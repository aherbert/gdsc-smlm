/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */

package uk.ac.sussex.gdsc.smlm.results;

import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Iterator;
import java.util.List;
import java.util.function.BiPredicate;
import java.util.function.Function;
import java.util.function.IntConsumer;
import java.util.function.IntPredicate;
import java.util.function.IntUnaryOperator;
import java.util.function.ToIntFunction;
import java.util.stream.Stream;
import uk.ac.sussex.gdsc.core.data.utils.ConversionException;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.logging.NullTrackProgress;
import uk.ac.sussex.gdsc.core.logging.TrackProgress;
import uk.ac.sussex.gdsc.core.match.Matchings;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.ValidationUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationHelper;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.data.config.PsfHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.results.count.Counter;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;

/**
 * Trace localisations through a time stack to identify single molecules.
 */
public class TraceManager {
  /**
   * Set the mode used to search backwards in time.
   */
  public enum TraceMode {
    //@formatter:off
    /**
     * Search the latest localisations first. This is equivalent to a downwards search.
     * When a localisation is found no more time points will be searched.
     */
    LATEST_FORERUNNER{ @Override
    public String getName() { return "Latest forerunner"; }},
    /**
     * Search the earliest localisations first. This is equivalent to a depth first search.
     * When a localisation is found no more time points will be searched.
     */
    EARLIEST_FORERUNNER{ @Override
    public String getName() { return "Earliest forerunner"; }},
    /**
     * Search all time points within the distance threshold. This is slower as all time points
     * are searched. It is equivalent to single-linkage clustering with a time constraint on
     * joined localisations.
     */
    SINGLE_LINKAGE{ @Override
    public String getName() { return "Single linkage"; }};
    //@formatter:on

    @Override
    public String toString() {
      return getName();
    }

    /**
     * Gets the name.
     *
     * @return the name
     */
    public abstract String getName();
  }

  private final MemoryPeakResults results;
  private final Localisation[] localisations;
  private Localisation[] endLocalisations;
  private int[] indexes;
  private int[] endIndexes;
  private int[] maxTime;
  private int totalTraces;
  private int totalFiltered;
  private TrackProgress tracker;
  private int activationFrameInterval;
  private int activationFrameWindow;
  private double distanceExclusion;
  private boolean filterActivationFrames;
  private TraceMode traceMode = TraceMode.LATEST_FORERUNNER;
  private int pulseInterval;

  /**
   * Contains trace localisation data.
   */
  private static class Localisation {
    int time;
    int endTime;
    int id;
    int trace;
    float x;
    float y;

    Localisation(int id, int time, int endTime, float x, float y) {
      ValidationUtils.checkArgument(endTime >= time, "End time (%d) is before the start time (%d)",
          endTime, time);

      this.time = time;
      this.endTime = endTime;
      this.id = id;
      this.x = x;
      this.y = y;
    }

    float distance2(Localisation other) {
      final float dx = x - other.x;
      final float dy = y - other.y;
      return dx * dx + dy * dy;
    }
  }

  /**
   * Contains a trajectory.
   */
  private static class Trajectory {
    int id;
    Localisation tail;

    Trajectory(int id, Localisation point) {
      this.id = id;
      add(point);
    }

    void add(Localisation point) {
      tail = point;
      point.trace = id;
    }

    double distance2(Localisation other) {
      final double dx = (double) tail.x - other.x;
      final double dy = (double) tail.y - other.y;
      return dx * dx + dy * dy;
    }

    int endTime() {
      return tail.endTime;
    }
  }

  /**
   * Iterator over a single list.
   */
  private static class SingletonIterator implements Iterator<List<Trajectory>> {
    private final List<Trajectory> traces;
    private boolean hasNext = true;

    /**
     * Create an instance.
     *
     * @param traces the traces
     */
    SingletonIterator(List<Trajectory> traces) {
      this.traces = traces;
    }

    @Override
    public boolean hasNext() {
      return hasNext;
    }

    @Override
    public List<Trajectory> next() {
      // Note: no check that hasNext is false and don't bother with NoSuchElementException
      hasNext = false;
      return traces;
    }
  }

  /**
   * Iterator over a list grouped by ascending end time order.
   */
  private static class EarliestIterator implements Iterator<List<Trajectory>> {
    private final List<Trajectory> traces;
    private int index;

    /**
     * Create an instance.
     *
     * @param traces the traces (assumed to not be empty)
     */
    EarliestIterator(List<Trajectory> traces) {
      this.traces = traces;
    }

    @Override
    public boolean hasNext() {
      return index < traces.size();
    }

    @Override
    public List<Trajectory> next() {
      // Note: no check that hasNext is false and don't bother with NoSuchElementException
      final int from = index;
      final int end = traces.size();
      final int t = traces.get(from).endTime();
      int i = from + 1;
      while (i < end && t == traces.get(i).endTime()) {
        i++;
      }
      index = i;
      return traces.subList(from, i);
    }
  }

  /**
   * Iterator over a list grouped by descending end time order.
   */
  private static class LatestIterator implements Iterator<List<Trajectory>> {
    private final List<Trajectory> traces;
    private int index;

    /**
     * Create an instance.
     *
     * @param traces the traces (assumed to not be empty)
     */
    LatestIterator(List<Trajectory> traces) {
      this.traces = traces;
      index = traces.size() - 1;
    }

    @Override
    public boolean hasNext() {
      return index >= 0;
    }

    @Override
    public List<Trajectory> next() {
      // Note: no check that hasNext is false and don't bother with NoSuchElementException
      final int to = index;
      final int t = traces.get(to).endTime();
      int i = to - 1;
      while (i >= 0 && t == traces.get(i).endTime()) {
        i--;
      }
      index = i;
      return traces.subList(i + 1, to + 1);
    }
  }

  /**
   * Instantiates a new trace manager.
   *
   * @param results the results
   * @throws IllegalArgumentException if results are null or empty
   */
  public TraceManager(final MemoryPeakResults results) {
    if (results == null || results.size() == 0) {
      throw new IllegalArgumentException("Results are null or empty");
    }
    this.results = results;

    // Assign localisations. We use the native result units to avoid conversion exceptions.
    localisations = new Localisation[results.size()];
    final Counter counter = new Counter();
    results.forEach((PeakResultProcedure) result -> {
      final int id = counter.getAndIncrement();
      localisations[id] = new Localisation(id, result.getFrame(), result.getEndFrame(),
          result.getXPosition(), result.getYPosition());
    });

    totalTraces = localisations.length;

    // Sort by start time
    Arrays.sort(localisations, (o1, o2) -> Integer.compare(o1.time, o2.time));

    // The algorithm assumes minT is positive
    if (localisations[0].time < 0) {
      throw new IllegalArgumentException("Lowest start time is negative");
    }

    // TODO - Assign a more efficient localisation representation using a grid
  }

  /**
   * Create a look-up table of the starting index in the localisations array for each possible time
   * point. This allows looping over all localisations for a given t using the index index[t] <= i <
   * index[t+1]
   *
   * @param localisations the localisations
   * @param timeFunction the time function
   * @return the int[]
   */
  private static int[] createTimeIndexes(Localisation[] localisations,
      ToIntFunction<Localisation> timeFunction) {
    final int maxT = timeFunction.applyAsInt(localisations[localisations.length - 1]);

    final int[] indexes = new int[maxT + 2];
    int time = -1;
    for (int i = 0; i < localisations.length; i++) {
      while (time < timeFunction.applyAsInt(localisations[i])) {
        indexes[++time] = i;
      }
    }
    indexes[maxT + 1] = localisations.length;
    return indexes;
  }

  private void createTimeIndexes() {
    if (endLocalisations == null) {
      // Build a second localisations list sorted by end time
      endLocalisations = localisations.clone();
      Arrays.sort(endLocalisations, (o1, o2) -> Integer.compare(o1.endTime, o2.endTime));

      // Create a look-up table of the starting index in each localisations array for each possible
      // time point. This allows looping over all localisations for a given t using the index
      // index[t] <= i < index[t+1]

      indexes = createTimeIndexes(localisations, point -> point.time);
      endIndexes = createTimeIndexes(endLocalisations, point -> point.endTime);
    }
  }

  /**
   * Trace localisations across frames that are the same molecule.
   *
   * <p>Any spot that occurred within time threshold and distance threshold of a previous spot is
   * grouped into the same trace as that previous spot. The resulting trace is assigned a spatial
   * position equal to the centroid position of all the spots included in the trace.
   *
   * <p>See Coltharp, et al. Accurate Construction of Photoactivated Localization Microscopy (PALM)
   * Images for Quantitative Measurements (2012). PLoS One. 7(12): e51725. DOI:
   * http://dx.doi.org/10.1371%2Fjournal.pone.0051725
   *
   * <p>Note: The actual traces representing molecules can be obtained by calling
   * {@link #getTraces()}
   *
   * @param distanceThreshold The distance threshold in the native units of the results
   * @param timeThreshold The time threshold in frames
   * @return The number of traces
   */
  public int traceMolecules(final double distanceThreshold, final int timeThreshold) {
    if (timeThreshold <= 0 || distanceThreshold < 0) {
      totalFiltered = 0;
      totalTraces = localisations.length;
      return totalTraces;
    }

    final TrackProgress tp = NullTrackProgress.createIfNull(this.tracker);
    tp.progress(0);

    final LocalList<Trajectory> traces = new LocalList<>(100);

    final IntPredicate ignoreActivation = createIgnoreActivationFilter();
    final IntConsumer updateTrajectories = createUpdateTrajectoriesFunction(traces, timeThreshold,
        pulseInterval, localisations[0].time);
    final BiPredicate<List<Trajectory>, Localisation> exclusionFilter =
        createExclusionFilter(distanceThreshold, distanceExclusion, ignoreActivation);
    final Function<List<Trajectory>, Iterator<List<Trajectory>>> toIterator =
        createIteratorFunction(timeThreshold, traceMode);

    int nextIndex = initialiseTrajectories(traces, ignoreActivation);

    // Now process the remaining frames, comparing them to previous frames
    LocalList<Localisation> candidates = new LocalList<>(100);
    LocalList<Localisation> targetCandidates = new LocalList<>(100);
    final double threshold = distanceThreshold * distanceThreshold;

    while (nextIndex < localisations.length) {
      tp.progress(nextIndex, localisations.length);

      // Identify the next localisations
      final int t = localisations[nextIndex].time;
      nextIndex = extractCandidates(nextIndex, candidates);

      updateTrajectories.accept(t);

      // Respect the exclusion distance. This removes any candidates that could be
      // matched a trace but have another trace within the exclusion distance.
      if (exclusionFilter != null) {
        candidates.removeIf(point -> exclusionFilter.test(traces, point));
        if (candidates.isEmpty()) {
          continue;
        }
      }

      // Match candidates to existing traces.
      // The active traces are iterated in batches based on the trace mode.
      // Traces are maintained in ascending end time order. Ignore any with
      // an end time effectively in the future.
      final int index = traces.findLastIndex(trace -> trace.endTime() < t);
      if (index >= 0) {
        final List<Trajectory> targetTraces =
            index + 1 == traces.size() ? traces : traces.subList(0, index + 1);
        Iterator<List<Trajectory>> iter = toIterator.apply(targetTraces);
        while (iter.hasNext()) {
          List<Trajectory> activeTraces = iter.next();

          // Rotate the candidate lists
          final LocalList<Localisation> tmp = targetCandidates;
          targetCandidates = candidates;
          candidates = tmp;
          candidates.clear();

          Matchings.nearestNeighbour(activeTraces, targetCandidates, Trajectory::distance2,
              threshold, (trace, point) -> {
                trace.add(point);
                // Ensure we maintain the time range for this trace
                maxTime[point.trace] = point.time;
              }, null, candidates::add);

          if (candidates.isEmpty()) {
            break;
          }
        }
      }

      // Create new traces from unmatched candidates
      if (!candidates.isEmpty()) {
        final boolean ignore = ignoreActivation.test(t);
        candidates.forEach(point -> addTrace(ignore, point, traces));
      }
    }

    tp.progress(1.0);

    return getTotalTraces();
  }

  private void addSingleTrace(boolean ignore, Localisation point) {
    totalTraces++;
    maxTime[totalTraces] = point.endTime;
    point.trace = totalTraces;
    if (ignore) {
      // Count the number of traces that will be filtered
      // (i.e. the time is not within an activation window).
      totalFiltered++;
    }
  }

  private void addTrace(boolean ignore, Localisation point, List<Trajectory> traces) {
    totalTraces++;
    maxTime[totalTraces] = point.endTime;
    if (ignore) {
      // Count the number of traces that will be filtered
      // (i.e. the time is not within an activation window).
      totalFiltered++;
      point.trace = totalTraces;
    } else {
      // Add to the active traces
      traces.add(new Trajectory(totalTraces, point));
    }
  }

  /**
   * Initialise the first traces using the first frame.
   *
   * @param traces the traces
   * @param ignoreActivation the ignore activation filter
   * @return the index in the localisations for the next time frame
   */
  private int initialiseTrajectories(LocalList<Trajectory> traces, IntPredicate ignoreActivation) {
    // Used to track the highest frame containing spots for a trace
    totalTraces = totalFiltered = 0;
    maxTime = new int[localisations.length + 1];

    // Initialise the first traces using the first frame
    final int t = localisations[0].time;
    final boolean ignore = ignoreActivation.test(t);
    addTrace(ignore, localisations[0], traces);
    // Add all with the same time
    for (int index = 1; index < localisations.length; index++) {
      if (localisations[index].time != t) {
        return index;
      }
      addTrace(ignore, localisations[index], traces);
    }
    return localisations.length;
  }

  /**
   * Creates the filter to ignore beginning a trajectory if the time does not occur in the
   * activation window. The predictate will return false if filtering is not enabled.
   *
   * @return the predicate
   */
  private IntPredicate createIgnoreActivationFilter() {
    if (filterActivationFrames) {
      final int interval = activationFrameInterval;
      final int window = activationFrameWindow;
      return time -> time % interval > window;
    }
    return point -> false;
  }

  /**
   * Creates a function to update the list of current trajectories. The list of trajectories is
   * modified in place based on the current time.
   *
   * @param traces the traces
   * @param timeThreshold the time threshold
   * @param pulseInterval the pulse interval
   * @param time the time of the first localisation
   * @return the function (for the current time)
   */
  private static IntConsumer createUpdateTrajectoriesFunction(LocalList<Trajectory> traces,
      int timeThreshold, int pulseInterval, int time) {
    // Support for splitting traces across pulse boundaries.
    if (pulseInterval > 0) {
      final IntUnaryOperator timeToPulse = t -> 1 + pulseInterval * ((t - 1) / pulseInterval);
      return new IntConsumer() {
        private int pulse = timeToPulse.applyAsInt(time);

        @Override
        public void accept(int t) {
          if (pulse != timeToPulse.applyAsInt(t)) {
            // Crossed a pulse boundary
            pulse = timeToPulse.applyAsInt(t);
            traces.clear();
          } else {
            updateTrajectories(timeThreshold, traces, t);
          }
        }
      };
    }
    return t -> updateTrajectories(timeThreshold, traces, t);
  }

  /**
   * Filter trajectories to remove those with an end time before the start of the time window.
   *
   * @param timeThreshold the time threshold
   * @param traces the traces
   * @param t the current time
   */
  private static void updateTrajectories(final int timeThreshold, LocalList<Trajectory> traces,
      final int t) {
    // Filter to those traces within the window
    final int pastT = t - timeThreshold;
    traces.removeIf(trace -> trace.endTime() < pastT);
    traces.sort((a, b) -> Integer.compare(a.endTime(), b.endTime()));
  }

  /**
   * Creates the exclusion filter. This will mark localisations to ignore if they match multiple
   * trajectories within the exclusion distance and at least 1 trajectory within the trace distance.
   *
   * <p>Filtered localisations are used to create a trace of a single localisation.
   *
   * @param distance the trace distance
   * @param upperDistance the exclusion distance
   * @return the predicate (or null)
   */
  private BiPredicate<List<Trajectory>, Localisation> createExclusionFilter(double distance,
      double upperDistance, IntPredicate ignoreActivation) {
    if (upperDistance >= distance) {
      final double lower = distance * distance;
      final double upper = upperDistance * upperDistance;
      return (traces, point) -> {
        boolean match = false;
        int exclusionCount = 0;
        for (final Trajectory trace : traces) {
          final double d = trace.distance2(point);
          if (d <= upper) {
            exclusionCount++;
            if (d <= lower) {
              match = true;
            }
            if (match && exclusionCount > 1) {
              addSingleTrace(ignoreActivation.test(point.time), point);
              return true;
            }
          }
        }
        return false;
      };
    }
    return null;
  }

  /**
   * Create a function to divide the active traces into subsets based on the trace mode. If the time
   * threshold is 1 then the trace mode is ignored and the entire set of active traces is returned.
   *
   * @param traces the traces
   * @return the iterator
   */
  private static Function<List<Trajectory>, Iterator<List<Trajectory>>>
      createIteratorFunction(int timeThreshold, TraceMode traceMode) {
    if (traceMode == TraceMode.SINGLE_LINKAGE || timeThreshold == 1) {
      return SingletonIterator::new;
    }
    // Return in chunks using the end time
    return traces -> iteratorByTime(traces, traceMode == TraceMode.LATEST_FORERUNNER);
  }

  /**
   * Create an iterable over the active traces (assumed to be sorted by end time).
   *
   * @param traces the traces
   * @param reversed true to return in reversed order
   * @return the iterable
   */
  private static Iterator<List<Trajectory>> iteratorByTime(List<Trajectory> traces,
      boolean reversed) {
    // If all the same then use a single list
    if (traces.isEmpty() || traces.get(0).endTime() == traces.get(traces.size() - 1).endTime()) {
      return new SingletonIterator(traces);
    }
    return reversed ? new LatestIterator(traces) : new EarliestIterator(traces);
  }

  /**
   * Extract the candidates into the provided list.
   *
   * @param currentIndex the current index
   * @param candidates the candidates
   * @return the index in the localisations for the next time frame
   */
  private int extractCandidates(int currentIndex, LocalList<Localisation> candidates) {
    candidates.clear();
    final int t = localisations[currentIndex].time;
    candidates.add(localisations[currentIndex]);
    for (int index = currentIndex + 1; index < localisations.length; index++) {
      if (localisations[index].time != t) {
        return index;
      }
      candidates.add(localisations[index]);
    }
    return localisations.length;
  }

  /**
   * Gets the traces that have been found using {@link #traceMolecules(double, int)}.
   *
   * @return The traces
   */
  public Trace[] getTraces() {
    final PeakResult[] peakResults = results.toArray();

    // No tracing yet performed or no thresholds
    if (totalTraces == localisations.length) {
      Stream<PeakResult> stream = Arrays.stream(peakResults);
      if (filterActivationFrames) {
        final IntPredicate ignoreActivation = createIgnoreActivationFilter();
        stream = stream.filter(peakResult -> !ignoreActivation.test(peakResult.getFrame()));
      }
      final int[] id = {0};
      return stream.map(r -> {
        final Trace trace = new Trace(r);
        trace.setId(++id[0]);
        return trace;
      }).toArray(Trace[]::new);
    }

    if (tracker != null) {
      tracker.progress(0);
    }

    // Build the list of traces
    final Trace[] traces = new Trace[getTotalTraces()];
    int count = 0;

    // for (int index = 0; index < localisations.length; index++)
    // if (localisations[index].trace == 0)
    // System.out.printf("error @ %d\n", index);

    // Since the trace numbers are allocated by processing the spots in frames, each frame can have
    // trace number out-of-order. This occurs if re-allocation has been performed,
    // e.g. [1,2,2,1,3] => [1,2,5,4,3] when spots in group 1 are reallocated before spots in group
    // 2.

    // TODO - Update localisation to support a single linked list of the trace using a nextIndex.
    // Created the traces by following the list.
    // This should also mark any traces that were filtered so we do not require the
    // filtering to be done here.

    final BitSet processedTraces = new BitSet(totalTraces + 1);
    final IntPredicate ignoreActivation = createIgnoreActivationFilter();
    for (int i = 0; i < localisations.length; i++) {
      if (tracker != null && (i & 0xff) == 0) {
        tracker.progress(i, localisations.length);
      }

      final int traceId = localisations[i].trace;

      if (processedTraces.get(traceId) || ignoreActivation.test(localisations[i].time)) {
        continue;
      }
      processedTraces.set(traceId);

      final PeakResult peakResult = peakResults[localisations[i].id];

      final Trace nextTrace = new Trace(peakResult);
      nextTrace.setId(traceId);
      final int tLimit = maxTime[traceId];

      // Check if the trace has later frames
      if (tLimit > localisations[i].time) {
        for (int j = i + 1; j < localisations.length; j++) {
          if (localisations[j].time > tLimit) {
            break;
          }
          if (localisations[j].trace == traceId) {
            nextTrace.add(peakResults[localisations[j].id]);
          }
        }
      }

      traces[count++] = nextTrace;
    }

    if (tracker != null) {
      tracker.progress(1.0);
    }

    return traces;
  }

  /**
   * Convert a list of traces into peak results. The centroid of each trace is used as the
   * coordinates. The standard deviation of positions from the centre is used as the width. The
   * amplitude is the average from all the peaks in the trace.
   *
   * @param traces the traces
   * @return the peak results
   */
  public static MemoryPeakResults toCentroidPeakResults(final Trace[] traces) {
    final int capacity = 1 + ((traces != null) ? traces.length : 0);
    final MemoryPeakResults results = new MemoryPeakResults(capacity);
    if (traces != null) {
      for (int i = 0; i < traces.length; i++) {
        final PeakResult result = traces[i].getHead();
        if (result == null) {
          continue;
        }
        final float[] centroid = traces[i].getCentroid();
        final float sd = traces[i].getStandardDeviation();
        final float background = 0;
        final float signal = (float) traces[i].getSignal();
        final float[] params = {background, signal, 0, centroid[0], centroid[1], sd, sd};
        final int endFrame = traces[i].getTail().getEndFrame();
        results.add(new ExtendedPeakResult(result.getFrame(), result.getOrigX(), result.getOrigY(),
            result.getOrigValue(), 0, 0, 0, params, null, endFrame, i + 1));
      }
    }
    return results;
  }

  /**
   * Convert a list of traces into peak results. The signal weighted centroid of each trace is used
   * as the coordinates. The weighted localisation precision is used as the precision. The signal is
   * the sum from all the peaks in the trace.
   *
   * <p>If the trace is empty it is ignored.
   *
   * <p>The results will only have the standard parameters set [Background,Intensity,X,Y,Z] so a
   * custom PSF is added to the results.
   *
   * @param traces the traces
   * @param calibration the calibration
   * @return the peak results
   * @see PSFType#CUSTOM
   */
  public static MemoryPeakResults toCentroidPeakResults(final Trace[] traces,
      final Calibration calibration) {
    final int capacity = 1 + ((traces != null) ? traces.length : 0);
    final MemoryPeakResults results = new MemoryPeakResults(capacity);
    results.setCalibration(calibration);
    results.setPsf(PsfHelper.create(PSFType.CUSTOM));
    if (traces != null) {
      TypeConverter<DistanceUnit> converter = null;
      if (calibration != null) {
        try {
          converter = CalibrationHelper.getDistanceConverter(calibration, DistanceUnit.NM);
        } catch (final ConversionException ignored) {
          // Ignore
        }
      }

      // Ensure all results are added as extended peak results with their trace ID.
      for (final Trace trace : traces) {
        if (trace == null || trace.size() == 0) {
          continue;
        }

        final PeakResult result = trace.getHead();
        if (trace.size() == 1) {
          final AttributePeakResult peakResult = new AttributePeakResult(result.getFrame(),
              result.getOrigX(), result.getOrigY(), result.getOrigValue(), 0, result.getNoise(),
              result.getMeanIntensity(), result.getParameters(), null);
          peakResult.setId(trace.getId());
          peakResult.setEndFrame(result.getEndFrame());
          if (converter != null) {
            peakResult.setPrecision(trace.getLocalisationPrecision(converter));
          }
          results.add(peakResult);
        } else {
          trace.sort();
          trace.resetCentroid();
          final float[] centroid = trace.getCentroid();
          float background = 0;
          double noise = 0;
          for (final PeakResult r : trace.getPoints().toArray()) {
            noise += r.getNoise() * r.getNoise();
            background += r.getBackground();
          }
          noise = Math.sqrt(noise);
          background /= trace.size();
          final double signal = trace.getSignal();
          final int endFrame = trace.getTail().getEndFrame();
          final AttributePeakResult peakResult =
              new AttributePeakResult(result.getFrame(), centroid[0], centroid[1], (float) signal);
          // Build standard peak data
          peakResult.setBackground(background);
          peakResult.setNoise((float) noise);
          // These could be weighted, at the moment we use the first peak
          peakResult.setOrigX(result.getOrigX());
          peakResult.setOrigY(result.getOrigY());
          peakResult.setOrigValue(result.getOrigValue());

          peakResult.setId(trace.getId());
          peakResult.setEndFrame(endFrame);
          if (converter != null) {
            peakResult.setPrecision(trace.getLocalisationPrecision(converter));
          }
          results.add(peakResult);
        }
      }
    }
    return results;
  }

  /**
   * Convert a list of traces into peak results setting the trace ID into the results.
   *
   * <p>If the trace is empty it is ignored.
   *
   * @param traces the traces
   * @param calibration the calibration
   * @return the peak results
   */
  public static MemoryPeakResults toPeakResults(final Trace[] traces,
      final Calibration calibration) {
    return toPeakResults(traces, calibration, false);
  }

  /**
   * Convert a list of traces into peak results setting the trace ID into the results.
   *
   * <p>If the trace is empty it is ignored.
   *
   * @param traces the traces
   * @param calibration the calibration
   * @param newId Set to true to use a new ID for each trace
   * @return the peak results
   */
  public static MemoryPeakResults toPeakResults(final Trace[] traces, final Calibration calibration,
      boolean newId) {
    final int capacity = 1 + ((traces != null) ? traces.length : 0);
    final MemoryPeakResults results = new MemoryPeakResults(capacity);
    results.setCalibration(calibration);
    if (traces != null) {
      // Ensure all results are added as extended peak results with their trace ID.
      int id = 0;
      for (final Trace trace : traces) {
        if (trace == null || trace.size() == 0) {
          continue;
        }

        final int traceId = (newId) ? ++id : trace.getId();
        trace.getPoints().forEach((PeakResultProcedure) result -> {
          final AttributePeakResult copy = new AttributePeakResult(result);
          copy.setId(traceId);
          results.add(copy);
        });
      }
    }
    return results;
  }

  /**
   * Convert a list of traces into peak results. The signal weighted centroid of each trace is used
   * as the coordinates. The weighted localisation precision is used as the width. The amplitude is
   * the average from all the peaks in the trace.
   *
   * <p>Uses the title and bounds from the constructor peak results. The title has the word 'Traced
   * Centroids' appended.
   *
   * @param traces the traces
   * @return the peak results
   */
  public MemoryPeakResults convertToCentroidPeakResults(final Trace[] traces) {
    return convertToCentroidPeakResults(results, traces);
  }

  /**
   * Convert a list of traces into peak results. The signal weighted centroid of each trace is used
   * as the coordinates. The weighted localisation precision is used as the width. The amplitude is
   * the average from all the peaks in the trace.
   *
   * <p>Uses the title and bounds from the provided peak results. The title has the word 'Traced
   * Centroids' appended.
   *
   * @param source the source
   * @param traces the traces
   * @return the peak results
   */
  public static MemoryPeakResults convertToCentroidPeakResults(MemoryPeakResults source,
      final Trace[] traces) {
    final MemoryPeakResults results = toCentroidPeakResults(traces, source.getCalibration());
    // Do not copy the PSF as it is invalid for centroids
    final PSF psf = results.getPsf();
    results.copySettings(source);
    results.setPsf(psf);
    // Change name
    results.setName(source.getSource() + " Traced Centroids");
    // TODO - Add the tracing settings
    return results;
  }

  /**
   * Convert a list of traces into peak results.
   *
   * <p>Uses the title and bounds from the constructor peak results. The title has the word 'Traced'
   * appended.
   *
   * @param traces the traces
   * @return the peak results
   */
  public MemoryPeakResults convertToPeakResults(final Trace[] traces) {
    return convertToPeakResults(results, traces);
  }

  /**
   * Convert a list of traces into peak results.
   *
   * <p>Uses the title and bounds from the provided peak results. The title has the word 'Traced'
   * appended.
   *
   * @param source the source
   * @param traces the traces
   * @return the peak results
   */
  public static MemoryPeakResults convertToPeakResults(MemoryPeakResults source,
      final Trace[] traces) {
    final MemoryPeakResults results = toPeakResults(traces, source.getCalibration());
    results.copySettings(source);
    // Change name
    results.setName(source.getSource() + " Traced");
    // TODO - Add the tracing settings
    return results;
  }

  /**
   * Gets the tracker.
   *
   * @return the tracker.
   */
  public TrackProgress getTracker() {
    return tracker;
  }

  /**
   * Sets the tracker.
   *
   * @param tracker the tracker to set
   */
  public void setTracker(TrackProgress tracker) {
    this.tracker = tracker;
  }

  /**
   * Gets the activation frame interval.
   *
   * @return the activation frame interval
   * @see #setActivationFrameInterval(int)
   */
  public int getActivationFrameInterval() {
    return activationFrameInterval;
  }

  /**
   * Set the interval at which the activation laser is used. These form staging points for the
   * traces.
   *
   * @param activationFrameInterval the activationFrameInterval to set
   */
  public void setActivationFrameInterval(int activationFrameInterval) {
    this.activationFrameInterval = activationFrameInterval;
    resetFilterActivationFramesFlag();
  }

  /**
   * Gets the activation frame window.
   *
   * @return the activation frame window
   * @see #setActivationFrameWindow(int)
   */
  public int getActivationFrameWindow() {
    return activationFrameWindow;
  }

  /**
   * Set the window after the activation pulse that will be used for traces. Any trace that does not
   * start within this window will be discarded.
   *
   * @param activationFrameWindow the activationFrameWindow to set
   */
  public void setActivationFrameWindow(int activationFrameWindow) {
    this.activationFrameWindow = activationFrameWindow;
    resetFilterActivationFramesFlag();
  }

  private void resetFilterActivationFramesFlag() {
    filterActivationFrames = (activationFrameInterval > 1 && activationFrameWindow > 0);
  }

  /**
   * Gets the trace mode.
   *
   * @return the trace mode
   */
  public TraceMode getTraceMode() {
    return traceMode;
  }

  /**
   * Sets the trace mode.
   *
   * @param traceMode the trace mode to set
   */
  public void setTraceMode(TraceMode traceMode) {
    this.traceMode = traceMode;
  }

  /**
   * Gets the pulse interval.
   *
   * @return the pulse interval.
   */
  public int getPulseInterval() {
    return pulseInterval;
  }

  /**
   * Set a pulse interval. Traces will only be created by joining localisations within each pulse.
   *
   * @param pulseInterval the pulse interval
   */
  public void setPulseInterval(int pulseInterval) {
    this.pulseInterval = Math.max(0, pulseInterval);
  }

  /**
   * Gets the distance exclusion.
   *
   * @return the distance exclusion
   * @see #setDistanceExclusion(double)
   */
  public double getDistanceExclusion() {
    return distanceExclusion;
  }

  /**
   * Set the minimum distance the next candidate spot must be in the same frame, i.e. choose
   * localisations closer than the distance threshold but no other spots are closer than this
   * distance exclusion
   *
   * <p>If less that the tracing distance threshold this value is ignored.
   *
   * @param distanceExclusion the distance exclusion
   */
  public void setDistanceExclusion(double distanceExclusion) {
    this.distanceExclusion = distanceExclusion;
  }

  /**
   * Gets the total traces from the last call of {@link #traceMolecules(double, int)}.
   *
   * @return the total traces
   */
  public int getTotalTraces() {
    return totalTraces - totalFiltered;
  }

  /**
   * Return the number of traces that were filtered since the trace was first activated outside the
   * configured activation window.
   *
   * @return the total filtered
   */
  public int getTotalFiltered() {
    return totalFiltered;
  }

  /**
   * Filter the traces that start during an activation frame.
   *
   * @param traces the traces
   * @param activationFrameInterval the interval at which the activation laser is used
   * @return the filtered traces
   */
  public Trace[] filterTraces(Trace[] traces, int activationFrameInterval) {
    return Arrays.stream(traces).filter(trace -> {
      final PeakResult r = trace.getHead();
      return r != null && r.getFrame() % activationFrameInterval == 1;
    }).toArray(Trace[]::new);
  }

  /**
   * Find the neighbour for each result within the given time and distance thresholds. The neighbour
   * with the strongest signal is selected.
   *
   * @param distanceThreshold The distance threshold in the native units of the results
   * @param timeThreshold The time threshold in frames
   * @return A list of traces containing the molecule and neighbour. If no neighbour is found then
   *         the trace will contain a single result
   */
  public Trace[] findNeighbours(final double distanceThreshold, final int timeThreshold) {
    if (distanceThreshold <= 0 || timeThreshold <= 0) {
      throw new IllegalArgumentException("Distance and time thresholds must be positive");
    }

    createTimeIndexes();

    final Localisation[] localisations = this.localisations;
    final Localisation[] endLocalisations = this.endLocalisations;
    final int[] indexes = this.indexes;
    final int[] endIndexes = this.endIndexes;

    final Trace[] neighbours = new Trace[results.size()];
    final PeakResult[] peakResults = results.toArray();

    final float dThresh2 = (float) (distanceThreshold * distanceThreshold);

    final TrackProgress tp = NullTrackProgress.createIfNull(this.tracker);
    tp.progress(0);

    // Initialise
    int nextIndex = 0;

    // Now process all the frames, comparing them to previous and future frames
    while (nextIndex < localisations.length) {
      tp.progress(nextIndex, localisations.length);

      final int currentIndex = nextIndex;
      final int t = localisations[currentIndex].time;
      // Look for the index for the next time-frame
      for (int tt = t + 1; tt < indexes.length; tt++) {
        nextIndex = indexes[tt];
        if (nextIndex != currentIndex) {
          break;
        }
      }
      final int pastEndIndex = endIndexes[Math.max(t - timeThreshold, 0)];
      final int currentEndIndex = endIndexes[t];
      final int futureIndex =
          Math.max(nextIndex, indexes[Math.min(t + 1 + timeThreshold, indexes.length - 1)]);

      // Process all spots from this frame.
      for (int index = currentIndex; index < nextIndex; index++) {
        final Localisation l = localisations[index];

        float maxSignal = 0;
        int neighbour = -1;

        // Look back
        for (int i = pastEndIndex; i < currentEndIndex; i++) {
          if (l.distance2(endLocalisations[i]) < dThresh2) {
            final float signal = peakResults[endLocalisations[i].id].getIntensity();
            if (maxSignal < signal) {
              maxSignal = signal;
              neighbour = endLocalisations[i].id;
            }
          }
        }

        // Look forward
        for (int i = nextIndex; i < futureIndex; i++) {
          if (l.distance2(localisations[i]) < dThresh2) {
            final float signal = peakResults[localisations[i].id].getIntensity();
            if (maxSignal < signal) {
              maxSignal = signal;
              neighbour = localisations[i].id;
            }
          }
        }

        // Assign
        final Trace trace = new Trace(peakResults[l.id]);
        if (neighbour > -1) {
          trace.add(peakResults[neighbour]);
        }
        neighbours[index] = trace;
      }
    }

    tp.progress(1.0);

    return neighbours;
  }

  /**
   * Collect all peak results with the same ID into traces. IDs of zero are ignored.
   *
   * @param results the results
   * @return The traces
   */
  public static Trace[] convert(MemoryPeakResults results) {
    if (results == null || results.size() == 0) {
      return new Trace[0];
    }

    final PeakResult[] list = results.toArray();
    // Find the max trace ID
    int max = 0;
    for (final PeakResult result : list) {
      if (max < result.getId()) {
        max = result.getId();
      }
    }
    if (max == 0) {
      return new Trace[0];
    }

    if (max < 10000) {
      // Small set of IDs so just use an array look-up table
      final Trace[] traces = new Trace[max + 1];
      for (final PeakResult result : list) {
        final int id = result.getId();
        if (id > 0) {
          if (traces[id] == null) {
            traces[id] = new Trace();
            traces[id].setId(id);
          }
          traces[id].add(result);
        }
      }

      // Consolidate to remove empty positions
      int count = 0;
      for (int i = 1; i < traces.length; i++) {
        if (traces[i] != null) {
          traces[count++] = traces[i];
        }
      }
      return Arrays.copyOf(traces, count);
    }
    // Use a map when the size is potentially large
    final Int2ObjectOpenHashMap<Trace> map = new Int2ObjectOpenHashMap<>();
    for (final PeakResult result : list) {
      final int id = result.getId();
      if (id > 0) {
        map.computeIfAbsent(id, i -> {
          final Trace trace = new Trace();
          trace.setId(i);
          return trace;
        }).add(result);
      }
    }

    // Extract the traces
    final Trace[] traces = map.values().toArray(new Trace[map.size()]);
    Arrays.sort(traces, (t1, t2) -> Integer.compare(t1.getId(), t2.getId()));
    return traces;
  }
}
