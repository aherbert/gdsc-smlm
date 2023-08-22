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
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import uk.ac.sussex.gdsc.core.data.utils.ConversionException;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.logging.TrackProgress;
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
  private final Localisation[] endLocalisations;
  private final int[] indexes;
  private final int[] endIndexes;
  private int[] maxTime;
  private int totalTraces;
  private int totalFiltered;
  private float distanceThreshSquared;
  private float distanceExclusionSquared;
  private TrackProgress tracker;
  private int activationFrameInterval;
  private int activationFrameWindow;
  private double distanceExclusion;
  private boolean filterActivationFrames;
  private TraceMode traceMode = TraceMode.LATEST_FORERUNNER;
  private int pulseInterval;

  /**
   * The distance between the localisation and its assigned forerunner.
   *
   * <p>Set in {@link #findForerunner(int, int, int)} and
   * {@link #findAlternativeForerunner(int, int, int, int, int[])}.
   */
  private float minD;

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
   * Contains trace assignment data.
   */
  private static class Assignment {
    int index;
    float distance;
    int traceId;

    Assignment(int index, float distance, int traceId) {
      this.index = index;
      this.distance = distance;
      this.traceId = traceId;
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

    // Build a second localisations list sorted by end time
    endLocalisations = Arrays.copyOf(localisations, totalTraces);
    Arrays.sort(endLocalisations, (o1, o2) -> Integer.compare(o1.endTime, o2.endTime));

    // Create a look-up table of the starting index in each localisations array for each possible
    // time point. This allows looping over all localisations for a given t using the index
    // index[t] <= i < index[t+1]

    int maxT = localisations[totalTraces - 1].time;

    indexes = new int[maxT + 2];
    int time = -1;
    for (int i = 0; i < localisations.length; i++) {
      while (time < localisations[i].time) {
        indexes[++time] = i;
      }
    }
    indexes[maxT + 1] = totalTraces;

    maxT = endLocalisations[totalTraces - 1].endTime;
    endIndexes = new int[maxT + 2];
    time = -1;
    for (int i = 0; i < endLocalisations.length; i++) {
      while (time < endLocalisations[i].endTime) {
        endIndexes[++time] = i;
      }
    }
    endIndexes[maxT + 1] = totalTraces;

    // TODO - Assign a more efficient localisation representation using a grid
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
      totalTraces = localisations.length;
      return totalTraces;
    }

    totalTraces = totalFiltered = 0;
    distanceThreshSquared = (float) (distanceThreshold * distanceThreshold);
    distanceExclusionSquared =
        (distanceExclusion >= distanceThreshold) ? (float) (distanceExclusion * distanceExclusion)
            : 0;

    if (tracker != null) {
      tracker.progress(0);
    }

    // Used to track the highest frame containing spots for a trace
    maxTime = new int[localisations.length + 1];
    final int[] traceIdToLocalisationsIndexMap = new int[localisations.length + 1];

    // Initialise the first traces using the first frame
    int nextIndex = indexes[localisations[0].time + 1]; // findNextStartTimeIndex(0);
    for (int index = 0; index < nextIndex; index++) {
      localisations[index].trace = addTrace(index, traceIdToLocalisationsIndexMap, maxTime);
    }

    Assignment[] assigned = new Assignment[10];

    // Now process the remaining frames, comparing them to previous frames
    while (nextIndex < localisations.length) {
      if (tracker != null) {
        tracker.progress(nextIndex, localisations.length);
      }

      final int currentIndex = nextIndex;
      final int t = localisations[currentIndex].time;
      nextIndex = indexes[t + 1];
      int pastT = Math.max(t - timeThreshold, 0);
      if (pulseInterval > 0) {
        // Support for splitting traces across pulse boundaries. Simply round the
        // previous timepoint to the next pulse boundary. Assume pulses start at t=1
        final int intervalBoundary = 1 + pulseInterval * ((t - 1) / pulseInterval);
        if (pastT < intervalBoundary) {
          pastT = intervalBoundary;
        }
      }
      final int pastEndIndex = endIndexes[pastT];
      final int currentEndIndex = endIndexes[t];

      // If no previous spots within the time threshold then create new traces
      if (pastEndIndex == currentEndIndex) {
        for (int index = currentIndex; index < nextIndex; index++) {
          localisations[index].trace = addTrace(index, traceIdToLocalisationsIndexMap, maxTime);
        }
        continue;
      }

      // Check the allocated buffer is larger enough
      if (assigned.length < nextIndex - currentIndex) {
        assigned = new Assignment[nextIndex - currentIndex];
      }

      // Process all spots from this frame. Note if a spot is allocated to an existing trace.
      int assignedToTrace = 0;
      for (int index = currentIndex; index < nextIndex; index++) {
        final int traceId = findForerunner(index, pastEndIndex, currentEndIndex);
        if (traceId == 0) {
          localisations[index].trace = addTrace(index, traceIdToLocalisationsIndexMap, maxTime);
        } else {
          // Tentatively assign
          assigned[assignedToTrace++] = new Assignment(index, minD, traceId);
        }
      }

      if (assignedToTrace > 1) {
        // Check if duplicate allocations are made. Each trace can only
        // be allocated one localisation so in the event of a multiple
        // allocation then only the closest spot should be allocated
        final int[] dualAllocation = new int[assignedToTrace];
        final int[] ignore = new int[assignedToTrace];
        int ignoreCount = 0;

        // Only check for duplicates if two assignments are remaining
        boolean reSort = true;
        for (int i = 0; i < assignedToTrace - 1; i++) {
          // If the distance is negative then this can be skipped as it was a new trace
          // (allocated in a previous loop).
          if (assigned[i].distance < 0) {
            continue;
          }

          // Sort the remaining allocations by their distance
          if (reSort) {
            reSort = false;
            Arrays.sort(assigned, i, assignedToTrace,
                (o1, o2) -> Double.compare(o1.distance, o2.distance));
            // Check for new traces (allocated in a previous loop). These have distance <0 so will
            // be sorted to the front.
            if (assigned[i].distance < 0) {
              continue;
            }
          }

          int dualAllocationCount = 0;

          for (int j = i + 1; j < assignedToTrace; j++) {
            // Dual allocation
            if (assigned[i].traceId == assigned[j].traceId) {
              dualAllocation[dualAllocationCount++] = j;
            }
          }

          // This trace has been taken so ignore when finding alternatives
          ignore[ignoreCount++] = assigned[i].traceId;

          if (dualAllocationCount > 0) {
            // Re-allocate the other spots
            for (int a = 0; a < dualAllocationCount; a++) {
              final int j = dualAllocation[a];
              final int index = assigned[j].index;
              int traceId = findAlternativeForerunner(index, pastEndIndex, currentEndIndex,
                  ignoreCount, ignore);
              if (traceId == 0) {
                traceId = addTrace(index, traceIdToLocalisationsIndexMap, maxTime);
                // Mark to ignore
                assigned[j].distance = -1;
              } else {
                // Indicate that the distances have changed and a re-sort is needed
                reSort = true;
                assigned[j].distance = minD;
              }
              assigned[j].traceId = traceId;
            }
          }
          // Ensure nothing can be sorted ahead of this trace assignment
          assigned[i].distance = -1;
        }
      }

      // Assign the localisations
      for (int i = 0; i < assignedToTrace; i++) {
        localisations[assigned[i].index].trace = assigned[i].traceId;
        maxTime[assigned[i].traceId] = localisations[assigned[i].index].endTime;
      }
    }

    if (tracker != null) {
      tracker.progress(1.0);
    }

    return getTotalTraces();
  }

  private int addTrace(int index, int[] traceIdToLocalisationsIndexMap, int[] maxT) {
    if (filterActivationFrames
        // Count the number of traces that will be filtered
        // (i.e. the time is not within an activation window)
        && outsideActivationWindow(localisations[index].time)) {
      totalFiltered++;
    }

    final int traceId = ++totalTraces;
    traceIdToLocalisationsIndexMap[traceId] = index;
    maxT[traceId] = localisations[index].endTime;
    return traceId;
  }

  private boolean outsideActivationWindow(int time) {
    return time % activationFrameInterval > activationFrameWindow;
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
      if (filterActivationFrames) {
        final ArrayList<Trace> traces = new ArrayList<>();
        for (int i = 0; i < totalTraces; i++) {
          final PeakResult peakResult = peakResults[localisations[i].id];
          if (!outsideActivationWindow(peakResult.getFrame())) {
            traces.add(new Trace(peakResult));
          }
        }
        return traces.toArray(new Trace[0]);
      }
      final Trace[] traces = new Trace[localisations.length];
      for (int i = 0; i < traces.length; i++) {
        traces[i] = new Trace(peakResults[localisations[i].id]);
      }
      return traces;
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

    final IntOpenHashSet processedTraces = new IntOpenHashSet(traces.length);
    for (int i = 0; i < localisations.length; i++) {
      if (tracker != null && (i & 0xff) == 0) {
        tracker.progress(i, localisations.length);
      }

      final int traceId = localisations[i].trace;

      if (!processedTraces.add(traceId)
          || (filterActivationFrames && outsideActivationWindow(localisations[i].time))) {
        continue;
      }

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
   * From the given index, move forward to a localisation with a new start time frame. If no more
   * frames return the number of localisations.
   *
   * @param index the index
   * @return The index of the next time frame
   */
  @SuppressWarnings("unused")
  private int findNextStartTimeIndex(int index) {
    final int t = localisations[index].time;
    while (index < localisations.length && localisations[index].time <= t) {
      index++;
    }
    return index;
  }

  /**
   * From the given index, move forward to a localisation with a start time beyond the time
   * threshold. If no more frames return the number of localisations.
   *
   * @param index the index
   * @param timeThreshold the time threshold
   * @return The index of the next time frame
   */
  @SuppressWarnings("unused")
  private int findNextStartTimeIndex(int index, final int timeThreshold) {
    final int t = localisations[index].time + timeThreshold;
    while (index < localisations.length && localisations[index].time <= t) {
      index++;
    }
    return index;
  }

  /**
   * From the given index, move backward to the earliest localisations within the time threshold.
   *
   * @param index the index
   * @param timeThreshold the time threshold
   * @return The index of the earliest localisation within the time threshold
   */
  @SuppressWarnings("unused")
  private int findPastTimeIndex(int index, final int timeThreshold) {
    final int t = localisations[index].time - timeThreshold;
    while (index > 0) {
      index--;
      if (localisations[index].time < t) {
        index++; // Set back to within the time threshold
        break;
      }
    }
    return index;
  }

  /**
   * Find the earliest forerunner spot (from pastIndex to currentIndex) that is within the distance
   * threshold of the given spot. In the event that multiple forerunner spots from the same frame
   * are within the distance, assign the closest spot.
   *
   * @param index The index of the spot
   * @param pastIndex The index of the earliest forerunner spot
   * @param currentIndex The index of the first spot in the same frame (i.e. end of forerunner
   *        spots)
   * @return The assigned trace
   */
  private int findForerunner(final int index, final int pastIndex, final int currentIndex) {
    if (distanceExclusionSquared == 0) {
      return findForerunnerNoExclusion(index, pastIndex, currentIndex);
    }
    return findForerunnerWithExclusion(index, pastIndex, currentIndex);
  }

  /**
   * Find the earliest forerunner spot (from pastIndex to currentIndex) that is within the distance
   * threshold of the given spot. In the event that multiple forerunner spots from the same frame
   * are within the distance, assign the closest spot.
   *
   * @param index The index of the spot
   * @param pastIndex The index of the earliest forerunner spot
   * @param currentIndex The index of the first spot in the same frame (i.e. end of forerunner
   *        spots)
   * @return The assigned trace
   */
  private int findForerunnerNoExclusion(final int index, final int pastIndex,
      final int currentIndex) {
    final Localisation spot = localisations[index];
    if (traceMode == TraceMode.EARLIEST_FORERUNNER) {
      for (int i = pastIndex; i < currentIndex; i++) {
        final float d2 = spot.distance2(endLocalisations[i]);
        if (d2 <= distanceThreshSquared) {
          minD = d2;
          int trace = endLocalisations[i].trace;

          // Search all remaining spots that end in this time frame and pick the closest
          final int nextIndex = endIndexes[endLocalisations[i].endTime + 1];
          for (int ii = i + 1; ii < nextIndex; ii++) {
            final float dd2 = spot.distance2(endLocalisations[ii]);
            if (dd2 < minD) {
              minD = dd2;
              trace = endLocalisations[ii].trace;
            }
          }

          return trace;
        }
      }
    } else if (traceMode == TraceMode.LATEST_FORERUNNER) {
      for (int i = currentIndex; i-- > pastIndex;) {
        final float d2 = spot.distance2(endLocalisations[i]);
        if (d2 <= distanceThreshSquared) {
          minD = d2;
          int trace = endLocalisations[i].trace;

          // Search all remaining spots in this time frame and pick the closest
          final int previousIndex = endIndexes[endLocalisations[i].endTime];
          //// DEBUG
          // int previousIndex = i;
          //// Look for the index for the previous time-frame
          // while (previousIndex > 0 && endLocalisations[previousIndex-1].t ==
          //// endLocalisations[i].t)
          // previousIndex--;
          // if (previousIndex != endIndex[endLocalisations[i].endT])
          // {
          // System.out.printf("Error when tracing: %d != %d\n", previousIndex,
          // endIndex[endLocalisations[i].endT]);
          // }
          for (int ii = i; ii-- > previousIndex;) {
            final float dd2 = spot.distance2(endLocalisations[ii]);
            if (dd2 < minD) {
              minD = dd2;
              trace = endLocalisations[ii].trace;
            }
          }

          return trace;
        }
      }
    } else {
      // traceMode == TraceMode.SINGLE_LINKAGE

      // Find the closest spot
      minD = distanceThreshSquared;
      int minI = -1;
      for (int i = pastIndex; i < currentIndex; i++) {
        final float d2 = spot.distance2(endLocalisations[i]);
        if (d2 <= minD) {
          minD = d2;
          minI = i;
        }
      }

      if (minI == -1) {
        return 0;
      }

      return endLocalisations[minI].trace;
    }
    return 0;
  }

  /**
   * Find the earliest forerunner spot (from pastIndex to currentIndex) that is within the distance
   * threshold of the given spot. In the event that multiple forerunner spots from the same frame
   * are within the distance, assign the closest spot.
   *
   * <p>This method respects the exclusion distance. No spot can be assigned if a the next closest
   * spot is within the exclusion distance.
   *
   * @param index The index of the spot
   * @param pastIndex The index of the earliest forerunner spot
   * @param currentIndex The index of the first spot in the same frame (i.e. end of forerunner
   *        spots)
   * @return The assigned trace
   */
  private int findForerunnerWithExclusion(final int index, final int pastIndex,
      final int currentIndex) {
    final Localisation spot = localisations[index];
    // Check that the next farthest spot is above the exclusion distance
    float nextMinD = Float.POSITIVE_INFINITY;
    int currentT;
    if (traceMode == TraceMode.EARLIEST_FORERUNNER) {
      currentT = endLocalisations[pastIndex].time;
      for (int i = pastIndex; i < currentIndex; i++) {
        final float d2 = spot.distance2(endLocalisations[i]);
        if (d2 <= distanceThreshSquared) {
          minD = d2;
          int trace = endLocalisations[i].trace;

          // Search all remaining spots that end in this time frame and pick the closest
          final int nextIndex = endIndexes[endLocalisations[i].endTime + 1];
          for (int ii = i + 1; ii < nextIndex; ii++) {
            final float dd2 = spot.distance2(endLocalisations[ii]);
            if (dd2 < minD) {
              nextMinD = minD;
              minD = dd2;
              trace = endLocalisations[ii].trace;
            }
          }

          return (nextMinD > distanceExclusionSquared) ? trace : 0;
        } else if (currentT == endLocalisations[i].time) {
          // Store the minimum distance to the next spot in the same frame
          if (d2 < nextMinD) {
            nextMinD = d2;
          }
        } else {
          // New time frame so reset the distance to the next spot in the same frame
          nextMinD = d2;
        }
        currentT = endLocalisations[i].time;
      }
    } else if (traceMode == TraceMode.LATEST_FORERUNNER) {
      currentT = endLocalisations[currentIndex].time;
      for (int i = currentIndex; i-- > pastIndex;) {
        final float d2 = spot.distance2(endLocalisations[i]);
        if (d2 <= distanceThreshSquared) {
          minD = d2;
          int trace = endLocalisations[i].trace;

          // Search all remaining spots in this time frame and pick the closest
          final int previousIndex = endIndexes[endLocalisations[i].endTime];
          // int previousIndex = i;
          //// Look for the index for the previous time-frame
          // while (previousIndex > 0 && endLocalisations[previousIndex-1].t ==
          // endLocalisations[i].t)
          // previousIndex--;
          // if (previousIndex != endIndex[endLocalisations[i].endT])
          // {
          // System.out.printf("Error when tracing: %d != %d\n", previousIndex,
          // endIndex[endLocalisations[i].endT]);
          // }
          for (int ii = i; ii-- > previousIndex;) {
            final float dd2 = spot.distance2(endLocalisations[ii]);
            if (dd2 < minD) {
              nextMinD = minD;
              minD = dd2;
              trace = endLocalisations[ii].trace;
            }
          }

          return (nextMinD > distanceExclusionSquared) ? trace : 0;
        } else if (currentT == endLocalisations[i].time) {
          // Store the minimum distance to the next spot in the same frame
          if (d2 < nextMinD) {
            nextMinD = d2;
          }
        } else {
          // New time frame so reset the distance to the next spot in the same frame
          nextMinD = d2;
        }
        currentT = endLocalisations[i].time;
      }
    } else {
      // traceMode == TraceMode.SINGLE_LINKAGE

      // Find the closest spot
      minD = distanceThreshSquared;
      int minI = -1;
      for (int i = pastIndex; i < currentIndex; i++) {
        final float d2 = spot.distance2(endLocalisations[i]);
        if (d2 <= minD) {
          minD = d2;
          minI = i;
        }
      }

      if (minI == -1) {
        return 0;
      }

      if (distanceExclusionSquared > 0) {
        // Check all spots in the same frame
        final int previousIndex = endIndexes[endLocalisations[minI].endTime];
        final int nextIndex = endIndexes[endLocalisations[minI].endTime + 1];

        for (int i = previousIndex; i < nextIndex; i++) {
          if (i == minI) {
            continue;
          }
          final float d2 = spot.distance2(endLocalisations[i]);
          if (d2 <= nextMinD) {
            nextMinD = d2;
          }
        }
      }

      return (nextMinD > distanceExclusionSquared) ? endLocalisations[minI].trace : 0;
    }
    return 0;
  }

  /**
   * Find the earliest forerunner spot (from pastIndex to currentIndex) that is within the distance
   * threshold of the given spot. In the event that multiple forerunner spots from the same frame
   * are within the distance, assign the closest spot.
   *
   * <p>Do not assigned to the specified trace to ignore.
   *
   * @param index The index of the spot
   * @param pastIndex The index of the earliest forerunner spot
   * @param currentIndex The index of the first spot in the same frame (i.e. end of forerunner
   *        spots)
   * @param ignoreCount The count of traces to ignore
   * @param ignore The traces to ignore
   * @return The assigned trace
   */
  private int findAlternativeForerunner(final int index, final int pastIndex,
      final int currentIndex, final int ignoreCount, final int[] ignore) {
    if (distanceExclusionSquared == 0) {
      return findAlternativeForerunnerNoExclusion(index, pastIndex, currentIndex, ignoreCount,
          ignore);
    }
    return findAlternativeForerunnerWithExclusion(index, pastIndex, currentIndex, ignoreCount,
        ignore);
  }

  /**
   * Find the earliest forerunner spot (from pastIndex to currentIndex) that is within the distance
   * threshold of the given spot. In the event that multiple forerunner spots from the same frame
   * are within the distance, assign the closest spot.
   *
   * <p>Do not assigned to the specified trace to ignore.
   *
   * @param index The index of the spot
   * @param pastIndex The index of the earliest forerunner spot
   * @param currentIndex The index of the first spot in the same frame (i.e. end of forerunner
   *        spots)
   * @param ignoreCount The count of traces to ignore
   * @param ignore The traces to ignore
   * @return The assigned trace
   */
  private int findAlternativeForerunnerNoExclusion(final int index, final int pastIndex,
      final int currentIndex, final int ignoreCount, final int[] ignore) {
    final Localisation spot = localisations[index];

    if (traceMode == TraceMode.EARLIEST_FORERUNNER) {
      for (int i = pastIndex; i < currentIndex; i++) {
        if (ignore(i, ignoreCount, ignore)) {
          continue;
        }

        final float d2 = spot.distance2(endLocalisations[i]);
        if (d2 <= distanceThreshSquared) {
          minD = d2;
          int trace = endLocalisations[i].trace;

          // Search all remaining spots in this time frame and pick the closest
          final int nextIndex = endIndexes[endLocalisations[i].endTime + 1];
          // int nextIndex = i;
          // // Look for the index for the next time-frame
          // for (int tt = endLocalisations[i].endT + 1; tt < endIndex.length; tt++)
          // {
          // nextIndex = endIndex[tt];
          // if (nextIndex != i)
          // break;
          // }
          for (int ii = i + 1; ii < nextIndex; ii++) {
            if (ignore(ii, ignoreCount, ignore)) {
              continue;
            }

            final float dd2 = spot.distance2(endLocalisations[ii]);
            if (dd2 < minD) {
              minD = dd2;
              trace = endLocalisations[ii].trace;
            }
          }

          return trace;
        }
      }
    } else if (traceMode == TraceMode.LATEST_FORERUNNER) {
      for (int i = currentIndex; i-- > pastIndex;) {
        if (ignore(i, ignoreCount, ignore)) {
          continue;
        }

        final float d2 = spot.distance2(endLocalisations[i]);
        if (d2 <= distanceThreshSquared) {
          minD = d2;
          int trace = endLocalisations[i].trace;

          // Search all remaining spots in this time frame and pick the closest
          final int previousIndex = endIndexes[endLocalisations[i].endTime];
          // int previousIndex = i;
          //// Look for the index for the previous time-frame
          // while (previousIndex > 0 && endLocalisations[previousIndex-1].t ==
          // endLocalisations[i].t)
          // previousIndex--;
          // if (previousIndex != endIndex[endLocalisations[i].endT])
          // {
          // System.out.printf("Error when tracing: %d != %d\n", previousIndex,
          // endIndex[endLocalisations[i].endT]);
          // }
          for (int ii = i; ii-- > previousIndex;) {
            if (ignore(ii, ignoreCount, ignore)) {
              continue;
            }

            final float dd2 = spot.distance2(endLocalisations[ii]);
            if (dd2 < minD) {
              minD = dd2;
              trace = endLocalisations[ii].trace;
            }
          }

          return trace;
        }
      }
    } else {
      // traceMode == TraceMode.SINGLE_LINKAGE

      // Find the closest spot
      minD = distanceThreshSquared;
      int minI = -1;
      for (int i = pastIndex; i < currentIndex; i++) {
        if (ignore(i, ignoreCount, ignore)) {
          continue;
        }

        final float d2 = spot.distance2(endLocalisations[i]);
        if (d2 <= minD) {
          minD = d2;
          minI = i;
        }
      }

      if (minI == -1) {
        return 0;
      }

      return endLocalisations[minI].trace;
    }
    return 0;
  }

  /**
   * Find the earliest forerunner spot (from pastIndex to currentIndex) that is within the distance
   * threshold of the given spot. In the event that multiple forerunner spots from the same frame
   * are within the distance, assign the closest spot.
   *
   * <p>This method respects the exclusion distance. No spot can be assigned if a the next closest
   * spot is within the exclusion distance.
   *
   * <p>Do not assigned to the specified trace to ignore.
   *
   * @param index The index of the spot
   * @param pastIndex The index of the earliest forerunner spot
   * @param currentIndex The index of the first spot in the same frame (i.e. end of forerunner
   *        spots)
   * @param ignoreCount The count of traces to ignore
   * @param ignore The traces to ignore
   * @return The assigned trace
   */
  private int findAlternativeForerunnerWithExclusion(final int index, final int pastIndex,
      final int currentIndex, final int ignoreCount, final int[] ignore) {
    final Localisation spot = localisations[index];

    // Check that the next farthest spot is above the exclusion distance.
    // Note: It is assumed that the spots to ignore have already been assigned following the
    // exclusion distance rules. So it should be impossible for any ignore spots to be closer than
    // the exclusion distance (otherwise they could not be assigned and ignored).
    float nextMinD = Float.POSITIVE_INFINITY;
    int currentT;

    if (traceMode == TraceMode.EARLIEST_FORERUNNER) {
      currentT = endLocalisations[pastIndex].time;
      for (int i = pastIndex; i < currentIndex; i++) {
        if (ignore(i, ignoreCount, ignore)) {
          continue;
        }

        final float d2 = spot.distance2(endLocalisations[i]);
        if (d2 <= distanceThreshSquared) {
          minD = d2;
          int trace = endLocalisations[i].trace;

          // Search all remaining spots in this time frame and pick the closest
          final int nextIndex = endIndexes[endLocalisations[i].endTime + 1];
          // int nextIndex = i;
          // // Look for the index for the next time-frame
          // for (int tt = endLocalisations[i].endT + 1; tt < endIndex.length; tt++)
          // {
          // nextIndex = endIndex[tt];
          // if (nextIndex != i)
          // break;
          // }
          for (int ii = i + 1; ii < nextIndex; ii++) {
            if (ignore(ii, ignoreCount, ignore)) {
              continue;
            }

            final float dd2 = spot.distance2(endLocalisations[ii]);
            if (dd2 < minD) {
              nextMinD = minD;
              minD = dd2;
              trace = endLocalisations[ii].trace;
            }
          }

          return (nextMinD > distanceExclusionSquared) ? trace : 0;
        } else if (currentT == endLocalisations[i].time) {
          // Store the minimum distance to the next spot in the same frame
          if (d2 < nextMinD) {
            nextMinD = d2;
          }
        } else {
          // New time frame so reset the distance to the next spot in the same frame
          nextMinD = d2;
        }
        currentT = endLocalisations[i].time;
      }
    } else if (traceMode == TraceMode.LATEST_FORERUNNER) {
      currentT = endLocalisations[currentIndex].time;
      for (int i = currentIndex; i-- > pastIndex;) {
        if (ignore(i, ignoreCount, ignore)) {
          continue;
        }

        final float d2 = spot.distance2(endLocalisations[i]);
        if (d2 <= distanceThreshSquared) {
          minD = d2;
          int trace = endLocalisations[i].trace;

          // Search all remaining spots in this time frame and pick the closest
          final int previousIndex = endIndexes[endLocalisations[i].endTime];
          // int previousIndex = i;
          //// Look for the index for the previous time-frame
          // while (previousIndex > 0 && endLocalisations[previousIndex-1].t ==
          // endLocalisations[i].t)
          // previousIndex--;
          // if (previousIndex != endIndex[endLocalisations[i].endT])
          // {
          // System.out.printf("Error when tracing: %d != %d\n", previousIndex,
          // endIndex[endLocalisations[i].endT]);
          // }
          for (int ii = i; ii-- > previousIndex;) {
            if (ignore(ii, ignoreCount, ignore)) {
              continue;
            }

            final float dd2 = spot.distance2(endLocalisations[ii]);
            if (dd2 < minD) {
              nextMinD = minD;
              minD = dd2;
              trace = endLocalisations[ii].trace;
            }
          }

          return (nextMinD > distanceExclusionSquared) ? trace : 0;
        } else if (currentT == endLocalisations[i].time) {
          // Store the minimum distance to the next spot in the same frame
          if (d2 < nextMinD) {
            nextMinD = d2;
          }
        } else {
          // New time frame so reset the distance to the next spot in the same frame
          nextMinD = d2;
        }
        currentT = endLocalisations[i].time;
      }
    } else {
      // traceMode == TraceMode.SINGLE_LINKAGE

      // Find the closest spot
      minD = distanceThreshSquared;
      int minI = -1;
      for (int i = pastIndex; i < currentIndex; i++) {
        if (ignore(i, ignoreCount, ignore)) {
          continue;
        }

        final float d2 = spot.distance2(endLocalisations[i]);
        if (d2 <= minD) {
          minD = d2;
          minI = i;
        }
      }

      if (minI == -1) {
        return 0;
      }

      if (distanceExclusionSquared > 0) {
        // Check all spots in the same frame
        final int previousIndex = endIndexes[endLocalisations[minI].endTime];
        final int nextIndex = endIndexes[endLocalisations[minI].endTime + 1];

        for (int i = previousIndex; i < nextIndex; i++) {
          if (i == minI || ignore(i, ignoreCount, ignore)) {
            continue;
          }

          final float d2 = spot.distance2(endLocalisations[i]);
          if (d2 <= nextMinD) {
            nextMinD = d2;
          }
        }
      }

      return (nextMinD > distanceExclusionSquared) ? endLocalisations[minI].trace : 0;
    }
    return 0;
  }

  /**
   * Check if the localisation at the specified index has a trace ID that matches any in the ignore
   * array.
   *
   * @param index the index
   * @param ignoreCount the ignore count
   * @param ignore the ignore
   * @return true, if successful
   */
  private boolean ignore(int index, int ignoreCount, int[] ignore) {
    for (int j = 0; j < ignoreCount; j++) {
      if (localisations[index].trace == ignore[j]) {
        return true;
      }
    }
    return false;
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
   * Filter the traces that start during an activation frame.
   *
   * @param traces the traces
   * @param activationFrameInterval the interval at which the activation laser is used
   * @return the filtered traces
   */
  public Trace[] filterTraces(Trace[] traces, int activationFrameInterval) {
    final Trace[] newTraces = new Trace[traces.length];
    int count = 0;
    for (final Trace trace : traces) {
      final PeakResult r = trace.getHead();
      if (r != null && (r.getFrame() % activationFrameInterval) == 1) {
        newTraces[count++] = trace;
      }
    }
    return Arrays.copyOf(newTraces, count);
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
      throw new IllegalArgumentException("Distancet and time thresholds must be positive");
    }

    final Trace[] neighbours = new Trace[results.size()];
    final PeakResult[] peakResults = results.toArray();

    final float dThresh2 = (float) (distanceThreshold * distanceThreshold);

    if (tracker != null) {
      tracker.progress(0);
    }

    // Initialise
    int nextIndex = 0;

    // Now process all the frames, comparing them to previous and future frames
    while (nextIndex < localisations.length) {
      if (tracker != null) {
        tracker.progress(nextIndex, localisations.length);
      }

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

    if (tracker != null) {
      tracker.progress(1.0);
    }

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
