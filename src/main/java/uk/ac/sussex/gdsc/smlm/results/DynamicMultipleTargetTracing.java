/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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

import gnu.trove.list.array.TIntArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Consumer;
import java.util.function.ToDoubleBiFunction;
import java.util.stream.Collectors;
import uk.ac.sussex.gdsc.core.data.VisibleForTesting;
import uk.ac.sussex.gdsc.core.data.utils.ConversionException;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.match.Matchings;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.ValidationUtils;
import uk.ac.sussex.gdsc.smlm.data.config.ConfigurationException;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.results.count.FrameCounter;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;

/**
 * Trace localisations through a time stack to identify diffusing molecules.
 *
 * <p>This class implements the reconnection method described in:
 *
 * <blockquote> Sergé, et al (2008) Dynamic multiple-target tracing to probe spatiotemporal
 * cartography of cell membranes. Nature Method 5, 687-694. </blockquote>
 *
 * <p>The reconnection of localisation points to existing trajectories is based on maximising the
 * likelihood of the available connection combinations. The connection of a point to a trajectory
 * uses the following three components:
 *
 * <ul>
 *
 * <li>Probability of position based on laws of diffusion
 *
 * <li>Probability of intensity
 *
 * <li>Probability of blink (disappearance)
 *
 * </ul>
 *
 * <p>Details of the probability model is in appendix 2 in the supplementary materials of Sergé, et
 * al (2008).
 */
public class DynamicMultipleTargetTracing {

  /** Reciprocal of the square root of 2 pi. */
  private static final double ONE_OVER_ROOT_2_PI = 1.0 / Math.sqrt(2 * Math.PI);

  /** The results. */
  private final MemoryPeakResults results;

  /** The frame results. */
  private final List<List<PeakResult>> frameResults;

  /** The mean intensity. */
  private final double meanI;

  /** The standard deviation of the intensity. */
  private final double sdI;

  /**
   * Configuration for Dynamic Multiple Target Tracing (DMMT).
   *
   * <p>This class is immutable.
   */
  public static class DmttConfiguration {
    /** The temporal window. */
    private final int temporalWindow;
    /** The local diffusion weight. */
    private final double localDiffusionWeight;
    /** The diffusion coefficient maximum. */
    private final double diffusionCoefficientMaximum;
    /** The on intensity weight. */
    private final double onIntensityWeight;
    /** The disappearance decay factor. */
    private final double disappearanceDecayFactor;
    /** The disappearance threshold. */
    private final int disappearanceThreshold;
    /** Flag to indicate that the intensity probability model should be disabled. */
    private final boolean disableIntensityModel;

    /**
     * Builder for the Dynamic Multiple Target Tracing (DMMT) configuration.
     */
    public static class Builder {
      /** The temporal window. */
      int temporalWindow = 5;
      /** The local diffusion weight. */
      double localDiffusionWeight = 0.9;
      /** The diffusion coefficient maximum. */
      double diffusionCoefficientMaximum;
      /** The on intensity weight. */
      double onIntensityWeight = 0.5;
      /** The disappearance decay factor. */
      double disappearanceDecayFactor = 5;
      /** The disappearance threshold. */
      int disappearanceThreshold = 15;
      /** Flag to indicate that the intensity probability model should be used. */
      boolean disableIntensityModel;

      /**
       * Create an instance.
       *
       * @param diffusionCoefficientMaximum the diffusion coefficient maximum (um<sup>2</sup>/s)
       * @throws IllegalArgumentException if the value is not strictly positive
       * @see #setDiffusionCoefficientMaximum(double)
       */
      Builder(double diffusionCoefficientMaximum) {
        setDiffusionCoefficientMaximum(diffusionCoefficientMaximum);
      }

      /**
       * Create an instance.
       *
       * @param source the source configuration
       */
      Builder(DmttConfiguration source) {
        this.temporalWindow = source.getTemporalWindow();
        this.localDiffusionWeight = source.getLocalDiffusionWeight();
        this.diffusionCoefficientMaximum = source.getDiffusionCoefficientMaximum();
        this.onIntensityWeight = source.getOnIntensityWeight();
        this.disappearanceDecayFactor = source.getDisappearanceDecayFactor();
        this.disappearanceThreshold = source.getDisappearanceThreshold();
        this.disableIntensityModel = source.isDisableIntensityModel();
      }

      /**
       * Sets the temporal window.
       *
       * <p>This is the window used to compute recent local diffusion and intensity statistics.
       *
       * @param temporalWindow the temporal window
       * @return this object
       * @throws IllegalArgumentException if the value is not {@code > 1}
       */
      public Builder setTemporalWindow(int temporalWindow) {
        ValidationUtils.checkArgument(temporalWindow > 1);
        this.temporalWindow = temporalWindow;
        return this;
      }

      /**
       * Sets the local diffusion weight.
       *
       * <p>This is the weighting given to the local diffusion radius verses the global diffusion
       * radius when computing the probability that a trajectory position can diffusion to a
       * candidate point.
       *
       * @param localDiffusionWeight the local diffusion weight
       * @return the local diffusion weight this object
       * @throws IllegalArgumentException if the value is not in the range {@code [0, 1]}
       */
      public Builder setLocalDiffusionWeight(double localDiffusionWeight) {
        checkWeight("localDiffusionWeight", localDiffusionWeight);
        this.localDiffusionWeight = localDiffusionWeight;
        return this;
      }

      /**
       * Sets the diffusion coefficient maximum.
       *
       * <p>This is the maximum expected diffusion coefficient for molecules. It is used to set an
       * upper radius for the circle of diffusion when connecting trajectories to new points.
       *
       * <p>Units are micrometre<sup>2</sup> / second (um<sup>2</sup>/s). The units will be
       * converted to appropriately match the distance and time calibration in the analysed results.
       *
       * @param diffusionCoefficientMaximum the diffusion coefficient maximum (um<sup>2</sup>/s)
       * @return this object
       * @throws IllegalArgumentException if the value is not strictly positive
       */
      public Builder setDiffusionCoefficientMaximum(double diffusionCoefficientMaximum) {
        ValidationUtils.checkStrictlyPositive(diffusionCoefficientMaximum);
        this.diffusionCoefficientMaximum = diffusionCoefficientMaximum;
        return this;
      }

      /**
       * Sets the on intensity weight.
       *
       * <p>This is the weighting {@code w} given to the probability that a trajectory intensity
       * matches a candidate point. {@code 1 - w} is the weighting given to the probability that the
       * trajectory blinked, i.e. went into the off state.
       *
       * @param onIntensityWeight the on intensity weight
       * @return this object
       * @throws IllegalArgumentException if the value is not in the range {@code [0, 1]}
       */
      public Builder setOnIntensityWeight(double onIntensityWeight) {
        checkWeight("onIntensityWeight", onIntensityWeight);
        this.onIntensityWeight = onIntensityWeight;
        return this;
      }

      /**
       * Sets the disappearance decay factor.
       *
       * <p>The factor {@code r_off} used to compute the exponential decay for the probability of
       * disappearance:
       *
       * <pre>
       * p_off(t) = exp(-(t - t_off) / r_off)
       * </pre>
       *
       * <p>for blinking, previously detected particles (thus in the full off state since time
       * t_off).
       *
       * @param disappearanceDecayFactor the disappearance decay factor
       * @return this object
       * @throws IllegalArgumentException if the value is not strictly positive
       */
      public Builder setDisappearanceDecayFactor(double disappearanceDecayFactor) {
        ValidationUtils.checkStrictlyPositive(disappearanceDecayFactor);
        this.disappearanceDecayFactor = disappearanceDecayFactor;
        return this;
      }

      /**
       * Sets the disappearance threshold.
       *
       * <p>This is the cut-off to end a trajectory that is in the full off state.
       *
       * @param disappearanceThreshold the disappearance threshold
       * @return this object
       * @throws IllegalArgumentException if the value is not positive
       */
      public Builder setDisappearanceThreshold(int disappearanceThreshold) {
        ValidationUtils.checkPositive(disappearanceThreshold);
        this.disappearanceThreshold = disappearanceThreshold;
        return this;
      }

      /**
       * Set if the intensity probability model is disabled.
       *
       * @param disableIntensityModel true if the intensity model is disabled
       * @return this object
       */
      public Builder setDisableIntensityModel(boolean disableIntensityModel) {
        this.disableIntensityModel = disableIntensityModel;
        return this;
      }

      /**
       * Builds the Dynamic Multiple Target Tracing (DMMT) configuration.
       *
       * @return the DMTT configuration
       */
      public DmttConfiguration build() {
        return new DmttConfiguration(this);
      }

      /**
       * Check the weight is in the range {@code [0, 1]}.
       *
       * @param name the name
       * @param weight the weight
       * @throws IllegalArgumentException if the value is not in the range {@code [0, 1]}
       */
      private static void checkWeight(String name, double weight) {
        ValidationUtils.checkArgument(weight >= 0 && weight <= 1,
            () -> name + " not in the range [0, 1]: " + weight);
      }
    }

    /**
     * Create an instance.
     *
     * @param source the source
     */
    DmttConfiguration(Builder source) {
      this.temporalWindow = source.temporalWindow;
      this.localDiffusionWeight = source.localDiffusionWeight;
      this.diffusionCoefficientMaximum = source.diffusionCoefficientMaximum;
      this.onIntensityWeight = source.onIntensityWeight;
      this.disappearanceDecayFactor = source.disappearanceDecayFactor;
      this.disappearanceThreshold = source.disappearanceThreshold;
      this.disableIntensityModel = source.disableIntensityModel;
    }

    /**
     * Create a new builder.
     *
     * @param diffusionCoefficientMaximum the diffusion coefficient maximum (um<sup>2</sup>/s)
     * @return the builder
     * @throws IllegalArgumentException if the value is not strictly positive
     */
    public static Builder newBuilder(double diffusionCoefficientMaximum) {
      return new Builder(diffusionCoefficientMaximum);
    }

    /**
     * Create a new builder from the current configuration.
     *
     * @return the builder
     */
    public Builder toBuilder() {
      return new Builder(this);
    }

    /**
     * Gets the temporal window.
     *
     * <p>This is the window used to compute recent local diffusion and intensity statistics.
     *
     * @return the temporal window
     */
    public int getTemporalWindow() {
      return temporalWindow;
    }

    /**
     * Gets the local diffusion weight.
     *
     * <p>This is the weighting given to the local diffusion radius verses the global diffusion
     * radius when computing the probability that a trajectory position can diffusion to a candidate
     * point.
     *
     * @return the local diffusion weight
     */
    public double getLocalDiffusionWeight() {
      return localDiffusionWeight;
    }

    /**
     * Gets the diffusion coefficient maximum.
     *
     * <p>This is the maximum expected diffusion coefficient for molecules. It is used to set an
     * upper radius for the circle of diffusion when connecting trajectories to new points.
     *
     * <p>Units are micrometre<sup>2</sup> / second (um<sup>2</sup>/s). The units will be converted
     * to appropriately match the distance and time calibration in the analysed results.
     *
     * @return the diffusion coefficient maximum
     */
    public double getDiffusionCoefficientMaximum() {
      return diffusionCoefficientMaximum;
    }

    /**
     * Gets the on intensity weight.
     *
     * <p>This is the weighting {@code w} given to the probability that a trajectory intensity
     * matches a candidate point. {@code 1 - w} is the weighting given to the probability that the
     * trajectory blinked, i.e. went into the off state.
     *
     * @return the on intensity weight
     */
    public double getOnIntensityWeight() {
      return onIntensityWeight;
    }

    /**
     * Gets the disappearance decay factor.
     *
     * <p>The factor {@code r_off} used to compute the exponential decay for the probability of
     * disappearance:
     *
     * <pre>
     * p_off(t) = exp(-(t - t_off) / r_off)
     * </pre>
     *
     * <p>for blinking, previously detected particles (thus in the full off state since time t_off).
     *
     * @return the disappearance decay factor
     */
    public double getDisappearanceDecayFactor() {
      return disappearanceDecayFactor;
    }

    /**
     * Gets the disappearance threshold.
     *
     * <p>This is the cut-off to end a trajectory that is in the full off state.
     *
     * @return the disappearance threshold
     */
    public int getDisappearanceThreshold() {
      return disappearanceThreshold;
    }

    /**
     * Checks if the intensity probability model is disabled.
     *
     * @return true if the intensity model is disabled
     */
    public boolean isDisableIntensityModel() {
      return disableIntensityModel;
    }
  }

  /**
   * Represent a trajectory. Contains the local diffusion coefficient (mean of the squared jump
   * distances) and mean and standard deviation of the intensity over the previous n frames.
   *
   * <p>This class is a state holder. State is not encapsulated and the update logic is performed
   * externally.
   */
  static class Trajectory {
    /** The id. */
    final int id;
    /** The results. */
    final PeakResultStoreList results = new ArrayPeakResultStore(11);
    /** The indices for frames when the trajectory was on. */
    final TIntArrayList onFrames;
    /**
     * The gap between the last frame in the trajectory and the current frame minus 1. A value of
     * zero indicates the current frame is adjacent.
     */
    int gap;
    /**
     * Flag to indicate the trajectory has a local diffusion coefficient. If false then the
     * probability for diffusion will use only the global model.
     */
    boolean isLocalDiffusion;
    /**
     * Flag to indicate the trajectory has a local intensity. If false then the probability for
     * intensity will use only the global model.
     */
    boolean isLocalIntensity;
    /**
     * The local mean squared jump distance (MSD or r^2). Corresponds to the diffusion coefficient
     * (D) multiplied by 4 (MSD = 4D).
     */
    double r2;
    /**
     * The local root mean squared jump distance normalised by the frame gap (sqrt(r^2 * t)).
     * Corresponds to the sqrt of r2 but adjusted using the frame gap to increase the local search
     * radius around the trajectory.
     */
    double rlocal;
    /** The local mean intensity. */
    double meanI;
    /** The local standard deviation of the intensity. */
    double sdI;

    /**
     * Create an instance.
     *
     * @param id the id
     * @param result the result
     * @param on the on flag
     */
    Trajectory(int id, PeakResult result, boolean on) {
      this.id = id;
      onFrames = new TIntArrayList();
      add(result, on);
    }

    /**
     * Create an instance. Specialised version when no on-frames are being tracked.
     *
     * @param id the id
     * @param result the result
     */
    Trajectory(int id, PeakResult result) {
      this.id = id;
      onFrames = null;
      add(result);
    }

    /**
     * Adds the result with a flag indicating if this is classified as an on-frame.
     *
     * @param result the result
     * @param on the on flag
     */
    void add(PeakResult result, boolean on) {
      if (on) {
        onFrames.add(results.size());
      }
      results.add(result);
    }

    /**
     * Adds the result. Specialised version when no on-frames are being tracked.
     *
     * @param result the result
     */
    void add(PeakResult result) {
      results.add(result);
    }

    /**
     * Gets the id.
     *
     * @return the id
     */
    int getId() {
      return id;
    }

    /**
     * Reset the trajectory for the current frame.
     *
     * @param frame the frame
     */
    void reset(int frame) {
      gap = frame - getLast(-1).getFrame() - 1;
    }

    /**
     * Gets the last n-th result from the end of the results identified using a negative index, e.g.
     * use -1 for the last result and -2 for the second to last.
     *
     * <p>Warning: Performs no checks that the results size is greater than or equal to
     * {@code |index|}.
     *
     * @param index the index (must be negative)
     * @return the last n-th result
     */
    PeakResult getLast(int index) {
      return results.get(results.size() + index);
    }

    /**
     * Get the size.
     *
     * @return the size
     */
    int size() {
      return results.size();
    }

    /**
     * Get the size of the on-frames.
     *
     * @return the size
     */
    int onSize() {
      return onFrames.size();
    }

    /**
     * Sets the local diffusion. If the mean squared jump distance is below the maximum then the
     * local diffusion flag is set to true.
     *
     * @param r2 the mean squared jump distance
     * @param r2max the maximum mean squared jump distance
     */
    void setLocalDiffusion(double r2, double r2max) {
      this.r2 = r2;
      // Local diffusion must be below the maximum
      isLocalDiffusion = r2 < r2max;
    }

    /**
     * Sets the local intensity. If the standard deviation is above 0 then the local intensity flag
     * is set to true.
     *
     * @param meanI the mean intensity
     * @param sdI the standard deviation of the intensity
     */
    void setLocalIntensity(double meanI, double sdI) {
      this.meanI = meanI;
      this.sdI = sdI;
      // Local intensity must have a standard deviation
      isLocalIntensity = sdI > 0;
    }

    /**
     * Pass the last specified number of results to the action.
     *
     * <p>Warning: Performs no checks that the results size is greater than or equal to n.
     *
     * @param n the number of results
     * @param action the action
     */
    void forLast(int n, Consumer<PeakResult> action) {
      int count = n;
      for (int i = results.size() - 1; count > 0; i--, count--) {
        action.accept(results.get(i));
      }
    }

    /**
     * Pass the last specified number of on results (i.e. non-blinking) to the action.
     *
     * <p>Warning: Performs no checks that the results size is greater than or equal to n.
     *
     * @param n the number of results
     * @param action the action
     */
    void forLastOn(int n, Consumer<PeakResult> action) {
      int count = n;
      for (int i = onFrames.size() - 1; count > 0; i--, count--) {
        action.accept(results.get(onFrames.getQuick(i)));
      }
    }

    /**
     * Convert the trajectory to a trace.
     *
     * @return the trace
     */
    Trace toTrace() {
      final Trace trace = new Trace(results);
      trace.setId(id);
      return trace;
    }
  }

  /**
   * Create an instance.
   *
   * @param results the results
   * @throws IllegalArgumentException if results are empty or have no exposure time
   * @throws ConversionException if results are not calibrated for distance
   * @throws ConfigurationException if results are not calibrated
   */
  public DynamicMultipleTargetTracing(final MemoryPeakResults results) {
    ValidationUtils.checkNotNull(results, "results");
    ValidationUtils.checkStrictlyPositive(results.size(), "results size");
    // Require conversion of diffusion coefficients specified in um^2/s
    results.getDistanceConverter(DistanceUnit.UM);
    ValidationUtils.checkArgument(results.getCalibrationReader().hasExposureTime(),
        "no exposure time");

    // Order by frame
    this.results = results.copy();
    this.results.sort();

    // Extract each frame
    final LocalList<PeakResult> list = new LocalList<>();
    frameResults = new LocalList<>();
    final FrameCounter counter = new FrameCounter(this.results.getFirstFrame() - 1);

    // Compute mean and SD of the intensity
    final Statistics stats = new Statistics();

    this.results.forEach((PeakResultProcedure) r -> {
      if (counter.advance(r.getFrame())) {
        frameResults.add(list.copy());
        list.clear();
      }
      list.add(r);
      stats.add(r.getIntensity());
    });
    frameResults.add(list.copy());
    // remove the first empty array
    frameResults.remove(0);

    // Allow a fixed intensity (i.e. no standard deviation).
    // This is done using an abstraction of the probability model for intensity
    // so that it returns a constant (p=1.0) if the intensity is fixed. This allows
    // tracing based only on diffusion and blinking.
    meanI = stats.getMean();
    sdI = stats.getStandardDeviation();
  }

  /**
   * Trace localisations across frames that are the same molecule.
   *
   * <p>Note: A trace may contain one or more peak results. All traces have unique identifiers
   * starting from 1.
   *
   * @param configuration the configuration
   * @return the traces
   */
  public List<Trace> traceMolecules(DmttConfiguration configuration) {
    ValidationUtils.checkNotNull(configuration, "configuration");

    // Initialise with the first frame
    final List<Trajectory> activeTrajectories = createTrajectories();
    final List<Trajectory> allTrajectories = new LocalList<>(activeTrajectories);

    // Maximum likelihood reconnection test.
    //
    // Performs a Kuhn-Munkres algorithm for maximising the likelihood of all-vs-all connections.
    // Complexity is O(kkl) where k = min(n, m) and l = max(n, m).
    // The connection test is all-vs-all for complexity O(nm).
    //
    // The likelihood is reversed for a minimum distance matching.
    // We can use double max value for the negative log likelihood threshold since
    // any invalid matching is computed as positive infinity.
    // Note: The actual scores for matches will be re-mapped using the range of valid values
    // by the Matchings class so the threshold value does not effect the scoring.
    // The resulting matchings are consumed directly to either extend a
    // trajectory or create a new one.

    // Create the functions.

    // Compute the negative log likelihood for the connection.
    final ToDoubleBiFunction<Trajectory, PeakResult> edges = createConnectionModel(configuration);

    // Extend each matched trajectory with the result.
    // This could be updated with different models for computing the local statistics.
    final BiConsumer<Trajectory, PeakResult> matched = createMatchedAction(configuration);

    // This could be used to mark the trajectory as entered into an off state.
    final Consumer<Trajectory> unmatchedA = null;

    // For any unmatched result create a new trajectory
    final Consumer<PeakResult> unmatchedB =
        createUnmatchedPeakAction(configuration, allTrajectories);

    // For each remaining frame connect trajectories
    for (int i = 1; i < frameResults.size(); i++) {
      final List<PeakResult> results = frameResults.get(i);
      final int frame = results.get(0).getFrame();

      updateTrajectories(configuration, activeTrajectories, frame);
      final int currentSize = allTrajectories.size();

      // Performs the reconnection test
      Matchings.minimumDistance(activeTrajectories, results, edges, Double.MAX_VALUE, matched,
          unmatchedA, unmatchedB);

      // Copy new trajectories to the currently active trajectories
      activeTrajectories.addAll(allTrajectories.subList(currentSize, allTrajectories.size()));
    }

    // Convert trajectories to traces
    return allTrajectories.stream().map(Trajectory::toTrace).collect(Collectors.toList());
  }

  /**
   * Checks the intensity model should be disabled.
   *
   * @param configuration the configuration
   * @return true to disable the intensity model
   */
  private boolean isDisableIntensityModel(DmttConfiguration configuration) {
    return configuration.isDisableIntensityModel() || sdI == 0;
  }

  /**
   * Creates the connection model. This computes the negative log-likelihood for the connection
   * between a trajectory and a localisation. No connection is positive infinity.
   *
   * @param configuration the configuration
   * @return the connection model
   */
  private ToDoubleBiFunction<Trajectory, PeakResult>
      createConnectionModel(DmttConfiguration configuration) {
    final double dMax = computeDMax(configuration);

    // Create the maximum radius for diffusion with different dark frames
    final double[] rMax = createRMax(dMax, configuration);
    // final double r2max = rMax[0] * rMax[0];
    // Pre-compute the search radius threshold
    // Use r < 3r_max => r^2 < 9 r_max^2.
    final double[] r2maxThreshold = Arrays.stream(rMax).map(d -> d * d * 9).toArray();
    final double[] logPOff = createLogPOff(configuration);

    final double ldw = configuration.getLocalDiffusionWeight();
    final double gdw = 1 - ldw;

    // Optionally disable the intensity model
    if (isDisableIntensityModel(configuration)) {
      return (t, r) -> {
        final double r2 = r.distance2(t.getLast(-1));
        // Must be within the maximum radius for the frame gap.
        if (r2 < r2maxThreshold[t.gap]) {
          // p(reconnection) = p(diffusion) * p(intensity) * p(off)
          // sum for log probability.
          // Use the local or global model as appropriate.
          // The p(off) can be pre-computed.
          double pDiff = probDiffusion(r2, rMax[t.gap]);
          if (t.isLocalDiffusion) {
            // Combined global and local diffusion
            pDiff = gdw * pDiff + ldw * probDiffusion(r2, t.rlocal);
          }
          // Negative log likelihood: p(Diffusion) * p(Off)
          return -Math.log(pDiff) - logPOff[t.gap];
        }

        // no connection probability
        return Double.POSITIVE_INFINITY;
      };
    }

    final double oiw = configuration.getOnIntensityWeight();
    final double biw = 1 - oiw;

    return (t, r) -> {
      final double r2 = r.distance2(t.getLast(-1));
      // Must be within the maximum radius for the frame gap.
      if (r2 < r2maxThreshold[t.gap]) {
        // p(reconnection) = p(diffusion) * p(intensity) * p(off)
        // sum for log probability.
        // Use the local or global model as appropriate.
        // The p(off) can be pre-computed.
        double pDiff = probDiffusion(r2, rMax[t.gap]);
        if (t.isLocalDiffusion) {
          // Combined global and local diffusion
          pDiff = gdw * pDiff + ldw * probDiffusion(r2, t.rlocal);
        }
        final double intensity = r.getIntensity();
        final double pInt =
            t.isLocalIntensity ? oiw * probOn(intensity, t.meanI, t.sdI) + biw * probBlink(t.meanI)
                : oiw * probOn(intensity, meanI, sdI) + biw * probBlink(meanI);
        // Negative log likelihood: p(Diffusion) * p(Intensity) * p(Off)
        return -Math.log(pDiff) - Math.log(pInt) - logPOff[t.gap];
      }

      // no connection probability
      return Double.POSITIVE_INFINITY;
    };
  }

  /**
   * Creates the action to process a connection between a trajectory and a localisation.
   *
   * <p>Currently the temporal window is occurrences and not an actual time span. Thus some local
   * statistics can be averaged over a long duration. However if the trajectory is for a blinking
   * particle with short on-times and long off-times the local statistics are no better or worse
   * than the global model.
   *
   * @param configuration the configuration
   * @return the matched action
   */
  private BiConsumer<Trajectory, PeakResult> createMatchedAction(DmttConfiguration configuration) {
    final double dMax = computeDMax(configuration);
    // D_max = r^2 / 4t ('t' is frame acquisition time)
    // dMax has been converted to frames
    final double r2max = dMax * 4;

    // Optionally disable the intensity model
    if (isDisableIntensityModel(configuration)) {
      return (t, r) -> {
        // Add the peak to the trajectory. Do not track on-frames for intensity.
        t.add(r);

        updateLocalDiffusion(configuration, r2max, t);
        // Ignore local intensity
      };
    }

    final Statistics stats = new Statistics();
    return (t, r) -> {
      // Add the peak to the trajectory tracking on-frames
      final double intensity = r.getIntensity();
      final boolean on =
          t.isLocalIntensity ? isOn(intensity, t.meanI, t.sdI) : isOn(intensity, meanI, sdI);
      t.add(r, on);

      updateLocalDiffusion(configuration, r2max, t);
      updateLocalIntensity(configuration, stats, t);
    };
  }

  /**
   * Update local diffusion. The maximum jump distance is used to switch global to local diffusion
   * (which must be slower than the global diffusion).
   *
   * @param configuration the configuration
   * @param r2max the maximum mean squared jump distance
   * @param t the trajectory
   */
  private static void updateLocalDiffusion(DmttConfiguration configuration, final double r2max,
      Trajectory t) {
    if (t.size() >= configuration.getTemporalWindow()) {
      PeakResult after = t.getLast(-1);
      double sum = 0;
      for (int j = 2; j <= configuration.getTemporalWindow(); j++) {
        final PeakResult before = t.getLast(-j);
        sum += after.distance2(before) / (after.getFrame() - before.getFrame());
        after = before;
      }
      t.setLocalDiffusion(sum / (configuration.getTemporalWindow() - 1), r2max);
    }
  }

  /**
   * Update local intensity.
   *
   * @param configuration the configuration
   * @param stats the working statistics instance
   * @param t the trajectory
   */
  private static void updateLocalIntensity(DmttConfiguration configuration, final Statistics stats,
      Trajectory t) {
    if (t.onSize() >= configuration.getTemporalWindow()) {
      stats.reset();
      t.forLastOn(configuration.getTemporalWindow(), peak -> stats.add(peak.getIntensity()));
      t.setLocalIntensity(stats.getMean(), stats.getStandardDeviation());
    }
  }

  /**
   * Creates the the action for an unmatched localisation. This will create a new trajectory and add
   * it to the list.
   *
   * @param configuration the configuration
   * @param allTrajectories the list of all trajectories
   * @return the action
   */
  private Consumer<PeakResult> createUnmatchedPeakAction(DmttConfiguration configuration,
      List<Trajectory> allTrajectories) {
    // Optionally disable the intensity model
    if (isDisableIntensityModel(configuration)) {
      return r -> {
        // Do not track on-frames for intensity.
        allTrajectories.add(new Trajectory(allTrajectories.size() + 1, r));
      };
    }
    return r -> {
      final double intensity = r.getIntensity();
      final boolean on = isOn(intensity, meanI, sdI);
      allTrajectories.add(new Trajectory(allTrajectories.size() + 1, r, on));
    };
  }

  /**
   * Compute the maximum diffusion coefficient (D_max) in the units of the data:
   * {@code um^2/s -> units^2/frame} (where units are the results native distance units).
   *
   * @param configuration the configuration
   * @return D_max
   */
  private double computeDMax(DmttConfiguration configuration) {
    // Convert the diffusion coefficient to the units of the data
    // um^2/s -> units^2/frame (where units are the results native distance units)
    final TypeConverter<DistanceUnit> dc = results.getDistanceConverter(DistanceUnit.UM);
    // um^2/s -> units^2/s
    double dMax = dc.convertBack(dc.convertBack(configuration.getDiffusionCoefficientMaximum()));
    // units^2/s -> units^2/frame
    final double exposureTimeSecondsPerFrame =
        results.getCalibrationReader().getExposureTime() / 1000;
    dMax = dMax * exposureTimeSecondsPerFrame;
    return dMax;
  }

  /**
   * Creates the maximum radius for diffusion (R_max) when the trajectory was last observed
   * {@code n+1} frames in the past. Thus {@code R_max[0]} is the previous frame, {@code R_max[1]}
   * is the frame before that all the way up to
   * {@link DmttConfiguration#getDisappearanceThreshold()}.
   *
   * <p>R_max is equivalent to the root mean squared distance expected for the maximum diffusion.
   * Typically 3 times this level sets an upper limit on the distance for diffusion assuming
   * Brownian motion is a Gaussian random walk with standard deviation R_max.
   *
   * @param dMax the maximum diffusion coefficient (units^2/frame)
   * @param configuration the configuration
   * @return R_max
   */
  private static double[] createRMax(double dMax, DmttConfiguration configuration) {
    final int n = configuration.getDisappearanceThreshold();
    final double[] rMax = new double[n + 1];
    // D_max = r^2 / 4t ('t' is frame acquisition time)
    // dMax has been converted to frames
    final double r2 = dMax * 4;
    // Increase maximum radius using:
    // R^2 = (3r_max)^2 * t_off = 9 * (r_max)^2 * t
    // (with t_off being the number of elapsed frames since the particle turned off)
    // The value 3 is based on Brownian diffusion and the 99% range of a Gaussian for
    // the Mean Squared Distance. It could be a configuration parameter to limit the
    // range for diffusion. The same effect is achieved by just increasing D_max although
    // the relationship is not linear.
    // Here we omit the factor 3 thus we are scaling the r_max^2 value by the frame gap.
    for (int i = 0; i <= n; i++) {
      rMax[i] = Math.sqrt(r2 * (i + 1));
    }
    return rMax;
  }

  /**
   * Creates the log of the probability for blink/disappearance (p(off)).
   *
   * @param configuration the configuration
   * @return log(p(off))
   */
  private static double[] createLogPOff(DmttConfiguration configuration) {
    final int n = configuration.getDisappearanceThreshold();
    final double[] logPOff = new double[n + 1];
    // p_off(t) = exp(-(t - t_off) / r_off)
    final double r_off = configuration.getDisappearanceDecayFactor();
    // Here t is the gap in frames between the current frame t and the frame when the trajectory
    // entered the full off state (t_off).
    // t=0 indicates no gap and the the trajectory was on in the last frame. This has a probability
    // of 1.0 and so log(p) of zero.
    for (int t = 1; t <= n; t++) {
      logPOff[t] = -t / r_off;
    }
    return logPOff;
  }

  /**
   * Creates the trajectories from the first frame.
   *
   * @return the trajectories
   */
  private List<Trajectory> createTrajectories() {
    final LocalList<Trajectory> trajectories = new LocalList<>();
    for (final PeakResult result : frameResults.get(0)) {
      // Compute if on.
      final boolean on = isOn(result.getIntensity(), meanI, sdI);
      trajectories.add(new Trajectory(trajectories.size() + 1, result, on));
    }
    return trajectories;
  }

  /**
   * Checks if the molecule is on. Assume it is on if the intensity is above the mean intensity.
   * Otherwise the on probability must be above the blink probability.
   *
   * @param intensity the intensity
   * @param meanI the mean intensity
   * @param sdI the standard deviation of the intensity
   * @return true if is on
   */
  @VisibleForTesting
  static boolean isOn(double intensity, double meanI, double sdI) {
    return (intensity > meanI) || probOn(intensity, meanI, sdI) > probBlink(meanI);
  }

  /**
   * Update the trajectories for the new frame. Removes old trajectories. Updates the local
   * diffusion radius for the frame gap.
   *
   * @param configuration the configuration
   * @param traj the trajectories
   * @param frame the frame
   */
  private static void updateTrajectories(DmttConfiguration configuration, List<Trajectory> traj,
      final int frame) {
    // Remove all old trajectories
    traj.removeIf(r -> {
      r.reset(frame);
      if (r.gap > configuration.getDisappearanceThreshold()) {
        return true;
      }
      // Update the local diffusion radius for the frame gap
      r.rlocal = Math.sqrt(r.r2 * (r.gap + 1));
      return false;
    });
  }

  /**
   * Compute the diffusion probability (p(diff)). This is a Gaussian distributed probability
   * function.
   *
   * @param r2 the squared distance between the points
   * @param rLocal the local root mean squared radius expected for diffusion.
   * @return p(diff)
   */
  private static double probDiffusion(double r2, double rLocal) {
    return ONE_OVER_ROOT_2_PI * Math.exp(-r2 / (2 * rLocal * rLocal));
  }

  /**
   * Compute the on probability (p(on)). This is a Gaussian distributed probability function.
   *
   * @param i the intensity
   * @param mean the mean
   * @param s the standard deviation
   * @return p(on)
   */
  private static double probOn(double i, double mean, double s) {
    final double di = i - mean;
    return ONE_OVER_ROOT_2_PI * Math.exp(-di * di / (2 * s * s));
  }

  /**
   * Compute the blink probability (p(blink)). This is a uniformly distributed probability function.
   *
   * @param mean the mean
   * @return p(blink)
   */
  private static double probBlink(double mean) {
    return 1.0 / mean;
  }
}
