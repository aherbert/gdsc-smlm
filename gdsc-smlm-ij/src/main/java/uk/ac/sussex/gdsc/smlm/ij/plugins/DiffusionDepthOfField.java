/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Package
 *
 * Software for single molecule localisation microscopy (SMLM) in ImageJ
 * %%
 * Copyright (C) 2011 - 2025 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.ij.plugins;

import ij.IJ;
import ij.Prefs;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.process.LUT;
import ij.text.TextWindow;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import java.awt.Color;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.InputMismatchException;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.Scanner;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.DoublePredicate;
import java.util.function.DoubleUnaryOperator;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.NonMonotonicSequenceException;
import org.apache.commons.math3.exception.NumberIsTooSmallException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.exception.TooManyIterationsException;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.OptimizationData;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.SimpleValueChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.rng.JumpableUniformRandomProvider;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.ContinuousSampler;
import org.apache.commons.rng.sampling.distribution.NormalizedGaussianSampler;
import uk.ac.sussex.gdsc.core.ij.ImageJPluginLoggerHelper;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper.LutColour;
import uk.ac.sussex.gdsc.core.logging.LoggerUtils;
import uk.ac.sussex.gdsc.core.logging.Ticker;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.rng.RandomGeneratorAdapter;
import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;
import uk.ac.sussex.gdsc.core.utils.rng.UniformRandomProviders;
import uk.ac.sussex.gdsc.smlm.fitting.DiffusionAnalysis;
import uk.ac.sussex.gdsc.smlm.math3.optim.nonlinear.scalar.noderiv.CustomPowellOptimizer;

/**
 * Move a set of molecules within a depth-of-field (DoF) to generate a distribution of the
 * probability of remaining within the DoF.
 *
 * <p>This empirical distribution can be used to fit correction factors {@code a} and {@code b} used
 * scale the depth of field {@code dz} to {@code dz_corr} for a given diffusion coefficient
 * {@code D}. The scaling factors depend on the DoF {@code dz}, the diffusion time step {@code dt},
 * and the gap size {@code g} allowed in tracks.
 *
 * <pre>
 * dz_corr = dz + a * sqrt(D) + b
 * </pre>
 *
 * <p>The corrected depth of field can be used to compute the probability that a diffusing molecule
 * remains within the depth-of-field after a given time. This is based on the Spot-On model
 * described in the methods section of the paper:
 *
 * <p>Hansen, A.S., Woringer, M., Grimm, J.B., Lavis, L.D., Tjian, R., and Darzacq, X. (2018) Robust
 * model-based analysis of single-particle tracking experiments with Spot-On. eLife 7, e33125.
 * doi:10.7554/eLife.33125.
 */
public class DiffusionDepthOfField implements PlugIn {
  private static final String TITLE = "Diffusion Depth of Field";
  private static final String[] OPTIMISER_MODES = {"Powell", "CMA-ES", "BOBYQA"};
  private static final int OPT_MODE_CMAES = 1;
  private static final int OPT_MODE_BOBYQA = 2;

  private static final AtomicReference<TextWindow> TABLE_REF = new AtomicReference<>();

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> INSTANCE = new AtomicReference<>(new Settings());

    double depthOfField;
    double exposureTime;
    int gap;

    int numberOfMolecules;
    int maxT;
    boolean useDetector;
    String detectorCurve;
    boolean allowRestarts;
    double halfLife;
    double minD;
    double maxD;
    int sampleD;

    int optimiserMode;
    int repeats;
    double minA;
    double maxA;
    double minB;
    double maxB;

    Settings() {
      depthOfField = 750;
      exposureTime = 10;
      gap = 1;
      // SpotOn paper used: n=50000, maxT=15; D in [1, 12]
      // Using a lower min D anchors the fit to slower moving molecules.
      numberOfMolecules = 100000;
      maxT = 15;
      minD = 0.2;
      maxD = 12;
      sampleD = 20;
      // Fitting
      optimiserMode = OPT_MODE_CMAES;
      repeats = 3;
      maxA = 0.3;
      maxB = 0.3;
    }

    Settings(Settings source) {
      depthOfField = source.depthOfField;
      exposureTime = source.exposureTime;
      gap = source.gap;

      numberOfMolecules = source.numberOfMolecules;
      maxT = source.maxT;
      useDetector = source.useDetector;
      detectorCurve = source.detectorCurve;
      halfLife = source.halfLife;
      allowRestarts = source.allowRestarts;
      minD = source.minD;
      maxD = source.maxD;
      sampleD = source.sampleD;

      optimiserMode = source.optimiserMode;
      repeats = source.repeats;
      minA = source.minA;
      maxA = source.maxA;
      minB = source.minB;
      maxB = source.maxB;
    }

    Settings copy() {
      return new Settings(this);
    }

    /**
     * Load a copy of the settings.
     *
     * @return the settings
     */
    static Settings load() {
      return INSTANCE.get().copy();
    }

    /**
     * Save the settings. This can be called only once as it saves via a reference.
     */
    void save() {
      INSTANCE.set(this);
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    settings = Settings.load();
    run(settings);
  }

  double[] run(double dz) {
    settings = Settings.load();
    if (dz > 0) {
      settings.depthOfField = dz;
    }
    return run(settings);
  }

  private double[] run(Settings settings) {
    settings.save();

    if (!showDialog()) {
      return null;
    }
    final int threadCount = Prefs.getThreads();
    final ExecutorService executor = Executors.newFixedThreadPool(threadCount);

    double[] ab = null;
    try {
      final double[][] probability = simulateRemaining(threadCount, executor);

      // Fit the observed probability
      ab = fitRemaining(probability, threadCount, executor);

      addToResultTable(ab);

      final double[][] fittedProbability = createProbability(ab);

      plotRemaining(probability, fittedProbability);

    } finally {
      executor.shutdown();
    }

    ImageJUtils.finished(TITLE + " done");
    return ab;
  }

  private boolean showDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

    gd.addNumericField("Depth_of_field", settings.depthOfField, 1, 6, "nm");
    gd.addNumericField("Exposure_time", settings.exposureTime, 1, 6, "ms");
    gd.addSlider("Gap", 1, 5, settings.gap);

    gd.addNumericField("Number_of_molecules", settings.numberOfMolecules);
    gd.addSlider("Max_t", 5, 15, settings.maxT);
    gd.addCheckbox("Use_dectector", settings.useDetector);
    gd.addFilenameField("Detector_curve", settings.detectorCurve);
    gd.addCheckbox("Allow_restarts", settings.allowRestarts);
    gd.addNumericField("Half_life", settings.halfLife);
    gd.addNumericField("Min_D", settings.minD, 3, 6, "um^2/s");
    gd.addNumericField("Max_D", settings.maxD, 3, 6, "um^2/s");
    gd.addSlider("Sample_D", 2, 30, settings.sampleD);

    gd.addChoice("Optimiser_mode", OPTIMISER_MODES, settings.optimiserMode);
    gd.addSlider("Repeats", 0, 5, settings.repeats);
    gd.addNumericField("Min_A", settings.minA, 3);
    gd.addNumericField("Max_A", settings.maxA, 3);
    gd.addNumericField("Min_B", settings.minB, 3);
    gd.addNumericField("Max_B", settings.maxB, 3);

    gd.addHelp(HelpUrls.getUrl("diffusion-depth-of-field"));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }

    settings.depthOfField = gd.getNextNumber();
    settings.exposureTime = gd.getNextNumber();
    settings.gap = (int) gd.getNextNumber();

    settings.numberOfMolecules = (int) gd.getNextNumber();
    settings.maxT = (int) gd.getNextNumber();
    settings.useDetector = gd.getNextBoolean();
    settings.detectorCurve = gd.getNextString();
    settings.allowRestarts = gd.getNextBoolean();
    settings.halfLife = gd.getNextNumber();
    settings.minD = gd.getNextNumber();
    settings.maxD = gd.getNextNumber();
    settings.sampleD = (int) gd.getNextNumber();

    settings.optimiserMode = gd.getNextChoiceIndex();
    settings.repeats = (int) gd.getNextNumber();
    settings.minA = gd.getNextNumber();
    settings.maxA = gd.getNextNumber();
    settings.minB = gd.getNextNumber();
    settings.maxB = gd.getNextNumber();

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Depth of Field", settings.depthOfField);
      ParameterUtils.isAboveZero("Exposure time", settings.exposureTime);
      ParameterUtils.isAboveZero("Gap", settings.gap);

      ParameterUtils.isAboveZero("Number of molecules", settings.numberOfMolecules);
      ParameterUtils.isAboveZero("Max T", settings.maxT);
      ParameterUtils.isAboveZero("Min D", settings.minD);
      ParameterUtils.isAboveZero("Max D", settings.maxD);
      ParameterUtils.isAbove("Sample D", settings.sampleD, 1);

      ParameterUtils.isPositive("Repeats", settings.repeats);
      ParameterUtils.isPositive("Min A", settings.minA);
      ParameterUtils.isAboveZero("Max A", settings.maxA);
      ParameterUtils.isPositive("Min B", settings.minB);
      ParameterUtils.isAboveZero("Max B", settings.maxB);

      ParameterUtils.isEqualOrAbove("Max D", settings.maxD, settings.minD);
      ParameterUtils.isEqualOrAbove("Max A", settings.maxA, settings.minA);
      ParameterUtils.isEqualOrAbove("Max B", settings.maxB, settings.minB);

      if (settings.allowRestarts) {
        ParameterUtils.isAboveZero("Half-life with restarts", settings.halfLife);
      }
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    if (gd.invalidNumber()) {
      return false;
    }

    return true;
  }

  /**
   * Simulate the probability of remaining in the depth-of-field.
   *
   * <p>Returns a probability for each diffusion coefficient over time.
   *
   * @param threadCount the thread count
   * @param executor the executor
   * @return the probability
   */
  private double[][] simulateRemaining(int threadCount, ExecutorService executor) {
    IJ.showStatus("Simulating tracks...");

    final double halfDz = settings.depthOfField / 2000;
    final double dt = settings.exposureTime / 1000;
    final int g = settings.gap;
    final int maxT = settings.maxT;
    final boolean allowRestarts = settings.allowRestarts;
    final double halfLife = settings.halfLife;

    // Create the depth of field detector curve
    final DoubleUnaryOperator detectorCurve =
        settings.useDetector ? loadDetectorCurve(settings.detectorCurve) : null;

    // Simulate tracks across the depth of field.
    // Each track is scaled using the diffusion coefficient and the molecule tested if
    // it remains in the depth-of-field.

    final double[] sampleD = createDiffusionCoefficients();
    // Create the Gaussian step for each diffusion coefficient
    final double[] step = Arrays.stream(sampleD).map(
        // Std.dev of Gaussian for the step size: sqrt(2D * dt)
        d -> Math.sqrt(2 * d * dt)).toArray();

    final List<Future<int[][]>> futures = new LinkedList<>();

    final int total = settings.numberOfMolecules;
    final Ticker ticker = ImageJUtils.createTicker(total, threadCount);

    final AtomicInteger position = new AtomicInteger(total);
    final int[] numberOfMolecules = new int[step.length];
    final JumpableUniformRandomProvider rng = UniformRandomProviders.createJumpable();
    for (int n = 0; n < threadCount; n++) {
      final NormalizedGaussianSampler sampler =
          SamplerUtils.createNormalizedGaussianSampler(rng.jump());
      // If no restarts then there is no requirement for a half-life.
      // The half-life parameter can be zero and the track will reach max T.
      final ContinuousSampler lifeSampler =
          halfLife <= 0 ? () -> maxT : SamplerUtils.createExponentialSampler(rng.jump(), halfLife);
      // Use an optional detector curve
      final DoublePredicate dof;
      if (settings.useDetector) {
        dof = createDetector(detectorCurve, rng.jump());
      } else {
        dof = z -> Math.abs(z) <= halfDz;
      }
      futures.add(executor.submit(() -> {
        final int[] tracks = new int[step.length];
        // Count of number of molecules inside the depth of field
        // for each diffusion coefficient and time frame
        final int[][] observed = new int[step.length][maxT];
        for (;;) {
          final int pos = position.decrementAndGet();
          if (pos < 0) {
            break;
          }
          // Simulate across the depth of field [-z/2, z/2]
          final double p = pos / (total - 1.0);
          final double originZ = (1 - p) * halfDz - p * halfDz;
          // for each diffusion rate test: |z| < z/2
          for (int j = 0; j < step.length; j++) {
            final double s = step[j];
            final int[] obs = observed[j];
            double z = originZ;
            if (!dof.test(z)) {
              // Not detected
              break;
            }
            // Create a lifetime within the depth-of-field.
            // Round lifetime to nearest integer using +0.5.
            final int lifetime = (int) Math.min(0.5 + lifeSampler.sample(), 20 * maxT);
            // Last frame the molecule was inside the DoF
            int lastT = 0;
            // Origin frame inside the DoF
            // Add 1 to the origin so that t - originT = steps - 1
            int originT = 1;
            tracks[j]++;
            // Diffuse until lifetime expires
            for (int t = 1; t <= lifetime; t++) {
              z += sampler.sample() * s;
              if (dof.test(z)) {
                if (originT < 0) {
                  // Start a track
                  // Add 1 to the origin so that t - originT = steps - 1
                  originT = t + 1;
                  tracks[j]++;
                } else if (t - originT < maxT) {
                  // Record the track
                  obs[t - originT]++;
                }
                lastT = t;
              } else if (t - lastT >= g) {
                // Above the configured gap size; molecule has been lost.
                if (allowRestarts) {
                  // Allow restarting a new track. Signal that the track is unstarted using -1.
                  // This can increase the activations around the edge of the
                  // depth of field.
                  originT = -1;
                } else {
                  break;
                }
              }
            }
          }
          ticker.tick();
        }
        // Final track counts
        synchronized (numberOfMolecules) {
          for (int j = 0; j < step.length; j++) {
            numberOfMolecules[j] += tracks[j];
          }
        }
        return observed;
      }));
    }

    // Collate results
    final double[][] probability = new double[step.length][maxT];
    for (final Future<int[][]> f : futures) {
      try {
        final int[][] observed = f.get();
        for (int j = 0; j < observed.length; j++) {
          for (int i = 0; i < maxT; i++) {
            probability[j][i] += observed[j][i];
          }
        }
      } catch (InterruptedException | ExecutionException e) {
        throw new RuntimeException(e);
      }
    }
    for (int j = 0; j < probability.length; j++) {
      final double totalTracks = numberOfMolecules[j];
      for (int i = 0; i < maxT; i++) {
        probability[j][i] /= totalTracks;
      }
    }

    return probability;
  }

  /**
   * Load a depth-of-field detector curve. The file should contains entries of: z (nm), probability.
   * The returned detector will by calibrated in micrometers.
   *
   * @param detectionCurve the detection curve
   * @return the detector curve
   * @throws UncheckedIOException if the curve cannot be read
   * @throws IllegalArgumentException if the curve is invalid
   */
  static DoubleUnaryOperator loadDetectorCurve(String detectionCurve) {
    try {
      final DoubleArrayList depth = new DoubleArrayList();
      final DoubleArrayList p = new DoubleArrayList();
      final Pattern delim = Pattern.compile("[\t ,]");
      Files.lines(Paths.get(detectionCurve)).forEach(line -> {
        try (Scanner scanner = new Scanner(line)) {
          scanner.useDelimiter(delim);
          scanner.useLocale(Locale.ROOT);
          final double pos = scanner.nextDouble();
          final double prob = scanner.nextDouble();
          if (!(prob >= 0 && prob <= 1)) {
            throw new IllegalArgumentException("Invalid probability: " + prob);
          }
          depth.add(pos * 1e-3);
          p.add(prob);
        } catch (final InputMismatchException e) {
          // Ignore header
          if (!depth.isEmpty()) {
            throw new IllegalArgumentException("Unrecognised curve entry: " + line);
          }
        }
      });
      try {
        final PolynomialSplineFunction fun =
            new SplineInterpolator().interpolate(depth.toDoubleArray(), p.toDoubleArray());
        final double minZ = depth.getDouble(0);
        final double maxZ = depth.getDouble(depth.size() - 1);
        return z -> {
          if (z < minZ | z > maxZ) {
            return 0;
          }
          return fun.value(z);
        };
      } catch (NonMonotonicSequenceException | NumberIsTooSmallException e) {
        throw new IllegalArgumentException("Invalid detector curve", e);
      }
    } catch (final IOException e) {
      throw new UncheckedIOException(e);
    }
  }

  /**
   * Creates the detector.
   *
   * @param detectorCurve the detector curve
   * @param rng the source of randomness
   * @return the detector
   */
  static DoublePredicate createDetector(DoubleUnaryOperator detectorCurve,
      UniformRandomProvider rng) {
    return z -> {
      final double p = detectorCurve.applyAsDouble(z);
      if (p <= 0) {
        return false;
      }
      return rng.nextDouble() < p;
    };
  }

  /**
   * Creates the diffusion coefficients.
   *
   * <p>Returns the configured number of samples log-uniformly distributed from min D to max D.
   *
   * @return the diffusion coefficients
   */
  private double[] createDiffusionCoefficients() {
    // Sample diffusion coefficients
    // Use logarithmic scale
    final double lnMin = Math.log(settings.minD);
    final double lnMax = Math.log(settings.maxD);
    final double[] sampleD = IntStream.rangeClosed(1, settings.sampleD).mapToDouble(i -> {
      // interpolate D in [min, max] using [1, n] as [0, 1]
      final double p = (i - 1.0) / (settings.sampleD - 1);
      return Math.exp((1 - p) * lnMin + p * lnMax);
    }).toArray();
    return sampleD;
  }

  /**
   * Plot the probability of remaining in the depth of field.
   *
   * @param probability the simulated probability
   * @param fitted the fitted probability
   */
  private void plotRemaining(final double[][] probability, double[][] fitted) {
    final double dt = settings.exposureTime / 1000;
    final int maxT = settings.maxT;

    // Plot the observed distribution.
    // For each D add a line for all time steps
    final float[] time = SimpleArrayUtils.newArray(maxT, 1.0f, 1.0f);
    SimpleArrayUtils.apply(time, x -> (float) (x * dt));
    final String title = TITLE + " probability";
    final Plot plot = new Plot(title, "Time (seconds)", "Probability");
    LUT lut = LutHelper.createLut(LutColour.RED);
    for (int i = 0; i < probability.length; i++) {
      plot.setColor(LutHelper.getColour(lut, 2 * (i + 1), 1, 2 * probability.length));
      plot.addPoints(time, SimpleArrayUtils.toFloat(probability[i]), Plot.LINE);
    }
    lut = LutHelper.createLut(LutColour.BLUE);
    for (int i = 0; i < fitted.length; i++) {
      plot.setColor(LutHelper.getColour(lut, 2 * (i + 1), 1, 2 * fitted.length));
      plot.addPoints(time, SimpleArrayUtils.toFloat(fitted[i]), Plot.CIRCLE);
    }
    plot.setColor(Color.black);
    plot.setLimits(dt * 0.5, dt * maxT + dt * 0.5, 0, 1);

    ImageJUtils.display(title, plot);
  }

  /**
   * Fit the probability of remaining in the depth of field.
   *
   * <p>Fit correction coefficient a and b.
   *
   * <pre>
   * dz_corr = dz + a * sqrt(D) + b
   * </pre>
   *
   * <p>dz_corr is used to compute the probability of remaining within the depth-of-field P(dt,
   * dz_corr, D) and fit to the empirical distribution from the simulation.
   *
   * @param probability the probability
   * @param threadCount the thread count
   * @param executor the executor
   * @return the coefficients
   */
  private double[] fitRemaining(double[][] probability, int threadCount, ExecutorService executor) {
    if (settings.repeats == 0) {
      // Allows arbitrary plotting of values. Typically these are zero.
      return new double[] {settings.minA, settings.minB};
    }
    IJ.showStatus("Fitting probability...");

    final double dt = settings.exposureTime / 1000;
    final double dz = settings.depthOfField / 1000;

    final double[] diffusionCoefficients = createDiffusionCoefficients();

    final MultivariateFunction function =
        new DepthOfFieldFunction(dt, dz, diffusionCoefficients, probability, threadCount, executor);

    final LocalList<OptimizationData> args = new LocalList<>();
    // CMA-ES can take 30 N to 300 N^2 evaluations for N parameters
    args.add(new MaxEval(3000));

    final Logger logger = ImageJPluginLoggerHelper.getLogger(getClass());

    logger.log(Level.INFO,
        () -> Arrays.stream(diffusionCoefficients).mapToObj(d -> MathUtils.rounded(d))
            .collect(Collectors.joining(", ", "Diffusion coefficients: ", "")));

    double best = Double.POSITIVE_INFINITY;
    double[] result = new double[2];

    args.add(new ObjectiveFunction(function));
    args.add(new SimpleBounds(new double[2],
        new double[] {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY}));
    args.add(GoalType.MINIMIZE);

    final UniformRandomProvider rng = UniformRandomProviders.create();

    MultivariateOptimizer optimizer;
    if (settings.optimiserMode == OPT_MODE_CMAES) {
      optimizer = createCmaesOptimizer(rng);
      // The sigma determines the search range for the variables.
      // It should be 1/3 of the initial search region.
      args.add(new CMAESOptimizer.Sigma(new double[] {0.02, 0.01}));
      args.add(new CMAESOptimizer.PopulationSize((int) (4 + Math.floor(3 * Math.log(2)))));
    } else if (settings.optimiserMode == OPT_MODE_BOBYQA) {
      optimizer = createBobyqaOptimizer(2);
    } else {
      optimizer = createCustomPowellOptimizer();
      args.add(new CustomPowellOptimizer.BasisStep(new double[] {0.02, 0.01}));
    }

    final Ticker ticker = ImageJUtils.createTicker(settings.repeats, 1);
    for (int n = 1; n <= settings.repeats; n++) {
      try {
        // Guess in range
        final double[] start = {rng.nextDouble(settings.minA, settings.maxA),
            rng.nextDouble(settings.minB, settings.maxB),};

        args.push(new InitialGuess(start));
        final PointValuePair solution =
            optimizer.optimize(args.toArray(new OptimizationData[args.size()]));
        args.pop();

        final double[] r = solution.getPointRef();

        LoggerUtils.log(logger, Level.INFO, "[%d] Fit [%.3f, %.3f]: SS = %s (%d evaluations)", n,
            r[0], r[1], solution.getValue(), optimizer.getEvaluations());
        if (solution.getValue() < best) {
          best = solution.getValue();
          result = solution.getPointRef();
        }
      } catch (final TooManyEvaluationsException ex) {
        LoggerUtils.log(logger, Level.INFO, "Optimiser [%d] failed : Too many evaluation (%d)", n,
            optimizer.getEvaluations());
      } catch (final TooManyIterationsException ex) {
        LoggerUtils.log(logger, Level.INFO, "Optimiser [%d] failed : Too many iterations (%d)", n,
            optimizer.getIterations());
      } catch (final ConvergenceException ex) {
        LoggerUtils.log(logger, Level.INFO, "Powell optimiser [%d] failed to fit : %s", n,
            ex.getMessage());
      }
      ticker.tick();
    }

    return result;
  }

  private void addToResultTable(double[] ab) {
    createTable().append(addResult(ab));
  }

  private TextWindow createTable() {
    return ImageJUtils.refresh(TABLE_REF, () -> {
      return new TextWindow(TITLE + " Analysis", createHeader(), "", 1200, 300);
    });
  }

  private String createHeader() {
    return Arrays
        .stream(new String[] {"dz (nm)", "dt (ms)", "g", "n", "max t", "detector", "restarts",
            "half-life", "min D (um^2/s)", "max D (um^2/s)", "sample D", "a", "b",})
        .collect(Collectors.joining("\t"));
  }

  private String addResult(double[] ab) {
    final StringBuilder sb = new StringBuilder();
    //@formatter:off
    sb.append(MathUtils.rounded(settings.depthOfField)).append('\t')
      .append(MathUtils.rounded(settings.exposureTime)).append('\t')
      .append(settings.gap).append('\t')
      .append(settings.numberOfMolecules).append('\t')
      .append(settings.maxT).append('\t')
      .append(settings.useDetector).append('\t')
      .append(settings.allowRestarts).append('\t')
      .append(MathUtils.rounded(settings.halfLife)).append('\t')
      .append(MathUtils.rounded(settings.minD)).append('\t')
      .append(MathUtils.rounded(settings.maxD)).append('\t')
      .append(settings.sampleD).append('\t')
      .append(MathUtils.rounded(ab[0])).append('\t')
      .append(MathUtils.rounded(ab[1]))
      ;
    //@formatter:on
    return sb.toString();
  }

  private static CustomPowellOptimizer createCustomPowellOptimizer() {
    final double rel = 1e-8;
    final double abs = 1e-100;
    final boolean basisConvergence = false;
    return new CustomPowellOptimizer(rel, abs, null, basisConvergence);
  }

  private static CMAESOptimizer createCmaesOptimizer(UniformRandomProvider rng) {
    final double rel = 1e-8;
    final double abs = 1e-100;
    final int maxIterations = 2000;
    final double stopFitness = 0;
    final boolean isActiveCma = true;
    final int diagonalOnly = 20;
    final int checkFeasableCount = 1;
    final RandomGenerator random = new RandomGeneratorAdapter(rng);
    final boolean generateStatistics = false;
    final ConvergenceChecker<PointValuePair> checker = new SimpleValueChecker(rel, abs);
    return new CMAESOptimizer(maxIterations, stopFitness, isActiveCma, diagonalOnly,
        checkFeasableCount, random, generateStatistics, checker);
  }

  private static BOBYQAOptimizer createBobyqaOptimizer(int n) {
    final int numberOfInterpolationpoints = n + 2;
    // Optional parameters: initialTrustRegionRadius; stoppingTrustRegionRadius.
    // It is unclear what the initial radius should be. The default is 10.
    // When the optimiser starts it sets the value to the minimum bound difference
    // [upper - lower] / 3. So when using a tight bounds this settings does not
    // matter unless we require it to be smaller.
    // The default stopping trust radius is 1e-8.
    return new BOBYQAOptimizer(numberOfInterpolationpoints);
  }

  private double[][] createProbability(double[] ab) {
    final double dt = settings.exposureTime / 1000;
    final double dz = settings.depthOfField / 1000;
    final double[] diffusionCoefficients = createDiffusionCoefficients();

    final double a = ab[0];
    final double b = ab[1];
    final double[][] probability = new double[diffusionCoefficients.length][settings.maxT];

    for (int i = 0; i < probability.length; i++) {
      final double d = diffusionCoefficients[i];
      for (int j = 0; j < probability[i].length; j++) {
        final double dz_corr = dz + a * Math.sqrt(d) + b;
        probability[i][j] = DiffusionAnalysis.remaining(dt * (j + 1), dz_corr, d);
      }
    }

    return probability;
  }

  /**
   * Function to evaluate the squared distance between the simulated probability and the computed
   * probability using P(dt, dz_corr, D) where z_corr is dz + a * sqrt(D) + b.
   *
   * <p>This function is slow and can use an {@link ExecutorService} for parallel evaluation.
   */
  static class DepthOfFieldFunction implements MultivariateFunction {
    double dt;
    double dz;
    double[] diffusionCoefficients;
    double[][] probability;
    int threadCount;
    ExecutorService executor;

    DepthOfFieldFunction(double dt, double dz, double[] diffusionCoefficients,
        double[][] probability, int threadCount, ExecutorService executor) {
      this.dt = dt;
      this.dz = dz;
      this.diffusionCoefficients = diffusionCoefficients;
      this.probability = probability;
      this.threadCount = threadCount;
      this.executor = executor;
    }

    @Override
    public double value(double[] point) {
      final double a = point[0];
      final double b = point[1];
      return value(dz, a, b);
    }

    private double value(double dz, double a, double b) {
      if (threadCount < 2 || executor == null) {
        double ss = 0;
        for (int i = 0; i < probability.length; i++) {
          final double d = diffusionCoefficients[i];
          for (int j = 0; j < probability[i].length; j++) {
            final double dz_corr = dz + a * Math.sqrt(d) + b;
            final double p = DiffusionAnalysis.remaining(dt * (j + 1), dz_corr, d);
            final double dp = probability[i][j] - p;
            ss += dp * dp;
          }
        }
        return ss;
      }

      final List<Future<Double>> futures = new LinkedList<>();
      final int maxj = probability[0].length;
      final AtomicInteger position = new AtomicInteger(probability.length * maxj);

      for (int n = 0; n < threadCount; n++) {
        futures.add(executor.submit(() -> {
          double ss = 0;
          for (;;) {
            final int p = position.decrementAndGet();
            if (p < 0) {
              break;
            }
            final int i = p / maxj;
            final int j = p % maxj;

            final double d = diffusionCoefficients[i];
            final double dz_corr = dz + a * Math.sqrt(d) + b;
            final double dp =
                probability[i][j] - DiffusionAnalysis.remaining(dt * (j + 1), dz_corr, d);
            ss += dp * dp;
          }
          return ss;
        }));
      }

      // Finish processing data
      double ss = 0;
      for (final Future<Double> f : futures) {
        try {
          ss += f.get();
        } catch (InterruptedException | ExecutionException e) {
          throw new RuntimeException(e);
        }
      }
      return ss;
    }
  }
}
