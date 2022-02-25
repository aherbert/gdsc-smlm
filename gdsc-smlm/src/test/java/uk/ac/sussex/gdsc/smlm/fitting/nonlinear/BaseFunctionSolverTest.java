/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.fitting.nonlinear;

import java.util.Arrays;
import java.util.concurrent.ConcurrentHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.stat.inference.TTest;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.ContinuousSampler;
import org.apache.commons.rng.sampling.distribution.NormalizedGaussianSampler;
import org.apache.commons.rng.sampling.distribution.SharedStateContinuousSampler;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import uk.ac.sussex.gdsc.core.math.SimpleArrayMoment;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.rng.MarsagliaTsangGammaSampler;
import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;
import uk.ac.sussex.gdsc.smlm.GdscSmlmTestUtils;
import uk.ac.sussex.gdsc.smlm.fitting.FisherInformationMatrix;
import uk.ac.sussex.gdsc.smlm.fitting.FitStatus;
import uk.ac.sussex.gdsc.smlm.fitting.FunctionSolver;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.PoissonGradientProcedure;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.PoissonGradientProcedureUtils;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;
import uk.ac.sussex.gdsc.smlm.function.Gradient2Function;
import uk.ac.sussex.gdsc.smlm.function.OffsetGradient2Function;
import uk.ac.sussex.gdsc.smlm.function.StandardValueProcedure;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils.TestLevel;
import uk.ac.sussex.gdsc.test.utils.functions.FunctionUtils;

/**
 * Base class for testing the function solvers.
 */
@SuppressWarnings({"javadoc"})
public abstract class BaseFunctionSolverTest {
  protected static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(BaseFunctionSolverTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  // Basic Gaussian
  static double[] params = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
  static double[] base = {0.8, 1, 1.2}; // Applied (*) to the background
  //@formatter:off
  static double[] signal = {
      1000, 2000, 5000, 10000
      //100, 200, 400, 800
    };
  //@formatter:on
  static double[] shift = {-0.5, 0, 0.5}; // Applied (+/-) to the x/y position
  static double[] factor = {0.9, 1, 1.1}; // Applied (*) to the width
  static int size = 11;

  static {
    // Keep SNR reasonable. This should be an "easy" test since the bounds
    // for a correct answer are strict.
    final double minSnr = 10;
    final double sd = 1.3;
    final double mean = Gaussian2DPeakResultHelper.getMeanSignalUsingP05(signal[0], sd, sd);
    // snr = mean/background => background = mean/snr
    params[Gaussian2DFunction.BACKGROUND] = mean / minSnr;
    params[Gaussian2DFunction.X_POSITION] = size / 2;
    params[Gaussian2DFunction.Y_POSITION] = size / 2;
    params[Gaussian2DFunction.X_SD] = sd;
  }

  static double[] defaultClampValues;

  static {
    defaultClampValues = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    // Taken from the 3D-DAO-STORM paper:
    // (Babcock et al. 2012) A high-density 3D localization algorithm for stochastic optical
    // reconstruction microscopy. Optical Nanoscopy. 2012 1:6
    // DOI: 10.1186/2192-2853-1-6
    // Page 3
    // Note: It is not clear if the background/signal are in ADUs or photons. I assume photons.
    // Note: The clamp value is the maximum permitted single step.
    // If the desired step is equal to the maximum step then the clamped step will be half.

    // This seems big for background in photons
    defaultClampValues[Gaussian2DFunction.BACKGROUND] = 100;
    // defaultClampValues[Gaussian2DFunction.BACKGROUND] = 20;
    defaultClampValues[Gaussian2DFunction.SIGNAL] = 1000;
    defaultClampValues[Gaussian2DFunction.X_POSITION] = 1;
    defaultClampValues[Gaussian2DFunction.Y_POSITION] = 1;
    defaultClampValues[Gaussian2DFunction.Z_POSITION] = 1;
    // This seems big for width changes
    defaultClampValues[Gaussian2DFunction.X_SD] = 3;
    defaultClampValues[Gaussian2DFunction.Y_SD] = 3;
    defaultClampValues[Gaussian2DFunction.ANGLE] = Math.PI;

    // More restrictive ...

    defaultClampValues[Gaussian2DFunction.BACKGROUND] = 5;
    // defaultClampValues[Gaussian2DFunction.SIGNAL] = 1000;
    // defaultClampValues[Gaussian2DFunction.X_POSITION] = 1;
    // defaultClampValues[Gaussian2DFunction.Y_POSITION] = 1;
    defaultClampValues[Gaussian2DFunction.X_SD] = 1;
    defaultClampValues[Gaussian2DFunction.Y_SD] = 1;
  }

  enum NoiseModel {
    NONE, CCD, EMCCD, SCMOS
  }

  // Based on Huang et al (2015) sCMOS per pixel noise.
  // Variance = Exponential (equivalent to chi-squared with k=1, i.e.
  // sum of the squares of 1 normal distribution).
  // We want 99.9% @ 400 ADU based on supplementary figure 1.a/1.b
  // cumul = 1 - e^-lx (with l = 1/mean)
  // => e^-lx = 1 - cumul
  // => -lx = log(1-0.999)
  // => l = -log(0.001) / 400 (since x==400)
  // => 1/l = 57.9
  private static double variance = 57.9; // SD = 7.6

  // Gain = Approximately Normal
  private static double gain = 2.2;
  private static double gainSD = 0.2;

  // Other noise models
  private static double noiseCCD = Math.sqrt(variance / (gain * gain)); // Same as sCMOS
  private static double emGain = 300;
  private static double noiseEMCCD = 0.02;

  private static double[][] weights = new double[NoiseModel.values().length][];
  private static double[][] noise = new double[weights.length][];

  private static ConcurrentHashMap<RandomSeed, double[][]> ConcurrentHashMap =
      new ConcurrentHashMap<>();

  private static double[] getWeights(RandomSeed seed, NoiseModel noiseModel) {
    final int index = noiseModel.ordinal();
    if (weights[index] == null) {
      if (noiseModel == NoiseModel.SCMOS) {
        // Special case of per-pixel random weights
        final double[][] data =
            ConcurrentHashMap.computeIfAbsent(seed, BaseFunctionSolverTest::createData);
        return data[0];
      }
      computeWeights(noiseModel, index);
    }
    return weights[index];
  }

  private static double[] getNoise(RandomSeed seed, NoiseModel noiseModel) {
    final int index = noiseModel.ordinal();
    if (noise[index] == null) {
      if (noiseModel == NoiseModel.SCMOS) {
        // Special case of per-pixel random noise
        final double[][] data =
            ConcurrentHashMap.computeIfAbsent(seed, BaseFunctionSolverTest::createData);
        return data[1];
      }
      computeWeights(noiseModel, index);
    }
    return noise[index];
  }

  private static void computeWeights(NoiseModel noiseModel, int index) {
    if (noiseModel == NoiseModel.SCMOS) {
      throw new IllegalStateException("This requires a random generator");
    }
    // The rest are fixed for all pixels
    final double[] w = new double[size * size];
    final double[] n = new double[size * size];
    switch (noiseModel) {
      case CCD:
        computeWeights(w, n, noiseCCD);
        break;
      case EMCCD:
        computeWeights(w, n, noiseEMCCD);
        break;
      case NONE:
      default:
        // Nothing to do
        break;
    }
    noise[index] = n;
    weights[index] = w;
  }

  private static void computeWeights(double[] weights, double[] noise, double sd) {
    Arrays.fill(weights, sd * sd);
    Arrays.fill(noise, sd);
  }

  private static double[][] createData(RandomSeed source) {
    // Per observation read noise.
    // This is generated once so create the randon generator here.
    final UniformRandomProvider rg = RngUtils.create(source.get());
    final ContinuousSampler ed = SamplerUtils.createExponentialSampler(rg, variance);
    final SharedStateContinuousSampler gs = SamplerUtils.createGaussianSampler(rg, gain, gainSD);
    final double[] w = new double[size * size];
    final double[] n = new double[size * size];
    for (int i = 0; i < w.length; i++) {
      final double pixelVariance = ed.sample();
      final double pixelGain = Math.max(0.1, gs.sample());
      // weights = var / g^2
      w[i] = pixelVariance / (pixelGain * pixelGain);

      // Convert to standard deviation for simulation
      n[i] = Math.sqrt(w[i]);
    }
    return new double[][] {w, n};
  }

  void canFitSingleGaussian(RandomSeed seed, FunctionSolver solver, boolean applyBounds) {
    // This is here to support the old solver tests which used to have fixed noise of 0.1, 0.5, 1.
    // Those levels are irrelevant for modern EM-CCD cameras which have effectively very low noise.
    // The test has been changed to better simulate the cameras encountered.
    canFitSingleGaussian(seed, solver, applyBounds, NoiseModel.NONE);
    // We are not interested in high noise CCD so this is commented out
    // canFitSingleGaussian(solver, applyBounds, NoiseModel.CCD);
    canFitSingleGaussian(seed, solver, applyBounds, NoiseModel.EMCCD);
    canFitSingleGaussian(seed, solver, applyBounds, NoiseModel.SCMOS);
  }

  void canFitSingleGaussian(RandomSeed seed, FunctionSolver solver, boolean applyBounds,
      NoiseModel noiseModel) {
    // Allow reporting the fit deviations
    final boolean report = logger.isLoggable(TestLevel.TEST_INFO);
    double[] crlb = null;
    SimpleArrayMoment moment = null;

    final double[] noise = getNoise(seed, noiseModel);
    if (solver.isWeighted()) {
      solver.setWeights(getWeights(seed, noiseModel));
    }

    final UniformRandomProvider rg = RngUtils.create(seed.get());

    for (final double s : signal) {
      final double[] expected = createParams(1, s, 0, 0, 1);
      final double[] lower = createParams(0, s * 0.5, -0.3, -0.3, 0.8);
      final double[] upper = createParams(3, s * 2, 0.3, 0.3, 1.2);
      if (applyBounds) {
        solver.setBounds(lower, upper);
      }
      if (report) {
        // Compute the CRLB for a Poisson process
        final PoissonGradientProcedure gp = PoissonGradientProcedureUtils
            .create((Gradient1Function) ((BaseFunctionSolver) solver).getGradientFunction());
        gp.computeFisherInformation(expected);
        final FisherInformationMatrix f =
            new FisherInformationMatrix(gp.getLinear(), gp.numberOfGradients);
        crlb = f.crlbSqrt();
        // Compute the deviations.
        // Note this is not the same as the CRLB as the fit is repeated
        // against the same data.
        // It should be repeated against different data generated with constant
        // parameters and variable noise.
        moment = new SimpleArrayMoment();
      }
      final double[] data = drawGaussian(expected, noise, noiseModel, rg);
      for (final double db : base) {
        for (final double dx : shift) {
          for (final double dy : shift) {
            for (final double dsx : factor) {
              final double[] p = createParams(db, s, dx, dy, dsx);
              final double[] fp = fitGaussian(solver, data, p, expected);
              for (int i = 0; i < expected.length; i++) {
                if (fp[i] < lower[i]) {
                  Assertions
                      .fail(FunctionUtils.getSupplier("Fit Failed: [%d] %.2f < %.2f: %s != %s", i,
                          fp[i], lower[i], Arrays.toString(fp), Arrays.toString(expected)));
                }
                if (fp[i] > upper[i]) {
                  Assertions
                      .fail(FunctionUtils.getSupplier("Fit Failed: [%d] %.2f > %.2f: %s != %s", i,
                          fp[i], upper[i], Arrays.toString(fp), Arrays.toString(expected)));
                }
                if (report) {
                  fp[i] = expected[i] - fp[i];
                }
              }
              // Store the deviations
              if (report) {
                moment.add(fp);
              }
            }
          }
        }
      }
      // Report
      if (report) {
        logger.log(TestLevel.TEST_INFO,
            FunctionUtils.getSupplier("%s %s %f : CRLB = %s, Deviations = %s",
                solver.getClass().getSimpleName(), noiseModel, s, Arrays.toString(crlb),
                Arrays.toString(moment.getStandardDeviation())));
      }
    }
  }

  void canFitSingleGaussianBetter(RandomSeed seed, FunctionSolver solver, boolean applyBounds,
      FunctionSolver solver2, boolean applyBounds2, String name, String name2) {
    // This is here to support the old solver tests which used to have fixed noise of 0.1, 0.5, 1.
    // Those levels are irrelevant for modern EM-CCD cameras which have effectively very low noise.
    // The test has been changed to better simulate the cameras encountered.
    canFitSingleGaussianBetter(seed, solver, applyBounds, solver2, applyBounds2, name, name2,
        NoiseModel.NONE);
    // We are not interested in high noise CCD so this is commented out
    // canFitSingleGaussianBetter(solver, applyBounds, solver2, applyBounds2, name, name2,
    // NoiseModel.CCD);
    canFitSingleGaussianBetter(seed, solver, applyBounds, solver2, applyBounds2, name, name2,
        NoiseModel.EMCCD);
  }

  void canFitSingleGaussianBetter(RandomSeed seed, FunctionSolver solver, boolean applyBounds,
      FunctionSolver solver2, boolean applyBounds2, String name, String name2,
      NoiseModel noiseModel) {
    final double[] noise = getNoise(seed, noiseModel);
    if (solver.isWeighted()) {
      solver.setWeights(getWeights(seed, noiseModel));
    }

    final int loops = 5;
    final UniformRandomProvider rg = RngUtils.create(seed.get());
    final StoredDataStatistics[] stats = new StoredDataStatistics[6];
    final String[] statName = {"Signal", "X", "Y"};

    final int[] betterPrecision = new int[3];
    final int[] totalPrecision = new int[3];
    final int[] betterAccuracy = new int[3];
    final int[] totalAccuracy = new int[3];

    final String msg = "%s vs %s : %.1f (%s) %s %f +/- %f vs %f +/- %f  (N=%d) %b %s";

    int i1 = 0;
    int i2 = 0;
    for (final double s : signal) {
      final double[] expected = createParams(1, s, 0, 0, 1);
      double[] lower = null;
      double[] upper = null;
      if (applyBounds || applyBounds2) {
        lower = createParams(0, s * 0.5, -0.3, -0.3, 0.8);
        upper = createParams(3, s * 2, 0.3, 0.3, 1.2);
      }
      if (applyBounds) {
        solver.setBounds(lower, upper);
      }
      if (applyBounds2) {
        solver2.setBounds(lower, upper);
      }

      for (int loop = loops; loop-- > 0;) {
        final double[] data = drawGaussian(expected, noise, noiseModel, rg);

        for (int i = 0; i < stats.length; i++) {
          stats[i] = new StoredDataStatistics();
        }

        for (final double db : base) {
          for (final double dx : shift) {
            for (final double dy : shift) {
              for (final double dsx : factor) {
                final double[] p = createParams(db, s, dx, dy, dsx);
                final double[] fp = fitGaussian(solver, data, p, expected);
                i1 += solver.getEvaluations();

                final double[] fp2 = fitGaussian(solver2, data, p, expected);
                i2 += solver2.getEvaluations();

                // Get the mean and sd (the fit precision)
                compare(fp, expected, fp2, expected, Gaussian2DFunction.SIGNAL, stats[0], stats[1]);

                compare(fp, expected, fp2, expected, Gaussian2DFunction.X_POSITION, stats[2],
                    stats[3]);
                compare(fp, expected, fp2, expected, Gaussian2DFunction.Y_POSITION, stats[4],
                    stats[5]);

                // Use the distance
                // stats[2].add(distance(fp, expected));
                // stats[3].add(distance(fp2, expected2));
              }
            }
          }
        }

        final double alpha = 0.05; // two sided
        for (int i = 0; i < stats.length; i += 2) {
          double u1 = stats[i].getMean();
          double u2 = stats[i + 1].getMean();
          final double sd1 = stats[i].getStandardDeviation();
          final double sd2 = stats[i + 1].getStandardDeviation();

          final TTest tt = new TTest();
          final boolean diff = tt.tTest(stats[i].getValues(), stats[i + 1].getValues(), alpha);

          final int index = i / 2;
          final Object[] args = new Object[] {name2, name, s, noiseModel, statName[index], u2, sd2,
              u1, sd1, stats[i].getN(), diff, ""};
          if (diff) {
            // Different means. Check they are roughly the same
            if (DoubleEquality.almostEqualRelativeOrAbsolute(u1, u2, 0.1, 0)) {
              // Basically the same. Check which is more precise
              if (!DoubleEquality.almostEqualRelativeOrAbsolute(sd1, sd2, 0.05, 0)) {
                if (sd2 < sd1) {
                  betterPrecision[index]++;
                  args[args.length - 1] = "P*";
                  logger.log(TestLogUtils.getRecord(TestLevel.TEST_DEBUG, msg, args));
                } else {
                  args[args.length - 1] = "P";
                  logger.log(TestLogUtils.getRecord(TestLevel.TEST_DEBUG, msg, args));
                }
                totalPrecision[index]++;
              }
            } else {
              // Check which is more accurate (closer to zero)
              u1 = Math.abs(u1);
              u2 = Math.abs(u2);
              if (u2 < u1) {
                betterAccuracy[index]++;
                args[args.length - 1] = "A*";
                logger.log(TestLogUtils.getRecord(TestLevel.TEST_DEBUG, msg, args));
              } else {
                args[args.length - 1] = "A";
                logger.log(TestLogUtils.getRecord(TestLevel.TEST_DEBUG, msg, args));
              }
              totalAccuracy[index]++;
            }

            // The same means. Check that it is more precise
          } else if (!DoubleEquality.almostEqualRelativeOrAbsolute(sd1, sd2, 0.05, 0)) {
            if (sd2 < sd1) {
              betterPrecision[index]++;
              args[args.length - 1] = "P*";
              logger.log(TestLogUtils.getRecord(TestLevel.TEST_DEBUG, msg, args));
            } else {
              args[args.length - 1] = "P";
              logger.log(TestLogUtils.getRecord(TestLevel.TEST_DEBUG, msg, args));
            }
            totalPrecision[index]++;
          }
        }
      }
    }

    int better = 0;
    int total = 0;
    for (int index = 0; index < statName.length; index++) {
      better += betterPrecision[index] + betterAccuracy[index];
      total += totalPrecision[index] + totalAccuracy[index];
      test(name2, name, statName[index] + " P", betterPrecision[index], totalPrecision[index],
          TestLevel.TEST_INFO);
      test(name2, name, statName[index] + " A", betterAccuracy[index], totalAccuracy[index],
          TestLevel.TEST_INFO);
    }
    test(name2, name, String.format("All (eval [%d] [%d]) : ", i1, i2), better, total,
        TestLevel.TEST_INFO);
  }

  private static void test(String name2, String name, String statName, int better, int total,
      Level logLevel) {
    final double p = (total == 0) ? 0 : 100.0 * better / total;
    logger.log(logLevel, FunctionUtils.getSupplier("%s vs %s : %s %d / %d  (%.1f)", name2, name,
        statName, better, total, p));
    // Do not test if we don't have many examples
    if (total <= 10) {
      return;
    }

    // Disable this for now so builds do not fail during the test phase

    // It seems that most of the time clamping and bounds improve things.
    // There are a few cases where Bounds or Clamping alone do not improve things.
    // Use of Dynamic Clamping is always better.
    // Use of Bounded Dynamic Clamping is always better.

    // The test may be unrealistic as the initial params are close to the actual answer.

    // if (p<50)
    // ExtraAssertions.fail(("%s vs %s : %s %d / %d (%.1f)", name2, name, statName, better, total,
    // p));
  }

  private static void compare(double[] o1, double[] e1, double[] o2, double[] e2, int index,
      Statistics stats1, Statistics stats2) {
    compare(o1[index], e1[index], o2[index], e2[index], stats1, stats2);
  }

  private static void compare(double o1, double e1, double o2, double e2, Statistics stats1,
      Statistics stats2) {
    stats1.add(o1 - e1);
    stats2.add(o2 - e2);
  }

  double[] createParams(double db, double signal, double dx, double dy, double dsx) {
    final double[] p = params.clone();
    p[Gaussian2DFunction.BACKGROUND] *= db;
    p[Gaussian2DFunction.SIGNAL] = signal;
    p[Gaussian2DFunction.X_POSITION] += dx;
    p[Gaussian2DFunction.Y_POSITION] += dy;
    p[Gaussian2DFunction.X_SD] *= dsx;
    return p;
  }

  double[] addBiasToParams(double[] params, double bias) {
    params = params.clone();
    params[Gaussian2DFunction.BACKGROUND] += bias;
    return params;
  }

  double[] fitGaussian(FunctionSolver solver, double[] data, double[] params, double[] expected) {
    // logger.log(TestLevel.TEST_INFO, FunctionUtils.getSupplier("%s : Expected %s",
    // solver.getClass().getSimpleName(), Arrays.toString(expected)));
    params = params.clone();
    final FitStatus status = solver.fit(data, null, params, null);
    if (status != FitStatus.OK) {
      Assertions.fail(FunctionUtils.getSupplier("Fit Failed: %s i=%d: %s != %s", status.toString(),
          solver.getIterations(), Arrays.toString(params), Arrays.toString(expected)));
    }
    return params;
  }

  static double[] drawGaussian(double[] params, UniformRandomProvider rg) {
    return drawGaussian(params, null, NoiseModel.NONE, rg);
  }

  static double[] drawGaussian(double[] params, double[] noise, UniformRandomProvider rg) {
    return drawGaussian(params, noise, NoiseModel.NONE, rg);
  }

  static final int flags = GaussianFunctionFactory.FIT_ERF_CIRCLE;

  /**
   * Draw a Gaussian with Poisson shot noise and Gaussian read noise.
   *
   * @param params The Gaussian parameters
   * @param noise The read noise
   * @param noiseModel the noise model
   * @param rg the random generator
   * @return The data
   */
  static double[] drawGaussian(double[] params, double[] noise, NoiseModel noiseModel,
      UniformRandomProvider rg) {
    final int n = params.length / Gaussian2DFunction.PARAMETERS_PER_PEAK;
    final Gaussian2DFunction f = GaussianFunctionFactory.create2D(n, size, size, flags, null);
    final double[] data = f.computeValues(params);

    // Poisson noise
    for (int i = 0; i < data.length; i++) {
      if (data[i] > 0) {
        data[i] = GdscSmlmTestUtils.createPoissonSampler(rg, data[i]).sample();
      }
    }

    // Simulate EM-gain
    if (noiseModel == NoiseModel.EMCCD) {
      // Use a gamma distribution
      // Since the call random.nextGamma(...) creates a Gamma distribution
      // which pre-calculates factors only using the scale parameter we
      // create a custom gamma distribution where the shape can be set as a property.
      final MarsagliaTsangGammaSampler gd = new MarsagliaTsangGammaSampler(rg, 1, emGain);

      for (int i = 0; i < data.length; i++) {
        if (data[i] > 0) {
          gd.setAlpha(data[i]);
          // The sample will amplify the signal so we remap to the original scale
          data[i] = gd.sample() / emGain;
        }
      }
    }

    // Read-noise
    final NormalizedGaussianSampler gs = SamplerUtils.createNormalizedGaussianSampler(rg);
    if (noise != null) {
      for (int i = 0; i < data.length; i++) {
        data[i] += gs.sample() * noise[i];
      }
    }

    // uk.ac.sussex.gdsc.core.ij.Utils.display("Spot", data, size, size);
    return data;
  }

  static double[] p1;
  static double[] p12;
  static double[] p2v;

  static {
    p1 = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    final double[] p2 = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    p12 = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK * 2];
    p1[Gaussian2DFunction.BACKGROUND] = 5;
    p1[Gaussian2DFunction.SIGNAL] = 1000;
    p1[Gaussian2DFunction.X_POSITION] = 3.1;
    p1[Gaussian2DFunction.Y_POSITION] = 4.2;
    p1[Gaussian2DFunction.X_SD] = 1.2;
    // p2[Gaussian2DFunction.BACKGROUND] = p1[Gaussian2DFunction.BACKGROUND];
    p2[Gaussian2DFunction.SIGNAL] = 600;
    p2[Gaussian2DFunction.X_POSITION] = 7.3;
    p2[Gaussian2DFunction.Y_POSITION] = 8.4;
    p2[Gaussian2DFunction.X_SD] = 1.1;
    System.arraycopy(p1, 0, p12, 0, p1.length);
    System.arraycopy(p2, 1, p12, p1.length, Gaussian2DFunction.PARAMETERS_PER_PEAK);

    final StandardValueProcedure p = new StandardValueProcedure();
    p2v = p.getValues(GaussianFunctionFactory.create2D(1, size, size, flags, null), p2);
  }

  /**
   * Check the fit and compute deviations match. The first solver will be used to do the fit. This
   * is initialised from the solution so the convergence criteria can be set to accept the first
   * step. The second solver is used to compute deviations (thus is not initialised for fitting).
   *
   * @param seed the seed
   * @param solver1 the solver 1
   * @param solver2 the solver 2
   * @param noiseModel the noise model
   * @param useWeights the use weights
   */
  void fitAndComputeDeviationsMatch(RandomSeed seed, BaseFunctionSolver solver1,
      BaseFunctionSolver solver2, NoiseModel noiseModel, boolean useWeights) {
    final double[] noise = getNoise(seed, noiseModel);
    if (solver1.isWeighted() && useWeights) {
      solver1.setWeights(getWeights(seed, noiseModel));
      solver2.setWeights(getWeights(seed, noiseModel));
    }

    // Draw target data
    final UniformRandomProvider rg = RngUtils.create(seed.get());
    final double[] data = drawGaussian(p12, noise, noiseModel, rg);

    // fit with 2 peaks using the known params.
    // compare to 2 peak deviation computation.
    final Gaussian2DFunction f2 = GaussianFunctionFactory.create2D(2, size, size, flags, null);
    solver1.setGradientFunction(f2);
    solver2.setGradientFunction(f2);
    double[] params = p12.clone();
    double[] expected = new double[params.length];
    double[] observed = new double[params.length];
    solver1.fit(data, null, params, expected);
    // System.out.TestLog.fine(logger,"a="+Arrays.toString(a));
    solver2.computeDeviations(data, params, observed);

    // System.out.TestLog.fine(logger,"e2="+Arrays.toString(e));
    // System.out.TestLog.fine(logger,"o2="+Arrays.toString(o));
    Assertions.assertArrayEquals(observed, expected,
        "Fit 2 peaks and deviations 2 peaks do not match");

    // Try again with y-fit values
    params = p12.clone();
    final double[] o1 = new double[f2.size()];
    final double[] o2 = new double[o1.length];
    solver1.fit(data, o1, params, expected);
    // System.out.TestLog.fine(logger,"a="+Arrays.toString(a));
    solver2.computeValue(data, o2, params);

    Assertions.assertArrayEquals(observed, expected,
        "Fit 2 peaks with yFit and deviations 2 peaks do not match");

    final StandardValueProcedure p = new StandardValueProcedure();
    double[] ev = p.getValues(f2, params);
    Assertions.assertArrayEquals(ev, o1, 1e-8, "Fit 2 peaks yFit");
    Assertions.assertArrayEquals(ev, o2, 1e-8, "computeValue 2 peaks yFit");

    if (solver1 instanceof SteppingFunctionSolver) {
      // fit with 1 peak + 1 precomputed using the known params.
      // compare to 2 peak deviation computation.
      final ErfGaussian2DFunction f1 =
          (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(1, size, size, flags, null);
      final Gradient2Function pf1 = OffsetGradient2Function.wrapGradient2Function(f1, p2v);
      solver1.setGradientFunction(pf1);
      params = p1.clone();
      expected = new double[params.length];
      solver1.fit(data, null, params, expected);

      final double[] a2 = p12.clone(); // To copy the second peak
      System.arraycopy(params, 0, a2, 0, params.length); // Add the same fitted first peak
      solver2.computeDeviations(data, a2, observed);
      // System.out.TestLog.fine(logger,"e1p1=" + Arrays.toString(e));
      // System.out.TestLog.fine(logger,"o2=" + Arrays.toString(o));

      // Deviation should be lower with only 1 peak.
      // Due to matrix inversion this may not be the case for all parameters so count.
      int ok = 0;
      int fail = 0;
      final StringBuilder sb = new StringBuilder();
      for (int i = 0; i < expected.length; i++) {
        if (expected[i] <= observed[i]) {
          ok++;
          continue;
        }
        fail++;
        TextUtils.formatTo(sb,
            "Fit 1 peak + 1 precomputed is higher than deviations 2 peaks %s: %s > %s",
            Gaussian2DFunction.getName(i), expected[i], observed[i]);
      }
      if (fail > ok) {
        Assertions.fail(sb.toString());
      }

      // Try again with y-fit values
      params = p1.clone();
      Arrays.fill(o1, 0);
      Arrays.fill(o2, 0);
      observed = new double[params.length];
      solver1.fit(data, o1, params, observed);
      solver2.computeValue(data, o2, a2);

      Assertions.assertArrayEquals(observed, expected, 1e-8,
          "Fit 1 peak + 1 precomputed with yFit and deviations 1 peak + "
              + "1 precomputed do not match");

      ev = p.getValues(pf1, params);
      Assertions.assertArrayEquals(ev, o1, 1e-8, "Fit 1 peak + 1 precomputed yFit");
      Assertions.assertArrayEquals(ev, o2, 1e-8, "computeValue 1 peak + 1 precomputed yFit");
    }
  }

  /**
   * Check the fit and compute values match. The first solver will be used to do the fit. This is
   * initialised from the solution so the convergence criteria can be set to accept the first step.
   * The second solver is used to compute values (thus is not initialised for fitting).
   *
   * @param seed the seed
   * @param solver1 the solver
   * @param solver2 the solver 2
   * @param noiseModel the noise model
   * @param useWeights the use weights
   */
  void fitAndComputeValueMatch(RandomSeed seed, BaseFunctionSolver solver1,
      BaseFunctionSolver solver2, NoiseModel noiseModel, boolean useWeights) {
    final double[] noise = getNoise(seed, noiseModel);
    if (solver1.isWeighted() && useWeights) {
      solver1.setWeights(getWeights(seed, noiseModel));
      solver2.setWeights(getWeights(seed, noiseModel));
    }

    // Draw target data
    final UniformRandomProvider rg = RngUtils.create(seed.get());
    final double[] data = drawGaussian(p12, noise, noiseModel, rg);

    // fit with 2 peaks using the known params.
    final Gaussian2DFunction f2 = GaussianFunctionFactory.create2D(2, size, size, flags, null);
    solver1.setGradientFunction(f2);
    solver2.setGradientFunction(f2);
    double[] params = p12.clone();
    solver1.fit(data, null, params, null);
    solver2.computeValue(data, null, params);

    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-10, 0);

    double v1 = solver1.getValue();
    double v2 = solver2.getValue();
    TestAssertions.assertTest(v1, v2, predicate, "Fit 2 peaks and computeValue");

    final double[] o1 = new double[f2.size()];
    final double[] o2 = new double[o1.length];

    solver1.fit(data, o1, params, null);
    solver2.computeValue(data, o2, params);

    v1 = solver1.getValue();
    v2 = solver2.getValue();
    TestAssertions.assertTest(v1, v2, predicate, "Fit 2 peaks and computeValue with yFit");

    final StandardValueProcedure p = new StandardValueProcedure();
    double[] expected = p.getValues(f2, params);
    Assertions.assertArrayEquals(expected, o1, 1e-8, "Fit 2 peaks yFit");
    Assertions.assertArrayEquals(expected, o2, 1e-8, "computeValue 2 peaks yFit");

    if (solver1 instanceof SteppingFunctionSolver) {
      // fit with 1 peak + 1 precomputed using the known params.
      // compare to 2 peak computation.
      final ErfGaussian2DFunction f1 =
          (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(1, size, size, flags, null);
      final Gradient2Function pf1 = OffsetGradient2Function.wrapGradient2Function(f1, p2v);
      solver1.setGradientFunction(pf1);
      solver2.setGradientFunction(pf1);

      params = p1.clone();
      solver1.fit(data, null, params, null);
      solver2.computeValue(data, null, params);

      v1 = solver1.getValue();
      v2 = solver2.getValue();
      TestAssertions.assertTest(v1, v2, predicate, "Fit 1 peak + 1 precomputed and computeValue");

      Arrays.fill(o1, 0);
      Arrays.fill(o2, 0);

      solver1.fit(data, o1, params, null);
      solver2.computeValue(data, o2, params);

      v1 = solver1.getValue();
      v2 = solver2.getValue();
      TestAssertions.assertTest(v1, v2, predicate,
          "Fit 1 peak + 1 precomputed and computeValue with yFit");

      expected = p.getValues(pf1, params);
      Assertions.assertArrayEquals(expected, o1, 1e-8, "Fit 1 peak + 1 precomputed yFit");
      Assertions.assertArrayEquals(expected, o2, 1e-8, "computeValue 1 peak + 1 precomputed yFit");
    }
  }
}
