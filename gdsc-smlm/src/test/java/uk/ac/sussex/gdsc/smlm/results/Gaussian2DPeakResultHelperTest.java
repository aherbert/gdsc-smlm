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

package uk.ac.sussex.gdsc.smlm.results;

import java.util.logging.Logger;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.data.config.PsfHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.test.api.Predicates;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngFactory;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogging.TestLevel;
import uk.ac.sussex.gdsc.test.utils.TestSettings;
import uk.ac.sussex.gdsc.test.utils.functions.FunctionUtils;

@SuppressWarnings({"javadoc"})
class Gaussian2DPeakResultHelperTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(Gaussian2DPeakResultHelperTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  double[] testA = {100, 130, 160};
  double[] testS = {80, 100, 140};
  double[] testN = {1, 10, 30, 100, 1000};
  double[] testB2 = {0, 1, 2, 4, 8};
  int minpoints = 3;
  int maxpoints = 20;

  @Test
  void canCalculateMaximumLikelihoodVariance() {
    int min = Gaussian2DPeakResultHelper.POINTS;
    int max = min;
    if (TestSettings.allow(TestComplexity.HIGH)) {
      min = 3;
      max = 20;
    }
    for (final double a : testA) {
      for (final double s : testS) {
        for (final double n : testN) {
          for (final double b2 : testB2) {
            for (int points = min; points <= max; points++) {
              Gaussian2DPeakResultHelper.getMLVarianceX(a, s, n, b2, true, points);
            }
          }
        }
      }
    }
  }

  @Test
  void lowerIntegrationpointsApproximateMaximumLikelihoodVariance() {
    Assumptions.assumeTrue(logger.isLoggable(TestLevel.TEST_INFO));
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.HIGH));

    final double[] sum = new double[maxpoints + 1];
    int count = 0;
    for (final double a : testA) {
      for (final double s : testS) {
        for (final double n : testN) {
          for (final double b2 : testB2) {
            count++;
            final double e = Gaussian2DPeakResultHelper.getMLVarianceX(a, s, n, b2, true, 30);
            for (int points = minpoints; points <= maxpoints; points++) {
              final double o = Gaussian2DPeakResultHelper.getMLVarianceX(a, s, n, b2, true, points);
              final double error = DoubleEquality.relativeError(e, o);
              sum[points] += error;
              if (error > 1e-2) {
                Assertions.fail(String.format("a=%f, s=%f, N=%f, b2=%f, points=%d : %f != %f : %f",
                    a, s, n, b2, points, e, o, error));
              }
            }
          }
        }
      }
    }

    for (int points = minpoints; points <= maxpoints; points++) {
      logger.log(TestLevel.TEST_INFO,
          FunctionUtils.getSupplier("Points = %d, Av error = %f", points, sum[points] / count));
    }
  }

  @Test
  void runSpeedTest() {
    Assumptions.assumeTrue(logger.isLoggable(TestLevel.TEST_INFO));
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.HIGH));

    // Test with realistic parameters

    // Warm-up
    for (final double a : new double[] {108}) {
      for (final double s : new double[] {120}) {
        for (final double n : new double[] {50, 100, 300}) {
          for (final double b2 : new double[] {0.5, 1, 2}) {
            for (int points = 3; points <= 20; points++) {
              Gaussian2DPeakResultHelper.getMLVarianceX(a, s, n, b2, true, points);
            }
          }
        }
      }
    }

    // Get average performance
    final double[] sum = new double[maxpoints + 1];
    final double[] sum2 = new double[sum.length];
    final long[] time = new long[sum.length];
    long count = 0;
    long count2 = 0;

    for (final double a : new double[] {108}) {
      for (final double s : new double[] {120}) {
        for (final double n : new double[] {50, 100, 300}) {
          for (final double b2 : new double[] {0.5, 1, 2}) {
            long min = Long.MAX_VALUE;
            for (int points = 3; points <= 20; points++) {
              long nanos = System.nanoTime();
              for (int i = 0; i < 1000; i++) {
                Gaussian2DPeakResultHelper.getMLVarianceX(a, s, n, b2, true, points);
              }
              nanos = time[points] = System.nanoTime() - nanos;
              if (min > nanos) {
                min = nanos;
              }
            }
            // Proportional weighting to the calculation that takes the longest
            count++;
            count2 += min;

            // Store relative performance
            final double factor = 1.0 / min;
            for (int points = 3; points <= 20; points++) {
              sum[points] += time[points] * factor;
              sum2[points] += time[points];
            }
          }
        }
      }
    }

    for (int points = minpoints; points <= maxpoints; points++) {
      logger.log(TestLevel.TEST_INFO,
          FunctionUtils.getSupplier("Points = %d, Av relative time = %f, Slow down factor = %f",
              points, sum[points] / count, sum2[points] / count2));
    }
  }

  @Test
  void canComputePixelAmplitude() {
    final float[] x = new float[] {0f, 0.1f, 0.3f, 0.5f, 0.7f, 1f};
    final float[] s = new float[] {0.8f, 1f, 1.5f, 2.2f};

    final float[] paramsf = new float[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    paramsf[Gaussian2DFunction.BACKGROUND] = 0;
    paramsf[Gaussian2DFunction.SIGNAL] = 105;

    final Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, 1, 1,
        GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);

    final SimpleRegression r = new SimpleRegression(false);

    for (final float tx : x) {
      for (final float ty : x) {
        for (final float sx : s) {
          for (final float sy : s) {
            paramsf[Gaussian2DFunction.X_POSITION] = tx;
            paramsf[Gaussian2DFunction.Y_POSITION] = ty;
            paramsf[Gaussian2DFunction.X_SD] = sx;
            paramsf[Gaussian2DFunction.Y_SD] = sy;

            // Get the answer using a single pixel image
            // Note the Gaussian2D functions set the centre of the pixel as 0,0 so offset
            final double[] params = SimpleArrayUtils.toDouble(paramsf);
            params[Gaussian2DFunction.X_POSITION] -= 0.5;
            params[Gaussian2DFunction.Y_POSITION] -= 0.5;
            f.initialise0(params);
            final double e = f.eval(0);

            final PSF psf = PsfHelper.create(PSFType.TWO_AXIS_GAUSSIAN_2D);
            final CalibrationWriter calibration = new CalibrationWriter();
            calibration.setCountPerPhoton(1);
            calibration.setIntensityUnit(IntensityUnit.PHOTON);
            calibration.setNmPerPixel(1);
            calibration.setDistanceUnit(DistanceUnit.PIXEL);
            final Gaussian2DPeakResultCalculator calc = Gaussian2DPeakResultHelper.create(psf,
                calibration,
                Gaussian2DPeakResultHelper.AMPLITUDE | Gaussian2DPeakResultHelper.PIXEL_AMPLITUDE);
            final double o1 = calc.getAmplitude(paramsf);
            final double o2 = calc.getPixelAmplitude(paramsf);

            // logger.fine(FunctionUtils.getSupplier("e=%f, o1=%f, o2=%f", e, o1, o2));
            Assertions.assertEquals(e, o2, 1e-3);
            r.addData(e, o1);
          }
        }
      }
    }

    // logger.fine(FunctionUtils.getSupplier("Regression: pixel amplitude vs amplitude = %f,
    // slope=%f, n=%d", r.getR(), r.getSlope(),
    // r.getN()));
    // The simple amplitude over estimates the actual pixel amplitude
    Assertions.assertTrue(r.getSlope() > 1);
  }

  @Test
  void canComputeCumulative() {
    Assertions.assertEquals(0, Gaussian2DPeakResultHelper.cumulative(0));
    Assertions.assertEquals(0.6827, Gaussian2DPeakResultHelper.cumulative(1), 1e-3);
    Assertions.assertEquals(0.9545, Gaussian2DPeakResultHelper.cumulative(2), 1e-3);
    Assertions.assertEquals(0.9974, Gaussian2DPeakResultHelper.cumulative(3), 1e-3);
    Assertions.assertTrue(1 == Gaussian2DPeakResultHelper.cumulative(Double.POSITIVE_INFINITY));
  }

  @Test
  void canComputeCumulative2DAndInverse() {
    Assertions.assertEquals(0, Gaussian2DPeakResultHelper.cumulative2D(0));
    Assertions.assertTrue(1 == Gaussian2DPeakResultHelper.cumulative2D(Double.POSITIVE_INFINITY));
    Assertions.assertEquals(0, Gaussian2DPeakResultHelper.inverseCumulative2D(0));
    Assertions
        .assertTrue(Double.POSITIVE_INFINITY == Gaussian2DPeakResultHelper.inverseCumulative2D(1));
    final DoubleDoubleBiPredicate predicate = Predicates.doublesAreClose(1e-8, 0);
    for (int i = 1; i <= 10; i++) {
      final double r = i / 10.0;
      final double p = Gaussian2DPeakResultHelper.cumulative2D(r);
      final double r2 = Gaussian2DPeakResultHelper.inverseCumulative2D(p);
      TestAssertions.assertTest(r, r2, predicate);
    }
  }

  @SeededTest
  void canComputeMeanSignalUsingR(RandomSeed seed) {
    final UniformRandomProvider rg = RngFactory.create(seed.get());

    final DoubleDoubleBiPredicate predicate = Predicates.doublesAreClose(1e-10, 0);
    for (int i = 0; i < 10; i++) {
      final double intensity = rg.nextDouble() * 100;
      final double sx = rg.nextDouble() * 2;
      final double sy = rg.nextDouble() * 2;
      final double r = rg.nextDouble() * 5;
      assertEquals(
          intensity * Gaussian2DPeakResultHelper.cumulative2D(r) / (Math.PI * r * r * sx * sy),
          Gaussian2DPeakResultHelper.getMeanSignalUsingR(intensity, sx, sy, r), predicate);

      // Test fixed versions verse dynamic
      assertEquals(Gaussian2DPeakResultHelper.getMeanSignalUsingR(intensity, sx, sy, 1),
          Gaussian2DPeakResultHelper.getMeanSignalUsingR1(intensity, sx, sy), predicate);
      assertEquals(Gaussian2DPeakResultHelper.getMeanSignalUsingR(intensity, sx, sy, 2),
          Gaussian2DPeakResultHelper.getMeanSignalUsingR2(intensity, sx, sy), predicate);
    }
  }

  private static void assertEquals(double expected, double observed,
      DoubleDoubleBiPredicate predicate) {
    TestAssertions.assertTest(expected, observed, predicate);
  }

  @SeededTest
  void canComputeMeanSignalUsingP(RandomSeed seed) {
    final UniformRandomProvider rg = RngFactory.create(seed.get());

    final DoubleDoubleBiPredicate predicate = Predicates.doublesAreClose(1e-10, 0);
    for (int i = 0; i < 10; i++) {
      final double intensity = rg.nextDouble() * 100;
      final double sx = rg.nextDouble() * 2;
      final double sy = rg.nextDouble() * 2;
      final double p = rg.nextDouble();
      double expected = intensity * p
          / (Math.PI * MathUtils.pow2(Gaussian2DPeakResultHelper.inverseCumulative2D(p)) * sx * sy);
      double observed = Gaussian2DPeakResultHelper.getMeanSignalUsingP(intensity, sx, sy, p);
      assertEquals(expected, observed, predicate);

      // Test fixed versions verse dynamic
      expected = Gaussian2DPeakResultHelper.getMeanSignalUsingP(intensity, sx, sy, 0.5);
      observed = Gaussian2DPeakResultHelper.getMeanSignalUsingP05(intensity, sx, sy);
      assertEquals(expected, observed, predicate);
    }
  }
}
