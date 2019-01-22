/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.function;

import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TestSettings;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.function.Supplier;
import java.util.logging.Level;
import java.util.logging.Logger;

@SuppressWarnings({"javadoc"})
public class PoissonGaussianConvolutionFunctionTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(PoissonGaussianConvolutionFunctionTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  double[] gain = PoissonGaussianFunctionTest.gain;
  double[] photons = PoissonGaussianFunctionTest.photons;
  double[] noise = PoissonGaussianFunctionTest.noise;

  @Test
  public void cumulativeProbabilityIsOneWithPdf() {
    for (final double g : gain) {
      for (final double p : photons) {
        for (final double s : noise) {
          cumulativeProbabilityIsOne(g, p, s, false);
        }
      }
    }
  }

  @Test
  public void cumulativeProbabilityIsOneWithPmf() {
    for (final double g : gain) {
      for (final double p : photons) {
        for (final double s : noise) {
          cumulativeProbabilityIsOne(g, p, s, true);
        }
      }
    }
  }

  @Test
  public void probabilityMatchesLogProbabilityWithPdf() {
    for (final double g : gain) {
      for (final double p : photons) {
        for (final double s : noise) {
          probabilityMatchesLogProbability(g, p, s, false);
        }
      }
    }
  }

  @Test
  public void probabilityMatchesLogProbabilityWithPmf() {
    for (final double g : gain) {
      for (final double p : photons) {
        for (final double s : noise) {
          probabilityMatchesLogProbability(g, p, s, true);
        }
      }
    }
  }

  private static void cumulativeProbabilityIsOne(final double gain, final double mu,
      final double sd, boolean computePmf) {
    final double p2 = cumulativeProbability(gain, mu, sd, computePmf);
    Assertions.assertEquals(1, p2, 0.02,
        () -> String.format("g=%f, mu=%f, s=%f, erf=%b", gain, mu, sd, computePmf));
  }

  private static double cumulativeProbability(final double gain, final double mu, final double sd,
      boolean computePmf) {
    // Note: The input s parameter is pre-gain.
    final PoissonGaussianConvolutionFunction f =
        PoissonGaussianConvolutionFunction.createWithStandardDeviation(1.0 / gain, sd * gain);
    f.setComputePmf(computePmf);

    // final PoissonGaussianConvolutionFunction f2 =
    // PoissonGaussianConvolutionFunction.createWithStandardDeviation(1.0 / gain, mu*gain, s *
    // gain);
    // f2.setUsePicardApproximation(usePicard);

    double pvalue = 0;
    int min = 1;
    int max = 0;

    // Note: The input mu parameter is pre-gain.
    final double e = mu;

    // Evaluate an initial range.
    // Gaussian should have >99% within +/- s
    // Poisson will have mean mu with a variance mu.
    // At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
    if (mu > 0) {
      final int[] range = PoissonGaussianFunctionTest.getRange(gain, mu, sd);
      min = range[0];
      max = range[1];
      for (int x = min; x <= max; x++) {
        final double pp = f.likelihood(x, e);
        // logger.fine(FunctionUtils.getSupplier("x=%d, p=%g", x, pp);
        // logger.fine(FunctionUtils.getSupplier("x=%d, p=%f %f", x, pp, f2.probability(x));
        pvalue += pp;
      }
      // if (p > 1.01)
      // Assertions.fail("P > 1: " + p);
    }

    // We have most of the likelihood density.
    // Now keep evaluating up and down until no difference
    final double changeTolerance = 1e-6;
    for (int x = min - 1;; x--) {
      min = x;
      final double pp = f.likelihood(x, e);
      // logger.fine(FunctionUtils.getSupplier("x=%d, p=%g", x, pp);
      pvalue += pp;
      if (pp == 0 || pp / pvalue < changeTolerance) {
        break;
      }
    }
    for (int x = max + 1;; x++) {
      max = x;
      final double pp = f.likelihood(x, e);
      // logger.fine(FunctionUtils.getSupplier("x=%d, p=%g", x, pp);
      pvalue += pp;
      if (pp == 0 || pp / pvalue < changeTolerance) {
        break;
      }
    }

    double p2 = pvalue;
    if (!computePmf) {
      // Do a formal integration if the PDF
      // if (p < 0.98 || p > 1.02)
      // logger.fine(FunctionUtils.getSupplier("g=%f, mu=%f, s=%f p=%f", gain, mu, s, p);
      final UnivariateIntegrator in =
          new SimpsonIntegrator(1e-4, 1e-6, 4, SimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT);
      p2 = in.integrate(Integer.MAX_VALUE, new UnivariateFunction() {
        @Override
        public double value(double x) {
          return f.likelihood(x, e);
        }
      }, min, max);
    }

    if (p2 < 0.98 || p2 > 1.02) {
      logger.log(TestLogUtils.getRecord(Level.INFO, "g=%f, mu=%f, s=%f p=%f  %f", gain, mu, sd,
          pvalue, p2));
    }

    return p2;
  }

  private static void probabilityMatchesLogProbability(final double gain, double mu,
      final double sd, boolean computePmf) {
    // Note: The input s parameter is pre-gain.
    final PoissonGaussianConvolutionFunction f =
        PoissonGaussianConvolutionFunction.createWithStandardDeviation(1.0 / gain, sd * gain);
    f.setComputePmf(computePmf);

    // Evaluate an initial range.
    // Gaussian should have >99% within +/- s
    // Poisson will have mean mu with a variance mu.
    // At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
    final int[] range = PoissonGaussianFunctionTest.getRange(gain, mu, sd);
    final int min = range[0];
    final int max = range[1];
    // Note: The input mu parameter is pre-gain.
    final double e = mu;
    final Supplier<String> msg =
        () -> String.format("g=%f, mu=%f, s=%f, erf=%b", gain, mu, sd, computePmf);
    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-3, 0);
    for (int x = min; x <= max; x++) {
      final double p = f.likelihood(x, e);
      if (p == 0) {
        continue;
      }
      final double logP = f.logLikelihood(x, e);
      TestAssertions.assertTest(Math.log(p), logP, predicate, msg);
    }
  }

  @SpeedTag
  @SeededTest
  public void pdfFasterThanPmf(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    // Realistic CCD parameters for speed test
    final double s = 7.16;
    final double g = 3.1;

    final PoissonGaussianConvolutionFunction f1 =
        PoissonGaussianConvolutionFunction.createWithStandardDeviation(1 / g, s);
    f1.setComputePmf(true);

    final PoissonGaussianConvolutionFunction f2 =
        PoissonGaussianConvolutionFunction.createWithStandardDeviation(1 / g, s);
    f2.setComputePmf(false);

    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());

    // Generate realistic data from the probability mass function
    final double[][] samples = new double[photons.length][];
    for (int j = 0; j < photons.length; j++) {
      final int start = (int) (4 * -s);
      int mu = start;
      final StoredDataStatistics stats = new StoredDataStatistics();
      while (stats.getSum() < 0.995) {
        final double p = f1.likelihood(mu, photons[j]);
        stats.add(p);
        if (mu > 10 && p / stats.getSum() < 1e-6) {
          break;
        }
        mu++;
      }

      // Generate cumulative probability
      final double[] data = stats.getValues();
      for (int i = 1; i < data.length; i++) {
        data[i] += data[i - 1];
      }
      // Normalise
      for (int i = 0, end = data.length - 1; i < data.length; i++) {
        data[i] /= data[end];
      }

      // Sample
      final double[] sample = new double[1000];
      for (int i = 0; i < sample.length; i++) {
        final double p = rg.nextDouble();
        int x = 0;
        while (x < data.length && data[x] < p) {
          x++;
        }
        sample[i] = start + x;
      }
      samples[j] = sample;
    }

    // Warm-up
    run(f1, samples, photons);
    run(f2, samples, photons);

    long t1 = 0;
    for (int i = 0; i < 5; i++) {
      t1 += run(f1, samples, photons);
    }

    long t2 = 0;
    for (int i = 0; i < 5; i++) {
      t2 += run(f2, samples, photons);
    }

    logger.log(TestLogUtils.getTimingRecord("cdf", t1, "pdf", t2));
  }

  private static long run(PoissonGaussianConvolutionFunction func, double[][] samples,
      double[] photons) {
    final long start = System.nanoTime();
    for (int j = 0; j < photons.length; j++) {
      final double p = photons[j];
      for (final double x : samples[j]) {
        func.likelihood(x, p);
      }
    }
    return System.nanoTime() - start;
  }
}
