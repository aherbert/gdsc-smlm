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

package uk.ac.sussex.gdsc.smlm.fitting;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Logger;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.SharedStateContinuousSampler;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;
import uk.ac.sussex.gdsc.smlm.fitting.JumpDistanceAnalysis.JumpDistanceCumulFunction;
import uk.ac.sussex.gdsc.smlm.fitting.JumpDistanceAnalysis.JumpDistanceFunction;
import uk.ac.sussex.gdsc.smlm.fitting.JumpDistanceAnalysis.MixedJumpDistanceCumulFunction;
import uk.ac.sussex.gdsc.smlm.fitting.JumpDistanceAnalysis.MixedJumpDistanceFunction;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestSettings;
import uk.ac.sussex.gdsc.test.utils.functions.FunctionUtils;

@SuppressWarnings({"javadoc"})
public class JumpDistanceAnalysisTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(JumpDistanceAnalysisTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  // Based on the paper: Weimann, L., Ganzinger, K.A., McColl, J., Irvine, K.L., Davis, S.J.,
  // Gay, N.J., Bryant, C.E., Klenerman, D. (2013) A Quantitative Comparison of Single-Dye
  // Tracking Analysis Tools Using Monte Carlo Simulations. PLoS One 8, Issue 5, e64287
  //
  // This paper simulated tracks using 150 particle over 30 frames with a SNR of 6. The spots
  // were fit and then MSD or Jump Distance analysis performed. This means that each image
  // could have 150*29 jumps = 4350. They compute the fit using 750 trajectories (21750 jumps)
  // and repeat this 10 times to get a mean fit. They indicate success if
  // D is within 10% of the true value (the threshold for F is not explicitly stated
  // but appears to be around 20% or less). 2 population sample used D = 0.1, 0.02.
  // Good results were obtained when F(mobile) > 20%.

  // TODO - revise this so that fitting always works and then reduce the sample size down gradually
  // to see if there is a fail limit.

  // Test MLE fitting and histogram fitting separately.
  // Is MLE fitting worth doing. Can it be made better?

  final DoubleDoubleBiPredicate deltaD = TestHelper.doublesAreClose(0.1, 0);
  final DoubleDoubleBiPredicate deltaF = TestHelper.doublesAreClose(0.2, 0);

  // Used for testing single populations
  // Used for testing dual populations:
  // 15-fold, 5-fold, 3-fold difference between pairs
  // double[] D = new double[] { 0.2, 3, 1 };
  // 5-fold difference between pairs
  // double[] D = new double[] { 0.2, 1 };
  // For proteins with mass 823 and 347 kDa the
  // difference using predicted diffusion coefficients is 3:1
  static final double[] D = new double[] {3, 1};

  @Disabled("Commented out as this test always passes")
  @Test
  public void canIntegrateProbabilityToCumulativeWithSinglePopulation() {
    final JumpDistanceAnalysis jd = new JumpDistanceAnalysis();
    jd.setMinD(0);
    jd.setMinFraction(0);
    final SimpsonIntegrator si =
        new SimpsonIntegrator(1e-3, 1e-8, 2, SimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT);
    final DoubleDoubleBiPredicate equality = TestHelper.doublesAreClose(1e-2, 0);
    for (final double d : D) {
      final double[] params = new double[] {d};
      final JumpDistanceFunction fp = new JumpDistanceFunction(null, d);
      final JumpDistanceCumulFunction fc = new JumpDistanceCumulFunction(null, null, d);
      double x = d / 8;
      final UnivariateFunction func = new UnivariateFunction() {
        @Override
        public double value(double x) {
          return fp.evaluate(x, params);
        }
      };
      for (int i = 1; i < 10; i++, x *= 2) {
        final double e = fc.evaluate(x, params);
        // Integrate
        final double o = si.integrate(10000, func, 0, x);
        // logger.info(FunctionUtils.getSupplier("Integrate d=%.1f : x=%.1f, e=%f, o=%f, iter=%d,
        // eval=%d", d, x, e, o, si.getIterations(),
        // si.getEvaluations());
        TestAssertions.assertTest(e, o, equality,
            FunctionUtils.getSupplier("Failed to integrate: x=%g", x));
      }
    }
  }

  @Disabled("Commented out as this test always passes")
  @Test
  public void canIntegrateProbabilityToCumulativeWithMixedPopulation() {
    final JumpDistanceAnalysis jd = new JumpDistanceAnalysis();
    jd.setMinD(0);
    jd.setMinFraction(0);
    final SimpsonIntegrator si =
        new SimpsonIntegrator(1e-3, 1e-8, 2, SimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT);
    final DoubleDoubleBiPredicate equality = TestHelper.doublesAreClose(1e-2, 0);
    for (final double d : D) {
      for (final double f : new double[] {0, 0.1, 0.2, 0.4, 0.7, 0.9, 1}) {
        final double[] params = new double[] {f, d, 1 - f, d * 0.1};
        final MixedJumpDistanceFunction fp = new MixedJumpDistanceFunction(null, d, 2);
        final MixedJumpDistanceCumulFunction fc =
            new MixedJumpDistanceCumulFunction(null, null, d, 2);
        final UnivariateFunction func = new UnivariateFunction() {
          @Override
          public double value(double x) {
            return fp.evaluate(x, params);
          }
        };
        double x = d / 8;
        for (int i = 1; i < 10; i++, x *= 2) {
          final double e = fc.evaluate(x, params);
          // Integrate
          final double o = si.integrate(10000, func, 0, x);
          // logger.info(FunctionUtils.getSupplier("Integrate d=%.1f, f=%.1f : x=%.1f, e=%f, o=%f,
          // iter=%d, eval=%d", d, f, x, e, o,
          // si.getIterations(), si.getEvaluations());
          TestAssertions.assertTest(e, o, equality,
              FunctionUtils.getSupplier("Failed to integrate: x=%g", x));
        }
      }
    }
  }

  // @formatter:off
  @SeededTest
  public void canFitSinglePopulationMLE(RandomSeed seed)  { fitSinglePopulation(seed, true);  }
  @SeededTest
  public void canFitSinglePopulation(RandomSeed seed) { fitSinglePopulation(seed, false); }
  // @formatter:on

  private void fitSinglePopulation(RandomSeed seed, boolean mle) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeed());
    final String title = String.format("%s Single  ", (mle) ? "MLE" : "LSQ");
    AssertionError error = null;
    NEXT_D: for (final double d : D) {
      for (int samples = 500, k = 0; k < 6; samples *= 2, k++) {
        try {
          fit(rg, title, samples, 0, new double[] {d}, new double[] {1}, mle);
          error = null;
          continue NEXT_D;
        } catch (final AssertionError ex) {
          error = ex;
        }
      }
      if (error != null) {
        throw error;
      }
    }
  }

  // @formatter:off
  @SeededTest
  public void canFitDual0_1PopulationMLE(RandomSeed seed) { fitDualPopulation(seed, true,  0.1); }
  @SeededTest
  public void canFitDual0_1Population(RandomSeed seed)    { fitDualPopulation(seed, false, 0.1); }
  @SeededTest
  public void canFitDual0_2PopulationMLE(RandomSeed seed) { fitDualPopulation(seed, true,  0.2); }
  @SeededTest
  public void canFitDual0_2Population(RandomSeed seed)    { fitDualPopulation(seed, false, 0.2); }
  @SeededTest
  public void canFitDual0_3PopulationMLE(RandomSeed seed) { fitDualPopulation(seed, true,  0.3); }
  @SeededTest
  public void canFitDual0_3Population(RandomSeed seed)    { fitDualPopulation(seed, false, 0.3); }
  @SeededTest
  public void canFitDual0_4PopulationMLE(RandomSeed seed) { fitDualPopulation(seed, true,  0.4); }
  @SeededTest
  public void canFitDual0_4Population(RandomSeed seed)    { fitDualPopulation(seed, false, 0.4); }
  @SeededTest
  public void canFitDual0_5PopulationMLE(RandomSeed seed) { fitDualPopulation(seed, true,  0.5); }
  @SeededTest
  public void canFitDual0_5Population(RandomSeed seed)    { fitDualPopulation(seed, false, 0.5); }
  @SeededTest
  public void canFitDual0_6PopulationMLE(RandomSeed seed) { fitDualPopulation(seed, true,  0.6); }
  @SeededTest
  public void canFitDual0_6Population(RandomSeed seed)    { fitDualPopulation(seed, false, 0.6); }
  @SeededTest
  public void canFitDual0_7PopulationMLE(RandomSeed seed) { fitDualPopulation(seed, true,  0.7); }
  @SeededTest
  public void canFitDual0_7Population(RandomSeed seed)    { fitDualPopulation(seed, false, 0.7); }
  @SeededTest
  public void canFitDual0_8PopulationMLE(RandomSeed seed) { fitDualPopulation(seed, true,  0.8); }
  @SeededTest
  public void canFitDual0_8Population(RandomSeed seed)    { fitDualPopulation(seed, false, 0.8); }
  @SeededTest
  public void canFitDual0_9PopulationMLE(RandomSeed seed) { fitDualPopulation(seed, true,  0.9); }
  @SeededTest
  public void canFitDual0_9Population(RandomSeed seed)    { fitDualPopulation(seed, false, 0.9); }
  // @formatter:on

  private void fitDualPopulation(RandomSeed seed, boolean mle, double fraction) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MAXIMUM));
    final UniformRandomProvider rg = RngUtils.create(seed.getSeed());

    final String title = String.format("%s Dual=%.1f", (mle) ? "MLE" : "LSQ", fraction);
    AssertionError error = null;
    for (int i = 0; i < D.length; i++) {
      NEXT_D: for (int j = i + 1; j < D.length; j++) {
        for (int samples = 500, k = 0; k < 6; samples *= 2, k++) {
          try {
            fit(rg, title, samples, 0, new double[] {D[i], D[j]},
                new double[] {fraction, 1 - fraction}, mle);
            error = null;
            continue NEXT_D;
          } catch (final AssertionError ex) {
            error = ex;
          }
        }
        if (error != null) {
          throw error;
        }
      }
    }
  }

  private OutputStreamWriter out = null;

  /**
   * This is not actually a test but runs the fitting algorithm many times to collect benchmark data
   * to file.
   */
  @SeededTest
  public void canDoBenchmark(RandomSeed seed) {
    // Skip this as it is slow
    Assumptions.assumeTrue(false);
    final UniformRandomProvider rg = RngUtils.create(seed.getSeed());

    out = null;
    try {
      final FileOutputStream fos = new FileOutputStream("JumpDistanceAnalysisTest.dat");
      out = new OutputStreamWriter(fos, "UTF-8");

      // Run the fitting to produce benchmark data for a mixed population of 2
      final int n = 2;
      writeHeader(n);
      for (int repeat = 10; repeat-- > 0;) {
        resetData();
        for (final boolean mle : new boolean[] {true, false}) {
          for (int f = 1; f <= 9; f++) {
            final double fraction = f / 10.0;
            final String title = String.format("%s Dual=%.1f", (mle) ? "MLE" : "LSQ", fraction);
            for (int samples = 500, k = 0; k < 6; samples *= 2, k++) {
              for (int i = 0; i < D.length; i++) {
                for (int j = i + 1; j < D.length; j++) {
                  try {
                    fit(rg, title, samples, 0, new double[] {D[i], D[j]},
                        new double[] {fraction, 1 - fraction}, mle);
                  } catch (final AssertionError ex) {
                    // Carry on with the benchmarking
                  }
                  // If the fit had the correct N then no need to repeat
                  if (fitN == n) {
                    continue;
                  }
                  try {
                    fit(rg, title + " Fixed", samples, n, new double[] {D[i], D[j]},
                        new double[] {fraction, 1 - fraction}, mle);
                  } catch (final AssertionError ex) {
                    // Carry on with the benchmarking
                  }
                }
              }
            }
          }
        }
      }
    } catch (final Exception ex) {
      throw new AssertionError("Failed to complete benchmark", ex);
    } finally {
      closeOutput();
    }
  }

  private void closeOutput() {
    if (out == null) {
      return;
    }
    try {
      out.close();
    } catch (final Exception ex) {
      // Ignore exception
    } finally {
      out = null;
    }
  }

  // Store the fitted N to allow repeat in the benchmark with fixed N if necessary
  int fitN = 0;

  private void fit(UniformRandomProvider rg, String title, int samples, int n, double[] dc,
      double[] fraction, boolean mle) {
    // Used for testing
    // @formatter:off
    //if (!mle) return;
    //if (mle) return;
    // For easy mixed populations
    //if (f.length == 2 && Math.min(f[0],f[1])/(f[0]+f[1]) <= 0.2) return;
    //n = 2;
    // @formatter:on

    JumpDistanceAnalysis.sort(dc, fraction);
    final double[] jumpDistances = createData(rg, samples, dc, fraction);
    final JumpDistanceAnalysis jd = new JumpDistanceAnalysis();
    jd.setFitRestarts(3);
    double[][] fit;
    if (n == 0) {
      jd.setMinFraction(0.05);
      jd.setMinDifference(2);
      jd.setMaxN(10);
      fit = (mle) ? jd.fitJumpDistancesMle(jumpDistances) : jd.fitJumpDistances(jumpDistances);
    } else {
      // No validation
      jd.setMinFraction(0);
      jd.setMinDifference(0);
      fit =
          (mle) ? jd.fitJumpDistancesMle(jumpDistances, n) : jd.fitJumpDistances(jumpDistances, n);
    }
    final double[] fitD = (fit == null) ? new double[0] : fit[0];
    final double[] fitF = (fit == null) ? new double[0] : fit[1];

    // Record results to file
    if (out != null) {
      writeResult(title, sample.dc, sample.fraction, samples, n, dc, fraction, mle, fitD, fitF);
    }

    fitN = fitD.length;
    AssertionError error = null;
    try {
      Assertions.assertEquals(dc.length, fitD.length, "Failed to fit n");
      TestAssertions.assertArrayTest(dc, fitD, deltaD, "Failed to fit d");
      TestAssertions.assertArrayTest(fraction, fitF, deltaF, "Failed to fit f");
    } catch (final AssertionError ex) {
      error = ex;
    } finally {
      final double[] e1 = getPercentError(dc, fitD);
      final double[] e2 = getPercentError(fraction, fitF);
      logger.info(
          FunctionUtils.getSupplier("%s %s N=%d sample=%d, n=%d : %s = %s [%s] : %s = %s [%s]",
              (error == null) ? "+++ Pass" : "--- Fail", title, dc.length, samples, n, toString(dc),
              toString(fitD), toString(e1), toString(fraction), toString(fitF), toString(e2)));
      if (error != null) {
        throw error;
      }
    }
  }

  private void writeHeader(int size) {
    final StringBuilder sb = new StringBuilder("title");
    sb.append('\t').append("repeat");
    for (int i = 0; i < size; i++) {
      sb.append('\t').append("D").append(i);
    }
    for (int i = 0; i < size; i++) {
      sb.append('\t').append("F").append(i);
    }
    sb.append('\t').append("samples");
    sb.append('\t').append("mle");
    sb.append('\t').append("n");
    sb.append('\t').append("size");
    sb.append('\t').append("fsize");
    for (int i = 0; i < size; i++) {
      sb.append('\t').append("d").append(i);
    }
    for (int i = 0; i < size; i++) {
      sb.append('\t').append("fd").append(i);
    }
    for (int i = 0; i < size; i++) {
      sb.append('\t').append("ed").append(i);
    }
    for (int i = 0; i < size; i++) {
      sb.append('\t').append("f").append(i);
    }
    for (int i = 0; i < size; i++) {
      sb.append('\t').append("ff").append(i);
    }
    for (int i = 0; i < size; i++) {
      sb.append('\t').append("ef").append(i);
    }
    sb.append(System.lineSeparator());
    try {
      out.write(sb.toString());
    } catch (final IOException ex) {
      throw new AssertionError("Failed to write result to file", ex);
    }
  }

  private void writeResult(String title, double[] actualD, double[] actualF, int samples, int n,
      double[] dc, double[] fraction, boolean mle, double[] fitDc, double[] fitFraction) {
    final int size = dc.length;
    final int fitSize = fitDc.length;
    // Pad results if they are too small
    if (fitSize < size) {
      fitDc = Arrays.copyOf(fitDc, size);
      fitFraction = Arrays.copyOf(fitFraction, size);
    }
    final double[] ed = getRelativeError(dc, fitDc);
    final double[] ef = getRelativeError(fraction, fitFraction);

    final StringBuilder sb = new StringBuilder(title);
    sb.append('\t').append(repeat);
    for (int i = 0; i < size; i++) {
      sb.append('\t').append(actualD[i]);
    }
    for (int i = 0; i < size; i++) {
      sb.append('\t').append(actualF[i]);
    }
    sb.append('\t').append(samples);
    sb.append('\t').append(mle);
    sb.append('\t').append(n);
    sb.append('\t').append(size);
    sb.append('\t').append(fitSize);
    for (int i = 0; i < size; i++) {
      sb.append('\t').append(dc[i]);
    }
    for (int i = 0; i < size; i++) {
      sb.append('\t').append(fitDc[i]);
    }
    for (int i = 0; i < size; i++) {
      sb.append('\t').append(ed[i]);
    }
    for (int i = 0; i < size; i++) {
      sb.append('\t').append(fraction[i]);
    }
    for (int i = 0; i < size; i++) {
      sb.append('\t').append(fitFraction[i]);
    }
    for (int i = 0; i < size; i++) {
      sb.append('\t').append(ef[i]);
    }
    sb.append(System.lineSeparator());
    try {
      out.write(sb.toString());
    } catch (final IOException ex) {
      throw new AssertionError("Failed to write result to file", ex);
    }
  }

  private static double[] getPercentError(double[] exp, double[] obs) {
    final double[] error = new double[Math.min(exp.length, obs.length)];
    for (int i = 0; i < error.length; i++) {
      error[i] = 100.0 * (obs[i] - exp[i]) / exp[i];
    }
    return error;
  }

  private static double[] getRelativeError(double[] exp, double[] obs) {
    final double[] error = new double[Math.min(exp.length, obs.length)];
    for (int i = 0; i < error.length; i++) {
      // As per the Weimann Plos One paper
      error[i] = Math.abs(obs[i] - exp[i]) / exp[i];
      // Use the relative error from the largest value
      // error[i] = uk.ac.sussex.gdsc.smlm.utils.DoubleEquality.relativeError(o[i], e[i]);
    }
    return error;
  }

  private static String toString(double[] data) {
    if (data.length == 0) {
      return "";
    }
    if (data.length == 1) {
      return format(data[0]);
    }
    final StringBuilder sb = new StringBuilder();
    sb.append(format(data[0]));
    for (int i = 1; i < data.length; i++) {
      sb.append(',').append(format(data[i]));
    }
    return sb.toString();
  }

  private static String format(double value) {
    return String.format("%.3f", value);
  }

  class DataSample {
    double[] dc;
    double[] fraction;
    double[] sigma;
    double[] data;
    int[] sample;

    DataSample(double[] dc, double[] fraction) {
      this.dc = dc.clone();

      // Convert diffusion co-efficient into the standard deviation for the random move in each
      // dimension.
      // For 1D diffusion: sigma^2 = 2D
      // sigma = sqrt(2D)
      // See: https://en.wikipedia.org/wiki/Brownian_motion#Einstein.27s_theory
      sigma = new double[dc.length];
      double sum = 0;
      for (int i = 0; i < sigma.length; i++) {
        sigma[i] = Math.sqrt(2 * dc[i]);
        sum += fraction[i];
      }
      this.fraction = new double[fraction.length];
      for (int i = 0; i < fraction.length; i++) {
        this.fraction[i] = fraction[i] / sum;
      }
    }

    void add(double[] data2, int[] sample2) {
      if (data == null) {
        data = data2;
        sample = sample2;
        return;
      }
      final int size = data.length;
      final int newSize = size + data2.length;
      data = Arrays.copyOf(data, newSize);
      sample = Arrays.copyOf(sample, newSize);
      System.arraycopy(data2, 0, data, size, data2.length);
      System.arraycopy(sample2, 0, sample, size, sample2.length);
    }

    int getSize() {
      return (data == null) ? 0 : data.length;
    }

    double[][] getSample(UniformRandomProvider random, int size) {
      if (size > getSize()) {
        final SharedStateContinuousSampler gs = SamplerUtils.createGaussianSampler(random, 0, 1);
        final int extra = size - getSize();

        // Get cumulative fraction
        final double[] c = new double[fraction.length];
        double sum = 0;
        for (int i = 0; i < fraction.length; i++) {
          sum += fraction[i];
          c[i] = sum;
        }

        final double[] data = new double[extra];

        // Pick the population using the fraction.
        // Do this before sampling since the nextGaussian function computes random variables
        // in pairs so we want to process all the same sample together
        final int[] sample = new int[extra];
        if (c.length > 1) {
          for (int i = 0; i < data.length; i++) {
            sample[i] = pick(c, random.nextDouble());
          }
        }
        Arrays.sort(sample);

        for (int i = 0; i < data.length; i++) {
          // Pick the population using the fraction
          final int j = sample[i];
          // Get the x/y shifts
          final double x = gs.sample() * sigma[j];
          final double y = gs.sample() * sigma[j];
          // Get the squared jump distance
          data[i] = x * x + y * y;
        }
        add(data, sample);
      }

      // Build the sample data and return the D and fractions
      final double[] data = Arrays.copyOf(this.data, size);
      final double[] d = new double[this.dc.length];
      final double[] f = new double[d.length];
      for (int i = 0; i < size; i++) {
        final int j = sample[i];
        d[j] += data[i];
        f[j]++;
      }
      for (int i = 0; i < d.length; i++) {
        // 4D = MSD
        // D = MSD / 4
        d[i] = (d[i] / f[i]) / 4;
        f[i] /= size;
      }

      return new double[][] {data, d, f};
    }

    @Override
    public boolean equals(Object obj) {
      if (!(obj instanceof DataSample)) {
        return super.equals(obj);
      }
      final DataSample that = (DataSample) obj;
      if (that.dc.length != this.dc.length) {
        return false;
      }
      for (int i = dc.length; i-- > 0;) {
        if (that.dc[i] != this.dc[i]) {
          return false;
        }
        if (that.fraction[i] != this.fraction[i]) {
          return false;
        }
      }
      return true;
    }

    @Override
    public int hashCode() {
      int hash = 1;
      for (int i = dc.length; i-- > 0;) {
        hash = hash * 31 + Double.hashCode(dc[i]);
        hash = hash * 31 + Double.hashCode(fraction[i]);
      }
      return hash;
    }
  }

  static ArrayList<DataSample> samples = new ArrayList<>();
  DataSample sample = null;
  private int repeat = 0;

  private void resetData() {
    samples.clear();
    repeat++;
  }

  /**
   * Create random jump distances.
   *
   * @param n Number of jump distances
   * @param dc Diffusion rate (should be ascending order of magnitude)
   * @param fraction Fraction of population (will be updated with the actual fractions, normalised
   *        to sum to 1)
   * @return The jump distances
   */
  private double[] createData(UniformRandomProvider rg, int n, double[] dc, double[] fraction) {
    // Cache the data so that if we run a second test with
    // the same d and f we use the same data
    sample = new DataSample(dc, fraction);
    final int index = samples.indexOf(sample);
    if (index != -1) {
      sample = samples.get(index);
    } else {
      samples.add(sample);
    }

    final double[][] dataSample = sample.getSample(rg, n);
    final double[] data = dataSample[0];
    final double[] dc2 = dataSample[1];
    final double[] f2 = dataSample[2];

    // Update with the real values
    for (int i = 0; i < dc.length; i++) {
      dc[i] = dc2[i];
      fraction[i] = f2[i];
    }

    // Debug
    // uk.ac.sussex.gdsc.core.utils.StoredDataStatistics stats = new
    // uk.ac.sussex.gdsc.core.utils.StoredDataStatistics(data);
    // uk.ac.sussex.gdsc.core.ij.new HistogramPlotBuilder(
    // "MSD",
    // stats,
    // "MSD",
    // 0,
    // 0,
    // n / 10,
    // true,
    // String.format("%s : %s : u=%f, sd=%f", toString(d), toString(f), stats.getMean(),
    // stats.getStandardDeviation()));

    return data;
  }

  private static int pick(double[] fraction, double nextDouble) {
    for (int i = 0; i < fraction.length; i++) {
      if (nextDouble < fraction[i]) {
        return i;
      }
    }
    return fraction.length - 1;
  }
}
