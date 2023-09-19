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

package uk.ac.sussex.gdsc.smlm.fitting;

import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.Supplier;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.AliasMethodDiscreteSampler;
import org.apache.commons.rng.sampling.distribution.SharedStateDiscreteSampler;
import org.apache.commons.rng.sampling.distribution.ZigguratSampler;
import org.apache.commons.statistics.distribution.ExponentialDistribution;
import org.apache.commons.statistics.distribution.GeometricDistribution;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.ValueSource;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.smlm.fitting.ResidenceTimeFitting.FitResult;
import uk.ac.sussex.gdsc.smlm.fitting.ResidenceTimeFitting.Function1;
import uk.ac.sussex.gdsc.smlm.fitting.ResidenceTimeFitting.Function1T;
import uk.ac.sussex.gdsc.smlm.fitting.ResidenceTimeFitting.Function2;
import uk.ac.sussex.gdsc.smlm.fitting.ResidenceTimeFitting.Function2T;
import uk.ac.sussex.gdsc.smlm.fitting.ResidenceTimeFitting.Function3;
import uk.ac.sussex.gdsc.smlm.fitting.ResidenceTimeFitting.Function3T;
import uk.ac.sussex.gdsc.smlm.fitting.ResidenceTimeFitting.Model;
import uk.ac.sussex.gdsc.test.api.Predicates;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngFactory;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogging.TestLevel;
import uk.ac.sussex.gdsc.test.utils.TestSettings;
import uk.ac.sussex.gdsc.test.utils.functions.FormatSupplier;

/**
 * Tests for {@link ResidenceTimeFitting}.
 */
@SuppressWarnings({"javadoc"})
class ResidenceTimeFittingTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(ResidenceTimeFittingTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
    samples.clear();
  }

  @Test
  void testCreateThrows() {
    Assertions.assertThrows(IllegalArgumentException.class,
        () -> ResidenceTimeFitting.of(0, new int[] {1, 1}));
    Assertions.assertThrows(IllegalArgumentException.class,
        () -> ResidenceTimeFitting.of(1, new int[] {0}));
    Assertions.assertThrows(IllegalArgumentException.class,
        () -> ResidenceTimeFitting.of(1, new int[] {0, 0}));
    Assertions.assertThrows(IllegalArgumentException.class,
        () -> ResidenceTimeFitting.of(1, new int[] {1, 0}));
    Assertions.assertThrows(IllegalArgumentException.class,
        () -> ResidenceTimeFitting.of(1, new int[] {0, 1}));
    Assertions.assertThrows(IllegalArgumentException.class,
        () -> ResidenceTimeFitting.of(1, new int[] {1, 1, -1}));
    Assertions.assertThrows(IllegalArgumentException.class,
        () -> ResidenceTimeFitting.of(1, new int[] {1 << 30, 1 << 30}));
  }

  @Test
  void testFitThrows() {
    // OK with two non-zero values
    final ResidenceTimeFitting rta = ResidenceTimeFitting.of(1, new int[] {6, 1});
    Assertions.assertThrows(IllegalArgumentException.class, () -> rta.fit(4));
  }

  @Test
  void testLogging() {
    final ResidenceTimeFitting rta = ResidenceTimeFitting.of(1, new int[] {16, 8, 4, 2, 1});
    final Logger logger = new Logger("Test", null) {
      @Override
      public void info(Supplier<String> msg) {
        Assertions.assertNotNull(msg.get());
      }
    };
    logger.setLevel(Level.INFO);
    for (int i = 1; i <= 3; i++) {
      final Pair<FitResult, Model> ra = rta.fit(i);
      final Pair<FitResult, Model> rb = rta.fit(i, logger);
      Assertions.assertEquals(ra.getLeft().getLogLikelihood(), rb.getLeft().getLogLikelihood());
    }
  }

  @Test
  void testMaxEval() {
    final ResidenceTimeFitting rta = ResidenceTimeFitting.of(1, new int[] {16, 8, 4, 2, 1});
    // Not enough evaluations to fit anything
    final MaxEval eval1 = new MaxEval(3);
    // Not enough evaluations for the n=2 Powell optimizer
    // This forces use of the CMAES optimizer
    final MaxEval eval2 = new MaxEval(185);
    for (int i = 1; i <= 2; i++) {
      Assertions.assertNull(rta.fit(i, eval1));
      final Pair<FitResult, Model> ra = rta.fit(i, eval2);
      Assertions.assertNotNull(ra);
    }
  }

  @ParameterizedTest
  @CsvSource({"1, 1", "1, 2", "1, 0.5", "0.5, 1", "0.5, 2", "0.5, 0.5", "0.75, 0.5", "0.75, 1.5",
      "1.5, 2"})
  void testRate(double k, double resolution) {
    // This is the computed assuming a geometric distribution histogram.
    final double mean = 1 / k;
    final ExponentialDistribution exp = ExponentialDistribution.of(mean);
    // 99.99 percentile of exponential is 9.21
    final double N = 1_000_000_000 * resolution;
    final int[] n = new int[(int) (12 / resolution)];
    for (int i = 0; i < n.length; i++) {
      n[i] = (int) Math.round(N * (1 - exp.cumulativeProbability((i + 1) * resolution)));
    }
    Assertions.assertEquals(k, ResidenceTimeFitting.of(resolution, n).getRate(), k * 0.05);
  }

  @ParameterizedTest
  @ValueSource(doubles = {0.1, 3})
  void testModel1Default(double k0) {
    final ExponentialDistribution d = ExponentialDistribution.of(1 / k0);
    final Model m = new Model() {
      @Override
      public int getSize() {
        return 1;
      }

      @Override
      public double getRate(int index) {
        return k0;
      }

      @Override
      public double getFraction(int index) {
        return 1;
      }
    };
    Assertions.assertEquals(1 / k0, m.getResidenceTime(0));
    Assertions.assertEquals(Double.POSITIVE_INFINITY, m.getUpperBound());
    final DoubleDoubleBiPredicate test = Predicates.doublesAreRelativelyClose(1e-10);
    DoubleStream.iterate(0, x -> x + 0.5).limit(10).forEach(t -> {
      final double e = d.survivalProbability(t);
      TestAssertions.assertTest(e, m.sf(t), test, () -> "sf: " + t);
    });
  }

  @ParameterizedTest
  @CsvSource({"0.1, Infinity", "0.3, Infinity", "0.1, 20", "0.1, 30"})
  void testModel1DefaultT(double k0, double upper) {
    final ExponentialDistribution d = ExponentialDistribution.of(1 / k0);
    final Model m = new Model() {
      @Override
      public int getSize() {
        return 1;
      }

      @Override
      public double getRate(int index) {
        return k0;
      }

      @Override
      public double getFraction(int index) {
        return 1;
      }

      @Override
      public double getUpperBound() {
        return upper;
      }
    };
    Assertions.assertEquals(1 / k0, m.getResidenceTime(0));
    Assertions.assertEquals(0, m.sf(upper), () -> "sf: " + upper);
    Assertions.assertEquals(0, m.sf(upper * 2), () -> "sf: " + upper * 2);
    final DoubleDoubleBiPredicate test = Predicates.doublesAreRelativelyClose(1e-10);
    final double cdfU = d.cumulativeProbability(upper);
    // This iteration (x < 5) is always below upper
    DoubleStream.iterate(0, x -> x + 0.5).limit(10).forEach(t -> {
      final double e = 1 - d.cumulativeProbability(t) / cdfU;
      TestAssertions.assertTest(e, m.sf(t), test, () -> "sf: " + t);
    });
  }

  @ParameterizedTest
  @CsvSource({"0.1, 3, 0.5", "2, 4, 0.25"})
  void testModel2Default(double k0, double k1, double f0) {
    final ExponentialDistribution d0 = ExponentialDistribution.of(1 / k0);
    final ExponentialDistribution d1 = ExponentialDistribution.of(1 / k1);
    final double f1 = 1 - f0;
    final Model m = new Model() {
      @Override
      public int getSize() {
        return 2;
      }

      @Override
      public double getRate(int index) {
        return index == 0 ? k0 : k1;
      }

      @Override
      public double getFraction(int index) {
        return index == 0 ? f0 : f1;
      }
    };
    Assertions.assertEquals(1 / k0, m.getResidenceTime(0));
    Assertions.assertEquals(1 / k1, m.getResidenceTime(1));
    Assertions.assertEquals(Double.POSITIVE_INFINITY, m.getUpperBound());
    final DoubleDoubleBiPredicate test = Predicates.doublesAreRelativelyClose(1e-10);
    DoubleStream.iterate(0, x -> x + 0.5).limit(10).forEach(t -> {
      final double e = f0 * d0.survivalProbability(t) + f1 * d1.survivalProbability(t);
      TestAssertions.assertTest(e, m.sf(t), test, () -> "sf: " + t);
    });
  }

  @ParameterizedTest
  @CsvSource({"0.1, 3, 0.5, Infinity", "2, 4, 0.25, 7", "2, 4, 0.25, 5",})
  void testModel2DefaultT(double k0, double k1, double f0, double upper) {
    final ExponentialDistribution d0 = ExponentialDistribution.of(1 / k0);
    final ExponentialDistribution d1 = ExponentialDistribution.of(1 / k1);
    final double f1 = 1 - f0;
    final Model m = new Model() {
      @Override
      public int getSize() {
        return 2;
      }

      @Override
      public double getRate(int index) {
        return index == 0 ? k0 : k1;
      }

      @Override
      public double getFraction(int index) {
        return index == 0 ? f0 : f1;
      }

      @Override
      public double getUpperBound() {
        return upper;
      }
    };
    Assertions.assertEquals(1 / k0, m.getResidenceTime(0));
    Assertions.assertEquals(1 / k1, m.getResidenceTime(1));
    Assertions.assertEquals(0, m.sf(upper), () -> "sf: " + upper);
    Assertions.assertEquals(0, m.sf(upper * 2), () -> "sf: " + upper * 2);
    final DoubleDoubleBiPredicate test = Predicates.doublesAreRelativelyClose(1e-10);
    final double cdfU0 = d0.cumulativeProbability(upper);
    final double cdfU1 = d1.cumulativeProbability(upper);
    // This iteration (x < 5) is always below upper
    DoubleStream.iterate(0, x -> x + 0.5).limit(10).forEach(t -> {
      final double e = f0 * (1 - d0.cumulativeProbability(t) / cdfU0)
          + f1 * (1 - d1.cumulativeProbability(t) / cdfU1);
      TestAssertions.assertTest(e, m.sf(t), test, () -> "sf: " + t);
    });
  }

  @ParameterizedTest
  @CsvSource({"0.1, 3, 5, 0.5, 0.25", "2, 4, 10, 0.25, 0.125"})
  void testModel3Default(double k0, double k1, double k2, double f0, double f1) {
    final ExponentialDistribution d0 = ExponentialDistribution.of(1 / k0);
    final ExponentialDistribution d1 = ExponentialDistribution.of(1 / k1);
    final ExponentialDistribution d2 = ExponentialDistribution.of(1 / k2);
    final double f2 = 1 - (f0 + f1);
    final Model m = new Model() {
      @Override
      public int getSize() {
        return 3;
      }

      @Override
      public double getRate(int index) {
        switch (index) {
          case 0:
            return k0;
          case 1:
            return k1;
          default:
            return k2;
        }
      }

      @Override
      public double getFraction(int index) {
        switch (index) {
          case 0:
            return f0;
          case 1:
            return f1;
          default:
            return f2;
        }
      }
    };
    Assertions.assertEquals(1 / k0, m.getResidenceTime(0));
    Assertions.assertEquals(1 / k1, m.getResidenceTime(1));
    Assertions.assertEquals(1 / k2, m.getResidenceTime(2));
    Assertions.assertEquals(Double.POSITIVE_INFINITY, m.getUpperBound());
    final DoubleDoubleBiPredicate test = Predicates.doublesAreRelativelyClose(1e-10);
    DoubleStream.iterate(0, x -> x + 0.5).limit(10).forEach(t -> {
      final double e = f0 * d0.survivalProbability(t) + f1 * d1.survivalProbability(t)
          + f2 * d2.survivalProbability(t);
      TestAssertions.assertTest(e, m.sf(t), test, () -> "sf: " + t);
    });
  }

  @ParameterizedTest
  @CsvSource({"0.1, 3, 5, 0.5, 0.25, Infinity", "2, 4, 10, 0.25, 0.125, 7",
      "2, 4, 10, 0.25, 0.125, 5",})
  void testModel3DefaultT(double k0, double k1, double k2, double f0, double f1, double upper) {
    final ExponentialDistribution d0 = ExponentialDistribution.of(1 / k0);
    final ExponentialDistribution d1 = ExponentialDistribution.of(1 / k1);
    final ExponentialDistribution d2 = ExponentialDistribution.of(1 / k2);
    final double f2 = 1 - (f0 + f1);
    final Model m = new Model() {
      @Override
      public int getSize() {
        return 3;
      }

      @Override
      public double getRate(int index) {
        switch (index) {
          case 0:
            return k0;
          case 1:
            return k1;
          default:
            return k2;
        }
      }

      @Override
      public double getFraction(int index) {
        switch (index) {
          case 0:
            return f0;
          case 1:
            return f1;
          default:
            return f2;
        }
      }

      @Override
      public double getUpperBound() {
        return upper;
      }
    };
    Assertions.assertEquals(1 / k0, m.getResidenceTime(0));
    Assertions.assertEquals(1 / k1, m.getResidenceTime(1));
    Assertions.assertEquals(1 / k2, m.getResidenceTime(2));
    Assertions.assertEquals(0, m.sf(upper), () -> "sf: " + upper);
    Assertions.assertEquals(0, m.sf(upper * 2), () -> "sf: " + upper * 2);
    final DoubleDoubleBiPredicate test = Predicates.doublesAreRelativelyClose(1e-10);
    final double cdfU0 = d0.cumulativeProbability(upper);
    final double cdfU1 = d1.cumulativeProbability(upper);
    final double cdfU2 = d2.cumulativeProbability(upper);
    // This iteration (x < 5) is always below upper
    DoubleStream.iterate(0, x -> x + 0.5).limit(10).forEach(t -> {
      final double e = f0 * (1 - d0.cumulativeProbability(t) / cdfU0)
          + f1 * (1 - d1.cumulativeProbability(t) / cdfU1)
          + f2 * (1 - d2.cumulativeProbability(t) / cdfU2);
      TestAssertions.assertTest(e, m.sf(t), test, () -> "sf: " + t);
    });
  }

  @ParameterizedTest
  @CsvSource({"0.1, Infinity", "0.3, Infinity", "0.1, 7", "0.1, 4",})
  void testModel1(double k0, double upper) {
    final ExponentialDistribution d = ExponentialDistribution.of(1 / k0);
    final Model m = ResidenceTimeFitting.createModel(k0, upper);
    Assertions.assertEquals(1 / k0, m.getResidenceTime(0));
    Assertions.assertEquals(1, m.getSize());
    Assertions.assertEquals(k0, m.getRate(0));
    Assertions.assertEquals(1, m.getFraction(0));
    Assertions.assertEquals(upper, m.getUpperBound());
    final DoubleDoubleBiPredicate test = Predicates.doublesAreRelativelyClose(1e-10);
    if (Double.isInfinite(upper)) {
      DoubleStream.iterate(0, x -> x + 0.5).limit(10).forEach(t -> {
        final double e = d.survivalProbability(t);
        TestAssertions.assertTest(e, m.sf(t), test, () -> "sf: " + t);
      });
    } else {
      final double cdfU = d.cumulativeProbability(upper);
      TestAssertions.assertTest(0, m.sf(upper), test, () -> "sf: " + upper);
      TestAssertions.assertTest(0, m.sf(2 * upper), test, () -> "sf: " + 2 * upper);
      final int n = 10;
      IntStream.rangeClosed(0, n).mapToDouble(i -> upper * i / n).forEach(t -> {
        final double e = 1 - d.cumulativeProbability(t) / cdfU;
        TestAssertions.assertTest(e, m.sf(t), test, () -> "sf: " + t);
      });
    }
  }

  @ParameterizedTest
  @CsvSource({"0.1, 3, 0.5, Infinity", "2, 4, 0.25, 0.35", "2, 4, 0.25, 0.5",})
  void testModel2(double k0, double k1, double f0, double upper) {
    final ExponentialDistribution d0 = ExponentialDistribution.of(1 / k0);
    final ExponentialDistribution d1 = ExponentialDistribution.of(1 / k1);
    final double f1 = 1 - f0;
    final Model m = ResidenceTimeFitting.createModel(k0, k1, f0, upper);
    Assertions.assertEquals(1 / k0, m.getResidenceTime(0));
    Assertions.assertEquals(1 / k1, m.getResidenceTime(1));
    Assertions.assertEquals(2, m.getSize());
    Assertions.assertEquals(k0, m.getRate(0));
    Assertions.assertEquals(k1, m.getRate(1));
    Assertions.assertEquals(f0, m.getFraction(0));
    Assertions.assertEquals(f1, m.getFraction(1));
    Assertions.assertEquals(upper, m.getUpperBound());
    final DoubleDoubleBiPredicate test = Predicates.doublesAreRelativelyClose(1e-10);
    if (Double.isInfinite(upper)) {
      DoubleStream.iterate(0, x -> x + 0.5).limit(10).forEach(t -> {
        final double e = f0 * d0.survivalProbability(t) + f1 * d1.survivalProbability(t);
        TestAssertions.assertTest(e, m.sf(t), test, () -> "sf: " + t);
      });
    } else {
      final double cdfU0 = d0.cumulativeProbability(upper);
      final double cdfU1 = d1.cumulativeProbability(upper);
      TestAssertions.assertTest(0, m.sf(upper), test, () -> "sf: " + upper);
      TestAssertions.assertTest(0, m.sf(2 * upper), test, () -> "sf: " + 2 * upper);
      final int n = 10;
      IntStream.rangeClosed(0, n).mapToDouble(i -> upper * i / n).forEach(t -> {
        final double e = f0 * (1 - d0.cumulativeProbability(t) / cdfU0)
            + f1 * (1 - d1.cumulativeProbability(t) / cdfU1);
        TestAssertions.assertTest(e, m.sf(t), test, () -> "sf: " + t);
      });
    }
  }

  @ParameterizedTest
  @CsvSource({"0.1, 3, 5, 0.5, 0.25, Infinity", "2, 4, 10, 0.25, 0.125, 0.5",
      "2, 4, 10, 0.25, 0.125, 0.375",})
  void testModel3(double k0, double k1, double k2, double f0, double f1, double upper) {
    final ExponentialDistribution d0 = ExponentialDistribution.of(1 / k0);
    final ExponentialDistribution d1 = ExponentialDistribution.of(1 / k1);
    final ExponentialDistribution d2 = ExponentialDistribution.of(1 / k2);
    final double f2 = 1 - (f0 + f1);
    final Model m = ResidenceTimeFitting.createModel(k0, k1, k2, f0, f1, upper);
    Assertions.assertEquals(1 / k0, m.getResidenceTime(0));
    Assertions.assertEquals(1 / k1, m.getResidenceTime(1));
    Assertions.assertEquals(1 / k2, m.getResidenceTime(2));
    Assertions.assertEquals(3, m.getSize());
    Assertions.assertEquals(k0, m.getRate(0));
    Assertions.assertEquals(k1, m.getRate(1));
    Assertions.assertEquals(k2, m.getRate(2));
    Assertions.assertEquals(f0, m.getFraction(0));
    Assertions.assertEquals(f1, m.getFraction(1));
    Assertions.assertEquals(f2, m.getFraction(2));
    Assertions.assertEquals(upper, m.getUpperBound());
    final DoubleDoubleBiPredicate test = Predicates.doublesAreRelativelyClose(1e-10);
    if (Double.isInfinite(upper)) {
      DoubleStream.iterate(0, x -> x + 0.5).limit(10).forEach(t -> {
        final double e = f0 * d0.survivalProbability(t) + f1 * d1.survivalProbability(t)
            + f2 * d2.survivalProbability(t);
        TestAssertions.assertTest(e, m.sf(t), test, () -> "sf: " + t);
      });
    } else {
      final double cdfU0 = d0.cumulativeProbability(upper);
      final double cdfU1 = d1.cumulativeProbability(upper);
      final double cdfU2 = d2.cumulativeProbability(upper);
      TestAssertions.assertTest(0, m.sf(upper), test, () -> "sf: " + upper);
      TestAssertions.assertTest(0, m.sf(2 * upper), test, () -> "sf: " + 2 * upper);
      final int n = 10;
      IntStream.rangeClosed(0, n).mapToDouble(i -> upper * i / n).forEach(t -> {
        final double e = f0 * (1 - d0.cumulativeProbability(t) / cdfU0)
            + f1 * (1 - d1.cumulativeProbability(t) / cdfU1)
            + f2 * (1 - d2.cumulativeProbability(t) / cdfU2);
        TestAssertions.assertTest(e, m.sf(t), test, () -> "sf: " + t);
      });
    }
  }

  @ParameterizedTest
  @CsvSource({"0.1, 50, 1", "1.5, 50, 0.125", "1.5, 100, 0.25"})
  void testFunction1(double k0, int size, double resolution) {
    // Create data
    final UniformRandomProvider rng = RngFactory.create(0xbeefc0de);
    final double[] x =
        new DataSample(new double[] {1 / k0}, new double[] {1}).getSample(rng, size)[0];
    final int[] n = DataSample.toBins(x, resolution);

    final Function1 f = new Function1(n, resolution);

    final DoubleDoubleBiPredicate testExp = Predicates.doublesAreRelativelyClose(0.05);
    final DoubleDoubleBiPredicate testExp2 = Predicates.doublesAreRelativelyClose(1e-6);
    final DoubleDoubleBiPredicate testGeo = Predicates.doublesAreRelativelyClose(1e-10);
    for (double s = 1.0 / 4; s <= 4; s *= 2) {
      final double k = k0 * s;
      final double ll = f.ll(k);
      final ExponentialDistribution de = ExponentialDistribution.of(1 / k);
      double e = 0;
      for (int i = 0; i < x.length; i++) {
        e += de.logDensity(x[i]);
      }
      TestAssertions.assertTest(e, ll, testExp, () -> "rate = " + k);
      final GeometricDistribution dg = GeometricDistribution.of(-Math.expm1(-k * resolution));
      e = 0;
      double e2 = 0;
      final double logr = Math.log(resolution);
      for (int i = 0; i < n.length; i++) {
        // Scaled probability maps each continuous interval to a discrete PMF value: p / resolution
        // log(p) - log(resolution)
        e += n[i] * Math.log(de.probability(i * resolution, (i + 1) * resolution) / resolution);
        e2 += n[i] * (dg.logProbability(i) - logr);
      }
      TestAssertions.assertTest(e, ll, testExp2, () -> "rate = " + k);
      TestAssertions.assertTest(e2, ll, testGeo, () -> "rate = " + k);
    }
  }

  @ParameterizedTest
  @CsvSource({"0.1, 50, 1, 7", "1.5, 50, 0.125, 0.5", "1.5, 100, 0.25, 0.75", "1.5, 100, 0.25, 1"})
  void testFunction1T(double k0, int size, double resolution, double upper) {
    // Create data
    final UniformRandomProvider rng = RngFactory.create(0xbeefc0de);
    final double[] x =
        new DataSample(new double[] {1 / k0}, new double[] {1}).getSample(rng, size)[0];
    final int[] n = DataSample.toBins(x, resolution, upper);

    final Function1T f = new Function1T(n, resolution, upper);

    final DoubleDoubleBiPredicate testExp = Predicates.doublesAreRelativelyClose(1e-6);
    final DoubleDoubleBiPredicate testGeo = Predicates.doublesAreRelativelyClose(1e-10);
    for (double s = 1.0 / 4; s <= 4; s *= 2) {
      final double k = k0 * s;
      final double ll = f.ll(k);
      final ExponentialDistribution de = ExponentialDistribution.of(1 / k);
      final double norme = Math.log(de.cumulativeProbability(upper));
      final GeometricDistribution dg = GeometricDistribution.of(-Math.expm1(-k * resolution));
      final double normg = Math.log(dg.cumulativeProbability((int) (upper / resolution) - 1));
      double e = 0;
      double e2 = 0;
      final double logr = Math.log(resolution);
      for (int i = 0; i < n.length; i++) {
        // Scaled probability maps each continuous interval to a discrete PMF value: p / resolution
        // log(p) - log(resolution)
        e += n[i]
            * (Math.log(de.probability(i * resolution, (i + 1) * resolution) / resolution) - norme);
        e2 += n[i] * (dg.logProbability(i) - logr - normg);
      }
      TestAssertions.assertTest(e, ll, testExp, () -> "rate = " + k);
      TestAssertions.assertTest(e2, ll, testGeo, () -> "rate = " + k);
    }
  }

  @ParameterizedTest
  @CsvSource({"0.1, 3, 0.5, 200, 0.125", "2, 4, 0.25, 200, 0.0625"})
  void testFunction2(double k0, double k1, double f0, int size, double resolution) {
    // Create data
    final double f1 = 1 - f0;
    final UniformRandomProvider rng = RngFactory.create(0xdeadbeef);
    final double[] x = new DataSample(new double[] {1 / k0, 1 / k1}, new double[] {f0, f1})
        .getSample(rng, size)[0];
    final int[] n = DataSample.toBins(x, resolution);

    final Function2 f = new Function2(n, resolution);

    // The LL values are not close when the expected approaches close to zero
    // so we lower the tolerance. At the optimum and far from the optimum the LL
    // is far from zero and here the values are close.
    final DoubleDoubleBiPredicate testExp = Predicates.doublesAreRelativelyClose(0.05);
    final DoubleDoubleBiPredicate testExp1 = Predicates.doublesAreRelativelyClose(0.3);
    final DoubleDoubleBiPredicate testExp2 = Predicates.doublesAreRelativelyClose(1e-6);
    final DoubleDoubleBiPredicate testGeo = Predicates.doublesAreRelativelyClose(1e-10);
    for (double s = 1.0 / 4; s <= 4; s *= 2) {
      final double ka = k0 * s;
      final double kb = k1 * s;
      final double ll = f.ll(ka, kb, f0);
      final ExponentialDistribution de0 = ExponentialDistribution.of(1 / ka);
      final ExponentialDistribution de1 = ExponentialDistribution.of(1 / kb);
      double e = 0;
      for (int i = 0; i < x.length; i++) {
        e += Math.log(f0 * de0.density(x[i]) + f1 * de1.density(x[i]));
      }
      TestAssertions.assertTest(e, ll, Math.abs(e) < 5 ? testExp1 : testExp);
      final GeometricDistribution dg0 = GeometricDistribution.of(-Math.expm1(-ka * resolution));
      final GeometricDistribution dg1 = GeometricDistribution.of(-Math.expm1(-kb * resolution));
      e = 0;
      double e2 = 0;
      for (int i = 0; i < n.length; i++) {
        // Scaled probability maps each continuous interval to a discrete PMF value: p / resolution
        e += n[i] * Math.log(f0 * de0.probability(i * resolution, (i + 1) * resolution) / resolution
            + f1 * de1.probability(i * resolution, (i + 1) * resolution) / resolution);
        e2 += n[i] * Math.log((f0 * dg0.probability(i) + f1 * dg1.probability(i)) / resolution);
      }
      TestAssertions.assertTest(e, ll, testExp2, () -> "rate = " + ka + " " + kb);
      TestAssertions.assertTest(e2, ll, testGeo, () -> "rate = " + ka + " " + kb);
    }
  }

  @ParameterizedTest
  @CsvSource({"0.1, 3, 0.5, 200, 0.125, 7", "2, 4, 0.25, 200, 0.0625, 0.75",
      "2, 4, 0.25, 200, 0.0625, 0.5",})
  void testFunction2T(double k0, double k1, double f0, int size, double resolution, double upper) {
    // Create data
    final double f1 = 1 - f0;
    final UniformRandomProvider rng = RngFactory.create(0xdeadbeef);
    final double[] x = new DataSample(new double[] {1 / k0, 1 / k1}, new double[] {f0, f1})
        .getSample(rng, size)[0];
    final int[] n = DataSample.toBins(x, resolution, upper);

    final Function2T f = new Function2T(n, resolution, upper);

    // The LL values are not close when the expected approaches close to zero
    // so we lower the tolerance. At the optimum and far from the optimum the LL
    // is far from zero and here the values are close.
    final DoubleDoubleBiPredicate testExp = Predicates.doublesAreRelativelyClose(1e-6);
    final DoubleDoubleBiPredicate testGeo = Predicates.doublesAreRelativelyClose(1e-10);
    for (double s = 1.0 / 4; s <= 4; s *= 2) {
      final double ka = k0 * s;
      final double kb = k1 * s;
      final double ll = f.ll(ka, kb, f0);
      final ExponentialDistribution de0 = ExponentialDistribution.of(1 / ka);
      final ExponentialDistribution de1 = ExponentialDistribution.of(1 / kb);
      final double norme0 = de0.cumulativeProbability(upper);
      final double norme1 = de1.cumulativeProbability(upper);
      final GeometricDistribution dg0 = GeometricDistribution.of(-Math.expm1(-ka * resolution));
      final GeometricDistribution dg1 = GeometricDistribution.of(-Math.expm1(-kb * resolution));
      final double normg0 = dg0.cumulativeProbability((int) (upper / resolution) - 1);
      final double normg1 = dg1.cumulativeProbability((int) (upper / resolution) - 1);
      double e = 0;
      double e2 = 0;
      for (int i = 0; i < n.length; i++) {
        // Scaled probability maps each continuous interval to a discrete PMF value: p / resolution
        e += n[i] * Math
            .log(f0 * de0.probability(i * resolution, (i + 1) * resolution) / resolution / norme0
                + f1 * de1.probability(i * resolution, (i + 1) * resolution) / resolution / norme1);
        e2 += n[i] * Math.log(
            (f0 * dg0.probability(i) / normg0 + f1 * dg1.probability(i) / normg1) / resolution);
      }
      TestAssertions.assertTest(e, ll, testExp, () -> "rate = " + ka + " " + kb);
      TestAssertions.assertTest(e2, ll, testGeo, () -> "rate = " + ka + " " + kb);
    }
  }

  @ParameterizedTest
  @CsvSource({"0.1, 3, 5, 0.5, 0.25, 200, 0.125", "2, 4, 10, 0.25, 0.5, 200, 0.0625"})
  void testFunction3(double k0, double k1, double k2, double f0, double f1, int size,
      double resolution) {
    // Create data
    final double f2 = 1 - (f0 + f1);
    final UniformRandomProvider rng = RngFactory.create(0xdeadbeef);
    final double[] x =
        new DataSample(new double[] {1 / k0, 1 / k1, 1 / k2}, new double[] {f0, f1, f2})
            .getSample(rng, size)[0];
    final int[] n = DataSample.toBins(x, resolution);

    final Function3 f = new Function3(n, resolution);

    // The LL values are not a close match for direct computation using the log density.
    final DoubleDoubleBiPredicate testExp = Predicates.doublesAreRelativelyClose(0.3);
    final DoubleDoubleBiPredicate testExp2 = Predicates.doublesAreRelativelyClose(1e-6);
    final DoubleDoubleBiPredicate testGeo = Predicates.doublesAreRelativelyClose(1e-10);
    for (double s = 1.0 / 4; s <= 4; s *= 2) {
      final double ka = k0 * s;
      final double kb = k1 * s;
      final double kc = k2 * s;
      final double ll = f.ll(ka, kb, kc, f0, f1, f2);
      final ExponentialDistribution de0 = ExponentialDistribution.of(1 / ka);
      final ExponentialDistribution de1 = ExponentialDistribution.of(1 / kb);
      final ExponentialDistribution de2 = ExponentialDistribution.of(1 / kc);
      double e = 0;
      for (int i = 0; i < x.length; i++) {
        e += Math.log(f0 * de0.density(x[i]) + f1 * de1.density(x[i]) + f2 * de2.density(x[i]));
      }
      TestAssertions.assertTest(e, ll, testExp);
      final GeometricDistribution dg0 = GeometricDistribution.of(-Math.expm1(-ka * resolution));
      final GeometricDistribution dg1 = GeometricDistribution.of(-Math.expm1(-kb * resolution));
      final GeometricDistribution dg2 = GeometricDistribution.of(-Math.expm1(-kc * resolution));
      e = 0;
      double e2 = 0;
      for (int i = 0; i < n.length; i++) {
        // Scaled probability maps each continuous interval to a discrete PMF value: p / resolution
        e += n[i] * Math.log(f0 * de0.probability(i * resolution, (i + 1) * resolution) / resolution
            + f1 * de1.probability(i * resolution, (i + 1) * resolution) / resolution
            + f2 * de2.probability(i * resolution, (i + 1) * resolution) / resolution);
        e2 += n[i]
            * Math.log((f0 * dg0.probability(i) + f1 * dg1.probability(i) + f2 * dg2.probability(i))
                / resolution);
      }
      TestAssertions.assertTest(e, ll, testExp2, () -> "rate = " + ka + " " + kb + " " + kc);
      TestAssertions.assertTest(e2, ll, testGeo, () -> "rate = " + ka + " " + kb + " " + kc);
    }
  }

  @ParameterizedTest
  @CsvSource({"0.1, 3, 5, 0.5, 0.25, 200, 0.125, 4", "2, 4, 10, 0.25, 0.5, 200, 0.0625, 0.5625",
      "2, 4, 10, 0.25, 0.5, 200, 0.0625, 0.5",})
  void testFunction3T(double k0, double k1, double k2, double f0, double f1, int size,
      double resolution, double upper) {
    // Create data
    final double f2 = 1 - (f0 + f1);
    final UniformRandomProvider rng = RngFactory.create(0xdeadbeef);
    final double[] x =
        new DataSample(new double[] {1 / k0, 1 / k1, 1 / k2}, new double[] {f0, f1, f2})
            .getSample(rng, size)[0];
    final int[] n = DataSample.toBins(x, resolution, upper);

    final Function3T f = new Function3T(n, resolution, upper);

    // The LL values are not a close match for direct computation using the log density.
    final DoubleDoubleBiPredicate testExp = Predicates.doublesAreRelativelyClose(1e-6);
    final DoubleDoubleBiPredicate testGeo = Predicates.doublesAreRelativelyClose(1e-10);
    for (double s = 1.0 / 4; s <= 4; s *= 2) {
      final double ka = k0 * s;
      final double kb = k1 * s;
      final double kc = k2 * s;
      final double ll = f.ll(ka, kb, kc, f0, f1, f2);
      final ExponentialDistribution de0 = ExponentialDistribution.of(1 / ka);
      final ExponentialDistribution de1 = ExponentialDistribution.of(1 / kb);
      final ExponentialDistribution de2 = ExponentialDistribution.of(1 / kc);
      final double norme0 = de0.cumulativeProbability(upper);
      final double norme1 = de1.cumulativeProbability(upper);
      final double norme2 = de2.cumulativeProbability(upper);
      final GeometricDistribution dg0 = GeometricDistribution.of(-Math.expm1(-ka * resolution));
      final GeometricDistribution dg1 = GeometricDistribution.of(-Math.expm1(-kb * resolution));
      final GeometricDistribution dg2 = GeometricDistribution.of(-Math.expm1(-kc * resolution));
      final double normg0 = dg0.cumulativeProbability((int) (upper / resolution) - 1);
      final double normg1 = dg1.cumulativeProbability((int) (upper / resolution) - 1);
      final double normg2 = dg2.cumulativeProbability((int) (upper / resolution) - 1);
      double e = 0;
      double e2 = 0;
      for (int i = 0; i < n.length; i++) {
        // Scaled probability maps each continuous interval to a discrete PMF value: p / resolution
        e += n[i] * Math
            .log(f0 * de0.probability(i * resolution, (i + 1) * resolution) / resolution / norme0
                + f1 * de1.probability(i * resolution, (i + 1) * resolution) / resolution / norme1
                + f2 * de2.probability(i * resolution, (i + 1) * resolution) / resolution / norme2);
        e2 += n[i] * Math.log((f0 * dg0.probability(i) / normg0 + f1 * dg1.probability(i) / normg1
            + f2 * dg2.probability(i) / normg2) / resolution);
      }
      TestAssertions.assertTest(e, ll, testExp, () -> "rate = " + ka + " " + kb + " " + kc);
      TestAssertions.assertTest(e2, ll, testGeo, () -> "rate = " + ka + " " + kb + " " + kc);
    }
  }

  // Based on the paper: Kim et al (2021)
  // Single-molecule imaging of chromatin remodelers reveals role of ATPase in promoting fast
  // kinetics of target search and dissociation from chromatin. eLife 10:e69387.
  // https://doi.org/10.7554/eLife.69387
  //
  // This paper fitted a model for stable- and transient- binding events using the CDF of
  // residence times:
  // 1 - CDF(t) = f_sb e^{-k_sb t} + f_tb e^{-k_tb t} ; f_sb + f_tb == 1
  // rho_sb = 1 / (k_sb + k_sb,H2B)
  // Example results show:
  // Sample, f_sb (%), rho_sb (s)
  // RSC, 27 +/- 2, 5.0 +/- 0.7
  // SWI/SNF, 24 +/- 6, 4.4 +/- 1.2
  // ISW2, 13 +/- 3, 4.9 +/- 2.2
  // CHD1, 15 +/- 3, 7.2 +/- 3.3
  // (Errors from bootstrap re-sampling 100 repeats)
  // No data is presented for the original k_sb and k_tb (rates)
  // Mean of exponential = 1/k_sb and 1/k_tb
  // => 1/k_sb ~ 5 s ; 1/k_tb is faster

  final DoubleDoubleBiPredicate deltaR = Predicates.doublesAreClose(0.1, 0);
  final DoubleDoubleBiPredicate deltaF = Predicates.doublesAreClose(0.2, 0);

  // Used for testing dual populations:
  // 10-fold, 50-fold, 500-fold difference between residence times
  // Note:
  // The exponential is mapped to a discrete time using floor(x). The 99th percentile
  // of an exponential distribution is at ~ 4.6. If the time resolution is much more
  // than 5 times the mean residence time then all observations will be in a single
  // time point and the fit value is not accurate with small n
  static final double[] RHO = new double[] {5, 0.5, 0.01};

  // @formatter:off
  @SeededTest
  void canFitSinglePopulation0_1(RandomSeed seed) { fitSinglePopulation(seed, 0.1); }
  @SeededTest
  void canFitSinglePopulation1(RandomSeed seed) { fitSinglePopulation(seed, 1); }
  @SeededTest
  void canFitSinglePopulation5(RandomSeed seed) { fitSinglePopulation(seed, 5); }
  // exp sf(7.5) * 16000 = 8.85
  @SeededTest
  void canFitSinglePopulation7_5(RandomSeed seed) { fitSinglePopulation(seed, 7.5); }
  // @formatter:on

  // @SeededTest
  void canFitSinglePopulation(RandomSeed seed) {
    // Debug limit
    // Reached at 11.39: exp sf(11.39) * 16000 = 0.179
    // Not enough non-zero counts in the next bin (i.e. all counts are in n[0])
    double s = 1;
    for (;;) {
      s *= 1.5;
      System.out.printf("scale %s%n", s);
      fitSinglePopulation(seed, s);
    }
  }

  private void fitSinglePopulation(RandomSeed seed, double scale) {
    final UniformRandomProvider rng = RngFactory.create(seed.get());
    // final UniformRandomProvider rng = RngFactory.create(0xdeadbeef);
    final String title = "Single  ";
    Throwable error = null;
    // Only need to fit 1 RHO as the resolution defines how hard it is to fit
    final double r = RHO[0];
    final double resolution = r * scale;
    for (int samples = 500, k = 0; k < 6; samples *= 2, k++) {
      try {
        fit(rng, title, samples, 1, new double[] {r}, new double[] {1}, resolution);
        return;
      } catch (final IllegalArgumentException | AssertionError ex) {
        error = ex;
      }
    }
    if (error instanceof RuntimeException) {
      throw (RuntimeException) error;
    }
    if (error instanceof Error) {
      throw (Error) error;
    }
    Assertions.assertNull(error);
  }

  // @formatter:off
  @SeededTest
  void canFitDual0_1Population(RandomSeed seed) { fitDualPopulation(seed, 0.25, 0.1); }
  @SeededTest
  void canFitDual0_2Population(RandomSeed seed) { fitDualPopulation(seed, 0.25, 0.2); }
  @SeededTest
  void canFitDual0_3Population(RandomSeed seed) { fitDualPopulation(seed, 0.25, 0.3); }
  @SeededTest
  void canFitDual0_4Population(RandomSeed seed) { fitDualPopulation(seed, 0.25, 0.4); }
  @SeededTest
  void canFitDual0_5Population(RandomSeed seed) { fitDualPopulation(seed, 0.25, 0.25); }
  @SeededTest
  void canFitDual0_6Population(RandomSeed seed) { fitDualPopulation(seed, 0.25, 0.6); }
  @SeededTest
  void canFitDual0_7Population(RandomSeed seed) { fitDualPopulation(seed, 0.25, 0.7); }
  @SeededTest
  void canFitDual0_8Population(RandomSeed seed) { fitDualPopulation(seed, 0.25, 0.8); }
  @SeededTest
  void canFitDual0_9Population(RandomSeed seed) { fitDualPopulation(seed, 0.25, 0.9); }
  // @formatter:on

  @Test
  void canFitDual0_2PopulationFixedSeed() {
    // Ensure the dual-population fit is performed at least once
    fitDualPopulation(null, 0.25, 0.2, 1236478126821268L);
  }

  // @SeededTest
  void canFitDualPopulation(RandomSeed seed) {
    // Debug limit
    // Reached at 1.90 or sometimes 2.85: exp sf(2.85) ~ 0.0578
    // exp sf(2.85) * 20000 * 0.8 = 925
    // Too few counts in the second bin to allow the smaller residence time to be detected ???
    double s = 0.25;
    for (;;) {
      s *= 1.5;
      System.out.printf("scale %s%n", s);
      fitDualPopulation(seed, s, 0.2, 1);
    }
  }

  private void fitDualPopulation(RandomSeed seed, double scale, double fraction) {
    fitDualPopulation(seed, scale, fraction, 0);
  }

  private void fitDualPopulation(RandomSeed seed, double scale, double fraction, long longSeed) {
    Assumptions.assumeTrue(longSeed != 0 || TestSettings.allow(TestComplexity.MAXIMUM));
    final UniformRandomProvider rng =
        seed != null ? RngFactory.create(seed.get()) : RngFactory.create(longSeed);

    final String title = String.format("Dual=%.2f", fraction);
    AssertionError error = null;
    for (int i = 0; i < RHO.length; i++) {
      NEXT_D: for (int j = i + 1; j < RHO.length; j++) {
        // Set resolution using the smallest residence time
        final double resolution = scale * Math.min(RHO[i], RHO[j]);
        for (int samples = 5000, k = 0; k < 3; samples *= 2, k++) {
          try {
            fit(rng, title, samples, 2, new double[] {RHO[i], RHO[j]},
                new double[] {fraction, 1 - fraction}, resolution);
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

  // @formatter:off
  @SeededTest
  void canFitTriple0_2_0_2Population(RandomSeed s) { fitTriplePopulation(s, 0.125, 0.2, 0.2); }
  @SeededTest
  void canFitTriple0_2_0_3Population(RandomSeed s) { fitTriplePopulation(s, 0.125, 0.2, 0.3); }
  @SeededTest
  void canFitTriple0_2_0_4Population(RandomSeed s) { fitTriplePopulation(s, 0.125, 0.2, 0.4); }
  @SeededTest
  void canFitTriple0_2_0_5Population(RandomSeed s) { fitTriplePopulation(s, 0.125, 0.2, 0.5); }
  // @formatter:on

  @Test
  void canFitTriple0_2PopulationFixedSeed() {
    // Ensure the triple-population fit is performed at least once
    fitTriplePopulation(null, 0.25, 0.2, 0.3, 1236478126821268L);
  }

  // @SeededTest
  void canFitTriplePopulation(RandomSeed seed) {
    // Debug limit
    // Triple population fitting is flakey and fails more than dual population
    double s = 0.125;
    for (;;) {
      s *= 1.5;
      System.out.printf("scale %s%n", s);
      fitTriplePopulation(seed, s, 0.2, 0.2, 1);
    }
  }

  private void fitTriplePopulation(RandomSeed seed, double scale, double f0, double f1) {
    fitTriplePopulation(seed, scale, f0, f1, 0);
  }

  private void fitTriplePopulation(RandomSeed seed, double scale, double f0, double f1,
      long longSeed) {
    Assumptions.assumeTrue(longSeed != 0 || TestSettings.allow(TestComplexity.MAXIMUM));
    final UniformRandomProvider rng =
        seed != null ? RngFactory.create(seed.get()) : RngFactory.create(longSeed);

    final double f2 = 1 - (f0 + f1);
    final String title = String.format("Triple=%.2f,%.2f", f0, f1);
    AssertionError error = null;
    for (int i = 0; i < RHO.length; i++) {
      for (int j = i + 1; j < RHO.length; j++) {
        NEXT_D: for (int k = j + 1; k < RHO.length; k++) {
          // Set resolution using the smallest residence time
          final double resolution = scale * MathUtils.min(RHO[i], RHO[j], RHO[k]);
          for (int samples = 5000, l = 0; l < 3; samples *= 2, l++) {
            try {
              fit(rng, title, samples, 3, new double[] {RHO[i], RHO[j], RHO[k]},
                  new double[] {f0, f1, f2}, resolution);
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
  }

  @Test
  void canFitTruncatedSinglePopulation() {
    final UniformRandomProvider rng = RngFactory.create(1236478126821268L);
    final double r = RHO[0];
    final double resolution = r * 0.125;
    fit(rng, "Single ", 5000, 1, new double[] {r}, new double[] {1}, resolution, r);
    fit(rng, "Single ", 5000, 1, new double[] {r}, new double[] {1}, resolution, r * 2);
  }

  @Test
  void canFitTruncatedDualPopulation() {
    final UniformRandomProvider rng = RngFactory.create(12364781926821268L);
    final double r0 = RHO[0];
    final double r1 = RHO[1];
    // Set resolution using the smallest residence time
    final double resolution = Math.min(r0, r1) * 0.125;
    // Truncate using the largest residence time.
    // This cannot remove too much of the the largest distribution.
    final double upper = Math.max(r0, r1);
    fit(rng, "Dual ", 5000, 2, new double[] {r0, r1}, new double[] {0.75, 0.25}, resolution,
        upper * 2);
    fit(rng, "Dual ", 10000, 2, new double[] {r0, r1}, new double[] {0.75, 0.25}, resolution,
        upper * 1.5);
  }

  @Test
  void canFitTruncatedTriplePopulation() {
    final UniformRandomProvider rng = RngFactory.create(1236478126821268L);
    final double r0 = RHO[0];
    final double r1 = RHO[1];
    final double r2 = RHO[2];
    // Set resolution using the smallest residence time
    final double resolution = MathUtils.min(r0, r1, r2) * 0.125;
    // Truncate using the largest residence time
    // This cannot remove too much of the the largest distribution.
    final double upper = MathUtils.max(r0, r1, r2);
    fit(rng, "Triple ", 5000, 3, new double[] {r0, r1, r2}, new double[] {0.5, 0.25, 0.25},
        resolution, upper * 2);
    fit(rng, "Triple ", 10000, 3, new double[] {r0, r1, r2}, new double[] {0.5, 0.25, 0.25},
        resolution, upper * 1.5);
  }

  private OutputStreamWriter out = null;

  // /**
  // * This is not actually a test but runs the fitting algorithm many times to collect benchmark
  // data
  // * to file.
  // */
  // @SeededTest
  // void canDoBenchmark(RandomSeed seed) {
  // // Skip this as it is slow
  // Assumptions.assumeTrue(false);
  // final UniformRandomProvider rng = RngFactory.create(seed.get());
  //
  // out = null;
  // try {
  // final FileOutputStream fos = new FileOutputStream("ResidenceTimeFittingTest.dat");
  // out = new OutputStreamWriter(fos, "UTF-8");
  //
  // // Run the fitting to produce benchmark data for a mixed population of 2
  // final int n = 2;
  // writeHeader(n);
  // for (int repeat = 10; repeat-- > 0;) {
  // resetData();
  // for (final boolean mle : new boolean[] {true, false}) {
  // for (int f = 1; f <= 9; f++) {
  // final double fraction = f / 10.0;
  // final String title = String.format("%s Dual=%.1f", (mle) ? "MLE" : "LSQ", fraction);
  // for (int samples = 500, k = 0; k < 6; samples *= 2, k++) {
  // for (int i = 0; i < RHO.length; i++) {
  // for (int j = i + 1; j < RHO.length; j++) {
  // try {
  // fit(rng, title, samples, 0, new double[] {RHO[i], RHO[j]},
  // new double[] {fraction, 1 - fraction}, mle);
  // } catch (final AssertionError ex) {
  // // Carry on with the benchmarking
  // }
  // // If the fit had the correct N then no need to repeat
  // if (fitN == n) {
  // continue;
  // }
  // try {
  // fit(rng, title + " Fixed", samples, n, new double[] {RHO[i], RHO[j]},
  // new double[] {fraction, 1 - fraction}, mle);
  // } catch (final AssertionError ex) {
  // // Carry on with the benchmarking
  // }
  // }
  // }
  // }
  // }
  // }
  // }
  // } catch (final Exception ex) {
  // throw new AssertionError("Failed to complete benchmark", ex);
  // } finally {
  // closeOutput();
  // }
  // }

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

  private void fit(UniformRandomProvider rng, String title, int samples, int n, double[] r,
      double[] fraction, double resolution) {
    fit(rng, title, samples, n, r, fraction, resolution, Double.POSITIVE_INFINITY);
  }

  private void fit(UniformRandomProvider rng, String title, int samples, int n, double[] r,
      double[] fraction, double resolution, double upper) {
    // Used for testing
    // @formatter:off
    // For easy mixed populations
    //if (f.length == 2 && Math.min(f[0],f[1])/(f[0]+f[1]) <= 0.2) return;
    // @formatter:on

    JumpDistanceAnalysis.sort(r, fraction);
    final double[] data = createData(rng, samples, r, fraction);
    final int[] h = DataSample.toBins(data, resolution, upper);

    final boolean truncated = upper < Double.POSITIVE_INFINITY;
    Assumptions.assumeTrue(h.length > 1,
        () -> String.format("Simulation did not generate enough timepoints: %d, %s, %s, %s, %s",
            samples, Arrays.toString(r), Arrays.toString(fraction), resolution, upper));

    final ResidenceTimeFitting rta = ResidenceTimeFitting.of(resolution, h, truncated);
    // rta.setLogger(logger);
    final Pair<FitResult, Model> fit = rta.fit(n);
    double[] fitR = new double[0];
    double[] fitF = fitR;
    if (fit != null) {
      final Model model = fit.getRight();
      Assertions.assertEquals(0, model.sf(upper), () -> "sf: " + upper);
      Assertions.assertEquals(0, model.sf(upper * 2), () -> "sf: " + upper * 2);
      Assertions.assertThrows(AssertionError.class, () -> model.getRate(n));
      Assertions.assertThrows(AssertionError.class, () -> model.getFraction(n));
      fitR = IntStream.range(0, n).mapToDouble(i -> 1 / model.getRate(i)).toArray();
      fitF = IntStream.range(0, n).mapToDouble(model::getFraction).toArray();
      // Check the fit result log-likelihood is consistent with the model
      final FitResult fr = fit.getLeft();
      Assertions.assertEquals(Arrays.stream(h).sum(), fr.getN(), "n");
      Assertions.assertEquals(2 * n - 1, fr.getP(), "parameters");
      final double ll = IntStream.range(0, h.length)
          .mapToDouble(i -> h[i]
              * Math.log((model.sf(i * resolution) - model.sf((i + 1) * resolution)) / resolution))
          .sum();
      Assertions.assertEquals(ll, fr.getLogLikelihood(), Math.abs(ll) * 0.05, "log-likelihood");
    }

    // Record results to file
    if (out != null) {
      writeResult(title, sample.mean, sample.fraction, samples, n, r, fraction, fitR, fitF);
    }

    AssertionError error = null;
    try {
      TestAssertions.assertArrayTest(r, fitR, deltaR, "Failed to fit residence time (r)");
      TestAssertions.assertArrayTest(fraction, fitF, deltaF, "Failed to fit f");
    } catch (final AssertionError ex) {
      error = ex;
    } finally {
      final double[] e1 = getPercentError(r, fitR);
      final double[] e2 = getPercentError(fraction, fitF);
      logger.log(TestLevel.TEST_INFO,
          FormatSupplier.getSupplier(
              "%s %s N=%d sample=%d, n=%d, res=%s : %s = %s [%s] : %s = %s [%s]",
              (error == null) ? "+++ Pass" : "--- Fail", title, r.length, samples, n, resolution,
              toString(r), toString(fitR), toString(e1), toString(fraction), toString(fitF),
              toString(e2)));
      if (error != null) {
        throw error;
      }
    }
  }

  private void writeHeader(int size) {
    final StringBuilder sb = new StringBuilder("title");
    sb.append('\t').append("repeat");
    for (int i = 0; i < size; i++) {
      sb.append('\t').append("R").append(i);
    }
    for (int i = 0; i < size; i++) {
      sb.append('\t').append("F").append(i);
    }
    sb.append('\t').append("samples");
    sb.append('\t').append("n");
    sb.append('\t').append("size");
    sb.append('\t').append("fsize");
    for (int i = 0; i < size; i++) {
      sb.append('\t').append("r").append(i);
    }
    for (int i = 0; i < size; i++) {
      sb.append('\t').append("fr").append(i);
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

  private void writeResult(String title, double[] actualR, double[] actualF, int samples, int n,
      double[] r, double[] fraction, double[] fitR, double[] fitFraction) {
    final int size = r.length;
    final int fitSize = fitR.length;
    // Pad results if they are too small
    if (fitSize < size) {
      fitR = Arrays.copyOf(fitR, size);
      fitFraction = Arrays.copyOf(fitFraction, size);
    }
    final double[] ed = getRelativeError(r, fitR);
    final double[] ef = getRelativeError(fraction, fitFraction);

    final StringBuilder sb = new StringBuilder(title);
    sb.append('\t').append(repeat);
    for (int i = 0; i < size; i++) {
      sb.append('\t').append(actualR[i]);
    }
    for (int i = 0; i < size; i++) {
      sb.append('\t').append(actualF[i]);
    }
    sb.append('\t').append(samples);
    sb.append('\t').append(n);
    sb.append('\t').append(size);
    sb.append('\t').append(fitSize);
    for (int i = 0; i < size; i++) {
      sb.append('\t').append(r[i]);
    }
    for (int i = 0; i < size; i++) {
      sb.append('\t').append(fitR[i]);
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
    return String.format("%.6f", value);
  }

  private static class DataSample {
    double[] mean;
    double[] fraction;
    double[] data;
    int[] sample;

    DataSample(double[] mean, double[] fraction) {
      this.mean = mean.clone();
      this.fraction = fraction.clone();
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

    double[][] getSample(UniformRandomProvider rng, int size) {
      if (size > getSize()) {
        final int extra = size - getSize();

        final double[] data = new double[extra];
        final int[] sample = new int[extra];

        // Residence times follow an exponential distribution with the specified mean
        if (fraction.length > 1) {
          // Pick the population using the fraction and apply the mean.
          final SharedStateDiscreteSampler select = AliasMethodDiscreteSampler.of(rng, fraction);
          final ZigguratSampler.Exponential exp = ZigguratSampler.Exponential.of(rng);
          for (int i = 0; i < data.length; i++) {
            sample[i] = select.sample();
            data[i] = exp.sample() * mean[sample[i]];
          }
        } else {
          // Fixed mean
          final ZigguratSampler.Exponential exp = ZigguratSampler.Exponential.of(rng, mean[0]);
          for (int i = 0; i < data.length; i++) {
            data[i] = exp.sample();
          }
        }
        add(data, sample);
      }

      // Build the sample data and return the mean residence time and fractions
      final double[] data = Arrays.copyOf(this.data, size);
      final double[] r = new double[this.mean.length];
      final double[] f = new double[r.length];
      for (int i = 0; i < size; i++) {
        final int j = sample[i];
        r[j] += data[i];
        f[j]++;
      }
      for (int i = 0; i < r.length; i++) {
        // Mean residence time
        r[i] /= f[i];
        // Normalise the fraction
        f[i] /= size;
      }

      // Sort by descending residence time
      JumpDistanceAnalysis.sort(r, f);

      return new double[][] {data, r, f};
    }

    /**
     * Convert continuous residence times to binned data.
     *
     * @param x the residence times
     * @param resolution the time resolution
     * @return the binned data
     */
    static int[] toBins(double[] x, double resolution) {
      // Get the max bin
      final int max = (int) (MathUtils.max(x) / resolution);
      final int[] h = new int[max + 1];
      // Assign to bins by rounding up to the time defined by the bin width.
      // (Avoids t=0 which errors in the ResidenceTimeFitting constructor.)
      Arrays.stream(x).mapToInt(t -> (int) (t / resolution)).forEach(i -> h[i]++);
      return h;
    }

    /**
     * Convert continuous residence times to binned data.
     *
     * @param x the residence times
     * @param resolution the time resolution
     * @param upper the upper limit of the doamin
     * @return the binned data
     */
    static int[] toBins(double[] x, double resolution, double upper) {
      // Get the max bin
      final int max = (int) (Arrays.stream(x).filter(d -> d < upper).max().orElse(0) / resolution);
      final int[] h = new int[max + 1];
      // Assign to bins by rounding up to the time defined by the bin width.
      // (Avoids t=0 which errors in the ResidenceTimeFitting constructor.)
      Arrays.stream(x).filter(d -> d < upper).mapToInt(t -> (int) (t / resolution))
          .forEach(i -> h[i]++);
      return h;
    }

    @Override
    public boolean equals(Object obj) {
      if (!(obj instanceof DataSample)) {
        return super.equals(obj);
      }
      final DataSample that = (DataSample) obj;
      if (that.mean.length != this.mean.length) {
        return false;
      }
      for (int i = mean.length; i-- > 0;) {
        if ((that.mean[i] != this.mean[i]) || (that.fraction[i] != this.fraction[i])) {
          return false;
        }
      }
      return true;
    }

    @Override
    public int hashCode() {
      int hash = 1;
      for (int i = mean.length; i-- > 0;) {
        hash = hash * 31 + Double.hashCode(mean[i]);
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
   * Create residence times.
   *
   * @param n Number of residence times
   * @param r Residence time
   * @param fraction Fraction of population (will be updated with the actual fractions, normalised
   *        to sum to 1)
   * @return The residence times
   */
  private double[] createData(UniformRandomProvider rng, int n, double[] r, double[] fraction) {
    // Cache the data so that if we run a second test with
    // the same r and f we use the same data
    sample = new DataSample(r, fraction);
    final int index = samples.indexOf(sample);
    if (index != -1) {
      sample = samples.get(index);
    } else {
      samples.add(sample);
    }

    final double[][] dataSample = sample.getSample(rng, n);
    final double[] data = dataSample[0];
    final double[] r2 = dataSample[1];
    final double[] f2 = dataSample[2];

    // System.out.printf("r: %s -> %s ; f: %s -> %s%n",
    // Arrays.toString(r), Arrays.toString(r2), Arrays.toString(fraction), Arrays.toString(f2));

    // Update with the real values
    System.arraycopy(r2, 0, r, 0, r.length);
    System.arraycopy(f2, 0, fraction, 0, fraction.length);

    return data;
  }
}
