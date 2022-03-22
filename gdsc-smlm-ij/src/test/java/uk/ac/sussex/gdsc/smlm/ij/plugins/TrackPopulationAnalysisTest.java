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

package uk.ac.sussex.gdsc.smlm.ij.plugins;

import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.Pair;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.rng.RandomGeneratorAdapter;
import uk.ac.sussex.gdsc.smlm.ij.plugins.TrackPopulationAnalysis.BrownianDiffusionFunction;
import uk.ac.sussex.gdsc.smlm.ij.plugins.TrackPopulationAnalysis.ExponentialDataFunction;
import uk.ac.sussex.gdsc.smlm.ij.plugins.TrackPopulationAnalysis.FbmDiffusionFunction;
import uk.ac.sussex.gdsc.test.api.Predicates;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.rng.RngFactory;

@SuppressWarnings({"javadoc"})
class TrackPopulationAnalysisTest {
  @Test
  void canComputeBrownianDiffusionFunction1() {
    final int size = 10;
    final double delta = 1e-6;
    final DoubleDoubleBiPredicate test = Predicates.doublesAreRelativelyClose(1e-5);
    for (final double t : new double[] {0.8, 1, 1.2}) {
      final MultivariateJacobianFunction f = new BrownianDiffusionFunction(size, t);
      for (final double d : new double[] {0.8, 0.9, 1, 1.1, 1.2}) {
        for (final double s : new double[] {2, 20}) {
          // Check the value and Jacobian
          final RealVector point = new ArrayRealVector(new double[] {d, s}, false);
          final Pair<RealVector, RealMatrix> p = f.value(point);
          final double[] value = p.getFirst().toArray();
          Assertions.assertEquals(size, value.length);
          for (int n = 1; n <= size; n++) {
            // MSD = 4Dt * (n - 1/3) + 4s^2
            final double msd = 4 * d * t * (n - 1.0 / 3) + 4 * s * s;
            TestAssertions.assertTest(msd, value[n - 1], test, "value");
          }
          // Columns of the Jacobian
          final double[] dfda1 = p.getSecond().getColumn(0);
          final double[] dfdb1 = p.getSecond().getColumn(1);
          point.setEntry(0, d - delta);
          RealVector v1 = f.value(point).getFirst();
          point.setEntry(0, d + delta);
          RealVector v2 = f.value(point).getFirst();
          final double[] dfda = v2.subtract(v1).mapDivide(2 * delta).toArray();
          point.setEntry(0, d);
          point.setEntry(1, s - delta);
          v1 = f.value(point).getFirst();
          point.setEntry(1, s + delta);
          v2 = f.value(point).getFirst();
          final double[] dfdb = v2.subtract(v1).mapDivide(2 * delta).toArray();
          // Element-by-element relative error
          TestAssertions.assertArrayTest(dfda, dfda1, test, "jacobian dfda");
          TestAssertions.assertArrayTest(dfdb, dfdb1, test, "jacobian dfdb");
        }
      }
    }
  }

  @Test
  void canComputeFbmDiffusionFunction() {
    final int size = 10;
    final double delta = 1e-6;
    final DoubleDoubleBiPredicate test = Predicates.doublesAreRelativelyClose(1e-5);
    for (final double t : new double[] {0.8, 1, 1.2}) {
      final MultivariateJacobianFunction f = new FbmDiffusionFunction(size, t);
      for (final double d : new double[] {0.8, 0.9, 1, 1.1, 1.2}) {
        for (final double s : new double[] {2, 20}) {
          for (final double alpha : new double[] {1e-2, 0.8, 0.9, 1, 1.1, 1.2, 2.0}) {
            // Check the value and Jacobian
            final RealVector point = new ArrayRealVector(new double[] {d, s, alpha}, false);
            final Pair<RealVector, RealMatrix> p = f.value(point);
            final double[] value = p.getFirst().toArray();
            Assertions.assertEquals(size, value.length);
            final double a1 = alpha + 1;
            final double a2 = alpha + 2;
            final double ta = Math.pow(t, alpha);
            for (int n = 1; n <= size; n++) {
              // MSD = [4Dt^a / (a+2)(a+1)] * [(n+1)^(a+2) + (n-1)^(a+2) - 2n^(a+2)]
              // - [8Dt^a / (a+2)(a+1)] + 4s^2
              // assume t=1.
              final double msd = (4 * d * ta / (a2 * a1))
                  * (Math.pow(n + 1, a2) + Math.pow(n - 1, a2) - 2 * Math.pow(n, a2))
                  - 8 * d * ta / (a2 * a1) + 4 * s * s;
              TestAssertions.assertTest(msd, value[n - 1], test, "value");
            }
            // Columns of the Jacobian
            final double[] dfda1 = p.getSecond().getColumn(0);
            final double[] dfdb1 = p.getSecond().getColumn(1);
            final double[] dfdc1 = p.getSecond().getColumn(2);
            point.setEntry(0, d - delta);
            RealVector v1 = f.value(point).getFirst();
            point.setEntry(0, d + delta);
            RealVector v2 = f.value(point).getFirst();
            final double[] dfda = v2.subtract(v1).mapDivide(2 * delta).toArray();
            point.setEntry(0, d);
            point.setEntry(1, s - delta);
            v1 = f.value(point).getFirst();
            point.setEntry(1, s + delta);
            v2 = f.value(point).getFirst();
            final double[] dfdb = v2.subtract(v1).mapDivide(2 * delta).toArray();
            point.setEntry(1, s);
            point.setEntry(2, alpha - delta);
            v1 = f.value(point).getFirst();
            point.setEntry(2, alpha + delta);
            v2 = f.value(point).getFirst();
            final double[] dfdc = v2.subtract(v1).mapDivide(2 * delta).toArray();
            // Element-by-element relative error
            TestAssertions.assertArrayTest(dfda, dfda1, test, "jacobian dfda");
            TestAssertions.assertArrayTest(dfdb, dfdb1, test, "jacobian dfdb");
            TestAssertions.assertArrayTest(dfdc, dfdc1, test, "jacobian dfdc");
          }
        }
      }
    }
  }

  @Test
  void canComputeBrownianModelUsingFbmFunction() {
    final int size = 10;
    final DoubleDoubleBiPredicate test = Predicates.doublesAreRelativelyClose(1e-5);
    final double alpha = 1.0;
    for (final double t : new double[] {0.8, 1, 1.2}) {
      final BrownianDiffusionFunction f1 = new BrownianDiffusionFunction(size, t);
      final FbmDiffusionFunction f2 = new FbmDiffusionFunction(size, t);
      for (final double d : new double[] {0.8, 0.9, 1, 1.1, 1.2}) {
        for (final double s : new double[] {2, 20}) {
          // Check the value and Jacobian
          final RealVector point1 = new ArrayRealVector(new double[] {d, s, alpha}, false);
          final Pair<RealVector, RealMatrix> p1 = f1.value(point1);
          final double[] value1 = p1.getFirst().toArray();
          final RealVector point2 = new ArrayRealVector(new double[] {d, s, alpha}, false);
          final Pair<RealVector, RealMatrix> p2 = f2.value(point2);
          final double[] value2 = p2.getFirst().toArray();
          TestAssertions.assertArrayTest(value1, value2, test, "value");
          final double[] dfda1 = p1.getSecond().getColumn(0);
          final double[] dfdb1 = p1.getSecond().getColumn(1);
          final double[] dfda2 = p2.getSecond().getColumn(0);
          final double[] dfdb2 = p2.getSecond().getColumn(1);
          TestAssertions.assertArrayTest(dfda1, dfda2, test, "jacobian dfda");
          TestAssertions.assertArrayTest(dfdb1, dfdb2, test, "jacobian dfdb");
        }
      }
    }
  }

  @Test
  void canComputeExponentialLogLikelihood() {
    final RandomGenerator rng = new RandomGeneratorAdapter(RngFactory.createWithFixedSeed());
    final double delta = 1e-6;
    for (final double mean : new double[] {2, 10}) {
      // Create exponential data
      ExponentialDistribution ed = new ExponentialDistribution(rng, mean);
      final double[] values = ed.sample(1000);
      // Discretise
      final double factor = mean / 10;
      for (int i = 0; i < values.length; i++) {
        values[i] = Math.round(values[i] / factor) * factor;
      }
      // Histogram these
      final double[][] h = MathUtils.cumulativeHistogram(values, false);
      final ExponentialDataFunction ef = ExponentialDataFunction.fromCdf(h);
      // Test with different means
      for (final double du : new double[] {-0.1, 0, 0.1}) {
        final double mu = mean + du;
        ed = new ExponentialDistribution(rng, mu);
        double ll = 0;
        for (final double x : values) {
          ll += ed.logDensity(x);
        }
        final double ll2 = ef.value(mu);
        Assertions.assertEquals(ll, ll2, 1e-10, "Log likelihood");
        // Check for gradient
        final double ll2a = ef.value(mu + delta);
        final double ll2b = ef.value(mu - delta);
        final double gradient = (ll2a - ll2b) / (2 * delta);
        Assertions.assertEquals(gradient, ef.gradient(mu), Math.abs(gradient) * 1e-4, "Gradient");
      }
    }
  }
}
