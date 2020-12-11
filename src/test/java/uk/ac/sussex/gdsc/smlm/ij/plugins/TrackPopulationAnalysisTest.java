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

package uk.ac.sussex.gdsc.smlm.ij.plugins;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.smlm.ij.plugins.TrackPopulationAnalysis.BrownianDiffusionFunction;
import uk.ac.sussex.gdsc.smlm.ij.plugins.TrackPopulationAnalysis.FbmDiffusionFunction;
import uk.ac.sussex.gdsc.smlm.ij.plugins.TrackPopulationAnalysis.OffsetPowerFunction;
import uk.ac.sussex.gdsc.smlm.ij.plugins.TrackPopulationAnalysis.OffsetPowerFunction1;
import uk.ac.sussex.gdsc.smlm.ij.plugins.TrackPopulationAnalysis.PowerFunction;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;

@SuppressWarnings({"javadoc"})
class TrackPopulationAnalysisTest {
  @Test
  void canComputeBrownianDiffusionFunction() {
    final int size = 10;
    final BrownianDiffusionFunction f = new BrownianDiffusionFunction(size);
    final double delta = 1e-6;
    final DoubleDoubleBiPredicate test = TestHelper.doublesAreClose(1e-5);
    for (final double d : new double[] {0.8, 0.9, 1, 1.1, 1.2}) {
      for (final double s : new double[] {2, 20}) {
        // Check the value and Jacobian
        final RealVector point = new ArrayRealVector(new double[] {d, s}, false);
        final Pair<RealVector, RealMatrix> p = f.value(point);
        final double[] value = p.getFirst().toArray();
        Assertions.assertEquals(size, value.length);
        for (int n = 1; n <= size; n++) {
          // MSD = 4Dt * (n - 1/3) + 4s^2
          Assertions.assertEquals(4 * (d * (n - 1.0 / 3) + s * s), value[n - 1], "value");
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

  @Test
  void canComputeFbmDiffusionFunction() {
    final int size = 10;
    final FbmDiffusionFunction f = new FbmDiffusionFunction(size);
    final double delta = 1e-6;
    final DoubleDoubleBiPredicate test = TestHelper.doublesAreClose(1e-5);
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
          for (int n = 1; n <= size; n++) {
            // MSD = [4Dt^a / (a+2)(a+1)] * [(n+1)^(a+2) + (n-1)^(a+2) - 2n^(a+2)]
            // - [8Dt^a / (a+2)(a+1)] + 4s^2
            // assume t=1.
            final double msd = (4 * d / (a2 * a1))
                * (Math.pow(n + 1, a2) + Math.pow(n - 1, a2) - 2 * Math.pow(n, a2))
                - 8 * d / (a2 * a1) + 4 * s * s;
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

  @Test
  void canComputeBrownianModelUsingFbmFunction() {
    final int size = 10;
    final BrownianDiffusionFunction f1 = new BrownianDiffusionFunction(size);
    final FbmDiffusionFunction f2 = new FbmDiffusionFunction(size);
    final DoubleDoubleBiPredicate test = TestHelper.doublesAreClose(1e-5);
    final double alpha = 1.0;
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

  @Test
  void canComputePowerFunction() {
    final int size = 10;
    final PowerFunction f = new PowerFunction(size);
    final double delta = 1e-6;
    final DoubleDoubleBiPredicate test = TestHelper.doublesAreClose(1e-5);
    for (final double alpha : new double[] {0.8, 0.9, 1, 1.1, 1.2}) {
      for (final double beta : new double[] {2, 20}) {
        // Check the value and Jacobian
        final RealVector point = new ArrayRealVector(new double[] {alpha, beta}, false);
        final Pair<RealVector, RealMatrix> p = f.value(point);
        final double[] value = p.getFirst().toArray();
        Assertions.assertEquals(size, value.length);
        for (int x = 0; x < size; x++) {
          Assertions.assertEquals(beta * Math.pow(x, alpha), value[x], "value");
        }
        // Columns of the Jacobian
        final double[] dfda1 = p.getSecond().getColumn(0);
        final double[] dfdb1 = p.getSecond().getColumn(1);
        point.setEntry(0, alpha - delta);
        RealVector v1 = f.value(point).getFirst();
        point.setEntry(0, alpha + delta);
        RealVector v2 = f.value(point).getFirst();
        final double[] dfda = v2.subtract(v1).mapDivide(2 * delta).toArray();
        point.setEntry(0, alpha);
        point.setEntry(1, beta - delta);
        v1 = f.value(point).getFirst();
        point.setEntry(1, beta + delta);
        v2 = f.value(point).getFirst();
        final double[] dfdb = v2.subtract(v1).mapDivide(2 * delta).toArray();
        // Element-by-element relative error
        TestAssertions.assertArrayTest(dfda, dfda1, test, "jacobian dfda");
        TestAssertions.assertArrayTest(dfdb, dfdb1, test, "jacobian dfdb");
      }
    }
  }

  @Test
  void canComputeOffsetPowerFunction() {
    final int size = 10;
    final OffsetPowerFunction f = new OffsetPowerFunction(size);
    final double delta = 1e-6;
    final DoubleDoubleBiPredicate test = TestHelper.doublesAreClose(1e-5);
    for (final double alpha : new double[] {0.8, 0.9, 1, 1.1, 1.2}) {
      for (final double beta : new double[] {2, 20}) {
        for (final double c : new double[] {0, 1.5}) {
          // Check the value and Jacobian
          final RealVector point = new ArrayRealVector(new double[] {alpha, beta, c}, false);
          final Pair<RealVector, RealMatrix> p = f.value(point);
          final double[] value = p.getFirst().toArray();
          Assertions.assertEquals(size, value.length);
          for (int x = 0; x < size; x++) {
            Assertions.assertEquals(c + beta * Math.pow(x, alpha), value[x], "value");
          }
          // Columns of the Jacobian
          final double[] dfda1 = p.getSecond().getColumn(0);
          final double[] dfdb1 = p.getSecond().getColumn(1);
          final double[] dfdc1 = p.getSecond().getColumn(2);
          point.setEntry(0, alpha - delta);
          RealVector v1 = f.value(point).getFirst();
          point.setEntry(0, alpha + delta);
          RealVector v2 = f.value(point).getFirst();
          final double[] dfda = v2.subtract(v1).mapDivide(2 * delta).toArray();
          point.setEntry(0, alpha);
          point.setEntry(1, beta - delta);
          v1 = f.value(point).getFirst();
          point.setEntry(1, beta + delta);
          v2 = f.value(point).getFirst();
          final double[] dfdb = v2.subtract(v1).mapDivide(2 * delta).toArray();
          point.setEntry(1, beta);
          point.setEntry(2, c - delta);
          v1 = f.value(point).getFirst();
          point.setEntry(2, c + delta);
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

  @Test
  void canComputeOffsetPowerFunction1() {
    final int size = 10;
    final OffsetPowerFunction1 f = new OffsetPowerFunction1(size);
    final double delta = 1e-6;
    final DoubleDoubleBiPredicate test = TestHelper.doublesAreClose(1e-5);
    for (final double alpha : new double[] {0.8, 0.9, 1, 1.1, 1.2}) {
      for (final double beta : new double[] {2, 20}) {
        for (final double c : new double[] {0, 1.5}) {
          // Check the value and Jacobian
          final RealVector point = new ArrayRealVector(new double[] {alpha, beta, c}, false);
          final Pair<RealVector, RealMatrix> p = f.value(point);
          final double[] value = p.getFirst().toArray();
          Assertions.assertEquals(size - 1, value.length);
          for (int x = 1; x < size; x++) {
            Assertions.assertEquals(c + beta * Math.pow(x, alpha), value[x - 1], "value");
          }
          // Columns of the Jacobian
          final double[] dfda1 = p.getSecond().getColumn(0);
          final double[] dfdb1 = p.getSecond().getColumn(1);
          final double[] dfdc1 = p.getSecond().getColumn(2);
          point.setEntry(0, alpha - delta);
          RealVector v1 = f.value(point).getFirst();
          point.setEntry(0, alpha + delta);
          RealVector v2 = f.value(point).getFirst();
          final double[] dfda = v2.subtract(v1).mapDivide(2 * delta).toArray();
          point.setEntry(0, alpha);
          point.setEntry(1, beta - delta);
          v1 = f.value(point).getFirst();
          point.setEntry(1, beta + delta);
          v2 = f.value(point).getFirst();
          final double[] dfdb = v2.subtract(v1).mapDivide(2 * delta).toArray();
          point.setEntry(1, beta);
          point.setEntry(2, c - delta);
          v1 = f.value(point).getFirst();
          point.setEntry(2, c + delta);
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
