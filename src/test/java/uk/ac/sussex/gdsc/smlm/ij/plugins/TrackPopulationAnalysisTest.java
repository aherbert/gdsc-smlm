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
import uk.ac.sussex.gdsc.smlm.ij.plugins.TrackPopulationAnalysis.OffsetPowerFunction;
import uk.ac.sussex.gdsc.smlm.ij.plugins.TrackPopulationAnalysis.OffsetPowerFunction1;
import uk.ac.sussex.gdsc.smlm.ij.plugins.TrackPopulationAnalysis.PowerFunction;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;

@SuppressWarnings({"javadoc"})
class TrackPopulationAnalysisTest {
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
