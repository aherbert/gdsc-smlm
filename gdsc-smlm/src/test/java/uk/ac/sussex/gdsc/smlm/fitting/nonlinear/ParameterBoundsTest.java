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

package uk.ac.sussex.gdsc.smlm.fitting.nonlinear;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.smlm.function.FakeGradientFunction;
import uk.ac.sussex.gdsc.test.utils.functions.IndexSupplier;

@SuppressWarnings({"javadoc"})
class ParameterBoundsTest {
  @Test
  void canUpperBoundParameter() {
    canBoundParameter(1);
  }

  @Test
  void canLowerBoundParameter() {
    canBoundParameter(-1);
  }

  private static void canBoundParameter(double value) {
    final ParameterBounds bounds = ParameterBounds.create(new FakeGradientFunction(1, 1, 1));
    if (value < 0) {
      bounds.setBounds(new double[] {2 * value}, null);
    } else {
      bounds.setBounds(null, new double[] {2 * value});
    }
    final double[] a1 = new double[1];
    final double[] a2 = new double[1];
    final double[] step = new double[] {value};
    bounds.applyBounds(a1, step, a2);
    Assertions.assertArrayEquals(a2, new double[] {1 * value}, "Step 1");
    bounds.applyBounds(a2, step, a1);
    Assertions.assertArrayEquals(a1, new double[] {2 * value}, "Step 2");
    bounds.applyBounds(a1, step, a2);
    // Should be bounded
    Assertions.assertArrayEquals(a2, new double[] {2 * value}, "Step 3");
  }

  @Test
  void canDoubleBoundParameter() {
    final ParameterBounds bounds = ParameterBounds.create(new FakeGradientFunction(1, 1, 1));
    final double s = 2;
    bounds.setBounds(new double[] {-s}, new double[] {s});
    final double[] a1 = new double[1];
    final double[] a2 = new double[1];
    bounds.applyBounds(a1, new double[] {10}, a2);
    Assertions.assertArrayEquals(a2, new double[] {s}, "Step 10");

    bounds.applyBounds(a1, new double[] {-10}, a2);
    Assertions.assertArrayEquals(a2, new double[] {-s}, "Step -10");
  }

  @Test
  void canStepParameter() {
    canStepParameter(0.5);
    canStepParameter(1);
    canStepParameter(-0.5);
    canStepParameter(-1);
  }

  private static void canStepParameter(double value) {
    final ParameterBounds bounds = ParameterBounds.create(new FakeGradientFunction(1, 1, 1));
    double[] a1 = new double[1];
    double[] a2 = new double[1];
    double[] tmp;
    final double[] step = new double[] {value};
    final IndexSupplier msg = new IndexSupplier(1, "Step ", null);
    for (int i = 1; i <= 10; i++) {
      bounds.applyBounds(a1, step, a2);
      Assertions.assertArrayEquals(a2, new double[] {i * value}, msg.set(0, i));
      tmp = a1;
      a1 = a2;
      a2 = tmp;
    }
  }
}
