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

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.opentest4j.AssertionFailedError;
import uk.ac.sussex.gdsc.core.utils.BitFlagUtils;

/**
 * Test the ToleranceChecker can converge as expected.
 */
@SuppressWarnings({"javadoc"})
class ToleranceCheckerTest {
  static final double NONE = ToleranceChecker.IGNORE_TOLERANCE;
  static final int IGNORE = ToleranceChecker.IGNORE_MAX_ITERATIONS;

  @Test
  void throwsIfCannotConverge() {
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      canConverge(false, NONE, NONE, NONE, NONE, IGNORE, 0);
    });
  }

  @Test
  void canConvergeOnMaximumRelativeValue() {
    canConverge(false, 1e-2, NONE, NONE, NONE, IGNORE, ToleranceChecker.STATUS_VALUE);
  }

  @Test
  void canConvergeOnMinimumRelativeValue() {
    canConverge(true, 1e-2, NONE, NONE, NONE, IGNORE, ToleranceChecker.STATUS_VALUE);
  }

  @Test
  void canConvergeOnMaximumAbsoluteValue() {
    canConverge(false, NONE, 1e-2, NONE, NONE, IGNORE, ToleranceChecker.STATUS_VALUE);
  }

  @Test
  void canConvergeOnMinimumAbsoluteValue() {
    canConverge(true, NONE, 1e-2, NONE, NONE, IGNORE, ToleranceChecker.STATUS_VALUE);
  }

  @Test
  void cannotConvergeOnMaximumRelativeValueIfMinimising() {
    Assertions.assertThrows(AssertionFailedError.class, () -> {
      canConverge(false, -1, 1e-2, NONE, NONE, NONE, IGNORE, 0);
    });
  }

  @Test
  void cannotConvergeOnMaximumAbsoluteValueIfMinimising() {
    Assertions.assertThrows(AssertionFailedError.class, () -> {
      canConverge(false, -1, NONE, 1e-2, NONE, NONE, IGNORE, 0);
    });
  }

  @Test
  void cannotConvergeOnMinimumRelativeValueIfMaximising() {
    Assertions.assertThrows(AssertionFailedError.class, () -> {
      canConverge(true, 1, 1e-2, NONE, NONE, NONE, IGNORE, 0);
    });
  }

  @Test
  void cannotConvergeOnMinimumAbsoluteValueIfMaximising() {
    Assertions.assertThrows(AssertionFailedError.class, () -> {
      canConverge(true, 1, NONE, 1e-2, NONE, NONE, IGNORE, 0);
    });
  }

  @Test
  void canConvergeOnRelativeParameters() {
    canConverge(true, NONE, NONE, 1e-2, NONE, IGNORE, ToleranceChecker.STATUS_PARAMETERS);
    canConverge(false, NONE, NONE, 1e-2, NONE, IGNORE, ToleranceChecker.STATUS_PARAMETERS);
  }

  @Test
  void canConvergeOnAbsoluteParameters() {
    canConverge(true, NONE, NONE, NONE, 1e-2, IGNORE, ToleranceChecker.STATUS_PARAMETERS);
    canConverge(false, NONE, NONE, NONE, 1e-2, IGNORE, ToleranceChecker.STATUS_PARAMETERS);
  }

  @Test
  void canConvergeOnIterations() {
    canConverge(true, NONE, NONE, NONE, NONE, -10, ToleranceChecker.STATUS_TARGET_ITERATIONS);
    canConverge(false, NONE, NONE, NONE, NONE, -10, ToleranceChecker.STATUS_TARGET_ITERATIONS);
  }

  @Test
  void canConvergeOnMaxIterations() {
    canConverge(true, NONE, NONE, NONE, NONE, 20, ToleranceChecker.STATUS_MAX_ITERATIONS);
    canConverge(false, NONE, NONE, NONE, NONE, 20, ToleranceChecker.STATUS_MAX_ITERATIONS);
  }

  private static void canConverge(boolean minimiseValue, double relativeValue, double absoluteValue,
      double relativeParameters, double absoluteParameters, int maxIterations, int expected) {
    final double dir = (minimiseValue) ? -1 : 1;
    canConverge(minimiseValue, dir, relativeValue, absoluteValue, relativeParameters,
        absoluteParameters, maxIterations, expected);
  }

  private static void canConverge(boolean minimiseValue, double dir, double relativeValue,
      double absoluteValue, double relativeParameters, double absoluteParameters, int maxIterations,
      int expected) {
    final ToleranceChecker tc = new ToleranceChecker(minimiseValue, relativeValue, absoluteValue,
        relativeParameters, absoluteParameters, maxIterations);

    double value = 10;
    double v2 = 1 + value;
    double param = 20;
    double[] p2 = new double[] {1 + param};
    for (int i = 0; i < 20; i++) {
      final double v1 = v2;
      final double[] p1 = p2;
      value *= 0.5;
      param *= 0.5;
      v2 = v1 + dir * value;
      // logger.fine(FormatSupplier.getSupplier("v2 = %f", v2);
      p2 = new double[] {p1[0] + dir * param};
      final int observed = tc.converged(v1, p1, v2, p2);
      if (observed != 0) {
        Assertions.assertEquals(expected, observed);
        if (BitFlagUtils.areSet(expected, ToleranceChecker.STATUS_TARGET_ITERATIONS)) {
          Assertions.assertEquals(-maxIterations, tc.getIterations());
        }
        if (BitFlagUtils.areSet(expected, ToleranceChecker.STATUS_MAX_ITERATIONS)) {
          Assertions.assertEquals(maxIterations, tc.getIterations());
        }
        return;
      }
    }

    Assertions.fail("Failed to converge");
  }

  @Test
  void canConvergeOnImprovedValueIfMaximising() {
    final double tolerance = 1e-2;
    final ToleranceChecker tc = new ToleranceChecker(false, NONE, tolerance, NONE, NONE, 100);
    Assertions.assertEquals(0, tc.converged(0, null, 1, null));
    Assertions.assertEquals(0, tc.converged(0, null, -1, null));
    Assertions.assertEquals(0, tc.converged(0, null, 2 * tolerance, null));
    Assertions.assertEquals(0, tc.converged(0, null, -2 * tolerance, null));
    Assertions.assertEquals(ToleranceChecker.STATUS_VALUE, tc.converged(0, null, tolerance, null));
    Assertions.assertEquals(0, tc.converged(0, null, -tolerance, null));
    Assertions.assertEquals(ToleranceChecker.STATUS_VALUE, tc.converged(0, null, 0, null));
  }

  @Test
  void canConvergeOnImprovedValueIfMinimising() {
    final double tolerance = 1e-2;
    final ToleranceChecker tc = new ToleranceChecker(true, NONE, tolerance, NONE, NONE, 100);
    Assertions.assertEquals(0, tc.converged(0, null, 1, null));
    Assertions.assertEquals(0, tc.converged(0, null, -1, null));
    Assertions.assertEquals(0, tc.converged(0, null, 2 * tolerance, null));
    Assertions.assertEquals(0, tc.converged(0, null, -2 * tolerance, null));
    Assertions.assertEquals(0, tc.converged(0, null, tolerance, null));
    Assertions.assertEquals(ToleranceChecker.STATUS_VALUE, tc.converged(0, null, -tolerance, null));
    Assertions.assertEquals(ToleranceChecker.STATUS_VALUE, tc.converged(0, null, 0, null));
  }

  @Test
  void canConvergeOnValueUsingZeroTolerance() {
    final double tolerance = 0;
    ToleranceChecker tc;

    tc = new ToleranceChecker(false, NONE, NONE, NONE, NONE, 100);
    Assertions.assertFalse(tc.checkValue);
    Assertions.assertFalse(tc.checkParameters);

    tc = new ToleranceChecker(true, tolerance, NONE, NONE, NONE, 100);
    Assertions.assertTrue(tc.checkValue);
    Assertions.assertFalse(tc.checkParameters);
    Assertions.assertEquals(0, tc.converged(0, null, -Double.MIN_VALUE, null));
    Assertions.assertEquals(0, tc.converged(0, null, Double.MIN_VALUE, null));
    Assertions.assertEquals(ToleranceChecker.STATUS_VALUE, tc.converged(0, null, 0, null));

    tc = new ToleranceChecker(false, tolerance, NONE, NONE, NONE, 100);
    Assertions.assertTrue(tc.checkValue);
    Assertions.assertFalse(tc.checkParameters);
    Assertions.assertEquals(0, tc.converged(0, null, -Double.MIN_VALUE, null));
    Assertions.assertEquals(0, tc.converged(0, null, Double.MIN_VALUE, null));
    Assertions.assertEquals(ToleranceChecker.STATUS_VALUE, tc.converged(0, null, 0, null));

    tc = new ToleranceChecker(true, NONE, tolerance, NONE, NONE, 100);
    Assertions.assertTrue(tc.checkValue);
    Assertions.assertFalse(tc.checkParameters);
    Assertions.assertEquals(0, tc.converged(0, null, -Double.MIN_VALUE, null));
    Assertions.assertEquals(0, tc.converged(0, null, Double.MIN_VALUE, null));
    Assertions.assertEquals(ToleranceChecker.STATUS_VALUE, tc.converged(0, null, 0, null));

    tc = new ToleranceChecker(false, NONE, tolerance, NONE, NONE, 100);
    Assertions.assertTrue(tc.checkValue);
    Assertions.assertFalse(tc.checkParameters);
    Assertions.assertEquals(0, tc.converged(0, null, -Double.MIN_VALUE, null));
    Assertions.assertEquals(0, tc.converged(0, null, Double.MIN_VALUE, null));
    Assertions.assertEquals(ToleranceChecker.STATUS_VALUE, tc.converged(0, null, 0, null));
  }

  @Test
  void canConvergeOnParametersUsingZeroTolerance() {
    final double tolerance = 0;
    ToleranceChecker tc;
    final double[] p = new double[1];

    tc = new ToleranceChecker(false, NONE, NONE, NONE, NONE, 100);
    Assertions.assertFalse(tc.checkValue);
    Assertions.assertFalse(tc.checkParameters);

    tc = new ToleranceChecker(true, NONE, NONE, tolerance, NONE, 100);
    Assertions.assertFalse(tc.checkValue);
    Assertions.assertTrue(tc.checkParameters);
    Assertions.assertEquals(0, tc.converged(0, p, 0, new double[] {-Double.MIN_VALUE}));
    Assertions.assertEquals(0, tc.converged(0, p, 0, new double[] {Double.MIN_VALUE}));
    Assertions.assertEquals(ToleranceChecker.STATUS_PARAMETERS, tc.converged(0, p, 0, p));

    tc = new ToleranceChecker(true, NONE, NONE, NONE, tolerance, 100);
    Assertions.assertFalse(tc.checkValue);
    Assertions.assertTrue(tc.checkParameters);
    Assertions.assertEquals(0, tc.converged(0, p, 0, new double[] {-Double.MIN_VALUE}));
    Assertions.assertEquals(0, tc.converged(0, p, 0, new double[] {Double.MIN_VALUE}));
    Assertions.assertEquals(ToleranceChecker.STATUS_PARAMETERS, tc.converged(0, p, 0, p));
  }
}
