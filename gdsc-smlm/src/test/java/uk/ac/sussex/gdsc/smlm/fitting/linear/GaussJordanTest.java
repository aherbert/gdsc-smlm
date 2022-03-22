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

package uk.ac.sussex.gdsc.smlm.fitting.linear;

import java.util.Arrays;
import java.util.logging.Logger;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.test.utils.TestLogging.TestLevel;

@SuppressWarnings({"javadoc"})
class GaussJordanTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(GaussJordanTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  @Test
  void canSolveLinearEquation() {
    final GaussJordan solver = new GaussJordan();

    // Solves (one) linear equation, a x = b, for x[n]

    // Example taken from http://eigen.tuxfamily.org/dox/TutorialLinearAlgebra.html
    final float[][] a =
        new float[][] {new float[] {1, 2, 3}, new float[] {4, 5, 6}, new float[] {7, 8, 10}};
    final float[] b = new float[] {3, 3, 4};
    final float[] expecteds = new float[] {-2, 1, 1};

    final boolean result = solver.solve(a, b);

    Assertions.assertTrue(result);
    Assertions.assertArrayEquals(expecteds, b, 1e-4f);

    if (logger.isLoggable(TestLevel.TEST_INFO)) {
      logger.log(TestLevel.TEST_INFO, () -> String.format("x = %s", Arrays.toString(b)));
      for (int i = 0; i < b.length; i++) {
        final int ii = i;
        logger.log(TestLevel.TEST_INFO,
            () -> String.format("a[%d] = %s", ii, Arrays.toString(a[ii])));
      }
    }
  }
}
