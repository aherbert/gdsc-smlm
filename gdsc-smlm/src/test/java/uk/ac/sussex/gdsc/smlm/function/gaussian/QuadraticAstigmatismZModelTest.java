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

package uk.ac.sussex.gdsc.smlm.function.gaussian;

import java.util.logging.Logger;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils.TestLevel;

@SuppressWarnings({"javadoc"})
class QuadraticAstigmatismZModelTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(QuadraticAstigmatismZModelTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  protected DoubleEquality eq = new DoubleEquality(1e-5, 1e-7);

  // Compute as per Numerical Recipes 5.7.
  // Approximate error accuracy in single precision: Ef
  // Step size for derivatives:
  // h ~ (Ef)^(1/3) * xc
  // xc is the characteristic scale over which x changes, assumed to be 1 (not x as per NR since x
  // is close to zero)
  protected double stepH = 0.0001; // (double) (Math.pow(1e-3f, 1.0 / 3));

  @Test
  void canStaticComputeGradient() {
    canStaticComputeGradient(1.2);
    canStaticComputeGradient(0.3);
  }

  private void canStaticComputeGradient(double zDepth) {
    final double[] ds_dz = new double[1];
    final double[] ds_dz2 = new double[2];
    final double[] ds_duz = new double[1];
    final double[] ds_dlz = new double[1];
    final boolean record = logger.isLoggable(TestLevel.TEST_INFO);
    for (double z = -0.5; z < 0.5; z += 0.01) {
      final double s0 = QuadraticAstigmatismZModel.getS(z, zDepth);
      final double s1 = QuadraticAstigmatismZModel.getS1(z, zDepth, ds_dz);
      final double s2 = QuadraticAstigmatismZModel.getS2(z, zDepth, ds_dz2);

      Assertions.assertEquals(s0, s1);
      Assertions.assertEquals(s0, s2);
      Assertions.assertEquals(ds_dz[0], ds_dz2[0]);

      final double uz = z + stepH;
      final double lz = z - stepH;
      final double upper = QuadraticAstigmatismZModel.getS1(uz, zDepth, ds_duz);
      final double lower = QuadraticAstigmatismZModel.getS1(lz, zDepth, ds_dlz);

      final double e1 = (upper - lower) / (uz - lz);
      final double o1 = ds_dz[0];

      // Second gradient
      final double e2 = (ds_duz[0] - ds_dlz[0]) / (uz - lz);
      final double o2 = ds_dz2[1];

      if (record) {
        logger.log(TestLogUtils.getRecord(TestLevel.TEST_INFO,
            "z=%f s=%f : ds_dz=%g  %g  (%g): d2s_dz2=%g   %g  (%g)", z, s0, e1, o1,
            DoubleEquality.relativeError(o1, e1), e2, o2, DoubleEquality.relativeError(o2, e2)));
      }

      // double error = DoubleEquality.relativeError(o, e);
      if (Math.abs(z) > 0.02) {
        Assertions.assertTrue((e1 * o1) >= 0, () -> String.format("%s sign != %s", e1, o1));
      }
      Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(e1, o1),
          () -> String.format("%s != %s", e1, o1));

      if (Math.abs(z) > 0.02) {
        Assertions.assertTrue((e2 * o2) >= 0, () -> String.format("%s sign != %s", e2, o2));
      }
      Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(e2, o2),
          () -> String.format("%s != %s", e2, o2));
    }
  }
}
