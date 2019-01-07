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

package uk.ac.sussex.gdsc.smlm.utils;

import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

@SuppressWarnings({"javadoc"})
public class TensorTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(TensorTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  @Test
  public void canComputeTensor3D() {
    //@formatter:off
    final float[][] data = new float[][] {
      { 2, 1, 0, 1, 2, 1, 0, 1, 2 },
      //{ 1, 0, 0, 0, 1, 0, 0, 0, 1 },
      //{ 1, 0, 0, 0, 1, 0, 0, 0, 1 },
    };
    //@formatter:on
    final Tensor3D t = new Tensor3D(data, 3, 3);
    Assertions.assertTrue(t.hasDecomposition());
    final double[] com = t.getCentreOfMass();
    Assertions.assertArrayEquals(new double[] {1, 1, 0}, com);
    final double[] v = t.getEigenValues();
    final double[][] vv = t.getEigenVectors();
    print(com, v, vv);
    for (int i = 1; i < v.length; i++) {
      Assertions.assertTrue(v[i - 1] >= v[i]);
    }
  }

  private static void print(double[] com, double[] v, double[][] vv) {
    if (logger.isLoggable(Level.INFO)) {
      final StringBuilder sb = new StringBuilder();
      final String newLine = System.getProperty("line.separator");
      sb.append(String.format("%scom = %s", newLine, Arrays.toString(com)));
      for (int i = 0; i < v.length; i++) {
        sb.append(String.format("%s[%d] %f = %s  %.2f", newLine, i, v[i],
            java.util.Arrays.toString(vv[i]), 180.0 * Math.atan2(vv[i][1], vv[i][0]) / Math.PI));
      }
      logger.info(sb.toString());
    }
  }

  @Test
  public void canComputeTensor2D() {
    //@formatter:off
    // Line through [0][0], [1][1], [2][2]
    // longest axis of object is -45 degrees
    final float[] data = new float[] {
        //1, 0, 0, 0, 1, 0, 0, 0, 1
        2, 1, 0, 1, 2, 1, 0, 1, 2
        //2, 0, 0, 0, 0, 0, 0, 0, 2
        };
    //@formatter:on
    final Tensor2D t = new Tensor2D(data, 3, 3);
    Assertions.assertTrue(t.hasDecomposition());
    final double[] com = t.getCentreOfMass();
    Assertions.assertArrayEquals(new double[] {1, 1}, com);
    final double[] v = t.getEigenValues();
    final double[][] vv = t.getEigenVectors();
    print(com, v, vv);
    for (int i = 1; i < v.length; i++) {
      Assertions.assertTrue(v[i - 1] >= v[i]);
    }
  }

  @SeededTest
  public void canComputeSameTensor(RandomSeed seed) {
    final UniformRandomProvider random = RngUtils.create(seed.getSeedAsLong());
    final int w = 3;
    final int h = 4;
    final float[] data = new float[w * h];
    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-6, 0);
    for (int i = 0; i < 10; i++) {
      for (int j = data.length; j-- > 0;) {
        data[j] = random.nextFloat();
      }

      final Tensor2D t2 = new Tensor2D(data, w, h);
      final Tensor3D t3 = new Tensor3D(new float[][] {data}, w, h);

      final double[] com2 = t2.getCentreOfMass();
      final double[] v2 = t2.getEigenValues();
      final double[][] vv2 = t2.getEigenVectors();
      final double[] com3 = t3.getCentreOfMass();
      final double[] v3 = t3.getEigenValues();
      final double[][] vv3 = t3.getEigenVectors();
      for (int k = 0; k < 2; k++) {
        Assertions.assertEquals(com2[k], com3[k]);
        TestAssertions.assertTest(v2[k], v3[k + 1], predicate);
        for (int kk = 0; kk < 2; kk++) {
          // Swap vector direction
          if (Math.signum(vv2[k][kk]) != Math.signum(vv3[k + 1][kk])) {
            vv2[k][0] = -vv2[k][0];
            vv2[k][1] = -vv2[k][1];
          }
          TestAssertions.assertTest(vv2[k][kk], vv3[k + 1][kk], predicate);
        }
      }
    }
  }
}
