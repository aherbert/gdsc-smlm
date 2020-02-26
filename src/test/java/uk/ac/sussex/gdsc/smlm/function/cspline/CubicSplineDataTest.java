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

package uk.ac.sussex.gdsc.smlm.function.cspline;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import uk.ac.sussex.gdsc.core.math.interpolation.CustomTricubicFunction;
import uk.ac.sussex.gdsc.core.math.interpolation.CustomTricubicFunctionUtils;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;

@SuppressWarnings({"javadoc"})
public class CubicSplineDataTest {
  @SeededTest
  public void canExternaliseDoubleFunction(RandomSeed seed) throws IOException {
    canExternaliseFunction(seed, false);
  }

  @SeededTest
  public void canExternaliseFloatFunction(RandomSeed seed) throws IOException {
    canExternaliseFunction(seed, true);
  }

  private static void canExternaliseFunction(RandomSeed seed, boolean singlePrecision)
      throws IOException {
    final UniformRandomProvider r = RngUtils.create(seed.getSeed());
    final int x = 6;
    final int y = 5;
    final int z = 4;

    final int size = x * y;
    final CustomTricubicFunction[][] splines = new CustomTricubicFunction[z][x * y];
    final double[] a = new double[64];
    for (int zz = 0; zz < z; zz++) {
      for (int i = 0; i < size; i++) {
        for (int j = 0; j < 64; j++) {
          a[j] = r.nextDouble();
        }
        splines[zz][i] = CustomTricubicFunctionUtils.create(a);
        if (singlePrecision) {
          splines[zz][i] = splines[zz][i].toSinglePrecision();
        }
      }
    }
    final CubicSplineData f1 = new CubicSplineData(x, y, splines);

    final ByteArrayOutputStream b = new ByteArrayOutputStream();
    f1.write(b);

    final byte[] bytes = b.toByteArray();
    final CubicSplineData f2 = CubicSplineData.read(new ByteArrayInputStream(bytes));

    final double[] exp = new double[64];
    final double[] obs = new double[64];
    for (int zz = 0; zz < z; zz++) {
      for (int i = 0; i < size; i++) {
        f1.splines[zz][i].getCoefficients(exp);
        f2.splines[zz][i].getCoefficients(obs);
        Assertions.assertArrayEquals(exp, obs);
      }
    }
  }
}
