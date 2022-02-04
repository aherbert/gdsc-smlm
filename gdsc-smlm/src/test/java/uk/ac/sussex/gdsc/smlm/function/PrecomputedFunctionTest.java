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

package uk.ac.sussex.gdsc.smlm.function;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import uk.ac.sussex.gdsc.smlm.GdscSmlmTestUtils;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;

@SuppressWarnings({"javadoc"})
class PrecomputedFunctionTest {
  @SeededTest
  void precomputedValueFunctionWrapsPrecomputedValues(RandomSeed seed) {
    final UniformRandomProvider r = RngUtils.create(seed.getSeed());
    final int size = 100;
    final double[] v = GdscSmlmTestUtils.generateDoubles(size, r);
    final ValueFunction func = new PrecomputedValueFunction(v);
    final double[] vo = evaluateValueFunction(func);
    Assertions.assertArrayEquals(v, vo, "values");
  }

  private static double[] evaluateValueFunction(ValueFunction func) {
    final double[] v = new double[func.size()];
    func.initialise0(null);
    func.forEach(new ValueProcedure() {
      int index = 0;

      @Override
      public void execute(double value) {
        v[index++] = value;
      }
    });
    return v;
  }

  @SeededTest
  void precomputedGradient1FunctionWrapsPrecomputedValues(RandomSeed seed) {
    final int n = 3;
    final UniformRandomProvider r = RngUtils.create(seed.getSeed());
    final int size = 100;
    final double[] v = GdscSmlmTestUtils.generateDoubles(size, r);
    final double[][] g1 = new double[size][];
    for (int i = 0; i < g1.length; i++) {
      g1[i] = GdscSmlmTestUtils.generateDoubles(n, r);
    }
    final Gradient1Function func = new PrecomputedGradient1Function(v, g1);
    final double[][] g1o = new double[size][];
    final double[] vo = evaluateGradient1Function(func, g1o);
    Assertions.assertArrayEquals(v, vo, "values");
    Assertions.assertArrayEquals(g1, g1o, "g1");
  }

  private static double[] evaluateGradient1Function(Gradient1Function func, final double[][] g1) {
    final double[] v = new double[func.size()];
    func.initialise1(null);
    func.forEach(new Gradient1Procedure() {
      int index = 0;

      @Override
      public void execute(double value, double[] dyDa) {
        v[index] = value;
        g1[index] = dyDa;
        index++;
      }
    });
    return v;
  }

  @SeededTest
  void precomputedGradient2FunctionWrapsPrecomputedValues(RandomSeed seed) {
    final int n = 3;
    final UniformRandomProvider r = RngUtils.create(seed.getSeed());
    final int size = 100;
    final double[] v = GdscSmlmTestUtils.generateDoubles(size, r);
    final double[][] g1 = new double[size][];
    final double[][] g2 = new double[size][];
    for (int i = 0; i < g1.length; i++) {
      g1[i] = GdscSmlmTestUtils.generateDoubles(n, r);
      g2[i] = GdscSmlmTestUtils.generateDoubles(n, r);
    }
    final Gradient2Function func = new PrecomputedGradient2Function(v, g1, g2);
    final double[][] g1o = new double[size][];
    final double[][] g2o = new double[size][];
    final double[] vo = evaluateGradient2Function(func, g1o, g2o);
    Assertions.assertArrayEquals(v, vo, "values");
    Assertions.assertArrayEquals(g1, g1o, "g1");
    Assertions.assertArrayEquals(g2, g2o, "g2");
  }

  private static double[] evaluateGradient2Function(Gradient2Function func, final double[][] g1,
      final double[][] g2) {
    final double[] v = new double[func.size()];
    func.initialise2(null);
    func.forEach(new Gradient2Procedure() {
      int index = 0;

      @Override
      public void execute(double value, double[] dyDa, double[] d2yDa2) {
        v[index] = value;
        g1[index] = dyDa;
        g2[index] = d2yDa2;
        index++;
      }
    });
    return v;
  }
}
