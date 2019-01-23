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

package uk.ac.sussex.gdsc.smlm.function;

import uk.ac.sussex.gdsc.core.utils.PseudoRandomGenerator;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;

@SuppressWarnings({"javadoc"})
public class OffsetFunctionTest {
  @SeededTest
  public void offsetValueFunctionWrapsPrecomputedValues(RandomSeed seed) {
    final int n = 3;
    final UniformRandomProvider r = RngUtils.create(seed.getSeedAsLong());
    final ValueFunction f0 = new FakeGradientFunction(3, n);
    final int size = f0.size();
    final double[] b1 = new PseudoRandomGenerator(size, r).getSequence();
    final double[] b2 = new PseudoRandomGenerator(size, r).getSequence();
    final ValueFunction f1 = OffsetValueFunction.wrapValueFunction(f0, b1);
    final ValueFunction f2 = OffsetValueFunction.wrapValueFunction(f1, b2);
    final double[] p = new double[n];
    for (int i = 0; i < n; i++) {
      p[i] = r.nextDouble();
    }
    final double[] v0 = evaluateValueFunction(f0, p);
    final double[] v1 = evaluateValueFunction(f1, p);
    final double[] v2 = evaluateValueFunction(f2, p);
    for (int i = 0; i < v0.length; i++) {
      final double e = v0[i] + b1[i] + b2[i];
      final double o1 = v1[i] + b2[i];
      final double o2 = v2[i];
      Assertions.assertEquals(e, o1, "o1");
      Assertions.assertEquals(e, o2, 1e-6, "o2");
    }
  }

  private static double[] evaluateValueFunction(ValueFunction func, double[] params) {
    func.initialise0(params);
    final double[] v = new double[func.size()];
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
  public void offsetGradient1FunctionWrapsPrecomputedValues(RandomSeed seed) {
    final int n = 3;
    final UniformRandomProvider r = RngUtils.create(seed.getSeedAsLong());
    final Gradient1Function f0 = new FakeGradientFunction(3, n);
    final int size = f0.size();
    final double[] b1 = new PseudoRandomGenerator(size, r).getSequence();
    final double[] b2 = new PseudoRandomGenerator(size, r).getSequence();
    final Gradient1Function f1 = OffsetGradient1Function.wrapGradient1Function(f0, b1);
    final Gradient1Function f2 = OffsetGradient1Function.wrapGradient1Function(f1, b2);
    final double[] p = new double[n];
    for (int i = 0; i < n; i++) {
      p[i] = r.nextDouble();
    }
    final double[] d0 = new double[n];
    final double[] d1 = new double[n];
    final double[] d2 = new double[n];
    final double[] v0 = evaluateGradient1Function(f0, p, d0);
    final double[] v1 = evaluateGradient1Function(f1, p, d1);
    final double[] v2 = evaluateGradient1Function(f2, p, d2);
    for (int i = 0; i < v0.length; i++) {
      final double e = v0[i] + b1[i] + b2[i];
      final double o1 = v1[i] + b2[i];
      final double o2 = v2[i];
      Assertions.assertEquals(e, o1, "o1");
      Assertions.assertEquals(e, o2, 1e-6, "o2");
    }
    Assertions.assertArrayEquals(d0, d1, "d1");
    Assertions.assertArrayEquals(d0, d2, "d2");
  }

  private static double[] evaluateGradient1Function(Gradient1Function func, double[] params,
      final double[] dyda) {
    func.initialise0(params);
    final double[] v = new double[func.size()];
    func.forEach(new Gradient1Procedure() {
      int index = 0;

      @Override
      public void execute(double value, double[] dyDa) {
        v[index++] = value;
        for (int j = 0; j < dyDa.length; j++) {
          dyda[j] += dyDa[j];
        }
      }
    });
    return v;
  }

  @SeededTest
  public void offsetGradient2FunctionWrapsPrecomputedValues(RandomSeed seed) {
    final int n = 3;
    final UniformRandomProvider r = RngUtils.create(seed.getSeedAsLong());
    final Gradient2Function f0 = new FakeGradientFunction(3, n);
    final int size = f0.size();
    final double[] b1 = new PseudoRandomGenerator(size, r).getSequence();
    final double[] b2 = new PseudoRandomGenerator(size, r).getSequence();
    final Gradient2Function f1 = OffsetGradient2Function.wrapGradient2Function(f0, b1);
    final Gradient2Function f2 = OffsetGradient2Function.wrapGradient2Function(f1, b2);
    final double[] p = new double[n];
    for (int i = 0; i < n; i++) {
      p[i] = r.nextDouble();
    }
    final double[] d0 = new double[n];
    final double[] d1 = new double[n];
    final double[] d2 = new double[n];
    final double[] d20 = new double[n];
    final double[] d21 = new double[n];
    final double[] d22 = new double[n];
    final double[] v0 = evaluateGradient2Function(f0, p, d0, d20);
    final double[] v1 = evaluateGradient2Function(f1, p, d1, d21);
    final double[] v2 = evaluateGradient2Function(f2, p, d2, d22);
    for (int i = 0; i < v0.length; i++) {
      final double e = v0[i] + b1[i] + b2[i];
      final double o1 = v1[i] + b2[i];
      final double o2 = v2[i];
      Assertions.assertEquals(e, o1, "o1");
      Assertions.assertEquals(e, o2, 1e-6, "o2");
    }
    Assertions.assertArrayEquals(d0, d1, "d1");
    Assertions.assertArrayEquals(d0, d2, "d2");
    Assertions.assertArrayEquals(d20, d21, "d21");
    Assertions.assertArrayEquals(d20, d22, "d22");
  }

  private static double[] evaluateGradient2Function(Gradient2Function func, double[] params,
      final double[] dyda, final double[] d2yda2) {
    func.initialise0(params);
    final double[] v = new double[func.size()];
    func.forEach(new Gradient2Procedure() {
      int index = 0;

      @Override
      public void execute(double value, double[] dyDa, double[] d2yDa2) {
        v[index++] = value;
        for (int j = 0; j < dyDa.length; j++) {
          dyda[j] += dyDa[j];
          d2yda2[j] += d2yDa2[j];
        }
      }
    });
    return v;
  }

  @SeededTest
  public void offsetValueFunctionCanWrapPrecomputed() {
    final double[] a = new double[] {3.2, 5.6};
    final FakeGradientFunction f = new FakeGradientFunction(10, a.length);
    final double[] b = SimpleArrayUtils.newArray(f.size(), 1.0, 0);
    final StandardValueProcedure sp = new StandardValueProcedure();
    final double[] e = sp.getValues(f, a).clone();
    ValueFunction vf = f;
    for (int n = 0; n < 3; n++) {
      vf = OffsetValueFunction.wrapValueFunction(vf, b);
      final double[] o = sp.getValues(vf, a);
      for (int i = 0; i < e.length; i++) {
        e[i] += b[i];
      }
      Assertions.assertArrayEquals(e, o);
      Assertions.assertTrue(((OffsetValueFunction) vf).getValueFunction() == f);
    }
  }

  @SeededTest
  public void offsetGradient1FunctionCanWrapPrecomputed() {
    final double[] a = new double[] {3.2, 5.6};
    final FakeGradientFunction f = new FakeGradientFunction(10, a.length);
    final double[] b = SimpleArrayUtils.newArray(f.size(), 1.0, 0);
    final StandardValueProcedure sp = new StandardValueProcedure();
    final double[] e = sp.getValues(f, a).clone();
    Gradient1Function vf = f;
    for (int n = 0; n < 3; n++) {
      vf = OffsetGradient1Function.wrapGradient1Function(vf, b);
      final double[] o = sp.getValues(vf, a);
      for (int i = 0; i < e.length; i++) {
        e[i] += b[i];
      }
      Assertions.assertArrayEquals(e, o);
      Assertions.assertTrue(((OffsetGradient1Function) vf).getGradient1Function() == f);
    }
  }

  @SeededTest
  public void offsetGradient2FunctionCanWrapPrecomputed() {
    final double[] a = new double[] {3.2, 5.6};
    final FakeGradientFunction f = new FakeGradientFunction(10, a.length);
    final double[] b = SimpleArrayUtils.newArray(f.size(), 1.0, 0);
    final StandardValueProcedure sp = new StandardValueProcedure();
    final double[] e = sp.getValues(f, a).clone();
    Gradient2Function vf = f;
    for (int n = 0; n < 3; n++) {
      vf = OffsetGradient2Function.wrapGradient2Function(vf, b);
      final double[] o = sp.getValues(vf, a);
      for (int i = 0; i < e.length; i++) {
        e[i] += b[i];
      }
      Assertions.assertArrayEquals(e, o);
      Assertions.assertTrue(((OffsetGradient2Function) vf).getGradient2Function() == f);
    }
  }

  @SeededTest
  public void offsetExtendedGradient2FunctionCanWrapPrecomputed() {
    final double[] a = new double[] {3.2, 5.6};
    final FakeGradientFunction f = new FakeGradientFunction(10, a.length);
    final double[] b = SimpleArrayUtils.newArray(f.size(), 1.0, 0);
    final StandardValueProcedure sp = new StandardValueProcedure();
    final double[] e = sp.getValues(f, a).clone();
    ExtendedGradient2Function vf = f;
    for (int n = 0; n < 3; n++) {
      vf = OffsetExtendedGradient2Function.wrapExtendedGradient2Function(vf, b);
      final double[] o = sp.getValues(vf, a);
      for (int i = 0; i < e.length; i++) {
        e[i] += b[i];
      }
      Assertions.assertArrayEquals(e, o);
      Assertions
          .assertTrue(((OffsetExtendedGradient2Function) vf).getExtendedGradient2Function() == f);
    }
  }
}
