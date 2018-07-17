/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
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
package gdsc.smlm.function;

import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.PseudoRandomGenerator;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.test.TestSettings;

@SuppressWarnings({ "javadoc" })
public class OffsetFunctionTest
{
	@Test
	public void offsetValueFunctionWrapsPrecomputedValues()
	{
		int n = 3;
		RandomGenerator r = TestSettings.getRandomGenerator();
		ValueFunction f0 = new FakeGradientFunction(3, n);
		int size = f0.size();
		double[] b1 = new PseudoRandomGenerator(size, r).getSequence();
		double[] b2 = new PseudoRandomGenerator(size, r).getSequence();
		ValueFunction f1 = OffsetValueFunction.wrapValueFunction(f0, b1);
		ValueFunction f2 = OffsetValueFunction.wrapValueFunction(f1, b2);
		double[] p = new double[n];
		for (int i = 0; i < n; i++)
			p[i] = r.nextDouble();
		double[] v0 = evaluateValueFunction(f0, p);
		double[] v1 = evaluateValueFunction(f1, p);
		double[] v2 = evaluateValueFunction(f2, p);
		for (int i = 0; i < v0.length; i++)
		{
			double e = v0[i] + b1[i] + b2[i];
			double o1 = v1[i] + b2[i];
			double o2 = v2[i];
			Assert.assertEquals("o1", e, o1, 0);
			Assert.assertEquals("o2", e, o2, 1e-6);
		}
	}

	private static double[] evaluateValueFunction(ValueFunction f, double[] p)
	{
		f.initialise0(p);
		final double[] v = new double[f.size()];
		f.forEach(new ValueProcedure()
		{
			int i = 0;

			@Override
			public void execute(double value)
			{
				v[i++] = value;
			}
		});
		return v;
	}

	@Test
	public void offsetGradient1FunctionWrapsPrecomputedValues()
	{
		int n = 3;
		RandomGenerator r = TestSettings.getRandomGenerator();
		Gradient1Function f0 = new FakeGradientFunction(3, n);
		int size = f0.size();
		double[] b1 = new PseudoRandomGenerator(size, r).getSequence();
		double[] b2 = new PseudoRandomGenerator(size, r).getSequence();
		Gradient1Function f1 = OffsetGradient1Function.wrapGradient1Function(f0, b1);
		Gradient1Function f2 = OffsetGradient1Function.wrapGradient1Function(f1, b2);
		double[] p = new double[n];
		for (int i = 0; i < n; i++)
			p[i] = r.nextDouble();
		double[] d0 = new double[n];
		double[] d1 = new double[n];
		double[] d2 = new double[n];
		double[] v0 = evaluateGradient1Function(f0, p, d0);
		double[] v1 = evaluateGradient1Function(f1, p, d1);
		double[] v2 = evaluateGradient1Function(f2, p, d2);
		for (int i = 0; i < v0.length; i++)
		{
			double e = v0[i] + b1[i] + b2[i];
			double o1 = v1[i] + b2[i];
			double o2 = v2[i];
			Assert.assertEquals("o1", e, o1, 0);
			Assert.assertEquals("o2", e, o2, 1e-6);
		}
		Assert.assertArrayEquals("d1", d0, d1, 0);
		Assert.assertArrayEquals("d2", d0, d2, 0);
	}

	private static double[] evaluateGradient1Function(Gradient1Function f, double[] p, final double[] dyda)
	{
		f.initialise0(p);
		final double[] v = new double[f.size()];
		f.forEach(new Gradient1Procedure()
		{
			int i = 0;

			@Override
			public void execute(double value, double[] dy_da)
			{
				v[i++] = value;
				for (int j = 0; j < dy_da.length; j++)
					dyda[j] += dy_da[j];
			}
		});
		return v;
	}

	@Test
	public void offsetGradient2FunctionWrapsPrecomputedValues()
	{
		int n = 3;
		RandomGenerator r = TestSettings.getRandomGenerator();
		Gradient2Function f0 = new FakeGradientFunction(3, n);
		int size = f0.size();
		double[] b1 = new PseudoRandomGenerator(size, r).getSequence();
		double[] b2 = new PseudoRandomGenerator(size, r).getSequence();
		Gradient2Function f1 = OffsetGradient2Function.wrapGradient2Function(f0, b1);
		Gradient2Function f2 = OffsetGradient2Function.wrapGradient2Function(f1, b2);
		double[] p = new double[n];
		for (int i = 0; i < n; i++)
			p[i] = r.nextDouble();
		double[] d0 = new double[n];
		double[] d1 = new double[n];
		double[] d2 = new double[n];
		double[] d20 = new double[n];
		double[] d21 = new double[n];
		double[] d22 = new double[n];
		double[] v0 = evaluateGradient2Function(f0, p, d0, d20);
		double[] v1 = evaluateGradient2Function(f1, p, d1, d21);
		double[] v2 = evaluateGradient2Function(f2, p, d2, d22);
		for (int i = 0; i < v0.length; i++)
		{
			double e = v0[i] + b1[i] + b2[i];
			double o1 = v1[i] + b2[i];
			double o2 = v2[i];
			Assert.assertEquals("o1", e, o1, 0);
			Assert.assertEquals("o2", e, o2, 1e-6);
		}
		Assert.assertArrayEquals("d1", d0, d1, 0);
		Assert.assertArrayEquals("d2", d0, d2, 0);
		Assert.assertArrayEquals("d21", d20, d21, 0);
		Assert.assertArrayEquals("d22", d20, d22, 0);
	}

	private static double[] evaluateGradient2Function(Gradient2Function f, double[] p, final double[] dyda,
			final double[] d2yda2)
	{
		f.initialise0(p);
		final double[] v = new double[f.size()];
		f.forEach(new Gradient2Procedure()
		{
			int i = 0;

			@Override
			public void execute(double value, double[] dy_da, double[] d2y_da2)
			{
				v[i++] = value;
				for (int j = 0; j < dy_da.length; j++)
				{
					dyda[j] += dy_da[j];
					d2yda2[j] += d2y_da2[j];
				}
			}
		});
		return v;
	}

	@Test
	public void offsetValueFunctionCanWrapPrecomputed()
	{
		double[] a = new double[] { 3.2, 5.6 };
		FakeGradientFunction f = new FakeGradientFunction(10, a.length);
		double[] b = SimpleArrayUtils.newArray(f.size(), 1.0, 0);
		StandardValueProcedure sp = new StandardValueProcedure();
		double[] e = sp.getValues(f, a).clone();
		ValueFunction vf = f;
		for (int n = 0; n < 3; n++)
		{
			vf = OffsetValueFunction.wrapValueFunction(vf, b);
			double[] o = sp.getValues(vf, a);
			for (int i = 0; i < e.length; i++)
				e[i] += b[i];
			Assert.assertArrayEquals(e, o, 0);
			Assert.assertTrue(((OffsetValueFunction) vf).getValueFunction() == f);
		}
	}

	@Test
	public void offsetGradient1FunctionCanWrapPrecomputed()
	{
		double[] a = new double[] { 3.2, 5.6 };
		FakeGradientFunction f = new FakeGradientFunction(10, a.length);
		double[] b = SimpleArrayUtils.newArray(f.size(), 1.0, 0);
		StandardValueProcedure sp = new StandardValueProcedure();
		double[] e = sp.getValues(f, a).clone();
		Gradient1Function vf = f;
		for (int n = 0; n < 3; n++)
		{
			vf = OffsetGradient1Function.wrapGradient1Function(vf, b);
			double[] o = sp.getValues(vf, a);
			for (int i = 0; i < e.length; i++)
				e[i] += b[i];
			Assert.assertArrayEquals(e, o, 0);
			Assert.assertTrue(((OffsetGradient1Function) vf).getGradient1Function() == f);
		}
	}

	@Test
	public void offsetGradient2FunctionCanWrapPrecomputed()
	{
		double[] a = new double[] { 3.2, 5.6 };
		FakeGradientFunction f = new FakeGradientFunction(10, a.length);
		double[] b = SimpleArrayUtils.newArray(f.size(), 1.0, 0);
		StandardValueProcedure sp = new StandardValueProcedure();
		double[] e = sp.getValues(f, a).clone();
		Gradient2Function vf = f;
		for (int n = 0; n < 3; n++)
		{
			vf = OffsetGradient2Function.wrapGradient2Function(vf, b);
			double[] o = sp.getValues(vf, a);
			for (int i = 0; i < e.length; i++)
				e[i] += b[i];
			Assert.assertArrayEquals(e, o, 0);
			Assert.assertTrue(((OffsetGradient2Function) vf).getGradient2Function() == f);
		}
	}

	@Test
	public void offsetExtendedGradient2FunctionCanWrapPrecomputed()
	{
		double[] a = new double[] { 3.2, 5.6 };
		FakeGradientFunction f = new FakeGradientFunction(10, a.length);
		double[] b = SimpleArrayUtils.newArray(f.size(), 1.0, 0);
		StandardValueProcedure sp = new StandardValueProcedure();
		double[] e = sp.getValues(f, a).clone();
		ExtendedGradient2Function vf = f;
		for (int n = 0; n < 3; n++)
		{
			vf = OffsetExtendedGradient2Function.wrapExtendedGradient2Function(vf, b);
			double[] o = sp.getValues(vf, a);
			for (int i = 0; i < e.length; i++)
				e[i] += b[i];
			Assert.assertArrayEquals(e, o, 0);
			Assert.assertTrue(((OffsetExtendedGradient2Function) vf).getExtendedGradient2Function() == f);
		}
	}
}
