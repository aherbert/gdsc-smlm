package gdsc.smlm.function;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.PseudoRandomGenerator;
import gdsc.core.utils.SimpleArrayUtils;

public class PrecomputedFunctionTest
{
	@Test
	public void precomputedValueFunctionWrapsPrecomputedValues()
	{
		int n = 3;
		RandomGenerator r = new Well19937c(30051977);
		ValueFunction f0 = new FakeGradientFunction(3, n);
		int size = f0.size();
		double[] b1 = new PseudoRandomGenerator(size, r).getSequence();
		double[] b2 = new PseudoRandomGenerator(size, r).getSequence();
		ValueFunction f1 = PrecomputedValueFunction.wrapValueFunction(f0, b1);
		ValueFunction f2 = PrecomputedValueFunction.wrapValueFunction(f1, b2);
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

	private double[] evaluateValueFunction(ValueFunction f, double[] p)
	{
		f.initialise0(p);
		final double[] v = new double[f.size()];
		f.forEach(new ValueProcedure()
		{
			int i = 0;

			public void execute(double value)
			{
				v[i++] = value;
			}
		});
		return v;
	}

	@Test
	public void precomputedGradient1FunctionWrapsPrecomputedValues()
	{
		int n = 3;
		RandomGenerator r = new Well19937c(30051977);
		Gradient1Function f0 = new FakeGradientFunction(3, n);
		int size = f0.size();
		double[] b1 = new PseudoRandomGenerator(size, r).getSequence();
		double[] b2 = new PseudoRandomGenerator(size, r).getSequence();
		Gradient1Function f1 = PrecomputedGradient1Function.wrapGradient1Function(f0, b1);
		Gradient1Function f2 = PrecomputedGradient1Function.wrapGradient1Function(f1, b2);
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

	private double[] evaluateGradient1Function(Gradient1Function f, double[] p, final double[] dyda)
	{
		f.initialise0(p);
		final double[] v = new double[f.size()];
		f.forEach(new Gradient1Procedure()
		{
			int i = 0;

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
	public void precomputedGradient2FunctionWrapsPrecomputedValues()
	{
		int n = 3;
		RandomGenerator r = new Well19937c(30051977);
		Gradient2Function f0 = new FakeGradientFunction(3, n);
		int size = f0.size();
		double[] b1 = new PseudoRandomGenerator(size, r).getSequence();
		double[] b2 = new PseudoRandomGenerator(size, r).getSequence();
		Gradient2Function f1 = PrecomputedGradient2Function.wrapGradient2Function(f0, b1);
		Gradient2Function f2 = PrecomputedGradient2Function.wrapGradient2Function(f1, b2);
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

	private double[] evaluateGradient2Function(Gradient2Function f, double[] p, final double[] dyda,
			final double[] d2yda2)
	{
		f.initialise0(p);
		final double[] v = new double[f.size()];
		f.forEach(new Gradient2Procedure()
		{
			int i = 0;

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
	public void precomputedValueFunctionCanWrapPrecomputed()
	{
		double[] a = new double[] { 3.2, 5.6 };
		FakeGradientFunction f = new FakeGradientFunction(10, a.length);
		double[] b = SimpleArrayUtils.newArray(f.size(), 1.0, 0);
		StandardValueProcedure sp = new StandardValueProcedure();
		double[] e = sp.getValues(f, a).clone();
		ValueFunction vf = f;
		for (int n = 0; n < 3; n++)
		{
			vf = PrecomputedValueFunction.wrapValueFunction(vf, b);
			double[] o = sp.getValues(vf, a);
			for (int i = 0; i < e.length; i++)
				e[i] += b[i];
			Assert.assertArrayEquals(e, o, 0);
			Assert.assertTrue(((PrecomputedValueFunction) vf).getValueFunction() == f);
		}
	}

	@Test
	public void precomputedGradient1FunctionCanWrapPrecomputed()
	{
		double[] a = new double[] { 3.2, 5.6 };
		FakeGradientFunction f = new FakeGradientFunction(10, a.length);
		double[] b = SimpleArrayUtils.newArray(f.size(), 1.0, 0);
		StandardValueProcedure sp = new StandardValueProcedure();
		double[] e = sp.getValues(f, a).clone();
		Gradient1Function vf = f;
		for (int n = 0; n < 3; n++)
		{
			vf = PrecomputedGradient1Function.wrapGradient1Function(vf, b);
			double[] o = sp.getValues(vf, a);
			for (int i = 0; i < e.length; i++)
				e[i] += b[i];
			Assert.assertArrayEquals(e, o, 0);
			Assert.assertTrue(((PrecomputedGradient1Function) vf).getGradient1Function() == f);
		}
	}

	@Test
	public void precomputedGradient2FunctionCanWrapPrecomputed()
	{
		double[] a = new double[] { 3.2, 5.6 };
		FakeGradientFunction f = new FakeGradientFunction(10, a.length);
		double[] b = SimpleArrayUtils.newArray(f.size(), 1.0, 0);
		StandardValueProcedure sp = new StandardValueProcedure();
		double[] e = sp.getValues(f, a).clone();
		Gradient2Function vf = f;
		for (int n = 0; n < 3; n++)
		{
			vf = PrecomputedGradient2Function.wrapGradient2Function(vf, b);
			double[] o = sp.getValues(vf, a);
			for (int i = 0; i < e.length; i++)
				e[i] += b[i];
			Assert.assertArrayEquals(e, o, 0);
			Assert.assertTrue(((PrecomputedGradient2Function) vf).getGradient2Function() == f);
		}
	}

	@Test
	public void precomputedExtendedGradient2FunctionCanWrapPrecomputed()
	{
		double[] a = new double[] { 3.2, 5.6 };
		FakeGradientFunction f = new FakeGradientFunction(10, a.length);
		double[] b = SimpleArrayUtils.newArray(f.size(), 1.0, 0);
		StandardValueProcedure sp = new StandardValueProcedure();
		double[] e = sp.getValues(f, a).clone();
		ExtendedGradient2Function vf = f;
		for (int n = 0; n < 3; n++)
		{
			vf = PrecomputedExtendedGradient2Function.wrapExtendedGradient2Function(vf, b);
			double[] o = sp.getValues(vf, a);
			for (int i = 0; i < e.length; i++)
				e[i] += b[i];
			Assert.assertArrayEquals(e, o, 0);
			Assert.assertTrue(((PrecomputedExtendedGradient2Function) vf).getExtendedGradient2Function() == f);
		}
	}
}
