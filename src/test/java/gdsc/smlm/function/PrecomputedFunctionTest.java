package gdsc.smlm.function;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.PseudoRandomGenerator;

public class PrecomputedFunctionTest
{
	@Test
	public void precomputedValueFunctionWrapsPrecomputedValues()
	{
		RandomGenerator r = new Well19937c(30051977);
		int size = 100;
		double[] v = new PseudoRandomGenerator(size, r).getSequence();
		ValueFunction f = new PrecomputedValueFunction(v);
		double[] vo = evaluateValueFunction(f);
		Assert.assertArrayEquals("values", v, vo, 0);
	}

	private double[] evaluateValueFunction(ValueFunction f)
	{
		final double[] v = new double[f.size()];
		f.initialise0(null);
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
		int size = 100;
		double[] v = new PseudoRandomGenerator(size, r).getSequence();
		double[][] g1 = new double[size][];
		for (int i = 0; i < g1.length; i++)
		{
			g1[i] = new PseudoRandomGenerator(n, r).getSequence();
		}
		Gradient1Function f = new PrecomputedGradient1Function(v, g1);
		double[][] g1o = new double[size][];
		double[] vo = evaluateGradient1Function(f, g1o);
		Assert.assertArrayEquals("values", v, vo, 0);
		for (int i = 0; i < vo.length; i++)
		{
			Assert.assertArrayEquals("g2", g1[i], g1o[i], 0);
		}
	}

	private double[] evaluateGradient1Function(Gradient1Function f, final double[][] g1)
	{
		final double[] v = new double[f.size()];
		f.initialise1(null);
		f.forEach(new Gradient1Procedure()
		{
			int i = 0;

			public void execute(double value, double[] dy_da)
			{
				v[i] = value;
				g1[i] = dy_da;
				i++;
			}
		});
		return v;
	}

	@Test
	public void precomputedGradient2FunctionWrapsPrecomputedValues()
	{
		int n = 3;
		RandomGenerator r = new Well19937c(30051977);
		int size = 100;
		double[] v = new PseudoRandomGenerator(size, r).getSequence();
		double[][] g1 = new double[size][];
		double[][] g2 = new double[size][];
		for (int i = 0; i < g1.length; i++)
		{
			g1[i] = new PseudoRandomGenerator(n, r).getSequence();
			g2[i] = new PseudoRandomGenerator(n, r).getSequence();
		}
		Gradient2Function f = new PrecomputedGradient2Function(v, g1, g2);
		double[][] g1o = new double[size][];
		double[][] g2o = new double[size][];
		double[] vo = evaluateGradient2Function(f, g1o, g2o);
		Assert.assertArrayEquals("values", v, vo, 0);
		for (int i = 0; i < vo.length; i++)
		{
			Assert.assertArrayEquals("g1", g1[i], g1o[i], 0);
			Assert.assertArrayEquals("g2", g1[i], g1o[i], 0);
		}
	}

	private double[] evaluateGradient2Function(Gradient2Function f, final double[][] g1, final double[][] g2)
	{
		final double[] v = new double[f.size()];
		f.initialise2(null);
		f.forEach(new Gradient2Procedure()
		{
			int i = 0;

			public void execute(double value, double[] dy_da, double[] d2y_da2)
			{
				v[i] = value;
				g1[i] = dy_da;
				g2[i] = d2y_da2;
				i++;
			}
		});
		return v;
	}
}
