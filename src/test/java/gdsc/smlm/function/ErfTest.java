package gdsc.smlm.function;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Precision;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.test.BaseTimingTask;
import gdsc.core.test.TimingService;
import gdsc.core.utils.DoubleEquality;

public class ErfTest
{
	//@formatter:off
	private abstract class BaseErf
	{
		String name;
		BaseErf(String name) { this.name = name; }
		abstract double erf(double x);
		abstract double erf(double x1, double x2);
	}
	private class ApacheErf extends BaseErf
	{
		ApacheErf() {	super("apache erf"); }
		double erf(double x) { return org.apache.commons.math3.special.Erf.erf(x); }
		double erf(double x1, double x2) { return org.apache.commons.math3.special.Erf.erf(x1, x2); }
	}
	private class Erf extends BaseErf
	{
		Erf() {	super("erf"); }
		double erf(double x) { return gdsc.smlm.function.Erf.erf(x); }
		double erf(double x1, double x2) { return gdsc.smlm.function.Erf.erf(x1, x2); }
	}
	private class Erf0 extends BaseErf
	{
		Erf0() { super("erf0"); }
		double erf(double x) { return gdsc.smlm.function.Erf.erf0(x); }
		double erf(double x1, double x2) { return gdsc.smlm.function.Erf.erf0(x1, x2); }
	}
	private class Erf2 extends BaseErf
	{
		Erf2() { super("erf2"); }
		double erf(double x) { return gdsc.smlm.function.Erf.erf2(x); }
		double erf(double x1, double x2) { return gdsc.smlm.function.Erf.erf2(x1, x2); }
	}
	//@formatter:on

	@Test
	public void erf0xHasLowError()
	{
		erfxHasLowError(new Erf0(), 5e-4);
	}

	@Test
	public void erfxHasLowError()
	{
		erfxHasLowError(new Erf(), 3e-7);
	}

	@Test
	public void erf2xHasLowError()
	{
		erfxHasLowError(new Erf2(), 1.3e-4);
	}

	private void erfxHasLowError(BaseErf erf, double expected)
	{
		RandomGenerator rg = new Well19937c(30051977);
		int range = 8;
		double max = 0;

		for (int xi = -range; xi <= range; xi++)
		{
			for (int i = 0; i < 5; i++)
			{
				double x = xi + rg.nextDouble();
				double o = erf.erf(x);
				double e = org.apache.commons.math3.special.Erf.erf(x);
				double error = Math.abs(o - e);
				if (max < error)
					max = error;
				//System.out.printf("x=%f, e=%f, o=%f, error=%f\n", x, e, o, error);
				Assert.assertTrue(error < expected);
			}
		}
		System.out.printf("erfx %s max error = %g\n", erf.name, max);
	}

	@Test
	public void erfApachexIndistinguishableFrom1()
	{
		erfxIndistinguishableFrom1(new ApacheErf());
	}

	@Test
	public void erf0xIndistinguishableFrom1()
	{
		erfxIndistinguishableFrom1(new Erf0());
	}

	@Test
	public void erfxIndistinguishableFrom1()
	{
		erfxIndistinguishableFrom1(new Erf());
	}

	@Test
	public void erf2xIndistinguishableFrom1()
	{
		erfxIndistinguishableFrom1(new Erf2());
	}

	private void erfxIndistinguishableFrom1(BaseErf erf)
	{
		// Find switch using a binary search
		double lower = 1;
		double upper = 40;
		while (DoubleEquality.complement(lower, upper) > 1)
		{
			double mid = (upper + lower) * 0.5;
			double o = erf.erf(mid);
			if (o == 1)
			{
				upper = mid;
			}
			else
			{
				lower = mid;
			}
		}

		System.out.printf("erfx %s indistinguishable from 1: x > %s, x >= %s\n", erf.name, Double.toString(lower),
				Double.toString(upper));
	}

	@Test
	public void erf0xxHasLowError()
	{
		erfxxHasLowError(new Erf0(), 4e-2);
	}

	@Test
	public void erfxxHasLowError()
	{
		erfxxHasLowError(new Erf(), 7e-4);
	}

	@Test
	public void erf2xxHasLowError()
	{
		erfxxHasLowError(new Erf2(), 1.1e-2);
	}

	private void erfxxHasLowError(BaseErf erf, double expected)
	{
		RandomGenerator rg = new Well19937c(30051977);

		int range = 3;
		double max = 0;

		for (int xi = -range; xi <= range; xi++)
		{
			for (int xi2 = -range; xi2 <= range; xi2++)
			{
				for (int i = 0; i < 5; i++)
				{
					double x = xi + rg.nextDouble();
					for (int j = 0; j < 5; j++)
					{
						double x2 = xi2 + rg.nextDouble();

						double o = erf.erf(x, x2);
						double e = org.apache.commons.math3.special.Erf.erf(x, x2);
						double error = Math.abs(o - e);
						if (max < error)
							max = error;
						//System.out.printf("x=%f, x2=%f, e=%f, o=%f, error=%f\n", x, x2, e, o, error);
						Assert.assertTrue(error < expected);
					}
				}
			}
		}

		System.out.printf("erfxx %s max error = %g\n", erf.name, max);
	}

	@Test
	public void erf0xxHasLowErrorForUnitBlocks()
	{
		erfxxHasLowErrorForUnitBlocks(new Erf0(), 5e-4);
	}

	@Test
	public void erfxxHasLowErrorForUnitBlocks()
	{
		erfxxHasLowErrorForUnitBlocks(new Erf(), 5e-7);
	}

	@Test
	public void erf2xxHasLowErrorForUnitBlocks()
	{
		erfxxHasLowErrorForUnitBlocks(new Erf2(), 1e-4);
	}

	private void erfxxHasLowErrorForUnitBlocks(BaseErf erf, double expected)
	{
		int range = 8;
		double max = 0;

		for (int xi = -range; xi <= range; xi++)
		{
			double x = xi;
			double x2 = xi + 1;
			double o = erf.erf(x, x2);
			double e = org.apache.commons.math3.special.Erf.erf(x, x2);
			double error = Math.abs(o - e);
			if (max < error)
				max = error;
			//System.out.printf("x=%f, x2=%f, e=%f, o=%f, error=%f\n", x, x2, e, o, error);
			Assert.assertTrue(error < expected);
		}

		System.out.printf("erfxx %s unit max error = %g\n", erf.name, max);
	}

	@Test
	public void erf0xxHasLowerErrorThanGaussianApproximationForUnitBlocks()
	{
		erfxxHasLowerErrorThanGaussianApproximationForUnitBlocks(new Erf0());
	}

	@Test
	public void erfxxHasLowerErrorThanGaussianApproximationForUnitBlocks()
	{
		erfxxHasLowerErrorThanGaussianApproximationForUnitBlocks(new Erf());
	}

	@Test
	public void erf2xxHasLowerErrorThanGaussianApproximationForUnitBlocks()
	{
		erfxxHasLowerErrorThanGaussianApproximationForUnitBlocks(new Erf2());
	}

	private void erfxxHasLowerErrorThanGaussianApproximationForUnitBlocks(BaseErf erf)
	{
		int range = 5;
		double max = 0, max2 = 0;

		// Standard deviation
		double s = 1.3;
		final double twos2 = 2 * s * s;
		double norm = 1 / (Math.PI * twos2);
		final double denom = 1.0 / (Math.sqrt(2.0) * s);

		double sum1 = 0, sum2 = 0, sum3 = 0;

		for (int x = -range; x <= range; x++)
		{
			double o1 = 0.5 * erf.erf((x - 0.5) * denom, (x + 0.5) * denom);
			double e1 = 0.5 * org.apache.commons.math3.special.Erf.erf((x - 0.5) * denom, (x + 0.5) * denom);
			for (int y = -range; y <= range; y++)
			{
				double o2 = 0.5 * erf.erf((y - 0.5) * denom, (y + 0.5) * denom);
				double e2 = 0.5 * org.apache.commons.math3.special.Erf.erf((y - 0.5) * denom, (y + 0.5) * denom);

				double o = o1 * o2;
				double e = e1 * e2;
				double oo = norm * FastMath.exp(-(x * x + y * y) / twos2);

				sum1 += e;
				sum2 += o;
				sum3 += oo;

				double absError = Math.abs(o - e);
				if (e < 1e-4 || absError < 1e-10)
					continue;
				double error = DoubleEquality.relativeError(o, e);
				double error2 = DoubleEquality.relativeError(oo, e);
				if (max < error)
					max = error;
				if (max2 < error2)
					max2 = error2;
				//System.out.printf("x=%d, y=%d, e=%g, o=%g, o2=%g, error=%f, error2=%f\n", x, y, e, o, oo, error, error2);
				Assert.assertTrue(error < error2);
			}
		}

		Assert.assertTrue(erf.name + " Gaussian 2D integral is not 1", sum1 > 0.999);
		Assert.assertTrue(erf.name + " Erf approx integral is incorrect",
				DoubleEquality.relativeError(sum1, sum2) < 1e-3);
		Assert.assertTrue(erf.name + " Gaussian approx integral is incorrect",
				DoubleEquality.relativeError(sum1, sum3) < 1e-3);

		System.out.printf(erf.name + " Erf approx pixel unit max error = %f\n", max);
		System.out.printf(erf.name + " Gaussian approx pixel unit max error = %f\n", max2);
	}

	private class ErfTimingTask extends BaseTimingTask
	{
		BaseErf erf;
		double[] x;

		public ErfTimingTask(BaseErf erf, double[] x)
		{
			super(erf.name);
			this.erf = erf;
			this.x = x;
		}

		public int getSize()
		{
			return 1;
		}

		public Object getData(int i)
		{
			return null;
		}

		public Object run(Object data)
		{
			for (int i = 0; i < x.length; i++)
				erf.erf(x[i]);
			return null;
		}
	}

	@Test
	public void erfApproxIsFaster()
	{
		int range = 5;
		int steps = 10000;
		final double[] x = new double[steps];
		double total = 2 * range;
		double step = total / steps;
		for (int i = 0; i < steps; i++)
			x[i] = -range + i * step;

		TimingService ts = new TimingService(5);
		ts.execute(new ErfTimingTask(new ApacheErf(), x));
		ts.execute(new ErfTimingTask(new Erf(), x));
		ts.execute(new ErfTimingTask(new Erf0(), x));
		ts.execute(new ErfTimingTask(new Erf2(), x));

		int size = ts.getSize();
		ts.repeat(size);
		ts.report();
	}

	@Test
	public void gaussianIntegralApproximatesErf()
	{
		double x = 1.3, y = 2.2, s = 1.14;
		int minx = (int) x;
		int miny = (int) y;
		int maxx = minx + 1;
		int maxy = miny + 1;

		// Full integration using the Erf
		// Note: The PSF of a 2D Gaussian is described in Smith et all using a denominator
		// of (2.0 * s * s) for both x and Y directions. This is wrong. We need the 
		// integral of the single Guassian in each dimension so the denomiator is (sqrt(2.0) * s). 
		// See: Smith et al, (2010). Fast, single-molecule localisation that achieves 
		// theoretically minimum uncertainty. Nature Methods 7, 373-375
		// (supplementary note).
		//final double denom = 1.0 / (2.0 * s * s); // As per Smith, etal (2010),

		final double denom = 1.0 / (Math.sqrt(2.0) * s);
		double e1 = 0.5 * org.apache.commons.math3.special.Erf.erf(minx * denom, maxx * denom);
		double e2 = 0.5 * org.apache.commons.math3.special.Erf.erf(miny * denom, maxy * denom);
		double e = e1 * e2;

		double o = 0;
		// Numeric integration
		final double twos2 = 2 * s * s;
		double norm = 1 / (Math.PI * twos2);
		for (int i = 0, steps = 1; i < 4; i++, steps = (int) FastMath.pow(10, i))
		{
			// Gaussian is: FastMath.exp(-(x * x + y * y) / twos2) over all x and y
			// But we can do this by separating x and y:
			// FastMath.exp(-(x * x) / twos2) * FastMath.exp(-(y * y) / twos2)

			// pre-compute
			double[] ex = new double[steps];
			double sumey = 0;
			if (steps == 1)
			{
				// Use the actual values for x and y
				ex[0] = FastMath.exp(-(x * x) / twos2);
				sumey = FastMath.exp(-(y * y) / twos2);
			}
			else
			{
				for (int j = 0; j < steps; j++)
				{
					double xx = minx + (double) j / steps;
					double yy = miny + (double) j / steps;
					ex[j] = FastMath.exp(-(xx * xx) / twos2);
					sumey += FastMath.exp(-(yy * yy) / twos2);
				}
			}

			double sum = 0;
			for (int j = 0; j < steps; j++)
			{
				sum += ex[j] * sumey;
			}

			//// Check
			//double sum2 = 0;
			//for (int j = 0; j <= steps; j++)
			//{
			//	double xx = minx + (double) j / steps;
			//	for (int k = 0; k <= steps; k++)
			//	{
			//		double yy = miny + (double) k / steps;
			//		sum2 += FastMath.exp(-(xx * xx + yy * yy) / twos2);
			//	}
			//}
			//System.out.printf("sum=%f, sum2=%f\n", sum, sum2);

			int n = steps * steps;
			o = norm * sum / n;
			System.out.printf("n=%d, e=%f, o=%f, error=%f\n", n, e, o, DoubleEquality.relativeError(e, o));
		}

		Assert.assertEquals(e, o, e * 1e-2);
	}

	@Test
	public void analyticErfGradientCorrectForErfApproximation()
	{
		BaseErf erf = new Erf();
		int range = 7;
		int steps = 10000;
		double step = (double) range / steps;
		double delta = 1e-3;
		DoubleEquality eq = new DoubleEquality(5e-4, 1e-6);
		for (int i = 0; i < steps; i++)
		{
			double x = i * step;
			double x1 = x + Precision.representableDelta(x, delta);
			double x2 = x - Precision.representableDelta(x, delta);
			double o1 = erf.erf(x1);
			double o2 = erf.erf(x2);
			double delta2 = x1 - x2;
			double g = (o1 - o2) / delta2;
			double e = gdsc.smlm.function.Erf.dErf_dx(x);
			if (!eq.almostEqualRelativeOrAbsolute(e, g))
				Assert.assertTrue(x + " : " + e + " != " + g, false);
		}
	}

	@Test
	public void canComputePower4()
	{
		for (int i = -10; i <= 10; i++)
		{
			for (double d : new double[] { 0, 0.1, 0.01, 0.001 })
			{
				double f = i + d;
				double e = Math.pow(f, 4);
				double o = gdsc.smlm.function.Erf.pow4(f);
				Assert.assertEquals("x=" + f, e, o, e * 1e-10);
			}
		}
	}

	@Test
	public void canComputePower16()
	{
		for (int i = -10; i <= 10; i++)
		{
			for (double d : new double[] { 0, 0.1, 0.01, 0.001 })
			{
				double f = i + d;
				double e = Math.pow(f, 16);
				double o = gdsc.smlm.function.Erf.pow16(f);
				Assert.assertEquals("x=" + f, e, o, e * 1e-10);
			}
		}
	}

	// See if power functions are faster

	//@formatter:off
	private abstract class BasePow
	{
		String name;
		BasePow(String name) { this.name = name; }
		abstract double pow(double x);
	}
	private class MathPow4 extends BasePow
	{
		MathPow4() {	super("Math pow4"); }
		double pow(double x) { return Math.pow(x, 4); }
	}
	private class FastMathPow4 extends BasePow
	{
		FastMathPow4() {	super("FastMath pow4"); }
		double pow(double x) { return FastMath.pow(x, 4L); }
	}
	private class Pow4 extends BasePow
	{
		Pow4() {	super("pow4"); }
		double pow(double x) { return gdsc.smlm.function.Erf.pow4(x); }
	}
	private class MathPow16 extends BasePow
	{
		MathPow16() {	super("Math pow16"); }
		double pow(double x) { return Math.pow(x, 16); }
	}
	private class FastMathPow16 extends BasePow
	{
		FastMathPow16() {	super("FastMath pow16"); }
		double pow(double x) { return FastMath.pow(x, 16); }
	}
	private class Pow16 extends BasePow
	{
		Pow16() {	super("pow16"); }
		double pow(double x) { return gdsc.smlm.function.Erf.pow16(x); }
	}
	//@formatter:on

	private class PowTimingTask extends BaseTimingTask
	{
		BasePow pow;
		double[] x;

		public PowTimingTask(BasePow pow, double[] x)
		{
			super(pow.name);
			this.pow = pow;
			this.x = x;
		}

		public int getSize()
		{
			return 1;
		}

		public Object getData(int i)
		{
			return null;
		}

		public Object run(Object data)
		{
			for (int i = 0; i < x.length; i++)
				pow.pow(x[i]);
			return null;
		}
	}

	@Test
	public void powerApproxIsFaster()
	{
		int range = 5000;
		int steps = 100000;
		final double[] x = new double[steps];
		double step = range / steps;
		for (int i = 0; i < steps; i++)
			x[i] = i * step;

		TimingService ts = new TimingService(5);
		ts.execute(new PowTimingTask(new MathPow4(), x));
		ts.execute(new PowTimingTask(new FastMathPow4(), x));
		ts.execute(new PowTimingTask(new Pow4(), x));
		ts.execute(new PowTimingTask(new MathPow16(), x));
		ts.execute(new PowTimingTask(new FastMathPow16(), x));
		ts.execute(new PowTimingTask(new Pow16(), x));

		int size = ts.getSize();
		ts.repeat(size);
		ts.report();

		for (int i = 0; i < 2; i++)
		{
			int j = -(1 + i * 3);
			Assert.assertTrue(ts.get(j).getMean() < ts.get(j - 1).getMean());
			Assert.assertTrue(ts.get(j).getMean() < ts.get(j - 2).getMean());
		}
	}
}
