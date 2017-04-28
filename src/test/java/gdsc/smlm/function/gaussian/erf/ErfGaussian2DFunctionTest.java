package gdsc.smlm.function.gaussian.erf;

import org.junit.Assert;
import org.junit.Test;

import gdsc.core.test.BaseTimingTask;
import gdsc.core.test.TimingService;
import gdsc.core.utils.BitFlags;
import gdsc.core.utils.TurboList;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.Gaussian2DFunctionTest;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;

public abstract class ErfGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	public ErfGaussian2DFunctionTest()
	{
		super();
		// Test fitting of second derivatives 
		flags |= GaussianFunctionFactory.FIT_2_DERIVATIVES;
	}

	@Test
	public void functionComputesSecondBackgroundGradient()
	{
		if (f1.evaluatesBackground())
			functionComputesSecondTargetGradient(Gaussian2DFunction.BACKGROUND);
	}

	@Test
	public void functionComputesSecondAmplitudeGradient()
	{
		if (f1.evaluatesSignal())
			functionComputesSecondTargetGradient(Gaussian2DFunction.SIGNAL);
	}

	//@Test
	//public void functionComputesSecondShapeGradient()
	//{
	//	if (f1.evaluatesShape())
	//		functionComputesSecondTargetGradient(Gaussian2DFunction.SHAPE);
	//}

	@Test
	public void functionComputesSecondXGradient()
	{
		functionComputesSecondTargetGradient(Gaussian2DFunction.X_POSITION);
	}

	@Test
	public void functionComputesSecondYGradient()
	{
		functionComputesSecondTargetGradient(Gaussian2DFunction.Y_POSITION);
	}

	@Test
	public void functionComputesSecondXWidthGradient()
	{
		if (f1.evaluatesSD0())
			functionComputesSecondTargetGradient(Gaussian2DFunction.X_SD);
	}

	@Test
	public void functionComputesSecondYWidthGradient()
	{
		if (f1.evaluatesSD1())
			functionComputesSecondTargetGradient(Gaussian2DFunction.Y_SD);
	}

	private void functionComputesSecondTargetGradient(int targetParameter)
	{
		int gradientIndex = findGradientIndex(f1, targetParameter);
		double[] dyda = new double[f1.gradientIndices().length];
		double[] dyda1 = new double[dyda.length];
		double[] dyda2 = new double[dyda.length];
		double[] a;

		ErfGaussian2DFunction f1a = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(1, maxx, maxx, flags);
		ErfGaussian2DFunction f1b = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(1, maxx, maxx, flags);
		System.out.printf("functionComputesSecondTargetGradient %s %s\n", f1.getClass().getSimpleName(),
				f1.getName(targetParameter));

		for (double background : testbackground)
			// Peak 1
			for (double amplitude1 : testamplitude1)
				for (double shape1 : testshape1)
					for (double cx1 : testcx1)
						for (double cy1 : testcy1)
							for (double[] w1 : testw1)
							{
								a = createParameters(background, amplitude1, shape1, cx1, cy1, w1[0], w1[1]);
								f1.initialise(a);

								// Numerically solve gradient. 
								// Calculate the step size h to be an exact numerical representation
								final double xx = a[targetParameter];

								// Get h to minimise roundoff error
								double h = h_; //((xx == 0) ? 1 : xx) * h_;
								final double temp = xx + h;
								doNothing(temp);
								h = temp - xx;

								// Evaluate at (x+h) and (x-h)
								a = createParameters(background, amplitude1, shape1, cx1, cy1, w1[0], w1[1]);
								a[targetParameter] = xx + h;
								f1a.initialise(a);

								a = createParameters(background, amplitude1, shape1, cx1, cy1, w1[0], w1[1]);
								a[targetParameter] = xx - h;
								f1b.initialise(a);

								for (int x : testx)
									for (int y : testy)
									{
										int i = y * maxx + x;
										f1a.eval(i, dyda1, dyda2);
										double value2 = dyda1[gradientIndex];
										f1b.eval(i, dyda1, dyda2);
										double value3 = dyda1[gradientIndex];
										f1.eval(i, dyda, dyda2);

										double gradient = (value2 - value3) / (2 * h);
										//System.out.printf("[%d,%d] %f == [%d] %f?\n", x, y, gradient, gradientIndex, dyda2[gradientIndex]);
										Assert.assertTrue(gradient + " != " + dyda2[gradientIndex],
												eq.almostEqualComplement(gradient, dyda2[gradientIndex]));
									}
							}
	}

	private class MyTimingTask extends BaseTimingTask
	{
		Gaussian2DFunction f;
		double[][] x;
		int order;
		final double[] dyda;
		final int n = maxx * maxx;

		public MyTimingTask(Gaussian2DFunction f, double[][] x, int order)
		{
			super(f.getClass().getSimpleName() + " " + order);
			this.f = f;
			this.x = x;
			this.order = order;
			dyda = new double[f.gradientIndices().length];
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
			double s = 0;
			if (order == 0)
			{
				for (int i = 0; i < x.length; i++)
				{
					f.initialise(x[i]);
					for (int j = 0; j < n; j++)
						s += f.eval(j);
				}
			}
			else
			{
				for (int i = 0; i < x.length; i++)
				{
					f.initialise(x[i]);
					for (int j = 0; j < n; j++)
						s += f.eval(j, dyda);
				}
			}
			return s;
		}
	}

	// Speed test verses equivalent Gaussian2DFunction
	@Test
	public void functionIsFasterThanEquivalentGaussian2DFunction()
	{
		final ErfGaussian2DFunction f2 = (ErfGaussian2DFunction) this.f1.create(2);
		final ErfGaussian2DFunction f1 = (ErfGaussian2DFunction) this.f1.create(1);
		final ErfGaussian2DFunction f0 = (ErfGaussian2DFunction) this.f1.create(0);
		int flags = BitFlags.unset(this.flags, GaussianFunctionFactory.FIT_ERF);
		final Gaussian2DFunction gf = GaussianFunctionFactory.create2D(1, maxx, maxx, flags);

		final TurboList<double[]> params = new TurboList<double[]>();
		for (double background : testbackground)
			// Peak 1
			for (double amplitude1 : testamplitude1)
				for (double shape1 : testshape1)
					for (double cx1 : testcx1)
						for (double cy1 : testcy1)
							for (double[] w1 : testw1)
								params.add(createParameters(background, amplitude1, shape1, cx1, cy1, w1[0], w1[1]));
		double[][] x = params.toArray(new double[params.size()][]);

		int runs = 10000 / x.length;
		TimingService ts = new TimingService(runs);
		ts.execute(new MyTimingTask(gf, x, 1));
		ts.execute(new MyTimingTask(gf, x, 0));
		ts.execute(new MyTimingTask(f2, x, 2));
		ts.execute(new MyTimingTask(f1, x, 1));
		ts.execute(new MyTimingTask(f0, x, 0));

		int size = ts.getSize();
		ts.repeat(size);
		ts.report();
	}

	// TODO - Computation of widths given z. This requires support for astigmatism.

	// 
}
