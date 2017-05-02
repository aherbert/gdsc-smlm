package gdsc.smlm.function.gaussian.erf;

import org.apache.commons.math3.util.Pair;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.ij.Utils;
import gdsc.core.test.BaseTimingTask;
import gdsc.core.test.TimingService;
import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.Statistics;
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
		//flags |= GaussianFunctionFactory.FIT_2_DERIVATIVES;

		// The derivative check can be tighter with the ERF since it is a true integration
		h_ = 0.0001;
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

	@Test
	public void functionComputesSecondShapeGradient()
	{
		if (f1.evaluatesShape())
			functionComputesSecondTargetGradient(Gaussian2DFunction.SHAPE);
	}

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
		double[] dyda2 = new double[dyda.length];
		double[] a;

		// Test fitting of second derivatives 
		int flags = this.flags | GaussianFunctionFactory.FIT_2_DERIVATIVES;
		ErfGaussian2DFunction f1a = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(1, maxx, maxx, flags,
				zModel);
		ErfGaussian2DFunction f1b = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(1, maxx, maxx, flags,
				zModel);
		Statistics s = new Statistics();

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
										f1a.eval(i, dyda, dyda2);
										double value2 = dyda[gradientIndex];
										f1b.eval(i, dyda, dyda2);
										double value3 = dyda[gradientIndex];
										f1.eval(i, dyda, dyda2);

										double gradient = (value2 - value3) / (2 * h);
										double error = DoubleEquality.relativeError(gradient, dyda2[gradientIndex]);
										s.add(error);
										Assert.assertTrue(gradient + " sign != " + dyda2[gradientIndex],
												(gradient * dyda2[gradientIndex]) >= 0);
										//System.out.printf("[%d,%d] %f == [%d] %f? (%g)\n", x, y, gradient,
										//		gradientIndex, dyda2[gradientIndex], error);
										Assert.assertTrue(gradient + " != " + dyda2[gradientIndex],
												eq.almostEqualComplement(gradient, dyda2[gradientIndex]));
									}
							}
		System.out.printf("functionComputesSecondTargetGradient %s %s (error %s +/- %s)\n",
				f1.getClass().getSimpleName(), f1.getName(targetParameter), Utils.rounded(s.getMean()),
				Utils.rounded(s.getStandardDeviation()));
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
		int flags = this.flags & ~GaussianFunctionFactory.FIT_ERF;
		final Gaussian2DFunction gf = GaussianFunctionFactory.create2D(1, maxx, maxx, flags, zModel);

		boolean zDepth = (flags & GaussianFunctionFactory.FIT_Z) != 0;

		final TurboList<double[]> params = new TurboList<double[]>();
		final TurboList<double[]> params2 = new TurboList<double[]>();
		for (double background : testbackground)
			// Peak 1
			for (double amplitude1 : testamplitude1)
				for (double shape1 : testshape1)
					for (double cx1 : testcx1)
						for (double cy1 : testcy1)
							for (double[] w1 : testw1)
							{
								double[] a = createParameters(background, amplitude1, shape1, cx1, cy1, w1[0], w1[1]);
								params.add(a);
								if (zDepth)
								{
									// Change to a standard free circular function
									a = a.clone();
									a[Gaussian2DFunction.X_SD] *= zModel.getSx(a[Gaussian2DFunction.SHAPE]);
									a[Gaussian2DFunction.Y_SD] *= zModel.getSy(a[Gaussian2DFunction.SHAPE]);
									a[Gaussian2DFunction.SHAPE] = 0;
									params2.add(a);
								}
							}
		double[][] x = params.toArray(new double[params.size()][]);
		double[][] x2 = (zDepth) ? params2.toArray(new double[params2.size()][]) : x;

		int runs = 10000 / x.length;
		TimingService ts = new TimingService(runs);
		ts.execute(new MyTimingTask(gf, x2, 1));
		ts.execute(new MyTimingTask(gf, x2, 0));
		ts.execute(new MyTimingTask(f2, x, 2));
		ts.execute(new MyTimingTask(f1, x, 1));
		ts.execute(new MyTimingTask(f0, x, 0));

		int size = ts.getSize();
		ts.repeat(size);
		ts.report();

		int n = ts.getSize() - 1;
		Assert.assertTrue("0 order", ts.get(n).getMean() < ts.get(n - 3).getMean());
		n--;
		Assert.assertTrue("1 order", ts.get(n).getMean() < ts.get(n - 3).getMean());
	}

	// Test that the value and jacobian is correct since this is re-implemented 
	@Test
	public void functionComputesValueAndJacobian()
	{
		final ErfGaussian2DFunction f1 = (ErfGaussian2DFunction) this.f1.create(1);
		
		final int n = maxx * maxx;
		double[] du_da = new double[f1.gradientIndices().length];
		
		for (double background : testbackground)
			// Peak 1
			for (double amplitude1 : testamplitude1)
				for (double shape1 : testshape1)
					for (double cx1 : testcx1)
						for (double cy1 : testcy1)
							for (double[] w1 : testw1)
							{
								double[] a = createParameters(background, amplitude1, shape1, cx1, cy1, w1[0], w1[1]);
								f1.initialise(a);
								double[] values = f1.computeValues(a);
								double[][] jacobian = f1.computeJacobian(a);
								Pair<double[], double[][]> pair = f1.computeValuesAndJacobian(a);
								Assert.assertArrayEquals("Values!=Values from ValuesAndJacobian", values, pair.getFirst(), 1e-10);
								double[][] jacobian2 = pair.getSecond();
								for (int i=0; i<n; i++)
								{
									Assert.assertArrayEquals("Jacobian!=Jacobian from ValuesAndJacobian", jacobian[i], jacobian2[i], 1e-10);
									double e = f1.eval(i, du_da);
									Assert.assertEquals("Value!=Values", e, values[i], 1e-10);
									Assert.assertArrayEquals("Jacobian!=Jacobians", jacobian[i], du_da, 1e-10);
								}
							}
	}
	
	// Speed test value and jacobian 
}
