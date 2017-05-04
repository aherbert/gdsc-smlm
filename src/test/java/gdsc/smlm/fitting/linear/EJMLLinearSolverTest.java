package gdsc.smlm.fitting.linear;

import java.util.Arrays;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.test.BaseTimingTask;
import gdsc.core.test.TimingService;
import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.TurboList;
import gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator;
import gdsc.smlm.fitting.nonlinear.gradient.GradientCalculatorFactory;
import gdsc.smlm.function.ValueProcedure;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;

public class EJMLLinearSolverTest
{
	@Test
	public void canSolveLinearEquation()
	{
		EJMLLinearSolver solver = new EJMLLinearSolver(3, 1e-6);

		// Solves (one) linear equation, a x = b, for x[n]

		// Taken from https://en.wikipedia.org/wiki/Positive-definite_matrix
		double[][] a = new double[][] { new double[] { 2, -1, 0 }, new double[] { -1, 2, -1 },
				new double[] { 0, -1, 2 } };
		double[] b = new double[] { 3, 3, 4 };

		// Expected solution
		double[] x = new double[] { 4.75, 6.5, 5.25 };
		double[][] a_inv = new double[][] { new double[] { 0.75, 0.5, 0.25 }, new double[] { 0.5, 1, 0.5 },
				new double[] { 0.25, 0.5, 0.75 } };

		boolean result = solver.solve(a, b);
		solver.invert(a);

		Assert.assertTrue("Failed to solve", result);
		Assert.assertArrayEquals("Bad solution", x, b, 1e-4f);

		log("x = %s\n", Arrays.toString(b));
		for (int i = 0; i < b.length; i++)
		{
			log("a[%d] = %s\n", i, Arrays.toString(a[i]));
			Assert.assertArrayEquals("Bad inversion", a_inv[i], a[i], 1e-4f);
		}
	}

	@Test
	public void canSolveLinearEquationWithZeroInB()
	{
		EJMLLinearSolver solver = new EJMLLinearSolver(3, 1e-6);

		// Solves (one) linear equation, a x = b, for x[n]

		// Taken from https://en.wikipedia.org/wiki/Positive-definite_matrix
		double[][] a = new double[][] { new double[] { 2, -1, 0 }, new double[] { -1, 2, -1 },
				new double[] { 0, -1, 2 } };
		double[] b = new double[] { 3, 0, 4 };

		// Expected solution
		double[] x = new double[] { 3.25, 3.5, 3.75 };
		double[][] a_inv = new double[][] { new double[] { 0.75, 0.5, 0.25 }, new double[] { 0.5, 1, 0.5 },
				new double[] { 0.25, 0.5, 0.75 } };

		boolean result = solver.solve(a, b);
		solver.invert(a);

		Assert.assertTrue("Failed to solve", result);
		Assert.assertArrayEquals("Bad solution", x, b, 1e-4f);

		log("x = %s\n", Arrays.toString(b));
		for (int i = 0; i < b.length; i++)
		{
			log("a[%d] = %s\n", i, Arrays.toString(a[i]));
			Assert.assertArrayEquals("Bad inversion", a_inv[i], a[i], 1e-4f);
		}
	}

	@Test
	public void canSolveLinearEquationWithZeroInA()
	{
		EJMLLinearSolver solver = new EJMLLinearSolver(3, 1e-6);

		// Solves (one) linear equation, a x = b, for x[n]

		double[][] a = new double[][] { new double[] { 2, 0, -1, 0 }, new double[] { 0, 0, 0, 0 },
				new double[] { -1, 0, 2, -1 }, new double[] { 0, 0, -1, 2 } };
		double[] b = new double[] { 3, 0, 3, 4 };

		// Expected solution
		double[] x = new double[] { 4.75, 0, 6.5, 5.25 };
		double[][] a_inv = new double[][] { new double[] { 0.75, 0, 0.5, 0.25 }, new double[] { 0, 0, 0, 0 },
				new double[] { 0.5, 0, 1, 0.5 }, new double[] { 0.25, 0, 0.5, 0.75 } };

		boolean result = solver.solve(a, b);
		solver.invert(a);

		Assert.assertTrue("Failed to solve", result);
		Assert.assertArrayEquals("Bad solution", x, b, 1e-4f);

		log("x = %s\n", Arrays.toString(b));
		for (int i = 0; i < b.length; i++)
		{
			log("a[%d] = %s\n", i, Arrays.toString(a[i]));
			Assert.assertArrayEquals("Bad inversion", a_inv[i], a[i], 1e-4f);
		}
	}

	@Test
	public void canSolveLinearEquationWithZerosInA()
	{
		EJMLLinearSolver solver = new EJMLLinearSolver();
		DoubleEquality eq = new DoubleEquality(3, 1e-16);
		solver.setEqual(eq);

		// Solves (one) linear equation, a x = b, for x[n]

		double[][] a = new double[][] { new double[] { 2, 0, -1, 0, 0, 0 }, new double[] { 0, 0, 0, 0, 0, 0 },
				new double[] { -1, 0, 2, 0, 0, -1 }, new double[] { 0, 0, 0, 0, 0, 0 },
				new double[] { 0, 0, 0, 0, 0, 0 }, new double[] { 0, 0, -1, 0, 0, 2 } };
		double[] b = new double[] { 3, 0, 3, 0, 0, 4 };

		// Expected solution
		double[] x = new double[] { 4.75, 0, 6.5, 0, 0, 5.25 };
		double[][] a_inv = new double[][] { new double[] { 0.75, 0, 0.5, 0, 0, 0.25 },
				new double[] { 0, 0, 0, 0, 0, 0 }, new double[] { 0.5, 0, 1, 0, 0, 0.5 },
				new double[] { 0, 0, 0, 0, 0, 0 }, new double[] { 0, 0, 0, 0, 0, 0 },
				new double[] { 0.25, 0, 0.5, 0, 0, 0.75 } };

		boolean result = solver.solve(a, b);
		solver.invert(a);

		Assert.assertTrue("Failed to solve", result);
		Assert.assertArrayEquals("Bad solution", x, b, 1e-4f);

		log("x = %s\n", Arrays.toString(b));
		for (int i = 0; i < b.length; i++)
		{
			log("a[%d] = %s\n", i, Arrays.toString(a[i]));
			Assert.assertArrayEquals("Bad inversion", a_inv[i], a[i], 1e-4f);
		}
	}

	private abstract class SolverTimingTask extends BaseTimingTask
	{
		double[][][] a;
		double[][] b;
		// No validation for a pure speed test
		EJMLLinearSolver solver = new EJMLLinearSolver();

		public SolverTimingTask(String name, double[][][] a, double[][] b)
		{
			super(name + " " + b[0].length);
			// Clone the data
			this.a = a;
			this.b = b;
			// Check the solver gets a good answer
			solver.setEqual(new DoubleEquality(3, 1e-6));
			Object data = getData(0);

			a = (double[][][]) ((Object[]) data)[0];
			b = (double[][]) ((Object[]) data)[1];
			for (int i = 0; i < a.length; i++)
			{
				if (!solve(a[i], b[i]))
				{
					throw new RuntimeException(getName() + " failed to solve");
				}
			}
			solver.setEqual(null);
		}

		public int getSize()
		{
			return 1;
		}

		public Object getData(int i)
		{
			// Clone
			int n = b.length;
			int m = b[0].length;
			double[][][] a = new double[n][][];
			double[][] b = new double[n][];
			while (n-- > 0)
			{
				a[n] = new double[m][];
				for (int j = m; j-- > 0;)
					a[n][j] = this.a[n][j].clone();
				b[n] = this.b[n].clone();
			}
			return new Object[] { a, b };
		}

		public Object run(Object data)
		{
			double[][][] a = (double[][][]) ((Object[]) data)[0];
			double[][] b = (double[][]) ((Object[]) data)[1];
			for (int i = 0; i < a.length; i++)
			{
				solve(a[i], b[i]);
			}
			return null;
		}

		abstract boolean solve(double[][] a, double[] b);
	}

	private class LinearSolverTimingTask extends SolverTimingTask
	{
		public LinearSolverTimingTask(double[][][] a, double[][] b)
		{
			super("Linear", a, b);
		}

		boolean solve(double[][] a, double[] b)
		{
			return solver.solveLinear(a, b);
		}
	}

	private class CholeskySolverTimingTask extends SolverTimingTask
	{
		public CholeskySolverTimingTask(double[][][] a, double[][] b)
		{
			super("Cholesky", a, b);
		}

		boolean solve(double[][] a, double[] b)
		{
			return solver.solveCholesky(a, b);
		}
	}

	private class CholeskyLDLTSolverTimingTask extends SolverTimingTask
	{
		public CholeskyLDLTSolverTimingTask(double[][][] a, double[][] b)
		{
			super("CholeskyLDLT", a, b);
		}

		boolean solve(double[][] a, double[] b)
		{
			return solver.solveCholeskyLDLT(a, b);
		}
	}

	private class PseudoInverseSolverTimingTask extends SolverTimingTask
	{
		public PseudoInverseSolverTimingTask(double[][][] a, double[][] b)
		{
			super("PseudoInverse", a, b);
		}

		boolean solve(double[][] a, double[] b)
		{
			return solver.solvePseudoInverse(a, b);
		}
	}

	private class DirectInversionSolverTimingTask extends SolverTimingTask
	{
		public DirectInversionSolverTimingTask(double[][][] a, double[][] b)
		{
			super("DirectInversion", a, b);
		}

		boolean solve(double[][] a, double[] b)
		{
			return solver.solveDirectInversion(a, b);
		}
	}

	// Create a speed test of the different methods
	@Test
	public void runSpeedTest5()
	{
		runSpeedTest(GaussianFunctionFactory.FIT_ERF_CIRCLE);
	}

	@Test
	public void runSpeedTest4()
	{
		runSpeedTest(GaussianFunctionFactory.FIT_ERF_FIXED);
	}

	@Test
	public void runSpeedTest3()
	{
		runSpeedTest(GaussianFunctionFactory.FIT_NB_FIXED);
	}

	@Test
	public void runSpeedTest2()
	{
		runSpeedTest(GaussianFunctionFactory.FIT_NS_NB_FIXED);
	}

	private void runSpeedTest(int flags)
	{
		final Gaussian2DFunction f0 = GaussianFunctionFactory.create2D(1, 10, 10, flags, null);
		int n = f0.size();
		final double[] y = new double[n];
		final TurboList<double[][]> aList = new TurboList<double[][]>();
		final TurboList<double[]> bList = new TurboList<double[]>();
		double[] testbackground = new double[] { 0.2, 0.7 };
		double[] testsignal1 = new double[] { 30, 100, 300 };
		double[] testcx1 = new double[] { 4.9, 5.3 };
		double[] testcy1 = new double[] { 4.8, 5.2 };
		double[] testw1 = new double[] { 1.1, 1.2, 1.5 };
		int np = f0.getNumberOfGradients();
		GradientCalculator calc = GradientCalculatorFactory.newCalculator(np);
		final RandomDataGenerator rdg = new RandomDataGenerator(new Well19937c(30051977));
		//double lambda = 10;
		for (double background : testbackground)
			// Peak 1
			for (double signal1 : testsignal1)
				for (double cx1 : testcx1)
					for (double cy1 : testcy1)
						for (double w1 : testw1)
						{
							double[] p = new double[] { background, signal1, 0, cx1, cy1, w1, w1 };
							f0.initialise(p);
							f0.forEach(new ValueProcedure()
							{
								int i = 0;

								public void execute(double value)
								{
									// Poisson data 
									y[i++] = rdg.nextPoisson(value);
								}
							});
							double[][] alpha = new double[np][np];
							double[] beta = new double[np];
							//double ss = 
							calc.findLinearised(n, y, p, alpha, beta, f0);
							//System.out.printf("SS = %f\n", ss);
							// As per the LVM algorithm
							//for (int i = 0; i < np; i++)
							//	alpha[i][i] *= lambda;
							aList.add(alpha);
							bList.add(beta);
						}

		double[][][] a = aList.toArray(new double[aList.size()][][]);
		double[][] b = bList.toArray(new double[bList.size()][]);
		int runs = 10000 / a.length;
		TimingService ts = new TimingService(runs);
		ts.execute(new LinearSolverTimingTask(a, b));
		ts.execute(new CholeskySolverTimingTask(a, b));
		ts.execute(new CholeskyLDLTSolverTimingTask(a, b));
		ts.execute(new PseudoInverseSolverTimingTask(a, b));
		ts.execute(new DirectInversionSolverTimingTask(a, b));
		ts.repeat();
		ts.report();
	}

	void log(String format, Object... args)
	{
		System.out.printf(format, args);
	}
}
