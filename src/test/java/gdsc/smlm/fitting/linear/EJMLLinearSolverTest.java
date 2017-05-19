package gdsc.smlm.fitting.linear;

import java.util.Arrays;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.ejml.data.DenseMatrix64F;
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
	//@formatter:off
	@Test
	public void canSolveLinearEquation()
	{
		EJMLLinearSolver solver = new EJMLLinearSolver(5e-3, 1e-6);

		// Solves (one) linear equation, a x = b, for x[n]

		// Taken from https://en.wikipedia.org/wiki/Positive-definite_matrix
		double[][] a = new double[][] { 
			new double[] { 2, -1, 0 }, 
			new double[] { -1, 2, -1 },
			new double[] { 0, -1, 2 } };
		double[] b = new double[] { 3, 3, 4 };

		// Expected solution
		double[] x = new double[] { 4.75, 6.5, 5.25 };
		double[][] a_inv = new double[][] { 
			new double[] { 0.75, 0.5, 0.25 }, 
			new double[] { 0.5, 1, 0.5 },
			new double[] { 0.25, 0.5, 0.75 } };

		boolean result = solver.solve(a, b);
		solver.invertLastA(a);

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
		EJMLLinearSolver solver = new EJMLLinearSolver(5e-3, 1e-6);

		// Solves (one) linear equation, a x = b, for x[n]

		// Taken from https://en.wikipedia.org/wiki/Positive-definite_matrix
		double[][] a = new double[][] { 
			new double[] { 2, -1, 0 }, 
			new double[] { -1, 2, -1 },
			new double[] { 0, -1, 2 } };
		double[] b = new double[] { 3, 0, 4 };

		// Expected solution
		double[] x = new double[] { 3.25, 3.5, 3.75 };
		double[][] a_inv = new double[][] { 
			new double[] { 0.75, 0.5, 0.25 }, 
			new double[] { 0.5, 1, 0.5 },
			new double[] { 0.25, 0.5, 0.75 } };

		boolean result = solver.solve(a, b);
		solver.invertLastA(a);

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
		EJMLLinearSolver solver = new EJMLLinearSolver(5e-3, 1e-6);

		// Solves (one) linear equation, a x = b, for x[n]

		double[][] a = new double[][] { 
			new double[] { 2, 0, -1, 0 }, 
			new double[] { 0, 0, 0, 0 },
			new double[] { -1, 0, 2, -1 }, 
			new double[] { 0, 0, -1, 2 } };
		double[] b = new double[] { 3, 0, 3, 4 };

		// Expected solution
		double[] x = new double[] { 4.75, 0, 6.5, 5.25 };
		double[][] a_inv = new double[][] { 
			new double[] { 0.75, 0, 0.5, 0.25 }, 
			new double[] { 0, 0, 0, 0 },
			new double[] { 0.5, 0, 1, 0.5 },
			new double[] { 0.25, 0, 0.5, 0.75 } };

		boolean result = solver.solve(a, b);
		solver.invertLastA(a);

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
		DoubleEquality eq = new DoubleEquality(5e-3, 1e-16);
		solver.setEqual(eq);

		// Solves (one) linear equation, a x = b, for x[n]

		double[][] a = new double[][] { 
			new double[] { 2, 0, -1, 0, 0, 0 }, 
			new double[] { 0, 0, 0, 0, 0, 0 },
			new double[] { -1, 0, 2, 0, 0, -1 }, 
			new double[] { 0, 0, 0, 0, 0, 0 },
			new double[] { 0, 0, 0, 0, 0, 0 }, 
			new double[] { 0, 0, -1, 0, 0, 2 } };
		double[] b = new double[] { 3, 0, 3, 0, 0, 4 };

		// Expected solution
		double[] x = new double[] { 4.75, 0, 6.5, 0, 0, 5.25 };
		double[][] a_inv = new double[][] { 
			new double[] { 0.75, 0, 0.5, 0, 0, 0.25 },
			new double[] { 0, 0, 0, 0, 0, 0 }, 
			new double[] { 0.5, 0, 1, 0, 0, 0.5 },
			new double[] { 0, 0, 0, 0, 0, 0 }, 
			new double[] { 0, 0, 0, 0, 0, 0 },
			new double[] { 0.25, 0, 0.5, 0, 0, 0.75 } };

		boolean result = solver.solve(a, b);
		solver.invertLastA(a);

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
	public void canInvert()
	{
		EJMLLinearSolver solver = EJMLLinearSolver.createForInversion(1e-2);

		// Solves (one) linear equation, a x = b, for x[n]

		// Taken from https://en.wikipedia.org/wiki/Positive-definite_matrix
		double[][] a = new double[][] { 
			new double[] { 2, -1, 0 }, 
			new double[] { -1, 2, -1 },
			new double[] { 0, -1, 2 } };

		// Expected solution
		double[][] a_inv = new double[][] { 
			new double[] { 0.75, 0.5, 0.25 }, 
			new double[] { 0.5, 1, 0.5 },
			new double[] { 0.25, 0.5, 0.75 } };

		boolean result = solver.invert(a);

		Assert.assertTrue("Failed to invert", result);

		for (int i = 0; i < a[0].length; i++)
		{
			log("a[%d] = %s\n", i, Arrays.toString(a[i]));
			Assert.assertArrayEquals("Bad inversion", a_inv[i], a[i], 1e-4f);
		}
	}	
	
	@Test
	public void canInvertWithZeros()
	{
		EJMLLinearSolver solver = EJMLLinearSolver.createForInversion(1e-2);

		// Solves (one) linear equation, a x = b, for x[n]

		// Taken from https://en.wikipedia.org/wiki/Positive-definite_matrix
		double[][] a = new double[][] { 
			new double[] { 2, 0, -1, 0, 0, 0 }, 
			new double[] { 0, 0, 0, 0, 0, 0 },
			new double[] { -1, 0, 2, 0, 0, -1 }, 
			new double[] { 0, 0, 0, 0, 0, 0 },
			new double[] { 0, 0, 0, 0, 0, 0 }, 
			new double[] { 0, 0, -1, 0, 0, 2 } };

		// Expected solution
		double[][] a_inv = new double[][] { 
			new double[] { 0.75, 0, 0.5, 0, 0, 0.25 },
			new double[] { 0, 0, 0, 0, 0, 0 }, 
			new double[] { 0.5, 0, 1, 0, 0, 0.5 },
			new double[] { 0, 0, 0, 0, 0, 0 }, 
			new double[] { 0, 0, 0, 0, 0, 0 },
			new double[] { 0.25, 0, 0.5, 0, 0, 0.75 } };

		boolean result = solver.invert(a);

		Assert.assertTrue("Failed to invert", result);

		for (int i = 0; i < a[0].length; i++)
		{
			log("a[%d] = %s\n", i, Arrays.toString(a[i]));
			Assert.assertArrayEquals("Bad inversion", a_inv[i], a[i], 1e-4f);
		}
	}	
	
	@Test
	public void canInvertDiagonal()
	{
		EJMLLinearSolver solver = EJMLLinearSolver.createForInversion(1e-2);

		// Solves (one) linear equation, a x = b, for x[n]

		// Taken from https://en.wikipedia.org/wiki/Positive-definite_matrix
		double[][] a = new double[][] { 
			new double[] { 2, -1, 0 }, 
			new double[] { -1, 2, -1 },
			new double[] { 0, -1, 2 } };

		// Expected solution
		double[] e = new double[] { 0.75, 1, 0.75 };

		double[] o = solver.invertDiagonal(a);

		Assert.assertNotNull("Failed to invert", o);

		log("a diagonal = %s\n", Arrays.toString(o));
		Assert.assertArrayEquals("Bad inversion", e, o, 1e-4);
	}		
	
	@Test
	public void canInvertDiagonalWithZeros()
	{
		EJMLLinearSolver solver = EJMLLinearSolver.createForInversion(1e-2);

		// Solves (one) linear equation, a x = b, for x[n]

		// Taken from https://en.wikipedia.org/wiki/Positive-definite_matrix
		double[][] a = new double[][] { 
			new double[] { 2, 0, -1, 0, 0, 0 }, 
			new double[] { 0, 0, 0, 0, 0, 0 },
			new double[] { -1, 0, 2, 0, 0, -1 }, 
			new double[] { 0, 0, 0, 0, 0, 0 },
			new double[] { 0, 0, 0, 0, 0, 0 }, 
			new double[] { 0, 0, -1, 0, 0, 2 } };

		// Expected solution
		double[] e = new double[] { 0.75, 0, 1, 0, 0, 0.75 };

		double[] o = solver.invertDiagonal(a);

		Assert.assertNotNull("Failed to invert", o);

		log("a diagonal = %s\n", Arrays.toString(o));
		Assert.assertArrayEquals("Bad inversion", e, o, 1e-4);
	}		
	//@formatter:on

	private abstract class SolverTimingTask extends BaseTimingTask
	{
		DenseMatrix64F[] a;
		DenseMatrix64F[] b;
		// No validation for a pure speed test
		EJMLLinearSolver solver = new EJMLLinearSolver();

		public SolverTimingTask(String name, DenseMatrix64F[] a, DenseMatrix64F[] b)
		{
			super(name + " " + a[0].numCols);
			// Clone the data
			this.a = a;
			this.b = b;
			// Check the solver gets a good answer
			solver.setEqual(new DoubleEquality(5e-3, 1e-6));
			Object data = getData(0);

			a = (DenseMatrix64F[]) ((Object[]) data)[0];
			b = (DenseMatrix64F[]) ((Object[]) data)[1];
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
			DenseMatrix64F[] a = new DenseMatrix64F[n];
			DenseMatrix64F[] b = new DenseMatrix64F[n];
			while (n-- > 0)
			{
				a[n] = this.a[n].copy();
				b[n] = this.b[n].copy();
			}
			return new Object[] { a, b };
		}

		public Object run(Object data)
		{
			DenseMatrix64F[] a = (DenseMatrix64F[]) ((Object[]) data)[0];
			DenseMatrix64F[] b = (DenseMatrix64F[]) ((Object[]) data)[1];
			for (int i = 0; i < a.length; i++)
			{
				solve(a[i], b[i]);
			}
			return null;
		}

		abstract boolean solve(DenseMatrix64F a, DenseMatrix64F b);
	}

	private class LinearSolverTimingTask extends SolverTimingTask
	{
		public LinearSolverTimingTask(DenseMatrix64F[] a, DenseMatrix64F[] b)
		{
			super("Linear Solver", a, b);
		}

		boolean solve(DenseMatrix64F a, DenseMatrix64F b)
		{
			return solver.solveLinear(a, b);
		}
	}

	private class CholeskySolverTimingTask extends SolverTimingTask
	{
		public CholeskySolverTimingTask(DenseMatrix64F[] a, DenseMatrix64F[] b)
		{
			super("Cholesky Solver", a, b);
		}

		boolean solve(DenseMatrix64F a, DenseMatrix64F b)
		{
			return solver.solveCholesky(a, b);
		}
	}

	private class CholeskyLDLTSolverTimingTask extends SolverTimingTask
	{
		public CholeskyLDLTSolverTimingTask(DenseMatrix64F[] a, DenseMatrix64F[] b)
		{
			super("CholeskyLDLT Solver", a, b);
		}

		boolean solve(DenseMatrix64F a, DenseMatrix64F b)
		{
			return solver.solveCholeskyLDLT(a, b);
		}
	}

	private class PseudoInverseSolverTimingTask extends SolverTimingTask
	{
		public PseudoInverseSolverTimingTask(DenseMatrix64F[] a, DenseMatrix64F[] b)
		{
			super("PseudoInverse Solver", a, b);
		}

		boolean solve(DenseMatrix64F a, DenseMatrix64F b)
		{
			return solver.solvePseudoInverse(a, b);
		}
	}

	private class DirectInversionSolverTimingTask extends SolverTimingTask
	{
		public DirectInversionSolverTimingTask(DenseMatrix64F[] a, DenseMatrix64F[] b)
		{
			super("DirectInversion Solver", a, b);
		}

		boolean solve(DenseMatrix64F a, DenseMatrix64F b)
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
		runSpeedTest(GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED);
	}

	@Test
	public void runSpeedTest2()
	{
		runSpeedTest(GaussianFunctionFactory.FIT_SIMPLE_NS_NB_FIXED);
	}

	private void runSpeedTest(int flags)
	{
		final Gaussian2DFunction f0 = GaussianFunctionFactory.create2D(1, 10, 10, flags, null);
		int n = f0.size();
		final double[] y = new double[n];
		final TurboList<DenseMatrix64F> aList = new TurboList<DenseMatrix64F>();
		final TurboList<DenseMatrix64F> bList = new TurboList<DenseMatrix64F>();
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
							aList.add(EJMLLinearSolver.toA(alpha));
							bList.add(EJMLLinearSolver.toB(beta));
						}

		DenseMatrix64F[] a = aList.toArray(new DenseMatrix64F[aList.size()]);
		DenseMatrix64F[] b = bList.toArray(new DenseMatrix64F[bList.size()]);
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
