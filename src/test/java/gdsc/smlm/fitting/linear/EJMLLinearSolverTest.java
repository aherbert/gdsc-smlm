package gdsc.smlm.fitting.linear;

import java.util.Arrays;

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import gdsc.smlm.fitting.linear.EJMLLinearSolver;

import org.junit.Assert;
import org.junit.Test;

public class EJMLLinearSolverTest
{
	@Test
	public void floatCanSolveLinearEquation()
	{
		EJMLLinearSolver solver = new EJMLLinearSolver();
		DoubleEquality eq = new DoubleEquality(3, 1e-16);
		solver.setEqual(eq);
		
		// Solves (one) linear equation, a x = b, for x[n]
		
		// Example taken from http://eigen.tuxfamily.org/dox/TutorialLinearAlgebra.html
		float[][] a = new float[][] { 
				new float[] { 1, 2, 3 }, 
				new float[] { 4, 5, 6 }, 
				new float[] { 7, 8, 10 } 
				};
		float[] b = new float[] { 3, 3, 4 };
		
		// Expected solution
		float[] x = new float[] { -2, 1, 1 };
		float[][] a_inv = new float[][] { 
				new float[] { -0.6666667f, -1.3333334f, 1.0f }, 
				new float[] { -0.6666667f, 3.6666667f, -2.0f }, 
				new float[] { 1.0f, -2.0f, 1.0f } 
				};
		
		boolean result = solver.solve(a, b);
		solver.invert(a);
		
		Assert.assertTrue("Failed to solve", result);
		Assert.assertArrayEquals("Bad solution", x, b, 1e-4f);
		
		log("x = %s\n", Arrays.toString(b));
		for (int i=0; i<b.length; i++)
		{
			Assert.assertArrayEquals("Bad inversion", a_inv[i], a[i], 1e-4f);
			log("a[%d] = %s\n", i, Arrays.toString(a[i]));
		}
	}
	
	@Test
	public void floatCanSolveLinearEquationWithZero()
	{
		EJMLLinearSolver solver = new EJMLLinearSolver();
		DoubleEquality eq = new DoubleEquality(3, 1e-16);
		solver.setEqual(eq);
		
		// Solves (one) linear equation, a x = b, for x[n]
		
		// Example taken from http://eigen.tuxfamily.org/dox/TutorialLinearAlgebra.html
		float[][] a = new float[][] { 
				new float[] { 1, 0, 2, 3 }, 
				new float[] { 0, 0, 0, 0 }, 
				new float[] { 4, 0, 5, 6 }, 
				new float[] { 7, 0, 8, 10 } 
				};
		float[] b = new float[] { 3, 0, 3, 4 };
		
		// Expected solution
		float[] x = new float[] { -2, 0, 1, 1 };
		float[][] a_inv = new float[][] { 
				new float[] { -0.6666667f, 0, -1.3333334f, 1.0f }, 
				new float[] { 0, 0, 0, 0 }, 
				new float[] { -0.6666667f, 0, 3.6666667f, -2.0f }, 
				new float[] { 1.0f, 0, -2.0f, 1.0f } 
				};
		
		boolean result = solver.solveWithZeros(a, b);
		solver.invert(a);
		
		Assert.assertTrue("Failed to solve", result);
		Assert.assertArrayEquals("Bad solution", x, b, 1e-4f);
		
		log("x = %s\n", Arrays.toString(b));
		for (int i=0; i<b.length; i++)
		{
			Assert.assertArrayEquals("Bad inversion", a_inv[i], a[i], 1e-4f);
			log("a[%d] = %s\n", i, Arrays.toString(a[i]));
		}
	}
	
	@Test
	public void floatCanSolveLinearEquationWithZeros()
	{
		EJMLLinearSolver solver = new EJMLLinearSolver();
		DoubleEquality eq = new DoubleEquality(3, 1e-16);
		solver.setEqual(eq);
		
		// Solves (one) linear equation, a x = b, for x[n]
		
		// Example taken from http://eigen.tuxfamily.org/dox/TutorialLinearAlgebra.html
		float[][] a = new float[][] { 
				new float[] { 1, 0, 2, 0, 0, 3 }, 
				new float[] { 0, 0, 0, 0, 0, 0 }, 
				new float[] { 4, 0, 5, 0, 0, 6 }, 
				new float[] { 0, 0, 0, 0, 0, 0 }, 
				new float[] { 0, 0, 0, 0, 0, 0 }, 
				new float[] { 7, 0, 8, 0, 0, 10 } 
				};
		float[] b = new float[] { 3, 0, 3, 0, 0, 4 };
		
		// Expected solution
		float[] x = new float[] { -2, 0, 1, 0, 0, 1 };
		float[][] a_inv = new float[][] { 
				new float[] { -0.6666667f, 0, -1.3333334f, 0, 0, 1.0f }, 
				new float[] { 0, 0, 0, 0, 0, 0 }, 
				new float[] { -0.6666667f, 0, 3.6666667f, 0, 0, -2.0f }, 
				new float[] { 0, 0, 0, 0, 0, 0 }, 
				new float[] { 0, 0, 0, 0, 0, 0 }, 
				new float[] { 1.0f, 0, -2.0f, 0, 0, 1.0f } 
				};
		
		boolean result = solver.solveWithZeros(a, b);
		solver.invert(a);
		
		Assert.assertTrue("Failed to solve", result);
		Assert.assertArrayEquals("Bad solution", x, b, 1e-4f);
		
		log("x = %s\n", Arrays.toString(b));
		for (int i=0; i<b.length; i++)
		{
			Assert.assertArrayEquals("Bad inversion", a_inv[i], a[i], 1e-4f);
			log("a[%d] = %s\n", i, Arrays.toString(a[i]));
		}
	}
	
	@Test
	public void doubleCanSolveLinearEquation()
	{
		EJMLLinearSolver solver = new EJMLLinearSolver();
		DoubleEquality eq = new DoubleEquality(DoubleEquality.getRelativeErrorTerm(3), 1e-16);
		solver.setEqual(eq);
		
		// Solves (one) linear equation, a x = b, for x[n]
		
		// Example taken from http://eigen.tuxfamily.org/dox/TutorialLinearAlgebra.html
		double[][] a = new double[][] { 
				new double[] { 1, 2, 3 }, 
				new double[] { 4, 5, 6 }, 
				new double[] { 7, 8, 10 } 
				};
		double[] b = new double[] { 3, 3, 4 };
		
		// Expected solution
		double[] x = new double[] { -2, 1, 1 };
		double[][] a_inv = new double[][] { 
				new double[] { -0.6666667f, -1.3333334f, 1.0f }, 
				new double[] { -0.6666667f, 3.6666667f, -2.0f }, 
				new double[] { 1.0f, -2.0f, 1.0f } 
				};
		
		boolean result = solver.solve(a, b);
		solver.invert(a);
		
		Assert.assertTrue("Failed to solve", result);
		Assert.assertArrayEquals("Bad solution", x, b, 1e-4f);
		
		log("x = %s\n", Arrays.toString(b));
		for (int i=0; i<b.length; i++)
		{
			Assert.assertArrayEquals("Bad inversion", a_inv[i], a[i], 1e-4f);
			log("a[%d] = %s\n", i, Arrays.toString(a[i]));
		}
	}
	
	@Test
	public void doubleCanSolveLinearEquationWithZero()
	{
		EJMLLinearSolver solver = new EJMLLinearSolver();
		DoubleEquality eq = new DoubleEquality(DoubleEquality.getRelativeErrorTerm(3), 1e-16);
		solver.setEqual(eq);
		
		// Solves (one) linear equation, a x = b, for x[n]
		
		// Example taken from http://eigen.tuxfamily.org/dox/TutorialLinearAlgebra.html
		double[][] a = new double[][] { 
				new double[] { 1, 0, 2, 3 }, 
				new double[] { 0, 0, 0, 0 }, 
				new double[] { 4, 0, 5, 6 }, 
				new double[] { 7, 0, 8, 10 } 
				};
		double[] b = new double[] { 3, 0, 3, 4 };
		
		// Expected solution
		double[] x = new double[] { -2, 0, 1, 1 };
		double[][] a_inv = new double[][] { 
				new double[] { -0.6666667f, 0, -1.3333334f, 1.0f }, 
				new double[] { 0, 0, 0, 0 }, 
				new double[] { -0.6666667f, 0, 3.6666667f, -2.0f }, 
				new double[] { 1.0f, 0, -2.0f, 1.0f } 
				};
		
		boolean result = solver.solveWithZeros(a, b);
		solver.invert(a);
		
		Assert.assertTrue("Failed to solve", result);
		Assert.assertArrayEquals("Bad solution", x, b, 1e-4f);
		
		log("x = %s\n", Arrays.toString(b));
		for (int i=0; i<b.length; i++)
		{
			Assert.assertArrayEquals("Bad inversion", a_inv[i], a[i], 1e-4f);
			log("a[%d] = %s\n", i, Arrays.toString(a[i]));
		}
	}
	
	@Test
	public void doubleCanSolveLinearEquationWithZeros()
	{
		EJMLLinearSolver solver = new EJMLLinearSolver();
		DoubleEquality eq = new DoubleEquality(3, 1e-16);
		solver.setEqual(eq);
		
		// Solves (one) linear equation, a x = b, for x[n]
		
		// Example taken from http://eigen.tuxfamily.org/dox/TutorialLinearAlgebra.html
		double[][] a = new double[][] { 
				new double[] { 1, 0, 2, 0, 0, 3 }, 
				new double[] { 0, 0, 0, 0, 0, 0 }, 
				new double[] { 4, 0, 5, 0, 0, 6 }, 
				new double[] { 0, 0, 0, 0, 0, 0 }, 
				new double[] { 0, 0, 0, 0, 0, 0 }, 
				new double[] { 7, 0, 8, 0, 0, 10 } 
				};
		double[] b = new double[] { 3, 0, 3, 0, 0, 4 };
		
		// Expected solution
		double[] x = new double[] { -2, 0, 1, 0, 0, 1 };
		double[][] a_inv = new double[][] { 
				new double[] { -0.6666667f, 0, -1.3333334f, 0, 0, 1.0f }, 
				new double[] { 0, 0, 0, 0, 0, 0 }, 
				new double[] { -0.6666667f, 0, 3.6666667f, 0, 0, -2.0f }, 
				new double[] { 0, 0, 0, 0, 0, 0 }, 
				new double[] { 0, 0, 0, 0, 0, 0 }, 
				new double[] { 1.0f, 0, -2.0f, 0, 0, 1.0f } 
				};
		
		boolean result = solver.solveWithZeros(a, b);
		solver.invert(a);
		
		Assert.assertTrue("Failed to solve", result);
		Assert.assertArrayEquals("Bad solution", x, b, 1e-4f);
		
		log("x = %s\n", Arrays.toString(b));
		for (int i=0; i<b.length; i++)
		{
			Assert.assertArrayEquals("Bad inversion", a_inv[i], a[i], 1e-4f);
			log("a[%d] = %s\n", i, Arrays.toString(a[i]));
		}
	}
	
	void log(String format, Object... args)
	{
		System.out.printf(format, args);
	}
}
