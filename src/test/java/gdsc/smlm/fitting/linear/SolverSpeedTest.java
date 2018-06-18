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
package gdsc.smlm.fitting.linear;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

import gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.SingleFreeCircularGaussian2DFunction;
import gdsc.test.TestSettings;

public class SolverSpeedTest
{
	private static ArrayList<float[][]> Adata = new ArrayList<float[][]>();
	private static ArrayList<float[]> Bdata = new ArrayList<float[]>();
	private static RandomGenerator rand = TestSettings.getRandomGenerator();

	private static synchronized void ensureData(int size)
	{
		while (Adata.size() < size)
		{
			float[][] a = new float[6][6];
			float[] b = new float[6];
			if (createData(rand, a, b, false))
			{
				Adata.add(a);
				Bdata.add(b);
			}
		}
	}

	@Test
	public void solveLinearAndGaussJordanReturnSameSolutionAndInversionResult()
	{
		int ITER = 100; ensureData(ITER);
		ArrayList<double[][]> A = copyAdouble(SolverSpeedTest.Adata, ITER);
		ArrayList<double[]> B = copyBdouble(SolverSpeedTest.Bdata, ITER);
		ArrayList<double[][]> A2 = copyAdouble(SolverSpeedTest.Adata, ITER);
		ArrayList<double[]> B2 = copyBdouble(SolverSpeedTest.Bdata, ITER);

		GaussJordan solver = new GaussJordan();
		EJMLLinearSolver solver2 = new EJMLLinearSolver();

		for (int i = 0; i < A.size(); i++)
		{
			double[][] a = A.get(i);
			double[] b = B.get(i);
			double[][] a2 = A2.get(i);
			double[] b2 = B2.get(i);
			boolean r1 = solver.solve(a, b);
			boolean r2 = solver2.solveLinear(a2, b2);
			solver2.invertLastA(a2);
			Assert.assertSame("Different solve result @ " + i, r1, r2);
			for (int j = 0; j < b.length; j++)
			{
				Assert.assertEquals("Different b result @ " + i, b[j] / b2[j], 1.0, 1e-2);
				String msg2 = "Different a[" + j + "] result @ " + i;
				for (int k = 0; k < b.length; k++)
					Assert.assertEquals(msg2, a[j][k] / a2[j][k], 1, 0.2);
			}
		}
	}

	@Test
	public void solveLinearAndGaussJordanReturnSameSolutionResult()
	{
		int ITER = 100; ensureData(ITER);
		ArrayList<double[][]> A = copyAdouble(SolverSpeedTest.Adata, ITER);
		ArrayList<double[]> B = copyBdouble(SolverSpeedTest.Bdata, ITER);
		ArrayList<double[][]> A2 = copyAdouble(SolverSpeedTest.Adata, ITER);
		ArrayList<double[]> B2 = copyBdouble(SolverSpeedTest.Bdata, ITER);

		GaussJordan solver = new GaussJordan();
		EJMLLinearSolver solver2 = new EJMLLinearSolver();

		for (int i = 0; i < ITER; i++)
		{
			double[][] a = A.get(i);
			double[] b = B.get(i);
			double[][] a2 = A2.get(i);
			double[] b2 = B2.get(i);
			boolean r1 = solver.solve(a, b);
			boolean r2 = solver2.solve(a2, b2);
			Assert.assertSame("Different solve result @ " + i, r1, r2);
			for (int j = 0; j < b.length; j++)
				Assert.assertEquals("Different b result @ " + i, b[j] / b2[j], 1.0, 1e-2);
		}
	}

	@Test
	public void gaussJordanFloatAndDoubleReturnSameSolutionAndInversionResult()
	{
		int ITER = 100; ensureData(ITER);
		ArrayList<float[][]> A = copyAfloat(SolverSpeedTest.Adata, ITER);
		ArrayList<float[]> B = copyBfloat(SolverSpeedTest.Bdata, ITER);
		ArrayList<double[][]> A2 = copyAdouble(SolverSpeedTest.Adata, ITER);
		ArrayList<double[]> B2 = copyBdouble(SolverSpeedTest.Bdata, ITER);

		GaussJordan solver = new GaussJordan();

		for (int i = 0; i < A.size(); i++)
		{
			float[][] a = A.get(i);
			float[] b = B.get(i);
			double[][] a2 = A2.get(i);
			double[] b2 = B2.get(i);
			boolean r1 = solver.solve(a, b);
			boolean r2 = solver.solve(a2, b2);
			Assert.assertSame("Different solve result @ " + i, r1, r2);
			for (int j = 0; j < b.length; j++)
			{
				Assert.assertEquals("Different b result @ " + i, b[j] / b2[j], 1.0, 1e-2);
				String msg2 = "Different a[" + j + "] result @ " + i;
				for (int k = 0; k < b.length; k++)
					Assert.assertEquals(msg2, a[j][k] / a2[j][k], 1, 0.2);
			}
		}
	}

	@Test
	public void solveLinearWithInversionIsFasterThanGaussJordanFloat()
	{
		TestSettings.assumeMediumComplexity();

		int ITER = 10000; ensureData(ITER);
		ArrayList<float[][]> A = copyAfloat(SolverSpeedTest.Adata, ITER);
		ArrayList<float[]> B = copyBfloat(SolverSpeedTest.Bdata, ITER);
		ArrayList<double[][]> A2 = copyAdouble(SolverSpeedTest.Adata, ITER);
		ArrayList<double[]> B2 = copyBdouble(SolverSpeedTest.Bdata, ITER);

		GaussJordan solver = new GaussJordan();
		EJMLLinearSolver solver2 = new EJMLLinearSolver();

		runFloat(copyAfloat(SolverSpeedTest.Adata, ITER), copyBfloat(SolverSpeedTest.Bdata, ITER), ITER, solver);
		solveLinearWithInversion(copyAdouble(SolverSpeedTest.Adata, ITER), copyBdouble(SolverSpeedTest.Bdata, ITER), ITER, solver2);

		long start1 = System.nanoTime();
		runFloat(A, B, ITER, solver);
		start1 = System.nanoTime() - start1;

		long start2 = System.nanoTime();
		solveLinearWithInversion(A2, B2, ITER, solver2);
		start2 = System.nanoTime() - start2;

		log("GaussJordanFloat = %d : LinearSolver.solveLinearWithInversion = %d : %fx\n", start1, start2,
				(1.0 * start1) / start2);
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(start2 < start1);
	}

	@Test
	public void solveLinearIsFasterThanGaussJordanFloat()
	{
		TestSettings.assumeMediumComplexity();

		int ITER = 10000; ensureData(ITER);
		ArrayList<float[][]> A = copyAfloat(SolverSpeedTest.Adata, ITER);
		ArrayList<float[]> B = copyBfloat(SolverSpeedTest.Bdata, ITER);
		ArrayList<double[][]> A2 = copyAdouble(SolverSpeedTest.Adata, ITER);
		ArrayList<double[]> B2 = copyBdouble(SolverSpeedTest.Bdata, ITER);

		GaussJordan solver = new GaussJordan();
		EJMLLinearSolver solver2 = new EJMLLinearSolver();

		runFloat(copyAfloat(SolverSpeedTest.Adata, ITER), copyBfloat(SolverSpeedTest.Bdata, ITER), ITER, solver);
		solveLinear(copyAdouble(SolverSpeedTest.Adata, ITER), copyBdouble(SolverSpeedTest.Bdata, ITER), ITER, solver2);

		long start1 = System.nanoTime();
		runFloat(A, B, ITER, solver);
		start1 = System.nanoTime() - start1;

		long start2 = System.nanoTime();
		solveLinear(A2, B2, ITER, solver2);
		start2 = System.nanoTime() - start2;

		log("GaussJordanFloat = %d : LinearSolver.solveLinear = %d : %fx\n", start1, start2, (1.0 * start1) / start2);
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(start2 < start1);
	}

	protected void runFloat(ArrayList<float[][]> A, ArrayList<float[]> B, int ITER, GaussJordan solver)
	{
		for (int i = 0; i < ITER; i++)
		{
			solver.solve(A.get(i), B.get(i));
		}
	}

	@Test
	public void solveLinearWithInversionIsFasterThanGaussJordanDouble()
	{
		TestSettings.assumeMediumComplexity();

		int ITER = 10000; ensureData(ITER);
		ArrayList<double[][]> A = copyAdouble(SolverSpeedTest.Adata, ITER);
		ArrayList<double[]> B = copyBdouble(SolverSpeedTest.Bdata, ITER);
		ArrayList<double[][]> A2 = copyAdouble(SolverSpeedTest.Adata, ITER);
		ArrayList<double[]> B2 = copyBdouble(SolverSpeedTest.Bdata, ITER);

		GaussJordan solver = new GaussJordan();
		EJMLLinearSolver solver2 = new EJMLLinearSolver();

		solveGaussJordan(copyAdouble(SolverSpeedTest.Adata, ITER), copyBdouble(SolverSpeedTest.Bdata, ITER), ITER, solver);
		solveLinearWithInversion(copyAdouble(SolverSpeedTest.Adata, ITER), copyBdouble(SolverSpeedTest.Bdata, ITER), ITER, solver2);

		long start1 = System.nanoTime();
		solveGaussJordan(A, B, ITER, solver);
		start1 = System.nanoTime() - start1;

		long start2 = System.nanoTime();
		solveLinearWithInversion(A2, B2, ITER, solver2);
		start2 = System.nanoTime() - start2;

		log("GaussJordanDouble = %d : LinearSolver.solveLinearWithInversion = %d : %fx\n", start1, start2,
				(1.0 * start1) / start2);
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(start2 < start1);
	}

	@Test
	public void solveLinearIsFasterThanGaussJordanDouble()
	{
		TestSettings.assumeMediumComplexity();

		int ITER = 10000; ensureData(ITER);
		ArrayList<double[][]> A = copyAdouble(SolverSpeedTest.Adata, ITER);
		ArrayList<double[]> B = copyBdouble(SolverSpeedTest.Bdata, ITER);
		ArrayList<double[][]> A2 = copyAdouble(SolverSpeedTest.Adata, ITER);
		ArrayList<double[]> B2 = copyBdouble(SolverSpeedTest.Bdata, ITER);

		GaussJordan solver = new GaussJordan();
		EJMLLinearSolver solver2 = new EJMLLinearSolver();

		solveGaussJordan(copyAdouble(SolverSpeedTest.Adata, ITER), copyBdouble(SolverSpeedTest.Bdata, ITER), ITER, solver);
		solveLinear(copyAdouble(SolverSpeedTest.Adata, ITER), copyBdouble(SolverSpeedTest.Bdata, ITER), ITER, solver2);

		long start1 = System.nanoTime();
		solveGaussJordan(A, B, ITER, solver);
		start1 = System.nanoTime() - start1;

		long start2 = System.nanoTime();
		solveLinear(A2, B2, ITER, solver2);
		start2 = System.nanoTime() - start2;

		log("GaussJordanDouble = %d : LinearSolver.solveLinear = %d : %fx\n", start1, start2, (1.0 * start1) / start2);
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(start2 < start1);
	}

	@Test
	public void solveCholeskyIsFasterThanGaussJordanDouble()
	{
		TestSettings.assumeMediumComplexity();

		int ITER = 10000; ensureData(ITER);
		ArrayList<double[][]> A = copyAdouble(SolverSpeedTest.Adata, ITER);
		ArrayList<double[]> B = copyBdouble(SolverSpeedTest.Bdata, ITER);
		ArrayList<double[][]> A2 = copyAdouble(SolverSpeedTest.Adata, ITER);
		ArrayList<double[]> B2 = copyBdouble(SolverSpeedTest.Bdata, ITER);

		GaussJordan solver = new GaussJordan();
		EJMLLinearSolver solver2 = new EJMLLinearSolver();

		solveGaussJordan(copyAdouble(SolverSpeedTest.Adata, ITER), copyBdouble(SolverSpeedTest.Bdata, ITER), ITER, solver);
		solveCholesky(copyAdouble(SolverSpeedTest.Adata, ITER), copyBdouble(SolverSpeedTest.Bdata, ITER), ITER, solver2);

		long start1 = System.nanoTime();
		solveGaussJordan(A, B, ITER, solver);
		start1 = System.nanoTime() - start1;

		long start2 = System.nanoTime();
		solveCholesky(A2, B2, ITER, solver2);
		start2 = System.nanoTime() - start2;

		log("GaussJordanDouble = %d : LinearSolver.solveCholesky = %d : %fx\n", start1, start2,
				(1.0 * start1) / start2);
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(start2 < start1);
	}

	@Test
	public void solveCholeskyLDLTIsFasterThanGaussJordanDouble()
	{
		TestSettings.assumeMediumComplexity();

		int ITER = 10000; ensureData(ITER);
		ArrayList<double[][]> A = copyAdouble(SolverSpeedTest.Adata, ITER);
		ArrayList<double[]> B = copyBdouble(SolverSpeedTest.Bdata, ITER);
		ArrayList<double[][]> A2 = copyAdouble(SolverSpeedTest.Adata, ITER);
		ArrayList<double[]> B2 = copyBdouble(SolverSpeedTest.Bdata, ITER);

		GaussJordan solver = new GaussJordan();
		EJMLLinearSolver solver2 = new EJMLLinearSolver();

		solveGaussJordan(copyAdouble(SolverSpeedTest.Adata, ITER), copyBdouble(SolverSpeedTest.Bdata, ITER), ITER, solver);
		solveCholeskyLDLT(copyAdouble(SolverSpeedTest.Adata, ITER), copyBdouble(SolverSpeedTest.Bdata, ITER), ITER, solver2);

		long start1 = System.nanoTime();
		solveGaussJordan(A, B, ITER, solver);
		start1 = System.nanoTime() - start1;

		long start2 = System.nanoTime();
		solveCholeskyLDLT(A2, B2, ITER, solver2);
		start2 = System.nanoTime() - start2;

		log("GaussJordanDouble = %d : LinearSolver.solveCholeskyLDLT = %d : %fx\n", start1, start2,
				(1.0 * start1) / start2);
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(start2 < start1);
	}

	@Test
	public void solveIsFasterThanGaussJordanDouble()
	{
		TestSettings.assumeMediumComplexity();

		int ITER = 10000; ensureData(ITER);
		ArrayList<double[][]> A = copyAdouble(SolverSpeedTest.Adata, ITER);
		ArrayList<double[]> B = copyBdouble(SolverSpeedTest.Bdata, ITER);
		ArrayList<double[][]> A2 = copyAdouble(SolverSpeedTest.Adata, ITER);
		ArrayList<double[]> B2 = copyBdouble(SolverSpeedTest.Bdata, ITER);

		GaussJordan solver = new GaussJordan();
		EJMLLinearSolver solver2 = new EJMLLinearSolver();

		solveGaussJordan(copyAdouble(SolverSpeedTest.Adata, ITER), copyBdouble(SolverSpeedTest.Bdata, ITER), ITER, solver);
		solve(copyAdouble(SolverSpeedTest.Adata, ITER), copyBdouble(SolverSpeedTest.Bdata, ITER), ITER, solver2);

		long start1 = System.nanoTime();
		solveGaussJordan(A, B, ITER, solver);
		start1 = System.nanoTime() - start1;

		long start2 = System.nanoTime();
		solve(A2, B2, ITER, solver2);
		start2 = System.nanoTime() - start2;

		log("GaussJordanDouble = %d : LinearSolver.solve = %d : %fx\n", start1, start2, (1.0 * start1) / start2);
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(start2 < start1);
	}

	private static boolean createData(RandomGenerator rand, float[][] alpha, float[] beta, boolean positiveDifinite)
	{
		// Generate a 2D Gaussian
		SingleFreeCircularGaussian2DFunction func = new SingleFreeCircularGaussian2DFunction(10, 10);
		double[] a = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
		a[Gaussian2DFunction.BACKGROUND] = 2 + rand.nextDouble() * 2;
		a[Gaussian2DFunction.SIGNAL] = 100 + rand.nextDouble() * 5;
		a[Gaussian2DFunction.X_POSITION] = 4.5 + rand.nextDouble();
		a[Gaussian2DFunction.Y_POSITION] = 4.5 + rand.nextDouble();
		a[Gaussian2DFunction.X_SD] = 1 + rand.nextDouble();
		a[Gaussian2DFunction.Y_SD] = 1 + rand.nextDouble();
		a[Gaussian2DFunction.ANGLE] = rand.nextDouble();

		int[] x = new int[100];
		double[] y = new double[100];
		func.initialise(a);
		for (int i = 0; i < x.length; i++)
		{
			// Add random noise
			y[i] = func.eval(i) + ((rand.nextDouble() < 0.5) ? -rand.nextDouble() * 5 : rand.nextDouble() * 5);
		}

		// Randomise parameters
		for (int i = 0; i < a.length; i++)
			a[i] += (rand.nextDouble() < 0.5) ? -rand.nextDouble() : rand.nextDouble();

		// Compute the Hessian and parameter gradient vector
		GradientCalculator calc = new GradientCalculator(6);
		double[][] alpha2 = new double[6][6];
		double[] beta2 = new double[6];
		calc.findLinearised(y.length, y, a, alpha2, beta2, func);

		// Update the Hessian using a lambda shift
		double lambda = 1.001;
		for (int i = 0; i < alpha2.length; i++)
			alpha2[i][i] *= lambda;

		// Copy back
		for (int i = 0; i < beta.length; i++)
		{
			beta[i] = (float) beta2[i];
			for (int j = 0; j < beta.length; j++)
			{
				alpha[i][j] = (float) alpha2[i][j];
			}
		}

		// Check for a positive definite matrix
		if (positiveDifinite)
		{
			EJMLLinearSolver solver = new EJMLLinearSolver();
			return solver.solveCholeskyLDLT(copydouble(alpha), copydouble(beta));
		}

		return true;
	}

	private static ArrayList<float[][]> copyAfloat(ArrayList<float[][]> a, int iter)
	{
		iter = FastMath.min(a.size(), iter);
		ArrayList<float[][]> a2 = new ArrayList<float[][]>(iter);
		for (int i = 0; i < iter; i++)
			a2.add(copyfloat(a.get(i)));
		return a2;
	}

	private static float[][] copyfloat(float[][] d)
	{
		float[][] d2 = new float[d.length][d.length];
		for (int i = 0; i < d.length; i++)
			for (int j = 0; j < d.length; j++)
				d2[i][j] = d[i][j];
		return d2;
	}

	private static ArrayList<float[]> copyBfloat(ArrayList<float[]> b, int iter)
	{
		iter = FastMath.min(b.size(), iter);
		ArrayList<float[]> b2 = new ArrayList<float[]>(iter);
		for (int i = 0; i < iter; i++)
			b2.add(Arrays.copyOf(b.get(i), b.get(i).length));
		return b2;
	}

	private static ArrayList<double[][]> copyAdouble(ArrayList<float[][]> a, int iter)
	{
		iter = FastMath.min(a.size(), iter);
		ArrayList<double[][]> a2 = new ArrayList<double[][]>(iter);
		for (int i = 0; i < iter; i++)
			a2.add(copydouble(a.get(i)));
		return a2;
	}

	private static double[][] copydouble(float[][] d)
	{
		double[][] d2 = new double[d.length][d.length];
		for (int i = 0; i < d.length; i++)
			for (int j = 0; j < d.length; j++)
				d2[i][j] = d[i][j];
		return d2;
	}

	private static ArrayList<double[]> copyBdouble(ArrayList<float[]> b, int iter)
	{
		iter = FastMath.min(b.size(), iter);
		ArrayList<double[]> b2 = new ArrayList<double[]>(iter);
		for (int i = 0; i < iter; i++)
			b2.add(copydouble(b.get(i)));
		return b2;
	}

	private static double[] copydouble(float[] d)
	{
		double[] d2 = new double[d.length];
		for (int i = 0; i < d.length; i++)
			d2[i] = d[i];
		return d2;
	}

	protected void solveGaussJordan(ArrayList<double[][]> A, ArrayList<double[]> B, int ITER, GaussJordan solver)
	{
		ITER = FastMath.min(ITER, A.size());
		for (int i = 0; i < ITER; i++)
		{
			solver.solve(A.get(i), B.get(i));
		}
	}

	protected void solveLinearWithInversion(ArrayList<double[][]> A, ArrayList<double[]> B, int ITER,
			EJMLLinearSolver solver)
	{
		ITER = FastMath.min(ITER, A.size());
		for (int i = 0; i < ITER; i++)
		{
			double[][] a = A.get(i);
			solver.solveLinear(a, B.get(i));
			solver.invertLastA(a);
		}
	}

	protected void solveLinear(ArrayList<double[][]> A, ArrayList<double[]> B, int ITER, EJMLLinearSolver solver)
	{
		ITER = FastMath.min(ITER, A.size());
		for (int i = 0; i < ITER; i++)
		{
			solver.solveLinear(A.get(i), B.get(i));
		}
	}

	protected void solveCholesky(ArrayList<double[][]> A, ArrayList<double[]> B, int ITER, EJMLLinearSolver solver)
	{
		ITER = FastMath.min(ITER, A.size());
		for (int i = 0; i < ITER; i++)
		{
			solver.solveCholesky(A.get(i), B.get(i));
		}
	}

	protected void solveCholeskyLDLT(ArrayList<double[][]> A, ArrayList<double[]> B, int ITER, EJMLLinearSolver solver)
	{
		ITER = FastMath.min(ITER, A.size());
		for (int i = 0; i < ITER; i++)
		{
			solver.solveCholeskyLDLT(A.get(i), B.get(i));
		}
	}

	protected void solve(ArrayList<double[][]> A, ArrayList<double[]> B, int ITER, EJMLLinearSolver solver)
	{
		ITER = FastMath.min(ITER, A.size());
		for (int i = 0; i < ITER; i++)
		{
			solver.solve(A.get(i), B.get(i));
		}
	}

	void log(String format, Object... args)
	{
		System.out.printf(format, args);
	}
}
