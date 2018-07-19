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
package uk.ac.sussex.gdsc.smlm.fitting.linear;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.FastMath;
import org.ejml.data.DenseMatrix64F;
import org.junit.Test;

import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.SingleFreeCircularGaussian2DFunction;
import uk.ac.sussex.gdsc.test.TestCounter;
import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.junit4.TestAssert;
import uk.ac.sussex.gdsc.test.junit4.TestAssume;

@SuppressWarnings({ "javadoc" })
public class SolverSpeedTest
{
	private static ArrayList<float[][]> Adata = new ArrayList<>();
	private static ArrayList<float[]> Bdata = new ArrayList<>();
	private static RandomGenerator rand = TestSettings.getRandomGenerator();

	private static synchronized void ensureData(int size)
	{
		while (Adata.size() < size)
		{
			final float[][] a = new float[6][6];
			final float[] b = new float[6];
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
		final int ITER = 100;
		ensureData(ITER);
		final ArrayList<double[][]> A = copyAdouble(SolverSpeedTest.Adata, ITER);
		final ArrayList<double[]> B = copyBdouble(SolverSpeedTest.Bdata, ITER);
		final ArrayList<double[][]> A2 = copyAdouble(SolverSpeedTest.Adata, ITER);
		final ArrayList<double[]> B2 = copyBdouble(SolverSpeedTest.Bdata, ITER);

		final GaussJordan solver = new GaussJordan();
		final EJMLLinearSolver solver2 = new EJMLLinearSolver();

		final int failureLimit = TestCounter.computeFailureLimit(ITER, 0.1);
		final TestCounter failCounter = new TestCounter(failureLimit, 2);

		int c = 0;
		for (int i = 0; i < A.size(); i++)
		{
			final double[][] a = A.get(i);
			final double[] b = B.get(i);
			final double[][] a2 = A2.get(i);
			final double[] b2 = B2.get(i);
			final boolean r1 = solver.solve(a, b);
			final boolean r2 = solver2.solveLinear(a2, b2);
			solver2.invertLastA(a2);
			//Assert.assertTrue("Different solve result @ " + i, r1 == r2);
			if (r1 && r2)
			{
				failCounter.run(0, () -> {
					TestAssert.assertArrayEqualsRelative("Different b result", b, b2, 1e-2);
				});
				failCounter.run(0, () -> {
					TestAssert.assertDoubleArrayEqualsRelative("Different a result", a, a2, 1e-2);
				});
			}
			else
				c++;
		}
		if (c > ITER / 2)
			TestAssert.fail("Failed to solve %d / %d", c, ITER);
	}

	@Test
	public void solveLinearAndGaussJordanReturnSameSolutionResult()
	{
		final int ITER = 100;
		ensureData(ITER);
		final ArrayList<double[][]> A = copyAdouble(SolverSpeedTest.Adata, ITER);
		final ArrayList<double[]> B = copyBdouble(SolverSpeedTest.Bdata, ITER);
		final ArrayList<double[][]> A2 = copyAdouble(SolverSpeedTest.Adata, ITER);
		final ArrayList<double[]> B2 = copyBdouble(SolverSpeedTest.Bdata, ITER);

		final GaussJordan solver = new GaussJordan();
		final EJMLLinearSolver solver2 = new EJMLLinearSolver();

		final int failureLimit = TestCounter.computeFailureLimit(ITER, 0.1);
		final TestCounter failCounter = new TestCounter(failureLimit);

		int c = 0;
		for (int i = 0; i < ITER; i++)
		{
			final double[][] a = A.get(i);
			final double[] b = B.get(i);
			final double[][] a2 = A2.get(i);
			final double[] b2 = B2.get(i);
			final boolean r1 = solver.solve(a, b);
			final boolean r2 = solver2.solve(a2, b2);
			//Assert.assertTrue("Different solve result @ " + i, r1 == r2);
			if (r1 && r2)
				failCounter.run(() -> {
					TestAssert.assertArrayEqualsRelative("Different b result", b, b2, 1e-2);
				});
			else
				c++;
		}
		if (c > ITER / 2)
			TestAssert.fail("Failed to solve %d / %d", c, ITER);
	}

	@Test
	public void gaussJordanFloatAndDoubleReturnSameSolutionAndInversionResult()
	{
		final int ITER = 100;
		ensureData(ITER);
		final ArrayList<float[][]> A = copyAfloat(SolverSpeedTest.Adata, ITER);
		final ArrayList<float[]> B = copyBfloat(SolverSpeedTest.Bdata, ITER);
		final ArrayList<double[][]> A2 = copyAdouble(SolverSpeedTest.Adata, ITER);
		final ArrayList<double[]> B2 = copyBdouble(SolverSpeedTest.Bdata, ITER);

		final GaussJordan solver = new GaussJordan();

		final int failureLimit = TestCounter.computeFailureLimit(ITER, 0.1);
		final TestCounter failCounter = new TestCounter(failureLimit, 2);

		int c = 0;
		for (int i = 0; i < A.size(); i++)
		{
			final float[][] a = A.get(i);
			final float[] b = B.get(i);
			final double[][] a2 = A2.get(i);
			final double[] b2 = B2.get(i);
			final boolean r1 = solver.solve(a, b);
			final boolean r2 = solver.solve(a2, b2);
			//Assert.assertTrue("Different solve result @ " + i, r1 == r2);
			if (r1 && r2)
			{
				final double[] b1 = SimpleArrayUtils.toDouble(b);
				final double[][] a1 = new double[a.length][];
				for (int j = a1.length; j-- > 0;)
					a1[j] = SimpleArrayUtils.toDouble(a[j]);
				failCounter.run(0, () -> {
					TestAssert.assertArrayEqualsRelative("Different b result", b1, b2, 1e-2);
				});
				failCounter.run(1, () -> {
					TestAssert.assertDoubleArrayEqualsRelative("Different a result", a1, a2, 1e-2);
				});
			}
			else
				c++;
		}
		if (c > ITER / 2)
			TestAssert.fail("Failed to solve %d / %d", c, ITER);
	}

	@Test
	public void solveLinearWithInversionIsNotFasterThanGaussJordanFloat()
	{
		TestAssume.assumeSpeedTest();

		final int ITER = 10000;
		ensureData(ITER);

		final GaussJordan solver = new GaussJordan();
		final EJMLLinearSolver solver2 = new EJMLLinearSolver();

		long t1 = Long.MAX_VALUE, t2 = Long.MAX_VALUE;

		for (int loops = 5; loops-- > 0;)
		{
			final ArrayList<float[][]> A = copyAfloat(SolverSpeedTest.Adata, ITER);
			final ArrayList<float[]> B = copyBfloat(SolverSpeedTest.Bdata, ITER);
			final ArrayList<double[]> A2 = copyA2double(SolverSpeedTest.Adata, ITER);
			final ArrayList<double[]> B2 = copyBdouble(SolverSpeedTest.Bdata, ITER);

			final long start1 = System.nanoTime();
			runFloat(A, B, ITER, solver);
			t1 = Math.min(t1, System.nanoTime() - start1);

			final long start2 = System.nanoTime();
			solveLinearWithInversion(A2, B2, ITER, solver2);
			t2 = Math.min(t2, System.nanoTime() - start2);
		}

		TestLog.logSpeedTestResult(t2 > t1,
				"GaussJordanFloat = %d : LinearSolver.solveLinearWithInversion = %d : %fx\n", t1, t2, (1.0 * t1) / t2);
	}

	@Test
	public void solveLinearIsFasterThanGaussJordanFloat()
	{
		TestAssume.assumeSpeedTest();

		final int ITER = 10000;
		ensureData(ITER);

		final GaussJordan solver = new GaussJordan();
		final EJMLLinearSolver solver2 = new EJMLLinearSolver();

		long t1 = Long.MAX_VALUE, t2 = Long.MAX_VALUE;

		for (int loops = 5; loops-- > 0;)
		{
			final ArrayList<float[][]> A = copyAfloat(SolverSpeedTest.Adata, ITER);
			final ArrayList<float[]> B = copyBfloat(SolverSpeedTest.Bdata, ITER);
			final ArrayList<double[]> A2 = copyA2double(SolverSpeedTest.Adata, ITER);
			final ArrayList<double[]> B2 = copyBdouble(SolverSpeedTest.Bdata, ITER);

			final long start1 = System.nanoTime();
			runFloat(A, B, ITER, solver);
			t1 = Math.min(t1, System.nanoTime() - start1);

			final long start2 = System.nanoTime();
			solveLinear(A2, B2, ITER, solver2);
			t2 = Math.min(t2, System.nanoTime() - start2);
		}

		TestLog.logSpeedTestResult(t2 < t1, "GaussJordanFloat = %d : LinearSolver.solveLinear = %d : %fx\n", t1,
				t2, (1.0 * t1) / t2);
	}

	protected void runFloat(ArrayList<float[][]> A, ArrayList<float[]> B, int ITER, GaussJordan solver)
	{
		for (int i = 0; i < ITER; i++)
			solver.solve(A.get(i), B.get(i));
	}

	@Test
	public void solveLinearWithInversionIsNotFasterThanGaussJordanDouble()
	{
		TestAssume.assumeSpeedTest();

		final int ITER = 10000;
		ensureData(ITER);

		final GaussJordan solver = new GaussJordan();
		final EJMLLinearSolver solver2 = new EJMLLinearSolver();

		long t1 = Long.MAX_VALUE, t2 = Long.MAX_VALUE;

		for (int loops = 5; loops-- > 0;)
		{
			final ArrayList<double[][]> A = copyAdouble(SolverSpeedTest.Adata, ITER);
			final ArrayList<double[]> B = copyBdouble(SolverSpeedTest.Bdata, ITER);
			final ArrayList<double[]> A2 = copyA2double(SolverSpeedTest.Adata, ITER);
			final ArrayList<double[]> B2 = copyBdouble(SolverSpeedTest.Bdata, ITER);

			final long start1 = System.nanoTime();
			solveGaussJordan(A, B, ITER, solver);
			t1 = Math.min(t1, System.nanoTime() - start1);

			final long start2 = System.nanoTime();
			solveLinearWithInversion(A2, B2, ITER, solver2);
			t2 = Math.min(t2, System.nanoTime() - start2);
		}

		TestLog.logSpeedTestResult(t2 > t1,
				"GaussJordanDouble = %d : LinearSolver.solveLinearWithInversion = %d : %fx\n", t1, t2, (1.0 * t1) / t2);
	}

	@Test
	public void solveLinearIsFasterThanGaussJordanDouble()
	{
		TestAssume.assumeSpeedTest();

		final int ITER = 10000;
		ensureData(ITER);

		final GaussJordan solver = new GaussJordan();
		final EJMLLinearSolver solver2 = new EJMLLinearSolver();

		long t1 = Long.MAX_VALUE, t2 = Long.MAX_VALUE;

		for (int loops = 5; loops-- > 0;)
		{
			final ArrayList<double[][]> A = copyAdouble(SolverSpeedTest.Adata, ITER);
			final ArrayList<double[]> B = copyBdouble(SolverSpeedTest.Bdata, ITER);
			final ArrayList<double[]> A2 = copyA2double(SolverSpeedTest.Adata, ITER);
			final ArrayList<double[]> B2 = copyBdouble(SolverSpeedTest.Bdata, ITER);

			final long start1 = System.nanoTime();
			solveGaussJordan(A, B, ITER, solver);
			t1 = Math.min(t1, System.nanoTime() - start1);

			final long start2 = System.nanoTime();
			solveLinear(A2, B2, ITER, solver2);
			t2 = Math.min(t2, System.nanoTime() - start2);
		}

		TestLog.logSpeedTestResult(t2 < t1, "GaussJordanDouble = %d : LinearSolver.solveLinear = %d : %fx\n", t1,
				t2, (1.0 * t1) / t2);
	}

	@Test
	public void solveCholeskyIsFasterThanGaussJordanDouble()
	{
		TestAssume.assumeSpeedTest();

		final int ITER = 10000;
		ensureData(ITER);

		final GaussJordan solver = new GaussJordan();
		final EJMLLinearSolver solver2 = new EJMLLinearSolver();

		long t1 = Long.MAX_VALUE, t2 = Long.MAX_VALUE;

		for (int loops = 5; loops-- > 0;)
		{
			final ArrayList<double[][]> A = copyAdouble(SolverSpeedTest.Adata, ITER);
			final ArrayList<double[]> B = copyBdouble(SolverSpeedTest.Bdata, ITER);
			final ArrayList<double[]> A2 = copyA2double(SolverSpeedTest.Adata, ITER);
			final ArrayList<double[]> B2 = copyBdouble(SolverSpeedTest.Bdata, ITER);

			final long start1 = System.nanoTime();
			solveGaussJordan(A, B, ITER, solver);
			t1 = Math.min(t1, System.nanoTime() - start1);

			final long start2 = System.nanoTime();
			solveCholesky(A2, B2, ITER, solver2);
			t2 = Math.min(t2, System.nanoTime() - start2);
		}

		TestLog.logSpeedTestResult(t2 < t1, "GaussJordanDouble = %d : LinearSolver.solveCholesky = %d : %fx\n", t1,
				t2, (1.0 * t1) / t2);
	}

	@Test
	public void solveCholeskyLDLTIsFasterThanGaussJordanDouble()
	{
		TestAssume.assumeSpeedTest();

		final int ITER = 10000;
		ensureData(ITER);

		final GaussJordan solver = new GaussJordan();
		final EJMLLinearSolver solver2 = new EJMLLinearSolver();

		long t1 = Long.MAX_VALUE, t2 = Long.MAX_VALUE;

		for (int loops = 5; loops-- > 0;)
		{
			final ArrayList<double[][]> A = copyAdouble(SolverSpeedTest.Adata, ITER);
			final ArrayList<double[]> B = copyBdouble(SolverSpeedTest.Bdata, ITER);
			final ArrayList<double[]> A2 = copyA2double(SolverSpeedTest.Adata, ITER);
			final ArrayList<double[]> B2 = copyBdouble(SolverSpeedTest.Bdata, ITER);

			final long start1 = System.nanoTime();
			solveGaussJordan(A, B, ITER, solver);
			t1 = Math.min(t1, System.nanoTime() - start1);

			final long start2 = System.nanoTime();
			solveCholeskyLDLT(A2, B2, ITER, solver2);
			t2 = Math.min(t2, System.nanoTime() - start2);
		}

		TestLog.logSpeedTestResult(t2 < t1, "GaussJordanDouble = %d : LinearSolver.solveCholeskyLDLT = %d : %fx\n",
				t1, t2, (1.0 * t1) / t2);
	}

	@Test
	public void solveIsFasterThanGaussJordanDouble()
	{
		TestAssume.assumeSpeedTest();

		final int ITER = 10000;
		ensureData(ITER);

		final GaussJordan solver = new GaussJordan();
		final EJMLLinearSolver solver2 = new EJMLLinearSolver();

		long t1 = Long.MAX_VALUE, t2 = Long.MAX_VALUE;

		for (int loops = 5; loops-- > 0;)
		{
			final ArrayList<double[][]> A = copyAdouble(SolverSpeedTest.Adata, ITER);
			final ArrayList<double[]> B = copyBdouble(SolverSpeedTest.Bdata, ITER);
			final ArrayList<double[]> A2 = copyA2double(SolverSpeedTest.Adata, ITER);
			final ArrayList<double[]> B2 = copyBdouble(SolverSpeedTest.Bdata, ITER);

			final long start1 = System.nanoTime();
			solveGaussJordan(A, B, ITER, solver);
			t1 = Math.min(t1, System.nanoTime() - start1);

			final long start2 = System.nanoTime();
			solve(A2, B2, ITER, solver2);
			t2 = Math.min(t2, System.nanoTime() - start2);
		}

		TestLog.logSpeedTestResult(t2 < t1, "GaussJordanDouble = %d : LinearSolver.solve = %d : %fx\n", t1, t2,
				(1.0 * t1) / t2);
	}

	private static boolean createData(RandomGenerator rand, float[][] alpha, float[] beta, boolean positiveDifinite)
	{
		// Generate a 2D Gaussian
		final SingleFreeCircularGaussian2DFunction func = new SingleFreeCircularGaussian2DFunction(10, 10);
		final double[] a = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
		a[Gaussian2DFunction.BACKGROUND] = 2 + rand.nextDouble() * 2;
		a[Gaussian2DFunction.SIGNAL] = 100 + rand.nextDouble() * 5;
		a[Gaussian2DFunction.X_POSITION] = 4.5 + rand.nextDouble();
		a[Gaussian2DFunction.Y_POSITION] = 4.5 + rand.nextDouble();
		a[Gaussian2DFunction.X_SD] = 1 + rand.nextDouble();
		a[Gaussian2DFunction.Y_SD] = 1 + rand.nextDouble();
		a[Gaussian2DFunction.ANGLE] = rand.nextDouble();

		final int[] x = new int[100];
		final double[] y = new double[100];
		func.initialise(a);
		for (int i = 0; i < x.length; i++)
			// Add random noise
			y[i] = func.eval(i) + ((rand.nextDouble() < 0.5) ? -rand.nextDouble() * 5 : rand.nextDouble() * 5);

		// Randomise parameters
		for (int i = 0; i < a.length; i++)
			a[i] += (rand.nextDouble() < 0.5) ? -rand.nextDouble() : rand.nextDouble();

		// Compute the Hessian and parameter gradient vector
		final GradientCalculator calc = new GradientCalculator(6);
		final double[][] alpha2 = new double[6][6];
		final double[] beta2 = new double[6];
		calc.findLinearised(y.length, y, a, alpha2, beta2, func);

		// Update the Hessian using a lambda shift
		final double lambda = 1.001;
		for (int i = 0; i < alpha2.length; i++)
			alpha2[i][i] *= lambda;

		// Copy back
		for (int i = 0; i < beta.length; i++)
		{
			beta[i] = (float) beta2[i];
			for (int j = 0; j < beta.length; j++)
				alpha[i][j] = (float) alpha2[i][j];
		}

		// Check for a positive definite matrix
		if (positiveDifinite)
		{
			final EJMLLinearSolver solver = new EJMLLinearSolver();
			return solver.solveCholeskyLDLT(copydouble(alpha), copydouble(beta));
		}

		return true;
	}

	private static ArrayList<float[][]> copyAfloat(ArrayList<float[][]> a, int iter)
	{
		iter = FastMath.min(a.size(), iter);
		final ArrayList<float[][]> a2 = new ArrayList<>(iter);
		for (int i = 0; i < iter; i++)
			a2.add(copyfloat(a.get(i)));
		return a2;
	}

	private static float[][] copyfloat(float[][] d)
	{
		final float[][] d2 = new float[d.length][d.length];
		for (int i = 0; i < d.length; i++)
			for (int j = 0; j < d.length; j++)
				d2[i][j] = d[i][j];
		return d2;
	}

	private static ArrayList<float[]> copyBfloat(ArrayList<float[]> b, int iter)
	{
		iter = FastMath.min(b.size(), iter);
		final ArrayList<float[]> b2 = new ArrayList<>(iter);
		for (int i = 0; i < iter; i++)
			b2.add(Arrays.copyOf(b.get(i), b.get(i).length));
		return b2;
	}

	private static ArrayList<double[][]> copyAdouble(ArrayList<float[][]> a, int iter)
	{
		iter = FastMath.min(a.size(), iter);
		final ArrayList<double[][]> a2 = new ArrayList<>(iter);
		for (int i = 0; i < iter; i++)
			a2.add(copydouble(a.get(i)));
		return a2;
	}

	private static ArrayList<double[]> copyA2double(ArrayList<float[][]> a, int iter)
	{
		iter = FastMath.min(a.size(), iter);
		final ArrayList<double[]> a2 = new ArrayList<>(iter);
		for (int i = 0; i < iter; i++)
			a2.add(new DenseMatrix64F(copydouble(a.get(i))).data);
		return a2;
	}

	private static double[][] copydouble(float[][] d)
	{
		final double[][] d2 = new double[d.length][d.length];
		for (int i = 0; i < d.length; i++)
			for (int j = 0; j < d.length; j++)
				d2[i][j] = d[i][j];
		return d2;
	}

	private static ArrayList<double[]> copyBdouble(ArrayList<float[]> b, int iter)
	{
		iter = FastMath.min(b.size(), iter);
		final ArrayList<double[]> b2 = new ArrayList<>(iter);
		for (int i = 0; i < iter; i++)
			b2.add(copydouble(b.get(i)));
		return b2;
	}

	private static double[] copydouble(float[] d)
	{
		final double[] d2 = new double[d.length];
		for (int i = 0; i < d.length; i++)
			d2[i] = d[i];
		return d2;
	}

	protected void solveGaussJordan(ArrayList<double[][]> A, ArrayList<double[]> B, int ITER, GaussJordan solver)
	{
		ITER = FastMath.min(ITER, A.size());
		for (int i = 0; i < ITER; i++)
			solver.solve(A.get(i), B.get(i));
	}

	protected void solveLinearWithInversion(ArrayList<double[]> A, ArrayList<double[]> B, int ITER,
			EJMLLinearSolver solver)
	{
		ITER = FastMath.min(ITER, A.size());
		for (int i = 0; i < ITER; i++)
		{
			final double[] a = A.get(i);
			solver.solveLinear(a, B.get(i));
			solver.invertLastA(a);
		}
	}

	protected void solveLinear(ArrayList<double[]> A, ArrayList<double[]> B, int ITER, EJMLLinearSolver solver)
	{
		ITER = FastMath.min(ITER, A.size());
		for (int i = 0; i < ITER; i++)
			solver.solveLinear(A.get(i), B.get(i));
	}

	protected void solveCholesky(ArrayList<double[]> A, ArrayList<double[]> B, int ITER, EJMLLinearSolver solver)
	{
		ITER = FastMath.min(ITER, A.size());
		for (int i = 0; i < ITER; i++)
			solver.solveCholesky(A.get(i), B.get(i));
	}

	protected void solveCholeskyLDLT(ArrayList<double[]> A, ArrayList<double[]> B, int ITER, EJMLLinearSolver solver)
	{
		ITER = FastMath.min(ITER, A.size());
		for (int i = 0; i < ITER; i++)
			solver.solveCholeskyLDLT(A.get(i), B.get(i));
	}

	protected void solve(ArrayList<double[]> A, ArrayList<double[]> B, int ITER, EJMLLinearSolver solver)
	{
		ITER = FastMath.min(ITER, A.size());
		for (int i = 0; i < ITER; i++)
			solver.solve(A.get(i), B.get(i));
	}

	void log(String format, Object... args)
	{
		System.out.printf(format, args);
	}
}
