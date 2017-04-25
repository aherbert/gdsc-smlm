package gdsc.smlm.fitting.linear;

import org.ejml.factory.LinearSolver;
import org.ejml.factory.LinearSolverFactory;
import org.ejml.ops.CommonOps;

import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.Maths;

import org.ejml.alg.dense.linsol.chol.LinearSolverCholLDL;
import org.ejml.data.DenseMatrix64F;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Solves (one) linear equation, a x = b.
 * <p>
 * Wraps the LinearSolver class from the EJML (Efficient Java Matrix Library) library.
 */
public class EJMLLinearSolver
{
	private LinearSolver<DenseMatrix64F> linearSolver;
	private LinearSolver<DenseMatrix64F> choleskySolver;
	private LinearSolver<DenseMatrix64F> choleskyLDLTSolver;
	private LinearSolver<DenseMatrix64F> lastSuccessfulSolver = null;
	private DenseMatrix64F X;
	private DenseMatrix64F A_inv;
	private int solverSize = 0;
	private boolean errorChecking = false;
	private DoubleEquality equal = null;

	// Variables to handle zero entries in input vector B:

	// The number of entries in vector B which are zero
	private int zeroCount = 0;
	// Array holding which columns in the input vector B are zero 
	private boolean[] zeros = new boolean[0];
	// Array mapping the solution index X back to the original index in B
	private int[] zerosMap = new int[0];

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On input have a[n][n], b[n]. On output these replaced by a_inverse[n][n], x[n].
	 * 
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solveLinearWithInversion(float[][] a, float[] b)
	{
		// Convert data to double
		double[] a2 = new double[b.length * b.length];
		double[] b2 = new double[b.length];
		for (int i = 0, index = 0; i < b.length; i++)
		{
			b2[i] = b[i];
			for (int j = 0; j < b.length; j++)
			{
				a2[index++] = a[i][j];
			}
		}

		DenseMatrix64F A = new DenseMatrix64F(b.length, b.length, true, a2);

		createSolver(b.length);

		if (!initialiseSolver(linearSolver, A))
			return false;

		DenseMatrix64F B = DenseMatrix64F.wrap(b.length, 1, b2);
		DenseMatrix64F x = (linearSolver.modifiesB()) ? X : B;
		linearSolver.solve(B, x);

		zeroCount = 0;
		invert(a);

		// Copy back result
		for (int i = 0; i < b.length; i++)
		{
			b[i] = (float) x.data[i];
		}

		return true;
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On input have a[n][n], b[n]. On output these replaced by a_inverse[n][n], x[n].
	 * 
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solveLinearWithInversion(double[][] a, double[] b)
	{
		DenseMatrix64F A = new DenseMatrix64F(a);

		createSolver(b.length);

		if (!initialiseSolver(linearSolver, A))
			return false;

		DenseMatrix64F B = DenseMatrix64F.wrap(b.length, 1, b);
		DenseMatrix64F x = (linearSolver.modifiesB()) ? X : B;
		linearSolver.solve(B, x);

		zeroCount = 0;
		invert(a);

		// Copy back result
		if (linearSolver.modifiesB())
		{
			for (int i = 0; i < b.length; i++)
			{
				b[i] = x.data[i];
			}
		}

		return true;
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On input have a[n][n], b[n]. On output b replaced by x[n].
	 * 
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solveLinear(float[][] a, float[] b)
	{
		// Convert data to double
		DenseMatrix64F[] data = wrapData(a, b);

		if (solve(linearSolver, data[0], data[1]))
		{
			for (int i = 0; i < b.length; i++)
				b[i] = (float) data[1].data[i];
			return true;
		}

		return false;
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On input have a[n][n], b[n]. On output b replaced by x[n].
	 * 
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solveCholesky(float[][] a, float[] b)
	{
		// Convert data to double
		DenseMatrix64F[] data = wrapData(a, b);

		if (solve(choleskySolver, data[0], data[1]))
		{
			for (int i = 0; i < b.length; i++)
				b[i] = (float) data[1].data[i];
			return true;
		}

		return false;
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On input have a[n][n], b[n]. On output b replaced by x[n].
	 * 
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solveCholeskyLDLT(float[][] a, float[] b)
	{
		// Convert data to double
		DenseMatrix64F[] data = wrapData(a, b);

		if (solve(choleskyLDLTSolver, data[0], data[1]))
		{
			for (int i = 0; i < b.length; i++)
				b[i] = (float) data[1].data[i];
			return true;
		}

		return false;
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On input have a[n][n], b[n]. On output b replaced by x[n].
	 * 
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solveLinear(double[][] a, double[] b)
	{
		createSolver(b.length);
		return solve(linearSolver, new DenseMatrix64F(a), DenseMatrix64F.wrap(b.length, 1, b));
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On input have a[n][n], b[n]. On output b replaced by x[n].
	 * 
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solveCholesky(double[][] a, double[] b)
	{
		createSolver(b.length);
		return solve(choleskySolver, new DenseMatrix64F(a), DenseMatrix64F.wrap(b.length, 1, b));
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On input have a[n][n], b[n]. On output b replaced by x[n].
	 * 
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solveCholeskyLDLT(double[][] a, double[] b)
	{
		createSolver(b.length);
		return solve(choleskyLDLTSolver, new DenseMatrix64F(a), DenseMatrix64F.wrap(b.length, 1, b));
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On input have a[n][n], b[n]. On output b replaced by x[n].
	 * 
	 * @return False if the equation is singular (no solution)
	 */
	private boolean solve(LinearSolver<DenseMatrix64F> solver, DenseMatrix64F A, DenseMatrix64F B)
	{
		zeroCount = 0;
		DenseMatrix64F x = (solver.modifiesB()) ? X : B;

		if (!solve(solver, A, B, x))
			return false;

		// Copy back result if necessary
		if (solver.modifiesB())
		{
			for (int i = 0; i < B.numRows; i++)
			{
				B.data[i] = x.data[i];
			}
		}

		return true;
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On input have a[n][n], b[n]. On output x[n].
	 * <p>
	 * If checking the solution then A and/or B will not be modified. If no solution checking is enabled then A and/or B
	 * could be modified.
	 * 
	 * @param solver
	 * @param A
	 * @param B
	 * @param X
	 * @return False if the equation is singular (no solution)
	 */
	private boolean solve(LinearSolver<DenseMatrix64F> solver, DenseMatrix64F A, DenseMatrix64F B, DenseMatrix64F X)
	{
		if (errorChecking)
		{
			if (!solveSafe(solver, A, B, X))
				return false;

			return validate(A, X, B);
		}
		else
			return solveUnsafe(solver, A, B, X);
	}

	/**
	 * Check that the solution for x satisfies A x = b within the error tolerance
	 * 
	 * @param A
	 * @param x
	 * @param b
	 * @return
	 */
	private boolean validate(DenseMatrix64F A, DenseMatrix64F x, DenseMatrix64F b)
	{
		if (errorChecking)
		{
			// Compute A x = b

			// Used for debugging
			//DenseMatrix64F b2 = new DenseMatrix64F(size, 1);
			//CommonOps.mult(A, x, b2);
			//return equal.almostEqualComplement(b.data, b2.data);

			for (int i = 0, index = 0; i < b.numRows; i++)
			{
				double bi = 0;
				for (int j = 0; j < b.numRows; j++)
				{
					bi += A.data[index++] * x.data[j];
				}
				if (!equal.almostEqualComplement(b.data[i], bi))
				{
					//System.out.printf("Bad solution: %g != %g (%g = %d)\n", b.data[i], bi, 
					//		equal.relativeError(b.data[i], bi), equal.complement(b.data[i], bi));
					return false;
				}
			}
		}
		return true;
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On input have a[n][n], b[n]. On output x[n].
	 * <p>
	 * A and/or B will not be modified. If you do not care then use
	 * {@link #solveUnsafe(LinearSolver, DenseMatrix64F, DenseMatrix64F, DenseMatrix64F)}
	 * 
	 * @param solver
	 * @param A
	 * @param B
	 * @param X
	 * @return False if the equation is singular (no solution)
	 */
	private boolean solveSafe(LinearSolver<DenseMatrix64F> solver, DenseMatrix64F A, DenseMatrix64F B, DenseMatrix64F X)
	{
		if (solver.modifiesA())
			A = A.copy();
		if (solver.modifiesB())
			B = B.copy();

		if (!initialiseSolver(solver, A))
			return false;

		solver.solve(B, X);

		return true;
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On input have a[n][n], b[n]. On output x[n].
	 * <p>
	 * A and/or B may be modified. Check the solver before calling or use
	 * {@link #solveSafe(LinearSolver, DenseMatrix64F, DenseMatrix64F, DenseMatrix64F)}
	 * 
	 * @param solver
	 * @param A
	 * @param B
	 * @param X
	 * @return False if the equation is singular (no solution)
	 */
	private boolean solveUnsafe(LinearSolver<DenseMatrix64F> solver, DenseMatrix64F A, DenseMatrix64F B,
			DenseMatrix64F X)
	{
		if (!initialiseSolver(solver, A))
			return false;

		solver.solve(B, X);

		return true;
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On input have a[n][n], b[n]. On output b replaced by x[n].
	 * <p>
	 * Solve using the CholeskyLDLT method or if that fails the Linear decomposition.
	 * 
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solve(double[][] a, double[] b)
	{
		zeroCount = 0;
		DenseMatrix64F A = new DenseMatrix64F(a);
		DenseMatrix64F B = DenseMatrix64F.wrap(b.length, 1, b);

		createSolver(b.length);

		// Note: Use solveSafe as we need the A and B matrix for the subsequent 
		// solve attempt if failure

		if (solveSafe(choleskyLDLTSolver, A, B, X))
		{
			if (validate(A, X, B))
			{
				for (int i = 0; i < b.length; i++)
					b[i] = X.data[i];
				return true;
			}
		}

		// No need for an explicit solveSafe this time. Use the solve() method which will include validate().
		if (solve(linearSolver, A, B, X))
		{
			for (int i = 0; i < b.length; i++)
				b[i] = X.data[i];
			return true;
		}

		return false;
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On input have a[n][n], b[n]. On output b replaced by x[n].
	 * <p>
	 * Solve using the CholeskyLDLT method or if that fails the Linear decomposition.
	 * 
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solve(float[][] a, float[] b)
	{
		zeroCount = 0;
		DenseMatrix64F[] data = wrapData(a, b);
		DenseMatrix64F A = data[0];
		DenseMatrix64F B = data[1];

		createSolver(b.length);

		// Note: Use solveSafe as we need the A and B matrix for the subsequent 
		// solve attempt if failure

		if (solveSafe(choleskyLDLTSolver, A, B, X))
		{
			if (validate(A, X, B))
			{
				for (int i = 0; i < b.length; i++)
					b[i] = (float) X.data[i];
				return true;
			}
		}

		// No need for an explicit solveSafe this time. Use the solve() method which will include validate().
		if (solve(linearSolver, A, B, X))
		{
			for (int i = 0; i < b.length; i++)
				b[i] = (float) X.data[i];
			return true;
		}

		return false;
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On input have a[n][n], b[n]. On output b replaced by x[n].
	 * <p>
	 * Solve using the CholeskyLDLT method or if that fails the Linear decomposition.
	 * <p>
	 * Note: Any zero elements in b are not solved.
	 * 
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solveWithZeros(double[][] a, double[] b)
	{
		checkForZeros(b);
		if (zeroCount == 0)
			return solve(a, b);

		if (b.length == zeroCount)
			// This can be solved using all zeros but it cannot be inverted
			return false;

		createSolver(b.length - zeroCount);

		DenseMatrix64F[] data = wrapDataWithZeros(a, b);
		DenseMatrix64F A = data[0];
		DenseMatrix64F B = data[1];
		DenseMatrix64F X = new DenseMatrix64F(B.numRows, 1);

		// Note: Use solveSafe as we need the A and B matrix for the subsequent 
		// solve attempt if failure

		if (solveSafe(choleskyLDLTSolver, A, B, X))
		{
			if (validate(A, X, B))
			{
				for (int i = 0; i < X.data.length; i++)
					b[zerosMap[i]] = X.data[i];
				return true;
			}
		}

		// No need for an explicit solveSafe this time. Use the solve() method which will include validate().
		if (solve(linearSolver, A, B, X))
		{
			for (int i = 0; i < X.data.length; i++)
				b[zerosMap[i]] = X.data[i];
			return true;
		}

		return false;
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On input have a[n][n], b[n]. On output b replaced by x[n].
	 * <p>
	 * Solve using the CholeskyLDLT method or if that fails the Linear decomposition.
	 * <p>
	 * Note: Any zero elements in b are not solved.
	 * 
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solveWithZeros(float[][] a, float[] b)
	{
		checkForZeros(b);
		if (zeroCount == 0)
			return solve(a, b);

		if (b.length == zeroCount)
			// This can be solved using all zeros but it cannot be inverted
			return false;

		createSolver(b.length - zeroCount);

		DenseMatrix64F[] data = wrapDataWithZeros(a, b);
		DenseMatrix64F A = data[0];
		DenseMatrix64F B = data[1];
		DenseMatrix64F X = new DenseMatrix64F(B.numRows, 1);

		// Note: Use solveSafe as we need the A and B matrix for the subsequent 
		// solve attempt if failure

		if (solveSafe(choleskyLDLTSolver, A, B, X))
		{
			if (validate(A, X, B))
			{
				for (int i = 0; i < X.data.length; i++)
					b[zerosMap[i]] = (float) X.data[i];
				return true;
			}
		}

		// No need for an explicit solveSafe this time. Use the solve() method which will include validate().
		if (solve(linearSolver, A, B, X))
		{
			for (int i = 0; i < X.data.length; i++)
				b[zerosMap[i]] = (float) X.data[i];
			return true;
		}

		return false;
	}

	private void checkForZeros(float[] b)
	{
		zeroCount = 0;
		for (int i = 0; i < b.length; i++)
			if (isZero(b[i]))
				zeroCount++;
	}

	private boolean isZero(float f)
	{
		return f == 0; // Q. Should this be an fEQ operation?
	}

	private void createZeroArrays(float[] b)
	{
		// Resize array if necessary
		if (zeros.length < b.length)
		{
			zeros = new boolean[b.length];
			zerosMap = new int[b.length];
		}
		for (int i = 0, j = 0; i < b.length; i++)
		{
			if (isZero(b[i]))
			{
				zeros[i] = true;
			}
			else
			{
				zeros[i] = false;
				// Map new data index back to original index
				zerosMap[j++] = i;
			}
		}
	}

	private void checkForZeros(double[] b)
	{
		zeroCount = 0;
		for (int i = 0; i < b.length; i++)
			if (isZero(b[i]))
				zeroCount++;
	}

	private boolean isZero(double f)
	{
		return f == 0; // Q. Should this be an fEQ operation?
	}

	private void createZeroArrays(double[] b)
	{
		// Resize array if necessary
		if (zeros.length < b.length)
		{
			zeros = new boolean[b.length];
			zerosMap = new int[b.length];
		}
		for (int i = 0, j = 0; i < b.length; i++)
		{
			if (isZero(b[i]))
			{
				zeros[i] = true;
			}
			else
			{
				zeros[i] = false;
				// Map new data index back to original index
				zerosMap[j++] = i;
			}
		}
	}

	private void checkForZerosAndCreateZeroArrays(double[][] a)
	{
		zeroCount = 0;
		// Resize array if necessary
		if (zeros.length < a.length)
		{
			zeros = new boolean[a.length];
			zerosMap = new int[a.length];
		}
		for (int i = 0, j = 0; i < a.length; i++)
		{
			// Check diagonal element first
			if (isZero(a[i][i]))
			{
				// Check row and column entries
				if (isZero(a, i))
				{
					zeros[i] = true;
					zeroCount++;
					continue;
				}
			}

			zeros[i] = false;
			// Map new data index back to original index
			zerosMap[j++] = i;
		}
	}

	private boolean isZero(double[][] a, int i)
	{
		final int length = a.length;
		for (int k = length; k-- > 0;)
		{
			if (!isZero(a[i][k]))
				return false;
			if (!isZero(a[k][i]))
				return false;
		}
		return true;
	}

	/**
	 * Computes the inverse of the 'A' matrix passed into the last successful solve method.
	 * <p>
	 * On output a[n][n] replaced by the inverse of the solved matrix a. If any column/row index was removed (as it was
	 * set to zero in the input matrix) then the resulting column/row index will be set to zero.
	 * 
	 * @return False if the last solve attempt failed, or inversion produces non finite values
	 */
	public boolean invert(double[][] a)
	{
		if (lastSuccessfulSolver == null)
			return false;

		lastSuccessfulSolver.invert(A_inv);

		// Check for NaN or Infinity
		for (int i = A_inv.data.length; i-- > 0;)
			if (!Maths.isFinite(A_inv.data[i]))
				return false;

		if (zeroCount > 0)
		{
			int size = this.solverSize + zeroCount;
			for (int i = 0, index = 0; i < size; i++)
			{
				if (zeros[i])
				{
					// Set all zero columns to zero.
					for (int j = 0; j < size; j++)
						a[i][j] = a[j][i] = 0;
					continue;
				}
				for (int j = 0; j < size; j++)
				{
					if (zeros[j])
						continue;
					a[i][j] = A_inv.data[index++];
				}
			}
		}
		else
		{
			for (int i = 0, index = 0; i < solverSize; i++)
				for (int j = 0; j < solverSize; j++)
				{
					a[i][j] = A_inv.data[index++];
				}
		}

		return true;
	}

	/**
	 * Computes the inverse of the 'A' matrix passed into the last successful solve method.
	 * <p>
	 * On output a[n][n] replaced by the inverse of the solved matrix a. If any column/row index was removed (as it was
	 * set to zero in the input matrix) then the resulting column/row index will be set to zero.
	 * 
	 * @return False if the last solve attempt failed, or inversion produces non finite values
	 */
	public boolean invert(float[][] a)
	{
		if (lastSuccessfulSolver == null)
			return false;

		lastSuccessfulSolver.invert(A_inv);

		// Check for NaN or Infinity
		for (int i = A_inv.data.length; i-- > 0;)
			if (!Maths.isFinite(A_inv.data[i]))
				return false;

		if (zeroCount > 0)
		{
			int size = this.solverSize + zeroCount;
			for (int i = 0, index = 0; i < size; i++)
			{
				if (zeros[i])
				{
					// Set all zero columns to zer0
					for (int j = 0; j < size; j++)
						a[i][j] = a[j][i] = 0;
					continue;
				}
				for (int j = 0; j < size; j++)
				{
					if (zeros[j])
						continue;
					a[i][j] = (float) A_inv.data[index++];
				}
			}
		}
		else
		{
			for (int i = 0, index = 0; i < solverSize; i++)
				for (int j = 0; j < solverSize; j++)
				{
					a[i][j] = (float) A_inv.data[index++];
				}
		}

		return true;
	}

	private boolean initialiseSolver(LinearSolver<DenseMatrix64F> solver, DenseMatrix64F A)
	{
		lastSuccessfulSolver = null;
		if (!solver.setA(A))
			return false;

		lastSuccessfulSolver = solver;

		// Note: 
		// An analysis of the quality measure shows that the LinearSolverCholLDL always returns 
		// a quality of 1 when it has successfully setA(). The quality is the product of the 
		// diagonal and the LDLT method must make this equal to 1 on all elements.
		// The other linear solver is a LU decomposition. The quality is the product of the diagonal
		// with each element divided by the maximum matrix element value. This can return a quality 
		// down in the 1e-20 range but the solution is still valid for A x = b.
		//
		// Thus we will check the solution to our equation instead for a maximum relative error.

		//double q = solver.quality(); 
		//System.out.println(solver.toString() + " quality = "+q);
		//if (q <= 1e-8)
		//{
		//System.out.println(solver.toString() + " poor quality = " + q); // + "\n"+A.toString());
		//return false;
		//}

		return true;
	}

	private void createSolver(int length)
	{
		if (solverSize != length)
		{
			solverSize = length;

			choleskyLDLTSolver = new LinearSolverCholLDL();

			// This is a Cholesky solver that only works on symmetric positive definite matrices
			choleskySolver = LinearSolverFactory.symmPosDef(solverSize);

			// This should work on any matrix
			linearSolver = LinearSolverFactory.linear(solverSize);

			X = new DenseMatrix64F(solverSize, 1);
			A_inv = new DenseMatrix64F(solverSize, solverSize);
		}
	}

	private DenseMatrix64F[] wrapData(float[][] a, float[] b)
	{
		double[] a2 = new double[b.length * b.length];
		double[] b2 = new double[b.length];
		for (int i = 0, index = 0; i < b.length; i++)
		{
			b2[i] = b[i];
			for (int j = 0; j < b.length; j++)
			{
				a2[index++] = a[i][j];
			}
		}
		return new DenseMatrix64F[] { DenseMatrix64F.wrap(b.length, b.length, a2),
				DenseMatrix64F.wrap(b.length, 1, b2) };
	}

	private DenseMatrix64F[] wrapDataWithZeros(float[][] a, float[] b)
	{
		createZeroArrays(b);

		int size = b.length - zeroCount;
		double[] a2 = new double[size * size];
		double[] b2 = new double[size];
		for (int i = 0, ii = 0, index = 0; i < b.length; i++)
		{
			if (zeros[i])
				continue;

			b2[ii++] = b[i];
			for (int j = 0; j < b.length; j++)
			{
				if (zeros[j])
					continue;

				a2[index++] = a[i][j];
			}
		}
		return new DenseMatrix64F[] { DenseMatrix64F.wrap(size, size, a2), DenseMatrix64F.wrap(size, 1, b2) };
	}

	private DenseMatrix64F[] wrapDataWithZeros(double[][] a, double[] b)
	{
		createZeroArrays(b);

		int size = b.length - zeroCount;
		double[] a2 = new double[size * size];
		double[] b2 = new double[size];
		for (int i = 0, ii = 0, index = 0; i < b.length; i++)
		{
			if (zeros[i])
				continue;

			b2[ii++] = b[i];
			for (int j = 0; j < b.length; j++)
			{
				if (zeros[j])
					continue;

				a2[index++] = a[i][j];
			}
		}
		return new DenseMatrix64F[] { DenseMatrix64F.wrap(size, size, a2), DenseMatrix64F.wrap(size, 1, b2) };
	}

	private DenseMatrix64F wrapDataWithZeros(double[][] a)
	{
		checkForZerosAndCreateZeroArrays(a);

		if (zeroCount == 0)
			return new DenseMatrix64F(a);

		final int length = a.length;

		int size = length - zeroCount;
		double[] a2 = new double[size * size];
		for (int i = 0, index = 0; i < length; i++)
		{
			if (zeros[i])
				continue;
			for (int j = 0; j < length; j++)
			{
				if (zeros[j])
					continue;
				a2[index++] = a[i][j];
			}
		}
		return DenseMatrix64F.wrap(size, size, a2);
	}

	/**
	 * @param equal
	 *            the equality class to compare that the solution x in A x = b is within tolerance
	 */
	public void setEqual(DoubleEquality equal)
	{
		this.equal = equal;
		errorChecking = (equal != null);
	}

	/**
	 * @return the equality class to compare that the solution x in A x = b is within tolerance
	 */
	public DoubleEquality getEqual()
	{
		return equal;
	}

	/**
	 * Computes the inverse of the symmetric positive definite matrix. On output a is replaced by the inverse of a.
	 * <p>
	 * Note: Rows and columns for element i that are entirely zero are removed before inversion and the result will have
	 * zero in the result.
	 *
	 * @param a
	 *            the matrix a
	 * @return False if there is no solution
	 */
	public boolean invertSymmPosDef(double[][] a)
	{
		DenseMatrix64F A = wrapDataWithZeros(a);
		DenseMatrix64F A2 = (errorChecking) ? A.copy() : null;
		createSolver(A.numCols);
		if (!initialiseSolver(choleskySolver, A))
		{
			if (choleskySolver.modifiesA())
				A = wrapDataWithZeros(a);
			if (!initialiseSolver(linearSolver, A))
				return false;
		}

		if (!invert(a))
			return false;

		if (errorChecking)
			return validateInversion(A2);

		return true;
	}

	private boolean validateInversion(DenseMatrix64F A)
	{
		// Compute A A_inv = I
		final int n = A.numCols;
		DenseMatrix64F I = new DenseMatrix64F(n, n);
		CommonOps.mult(A, A_inv, I);

		for (int i = n, index = I.data.length; i-- > 0;)
		{
			for (int j = n; j-- > 0;)
			{
				double target = (j == i) ? 1 : 0;
				if (!equal.almostEqualComplement(I.data[--index], target))
				{
					System.out.printf("Bad solution: %g != %g (%g = %d)\n", I.data[index], target,
							DoubleEquality.relativeError(I.data[index], target),
							DoubleEquality.complement(I.data[index], target));
					return false;
				}
			}
		}
		return true;
	}

	/**
	 * Computes the inverse of the matrix. This should work on any square matrix. On output a is replaced by the
	 * inverse of a.
	 * <p>
	 * Note: Rows and columns for element i that are entirely zero are removed before inversion and the result will have
	 * zero in the result.
	 * 
	 * @param a
	 *            the matrix a
	 * @return False if there is no solution
	 */
	public boolean invertLinear(double[][] a)
	{
		DenseMatrix64F A = wrapDataWithZeros(a);
		DenseMatrix64F A2 = (errorChecking) ? A.copy() : null;
		createSolver(A.numCols);
		if (!initialiseSolver(linearSolver, A))
			return false;

		if (!invert(a))
			return false;

		if (errorChecking)
			return validateInversion(A2);

		return true;
	}
}
