package gdsc.smlm.fitting.linear;

import org.ejml.factory.LinearSolver;
import org.ejml.factory.LinearSolverFactory;
import org.ejml.ops.CommonOps;

import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.Maths;

import org.ejml.alg.dense.linsol.chol.LinearSolverCholLDL;
import org.ejml.alg.dense.misc.UnrolledInverseFromMinor;
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
 * <p>
 * This class assumes that the input matrix is symmetric positive definite matrices, for instance it has been generated
 * in the numerical solution of partial differential equations (for example a Hessian matrix).
 * <p>
 * If zeros occur on the diagonal of the matrix, for example if no gradient is available for a parameter, then the index
 * can be excised from the matrix before solving.
 */
public class EJMLLinearSolver
{
	// TODO - 
	// Rewrite to accept DenseMatrix64 as the primary input with wrapper 
	// functions to accept double[][].
	// Do not worry about double[] since this is easily wrapped.

	/**
	 * Solve the matrix using direct inversion
	 */
	private class InversionSolver implements LinearSolver<DenseMatrix64F>
	{
		private DenseMatrix64F A;

		public boolean setA(DenseMatrix64F A)
		{
			if (A.numCols <= UnrolledInverseFromMinor.MAX)
			{
				// Direct inversion using the determinant 
				if (A.numCols >= 2)
				{
					UnrolledInverseFromMinor.inv(A, A);
				}
				else
				{
					A.set(0, 1.0 / A.get(0));
				}

				// Check for NaN or Infinity
				for (int i = A.data.length; i-- > 0;)
					if (!Maths.isFinite(A.data[i]))
						return false;

				this.A = A;
				return true;
			}
			return false;
		}

		public double quality()
		{
			return 0;
		}

		public void solve(DenseMatrix64F B, DenseMatrix64F X)
		{
			CommonOps.mult(A, B, X);
		}

		public void invert(DenseMatrix64F A_inv)
		{
			System.arraycopy(A.data, 0, A_inv.data, 0, A.data.length);
		}

		public boolean modifiesA()
		{
			return true;
		}

		public boolean modifiesB()
		{
			return false;
		}
	}

	private LinearSolver<DenseMatrix64F> linearSolver;
	private LinearSolver<DenseMatrix64F> pseudoInverseSolver;
	private LinearSolver<DenseMatrix64F> choleskySolver;
	private LinearSolver<DenseMatrix64F> choleskyLDLTSolver;
	private LinearSolver<DenseMatrix64F> inversionSolver;
	private LinearSolver<DenseMatrix64F> lastSuccessfulSolver = null;
	private DenseMatrix64F X;
	private DenseMatrix64F A_inv;
	private int solverSize = 0;
	private boolean errorChecking = false;
	private DoubleEquality equal = null;

	/**
	 * Instantiates a new EJML linear solver.
	 */
	public EJMLLinearSolver()
	{
	}

	/**
	 * Instantiates a new EJML linear solver with tolerance for the linear solution.
	 *
	 * @param equal
	 *            the object for equality
	 */
	public EJMLLinearSolver(DoubleEquality equal)
	{
		setEqual(equal);
	}

	/**
	 * Instantiates a new EJML linear solver with tolerance for the linear solution
	 *
	 * @param significantDigits
	 *            the significant digits for equality
	 * @param maxAbsoluteError
	 *            the max absolute error for equality
	 */
	public EJMLLinearSolver(int significantDigits, double maxAbsoluteError)
	{
		this(new DoubleEquality(significantDigits, maxAbsoluteError));
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
		return solve(getLinearSolver(), new DenseMatrix64F(a), DenseMatrix64F.wrap(b.length, 1, b));
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
		return solve(getCholeskySolver(), new DenseMatrix64F(a), DenseMatrix64F.wrap(b.length, 1, b));
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
		return solve(getCholeskyLDLTSolver(), new DenseMatrix64F(a), DenseMatrix64F.wrap(b.length, 1, b));
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On input have a[n][n], b[n]. On output b replaced by x[n].
	 * 
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solvePseudoInverse(double[][] a, double[] b)
	{
		createSolver(b.length);
		return solve(getPseudoInverseSolver(), new DenseMatrix64F(a), DenseMatrix64F.wrap(b.length, 1, b));
	}

	/**
	 * Solves (one) linear equation, a x = b by direct inversion. Works on small matrices up to size 5.
	 * <p>
	 * On input have a[n][n], b[n]. On output b replaced by x[n].
	 * 
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solveDirectInversion(double[][] a, double[] b)
	{
		createSolver(b.length);
		return solve(getInversionSolver(), new DenseMatrix64F(a), DenseMatrix64F.wrap(b.length, 1, b));
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
		boolean copy = (errorChecking || solver.modifiesB());
		DenseMatrix64F x = (copy) ? getX() : B;

		if (!solve(solver, A, B, x))
			return false;

		// Copy back result if necessary
		if (copy)
			System.arraycopy(x.data, 0, B.data, 0, B.numRows);

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
					//		DoubleEquality.relativeError(b.data[i], bi), DoubleEquality.complement(b.data[i], bi));
					return false;
				}
			}
			//System.out.println("OK");
		}
		return true;
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On input have a[n][n], b[n]. On output b replaced by x[n].
	 * <p>
	 * Solve using the CholeskyLDLT method or, if that fails (due to a singular matrix), the PseudoInversion
	 * decomposition.
	 * 
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solve(double[][] a, double[] b)
	{
		DenseMatrix64F A = new DenseMatrix64F(a);
		DenseMatrix64F B = DenseMatrix64F.wrap(b.length, 1, b);

		createSolver(b.length);

		// Note: Use solveSafe as we need the A and B matrix for the subsequent 
		// solve attempt if failure

		if (solveSafe(getCholeskyLDLTSolver(), A, B, getX()))
		{
			if (validate(A, X, B))
			{
				System.arraycopy(X.data, 0, b, 0, b.length);
				return true;
			}
		}

		// TODO - Count how often the primary method fails on a realistic set of fitting data
		// since the PseudoInverse method is slow. We may want to try a different solver first,
		// e.g. LinearSolver.

		// No need for an explicit solveSafe this time. 
		// Use the solve() method which will include validate().
		if (solve(getPseudoInverseSolver(), A, B, X))
		{
			System.arraycopy(X.data, 0, b, 0, b.length);
			return true;
		}

		return false;
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

		lastSuccessfulSolver.invert(getA_inv());

		return checkA_inv(a);
	}

	private boolean checkA_inv(double[][] a)
	{
		// Check for NaN or Infinity
		for (int i = A_inv.data.length; i-- > 0;)
			if (!Maths.isFinite(A_inv.data[i]))
				return false;

		for (int i = 0, index = 0; i < solverSize; i++)
			for (int j = 0; j < solverSize; j++)
			{
				a[i][j] = A_inv.data[index++];
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

			// Reset size dependent objects
			choleskySolver = null;
			linearSolver = null;
			X = null;
			A_inv = null;
		}
	}

	private DenseMatrix64F getX()
	{
		if (X == null)
			X = new DenseMatrix64F(solverSize, 1);
		return X;
	}

	private DenseMatrix64F getA_inv()
	{
		if (A_inv == null)
			A_inv = new DenseMatrix64F(solverSize, solverSize);
		return A_inv;
	}

	private LinearSolver<DenseMatrix64F> getCholeskyLDLTSolver()
	{
		// This is a Cholesky LDLT solver that should be faster than the Cholesky solver.
		// It only works on symmetric positive definite matrices.
		if (choleskyLDLTSolver == null)
			choleskyLDLTSolver = new LinearSolverCholLDL();
		return choleskyLDLTSolver;
	}

	private LinearSolver<DenseMatrix64F> getCholeskySolver()
	{
		if (choleskySolver == null)
			// This is a Cholesky solver that only works on symmetric positive definite matrices
			choleskySolver = LinearSolverFactory.symmPosDef(solverSize);
		return choleskySolver;
	}

	private LinearSolver<DenseMatrix64F> getLinearSolver()
	{
		if (linearSolver == null)
			// This should work on any matrix
			linearSolver = LinearSolverFactory.linear(solverSize);
		return linearSolver;
	}

	private LinearSolver<DenseMatrix64F> getPseudoInverseSolver()
	{
		// The pseudo inverse is constructed using the non-singular sub matrix of A

		if (pseudoInverseSolver == null)
			pseudoInverseSolver = LinearSolverFactory.pseudoInverse(false);
		return pseudoInverseSolver;
	}

	private LinearSolver<DenseMatrix64F> getInversionSolver()
	{
		if (inversionSolver == null)
			// This should work on any matrix
			inversionSolver = new InversionSolver();
		return inversionSolver;
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
	 * Note: If the matrix is singular then a pseudo inverse will be computed.
	 *
	 * @param a
	 *            the matrix a
	 * @return False if there is no solution
	 */
	public boolean invertSymmPosDef(double[][] a)
	{
		DenseMatrix64F A = new DenseMatrix64F(a);

		createSolver(A.numCols);
		DenseMatrix64F A2 = (errorChecking) ? A.copy() : null;

		LinearSolver<DenseMatrix64F> primarySolver = (A.numCols <= UnrolledInverseFromMinor.MAX) ? getInversionSolver()
				: getCholeskyLDLTSolver();

		if (invert(primarySolver, A, A2, a, false))
			return true;
		if (primarySolver.modifiesA())
			A = (errorChecking) ? A2.copy() : new DenseMatrix64F(a);
		return invert(getPseudoInverseSolver(), A, A2, a, true);
	}

	private boolean invert(LinearSolver<DenseMatrix64F> solver, DenseMatrix64F A, DenseMatrix64F A2, double[][] a,
			boolean pseudoInverse)
	{
		if (!initialiseSolver(solver, A))
			return false;
		if (!invert(a))
			return false;
		if (A2 != null)
			return validateInversion(A2, pseudoInverse);
		return true;
	}

	private boolean validateInversion(DenseMatrix64F A, boolean pseudoInverse)
	{
		// Check for the identity matrix:
		// Compute A A_inv = I
		final int n = A.numCols;
		DenseMatrix64F I = new DenseMatrix64F(n, n);
		CommonOps.mult(A, A_inv, I);

		if (pseudoInverse)
		{
			for (int i = n, index = I.data.length; i-- > 0;)
			{
				for (int j = n; j-- > 0;)
				{
					if (j == i)
					{
						--index;
						// If using the pseudo inverse then the diagonal can be zero or 1
						if (invalid(I.data[index], 1) && invalid(I.data[index], 0))
							return false;
					}
					else
					{
						if (invalid(I.data[--index], 0))
							return false;
					}
				}
			}
		}
		else
		{
			for (int i = n, index = I.data.length; i-- > 0;)
			{
				for (int j = n; j-- > 0;)
				{
					if (invalid(I.data[--index], (j == i) ? 1 : 0))
						return false;
				}
			}
		}
		return true;
	}

	private boolean invalid(double e, double o)
	{
		if (equal.almostEqualComplement(e, o))
			return false;
		//System.out.printf("Bad solution: %g != %g (%g = %d)\n", e, o, DoubleEquality.relativeError(e, o),
		//		DoubleEquality.complement(e, o));
		return true;
	}

	/**
	 * Computes the inverse of the symmetric positive definite matrix and returns only the diagonal.
	 * <p>
	 * Note: If the matrix is singular then a pseudo inverse will be computed.
	 *
	 * @param a
	 *            the matrix a
	 * @return The diagonal of the inverted matrix (or null)
	 */
	public double[] invertSymmPosDefDiagonal(double[][] a)
	{
		DenseMatrix64F A = new DenseMatrix64F(a);

		// Try a fast inversion of the diagonal
		if (A.numCols <= UnrolledInverseFromMinorExt.MAX)
		{
			double[] d = UnrolledInverseFromMinorExt.inv(A);
			if (d != null)
				return d;
		}

		createSolver(A.numCols);
		DenseMatrix64F A2 = (errorChecking) ? A.copy() : null;

		if (A.numCols <= UnrolledInverseFromMinorExt.MAX)
		{
			// The fast inversion failed so try a pseudo-inverse
			if (!invert(getPseudoInverseSolver(), A, A2, a, true))
				return null;
		}
		else
		{
			// The matrix was too big for fast inversion so try linear algebra
			if (!invert(getCholeskyLDLTSolver(), A, A2, a, false))
			{
				if (choleskyLDLTSolver.modifiesA())
					A = (errorChecking) ? A2.copy() : new DenseMatrix64F(a);
				if (!invert(getPseudoInverseSolver(), A, A2, a, true))
					return null;
			}
		}

		// We reach here when 'a' has been inverted
		double[] d = new double[A.numCols];
		for (int i = 0; i < d.length; i++)
			d[i] = a[i][i];
		return d;
	}

	/**
	 * Convert a dense matrix to a row/column format
	 *
	 * @param A
	 *            the matrix
	 * @return the row/column format
	 */
	public static double[][] toData(DenseMatrix64F A)
	{
		final int numRows = A.numRows, numCols = A.numCols;
		final double[][] out = new double[numRows][];
		for (int i = 0, pos = 0; i < numRows; i++, pos += numRows)
		{
			out[i] = new double[numCols];
			System.arraycopy(A.data, pos, out[i], 0, numCols);
		}
		return out;
	}

	/**
	 * Convert a dense matrix to a row/column format
	 *
	 * @param A
	 *            the matrix
	 * @param out
	 *            the row/column format
	 */
	private static void toSquareData(DenseMatrix64F A, double[][] out)
	{
		final int numRows = A.numRows;
		for (int i = 0, pos = 0; i < numRows; i++, pos += numRows)
			System.arraycopy(A.data, pos, out[i], 0, numRows);
	}
}
