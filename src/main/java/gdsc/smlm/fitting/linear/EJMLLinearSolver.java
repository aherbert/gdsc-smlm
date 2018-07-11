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

import org.ejml.alg.dense.linsol.chol.LinearSolverCholLDL;
import org.ejml.alg.dense.misc.UnrolledInverseFromMinor;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.LinearSolver;
import org.ejml.factory.LinearSolverFactory;
import org.ejml.ops.CommonOps;

import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.Maths;

/**
 * Solves (one) linear equation, a x = b.
 * <p>
 * Wraps the LinearSolver class from the EJML (Efficient Java Matrix Library) library.
 * <p>
 * This class assumes that the input matrix is symmetric positive definite matrices, for instance it has been generated
 * in the numerical solution of partial differential equations (for example a Hessian matrix).
 * <p>
 * The solve and invert methods uses two solvers, the first is an optimised solver chosen for speed. If the solver fails
 * then a pseudo inverse solver is used so that a solution can be found. For example the first solver may fail if zeros
 * occur on the diagonal of the matrix, for example if no gradient is available for a parameter.
 */
public class EJMLLinearSolver
{
	/**
	 * Solve the matrix using direct inversion
	 */
	private class InversionSolver implements LinearSolver<DenseMatrix64F>
	{
		private DenseMatrix64F A;

		@Override
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

		@Override
		public double quality()
		{
			return 0;
		}

		@Override
		public void solve(DenseMatrix64F B, DenseMatrix64F X)
		{
			CommonOps.mult(A, B, X);
		}

		@Override
		public void invert(DenseMatrix64F A_inv)
		{
			System.arraycopy(A.data, 0, A_inv.data, 0, A.data.length);
		}

		@Override
		public boolean modifiesA()
		{
			return true;
		}

		@Override
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
	private double inversionTolerance = 0;

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
	 * Instantiates a new EJML linear solver with tolerance for the linear solution.
	 *
	 * @param equal
	 *            the object for equality
	 * @param inversionTolerance
	 *            the inversion tolerance
	 */
	public EJMLLinearSolver(DoubleEquality equal, double inversionTolerance)
	{
		setEqual(equal);
		setInversionTolerance(inversionTolerance);
	}

	/**
	 * Instantiates a new EJML linear solver with tolerance for the linear solution.
	 *
	 * @param maxRelativeError
	 *            the max relative error for equality
	 * @param maxAbsoluteError
	 *            the max absolute error for equality
	 */
	public EJMLLinearSolver(double maxRelativeError, double maxAbsoluteError)
	{
		this(new DoubleEquality(maxRelativeError, maxAbsoluteError));
	}

	/**
	 * Instantiates a new EJML linear solver with tolerance for the linear solution.
	 *
	 * @param maxRelativeError
	 *            the max relative error for equality
	 * @param maxAbsoluteError
	 *            the max absolute error for equality
	 * @param inversionTolerance
	 *            the inversion tolerance
	 */
	public EJMLLinearSolver(double maxRelativeError, double maxAbsoluteError, double inversionTolerance)
	{
		this(new DoubleEquality(maxRelativeError, maxAbsoluteError), inversionTolerance);
	}

	/**
	 * Instantiates a new EJML linear solver with tolerance for the linear solution
	 *
	 * @param significantDigits
	 *            the significant digits for equality
	 * @param maxAbsoluteError
	 *            the max absolute error for equality
	 * @deprecated
	 */
	@Deprecated
	public EJMLLinearSolver(int significantDigits, double maxAbsoluteError)
	{
		this(new DoubleEquality(significantDigits, maxAbsoluteError));
	}

	/**
	 * Creates the solver for inversion.
	 *
	 * @param inversionTolerance
	 *            the inversion tolerance
	 * @return the EJML linear solver
	 */
	public static EJMLLinearSolver createForInversion(double inversionTolerance)
	{
		final EJMLLinearSolver s = new EJMLLinearSolver();
		s.setInversionTolerance(inversionTolerance);
		return s;
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On output b replaced by x. Matrix a may be modified.
	 *
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solveLinear(DenseMatrix64F A, DenseMatrix64F B)
	{
		createSolver(A.numCols);
		return solve(getLinearSolver(), A, B);
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On output b replaced by x. Matrix a may be modified.
	 *
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solveCholesky(DenseMatrix64F A, DenseMatrix64F B)
	{
		createSolver(A.numCols);
		return solve(getCholeskySolver(), A, B);
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On output b replaced by x. Matrix a may be modified.
	 *
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solveCholeskyLDLT(DenseMatrix64F A, DenseMatrix64F B)
	{
		createSolver(A.numCols);
		return solve(getCholeskyLDLTSolver(), A, B);
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On output b replaced by x. Matrix a may be modified.
	 *
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solvePseudoInverse(DenseMatrix64F A, DenseMatrix64F B)
	{
		createSolver(A.numCols);
		return solve(getPseudoInverseSolver(), A, B);
	}

	/**
	 * Solves (one) linear equation, a x = b by direct inversion. Works on small matrices up to size 5.
	 * <p>
	 * On output b replaced by x. Matrix a may be modified.
	 *
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solveDirectInversion(DenseMatrix64F A, DenseMatrix64F B)
	{
		createSolver(A.numCols);
		return solve(getInversionSolver(), A, B);
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On output b replaced by x. Matrix a may be modified.
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
	 * Output is written to x.
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
			// We need A+B for error checking so solve without modification
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
	 * Output is written to x.
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
	 * Output is written to x.
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
			if (!equal.almostEqualRelativeOrAbsolute(b.data[i], bi))
			{
				//System.out.printf("Bad solution: %g != %g (%g = %d)\n", b.data[i], bi,
				//		DoubleEquality.relativeError(b.data[i], bi), DoubleEquality.complement(b.data[i], bi));
				return false;
			}
		}
		//System.out.println("OK");
		return true;
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On output b replaced by x. Matrix a may be modified.
	 * <p>
	 * Solve using the CholeskyLDLT method or, if that fails (due to a singular matrix), the PseudoInverse
	 * decomposition.
	 *
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solve(DenseMatrix64F A, DenseMatrix64F B)
	{
		createSolver(A.numCols);

		// Speed tests show the Cholesky solver marginally outperforms the
		// CholeskyLDLT solver as the size increases. Note: The EJML factory
		// returns a Cholesky solver for symmetric positive definite matrices
		// so we use this.

		// Note: Use solveSafe as we need the A and B matrix for the subsequent
		// solve attempt if failure

		if (solveSafe(getCholeskySolver(), A, B, getX()))
		{
			if (!errorChecking || validate(A, X, B))
			{
				System.arraycopy(X.data, 0, B.data, 0, A.numCols);
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
			System.arraycopy(X.data, 0, B.data, 0, A.numCols);
			return true;
		}

		return false;
	}

	/**
	 * Checks if a solve may modify A.
	 *
	 * @return true, if a solve may modify A
	 * @see #solve(DenseMatrix64F, DenseMatrix64F)
	 * @see #solve(double[], double[])
	 */
	public boolean solveModifiesA()
	{
		if (errorChecking)
			return false;
		return getPseudoInverseSolver().modifiesA();
	}

	/**
	 * Computes the inverse of the 'A' matrix passed into the last successful solve method.
	 * <p>
	 * On output a[n][n] replaced by the inverse of the solved matrix a.
	 *
	 * @param A
	 *            the matrix a
	 *
	 * @return False if the last solve attempt failed, or inversion produces non finite values
	 */
	public boolean invertLastA(DenseMatrix64F A)
	{
		if (lastSuccessfulSolver == null)
			return false;

		lastSuccessfulSolver.invert(getA_inv());

		// Check for NaN or Infinity
		double[] a_inv = A_inv.data;
		for (int i = a_inv.length; i-- > 0;)
			if (!Maths.isFinite(a_inv[i]))
				return false;

		// Q. Should we check the product is the identity matrix?
		// This will require that we have the original matrix A used to initialise the solver.

		System.arraycopy(a_inv, 0, A.data, 0, a_inv.length);

		return true;
	}

	private boolean initialiseSolver(LinearSolver<DenseMatrix64F> solver, DenseMatrix64F A)
	{
		lastSuccessfulSolver = null;
		if (!solver.setA(A))
			return false;
		lastSuccessfulSolver = solver;
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

	/**
	 * Gets the working space to store the output solution x.
	 *
	 * @return the x
	 */
	private DenseMatrix64F getX()
	{
		if (X == null)
			X = new DenseMatrix64F(solverSize, 1);
		return X;
	}

	/**
	 * Gets the working space to store the inverse of A.
	 *
	 * @return the space for A^-1
	 */
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
	 * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
	 *
	 * @return False if the matrix is singular (no solution)
	 */
	public boolean invertLinear(DenseMatrix64F A)
	{
		createSolver(A.numCols);
		return invertUnsafe(getLinearSolver(), A, false);
	}

	/**
	 * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
	 *
	 * @return False if the matrix is singular (no solution)
	 */
	public boolean invertCholesky(DenseMatrix64F A)
	{
		createSolver(A.numCols);
		return invertUnsafe(getCholeskySolver(), A, false);
	}

	/**
	 * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
	 *
	 * @return False if the matrix is singular (no solution)
	 */
	public boolean invertCholeskyLDLT(DenseMatrix64F A)
	{
		createSolver(A.numCols);
		return invertUnsafe(getCholeskyLDLTSolver(), A, false);
	}

	/**
	 * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
	 *
	 * @return False if the matrix is singular (no solution)
	 */
	public boolean invertPseudoInverse(DenseMatrix64F A)
	{
		createSolver(A.numCols);
		return invertUnsafe(getPseudoInverseSolver(), A, true);
	}

	/**
	 * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
	 *
	 * @return False if the matrix is singular (no solution)
	 */
	public boolean invertDirectInversion(DenseMatrix64F A)
	{
		createSolver(A.numCols);
		return invertUnsafe(getInversionSolver(), A, false);
	}

	/**
	 * Invert symmetric positive definite matrix A and returns only the diagonal.
	 * <p>
	 * Only supports matrix size up to 5.
	 *
	 * @return The diagonal of A^-1 (or null if the matrix is singular (no solution) or too large)
	 */
	public double[] invertDiagonalDirectInversion(DenseMatrix64F A)
	{
		return (A.numCols <= UnrolledInverseFromMinorExt.MAX) ? UnrolledInverseFromMinorExt.inv(A) : null;
	}

	/**
	 * Computes the inverse of the symmetric positive definite matrix. On output a is replaced by the inverse of a.
	 * <p>
	 * Note: If the matrix is singular then a pseudo inverse will be computed.
	 *
	 * @param A
	 *            the matrix a
	 * @return False if there is no solution
	 */
	public boolean invert(DenseMatrix64F A)
	{
		createSolver(A.numCols);

		// Speed tests show the Cholesky solver marginally outperforms the
		// CholeskyLDLT solver as the size increases.
		// The DirectInversion solver is faster when the size < 5.
		// Note: The EJML factory returns a Cholesky solver for symmetric
		// positive definite matrices so we use this in preference to a CholeskyLDLT.

		LinearSolver<DenseMatrix64F> primarySolver = (A.numCols < 5) ? getInversionSolver() : getCholeskySolver();
		if (invertSafe(primarySolver, A, false))
			return true;

		return invertUnsafe(getPseudoInverseSolver(), A, true);
	}

	/**
	 * Checks if an inversion may modify A.
	 *
	 * @return true, if an inversion may modify A
	 * @see #invert(DenseMatrix64F)
	 * @see #invert(double[], int)
	 */
	public boolean invertModifiesA()
	{
		if (isInversionTolerance())
			return false;
		return getPseudoInverseSolver().modifiesA();
	}

	private boolean invertSafe(LinearSolver<DenseMatrix64F> solver, DenseMatrix64F A, boolean pseudoInverse)
	{
		DenseMatrix64F Ain = (solver.modifiesA() || isInversionTolerance()) ? A.copy() : A;
		if (!initialiseSolver(solver, Ain))
			return false;
		solver.invert(getA_inv());

		// Check for NaN or Infinity
		double[] a_inv = A_inv.data;
		for (int i = a_inv.length; i-- > 0;)
			if (!Maths.isFinite(a_inv[i]))
				return false;

		if (isInversionTolerance() && invalidInversion(A, pseudoInverse))
			return false;

		System.arraycopy(a_inv, 0, A.data, 0, a_inv.length);

		return true;
	}

	private boolean invertUnsafe(LinearSolver<DenseMatrix64F> solver, DenseMatrix64F A, boolean pseudoInverse)
	{
		DenseMatrix64F Ain = (isInversionTolerance()) ? A.copy() : A;
		if (!initialiseSolver(solver, Ain))
			return false;
		solver.invert(getA_inv());

		// Check for NaN or Infinity
		double[] a_inv = A_inv.data;
		for (int i = a_inv.length; i-- > 0;)
			if (!Maths.isFinite(a_inv[i]))
			{
				System.out.printf("Not finite\n");
				return false;
			}

		if (isInversionTolerance() && invalidInversion(A, pseudoInverse))
			return false;

		System.arraycopy(a_inv, 0, A.data, 0, a_inv.length);

		return true;
	}

	private boolean invalidInversion(DenseMatrix64F A, boolean pseudoInverse)
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
							return true;
					}
					else
					{
						if (invalid(I.data[--index], 0))
							return true;
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
						return true;
				}
			}
		}
		return false;
	}

	private boolean invalid(double e, double o)
	{
		//if (Math.abs(e - o) > inversionTolerance)
		//	System.out.printf("Bad solution: %g != %g (%g = %d)\n", e, o, DoubleEquality.relativeError(e, o),
		//			DoubleEquality.complement(e, o));
		return (Math.abs(e - o) > inversionTolerance);
	}

	/**
	 * Computes the inverse of the symmetric positive definite matrix and returns only the diagonal. Will not modify the
	 * matrix.
	 * <p>
	 * Note: If the matrix is singular then a pseudo inverse will be computed.
	 *
	 * @param A
	 *            the matrix a
	 * @return The diagonal of the inverted matrix (or null)
	 */
	public double[] invertDiagonal(DenseMatrix64F A)
	{
		// Try a fast inversion of the diagonal
		if (A.numCols <= UnrolledInverseFromMinorExt.MAX)
		{
			double[] d = UnrolledInverseFromMinorExt.inv(A);
			if (d != null)
				return d;
		}

		createSolver(A.numCols);

		if (A.numCols <= UnrolledInverseFromMinorExt.MAX)
		{
			// The fast inversion failed so try a pseudo-inverse
			if (!invertSafe(getPseudoInverseSolver(), A, true))
				return null;
		}
		else
		{
			// The matrix was too big for fast inversion so try linear algebra
			if (!invertSafe(getCholeskyLDLTSolver(), A, false))
			{
				if (!invertSafe(getPseudoInverseSolver(), A, true))
					return null;
			}
		}

		// We reach here when 'a' has been inverted
		double[] d = new double[A.numCols];
		for (int i = 0, j = 0; i < d.length; i++, j += A.numCols + 1)
			d[i] = A.get(j);
		return d;
	}

	/**
	 * Computes the inverse of the symmetric positive definite matrix and returns only the diagonal. May modify the
	 * matrix.
	 * <p>
	 * Note: If the matrix is singular then a pseudo inverse will be computed.
	 *
	 * @param A
	 *            the matrix a
	 * @return The diagonal of the inverted matrix (or null)
	 */
	private double[] invertDiagonalUnsafe(DenseMatrix64F A)
	{
		// Try a fast inversion of the diagonal
		if (A.numCols <= UnrolledInverseFromMinorExt.MAX)
		{
			double[] d = UnrolledInverseFromMinorExt.inv(A);
			if (d != null)
				return d;
		}

		createSolver(A.numCols);

		if (A.numCols <= UnrolledInverseFromMinorExt.MAX)
		{
			// The fast inversion failed so try a pseudo-inverse
			if (!invertUnsafe(getPseudoInverseSolver(), A, true))
				return null;
		}
		else
		{
			// The matrix was too big for fast inversion so try linear algebra
			if (!invertSafe(getCholeskyLDLTSolver(), A, false))
			{
				if (!invertUnsafe(getPseudoInverseSolver(), A, true))
					return null;
			}
		}

		// We reach here when 'a' has been inverted
		double[] d = new double[A.numCols];
		for (int i = 0, j = 0; i < d.length; i++, j += A.numCols + 1)
			d[i] = A.get(j);
		return d;
	}

	/**
	 * Convert a dense matrix to a row/column format
	 *
	 * @param A
	 *            the matrix
	 * @return the row/column format
	 */
	public static double[][] toSquareData(DenseMatrix64F A)
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
	 * Convert a dense matrix to a row/column format. The output arrays must be the correct size.
	 *
	 * @param A
	 *            the matrix
	 * @param out
	 *            the row/column format
	 */
	public static void toSquareData(DenseMatrix64F A, double[][] out)
	{
		final int numRows = A.numRows;
		for (int i = 0, pos = 0; i < numRows; i++, pos += numRows)
			System.arraycopy(A.data, pos, out[i], 0, numRows);
	}

	/**
	 * Create a new DenseMatrix from the input matrix a. Modifications to the matrix are not passed through to the input
	 * array! The matrix can be converted back using {@link #toSquareData(DenseMatrix64F, double[][])}.
	 * <p>
	 * This is provided as a bridge method between the functions that accept primitive arrays and those that accept
	 * DenseMatrix.
	 *
	 * @param a
	 *            the matrix
	 * @return the dense matrix
	 */
	public static DenseMatrix64F toA(double[][] a)
	{
		return new DenseMatrix64F(a);
	}

	/**
	 * Create a new DenseMatrix from the input matrix a. Modifications to the matrix are passed through to the input
	 * array!
	 * <p>
	 * This is provided as a bridge method between the functions that accept primitive arrays and those that accept
	 * DenseMatrix.
	 *
	 * @param a
	 *            the matrix
	 * @param n
	 *            the number of columns/rows
	 * @return the dense matrix
	 */
	public static DenseMatrix64F toA(double[] a, int n)
	{
		return DenseMatrix64F.wrap(n, n, a);
	}

	/**
	 * Wrap the input array b in a DenseMatrix. Modifications to the matrix are passed through to the input array.
	 * <p>
	 * This is provided as a bridge method between the functions that accept primitive arrays and those that accept
	 * DenseMatrix.
	 *
	 * @param b
	 *            the array
	 * @return the dense matrix
	 */
	public static DenseMatrix64F toB(double[] b)
	{
		return DenseMatrix64F.wrap(b.length, 1, b);
	}

	// Methods for input of primitive arrays

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On output b replaced by x.
	 *
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solveLinear(double[][] a, double[] b)
	{
		return solveLinear(toA(a), toB(b));
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On output b replaced by x.
	 *
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solveCholesky(double[][] a, double[] b)
	{
		return solveCholesky(toA(a), toB(b));
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On output b replaced by x.
	 *
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solveCholeskyLDLT(double[][] a, double[] b)
	{
		return solveCholeskyLDLT(toA(a), toB(b));
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On output b replaced by x.
	 *
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solvePseudoInverse(double[][] a, double[] b)
	{
		return solvePseudoInverse(toA(a), toB(b));
	}

	/**
	 * Solves (one) linear equation, a x = b by direct inversion. Works on small matrices up to size 5.
	 * <p>
	 * On output b replaced by x.
	 *
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solveDirectInversion(double[][] a, double[] b)
	{
		return solveDirectInversion(toA(a), toB(b));
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On output b replaced by x.
	 * <p>
	 * Solve using the CholeskyLDLT method or, if that fails (due to a singular matrix), the PseudoInversion
	 * decomposition.
	 *
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solve(double[][] a, double[] b)
	{
		return solve(toA(a), toB(b));
	}

	/**
	 * Computes the inverse of the 'A' matrix passed into the last successful solve method.
	 * <p>
	 * On output a[n][n] replaced by the inverse of the solved matrix a. If any column/row index was removed (as it was
	 * set to zero in the input matrix) then the resulting column/row index will be set to zero.
	 *
	 * @param a
	 *            the matrix a
	 * @return False if the last solve attempt failed, or inversion produces non finite values
	 */
	public boolean invertLastA(double[][] a)
	{
		if (lastSuccessfulSolver == null)
			return false;

		lastSuccessfulSolver.invert(getA_inv());

		// Check for NaN or Infinity
		double[] a_inv = A_inv.data;
		for (int i = a_inv.length; i-- > 0;)
			if (!Maths.isFinite(a_inv[i]))
				return false;

		// Q. Should we check the product is the identity matrix?
		// This will require that we have the original matrix A used to initialise the solver.

		toSquareData(A_inv, a);
		return true;
	}

	/**
	 * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
	 *
	 * @param a
	 *            the matrix a
	 * @return False if the matrix is singular (no solution)
	 */
	public boolean invertLinear(double[][] a)
	{
		DenseMatrix64F A = toA(a);
		if (!invertLinear(A))
			return false;
		toSquareData(A, a);
		return true;
	}

	/**
	 * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
	 *
	 * @param a
	 *            the matrix a
	 * @return False if the matrix is singular (no solution)
	 */
	public boolean invertCholesky(double[][] a)
	{
		DenseMatrix64F A = toA(a);
		if (!invertCholesky(A))
			return false;
		toSquareData(A, a);
		return true;
	}

	/**
	 * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
	 *
	 * @param a
	 *            the matrix a
	 * @return False if the matrix is singular (no solution)
	 */
	public boolean invertCholeskyLDLT(double[][] a)
	{
		DenseMatrix64F A = toA(a);
		if (!invertCholeskyLDLT(A))
			return false;
		toSquareData(A, a);
		return true;
	}

	/**
	 * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
	 *
	 * @param a
	 *            the matrix a
	 * @return False if the matrix is singular (no solution)
	 */
	public boolean invertPseudoInverse(double[][] a)
	{
		DenseMatrix64F A = toA(a);
		if (!invertPseudoInverse(A))
			return false;
		toSquareData(A, a);
		return true;
	}

	/**
	 * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
	 *
	 * @param a
	 *            the matrix a
	 * @return False if the matrix is singular (no solution)
	 */
	public boolean invertDirectInversion(double[][] a)
	{
		DenseMatrix64F A = toA(a);
		if (!invertDirectInversion(A))
			return false;
		toSquareData(A, a);
		return true;
	}

	/**
	 * Invert symmetric positive definite matrix A and returns only the diagonal.
	 * <p>
	 * Only supports matrix size up to 5.
	 *
	 * @param a
	 *            the matrix a
	 * @return The diagonal of A^-1 (or null if the matrix is singular (no solution) or too large)
	 */
	public double[] invertDiagonalDirectInversion(double[][] a)
	{
		if (a.length > UnrolledInverseFromMinorExt.MAX)
			return null;
		return UnrolledInverseFromMinorExt.inv(toA(a));
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
	public boolean invert(double[][] a)
	{
		DenseMatrix64F A = toA(a);
		if (!invert(A))
			return false;
		toSquareData(A, a);
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
	public double[] invertDiagonal(double[][] a)
	{
		// Use the unsafe method as the matrix has been converted so can be modified
		return invertDiagonalUnsafe(new DenseMatrix64F(a));
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On output b replaced by x. Matrix a may be modified.
	 *
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solveLinear(double[] a, double[] b)
	{
		return solveLinear(toA(a, b.length), toB(b));
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On output b replaced by x. Matrix a may be modified.
	 *
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solveCholesky(double[] a, double[] b)
	{
		return solveCholesky(toA(a, b.length), toB(b));
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On output b replaced by x. Matrix a may be modified.
	 *
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solveCholeskyLDLT(double[] a, double[] b)
	{
		return solveCholeskyLDLT(toA(a, b.length), toB(b));
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On output b replaced by x. Matrix a may be modified.
	 *
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solvePseudoInverse(double[] a, double[] b)
	{
		return solvePseudoInverse(toA(a, b.length), toB(b));
	}

	/**
	 * Solves (one) linear equation, a x = b by direct inversion. Works on small matrices up to size 5.
	 * <p>
	 * On output b replaced by x. Matrix a may be modified.
	 *
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solveDirectInversion(double[] a, double[] b)
	{
		return solveDirectInversion(toA(a, b.length), toB(b));
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On output b replaced by x. Matrix a may be modified.
	 * <p>
	 * Solve using the CholeskyLDLT method or, if that fails (due to a singular matrix), the PseudoInversion
	 * decomposition.
	 *
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solve(double[] a, double[] b)
	{
		return solve(toA(a, b.length), toB(b));
	}

	/**
	 * Computes the inverse of the 'A' matrix passed into the last successful solve method.
	 * <p>
	 * On output a[n][n] replaced by the inverse of the solved matrix a. If any column/row index was removed (as it was
	 * set to zero in the input matrix) then the resulting column/row index will be set to zero.
	 *
	 * @param a
	 *            the matrix a
	 * @return False if the last solve attempt failed, or inversion produces non finite values
	 */
	public boolean invertLastA(double[] a)
	{
		if (lastSuccessfulSolver == null)
			return false;

		lastSuccessfulSolver.invert(getA_inv());

		// Check for NaN or Infinity
		double[] a_inv = A_inv.data;
		for (int i = a_inv.length; i-- > 0;)
			if (!Maths.isFinite(a_inv[i]))
				return false;

		// Q. Should we check the product is the identity matrix?
		// This will require that we have the original matrix A used to initialise the solver.

		System.arraycopy(a_inv, 0, a, 0, a_inv.length);

		return true;
	}

	/**
	 * Computes the inverse of the symmetric positive definite matrix. On output a is replaced by the inverse of a.
	 * <p>
	 * Note: If the matrix is singular then a pseudo inverse will be computed.
	 *
	 * @param a
	 *            the matrix a
	 * @param n
	 *            the number of columns/rows
	 * @return False if there is no solution
	 */
	public boolean invert(double[] a, int n)
	{
		return invert(toA(a, n));
	}

	/**
	 * Computes the inverse of the symmetric positive definite matrix and returns only the diagonal.
	 * <p>
	 * Note: If the matrix is singular then a pseudo inverse will be computed.
	 *
	 * @param a
	 *            the matrix a
	 * @param n
	 *            the number of columns/rows
	 * @return The diagonal of the inverted matrix (or null)
	 */
	public double[] invertDiagonal(double[] a, int n)
	{
		return invertDiagonal(toA(a, n));
	}

	/**
	 * Gets the inversion tolerance. Inversions are checked by ensuring that the product matches the identity matrix: A
	 * * A^-1 = I. Elements must be within the tolerance or else the inversion is rejected. Set to zero to disable.
	 *
	 * @return the inversion tolerance
	 */
	public double getInversionTolerance()
	{
		return inversionTolerance;
	}

	/**
	 * Sets the inversion tolerance. Inversions are checked by ensuring that the product matches the identity matrix: A
	 * * A^-1 = I. Elements must be within the tolerance or else the inversion is rejected. Set to zero to disable.
	 *
	 * @param inversionTolerance
	 *            the new inversion tolerance
	 */
	public void setInversionTolerance(double inversionTolerance)
	{
		this.inversionTolerance = inversionTolerance;
	}

	/**
	 * Checks if there is an inversion tolerance.
	 *
	 * @return true, if there is an inversion tolerance
	 */
	private boolean isInversionTolerance()
	{
		return inversionTolerance > 0;
	}
}
