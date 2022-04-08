/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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

import java.util.logging.Logger;
import org.ejml.alg.dense.linsol.chol.LinearSolverCholLDL;
import org.ejml.alg.dense.misc.UnrolledInverseFromMinor;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.LinearSolver;
import org.ejml.factory.LinearSolverFactory;
import org.ejml.ops.CommonOps;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;

/**
 * Solves (one) linear equation, A x = b.
 *
 * <p>Wraps the LinearSolver class from the EJML (Efficient Java Matrix Library) library.
 *
 * <p>This class assumes that the input matrix is symmetric positive definite matrices, for instance
 * it has been generated in the numerical solution of partial differential equations (for example a
 * Hessian matrix).
 *
 * <p>The solve and invert methods uses two solvers, the first is an optimised solver chosen for
 * speed. If the solver fails then a pseudo inverse solver is used so that a solution can be found.
 * For example the first solver may fail if zeros occur on the diagonal of the matrix, for example
 * if no gradient is available for a parameter.
 */
public class EjmlLinearSolver {
  // Allow matrix name A and vector names b & x for equation A x = b
  // CHECKSTYLE.OFF: MemberName
  // CHECKSTYLE.OFF: ParameterName

  /**
   * Solve the matrix using direct inversion.
   */
  private static class InversionSolver implements LinearSolver<DenseMatrix64F> {
    /** The matrix A. */
    private DenseMatrix64F A;

    /**
     * Instantiates a new inversion solver.
     */
    InversionSolver() {}

    @Override
    public boolean setA(DenseMatrix64F A) {
      if (A.numCols <= UnrolledInverseFromMinor.MAX) {
        // Direct inversion using the determinant
        if (A.numCols >= 2) {
          UnrolledInverseFromMinor.inv(A, A);
        } else {
          A.set(0, 1.0 / A.get(0));
        }

        // Check for NaN or Infinity
        for (int i = A.data.length; i-- > 0;) {
          if (!Double.isFinite(A.data[i])) {
            return false;
          }
        }

        this.A = A;
        return true;
      }
      return false;
    }

    @Override
    public double quality() {
      return 0;
    }

    @Override
    public void solve(DenseMatrix64F b, DenseMatrix64F x) {
      CommonOps.mult(A, b, x);
    }

    @Override
    public void invert(DenseMatrix64F Ainv) {
      System.arraycopy(A.data, 0, Ainv.data, 0, A.data.length);
    }

    @Override
    public boolean modifiesA() {
      return true;
    }

    @Override
    public boolean modifiesB() {
      return false;
    }
  }

  /** The linear solver. */
  private LinearSolver<DenseMatrix64F> linearSolver;

  /** The pseudo inverse solver. */
  private LinearSolver<DenseMatrix64F> pseudoInverseSolver;

  /** The cholesky solver. */
  private LinearSolver<DenseMatrix64F> choleskySolver;

  /** The cholesky LDLT solver. */
  private LinearSolver<DenseMatrix64F> choleskyLdlTSolver;

  /** The inversion solver. */
  private LinearSolver<DenseMatrix64F> inversionSolver;

  /** The last successful solver. */
  private LinearSolver<DenseMatrix64F> lastSuccessfulSolver;

  /** The vector x. */
  private DenseMatrix64F x;

  /** The inverse of matrix A. */
  private DenseMatrix64F invA;

  /** The solver size. */
  private int solverSize;

  /** The error checking flag. Set to true to check the solution x to linear equations A x = b. */
  private boolean errorChecking;

  /** The class used to check equality (with zero). */
  private DoubleEquality equal;

  /** The inversion tolerance. */
  private double inversionTolerance;

  /**
   * Instantiates a new EJML linear solver.
   */
  public EjmlLinearSolver() {}

  /**
   * Instantiates a new EJML linear solver with tolerance for the linear solution.
   *
   * @param equal the object for equality
   */
  public EjmlLinearSolver(DoubleEquality equal) {
    this(equal, 0);
  }

  /**
   * Instantiates a new EJML linear solver with tolerance for the linear solution.
   *
   * @param equal the object for equality
   * @param inversionTolerance the inversion tolerance
   */
  public EjmlLinearSolver(DoubleEquality equal, double inversionTolerance) {
    this.equal = equal;
    errorChecking = (equal != null);
    this.inversionTolerance = inversionTolerance;
  }

  /**
   * Instantiates a new EJML linear solver with tolerance for the linear solution.
   *
   * @param maxRelativeError the max relative error for equality
   * @param maxAbsoluteError the max absolute error for equality
   */
  public EjmlLinearSolver(double maxRelativeError, double maxAbsoluteError) {
    this(new DoubleEquality(maxRelativeError, maxAbsoluteError));
  }

  /**
   * Instantiates a new EJML linear solver with tolerance for the linear solution.
   *
   * @param maxRelativeError the max relative error for equality
   * @param maxAbsoluteError the max absolute error for equality
   * @param inversionTolerance the inversion tolerance
   */
  public EjmlLinearSolver(double maxRelativeError, double maxAbsoluteError,
      double inversionTolerance) {
    this(new DoubleEquality(maxRelativeError, maxAbsoluteError), inversionTolerance);
  }

  /**
   * Creates the solver for inversion.
   *
   * @param inversionTolerance the inversion tolerance
   * @return the EJML linear solver
   */
  public static EjmlLinearSolver createForInversion(double inversionTolerance) {
    final EjmlLinearSolver s = new EjmlLinearSolver();
    s.setInversionTolerance(inversionTolerance);
    return s;
  }

  /**
   * Gets the inversion tolerance. Inversions are checked by ensuring that the product matches the
   * identity matrix: A * A^-1 = I. Elements must be within the tolerance or else the inversion is
   * rejected. Set to zero to disable.
   *
   * @return the inversion tolerance
   */
  public double getInversionTolerance() {
    return inversionTolerance;
  }

  /**
   * Sets the inversion tolerance. Inversions are checked by ensuring that the product matches the
   * identity matrix: A * A^-1 = I. Elements must be within the tolerance or else the inversion is
   * rejected. Set to zero to disable.
   *
   * @param inversionTolerance the new inversion tolerance
   */
  public void setInversionTolerance(double inversionTolerance) {
    this.inversionTolerance = inversionTolerance;
  }

  /**
   * Checks if there is an inversion tolerance.
   *
   * @return true, if there is an inversion tolerance
   */
  private boolean isInversionTolerance() {
    return inversionTolerance > 0;
  }

  /**
   * Solves (one) linear equation, A x = b.
   *
   * <p>On output b replaced by x. Matrix a may be modified.
   *
   * @param A the matrix A
   * @param b the vector b
   * @return False if the equation is singular (no solution)
   */
  public boolean solveLinear(DenseMatrix64F A, DenseMatrix64F b) {
    createSolver(A.numCols);
    return solveEquation(getLinearSolver(), A, b);
  }

  /**
   * Solves (one) linear equation, A x = b.
   *
   * <p>On output b replaced by x. Matrix a may be modified.
   *
   * @param A the matrix A
   * @param b the vector b
   * @return False if the equation is singular (no solution)
   */
  public boolean solveCholesky(DenseMatrix64F A, DenseMatrix64F b) {
    createSolver(A.numCols);
    return solveEquation(getCholeskySolver(), A, b);
  }

  /**
   * Solves (one) linear equation, A x = b.
   *
   * <p>On output b replaced by x. Matrix a may be modified.
   *
   * @param A the matrix A
   * @param b the vector b
   * @return False if the equation is singular (no solution)
   */
  public boolean solveCholeskyLdlT(DenseMatrix64F A, DenseMatrix64F b) {
    createSolver(A.numCols);
    return solveEquation(getCholeskyLdlTSolver(), A, b);
  }

  /**
   * Solves (one) linear equation, A x = b.
   *
   * <p>On output b replaced by x. Matrix a may be modified.
   *
   * @param A the matrix A
   * @param b the vector b
   * @return False if the equation is singular (no solution)
   */
  public boolean solvePseudoInverse(DenseMatrix64F A, DenseMatrix64F b) {
    createSolver(A.numCols);
    return solveEquation(getPseudoInverseSolver(), A, b);
  }

  /**
   * Solves (one) linear equation, A x = b by direct inversion. Works on small matrices up to size
   * 5.
   *
   * <p>On output b replaced by x. Matrix a may be modified.
   *
   * @param A the matrix A
   * @param b the vector b
   * @return False if the equation is singular (no solution)
   */
  public boolean solveDirectInversion(DenseMatrix64F A, DenseMatrix64F b) {
    createSolver(A.numCols);
    return solveEquation(getInversionSolver(), A, b);
  }

  /**
   * Solves (one) linear equation, A x = b.
   *
   * <p>On output b replaced by x. Matrix a may be modified.
   *
   * @param solver the solver
   * @param A the matrix A
   * @param b the vector b
   * @return False if the equation is singular (no solution)
   */
  private boolean solveEquation(LinearSolver<DenseMatrix64F> solver, DenseMatrix64F A,
      DenseMatrix64F b) {
    final boolean copy = (errorChecking || solver.modifiesB());
    final DenseMatrix64F x = (copy) ? getX() : b;

    if (!solveEquation(solver, A, b, x)) {
      return false;
    }

    // Copy back result if necessary
    if (copy) {
      System.arraycopy(x.data, 0, b.data, 0, b.numRows);
    }

    return true;
  }

  /**
   * Solves (one) linear equation, A x = b.
   *
   * <p>Output is written to x.
   *
   * <p>If checking the solution then A and/or b will not be modified. If no solution checking is
   * enabled then A and/or b could be modified.
   *
   * @param solver the solver
   * @param A the matrix A
   * @param b the vector b
   * @param x the vector x
   * @return False if the equation is singular (no solution)
   */
  private boolean solveEquation(LinearSolver<DenseMatrix64F> solver, DenseMatrix64F A,
      DenseMatrix64F b, DenseMatrix64F x) {
    if (errorChecking) {
      // We need A+b for error checking so solve without modification
      if (!solveEquationSafe(solver, A, b, x)) {
        return false;
      }
      return validate(A, x, b);
    }
    return solveEquationUnsafe(solver, A, b, x);
  }

  /**
   * Solves (one) linear equation, A x = b.
   *
   * <p>Output is written to x.
   *
   * <p>A and/or b will not be modified. If you do not care then use
   * {@link #solveEquationUnsafe(LinearSolver, DenseMatrix64F, DenseMatrix64F, DenseMatrix64F)}
   *
   * @param solver the solver
   * @param A the matrix A
   * @param b the vector b
   * @param x the vector x
   * @return False if the equation is singular (no solution)
   */
  private boolean solveEquationSafe(LinearSolver<DenseMatrix64F> solver, DenseMatrix64F A,
      DenseMatrix64F b, DenseMatrix64F x) {
    if (solver.modifiesA()) {
      A = A.copy();
    }
    if (solver.modifiesB()) {
      b = b.copy();
    }

    if (!initialiseSolver(solver, A)) {
      return false;
    }

    solver.solve(b, x);

    return true;
  }

  /**
   * Solves (one) linear equation, A x = b.
   *
   * <p>Output is written to x.
   *
   * <p>A and/or b may be modified. Check the solver before calling or use
   * {@link #solveEquationSafe(LinearSolver, DenseMatrix64F, DenseMatrix64F, DenseMatrix64F)}
   *
   * @param solver the solver
   * @param A the matrix A
   * @param b the vector b
   * @param x the vector x
   * @return False if the equation is singular (no solution)
   */
  private boolean solveEquationUnsafe(LinearSolver<DenseMatrix64F> solver, DenseMatrix64F A,
      DenseMatrix64F b, DenseMatrix64F x) {
    if (!initialiseSolver(solver, A)) {
      return false;
    }

    solver.solve(b, x);

    return true;
  }

  /**
   * Check that the solution for x satisfies A x = b within the error tolerance.
   *
   * @param A the matrix A
   * @param x the x
   * @param b the b
   * @return true, if successful
   */
  private boolean validate(DenseMatrix64F A, DenseMatrix64F x, DenseMatrix64F b) {
    // Compute A x = b
    for (int i = 0, index = 0; i < b.numRows; i++) {
      double bi = 0;
      for (int j = 0; j < b.numRows; j++) {
        bi += A.data[index++] * x.data[j];
      }
      if (!equal.almostEqualRelativeOrAbsolute(b.data[i], bi)) {
        return false;
      }
    }
    return true;
  }

  /**
   * Solves (one) linear equation, A x = b.
   *
   * <p>On output b replaced by x. Matrix a may be modified.
   *
   * <p>Solve using the CholeskyLDLT method or, if that fails (due to a singular matrix), the
   * PseudoInverse decomposition.
   *
   * @param A the matrix A
   * @param b the vector b
   * @return False if the equation is singular (no solution)
   */
  public boolean solve(DenseMatrix64F A, DenseMatrix64F b) {
    createSolver(A.numCols);

    // Speed tests show the Cholesky solver marginally outperforms the
    // CholeskyLDLT solver as the size increases. Note: The EJML factory
    // returns a Cholesky solver for symmetric positive definite matrices
    // so we use this.

    // Note: Use solveSafe as we need the A and b matrix for the subsequent
    // solve attempt if failure

    if (solveEquationSafe(getCholeskySolver(), A, b, getX())
        && (!errorChecking || validate(A, x, b))) {
      System.arraycopy(x.data, 0, b.data, 0, A.numCols);
      return true;
    }

    // TODO - Count how often the primary method fails on a realistic set of fitting data
    // since the PseudoInverse method is slow. We may want to try a different solver first,
    // e.g. LinearSolver.

    // No need for an explicit solveSafe this time.
    // Use the solve() method which will include validate().
    if (solveEquation(getPseudoInverseSolver(), A, b, x)) {
      System.arraycopy(x.data, 0, b.data, 0, A.numCols);
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
  public boolean solveModifiesA() {
    if (errorChecking) {
      return false;
    }
    return getPseudoInverseSolver().modifiesA();
  }

  /**
   * Computes the inverse of the 'A' matrix passed into the last successful solve method.
   *
   * <p>On output a[n][n] replaced by the inverse of the solved matrix a.
   *
   * @param A the matrix a
   *
   * @return False if the last solve attempt failed, or inversion produces non finite values
   */
  public boolean invertLastA(DenseMatrix64F A) {
    if (lastSuccessfulSolver == null) {
      return false;
    }

    lastSuccessfulSolver.invert(getAinv());

    // Check for NaN or Infinity
    final double[] a_inv = invA.data;
    for (int i = a_inv.length; i-- > 0;) {
      if (!Double.isFinite(a_inv[i])) {
        return false;
      }
    }

    // Q. Should we check the product is the identity matrix?
    // This will require that we have the original matrix A used to initialise the solver.

    System.arraycopy(a_inv, 0, A.data, 0, a_inv.length);

    return true;
  }

  /**
   * Initialise solver.
   *
   * @param solver the solver
   * @param A the matrix A
   * @return true, if successful
   */
  private boolean initialiseSolver(LinearSolver<DenseMatrix64F> solver, DenseMatrix64F A) {
    lastSuccessfulSolver = null;
    if (!solver.setA(A)) {
      return false;
    }
    lastSuccessfulSolver = solver;
    return true;
  }

  /**
   * Creates the solver.
   *
   * @param length the length
   */
  private void createSolver(int length) {
    if (solverSize != length) {
      solverSize = length;

      // Reset size dependent objects
      choleskySolver = null;
      linearSolver = null;
      x = null;
      invA = null;
    }
  }

  /**
   * Gets the working space to store the output solution x.
   *
   * @return the x
   */
  private DenseMatrix64F getX() {
    if (x == null) {
      x = new DenseMatrix64F(solverSize, 1);
    }
    return x;
  }

  /**
   * Gets the working space to store the inverse of A.
   *
   * @return the space for A^-1
   */
  private DenseMatrix64F getAinv() {
    if (invA == null) {
      invA = new DenseMatrix64F(solverSize, solverSize);
    }
    return invA;
  }

  /**
   * Gets the cholesky LDLT solver.
   *
   * @return the cholesky LDLT solver
   */
  private LinearSolver<DenseMatrix64F> getCholeskyLdlTSolver() {
    // This is a Cholesky LDLT solver that should be faster than the Cholesky solver.
    // It only works on symmetric positive definite matrices.
    if (choleskyLdlTSolver == null) {
      choleskyLdlTSolver = new LinearSolverCholLDL();
    }
    return choleskyLdlTSolver;
  }

  /**
   * Gets the cholesky solver.
   *
   * @return the cholesky solver
   */
  private LinearSolver<DenseMatrix64F> getCholeskySolver() {
    if (choleskySolver == null) {
      // This is a Cholesky solver that only works on symmetric positive definite matrices
      choleskySolver = LinearSolverFactory.symmPosDef(solverSize);
    }
    return choleskySolver;
  }

  /**
   * Gets the linear solver.
   *
   * @return the linear solver
   */
  private LinearSolver<DenseMatrix64F> getLinearSolver() {
    if (linearSolver == null) {
      // This should work on any matrix
      linearSolver = LinearSolverFactory.linear(solverSize);
    }
    return linearSolver;
  }

  /**
   * Gets the pseudo inverse solver.
   *
   * @return the pseudo inverse solver
   */
  private LinearSolver<DenseMatrix64F> getPseudoInverseSolver() {
    // The pseudo inverse is constructed using the non-singular sub matrix of A

    if (pseudoInverseSolver == null) {
      pseudoInverseSolver = LinearSolverFactory.pseudoInverse(false);
    }
    return pseudoInverseSolver;
  }

  /**
   * Gets the inversion solver.
   *
   * @return the inversion solver
   */
  private LinearSolver<DenseMatrix64F> getInversionSolver() {
    if (inversionSolver == null) {
      // Supports any matrix up to size 5
      inversionSolver = new InversionSolver();
    }
    return inversionSolver;
  }

  /**
   * Sets the equal.
   *
   * @param equal the equality class to compare that the solution x in A x = b is within tolerance
   */
  public void setEqual(DoubleEquality equal) {
    this.equal = equal;
    errorChecking = (equal != null);
  }

  /**
   * Gets the equal.
   *
   * @return the equality class to compare that the solution x in A x = b is within tolerance
   */
  public DoubleEquality getEqual() {
    return equal;
  }

  /**
   * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
   *
   * @param A the matrix A
   * @return False if the matrix is singular (no solution)
   */
  public boolean invertLinear(DenseMatrix64F A) {
    createSolver(A.numCols);
    return invertUnsafe(getLinearSolver(), A, false);
  }

  /**
   * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
   *
   * @param A the matrix A
   * @return False if the matrix is singular (no solution)
   */
  public boolean invertCholesky(DenseMatrix64F A) {
    createSolver(A.numCols);
    return invertUnsafe(getCholeskySolver(), A, false);
  }

  /**
   * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
   *
   * @param A the matrix A
   * @return False if the matrix is singular (no solution)
   */
  public boolean invertCholeskyLdlT(DenseMatrix64F A) {
    createSolver(A.numCols);
    return invertUnsafe(getCholeskyLdlTSolver(), A, false);
  }

  /**
   * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
   *
   * @param A the matrix A
   * @return False if the matrix is singular (no solution)
   */
  public boolean invertPseudoInverse(DenseMatrix64F A) {
    createSolver(A.numCols);
    return invertUnsafe(getPseudoInverseSolver(), A, true);
  }

  /**
   * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
   *
   * @param A the matrix A
   * @return False if the matrix is singular (no solution)
   */
  public boolean invertDirectInversion(DenseMatrix64F A) {
    createSolver(A.numCols);
    return invertUnsafe(getInversionSolver(), A, false);
  }

  /**
   * Invert symmetric positive definite matrix A and returns only the diagonal.
   *
   * <p>Only supports matrix size up to 5.
   *
   * @param A the matrix A
   * @return The diagonal of A^-1 (or null if the matrix is singular (no solution) or too large)
   */
  public static double[] invertDiagonalDirectInversion(DenseMatrix64F A) {
    return (A.numCols <= UnrolledInverseFromMinorExt.MAX) ? UnrolledInverseFromMinorExt.inv(A)
        : null;
  }

  /**
   * Computes the inverse of the symmetric positive definite matrix. On output a is replaced by the
   * inverse of a.
   *
   * <p>Note: If the matrix is singular then a pseudo inverse will be computed.
   *
   * @param A the matrix a
   * @return False if there is no solution
   */
  public boolean invert(DenseMatrix64F A) {
    createSolver(A.numCols);

    // Speed tests show the Cholesky solver marginally outperforms the
    // CholeskyLDLT solver as the size increases.
    // The DirectInversion solver is faster when the size < 5.
    // Note: The EJML factory returns a Cholesky solver for symmetric
    // positive definite matrices so we use this in preference to a CholeskyLDLT.

    final LinearSolver<DenseMatrix64F> primarySolver =
        (A.numCols < 5) ? getInversionSolver() : getCholeskySolver();
    if (invertSafe(primarySolver, A, false)) {
      return true;
    }

    return invertUnsafe(getPseudoInverseSolver(), A, true);
  }

  /**
   * Checks if an inversion may modify A.
   *
   * @return true, if an inversion may modify A
   * @see #invert(DenseMatrix64F)
   * @see #invert(double[], int)
   */
  public boolean invertModifiesA() {
    if (isInversionTolerance()) {
      return false;
    }
    return getPseudoInverseSolver().modifiesA();
  }

  /**
   * Invert safe.
   *
   * @param solver the solver
   * @param A the matrix A
   * @param pseudoInverse the pseudo inverse
   * @return true, if successful
   */
  private boolean invertSafe(LinearSolver<DenseMatrix64F> solver, DenseMatrix64F A,
      boolean pseudoInverse) {
    final DenseMatrix64F Ain = (solver.modifiesA() || isInversionTolerance()) ? A.copy() : A;
    if (!initialiseSolver(solver, Ain)) {
      return false;
    }
    solver.invert(getAinv());

    // Check for NaN or Infinity
    final double[] a_inv = invA.data;
    for (int i = a_inv.length; i-- > 0;) {
      if (!Double.isFinite(a_inv[i])) {
        return false;
      }
    }

    if (isInversionTolerance() && invalidInversion(A, pseudoInverse)) {
      return false;
    }

    System.arraycopy(a_inv, 0, A.data, 0, a_inv.length);

    return true;
  }

  /**
   * Invert unsafe.
   *
   * @param solver the solver
   * @param A the matrix A
   * @param pseudoInverse the pseudo inverse
   * @return true, if successful
   */
  private boolean invertUnsafe(LinearSolver<DenseMatrix64F> solver, DenseMatrix64F A,
      boolean pseudoInverse) {
    final DenseMatrix64F Ain = (isInversionTolerance()) ? A.copy() : A;
    if (!initialiseSolver(solver, Ain)) {
      return false;
    }
    solver.invert(getAinv());

    // Check for NaN or Infinity
    final double[] inverse = invA.data;
    for (int i = inverse.length; i-- > 0;) {
      if (!Double.isFinite(inverse[i])) {
        Logger.getLogger(getClass().getName()).warning("Inversion not finite");
        return false;
      }
    }

    if (isInversionTolerance() && invalidInversion(A, pseudoInverse)) {
      return false;
    }

    System.arraycopy(inverse, 0, A.data, 0, inverse.length);

    return true;
  }

  /**
   * Invalid inversion.
   *
   * @param A the matrix A
   * @param pseudoInverse the pseudo inverse
   * @return true, if successful
   */
  private boolean invalidInversion(DenseMatrix64F A, boolean pseudoInverse) {
    // Check for the identity matrix:
    // Compute A Ainv = I
    final int n = A.numCols;
    final DenseMatrix64F I = new DenseMatrix64F(n, n);
    CommonOps.mult(A, invA, I);

    if (pseudoInverse) {
      for (int i = n, index = I.data.length; i-- > 0;) {
        for (int j = n; j-- > 0;) {
          if (j == i) {
            --index;
            // If using the pseudo inverse then the diagonal can be zero or 1
            if (invalid(I.data[index], 1) && invalid(I.data[index], 0)) {
              return true;
            }
          } else if (invalid(I.data[--index], 0)) {
            return true;
          }
        }
      }
    } else {
      for (int i = n, index = I.data.length; i-- > 0;) {
        for (int j = n; j-- > 0;) {
          if (invalid(I.data[--index], (j == i) ? 1 : 0)) {
            return true;
          }
        }
      }
    }
    return false;
  }

  /**
   * Invalid.
   *
   * @param e the e
   * @param o the o
   * @return true, if successful
   */
  private boolean invalid(double e, double o) {
    return (Math.abs(e - o) > inversionTolerance);
  }

  /**
   * Computes the inverse of the symmetric positive definite matrix and returns only the diagonal.
   * Will not modify the matrix.
   *
   * <p>Note: If the matrix is singular then a pseudo inverse will be computed.
   *
   * @param A the matrix a
   * @return The diagonal of the inverted matrix (or null)
   */
  public double[] invertDiagonal(DenseMatrix64F A) {
    // Try a fast inversion of the diagonal
    if (A.numCols <= UnrolledInverseFromMinorExt.MAX) {
      final double[] d = UnrolledInverseFromMinorExt.inv(A);
      if (d != null) {
        return d;
      }
    }

    createSolver(A.numCols);

    // Do the LdlT solver only if the fast inversion failed
    if ((A.numCols <= UnrolledInverseFromMinorExt.MAX
        || !invertSafe(getCholeskyLdlTSolver(), A, false))
        // The first inversion failed so try a pseudo-inverse
        && !invertSafe(getPseudoInverseSolver(), A, true)) {
      return null;
    }

    // We reach here when 'a' has been inverted
    final double[] d = new double[A.numCols];
    for (int i = 0, j = 0; i < d.length; i++, j += A.numCols + 1) {
      d[i] = A.get(j);
    }
    return d;
  }

  /**
   * Computes the inverse of the symmetric positive definite matrix and returns only the diagonal.
   * May modify the matrix.
   *
   * <p>Note: If the matrix is singular then a pseudo inverse will be computed.
   *
   * @param A the matrix a
   * @return The diagonal of the inverted matrix (or null)
   */
  private double[] invertDiagonalUnsafe(DenseMatrix64F A) {
    // Try a fast inversion of the diagonal
    if (A.numCols <= UnrolledInverseFromMinorExt.MAX) {
      final double[] d = UnrolledInverseFromMinorExt.inv(A);
      if (d != null) {
        return d;
      }
    }

    createSolver(A.numCols);

    // Do the LdlT solver only if the fast inversion failed
    if ((A.numCols <= UnrolledInverseFromMinorExt.MAX
        || !invertSafe(getCholeskyLdlTSolver(), A, false))
        // The first inversion failed so try a pseudo-inverse
        && !invertUnsafe(getPseudoInverseSolver(), A, true)) {
      return null;
    }

    // We reach here when 'a' has been inverted
    final double[] d = new double[A.numCols];
    for (int i = 0, j = 0; i < d.length; i++, j += A.numCols + 1) {
      d[i] = A.get(j);
    }
    return d;
  }

  /**
   * Convert a dense matrix to a row/column format.
   *
   * @param A the matrix
   * @return the row/column format
   */
  public static double[][] toSquareData(DenseMatrix64F A) {
    final int numRows = A.numRows;
    final int numCols = A.numCols;
    final double[][] out = new double[numRows][];
    for (int i = 0, pos = 0; i < numRows; i++, pos += numRows) {
      out[i] = new double[numCols];
      System.arraycopy(A.data, pos, out[i], 0, numCols);
    }
    return out;
  }

  /**
   * Convert a dense matrix to a row/column format. The output arrays must be the correct size.
   *
   * @param A the matrix
   * @param out the row/column format
   */
  public static void toSquareData(DenseMatrix64F A, double[][] out) {
    final int numRows = A.numRows;
    for (int i = 0, pos = 0; i < numRows; i++, pos += numRows) {
      System.arraycopy(A.data, pos, out[i], 0, numRows);
    }
  }

  /**
   * Create a new DenseMatrix from the input matrix a. Modifications to the matrix are not passed
   * through to the input array! The matrix can be converted back using
   * {@link #toSquareData(DenseMatrix64F, double[][])}.
   *
   * <p>This is provided as a bridge method between the functions that accept primitive arrays and
   * those that accept DenseMatrix.
   *
   * @param a the matrix
   * @return the dense matrix
   */
  public static DenseMatrix64F toA(double[][] a) {
    return new DenseMatrix64F(a);
  }

  /**
   * Create a new DenseMatrix from the input matrix a. Modifications to the matrix are passed
   * through to the input array!
   *
   * <p>This is provided as a bridge method between the functions that accept primitive arrays and
   * those that accept DenseMatrix.
   *
   * @param a the matrix
   * @param n the number of columns/rows
   * @return the dense matrix
   */
  public static DenseMatrix64F toA(double[] a, int n) {
    return DenseMatrix64F.wrap(n, n, a);
  }

  /**
   * Wrap the input array b in a DenseMatrix. Modifications to the matrix are passed through to the
   * input array.
   *
   * <p>This is provided as a bridge method between the functions that accept primitive arrays and
   * those that accept DenseMatrix.
   *
   * @param b the array
   * @return the dense matrix
   */
  public static DenseMatrix64F toB(double[] b) {
    return DenseMatrix64F.wrap(b.length, 1, b);
  }

  // Methods for input of primitive arrays
  // CHECKSTYLE.OFF: OverloadMethodsDeclarationOrder

  /**
   * Solves (one) linear equation, A x = b.
   *
   * <p>On output b replaced by x.
   *
   * @param a the a
   * @param b the b
   * @return False if the equation is singular (no solution)
   */
  public boolean solveLinear(double[][] a, double[] b) {
    return solveLinear(toA(a), toB(b));
  }

  /**
   * Solves (one) linear equation, A x = b.
   *
   * <p>On output b replaced by x.
   *
   * @param a the a
   * @param b the b
   * @return False if the equation is singular (no solution)
   */
  public boolean solveCholesky(double[][] a, double[] b) {
    return solveCholesky(toA(a), toB(b));
  }

  /**
   * Solves (one) linear equation, A x = b.
   *
   * <p>On output b replaced by x.
   *
   * @param a the a
   * @param b the b
   * @return False if the equation is singular (no solution)
   */
  public boolean solveCholeskyLdlT(double[][] a, double[] b) {
    return solveCholeskyLdlT(toA(a), toB(b));
  }

  /**
   * Solves (one) linear equation, A x = b.
   *
   * <p>On output b replaced by x.
   *
   * @param a the a
   * @param b the b
   * @return False if the equation is singular (no solution)
   */
  public boolean solvePseudoInverse(double[][] a, double[] b) {
    return solvePseudoInverse(toA(a), toB(b));
  }

  /**
   * Solves (one) linear equation, A x = b by direct inversion. Works on small matrices up to size
   * 5.
   *
   * <p>On output b replaced by x.
   *
   * @param a the a
   * @param b the b
   * @return False if the equation is singular (no solution)
   */
  public boolean solveDirectInversion(double[][] a, double[] b) {
    return solveDirectInversion(toA(a), toB(b));
  }

  /**
   * Solves (one) linear equation, A x = b.
   *
   * <p>On output b replaced by x.
   *
   * <p>Solve using the CholeskyLDLT method or, if that fails (due to a singular matrix), the
   * PseudoInversion decomposition.
   *
   * @param a the a
   * @param b the b
   * @return False if the equation is singular (no solution)
   */
  public boolean solve(double[][] a, double[] b) {
    return solve(toA(a), toB(b));
  }

  /**
   * Computes the inverse of the 'A' matrix passed into the last successful solve method.
   *
   * <p>On output a[n][n] replaced by the inverse of the solved matrix a. If any column/row index
   * was removed (as it was set to zero in the input matrix) then the resulting column/row index
   * will be set to zero.
   *
   * @param a the matrix a
   * @return False if the last solve attempt failed, or inversion produces non finite values
   */
  public boolean invertLastA(double[][] a) {
    if (lastSuccessfulSolver == null) {
      return false;
    }

    lastSuccessfulSolver.invert(getAinv());

    // Check for NaN or Infinity
    final double[] a_inv = invA.data;
    for (int i = a_inv.length; i-- > 0;) {
      if (!Double.isFinite(a_inv[i])) {
        return false;
      }
    }

    // Q. Should we check the product is the identity matrix?
    // This will require that we have the original matrix A used to initialise the solver.

    toSquareData(invA, a);
    return true;
  }

  /**
   * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
   *
   * @param a the matrix a
   * @return False if the matrix is singular (no solution)
   */
  public boolean invertLinear(double[][] a) {
    final DenseMatrix64F A = toA(a);
    if (!invertLinear(A)) {
      return false;
    }
    toSquareData(A, a);
    return true;
  }

  /**
   * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
   *
   * @param a the matrix a
   * @return False if the matrix is singular (no solution)
   */
  public boolean invertCholesky(double[][] a) {
    final DenseMatrix64F A = toA(a);
    if (!invertCholesky(A)) {
      return false;
    }
    toSquareData(A, a);
    return true;
  }

  /**
   * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
   *
   * @param a the matrix a
   * @return False if the matrix is singular (no solution)
   */
  public boolean invertCholeskyLdlT(double[][] a) {
    final DenseMatrix64F A = toA(a);
    if (!invertCholeskyLdlT(A)) {
      return false;
    }
    toSquareData(A, a);
    return true;
  }

  /**
   * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
   *
   * @param a the matrix a
   * @return False if the matrix is singular (no solution)
   */
  public boolean invertPseudoInverse(double[][] a) {
    final DenseMatrix64F A = toA(a);
    if (!invertPseudoInverse(A)) {
      return false;
    }
    toSquareData(A, a);
    return true;
  }

  /**
   * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
   *
   * @param a the matrix a
   * @return False if the matrix is singular (no solution)
   */
  public boolean invertDirectInversion(double[][] a) {
    final DenseMatrix64F A = toA(a);
    if (!invertDirectInversion(A)) {
      return false;
    }
    toSquareData(A, a);
    return true;
  }

  /**
   * Invert symmetric positive definite matrix A and returns only the diagonal.
   *
   * <p>Only supports matrix size up to 5.
   *
   * @param a the matrix a
   * @return The diagonal of A^-1 (or null if the matrix is singular (no solution) or too large)
   */
  public static double[] invertDiagonalDirectInversion(double[][] a) {
    if (a.length > UnrolledInverseFromMinorExt.MAX) {
      return null;
    }
    return UnrolledInverseFromMinorExt.inv(toA(a));
  }

  /**
   * Computes the inverse of the symmetric positive definite matrix. On output a is replaced by the
   * inverse of a.
   *
   * <p>Note: If the matrix is singular then a pseudo inverse will be computed.
   *
   * @param a the matrix a
   * @return False if there is no solution
   */
  public boolean invert(double[][] a) {
    final DenseMatrix64F A = toA(a);
    if (!invert(A)) {
      return false;
    }
    toSquareData(A, a);
    return true;
  }

  /**
   * Computes the inverse of the symmetric positive definite matrix and returns only the diagonal.
   *
   * <p>Note: If the matrix is singular then a pseudo inverse will be computed.
   *
   * @param a the matrix a
   * @return The diagonal of the inverted matrix (or null)
   */
  public double[] invertDiagonal(double[][] a) {
    // Use the unsafe method as the matrix has been converted so can be modified
    return invertDiagonalUnsafe(new DenseMatrix64F(a));
  }

  /**
   * Solves (one) linear equation, A x = b.
   *
   * <p>On output b replaced by x. Matrix a may be modified.
   *
   * @param a the a
   * @param b the b
   * @return False if the equation is singular (no solution)
   */
  public boolean solveLinear(double[] a, double[] b) {
    return solveLinear(toA(a, b.length), toB(b));
  }

  /**
   * Solves (one) linear equation, A x = b.
   *
   * <p>On output b replaced by x. Matrix a may be modified.
   *
   * @param a the a
   * @param b the b
   * @return False if the equation is singular (no solution)
   */
  public boolean solveCholesky(double[] a, double[] b) {
    return solveCholesky(toA(a, b.length), toB(b));
  }

  /**
   * Solves (one) linear equation, A x = b.
   *
   * <p>On output b replaced by x. Matrix a may be modified.
   *
   * @param a the a
   * @param b the b
   * @return False if the equation is singular (no solution)
   */
  public boolean solveCholeskyLdlT(double[] a, double[] b) {
    return solveCholeskyLdlT(toA(a, b.length), toB(b));
  }

  /**
   * Solves (one) linear equation, A x = b.
   *
   * <p>On output b replaced by x. Matrix a may be modified.
   *
   * @param a the a
   * @param b the b
   * @return False if the equation is singular (no solution)
   */
  public boolean solvePseudoInverse(double[] a, double[] b) {
    return solvePseudoInverse(toA(a, b.length), toB(b));
  }

  /**
   * Solves (one) linear equation, A x = b by direct inversion. Works on small matrices up to size
   * 5.
   *
   * <p>On output b replaced by x. Matrix a may be modified.
   *
   * @param a the a
   * @param b the b
   * @return False if the equation is singular (no solution)
   */
  public boolean solveDirectInversion(double[] a, double[] b) {
    return solveDirectInversion(toA(a, b.length), toB(b));
  }

  /**
   * Solves (one) linear equation, A x = b.
   *
   * <p>On output b replaced by x. Matrix a may be modified.
   *
   * <p>Solve using the CholeskyLDLT method or, if that fails (due to a singular matrix), the
   * PseudoInversion decomposition.
   *
   * @param a the a
   * @param b the b
   * @return False if the equation is singular (no solution)
   */
  public boolean solve(double[] a, double[] b) {
    return solve(toA(a, b.length), toB(b));
  }

  /**
   * Computes the inverse of the 'A' matrix passed into the last successful solve method.
   *
   * <p>On output a[n][n] replaced by the inverse of the solved matrix a. If any column/row index
   * was removed (as it was set to zero in the input matrix) then the resulting column/row index
   * will be set to zero.
   *
   * @param a the matrix a
   * @return False if the last solve attempt failed, or inversion produces non finite values
   */
  public boolean invertLastA(double[] a) {
    if (lastSuccessfulSolver == null) {
      return false;
    }

    lastSuccessfulSolver.invert(getAinv());

    // Check for NaN or Infinity
    final double[] inv = invA.data;
    for (int i = inv.length; i-- > 0;) {
      if (!Double.isFinite(inv[i])) {
        return false;
      }
    }

    // Q. Should we check the product is the identity matrix?
    // This will require that we have the original matrix A used to initialise the solver.

    System.arraycopy(inv, 0, a, 0, inv.length);

    return true;
  }

  /**
   * Computes the inverse of the symmetric positive definite matrix. On output a is replaced by the
   * inverse of a.
   *
   * <p>Note: If the matrix is singular then a pseudo inverse will be computed.
   *
   * @param a the matrix a
   * @param n the number of columns/rows
   * @return False if there is no solution
   */
  public boolean invert(double[] a, int n) {
    return invert(toA(a, n));
  }

  /**
   * Computes the inverse of the symmetric positive definite matrix and returns only the diagonal.
   *
   * <p>Note: If the matrix is singular then a pseudo inverse will be computed.
   *
   * @param a the matrix a
   * @param n the number of columns/rows
   * @return The diagonal of the inverted matrix (or null)
   */
  public double[] invertDiagonal(double[] a, int n) {
    return invertDiagonal(toA(a, n));
  }
}
