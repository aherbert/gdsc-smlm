/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.factory.LinearSolverFactory_DDRM;
import org.ejml.dense.row.linsol.chol.LinearSolverCholLDL_DDRM;
import org.ejml.dense.row.misc.UnrolledInverseFromMinor_DDRM;
import org.ejml.interfaces.decomposition.DecompositionInterface;
import org.ejml.interfaces.linsol.LinearSolver;
import org.ejml.interfaces.linsol.LinearSolverDense;
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

  /** The linear solver. */
  private LinearSolverDense<DMatrixRMaj> linearSolver;

  /** The pseudo inverse solver. */
  private LinearSolverDense<DMatrixRMaj> pseudoInverseSolver;

  /** The cholesky solver. */
  private LinearSolverDense<DMatrixRMaj> choleskySolver;

  /** The cholesky LDLT solver. */
  private LinearSolverDense<DMatrixRMaj> choleskyLdlTSolver;

  /** The inversion solver. */
  private LinearSolverDense<DMatrixRMaj> inversionSolver;

  /** The last successful solver. */
  private LinearSolverDense<DMatrixRMaj> lastSuccessfulSolver;

  /** The vector x. */
  private DMatrixRMaj x;

  /** The inverse of matrix A. */
  private DMatrixRMaj invA;

  /** The solver size. */
  private int solverSize;

  /** The error checking flag. Set to true to check the solution x to linear equations A x = b. */
  private boolean errorChecking;

  /** The class used to check equality (with zero). */
  private DoubleEquality equal;

  /** The inversion tolerance. */
  private double inversionTolerance;

  /**
   * Solve the matrix using direct inversion.
   */
  private static class InversionSolver implements LinearSolverDense<DMatrixRMaj> {
    /** The matrix A. */
    private DMatrixRMaj a;

    /**
     * Instantiates a new inversion solver.
     */
    InversionSolver() {
      // Intentionally empty
    }

    @Override
    public boolean setA(DMatrixRMaj a) {
      if (a.numCols <= UnrolledInverseFromMinor_DDRM.MAX) {
        // Direct inversion using the determinant
        if (a.numCols >= 2) {
          UnrolledInverseFromMinor_DDRM.inv(a, a);
        } else {
          a.set(0, 1.0 / a.get(0));
        }

        // Check for NaN or Infinity
        for (int i = a.data.length; i-- > 0;) {
          if (!Double.isFinite(a.data[i])) {
            return false;
          }
        }

        this.a = a;
        return true;
      }
      return false;
    }

    @Override
    public double quality() {
      return 0;
    }

    @Override
    public void solve(DMatrixRMaj b, DMatrixRMaj x) {
      CommonOps_DDRM.mult(a, b, x);
    }

    @Override
    public void invert(DMatrixRMaj aInv) {
      System.arraycopy(a.data, 0, aInv.data, 0, a.data.length);
    }

    @Override
    public boolean modifiesA() {
      return true;
    }

    @Override
    public boolean modifiesB() {
      return false;
    }

    @Override
    public <Decomposition extends DecompositionInterface> Decomposition getDecomposition() {
      return null;
    }
  }

  /**
   * Instantiates a new EJML linear solver.
   */
  public EjmlLinearSolver() {
    // Intentionally empty
  }

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
   * @param a the matrix A
   * @param b the vector b
   * @return False if the equation is singular (no solution)
   */
  public boolean solveLinear(DMatrixRMaj a, DMatrixRMaj b) {
    createSolver(a.numCols);
    return solveEquation(getLinearSolver(), a, b);
  }

  /**
   * Solves (one) linear equation, A x = b.
   *
   * <p>On output b replaced by x. Matrix a may be modified.
   *
   * @param a the matrix A
   * @param b the vector b
   * @return False if the equation is singular (no solution)
   */
  public boolean solveCholesky(DMatrixRMaj a, DMatrixRMaj b) {
    createSolver(a.numCols);
    return solveEquation(getCholeskySolver(), a, b);
  }

  /**
   * Solves (one) linear equation, A x = b.
   *
   * <p>On output b replaced by x. Matrix a may be modified.
   *
   * @param a the matrix A
   * @param b the vector b
   * @return False if the equation is singular (no solution)
   */
  public boolean solveCholeskyLdlT(DMatrixRMaj a, DMatrixRMaj b) {
    createSolver(a.numCols);
    return solveEquation(getCholeskyLdlTSolver(), a, b);
  }

  /**
   * Solves (one) linear equation, A x = b.
   *
   * <p>On output b replaced by x. Matrix a may be modified.
   *
   * @param a the matrix A
   * @param b the vector b
   * @return False if the equation is singular (no solution)
   */
  public boolean solvePseudoInverse(DMatrixRMaj a, DMatrixRMaj b) {
    createSolver(a.numCols);
    return solveEquation(getPseudoInverseSolver(), a, b);
  }

  /**
   * Solves (one) linear equation, A x = b by direct inversion. Works on small matrices up to size
   * 5.
   *
   * <p>On output b replaced by x. Matrix a may be modified.
   *
   * @param a the matrix A
   * @param b the vector b
   * @return False if the equation is singular (no solution)
   */
  public boolean solveDirectInversion(DMatrixRMaj a, DMatrixRMaj b) {
    createSolver(a.numCols);
    return solveEquation(getInversionSolver(), a, b);
  }

  /**
   * Solves (one) linear equation, A x = b.
   *
   * <p>On output b replaced by x. Matrix a may be modified.
   *
   * @param solver the solver
   * @param a the matrix A
   * @param b the vector b
   * @return False if the equation is singular (no solution)
   */
  private boolean solveEquation(LinearSolverDense<DMatrixRMaj> solver, DMatrixRMaj a,
      DMatrixRMaj b) {
    final boolean copy = (errorChecking || solver.modifiesB());
    final DMatrixRMaj x = (copy) ? getX() : b;

    if (!solveEquation(solver, a, b, x)) {
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
   * @param a the matrix A
   * @param b the vector b
   * @param x the vector x
   * @return False if the equation is singular (no solution)
   */
  private boolean solveEquation(LinearSolverDense<DMatrixRMaj> solver, DMatrixRMaj a,
      DMatrixRMaj b, DMatrixRMaj x) {
    if (errorChecking) {
      // We need A+b for error checking so solve without modification
      if (!solveEquationSafe(solver, a, b, x)) {
        return false;
      }
      return validate(a, x, b);
    }
    return solveEquationUnsafe(solver, a, b, x);
  }

  /**
   * Solves (one) linear equation, A x = b.
   *
   * <p>Output is written to x.
   *
   * <p>A and/or b will not be modified. If you do not care then use
   * {@link #solveEquationUnsafe(LinearSolver, DMatrixRMaj, DMatrixRMaj, DMatrixRMaj)}
   *
   * @param solver the solver
   * @param a the matrix A
   * @param b the vector b
   * @param x the vector x
   * @return False if the equation is singular (no solution)
   */
  private boolean solveEquationSafe(LinearSolverDense<DMatrixRMaj> solver, DMatrixRMaj a,
      DMatrixRMaj b, DMatrixRMaj x) {
    if (solver.modifiesA()) {
      a = a.copy();
    }
    if (solver.modifiesB()) {
      b = b.copy();
    }

    if (!initialiseSolver(solver, a)) {
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
   * {@link #solveEquationSafe(LinearSolver, DMatrixRMaj, DMatrixRMaj, DMatrixRMaj)}
   *
   * @param solver the solver
   * @param a the matrix A
   * @param b the vector b
   * @param x the vector x
   * @return False if the equation is singular (no solution)
   */
  private boolean solveEquationUnsafe(LinearSolverDense<DMatrixRMaj> solver, DMatrixRMaj a,
      DMatrixRMaj b, DMatrixRMaj x) {
    if (!initialiseSolver(solver, a)) {
      return false;
    }

    solver.solve(b, x);

    return true;
  }

  /**
   * Check that the solution for x satisfies A x = b within the error tolerance.
   *
   * @param a the matrix A
   * @param x the x
   * @param b the b
   * @return true, if successful
   */
  private boolean validate(DMatrixRMaj a, DMatrixRMaj x, DMatrixRMaj b) {
    // Compute A x = b
    for (int i = 0, index = 0; i < b.numRows; i++) {
      double bi = 0;
      for (int j = 0; j < b.numRows; j++) {
        bi += a.data[index++] * x.data[j];
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
   * @param a the matrix A
   * @param b the vector b
   * @return False if the equation is singular (no solution)
   */
  public boolean solve(DMatrixRMaj a, DMatrixRMaj b) {
    createSolver(a.numCols);

    // Speed tests show the Cholesky solver marginally outperforms the
    // CholeskyLDLT solver as the size increases. Note: The EJML factory
    // returns a Cholesky solver for symmetric positive definite matrices
    // so we use this.

    // Note: Use solveSafe as we need the A and b matrix for the subsequent
    // solve attempt if failure

    if (solveEquationSafe(getCholeskySolver(), a, b, getX())
        && (!errorChecking || validate(a, x, b))) {
      System.arraycopy(x.data, 0, b.data, 0, a.numCols);
      return true;
    }

    // TODO - Count how often the primary method fails on a realistic set of fitting data
    // since the PseudoInverse method is slow. We may want to try a different solver first,
    // e.g. LinearSolver.

    // No need for an explicit solveSafe this time.
    // Use the solve() method which will include validate().
    if (solveEquation(getPseudoInverseSolver(), a, b, x)) {
      System.arraycopy(x.data, 0, b.data, 0, a.numCols);
      return true;
    }

    return false;
  }

  /**
   * Checks if a solve may modify A.
   *
   * @return true, if a solve may modify A
   * @see #solve(DMatrixRMaj, DMatrixRMaj)
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
   * @param a the matrix A
   *
   * @return False if the last solve attempt failed, or inversion produces non finite values
   */
  public boolean invertLastA(DMatrixRMaj a) {
    if (lastSuccessfulSolver == null) {
      return false;
    }

    lastSuccessfulSolver.invert(getAinv());

    // Check for NaN or Infinity
    final double[] aInv = invA.data;
    for (int i = aInv.length; i-- > 0;) {
      if (!Double.isFinite(aInv[i])) {
        return false;
      }
    }

    // Q. Should we check the product is the identity matrix?
    // This will require that we have the original matrix A used to initialise the solver.

    System.arraycopy(aInv, 0, a.data, 0, aInv.length);

    return true;
  }

  /**
   * Initialise solver.
   *
   * @param solver the solver
   * @param a the matrix A
   * @return true, if successful
   */
  private boolean initialiseSolver(LinearSolverDense<DMatrixRMaj> solver, DMatrixRMaj a) {
    lastSuccessfulSolver = null;
    if (!solver.setA(a)) {
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
  private DMatrixRMaj getX() {
    if (x == null) {
      x = new DMatrixRMaj(solverSize, 1);
    }
    return x;
  }

  /**
   * Gets the working space to store the inverse of A.
   *
   * @return the space for A^-1
   */
  private DMatrixRMaj getAinv() {
    if (invA == null) {
      invA = new DMatrixRMaj(solverSize, solverSize);
    }
    return invA;
  }

  /**
   * Gets the cholesky LDLT solver.
   *
   * @return the cholesky LDLT solver
   */
  private LinearSolverDense<DMatrixRMaj> getCholeskyLdlTSolver() {
    // This is a Cholesky LDLT solver that should be faster than the Cholesky solver.
    // It only works on symmetric positive definite matrices.
    if (choleskyLdlTSolver == null) {
      choleskyLdlTSolver = new LinearSolverCholLDL_DDRM();
    }
    return choleskyLdlTSolver;
  }

  /**
   * Gets the cholesky solver.
   *
   * @return the cholesky solver
   */
  private LinearSolverDense<DMatrixRMaj> getCholeskySolver() {
    if (choleskySolver == null) {
      // This is a Cholesky solver that only works on symmetric positive definite matrices
      choleskySolver = LinearSolverFactory_DDRM.symmPosDef(solverSize);
    }
    return choleskySolver;
  }

  /**
   * Gets the linear solver.
   *
   * @return the linear solver
   */
  private LinearSolverDense<DMatrixRMaj> getLinearSolver() {
    if (linearSolver == null) {
      // This should work on any matrix
      linearSolver = LinearSolverFactory_DDRM.linear(solverSize);
    }
    return linearSolver;
  }

  /**
   * Gets the pseudo inverse solver.
   *
   * @return the pseudo inverse solver
   */
  private LinearSolverDense<DMatrixRMaj> getPseudoInverseSolver() {
    // The pseudo inverse is constructed using the non-singular sub matrix of A

    if (pseudoInverseSolver == null) {
      pseudoInverseSolver = LinearSolverFactory_DDRM.pseudoInverse(false);
    }
    return pseudoInverseSolver;
  }

  /**
   * Gets the inversion solver.
   *
   * @return the inversion solver
   */
  private LinearSolverDense<DMatrixRMaj> getInversionSolver() {
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
   * @param a the matrix A
   * @return False if the matrix is singular (no solution)
   */
  public boolean invertLinear(DMatrixRMaj a) {
    createSolver(a.numCols);
    return invertUnsafe(getLinearSolver(), a, false);
  }

  /**
   * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
   *
   * @param a the matrix A
   * @return False if the matrix is singular (no solution)
   */
  public boolean invertCholesky(DMatrixRMaj a) {
    createSolver(a.numCols);
    return invertUnsafe(getCholeskySolver(), a, false);
  }

  /**
   * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
   *
   * @param a the matrix A
   * @return False if the matrix is singular (no solution)
   */
  public boolean invertCholeskyLdlT(DMatrixRMaj a) {
    createSolver(a.numCols);
    return invertUnsafe(getCholeskyLdlTSolver(), a, false);
  }

  /**
   * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
   *
   * @param a the matrix A
   * @return False if the matrix is singular (no solution)
   */
  public boolean invertPseudoInverse(DMatrixRMaj a) {
    createSolver(a.numCols);
    return invertUnsafe(getPseudoInverseSolver(), a, true);
  }

  /**
   * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
   *
   * @param a the matrix A
   * @return False if the matrix is singular (no solution)
   */
  public boolean invertDirectInversion(DMatrixRMaj a) {
    createSolver(a.numCols);
    return invertUnsafe(getInversionSolver(), a, false);
  }

  /**
   * Invert symmetric positive definite matrix A and returns only the diagonal.
   *
   * <p>Only supports matrix size up to 5.
   *
   * @param a the matrix A
   * @return The diagonal of A^-1 (or null if the matrix is singular (no solution) or too large)
   */
  public static double[] invertDiagonalDirectInversion(DMatrixRMaj a) {
    return (a.numCols <= UnrolledInverseFromMinorExt.MAX) ? UnrolledInverseFromMinorExt.inv(a)
        : null;
  }

  /**
   * Computes the inverse of the symmetric positive definite matrix. On output a is replaced by the
   * inverse of a.
   *
   * <p>Note: If the matrix is singular then a pseudo inverse will be computed.
   *
   * @param a the matrix A
   * @return False if there is no solution
   */
  public boolean invert(DMatrixRMaj a) {
    createSolver(a.numCols);

    // Speed tests show the Cholesky solver marginally outperforms the
    // CholeskyLDLT solver as the size increases.
    // The DirectInversion solver is faster when the size < 5.
    // Note: The EJML factory returns a Cholesky solver for symmetric
    // positive definite matrices so we use this in preference to a CholeskyLDLT.

    final LinearSolverDense<DMatrixRMaj> primarySolver =
        (a.numCols < 5) ? getInversionSolver() : getCholeskySolver();
    if (invertSafe(primarySolver, a, false)) {
      return true;
    }

    return invertUnsafe(getPseudoInverseSolver(), a, true);
  }

  /**
   * Checks if an inversion may modify A.
   *
   * @return true, if an inversion may modify A
   * @see #invert(DMatrixRMaj)
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
   * @param a the matrix A
   * @param pseudoInverse the pseudo inverse
   * @return true, if successful
   */
  private boolean invertSafe(LinearSolverDense<DMatrixRMaj> solver, DMatrixRMaj a,
      boolean pseudoInverse) {
    final DMatrixRMaj aa = (solver.modifiesA() || isInversionTolerance()) ? a.copy() : a;
    if (!initialiseSolver(solver, aa)) {
      return false;
    }
    solver.invert(getAinv());

    // Check for NaN or Infinity
    final double[] aInv = invA.data;
    for (int i = aInv.length; i-- > 0;) {
      if (!Double.isFinite(aInv[i])) {
        return false;
      }
    }

    if (isInversionTolerance() && invalidInversion(a, pseudoInverse)) {
      return false;
    }

    System.arraycopy(aInv, 0, a.data, 0, aInv.length);

    return true;
  }

  /**
   * Invert unsafe.
   *
   * @param solver the solver
   * @param a the matrix A
   * @param pseudoInverse the pseudo inverse
   * @return true, if successful
   */
  private boolean invertUnsafe(LinearSolverDense<DMatrixRMaj> solver, DMatrixRMaj a,
      boolean pseudoInverse) {
    final DMatrixRMaj aa = (isInversionTolerance()) ? a.copy() : a;
    if (!initialiseSolver(solver, aa)) {
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

    if (isInversionTolerance() && invalidInversion(a, pseudoInverse)) {
      return false;
    }

    System.arraycopy(inverse, 0, a.data, 0, inverse.length);

    return true;
  }

  /**
   * Invalid inversion.
   *
   * @param a the matrix A
   * @param pseudoInverse the pseudo inverse
   * @return true, if successful
   */
  private boolean invalidInversion(DMatrixRMaj a, boolean pseudoInverse) {
    // Check for the identity matrix:
    // Compute A Ainv = I
    final int n = a.numCols;
    final DMatrixRMaj eye = new DMatrixRMaj(n, n);
    CommonOps_DDRM.mult(a, invA, eye);

    if (pseudoInverse) {
      for (int i = n, index = eye.data.length; i-- > 0;) {
        for (int j = n; j-- > 0;) {
          if (j == i) {
            --index;
            // If using the pseudo inverse then the diagonal can be zero or 1
            if (invalid(eye.data[index], 1) && invalid(eye.data[index], 0)) {
              return true;
            }
          } else if (invalid(eye.data[--index], 0)) {
            return true;
          }
        }
      }
    } else {
      for (int i = n, index = eye.data.length; i-- > 0;) {
        for (int j = n; j-- > 0;) {
          if (invalid(eye.data[--index], (j == i) ? 1 : 0)) {
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
   * @param a the matrix A
   * @return The diagonal of the inverted matrix (or null)
   */
  public double[] invertDiagonal(DMatrixRMaj a) {
    // Try a fast inversion of the diagonal
    if (a.numCols <= UnrolledInverseFromMinorExt.MAX) {
      final double[] d = UnrolledInverseFromMinorExt.inv(a);
      if (d != null) {
        return d;
      }
    }

    createSolver(a.numCols);

    // Do the LdlT solver only if the fast inversion failed
    if ((a.numCols <= UnrolledInverseFromMinorExt.MAX
        || !invertSafe(getCholeskyLdlTSolver(), a, false))
        // The first inversion failed so try a pseudo-inverse
        && !invertSafe(getPseudoInverseSolver(), a, true)) {
      return null;
    }

    // We reach here when 'a' has been inverted
    final double[] d = new double[a.numCols];
    for (int i = 0, j = 0; i < d.length; i++, j += a.numCols + 1) {
      d[i] = a.get(j);
    }
    return d;
  }

  /**
   * Computes the inverse of the symmetric positive definite matrix and returns only the diagonal.
   * May modify the matrix.
   *
   * <p>Note: If the matrix is singular then a pseudo inverse will be computed.
   *
   * @param a the matrix A
   * @return The diagonal of the inverted matrix (or null)
   */
  private double[] invertDiagonalUnsafe(DMatrixRMaj a) {
    // Try a fast inversion of the diagonal
    if (a.numCols <= UnrolledInverseFromMinorExt.MAX) {
      final double[] d = UnrolledInverseFromMinorExt.inv(a);
      if (d != null) {
        return d;
      }
    }

    createSolver(a.numCols);

    // Do the LdlT solver only if the fast inversion failed
    if ((a.numCols <= UnrolledInverseFromMinorExt.MAX
        || !invertSafe(getCholeskyLdlTSolver(), a, false))
        // The first inversion failed so try a pseudo-inverse
        && !invertUnsafe(getPseudoInverseSolver(), a, true)) {
      return null;
    }

    // We reach here when 'a' has been inverted
    final double[] d = new double[a.numCols];
    for (int i = 0, j = 0; i < d.length; i++, j += a.numCols + 1) {
      d[i] = a.get(j);
    }
    return d;
  }

  /**
   * Convert a dense matrix to a row/column format.
   *
   * @param a the matrix
   * @return the row/column format
   */
  public static double[][] toSquareData(DMatrixRMaj a) {
    final int numRows = a.numRows;
    final int numCols = a.numCols;
    final double[][] out = new double[numRows][];
    for (int i = 0, pos = 0; i < numRows; i++, pos += numRows) {
      out[i] = new double[numCols];
      System.arraycopy(a.data, pos, out[i], 0, numCols);
    }
    return out;
  }

  /**
   * Convert a dense matrix to a row/column format. The output arrays must be the correct size.
   *
   * @param a the matrix
   * @param out the row/column format
   */
  public static void toSquareData(DMatrixRMaj a, double[][] out) {
    final int numRows = a.numRows;
    for (int i = 0, pos = 0; i < numRows; i++, pos += numRows) {
      System.arraycopy(a.data, pos, out[i], 0, numRows);
    }
  }

  /**
   * Create a new DenseMatrix from the input matrix a. Modifications to the matrix Are not passed
   * through to the input array! The matrix can be converted back using
   * {@link #toSquareData(DMatrixRMaj, double[][])}.
   *
   * <p>This is provided as a bridge method between the functions that accept primitive arrays and
   * those that accept DenseMatrix.
   *
   * @param a the matrix
   * @return the dense matrix
   */
  public static DMatrixRMaj toA(double[][] a) {
    return new DMatrixRMaj(a);
  }

  /**
   * Create a new DenseMatrix from the input matrix a. Modifications to the matrix Are passed
   * through to the input array!
   *
   * <p>This is provided as a bridge method between the functions that accept primitive arrays and
   * those that accept DenseMatrix.
   *
   * @param a the matrix
   * @param n the number of columns/rows
   * @return the dense matrix
   */
  public static DMatrixRMaj toA(double[] a, int n) {
    return DMatrixRMaj.wrap(n, n, a);
  }

  /**
   * Wrap the input array b in a DenseMatrix. Modifications to the matrix Are passed through to the
   * input array.
   *
   * <p>This is provided as a bridge method between the functions that accept primitive arrays and
   * those that accept DenseMatrix.
   *
   * @param b the array
   * @return the dense matrix
   */
  public static DMatrixRMaj toB(double[] b) {
    return DMatrixRMaj.wrap(b.length, 1, b);
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
   * @param a the matrix A
   * @return False if the last solve attempt failed, or inversion produces non finite values
   */
  public boolean invertLastA(double[][] a) {
    if (lastSuccessfulSolver == null) {
      return false;
    }

    lastSuccessfulSolver.invert(getAinv());

    // Check for NaN or Infinity
    final double[] aInv = invA.data;
    for (int i = aInv.length; i-- > 0;) {
      if (!Double.isFinite(aInv[i])) {
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
   * @param a the matrix A
   * @return False if the matrix is singular (no solution)
   */
  public boolean invertLinear(double[][] a) {
    final DMatrixRMaj aa = toA(a);
    if (!invertLinear(aa)) {
      return false;
    }
    toSquareData(aa, a);
    return true;
  }

  /**
   * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
   *
   * @param a the matrix A
   * @return False if the matrix is singular (no solution)
   */
  public boolean invertCholesky(double[][] a) {
    final DMatrixRMaj aa = toA(a);
    if (!invertCholesky(aa)) {
      return false;
    }
    toSquareData(aa, a);
    return true;
  }

  /**
   * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
   *
   * @param a the matrix A
   * @return False if the matrix is singular (no solution)
   */
  public boolean invertCholeskyLdlT(double[][] a) {
    final DMatrixRMaj aa = toA(a);
    if (!invertCholeskyLdlT(aa)) {
      return false;
    }
    toSquareData(aa, a);
    return true;
  }

  /**
   * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
   *
   * @param a the matrix A
   * @return False if the matrix is singular (no solution)
   */
  public boolean invertPseudoInverse(double[][] a) {
    final DMatrixRMaj aa = toA(a);
    if (!invertPseudoInverse(aa)) {
      return false;
    }
    toSquareData(aa, a);
    return true;
  }

  /**
   * Invert symmetric positive definite matrix A. On output a replaced by A^-1.
   *
   * @param a the matrix A
   * @return False if the matrix is singular (no solution)
   */
  public boolean invertDirectInversion(double[][] a) {
    final DMatrixRMaj aa = toA(a);
    if (!invertDirectInversion(aa)) {
      return false;
    }
    toSquareData(aa, a);
    return true;
  }

  /**
   * Invert symmetric positive definite matrix A and returns only the diagonal.
   *
   * <p>Only supports matrix size up to 5.
   *
   * @param a the matrix A
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
   * @param a the matrix A
   * @return False if there is no solution
   */
  public boolean invert(double[][] a) {
    final DMatrixRMaj aa = toA(a);
    if (!invert(aa)) {
      return false;
    }
    toSquareData(aa, a);
    return true;
  }

  /**
   * Computes the inverse of the symmetric positive definite matrix and returns only the diagonal.
   *
   * <p>Note: If the matrix is singular then a pseudo inverse will be computed.
   *
   * @param a the matrix A
   * @return The diagonal of the inverted matrix (or null)
   */
  public double[] invertDiagonal(double[][] a) {
    // Use the unsafe method as the matrix has been converted so can be modified
    return invertDiagonalUnsafe(new DMatrixRMaj(a));
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
   * @param a the matrix A
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
   * @param a the matrix A
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
   * @param a the matrix A
   * @param n the number of columns/rows
   * @return The diagonal of the inverted matrix (or null)
   */
  public double[] invertDiagonal(double[] a, int n) {
    return invertDiagonal(toA(a, n));
  }
}
