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

package uk.ac.sussex.gdsc.smlm.fitting;

import uk.ac.sussex.gdsc.smlm.fitting.linear.EJMLLinearSolver;

import org.ejml.data.DenseMatrix64F;

import java.util.Arrays;

/**
 * Container for the Fisher information, a symmetric positive definite matrix containing the amount
 * of information that an observable random variable X carries about an unknown parameter θ of a
 * distribution that models X.
 */
public class FisherInformationMatrix {
  /** The default inversion tolerance. */
  public static final double DEFAULT_INVERSION_TOLERANCE = 1e-2;

  private static final byte YES = 1;
  private static final byte UNKNOWN = 0;
  private static final byte NO = -1;

  private final DenseMatrix64F m;
  private double[] crlb = null;
  private byte inverted = UNKNOWN;
  private double inversionTolerance = 0;

  /**
   * Instantiates a new fisher information matrix.
   *
   * @param m the fisher information matrix
   * @param inversionTolerance the inversion tolerance
   */
  public FisherInformationMatrix(double[][] m, double inversionTolerance) {
    this(EJMLLinearSolver.toA(m), inversionTolerance);
  }

  /**
   * Instantiates a new fisher information matrix.
   *
   * @param m the fisher information matrix
   * @param n the number of columns/rows
   * @param inversionTolerance the inversion tolerance
   */
  public FisherInformationMatrix(double[] m, int n, double inversionTolerance) {
    this(EJMLLinearSolver.toA(m, n), inversionTolerance);
  }

  /**
   * Instantiates a new fisher information matrix.
   *
   * @param m the fisher information matrix
   * @param inversionTolerance the inversion tolerance
   */
  public FisherInformationMatrix(DenseMatrix64F m, double inversionTolerance) {
    this.m = m;
    setInversionTolerance(inversionTolerance);
  }

  /**
   * Instantiates a new fisher information matrix.
   *
   * @param m the fisher information matrix
   */
  public FisherInformationMatrix(double[][] m) {
    this(EJMLLinearSolver.toA(m));
  }

  /**
   * Instantiates a new fisher information matrix.
   *
   * @param m the fisher information matrix
   * @param n the number of columns/rows
   */
  public FisherInformationMatrix(double[] m, int n) {
    this(EJMLLinearSolver.toA(m, n));
  }

  /**
   * Instantiates a new fisher information matrix.
   *
   * @param m the fisher information matrix
   */
  public FisherInformationMatrix(DenseMatrix64F m) {
    this(m, DEFAULT_INVERSION_TOLERANCE);
  }

  /**
   * Create a subset of this matrix using the specified row/column indices.
   *
   * @param indices the indices
   * @return the subset fisher information matrix
   */
  public FisherInformationMatrix subset(int[] indices) {
    final int n = m.getNumCols();

    // Check the indices are within the matrix
    for (int i = 0; i < indices.length; i++) {
      if (indices[i] < 0 || indices[i] >= n) {
        throw new IllegalArgumentException("Indices must be >=0 and <" + n + ": " + indices[i]);
      }
    }

    final double[] in = m.getData();
    final int m = indices.length;
    final double[] out = new double[m * m];
    for (int i = 0, ii = 0; i < indices.length; i++) {
      final int index = indices[i] * n;
      for (int j = 0; j < indices.length; j++, ii++) {
        out[ii] = in[index + indices[j]];
      }
    }

    return new FisherInformationMatrix(out, m);
  }

  private void invert() {
    if (inverted != UNKNOWN) {
      return;
    }

    if (m.numCols == 0) {
      // Nothing to do
      crlb = new double[0];
      inverted = YES;
      return;
    }

    inverted = NO;

    // Matrix inversion
    final EJMLLinearSolver solver = EJMLLinearSolver.createForInversion(inversionTolerance);
    final double[] crlb = solver.invertDiagonal(m); // Does not modify the matrix
    if (crlb != null) {
      // Check all diagonal values are zero or above
      if (inversionTolerance > 0) {
        // Already checked so just ignore values just below zero
        for (int i = m.numCols; i-- > 0;) {
          if (crlb[i] < 0) {
            crlb[i] = 0;
          }
        }
      } else {
        // A small error is OK
        for (int i = m.numCols; i-- > 0;) {
          if (crlb[i] < 0) {
            if (crlb[i] > -DEFAULT_INVERSION_TOLERANCE) {
              crlb[i] = 0;
              continue;
            }
            return;
          }
        }
      }

      // Check all diagonal values are zero or above

      inverted = YES;
      this.crlb = crlb;
    }
  }

  /**
   * Compute the Cramér–Rao Lower Bound (CRLB) variance for fitted variables using the central
   * diagonal of the inverted Fisher information matrix.
   *
   * <p>The information matrix is inverted and the central diagonal returned.
   *
   * @return CRLB (or null if inversion failed)
   */
  public double[] crlb() {
    return crlb(false);
  }

  /**
   * Compute the Cramér–Rao Lower Bound (CRLB) variance for fitted variables using the central
   * diagonal of the inverted Fisher information matrix.
   *
   * <p>The information matrix is inverted and the central diagonal returned. If the inversion fails
   * then the routine optionally returns the reciprocal of the diagonal element to find a (possibly
   * loose) lower bound.
   *
   * @param allowReciprocal the allow reciprocal flag
   * @return CRLB (or null if inversion failed and the reciprocal is not used)
   */
  public double[] crlb(boolean allowReciprocal) {
    invert();

    if (inverted == YES) {
      return crlb;
    }

    if (allowReciprocal) {
      return crlbReciprocal();
    }

    return null;
  }

  /**
   * Compute the Cramér–Rao Lower Bound (CRLB) variance for fitted variables using the reciprocal of
   * the central diagonal of the Fisher information matrix.
   *
   * <p>The information matrix is NOT inverted. The reciprocal of the central diagonal returned for
   * a (possibly loose) lower bound.
   *
   * @return CRLB (or null if inversion failed)
   */
  public double[] crlbReciprocal() {
    final double[] crlb = new double[m.numCols];
    for (int i = 0, j = 0, n = m.numCols; i < n; i++, j += n + 1) {
      crlb[i] = reciprocal(m.data[j]);
    }
    return crlb;
  }

  /**
   * Compute the Cramér–Rao Lower Bound (CRLB) for fitted variables using the central diagonal of
   * the inverted Fisher information matrix.
   *
   * <p>The information matrix is inverted and the central diagonal returned.
   *
   * @return CRLB (or null if inversion failed)
   */
  public double[] crlbSqrt() {
    return crlbSqrt(false);
  }

  /**
   * Compute the Cramér–Rao Lower Bound (CRLB) for fitted variables using the central diagonal of
   * the inverted Fisher information matrix.
   *
   * <p>The information matrix is inverted and the square root of the central diagonal returned. If
   * the inversion fails then the routine optionally returns the square root of the reciprocal of
   * the diagonal element to find a (possibly loose) lower bound.
   *
   * @param allowReciprocal the allow reciprocal flag
   * @return CRLB (or null if inversion failed and the reciprocal is not used)
   */
  public double[] crlbSqrt(boolean allowReciprocal) {
    invert();

    if (inverted == YES) {
      final double[] crlb = new double[this.crlb.length];
      for (int i = crlb.length; i-- > 0;) {
        crlb[i] = Math.sqrt(this.crlb[i]);
      }
      return crlb;
    }

    if (allowReciprocal) {
      return crlbReciprocalSqrt();
    }

    return null;
  }

  /**
   * Compute the Cramér–Rao Lower Bound (CRLB) for fitted variables using the reciprocal of the
   * central diagonal of the Fisher information matrix.
   *
   * <p>The information matrix is NOT inverted. Uses the square root of the reciprocal of the
   * central diagonal returned for a (possibly loose) lower bound.
   *
   * @return CRLB (or null if inversion failed)
   */
  public double[] crlbReciprocalSqrt() {
    final double[] crlb = new double[m.numCols];
    for (int i = 0, j = 0, n = m.numCols; i < n; i++, j += n + 1) {
      crlb[i] = reciprocalSqrt(m.data[j]);
    }
    return crlb;
  }

  /**
   * Compute the Cramér–Rao Lower Bound (CRLB) variance for fitted variables using the reciprocal of
   * the central diagonal of the Fisher information matrix.
   *
   * <p>The information matrix is NOT inverted. Uses the reciprocal of the central diagonal returned
   * for a (possibly loose) lower bound.
   *
   * @param m the fisher information matrix
   * @return CRLB
   */
  public static double[] crlbReciprocal(double[][] m) {
    final int n = m.length;
    final double[] crlb = new double[n];
    for (int i = 0; i < n; i++) {
      crlb[i] = reciprocal(m[i][i]);
    }
    return crlb;
  }

  /**
   * Compute the Cramér–Rao Lower Bound (CRLB) for fitted variables using the reciprocal of the
   * central diagonal of the Fisher information matrix.
   *
   * <p>The information matrix is NOT inverted. Uses the square root of the reciprocal of the
   * central diagonal returned for a (possibly loose) lower bound.
   *
   * @param m the fisher information matrix
   * @return CRLB
   */
  public static double[] crlbReciprocalSqrt(double[][] m) {
    final int n = m.length;
    final double[] crlb = new double[n];
    for (int i = 0; i < n; i++) {
      crlb[i] = reciprocalSqrt(m[i][i]);
    }
    return crlb;
  }

  /**
   * Return the reciprocal of the input. If the number is not strictly positive then zero is
   * returned. Note that zero can only be returned if there was no Fisher information. This is done
   * to match the return value from matrix inversion when there is no Fisher information for a
   * parameter i within the matrix. In that case the zero column and row is removed from the matrix
   * before inversion and the inverted matrix contains zeros.
   *
   * <p>The reciprocal of the diagonal element of the Fisher information matrix is a (possibly
   * loose) lower bound.
   *
   * @param d the input value
   * @return the reciprocal of the square root of the input value
   */
  public static double reciprocal(double d) {
    return (d > 0) ? 1.0 / d : 0;
  }

  /**
   * Return the reciprocal of the square root of the input. If the number is not strictly positive
   * then zero is returned. Note that zero can only be returned if there was no Fisher information.
   * This is done to match the return value from matrix inversion when there is no Fisher
   * information for a parameter i within the matrix. In that case the zero column and row is
   * removed from the matrix before inversion and the inverted matrix contains zeros.
   *
   * <p>The square root of the reciprocal of the diagonal element of the Fisher information matrix
   * is a (possibly loose) lower bound.
   *
   * @param d the input value
   * @return the reciprocal of the square root of the input value
   */
  public static double reciprocalSqrt(double d) {
    return (d > 0) ? 1.0 / Math.sqrt(d) : 0;
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
   * Gets a copy of the matrix.
   *
   * @return the matrix
   */
  public double[][] getSquareMatrix() {
    return EJMLLinearSolver.toSquareData(m);
  }

  /**
   * Gets a reference to the matrix.
   *
   * @return the matrix
   */
  public DenseMatrix64F getMatrix() {
    return m;
  }

  @Override
  public String toString() {
    return "CRLB=" + Arrays.toString(crlb) + "\nM=" + m.toString();
  }
}
