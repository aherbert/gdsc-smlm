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

/**
 * Solves (one) linear equation, a x = b.
 */
public class GaussJordan {
  private int maxRow;
  private int maxCol;

  // Allow use of parameter names a and b
  // @CHECKSTYLE.OFF: ParameterName

  private boolean findPivot(float[][] a, int[] piv) {
    float max = 0;

    for (int i = 0; i < piv.length; i++) {
      if (piv[i] != 1) {
        for (int j = 0; j < piv.length; j++) {
          if (piv[j] == 0) {
            if (Math.abs(a[i][j]) >= max) {
              max = Math.abs(a[i][j]);
              maxRow = i;
              maxCol = j;
            }
          } else if (piv[j] > 1) {
            // This should not happen, i.e. a second pivot around a column
            return false;
          }
        }
      }
    }

    piv[maxCol]++;

    return true;
  }

  private void interchangeRowsVector(float[][] a, float[] b) {
    for (int j = a[maxRow].length; j-- > 0;) {
      final float tmp = a[maxRow][j];
      a[maxRow][j] = a[maxCol][j];
      a[maxCol][j] = tmp;
    }

    final float tmp = b[maxRow];
    b[maxRow] = b[maxCol];
    b[maxCol] = tmp;
  }

  private boolean pivotVector(float[][] a, float[] b) {
    if (a[maxCol][maxCol] == 0) {
      return false;
    }

    final float piv_inv = 1 / a[maxCol][maxCol];

    a[maxCol][maxCol] = 1;

    for (int i = 0; i < a[maxCol].length; i++) {
      a[maxCol][i] *= piv_inv;
    }
    b[maxCol] *= piv_inv;

    for (int i = 0; i < a[maxCol].length; i++) {
      if (i != maxCol) {
        final float x = a[i][maxCol];
        a[i][maxCol] = 0;

        for (int j = 0; j < a[maxCol].length; j++) {
          a[i][j] -= x * a[maxCol][j];
        }

        b[i] -= x * b[maxCol];
      }
    }

    return true;
  }

  private static void unscrambleVector(float[][] a, int[] row, int[] col) {
    for (int j = row.length; j-- > 0;) {
      if (row[j] != col[j]) {
        for (int i = row.length; i-- > 0;) {
          final float tmp = a[i][row[j]];
          a[i][row[j]] = a[i][col[j]];
          a[i][col[j]] = tmp;
        }
      }
    }
  }

  /**
   * Solves (one) linear equation, a x = b, for x[n].
   *
   * <p>On input have a[n][n], b[n]. On output these replaced by a_inverse[n][n], x[n].
   *
   * @param a the a
   * @param b the b
   * @return False if the equation is singular (no solution)
   */
  public boolean solve(float[][] a, float[] b) {
    final int[] piv = new int[b.length];
    final int[] row = new int[b.length];
    final int[] col = new int[b.length];
    return solve(a, b, piv, row, col);
  }

  /**
   * Solves (one) linear equation, a x = b, for x[n].
   *
   * <p>On input have a[n][n], b[n]. On output these replaced by a_inverse[n][n], x[n].
   *
   * <p>piv[n], row[n], col[n] (all ints) are used for storage
   *
   * @param a the a
   * @param b the b
   * @param piv the pivot storage
   * @param row the row storage
   * @param col the column storage
   * @return False if the equation is singular (no solution)
   */
  private boolean solve(float[][] a, float[] b, int[] piv, int[] row, int[] col) {
    maxRow = 0;
    maxCol = 0;

    for (int i = 0; i < piv.length; i++) {
      piv[i] = 0;
    }

    for (int i = 0; i < piv.length; i++) {
      if (!findPivot(a, piv)) {
        return false;
      }

      if (maxRow != maxCol) {
        interchangeRowsVector(a, b);
      }

      row[i] = maxRow;
      col[i] = maxCol;

      if (!pivotVector(a, b)) {
        return false;
      }
    }

    unscrambleVector(a, row, col);
    return true;
  }

  // The above code is repeated for <double>
  // @CHECKSTYLE.OFF: OverloadMethodsDeclarationOrder

  private boolean findPivot(double[][] a, int[] piv) {
    double max = 0;

    for (int i = 0; i < piv.length; i++) {
      if (piv[i] != 1) {
        for (int j = 0; j < piv.length; j++) {
          if (piv[j] == 0) {
            if (Math.abs(a[i][j]) >= max) {
              max = Math.abs(a[i][j]);
              maxRow = i;
              maxCol = j;
            }
          } else if (piv[j] > 1) {
            // This should not happen, i.e. a second pivot around a column
            return false;
          }
        }
      }
    }

    piv[maxCol]++;

    return true;
  }

  private void interchangeRowsVector(double[][] a, double[] b) {
    for (int j = a[maxRow].length; j-- > 0;) {
      final double tmp = a[maxRow][j];
      a[maxRow][j] = a[maxCol][j];
      a[maxCol][j] = tmp;
    }

    final double tmp = b[maxRow];
    b[maxRow] = b[maxCol];
    b[maxCol] = tmp;
  }

  private boolean pivotVector(double[][] a, double[] b) {
    if (a[maxCol][maxCol] == 0) {
      return false;
    }

    final double piv_inv = 1 / a[maxCol][maxCol];

    a[maxCol][maxCol] = 1;

    for (int i = 0; i < a[maxCol].length; i++) {
      a[maxCol][i] *= piv_inv;
    }
    b[maxCol] *= piv_inv;

    for (int i = 0; i < a[maxCol].length; i++) {
      if (i != maxCol) {
        final double x = a[i][maxCol];
        a[i][maxCol] = 0;

        for (int j = 0; j < a[maxCol].length; j++) {
          a[i][j] -= x * a[maxCol][j];
        }

        b[i] -= x * b[maxCol];
      }
    }

    return true;
  }

  private static void unscrambleVector(double[][] a, int[] row, int[] col) {
    for (int j = row.length; j-- > 0;) {
      if (row[j] != col[j]) {
        for (int i = row.length; i-- > 0;) {
          final double tmp = a[i][row[j]];
          a[i][row[j]] = a[i][col[j]];
          a[i][col[j]] = tmp;
        }
      }
    }
  }

  /**
   * Solves (one) linear equation, a x = b, for x[n].
   *
   * <p>On input have a[n][n], b[n]. On output these replaced by a_inverse[n][n], x[n].
   *
   * @param a the a
   * @param b the b
   * @return False if the equation is singular (no solution)
   */
  public boolean solve(double[][] a, double[] b) {
    final int[] piv = new int[b.length];
    final int[] row = new int[b.length];
    final int[] col = new int[b.length];
    return solve(a, b, piv, row, col);
  }

  /**
   * Solves (one) linear equation, a x = b, for x[n].
   *
   * <p>On input have a[n][n], b[n]. On output these replaced by a_inverse[n][n], x[n].
   *
   * <p>piv[n], row[n], col[n] (all ints) are used for storage
   *
   * @param a the a
   * @param b the b
   * @param piv the pivot storage
   * @param row the row storage
   * @param col the column storage
   * @return False if the equation is singular (no solution)
   */
  public boolean solve(double[][] a, double[] b, int[] piv, int[] row, int[] col) {
    maxRow = 0;
    maxCol = 0;

    for (int i = 0; i < piv.length; i++) {
      piv[i] = 0;
    }

    for (int i = 0; i < piv.length; i++) {
      if (!findPivot(a, piv)) {
        return false;
      }

      if (maxRow != maxCol) {
        interchangeRowsVector(a, b);
      }

      row[i] = maxRow;
      col[i] = maxCol;

      if (!pivotVector(a, b)) {
        return false;
      }
    }

    unscrambleVector(a, row, col);
    return true;
  }
}
