/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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

/*
 * Copyright (c) 2009-2012, Peter Abeles. All Rights Reserved.
 *
 * This file is part of Efficient Java Matrix Library (EJML).
 *
 * EJML is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
 * General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * EJML is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
 * General Public License for more det_recipails.
 *
 * You should have received a copy of the GNU Lesser General Public License along with EJML. If not,
 * see <http://www.gnu.org/licenses/>.
 */

import org.ejml.data.DenseMatrix64F;

/**
 * Alex Herbert created this extension of {@link org.ejml.alg.dense.misc.UnrolledInverseFromMinor}
 * to compute only the diagonal of the inverse.
 */
public class UnrolledInverseFromMinorExt {
  /** The maximum supported matrix size. */
  public static final int MAX = 5;

  /**
   * Invert the matrix and return the diagonal.
   *
   * @param mat the matrix
   * @return the diagonal of the inverse (or null if there is no det_reciperminant)
   */
  public static double[] inv(DenseMatrix64F mat) {
    double max = Math.abs(mat.data[0]);
    final int N = mat.getNumElements();

    for (int i = 1; i < N; i++) {
      final double a = Math.abs(mat.data[i]);
      if (a > max) {
        max = a;
      }
    }

    if (mat.numRows == 5) {
      return inv5(mat, 1.0 / max);
    } else if (mat.numRows == 4) {
      return inv4(mat, 1.0 / max);
    } else if (mat.numRows == 3) {
      return inv3(mat, 1.0 / max);
    } else if (mat.numRows == 2) {
      return inv2(mat, 1.0 / max);
    } else if (mat.numRows == 1) {
      return new double[] {1.0 / mat.data[0]};
    } else {
      throw new IllegalArgumentException("Not supported");
    }
  }

  /**
   * Invert the 2x2 matrix and return the diagonal.
   *
   * @param mat the 2x2 matrix
   * @param scale the scale (set to the maximum of the matrix data, or 1)
   * @return the diagonal of the inverse (or null if there is no det_reciperminant)
   */
  public static double[] inv2(DenseMatrix64F mat, double scale) {
    final double[] data = mat.data;

    final double a11 = data[0] * scale;
    final double a12 = data[1] * scale;
    final double a21 = data[2] * scale;
    final double a22 = data[3] * scale;

    final double m11 = a22;
    final double m12 = -(a21);

    final double det_recip = scale / (a11 * m11 + a12 * m12);
    if (!Double.isFinite(det_recip)) {
      return null;
    }

    final double m22 = a11;

    return new double[] {m11 * det_recip, m22 * det_recip};
  }

  /**
   * Invert the 3x3 matrix and return the diagonal.
   *
   * @param mat the 3x3 matrix
   * @param scale the scale (set to the maximum of the matrix data, or 1)
   * @return the diagonal of the inverse (or null if there is no det_reciperminant)
   */
  public static double[] inv3(DenseMatrix64F mat, double scale) {
    final double[] data = mat.data;

    final double a11 = data[0] * scale;
    final double a12 = data[1] * scale;
    final double a13 = data[2] * scale;
    final double a21 = data[3] * scale;
    final double a22 = data[4] * scale;
    final double a23 = data[5] * scale;
    final double a31 = data[6] * scale;
    final double a32 = data[7] * scale;
    final double a33 = data[8] * scale;

    final double m11 = a22 * a33 - a23 * a32;
    final double m12 = -(a21 * a33 - a23 * a31);
    final double m13 = a21 * a32 - a22 * a31;

    final double det_recip = scale / (a11 * m11 + a12 * m12 + a13 * m13);
    if (!Double.isFinite(det_recip)) {
      return null;
    }

    final double m22 = a11 * a33 - a13 * a31;
    final double m33 = a11 * a22 - a12 * a21;

    return new double[] {m11 * det_recip, m22 * det_recip, m33 * det_recip};
  }

  /**
   * Invert the 4x4 matrix and return the diagonal.
   *
   * @param mat the 4x4 matrix
   * @param scale the scale (set to the maximum of the matrix data, or 1)
   * @return the diagonal of the inverse (or null if there is no det_reciperminant)
   */
  public static double[] inv4(DenseMatrix64F mat, double scale) {
    final double[] data = mat.data;

    final double a11 = data[0] * scale;
    final double a12 = data[1] * scale;
    final double a13 = data[2] * scale;
    final double a14 = data[3] * scale;
    final double a21 = data[4] * scale;
    final double a22 = data[5] * scale;
    final double a23 = data[6] * scale;
    final double a24 = data[7] * scale;
    final double a31 = data[8] * scale;
    final double a32 = data[9] * scale;
    final double a33 = data[10] * scale;
    final double a34 = data[11] * scale;
    final double a41 = data[12] * scale;
    final double a42 = data[13] * scale;
    final double a43 = data[14] * scale;
    final double a44 = data[15] * scale;

    final double m11 = +a22 * (a33 * a44 - a34 * a43) - a23 * (a32 * a44 - a34 * a42)
        + a24 * (a32 * a43 - a33 * a42);
    final double m12 = -(+a21 * (a33 * a44 - a34 * a43) - a23 * (a31 * a44 - a34 * a41)
        + a24 * (a31 * a43 - a33 * a41));
    final double m13 = +a21 * (a32 * a44 - a34 * a42) - a22 * (a31 * a44 - a34 * a41)
        + a24 * (a31 * a42 - a32 * a41);
    final double m14 = -(+a21 * (a32 * a43 - a33 * a42) - a22 * (a31 * a43 - a33 * a41)
        + a23 * (a31 * a42 - a32 * a41));

    final double det_recip = scale / (a11 * m11 + a12 * m12 + a13 * m13 + a14 * m14);
    if (!Double.isFinite(det_recip)) {
      return null;
    }

    final double m22 = +a11 * (a33 * a44 - a34 * a43) - a13 * (a31 * a44 - a34 * a41)
        + a14 * (a31 * a43 - a33 * a41);
    final double m33 = +a11 * (a22 * a44 - a24 * a42) - a12 * (a21 * a44 - a24 * a41)
        + a14 * (a21 * a42 - a22 * a41);
    final double m44 = +a11 * (a22 * a33 - a23 * a32) - a12 * (a21 * a33 - a23 * a31)
        + a13 * (a21 * a32 - a22 * a31);

    return new double[] {m11 * det_recip, m22 * det_recip, m33 * det_recip, m44 * det_recip};
  }

  /**
   * Invert the 5x5 matrix and return the diagonal.
   *
   * @param mat the 5x5 matrix
   * @param scale the scale (set to the maximum of the matrix data, or 1)
   * @return the diagonal of the inverse (or null if there is no det_reciperminant)
   */
  public static double[] inv5(DenseMatrix64F mat, double scale) {
    final double[] data = mat.data;

    final double a11 = data[0] * scale;
    final double a12 = data[1] * scale;
    final double a13 = data[2] * scale;
    final double a14 = data[3] * scale;
    final double a15 = data[4] * scale;
    final double a21 = data[5] * scale;
    final double a22 = data[6] * scale;
    final double a23 = data[7] * scale;
    final double a24 = data[8] * scale;
    final double a25 = data[9] * scale;
    final double a31 = data[10] * scale;
    final double a32 = data[11] * scale;
    final double a33 = data[12] * scale;
    final double a34 = data[13] * scale;
    final double a35 = data[14] * scale;
    final double a41 = data[15] * scale;
    final double a42 = data[16] * scale;
    final double a43 = data[17] * scale;
    final double a44 = data[18] * scale;
    final double a45 = data[19] * scale;
    final double a51 = data[20] * scale;
    final double a52 = data[21] * scale;
    final double a53 = data[22] * scale;
    final double a54 = data[23] * scale;
    final double a55 = data[24] * scale;

    final double m11 = +a22
        * (+a33 * (a44 * a55 - a45 * a54) - a34 * (a43 * a55 - a45 * a53)
            + a35 * (a43 * a54 - a44 * a53))
        - a23 * (+a32 * (a44 * a55 - a45 * a54) - a34 * (a42 * a55 - a45 * a52)
            + a35 * (a42 * a54 - a44 * a52))
        + a24 * (+a32 * (a43 * a55 - a45 * a53) - a33 * (a42 * a55 - a45 * a52)
            + a35 * (a42 * a53 - a43 * a52))
        - a25 * (+a32 * (a43 * a54 - a44 * a53) - a33 * (a42 * a54 - a44 * a52)
            + a34 * (a42 * a53 - a43 * a52));
    final double m12 = -(+a21
        * (+a33 * (a44 * a55 - a45 * a54) - a34 * (a43 * a55 - a45 * a53)
            + a35 * (a43 * a54 - a44 * a53))
        - a23 * (+a31 * (a44 * a55 - a45 * a54) - a34 * (a41 * a55 - a45 * a51)
            + a35 * (a41 * a54 - a44 * a51))
        + a24 * (+a31 * (a43 * a55 - a45 * a53) - a33 * (a41 * a55 - a45 * a51)
            + a35 * (a41 * a53 - a43 * a51))
        - a25 * (+a31 * (a43 * a54 - a44 * a53) - a33 * (a41 * a54 - a44 * a51)
            + a34 * (a41 * a53 - a43 * a51)));
    final double m13 = +a21
        * (+a32 * (a44 * a55 - a45 * a54) - a34 * (a42 * a55 - a45 * a52)
            + a35 * (a42 * a54 - a44 * a52))
        - a22 * (+a31 * (a44 * a55 - a45 * a54) - a34 * (a41 * a55 - a45 * a51)
            + a35 * (a41 * a54 - a44 * a51))
        + a24 * (+a31 * (a42 * a55 - a45 * a52) - a32 * (a41 * a55 - a45 * a51)
            + a35 * (a41 * a52 - a42 * a51))
        - a25 * (+a31 * (a42 * a54 - a44 * a52) - a32 * (a41 * a54 - a44 * a51)
            + a34 * (a41 * a52 - a42 * a51));
    final double m14 = -(+a21
        * (+a32 * (a43 * a55 - a45 * a53) - a33 * (a42 * a55 - a45 * a52)
            + a35 * (a42 * a53 - a43 * a52))
        - a22 * (+a31 * (a43 * a55 - a45 * a53) - a33 * (a41 * a55 - a45 * a51)
            + a35 * (a41 * a53 - a43 * a51))
        + a23 * (+a31 * (a42 * a55 - a45 * a52) - a32 * (a41 * a55 - a45 * a51)
            + a35 * (a41 * a52 - a42 * a51))
        - a25 * (+a31 * (a42 * a53 - a43 * a52) - a32 * (a41 * a53 - a43 * a51)
            + a33 * (a41 * a52 - a42 * a51)));
    final double m15 = +a21
        * (+a32 * (a43 * a54 - a44 * a53) - a33 * (a42 * a54 - a44 * a52)
            + a34 * (a42 * a53 - a43 * a52))
        - a22 * (+a31 * (a43 * a54 - a44 * a53) - a33 * (a41 * a54 - a44 * a51)
            + a34 * (a41 * a53 - a43 * a51))
        + a23 * (+a31 * (a42 * a54 - a44 * a52) - a32 * (a41 * a54 - a44 * a51)
            + a34 * (a41 * a52 - a42 * a51))
        - a24 * (+a31 * (a42 * a53 - a43 * a52) - a32 * (a41 * a53 - a43 * a51)
            + a33 * (a41 * a52 - a42 * a51));

    final double det_recip = scale / (a11 * m11 + a12 * m12 + a13 * m13 + a14 * m14 + a15 * m15);
    if (!Double.isFinite(det_recip)) {
      return null;
    }

    final double m22 = +a11
        * (+a33 * (a44 * a55 - a45 * a54) - a34 * (a43 * a55 - a45 * a53)
            + a35 * (a43 * a54 - a44 * a53))
        - a13 * (+a31 * (a44 * a55 - a45 * a54) - a34 * (a41 * a55 - a45 * a51)
            + a35 * (a41 * a54 - a44 * a51))
        + a14 * (+a31 * (a43 * a55 - a45 * a53) - a33 * (a41 * a55 - a45 * a51)
            + a35 * (a41 * a53 - a43 * a51))
        - a15 * (+a31 * (a43 * a54 - a44 * a53) - a33 * (a41 * a54 - a44 * a51)
            + a34 * (a41 * a53 - a43 * a51));
    final double m33 = +a11
        * (+a22 * (a44 * a55 - a45 * a54) - a24 * (a42 * a55 - a45 * a52)
            + a25 * (a42 * a54 - a44 * a52))
        - a12 * (+a21 * (a44 * a55 - a45 * a54) - a24 * (a41 * a55 - a45 * a51)
            + a25 * (a41 * a54 - a44 * a51))
        + a14 * (+a21 * (a42 * a55 - a45 * a52) - a22 * (a41 * a55 - a45 * a51)
            + a25 * (a41 * a52 - a42 * a51))
        - a15 * (+a21 * (a42 * a54 - a44 * a52) - a22 * (a41 * a54 - a44 * a51)
            + a24 * (a41 * a52 - a42 * a51));
    final double m44 = +a11
        * (+a22 * (a33 * a55 - a35 * a53) - a23 * (a32 * a55 - a35 * a52)
            + a25 * (a32 * a53 - a33 * a52))
        - a12 * (+a21 * (a33 * a55 - a35 * a53) - a23 * (a31 * a55 - a35 * a51)
            + a25 * (a31 * a53 - a33 * a51))
        + a13 * (+a21 * (a32 * a55 - a35 * a52) - a22 * (a31 * a55 - a35 * a51)
            + a25 * (a31 * a52 - a32 * a51))
        - a15 * (+a21 * (a32 * a53 - a33 * a52) - a22 * (a31 * a53 - a33 * a51)
            + a23 * (a31 * a52 - a32 * a51));
    final double m55 = +a11
        * (+a22 * (a33 * a44 - a34 * a43) - a23 * (a32 * a44 - a34 * a42)
            + a24 * (a32 * a43 - a33 * a42))
        - a12 * (+a21 * (a33 * a44 - a34 * a43) - a23 * (a31 * a44 - a34 * a41)
            + a24 * (a31 * a43 - a33 * a41))
        + a13 * (+a21 * (a32 * a44 - a34 * a42) - a22 * (a31 * a44 - a34 * a41)
            + a24 * (a31 * a42 - a32 * a41))
        - a14 * (+a21 * (a32 * a43 - a33 * a42) - a22 * (a31 * a43 - a33 * a41)
            + a23 * (a31 * a42 - a32 * a41));

    return new double[] {m11 * det_recip, m22 * det_recip, m33 * det_recip, m44 * det_recip,
        m55 * det_recip};
  }
}
