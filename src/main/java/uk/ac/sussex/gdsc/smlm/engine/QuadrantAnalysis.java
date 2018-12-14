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

package uk.ac.sussex.gdsc.smlm.engine;

import uk.ac.sussex.gdsc.core.utils.MathUtils;

/**
 * Performs quadrant analysis on fit residuals to look for asymmetry.
 */
public class QuadrantAnalysis {
  // Make these public for simplicity

  /** The sum of the 4 quadrants using the X dividing lines (two diagonals through the centre). */
  public double sumABCD;
  /** The sum of the A quadrant using the X dividing lines (two diagonals through the centre). */
  public double sumA;
  /** The sum of the B quadrant using the X dividing lines (two diagonals through the centre). */
  public double sumB;
  /** The sum of the C quadrant using the X dividing lines (two diagonals through the centre). */
  public double sumC;
  /** The sum of the D quadrant using the X dividing lines (two diagonals through the centre). */
  public double sumD;
  /**
   * The sum of the 4 quadrants using the + dividing lines (horizontal and vertical through the
   * centre).
   */
  public double sumABCD2;
  /**
   * The sum of the A quadrant using the + dividing lines (horizontal and vertical through the
   * centre).
   */
  public double sumA2;
  /**
   * The sum of the B quadrant using the + dividing lines (horizontal and vertical through the
   * centre).
   */
  public double sumB2;
  /**
   * The sum of the C quadrant using the + dividing lines (horizontal and vertical through the
   * centre).
   */
  public double sumC2;
  /**
   * The sum of the D quadrant using the + dividing lines (horizontal and vertical through the
   * centre).
   */
  public double sumD2;

  /**
   * {@link #sumA} + {@link #sumC}.
   */
  public double sumAC;
  /**
   * {@link #sumB} + {@link #sumD}.
   */
  public double sumBD;
  /**
   * The asymmetry score for the + dividing lines. Math.abs({@link #sumAC} - {@link #sumBD}) /
   * {@link #sumABCD}.
   */
  public double score1;

  /**
   * {@link #sumA2} + {@link #sumC2}.
   */
  public double sumAC2;
  /**
   * {@link #sumB2} + {@link #sumD2}.
   */
  public double sumBD2;
  /**
   * The asymmetry score for the + dividing lines. Math.abs({@link #sumAC2} - {@link #sumBD2}) /
   * {@link #sumABCD2}.
   */
  public double score2;

  /** The vector direction defined by the assymetry in the quadrants. */
  public int[] vector;
  /**
   * The maximum asymmetry score for quadrant analysis. The max of {@link #score1} and
   * {@link #score2}.
   */
  public double score;

  /**
   * Proposed x coordinate for centre 1 created by
   * {@link #computeDoubletCentres(int, int, int, int, double, double)}.
   */
  public double x1;
  /**
   * Proposed y coordinate for centre 1 created by
   * {@link #computeDoubletCentres(int, int, int, int, double, double)}.
   */
  public double y1;
  /**
   * Proposed x coordinate for centre 2 created by
   * {@link #computeDoubletCentres(int, int, int, int, double, double)}.
   */
  public double x2;
  /**
   * Proposed y coordinate for centre 2 created by
   * {@link #computeDoubletCentres(int, int, int, int, double, double)}.
   */
  public double y2;

  /**
   * Proposed integer x coordinate for centre 1 created by
   * {@link #computeDoubletCentres(int, int, int, int, double, double)}.
   */
  public int xi1;
  /**
   * Proposed integer y coordinate for centre 1 created by
   * {@link #computeDoubletCentres(int, int, int, int, double, double)}.
   */
  public int yi1;
  /**
   * Proposed integer x coordinate for centre 2 created by
   * {@link #computeDoubletCentres(int, int, int, int, double, double)}.
   */
  public int xi2;
  /**
   * Proposed integer y coordinate for centre 2 created by
   * {@link #computeDoubletCentres(int, int, int, int, double, double)}.
   */
  public int yi2;

  /**
   * Perform quadrant analysis as per rapidSTORM.
   *
   * <p>When two fluorophores emit close to each other, typically the nonlinear fit will result in a
   * suspected fluorophore position midway between the two fluorophores and with a high amplitude.
   * In this case, the fit results show a characteristic handle structure: The two true fluorophore
   * emissions leave slightly positive residues, while there are negative residues on an axis
   * perpendicular to the one connecting the fluorophores.
   *
   * <p>This condition is detected well by quadrant-differential residue analysis: The residue
   * matrix is divided into quadrants, with the pixels above both diagonals forming the upper
   * quadrant, the pixels above the main and below the off diagonal forming the right quadrants and
   * so on. Pixels right on the diagonals are ignored. Then, the opposing quadrants are summed, and
   * these sums substracted from another, resulting in two quadrant differences: upper and lower
   * minus right and left and right and left minus upper and lower. This process is repeated for the
   * quadrants defined by the central row and the central column.
   *
   * <p>The maximum sum obtained in this way divided by the sum of the absolute quadrant
   * contributions is an observable correlating highly with the true presence of double emitters.
   * Also, the quadrants containing the positive contribution in the highest sum indicate along
   * which axis the double emission happened.
   *
   * @param residuals the residuals
   * @param width the width
   * @param height the height
   * @param cx the centre in x
   * @param cy the centre in y
   * @return true, if successful
   */
  public boolean quadrantAnalysis(final double[] residuals, final int width, final int height,
      final int cx, final int cy) {
    vector = null;

    if (cx < 0 || cx >= width || cy < 0 || cy >= height) {
      return false;
    }

    // Compute quadrants

    // X quadrant:
    // .AAA.
    // D.A.B
    // DD.BB
    // D.C.B
    // .CCC.
    sumABCD = 0;
    sumA = 0;
    sumB = 0;
    sumC = 0;
    sumD = 0;
    for (int y = cy, xa = cx, xb = cx; y < height; y++, xa--, xb++) {
      for (int x = 0, index = y * width; x < width; x++, index++) {
        sumABCD += Math.abs(residuals[index]);
        if (x < xa) {
          sumD += residuals[index];
        } else if (x < xb && x > xa) {
          sumC += residuals[index];
        } else if (x > xb) {
          sumB += residuals[index];
        } else {
          sumABCD -= Math.abs(residuals[index]);
        }
      }
    }
    for (int y = cy - 1, xa = cx - 1, xb = cx + 1; y >= 0; y--, xa--, xb++) {
      for (int x = 0, index = y * width; x < width; x++, index++) {
        sumABCD += Math.abs(residuals[index]);
        if (x < xa) {
          sumD += residuals[index];
        } else if (x < xb && x > xa) {
          sumA += residuals[index];
        } else if (x > xb) {
          sumB += residuals[index];
        } else {
          sumABCD -= Math.abs(residuals[index]);
        }
      }
    }

    // Similar for + quadrants:
    // AA.BB
    // AA.BB
    // .....
    // DD.CC
    // DD.CC
    sumABCD2 = 0;
    sumA2 = 0;
    sumB2 = 0;
    sumC2 = 0;
    sumD2 = 0;
    for (int y = cy + 1; y < height; y++) {
      for (int x = 0, index = y * width; x < width; x++, index++) {
        sumABCD2 += Math.abs(residuals[index]);
        if (x < cx) {
          sumD2 += residuals[index];
        } else if (x > cx) {
          sumC2 += residuals[index];
        }
      }
    }
    for (int y = cy - 1; y >= 0; y--) {
      for (int x = 0, index = y * width; x < width; x++, index++) {
        sumABCD2 += Math.abs(residuals[index]);
        if (x < cx) {
          sumA2 += residuals[index];
        } else if (x > cx) {
          sumB2 += residuals[index];
        }
      }
    }

    // X quadrant:
    // .AAA.
    // D.A.B
    // DD.BB
    // D.C.B
    // .CCC.
    sumAC = sumA + sumC;
    sumBD = sumB + sumD;
    score1 = Math.abs(sumAC - sumBD) / sumABCD;

    // + quadrant:
    // AA.BB
    // AA.BB
    // .....
    // DD.CC
    // DD.CC
    sumAC2 = sumA2 + sumC2;
    sumBD2 = sumB2 + sumD2;
    score2 = Math.abs(sumAC2 - sumBD2) / sumABCD2;

    if (score1 > score2) {
      vector = (sumAC > sumBD) ? new int[] {0, 1} : new int[] {1, 0};
      score = score1;
    } else {
      vector = (sumAC2 > sumBD2) ? new int[] {1, 1} : new int[] {1, -1};
      score = score2;
    }

    return true;
  }

  /**
   * Locate the 2 new centres by moving out into the quadrant defined by the computed vector by the
   * defined shift.
   *
   * <p>Requires a valid call to {@link #quadrantAnalysis(double[], int, int, int, int)} to create
   * the vector
   *
   * @param width the width
   * @param height the height
   * @param cx the centre in x
   * @param cy the centre in y
   * @param shiftx the shiftx
   * @param shifty the shifty
   * @return true, if successful
   */
  public boolean computeDoubletCentres(final int width, final int height, final int cx,
      final int cy, double shiftx, double shifty) {
    if (vector == null || cx < 0 || cx >= width || cy < 0 || cy >= height) {
      return false;
    }

    // Pick double coords. The input centres must be shifted to the centre of the pixel by 0.5.
    x1 = cx + 0.5 + vector[0] * shiftx;
    y1 = cy + 0.5 + vector[1] * shifty;
    x2 = cx + 0.5 - vector[0] * shiftx;
    y2 = cy + 0.5 - vector[1] * shifty;

    final double maxWidth = width - 0.01;
    final double maxHeight = height - 0.01;

    // Check bounds
    x1 = MathUtils.clip(0, maxWidth, x1);
    y1 = MathUtils.clip(0, maxHeight, y1);
    x2 = MathUtils.clip(0, maxWidth, x2);
    y2 = MathUtils.clip(0, maxHeight, y2);

    xi1 = (int) (x1);
    yi1 = (int) (y1);
    xi2 = (int) (x2);
    yi2 = (int) (y2);

    // Check the two points are not the same pixel.
    // This is an edge case where the shift is small.
    if (xi1 == xi2 && yi1 == yi2) {
      // This can only happen when the shift is zero after rounding.
      // If they are the same then the value (xi1,yi1) should be cx,cy
      // and we can move along the vector.
      x1 = cx + 0.5 + vector[0];
      y1 = cy + 0.5 + vector[1];
      x2 = cx + 0.5 - vector[0];
      y2 = cy + 0.5 - vector[1];

      // Check bounds again
      x1 = MathUtils.clip(0, maxWidth, x1);
      y1 = MathUtils.clip(0, maxHeight, y1);
      x2 = MathUtils.clip(0, maxWidth, x2);
      y2 = MathUtils.clip(0, maxHeight, y2);

      xi1 = (int) (x1);
      yi1 = (int) (y1);
      xi2 = (int) (x2);
      yi2 = (int) (y2);

      return (xi1 != xi2 || yi1 != yi2);
    }

    return true;
  }

  /**
   * Gets the angle between two vectors.
   *
   * @param v1 the first vector
   * @param v2 the second vector
   * @return the angle (in radians)
   */
  public static double getAngle(int[] v1, double[] v2) {
    double d1 = (double) v1[0] * v1[0] + v1[1] * v1[1];
    double d2 = v2[0] * v2[0] + v2[1] * v2[1];
    if (d1 > 0.0 && d2 > 0.0) {
      d1 = Math.sqrt(d1);
      d2 = Math.sqrt(d2);
      final double sum = v1[0] * v2[0] + v1[1] * v2[1];
      final double cosang = sum / (d1 * d2);

      if (cosang > 1.0) {
        return 0;
      } else if (cosang < -1.0) {
        return Math.PI;
      }

      return Math.acos(cosang);
    }
    return 999;
  }

  /**
   * Gets the angle between two vectors.
   *
   * @param v1 the first vector
   * @param v2 the second vector
   * @return the angle (in radians)
   */
  public static double getAngle(double[] v1, double[] v2) {
    double d1 = v1[0] * v1[0] + v1[1] * v1[1];
    double d2 = v2[0] * v2[0] + v2[1] * v2[1];
    if (d1 > 0.0 && d2 > 0.0) {
      d1 = Math.sqrt(d1);
      d2 = Math.sqrt(d2);
      final double sum = v1[0] * v2[0] + v1[1] * v2[1];
      final double cosang = sum / (d1 * d2);

      if (cosang > 1.0) {
        return 0;
      } else if (cosang < -1.0) {
        return Math.PI;
      }

      return Math.acos(cosang);
    }
    return 999;
  }
}
