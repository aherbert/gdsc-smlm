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

package uk.ac.sussex.gdsc.smlm.engine;

import java.util.Arrays;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.utils.MathUtils;

@SuppressWarnings({"javadoc"})
class QuadrantAnalysisTest {
  @Test
  void testQuadrantAnalysisEdgeCases() {
    final QuadrantAnalysis qa = new QuadrantAnalysis();
    final double[] residuals = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    final int width = 3;
    final int height = 3;
    final int cx = 1;
    final int cy = 1;
    Assertions.assertTrue(qa.quadrantAnalysis(residuals, width, height, cx, cy));
    Assertions.assertNotNull(qa.vector);
    Assertions.assertFalse(qa.quadrantAnalysis(residuals, width, height, -1, cy));
    Assertions.assertNull(qa.vector);
    Assertions.assertFalse(qa.quadrantAnalysis(residuals, width, height, cx, -1));
    Assertions.assertFalse(qa.quadrantAnalysis(residuals, width, height, width, cy));
    Assertions.assertFalse(qa.quadrantAnalysis(residuals, width, height, cx, height));
  }

  @Test
  void testQuadrantAnalysisEmpty() {
    final QuadrantAnalysis qa = new QuadrantAnalysis();
    final int width = 7;
    final int height = 5;
    final int cx = 3;
    final int cy = 2;
    Assertions.assertTrue(qa.quadrantAnalysis(new double[width * height], width, height, cx, cy));
    assertQa(qa, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, new int[] {1, -1});
  }

  @Test
  void testQuadrantAnalysis3x3() {
    final QuadrantAnalysis qa = new QuadrantAnalysis();
    //@formatter:off
    final double[] residuals = {
        2, 1, 3,
       -1,-2, 1,
       -5,-3, 1};
    //@formatter:on
    final int width = 3;
    final int height = 3;
    final int cx = 1;
    final int cy = 1;
    Assertions.assertTrue(qa.quadrantAnalysis(residuals, width, height, cx, cy));
    final double sum = Arrays.stream(residuals).map(Math::abs).sum();
    final double sumAbcd = sum - (2 + 2 + 1 + 5 + 3);
    final double sumA = 1;
    final double sumB = 1;
    final double sumC = -3;
    final double sumD = -1;
    final double sumAbcd2 = sum - (1 + 2 + 1 + 1 + 3);
    final double sumA2 = 2;
    final double sumB2 = 3;
    final double sumC2 = 1;
    final double sumD2 = -5;
    assertQa(qa, sumAbcd, sumA, sumB, sumC, sumD, sumAbcd2, sumA2, sumB2, sumC2, sumD2,
        new int[] {1, 1});
  }

  @Test
  void testQuadrantAnalysis7x5() {
    final QuadrantAnalysis qa = new QuadrantAnalysis();
    //@formatter:off
    final double[] residuals = {
        2, 1, 3, 1, 4,-1, 4,
       -1,-2, 1, 1, 2, 0, 1,
       -5,-3, 1, 7, 1,-1,-3,
        2, 2,-1,-1, 2,-2,-2,
        1,-1, 1, 4, 3, 1,-1};
    //@formatter:on
    final int width = 7;
    final int height = 5;
    final int cx = 3;
    final int cy = 2;
    Assertions.assertTrue(qa.quadrantAnalysis(residuals, width, height, cx, cy));
    final double sum = Arrays.stream(residuals).map(Math::abs).sum();
    final double sumAbcd = sum - (1 + 1 + 7 + 2 + 1 + 1 + 1 + 2 + 1);
    final double sumA = 3 + 1 + 4 + 1;
    final double sumB = 4 + 1 - 3 - 2 - 1 + 0 - 1 - 2 + 1;
    final double sumC = 3 + 4 + 1 - 1;
    final double sumD = 2 - 1 - 5 + 2 + 1 - 2 - 3 + 2 + 1;
    final double sumAbcd2 = sum - (5 + 3 + 1 + 7 + 1 + 1 + 3 + 1 + 1 + 1 + 4);
    final double sumA2 = 2 + 1 + 3 - 1 - 2 + 1;
    final double sumB2 = 4 - 1 + 4 + 2 + 0 + 1;
    final double sumC2 = 2 - 2 - 2 + 3 + 1 - 1;
    final double sumD2 = 2 + 2 - 1 + 1 - 1 + 1;
    assertQa(qa, sumAbcd, sumA, sumB, sumC, sumD, sumAbcd2, sumA2, sumB2, sumC2, sumD2,
        new int[] {0, 1});
  }

  @Test
  void testQuadrantAnalysis3x3a() {
    final QuadrantAnalysis qa = new QuadrantAnalysis();
    //@formatter:off
    final double[] residuals = {
        0, 1, 0,
        0, 0, 0,
        0, 2, 0};
    //@formatter:on
    final int width = 3;
    final int height = 3;
    final int cx = 1;
    final int cy = 1;
    Assertions.assertTrue(qa.quadrantAnalysis(residuals, width, height, cx, cy));
    final double sum = Arrays.stream(residuals).map(Math::abs).sum();
    final double sumAbcd = sum;
    final double sumA = 1;
    final double sumB = 0;
    final double sumC = 2;
    final double sumD = 0;
    final double sumAbcd2 = 0;
    final double sumA2 = 0;
    final double sumB2 = 0;
    final double sumC2 = 0;
    final double sumD2 = 0;
    assertQa(qa, sumAbcd, sumA, sumB, sumC, sumD, sumAbcd2, sumA2, sumB2, sumC2, sumD2,
        new int[] {0, 1});
  }

  @Test
  void testQuadrantAnalysis3x3b() {
    final QuadrantAnalysis qa = new QuadrantAnalysis();
    //@formatter:off
    final double[] residuals = {
        0,-1, 0,
        0, 0, 0,
        0,-2, 0};
    //@formatter:on
    final int width = 3;
    final int height = 3;
    final int cx = 1;
    final int cy = 1;
    Assertions.assertTrue(qa.quadrantAnalysis(residuals, width, height, cx, cy));
    final double sum = Arrays.stream(residuals).map(Math::abs).sum();
    final double sumAbcd = sum;
    final double sumA = -1;
    final double sumB = 0;
    final double sumC = -2;
    final double sumD = 0;
    final double sumAbcd2 = 0;
    final double sumA2 = 0;
    final double sumB2 = 0;
    final double sumC2 = 0;
    final double sumD2 = 0;
    assertQa(qa, sumAbcd, sumA, sumB, sumC, sumD, sumAbcd2, sumA2, sumB2, sumC2, sumD2,
        new int[] {1, 0});
  }

  @Test
  void testQuadrantAnalysis3x3c() {
    final QuadrantAnalysis qa = new QuadrantAnalysis();
    //@formatter:off
    final double[] residuals = {
        1, 0, 0,
        0, 0, 0,
        0, 0, 2};
    //@formatter:on
    final int width = 3;
    final int height = 3;
    final int cx = 1;
    final int cy = 1;
    Assertions.assertTrue(qa.quadrantAnalysis(residuals, width, height, cx, cy));
    final double sum = Arrays.stream(residuals).map(Math::abs).sum();
    final double sumAbcd = 0;
    final double sumA = 0;
    final double sumB = 0;
    final double sumC = 0;
    final double sumD = 0;
    final double sumAbcd2 = sum;
    final double sumA2 = 1;
    final double sumB2 = 0;
    final double sumC2 = 2;
    final double sumD2 = 0;
    assertQa(qa, sumAbcd, sumA, sumB, sumC, sumD, sumAbcd2, sumA2, sumB2, sumC2, sumD2,
        new int[] {1, 1});
  }

  @Test
  void testQuadrantAnalysis3x3d() {
    final QuadrantAnalysis qa = new QuadrantAnalysis();
    //@formatter:off
    final double[] residuals = {
       -1, 0, 0,
        0, 0, 0,
        0, 0,-2};
    //@formatter:on
    final int width = 3;
    final int height = 3;
    final int cx = 1;
    final int cy = 1;
    Assertions.assertTrue(qa.quadrantAnalysis(residuals, width, height, cx, cy));
    final double sum = Arrays.stream(residuals).map(Math::abs).sum();
    final double sumAbcd = 0;
    final double sumA = 0;
    final double sumB = 0;
    final double sumC = 0;
    final double sumD = 0;
    final double sumAbcd2 = sum;
    final double sumA2 = -1;
    final double sumB2 = 0;
    final double sumC2 = -2;
    final double sumD2 = 0;
    assertQa(qa, sumAbcd, sumA, sumB, sumC, sumD, sumAbcd2, sumA2, sumB2, sumC2, sumD2,
        new int[] {1, -1});
  }

  private static void assertQa(QuadrantAnalysis qa, double sumAbcd, double sumA, double sumB,
      double sumC, double sumD, double sumAbcd2, double sumA2, double sumB2, double sumC2,
      double sumD2, int[] vector) {
    Assertions.assertEquals(sumAbcd, qa.sumAbcd);
    Assertions.assertEquals(sumA, qa.sumA);
    Assertions.assertEquals(sumB, qa.sumB);
    Assertions.assertEquals(sumC, qa.sumC);
    Assertions.assertEquals(sumD, qa.sumD);
    Assertions.assertEquals(sumAbcd2, qa.sumAbcd2);
    Assertions.assertEquals(sumA2, qa.sumA2);
    Assertions.assertEquals(sumB2, qa.sumB2);
    Assertions.assertEquals(sumC2, qa.sumC2);
    Assertions.assertEquals(sumD2, qa.sumD2);
    Assertions.assertEquals(qa.sumAc, qa.sumA + qa.sumC);
    Assertions.assertEquals(qa.sumBd, qa.sumB + qa.sumD);
    Assertions.assertEquals(qa.sumAc2, qa.sumA2 + qa.sumC2);
    Assertions.assertEquals(qa.sumBd2, qa.sumB2 + qa.sumD2);
    Assertions.assertEquals(MathUtils.div0(Math.abs(sumA + sumC - sumB - sumD), sumAbcd),
        qa.score1);
    Assertions.assertEquals(MathUtils.div0(Math.abs(sumA2 + sumC2 - sumB2 - sumD2), sumAbcd2),
        qa.score2);
    Assertions.assertEquals(Math.max(qa.score1, qa.score2), qa.score);
    Assertions.assertArrayEquals(vector, qa.vector);
  }

  @Test
  void testComputeDoubleCentresEdgeCases() {
    final QuadrantAnalysis qa = new QuadrantAnalysis();
    final int width = 7;
    final int height = 5;
    final int cx = 3;
    final int cy = 2;
    final double shiftx = 0.5;
    final double shifty = 0.75;
    Assertions.assertFalse(qa.computeDoubletCentres(width, height, cx, cy, shiftx, shifty));
    // Set fake vector
    qa.vector = new int[] {0, 1};
    Assertions.assertTrue(qa.computeDoubletCentres(width, height, cx, cy, shiftx, shifty));
    Assertions.assertFalse(qa.computeDoubletCentres(width, height, -1, cy, shiftx, shifty));
    Assertions.assertFalse(qa.computeDoubletCentres(width, height, width, cy, shiftx, shifty));
    Assertions.assertFalse(qa.computeDoubletCentres(width, height, cx, -1, shiftx, shifty));
    Assertions.assertFalse(qa.computeDoubletCentres(width, height, cx, height, shiftx, shifty));
  }

  @Test
  void testComputeDoubleCentres() {
    assertComputeDoubleCentres(0, 1, 7, 5, 3, 2, 1.25, 1, 3.5, 3.5, 3.5, 1.5, 3, 3, 3, 1);
    assertComputeDoubleCentres(0, 1, 7, 5, 3, 2, 1.25, 2.75, 3.5, 4.99, 3.5, 0, 3, 4, 3, 0);
    assertComputeDoubleCentres(0, -1, 7, 5, 3, 2, 1.25, 1, 3.5, 1.5, 3.5, 3.5, 3, 1, 3, 3);
    assertComputeDoubleCentres(0, -1, 7, 5, 3, 2, 1.25, 2.75, 3.5, 0, 3.5, 4.99, 3, 0, 3, 4);
    assertComputeDoubleCentres(1, 0, 7, 5, 3, 2, 1.25, 1, 4.75, 2.5, 2.25, 2.5, 4, 2, 2, 2);
    assertComputeDoubleCentres(1, 0, 7, 5, 3, 2, 3.75, 1, 6.99, 2.5, 0, 2.5, 6, 2, 0, 2);
    assertComputeDoubleCentres(-1, 0, 7, 5, 3, 2, 1.25, 1, 2.25, 2.5, 4.75, 2.5, 2, 2, 4, 2);
    assertComputeDoubleCentres(-1, 0, 7, 5, 3, 2, 3.75, 1, 0, 2.5, 6.99, 2.5, 0, 2, 6, 2);
    // Small shift should be moved to a different pixel using a shift of 1 along the vector
    assertComputeDoubleCentres(0, 1, 7, 5, 3, 2, 0.25, 0.1, 3.5, 3.5, 3.5, 1.5, 3, 3, 3, 1);
    assertComputeDoubleCentres(0, -1, 7, 5, 3, 2, 0.25, 0.1, 3.5, 1.5, 3.5, 3.5, 3, 1, 3, 3);
    assertComputeDoubleCentres(1, 0, 7, 5, 3, 2, 0.25, 0.1, 4.5, 2.5, 2.5, 2.5, 4, 2, 2, 2);
    assertComputeDoubleCentres(-1, 0, 7, 5, 3, 2, 0.25, 0.1, 2.5, 2.5, 4.5, 2.5, 2, 2, 4, 2);
    // No vector results in colocated centres.
    assertComputeDoubleCentres(false, 0, 0, 7, 5, 3, 2, 0.25, 0.1, 3.5, 2.5, 3.5, 2.5, 3, 2, 3, 2);
  }

  private static void assertComputeDoubleCentres(int vectorx, int vectory, int width, int height,
      int cx, int cy, double shiftx, double shifty, double x1, double y1, double x2, double y2,
      int xi1, int yi1, int xi2, int yi2) {
    assertComputeDoubleCentres(true, vectorx, vectory, width, height, cx, cy, shiftx, shifty, x1,
        y1, x2, y2, xi1, yi1, xi2, yi2);
  }

  private static void assertComputeDoubleCentres(boolean result, int vectorx, int vectory,
      int width, int height, int cx, int cy, double shiftx, double shifty, double x1, double y1,
      double x2, double y2, int xi1, int yi1, int xi2, int yi2) {
    final QuadrantAnalysis qa = new QuadrantAnalysis();
    qa.vector = new int[] {vectorx, vectory};
    Assertions.assertEquals(result,
        qa.computeDoubletCentres(width, height, cx, cy, shiftx, shifty));
    Assertions.assertEquals(x1, qa.x1);
    Assertions.assertEquals(y1, qa.y1);
    Assertions.assertEquals(x2, qa.x2);
    Assertions.assertEquals(y2, qa.y2);
    Assertions.assertEquals(xi1, qa.xi1);
    Assertions.assertEquals(yi1, qa.yi1);
    Assertions.assertEquals(xi2, qa.xi2);
    Assertions.assertEquals(yi2, qa.yi2);
  }

  @Test
  void canGetAngle() {
    Assertions.assertEquals(0, QuadrantAnalysis.getAngle(new int[] {0, 1}, new double[] {0, 1}));
    Assertions.assertEquals(Math.PI,
        QuadrantAnalysis.getAngle(new int[] {0, 1}, new double[] {0, -1}));
    Assertions.assertEquals(999.0,
        QuadrantAnalysis.getAngle(new int[] {0, 0}, new double[] {0, 1}));
    Assertions.assertEquals(999.0,
        QuadrantAnalysis.getAngle(new int[] {0, 1}, new double[] {0, 0}));

    Assertions.assertEquals(Math.toRadians(90),
        QuadrantAnalysis.getAngle(new int[] {0, 1}, new double[] {1, 0}), 1e-10);
    Assertions.assertEquals(Math.toRadians(90),
        QuadrantAnalysis.getAngle(new int[] {0, 1}, new double[] {-1, 0}), 1e-10);
    Assertions.assertEquals(Math.toRadians(45),
        QuadrantAnalysis.getAngle(new int[] {0, 1}, new double[] {1, 1}), 1e-10);
    Assertions.assertEquals(Math.toRadians(135),
        QuadrantAnalysis.getAngle(new int[] {0, 1}, new double[] {1, -1}), 1e-10);
  }
}
