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

package uk.ac.sussex.gdsc.smlm.ij.results;

import java.awt.Rectangle;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageSettings;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageSizeMode;

/**
 * Test the {@link ImagePeakResultsFactory} functionality.
 */
@SuppressWarnings({"javadoc"})
class ImagePeakResultsFactoryTest {
  @Test
  void testGetScaleImageSize() {
    assertGetScaleImageSize(200, 300, 0, 1);
    assertGetScaleImageSize(0, 0, 600, 1);
    assertGetScaleImageSize(200, 300, 600, 2);
    assertGetScaleImageSize(300, 200, 600, 2);
    // Rounding errors
    assertGetScaleImageSize(191, 1956);
    assertGetScaleImageSize(631, 796);
    assertGetScaleImageSize(193, 1892);
    assertGetScaleImageSize(149, 1541);
    assertGetScaleImageSize(191, 1956);
    assertGetScaleImageSize(883, 896);
    assertGetScaleImageSize(823, 842);
    assertGetScaleImageSize(911, 1848);
  }

  private static void assertGetScaleImageSize(int width, int height, int imageSize,
      float expected) {
    final ResultsImageSettings.Builder b = ResultsImageSettings.newBuilder()
        .setImageSizeMode(ResultsImageSizeMode.IMAGE_SIZE).setImageSize(imageSize);
    final Rectangle bounds = new Rectangle(width, height);
    Assertions.assertEquals(expected, ImagePeakResultsFactory.getScale(b, bounds, 0));
  }

  private static void assertGetScaleImageSize(int max, int imageSize) {
    final ResultsImageSettings.Builder b = ResultsImageSettings.newBuilder()
        .setImageSizeMode(ResultsImageSizeMode.IMAGE_SIZE).setImageSize(imageSize);
    final Rectangle bounds = new Rectangle(max, max);
    float scale = ImagePeakResultsFactory.getScale(b, bounds, 0);
    Assertions.assertEquals(imageSize, (int) Math.ceil(max * scale));
  }

  @Test
  void testGetScalePixelSize() {
    assertGetScalePixelSize(7.5f, 0, 1);
    assertGetScalePixelSize(0, 75, 1);
    assertGetScalePixelSize(5, 75, 15);
    assertGetScalePixelSize(5, 100, 20);
  }

  private static void assertGetScalePixelSize(double pixelSize, double nmPerPixel, float expected) {
    final ResultsImageSettings.Builder b = ResultsImageSettings.newBuilder()
        .setImageSizeMode(ResultsImageSizeMode.PIXEL_SIZE).setPixelSize(pixelSize);
    Assertions.assertEquals(expected, ImagePeakResultsFactory.getScale(b, null, nmPerPixel));
  }

  @Test
  void testGetScaleScaled() {
    assertGetScaleScaled(0, 1);
    assertGetScaleScaled(-1, 1);
    assertGetScaleScaled(7.5, 7.5f);
    assertGetScaleScaled(5, 5);
    assertGetScaleScaled(2, 2);
  }

  private static void assertGetScaleScaled(double scale, float expected) {
    // Check this is the default (no images size mode is set)
    final ResultsImageSettings.Builder b = ResultsImageSettings.newBuilder().setScale(scale);
    Assertions.assertEquals(expected, ImagePeakResultsFactory.getScale(b, null, 0));
    // Explicitly set the size mode
    b.setImageSizeMode(ResultsImageSizeMode.SCALED);
    Assertions.assertEquals(expected, ImagePeakResultsFactory.getScale(b, null, 0));
  }
}
