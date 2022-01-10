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
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageMode;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageSettings;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageSettingsOrBuilder;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageSizeMode;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageType;

/**
 * A factory for creating ImagePeakResults objects.
 */
public class ImagePeakResultsFactory {

  /**
   * Create a PeakResults image using the specified parameters.
   *
   * @param resultsImage The type of image
   * @param weighted Flag to indicate is the values should be bilinearly weighted on surrounding 4
   *        pixels (applied when plotting localisations at a single point)
   * @param equalised Flag to indicate if the image should have histogram equalisation applied
   * @param title The title of the image
   * @param bounds Define the bounding rectangle of the result coordinates
   * @param nmPerPixel The results scale in nanometers per pixel
   * @param imageScale Define the scale of the image relative to the bounding rectangle
   * @param precision For average precision plots this parameter specifies the fixed width of the
   *        PSF (in nm). If less than zero then defaults to the nmPerPixel value.
   * @param mode The mode for showing consecutive results in the same pixel location
   * @return The PeakResults image
   * @see #createPeakResultsImage(ResultsImageSettingsOrBuilder, String, Rectangle, double)
   */
  public static ImageJImagePeakResults createPeakResultsImage(ResultsImageType resultsImage,
      boolean weighted, boolean equalised, String title, Rectangle bounds, double nmPerPixel,
      double imageScale, double precision, ResultsImageMode mode) {
    // Delegate to the method using ResultsImageSettings
    final ResultsImageSettings.Builder builder =
        ResultsImageSettings.newBuilder().setImageType(resultsImage).setWeighted(weighted)
            .setEqualised(equalised).setScale(imageScale).setAveragePrecision(precision)
            .setImageMode(mode).setImageSizeMode(ResultsImageSizeMode.SCALED);
    final ImageJImagePeakResults image = createPeakResultsImage(builder, title, bounds, nmPerPixel);
    // For backwards compatibility the calibration is set using the nm/pixel argument
    image.setCalibration(nmPerPixel, 1);
    return image;
  }

  /**
   * Create a PeakResults image using the specified parameters.
   *
   * <p>The image scale (coordinate-to-pixel mapping) is determined based on the configured
   * settings. If the image scale cannot be set then the output image will be rendered using a
   * 1-to-1 coordinates-to-pixel mapping. Note: The nm-per-pixel calibration of the results is
   * required to correctly scale the image when rendering a fixed size output pixel. If not using a
   * fixed size output pixel then this parameter is ignored.
   *
   * <p>A correct calibration can be added to the returned image using
   * {@link uk.ac.sussex.gdsc.smlm.results.PeakResults#setCalibration(uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration)
   * PeakResults.setCalibration} before calling the
   * {@link uk.ac.sussex.gdsc.smlm.results.PeakResults#begin()} method. Alternatively use the
   * {@link ImageJImagePeakResults#setUncalibrated(boolean)} method to mark the image as
   * uncalibrated.
   *
   * <p>No checks are made on the dimensions of the image. The {@link ImageJImagePeakResults} will
   * throw an exception in the {@code begin()} method if the output image is too large to fit in a
   * 2D array or there is not enough memory.
   *
   * @param resultsSettings The results image settings
   * @param title The title of the image
   * @param bounds Define the bounding rectangle of the result coordinates
   * @param nmPerPixel The results scale in nanometers per pixel (ignored unless image type uses
   *        fixed size pixels)
   * @return The PeakResults image
   */
  public static ImageJImagePeakResults createPeakResultsImage(
      ResultsImageSettingsOrBuilder resultsSettings, String title, Rectangle bounds,
      double nmPerPixel) {
    final ResultsImageType resultsImage = resultsSettings.getImageType();
    final float imageScale = getScale(resultsSettings, bounds, nmPerPixel);
    ImageJImagePeakResults image;
    switch (resultsImage) {
      case DRAW_FITTED_PSF:
      case DRAW_INTENSITY_PRECISION:
      case DRAW_LOCALISATIONS_PRECISION:
      case DRAW_INTENSITY_AVERAGE_PRECISION:
      case DRAW_LOCALISATIONS_AVERAGE_PRECISION:
        // Special case for full PSF image
        final PsfImagePeakResults image2 = new PsfImagePeakResults(title, bounds, imageScale);
        if (resultsImage == ResultsImageType.DRAW_INTENSITY_AVERAGE_PRECISION
            || resultsImage == ResultsImageType.DRAW_LOCALISATIONS_AVERAGE_PRECISION) {
          // Fixed width display (in pixels)
          image2.setWidth((float) (((resultsSettings.getAveragePrecision() <= 0) ? nmPerPixel
              : resultsSettings.getAveragePrecision()) / nmPerPixel));
        } else if (resultsImage == ResultsImageType.DRAW_INTENSITY_PRECISION
            || resultsImage == ResultsImageType.DRAW_LOCALISATIONS_PRECISION) {
          image2.setCalculatedPrecision(true);
        }
        image = image2;
        break;

      default:
        image = new ImageJImagePeakResults(title, bounds, imageScale);
    }
    int flags = 0;

    switch (resultsImage) {
      case DRAW_INTENSITY:
      case DRAW_INTENSITY_PRECISION:
      case DRAW_INTENSITY_AVERAGE_PRECISION:
      case DRAW_FITTED_PSF:
        flags |= ImageJImagePeakResults.DISPLAY_SIGNAL;
        break;
      case DRAW_FRAME_NUMBER:
        flags |= ImageJImagePeakResults.DISPLAY_PEAK;
        break;
      case DRAW_FIT_ERROR:
        flags |= ImageJImagePeakResults.DISPLAY_ERROR;
        break;
      case DRAW_Z_POSITION:
        flags |= ImageJImagePeakResults.DISPLAY_Z_POSITION;
        break;
      case DRAW_ID:
        flags |= ImageJImagePeakResults.DISPLAY_ID;
        break;
      default:
        // Nothing to do for the other cases
        break;
    }

    switch (resultsSettings.getImageMode()) {
      case IMAGE_MAX:
        flags |= ImageJImagePeakResults.DISPLAY_MAX;
        break;
      case IMAGE_REPLACE:
        flags |= ImageJImagePeakResults.DISPLAY_REPLACE;
        break;
      default:
        // Nothing to do for the other cases
        break;
    }

    if (resultsSettings.getWeighted()) {
      flags |= ImageJImagePeakResults.DISPLAY_WEIGHTED;
    }
    if (resultsSettings.getEqualised()) {
      flags |= ImageJImagePeakResults.DISPLAY_EQUALIZED;
    }
    image.setDisplayFlags(flags);
    image.setRollingWindowSize(resultsSettings.getRollingWindowSize());
    image.setLutName(resultsSettings.getLutName());
    return image;
  }

  /**
   * Gets the image scale based on the image rendering settings.
   *
   * @param resultsSettings The results image settings
   * @param bounds Define the bounding rectangle of the result coordinates
   * @param nmPerPixel The results scale in nanometers per pixel (ignored unless image type uses
   *        fixed size pixels)
   * @return the scale
   */
  public static float getScale(ResultsImageSettingsOrBuilder resultsSettings, Rectangle bounds,
      double nmPerPixel) {
    float imageScale;
    switch (resultsSettings.getImageSizeMode()) {
      case IMAGE_SIZE:
        // Note: This may create an output image of the wrong size due to
        // floating-point rounding error.
        final int imageSize = resultsSettings.getImageSize();
        final int max = Math.max(bounds.width, bounds.height);
        if (imageSize > 0 && max > 0) {
          imageScale = (float) (imageSize / (double) max);
          // Check the output size and correct rounding errors.
          // The scale is likely to be too large due to the use of ceil.
          // No cases of the scale being too small are known.
          // This would occur if imageSize / max is exact as a double but has
          // more than 24 bits in the mantissa causing rounding down.
          // Ignore this possibility given that imageSize and max are effectively
          // limited to much less than 31-bits due to the limit imposed by
          // 2D array images. A single dimension is likely limited to 2^16.
          final int size = (int) Math.ceil(max * imageScale);
          if (size > imageSize) {
            // Use a single offset correction of half a pixel
            final float delta = Math.max(0.5f / max, Math.ulp(imageScale));
            imageScale -= delta;
          }
        } else {
          imageScale = 0;
        }
        break;
      case PIXEL_SIZE:
        imageScale = (float) (nmPerPixel / resultsSettings.getPixelSize());
        break;
      case SCALED:
      default:
        imageScale = (float) resultsSettings.getScale();
        break;
    }
    // Handle bad scaling
    return imageScale > 0 && imageScale < Double.POSITIVE_INFINITY ? imageScale : 1;
  }
}
