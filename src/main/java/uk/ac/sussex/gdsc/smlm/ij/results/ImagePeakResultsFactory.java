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
   */
  public static ImageJImagePeakResults createPeakResultsImage(ResultsImageType resultsImage,
      boolean weighted, boolean equalised, String title, Rectangle bounds, double nmPerPixel,
      double imageScale, double precision, ResultsImageMode mode) {
    ImageJImagePeakResults image;
    switch (resultsImage) {
      case DRAW_FITTED_PSF:
      case DRAW_INTENSITY_PRECISION:
      case DRAW_LOCALISATIONS_PRECISION:
      case DRAW_INTENSITY_AVERAGE_PRECISION:
      case DRAW_LOCALISATIONS_AVERAGE_PRECISION:
        // Special case for full PSF image
        final PsfImagePeakResults image2 =
            new PsfImagePeakResults(title, bounds, (float) imageScale);
        if (resultsImage == ResultsImageType.DRAW_INTENSITY_AVERAGE_PRECISION
            || resultsImage == ResultsImageType.DRAW_LOCALISATIONS_AVERAGE_PRECISION) {
          // Fixed width display (in pixels)
          image2.setWidth((float) (((precision <= 0) ? nmPerPixel : precision) / nmPerPixel));
        } else if (resultsImage == ResultsImageType.DRAW_INTENSITY_PRECISION
            || resultsImage == ResultsImageType.DRAW_LOCALISATIONS_PRECISION) {
          image2.setCalculatedPrecision(true);
        }
        image = image2;
        break;

      default:
        image = new ImageJImagePeakResults(title, bounds, (float) imageScale);
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

    switch (mode) {
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

    if (weighted) {
      flags |= ImageJImagePeakResults.DISPLAY_WEIGHTED;
    }
    if (equalised) {
      flags |= ImageJImagePeakResults.DISPLAY_EQUALIZED;
    }
    image.setDisplayFlags(flags);
    image.setCalibration(nmPerPixel, 1);
    return image;
  }
}
