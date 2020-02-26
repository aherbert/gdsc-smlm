/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.ij.settings;

import java.util.Map;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.ImagePSF;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.ImagePSFOrBuilder;

/**
 * Contain helper functions for the ImagePOSF class.
 */
public final class ImagePsfHelper {

  /** No public constructor. */
  private ImagePsfHelper() {}

  /**
   * Convert the ImagePSF to a string.
   *
   * @param imagePsf the image psf
   * @return the string
   */
  public static String toString(ImagePSFOrBuilder imagePsf) {
    return SettingsManager.toJson(imagePsf, SettingsManager.FLAG_JSON_WHITESPACE);
  }

  /**
   * Get the ImagePSF from a string.
   *
   * @param string the string
   * @return the image PSF
   */
  public static ImagePSF fromString(String string) {
    final ImagePSF.Builder builder = ImagePSF.newBuilder();
    if (SettingsManager.fromJson(string, builder)) {
      return builder.build();
    }
    return null;
  }

  /**
   * Creates the ImagePSF.
   *
   * @param centreImage the centre image
   * @param pixelSize the pixel size
   * @param pixelDepth the pixel depth
   * @param imageCount the image count
   * @param fwhm the fwhm
   * @return the image PSF
   */
  public static ImagePSF create(int centreImage, double pixelSize, double pixelDepth,
      int imageCount, double fwhm) {
    return create(centreImage, pixelSize, pixelDepth, imageCount, fwhm, null);
  }

  /**
   * Creates the ImagePSF.
   *
   * @param centreImage the centre image
   * @param pixelSize the pixel size
   * @param pixelDepth the pixel depth
   * @param imageCount the image count
   * @param fwhm the fwhm
   * @param notes the notes
   * @return the image PSF
   */
  public static ImagePSF create(int centreImage, double pixelSize, double pixelDepth,
      int imageCount, double fwhm, Map<String, String> notes) {
    final ImagePSF.Builder builder = ImagePSF.newBuilder();
    builder.setCentreImage(centreImage);
    builder.setPixelSize(pixelSize);
    builder.setPixelDepth(pixelDepth);
    builder.setImageCount(imageCount);
    builder.setFwhm(fwhm);
    if (notes != null) {
      builder.putAllNotes(notes);
    }
    return builder.build();
  }
}
