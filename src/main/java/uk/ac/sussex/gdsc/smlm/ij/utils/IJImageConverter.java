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

package uk.ac.sussex.gdsc.smlm.ij.utils;

import uk.ac.sussex.gdsc.core.annotation.Nullable;
import uk.ac.sussex.gdsc.smlm.utils.ImageConverter;

import ij.process.ImageProcessor;

import java.awt.Rectangle;

/**
 * Contains methods for converting an image to float data. Simple wrapper around
 * uk.ac.sussex.gdsc.smlm.utils.ImageConverter;
 */
public final class IJImageConverter {
  /**
   * The instance that does the conversion. This can be manipulated for different RGB colour
   * conversions.
   */
  public static final ImageConverter IMAGE_CONVERTER = new ImageConverter();

  /**
   * No public construction.
   */
  private IJImageConverter() {}

  /**
   * Get the data from the image processor as a float array (include cropping to the ROI).
   *
   * @param ip the image
   * @return The float array data
   */
  public static float[] getData(ImageProcessor ip) {
    return getData(ip, null);
  }

  /**
   * Get the data from the image processor as a float array (include cropping to the ROI). Data is
   * duplicated if the InputImage is a FloatProcessor.
   *
   * <p>Allows reuse of an existing buffer if provided. This will not be truncated if it is larger
   * than the ImageProcessor ROI bounds. If smaller then a new buffer will be created.
   *
   * @param ip the image
   * @param buffer the buffer
   * @return The float array data
   */
  public static float[] getData(ImageProcessor ip, float[] buffer) {
    return getData(ip.getPixels(), ip.getWidth(), ip.getHeight(), ip.getRoi(), buffer);
  }

  /**
   * Get the data from the image as a float array (include cropping to the ROI). Data is duplicated
   * if the input is already a float array.
   *
   * <p>Allows reuse of an existing buffer if provided. This will not be truncated if it is larger
   * than the ImageProcessor ROI bounds. If smaller then a new buffer will be created.
   *
   * <p>If the object pixels array is incorrect size (it should be width*height) then null will be
   * returned.
   *
   * @param pixels the pixels object
   * @param width the width
   * @param height the height
   * @param bounds the bounds
   * @param buffer the buffer
   * @return The float array data
   */
  public static float[] getData(final Object pixels, final int width, final int height,
      final Rectangle bounds, float[] buffer) {
    return IMAGE_CONVERTER.getData(pixels, width, height, bounds, buffer);
  }

  /**
   * Get the data from the image pixels as a float array. Data is not duplicated if the input is
   * already a float array unless a buffer is provided.
   *
   * <p>Allows reuse of an existing buffer if provided. This will not be truncated if it is larger
   * than the pixels array. If smaller then a new buffer will be created.
   *
   * @param pixels the pixels
   * @param buffer the buffer
   * @return The float array data
   */
  public static float[] getData(final Object pixels, float[] buffer) {
    return IMAGE_CONVERTER.getData(pixels, buffer);
  }

  /**
   * Get the data from the image processor as a double array (include cropping to the ROI). Data is
   * duplicated if the InputImage is a FloatProcessor.
   *
   * <p>Allows reuse of an existing buffer if provided. This will not be truncated if it is larger
   * than the ImageProcessor ROI bounds. If smaller then a new buffer will be created.
   *
   * @param ip the image
   * @param buffer the buffer
   * @return The double array data
   */
  public static @Nullable double[] getDoubleData(ImageProcessor ip, double[] buffer) {
    return getDoubleData(ip.getPixels(), ip.getWidth(), ip.getHeight(), ip.getRoi(), buffer);
  }

  /**
   * Get the data from the image as a double array (include cropping to the ROI). Data is duplicated
   * if the input is already a double array.
   *
   * <p>Allows reuse of an existing buffer if provided. This will not be truncated if it is larger
   * than the ImageProcessor ROI bounds. If smaller then a new buffer will be created.
   *
   * <p>If the object pixels array is incorrect size (it should be width*height) then null will be
   * returned.
   *
   * @param pixels the pixels object
   * @param width the width
   * @param height the height
   * @param bounds the bounds
   * @param buffer the buffer
   * @return The double array data
   */
  public static double[] getDoubleData(final Object pixels, final int width, final int height,
      final Rectangle bounds, double[] buffer) {
    return IMAGE_CONVERTER.getDoubleData(pixels, width, height, bounds, buffer);
  }
}
