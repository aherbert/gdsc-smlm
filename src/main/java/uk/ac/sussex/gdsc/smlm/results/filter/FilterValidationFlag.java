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

package uk.ac.sussex.gdsc.smlm.results.filter;

/**
 * Support direct filtering of PreprocessedPeakResult objects.
 *
 * <p>The decision to support for filtering as both a DirectFilter and Filter concurrently is left
 * to the implementing class. It is not a requirement.
 */
public final class FilterValidationFlag {
  /**
   * Validation flag for the signal in photons.
   */
  public static final int PHOTONS = 0x000000001;

  /**
   * Validation flag for the SNR.
   */
  public static final int SNR = 0x000000002;

  /**
   * Validation flag for the noise.
   */
  public static final int NOISE = 0x000000004;

  /**
   * Validation flag for the location variance.
   */
  public static final int LOCATION_VARIANCE = 0x000000008;

  /**
   * Validation flag for the location variance using the local background.
   */
  public static final int LOCATION_VARIANCE2 = 0x000000010;

  /**
   * Validation flag for the average peak standard deviation in the X and Y dimension.
   */
  public static final int SD = 0x000000020;

  /**
   * Validation flag for the background.
   */
  public static final int BACKGROUND = 0x000000040;

  /**
   * Validation flag for the amplitude.
   */
  public static final int AMPLITUDE = 0x000000080;

  /**
   * Validation flag for the angle (for an elliptical Gaussian peak).
   */
  public static final int ANGLE = 0x000000100;

  /**
   * Validation flag for the x position.
   */
  public static final int X = 0x000000200;

  /**
   * Validation flag for the y position.
   */
  public static final int Y = 0x000000400;

  /**
   * Validation flag for the relative x position shift squared.
   */
  public static final int X_RELATIVE_SHIFT = 0x000000800;

  /**
   * Validation flag for the relative y position shift squared.
   */
  public static final int Y_RELATIVE_SHIFT = 0x000001000;

  /**
   * Validation flag for the x-dimension standard deviation.
   */
  public static final int X_SD = 0x000002000;

  /**
   * Validation flag for the y-dimension standard deviation.
   */
  public static final int Y_SD = 0x000004000;

  /**
   * Validation flag for the x-dimension width factor.
   */
  public static final int X_SD_FACTOR = 0x000008000;

  /**
   * Validation flag for the y-dimension width factor.
   */
  public static final int Y_SD_FACTOR = 0x000010000;

  /**
   * Validation flag for the location variance using the fitted x/y parameter Cramér-Rao lower
   * bound.
   */
  public static final int LOCATION_VARIANCE_CRLB = 0x000020000;

  /**
   * Validation flag for the z position.
   */
  public static final int Z = 0x000040000;

  /** No public construction. */
  private FilterValidationFlag() {}
}
