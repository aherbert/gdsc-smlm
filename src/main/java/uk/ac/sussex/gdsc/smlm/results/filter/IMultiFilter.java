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

package uk.ac.sussex.gdsc.smlm.results.filter;

/**
 * Returns the settings of the filter for multiple thresholds: Signal, SNR, width, coordinate shift,
 * precision and z-depth.
 */
public interface IMultiFilter {
  /**
   * Gets the signal.
   *
   * @return the signal
   */
  double getSignal();

  /**
   * Gets the Signal-to-noise ratio.
   *
   * @return the Signal-to-noise ratio (SNR)
   */
  double getSNR();

  /**
   * Gets the min width.
   *
   * @return the min width
   */
  double getMinWidth();

  /**
   * Gets the max width.
   *
   * @return the max width
   */
  double getMaxWidth();

  /**
   * Gets the coordinate shift.
   *
   * @return the shift
   */
  double getShift();

  /**
   * Gets the Euclidian shift.
   *
   * @return the Euclidian shift
   */
  double getEShift();

  /**
   * Gets the precision.
   *
   * @return the precision
   */
  double getPrecision();

  /**
   * Gets the precision type.
   *
   * @return the precision type
   */
  PrecisionType getPrecisionType();

  /**
   * Gets the min allowed z depth. Note z==0 is the focal plane. If both the min and max z depth are
   * zero then it is assumed that depth filtering is disabled.
   *
   * @return the min z depth
   */
  double getMinZ();

  /**
   * Gets the max allowed z depth. Note z==0 is the focal plane. If both the min and max z depth are
   * zero then it is assumed that depth filtering is disabled.
   *
   * @return the max z depth
   */
  double getMaxZ();
}
