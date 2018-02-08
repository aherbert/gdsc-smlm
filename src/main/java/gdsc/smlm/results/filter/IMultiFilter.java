package gdsc.smlm.results.filter;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Returns the settings of the filter for multiple thresholds: Signal, SNR, width, coordinate shift, precision and
 * z-depth.
 */
public interface IMultiFilter
{
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
	 * Gets the min allowed z depth. Note z==0 is the focal plane. If both the min and max z depth are zero then it is
	 * assumed that depth filtering is disabled.
	 *
	 * @return the min z depth
	 */
	double getMinZ();

	/**
	 * Gets the max allowed z depth. Note z==0 is the focal plane. If both the min and max z depth are zero then it is
	 * assumed that depth filtering is disabled.
	 *
	 * @return the max z depth
	 */
	double getMaxZ();
}