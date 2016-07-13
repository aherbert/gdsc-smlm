package gdsc.smlm.results.filter;

// TODO: Auto-generated Javadoc
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
 * Returns the settings of the filter for multiple thresholds: Signal, SNR, width, coordinate shift and precision.
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
	 * Gets the Signa-to-noise ratio.
	 *
	 * @return the snr
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
	 * Checks if precision uses local background.
	 *
	 * @return true, if precision uses local background
	 */
	boolean isPrecisionUsesLocalBackground();
}