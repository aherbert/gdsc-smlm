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
 * Support direct filtering of PreprocessedPeakResult objects.
 * <p>
 * The decision to support for filtering as both a DirectFilter and Filter concurrently is left to the implementing
 * class. It is not a requirement.
 */
public interface IDirectFilter
{
	/**
	 * Validation flag for the signal in photons
	 */
	final static int V_PHOTONS = 1;

	/**
	 * Validation flag for the SNR
	 */
	final static int V_SNR = 2;

	/**
	 * Validation flag for the noise
	 */
	final static int V_NOISE = 4;

	/**
	 * Validation flag for the location variance
	 */
	final static int V_LOCATION_VARIANCE = 8;

	/**
	 * Validation flag for the location variance using the local background
	 */
	final static int V_LOCATION_VARIANCE2 = 16;

	/**
	 * Validation flag for the average peak standard deviation in the X and Y dimension
	 */
	final static int V_SD = 32;

	/**
	 * Validation flag for the background
	 */
	final static int V_BACKGROUND = 64;

	/**
	 * Validation flag for the amplitude
	 */
	final static int V_AMPLITUDE = 128;

	/**
	 * Validation flag for the angle (for an elliptical Gaussian peak)
	 */
	final static int V_ANGLE = 256;

	/**
	 * Validation flag for the x position
	 */
	final static int V_X = 512;

	/**
	 * Validation flag for the y position
	 */
	final static int V_Y = 1024;

	/**
	 * Validation flag for the relative x position shift squared
	 */
	final static int V_X_RELATIVE_SHIFT = 2048;

	/**
	 * Validation flag for the relative y position shift squared
	 */
	final static int V_Y_RELATIVE_SHIFT = 4096;

	/**
	 * Validation flag for the x-dimension standard deviation
	 */
	final static int V_X_SD = 8192;

	/**
	 * Validation flag for the y-dimension standard deviation
	 */
	final static int V_Y_SD = 16384;

	/**
	 * Validation flag for the x-dimension width factor
	 */
	final static int V_X_SD_FACTOR = 32768;

	/**
	 * Validation flag for the y-dimension width factor
	 */
	final static int V_Y_SD_FACTOR = 65536;

	/**
	 * Disable filtering using the width of the result
	 */
	final static int NO_WIDTH = 1;

	/**
	 * Called before the accept method is called for PreprocessedPeakResult
	 * <p>
	 * This should be called once to initialise the filter before processing a batch of results.
	 * 
	 * @see #validate(PreprocessedPeakResult)
	 */
	void setup();

	/**
	 * Called before the accept method is called for PreprocessedPeakResult. the flags can control the type of filtering
	 * requested. Filters are asked to respect the flags defined in this class.
	 * <p>
	 * This should be called once to initialise the filter before processing a batch of results.
	 * 
	 * @param flags
	 *            Flags used to control the filter
	 * @see #validate(PreprocessedPeakResult)
	 */
	void setup(final int flags);

	/**
	 * Filter the peak result.
	 * <p>
	 * Calls {@link #validate(PreprocessedPeakResult)} and stores the result. This can be obtained using
	 * {@link #getResult()}.
	 * 
	 * @param peak
	 *            The peak result
	 * @return true if the peak should be accepted
	 */
	boolean accept(final PreprocessedPeakResult peak);

	/**
	 * Filter the peak result.
	 * 
	 * @param peak
	 *            The peak result
	 * @return zero if the peak should be accepted, otherwise set to flags indicating the field that failed validation.
	 */
	int validate(final PreprocessedPeakResult peak);

	/**
	 * Return the type of filter. This should be a DirectFilter.
	 * 
	 * @return Should return DirectFilter
	 */
	FilterType getFilterType();

	/**
	 * Return the result flag generated during the last call to {@link #accept(PreprocessedPeakResult)}.
	 * 
	 * @return the validation result from the last call to {@link #accept(PreprocessedPeakResult)}
	 */
	int getResult();
	
	/**
	 * Copy this filter.
	 *
	 * @return the copy
	 */
	IDirectFilter copy();
}