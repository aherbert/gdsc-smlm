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
 * Specifies a peak fitting result for use in filtering. Any result implementing this interface can be directly filtered
 * without the requiring the filter to be initialised with calibration data.
 */
public interface PreprocessedPeakResult
{
	/**
	 * Get the frame containing the result
	 * 
	 * @return The frame
	 */
	int getFrame();

	/**
	 * Get the signal strength (i.e. the volume under the Gaussian peak, amplitude * 2 * pi * sx * sy)
	 * 
	 * @return The signal (non calibrated)
	 */
	float getSignal();

	/**
	 * Get the signal strength (i.e. the volume under the Gaussian peak, amplitude * 2 * pi * sx * sy) calibrated to
	 * photons
	 * 
	 * @return The signal in photons
	 */
	float getPhotons();

	/**
	 * Get the signal-to-noise ratio (SNR)
	 * 
	 * @return The SNR
	 */
	float getSNR();

	/**
	 * Get the background noise
	 * 
	 * @return The noise
	 */
	float getNoise();

	/**
	 * The localisation variance
	 * 
	 * @return The location variance in nm
	 */
	double getLocationVariance();

	/**
	 * @return The average peak standard deviation in the X and Y dimension
	 */
	float getSD();

	/**
	 * @return The background
	 */
	float getBackground();

	/**
	 * Get the amplitude. Amplitude = Signal / (2*pi*sd0*sd1).
	 * 
	 * @return The amplitude
	 */
	float getAmplitude();

	/**
	 * @return The angle (for an elliptical Gaussian peak)
	 */
	float getAngle();

	/**
	 * @return The x position
	 */
	float getX();

	/**
	 * @return The y position
	 */
	float getY();

	/**
	 * Return the squared shift between the initial and the final x-position, i.e. how much the position moved during
	 * fitting, relative to the initial peak SD
	 * 
	 * @return The relative x position shift squared
	 */
	float getXRelativeShift2();

	/**
	 * Return the squared shift between the initial and the final y-position, i.e. how much the position moved during
	 * fitting, relative to the initial peak SD.
	 * 
	 * @return The relative y position shift squared
	 */
	float getYRelativeShift2();

	/**
	 * @return The x-dimension standard deviation
	 */
	float getXSD();

	/**
	 * @return The y-dimension standard deviation
	 */
	float getYSD();

	/**
	 * Return the ratio between the initial and the final x-dimension standard deviation, i.e. how much the width
	 * changed during fitting
	 * 
	 * @return The x-dimension width factor
	 */
	float getXSDFactor();

	/**
	 * Return the ratio between the initial and the final y-dimension standard deviation, i.e. how much the width
	 * changed during fitting
	 * 
	 * @return The y-dimension width factor
	 */
	float getYSDFactor();
}
