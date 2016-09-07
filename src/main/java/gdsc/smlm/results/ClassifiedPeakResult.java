package gdsc.smlm.results;

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
 * Specifies a fitted peak result that has been scored for classification analysis 
 */
public interface ClassifiedPeakResult
{
	/**
	 * Get the signal strength (e.g. the volume under the Gaussian peak, amplitude * 2 * pi * sx * sy)
	 * 
	 * @return The signal
	 */
	public float getSignal();
	
	/**
	 * Gets the noise.
	 *
	 * @return the noise
	 */
	public float getNoise();

	/**
	 * @return The average peak standard deviation in the X and Y dimension
	 */
	public float getSD();

	/**
	 * @return The background
	 */
	public float getBackground();

	/**
	 * Get the amplitude
	 * 
	 * @return The amplitude
	 */
	public float getAmplitude();

	/**
	 * @return The angle
	 */
	public float getAngle();

	/**
	 * @return The x position
	 */
	public float getXPosition();

	/**
	 * @return The y position
	 */
	public float getYPosition();

	/**
	 * @return The x shift
	 */
	public float getXShift();

	/**
	 * @return The y shift
	 */
	public float getYShift();

	/**
	 * @return The x-dimension standard deviation
	 */
	public float getXSD();

	/**
	 * @return The y-dimension standard deviation
	 */
	public float getYSD();

	/**
	 * @return The results identifier
	 */
	public int getId();

	/**
	 * Return the true positive score for use in classification analysis
	 * 
	 * @return The true positive score
	 */
	public double getTruePositiveScore();

	/**
	 * Return the false positive score for use in classification analysis
	 * 
	 * @return The false positive score
	 */
	public double getFalsePositiveScore();

	/**
	 * Return the true negative score for use in classification analysis
	 * 
	 * @return The true negative score
	 */
	public double getTrueNegativeScore();

	/**
	 * Return the false negative score for use in classification analysis
	 * 
	 * @return The false negative score
	 */
	public double getFalseNegativeScore();
}
