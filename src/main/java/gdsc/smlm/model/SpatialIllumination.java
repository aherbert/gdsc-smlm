package gdsc.smlm.model;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Contains methods for specifying the illumination photons of a spatial position.
 */
public interface SpatialIllumination
{
	/**
	 * Get the number of photons for the position
	 * 
	 * @param xyz
	 * @return The photons
	 */
	double getPhotons(double[] xyz);

	/**
	 * Get the number of photons for the position at the specified time.
	 * <p>
	 * The return value is an array containing the number of photons that occurred before the time frame and then the
	 * number of photons during the time frame. This allows simulation of a pulsed illumination source where the pulse
	 * is modelled as a zero time event.
	 * 
	 * @param xyz
	 * @return The photons [before,during]
	 */
	double[] getPulsedPhotons(double[] xyz, int t);
	
	
	/**
	 * @return An estimate of the average photons per time frame
	 */
	double getAveragePhotons();
}
