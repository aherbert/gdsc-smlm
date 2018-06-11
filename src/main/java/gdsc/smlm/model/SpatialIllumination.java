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
package gdsc.smlm.model;

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
