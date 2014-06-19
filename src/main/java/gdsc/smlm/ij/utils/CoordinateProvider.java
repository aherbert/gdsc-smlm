package gdsc.smlm.ij.utils;

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
 * Provide methods to retrieve the coordinates from a text line
 */
public interface CoordinateProvider
{
	/**
	 * Get the coordinates from the line
	 * @param line
	 * @return The coordinates (or null)
	 */
	double[] getCoordinates(String line);
}
