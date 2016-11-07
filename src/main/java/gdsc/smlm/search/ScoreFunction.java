package gdsc.smlm.search;

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
 * Calculate the optimum point within a search space
 */
public interface ScoreFunction
{
	/**
	 * Find the optimum point.
	 *
	 * @param points
	 *            the points
	 * @return the optimum
	 */
	double[] findOptimum(double[][] points);
}
