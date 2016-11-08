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
public interface ScoreFunction<T extends Comparable<T>>
{
	/**
	 * Find the optimum point. Return the best point from the input points with a score that can be compared to other
	 * results.
	 *
	 * @param points
	 *            the points
	 * @return the result for the optimum of the points
	 */
	SearchResult<T> findOptimum(double[][] points);
}
