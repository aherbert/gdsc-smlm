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
 * Calculate the score of points within a search space
 */
public interface FullScoreFunction<T extends Comparable<T>> extends ScoreFunction<T>
{
	/**
	 * Return the score of the input points.
	 *
	 * @param points
	 *            the points
	 * @return the scores
	 */
	SearchResult<T>[] score(double[][] points);
}
