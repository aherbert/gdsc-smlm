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
 * Store the result of scoring a point within a search space. Allows the scores to be compared.
 */
public class SearchResult<T extends Comparable<T>> implements Comparable<SearchResult<T>>
{
	public final double[] point;
	public final T score;

	public SearchResult(double[] point, T score)
	{
		if (point == null)
			throw new IllegalArgumentException("Point is null");
		if (score == null)
			throw new IllegalArgumentException("Score is null");
		this.point = point;
		this.score = score;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo(SearchResult<T> o)
	{
		return score.compareTo(o.score);
	}
}
