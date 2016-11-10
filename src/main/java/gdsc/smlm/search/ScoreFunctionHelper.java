package gdsc.smlm.search;

import java.util.Arrays;

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
public class ScoreFunctionHelper<T extends Comparable<T>>
{
	/**
	 * Cut the list of scores down to the given size by selecting only the best results. The input list may not be
	 * sorted. The results should contain the best result at position 0 in the output array.
	 * <p>
	 * Helper implementation of the FullScoreFunction.cut(...) method. Uses a full sort then truncation to the given
	 * size.
	 * 
	 * @param scores
	 *            The scores
	 * @param size
	 *            The size
	 * @return The reduced list
	 */
	public static <T extends Comparable<T>> SearchResult<T>[] cut(SearchResult<T>[] scores, int size)
	{
		if (scores == null || scores.length == 1)
			return scores;
		Arrays.sort(scores);
		return (size < scores.length) ? Arrays.copyOf(scores, size) : scores;
	}
}
