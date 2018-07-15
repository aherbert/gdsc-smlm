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
package gdsc.smlm.search;

import java.util.Arrays;

/**
 * Calculate the score of points within a search space.
 *
 * @param <T>
 *            the generic type
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
	 * @param <T>
	 *            the generic type
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
