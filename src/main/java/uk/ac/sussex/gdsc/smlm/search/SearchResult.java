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
package uk.ac.sussex.gdsc.smlm.search;

/**
 * Store the result of scoring a point within a search space. Allows the scores to be compared.
 *
 * @param <T>
 *            the generic type
 */
public class SearchResult<T extends Comparable<T>> implements Comparable<SearchResult<T>>
{
	private final double[] point;
	private final T score;

	/**
	 * Instantiates a new search result.
	 *
	 * @param point
	 *            the point
	 * @param score
	 *            the score
	 */
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
	@Override
	public int compareTo(SearchResult<T> o)
	{
		if (o == null)
			return -1;
		return getScore().compareTo(o.getScore());
	}

	/**
	 * Gets the point.
	 *
	 * @return the point
	 */
	public double[] getPoint()
	{
		return point;
	}

	/**
	 * Gets the score.
	 *
	 * @return the score
	 */
	public T getScore()
	{
		return score;
	}
}
