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
package gdsc.smlm.results.filter;

import gdsc.core.match.FractionalAssignment;

/**
 * Stores a base assignment between two identified points, the distance between them and the score for the match.
 * <p>
 * This class is used to store the data for a fractional assignment between a single predicted point and many actual
 * target points. It can be converted to a fractional assignment by specifying a predicted Id.
 */
public class ResultAssignment implements Comparable<ResultAssignment>
{
	/**
	 * The ID of the result that this assignment matches
	 */
	public final int targetId;
	/**
	 * The distance to the target
	 */
	public final double distance;
	/**
	 * The score for matching the target
	 */
	public final double score;

	/**
	 * Instantiates a new base assignment.
	 *
	 * @param targetId
	 *            the target id
	 * @param distance
	 *            the distance (zero is perfect match)
	 * @param score
	 *            The true positive score (must be 0-1)
	 */
	public ResultAssignment(int targetId, double distance, double score)
	{
		this.targetId = targetId;
		this.distance = distance;
		this.score = score;
	}

	/**
	 * Create a FractionalAssignment using the specified predicted id and a reference the result.
	 *
	 * @param predictedId
	 *            the predicted id
	 * @param peakResult
	 *            the peak result
	 * @return the fractional assignment
	 */
	public FractionalAssignment toFractionalAssignment(final int predictedId, final PreprocessedPeakResult peakResult)
	{
		return new PeakFractionalAssignment(targetId, predictedId, distance, score, peakResult);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	@Override
	public int compareTo(ResultAssignment that)
	{
		if (this.distance < that.distance)
			return -1;
		if (this.distance > that.distance)
			return 1;
		return 0;
	}
}
