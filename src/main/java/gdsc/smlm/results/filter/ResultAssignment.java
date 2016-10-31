package gdsc.smlm.results.filter;

import gdsc.core.match.FractionalAssignment;

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
	public int compareTo(ResultAssignment that)
	{
		if (this.distance < that.distance)
			return -1;
		if (this.distance > that.distance)
			return 1;
		return 0;
	}
}