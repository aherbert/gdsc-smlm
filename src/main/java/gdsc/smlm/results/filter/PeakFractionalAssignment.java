package gdsc.smlm.results.filter;

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

import gdsc.core.match.ImmutableFractionalAssignment;

/**
 * Extends the fractional assignment to add a reference to the peak result
 */
public class PeakFractionalAssignment extends ImmutableFractionalAssignment
{
	public final PreprocessedPeakResult peakResult;

	/**
	 * Instantiates a new custom fractional assignment.
	 *
	 * @param targetId
	 *            the target id
	 * @param predictedId
	 *            the predicted id
	 * @param distance
	 *            the distance
	 * @param score
	 *            the score
	 * @param peakResult
	 *            the peak result
	 */
	public PeakFractionalAssignment(int targetId, int predictedId, double distance, double score,
			PreprocessedPeakResult peakResult)
	{
		super(targetId, predictedId, distance, score);
		this.peakResult = peakResult;
	}
}
