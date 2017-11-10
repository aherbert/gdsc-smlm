/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2014 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

package org.apache.commons.math3.optim;

/**
 * Check if the position has converged
 */
public class PositionChecker extends SimplePointChecker<PointValuePair>
		implements OptimizationData
{
	/**
	 * Build an instance with specified thresholds.
	 *
	 * In order to perform only relative checks, the absolute tolerance
	 * must be set to a negative value. In order to perform only absolute
	 * checks, the relative tolerance must be set to a negative value.
	 *
	 * @param relativeThreshold
	 *            relative tolerance threshold
	 * @param absoluteThreshold
	 *            absolute tolerance threshold
	 */
	public PositionChecker(double relative, double absolute)
	{
		super(relative, absolute);
	}
}
