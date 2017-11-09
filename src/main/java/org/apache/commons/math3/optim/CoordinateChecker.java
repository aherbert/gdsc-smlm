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
 * 
 * This code is based on the ideas expressed in Numerical Recipes in C++, 
 * The Art of Scientific Computing, Second Edition, W.H. Press, 
 * S.A. Teukolsky, W.T. Vetterling, B.P. Flannery (Cambridge University Press, 
 * Cambridge, 2002).
 *---------------------------------------------------------------------------*/

package org.apache.commons.math3.optim;

import org.apache.commons.math3.util.FastMath;

/**
 * Check if the coordinates have converged
 */
public abstract class CoordinateChecker implements OptimizationData, ConvergenceChecker<PointValuePair>
{
	final double relative, absolute;
	final int fixedIterations;

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
	public CoordinateChecker(double relative, double absolute)
	{
		this(relative, absolute, Integer.MAX_VALUE);
	}

	/**
	 * Build an instance with specified thresholds.
	 * 
	 * In order to perform only relative checks, the absolute tolerance
	 * must be set to a negative value. In order to perform only absolute
	 * checks, the relative tolerance must be set to a negative value.
	 *
	 * @param relative
	 *            the relative
	 * @param absolute
	 *            the absolute
	 * @param fixedIterations
	 *            the fixed number of iterations to signal convergence
	 */
	public CoordinateChecker(double relative, double absolute, int fixedIterations)
	{
		this.relative = relative;
		this.absolute = absolute;
		this.fixedIterations = fixedIterations;
	}

	/**
	 * Check if the coordinates have converged
	 * 
	 * @param p
	 *            Previous
	 * @param c
	 *            Current
	 * @return True if converged
	 */
	public boolean converged(final double[] p, final double[] c)
	{
		for (int i = 0; i < p.length; ++i)
		{
			final double pi = p[i];
			final double ci = c[i];
			final double difference = Math.abs(pi - ci);
			final double size = FastMath.max(Math.abs(pi), Math.abs(ci));
			if (difference > size * relative && difference > absolute)
			{
				return false;
			}
		}
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.apache.commons.math3.optim.ConvergenceChecker#converged(int, java.lang.Object, java.lang.Object)
	 */
	public boolean converged(int iteration, PointValuePair previous, PointValuePair current)
	{
		if (iteration >= fixedIterations)
			return true;
		return converged(previous.getPointRef(), current.getPointRef());
	}
}
