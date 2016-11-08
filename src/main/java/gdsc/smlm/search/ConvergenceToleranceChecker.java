package gdsc.smlm.search;

import org.apache.commons.math3.util.FastMath;

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
 * 
 * This code is based on the ideas expressed in Numerical Recipes in C++, 
 * The Art of Scientific Computing, Second Edition, W.H. Press, 
 * S.A. Teukolsky, W.T. Vetterling, B.P. Flannery (Cambridge University Press, 
 * Cambridge, 2002).
 *---------------------------------------------------------------------------*/

/**
 * Check if converged using a tolerance on the score and/or position change, and the number of iterations
 */
public class ConvergenceToleranceChecker<T extends Comparable<T>> implements ConvergenceChecker<T>
{
	final double relative, absolute;
	final boolean checkScore, checkSequence;
	final int maxIterations;

	private int iterations = 0;

	/**
	 * Build an instance with specified thresholds. This only check convergence using the score.
	 *
	 * In order to perform only relative checks, the absolute tolerance
	 * must be set to a negative value. In order to perform only absolute
	 * checks, the relative tolerance must be set to a negative value.
	 *
	 * @param relativeThreshold
	 *            relative tolerance threshold
	 * @param absoluteThreshold
	 *            absolute tolerance threshold
	 * @throws IllegalArgumentException
	 *             if none of the convergence criteria are valid
	 */
	public ConvergenceToleranceChecker(double relative, double absolute)
	{
		this(relative, absolute, true, false, 0);
	}

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
	 * @param checkScore
	 *            Set to true to check the score
	 * @param checkSequence
	 *            Set to true to check the position
	 * @param maxIterations
	 *            Set above zero to check the iterations (number of time {@link #converged(Chromosome, Chromosome)} is
	 *            called)
	 * @throws IllegalArgumentException
	 *             if none of the convergence criteria are valid
	 */
	public ConvergenceToleranceChecker(double relative, double absolute, boolean checkScore, boolean checkSequence,
			int maxIterations)
	{
		if (maxIterations < 0)
			maxIterations = 0;
		boolean canConverge = maxIterations != 0;

		if (checkScore || checkSequence)
			canConverge |= (relative > 0 || absolute > 0);

		if (!canConverge)
			noConvergenceCriteria();

		this.relative = relative;
		this.absolute = absolute;
		this.checkScore = checkScore;
		this.checkSequence = checkSequence;
		this.maxIterations = maxIterations;
	}
	
	protected void noConvergenceCriteria()
	{
		throw new IllegalArgumentException("No valid convergence criteria");
	}

	/**
	 * Check if the position has converged
	 * 
	 * @param p
	 *            Previous
	 * @param c
	 *            Current
	 * @return True if converged
	 */
	private boolean converged(final double[] p, final double[] c)
	{
		for (int i = 0; i < p.length; ++i)
		{
			if (!converged(p[i], c[i]))
			{
				return false;
			}
		}
		return true;
	}

	/**
	 * Check if the value has converged
	 * 
	 * @param p
	 *            Previous
	 * @param c
	 *            Current
	 * @return True if converged
	 */
	private boolean converged(final double p, final double c)
	{
		final double difference = Math.abs(p - c);
		final double size = FastMath.max(Math.abs(p), Math.abs(c));
		if (difference > size * relative && difference > absolute)
		{
			return false;
		}
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.search.ConvergenceChecker#converged(gdsc.smlm.search.ScoreResult, gdsc.smlm.search.ScoreResult)
	 */
	public boolean converged(SearchResult<T> previous, SearchResult<T> current)
	{
		iterations++;
		if (maxIterations != 0 && iterations >= maxIterations)
			return true;
		if (checkScore && converged(previous.score, current.score))
			return true;
		if (checkSequence && converged(previous.point, current.point))
			return true;
		return false;
	}

	private boolean converged(T score, T score2)
	{
		return score.compareTo(score2) == 0;
	}

	/**
	 * @return the iterations
	 */
	public int getIterations()
	{
		return iterations;
	}
}
