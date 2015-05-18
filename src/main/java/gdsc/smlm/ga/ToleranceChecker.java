package gdsc.smlm.ga;

import org.apache.commons.math3.util.FastMath;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2015 Alex Herbert
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
 * Check if converged using a tolerance on the fitness and/or sequence change, and the number of iterations
 */
public class ToleranceChecker implements ConvergenceChecker
{
	final double relative, absolute;
	final boolean checkFitness, checkSequence;
	final int maxIterations;

	private int iterations = 0;

	/**
	 * Build an instance with specified thresholds. This only check convergence using the fitness.
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
	public ToleranceChecker(double relative, double absolute)
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
	 * @param checkFitness
	 *            Set to true to check the fitness
	 * @param checkSequence
	 *            Set to true to check the sequence
	 * @param maxIterations
	 *            Set above zero to check the iterations (number of time {@link #converged(Chromosome, Chromosome)} is
	 *            called)
	 * @throws IllegalArgumentException
	 *             if none of the convergence criteria are valid
	 */
	public ToleranceChecker(double relative, double absolute, boolean checkFitness, boolean checkSequence,
			int maxIterations)
	{
		if (maxIterations < 0)
			maxIterations = 0;
		boolean canConverge = maxIterations != 0;

		if (checkFitness || checkSequence)
			canConverge |= (relative > 0 || absolute > 0);

		if (!canConverge)
			throw new IllegalArgumentException("No valid convergence criteria");

		this.relative = relative;
		this.absolute = absolute;
		this.checkFitness = checkFitness;
		this.checkSequence = checkSequence;
		this.maxIterations = maxIterations;
	}

	/**
	 * Check if the sequence has converged
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
	private boolean converged(final double pi, final double ci)
	{
		final double difference = Math.abs(pi - ci);
		final double size = FastMath.max(Math.abs(pi), Math.abs(ci));
		if (difference > size * relative && difference > absolute)
		{
			return false;
		}
		return true;
	}

	public boolean converged(Chromosome previous, Chromosome current)
	{
		iterations++;
		if (maxIterations != 0 && iterations >= maxIterations)
			return true;
		if (checkFitness && converged(previous.getFitness(), current.getFitness()))
			return true;
		if (checkSequence && converged(previous.sequence(), current.sequence()))
			return true;
		return false;
	}

	/**
	 * @return the iterations
	 */
	public int getIterations()
	{
		return iterations;
	}
}
