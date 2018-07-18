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
package uk.ac.sussex.gdsc.smlm.ga;

import org.apache.commons.math3.util.FastMath;

/**
 * Check if converged using a tolerance on the fitness and/or sequence change, and the number of iterations.
 *
 * @param <T>
 *            the generic type
 */
public abstract class ToleranceChecker<T extends Comparable<T>> implements ConvergenceChecker<T>
{
	/** The relative tolerance threshold. */
	final double relative;
	/** The absolute tolerance threshold. */
	final double  absolute;
	/** Set to true to check the fitness */
	final boolean checkFitness;
	/** Set to true to check the sequence */
	final boolean checkSequence;
	/** The maximum iterations. */
	final int maxIterations;

	/** The iterations. */
	private int iterations = 0;

	/**
	 * Build an instance with specified thresholds. This only check convergence using the sequence.
	 *
	 * In order to perform only relative checks, the absolute tolerance
	 * must be set to a negative value. In order to perform only absolute
	 * checks, the relative tolerance must be set to a negative value.
	 *
	 * @param relative
	 *            relative tolerance threshold
	 * @param absolute
	 *            absolute tolerance threshold
	 * @throws IllegalArgumentException
	 *             if none of the convergence criteria are valid
	 */
	public ToleranceChecker(double relative, double absolute)
	{
		this(relative, absolute, false, true, 0);
	}

	/**
	 * Build an instance with specified thresholds.
	 *
	 * In order to perform only relative checks, the absolute tolerance
	 * must be set to a negative value. In order to perform only absolute
	 * checks, the relative tolerance must be set to a negative value.
	 *
	 * @param relative
	 *            relative tolerance threshold
	 * @param absolute
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
			canConverge |= (relative >= 0 || absolute >= 0);

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
			if (!converged(p[i], c[i]))
				return false;
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
	public boolean converged(final double p, final double c)
	{
		final double difference = Math.abs(p - c);
		final double size = FastMath.max(Math.abs(p), Math.abs(c));
		if (difference > size * relative && difference > absolute)
			return false;
		return true;
	}

	@Override
	public boolean converged(Chromosome<T> previous, Chromosome<T> current)
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
	 * Check if the fitness has converged.
	 *
	 * @param previous
	 *            the previous
	 * @param current
	 *            the current
	 * @return true, if successful
	 */
	protected abstract boolean converged(T previous, T current);

	/**
	 * @return the iterations
	 */
	public int getIterations()
	{
		return iterations;
	}
}
