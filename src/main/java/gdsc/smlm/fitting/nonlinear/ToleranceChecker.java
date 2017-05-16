package gdsc.smlm.fitting.nonlinear;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
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
 * Check if converged using a tolerance on the value, parameters, and/or the number of iterations
 */
public class ToleranceChecker
{
	final double relativeValue, absoluteValue;
	final double relativeParameters, absoluteParameters;
	final boolean checkValue, checkParameters, minimiseValue;
	final int maxIterations;

	private int iterations = 0;

	/**
	 * Build an instance with specified thresholds. This only checks convergence using the parameters.
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
		this(false, false, 0, 0, true, relative, absolute, 0);
	}

	/**
	 * Build an instance with specified thresholds.
	 * 
	 * In order to perform only relative checks, the absolute tolerance
	 * must be set to a negative value. In order to perform only absolute
	 * checks, the relative tolerance must be set to a negative value.
	 *
	 * @param checkValue
	 *            Set to true to check the value
	 * @param minimiseValue
	 *            Set to true to ensure the value is minimised at converge (otherwise it is maximised)
	 * @param relativeValue
	 *            relative tolerance threshold on the value
	 * @param absoluteValue
	 *            absolute tolerance threshold on the value
	 * @param checkParameters
	 *            Set to true to check the parameters
	 * @param relativeParameters
	 *            relative tolerance threshold on the parameters
	 * @param absoluteParameters
	 *            absolute tolerance threshold on the parameters
	 * @param maxIterations
	 *            Set above zero to limit the iterations. Set below zero to define a specific number of iterations.
	 * @throws IllegalArgumentException
	 *             if none of the convergence criteria are valid (i.e. convergence is not possible)
	 */
	public ToleranceChecker(boolean checkValue, boolean minimiseValue, double relativeValue, double absoluteValue,
			boolean checkParameters, double relativeParameters, double absoluteParameters, int maxIterations)
	{
		boolean canConverge = maxIterations != 0;

		if (!canConverge && checkValue)
			canConverge |= (relativeValue >= 0 || absoluteValue >= 0);
		if (!canConverge && checkParameters)
			canConverge |= (relativeParameters >= 0 || absoluteParameters >= 0);

		if (!canConverge)
			throw new IllegalArgumentException("No valid convergence criteria");

		this.checkValue = checkValue;
		this.minimiseValue = minimiseValue;
		this.relativeValue = relativeValue;
		this.absoluteValue = absoluteValue;
		this.checkParameters = checkParameters;
		this.relativeParameters = relativeParameters;
		this.absoluteParameters = absoluteParameters;
		this.maxIterations = maxIterations;
	}

	/**
	 * Check if the parameters have converged
	 * 
	 * @param p
	 *            Previous
	 * @param c
	 *            Current
	 * @param relative
	 *            relative tolerance threshold (set to negative to ignore)
	 * @param absolute
	 *            absolute tolerance threshold (set to negative to ignore)
	 * @return True if converged
	 */
	public boolean converged(final double[] p, final double[] c, double absolute, double relative)
	{
		for (int i = 0; i < p.length; ++i)
		{
			if (!converged(p[i], c[i], absolute, relative))
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
	 * @param relative
	 *            relative tolerance threshold (set to negative to ignore)
	 * @param absolute
	 *            absolute tolerance threshold (set to negative to ignore)
	 * @return True if converged
	 */
	public static boolean converged(final double p, final double c, double absolute, double relative)
	{
		final double difference = Math.abs(p - c);
		if (difference <= absolute)
			return true;
		final double size = max(Math.abs(p), Math.abs(c));
		return (difference <= size * relative);
	}

	private static double max(final double a, final double b)
	{
		// Ignore NaN
		return (a > b) ? a : b;
	}

	/** Flag to indicate convergence on the max iterations. */
	public static final int STATUS_MAX_ITERATIONS = 0x00000001;
	/** Flag to indicate convergence on the value. */
	public static final int STATUS_VALUE = 0x00000002;
	/** Flag to indicate convergence on the parameters. */
	public static final int STATUS_PARAMETERS = 0x00000004;
	/** Flag to indicate convergence on the target number of iterations. */
	public static final int STATUS_TARGET_ITERATIONS = 0x00000008;
	/** Flag to indicate all valid convergence flags. */
	public static final int STATUS_CONVERGED = STATUS_VALUE | STATUS_PARAMETERS | STATUS_TARGET_ITERATIONS;

	/**
	 * Check if converged
	 *
	 * @param previousValue
	 *            the previous value
	 * @param previousParameters
	 *            the previous parameters
	 * @param currentValue
	 *            the current value
	 * @param currentParameters
	 *            the current parameters
	 * @return The status flag
	 */
	public int converged(double previousValue, double[] previousParameters, double currentValue,
			double[] currentParameters)
	{
		iterations++;
		int status = 0;
		if (checkValue && correctDirection(previousValue, currentValue) &&
				converged(previousValue, currentValue, absoluteValue, relativeValue))
			status |= STATUS_VALUE;
		if (checkParameters && converged(previousParameters, currentParameters, absoluteParameters, relativeParameters))
			status |= STATUS_PARAMETERS;
		if (maxIterations != 0 && iterations >= Math.abs(maxIterations))
			status |= (maxIterations < 0) ? STATUS_TARGET_ITERATIONS : STATUS_MAX_ITERATIONS;
		return status;
	}

	private boolean correctDirection(double previousValue, double currentValue)
	{
		return (minimiseValue) ? currentValue <= previousValue : currentValue >= previousValue;
	}

	/**
	 * @return the iterations
	 */
	public int getIterations()
	{
		return iterations;
	}
}
