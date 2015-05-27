package gdsc.smlm.utils;

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
 *---------------------------------------------------------------------------*/

/**
 * Provide a score function that ramps smoothly between the configured limits
 */
public class RampedScore
{
	public final double lower, upper;
	private final double range;

	/**
	 * Create the score function with the specified limits
	 * 
	 * @param lower
	 * @param upper
	 */
	public RampedScore(double lower, double upper)
	{
		if (lower > upper)
		{
			final double tmp = lower;
			lower = upper;
			upper = tmp;
		}
		this.upper = upper;
		this.lower = lower;
		this.range = upper - lower;
	}

	/**
	 * Provide a score between 0 and 1 for the value. Return 1 if below the lower limit, 0 if above the upper limit,
	 * otherwise ramp smoothly from 1 to 0.
	 * 
	 * @param value
	 * @return the score
	 */
	public double score(double value)
	{
		double score = 0;
		if (value <= upper)
		{
			if (value <= lower)
			{
				score = 1;
			}
			else
			{
				// Interpolate from the minimum to the maximum match distance:
				// Cosine
				score = (0.5 * (1 + Math.cos(((value - lower) / range) * Math.PI)));
			}
		}
		return score;
	}

	/**
	 * Provide a score between 0 and 1 for the value. Return 1 if below the lower limit, 0 if above the upper limit,
	 * otherwise ramp smoothly from 1 to 0. Flatten the score to a new score that will have a maximum number of steps
	 * between 0 and 1.
	 * 
	 * @param value
	 * @param steps
	 * @return the score
	 */
	public double scoreAndFlatten(double value, int steps)
	{
		return flatten(score(value), steps);
	}

	/**
	 * Flatten the score to a new score that will have a maximum number of steps between 0 and 1.
	 * 
	 * @param score
	 * @param steps
	 * @return The new score
	 */
	public static double flatten(double score, int steps)
	{
		return (Math.round(score * steps)) / (double) steps;
	}

	/**
	 * Flatten the score to a new score that will have a maximum number of steps between 0 and 1.
	 * 
	 * @param score
	 * @param steps
	 * @return The new score
	 */
	public static float flatten(float score, int steps)
	{
		return (Math.round(score * steps)) / (float) steps;
	}
}
