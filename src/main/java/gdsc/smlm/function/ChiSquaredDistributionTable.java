package gdsc.smlm.function;

import java.util.Arrays;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.special.Gamma;

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
 *---------------------------------------------------------------------------*/

/**
 * Computes probability values from the Chi-squared distribution
 * <p>
 * Note that Wilk's theorem states that the log-likelihood ratio asymptotically approaches a
 * Chi-squared distribution with degrees of freedom equal to the difference in dimensionality of the two models
 * (alternative and null).
 * 
 * @see https://en.wikipedia.org/wiki/Likelihood-ratio_test#Wilks.27_theorem
 */
public class ChiSquaredDistributionTable
{
	/**
	 * The p-value for computing the value using the cumulative probability (p=1-q)
	 */
	private final double p;

	final double[] chiSquared;

	/**
	 * Instantiates a new chi squared distribution table. Q is 1-p, with p the cumulative probability of the chi-squared
	 * distribution.
	 * <p>
	 * Setting q will cause {@link #getChiSquared(int)} to generate a value where the probability of a random
	 * chi-squared observation above that value is equal to q.
	 *
	 * @param q
	 *            the q-value for significance (lower will produce a higher allowed chi-squared value)
	 * @param df
	 *            the maximum degrees of freedom
	 */
	public ChiSquaredDistributionTable(double q, int df)
	{
		this.p = 1 - q;
		chiSquared = new double[df + 1];
		Arrays.fill(chiSquared, Double.NaN);
	}

	/**
	 * Get the q-value for significance.
	 * <p>
	 * Q refers to the probability of Chi squared being equal to or above a given value by chance.
	 *
	 * @return the q value
	 */
	public double getQValue()
	{
		return 1 - p;
	}

	/**
	 * Gets the chi squared value for the configured significance level and degrees of freedom.
	 * <p>
	 * This is the value of Chi squared where the chance of obtaining a value more extreme than this point is equal to
	 * the configured significance level (q-value).
	 *
	 * @param degreesOfFreedom
	 *            the degrees of freedom
	 * @return the chi squared
	 */
	public double getChiSquared(int degreesOfFreedom)
	{
		if (degreesOfFreedom >= chiSquared.length)
			throw new IllegalStateException("Maximum degrees of freedom = " + (chiSquared.length - 1));
		if (Double.isNaN(chiSquared[degreesOfFreedom]))
			chiSquared[degreesOfFreedom] = new ChiSquaredDistribution(null, degreesOfFreedom)
					.inverseCumulativeProbability(p);
		return chiSquared[degreesOfFreedom];
	}

	/**
	 * Checks if the value of chi-squared is significant at the configured significance level. Tests if Chi squared is
	 * below {@link #getChiSquared(int)}.
	 *
	 * @param chiSquared
	 *            the chi squared
	 * @param degreesOfFreedom
	 *            the degrees of freedom
	 * @return true, if is significant
	 */
	public boolean isSignificant(double chiSquared, int degreesOfFreedom)
	{
		double level = getChiSquared(degreesOfFreedom);
		return chiSquared < level;
	}

	/**
	 * Gets the chi squared value for the configured significance level and degrees of freedom.
	 *
	 * @param p
	 *            the p-value for significance (lower will produce a higher chi-squared value)
	 * @param degreesOfFreedom
	 *            the degrees of freedom
	 * @return the chi squared
	 */
	public static double getChiSquared(double p, int degreesOfFreedom)
	{
		return new ChiSquaredDistribution(null, degreesOfFreedom).inverseCumulativeProbability(1 - p);
	}

	/**
	 * Compute the q-value of the Chi-squared distribution.
	 * <p>
	 * This is the probability of obtaining a value more extreme than this point by chance. A chi-squared value is
	 * significant if q is higher than the p-value for significance (e.g. 0.05).
	 *
	 * @param chiSquared
	 *            the chi squared
	 * @param degreesOfFreedom
	 *            the degrees of freedom
	 * @return the q-value
	 */
	public static double computeQValue(double chiSquared, int degreesOfFreedom)
	{
		if (chiSquared <= 0)
		{
			return 1;
		}
		else
		{
			return Gamma.regularizedGammaQ(degreesOfFreedom / 2.0, chiSquared / 2.0);
		}
	}

	/**
	 * Compute the p-value of the Chi-squared distribution. This is the cumulative probability of the chi-squared
	 * distribution.
	 * <p>
	 * This is the probability of obtaining a value less extreme than this point by chance.
	 *
	 * @param chiSquared
	 *            the chi squared
	 * @param degreesOfFreedom
	 *            the degrees of freedom
	 * @return the p-value
	 */
	public static double computePValue(double chiSquared, int degreesOfFreedom)
	{
		if (chiSquared <= 0)
		{
			return 0;
		}
		else
		{
			return Gamma.regularizedGammaP(degreesOfFreedom / 2.0, chiSquared / 2.0);
		}
	}
}