package gdsc.smlm.utils;

import java.util.Arrays;

import org.apache.commons.math3.util.FastMath;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Simple class to calculate statistics of data
 */
public class Maths
{
	public static double min(double... data)
	{
		if (data == null || data.length == 0)
			return Double.NaN;
		return minDefault(Double.POSITIVE_INFINITY, data);
	}

	public static double minDefault(double min, double... data)
	{
		if (data == null || data.length == 0)
			return min;
		for (double d : data)
			if (min > d)
				min = d;
		return min;
	}

	public static double max(double... data)
	{
		if (data == null || data.length == 0)
			return Double.NaN;
		return maxDefault(Double.NEGATIVE_INFINITY, data);
	}

	public static double maxDefault(double max, double... data)
	{
		if (data == null || data.length == 0)
			return max;
		for (double d : data)
			if (max < d)
				max = d;
		return max;
	}

	/**
	 * Compute the min and max of the data
	 * 
	 * @param data
	 * @return [min, max]
	 */
	public static double[] limits(double... data)
	{
		if (data == null || data.length == 0)
			return noDoubleLimits();
		return limits(null, data);
	}

	private static double[] noDoubleLimits()
	{
		return new double[] { Double.NaN, Double.NaN };
	}

	/**
	 * Compute the min and max of the data
	 * 
	 * @param limits
	 *            The current [min, max]
	 * @param data
	 * @return [min, max]
	 */
	public static double[] limits(double[] limits, double... data)
	{
		if (data == null || data.length == 0)
			return (limits == null || limits.length < 2) ? noDoubleLimits() : limits;
		if (limits == null || limits.length < 2)
			limits = new double[] { Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY };
		double min = limits[0];
		double max = limits[1];
		for (double d : data)
		{
			if (min > d)
				min = d;
			if (max < d)
				max = d;
		}
		limits[0] = min;
		limits[1] = max;
		return limits;
	}

	public static float min(float... data)
	{
		if (data == null || data.length == 0)
			return Float.NaN;
		return minDefault(Float.POSITIVE_INFINITY, data);
	}

	public static float minDefault(float min, float... data)
	{
		if (data == null || data.length == 0)
			return min;
		for (float d : data)
			if (min > d)
				min = d;
		return min;
	}

	public static float max(float... data)
	{
		if (data == null || data.length == 0)
			return Float.NaN;
		return maxDefault(Float.NEGATIVE_INFINITY, data);
	}

	public static float maxDefault(float max, float... data)
	{
		if (data == null || data.length == 0)
			return max;
		for (float d : data)
			if (max < d)
				max = d;
		return max;
	}

	/**
	 * Compute the min and max of the data
	 * 
	 * @param data
	 * @return [min, max]
	 */
	public static float[] limits(float... data)
	{
		if (data == null || data.length == 0)
			return noFloatLimits();
		return limits(null, data);
	}

	private static float[] noFloatLimits()
	{
		return new float[] { Float.NaN, Float.NaN };
	}

	/**
	 * Compute the min and max of the data
	 * 
	 * @param limits
	 *            The current [min, max]
	 * @param data
	 * @return [min, max]
	 */
	public static float[] limits(float[] limits, float... data)
	{
		if (data == null || data.length == 0)
			return (limits == null || limits.length < 2) ? noFloatLimits() : limits;
		if (limits == null || limits.length < 2)
			limits = new float[] { Float.POSITIVE_INFINITY, Float.NEGATIVE_INFINITY };
		float min = limits[0];
		float max = limits[1];
		for (float d : data)
		{
			if (min > d)
				min = d;
			if (max < d)
				max = d;
		}
		limits[0] = min;
		limits[1] = max;
		return limits;
	}

	public static int min(int... data)
	{
		return minDefault(Integer.MAX_VALUE, data);
	}

	public static int minDefault(int min, int... data)
	{
		if (data == null || data.length == 0)
			return min;
		for (int d : data)
			if (min > d)
				min = d;
		return min;
	}

	public static int max(int... data)
	{
		return maxDefault(Integer.MIN_VALUE, data);
	}

	public static int maxDefault(int max, int... data)
	{
		if (data == null || data.length == 0)
			return max;
		for (int d : data)
			if (max < d)
				max = d;
		return max;
	}

	/**
	 * Compute the min and max of the data
	 * 
	 * @param data
	 * @return [min, max]
	 */
	public static int[] limits(int... data)
	{
		if (data == null || data.length == 0)
			return noIntegerLimits();
		return limits(null, data);
	}

	private static int[] noIntegerLimits()
	{
		return new int[] { 0, 0 };
	}

	/**
	 * Compute the min and max of the data
	 * 
	 * @param limits
	 *            The current [min, max]
	 * @param data
	 * @return [min, max]
	 */
	public static int[] limits(int[] limits, int... data)
	{
		if (data == null || data.length == 0)
			return (limits == null || limits.length < 2) ? noIntegerLimits() : limits;
		if (limits == null || limits.length < 2)
			limits = new int[] { Integer.MAX_VALUE, Integer.MIN_VALUE };
		int min = limits[0];
		int max = limits[1];
		for (int d : data)
		{
			if (min > d)
				min = d;
			if (max < d)
				max = d;
		}
		limits[0] = min;
		limits[1] = max;
		return limits;
	}

	public static short min(short... data)
	{
		return minDefault(Short.MAX_VALUE, data);
	}

	public static short minDefault(short min, short... data)
	{
		if (data == null || data.length == 0)
			return min;
		for (short d : data)
			if (min > d)
				min = d;
		return min;
	}

	public static short max(short... data)
	{
		return maxDefault(Short.MIN_VALUE, data);
	}

	public static short maxDefault(short max, short... data)
	{
		if (data == null || data.length == 0)
			return max;
		for (short d : data)
			if (max < d)
				max = d;
		return max;
	}

	/**
	 * Compute the min and max of the data
	 * 
	 * @param data
	 * @return [min, max]
	 */
	public static short[] limits(short... data)
	{
		if (data == null || data.length == 0)
			return noShortLimits();
		return limits(null, data);
	}

	private static short[] noShortLimits()
	{
		return new short[] { 0, 0 };
	}

	/**
	 * Compute the min and max of the data
	 * 
	 * @param limits
	 *            The current [min, max]
	 * @param data
	 * @return [min, max]
	 */
	public static short[] limits(short[] limits, short... data)
	{
		if (data == null || data.length == 0)
			return (limits == null || limits.length < 2) ? noShortLimits() : limits;
		if (limits == null || limits.length < 2)
			limits = new short[] { Short.MAX_VALUE, Short.MIN_VALUE };
		short min = limits[0];
		short max = limits[1];
		for (short d : data)
		{
			if (min > d)
				min = d;
			if (max < d)
				max = d;
		}
		limits[0] = min;
		limits[1] = max;
		return limits;
	}

	/**
	 * Calculate a cumulative histogram of the input values. The data is sorted and the first value in the returned
	 * values array will be the lowest value. NaN are ignored.
	 * 
	 * @param values
	 * @param normalise
	 *            Normalise so the total is 1
	 * @return Histogram values and cumulative total
	 */
	public static double[][] cumulativeHistogram(double[] values, boolean normalise)
	{
		if (values == null || values.length == 0)
			return new double[2][0];

		values = Arrays.copyOf(values, values.length);
		Arrays.sort(values);

		// Arrays.sort() put the NaN values higher than all others. If this is the first value then stop
		if (Double.isNaN(values[0]))
			return new double[2][0];

		double[] sum = new double[values.length];
		double lastValue = values[0];
		int position = 0, count = 0;
		for (int i = 0; i < values.length; i++)
		{
			// Arrays.sort() put the NaN values higher than all others so this should occur at the end
			if (Double.isNaN(values[i]))
				break;

			// When a new value is reached, store the cumulative total for the previous value
			if (lastValue != values[i])
			{
				values[position] = lastValue;
				sum[position] = count;
				lastValue = values[i];
				position++;
			}
			count++;
		}

		// Record the final value
		values[position] = lastValue;
		sum[position] = count;
		position++;

		// Truncate if necessary
		if (position < values.length)
		{
			values = Arrays.copyOf(values, position);
			sum = Arrays.copyOf(sum, position);
		}

		// Normalise
		if (normalise)
		{
			for (int i = 0; i < sum.length; i++)
			{
				sum[i] /= count;
			}
		}

		return new double[][] { values, sum };
	}

	/**
	 * @param sumOfSquaredResiduals
	 *            the sum of squared residuals from the nonlinear least-squares fit
	 * @param n
	 *            The number of data points
	 * @param p
	 *            The number of fitted parameters
	 * @return The Information Criterion
	 */
	public static double getInformationCriterion(double sumOfSquaredResiduals, int n, int p)
	{
		final double logLikelihood = 0.5 * (-n * (Math.log(2 * Math.PI) + 1 - Math.log(n) + Math
				.log(sumOfSquaredResiduals)));

		// Note: The true bias corrected AIC is derived from the 2nd, 3rd and 4th derivatives of the 
		// negative log-likelihood function. This is complex and so is not implemented.
		// See: 
		// http://www.math.sci.hiroshima-u.ac.jp/stat/TR/TR11/TR11-06.pdf
		// http://www.sciencedirect.com/science/article/pii/S0024379512000821#

		// This paper explains that the AIC or BIC are much better than the Adjusted coefficient of determination
		// for model selection:
		// http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2892436/

		//double aic = 2.0 * p - 2.0 * logLikelihood;

		// The Bias Corrected Akaike Information Criterion (AICc)
		// http://en.wikipedia.org/wiki/Akaike_information_criterion#AICc
		// Assumes a univariate linear model.
		//aic = aic + (2.0 * p * (p + 1)) / (n - p - 1);

		// Optimised 
		final double ic = 2.0 * (p - logLikelihood) + (2.0 * p * (p + 1)) / (n - p - 1);

		// Bayesian Information Criterion (BIC), which gives a higher penalty on the number of parameters
		// http://en.wikipedia.org/wiki/Bayesian_information_criterion
		//final double ic = p * Math.log(n) - 2.0 * logLikelihood;

		return ic;
	}
}