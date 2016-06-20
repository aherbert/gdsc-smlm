package gdsc.smlm.function;

import java.util.Arrays;

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

/**
 * This is a wrapper for any function to compute the negative log-likelihood assuming a Poisson distribution:<br/>
 * f(x) = l(x) - k * ln(l(x)) + log(k!)<br/>
 * Where:<br/>
 * l(x) is the function (expected) value<br/>
 * k is the observed value
 * <p>
 * The negative log-likelihood (and gradient) can be evaluated over the entire set of observed values or for a chosen
 * observed value.
 */
public class PoissonLikelihoodWrapper extends LikelihoodWrapper
{
	private static double[] logFactorial;
	private final double[] logFactorialK;
	private final double sumLogFactorialK;

	/** All long-representable factorials */
	static final long[] FACTORIALS = new long[] { 1l, 1l, 2l, 6l, 24l, 120l, 720l, 5040l, 40320l, 362880l, 3628800l,
			39916800l, 479001600l, 6227020800l, 87178291200l, 1307674368000l, 20922789888000l, 355687428096000l,
			6402373705728000l, 121645100408832000l, 2432902008176640000l };

	static
	{
		logFactorial = new double[FACTORIALS.length];
		for (int k = 0; k < FACTORIALS.length; k++)
			logFactorial[k] = Math.log(FACTORIALS[k]);
	}

	private static double[] initialiseFactorial(double[] data)
	{
		int max = 0;
		for (double d : data)
		{
			final int i = (int) d;
			if (i != d)
				throw new IllegalArgumentException("Input observed values must be integers: " + d);
			if (max < i)
				max = i;
		}

		if (logFactorial.length <= max)
			populate(max);

		final double[] f = new double[data.length];
		for (int i = 0; i < data.length; i++)
		{
			f[i] = logFactorial[(int) data[i]];
		}
		return f;
	}

	private static synchronized void populate(int n)
	{
		if (logFactorial.length <= n)
		{
			int k = logFactorial.length - 1;
			double logSum = logFactorial[k];

			logFactorial = Arrays.copyOf(logFactorial, n + 1);
			while (k < n)
			{
				k++;
				logSum += Math.log(k);
				logFactorial[k] = logSum;
			}
		}
	}

	/**
	 * Initialise the function.
	 * <p>
	 * The input parameters must be the full parameters for the non-linear function. Only those parameters with gradient
	 * indices should be passed in to the functions to obtain the value (and gradient).
	 * 
	 * @param f
	 *            The function to be used to calculated the expected values
	 * @param a
	 *            The initial parameters for the function
	 * @param k
	 *            The observed values
	 * @param n
	 *            The number of observed values
	 * @throws IllegalArgumentException
	 *             if the input observed values are not integers
	 */
	public PoissonLikelihoodWrapper(NonLinearFunction f, double[] a, double[] k, int n)
	{
		super(f, a, k, n);
		// Initialise the factorial table
		logFactorialK = initialiseFactorial(k);
		double sum = 0;
		for (double d : logFactorialK)
			sum += d;
		sumLogFactorialK = sum;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.LikelihoodWrapper#computeLikelihood()
	 */
	public double computeLikelihood()
	{
		// Compute the negative log-likelihood to be minimised
		// f(x) = l(x) - k * ln(l(x)) + log(k!)
		double ll = sumLogFactorialK;
		for (int i = 0; i < n; i++)
		{
			final double l = f.eval(i);

			// Check for zero and return the worst likelihood score
			if (l <= 0)
			{
				// Since ln(0) -> -Infinity
				return Double.POSITIVE_INFINITY;
			}

			final double k = data[i];
			ll += l - k * Math.log(l);
		}
		return ll;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.LikelihoodWrapper#computeLikelihood(double[])
	 */
	public double computeLikelihood(double[] gradient)
	{
		// Compute the negative log-likelihood to be minimised
		// f(x) = l(x) - k * ln(l(x)) + log(k!)
		// 
		// Since (k * ln(l(x)))' = (k * ln(l(x))') * l'(x)
		//                       = (k / l(x)) * l'(x)

		// f'(x) = l'(x) - (k/l(x) * l'(x))
		// f'(x) = l'(x) * (1 - k/l(x))

		double ll = sumLogFactorialK;
		for (int j = 0; j < nVariables; j++)
			gradient[j] = 0;
		double[] dl_da = new double[nVariables];
		for (int i = 0; i < n; i++)
		{
			final double l = f.eval(i, dl_da);

			final double k = data[i];

			// Check for zero and return the worst likelihood score
			if (l <= 0)
			{
				// Since ln(0) -> -Infinity
				return Double.POSITIVE_INFINITY;
			}
			ll += l - k * Math.log(l);

			// Continue to work out the gradient since this does not involve logs.
			// Note: if l==0 then we get divide by zero and a NaN value
			final double factor = 1 - k / l;
			for (int j = 0; j < gradient.length; j++)
			{
				//gradient[j] += dl_da[j] - (dl_da[j] * k / l);
				//gradient[j] += dl_da[j] * (1 - k / l);
				gradient[j] += dl_da[j] * factor;
			}
		}
		return ll;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.LikelihoodWrapper#computeLikelihood(int)
	 */
	public double computeLikelihood(int i)
	{
		final double l = f.eval(i);

		// Check for zero and return the worst likelihood score
		if (l <= 0)
		{
			// Since ln(0) -> -Infinity
			return Double.POSITIVE_INFINITY;
		}

		final double k = data[i];
		return l - k * Math.log(l) + logFactorialK[i];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.LikelihoodWrapper#computeLikelihood(double[], int)
	 */
	public double computeLikelihood(double[] gradient, int i)
	{
		for (int j = 0; j < nVariables; j++)
			gradient[j] = 0;
		double[] dl_da = new double[nVariables];
		final double l = f.eval(i, dl_da);

		// Check for zero and return the worst likelihood score
		if (l <= 0)
		{
			// Since ln(0) -> -Infinity
			return Double.POSITIVE_INFINITY;
		}

		final double k = data[i];
		final double factor = 1 - k / l;
		for (int j = 0; j < gradient.length; j++)
		{
			//gradient[j] = dl_da[j] - (dl_da[j] * k / l);
			//gradient[j] = dl_da[j] * (1 - k / l);
			gradient[j] = dl_da[j] * factor;
		}
		return l - k * Math.log(l) + logFactorialK[i];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.LikelihoodWrapper#canComputeGradient()
	 */
	@Override
	public boolean canComputeGradient()
	{
		return true;
	}
}