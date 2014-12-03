package gdsc.smlm.function;

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
 * f(x) = l(x) - k * ln(l(x))<br/>
 * Where:<br/>
 * l(x) is the function (expected) value<br/>
 * k is the observed value
 * <p>
 * The negative log-likelihood (and gradient) can be evaluated over the entire set of observed values or for a chosen
 * observed value.
 */
public class PoissonLikelihoodWrapper extends LikelihoodWrapper
{
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
	 */
	public PoissonLikelihoodWrapper(NonLinearFunction f, double[] a, double[] k, int n)
	{
		super(f, a, k, n);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.LikelihoodWrapper#computeLikelihood()
	 */
	public double computeLikelihood()
	{
		// Compute the negative log-likelihood to be minimised
		double ll = 0;
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
		// f(x) = l(x) - k * ln(l(x))
		// 
		// Since (k * ln(l(x)))' = (k * ln(l(x))') * l'(x)
		//                       = (k / l(x)) * l'(x)

		// f'(x) = l'(x) - (k/l(x) * l'(x))
		// f'(x) = l'(x) * (1 - k/l(x))

		double ll = 0;
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
		return l - k * Math.log(l);
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
		return l - k * Math.log(l);
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