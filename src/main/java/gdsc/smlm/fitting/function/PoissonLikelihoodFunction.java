package gdsc.smlm.fitting.function;

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

import java.util.Arrays;

import org.apache.commons.math3.analysis.DifferentiableMultivariateFunction;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.MultivariateVectorFunction;

/**
 * This is a wrapper for any function to compute the negative log-likelihood assuming a Poisson distribution:<br/>
 * f(x) = l(x) - k * ln(l(x))<br/>
 * Where:<br/>
 * l(x) is the function (expected) value<br/>
 * k is the observed value
 * <p>
 * The negative log-likelihood (and gradient) can be evaluated over the entire set of observed values or for a chosen
 * observed value.
 * <p>
 * This implements the Apache Commons Math API so that optimisers can be used for Maximum Likelihood Estimation. Uses
 * the deprecated API since the new API for version 4.0 is not a fully documented final release.
 */
public class PoissonLikelihoodFunction implements DifferentiableMultivariateFunction
{
	private NonLinearFunction f;
	private float[] a, data;
	private int n;

	/**
	 * Initialise the function.
	 * <p>
	 * The input parameters must be the full parameters for the non-linear function. Only those
	 * parameters with gradient indices should be passed in to the functions to obtain the value (and gradient).
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
	public PoissonLikelihoodFunction(NonLinearFunction f, float[] a, float[] k, int n)
	{
		this.f = f;
		this.a = Arrays.copyOf(a, a.length);
		this.data = k;
		this.n = n;
	}

	/**
	 * Copy the variables into the appropriate parameter positions for the NonLinearFunction
	 * 
	 * @param variables
	 */
	private void initialiseFunction(double[] variables)
	{
		int[] gradientIndices = f.gradientIndices();
		for (int i = 0; i < gradientIndices.length; i++)
			a[gradientIndices[i]] = (float) variables[i];
		f.initialise(a);
	}

	/**
	 * Compute the value. Returns positive infinity if the function evaluates to zero (or below) at any point in the
	 * observed values.
	 * 
	 * @param variables
	 *            The variables of the function
	 * @return The negative log likelihood
	 */
	public double value(double[] variables)
	{
		initialiseFunction(variables);

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
	 * @see org.apache.commons.math3.analysis.DifferentiableMultivariateFunction#partialDerivative(int)
	 */
	@Override
	public MultivariateFunction partialDerivative(int k)
	{
		// This is not required for the AbstractDifferentiableOptimizer classes
		return null;
	}

	/**
	 * Compute the value and gradient of the function. Returns positive infinity if the function evaluates to zero (or
	 * below) at any point in the observed values. In this case the gradient computed so far will be invalid.
	 * 
	 * @param variables
	 *            The variables of the function
	 * @param gradient
	 *            The gradient (must be equal length to the variables array)
	 * @return The negative log likelihood
	 */
	public double value(double[] variables, double[] gradient)
	{
		initialiseFunction(variables);

		// Compute the negative log-likelihood to be minimised
		// f(x) = l(x) - k * ln(l(x))
		// 
		// Since (k * ln(l(x)))' = (k * ln(l(x)))' * l'(x) 

		// f'(x) = l'(x) - (k/l(x) * l'(x))
		// f'(x) = l'(x) * (1 - k/l(x))

		double ll = 0;
		for (int j = 0; j < variables.length; j++)
			gradient[j] = 0;
		float[] dl_da = new float[variables.length];
		for (int i = 0; i < n; i++)
		{
			final double l = f.eval(i, dl_da);

			// Check for zero and return the worst likelihood score
			if (l <= 0)
			{
				// Since ln(0) -> -Infinity
				return Double.POSITIVE_INFINITY;
			}

			final double k = data[i];
			ll += l - k * Math.log(l);
			for (int j = 0; j < gradient.length; j++)
				gradient[j] += dl_da[j] - (dl_da[j] * k / l);
		}
		return ll;
	}

	/**
	 * Compute the value and gradient of the function at observed value i. Returns positive infinity if the function
	 * evaluates to zero (or below) at the observed value. In this case the gradient computed so far will
	 * be invalid.
	 * 
	 * @param variables
	 *            The variables of the function
	 * @param gradient
	 *            The gradient (must be equal length to the variables array)
	 * @param i
	 *            Observed value i
	 * @return The negative log likelihood
	 */
	public double value(double[] variables, double[] gradient, int i)
	{
		initialiseFunction(variables);

		for (int j = 0; j < variables.length; j++)
			gradient[j] = 0;
		float[] dl_da = new float[variables.length];
		final double l = f.eval(i, dl_da);

		// Check for zero and return the worst likelihood score
		if (l <= 0)
		{
			// Since ln(0) -> -Infinity
			return Double.POSITIVE_INFINITY;
		}

		final double k = data[i];
		for (int j = 0; j < gradient.length; j++)
			gradient[j] = dl_da[j] - (dl_da[j] * k / l);
		return l - k * Math.log(l);

	}

	/**
	 * Compute the value of the function at observed value i. Returns positive infinity if the function
	 * evaluates to zero (or below) at the observed value.
	 * 
	 * @param variables
	 *            The variables of the function
	 * @param i
	 *            Observed value i
	 * @return The negative log likelihood
	 */
	public double value(double[] variables, int i)
	{
		initialiseFunction(variables);

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
	 * @see org.apache.commons.math3.analysis.DifferentiableMultivariateFunction#gradient()
	 */
	@Override
	public MultivariateVectorFunction gradient()
	{
		return new MultivariateVectorFunction()
		{
			public double[] value(double[] point) throws IllegalArgumentException
			{
				double[] gradient = new double[point.length];
				if (PoissonLikelihoodFunction.this.value(point, gradient) == Double.POSITIVE_INFINITY)
					return new double[point.length]; // No gradient if the function is invalid
				return gradient;
			}
		};
	}
}