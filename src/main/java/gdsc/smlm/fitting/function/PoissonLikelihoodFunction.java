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
 * 
 * This is an adaption of the C-code contained in the CcpNmr Analysis Program:
 *   CCPN website (http://www.ccpn.ac.uk/). 
 * The CCPN code was based on Numerical Recipes. 
 *---------------------------------------------------------------------------*/

import java.util.Arrays;

import org.apache.commons.math3.analysis.DifferentiableMultivariateFunction;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.MultivariateVectorFunction;

/**
 * Wrapper for any function to allow use of the Apache Commons optimiser for Maximum Likelihood Estimation.
 * <p>
 * Uses the deprecated API since the new API for version 4.0 is not a fully documented final release.
 */
public class PoissonLikelihoodFunction implements DifferentiableMultivariateFunction
{
	private NonLinearFunction f;
	private float[] a, y;
	private int n;

	/**
	 * @param f
	 *            The function to be used to calculated the expected values
	 * @param a
	 *            The initial parameters for the function
	 * @param y
	 *            The observed values
	 * @param n
	 *            The number of observed values
	 */
	public PoissonLikelihoodFunction(NonLinearFunction f, float[] a, float[] y, int n)
	{
		this.f = f;
		this.a = Arrays.copyOf(a, a.length);
		this.y = y;
		this.n = n;
	}

	/**
	 * Copy the variables into the appropriate parameter positions for the NonLinearFunction
	 * 
	 * @param variables
	 * @return true if the function could evaluate zero
	 */
	private boolean initialiseFunction(double[] variables)
	{
		int[] gradientIndices = f.gradientIndices();
		for (int i = 0; i < gradientIndices.length; i++)
			a[gradientIndices[i]] = (float) variables[i];
		// Do not allow negative values to be computed (since log(-x) is undefined)
		boolean zero = true;

		//		// TODO - Remove dependency on Gaussian2DFunction
		//		// Pass in bounds to the function parameters that require restriction
		//		if (a[Gaussian2DFunction.AMPLITUDE] < 0)
		//		{
		//			zero = true;
		//			a[Gaussian2DFunction.AMPLITUDE] = 0;
		//		}
		//		if (a[Gaussian2DFunction.BACKGROUND] < 0)
		//		{
		//			zero = true;
		//			a[Gaussian2DFunction.BACKGROUND] = 0;
		//		}
		f.initialise(a);
		return zero;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.apache.commons.math3.analysis.MultivariateFunction#value(double[])
	 */
	public double value(double[] variables)
	{
		boolean zero = initialiseFunction(variables);

		// Compute the negative log-likelihood to be minimised
		double ll = 0;
		for (int i = 0; i < n; i++)
		{
			final double l = f.eval(i);

			// If the function was called with arguments that could be zero, check for zero
			// and return the worst likelihood score
			if (zero && l <= 0)
				// Since ln(0) -> -Infinity
				return Double.POSITIVE_INFINITY;

			final double k = y[i];
			ll += l - k * Math.log(l);
		}
		return ll;
	}

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
		boolean zero = initialiseFunction(variables);

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

			// If the function was called with arguments that could be zero, check for zero
			// and return the worst likelihood score
			if (zero && l <= 0)
			{
				// Since ln(0) -> -Infinity
				return Double.POSITIVE_INFINITY;
			}

			final double k = y[i];
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
		boolean zero = initialiseFunction(variables);

		// Compute the negative log-likelihood to be minimised
		// f(x) = l(x) - k * ln(l(x))
		// 
		// Since (k * ln(l(x)))' = (k * ln(l(x)))' * l'(x) 

		// f'(x) = l'(x) - (k/l(x) * l'(x))
		// f'(x) = l'(x) * (1 - k/l(x))

		double ll = 0;
		if (gradient == null)
			gradient = new double[variables.length];
		float[] dl_da = new float[variables.length];
		final double l = f.eval(i, dl_da);

		// If the function was called with arguments that could be zero, check for zero
		// and return the worst likelihood score
		if (zero && l <= 0)
		{
			// Since ln(0) -> -Infinity
			return Double.POSITIVE_INFINITY;
		}

		final double k = y[i];
		ll = l - k * Math.log(l);
		for (int j = 0; j < gradient.length; j++)
			gradient[j] = dl_da[j] - (dl_da[j] * k / l);
		return ll;

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
		boolean zero = initialiseFunction(variables);

		// Compute the negative log-likelihood to be minimised
		// TODO - Check if this is valid:
		// f(x) = l(x) - k * ln(l(x))
		// 
		// Since (k * ln(l(x)))' = (k * ln(l(x)))' * l'(x) 

		// f'(x) = l'(x) - (k/l(x) * l'(x))
		// f'(x) = l'(x) * (1 - k/l(x))

		double ll = 0;
		final double l = f.eval(i);

		// If the function was called with arguments that could be zero, check for zero
		// and return the worst likelihood score
		if (zero && l <= 0)
		{
			// Since ln(0) -> -Infinity
			return Double.POSITIVE_INFINITY;
		}

		final double k = y[i];
		ll += l - k * Math.log(l);
		return ll;

	}

	@Override
	public MultivariateVectorFunction gradient()
	{
		return new MultivariateVectorFunction()
		{
			public double[] value(double[] point) throws IllegalArgumentException
			{
				double[] gradient = new double[point.length];
				if (PoissonLikelihoodFunction.this.value(point, gradient) == Double.POSITIVE_INFINITY)
					return new double[point.length];
				return gradient;
			}
		};
	}
}