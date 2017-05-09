package gdsc.smlm.fitting.nonlinear.gradient;

import java.util.Arrays;

import org.apache.commons.math3.special.Gamma;

import gdsc.smlm.function.Gradient1Procedure;
import gdsc.smlm.function.Gradient2Function;
import gdsc.smlm.function.Gradient2Procedure;
import gdsc.smlm.function.ValueProcedure;

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
 * Calculates the Newton-Raphson update vector for a Poisson process using the first and second partial derivatives.
 * 
 * @see Smith et al, (2010). Fast, single-molecule localisation that achieves theoretically minimum uncertainty.
 *      Nature Methods 7, 373-375 (supplementary note), Eq. 12.
 */
public class NewtonRaphsonGradient2Procedure implements ValueProcedure, Gradient1Procedure, Gradient2Procedure
{
	protected final double[] x;
	protected final Gradient2Function func;

	/**
	 * The number of gradients
	 */
	public final int n;
	/**
	 * The first derivative of the Poisson log likelihood with respect to each parameter
	 */
	public final double[] d1;
	/**
	 * The second derivative of the Poisson log likelihood with respect to each parameter
	 */
	public final double[] d2;

	/** Counter */
	protected int k;

	/** The log likelihood. */
	protected double ll;

	/**
	 * @param x
	 *            Data to fit (must be positive, i.e. the value of a Poisson process)
	 * @param func
	 *            Gradient function (must produce a strictly positive value, i.e. the mean of a Poisson process)
	 */
	public NewtonRaphsonGradient2Procedure(final double[] x, final Gradient2Function func)
	{
		this.x = x;
		this.func = func;
		this.n = func.getNumberOfGradients();
		d1 = new double[n];
		d2 = new double[n];
	}

	/**
	 * Calculates the Newton-Raphson update vector for a Poisson process.
	 *
	 * @param a
	 *            Set of coefficients for the function
	 * @return the update vector of the function's parameters
	 */
	public double[] computeUpdate(final double[] a)
	{
		k = 0;
		reset2();
		func.initialise2(a);
		func.forEach((Gradient2Procedure) this);
		return getUpdate();
	}

	/**
	 * Reset the first and second derivative vectors
	 */
	protected void reset2()
	{
		Arrays.fill(d1, 0);
		Arrays.fill(d2, 0);
	}

	/**
	 * Calculates the Newton-Raphson update vector for a Poisson process. Variables are named as per the Smith, et al
	 * (2010) paper.
	 * 
	 * {@inheritDoc}
	 * 
	 * @see gdsc.smlm.function.Gradient2Procedure#execute(double, double[], double[])
	 */
	public void execute(double uk, double[] duk_dt, double[] d2uk_dt2)
	{
		final double xk = x[k++];
		final double xk_uk_minus1 = xk / uk - 1.0;
		final double xk_uk2 = xk / (uk * uk);
		for (int i = 0; i < n; i++)
		{
			d1[i] += duk_dt[i] * xk_uk_minus1;
			d2[i] += d2uk_dt2[i] * xk_uk_minus1 - duk_dt[i] * duk_dt[i] * xk_uk2;
		}
	}

	/**
	 * Gets the update vector of the function's parameters (size n).
	 *
	 * @return the update vector
	 */
	public double[] getUpdate()
	{
		double[] update = new double[n];
		for (int i = 0; i < n; i++)
			update[i] = d1[i] / d2[i];
		return update;
	}

	/**
	 * Calculates the first derivative of the Poisson log likelihood with respect to each parameter.
	 *
	 * @param a
	 *            Set of coefficients for the function
	 * @return the first derivative of the Poisson log likelihood with respect to each parameter
	 */
	public double[] computeFirstDerivative(final double[] a)
	{
		k = 0;
		reset1();
		func.initialise1(a);
		func.forEach((Gradient1Procedure) this);
		return d1;
	}

	/**
	 * Reset the first derivative vector
	 */
	protected void reset1()
	{
		Arrays.fill(d1, 0);
	}

	/**
	 * Variables are named as per the Smith, et al (2010) paper.
	 * 
	 * {@inheritDoc}
	 * 
	 * @see gdsc.smlm.function.Gradient1Procedure#execute(double, double[])
	 */
	public void execute(double uk, double[] duk_dt)
	{
		final double xk = x[k++];
		final double xk_uk_minus1 = xk / uk - 1.0;
		for (int i = 0; i < n; i++)
		{
			d1[i] += duk_dt[i] * xk_uk_minus1;
		}
	}

	/**
	 * Calculates the Poisson log likelihood.
	 *
	 * @param a
	 *            Set of coefficients for the function
	 * @return the Poisson log likelihood
	 */
	public double computeLogLikelihood(final double[] a)
	{
		ll = 0;
		k = 0;
		func.initialise0(a);
		func.forEach((ValueProcedure) this);
		return ll;
	}

	/**
	 * Variables are named as per the Smith, et al (2010) paper.
	 * 
	 * {@inheritDoc}
	 * 
	 * @see gdsc.smlm.function.ValueProcedure#execute(double)
	 */
	public void execute(double uk)
	{
		ll += logLikelihood(uk, x[k++]);
	}

	/**
	 * Get the Poisson log likelihood of value x given the mean. The mean must be strictly positive. x must be positive.
	 *
	 * @param u
	 *            the mean
	 * @param x
	 *            the x
	 * @return the log likelihood
	 */
	public static double logLikelihood(double u, double x)
	{
		// X is not zero very often so skip this ...
		//if (x == 0)
		//	return -u;
		return x * Math.log(u) - u - logFactorial(x);
	}

	private static double logFactorial(double k)
	{
		if (k <= 1)
			return 0;
		return Gamma.logGamma(k + 1);
	}
}
