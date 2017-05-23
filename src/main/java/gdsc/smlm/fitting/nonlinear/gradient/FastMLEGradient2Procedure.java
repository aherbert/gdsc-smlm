package gdsc.smlm.fitting.nonlinear.gradient;

import java.util.Arrays;

import gdsc.smlm.function.Gradient1Procedure;
import gdsc.smlm.function.Gradient2Function;
import gdsc.smlm.function.Gradient2Procedure;
import gdsc.smlm.function.PoissonCalculator;
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
 * <p>
 * Ref: Smith et al, (2010). Fast, single-molecule localisation that achieves theoretically minimum uncertainty.
 * Nature Methods 7, 373-375 (supplementary note), Eq. 12.
 */
public class FastMLEGradient2Procedure implements ValueProcedure, Gradient1Procedure, Gradient2Procedure
{
	protected final double[] x;
	protected final Gradient2Function func;
	protected PoissonCalculator poissonCalculator = null;

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

	/**
	 * The value of the function. This is updated by calls to {@link #computeValue(double[])},
	 * {@link #computeFirstDerivative(double[])}, {@link #computeSecondDerivative(double[])}
	 */
	public final double[] u;

	/**
	 * @param x
	 *            Data to fit (must be positive, i.e. the value of a Poisson process)
	 * @param func
	 *            Gradient function (must produce a strictly positive value, i.e. the mean of a Poisson process)
	 */
	public FastMLEGradient2Procedure(final double[] x, final Gradient2Function func)
	{
		this.x = x;
		this.u = new double[x.length];
		this.func = func;
		this.n = func.getNumberOfGradients();
		d1 = new double[n];
		d2 = new double[n];
	}

	/**
	 * Calculates the first and second derivative of the Poisson log likelihood with respect to each parameter.
	 *
	 * @param a
	 *            Set of coefficients for the function
	 */
	public void computeSecondDerivative(final double[] a)
	{
		k = 0;
		reset2();
		func.initialise2(a);
		func.forEach((Gradient2Procedure) this);
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
		u[k] = uk;
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
		u[k] = uk;
		final double xk = x[k++];
		final double xk_uk_minus1 = xk / uk - 1.0;
		for (int i = 0; i < n; i++)
		{
			d1[i] += duk_dt[i] * xk_uk_minus1;
		}
	}

	/**
	 * Compute the value of the function.
	 *
	 * @param a
	 *            the a
	 * @return the double[]
	 */
	public double[] computeValue(final double[] a)
	{
		k = 0;
		func.initialise0(a);
		func.forEach((ValueProcedure) this);
		return u;
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
		u[k++] = uk;
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
		computeValue(a);
		return computeLogLikelihood();
	}

	/**
	 * Calculates the Poisson log likelihood using the last value of the function.
	 *
	 * @return the Poisson log likelihood
	 */
	public double computeLogLikelihood()
	{
		return getPoissonCalculator().logLikelihood(u);
	}

	/**
	 * Gets the poisson calculator, creating using the values x if necessary.
	 *
	 * @return the poisson calculator
	 */
	private PoissonCalculator getPoissonCalculator()
	{
		if (poissonCalculator == null)
			poissonCalculator = new PoissonCalculator(x);
		return poissonCalculator;
	}

	/**
	 * Calculates the Poisson log likelihood ratio.
	 *
	 * @param a
	 *            Set of coefficients for the function
	 * @return the Poisson log likelihood ratio
	 */
	public double computeLogLikelihoodRatio(final double[] a)
	{
		computeValue(a);
		return computeLogLikelihoodRatio();
	}

	/**
	 * Calculates the Poisson log likelihood ratio using the last value of the function.
	 *
	 * @return the Poisson log likelihood ratio
	 */
	public double computeLogLikelihoodRatio()
	{
		return getPoissonCalculator().getLogLikelihoodRatio(computeLogLikelihood());
	}

	/**
	 * @return True if the last update calculation produced gradients with NaN values
	 */
	public boolean isNaNGradients()
	{
		for (int i = n; i-- > 0;)
		{
			if (Double.isNaN(d1[i]))
				return true;
			if (Double.isNaN(d2[i]))
				return true;
		}
		return false;
	}
}
