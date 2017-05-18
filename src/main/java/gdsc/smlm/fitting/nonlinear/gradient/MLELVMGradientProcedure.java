package gdsc.smlm.fitting.nonlinear.gradient;

import gdsc.smlm.function.Gradient1Function;

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
 * Calculates the scaled Hessian matrix (the square matrix of second-order partial derivatives of a function)
 * and the scaled gradient vector of the function's partial first derivatives with respect to the parameters.
 * This is used within the Levenberg-Marquardt method to fit a nonlinear model with coefficients (a) for a
 * set of data points (x, y).
 * <p>
 * This calculator computes a modified Chi-squared expression to perform Maximum Likelihood Estimation assuming Poisson
 * model. See Laurence & Chromy (2010) Efficient maximum likelihood estimator. Nature Methods 7, 338-339. The input data
 * must be Poisson distributed for this to be relevant.
 */
public class MLELVMGradientProcedure extends LSQLVMGradientProcedure
{
	/**
	 * @param y
	 *            Data to fit (must be positive)
	 * @param func
	 *            Gradient function
	 */
	public MLELVMGradientProcedure(final double[] y, final Gradient1Function func)
	{
		super(y, func);
		// We could check that y is positive ...
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient1Procedure#execute(double, double[])
	 */
	public void execute(double fi, double[] dfi_da)
	{
		++yi;
		// Function must produce a strictly positive output.
		// ---
		// The code provided in Laurence & Chromy (2010) Nature Methods 7, 338-339, SI
		// effectively ignores any function value below zero. This could lead to a 
		// situation where the best chisq value can be achieved by setting the output
		// function to produce 0 for all evaluations.
		// Optimally the function should be bounded to always produce a positive number.
		// ---
		if (fi > 0.0)
		{
			final double xi = y[yi];

			// We assume y[i] is positive
			if (xi <= 0.0)
			{
				value += fi;
				for (int k = 0; k < n; k++)
				{
					beta[k] -= dfi_da[k];
				}
			}
			else
			{
				value += (fi - xi - xi * Math.log(fi / xi));
				final double xi_fi2 = xi / fi / fi;
				final double e = 1 - (xi / fi);
				for (int k = 0, i = 0; k < n; k++)
				{
					beta[k] -= e * dfi_da[k];
					final double w = dfi_da[k] * xi_fi2;
					for (int l = 0; l <= k; l++)
						alpha[i++] += w * dfi_da[l];
				}
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.ValueProcedure#execute(double)
	 */
	public void execute(double fi)
	{
		++yi;
		// Function must produce a strictly positive output.
		if (fi > 0.0)
		{
			final double xi = y[yi];

			// We assume y[i] is positive
			if (xi <= 0.0)
			{
				value += fi;
			}
			else
			{
				value += (fi - xi - xi * Math.log(fi / xi));
			}
		}
	}

	@Override
	protected void finishGradient()
	{
		// Move the factor of 2 to the end
		value *= 2;
	}

	@Override
	protected void finishValue()
	{
		// Move the factor of 2 to the end
		value *= 2;
	}
}
