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
 * Note that the Hessian matrix is scaled by 1/2 and the gradient vector is scaled by -1/2 for convenience in solving
 * the non-linear model. See Numerical Recipes in C++, 2nd Ed. Equation 15.5.8 for Nonlinear Models.
 */
public abstract class BaseLSQLVMGradientProcedure extends LVMGradientProcedure
{
	/**
	 * @param y
	 *            Data to fit
	 * @param func
	 *            Gradient function
	 */
	public BaseLSQLVMGradientProcedure(final double[] y, final Gradient1Function func)
	{
		super(y, func);
	}

	/**
	 * @param y
	 *            Data to fit
	 * @param b
	 *            Baseline pre-computed y-values
	 * @param func
	 *            Gradient function
	 */
	public BaseLSQLVMGradientProcedure(final double[] y, final double[] b, final Gradient1Function func)
	{
		super(y, b, func);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.ValueProcedure#execute(double)
	 */
	public void execute(double value)
	{
		// Produce a sum-of-squares
		final double dy = y[++yi] - value;
		this.value += dy * dy;
	}

	@Override
	protected void initialiseValue()
	{
		// Do nothing		
	}

	@Override
	protected void finishValue()
	{
		// Do nothing		
	}
}
