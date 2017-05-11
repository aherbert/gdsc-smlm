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
public class MLELVMGradientProcedureB4 extends MLELVMGradientProcedure4
{
	protected final double[] b;

	/**
	 * @param y
	 *            Data to fit (must be positive)
	 * @param b
	 *            Baseline pre-computed y-values
	 * @param func
	 *            Gradient function
	 */
	public MLELVMGradientProcedureB4(final double[] y, final double[] b, final Gradient1Function func)
	{
		super(y, func);
		this.b = b;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient1Procedure#execute(double, double[])
	 */
	public void execute(double fi, double[] dfi_da)
	{
		// Add the baseline to the function value
		super.execute(fi + b[yi + 1], dfi_da);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.ValueProcedure#execute(double)
	 */
	public void execute(double fi)
	{
		// Add the baseline to the function value
		super.execute(fi + b[yi + 1]);
	}
}
