package gdsc.smlm.fitting.nonlinear.gradient;

import gdsc.smlm.function.Gradient1Procedure;
import gdsc.smlm.function.GradientFunction;

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
 * Calculates the Hessian matrix (the square matrix of second-order partial derivatives of a function)
 * and the scaled gradient vector of the function's partial first derivatives with respect to the parameters.
 * This is used within the Levenberg-Marquardt method to fit a nonlinear model with coefficients (a) for a
 * set of data points (x, y).
 * <p>
 * Note that the Hessian matrix is scaled by 1/2 and the gradient vector is scaled by -1/2 for convenience in solving
 * the non-linear model. See Numerical Recipes in C++, 2nd Ed. Equation 15.5.8 for Nonlinear Models.
 */
public class LSQGradientProcedureMatrix implements Gradient1Procedure
{
	protected final double[] y;
	protected final GradientFunction func;
	
	/**
	 * The number of gradients
	 */
	public final int n;
	/**
	 * The scaled Hessian curvature matrix (size n*n)
	 */
	public final double[][] alpha;
	/**
	 * The scaled gradient vector of the function's partial first derivatives with respect to the parameters
	 * (size n)
	 */
	public final double[] beta;
	/**
	 * The sum-of-squares value for the fit
	 */
	public double ssx;

	protected int yi;
	protected boolean isNanGradients;

	/**
	 * @param y
	 *            Data to fit
	 * @param func
	 *            Gradient function
	 */
	public LSQGradientProcedureMatrix(final double[] y, final GradientFunction func)
	{
		this.y = y;
		this.func = func;
		this.n = func.getNumberOfGradients();
		alpha = new double[n][n];
		beta = new double[n];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient1Procedure#execute(double, double[])
	 */
	public void execute(double value, double[] dy_da)
	{
		final double dy = y[yi++] - value;

		// Compute:
		// - the scaled Hessian matrix (the square matrix of second-order partial derivatives of a function; 
		//   that is, it describes the local curvature of a function of many variables.)
		// - the scaled gradient vector of the function's partial first derivatives with respect to the parameters

		for (int j = 0; j < n; j++)
		{
			final double wgt = dy_da[j];

			for (int k = 0; k <= j; k++)
				alpha[j][k] += wgt * dy_da[k];

			beta[j] += wgt * dy;
		}

		ssx += dy * dy;
	}

	/**
	 * Evaluate the function and compute the sum-of-squares and the curvature matrix.
	 * <p>
	 * A call to {@link #isNaNGradients()} will indicate if the gradients were invalid.
	 *
	 * @param a
	 *            Set of coefficients for the function
	 */
	public void run(final double[] a)
	{
		initialise();
		func.initialise(a);
		func.forEach(this);
		finish();
	}

	protected void initialise()
	{
		ssx = 0;
		yi = 0;
		for (int i = 0; i < n; i++)
		{
			beta[i] = 0;
			for (int j = 0; j <= i; j++)
				alpha[i][j] = 0;
		}
	}

	protected void finish()
	{
		// Generate symmetric matrix
		for (int i = 0; i < n - 1; i++)
			for (int j = i + 1; j < n; j++)
				alpha[i][j] = alpha[j][i];
		isNanGradients = checkGradients();
	}

	protected boolean checkGradients()
	{
		for (int i = 0; i < n; i++)
		{
			if (Double.isNaN(beta[i]))
				return true;
			for (int j = 0; j <= i; j++)
				if (Double.isNaN(alpha[i][j]))
					return true;
		}
		return false;
	}

	/**
	 * @return True if the last calculation produced gradients with NaN values
	 */
	public boolean isNaNGradients()
	{
		return isNanGradients;
	}
}
