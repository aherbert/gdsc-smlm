package gdsc.smlm.fitting.nonlinear.gradient;

import gdsc.smlm.function.Gradient2Function;

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
 * Computes the Jacobian matrix of the partial derivatives, dFi/dxj, for all n parameters. dFi is the first partial
 * derivative of the log likelihood function with respect to parameter i. dFi/dxj is the first partial derivative of dFi
 * with respect to parameter j.
 * <p>
 * Ref: Smith et al, (2010). Fast, single-molecule localisation that achieves theoretically minimum uncertainty.
 * Nature Methods 7, 373-375 (supplementary note), Eq. 12.
 */
public class FastMLEJacobianGradient2Procedure extends FastMLEGradient2Procedure
{
	/**
	 * The Jacobian matrix of the partial derivatives, dFi/dxj, for all n parameters. dFi is the first partial
	 * derivative of the log likelihood function with respect to parameter i. dFi/dxj is the first partial derivative of
	 * dFi with respect to parameter j.
	 */
	public final double[] J;

	/**
	 * @param x
	 *            Data to fit (must be positive, i.e. the value of a Poisson process)
	 * @param func
	 *            Gradient function (must produce a strictly positive value, i.e. the mean of a Poisson process)
	 */
	public FastMLEJacobianGradient2Procedure(final double[] x, final Gradient2Function func)
	{
		super(x, func);
		J = new double[n * n];
	}

	@Override
	public void computeSecondDerivative(double[] a)
	{
		super.computeSecondDerivative(a);
		
		// Note:
		// This class is an attempt to extend the method to implement correct multi-dimensional 
		// Newton-Raphson root finding (see Numerical Recipes in C++, 2nd Ed, page 385, function mnewt).
		// This involves computing the Jacobian matrix from the vector of gradient functions.
		// However JUnit gradient tests show the Jacobian is just a copy of the second partial 
		// derivative for each parameter and so LU decomposition fails with a singular matrix. 
		
		// Finish the Jacobian:
		// This is just a copy of the second partial derivative for each parameter
		// This passes JUnit gradient tests.  
		for (int i = 0, index = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				J[index++] = d2[j];
	}
}
