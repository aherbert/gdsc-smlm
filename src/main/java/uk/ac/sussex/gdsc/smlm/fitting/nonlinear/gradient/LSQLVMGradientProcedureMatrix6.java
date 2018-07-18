/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient;

import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;

/**
 * Calculates the Hessian matrix (the square matrix of second-order partial derivatives of a function)
 * and the scaled gradient vector of the function's partial first derivatives with respect to the parameters.
 * This is used within the Levenberg-Marquardt method to fit a nonlinear model with coefficients (a) for a
 * set of data points (x, y).
 * <p>
 * Note that the Hessian matrix is scaled by 1/2 and the gradient vector is scaled by -1/2 for convenience in solving
 * the non-linear model. See Numerical Recipes in C++, 2nd Ed. Equation 15.5.8 for Nonlinear Models.
 */
public class LSQLVMGradientProcedureMatrix6 extends LSQLVMGradientProcedureMatrix
{
	/**
	 * @param y
	 *            Data to fit
	 * @param b
	 *            Baseline pre-computed y-values
	 * @param func
	 *            Gradient function
	 */
	public LSQLVMGradientProcedureMatrix6(final double[] y, final double[] b, final Gradient1Function func)
	{
		super(y, b, func);
		if (n != 6)
			throw new IllegalArgumentException("Function must compute 6 gradients");
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.function.Gradient1Procedure#execute(double, double[])
	 */
	@Override
	public void execute(double value, double[] dy_da)
	{
		final double dy = y[++yi] - value;

		alpha[0][0] += dy_da[0] * dy_da[0];
		alpha[1][0] += dy_da[1] * dy_da[0];
		alpha[1][1] += dy_da[1] * dy_da[1];
		alpha[2][0] += dy_da[2] * dy_da[0];
		alpha[2][1] += dy_da[2] * dy_da[1];
		alpha[2][2] += dy_da[2] * dy_da[2];
		alpha[3][0] += dy_da[3] * dy_da[0];
		alpha[3][1] += dy_da[3] * dy_da[1];
		alpha[3][2] += dy_da[3] * dy_da[2];
		alpha[3][3] += dy_da[3] * dy_da[3];
		alpha[4][0] += dy_da[4] * dy_da[0];
		alpha[4][1] += dy_da[4] * dy_da[1];
		alpha[4][2] += dy_da[4] * dy_da[2];
		alpha[4][3] += dy_da[4] * dy_da[3];
		alpha[4][4] += dy_da[4] * dy_da[4];
		alpha[5][0] += dy_da[5] * dy_da[0];
		alpha[5][1] += dy_da[5] * dy_da[1];
		alpha[5][2] += dy_da[5] * dy_da[2];
		alpha[5][3] += dy_da[5] * dy_da[3];
		alpha[5][4] += dy_da[5] * dy_da[4];
		alpha[5][5] += dy_da[5] * dy_da[5];

		beta[0] += dy_da[0] * dy;
		beta[1] += dy_da[1] * dy;
		beta[2] += dy_da[2] * dy;
		beta[3] += dy_da[3] * dy;
		beta[4] += dy_da[4] * dy;
		beta[5] += dy_da[5] * dy;

		this.value += dy * dy;
	}

	@Override
	protected void initialiseGradient()
	{
		alpha[0][0] = 0;
		alpha[1][0] = 0;
		alpha[1][1] = 0;
		alpha[2][0] = 0;
		alpha[2][1] = 0;
		alpha[2][2] = 0;
		alpha[3][0] = 0;
		alpha[3][1] = 0;
		alpha[3][2] = 0;
		alpha[3][3] = 0;
		alpha[4][0] = 0;
		alpha[4][1] = 0;
		alpha[4][2] = 0;
		alpha[4][3] = 0;
		alpha[4][4] = 0;
		alpha[5][0] = 0;
		alpha[5][1] = 0;
		alpha[5][2] = 0;
		alpha[5][3] = 0;
		alpha[5][4] = 0;
		alpha[5][5] = 0;
		beta[0] = 0;
		beta[1] = 0;
		beta[2] = 0;
		beta[3] = 0;
		beta[4] = 0;
		beta[5] = 0;
	}

	@Override
	protected void finishGradient()
	{
		alpha[0][1] = alpha[1][0];
		alpha[0][2] = alpha[2][0];
		alpha[0][3] = alpha[3][0];
		alpha[0][4] = alpha[4][0];
		alpha[0][5] = alpha[5][0];
		alpha[1][2] = alpha[2][1];
		alpha[1][3] = alpha[3][1];
		alpha[1][4] = alpha[4][1];
		alpha[1][5] = alpha[5][1];
		alpha[2][3] = alpha[3][2];
		alpha[2][4] = alpha[4][2];
		alpha[2][5] = alpha[5][2];
		alpha[3][4] = alpha[4][3];
		alpha[3][5] = alpha[5][3];
		alpha[4][5] = alpha[5][4];
	}
}
