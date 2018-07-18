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
package gdsc.smlm.fitting.nonlinear.gradient;

import gdsc.smlm.function.NonLinearFunction;

/**
 * Calculates the Hessian matrix (the square matrix of second-order partial derivatives of a function)
 * and the gradient vector of the function's partial first derivatives with respect to the parameters.
 * This is used within the Levenberg-Marquardt method to fit a nonlinear model with coefficients (a) for a
 * set of data points (x, y).
 */
public class GradientCalculator4 extends GradientCalculator
{
	/**
	 * Instantiates a new gradient calculator for a 4x4 matrix.
	 */
	public GradientCalculator4()
	{
		super(4);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.fitting.model.GradientCalculator#findLinearised(int[], double[] double[], double[][], double[],
	 * gdsc.fitting.function.NonLinearFunction)
	 */
	@Override
	public double findLinearised(int[] x, double[] y, double[] a, double[][] alpha, double[] beta,
			NonLinearFunction func)
	{
		double ssx = 0;
		final double[] dy_da = new double[4];

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

		beta[0] = 0;
		beta[1] = 0;
		beta[2] = 0;
		beta[3] = 0;

		func.initialise(a);

		if (func.canComputeWeights())
		{
			final double[] w = new double[1];
			for (int i = 0; i < x.length; i++)
			{
				final double dy = y[i] - func.eval(x[i], dy_da, w);
				final double weight = getWeight(w[0]);

				alpha[0][0] += dy_da[0] * weight * dy_da[0];
				alpha[1][0] += dy_da[1] * weight * dy_da[0];
				alpha[1][1] += dy_da[1] * weight * dy_da[1];
				alpha[2][0] += dy_da[2] * weight * dy_da[0];
				alpha[2][1] += dy_da[2] * weight * dy_da[1];
				alpha[2][2] += dy_da[2] * weight * dy_da[2];
				alpha[3][0] += dy_da[3] * weight * dy_da[0];
				alpha[3][1] += dy_da[3] * weight * dy_da[1];
				alpha[3][2] += dy_da[3] * weight * dy_da[2];
				alpha[3][3] += dy_da[3] * weight * dy_da[3];

				beta[0] += dy_da[0] * weight * dy;
				beta[1] += dy_da[1] * weight * dy;
				beta[2] += dy_da[2] * weight * dy;
				beta[3] += dy_da[3] * weight * dy;

				ssx += dy * dy * weight;
			}
		}
		else
			for (int i = 0; i < x.length; i++)
			{
				final double dy = y[i] - func.eval(x[i], dy_da);

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

				beta[0] += dy_da[0] * dy;
				beta[1] += dy_da[1] * dy;
				beta[2] += dy_da[2] * dy;
				beta[3] += dy_da[3] * dy;

				ssx += dy * dy;
			}

		// Generate symmetric matrix
		alpha[0][1] = alpha[1][0];
		alpha[0][2] = alpha[2][0];
		alpha[0][3] = alpha[3][0];
		alpha[1][2] = alpha[2][1];
		alpha[1][3] = alpha[3][1];
		alpha[2][3] = alpha[3][2];

		return checkGradients(alpha, beta, nparams, ssx);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.fitting.nonlinear.gradient.GradientCalculator#findLinearised(int, double[] double[], double[][],
	 * double[], gdsc.fitting.function.NonLinearFunction)
	 */
	@Override
	public double findLinearised(int n, double[] y, double[] a, double[][] alpha, double[] beta, NonLinearFunction func)
	{
		double ssx = 0;
		final double[] dy_da = new double[4];

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

		beta[0] = 0;
		beta[1] = 0;
		beta[2] = 0;
		beta[3] = 0;

		func.initialise(a);

		if (func.canComputeWeights())
		{
			final double[] w = new double[1];
			for (int i = 0; i < n; i++)
			{
				final double dy = y[i] - func.eval(i, dy_da, w);
				final double weight = getWeight(w[0]);

				alpha[0][0] += dy_da[0] * weight * dy_da[0];
				alpha[1][0] += dy_da[1] * weight * dy_da[0];
				alpha[1][1] += dy_da[1] * weight * dy_da[1];
				alpha[2][0] += dy_da[2] * weight * dy_da[0];
				alpha[2][1] += dy_da[2] * weight * dy_da[1];
				alpha[2][2] += dy_da[2] * weight * dy_da[2];
				alpha[3][0] += dy_da[3] * weight * dy_da[0];
				alpha[3][1] += dy_da[3] * weight * dy_da[1];
				alpha[3][2] += dy_da[3] * weight * dy_da[2];
				alpha[3][3] += dy_da[3] * weight * dy_da[3];

				beta[0] += dy_da[0] * weight * dy;
				beta[1] += dy_da[1] * weight * dy;
				beta[2] += dy_da[2] * weight * dy;
				beta[3] += dy_da[3] * weight * dy;

				ssx += dy * dy * weight;
			}
		}
		else
			for (int i = 0; i < n; i++)
			{
				final double dy = y[i] - func.eval(i, dy_da);

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

				beta[0] += dy_da[0] * dy;
				beta[1] += dy_da[1] * dy;
				beta[2] += dy_da[2] * dy;
				beta[3] += dy_da[3] * dy;

				ssx += dy * dy;
			}

		// Generate symmetric matrix
		alpha[0][1] = alpha[1][0];
		alpha[0][2] = alpha[2][0];
		alpha[0][3] = alpha[3][0];
		alpha[1][2] = alpha[2][1];
		alpha[1][3] = alpha[3][1];
		alpha[2][3] = alpha[3][2];

		return checkGradients(alpha, beta, nparams, ssx);
	}
}
