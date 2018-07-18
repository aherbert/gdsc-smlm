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

/**
 * Calculates the Hessian matrix (the square matrix of second-order partial derivatives of a function)
 * and the gradient vector of the function's partial first derivatives with respect to the parameters.
 * This is used within the Levenberg-Marquardt method to fit a nonlinear model with coefficients (a) for a
 * set of data points (x, y).
 * <p>
 * This calculator computes a modified Chi-squared expression to perform Maximum Likelihood Estimation assuming Poisson
 * model. See Laurence & Chromy (2010) Efficient maximum likelihood estimator. Nature Methods 7, 338-339. The input data
 * must be Poisson distributed for this to be relevant.
 */
public class MLEGradientCalculator3 extends MLEGradientCalculator
{
	/**
	 * Instantiates a new MLE gradient calculator.
	 */
	public MLEGradientCalculator3()
	{
		super(3);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.MLEGradientCalculator#zero(double[][], double[])
	 */
	@Override
	protected void zero(final double[][] alpha, final double[] beta)
	{
		alpha[0][0] = 0;
		alpha[1][0] = 0;
		alpha[1][1] = 0;
		alpha[2][0] = 0;
		alpha[2][1] = 0;
		alpha[2][2] = 0;

		beta[0] = 0;
		beta[1] = 0;
		beta[2] = 0;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.MLEGradientCalculator#compute(double[][], double[], double[], double,
	 * double)
	 */
	@Override
	protected void compute(final double[][] alpha, final double[] beta, final double[] dfi_da, final double fi,
			final double xi)
	{
		final double xi_fi = xi / fi;
		final double xi_fi2 = xi_fi / fi;
		final double e = 1 - (xi_fi);

		alpha[0][0] += dfi_da[0] * xi_fi2 * dfi_da[0];
		double w;
		w = dfi_da[1] * xi_fi2;
		alpha[1][0] += w * dfi_da[0];
		alpha[1][1] += w * dfi_da[1];
		w = dfi_da[2] * xi_fi2;
		alpha[2][0] += w * dfi_da[0];
		alpha[2][1] += w * dfi_da[1];
		alpha[2][2] += w * dfi_da[2];

		beta[0] -= e * dfi_da[0];
		beta[1] -= e * dfi_da[1];
		beta[2] -= e * dfi_da[2];
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.MLEGradientCalculator#symmetric(double[][])
	 */
	@Override
	protected void symmetric(final double[][] alpha)
	{
		alpha[0][1] = alpha[1][0];
		alpha[0][2] = alpha[2][0];
		alpha[1][2] = alpha[2][1];
	}
}
