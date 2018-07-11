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
package gdsc.smlm.fitting.nonlinear;

import gdsc.core.utils.Maths;
import gdsc.smlm.fitting.FunctionSolverType;
import gdsc.smlm.fitting.LSEFunctionSolver;
import gdsc.smlm.fitting.linear.EJMLLinearSolver;
import gdsc.smlm.function.GradientFunction;

/**
 * Abstract class with utility methods for the LSEFunctionSolver interface.
 */
public abstract class LSEBaseFunctionSolver extends BaseFunctionSolver implements LSEFunctionSolver
{
	protected double totalSumOfSquares = Double.NaN;

	/**
	 * Default constructor
	 *
	 * @throws NullPointerException
	 *             if the function is null
	 */
	public LSEBaseFunctionSolver(GradientFunction f)
	{
		super(FunctionSolverType.LSE, f);
	}

	@Override
	protected void preProcess()
	{
		totalSumOfSquares = Double.NaN;
	}

	/**
	 * Gets the total sum of squares.
	 *
	 * @param y
	 *            the y
	 * @return the total sum of squares
	 */
	public static double getTotalSumOfSquares(double[] y)
	{
		double sx = 0, ssx = 0;
		for (int i = y.length; i-- > 0;)
		{
			sx += y[i];
			ssx += y[i] * y[i];
		}
		final double sumOfSquares = ssx - (sx * sx) / (y.length);
		return sumOfSquares;
	}

	/**
	 * Compute the error
	 *
	 * @param value
	 * @param noise
	 * @param numberOfFittedPoints
	 * @param numberOfFittedParameters
	 * @return the error
	 */
	public static double getError(double value, double noise, int numberOfFittedPoints, int numberOfFittedParameters)
	{
		double error = value;

		// Divide by the uncertainty in the individual measurements yi to get the chi-squared
		if (noise > 0)
		{
			error /= numberOfFittedPoints * noise * noise;
		}

		// This updates the chi-squared value to the average error for a single fitted
		// point using the degrees of freedom (N-M)?
		// Note: This matches the mean squared error output from the MatLab fitting code.
		// If a noise estimate was provided for individual measurements then this will be the
		// reduced chi-square (see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2892436/)
		if (numberOfFittedPoints > numberOfFittedParameters)
			error /= (numberOfFittedPoints - numberOfFittedParameters);
		else
			error = 0;

		return error;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.fitting.LSEFunctionSolver#getTotalSumOfSquares()
	 */
	@Override
	public double getTotalSumOfSquares()
	{
		if (Double.isNaN(totalSumOfSquares) && lastY != null)
		{
			totalSumOfSquares = getTotalSumOfSquares(lastY);
		}
		return totalSumOfSquares;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.fitting.LSEFunctionSolver#getResidualSumOfSquares()
	 */
	@Override
	public double getResidualSumOfSquares()
	{
		return value;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.fitting.LSEFunctionSolver#getCoefficientOfDetermination()
	 */
	@Override
	public double getCoefficientOfDetermination()
	{
		return 1.0 - (value / getTotalSumOfSquares());
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.fitting.LSEFunctionSolver#getAdjustedCoefficientOfDetermination()
	 */
	@Override
	public double getAdjustedCoefficientOfDetermination()
	{
		return Maths.getAdjustedCoefficientOfDetermination(getResidualSumOfSquares(), getTotalSumOfSquares(),
				getNumberOfFittedPoints(), getNumberOfFittedParameters());
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.fitting.LSEFunctionSolver#getMeanSquaredError()
	 */
	@Override
	public double getMeanSquaredError()
	{
		return getResidualSumOfSquares() / (getNumberOfFittedPoints() - getNumberOfFittedParameters());
	}

	/**
	 * Compute the covariance matrix of the parameters of the function assuming a least squares fit of a Poisson
	 * process.
	 * <p>
	 * Uses the Mortensen formula (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 25.
	 * <p>
	 * The method involves inversion of a matrix and may fail.
	 *
	 * <pre>
	 * I = sum_i { Ei,a * Ei,b }
	 * E = sum_i { Ei * Ei,a * Ei,b }
	 *
	 * with
	 * i the number of data points fit using least squares using a function of n variable parameters
	 * Ei the expected value of the function at i
	 * Ei,a the gradient the function at i with respect to parameter a
	 * Ei,b the gradient the function at i with respect to parameter b
	 * </pre>
	 *
	 * @param I
	 *            the Iab matrix
	 * @param E
	 *            the Ei_Eia_Eib matrix
	 * @return the covariance matrix (or null)
	 */
	public static double[][] covariance(double[][] I, double[][] E)
	{
		int n = I.length;

		// Invert the matrix
		EJMLLinearSolver solver = EJMLLinearSolver.createForInversion(1e-2);
		if (!solver.invert(I))
			return null;

		// Note that I now refers to I^-1 in the Mortensen notation

		double[][] covar = new double[n][n];
		for (int a = 0; a < n; a++)
		{
			for (int b = 0; b < n; b++)
			{
				double v = 0;
				for (int ap = 0; ap < n; ap++)
				{
					for (int bp = 0; bp < n; bp++)
					{
						v += I[a][ap] * E[ap][bp] * I[bp][b];
					}
				}
				covar[a][b] = v;
			}
		}

		return covar;
	}

	/**
	 * Compute the variance of the parameters of the function assuming a least squares fit of a Poisson process.
	 * <p>
	 * Uses the Mortensen formula (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 25.
	 * <p>
	 * The method involves inversion of a matrix and may fail.
	 *
	 * <pre>
	 * I = sum_i { Ei,a * Ei,b }
	 * E = sum_i { Ei * Ei,a * Ei,b }
	 *
	 * with
	 * i the number of data points fit using least squares using a function of n variable parameters
	 * Ei the expected value of the function at i
	 * Ei,a the gradient the function at i with respect to parameter a
	 * Ei,b the gradient the function at i with respect to parameter b
	 * </pre>
	 *
	 * @param I
	 *            the Iab matrix
	 * @param E
	 *            the Ei_Eia_Eib matrix
	 * @return the variance (or null)
	 */
	public static double[] variance(double[][] I, double[][] E)
	{
		int n = I.length;

		// Invert the matrix
		EJMLLinearSolver solver = EJMLLinearSolver.createForInversion(1e-2);
		if (!solver.invert(I))
			return null;

		// Note that I now refers to I^-1 in the Mortensen notation

		double[] covar = new double[n];
		for (int a = 0; a < n; a++)
		{
			// Note: b==a as we only do the diagonal
			double v = 0;
			for (int ap = 0; ap < n; ap++)
			{
				for (int bp = 0; bp < n; bp++)
				{
					v += I[a][ap] * E[ap][bp] * I[bp][a];
				}
			}
			covar[a] = v;
		}

		return covar;
	}
}
