package gdsc.smlm.fitting.nonlinear.gradient;

import java.util.Arrays;

import gdsc.smlm.function.Gradient1Function;
import gdsc.smlm.function.Gradient1Procedure;

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
 * Calculates the Fisher information matrix for a Poisson process.
 * 
 * @see Smith et al, (2010). Fast, single-molecule localisation that achieves theoretically minimum uncertainty.
 *      Nature Methods 7, 373-375 (supplementary note), Eq. 9.
 */
public class PoissonGradientProcedure implements Gradient1Procedure
{
	protected final Gradient1Function func;

	/**
	 * The number of gradients
	 */
	public final int n;

	/**
	 * Working space for the Fisher information matrix (size n*n)
	 */
	protected double[] data;

	/**
	 * @param func
	 *            Gradient function
	 */
	public PoissonGradientProcedure(final Gradient1Function func)
	{
		this.func = func;
		this.n = func.getNumberOfGradients();
	}

	/**
	 * Compute the Fisher information matrix
	 * 
	 * <pre>
	 * Iab = E [ ( d ln(L(x|p)) / da ) * ( d ln(L(x|p)) / db ) ]
	 * p = parameters
	 * x = obsereved values
	 * L(x|p) = likelihood of X given p
	 * E = expected value
	 * </pre>
	 * 
	 * Note that this is only a true Fisher information diagonal if the function returns the expected value for a
	 * Poisson process. In this case the equation reduces to:
	 * 
	 * <pre>
	 * Iaa = sum(i) (dYi da) * (dYi da) / Yi
	 * </pre>
	 * 
	 * See Smith et al, (2010). Fast, single-molecule localisation that achieves theoretically minimum uncertainty.
	 * Nature Methods 7, 373-375 (supplementary note), Eq. 9.
	 * <p>
	 * A call to {@link #isNaNGradients()} will indicate if the gradients were invalid.
	 *
	 * @param a
	 *            Set of coefficients for the function
	 */
	public void computeFisherInformation(final double[] a)
	{
		if (data == null)
			data = new double[n * (n + 1) / 2];
		else
			initialiseWorkingMatrix();
		func.initialise1(a);
		func.forEach((Gradient1Procedure) this);
	}

	/**
	 * Initialise for the computation of the Fisher information matrix. Zero the working matrix.
	 */
	protected void initialiseWorkingMatrix()
	{
		Arrays.fill(data, 0);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient1Procedure#execute(double, double[])
	 */
	public void execute(double value, double[] dy_da)
	{
		if (value > 0)
		{
			final double f = 1.0 / value;
			for (int j = 0, i = 0; j < n; j++)
			{
				final double wgt = f * dy_da[j];
				for (int k = 0; k <= j; k++)
				{
					data[i++] += wgt * dy_da[k];
				}
			}
		}
	}

	protected boolean checkGradients()
	{
		for (int i = 0, len = data.length; i < len; i++)
			if (Double.isNaN(data[i]))
				return true;
		return false;
	}

	/**
	 * Get the scaled Fisher Information matrix (size n*n).
	 *
	 * @return the alpha
	 */
	public double[][] getMatrix()
	{
		double[][] a = new double[n][n];
		getMatrix(a);
		return a;
	}

	/**
	 * Get the scaled Fisher Information matrix (size n*n) into the provided storage.
	 *
	 * @param matrix
	 *            the matrix
	 * @return the matrix
	 */
	public void getMatrix(double[][] matrix)
	{
		GradientProcedureHelper.getMatrix(data, matrix, n);
	}

	/**
	 * Get the scaled Fisher Information matrix (size n*n).
	 *
	 * @return the alpha
	 */
	public double[] getLinear()
	{
		double[] a = new double[n * n];
		getLinear(a);
		return a;
	}

	/**
	 * Get the scaled Fisher Information matrix (size n*n) into the provided storage.
	 *
	 * @param matrix
	 *            the matrix
	 * @return the linear matrix
	 */
	public void getLinear(double[] matrix)
	{
		GradientProcedureHelper.getMatrix(data, matrix, n);
	}

	/**
	 * @return True if the last calculation produced gradients with NaN values
	 */
	public boolean isNaNGradients()
	{
		return checkGradients();
	}
}
