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
package gdsc.smlm.function;

import gdsc.smlm.utils.Pair;

/**
 * Wrap the NonLinearFunction to remove the parameters that are fixed from the evaluation methods
 */
public class NonLinearFunctionWrapper implements ExtendedNonLinearFunction
{
	private final NonLinearFunction fun;
	private final double[] a;
	private final int n;
	private final int[] gradientIndices;

	/**
	 * Create a new instance using the full set of parameters for the function and the number of points the function
	 * evaluates. The parameters that are not within the function gradient indices array will be fixed.
	 *
	 * @param fun
	 *            The function
	 * @param a
	 *            The parameters
	 * @param n
	 *            The number of data points to evaluate
	 */
	public NonLinearFunctionWrapper(NonLinearFunction fun, double[] a, int n)
	{
		this.fun = fun;
		this.a = a.clone();
		this.n = n;
		// This wrapper will evaluate all the indices that are not fixed
		gradientIndices = new int[fun.getNumberOfGradients()];
		for (int i = 0; i < gradientIndices.length; i++)
			gradientIndices[i] = i;
	}

	/**
	 * Set the predictor coefficients for the function that are not fixed (i.e. those corresponding to the gradient
	 * indices in the wrapped function). The fixed coefficients are set in the constructor.
	 *
	 * @see gdsc.smlm.function.NonLinearFunction#initialise(double[])
	 */
	@Override
	public void initialise(double[] variables)
	{
		int[] gradientIndices = fun.gradientIndices();
		for (int i = 0; i < gradientIndices.length; i++)
			a[gradientIndices[i]] = variables[i];
		fun.initialise(a);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.fitting.function.NonLinearFunction#gradientIndices()
	 */
	@Override
	public int[] gradientIndices()
	{
		return gradientIndices;
	}

	@Override
	public int getNumberOfGradients()
	{
		return gradientIndices.length;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.fitting.function.NonLinearFunction#eval(int, double[])
	 */
	@Override
	public double eval(int x, double[] dyda)
	{
		return fun.eval(x, dyda);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.fitting.function.NonLinearFunction#eval(int)
	 */
	@Override
	public double eval(int x)
	{
		return fun.eval(x);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.fitting.function.NonLinearFunction#eval(int, double[], double[])
	 */
	@Override
	public double eval(int x, double[] dyda, double[] w)
	{
		return fun.eval(x, dyda, w);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.function.NonLinearFunction#evalw(int, double[])
	 */
	@Override
	public double evalw(int x, double[] w)
	{
		return fun.eval(x, w);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.fitting.function.NonLinearFunction#canComputeWeights()
	 */
	@Override
	public boolean canComputeWeights()
	{
		return fun.canComputeWeights();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.function.ExtendedNonLinearFunction#computeValues(double[])
	 */
	@Override
	public double[] computeValues(double[] variables)
	{
		initialise(variables);
		final double[] values = new double[n];
		for (int i = 0; i < values.length; i++)
		{
			// Assume linear X from 0..N
			values[i] = fun.eval(i);
		}
		return values;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.function.ExtendedNonLinearFunction#computeJacobian(double[])
	 */
	@Override
	public double[][] computeJacobian(double[] variables)
	{
		initialise(variables);

		final double[][] jacobian = new double[n][];

		for (int i = 0; i < n; ++i)
		{
			// Assume linear X from 0..N
			final double[] dyda = new double[variables.length];
			fun.eval(i, dyda);
			jacobian[i] = dyda;
		}

		return jacobian;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.function.ExtendedNonLinearFunction#canComputeValuesAndJacobian()
	 */
	@Override
	public boolean canComputeValuesAndJacobian()
	{
		return true;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.function.ExtendedNonLinearFunction#computeValuesAndJacobian(double[])
	 */
	@Override
	public Pair<double[], double[][]> computeValuesAndJacobian(double[] variables)
	{
		initialise(variables);

		final double[][] jacobian = new double[n][];
		final double[] values = new double[n];

		for (int i = 0; i < n; ++i)
		{
			// Assume linear X from 0..N
			final double[] dyda = new double[variables.length];
			values[i] = fun.eval(i, dyda);
			jacobian[i] = dyda;
		}

		return new Pair<>(values, jacobian);
	}
}
