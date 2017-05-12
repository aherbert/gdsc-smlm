package gdsc.smlm.function;

import org.apache.commons.math3.util.Pair;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

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
	public int[] gradientIndices()
	{
		return gradientIndices;
	}


	public int getNumberOfGradients()
	{
		return gradientIndices.length;
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.function.NonLinearFunction#eval(int, double[])
	 */
	public double eval(int x, double[] dyda)
	{
		return fun.eval(x, dyda);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.function.NonLinearFunction#eval(int)
	 */
	public double eval(int x)
	{
		return fun.eval(x);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.function.NonLinearFunction#eval(int, double[], double[])
	 */
	public double eval(int x, double[] dyda, double[] w)
	{
		return fun.eval(x, dyda, w);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.NonLinearFunction#evalw(int, double[])
	 */
	public double evalw(int x, double[] w)
	{
		return fun.eval(x, w);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.function.NonLinearFunction#canComputeWeights()
	 */
	public boolean canComputeWeights()
	{
		return fun.canComputeWeights();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.ExtendedNonLinearFunction#computeValues(double[])
	 */
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
	public boolean canComputeValuesAndJacobian()
	{
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.ExtendedNonLinearFunction#computeValuesAndJacobian(double[])
	 */
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

		return new Pair<double[], double[][]>(values, jacobian);
	}
}
