package gdsc.smlm.function;

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
public class NonLinearFunctionWrapper implements NonLinearFunction
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
		gradientIndices = new int[fun.gradientIndices().length];
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
	 * @see gdsc.smlm.fitting.function.NonLinearFunction#canComputeWeights()
	 */
	public boolean canComputeWeights()
	{
		return fun.canComputeWeights();
	}

	/**
	 * Helper method to return the value of all the data points
	 * 
	 * @param variables
	 * @return
	 */
	public double[] computeValue(double[] variables)
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

	/**
	 * Helper method to return the jacobian of the gradients of all the data points
	 * 
	 * @param variables
	 * @return
	 */
	public double[][] computeJacobian(double[] variables)
	{
		initialise(variables);

		final double[][] jacobian = new double[n][variables.length];
		final double[] dyda = new double[variables.length];

		for (int i = 0; i < jacobian.length; ++i)
		{
			//float y = gf.eval(x.get(i).intValue());
			// Assume linear X from 0..N
			fun.eval(i, dyda);

			// Differentiate with respect to each parameter:
			for (int j = 0; j < dyda.length; j++)
				jacobian[i][j] = dyda[j];
		}

		return jacobian;
	}
}