package gdsc.smlm.fitting.function;

import java.util.Arrays;

import org.apache.commons.math3.analysis.DifferentiableMultivariateVectorFunction;
import org.apache.commons.math3.analysis.MultivariateMatrixFunction;

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
 * 
 * This is an adaption of the C-code contained in the CcpNmr Analysis Program:
 * CCPN website (http://www.ccpn.ac.uk/)
 *---------------------------------------------------------------------------*/

/**
 * Wrapper for the 2-dimensional Gaussian function for a configured number of peaks to allow use of the Apache Commons
 * optimiser.
 * <p>
 * Uses the deprecated API since the new API for version 4.0 is not a fully documented final release.
 */
public class DifferentiableGaussian2DFunction implements DifferentiableMultivariateVectorFunction
{
	private Gaussian2DFunction gf;
	private double[] y = null;
	private float[] a;

	/**
	 * @param gf The gaussian function to be use to calculate the gradients
	 * @param a The initial parameters for the gaussian function
	 */
	public DifferentiableGaussian2DFunction(Gaussian2DFunction gf, float[] a)
	{
		this.gf = gf;
		this.a = Arrays.copyOf(a, a.length);
	}

	public void addData(int n, float[] y)
	{
		this.y = toDouble(n, y);
	}

	public double[] getY()
	{
		return y;
	}

	public double[] getWeights()
	{
		double[] w = new double[y.length];
		Arrays.fill(w, 1);
		return w;
	}

	// Adapted from http://commons.apache.org/proper/commons-math/userguide/optimization.html
	// Use the deprecated API since the new one is not yet documented.

	private double[][] jacobian(double[] variables)
	{
		initialiseFunction(variables);
		
		double[][] jacobian = new double[y.length][variables.length];
		float[] dyda = new float[variables.length];

		for (int i = 0; i < jacobian.length; ++i)
		{
			//float y = gf.eval(x.get(i).intValue());
			// Assume linear X from 0..N
			gf.eval(i, dyda);

			// Differentiate with respect to N:
			for (int j = 0; j < dyda.length; j++)
				jacobian[i][j] = dyda[j];
		}

		return jacobian;
	}

	private void initialiseFunction(double[] variables)
	{
		int[] gradientIndices = gf.gradientIndices();
		for (int i = 0; i < gradientIndices.length; i++)
			a[gradientIndices[i]] = (float)variables[i];
		gf.initialise(a);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.apache.commons.math3.analysis.MultivariateVectorFunction#value(double[])
	 */
	public double[] value(double[] variables)
	{
		initialiseFunction(variables);

		double[] values = new double[y.length];
		for (int i = 0; i < values.length; i++)
		{
			//values[i] = gf.eval(x.get(i).intValue());
			// Assume linear X from 0..N
			values[i] = gf.eval(i);
		}
		return values;
	}

	public static double[] toDouble(int n, float[] in)
	{
		double[] out = new double[n];
		for (int i = 0; i < out.length; i++)
			out[i] = in[i];
		return out;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.apache.commons.math3.analysis.DifferentiableMultivariateVectorFunction#jacobian()
	 */
	public MultivariateMatrixFunction jacobian()
	{
		return new MultivariateMatrixFunction()
		{
			public double[][] value(double[] variables)
			{
				return jacobian(variables);
			}
		};
	}
}
