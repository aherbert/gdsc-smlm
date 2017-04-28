package gdsc.smlm.function.gaussian.erf;

import org.apache.commons.math3.util.Pair;

import gdsc.smlm.function.gaussian.Gaussian2DFunction;

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
 * Evaluates a 2-dimensional Gaussian function for a single peak.
 */
public class SingleCircularErfGaussian2DFunction extends SingleFreeCircularErfGaussian2DFunction
{
	private static final int[] gradientIndices;
	static
	{
		gradientIndices = createGradientIndices(1, new SingleCircularErfGaussian2DFunction(1, 1, 0));
	}

	/**
	 * Constructor.
	 *
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 * @param maxy
	 *            The maximum y value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 * @param derivativeOrder
	 *            Set to the order of partial derivatives required
	 */
	public SingleCircularErfGaussian2DFunction(int maxx, int maxy, int derivativeOrder)
	{
		super(maxx, maxy, derivativeOrder);
	}

	@Override
	public Gaussian2DFunction create(int derivativeOrder)
	{
		if (derivativeOrder == this.derivativeOrder)
			return this;
		return new SingleCircularErfGaussian2DFunction(maxx, maxy, derivativeOrder);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.fitting.function.NonLinearFunction#initialise(double[])
	 */
	public void initialise(double[] a)
	{
		I0 = a[Gaussian2DFunction.SIGNAL];
		bg = a[Gaussian2DFunction.BACKGROUND];
		// Pre-compute the offset by 0.5
		x0 = a[Gaussian2DFunction.X_POSITION] + 0.5;
		y0 = a[Gaussian2DFunction.Y_POSITION] + 0.5;
		sx0 = sy0 = a[Gaussian2DFunction.X_SD];

		// We can pre-compute part of the derivatives for position and sd in arrays 
		// since the Gaussian is XY separable

		if (derivativeOrder == (byte) 2)
		{
			createSecondOrderTables(deltaEx, dudx, dudsx, d2udx2, d2udsx2, x0, sx0);
			createSecondOrderTables(deltaEy, dudy, dudsy, d2udy2, d2udsy2, y0, sy0);
		}
		else if (derivativeOrder == (byte) 1)
		{
			createFirstOrderTables(deltaEx, dudx, dudsx, x0, sx0);
			createFirstOrderTables(deltaEy, dudy, dudsy, y0, sy0);
		}
		else
		{
			createDeltaETable(deltaEx, x0, sx0);
			createDeltaETable(deltaEy, y0, sy0);
		}
	}

	/**
	 * Evaluates an 2-dimensional Gaussian function for a single peak.
	 * 
	 * @param i
	 *            Input predictor
	 * @return The Gaussian value
	 * 
	 * @see gdsc.fitting.function.NonLinearFunction#eval(int)
	 */
	public double eval(final int i)
	{
		// Unpack the predictor into the dimensions
		final int y = i / maxx;
		final int x = i % maxx;

		return bg + I0 * deltaEx[x] * deltaEy[y];
	}

	/**
	 * Evaluates an 2-dimensional Gaussian function for a single peak.
	 * 
	 * @param i
	 *            Input predictor
	 * @param duda
	 *            Partial gradient of function with respect to each coefficient
	 * @return The predicted value
	 * 
	 * @see gdsc.smlm.function.NonLinearFunction#eval(int, double[])
	 */
	public double eval(final int i, final double[] duda)
	{
		// Unpack the predictor into the dimensions
		final int y = i / maxx;
		final int x = i % maxx;

		// Return in order of Gaussian2DFunction.createGradientIndices().
		// Use pre-computed gradients
		duda[0] = 1.0;
		duda[1] = deltaEx[x] * deltaEy[y];
		duda[2] = dudx[x] * deltaEy[y];
		duda[3] = dudy[y] * deltaEx[x];
		duda[4] = dudsx[x] * deltaEy[y] + dudsy[y] * deltaEx[x];

		return bg + I0 * duda[1];
	}

	/**
	 * Evaluates an 2-dimensional Gaussian function for a single peak.
	 * 
	 * @param i
	 *            Input predictor
	 * @param duda
	 *            Partial first gradient of function with respect to each coefficient
	 * @param d2uda2
	 *            Partial second gradient of function with respect to each coefficient
	 * @return The predicted value
	 */
	public double eval(final int i, final double[] duda, final double[] d2uda2)
	{
		// Unpack the predictor into the dimensions
		final int y = i / maxx;
		final int x = i % maxx;

		// Return in order of Gaussian2DFunction.createGradientIndices().
		// Use pre-computed gradients
		duda[0] = 1.0;
		duda[1] = deltaEx[x] * deltaEy[y];
		duda[2] = dudx[x] * deltaEy[y];
		duda[3] = dudy[y] * deltaEx[x];
		duda[4] = dudsx[x] * deltaEy[y] + dudsy[y] * deltaEx[x];
		d2uda2[0] = 0;
		d2uda2[1] = 0;
		d2uda2[2] = d2udx2[x] * deltaEy[y];
		d2uda2[3] = d2udy2[y] * deltaEx[x];
		d2uda2[4] = d2udsx2[x] * deltaEy[y] + d2udsy2[y] * deltaEx[x];

		return bg + I0 * duda[1];
	}

	// Support for ExtendedNonLinear Function. This can take advantage of x/y iteration.
	@Override
	public double[][] computeJacobian(double[] variables)
	{
		initialise(variables);

		final int n = maxx * maxy;
		final double[][] jacobian = new double[n][];

		for (int y = 0, i = 0; y < maxy; y++)
		{
			final double dudy = this.dudy[y];
			final double deltaEy = this.deltaEy[y];
			final double dudsy = this.dudsy[y];
			for (int x = 0; x < maxx; x++, i++)
			{
				final double[] duda = new double[variables.length];
				duda[0] = 1.0;
				duda[1] = deltaEx[x] * deltaEy;
				duda[2] = dudx[x] * deltaEy;
				duda[3] = dudy * deltaEx[x];
				duda[4] = dudsx[x] * deltaEy + dudsy * deltaEx[x];
				jacobian[i] = duda;
			}
		}

		return jacobian;
	}

	// Support for ExtendedNonLinear Function. This can take advantage of x/y iteration.
	@Override
	public Pair<double[], double[][]> computeValuesAndJacobian(double[] variables)
	{
		initialise(variables);

		final int n = maxx * maxy;
		final double[][] jacobian = new double[n][];
		final double[] values = new double[n];

		for (int y = 0, i = 0; y < maxy; y++)
		{
			final double dudy = this.dudy[y];
			final double deltaEy = this.deltaEy[y];
			final double dudsy = this.dudsy[y];
			for (int x = 0; x < maxx; x++, i++)
			{
				final double[] duda = new double[variables.length];
				duda[0] = 1.0;
				duda[1] = deltaEx[x] * deltaEy;
				duda[2] = dudx[x] * deltaEy;
				duda[3] = dudy * deltaEx[x];
				duda[4] = dudsx[x] * deltaEy + dudsy * deltaEx[x];
				values[i] = bg + I0 * duda[1];
				jacobian[i] = duda;
			}
		}

		return new Pair<double[], double[][]>(values, jacobian);
	}

	@Override
	public int getNPeaks()
	{
		return 1;
	}

	@Override
	public boolean evaluatesBackground()
	{
		return true;
	}

	@Override
	public boolean evaluatesSignal()
	{
		return true;
	}

	@Override
	public boolean evaluatesShape()
	{
		return false;
	}

	@Override
	public boolean evaluatesPosition()
	{
		return true;
	}

	@Override
	public boolean evaluatesSD0()
	{
		return true;
	}

	@Override
	public boolean evaluatesSD1()
	{
		return false;
	}

	@Override
	public int getParametersPerPeak()
	{
		return 4;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.fitting.function.NonLinearFunction#gradientIndices()
	 */
	public int[] gradientIndices()
	{
		return gradientIndices;
	}
}
