package gdsc.smlm.function.gaussian.erf;

//import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Pair;

import gdsc.smlm.function.Erf;
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
public class SingleFreeCircularErfGaussian2DFunction extends ErfGaussian2DFunction
{
	private static final int[] gradientIndices;
	static
	{
		gradientIndices = createGradientIndices(1, new SingleFreeCircularErfGaussian2DFunction(1, 1, 0));
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
	public SingleFreeCircularErfGaussian2DFunction(int maxx, int maxy, int derivativeOrder)
	{
		super(maxx, maxy, derivativeOrder);
	}

	@Override
	protected void createArrays()
	{
		if (derivativeOrder <= 0)
			return;
		dudx = new double[this.maxx];
		dudy = new double[this.maxy];
		dudsx = new double[this.maxx];
		dudsy = new double[this.maxy];
		if (derivativeOrder <= 1)
			return;
		d2udx2 = new double[this.maxx];
		d2udy2 = new double[this.maxy];
		d2udsx2 = new double[this.maxx];
		d2udsy2 = new double[this.maxy];
	}

	@Override
	public Gaussian2DFunction create(int derivativeOrder)
	{
		if (derivativeOrder == this.derivativeOrder)
			return this;
		return new SingleFreeCircularErfGaussian2DFunction(maxx, maxy, derivativeOrder);
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
		sx0 = a[Gaussian2DFunction.X_SD];
		sy0 = a[Gaussian2DFunction.Y_SD];

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
	 * Creates the delta E array. This is the sum of the Gaussian function using the error function for each of the
	 * pixels from 0 to n.
	 *
	 * @param deltaE0
	 *            the delta E (difference between the error function at the start and end of each pixel)
	 * @param u0
	 *            the mean of the Gaussian for dimension 0
	 * @param s0
	 *            the standard deviation of the Gaussian for dimension 0
	 */
	protected void createDeltaETable(double[] deltaE0, double u0, double s0)
	{
		final double one_oversSqrt2 = ONE_OVER_ROOT2 / s0;

		// Note: The paper by Smith, et al computes the integral for the kth pixel centred at (x,y)
		// If x=u then the Erf will be evaluated at x-u+0.5 - x-u-0.5 => integral from -0.5 to 0.5.
		// This code sets the first pixel at (0,0).

		// All computations for pixel k (=(x,y)) that require the exponential can use x,y indices for the 
		// lower boundary value and x+1,y+1 indices for the upper value.

		// The first position:
		// Offset x by the position and get the pixel lower bound.
		// (x - u - 0.5) with x=0 and u offset by +0.5
		double x = -u0;
		double erf = 0.5 * Erf.erf(x * one_oversSqrt2);
		for (int i = 0, n = deltaE0.length; i < n; i++)
		{
			x += 1.0;
			final double erf2 = 0.5 * Erf.erf(x * one_oversSqrt2);
			deltaE0[i] = erf2 - erf;
			erf = erf2;
		}
	}

	/**
	 * Creates the first order derivatives.
	 *
	 * @param deltaE0
	 *            the delta E for dimension 0 (difference between the error function at the start and end of each pixel)
	 * @param dudx0
	 *            the first order x derivative for dimension 0
	 * @param duds0
	 *            the first order s derivative for dimension 0
	 * @param u0
	 *            the mean of the Gaussian for dimension 0
	 * @param s0
	 *            the standard deviation of the Gaussian for dimension 0
	 */
	protected void createFirstOrderTables(double[] deltaE0, double[] dudx0, double[] duds0, double u0, double s0)
	{
		final double one_oversSqrt2 = ONE_OVER_ROOT2 / s0;
		final double one_2ss = 0.5 / (s0 * s0);
		final double I0_sSqrt2pi = I0 * ONE_OVER_ROOT2PI / s0;

		// Note: The paper by Smith, et al computes the integral for the kth pixel centred at (x,y)
		// If x=u then the Erf will be evaluated at x-u+0.5 - x-u-0.5 => integral from -0.5 to 0.5.
		// This code sets the first pixel at (0,0).

		// All computations for pixel k (=(x,y)) that require the exponential can use x,y indices for the 
		// lower boundary value and x+1,y+1 indices for the upper value.

		// The first position:
		// Offset x by the position and get the pixel lower bound.
		// (x - u - 0.5) with x=0 and u offset by +0.5
		double x = -u0;
		double erf = 0.5 * Erf.erf(x * one_oversSqrt2);
		double exp = FastMath.exp(-(x * x * one_2ss));
		for (int i = 0, n = deltaE0.length; i < n; i++)
		{
			x += 1.0;
			final double erf2 = 0.5 * Erf.erf(x * one_oversSqrt2);
			deltaE0[i] = erf2 - erf;
			erf = erf2;

			final double exp2 = FastMath.exp(-(x * x * one_2ss));
			final double pre = I0_sSqrt2pi;
			dudx0[i] = pre * (exp - exp2);
			// Compute: I0 * G21(xk)
			duds0[i] = (pre / s0) * ((x - 1) * exp - x * exp2);

			exp = exp2;
		}
	}

	/**
	 * Creates the first and second order derivatives.
	 *
	 * @param deltaE0
	 *            the delta E for dimension 0 (difference between the error function at the start and end of each pixel)
	 * @param dudx0
	 *            the first order x derivative for dimension 0
	 * @param duds0
	 *            the first order s derivative for dimension 0
	 * @param d2udx02
	 *            the second order x derivative for dimension 0
	 * @param d2uds02
	 *            the second order s derivative for dimension 0
	 * @param u0
	 *            the mean of the Gaussian for dimension 0
	 * @param s0
	 *            the standard deviation of the Gaussian for dimension 0
	 */
	protected void createSecondOrderTables(double[] deltaE0, double[] dudx0, double[] duds0, double[] d2udx02,
			double[] d2uds02, double u0, double s0)
	{
		final double ss = s0 * s0;
		final double one_oversSqrt2 = ONE_OVER_ROOT2 / s0;
		final double one_2ss = 0.5 / ss;
		final double one_sSqrt2pi = ONE_OVER_ROOT2PI / s0;
		final double I0_sSqrt2pi = I0 * one_sSqrt2pi;
		final double one_sssSqrt2pi = one_sSqrt2pi / ss;
		final double one_sssssSqrt2pi = one_sssSqrt2pi / ss;

		// Note: The paper by Smith, et al computes the integral for the kth pixel centred at (x,y)
		// If x=u then the Erf will be evaluated at x-u+0.5 - x-u-0.5 => integral from -0.5 to 0.5.
		// This code sets the first pixel at (0,0).

		// All computations for pixel k (=(x,y)) that require the exponential can use x,y indices for the 
		// lower boundary value and x+1,y+1 indices for the upper value.

		// The first position:
		// Offset x by the position and get the pixel lower bound.
		// (x - u - 0.5) with x=0 and u offset by +0.5
		double x = -u0;
		double erf = 0.5 * Erf.erf(x * one_oversSqrt2);
		double exp = FastMath.exp(-(x * x * one_2ss));
		for (int i = 0, n = deltaE0.length; i < n; i++)
		{
			x += 1.0;
			final double erf2 = 0.5 * Erf.erf(x * one_oversSqrt2);
			deltaE0[i] = erf2 - erf;
			erf = erf2;

			final double exp2 = FastMath.exp(-(x * x * one_2ss));
			double lx = x - 1;
			final double pre = I0_sSqrt2pi; // * deltaE1[i];
			dudx0[i] = pre * (exp - exp2);
			// Compute: I0 * G21(xk)
			final double pre2 = (lx * exp - x * exp2);
			duds0[i] = (pre / s0) * pre2;

			// Second derivatives
			d2udx02[i] = (pre / ss) * pre2;

			// Compute G31(xk)
			final double G31 = one_sssSqrt2pi * pre2;

			// Compute G53(xk)
			lx = lx * lx * lx;
			final double ux = x * x * x;
			final double G53 = one_sssssSqrt2pi * (lx * exp - ux * exp2);
			d2uds02[i] = I0 * (G53 - 2 * G31);

			exp = exp2;
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
		duda[4] = dudsx[x] * deltaEy[y];
		duda[5] = dudsy[y] * deltaEx[x];

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
		duda[4] = dudsx[x] * deltaEy[y];
		duda[5] = dudsy[y] * deltaEx[x];
		d2uda2[0] = 0;
		d2uda2[1] = 0;
		d2uda2[2] = d2udx2[x] * deltaEy[y];
		d2uda2[3] = d2udy2[y] * deltaEx[x];
		d2uda2[4] = d2udsx2[x] * deltaEy[y];
		d2uda2[5] = d2udsy2[y] * deltaEx[x];

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
				duda[4] = dudsx[x] * deltaEy;
				duda[5] = dudsy * deltaEx[x];
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
				duda[4] = dudsx[x] * deltaEy;
				duda[5] = dudsy * deltaEx[x];
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
		return true;
	}

	@Override
	public int getParametersPerPeak()
	{
		return 5;
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
