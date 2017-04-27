package gdsc.smlm.function.gaussian.erf;

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
	protected final static double ONE_OVER_ROOT2 = 1.0 / Math.sqrt(2);
	protected final static double ONE_OVER_ROOT2PI = 1.0 / Math.sqrt(2 * Math.PI);

	private static final int[] gradientIndices;
	static
	{
		gradientIndices = createGradientIndices(1, new SingleFreeCircularErfGaussian2DFunction(1, 1, true));
	}

	// Required for the PSF
	protected final double[] deltaEx, deltaEy;
	protected double bg, I0, x0, y0, sx0, sy0;

	// Required for the gradients
	protected double[] expx, expy;
	protected double oneOverSxRoot2pi, oneOverSyRoot2pi;

	/**
	 * Constructor.
	 *
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 * @param maxy
	 *            The maximum y value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 * @param noGradients
	 *            Set to true if the gradients are not required
	 */
	public SingleFreeCircularErfGaussian2DFunction(int maxx, int maxy, boolean noGradients)
	{
		super(maxx, maxy, noGradients);
		deltaEx = new double[this.maxx];
		deltaEy = new double[this.maxy];
		if (noGradients)
			return;
		expx = new double[this.maxx + 1];
		expy = new double[this.maxy + 1];
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
		
		// TODO - initialise for 0,1,2 derivatives
		
		// We can store the derivatives for position and sd in arrays 
		// since the Gaussian is XY separable
		
		if (noGradients)
		{
			// Only create the tables for the PSF
			createDeltaETable(deltaEx, x0, sx0);
			createDeltaETable(deltaEy, y0, sy0);
		}
		else
		{
			// Create the tables for the PSF and the gradients
			createTables(deltaEx, expx, x0, sx0);
			createTables(deltaEy, expy, y0, sy0);
			oneOverSxRoot2pi = ONE_OVER_ROOT2PI / sx0;
			oneOverSyRoot2pi = ONE_OVER_ROOT2PI / sy0;
		}
	}

	/**
	 * Creates the delta E array. This is the sum of the Gaussian function using the error function for each of the
	 * pixels from 0 to n. Also compute the value of the Gaussian function at the pixel boundaries.
	 *
	 * @param deltaE
	 *            the delta E (difference between the error function at the start and end of each pixel)
	 * @param exp
	 *            the Gaussian exponential function at each pixel boundary
	 * @param u
	 *            the mean of the Gaussian
	 * @param s
	 *            the standard deviation of the Gaussian
	 * @return the double[]
	 */
	private void createDeltaETable(double[] deltaE, double u, double s)
	{
		// Used for the denominator: sqrt(2) * s
		final double f = ONE_OVER_ROOT2 / s;

		// Note: The paper by Smith, et al computes the integral for the kth pixel centred at (x,y)
		// If x=u then the Erf will be evaluated at x-u+0.5 - x-u-0.5 => integral from -0.5 to 0.5.
		// This code sets the first pixel at (0,0).

		// All computations for pixel k (=(x,y)) that require the exponential can use x,y indices for the 
		// lower boundary value and x+1,y+1 indices for the upper value.

		// The first position:
		// Offset x by the position and get the pixel lower bound.
		// (x - u - 0.5) with x=0 and u offset by +0.5
		double x = -u;
		double erf = 0.5 * Erf.erf(x * f);
		for (int i = 0, n = deltaE.length; i < n; i++)
		{
			x += 1.0;
			final double erf2 = 0.5 * Erf.erf(x * f);
			deltaE[i] = erf2 - erf;
			erf = erf2;
		}
	}

	/**
	 * Creates the delta E array. This is the sum of the Gaussian function using the error function for each of the
	 * pixels from 0 to n. Also compute the value of the Gaussian function at the pixel boundaries.
	 *
	 * @param deltaE
	 *            the delta E (difference between the error function at the start and end of each pixel)
	 * @param exp
	 *            the Gaussian exponential function at each pixel boundary
	 * @param u
	 *            the mean of the Gaussian
	 * @param s
	 *            the standard deviation of the Gaussian
	 * @return the double[]
	 */
	private void createTables(double[] deltaE, double[] exp, double u, double s)
	{
		// Used for the denominator: sqrt(2) * s
		final double f = ONE_OVER_ROOT2 / s;
		// Used for the denominator: 2 * s^2
		final double f2 = 0.5 / (s * s);

		// Note: The paper by Smith, et al computes the integral for the kth pixel centred at (x,y)
		// If x=u then the Erf will be evaluated at x-u+0.5 - x-u-0.5 => integral from -0.5 to 0.5.
		// This code sets the first pixel at (0,0).

		// All computations for pixel k (=(x,y)) that require the exponential can use x,y indices for the 
		// lower boundary value and x+1,y+1 indices for the upper value.

		// The first position:
		// Offset x by the position and get the pixel lower bound.
		// (x - u - 0.5) with x=0 and u offset by +0.5
		double x = -u;
		double erf = 0.5 * Erf.erf(x * f);
		exp[0] = FastMath.exp(-(x * x * f2));
		for (int i = 0, n = deltaE.length; i < n;)
		{
			x += 1.0;
			final double erf2 = 0.5 * Erf.erf(x * f);
			deltaE[i] = erf2 - erf;
			erf = erf2;
			exp[++i] = FastMath.exp(-(x * x * f2));
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

		// Use pre-computed gradients
		
		duda[Gaussian2DFunction.X_POSITION] = I0 * oneOverSxRoot2pi * (expx[x + 1] - expx[x]) * deltaEy[y];
		duda[Gaussian2DFunction.Y_POSITION] = I0 * oneOverSyRoot2pi * (expy[y + 1] - expy[y]) * deltaEx[x];
		duda[Gaussian2DFunction.SIGNAL] = deltaEx[x] * deltaEy[y];
		duda[Gaussian2DFunction.BACKGROUND] = 1.0;

		duda[Gaussian2DFunction.X_SD] = I0 * du_ds(deltaEx, oneOverSxRoot2pi, x, x - x0, sx0, expx);
		duda[Gaussian2DFunction.Y_SD] = I0 * du_ds(deltaEy, oneOverSyRoot2pi, y, y - y0, sy0, expy);

		return bg + I0 * duda[Gaussian2DFunction.SIGNAL];
	}

	/**
	 * Compute the first derivative with respect to the standard deviation. This does not include the intensity of the
	 * Gaussian.
	 *
	 * @param deltaE
	 *            the delta E array
	 * @param oneOverSRoot2pi
	 *            the one over S root 2 pi constant
	 * @param x
	 *            the pixel index
	 * @param lx
	 *            the lower x offset (x-u-0.5)
	 * @param s
	 *            the standard deviation
	 * @param exp
	 *            the Gaussian exponential array
	 * @return the first derivative
	 */
	protected static double du_ds(double[] deltaE, double oneOverSRoot2pi, int x, double lx, double s, double[] exp)
	{
		// Static to allow JVM optimisation
		final double ux = lx + 1;
		// Compute G21(xk)
		return deltaE[x] * (oneOverSRoot2pi / s) * (lx * exp[x] - ux * exp[x + 1]);
	}

	/**
	 * Compute the second derivative with respect to the standard deviation. This does not include the intensity of the
	 * Gaussian.
	 *
	 * @param deltaE
	 *            the delta E array
	 * @param oneOverSRoot2pi
	 *            the one over S root 2 pi constant
	 * @param x
	 *            the pixel index
	 * @param lx
	 *            the lower x offset (x-u-0.5)
	 * @param s
	 *            the standard deviation
	 * @param exp
	 *            the Gaussian exponential array
	 * @return the second derivative
	 */
	protected static double d2u_ds2(double[] deltaE, double oneOverSRoot2pi, int x, double lx, double s, double[] exp)
	{
		// Static to allow JVM optimisation
		double ux = lx + 1;
		// Compute G31(xk)
		final double s2 = s * s;
		final double G31 = (oneOverSRoot2pi / s2) * (lx * exp[x] - ux * exp[x + 1]);
		// Compute G53(xk)
		lx = lx * lx * lx;
		ux = ux * ux * ux;
		final double G53 = (oneOverSRoot2pi / (s2 * s2)) * (lx * exp[x] - ux * exp[x + 1]);
		return deltaE[x] * (G53 - 2 * G31);
	}

	// Support for ExtendedNonLinear Function. This can take advantage of x/y iteration.
	@Override
	public Pair<double[], double[][]> computeValuesAndJacobian(double[] variables)
	{
		initialise(variables);

		final int n = maxx * maxy;
		final double[][] jacobian = new double[n][];
		final double[] values = new double[n];

		for (int i = 0; i < n; ++i)
		{
			// Assume linear X from 0..N
			final double[] dyda = new double[variables.length];
			values[i] = eval(i, dyda);
			jacobian[i] = dyda;
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
