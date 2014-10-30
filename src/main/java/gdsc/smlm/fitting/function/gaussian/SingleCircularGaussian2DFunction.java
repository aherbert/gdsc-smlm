package gdsc.smlm.fitting.function.gaussian;

import org.apache.commons.math3.util.FastMath;

import gdsc.smlm.fitting.function.Gaussian2DFunction;

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
 * Evaluates an 2-dimensional Gaussian function for a single peak.
 * <p>
 * The single parameter x in the {@link #eval(int, double[])} function is assumed to be a linear index into 2-dimensional
 * data. The dimensions of the data must be specified to allow unpacking to coordinates.
 * <p>
 * Data should be packed in descending dimension order, e.g. Y,X : Index for [x,y] = MaxX*y + x.
 */
public class SingleCircularGaussian2DFunction extends Gaussian2DFunction
{
	private static int[] gradientIndices;
	static
	{
		gradientIndices = createGradientIndices(1, new SingleCircularGaussian2DFunction(1));
	}

	protected double background;
	protected double amplitude;
	protected double x0pos;
	protected double x1pos;

	protected double aa;
	protected double aa2;
	protected double ax;

	/**
	 * Constructor
	 * 
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 */
	public SingleCircularGaussian2DFunction(int maxx)
	{
		super(maxx);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.fitting.function.NonLinearFunction#initialise(double[])
	 */
	public void initialise(double[] a)
	{
		background = a[BACKGROUND];
		amplitude = a[AMPLITUDE];
		x0pos = a[X_POSITION];
		x1pos = a[Y_POSITION];

		final double sx = a[X_SD];
		final double sx2 = sx * sx;
		final double sx3 = sx2 * sx;

		// All prefactors are negated since the Gaussian uses the exponential to the negative:
		// A * exp( -( a(x-x0)^2 + 2b(x-x0)(y-y0) + c(y-y0)^2 ) )

		aa = -0.5 / sx2;
		aa2 = -2.0 * aa;

		// For the x-width gradient
		ax = 1.0 / sx3;
	}

	/**
	 * Produce an output predicted value for a given set of input
	 * predictors (x) and coefficients (a).
	 * <p>
	 * Evaluates an 2-dimensional elliptical Gaussian function for a single peak.
	 * <p>
	 * The first coefficient is the Gaussian background level (B). The coefficients are then packed for each peak:
	 * Amplitude; Angle; position[N]; sd[N]. Amplitude (A) is the height of the Gaussian. Angle (r) is the rotation
	 * angle of the ellipse. Position (x,y) is the position of the Gaussian in each of the N-dimensions. SD (sx,sy) is
	 * the standard deviation in each of the N-dimensions.
	 * <p>
	 * The equation per peak is:<br/>
	 * y_peak = A * exp( -( a(x-x0)^2 + c(y-y0)^2 ) )<br/>
	 * Where: <br/>
	 * a = 1/(2*sx^2) <br/>
	 * c = 1/(2*sy^2)
	 * 
	 * @param x
	 *            Input predictor
	 * @param dyda
	 *            Partial gradient of function with respect to each coefficient
	 * @return The predicted value
	 * 
	 * @see gdsc.smlm.fitting.function.NonLinearFunction#eval(int, double[])
	 */
	public double eval(final int x, final double[] dyda)
	{
		// First parameter is the background level 
		dyda[0] = 1; // Gradient for a constant background is 1

		// Unpack the predictor into the dimensions
		final int x1 = x / maxx;
		final int x0 = x % maxx;

		return background + gaussian(x0, x1, dyda);
	}

	private double gaussian(final int x0, final int x1, final double[] dy_da)
	{
		final double h = amplitude;

		final double dx = x0 - x0pos;
		final double dy = x1 - x1pos;
		final double dx2dy2 = dx * dx + dy * dy;

		//final double y = (double) (h * FastMath.exp(aa * (dx2dy2)));

		// Calculate gradients
		//dy_da[1] = y / h;

		dy_da[1] = FastMath.exp(aa * (dx2dy2));
		final double y = h * dy_da[1];
		final double yaa2 = y * aa2;
		dy_da[2] = yaa2 * dx;
		dy_da[3] = yaa2 * dy;

		dy_da[4] = y * (ax * (dx2dy2));

		return y;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.fitting.function.NonLinearFunction#eval(int)
	 */
	public double eval(final int x)
	{
		// Unpack the predictor into the dimensions
		final int x1 = x / maxx;
		final int x0 = x % maxx;

		final double dx = x0 - x0pos;
		final double dy = x1 - x1pos;

		return background + amplitude * FastMath.exp(aa * (dx * dx + dy * dy));
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
	public boolean evaluatesAmplitude()
	{
		return true;
	}

	@Override
	public boolean evaluatesAngle()
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
