package gdsc.smlm.function.gaussian;

import org.apache.commons.math3.util.FastMath;

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
 * Evaluates an 2-dimensional elliptical Gaussian function for a single peak.
 * <p>
 * The single parameter x in the {@link #eval(int, double[])} function is assumed to be a linear index into
 * 2-dimensional data. The dimensions of the data must be specified to allow unpacking to coordinates.
 * <p>
 * Data should be packed in descending dimension order, e.g. Y,X : Index for [x,y] = MaxX*y + x.
 */
public class SingleEllipticalGaussian2DFunction extends Gaussian2DFunction
{
	private static int[] gradientIndices;
	static
	{
		gradientIndices = createGradientIndices(1, new SingleEllipticalGaussian2DFunction(1));
	}

	protected double background;
	protected double x0pos;
	protected double x1pos;

	protected double n;
	protected double height;
	protected double aa;
	protected double bb;
	protected double cc;
	protected double aa2;
	protected double bb2;
	protected double cc2;
	protected double nx;
	protected double ax;
	protected double bx;
	protected double cx;
	protected double ny;
	protected double ay;
	protected double by;
	protected double cy;

	/**
	 * Constructor
	 * 
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 */
	public SingleEllipticalGaussian2DFunction(int maxx)
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
		x0pos = a[X_POSITION];
		x1pos = a[Y_POSITION];

		// Precalculate multiplication factors
		final double theta = a[ANGLE];
		final double sx = a[X_SD];
		final double sy = a[Y_SD];
		final double sx2 = sx * sx;
		final double sy2 = sy * sy;
		final double sx3 = sx2 * sx;
		final double sy3 = sy2 * sy;
		final double cosSqt = Math.cos(theta) * Math.cos(theta);
		final double sinSqt = Math.sin(theta) * Math.sin(theta);
		final double sincost = Math.sin(theta) * Math.cos(theta);
		final double sin2t = Math.sin(2 * theta);
		final double cos2t = Math.cos(2 * theta);

		n = ONE_OVER_TWO_PI / (sx * sy);
		height = a[SIGNAL] * n;

		// All prefactors are negated since the Gaussian uses the exponential to the negative:
		// (A/2*pi*sx*sy) * exp( -( a(x-x0)^2 + 2b(x-x0)(y-y0) + c(y-y0)^2 ) )

		aa = -0.5 * (cosSqt / sx2 + sinSqt / sy2);
		bb = -0.25 * (-sin2t / sx2 + sin2t / sy2);
		cc = -0.5 * (sinSqt / sx2 + cosSqt / sy2);

		// For the angle gradient
		aa2 = -(-sincost / sx2 + sincost / sy2);
		bb2 = -0.5 * (-cos2t / sx2 + cos2t / sy2);
		cc2 = -(sincost / sx2 - sincost / sy2);

		// For the x-width gradient
		nx = -1 / sx;
		ax = cosSqt / sx3;
		bx = -0.5 * sin2t / sx3;
		cx = sinSqt / sx3;

		// For the y-width gradient
		ny = -1 / sy;
		ay = sinSqt / sy3;
		by = 0.5 * sin2t / sy3;
		cy = cosSqt / sy3;
	}

	/**
	 * Produce an output predicted value for a given set of input
	 * predictors (x) and coefficients (a).
	 * <p>
	 * Evaluates an 2-dimensional elliptical Gaussian function for a single peak.
	 * <p>
	 * The first coefficient is the Gaussian background level (B). The coefficients are then packed for each peak:
	 * Amplitude; Angle; position[N]; sd[N]. Amplitude (A) is the volume of the Gaussian. Angle (r) is the rotation
	 * angle of the ellipse. Position (x,y) is the position of the Gaussian in each of the N-dimensions. SD (sx,sy) is
	 * the standard deviation in each of the N-dimensions.
	 * <p>
	 * The equation per peak is:<br/>
	 * y_peak = A/(2*pi*sx*sy) * exp( -( a(x-x0)^2 + 2b(x-x0)(y-y0) + c(y-y0)^2 ) )<br/>
	 * Where: <br/>
	 * a = cos(r)^2/(2*sx^2) + sin(r)^2 /(2*sy^2) <br/>
	 * b = -sin(2r)^2/(4*sx^2) + sin(2r)^2/(4*sy^2) <br/>
	 * c = sin(r)^2/(2*sx^2) + cos(r)^2/(2*sy^2)
	 * 
	 * @param x
	 *            Input predictor
	 * @param dyda
	 *            Partial gradient of function with respect to each coefficient
	 * @return The predicted value
	 * 
	 * @see gdsc.smlm.function.NonLinearFunction#eval(int, double[])
	 */
	public double eval(final int x, final double[] dyda)
	{
		// First parameter is the background level 
		dyda[0] = 1.0; // Gradient for a constant background is 1

		// Unpack the predictor into the dimensions
		final int x1 = x / maxx;
		final int x0 = x % maxx;

		return background + gaussian(x0, x1, dyda);
	}

	private double gaussian(final int x0, final int x1, final double[] dy_da)
	{
		final double dx = x0 - x0pos;
		final double dy = x1 - x1pos;
		final double dx2 = dx * dx;
		final double dxy = dx * dy;
		final double dy2 = dy * dy;

		// Calculate gradients

		final double exp = FastMath.exp(aa * dx2 + bb * dxy + cc * dy2);
		dy_da[1] = n * exp;
		final double y = height * exp;
		dy_da[2] = y * (aa2 * dx2 + bb2 * dxy + cc2 * dy2);

		dy_da[3] = y * (-2.0 * aa * dx - bb * dy);
		dy_da[4] = y * (-2.0 * cc * dy - bb * dx);

		dy_da[5] = y * (nx + ax * dx2 + bx * dxy + cx * dy2);
		dy_da[6] = y * (ny + ay * dx2 + by * dxy + cy * dy2);

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

		return background + height * FastMath.exp(aa * dx * dx + bb * dx * dy + cc * dy * dy);
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
	public boolean evaluatesAngle()
	{
		return true;
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
		return 6;
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
