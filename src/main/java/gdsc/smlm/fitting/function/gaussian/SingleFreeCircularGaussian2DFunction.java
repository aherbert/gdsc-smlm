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
 * The single parameter x in the {@link #eval(int, float[])} function is assumed to be a linear index into 2-dimensional
 * data. The dimensions of the data must be specified to allow unpacking to coordinates.
 * <p>
 * Data should be packed in descending dimension order, e.g. Y,X : Index for [x,y] = MaxX*y + x.
 */
public class SingleFreeCircularGaussian2DFunction extends Gaussian2DFunction
{
	private static int[] gradientIndices;
	static
	{
		gradientIndices = createGradientIndices(1, new SingleFreeCircularGaussian2DFunction(1));
	}

	protected float background;
	protected float amplitude;
	protected float x0pos;
	protected float x1pos;

	protected float aa;
	protected float bb;
	protected float cc;
	protected float ax;
	protected float bx;
	protected float cx;
	protected float ay;
	protected float by;
	protected float cy;

	/**
	 * Constructor
	 * 
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 */
	public SingleFreeCircularGaussian2DFunction(int maxx)
	{
		super(maxx);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.fitting.function.NonLinearFunction#initialise(float[])
	 */
	public void initialise(float[] a)
	{
		background = a[BACKGROUND];
		amplitude = a[AMPLITUDE];
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
		final double sin2t = Math.sin(2 * theta);

		// All prefactors are negated since the Gaussian uses the exponential to the negative:
		// A * exp( -( a(x-x0)^2 + 2b(x-x0)(y-y0) + c(y-y0)^2 ) )
		
		aa = (float) (-0.5 * (cosSqt / sx2 + sinSqt / sy2));
		bb = (float) (-0.25 * (-sin2t / sx2 + sin2t / sy2));
		cc = (float) (-0.5 * (sinSqt / sx2 + cosSqt / sy2));;

		// For the x-width gradient
		ax = (float) (cosSqt / sx3);
		bx = (float) (-0.5 * sin2t / sx3);
		cx = (float) (sinSqt / sx3);

		// For the y-width gradient
		ay = (float) (sinSqt / sy3);
		by = (float) (0.5 * sin2t / sy3);
		cy = (float) (cosSqt / sy3);
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
	 * y_peak = A * exp( -( a(x-x0)^2 + 2b(x-x0)(y-y0) + c(y-y0)^2 ) )<br/>
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
	 * @see gdsc.smlm.fitting.function.NonLinearFunction#eval(int, float[])
	 */
	public float eval(final int x, final float[] dyda)
	{
		// First parameter is the background level 
		dyda[0] = 1; // Gradient for a constant background is 1

		// Unpack the predictor into the dimensions
		final int x1 = x / maxx;
		final int x0 = x % maxx;

		return background + gaussian(x0, x1, dyda);
	}

	private float gaussian(final int x0, final int x1, final float[] dy_da)
	{
		final float h = amplitude;

		final float dx = x0 - x0pos;
		final float dy = x1 - x1pos;
		final float dx2 = dx * dx;
		final float dxy = dx * dy;
		final float dy2 = dy * dy;

		final float y = (float) (h * FastMath.exp(aa * dx2 + bb * dxy + cc * dy2));

		// Calculate gradients
		dy_da[1] = y / h;

		dy_da[2] = y * (-2.0f * aa * dx - bb * dy);
		dy_da[3] = y * (-2.0f * cc * dy - bb * dx);

		dy_da[4] = y * (ax * dx2 + bx * dxy + cx * dy2);
		dy_da[5] = y * (ay * dx2 + by * dxy + cy * dy2);

		return y;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.fitting.function.NonLinearFunction#eval(int)
	 */
	public float eval(final int x)
	{
		// Unpack the predictor into the dimensions
		final int x1 = x / maxx;
		final int x0 = x % maxx;

		final float dx = x0 - x0pos;
		final float dy = x1 - x1pos;

		return background + amplitude * (float) FastMath.exp(aa * dx * dx + bb * dx * dy + cc * dy * dy);
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
