package gdsc.smlm.fitting.function.gaussian;

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
 * 
 * This is an adaption of the C-code contained in the CcpNmr Analysis Program:
 * CCPN website (http://www.ccpn.ac.uk/)
 *---------------------------------------------------------------------------*/

/**
 * Evaluates an 2-dimensional Gaussian function for a single peak.
 * <p>
 * The single parameter x in the {@link #eval(int, float[])} function is assumed to be a linear index into 2-dimensional
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

	protected float background;
	protected float amplitude;
	protected float x0pos;
	protected float x1pos;

	protected float aa;
	protected float ax;

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

	protected static final double MinusFourLog2 = (double) (-4.0 * Math.log(2.0));

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
		final double sigma_x = a[X_WIDTH];

		aa = (float) (MinusFourLog2 / (sigma_x * sigma_x));

		// For the x-width gradient
		ax = (float) (-2.0 * MinusFourLog2 / (sigma_x * sigma_x * sigma_x));
	}

	/**
	 * Produce an output predicted value for a given set of input
	 * predictors (x) and coefficients (a).
	 * <p>
	 * Evaluates an 2-dimensional Gaussian function for a single peak.
	 * <p>
	 * The first coefficient is the Gaussian background level (B). The coefficients are then packed for each peak:
	 * Amplitude; position[N]; width[N]. Amplitude (A) is the height of the Gaussian. Position (p) is the position of
	 * the Gaussian in each of the N-dimensions. Width (w) is the peak width at half-max in each of the N-dimensions.
	 * This produces an additional 1+2N coefficents per peak.
	 * <p>
	 * The equation per peak is:<br/>
	 * y_peak = A * product(N) [ exp(-4 * log(2) * (x[i]-p[i]) * (x[i]-p[i]) / (w[i] * w[i])) ]
	 * <p>
	 * Note the use of the width as the peak width at half max. This can be converted to the Gaussian peak standard
	 * deviation using: SD = pwhm / 2*sqrt(2*log(2)); SD = pwhm / 2.35482
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
		;
		final int x1 = x / maxx;
		final int x0 = x % maxx;

		return background + gaussian(x0, x1, dyda);
	}

	private float gaussian(final int x0, final int x1, final float[] dy_da)
	{
		final float h = amplitude;

		final float dx = x0 - x0pos;
		final float dy = x1 - x1pos;

		final float y = (float) (h * Math.exp(aa * dx * dx + aa * dy * dy));

		// Calculate gradients
		dy_da[1] = y / h;

		dy_da[2] = y * (-2.0f * aa * dx);
		dy_da[3] = y * (-2.0f * aa * dy);

		dy_da[4] = y * (ax * dx * dx + ax * dy * dy);

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

		return background + amplitude * (float) Math.exp(aa * (dx * dx + dy * dy));
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
	public boolean evaluatesWidth0()
	{
		return true;
	}

	@Override
	public boolean evaluatesWidth1()
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
