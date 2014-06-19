package gdsc.smlm.fitting.function.gaussian;

import gdsc.smlm.fitting.function.MultiPeakGaussian2DFunction;

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
 * Evaluates an 2-dimensional elliptical Gaussian function for a configured number of peaks.
 * <p>
 * The single parameter x in the {@link #eval(int, float[])} function is assumed to be a linear index into 2-dimensional
 * data. The dimensions of the data must be specified to allow unpacking to coordinates.
 * <p>
 * Data should be packed in descending dimension order, e.g. Y,X : Index for [x,y] = MaxX*y + x.
 */
public class EllipticalGaussian2DFunction extends MultiPeakGaussian2DFunction
{
	protected static final int PARAMETERS_PER_PEAK = 6;

	protected float[][] peakFactors;
	protected float[] a;

	/**
	 * Constructor
	 * 
	 * @param npeaks
	 *            The number of peaks
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 */
	public EllipticalGaussian2DFunction(int npeaks, int maxx)
	{
		super(npeaks, maxx);
	}

	protected static final double MinusFourLog2 = (double) (-4.0 * Math.log(2.0));
	protected static final int AA = 0;
	protected static final int BB = 1;
	protected static final int CC = 2;
	protected static final int AA2 = 3;
	protected static final int BB2 = 4;
	protected static final int CC2 = 5;
	protected static final int AX = 6;
	protected static final int BX = 7;
	protected static final int CX = 8;
	protected static final int AY = 9;
	protected static final int BY = 10;
	protected static final int CY = 11;

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.fitting.function.NonLinearFunction#initialise(float[])
	 */
	public void initialise(float[] a)
	{
		this.a = a;
		// Precalculate multiplication factors
		peakFactors = new float[npeaks][12];
		for (int j = 0; j < npeaks; j++)
		{
			final double theta = a[j * 6 + ANGLE];
			final double sigma_x = a[j * 6 + X_WIDTH];
			final double sigma_y = a[j * 6 + Y_WIDTH];

			peakFactors[j][AA] = (float) (MinusFourLog2 * Math.cos(theta) * Math.cos(theta) / (sigma_x * sigma_x) + MinusFourLog2 *
					Math.sin(theta) * Math.sin(theta) / (sigma_y * sigma_y));
			peakFactors[j][BB] = (float) (-MinusFourLog2 * Math.sin(2 * theta) / (sigma_x * sigma_x) + MinusFourLog2 *
					Math.sin(2 * theta) / (sigma_y * sigma_y));
			peakFactors[j][CC] = (float) (MinusFourLog2 * Math.sin(theta) * Math.sin(theta) / (sigma_x * sigma_x) + MinusFourLog2 *
					Math.cos(theta) * Math.cos(theta) / (sigma_y * sigma_y));

			// For the angle gradient
			peakFactors[j][AA2] = (float) (MinusFourLog2 * -2.0 * Math.sin(theta) * Math.cos(theta) /
					(sigma_x * sigma_x) + MinusFourLog2 * 2.0 * Math.sin(theta) * Math.cos(theta) / (sigma_y * sigma_y));
			peakFactors[j][BB2] = (float) (MinusFourLog2 * -2.0 * Math.cos(2 * theta) / (sigma_x * sigma_x) + MinusFourLog2 *
					2.0 * Math.cos(2 * theta) / (sigma_y * sigma_y));
			peakFactors[j][CC2] = (float) (MinusFourLog2 * 2.0 * Math.sin(theta) * Math.cos(theta) /
					(sigma_x * sigma_x) + MinusFourLog2 * -2.0 * Math.sin(theta) * Math.cos(theta) /
					(sigma_y * sigma_y));

			// For the x-width gradient
			peakFactors[j][AX] = (float) (-2.0 * MinusFourLog2 * Math.cos(theta) * Math.cos(theta) / (sigma_x * sigma_x * sigma_x));
			peakFactors[j][BX] = (float) (2.0 * MinusFourLog2 * Math.sin(2 * theta) / (sigma_x * sigma_x * sigma_x));
			peakFactors[j][CX] = (float) (-2.0 * MinusFourLog2 * Math.sin(theta) * Math.sin(theta) / (sigma_x * sigma_x * sigma_x));

			// For the y-width gradient
			peakFactors[j][AY] = (float) (-2.0 * MinusFourLog2 * Math.sin(theta) * Math.sin(theta) / (sigma_y * sigma_y * sigma_y));
			peakFactors[j][BY] = (float) (-2.0 * MinusFourLog2 * Math.sin(2 * theta) / (sigma_y * sigma_y * sigma_y));
			peakFactors[j][CY] = (float) (-2.0 * MinusFourLog2 * Math.cos(theta) * Math.cos(theta) / (sigma_y * sigma_y * sigma_y));
		}
	}

	/**
	 * Produce an output predicted value for a given set of input
	 * predictors (x) and coefficients (a).
	 * <p>
	 * Evaluates an 2-dimensional elliptical Gaussian function for a single peak.
	 * <p>
	 * The first coefficient is the Gaussian background level (B). The coefficients are then packed for each peak:
	 * Amplitude; Angle; position[N]; width[N]. Amplitude (A) is the height of the Gaussian. Angle (r) is the rotation
	 * angle of the ellipse. Position (p) is the position of the Gaussian in each of the N-dimensions. Width (w) is the
	 * peak width at half-max in each of the N-dimensions. This produces an additional 1+2N coefficients per peak.
	 * <p>
	 * The equation per peak is:<br/>
	 * y_peak = A * exp( -( a(x-x0)^2 + 2b(x-x0)(y-y0) + c(y-y0)^2 ) )<br/>
	 * Where: <br/>
	 * a = 4*ln(2)*cos(r)^2 /lwX^2 + 4*ln(2)*sin(r)^2 /lwY^2 <br/>
	 * b = -2*ln(2)*sin(2r)^2/lwX^2 + 2*ln(2)*sin(2r)^2/lwY^2 <br/>
	 * c = 4*ln(2)*sin(r)^2 /lwX^2 + 4*ln(2)*cos(r)^2 /lwY^2
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
		// Track the position of the parameters
		int apos = 0;
		int dydapos = 0;

		// First parameter is the background level 
		float y_fit = a[BACKGROUND];
		dyda[dydapos++] = 1; // Gradient for a constant background is 1

		// Unpack the predictor into the dimensions
		final int x1 = x / maxx;
		final int x0 = x % maxx;

		for (int j = 0; j < npeaks; j++)
		{
			y_fit += gaussian(x0, x1, dyda, apos, dydapos, peakFactors[j]);
			apos += 6;
			dydapos += PARAMETERS_PER_PEAK;
		}

		return y_fit;
	}

	protected float gaussian(final int x0, final int x1, final float[] dy_da, final int apos, final int dydapos,
			final float[] factors)
	{
		final float h = a[apos + AMPLITUDE];

		final float dx = x0 - a[apos + X_POSITION];
		final float dy = x1 - a[apos + Y_POSITION];

		final float aa = factors[AA];
		final float bb = factors[BB];
		final float cc = factors[CC];
		final float aa2 = factors[AA2];
		final float bb2 = factors[BB2];
		final float cc2 = factors[CC2];
		final float ax = factors[AX];
		final float bx = factors[BX];
		final float cx = factors[CX];
		final float ay = factors[AY];
		final float by = factors[BY];
		final float cy = factors[CY];

		final float y = (float) (h * Math.exp(aa * dx * dx + bb * dx * dy + cc * dy * dy));

		// Calculate gradients
		dy_da[dydapos] = y / h;
		dy_da[dydapos + 1] = y * (aa2 * dx * dx + bb2 * dx * dy + cc2 * dy * dy);

		dy_da[dydapos + 2] = y * (-2.0f * aa * dx - bb * dy);
		dy_da[dydapos + 3] = y * (-2.0f * cc * dy - bb * dx);

		dy_da[dydapos + 4] = y * (ax * dx * dx + bx * dx * dy + cx * dy * dy);
		dy_da[dydapos + 5] = y * (ay * dx * dx + by * dx * dy + cy * dy * dy);

		return y;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.fitting.function.NonLinearFunction#eval(int)
	 */
	public float eval(final int x)
	{
		// Track the position of the parameters
		int apos = 0;

		// First parameter is the background level 
		float y_fit = a[BACKGROUND];

		// Unpack the predictor into the dimensions
		final int x1 = x / maxx;
		final int x0 = x % maxx;

		for (int j = 0; j < npeaks; j++, apos += 6)
		{
			y_fit += gaussian(x0, x1, apos, peakFactors[j]);
		}

		return y_fit;
	}

	protected float gaussian(final int x0, final int x1, final int apos, final float[] factors)
	{
		final float h = a[apos + AMPLITUDE];

		final float dx = x0 - a[apos + X_POSITION];
		final float dy = x1 - a[apos + Y_POSITION];

		final float aa = factors[AA];
		final float bb = factors[BB];
		final float cc = factors[CC];

		return h * (float) Math.exp(aa * dx * dx + bb * dx * dy + cc * dy * dy);
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
		return true;
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
		return true;
	}

	@Override
	public int getParametersPerPeak()
	{
		return PARAMETERS_PER_PEAK;
	}
}
