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
 * Evaluates an 2-dimensional Gaussian function for a configured number of peaks.
 * <p>
 * The single parameter x in the {@link #eval(int, double[])} function is assumed to be a linear index into
 * 2-dimensional
 * data. The dimensions of the data must be specified to allow unpacking to coordinates.
 * <p>
 * Data should be packed in descending dimension order, e.g. Y,X : Index for [x,y] = MaxX*y + x.
 */
public class FreeCircularGaussian2DFunction extends MultiPeakGaussian2DFunction
{
	protected static final int PARAMETERS_PER_PEAK = 5;

	protected boolean[] zeroAngle;
	protected final double[][] peakFactors;
	protected double[] a;

	/**
	 * Constructor
	 * 
	 * @param npeaks
	 *            The number of peaks
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 * @param maxy
	 *            The maximum y value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 */
	public FreeCircularGaussian2DFunction(int npeaks, int maxx, int maxy)
	{
		super(npeaks, maxx, maxy);
		zeroAngle = new boolean[npeaks];
		peakFactors = new double[npeaks][13];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.gaussian.Gaussian2DFunction#copy()
	 */
	@Override
	public Gaussian2DFunction copy()
	{
		return new FreeCircularGaussian2DFunction(npeaks, maxx, maxy);
	}

	protected static final int N = 0;
	protected static final int HEIGHT = 1;
	protected static final int AA = 2;
	protected static final int BB = 3;
	protected static final int CC = 4;
	protected static final int NX = 5;
	protected static final int AX = 6;
	protected static final int BX = 7;
	protected static final int CX = 8;
	protected static final int NY = 9;
	protected static final int AY = 10;
	protected static final int BY = 11;
	protected static final int CY = 12;

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.fitting.function.NonLinearFunction#initialise(double[])
	 */
	public void initialise(double[] a)
	{
		this.a = a;
		// Precalculate multiplication factors
		for (int j = 0; j < npeaks; j++)
		{
			final double theta = a[j * 6 + SHAPE];
			final double sx = a[j * 6 + X_SD];
			final double sy = a[j * 6 + Y_SD];
			final double sx2 = sx * sx;
			final double sy2 = sy * sy;
			final double sx3 = sx2 * sx;
			final double sy3 = sy2 * sy;
			peakFactors[j][N] = ONE_OVER_TWO_PI / (sx * sy);
			peakFactors[j][HEIGHT] = a[j * 6 + SIGNAL] * peakFactors[j][N];

			// All prefactors are negated since the Gaussian uses the exponential to the negative:
			// (A/2*pi*sx*sy) * exp( -( a(x-x0)^2 + 2b(x-x0)(y-y0) + c(y-y0)^2 ) )

			if (theta == 0)
			{
				zeroAngle[j] = true;

				// cosSqt = 1
				// sinSqt = 0
				// sin2t = 0

				peakFactors[j][AA] = -0.5 / sx2;
				peakFactors[j][BB] = 0;
				peakFactors[j][CC] = -0.5 / sy2;

				// For the x-width gradient
				peakFactors[j][NX] = -1.0 / sx;
				peakFactors[j][AX] = 1.0 / sx3;
				peakFactors[j][BX] = 0;
				peakFactors[j][CX] = 0;

				// For the y-width gradient
				peakFactors[j][NY] = -1.0 / sy;
				peakFactors[j][AY] = 0;
				peakFactors[j][BY] = 0;
				peakFactors[j][CY] = 1.0 / sy3;
			}
			else
			{
				zeroAngle[j] = false;

				final double cosSqt = Math.cos(theta) * Math.cos(theta);
				final double sinSqt = Math.sin(theta) * Math.sin(theta);
				final double sin2t = Math.sin(2.0 * theta);

				peakFactors[j][AA] = -0.5 * (cosSqt / sx2 + sinSqt / sy2);
				peakFactors[j][BB] = -0.25 * (-sin2t / sx2 + sin2t / sy2);
				peakFactors[j][CC] = -0.5 * (sinSqt / sx2 + cosSqt / sy2);

				// For the x-width gradient
				peakFactors[j][NX] = -1.0 / sx;
				peakFactors[j][AX] = cosSqt / sx3;
				peakFactors[j][BX] = -0.5 * sin2t / sx3;
				peakFactors[j][CX] = sinSqt / sx3;

				// For the y-width gradient
				peakFactors[j][NY] = -1.0 / sy;
				peakFactors[j][AY] = sinSqt / sy3;
				peakFactors[j][BY] = 0.5 * sin2t / sy3;
				peakFactors[j][CY] = cosSqt / sy3;
			}
		}
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
	 * the standard deviation in each of the N-dimensions. This produces an additional 1+2N coefficients per peak.
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
		// Track the position of the parameters
		int apos = 0;
		int dydapos = 0;

		// First parameter is the background level 
		double y_fit = a[BACKGROUND];
		dyda[dydapos++] = 1.0; // Gradient for a constant background is 1

		// Unpack the predictor into the dimensions
		final int x1 = x / maxx;
		final int x0 = x % maxx;

		for (int j = 0; j < npeaks; j++)
		{
			y_fit += gaussian(x0, x1, dyda, apos, dydapos, zeroAngle[j], peakFactors[j]);
			apos += 6;
			dydapos += PARAMETERS_PER_PEAK;
		}

		return y_fit;
	}

	protected double gaussian(final int x0, final int x1, final double[] dy_da, final int apos, final int dydapos,
			boolean zeroAngle, final double[] factors)
	{
		final double dx = x0 - a[apos + X_POSITION];
		final double dy = x1 - a[apos + Y_POSITION];
		final double dx2 = dx * dx;
		final double dxy = dx * dy;
		final double dy2 = dy * dy;

		// Calculate gradients

		final double aa = factors[AA];
		final double cc = factors[CC];

		if (zeroAngle)
		{
			final double exp = FastMath.exp(aa * dx2 + cc * dy2);
			dy_da[dydapos] = factors[N] * exp;
			final double y = factors[HEIGHT] * exp;

			dy_da[dydapos + 1] = y * (-2.0 * aa * dx);
			dy_da[dydapos + 2] = y * (-2.0 * cc * dy);

			dy_da[dydapos + 3] = y * (factors[NX] + factors[AX] * dx2);
			dy_da[dydapos + 4] = y * (factors[NY] + factors[CY] * dy2);

			return y;
		}
		else
		{
			final double bb = factors[BB];

			final double exp = FastMath.exp(aa * dx2 + bb * dxy + cc * dy2);
			dy_da[dydapos] = factors[N] * exp;
			final double y = factors[HEIGHT] * exp;

			dy_da[dydapos + 1] = y * (-2.0 * aa * dx - bb * dy);
			dy_da[dydapos + 2] = y * (-2.0 * cc * dy - bb * dx);

			dy_da[dydapos + 3] = y * (factors[NX] + factors[AX] * dx2 + factors[BX] * dxy + factors[CX] * dy2);
			dy_da[dydapos + 4] = y * (factors[NY] + factors[AY] * dx2 + factors[BY] * dxy + factors[CY] * dy2);

			return y;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.fitting.function.NonLinearFunction#eval(int)
	 */
	public double eval(final int x)
	{
		// Track the position of the parameters
		int apos = 0;

		// First parameter is the background level 
		double y_fit = a[BACKGROUND];

		// Unpack the predictor into the dimensions
		final int x1 = x / maxx;
		final int x0 = x % maxx;

		for (int j = 0; j < npeaks; j++, apos += 6)
		{
			y_fit += gaussian(x0, x1, apos, zeroAngle[j], peakFactors[j]);
		}

		return y_fit;
	}

	protected double gaussian(final int x0, final int x1, final int apos, boolean zeroAngle, final double[] factors)
	{
		final double dx = x0 - a[apos + X_POSITION];
		final double dy = x1 - a[apos + Y_POSITION];

		if (zeroAngle)
			return factors[HEIGHT] * FastMath.exp(factors[AA] * dx * dx + factors[CC] * dy * dy);
		else
			return factors[HEIGHT] *
					FastMath.exp(factors[AA] * dx * dx + factors[BB] * dx * dy + factors[CC] * dy * dy);
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
		return PARAMETERS_PER_PEAK;
	}
}
