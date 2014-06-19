package gdsc.smlm.fitting.function.gaussian;


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
public class NBEllipticalGaussian2DFunction extends EllipticalGaussian2DFunction
{
	/**
	 * Constructor
	 * 
	 * @param npeaks
	 *            The number of peaks
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 */
	public NBEllipticalGaussian2DFunction(int npeaks, int maxx)
	{
		super(npeaks, maxx);
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
		
		// Unpack the predictor into the dimensions
		;
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
	
	@Override
	public boolean evaluatesBackground()
	{
		return false;
	}
}
