package gdsc.smlm.fitting.function;

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
 * Abstract base class for an 2-dimensional Gaussian function for a configured number of peaks.
 * <p>
 * The function will calculate the value of the Gaussian and evaluate the gradient of a set of parameters. The class can
 * specify which of the following parameters the function will evaluate:<br/>
 * background, amplitude, angle, position0, position1, width0, width1
 * <p>
 * The class provides an index of the position in the parameter array where the parameter is expected.
 */
public abstract class Gaussian2DFunction extends GaussianFunction
{
	public static final int BACKGROUND = 0;
	public static final int AMPLITUDE = 1;
	public static final int ANGLE = 2;
	public static final int X_POSITION = 3;
	public static final int Y_POSITION = 4;
	public static final int X_WIDTH = 5;
	public static final int Y_WIDTH = 6;
	
	protected int maxx;
	
	public Gaussian2DFunction(int maxx)
	{
		this.maxx = maxx;
	}
	
	/**
	 * @return the number of dimensions
	 */
	public int getNDimensions()
	{
		return 2;
	}

	/**
	 * @return the dimensions
	 */
	public int[] getDimensions()
	{
		return new int[] { maxx, 0 };
	}

	/**
	 * Build the index array that maps the gradient index back to the original parameter index so that:<br/>
	 * a[indices[i]] += dy_da[i]
	 * 
	 * @param nPeaks
	 * @return The indices
	 */
	protected int[] createGradientIndices(int nPeaks)
	{
		return createGradientIndices(nPeaks, this);
	}

	protected static int[] createGradientIndices(int nPeaks, GaussianFunction gf)
	{
		// Parameters are: 
		// Background + n * { Amplitude, Angle, Xpos, Ypos, Xwidth, yWidth }
		int nparams = (gf.evaluatesBackground() ? 1 : 0) + nPeaks * gf.getParametersPerPeak();
		int[] indices = new int[nparams];

		int p = 0;
		if (gf.evaluatesBackground())
			indices[p++] = 0;
		for (int n = 0, i = 0; n < nPeaks; n++, i += 6)
		{
			indices[p++] = i + AMPLITUDE;
			if (gf.evaluatesAngle())
				indices[p++] = i + ANGLE;
			indices[p++] = i + X_POSITION;
			indices[p++] = i + Y_POSITION;
			if (gf.evaluatesWidth0())
				indices[p++] = i + X_WIDTH;
			if (gf.evaluatesWidth1())
				indices[p++] = i + Y_WIDTH;
		}

		return indices;
	}
	
	/**
	 * Locate the index within the gradient indices for the specified parameter
	 * @param parameterIndex
	 * @return the gradient index (or -1 if not present)
	 */
	public int findGradientIndex(int parameterIndex)
	{
		int[] gradientIndices = gradientIndices();
		for (int i=0; i<gradientIndices.length; i++)
			if (gradientIndices[i] == parameterIndex)
				return i;
		return -1;
	}
}
