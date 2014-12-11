package gdsc.smlm.function.gaussian;

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
 * Abstract base class for an 2-dimensional Gaussian function for a configured number of peaks.
 * <p>
 * The function will calculate the value of the Gaussian and evaluate the gradient of a set of parameters. The class can
 * specify which of the following parameters the function will evaluate:<br/>
 * background, signal, angle, position0, position1, sd0, sd1
 * <p>
 * The class provides an index of the position in the parameter array where the parameter is expected.
 */
public abstract class Gaussian2DFunction extends GaussianFunction
{
	public static final double ONE_OVER_TWO_PI = 0.5 / Math.PI; 

	public static final int BACKGROUND = 0;
	public static final int SIGNAL = 1;
	public static final int ANGLE = 2;
	public static final int X_POSITION = 3;
	public static final int Y_POSITION = 4;
	public static final int X_SD = 5;
	public static final int Y_SD = 6;
	
	protected int maxx;
	
	public Gaussian2DFunction(int maxx)
	{
		setMaxX(maxx);
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
	 * @maxx the maximum size in the first dimension. Default to 1 if not positive.
	 */
	public void setMaxX(int maxx)
	{
		if (maxx < 1)
			maxx = 1;
		this.maxx = maxx;
	}

	/**
	 * @return the maximum size in the first dimension
	 */
	public int getMaxX()
	{
		return maxx;
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
		// Background + n * { Signal, Angle, Xpos, Ypos, Xsd, Ysd }
		int nparams = (gf.evaluatesBackground() ? 1 : 0) + nPeaks * gf.getParametersPerPeak();
		int[] indices = new int[nparams];

		int p = 0;
		if (gf.evaluatesBackground())
			indices[p++] = 0;
		for (int n = 0, i = 0; n < nPeaks; n++, i += 6)
		{
			if (gf.evaluatesSignal())
				indices[p++] = i + SIGNAL;
			if (gf.evaluatesAngle())
				indices[p++] = i + ANGLE;
			// All functions evaluate the position gradient
			indices[p++] = i + X_POSITION;
			indices[p++] = i + Y_POSITION;
			if (gf.evaluatesSD0())
				indices[p++] = i + X_SD;
			if (gf.evaluatesSD1())
				indices[p++] = i + Y_SD;
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
