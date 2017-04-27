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
 * Abstract base class for an N-dimensional Gaussian function for a configured number of peaks.
 * <p>
 * The function will calculate the value of the Gaussian and evaluate the gradient of a set of parameters. The class can
 * specify which of the following parameters the function will evaluate:<br/>
 * background, amplitude, angle[N-1], position[N], sd[N]
 * <p>
 * The class provides the number of peaks and the gradient indices.
 */
public abstract class MultiPeakGaussian2DFunction extends Gaussian2DFunction
{
	protected final int npeaks;
	protected final int[] gradientIndices;

	/**
	 * Instantiates a new multi peak gaussian 2D function.
	 *
	 * @param npeaks
	 *            The number of peaks
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 * @param maxy
	 *            The maximum y value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 */
	public MultiPeakGaussian2DFunction(int npeaks, int maxx, int maxy)
	{
		super(maxx, maxy);
		this.npeaks = npeaks;
		this.gradientIndices = createGradientIndices(npeaks);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.fitting.function.GaussianFunction#getNPeaks()
	 */
	@Override
	public int getNPeaks()
	{
		return npeaks;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.function.NonLinearFunction#gradientIndices()
	 */
	public int[] gradientIndices()
	{
		return gradientIndices;
	}
}
