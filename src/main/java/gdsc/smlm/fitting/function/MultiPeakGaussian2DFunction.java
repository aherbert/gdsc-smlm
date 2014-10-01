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
	protected int npeaks;
	protected int[] gradientIndices;
	
	/**
	 * @param npeaks The number of peaks
	 */
	public MultiPeakGaussian2DFunction(int npeaks, int maxx)
	{
		super(maxx);
		this.npeaks = npeaks;
		this.gradientIndices = createGradientIndices(npeaks);
	}
	
	/* (non-Javadoc)
	 * @see gdsc.fitting.function.GaussianFunction#getNPeaks()
	 */
	@Override
	public int getNPeaks()
	{
		return npeaks;
	}
	
	/* (non-Javadoc)
	 * @see gdsc.smlm.fitting.function.NonLinearFunction#gradientIndices()
	 */
	public int[] gradientIndices()
	{
		return gradientIndices;
	}
}
