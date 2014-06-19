package gdsc.smlm.fitting;

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
 * Contains the fitting result
 */
public class FitResult
{
	private FitStatus status;
	private int degreesOfFreedom;
	private double error;
	private float[] initialParameters;
	private float[] parameters;
	private float[] parametersDev;
	private int nPeaks;
	private int nFittedParameters;
	private Object data;

	/**
	 * Constructor
	 * 
	 * @param status
	 * @param degreesOfFreedom
	 * @param error
	 * @param initialParameters
	 * @param parameters
	 * @param parametersDev
	 * @param nPeaks
	 * @param nFittedParameters
	 * @param data
	 */
	public FitResult(FitStatus status, int degreesOfFreedom, double error, float[] initialParameters,
			float[] parameters, float[] parametersDev, int nPeaks, int nFittedParameters, Object data)
	{
		this.status = status;
		this.degreesOfFreedom = degreesOfFreedom;
		this.error = error;
		this.initialParameters = initialParameters;
		this.parameters = parameters;
		this.parametersDev = parametersDev;
		this.nPeaks = nPeaks;
		this.nFittedParameters = nFittedParameters;
		this.data = data;
	}

	/**
	 * Constructor
	 * 
	 * @param result
	 */
	public FitResult(FitStatus result)
	{
		this.status = result;
		this.degreesOfFreedom = 0;
		this.error = 0;
		this.initialParameters = null;
		this.parameters = null;
		this.parametersDev = null;
		this.nPeaks = 0;
		this.data = null;
	}

	/**
	 * @return the fit status
	 */
	public FitStatus getStatus()
	{
		return status;
	}

	/**
	 * @return the degrees of freedom for the fit
	 */
	public int getDegreesOfFreedom()
	{
		return degreesOfFreedom;
	}

	/**
	 * @return the error value for the fit (e.g. mean-squared error or the reduced Chi-squared measure)
	 */
	public double getError()
	{
		return error;
	}

	/**
	 * @return the initial parameters
	 */
	public float[] getInitialParameters()
	{
		return initialParameters;
	}

	/**
	 * Get the fitted parameters
	 * <p>
	 * The first coefficient is the Gaussian background level (B). The coefficients are then packed for each peak:
	 * Amplitude; angle[N-1]; position[N]; width[N]. Amplitude (A) is the height of the Gaussian. Angle is the rotation
	 * in the i to i+1 axis (in degrees). Position (p) is the position of the Gaussian in each of the N-dimensions.
	 * Width (w) is the peak width at half-max in each of the N-dimensions. This produces an additional 3N coefficients
	 * per peak.
	 * 
	 * @return the fitted parameters
	 */
	public float[] getParameters()
	{
		return parameters;
	}

	/**
	 * Get the standard deviations for the fitted parameters
	 * <p>
	 * The first coefficient is the Gaussian background level (B). The coefficients are then packed for each peak:
	 * Amplitude; angle[N-1]; position[N]; width[N]. Amplitude (A) is the height of the Gaussian. Angle is the rotation
	 * in the i to i+1 axis (in degrees). Position (p) is the position of the Gaussian in each of the N-dimensions.
	 * Width (w) is the peak width at half-max in each of the N-dimensions. This produces an additional 3N coefficients
	 * per peak.
	 * 
	 * @return the fitted parameters
	 */
	public float[] getParameterStdDev()
	{
		return parametersDev;
	}

	/**
	 * @return the number of fitted peaks
	 */
	public int getNumberOfPeaks()
	{
		return nPeaks;
	}

	/**
	 * @return the number of fitted parameters (e.g. fixed width fitting will reduce the number of fitted parameters)
	 */
	public int getNumberOfFittedParameters()
	{
		return nFittedParameters;
	}

	/**
	 * Returns an object containing data about the fit status. This is used to pass additional information depening on
	 * the fitting status.
	 * 
	 * @return the data
	 */
	public Object getStatusData()
	{
		return data;
	}
}