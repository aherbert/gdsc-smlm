package gdsc.smlm.function.gaussian;

import gdsc.smlm.function.NoiseModel;
import gdsc.smlm.function.NonLinearFunction;

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
 * background, signal, angle[N-1], position[N], sd[N]
 * <p>
 * The class provides an index of the position in the parameter array where the parameter is expected.
 */
public abstract class GaussianFunction implements NonLinearFunction
{
	/**
	 * The factor for converting a Gaussian standard deviation to Full Width at Half Maxima (FWHM)
	 */
	public static final double SD_TO_FWHM_FACTOR = (2.0 * Math.sqrt(2.0 * Math.log(2.0)));

	/**
	 * The factor for converting a Gaussian standard deviation to Half Width at Half Maxima (FWHM)
	 */
	public static final double SD_TO_HWHM_FACTOR = (Math.sqrt(2.0 * Math.log(2.0)));
	
	private NoiseModel noiseModel = null;

	/**
	 * @return the number of dimensions
	 */
	public abstract int getNDimensions();

	/**
	 * @return the dimensions
	 */
	public abstract int[] getDimensions();

	/**
	 * @return the number of peaks
	 */
	public abstract int getNPeaks();

	/**
	 * @return True if the function can evaluate the background gradient
	 */
	public abstract boolean evaluatesBackground();

	/**
	 * @return True if the function can evaluate the signal gradient
	 */
	public abstract boolean evaluatesSignal();

	/**
	 * @return True if the function can evaluate the angle gradient
	 */
	public abstract boolean evaluatesAngle();

	/**
	 * @return True if the function can evaluate the position gradient
	 */
	public abstract boolean evaluatesPosition();

	/**
	 * @return True if the function can evaluate the standard deviation gradient for the 1st dimension
	 */
	public abstract boolean evaluatesSD0();

	/**
	 * @return True if the function can evaluate the standard deviation gradient for the 2nd dimension
	 */
	public abstract boolean evaluatesSD1();

	/**
	 * @return The number of parameters per peak
	 */
	public abstract int getParametersPerPeak();

	/**
	 * Execute the {@link #eval(int, float[])} method and set the expected variance using the noise model
	 * 
	 * @throws NullPointerException
	 *             if the noise model is null
	 * @see gdsc.smlm.function.NonLinearFunction#eval(int, float[], float[])
	 */
	public double eval(final int x, final double[] dyda, final double[] w) throws NullPointerException
	{
		final double value = eval(x, dyda);
		//w[0] = (noiseModel == null) ? 1 : noiseModel.variance(value);
		// Just throw a null pointer exception if noiseModel is null
		w[0] = noiseModel.variance(value);
		return value;
	}

	/**
	 * @return the noise model
	 */
	public NoiseModel getNoiseModel()
	{
		return noiseModel;
	}

	/**
	 * Set the noise model used in {@link #eval(int, float[], float[])}.
	 * 
	 * @param noiseModel
	 *            the noise model to set
	 */
	public void setNoiseModel(NoiseModel noiseModel)
	{
		this.noiseModel = noiseModel;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.function.NonLinearFunction#canComputeWeights()
	 */
	public boolean canComputeWeights()
	{
		return (noiseModel != null);
	}
}
