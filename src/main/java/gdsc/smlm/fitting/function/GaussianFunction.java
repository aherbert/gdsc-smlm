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
 * background, amplitude, angle[N-1], position[N], width[N]
 * <p>
 * The class provides an index of the position in the parameter array where the parameter is expected.
 */
public abstract class GaussianFunction implements NonLinearFunction
{
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
	 * @return True if the function can evaluate the amplitude gradient
	 */
	public abstract boolean evaluatesAmplitude();

	/**
	 * @return True if the function can evaluate the angle gradient
	 */
	public abstract boolean evaluatesAngle();

	/**
	 * @return True if the function can evaluate the position gradient
	 */
	public abstract boolean evaluatesPosition();

	/**
	 * @return True if the function can evaluate the width gradient for the 1st dimension
	 */
	public abstract boolean evaluatesWidth0();

	/**
	 * @return True if the function can evaluate the width gradient for the 2nd dimension
	 */
	public abstract boolean evaluatesWidth1();

	/**
	 * @return The number of parameters per peak
	 */
	public abstract int getParametersPerPeak();

	/**
	 * Execute the {@link #eval(int, float[])} method and set the expected variance using the noise model
	 * 
	 * @throws NullPointerException
	 *             if the noise model is null
	 * @see gdsc.smlm.fitting.function.NonLinearFunction#eval(int, float[], float[])
	 */
	public float eval(final int x, final float[] dyda, final float[] w) throws NullPointerException
	{
		final float value = eval(x, dyda);
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
