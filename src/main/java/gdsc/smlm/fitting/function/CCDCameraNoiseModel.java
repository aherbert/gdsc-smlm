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
 *---------------------------------------------------------------------------*/

/**
 * Defines the expected variance of a signal recorded on a CCD or EM-CCD Camera. The model assumes a Gaussian read
 * noise, photon shot noise and an EM-gain noise factor.
 */
public class CCDCameraNoiseModel implements NoiseModel
{
	private float bias = 0;
	private float readNoise2 = 0;
	private boolean emCCD = true;

	public CCDCameraNoiseModel(final float readNoise, final boolean emCCD)
	{
		setReadNoise(readNoise);
		setEmCCD(emCCD);
	}

	public CCDCameraNoiseModel(final float readNoise, final float bias, final boolean emCCD)
	{
		setReadNoise(readNoise);
		setBias(bias);
		setEmCCD(emCCD);
	}

	/**
	 * Compute the expected variance of the signal from a CCD camera.
	 * <p>
	 * The noise model assumes that the camera may have a bias offset. The signal is computed as the input value minus
	 * the bias. The variance is computed using:
	 * 
	 * <pre>
	 * variance = read_noise^2 + shot_noise^2 x em-ccd noise factor
	 *          = read_noise^2 + signal x (emCCD) ? 2 : 1
	 * </pre>
	 * 
	 * The read noise is Gaussian read noise of the CCD camera.
	 * <p>
	 * The shot noise is Poisson noise of the signal. Since the variance of the Poisson distribution is the mean so we
	 * can use the signal directly.
	 * <p>
	 * The em-ccd noise factor is sqrt(2) for EM CCD cameras, otherwise it is 1. This is only applied to the signal
	 * noise standard deviation. Applying it directly to the signal variance uses a factor of sqrt(2)^2 = 2.
	 * 
	 * @see gdsc.smlm.fitting.function.NoiseModel#variance(float)
	 */
	public float variance(final float value)
	{
		return readNoise2 + Math.max(value - bias, 0f) * ((emCCD) ? 2f : 1f);
	}

	/**
	 * @return the bias
	 */
	public float getBias()
	{
		return bias;
	}

	/**
	 * @param bias
	 *            the bias to set
	 */
	public void setBias(float bias)
	{
		this.bias = bias;
	}

	/**
	 * @return the read noise
	 */
	public float getReadNoise()
	{
		return (float) Math.sqrt(readNoise2);
	}

	/**
	 * @param readNoise
	 *            the read noise to set
	 */
	public void setReadNoise(float readNoise)
	{
		readNoise2 = readNoise * readNoise;
	}

	/**
	 * @return true if the camera is an Electron Multiplying CCD camera
	 */
	public boolean isEmCCD()
	{
		return emCCD;
	}

	/**
	 * @param emCCD
	 *            Set to true if the camera is an Electron Multiplying CCD camera
	 */
	public void setEmCCD(boolean emCCD)
	{
		this.emCCD = emCCD;
	}
}
