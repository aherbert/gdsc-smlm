/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package gdsc.smlm.function;

/**
 * Defines the expected variance of a signal recorded on a CCD or EM-CCD Camera. The model assumes a Gaussian read
 * noise, photon shot noise and an EM-gain noise factor.
 */
public abstract class CameraNoiseModel implements NoiseModel
{
	/** The bias. */
	protected double bias = 0;
	
	/** The read noise squared (read variance). */
	protected double readNoise2 = 0;

	/**
	 * Instantiates a new camera noise model.
	 *
	 * @param readNoise
	 *            the read noise
	 */
	protected CameraNoiseModel(final double readNoise)
	{
		setReadNoise(readNoise);
	}

	/**
	 * Instantiates a new camera noise model.
	 *
	 * @param readNoise
	 *            the read noise
	 * @param bias
	 *            the bias
	 */
	protected CameraNoiseModel(final double readNoise, final double bias)
	{
		setReadNoise(readNoise);
		setBias(bias);
	}

	/**
	 * Factory method for creating camera noise models from the sub-classes.
	 *
	 * @param readNoise
	 *            the read noise
	 * @param bias
	 *            the bias
	 * @param emCCD
	 *            the EM-CCD flag
	 * @return the camera noise model
	 */
	public static CameraNoiseModel createNoiseModel(final double readNoise, final double bias, final boolean emCCD)
	{
		return (emCCD) ? new EMCCDCameraNoiseModel(readNoise, bias) : new CCDCameraNoiseModel(readNoise, bias);
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
	 * @see gdsc.smlm.function.NoiseModel#variance(double)
	 */
	@Override
	public abstract double variance(final double value);

	/**
	 * @return the bias
	 */
	public double getBias()
	{
		return bias;
	}

	/**
	 * @param bias
	 *            the bias to set
	 */
	public void setBias(double bias)
	{
		this.bias = bias;
	}

	/**
	 * @return the read noise
	 */
	public double getReadNoise()
	{
		return Math.sqrt(readNoise2);
	}

	/**
	 * @param readNoise
	 *            the read noise to set
	 */
	public void setReadNoise(double readNoise)
	{
		readNoise2 = readNoise * readNoise;
	}

	/**
	 * @return true if the camera is an Electron Multiplying CCD camera
	 */
	public abstract boolean isEmCCD();
}
