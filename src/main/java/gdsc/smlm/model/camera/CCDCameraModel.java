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
package gdsc.smlm.model.camera;

import java.awt.Rectangle;

/**
 * An CCD camera model with all pixels treated equally.
 *
 * @author Alex Herbert
 */
public class CCDCameraModel extends FixedPixelCameraModel
{
	/**
	 * Instantiates a new CCD camera model.
	 *
	 * @param bias
	 *            the bias (in counts)
	 * @param gain
	 *            the total gain (count/photon)
	 * @param emGain
	 *            the EM-gain
	 */
	public CCDCameraModel(float bias, float gain)
	{
		this(bias, gain, 0f);
	}

	/**
	 * Instantiates a new CCD camera model.
	 *
	 * @param bias
	 *            the bias (in counts)
	 * @param gain
	 *            the total gain (count/photon)
	 * @param emGain
	 *            the EM-gain
	 */
	public CCDCameraModel(double bias, double gain)
	{
		this(bias, gain, 0d);
	}

	/**
	 * Instantiates a new CCD camera model.
	 *
	 * @param bias
	 *            the bias (in counts)
	 * @param gain
	 *            the total gain (count/photon)
	 * @param emGain
	 *            the EM-gain
	 * @param variance
	 *            the variance (in counts)
	 */
	public CCDCameraModel(float bias, float gain, float variance)
	{
		super(bias, gain, variance);
	}

	/**
	 * Instantiates a new CCD camera model.
	 *
	 * @param bias
	 *            the bias (in counts)
	 * @param gain
	 *            the total gain (count/photon)
	 * @param emGain
	 *            the EM-gain
	 * @param variance
	 *            the variance (in counts)
	 */
	public CCDCameraModel(double bias, double gain, double variance)
	{
		super(bias, gain, variance);
	}

	/**
	 * {@inheritDoc}
	 * <p>
	 * Note: This is an CCD camera model. The normalised variance represents the effective read noise in incident
	 * photons (i.e. before gain). This can be combined with the expected shot variance of a Poisson distribution
	 * (mean) to obtain the total variance in photon units:
	 * 
	 * <pre>
	 * Total variance (photons) = [Poisson mean] + [normalised variance]
	 * </pre>
	 * 
	 * This value multiplied by the [gain]^2 is the variance in counts.
	 * 
	 * @see gdsc.smlm.model.camera.FixedPixelCameraModel#getNormalisedVariance(java.awt.Rectangle)
	 */
	@Override
	public float[] getNormalisedVariance(Rectangle bounds)
	{
		return super.getNormalisedVariance(bounds);
	}

	/**
	 * {@inheritDoc}
	 * <p>
	 * Note: This is an CCD camera model. The normalised variance represents the effective read noise in incident
	 * photons (i.e. before gain). This can be combined with the expected shot variance of a Poisson distribution
	 * (mean) to obtain the total variance in photon units:
	 * 
	 * <pre>
	 * Total variance (photons) = [Poisson mean] + [normalised variance]
	 * </pre>
	 * 
	 * This value multiplied by the [gain]^2 is the variance in counts.
	 * 
	 * @see gdsc.smlm.model.camera.FixedPixelCameraModel#getMeanNormalisedVariance(java.awt.Rectangle)
	 */
	@Override
	public double getMeanNormalisedVariance(Rectangle bounds)
	{
		return super.getMeanNormalisedVariance(bounds);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#copy()
	 */
	@Override
	public CCDCameraModel copy()
	{
		return clone();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	protected CCDCameraModel clone()
	{
		return (CCDCameraModel) super.clone();
	}
}
