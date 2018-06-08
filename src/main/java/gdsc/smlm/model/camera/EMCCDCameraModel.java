package gdsc.smlm.model.camera;

import java.awt.Rectangle;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * An EM-CCD camera model with all pixels treated equally.
 *
 * @author Alex Herbert
 */
public class EMCCDCameraModel extends FixedPixelCameraModel
{
	/**
	 * Instantiates a new EM-CCD camera model.
	 *
	 * @param bias
	 *            the bias (in counts)
	 * @param gain
	 *            the total gain (count/photon)
	 * @param emGain
	 *            the EM-gain
	 */
	public EMCCDCameraModel(float bias, float gain)
	{
		this(bias, gain, 0f);
	}

	/**
	 * Instantiates a new EM-CCD camera model.
	 *
	 * @param bias
	 *            the bias (in counts)
	 * @param gain
	 *            the total gain (count/photon)
	 * @param emGain
	 *            the EM-gain
	 */
	public EMCCDCameraModel(double bias, double gain)
	{
		this(bias, gain, 0d);
	}

	/**
	 * Instantiates a new EM-CCD camera model.
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
	public EMCCDCameraModel(float bias, float gain, float variance)
	{
		super(bias, gain, variance);
	}

	/**
	 * Instantiates a new EM-CCD camera model.
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
	public EMCCDCameraModel(double bias, double gain, double variance)
	{
		super(bias, gain, variance);
	}

	/**
	 * {@inheritDoc}
	 * <p>
	 * Note: This is an EM-CCD camera model. The normalised variance represents the effective read noise in incident
	 * photons (i.e. before EM-gain). This can be combined with the expected shot variance of a Poisson distribution
	 * (mean) scaled by the EM-amplification noise factor (2) to obtain the total variance in photon units:
	 * 
	 * <pre>
	 * Total variance (photons) = [Poisson mean] * 2 + [normalised variance]
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
	 * Note: This is an EM-CCD camera model. The normalised variance represents the effective read noise in incident
	 * photons (i.e. before EM-gain). This can be combined with the expected shot variance of a Poisson distribution
	 * (mean) scaled by the EM-amplification noise factor (2) to obtain the total variance in photon units:
	 * 
	 * <pre>
	 * Total variance (photons) = [Poisson mean] * 2 + [normalised variance]
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
	public EMCCDCameraModel copy()
	{
		return clone();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	protected EMCCDCameraModel clone()
	{
		return (EMCCDCameraModel) super.clone();
	}
}
