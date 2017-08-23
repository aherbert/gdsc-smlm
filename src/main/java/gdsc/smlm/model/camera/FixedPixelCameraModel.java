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
 * A camera model with all pixels treated equally.
 *
 * @author Alex Herbert
 */
public class FixedPixelCameraModel extends BaseCameraModel
{
	private final float bias, gain, variance, var_g2;

	/**
	 * Instantiates a new fixed pixel camera model.
	 *
	 * @param bias
	 *            the bias (in counts)
	 * @param gain
	 *            the gain (count/photon)
	 */
	public FixedPixelCameraModel(float bias, float gain)
	{
		this(bias, gain, 0);
	}

	/**
	 * Instantiates a new fixed pixel camera model.
	 *
	 * @param bias
	 *            the bias (in counts)
	 * @param gain
	 *            the gain (count/photon)
	 * @param variance
	 *            the variance (in counts)
	 */
	public FixedPixelCameraModel(float bias, float gain, float variance)
	{
		checkBias(bias);
		checkGain(gain);
		checkVariance(variance);
		this.bias = bias;
		this.gain = gain;
		this.variance = variance;
		this.var_g2 = variance / (gain * gain);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#getBounds()
	 */
	public Rectangle getBounds()
	{
		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#crop(java.awt.Rectangle)
	 */
	public CameraModel crop(Rectangle bounds)
	{
		return this;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#isPerPixelModel()
	 */
	public boolean isPerPixelModel()
	{
		return false;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#getBias(java.awt.Rectangle)
	 */
	public float[] getBias(Rectangle bounds)
	{
		return newArray(bounds, bias);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#getGain(java.awt.Rectangle)
	 */
	public float[] getGain(Rectangle bounds)
	{
		return newArray(bounds, gain);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#getVariance(java.awt.Rectangle)
	 */
	public float[] getVariance(Rectangle bounds)
	{
		return newArray(bounds, variance);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#getNormalisedVariance(java.awt.Rectangle)
	 */
	public float[] getNormalisedVariance(Rectangle bounds)
	{
		return newArray(bounds, var_g2);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#getWeights(java.awt.Rectangle)
	 */
	public float[] getWeights(Rectangle bounds)
	{
		return newArray(bounds, 1f);
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#removeBias(java.awt.Rectangle, float[])
	 */
	public void removeBias(Rectangle bounds, float[] data)
	{
		for (int i = 0; i < data.length; i++)
			data[i] -= bias;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#removeGain(java.awt.Rectangle, float[])
	 */
	public void removeGain(Rectangle bounds, float[] data)
	{
		for (int i = 0; i < data.length; i++)
			data[i] /= gain;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#removeBiasAndRemoveGain(java.awt.Rectangle, float[])
	 */
	public void removeBiasAndGain(Rectangle bounds, float[] data)
	{
		for (int i = 0; i < data.length; i++)
			data[i] = (data[i] - bias) / gain;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.BaseCameraModel#copy()
	 */
	@Override
	public FixedPixelCameraModel copy()
	{
		return clone();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	protected FixedPixelCameraModel clone()
	{
		try
		{
			return (FixedPixelCameraModel) super.clone();
		}
		catch (CloneNotSupportedException e)
		{
			return null;
		}
	}
}
