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
	 */
	public FixedPixelCameraModel(double bias, double gain)
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
	public FixedPixelCameraModel(double bias, double gain, double variance)
	{
		// Re-implement the constructor (rather than chaining)
		// to take advantage of double precision computation of var_g2.

		// Cast to float then check
		this.bias = (float) bias;
		checkBias(this.bias);
		this.gain = (float) gain;
		checkGain(this.gain);
		this.variance = (float) variance;
		checkVariance(this.variance);
		this.var_g2 = (float) (variance / (gain * gain));
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
	 * @see gdsc.smlm.model.camera.CameraModel#setOrigin(int, int)
	 */
	public void setOrigin(int x, int y)
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#crop(java.awt.Rectangle, boolean)
	 */
	public CameraModel crop(Rectangle bounds, boolean resetOrigin)
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
	 * @see gdsc.smlm.model.camera.CameraModel#getNormalisedWeights(java.awt.Rectangle)
	 */
	public float[] getNormalisedWeights(Rectangle bounds)
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
		removeBias(data);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#removeGain(java.awt.Rectangle, float[])
	 */
	public void removeGain(Rectangle bounds, float[] data)
	{
		removeGain(data);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#removeBiasAndRemoveGain(java.awt.Rectangle, float[])
	 */
	public void removeBiasAndGain(Rectangle bounds, float[] data)
	{
		removeBiasAndGain(data);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#applyBias(java.awt.Rectangle, float[])
	 */
	public void applyBias(Rectangle bounds, float[] data)
	{
		applyBias(data);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#applyGain(java.awt.Rectangle, float[])
	 */
	public void applyGain(Rectangle bounds, float[] data)
	{
		applyGain(data);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#applyGainAndBias(java.awt.Rectangle, float[])
	 */
	public void applyGainAndBias(Rectangle bounds, float[] data)
	{
		applyGainAndBias(data);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#removeBias(float[])
	 */
	public void removeBias(float[] data)
	{
		if (data == null)
			return;
		for (int i = 0; i < data.length; i++)
			data[i] -= bias;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#removeGain(float[])
	 */
	public void removeGain(float[] data)
	{
		if (data == null)
			return;
		for (int i = 0; i < data.length; i++)
			data[i] /= gain;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#removeBiasAndRemoveGain(float[])
	 */
	public void removeBiasAndGain(float[] data)
	{
		if (data == null)
			return;
		for (int i = 0; i < data.length; i++)
			data[i] = (data[i] - bias) / gain;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#applyBias(float[])
	 */
	public void applyBias(float[] data)
	{
		if (data == null)
			return;
		for (int i = 0; i < data.length; i++)
			data[i] += bias;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#applyGain(float[])
	 */
	public void applyGain(float[] data)
	{
		if (data == null)
			return;
		for (int i = 0; i < data.length; i++)
			data[i] *= gain;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#applyGainAndBias(float[])
	 */
	public void applyGainAndBias(float[] data)
	{
		if (data == null)
			return;
		for (int i = 0; i < data.length; i++)
			data[i] = data[i] * gain + bias;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#copy()
	 */
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
