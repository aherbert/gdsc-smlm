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
public class NullCameraModel extends BaseCameraModel
{
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
		return newArray(bounds, 0);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#getGain(java.awt.Rectangle)
	 */
	public float[] getGain(Rectangle bounds)
	{
		return newArray(bounds, 1f);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#getVariance(java.awt.Rectangle)
	 */
	public float[] getVariance(Rectangle bounds)
	{
		return newArray(bounds, 0);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#getNormalisedVariance(java.awt.Rectangle)
	 */
	public float[] getNormalisedVariance(Rectangle bounds)
	{
		return newArray(bounds, 0);
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
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#removeGain(java.awt.Rectangle, float[])
	 */
	public void removeGain(Rectangle bounds, float[] data)
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#removeBiasAndRemoveGain(java.awt.Rectangle, float[])
	 */
	public void removeBiasAndGain(Rectangle bounds, float[] data)
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#applyBias(java.awt.Rectangle, float[])
	 */
	public void applyBias(Rectangle bounds, float[] data)
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#applyGain(java.awt.Rectangle, float[])
	 */
	public void applyGain(Rectangle bounds, float[] data)
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#applyGainAndBias(java.awt.Rectangle, float[])
	 */
	public void applyGainAndBias(Rectangle bounds, float[] data)
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#removeBias(float[])
	 */
	public void removeBias(float[] data)
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#removeGain(float[])
	 */
	public void removeGain(float[] data)
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#removeBiasAndGain(float[])
	 */
	public void removeBiasAndGain(float[] data)
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#applyBias(float[])
	 */
	public void applyBias(float[] data)
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#applyGain(float[])
	 */
	public void applyGain(float[] data)
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#applyGainAndBias(float[])
	 */
	public void applyGainAndBias(float[] data)
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#copy()
	 */
	public NullCameraModel copy()
	{
		return this; // no state so no need to clone()
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	protected NullCameraModel clone()
	{
		try
		{
			return (NullCameraModel) super.clone();
		}
		catch (CloneNotSupportedException e)
		{
			return null;
		}
	}
}
