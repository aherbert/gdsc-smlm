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
package uk.ac.sussex.gdsc.smlm.model.camera;

import java.awt.Rectangle;

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
	 * @see uk.ac.sussex.gdsc.smlm.model.camera.CameraModel#getBounds()
	 */
	@Override
	public Rectangle getBounds()
	{
		return null;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.camera.CameraModel#setOrigin(int, int)
	 */
	@Override
	public void setOrigin(int x, int y)
	{
		// Ignore
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.camera.CameraModel#crop(java.awt.Rectangle, boolean)
	 */
	@Override
	public CameraModel crop(Rectangle bounds, boolean resetOrigin)
	{
		return this;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.camera.CameraModel#isPerPixelModel()
	 */
	@Override
	public boolean isPerPixelModel()
	{
		return false;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.camera.CameraModel#getBias(java.awt.Rectangle)
	 */
	@Override
	public float[] getBias(Rectangle bounds)
	{
		return newArray(bounds, 0);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.camera.CameraModel#getGain(java.awt.Rectangle)
	 */
	@Override
	public float[] getGain(Rectangle bounds)
	{
		return newArray(bounds, 1f);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.camera.CameraModel#getVariance(java.awt.Rectangle)
	 */
	@Override
	public float[] getVariance(Rectangle bounds)
	{
		return newArray(bounds, 0);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.camera.CameraModel#getNormalisedVariance(java.awt.Rectangle)
	 */
	@Override
	public float[] getNormalisedVariance(Rectangle bounds)
	{
		return newArray(bounds, 0);
	}

	@Override
	public float getBias(int x, int y)
	{
		return 0f;
	}

	@Override
	public float getGain(int x, int y)
	{
		return 1f;
	}

	@Override
	public float getVariance(int x, int y)
	{
		return 0f;
	}

	@Override
	public float getNormalisedVariance(int x, int y)
	{
		return 0f;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.camera.CameraModel#getMeanVariance(java.awt.Rectangle)
	 */
	@Override
	public double getMeanVariance(Rectangle bounds)
	{
		return 0d;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.camera.CameraModel#getMeanNormalisedVariance(java.awt.Rectangle)
	 */
	@Override
	public double getMeanNormalisedVariance(Rectangle bounds)
	{
		return 0d;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.camera.CameraModel#getWeights(java.awt.Rectangle)
	 */
	@Override
	public float[] getWeights(Rectangle bounds)
	{
		return newArray(bounds, 1f);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.camera.CameraModel#getNormalisedWeights(java.awt.Rectangle)
	 */
	@Override
	public float[] getNormalisedWeights(Rectangle bounds)
	{
		return newArray(bounds, 1f);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.camera.CameraModel#removeBias(java.awt.Rectangle, float[])
	 */
	@Override
	public void removeBias(Rectangle bounds, float[] data)
	{
		// Ignore
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.camera.CameraModel#removeGain(java.awt.Rectangle, float[])
	 */
	@Override
	public void removeGain(Rectangle bounds, float[] data)
	{
		// Ignore
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.camera.CameraModel#removeBiasAndRemoveGain(java.awt.Rectangle, float[])
	 */
	@Override
	public void removeBiasAndGain(Rectangle bounds, float[] data)
	{
		// Ignore
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.camera.CameraModel#applyBias(java.awt.Rectangle, float[])
	 */
	@Override
	public void applyBias(Rectangle bounds, float[] data)
	{
		// Ignore
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.camera.CameraModel#applyGain(java.awt.Rectangle, float[])
	 */
	@Override
	public void applyGain(Rectangle bounds, float[] data)
	{
		// Ignore
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.camera.CameraModel#applyGainAndBias(java.awt.Rectangle, float[])
	 */
	@Override
	public void applyGainAndBias(Rectangle bounds, float[] data)
	{
		// Ignore
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.camera.CameraModel#removeBias(float[])
	 */
	@Override
	public void removeBias(float[] data)
	{
		// Ignore
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.camera.CameraModel#removeGain(float[])
	 */
	@Override
	public void removeGain(float[] data)
	{
		// Ignore
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.camera.CameraModel#removeBiasAndGain(float[])
	 */
	@Override
	public void removeBiasAndGain(float[] data)
	{
		// Ignore
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.camera.CameraModel#applyBias(float[])
	 */
	@Override
	public void applyBias(float[] data)
	{
		// Ignore
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.camera.CameraModel#applyGain(float[])
	 */
	@Override
	public void applyGain(float[] data)
	{
		// Ignore
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.camera.CameraModel#applyGainAndBias(float[])
	 */
	@Override
	public void applyGainAndBias(float[] data)
	{
		// Ignore
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.model.camera.CameraModel#copy()
	 */
	@Override
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
		catch (final CloneNotSupportedException e)
		{
			return null;
		}
	}
}
