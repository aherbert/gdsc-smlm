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
 * Define the methods for manipulating camera pixel data.
 *
 * @author Alex Herbert
 */
public class PerPixelCameraModel extends BaseCameraModel
{
	private final Rectangle cameraBounds;

	private final float[] bias, gain, var_g2;

	/**
	 * Instantiates a new per pixel camera model.
	 *
	 * @param width
	 *            the width
	 * @param height
	 *            the height
	 * @param bias
	 *            the bias (in counts)
	 * @param gain
	 *            the gain (count/photon)
	 * @param variance
	 *            the variance (in counts)
	 * @throws IllegalArgumentException
	 *             If the data is not valid
	 */
	public PerPixelCameraModel(int width, int height, float[] bias, float[] gain, float[] variance)
			throws IllegalArgumentException
	{
		this(0, 0, width, height, bias, gain, variance);
	}

	/**
	 * Instantiates a new per pixel camera model.
	 *
	 * @param xorigin
	 *            the xorigin
	 * @param yorigin
	 *            the yorigin
	 * @param width
	 *            the width
	 * @param height
	 *            the height
	 * @param bias
	 *            the bias (in counts)
	 * @param gain
	 *            the gain (count/photon)
	 * @param variance
	 *            the variance (in counts)
	 * @throws IllegalArgumentException
	 *             If the data is not valid
	 */
	public PerPixelCameraModel(int xorigin, int yorigin, int width, int height, float[] bias, float[] gain,
			float[] variance) throws IllegalArgumentException
	{
		this(new Rectangle(xorigin, yorigin, width, height), bias, gain, variance, false);
	}

	/**
	 * Instantiates a new per pixel camera model.
	 *
	 * @param bounds
	 *            the bounds
	 * @param bias
	 *            the bias (in counts)
	 * @param gain
	 *            the gain (count/photon)
	 * @param variance
	 *            the variance (in counts)
	 * @throws IllegalArgumentException
	 *             If the data is not valid
	 */
	public PerPixelCameraModel(Rectangle bounds, float[] bias, float[] gain, float[] variance)
			throws IllegalArgumentException
	{
		this(bounds, bias, gain, variance, true);

	}

	/**
	 * Instantiates a new per pixel camera model.
	 *
	 * @param bounds
	 *            the bounds
	 * @param bias
	 *            the bias (in counts)
	 * @param gain
	 *            the gain (count/photon)
	 * @param variance
	 *            the variance (in counts)
	 * @param cloneBounds
	 *            Set to true to clone the bounds
	 * @throws IllegalArgumentException
	 *             If the data is not valid
	 */
	private PerPixelCameraModel(Rectangle bounds, float[] bias, float[] gain, float[] variance, boolean cloneBounds)
			throws IllegalArgumentException
	{
		if (bounds == null)
			throw new IllegalArgumentException("Bounds must not be null");
		checkBounds(bounds);
		cameraBounds = (cloneBounds) ? new Rectangle(bounds) : bounds;
		int size = bounds.width * bounds.height;
		checkArray(bias, size);
		checkArray(gain, size);
		checkArray(variance, size);
		this.bias = bias.clone();
		this.gain = gain.clone();
		this.var_g2 = new float[size];
		for (int i = 0; i < size; i++)
		{
			checkBias(bias[i]);
			checkGain(gain[i]);
			checkVariance(variance[i]);
			var_g2[i] = variance[i] / (gain[i] * gain[i]);
		}
	}

	/**
	 * Check the bounds have positive coordinates and widths.
	 *
	 * @param bounds
	 *            the bounds
	 * @throws IllegalArgumentException
	 *             If the data is not valid
	 */
	private void checkBounds(Rectangle bounds) throws IllegalArgumentException
	{
		if (bounds.x < 0)
			throw new IllegalArgumentException("Bounds must have positive x origin");
		if (bounds.y < 0)
			throw new IllegalArgumentException("Bounds must have positive y origin");
		if (bounds.width < 0)
			throw new IllegalArgumentException("Bounds must have positive width");
		if (bounds.height < 0)
			throw new IllegalArgumentException("Bounds must have positive height");
	}

	/**
	 * Instantiates a new per pixel camera model, copying all input fields.
	 *
	 * @param duplicate
	 *            a flag to indicate the data should be duplicated
	 * @param bounds
	 *            the bounds
	 * @param bias
	 *            the bias
	 * @param gain
	 *            the gain
	 * @param var_g2
	 *            the var_g2 array
	 */
	private PerPixelCameraModel(boolean duplicate, Rectangle bounds, float[] bias, float[] gain, float[] var_g2)
	{
		if (duplicate)
		{
			// This is used privately so we just copy the input
			cameraBounds = new Rectangle(bounds);
			this.bias = bias.clone();
			this.gain = gain.clone();
			this.var_g2 = var_g2.clone();
		}
		else
		{
			// This is used publicly so we check the input data
			checkBounds(bounds);
			cameraBounds = bounds;
			int size = bounds.width * bounds.height;
			checkArray(bias, size);
			checkArray(gain, size);
			checkArray(var_g2, size);
			this.bias = bias;
			this.gain = gain;
			this.var_g2 = var_g2;
			for (int i = 0; i < size; i++)
			{
				checkBias(bias[i]);
				checkGain(gain[i]);
				checkVariance(var_g2[i]);
			}
		}
	}

	/**
	 * Creates a new per pixel camera model.
	 *
	 * @param bounds
	 *            the bounds
	 * @param bias
	 *            the bias (in counts)
	 * @param gain
	 *            the gain (count/photon)
	 * @param normalisedVariance
	 *            the normalised variance (in photons)
	 * @return the per pixel camera model
	 * @throws IllegalArgumentException
	 *             If the data is not valid
	 */
	public static PerPixelCameraModel create(Rectangle bounds, float[] bias, float[] gain, float[] normalisedVariance)
			throws IllegalArgumentException
	{
		return new PerPixelCameraModel(false, bounds, bias, gain, normalisedVariance);
	}

	/**
	 * Check array.
	 *
	 * @param array
	 *            the array
	 * @param size
	 *            the size
	 */
	private void checkArray(float[] array, int size)
	{
		if (array == null || array.length != size)
			throw new IllegalArgumentException("Input array must match the size of the input bounds");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#getBounds()
	 */
	public Rectangle getBounds()
	{
		return new Rectangle(cameraBounds);
	}

	/**
	 * Gets a copy of the bias for the current bounds.
	 *
	 * @return the bias
	 */
	public float[] getBias()
	{
		return bias.clone();
	}

	/**
	 * Gets a copy of the gain for the current bounds.
	 *
	 * @return the gain
	 */
	public float[] getGain()
	{
		return gain.clone();
	}

	/**
	 * Gets a copy of the normalised variance for the current bounds.
	 *
	 * @return the normalised variance
	 */
	public float[] getNormalisedVariance()
	{
		return var_g2.clone();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#isPerPixelModel()
	 */
	public boolean isPerPixelModel()
	{
		return true;
	}

	/**
	 * Gets the x origin.
	 *
	 * @return the x origin
	 */
	public int getXOrigin()
	{
		return cameraBounds.x;
	}

	/**
	 * Gets the y origin.
	 *
	 * @return the y origin
	 */
	public int getYOrigin()
	{
		return cameraBounds.y;
	}

	/**
	 * Gets the width.
	 *
	 * @return the width
	 */
	public int getWidth()
	{
		return cameraBounds.width;
	}

	/**
	 * Gets the height.
	 *
	 * @return the height
	 */
	public int getHeight()
	{
		return cameraBounds.height;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#getBias(java.awt.Rectangle)
	 */
	public float[] getBias(Rectangle bounds)
	{
		return getData(bounds, bias);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#getGain(java.awt.Rectangle)
	 */
	public float[] getGain(Rectangle bounds)
	{
		return getData(bounds, gain);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#getNormalisedVariance(java.awt.Rectangle)
	 */
	public float[] getNormalisedVariance(Rectangle bounds)
	{
		return getData(bounds, var_g2);
	}

	/**
	 * Gets the data from the values using the intersection of the bounds.
	 *
	 * @param bounds
	 *            the bounds
	 * @param values
	 *            the values
	 * @return the data
	 */
	private float[] getData(Rectangle bounds, float[] values)
	{
		Rectangle intersection = getIntersection(bounds);
		return getData(values, intersection, true);
	}

	/**
	 * Gets the intersection between the target bounds and the camera bounds. The results is offset by the camera bounds
	 * origin so can be used to crop the camera per-pixel data.
	 *
	 * @param bounds
	 *            the bounds
	 * @return the intersection
	 */
	private Rectangle getIntersection(Rectangle bounds)
	{
		if (bounds == null)
			throw new IllegalArgumentException("Bounds are null");
		if (equal(bounds))
			return new Rectangle(getWidth(), getHeight()); // Offset to the origin
		checkBounds(bounds);
		// We avoid over/underflow since we have checked the bounds are positive integers
		int minx = bounds.x - cameraBounds.x;
		int miny = bounds.y - cameraBounds.y;
		int maxx = minx + bounds.width;
		int maxy = miny + bounds.height;
		if (minx < 0 || miny < 0 || maxx > cameraBounds.width || maxy > cameraBounds.height)
			throw new IllegalArgumentException("Bounds must be within the camera bounds");
		return new Rectangle(minx, miny, bounds.width, bounds.height);

		//			// Using java.awt.Rectangle.intesection function
		//			Rectangle intersection = cameraBounds.intersection(bounds);
		//			if (intersection.width != bounds.width || intersection.height != bounds.height)
		//				throw new IllegalArgumentException("Bounds must be within the camera bounds");
		//			intersection.x -= cameraBounds.x;
		//			intersection.y -= cameraBounds.y;
		//			return intersection;
	}

	/**
	 * Check if the bounds are equal to the camera bounds.
	 *
	 * @param bounds
	 *            the bounds
	 * @return true, if successful
	 */
	private boolean equal(Rectangle bounds)
	{
		//@formatter:off
		return cameraBounds.x == bounds.x && 
			   cameraBounds.y == bounds.y && 
			   cameraBounds.width == bounds.width &&
			   cameraBounds.height == bounds.height;
		//@formatter:on
	}

	/**
	 * Crop the data from the per-pixel data using the given bounds.
	 *
	 * @param pixels
	 *            the pixels
	 * @param bounds
	 *            the bounds
	 * @param copy
	 *            Set to true to copy the values (if the bounds cover all the pixel data)
	 * @return the data
	 */
	private float[] getData(float[] pixels, final Rectangle bounds, boolean copy)
	{
		if (bounds.x != 0 || bounds.y != 0 || bounds.width != cameraBounds.width ||
				bounds.height != cameraBounds.height)
		{
			float[] pixels2 = new float[bounds.width * bounds.height];
			final int width = cameraBounds.width;
			for (int ys = 0, offset1 = 0; ys < bounds.height; ys++)
			{
				for (int xs = 0, offset2 = (ys + bounds.y) * width + bounds.x; xs < bounds.width; xs++)
					pixels2[offset1++] = pixels[offset2++];
			}
			return pixels2;
		}
		else
		{
			return (copy) ? pixels.clone() : pixels;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#removeBias(java.awt.Rectangle, float[])
	 */
	public void removeBias(Rectangle bounds, float[] data)
	{
		if (data == null)
			return;
		Rectangle intersection = getIntersection(bounds);
		float[] bias = getData(this.bias, intersection, false);
		if (data.length != bias.length)
			throw new IllegalArgumentException("Bounds must match the data size");
		for (int i = 0; i < data.length; i++)
			data[i] -= bias[i];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#removeGain(java.awt.Rectangle, float[])
	 */
	public void removeGain(Rectangle bounds, float[] data)
	{
		if (data == null)
			return;
		Rectangle intersection = getIntersection(bounds);
		float[] gain = getData(this.gain, intersection, false);
		if (data.length != gain.length)
			throw new IllegalArgumentException("Bounds must match the data size");
		for (int i = 0; i < data.length; i++)
			data[i] /= gain[i];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#removeBiasAndRemoveGain(java.awt.Rectangle, float[])
	 */
	public void removeBiasAndGain(Rectangle bounds, float[] data)
	{
		if (data == null)
			return;
		Rectangle intersection = getIntersection(bounds);
		float[] bias = getData(this.bias, intersection, false);
		if (data.length != bias.length)
			throw new IllegalArgumentException("Bounds must match the data size");
		float[] gain = getData(this.gain, intersection, false);
		for (int i = 0; i < data.length; i++)
			data[i] = (data[i] - bias[i]) / gain[i];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.CameraModel#crop(java.awt.Rectangle)
	 */
	public CameraModel crop(Rectangle bounds)
	{
		if (bounds == null)
			throw new IllegalArgumentException("Bounds are null");
		if (equal(bounds))
			return this;
		Rectangle intersection = getIntersection(bounds);
		float[] bias = getData(this.bias, intersection, true);
		float[] gain = getData(this.gain, intersection, true);
		float[] var_g2 = getData(this.var_g2, intersection, true);
		return new PerPixelCameraModel(true, bounds, bias, gain, var_g2);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.camera.BaseCameraModel#copy()
	 */
	public PerPixelCameraModel copy()
	{
		return new PerPixelCameraModel(true, cameraBounds, bias, gain, var_g2);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	protected PerPixelCameraModel clone()
	{
		try
		{
			return (PerPixelCameraModel) super.clone();
		}
		catch (CloneNotSupportedException e)
		{
			return null;
		}
	}
}
