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

import gdsc.core.utils.SimpleArrayUtils;

/**
 * Define the methods for manipulating camera pixel data.
 *
 * @author Alex Herbert
 */
public class PerPixelCameraModel extends BaseCameraModel
{
	private final Rectangle cameraBounds;

	private final float[] bias, gain, variance;
	// This is computed when required
	private float[] var_g2;

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
		this(new Rectangle(xorigin, yorigin, width, height), bias, gain, variance, false, true);
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
		this(bounds, bias, gain, variance, true, true);
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
	 * @param cloneData
	 *            Set to true to clone the data
	 * @throws IllegalArgumentException
	 *             If the data is not valid
	 */
	private PerPixelCameraModel(Rectangle bounds, float[] bias, float[] gain, float[] variance, boolean cloneBounds,
			boolean cloneData) throws IllegalArgumentException
	{
		if (bounds == null)
			throw new IllegalArgumentException("Bounds must not be null");
		checkBounds(bounds);
		cameraBounds = (cloneBounds) ? new Rectangle(bounds) : bounds;
		int size = SimpleArrayUtils.check2DSize(bounds.width, bounds.height);
		checkArray(bias, size);
		checkArray(gain, size);
		checkArray(variance, size);
		if (cloneData)
		{
			this.bias = bias.clone();
			this.gain = gain.clone();
			this.variance = variance.clone();
		}
		else
		{
			this.bias = bias;
			this.gain = gain;
			this.variance = variance;
		}
		for (int i = 0; i < size; i++)
		{
			checkBias(bias[i]);
			checkGain(gain[i]);
			checkVariance(variance[i]);
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
	private static void checkBounds(Rectangle bounds) throws IllegalArgumentException
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
	 * <p>
	 * This is an internally used copy constructor.
	 *
	 * @param duplicate
	 *            a flag to indicate the data should be duplicated
	 * @param bounds
	 *            the bounds
	 * @param bias
	 *            the bias
	 * @param gain
	 *            the gain
	 * @param variance
	 *            the variance array
	 * @param var_g2
	 *            the normalised variance array
	 */
	private PerPixelCameraModel(boolean duplicate, Rectangle bounds, float[] bias, float[] gain, float[] variance,
			float[] var_g2)
	{
		if (duplicate)
		{
			cameraBounds = new Rectangle(bounds);
			this.bias = bias.clone();
			this.gain = gain.clone();
			this.variance = variance.clone();
			this.var_g2 = (var_g2 == null) ? null : var_g2.clone();
		}
		else
		{
			cameraBounds = bounds;
			this.bias = bias;
			this.gain = gain;
			this.variance = variance;
			this.var_g2 = var_g2;
		}
	}

	/**
	 * Creates a new per pixel camera model. The input arguments are not cloned.
	 *
	 * @param bounds
	 *            the bounds
	 * @param bias
	 *            the bias (in counts)
	 * @param gain
	 *            the gain (count/photon)
	 * @param variance
	 *            the variance (in counts)
	 * @return the per pixel camera model
	 * @throws IllegalArgumentException
	 *             If the data is not valid
	 */
	public static PerPixelCameraModel create(Rectangle bounds, float[] bias, float[] gain, float[] variance)
			throws IllegalArgumentException
	{
		return new PerPixelCameraModel(bounds, bias, gain, variance, false, false);
	}

	/**
	 * Check array.
	 *
	 * @param array
	 *            the array
	 * @param size
	 *            the size
	 */
	private static void checkArray(float[] array, int size)
	{
		if (array == null || array.length != size)
			throw new IllegalArgumentException("Input array must match the size of the input bounds");
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.model.camera.CameraModel#getBounds()
	 */
	@Override
	public Rectangle getBounds()
	{
		return new Rectangle(cameraBounds);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.model.camera.CameraModel#setOrigin(int, int)
	 */
	@Override
	public void setOrigin(int x, int y)
	{
		cameraBounds.x = x;
		cameraBounds.y = y;
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
	 * Gets a copy of the variance for the current bounds.
	 *
	 * @return the variance
	 */
	public float[] getVariance()
	{
		return variance.clone();
	}

	/**
	 * Gets a copy of the normalised variance for the current bounds.
	 *
	 * @return the normalised variance
	 */
	public float[] getNormalisedVariance()
	{
		return getNormalisedVarianceInternal().clone();
	}

	/**
	 * Initialise the model. This allows caching of precomputed values but it is not required to call this method before
	 * using the model.
	 */
	public void initialise()
	{
		getNormalisedVarianceInternal();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.model.camera.CameraModel#isPerPixelModel()
	 */
	@Override
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
	@Override
	public float[] getBias(Rectangle bounds)
	{
		return getData(bounds, bias);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.model.camera.CameraModel#getGain(java.awt.Rectangle)
	 */
	@Override
	public float[] getGain(Rectangle bounds)
	{
		return getData(bounds, gain);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.model.camera.CameraModel#getVariance(java.awt.Rectangle)
	 */
	@Override
	public float[] getVariance(Rectangle bounds)
	{
		return getData(bounds, variance);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.model.camera.CameraModel#getNormalisedVariance(java.awt.Rectangle)
	 */
	@Override
	public float[] getNormalisedVariance(Rectangle bounds)
	{
		return getData(bounds, getNormalisedVarianceInternal());
	}

	private float[] getNormalisedVarianceInternal()
	{
		if (var_g2 == null)
		{
			createNormalisedVariance();
		}
		return var_g2;
	}

	private synchronized void createNormalisedVariance()
	{
		if (var_g2 == null)
		{
			int size = variance.length;
			var_g2 = new float[size];
			for (int i = 0; i < size; i++)
			{
				var_g2[i] = variance[i] / (gain[i] * gain[i]);
			}
		}
	}

	@Override
	public float getBias(int x, int y)
	{
		return getData(x, y, bias);
	}

	@Override
	public float getGain(int x, int y)
	{
		return getData(x, y, gain);
	}

	@Override
	public float getVariance(int x, int y)
	{
		return getData(x, y, variance);
	}

	@Override
	public float getNormalisedVariance(int x, int y)
	{
		return getData(x, y, getNormalisedVarianceInternal());
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.model.camera.CameraModel#getMeanVariance(java.awt.Rectangle)
	 */
	@Override
	public double getMeanVariance(Rectangle bounds)
	{
		return getMean(bounds, variance);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.model.camera.CameraModel#getMeanNormalisedVariance(java.awt.Rectangle)
	 */
	@Override
	public double getMeanNormalisedVariance(Rectangle bounds)
	{
		return getMean(bounds, getNormalisedVarianceInternal());
	}

	/**
	 * Return the weights as 1/variance.
	 *
	 * @see gdsc.smlm.model.camera.CameraModel#getWeights(java.awt.Rectangle)
	 */
	@Override
	public float[] getWeights(Rectangle bounds)
	{
		return toWeights(getVariance(bounds));
	}

	/**
	 * Return the weights as 1/[normalised variance].
	 *
	 * @see gdsc.smlm.model.camera.CameraModel#getNormalisedWeights(java.awt.Rectangle)
	 */
	@Override
	public float[] getNormalisedWeights(Rectangle bounds)
	{
		return toWeights(getNormalisedVariance(bounds));
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
	 * Gets the mean of the data from the values using the intersection of the bounds.
	 *
	 * @param bounds
	 *            the bounds
	 * @param values
	 *            the values
	 * @return the mean of the data
	 */
	private double getMean(Rectangle bounds, float[] values)
	{
		Rectangle intersection = getIntersection(bounds);
		return getMean(values, intersection);
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
		// We avoid underflow since we have checked the bounds are positive integers
		int minx = bounds.x - cameraBounds.x;
		int miny = bounds.y - cameraBounds.y;
		if (minx < 0 || miny < 0)
			throw new IllegalArgumentException("Bounds must be within the camera bounds");
		// Avoid overflow using a long result
		long maxx = (long) minx + bounds.width;
		long maxy = (long) miny + bounds.height;
		if (maxx > cameraBounds.width || maxy > cameraBounds.height)
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

	/**
	 * Get the mean of the data from the per-pixel data using the given bounds.
	 *
	 * @param pixels
	 *            the pixels
	 * @param bounds
	 *            the bounds
	 * @return the mean of the data
	 */
	private double getMean(float[] pixels, final Rectangle bounds)
	{
		double sum = 0;
		if (bounds.x != 0 || bounds.y != 0 || bounds.width != cameraBounds.width ||
				bounds.height != cameraBounds.height)
		{
			final int width = cameraBounds.width;
			for (int ys = 0; ys < bounds.height; ys++)
			{
				for (int xs = 0, offset2 = (ys + bounds.y) * width + bounds.x; xs < bounds.width; xs++)
					sum += pixels[offset2++];
			}
			return sum / (bounds.height * bounds.width);
		}
		else
		{
			for (int i = pixels.length; i-- > 0;)
				sum += pixels[i];
			return sum / pixels.length;
		}
	}

	/**
	 * Gets the data value.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param data
	 *            the data
	 * @return the data value
	 * @throws IllegalArgumentException
	 *             If the coordinates are not inside the bounds
	 */
	private float getData(int x, int y, float[] data) throws IllegalArgumentException
	{
		x -= cameraBounds.x;
		y -= cameraBounds.y;
		if (x < 0 || y < 0 || x >= cameraBounds.width || y >= cameraBounds.height)
			throw new IllegalArgumentException("Coordinates must be within the camera bounds");
		return data[y * cameraBounds.width + x];
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.model.camera.CameraModel#removeBias(java.awt.Rectangle, float[])
	 */
	@Override
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
	@Override
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
	@Override
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
	 * @see gdsc.smlm.model.camera.CameraModel#applyBias(java.awt.Rectangle, float[])
	 */
	@Override
	public void applyBias(Rectangle bounds, float[] data)
	{
		if (data == null)
			return;
		Rectangle intersection = getIntersection(bounds);
		float[] bias = getData(this.bias, intersection, false);
		if (data.length != bias.length)
			throw new IllegalArgumentException("Bounds must match the data size");
		for (int i = 0; i < data.length; i++)
			data[i] += bias[i];
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.model.camera.CameraModel#applyGain(java.awt.Rectangle, float[])
	 */
	@Override
	public void applyGain(Rectangle bounds, float[] data)
	{
		if (data == null)
			return;
		Rectangle intersection = getIntersection(bounds);
		float[] gain = getData(this.gain, intersection, false);
		if (data.length != gain.length)
			throw new IllegalArgumentException("Bounds must match the data size");
		for (int i = 0; i < data.length; i++)
			data[i] *= gain[i];
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.model.camera.CameraModel#applyGainAndBias(java.awt.Rectangle, float[])
	 */
	@Override
	public void applyGainAndBias(Rectangle bounds, float[] data)
	{
		if (data == null)
			return;
		Rectangle intersection = getIntersection(bounds);
		float[] bias = getData(this.bias, intersection, false);
		if (data.length != bias.length)
			throw new IllegalArgumentException("Bounds must match the data size");
		float[] gain = getData(this.gain, intersection, false);
		for (int i = 0; i < data.length; i++)
			data[i] = data[i] * gain[i] + bias[i];
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.model.camera.CameraModel#removeBias(float[])
	 */
	@Override
	public void removeBias(float[] data)
	{
		if (data == null)
			return;
		if (data.length != bias.length)
			throw new IllegalArgumentException("Data length must match the bounds (width x height)");
		for (int i = 0; i < data.length; i++)
			data[i] -= bias[i];
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.model.camera.CameraModel#removeGain(float[])
	 */
	@Override
	public void removeGain(float[] data)
	{
		if (data == null)
			return;
		if (data.length != gain.length)
			throw new IllegalArgumentException("Data length must match the bounds (width x height)");
		for (int i = 0; i < data.length; i++)
			data[i] /= gain[i];
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.model.camera.CameraModel#removeBiasAndRemoveGain(float[])
	 */
	@Override
	public void removeBiasAndGain(float[] data)
	{
		if (data == null)
			return;
		if (data.length != bias.length)
			throw new IllegalArgumentException("Data length must match the bounds (width x height)");
		for (int i = 0; i < data.length; i++)
			data[i] = (data[i] - bias[i]) / gain[i];
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.model.camera.CameraModel#applyBias(float[])
	 */
	@Override
	public void applyBias(float[] data)
	{
		if (data == null)
			return;
		if (data.length != bias.length)
			throw new IllegalArgumentException("Data length must match the bounds (width x height)");
		for (int i = 0; i < data.length; i++)
			data[i] += bias[i];
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.model.camera.CameraModel#applyGain(float[])
	 */
	@Override
	public void applyGain(float[] data)
	{
		if (data == null)
			return;
		if (data.length != gain.length)
			throw new IllegalArgumentException("Data length must match the bounds (width x height)");
		for (int i = 0; i < data.length; i++)
			data[i] *= gain[i];
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.model.camera.CameraModel#applyGainAndBias(float[])
	 */
	@Override
	public void applyGainAndBias(float[] data)
	{
		if (data == null)
			return;
		if (data.length != bias.length)
			throw new IllegalArgumentException("Data length must match the bounds (width x height)");
		for (int i = 0; i < data.length; i++)
			data[i] = data[i] * gain[i] + bias[i];
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.model.camera.CameraModel#crop(java.awt.Rectangle, boolean)
	 */
	@Override
	public CameraModel crop(Rectangle bounds, boolean resetOrigin)
	{
		if (bounds == null)
			throw new IllegalArgumentException("Bounds are null");
		if (equal(bounds))
		{
			if (resetOrigin)
			{
				PerPixelCameraModel model = copy();
				model.setOrigin(0, 0);
				return model;
			}
			return this;
		}
		Rectangle intersection = getIntersection(bounds);
		float[] bias = getData(this.bias, intersection, true);
		float[] gain = getData(this.gain, intersection, true);
		float[] variance = getData(this.variance, intersection, true);
		float[] var_g2 = (this.var_g2 == null) ? null : getData(this.var_g2, intersection, true);
		return new PerPixelCameraModel(false, bounds, bias, gain, variance, var_g2);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.model.camera.CameraModel#copy()
	 */
	@Override
	public PerPixelCameraModel copy()
	{
		// Deep copy
		return new PerPixelCameraModel(true, cameraBounds, bias, gain, variance, var_g2);
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
