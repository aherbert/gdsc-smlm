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
 * Define the methods for manipulating camera pixel data
 * 
 * @author Alex Herbert
 */
public interface CameraModel
{
	/**
	 * Gets the bounds of the camera pixel data. This could be null if the camera has infinite bounds, e.g. all pixels
	 * are treated equally.
	 *
	 * @return the bounds
	 */
	public Rectangle getBounds();

	/**
	 * Crop the camera to the given bounds. The bounds are expected to fit within the camera bounds.
	 * <p>
	 * This can be used to create a more efficient representation if no data outside the bounds are required.
	 * <p>
	 * Note: If the bounds match the current bounds then the returned model may not be a copy.
	 *
	 * @param bounds
	 *            the bounds
	 * @return the camera model
	 */
	public CameraModel crop(Rectangle bounds);

	/**
	 * Checks if is per pixel model. If false then all pixels are treated equally.
	 *
	 * @return true, if is per pixel model
	 */
	public boolean isPerPixelModel();

	/**
	 * Gets the per-pixel camera bias (offset). The bounds are expected to fit within the camera bounds.
	 *
	 * @param bounds
	 *            the bounds
	 * @return the bias
	 */
	public float[] getBias(Rectangle bounds);

	/**
	 * Gets the per-pixel camera gain (in count/photon). The bounds are expected to fit within the camera bounds.
	 *
	 * @param bounds
	 *            the bounds
	 * @return the gain
	 */
	public float[] getGain(Rectangle bounds);

	/**
	 * Gets the per-pixel variance. This is the variance of the pixel in camera counts. The bounds are expected to fit
	 * within the camera bounds.
	 *
	 * @param bounds
	 *            the bounds
	 * @return the variance
	 */
	public float[] getVariance(Rectangle bounds);

	/**
	 * Gets the per-pixel normalised variance. This is the variance of the pixel in camera counts divided by the squared
	 * gain, i.e. the variance in photon units. The bounds are expected to fit within the camera bounds.
	 *
	 * @param bounds
	 *            the bounds
	 * @return the normalised variance
	 */
	public float[] getNormalisedVariance(Rectangle bounds);

	/**
	 * Remove the per-pixel camera bias (offset) from the crop of the camera data. The bounds are expected to fit
	 * within the camera bounds.
	 *
	 * @param bounds
	 *            the bounds of the data.
	 * @param data
	 *            the data
	 */
	public void removeBias(Rectangle bounds, float[] data);

	/**
	 * Remove the per-pixel gain from the crop of the camera data. The bounds are expected to fit within the camera
	 * bounds.
	 *
	 * @param bounds
	 *            the bounds of the data.
	 * @param data
	 *            the data
	 */
	public void removeGain(Rectangle bounds, float[] data);

	/**
	 * Remove the per-pixel camera bias (offset) and gain from the crop of the camera data. The bounds are expected to
	 * fit within the camera bounds.
	 *
	 * @param bounds
	 *            the bounds of the data.
	 * @param data
	 *            the data
	 */
	public void removeBiasAndGain(Rectangle bounds, float[] data);
}
