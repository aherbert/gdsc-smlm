package gdsc.smlm.model.camera;

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
 * A per-pixel camera model with all pixels treated equally. Note that this concept is invalid since this model reports
 * itself as a per-pixel model even though all pixels are treated equally. This allows testing algorithms that require a
 * per-pixel model with a fixed-pixel model, e.g. for performance comparison.
 *
 * @author Alex Herbert
 */
public class FakePerPixelCameraModel extends FixedPixelCameraModel
{
	/**
	 * Instantiates a new fake per pixel camera model.
	 *
	 * @param bias
	 *            the bias (in counts)
	 * @param gain
	 *            the gain (count/photon)
	 */
	public FakePerPixelCameraModel(float bias, float gain)
	{
		super(bias, gain);
	}

	/**
	 * Instantiates a new fake per pixel camera model.
	 *
	 * @param bias
	 *            the bias (in counts)
	 * @param gain
	 *            the gain (count/photon)
	 * @param variance
	 *            the variance (in counts)
	 */
	public FakePerPixelCameraModel(float bias, float gain, float variance)
	{
		super(bias, gain, variance);
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
}
