package gdsc.smlm.model.camera;

import gdsc.core.utils.Maths;

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
 * Base class for the camera model
 *
 * @author Alex Herbert
 */
public abstract class BaseCameraModel implements CameraModel, Cloneable
{
	public void checkBias(float bias)
	{
		if (!Maths.isFinite(bias))
			throw new IllegalArgumentException("Bias must be a finite number");
	}
	
	public void checkGain(float gain)
	{
		if (!(gain <= Double.MAX_VALUE && gain > 0))
			throw new IllegalArgumentException("Gain must be strictly positive");
	}
	
	public void checkVariance(float variance)
	{
		if (!(variance <= Double.MAX_VALUE && variance >= 0))
			throw new IllegalArgumentException("Variance must be positive");
	}
	
	/**
	 * Copy this camera model. This is a deep copy of any structures.
	 *
	 * @return the base camera model
	 */
	public abstract BaseCameraModel copy();
}
