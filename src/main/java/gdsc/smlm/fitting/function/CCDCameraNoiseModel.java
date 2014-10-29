package gdsc.smlm.fitting.function;

import org.apache.commons.math3.util.FastMath;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Defines the expected variance of a signal recorded on a CCD or EM-CCD Camera. The model assumes a Gaussian read
 * noise, photon shot noise and an EM-gain noise factor.
 */
public class CCDCameraNoiseModel extends CameraNoiseModel
{
	public CCDCameraNoiseModel(final float readNoise)
	{
		super(readNoise);
	}

	public CCDCameraNoiseModel(final float readNoise, final float bias)
	{
		super(readNoise, bias);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.function.CameraNoiseModel#variance(float)
	 */
	public float variance(final float value)
	{
		return readNoise2 + FastMath.max(value - bias, 0f);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.function.CameraNoiseModel#isEmCCD()
	 */
	public boolean isEmCCD()
	{
		return false;
	}
}
