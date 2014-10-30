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
public class EMCCDCameraNoiseModel extends CameraNoiseModel
{
	public EMCCDCameraNoiseModel(final double readNoise)
	{
		super(readNoise);
	}

	public EMCCDCameraNoiseModel(final double readNoise, final double bias)
	{
		super(readNoise, bias);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.function.CameraNoiseModel#variance(double)
	 */
	public double variance(final double value)
	{
		return readNoise2 + FastMath.max(value - bias, 0.0) * 2.0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.function.CameraNoiseModel#isEmCCD()
	 */
	public boolean isEmCCD()
	{
		return true;
	}
}