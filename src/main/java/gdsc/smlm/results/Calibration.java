package gdsc.smlm.results;

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
 * Contain the calibration settings for the microscope
 */
public class Calibration
{
	/**
	 * The image pixel size in nanometers
	 */
	public double nmPerPixel = 107;
	/**
	 * The gain (ADUs/photon). Can be used to convert the signal in Analogue-to-Digital units
	 * (ADUs) to photons.
	 */
	public float gain = 37.7f;
	/**
	 * The exposure time in milliseconds per frame
	 */
	public double exposureTime;
	
	/**
	 * The CCD camera Gaussian read noise (in ADUs)
	 */
	public float readNoise;
	/**
	 * The CCD camera bias (in ADUs)
	 */
	public float bias;
	/**
	 * True if the CCD camera was run in Electron Multiplying (EM) mode
	 */
	public boolean emCCD;

	/**
	 * Default constructor
	 */
	public Calibration()
	{
	}

	/**
	 * Copy constructor
	 * 
	 * @param calibration
	 */
	public Calibration(Calibration calibration)
	{
		this.nmPerPixel = calibration.nmPerPixel;
		this.gain = calibration.gain;
		this.exposureTime = calibration.exposureTime;
	}

	/**
	 * Parameterised constructor.
	 * <p>
	 * If gain is zero then set to 1.
	 * 
	 * @param nmPerPixel
	 * @param gain
	 * @param exposureTime
	 */
	public Calibration(double nmPerPixel, float gain, double exposureTime)
	{
		this.nmPerPixel = nmPerPixel;
		if (gain <= 0)
			gain = 1;
		this.gain = gain;
		this.exposureTime = exposureTime;
	}
}
