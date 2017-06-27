package gdsc.smlm.data.config;

import gdsc.smlm.data.config.CalibrationConfig.Calibration;
import gdsc.smlm.data.config.CalibrationConfig.CameraType;

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
 * Contains helper functions for the CalibrationConfig class.
 */
public class CalibrationConfigHelper
{
	/** The default Calibration */
	public static final Calibration defaultCalibration;
	static
	{
		Calibration.Builder builder = Calibration.newBuilder();
		defaultCalibration = builder.build();
	}

	/**
	 * Gets the name.
	 *
	 * @param value
	 *            the value
	 * @return the name
	 */
	public static String getName(CameraType value)
	{
		switch (value)
		{
			case CAMERA_TYPE_NA:
				return "NA";
			case CCD:
				return "CCD";
			case EMCCD:
				return "EMCCD";
			case SCMOS:
				return "sCMOS";
			case UNRECOGNIZED:
				return "Unknown";
			default:
				throw new IllegalStateException("Unknown name: " + value);
		}
	}

}
