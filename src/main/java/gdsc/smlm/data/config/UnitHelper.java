package gdsc.smlm.data.config;

import java.util.ArrayList;

import gdsc.smlm.data.config.SMLMSettings.AngleUnit;
import gdsc.smlm.data.config.SMLMSettings.Calibration;
import gdsc.smlm.data.config.SMLMSettings.DistanceCalibration;
import gdsc.smlm.data.config.SMLMSettings.DistanceUnit;
import gdsc.smlm.data.config.SMLMSettings.IntensityUnit;
import gdsc.smlm.data.utils.ConversionException;
import gdsc.smlm.data.utils.IdentityTypeConverter;
import gdsc.smlm.data.utils.TypeConverter;

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
 * Contains helper functions for the units.
 */
public class UnitHelper
{
	/**
	 * Gets the name.
	 *
	 * @param unit
	 *            the unit
	 * @return the name
	 */
	public static String getName(DistanceUnit unit)
	{
		switch (unit)
		{
			case NM:
				return "nanometer";
			case PIXEL:
				return "pixel";
			case UM:
				return "micrometer";
			default:
				return "unknown";
		}
	}

	/**
	 * Gets the name.
	 *
	 * @param unit
	 *            the unit
	 * @return the name
	 */
	public static String getName(IntensityUnit unit)
	{
		switch (unit)
		{
			case COUNT:
				return "count";
			case PHOTON:
				return "photon";
			default:
				return "unknown";
		}
	}

	/**
	 * Gets the name.
	 *
	 * @param unit
	 *            the unit
	 * @return the name
	 */
	public static String getName(AngleUnit unit)
	{
		switch (unit)
		{
			case DEGREE:
				return "degree";
			case RADIAN:
				return "radian";
			default:
				return "unknown";
		}
	}

	/**
	 * Gets the short name.
	 *
	 * @param unit
	 *            the unit
	 * @return the short name
	 */
	public static String getShortName(DistanceUnit unit)
	{
		switch (unit)
		{
			case NM:
				return "nm";
			case PIXEL:
				return "px";
			case UM:
				return "um";
			default:
				return "na";
		}
	}

	/**
	 * Gets the short name.
	 *
	 * @param unit
	 *            the unit
	 * @return the short name
	 */
	public static String getShortName(IntensityUnit unit)
	{
		switch (unit)
		{
			case COUNT:
				return "count";
			case PHOTON:
				return "photon";
			default:
				return "na";
		}
	}

	/**
	 * Gets the short name.
	 *
	 * @param unit
	 *            the unit
	 * @return the short name
	 */
	public static String getShortName(AngleUnit unit)
	{
		switch (unit)
		{
			case DEGREE:
				return "°";
			case RADIAN:
				return "rad";
			default:
				return "na";
		}
	}
}
