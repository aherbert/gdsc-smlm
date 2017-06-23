package gdsc.smlm.data.config;

import gdsc.smlm.data.config.UnitConfig.AngleUnit;
import gdsc.smlm.data.config.UnitConfig.DistanceUnit;
import gdsc.smlm.data.config.UnitConfig.IntensityUnit;
import gdsc.smlm.data.config.UnitConfig.TimeUnit;

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
 * Contains helper functions for the units. Adds functionality to the protocol buffer unit enums.
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
	 * Guess distance unit from short name.
	 *
	 * @param name
	 *            the name
	 * @return the distance unit
	 */
	public static DistanceUnit guessDistanceUnitFromShortName(String name)
	{
		if (name != null)
		{
			if (name.equalsIgnoreCase("px"))
				return DistanceUnit.PIXEL;
			if (name.equalsIgnoreCase("um"))
				return DistanceUnit.UM;
			if (name.equalsIgnoreCase("nm"))
				return DistanceUnit.NM;
		}
		return null;
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
	 * Guess intensity unit from short name.
	 *
	 * @param name
	 *            the name
	 * @return the intensity unit
	 */
	public static IntensityUnit guessIntensityUnitFromShortName(String name)
	{
		if (name != null)
		{
			if (name.equalsIgnoreCase("photon"))
				return IntensityUnit.PHOTON;
			if (name.equalsIgnoreCase("count"))
				return IntensityUnit.COUNT;
		}
		return null;
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

	/**
	 * Guess angle unit from short name.
	 *
	 * @param name
	 *            the name
	 * @return the angle unit
	 */
	public static AngleUnit guessAngleUnitFromShortName(String name)
	{
		if (name != null)
		{
			if (name.equalsIgnoreCase("°"))
				return AngleUnit.DEGREE;
			if (name.equalsIgnoreCase("rad"))
				return AngleUnit.RADIAN;
		}
		return null;
	}

	/**
	 * Gets the name.
	 *
	 * @param unit
	 *            the unit
	 * @return the name
	 */
	public static String getName(TimeUnit unit)
	{
		switch (unit)
		{
			case FRAME:
				return "frame";
			case MILLISECOND:
				return "millisecond";
			case SECOND:
				return "second";
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
	public static String getShortName(TimeUnit unit)
	{
		switch (unit)
		{
			case FRAME:
				return "frame";
			case MILLISECOND:
				return "ms";
			case SECOND:
				return "s";
			default:
				return "na";
		}
	}

	/**
	 * Guess time unit from short name.
	 *
	 * @param name
	 *            the name
	 * @return the time unit
	 */
	public static TimeUnit guessTimeUnitFromShortName(String name)
	{
		if (name != null)
		{
			if (name.equalsIgnoreCase("frame"))
				return TimeUnit.FRAME;
			if (name.equalsIgnoreCase("s"))
				return TimeUnit.SECOND;
			if (name.equalsIgnoreCase("ms"))
				return TimeUnit.MILLISECOND;
		}
		return null;
	}
}
