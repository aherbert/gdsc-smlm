package gdsc.smlm.data.config;

import gdsc.core.utils.Maths;
import gdsc.smlm.data.config.SMLMSettings.AngleUnit;
import gdsc.smlm.data.config.SMLMSettings.DistanceUnit;
import gdsc.smlm.data.config.SMLMSettings.IntensityUnit;
import gdsc.smlm.data.config.SMLMSettings.TimeUnit;
import gdsc.smlm.data.utils.AddMultiplyTypeConverter;
import gdsc.smlm.data.utils.ConversionException;
import gdsc.smlm.data.utils.IdentityTypeConverter;
import gdsc.smlm.data.utils.MultiplyAddTypeConverter;
import gdsc.smlm.data.utils.MultiplyTypeConverter;
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
 * Factory for creating unit converters
 */
public class UnitConverterFactory
{
	/**
	 * Creates the converter.
	 *
	 * @param from
	 *            the from
	 * @param to
	 *            the to
	 * @param nmPerPixel
	 *            the nm per pixel
	 * @return the unit converter
	 * @throws ConversionException
	 *             if a converter cannot be created
	 */
	public static TypeConverter<DistanceUnit> createConverter(DistanceUnit from, DistanceUnit to, double nmPerPixel)
			throws ConversionException
	{
		if (from == to)
			return new IdentityTypeConverter<DistanceUnit>(from);
		if (from == null)
			throw new ConversionException("from must not be null");
		if (to == null)
			throw new ConversionException("to must not be null");

		switch (from)
		{
			case NM:
				switch (to)
				{
					case PIXEL:
						return new MultiplyTypeConverter<DistanceUnit>(from, to, 1.0 / checkNmPerPixel(nmPerPixel));
					case UM:
						return new MultiplyTypeConverter<DistanceUnit>(from, to, 1e-3);
					default:
						throw new ConversionException(from + " to " + to);
				}

			case PIXEL:
				switch (to)
				{
					case NM:
						return new MultiplyTypeConverter<DistanceUnit>(from, to, checkNmPerPixel(nmPerPixel));
					case UM:
						return new MultiplyTypeConverter<DistanceUnit>(from, to, checkNmPerPixel(nmPerPixel) / 1e3);
					default:
						throw new ConversionException(from + " to " + to);
				}

			case UM:
				switch (to)
				{
					case NM:
						return new MultiplyTypeConverter<DistanceUnit>(from, to, 1e3);
					case PIXEL:
						return new MultiplyTypeConverter<DistanceUnit>(from, to, 1e3 / checkNmPerPixel(nmPerPixel));
					default:
						throw new ConversionException(from + " to " + to);
				}

			default:
				throw new ConversionException(from + " to " + to);
		}
	}

	private static double checkNmPerPixel(double nmPerPixel)
	{
		if (!(nmPerPixel > 0 && nmPerPixel <= java.lang.Double.MAX_VALUE))
			throw new ConversionException("nm/pixel must be positive");
		return nmPerPixel;
	}

	/**
	 * Creates the converter.
	 *
	 * @param from
	 *            the from
	 * @param to
	 *            the to
	 * @param msPerFrame
	 *            the ms per frame
	 * @return the unit converter
	 * @throws ConversionException
	 *             if a converter cannot be created
	 */
	public static TypeConverter<TimeUnit> createConverter(TimeUnit from, TimeUnit to, double msPerFrame)
			throws ConversionException
	{
		if (from == to)
			return new IdentityTypeConverter<TimeUnit>(from);
		if (from == null)
			throw new ConversionException("from must not be null");
		if (to == null)
			throw new ConversionException("to must not be null");

		switch (from)
		{
			case MILLISECOND:
				switch (to)
				{
					case FRAME:
						return new MultiplyTypeConverter<TimeUnit>(from, to, 1.0 / checkMsPerFrame(msPerFrame));
					case SECOND:
						return new MultiplyTypeConverter<TimeUnit>(from, to, 1e-3);
					default:
						throw new ConversionException(from + " to " + to);
				}

			case FRAME:
				switch (to)
				{
					case MILLISECOND:
						return new MultiplyTypeConverter<TimeUnit>(from, to, checkMsPerFrame(msPerFrame));
					case SECOND:
						return new MultiplyTypeConverter<TimeUnit>(from, to, checkMsPerFrame(msPerFrame) / 1e3);
					default:
						throw new ConversionException(from + " to " + to);
				}

			case SECOND:
				switch (to)
				{
					case MILLISECOND:
						return new MultiplyTypeConverter<TimeUnit>(from, to, 1e3);
					case FRAME:
						return new MultiplyTypeConverter<TimeUnit>(from, to, 1e3 / checkMsPerFrame(msPerFrame));
					default:
						throw new ConversionException(from + " to " + to);
				}

			default:
				throw new ConversionException(from + " to " + to);
		}
	}

	private static double checkMsPerFrame(double msPerFrame)
	{
		if (!(msPerFrame > 0 && msPerFrame <= java.lang.Double.MAX_VALUE))
			throw new ConversionException("ms/frame must be positive");
		return msPerFrame;
	}

	/**
	 * Creates the converter.
	 *
	 * @param from
	 *            the from
	 * @param to
	 *            the to
	 * @return the unit converter
	 * @throws ConversionException
	 *             if a converter cannot be created
	 */
	public static TypeConverter<AngleUnit> createConverter(AngleUnit from, AngleUnit to) throws ConversionException
	{
		if (from == to)
			return new IdentityTypeConverter<AngleUnit>(from);
		if (from == null)
			throw new ConversionException("from must not be null");
		if (to == null)
			throw new ConversionException("to must not be null");

		switch (from)
		{
			case RADIAN:
				switch (to)
				{
					case DEGREE:
						return new MultiplyTypeConverter<AngleUnit>(from, to, 180.0 / Math.PI);
					default:
						throw new ConversionException(from + " to " + to);
				}

			case DEGREE:
				switch (to)
				{
					case RADIAN:
						return new MultiplyTypeConverter<AngleUnit>(from, to, Math.PI / 180.0);
					default:
						throw new ConversionException(from + " to " + to);
				}

			default:
				throw new ConversionException(from + " to " + to);
		}
	}

	/**
	 * Creates the converter.
	 *
	 * @param from
	 *            the from
	 * @param to
	 *            the to
	 * @param offset
	 *            the offset
	 * @param countPerPhoton
	 *            the count per photon
	 * @return the unit converter
	 * @throws ConversionException
	 *             if a converter cannot be created
	 */
	public static TypeConverter<IntensityUnit> createConverter(IntensityUnit from, IntensityUnit to, double offset,
			double countPerPhoton) throws ConversionException
	{
		if (from == to)
			return new IdentityTypeConverter<IntensityUnit>(from);
		if (from == null)
			throw new ConversionException("from must not be null");
		if (to == null)
			throw new ConversionException("to must not be null");

		switch (from)
		{
			case COUNT:
				switch (to)
				{
					case PHOTON:
						return new AddMultiplyTypeConverter<IntensityUnit>(from, to, checkOffset(-offset),
								1.0 / checkCountPerPhoton(countPerPhoton));
					default:
						throw new ConversionException(from + " to " + to);
				}

			case PHOTON:
				switch (to)
				{
					case COUNT:
						return new MultiplyAddTypeConverter<IntensityUnit>(from, to,
								checkCountPerPhoton(countPerPhoton), checkOffset(offset));
					default:
						throw new ConversionException(from + " to " + to);
				}

			default:
				throw new ConversionException(from + " to " + to);
		}
	}

	private static double checkOffset(double offset)
	{
		if (!Maths.isFinite(offset))
			throw new ConversionException("offset must be finite");
		return offset;
	}

	private static double checkCountPerPhoton(double countPerPhoton)
	{
		if (!(countPerPhoton > 0 && countPerPhoton <= java.lang.Double.MAX_VALUE))
			throw new ConversionException("count/photon must be positive");
		return countPerPhoton;
	}

	/**
	 * Creates the converter.
	 *
	 * @param from
	 *            the from
	 * @param to
	 *            the to
	 * @param countPerPhoton
	 *            the count per photon
	 * @return the unit converter
	 * @throws ConversionException
	 *             if a converter cannot be created
	 */
	public static TypeConverter<IntensityUnit> createConverter(IntensityUnit from, IntensityUnit to,
			double countPerPhoton) throws ConversionException
	{
		if (from == to)
			return new IdentityTypeConverter<IntensityUnit>(from);
		if (from == null)
			throw new ConversionException("from must not be null");
		if (to == null)
			throw new ConversionException("to must not be null");

		switch (from)
		{
			case COUNT:
				switch (to)
				{
					case PHOTON:
						return new MultiplyTypeConverter<IntensityUnit>(from, to,
								1.0 / checkCountPerPhoton(countPerPhoton));
					default:
						throw new ConversionException(from + " to " + to);
				}

			case PHOTON:
				switch (to)
				{
					case COUNT:
						return new MultiplyTypeConverter<IntensityUnit>(from, to, checkCountPerPhoton(countPerPhoton));
					default:
						throw new ConversionException(from + " to " + to);
				}

			default:
				throw new ConversionException(from + " to " + to);
		}
	}
}
