package gdsc.smlm.data.config;

import gdsc.core.data.utils.ConversionException;
import gdsc.core.data.utils.IdentityTypeConverter;
import gdsc.core.data.utils.TypeConverter;
import gdsc.smlm.data.config.UnitConfig.AngleUnit;
import gdsc.smlm.data.config.CalibrationConfig.Calibration;
import gdsc.smlm.data.config.CalibrationConfig.CalibrationOrBuilder;
import gdsc.smlm.data.config.CalibrationConfig.DistanceCalibration;
import gdsc.smlm.data.config.UnitConfig.DistanceUnit;
import gdsc.smlm.data.config.CalibrationConfig.IntensityCalibration;
import gdsc.smlm.data.config.UnitConfig.IntensityUnit;
import gdsc.smlm.data.config.CalibrationConfig.PSFCalibration;

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
 * Contains helper functions for the Calibration class.
 */
public class CalibrationHelper
{
	/**
	 * Gets a distance converter to update values.
	 * <p>
	 * If the conversion is not possible then an exception is thrown.
	 *
	 * @param calibration
	 *            the calibration
	 * @param toDistanceUnit
	 *            the distance unit
	 * @return the distance converter
	 * @throws ConversionException
	 *             if a converter cannot be created
	 */
	public static TypeConverter<DistanceUnit> getDistanceConverter(CalibrationOrBuilder calibration,
			DistanceUnit toDistanceUnit) throws ConversionException
	{
		if (calibration != null && toDistanceUnit != null && calibration.hasDistanceCalibration())
		{
			DistanceCalibration distanceCalibration = calibration.getDistanceCalibration();
			return UnitConverterFactory.createConverter(distanceCalibration.getDistanceUnit(), toDistanceUnit,
					distanceCalibration.getNmPerPixel());
		}
		throw new ConversionException();
	}

	/**
	 * Gets a intensity converter to update values.
	 * <p>
	 * If the conversion is not possible then an exception is thrown.
	 *
	 * @param calibration
	 *            the calibration
	 * @param toIntensityUnit
	 *            the intensity unit
	 * @return the intensity converter
	 * @throws ConversionException
	 *             if a converter cannot be created
	 */
	public static TypeConverter<IntensityUnit> getIntensityConverter(CalibrationOrBuilder calibration,
			IntensityUnit toIntensityUnit) throws ConversionException
	{
		if (calibration != null && toIntensityUnit != null && calibration.hasIntensityCalibration())
		{
			IntensityCalibration intensityCalibration = calibration.getIntensityCalibration();
			return UnitConverterFactory.createConverter(intensityCalibration.getIntensityUnit(), toIntensityUnit,
					intensityCalibration.getGain());
		}
		throw new ConversionException();
	}

	/**
	 * Gets an angle converter to update values.
	 * <p>
	 * If the conversion is not possible then an exception is thrown.
	 *
	 * @param calibration
	 *            the calibration
	 * @param toAngleUnit
	 *            the angle unit
	 * @return the angle converter
	 * @throws ConversionException
	 *             if a converter cannot be created
	 */
	public static TypeConverter<AngleUnit> getAngleConverter(CalibrationOrBuilder calibration, AngleUnit toAngleUnit)
			throws ConversionException
	{
		if (calibration != null && toAngleUnit != null && calibration.hasPsfCalibration())
		{
			PSFCalibration psfCalibration = calibration.getPsfCalibration();
			return UnitConverterFactory.createConverter(psfCalibration.getAngleUnit(), toAngleUnit);
		}
		throw new ConversionException();
	}

	/**
	 * Gets a distance converter to update values.
	 * <p>
	 * If the calibration is already in the given units or conversion is not possible
	 * then an identity converter will be returned.
	 *
	 * @param calibration
	 *            the calibration
	 * @param toDistanceUnit
	 *            the distance unit
	 * @return the distance converter
	 */
	public static TypeConverter<DistanceUnit> getDistanceConverterSafe(CalibrationOrBuilder calibration,
			DistanceUnit toDistanceUnit)
	{
		try
		{
			return getDistanceConverter(calibration, toDistanceUnit);
		}
		catch (ConversionException e)
		{
			return new IdentityTypeConverter<DistanceUnit>(null);
		}
	}

	/**
	 * Gets an intensity converter to update values.
	 * <p>
	 * If the calibration is already in the given units or conversion is not possible
	 * then an identity converter will be returned.
	 *
	 * @param calibration
	 *            the calibration
	 * @param toIntensityUnit
	 *            the intensity unit
	 * @return the intensity converter
	 */
	public static TypeConverter<IntensityUnit> getIntensityConverterSafe(CalibrationOrBuilder calibration,
			IntensityUnit toIntensityUnit)
	{
		try
		{
			return getIntensityConverter(calibration, toIntensityUnit);
		}
		catch (ConversionException e)
		{
			return new IdentityTypeConverter<IntensityUnit>(null);
		}
	}

	/**
	 * Gets an angle converter to update values.
	 * <p>
	 * If the calibration is already in the given units or conversion is not possible
	 * then an identity converter will be returned.
	 *
	 * @param calibration
	 *            the calibration
	 * @param toAngleUnit
	 *            the angle unit
	 * @return the angle converter
	 */
	public static TypeConverter<AngleUnit> getAngleConverterSafe(CalibrationOrBuilder calibration,
			AngleUnit toAngleUnit)
	{
		try
		{
			return getAngleConverter(calibration, toAngleUnit);
		}
		catch (ConversionException e)
		{
			return new IdentityTypeConverter<AngleUnit>(null);
		}
	}

	/**
	 * Create a new calibration using the given properties.
	 *
	 * @param nmPerPixel
	 *            the nm per pixel
	 * @param gain
	 *            the gain
	 * @param exposureTime
	 *            the exposure time
	 * @return the calibration
	 */
	public static Calibration create(double nmPerPixel, double gain, double exposureTime)
	{
		CalibrationWriter cw = new CalibrationWriter();
		if (nmPerPixel > 0)
			cw.setNmPerPixel(nmPerPixel);
		if (gain > 0)
			cw.setGain(gain);
		if (exposureTime > 0)
			cw.setExposureTime(exposureTime);
		return cw.getCalibration();
	}
}
