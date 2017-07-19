package gdsc.smlm.data.config;

import gdsc.core.data.utils.ConversionException;
import gdsc.core.data.utils.IdentityTypeConverter;
import gdsc.core.data.utils.TypeConverter;
import gdsc.smlm.data.config.CalibrationProtos.AngleCalibrationOrBuilder;
import gdsc.smlm.data.config.CalibrationProtos.Calibration;
import gdsc.smlm.data.config.CalibrationProtos.CalibrationOrBuilder;
import gdsc.smlm.data.config.CalibrationProtos.DistanceCalibrationOrBuilder;
import gdsc.smlm.data.config.CalibrationProtos.IntensityCalibrationOrBuilder;
import gdsc.smlm.data.config.UnitProtos.AngleUnit;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.data.config.UnitProtos.IntensityUnit;

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
			DistanceCalibrationOrBuilder distanceCalibration = calibration.getDistanceCalibrationOrBuilder();
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
			IntensityCalibrationOrBuilder intensityCalibration = calibration.getIntensityCalibrationOrBuilder();
			return UnitConverterFactory.createConverter(intensityCalibration.getIntensityUnit(), toIntensityUnit,
					intensityCalibration.getCountPerPhoton());
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
		if (calibration != null && toAngleUnit != null && calibration.hasAngleCalibration())
		{
			AngleCalibrationOrBuilder psfCalibration = calibration.getAngleCalibrationOrBuilder();
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
			if (calibration != null && calibration.hasDistanceCalibration())
				return new IdentityTypeConverter<DistanceUnit>(
						calibration.getDistanceCalibrationOrBuilder().getDistanceUnit());
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
			if (calibration != null && calibration.hasIntensityCalibration())
				return new IdentityTypeConverter<IntensityUnit>(
						calibration.getIntensityCalibrationOrBuilder().getIntensityUnit());
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
			if (calibration != null && calibration.hasAngleCalibration())
				return new IdentityTypeConverter<AngleUnit>(calibration.getAngleCalibrationOrBuilder().getAngleUnit());
			return new IdentityTypeConverter<AngleUnit>(null);
		}
	}

	/**
	 * Create a new calibration using the given properties. Distance unit is set to pixel and intensity unit to photon.
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
		cw.setIntensityUnit(IntensityUnit.PHOTON);
		cw.setDistanceUnit(DistanceUnit.PIXEL);
		if (nmPerPixel > 0)
			cw.setNmPerPixel(nmPerPixel);
		if (gain > 0)
			cw.setGain(gain);
		if (exposureTime > 0)
			cw.setExposureTime(exposureTime);
		return cw.getCalibration();
	}
}
