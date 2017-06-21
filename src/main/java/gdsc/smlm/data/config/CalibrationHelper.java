package gdsc.smlm.data.config;

import java.util.ArrayList;

import gdsc.core.data.utils.ConversionException;
import gdsc.core.data.utils.IdentityTypeConverter;
import gdsc.core.data.utils.TypeConverter;
import gdsc.smlm.data.config.SMLMSettings.AngleUnit;
import gdsc.smlm.data.config.SMLMSettings.Calibration;
import gdsc.smlm.data.config.SMLMSettings.DistanceCalibration;
import gdsc.smlm.data.config.SMLMSettings.DistanceUnit;
import gdsc.smlm.data.config.SMLMSettings.IntensityCalibration;
import gdsc.smlm.data.config.SMLMSettings.IntensityUnit;
import gdsc.smlm.data.config.SMLMSettings.PSFCalibration;

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
	private final Calibration.Builder calibrationBuilder;

	/**
	 * Instantiates a new calibration helper.
	 *
	 * @param calibration
	 *            the calibration
	 * @throws IllegalArgumentException
	 *             if the calibration is null
	 */
	public CalibrationHelper(Calibration calibration) throws IllegalArgumentException
	{
		if (calibration == null)
			throw new IllegalArgumentException("Calibration is null");
		calibrationBuilder = calibration.toBuilder();
	}

	/**
	 * Gets the calibration.
	 *
	 * @return the calibration
	 */
	public Calibration getCalibration()
	{
		return calibrationBuilder.build();
	}

	/**
	 * Gets a distance converter to update values.
	 * <p>
	 * If the conversion is not possible then an exception is thrown.
	 * <p>
	 * If a converter if successfully created then the instance calibration will be updated.
	 *
	 * @param toDistanceUnit
	 *            the distance unit
	 * @return the distance converter
	 * @throws ConversionException
	 *             if a converter cannot be created
	 */
	public TypeConverter<DistanceUnit> getDistanceConverter(DistanceUnit toDistanceUnit)
	{
		if (toDistanceUnit != null && calibrationBuilder.hasDistanceCalibration())
		{
			DistanceCalibration.Builder distanceCalibration = calibrationBuilder.getDistanceCalibrationBuilder();
			TypeConverter<DistanceUnit> c = UnitConverterFactory.createConverter(distanceCalibration.getUnit(),
					toDistanceUnit, distanceCalibration.getNmPerPixel());
			distanceCalibration.setUnit(toDistanceUnit);
			return c;
		}
		throw new ConversionException();
	}

	/**
	 * Gets an intensity converters to update values.
	 * <p>
	 * If the conversion is not possible then an exception is thrown.
	 * <p>
	 * If a converter if successfully created then the instance calibration will be updated.
	 *
	 * @param toIntensityUnit
	 *            the intensity unit
	 * @return the intensity converter
	 * @throws ConversionException
	 *             if a converter cannot be created
	 */
	public TypeConverter<IntensityUnit> getIntensityConverter(IntensityUnit toIntensityUnit)
	{
		if (toIntensityUnit != null && calibrationBuilder.hasIntensityCalibration())
		{
			IntensityCalibration.Builder intensityCalibration = calibrationBuilder.getIntensityCalibrationBuilder();
			TypeConverter<IntensityUnit> c = UnitConverterFactory.createConverter(intensityCalibration.getUnit(),
					toIntensityUnit, intensityCalibration.getGain());
			intensityCalibration.setUnit(toIntensityUnit);
			return c;
		}
		throw new ConversionException();
	}

	/**
	 * Gets a angle converter to update values.
	 * <p>
	 * If the conversion is not possible then an exception is thrown.
	 * <p>
	 * If a converter if successfully created then the instance calibration will be updated.
	 *
	 * @param toAngleUnit
	 *            the angle unit
	 * @return the angle converter
	 * @throws ConversionException
	 *             if a converter cannot be created
	 */
	public TypeConverter<AngleUnit> getAngleConverter(AngleUnit toAngleUnit)
	{
		if (toAngleUnit != null && calibrationBuilder.hasPsfCalibration())
		{
			PSFCalibration.Builder psfCalibration = calibrationBuilder.getPsfCalibrationBuilder();
			TypeConverter<AngleUnit> c = UnitConverterFactory.createConverter(psfCalibration.getAngleUnit(),
					toAngleUnit);
			// Update the units
			psfCalibration.setAngleUnit(toAngleUnit);
			return c;
		}
		throw new ConversionException();
	}

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
	public static TypeConverter<DistanceUnit> getDistanceConverter(Calibration calibration, DistanceUnit toDistanceUnit)
	{
		if (toDistanceUnit != null && calibration.hasDistanceCalibration())
		{
			DistanceCalibration distanceCalibration = calibration.getDistanceCalibration();
			return UnitConverterFactory.createConverter(distanceCalibration.getUnit(), toDistanceUnit,
					distanceCalibration.getNmPerPixel());
		}
		throw new ConversionException();
	}

	/**
	 * Gets an intensity converters to update values.
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
	public static TypeConverter<IntensityUnit> getIntensityConverter(Calibration calibration,
			IntensityUnit toIntensityUnit)
	{
		if (toIntensityUnit != null && calibration.hasIntensityCalibration())
		{
			IntensityCalibration intensityCalibration = calibration.getIntensityCalibration();
			return UnitConverterFactory.createConverter(intensityCalibration.getUnit(), toIntensityUnit,
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
	public static TypeConverter<AngleUnit> getAngleConverter(Calibration calibration, AngleUnit toAngleUnit)
	{
		if (toAngleUnit != null && calibration.hasPsfCalibration())
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
	 * <p>
	 * If a converter if successfully created then the instance calibration will be updated.
	 *
	 * @param toDistanceUnit
	 *            the distance unit
	 * @return the distance converter
	 */
	public TypeConverter<DistanceUnit> getDistanceConverterSafe(DistanceUnit toDistanceUnit)
	{
		try
		{
			return getDistanceConverter(toDistanceUnit);
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
	 * <p>
	 * If a converter if successfully created then the instance calibration will be updated.
	 *
	 * @param toIntensityUnit
	 *            the intensity unit
	 * @return the intensity converter
	 */
	public TypeConverter<IntensityUnit> getIntensityConverterSafe(IntensityUnit toIntensityUnit)
	{
		try
		{
			return getIntensityConverter(toIntensityUnit);
		}
		catch (ConversionException e)
		{
			return new IdentityTypeConverter<IntensityUnit>(null);
		}
	}

	/**
	 * Gets a angle converter to update values.
	 * <p>
	 * If the calibration is already in the given units or conversion is not possible
	 * then an identity converter will be returned.
	 * <p>
	 * If a converter if successfully created then the instance calibration will be updated.
	 *
	 * @param toAngleUnit
	 *            the angle unit
	 * @return the angle converter
	 */
	public TypeConverter<AngleUnit> getAngleConverterSafe(AngleUnit toAngleUnit)
	{
		try
		{
			return getAngleConverter(toAngleUnit);
		}
		catch (ConversionException e)
		{
			return new IdentityTypeConverter<AngleUnit>(null);
		}
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
	public static TypeConverter<DistanceUnit> getDistanceConverterSafe(Calibration calibration,
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
	public static TypeConverter<IntensityUnit> getIntensityConverterSafe(Calibration calibration,
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
	 * Gets a angle converter to update values.
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
	public static TypeConverter<AngleUnit> getAngleConverterSafe(Calibration calibration, AngleUnit toAngleUnit)
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

}
