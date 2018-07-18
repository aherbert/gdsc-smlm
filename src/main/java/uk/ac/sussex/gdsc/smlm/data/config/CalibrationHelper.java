/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package uk.ac.sussex.gdsc.smlm.data.config;

import uk.ac.sussex.gdsc.core.data.utils.ConversionException;
import uk.ac.sussex.gdsc.core.data.utils.IdentityTypeConverter;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.AngleCalibrationOrBuilder;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CalibrationOrBuilder;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.DistanceCalibrationOrBuilder;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.IntensityCalibrationOrBuilder;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.TimeCalibrationOrBuilder;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.TimeUnit;

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
			final DistanceCalibrationOrBuilder distanceCalibration = calibration.getDistanceCalibrationOrBuilder();
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
			final IntensityCalibrationOrBuilder intensityCalibration = calibration.getIntensityCalibrationOrBuilder();
			return UnitConverterFactory.createConverter(intensityCalibration.getIntensityUnit(), toIntensityUnit,
					intensityCalibration.getCountPerPhoton());
		}
		throw new ConversionException();
	}

	/**
	 * Gets a time converter to update values.
	 * <p>
	 * If the conversion is not possible then an exception is thrown.
	 *
	 * @param calibration
	 *            the calibration
	 * @param toTimeUnit
	 *            the time unit
	 * @return the time converter
	 * @throws ConversionException
	 *             if a converter cannot be created
	 */
	public static TypeConverter<TimeUnit> getTimeConverter(CalibrationOrBuilder calibration, TimeUnit toTimeUnit)
			throws ConversionException
	{
		if (calibration != null && toTimeUnit != null && calibration.hasTimeCalibration())
		{
			final TimeCalibrationOrBuilder timeCalibration = calibration.getTimeCalibrationOrBuilder();
			// Assume time is in frames
			final TimeUnit timeUnit = TimeUnit.FRAME; // timeCalibration.getTimeUnit()
			return UnitConverterFactory.createConverter(timeUnit, toTimeUnit, timeCalibration.getExposureTime());
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
			final AngleCalibrationOrBuilder psfCalibration = calibration.getAngleCalibrationOrBuilder();
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
		catch (final ConversionException e)
		{
			if (calibration != null && calibration.hasDistanceCalibration())
				return new IdentityTypeConverter<>(
						calibration.getDistanceCalibrationOrBuilder().getDistanceUnit());
			return new IdentityTypeConverter<>(null);
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
		catch (final ConversionException e)
		{
			if (calibration != null && calibration.hasIntensityCalibration())
				return new IdentityTypeConverter<>(
						calibration.getIntensityCalibrationOrBuilder().getIntensityUnit());
			return new IdentityTypeConverter<>(null);
		}
	}

	/**
	 * Gets an time converter to update values.
	 * <p>
	 * If the calibration is already in the given units or conversion is not possible
	 * then an identity converter will be returned.
	 *
	 * @param calibration
	 *            the calibration
	 * @param toTimeUnit
	 *            the time unit
	 * @return the time converter
	 */
	public static TypeConverter<TimeUnit> getTimeConverterSafe(CalibrationOrBuilder calibration, TimeUnit toTimeUnit)
	{
		try
		{
			return getTimeConverter(calibration, toTimeUnit);
		}
		catch (final ConversionException e)
		{
			// Calibration is assumed to be in frames
			return new IdentityTypeConverter<>(TimeUnit.FRAME);
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
		catch (final ConversionException e)
		{
			if (calibration != null && calibration.hasAngleCalibration())
				return new IdentityTypeConverter<>(calibration.getAngleCalibrationOrBuilder().getAngleUnit());
			return new IdentityTypeConverter<>(null);
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
		final CalibrationWriter cw = new CalibrationWriter();
		cw.setIntensityUnit(IntensityUnit.PHOTON);
		cw.setDistanceUnit(DistanceUnit.PIXEL);
		if (nmPerPixel > 0)
			cw.setNmPerPixel(nmPerPixel);
		if (gain > 0)
			cw.setCountPerPhoton(gain);
		if (exposureTime > 0)
			cw.setExposureTime(exposureTime);
		return cw.getCalibration();
	}
}
