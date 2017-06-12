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
 * Contains helper functions for the Calibration class.
 */
public class CalibrationHelper
{
	// TODO - Create an instance version that can build an updated calibration 
	// if the converter was successfully created.
	
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
	public static TypeConverter<DistanceUnit> getDistanceConverter(Calibration calibration, DistanceUnit toDistanceUnit)
	{
		TypeConverter<DistanceUnit> c = null;
		if (toDistanceUnit != null && calibration.hasDistanceCalibration())
		{
			DistanceCalibration distanceCalibration = calibration.getDistanceCalibration();
			try
			{
				c = UnitConverterFactory.createConverter(distanceCalibration.getUnit(), toDistanceUnit,
						distanceCalibration.getNmPerPixel());
			}
			catch (ConversionException e)
			{
				// Ignore this
			}
		}
		if (c == null)
			c = new IdentityTypeConverter<DistanceUnit>(null);
		return c;
	}
	
	// Fix these

//	/**
//	 * Gets intensity converters to update values. If the converter can be created then the current intensity unit in
//	 * this instance is updated.
//	 * <p>
//	 * If the calibration is already in the given units or conversion is not possible
//	 * then an identity converter will be returned.
//	 * <p>
//	 * The returned list calibration.has a converter with only the gain, and a second converter with the gain and bias.
//	 * <p>
//	 * It is recommended to clone the calibration before invoking this method as the state may change.
//	 *
//	 * @param calibration
//	 *            the calibration
//	 * @param toIntensityUnit
//	 *            the intensity unit
//	 * @return the intensity converters (gain, gain + bias)
//	 */
//	public ArrayList<TypeConverter<IntensityUnit>> getIntensityConverter(Calibration calibration,
//			IntensityUnit toIntensityUnit)
//	{
//		ArrayList<TypeConverter<IntensityUnit>> list = new ArrayList<TypeConverter<IntensityUnit>>(2);
//		if (toIntensityUnit != null && calibration.hasIntensityUnit() && intensityUnit != toIntensityUnit &&
//				calibration.hasGain() && calibration.hasBias())
//		{
//			try
//			{
//				list.add(UnitConverterFactory.createConverter(intensityUnit, toIntensityUnit, gain));
//				list.add(UnitConverterFactory.createConverter(intensityUnit, toIntensityUnit, bias, gain));
//				setIntensityUnit(toIntensityUnit);
//			}
//			catch (ConversionException e)
//			{
//				// Ignore this
//			}
//		}
//		if (list.size() != 2)
//		{
//			list.clear();
//			TypeConverter<IntensityUnit> c = new IdentityTypeConverter<IntensityUnit>(intensityUnit);
//			list.add(c);
//			list.add(c);
//		}
//		return list;
//	}
//
//	/**
//	 * Gets a angle converter to update values. If the converter can be created then the current angle unit in
//	 * this instance is updated.
//	 * <p>
//	 * If the calibration is already in the given units or conversion is not possible
//	 * then an identity converter will be returned.
//	 * <p>
//	 * It is recommended to clone the calibration before invoking this method as the state may change.
//	 *
//	 * @param calibration
//	 *            the calibration
//	 * @param toAngleUnit
//	 *            the angle unit
//	 * @return the angle converter
//	 */
//	public TypeConverter<AngleUnit> getAngleConverter(Calibration calibration, AngleUnit toAngleUnit)
//	{
//		TypeConverter<AngleUnit> c = null;
//		if (toAngleUnit != null && calibration.hasAngleUnit() && angleUnit != toAngleUnit)
//		{
//			try
//			{
//				c = UnitConverterFactory.createConverter(angleUnit, toAngleUnit);
//				setAngleUnit(toAngleUnit);
//			}
//			catch (ConversionException e)
//			{
//				// Ignore this
//			}
//		}
//		if (c == null)
//			c = new IdentityTypeConverter<AngleUnit>(angleUnit);
//		return c;
//	}
}
