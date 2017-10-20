/*
 * 
 */
package gdsc.smlm.data.config;

import gdsc.core.data.utils.ConversionException;
import gdsc.core.data.utils.TypeConverter;
import gdsc.core.utils.Maths;
import gdsc.core.utils.TextUtils;
import gdsc.smlm.data.config.UnitProtos.AngleUnit;
import gdsc.smlm.data.config.CalibrationProtos.CalibrationOrBuilder;
import gdsc.smlm.data.config.CalibrationProtos.CameraType;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import gdsc.smlm.data.config.UnitProtos.TimeUnit;

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
 * Contains helper functions for reading the Calibration class.
 */
public class CalibrationReader
{
	/** The calibration or builder. */
	private CalibrationOrBuilder calibrationOrBuilder;

	/**
	 * Instantiates a new calibration reader with no calibration.
	 */
	protected CalibrationReader()
	{
	}

	/**
	 * Instantiates a new calibration reader.
	 *
	 * @param calibration
	 *            the calibration
	 * @throws IllegalArgumentException
	 *             if the calibration is null
	 */
	public CalibrationReader(CalibrationOrBuilder calibration) throws IllegalArgumentException
	{
		if (calibration == null)
			throw new IllegalArgumentException("Calibration is null");
		this.calibrationOrBuilder = calibration;
	}

	/**
	 * Gets the calibration or builder with the latest changes.
	 *
	 * @return the calibration or builder
	 */
	public CalibrationOrBuilder getCalibrationOrBuilder()
	{
		return calibrationOrBuilder;
	}

	/**
	 * Gets a distance converter to update values.
	 * <p>
	 * If the conversion is not possible then an exception is thrown.
	 *
	 * @param toDistanceUnit
	 *            the distance unit
	 * @return the distance converter
	 * @throws ConversionException
	 *             the conversion exception
	 */
	public TypeConverter<DistanceUnit> getDistanceConverter(DistanceUnit toDistanceUnit) throws ConversionException
	{
		return CalibrationHelper.getDistanceConverter(getCalibrationOrBuilder(), toDistanceUnit);
	}

	/**
	 * Gets an intensity converter to update values.
	 * <p>
	 * If the conversion is not possible then an exception is thrown.
	 *
	 * @param toIntensityUnit
	 *            the intensity unit
	 * @return the intensity converter
	 * @throws ConversionException
	 *             the conversion exception
	 */
	public TypeConverter<IntensityUnit> getIntensityConverter(IntensityUnit toIntensityUnit) throws ConversionException
	{
		return CalibrationHelper.getIntensityConverter(getCalibrationOrBuilder(), toIntensityUnit);
	}

	/**
	 * Gets an time converter to update values.
	 * <p>
	 * If the conversion is not possible then an exception is thrown.
	 *
	 * @param toTimeUnit
	 *            the time unit
	 * @return the time converter
	 * @throws ConversionException
	 *             the conversion exception
	 */
	public TypeConverter<TimeUnit> getTimeConverter(TimeUnit toTimeUnit) throws ConversionException
	{
		return CalibrationHelper.getTimeConverter(getCalibrationOrBuilder(), toTimeUnit);
	}

	/**
	 * Gets an angle converter to update values.
	 * <p>
	 * If the conversion is not possible then an exception is thrown.
	 *
	 * @param toAngleUnit
	 *            the angle unit
	 * @return the angle converter
	 * @throws ConversionException
	 *             the conversion exception
	 */
	public TypeConverter<AngleUnit> getAngleConverter(AngleUnit toAngleUnit) throws ConversionException
	{
		return CalibrationHelper.getAngleConverter(getCalibrationOrBuilder(), toAngleUnit);
	}

	/**
	 * Gets a distance converter to update values.
	 * <p>
	 * If the calibration is already in the given units or conversion is not possible
	 * then an identity converter will be returned.
	 *
	 * @param toDistanceUnit
	 *            the distance unit
	 * @return CalibrationHelper.the distance converter
	 */
	public TypeConverter<DistanceUnit> getDistanceConverterSafe(DistanceUnit toDistanceUnit)
	{
		return CalibrationHelper.getDistanceConverterSafe(getCalibrationOrBuilder(), toDistanceUnit);
	}

	/**
	 * Gets an intensity converter to update values.
	 * <p>
	 * If the calibration is already in the given units or conversion is not possible
	 * then an identity converter will be returned.
	 *
	 * @param toIntensityUnit
	 *            the intensity unit
	 * @return CalibrationHelper.the intensity converter
	 */
	public TypeConverter<IntensityUnit> getIntensityConverterSafe(IntensityUnit toIntensityUnit)
	{
		return CalibrationHelper.getIntensityConverterSafe(getCalibrationOrBuilder(), toIntensityUnit);
	}

	/**
	 * Gets an time converter to update values.
	 * <p>
	 * If the calibration is already in the given units or conversion is not possible
	 * then an identity converter will be returned.
	 *
	 * @param toTimeUnit
	 *            the time unit
	 * @return CalibrationHelper.the time converter
	 */
	public TypeConverter<TimeUnit> getTimeConverterSafe(TimeUnit toTimeUnit)
	{
		return CalibrationHelper.getTimeConverterSafe(getCalibrationOrBuilder(), toTimeUnit);
	}

	/**
	 * Gets an angle converter to update values.
	 * <p>
	 * If the calibration is already in the given units or conversion is not possible
	 * then an identity converter will be returned.
	 *
	 * @param toAngleUnit
	 *            the angle unit
	 * @return CalibrationHelper.the angle converter
	 */
	public TypeConverter<AngleUnit> getAngleConverterSafe(AngleUnit toAngleUnit)
	{
		return CalibrationHelper.getAngleConverterSafe(getCalibrationOrBuilder(), toAngleUnit);
	}

	/**
	 * Gets image pixel size in nanometers.
	 *
	 * @return image pixel size in nanometers
	 */
	public double getNmPerPixel()
	{
		CalibrationOrBuilder c = getCalibrationOrBuilder();
		return (c.hasDistanceCalibration()) ? c.getDistanceCalibration().getNmPerPixel() : 0;
	}

	/**
	 * Checks for nm per pixel.
	 *
	 * @return true, if successful
	 */
	public boolean hasNmPerPixel()
	{
		return getNmPerPixel() > 0;
	}

	/**
	 * Gets the gain (Count/photon). Can be used to convert the signal in count units to photons.
	 *
	 * @return the gain
	 */
	public double getCountPerPhoton()
	{
		CalibrationOrBuilder c = getCalibrationOrBuilder();
		return (c.hasIntensityCalibration()) ? c.getIntensityCalibration().getCountPerPhoton() : 0;
	}

	/**
	 * Checks for gain (Count/photon).
	 *
	 * @return true, if successful
	 */
	public boolean hasCountPerPhoton()
	{
		return getCountPerPhoton() > 0;
	}

	/**
	 * Gets the exposure time in milliseconds per frame.
	 *
	 * @return the exposure time
	 */
	public double getExposureTime()
	{
		CalibrationOrBuilder c = getCalibrationOrBuilder();
		return (c.hasTimeCalibration()) ? c.getTimeCalibration().getExposureTime() : 0;
	}

	/**
	 * Checks for exposure time.
	 *
	 * @return true, if successful
	 */
	public boolean hasExposureTime()
	{
		return getExposureTime() > 0;
	}

	/**
	 * Gets the camera Gaussian read noise (in Count units).
	 *
	 * @return the camera Gaussian read noise (in Count units)
	 */
	public double getReadNoise()
	{
		CalibrationOrBuilder c = getCalibrationOrBuilder();
		return (c.hasCameraCalibration()) ? c.getCameraCalibration().getReadNoise() : 0;
	}

	/**
	 * Checks for read noise.
	 *
	 * @return true, if successful
	 */
	public boolean hasReadNoise()
	{
		return getReadNoise() > 0;
	}

	/**
	 * Gets camera bias (in Count units).
	 *
	 * @return camera bias (in Count units)
	 */
	public double getBias()
	{
		CalibrationOrBuilder c = getCalibrationOrBuilder();
		return (c.hasCameraCalibration()) ? c.getCameraCalibration().getBias() : 0;
	}

	/**
	 * Checks for bias.
	 *
	 * @return true, if successful
	 */
	public boolean hasBias()
	{
		// Bias can be zero (the default) so check for a camera calibration first
		CalibrationOrBuilder c = getCalibrationOrBuilder();
		return c.hasCameraCalibration() && c.getCameraCalibration().getBias() >= 0;
	}

	/**
	 * Get the camera type.
	 *
	 * @return the camera type
	 */
	public CameraType getCameraType()
	{
		CalibrationOrBuilder c = getCalibrationOrBuilder();
		return (c.hasCameraCalibration()) ? c.getCameraCalibration().getCameraType() : CameraType.CAMERA_TYPE_NA;
	}

	/**
	 * Checks for camera type.
	 *
	 * @return true, if successful
	 */
	public boolean hasCameraType()
	{
		return getCameraType().getNumber() > 0;
	}

	/**
	 * Checks for camera calibration.
	 *
	 * @return true, if successful
	 */
	public boolean hasCameraCalibration()
	{
		return getCalibrationOrBuilder().hasCameraCalibration();
	}

	/**
	 * Checks for a CCD camera.
	 *
	 * @return true, if successful
	 */
	public boolean isCCDCamera()
	{
		CalibrationOrBuilder c = getCalibrationOrBuilder();
		if (c.hasCameraCalibration())
		{
			return CalibrationProtosHelper.isCCDCameraType(c.getCameraCalibration().getCameraType());
		}
		return false;
	}

	/**
	 * Checks if the camera type was an Electron Multiplying (EM) CCD.
	 *
	 * @return true, if the camera type was an Electron Multiplying (EM) CCD
	 */
	public boolean isEMCCD()
	{
		return getCameraType() == CameraType.EMCCD;
	}

	/**
	 * Checks if the camera type was a standard CCD.
	 *
	 * @return true, if the camera type was a standard CCD.
	 */
	public boolean isCCD()
	{
		return getCameraType() == CameraType.CCD;
	}

	/**
	 * Checks if the camera type was a sCMOS.
	 *
	 * @return true, if the camera type was a sCMOS.
	 */
	public boolean isSCMOS()
	{
		return getCameraType() == CameraType.SCMOS;
	}

	/**
	 * Get the camera quantum efficiency (e-/photon) used when modelling a microscope camera.
	 * <p>
	 * Note that the camera noise model assumes that photons are converted to counts by
	 * a process that is not perfect (i.e. it has noise). The underlying process is
	 * photons converted to electrons in the camera chip and then amplification
	 * (count/electron) occurring in the camera hardware. Ideally this should be recorded
	 * by storing the QE and the amplification. However the total gain (Count/photon)
	 * is already stored with the results. Thus the amplification can be inferred by
	 * dividing the total gain by the quantum efficiency which should be in the range 0-1.
	 *
	 * @return the quantum efficiency
	 */
	public double getQuantumEfficiency()
	{
		CalibrationOrBuilder c = getCalibrationOrBuilder();
		return (c.hasCameraCalibration()) ? Maths.clip(0, 1, c.getCameraCalibration().getQuantumEfficiency()) : 0;
	}

	/**
	 * Checks for quantum efficiency (e-/photon).
	 *
	 * @return true, if successful (QE is in the range 0 (exclusive) to 1 (inclusive)).
	 */
	public boolean hasQuantumEfficiency()
	{
		double qe = getQuantumEfficiency();
		return qe > 0 && qe <= 1;
	}

	/**
	 * Get the camera amplification (Count/e-) used when modelling a microscope camera.
	 * <p>
	 * Note that the camera noise model assumes that electrons are converted to Count units by amplification that is not
	 * perfect (i.e. it has noise). The amplification is equal to the gain (Count/photon) divided by the quantum
	 * efficiency [QE] (e-/photon).
	 * <p>
	 * If the QE is not set then it is assumed to be 1.
	 *
	 * @return the amplification
	 */
	public double getCountPerElectron()
	{
		double countPerPhoton = getCountPerPhoton();
		if (countPerPhoton > 0)
		{
			double qe = getQuantumEfficiency();
			if (qe == 0)
				qe = 1;
			return countPerPhoton / qe;
		}
		return 0;
	}

	/**
	 * Checks for amplification (Count/e-). This is true if there is a count/photon since the quantum efficiency is
	 * assumed to be 1 if absent.
	 *
	 * @return true, if successful
	 */
	public boolean hasCountPerElectron()
	{
		return hasCountPerPhoton(); // && hasQuantumEfficiency();
	}

	/**
	 * Gets the camera model name. This should contain all the information required to load the camera model, e.g. in
	 * the case of a per-pixel camera model for sCMOS cameras.
	 *
	 * @return the camera model name
	 */
	public String getCameraModelName()
	{
		CalibrationOrBuilder c = getCalibrationOrBuilder();
		return (c.hasCameraCalibration()) ? c.getCameraCalibration().getCameraModelName() : null;
	}

	/**
	 * Checks for camera model name.
	 *
	 * @return true, if successful
	 */
	public boolean hasCameraModelName()
	{
		return !TextUtils.isNullOrEmpty(getCameraModelName());
	}

	/**
	 * Get the distance unit used for the results.
	 *
	 * @return the distanceUnit
	 */
	public DistanceUnit getDistanceUnit()
	{
		CalibrationOrBuilder c = getCalibrationOrBuilder();
		return (c.hasDistanceCalibration()) ? c.getDistanceCalibration().getDistanceUnit()
				: DistanceUnit.DISTANCE_UNIT_NA;
	}

	/**
	 * Checks for distance unit.
	 *
	 * @return true, if successful
	 */
	public boolean hasDistanceUnit()
	{
		return getDistanceUnit().getNumber() > 0;
	}

	/**
	 * Get the intensity unit used for the results.
	 *
	 * @return the intensityUnit
	 */
	public IntensityUnit getIntensityUnit()
	{
		CalibrationOrBuilder c = getCalibrationOrBuilder();
		return (c.hasIntensityCalibration()) ? c.getIntensityCalibration().getIntensityUnit()
				: IntensityUnit.INTENSITY_UNIT_NA;
	}

	/**
	 * Checks for intensity unit.
	 *
	 * @return true, if successful
	 */
	public boolean hasIntensityUnit()
	{
		return getIntensityUnit().getNumber() > 0;
	}

	/**
	 * Get the angle unit used for the results.
	 *
	 * @return the angleUnit
	 */
	public AngleUnit getAngleUnit()
	{
		CalibrationOrBuilder c = getCalibrationOrBuilder();
		return (c.hasAngleCalibration()) ? c.getAngleCalibration().getAngleUnit() : AngleUnit.ANGLE_UNIT_NA;
	}

	/**
	 * Checks for angle unit.
	 *
	 * @return true, if successful
	 */
	public boolean hasAngleUnit()
	{
		return getAngleUnit().getNumber() > 0;
	}
}
