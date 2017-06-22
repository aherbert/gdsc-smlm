package gdsc.smlm.data.config;

import gdsc.core.data.utils.ConversionException;
import gdsc.core.data.utils.IdentityTypeConverter;
import gdsc.core.data.utils.TypeConverter;
import gdsc.smlm.data.config.SMLMSettings.AngleUnit;
import gdsc.smlm.data.config.SMLMSettings.Calibration;
import gdsc.smlm.data.config.SMLMSettings.CalibrationOrBuilder;
import gdsc.smlm.data.config.SMLMSettings.CameraType;
import gdsc.smlm.data.config.SMLMSettings.DistanceCalibration;
import gdsc.smlm.data.config.SMLMSettings.DistanceUnit;
import gdsc.smlm.data.config.SMLMSettings.IntensityCalibration;
import gdsc.smlm.data.config.SMLMSettings.IntensityUnit;
import gdsc.smlm.data.config.SMLMSettings.PSFCalibration;

// TODO: Auto-generated Javadoc
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
	/** The original calibration. This is immutable. */
	private Calibration calibration;

	/** The calibration builder holding the latest calibration if changes have been made. */
	private Calibration.Builder calibrationBuilder;

	/**
	 * Instantiates a new calibration helper with a default calibration.
	 */
	public CalibrationHelper()
	{
		this.calibration = Calibration.getDefaultInstance();
	}

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
		this.calibration = calibration;
	}

	/**
	 * Instantiates a new calibration helper.
	 *
	 * @param calibration
	 *            the calibration
	 * @throws IllegalArgumentException
	 *             if the calibration is null
	 */
	public CalibrationHelper(Calibration.Builder calibration) throws IllegalArgumentException
	{
		if (calibration == null)
			throw new IllegalArgumentException("Calibration is null");
		this.calibrationBuilder = calibration;
	}

	/**
	 * Gets the latest calibration.
	 *
	 * @return the calibration
	 */
	public Calibration getCalibration()
	{
		return (calibrationBuilder != null) ? calibrationBuilder.build() : calibration;
	}

	/**
	 * Gets the builder containing the latest calibration.
	 *
	 * @return the builder
	 */
	public Calibration.Builder getBuilder()
	{
		if (calibrationBuilder == null)
			calibrationBuilder = calibration.toBuilder();
		return calibrationBuilder;
	}

	/**
	 * Gets the calibration or builder with the latest changes.
	 *
	 * @return the calibration or builder
	 */
	public CalibrationOrBuilder getCalibrationOrBuilder()
	{
		return (calibrationBuilder != null) ? calibrationBuilder : calibration;
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
	public static TypeConverter<DistanceUnit> getDistanceConverter(CalibrationOrBuilder calibration,
			DistanceUnit toDistanceUnit) throws ConversionException
	{
		if (calibration != null && toDistanceUnit != null && calibration.hasDistanceCalibration())
		{
			DistanceCalibration distanceCalibration = calibration.getDistanceCalibration();
			return UnitConverterFactory.createConverter(distanceCalibration.getUnit(), toDistanceUnit,
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
		return getDistanceConverter(getCalibrationOrBuilder(), toDistanceUnit);
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
		return getIntensityConverter(getCalibrationOrBuilder(), toIntensityUnit);
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
		return getAngleConverter(getCalibrationOrBuilder(), toAngleUnit);
	}

	/**
	 * Gets a distance converter to update values.
	 * <p>
	 * If the calibration is already in the given units or conversion is not possible
	 * then an identity converter will be returned.
	 *
	 * @param toDistanceUnit
	 *            the distance unit
	 * @return the distance converter
	 */
	public TypeConverter<DistanceUnit> getDistanceConverterSafe(DistanceUnit toDistanceUnit)
	{
		return getDistanceConverterSafe(getCalibrationOrBuilder(), toDistanceUnit);
	}

	/**
	 * Gets an intensity converter to update values.
	 * <p>
	 * If the calibration is already in the given units or conversion is not possible
	 * then an identity converter will be returned.
	 *
	 * @param toIntensityUnit
	 *            the intensity unit
	 * @return the intensity converter
	 */
	public TypeConverter<IntensityUnit> getIntensityConverterSafe(IntensityUnit toIntensityUnit)
	{
		return getIntensityConverterSafe(getCalibrationOrBuilder(), toIntensityUnit);
	}

	/**
	 * Gets an angle converter to update values.
	 * <p>
	 * If the calibration is already in the given units or conversion is not possible
	 * then an identity converter will be returned.
	 *
	 * @param toAngleUnit
	 *            the angle unit
	 * @return the angle converter
	 */
	public TypeConverter<AngleUnit> getAngleConverterSafe(AngleUnit toAngleUnit)
	{
		return getAngleConverterSafe(getCalibrationOrBuilder(), toAngleUnit);
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
	 * Sets the image pixel size in nanometers.
	 *
	 * @param nmPerPixel
	 *            image pixel size in nanometers
	 */
	public void setNmPerPixel(double nmPerPixel)
	{
		getBuilder().getDistanceCalibrationBuilder().setNmPerPixel(nmPerPixel);
	}

	/**
	 * Gets the gain (Count/photon). Can be used to convert the signal in count units to photons.
	 *
	 * @return the gain
	 */
	public double getGain()
	{
		CalibrationOrBuilder c = getCalibrationOrBuilder();
		return (c.hasIntensityCalibration()) ? c.getIntensityCalibration().getGain() : 0;
	}

	/**
	 * Checks for gain.
	 *
	 * @return true, if successful
	 */
	public boolean hasGain()
	{
		return getGain() > 0;
	}

	/**
	 * Sets the gain (Count/photon). Can be used to convert the signal in count units to photons.
	 *
	 * @param gain
	 *            the new gain
	 */
	public void setGain(double gain)
	{
		getBuilder().getIntensityCalibrationBuilder().setGain(gain);
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
	 * Sets the exposure time in milliseconds per frame.
	 *
	 * @param exposureTime
	 *            the new exposure time
	 */
	public void setExposureTime(double exposureTime)
	{
		getBuilder().getTimeCalibrationBuilder().setExposureTime(exposureTime);
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
	 * Sets camera Gaussian read noise (in Count units).
	 *
	 * @param readNoise
	 *            the new read noise
	 */
	public void setReadNoise(double readNoise)
	{
		getBuilder().getCameraCalibrationBuilder().setReadNoise(readNoise);
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

	public boolean hasBias()
	{
		CalibrationOrBuilder c = getCalibrationOrBuilder();
		return c.hasCameraCalibration();
	}

	/**
	 * Sets camera bias (in Count units).
	 *
	 * @param bias
	 *            the new bias
	 */
	public void setBias(double bias)
	{
		getBuilder().getCameraCalibrationBuilder().setReadNoise(bias);
	}

	/**
	 * Get the camera type.
	 *
	 * @return the camera type
	 */
	public CameraType getCameraType()
	{
		CalibrationOrBuilder c = getCalibrationOrBuilder();
		return (c.hasCameraCalibration()) ? c.getCameraCalibration().getCameraType() : null;
	}

	/**
	 * Checks for camera type.
	 *
	 * @return true, if successful
	 */
	public boolean hasCameraType()
	{
		return getCameraType() != null;
	}

	/**
	 * Set the camera type.
	 *
	 * @param cameraType
	 *            the new camera type
	 */
	public void setCameraType(CameraType cameraType)
	{
		if (cameraType == null)
			getBuilder().getCameraCalibrationBuilder().clearCameraType();
		else
			getBuilder().getCameraCalibrationBuilder().setCameraType(cameraType);
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
			switch (c.getCameraCalibration().getCameraType())
			{
				case CCD:
				case EMCCD:
					return true;
				default:
					break;
			}
		}
		return false;
	}

	/**
	 * Checks if the camera type was an Electron Multiplying (EM) CCD.
	 *
	 * @return true, if the camera type was an Electron Multiplying (EM) CCD
	 */
	public boolean isEmCCD()
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
	 * Get the camera amplification (Count/e-) used when modelling a microscope camera.
	 * <p>
	 * Note that the camera noise model assumes that electrons are converted to Count units by amplification that is not
	 * perfect (i.e. it has noise). The amplification is equal to the gain (Count/photon) divided by the quantum
	 * efficiency (e-/photon).
	 *
	 * @return the amplification
	 */
	public double getAmplification()
	{
		CalibrationOrBuilder c = getCalibrationOrBuilder();
		return (c.hasCameraCalibration()) ? c.getCameraCalibration().getAmplification() : 0;
	}

	/**
	 * Checks for amplification.
	 *
	 * @return true, if successful
	 */
	public boolean hasAmplification()
	{
		return getAmplification() > 0;
	}

	/**
	 * Set the camera amplification (Count/e-) used when modelling a microscope camera.
	 * <p>
	 * Note that the camera noise model assumes that electrons are converted to Count units by amplification that is not
	 * perfect (i.e. it has noise). The amplification is equal to the gain (Count/photon) divided by the quantum
	 * efficiency (e-/photon).
	 *
	 * @param amplification
	 *            the new amplification
	 */
	public void setAmplification(double amplification)
	{
		getBuilder().getCameraCalibrationBuilder().setAmplification(amplification);
	}

	/**
	 * Get the distance unit used for the results.
	 *
	 * @return the distanceUnit
	 */
	public DistanceUnit getDistanceUnit()
	{
		CalibrationOrBuilder c = getCalibrationOrBuilder();
		return (c.hasDistanceCalibration()) ? c.getDistanceCalibration().getUnit() : null;
	}

	/**
	 * Checks for distance unit.
	 *
	 * @return true, if successful
	 */
	public boolean hasDistanceUnit()
	{
		return getDistanceUnit() != null;
	}

	/**
	 * Set the distance unit used for the results.
	 *
	 * @param distanceUnit
	 *            the new distanceUnit
	 */
	public void setDistanceUnit(DistanceUnit distanceUnit)
	{
		if (distanceUnit == null)
			getBuilder().getDistanceCalibrationBuilder().clearUnit();
		else
			getBuilder().getDistanceCalibrationBuilder().setUnit(distanceUnit);
	}

	/**
	 * Get the intensity unit used for the results.
	 *
	 * @return the intensityUnit
	 */
	public IntensityUnit getIntensityUnit()
	{
		CalibrationOrBuilder c = getCalibrationOrBuilder();
		return (c.hasIntensityCalibration()) ? c.getIntensityCalibration().getUnit() : null;
	}

	/**
	 * Checks for intensity unit.
	 *
	 * @return true, if successful
	 */
	public boolean hasIntensityUnit()
	{
		return getIntensityUnit() != null;
	}

	/**
	 * Set the intensity unit used for the results.
	 *
	 * @param intensityUnit
	 *            the new intensityUnit
	 */
	public void setIntensityUnit(IntensityUnit intensityUnit)
	{
		if (intensityUnit == null)
			getBuilder().getIntensityCalibrationBuilder().clearUnit();
		else
			getBuilder().getIntensityCalibrationBuilder().setUnit(intensityUnit);
	}

	/**
	 * Get the angle unit used for the results.
	 *
	 * @return the angleUnit
	 */
	public AngleUnit getAngleUnit()
	{
		CalibrationOrBuilder c = getCalibrationOrBuilder();
		return (c.hasPsfCalibration()) ? c.getPsfCalibration().getAngleUnit() : null;
	}

	/**
	 * Checks for angle unit.
	 *
	 * @return true, if successful
	 */
	public boolean hasAngleUnit()
	{
		return getAngleUnit() != null;
	}

	/**
	 * Set the angle unit used for the results.
	 *
	 * @param angleUnit
	 *            The new angleUnit
	 */
	public void setAngleUnit(AngleUnit angleUnit)
	{
		if (angleUnit == null)
			getBuilder().getPsfCalibrationBuilder().clearAngleUnit();
		else
			getBuilder().getPsfCalibrationBuilder().setAngleUnit(angleUnit);
	}
}
