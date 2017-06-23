package gdsc.smlm.data.config;

import gdsc.smlm.data.config.UnitConfig.AngleUnit;
import gdsc.smlm.data.config.CalibrationConfig.Calibration;
import gdsc.smlm.data.config.CalibrationConfig.CalibrationOrBuilder;
import gdsc.smlm.data.config.CalibrationConfig.CameraType;
import gdsc.smlm.data.config.UnitConfig.DistanceUnit;
import gdsc.smlm.data.config.UnitConfig.IntensityUnit;

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
 * Contains helper functions for the writing a Calibration class.
 */
public class CalibrationWriter extends CalibrationReader
{
	/** The original calibration. This is immutable. */
	private Calibration calibration;

	/** The calibration builder holding the latest calibration if changes have been made. */
	private Calibration.Builder calibrationBuilder;

	/**
	 * Instantiates a new calibration writer with a default calibration.
	 */
	public CalibrationWriter()
	{
		this.calibration = Calibration.getDefaultInstance();
	}

	/**
	 * Instantiates a new calibration writer.
	 *
	 * @param calibration
	 *            the calibration
	 * @throws IllegalArgumentException
	 *             if the calibration is null
	 */
	public CalibrationWriter(Calibration calibration) throws IllegalArgumentException
	{
		if (calibration == null)
			throw new IllegalArgumentException("Calibration is null");
		this.calibration = calibration;
	}

	/**
	 * Instantiates a new calibration writer.
	 *
	 * @param calibration
	 *            the calibration
	 * @throws IllegalArgumentException
	 *             if the calibration is null
	 */
	public CalibrationWriter(Calibration.Builder calibration) throws IllegalArgumentException
	{
		if (calibration == null)
			throw new IllegalArgumentException("Calibration is null");
		this.calibrationBuilder = calibration;
	}

	/**
	 * Creates the calibration writer.
	 *
	 * @param calibration
	 *            the calibration (may be null)
	 * @return the calibration writer
	 */
	public static CalibrationWriter create(Calibration calibration)
	{
		return (calibration != null) ? new CalibrationWriter(calibration) : new CalibrationWriter();
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

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.data.config.CalibrationReader#getCalibrationOrBuilder()
	 */
	public CalibrationOrBuilder getCalibrationOrBuilder()
	{
		// If changes have been made then return the builder
		return (calibrationBuilder != null) ? calibrationBuilder : calibration;
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
	 * Sets camera bias (in Count units).
	 *
	 * @param bias
	 *            the new bias
	 */
	public void setBias(double bias)
	{
		getBuilder().getCameraCalibrationBuilder().setBias(bias);
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
	 * Sets the camera type to EMCCD or CCD.
	 *
	 * @param isEMCCD true if an EMCCD
	 * @deprecated This has been replaced by {@link #setCameraType(CameraType)}
	 */
	@Deprecated
	public void setEmCCD(boolean isEMCCD)
	{
		setCameraType((isEMCCD) ? CameraType.EMCCD : CameraType.CCD);
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
	 * Set the distance unit used for the results.
	 *
	 * @param distanceUnit
	 *            the new distanceUnit
	 */
	public void setDistanceUnit(DistanceUnit distanceUnit)
	{
		if (distanceUnit == null)
			getBuilder().getDistanceCalibrationBuilder().clearDistanceUnit();
		else
			getBuilder().getDistanceCalibrationBuilder().setDistanceUnit(distanceUnit);
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
			getBuilder().getIntensityCalibrationBuilder().clearIntensityUnit();
		else
			getBuilder().getIntensityCalibrationBuilder().setIntensityUnit(intensityUnit);
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
