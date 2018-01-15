package gdsc.smlm.data.config;

import gdsc.smlm.data.config.UnitProtos.AngleUnit;
import gdsc.core.utils.TextUtils;
import gdsc.smlm.data.config.CalibrationProtos.Calibration;
import gdsc.smlm.data.config.CalibrationProtos.CalibrationOrBuilder;
import gdsc.smlm.data.config.CalibrationProtos.CameraCalibration;
import gdsc.smlm.data.config.CalibrationProtos.CameraType;
import gdsc.smlm.data.config.FitProtos.PrecisionMethod;
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
	 * Merge the calibration.
	 *
	 * @param calibration
	 *            the new calibration
	 */
	public void mergeCalibration(Calibration calibration)
	{
		getBuilder().mergeFrom(calibration);
	}

	/**
	 * Sets the calibration.
	 *
	 * @param calibration
	 *            the new calibration
	 */
	public void setCalibration(Calibration calibration)
	{
		getBuilder().clear().mergeFrom(calibration);
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
		if (nmPerPixel < 0)
			throw new IllegalArgumentException();
		getBuilder().getDistanceCalibrationBuilder().setNmPerPixel(nmPerPixel);
	}

	/**
	 * Sets the gain (Count/photon). Can be used to convert the signal in count units to photons.
	 *
	 * @param gain
	 *            the new gain
	 */
	public void setCountPerPhoton(double gain)
	{
		if (gain < 0)
			throw new IllegalArgumentException();
		getBuilder().getIntensityCalibrationBuilder().setCountPerPhoton(gain);
	}

	/**
	 * Sets the exposure time in milliseconds per frame.
	 *
	 * @param exposureTime
	 *            the new exposure time
	 */
	public void setExposureTime(double exposureTime)
	{
		if (exposureTime < 0)
			throw new IllegalArgumentException();
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
		if (readNoise < 0)
			throw new IllegalArgumentException();
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
		if (bias < 0)
			throw new IllegalArgumentException();
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
			getBuilder().getCameraCalibrationBuilder().setCameraTypeValue(cameraType.getNumber());
	}

	/**
	 * Sets the camera type to EMCCD or CCD.
	 *
	 * @param isEMCCD
	 *            true if an EMCCD
	 * @deprecated This has been replaced by {@link #setCameraType(CameraType)}
	 */
	@Deprecated
	public void setEmCCD(boolean isEMCCD)
	{
		setCameraType((isEMCCD) ? CameraType.EMCCD : CameraType.CCD);
	}

	/**
	 * Sets the camera quantum efficiency (e-/photon) used when modelling a microscope camera.
	 * <p>
	 * Note that the camera noise model assumes that photons are converted to counts by
	 * a process that is not perfect (i.e. it has noise). The underlying process is
	 * photons converted to electrons in the camera chip and then amplification
	 * (count/electron) occurring in the camera hardware. Ideally this should be recorded
	 * by storing the QE and the amplification. However the total gain (Count/photon)
	 * is already stored with the results. Thus the amplification can be inferred by
	 * dividing the total gain by the quantum efficiency which should be in the range 0-1.
	 *
	 * @param quantumEfficiency
	 *            the new quantum efficiency
	 */
	public void setQuantumEfficiency(double quantumEfficiency)
	{
		if (quantumEfficiency < 0 || quantumEfficiency > 1)
			throw new IllegalArgumentException();
		getBuilder().getCameraCalibrationBuilder().setQuantumEfficiency(quantumEfficiency);
	}

	/**
	 * Sets the camera model name. This should contain all the information required to load the camera model, e.g. in
	 * the case of a per-pixel camera model for sCMOS cameras.
	 *
	 * @param cameraModelName
	 *            the new camera model name
	 */
	public void setCameraModelName(String cameraModelName)
	{
		if (TextUtils.isNullOrEmpty(cameraModelName))
			getBuilder().getCameraCalibrationBuilder().clearCameraModelName();
		else
			getBuilder().getCameraCalibrationBuilder().setCameraModelName(cameraModelName);
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
			getBuilder().getDistanceCalibrationBuilder().setDistanceUnitValue(distanceUnit.getNumber());
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
			getBuilder().getIntensityCalibrationBuilder().setIntensityUnitValue(intensityUnit.getNumber());
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
			getBuilder().getAngleCalibrationBuilder().clearAngleUnit();
		else
			getBuilder().getAngleCalibrationBuilder().setAngleUnitValue(angleUnit.getNumber());
	}

	/**
	 * Set the precision method used for the results.
	 *
	 * @param precision method
	 *            The new precision method
	 */
	public void setPrecisionMethod(PrecisionMethod precisionMethod)
	{
		if (precisionMethod == null)
			getBuilder().getResultDataCalibrationBuilder().clearPrecisionMethod();
		else
			getBuilder().getResultDataCalibrationBuilder().setPrecisionMethodValue(precisionMethod.getNumber());
	}

	/**
	 * Clear global camera settings. This should be used to remove the global camera settings, e.g. bias, read noise,
	 * count/photon, for example when a per-pixel camera model is used such as sCMOS camera type.
	 */
	public void clearGlobalCameraSettings()
	{
		CameraCalibration.Builder b = getBuilder().getCameraCalibrationBuilder();
		b.setBias(0);
		b.setReadNoise(0);
		setCountPerPhoton(0); // gain
		
		// This is still relevant to per-pixel camera settings as it is related to the chip sensitivity
		//b.setQuantumEfficiency(0);   
	}
}
