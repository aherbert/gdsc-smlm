package gdsc.smlm.results;

import java.util.ArrayList;

import gdsc.smlm.units.DistanceUnit;
import gdsc.smlm.units.IdentityUnitConverter;
import gdsc.smlm.units.IntensityUnit;
import gdsc.smlm.units.UnitConversionException;
import gdsc.smlm.units.UnitConverter;

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
 * Contain the calibration settings for the microscope
 * <p>
 * The calibration has flags to indicate that a valid value has been set for each property. If these are false then the
 * property get method can optionally throw an exception.
 */
public class Calibration implements Cloneable
{
	/**
	 * The camera type.
	 */
	public enum CameraType
	{

		/** The ccd. */
		//@formatter:off
		CCD { String getName() {return "CCD"; } }, 
		
		/** The em ccd. */
		EM_CCD { String getName() {return "EM-CCD"; } },
		
		/** The scmos. */
		SCMOS { String getName() {return "sCMOS"; } },
		;
		//@formatter:on

		/**
		 * Gets the name.
		 *
		 * @return the name
		 */
		abstract String getName();

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Enum#toString()
		 */
		public String toString()
		{
			return getName();
		}
	}

	/** The field missing exception. */
	// State flags
	private static int FIELD_MISSING_EXCEPTION = 0x00000001;

	/** The field nm per pixel. */
	private static int FIELD_NM_PER_PIXEL = 0x00000002;

	/** The field gain. */
	private static int FIELD_GAIN = 0x00000004;

	/** The field exposure time. */
	private static int FIELD_EXPOSURE_TIME = 0x00000008;

	/** The field read noise. */
	private static int FIELD_READ_NOISE = 0x00000010;

	/** The field bias. */
	private static int FIELD_BIAS = 0x00000020;

	/** The field camera type. */
	private static int FIELD_CAMERA_TYPE = 0x00000040;

	/** The field amplification. */
	private static int FIELD_AMPLIFICATION = 0x00000080;

	/** The field distance unit. */
	private static int FIELD_DISTANCE_UNIT = 0x00000100;

	/** The field intensity unit. */
	private static int FIELD_INTENSITY_UNIT = 0x00000200;

	/** The fields. */
	private int fields = 0;

	/**
	 * Checks if an exception will be thrown when accessing a field that is missing.
	 *
	 * @return true, if is field missing exceptions are enabled
	 */
	public boolean isFieldMissingException()
	{
		return ((fields & FIELD_MISSING_EXCEPTION) == FIELD_MISSING_EXCEPTION);
	}

	/**
	 * Sets the field missing exception enabled flag. If true an exception will be thrown when accessing a field that is
	 * missing.
	 *
	 * @param enabled
	 *            the new field missing exception enabled flag
	 */
	public void setFieldMissingException(boolean enabled)
	{
		if (enabled)
			fields |= FIELD_MISSING_EXCEPTION;
		else
			fields = fields & ~FIELD_MISSING_EXCEPTION;
	}

	/**
	 * Checks for nm per pixel.
	 *
	 * @return true, if successful
	 */
	public boolean hasNmPerPixel()
	{
		return ((fields & FIELD_NM_PER_PIXEL) == FIELD_NM_PER_PIXEL);
	}

	/**
	 * Checks for gain.
	 *
	 * @return true, if successful
	 */
	public boolean hasGain()
	{
		return ((fields & FIELD_GAIN) == FIELD_GAIN);
	}

	/**
	 * Checks for exposure time.
	 *
	 * @return true, if successful
	 */
	public boolean hasExposureTime()
	{
		return ((fields & FIELD_EXPOSURE_TIME) == FIELD_EXPOSURE_TIME);
	}

	/**
	 * Checks for read noise.
	 *
	 * @return true, if successful
	 */
	public boolean hasReadNoise()
	{
		return ((fields & FIELD_READ_NOISE) == FIELD_READ_NOISE);
	}

	/**
	 * Checks for bias.
	 *
	 * @return true, if successful
	 */
	public boolean hasBias()
	{
		return ((fields & FIELD_BIAS) == FIELD_BIAS);
	}

	/**
	 * Checks for camera type.
	 *
	 * @return true, if successful
	 */
	public boolean hasCameraType()
	{
		return ((fields & FIELD_CAMERA_TYPE) == FIELD_CAMERA_TYPE);
	}

	/**
	 * Checks for amplification.
	 *
	 * @return true, if successful
	 */
	public boolean hasAmplification()
	{
		return ((fields & FIELD_AMPLIFICATION) == FIELD_AMPLIFICATION);
	}

	/**
	 * Checks for distance unit.
	 *
	 * @return true, if successful
	 */
	public boolean hasDistanceUnit()
	{
		return ((fields & FIELD_DISTANCE_UNIT) == FIELD_DISTANCE_UNIT);
	}

	/**
	 * Checks for intensity unit.
	 *
	 * @return true, if successful
	 */
	public boolean hasIntensityUnit()
	{
		return ((fields & FIELD_INTENSITY_UNIT) == FIELD_INTENSITY_UNIT);
	}

	/**
	 * Sets the has nm per pixel.
	 */
	private void setHasNmPerPixel()
	{
		fields |= FIELD_NM_PER_PIXEL;
	}

	/**
	 * Sets the has gain.
	 */
	private void setHasGain()
	{
		fields |= FIELD_GAIN;
	}

	/**
	 * Sets the has exposure time.
	 */
	private void setHasExposureTime()
	{
		fields |= FIELD_EXPOSURE_TIME;
	}

	/**
	 * Sets the has read noise.
	 */
	private void setHasReadNoise()
	{
		fields |= FIELD_READ_NOISE;
	}

	/**
	 * Sets the has bias.
	 */
	private void setHasBias()
	{
		fields |= FIELD_BIAS;
	}

	/**
	 * Sets the has camera type.
	 */
	private void setHasCameraType()
	{
		fields |= FIELD_CAMERA_TYPE;
	}

	/**
	 * Sets the has amplification.
	 */
	private void setHasAmplification()
	{
		fields |= FIELD_AMPLIFICATION;
	}

	/**
	 * Sets the has distance unit.
	 */
	private void setHasDistanceUnit()
	{
		fields |= FIELD_DISTANCE_UNIT;
	}

	/**
	 * Sets the has intensity unit.
	 */
	private void setHasIntensityUnit()
	{
		fields |= FIELD_INTENSITY_UNIT;
	}

	/**
	 * Clear has nm per pixel.
	 */
	public void clearHasNmPerPixel()
	{
		fields = fields & ~FIELD_NM_PER_PIXEL;
		nmPerPixel = 0;
	}

	/**
	 * Clear has gain.
	 */
	public void clearHasGain()
	{
		fields = fields & ~FIELD_GAIN;
		gain = 0;
	}

	/**
	 * Clear has exposure time.
	 */
	public void clearHasExposureTime()
	{
		fields = fields & ~FIELD_EXPOSURE_TIME;
		exposureTime = 0;
	}

	/**
	 * Clear has read noise.
	 */
	public void clearHasReadNoise()
	{
		fields = fields & ~FIELD_READ_NOISE;
		readNoise = -1;
	}

	/**
	 * Clear has bias.
	 */
	public void clearHasBias()
	{
		fields = fields & ~FIELD_BIAS;
		bias = -1;
	}

	/**
	 * Clear has camera type.
	 */
	public void clearHasCameraType()
	{
		fields = fields & ~FIELD_CAMERA_TYPE;
		cameraType = null;
	}

	/**
	 * Clear has EMCCD.
	 */
	@Deprecated
	public void clearHasEMCCD()
	{
		clearHasCameraType();
	}

	/**
	 * Clear has amplification.
	 */
	public void clearHasAmplification()
	{
		fields = fields & ~FIELD_AMPLIFICATION;
		amplification = 0;
	}

	/**
	 * Clear has distance unit.
	 */
	public void clearHasDistanceUnit()
	{
		fields = fields & ~FIELD_DISTANCE_UNIT;
		distanceUnit = null;
	}

	/**
	 * Clear has intensity unit.
	 */
	public void clearHasIntensityUnit()
	{
		fields = fields & ~FIELD_INTENSITY_UNIT;
		intensityUnit = null;
	}

	/** The image pixel size in nanometers. */
	private double nmPerPixel = 0;
	/**
	 * The gain (ADUs/photon). Can be used to convert the signal in Analogue-to-Digital units
	 * (ADUs) to photons.
	 */
	private double gain = 0;

	/** The exposure time in milliseconds per frame. */
	private double exposureTime = 0;

	/** The camera Gaussian read noise (in ADUs). */
	private double readNoise = -1;

	/** The camera bias (in ADUs). */
	private double bias = -1;

	/** The camera type. */
	private CameraType cameraType = null;

	/**
	 * True if the camera was run in Electron Multiplying (EM) mode.
	 *
	 * @deprecated This has been replaced by camaraType. It is left to enable XStream deserialisation of old
	 *             configuration.
	 */
	@Deprecated
	private boolean emCCD;
	/**
	 * The camera amplification (ADUs/e-) used when modelling a microscope camera.
	 * <p>
	 * Note that the camera noise model assumes that electrons are converted to ADUs by amplification that is not
	 * perfect (i.e. it has noise). The amplification is equal to the gain (ADUs/photon) divided by the quantum
	 * efficiency (e-/photon).
	 */
	private double amplification = 0;

	/** The distance unit. */
	private DistanceUnit distanceUnit = null;

	/** The intensity unit. */
	private IntensityUnit intensityUnit = null;

	/**
	 * Default constructor. All properties are set to be invalid but field missing exceptions are disabled.
	 */
	public Calibration()
	{
	}

	/**
	 * All properties are set to be invalid but missing exceptions can be enabled.
	 *
	 * @param fieldMissingExceptionEnabled
	 *            Set to true to enable the field missing exceptions in the get methods
	 */
	public Calibration(boolean fieldMissingExceptionEnabled)
	{
		setFieldMissingException(fieldMissingExceptionEnabled);
	}

	/**
	 * Parameterised constructor for essential settings.
	 *
	 * @param nmPerPixel
	 *            the nm per pixel
	 * @param gain
	 *            the gain
	 * @param exposureTime
	 *            the exposure time
	 */
	public Calibration(double nmPerPixel, double gain, double exposureTime)
	{
		this.setNmPerPixel(nmPerPixel);
		this.setGain(gain);
		this.setExposureTime(exposureTime);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public Calibration clone()
	{
		Calibration c;
		try
		{
			c = (Calibration) super.clone();
		}
		catch (CloneNotSupportedException e)
		{
			return null;
		}
		return c;
	}

	/**
	 * Clear the calibration (set all valid flags to false).
	 */
	public void clear()
	{
		// Clear all but the field missing flag
		fields = fields & FIELD_MISSING_EXCEPTION;
	}

	/**
	 * Validate the calibration by calling each property setter with the current value. This will set the valid flags to
	 * false If the current value is not valid.
	 */
	public void validate()
	{
		setNmPerPixel(nmPerPixel);
		setGain(gain);
		setExposureTime(exposureTime);
		setReadNoise(readNoise);
		setBias(bias);
		setCameraType(cameraType);
		setAmplification(amplification);
		setDistanceUnit(distanceUnit);
		setIntensityUnit(intensityUnit);
	}

	/**
	 * Gets image pixel size in nanometers.
	 *
	 * @return image pixel size in nanometers
	 * @throws IllegalStateException
	 *             if the missing field exceptions is enabled and the field has not been set to a valid value
	 */
	public double getNmPerPixel()
	{
		if (isFieldMissingException() && !hasNmPerPixel())
			throw new IllegalStateException();
		return nmPerPixel;
	}

	/**
	 * Sets the image pixel size in nanometers.
	 *
	 * @param nmPerPixel
	 *            image pixel size in nanometers
	 */
	public void setNmPerPixel(double nmPerPixel)
	{
		if (nmPerPixel > 0)
		{
			setHasNmPerPixel();
			this.nmPerPixel = nmPerPixel;
		}
		else
			clearHasNmPerPixel();
	}

	/**
	 * Gets the gain (ADUs/photon). Can be used to convert the signal in Analogue-to-Digital units
	 * (ADUs) to photons.
	 *
	 * @return the gain
	 * @throws IllegalStateException
	 *             if the missing field exceptions is enabled and the field has not been set to a valid value
	 */
	public double getGain()
	{
		if (isFieldMissingException() && !hasGain())
			throw new IllegalStateException();
		return gain;
	}

	/**
	 * Sets the (ADUs/photon). Can be used to convert the signal in Analogue-to-Digital units
	 * (ADUs) to photons.
	 *
	 * @param gain
	 *            the new gain
	 */
	public void setGain(double gain)
	{
		if (gain > 0)
		{
			setHasGain();
			this.gain = gain;
		}
		else
			clearHasGain();
	}

	/**
	 * Gets the exposure time in milliseconds per frame.
	 *
	 * @return the exposure time
	 * @throws IllegalStateException
	 *             if the missing field exceptions is enabled and the field has not been set to a valid value
	 */
	public double getExposureTime()
	{
		if (isFieldMissingException() && !hasExposureTime())
			throw new IllegalStateException();
		return exposureTime;
	}

	/**
	 * Sets the exposure time in milliseconds per frame.
	 *
	 * @param exposureTime
	 *            the new exposure time
	 */
	public void setExposureTime(double exposureTime)
	{
		if (exposureTime > 0)
		{
			setHasExposureTime();
			this.exposureTime = exposureTime;
		}
		else
			clearHasExposureTime();
	}

	/**
	 * Gets the camera Gaussian read noise (in ADUs).
	 *
	 * @return the camera Gaussian read noise (in ADUs)
	 * @throws IllegalStateException
	 *             if the missing field exceptions is enabled and the field has not been set to a valid value
	 */
	public double getReadNoise()
	{
		if (isFieldMissingException() && !hasReadNoise())
			throw new IllegalStateException();
		return readNoise;
	}

	/**
	 * Sets camera Gaussian read noise (in ADUs).
	 *
	 * @param readNoise
	 *            the new read noise
	 */
	public void setReadNoise(double readNoise)
	{
		if (readNoise >= 0)
		{
			setHasReadNoise();
			this.readNoise = readNoise;
		}
		else
			clearHasReadNoise();
	}

	/**
	 * Gets camera bias (in ADUs).
	 *
	 * @return camera bias (in ADUs)
	 * @throws IllegalStateException
	 *             if the missing field exceptions is enabled and the field has not been set to a valid value
	 */
	public double getBias()
	{
		if (isFieldMissingException() && !hasBias())
			throw new IllegalStateException();
		return bias;
	}

	/**
	 * Sets camera bias (in ADUs).
	 *
	 * @param bias
	 *            the new bias
	 */
	public void setBias(double bias)
	{
		if (bias >= 0)
		{
			setHasBias();
			this.bias = bias;
		}
		else
			clearHasBias();
	}

	/**
	 * Get the camera type.
	 *
	 * @return the camera type
	 * @throws IllegalStateException
	 *             if the missing field exceptions is enabled and the field has not been set to a valid value
	 */
	public CameraType getCameraType()
	{
		if (isFieldMissingException() && !hasCameraType())
			throw new IllegalStateException();
		return cameraType;
	}

	/**
	 * Set the camera type.
	 *
	 * @param cameraType
	 *            the new camera type
	 */
	public void setCameraType(CameraType cameraType)
	{
		if (cameraType != null)
		{
			setHasCameraType();
			this.cameraType = cameraType;
		}
		else
			clearHasCameraType();
	}

	/**
	 * Checks for a CCD camera.
	 *
	 * @return true, if successful
	 */
	public boolean isCCDCamera()
	{
		return hasCameraType() && (cameraType == CameraType.CCD || cameraType == CameraType.EM_CCD);
	}
	
	/**
	 * Checks if the camera type was an Electron Multiplying (EM) CCD.
	 *
	 * @return true, if the camera type was an Electron Multiplying (EM) CCD
	 * @throws IllegalStateException
	 *             if the missing field exceptions is enabled and the camera type field has not been set
	 */
	public boolean isEmCCD()
	{
		return getCameraType() == CameraType.EM_CCD;
	}
	
	/**
	 * Checks if the camera type was a standard CCD.
	 *
	 * @return true, if the camera type was a standard CCD.
	 * @throws IllegalStateException
	 *             if the missing field exceptions is enabled and the camera type field has not been set
	 */
	public boolean isCCD()
	{
		return getCameraType() == CameraType.CCD;
	}

	/**
	 * Set if the CCD camera was run in Electron Multiplying (EM) mode, otherwise assume a standard CCD.
	 *
	 * @param emCCD
	 *            true, if the CCD camera was run in Electron Multiplying (EM) mode, otherwise set to CCD
	 * @deprecated Replaced by camera type
	 */
	@Deprecated
	public void setEmCCD(boolean emCCD)
	{
		setCameraType((emCCD) ? CameraType.EM_CCD : CameraType.CCD);
	}

	/**
	 * Sets the camera type from the deprecated EM CCD field.
	 */
	void setCameraTypeFromEmCCDField()
	{
		setCameraType((emCCD) ? CameraType.EM_CCD : CameraType.CCD);
	}

	/**
	 * Get the camera amplification (ADUs/e-) used when modelling a microscope camera.
	 * <p>
	 * Note that the camera noise model assumes that electrons are converted to ADUs by amplification that is not
	 * perfect (i.e. it has noise). The amplification is equal to the gain (ADUs/photon) divided by the quantum
	 * efficiency (e-/photon).
	 *
	 * @return the amplification
	 * @throws IllegalStateException
	 *             if the missing field exceptions is enabled and the field has not been set to a valid value
	 */
	public double getAmplification()
	{
		if (isFieldMissingException() && !hasAmplification())
			throw new IllegalStateException();
		return amplification;
	}

	/**
	 * Set the camera amplification (ADUs/e-) used when modelling a microscope camera.
	 * <p>
	 * Note that the camera noise model assumes that electrons are converted to ADUs by amplification that is not
	 * perfect (i.e. it has noise). The amplification is equal to the gain (ADUs/photon) divided by the quantum
	 * efficiency (e-/photon).
	 *
	 * @param amplification
	 *            the new amplification
	 */
	public void setAmplification(double amplification)
	{
		if (amplification > 0)
		{
			setHasAmplification();
			this.amplification = amplification;
		}
		else
			clearHasAmplification();
	}

	/**
	 * Get the distance unit used for the results.
	 *
	 * @return the distanceUnit
	 * @throws IllegalStateException
	 *             if the missing field exceptions is enabled and the field has not been set to a valid value
	 */
	public DistanceUnit getDistanceUnit()
	{
		if (isFieldMissingException() && !hasDistanceUnit())
			throw new IllegalStateException();
		return distanceUnit;
	}

	/**
	 * Set the distance unit used for the results.
	 *
	 * @param distanceUnit
	 *            the new distanceUnit
	 */
	public void setDistanceUnit(DistanceUnit distanceUnit)
	{
		if (distanceUnit != null)
		{
			setHasDistanceUnit();
			this.distanceUnit = distanceUnit;
		}
		else
			clearHasDistanceUnit();
	}

	/**
	 * Get the intensity unit used for the results.
	 *
	 * @return the intensityUnit
	 * @throws IllegalStateException
	 *             if the missing field exceptions is enabled and the field has not been set to a valid value
	 */
	public IntensityUnit getIntensityUnit()
	{
		if (isFieldMissingException() && !hasIntensityUnit())
			throw new IllegalStateException();
		return intensityUnit;
	}

	/**
	 * Set the intensity unit used for the results.
	 *
	 * @param intensityUnit
	 *            The new intensityUnit
	 */
	public void setIntensityUnit(IntensityUnit intensityUnit)
	{
		if (intensityUnit != null)
		{
			setHasIntensityUnit();
			this.intensityUnit = intensityUnit;
		}
		else
			clearHasIntensityUnit();
	}

	/**
	 * Gets a distance converter to update values. If the converter can be created then the current distance unit in
	 * this instance is updated.
	 * <p>
	 * If the calibration is already in the given units or conversion is not possible
	 * then an identity converter will be returned.
	 * <p>
	 * It is recommended to clone the calibration before invoking this method as the state may change.
	 *
	 * @param toDistanceUnit
	 *            the distance unit
	 * @return the distance converter
	 */
	public UnitConverter<DistanceUnit> getDistanceConverter(DistanceUnit toDistanceUnit)
	{
		UnitConverter<DistanceUnit> c = null;
		if (toDistanceUnit != null && hasDistanceUnit() && distanceUnit != toDistanceUnit && hasNmPerPixel())
		{
			try
			{
				c = distanceUnit.createConverter(toDistanceUnit, nmPerPixel);
				setDistanceUnit(toDistanceUnit);
			}
			catch (UnitConversionException e)
			{
				// Ignore this
			}
		}
		if (c == null)
			c = new IdentityUnitConverter<DistanceUnit>(distanceUnit);
		return c;
	}

	/**
	 * Gets intensity converters to update values. If the converter can be created then the current intensity unit in
	 * this instance is updated.
	 * <p>
	 * If the calibration is already in the given units or conversion is not possible
	 * then an identity converter will be returned.
	 * <p>
	 * The returned list has a converter with only the gain, and a second converter with the gain and bias.
	 * <p>
	 * It is recommended to clone the calibration before invoking this method as the state may change.
	 *
	 * @param toIntensityUnit
	 *            the intensity unit
	 * @return the intensity converters (gain, gain + bias)
	 */
	public ArrayList<UnitConverter<IntensityUnit>> getIntensityConverter(IntensityUnit toIntensityUnit)
	{
		ArrayList<UnitConverter<IntensityUnit>> list = new ArrayList<UnitConverter<IntensityUnit>>(2);
		if (toIntensityUnit != null && hasIntensityUnit() && intensityUnit != toIntensityUnit && hasGain() && hasBias())
		{
			try
			{
				list.add(intensityUnit.createConverter(toIntensityUnit, gain));
				list.add(intensityUnit.createConverter(toIntensityUnit, bias, gain));
				setIntensityUnit(toIntensityUnit);
			}
			catch (UnitConversionException e)
			{
				// Ignore this
			}
		}
		if (list.size() != 2)
		{
			list.clear();
			UnitConverter<IntensityUnit> c = new IdentityUnitConverter<IntensityUnit>(intensityUnit);
			list.add(c);
			list.add(c);
		}
		return list;
	}
}
