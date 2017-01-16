package gdsc.smlm.results;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
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
	// State flags
	private static int FIELD_MISSING_EXCEPTION = 0x00000001;
	private static int FIELD_NM_PER_PIXEL = 0x00000002;
	private static int FIELD_GAIN = 0x00000004;
	private static int FIELD_EXPOSURE_TIME = 0x00000010;
	private static int FIELD_READ_NOISE = 0x00000020;
	private static int FIELD_BIAS = 0x00000040;
	private static int FIELD_EM_CCD = 0x00000100;
	private static int FIELD_AMPLIFICATION = 0x00000200;
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

	public boolean hasNmPerPixel()
	{
		return ((fields & FIELD_NM_PER_PIXEL) == FIELD_NM_PER_PIXEL);
	}

	public boolean hasGain()
	{
		return ((fields & FIELD_GAIN) == FIELD_GAIN);
	}

	public boolean hasExposureTime()
	{
		return ((fields & FIELD_EXPOSURE_TIME) == FIELD_EXPOSURE_TIME);
	}

	public boolean hasReadNoise()
	{
		return ((fields & FIELD_READ_NOISE) == FIELD_READ_NOISE);
	}

	public boolean hasBias()
	{
		return ((fields & FIELD_BIAS) == FIELD_BIAS);
	}

	public boolean hasEMCCD()
	{
		return ((fields & FIELD_EM_CCD) == FIELD_EM_CCD);
	}

	public boolean hasAmplification()
	{
		return ((fields & FIELD_AMPLIFICATION) == FIELD_AMPLIFICATION);
	}

	private void setHasNmPerPixel()
	{
		fields |= FIELD_NM_PER_PIXEL;
	}

	private void setHasGain()
	{
		fields |= FIELD_GAIN;
	}

	private void setHasExposureTime()
	{
		fields |= FIELD_EXPOSURE_TIME;
	}

	private void setHasReadNoise()
	{
		fields |= FIELD_READ_NOISE;
	}

	private void setHasBias()
	{
		fields |= FIELD_BIAS;
	}

	private void setHasEMCCD()
	{
		fields |= FIELD_EM_CCD;
	}

	private void setHasAmplification()
	{
		fields |= FIELD_AMPLIFICATION;
	}

	public void clearHasNmPerPixel()
	{
		fields = fields & ~FIELD_NM_PER_PIXEL;
	}

	public void clearHasGain()
	{
		fields = fields & ~FIELD_GAIN;
	}

	public void clearHasExposureTime()
	{
		fields = fields & ~FIELD_EXPOSURE_TIME;
	}

	public void clearHasReadNoise()
	{
		fields = fields & ~FIELD_READ_NOISE;
	}

	public void clearHasBias()
	{
		fields = fields & ~FIELD_BIAS;
	}

	public void clearHasEMCCD()
	{
		fields = fields & ~FIELD_EM_CCD;
	}

	public void clearHasAmplification()
	{
		fields = fields & ~FIELD_AMPLIFICATION;
	}

	/**
	 * The image pixel size in nanometers
	 */
	private double nmPerPixel = 107;
	/**
	 * The gain (ADUs/photon). Can be used to convert the signal in Analogue-to-Digital units
	 * (ADUs) to photons.
	 */
	private double gain = 37.7;
	/**
	 * The exposure time in milliseconds per frame
	 */
	private double exposureTime;

	/**
	 * The CCD camera Gaussian read noise (in ADUs)
	 */
	private double readNoise;
	/**
	 * The CCD camera bias (in ADUs)
	 */
	private double bias;
	/**
	 * True if the CCD camera was run in Electron Multiplying (EM) mode
	 */
	private boolean emCCD;
	/**
	 * The camera amplification (ADUs/e-) used when modelling a microscope camera.
	 * <p>
	 * Note that the camera noise model assumes that electrons are converted to ADUs by amplification that is not
	 * perfect (i.e. it has noise). The amplification is equal to the gain (ADUs/photon) divided by the quantum
	 * efficiency (e-/photon).
	 */
	private double amplification;

	/**
	 * Default constructor. All properties are set to be invalid but missing exception are disabled.
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
	 * Parameterised constructor.
	 * <p>
	 * If gain is zero then set to 1.
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
	 * false if the current value is not valid.
	 */
	public void validate()
	{
		setNmPerPixel(nmPerPixel);
		setGain(gain);
		setExposureTime(exposureTime);
		setReadNoise(readNoise);
		setBias(bias);
		setEmCCD(emCCD);
		setAmplification(amplification);
	}

	/**
	 * Gets image pixel size in nanometers.
	 *
	 * @return image pixel size in nanometers
	 * @throws IllegalStateException if the missing field exceptions is enabled and the field has not been set to a valid value
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
	 * @throws IllegalStateException if the missing field exceptions is enabled and the field has not been set to a valid value
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
	 * @throws IllegalStateException if the missing field exceptions is enabled and the field has not been set to a valid value
	 */
	public double getExposureTime()
	{
		if (isFieldMissingException() && !hasExposureTime())
			throw new IllegalStateException();
		return exposureTime;
	}

	/**
	 * Sets the exposure time in milliseconds per frame
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
	 * Gets the CCD camera Gaussian read noise (in ADUs).
	 *
	 * @return the CCD camera Gaussian read noise (in ADUs)
	 * @throws IllegalStateException if the missing field exceptions is enabled and the field has not been set to a valid value
	 */
	public double getReadNoise()
	{
		if (isFieldMissingException() && !hasReadNoise())
			throw new IllegalStateException();
		return readNoise;
	}

	/**
	 * Sets CCD camera Gaussian read noise (in ADUs).
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
	 * Gets CCD camera bias (in ADUs).
	 *
	 * @return CCD camera bias (in ADUs)
	 * @throws IllegalStateException if the missing field exceptions is enabled and the field has not been set to a valid value
	 */
	public double getBias()
	{
		if (isFieldMissingException() && !hasBias())
			throw new IllegalStateException();
		return bias;
	}

	/**
	 * Sets CCD camera bias (in ADUs).
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
	 * Checks if the CCD camera was run in Electron Multiplying (EM) mode.
	 *
	 * @return true, if the CCD camera was run in Electron Multiplying (EM) mode.
	 * @throws IllegalStateException if the missing field exceptions is enabled and the field has not been set
	 */
	public boolean isEmCCD()
	{
		if (isFieldMissingException() && !hasEMCCD())
			throw new IllegalStateException();
		return emCCD;
	}

	/**
	 * Set if the CCD camera was run in Electron Multiplying (EM) mode.
	 *
	 * @param emCCD
	 *            true, if the CCD camera was run in Electron Multiplying (EM) mode.
	 */
	public void setEmCCD(boolean emCCD)
	{
		setHasEMCCD();
		this.emCCD = emCCD;
	}

	/**
	 * Get the camera amplification (ADUs/e-) used when modelling a microscope camera.
	 * <p>
	 * Note that the camera noise model assumes that electrons are converted to ADUs by amplification that is not
	 * perfect (i.e. it has noise). The amplification is equal to the gain (ADUs/photon) divided by the quantum
	 * efficiency (e-/photon).
	 *
	 * @return the amplification
	 * @throws IllegalStateException if the missing field exceptions is enabled and the field has not been set to a valid value
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
}
