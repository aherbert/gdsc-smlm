package gdsc.smlm.data.config;

import gdsc.smlm.data.config.GUIProtos.CreateDataSettingsOrBuilder;
import gdsc.smlm.model.DiffusionType;

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
 * Contain the helper functionality for the CreateDataSettings
 */
public class CreateDataSettingsHelper
{
	CreateDataSettingsOrBuilder createDataSettings;

	private double totalGain = 0;

	/**
	 * Instantiates a new creates the data createDataSettings helper.
	 *
	 * @param createDataSettings
	 *            the create data createDataSettings
	 */
	public CreateDataSettingsHelper(CreateDataSettingsOrBuilder createDataSettings)
	{
		if (createDataSettings == null)
			throw new IllegalArgumentException("CreateDataSettings must not be null");
		if (!CalibrationProtosHelper.isCCDCameraType(createDataSettings.getCameraType()))
			throw new IllegalArgumentException("Helper instance must be used for a CCD-type camera");
		this.createDataSettings = createDataSettings;
		updateTotalGain();
	}

	/**
	 * Get the amplification (Count/electron). This is equal to the EM-gain multiplied by the camera gain. If either
	 * gain or QE are disabled then they will be ignored. An amplification of zero means no amplification is applied.
	 * 
	 * @return
	 */
	public double getAmplification()
	{
		return totalGain / getQuantumEfficiency();
	}

	/**
	 * EM-gain cannot be below 1. If so it is set to zero and disabled.
	 * 
	 * @return the emGain
	 */
	public double getEmGain()
	{
		double emGain = createDataSettings.getEmGain();
		if (emGain < 1)
			emGain = 0;
		return emGain;
	}

	/**
	 * Camera gain cannot be below 0. If so it is set to zero and disabled.
	 * 
	 * @return the cameraGain (Count/electron)
	 */
	public double getCameraGain()
	{
		double cameraGain = createDataSettings.getCameraGain();
		if (cameraGain < 0)
			cameraGain = 0;
		return cameraGain;
	}

	/**
	 * QE cannot be below 0 or above 1. If so it is set to one and disabled.
	 * 
	 * @return the quantumEfficiency
	 */
	public double getQuantumEfficiency()
	{
		double quantumEfficiency = createDataSettings.getQuantumEfficiency();
		if (quantumEfficiency < 0 || quantumEfficiency > 1)
			quantumEfficiency = 1;
		return quantumEfficiency;
	}

	/**
	 * Set the total gain (the EM-gain multiplied by the camera gain multiplied by the quantum efficiency). If either
	 * gain is disabled then it will be ignored. A total gain of zero means no gain is applied.
	 */
	private void updateTotalGain()
	{
		totalGain = 0;
		double emGain = getEmGain();
		double cameraGain = getCameraGain();
		if (cameraGain > 0)
		{
			totalGain = (emGain > 0) ? emGain * cameraGain : cameraGain;
		}
		else if (emGain > 0)
			totalGain = emGain;
		totalGain *= getQuantumEfficiency();
	}

	/**
	 * Get the total gain (Count/Photon). This is equal to the EM-gain multiplied by the camera gain multiplied by the
	 * quantum efficiency. If either gain is disabled then they will be ignored. A total gain of zero means no gain is
	 * applied.
	 * 
	 * @return
	 */
	public double getTotalGain()
	{
		return totalGain;
	}

	/**
	 * Get the total gain (Count/Photon). This is equal to the EM-gain multiplied by the camera gain multiplied by the
	 * quantum efficiency. If either gain is disabled then they will be ignored. If the total gain is not above zero
	 * (due to invalid configuration) then 1 will be returned.
	 * 
	 * @return
	 */
	public double getTotalGainSafe()
	{
		return (totalGain > 0) ? totalGain : 1;
	}

	/**
	 * Gets the diffusion type.
	 *
	 * @param diffusionType
	 *            the diffusion type
	 * @return the diffusion type
	 */
	public static DiffusionType getDiffusionType(int diffusionType)
	{
		if (diffusionType >= 0 && diffusionType < DiffusionType.values().length)
		{
			return DiffusionType.values()[diffusionType];
		}
		// Set a default
		return DiffusionType.RANDOM_WALK;
	}

	/**
	 * Gets the read noise in Counts. This is the read noise in electrons multiplied by the camera gain.
	 *
	 * @return the read noise in Counts
	 */
	public double getReadNoiseInCounts()
	{
		double readNoise = 0;
		if (createDataSettings.getReadNoise() > 0)
		{
			readNoise = createDataSettings.getReadNoise();
			// Read noise is in electrons. Apply camera gain to get the noise in ADUs.
			if (createDataSettings.getCameraGain() != 0)
				readNoise *= createDataSettings.getCameraGain();
		}
		return readNoise;
	}
}
