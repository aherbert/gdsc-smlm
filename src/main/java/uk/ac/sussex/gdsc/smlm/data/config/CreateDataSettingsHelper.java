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

import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraType;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.CreateDataSettingsOrBuilder;
import uk.ac.sussex.gdsc.smlm.model.DiffusionType;

/**
 * Contain the helper functionality for the CreateDataSettings
 */
public class CreateDataSettingsHelper
{
	/** The create data settings. */
	CreateDataSettingsOrBuilder createDataSettings;

	/** Set to true if the camera type is {@link CameraType#EMCCD} */
	final public boolean isEMCCD;
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
		isEMCCD = createDataSettings.getCameraTypeValue() == CameraType.EMCCD_VALUE;
		updateTotalGain();
	}

	/**
	 * Get the amplification (Count/electron). This is equal to the EM-gain multiplied by the camera gain. If either
	 * gain or QE are disabled then they will be ignored. An amplification of zero means no amplification is applied.
	 *
	 * @return the amplification
	 */
	public double getAmplification()
	{
		return totalGain / getQuantumEfficiency();
	}

	/**
	 * EM-gain cannot be below 1. If so it is set to zero and disabled.
	 * <p>
	 * This is also zero for a non EM-CCD camera.
	 *
	 * @return the emGain
	 */
	public double getEmGain()
	{
		if (isEMCCD)
		{
			double emGain = createDataSettings.getEmGain();
			if (emGain < 1)
				emGain = 0;
			return emGain;
		}
		return 0;
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
		final double emGain = getEmGain();
		final double cameraGain = getCameraGain();
		if (cameraGain > 0)
			totalGain = (emGain > 0) ? emGain * cameraGain : cameraGain;
		else if (emGain > 0)
			totalGain = emGain;
		totalGain *= getQuantumEfficiency();
	}

	/**
	 * Get the total gain (Count/Photon). This is equal to the EM-gain multiplied by the camera gain multiplied by the
	 * quantum efficiency. If either gain is disabled then they will be ignored. A total gain of zero means no gain is
	 * applied.
	 *
	 * @return the total gain
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
	 * @return the total gain safe
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
			return DiffusionType.values()[diffusionType];
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
