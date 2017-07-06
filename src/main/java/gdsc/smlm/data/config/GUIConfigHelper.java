package gdsc.smlm.data.config;

import gdsc.smlm.data.config.GUIConfig.*;

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
 * Contains helper functions for the GUIConfig class.
 */
public class GUIConfigHelper
{
	/** The default GUIFilterSettings */
	public static final GUIFilterSettings defaultGUIFilterSettings = GUIFilterSettings.getDefaultInstance();

	/** The default PSFCalculatorSettings */
	public static final PSFCalculatorSettings defaultPSFCalculatorSettings;
	static
	{
		PSFCalculatorSettings.Builder builder = PSFCalculatorSettings.newBuilder();
		builder.setPixelPitch(6.45);
		builder.setMagnification(63);
		builder.setBeamExpander(1);
		builder.setWavelength(500);
		builder.setNumericalAperture(1.4);
		builder.setAdjustForSquarePixels(true);
		builder.setProportionalityFactor(1.52);
		defaultPSFCalculatorSettings = builder.build();
	}

	/** The default PSFEstimatorSettings */
	public static final PSFEstimatorSettings defaultPSFEstimatorSettings;
	static
	{
		PSFEstimatorSettings.Builder builder = PSFEstimatorSettings.newBuilder();
		builder.setNumberOfPeaks(1000);
		builder.setPValue(0.01);
		builder.setUpdatePreferences(true);
		builder.setDebugPsfEstimator(false);
		builder.setIterate(true);
		builder.setShowHistograms(false);
		builder.setHistogramBins(100);
		defaultPSFEstimatorSettings = builder.build();
	}

	/** The default CreateDataSettings */
	public static final CreateDataSettings defaultCreateDataSettings;
	static
	{
		CreateDataSettings.Builder builder = CreateDataSettings.newBuilder();
		builder.setSize(512);
		builder.setDepth(3000);
		builder.setSeconds(100);
		builder.setExposureTime(100);
		builder.setStepsPerSecond(10);
		builder.setDistributionMaskSliceDepth(25);
		builder.setPoissonNoise(true);
		builder.setBackground(1);
		builder.setEmGain(255);
		builder.setCameraGain(0.1557);
		builder.setQuantumEfficiency(0.95);
		builder.setReadNoise(46);
		builder.setBias(500);
		builder.setParticles(300);
		builder.setSamplePerFrame(true);
		builder.setPhotonsPerSecond(1000);
		builder.setPhotonsPerSecondMaximum(2000);
		builder.setPhotonShape(2.5);
		builder.setCorrelation(-0.35);
		builder.setWavelength(561);
		builder.setNumericalAperture(1.4);
		builder.setPsfSd(130);
		builder.setPixelPitch(107);
		builder.setDensity(1);
		builder.setConfinementMaskSliceDepth(25);
		builder.setConfinementRadius(10);
		builder.setPulseRatio(100);
		builder.setTOn(40);
		builder.setTOffShort(25);
		builder.setTOffLong(631);
		builder.setNBlinksShort(6.3);
		builder.setNBlinksLong(1.8);
		builder.setMinPhotons(30);
		builder.setCellSize(32);
		builder.setProbabilityBinary(0.1);
		builder.setMaxBinaryDistance(30);
		builder.setDensityRadius(3);
		builder.setDepthOfField(250);
		defaultCreateDataSettings = builder.build();
	}

	/** The default LoadLocalisationsSettings */
	public static final LoadLocalisationsSettings defaultLoadLocalisationsSettings = LoadLocalisationsSettings
			.getDefaultInstance();
}
