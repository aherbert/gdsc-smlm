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
}
