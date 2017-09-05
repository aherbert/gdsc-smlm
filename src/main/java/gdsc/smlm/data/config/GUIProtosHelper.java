package gdsc.smlm.data.config;

import java.io.File;

import gdsc.core.clustering.ClusteringAlgorithm;
import gdsc.core.clustering.optics.SampleMode;
import gdsc.smlm.data.config.GUIProtos.*;
import gdsc.smlm.data.config.UnitProtos.TimeUnit;
import gdsc.smlm.ij.plugins.OPTICS.ClusteringMode;
import gdsc.smlm.ij.plugins.OPTICS.ImageMode;
import gdsc.smlm.ij.plugins.OPTICS.OpticsMode;
import gdsc.smlm.ij.plugins.OPTICS.OutlineMode;
import gdsc.smlm.ij.plugins.OPTICS.PlotMode;
import gdsc.smlm.ij.plugins.OPTICS.SpanningTreeMode;
import gdsc.smlm.results.TraceManager.TraceMode;

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
 * Contains helper functions for the GUIProtos class.
 */
public class GUIProtosHelper
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
	
	/** The default ClusteringSettings */
	public static final ClusteringSettings defaultClusteringSettings;
	static
	{
		ClusteringSettings.Builder builder = ClusteringSettings.newBuilder();
		builder.setDistanceThreshold(50);
		builder.setDistanceExclusion(0);
		builder.setTimeThreshold(1);
		builder.setTimeUnit(TimeUnit.SECOND);
		builder.setTraceMode(TraceMode.LATEST_FORERUNNER.ordinal());
		builder.setClusteringAlgorithm(ClusteringAlgorithm.PAIRWISE.ordinal());
		builder.setBlinkingRate(1);
		builder.setMaxDistanceThreshold(500);
		builder.setMaxTimeThreshold(20);
		builder.setOptimiserSteps(10);
		builder.setOptimiserPlot(2);
		builder.setMinimumTraceLength(6);
		builder.setInternalDistances(true);
		builder.setIgnoreEnds(true);
		builder.setPrecisionCorrection(true);
		builder.setMsdCorrection(true);
		builder.setMle(true);
		builder.setFitLength(6);
		builder.setFitRestarts(1);
		builder.setJumpDistance(1);
		defaultClusteringSettings = builder.build();
	}
	
	/** The default OpticsSettings */
	public static final OpticsSettings defaultOpticsSettings;
	static
	{
		OpticsSettings.Builder builder = OpticsSettings.newBuilder();
		builder.setOpticsMode(OpticsMode.FAST_OPTICS.ordinal());
		builder.setSampleMode(SampleMode.RANDOM.ordinal());
		builder.setMinPoints(4);
		builder.setClusteringMode(ClusteringMode.XI.ordinal());
		builder.setXi(0.03);
		builder.setSamples(100);
		builder.setSampleFraction(0.05);
		builder.setFractionNoise(0.05);
		builder.setImageScale(2);
		builder.setImageMode(ImageMode.VALUE.ordinal());
		builder.setWeighted(true);
		builder.setEqualised(true);
		builder.setPlotMode(PlotMode.COLOURED_BY_DEPTH_WITH_CLUSTERS.ordinal());
		builder.setOutlineMode(OutlineMode.COLOURED_BY_CLUSTER.ordinal());
		builder.setSpanningTreeMode(SpanningTreeMode.OFF.ordinal());
		builder.setLambda(3);
		defaultOpticsSettings = builder.build();
	}
	
	/** The default ConfigurationTemplateSettings */
	public static final ConfigurationTemplateSettings defaultConfigurationTemplateSettings;
	static
	{
		ConfigurationTemplateSettings.Builder builder = ConfigurationTemplateSettings.newBuilder();
		builder.setSelectStandardTemplates(true);
		builder.setSelectCustomDirectory(false);
		builder.setConfigurationDirectory(System.getProperty("user.home") + File.separator + "gdsc.smlm");
		defaultConfigurationTemplateSettings = builder.build();
	}
}
