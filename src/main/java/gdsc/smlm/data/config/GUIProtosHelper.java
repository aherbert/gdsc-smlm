package gdsc.smlm.data.config;

import java.io.File;

import gdsc.core.clustering.ClusteringAlgorithm;
import gdsc.core.clustering.optics.SampleMode;
import gdsc.smlm.data.config.CalibrationProtos.CameraType;
import gdsc.smlm.data.config.FitProtos.FilterSettings;
import gdsc.smlm.data.config.FitProtos.FitEngineSettings;
import gdsc.smlm.data.config.FitProtos.PrecisionMethod;
import gdsc.smlm.data.config.GUIProtos.CameraModelManagerSettings;
import gdsc.smlm.data.config.GUIProtos.ClusteringSettings;
import gdsc.smlm.data.config.GUIProtos.ConfigurationTemplateSettings;
import gdsc.smlm.data.config.GUIProtos.CreateDataSettings;
import gdsc.smlm.data.config.GUIProtos.CubicSplineManagerSettings;
import gdsc.smlm.data.config.GUIProtos.FailCountManagerSettings;
import gdsc.smlm.data.config.GUIProtos.GUIFilterSettings;
import gdsc.smlm.data.config.GUIProtos.Image3DDrawingMode;
import gdsc.smlm.data.config.GUIProtos.ImageJ3DResultsViewerSettings;
import gdsc.smlm.data.config.GUIProtos.LoadLocalisationsSettings;
import gdsc.smlm.data.config.GUIProtos.NucleusMaskSettings;
import gdsc.smlm.data.config.GUIProtos.OpticsEventSettings;
import gdsc.smlm.data.config.GUIProtos.OpticsSettings;
import gdsc.smlm.data.config.GUIProtos.AstigmatismModelManagerSettings;
import gdsc.smlm.data.config.GUIProtos.PSFCalculatorSettings;
import gdsc.smlm.data.config.GUIProtos.PSFCreatorSettings;
import gdsc.smlm.data.config.GUIProtos.PSFEstimatorSettings;
import gdsc.smlm.data.config.UnitProtos.TimeUnit;
import gdsc.smlm.ij.plugins.OPTICS.ClusteringMode;
import gdsc.smlm.ij.plugins.OPTICS.ImageMode;
import gdsc.smlm.ij.plugins.OPTICS.OpticsMode;
import gdsc.smlm.ij.plugins.OPTICS.OutlineMode;
import gdsc.smlm.ij.plugins.OPTICS.PlotMode;
import gdsc.smlm.ij.plugins.OPTICS.SpanningTreeMode;
import gdsc.smlm.results.TraceManager.TraceMode;
import ij.process.LUTHelper.LutColour;

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
		builder.setCameraType(CameraType.EMCCD);
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
		builder.setDepthOfFocus(450);
		defaultCreateDataSettings = builder.build();
	}

	/** The default LoadLocalisationsSettings */
	public static final LoadLocalisationsSettings defaultLoadLocalisationsSettings;
	static {
		LoadLocalisationsSettings.Builder builder = LoadLocalisationsSettings.newBuilder();
		builder.setFieldT(0);
		builder.setFieldId(-1);
		builder.setFieldX(1);
		builder.setFieldY(2);
		builder.setFieldZ(-1);
		builder.setFieldI(3);
		builder.setFieldSx(-1);
		builder.setFieldSy(-1);
		builder.setFieldPrecision(-1);
		builder.setHeaderLines(1);
		builder.setComment("#");
		builder.setDelimiter("\\s+");
		builder.setName("Localisations");
		defaultLoadLocalisationsSettings = builder.build();
	}

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
		OpticsEventSettings.Builder b = builder.getOpticsEventSettingsBuilder();
		b.setShowSelectionTable(true);
		b.setTableCreateSelection(true);
		b.setTableShowSelection(true);
		b.setImageCreateSelection(true);
		b.setImageShowSelection(true);
		b.setPlotCreateSelection(true);
		b.setPlotShowSelection(true);
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

	/** The default NucleusMaskSettings */
	public static final NucleusMaskSettings defaultNucleusMaskSettings;
	static
	{
		NucleusMaskSettings.Builder builder = NucleusMaskSettings.newBuilder();
		builder.setMode(1);
		builder.setFieldWidth(512);
		builder.setYDither(4);
		builder.setZDither(1);
		builder.setNmPerPixel(100);
		builder.setNmPerSlice(20);
		builder.setDiameter(2);
		defaultNucleusMaskSettings = builder.build();
	}

	/** The default PSFCreatorSettings */
	public static final PSFCreatorSettings defaultPSFCreatorSettings;
	static
	{
		PSFCreatorSettings.Builder builder = PSFCreatorSettings.newBuilder();
		builder.setRadius(10);
		builder.setAmplitudeFraction(0.2);
		builder.setStartBackgroundFrames(5);
		builder.setEndBackgroundFrames(5);
		builder.setMagnification(10);
		builder.setSmoothing(0.1);
		builder.setCentreEachSlice(false);
		builder.setComCutOff(5e-2);
		builder.setInteractiveMode(false);
		builder.setInterpolationMethod(2); // ImageProcessor.BICUBIC
		builder.setPsfType(0);
		builder.setAnalysisWindow(1);
		builder.setComWindow(2);
		builder.setAlignmentMagnification(2);
		builder.setMaxIterations(20);
		builder.setCheckAlignments(false);
		builder.setPsfMagnification(4);
		builder.setWindow(3);
		builder.setSmoothStackSignal(true);
		builder.setComBorder(0.40); // Right in the middle of a spot
		builder.setAlignmentMode(0); // 2D Projections (3D cross correlation has normalisations stability issues)
		builder.setAlignmentZRadius(0); // All of the z-stack
		builder.setSubPixelPrecision(0.01);
		builder.setRmsdXyThreshold(0.01);
		builder.setRmsdZThreshold(0.05);
		builder.setComShiftThreshold(0.01);
		defaultPSFCreatorSettings = builder.build();
	}

	/** The default CameraModelManagerSettings */
	public static final CameraModelManagerSettings defaultCameraModelManagerSettings;
	static
	{
		defaultCameraModelManagerSettings = CameraModelManagerSettings.getDefaultInstance();
	}

	/** The default CubicSplineManagerSettings */
	public static final CubicSplineManagerSettings defaultCubicSplineManagerSettings;
	static
	{
		CubicSplineManagerSettings.Builder builder = CubicSplineManagerSettings.newBuilder();
		builder.setMagnification(3);
		builder.setScale(2);
		defaultCubicSplineManagerSettings = builder.build();
	}

	/** The default FailCountManagerSettings */
	public static final FailCountManagerSettings defaultFailCountManagerSettings;
	static
	{
		FailCountManagerSettings.Builder builder = FailCountManagerSettings.newBuilder();
		builder.setMaxFrames(100);
		builder.setFailCountLimit(50);
		builder.setTargetPassFraction(0.9);
		builder.setPlotRollingWindow(10);
		builder.setPlotPassWeight(2);
		builder.setPlotFailWeight(1);
		builder.setPlotResetFraction(0.5);
		builder.setPlotFixedXAxis(true);
		builder.setTableTopN(1);
		builder.setRollingCounterMinAllowedFailures(1);
		builder.setRollingCounterMaxAllowedFailures(100);
		builder.setRollingCounterMinWindow(2);
		builder.setRollingCounterMaxWindow(200);
		builder.setWeightedCounterMinAllowedFailures(1);
		builder.setWeightedCounterMaxAllowedFailures(200);
		builder.setWeightedCounterMinPassDecrement(1);
		builder.setWeightedCounterMaxPassDecrement(10);
		builder.setResettingCounterMinAllowedFailures(1);
		builder.setResettingCounterMaxAllowedFailures(200);
		builder.setResettingCounterMinResetFraction(0.5);
		builder.setResettingCounterMaxResetFraction(0.95);
		builder.setResettingCounterIncResetFraction(0.05);
		builder.setPassRateCounterMinAllowedCounts(0);
		builder.setPassRateCounterMaxAllowedCounts(5);
		builder.setPassRateCounterMinPassRate(0.05);
		builder.setPassRateCounterMaxPassRate(0.2);
		builder.setPassRateCounterIncPassRate(0.01);
		defaultFailCountManagerSettings = builder.build();
	}

	/** The default AstigmatismModelManagerSettings */
	public static final AstigmatismModelManagerSettings defaultAstigmatismModelManagerSettings;
	static
	{
		AstigmatismModelManagerSettings.Builder builder = AstigmatismModelManagerSettings.newBuilder();
		builder.setSmoothing(0.2);
		builder.setWeightedFit(true);
		builder.setSaveFitWidth(true);
		builder.setSaveModel(true);
		FitEngineSettings.Builder b = FitProtosHelper.defaultFitEngineSettings.toBuilder();
		
		// Adjust for a wider fit range
		b.getFittingBuilder().setValue(10).setAbsolute(true);
		
		// Simple filter
		FilterSettings.Builder fb = b.getFitSettingsBuilder().getFilterSettingsBuilder();
		fb.setSmartFilter(false);
		fb.setDisableSimpleFilter(false);
		fb.setShiftFactor(2);
		fb.setSignalStrength(0);
		fb.setMinPhotons(50);
		fb.setMinWidthFactor(0.5);
		fb.setMaxWidthFactor(5);
		fb.setPrecisionThreshold(50);
		fb.setPrecisionMethodValue(PrecisionMethod.POISSON_CRLB_VALUE);
		
		builder.setFitEngineSettings(b);
		builder.setPsf(PSFProtosHelper.defaultTwoAxisGaussian2DPSF);

		defaultAstigmatismModelManagerSettings = builder.build();
	}
	
	public static String getName(Image3DDrawingMode value)
	{
		switch (value)
		{
			case DRAW_3D_FIXED_SIZE:
				return "Fixed size";
			case DRAW_3D_XYZ_DEVIATIONS:
				return "XYZ Deviations";
			case DRAW_3D_XY_PRECISION:
				return "XY Precision";
			case UNRECOGNIZED:
				return "Unknown";
			default:
				throw new IllegalStateException("Unknown name: " + value);
		}
	}	

	/** The default ImageJ3DResultsViewerSettings */
	public static final ImageJ3DResultsViewerSettings defaultImageJ3DResultsViewerSettings;
	static
	{
		ImageJ3DResultsViewerSettings.Builder builder = ImageJ3DResultsViewerSettings.newBuilder();
		builder.setSize(10); // NM units
		builder.setTransparency(0.5);
		builder.setLut(LutColour.FIRE.ordinal());
		builder.setRendering(4); // Octahedron
		builder.setShaded(true);
		builder.setDrawingModeValue(Image3DDrawingMode.DRAW_3D_FIXED_SIZE_VALUE);
		builder.setPixelSize(2); // Like weighted 2D pixel rendering using 4 pixels
		defaultImageJ3DResultsViewerSettings = builder.build();
	}
}
