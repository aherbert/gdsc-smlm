/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.ij.settings;

import java.io.File;
import uk.ac.sussex.gdsc.core.clustering.ClusteringAlgorithm;
import uk.ac.sussex.gdsc.core.clustering.optics.SampleMode;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper.LutColour;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraType;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FilterSettings;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FitEngineSettings;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.PrecisionMethod;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.data.config.PsfProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageSettings;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageSizeMode;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageType;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsTableSettings;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.TimeUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.Optics.ClusteringMode;
import uk.ac.sussex.gdsc.smlm.ij.plugins.Optics.ImageMode;
import uk.ac.sussex.gdsc.smlm.ij.plugins.Optics.OpticsMode;
import uk.ac.sussex.gdsc.smlm.ij.plugins.Optics.OutlineMode;
import uk.ac.sussex.gdsc.smlm.ij.plugins.Optics.PlotMode;
import uk.ac.sussex.gdsc.smlm.ij.plugins.Optics.SpanningTreeMode;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.AstigmatismModelManagerSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.CameraModelAnalysisSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.CameraModelFisherInformationAnalysisSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.CameraModelManagerSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.ClusteringSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.ConfigurationTemplateSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.CreateDataSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.CubicSplineManagerSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.FailCountManagerSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.GUIFilterSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.ImageJ3DResultsViewerSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.LoadLocalisationsSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.NucleusMaskSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.OpticsEventSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.OpticsSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.PSFCalculatorSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.PSFCreatorSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.PSFEstimatorSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.SpotFitSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.TcPalmAnalysisSettings;
import uk.ac.sussex.gdsc.smlm.results.DynamicMultipleTargetTracing;
import uk.ac.sussex.gdsc.smlm.results.TraceManager.TraceMode;

/**
 * Contains default graphic user interface (GUI) settings.
 */
public final class GuiSettings {
  /** The default GUIFilterSettings. */
  public static class DefaultGUIFilterSettings {
    /** Default settings instance. */
    public static final GUIFilterSettings INSTANCE;

    static {
      final GUIFilterSettings.Builder builder = GUIFilterSettings.newBuilder();
      INSTANCE = builder.build();
    }
  }

  /** The default PSFCalculatorSettings. */
  public static class DefaultPSFCalculatorSettings {
    /** Default settings instance. */
    public static final PSFCalculatorSettings INSTANCE;

    static {
      final PSFCalculatorSettings.Builder builder = PSFCalculatorSettings.newBuilder();
      builder.setPixelPitch(6.45);
      builder.setMagnification(63);
      builder.setBeamExpander(1);
      builder.setWavelength(500);
      builder.setNumericalAperture(1.4);
      builder.setAdjustForSquarePixels(true);
      builder.setProportionalityFactor(1.52);
      INSTANCE = builder.build();
    }
  }

  /** The default PSFEstimatorSettings. */
  public static class DefaultPSFEstimatorSettings {
    /** Default settings instance. */
    public static final PSFEstimatorSettings INSTANCE;

    static {
      final PSFEstimatorSettings.Builder builder = PSFEstimatorSettings.newBuilder();
      builder.setNumberOfPeaks(1000);
      builder.setPValue(0.01);
      builder.setUpdatePreferences(true);
      builder.setDebugPsfEstimator(false);
      builder.setIterate(true);
      builder.setShowHistograms(false);
      builder.setHistogramBins(100);
      INSTANCE = builder.build();
    }
  }

  /** The default CreateDataSettings. */
  public static class DefaultCreateDataSettings {
    /** Default settings instance. */
    public static final CreateDataSettings INSTANCE;

    static {
      final CreateDataSettings.Builder builder = CreateDataSettings.newBuilder();
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
      INSTANCE = builder.build();
    }
  }

  /** The default LoadLocalisationsSettings. */
  public static class DefaultLoadLocalisationsSettings {
    /** Default settings instance. */
    public static final LoadLocalisationsSettings INSTANCE;

    static {
      final LoadLocalisationsSettings.Builder builder = LoadLocalisationsSettings.newBuilder();
      builder.setFieldT(0);
      builder.setFieldId(-1);
      builder.setFieldCategory(-1);
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
      INSTANCE = builder.build();
    }
  }
  /** The default ClusteringSettings. */
  public static class DefaultClusteringSettings {
    /** Default settings instance. */
    public static final ClusteringSettings INSTANCE;

    static {
      final ClusteringSettings.Builder builder = ClusteringSettings.newBuilder();
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
      // DynamicMultipleTargetTracing
      final DynamicMultipleTargetTracing.DmttConfiguration config =
          DynamicMultipleTargetTracing.DmttConfiguration.newBuilder(1).build();
      builder.setTemporalWindow(config.getTemporalWindow());
      builder.setLocalDiffusionWeight(config.getLocalDiffusionWeight());
      builder.setOnIntensityWeight(config.getLocalDiffusionWeight());
      builder.setDisappearanceDecayFactor(config.getDisappearanceDecayFactor());
      builder.setDisappearanceThreshold(config.getDisappearanceThreshold());
      INSTANCE = builder.build();
    }
  }

  /** The default OpticsSettings. */
  public static class DefaultOpticsSettings {
    /** Default settings instance. */
    public static final OpticsSettings INSTANCE;

    static {
      final OpticsSettings.Builder builder = OpticsSettings.newBuilder();
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
      builder.setDiggingThreshold(2.0);
      final OpticsEventSettings.Builder b = builder.getOpticsEventSettingsBuilder();
      b.setShowSelectionTable(true);
      b.setTableCreateSelection(true);
      b.setTableShowSelection(true);
      b.setImageCreateSelection(true);
      b.setImageShowSelection(true);
      b.setPlotCreateSelection(true);
      b.setPlotShowSelection(true);
      INSTANCE = builder.build();
    }
  }

  /** The default ConfigurationTemplateSettings. */
  public static class DefaultConfigurationTemplateSettings {
    /** Default settings instance. */
    public static final ConfigurationTemplateSettings INSTANCE;

    static {
      final ConfigurationTemplateSettings.Builder builder =
          ConfigurationTemplateSettings.newBuilder();
      builder.setSelectStandardTemplates(true);
      builder.setSelectCustomDirectory(false);
      builder.setConfigurationDirectory(
          System.getProperty("user.home") + File.separator + "uk.ac.sussex.gdsc.smlm");
      INSTANCE = builder.build();
    }
  }

  /** The default NucleusMaskSettings. */
  public static class DefaultNucleusMaskSettings {
    /** Default settings instance. */
    public static final NucleusMaskSettings INSTANCE;

    static {
      final NucleusMaskSettings.Builder builder = NucleusMaskSettings.newBuilder();
      builder.setMode(1);
      builder.setFieldWidth(512);
      builder.setYDither(4);
      builder.setZDither(1);
      builder.setNmPerPixel(100);
      builder.setNmPerSlice(20);
      builder.setDiameter(2);
      INSTANCE = builder.build();
    }
  }

  /** The default PSFCreatorSettings. */
  public static class DefaultPSFCreatorSettings {
    /** Default settings instance. */
    public static final PSFCreatorSettings INSTANCE;

    static {
      final PSFCreatorSettings.Builder builder = PSFCreatorSettings.newBuilder();
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
      // Use 2D Projections (3D cross correlation has normalisation stability issues)
      builder.setAlignmentMode(0);
      builder.setAlignmentZRadius(0); // All of the z-stack
      builder.setSubPixelPrecision(0.01);
      builder.setRmsdXyThreshold(0.01);
      builder.setRmsdZThreshold(0.05);
      builder.setComShiftThreshold(0.01);
      // Copy default fit settings but modify for fitting standalone PSFs
      final FitEngineSettings.Builder febuilder =
          FitProtosHelper.DefaultFitEngineSettings.INSTANCE.toBuilder();
      febuilder.setIncludeNeighbours(false);
      builder.setFitEngineSettings(febuilder);
      builder.setPsf(PsfProtosHelper.getDefaultPsf(PSFType.TWO_AXIS_GAUSSIAN_2D));
      INSTANCE = builder.build();
    }
  }

  /** The default CameraModelManagerSettings. */
  public static class DefaultCameraModelManagerSettings {
    /** Default settings instance. */
    public static final CameraModelManagerSettings INSTANCE;

    static {
      final CameraModelManagerSettings.Builder builder = CameraModelManagerSettings.newBuilder();
      INSTANCE = builder.build();
    }
  }

  /** The default CameraModelAnalysisSettings. */
  public static class DefaultCameraModelAnalysisSettings {
    /** Default settings instance. */
    public static final CameraModelAnalysisSettings INSTANCE;

    static {
      final CameraModelAnalysisSettings.Builder builder = CameraModelAnalysisSettings.newBuilder();
      builder.setPhotons(10);
      // Note that the total gain is likely to be very different if using EM-CCD/CCD/sCMOS
      // so these are separate.
      // Use counts as it is simpler to understand (and measure) than electrons on manufacturers
      // specification sheets.

      // Need a better CCD estimate. This is Photometrics CoolSNAP HQ2
      builder.setGain(1); // Count/electron
      builder.setNoise(5.5); // Count

      // EM-CCD: Photometrics Evolve 512
      builder.setEmGain(40); // This is the total gain in Count/electron
      builder.setEmNoise(13); // Count

      // sCMOS: Photometrics Prime95b
      builder.setCmosGain(1.7); // Count/electron
      builder.setCmosNoise(3.4); // Count

      builder.setSamples(20000);
      builder.setNoiseSamples(10);
      builder.setEmSamples(5);
      INSTANCE = builder.build();
    }
  }

  /** The default CameraModelFisherInformationAnalysisSettings. */
  public static class DefaultCameraModelFisherInformationAnalysisSettings {
    /** Default settings instance. */
    public static final CameraModelFisherInformationAnalysisSettings INSTANCE;

    static {
      final CameraModelFisherInformationAnalysisSettings.Builder builder =
          CameraModelFisherInformationAnalysisSettings.newBuilder();
      builder.setMinExponent(-6);
      builder.setMaxExponent(2);
      builder.setSubDivisions(1);
      builder.setCamera2Type(1);
      builder.setCamera1Gain(1);
      builder.setCamera1Noise(4);
      builder.setCamera2Type(3);
      builder.setCamera2Gain(20);
      builder.setCamera2Noise(8);
      INSTANCE = builder.build();
    }
  }
  /** The default CubicSplineManagerSettings. */
  public static class DefaultCubicSplineManagerSettings {
    /** Default settings instance. */
    public static final CubicSplineManagerSettings INSTANCE;

    static {
      final CubicSplineManagerSettings.Builder builder = CubicSplineManagerSettings.newBuilder();
      builder.setMagnification(3);
      builder.setScale(2);
      INSTANCE = builder.build();
    }
  }

  /** The default FailCountManagerSettings. */
  public static class DefaultFailCountManagerSettings {
    /** Default settings instance. */
    public static final FailCountManagerSettings INSTANCE;

    static {
      final FailCountManagerSettings.Builder builder = FailCountManagerSettings.newBuilder();
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
      INSTANCE = builder.build();
    }
  }

  /** The default AstigmatismModelManagerSettings. */
  public static class DefaultAstigmatismModelManagerSettings {
    /** Default settings instance. */
    public static final AstigmatismModelManagerSettings INSTANCE;

    static {
      final AstigmatismModelManagerSettings.Builder builder =
          AstigmatismModelManagerSettings.newBuilder();
      builder.setSmoothing(0.2);
      builder.setWeightedFit(true);
      builder.setSaveFitWidth(true);
      builder.setSaveModel(true);
      final FitEngineSettings.Builder b =
          FitProtosHelper.DefaultFitEngineSettings.INSTANCE.toBuilder();

      // Adjust for a wider fit range
      b.getFittingBuilder().setValue(10).setAbsolute(true);

      // Simple filter
      final FilterSettings.Builder fb = b.getFitSettingsBuilder().getFilterSettingsBuilder();
      fb.setSmartFilter(false);
      fb.setDisableSimpleFilter(false);
      fb.setShiftFactor(2);
      fb.setSignalStrength(5);
      fb.setMinPhotons(50);
      fb.setMinWidthFactor(0.5);
      fb.setMaxWidthFactor(5);
      fb.setPrecisionThreshold(50);
      fb.setPrecisionMethodValue(PrecisionMethod.POISSON_CRLB_VALUE);

      builder.setFitEngineSettings(b);
      builder.setPsf(PsfProtosHelper.DefaultTwoAxisGaussian2dPsf.INSTANCE);

      INSTANCE = builder.build();
    }
  }

  /** The default ImageJ3DResultsViewerSettings. */
  public static class DefaultImageJ3DResultsViewerSettings {
    /** Default settings instance. */
    public static final ImageJ3DResultsViewerSettings INSTANCE;

    static {
      final ImageJ3DResultsViewerSettings.Builder builder =
          ImageJ3DResultsViewerSettings.newBuilder();
      builder.setNewWindow(false);
      builder.setTransparency(0.65);
      builder.setLut(LutColour.FIRE_GLOW.ordinal());
      builder.setRendering(7); // Octahedron
      builder.setPixelSize(2); // Like weighted 2D pixel rendering using 4 pixels
      builder.setShaded(true);
      builder.setSizeMode(0); // Fixed
      builder.setSize(10); // NM units
      builder.setSortMode(0); // None
      builder.setTransparencyMode(1); // From the size
      builder.setMinTransparency(0);
      builder.setMaxTransparency(0.95);
      builder.setDepthMode(1);
      builder.setDepthRange(20);
      builder.setDitherSeed(123456789);
      builder.setNameOption(1);
      builder.setNameSuffix(" Cropped");
      final ResultsTableSettings.Builder resultsTableSettings =
          builder.getResultsTableSettingsBuilder();
      resultsTableSettings.setDistanceUnit(DistanceUnit.NM);
      resultsTableSettings.setShowTable(true);
      INSTANCE = builder.build();
    }
  }

  /** The default SpotFitSettings. */
  public static class DefaultSpotFitSettings {
    /** Default settings instance. */
    public static final SpotFitSettings INSTANCE;

    static {
      final SpotFitSettings.Builder builder = SpotFitSettings.newBuilder();
      builder.setChannel(1);
      builder.setSearchRadius(3);
      builder.setFitRadius(10);
      builder.setShowOverlay(true);
      builder.setAnalysisRadius(5);
      INSTANCE = builder.build();
    }
  }

  /** The default TcPalmAnalysisSettings. */
  public static class DefaultTcPalmAnalysisSettings {
    /** Default settings instance. */
    public static final TcPalmAnalysisSettings INSTANCE;

    static {
      final TcPalmAnalysisSettings.Builder builder = TcPalmAnalysisSettings.newBuilder();
      ResultsImageSettings.Builder resultsImageSettings = builder.getResultsImageSettingsBuilder();
      resultsImageSettings.setImageType(ResultsImageType.DRAW_LOCALISATIONS);
      resultsImageSettings.setLutName(LutColour.FIRE.getName());
      builder.setLoopSize(512);
      // @formatter:off
      builder.getLoopImageSettingsBuilder()
        .setImageType(ResultsImageType.DRAW_LOCALISATIONS)
        .setImageSizeMode(ResultsImageSizeMode.IMAGE_SIZE)
        .setScale(10)
        .setImageSize(512)
        .setPixelSize(5);
      // @formatter:on
      // No LUT to use the grey default. This allows a colour overlay to be distinct.
      INSTANCE = builder.build();
    }
  }

  /** No public constructor. */
  private GuiSettings() {}
}
