/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.ij.plugins;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.gui.YesNoCancelDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.process.LUT;
import java.awt.Checkbox;
import java.awt.Choice;
import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.Insets;
import java.awt.Label;
import java.awt.Panel;
import java.awt.Rectangle;
import java.awt.Scrollbar;
import java.awt.SystemColor;
import java.awt.TextField;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.IntFunction;
import java.util.function.Supplier;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import org.apache.commons.lang3.StringUtils;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.ij.ImageJPluginLoggerHelper;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.SeriesOpener;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog.OptionListener;
import uk.ac.sussex.gdsc.core.ij.gui.OffsetPointRoi;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper.LutColour;
import uk.ac.sussex.gdsc.core.logging.TrackProgressAdaptor;
import uk.ac.sussex.gdsc.core.utils.BitFlagUtils;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraType;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.DataFilterMethod;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FitEngineSettings;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FitSolver;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.NoiseEstimatorMethod;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.PrecisionMethod;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.PSFCalculatorSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GuiProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.AstigmatismModel;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFParameter;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.data.config.PsfHelper;
import uk.ac.sussex.gdsc.smlm.data.config.PsfProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsFileSettings;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageSettings;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageType;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsSettings;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsTableSettings;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.TemplateProtos.TemplateSettings;
import uk.ac.sussex.gdsc.smlm.data.config.UnitConverterUtils;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.TimeUnit;
import uk.ac.sussex.gdsc.smlm.engine.FitConfiguration;
import uk.ac.sussex.gdsc.smlm.engine.FitEngine;
import uk.ac.sussex.gdsc.smlm.engine.FitEngineConfiguration;
import uk.ac.sussex.gdsc.smlm.engine.FitJob;
import uk.ac.sussex.gdsc.smlm.engine.FitParameters;
import uk.ac.sussex.gdsc.smlm.engine.FitParameters.FitTask;
import uk.ac.sussex.gdsc.smlm.engine.FitQueue;
import uk.ac.sussex.gdsc.smlm.engine.FitWorker;
import uk.ac.sussex.gdsc.smlm.engine.ParameterisedFitJob;
import uk.ac.sussex.gdsc.smlm.filters.SpotFilter;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.FastMleSteppingFunctionSolver;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.MaximumLikelihoodFitter;
import uk.ac.sussex.gdsc.smlm.ij.IJImageSource;
import uk.ac.sussex.gdsc.smlm.ij.SeriesImageSource;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.results.ImageJImagePeakResults;
import uk.ac.sussex.gdsc.smlm.ij.results.ImageJTablePeakResults;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.ij.utils.ImageJImageConverter;
import uk.ac.sussex.gdsc.smlm.model.camera.CameraModel;
import uk.ac.sussex.gdsc.smlm.model.camera.PerPixelCameraModel;
import uk.ac.sussex.gdsc.smlm.results.AggregatedImageSource;
import uk.ac.sussex.gdsc.smlm.results.FilePeakResults;
import uk.ac.sussex.gdsc.smlm.results.ImageSource;
import uk.ac.sussex.gdsc.smlm.results.InterlacedImageSource;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.PeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResultsList;
import uk.ac.sussex.gdsc.smlm.results.count.Counter;
import uk.ac.sussex.gdsc.smlm.results.count.FrameCounter;
import uk.ac.sussex.gdsc.smlm.results.filter.DirectFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.Filter;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedureX;
import uk.ac.sussex.gdsc.smlm.results.procedures.XyResultProcedure;

/**
 * Fits local maxima using a 2D Gaussian. Process each frame until a successive number of fits fail
 * to meet the fit criteria.
 */
public class PeakFit implements PlugInFilter {
  private static final String TITLE = "PeakFit";
  private static final String LOG_SPACER = "-=-=-=-";

  /** Flag to indicate that additional options can be configured. */
  public static final int FLAG_EXTRA_OPTIONS = 0x01;
  /** Flag to indicate that the calibration should not be configured. */
  public static final int FLAG_IGNORE_CALIBRATION = 0x02;
  /** Flag to indicate that configuration should not be saved. */
  public static final int FLAG_NO_SAVE = 0x04;

  private static final int FLAGS = DOES_16 | DOES_8G | DOES_32 | NO_CHANGES;

  private int pluginFlags;
  private int singleFrame;
  private ImagePlus imp;

  // Used within the run(ImageProcessor) method
  private Rectangle bounds;
  private boolean ignoreBoundsForNoise = true;
  private Logger logger;

  private ImageSource source;
  private String resultsSuffix;
  private PeakResultsList results;
  private long time;
  private long runTime;
  private FitEngineConfiguration config;
  private FitConfiguration fitConfig;
  private ResultsSettings.Builder resultsSettings;
  private boolean silent;

  /** Flag for extra options in the dialog (shift-key down). */
  private boolean extraOptions;
  /** True if running in maxima identification mode. */
  private boolean maximaIdentification;
  /** True if running in fit maxima mode. */
  private boolean fitMaxima;
  /** True if running in simple fit mode. */
  private boolean simpleFit;

  /** The calculator setting used for the Simple Fit wizard. This is immutable. */
  private static PSFCalculatorSettings calculatorSettings =
      GuiProtosHelper.defaultPSFCalculatorSettings;

  /** The plugin settings. Loaded when a dialog is shown. */
  private Settings settings;
  /**
   * The extra settings for uncommon processing options. This is always loaded as the class can be
   * used via API method calls.
   */
  private ExtraSettings extraSettings = new ExtraSettings();

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    double fractionOfThreads;
    String inputOption;
    boolean showTable;
    boolean showImage;
    boolean fitAcrossAllFrames;

    Settings() {
      // Allow 1 thread free.
      fractionOfThreads = 0.99;
      inputOption = "";
      showTable = true;
      showImage = true;
    }

    Settings(Settings source) {
      fractionOfThreads = source.fractionOfThreads;
      inputOption = source.inputOption;
      showTable = source.showTable;
      showImage = source.showImage;
      fitAcrossAllFrames = source.fitAcrossAllFrames;
    }

    Settings copy() {
      return new Settings(this);
    }

    /**
     * Load a copy of the settings.
     *
     * @return the settings
     */
    static Settings load() {
      return lastSettings.get().copy();
    }

    /**
     * Save the settings.
     */
    void save() {
      lastSettings.set(this);
    }
  }

  /**
   * Contains the extra settings that are the re-usable state of the plugin.
   *
   * <p>These are not commonly used settings. They can be set in the plugin dialog by holding down
   * the shift key.
   */
  private static class ExtraSettings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<ExtraSettings> lastSettings =
        new AtomicReference<>(new ExtraSettings());

    /**
     * This option is only required for the dialog when the input image has a crop. Otherwise the
     * class level {@link PeakFit#ignoreBoundsForNoise} will be set via the API method
     * {@link PeakFit#initialiseImage(ImageSource, Rectangle, boolean)}.
     */
    boolean optionIgnoreBoundsForNoise;
    int integrateFrames;
    boolean interlacedData;
    int dataStart;
    int dataBlock;
    int dataSkip;
    boolean showProcessedFrames;

    ExtraSettings() {
      // Set defaults
      optionIgnoreBoundsForNoise = true;
      integrateFrames = 1;
      dataStart = 1;
      dataBlock = 1;
    }

    ExtraSettings(ExtraSettings source) {
      optionIgnoreBoundsForNoise = source.optionIgnoreBoundsForNoise;
      integrateFrames = source.integrateFrames;
      interlacedData = source.interlacedData;
      dataStart = source.dataStart;
      dataBlock = source.dataBlock;
      dataSkip = source.dataSkip;
      showProcessedFrames = source.showProcessedFrames;
    }

    ExtraSettings copy() {
      return new ExtraSettings(this);
    }

    /**
     * Load a copy of the settings.
     *
     * @return the settings
     */
    static ExtraSettings load() {
      return lastSettings.get().copy();
    }

    /**
     * Save the settings.
     */
    void save() {
      lastSettings.set(this);
    }
  }

  /**
   * Class to control listening to the main dialog and updating fields.
   *
   * <p>Supports changing all field values using the selected template.
   */
  private class ItemDialogListener implements ItemListener {

    // All the fields that will be updated when reloading the configuration file
    private Choice textCameraType;
    private TextField textNmPerPixel;
    private TextField textExposure;
    private Choice textPsf;
    private Choice textDataFilterType;
    private Choice textDataFilterMethod;
    private TextField textSmooth;
    private TextField textSearch;
    private TextField textBorder;
    private TextField textFitting;
    private Choice textFitSolver;
    private Checkbox textFitBackground;
    private TextField textFailuresLimit;
    private TextField textPassRate;
    private Checkbox textIncludeNeighbours;
    private TextField textNeighbourHeightThreshold;
    private TextField textResidualsThreshold;
    private TextField textDuplicateDistance;
    // Cannot be accessed through the GenericDialog API
    private final Scrollbar sliderCoordinateShiftFactor;
    private TextField textCoordinateShiftFactor;
    private TextField textSignalStrength;
    private TextField textMinPhotons;
    private TextField textPrecisionThreshold;
    private Checkbox textSmartFilter;
    private Checkbox textDisableSimpleFilter;
    private TextField textNoise;
    private Choice textNoiseMethod;
    private TextField textMinWidthFactor;
    private TextField textWidthFactor;
    private Checkbox textLogProgress;
    private Checkbox textShowDeviations;
    private Checkbox textResultsTable;
    private Choice textResultsImage;
    private TextField textResultsDirectory;
    private Choice textFileFormat;
    private Checkbox textResultsInMemory;

    ItemDialogListener(Scrollbar sliderCoordinateShiftFactor) {
      // Special case of a Scrollbar that cannot be accessed easily through the GenericDialog API.
      this.sliderCoordinateShiftFactor = sliderCoordinateShiftFactor;
    }

    void attach(GenericDialog gd, boolean isCrop) {

      final Vector<TextField> texts = gd.getStringFields();
      final Vector<TextField> numerics = gd.getNumericFields();
      final Vector<Checkbox> checkboxes = gd.getCheckboxes();
      final Vector<Choice> choices = gd.getChoices();

      final Iterator<TextField> te = texts.iterator();
      final Iterator<TextField> nu = numerics.iterator();
      final Iterator<Checkbox> cb = checkboxes.iterator();
      final Iterator<Choice> ch = choices.iterator();

      final Choice textTemplate = ch.next();
      textTemplate.addItemListener(this);

      textCameraType = ch.next();
      textNmPerPixel = nu.next();
      textExposure = nu.next();
      if (isCrop) {
        cb.next();
      }
      textPsf = ch.next();
      textDataFilterType = ch.next();
      textDataFilterMethod = ch.next();
      textSmooth = nu.next();
      textSearch = nu.next();
      textBorder = nu.next();
      textFitting = nu.next();
      if (extraOptions && !fitMaxima) {
        cb.next(); // Skip over the interlaced data option
        nu.next(); // Skip over the integrate frames option
      }
      if (!maximaIdentification) {
        textFitSolver = ch.next();
        if (extraOptions) {
          textFitBackground = cb.next();
        }
        textFailuresLimit = nu.next();
        textPassRate = nu.next();
        textIncludeNeighbours = cb.next();
        textNeighbourHeightThreshold = nu.next();
        textResidualsThreshold = nu.next();
        textDuplicateDistance = nu.next();
        textSmartFilter = cb.next();
        textDisableSimpleFilter = cb.next();
        textCoordinateShiftFactor = nu.next();
        textSignalStrength = nu.next();
        textMinPhotons = nu.next();
        if (extraOptions) {
          textNoise = nu.next();
          textNoiseMethod = ch.next();
        }
        textMinWidthFactor = nu.next();
        textWidthFactor = nu.next();
        textPrecisionThreshold = nu.next();

        updateFilterInput();
        textSmartFilter.addItemListener(this);
        textDisableSimpleFilter.addItemListener(this);
      }
      textLogProgress = cb.next();
      if (!maximaIdentification) {
        textShowDeviations = cb.next();
      }
      textResultsTable = cb.next();
      textResultsImage = ch.next();
      if (extraOptions) {
        cb.next(); // Skip over show processed frames option
      }
      textResultsDirectory = te.next();

      textFileFormat = ch.next();
      textResultsInMemory = cb.next();
    }

    @Override
    public void itemStateChanged(ItemEvent event) {
      if (event.getSource() instanceof Choice) {
        // Update the settings from the template
        final Choice choice = (Choice) event.getSource();
        final String templateName = choice.getSelectedItem();

        // Get the configuration template
        final TemplateSettings template = ConfigurationTemplate.getTemplate(templateName);

        if (template != null) {
          IJ.log("Applying template: " + templateName);

          for (final String note : template.getNotesList()) {
            IJ.log(note);
          }

          final boolean custom = ConfigurationTemplate.isCustomTemplate(templateName);
          if (template.hasCalibration()) {
            refreshSettings(template.getCalibration());
          }
          if (template.hasPsf()) {
            refreshSettings(template.getPsf(), custom);
          }
          if (template.hasFitEngineSettings()) {
            refreshSettings(template.getFitEngineSettings(), custom);
          }
          if (template.hasResultsSettings()) {
            refreshSettings(template.getResultsSettings());
          }
        }
      } else if (event.getSource() instanceof Checkbox) {
        if (event.getSource() == textSmartFilter) {
          // Prevent both filters being enabled
          textDisableSimpleFilter.setState(textSmartFilter.getState());
          updateFilterInput();
        } else if (event.getSource() == textDisableSimpleFilter) {
          updateFilterInput();
        }
      }
    }

    private void updateFilterInput() {
      if (textDisableSimpleFilter.getState()) {
        sliderCoordinateShiftFactor.setEnabled(false);
        disableEditing(textCoordinateShiftFactor);
        disableEditing(textSignalStrength);
        disableEditing(textMinPhotons);
        // These are used to set bounds
        // disableEditing(textMinWidthFactor)
        // disableEditing(textWidthFactor)
        disableEditing(textPrecisionThreshold);
      } else {
        sliderCoordinateShiftFactor.setEnabled(true);
        enableEditing(textCoordinateShiftFactor);
        enableEditing(textSignalStrength);
        enableEditing(textMinPhotons);
        // enableEditing(textMinWidthFactor)
        // enableEditing(textWidthFactor)
        enableEditing(textPrecisionThreshold);
      }
    }

    private void disableEditing(TextField textField) {
      textField.setEditable(false);
      textField.setBackground(SystemColor.control);
    }

    private void enableEditing(TextField textField) {
      textField.setEditable(true);
      textField.setBackground(Color.white);
    }

    private void refreshSettings(Calibration cal) {
      if (cal == null) {
        return;
      }

      // Do not use set() as we support merging a partial calibration
      fitConfig.mergeCalibration(cal);
      final CalibrationReader calibration = fitConfig.getCalibrationReader();

      textCameraType.select(CalibrationProtosHelper.getName(calibration.getCameraType()));
      if (calibration.hasNmPerPixel()) {
        textNmPerPixel.setText("" + calibration.getNmPerPixel());
      }
      if (calibration.hasExposureTime()) {
        textExposure.setText("" + calibration.getExposureTime());
      }
    }

    private void refreshSettings(PSF psf, boolean isCustomTemplate) {
      if (!isCustomTemplate || psf == null) {
        return;
      }

      // Do not use set() as we support merging a partial PSF
      fitConfig.mergePsf(psf);

      textPsf.select(PsfProtosHelper.getName(fitConfig.getPsfType()));
    }

    /**
     * Refresh settings.
     *
     * <p>If this is a custom template then use all the settings. If a default template then leave
     * some existing spot settings untouched as the user may have updated them (e.g. PSF width).
     *
     * @param fitEngineSettings the config
     * @param isCustomTemplate True if a custom template.
     */
    private void refreshSettings(FitEngineSettings fitEngineSettings, boolean isCustomTemplate) {
      // Set the configuration
      // This will clear everything and merge the configuration
      PeakFit.this.config.setFitEngineSettings(fitEngineSettings);
      fitConfig = PeakFit.this.config.getFitConfiguration();

      textDataFilterType
          .select(SettingsManager.getDataFilterTypeNames()[config.getDataFilterType().ordinal()]);
      textDataFilterMethod.select(
          SettingsManager.getDataFilterMethodNames()[config.getDataFilterMethod(0).ordinal()]);
      textSmooth.setText("" + config.getDataFilterParameterValue(0));
      textSearch.setText("" + config.getSearch());
      textBorder.setText("" + config.getBorder());
      textFitting.setText("" + config.getFitting());
      if (!maximaIdentification) {
        textFitSolver
            .select(SettingsManager.getFitSolverNames()[fitConfig.getFitSolver().ordinal()]);
        if (extraOptions) {
          textFitBackground.setState(fitConfig.isBackgroundFitting());
        }
        textFailuresLimit.setText("" + config.getFailuresLimit());
        textPassRate.setText("" + config.getPassRate());
        textIncludeNeighbours.setState(config.isIncludeNeighbours());
        textNeighbourHeightThreshold.setText("" + config.getNeighbourHeightThreshold());
        textResidualsThreshold.setText("" + config.getResidualsThreshold());
        textDuplicateDistance.setText("" + config.getDuplicateDistance());

        // Filtering
        textSmartFilter.setState(fitConfig.isSmartFilter());
        textDisableSimpleFilter.setState(fitConfig.isDisableSimpleFilter());
        if (!fitConfig.isDisableSimpleFilter()) {
          textCoordinateShiftFactor.setText("" + fitConfig.getCoordinateShiftFactor());
          textSignalStrength.setText("" + fitConfig.getSignalStrength());
          textWidthFactor.setText("" + fitConfig.getMaxWidthFactor());
          textPrecisionThreshold.setText("" + fitConfig.getPrecisionThreshold());
        }
        // These are used for settings the bounds so they are included
        textMinPhotons.setText("" + fitConfig.getMinPhotons());
        textMinWidthFactor.setText("" + fitConfig.getMinWidthFactor());
        updateFilterInput();

        if (extraOptions) {
          textNoise.setText("" + fitConfig.getNoise());
          textNoiseMethod.select(config.getNoiseMethod().ordinal());
        }
      }
    }

    private void refreshSettings(ResultsSettings resultsSettings) {
      PeakFit.this.resultsSettings = resultsSettings.toBuilder();
      textLogProgress.setState(resultsSettings.getLogProgress());
      if (!maximaIdentification) {
        textShowDeviations.setState(resultsSettings.getShowDeviations());
      }
      textResultsTable.setState(resultsSettings.getResultsTableSettings().getShowTable());
      textResultsImage.select(resultsSettings.getResultsImageSettings().getImageTypeValue());
      textResultsDirectory
          .setText("" + resultsSettings.getResultsFileSettings().getResultsDirectory());
      textFileFormat.select(resultsSettings.getResultsFileSettings().getFileFormatValue());
      textResultsInMemory.setState(resultsSettings.getResultsInMemorySettings().getInMemory());
    }
  }

  /**
   * Lazy load the {@link PSFType} values and names.
   */
  private static class PsfTypeLoader {
    private static final PSFType[] psfTypeValues;
    private static final String[] psfTypeNames;

    static {
      //@formatter:off
      final EnumSet<PSFType> d = EnumSet.of(
          PSFType.ONE_AXIS_GAUSSIAN_2D,
          PSFType.TWO_AXIS_GAUSSIAN_2D,
          PSFType.TWO_AXIS_AND_THETA_GAUSSIAN_2D,
          PSFType.ASTIGMATIC_GAUSSIAN_2D);
      //@formatter:on
      psfTypeValues = d.toArray(new PSFType[0]);
      psfTypeNames = new String[psfTypeValues.length];
      for (int i = 0; i < psfTypeValues.length; i++) {
        psfTypeNames[i] = PsfProtosHelper.getName(psfTypeValues[i]);
      }
    }
  }

  /**
   * Instantiates a new peak fit.
   */
  public PeakFit() {
    init(null, null);
  }

  /**
   * Instantiates a new peak fit.
   *
   * @param config the config
   */
  public PeakFit(FitEngineConfiguration config) {
    init(config, null);
  }

  /**
   * Instantiates a new peak fit.
   *
   * @param config the config
   * @param resultsSettings the results settings
   */
  public PeakFit(FitEngineConfiguration config, ResultsSettings resultsSettings) {
    init(config, resultsSettings);
  }

  private void init(FitEngineConfiguration config, ResultsSettings resultsSettings) {
    this.config = (config != null) ? config : new FitEngineConfiguration();
    fitConfig = this.config.getFitConfiguration();
    this.resultsSettings = (resultsSettings != null) ? ResultsSettings.newBuilder(resultsSettings)
        : ResultsSettings.newBuilder();
  }

  @Override
  public int setup(String arg, ImagePlus imp) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    pluginFlags = FLAGS;
    extraOptions = ImageJUtils.isExtraOptions();

    maximaIdentification = StringUtils.contains(arg, "spot");
    fitMaxima = StringUtils.contains(arg, "maxima");
    simpleFit = StringUtils.contains(arg, "simple");
    final boolean runSeries = StringUtils.contains(arg, "series");

    ImageSource imageSource = null;
    if (fitMaxima) {
      // The maxima will have been identified already.
      // The image source will be found from the peak results.
      if (!showMaximaDialog()) {
        return DONE;
      }

      final MemoryPeakResults localResults =
          ResultsManager.loadInputResults(settings.inputOption, false, DistanceUnit.PIXEL);
      if (localResults == null || localResults.size() == 0) {
        IJ.error(TITLE, "No results could be loaded");
        return DONE;
      }

      if (settings.fitAcrossAllFrames) {
        // Allow the user to select a different image. The source will be set as per the
        // main fit routine from the image (imp).
        singleFrame = 0;
      } else {
        // Check for single frame
        singleFrame = getSingleFrame(localResults);
        // Forces the maxima to be used with their original source.
        imp = null;
        imageSource = localResults.getSource();
        pluginFlags |= NO_IMAGE_REQUIRED;
      }
    } else if (runSeries) {
      imp = null;
      // Select input folder
      final String inputDirectory = IJ.getDirectory("Select image series ...");
      if (inputDirectory == null) {
        return DONE;
      }

      // Load input series ...
      SeriesOpener series;
      if (extraOptions) {
        final String helpKey = maximaIdentification ? "spot-finder-series" : "peak-fit-series";
        series = SeriesOpener.create(inputDirectory, true, HelpUrls.getUrl(helpKey));
      } else {
        series = new SeriesOpener(inputDirectory);
      }
      if (series.getNumberOfImages() == 0) {
        IJ.error(TITLE, "No images in the selected directory:\n" + inputDirectory);
        return DONE;
      }

      final SeriesImageSource seriesImageSource =
          new SeriesImageSource(getName(series.getImageList()), series);
      // TrackProgress logging is very verbose if the series has many images
      // Status is used only when reading TIFF info.
      // seriesImageSource.setTrackProgress(SimpleImageJTrackProgress.getInstance());
      seriesImageSource.setTrackProgress(new TrackProgressAdaptor() {
        @Override
        public void status(String format, Object... args) {
          ImageJUtils.showStatus(() -> String.format(format, args));
        }
      });
      imageSource = seriesImageSource;

      pluginFlags |= NO_IMAGE_REQUIRED;
    }

    // If the image source has not been set then use the input image
    if (imageSource == null) {
      if (imp == null) {
        IJ.noImage();
        return DONE;
      }

      // Check it is not a previous result
      if (imp.getTitle().endsWith(ImageJImagePeakResults.IMAGE_SUFFIX)) {
        IJImageSource ijImageSource = null;

        // Check the image to see if it has an image source XML structure in the info property
        final Object o = imp.getProperty("Info");
        final Pattern pattern =
            Pattern.compile("Source: (<.*IJImageSource>.*<.*IJImageSource>)", Pattern.DOTALL);
        final Matcher match = pattern.matcher((o == null) ? "" : o.toString());
        if (match.find()) {
          final ImageSource tmpSource = ImageSource.fromXml(match.group(1));
          if (tmpSource instanceof IJImageSource) {
            ijImageSource = (IJImageSource) tmpSource;
            if (!ijImageSource.open()) {
              ijImageSource = null;
            } else {
              imp = WindowManager.getImage(ijImageSource.getName());
            }
          }
        }

        if (ijImageSource == null) {
          // Look for a parent using the title
          final String parentTitle = imp.getTitle().substring(0,
              imp.getTitle().length() - ImageJImagePeakResults.IMAGE_SUFFIX.length() - 1);
          final ImagePlus parentImp = WindowManager.getImage(parentTitle);
          if (parentImp != null) {
            ijImageSource = new IJImageSource(parentImp);
            imp = parentImp;
          }
        }
        String message = "The selected image may be a previous fit result";
        if (ijImageSource != null) {
          if (!TextUtils.isNullOrEmpty(ijImageSource.getName())) {
            message += " of: \n \n" + ijImageSource.getName();
          }
          message += " \n \nFit the parent?";
        } else {
          message += " \n \nDo you want to continue?";
        }

        final YesNoCancelDialog d = new YesNoCancelDialog(null, TITLE, message);
        if (ijImageSource == null) {
          if (!d.yesPressed()) {
            return DONE;
          }
        } else {
          if (d.yesPressed()) {
            imageSource = ijImageSource;
          }
          if (d.cancelPressed()) {
            return DONE;
          }
        }
      }

      if (imageSource == null) {
        try {
          imageSource = new IJImageSource(imp);
        } catch (final IllegalArgumentException ex) {
          // This can happen if the image has an origin not in integer pixels
          // e.g. the plugin is run on a plot
          IJ.error(TITLE, "Error using image: " + imp.getTitle() + "\n \n" + ex.getMessage());
          return DONE;
        }
      }
    }

    time = -1;

    if (!initialiseImage(imageSource, getBounds(imp), false)) {
      IJ.error(TITLE, "Failed to initialise the source image: " + imageSource.getName());
      return DONE;
    }

    final int flags = showDialog(imp);
    if ((flags & DONE) == 0) {
      initialiseFitting();
    }
    return flags;
  }

  /**
   * Gets the single frame containing all the results (if they are all in a single frame), else 0.
   *
   * @param results the results (must not be empty)
   * @return the single frame (or zero)
   */
  private static int getSingleFrame(MemoryPeakResults results) {
    final FrameCounter counter = new FrameCounter(results.getFirstFrame());
    // The counter will return true (stop execution) if a new frame
    results.forEach((PeakResultProcedureX) peakResult -> counter.advance(peakResult.getFrame()));
    if (counter.currentFrame() != counter.previousFrame()) {
      return 0;
    }
    return counter.currentFrame();
  }

  /**
   * Gets the name.
   *
   * @param imageList the image list
   * @return the name
   */
  static String getName(String[] imageList) {
    String name = imageList[0];
    // Remove directory
    final int index = name.lastIndexOf(File.separatorChar);
    if (index > -1) {
      name = name.substring(index + 1);
    }
    return "Series " + name;
  }

  private static Rectangle getBounds(ImagePlus imp) {
    if (imp == null) {
      return null;
    }
    final Roi roi = imp.getRoi();
    if (roi != null && roi.isArea()) {
      return roi.getBounds();
    }
    return null;
  }

  /**
   * Initialise a new image for fitting and prepare the output results.
   *
   * <p>Calls {@link #initialise(ImageSource, Rectangle, boolean)} then
   * {@link #initialiseFitting()}.
   *
   * @param imageSource The image source
   * @param bounds The region to process from the image
   * @param ignoreBoundsForNoise Set to true if the bounds should be ignored when computing the
   *        noise estimate for each frame
   * @return True if the image was valid and the initialisation was successful
   */
  public boolean initialise(ImageSource imageSource, Rectangle bounds,
      boolean ignoreBoundsForNoise) {
    if (!initialiseImage(imageSource, bounds, ignoreBoundsForNoise)) {
      return false;
    }
    return initialiseFitting();
  }

  /**
   * Initialise a new image.
   *
   * <p>Does not set-up for fitting. This can be done using a subsequent call to
   * {@link #initialiseFitting()}.
   *
   * <p>This mechanism allows additional result outputs to be added after initialisation using
   * {@link #addPeakResults(PeakResults)}.
   *
   * @param imageSource The image source
   * @param bounds The region to process from the image
   * @param ignoreBoundsForNoise Set to true if the bounds should be ignored when computing the
   *        noise estimate for each frame
   * @return True if the image was valid and the initialisation was successful
   */
  public boolean initialiseImage(ImageSource imageSource, Rectangle bounds,
      boolean ignoreBoundsForNoise) {
    // Initialise for image processing
    if (!setSource(imageSource)) {
      return false;
    }

    if (bounds == null) {
      // No region so no need to ignore the bounds.
      this.bounds = new Rectangle(0, 0, source.getWidth(), source.getHeight());
      this.ignoreBoundsForNoise = false;
    } else {
      // The bounds must fit in the image
      try {
        imageSource.checkBounds(bounds);
      } catch (final IllegalArgumentException ex) {
        return false;
      }
      this.bounds = bounds;
      this.ignoreBoundsForNoise = ignoreBoundsForNoise;
    }

    results = new PeakResultsList();

    time = 0;

    return true;
  }

  private boolean setSource(ImageSource imageSource) {
    // Reset
    this.source = null;
    if (imageSource == null) {
      return false;
    }

    // Open the image to ensure it is accessible and the width/height are known
    if (!imageSource.open()) {
      return false;
    }

    this.source = imageSource;
    return true;
  }

  /**
   * Sets the results suffix.
   *
   * @param resultsSuffix the new results suffix
   */
  public void setResultsSuffix(String resultsSuffix) {
    this.resultsSuffix = resultsSuffix;
  }

  /**
   * Set-up the fitting using all the configured properties. Prepare the output results.
   *
   * @return true, if successful
   */
  public boolean initialiseFitting() {
    if (source == null) {
      return false;
    }

    // Do this to ensure the serialised configuration is correct
    updateFitConfiguration(config);

    results.setSource(source);
    String name = source.getName();
    if (!TextUtils.isNullOrEmpty(resultsSuffix)) {
      name += " " + resultsSuffix;
    }
    if (maximaIdentification) {
      name += " (Maxima)";
    } else if (fitMaxima) {
      name += " (" + getSolverName() + " Fit Maxima)";
    } else {
      name += " (" + getSolverName() + ")";
    }
    results.setName(name);
    results.setBounds(bounds);

    // Calibration cal = calibration.clone();
    // Account for the frame integration
    // TODO - Should we change this so that if integrate frames is used then the data
    // are converted to ExtendedPeakResult with a start and end frame
    // cal.exposureTime *= integrateFrames;
    // if (interlacedData)
    // {
    // cal.exposureTime *= ((double)dataBlock / (dataBlock + dataSkip));
    // }

    final CalibrationWriter cal = fitConfig.getCalibrationWriter();
    cal.setTimeUnit(TimeUnit.FRAME);
    results.setCalibration(cal.getCalibration());
    results.setPsf(fitConfig.getPsf());
    results.setConfiguration(SettingsManager.toJson(config.getFitEngineSettings()));

    // This is added first as it cannot be closed. If the table is closed then the
    // number of results at the end is reported incorrectly.
    addMemoryResults(results, false);

    addTableResults(results);
    ResultsManager.addImageResults(results, resultsSettings.getResultsImageSettings(), bounds,
        (extraOptions) ? ResultsManager.FLAG_EXTRA_OPTIONS : 0);
    addFileResults(results);
    addDefaultResults(results);

    results.begin();

    if (simpleFit && settings.showImage) {
      for (final PeakResults r : results.toArray()) {
        if (r instanceof ImageJImagePeakResults) {
          final ImagePlus i = ((ImageJImagePeakResults) r).getImagePlus();
          ImageJUtils.log("Super-resolution image title = " + i.getTitle());
          WindowManager.toFront(i.getWindow());
        }
      }
    }

    return true;
  }

  private String getSolverName() {
    return getSolverName(config.getFitConfiguration());
  }

  /**
   * Gets the solver name.
   *
   * @param fitConfig the fit config
   * @return the solver name
   */
  public static String getSolverName(FitConfiguration fitConfig) {
    final FitSolver solver = fitConfig.getFitSolver();
    String name = FitProtosHelper.getName(solver);
    if (solver == FitSolver.MLE) {
      name += " " + FitProtosHelper.getName(fitConfig.getSearchMethod());
    }
    return name;
  }

  /**
   * Show results.
   */
  protected void showResults() {
    IJ.showProgress(1.0);
    if (time >= 0) {
      if (silent) {
        results.end();
        return;
      }

      // Check if we are sorting
      IJ.showStatus("Finalising results ...");
      for (final PeakResults r : results.toArray()) {
        if (r instanceof MemoryPeakResults) {
          if (((MemoryPeakResults) r).isSortAfterEnd()) {
            IJ.showStatus("Sorting " + r.size() + " results ...");
          }
          break;
        }
      }

      results.end();

      final String textTime = TextUtils.nanosToString(time);
      final String textRunTime = TextUtils.nanosToString(runTime);

      final int size = getSize();
      final String message = String.format("%s. Total fitting time = %s. Run time = %s",
          TextUtils.pleural(size, "localisation"), textTime, textRunTime);
      if (resultsSettings.getLogProgress()) {
        IJ.log(LOG_SPACER);
      }
      IJ.log(message);
      IJ.showStatus(message);
    } else {
      IJ.showStatus("");
    }
  }

  private boolean showMaximaDialog() {
    final int size = MemoryPeakResults.countMemorySize();
    if (size == 0) {
      IJ.error(TITLE, "There are no fitting results in memory");
      return false;
    }

    settings = Settings.load();
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("fit-maxima"));
    gd.addMessage("Select identified maxima for fitting");

    ResultsManager.addInput(gd, settings.inputOption, InputSource.MEMORY);
    gd.addCheckbox("Fit_across_all_frames", settings.fitAcrossAllFrames);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.inputOption = ResultsManager.getInputSource(gd);
    settings.fitAcrossAllFrames = gd.getNextBoolean();
    settings.save();

    return true;
  }

  /**
   * Gets the PSF type values.
   *
   * @return the PSF type values
   */
  public static PSFType[] getPsfTypeValues() {
    return PsfTypeLoader.psfTypeValues;
  }

  /**
   * Gets the PSF type names.
   *
   * @return the PSF type names
   */
  public static String[] getPsfTypeNames() {
    return PsfTypeLoader.psfTypeNames;
  }

  private int showDialog(ImagePlus imp) {
    // Executing as an ImageJ plugin.

    // Load the settings
    resultsSettings = SettingsManager.readResultsSettings(0).toBuilder();
    // Settings are within the FitEngineSettings
    config = SettingsManager.readFitEngineConfiguration(0);
    fitConfig = config.getFitConfiguration();
    settings = Settings.load();

    if (simpleFit) {
      return showSimpleDialog();
    }

    // Note: The bounds are not set when running in the fit maxima option (since all candidates
    // have been identified already the crop is not required).
    final boolean isCrop = (bounds != null && imp != null
        && (bounds.width < imp.getWidth() || bounds.height < imp.getHeight()));

    // Some options are not always needed
    if (extraOptions || isCrop) {
      extraSettings = ExtraSettings.load();
      ignoreBoundsForNoise = extraSettings.optionIgnoreBoundsForNoise;
    }

    if (!extraOptions) {
      resultsSettings.getResultsImageSettingsBuilder().setRollingWindowSize(0);
      fitConfig.setBackgroundFitting(true);
      fitConfig.setNoise(0);
      config.setNoiseMethod(NoiseEstimatorMethod.QUICK_RESIDUALS_LEAST_MEAN_OF_SQUARES);
    }

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    String helpKey;
    if (maximaIdentification) {
      helpKey = "spot-finder";
      gd.addMessage("Identify candidate maxima");
    } else {
      helpKey = "peak-fit";
      gd.addMessage("Fit 2D Gaussian to identified maxima");
    }
    // Note: Currently is is not useful to append "-series" when running for a series image
    // source since the params in this dialog do not concern the image input.
    gd.addHelp(HelpUrls.getUrl(helpKey));

    final String[] templates = ConfigurationTemplate.getTemplateNames(true);
    gd.addChoice("Template", templates, templates[0]);

    final CalibrationReader calibration = fitConfig.getCalibrationReader();
    addCameraOptions(gd, 0, fitConfig);
    gd.addNumericField("Calibration", calibration.getNmPerPixel(), 2, 6, "nm/px");
    gd.addNumericField("Exposure_time", calibration.getExposureTime(), 2, 6, "ms");

    if (isCrop) {
      gd.addCheckbox("Ignore_bounds_for_noise", ignoreBoundsForNoise);
    }

    final FitConfigurationProvider fitConfigurationProvider = () -> fitConfig;
    final FitEngineConfigurationProvider fitEngineConfigurationProvider = () -> config;

    addPsfOptions(gd, fitConfigurationProvider);
    addDataFilterOptions(gd, fitEngineConfigurationProvider);
    addSearchOptions(gd, fitEngineConfigurationProvider);
    addBorderOptions(gd, fitEngineConfigurationProvider);
    addFittingOptions(gd, fitEngineConfigurationProvider);
    if (extraOptions && !fitMaxima) {
      gd.addCheckbox("Interlaced_data", extraSettings.interlacedData);
      gd.addSlider("Integrate_frames", 1, 5, extraSettings.integrateFrames);
    }

    // Special case top get the slider since the GenericDialog does not provide access to this.
    Scrollbar sliderCoordinateShiftFactor = null;
    final boolean isShowGenericDialog = ImageJUtils.isShowGenericDialog();

    if (!maximaIdentification) {
      gd.addMessage("--- Gaussian fitting ---");
      gd.addChoice("Fit_solver", SettingsManager.getFitSolverNames(),
          fitConfig.getFitSolver().ordinal());
      if (extraOptions) {
        gd.addCheckbox("Fit_background", fitConfig.isBackgroundFitting());
      }

      // Parameters specific to each Fit solver are collected in a second dialog

      gd.addNumericField("Fail_limit", config.getFailuresLimit(), 0);
      gd.addNumericField("Pass_rate", config.getPassRate(), 2);
      gd.addCheckbox("Include_neighbours", config.isIncludeNeighbours());
      gd.addSlider("Neighbour_height", 0.01, 1, config.getNeighbourHeightThreshold());
      gd.addSlider("Residuals_threshold", 0.01, 1, config.getResidualsThreshold());

      addDuplicateDistanceOptions(gd, fitEngineConfigurationProvider);

      gd.addMessage(
          "--- Peak filtering ---\nDiscard fits that shift; are too low; or expand/contract");

      gd.addCheckbox("Smart_filter", fitConfig.isSmartFilter());
      gd.addCheckbox("Disable_simple_filter", fitConfig.isDisableSimpleFilter());
      gd.addSlider("Shift_factor", 0.0, 2.5, fitConfig.getCoordinateShiftFactor());
      if (isShowGenericDialog) {
        sliderCoordinateShiftFactor = gd.getLastScrollbar();
      }
      gd.addNumericField("Signal_strength", fitConfig.getSignalStrength(), 2);
      gd.addNumericField("Min_photons", fitConfig.getMinPhotons(), 0);
      if (extraOptions) {
        gd.addNumericField("Noise", fitConfig.getNoise(), 2);
        gd.addChoice("Noise_method", SettingsManager.getNoiseEstimatorMethodNames(),
            config.getNoiseMethod().ordinal());
      }
      gd.addSlider("Min_width_factor", 0, 0.99, fitConfig.getMinWidthFactor());
      gd.addSlider("Width_factor", 1, 4.5, fitConfig.getMaxWidthFactor());
      addPrecisionOptions(gd, fitConfigurationProvider);
      // Q. Add dynamically displayed options for z-filtering here?
    }

    gd.addMessage("--- Results ---");
    gd.addCheckbox("Log_progress", resultsSettings.getLogProgress());
    if (!maximaIdentification) {
      gd.addCheckbox("Show_deviations", resultsSettings.getShowDeviations());
    }
    ResultsManager.addTableResultsOptions(gd, resultsSettings);
    ResultsManager.addImageResultsOptions(gd, resultsSettings,
        (extraOptions) ? ResultsManager.FLAG_EXTRA_OPTIONS : 0);
    if (extraOptions) {
      gd.addCheckbox("Show_processed_frames", extraSettings.showProcessedFrames);
    }
    ResultsManager.addFileResultsOptions(gd, resultsSettings,
        ResultsManager.FLAG_RESULTS_DIRECTORY);
    ResultsManager.addInMemoryResultsOptions(gd, resultsSettings);

    if (extraOptions) {
      gd.addMessage("--- Misc ---");
      gd.addSlider("Fraction_of_threads", 0.1, 1, settings.fractionOfThreads);
    }

    // Add a mouse listener to the config file field
    if (isShowGenericDialog) {
      new ItemDialogListener(sliderCoordinateShiftFactor).attach(gd, isCrop);
    }

    gd.showDialog();

    if (gd.wasCanceled() || !readDialog(gd, isCrop)) {
      return DONE;
    }

    if (imp != null) {
      // Store whether the user selected to process all the images.
      final int flags = IJ.setupDialog(imp, pluginFlags);

      // Check if cancelled
      if ((flags & DONE) != 0) {
        return DONE;
      }

      if ((flags & DOES_STACKS) == 0) {
        // Save the slice number for the overlay
        singleFrame = imp.getCurrentSlice();

        // Account for interlaced data
        if (extraSettings.interlacedData) {
          int start = singleFrame;

          // Calculate the first frame that is not skipped
          while (ignoreFrame(start) && start > extraSettings.dataStart) {
            start--;
          }
          if (start < extraSettings.dataStart) {
            log("The current frame (%d) is before the start of the interlaced data", singleFrame);
            return DONE;
          }
          if (start != singleFrame) {
            log("Updated the current frame (%d) to a valid interlaced data frame (%d)", singleFrame,
                start);
          }
          singleFrame = start;
        }

        // Account for integrated frames
        int endFrame = singleFrame;
        if (extraSettings.integrateFrames > 1) {
          int totalFrames = 1;
          while (totalFrames < extraSettings.integrateFrames) {
            endFrame++;
            if (!ignoreFrame(endFrame)) {
              totalFrames++;
            }
          }
          log("Updated the image end frame (%d) to %d allow %d integrated frames", singleFrame,
              endFrame, extraSettings.integrateFrames);
        }

        // Create a new image source with the correct frames
        setSource(new IJImageSource(imp, singleFrame, endFrame - singleFrame));

        // Store the image so the results can be added as an overlay
        this.imp = imp;
        this.imp.setOverlay(null);
      }
    }

    // Allow interlaced data by wrapping the image source
    if (extraSettings.interlacedData) {
      setSource(new InterlacedImageSource(this.source, extraSettings.dataStart,
          extraSettings.dataBlock, extraSettings.dataSkip));
    }

    // Allow frame aggregation by wrapping the image source
    if (extraSettings.integrateFrames > 1) {
      setSource(new AggregatedImageSource(this.source, extraSettings.integrateFrames));
    }

    // Ask if the user wants to log progress on multiple frame images
    if (resultsSettings.getLogProgress() && source.getFrames() > 1) {
      final ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE);
      egd.addMessage("Warning: Log progress on multiple-frame image will be slow");
      egd.addCheckbox("Log_progress", true);
      egd.showDialog();
      if (egd.wasCanceled()) {
        return DONE;
      }
      if (!egd.getNextBoolean()) {
        resultsSettings.setLogProgress(false);
        SettingsManager.writeSettings(resultsSettings.build());
      }
    }

    // Return the plugin flags (without the DOES_STACKS flag).
    // The call to run(ImageProcessor) will process the image in 'this.imp' so we only want a
    // single call to be made.
    return pluginFlags;
  }

  /**
   * Allow the latest calibration to be provided for update.
   */
  public interface CalibrationProvider {
    /**
     * Gets the calibration.
     *
     * @return the calibration
     */
    Calibration getCalibration();

    /**
     * Save the calibration. This is used to save changes to the calibration back to the provider.
     *
     * @param calibration the calibration
     */
    void saveCalibration(Calibration calibration);
  }

  /**
   * Adds the camera options.
   *
   * @param gd the dialog
   * @param fitConfig the fit config
   */
  public static void addCameraOptions(final ExtendedGenericDialog gd,
      final FitConfiguration fitConfig) {
    addCameraOptions(gd, 0, fitConfig);
  }

  /**
   * Adds the camera options.
   *
   * @param gd the dialog
   * @param options the options
   * @param fitConfig the fit config
   */
  public static void addCameraOptions(final ExtendedGenericDialog gd, final int options,
      final FitConfiguration fitConfig) {
    addCameraOptions(gd, options, new CalibrationProvider() {
      @Override
      public Calibration getCalibration() {
        return fitConfig.getCalibration();
      }

      @Override
      public void saveCalibration(Calibration calibration) {
        fitConfig.setCalibration(calibration);
      }
    });
  }

  /**
   * Adds the camera options.
   *
   * @param gd the dialog
   * @param calibrationWriter the calibration writer
   */
  public static void addCameraOptions(final ExtendedGenericDialog gd,
      final CalibrationWriter calibrationWriter) {
    addCameraOptions(gd, 0, calibrationWriter);
  }

  /**
   * Adds the camera options.
   *
   * @param gd the dialog
   * @param options the options
   * @param calibrationWriter the calibration writer
   */
  public static void addCameraOptions(final ExtendedGenericDialog gd, final int options,
      final CalibrationWriter calibrationWriter) {
    addCameraOptions(gd, options, new PeakFit.CalibrationProvider() {
      @Override
      public Calibration getCalibration() {
        return calibrationWriter.getCalibration();
      }

      @Override
      public void saveCalibration(Calibration calibration) {
        calibrationWriter.setCalibration(calibration);
      }
    });
  }

  /** Flag to indicate that read noise should not be configured. */
  public static final int FLAG_NO_READ_NOISE = 0x00000001;
  /** Flag to indicate that quantum efficiency should be configured. */
  public static final int FLAG_QUANTUM_EFFICIENCY = 0x00000002;
  /** Flag to indicate that gain should not be configured. */
  public static final int FLAG_NO_GAIN = 0x00000004;

  /**
   * Adds the camera options.
   *
   * @param gd the dialog
   * @param options the options
   * @param calibrationProvider the calibration provider
   */
  public static void addCameraOptions(final ExtendedGenericDialog gd, final int options,
      final CalibrationProvider calibrationProvider) {
    final CalibrationReader calibration =
        new CalibrationReader(calibrationProvider.getCalibration());

    gd.addChoice("Camera_type", SettingsManager.getCameraTypeNames(),
        CalibrationProtosHelper.getName(calibration.getCameraType()),
        new OptionListener<Integer>() {
          @Override
          public boolean collectOptions(Integer field) {
            final CalibrationWriter calibration =
                new CalibrationWriter(calibrationProvider.getCalibration());
            final CameraType t = SettingsManager.getCameraTypeValues()[field];
            if (calibration.getCameraType() != t) {
              calibration.setCameraType(t);
              calibrationProvider.saveCalibration(calibration.getCalibration());
            }
            return collectOptions(false);
          }

          @Override
          public boolean collectOptions() {
            return collectOptions(true);
          }

          private boolean collectOptions(boolean silent) {
            final CalibrationWriter calibration =
                new CalibrationWriter(calibrationProvider.getCalibration());
            final ExtendedGenericDialog egd =
                new ExtendedGenericDialog("Camera type options", null);
            if (calibration.isCcdCamera()) {
              egd.addNumericField("Camera_bias", calibration.getBias(), 2, 6, "Count");
              if (BitFlagUtils.anyNotSet(options, FLAG_NO_GAIN)) {
                egd.addNumericField("Gain", calibration.getCountPerPhoton(), 4, 6, "Count/photon");
              }
              if (BitFlagUtils.anyNotSet(options, FLAG_NO_READ_NOISE)) {
                egd.addNumericField("Read_noise", calibration.getReadNoise(), 4, 6, "Count");
              }
              if (BitFlagUtils.areSet(options, FLAG_QUANTUM_EFFICIENCY)) {
                egd.addNumericField("Quantum_efficiency", calibration.getQuantumEfficiency(), 4, 6,
                    "electron/photon");
              }
            } else if (calibration.isScmos()) {
              final String[] models = CameraModelManager.listCameraModels(true);
              egd.addChoice("Camera_model_name", models, calibration.getCameraModelName());
              if (BitFlagUtils.areSet(options, FLAG_QUANTUM_EFFICIENCY)) {
                egd.addNumericField("Quantum_efficiency", calibration.getQuantumEfficiency(), 4, 6,
                    "electron/photon");
              }
            } else {
              IJ.error("Unsupported camera type "
                  + CalibrationProtosHelper.getName(calibration.getCameraType()));
              return false;
            }
            egd.setSilent(silent);
            egd.showDialog(true, gd);
            if (egd.wasCanceled()) {
              return false;
            }
            final Calibration old = calibration.getCalibration();
            if (calibration.isCcdCamera()) {
              calibration.setBias(Math.abs(egd.getNextNumber()));
              if (BitFlagUtils.anyNotSet(options, FLAG_NO_GAIN)) {
                calibration.setCountPerPhoton(Math.abs(egd.getNextNumber()));
              }
              if (BitFlagUtils.anyNotSet(options, FLAG_NO_READ_NOISE)) {
                calibration.setReadNoise(Math.abs(egd.getNextNumber()));
              }
              if (BitFlagUtils.areSet(options, FLAG_QUANTUM_EFFICIENCY)) {
                calibration.setQuantumEfficiency(Math.abs(egd.getNextNumber()));
              }
            } else if (calibration.isScmos()) {
              // Note: Since this does not go through the FitConfiguration object the
              // camera model is not invalidated. However any code using this function
              // should later call configureFitSolver(...) which will set the camera model
              // using the camera model name.
              calibration.setCameraModelName(egd.getNextChoice());
              if (BitFlagUtils.areSet(options, FLAG_QUANTUM_EFFICIENCY)) {
                calibration.setQuantumEfficiency(Math.abs(egd.getNextNumber()));
              }
            }
            final Calibration current = calibration.getCalibration();
            final boolean changed = !old.equals(current);
            if (changed) {
              calibrationProvider.saveCalibration(current);
            }
            return changed;
          }
        });
  }

  /**
   * Allow the latest fitEngineConfiguration to be provided for update.
   */
  @FunctionalInterface
  public interface FitEngineConfigurationProvider {
    /**
     * Gets the fitEngineConfiguration.
     *
     * @return the fitEngineConfiguration
     */
    FitEngineConfiguration getFitEngineConfiguration();
  }

  /**
   * Allow the latest fitConfiguration to be provided for update.
   */
  @FunctionalInterface
  public interface FitConfigurationProvider {
    /**
     * Gets the fitConfiguration.
     *
     * @return the fitConfiguration
     */
    FitConfiguration getFitConfiguration();
  }

  /**
   * Simple implementation of {@link FitEngineConfigurationProvider}.
   */
  public static class SimpleFitEngineConfigurationProvider
      implements FitEngineConfigurationProvider {
    /** The fit engine configuration. */
    private final FitEngineConfiguration config;

    /**
     * Instantiates a new simple fit engine configuration provider.
     *
     * @param config the configuration
     */
    public SimpleFitEngineConfigurationProvider(FitEngineConfiguration config) {
      this.config = config;
    }

    @Override
    public FitEngineConfiguration getFitEngineConfiguration() {
      return config;
    }
  }

  /**
   * Simple implementation of {@link FitConfigurationProvider}.
   */
  public static class SimpleFitConfigurationProvider implements FitConfigurationProvider {
    /** The configuration. */
    private final FitConfiguration fitConfig;

    /**
     * Instantiates a new simple fit configuration provider.
     *
     * @param fitConfig the configuration
     */
    public SimpleFitConfigurationProvider(FitConfiguration fitConfig) {
      this.fitConfig = fitConfig;
    }

    @Override
    public FitConfiguration getFitConfiguration() {
      return fitConfig;
    }
  }

  /**
   * Adds the PSF options.
   *
   * <p>Note that if an astigmatic PSF is selected then the model must be created with
   * {@link #configurePsfModel(FitEngineConfiguration, int)}.
   *
   * @param gd the dialog
   * @param fitConfiguration the fit configuration
   */
  public static void addPsfOptions(final ExtendedGenericDialog gd,
      final FitConfiguration fitConfiguration) {
    addPsfOptions(gd, new SimpleFitConfigurationProvider(fitConfiguration));
  }

  /**
   * Adds the PSF options.
   *
   * <p>Note that if an astigmatic PSF is selected then the model must be created with
   * {@link #configurePsfModel(FitEngineConfiguration, int)}.
   *
   * @param gd the dialog
   * @param fitConfigurationProvider the fit configuration provider
   */
  public static void addPsfOptions(final ExtendedGenericDialog gd,
      final FitConfigurationProvider fitConfigurationProvider) {
    final FitConfiguration fitConfig = fitConfigurationProvider.getFitConfiguration();
    gd.addChoice("PSF", getPsfTypeNames(), PsfProtosHelper.getName(fitConfig.getPsfType()),
        new OptionListener<Integer>() {
          @Override
          public boolean collectOptions(Integer field) {
            final FitConfiguration fitConfig = fitConfigurationProvider.getFitConfiguration();
            fitConfig.setPsfType(PeakFit.getPsfTypeValues()[field]);
            return collectOptions(false);
          }

          @Override
          public boolean collectOptions() {
            return collectOptions(true);
          }

          private boolean collectOptions(boolean silent) {
            final FitConfiguration localFitConfig = fitConfigurationProvider.getFitConfiguration();
            final PSFType psfType = localFitConfig.getPsfType();
            final ExtendedGenericDialog egd = new ExtendedGenericDialog("PSF Options", null);
            PSF oldPsf = null;
            if (psfType == PSFType.ASTIGMATIC_GAUSSIAN_2D) {
              // The PSF is entirely defined in the model
              String[] list = AstigmatismModelManager.listAstigmatismModels(false,
                  localFitConfig.getCalibrationReader().getNmPerPixel(), 0.1);
              // In case the calibration has not been updated
              if (list.length == 0) {
                list = AstigmatismModelManager.listAstigmatismModels(false, true);
              }
              egd.addChoice("Z-model", list, localFitConfig.getPsfModelName());
            } else {
              // Collect the PSF parameters
              oldPsf = localFitConfig.getPsf();
              for (int i = 0; i < oldPsf.getParametersCount(); i++) {
                final PSFParameter p = oldPsf.getParameters(i);
                egd.addNumericField(String.format("PSF_parameter_%d (%s)", i + 1, p.getName()),
                    p.getValue(), 3);
              }
              if (psfType == PSFType.ONE_AXIS_GAUSSIAN_2D) {
                egd.addCheckbox("Fixed", localFitConfig.isFixedPsf());
              }
            }

            egd.setSilent(silent);
            egd.showDialog(true, gd);
            if (egd.wasCanceled()) {
              return false;
            }
            if (psfType == PSFType.ASTIGMATIC_GAUSSIAN_2D) {
              // The PSF is entirely defined in the model
              localFitConfig.setPsfModelName(egd.getNextChoice());
              return true;
            }
            @SuppressWarnings("null")
            final PSF.Builder b = oldPsf.toBuilder();
            final int n = b.getParametersCount();
            for (int i = 0; i < n; i++) {
              b.getParametersBuilder(i).setValue(egd.getNextNumber());
            }
            final PSF newPsf = b.build();
            localFitConfig.setPsf(newPsf);
            boolean changed = !oldPsf.equals(newPsf);
            if (psfType == PSFType.ONE_AXIS_GAUSSIAN_2D) {
              final boolean newFixed = egd.getNextBoolean();
              changed = changed || (newFixed != localFitConfig.isFixedPsf());
              localFitConfig.setFixedPsf(newFixed);
            }
            return changed;
          }
        });
  }

  /**
   * Used for relative parameters.
   */
  abstract static class RelativeParameterProvider {
    /** The min. */
    double min;
    /** The max. */
    double max;
    /** The name. */
    String name;
    /** The fit engine configuration provider. */
    FitEngineConfigurationProvider fitEngineConfigurationProvider;
    /**
     * Set to true to include the value in {@link #getMin()} and {@link #getMax()}.
     */
    boolean includeValue;

    /**
     * Instantiates a new relative parameter provider.
     *
     * @param min the min
     * @param max the max
     * @param name the name
     * @param fitEngineConfigurationProvider the fit engine configuration provider
     */
    public RelativeParameterProvider(double min, double max, String name,
        FitEngineConfigurationProvider fitEngineConfigurationProvider) {
      this(min, max, name, fitEngineConfigurationProvider, false);
    }

    /**
     * Instantiates a new relative parameter provider.
     *
     * @param min the min
     * @param max the max
     * @param name the name
     * @param fitEngineConfigurationProvider the fit engine configuration provider
     * @param includeValue the include value
     */
    public RelativeParameterProvider(double min, double max, String name,
        FitEngineConfigurationProvider fitEngineConfigurationProvider, boolean includeValue) {
      this.min = min;
      this.max = max;
      this.name = name;
      this.fitEngineConfigurationProvider = fitEngineConfigurationProvider;
      this.includeValue = includeValue;
    }

    /**
     * Gets the min.
     *
     * @return the min
     */
    double getMin() {
      return (includeValue) ? Math.min(min, getValue()) : min;
    }

    /**
     * Gets the max.
     *
     * @return the max
     */
    double getMax() {
      return (includeValue) ? Math.max(max, getValue()) : max;
    }

    /**
     * Gets the value.
     *
     * @return the value
     */
    abstract double getValue();

    /**
     * Checks if is absolute.
     *
     * @return true, if is absolute
     */
    abstract boolean isAbsolute();

    /**
     * Sets the absolute.
     *
     * @param absolute the new absolute
     */
    abstract void setAbsolute(boolean absolute);

    /**
     * Gets the dialog name.
     *
     * @return the dialog name
     */
    String getDialogName() {
      return name.replace(' ', '_');
    }
  }

  /**
   * Adds the relative parameter options.
   *
   * @param gd the dialog
   * @param rp the relative parameter
   */
  static void addRelativeParameterOptions(final ExtendedGenericDialog gd,
      final RelativeParameterProvider rp) {
    final String label = rp.getDialogName();
    gd.addSlider(label, rp.getMin(), rp.getMax(), rp.getValue(), new OptionListener<Double>() {
      @Override
      public boolean collectOptions(Double value) {
        // Nothing depends on the input double value so just collect the options
        return collectOptions(false);
      }

      @Override
      public boolean collectOptions() {
        return collectOptions(true);
      }

      private boolean collectOptions(boolean silent) {
        final ExtendedGenericDialog egd = new ExtendedGenericDialog(rp.name + " Options", null);
        final boolean oldValue = rp.isAbsolute();
        egd.addCheckbox(rp.getDialogName() + "_absolute", oldValue);
        egd.setSilent(silent);
        egd.showDialog(true, gd);
        if (egd.wasCanceled()) {
          return false;
        }
        final boolean newValue = egd.getNextBoolean();
        rp.setAbsolute(newValue);
        return oldValue != newValue;
      }
    });

    // Add a label after the button.
    // The button is added in a panel with a GridBagLayout.
    final Panel p = gd.getLastPanel();
    final GridBagConstraints pc = new GridBagConstraints();
    pc.gridy = 0;
    pc.gridx = 3;
    pc.insets = new Insets(5, 5, 0, 0);
    pc.anchor = GridBagConstraints.EAST;
    final Label flagLabel = new Label();
    updateFlag(flagLabel, rp.isAbsolute());
    p.add(flagLabel, pc);

    gd.addOptionCollectedListener(event -> {
      if (label.equals(event.getLabel())) {
        updateFlag(flagLabel, rp.isAbsolute());
      }
    });
  }

  private static void updateFlag(Label label, boolean absolute) {
    // Note: Add spaces so the label has extra width for font width differences
    if (absolute) {
      label.setText("Absolute   ");
    } else {
      label.setText("Relative   ");
    }
  }

  /**
   * Adds the data filter options for the first filter. Adds to the dialog: <ul> <li>a choice of
   * filter type (e.g. single, difference, etc)</li> <li>a choice of primary filter (e.g. mean,
   * Gaussian, etc)</li> <li>a single slider for the primary filter parameter</li> </ul>
   *
   * @param gd the dialog
   * @param fitEngineConfigurationProvider the fit engine configuration provider
   */
  public static void addDataFilterOptions(final ExtendedGenericDialog gd,
      final FitEngineConfigurationProvider fitEngineConfigurationProvider) {
    final int n = 0;
    final DataFilterMethod defaultFilterMethod = DataFilterMethod.GAUSSIAN;
    final double defaultFilterSmoothing = 0.5;
    final FitEngineConfiguration config =
        fitEngineConfigurationProvider.getFitEngineConfiguration();
    gd.addChoice("Spot_filter_type", SettingsManager.getDataFilterTypeNames(),
        config.getDataFilterType().ordinal());
    gd.addChoice("Spot_filter", SettingsManager.getDataFilterMethodNames(),
        config.getDataFilterMethod(n, defaultFilterMethod).ordinal());
    addRelativeParameterOptions(gd,
        new RelativeParameterProvider(0, 2.5, "Smoothing", fitEngineConfigurationProvider, true) {
          @Override
          void setAbsolute(boolean absolute) {
            final FitEngineConfiguration c =
                fitEngineConfigurationProvider.getFitEngineConfiguration();
            final DataFilterMethod m = c.getDataFilterMethod(n, defaultFilterMethod);
            final double smooth = c.getDataFilterParameterValue(n, defaultFilterSmoothing);
            c.setDataFilter(m, smooth, absolute, n);
          }

          @Override
          boolean isAbsolute() {
            return fitEngineConfigurationProvider.getFitEngineConfiguration()
                .getDataFilterParameterAbsolute(n, false);
          }

          @Override
          double getValue() {
            return fitEngineConfigurationProvider.getFitEngineConfiguration()
                .getDataFilterParameterValue(n, defaultFilterSmoothing);
          }
        });
  }

  /**
   * Adds the search options. A single slider for the search parameter is added to the dialog.
   *
   * @param gd the dialog
   * @param fitEngineConfigurationProvider the fit engine configuration provider
   */
  public static void addSearchOptions(final ExtendedGenericDialog gd,
      final FitEngineConfigurationProvider fitEngineConfigurationProvider) {
    addRelativeParameterOptions(gd, new RelativeParameterProvider(0.5, 2.5, "Search Width",
        fitEngineConfigurationProvider, true) {
      @Override
      void setAbsolute(boolean absolute) {
        fitEngineConfigurationProvider.getFitEngineConfiguration().setSearchAbsolute(absolute);
      }

      @Override
      boolean isAbsolute() {
        return fitEngineConfigurationProvider.getFitEngineConfiguration().getSearchAbsolute();
      }

      @Override
      double getValue() {
        return fitEngineConfigurationProvider.getFitEngineConfiguration().getSearch();
      }
    });
  }

  /**
   * Adds the border options. A single slider for the border parameter is added to the dialog.
   *
   * @param gd the dialog
   * @param fitEngineConfigurationProvider the fit engine configuration provider
   */
  public static void addBorderOptions(final ExtendedGenericDialog gd,
      final FitEngineConfigurationProvider fitEngineConfigurationProvider) {
    addRelativeParameterOptions(gd, new RelativeParameterProvider(0.5, 2.5, "Border Width",
        fitEngineConfigurationProvider, true) {
      @Override
      void setAbsolute(boolean absolute) {
        fitEngineConfigurationProvider.getFitEngineConfiguration().setBorderAbsolute(absolute);
      }

      @Override
      boolean isAbsolute() {
        return fitEngineConfigurationProvider.getFitEngineConfiguration().getBorderAbsolute();
      }

      @Override
      double getValue() {
        return fitEngineConfigurationProvider.getFitEngineConfiguration().getBorder();
      }
    });
  }

  /**
   * Adds the fitting options. A single slider for the fitting parameter is added to the dialog.
   *
   * @param gd the dialog
   * @param fitEngineConfigurationProvider the fit engine configuration provider
   */
  public static void addFittingOptions(final ExtendedGenericDialog gd,
      final FitEngineConfigurationProvider fitEngineConfigurationProvider) {
    // For this we allow the slider range to increase as the user may have a large fit width
    addRelativeParameterOptions(gd, new RelativeParameterProvider(2, 4.5, "Fitting Width",
        fitEngineConfigurationProvider, true) {
      @Override
      void setAbsolute(boolean absolute) {
        fitEngineConfigurationProvider.getFitEngineConfiguration().setFittingAbsolute(absolute);
      }

      @Override
      boolean isAbsolute() {
        return fitEngineConfigurationProvider.getFitEngineConfiguration().getFittingAbsolute();
      }

      @Override
      double getValue() {
        return fitEngineConfigurationProvider.getFitEngineConfiguration().getFitting();
      }
    });
  }

  /**
   * Adds the duplicate distance options. A single slider for the duplicate distance parameter is
   * added to the dialog.
   *
   * @param gd the dialog
   * @param fitEngineConfigurationProvider the fit engine configuration provider
   */
  public static void addDuplicateDistanceOptions(final ExtendedGenericDialog gd,
      final FitEngineConfigurationProvider fitEngineConfigurationProvider) {
    addRelativeParameterOptions(gd, new RelativeParameterProvider(0, 1.5, "Duplicate Distance",
        fitEngineConfigurationProvider) {
      @Override
      void setAbsolute(boolean absolute) {
        fitEngineConfigurationProvider.getFitEngineConfiguration()
            .setDuplicateDistanceAbsolute(absolute);
      }

      @Override
      boolean isAbsolute() {
        return fitEngineConfigurationProvider.getFitEngineConfiguration()
            .getDuplicateDistanceAbsolute();
      }

      @Override
      double getValue() {
        return fitEngineConfigurationProvider.getFitEngineConfiguration().getDuplicateDistance();
      }
    });
  }

  /**
   * Adds the precision options. A single numeric field for the precision is added. A pop-up is
   * added to allow the precision method to be configured.
   *
   * @param gd the dialog
   * @param fitConfigurationProvider the fit configuration provider
   */
  public static void addPrecisionOptions(final ExtendedGenericDialog gd,
      final FitConfigurationProvider fitConfigurationProvider) {
    gd.addNumericField("Precision",
        fitConfigurationProvider.getFitConfiguration().getPrecisionThreshold(), 2,
        new OptionListener<Double>() {
          @Override
          public boolean collectOptions(Double field) {
            fitConfigurationProvider.getFitConfiguration().setPrecisionThreshold(field);
            return collectOptions(false);
          }

          @Override
          public boolean collectOptions() {
            return collectOptions(true);
          }

          private boolean collectOptions(boolean silent) {
            final FitConfiguration localFitConfig = fitConfigurationProvider.getFitConfiguration();
            final ExtendedGenericDialog egd = new ExtendedGenericDialog("Precision Options", null);
            final int oldIndex = localFitConfig.getPrecisionMethod().ordinal();
            egd.addChoice("Precision_method", SettingsManager.getPrecisionMethodNames(), oldIndex);
            egd.setSilent(silent);
            egd.showDialog(true, gd);
            if (egd.wasCanceled()) {
              return false;
            }
            final int newIndex = egd.getNextChoiceIndex();
            localFitConfig.setPrecisionMethod(newIndex);
            return oldIndex != newIndex;
          }
        });
  }

  private void log(String format, Object... args) {
    if (!silent) {
      ImageJUtils.log(format, args);
    }
  }

  private int showSimpleDialog() {
    // Just support circular fitting
    fitConfig.setPsf(PsfProtosHelper.defaultOneAxisGaussian2DPSF);
    fitConfig.setFixedPsf(false);

    // TODO - Support sCMOS camera. This may be 'too difficult' as the
    // user will need to have created a per-pixel calibration image

    final CalibrationWriter calibration = fitConfig.getCalibrationWriter();
    final boolean requireCalibration = requireCalibration(calibration);
    if (requireCalibration && !showCalibrationWizard(calibration, true)) {
      return DONE;
    }

    // Present dialog with simple output options: Image, Table
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("simple-fit"));
    gd.addMessage("Fit single-molecule localisations");

    if (!requireCalibration) {
      gd.addCheckbox("Use_current_calibration", true);
    }
    gd.addCheckbox("Show_table", settings.showTable);
    gd.addCheckbox("Show_image", settings.showImage);

    gd.showDialog();
    if (gd.wasCanceled()) {
      return DONE;
    }

    boolean useCurrentCalibration = true;
    if (!requireCalibration) {
      useCurrentCalibration = gd.getNextBoolean();
    }
    settings.showTable = gd.getNextBoolean();
    settings.showImage = gd.getNextBoolean();

    if (!useCurrentCalibration && !showCalibrationWizard(calibration, false)) {
      return DONE;
    }

    // Restore fitting to default settings but maintain the calibrated width
    final double sd = fitConfig.getInitialXSd();
    config = new FitEngineConfiguration();
    fitConfig = config.getFitConfiguration();
    fitConfig.setInitialPeakStdDev(sd);
    // Allow to move 1 SD
    fitConfig.setCoordinateShiftFactor(1);
    resultsSettings = ResultsSettings.newBuilder();

    // Do simple results output. We only need to set non-default values.
    resultsSettings.getResultsInMemorySettingsBuilder().setInMemory(true);
    if (settings.showTable) {
      final ResultsTableSettings.Builder tableSettings =
          resultsSettings.getResultsTableSettingsBuilder();
      tableSettings.setShowTable(true);
    }
    if (settings.showImage) {
      final ResultsImageSettings.Builder imageSettings =
          resultsSettings.getResultsImageSettingsBuilder();
      imageSettings.setImageType(ResultsImageType.DRAW_INTENSITY);
      imageSettings.setScale(Math.ceil(1024.0 / Math.max(bounds.width, bounds.height)));
      imageSettings.setWeighted(true);
      imageSettings.setEqualised(true);
    }

    // Log the settings we care about:
    IJ.log(LOG_SPACER);
    IJ.log("Peak Fit");
    IJ.log(LOG_SPACER);
    ImageJUtils.log("Pixel pitch = %s", MathUtils.rounded(calibration.getNmPerPixel(), 4));
    ImageJUtils.log("Exposure Time = %s", MathUtils.rounded(calibration.getExposureTime(), 4));
    ImageJUtils.log("PSF width = %s", MathUtils.rounded(fitConfig.getInitialXSd(), 4));

    // Save
    fitConfig.setCalibration(calibration.getCalibration());
    saveFitEngineSettings();
    SettingsManager.writeSettings(resultsSettings.build());

    return FLAGS;
  }

  private boolean saveFitEngineSettings() {
    return saveFitEngineSettings(config);
  }

  private static boolean saveFitEngineSettings(FitEngineConfiguration config) {
    return SettingsManager.writeSettings(config, 0);
  }

  /**
   * Check if additional calibration information is required.
   *
   * <p>Check the calibration is valid for fitting.
   *
   * @param calibration the calibration
   * @return True if additional calibration information is required, false if the system is
   *         calibrated
   */
  private boolean requireCalibration(CalibrationWriter calibration) {
    // Check for a supported camera
    if (!calibration.isCcdCamera()) {
      return true;
    }
    // Check if the calibration contains: Pixel pitch, Gain (can be 1), Exposure time
    if (!calibration.hasNmPerPixel()) {
      return true;
    }
    // Bias can be zero (but this is unlikely)
    if (!calibration.hasCountPerPhoton() || calibration.getBias() <= 0) {
      return true;
    }
    if (!calibration.hasExposureTime()) {
      return true;
    }

    // Check for a PSF width
    return fitConfig.getInitialXSd() <= 0;
  }

  private boolean showCalibrationWizard(CalibrationWriter calibration, boolean showIntroduction) {
    if (showIntroduction) {
      // Currently we do not support a sCMOS camera.
      // If this is set then the settings must have been created in Peak Fit and not
      // the calibration wizard.
      final ArrayList<String> msgs = new ArrayList<>();
      if (calibration.isScmos()) {
        msgs.add("The configuration file specifies a sCMOS camera. "
            + "This is only supported in the Peak Fit plugin.");
        msgs.add("Continuing with the wizard will overwrite the current cofiguration.");
        msgs.add("Press Cancel to close the wizard. "
            + "You can then run Peak Fit with the current configuration.");
      } else {
        msgs.add("No configuration file could be loaded.");
        msgs.add("Please follow the configuration wizard to calibrate.");
      }

      final ExtendedGenericDialog gd = newWizardDialog(msgs.toArray(new String[0]));
      gd.showDialog();
      if (gd.wasCanceled()) {
        return false;
      }
    }

    if (!getCameraType(calibration)) {
      return false;
    }
    if (!getPixelPitch(calibration)) {
      return false;
    }
    if (!getGain(calibration)) {
      return false;
    }
    if (!getExposureTime(calibration)) {
      return false;
    }
    if (!getPeakWidth(calibration)) {
      return false;
    }

    // Check parameters
    try {
      ParameterUtils.isAboveZero("nm per pixel", calibration.getNmPerPixel());
      // We allow bias to be zero
      ParameterUtils.isAboveZero("Gain", calibration.getCountPerPhoton());
      ParameterUtils.isAboveZero("Exposure time", calibration.getExposureTime());
      ParameterUtils.isAboveZero("Initial SD", fitConfig.getInitialXSd());
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return true;
  }

  private static ExtendedGenericDialog newWizardDialog(String... messages) {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("simple-fit"));
    final String header = "-=-";
    gd.addMessage(header + " " + TITLE + " Configuration Wizard " + header);
    for (final String message : messages) {
      gd.addMessage(TextUtils.wrap(message, 80));
    }
    return gd;
  }

  private static boolean getCameraType(CalibrationWriter calibration) {
    final ExtendedGenericDialog gd = newWizardDialog("Enter the type of camera.");
    gd.addChoice("Camera_type", SettingsManager.getCameraTypeNames(),
        CalibrationProtosHelper.getName(calibration.getCameraType()));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }
    calibration.setCameraType(SettingsManager.getCameraTypeValues()[gd.getNextChoiceIndex()]);
    if (!calibration.isCcdCamera()) {
      // TODO - Support sCMOS camera

      IJ.error("Unsupported camera type "
          + CalibrationProtosHelper.getName(calibration.getCameraType()));
      return false;
    }
    return true;
  }

  private static boolean getPixelPitch(CalibrationWriter calibration) {
    final ExtendedGenericDialog gd = newWizardDialog(
        "Enter the size of each pixel. This is required to ensure the dimensions of the image"
            + " are calibrated.",
        "E.g. a camera with a 6.45um pixel size and a 60x objective will have a pitch of"
            + " 6450/60 = 107.5nm.");
    // TODO - Add a pop-up calculator...
    gd.addNumericField("Calibration", calibration.getNmPerPixel(), 2, 6, "nm/px");
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }
    calibration.setNmPerPixel(Math.abs(gd.getNextNumber()));
    return true;
  }

  private static boolean getGain(CalibrationWriter calibration) {
    final ExtendedGenericDialog gd = newWizardDialog("Enter the bias and total gain.",
        "This is usually supplied with your camera certificate. The bias is a fixed offset added"
            + " to the camera counts. The gain indicates how many Analogue-to-Digital-Units (Count)"
            + " are recorded at the pixel for each photon registered on the sensor.",
        "The gain is usually expressed using the product of the EM-gain (if applicable), the"
            + " camera gain and the sensor quantum efficiency.",
        "A value of 1 means no conversion to photons will occur.");
    // TODO - Add a wizard to allow calculation of total gain from EM-gain, camera gain and QE
    gd.addNumericField("Camera_bias", calibration.getBias(), 2, 6, "Count");
    gd.addNumericField("Gain", calibration.getCountPerPhoton(), 2, 6, "Count/photon");
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }
    calibration.setBias(Math.abs(gd.getNextNumber()));
    calibration.setCountPerPhoton(Math.abs(gd.getNextNumber()));
    return true;
  }

  private static boolean getExposureTime(CalibrationWriter calibration) {
    final ExtendedGenericDialog gd = newWizardDialog(
        "Enter the exposure time. Calibration of the exposure time allows correct reporting of on"
            + " and off times.",
        "This is the length of time for each frame in the image.");
    gd.addNumericField("Exposure_time", calibration.getExposureTime(), 2, 6, "ms");
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }
    calibration.setExposureTime(Math.abs(gd.getNextNumber()));
    return true;
  }

  private boolean getPeakWidth(final CalibrationWriter calibration) {
    final ExtendedGenericDialog gd = newWizardDialog("Enter the expected peak width in pixels.",
        "A point source of light will not be focussed perfectly by the microscope but will appear"
            + " as a spread out peak. This Point Spread Function (PSF) can be modelled using a"
            + " 2D Gaussian curve.",
        "An optimised optical system (lens and camera sensor) should have a peak standard"
            + " deviation of approximately 1 pixel when in focus. This allows the fitting routine"
            + " to have enough data to identify the centre of the peak without spreading the light"
            + " over too many pixels (which increases noise).",
        "The peak width can be estimated using the wavelength of light emitted by the single"
            + " molecules and the parameters of the microscope. Use a PSF calculator by clicking"
            + " the checkbox below:");
    gd.addNumericField("Gaussian_SD", fitConfig.getInitialXSd(), 3);
    // Add ability to run the PSF Calculator to get the width
    if (ImageJUtils.isShowGenericDialog()) {
      final TextField textInitialPeakStdDev0 = (TextField) gd.getNumericFields().get(0);
      gd.addAndGetButton("Run PSF calculator", event -> {
        // Run the PSF Calculator
        calculatorSettings =
            calculatorSettings.toBuilder().setPixelPitch(calibration.getNmPerPixel() / 1000.0)
                .setMagnification(1).setBeamExpander(1).build();
        final double sd = new PsfCalculator().calculate(calculatorSettings, true);
        if (sd > 0) {
          textInitialPeakStdDev0.setText(Double.toString(sd));
        }
      });
    }
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }
    fitConfig.setInitialPeakStdDev(Math.abs(gd.getNextNumber()));
    return true;
  }

  private boolean readDialog(ExtendedGenericDialog gd, boolean isCrop) {
    // Ignore the template
    gd.getNextChoice();

    final CalibrationWriter calibration = fitConfig.getCalibrationWriter();
    calibration.setCameraType(SettingsManager.getCameraTypeValues()[gd.getNextChoiceIndex()]);
    calibration.setNmPerPixel(Math.abs(gd.getNextNumber()));
    calibration.setExposureTime(Math.abs(gd.getNextNumber()));
    fitConfig.setCalibration(calibration.getCalibration());

    // Note: The bias and read noise will just end up being what was in the configuration file.
    // One fix for this is to save/load only the settings that are required from the configuration
    // file (the others will remain unchanged). This will require a big refactor of the settings
    // save/load.
    // The simple fix is to create a plugin to allow the configuration to be changed for results.
    if (isCrop) {
      ignoreBoundsForNoise = extraSettings.optionIgnoreBoundsForNoise = gd.getNextBoolean();
    }

    fitConfig.setPsfType(PeakFit.getPsfTypeValues()[gd.getNextChoiceIndex()]);
    config.setDataFilterType(gd.getNextChoiceIndex());
    // Note: The absolute flag is set in extra options
    config.setDataFilter(gd.getNextChoiceIndex(), Math.abs(gd.getNextNumber()), 0);
    config.setSearch(gd.getNextNumber());
    config.setBorder(gd.getNextNumber());
    config.setFitting(gd.getNextNumber());
    if (extraOptions && !fitMaxima) {
      extraSettings.interlacedData = gd.getNextBoolean();
      extraSettings.integrateFrames = (int) gd.getNextNumber();
    }

    if (!maximaIdentification) {
      fitConfig.setFitSolver(gd.getNextChoiceIndex());
      if (extraOptions) {
        fitConfig.setBackgroundFitting(gd.getNextBoolean());
      }
      config.setFailuresLimit((int) gd.getNextNumber());
      config.setPassRate(gd.getNextNumber());
      config.setIncludeNeighbours(gd.getNextBoolean());
      config.setNeighbourHeightThreshold(gd.getNextNumber());
      config.setResidualsThreshold(gd.getNextNumber());

      config.setDuplicateDistance(gd.getNextNumber());

      fitConfig.setSmartFilter(gd.getNextBoolean());
      fitConfig.setDisableSimpleFilter(gd.getNextBoolean());
      fitConfig.setCoordinateShiftFactor(gd.getNextNumber());
      fitConfig.setSignalStrength(gd.getNextNumber());
      fitConfig.setMinPhotons(gd.getNextNumber());
      if (extraOptions) {
        fitConfig.setNoise(gd.getNextNumber());
        config.setNoiseMethod(gd.getNextChoiceIndex());
      }
      fitConfig.setMinWidthFactor(gd.getNextNumber());
      fitConfig.setMaxWidthFactor(gd.getNextNumber());
      fitConfig.setPrecisionThreshold(gd.getNextNumber());
    }

    resultsSettings.setLogProgress(gd.getNextBoolean());
    if (!maximaIdentification) {
      resultsSettings.setShowDeviations(gd.getNextBoolean());
    }

    resultsSettings.getResultsTableSettingsBuilder().setShowTable(gd.getNextBoolean());
    resultsSettings.getResultsImageSettingsBuilder().setImageTypeValue(gd.getNextChoiceIndex());
    if (extraOptions) {
      extraSettings.showProcessedFrames = gd.getNextBoolean();
    }
    resultsSettings.getResultsFileSettingsBuilder().setFileFormatValue(gd.getNextChoiceIndex());
    resultsSettings.getResultsFileSettingsBuilder().setResultsDirectory(gd.getNextString());
    resultsSettings.getResultsInMemorySettingsBuilder().setInMemory(gd.getNextBoolean());
    if (extraOptions) {
      settings.fractionOfThreads = Math.abs(gd.getNextNumber());
    }

    gd.collectOptions();

    // Save to allow dialog state to be maintained even with invalid parameters
    saveFitEngineSettings();
    SettingsManager.writeSettings(resultsSettings.build());

    if (gd.invalidNumber()) {
      return false;
    }

    // Check arguments
    try {
      // No check on camera calibration. This is left to the FitConfiguration to
      // error if the settings are incorrect

      ParameterUtils.isAboveZero("nm per pixel", calibration.getNmPerPixel());
      ParameterUtils.isAboveZero("Exposure time", calibration.getExposureTime());
      if (fitConfig.getPsfTypeValue() != PSFType.ASTIGMATIC_GAUSSIAN_2D_VALUE) {
        ParameterUtils.isAboveZero("Initial SD0", fitConfig.getInitialXSd());
        if (fitConfig.getPsf().getParametersCount() > 1) {
          ParameterUtils.isAboveZero("Initial SD1", fitConfig.getInitialYSd());
        }
      }
      ParameterUtils.isAboveZero("Search_width", config.getSearch());
      ParameterUtils.isAboveZero("Fitting_width", config.getFitting());
      if (extraOptions && !fitMaxima) {
        ParameterUtils.isPositive("Integrate frames", extraSettings.integrateFrames);
      }
      if (!maximaIdentification) {
        // This can be negative to disable, i.e. fit everything
        // Parameters.isPositive("Failures limit", config.getFailuresLimit())
        ParameterUtils.isPositive("Neighbour height threshold",
            config.getNeighbourHeightThreshold());
        ParameterUtils.isPositive("Residuals threshold", config.getResidualsThreshold());
        ParameterUtils.isPositive("Duplicate distance", config.getDuplicateDistance());

        if (!fitConfig.isSmartFilter()) {
          ParameterUtils.isPositive("Coordinate Shift factor",
              fitConfig.getCoordinateShiftFactor());
          ParameterUtils.isPositive("Signal strength", fitConfig.getSignalStrength());
          ParameterUtils.isPositive("Min photons", fitConfig.getMinPhotons());
        }
        if (extraOptions) {
          ParameterUtils.isPositive("Noise", fitConfig.getNoise());
        }
        if (!fitConfig.isSmartFilter()) {
          ParameterUtils.isPositive("Min width factor", fitConfig.getMinWidthFactor());
          ParameterUtils.isPositive("Width factor", fitConfig.getMaxWidthFactor());
          ParameterUtils.isPositive("Precision threshold", fitConfig.getPrecisionThreshold());
          if (fitConfig.getPrecisionThreshold() > 0) {
            if (fitConfig.getPrecisionMethod() == PrecisionMethod.PRECISION_METHOD_NA) {
              throw new IllegalArgumentException("Precision filter requires a precision method");
            }
            if (fitConfig.isPrecisionUsingBackground() && calibration.isCcdCamera()
                && (calibration.getBias() == 0 || !calibration.hasCountPerPhoton())) {
              throw new IllegalArgumentException(
                  "Precision using the local background requires the camera bias");
            }
          }
        }
      }
      final ResultsImageSettings.Builder imageSettings =
          resultsSettings.getResultsImageSettingsBuilder();
      if (imageSettings.getImageType() == ResultsImageType.DRAW_INTENSITY_AVERAGE_PRECISION
          || imageSettings
              .getImageType() == ResultsImageType.DRAW_LOCALISATIONS_AVERAGE_PRECISION) {
        ParameterUtils.isAboveZero("Image precision", imageSettings.getAveragePrecision());
      }
      ParameterUtils.isAboveZero("Image scale", imageSettings.getScale());
      if (extraOptions) {
        ParameterUtils.isPositive("Image rolling window", imageSettings.getRollingWindowSize());
      }
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    final int flags = (extraOptions) ? FLAG_EXTRA_OPTIONS : 0;

    // If precision filtering then we need the camera bias
    if (!maximaIdentification) {
      if (!configurePsfModel(config, flags)) {
        return false;
      }

      if (!configureResultsFilter(config, flags)) {
        return false;
      }
    }

    if (!configureDataFilter(config, flags)) {
      return false;
    }

    // Second dialog for solver dependent parameters
    if (!maximaIdentification && !configureFitSolver(config, source.getBounds(), bounds, flags)) {
      return false;
    }

    // Extra parameters are needed for interlaced data
    if (extraSettings.interlacedData) {
      gd = new ExtendedGenericDialog(TITLE);
      gd.addMessage("Interlaced data requires a repeating pattern of frames to process.\n"
          + "Describe the regular repeat of the data:\n \n"
          + "Start = The first frame that contains data\n"
          + "Block = The number of continuous frames containing data\n"
          + "Skip = The number of continuous frames to ignore before the next data\n \n"
          + "E.G. 2:9:1 = Data was imaged from frame 2 for 9 frames, 1 frame to ignore,"
          + " then repeat.");
      gd.addNumericField("Start", extraSettings.dataStart, 0);
      gd.addNumericField("Block", extraSettings.dataBlock, 0);
      gd.addNumericField("Skip", extraSettings.dataSkip, 0);

      gd.showDialog();
      if (gd.wasCanceled()) {
        return false;
      }

      if (!gd.wasCanceled()) {
        extraSettings.dataStart = (int) gd.getNextNumber();
        extraSettings.dataBlock = (int) gd.getNextNumber();
        extraSettings.dataSkip = (int) gd.getNextNumber();

        if (extraSettings.dataStart > 0 && extraSettings.dataBlock > 0
            && extraSettings.dataSkip > 0) {
          // Store options for next time
          extraSettings.save();
        }
      } else {
        extraSettings.interlacedData = false;
      }
    }

    final boolean result = saveFitEngineSettings();
    if (!result) {
      IJ.error(TITLE, "Failed to save settings");
    }

    return result;
  }

  /**
   * Show a dialog to configure the PSF model. The updated settings are saved to the settings file.
   *
   * <p>If the configuration is for a 3D PSF then a dialog to configure the z model is shown.
   *
   * @param config the config
   * @return true, if successful
   */
  public static boolean configurePsfModel(FitEngineConfiguration config) {
    return configurePsfModel(config, FLAG_NO_SAVE);
  }

  /**
   * Show a dialog to configure the PSF model. The updated settings are saved to the settings file.
   *
   * <p>If the configuration is for a 3D PSF then a dialog to configure the z model is shown.
   *
   * @param config the config
   * @param flags the flags
   * @return true, if successful
   */
  public static boolean configurePsfModel(FitEngineConfiguration config, int flags) {
    final FitConfiguration fitConfig = config.getFitConfiguration();
    if (fitConfig.getPsfTypeValue() != PSFType.ASTIGMATIC_GAUSSIAN_2D_VALUE) {
      return true;
    }

    // Get the astigmatism z-model
    final AstigmatismModel model = AstigmatismModelManager.getModel(fitConfig.getPsfModelName());
    if (model == null) {
      IJ.error(TITLE, "Failed to load the model: " + fitConfig.getPsfModelName());
      return false;
    }

    // Conversion to the correct units in pixels is done within the FitConfiguration object.
    fitConfig.setAstigmatismModel(model);

    if (BitFlagUtils.anyNotSet(flags, FLAG_NO_SAVE)) {
      SettingsManager.writeSettings(config, 0);
    }

    return true;
  }

  /**
   * Show a dialog to configure the results filter. The updated settings are saved to the settings
   * file.
   *
   * <p>If the configuration is for a 3D PSF then a dialog to configure the z range for the results
   * is shown (see {@link #configureZFilter(FitEngineConfiguration, int)}).
   *
   * <p>If the configuration is for a smart filter then a dialog to configure the smart filter is
   * shown (see {@link #configureSmartFilter(FitEngineConfiguration, int)}).
   *
   * @param config the config
   * @param flags the flags
   * @return true, if successful
   */
  public static boolean configureResultsFilter(FitEngineConfiguration config, int flags) {
    boolean result = configureZFilter(config, flags);
    result = result && configureSmartFilter(config, flags);
    return result;
  }

  /**
   * Show a dialog to configure the results z filter. The updated settings are saved to the settings
   * file.
   *
   * <p>If the fit configuration PSF is not 3D or the simple filter is disabled then this method
   * returns true. If it is enabled then a dialog is shown to input the configuration for the z
   * filter.
   *
   * <p>Note: The PSF and any z-model must be correctly configured for fitting in pixel units.
   *
   * @param config the config
   * @param flags the flags
   * @return true, if successful
   */
  public static boolean configureZFilter(FitEngineConfiguration config, int flags) {
    final FitConfiguration fitConfig = config.getFitConfiguration();
    if (fitConfig.isDisableSimpleFilter() || !fitConfig.is3D()) {
      return true;
    }

    // Create a converter to map the model units in pixels to nm for the dialog.
    // Note the output units of pixels may not yet be set in the calibration so we assume it is
    // pixels.
    final TypeConverter<DistanceUnit> c = UnitConverterUtils.createConverter(DistanceUnit.PIXEL,
        DistanceUnit.NM, fitConfig.getCalibrationReader().getNmPerPixel());

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

    gd.addMessage("3D filter");
    gd.addNumericField("Min_z", c.convert(fitConfig.getMinZ()), 0, 6, "nm");
    gd.addNumericField("Max_z", c.convert(fitConfig.getMaxZ()), 0, 6, "nm");

    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }

    final double minZ = gd.getNextNumber();
    final double maxZ = gd.getNextNumber();

    if (gd.invalidNumber() || minZ > maxZ) {
      IJ.error(TITLE, "Min Z must be equal or below the max Z");
      return false;
    }

    // Map back
    fitConfig.setMinZ(c.convertBack(minZ));
    fitConfig.setMaxZ(c.convertBack(maxZ));

    if (BitFlagUtils.anyNotSet(flags, FLAG_NO_SAVE)) {
      SettingsManager.writeSettings(config, 0);
    }

    return true;
  }

  /**
   * Show a dialog to configure the smart filter. The updated settings are optionally saved to the
   * settings file.
   *
   * <p>If the fit configuration isSmartFilter is not enabled then this method returns true. If it
   * is enabled then a dialog is shown to input the configuration for a smart filter. If no valid
   * filter can be created from the input then the method returns false.
   *
   * <p>Note: If the smart filter is successfully configured then the user may want to disable the
   * standard fit validation.
   *
   * @param config the config
   * @param flags the flags
   * @return true, if successful
   */
  public static boolean configureSmartFilter(FitEngineConfiguration config, int flags) {
    final FitConfiguration fitConfig = config.getFitConfiguration();
    if (!fitConfig.isSmartFilter()) {
      return true;
    }
    final boolean result = configureSmartFilter(fitConfig);
    if (result) {
      if (BitFlagUtils.anyNotSet(flags, FLAG_NO_SAVE)) {
        SettingsManager.writeSettings(config, 0);
      }
    }
    return result;
  }

  /**
   * Show a dialog to configure the smart filter.
   *
   * <p>If the fit configuration isSmartFilter is not enabled then this method returns true. If it
   * is enabled then a dialog is shown to input the configuration for a smart filter. If no valid
   * filter can be created from the input then the method returns false.
   *
   * <p>Note: If the smart filter is successfully configured then the user may want to disable the
   * standard fit validation.
   *
   * @param fitConfig the fit config
   * @return true, if successful
   */
  public static boolean configureSmartFilter(final FitConfiguration fitConfig) {
    if (!fitConfig.isSmartFilter()) {
      return true;
    }

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

    String xml = fitConfig.getSmartFilterString();
    if (TextUtils.isNullOrEmpty(xml)) {
      xml = fitConfig.getDefaultSmartFilterXml();
    }

    gd.addMessage("Smart filter (used to pick optimum results during fitting)");
    gd.addTextAreas(uk.ac.sussex.gdsc.core.utils.XmlUtils.convertQuotes(xml), null, 8, 60);
    // Add message about precision filtering
    gd.addMessage(
        TextUtils.wrap("Note: Smart filters using precision may require a local background level. "
            + "Ensure the camera calibration is correct including any bias.", 80));

    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }

    xml = gd.getNextText();
    final Filter f = Filter.fromXml(xml);
    if (!(f instanceof DirectFilter)) {
      return false;
    }

    fitConfig.setDirectFilter((DirectFilter) f);

    return true;
  }

  /**
   * Show a dialog to configure the data filter. The data filter type and the first data filter must
   * ALREADY be set in the configuration. The subsequent filters are then configured, e.g. for
   * difference and jury filters.
   *
   * <p>The updated settings are saved to the settings file. An error message is shown if the dialog
   * is cancelled or the configuration is invalid.
   *
   * <p>If the configuration is for a per-pixel camera type (e.g. sCMOS) then the camera model will
   * be loaded using the configured camera model name. This will be used to validate the filter to
   * check the filter supports the per-pixel camera type.
   *
   * @param config the config
   * @param flags the flags
   * @return True if the configuration succeeded
   */
  public static boolean configureDataFilter(final FitEngineConfiguration config, int flags) {
    int numberOfFilters = 1;
    int filterCount;
    switch (config.getDataFilterType()) {
      case JURY:
        filterCount = Integer.MAX_VALUE;
        break;

      case DIFFERENCE:
        filterCount = 2;
        break;

      case SINGLE:
      default:
        filterCount = 1;
    }

    final String[] filterNames = SettingsManager.getDataFilterMethodNames();
    final DataFilterMethod[] filterValues = SettingsManager.getDataFilterMethodValues();

    // We use the previous value in the event the configuration does not have any current values.
    // Check we have at least the first filter.
    if (config.getDataFiltersCount() == 0) {
      throw new IllegalStateException("No primary filter is configured");
    }

    final FitEngineConfigurationProvider fitEngineConfigurationProvider = () -> config;

    for (int i = 1; i < filterCount; i++) {
      final int filter = i + 1;
      final int ii = i;

      final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
      if (filter == filterCount) {
        // This is maximum filter count so no continue option
        ImageJUtils.addMessage(gd, "Configure the %s filter.",
            FitProtosHelper.getName(config.getDataFilterType()));
      } else {
        gd.enableYesNoCancel("Add", "Continue");
        ImageJUtils.addMessage(gd,
            "Configure the %s filter.\nClick continue to proceed with the current set of %d.",
            FitProtosHelper.getName(config.getDataFilterType()), i);
      }

      final String fieldName = "Spot_filter" + filter;
      if (IJ.isMacro()) {
        // Use blank default value so bad macro parameters return nothing
        gd.addStringField(fieldName, "");
      } else {
        gd.addChoice(fieldName, filterNames, filterNames[config
            .getDataFilterMethod(ii, config.getDataFilterMethod(ii - 1)).ordinal()]);
      }
      addRelativeParameterOptions(gd, new RelativeParameterProvider(0, 4.5, "Smoothing" + filter,
          fitEngineConfigurationProvider, true) {
        @Override
        void setAbsolute(boolean absolute) {
          // Get the current settings
          final FitEngineConfiguration c =
              fitEngineConfigurationProvider.getFitEngineConfiguration();
          final DataFilterMethod m = c.getDataFilterMethod(ii);
          final double smooth = c.getDataFilterParameter(ii).getValue();
          // Reset with the new absolute value
          c.setDataFilter(m, smooth, absolute, ii);
        }

        @Override
        boolean isAbsolute() {
          final FitEngineConfiguration c =
              fitEngineConfigurationProvider.getFitEngineConfiguration();
          return c.getDataFilterParameterAbsolute(ii, c.getDataFilterParameterAbsolute(ii - 1));
        }

        @Override
        double getValue() {
          final FitEngineConfiguration c =
              fitEngineConfigurationProvider.getFitEngineConfiguration();
          return c.getDataFilterParameterValue(ii, c.getDataFilterParameterValue(ii - 1));
        }
      });
      gd.showDialog();
      if (gd.wasCanceled()) {
        return false;
      }
      if (gd.wasOKed()) {
        int filterIndex = -1;
        if (IJ.isMacro()) {
          final String filterName = gd.getNextString();
          for (int j = 0; j < filterNames.length; j++) {
            if (filterNames[j].equals(filterName)) {
              filterIndex = j;
              break;
            }
          }

          if (filterIndex < 0) {
            break;
          }
        } else {
          filterIndex = gd.getNextChoiceIndex();
        }
        // Note: The absolute flag is set in extra options
        config.setDataFilter(filterValues[filterIndex], Math.abs(gd.getNextNumber()), i);
        gd.collectOptions();
        numberOfFilters++;
      } else {
        break;
      }
    }
    config.setNumberOfFilters(numberOfFilters);

    if (BitFlagUtils.anyNotSet(flags, FLAG_NO_SAVE)) {
      saveFitEngineSettings(config);
    }

    final FitConfiguration fitConfig = config.getFitConfiguration();
    final CalibrationReader calibration = fitConfig.getCalibrationReader();
    if (calibration.isScmos()) {
      fitConfig.setCameraModel(CameraModelManager.load(fitConfig.getCameraModelName()));
    }

    try {
      config.createSpotFilter();
    } catch (final IllegalStateException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return true;
  }

  /**
   * Show a dialog to configure the fit solver. The updated settings are saved to the settings file.
   * An error message is shown if the dialog is cancelled or the configuration is invalid.
   *
   * <p>The bounds are used to validate the camera model. The camera model must be large enough to
   * cover the source bounds. If larger then it will be cropped. Optionally an internal region of
   * the input image can be specifed. This is relative to the width and height of the input image.
   * If no camera model is present then the bounds can be null.
   *
   * @param config the config
   * @param sourceBounds the source image bounds (used to validate the camera model dimensions)
   * @param bounds the crop bounds (relative to the input image, used to validate the camera model
   *        dimensions)
   * @param flags the flags
   * @return True if the configuration succeeded
   */
  public static boolean configureFitSolver(FitEngineConfiguration config, Rectangle sourceBounds,
      Rectangle bounds, int flags) {
    final boolean extraOptions = BitFlagUtils.anySet(flags, FLAG_EXTRA_OPTIONS);
    final boolean ignoreCalibration = BitFlagUtils.anySet(flags, FLAG_IGNORE_CALIBRATION);
    final boolean saveSettings = BitFlagUtils.anyNotSet(flags, FLAG_NO_SAVE);

    final FitConfiguration fitConfig = config.getFitConfiguration();
    final CalibrationWriter calibration = fitConfig.getCalibrationWriter();

    final FitSolver fitSolver = fitConfig.getFitSolver();

    final boolean isLvm = fitSolver == FitSolver.LVM_LSE || fitSolver == FitSolver.LVM_WLSE
        || fitSolver == FitSolver.LVM_MLE;
    final boolean isFastMml =
        fitSolver == FitSolver.FAST_MLE || fitSolver == FitSolver.BACKTRACKING_FAST_MLE;
    final boolean isSteppingFunctionSolver = isLvm || isFastMml;

    if (fitSolver == FitSolver.MLE) {
      final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
      if (!ignoreCalibration) {
        gd.addMessage("Maximum Likelihood Estimation requires CCD-type camera parameters");
        gd.addNumericField("Camera_bias", calibration.getBias(), 2, 6, "count");
        gd.addCheckbox("Model_camera_noise", fitConfig.isModelCamera());
        gd.addNumericField("Read_noise", calibration.getReadNoise(), 2, 6, "count");
        gd.addNumericField("Quantum_efficiency", calibration.getQuantumEfficiency(), 2, 6,
            "electron/photon");
        gd.addCheckbox("EM-CCD", calibration.isEmCcd());
      } else {
        gd.addMessage("Maximum Likelihood Estimation requires additional parameters");
      }
      // This works because the proto configuration enum matches the named enum
      final String[] searchNames =
          SettingsManager.getNames((Object[]) MaximumLikelihoodFitter.SearchMethod.values());
      gd.addChoice("Search_method", searchNames,
          searchNames[fitConfig.getSearchMethod().getNumber()]);
      gd.addStringField("Relative_threshold", MathUtils.rounded(fitConfig.getRelativeThreshold()));
      gd.addStringField("Absolute_threshold", MathUtils.rounded(fitConfig.getAbsoluteThreshold()));
      gd.addNumericField("Max_iterations", fitConfig.getMaxIterations(), 0);
      gd.addNumericField("Max_function_evaluations", fitConfig.getMaxFunctionEvaluations(), 0);
      if (extraOptions) {
        gd.addCheckbox("Gradient_line_minimisation", fitConfig.isGradientLineMinimisation());
      }
      gd.showDialog();
      if (gd.wasCanceled()) {
        return false;
      }
      if (!ignoreCalibration) {
        calibration.setBias(Math.abs(gd.getNextNumber()));
        fitConfig.setModelCamera(gd.getNextBoolean());
        calibration.setReadNoise(Math.abs(gd.getNextNumber()));
        calibration.setQuantumEfficiency(Math.abs(gd.getNextNumber()));
        calibration.setCameraType((gd.getNextBoolean()) ? CameraType.EMCCD : CameraType.CCD);
        fitConfig.setCalibration(calibration.getCalibration());
      }
      fitConfig.setSearchMethod(gd.getNextChoiceIndex());
      fitConfig.setRelativeThreshold(getThresholdNumber(gd));
      fitConfig.setAbsoluteThreshold(getThresholdNumber(gd));
      fitConfig.setMaxIterations((int) gd.getNextNumber());
      fitConfig.setMaxFunctionEvaluations((int) gd.getNextNumber());
      if (extraOptions) {
        fitConfig.setGradientLineMinimisation(gd.getNextBoolean());
      } else {
        // This option is for the Conjugate Gradient optimiser and makes it less stable
        fitConfig.setGradientLineMinimisation(false);
      }

      if (saveSettings) {
        saveFitEngineSettings(config);
      }

      try {
        ParameterUtils.isAboveZero("Relative threshold", fitConfig.getRelativeThreshold());
        ParameterUtils.isAboveZero("Absolute threshold", fitConfig.getAbsoluteThreshold());
        ParameterUtils.isAboveZero("Max iterations", fitConfig.getMaxIterations());
        ParameterUtils.isAboveZero("Max function evaluations",
            fitConfig.getMaxFunctionEvaluations());
        fitConfig.getFunctionSolver();
      } catch (final IllegalArgumentException | IllegalStateException ex) {
        IJ.error(TITLE, ex.getMessage());
        return false;
      }
    } else if (isSteppingFunctionSolver) {
      final boolean requireCalibration = !ignoreCalibration && fitSolver != FitSolver.LVM_LSE;

      // Collect options for LVM fitting
      final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
      final String fitSolverName = FitProtosHelper.getName(fitSolver);
      gd.addMessage(fitSolverName + " requires additional parameters");
      gd.addStringField("Relative_threshold", MathUtils.rounded(fitConfig.getRelativeThreshold()));
      gd.addStringField("Absolute_threshold", MathUtils.rounded(fitConfig.getAbsoluteThreshold()));
      gd.addStringField("Parameter_relative_threshold",
          MathUtils.rounded(fitConfig.getParameterRelativeThreshold()));
      gd.addStringField("Parameter_absolute_threshold",
          MathUtils.rounded(fitConfig.getParameterAbsoluteThreshold()));
      gd.addNumericField("Max_iterations", fitConfig.getMaxIterations(), 0);
      if (isLvm) {
        gd.addNumericField("Lambda", fitConfig.getLambda(), 4);
      }
      if (isFastMml) {
        gd.addCheckbox("Fixed_iterations", fitConfig.isFixedIterations());
        // This works because the proto configuration enum matches the named enum
        final String[] lineSearchNames = SettingsManager
            .getNames((Object[]) FastMleSteppingFunctionSolver.LineSearchMethod.values());
        gd.addChoice("Line_search_method", lineSearchNames,
            lineSearchNames[fitConfig.getLineSearchMethod().getNumber()]);
      }

      gd.addCheckbox("Use_clamping", fitConfig.isUseClamping());
      gd.addCheckbox("Dynamic_clamping", fitConfig.isUseDynamicClamping());
      final PSF psf = fitConfig.getPsf();
      final boolean isAstigmatism = psf.getPsfType() == PSFType.ASTIGMATIC_GAUSSIAN_2D;
      final int nParams = PsfHelper.getParameterCount(psf);
      if (extraOptions) {
        gd.addNumericField("Clamp_background", fitConfig.getClampBackground(), 2);
        gd.addNumericField("Clamp_signal", fitConfig.getClampSignal(), 2);
        gd.addNumericField("Clamp_x", fitConfig.getClampX(), 2);
        gd.addNumericField("Clamp_y", fitConfig.getClampY(), 2);
        if (isAstigmatism) {
          gd.addNumericField("Clamp_z", fitConfig.getClampZ(), 2);
        } else {
          if (nParams > 1 || !fitConfig.isFixedPsf()) {
            gd.addNumericField("Clamp_sx", fitConfig.getClampXSd(), 2);
          }
          if (nParams > 1) {
            gd.addNumericField("Clamp_sy", fitConfig.getClampYSd(), 2);
          }
          if (nParams > 2) {
            gd.addNumericField("Clamp_angle", fitConfig.getClampAngle(), 2);
          }
        }
      }

      // Extra parameters are needed for calibrated fit solvers
      if (requireCalibration) {
        switch (calibration.getCameraType()) {
          case CCD:
          case EMCCD:
          case SCMOS:
            break;
          default:
            IJ.error(TITLE, fitSolverName + " requires camera calibration");
            return false;
        }

        gd.addMessage(fitSolverName + " requires calibration for camera: "
            + CalibrationProtosHelper.getName(calibration.getCameraType()));
        if (calibration.isScmos()) {
          final String[] models = CameraModelManager.listCameraModels(true);
          gd.addChoice("Camera_model_name", models, fitConfig.getCameraModelName());
        } else {
          gd.addNumericField("Camera_bias", calibration.getBias(), 2, 6, "Count");
          gd.addNumericField("Gain", calibration.getCountPerPhoton(), 2, 6, "Count/photon");
          gd.addNumericField("Read_noise", calibration.getReadNoise(), 2, 6, "Count");
        }
      }

      gd.showDialog();
      if (gd.wasCanceled()) {
        return false;
      }

      fitConfig.setRelativeThreshold(getThresholdNumber(gd));
      fitConfig.setAbsoluteThreshold(getThresholdNumber(gd));
      fitConfig.setParameterRelativeThreshold(getThresholdNumber(gd));
      fitConfig.setParameterAbsoluteThreshold(getThresholdNumber(gd));
      fitConfig.setMaxIterations((int) gd.getNextNumber());
      if (isLvm) {
        fitConfig.setLambda(gd.getNextNumber());
      }
      if (isFastMml) {
        fitConfig.setFixedIterations(gd.getNextBoolean());
        fitConfig.setLineSearchMethod(gd.getNextChoiceIndex());
      }

      fitConfig.setUseClamping(gd.getNextBoolean());
      fitConfig.setUseDynamicClamping(gd.getNextBoolean());
      if (extraOptions) {
        fitConfig.setClampBackground(Math.abs(gd.getNextNumber()));
        fitConfig.setClampSignal(Math.abs(gd.getNextNumber()));
        fitConfig.setClampX(Math.abs(gd.getNextNumber()));
        fitConfig.setClampY(Math.abs(gd.getNextNumber()));
        if (isAstigmatism) {
          fitConfig.setClampZ(Math.abs(gd.getNextNumber()));
        } else {
          if (nParams > 1 || !fitConfig.isFixedPsf()) {
            fitConfig.setClampXSd(Math.abs(gd.getNextNumber()));
          }
          if (nParams > 1) {
            fitConfig.setClampYSd(Math.abs(gd.getNextNumber()));
          }
          if (nParams > 2) {
            fitConfig.setClampAngle(Math.abs(gd.getNextNumber()));
          }
        }
      }

      if (requireCalibration) {
        if (calibration.isScmos()) {
          fitConfig.setCameraModelName(gd.getNextChoice());
        } else {
          calibration.setBias(Math.abs(gd.getNextNumber()));
          calibration.setCountPerPhoton(Math.abs(gd.getNextNumber()));
          calibration.setReadNoise(Math.abs(gd.getNextNumber()));
          fitConfig.setCalibration(calibration.getCalibration());
        }
      }

      // Do this even if collection of calibration settings was ignored. This ensures the
      // camera model is set.
      if (calibration.isScmos()) {
        fitConfig.setCameraModel(CameraModelManager.load(fitConfig.getCameraModelName()));
        if (!checkCameraModel(fitConfig, sourceBounds, bounds, true)) {
          return false;
        }
      }

      if (saveSettings) {
        saveFitEngineSettings(config);
      }

      try {
        if (isLvm) {
          ParameterUtils.isAboveZero("Lambda", fitConfig.getLambda());
        }
        // This call will check if the configuration is OK (including convergence criteria)
        fitConfig.getFunctionSolver();
      } catch (final IllegalArgumentException | IllegalStateException ex) {
        IJ.error(TITLE, ex.getMessage());
        return false;
      }
    } else {
      IJ.error(TITLE, "Unknown fit solver: " + fitSolver);
      return false;
    }

    if (config.isIncludeNeighbours() && !fitConfig.getFunctionSolver().isBounded()) {
      IJ.error(TITLE, "Including neighbours requires a bounded fit solver");
      return false;
    }

    return true;
  }

  /**
   * Gets the threshold number from the next string field. If invalid then -1 is returned (which
   * deactivates the threshold).
   *
   * @param gd the generic dialog
   * @return the number
   */
  private static double getThresholdNumber(ExtendedGenericDialog gd) {
    try {
      return Double.parseDouble(gd.getNextString());
    } catch (final NumberFormatException ex) {
      return -1;
    }
  }

  /**
   * Check the camera model covers the region of the source.
   *
   * @param fitConfig the fit config
   * @param sourceBounds the source bounds of the input image
   * @param cropBounds the crop bounds (relative to the input image)
   * @param initialise the initialise flag
   * @return true, if successful
   * @throws IllegalStateException if no camera model exists for the camera type
   */
  private static boolean checkCameraModel(FitConfiguration fitConfig, Rectangle sourceBounds,
      Rectangle cropBounds, boolean initialise) {
    final CalibrationReader calibration = fitConfig.getCalibrationReader();
    if (calibration.isScmos() && sourceBounds != null) {
      CameraModel cameraModel = fitConfig.getCameraModel();

      // The camera model origin must be reset to be relative to the source bounds origin
      cameraModel = cropCameraModel(cameraModel, sourceBounds, cropBounds, true);
      if (cameraModel == null) {
        return false;
      }

      if (initialise && cameraModel instanceof PerPixelCameraModel) {
        ((PerPixelCameraModel) cameraModel).initialise();
      }
      fitConfig.setCameraModel(cameraModel);
    }
    return true;
  }

  /**
   * Combine the source bounds with an internal crop bounds to create the region required from the
   * camera model.
   *
   * @param sourceBounds the source bounds
   * @param cropBounds the crop bounds
   * @return the rectangle
   */
  public static Rectangle combineBounds(Rectangle sourceBounds, Rectangle cropBounds) {
    if (sourceBounds == null) {
      throw new NullPointerException("No source bounds");
    }
    if (cropBounds == null) {
      return sourceBounds;
    }

    // Check the crop is inside the source
    if (!new Rectangle(sourceBounds.width, sourceBounds.height).contains(cropBounds)) {
      throw new IllegalArgumentException(
          "Crop bounds does not fit within the source width x height");
    }

    // Make relative
    if (sourceBounds.x != 0 || sourceBounds.y != 0) {
      cropBounds = (Rectangle) cropBounds.clone();
      cropBounds.x += sourceBounds.x;
      cropBounds.y += sourceBounds.y;
    }

    return cropBounds;
  }

  /**
   * Crop a camera model for processing data from a cropped image frame of the given bounds. The
   * target bounds are created by combining the crop with the source bounds. The camera model bounds
   * will be checked to verify that the target fits within the model.
   *
   * <p>If the model is smaller then an error is thrown.
   *
   * <p>If the model is larger then a crop will be made using a dialog to select the crop.
   *
   * <p>If the model is the same size then no crop is made, even if the origin is incorrect.
   *
   * <p>Optionally the model can be updated so that the origin is relative to the source bounds. If
   * no crop is used the origin will be 0,0. Otherwise it will be equal to the crop origin.
   *
   * <p>This method can be used to prepare a camera model for processing images frames of crop width
   * x height.
   *
   * @param cameraModel the camera model
   * @param sourceBounds the source bounds
   * @param cropBounds the crop bounds (relative to the input image). If null then the full width x
   *        height of the source is used.
   * @param resetOrigin the reset origin flag (to set the output camera model origin to match the
   *        source bounds)
   * @return the camera model (or null if the dialog was cancelled)
   * @throws IllegalArgumentException If the model is null or the crop cannot be done
   */
  public static CameraModel cropCameraModel(CameraModel cameraModel, Rectangle sourceBounds,
      Rectangle cropBounds, boolean resetOrigin) {
    if (cameraModel == null) {
      throw new IllegalArgumentException("No camera model");
    }
    Rectangle modelBounds = cameraModel.getBounds();
    if (modelBounds == null || sourceBounds == null) {
      return cameraModel;
    }

    // Combine the source bounds with the crop
    sourceBounds = combineBounds(sourceBounds, cropBounds);

    final int width = sourceBounds.width;
    final int height = sourceBounds.height;
    if (modelBounds.width < width || modelBounds.height < height) {
      throw new IllegalArgumentException(String.format(
          "Camera model bounds [x=%d,y=%d,width=%d,height=%d] is smaller than image size [%dx%d]",
          modelBounds.x, modelBounds.y, modelBounds.width, modelBounds.height, width, height));
    } else if (modelBounds.width > width || modelBounds.height > height) {
      final GenericDialog gd2 = new GenericDialog("Crop Camera Model");
      //@formatter:off
      ImageJUtils.addMessage(gd2,
          "WARNING:\n \nCamera model bounds\n[x=%d,y=%d,width=%d,height=%d]\n"
              + "are larger than the image size [%dx%d].\n \nCrop the model?",
          modelBounds.x, modelBounds.y, modelBounds.width, modelBounds.height,
          width, height
          );
      //@formatter:on
      final int upperx = modelBounds.x + modelBounds.width - width;
      final int uppery = modelBounds.y + modelBounds.height - height;
      int ox = sourceBounds.x;
      int oy = sourceBounds.y;
      gd2.addSlider("Origin_x", modelBounds.x, upperx, MathUtils.clip(modelBounds.x, upperx, ox));
      gd2.addSlider("Origin_y", modelBounds.y, uppery, MathUtils.clip(modelBounds.y, uppery, oy));
      gd2.showDialog();
      if (gd2.wasCanceled()) {
        return null;
      }
      ox = (int) gd2.getNextNumber();
      oy = (int) gd2.getNextNumber();

      final Rectangle bounds = new Rectangle(ox, oy, width, height);
      cameraModel = cameraModel.crop(bounds, false);
      modelBounds = cameraModel.getBounds();
      if (modelBounds.width != bounds.width || modelBounds.height != bounds.height) {
        throw new IllegalArgumentException("Failed to crop camera model using bounds: " + bounds);
      }
    }

    if (resetOrigin) {
      // Reset origin to the source origin
      cameraModel = cameraModel.copy();
      if (cropBounds == null) {
        cameraModel.setOrigin(0, 0);
      } else {
        cameraModel.setOrigin(cropBounds.x, cropBounds.y);
      }
    }

    return cameraModel;
  }

  /**
   * Add a result output.
   *
   * <p>This can be called after {@link #initialiseImage(ImageSource, Rectangle, boolean)} and
   * before {@link #initialiseFitting()} to add to the configured result outputs.
   *
   * @param peakResults the peak results
   */
  public void addPeakResults(PeakResults peakResults) {
    if (results != null) {
      results.addOutput(peakResults);
    }
  }

  private void addTableResults(PeakResultsList resultsList) {
    final ImageJTablePeakResults peakResults =
        ResultsManager.addTableResults(resultsList, resultsSettings.getResultsTableSettings(),
            resultsSettings.getShowDeviations(), false, false, false);
    if (peakResults != null) {
      peakResults.setShowZ(PsfHelper.is3D(resultsList.getPsf()));
      peakResults.setClearAtStart(simpleFit);
      peakResults.setShowEndFrame(getShowEndFrame());
    }
  }

  private boolean getShowEndFrame() {
    return extraSettings != null && extraSettings.integrateFrames > 1;
  }

  private void addFileResults(PeakResultsList resultsList) {
    final ResultsFileSettings resultsFileSettings = this.resultsSettings.getResultsFileSettings();
    if (resultsFileSettings.getFileFormat().getNumber() > 0) {
      String resultsFilename = null;
      if (resultsFileSettings.getResultsDirectory() != null
          && new File(resultsFileSettings.getResultsDirectory()).exists()) {
        resultsFilename = resultsFileSettings.getResultsDirectory() + File.separatorChar
            + source.getName() + ".results."
            + ResultsProtosHelper.getExtension(resultsFileSettings.getFileFormat());

        // This is used for running via other code calling PeakFit methods,
        // i.e. not as an ImageJ plugin.
      } else if (pluginFlags == 0) {
        resultsFilename = resultsFileSettings.getResultsFilename();
      }
      final PeakResults r = ResultsManager.addFileResults(resultsList, resultsFileSettings,
          resultsFilename, this.resultsSettings.getShowDeviations(), getShowEndFrame(), false);
      if (r instanceof FilePeakResults) {
        final FilePeakResults fr = (FilePeakResults) r;
        fr.setSortAfterEnd(Prefs.getThreads() > 1);
      }
    }
  }

  private void addMemoryResults(PeakResultsList resultsList, boolean force) {
    if (resultsSettings.getResultsInMemorySettings().getInMemory() || force) {
      final MemoryPeakResults peakResults = new MemoryPeakResults();
      peakResults.setSortAfterEnd(Prefs.getThreads() > 1);
      resultsList.addOutput(peakResults);
      MemoryPeakResults.addResults(peakResults);
    }
  }

  private void addDefaultResults(PeakResultsList resultsList) {
    if (resultsList.numberOfOutputs() == 0) {
      if (logger != null) {
        logger.info("No results output configured. Defaulting to memory");
      }
      addMemoryResults(resultsList, true);
    }
  }

  @Override
  public void run(ImageProcessor ip) {
    if (source == null) {
      IJ.error(TITLE, "No valid image source configured");
      return;
    }
    if (fitMaxima) {
      runMaximaFitting();
    } else {
      // All setup is already done so simply run
      run();
    }

    addSingleFrameOverlay();
  }

  /**
   * Process the image. The current ROI will be used to define the region processed. The noise can
   * be estimated using the entire frame or the ROI region.
   *
   * @param imp the imp
   * @param ignoreBoundsForNoise If true estimate the noise from the entire frame, otherwise use
   *        only the ROI bounds
   */
  public void run(ImagePlus imp, boolean ignoreBoundsForNoise) {
    run(new IJImageSource(imp), getBounds(imp), ignoreBoundsForNoise);
  }

  /**
   * Process the image.
   *
   * @param imageSource the image source
   * @param bounds the bounds
   * @param ignoreBoundsForNoise the ignore bounds for noise
   */
  public void run(ImageSource imageSource, Rectangle bounds, boolean ignoreBoundsForNoise) {
    if (initialise(imageSource, bounds, ignoreBoundsForNoise)) {
      run();
    }
  }

  /**
   * Locate the peaks in the configured image source. Results are saved to the configured output.
   *
   * <p>This must be called after initialisation with an image source. Note that each call to this
   * method must be preceded with initialisation to prepare the image and output options.
   */
  @SuppressWarnings("null")
  public void run() {
    if (source == null) {
      return;
    }

    final int totalFrames = source.getFrames();

    final ImageStack stack =
        (extraSettings.showProcessedFrames) ? new ImageStack(bounds.width, bounds.height) : null;

    // Do not crop the region from the source if the bounds match the source dimensions
    final Rectangle cropBounds =
        (bounds.x == 0 && bounds.y == 0 && bounds.width == source.getWidth()
            && bounds.height == source.getHeight()) ? null : bounds;

    // Use the FitEngine to allow multi-threading.
    final FitEngine engine = createFitEngine(getNumberOfThreads(totalFrames));
    if (engine == null) {
      return;
    }

    final int step = ImageJUtils.getProgressInterval(totalFrames);

    // To pre-process data for noise estimation
    boolean isFitCameraCounts = false;
    CameraModel cameraModel = null;
    if (ignoreBoundsForNoise) {
      isFitCameraCounts = fitConfig.isFitCameraCounts();
      cameraModel = fitConfig.getCameraModel();
    }

    runTime = System.nanoTime();
    boolean shutdown = false;
    int slice = 0;
    final String format = String.format("Slice: %%d / %d (Results=%%d)", totalFrames);
    while (!shutdown) {
      // Noise can optionally be estimated from the entire frame
      float[] data = (ignoreBoundsForNoise) ? source.next() : source.next(cropBounds);
      if (data == null) {
        break;
      }

      if (++slice % step == 0) {
        final int frames = slice;
        if (ImageJUtils.showStatus(() -> String.format(format, frames, results.size()))) {
          IJ.showProgress(slice, totalFrames);
        }
      }

      float noise = Float.NaN;
      if (ignoreBoundsForNoise) {
        // We must pre-process the data before noise estimation
        final float[] data2 = data.clone();
        if (isFitCameraCounts) {
          cameraModel.removeBias(data2);
        } else {
          cameraModel.removeBiasAndGain(data2);
        }

        noise = FitWorker.estimateNoise(data2, source.getWidth(), source.getHeight(),
            config.getNoiseMethod());

        // Crop the data to the region
        data =
            ImageJImageConverter.getData(data, source.getWidth(), source.getHeight(), bounds, null);
      }

      if (stack != null) {
        stack.addSlice(String.format("Frame %d - %d", source.getStartFrameNumber(),
            source.getEndFrameNumber()), data);
      }

      // Get the frame number from the source to allow for interlaced and aggregated data
      engine.run(
          createJob(source.getStartFrameNumber(), source.getEndFrameNumber(), data, bounds, noise));

      shutdown = escapePressed();
    }

    engine.end(shutdown);
    time = engine.getTime();
    runTime = System.nanoTime() - runTime;

    if (stack != null) {
      ImageJUtils.display("Processed frames", stack);
    }

    showResults();

    source.close();
  }

  private void addSingleFrameOverlay() {
    // If a single frame was processed add the peaks as an overlay if they are in memory
    ImagePlus localImp = this.imp;

    if (fitMaxima && singleFrame > 0 && source instanceof IJImageSource) {
      final String title = source.getName();
      localImp = WindowManager.getImage(title);
    }

    if (singleFrame > 0 && localImp != null) {
      MemoryPeakResults memoryResults = null;
      for (final PeakResults r : this.results.toArray()) {
        if (r instanceof MemoryPeakResults) {
          memoryResults = (MemoryPeakResults) r;
          break;
        }
      }
      if (memoryResults == null || memoryResults.size() == 0) {
        return;
      }

      final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
      gd.enableYesNoCancel();
      gd.hideCancelButton();
      gd.addMessage("Add the fitted localisations as an overlay?");
      gd.showDialog();
      if (!gd.wasOKed()) {
        return;
      }

      final LUT lut = LutHelper.createLut(LutColour.ICE);
      final Overlay o = new Overlay();
      final int size = memoryResults.size();
      final Counter j = new Counter(size);
      final ImagePlus finalImp = localImp;
      memoryResults.forEach(DistanceUnit.PIXEL, (XyResultProcedure) (x, y) -> {
        final PointRoi roi = new OffsetPointRoi(x, y);
        final Color c = LutHelper.getColour(lut, j.decrementAndGet(), size);
        roi.setStrokeColor(c);
        roi.setFillColor(c);
        if (finalImp.getStackSize() > 1) {
          roi.setPosition(singleFrame);
        }
        o.add(roi);
      });
      localImp.setOverlay(o);
      localImp.getWindow().toFront();
    }
  }

  /**
   * Gets the number of threads.
   *
   * @param totalFrames the total frames
   * @return the number of threads
   */
  private int getNumberOfThreads(int totalFrames) {
    final int t = Math.max(1, (int) (settings.fractionOfThreads * Prefs.getThreads()));
    return Math.min(totalFrames, t);
  }

  /**
   * Check if the frame should be ignored (relevant when using interlaced data).
   *
   * @param frame the frame
   * @return True if the frame should be ignored
   */
  private boolean ignoreFrame(int frame) {
    if (!extraSettings.interlacedData) {
      return false;
    }

    // Check if the frame is allowed:
    // Start
    // |
    // |----|Block|Skip|Block|Skip|Block|Skip
    // Note the source data is still read so that the source is incremented.
    if (frame < extraSettings.dataStart) {
      return true;
    }
    final int frameInBlock =
        (frame - extraSettings.dataStart) % (extraSettings.dataBlock + extraSettings.dataSkip);
    return frameInBlock >= extraSettings.dataBlock;
  }

  private FitJob createJob(int startFrame, int endFrame, float[] data, Rectangle bounds,
      float noise) {
    FitParameters fitParams = null;
    if (startFrame != endFrame) {
      fitParams = new FitParameters();
      fitParams.endT = endFrame;
    }

    if (maximaIdentification) {
      if (fitParams == null) {
        fitParams = new FitParameters();
      }
      fitParams.fitTask = FitTask.MAXIMA_IDENITIFICATION;
      fitParams.noise = noise;
    } else if (!Float.isNaN(noise)) {
      if (fitParams == null) {
        fitParams = new FitParameters();
      }
      fitParams.fitTask = FitTask.PSF_FITTING;
      fitParams.noise = noise;
    }

    if (fitParams != null) {
      return new ParameterisedFitJob(fitParams, startFrame, data, bounds);
    }
    return new FitJob(startFrame, data, bounds);
  }

  /**
   * Creates the fit job to fit all the given candidate maxima.
   *
   * @param sliceCandidates the slice candidates
   * @param startFrame the start frame
   * @param endFrame the end frame
   * @param data the data
   * @param bounds the bounds
   * @param noise the noise
   * @return the fit job
   */
  private static FitJob createMaximaFitJob(int[] maxIndices, int startFrame, int endFrame,
      float[] data, Rectangle bounds, float noise) {
    final FitParameters fitParams = new FitParameters();
    fitParams.maxIndices = maxIndices;
    if (startFrame != endFrame) {
      fitParams.endT = endFrame;
    }
    fitParams.fitTask = FitTask.PSF_FITTING;
    fitParams.noise = noise;
    return new ParameterisedFitJob(fitParams, startFrame, data, bounds);
  }

  private static boolean escapePressed() {
    if (IJ.escapePressed()) {
      IJ.log(TITLE + " stopping ...");
      IJ.beep();
      return true;
    }
    return false;
  }

  /**
   * Creates a fitting engine using the current configuration.
   *
   * @return The fitting engine
   */
  public FitEngine createFitEngine() {
    return createFitEngine(Prefs.getThreads());
  }

  /**
   * Creates a fitting engine using the current configuration.
   *
   * @param numberOfThreads the number of threads
   * @return The fitting engine
   */
  public FitEngine createFitEngine(int numberOfThreads) {
    // Use a blocking queue to enable progress tracking on the IJ progress bar.
    // Use a large queue size to allow images read from disk to be pre-loaded.
    return createFitEngine(numberOfThreads, FitQueue.BLOCKING, numberOfThreads * 10);
  }

  /**
   * Creates a fitting engine using the current configuration.
   *
   * @param numberOfThreads the number of threads
   * @param queue the queue
   * @param queueSize the queue size
   * @return The fiting engine
   */
  public FitEngine createFitEngine(int numberOfThreads, FitQueue queue, int queueSize) {
    // Ensure thread safety
    final PeakResultsList list = (numberOfThreads > 1) ? results.getThreadSafeList() : results;

    // Reduce to single object for speed
    final PeakResults r = (results.numberOfOutputs() == 1) ? list.toArray()[0] : list;

    // Update the configuration
    if (!updateFitConfiguration(config)) {
      return null;
    }

    final FitEngine engine = FitEngine.create(config, r, numberOfThreads, queue, queueSize);

    // Write settings out to the IJ log
    if (resultsSettings.getLogProgress()) {
      IJ.log(LOG_SPACER);
      IJ.log("Peak Fit");
      IJ.log(LOG_SPACER);
      ImageJUtils.log("Initial Peak SD = %s,%s", MathUtils.rounded(fitConfig.getInitialXSd()),
          MathUtils.rounded(fitConfig.getInitialYSd()));
      final SpotFilter spotFilter = engine.getSpotFilter();
      IJ.log("Spot Filter = " + spotFilter.getDescription());
      final int w = 2 * engine.getFitting() + 1;
      ImageJUtils.log("Fit window = %d x %d", w, w);
      if (!fitConfig.isDisableSimpleFilter()) {
        IJ.log("Coordinate shift = "
            + MathUtils.rounded(config.getFitConfiguration().getCoordinateShift()));
        IJ.log("Signal strength = " + MathUtils.rounded(fitConfig.getSignalStrength()));
      }
      if (fitConfig.isDirectFilter()) {
        IJ.log("Smart filter = " + fitConfig.getSmartFilter().getDescription());
      }
      if (extraOptions) {
        IJ.log("Noise = " + MathUtils.rounded(fitConfig.getNoise()));
      }
      IJ.log("Width factor = " + MathUtils.rounded(fitConfig.getMaxWidthFactor()));
      IJ.log(LOG_SPACER);
    }

    return engine;
  }

  /**
   * Updates the configuration for peak fitting. Configures the calculation of residuals, logging
   * and peak validation.
   *
   * @param config the config
   * @return true, if successful
   */
  private boolean updateFitConfiguration(FitEngineConfiguration config) {
    final FitConfiguration localFitConfig = config.getFitConfiguration();

    // Adjust the settings that are relevant within the fitting configuration.
    localFitConfig.setComputeResiduals(config.getResidualsThreshold() < 1);
    logger =
        (resultsSettings.getLogProgress()) ? ImageJPluginLoggerHelper.getLogger(getClass()) : null;
    localFitConfig.setLog(logger);

    if (resultsSettings.getShowDeviations()) {
      // Note: This may already by true if the deviations are needed for the smart filter
      localFitConfig.setComputeDeviations(resultsSettings.getShowDeviations());
    }

    config.configureOutputUnits();

    return checkCameraModel(localFitConfig, source.getBounds(), bounds, true);
  }

  /**
   * Load the selected results from memory. All multiple frame results are added directly to the
   * results. All single frame results are added to a list of candidate maxima per frame and fitted
   * using the configured parameters.
   */
  private void runMaximaFitting() {
    final MemoryPeakResults memoryResults =
        ResultsManager.loadInputResults(settings.inputOption, false, DistanceUnit.PIXEL);
    if (memoryResults == null || memoryResults.size() == 0) {
      log("No results for maxima fitting");
      return;
    }

    // The total frames (for progress reporting)
    int totalFrames;
    // A function that can convert a frame into a set of candidate indices
    final IntFunction<int[]> frameToMaxIndices;
    // The frames to process (should be sorted ascending)
    Supplier<IntStream> frames;

    // Support fitting all time frames with the same results.
    if (settings.fitAcrossAllFrames) {
      // All results are candidates
      // Check if the input spans multiple frames
      if (getSingleFrame(memoryResults) == 0) {
        final int min = memoryResults.getMinFrame();
        final int max = memoryResults.getMaxFrame();
        final GenericDialog gd = new GenericDialog(TITLE);
        gd.enableYesNoCancel();
        gd.hideCancelButton();
        ImageJUtils.addMessage(gd,
            "Candidate maxima for fitting span multiple frames (%d-%d).\n \n"
                + "Please confirm the %s are correct.",
            min, max, TextUtils.pleural(memoryResults.size(), "candidate"));
        gd.showDialog();
        if (!gd.wasOKed()) {
          return;
        }
      }

      final int[] maxIndices = getMaxIndices(Arrays.asList(memoryResults.toArray()));

      // This may not work correctly if using for example a series image source that
      // incorrectly estimates the number of frames
      totalFrames = source.getFrames();
      frameToMaxIndices = frame -> maxIndices;
      frames = () -> IntStream.rangeClosed(1, totalFrames);
    } else {

      // Build a map between the time-frame and the results in that frame.
      final Map<Integer, List<PeakResult>> map = Arrays.stream(memoryResults.toArray()).parallel()
          // Ensure only single frame results are used to pick candidates
          .filter(peakResult -> peakResult.getFrame() == peakResult.getEndFrame())
          .collect(Collectors.groupingBy(PeakResult::getFrame));

      totalFrames = map.size();

      // Build a function that can convert a frame into a set of candidate indices
      frameToMaxIndices = frame -> getMaxIndices(map.get(frame));

      frames = () -> map.keySet().stream().mapToInt(Integer::intValue).sorted();
    }

    final ImageStack stack =
        (extraSettings.showProcessedFrames) ? new ImageStack(bounds.width, bounds.height) : null;

    // Use the FitEngine to allow multi-threading.
    final FitEngine engine = createFitEngine(getNumberOfThreads(totalFrames));
    if (engine == null) {
      return;
    }

    final int step = ImageJUtils.getProgressInterval(totalFrames);

    // No crop bounds are supported.
    // To pre-process data for noise estimation
    final boolean isFitCameraCounts = fitConfig.isFitCameraCounts();
    final CameraModel cameraModel = fitConfig.getCameraModel();

    runTime = System.nanoTime();
    final AtomicBoolean shutdown = new AtomicBoolean();
    final String format = String.format("Slice: %%d / %d (Results=%%d)", totalFrames);
    frames.get().forEachOrdered(slice -> {
      if (shutdown.get() || escapePressed()) {
        shutdown.set(true);
        return;
      }
      final float[] data = source.get(slice);
      if (data == null) {
        shutdown.set(true);
        return;
      }

      if (slice % step == 0) {
        if (ImageJUtils.showStatus(() -> String.format(format, slice, results.size()))) {
          IJ.showProgress(slice, totalFrames);
        }
      }

      // We must pre-process the data before noise estimation
      final float[] data2 = data.clone();
      if (isFitCameraCounts) {
        cameraModel.removeBias(data2);
      } else {
        cameraModel.removeBiasAndGain(data2);
      }

      final float noise = FitWorker.estimateNoise(data2, source.getWidth(), source.getHeight(),
          config.getNoiseMethod());

      if (stack != null) {
        stack.addSlice(String.format("Frame %d - %d", source.getStartFrameNumber(),
            source.getEndFrameNumber()), data);
      }

      // Get the frame number from the source to allow for interlaced and aggregated data
      engine.run(createMaximaFitJob(frameToMaxIndices.apply(slice), source.getStartFrameNumber(),
          source.getEndFrameNumber(), data, bounds, noise));
    });

    engine.end(shutdown.get());
    time = engine.getTime();
    runTime = System.nanoTime() - runTime;

    if (stack != null) {
      ImageJUtils.display("Processed frames", stack);
    }

    showResults();

    source.close();
  }

  private int[] getMaxIndices(List<PeakResult> sliceCandidates) {
    final int[] maxIndices = new int[sliceCandidates.size()];
    int count = 0;
    for (final PeakResult result : sliceCandidates) {
      maxIndices[count++] = result.getOrigX() + bounds.width * result.getOrigY();
    }
    return maxIndices;
  }

  /**
   * Gets the total fitting time in nanoseconds.
   *
   * @return The total fitting time in nanoseconds.
   */
  public long getTime() {
    return time;
  }

  /**
   * Gets the total number of localisations.
   *
   * @return The total number of localisations.
   */
  public int getSize() {
    // The list returns the size using the first entry.
    return results.size();
  }

  /**
   * Checks if is silent. If true, do not output any log messages (e.g. the final time and count).
   *
   * @return true, if is silent
   */
  public boolean isSilent() {
    return silent;
  }

  /**
   * Sets the silent option. If true, do not output any log messages (e.g. the final time and
   * count).
   *
   * @param silent the new silent
   */
  public void setSilent(boolean silent) {
    this.silent = silent;
  }
}
