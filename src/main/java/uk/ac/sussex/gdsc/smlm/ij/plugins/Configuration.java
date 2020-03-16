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
import ij.plugin.PlugIn;
import java.awt.Checkbox;
import java.awt.Choice;
import java.awt.Color;
import java.awt.SystemColor;
import java.awt.TextField;
import java.awt.event.ItemEvent;
import java.io.File;
import java.util.Iterator;
import java.util.Vector;
import java.util.concurrent.atomic.AtomicReference;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.utils.FileUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FitEngineSettings;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.data.config.PsfProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.TemplateProtos.TemplateSettings;
import uk.ac.sussex.gdsc.smlm.engine.FitConfiguration;
import uk.ac.sussex.gdsc.smlm.engine.FitEngineConfiguration;
import uk.ac.sussex.gdsc.smlm.ij.plugins.PeakFit.FitEngineConfigurationProvider;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;

/**
 * Adjust the configuration used for fitting.
 */
public class Configuration implements PlugIn {
  private static final String TITLE = "Fit Configuration";

  private FitEngineConfiguration config;
  private FitConfiguration fitConfig;

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
  private TextField textFailuresLimit;
  private TextField textPassRate;
  private Checkbox textIncludeNeighbours;
  private TextField textNeighbourHeightThreshold;
  private TextField textResidualsThreshold;
  private TextField textDuplicateDistance;
  private Checkbox textSmartFilter;
  private Checkbox textDisableSimpleFilter;
  private TextField textCoordinateShiftFactor;
  private TextField textSignalStrength;
  private TextField textMinPhotons;
  private TextField textPrecisionThreshold;
  private TextField textMinWidthFactor;
  private TextField textWidthFactor;

  /** The plugin settings. */
  private Settings pluginSettings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    String templateFilename = "";
    String notes = "";

    Settings() {
      // Set defaults
      templateFilename = "";
      notes = "";
    }

    Settings(Settings source) {
      templateFilename = source.templateFilename;
      notes = source.notes;
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
     * Save the settings. This can be called only once as it saves via a reference.
     */
    void save() {
      lastSettings.set(this);
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    showDialog(true);
  }

  /**
   * Show the current properties.
   *
   * @param save the save
   * @return true, if successful
   */
  public boolean showDialog(boolean save) {
    config = SettingsManager.readFitEngineConfiguration(0);
    return showDialog(config, save);
  }

  /**
   * Show the current properties.
   *
   * @param fitEngineConfiguration the fit engine configuration
   * @param save the save
   * @return true, if successful
   */
  public boolean showDialog(FitEngineConfiguration fitEngineConfiguration, boolean save) {
    this.config = fitEngineConfiguration;
    fitConfig = config.getFitConfiguration();
    pluginSettings = Settings.load();
    pluginSettings.save();

    final CalibrationReader calibrationReader = fitConfig.getCalibrationReader();

    ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("fit-configuration"));
    gd.addMessage("Configuration settings for the single-molecule localisation microscopy plugins");

    final String[] templates = ConfigurationTemplate.getTemplateNames(true);
    gd.addChoice("Template", templates, templates[0]);

    PeakFit.addCameraOptions(gd, fitConfig);
    gd.addNumericField("Calibration (nm/px)", calibrationReader.getNmPerPixel(), 2);
    gd.addNumericField("Exposure_time (ms)", calibrationReader.getExposureTime(), 2);

    gd.addMessage("--- Gaussian parameters ---");
    PeakFit.addPsfOptions(gd, fitConfig);

    gd.addMessage("--- Maxima identification ---");
    final FitEngineConfigurationProvider provider = this::getFitEngineConfiguration;
    PeakFit.addDataFilterOptions(gd, provider);
    PeakFit.addSearchOptions(gd, provider);
    PeakFit.addBorderOptions(gd, provider);
    PeakFit.addFittingOptions(gd, provider);

    gd.addMessage("--- Gaussian fitting ---");

    gd.addChoice("Fit_solver", SettingsManager.getFitSolverNames(),
        fitConfig.getFitSolver().ordinal());

    // Parameters specific to each Fit solver are collected in a second dialog

    gd.addNumericField("Fail_limit", config.getFailuresLimit(), 0);
    gd.addNumericField("Pass_rate", config.getPassRate(), 2);
    gd.addCheckbox("Include_neighbours", config.isIncludeNeighbours());
    gd.addSlider("Neighbour_height", 0.01, 1, config.getNeighbourHeightThreshold());
    gd.addSlider("Residuals_threshold", 0.01, 1, config.getResidualsThreshold());

    PeakFit.addDuplicateDistanceOptions(gd, provider);

    gd.addMessage(
        "--- Peak filtering ---\nDiscard fits that shift; are too low; or expand/contract");

    gd.addCheckbox("Smart_filter", fitConfig.isSmartFilter());
    gd.addCheckbox("Disable_simple_filter", fitConfig.isDisableSimpleFilter());
    gd.addSlider("Shift_factor", 0.01, 2, fitConfig.getCoordinateShiftFactor());
    gd.addNumericField("Signal_strength", fitConfig.getSignalStrength(), 2);
    gd.addNumericField("Min_photons", fitConfig.getMinPhotons(), 0);
    gd.addSlider("Min_width_factor", 0, 0.99, fitConfig.getMinWidthFactor());
    gd.addSlider("Width_factor", 1, 4.5, fitConfig.getMaxWidthFactor());
    PeakFit.addPrecisionOptions(gd, this::getFitConfiguration);

    // Add a mouse listener to the config file field
    if (ImageJUtils.isShowGenericDialog()) {
      final Vector<TextField> numerics = gd.getNumericFields();
      final Vector<Checkbox> checkboxes = gd.getCheckboxes();
      final Vector<Choice> choices = gd.getChoices();

      final Iterator<TextField> nu = numerics.iterator();
      final Iterator<Checkbox> cb = checkboxes.iterator();
      final Iterator<Choice> ch = choices.iterator();

      final Choice textTemplate = ch.next();
      textTemplate.addItemListener(this::itemStateChanged);

      textCameraType = ch.next();
      textNmPerPixel = nu.next();
      textExposure = nu.next();
      textPsf = ch.next();
      textDataFilterType = ch.next();
      textDataFilterMethod = ch.next();
      textSmooth = nu.next();
      textSearch = nu.next();
      textBorder = nu.next();
      textFitting = nu.next();
      textFitSolver = ch.next();
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
      textMinWidthFactor = nu.next();
      textWidthFactor = nu.next();
      textPrecisionThreshold = nu.next();

      updateFilterInput();
      textSmartFilter.addItemListener(this::itemStateChanged);
      textDisableSimpleFilter.addItemListener(this::itemStateChanged);
    }

    if (save) {
      gd.enableYesNoCancel("Save", "Save template");
    }

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    // In case a template updated the calibration
    final CalibrationWriter calibrationWriter = fitConfig.getCalibrationWriter();

    // Ignore the template
    gd.getNextChoice();

    calibrationWriter.setCameraType(SettingsManager.getCameraTypeValues()[gd.getNextChoiceIndex()]);
    calibrationWriter.setNmPerPixel(gd.getNextNumber());
    calibrationWriter.setExposureTime(gd.getNextNumber());
    fitConfig.setCalibration(calibrationWriter.getCalibration());
    fitConfig.setPsfType(PeakFit.getPsfTypeValues()[gd.getNextChoiceIndex()]);
    config.setDataFilterType(gd.getNextChoiceIndex());
    config.setDataFilter(gd.getNextChoiceIndex(), Math.abs(gd.getNextNumber()), false, 0);
    config.setSearch(gd.getNextNumber());
    config.setBorder(gd.getNextNumber());
    config.setFitting(gd.getNextNumber());

    fitConfig.setFitSolver(gd.getNextChoiceIndex());

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
    fitConfig.setMinWidthFactor(gd.getNextNumber());
    fitConfig.setMaxWidthFactor(gd.getNextNumber());
    fitConfig.setPrecisionThreshold(gd.getNextNumber());

    gd.collectOptions();

    // Check arguments
    try {
      ParameterUtils.isAboveZero("nm per pixel", calibrationWriter.getNmPerPixel());
      ParameterUtils.isAboveZero("Exposure time", calibrationWriter.getExposureTime());
      if (fitConfig.getPsfTypeValue() != PSFType.ASTIGMATIC_GAUSSIAN_2D_VALUE) {
        ParameterUtils.isAboveZero("Initial SD0", fitConfig.getInitialXSd());
        if (fitConfig.getPsf().getParametersCount() > 1) {
          ParameterUtils.isAboveZero("Initial SD1", fitConfig.getInitialYSd());
        }
      }
      ParameterUtils.isAboveZero("Search_width", config.getSearch());
      ParameterUtils.isAboveZero("Fitting_width", config.getFitting());
      ParameterUtils.isPositive("Neighbour height threshold", config.getNeighbourHeightThreshold());
      ParameterUtils.isPositive("Residuals threshold", config.getResidualsThreshold());
      ParameterUtils.isPositive("Duplicate distance", config.getDuplicateDistance());

      if (!fitConfig.isSmartFilter()) {
        ParameterUtils.isPositive("Coordinate Shift factor", fitConfig.getCoordinateShiftFactor());
        ParameterUtils.isPositive("Signal strength", fitConfig.getSignalStrength());
        ParameterUtils.isPositive("Min photons", fitConfig.getMinPhotons());
        ParameterUtils.isPositive("Min width factor", fitConfig.getMinWidthFactor());
        ParameterUtils.isPositive("Width factor", fitConfig.getMaxWidthFactor());
        ParameterUtils.isPositive("Precision threshold", fitConfig.getPrecisionThreshold());
      }
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    if (gd.invalidNumber()) {
      return false;
    }

    final int flags = PeakFit.FLAG_NO_SAVE;
    if (!PeakFit.configurePsfModel(config, flags)) {
      return false;
    }
    if (!PeakFit.configureResultsFilter(config, flags)) {
      return false;
    }
    if (!PeakFit.configureDataFilter(config, flags)) {
      return false;
    }
    PeakFit.configureFitSolver(config, null, null, flags);

    if (save) {
      final boolean saveToFile = !gd.wasOKed();
      if (saveToFile) {
        gd = new ExtendedGenericDialog(TITLE);
        gd.addFilenameField("Template_filename", pluginSettings.templateFilename);
        gd.addMessage("Add notes to the template ...");
        gd.addTextAreas(pluginSettings.notes, null, 10, 60);
        gd.showDialog();
        if (gd.wasCanceled()) {
          return false;
        }
        pluginSettings.save();
        final String filename = gd.getNextString();
        pluginSettings.notes = gd.getNextText();

        if (filename != null) {
          pluginSettings.templateFilename = FileUtils.replaceExtension(filename, ".txt");
          final File file = new File(pluginSettings.templateFilename);
          final String name = FileUtils.removeExtension(file.getName());
          final TemplateSettings.Builder settings = TemplateSettings.newBuilder();
          settings.addNotes(pluginSettings.notes);
          settings.setCalibration(fitConfig.getCalibration());
          settings.setFitEngineSettings(config.getFitEngineSettings());
          // Note: No results settings are currently supported
          settings.setPsf(fitConfig.getPsf());
          if (!ConfigurationTemplate.saveTemplate(name, settings.build(), file)) {
            IJ.error(TITLE, "Failed to save to file: " + pluginSettings.templateFilename);
          }
        }
      } else {
        SettingsManager.writeSettings(config, 0);
      }
    }

    return true;
  }

  /**
   * Gets the fit engine configuration.
   *
   * @return the fit engine configuration
   */
  FitEngineConfiguration getFitEngineConfiguration() {
    return config;
  }

  private FitConfiguration getFitConfiguration() {
    return config.getFitConfiguration();
  }

  private void itemStateChanged(ItemEvent event) {
    if (event.getSource() instanceof Choice) {
      // Update the settings from the template
      final Choice choice = (Choice) event.getSource();
      final String templateName = choice.getSelectedItem();

      // Get the configuration template
      final TemplateSettings template = ConfigurationTemplate.getTemplate(templateName);

      if (template != null) {
        IJ.log("Applying template: " + templateName);
        final File file = ConfigurationTemplate.getTemplateFile(templateName);
        if (file != null) {
          pluginSettings.templateFilename = file.getPath();
        }

        final StringBuilder sb = new StringBuilder();
        for (final String note : template.getNotesList()) {
          sb.append(note);
          if (!note.endsWith("\n")) {
            sb.append("\n");
          }
          IJ.log(note);
        }
        pluginSettings.notes = sb.toString();

        final boolean custom = ConfigurationTemplate.isCustomTemplate(templateName);
        if (template.hasCalibration()) {
          refreshSettings(template.getCalibration());
        }
        if (template.hasPsf()) {
          refreshSettings(template.getPsf(), custom);
        }
        if (template.hasFitEngineSettings()) {
          refreshSettings(template.getFitEngineSettings());
        }
      }
    } else if (event.getSource() instanceof Checkbox) {
      if (event.getSource() == textSmartFilter) {
        textDisableSimpleFilter.setState(textSmartFilter.getState());
      }
      updateFilterInput();
    }
  }

  private void updateFilterInput() {
    if (textDisableSimpleFilter.getState()) {
      disableEditing(textCoordinateShiftFactor);
      disableEditing(textSignalStrength);
      disableEditing(textMinPhotons);
      // These are used to set bounds:
      // - textMinWidthFactor
      // - textWidthFactor
      disableEditing(textPrecisionThreshold);
    } else {
      enableEditing(textCoordinateShiftFactor);
      enableEditing(textSignalStrength);
      enableEditing(textMinPhotons);
      enableEditing(textPrecisionThreshold);
    }
  }

  private static void disableEditing(TextField textField) {
    textField.setEditable(false);
    textField.setBackground(SystemColor.control);
  }

  private static void enableEditing(TextField textField) {
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
   */
  private void refreshSettings(FitEngineSettings fitEngineSettings) {
    // Set the configuration
    // This will clear everything and merge the configuration
    this.config.setFitEngineSettings(fitEngineSettings);
    fitConfig = this.config.getFitConfiguration();

    textDataFilterType
        .select(SettingsManager.getDataFilterTypeNames()[config.getDataFilterType().ordinal()]);
    textDataFilterMethod.select(
        SettingsManager.getDataFilterMethodNames()[config.getDataFilterMethod(0).ordinal()]);
    textSmooth.setText("" + config.getDataFilterParameterValue(0));
    textSearch.setText("" + config.getSearch());
    textBorder.setText("" + config.getBorder());
    textFitting.setText("" + config.getFitting());
    textFitSolver.select(SettingsManager.getFitSolverNames()[fitConfig.getFitSolver().ordinal()]);
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
  }
}
