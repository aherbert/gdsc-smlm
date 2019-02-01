/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

import uk.ac.sussex.gdsc.core.ij.ImageJTrackProgress;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog.OptionListener;
import uk.ac.sussex.gdsc.core.utils.BitFlagUtils;
import uk.ac.sussex.gdsc.core.utils.FileUtils;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.ValidationUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraType;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsFileFormat;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsFileSettings;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageMode;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageSettings;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageType;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsInMemorySettings;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsSettings;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsSettings.Builder;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsTableFormat;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsTableSettings;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.ij.IJImageSource;
import uk.ac.sussex.gdsc.smlm.ij.gui.PeakResultTableModel;
import uk.ac.sussex.gdsc.smlm.ij.gui.PeakResultTableModelFrame;
import uk.ac.sussex.gdsc.smlm.ij.results.ImageJImagePeakResults;
import uk.ac.sussex.gdsc.smlm.ij.results.ImageJTablePeakResults;
import uk.ac.sussex.gdsc.smlm.ij.results.ImagePeakResultsFactory;
import uk.ac.sussex.gdsc.smlm.ij.settings.Constants;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.BinaryFilePeakResults;
import uk.ac.sussex.gdsc.smlm.results.ExtendedPeakResult;
import uk.ac.sussex.gdsc.smlm.results.FixedPeakResultList;
import uk.ac.sussex.gdsc.smlm.results.MalkFilePeakResults;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResultsList;
import uk.ac.sussex.gdsc.smlm.results.PeakResultsReader;
import uk.ac.sussex.gdsc.smlm.results.ResultOption;
import uk.ac.sussex.gdsc.smlm.results.TextFilePeakResults;
import uk.ac.sussex.gdsc.smlm.results.TsfPeakResultsWriter;
import uk.ac.sussex.gdsc.smlm.results.count.Counter;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedureX;

import ij.IJ;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.YesNoCancelDialog;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import ij.plugin.frame.Recorder;
import ij.util.Java2;

import java.awt.Checkbox;
import java.awt.Choice;
import java.awt.EventQueue;
import java.awt.Label;
import java.awt.Panel;
import java.awt.Rectangle;
import java.awt.event.ItemListener;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.EnumSet;
import java.util.List;

import javax.swing.JFileChooser;

/**
 * Opens peaks results and displays/converts them.
 */
public class ResultsManager implements PlugIn {

  private static final String TITLE = "Peak Results Manager";

  /** The input file. */
  static final String INPUT_FILE = "File";

  /** The input memory. */
  static final String INPUT_MEMORY = "Memory";

  /** The input none. */
  static final String INPUT_NONE = "[None]";

  private static String inputOption = "";
  private static String inputFilename = Prefs.get(Constants.inputFilename, "");

  private ResultsSettings.Builder resultsSettings = ResultsSettings.newBuilder();
  private boolean extraOptions;

  private boolean fileInput;

  private static double inputNmPerPixel = Prefs.get(Constants.inputNmPerPixel, 0);
  private static double inputGain = Prefs.get(Constants.inputGain, 1);
  private static double inputExposureTime = Prefs.get(Constants.inputExposureTime, 0);
  private static List<String> selected;

  /**
   * The input source.
   */
  public enum InputSource {
    //@formatter:off
    /** File input. */
    FILE{ @Override
    public String getName() { return "File"; }},

    /** Memory input. */
    MEMORY{ @Override
    public String getName() { return "Memory"; }},

    /** Memory input with at least one result spanning frames (linked using ID). */
    MEMORY_MULTI_FRAME{ @Override
    public String getName() { return "Memory (Multi-Frame)"; }},

    /** Memory input with no results spanning frames (linked using ID). */
    MEMORY_SINGLE_FRAME{ @Override
    public String getName() { return "Memory (Single-Frame)"; }},

    /** Memory input with clustered results (ID above zero). */
    MEMORY_CLUSTERED{ @Override
    public String getName() { return "Memory (Clustered)"; }},

    /** No input. */
    NONE{ @Override
    public String getName() { return "None"; }};
    //@formatter:on

    @Override
    public String toString() {
      return getName();
    }

    /**
     * Gets the name.
     *
     * @return the name
     */
    public abstract String getName();
  }

  @Override
  public void run(String arg) {
    extraOptions = ImageJUtils.isExtraOptions();
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if ("load".equals(arg)) {
      batchLoad();
      return;
    }

    if ("save".equals(arg)) {
      batchSave();
      return;
    }

    if (arg != null && arg.startsWith("clear")) {
      if (MemoryPeakResults.isMemoryEmpty()) {
        IJ.error(TITLE, "There are no fitting results in memory");
        IJ.showStatus("");
        return;
      }

      Collection<MemoryPeakResults> allResults;
      boolean removeAll = false;
      if (arg.contains("multi")) {
        final MultiDialog md = new MultiDialog(TITLE, new MultiDialog.MemoryResultsItems());
        md.addSelected(selected);
        md.showDialog();
        if (md.wasCancelled()) {
          return;
        }
        selected = md.getSelectedResults();
        if (selected.isEmpty()) {
          return;
        }
        allResults = new ArrayList<>(selected.size());
        for (final String name : selected) {
          final MemoryPeakResults r = MemoryPeakResults.getResults(name);
          if (r != null) {
            allResults.add(r);
          }
        }
      } else {
        removeAll = true;
        allResults = MemoryPeakResults.getAllResults();
      }
      if (allResults.isEmpty()) {
        return;
      }

      long memorySize = 0;
      int size = 0;
      for (final MemoryPeakResults results : allResults) {
        memorySize += MemoryPeakResults.estimateMemorySize(results);
        size += results.size();
      }
      final String memory = TextUtils.bytesToString(memorySize);
      final String count = TextUtils.pleural(size, "result");
      final String sets = TextUtils.pleural(allResults.size(), "set");

      final GenericDialog gd = new GenericDialog(TITLE);

      gd.addMessage(
          String.format("Do you want to remove %s from memory (%s, %s)?", count, sets, memory));
      gd.showDialog();
      if (gd.wasCanceled()) {
        return;
      }

      if (removeAll) {
        MemoryPeakResults.clearMemory();
      } else {
        for (final MemoryPeakResults results : allResults) {
          MemoryPeakResults.removeResults(results.getName());
        }
      }

      SummariseResults.clearSummaryTable();
      ImageJUtils.log("Cleared %s (%s, %s)", count, sets, memory);
      return;
    }

    if (!showDialog()) {
      return;
    }

    final MemoryPeakResults results = loadResults(inputOption);

    if (results == null || results.size() == 0) {
      IJ.error(TITLE, "No results could be loaded");
      IJ.showStatus("");
      return;
    }

    IJ.showStatus("Loaded " + TextUtils.pleural(results.size(), "result"));

    boolean saved = false;
    if (resultsSettings.getResultsInMemorySettings().getInMemory() && fileInput) {
      MemoryPeakResults.addResults(results);
      saved = true;
    }

    if (resultsSettings.getResultsTableSettings().getResultsTableFormatValue() <= 0
        && resultsSettings.getResultsImageSettings().getImageTypeValue() <= 0
        && TextUtils.isNullOrEmpty(resultsSettings.getResultsFileSettings().getResultsFilename())) {
      // No outputs. Error if results were not saved to memory
      if (!saved) {
        IJ.error(TITLE, "No output selected");
      }
      return;
    }

    final Rectangle bounds = results.getBounds(true);
    final boolean showDeviations =
        resultsSettings.getShowDeviations() && canShowDeviations(results);
    final boolean showEndFrame = canShowEndFrame(results);
    final boolean showId = canShowId(results);

    // Display the configured output
    final PeakResultsList outputList = new PeakResultsList();

    outputList.copySettings(results);
    // String title = results.getSource();
    // if (title == null || title.length() == 0)
    // output.setSource(TITLE);

    final int tableFormat = resultsSettings.getResultsTableSettings().getResultsTableFormatValue();
    if (tableFormat == ResultsTableFormat.IMAGEJ_VALUE) {
      addImageJTableResults(outputList, resultsSettings.getResultsTableSettings(), showDeviations,
          showEndFrame, results.is3D(), showId);
    } else if (tableFormat == ResultsTableFormat.INTERACTIVE_VALUE) {
      showInteractiveTable(results, resultsSettings.getResultsTableSettings());
    }

    addImageResults(outputList, resultsSettings.getResultsImageSettings(), bounds,
        (extraOptions) ? FLAG_EXTRA_OPTIONS : 0);
    addFileResults(outputList, showDeviations, showEndFrame, showId);

    if (outputList.numberOfOutputs() == 0) {
      // This occurs when only using the interactive table
      IJ.showStatus("Processed " + TextUtils.pleural(results.size(), "result"));
      return;
    }

    IJ.showStatus("Processing outputs ...");

    // Reduce to single object for speed
    final PeakResults output =
        (outputList.numberOfOutputs() == 1) ? outputList.toArray()[0] : outputList;

    output.begin();

    // Note: We could add a batch iterator to the MemoryPeakResults.
    // However the speed increase will be marginal as the main time
    // taken is in processing the outputs.

    // Process in batches to provide progress
    final Counter progress = new Counter();
    final int totalProgress = results.size();
    final int batchSize = Math.max(100, totalProgress / 10);
    final FixedPeakResultList batch = new FixedPeakResultList(batchSize);
    IJ.showProgress(0);
    results.forEach((PeakResultProcedureX) result -> {
      batch.add(result);
      if (batch.size == batchSize) {
        if (IJ.escapePressed()) {
          batch.clear();
          return true;
        }
        output.addAll(batch.results);
        batch.clear();
        IJ.showProgress(progress.incrementAndGet(batchSize), totalProgress);
      }
      return false;
    });

    // Will be empty if interrupted
    if (batch.isNotEmpty()) {
      output.addAll(batch.toArray());
    }

    IJ.showProgress(1);
    output.end();

    if (output.size() == results.size()) {
      IJ.showStatus("Processed " + TextUtils.pleural(results.size(), "result"));
    } else {
      IJ.showStatus(
          String.format("A %d/%s", output.size(), TextUtils.pleural(results.size(), "result")));
    }
  }

  private static boolean canShowDeviations(MemoryPeakResults results) {
    return results.hasDeviations();
  }

  private static boolean canShowEndFrame(MemoryPeakResults results) {
    return results.hasEndFrame();
  }

  private static boolean canShowId(MemoryPeakResults results) {
    return results.hasId();
  }

  /**
   * Adds the table results.
   *
   * @param resultsList the results list
   * @param resultsSettings the results settings
   * @param showDeviations the show deviations
   * @param showEndFrame the show end frame
   * @param showZ the show Z
   * @param showId the show id
   * @return the IJ table peak results
   */
  public static ImageJTablePeakResults addTableResults(PeakResultsList resultsList,
      ResultsTableSettings resultsSettings, boolean showDeviations, boolean showEndFrame,
      boolean showZ, boolean showId) {
    if (resultsSettings.getShowTable()) {
      return addImageJTableResults(resultsList, resultsSettings, showDeviations, showEndFrame,
          showZ, showId);
    }
    return null;
  }

  private static ImageJTablePeakResults addImageJTableResults(PeakResultsList resultsList,
      ResultsTableSettings resultsSettings, boolean showDeviations, boolean showEndFrame,
      boolean showZ, boolean showId) {
    final ImageJTablePeakResults r = new ImageJTablePeakResults(showDeviations);
    r.setDistanceUnit(resultsSettings.getDistanceUnit());
    r.setIntensityUnit(resultsSettings.getIntensityUnit());
    r.setAngleUnit(resultsSettings.getAngleUnit());
    r.setShowPrecision(resultsSettings.getShowPrecision());
    if (resultsSettings.getShowPrecision()) {
      r.setComputePrecision(true);
    }
    r.setShowEndFrame(showEndFrame);
    r.setRoundingPrecision(resultsSettings.getRoundingPrecision());
    r.setShowZ(showZ);
    r.setShowFittingData(resultsSettings.getShowFittingData());
    r.setShowNoiseData(resultsSettings.getShowNoiseData());
    r.setShowId(showId);
    resultsList.addOutput(r);
    return r;
  }

  /**
   * Show interactive table.
   *
   * @param results the results
   * @param resultsTableSettings the results table settings
   */
  public static void showInteractiveTable(MemoryPeakResults results,
      ResultsTableSettings resultsTableSettings) {
    final PeakResultTableModel model =
        new PeakResultTableModel(results, true, resultsTableSettings);
    final PeakResultTableModelFrame frame = new PeakResultTableModelFrame(model);
    frame.setTitle(results.getName());
    frame.setVisible(true);
  }

  /**
   * Adds the image results.
   *
   * @param resultsList the results list
   * @param resultsSettings the results settings
   * @param bounds the bounds
   * @param flags the flags
   */
  public static void addImageResults(PeakResultsList resultsList,
      ResultsImageSettings resultsSettings, Rectangle bounds, int flags) {
    if (resultsSettings.getImageTypeValue() > 0) {
      final ImageJImagePeakResults image = ImagePeakResultsFactory.createPeakResultsImage(
          resultsSettings.getImageType(), resultsSettings.getWeighted(),
          resultsSettings.getEqualised(), resultsList.getName(), bounds,
          resultsList.getNmPerPixel(), resultsList.getGain(), resultsSettings.getScale(),
          resultsSettings.getAveragePrecision(), ResultsImageMode.IMAGE_ADD);
      if (BitFlagUtils.anySet(flags, FLAG_EXTRA_OPTIONS)) {
        image.setRollingWindowSize(resultsSettings.getRollingWindowSize());
      }
      image.setRepaintDelay(2000);
      resultsList.addOutput(image);
    }
  }

  private void addFileResults(PeakResultsList resultsList, boolean showDeviations,
      boolean showEndFrame, boolean showId) {
    final ResultsFileSettings resultsFileSettings = this.resultsSettings.getResultsFileSettings();
    if (!TextUtils.isNullOrEmpty(resultsFileSettings.getResultsFilename())) {
      // Remove extension
      final String resultsFilename =
          FileUtils.replaceExtension(resultsFileSettings.getResultsFilename(),
              ResultsProtosHelper.getExtension(resultsFileSettings.getFileFormat()));

      if (fileInput && inputFilename.equals(resultsFilename)) {
        IJ.log(TITLE + ": Input and output files are the same, skipping output ...");
        return;
      }

      // Update
      this.resultsSettings.getResultsFileSettingsBuilder().setResultsFilename(resultsFilename);
      SettingsManager.writeSettings(this.resultsSettings.build());

      // Check if file exists
      final File file = new File(resultsFilename);
      if (file.exists()) {
        final YesNoCancelDialog d = new YesNoCancelDialog(IJ.getInstance(), TITLE,
            "Overwrite existing file?\n" + resultsFilename);
        if (!d.yesPressed()) {
          return;
        }
      }

      addFileResults(resultsList, resultsFileSettings, resultsFilename, showDeviations,
          showEndFrame, showId);
    }
  }

  /**
   * Adds the file results.
   *
   * @param resultsList the results list
   * @param resultsSettings the results settings
   * @param resultsFilename the results filename
   * @param showDeviations the show deviations
   * @param showEndFrame the show end frame
   * @param showId the show id
   * @return the peak results
   */
  public static PeakResults addFileResults(PeakResultsList resultsList,
      ResultsFileSettings resultsSettings, String resultsFilename, boolean showDeviations,
      boolean showEndFrame, boolean showId) {
    if (resultsSettings.getFileFormatValue() > 0 && resultsFilename != null) {
      final File file = new File(resultsFilename);
      final File parent = file.getParentFile();
      if (parent != null && parent.exists()) {
        PeakResults results;
        switch (resultsSettings.getFileFormat()) {
          case BINARY:
            results = new BinaryFilePeakResults(resultsFilename, showDeviations, showEndFrame,
                showId, resultsSettings.getShowPrecision());
            break;
          case TEXT:
            final TextFilePeakResults f = new TextFilePeakResults(resultsFilename, showDeviations,
                showEndFrame, showId, resultsSettings.getShowPrecision());
            f.setDistanceUnit(resultsSettings.getDistanceUnit());
            f.setIntensityUnit(resultsSettings.getIntensityUnit());
            f.setAngleUnit(resultsSettings.getAngleUnit());
            f.setComputePrecision(true);
            results = f;
            break;
          case MALK:
            results = new MalkFilePeakResults(resultsFilename);
            break;
          case TSF:
            results = new TsfPeakResultsWriter(resultsFilename);
            break;
          default:
            throw new IllegalArgumentException(
                "Unsupported file format: " + resultsSettings.getFileFormat());
        }
        resultsList.addOutput(results);
        return results;
      }
    }
    return null;
  }

  private boolean showDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);

    resultsSettings = SettingsManager.readResultsSettings(0).toBuilder();

    gd.addMessage("Read the Peak Results and output to a new format");

    gd.addMessage("Select the Peak Results");
    addInput(gd, inputOption, InputSource.MEMORY, InputSource.FILE);

    final Choice inputChoice = gd.getLastChoice();

    addTableResultsOptions(gd, resultsSettings, FLAG_TABLE_FORMAT);
    addImageResultsOptions(gd, resultsSettings);
    addFileResultsOptions(gd, resultsSettings, 0);
    addInMemoryResultsOptions(gd, resultsSettings);

    final Label messageLabel = (Label) gd.getMessage();
    final Checkbox saveCheckbox = gd.getLastCheckbox();

    // Hide the in-memory settings if the input is not a file
    if (ImageJUtils.isShowGenericDialog()) {
      final Label saveLabel = gd.getLastLabel();
      final ItemListener listener = event -> {
        final boolean enable = INPUT_FILE.equals(inputChoice.getSelectedItem());
        if (enable != messageLabel.isVisible()) {
          messageLabel.setVisible(enable);
          saveCheckbox.setVisible(enable);
          saveLabel.setVisible(enable);
          gd.pack();
        }
      };

      // Run once to set up
      listener.itemStateChanged(null);

      inputChoice.addItemListener(listener);
    }

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    inputOption = ResultsManager.getInputSource(gd);
    inputFilename = gd.getNextString();
    resultsSettings.getResultsTableSettingsBuilder()
        .setResultsTableFormatValue(gd.getNextChoiceIndex());
    resultsSettings.getResultsImageSettingsBuilder().setImageTypeValue(gd.getNextChoiceIndex());
    resultsSettings.getResultsFileSettingsBuilder().setFileFormatValue(gd.getNextChoiceIndex());
    resultsSettings.getResultsFileSettingsBuilder().setResultsFilename(gd.getNextString());
    resultsSettings.getResultsInMemorySettingsBuilder().setInMemory(gd.getNextBoolean());

    gd.collectOptions();

    // Check arguments
    try {
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

    Prefs.set(Constants.inputFilename, inputFilename);

    SettingsManager.writeSettings(resultsSettings.build());

    return true;
  }

  /** Use this to add extra options to the dialog. */
  public static final int FLAG_EXTRA_OPTIONS = 0x00000001;
  /** Use this to add the results directory to the file results dialog. */
  public static final int FLAG_RESULTS_DIRECTORY = 0x00000002;
  /** Use this to add the results file to the file results dialog. */
  public static final int FLAG_RESULTS_FILE = 0x00000004;
  /** Use this to avoid adding the section header to the dialog. */
  public static final int FLAG_NO_SECTION_HEADER = 0x00000008;
  /** Use this to add a choice of table format to the dialog. */
  public static final int FLAG_TABLE_FORMAT = 0x00000010;

  /**
   * Adds the table results options.
   *
   * @param gd the dialog
   * @param resultsSettings the results settings
   */
  public static void addTableResultsOptions(final ExtendedGenericDialog gd,
      final Builder resultsSettings) {
    addTableResultsOptions(gd, resultsSettings, 0);
  }

  /**
   * Adds the table results options.
   *
   * @param gd the dialog
   * @param resultsSettings the results settings
   * @param flags the flags
   */
  public static void addTableResultsOptions(final ExtendedGenericDialog gd,
      final Builder resultsSettings, final int flags) {
    if (BitFlagUtils.anyNotSet(flags, FLAG_NO_SECTION_HEADER)) {
      gd.addMessage("--- Table output ---");
    }
    final ResultsTableSettings.Builder tableSettings =
        resultsSettings.getResultsTableSettingsBuilder();
    if (BitFlagUtils.anySet(flags, FLAG_TABLE_FORMAT)) {
      gd.addChoice("Table", SettingsManager.getResultsTableFormatNames(),
          tableSettings.getResultsTableFormatValue(), new OptionListener<Integer>() {
            @Override
            public boolean collectOptions(Integer field) {
              tableSettings.setResultsTableFormatValue(field);
              return collectOptions(false);
            }

            @Override
            public boolean collectOptions() {
              return collectOptions(true);
            }

            private boolean collectOptions(boolean silent) {
              if (tableSettings.getResultsTableFormatValue() <= 0) {
                return false;
              }
              final ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE, null);
              egd.addChoice("Table_distance_unit", SettingsManager.getDistanceUnitNames(),
                  tableSettings.getDistanceUnitValue());
              egd.addChoice("Table_intensity_unit", SettingsManager.getIntensityUnitNames(),
                  tableSettings.getIntensityUnitValue());
              egd.addChoice("Table_angle_unit", SettingsManager.getAngleUnitNames(),
                  tableSettings.getAngleUnitValue());
              egd.addCheckbox("Table_show_fitting_data", tableSettings.getShowFittingData());
              egd.addCheckbox("Table_show_noise_data", tableSettings.getShowNoiseData());
              egd.addCheckbox("Table_show_precision", tableSettings.getShowPrecision());
              egd.addSlider("Table_precision", 0, 10, tableSettings.getRoundingPrecision());
              egd.setSilent(silent);
              egd.showDialog(true, gd);
              if (egd.wasCanceled()) {
                return false;
              }
              tableSettings.setDistanceUnitValue(egd.getNextChoiceIndex());
              tableSettings.setIntensityUnitValue(egd.getNextChoiceIndex());
              tableSettings.setAngleUnitValue(egd.getNextChoiceIndex());
              tableSettings.setShowFittingData(egd.getNextBoolean());
              tableSettings.setShowNoiseData(egd.getNextBoolean());
              tableSettings.setShowPrecision(egd.getNextBoolean());
              tableSettings.setRoundingPrecision((int) egd.getNextNumber());
              return true;
            }
          });
    } else {
      gd.addCheckbox("Show_results_table", tableSettings.getShowTable(),
          new OptionListener<Boolean>() {

            @Override
            public boolean collectOptions(Boolean field) {
              tableSettings.setShowTable(field);
              return collectOptions(false);
            }

            @Override
            public boolean collectOptions() {
              return collectOptions(true);
            }

            private boolean collectOptions(boolean silent) {
              if (!tableSettings.getShowTable()) {
                return false;
              }
              final ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE, null);
              egd.addChoice("Table_distance_unit", SettingsManager.getDistanceUnitNames(),
                  tableSettings.getDistanceUnitValue());
              egd.addChoice("Table_intensity_unit", SettingsManager.getIntensityUnitNames(),
                  tableSettings.getIntensityUnitValue());
              egd.addChoice("Table_angle_unit", SettingsManager.getAngleUnitNames(),
                  tableSettings.getAngleUnitValue());
              egd.addCheckbox("Table_show_fitting_data", tableSettings.getShowFittingData());
              egd.addCheckbox("Table_show_noise_data", tableSettings.getShowNoiseData());
              egd.addCheckbox("Table_show_precision", tableSettings.getShowPrecision());
              egd.addSlider("Table_precision", 0, 10, tableSettings.getRoundingPrecision());
              egd.setSilent(silent);
              egd.showDialog(true, gd);
              if (egd.wasCanceled()) {
                return false;
              }
              tableSettings.setDistanceUnitValue(egd.getNextChoiceIndex());
              tableSettings.setIntensityUnitValue(egd.getNextChoiceIndex());
              tableSettings.setAngleUnitValue(egd.getNextChoiceIndex());
              tableSettings.setShowFittingData(egd.getNextBoolean());
              tableSettings.setShowNoiseData(egd.getNextBoolean());
              tableSettings.setShowPrecision(egd.getNextBoolean());
              tableSettings.setRoundingPrecision((int) egd.getNextNumber());
              return true;
            }

          });
    }
  }

  private void addImageResultsOptions(final ExtendedGenericDialog gd,
      final Builder resultsSettings) {
    addImageResultsOptions(gd, resultsSettings, (extraOptions) ? FLAG_EXTRA_OPTIONS : 0);
  }

  /**
   * Adds the image results options.
   *
   * @param gd the dialog
   * @param resultsSettings the results settings
   * @param flags the flags
   */
  public static void addImageResultsOptions(final ExtendedGenericDialog gd,
      final Builder resultsSettings, final int flags) {
    if (BitFlagUtils.anyNotSet(flags, FLAG_NO_SECTION_HEADER)) {
      gd.addMessage("--- Image output ---");
    }
    final ResultsImageSettings.Builder imageSettings =
        resultsSettings.getResultsImageSettingsBuilder();
    final EnumSet<ResultsImageType> requirePrecision =
        EnumSet.of(ResultsImageType.DRAW_LOCALISATIONS_AVERAGE_PRECISION,
            ResultsImageType.DRAW_INTENSITY_AVERAGE_PRECISION);
    final EnumSet<ResultsImageType> requireWeighted =
        EnumSet.of(ResultsImageType.DRAW_LOCALISATIONS, ResultsImageType.DRAW_INTENSITY,
            ResultsImageType.DRAW_FRAME_NUMBER, ResultsImageType.DRAW_FIT_ERROR);
    gd.addChoice("Image", SettingsManager.getResultsImageTypeNames(),
        imageSettings.getImageTypeValue(), new OptionListener<Integer>() {
          @Override
          public boolean collectOptions(Integer field) {
            imageSettings.setImageTypeValue(field);
            return collectOptions(false);
          }

          @Override
          public boolean collectOptions() {
            return collectOptions(true);
          }

          private boolean collectOptions(boolean silent) {
            final ResultsImageType resultsImage = imageSettings.getImageType();
            if (resultsImage.getNumber() <= 0) {
              return false;
            }
            final boolean extraOptions = BitFlagUtils.anySet(flags, FLAG_EXTRA_OPTIONS);
            final ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE, null);
            if (requireWeighted.contains(resultsImage)) {
              egd.addCheckbox("Weighted", imageSettings.getWeighted());
            }
            egd.addCheckbox("Equalised", imageSettings.getEqualised());
            if (requirePrecision.contains(resultsImage)) {
              egd.addSlider("Image_Precision (nm)", 5, 30, imageSettings.getAveragePrecision());
            }
            egd.addSlider("Image_Scale", 1, 15, imageSettings.getScale());
            if (extraOptions) {
              egd.addNumericField("Image_Window", imageSettings.getRollingWindowSize(), 0);
            }
            egd.setSilent(silent);
            egd.showDialog(true, gd);
            if (egd.wasCanceled()) {
              return false;
            }
            if (requireWeighted.contains(resultsImage)) {
              imageSettings.setWeighted(egd.getNextBoolean());
            }
            imageSettings.setEqualised(egd.getNextBoolean());
            if (requirePrecision.contains(resultsImage)) {
              imageSettings.setAveragePrecision(egd.getNextNumber());
            }
            imageSettings.setScale(egd.getNextNumber());
            if (extraOptions) {
              imageSettings.setRollingWindowSize((int) egd.getNextNumber());
            }
            return true;
          }
        });
  }

  /**
   * Adds the file results options.
   *
   * @param gd the dialog
   * @param resultsSettings the results settings
   * @param flags the flags
   */
  public static void addFileResultsOptions(final ExtendedGenericDialog gd,
      final Builder resultsSettings, final int flags) {
    if (BitFlagUtils.anyNotSet(flags, FLAG_NO_SECTION_HEADER)) {
      gd.addMessage("--- File output ---");
    }
    final ResultsFileSettings.Builder fileSettings =
        resultsSettings.getResultsFileSettingsBuilder();
    gd.addChoice("Results_format", SettingsManager.getResultsFileFormatNames(),
        fileSettings.getFileFormatValue(), new OptionListener<Integer>() {
          @Override
          public boolean collectOptions(Integer field) {
            fileSettings.setFileFormatValue(field);
            return collectOptions(false);
          }

          @Override
          public boolean collectOptions() {
            return collectOptions(true);
          }

          private boolean collectOptions(boolean silent) {
            final ResultsFileFormat resultsFileFormat = fileSettings.getFileFormat();
            if (!ResultsProtosHelper.isGdsc(resultsFileFormat)) {
              return false;
            }
            final ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE, null);
            if (resultsFileFormat == ResultsFileFormat.TEXT) {
              egd.addChoice("File_distance_unit", SettingsManager.getDistanceUnitNames(),
                  fileSettings.getDistanceUnit().ordinal());
              egd.addChoice("File_intensity_unit", SettingsManager.getIntensityUnitNames(),
                  fileSettings.getIntensityUnit().ordinal());
              egd.addChoice("File_angle_unit", SettingsManager.getAngleUnitNames(),
                  fileSettings.getAngleUnit().ordinal());
              egd.addCheckbox("File_show_precision", fileSettings.getShowPrecision());
            }
            egd.addCheckbox("Show_deviations", resultsSettings.getShowDeviations());
            egd.setSilent(silent);
            egd.showDialog(true, gd);
            if (egd.wasCanceled()) {
              return false;
            }
            fileSettings.setDistanceUnitValue(egd.getNextChoiceIndex());
            fileSettings.setIntensityUnitValue(egd.getNextChoiceIndex());
            fileSettings.setAngleUnitValue(egd.getNextChoiceIndex());
            fileSettings.setShowPrecision(egd.getNextBoolean());
            resultsSettings.setShowDeviations(egd.getNextBoolean());
            return true;
          }
        });
    if (BitFlagUtils.anySet(flags, FLAG_RESULTS_DIRECTORY)) {
      gd.addDirectoryField("Results_directory", fileSettings.getResultsDirectory());
    } else if (BitFlagUtils.anySet(flags, FLAG_RESULTS_FILE)) {
      gd.addFilenameField("Results_file", fileSettings.getResultsFilename());
    } else {
      // Do not add a selected results file to prevent constant overwrite messages.
      // However we can set the initial directory if we have a results file.
      if (fileSettings.getResultsFilename() != null) {
        final File dir = new File(fileSettings.getResultsFilename()).getParentFile();
        if (dir != null) {
          OpenDialog.setDefaultDirectory(dir.getPath());
        }
      }
      gd.addFilenameField("Results_file", "");
    }
  }

  /**
   * Adds the in memory results options.
   *
   * @param gd the dialog
   * @param resultsSettings the results settings
   */
  public static void addInMemoryResultsOptions(final ExtendedGenericDialog gd,
      final Builder resultsSettings) {
    addInMemoryResultsOptions(gd, resultsSettings, 0);
  }

  /**
   * Adds the in memory results options.
   *
   * @param gd the dialog
   * @param resultsSettings the results settings
   * @param flags the flags
   */
  public static void addInMemoryResultsOptions(final ExtendedGenericDialog gd,
      final Builder resultsSettings, int flags) {
    if (BitFlagUtils.anyNotSet(flags, FLAG_NO_SECTION_HEADER)) {
      gd.addMessage("--- Memory output ---");
    }
    final ResultsInMemorySettings.Builder memorySettings =
        resultsSettings.getResultsInMemorySettingsBuilder();
    gd.addCheckbox("Save_to_memory", memorySettings.getInMemory());
  }

  /**
   * Add a list of input sources to the generic dialog. The choice field will be named 'input'. If a
   * file input option is selected then a field will be added name 'Input_file'.
   *
   * <p>If the source is a memory source then it will not be added if it is empty. If not empty then
   * a summary of the number of localisation is added as a message to the dialog.
   *
   * @param gd the dialog
   * @param inputOption the input option
   * @param inputs the inputs
   */
  public static void addInput(ExtendedGenericDialog gd, String inputOption, InputSource... inputs) {
    addInput(gd, "Input", inputOption, inputs);
  }

  /**
   * Add a list of input sources to the generic dialog. The choice field will be named inputName. If
   * a file input option is selected then a field will be added name 'Input_file'.
   *
   * <p>If the source is a memory source then it will not be added if it is empty. If not empty then
   * a summary of the number of localisation is added as a message to the dialog.
   *
   * @param gd the dialog
   * @param inputName the input name
   * @param inputOption the input option
   * @param inputs the inputs
   */
  public static void addInput(ExtendedGenericDialog gd, String inputName, String inputOption,
      InputSource... inputs) {
    final ArrayList<String> source = new ArrayList<>(3);
    boolean fileInput = false;
    for (final InputSource input : inputs) {
      ResultsManager.addInputSource(source, input);
      if (input == InputSource.FILE) {
        fileInput = true;
      }
    }
    if (source.isEmpty()) {
      addInputSource(source, InputSource.NONE);
    }

    addInputSourceToDialog(gd, inputName, inputOption, source, fileInput);
  }

  /**
   * Add a list of input sources to the generic dialog. The choice field will be named inputName. If
   * the file input option is true then a field will be added name 'Input_file'.
   *
   * @param gd the dialog
   * @param inputName the input name
   * @param inputOption The option to select by default
   * @param source the source
   * @param fileInput the file input
   */
  @SuppressWarnings("unused")
  public static void addInputSourceToDialog(final ExtendedGenericDialog gd, String inputName,
      String inputOption, List<String> source, boolean fileInput) {
    final String[] options = source.toArray(new String[source.size()]);
    // Find the option
    inputOption = removeFormatting(inputOption);

    int optionIndex = 0;
    for (int i = 0; i < options.length; i++) {
      final String name = removeFormatting(options[i]);
      if (name.equals(inputOption)) {
        optionIndex = i;
        break;
      }
    }
    final Choice c = gd.addAndGetChoice(inputName, options, options[optionIndex]);
    if (fileInput) {
      // final TextField tf = gd.addFilenameField("Input_file", inputFilename)
      // final JButton b = gd.getLastOptionButton()

      // Add a listener to the choice to enable the file input field.
      // Currently we hide the filename field and pack the dialog.
      // We may wish to just disable the fields and leave them there.
      // This could be a user configured option in a global GDSC settings class.
      if (ImageJUtils.isShowGenericDialog()) {
        final Label l = gd.getLastLabel();
        final Panel p = gd.getLastPanel();
        final ItemListener listener = event -> {
          final boolean enable = INPUT_FILE.equals(c.getSelectedItem());
          if (enable != l.isVisible()) {
            l.setVisible(enable);
            p.setVisible(enable);
            // tf.setVisible(enable);
            // b.setVisible(enable);
            gd.pack();
          }
        };

        // Run once to set up
        listener.itemStateChanged(null);

        c.addItemListener(listener);
      }
    }
  }

  /**
   * Remove the extra information added to a name for use in dialogs.
   *
   * @param name the formatted name
   * @return The name
   */
  public static String removeFormatting(String name) {
    final int index = name.lastIndexOf('[');
    if (index > 0) {
      name = name.substring(0, index - 1);
    }
    return name;
  }

  /**
   * Add an input source the list. If the source is a memory source then it will not be added if it
   * is empty. If not empty then a summary of the number of localisation is added as a message to
   * the dialog.
   *
   * @param source the source
   * @param input the input
   */
  public static void addInputSource(List<String> source, InputSource input) {
    switch (input) {
      case FILE:
        source.add(INPUT_FILE);
        break;

      case MEMORY:
      case MEMORY_MULTI_FRAME:
      case MEMORY_SINGLE_FRAME:
      case MEMORY_CLUSTERED:
        for (final String name : MemoryPeakResults.getResultNames()) {
          addInputSource(source, MemoryPeakResults.getResults(name), input);
        }
        break;

      default:
        ValidationUtils.checkArgument(input == InputSource.NONE, input);
        source.add(INPUT_NONE);
        break;
    }
  }

  /**
   * Add a memory input source to the list.
   *
   * @param source the source
   * @param memoryResults the memory results
   * @param input MEMORY_MULTI_FRAME : Select only those results with at least one result spanning
   *        frames, MEMORY_CLUSTERED : Select only those results which have at least some IDs above
   *        zero (allowing zero to be a valid cluster Id for no cluster)
   */
  public static void addInputSource(List<String> source, MemoryPeakResults memoryResults,
      InputSource input) {
    if (memoryResults.size() > 0) {
      switch (input) {
        case MEMORY_MULTI_FRAME:
          if (!isMultiFrame(memoryResults)) {
            return;
          }
          break;
        case MEMORY_SINGLE_FRAME:
          if (isMultiFrame(memoryResults)) {
            return;
          }
          break;
        case MEMORY_CLUSTERED:
          if (!hasId(memoryResults)) {
            return;
          }
          break;
        default:
      }
      source.add(getName(memoryResults));
    }
  }

  /**
   * Get the name of the results for use in dialogs.
   *
   * @param memoryResults the memory results
   * @return The name
   */
  public static String getName(MemoryPeakResults memoryResults) {
    return memoryResults.getName() + " [" + memoryResults.size() + "]";
  }

  /**
   * Load results from memory using a name from a dialog.
   *
   * @param name the name
   * @return The results
   */
  public static MemoryPeakResults loadMemoryResults(String name) {
    return MemoryPeakResults.getResults(removeFormatting(name));
  }

  /**
   * Check for multi-frame results.
   *
   * @param memoryResults the memory results
   * @return True if at least one result spanning frames
   */
  public static boolean isMultiFrame(MemoryPeakResults memoryResults) {
    return memoryResults
        .forEach((PeakResultProcedureX) result -> result.getFrame() < result.getEndFrame());
  }

  /**
   * Check for any IDs above zero.
   *
   * @param memoryResults the memory results
   * @return True if any results have IDs above zero
   */
  public static boolean hasId(MemoryPeakResults memoryResults) {
    return memoryResults.forEach((PeakResultProcedureX) result -> result.getId() > 0);
  }

  /**
   * Check for all IDs above zero.
   *
   * @param memoryResults the memory results
   * @return True if all results have IDs above zero
   */
  public static boolean isId(MemoryPeakResults memoryResults) {
    return memoryResults.forEach((PeakResultProcedureX) result -> result.getId() <= 0);
  }

  /**
   * All results must be an ExtendedPeakResult.
   *
   * @param memoryResults the memory results
   * @return True if all are an ExtendedPeakResult
   */
  public static boolean isExtended(MemoryPeakResults memoryResults) {
    return memoryResults
        .forEach((PeakResultProcedureX) result -> !(result instanceof ExtendedPeakResult));
  }

  /**
   * Gets the name of the next input source from the dialog.
   *
   * @param gd the dialog
   * @return the input source
   */
  public static String getInputSource(GenericDialog gd) {
    final String source = gd.getNextChoice();
    return removeFormatting(source);
  }

  /**
   * Load the results from the named input option. If the results are not empty then a check can be
   * made for calibration, and data using the legacy standard units (distance in Pixel and intensity
   * in Count).
   *
   * <p>If the calibration cannot be obtained or the units are incorrect then the null will be
   * returned.
   *
   * @param inputOption the input option
   * @param checkCalibration Set to true to ensure the results have a valid calibration
   * @return the results
   * @deprecated Plugins should specify if they require results in a specific unit
   */
  @Deprecated
  public static MemoryPeakResults loadInputResults(String inputOption, boolean checkCalibration) {
    return loadInputResults(inputOption, checkCalibration, DistanceUnit.PIXEL,
        IntensityUnit.PHOTON);
  }

  /**
   * Load the results from the named input option. If the results are not empty then a check can be
   * made for calibration, and data using the specified units. If the calibration cannot be obtained
   * or the units are incorrect then the null will be returned.
   *
   * @param inputOption the input option
   * @param checkCalibration Set to true to ensure the results have a valid calibration
   * @param distanceUnit the required distance unit for the results
   * @return the results
   */
  public static MemoryPeakResults loadInputResults(String inputOption, boolean checkCalibration,
      DistanceUnit distanceUnit) {
    return loadInputResults(inputOption, checkCalibration, distanceUnit, null);
  }

  /**
   * Load the results from the named input option. If the results are not empty then a check can be
   * made for calibration, and data using the specified units. If the calibration cannot be obtained
   * or the units are incorrect then the null will be returned.
   *
   * @param inputOption the input option
   * @param checkCalibration Set to true to ensure the results have a valid calibration
   * @param intensityUnit the required intensity unit for the results
   * @return the results
   */
  public static MemoryPeakResults loadInputResults(String inputOption, boolean checkCalibration,
      IntensityUnit intensityUnit) {
    return loadInputResults(inputOption, checkCalibration, null, intensityUnit);
  }

  /**
   * Load the results from the named input option. If the results are not empty then a check can be
   * made for calibration, and data using the specified units. If the calibration cannot be obtained
   * or the units are incorrect then the null will be returned.
   *
   * @param inputOption the input option
   * @param checkCalibration Set to true to ensure the results have a valid calibration
   * @param distanceUnit the required distance unit for the results
   * @param intensityUnit the required intensity unit for the results
   * @return the results
   */
  public static MemoryPeakResults loadInputResults(String inputOption, boolean checkCalibration,
      DistanceUnit distanceUnit, IntensityUnit intensityUnit) {
    MemoryPeakResults results = null;
    PeakResultsReader reader = null;
    if (inputOption.equals(INPUT_NONE)) {
      // Do nothing
    } else if (inputOption.equals(INPUT_FILE)) {
      IJ.showStatus("Reading results file ...");
      reader = new PeakResultsReader(inputFilename);
      IJ.showStatus("Reading " + reader.getFormat() + " results file ...");
      final ResultOption[] options = reader.getOptions();
      if (options != null) {
        collectOptions(reader, options);
      }
      reader.setTracker(new ImageJTrackProgress());
      results = reader.getResults();
      reader.getTracker().progress(1.0);

      if (results != null && results.size() > 0) {
        // If the name contains a .tif suffix then create an image source
        if (results.getName() != null && results.getName().contains(".tif")
            && results.getSource() == null) {
          final int index = results.getName().indexOf(".tif");
          results.setSource(new IJImageSource(results.getName().substring(0, index)));
        }
      }
    } else {
      results = loadMemoryResults(inputOption);
    }

    try {
      if (results == null) {
        return null;
      }
      if (results.isEmpty()) {
        return results;
      }

      if (checkCalibration && !checkCalibration(results, reader)) {
        return null;
      }
      if (distanceUnit != null && results.getDistanceUnit() != distanceUnit) {
        ImageJUtils.log("Incorrect distance unit: " + results.getDistanceUnit());
        return null;
      }
      if (intensityUnit != null && results.getIntensityUnit() != intensityUnit) {
        ImageJUtils.log("Incorrect intensity unit: " + results.getDistanceUnit());
        return null;
      }
    } finally {
      IJ.showStatus("");
    }
    return results;
  }

  private static void collectOptions(PeakResultsReader reader, ResultOption[] options) {
    final GenericDialog gd = new GenericDialog(TITLE);
    gd.addMessage("Options required for file format: " + reader.getFormat().getName());
    for (final ResultOption option : options) {
      if (option.hasValues()) {
        final String[] items = new String[option.values.length];
        for (int i = 0; i < items.length; i++) {
          items[i] = option.values[i].toString();
        }
        gd.addChoice(getOptionName(option), items, option.getValue().toString());
      } else if (option.getValue() instanceof Number) {
        final Number n = (Number) option.getValue();
        if (n.doubleValue() == n.intValue()) {
          gd.addNumericField(getOptionName(option), n.intValue(), 0);
        } else {
          final String value = n.toString();
          int sig = 0;
          int index = value.indexOf('.');
          if (index != -1) {
            // There is a decimal point. Count the digits after it
            while (++index < value.length()) {
              if (!Character.isDigit(value.charAt(index))) {
                // A non-digit character after the decimal point is for scientific notation
                sig = -sig;
                break;
              }
              sig++;
            }
          }
          gd.addNumericField(getOptionName(option), n.doubleValue(), sig);
        }
      } else if (option.getValue() instanceof String) {
        gd.addStringField(getOptionName(option), (String) option.getValue());
      } else if (option.getValue() instanceof Boolean) {
        gd.addCheckbox(getOptionName(option), (Boolean) option.getValue());
      } else {
        IJ.log(TITLE + ": Unsupported reader option: " + option.name + "="
            + option.getValue().toString());
      }
    }
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    try {
      for (final ResultOption option : options) {
        if (option.hasValues()) {
          option.setValue(option.values[gd.getNextChoiceIndex()]);
        } else if (option.getValue() instanceof Number) {
          final double d = gd.getNextNumber();
          // Convert to the correct type using the String value constructor for the number
          option.setValue(option.getValue().getClass().getConstructor(String.class)
              .newInstance(Double.toString(d)));
        } else if (option.getValue() instanceof String) {
          option.setValue(gd.getNextString());
        } else if (option.getValue() instanceof Boolean) {
          option.setValue(gd.getNextBoolean());
        }
      }
      reader.setOptions(options);
    } catch (final Exception ex) {
      // This can occur if the options are not valid
      IJ.log(TITLE + ": Failed to configure reader options: " + ex.getMessage());
    }
  }

  private static String getOptionName(ResultOption option) {
    return option.name.replace(' ', '_');
  }

  /**
   * Check the calibration of the results exists, if not then prompt for it with a dialog.
   *
   * @param results The results
   * @return True if OK; false if calibration dialog cancelled
   */
  public static boolean checkCalibration(MemoryPeakResults results) {
    return checkCalibration(results, null);
  }

  /**
   * Check the calibration of the results exists, if not then prompt for it with a dialog.
   *
   * <p>The calibration is rechecked after the dialog is shown.
   *
   * @param results The results
   * @param reader Used to determine the file type
   * @return True if OK; false if calibration is missing
   */
  private static boolean checkCalibration(MemoryPeakResults results, PeakResultsReader reader) {
    // Check for Calibration
    final CalibrationWriter calibration = results.getCalibrationWriterSafe();
    String msg = "partially calibrated";
    if (!results.hasCalibration()) {
      // Make sure the user knows all the values have not been set
      msg = "uncalibrated";
    }

    boolean missing = isEssentialCalibrationMissing(calibration);

    if (missing) {
      final Rectangle2D.Float dataBounds = results.getDataBounds(null);

      final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
      gd.addMessage(String.format("Results are %s.\nData bounds = (%s,%s) to (%s,%s)", msg,
          MathUtils.rounded(dataBounds.x), MathUtils.rounded(dataBounds.y),
          MathUtils.rounded(dataBounds.y + dataBounds.getWidth()),
          MathUtils.rounded(dataBounds.x + dataBounds.getHeight())));
      gd.addChoice("Calibration_distance_unit", SettingsManager.getDistanceUnitNames(),
          calibration.getDistanceUnitValue());
      gd.addChoice("Calibration_intensity_unit", SettingsManager.getIntensityUnitNames(),
          calibration.getIntensityUnitValue());
      gd.addNumericField("Calibration (nm/px)", calibration.getNmPerPixel(), 2);
      gd.addNumericField("Exposure_time (ms)", calibration.getExposureTime(), 2);
      PeakFit.addCameraOptions(gd, calibration);
      gd.showDialog();
      if (gd.wasCanceled()) {
        return false;
      }
      calibration.setDistanceUnit(SettingsManager.getDistanceUnitValues()[gd.getNextChoiceIndex()]);
      calibration
          .setIntensityUnit(SettingsManager.getIntensityUnitValues()[gd.getNextChoiceIndex()]);
      calibration.setNmPerPixel(Math.abs(gd.getNextNumber()));
      calibration.setExposureTime(Math.abs(gd.getNextNumber()));
      calibration.setCameraType(SettingsManager.getCameraTypeValues()[gd.getNextChoiceIndex()]);

      gd.collectOptions();

      if (calibration.getCameraType() == CameraType.SCMOS) {
        calibration.clearGlobalCameraSettings();
      }

      missing = isEssentialCalibrationMissing(calibration);

      // Save for next time ...
      inputNmPerPixel = calibration.getNmPerPixel();
      inputExposureTime = calibration.getExposureTime();
      Prefs.set(Constants.inputNmPerPixel, inputNmPerPixel);
      Prefs.set(Constants.inputExposureTime, inputExposureTime);
      if (calibration.isCcdCamera()) {
        inputGain = calibration.getCountPerPhoton();
        Prefs.set(Constants.inputGain, inputGain);
      }

      results.setCalibration(calibration.getCalibration());
    }

    // Only OK if nothing is missing
    return !missing;
  }

  /**
   * Check for essential calibration settings (i.e. not readNoise, bias, emCCD, amplification).
   *
   * @param calibration the calibration
   * @return true, if calibration is missing
   */
  private static boolean isEssentialCalibrationMissing(final CalibrationWriter calibration) {
    boolean missing = false;
    if (!calibration.hasNmPerPixel()) {
      missing = true;
      calibration.setNmPerPixel(inputNmPerPixel);
    }
    if (!calibration.hasExposureTime()) {
      missing = true;
      calibration.setExposureTime(inputExposureTime);
    }
    if (!calibration.hasDistanceUnit()) {
      missing = true;
    }
    if (!calibration.hasIntensityUnit()) {
      missing = true;
    }

    switch (calibration.getCameraType()) {
      case CCD:
      case EMCCD:
        // Count-per-photon is required for CCD camera types
        missing |= !calibration.hasCountPerPhoton();
        calibration.setCountPerPhoton(inputGain);
        break;
      case SCMOS:
        break;
      case CAMERA_TYPE_NA:
      case UNRECOGNIZED:
      default:
        missing = true;
    }
    return missing;
  }

  /**
   * Load the results from the named input option.
   *
   * @param inputOption the input option
   * @return the memory peak results
   */
  private MemoryPeakResults loadResults(String inputOption) {
    if (inputOption.equals(INPUT_FILE)) {
      fileInput = true;
    }
    return loadInputResults(inputOption, true, null, null);
  }

  /**
   * Gets the input filename that will be used in
   * {@link ResultsManager#loadInputResults(String, boolean)}.
   *
   * @return the input filename
   */
  static String getInputFilename() {
    return inputFilename;
  }

  /**
   * Sets the input filename that will be used in
   * {@link ResultsManager#loadInputResults(String, boolean)}.
   *
   * @param inputFilename the new input filename
   */
  static void setInputFilename(String inputFilename) {
    ResultsManager.inputFilename = inputFilename;
  }

  private String omDirectory;
  private File[] omFiles;

  /**
   * Batch load a set of results files.
   */
  private void batchLoad() {
    // Adapted from ij.io.Opener.openMultiple

    final String resetInputFilename = inputFilename;

    Java2.setSystemLookAndFeel();
    // run JFileChooser in a separate thread to avoid possible thread deadlocks
    try {
      EventQueue.invokeAndWait(() -> {
        final JFileChooser fc = new JFileChooser();
        fc.setMultiSelectionEnabled(true);
        File dir = null;
        final String sdir = OpenDialog.getDefaultDirectory();
        if (sdir != null) {
          dir = new File(sdir);
        }
        if (dir != null) {
          fc.setCurrentDirectory(dir);
        }
        final int returnVal = fc.showOpenDialog(IJ.getInstance());
        if (returnVal != JFileChooser.APPROVE_OPTION) {
          return;
        }
        omFiles = fc.getSelectedFiles();
        if (omFiles.length == 0) { // getSelectedFiles does not work on some JVMs
          omFiles = new File[1];
          omFiles[0] = fc.getSelectedFile();
        }
        omDirectory = fc.getCurrentDirectory().getPath() + File.separator;
      });
    } catch (final Exception ex) {
      // Ignore
    }
    if (omDirectory == null) {
      return;
    }
    OpenDialog.setDefaultDirectory(omDirectory);
    for (int i = 0; i < omFiles.length; i++) {
      final String path = omDirectory + omFiles[i].getName();
      load(path);
    }

    inputFilename = resetInputFilename;
  }

  private static void load(String path) {
    inputFilename = path;
    // Record this as a single load of the results manager.
    // This should support any dialogs that are presented in loadInputResults(...)
    // to get the calibration.
    if (Recorder.record) {
      Recorder.setCommand("Results Manager");
      Recorder.recordOption("input", INPUT_FILE);
      Recorder.recordOption("input_file", path);
      Recorder.recordOption("image",
          SettingsManager.getResultsImageTypeNames()[ResultsImageType.DRAW_NONE_VALUE]);
      Recorder.recordOption("results_file", "[]");
      Recorder.recordOption("save_to_memory");
    }
    final MemoryPeakResults results = loadInputResults(INPUT_FILE, true, null, null);
    if (results == null || results.size() == 0) {
      IJ.error(TITLE, "No results could be loaded from " + path);
    } else {
      if (Recorder.record) {
        Recorder.saveCommand();
      }
      MemoryPeakResults.addResults(results);
    }
  }

  /**
   * Batch save a set of results files.
   */
  private void batchSave() {
    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "No localisations in memory");
      return;
    }
    final MultiDialog md = new MultiDialog(TITLE, new MultiDialog.MemoryResultsItems());
    md.addSelected(selected);
    md.showDialog();
    if (md.wasCancelled()) {
      return;
    }
    selected = md.getSelectedResults();
    if (selected.isEmpty()) {
      return;
    }
    resultsSettings = SettingsManager.readResultsSettings(0).toBuilder();
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    addFileResultsOptions(gd, resultsSettings, FLAG_RESULTS_DIRECTORY);
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    gd.collectOptions();
    final ResultsFileSettings.Builder b = resultsSettings.getResultsFileSettingsBuilder();
    final String dir = gd.getNextString();
    b.setFileFormatValue(gd.getNextChoiceIndex());
    b.setResultsDirectory(dir);
    SettingsManager.writeSettings(resultsSettings);
    if (dir == null || !new File(dir).exists()) {
      IJ.error(TITLE, "Output directory does not exist");
      return;
    }
    final ResultsFileSettings resultsFileSettings = resultsSettings.getResultsFileSettings();
    if (resultsFileSettings.getFileFormatValue() <= 0) {
      IJ.error(TITLE, "No output file format");
      return;
    }
    int count = 0;
    for (final String name : selected) {
      final MemoryPeakResults r = MemoryPeakResults.getResults(name);
      if (r != null && save(resultsFileSettings, r)) {
        count++;
      }
    }
    IJ.showStatus("Saved " + TextUtils.pleural(count, "dataset"));
  }

  private static boolean save(ResultsFileSettings resultsSettings, MemoryPeakResults source) {
    // Assume the directory exists
    String resultsFilename;
    try {
      resultsFilename = new File(resultsSettings.getResultsDirectory(),
          source.getName() + ".results."
              + ResultsProtosHelper.getExtension(resultsSettings.getFileFormat()))
                  .getCanonicalPath();
    } catch (final IOException ex) {
      return false;
    }
    PeakResults results;
    switch (resultsSettings.getFileFormat()) {
      case BINARY:
        results = new BinaryFilePeakResults(resultsFilename, source.hasDeviations(),
            source.hasEndFrame(), source.hasId(), resultsSettings.getShowPrecision());
        break;
      case TEXT:
        final TextFilePeakResults f =
            new TextFilePeakResults(resultsFilename, source.hasDeviations(), source.hasEndFrame(),
                source.hasId(), resultsSettings.getShowPrecision());
        f.setDistanceUnit(resultsSettings.getDistanceUnit());
        f.setIntensityUnit(resultsSettings.getIntensityUnit());
        f.setAngleUnit(resultsSettings.getAngleUnit());
        f.setComputePrecision(true);
        results = f;
        break;
      case MALK:
        results = new MalkFilePeakResults(resultsFilename);
        break;
      case TSF:
        results = new TsfPeakResultsWriter(resultsFilename);
        break;
      default:
        throw new IllegalArgumentException(
            "Unsupported file format: " + resultsSettings.getFileFormat());
    }
    results.copySettings(source);
    results.begin();
    results.addAll(source.toArray());
    results.end();
    ImageJUtils.log("Saved %s to %s", source.getName(), resultsFilename);
    return true;
  }
}
