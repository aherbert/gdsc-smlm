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

import uk.ac.sussex.gdsc.core.annotation.Nullable;
import uk.ac.sussex.gdsc.core.data.utils.Rounder;
import uk.ac.sussex.gdsc.core.data.utils.RounderUtils;
import uk.ac.sussex.gdsc.core.ij.BufferedTextWindow;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.logging.Ticker;
import uk.ac.sussex.gdsc.core.match.Coordinate;
import uk.ac.sussex.gdsc.core.match.MatchCalculator;
import uk.ac.sussex.gdsc.core.match.MatchResult;
import uk.ac.sussex.gdsc.core.match.PointPair;
import uk.ac.sussex.gdsc.core.utils.BitFlagUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.utils.ImageRoiPainter;
import uk.ac.sussex.gdsc.smlm.results.ImageSource;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.PeakResultPoint;
import uk.ac.sussex.gdsc.smlm.results.PeakResults;
import uk.ac.sussex.gdsc.smlm.results.TextFilePeakResults;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;

import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.procedure.TIntProcedure;
import gnu.trove.set.hash.TIntHashSet;

import ij.IJ;
import ij.Prefs;
import ij.WindowManager;
import ij.plugin.PlugIn;
import ij.text.TextWindow;

import java.awt.Point;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Consumer;
import java.util.stream.Stream;

/**
 * Compares the coordinates in two sets of results and computes the match statistics.
 */
public class ResultsMatchCalculator implements PlugIn {
  private static final String TITLE = "Results Match Calculator";

  private static AtomicReference<TextWindow> resultsWindowRef = new AtomicReference<>();
  private static AtomicReference<TextWindow> pairsWindowRef = new AtomicReference<>();

  /** The write header flag used in headless mode to ensure the header is only written once. */
  private static AtomicBoolean writeHeader = new AtomicBoolean(true);

  /** Used to draw ROIs on an image when a line is clicked in the results table. */
  private static AtomicReference<ImageRoiPainter> pairPainterRef = new AtomicReference<>();

  /** The rounder for the text output. */
  private static final Rounder rounder = RounderUtils.create(4);

  /** The plugin settings. */
  private Settings settings;

  /**
   * Specify the method used to convert a {@link PeakResult} into coordinates for matching.
   */
  public enum CoordinateMethod {
    /**
     * Create coordinates in each frame from the {@link PeakResult}.
     */
    ALL("All"),
    /**
     * Create coordinates in the first frame from the {@link PeakResult}.
     */
    FIRST_FRAME("First frame"),
    /**
     * Create coordinates in the last frame from the {@link PeakResult}.
     */
    LAST_FRAME("Last frame");

    private static final CoordinateMethod[] values = CoordinateMethod.values();

    /** The description. */
    private final String description;

    private CoordinateMethod(String description) {
      this.description = description;
    }

    /**
     * Gets the description.
     *
     * @return the description
     */
    public String getDescription() {
      return description;
    }

    /**
     * Get the method from the description.
     *
     * @param description the description
     * @return the coordinate method (or null)
     */
    public static CoordinateMethod fromDescription(String description) {
      return fromDescription(description, null);
    }

    /**
     * Get the method from the description.
     *
     * @param description the description
     * @param defaultValue the default value
     * @return the coordinate method
     */
    public static CoordinateMethod fromDescription(String description,
        CoordinateMethod defaultValue) {
      for (final CoordinateMethod value : values) {
        if (value.description.equals(description)) {
          return value;
        }
      }
      return defaultValue;
    }
  }

  /**
   * Lazy load the description for the dialog.
   */
  private static class CoordinateMethodDescriptions {
    /** The descriptions. */
    static final String[] descriptions;

    static {
      descriptions =
          Stream.of(CoordinateMethod.values()).map(m -> m.getDescription()).toArray(String[]::new);
    }
  }

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    // @formatter:off
    private static final String KEY_INPUT_OPTION1 = "gdsc.smlm.resultsmatchcalculator.inputOption1";
    private static final String KEY_INPUT_OPTION2 = "gdsc.smlm.resultsmatchcalculator.inputOption2";
    private static final String KEY_COORD_METHOD1 = "gdsc.smlm.resultsmatchcalculator.coordinateMethod1";
    private static final String KEY_COORD_METHOD2 = "gdsc.smlm.resultsmatchcalculator.coordinateMethod2";
    private static final String KEY_DISTANCE_THRESHOLD = "gdsc.smlm.resultsmatchcalculator.distanceThreshold";
    private static final String KEY_INCREMENTS = "gdsc.smlm.resultsmatchcalculator.increments";
    private static final String KEY_DELTA = "gdsc.smlm.resultsmatchcalculator.delta";
    private static final String KEY_BETA = "gdsc.smlm.resultsmatchcalculator.beta";
    private static final String KEY_SHOW_TABLE = "gdsc.smlm.resultsmatchcalculator.showTable";
    private static final String KEY_SHOW_PAIRS = "gdsc.smlm.resultsmatchcalculator.showPairs";
    private static final String KEY_SAVE_CLASS = "gdsc.smlm.resultsmatchcalculator.saveClassifications";
    private static final String KEY_CLASS_FILE = "gdsc.smlm.resultsmatchcalculator.classificationsFile";
    private static final String KEY_ID_ANALYSIS = "gdsc.smlm.resultsmatchcalculator.idAnalysis";
    private static final String KEY_SAVE_PAIRS = "gdsc.smlm.resultsmatchcalculator.savePairs";
    private static final String KEY_PAIRS_DIR = "gdsc.smlm.resultsmatchcalculator.pairsDirectory";
    private static final String KEY_OUTPUT_END_FRAME = "gdsc.smlm.resultsmatchcalculator.outputEndFrame";
    // @formatter:on

    /** The options to save the classifications. */
    private static final String[] saveClassificationsOptions =
        {"None", "Matched", "Unmatched", "All"};
    private static final int SAVE_MATCHED = 0x01;
    private static final int SAVE_UNMATCHED = 0x02;

    String inputOption1 = "";
    String inputOption2 = "";
    CoordinateMethod coordinateMethod1 = CoordinateMethod.ALL;
    CoordinateMethod coordinateMethod2 = CoordinateMethod.ALL;
    double distanceThreshold = 0.5;
    int increments = 5;
    double delta = 0.1;
    double beta = 4;
    boolean showTable = true;
    boolean showPairs;
    int saveClassificationsOption;
    String classificationsFile = "";
    boolean idAnalysis;
    boolean savePairs;
    String pairsDirectory = "";
    boolean outputEndFrame;

    Settings() {
      inputOption1 = Prefs.get(KEY_INPUT_OPTION1, "");
      inputOption2 = Prefs.get(KEY_INPUT_OPTION2, "");
      coordinateMethod1 = CoordinateMethod
          .fromDescription(Prefs.get(KEY_COORD_METHOD1, CoordinateMethod.ALL.getDescription()));
      coordinateMethod2 = CoordinateMethod
          .fromDescription(Prefs.get(KEY_COORD_METHOD2, CoordinateMethod.ALL.getDescription()));
      distanceThreshold = Prefs.get(KEY_DISTANCE_THRESHOLD, 0.5);
      increments = Prefs.getInt(KEY_INCREMENTS, 5);
      delta = Prefs.get(KEY_DELTA, 0.1);
      beta = Prefs.get(KEY_BETA, 4);
      showTable = Prefs.get(KEY_SHOW_TABLE, true);
      showPairs = Prefs.get(KEY_SHOW_PAIRS, false);
      saveClassificationsOption = Prefs.getInt(KEY_SAVE_CLASS, 0);
      classificationsFile = Prefs.get(KEY_CLASS_FILE, "");
      idAnalysis = Prefs.get(KEY_ID_ANALYSIS, false);
      savePairs = Prefs.get(KEY_SAVE_PAIRS, false);
      pairsDirectory = Prefs.get(KEY_PAIRS_DIR, "");
      outputEndFrame = Prefs.get(KEY_OUTPUT_END_FRAME, false);
    }

    Settings(Settings source) {
      this.inputOption1 = source.inputOption1;
      this.inputOption2 = source.inputOption2;
      this.coordinateMethod1 = source.coordinateMethod1;
      this.coordinateMethod2 = source.coordinateMethod2;
      this.distanceThreshold = source.distanceThreshold;
      this.increments = source.increments;
      this.delta = source.delta;
      this.beta = source.beta;
      this.showTable = source.showTable;
      this.showPairs = source.showPairs;
      this.saveClassificationsOption = source.saveClassificationsOption;
      this.classificationsFile = source.classificationsFile;
      this.idAnalysis = source.idAnalysis;
      this.savePairs = source.savePairs;
      this.pairsDirectory = source.pairsDirectory;
      this.outputEndFrame = source.outputEndFrame;
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
      Prefs.set(KEY_INPUT_OPTION1, inputOption1);
      Prefs.set(KEY_INPUT_OPTION2, inputOption2);
      Prefs.set(KEY_COORD_METHOD1, coordinateMethod1.getDescription());
      Prefs.set(KEY_COORD_METHOD2, coordinateMethod2.getDescription());
      Prefs.set(KEY_DISTANCE_THRESHOLD, distanceThreshold);
      Prefs.set(KEY_INCREMENTS, increments);
      Prefs.set(KEY_DELTA, delta);
      Prefs.set(KEY_BETA, beta);
      Prefs.set(KEY_SHOW_TABLE, showTable);
      Prefs.set(KEY_SHOW_PAIRS, showPairs);
      Prefs.set(KEY_SAVE_CLASS, saveClassificationsOption);
      Prefs.set(KEY_CLASS_FILE, classificationsFile);
      Prefs.set(KEY_ID_ANALYSIS, idAnalysis);
      Prefs.set(KEY_SAVE_PAIRS, savePairs);
      Prefs.set(KEY_PAIRS_DIR, pairsDirectory);
      Prefs.set(KEY_OUTPUT_END_FRAME, outputEndFrame);
    }

    /**
     * Check if the save matched flag is set.
     *
     * @return true if set
     */
    boolean isSaveMatched() {
      return BitFlagUtils.anySet(saveClassificationsOption, SAVE_MATCHED);
    }

    /**
     * Check if the save unmatched flag is set.
     *
     * @return true if set
     */
    boolean isSaveUnmatched() {
      return BitFlagUtils.anySet(saveClassificationsOption, SAVE_UNMATCHED);
    }

    /**
     * Check if any save classifications options are set.
     *
     * @return true if set
     */
    boolean isSaveClassifications() {
      return BitFlagUtils.anySet(saveClassificationsOption, SAVE_MATCHED | SAVE_UNMATCHED);
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "No localisations in memory");
      return;
    }

    if (!showDialog()) {
      return;
    }

    // Load the results
    final MemoryPeakResults results1 =
        ResultsManager.loadInputResults(settings.inputOption1, false, null, null);
    final MemoryPeakResults results2 =
        ResultsManager.loadInputResults(settings.inputOption2, false, null, null);
    IJ.showStatus("");
    if (results1 == null || results1.size() == 0) {
      IJ.error(TITLE, "No results 1 could be loaded");
      return;
    }
    if (results2 == null || results2.size() == 0) {
      IJ.error(TITLE, "No results 2 could be loaded");
      return;
    }
    if (results1.getDistanceUnit() != results2.getDistanceUnit()) {
      IJ.error(TITLE, "Distance unit should be the same for the results");
      return;
    }

    final long start = System.nanoTime();
    runCompareCoordinates(results1, results2);
    final long nanoseconds = System.nanoTime() - start;

    IJ.showStatus(String.format("%s = %s", TITLE, TextUtils.nanosToString(nanoseconds)));
  }

  private boolean showDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

    settings = Settings.load();
    gd.addMessage("Compare the points in two results sets\nand compute the match statistics");
    ResultsManager.addInput(gd, "Results1", settings.inputOption1, InputSource.MEMORY);
    ResultsManager.addInput(gd, "Results2", settings.inputOption2, InputSource.MEMORY);
    gd.addChoice("Coordinate_method1", CoordinateMethodDescriptions.descriptions,
        settings.coordinateMethod1.getDescription());
    gd.addChoice("Coordinate_method2", CoordinateMethodDescriptions.descriptions,
        settings.coordinateMethod2.getDescription());
    gd.addNumericField("Distance", settings.distanceThreshold, 2);

    gd.addSlider("Increments", 0, 10, settings.increments);
    gd.addNumericField("Delta", settings.delta, 2);
    gd.addNumericField("Beta", settings.beta, 2);
    gd.addCheckbox("Show_table", settings.showTable);
    gd.addCheckbox("Show_pairs", settings.showPairs);
    gd.addChoice("Save_classifications", Settings.saveClassificationsOptions,
        settings.saveClassificationsOption);
    gd.addCheckbox("Id_analysis", settings.idAnalysis);
    gd.addCheckbox("Save_pairs", settings.savePairs);
    gd.addCheckbox("Output_end_frame", settings.outputEndFrame);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.inputOption1 = gd.getNextChoice();
    settings.inputOption2 = gd.getNextChoice();
    settings.coordinateMethod1 =
        CoordinateMethod.fromDescription(gd.getNextChoice(), CoordinateMethod.ALL);
    settings.coordinateMethod2 =
        CoordinateMethod.fromDescription(gd.getNextChoice(), CoordinateMethod.ALL);
    settings.distanceThreshold = gd.getNextNumber();
    settings.increments = (int) gd.getNextNumber();
    settings.delta = gd.getNextNumber();
    settings.beta = gd.getNextNumber();
    settings.showTable = gd.getNextBoolean();
    settings.showPairs = gd.getNextBoolean();
    settings.saveClassificationsOption = gd.getNextChoiceIndex();
    settings.idAnalysis = gd.getNextBoolean();
    settings.savePairs = gd.getNextBoolean();
    settings.outputEndFrame = gd.getNextBoolean();
    settings.save();

    if (!(settings.showTable || settings.showPairs || !settings.isSaveClassifications())) {
      IJ.error(TITLE, "No outputs specified");
      return false;
    }

    // Check arguments
    try {
      ParameterUtils.isPositive("Distance threshold", settings.distanceThreshold);
      ParameterUtils.isPositive("Increments", settings.increments);
      ParameterUtils.isAboveZero("Delta", settings.delta);
      ParameterUtils.isPositive("Beta", settings.beta);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return true;
  }

  @SuppressWarnings("null")
  private void runCompareCoordinates(MemoryPeakResults results1, MemoryPeakResults results2) {
    final boolean requirePairs = settings.showPairs || settings.isSaveClassifications();
    final boolean saveMatched = settings.isSaveMatched();
    final boolean saveUnmatched = settings.isSaveUnmatched();

    final TextFilePeakResults fileResults = createFilePeakResults(results2);

    final List<PointPair> allMatches = new LinkedList<>();
    final List<PointPair> pairs = (requirePairs) ? new LinkedList<>() : null;

    final double maxDistance = settings.distanceThreshold + settings.increments * settings.delta;

    // Divide the results into time points
    final TIntObjectHashMap<ArrayList<Coordinate>> actualCoordinates =
        getCoordinates(results1, settings.coordinateMethod1);
    final TIntObjectHashMap<ArrayList<Coordinate>> predictedCoordinates =
        getCoordinates(results2, settings.coordinateMethod2);

    int n1 = 0;
    int n2 = 0;

    // Process each time point
    for (final Integer t : getTimepoints(actualCoordinates, predictedCoordinates)) {
      final Coordinate[] actual = getCoordinates(actualCoordinates, t);
      final Coordinate[] predicted = getCoordinates(predictedCoordinates, t);

      final List<Coordinate> tp = null;
      List<Coordinate> fp = null;
      List<Coordinate> fn = null;
      final List<PointPair> matches = new LinkedList<>();
      if (requirePairs) {
        fp = new LinkedList<>();
        fn = new LinkedList<>();
      }

      MatchCalculator.analyseResults2D(actual, predicted, maxDistance, tp, fp, fn, matches);

      // Aggregate
      n1 += actual.length;
      n2 += predicted.length;

      allMatches.addAll(matches);
      if (settings.showPairs) {
        pairs.addAll(matches);
        for (final Coordinate c : fn) {
          pairs.add(new PointPair(c, null));
        }
        for (final Coordinate c : fp) {
          pairs.add(new PointPair(null, c));
        }
      }
      if (fileResults != null) {
        // Matches are marked in the original value with 1 for true, 0 for false
        if (saveMatched) {
          for (final PointPair pair : matches) {
            PeakResult result = ((PeakResultPoint) pair.getPoint2()).getPeakResult();
            result = result.copy();
            result.setOrigValue(1);
            fileResults.add(result);
          }
        }
        if (saveUnmatched) {
          for (final Coordinate c : fp) {
            PeakResult result = ((PeakResultPoint) c).getPeakResult();
            result = result.copy();
            result.setOrigValue(0);
            fileResults.add(result);
          }
        }
      }
    }

    if (fileResults != null) {
      fileResults.end();
    }

    final boolean doIdAnalysis1 = settings.idAnalysis && haveIds(results1);
    final boolean doIdAnalysis2 = settings.idAnalysis && haveIds(results2);

    // Create output.
    // This supports headless mode with just the results table
    // or graphical mode with a results table and pairs window.
    final boolean headless = java.awt.GraphicsEnvironment.isHeadless();

    TextWindow resultsWindow = null;
    if (!headless) {
      resultsWindow =
          (settings.showTable) ? createResultsWindow(doIdAnalysis1 || doIdAnalysis2) : null;
      showPairs(results1, pairs, resultsWindow);
    }

    showResults(results1, results2, allMatches, n1, n2, doIdAnalysis1, doIdAnalysis2,
        resultsWindow);

    savePairs(results1, results2, allMatches);
  }

  /**
   * Show the match pairs in a results table.
   *
   * <p>Adds an ROI painter for the original image source of results set 1 (if it is visible) to the
   * table.
   *
   * @param results1 the first set of results
   * @param pairs the pairs
   * @param resultsWindow the results window
   */
  private void showPairs(MemoryPeakResults results1, final List<PointPair> pairs,
      final TextWindow resultsWindow) {
    if (!settings.showPairs) {
      return;
    }
    final TextWindow pairsWindow = createPairsWindow(resultsWindow);
    pairsWindow.getTextPanel().clear();
    createPairPainter(pairsWindow, results1.getSource());
    final Ticker ticker = ImageJUtils.createTicker(pairs.size(), 0, "Writing pairs table");
    try (final BufferedTextWindow bw = new BufferedTextWindow(pairsWindow)) {
      final StringBuilder sb = new StringBuilder();
      for (final PointPair pair : pairs) {
        bw.append(addPairResult(sb, pair));
        ticker.tick();
      }
    }
    ImageJUtils.finished();
  }

  @SuppressWarnings("null")
  private void showResults(MemoryPeakResults results1, MemoryPeakResults results2,
      final List<PointPair> allMatches, int n1, int n2, final boolean doIdAnalysis1,
      final boolean doIdAnalysis2, TextWindow resultsWindow) {
    if (!settings.showTable) {
      return;
    }

    // Output the results
    Consumer<String> output;
    if (resultsWindow != null) {
      output = resultsWindow::append;
    } else {
      // Headless mode
      output = IJ::log;
      if (writeHeader.get()) {
        writeHeader.set(false);
        IJ.log(createResultsHeader(settings.idAnalysis));
      }
    }

    // We have the results for the largest distance.
    // Now reduce the distance threshold and recalculate the results
    final double[] distanceThresholds =
        getDistances(settings.distanceThreshold, settings.increments, settings.delta);
    final double[] pairDistances = getPairDistances(allMatches);
    // Re-use storage for the ID analysis
    TIntHashSet id1 = null;
    TIntHashSet id2 = null;
    TIntHashSet matchId1 = null;
    TIntHashSet matchId2 = null;
    final boolean doIdAnalysis = doIdAnalysis1 || doIdAnalysis2;
    if (doIdAnalysis) {
      if (doIdAnalysis1) {
        id1 = getIds(results1);
        matchId1 = new TIntHashSet(id1.size());
      }
      if (doIdAnalysis2) {
        id2 = getIds(results2);
        matchId2 = new TIntHashSet(id2.size());
      }
    }
    final StringBuilder sb = new StringBuilder();
    for (final double distanceThreshold : distanceThresholds) {
      double rms = 0;
      int tp2 = 0;
      final double d2 = distanceThreshold * distanceThreshold;
      for (final double d : pairDistances) {
        if (d <= d2) {
          rms += d;
          tp2++;
        }
      }
      // All non-true positives must be added to the false totals.
      final int fp2 = n2 - tp2;
      final int fn2 = n1 - tp2;

      final MatchResult result =
          new MatchResult(tp2, fp2, fn2, (tp2 > 0) ? Math.sqrt(rms / tp2) : 0);

      MatchResult idResult1 = null;
      MatchResult idResult2 = null;
      if (doIdAnalysis) {
        if (doIdAnalysis1) {
          matchId1.clear();
        }
        if (doIdAnalysis2) {
          matchId2.clear();
        }
        int index = 0;
        for (final PointPair pair : allMatches) {
          if (pairDistances[index++] <= d2) {
            if (doIdAnalysis1) {
              matchId1.add(((PeakResultPoint) pair.getPoint1()).getPeakResult().getId());
            }
            if (doIdAnalysis2) {
              matchId2.add(((PeakResultPoint) pair.getPoint2()).getPeakResult().getId());
            }
          }
        }
        // Only the actual points are checked for Ids. For example these could be from the
        // Create Data plugin with actual fluorophore Ids.
        // => Only the recall will be valid: tp / (tp + fn)
        if (doIdAnalysis1) {
          idResult1 = new MatchResult(matchId1.size(), 0, id1.size() - matchId1.size(), 0);
        }
        if (doIdAnalysis2) {
          idResult2 = new MatchResult(matchId2.size(), 0, id2.size() - matchId2.size(), 0);
        }
      }

      addResult(sb, settings.inputOption1, settings.inputOption2, distanceThreshold, result,
          idResult1, idResult2);

      output.accept(sb.toString());
    }
  }

  private static boolean haveIds(MemoryPeakResults results) {
    return results.hasId();
  }

  private TextFilePeakResults createFilePeakResults(MemoryPeakResults results) {
    if (!settings.isSaveClassifications()) {
      return null;
    }
    final String filename =
        ImageJUtils.getFilename("Classifications_File", settings.classificationsFile);
    if (filename != null) {
      settings.classificationsFile = filename;
      final TextFilePeakResults r =
          new TextFilePeakResults(filename, false, outputEndFrame(results));
      r.copySettings(results);
      r.begin();
      return r;
    }
    return null;
  }

  private TextFilePeakResults createFilePeakResults(String directory, int set,
      MemoryPeakResults results, double distanceLow, double distanceHigh) {
    final String filename = directory + String.format("Match%d_%s_%s_%s.txt", set,
        results.getName(), rounder.toString(distanceLow), rounder.toString(distanceHigh));
    final TextFilePeakResults r = new TextFilePeakResults(filename, false, outputEndFrame(results));
    r.copySettings(results);
    r.begin();
    return r;
  }

  /**
   * Checks if the end frame is required in the file output.
   *
   * @param results the results
   * @return true, if successful
   */
  private boolean outputEndFrame(MemoryPeakResults results) {
    return settings.outputEndFrame || results.hasEndFrame();
  }

  /**
   * Build a map between the peak id (time point) and a list of coordinates.
   *
   * @param results the results
   * @return the coordinates
   */
  public static TIntObjectHashMap<ArrayList<Coordinate>> getCoordinates(MemoryPeakResults results) {
    return getCoordinates(results, CoordinateMethod.ALL, false);
  }

  /**
   * Build a map between the peak id (time point) and a list of coordinates.
   *
   * @param results the results
   * @param coordinateMethod the coordinate method
   * @return the coordinates
   */
  public static TIntObjectHashMap<ArrayList<Coordinate>> getCoordinates(MemoryPeakResults results,
      CoordinateMethod coordinateMethod) {
    return getCoordinates(results, coordinateMethod, false);
  }


  /**
   * Build a map between the peak id (time point) and a list of coordinates.
   *
   * @param results the results
   * @param integerCoordinates True if the values should be rounded down to integers
   * @return the coordinates
   */
  public static TIntObjectHashMap<ArrayList<Coordinate>> getCoordinates(MemoryPeakResults results,
      final boolean integerCoordinates) {
    return getCoordinates(results, CoordinateMethod.ALL, integerCoordinates);
  }

  /**
   * Build a map between the peak id (time point) and a list of coordinates.
   *
   * @param results the results
   * @param coordinateMethod the coordinate method
   * @param integerCoordinates True if the values should be rounded down to integers
   * @return the coordinates
   */
  public static TIntObjectHashMap<ArrayList<Coordinate>> getCoordinates(MemoryPeakResults results,
      CoordinateMethod coordinateMethod, final boolean integerCoordinates) {
    final TIntObjectHashMap<ArrayList<Coordinate>> coords = new TIntObjectHashMap<>();
    if (results.size() > 0) {
      // Do not use HashMap directly to build the coords object since there
      // will be many calls to getEntry(). Instead sort the results and use
      // a new list for each time point
      results.sort();
      final int minT = results.getFirstFrame();
      final int maxT = results.getLastFrame();

      // Create lists
      final ArrayList<ArrayList<Coordinate>> tmpCoords = new ArrayList<>(maxT - minT + 1);
      for (int t = minT; t <= maxT; t++) {
        tmpCoords.add(new ArrayList<Coordinate>());
      }

      // Add the results to the lists
      results.forEach((PeakResultProcedure) result -> {
        final float x;
        final float y;
        final float z;
        if (integerCoordinates) {
          x = (int) result.getXPosition();
          y = (int) result.getYPosition();
          z = (int) result.getZPosition();
        } else {
          x = result.getXPosition();
          y = result.getYPosition();
          z = result.getZPosition();
        }

        final int startFrame = getStartFrame(result, coordinateMethod);
        final int endFrame = getEndFrame(result, coordinateMethod);
        for (int t = startFrame - minT, i = endFrame - startFrame + 1; i-- > 0; t++) {
          tmpCoords.get(t).add(new PeakResultPoint(t + minT, x, y, z, result));
        }
      });

      // Put in the map
      for (int t = minT, i = 0; t <= maxT; t++, i++) {
        coords.put(t, tmpCoords.get(i));
      }
    }
    return coords;
  }

  /**
   * Return an array of coordinates for the given time point. Returns an empty array if there are no
   * coordinates.
   *
   * @param coords the coords
   * @param time the time
   * @return the coordinates
   */
  public static Coordinate[] getCoordinates(TIntObjectHashMap<ArrayList<Coordinate>> coords,
      int time) {
    final ArrayList<Coordinate> tmp = coords.get(time);
    if (tmp != null) {
      return tmp.toArray(new Coordinate[tmp.size()]);
    }
    return new Coordinate[0];
  }

  private static int getStartFrame(PeakResult result, CoordinateMethod coordinateMethod) {
    if (coordinateMethod == CoordinateMethod.LAST_FRAME) {
      return result.getEndFrame();
    }
    return result.getFrame();
  }

  private static int getEndFrame(PeakResult result, CoordinateMethod coordinateMethod) {
    if (coordinateMethod == CoordinateMethod.FIRST_FRAME) {
      return result.getFrame();
    }
    return result.getEndFrame();
  }

  /**
   * Merge the time points from each map into a single sorted list of unique time points.
   *
   * @param actualCoordinates the actual coordinates
   * @param predictedCoordinates the predicted coordinates
   * @return a list of time points
   */
  private static int[] getTimepoints(TIntObjectHashMap<ArrayList<Coordinate>> actualCoordinates,
      TIntObjectHashMap<ArrayList<Coordinate>> predictedCoordinates) {

    // Do inline to avoid materialising the keys arrays
    final TIntHashSet hashset =
        new TIntHashSet(Math.max(actualCoordinates.size(), predictedCoordinates.size()));
    final TIntProcedure p = value -> {
      hashset.add(value);
      return true;
    };
    actualCoordinates.forEachKey(p);
    predictedCoordinates.forEachKey(p);
    final int[] set = hashset.toArray();

    Arrays.sort(set);
    return set;
  }

  private static TextWindow createResultsWindow(boolean doIdAnalysis) {
    final String header = createResultsHeader(doIdAnalysis);
    final TextWindow resultsWindow = ImageJUtils.refresh(resultsWindowRef,
        () -> new TextWindow(TITLE + " Results", header, "", 900, 300));
    ImageJUtils.refreshHeadings(resultsWindow, header, true);
    return resultsWindow;
  }

  private static TextWindow createPairsWindow(TextWindow resultsWindow) {
    return ImageJUtils.refresh(pairsWindowRef, () -> {
      final TextWindow window = new TextWindow(TITLE + " Pairs", createPairsHeader(), "", 900, 300);
      // Position relative to results window
      if (resultsWindow != null) {
        final Point p = resultsWindow.getLocation();
        p.y += resultsWindow.getHeight();
        window.setLocation(p);
      }
      return window;
    });
  }

  private static String createResultsHeader(boolean idAnalysis) {
    final StringBuilder sb = new StringBuilder();
    sb.append("Image 1\t");
    sb.append("Image 2\t");
    sb.append("Distance (px)\t");
    sb.append("N\t");
    sb.append("TP\t");
    sb.append("FP\t");
    sb.append("FN\t");
    sb.append("Jaccard\t");
    sb.append("RMSD\t");
    sb.append("Precision\t");
    sb.append("Recall\t");
    sb.append("F0.5\t");
    sb.append("F1\t");
    sb.append("F2\t");
    sb.append("F-beta");
    if (idAnalysis) {
      sb.append("\tId1-N");
      sb.append("\tId1-TP");
      sb.append("\tId1-Recall");
      sb.append("\tId2-N");
      sb.append("\tId2-TP");
      sb.append("\tId2-Recall");
    }
    return sb.toString();
  }

  private void addResult(StringBuilder sb, String i1, String i2, double distanceThrehsold,
      MatchResult result, MatchResult idResult1, MatchResult idResult2) {
    sb.setLength(0);
    sb.append(i1).append('\t');
    sb.append(i2).append('\t');
    sb.append(rounder.round(distanceThrehsold)).append('\t');
    sb.append(result.getNumberPredicted()).append('\t');
    sb.append(result.getTruePositives()).append('\t');
    sb.append(result.getFalsePositives()).append('\t');
    sb.append(result.getFalseNegatives()).append('\t');
    sb.append(rounder.round(result.getJaccard())).append('\t');
    sb.append(rounder.round(result.getRmsd())).append('\t');
    sb.append(rounder.round(result.getPrecision())).append('\t');
    sb.append(rounder.round(result.getRecall())).append('\t');
    sb.append(rounder.round(result.getFScore(0.5))).append('\t');
    sb.append(rounder.round(result.getFScore(1.0))).append('\t');
    sb.append(rounder.round(result.getFScore(2.0))).append('\t');
    sb.append(rounder.round(result.getFScore(settings.beta)));
    if (idResult1 != null) {
      sb.append('\t').append(idResult1.getNumberPredicted());
      sb.append('\t').append(idResult1.getTruePositives());
      sb.append('\t').append(rounder.round(idResult1.getRecall()));
    } else if (idResult2 != null) {
      sb.append("\t-\t-\t-");
    }
    if (idResult2 != null) {
      sb.append('\t').append(idResult2.getNumberPredicted());
      sb.append('\t').append(idResult2.getTruePositives());
      sb.append('\t').append(rounder.round(idResult2.getRecall()));
    } else if (idResult1 != null) {
      sb.append("\t-\t-\t-");
    }
  }

  private static String createPairsHeader() {
    final StringBuilder sb = new StringBuilder();
    sb.append("T\t");
    sb.append("X1\t");
    sb.append("Y1\t");
    sb.append("Z1\t");
    sb.append("X2\t");
    sb.append("Y2\t");
    sb.append("Z2\t");
    sb.append("Distance\t");
    return sb.toString();
  }

  private static String addPairResult(StringBuilder sb, PointPair pair) {
    sb.setLength(0);
    final PeakResultPoint p1 = (PeakResultPoint) pair.getPoint1();
    final PeakResultPoint p2 = (PeakResultPoint) pair.getPoint2();
    final int t = (p1 != null) ? p1.getTime() : p2.getTime();
    sb.append(t).append('\t');
    addPoint(sb, p1);
    addPoint(sb, p2);
    final double d = pair.getXyDistance();
    if (d >= 0) {
      sb.append(rounder.round(d)).append('\t');
    } else {
      sb.append("-\t");
    }
    return sb.toString();
  }

  private static void addPoint(StringBuilder sb, PeakResultPoint result) {
    if (result == null) {
      sb.append("-\t-\t-\t");
    } else {
      sb.append(rounder.round(result.getX())).append('\t');
      sb.append(rounder.round(result.getY())).append('\t');
      sb.append(rounder.round(result.getZ())).append('\t');
    }
  }

  private static void createPairPainter(TextWindow pairsWindow, ImageSource source) {
    // Find if the image is open
    if (source != null && TextUtils.isNotEmpty(source.getOriginal().getName())) {
      final String title = source.getOriginal().getName();

      if (WindowManager.getImage(title) == null) {
        return;
      }

      // Create the pair painter if necessary
      ImageRoiPainter painter = pairPainterRef.get();
      if (painter == null) {
        painter = new ImageRoiPainter(null, null, ResultsMatchCalculator::getRoiCoordinates);
        pairPainterRef.set(painter);
      }

      // Update the painter
      painter.setTextPanel(pairsWindow.getTextPanel());
      painter.setTitle(title);
    }
  }

  @Nullable
  private static double[] getRoiCoordinates(String line) {
    // Extract the startT and x,y coordinates from the first available point in the pair
    final int[] index = {1, 4};
    final String[] fields = line.split("\t");
    final int startT = Integer.parseInt(fields[0]);
    for (final int i : index) {
      if (i < fields.length) {
        if (fields[i].equals("-")) {
          continue;
        }
        final double x = Double.parseDouble(fields[i]);
        final double y = Double.parseDouble(fields[i + 1]);
        return new double[] {startT, x, y};
      }
    }
    return null;
  }

  private static TIntHashSet getIds(MemoryPeakResults results) {
    final TIntHashSet ids = new TIntHashSet(results.size());
    results.forEach((PeakResultProcedure) result -> ids.add(result.getId()));
    return ids;
  }

  private static double[] getDistances(double distanceThreshold, int increments, double delta) {
    final double[] d = new double[increments + 1];
    for (int i = 0; i <= increments; i++) {
      d[i] = distanceThreshold + i * delta;
    }
    return d;
  }

  private static double[] getPairDistances(List<PointPair> pairs) {
    final double[] d = new double[pairs.size()];
    int index = 0;
    for (final PointPair pair : pairs) {
      d[index++] = pair.getXyDistanceSquared();
    }
    return d;
  }

  private void savePairs(MemoryPeakResults results1, MemoryPeakResults results2,
      List<PointPair> allMatches) {
    if (!settings.savePairs) {
      return;
    }

    // Get the directory
    final String directory = ImageJUtils.getDirectory("Pairs_directory", settings.pairsDirectory);
    if (directory == null) {
      return;
    }
    settings.pairsDirectory = directory;

    final double[] distanceThresholds =
        getDistances(settings.distanceThreshold, settings.increments, settings.delta);
    // Create output files for each distance band
    final PeakResults[] output1 = new PeakResults[distanceThresholds.length];
    final PeakResults[] output2 = new PeakResults[distanceThresholds.length];
    double high = 0;
    for (int i = 0; i < distanceThresholds.length; i++) {
      final double low = high;
      high = distanceThresholds[i];
      output1[i] = createFilePeakResults(directory, 1, results1, low, high);
      output2[i] = createFilePeakResults(directory, 2, results2, low, high);
    }

    // Square the thresholds
    SimpleArrayUtils.apply(distanceThresholds, v -> v * v);
    final double[] pairDistances = getPairDistances(allMatches);
    int index = 0;
    for (final PointPair pair : allMatches) {
      final int insert = search(distanceThresholds, pairDistances[index++]);
      if (insert != -1) {
        final PeakResult r1 = ((PeakResultPoint) pair.getPoint1()).getPeakResult();
        final PeakResult r2 = ((PeakResultPoint) pair.getPoint2()).getPeakResult();
        output1[insert].add(r1);
        output2[insert].add(r2);
      }
    }

    for (int i = 0; i < output1.length; i++) {
      output1[i].end();
      output2[i].end();
    }
  }

  private static int search(double[] distanceThresholds, double distanceSquared) {
    for (int i = 0; i < distanceThresholds.length; i++) {
      if (distanceSquared <= distanceThresholds[i]) {
        return i;
      }
    }
    // This should not happen since all pairs should be under the max distance threshold
    return -1;
  }

  /**
   * Compare the coordinates in two results sets.
   *
   * @param results1 the results 1
   * @param results2 the results 2
   * @param distance the distance
   * @return the match result
   */
  public static MatchResult compareCoordinates(MemoryPeakResults results1,
      MemoryPeakResults results2, double distance) {
    // Divide the results into time points
    final TIntObjectHashMap<ArrayList<Coordinate>> actualCoordinates = getCoordinates(results1);
    final TIntObjectHashMap<ArrayList<Coordinate>> predictedCoordinates = getCoordinates(results2);

    return compareCoordinates(actualCoordinates, predictedCoordinates, distance);
  }

  /**
   * Compare the coordinates on a frame-by-frame basis.
   *
   * @param actualCoordinates the actual coordinates
   * @param predictedCoordinates the predicted coordinates
   * @param distance the distance
   * @return the match result
   */
  public static MatchResult compareCoordinates(
      TIntObjectHashMap<ArrayList<Coordinate>> actualCoordinates,
      TIntObjectHashMap<ArrayList<Coordinate>> predictedCoordinates, double distance) {
    int tp = 0;
    int fp = 0;
    int fn = 0;

    // Process each time point
    for (final Integer t : getTimepoints(actualCoordinates, predictedCoordinates)) {
      final Coordinate[] actual = getCoordinates(actualCoordinates, t);
      final Coordinate[] predicted = getCoordinates(predictedCoordinates, t);

      final MatchResult r = MatchCalculator.analyseResults2D(actual, predicted, distance);

      // Aggregate
      tp += r.getTruePositives();
      fp += r.getFalsePositives();
      fn += r.getFalseNegatives();
    }

    return new MatchResult(tp, fp, fn, 0);
  }
}
