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

package uk.ac.sussex.gdsc.smlm.ij.plugins;

import ij.IJ;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import ij.text.TextWindow;
import java.awt.Color;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Consumer;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.match.ClassificationResult;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.GUIFilterSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResultsReader;
import uk.ac.sussex.gdsc.smlm.results.count.Counter;
import uk.ac.sussex.gdsc.smlm.results.filter.AndFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.Filter;
import uk.ac.sussex.gdsc.smlm.results.filter.FilterSet;
import uk.ac.sussex.gdsc.smlm.results.filter.FilterXStreamUtils;
import uk.ac.sussex.gdsc.smlm.results.filter.OrFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.PrecisionFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.PrecisionHysteresisFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.SnrFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.SnrHysteresisFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.TraceFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.WidthFilter;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;

/**
 * Run different filtering methods on a set of labelled peak results outputting performance
 * statistics on the success of the filter.
 *
 * <p>All results files in a specified directory are read. If the peak result original value is set
 * to 1 it is considered a true peak, 0 for a false peak. Filtering is done using e.g. SNR
 * threshold, Precision thresholds, etc. The statistics reported are shown in a table, e.g.
 * precision, Jaccard, F-score.
 */
public class FilterAnalysis implements PlugIn {
  private static final String TITLE = "Filter Analysis";
  private static AtomicReference<TextWindow> resultsWindowRef = new AtomicReference<>();
  private static AtomicReference<TextWindow> sensitivityWindowRef = new AtomicReference<>();

  private static final AtomicReference<LastResults> lastResults = new AtomicReference<>();

  private ArrayList<NamedPlot> plots;
  private HashMap<String, FilterScore> bestFilter;
  private LinkedList<String> bestFilterOrder;

  private final boolean isHeadless;

  private static class LastResults {
    List<MemoryPeakResults> resultsList;
    String inputDirectory;

    LastResults(List<MemoryPeakResults> resultsList, String inputDirectory) {
      this.resultsList = resultsList;
      this.inputDirectory = inputDirectory;
    }
  }

  private static class FilterScore {
    Filter filter;
    double score;

    FilterScore(Filter filter, double score) {
      update(filter, score);
    }

    void update(Filter filter, double score) {
      this.filter = filter;
      this.score = score;
    }
  }

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    boolean saveFilterSets;
    boolean showResultsTable;
    int plotTopN;
    boolean calculateSensitivity;
    double delta;

    boolean snrFilter;
    int minSnr;
    int maxSnr;
    double minWidth;
    double maxWidth;
    double incWidth;

    boolean precisionFilter;
    int minPrecision;
    int maxPrecision;

    boolean traceFilter;
    double minDistance;
    double maxDistance;
    double incDistance;
    int minTime;
    int maxTime;
    int incTime;

    boolean hysteresisSnrFilter;
    int minSnrGap;
    int maxSnrGap;
    int incSnrGap;

    boolean hysteresisPrecisionFilter;
    int minPrecisionGap;
    int maxPrecisionGap;
    int incPrecisionGap;

    Settings() {
      // Set defaults
      showResultsTable = true;
      delta = 0.1;
      snrFilter = true;
      minSnr = 20;
      maxSnr = 80;
      minWidth = 1.5;
      maxWidth = 2.0;
      incWidth = 0.5;
      precisionFilter = true;
      minPrecision = 20;
      maxPrecision = 80;
      minDistance = 0.3;
      maxDistance = 1.2;
      incDistance = 0.3;
      minTime = 1;
      maxTime = 80;
      incTime = 10;
      hysteresisSnrFilter = true;
      minSnrGap = 10;
      maxSnrGap = 40;
      incSnrGap = 10;
      hysteresisPrecisionFilter = true;
      minPrecisionGap = 10;
      maxPrecisionGap = 40;
      incPrecisionGap = 10;
    }

    Settings(Settings source) {
      saveFilterSets = source.saveFilterSets;
      showResultsTable = source.showResultsTable;
      plotTopN = source.plotTopN;
      calculateSensitivity = source.calculateSensitivity;
      delta = source.delta;
      snrFilter = source.snrFilter;
      minSnr = source.minSnr;
      maxSnr = source.maxSnr;
      minWidth = source.minWidth;
      maxWidth = source.maxWidth;
      incWidth = source.incWidth;
      precisionFilter = source.precisionFilter;
      minPrecision = source.minPrecision;
      maxPrecision = source.maxPrecision;
      traceFilter = source.traceFilter;
      minDistance = source.minDistance;
      maxDistance = source.maxDistance;
      incDistance = source.incDistance;
      minTime = source.minTime;
      maxTime = source.maxTime;
      incTime = source.incTime;
      hysteresisSnrFilter = source.hysteresisSnrFilter;
      minSnrGap = source.minSnrGap;
      maxSnrGap = source.maxSnrGap;
      incSnrGap = source.incSnrGap;
      hysteresisPrecisionFilter = source.hysteresisPrecisionFilter;
      minPrecisionGap = source.minPrecisionGap;
      maxPrecisionGap = source.maxPrecisionGap;
      incPrecisionGap = source.incPrecisionGap;
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

  private static class NamedPlot {
    String name;
    String xAxisName;
    double[] xValues;
    double[] yValues;
    double score;

    NamedPlot(String name, String xAxisName, double[] xValues, double[] yValues) {
      this.name = name;
      updateValues(xAxisName, xValues, yValues);
    }

    void updateValues(String xAxisName, double[] xValues, double[] yValues) {
      this.xAxisName = xAxisName;
      this.xValues = xValues;
      this.yValues = yValues;
      this.score = MathUtils.max(yValues);
    }

    static int compare(NamedPlot r1, NamedPlot r2) {
      return Double.compare(r2.score, r1.score);
    }
  }

  /**
   * Instantiates a new filter analysis.
   */
  public FilterAnalysis() {
    isHeadless = java.awt.GraphicsEnvironment.isHeadless();
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    final String inputDirectory = getInputDirectory();
    if (inputDirectory == null) {
      return;
    }

    final List<MemoryPeakResults> resultsList = readResults(inputDirectory);
    if (resultsList.isEmpty()) {
      IJ.error(TITLE, "No results could be loaded. Check files have the suffix .xls, .csv or .bin");
      return;
    }

    // Load filters from file or generate from dialog input
    List<FilterSet> filterSets = null;
    final boolean fileInput = (arg != null && arg.contains("file"));

    if (!showDialog(resultsList, fileInput)) {
      return;
    }

    if (fileInput) {
      filterSets = readFilterSets();
      if (filterSets.isEmpty()) {
        return;
      }
    } else {
      filterSets = createFilters();
    }

    if (filterSets.isEmpty()) {
      IJ.error(TITLE, "No filters specified");
      return;
    }

    if (!fileInput && settings.saveFilterSets) {
      saveFilterSets(filterSets);
    }

    analyse(resultsList, filterSets);
  }

  private static String getInputDirectory() {
    final GUIFilterSettings.Builder filterSettings =
        SettingsManager.readGuiFilterSettings(0).toBuilder();

    if (filterSettings.getFilterAnalysisDirectory() != null) {
      OpenDialog.setDefaultDirectory(filterSettings.getFilterAnalysisDirectory());
    }
    filterSettings.setFilterAnalysisDirectory(IJ.getDirectory("Select results directory ..."));
    if (filterSettings.getFilterAnalysisDirectory() == null) {
      return null;
    }

    SettingsManager.writeSettings(filterSettings.build());

    return filterSettings.getFilterAnalysisDirectory();
  }

  @SuppressWarnings("unchecked")
  private static List<FilterSet> readFilterSets() {
    final GUIFilterSettings.Builder filterSettings =
        SettingsManager.readGuiFilterSettings(0).toBuilder();

    final String[] path = ImageJUtils.decodePath(filterSettings.getFilterSetFilename());
    final OpenDialog chooser = new OpenDialog("Filter_File", path[0], path[1]);
    if (chooser.getFileName() != null) {
      IJ.showStatus("Reading filters ...");
      filterSettings.setFilterSetFilename(chooser.getDirectory() + chooser.getFileName());

      try (BufferedReader input =
          Files.newBufferedReader(Paths.get(filterSettings.getFilterSetFilename()))) {
        // Use the instance so we can catch the exception
        final Object o = FilterXStreamUtils.getXStreamInstance().fromXML(input);
        if (o instanceof List<?>) {
          SettingsManager.writeSettings(filterSettings.build());
          return (List<FilterSet>) o;
        }
        IJ.log("No filter sets defined in the specified file: "
            + filterSettings.getFilterSetFilename());
      } catch (final Exception ex) {
        IJ.log("Unable to load the filter sets from file: " + ex.getMessage());
      } finally {
        IJ.showStatus("");
      }
    }
    return Collections.emptyList();
  }

  private static void saveFilterSets(List<FilterSet> filterSets) {
    final GUIFilterSettings.Builder filterSettings =
        SettingsManager.readGuiFilterSettings(0).toBuilder();

    final String[] path = ImageJUtils.decodePath(filterSettings.getFilterSetFilename());
    final OpenDialog chooser = new OpenDialog("Filter_File", path[0], path[1]);
    if (chooser.getFileName() != null) {
      filterSettings.setFilterSetFilename(chooser.getDirectory() + chooser.getFileName());
      try (BufferedWriter out =
          Files.newBufferedWriter(Paths.get(filterSettings.getFilterSetFilename()))) {
        // Use the instance so we can catch the exception
        FilterXStreamUtils.getXStreamInstance().toXML(filterSets, out);
        SettingsManager.writeSettings(filterSettings.build());
      } catch (final Exception ex) {
        IJ.log("Unable to save the filter sets to file: " + ex.getMessage());
      }
    }
  }

  private static List<MemoryPeakResults> readResults(String inputDirectory) {
    final LastResults last = lastResults.get();
    if (last != null && last.inputDirectory.equals(inputDirectory)) {
      final GenericDialog gd = new GenericDialog(TITLE);
      gd.addMessage("Re-use results from the same directory (no to refresh)?");
      gd.enableYesNoCancel();
      gd.hideCancelButton();
      gd.showDialog();
      if (gd.wasOKed()) {
        return last.resultsList;
      }
    }

    final List<MemoryPeakResults> list = new LinkedList<>();
    final File[] fileList = (new File(inputDirectory)).listFiles(
        (dir, name) -> name.endsWith(".xls") || name.endsWith(".csv") || name.endsWith(".bin"));
    if (fileList != null) {
      // Exclude directories
      for (int i = 0; i < fileList.length; i++) {
        if (fileList[i].isFile()) {
          IJ.showStatus(String.format("Reading results ... %d/%d", i + 1, fileList.length));
          IJ.showProgress(i, fileList.length);
          final PeakResultsReader reader = new PeakResultsReader(fileList[i].getPath());
          final MemoryPeakResults results = reader.getResults();
          if (results != null && results.size() > 0) {
            list.add(results);
          }
        }
      }
    }
    ImageJUtils.finished();
    lastResults.set(new LastResults(list, inputDirectory));

    return list;
  }

  private boolean showDialog(List<MemoryPeakResults> resultsList, boolean fileInput) {
    final GenericDialog gd = new GenericDialog(TITLE);
    String helpKey = "filter-analysis";
    if (fileInput) {
      helpKey += "-file";
    }
    gd.addHelp(HelpUrls.getUrl(helpKey));

    int total = 0;
    final Counter tp = new Counter();
    for (final MemoryPeakResults r : resultsList) {
      total += r.size();
      r.forEach((PeakResultProcedure) result -> {
        if (result.getOrigValue() != 0) {
          tp.increment();
        }
      });
    }
    ImageJUtils.addMessage(gd, "%d files, %d results, %d True-Positives", resultsList.size(), total,
        tp.getCount());

    settings = Settings.load();
    settings.save();

    if (!fileInput) {
      gd.addCheckbox("SNR_filter", settings.snrFilter);
      gd.addNumericField("Min_SNR", settings.minSnr, 0);
      gd.addNumericField("Max_SNR", settings.maxSnr, 0);
      gd.addNumericField("Min_Width", settings.minWidth, 2);
      gd.addNumericField("Max_Width", settings.maxWidth, 2);
      gd.addNumericField("Increment_Width", settings.incWidth, 2);

      gd.addCheckbox("Precision_filter", settings.precisionFilter);
      gd.addNumericField("Min_Precision", settings.minPrecision, 0);
      gd.addNumericField("Max_Precision", settings.maxPrecision, 0);

      gd.addCheckbox("Trace_filter", settings.traceFilter);
      gd.addNumericField("Min_distance", settings.minDistance, 2);
      gd.addNumericField("Max_distance", settings.maxDistance, 2);
      gd.addNumericField("Increment_distance", settings.incDistance, 2);
      gd.addNumericField("Min_time", settings.minTime, 0);
      gd.addNumericField("Max_time", settings.maxTime, 0);
      gd.addNumericField("Increment_time", settings.incTime, 0);

      gd.addCheckbox("Hysteresis_SNR_filter", settings.hysteresisSnrFilter);
      gd.addNumericField("Min_SNR_gap", settings.minSnrGap, 0);
      gd.addNumericField("Max_SNR_gap", settings.maxSnrGap, 0);
      gd.addNumericField("Increment_SNR_gap", settings.incSnrGap, 0);

      gd.addCheckbox("Hysteresis_Precision_filter", settings.hysteresisPrecisionFilter);
      gd.addNumericField("Min_Precision_gap", settings.minPrecisionGap, 0);
      gd.addNumericField("Max_Precision_gap", settings.maxPrecisionGap, 0);
      gd.addNumericField("Increment_Precision_gap", settings.incPrecisionGap, 0);

      gd.addCheckbox("Save_filters", settings.saveFilterSets);
    }
    gd.addCheckbox("Show_table", settings.showResultsTable);
    gd.addSlider("Plot_top_n", 0, 20, settings.plotTopN);
    gd.addCheckbox("Calculate_sensitivity", settings.calculateSensitivity);
    gd.addSlider("Delta", 0.01, 1, settings.delta);

    gd.showDialog();

    return !gd.wasCanceled() && readDialog(gd, fileInput);
  }

  private boolean readDialog(GenericDialog gd, boolean fileInput) {
    if (!fileInput) {
      settings.snrFilter = gd.getNextBoolean();
      settings.minSnr = (int) gd.getNextNumber();
      settings.maxSnr = (int) gd.getNextNumber();
      settings.minWidth = gd.getNextNumber();
      settings.maxWidth = gd.getNextNumber();
      settings.incWidth = gd.getNextNumber();

      settings.precisionFilter = gd.getNextBoolean();
      settings.minPrecision = (int) gd.getNextNumber();
      settings.maxPrecision = (int) gd.getNextNumber();

      settings.traceFilter = gd.getNextBoolean();
      settings.minDistance = gd.getNextNumber();
      settings.maxDistance = gd.getNextNumber();
      settings.incDistance = gd.getNextNumber();
      settings.minTime = (int) gd.getNextNumber();
      settings.maxTime = (int) gd.getNextNumber();
      settings.incTime = (int) gd.getNextNumber();

      settings.hysteresisSnrFilter = gd.getNextBoolean();
      settings.minSnrGap = (int) gd.getNextNumber();
      settings.maxSnrGap = (int) gd.getNextNumber();
      settings.incSnrGap = (int) gd.getNextNumber();

      settings.hysteresisPrecisionFilter = gd.getNextBoolean();
      settings.minPrecisionGap = (int) gd.getNextNumber();
      settings.maxPrecisionGap = (int) gd.getNextNumber();
      settings.incPrecisionGap = (int) gd.getNextNumber();

      settings.saveFilterSets = gd.getNextBoolean();
    }
    settings.showResultsTable = gd.getNextBoolean();
    settings.plotTopN = (int) Math.abs(gd.getNextNumber());
    settings.calculateSensitivity = gd.getNextBoolean();
    settings.delta = gd.getNextNumber();

    // Check there is one output
    if (!settings.showResultsTable && !settings.calculateSensitivity && settings.plotTopN < 1) {
      IJ.error(TITLE, "No output selected");
      return false;
    }

    // Check arguments
    try {
      if (!fileInput) {
        ParameterUtils.isPositive("Min SNR", settings.minSnr);
        ParameterUtils.isAboveZero("Max SNR", settings.maxSnr);
        ParameterUtils.isPositive("Min width", settings.minWidth);
        ParameterUtils.isAboveZero("Max width", settings.maxWidth);
        ParameterUtils.isAboveZero("Increment width", settings.incWidth);
        ParameterUtils.isPositive("Min precision", settings.minPrecision);
        ParameterUtils.isAboveZero("Max precision", settings.maxPrecision);
        ParameterUtils.isPositive("Min Distance", settings.minDistance);
        ParameterUtils.isAboveZero("Max Distance", settings.maxDistance);
        ParameterUtils.isAboveZero("Increment Distance", settings.incDistance);
        ParameterUtils.isPositive("Min Time", settings.minTime);
        ParameterUtils.isAboveZero("Max Time", settings.maxTime);
        ParameterUtils.isAboveZero("Increment Time", settings.incTime);
        ParameterUtils.isPositive("Min Snr Gap", settings.minSnrGap);
        ParameterUtils.isAboveZero("Max Snr Gap", settings.maxSnrGap);
        ParameterUtils.isAboveZero("Increment Snr Gap", settings.incSnrGap);
        ParameterUtils.isPositive("Min Precision Gap", settings.minPrecisionGap);
        ParameterUtils.isAboveZero("Max Precision Gap", settings.maxPrecisionGap);
        ParameterUtils.isAboveZero("Increment Precision Gap", settings.incPrecisionGap);
      }

      ParameterUtils.isAboveZero("Delta", settings.delta);
      ParameterUtils.isBelow("Delta", settings.delta, 1);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return !gd.invalidNumber();
  }

  private List<FilterSet> createFilters() {
    IJ.showStatus("Creating filters ...");
    final List<FilterSet> filterSets = new LinkedList<>();
    addSnrFilters(filterSets);
    addPrecisionFilters(filterSets);
    addTraceFilters(filterSets);
    addSnrHysteresisFilters(filterSets);
    addPrecisionHysteresisFilters(filterSets);
    IJ.showStatus("");
    return filterSets;
  }

  private void addSnrFilters(List<FilterSet> filterSets) {
    if (!settings.snrFilter) {
      return;
    }
    for (double w = settings.minWidth; w <= settings.maxWidth; w += settings.incWidth) {
      final WidthFilter wf = new WidthFilter((float) w);
      final List<Filter> filters = new LinkedList<>();
      for (int snr = settings.minSnr; snr <= settings.maxSnr; snr++) {
        filters.add(new AndFilter(wf, new SnrFilter(snr)));
      }
      filterSets.add(new FilterSet(filters));
    }
  }

  private void addPrecisionFilters(List<FilterSet> filterSets) {
    if (!settings.precisionFilter) {
      return;
    }
    final List<Filter> filters = new LinkedList<>();
    for (int p = settings.minPrecision; p <= settings.maxPrecision; p++) {
      filters.add(new PrecisionFilter(p));
    }
    filterSets.add(new FilterSet(filters));
  }

  private void addTraceFilters(List<FilterSet> filterSets) {
    if (!settings.traceFilter) {
      return;
    }
    for (double d = settings.minDistance; d <= settings.maxDistance; d += settings.incDistance) {
      final SnrFilter snr = new SnrFilter(settings.maxSnr);
      final List<Filter> filters = new LinkedList<>();
      for (int t = settings.minTime; t <= settings.maxTime; t += settings.incTime) {
        filters.add(new OrFilter(snr, new TraceFilter(d, t)));
      }
      filterSets.add(new FilterSet(filters));
    }
  }

  private void addSnrHysteresisFilters(List<FilterSet> filterSets) {
    if (!settings.hysteresisSnrFilter) {
      return;
    }
    for (double w = settings.minWidth; w <= settings.maxWidth; w += settings.incWidth) {
      final WidthFilter wf = new WidthFilter((float) w);
      for (int snrGap = settings.minSnrGap; snrGap <= settings.maxSnrGap;
          snrGap += settings.incSnrGap) {
        final List<Filter> filters = new LinkedList<>();
        for (int snr = settings.minSnr; snr <= settings.maxSnr; snr++) {
          filters.add(new AndFilter(wf, new SnrHysteresisFilter(2, 0, 1, 0, snr, snrGap)));
        }
        filterSets.add(new FilterSet(filters));
      }
    }
  }

  private void addPrecisionHysteresisFilters(List<FilterSet> filterSets) {
    if (!settings.hysteresisPrecisionFilter) {
      return;
    }
    for (int precisionGap = settings.minPrecisionGap; precisionGap <= settings.maxPrecisionGap;
        precisionGap += settings.incPrecisionGap) {
      final List<Filter> filters = new LinkedList<>();
      for (int precision = settings.minPrecision; precision <= settings.maxPrecision; precision++) {
        filters.add(new PrecisionHysteresisFilter(2, 0, 1, 0, precision, precisionGap));
      }
      filterSets.add(new FilterSet(filters));
    }
  }

  /**
   * Run different filtering methods on a set of labelled peak results outputting performance
   * statistics on the success of the filter to an ImageJ table.
   *
   * <p>If the peak result original value is set to 1 it is considered a true peak, 0 for a false
   * peak. Filtering is done using e.g. SNR threshold, Precision thresholds, etc. The statistics
   * reported are shown in a table, e.g. precision, Jaccard, F-score.
   *
   * <p>For each filter set a plot is shown of the Jaccard score verses the filter value, thus
   * filters should be provided in ascending numerical order otherwise they are sorted.
   *
   * @param resultsList the results list
   * @param filterSets the filter sets
   */
  public void analyse(List<MemoryPeakResults> resultsList, List<FilterSet> filterSets) {
    final Consumer<String> output = createResultsWindow();
    plots = new ArrayList<>(settings.plotTopN);
    bestFilter = new HashMap<>();
    bestFilterOrder = new LinkedList<>();

    IJ.showStatus("Analysing filters ...");
    final int total = countFilters(filterSets);
    int count = 0;
    for (final FilterSet filterSet : filterSets) {
      IJ.showStatus("Analysing " + filterSet.getName() + " ...");
      count = runAnalysis(output, filterSet, resultsList, count, total);
    }
    ImageJUtils.finished();

    showPlots();
    calculateSensitivity(resultsList);
  }

  private static int countFilters(List<FilterSet> filterSets) {
    int count = 0;
    for (final FilterSet filterSet : filterSets) {
      count += filterSet.size();
    }
    return count;
  }

  private void showPlots() {
    if (plots.isEmpty()) {
      return;
    }

    // Display the top N plots
    final int[] list = new int[plots.size()];
    int index = 0;
    for (final NamedPlot p : plots) {
      final Plot plot = new Plot(p.name, p.xAxisName, "Jaccard");
      plot.addPoints(p.xValues, p.yValues, Plot.LINE);
      plot.setLimits(p.xValues[0], p.xValues[p.xValues.length - 1], 0, 1);
      plot.setColor(Color.RED);
      plot.draw();
      plot.setColor(Color.BLUE);
      plot.addPoints(p.xValues, p.yValues, Plot.CROSS);
      final PlotWindow plotWindow = ImageJUtils.display(p.name, plot);
      list[index++] = plotWindow.getImagePlus().getID();
    }
    WindowOrganiser.tileWindows(list);
  }

  private void calculateSensitivity(List<MemoryPeakResults> resultsList) {
    if (!settings.calculateSensitivity) {
      return;
    }
    if (!bestFilter.isEmpty()) {
      IJ.showStatus("Calculating sensitivity ...");
      final Consumer<String> output = createSensitivityWindow();

      int currentIndex = 0;
      for (final String type : bestFilterOrder) {
        IJ.showProgress(currentIndex++, bestFilter.size());

        final Filter filter = bestFilter.get(type).filter;

        final ClassificationResult s = filter.score(resultsList);

        final String message = type + "\t\t\t" + MathUtils.rounded(s.getJaccard(), 4) + "\t\t"
            + MathUtils.rounded(s.getPrecision(), 4) + "\t\t" + MathUtils.rounded(s.getRecall(), 4);

        output.accept(message);

        // List all the parameters that can be adjusted.
        final int parameters = filter.getNumberOfParameters();
        for (int index = 0; index < parameters; index++) {
          // For each parameter compute as upward + downward delta and get the average gradient
          final Filter higher = filter.adjustParameter(index, settings.delta);
          final Filter lower = filter.adjustParameter(index, -settings.delta);

          final ClassificationResult sHigher = higher.score(resultsList);
          final ClassificationResult sLower = lower.score(resultsList);

          final StringBuilder sb = new StringBuilder();
          sb.append('\t').append(filter.getParameterName(index)).append('\t');
          sb.append(MathUtils.rounded(filter.getParameterValue(index), 4)).append('\t');

          final double dx1 = higher.getParameterValue(index) - filter.getParameterValue(index);
          final double dx2 = filter.getParameterValue(index) - lower.getParameterValue(index);
          addSensitivityScore(sb, s.getJaccard(), sHigher.getJaccard(), sLower.getJaccard(), dx1,
              dx2);
          addSensitivityScore(sb, s.getPrecision(), sHigher.getPrecision(), sLower.getPrecision(),
              dx1, dx2);
          addSensitivityScore(sb, s.getRecall(), sHigher.getRecall(), sLower.getRecall(), dx1, dx2);

          output.accept(sb.toString());
        }
      }

      output.accept("-=-=-=-");

      ImageJUtils.finished();
    }
  }

  private static void addSensitivityScore(StringBuilder sb, double score, double s1, double s2,
      double dx1, double dx2) {
    // Use absolute in case this is not a local maximum. We are mainly interested in how
    // flat the curve is at this point in relation to parameter changes.
    final double abs1 = Math.abs(score - s1);
    final double abs2 = Math.abs(score - s2);
    final double dydx1 = (abs1) / dx1;
    final double dydx2 = (abs2) / dx2;
    final double relativeSensitivity = (abs1 + abs2) * 0.5;
    final double sensitivity = (dydx1 + dydx2) * 0.5;

    sb.append(MathUtils.rounded(relativeSensitivity, 4)).append('\t');
    sb.append(MathUtils.rounded(sensitivity, 4)).append('\t');
  }

  private Consumer<String> createResultsWindow() {
    if (!settings.showResultsTable) {
      return null;
    }

    if (isHeadless) {
      IJ.log(createResultsHeader());
      return IJ::log;
    }

    return ImageJUtils.refresh(resultsWindowRef, () -> {
      final String header = createResultsHeader();
      return new TextWindow(TITLE + " Results", header, "", 900, 300);
    })::append;
  }

  private static String createResultsHeader() {
    final StringBuilder sb = new StringBuilder();
    sb.append("Name\t");
    sb.append("N\t");
    sb.append("TP\t");
    sb.append("FP\t");
    sb.append("TN\t");
    sb.append("FN\t");
    sb.append("Jaccard\t");
    sb.append("Precision\t");
    sb.append("Recall\t");
    sb.append("F1");
    return sb.toString();
  }

  private Consumer<String> createSensitivityWindow() {
    if (isHeadless) {
      IJ.log(createSensitivityHeader());
      return IJ::log;
    }
    return ImageJUtils.refresh(sensitivityWindowRef, () -> {
      final String header = createSensitivityHeader();
      return new TextWindow(TITLE + " Sensitivity", header, "", 900, 300);
    })::append;
  }

  private static String createSensitivityHeader() {
    final StringBuilder sb = new StringBuilder();
    sb.append("Filter\t");
    sb.append("Param\t");
    sb.append("Value\t");
    sb.append("J Sensitivity (delta)\t");
    sb.append("J Sensitivity (unit)\t");
    sb.append("P Sensitivity (delta)\t");
    sb.append("P Sensitivity (unit)\t");
    sb.append("R Sensitivity (delta)\t");
    sb.append("R Sensitivity (unit)\t");
    return sb.toString();
  }

  private int runAnalysis(Consumer<String> output, FilterSet filterSet,
      List<MemoryPeakResults> resultsList, int count, final int total) {
    final double[] xValues = (isHeadless) ? null : new double[filterSet.size()];
    final double[] yValues = (isHeadless) ? null : new double[filterSet.size()];
    int index = 0;

    filterSet.sort();

    // Track if all the filters are the same type. If so then we can calculate the sensitivity of
    // each parameter.
    String type = null;
    boolean allSameType = true;
    Filter maxFilter = null;
    double maxScore = -1;

    for (final Filter filter : filterSet.getFilters()) {
      if (count++ % 16 == 0) {
        IJ.showProgress(count, total);
      }

      final ClassificationResult s = runFilter(output, filter, resultsList);

      if (type == null) {
        type = filter.getType();
      } else if (!type.equals(filter.getType())) {
        allSameType = false;
      }

      final double jaccard = s.getJaccard();
      if (maxScore < jaccard) {
        maxScore = jaccard;
        maxFilter = filter;
      }

      if (xValues != null && yValues != null) {
        xValues[index] = filter.getNumericalValue();
        yValues[index++] = jaccard;
      }
    }

    if (allSameType && settings.calculateSensitivity) {
      final FilterScore filterScore = bestFilter.get(type);
      if (filterScore != null) {
        if (filterScore.score < maxScore) {
          filterScore.update(maxFilter, maxScore);
        }
      } else {
        bestFilter.put(type, new FilterScore(maxFilter, maxScore));
        bestFilterOrder.add(type);
      }
    }

    // Add spacer at end of each result set
    if (settings.showResultsTable) {
      output.accept("");
    }

    if (!isHeadless && settings.plotTopN > 0 && xValues != null) {
      // Check the xValues are unique. Since the filters have been sorted by their
      // numeric value we only need to compare adjacent entries.
      boolean unique = true;
      for (int ii = 0; ii < xValues.length - 1; ii++) {
        if (xValues[ii] == xValues[ii + 1]) {
          unique = false;
          break;
        }
      }
      String xAxisName = filterSet.getValueName();
      // Check the values all refer to the same property
      for (final Filter filter : filterSet.getFilters()) {
        if (!xAxisName.equals(filter.getNumericalValueName())) {
          unique = false;
          break;
        }
      }
      if (!unique) {
        // If not unique then renumber them and use an arbitrary label
        xAxisName = "Filter";
        for (int ii = 0; ii < xValues.length; ii++) {
          xValues[ii] = ii + 1.0;
        }
      }

      final String title = filterSet.getName();

      // Check if a previous filter set had the same name, update if necessary
      final NamedPlot plot = getNamedPlot(title);
      if (plot == null) {
        plots.add(new NamedPlot(title, xAxisName, xValues, yValues));
      } else {
        plot.updateValues(xAxisName, xValues, yValues);
      }

      if (plots.size() > settings.plotTopN) {
        Collections.sort(plots, NamedPlot::compare);
        plots.remove(plots.size() - 1);
      }
    }

    return count;
  }

  private NamedPlot getNamedPlot(String title) {
    for (final NamedPlot p : plots) {
      if (p.name.equals(title)) {
        return p;
      }
    }
    return null;
  }

  private ClassificationResult runFilter(Consumer<String> output, Filter filter,
      List<MemoryPeakResults> resultsList) {
    final ClassificationResult s = filter.score(resultsList);

    if (settings.showResultsTable) {
      final StringBuilder sb = new StringBuilder();
      sb.append(filter.getName()).append('\t');
      sb.append(s.getTruePositives() + s.getFalsePositives()).append('\t');
      sb.append(s.getTruePositives()).append('\t');
      sb.append(s.getFalsePositives()).append('\t');
      sb.append(s.getTrueNegatives()).append('\t');
      sb.append(s.getFalseNegatives()).append('\t');
      sb.append(s.getJaccard()).append('\t');
      sb.append(s.getPrecision()).append('\t');
      sb.append(s.getRecall()).append('\t');
      sb.append(s.getF1Score());
      output.accept(sb.toString());
    }
    return s;
  }
}
