/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
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

import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.Plot2;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.match.ClassificationResult;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.UnicodeReader;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.GUIFilterSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.PeakResultsReader;
import uk.ac.sussex.gdsc.smlm.results.count.Counter;
import uk.ac.sussex.gdsc.smlm.results.filter.AndFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.Filter;
import uk.ac.sussex.gdsc.smlm.results.filter.FilterSet;
import uk.ac.sussex.gdsc.smlm.results.filter.OrFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.PrecisionFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.PrecisionHysteresisFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.SNRFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.SNRHysteresisFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.TraceFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.WidthFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.XStreamWrapper;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;

import ij.IJ;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import ij.text.TextWindow;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

/**
 * Run different filtering methods on a set of labelled peak results outputting performance
 * statistics on the success of the filter. <p> All results files in a specified directory are read.
 * If the peak result original value is set to 1 it is considered a true peak, 0 for a false peak.
 * Filtering is done using e.g. SNR threshold, Precision thresholds, etc. The statistics reported
 * are shown in a table, e.g. precision, Jaccard, F-score.
 */
public class FilterAnalysis implements PlugIn {
  private static final String TITLE = "Filter Analysis";
  private static TextWindow resultsWindow = null;
  private static TextWindow sensitivityWindow = null;

  private static boolean saveFilterSets = false;
  private static boolean showResultsTable = true;
  private static int plotTopN = 0;
  private ArrayList<NamedPlot> plots;
  private static boolean calculateSensitivity = false;
  private static double delta = 0.1;
  private HashMap<String, FilterScore> bestFilter;
  private LinkedList<String> bestFilterOrder;

  private static boolean snrFilter = true;
  private static int minSnr = 20;
  private static int maxSnr = 80;
  private static double minWidth = 1.5;
  private static double maxWidth = 2.0;
  private static double incWidth = 0.5;

  private static boolean precisionFilter = true;
  private static int minPrecision = 20;
  private static int maxPrecision = 80;

  private static boolean traceFilter = false;
  private static double minDistance = 0.3;
  private static double maxDistance = 1.2;
  private static double incDistance = 0.3;
  private static int minTime = 1;
  private static int maxTime = 80;
  private static int incTime = 10;

  private static boolean hysteresisSnrFilter = true;
  private static int minSnrGap = 10;
  private static int maxSnrGap = 40;
  private static int incSnrGap = 10;

  private static boolean hysteresisPrecisionFilter = true;
  private static int minPrecisionGap = 10;
  private static int maxPrecisionGap = 40;
  private static int incPrecisionGap = 10;

  private static List<MemoryPeakResults> resultsList = null;
  private String inputDirectory;
  private static String lastInputDirectory = "";

  private final boolean isHeadless;

  /**
   * Instantiates a new filter analysis.
   */
  public FilterAnalysis() {
    isHeadless = java.awt.GraphicsEnvironment.isHeadless();
  }

  /** {@inheritDoc} */
  @Override
  public void run(String arg) {
    SMLMUsageTracker.recordPlugin(this.getClass(), arg);

    if (getInputDirectory() == null) {
      return;
    }

    resultsList = readResults();
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
      if (filterSets == null) {
        return;
      }
    } else {
      filterSets = createFilters();
    }

    if (filterSets == null || filterSets.isEmpty()) {
      IJ.error(TITLE, "No filters specified");
      return;
    }

    if (!fileInput && saveFilterSets) {
      saveFilterSets(filterSets);
    }

    analyse(resultsList, filterSets);
  }

  private String getInputDirectory() {
    final GUIFilterSettings.Builder filterSettings =
        SettingsManager.readGUIFilterSettings(0).toBuilder();

    if (filterSettings.getFilterAnalysisDirectory() != null) {
      OpenDialog.setDefaultDirectory(filterSettings.getFilterAnalysisDirectory());
    }
    filterSettings.setFilterAnalysisDirectory(IJ.getDirectory("Select results directory ..."));
    if (filterSettings.getFilterAnalysisDirectory() == null) {
      return null;
    }

    SettingsManager.writeSettings(filterSettings.build());

    return inputDirectory = filterSettings.getFilterAnalysisDirectory();
  }

  @SuppressWarnings("unchecked")
  private static List<FilterSet> readFilterSets() {
    final GUIFilterSettings.Builder filterSettings =
        SettingsManager.readGUIFilterSettings(0).toBuilder();

    final String[] path = ImageJUtils.decodePath(filterSettings.getFilterSetFilename());
    final OpenDialog chooser = new OpenDialog("Filter_File", path[0], path[1]);
    if (chooser.getFileName() != null) {
      IJ.showStatus("Reading filters ...");
      filterSettings.setFilterSetFilename(chooser.getDirectory() + chooser.getFileName());

      try (BufferedReader input = new BufferedReader(
          new UnicodeReader(new FileInputStream(filterSettings.getFilterSetFilename()), null))) {
        final Object o = XStreamWrapper.getInstance().fromXML(input);
        if (o != null && o instanceof List<?>) {
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
    return null;
  }

  private static void saveFilterSets(List<FilterSet> filterSets) {
    final GUIFilterSettings.Builder filterSettings =
        SettingsManager.readGUIFilterSettings(0).toBuilder();

    final String[] path = ImageJUtils.decodePath(filterSettings.getFilterSetFilename());
    final OpenDialog chooser = new OpenDialog("Filter_File", path[0], path[1]);
    if (chooser.getFileName() != null) {
      filterSettings.setFilterSetFilename(chooser.getDirectory() + chooser.getFileName());
      try (OutputStreamWriter out = new OutputStreamWriter(
          new FileOutputStream(filterSettings.getFilterSetFilename()), "UTF-8")) {
        XStreamWrapper.getInstance().toXML(filterSets, out);
        SettingsManager.writeSettings(filterSettings.build());
      } catch (final Exception ex) {
        IJ.log("Unable to save the filter sets to file: " + ex.getMessage());
      }
    }
  }

  private List<MemoryPeakResults> readResults() {
    if (resultsList != null && inputDirectory.equals(lastInputDirectory)) {
      final GenericDialog gd = new GenericDialog(TITLE);
      gd.addMessage("Re-use results from the same directory (no to refresh)?");
      gd.enableYesNoCancel();
      gd.hideCancelButton();
      gd.showDialog();
      if (gd.wasOKed()) {
        return resultsList;
      }
    }

    final List<MemoryPeakResults> resultsList = new LinkedList<>();
    final File[] fileList = (new File(inputDirectory)).listFiles(new FilenameFilter() {
      @Override
      public boolean accept(File dir, String name) {
        return (name.endsWith(".xls") || name.endsWith(".csv") || name.endsWith(".bin"));
      }
    });
    if (fileList != null) {
      // Exclude directories
      for (int i = 0; i < fileList.length; i++) {
        if (fileList[i].isFile()) {
          IJ.showStatus(String.format("Reading results ... %d/%d", i + 1, fileList.length));
          IJ.showProgress(i, fileList.length);
          final PeakResultsReader reader = new PeakResultsReader(fileList[i].getPath());
          final MemoryPeakResults results = reader.getResults();
          if (results != null && results.size() > 0) {
            resultsList.add(results);
          }
        }
      }
    }
    IJ.showStatus("");
    IJ.showProgress(1);
    lastInputDirectory = inputDirectory;

    return resultsList;
  }

  private static boolean showDialog(List<MemoryPeakResults> resultsList, boolean fileInput) {
    final GenericDialog gd = new GenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);

    int total = 0;
    final Counter tp = new Counter();
    for (final MemoryPeakResults r : resultsList) {
      total += r.size();
      r.forEach(new PeakResultProcedure() {
        @Override
        public void execute(PeakResult p) {
          if (p.getOrigValue() != 0) {
            tp.increment();
          }
        }
      });
    }
    gd.addMessage(
        String.format("%d files, %d results, %d True-Positives", resultsList.size(), total, tp));

    if (!fileInput) {
      gd.addCheckbox("SNR_filter", snrFilter);
      gd.addNumericField("Min_SNR", minSnr, 0);
      gd.addNumericField("Max_SNR", maxSnr, 0);
      gd.addNumericField("Min_Width", minWidth, 2);
      gd.addNumericField("Max_Width", maxWidth, 2);
      gd.addNumericField("Increment_Width", incWidth, 2);

      gd.addCheckbox("Precision_filter", precisionFilter);
      gd.addNumericField("Min_Precision", minPrecision, 0);
      gd.addNumericField("Max_Precision", maxPrecision, 0);

      gd.addCheckbox("Trace_filter", traceFilter);
      gd.addNumericField("Min_distance", minDistance, 2);
      gd.addNumericField("Max_distance", maxDistance, 2);
      gd.addNumericField("Increment_distance", incDistance, 2);
      gd.addNumericField("Min_time", minTime, 0);
      gd.addNumericField("Max_time", maxTime, 0);
      gd.addNumericField("Increment_time", incTime, 0);

      gd.addCheckbox("Hysteresis_SNR_filter", hysteresisSnrFilter);
      gd.addNumericField("Min_SNR_gap", minSnrGap, 0);
      gd.addNumericField("Max_SNR_gap", maxSnrGap, 0);
      gd.addNumericField("Increment_SNR_gap", incSnrGap, 0);

      gd.addCheckbox("Hysteresis_Precision_filter", hysteresisPrecisionFilter);
      gd.addNumericField("Min_Precision_gap", minPrecisionGap, 0);
      gd.addNumericField("Max_Precision_gap", maxPrecisionGap, 0);
      gd.addNumericField("Increment_Precision_gap", incPrecisionGap, 0);

      gd.addCheckbox("Save_filters", saveFilterSets);
    }
    gd.addCheckbox("Show_table", showResultsTable);
    gd.addSlider("Plot_top_n", 0, 20, plotTopN);
    gd.addCheckbox("Calculate_sensitivity", calculateSensitivity);
    gd.addSlider("Delta", 0.01, 1, delta);

    gd.showDialog();

    if (gd.wasCanceled() || !readDialog(gd, fileInput)) {
      return false;
    }

    return true;
  }

  private static boolean readDialog(GenericDialog gd, boolean fileInput) {
    if (!fileInput) {
      snrFilter = gd.getNextBoolean();
      minSnr = (int) gd.getNextNumber();
      maxSnr = (int) gd.getNextNumber();
      minWidth = gd.getNextNumber();
      maxWidth = gd.getNextNumber();
      incWidth = gd.getNextNumber();

      precisionFilter = gd.getNextBoolean();
      minPrecision = (int) gd.getNextNumber();
      maxPrecision = (int) gd.getNextNumber();

      traceFilter = gd.getNextBoolean();
      minDistance = gd.getNextNumber();
      maxDistance = gd.getNextNumber();
      incDistance = gd.getNextNumber();
      minTime = (int) gd.getNextNumber();
      maxTime = (int) gd.getNextNumber();
      incTime = (int) gd.getNextNumber();

      hysteresisSnrFilter = gd.getNextBoolean();
      minSnrGap = (int) gd.getNextNumber();
      maxSnrGap = (int) gd.getNextNumber();
      incSnrGap = (int) gd.getNextNumber();

      hysteresisPrecisionFilter = gd.getNextBoolean();
      minPrecisionGap = (int) gd.getNextNumber();
      maxPrecisionGap = (int) gd.getNextNumber();
      incPrecisionGap = (int) gd.getNextNumber();

      saveFilterSets = gd.getNextBoolean();
    }
    showResultsTable = gd.getNextBoolean();
    plotTopN = (int) Math.abs(gd.getNextNumber());
    calculateSensitivity = gd.getNextBoolean();
    delta = gd.getNextNumber();

    // Check there is one output
    if (!showResultsTable && !calculateSensitivity && plotTopN < 1) {
      IJ.error(TITLE, "No output selected");
      return false;
    }

    // Check arguments
    try {
      if (!fileInput) {
        Parameters.isPositive("Min SNR", minSnr);
        Parameters.isAboveZero("Max SNR", maxSnr);
        Parameters.isPositive("Min width", minWidth);
        Parameters.isAboveZero("Max width", maxWidth);
        Parameters.isAboveZero("Increment width", incWidth);
        Parameters.isPositive("Min precision", minPrecision);
        Parameters.isAboveZero("Max precision", maxPrecision);
        Parameters.isPositive("Min Distance", minDistance);
        Parameters.isAboveZero("Max Distance", maxDistance);
        Parameters.isAboveZero("Increment Distance", incDistance);
        Parameters.isPositive("Min Time", minTime);
        Parameters.isAboveZero("Max Time", maxTime);
        Parameters.isAboveZero("Increment Time", incTime);
        Parameters.isPositive("Min Snr Gap", minSnrGap);
        Parameters.isAboveZero("Max Snr Gap", maxSnrGap);
        Parameters.isAboveZero("Increment Snr Gap", incSnrGap);
        Parameters.isPositive("Min Precision Gap", minPrecisionGap);
        Parameters.isAboveZero("Max Precision Gap", maxPrecisionGap);
        Parameters.isAboveZero("Increment Precision Gap", incPrecisionGap);
      }

      Parameters.isAboveZero("Delta", delta);
      Parameters.isBelow("Delta", delta, 1);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return !gd.invalidNumber();
  }

  private static List<FilterSet> createFilters() {
    IJ.showStatus("Creating filters ...");
    final List<FilterSet> filterSets = new LinkedList<>();
    addSNRFilters(filterSets);
    addPrecisionFilters(filterSets);
    addTraceFilters(filterSets);
    addSNRHysteresisFilters(filterSets);
    addPrecisionHysteresisFilters(filterSets);
    IJ.showStatus("");
    return filterSets;
  }

  private static void addSNRFilters(List<FilterSet> filterSets) {
    if (!snrFilter) {
      return;
    }
    for (double w = minWidth; w <= maxWidth; w += incWidth) {
      final WidthFilter wf = new WidthFilter((float) w);
      final List<Filter> filters = new LinkedList<>();
      for (int snr = minSnr; snr <= maxSnr; snr++) {
        filters.add(new AndFilter(wf, new SNRFilter(snr)));
      }
      filterSets.add(new FilterSet(filters));
    }
  }

  private static void addPrecisionFilters(List<FilterSet> filterSets) {
    if (!precisionFilter) {
      return;
    }
    final List<Filter> filters = new LinkedList<>();
    for (int p = minPrecision; p <= maxPrecision; p++) {
      filters.add(new PrecisionFilter(p));
    }
    filterSets.add(new FilterSet(filters));
  }

  private static void addTraceFilters(List<FilterSet> filterSets) {
    if (!traceFilter) {
      return;
    }
    for (double d = minDistance; d <= maxDistance; d += incDistance) {
      final SNRFilter snr = new SNRFilter(maxSnr);
      final List<Filter> filters = new LinkedList<>();
      for (int t = minTime; t <= maxTime; t += incTime) {
        filters.add(new OrFilter(snr, new TraceFilter(d, t)));
      }
      filterSets.add(new FilterSet(filters));
    }
  }

  private static void addSNRHysteresisFilters(List<FilterSet> filterSets) {
    if (!hysteresisSnrFilter) {
      return;
    }
    for (double w = minWidth; w <= maxWidth; w += incWidth) {
      final WidthFilter wf = new WidthFilter((float) w);
      for (int snrGap = minSnrGap; snrGap <= maxSnrGap; snrGap += incSnrGap) {
        final List<Filter> filters = new LinkedList<>();
        for (int snr = minSnr; snr <= maxSnr; snr++) {
          filters.add(new AndFilter(wf, new SNRHysteresisFilter(2, 0, 1, 0, snr, snrGap)));
        }
        filterSets.add(new FilterSet(filters));
      }
    }
  }

  private static void addPrecisionHysteresisFilters(List<FilterSet> filterSets) {
    if (!hysteresisPrecisionFilter) {
      return;
    }
    for (int precisionGap = minPrecisionGap; precisionGap <= maxPrecisionGap; precisionGap +=
        incPrecisionGap) {
      final List<Filter> filters = new LinkedList<>();
      for (int precision = minPrecision; precision <= maxPrecision; precision++) {
        filters.add(new PrecisionHysteresisFilter(2, 0, 1, 0, precision, precisionGap));
      }
      filterSets.add(new FilterSet(filters));
    }
  }

  /**
   * Run different filtering methods on a set of labelled peak results outputting performance
   * statistics on the success of the filter to an ImageJ table. <p> If the peak result original
   * value is set to 1 it is considered a true peak, 0 for a false peak. Filtering is done using
   * e.g. SNR threshold, Precision thresholds, etc. The statistics reported are shown in a table,
   * e.g. precision, Jaccard, F-score. <p> For each filter set a plot is shown of the Jaccard score
   * verses the filter value, thus filters should be provided in ascending numerical order otherwise
   * they are sorted.
   *
   * @param resultsList the results list
   * @param filterSets the filter sets
   */
  public void analyse(List<MemoryPeakResults> resultsList, List<FilterSet> filterSets) {
    createResultsWindow();
    plots = new ArrayList<>(plotTopN);
    bestFilter = new HashMap<>();
    bestFilterOrder = new LinkedList<>();

    IJ.showStatus("Analysing filters ...");
    final int total = countFilters(filterSets);
    int count = 0;
    for (final FilterSet filterSet : filterSets) {
      IJ.showStatus("Analysing " + filterSet.getName() + " ...");
      count = run(filterSet, resultsList, count, total);
    }
    IJ.showProgress(1);
    IJ.showStatus("");

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
    int i = 0;
    for (final NamedPlot p : plots) {
      final Plot2 plot = new Plot2(p.name, p.xAxisName, "Jaccard", p.xValues, p.yValues);
      plot.setLimits(p.xValues[0], p.xValues[p.xValues.length - 1], 0, 1);
      plot.setColor(Color.RED);
      plot.draw();
      plot.setColor(Color.BLUE);
      plot.addPoints(p.xValues, p.yValues, Plot.CROSS);
      final PlotWindow plotWindow = ImageJUtils.display(p.name, plot);
      list[i++] = plotWindow.getImagePlus().getID();
    }
    WindowOrganiser.tileWindows(list);
  }

  private void calculateSensitivity(List<MemoryPeakResults> resultsList) {
    if (!calculateSensitivity) {
      return;
    }
    if (!bestFilter.isEmpty()) {
      IJ.showStatus("Calculating sensitivity ...");
      createSensitivityWindow();

      int currentIndex = 0;
      for (final String type : bestFilterOrder) {
        IJ.showProgress(currentIndex++, bestFilter.size());

        final Filter filter = bestFilter.get(type).filter;

        final ClassificationResult s = filter.score(resultsList);

        final String message = type + "\t\t\t" + MathUtils.rounded(s.getJaccard(), 4) + "\t\t"
            + MathUtils.rounded(s.getPrecision(), 4) + "\t\t" + MathUtils.rounded(s.getRecall(), 4);

        if (isHeadless) {
          IJ.log(message);
        } else {
          sensitivityWindow.append(message);
        }

        // List all the parameters that can be adjusted.
        final int parameters = filter.getNumberOfParameters();
        for (int index = 0; index < parameters; index++) {
          // For each parameter compute as upward + downward delta and get the average gradient
          final Filter higher = filter.adjustParameter(index, delta);
          final Filter lower = filter.adjustParameter(index, -delta);

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

          if (isHeadless) {
            IJ.log(sb.toString());
          } else {
            sensitivityWindow.append(sb.toString());
          }
        }
      }

      final String message = "-=-=-=-";
      if (isHeadless) {
        IJ.log(message);
      } else {
        sensitivityWindow.append(message);
      }

      IJ.showProgress(1);
      IJ.showStatus("");
    }
  }

  private static void addSensitivityScore(StringBuilder sb, double s, double s1, double s2,
      double dx1, double dx2) {
    // Use absolute in case this is not a local maximum. We are mainly interested in how
    // flat the curve is at this point in relation to parameter changes.
    final double abs1 = Math.abs(s - s1);
    final double abs2 = Math.abs(s - s2);
    final double dydx1 = (abs1) / dx1;
    final double dydx2 = (abs2) / dx2;
    final double relativeSensitivity = (abs1 + abs2) * 0.5;
    final double sensitivity = (dydx1 + dydx2) * 0.5;

    sb.append(MathUtils.rounded(relativeSensitivity, 4)).append('\t');
    sb.append(MathUtils.rounded(sensitivity, 4)).append('\t');
  }

  private void createResultsWindow() {
    if (!showResultsTable) {
      return;
    }

    if (isHeadless) {
      IJ.log(createResultsHeader());
    } else if (resultsWindow == null || !resultsWindow.isShowing()) {
      final String header = createResultsHeader();
      resultsWindow = new TextWindow(TITLE + " Results", header, "", 900, 300);
    }
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

  private void createSensitivityWindow() {
    if (isHeadless) {
      IJ.log(createSensitivityHeader());
    } else if (sensitivityWindow == null || !sensitivityWindow.isShowing()) {
      final String header = createSensitivityHeader();
      sensitivityWindow = new TextWindow(TITLE + " Sensitivity", header, "", 900, 300);
    }
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

  private int run(FilterSet filterSet, List<MemoryPeakResults> resultsList, int count,
      final int total) {
    final double[] xValues = (isHeadless) ? null : new double[filterSet.size()];
    final double[] yValues = (isHeadless) ? null : new double[filterSet.size()];
    int i = 0;

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

      final ClassificationResult s = run(filter, resultsList);

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
        xValues[i] = filter.getNumericalValue();
        yValues[i++] = jaccard;
      }
    }

    if (allSameType && calculateSensitivity) {
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
    if (isHeadless) {
      if (showResultsTable) {
        IJ.log("");
      }
    } else {
      if (showResultsTable) {
        resultsWindow.append("");
      }

      if (plotTopN > 0 && xValues != null) {
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
            xValues[ii] = ii + 1;
          }
        }

        final String title = filterSet.getName();

        // Check if a previous filter set had the same name, update if necessary
        NamedPlot p = getNamedPlot(title);
        if (p == null) {
          plots.add(new NamedPlot(title, xAxisName, xValues, yValues));
        } else {
          p.updateValues(xAxisName, xValues, yValues);
        }

        if (plots.size() > plotTopN) {
          Collections.sort(plots);
          p = plots.remove(plots.size() - 1);
        }
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

  private static double getMaximum(double[] values) {
    double max = values[0];
    for (int i = 1; i < values.length; i++) {
      if (values[i] > max) {
        max = values[i];
      }
    }
    return max;
  }

  private ClassificationResult run(Filter filter, List<MemoryPeakResults> resultsList) {
    final ClassificationResult s = filter.score(resultsList);

    if (showResultsTable) {
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
      if (isHeadless) {
        IJ.log(sb.toString());
      } else {
        resultsWindow.append(sb.toString());
      }
    }
    return s;
  }

  private class FilterScore {
    Filter filter;
    double score;

    public FilterScore(Filter filter, double score) {
      update(filter, score);
    }

    public void update(Filter filter, double score) {
      this.filter = filter;
      this.score = score;
    }
  }

  private class NamedPlot implements Comparable<NamedPlot> {
    String name, xAxisName;
    double[] xValues, yValues;
    double score;

    public NamedPlot(String name, String xAxisName, double[] xValues, double[] yValues) {
      this.name = name;
      updateValues(xAxisName, xValues, yValues);
    }

    public void updateValues(String xAxisName, double[] xValues, double[] yValues) {
      this.xAxisName = xAxisName;
      this.xValues = xValues;
      this.yValues = yValues;
      this.score = getMaximum(yValues);
    }

    @Override
    public int compareTo(NamedPlot o) {
      return Double.compare(o.score, score);
    }
  }
}
