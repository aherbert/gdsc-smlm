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

import uk.ac.sussex.gdsc.core.ags.utils.data.trees.gen2.IntResultHeap;
import uk.ac.sussex.gdsc.core.ij.BufferedTextWindow;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.NonBlockingExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.logging.Ticker;
import uk.ac.sussex.gdsc.core.utils.BooleanArray;
import uk.ac.sussex.gdsc.core.utils.BooleanRollingArray;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.SortUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.TurboList;
import uk.ac.sussex.gdsc.core.utils.concurrent.ConcurrencyUtils;
import uk.ac.sussex.gdsc.core.utils.concurrent.ConcurrentMonoStack;
import uk.ac.sussex.gdsc.smlm.data.NamedObject;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.FailCountManagerSettings;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsSettings;
import uk.ac.sussex.gdsc.smlm.engine.FitEngine;
import uk.ac.sussex.gdsc.smlm.engine.FitEngineConfiguration;
import uk.ac.sussex.gdsc.smlm.engine.FitJob.Status;
import uk.ac.sussex.gdsc.smlm.engine.FitParameters;
import uk.ac.sussex.gdsc.smlm.engine.FitParameters.FitTask;
import uk.ac.sussex.gdsc.smlm.engine.ParameterisedFitJob;
import uk.ac.sussex.gdsc.smlm.ij.IJImageSource;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.ImageSource;
import uk.ac.sussex.gdsc.smlm.results.count.ConsecutiveFailCounter;
import uk.ac.sussex.gdsc.smlm.results.count.FailCounter;
import uk.ac.sussex.gdsc.smlm.results.count.PassRateFailCounter;
import uk.ac.sussex.gdsc.smlm.results.count.ResettingFailCounter;
import uk.ac.sussex.gdsc.smlm.results.count.RollingWindowFailCounter;
import uk.ac.sussex.gdsc.smlm.results.count.WeightedFailCounter;

import gnu.trove.list.array.TByteArrayList;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.text.TextWindow;

import java.awt.Rectangle;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicReference;
import java.util.regex.Pattern;

/**
 * This plugin handles generation and analysis of fail counts to optimise the stopping criteria for
 * a sequential analysis routine.
 */
public class FailCountManager implements PlugIn {
  private static final String TITLE = "Fail Count Manager";

  private static int maxCounters = 200000;

  private static final String[] OPTIONS =
      SettingsManager.getNames((Object[]) FailCountOption.values());

  private static TurboList<FailCountData> failCountData = new TurboList<>(1);
  private static AtomicReference<TextWindow> resultsWindowRef = new AtomicReference<>();

  private FailCountManagerSettings.Builder settings;

  //@formatter:off
  private enum FailCountOption implements NamedObject
  {
    CREATE_DATA { @Override
    public String getName() { return "Create Data"; } },
    LOAD_DATA { @Override
    public String getName() { return "Load Data"; } },
    SAVE_DATA { @Override
    public String getName() { return "Save Data"; } },
    PLOT_DATA { @Override
    public String getName() { return "Plot Data"; } },
    ANALYSE_DATA { @Override
    public String getName() { return "Analyse Data"; } },
    ;
    //@formatter:on

    @Override
    public String getShortName() {
      return getName();
    }

    public static FailCountOption forOrdinal(int ordinal) {
      final FailCountOption[] values = FailCountOption.values();
      if (ordinal < 0 || ordinal >= values.length) {
        ordinal = 0;
      }
      return values[ordinal];
    }
  }

  /**
   * Hold the fail count data for a single sequential analysis routine.
   */
  private static class FailCountData {
    /** The id of the data. */
    final int id;

    /** The results (pass./fail). */
    private final boolean[] results;

    private int maxConsFailCount = -1;

    // These are for plotting so we use float not int
    private float[] candidate;
    private float[] consFailCount;
    private float[] passCount;
    private float[] passRate;

    /**
     * The number of results to process before a fail counter is not OK. Used to score a fail
     * counter
     */
    private int target;
    private float targetPassCount;
    private double maxScore;

    FailCountData(int id, boolean[] results) {
      this.id = id;
      this.results = results;
    }

    int getMaxConsecutiveFailCount() {
      if (maxConsFailCount == -1) {
        int consFail = 0;
        int max = 0;
        for (int i = 0; i < results.length; i++) {
          if (results[i]) {
            consFail = 0;
          } else {
            consFail++;
            if (max < consFail) {
              max = consFail;
            }
          }
        }
        maxConsFailCount = max;
      }
      return maxConsFailCount;
    }

    int getPassCount() {
      createData();
      return (int) passCount[passCount.length - 1];
    }

    int getFailCount() {
      return results.length - getPassCount();
    }

    void createData() {
      if (candidate == null) {
        initialiseData();
      }
    }

    private synchronized void initialiseData() {
      if (candidate != null) {
        return;
      }
      int pass = 0;
      int consFail = 0;

      final int size = results.length;
      candidate = new float[size];
      passCount = new float[size];
      passRate = new float[size];
      consFailCount = new float[size];
      for (int i = 0; i < size; i++) {
        if (results[i]) {
          pass++;
          consFail = 0;
        } else {
          consFail++;
          // Only set this when non-zero
          consFailCount[i] = consFail;
        }
        candidate[i] = i + 1;
        passCount[i] = pass;
        passRate[i] = (float) pass / (i + 1);
      }
    }

    float[] getRollingFailCount(int rollingWindow) {
      final BooleanRollingArray c = new BooleanRollingArray(rollingWindow);
      final int size = results.length;
      final float[] failCount = new float[size];
      for (int i = 0; i < size; i++) {
        c.add(results[i]);
        failCount[i] = c.getFalseCount();
      }
      return failCount;
    }

    float[] getWeightedFailCount(int passWeight, int failWeight) {
      final WeightedFailCounter c =
          WeightedFailCounter.create(Integer.MAX_VALUE, failWeight, passWeight);
      final int size = results.length;
      final float[] failCount = new float[size];
      for (int i = 0; i < size; i++) {
        c.addResult(results[i]);
        failCount[i] = c.getFailCount();
      }
      return failCount;
    }

    float[] getResettingFailCount(double resetFraction) {
      final ResettingFailCounter c = ResettingFailCounter.create(Integer.MAX_VALUE, resetFraction);
      final int size = results.length;
      final float[] failCount = new float[size];
      for (int i = 0; i < size; i++) {
        c.addResult(results[i]);
        failCount[i] = c.getFailCount();
      }
      return failCount;
    }

    void initialiseAnalysis(double targetPassFraction) {
      initialiseData();
      targetPassCount = Math.round(passCount[passCount.length - 1] * targetPassFraction);
      final int size = results.length;
      target = 1;
      while (target <= size) {
        if (passCount[target - 1] >= targetPassCount) {
          break;
        }
        target++;
      }
      maxScore = score(size);
    }

    double score(FailCounter counter) {
      final int size = results.length;
      int index = 0;
      while (index < size) {
        if (results[index]) {
          counter.pass();
        } else {
          counter.fail();
        }
        index++;
        if (!counter.isOk()) {
          return score(index);
        }
      }
      return maxScore;
    }

    double score(int n) {
      if (n == target) {
        return 0; // Perfect
      }
      if (n < target) {
        // Penalise stopping too early.
        final float remaining = (targetPassCount - passCount[n - 1]) / targetPassCount;
        // This has a score from 0 to 1.
        // return remaining * remaining;
        return remaining;
      }

      // Penalise running too long.
      // Overrun will be above 0.
      final float overrun = ((float) (n - target)) / target;

      // This has a score from 0 to Infinity.
      // return overrun * overrun;
      return overrun;

      // This has a score from 0 to Infinity but does not heavily penalise large overrun
      // if (overrun < 1)
      // // So the gradient is 1 at x=1, f(x)=0.5*x^2, f'(x)=x
      // return overrun * overrun / 2.0;
      // else
      // return overrun;

      // This has a score from 0 to 1. This equally weights under/overrun.
      // However very long overrun is not penalised due to the exponential.
      // 0 0
      // 0.5 0.3934693403
      // 1 0.6321205588
      // 2 0.8646647168
      // 3 0.9502129316
      // 4 0.9816843611
      // 5 0.993262053
      // 10 0.9999546001
      // return 1.0 - FastMath.exp(-overrun);
    }
  }

  private enum CounterStatus {
    CONTINUE, ANALYSE, RETURN;
  }

  private static class PlotData {
    final int item;
    final int rollingWindow;
    final int passWeight;
    final int failWeight;
    final double resetFraction;
    final boolean fixedXAxis;

    PlotData(int item, boolean fixedXAxis, int rollingWindow, int passWeight, int failWeight,
        double resetFraction) {
      this.item = item;
      this.fixedXAxis = fixedXAxis;
      this.rollingWindow = rollingWindow;
      this.passWeight = passWeight;
      this.failWeight = failWeight;
      this.resetFraction = resetFraction;
    }

    /**
     * Test if this equals the other object.
     *
     * @param that the other object
     * @return true, if successful
     */
    public boolean equals(PlotData that) {
      //@formatter:off
      return that != null &&
          this.item == that.item &&
          this.fixedXAxis == that.fixedXAxis &&
          this.rollingWindow == that.rollingWindow &&
          this.passWeight == that.passWeight &&
          this.failWeight == that.failWeight &&
          this.resetFraction == that.resetFraction
          ;
      //@formatter:on
    }

    /**
     * Checks if is a new item.
     *
     * @param that the other object
     * @return true, if is new item
     */
    public boolean isNewItem(PlotData that) {
      if (that == null) {
        return true;
      }
      return this.item != that.item;
    }
  }

  private static class PlotWorker implements Runnable {
    final ConcurrentMonoStack<PlotData> stack;
    final TurboList<FailCountData> failCountData;
    PlotData lastPlotData;
    int maxSize;

    PlotWorker(ConcurrentMonoStack<PlotData> stack, TurboList<FailCountData> failCountData) {
      this.stack = stack;
      this.failCountData = failCountData;
      for (int i = 0; i < failCountData.size(); i++) {
        maxSize = Math.max(maxSize, failCountData.getf(i).results.length);
      }
    }

    @Override
    public void run() {
      // while (!Thread.interrupted())
      while (true) {
        try {
          final PlotData plotData = stack.pop();
          if (plotData == null) {
            break;
          }
          if (plotData.equals(lastPlotData)) {
            continue;
          }
          run(plotData);
          lastPlotData = plotData;
        } catch (final InterruptedException ex) {
          // Thread.currentThread().interrupt();
          break;
        }
      }
    }

    private void run(PlotData plotData) {
      final boolean isNew =
          plotData.isNewItem(lastPlotData) || plotData.fixedXAxis != lastPlotData.fixedXAxis;
      final int item = plotData.item - 1; // 0-based index
      if (item < 0 || item >= failCountData.size()) {
        return;
      }
      final FailCountData data = failCountData.get(item);

      // Ensure consistent synchronisation
      data.createData();
      final WindowOrganiser wo = new WindowOrganiser();
      if (isNew) {
        display(wo, "Pass Count", data.candidate, data.passCount, plotData.fixedXAxis);
        display(wo, "Pass Rate", data.candidate, data.passRate, plotData.fixedXAxis);
        display(wo, "Consecutive Fail Count", data.candidate, data.consFailCount,
            plotData.fixedXAxis);
      }

      // Only rebuild if changed
      if (isNew || plotData.rollingWindow != lastPlotData.rollingWindow) {
        display(wo, "Rolling Fail Count", data.candidate,
            data.getRollingFailCount(plotData.rollingWindow), plotData.fixedXAxis);
      }
      if (isNew || plotData.passWeight != lastPlotData.passWeight
          || plotData.failWeight != lastPlotData.failWeight) {
        display(wo, "Weighted Fail Count", data.candidate,
            data.getWeightedFailCount(plotData.passWeight, plotData.failWeight),
            plotData.fixedXAxis);
      }
      if (isNew || plotData.resetFraction != lastPlotData.resetFraction) {
        display(wo, "Resetting Fail Count", data.candidate,
            data.getResettingFailCount(plotData.resetFraction), plotData.fixedXAxis);
      }

      wo.tile();
    }

    private void display(WindowOrganiser wo, String string, float[] x, float[] y,
        boolean fixedXAxis) {
      final String title = TITLE + " " + string;
      final Plot plot = new Plot(title, "Candidate", string);
      final double maxx = (fixedXAxis) ? maxSize : x[x.length - 1];
      final double maxy = MathUtils.max(y);
      plot.setLimits(x[0], maxx, 0, maxy * 1.05);
      plot.addPoints(x, y, Plot.LINE);
      plot.addLabel(0, 0, "Max = " + maxy);
      ImageJUtils.display(title, plot, 0, wo);
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    settings = SettingsManager.readFailCountManagerSettings(0).toBuilder();

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addChoice("Option", OPTIONS, settings.getOption());
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    settings.setOption(gd.getNextChoiceIndex());

    final FailCountOption option = FailCountOption.forOrdinal(settings.getOption());
    switch (option) {
      case CREATE_DATA:
        createData();
        break;
      case LOAD_DATA:
        loadData();
        break;
      case SAVE_DATA:
        saveData();
        break;
      case PLOT_DATA:
        plotData();
        break;
      case ANALYSE_DATA:
        analyseData();
        break;
      default:
        throw new IllegalStateException("Unknown option: " + option);
    }

    SettingsManager.writeSettings(settings);
  }

  /**
   * Creates the fail count data by running fitting on the current image.
   */
  private void createData() {
    final ImagePlus imp = WindowManager.getCurrentImage();
    if (imp == null) {
      IJ.error(TITLE, "No image for fitting");
      return;
    }

    if (!showCreateDataDialog(imp)) {
      return;
    }

    // Get the current fit configuration
    final Configuration c = new Configuration();
    if (!c.showDialog(false)) {
      return;
    }

    final FitEngineConfiguration fitConfig = c.getFitEngineConfiguration();
    // Update stopping criteria.
    fitConfig.resetFailCounter();
    fitConfig.setFailuresLimit(settings.getFailCountLimit());

    final ImageSource source = new IJImageSource(imp);
    final PeakFit peakFit = new PeakFit(fitConfig, ResultsSettings.getDefaultInstance());
    peakFit.setResultsSuffix("(FailCountAnalysis)");
    if (!peakFit.initialise(source, null, false)) {
      IJ.error(TITLE, "Failed to initialise the fit engine");
      return;
    }
    final FitEngine engine = peakFit.createFitEngine();

    final Rectangle bounds = new Rectangle(source.getWidth(), source.getHeight());

    // Run
    final int totalFrames = Math.min(source.getFrames(), settings.getMaxFrames());
    final int step = ImageJUtils.getProgressInterval(totalFrames);
    IJ.showProgress(0);
    boolean shutdown = false;
    int slice = 0;
    final TurboList<ParameterisedFitJob> jobs = new TurboList<>(totalFrames);
    while (!shutdown && slice < totalFrames) {
      final float[] data = source.next();
      if (data == null) {
        break;
      }

      if (slice++ % step == 0) {
        if (ImageJUtils.showStatus("Fitting slice: " + slice + " / " + totalFrames)) {
          IJ.showProgress(slice, totalFrames);
        }
      }

      final ParameterisedFitJob job = createJob(source.getStartFrameNumber(), data, bounds);
      jobs.addf(job);
      engine.run(job);

      shutdown = escapePressed();
    }

    ImageJUtils.showStatus("Extracting fail count data");
    engine.end(shutdown);
    IJ.showProgress(1);
    source.close();

    // Extract the fail count data
    final TurboList<FailCountData> failCountData = new TurboList<>(jobs.size());
    for (int i = 0; i < jobs.size(); i++) {
      final ParameterisedFitJob job = jobs.getf(i);
      if (job.getStatus() == Status.FINISHED) {
        final FitParameters fitParams = job.getFitParameters();

        // Find the last success
        boolean[] results = fitParams.pass;
        int end = results.length - 1;
        while (end > 0 && !results[end]) {
          end--;
        }
        // Add on the configured fail count limit
        end = Math.min(end + 1 + settings.getFailCountLimit(), results.length);
        results = Arrays.copyOf(results, end);
        failCountData.add(new FailCountData(job.getSlice(), results));
      }
    }
    FailCountManager.failCountData = failCountData;
    ImageJUtils.showStatus("");

    // Save for the future
    if (settings.getSaveAfterFitting()) {
      saveData();
    }
  }

  private boolean showCreateDataDialog(ImagePlus imp) {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addMessage(TextUtils.wrap("Run the fit engine on the current image to generate "
        + "pass/fail data for sequential candidates in each frame. A second dialog will "
        + "be shown to check the fit settings.", 80));
    gd.addSlider("Max_frames", 1, imp.getStackSize(), settings.getMaxFrames());
    gd.addNumericField("Fail_count_limit", settings.getFailCountLimit(), 0);
    gd.addCheckbox("Save", settings.getSaveAfterFitting());
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }
    settings.setMaxFrames((int) gd.getNextNumber());
    settings.setFailCountLimit((int) gd.getNextNumber());
    settings.setSaveAfterFitting(gd.getNextBoolean());
    try {
      ParameterUtils.isAboveZero("Max frames", settings.getMaxFrames());
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }
    return true;
  }

  private static ParameterisedFitJob createJob(int startFrame, float[] data, Rectangle bounds) {
    final FitParameters fitParams = new FitParameters();
    fitParams.fitTask = FitTask.PSF_FITTING;
    // Signal that the fail count should be recorded
    fitParams.pass = new boolean[0];
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
   * Load the data from a file.
   */
  private void loadData() {
    final String filename =
        ImageJUtils.getFilename("Fail_count_data_filename", settings.getFilename());
    if (filename == null) {
      return;
    }
    settings.setFilename(filename);
    final TurboList<FailCountData> countData = new TurboList<>();

    try (BufferedReader br = Files.newBufferedReader(Paths.get(filename))) {
      final Pattern pattern = Pattern.compile("[\t, ]+");
      // Ignore the first line
      String line = br.readLine();
      final BooleanArray array = new BooleanArray(100);
      int lastId = 0;
      int lastCandidate = 0;
      while ((line = br.readLine()) != null) {
        final String[] data = pattern.split(line);
        if (data.length != 3) {
          throw new IOException("Require 3 fields in the data");
        }

        final int id = Integer.parseInt(data[0]);
        if (id < 1) {
          throw new IOException("ID must be strictly positive");
        }
        final int candidate = Integer.parseInt(data[1]);
        if (candidate < 1) {
          throw new IOException("Candidate must be strictly positive");
        }
        final boolean ok = guessStatus(data[2]);

        if (lastId != id) {
          if (array.size() > 0) {
            countData.add(new FailCountData(lastId, array.toArray()));
            array.clear();
          }
          if (candidate != 1) {
            throw new IOException("Candidate must start at 1");
          }
          lastId = id;
          lastCandidate = candidate - 1; // Ensure continuous
        }
        // Require continuous sequence
        if (candidate - lastCandidate == 1) {
          array.add(ok);
          lastCandidate = candidate;
        } else {
          // Make impossible to add any more for this ID
          lastCandidate = -1;
        }
      }
      // Final ID
      if (array.size() > 0) {
        countData.add(new FailCountData(lastId, array.toArray()));
      }

      IJ.showMessage(TITLE, "Loaded " + TextUtils.pleural(countData.size(), "sequence"));
      FailCountManager.failCountData = countData;
    } catch (final NumberFormatException | IOException ex) {
      IJ.error(TITLE, "Failed to load data:\n" + ex.getMessage());
    }
  }

  /**
   * Guess the pass/fail status.
   *
   * @param string the string
   * @return true, if successful
   * @throws IOException Signals that an I/O exception has occurred.
   */
  private static boolean guessStatus(String string) throws IOException {
    final int len = string.length();
    if (len < 1) {
      return false;
    }

    if (len == 1) {
      final char c = string.charAt(0);
      if (c == 'y' || c == '1') {
        return true;
      }
      if (c == 'n' || c == '0') {
        return false;
      }
    } else {
      string = string.toLowerCase(Locale.US);
      if (string.equals("pass")) {
        return true;
      }
      if (string.equals("fail")) {
        return false;
      }
      if (string.equals("ok")) {
        return true;
      }
    }
    throw new IOException("Unrecognised status: " + string);
  }

  /**
   * Save the data in memory to file.
   */
  private void saveData() {
    final TurboList<FailCountData> failCountData = FailCountManager.failCountData;
    if (failCountData.isEmpty()) {
      IJ.error(TITLE, "No fail count data in memory");
      return;
    }
    final String filename =
        ImageJUtils.getFilename("Fail_count_data_filename", settings.getFilename());
    if (filename == null) {
      return;
    }
    settings.setFilename(filename);

    try (BufferedWriter bw = Files.newBufferedWriter(Paths.get(filename))) {
      bw.write("ID,Candidate,Status");
      bw.newLine();
      for (int i = 0; i < failCountData.size(); i++) {
        final FailCountData d = failCountData.get(i);
        final String prefix = d.id + ",";
        final boolean[] pass = d.results;
        for (int j = 0; j < pass.length; j++) {
          bw.write(prefix);
          bw.write(Integer.toString(j + 1));
          bw.write(',');
          if (pass[j]) {
            bw.write('y');
          } else {
            bw.write('n');
          }
          bw.newLine();
        }
      }
    } catch (final IOException ex) {
      IJ.error(TITLE, "Failed to save data:\n" + ex.getMessage());
    }
  }

  /**
   * Show an interactive plot of the fail count data.
   */
  private void plotData() {
    final TurboList<FailCountData> failCountData = FailCountManager.failCountData;
    if (failCountData.isEmpty()) {
      IJ.error(TITLE, "No fail count data in memory");
      return;
    }

    // Find max fail count size
    final int max = getMaxConsecutiveFailCount(failCountData);

    final ConcurrentMonoStack<PlotData> stack = new ConcurrentMonoStack<>();
    new Thread(new PlotWorker(stack, failCountData)).start();

    final NonBlockingExtendedGenericDialog gd = new NonBlockingExtendedGenericDialog(TITLE);
    gd.addSlider("Item", 1, failCountData.size(), settings.getPlotItem());
    gd.addCheckbox("Fixed_x_axis", settings.getPlotFixedXAxis());
    gd.addMessage("Rolling Window Fail Count");
    gd.addSlider("Rolling_window", 1, 3 * max, settings.getPlotRollingWindow());
    gd.addMessage("Weighted Fail Count");
    gd.addSlider("Pass_weight", 1, 20, settings.getPlotPassWeight());
    gd.addSlider("Fail_weight", 1, 20, settings.getPlotFailWeight());
    gd.addMessage("Resetting Fail Count");
    gd.addSlider("Reset_fraction", 0.05, 0.95, settings.getPlotResetFraction());
    gd.addDialogListener((gd2, event) -> {
      final int item = (int) gd2.getNextNumber();
      final boolean fixedXAxis = gd2.getNextBoolean();
      final int rollingWindow = (int) gd2.getNextNumber();
      final int passWeight = (int) gd2.getNextNumber();
      final int failWeight = (int) gd2.getNextNumber();
      final double resetFraction = gd2.getNextNumber();
      settings.setPlotItem(item);
      settings.setPlotRollingWindow(rollingWindow);
      settings.setPlotPassWeight(passWeight);
      settings.setPlotFailWeight(failWeight);
      settings.setPlotResetFraction(resetFraction);
      stack.insert(
          new PlotData(item, fixedXAxis, rollingWindow, passWeight, failWeight, resetFraction));
      return true;
    });

    gd.hideCancelButton();
    gd.setOKLabel("Close");
    gd.showDialog();

    stack.close(gd.wasCanceled());
  }

  private static int getMaxConsecutiveFailCount(TurboList<FailCountData> failCountData) {
    int max = 1;
    for (int i = 0; i < failCountData.size(); i++) {
      max = Math.max(max, failCountData.getf(i).getMaxConsecutiveFailCount());
    }
    return max;
  }

  private static int getMaxFailCount(TurboList<FailCountData> failCountData) {
    int max = 1;
    for (int i = 0; i < failCountData.size(); i++) {
      max = Math.max(max, failCountData.getf(i).getFailCount());
    }
    return max;
  }

  @SuppressWarnings("unused")
  private static int getMaxPassCount(TurboList<FailCountData> failCountData) {
    int max = 1;
    for (int i = 0; i < failCountData.size(); i++) {
      max = Math.max(max, failCountData.getf(i).getPassCount());
    }
    return max;
  }

  private void analyseData() {
    final TurboList<FailCountData> failCountData = FailCountManager.failCountData;
    if (failCountData.isEmpty()) {
      IJ.error(TITLE, "No fail count data in memory");
      return;
    }

    if (!showAnalysisDialog()) {
      return;
    }

    final int maxCons = getMaxConsecutiveFailCount(failCountData);
    final int maxFail = getMaxFailCount(failCountData);
    // final int maxPass = getMaxPassCount(failCountData);

    // Create a set of fail counters
    final TurboList<FailCounter> counters = new TurboList<>();
    final TByteArrayList type = new TByteArrayList();
    for (int i = 0; i <= maxCons; i++) {
      counters.add(ConsecutiveFailCounter.create(i));
    }
    fill(type, counters, 0);

    // The other counters are user configured.
    // Ideally this would be a search to optimise the best parameters
    // for each counter as any enumeration may be way off the mark.

    // Note that 0 failures in a window can be scored using the consecutive fail counter.
    int max = Math.min(maxFail, settings.getRollingCounterMaxAllowedFailures());
    for (int fail = MathUtils.clip(1, maxFail, settings.getRollingCounterMinAllowedFailures());
        fail <= max; fail++) {
      // Note that n-1 failures in window n can be scored using the consecutive fail counter.
      for (int window = Math.max(fail + 2, settings.getRollingCounterMinWindow());
          window <= settings.getRollingCounterMaxWindow(); window++) {
        counters.add(RollingWindowFailCounter.create(fail, window));
      }
      switch (checkCounters(counters)) {
        case ANALYSE:
          break;
        case CONTINUE:
          break;
        case RETURN:
          return;
        default:
          throw new IllegalStateException();
      }
    }
    fill(type, counters, 1);

    max = Math.min(maxFail, settings.getWeightedCounterMaxAllowedFailures());
    for (int fail = MathUtils.min(maxFail, settings.getWeightedCounterMinAllowedFailures());
        fail <= max; fail++) {
      for (int w = settings.getWeightedCounterMinPassDecrement();
          w <= settings.getWeightedCounterMaxPassDecrement(); w++) {
        counters.add(WeightedFailCounter.create(fail, 1, w));
      }
      switch (checkCounters(counters)) {
        case ANALYSE:
          break;
        case CONTINUE:
          break;
        case RETURN:
          return;
        default:
          throw new IllegalStateException();
      }
    }
    fill(type, counters, 2);

    max = Math.min(maxFail, settings.getResettingCounterMaxAllowedFailures());
    for (int fail = MathUtils.min(maxFail, settings.getResettingCounterMinAllowedFailures());
        fail <= max; fail++) {
      for (double f = settings.getResettingCounterMinResetFraction();
          f <= settings.getResettingCounterMaxResetFraction();
          f += settings.getResettingCounterIncResetFraction()) {
        counters.add(ResettingFailCounter.create(fail, f));
      }
      switch (checkCounters(counters)) {
        case ANALYSE:
          break;
        case CONTINUE:
          break;
        case RETURN:
          return;
        default:
          throw new IllegalStateException();
      }
    }
    fill(type, counters, 3);

    for (int count = settings.getPassRateCounterMinAllowedCounts();
        count <= settings.getPassRateCounterMaxAllowedCounts(); count++) {
      for (double f = settings.getPassRateCounterMinPassRate();
          f <= settings.getPassRateCounterMaxPassRate();
          f += settings.getPassRateCounterIncPassRate()) {
        counters.add(PassRateFailCounter.create(count, f));
      }
      switch (checkCounters(counters)) {
        case ANALYSE:
          break;
        case CONTINUE:
          break;
        case RETURN:
          return;
        default:
          throw new IllegalStateException();
      }
    }
    fill(type, counters, 4);

    counters.trimToSize();

    // Score each of a set of standard fail counters against each frame using how
    // close they are to the target.
    final double[] score = new double[counters.size()];
    final double targetPassFraction = settings.getTargetPassFraction();

    final int nThreads = Prefs.getThreads();
    final ExecutorService executor = Executors.newFixedThreadPool(nThreads);
    final TurboList<Future<?>> futures = new TurboList<>(nThreads);

    final Ticker ticker = ImageJUtils.createTicker(failCountData.size(), nThreads);
    IJ.showStatus("Analysing " + TextUtils.pleural(counters.size(), "counter"));
    for (int i = 0; i < failCountData.size(); i++) {
      final FailCountData data = failCountData.getf(i);
      futures.add(executor.submit(() -> {
        if (IJ.escapePressed()) {
          return;
        }

        // TODO - Ideally this plugin should be run on benchmark data with ground truth.
        // The target could be to ensure all all the correct results are fit
        // and false positives are excluded from incrementing the pass counter.
        // This could be done by saving the results from a benchmarking scoring
        // plugin to memory as the current dataset.
        data.initialiseAnalysis(targetPassFraction);

        // Score in blocks and then do a synchronized write to the combined score
        final Thread t = Thread.currentThread();
        final double[] s = new double[8192];
        int index = 0;
        while (index < counters.size()) {
          if (t.isInterrupted()) {
            break;
          }
          final int block = Math.min(8192, counters.size() - index);
          for (int j = 0; j < block; j++) {
            final FailCounter counter = counters.getf(index + j).newCounter();
            s[j] = data.score(counter);
          }
          // Write to the combined score
          synchronized (score) {
            for (int j = 0; j < block; j++) {
              score[index + j] += s[j];
            }
          }
          index += block;
        }
        ticker.tick();
      }));
    }

    ConcurrencyUtils.waitForCompletionUnchecked(futures);
    executor.shutdown();
    IJ.showProgress(1);
    if (IJ.escapePressed()) {
      IJ.showStatus("");
      IJ.error(TITLE, "Cancelled analysis");
      return;
    }
    IJ.showStatus("Summarising results ...");

    // TODO - check if the top filter is at the bounds of the range
    final int minIndex = SimpleArrayUtils.findMinIndex(score);
    ImageJUtils.log(TITLE + " Analysis : Best counter = %s (Score = %f)",
        counters.getf(minIndex).getDescription(), score[minIndex]);

    // Show a table of results for the top N for each type
    final int topN = Math.min(settings.getTableTopN(), score.length);
    if (topN > 0) {
      final byte[] types = type.toArray();
      final byte maxType = types[types.length - 1];
      final TextWindow resultsWindow = createTable();
      for (byte b = 0; b <= maxType; b++) {
        int[] indices;
        // Use a heap to avoid a full sort
        final IntResultHeap heap = new IntResultHeap(topN);
        for (int i = 0; i < score.length; i++) {
          if (types[i] == b) {
            heap.addValue(score[i], i);
          }
        }
        if (heap.getSize() == 0) {
          continue;
        }
        indices = heap.getData();
        // Ensure sorted
        SortUtils.sortIndices(indices, score, false);

        final StringBuilder sb = new StringBuilder();
        try (final BufferedTextWindow tw = new BufferedTextWindow(resultsWindow)) {
          for (int i = 0; i < topN; i++) {
            sb.setLength(0);
            final int j = indices[i];
            sb.append(i + 1).append('\t');
            sb.append(counters.getf(j).getDescription()).append('\t');
            sb.append(score[j]);
            tw.append(sb.toString());
          }
        }
      }
    }

    // TODO - Save the best fail counter to the current fit configuration.

    IJ.showStatus("");
  }

  private static void fill(TByteArrayList type, TurboList<FailCounter> counters, int value) {
    final int n = counters.size() - type.size();
    ImageJUtils.log("Type %d = %d", value, n);
    type.fill(type.size(), counters.size(), (byte) value);
  }

  private static CounterStatus checkCounters(TurboList<FailCounter> counters) {
    if (counters.size() > maxCounters) {
      final GenericDialog gd = new GenericDialog(TITLE);
      gd.addMessage("Too many counters to analyse: " + counters.size());
      gd.addNumericField("Max_counters", maxCounters, 0);
      gd.enableYesNoCancel(" Analyse ", " Continue ");
      gd.showDialog();
      if (gd.wasCanceled()) {
        return CounterStatus.RETURN;
      }
      if (gd.wasOKed()) {
        return CounterStatus.ANALYSE;
      }
      final int newMaxCounters = (int) gd.getNextNumber();
      if (newMaxCounters <= maxCounters) {
        IJ.error(TITLE, "The max counters has not been increased, unable to continue");
        return CounterStatus.RETURN;
      }
      maxCounters = newMaxCounters;
    }
    return CounterStatus.CONTINUE;
  }

  private static TextWindow createTable() {
    return ImageJUtils.refresh(resultsWindowRef, () -> new TextWindow(TITLE + " Analysis Results",
        "Rank\tFail Counter\tScore", "", 600, 400));
  }

  private boolean showAnalysisDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addMessage(
        TextUtils.wrap("Analysis a set of fail counters on the current pass/fail data.", 80));
    gd.addSlider("Target_pass_fraction", 0.1, 1, settings.getTargetPassFraction());
    gd.addSliderIncludeDefault("Table_top_n", 0, 100, settings.getTableTopN());
    gd.addNumericField("Rolling_counter_min_allowed_failures",
        settings.getRollingCounterMinAllowedFailures(), 0);
    gd.addNumericField("Rolling_counter_max_allowed_failures",
        settings.getRollingCounterMaxAllowedFailures(), 0);
    gd.addNumericField("Rolling_counter_min_window", settings.getRollingCounterMinWindow(), 0);
    gd.addNumericField("Rolling_counter_max_window", settings.getRollingCounterMaxWindow(), 0);
    gd.addNumericField("Weighted_counter_min_allowed_failures",
        settings.getWeightedCounterMinAllowedFailures(), 0);
    gd.addNumericField("Weighted_counter_max_allowed_failures",
        settings.getWeightedCounterMaxAllowedFailures(), 0);
    gd.addNumericField("Weighted_counter_min_pass_decrement",
        settings.getWeightedCounterMinPassDecrement(), 0);
    gd.addNumericField("Weighted_counter_max_pass_decrement",
        settings.getWeightedCounterMaxPassDecrement(), 0);
    gd.addNumericField("Resetting_counter_min_allowed_failures",
        settings.getResettingCounterMinAllowedFailures(), 0);
    gd.addNumericField("Resetting_counter_max_allowed_failures",
        settings.getResettingCounterMaxAllowedFailures(), 0);
    gd.addNumericField("Resetting_counter_min_pass_decrement",
        settings.getResettingCounterMinResetFraction(), 2);
    gd.addNumericField("Resetting_counter_max_pass_decrement",
        settings.getResettingCounterMaxResetFraction(), 2);
    gd.addNumericField("Resetting_counter_inc_pass_decrement",
        settings.getResettingCounterIncResetFraction(), 2);
    gd.addNumericField("Pass_rate_counter_min_allowed_failures",
        settings.getPassRateCounterMinAllowedCounts(), 0);
    gd.addNumericField("Pass_rate_counter_max_allowed_failures",
        settings.getPassRateCounterMaxAllowedCounts(), 0);
    gd.addNumericField("Pass_rate_counter_min_pass_rate", settings.getPassRateCounterMinPassRate(),
        3);
    gd.addNumericField("Pass_rate_counter_max_pass_rate", settings.getPassRateCounterMaxPassRate(),
        3);
    gd.addNumericField("Pass_rate_counter_inc_pass_rate", settings.getPassRateCounterIncPassRate(),
        3);
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }
    settings.setTargetPassFraction(gd.getNextNumber());
    settings.setTableTopN((int) gd.getNextNumber());
    settings.setRollingCounterMinAllowedFailures((int) gd.getNextNumber());
    settings.setRollingCounterMaxAllowedFailures((int) gd.getNextNumber());
    settings.setRollingCounterMinWindow((int) gd.getNextNumber());
    settings.setRollingCounterMaxWindow((int) gd.getNextNumber());
    settings.setWeightedCounterMinAllowedFailures((int) gd.getNextNumber());
    settings.setWeightedCounterMaxAllowedFailures((int) gd.getNextNumber());
    settings.setWeightedCounterMinPassDecrement((int) gd.getNextNumber());
    settings.setWeightedCounterMaxPassDecrement((int) gd.getNextNumber());
    settings.setResettingCounterMinAllowedFailures((int) gd.getNextNumber());
    settings.setResettingCounterMaxAllowedFailures((int) gd.getNextNumber());
    settings.setResettingCounterMinResetFraction(gd.getNextNumber());
    settings.setResettingCounterMaxResetFraction(gd.getNextNumber());
    settings.setResettingCounterIncResetFraction(gd.getNextNumber());
    settings.setPassRateCounterMinAllowedCounts((int) gd.getNextNumber());
    settings.setPassRateCounterMaxAllowedCounts((int) gd.getNextNumber());
    settings.setPassRateCounterMinPassRate(gd.getNextNumber());
    settings.setPassRateCounterMaxPassRate(gd.getNextNumber());
    settings.setPassRateCounterIncPassRate(gd.getNextNumber());
    try {
      ParameterUtils.isAboveZero("Target pass fraction", settings.getTargetPassFraction());
      ParameterUtils.isAboveZero("Resetting counter inc pass decrement",
          settings.getResettingCounterIncResetFraction());
      ParameterUtils.isAboveZero("Pass rate counter inc pass rate",
          settings.getPassRateCounterIncPassRate());
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }
    return true;
  }
}
