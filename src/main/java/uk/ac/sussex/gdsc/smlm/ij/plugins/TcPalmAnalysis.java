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

import gnu.trove.map.hash.TIntIntHashMap;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.text.TextPanel;
import ij.text.TextWindow;
import java.awt.AWTEvent;
import java.awt.Color;
import java.awt.Rectangle;
import java.awt.Window;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.geom.Rectangle2D;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.BiPredicate;
import java.util.function.UnaryOperator;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.tuple.Pair;
import uk.ac.sussex.gdsc.core.ij.BufferedTextWindow;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.NonBlockingExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper.LutColour;
import uk.ac.sussex.gdsc.core.match.HopcroftKarpMatching.IntBiConsumer;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.SoftLock;
import uk.ac.sussex.gdsc.core.utils.SortUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.concurrent.ConcurrentMonoStack;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.TcPalmAnalysisSettings;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageSettings;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsSettings;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.results.ImageJImagePeakResults;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.ij.utils.TextPanelMouseListener;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.PeakResultsList;
import uk.ac.sussex.gdsc.smlm.results.count.FrameCounter;
import uk.ac.sussex.gdsc.smlm.results.sort.IdFramePeakResultComparator;

/**
 * Analyses the time-correlated activation of traced localisation data.
 */
public class TcPalmAnalysis implements PlugIn {

  /** The plugin title. */
  private static final String TITLE = "TC PALM Analysis";

  /** Text window showing the current clusters. */
  private static AtomicReference<CurrentClustersTextWindow> resultsWindowRef =
      new AtomicReference<>();

  /**
   * The instance lock used to prevent multiple instances running. This is because the output
   * results table and plots are reused and multiple instances cannot share the same output.
   */
  private static SoftLock instanceLock = new SoftLock();

  /** The plugin settings. */
  private TcPalmAnalysisSettings.Builder settings;

  /** The results. */
  private MemoryPeakResults results;

  /** The image. */
  private ImageJImagePeakResults image;

  /** The lock (must be held when processing the work queue). */
  private SoftLock lock;

  /** The work queue. */
  private ConcurrentMonoStack<Pair<Rectangle, TcPalmAnalysisSettings>> workQueue;

  /** The previous work. */
  private Pair<Rectangle, TcPalmAnalysisSettings> previous;

  /** The executor. */
  private ExecutorService executor;

  /** The cluster data. */
  private LocalList<ClusterData> clusterData;

  /** The minimum time frame. */
  private int minT;

  /** The maximum time frame. */
  private int maxT;

  /** The cluster selected listener. */
  private ClusterSelectedListener clusterSelectedListener;

  /**
   * The total activations plot data from the currently clusters. ALlows adding more lines to the
   * plot when clusters are selected.
   */
  private ActivationsPlotData activationsPlotData;

  /**
   * Store data on each cluster.
   */
  private static class ClusterData {
    final int id;
    int index;
    LocalList<PeakResult> results;
    float x;
    float y;
    float width;
    float height;
    float[] frames;
    int start;
    int end;

    ClusterData(PeakResult result) {
      id = result.getId();
      x = width = result.getXPosition();
      y = height = result.getYPosition();
      results = new LocalList<>();
      results.add(result);
    }

    void add(PeakResult result) {
      float value = result.getXPosition();
      if (value < x) {
        x = value;
      } else if (value > width) {
        width = value;
      }
      value = result.getYPosition();
      if (value < y) {
        y = value;
      } else if (value > height) {
        height = value;
      }
      results.add(result);
    }

    void finalise() {
      results.trimToSize();
      width -= x;
      height -= y;
      start = results.unsafeGet(0).getFrame();
      end = results.unsafeGet(results.size() - 1).getFrame();
    }

    /**
     * True if all points are within the rectangle.
     *
     * @param r the rectangle
     * @return true if within
     */
    boolean isWithin(Rectangle2D r) {
      return r.contains(x, y, width, height);
    }

    /**
     * True if any points are within the rectangle.
     *
     * @param r the rectangle
     * @return true if within
     */
    boolean intersects(Rectangle2D r) {
      if (r.intersects(x, y, width, height)) {
        final int size = results.size();
        for (int i = 0; i < size; i++) {
          final PeakResult result = results.unsafeGet(i);
          if (r.contains(result.getXPosition(), result.getYPosition())) {
            return true;
          }
        }
      }
      return false;
    }

    /**
     * Gets the frames as a float[] for convenience when plotting.
     *
     * @return the frames
     */
    float[] getFrames() {
      float[] frames = this.frames;
      if (frames == null) {
        frames = new float[results.size()];
        final int size = results.size();
        for (int i = 0; i < size; i++) {
          frames[i] = results.unsafeGet(i).getFrame();
        }
        this.frames = frames;
      }
      return frames;
    }

    double[] getXyz() {
      final double[] xyz = new double[3];
      results.forEach(r -> {
        xyz[0] += r.getXPosition();
        xyz[1] += r.getYPosition();
        xyz[2] += r.getZPosition();
      });
      final int size = results.size();
      for (int i = 0; i < 3; i++) {
        xyz[i] /= size;
      }
      return xyz;
    }
  }

  /**
   * Class to generate a cumulative count for plotting.
   */
  private static class CumulativeCount {
    private float[] count;

    CumulativeCount() {
      count = SimpleArrayUtils.newArray(20, 1f, 1f);
    }

    float[] getCount(int size) {
      float[] count = this.count;
      if (count.length < size) {
        count = this.count = SimpleArrayUtils.newArray(size, 1f, 1f);
      }
      return Arrays.copyOf(count, size);
    }
  }

  /**
   * Custom TextWindow to allow access to the associated TextPanelMouseListener.
   */
  private static class CurrentClustersTextWindow extends TextWindow {
    private static final long serialVersionUID = 1L;

    /**
     * Class to allow the events raised by the single instance of the current clusters table to pass
     * to a custom handler.
     */
    private static class ClusterSelectedListenerHolder extends TextPanelMouseListener {
      static final IntBiConsumer NOOP = (u, v) -> {
        // noop
      };

      IntBiConsumer action = NOOP;

      ClusterSelectedListenerHolder(TextPanel textPanel) {
        super(textPanel);
      }

      @Override
      public void selected(int selectionStart, int selectionEnd) {
        action.accept(selectionStart, selectionEnd);
      }
    }

    private final ClusterSelectedListenerHolder listener;

    /**
     * Create a new instance.
     *
     * @param title the title
     * @param headings the headings
     * @param width the width
     * @param height the height
     */
    public CurrentClustersTextWindow(String title, String headings, int width, int height) {
      super(title, headings, "", width, height);
      // Only one instance can bind to the text panel as it cannot be removed.
      // So instead allow a swap of the action performed by the listener.
      listener = new ClusterSelectedListenerHolder(getTextPanel());
    }

    /**
     * Sets the action for the TextPanelMouseListener.
     *
     * @param action the new action
     */
    void setAction(ClusterSelectedListener action) {
      if (action == null) {
        listener.action = ClusterSelectedListenerHolder.NOOP;
      } else {
        action.textPanel = getTextPanel();
        listener.action = action;
      }
    }
  }

  /**
   * Class to display data from the selected clusters.
   */
  private class ClusterSelectedListener implements IntBiConsumer {
    TextPanel textPanel;

    @Override
    public void accept(int selectionStart, int selectionEnd) {
      // Overlay the clusters on the image.
      final Overlay overlay = new Overlay();
      // For each cluster create a track and add to the overlay.
      int start = maxT;
      int end = minT;
      for (int i = selectionStart; i <= selectionEnd; i++) {
        final String line = textPanel.getLine(i);
        final int indexOf = line.indexOf('\t');
        final int index = Integer.parseInt(line.substring(0, indexOf));
        final ClusterData data = clusterData.get(index);
        appendRoi(overlay, data);
        start = Math.min(data.start, start);
        end = Math.max(data.end, end);
      }
      image.getImagePlus().setOverlay(overlay);

      // Draw an activation plot of the selected clusters with analysis of bursts.
      final ActivationsPlotData activationsPlotData = TcPalmAnalysis.this.activationsPlotData;
      // Find the start and end
      final int is = ArrayUtils.indexOf(activationsPlotData.frames, start);
      final int ie = ArrayUtils.indexOf(activationsPlotData.frames, end);
      // Extract the lines.
      final int[] frames = Arrays.copyOfRange(activationsPlotData.frames, is, ie + 1);
      final int[] counts = Arrays.copyOfRange(activationsPlotData.counts, is, ie + 1);

      // Add to the plot.
      final Plot plot2 = activationsPlotData.plot;
      plot2.restorePlotObjects();
      plot2.setColor(Color.red);
      plot2.addPoints(activationsPlotData.timeConverter.apply(SimpleArrayUtils.toFloat(frames)),
          SimpleArrayUtils.toFloat(counts), Plot.BAR);
      plot2.updateImage();
    }

    private void appendRoi(Overlay overlay, ClusterData data) {
      final int np = data.results.size();
      final float[] xp = new float[np];
      final float[] yp = new float[np];
      for (int i = 0; i < np; i++) {
        final PeakResult r = data.results.unsafeGet(i);
        xp[i] = image.mapX(r.getXPosition());
        yp[i] = image.mapY(r.getYPosition());
      }
      final PointRoi roi = new PointRoi(xp, yp, np);
      roi.setShowLabels(false);
      overlay.add(roi);
    }
  }

  /**
   * Class containing the data for the most recent total activations plot for current clusters.
   */
  private static class ActivationsPlotData {
    Plot plot;
    UnaryOperator<float[]> timeConverter;
    int[] frames;
    int[] counts;

    ActivationsPlotData(Plot plot, UnaryOperator<float[]> timeConverter, int[] frames,
        int[] counts) {
      this.plot = plot;
      this.timeConverter = timeConverter;
      this.frames = frames;
      this.counts = counts;
      plot.savePlotObjects();
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "No localisations in memory");
      return;
    }

    // Only allow 1 instance to run
    if (!instanceLock.acquire()) {
      final Window w = WindowManager.getWindow(TITLE);
      if (w != null) {
        w.toFront();
        return;
      }
      // Fall through to allow the plugin to run. This may still have concurrency issues if
      // another version is running but the window is not currently showing/registered with
      // the window manager. Perhaps show a dialog asking to continue.
    }

    if (!showDialog()) {
      instanceLock.release();
      return;
    }

    // Load the results
    results = ResultsManager.loadInputResults(settings.getInputOption(), false, null, null);
    if (MemoryPeakResults.isEmpty(results)) {
      IJ.error(TITLE, "No results could be loaded");
      instanceLock.release();
      return;
    }

    // Map all non-zero IDs to a natural series. This avoids issues with sparse cluster Ids.

    // Show a super-resolution image where clusters can be selected.
    final Rectangle bounds = results.getBounds();
    final PeakResultsList resultsList = new PeakResultsList();
    ResultsManager.addImageResults(resultsList, settings.getResultsImageSettings(), bounds, 0);
    resultsList.copySettings(results);
    resultsList.begin();
    resultsList.addAll(results.toArray());
    resultsList.end();
    image = (ImageJImagePeakResults) resultsList.getOutput(0);

    // Note: Setting the lut name in the image only has an effect if the image is not showing
    // thus the lut is applied afterwards.
    final ImagePlus imp = image.getImagePlus();
    if (TextUtils.isNotEmpty(image.getLutName())) {
      imp.setLut(LutHelper.createLut(LutColour.forName(image.getLutName()), true));
    }

    final ImageCanvas canvas = imp.getCanvas();
    if (canvas == null) {
      instanceLock.release();
      return;
    }

    // Set-up analysis processing:
    // Store latest image ROI bounds and analysis settings.
    // ConcurrentMonoStack to store next image ROI bounds and analysis settings.
    // When image is clicked submit for analysis.
    // When settings are changed submit for analysis.
    // Submit for analysis checks if ROI is area ROI. if so it:
    // - adds current settings to the next analysis monostack
    // - acquires a softlock and if available submits a runnable to do the analysis until
    // the monostack is empty
    lock = new SoftLock();
    workQueue = new ConcurrentMonoStack<>();
    previous = Pair.of(null, settings.build());
    executor = Executors.newSingleThreadExecutor();

    // Create the bounds and activation times for each cluster
    clusterData = createClusterData(results);

    // Add interactive monitor to the image where clusters can be selected.
    // For all selected clusters show on an Activations-vs-Time plot.
    final MouseListener imageListener = new MouseAdapter() {
      @Override
      public void mouseReleased(MouseEvent e) {
        addWork(imp.getRoi());
      }
    };
    canvas.addMouseListener(imageListener);

    // Add monitor for the selection of clusters in the current clusters table
    clusterSelectedListener = new ClusterSelectedListener();

    // Allow analysis of the activations-vs-time data for start/end of bursts based on a
    // steepness parameter: local window size (sec) and activation rate (per sec)
    try {
      showAnalysisDialog();
    } finally {
      canvas.removeMouseListener(imageListener);
      final CurrentClustersTextWindow tw = resultsWindowRef.get();
      if (tw != null) {
        tw.setAction(null);
      }
      executor.shutdown();
      instanceLock.release();
    }
  }

  /**
   * Show dialog.
   *
   * @return true, if successful
   */
  private boolean showDialog() {
    settings = SettingsManager.readTcPalmAnalysisSettings(0).toBuilder();
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addMessage("Analyse the time-correlated activation of traced data");
    ResultsManager.addInput(gd, "Input", settings.getInputOption(), InputSource.MEMORY_CLUSTERED);
    // Require results settings to use the standard ResultsManager image options
    final ResultsSettings.Builder tmp = ResultsSettings.newBuilder();
    tmp.setResultsImageSettings(settings.getResultsImageSettingsBuilder());
    final int flags = ResultsManager.FLAG_NO_SECTION_HEADER | ResultsManager.FLAG_IMAGE_REMOVE_NONE;
    ResultsManager.addImageResultsOptions(gd, tmp, flags);
    gd.addHelp(HelpUrls.getUrl("tc-palm-analysis"));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }
    settings.setInputOption(ResultsManager.getInputSource(gd));
    final ResultsImageSettings.Builder newImgSettings = tmp.getResultsImageSettingsBuilder();
    // Note: The initial none option was removed
    newImgSettings.setImageTypeValue(gd.getNextChoiceIndex() + 1);
    settings.setResultsImageSettings(newImgSettings);
    SettingsManager.writeSettings(settings);
    return true;
  }

  @SuppressWarnings("null")
  private static LocalList<ClusterData> createClusterData(MemoryPeakResults results) {
    results = results.copy();
    results.sort(IdFramePeakResultComparator.INSTANCE);
    final LocalList<ClusterData> clusterData = new LocalList<>();
    final FrameCounter counter = new FrameCounter(results.getFirst().getId() - 1);
    ClusterData data = null;
    final int size = results.size();
    for (int i = 0; i < size; i++) {
      final PeakResult result = results.get(i);
      if (counter.advance(result.getId())) {
        clusterData.add(data);
        data = new ClusterData(result);
      } else {
        data.add(result);
      }
    }
    // Final cluster
    clusterData.add(data);
    // Remove the first null object and compact the frame arrays
    clusterData.remove(0);
    clusterData.forEach(ClusterData::finalise);
    // Sort by time then cluster ID
    clusterData.sort((c1, c2) -> {
      final int result = Integer.compare(c1.start, c2.start);
      if (result != 0) {
        return result;
      }
      // result = Integer.compare(c1.end, c2.end);
      // if (result != 0) {
      // return result;
      // }
      return Integer.compare(c1.id, c2.id);
    });
    final int[] index = {0};
    clusterData.forEach(c -> c.index = index[0]++);
    return clusterData;
  }

  /**
   * Show analysis dialog.
   *
   * @return true, if successful
   */
  private boolean showAnalysisDialog() {
    final NonBlockingExtendedGenericDialog gd = new NonBlockingExtendedGenericDialog(TITLE);
    gd.addMessage("Analyse the time-correlated activation of traced data");
    gd.addCheckbox("ROI_intersects", settings.getIntersects());
    gd.addCheckbox("Time_in_seconds", settings.getTimeInSeconds());
    minT = results.getMinFrame();
    maxT = results.getMaxFrame();
    // No need to check other end of the range as the dialog slider will clip to the range.
    if (settings.getMinFrame() > maxT) {
      settings.setMinFrame(minT);
    }
    if (settings.getMaxFrame() < minT) {
      settings.setMaxFrame(maxT);
    }
    gd.addSlider("Min_frame", minT, maxT, settings.getMinFrame());
    gd.addSlider("Max_frame", minT, maxT, settings.getMaxFrame());
    gd.addCheckbox("Fixed_time_axis", settings.getFixedTimeAxis());
    gd.addDialogListener(this::readDialog);
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }
    return true;
  }

  private boolean readDialog(GenericDialog gd, @SuppressWarnings("unused") AWTEvent e) {
    settings.setIntersects(gd.getNextBoolean());
    settings.setTimeInSeconds(gd.getNextBoolean());
    settings.setMinFrame((int) gd.getNextNumber());
    settings.setMaxFrame((int) gd.getNextNumber());
    settings.setFixedTimeAxis(gd.getNextBoolean());
    addWork(previous.getLeft());
    return true;
  }

  /**
   * Add work to the queue and submit a job to process the queue if one is not already running.
   *
   * @param roi the roi
   */
  private void addWork(Roi roi) {
    if (roi == null || !roi.isArea()) {
      return;
    }
    addWork(roi.getBounds());
  }

  /**
   * Add work to the queue and submit a job to process the queue if one is not already running.
   *
   * @param bounds the bounds
   */
  private void addWork(Rectangle bounds) {
    if (bounds == null || executor.isShutdown()) {
      return;
    }
    workQueue.insert(Pair.of(new Rectangle(bounds), settings.build()));
    if (lock.acquire()) {
      executor.submit(() -> {
        Pair<Rectangle, TcPalmAnalysisSettings> current = previous;
        try {
          for (;;) {
            final Pair<Rectangle, TcPalmAnalysisSettings> next = workQueue.poll();
            if (next == null) {
              // queue is empty
              break;
            }
            if (!current.equals(next)) {
              // Settings have changed
              runAnalysis(next.getLeft(), next.getRight());
              current = next;
            }
          }
        } finally {
          previous = current;
          lock.release();
        }
      });
    }
  }

  /**
   * Run the analysis of clusters inside the bounds with the provided analysis settings.
   *
   * @param bounds the bounds
   * @param settings the settings
   */
  private void runAnalysis(Rectangle bounds, TcPalmAnalysisSettings settings) {
    // Map the bounds to the data bounds
    final double x = image.inverseMapX(bounds.x);
    final double y = image.inverseMapY(bounds.y);
    final double w = image.inverseMapX(bounds.x + bounds.width) - x;
    final double h = image.inverseMapY(bounds.y + bounds.height) - y;
    final Rectangle2D scaledBounds = new Rectangle2D.Double(x, y, w, h);

    // Identify all clusters using a custom filter
    BiPredicate<ClusterData, Rectangle2D> test = null;
    if (settings.getMinFrame() > minT) {
      final int min = settings.getMinFrame();
      test = (c, r) -> c.start >= min;
    }
    if (settings.getMaxFrame() < maxT) {
      final int max = settings.getMaxFrame();
      test = and(test, (c, r) -> c.end <= max);
    }
    final BiPredicate<ClusterData, Rectangle2D> filter =
        and(test, settings.getIntersects() ? ClusterData::intersects : ClusterData::isWithin);

    final LocalList<ClusterData> clusters = new LocalList<>();
    clusterData.forEach(c -> {
      if (filter.test(c, scaledBounds)) {
        clusters.add(c);
      }
    });
    if (clusters.isEmpty()) {
      return;
    }

    final WindowOrganiser wo = new WindowOrganiser();

    String timeLabel = "Time";
    UnaryOperator<float[]> timeConverter;
    if (settings.getTimeInSeconds()) {
      timeLabel += " (s)";
      final double scale = 1.0 / results.getCalibration().getTimeCalibration().getExposureTime();
      timeConverter = frames -> {
        final float[] updated = frames.clone();
        SimpleArrayUtils.apply(updated, f -> f *= scale);
        return updated;
      };
    } else {
      timeLabel += " (frame)";
      timeConverter = UnaryOperator.identity();
    }

    String title = TITLE + " Cluster Activations vs Time";
    final Plot plot = new Plot(title, timeLabel, "Cumulative count");
    plot.addLabel(0, 0, TextUtils.pleural(clusters.size(), "cluster"));
    final CumulativeCount count = new CumulativeCount();
    final TIntIntHashMap all = new TIntIntHashMap();
    clusters.forEach(c -> {
      final float[] frames = c.getFrames();
      for (final float t : frames) {
        all.adjustOrPutValue((int) t, 1, 1);
      }
      plot.addPoints(timeConverter.apply(frames), count.getCount(c.results.size()), Plot.LINE);
    });
    plot.draw();
    plot.setLimitsToFit(true);
    if (settings.getFixedTimeAxis()) {
      final double[] limits = plot.getLimits();
      limits[0] = minT - 1;
      limits[1] = maxT + 1;
      plot.setLimits(limits);
      plot.updateImage();
    }
    ImageJUtils.display(title, plot, wo);

    title = TITLE + " Total Activations vs Time";
    final Plot plot2 = new Plot(title, timeLabel, "Cumulative count");
    final int[] frames = all.keys();
    final int[] counts = all.values();
    SortUtils.sortData(counts, frames, true, false);
    final int localisations = (int) MathUtils.sum(counts);
    final int clashes = localisations - all.size();
    // Make counts cumulative
    for (int i = 1; i < counts.length; i++) {
      counts[i] += counts[i - 1];
    }
    plot2.addLabel(0, 0, TextUtils.pleural(localisations, "localisation") + " : " + clashes
        + TextUtils.pleuralise(clashes, " clash", " clashes"));
    plot2.addPoints(timeConverter.apply(SimpleArrayUtils.toFloat(frames)),
        SimpleArrayUtils.toFloat(counts), Plot.LINE);
    if (settings.getFixedTimeAxis()) {
      plot2.setLimits(minT - 1, maxT + 1, Double.NaN, Double.NaN);
    }
    ImageJUtils.display(title, plot2, wo);

    activationsPlotData = new ActivationsPlotData(plot2, timeConverter, frames, counts);

    wo.tile();

    // Add a table of the clusters (these should already be sorted by time then id).
    final CurrentClustersTextWindow resultsWindow = createTable();
    final TextPanel tp = resultsWindow.getTextPanel();
    tp.clear();
    try (BufferedTextWindow tw = new BufferedTextWindow(resultsWindow)) {
      final StringBuilder sb = new StringBuilder();
      clusters.forEach(c -> {
        sb.setLength(0);
        sb.append(c.index);
        sb.append('\t').append(c.id);
        final double[] xyz = c.getXyz();
        for (int i = 0; i < 3; i++) {
          sb.append('\t').append(MathUtils.round(xyz[i]));
        }
        sb.append('\t').append(c.results.size());
        sb.append('\t').append(c.start);
        sb.append('\t').append(c.end);
        tw.append(sb.toString());
      });
    }

    // Add a selected listener to allow selected clusters to be drawn on the
    // super-resolution image.
    resultsWindow.setAction(clusterSelectedListener);
  }

  /**
   * Create a combined {@code AND} predicate. If the first predicate is null the second is returned.
   *
   * @param <T> the generic type
   * @param <U> the generic type
   * @param first the first predicate (can be null)
   * @param second the second predicate
   * @return the AND predicate
   */
  private static <T, U> BiPredicate<T, U> and(BiPredicate<T, U> first, BiPredicate<T, U> second) {
    if (first == null) {
      return second;
    }
    return (T t, U u) -> first.test(t, u) && second.test(t, u);
  }

  /**
   * Creates the table for the current clusters.
   *
   * @return the text window
   */
  private static CurrentClustersTextWindow createTable() {
    return ImageJUtils.refresh(resultsWindowRef,
        () -> new CurrentClustersTextWindow(TITLE + " Current Clusters",
            "#\tID\tX\tY\tZ\tSize\tStart\tEnd", 600, 400));
  }
}
