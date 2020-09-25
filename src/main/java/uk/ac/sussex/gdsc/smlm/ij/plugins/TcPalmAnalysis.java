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
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.Plot;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.text.TextWindow;
import java.awt.AWTEvent;
import java.awt.Rectangle;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.geom.Rectangle2D;
import java.util.Arrays;
import java.util.Objects;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.BiPredicate;
import java.util.function.UnaryOperator;
import org.apache.commons.lang3.tuple.Pair;
import uk.ac.sussex.gdsc.core.ij.BufferedTextWindow;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.NonBlockingExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper.LutColour;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.SoftLock;
import uk.ac.sussex.gdsc.core.utils.SortUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.concurrent.ConcurrentMonoStack;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageMode;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageType;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.results.ImageJImagePeakResults;
import uk.ac.sussex.gdsc.smlm.ij.results.ImagePeakResultsFactory;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.count.FrameCounter;
import uk.ac.sussex.gdsc.smlm.results.sort.IdFramePeakResultComparator;

/**
 * Analyses the time-correlated activation of traced localisation data.
 */
public class TcPalmAnalysis implements PlugIn {

  /** The plugin title. */
  private static final String TITLE = "TC PALM Analysis";

  /** Text window showing the current clusters. */
  private static AtomicReference<TextWindow> resultsWindowRef = new AtomicReference<>();

  /** The plugin settings. */
  private Settings settings;

  /** The results. */
  private MemoryPeakResults results;

  /** The image. */
  private ImagePlus imp;

  /** The lock (must be held when processing the work queue). */
  private SoftLock lock;

  /** The work queue. */
  private ConcurrentMonoStack<Pair<Rectangle, Settings>> workQueue;

  /** The previous work. */
  private Pair<Rectangle, Settings> previous;

  /** The executor. */
  private ExecutorService executor;

  /** The cluster data. */
  private LocalList<ClusterData> clusterData;

  /** The minimum time frame. */
  private int minT;

  /** The maximum time frame. */
  private int maxT;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    String inputOption;
    int imageSize;
    int lut;
    boolean intersects;
    boolean timeInSeconds;
    int minFrame;
    int maxFrame;

    /**
     * Instantiates a new settings.
     */
    Settings() {
      // Set defaults
      inputOption = "";
      imageSize = 1024;
    }

    /**
     * Instantiates a new settings.
     *
     * @param source the source
     */
    Settings(Settings source) {
      inputOption = source.inputOption;
      imageSize = source.imageSize;
      lut = source.lut;
      intersects = source.intersects;
      timeInSeconds = source.timeInSeconds;
      minFrame = source.minFrame;
      maxFrame = source.maxFrame;
    }

    /**
     * Copy.
     *
     * @return the settings
     */
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

    @Override
    public boolean equals(Object obj) {
      if (obj instanceof Settings) {
        final Settings other = (Settings) obj;
        return inputOption.equals(other.inputOption) && imageSize == other.imageSize
            && lut == other.lut && intersects == other.intersects
            && timeInSeconds == other.timeInSeconds && minFrame == other.minFrame
            && maxFrame == other.maxFrame;
      }
      return super.equals(obj);
    }

    @Override
    public int hashCode() {
      return Objects.hash(inputOption, imageSize, lut, intersects, timeInSeconds, minFrame,
          maxFrame);
    }
  }

  /**
   * Store data on each cluster.
   */
  private static class ClusterData {
    final int id;
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
      double[] xyz = new double[3];
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
    results = ResultsManager.loadInputResults(settings.inputOption, false, null, null);
    if (MemoryPeakResults.isEmpty(results)) {
      IJ.error(TITLE, "No results could be loaded");
      return;
    }

    // Map all non-zero IDs to a natural series. This avoids issues with sparse cluster Ids.

    // Show a super-resolution image where clusters can be selected.
    final Rectangle bounds = results.getBounds();
    final double scale = settings.imageSize / Math.max(bounds.width, bounds.height);
    final ImageJImagePeakResults image =
        ImagePeakResultsFactory.createPeakResultsImage(ResultsImageType.DRAW_ID, false, false,
            settings.inputOption, bounds, 0, 0, scale, 0, ResultsImageMode.IMAGE_MAX);
    image.copySettings(results);
    image.begin();
    image.addAll(results.toArray());
    image.end();

    // Note: Setting the lut name in the image only has an effect if the image is not showing
    // thus the lut is applied afterwards.
    imp = image.getImagePlus();
    imp.setLut(LutHelper.createLut(LutColour.forNumber(settings.lut), true));

    final ImageCanvas canvas = imp.getCanvas();
    if (canvas == null) {
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
    previous = Pair.of(null, settings);
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

    // Allow analysis of the activations-vs-time data for start/end of bursts based on a
    // steepness parameter: local window size (sec) and activation rate (per sec)
    try {
      showAnalysisDialog();
    } finally {
      canvas.removeMouseListener(imageListener);
      executor.shutdown();
    }
  }

  /**
   * Show dialog.
   *
   * @return true, if successful
   */
  private boolean showDialog() {
    settings = Settings.load();
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addMessage("Analyse the time-correlated activation of traced data");
    ResultsManager.addInput(gd, "Input", settings.inputOption, InputSource.MEMORY_CLUSTERED);
    gd.addNumericField("Image_size", settings.imageSize, 0);
    gd.addChoice("LUT", LutHelper.getLutNames(), settings.lut);
    gd.addHelp(HelpUrls.getUrl("tc-palm-analysis"));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }
    settings.inputOption = ResultsManager.getInputSource(gd);
    settings.imageSize = MathUtils.clip(20, 1 << 16, (int) gd.getNextNumber());
    settings.lut = gd.getNextChoiceIndex();
    settings.save();
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
    gd.addCheckbox("ROI_intersects", settings.intersects);
    gd.addCheckbox("Time_in_seconds", settings.timeInSeconds);
    minT = results.getMinFrame();
    maxT = results.getMaxFrame();
    // No need to check other end of the range as the dialog slider will clip to the range.
    if (settings.minFrame > maxT) {
      settings.minFrame = minT;
    }
    if (settings.maxFrame < minT) {
      settings.maxFrame = maxT;
    }
    gd.addSlider("Min_frame", minT, maxT, settings.minFrame);
    gd.addSlider("Max_frame", minT, maxT, settings.maxFrame);
    gd.addDialogListener(this::readDialog);
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }
    return true;
  }

  private boolean readDialog(GenericDialog gd, @SuppressWarnings("unused") AWTEvent e) {
    settings.intersects = gd.getNextBoolean();
    settings.timeInSeconds = gd.getNextBoolean();
    settings.minFrame = (int) gd.getNextNumber();
    settings.maxFrame = (int) gd.getNextNumber();
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
    if (bounds == null) {
      return;
    }
    workQueue.insert(Pair.of(new Rectangle(bounds), settings.copy()));
    if (lock.acquire()) {
      executor.submit(() -> {
        Pair<Rectangle, Settings> current = previous;
        try {
          for (;;) {
            final Pair<Rectangle, Settings> next = workQueue.poll();
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
  private void runAnalysis(Rectangle bounds, Settings settings) {
    // Map the bounds to the data bounds
    final Rectangle dataBounds = results.getBounds();
    final double scalex = (double) dataBounds.width / imp.getWidth();
    final double scaley = (double) dataBounds.height / imp.getHeight();
    final double x = dataBounds.x + scalex * bounds.x;
    final double y = dataBounds.y + scaley * bounds.y;
    final double w = scalex * bounds.width;
    final double h = scaley * bounds.height;
    final Rectangle2D scaledBounds = new Rectangle2D.Double(x, y, w, h);

    // Identify all clusters using a custom filter
    BiPredicate<ClusterData, Rectangle2D> test = null;
    if (settings.minFrame > minT) {
      final int min = settings.minFrame;
      test = (c, r) -> c.start >= min;
    }
    if (settings.maxFrame < maxT) {
      final int max = settings.maxFrame;
      test = and(test, (c, r) -> c.end <= max);
    }
    final BiPredicate<ClusterData, Rectangle2D> filter =
        and(test, settings.intersects ? ClusterData::intersects : ClusterData::isWithin);

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
    if (settings.timeInSeconds) {
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
    ImageJUtils.display(title, plot2, wo);

    wo.tile();

    // Add a table of the clusters (these should already be sorted by time then id).
    final TextWindow resultsWindow = createTable();
    resultsWindow.getTextPanel().clear();
    try (BufferedTextWindow tw = new BufferedTextWindow(resultsWindow)) {
      StringBuilder sb = new StringBuilder();
      clusters.forEach(c -> {
        sb.setLength(0);
        sb.append(c.id);
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

    // Add a selected listener to allow
    // selected clusters to be drawn on the super-resolution image.

    // Q. can the selection be individual lines selected with the ctrl key?
    // IIRC this cannot be done for an ImageJ table.
    // Look at using an interactive ImageJ table instead.

    // And draw an activation plot of the selected clusters with analysis of bursts...
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
  private static TextWindow createTable() {
    return ImageJUtils.refresh(resultsWindowRef, () -> new TextWindow(TITLE + " Analysis Results",
        "ID\tX\tY\tZ\tSize\tStart\tEnd", "", 600, 400));
  }
}
