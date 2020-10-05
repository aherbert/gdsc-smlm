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

import gnu.trove.list.array.TIntArrayList;
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
import java.awt.AWTEvent;
import java.awt.Color;
import java.awt.Rectangle;
import java.awt.Window;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.geom.Rectangle2D;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.BiPredicate;
import java.util.function.Consumer;
import java.util.function.IntUnaryOperator;
import java.util.function.UnaryOperator;
import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.RowSorter;
import javax.swing.SwingConstants;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.DefaultTableCellRenderer;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import uk.ac.sussex.gdsc.core.data.utils.Rounder;
import uk.ac.sussex.gdsc.core.data.utils.RounderUtils;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.NonBlockingExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.ScreenDimensionHelper;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper.LutColour;
import uk.ac.sussex.gdsc.core.ij.roi.CoordinatePredicate;
import uk.ac.sussex.gdsc.core.ij.roi.CoordinatePredicateUtils;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.SoftLock;
import uk.ac.sussex.gdsc.core.utils.SortUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.concurrent.ConcurrencyUtils;
import uk.ac.sussex.gdsc.core.utils.concurrent.ConcurrentMonoStack;
import uk.ac.sussex.gdsc.core.utils.function.FloatUnaryOperator;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.TcPalmAnalysisSettings;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageSettings;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsSettings;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.results.ImageJImagePeakResults;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
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
  private static AtomicReference<ClusterDataTableModelFrame> currentClustersTable =
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
  private ConcurrentMonoStack<Pair<Roi, TcPalmAnalysisSettings>> workQueue;

  /** The previous work. */
  private Pair<Roi, TcPalmAnalysisSettings> previous;

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
    LocalList<PeakResult> results;
    float x;
    float y;
    float width;
    float height;
    float[] frames;
    int start;
    int end;
    double[] xyz;

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
      // Set dimensions non-zero so the rectangle contains/intersects methods do not ignore
      // the zero sized cluster.
      width = Math.max(width - x, Float.MIN_VALUE);
      height = Math.max(height - y, Float.MIN_VALUE);
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
     * True if all points are within the coordinate predicate. The result coordinates are first
     * mapped to the domain of the predicate.
     *
     * @param r the predicate
     * @param mapX the function to map X
     * @param mapY the function to map Y
     * @return true if within
     */
    boolean isWithin(CoordinatePredicate r, FloatUnaryOperator mapX, FloatUnaryOperator mapY) {
      final int size = results.size();
      for (int i = 0; i < size; i++) {
        final PeakResult result = results.unsafeGet(i);
        if (!r.test(mapX.applyAsFloat(result.getXPosition()),
            mapY.applyAsFloat(result.getYPosition()))) {
          return false;
        }
      }
      // Assume size is above zero, i.e. return size != 0
      return true;
    }

    /**
     * True if any points are within the rectangle.
     *
     * @param r the rectangle
     * @return true if intersects
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
     * True if any points are within the coordinate predicate. The result coordinates are first
     * mapped to the domain of the predicate.
     *
     * @param r the predicate
     * @param mapX the function to map X
     * @param mapY the function to map Y
     * @return true if intersects
     */
    boolean intersects(CoordinatePredicate r, FloatUnaryOperator mapX, FloatUnaryOperator mapY) {
      final int size = results.size();
      for (int i = 0; i < size; i++) {
        final PeakResult result = results.unsafeGet(i);
        if (r.test(mapX.applyAsFloat(result.getXPosition()),
            mapY.applyAsFloat(result.getYPosition()))) {
          return true;
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
      double[] xyz = this.xyz;
      if (xyz == null) {
        final double[] tmp = new double[3];
        results.forEach(r -> {
          tmp[0] += r.getXPosition();
          tmp[1] += r.getYPosition();
          tmp[2] += r.getZPosition();
        });
        final int size = results.size();
        for (int i = 0; i < 3; i++) {
          tmp[i] /= size;
        }
        xyz = this.xyz = tmp;
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
   * Class to display data from the selected clusters.
   */
  private class ClusterSelectedListener implements Consumer<List<ClusterData>> {
    @Override
    public void accept(List<ClusterData> clusters) {
      final ActivationsPlotData activationsPlotData = TcPalmAnalysis.this.activationsPlotData;
      final Plot plot2 = activationsPlotData.plot;
      plot2.restorePlotObjects();

      // In case no cluster data is provided
      if (clusters.isEmpty()) {
        image.getImagePlus().setOverlay(null);
        plot2.updateImage();
        return;
      }

      // TODO
      // Time linkage clustering to find bursts. Assume the distance is close enough.
      // Note: A burst may be more than a single molecule as it can represent recruitment
      // of multiple proteins to the same position.
      // - Sort by start time.
      // - Maintain start and end of current cluster
      // - Join to next if the gap is below a threshold.
      // - Else finish and start a new cluster.
      // - Add all ranges to the plot.

      // Overlay the clusters on the image.
      final Overlay overlay = new Overlay();

      // Find bursts of continuous time
      final LocalList<int[]> bursts = new LocalList<>();
      final int gap = settings.getBurstMaximumGap();
      clusters.sort((c1, c2) -> Integer.compare(c1.start, c2.start));

      int start = clusters.get(0).start;
      int end = clusters.get(0).end;
      for (int i = 1; i < clusters.size(); i++) {
        // For each cluster add to the overlay.
        final ClusterData data = clusters.get(i);
        appendRoi(overlay, data);
        if (data.start - end > gap) {
          // Save the burst
          bursts.add(new int[] {start, end});
          start = data.start;
          end = data.end;
        } else {
          // extend the burst
          end = Math.max(data.end, end);
        }
      }
      bursts.add(new int[] {start, end});
      image.getImagePlus().setOverlay(overlay);

      // For each cluster add to the plot.
      plot2.setColor(Color.red);
      final TIntArrayList tmpFrames = new TIntArrayList();
      final TIntArrayList tmpCounts = new TIntArrayList();
      bursts.forEach(range -> {
        // Find the start and end points on the plotted data.
        final int is = ArrayUtils.indexOf(activationsPlotData.frames, range[0]);
        final int ie = ArrayUtils.indexOf(activationsPlotData.frames, range[1]);

        // Pad the line with zeros at the end
        tmpFrames.resetQuick();
        tmpCounts.resetQuick();
        int[] frames = Arrays.copyOfRange(activationsPlotData.frames, is, ie + 1);
        tmpFrames.add(frames[0]);
        tmpFrames.add(frames);
        tmpFrames.add(frames[frames.length - 1]);
        tmpCounts.add(0);
        int[] counts = Arrays.copyOfRange(activationsPlotData.counts, is, ie + 1);
        tmpCounts.add(counts);
        tmpCounts.add(0);
        frames = tmpFrames.toArray();
        counts = tmpCounts.toArray();

        // Add to the plot.
        plot2.addPoints(activationsPlotData.timeConverter.apply(SimpleArrayUtils.toFloat(frames)),
            SimpleArrayUtils.toFloat(counts), Plot.LINE);
      });
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

  /**
   * Class to show the ClusterData in a JTable.
   */
  private static class ClusterDataTableModel extends AbstractTableModel {
    private static final long serialVersionUID = 1L;

    List<ClusterData> data = Collections.emptyList();

    /**
     * Sets the data.
     *
     * @param data the new data
     */
    void setData(LocalList<ClusterData> data) {
      // Only update if the data is different.
      // This expects the selected clusters to not change very often.
      if (!data.equals(this.data)) {
        this.data = data;
        fireTableDataChanged();
      }
    }

    @Override
    public int getRowCount() {
      return data.size();
    }

    @Override
    public int getColumnCount() {
      // ID X Y Z Size Start End
      return 7;
    }

    @Override
    public String getColumnName(int columnIndex) {
      switch (columnIndex) {
        // @formatter:off
        case 0: return "ID";
        case 1: return "X";
        case 2: return "Y";
        case 3: return "Z";
        case 4: return "Size";
        case 5: return "Start";
        case 6: return "End";
        // @formatter:on
        default:
          throw new IndexOutOfBoundsException("Bad column: " + columnIndex);
      }
    }

    @Override
    public Class<?> getColumnClass(int columnIndex) {
      switch (columnIndex) {
        case 1:
        case 2:
        case 3:
          return Double.class;
        default:
          return Integer.class;
      }
    }

    @Override
    public Object getValueAt(int rowIndex, int columnIndex) {
      final ClusterData c = data.get(rowIndex);
      switch (columnIndex) {
        // @formatter:off
        case 0: return c.id;
        case 1: return c.getXyz()[0];
        case 2: return c.getXyz()[1];
        case 3: return c.getXyz()[2];
        case 4: return c.results.size();
        case 5: return c.start;
        case 6: return c.end;
        // @formatter:on
        default:
          throw new IndexOutOfBoundsException("Bad column: " + columnIndex);
      }
    }
  }

  /**
   * Class to display ClusterDataTableModel in a JTable.
   */
  private static class ClusterDataJTable extends JTable {
    private static final long serialVersionUID = 1L;

    /**
     * Create a new instance.
     *
     * @param model the model
     */
    ClusterDataJTable(ClusterDataTableModel model) {
      super(model);
      final Rounder rounder = RounderUtils.create(4);
      final DefaultTableCellRenderer renderer = new DefaultTableCellRenderer() {
        private static final long serialVersionUID = 1L;

        @Override
        protected void setValue(Object value) {
          // Boxed primitives should never be null
          setText(rounder.toString(((Number) value).doubleValue()));
        }
      };
      renderer.setHorizontalAlignment(SwingConstants.TRAILING);
      setDefaultRenderer(Float.class, renderer);
      setDefaultRenderer(Double.class, renderer);
      setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
      setAutoCreateRowSorter(true);
    }

    /**
     * Returns the data of all selected rows. This maps the indices from the view to the data model.
     *
     * @return an array containing the data of all selected rows, or an empty array if no row is
     *         selected
     * @see #getSelectedRow
     */
    List<ClusterData> getSelectedData() {
      final ClusterDataTableModel model = (ClusterDataTableModel) dataModel;
      final List<ClusterData> data = model.data;
      final int iMin = selectionModel.getMinSelectionIndex();
      final int iMax = selectionModel.getMaxSelectionIndex();

      // Any negative
      if ((iMin | iMax) < 0) {
        return Collections.emptyList();
      }

      final LocalList<ClusterData> rvTmp = new LocalList<>(1 + (iMax - iMin));

      final RowSorter<?> sorter = getRowSorter();
      final IntUnaryOperator map = sorter == null ? i -> i : sorter::convertRowIndexToModel;
      for (int i = iMin; i <= iMax; i++) {
        if (selectionModel.isSelectedIndex(i)) {
          rvTmp.add(data.get(map.applyAsInt(i)));
        }
      }
      return rvTmp;
    }
  }

  /**
   * Class to display ClusterDataTableModel in a window frame.
   */
  private static class ClusterDataTableModelFrame extends JFrame {
    private static final long serialVersionUID = 1L;

    final ClusterDataJTable table;
    Consumer<List<ClusterData>> selectedAction;

    /**
     * Create a new instance.
     *
     * @param model the model
     */
    ClusterDataTableModelFrame(ClusterDataTableModel model) {
      table = new ClusterDataJTable(model);
      final JScrollPane scroll = new JScrollPane(table);

      final ScreenDimensionHelper helper = new ScreenDimensionHelper();
      helper.setMinHeight(300);
      helper.setup(scroll);

      add(scroll);
      pack();

      addWindowListener(new WindowAdapter() {
        @Override
        public void windowOpened(WindowEvent event) {
          WindowManager.addWindow(ClusterDataTableModelFrame.this);
        }

        @Override
        public void windowClosing(WindowEvent event) {
          WindowManager.removeWindow(ClusterDataTableModelFrame.this);
        }
      });

      table.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
        @Override
        public void valueChanged(ListSelectionEvent e) {
          if (e.getValueIsAdjusting() || selectedAction == null) {
            return;
          }
          // Notify the action of the current selection
          selectedAction.accept(table.getSelectedData());
        }
      });
    }

    /**
     * Gets the model.
     *
     * @return the model
     */
    ClusterDataTableModel getModel() {
      return (ClusterDataTableModel) table.getModel();
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

    // Optionally group singles (id=0) together. The default is to allocated them an Id
    // so noise localisations are included in counts from selected regions.
    results = settings.getGroupSingles() ? results.copy() : results.copyAndAssignZeroIds();

    // Show a super-resolution image where clusters can be selected.
    final Rectangle bounds = results.getBounds();
    final PeakResultsList resultsList = new PeakResultsList();
    resultsList.copySettings(results);
    ResultsManager.addImageResults(resultsList, settings.getResultsImageSettings(), bounds, 0);
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
      // Remove the action from the single instance of the current clusters table
      final ClusterDataTableModelFrame frame = currentClustersTable.get();
      if (frame != null) {
        frame.selectedAction = null;
      }
      executor.shutdown();
      instanceLock.release();
      SettingsManager.writeSettings(settings);
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
    gd.addCheckbox("Group_singles", settings.getGroupSingles());
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
    settings.setGroupSingles(gd.getNextBoolean());
    final ResultsImageSettings.Builder newImgSettings = tmp.getResultsImageSettingsBuilder();
    // Note: The initial none option was removed
    newImgSettings.setImageTypeValue(gd.getNextChoiceIndex() + 1);
    settings.setResultsImageSettings(newImgSettings);
    SettingsManager.writeSettings(settings);
    return true;
  }

  @SuppressWarnings("null")
  private static LocalList<ClusterData> createClusterData(MemoryPeakResults results) {
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
    gd.addSlider("Rate_window", 0, 100, settings.getRateWindow());
    gd.addSlider("Burst_maximum_gap", 0, 100, settings.getBurstMaximumGap());
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
    settings.setRateWindow((int) gd.getNextNumber());
    settings.setBurstMaximumGap((int) gd.getNextNumber());
    addWork(previous.getLeft());
    return true;
  }

  /**
   * Add work to the queue and submit a job to process the queue if one is not already running.
   *
   * @param roi the roi
   */
  private void addWork(Roi roi) {
    if (roi == null || !roi.isArea() || executor.isShutdown()) {
      return;
    }
    workQueue.insert(Pair.of((Roi) roi.clone(), settings.build()));
    if (lock.acquire()) {
      executor.submit(() -> {
        Pair<Roi, TcPalmAnalysisSettings> current = previous;
        try {
          for (;;) {
            final Pair<Roi, TcPalmAnalysisSettings> next = workQueue.poll();
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
   * @param roi the roi
   * @param settings the settings
   */
  private void runAnalysis(Roi roi, TcPalmAnalysisSettings settings) {
    // Support square ROI
    // - Map image ROI bounds to the data
    // - Check if the cluster is inside the rectangle bounds

    final Rectangle bounds = roi.getBounds();
    final double x = image.inverseMapX(bounds.x);
    final double y = image.inverseMapY(bounds.y);
    final double w = image.inverseMapX(bounds.x + bounds.width) - x;
    final double h = image.inverseMapY(bounds.y + bounds.height) - y;
    final Rectangle2D scaledBounds = new Rectangle2D.Double(x, y, w, h);

    // Identify all clusters using a custom filter. Filter on time first as this is simple.
    BiPredicate<ClusterData, Rectangle2D> test = null;
    if (settings.getMinFrame() > minT) {
      final int min = settings.getMinFrame();
      test = (c, r) -> c.start >= min;
    }
    if (settings.getMaxFrame() < maxT) {
      final int max = settings.getMaxFrame();
      test = and(test, (c, r) -> c.end <= max);
    }
    // Add square bounds check
    test = and(test, settings.getIntersects() ? ClusterData::intersects : ClusterData::isWithin);

    // - If non-square ROI:
    // -- Map cluster coordinates back to image bounds
    // -- Check if the cluster is inside the ROI using a CoordinatePredicate
    if (roi.getType() != Roi.RECTANGLE || roi.getCornerDiameter() != 0) {
      final CoordinatePredicate pred = CoordinatePredicateUtils.createContainsPredicate(roi);
      if (settings.getIntersects()) {
        test = and(test, (c, r) -> c.intersects(pred, image::mapX, image::mapY));
      } else {
        test = and(test, (c, r) -> c.isWithin(pred, image::mapX, image::mapY));
      }
    }

    final LocalList<ClusterData> clusters = new LocalList<>();
    final BiPredicate<ClusterData, Rectangle2D> filter = test;
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
    double timeScale;
    if (settings.getTimeInSeconds()) {
      timeLabel += " (s)";
      timeScale = 1.0 / results.getCalibration().getTimeCalibration().getExposureTime();
      timeConverter = frames -> {
        final float[] updated = frames.clone();
        SimpleArrayUtils.apply(updated, f -> f *= timeScale);
        return updated;
      };
    } else {
      timeLabel += " (frame)";
      timeConverter = UnaryOperator.identity();
      timeScale = 1;
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
      limits[0] = timeScale * (minT - 1);
      limits[1] = timeScale * (maxT + 1);
      plot.setLimits(limits);
      plot.updateImage();
    }
    ImageJUtils.display(title, plot, wo);

    title = TITLE + " Total Activations vs Time";
    final Plot plot2 = new Plot(title, timeLabel, "Cumulative count");
    int[] frames = all.keys();
    int[] counts = all.values();
    SortUtils.sortData(counts, frames, true, false);
    final int localisations = (int) MathUtils.sum(counts);
    final int clashes = localisations - all.size();
    // Make counts cumulative
    for (int i = 1; i < counts.length; i++) {
      counts[i] += counts[i - 1];
    }

    // Expand the total activations with extra frames to allow large spans to plot as horizontal
    final TIntArrayList frames2 = new TIntArrayList(frames.length);
    final TIntArrayList counts2 = new TIntArrayList(frames.length);
    int previousT = Integer.MIN_VALUE;
    int previousC = 0;
    for (int i = 0; i < frames.length; i++) {
      if (frames[i] > previousT + 1) {
        frames2.add(frames[i] - 1);
        counts2.add(previousC);
      }
      previousT = frames[i];
      previousC = counts[i];
      frames2.add(previousT);
      counts2.add(previousC);
    }
    frames = frames2.toArray();
    counts = counts2.toArray();

    plot2.addLabel(0, 0, TextUtils.pleural(localisations, "localisation") + " : " + clashes
        + TextUtils.pleuralise(clashes, " clash", " clashes"));
    plot2.addPoints(timeConverter.apply(SimpleArrayUtils.toFloat(frames)),
        SimpleArrayUtils.toFloat(counts), Plot.LINE);
    if (settings.getFixedTimeAxis()) {
      plot2.setLimits(timeScale * (minT - 1), timeScale * (maxT + 1), Double.NaN, Double.NaN);
    }
    ImageJUtils.display(title, plot2, wo);

    activationsPlotData = new ActivationsPlotData(plot2, timeConverter, frames, counts);

    // Plot the gradient (activation rate) using a rolling window
    final PolynomialSplineFunction fun = new LinearInterpolator()
        .interpolate(SimpleArrayUtils.toDouble(frames), SimpleArrayUtils.toDouble(counts));
    title = TITLE + " Activation rate vs Time";
    final Plot plot4 = new Plot(title, timeLabel, "Activation rate (count/frame)");
    // Clip to the range
    final int min = frames[0];
    final int max = frames[frames.length - 1];
    final float[] frames3 = SimpleArrayUtils.newArray(max - min + 1, min, 1f);
    // Configurable window. This works even for a window of zero. Here interpolation is not
    // necessary but we leave it in place for simplicity.
    final int window = settings.getRateWindow();
    final float[] values = new float[frames3.length];
    for (int i = 0; i < values.length; i++) {
      final int low = Math.max(min, (int) frames3[i] - window - 1);
      final int high = Math.min(max, (int) frames3[i] + window);
      values[i] = (float) ((fun.value(high) - fun.value(low)) / (high - low));
    }
    // Note the initial frame will be an added frame with a count of zero. If the window is zero
    // then this will be NaN. Replace this with no rate.
    if (window == 0) {
      values[0] = 0;
    }
    plot4.addPoints(timeConverter.apply(frames3), values, Plot.LINE);
    if (settings.getFixedTimeAxis()) {
      plot4.setLimits(timeScale * (minT - 1), timeScale * (maxT + 1), Double.NaN, Double.NaN);
    }
    ImageJUtils.display(title, plot4, wo);

    wo.tile();

    // Add a table of the clusters.
    final ClusterDataTableModelFrame clustersTable = createClustersTable();
    clustersTable.getModel().setData(clusters);

    // Allow a configurable action that accepts the array of ClusterData that is selected.
    clustersTable.selectedAction = clusterSelectedListener;

    // In case the clusters are the same then ensure the selected clusters are redisplayed.
    clusterSelectedListener.accept(clustersTable.table.getSelectedData());
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
  private static ClusterDataTableModelFrame createClustersTable() {
    return ConcurrencyUtils.refresh(currentClustersTable, JFrame::isShowing, () -> {
      final ClusterDataTableModelFrame frame =
          new ClusterDataTableModelFrame(new ClusterDataTableModel());
      frame.setTitle(TITLE + " Selected Clusters");
      frame.setVisible(true);
      return frame;
    });
  }
}
