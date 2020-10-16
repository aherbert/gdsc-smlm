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
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.gui.PointRoi;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.RoiListener;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import ij.process.LUT;
import java.awt.AWTEvent;
import java.awt.Color;
import java.awt.Component;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.TextField;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.geom.Rectangle2D;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.BiPredicate;
import java.util.function.Consumer;
import java.util.function.DoublePredicate;
import java.util.function.IntUnaryOperator;
import java.util.function.ToDoubleFunction;
import java.util.function.UnaryOperator;
import java.util.logging.Level;
import java.util.stream.DoubleStream;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.KeyStroke;
import javax.swing.ListSelectionModel;
import javax.swing.RowSorter;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.TableModelEvent;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableCellRenderer;
import uk.ac.sussex.gdsc.core.data.utils.IdentityTypeConverter;
import uk.ac.sussex.gdsc.core.data.utils.Rounder;
import uk.ac.sussex.gdsc.core.data.utils.RounderUtils;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot.HistogramPlotBuilder;
import uk.ac.sussex.gdsc.core.ij.ImageJPluginLoggerHelper;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.NonBlockingExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.ScreenDimensionHelper;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper.LutColour;
import uk.ac.sussex.gdsc.core.ij.roi.CoordinatePredicate;
import uk.ac.sussex.gdsc.core.ij.roi.CoordinatePredicateUtils;
import uk.ac.sussex.gdsc.core.math.hull.ConvexHull2d;
import uk.ac.sussex.gdsc.core.math.hull.Hull;
import uk.ac.sussex.gdsc.core.math.hull.Hull2d;
import uk.ac.sussex.gdsc.core.utils.FileUtils;
import uk.ac.sussex.gdsc.core.utils.LocalCollectors;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.SoftLock;
import uk.ac.sussex.gdsc.core.utils.SortUtils;
import uk.ac.sussex.gdsc.core.utils.StoredData;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.concurrent.ConcurrencyUtils;
import uk.ac.sussex.gdsc.core.utils.concurrent.ConcurrentMonoStack;
import uk.ac.sussex.gdsc.core.utils.function.FloatUnaryOperator;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationHelper;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.TcPalmAnalysisSettings;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageSettings;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsSettings;
import uk.ac.sussex.gdsc.smlm.data.config.UnitHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.TimeUnit;
import uk.ac.sussex.gdsc.smlm.ij.gui.TableColumnAdjuster;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.results.ImageJImagePeakResults;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.PeakResultsList;
import uk.ac.sussex.gdsc.smlm.results.Trace;
import uk.ac.sussex.gdsc.smlm.results.count.Counter;
import uk.ac.sussex.gdsc.smlm.results.count.FrameCounter;
import uk.ac.sussex.gdsc.smlm.results.sort.IdFramePeakResultComparator;

/**
 * Analyses the time-correlated activation of traced localisation data.
 */
public class TcPalmAnalysis implements PlugIn {

  /** The plugin title. */
  private static final String TITLE = "TC PALM Analysis";

  /** Text window showing the current localisation groups. */
  private static AtomicReference<ClusterDataTableModelFrame> currentGroupsTable =
      new AtomicReference<>();

  /** Text window showing the current clusters. */
  private static AtomicReference<ClusterDataTableModelFrame> currentClustersTable =
      new AtomicReference<>();

  /** Text window showing the clusters from all ROIs in the ROI Manager. */
  private static AtomicReference<ClusterDataTableModelFrame> allClustersTable =
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

  /** The data calibration. */
  private DataCalibration dataCalibration;

  /** The image. */
  private ImageJImagePeakResults image;

  /** The loop image of the selected region. */
  private LoopImage loopImage;

  /** The lock (must be held when processing the work queue). */
  private SoftLock lock;

  /** The work queue. */
  private ConcurrentMonoStack<Work> workQueue;

  /** The previous work. */
  private Work previous;

  /** The executor. */
  private ExecutorService executor;

  /** The cluster data. */
  private LocalList<ClusterData> clusterData;

  /** The currently selected clusters. */
  private LocalList<ClusterData> clusters;

  /** The count data of the selected clusters. */
  private CumulativeCountData countData;

  /** The current bursts of continuous time. */
  private LocalList<int[]> bursts;

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

  /** The colour map used to colour selected bursts. */
  private final ColourMap colourMap = new ColourMap();

  /** The dark time tolerance text field. */
  private TextField darkTimeToleranceTextField;

  /** The min cluster size text field. */
  private TextField minClusterSizeTextField;

  /**
   * Class to store the work.
   */
  private static class Work {
    final long timeout;
    final Roi roi;
    final TcPalmAnalysisSettings settings;

    /**
     * Create an instance.
     *
     * @param time the time
     * @param roi the roi
     * @param settings the settings
     */
    Work(long time, Roi roi, TcPalmAnalysisSettings settings) {
      this.timeout = time;
      this.roi = roi;
      this.settings = settings;
    }

    @Override
    public boolean equals(Object obj) {
      if (obj instanceof Work) {
        final Work other = (Work) obj;
        return Objects.equals(roi, other.roi) && Objects.equals(settings, other.settings);
      }
      return super.equals(obj);
    }

    @Override
    public int hashCode() {
      return Objects.hashCode(roi) ^ Objects.hashCode(settings);
    }
  }

  /**
   * Store data on each cluster.
   */
  private static class ClusterData {
    int id;
    LocalList<PeakResult> results;
    float x;
    float y;
    float width;
    float height;
    int start;
    int end;
    double[] xyz;
    Hull2d hull;
    double area = -1;
    Roi sourceRoi;

    ClusterData(PeakResult result) {
      id = result.getId();
      x = width = result.getXPosition();
      y = height = result.getYPosition();
      results = new LocalList<>();
      results.add(result);
    }

    ClusterData(int id, List<PeakResult> list) {
      this.id = id;
      results = new LocalList<>(list);
      results.sort((p1, p2) -> Integer.compare(p1.getFrame(), p2.getFrame()));
      final PeakResult result = results.unsafeGet(0);
      x = width = result.getXPosition();
      y = height = result.getYPosition();
      final int size = results.size();
      for (int i = 1; i < size; i++) {
        updateBounds(results.unsafeGet(i));
      }
      finish();
    }

    void add(PeakResult result) {
      updateBounds(result);
      results.add(result);
    }

    void updateBounds(PeakResult result) {
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
    }

    /**
     * Called when no more results will be added to the cluster data.
     */
    void finish() {
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

    Hull2d getHull() {
      Hull2d hull = this.hull;
      if (hull == null) {
        final Hull.Builder builder = ConvexHull2d.newBuilder();
        results.forEach(r -> {
          builder.add(r.getXPosition(), r.getYPosition());
        });
        hull = this.hull = (Hull2d) builder.build();
      }
      return hull;
    }

    double getArea() {
      double area = this.area;
      if (area < 0) {
        final Hull2d hull = getHull();
        area = this.area = hull == null ? 0 : hull.getArea();
      }
      return area;
    }

    int getDuration() {
      return end - start + 1;
    }
  }

  /**
   * Store the cumulative counts.
   */
  private static class CumulativeCountData {
    final int[] frames;
    final int[] counts;

    /** The frames padded for plotting. */
    final int[] plotFrames;
    /** The cumulative counts padded for plotting. */
    final int[] plotCounts;

    CumulativeCountData(int[] frames, int[] counts, boolean createPlotData) {
      this.frames = frames;
      this.counts = counts;

      if (createPlotData) {
        // Expand the total activations with extra frames to allow large spans to plot as horizontal
        final TIntArrayList tmpFrames = new TIntArrayList(frames.length);
        final TIntArrayList tmpCounts = new TIntArrayList(frames.length);

        int previousT = Integer.MIN_VALUE;
        int previousC = 0;
        for (int i = 0; i < frames.length; i++) {
          if (frames[i] > previousT + 1) {
            tmpFrames.add(frames[i] - 1);
            tmpCounts.add(previousC);
          }
          previousT = frames[i];
          previousC += counts[i];
          tmpFrames.add(previousT);
          tmpCounts.add(previousC);
        }
        plotFrames = tmpFrames.toArray();
        plotCounts = tmpCounts.toArray();
      } else {
        plotFrames = null;
        plotCounts = null;
      }
    }
  }

  /**
   * Class to display data from the selected clusters.
   */
  private class ClusterSelectedListener implements Consumer<List<ClusterData>> {
    @Override
    public void accept(List<ClusterData> clusters) {
      final LocalList<LocalList<PeakResult>> burstLocalisations = new LocalList<>();
      final LocalList<int[]> bursts = new LocalList<>();

      // In case no cluster data is provided
      if (clusters.isEmpty()) {
        runBurstOverlay(burstLocalisations);
        runBurstPlotSelection(bursts);
        return;
      }

      // Time linkage clustering to find bursts. Assume the distance is close enough.
      // Note: A burst may be more than a single molecule as it can represent recruitment
      // of multiple proteins to the same position.
      // - Sort by start time.
      // - Maintain start and end of current cluster
      // - Join to next if the gap is below a threshold.
      // - Else finish and start a new cluster.
      // - Add all ranges to the plot.

      // Find bursts of continuous time
      final int gap = settings.getDarkTimeTolerance() + 1;
      clusters.sort((c1, c2) -> Integer.compare(c1.start, c2.start));

      final ClusterData first = clusters.get(0);
      final LocalList<PeakResult> results = new LocalList<>();
      results.addAll(first.results);
      int start = first.start;
      int end = first.end;
      for (int i = 1; i < clusters.size(); i++) {
        // For each cluster add to the overlay.
        final ClusterData data = clusters.get(i);
        if (data.start - end > gap) {
          // Save the burst
          bursts.add(new int[] {start, end});
          burstLocalisations.add(results.copy());
          start = data.start;
          end = data.end;
          results.clear();
        } else {
          // extend the burst
          end = Math.max(data.end, end);
        }
        results.addAll(data.results);
      }
      bursts.add(new int[] {start, end});
      burstLocalisations.add(results);

      runBurstOverlay(burstLocalisations);
      runBurstPlotSelection(bursts);
    }
  }

  /**
   * Class containing the data for the most recent total activations plot for current clusters.
   */
  private static class ActivationsPlotData {
    final Plot plot;
    final UnaryOperator<float[]> timeConverter;
    final CumulativeCountData data;

    ActivationsPlotData(Plot plot, UnaryOperator<float[]> timeConverter, CumulativeCountData data) {
      this.plot = plot;
      this.timeConverter = timeConverter;
      this.data = data;
      plot.savePlotObjects();
    }
  }

  /**
   * Class to calibration data from frames and pixels to a time unit and distance unit.
   */
  private static class DataCalibration {
    static final DataCalibration DEFAULT = new DataCalibration();
    final TypeConverter<TimeUnit> timeConverter;
    final TypeConverter<DistanceUnit> distanceConverter;

    DataCalibration() {
      timeConverter = new IdentityTypeConverter<>(TimeUnit.TIME_UNIT_NA);
      distanceConverter = new IdentityTypeConverter<>(DistanceUnit.DISTANCE_UNIT_NA);
    }

    DataCalibration(Calibration cal) {
      timeConverter = CalibrationHelper.getTimeConverterSafe(cal, TimeUnit.SECOND);
      distanceConverter = CalibrationHelper.getDistanceConverterSafe(cal, DistanceUnit.UM);
    }

    String getTimeUnitName() {
      return UnitHelper.getShortName(timeConverter.to());
    }

    String getDistanceUnitName() {
      return UnitHelper.getShortName(distanceConverter.to());
    }

    double convertArea(double area) {
      return distanceConverter.convert(distanceConverter.convert(area));
    }
  }

  /**
   * Class to show the ClusterData in a JTable.
   */
  private static class ClusterDataTableModel extends AbstractTableModel {
    private static final long serialVersionUID = 1L;

    /** Flag to indicate the cluster statistics (area, density) should be displayed. */
    final boolean withStatistics;

    List<ClusterData> data = Collections.emptyList();
    DataCalibration calibration = DataCalibration.DEFAULT;

    ClusterDataTableModel(boolean withStatistics) {
      this.withStatistics = withStatistics;
    }

    /**
     * Sets the data.
     *
     * @param data the new data
     * @param calibration the calibration
     */
    void setData(List<ClusterData> data, DataCalibration calibration) {
      // Only update if the data is different.
      // This expects the selected clusters to not change very often.
      if (this.calibration != calibration) {
        this.calibration = calibration;
        fireTableStructureChanged();
      }
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
      return withStatistics ? 11 : 9;
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
        case 7: return "Duration";
        case 8: return "Duration (" + calibration.getTimeUnitName() + ")";
        case 9: return "Area (" + calibration.getDistanceUnitName() + "^2)";
        case 10: return "Density (" + calibration.getDistanceUnitName() + "^-2)";
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
        case 8:
        case 9:
        case 10:
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
        case 7: return c.getDuration();
        case 8: return calibration.timeConverter.convert(c.getDuration());
        case 9: return calibration.convertArea(c.getArea());
        case 10: return c.results.size() / calibration.convertArea(c.getArea());
        // @formatter:on
        default:
          throw new IndexOutOfBoundsException("Bad column: " + columnIndex);
      }
    }
  }

  /**
   * A simple renderer to colour the ID column.
   */
  private static class ForegroundColouredTableCellRenderer extends DefaultTableCellRenderer {
    private static final long serialVersionUID = 1L;

    ForegroundColouredTableCellRenderer(Color colour) {
      setBackground(colour);
      setHorizontalAlignment(SwingConstants.TRAILING);
    }

    @Override
    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected,
        boolean hasFocus, int row, int column) {
      // Ignore selected status of the row. The cell will paint the same if selected or not.
      return super.getTableCellRendererComponent(table, value, false, hasFocus, row, column);
    }
  }

  /**
   * Class to display ClusterDataTableModel in a JTable.
   */
  private static class ClusterDataJTable extends JTable {
    private static final long serialVersionUID = 1L;

    private boolean sized;

    /** The colour map used for the ID column. */
    ColourMap colourMap;

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

    @Override
    public TableCellRenderer getCellRenderer(int row, int column) {
      // Override to selectively colour the first column (using the cluster ID)
      if (column == 0 && colourMap != null) {
        final IntUnaryOperator map = getRowIndexToModelFunction();
        final ClusterDataTableModel model = (ClusterDataTableModel) dataModel;
        final List<ClusterData> data = model.data;
        final ClusterData c = data.get(map.applyAsInt(row));
        final Color colour = colourMap.getColour(c.id - 1);
        return new ForegroundColouredTableCellRenderer(colour);
      }
      return super.getCellRenderer(row, column);
    }

    @Override
    public void tableChanged(final TableModelEvent event) {
      super.tableChanged(event);

      if (!sized && getModel().getRowCount() != 0) {
        sized = true;
        // The whole thing changed so resize the columns
        SwingUtilities.invokeLater(() -> {
          final TableColumnAdjuster tca = new TableColumnAdjuster(this, 6, false);
          // Only process 10 rows (5 at start, 5 at end).
          tca.setMaxRows(5);
          tca.setOnlyAdjustLarger(true);
          tca.adjustColumns();
        });
      }
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

      final IntUnaryOperator map = getRowIndexToModelFunction();
      for (int i = iMin; i <= iMax; i++) {
        if (selectionModel.isSelectedIndex(i)) {
          rvTmp.add(data.get(map.applyAsInt(i)));
        }
      }
      return rvTmp;
    }

    IntUnaryOperator getRowIndexToModelFunction() {
      final RowSorter<?> sorter = getRowSorter();
      final IntUnaryOperator map = sorter == null ? i -> i : sorter::convertRowIndexToModel;
      return map;
    }

    IntUnaryOperator getModelIndexToRowFunction() {
      final RowSorter<?> sorter = getRowSorter();
      final IntUnaryOperator map = sorter == null ? i -> i : sorter::convertRowIndexToView;
      return map;
    }
  }

  /**
   * Class to display ClusterDataTableModel in a window frame.
   */
  private static class ClusterDataTableModelFrame extends JFrame implements ActionListener {
    private static final long serialVersionUID = 1L;

    private JMenuItem fileSave;
    final ClusterDataJTable table;
    Consumer<List<ClusterData>> selectedAction;
    private String saveName;

    /**
     * Create a new instance.
     *
     * @param model the model
     */
    ClusterDataTableModelFrame(ClusterDataTableModel model) {
      if (model.withStatistics) {
        setJMenuBar(createMenuBar());
      }

      table = new ClusterDataJTable(model);
      final JScrollPane scroll = new JScrollPane(table);

      final ScreenDimensionHelper helper = new ScreenDimensionHelper();
      helper.setMinHeight(250);
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

    private JMenuBar createMenuBar() {
      final JMenuBar menubar = new JMenuBar();
      menubar.add(createFileMenu());
      return menubar;
    }

    private JMenu createFileMenu() {
      final JMenu menu = new JMenu("File");
      menu.setMnemonic(KeyEvent.VK_F);
      menu.add(fileSave = add("Save ...", KeyEvent.VK_S, "ctrl pressed S"));
      return menu;
    }

    private JMenuItem add(String text, int mnemonic, String keyStroke) {
      final JMenuItem item = new JMenuItem(text, mnemonic);
      if (keyStroke != null) {
        item.setAccelerator(KeyStroke.getKeyStroke(keyStroke));
      }
      item.addActionListener(this);
      return item;
    }

    @Override
    public void actionPerformed(ActionEvent event) {
      final Object src = event.getSource();
      if (src == fileSave) {
        doFileSave();
      }
    }

    private void doFileSave() {
      final ClusterDataTableModel model = getModel();
      if (model == null || model.getRowCount() == 0) {
        return;
      }
      final String filename = ImageJUtils.getFilename("Results_set_name", saveName);
      if (filename == null) {
        return;
      }
      saveName = FileUtils.replaceExtension(filename, ".csv");
      final String newLine = System.lineSeparator();
      try (BufferedWriter out = Files.newBufferedWriter(Paths.get(saveName))) {
        final StringBuilder sb = new StringBuilder();
        // Use the same column names
        for (int i = 0, cols = model.getColumnCount(); i < cols; i++) {
          sb.append(model.getColumnName(i)).append(',');
        }
        write(newLine, out, sb);
        // Return the same order as the table
        final IntUnaryOperator map = table.getRowIndexToModelFunction();
        for (int j = 0, rows = model.getRowCount(); j < rows; j++) {
          final int row = map.applyAsInt(j);
          for (int i = 0, cols = model.getColumnCount(); i < cols; i++) {
            sb.append(model.getValueAt(row, i)).append(',');
          }
          write(newLine, out, sb);
        }
      } catch (final IOException e) {
        // Let ImageJ display the exception
        throw new RuntimeException(e);
      }
    }

    /**
     * Clips the last character from the StringBuilder then writes to the output with a newline.
     * Used to output delimited text to file which has a trailing delimiter. Resets the length to
     * zero for fresh usage.
     *
     * @param newLine the new line
     * @param out the out
     * @param sb the StringBuilder
     * @throws IOException Signals that an I/O exception has occurred.
     */
    private static void write(final String newLine, BufferedWriter out, final StringBuilder sb)
        throws IOException {
      sb.setLength(sb.length() - 1);
      sb.append(newLine);
      out.write(sb.toString());
      sb.setLength(0);
    }

    /**
     * Gets the model.
     *
     * @return the model
     */
    ClusterDataTableModel getModel() {
      return (ClusterDataTableModel) table.getModel();
    }

    /**
     * Select the cluster in the table that matches the given cluster data. The match uses the start
     * and end points which distinguish clusters.
     *
     * @param clusterData the cluster data
     */
    void select(ClusterData clusterData) {
      final ClusterDataTableModel model = getModel();
      final int rows = model.data.size();
      for (int i = 0; i < rows; i++) {
        final ClusterData tmp = model.data.get(i);
        if (tmp.start == clusterData.start && tmp.end == clusterData.end) {
          final int index = table.getModelIndexToRowFunction().applyAsInt(i);
          SwingUtilities.invokeLater(() -> table.setRowSelectionInterval(index, index));
          return;
        }
      }
    }
  }

  /**
   * Map an index to a colour.
   */
  private static class ColourMap {
    /** This should be a 256 colour LUT with no black. */
    final LUT lut = LutHelper.createLut(LutColour.PIMP);

    /**
     * Gets the colour.
     *
     * @param index the index
     * @return the colour
     */
    Color getColour(int index) {
      return LutHelper.getColour(lut, index & 0xff);
    }
  }

  /**
   * Encapsulate the high resolution loop image.
   */
  private static class LoopImage implements Consumer<List<ClusterData>> {
    private static final int STATE_IMAGE = 0;
    private static final int STATE_OVERLAY = 1;
    private static final int STATE_ROI = 2;
    private static final int STATE_DONE = 3;

    private int state;
    private MemoryPeakResults results;
    private TcPalmAnalysisSettings settings;
    private ColourMap colourMap;
    private LocalList<LocalList<PeakResult>> bursts;
    private List<ClusterData> clusters;
    private ImageJImagePeakResults image;

    LoopImage setResults(MemoryPeakResults results) {
      this.results = results;
      setBursts(new LocalList<>());
      state = STATE_IMAGE;
      return this;
    }

    LoopImage setSettings(TcPalmAnalysisSettings settings) {
      this.settings = settings;
      state = STATE_IMAGE;
      return this;
    }

    LoopImage setColourMap(ColourMap colourMap) {
      this.colourMap = colourMap;
      state = STATE_IMAGE;
      return this;
    }

    LoopImage setBursts(LocalList<LocalList<PeakResult>> bursts) {
      this.bursts = bursts;
      setClusters(Collections.emptyList());
      state = Math.min(state, STATE_OVERLAY);
      return this;
    }

    LoopImage setClusters(List<ClusterData> clusters) {
      this.clusters = clusters;
      state = Math.min(state, STATE_ROI);
      return this;
    }

    /**
     * Update the image.
     *
     * @return the loop image
     */
    LoopImage update() {
      if (state <= STATE_IMAGE) {
        // Allow for no loop image
        final ResultsImageSettings imageSettings = settings.getLoopImageSettings();
        if (imageSettings.getImageTypeValue() <= 0) {
          // No loop image
          if (image != null) {
            image.getImagePlus().close();
            image = null;
          }
          return this;
        }
        drawImage(imageSettings);
      }
      if (state <= STATE_OVERLAY) {
        drawOverlay();
      }
      if (state <= STATE_ROI) {
        drawRoi();
      }
      state = STATE_DONE;
      return this;
    }

    private void drawImage(final ResultsImageSettings imageSettings) {
      final Rectangle bounds = results.getBounds(true);
      final PeakResultsList resultsList = new PeakResultsList();
      resultsList.copySettings(results);
      // In case the selection is empty
      bounds.width = Math.max(bounds.width, 1);
      bounds.height = Math.max(bounds.height, 1);
      resultsList.setBounds(bounds);
      resultsList.setName(TITLE + " Loop");
      // Compute the scale to generate a fixed size loop image
      final ResultsImageSettings.Builder builder = imageSettings.toBuilder();
      builder.setScale(settings.getLoopSize() / Math.max(bounds.width, bounds.height));
      ResultsManager.addImageResults(resultsList, builder.build(), bounds, 0);
      resultsList.begin();
      resultsList.addAll(results.toArray());
      resultsList.end();
      image = (ImageJImagePeakResults) resultsList.getOutput(0);
      final ImagePlus imp = image.getImagePlus();
      if (TextUtils.isNotEmpty(image.getLutName())) {
        imp.setLut(LutHelper.createLut(LutColour.forName(image.getLutName()), true));
      }
    }

    private void drawOverlay() {
      runBurstOverlay(bursts, image, colourMap);
    }

    private void drawRoi() {
      // We are expecting only a single cluster
      Hull2d hull;
      if (clusters.isEmpty() || (hull = clusters.get(0).getHull()) == null) {
        image.getImagePlus().deleteRoi();
        return;
      }
      final double[][] vertices = hull.getVertices();
      final int np = vertices.length;
      final float[] xp = new float[np];
      final float[] yp = new float[np];
      for (int j = 0; j < np; j++) {
        final double[] r = vertices[j];
        xp[j] = image.mapX((float) r[0]);
        yp[j] = image.mapY((float) r[1]);
      }
      final PolygonRoi roi = new PolygonRoi(xp, yp, Roi.POLYGON);
      image.getImagePlus().setRoi(roi);
    }

    @Override
    public void accept(List<ClusterData> clusters) {
      setClusters(clusters);
      update();
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

    // Allocate singles an id for analysis.
    results = results.copyAndAssignZeroIds();
    dataCalibration = new DataCalibration(results.getCalibration());

    // Show a super-resolution image where clusters can be selected.
    final Rectangle bounds = results.getBounds(true);
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
    // Use the current ROI (which may remain from previous plugin execution)
    previous = new Work(0, imp.getRoi(), settings.build());
    executor = Executors.newSingleThreadExecutor();

    // Create the bounds and activation times for each cluster
    clusterData = createClusterData(results);

    // Add interactive monitor to the image where clusters can be selected.
    // For all selected clusters show on an Activations-vs-Time plot.
    final RoiListener roiListener = new RoiListener() {
      @Override
      public void roiModified(ImagePlus imp2, int id) {
        if (imp2 != null && imp.getID() == imp2.getID()) {
          addWork(imp.getRoi());
        }
      }
    };
    Roi.addRoiListener(roiListener);

    // Add monitor for the selection of clusters in the current clusters table
    clusterSelectedListener = new ClusterSelectedListener();

    // Initialise the loop view
    loopImage = new LoopImage().setSettings(previous.settings).setColourMap(colourMap);

    // Allow analysis of the activations-vs-time data for start/end of bursts based on a
    // steepness parameter: local window size (sec) and activation rate (per sec)
    try {
      showAnalysisDialog();
    } finally {
      Roi.removeRoiListener(roiListener);
      // Remove the action from the single instance of the current clusters table
      removeListener(currentGroupsTable.get());
      removeListener(currentClustersTable.get());
      removeListener(allClustersTable.get());
      executor.shutdown();
      instanceLock.release();
      SettingsManager.writeSettings(settings);
    }
  }

  private static void removeListener(ClusterDataTableModelFrame frame) {
    if (frame != null) {
      frame.selectedAction = null;
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
    ResultsManager.addInput(gd, "Input", settings.getInputOption(), InputSource.MEMORY);
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
    clusterData.forEach(ClusterData::finish);
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
    gd.addCheckbox("Time_in_seconds", settings.getTimeInSeconds());
    darkTimeToleranceTextField =
        gd.addAndGetSlider("Dark_time_tolerance", 0, 100, settings.getDarkTimeTolerance());
    minClusterSizeTextField =
        gd.addAndGetSlider("Min_cluster_size", 0, 100, settings.getMinClusterSize());
    gd.addAndGetButton("Loop settings", this::showLoopSettingsDialog);
    gd.addAndGetButton("Analysis settings", this::showAnalysisSettingsDialog);
    gd.addAndGetButton("Analyse ROIs", this::analyseRois);
    gd.addDialogListener(this::readDialog);
    gd.hideCancelButton();
    gd.setOKLabel("Close");
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }
    return true;
  }

  /**
   * Show a dialog to change the loop settings.
   *
   * @param event the event
   */
  private void showLoopSettingsDialog(ActionEvent event) {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    // Require results settings to use the standard ResultsManager image options
    final ResultsSettings.Builder tmp = ResultsSettings.newBuilder();
    tmp.setResultsImageSettings(settings.getLoopImageSettingsBuilder());
    final int flags = ResultsManager.FLAG_NO_SECTION_HEADER;
    gd.addSlider("Loop_size", 128, 1024, settings.getLoopSize());
    ResultsManager.addImageResultsOptions(gd, tmp, flags);
    gd.addHelp(HelpUrls.getUrl("tc-palm-analysis"));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    // Restrict to a sensible range
    settings.setLoopSize(MathUtils.clip(32, 4096, (int) gd.getNextNumber()));
    final ResultsImageSettings.Builder newImgSettings = tmp.getResultsImageSettingsBuilder();
    // Note: The initial none option was removed
    newImgSettings.setImageTypeValue(gd.getNextChoiceIndex());
    settings.setLoopImageSettings(newImgSettings);
    SettingsManager.writeSettings(settings);
    addWork(previous.roi);
  }

  /**
   * Show a dialog to change the analysis settings.
   *
   * @param event the event
   */
  private void showAnalysisSettingsDialog(ActionEvent event) {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addCheckbox("Show_size_histogram", settings.getShowSizeHistogram());
    gd.addCheckbox("Show_duration_histogram", settings.getShowDurationHistogram());
    gd.addCheckbox("Show_area_histogram", settings.getShowAreaHistogram());
    gd.addCheckbox("Show_density_histogram", settings.getShowDensityHistogram());
    gd.addHelp(HelpUrls.getUrl("tc-palm-analysis"));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    settings.setShowSizeHistogram(gd.getNextBoolean());
    settings.setShowDurationHistogram(gd.getNextBoolean());
    settings.setShowAreaHistogram(gd.getNextBoolean());
    settings.setShowDensityHistogram(gd.getNextBoolean());
  }

  private boolean readDialog(GenericDialog gd, @SuppressWarnings("unused") AWTEvent e) {
    settings.setIntersects(gd.getNextBoolean());
    settings.setMinFrame((int) gd.getNextNumber());
    settings.setMaxFrame((int) gd.getNextNumber());
    settings.setFixedTimeAxis(gd.getNextBoolean());
    settings.setTimeInSeconds(gd.getNextBoolean());
    settings.setDarkTimeTolerance((int) gd.getNextNumber());
    settings.setMinClusterSize((int) gd.getNextNumber());
    addWork(previous.roi);
    return true;
  }

  /**
   * Add work to the queue and submit a job to process the queue if one is not already running. This
   * method has no action if the ROI is not an area or the flag to ignore analysis events is set.
   *
   * @param roi the roi
   */
  private void addWork(Roi roi) {
    if (roi == null || !roi.isArea()) {
      return;
    }
    addWork(System.currentTimeMillis() + 100, roi, settings.build(), null);
  }

  /**
   * Add work to the queue and submit a job to process the queue if one is not already running.
   *
   * <p>Optionally supply an action to perform after analysis. This will execute when no further
   * work is queued.
   *
   * @param timestamp the timeout timestamp
   * @param roi the roi
   * @param settings the settings
   * @param postAnalysisAction the post analysis action (can be null)
   */
  private void addWork(long timestamp, Roi roi, TcPalmAnalysisSettings settings,
      Runnable postAnalysisAction) {
    if (executor.isShutdown()) {
      return;
    }
    workQueue.insert(new Work(timestamp, (Roi) roi.clone(), settings));
    if (lock.acquire()) {
      executor.submit(() -> {
        Work current = previous;
        try {
          for (;;) {
            Work next = workQueue.poll();
            if (next == null) {
              // queue is empty
              break;
            }

            // Respect the delay
            if (next.timeout != 0) {
              long timeout = next.timeout;
              while (System.currentTimeMillis() < timeout) {
                Thread.sleep(50);
                // Assume new work can be added to the inbox.
                final Work newWork = workQueue.poll();
                if (newWork != null) {
                  timeout = newWork.timeout;
                  next = newWork;
                }
              }
            }

            if (!current.equals(next)) {
              // Settings have changed
              runAnalysis(current, next);
              current = next;
            }
          }
          if (postAnalysisAction != null) {
            // Set the latest work so any post analysis action uses the most recent settings
            previous = current;
            postAnalysisAction.run();
          }
        } catch (final InterruptedException e) {
          // Signal this was interrupted while waiting.
          // Note: This currently does not matter and is ignored.
          Thread.currentThread().interrupt();
        } catch (final RuntimeException ex) {
          ImageJPluginLoggerHelper.getLogger(TcPalmAnalysis.class).log(Level.WARNING,
              "Failed to execute work", ex);
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
   * @param last the last work
   * @param current the current work
   */
  private void runAnalysis(Work last, Work current) {
    final TcPalmAnalysisSettings lastSettings = last.settings;
    final TcPalmAnalysisSettings settings = current.settings;

    final boolean selectionChanged = selectionChanged(last, current);

    if (selectionChanged) {
      runSelection(current);
    } else if (loopSettingsChanged(lastSettings, settings)) {
      loopImage.setSettings(settings).update();
    }

    final boolean plotChanged = selectionChanged || plotSettingsChanged(lastSettings, settings);

    if (plotChanged) {
      runPlot(settings);
    }

    final boolean analysisChanged =
        selectionChanged || analysisSettingsChanged(lastSettings, settings);

    if (analysisChanged) {
      clearClustersTableSelection();
      bursts = runBurstAnalysis(settings, countData);
      runBurstOverlay(createBurstLocalisations(clusters, bursts));
    }

    if (plotChanged || analysisChanged) {
      runBurstPlotSelection(bursts);
    }
  }

  private static boolean selectionChanged(Work first, Work other) {
    return !Objects.equals(first.roi, other.roi)
        || selectionSettingsChanged(first.settings, other.settings);
  }

  private static boolean selectionSettingsChanged(TcPalmAnalysisSettings first,
      TcPalmAnalysisSettings second) {
    boolean result = (first.getIntersects() != second.getIntersects());
    result = result || (first.getMinFrame() != second.getMinFrame());
    result = result || (first.getMaxFrame() != second.getMaxFrame());
    return result;
  }

  private static boolean loopSettingsChanged(TcPalmAnalysisSettings first,
      TcPalmAnalysisSettings second) {
    boolean result = (first.getLoopSize() != second.getLoopSize());
    result = result || !first.getLoopImageSettings().equals(second.getLoopImageSettings());
    return result;
  }

  private static boolean plotSettingsChanged(TcPalmAnalysisSettings first,
      TcPalmAnalysisSettings second) {
    boolean result = (first.getTimeInSeconds() != second.getTimeInSeconds());
    result = result || (first.getFixedTimeAxis() != second.getFixedTimeAxis());
    return result;
  }

  private static boolean analysisSettingsChanged(TcPalmAnalysisSettings first,
      TcPalmAnalysisSettings second) {
    boolean result = (first.getDarkTimeTolerance() != second.getDarkTimeTolerance());
    result = result || (first.getMinClusterSize() != second.getMinClusterSize());
    return result;
  }

  /**
   * Run selection of the current clusters.
   *
   * @param work the current work
   */
  private void runSelection(Work work) {
    // Support square ROI
    // - Map image ROI bounds to the data
    // - Check if the cluster is inside the rectangle bounds

    final Roi roi = work.roi;
    final TcPalmAnalysisSettings settings = work.settings;
    final Rectangle2D scaledBounds = createScaledBounds(roi);
    final BiPredicate<ClusterData, Rectangle2D> filter = createSelectionFilter(roi, settings);

    clusters = new LocalList<>();
    clusterData.forEach(c -> {
      if (filter.test(c, scaledBounds)) {
        clusters.add(c);
      }
    });

    // Build total activations data
    countData = createCumulativeCountData(clusters, true);

    // Add a table of the clusters.
    final ClusterDataTableModelFrame clustersTable = createGroupsTable();
    clustersTable.getModel().setData(clusters, dataCalibration);

    // Allow a configurable action that accepts the array of ClusterData that is selected.
    clustersTable.selectedAction = clusterSelectedListener;

    // Extract the localisations into a loop view
    final MemoryPeakResults subset = new MemoryPeakResults();
    subset.copySettings(results);
    clusters.forEach(c -> {
      c.results.forEach(peak -> {
        if (scaledBounds.contains(peak.getXPosition(), peak.getYPosition())) {
          subset.add(peak);
        }
      });
    });
    // Clear bounds to force a recompute
    subset.setBounds(null);
    loopImage.setSettings(settings).setResults(subset).update();
  }

  /**
   * Creates the scaled bounds.
   *
   * @param roi the roi from the super-resolution image
   * @return the bounds scaled to the original data
   */
  private Rectangle2D createScaledBounds(final Roi roi) {
    final Rectangle bounds = roi.getBounds();
    final double x = image.inverseMapX(bounds.x);
    final double y = image.inverseMapY(bounds.y);
    final double w = image.inverseMapX(bounds.x + bounds.width) - x;
    final double h = image.inverseMapY(bounds.y + bounds.height) - y;
    return new Rectangle2D.Double(x, y, w, h);
  }

  /**
   * Creates the selection filter to identify all cluster groups using a custom filter.
   *
   * @param roi the roi
   * @param settings the settings
   * @return the bi predicate
   */
  private BiPredicate<ClusterData, Rectangle2D> createSelectionFilter(final Roi roi,
      final TcPalmAnalysisSettings settings) {
    // Filter on time first as this is simple.
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

    return test;
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
   * Creates the cumulative count data.
   *
   * @param createPlotData set to true to create the plot data arrays
   * @return the cumulative count data
   */
  private CumulativeCountData createCumulativeCountData(LocalList<ClusterData> clusters,
      boolean createPlotData) {
    final TIntIntHashMap all = new TIntIntHashMap(maxT - minT + 1);
    clusters.forEach(c -> c.results.forEach(peak -> all.adjustOrPutValue(peak.getFrame(), 1, 1)));
    final int[] frames = all.keys();
    final int[] counts = all.values();
    SortUtils.sortData(counts, frames, true, false);
    return new CumulativeCountData(frames, counts, createPlotData);
  }

  /**
   * Creates the table for the current groups.
   *
   * @return the text window
   */
  private static ClusterDataTableModelFrame createGroupsTable() {
    return ConcurrencyUtils.refresh(currentGroupsTable, JFrame::isShowing, () -> {
      final ClusterDataTableModelFrame frame =
          new ClusterDataTableModelFrame(new ClusterDataTableModel(false));
      frame.setTitle(TITLE + " Current Localisation Groups");
      frame.setVisible(true);
      return frame;
    });
  }

  /**
   * Creates the table for the current clusters and position it below the parent when created.
   *
   * @param parent the parent
   * @return the text window
   */
  private static ClusterDataTableModelFrame createClustersTable(ClusterDataTableModelFrame parent) {
    return ConcurrencyUtils.refresh(currentClustersTable, JFrame::isShowing, () -> {
      final ClusterDataTableModelFrame frame =
          new ClusterDataTableModelFrame(new ClusterDataTableModel(true));
      frame.setTitle(TITLE + " Current Clusters");
      frame.table.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
      if (parent != null) {
        final Point p = parent.getLocation();
        p.y += parent.getHeight();
        frame.setLocation(p);
      }
      frame.setVisible(true);
      return frame;
    });
  }

  /**
   * Creates the table for the clusters from batch analysis of all ROIs in the ROI manager.
   *
   * @return the text window
   */
  private static ClusterDataTableModelFrame createAllClustersTable() {
    return ConcurrencyUtils.refresh(allClustersTable, JFrame::isShowing, () -> {
      final ClusterDataTableModelFrame frame =
          new ClusterDataTableModelFrame(new ClusterDataTableModel(true));
      frame.setTitle(TITLE + " All Clusters");
      frame.table.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
      frame.setVisible(true);
      return frame;
    });
  }

  /**
   * Clear the manually created selection from the current clusters table.
   */
  private static void clearClustersTableSelection() {
    final ClusterDataTableModelFrame frame = currentGroupsTable.get();
    if (frame != null) {
      frame.table.clearSelection();
    }
  }

  /**
   * Create the activations and total activations plots.
   *
   * @param settings the settings
   */
  private void runPlot(TcPalmAnalysisSettings settings) {
    final WindowOrganiser wo = new WindowOrganiser();

    String timeLabel = "Time";
    UnaryOperator<float[]> timeConverter;
    double timeScale;
    if (settings.getTimeInSeconds()) {
      timeLabel += " (s)";
      timeScale = 1.0 / results.getCalibration().getTimeCalibration().getExposureTime();
      timeConverter = frames -> {
        SimpleArrayUtils.apply(frames, f -> f *= timeScale);
        return frames;
      };
    } else {
      timeLabel += " (frame)";
      timeConverter = UnaryOperator.identity();
      timeScale = 1;
    }

    final CumulativeCountData data = this.countData;

    String title = TITLE + " Activations vs Time";
    final Plot plot1 = new Plot(title, timeLabel, "Count");
    plot1.addLabel(0, 0, TextUtils.pleural(clusters.size(), "cluster"));
    for (int i = 0; i < data.frames.length; i++) {
      final int t = data.frames[i];
      final int c = data.counts[i];
      plot1.addPoints(timeConverter.apply(new float[] {t, t}), new float[] {0, c}, Plot.LINE);
    }
    plot1.draw();
    plot1.setLimitsToFit(true);
    if (settings.getFixedTimeAxis()) {
      final double[] limits = plot1.getLimits();
      limits[0] = timeScale * (minT - 1);
      limits[1] = timeScale * (maxT + 1);
      plot1.setLimits(limits);
      plot1.updateImage();
    }
    final PlotWindow pw1 = ImageJUtils.display(title, plot1, ImageJUtils.NO_TO_FRONT, wo);

    title = TITLE + " Total Activations vs Time";
    final Plot plot2 = new Plot(title, timeLabel, "Cumulative count");
    final int localisations =
        data.counts.length == 0 ? 0 : data.plotCounts[data.plotCounts.length - 1];
    final int clashes = localisations - data.counts.length;
    plot2.addLabel(0, 0, TextUtils.pleural(localisations, "localisation") + " : " + clashes
        + TextUtils.pleuralise(clashes, " clash", " clashes"));
    plot2.addPoints(timeConverter.apply(SimpleArrayUtils.toFloat(data.plotFrames)),
        SimpleArrayUtils.toFloat(data.plotCounts), Plot.LINE);
    if (settings.getFixedTimeAxis()) {
      plot2.setLimits(timeScale * (minT - 1), timeScale * (maxT + 1), Double.NaN, Double.NaN);
    }
    final PlotWindow pw2 = ImageJUtils.display(title, plot2, ImageJUtils.NO_TO_FRONT, wo);

    activationsPlotData = new ActivationsPlotData(plot2, timeConverter, data);

    // Simple tile one window above the other only if both are new.
    if (wo.size() == 2) {
      final Point p = pw1.getLocation();
      p.y += pw1.getHeight();
      pw2.setLocation(p);
    }
  }

  /**
   * Run the burst analysis on the activation counts. Identifies continuous bursts of activations
   * that are connected within the dark time tolerance and above the minimum cluster size.
   *
   * @param settings the settings
   */
  private static LocalList<int[]> runBurstAnalysis(TcPalmAnalysisSettings settings,
      CumulativeCountData data) {
    final LocalList<int[]> bursts = new LocalList<>();
    final int[] frames = data.frames;
    final int[] counts = data.counts;

    if (frames.length != 0) {
      final int gap = settings.getDarkTimeTolerance() + 1;
      final int size = settings.getMinClusterSize();
      int start = frames[0];
      int end = start;
      int count = counts[0];
      for (int i = 1; i < frames.length; i++) {
        final int t = frames[i];
        if (t - end > gap) {
          // Save the burst
          if (count >= size) {
            bursts.add(new int[] {start, end});
          }
          start = end = t;
          count = counts[i];
        } else {
          end = t;
          count += counts[i];
        }
      }
      if (count >= size) {
        bursts.add(new int[] {start, end});
      }
    }

    return bursts;
  }

  /**
   * Creates the burst localisations using the ranges of the current bursts.
   *
   * @param clusters the clusters
   * @param bursts the bursts
   * @return the burst localisations
   */
  private static LocalList<LocalList<PeakResult>>
      createBurstLocalisations(LocalList<ClusterData> clusters, LocalList<int[]> bursts) {
    final LocalList<LocalList<PeakResult>> burstsLocalisations = new LocalList<>();
    // TODO: Make more efficient. Order clusters by frame and burst by frame
    // Do a single sweep over the clusters and extract the peaks within the burst ranges.
    bursts.forEach(range -> {
      final int min = range[0];
      final int max = range[1];
      final LocalList<PeakResult> list = new LocalList<>();
      // Build from the current selection
      clusters.forEach(cluster -> {
        cluster.results.forEach(peak -> {
          final int t = peak.getFrame();
          if (t >= min && t <= max) {
            list.add(peak);
          }
        });
      });
      burstsLocalisations.add(list);
    });
    return burstsLocalisations;
  }

  /**
   * Select the activation bursts on the cumulative activations plot. Also stores the bursts in the
   * class property for future plot updates with the same selection.
   *
   * <p>This is synchronized so that the table selection events do not clash with the background
   * analysis jobs.
   *
   * @param bursts the bursts
   */
  private synchronized void runBurstPlotSelection(LocalList<int[]> bursts) {
    this.bursts = bursts;
    final Plot plot = activationsPlotData.plot;
    plot.restorePlotObjects();

    if (!bursts.isEmpty()) {
      final int[] plotFrames = activationsPlotData.data.plotFrames;
      final int[] plotCounts = activationsPlotData.data.plotCounts;

      // For each cluster add to the plot.
      final TIntArrayList tmpFrames = new TIntArrayList();
      final TIntArrayList tmpCounts = new TIntArrayList();
      for (int i = 0; i < bursts.size(); i++) {
        // Find the start and end points on the plotted data.
        final int[] range = bursts.unsafeGet(i);
        final int is = Arrays.binarySearch(plotFrames, range[0]);
        final int ie = Arrays.binarySearch(plotFrames, range[1]);

        // Pad the line with zeros at the end
        tmpFrames.resetQuick();
        tmpCounts.resetQuick();
        int[] frames = Arrays.copyOfRange(plotFrames, is, ie + 1);
        tmpFrames.add(frames[0]);
        tmpFrames.add(frames);
        tmpFrames.add(frames[frames.length - 1]);
        tmpCounts.add(0);
        int[] counts = Arrays.copyOfRange(plotCounts, is, ie + 1);
        tmpCounts.add(counts);
        tmpCounts.add(0);
        frames = tmpFrames.toArray();
        counts = tmpCounts.toArray();

        // Add to the plot.
        plot.setColor(colourMap.getColour(i));
        plot.addPoints(activationsPlotData.timeConverter.apply(SimpleArrayUtils.toFloat(frames)),
            SimpleArrayUtils.toFloat(counts), Plot.LINE);
      }
    }
    plot.updateImage();
  }

  /**
   * Select the data for the first cluster in the provided list on the cumulative activations plot.
   *
   * @param clusters the clusters
   */
  private void runClusterPlotSelection(List<ClusterData> clusters) {
    final Plot plot = activationsPlotData.plot;
    // We are expecting only a single cluster
    if (clusters.isEmpty()) {
      plot.getImagePlus().deleteRoi();
      return;
    }
    final int[] plotFrames = activationsPlotData.data.plotFrames;
    final int[] plotCounts = activationsPlotData.data.plotCounts;
    final ClusterData data = clusters.get(0);
    final int is = Arrays.binarySearch(plotFrames, data.start);
    final int ie = Arrays.binarySearch(plotFrames, data.end);
    final int low = plotCounts[is];
    final int high = plotCounts[ie];
    // Pad slightly. This handles the case where start == end to create a width and height
    final float[] x =
        activationsPlotData.timeConverter.apply(new float[] {data.start - 1, data.end + 1});
    final double x1 = plot.scaleXtoPxl(x[0]);
    final double x2 = plot.scaleXtoPxl(x[1]);
    final double y1 = plot.scaleYtoPxl(high + 1);
    final double y2 = plot.scaleYtoPxl(low - 1);
    plot.getImagePlus().setRoi(new Roi(x1, y1, x2 - x1, y2 - y1));
  }

  /**
   * Select the activation burst localisations on the super-resolution image. Create a table of the
   * clusters.
   *
   * @param bursts the bursts
   */
  private void runBurstOverlay(LocalList<LocalList<PeakResult>> bursts) {
    // Overlay the clusters on the image.
    runBurstOverlay(bursts, image, colourMap);
    loopImage.setBursts(bursts).update();

    // Add to an activation bursts (clusters) table
    final LocalList<ClusterData> clusters = new LocalList<>(bursts.size());
    bursts.forEach(list -> clusters.add(new ClusterData(clusters.size() + 1, list)));
    final ClusterDataTableModelFrame clustersTable = createClustersTable(currentGroupsTable.get());
    clustersTable.selectedAction = clusterData -> {
      loopImage.accept(clusterData);
      runClusterPlotSelection(clusterData);
    };
    clustersTable.table.colourMap = this.colourMap;
    clustersTable.getModel().setData(clusters, dataCalibration);
  }

  /**
   * Overlay the clusters on the image.
   *
   * @param bursts the bursts
   * @param image the image
   * @param colourMap the colour map
   */
  private static void runBurstOverlay(LocalList<LocalList<PeakResult>> bursts,
      ImageJImagePeakResults image, ColourMap colourMap) {
    if (bursts.isEmpty()) {
      image.getImagePlus().setOverlay(null);
      return;
    }
    final Overlay overlay = new Overlay();
    for (int i = 0; i < bursts.size(); i++) {
      final LocalList<PeakResult> results = bursts.unsafeGet(i);
      final int np = results.size();
      final float[] xp = new float[np];
      final float[] yp = new float[np];
      for (int j = 0; j < np; j++) {
        final PeakResult r = results.unsafeGet(j);
        xp[j] = image.mapX(r.getXPosition());
        yp[j] = image.mapY(r.getYPosition());
      }
      final PointRoi roi = new PointRoi(xp, yp, np);
      roi.setShowLabels(false);
      roi.setPointType(3);
      final Color c = colourMap.getColour(i);
      roi.setFillColor(c);
      roi.setStrokeColor(c);
      overlay.add(roi);
    }
    image.getImagePlus().setOverlay(overlay);
  }

  /**
   * Analyses all the ROIs in the ROI manager.
   *
   * @param event the event
   */
  private void analyseRois(ActionEvent event) {
    final RoiManager manager = RoiManager.getInstance();
    if (manager == null) {
      IJ.error(TITLE, "ROI manager is not open");
      return;
    }
    final LocalList<Roi> rois = Arrays.stream(manager.getRoisAsArray()).filter(Roi::isArea)
        .collect(LocalCollectors.toLocalList());
    if (rois.isEmpty()) {
      IJ.error(TITLE, "No area ROIs");
      return;
    }

    // Check for overlaps.
    if (anyOverlap(rois)) {
      final GenericDialog gd = new GenericDialog(TITLE);
      gd.addMessage(TextUtils.wrap("WARNING - Bounding rectangles of ROIs overlap. You can verify "
          + "the ROIs on the image using the ROI manager 'Show all' function.", 80));
      gd.setOKLabel("Continue");
      gd.showDialog();
      if (gd.wasCanceled()) {
        return;
      }
    }

    // For each ROI:
    // - Extract the current groups
    // - Build the cumulative count plot
    // - Identify the bursts
    // - Extract ClusterData for each burst
    final TcPalmAnalysisSettings settings = this.settings.build();
    final LocalList<ClusterData> allClusters = rois.parallelStream().map(roi -> {
      final Rectangle2D scaledBounds = createScaledBounds(roi);
      final BiPredicate<ClusterData, Rectangle2D> filter = createSelectionFilter(roi, settings);

      // Filter all cluster groups
      final LocalList<ClusterData> clusters = new LocalList<>();
      clusterData.forEach(c -> {
        if (filter.test(c, scaledBounds)) {
          clusters.add(c);
        }
      });

      // Extract activation bursts
      final CumulativeCountData countData = createCumulativeCountData(clusters, false);
      final LocalList<int[]> bursts = runBurstAnalysis(settings, countData);
      final LocalList<LocalList<PeakResult>> burstLocalisations =
          createBurstLocalisations(clusters, bursts);
      clusters.clear();
      burstLocalisations.forEach(list -> {
        final ClusterData d = new ClusterData(clusters.size() + 1, list);
        // Save this for analysis
        d.sourceRoi = roi;
        d.getArea();
        clusters.add(d);
      });
      return clusters;
    }).collect(LocalList::new, LocalList::addAll, LocalList::addAll);

    // Reorder
    final Counter count = new Counter();
    allClusters.forEach(c -> c.id = count.incrementAndGet());

    // Display in a table
    final ClusterDataTableModelFrame frame = createAllClustersTable();
    frame.getModel().setData(allClusters, dataCalibration);

    // Allow the results to be repeated
    frame.selectedAction = clusters -> {
      // Expecting a single cluster.
      final ClusterData c = clusters.get(0);
      // Push the correct ROI and settings to the analysis action.
      // We do not directly update the ROI or dialog settings as
      // these trigger events that are processed to add work with a delay.
      // Updating them at the end should generate events that are
      // ignored when finally executed as the ROI/settings should be the same.
      addWork(0, c.sourceRoi, settings, () -> {
        // When analysis has finished update the settings and image ROI.
        image.getImagePlus().setRoi(c.sourceRoi);
        darkTimeToleranceTextField.setText(Integer.toString(settings.getDarkTimeTolerance()));
        minClusterSizeTextField.setText(Integer.toString(settings.getMinClusterSize()));
        // When analysis has finished the cluster should be selected in the
        // current clusters table.
        final ClusterDataTableModelFrame currentClusters = currentClustersTable.get();
        if (currentClusters != null) {
          currentClusters.select(c);
        }
      });
    };

    // Show histogram of cluster size/duration
    reportAnalysis(settings, allClusters, dataCalibration);

    // Save clusters to memory
    final Trace[] traces = allClusters.stream().map(c -> {
      final Trace t = new Trace();
      t.setId(c.id);
      c.results.forEach(t::add);
      return t;
    }).toArray(Trace[]::new);
    TraceMolecules.saveResults(results, traces, "TC PALM");

    IJ.showStatus(TITLE + ": " + TextUtils.pleural(allClusters.size(), "cluster"));
  }

  /**
   * Determine if any areas overlap. Uses only the bounding rectangles
   *
   * @param rois the rois
   * @return true if an overlap exists
   */
  private static boolean anyOverlap(LocalList<Roi> rois) {
    for (int i = 0; i < rois.size(); i++) {
      final Rectangle2D.Double bounds = rois.get(i).getFloatBounds();
      for (int j = i + 1; j < rois.size(); j++) {
        if (bounds.intersects(rois.get(j).getFloatBounds())) {
          return true;
        }
      }
    }
    return false;
  }

  /**
   * Report statistics on the analysis results.
   *
   * @param settings the settings
   * @param clusters the clusters
   * @param calibration the data calibration
   */
  private static void reportAnalysis(TcPalmAnalysisSettings settings,
      LocalList<ClusterData> clusters, DataCalibration calibration) {
    final WindowOrganiser wo = new WindowOrganiser();
    final Consumer<HistogramPlotBuilder> action = builder -> {
      /* noop. */
    };
    plotHistogram(settings.getShowSizeHistogram(), wo, clusters, "Size", c -> c.results.size(),
        null, action);
    plotHistogram(settings.getShowDurationHistogram(), wo, clusters,
        "Duration (" + calibration.getTimeUnitName() + ")",
        c -> calibration.timeConverter.convert(c.getDuration()), null, action);
    plotHistogram(settings.getShowAreaHistogram(), wo, clusters,
        "Area (" + calibration.getDistanceUnitName() + "^2)",
        c -> calibration.convertArea(c.getArea()), null, action);
    plotHistogram(settings.getShowDensityHistogram(), wo, clusters,
        "Density (" + calibration.getDistanceUnitName() + "^-2)",
        c -> c.results.size() / calibration.convertArea(c.getArea()), Double::isFinite, action);
    wo.tile();
  }

  /**
   * Plot a histogram of the extracted statistic.
   *
   * @param wo the window organiser
   * @param clusters the clusters
   * @param name the name
   * @param function the function to extract the plotted statistic
   * @param predicate the predicate to filter the stream of data
   * @param action the action to use on the histogram plot (to set non-standard options)
   */
  private static void plotHistogram(boolean show, WindowOrganiser wo,
      LocalList<ClusterData> clusters, String name, ToDoubleFunction<ClusterData> function,
      DoublePredicate predicate, Consumer<HistogramPlotBuilder> action) {
    if (!show) {
      return;
    }
    final StoredData data = new StoredData(clusters.size());
    DoubleStream stream = clusters.stream().mapToDouble(function);
    if (predicate != null) {
      stream = stream.filter(predicate);
    }
    stream.forEach(data::add);
    final HistogramPlotBuilder builder = new HistogramPlotBuilder(TITLE, data, name);
    action.accept(builder);
    builder.show(wo);
  }
}
