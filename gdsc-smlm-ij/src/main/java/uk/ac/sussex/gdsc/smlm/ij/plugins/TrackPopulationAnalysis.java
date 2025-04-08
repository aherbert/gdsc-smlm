/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Package
 *
 * Software for single molecule localisation microscopy (SMLM) in ImageJ
 * %%
 * Copyright (C) 2011 - 2025 Alex Herbert
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
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.ImageWindow;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PointRoi;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.LUT;
import ij.text.TextWindow;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.Serializable;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Formatter;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Consumer;
import java.util.function.DoubleUnaryOperator;
import java.util.function.IntConsumer;
import java.util.function.IntFunction;
import java.util.function.IntUnaryOperator;
import java.util.function.Supplier;
import java.util.function.UnaryOperator;
import java.util.stream.IntStream;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.KeyStroke;
import javax.swing.RowSorter;
import javax.swing.SwingConstants;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.DefaultTableCellRenderer;
import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyIterationsException;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresFactory;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer.Optimum;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem.Evaluation;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.fitting.leastsquares.ParameterValidator;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.NonPositiveDefiniteMatrixException;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.commons.math3.util.Pair;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.UnitSphereSampler;
import org.apache.commons.statistics.distribution.FDistribution;
import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.core.data.VisibleForTesting;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.ij.BufferedTextWindow;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot.BinMethod;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot.HistogramPlotBuilder;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.MultiDialog;
import uk.ac.sussex.gdsc.core.ij.gui.NonBlockingExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.ScreenDimensionHelper;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper.LutColour;
import uk.ac.sussex.gdsc.core.logging.Ticker;
import uk.ac.sussex.gdsc.core.math.Mean;
import uk.ac.sussex.gdsc.core.math.SumOfSquaredDeviations;
import uk.ac.sussex.gdsc.core.utils.DoubleData;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.MemoryUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.SortUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.StoredData;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.ValidationUtils;
import uk.ac.sussex.gdsc.core.utils.rng.Mixers;
import uk.ac.sussex.gdsc.core.utils.rng.UniformRandomProviders;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.function.ChiSquaredDistributionTable;
import uk.ac.sussex.gdsc.smlm.ij.gui.TableColumnAdjuster;
import uk.ac.sussex.gdsc.smlm.math3.distribution.fitting.MultivariateGaussianMixtureExpectationMaximization;
import uk.ac.sussex.gdsc.smlm.math3.distribution.fitting.MultivariateGaussianMixtureExpectationMaximization.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.smlm.math3.distribution.fitting.MultivariateGaussianMixtureExpectationMaximization.MixtureMultivariateGaussianDistribution;
import uk.ac.sussex.gdsc.smlm.math3.distribution.fitting.MultivariateGaussianMixtureExpectationMaximization.MixtureMultivariateGaussianDistribution.MultivariateGaussianDistribution;
import uk.ac.sussex.gdsc.smlm.results.AttributePeakResult;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.PeakResultStoreList;
import uk.ac.sussex.gdsc.smlm.results.Trace;
import uk.ac.sussex.gdsc.smlm.results.procedures.PrecisionResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.sort.IdFramePeakResultComparator;

/**
 * Analyse tracks using a local sliding window to extract parameters that characterise the current
 * diffusion. Fit the parameters using a multi-variate Gaussian mixture model to identify
 * sub-populations such as bound and unbound particles. Output analysis on the sub-populations.
 *
 * <blockquote>Basu, et al (2020) Live-cell 3D single-molecule tracking reveals how NuRD modulates
 * enhancer dynamics. doi: <a href="https://doi.org/10.1101/2020.04.03.003178">bioRxiv
 * 2020.04.03.003178</a> </blockquote>
 */
public class TrackPopulationAnalysis implements PlugIn {
  private static final String TITLE = "Track Population Analysis";
  /** The feature names. */
  static final String[] FEATURE_NAMES = {"Anomalous exponent", "Effective diffusion coefficient",
      "Length of confinement", "Drift vector magnitude"};
  private static final String[] FEATURE_UNITS = {null, "μm^2/s", "μm", "μm"};
  /** The diffusion coefficient feature. */
  static final int FEATURE_D = 1;
  private static final int SORT_DIMENSION = 1;

  // Limits for fitting the MSD
  /** The minimum diffusion coefficient during fitting. */
  static final double MIN_D = Double.MIN_NORMAL;
  /** The minimum sigma during fitting. */
  static final double MIN_SIGMA = Double.MIN_NORMAL;
  /** The minimum alpha during fitting. */
  static final double MIN_ALPHA = Math.ulp(1.0);

  private static final AtomicReference<TextWindow> MODEL_TABLE_REF = new AtomicReference<>();

  /** The plugin settings. */
  private Settings settings;
  private boolean extraOptions;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> INSTANCE = new AtomicReference<>(new Settings());

    List<String> input;
    int window;
    int minTrackLength;
    int maxIterations;
    double relativeError;
    int repeats;
    int seed;
    boolean debug;
    int maxComponents;
    double minWeight;
    boolean simpleAlpha;
    boolean useT1Fraction;
    double minAlpha;
    double maxAlpha;
    double significance;
    boolean ignoreAlpha;
    int lutIndex;
    int histogramBins;
    // Track Data table
    boolean showTrackImage;
    boolean showTrackProbabilityPlot;
    boolean showSingleJumpPlot;
    boolean[] showFeaturePlot;
    double nmPerPixel;
    int minDisplaySize;
    boolean showMsdOverT;
    String dataSaveName;
    boolean dataSaveSelected;
    double angleMinJumpDistance;

    Settings() {
      // Set defaults
      window = 11;
      minTrackLength = 5;
      maxIterations = 1000;
      relativeError = 1e-6;
      repeats = 30;
      seed = 42;
      maxComponents = 2;
      minWeight = 0.1;
      maxAlpha = 2;
      significance = 0.05;
      ignoreAlpha = true;
      lutIndex = LutColour.RED_BLUE.ordinal();
      showTrackImage = true;
      showTrackProbabilityPlot = true;
      showSingleJumpPlot = true;
      showFeaturePlot = new boolean[FEATURE_NAMES.length];
      Arrays.fill(showFeaturePlot, true);
      nmPerPixel = 10;
      minDisplaySize = 400;
      angleMinJumpDistance = 125;
    }

    Settings(Settings source) {
      this.input = source.input;
      this.window = source.window;
      this.minTrackLength = source.minTrackLength;
      this.maxIterations = source.maxIterations;
      this.relativeError = source.relativeError;
      this.repeats = source.repeats;
      this.seed = source.seed;
      this.debug = source.debug;
      this.maxComponents = source.maxComponents;
      this.minWeight = source.minWeight;
      this.simpleAlpha = source.simpleAlpha;
      this.useT1Fraction = source.useT1Fraction;
      this.minAlpha = source.minAlpha;
      this.maxAlpha = source.maxAlpha;
      this.significance = source.significance;
      this.ignoreAlpha = source.ignoreAlpha;
      this.histogramBins = source.histogramBins;
      this.lutIndex = source.lutIndex;
      this.showTrackImage = source.showTrackImage;
      this.showTrackProbabilityPlot = source.showTrackProbabilityPlot;
      this.showSingleJumpPlot = source.showSingleJumpPlot;
      this.showFeaturePlot = source.showFeaturePlot.clone();
      this.nmPerPixel = source.nmPerPixel;
      this.minDisplaySize = source.minDisplaySize;
      this.showMsdOverT = source.showMsdOverT;
      this.dataSaveName = source.dataSaveName;
      this.dataSaveSelected = source.dataSaveSelected;
      this.angleMinJumpDistance = source.angleMinJumpDistance;
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
      return INSTANCE.get().copy();
    }

    /**
     * Save the settings.
     */
    void save() {
      INSTANCE.set(this);
    }
  }

  /**
   * Class to hold data for each track.
   */
  private static class TrackData {
    final int id;
    final int[] component;
    final double[][] data;
    final double[][] fitData;
    final Trace trace;
    BitSet bits;
    int transitionsCount = -1;

    /**
     * Create an instance.
     *
     * @param id the id
     * @param component the component
     * @param data the data
     * @param fitData the fit data (can be null if no fitting was done)
     * @param trace the trace
     */
    TrackData(int id, int[] component, double[][] data, double[][] fitData, Trace trace) {
      this.id = id;
      this.component = component;
      this.data = data;
      this.fitData = fitData;
      this.trace = trace;
    }

    /**
     * Gets the count of the different components.
     *
     * @return the component count
     */
    int getComponentCount() {
      return getComponents().cardinality();
    }

    /**
     * Gets the count of transitions between components.
     *
     * @return the transition count
     */
    int getTransitionCount() {
      int change = transitionsCount;
      if (change < 0) {
        change = 0;
        int current = component[0];
        for (int i = 1; i < component.length; i++) {
          if (current != component[i]) {
            current = component[i];
            change++;
          }
        }
        transitionsCount = change;
      }
      return change;
    }

    BitSet getComponents() {
      BitSet bits = this.bits;
      if (bits == null) {
        bits = new BitSet();
        for (final int c : component) {
          bits.set(c);
        }
        this.bits = bits;
      }
      return bits;
    }
  }

  /**
   * Class to show the track data in a JTable.
   */
  private static class TrackDataTableModel extends AbstractTableModel {
    private static final long serialVersionUID = 1L;

    List<TrackData> data;

    /**
     * Create an instance.
     *
     * @param trackData the track data
     */
    TrackDataTableModel(List<TrackData> trackData) {
      this.data = trackData;
    }

    @Override
    public int getRowCount() {
      return data.size();
    }

    @Override
    public int getColumnCount() {
      return 5;
    }

    @Override
    public String getColumnName(int columnIndex) {
      switch (columnIndex) {
        // @formatter:off
        case 0: return "ID";
        case 1: return "Length";
        case 2: return "Component count";
        case 3: return "Transitions";
        case 4: return "Components";
        // @formatter:on
        default:
          throw new IndexOutOfBoundsException("Bad column: " + columnIndex);
      }
    }

    @Override
    public Class<?> getColumnClass(int columnIndex) {
      switch (columnIndex) {
        case 0:
        case 1:
        case 2:
        case 3:
          return Integer.class;
        default:
          return String.class;
      }
    }

    @Override
    public Object getValueAt(int rowIndex, int columnIndex) {
      final TrackData c = data.get(rowIndex);
      switch (columnIndex) {
        // @formatter:off
        case 0: return c.id;
        case 1: return c.component.length;
        case 2: return c.getComponentCount();
        case 3: return c.getTransitionCount();
        case 4: return c.getComponents().toString();
        // @formatter:on
        default:
          throw new IndexOutOfBoundsException("Bad column: " + columnIndex);
      }
    }
  }

  /**
   * Class to display TrackDataTableModel in a JTable.
   */
  private static class TrackDataJTable extends JTable {
    private static final long serialVersionUID = 1L;

    /**
     * Create a new instance.
     *
     * @param model the model
     */
    TrackDataJTable(TrackDataTableModel model) {
      super(model);
      ((DefaultTableCellRenderer) getDefaultRenderer(Number.class))
          .setHorizontalAlignment(SwingConstants.TRAILING);
      setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
      setAutoCreateRowSorter(true);

      final TableColumnAdjuster tca = new TableColumnAdjuster(this, 6, false);
      // Only process 10 rows (5 at start, 5 at end).
      tca.setMaxRows(5);
      tca.setOnlyAdjustLarger(true);
      tca.adjustColumns();
    }

    /**
     * Returns the data of all selected rows. This maps the indices from the view to the data model.
     *
     * @return an array containing the data of all selected rows, or an empty array if no row is
     *         selected
     * @see #getSelectedRow
     */
    List<TrackData> getSelectedData() {
      final int iMin = selectionModel.getMinSelectionIndex();
      final int iMax = selectionModel.getMaxSelectionIndex();

      // Any negative
      if ((iMin | iMax) < 0) {
        return Collections.emptyList();
      }

      final LocalList<TrackData> rvTmp = new LocalList<>(1 + (iMax - iMin));

      final List<TrackData> data = ((TrackDataTableModel) dataModel).data;
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
      return sorter == null ? i -> i : sorter::convertRowIndexToModel;
    }
  }

  /**
   * Class to display TrackDataTableModel in a window frame.
   *
   * <p>Although this class extends {@link java.awt.Component} it is not intended to be
   * {@link Serializable}.
   */
  private static class TrackDataTableModelFrame extends JFrame implements ActionListener {
    private static final long serialVersionUID = 1L;
    /** Used to convert an angle in radians to turns. */
    private static final double FROM_RADIANS = 1.0 / (2 * Math.PI);

    private final Settings settings;
    /** The table. */
    final TrackDataJTable table;
    private final MixtureMultivariateGaussianDistribution model;
    private final Calibration cal;
    private final TypeConverter<DistanceUnit> distanceConverter;
    private final double deltaT;
    private final IntFunction<Color> colourMap;
    private JMenuItem dataSave;
    private JMenuItem analysisAnomalousExponentView;
    private JMenuItem analysisJumpAngles;
    private JMenuItem analysisFitJumpDistances;
    private JMenuItem analysisResidenceTime;
    private JMenuItem optionsTrackData;
    /** The selected action. */
    transient Consumer<List<TrackData>> selectedAction;

    /**
     * Create a new instance.
     *
     * @param settings the settings
     * @param tableModel the model
     * @param model the model (can be null if no fitting was done)
     * @param calibration the calibration
     * @param colourMap the colour map
     */
    TrackDataTableModelFrame(Settings settings, TrackDataTableModel tableModel,
        MixtureMultivariateGaussianDistribution model, Calibration calibration,
        final IntFunction<Color> colourMap) {
      this.settings = settings;
      this.model = model;
      cal = calibration;
      final CalibrationReader cr = new CalibrationReader(cal);
      distanceConverter = cr.getDistanceConverter(DistanceUnit.UM);
      deltaT = cr.getExposureTime() / 1000.0;
      this.colourMap = colourMap;

      setJMenuBar(createMenuBar());
      createSelectedAction();

      table = new TrackDataJTable(tableModel);
      final JScrollPane scroll = new JScrollPane(table);

      final ScreenDimensionHelper helper = new ScreenDimensionHelper();
      helper.setMinHeight(250);
      helper.setup(scroll);

      add(scroll);
      pack();

      addWindowListener(new WindowAdapter() {
        @Override
        public void windowOpened(WindowEvent event) {
          WindowManager.addWindow(TrackDataTableModelFrame.this);
        }

        @Override
        public void windowClosing(WindowEvent event) {
          WindowManager.removeWindow(TrackDataTableModelFrame.this);
        }
      });

      table.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
        @Override
        public void valueChanged(ListSelectionEvent e) {
          if (e.getValueIsAdjusting()) {
            return;
          }
          // Notify the action of the current selection
          selectedAction.accept(table.getSelectedData());
        }
      });
    }

    private JMenuBar createMenuBar() {
      final JMenuBar menubar = new JMenuBar();
      menubar.add(createDataMenu());
      menubar.add(createAnalysisMenu());
      menubar.add(createOptionsMenu());
      return menubar;
    }

    private JMenu createDataMenu() {
      final JMenu menu = new JMenu("Data");
      menu.setMnemonic(KeyEvent.VK_D);
      menu.add(dataSave = add("Save ...", KeyEvent.VK_S, "ctrl pressed S"));
      return menu;
    }

    private JMenu createAnalysisMenu() {
      final JMenu menu = new JMenu("Analysis");
      menu.setMnemonic(KeyEvent.VK_A);
      menu.add(analysisAnomalousExponentView =
          add("Anomalous exponent view", KeyEvent.VK_E, "ctrl pressed E"));
      menu.add(analysisJumpAngles = add("Jump angles", KeyEvent.VK_J, "ctrl pressed J"));
      menu.add(
          analysisFitJumpDistances = add("Fit jump distances", KeyEvent.VK_D, "ctrl pressed D"));
      menu.add(analysisResidenceTime = add("Residence Time", KeyEvent.VK_R, "ctrl pressed R"));
      return menu;
    }

    private JMenu createOptionsMenu() {
      final JMenu menu = new JMenu("Options");
      menu.setMnemonic(KeyEvent.VK_O);
      menu.add(optionsTrackData = add("Track Data ...", KeyEvent.VK_T, "ctrl pressed T"));
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
      Runnable action = null;
      if (src == optionsTrackData) {
        action = this::doOptionsTrackData;
      } else if (src == analysisAnomalousExponentView) {
        action = this::doAnalysisAnomalousExponentView;
      } else if (src == analysisJumpAngles) {
        action = this::doAnalysisJumpAngles;
      } else if (src == analysisFitJumpDistances) {
        action = this::doAnalysisFitJumpDistances;
      } else if (src == analysisResidenceTime) {
        action = this::doAnalysisResidenceTime;
      } else if (src == dataSave) {
        action = this::doDataSave;
      }
      // This will block the executing thread which may be the event dispatch thread.
      // Run on a new thread so events are still processed.
      if (action != null) {
        new Thread(action).start();
      }
    }

    /**
     * Show a dialog to select what to display for a selected track.
     */
    private void doOptionsTrackData() {
      // Show dialog to choose the track data
      final ExtendedGenericDialog gd = new ExtendedGenericDialog(getTitle() + " options");
      gd.addCheckbox("Show_track_image", settings.showTrackImage);
      gd.addNumericField("nm_per_pixel", settings.nmPerPixel);
      gd.addNumericField("Min_display_size", settings.minDisplaySize, 0);
      gd.addCheckbox("Show_track_probability_plot", settings.showTrackProbabilityPlot);
      for (int i = 0; i < FEATURE_NAMES.length; i++) {
        gd.addCheckbox("Show_" + FEATURE_NAMES[i] + "_plot", settings.showFeaturePlot[i]);
      }
      gd.addCheckbox("Show_single_jump_plot", settings.showSingleJumpPlot);
      gd.showDialog();
      if (gd.wasCanceled()) {
        return;
      }

      // So we know what to close
      final Settings oldSettings = settings.copy();

      settings.showTrackImage = gd.getNextBoolean();
      final double nmPerPixel = gd.getNextNumber();
      settings.minDisplaySize = (int) gd.getNextNumber();
      settings.showTrackProbabilityPlot = gd.getNextBoolean();
      for (int i = 0; i < FEATURE_NAMES.length; i++) {
        settings.showFeaturePlot[i] = gd.getNextBoolean();
      }
      settings.showSingleJumpPlot = gd.getNextBoolean();

      // Check arguments
      try {
        ParameterUtils.isAboveZero("nm/pixel", nmPerPixel);
      } catch (final IllegalArgumentException ex) {
        IJ.error(TITLE, ex.getMessage());
        return;
      }
      settings.nmPerPixel = nmPerPixel;

      // Close data no longer required
      close(oldSettings.showTrackImage, settings.showTrackImage, () -> TrackImage.TITLE);
      close(oldSettings.showTrackProbabilityPlot, settings.showTrackProbabilityPlot,
          () -> TrackProbabilityPlot.TITLE);
      for (int i = 0; i < FEATURE_NAMES.length; i++) {
        final int index = i;
        close(oldSettings.showFeaturePlot[i], settings.showFeaturePlot[i],
            () -> "Track " + FEATURE_NAMES[index]);
      }

      createSelectedAction();
    }

    /**
     * Close a window of the given title if it was previously visible and has been set to not
     * visible.
     *
     * @param wasVisible the previous visible flag
     * @param visible the current visible flag
     * @param title the title of the window
     */
    private static void close(boolean wasVisible, boolean visible, Supplier<String> title) {
      if (wasVisible && !visible) {
        final Window window = WindowManager.getWindow(title.get());
        if (window != null) {
          window.setVisible(false);
        }
      }
    }

    /**
     * Creates the action to perform on a selected track.
     */
    private void createSelectedAction() {
      // Avoid null pointer by always having an action
      selectedAction = NoopAction.INSTANCE;

      final WindowOrganiser wo = new WindowOrganiser();

      // Show probability of each component verses time.
      if (settings.showTrackProbabilityPlot && model != null) {
        selectedAction =
            selectedAction.andThen(new TrackProbabilityPlot(model, deltaT, colourMap, wo));
      }

      // Show plot of feature verses time.
      for (int i = 0; i < FEATURE_NAMES.length; i++) {
        if (settings.showFeaturePlot[i]) {
          selectedAction = selectedAction.andThen(new TrackFeaturePlot(deltaT, i, colourMap, wo));
        }
      }

      if (settings.showSingleJumpPlot) {
        selectedAction =
            selectedAction.andThen(new TrackDPlot(deltaT, distanceConverter, colourMap, wo));
      }

      // Show an image of the track coloured by the component.
      // Do this last as all the plot are the same size and will tile nicely.
      if (settings.showTrackImage) {
        selectedAction = selectedAction.andThen(new TrackImage(distanceConverter,
            settings.nmPerPixel, settings.minDisplaySize, colourMap, wo));
      }

      if (selectedAction != NoopAction.INSTANCE) {
        // Tile only new windows
        selectedAction = selectedAction.andThen(selected -> {
          wo.tile();
          wo.clear();
        });
      }
    }

    /**
     * For the currently selected track, extract the coordinates and show a dialog that allows the
     * mean squared jump distances to be plotted for a selected window from the track. The plot is
     * fit with the Brownian and Fractional Brownian Motion diffusion models.
     */
    private void doAnalysisAnomalousExponentView() {
      // Get the current selected track
      final List<TrackData> selected = table.getSelectedData();
      if (selected.isEmpty()) {
        return;
      }
      final TrackData trackData = selected.get(0);

      // Extract the track coordinates
      final Trace track = trackData.trace;
      final int size = track.size();
      final double[] x = new double[size];
      final double[] y = new double[size];
      for (int i = 0; i < size; i++) {
        final PeakResult peak = track.get(i);
        x[i] = distanceConverter.convert(peak.getXPosition());
        y[i] = distanceConverter.convert(peak.getYPosition());
      }

      final boolean extraOptions = ImageJUtils.isExtraOptions();
      final MsdFitter msdFitter = new MsdFitter(settings, deltaT);
      msdFitter.setData(x, y);
      final double[] t = SimpleArrayUtils.newArray(msdFitter.s.length, deltaT, deltaT);
      final SumOfSquaredDeviations[] msds = new SumOfSquaredDeviations[t.length];

      // For the selected anomalous exponent data point:
      final IntConsumer action = frame -> {
        // Extract the jump distances
        // Fit with regression, Brownian and FBM models
        try {
          msdFitter.fit(frame - 1, msds, false);
        } catch (TooManyIterationsException | ConvergenceException ex) {
          if (settings.debug) {
            ImageJUtils.log("Failed to fit anomalous exponent: " + ex.getMessage());
          }
          return;
        }

        // Plot the raw data and the fit results.
        final double d = msdFitter.getRegressionD();
        final double alpha = msdFitter.getRegressionAlpha();
        final RealVector p1 = msdFitter.lvmSolution1.getPoint();
        final RealVector p2 = msdFitter.lvmSolution2.getPoint();
        final String title = TITLE + " MSD vs Time";
        final Plot plot = new Plot(title, "time (s)", "MSD (um^^2^^)");
        final StringBuilder msg = new StringBuilder(512);
        try (Formatter formatter = new Formatter(msg)) {
          formatter.format("Brownian D=%s, s=%s : FBM D=%s, s=%s, alpha=%s (p-value=%s)",
              MathUtils.rounded(p1.getEntry(0)), MathUtils.rounded(p1.getEntry(1)),
              MathUtils.rounded(p2.getEntry(0)), MathUtils.rounded(p2.getEntry(1)),
              MathUtils.rounded(p2.getEntry(2)), MathUtils.rounded(msdFitter.pValue));
          if (extraOptions) {
            formatter.format(" : Simple D=%s, alpha=%s", MathUtils.rounded(d),
                MathUtils.rounded(alpha));
          }
        }
        plot.addLabel(0, 0, msg.toString());

        // Add error bars
        final double[] error = Arrays.stream(msds)
            .mapToDouble(ss -> ss.getStandardDeviation() / Math.sqrt(ss.getN())).toArray();

        // Convert the plot y axis by optionally dividing by the time and plot on a log scale.
        UnaryOperator<double[]> converter;
        if (settings.showMsdOverT) {
          plot.setLogScaleX();
          plot.setLogScaleY();
          plot.setXYLabels("time (s)", "MSD/t (um^^2^^/s)");
          converter = yy -> {
            for (int i = 0; i < yy.length; i++) {
              yy[i] /= t[i];
            }
            return yy;
          };
        } else {
          converter = UnaryOperator.identity();
        }

        plot.addPoints(t, converter.apply(msdFitter.s), converter.apply(error), Plot.CROSS);
        plot.addPoints(t, converter.apply(msdFitter.model2.value(p2).getFirst().toArray()),
            Plot.LINE);
        plot.setColor(Color.BLUE);
        plot.addPoints(t, converter.apply(msdFitter.model1.value(p1).getFirst().toArray()),
            Plot.LINE);
        plot.setColor(Color.RED);
        plot.addPoints(t,
            converter.apply(Arrays.stream(t).map(tt -> 4 * d * Math.pow(tt, alpha)).toArray()),
            Plot.LINE);
        // plot.setColor(Color.RED);
        // final double[] yy = Arrays.stream(t).map(msdFitter.reg::predict).toArray();
        // plot.addPoints(t, yy, Plot.CIRCLE);
        ImageJUtils.display(title, plot, ImageJUtils.NO_TO_FRONT);
      };
      action.accept(1);

      // Show a dialog with the number of anomalous exponents
      final NonBlockingExtendedGenericDialog gd =
          new NonBlockingExtendedGenericDialog(TITLE + " Anomalous Exponent View");
      gd.addMessage("Track " + trackData.id + " Anomalous Exponent View");
      gd.addSlider("Data window", 1, trackData.data.length, 1);
      gd.addCheckbox("Show MSD over time", settings.showMsdOverT);
      if (extraOptions) {
        gd.addCheckbox("Use_t1_fraction", settings.useT1Fraction);
      }
      gd.addDialogListener((egd, e) -> {
        // Read the frame
        final int frame = (int) egd.getNextNumber();
        settings.showMsdOverT = egd.getNextBoolean();
        if (extraOptions) {
          msdFitter.useT1Fraction = settings.useT1Fraction = egd.getNextBoolean();
        }
        // Show the data point
        action.accept(frame);
        return true;
      });
      gd.hideCancelButton();
      gd.setOKLabel("Close");
      gd.showDialog();
    }

    /**
     * For each component extract the jump angles between successive jumps using 3 consecutive
     * localisation positions. The component label before and after the central position is not
     * required to match.
     */
    private void doAnalysisJumpAngles() {
      final List<TrackData> data = getMultiSelectedDataOrAll();
      if (data.isEmpty()) {
        return;
      }
      // Exclude angles between jumps under a threshold. This excludes jumps that are
      // within the localisation precision for static molecules.
      // Show dialog to choose the track data
      final ExtendedGenericDialog gd = new ExtendedGenericDialog(getTitle() + " options");
      gd.addMessage("Specify the minimum distance for both jumps to include their angle.");
      gd.addSlider("Min_jump_distance (nm)", 50, 200, settings.angleMinJumpDistance);
      gd.showDialog();
      if (gd.wasCanceled()) {
        return;
      }
      settings.angleMinJumpDistance = gd.getNextNumber();

      final Int2ObjectOpenHashMap<DoubleArrayList> angleMap =
          getAngleMap(data, distanceConverter, settings.angleMinJumpDistance);
      final WindowOrganiser wo = new WindowOrganiser();
      // Histogram each angle bin from 0 to 360 in 15 degree increments (24 bins).
      // The limits max is inclusive so subtract a single bin width to ensure there is no empty
      // bin at [360, 375)
      final double binWidth = 15;
      final HistogramPlotBuilder builder = new HistogramPlotBuilder("title").setName("Jump Angle")
          .setLimits(new double[] {0, 360 - binWidth}).setMinBinWidth(binWidth)
          // This is a higher number than 360/binWidth to force use of the bin width
          .setNumberOfBins((int) (2 * 360 / binWidth));
      final int[] comps = angleMap.keySet().toIntArray();
      Arrays.sort(comps);
      for (final int comp : comps) {
        final DoubleArrayList angles = angleMap.get(comp);
        SimpleArrayUtils.apply(angles.elements(), 0, angles.size(), d -> d * 360);
        // @formatter:off
        final HistogramPlot plot = builder.setTitle(String.format("Component %d", comp))
            .setData(DoubleData.wrap(angles.toDoubleArray()))
            .build();
        // @formatter:on
        if (plot.show(wo) == null) {
          // Not enough data
          continue;
        }
        // Get the measure of asymmetry.
        // See Izeddin, et al (2014)
        // eLife 2014;3:e02230 DOI: 10.7554/eLife.02230
        // Figure 4
        final double[] y = plot.getPlotYValues();
        final double fwd = y[0] + y[1] + y[y.length - 1] + y[y.length - 2];
        final int mid = y.length / 2;
        final double bwd = y[mid] + y[mid + 1] + y[mid - 1] + y[mid - 2];
        // log2( (0 +/- 30) / (180 +/- 30) )
        final double asymmetryCoefficient = Math.log(fwd / bwd) / Math.log(2);
        final Plot p = plot.getPlot();
        p.addLabel(0, 0, String.format("Asymmetry Coefficient = %s (n = %d)",
            MathUtils.rounded(asymmetryCoefficient), angles.size()));
        p.setOptions("xinterval=45");
        p.updateImage();
      }
      wo.tile();
    }

    /**
     * Gets the angle map containing all the clockwise turns for each component. The turns are in
     * the range {@code [0, 1]}.
     *
     * @param data the data
     * @param distanceConverter the distance converter (converts native units to um)
     * @param minJumpDistance the minimum distance for each jump to include the angle (in nm)
     * @return the angle map
     */
    private static Int2ObjectOpenHashMap<DoubleArrayList> getAngleMap(List<TrackData> data,
        TypeConverter<DistanceUnit> distanceConverter, double minJumpDistance) {
      final Int2ObjectOpenHashMap<DoubleArrayList> angleMap = new Int2ObjectOpenHashMap<>();
      // Set the small jump distance using a squared threshold.
      // Convert min jump distance in nm to the native units.
      // Note: The distance converter is for micrometres.
      final double min2 = minJumpDistance <= 0 ? 0
          : MathUtils.pow2(distanceConverter.convertBack(minJumpDistance / 1000));
      final Int2IntOpenHashMap excluded = new Int2IntOpenHashMap();
      data.forEach(track -> {
        final int[] component = track.component;
        final Trace trace = track.trace;
        final int offset = (trace.size() - component.length) / 2;
        for (int i = 0; i < component.length; i++) {
          final PeakResult p1 = trace.get(i + offset - 1);
          final PeakResult p2 = trace.get(i + offset);
          final PeakResult p3 = trace.get(i + offset + 1);
          // Exclude small jump distances
          if (min2 > 0 && (p1.distance2(p2) < min2 || p2.distance2(p3) < min2)) {
            // Count the number of excluded angles
            excluded.addTo(component[i], 1);
            continue;
          }
          final double x1 = p1.getXPosition();
          final double y1 = p1.getYPosition();
          final double x2 = p2.getXPosition();
          final double y2 = p2.getYPosition();
          final double x3 = p3.getXPosition();
          final double y3 = p3.getYPosition();
          DoubleArrayList list = angleMap.get(component[i]);
          if (list == null) {
            angleMap.put(component[i], list = new DoubleArrayList());
          }
          // Turn angle between vector p1 -> p2 and p2 -> p3.
          // An angle of 0 is the same direction. 0.5 is a reverse.
          list.add(clockwiseTurns(angle(x1, y1, x2, y2), angle(x2, y2, x3, y3)));
        }
      });

      // Report the number excluded
      if (min2 > 0) {
        // Create the entire key set
        angleMap.keySet().forEach(key -> excluded.putIfAbsent(key, 0));
        final int[] comp = excluded.keySet().toIntArray();
        Arrays.sort(comp);
        ImageJUtils.log("Jump angle analysis: Min distance=%.1fnm", minJumpDistance);
        for (final int c : comp) {
          final int ex = excluded.get(c);
          final DoubleArrayList list = angleMap.get(c);
          final int total = ex + (list == null ? 0 : list.size());
          ImageJUtils.log("Component %d: Excluded %d / %d angles (%.2f%%)", c, ex, total,
              100.0 * ex / total);
        }
      }

      return angleMap;
    }

    /**
     * Compute the angle of the vector between the point and the origin.
     *
     * <p>Uses {@link Math#atan2(double, double)} thus increasing angle is for the
     * <em>counter-clockwise</em> direction.
     *
     * @param ox the origin x
     * @param oy the origin y
     * @param px the point x
     * @param py the point y
     * @return the angle
     */
    private static double angle(double ox, double oy, double px, double py) {
      final double dx = px - ox;
      final double dy = py - oy;
      return Math.atan2(dy, dx);
    }

    /**
     * Compute the number of turns <em>clockwise</em> to move from angle 1 to 2. The angles should
     * be in radians in the range -pi to pi with increasing angle for the <em>counter-clockwise</em>
     * direction (as computed by {@link Math#atan2(double, double)}).
     *
     * <p>The returned number of turns is in the range {@code [0, 1]}.
     *
     * @param angle1 the first angle
     * @param angle2 the second angle
     * @return the number of turns clockwise
     */
    private static double clockwiseTurns(double angle1, double angle2) {
      // Adapted from org.apache.commons.numbers.angle.PlaneAngle.normalize(...)
      // Convert radians to turns.
      // Expected range is [-0.5, 0.5] turns for each angle.
      // Result is in the range [-1, 1].
      final double turns = (angle1 - angle2) * FROM_RADIANS;
      // Normalise into the range [0, 1]
      return turns - Math.floor(turns);
    }

    /**
     * For each component extract the mean squared jump distances using all jump intervals up to the
     * configured window size. The component label must be the same for all coordinates used to
     * extract the jumps. A track must be at least the length of the window to be included in the
     * analysis. This reduces jump errors when the state of the component switches by reducing
     * exposure to transitions and eliminating short tracks.
     *
     * <p>Fit the combined MSD curve with a Brownian and FBM model.
     */
    private void doAnalysisFitJumpDistances() {
      final List<TrackData> data = getMultiSelectedDataOrAll();
      if (data.isEmpty()) {
        return;
      }

      // Configure options
      final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
      gd.addCheckbox("Show MSD over time", settings.showMsdOverT);
      gd.showDialog();
      if (gd.wasCanceled()) {
        return;
      }
      settings.showMsdOverT = gd.getNextBoolean();

      final Int2ObjectOpenHashMap<SumOfSquaredDeviations[]> msdMap = getMsdMap(data);
      final WindowOrganiser wo = new WindowOrganiser();
      final MsdFitter msdFitter = new MsdFitter(settings, deltaT);
      final double[] t = SimpleArrayUtils.newArray(msdFitter.s.length, deltaT, deltaT);
      final int[] comps = msdMap.keySet().toIntArray();
      Arrays.sort(comps);
      for (final int comp : comps) {
        final SumOfSquaredDeviations[] msds = msdMap.get(comp);
        // Do fitting
        try {
          msdFitter.fit(msds);
        } catch (TooManyIterationsException | ConvergenceException ex) {
          if (settings.debug) {
            ImageJUtils.log("Failed to fit anomalous exponent for component %d: ", comp,
                ex.getMessage());
          }
          return;
        }

        // Show the result
        final RealVector p1 = msdFitter.lvmSolution1.getPoint();
        final RealVector p2 = msdFitter.lvmSolution2.getPoint();
        final String title = String.format("Component %d MSD vs Time", comp);
        final Plot plot = new Plot(title, "time (s)", "MSD (um^^2^^)");
        // Note:
        // The number of windows is the count of samples of the longest jump.
        final int numberOfWindows = (int) msds[msds.length - 1].getN();
        // The number of tracks is the number of additional samples that can be obtained with
        // a jump size 1 smaller.
        final int numberOfTracks = (int) (msds[0].getN() - msds[1].getN());
        // This is the total number of data points minus the 3 fit parameters
        final int denominatorDegreesOfFreedom =
            (int) (Arrays.stream(msds).mapToLong(SumOfSquaredDeviations::getN).sum() - 3);
        plot.addLabel(0, 0,
            String.format(
                "Sub-tracks=%d, N=%d, DF=%d. Brownian D=%s, s=%s : "
                    + "FBM D=%s, s=%s, alpha=%s (p-value=%s)",
                numberOfTracks, numberOfWindows, denominatorDegreesOfFreedom,
                MathUtils.rounded(p1.getEntry(0)), MathUtils.rounded(p1.getEntry(1)),
                MathUtils.rounded(p2.getEntry(0)), MathUtils.rounded(p2.getEntry(1)),
                MathUtils.rounded(p2.getEntry(2)), MathUtils.rounded(msdFitter.pValue)));

        // Add error bars
        final double[] error = Arrays.stream(msds)
            .mapToDouble(ss -> ss.getStandardDeviation() / Math.sqrt(ss.getN())).toArray();

        // Convert the plot y axis by optionally dividing by the time and plot on a log scale.
        UnaryOperator<double[]> converter;
        if (settings.showMsdOverT) {
          plot.setLogScaleX();
          plot.setLogScaleY();
          plot.setXYLabels("time (s)", "MSD/t (um^^2^^/s)");
          converter = yy -> {
            for (int i = 0; i < yy.length; i++) {
              yy[i] /= t[i];
            }
            return yy;
          };
        } else {
          converter = UnaryOperator.identity();
        }

        plot.addPoints(t, converter.apply(msdFitter.s), converter.apply(error), Plot.CROSS);
        plot.addPoints(t, converter.apply(msdFitter.model2.value(p2).getFirst().toArray()),
            Plot.LINE);
        plot.setColor(Color.BLUE);
        plot.addPoints(t, converter.apply(msdFitter.model1.value(p1).getFirst().toArray()),
            Plot.LINE);
        ImageJUtils.display(title, plot, wo);
      }
      wo.tile();
    }

    /**
     * Gets the map containing the MSD for each component for each jump distance up to the window
     * size minus 1.
     *
     * @param data the data
     * @return the MSD map
     */
    private Int2ObjectOpenHashMap<SumOfSquaredDeviations[]> getMsdMap(List<TrackData> data) {
      final Int2ObjectOpenHashMap<SumOfSquaredDeviations[]> msdMap = new Int2ObjectOpenHashMap<>();
      data.forEach(track -> {
        final int[] component = track.component;
        final Trace trace = track.trace;
        // Window - 1
        final int wm1 = trace.size() - component.length;
        // Extract coordinates
        final double[] x = new double[component.length];
        final double[] y = new double[component.length];
        final int offset = wm1 / 2;
        for (int i = 0; i < x.length; i++) {
          final PeakResult peak = trace.get(i + offset);
          x[i] = distanceConverter.convert(peak.getXPosition());
          y[i] = distanceConverter.convert(peak.getYPosition());
        }
        // Find contiguous regions of the same component at least the length of the window
        int start = 0;
        while (start < component.length) {
          int end = start + 1;
          while (end < component.length && component[start] == component[end]) {
            end++;
          }
          addMsd(msdMap, component[start], start, end - 1, wm1, x, y);
          start = end;
        }
      });

      return msdMap;
    }

    /**
     * Compute the squared jump distances and add them to the means for the component.
     *
     * @param msdMap the MSD map
     * @param component the component
     * @param start the start (inclusive)
     * @param end the end (inclusive)
     * @param wm1 the window size minus 1
     * @param x the x coords
     * @param y the y coords
     */
    private static void addMsd(Int2ObjectOpenHashMap<SumOfSquaredDeviations[]> msdMap,
        int component, int start, int end, int wm1, double[] x, double[] y) {
      final int length = end - start;
      if (length < wm1) {
        return;
      }
      SumOfSquaredDeviations[] msds = msdMap.get(component);
      if (msds == null) {
        msds = SimpleArrayUtils.fill(new SumOfSquaredDeviations[wm1],
            (Supplier<SumOfSquaredDeviations>) SumOfSquaredDeviations::new);
        msdMap.put(component, msds);
      }
      for (int m = 1; m <= wm1; m++) {
        final SumOfSquaredDeviations m1 = msds[m - 1];
        for (int i = end - m; i >= start; i--) {
          final double d = MathUtils.distance2(x[i], y[i], x[i + m], y[i + m]);
          m1.add(d);
        }
      }
    }

    /**
     * For each component extract a histogram of the length of time a molecule is in each state.
     *
     * <p>Fit the combined MSD curve with a Brownian and FBM model.
     */
    private void doAnalysisResidenceTime() {
      final List<TrackData> data = getMultiSelectedDataOrAll();
      if (data.isEmpty()) {
        return;
      }

      final Int2ObjectOpenHashMap<IntArrayList> timeMap = getTimeMap(data);
      final WindowOrganiser wo = new WindowOrganiser();
      final int[] comps = timeMap.keySet().toIntArray();
      Arrays.sort(comps);
      final HistogramPlotBuilder builder =
          new HistogramPlotBuilder("title").setName("Residence Time").setMinBinWidth(deltaT);
      final UniformRandomProvider rng = UniformRandomProviders.create();
      for (final int comp : comps) {
        final int[] times = timeMap.get(comp).toIntArray();
        final int total = times.length;
        final double[] dtimes = Arrays.stream(times).mapToDouble(t -> t * deltaT).toArray();
        // Build histogram
        final String compName = String.format("Component %d", comp);
        final HistogramPlot plot =
            builder.setTitle(compName).setData(DoubleData.wrap(dtimes)).build();

        // Cumulative histogram
        final double[][] h = MathUtils.cumulativeHistogram(dtimes, false);
        final ExponentialDataFunction ef = ExponentialDataFunction.fromCdf(h);

        final double mean = ef.sumx / ef.nx;

        // Fit the residence time to: y = A exp(-t/Td)
        // with A = scaling factor, t = time and Td = Dissociation time
        // Using the histogram as a PMF approximation is error prone. Instead use
        // maximum likelihood fitting of the observations to an exponential.
        // Note:
        // Actual iterative fitting is not required as the log-likelihood gradient is zero
        // when the mean (t) is equal to sum(x) / n. See ExponentialDataFunction#gradient(double).
        // TODO: Try an alternative to maximum likelihood that optimises a different
        // goodness-of-fit function.
        // try {
        // final BrentSolver bs = new BrentSolver(1e-9, 1e-20);
        // final BracketFinder bf = new BracketFinder();
        // final int maxEval = 3000;
        // final double initialMean = mean;
        // bf.search(ef::value, GoalType.MAXIMIZE, initialMean * 0.5, initialMean * 2);
        // mean = bs.solve(maxEval, ef::gradient, bf.getLo(), bf.getHi(), bf.getMid());
        // } catch (TooManyEvaluationsException | NoBracketingException ex) {
        // ImageJUtils.log("Failed to fit exponential distribution: " + ex.getMessage());
        // }

        plot.show(wo);

        // Plot CDF
        SimpleArrayUtils.multiply(h[1], 1.0 / total);
        final String title = compName + " Cumulative Residence Time";
        final Plot plot2 = new Plot(title, "Residence Time (s)", "CDF");
        plot2.addPoints(h[0], h[1], Plot.BAR);

        if (mean != 0) {
          // Bootstrap the data to get an estimate:
          // 1. arrange all values in an array (i.e. times)
          // 2. sample the same number of values from the array with replacement
          // 3. compute parameters of exponential (i.e. the mean)
          // 4. repeat n times and compute confidence interval from the mean and SD of parameter
          final Statistics ss = new Statistics();
          for (int i = 0; i < 100; i++) {
            long sum = 0;
            for (int n = 0; n < total; n++) {
              sum += times[rng.nextInt(total)];
            }
            // maximum likelihood estimate of the mean
            ss.add(deltaT * sum / total);
          }
          final double ci = ss.getConfidenceInterval(0.95);

          final String label =
              String.format("Td = %s (n = %d; %s +/- %s (95%% ci))", MathUtils.rounded(mean), total,
                  MathUtils.rounded(ss.getMean()), MathUtils.rounded(ci));
          plot.getPlot().addLabel(0, 0, label);
          plot2.addLabel(0, 0, label);

          // Create the exponential CDF
          final double meanT = mean;
          final DoubleUnaryOperator cdf = x -> 1 - Math.exp(-x / meanT);

          // Show exponential on histogram:
          // y = 1/T exp(-x/T)
          // This must be rescaled by the number of observations. This does not work as the
          // histogram is a poor representation of the PMF. So instead we use the CDF to compute
          // the expected number of observations in the range of each bin.
          final double[] x = plot.getPlotXValues();
          final double dx2 = (x[1] - x[0]) / 2;
          double[] y = Arrays.stream(x)
              .map(xi -> total * (cdf.applyAsDouble(xi + dx2) - cdf.applyAsDouble(xi - dx2)))
              .toArray();
          plot.getPlot().setColor(Color.red);
          plot.getPlot().addPoints(x, y, Plot.LINE);
          plot.getPlot().updateImage();

          // Show exponential on cumulative histogram
          // y = 1 - exp(-x/T)
          final double x0 = h[0][0];
          final double dxx = (h[0][h[0].length - 1] - x0) / 100;
          final double[] xx =
              IntStream.rangeClosed(0, 100).mapToDouble(i -> x0 + i * dxx).toArray();
          y = Arrays.stream(xx).map(cdf).toArray();
          plot2.setColor(Color.red);
          plot2.addPoints(xx, y, Plot.LINE);
        }

        ImageJUtils.display(title, plot2, wo);
      }
      wo.tile();
    }

    /**
     * Gets the map containing the residence time in frames for each component.
     *
     * <p>Ignore the first and last component in each track, thus only processes lengths with a
     * defined start and end.
     *
     * @param data the data
     * @return the time map
     */
    private static Int2ObjectOpenHashMap<IntArrayList> getTimeMap(List<TrackData> data) {
      final Int2ObjectOpenHashMap<IntArrayList> timeMap = new Int2ObjectOpenHashMap<>();
      data.forEach(track -> {
        final int[] component = track.component;
        int current = component[0];
        int start = 0;
        for (int i = 1; i < component.length; i++) {
          if (current != component[i]) {
            // Ignore first component as it does not have a defined start point
            if (start != 0) {
              addTime(timeMap, current, start, i);
            }
            current = component[i];
            start = i;
          }
        }
        // Ignore last component as it does not have a defined end point
      });

      return timeMap;
    }

    /**
     * Add the residence time in the component to the time map.
     *
     * @param timeMap the time map
     * @param component the component
     * @param start the start (inclusive)
     * @param end the end (exclusive)
     */
    private static void addTime(Int2ObjectOpenHashMap<IntArrayList> timeMap, int component,
        int start, int end) {
      IntArrayList times = timeMap.get(component);
      if (times == null) {
        times = new IntArrayList();
        timeMap.put(component, times);
      }
      times.add(end - start);
    }

    /**
     * Save the selected data to a new dataset.
     */
    private void doDataSave() {
      // Get the currently selected track(s), or else all the data
      final List<TrackData> selected = table.getSelectedData();
      // Get the name for the dataset
      final GenericDialog gd = new GenericDialog(TITLE);
      gd.addMessage("Save track data to a results set");
      gd.addStringField("Name", settings.dataSaveName, 30);
      if (!selected.isEmpty()) {
        gd.addCheckbox("Save_selected", settings.dataSaveSelected);
      }
      gd.showDialog();
      if (gd.wasCanceled()) {
        return;
      }
      final String name = gd.getNextString();
      boolean saveSelected = false;
      if (!selected.isEmpty()) {
        settings.dataSaveSelected = saveSelected = gd.getNextBoolean();
      }
      if (TextUtils.isNullOrEmpty(name)) {
        return;
      }
      settings.dataSaveName = name;
      final List<TrackData> data =
          saveSelected ? selected : ((TrackDataTableModel) table.getModel()).data;

      // This may be a combined dataset with overlapping traces from the same frames.
      // Save with simple calibration (i.e. no PSF or full calibration).
      final MemoryPeakResults newResults = new MemoryPeakResults();
      newResults.setName(name);
      // This will pass through the camera calibration:
      // newResults.setCalibration(cal);
      // The input may be multiple datasets so we do not assume the same camera model.
      // The only requirement was the same exposure time, distance unit and nm/pixel.
      final CalibrationWriter cw = new CalibrationWriter();
      final CalibrationReader cr = new CalibrationReader(cal);
      cw.setDistanceUnit(cr.getDistanceUnit());
      cw.setNmPerPixel(cr.getNmPerPixel());
      cw.setTimeUnit(cr.getTimeUnit());
      cw.setExposureTime(cr.getExposureTime());
      newResults.setCalibration(cw.getCalibration());
      data.forEach(track -> {
        // Add categories using:
        // track.component
        // The window means that some of the coordinates at the start and end are not classified.
        // Find the offset to the first classified point.
        final PeakResultStoreList points = track.trace.getPoints();
        final int[] component = track.component;
        final int offset = (points.size() - component.length) / 2;
        int i = 0;
        // No category before the component categories
        while (i < offset) {
          newResults.add(new AttributePeakResult(points.get(i++)));
        }
        // Set the component using 1-based indexing
        for (final int c : component) {
          final AttributePeakResult r = new AttributePeakResult(points.get(i++));
          r.setCategory(c + 1);
          newResults.add(r);
        }
        // No category after the component categories
        while (i < points.size()) {
          newResults.add(new AttributePeakResult(points.get(i++)));
        }
      });
      MemoryPeakResults.addResults(newResults);
    }

    /**
     * Gets the selected tracks if a multi-selection, or else all tracks.
     *
     * @return the tracks
     */
    private List<TrackData> getMultiSelectedDataOrAll() {
      final List<TrackData> selected = table.getSelectedData();
      if (selected.size() > 1) {
        return selected;
      }
      return ((TrackDataTableModel) table.getModel()).data;
    }
  }

  /**
   * Class use to compute the log likelihood of exponentially distributed data.
   */
  @VisibleForTesting
  static class ExponentialDataFunction {
    /** The sum of all observations x. */
    final double sumx;
    /** The total number of observations. */
    final int nx;

    /**
     * Create an instance.
     *
     * @param x the x observations
     * @param count the count of each observation
     */
    ExponentialDataFunction(double[] x, int[] count) {
      // y = 1/T exp(-x/T)
      // log(y) = -x/T - log(T)
      // Precompute sum of x and the total count
      double sumx = 0;
      int n = 0;
      for (int i = 0; i < x.length; i++) {
        sumx += x[i] * count[i];
        n += count[i];
      }
      this.sumx = sumx;
      this.nx = n;
    }

    /**
     * Create from a cumulative distribution function. The CDF must not be normalised.
     *
     * @param h the CDF ([x, count])
     * @return the exponential function
     * @see MathUtils#cumulativeHistogram(double[], boolean)
     */
    static ExponentialDataFunction fromCdf(double[][] h) {
      final int[] counts = new int[h[1].length];
      // Convert cumulative frequency to counts
      int previous = 0;
      for (int i = 0; i < counts.length; i++) {
        final int value = (int) h[1][i];
        counts[i] = value - previous;
        previous = value;
      }
      return new ExponentialDataFunction(h[0], counts);
    }

    /**
     * Compute the log-likelihood for the given mean.
     *
     * @param t the mean
     * @return the log-likelihood
     */
    double value(double t) {
      // y = 1/T exp(-x/T)
      // log(y) = -x/T - log(T)
      // This can be partially precomputed
      return -sumx / t - nx * Math.log(t);
    }

    /**
     * Compute the first derivative (gradient) of the log-likelihood for the given mean.
     *
     * @param t the mean
     * @return the gradient of the log-likelihood
     */
    double gradient(double t) {
      // Note:
      // The gradient will be zero when t is the mean: sum(x) / n
      // sumx / (t * t) - n / t
      // = sumx / (sumx^2 / n^2) - n / (sumx / n)
      // = sumx * n^2 / sumx^2 - n^2 / sumx
      // = n^2 / sumx - n^2 / sumx
      // = 0

      return sumx / (t * t) - nx / t;
    }
  }

  /**
   * No operation implementation of the consumer interface.
   */
  private static class NoopAction implements Consumer<List<TrackData>> {
    static final NoopAction INSTANCE = new NoopAction();

    @Override
    public void accept(List<TrackData> t) {
      // Do nothing
    }

    @SuppressWarnings("unchecked")
    @Override
    public Consumer<List<TrackData>> andThen(Consumer<? super List<TrackData>> after) {
      return (Consumer<List<TrackData>>) after;
    }
  }

  /**
   * Show a plot of the selected track feature.
   */
  private static class TrackImage implements Consumer<List<TrackData>> {
    static final String TITLE = "Track Image";
    private static final BasicStroke DASHED = new BasicStroke(1, BasicStroke.CAP_SQUARE,
        BasicStroke.JOIN_MITER, 10.0f, new float[] {3, 3}, 0.0f);
    final double nmPerPixel;
    final int minDisplaySize;
    final double scale;
    final IntFunction<Color> colourMap;
    final WindowOrganiser wo;

    /**
     * Create an instance.
     *
     * @param distanceConverter the distance converter to get the coordinates in micrometers
     * @param nmPerPixel the nm per pixel for the output track image
     * @param minDisplaySize the min display size for the image window (in pixels)
     * @param colourMap the colour map
     * @param wo the window organiser
     */
    TrackImage(TypeConverter<DistanceUnit> distanceConverter, double nmPerPixel, int minDisplaySize,
        IntFunction<Color> colourMap, WindowOrganiser wo) {
      this.nmPerPixel = nmPerPixel;
      this.minDisplaySize = minDisplaySize;
      // Get the scale factor
      final double unitsToUm = distanceConverter.convert(1);
      // Create a distance converter to convert the units to the desired nm per pixel
      this.scale = unitsToUm * 1000 / nmPerPixel;
      this.colourMap = colourMap;
      this.wo = wo;
    }

    @Override
    public void accept(List<TrackData> tracks) {
      if (tracks.isEmpty()) {
        return;
      }
      final TrackData track = tracks.get(0);
      final Trace trace = track.trace;
      final int length = track.component.length;
      final int points = trace.size();
      final double[] x = new double[points];
      final double[] y = new double[points];
      for (int i = 0; i < points; i++) {
        x[i] = (float) (trace.get(i).getXPosition() * scale);
        y[i] = (float) (trace.get(i).getYPosition() * scale);
      }
      final double[] xlimits = MathUtils.limits(x);
      final double[] ylimits = MathUtils.limits(y);
      // Add border
      final float border = 1;
      xlimits[0] = Math.floor(xlimits[0] - border);
      ylimits[0] = Math.floor(ylimits[0] - border);
      xlimits[1] = Math.ceil(xlimits[1] + border);
      ylimits[1] = Math.ceil(ylimits[1] + border);
      final int width = (int) (xlimits[1] - xlimits[0]);
      final int height = (int) (ylimits[1] - ylimits[0]);
      // Blank white image
      final ByteProcessor bp = new ByteProcessor(width, height);
      bp.setValue(255);
      bp.fill();
      // Offset the coordinates
      for (int i = 0; i < points; i++) {
        x[i] -= xlimits[0];
        y[i] -= ylimits[0];
      }
      // Create a track overlay
      final Overlay overlay = new Overlay();
      // The window means that some of the coordinates at the start and end are not classified.
      // Find the offset to the first classified point.
      final int offset = (points - length) / 2;
      // Add the unclassified start/end as dashed black lines.
      int n = 0;
      final float[] xx = new float[x.length];
      final float[] yy = new float[x.length];
      for (int i = 0; i <= offset; i++) {
        xx[n] = (float) x[i];
        yy[n] = (float) y[i];
        n++;
      }
      Roi roi = createRoi(xx, yy, n, Color.LIGHT_GRAY);
      roi.setStroke(DASHED);
      overlay.add(roi);
      n = 0;
      for (int i = offset + length - 1; i < points; i++) {
        xx[n] = (float) x[i];
        yy[n] = (float) y[i];
        n++;
      }
      roi = createRoi(xx, yy, n, Color.GRAY);
      roi.setStroke(DASHED);
      overlay.add(roi);
      // Add the track end
      final PointRoi end = new PointRoi(x[points - 1], y[points - 1]);
      end.setStrokeColor(Color.GRAY);
      end.setPointType(PointRoi.CIRCLE);
      overlay.add(end);

      // Add the classified points using the correct colour
      final int[] component = track.component;
      int current = component[0];
      xx[0] = (float) x[offset];
      yy[0] = (float) y[offset];
      n = 1;
      for (int i = 1; i < length; i++) {
        if (current != component[i]) {
          // Add ROI
          overlay.add(createRoi(xx, yy, n, current));
          current = component[i];
          xx[0] = (float) x[i + offset - 1];
          yy[0] = (float) y[i + offset - 1];
          n = 1;
        }
        xx[n] = (float) x[i + offset];
        yy[n] = (float) y[i + offset];
        n++;
      }
      // Add final ROI
      overlay.add(createRoi(xx, yy, n, current));

      // Show the image.
      final ImagePlus imp = ImageJUtils.display(TITLE, bp, ImageJUtils.NO_TO_FRONT, wo);
      imp.setIgnoreGlobalCalibration(true);
      final ij.measure.Calibration cal = imp.getCalibration();
      cal.pixelWidth = cal.pixelHeight = nmPerPixel;
      cal.setUnit("nm");
      // Use a pixel offset to output the correct coordinate values
      cal.xOrigin = -xlimits[0] / nmPerPixel;
      cal.yOrigin = -ylimits[0] / nmPerPixel;
      imp.setOverlay(overlay);
      // Zoom in
      final ImageWindow iw = imp.getWindow();
      for (int i = 9; i-- > 0 && Math.max(iw.getWidth(), iw.getHeight()) < minDisplaySize;) {
        iw.getCanvas().zoomIn(imp.getWidth() / 2, imp.getHeight() / 2);
      }
    }

    private Roi createRoi(float[] x, float[] y, int n, int current) {
      return createRoi(x, y, n, colourMap.apply(current));
    }

    private static Roi createRoi(float[] x, float[] y, int n, Color c) {
      final PolygonRoi roi = new PolygonRoi(x, y, n, Roi.POLYLINE);
      roi.setStrokeColor(c);
      return roi;
    }
  }

  /**
   * Show a plot of the selected track feature.
   */
  private static class TrackProbabilityPlot implements Consumer<List<TrackData>> {
    static final String TITLE = "Track Probability";
    final MixtureMultivariateGaussianDistribution model;
    final double deltaT;
    final IntFunction<Color> colourMap;
    final WindowOrganiser wo;

    TrackProbabilityPlot(MixtureMultivariateGaussianDistribution model, double deltaT,
        IntFunction<Color> colourMap, WindowOrganiser wo) {
      this.model = model;
      this.deltaT = deltaT;
      this.colourMap = colourMap;
      this.wo = wo;
    }

    @Override
    public void accept(List<TrackData> tracks) {
      if (tracks.isEmpty()) {
        return;
      }
      final TrackData track = tracks.get(0);
      final int length = track.component.length;
      final float[] x = new float[length];
      for (int i = 0; i < length; i++) {
        x[i] = (float) ((i + 1) * deltaT);
      }
      final Plot plot = new Plot(TITLE, "Time (s)", "Probability");
      final double[] weights = model.getWeights();
      final MultivariateGaussianDistribution[] distributions = model.getDistributions();
      final float[][] y = new float[weights.length][x.length];
      final double[] p = new double[weights.length];
      for (int i = 0; i < x.length; i++) {
        double sum = 0;
        for (int n = 0; n < weights.length; n++) {
          p[n] = weights[n] * distributions[n].density(track.fitData[i]);
          sum += p[n];
        }
        for (int n = 0; n < weights.length; n++) {
          y[n][i] = (float) (p[n] / sum);
        }
      }
      for (int n = 0; n < weights.length; n++) {
        plot.setColor(colourMap.apply(n));
        plot.addPoints(x, y[n], Plot.LINE);
      }
      ImageJUtils.display(plot.getTitle(), plot, ImageJUtils.NO_TO_FRONT, wo).getPlot()
          .setLimitsToFit(true);
    }
  }

  /**
   * Store float coordinates.
   */
  private static final class Coords {
    float[] x;
    float[] y;
    int size;

    Coords(int length) {
      x = new float[length];
      y = new float[length];
    }

    void add(float v, float w) {
      if (x.length == size) {
        x = Arrays.copyOf(x, MemoryUtils.createNewCapacity(size + 1, size));
        y = Arrays.copyOf(y, x.length);
      }
      x[size] = v;
      y[size] = w;
      size++;
    }

    void clear() {
      size = 0;
    }

    float[] getX() {
      return Arrays.copyOf(x, size);
    }

    float[] getY() {
      return Arrays.copyOf(y, size);
    }
  }

  /**
   * Show a plot of the selected track feature.
   */
  private static class TrackFeaturePlot implements Consumer<List<TrackData>> {
    final double deltaT;
    final int featureIndex;
    final IntFunction<Color> colourMap;
    final WindowOrganiser wo;

    TrackFeaturePlot(double deltaT, int featureIndex, IntFunction<Color> colourMap,
        WindowOrganiser wo) {
      this.deltaT = deltaT;
      this.featureIndex = featureIndex;
      this.colourMap = colourMap;
      this.wo = wo;
    }

    @Override
    public void accept(List<TrackData> tracks) {
      if (tracks.isEmpty()) {
        return;
      }
      final TrackData track = tracks.get(0);
      final int length = track.component.length;
      final float[] x = new float[length];
      final float[] y = new float[length];
      for (int i = 0; i < length; i++) {
        x[i] = (float) ((i + 1) * deltaT);
        y[i] = (float) track.data[i][featureIndex];
      }
      final Plot plot = new Plot("Track " + FEATURE_NAMES[featureIndex], "Time (s)",
          getFeatureLabel(featureIndex));
      // Get the number of different components.
      final BitSet bits = track.getComponents();
      final int n = bits.cardinality();
      final int[] component = track.component;
      // Add each component using the correct colour.
      if (n == 1) {
        plot.setColor(colourMap.apply(component[0]));
        plot.addPoints(x, y, Plot.CIRCLE);
        plot.addPoints(x, y, Plot.LINE);
      } else {
        final Coords c = new Coords(11);
        bits.stream().forEach(com -> {
          plot.setColor(colourMap.apply(com));
          c.clear();
          for (int i = 0; i < x.length; i++) {
            if (component[i] == com) {
              c.add(x[i], y[i]);
              // Add lines using the colour of the second component.
              if (i != 0) {
                plot.drawLine(x[i - 1], y[i - 1], x[i], y[i]);
              }
            }
          }
          plot.addPoints(c.getX(), c.getY(), Plot.CIRCLE);
        });
      }
      ImageJUtils.display(plot.getTitle(), plot, ImageJUtils.NO_TO_FRONT, wo).getPlot()
          .setLimitsToFit(true);
    }
  }

  /**
   * Show a plot of the effective diffusion coefficient from a single jump, i.e. not over the
   * analysis window.
   */
  private static class TrackDPlot implements Consumer<List<TrackData>> {
    final double deltaT;
    final double scale;
    final IntFunction<Color> colourMap;
    final WindowOrganiser wo;

    TrackDPlot(double deltaT, TypeConverter<DistanceUnit> distanceConverter,
        IntFunction<Color> colourMap, WindowOrganiser wo) {
      this.deltaT = deltaT;
      // Extract effective diffusion coefficient in um^2/s:
      // d = MSD / (4 * deltaT)
      scale = MathUtils.pow2(distanceConverter.convert(1)) / (4 * deltaT);
      this.colourMap = colourMap;
      this.wo = wo;
    }

    @Override
    public void accept(List<TrackData> tracks) {
      if (tracks.isEmpty()) {
        return;
      }
      final TrackData track = tracks.get(0);
      // Extract effective diffusion coefficient in um^2/s:
      // d = MSD / (4 * deltaT)
      final Trace trace = track.trace;
      final int length = track.component.length;
      final int points = trace.size();
      // The window means that some of the coordinates at the start and end are not classified.
      // Find the offset to the first classified point.
      // Reduce by 1 as we have jumps and not raw XY coordinates.
      // jump[i + offset] == component[i]
      final int offset = ((points - length) / 2) - 1;
      final float[] t = new float[points - 1];
      final float[] d = new float[points - 1];
      double x = trace.get(0).getXPosition();
      double y = trace.get(0).getYPosition();
      for (int i = 1; i < points; i++) {
        final double x1 = trace.get(i).getXPosition();
        final double y1 = trace.get(i).getYPosition();
        // Offset the time so the classified points match with the other plots
        t[i - 1] = (float) ((i - offset) * deltaT);
        d[i - 1] = (float) (scale * (MathUtils.pow2(x1 - x) + MathUtils.pow2(y1 - y)));
        x = x1;
        y = y1;
      }

      final Plot plot = new Plot("Track " + FEATURE_NAMES[FEATURE_D] + " (single jump)", "Time (s)",
          getFeatureLabel(FEATURE_D));
      // Get colours for the different components.
      final int[] component = track.component;
      final Color[] map = new Color[MathUtils.max(component) + 1];
      final Color[] colours = new Color[t.length];
      // Grey for no classification before the components
      for (int i = 0; i < offset; i++) {
        colours[i] = Color.GRAY;
      }
      // Colour the classified points
      for (int i = 0; i < component.length; i++) {
        Color c = map[component[i]];
        if (c == null) {
          c = map[component[i]] = colourMap.apply(component[i]);
        }
        colours[i + offset] = c;
      }
      // Grey for no classification after the components
      for (int i = component.length + offset; i < colours.length; i++) {
        colours[i] = Color.GRAY;
      }
      // Add lines using the colour of the second component.
      for (int i = 1; i < t.length; i++) {
        plot.setColor(colours[i]);
        plot.drawLine(t[i - 1], d[i - 1], t[i], d[i]);
      }
      // Add points of the same colour
      final Coords coords = new Coords(t.length);
      for (final Color c : colours) {
        if (c == null) {
          continue;
        }
        coords.clear();
        plot.setColor(c);
        for (int i = 0; i < component.length; i++) {
          final int j = i + offset;
          if (colours[j] == c) {
            coords.add(t[j], d[j]);
          }
        }
        plot.addPoints(coords.getX(), coords.getY(), Plot.CIRCLE);
      }
      ImageJUtils.display(plot.getTitle(), plot, ImageJUtils.NO_TO_FRONT, wo).getPlot()
          .setLimitsToFit(true);
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    extraOptions = ImageJUtils.isExtraOptions();

    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "No localisations in memory");
      return;
    }

    settings = Settings.load();
    // Saved by reference so just save now
    settings.save();

    // Read in multiple traced datasets
    // All datasets must have the same pixel pitch and exposure time
    // Get parameters
    // Convert datasets to tracks
    // For each track compute the 4 local track features using the configured window
    //
    // Optional:
    // Fit a multi-variate Gaussian mixture model to the data
    // (using the configured number of components/populations)
    // Assign each point in the track using the model.
    // Smooth the assignments.
    //
    // The alternative is to use the localisation category to assign populations.
    //
    // Plot histograms of each track parameter, coloured by component

    final List<MemoryPeakResults> combinedResults = new LocalList<>();

    if (!showInputDialog(combinedResults)) {
      return;
    }

    final boolean hasCategory = showHasCategoryDialog(combinedResults);

    if (!showDialog(hasCategory)) {
      return;
    }

    ImageJUtils.log(TITLE + "...");

    final List<Trace> tracks = getTracks(combinedResults, settings.window, settings.minTrackLength);
    if (tracks.isEmpty()) {
      IJ.error(TITLE, "No tracks. Please check the input data and min track length setting.");
      return;
    }

    final Calibration cal = combinedResults.get(0).getCalibration();
    final CalibrationReader cr = new CalibrationReader(cal);
    // Use micrometer / second
    final TypeConverter<DistanceUnit> distanceConverter = cr.getDistanceConverter(DistanceUnit.UM);
    final double exposureTime = cr.getExposureTime() / 1000.0;
    final Pair<int[], double[][]> trackData =
        extractTrackData(tracks, distanceConverter, exposureTime, hasCategory);
    final double[][] data = trackData.getValue();

    // Histogram the raw data.
    final Array2DRowRealMatrix raw = new Array2DRowRealMatrix(data, false);
    final WindowOrganiser wo = new WindowOrganiser();
    // Store the histogram data for plotting the components
    final double[][] columns = new double[FEATURE_NAMES.length][];
    final double[][] limits = new double[FEATURE_NAMES.length][];
    // Get column data
    for (int i = 0; i < FEATURE_NAMES.length; i++) {
      columns[i] = raw.getColumn(i);
      if (i == FEATURE_D) {
        // Plot using a logarithmic scale
        SimpleArrayUtils.apply(columns[i], Math::log10);
      }
      limits[i] = MathUtils.limits(columns[i]);
    }
    // Compute histogram bins
    final int[] bins = new int[FEATURE_NAMES.length];
    if (settings.histogramBins > 0) {
      Arrays.fill(bins, settings.histogramBins);
    } else {
      for (int i = 0; i < FEATURE_NAMES.length; i++) {
        bins[i] = HistogramPlot.getBins(StoredData.create(columns[i]), BinMethod.FD);
      }
      // Use the maximum so all histograms look the same
      Arrays.fill(bins, MathUtils.max(bins));
    }
    // Compute plots
    final Plot[] plots = new Plot[FEATURE_NAMES.length];
    for (int i = 0; i < FEATURE_NAMES.length; i++) {
      final double[][] hist =
          HistogramPlot.calcHistogram(columns[i], limits[i][0], limits[i][1], bins[i]);
      plots[i] =
          new Plot(TITLE + " " + FEATURE_NAMES[i], getFeatureLabel(i, i == FEATURE_D), "Frequency");
      plots[i].addPoints(hist[0], hist[1], Plot.BAR);
      ImageJUtils.display(plots[i].getTitle(), plots[i], ImageJUtils.NO_TO_FRONT, wo);
    }
    wo.tile();

    // The component for each data point
    int[] component;
    // The number of components
    int numComponents;

    // Data used to fit the Gaussian mixture model
    double[][] fitData;
    // The fitted model
    MixtureMultivariateGaussianDistribution model;

    if (hasCategory) {
      // Use the category as the component.
      // No fit data and no output model
      fitData = null;
      model = null;

      // The component is stored at the end of the raw track data.
      final int end = data[0].length - 1;
      component = Arrays.stream(data).mapToInt(d -> (int) d[end]).toArray();
      numComponents = MathUtils.max(component) + 1;

      // In the EM algorithm the probability of each data point is computed and normalised to
      // sum to 1. The normalised probabilities are averaged to create the weights.
      // Note the probability of each data point uses the previous weight and the algorithm
      // iterates.
      // This is not a fitted model but the input model so use
      // zero weights to indicate no fitting was performed.
      final double[] weights = new double[numComponents];
      // Remove the trailing component to show the 'model' in a table.
      createModelTable(Arrays.stream(data).map(d -> Arrays.copyOf(d, end)).toArray(double[][]::new),
          weights, component);
    } else {
      // Multivariate Gaussian mixture EM

      // Provide option to not use the anomalous exponent in the population mix.
      int sortDimension = SORT_DIMENSION;
      if (settings.ignoreAlpha) {
        // Remove index 0. This shifts the sort dimension.
        sortDimension--;
        fitData = Arrays.stream(data).map(d -> Arrays.copyOfRange(d, 1, d.length))
            .toArray(double[][]::new);
      } else {
        fitData = SimpleArrayUtils.deepCopy(data);
      }

      final MultivariateGaussianMixtureExpectationMaximization mixed =
          fitGaussianMixture(fitData, sortDimension);

      if (mixed == null) {
        IJ.error(TITLE, "Failed to fit a mixture model");
        return;
      }

      model = sortComponents(mixed.getFittedModel(), sortDimension);

      // For the best model, assign to the most likely population.
      component = assignData(fitData, model);

      // Table of the final model using the original data (i.e. not normalised)
      final double[] weights = model.getWeights();
      numComponents = weights.length;
      createModelTable(data, weights, component);
    }

    // Output coloured histograms of the populations.
    final LUT lut = LutHelper.createLut(settings.lutIndex);
    IntFunction<Color> colourMap;
    if (LutHelper.getColour(lut, 0).equals(Color.BLACK)) {
      colourMap = i -> LutHelper.getNonZeroColour(lut, i, 0, numComponents - 1);
    } else {
      colourMap = i -> LutHelper.getColour(lut, i, 0, numComponents - 1);
    }
    for (int i = 0; i < FEATURE_NAMES.length; i++) {
      // Extract the data for each component
      final double[] col = columns[i];
      final Plot plot = plots[i];
      for (int n = 0; n < numComponents; n++) {
        final StoredData feature = new StoredData();
        for (int j = 0; j < component.length; j++) {
          if (component[j] == n) {
            feature.add(col[j]);
          }
        }
        if (feature.size() == 0) {
          continue;
        }
        final double[][] hist =
            HistogramPlot.calcHistogram(feature.values(), limits[i][0], limits[i][1], bins[i]);
        // Colour the points
        plot.setColor(colourMap.apply(n));
        plot.addPoints(hist[0], hist[1], Plot.BAR);
      }
      plot.updateImage();
    }

    createTrackDataTable(tracks, trackData, fitData, model, component, cal, colourMap);

    // Analysis.
    // Assign the original localisations to their track component.
    // Q. What about the start/end not covered by the window?
    // Save tracks as a dataset labelled with the sub-track ID?

    // Output for the bound component and free components track parameters.
    // Compute dwell times.
    // Other ...

    // Track analysis plugin:
    // Extract all continuous segments of the same component.
    // Produce MSD plot with error bars.
    // Fit using FBM model.
  }

  /**
   * Creates a table to show the final model. This uses the assignments to create a mixture model
   * from the original data.
   *
   * @param data the data
   * @param weights the weights for each component
   * @param component the component
   */
  private static void createModelTable(double[][] data, double[] weights, int[] component) {
    final MixtureMultivariateGaussianDistribution model =
        MultivariateGaussianMixtureExpectationMaximization.createMixed(data, component);
    final MultivariateGaussianDistribution[] distributions = model.getDistributions();

    // Get the fraction of each component
    final int[] count = new int[MathUtils.max(component) + 1];
    Arrays.stream(component).forEach(c -> count[c]++);

    try (BufferedTextWindow tw = new BufferedTextWindow(ImageJUtils.refresh(MODEL_TABLE_REF,
        () -> new TextWindow("Track Population Model", createHeader(), "", 800, 300)))) {
      final StringBuilder sb = new StringBuilder();
      for (int i = 0; i < weights.length; i++) {
        sb.setLength(0);
        sb.append(i).append('\t');
        sb.append(MathUtils.rounded((double) count[i] / component.length)).append('\t');
        sb.append(MathUtils.rounded(weights[i]));
        final double[] means = distributions[i].getMeans();
        final double[] sd = distributions[i].getStandardDeviations();
        for (int j = 0; j < means.length; j++) {
          sb.append('\t').append(MathUtils.rounded(means[j])).append('\t')
              .append(MathUtils.rounded(sd[j]));
        }
        tw.append(sb.toString());
      }
    }
  }

  private static String createHeader() {
    final StringBuilder sb = new StringBuilder("Component\tFraction\tGaussian Weight");
    for (int i = 0; i < FEATURE_NAMES.length; i++) {
      sb.append('\t').append(getFeatureLabel(i)).append("\t+/-");
    }
    return sb.toString();
  }

  /**
   * Creates the track data table.
   *
   * @param tracks the tracks
   * @param trackData the track data
   * @param fitData the fit data (can be null if no fitting was done)
   * @param model the model (can be null if no fitting was done)
   * @param component the component
   * @param calibration the calibration
   * @param distanceConverter the distance converter (to get the coordinates in micrometers)
   * @param deltaT the time step of each frame in seconds (delta T)
   * @param colourMap the colour map
   */
  private void createTrackDataTable(List<Trace> tracks, Pair<int[], double[][]> trackData,
      double[][] fitData, MixtureMultivariateGaussianDistribution model, int[] component,
      Calibration calibration, IntFunction<Color> colourMap) {
    // Get the track data
    final LocalList<TrackData> list = new LocalList<>(tracks.size());
    final int[] lengths = trackData.getFirst();
    final double[][] data = trackData.getSecond();
    int from = 0;
    for (int i = 0; i < lengths.length; i++) {
      final int to = from + lengths[i];
      list.add(new TrackData(i, Arrays.copyOfRange(component, from, to),
          Arrays.copyOfRange(data, from, to),
          fitData == null ? null : Arrays.copyOfRange(fitData, from, to), tracks.get(i)));
      from = to;
    }

    // Track viewer:
    // Show table of tracks. Summarise:
    // length, set of components, count of component switches
    final TrackDataTableModelFrame frame = new TrackDataTableModelFrame(settings,
        new TrackDataTableModel(list), model, calibration, colourMap);
    frame.setTitle(TITLE + " Track Data");
    // Single selection is required for plotting. But some actions can work on multiple selection.
    // frame.table.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    frame.setVisible(true);
  }

  private boolean showInputDialog(List<MemoryPeakResults> combinedResults) {
    // Show a list box containing all the clustered results.
    // This should remember the last set of chosen items.
    final MultiDialog md = ResultsManager.createMultiDialog(TITLE, MemoryPeakResults::hasId);
    md.setSelected(settings.input);
    md.setHelpUrl(HelpUrls.getUrl("track-population-analysis"));

    md.showDialog();

    if (md.wasCancelled()) {
      return false;
    }

    final List<String> selected = md.getSelectedResults();
    if (selected.isEmpty()) {
      IJ.error(TITLE, "No results were selected");
      return false;
    }
    settings.input = selected;

    for (final String name : selected) {
      final MemoryPeakResults r = MemoryPeakResults.getResults(name);
      if (r != null) {
        combinedResults.add(r);
      }
    }

    // Check calibration exists for the first set of results
    if (combinedResults.isEmpty() || !checkCalibration(combinedResults.get(0))) {
      return false;
    }

    // Check the calibration is the same for the rest
    final CalibrationReader cal = combinedResults.get(0).getCalibrationReader();
    final double nmPerPixel = cal.getNmPerPixel();
    final double exposureTime = cal.getExposureTime();
    final DistanceUnit distanceUnit = cal.getDistanceUnit();
    for (int i = 1; i < combinedResults.size(); i++) {
      final MemoryPeakResults results = combinedResults.get(i);

      if (!results.hasCalibration()
          || results.getCalibrationReader().getExposureTime() != exposureTime
          || results.getNmPerPixel() != nmPerPixel || results.getDistanceUnit() != distanceUnit) {
        IJ.error(TITLE,
            "The exposure time, pixel pitch and distance unit must match across all the results");
        return false;
      }
    }

    return true;
  }

  /**
   * Check the results have a calibrated exposure time and pixel pitch. If not then show a dialog to
   * collect the calibration.
   *
   * @param results the results
   * @return True if calibrated
   */
  private static boolean checkCalibration(MemoryPeakResults results) {
    if (results.getCalibration() == null || !results.getCalibrationReader().hasExposureTime()
        || !results.getCalibrationReader().hasNmPerPixel()) {
      final CalibrationWriter cal = results.getCalibrationWriterSafe();

      final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
      gd.addMessage("Uncalibrated results! Please enter the calibration:");
      gd.addNumericField("Exposure_time", cal.getExposureTime(), 2, 6, "ms");
      gd.addNumericField("Pixel_pitch", cal.getNmPerPixel(), 2, 6, "nm");
      gd.showDialog();
      if (gd.wasCanceled() || gd.invalidNumber()) {
        return false;
      }
      cal.setExposureTime(gd.getNextNumber());
      cal.setNmPerPixel(gd.getNextNumber());
      if (cal.getExposureTime() <= 0 || cal.getNmPerPixel() <= 0) {
        return false;
      }
      results.setCalibration(cal.getCalibration());
    }
    return true;
  }

  private static boolean showHasCategoryDialog(List<MemoryPeakResults> combinedResults) {
    // Check if all results have a category
    if (combinedResults.stream().filter(MemoryPeakResults::hasCategory).count() == combinedResults
        .size()) {
      final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
      gd.addHelp(HelpUrls.getUrl("track-population-analysis"));
      gd.addMessage(
          "Input results contain categories.\n \nUse the categories as the track sub-populations?");
      gd.enableYesNoCancel();
      gd.hideCancelButton();
      gd.showDialog();
      return gd.wasOKed();
    }
    return false;
  }

  private boolean showDialog(boolean hasCategory) {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("track-population-analysis"));

    gd.addSlider("Window", 5, 20, settings.window);
    gd.addSlider("Min_track_length", 2, 20, settings.minTrackLength);
    gd.addMessage("Anomalous diffusion coefficient");
    if (extraOptions) {
      gd.addCheckbox("Simple_alpha", settings.simpleAlpha);
      gd.addCheckbox("Use_t1_fraction", settings.useT1Fraction);
    }
    gd.addNumericField("Fit_significance", settings.significance, -2);
    gd.addNumericField("Min_alpha", settings.minAlpha, -3);
    gd.addNumericField("Max_alpha", settings.maxAlpha, -3);
    if (!hasCategory) {
      gd.addMessage("Multi-variate Gaussian mixture Expectation-Maximisation");
      gd.addCheckbox("Ignore_alpha", settings.ignoreAlpha);
      gd.addSlider("Max_components", 2, 10, settings.maxComponents);
      gd.addSlider("Min_weight", 0.01, 1, settings.minWeight);
      gd.addNumericField("Max_iterations", settings.maxIterations, 0);
      gd.addNumericField("Relative_error", settings.relativeError, -1);
      gd.addNumericField("Repeats", settings.repeats, 0);
      gd.addNumericField("Seed", settings.seed, 0);
      gd.addCheckbox("Debug", settings.debug);
    }
    gd.addMessage("Output options");
    gd.addNumericField("Histogram_bins", settings.histogramBins, 0);
    gd.addChoice("LUT", LutHelper.getLutNames(), settings.lutIndex);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.window = (int) gd.getNextNumber();
    settings.minTrackLength = (int) gd.getNextNumber();
    if (extraOptions) {
      settings.simpleAlpha = gd.getNextBoolean();
      settings.useT1Fraction = gd.getNextBoolean();
    }
    settings.significance = gd.getNextNumber();
    settings.minAlpha = gd.getNextNumber();
    settings.maxAlpha = gd.getNextNumber();
    if (!hasCategory) {
      settings.ignoreAlpha = gd.getNextBoolean();
      settings.maxComponents = (int) gd.getNextNumber();
      settings.minWeight = gd.getNextNumber();
      settings.maxIterations = (int) gd.getNextNumber();
      settings.relativeError = gd.getNextNumber();
      settings.repeats = (int) gd.getNextNumber();
      settings.seed = (int) gd.getNextNumber();
      settings.debug = gd.getNextBoolean();
    }
    settings.histogramBins = (int) gd.getNextNumber();
    settings.lutIndex = gd.getNextChoiceIndex();

    if (gd.invalidNumber()) {
      return false;
    }

    // Check arguments
    try {
      // For fitting the alpha coefficient there should be:
      // (number of parameters) < (number of points)
      // where the number of points is the window size minus 1.
      ParameterUtils.isEqualOrAbove("Window", settings.window, 5);
      ParameterUtils.isEqualOrAbove("Min track length", settings.minTrackLength, 2);
      if (!hasCategory) {
        ParameterUtils.isEqualOrAbove("Max components", settings.maxComponents, 2);
        ParameterUtils.isPositive("Min weight", settings.minWeight);
        ParameterUtils.isAboveZero("Max iterations", settings.maxIterations);
        ParameterUtils.isAboveZero("Maximum N", settings.relativeError);
        ParameterUtils.isAboveZero("Repeats", settings.repeats);
      }
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return true;
  }

  /**
   * Gets the tracks. Each track has contiguous frames and the length is enough to fit
   * {@code minTrackLength} overlapping windows of the specified size:
   *
   * <pre>
   * length >= window + minTrackLength - 1
   * </pre>
   *
   * @param combinedResults the combined results
   * @param window the window size
   * @param minTrackLength the minimum track length (assumed to be {@code >= 1})
   * @return the tracks
   */
  private static List<Trace> getTracks(List<MemoryPeakResults> combinedResults, int window,
      int minTrackLength) {
    final LocalList<Trace> tracks = new LocalList<>();
    final Statistics stats = new Statistics();
    final int minSize = window + Math.max(minTrackLength, 1) - 1;
    combinedResults.forEach(results -> {
      final int start = tracks.size();
      // Sort by id then frame
      results = results.copy();
      results.sort(IdFramePeakResultComparator.INSTANCE);
      final int size = results.size();
      // Skip IDs not associated with clustering
      int index = 0;
      while (index < size && results.get(index).getId() < 1) {
        index++;
      }
      // Initialise current id and frame
      int id = results.get(index).getId() - 1;
      int frame = results.get(index).getFrame();
      Trace track = new Trace();
      for (; index < size; index++) {
        final PeakResult result = results.get(index);
        // Same ID and contiguous frames
        if (result.getId() != id || result.getFrame() != frame + 1) {
          addTrack(minSize, tracks, track);
          track = new Trace();
        }
        id = result.getId();
        frame = result.getFrame();
        track.add(result);
      }
      addTrack(minSize, tracks, track);

      stats.reset();
      for (int i = start; i < tracks.size(); i++) {
        stats.add(tracks.unsafeGet(i).size());
      }

      final StringBuilder sb = new StringBuilder(256);
      TextUtils.formatTo(sb, "%s tracks=%d, length=%s +/- %s", results.getName(), stats.getN(),
          MathUtils.rounded(stats.getMean(), 3),
          MathUtils.rounded(stats.getStandardDeviation(), 3));

      // Limit of diffusion coefficient from the localisation precision.
      // Just use the entire dataset for simplicity (i.e. not the tracks of min length).
      final PrecisionResultProcedure pp = new PrecisionResultProcedure(results);
      try {
        pp.getPrecision();
        final Mean mean = new Mean();
        for (final double p : pp.precisions) {
          mean.add(p);
        }
        // 2nDt = MSD (n = number of dimensions)
        // D = MSD / 2nt
        final CalibrationReader reader = results.getCalibrationReader();
        final double t = reader.getExposureTime() / 1000.0;
        // Assume computed in nm. Convert to um.
        final double x = mean.getMean() / 1000;
        final double d = x * x / (2 * t);
        TextUtils.formatTo(sb, ", precision=%s nm, D limit=%s um^2/s",
            MathUtils.rounded(x * 1000, 4), MathUtils.rounded(d, 4));
      } catch (final DataException ignored) {
        // No precision
      }
      IJ.log(sb.toString());
    });
    return tracks;
  }

  /**
   * Adds the track to the list of tracks if it is {@code >= minimumSize}.
   *
   * @param minimumSize the minimum size
   * @param tracks the tracks
   * @param track the track
   */
  private static void addTrack(int minimumSize, final List<Trace> tracks, Trace track) {
    if (track.size() >= minimumSize) {
      tracks.add(track);
    }
  }

  /**
   * Extract the track data. This extracts different descriptors of the track using a rolling local
   * window.
   *
   * <p>Distances are converted to {@code unit} using the provided converter and time units are
   * converted from frame to seconds (s). The diffusion coefficients is in unit^2/s.
   *
   * <p>If categories are added they are remapped to be a natural sequence starting from 0.
   *
   * @param tracks the tracks
   * @param distanceConverter the distance converter
   * @param deltaT the time step of each frame in seconds (delta T)
   * @param hasCategory if true add the category from the result to the end of the data
   * @return the track data (track lengths and descriptors)
   */
  private Pair<int[], double[][]> extractTrackData(List<Trace> tracks,
      TypeConverter<DistanceUnit> distanceConverter, double deltaT, boolean hasCategory) {
    final List<double[]> data = new LocalList<>(tracks.size());
    double[] x = new double[0];
    double[] y = new double[0];
    final int wm1 = settings.window - 1;
    final int valuesLength = hasCategory ? 5 : 4;
    final int mid = settings.window / 2;

    final MsdFitter msdFitter = new MsdFitter(settings, deltaT);
    final double significance = settings.significance;
    final int[] fitResult = new int[4];
    final boolean regressionOnly = extraOptions && settings.simpleAlpha;

    // Factor for the diffusion coefficient: 1/N * 1/(2*dimensions*deltaT) = 1 / 4Nt
    // with N the number of points to average.
    final double diffusionCoefficientFactor = 1.0 / (4 * wm1 * deltaT);
    // Drift coefficient: 1/N * 1/deltaT
    final double driftMagnitudeFactor = 1.0 / (wm1 * deltaT);

    // Used for the standard deviations
    final Statistics statsX = new Statistics();
    final Statistics statsY = new Statistics();
    final Ticker ticker = ImageJUtils.createTicker(tracks.size(), 1, "Computing track features...");

    // Collect the track categories. Later these are renumbered to a natural sequence from 0.
    final IntOpenHashSet catSet = new IntOpenHashSet();

    // final StoredDataStatistics statsAlpha = new StoredDataStatistics(tracks.size());

    // Process each track
    final IntArrayList lengths = new IntArrayList(tracks.size());
    for (final Trace track : tracks) {
      // Get xy coordinates
      final int size = track.size();
      if (x.length < size) {
        x = new double[size];
        y = new double[size];
      }
      for (int i = 0; i < size; i++) {
        final PeakResult peak = track.get(i);
        x[i] = distanceConverter.convert(peak.getXPosition());
        y[i] = distanceConverter.convert(peak.getYPosition());
      }
      final int smwm1 = size - wm1;
      final int previousSize = data.size();
      for (int k = 0; k < smwm1; k++) {
        final double[] values = new double[valuesLength];
        data.add(values);

        // Append the category to the data. This works best for odd sized windows as the
        // middle of even sized windows is between two localisations.
        if (hasCategory) {
          final int cat = track.get(k + mid).getCategory();
          values[4] = cat;
          catSet.add(cat);
        }

        // First point in window = k
        // Last point in window = k + w - 1 = k + wm1
        final int end = k + wm1;

        // 1. Anomalous exponent.
        msdFitter.setData(x, y);
        try {
          msdFitter.fit(k, null, regressionOnly);
          if (regressionOnly) {
            values[0] = msdFitter.getRegressionAlpha();
          } else {

            int fitIndex = msdFitter.positiveSlope ? 0 : 2;

            // If better then this is anomalous diffusion
            final double alpha = msdFitter.lvmSolution2.getPoint().getEntry(2);
            values[0] = alpha;
            if (msdFitter.pValue > significance) {
              fitIndex++;
            }
            fitResult[fitIndex]++;

            // // Debug
            // if (
            // // msdFitter.pValue < 0.2
            // msdFitter.pValue < 0.2 && values[0] < 0
            // // !msdFitter.positiveSlope
            // ) {
            // final RealVector p = msdFitter.lvmSolution2.getPoint();
            // final String title = "anomalous exponent";
            // final Plot plot = new Plot(title, "time (s)", "MSD (um^^2^^)");
            // final double[] t = SimpleArrayUtils.newArray(msdFitter.s.length, deltaT, deltaT);
            // plot.addLabel(0, 0, msdFitter.lvmSolution2.getPoint().toString() + " p="
            // + msdFitter.pValue + ". " + msdFitter.lvmSolution1.getPoint().toString());
            // plot.addPoints(t, msdFitter.s, Plot.CROSS);
            // plot.addPoints(t, msdFitter.model2.value(p).getFirst().toArray(), Plot.LINE);
            // plot.setColor(Color.BLUE);
            // plot.addPoints(t,
            // msdFitter.model1.value(msdFitter.lvmSolution1.getPoint()).getFirst().toArray(),
            // Plot.LINE);
            // plot.setColor(Color.RED);
            // final double[] yy = Arrays.stream(t).map(msdFitter.reg::predict).toArray();
            // plot.addPoints(t, yy, Plot.CIRCLE);
            // System.out.printf("%s : %s", msdFitter.lvmSolution2.getPoint(), values[0]);
            // ImageJUtils.display(title, plot, ImageJUtils.NO_TO_FRONT);
            // System.out.println();
            // }
          }
        } catch (TooManyIterationsException | ConvergenceException ex) {
          if (settings.debug) {
            ImageJUtils.log("Failed to fit anomalous exponent: " + ex.getMessage());
          }
          // Ignore this and leave as Brownian motion
          values[0] = 1.0;
        }

        // Referenced papers:
        // Hozé, N. H., D. (2017) Statistical methods for large ensembles of super-resolution
        // stochastic single particle trajectories in cell biology.
        // Annual Review of Statistics and Its Application 4, 189-223
        //
        // Amitai, A., Seeber, A., Gasser, S. M. & Holcman, D. (2017) Visualization of Chromatin
        // Decompaction and Break Site Extrusion as Predicted by Statistical Polymer
        // Modeling of Single-Locus Trajectories. Cell reports 18, 1200-1214

        // 2. Effective diffusion coefficient (Hozé, eq 10).
        // This is the average squared jump distance between successive points
        // multiplied by 1 / (2 * dimensions * deltaT), i.e. 1 / 4t.
        double sum = 0;
        for (int i = k; i < end; i++) {
          sum += MathUtils.distance2(x[i], y[i], x[i + 1], y[i + 1]);
        }
        values[1] = sum * diffusionCoefficientFactor;

        // 3. Length of confinement (Amitai et al, eq 1).
        // Compute the average of the standard deviation of the position in each dimension.
        statsX.reset();
        statsY.reset();
        for (int i = k; i <= end; i++) {
          statsX.add(x[i]);
          statsY.add(y[i]);
        }
        values[2] = (statsX.getStandardDeviation() + statsY.getStandardDeviation()) / 2;

        // 4. Magnitude of drift vector (Hozé, eq 9).
        // Note: The drift field is given as the expected distance between successive points, i.e.
        // the average step. Since all track windows are the same length this is the same
        // as the distance between the first and last point divided by the number of points.
        // The drift field in each dimension is combined to create a drift norm, i.e. Euclidean
        // distance.
        // double sumx = 0;
        // double sumy = 0;
        // for (int i = k; i < end; i++) {
        // sumx += x[i + 1] - x[i];
        // sumy += y[i + 1] - y[i];
        // }
        // System.out.println(DoubleEquality.relativeError(sumx, x[end] - x[k]));
        // System.out.println(DoubleEquality.relativeError(sumy, y[end] - y[k]));
        values[3] = MathUtils.distance(x[k], y[k], x[end], y[end]) * driftMagnitudeFactor;
      }
      lengths.add(data.size() - previousSize);
      ticker.tick();
    }
    ImageJUtils.finished();
    if (settings.debug) {
      ImageJUtils.log("  +Slope, significant:   %d", fitResult[0]);
      ImageJUtils.log("  +Slope, insignificant: %d", fitResult[1]);
      ImageJUtils.log("  -Slope, significant:   %d", fitResult[2]);
      ImageJUtils.log("  -Slope, insignificant: %d", fitResult[3]);
    }
    ImageJUtils.log("Insignificant anomalous exponents: %d / %d", fitResult[1] + fitResult[3],
        data.size());
    // System.out.println(statsAlpha.getStatistics().toString());

    final double[][] trackData = data.toArray(new double[0][0]);
    if (hasCategory) {
      final int[] categories = catSet.toIntArray();
      Arrays.sort(categories);
      // Only remap if non-compact (i.e. not 0 to n)
      final int max = categories[categories.length - 1];
      if (categories[0] != 0 || max + 1 != categories.length) {
        final int[] remap = new int[max + 1];
        for (int i = 0; i < categories.length; i++) {
          remap[categories[i]] = i;
        }
        final int end = valuesLength - 1;
        for (final double[] values : trackData) {
          values[end] = remap[(int) values[end]];
        }
      }
    }
    return Pair.create(lengths.toIntArray(), trackData);
  }

  /**
   * Gets the residual sum of squares.
   *
   * @param lvmSolution1 the lvm solution
   * @return the residual sum of squares
   */
  static double getResidualSumOfSquares(final Optimum lvmSolution1) {
    final RealVector res = lvmSolution1.getResiduals();
    return res.dotProduct(res);
  }

  /**
   * Fit the Gaussian mixture to the data. The fitter with the highest likelihood from a number of
   * repeats is returned.
   *
   * @param data the data
   * @param sortDimension the sort dimension
   * @return the multivariate gaussian mixture
   */
  private MultivariateGaussianMixtureExpectationMaximization
      fitGaussianMixture(final double[][] data, int sortDimension) {
    // Get the unmixed multivariate Guassian.
    MultivariateGaussianDistribution unmixed =
        MultivariateGaussianMixtureExpectationMaximization.createUnmixed(data);

    // Normalise the columns of the data
    // Get means and SD of each column
    final double[] means = unmixed.getMeans();
    final double[] sd = unmixed.getStandardDeviations();
    final int dimensions = means.length;
    for (final double[] value : data) {
      for (int i = 0; i < dimensions; i++) {
        value[i] = (value[i] - means[i]) / sd[i];
      }
    }

    // Repeat. The mean should be approximately 0 and std.dev. 1.
    unmixed = MultivariateGaussianMixtureExpectationMaximization.createUnmixed(data);

    // Record the likelihood of the unmixed model
    double logLikelihood = Arrays.stream(data).mapToDouble(unmixed::density).map(Math::log).sum();
    // x means, x*x covariances
    final int parametersPerGaussian = dimensions + dimensions * dimensions;
    double aic = MathUtils.getAkaikeInformationCriterion(logLikelihood, parametersPerGaussian);
    double bic = MathUtils.getBayesianInformationCriterion(logLikelihood, data.length,
        parametersPerGaussian);
    ImageJUtils.log("1 component log-likelihood=%s. AIC=%s. BIC=%s", logLikelihood, aic, bic);

    // Fit a mixed component model.
    // Increment the number of components up to a maximim or when the model does not improve.
    MultivariateGaussianMixtureExpectationMaximization mixed = null;
    for (int numComponents = 2; numComponents <= settings.maxComponents; numComponents++) {
      final MultivariateGaussianMixtureExpectationMaximization mixed2 =
          createMixed(data, dimensions, numComponents);
      if (mixed2 == null) {
        ImageJUtils.log("Failed to fit a %d component mixture model", numComponents);
        break;
      }
      final double logLikelihood2 = mixed2.getLogLikelihood();
      // n * (means, covariances, 1 weight) - 1
      // (Note: subtract 1 as the weights are constrained by summing to 1)
      final int param2 = numComponents * (parametersPerGaussian + 1) - 1;
      final double aic2 = MathUtils.getAkaikeInformationCriterion(logLikelihood2, param2);
      final double bic2 =
          MathUtils.getBayesianInformationCriterion(logLikelihood2, data.length, param2);

      // Log-likelihood ratio test statistic
      final double lambdaLr = -2 * (logLikelihood - logLikelihood2);
      // DF = difference in dimensionality from previous number of components
      // means, covariances, 1 weight
      final int degreesOfFreedom = parametersPerGaussian + 1;
      final double q = ChiSquaredDistributionTable.computeQValue(lambdaLr, degreesOfFreedom);
      ImageJUtils.log("%d component log-likelihood=%s. AIC=%s. BIC=%s. LLR significance=%s.",
          numComponents, logLikelihood2, aic2, bic2, MathUtils.rounded(q));
      final double[] weights = mixed2.getFittedModel().getWeights();

      // For consistency sort the mixture by the mean of the diffusion coefficient
      final double[] values = Arrays.stream(mixed2.getFittedModel().getDistributions())
          .mapToDouble(d -> d.getMeans()[sortDimension]).toArray();
      SortUtils.sortData(weights, values, false, false);
      ImageJUtils.log("Population weights: " + Arrays.toString(weights));

      if (MathUtils.min(weights) < settings.minWeight) {
        ImageJUtils.log("%d component model has population weight %s under minimum level %s",
            numComponents, MathUtils.min(weights), settings.minWeight);
        break;
      }
      if (aic <= aic2 || bic <= bic2 || q > 0.001) {
        ImageJUtils.log("%d component model is not significant", numComponents);
        break;
      }
      aic = aic2;
      bic = bic2;
      logLikelihood = logLikelihood2;
      mixed = mixed2;
    }

    return mixed;
  }

  /**
   * Creates the multivariate gaussian mixture as the best of many repeats of the expectation
   * maximisation algorithm.
   *
   * @param data the data
   * @param dimensions the dimensions
   * @param numComponents the number of components
   * @return the multivariate gaussian mixture expectation maximization
   */
  private MultivariateGaussianMixtureExpectationMaximization createMixed(final double[][] data,
      int dimensions, int numComponents) {
    // Fit a mixed multivariate Gaussian with different repeats.
    final UnitSphereSampler sampler = UnitSphereSampler
        .of(UniformRandomProviders.create(Mixers.stafford13(settings.seed++)), dimensions);
    final LocalList<CompletableFuture<MultivariateGaussianMixtureExpectationMaximization>> results =
        new LocalList<>(settings.repeats);
    final DoubleDoubleBiPredicate test = createConvergenceTest(settings.relativeError);
    if (settings.debug) {
      ImageJUtils.log("  Fitting %d components", numComponents);
    }
    final Ticker ticker = ImageJUtils.createTicker(settings.repeats, 2, "Fitting...");
    final AtomicInteger failures = new AtomicInteger();
    for (int i = 0; i < settings.repeats; i++) {
      final double[] vector = sampler.sample();
      results.add(CompletableFuture.supplyAsync(() -> {
        final MultivariateGaussianMixtureExpectationMaximization fitter =
            new MultivariateGaussianMixtureExpectationMaximization(data);
        try {
          // This may also throw the same exceptions due to inversion of the covariance matrix
          final MixtureMultivariateGaussianDistribution initialMixture =
              MultivariateGaussianMixtureExpectationMaximization.estimate(data, numComponents,
                  point -> {
                    double dot = 0;
                    for (int j = 0; j < dimensions; j++) {
                      dot += vector[j] * point[j];
                    }
                    return dot;
                  });
          final boolean result = fitter.fit(initialMixture, settings.maxIterations, test);
          // Log the result. Note: The ImageJ log is synchronized.
          if (settings.debug) {
            ImageJUtils.log("  Fit: log-likelihood=%s, iter=%d, converged=%b",
                fitter.getLogLikelihood(), fitter.getIterations(), result);
          }
          return result ? fitter : null;
        } catch (NonPositiveDefiniteMatrixException | SingularMatrixException ex) {
          failures.getAndIncrement();
          if (settings.debug) {
            ImageJUtils.log("  Fit failed during iteration %d. No variance in a sub-population "
                + "component (check alpha is not always 1.0).", fitter.getIterations());
          }
        } finally {
          ticker.tick();
        }
        return null;
      }));
    }
    ImageJUtils.finished();
    if (failures.get() != 0 && settings.debug) {
      ImageJUtils.log("  %d component fit failed %d/%d", numComponents, failures.get(),
          settings.repeats);
    }
    // Collect results and return the best model.
    return results.stream().map(f -> f.join()).filter(f -> f != null)
        .sorted((f1, f2) -> Double.compare(f2.getLogLikelihood(), f1.getLogLikelihood()))
        .findFirst().orElse(null);
  }

  /**
   * Creates the convergence test.
   *
   * @param relativeError the relative error
   * @return the predicate
   */
  private static DoubleDoubleBiPredicate createConvergenceTest(double relativeError) {
    return (v1, v2) -> DoubleEquality.relativeError(v1, v2) < relativeError;
  }

  /**
   * Sort the components by the mean of the given dimension.
   *
   * @param model the model
   * @param dimension the dimension
   * @return the mixture multivariate gaussian distribution
   */
  private static MixtureMultivariateGaussianDistribution
      sortComponents(MixtureMultivariateGaussianDistribution model, int dimension) {
    final double[] weights = model.getWeights();
    final MultivariateGaussianDistribution[] distributions = model.getDistributions();
    final LocalList<Pair<Double, MultivariateGaussianDistribution>> list =
        new LocalList<>(weights.length);
    for (int i = 0; i < weights.length; i++) {
      list.add(Pair.create(weights[i], distributions[i]));
    }
    list.sort((o1, o2) -> Double.compare(o1.getSecond().getMeans()[dimension],
        o2.getSecond().getMeans()[dimension]));
    for (int i = 0; i < weights.length; i++) {
      weights[i] = list.unsafeGet(i).getFirst();
      distributions[i] = list.unsafeGet(i).getSecond();
    }
    return MixtureMultivariateGaussianDistribution.create(weights, distributions);
  }

  /**
   * Gets the feature label.
   *
   * @param feature the feature
   * @return the feature label
   */
  static String getFeatureLabel(int feature) {
    return getFeatureLabel(feature, false);
  }

  /**
   * Gets the feature label.
   *
   * @param feature the feature
   * @param log10 true if the units are logarithmic using log10
   * @return the feature label
   */
  private static String getFeatureLabel(int feature, boolean log10) {
    if (FEATURE_UNITS[feature] == null) {
      return FEATURE_NAMES[feature];
    }
    final StringBuilder sb = new StringBuilder(FEATURE_NAMES[feature]);
    if (log10) {
      sb.append(" (log(");
    } else {
      sb.append(" (");
    }
    sb.append(FEATURE_UNITS[feature]);
    if (log10) {
      sb.append("))");
    } else {
      sb.append(')');
    }
    return sb.toString();
  }

  /**
   * Assign the data using the mixture model.
   *
   * @param data the data
   * @param model the model
   * @return the assignments
   */
  private static int[] assignData(double[][] data, MixtureMultivariateGaussianDistribution model) {
    final double[] weights = model.getWeights();
    final MultivariateGaussianDistribution[] distributions = model.getDistributions();
    // All initialised as component 0
    final int[] comp = new int[data.length];
    for (int i = 0; i < comp.length; i++) {
      // Assign using the highest probability.
      final double[] x = data[i];
      double max = weights[0] * distributions[0].density(x);
      for (int j = 1; j < weights.length; j++) {
        final double p = weights[j] * distributions[j].density(x);
        if (max < p) {
          max = p;
          comp[i] = j;
        }
      }
    }
    // Neighbour smoothing window of 3 eliminates isolated assignments, e.g. C in between U:
    // U C U => U U U
    // Note: For more than two components this does nothing to isolated assignments
    // between different neighbours:
    // U C X => U C X
    int n = 0;
    if (comp.length >= 3) {
      int p0 = comp[0];
      int p1 = comp[1];
      for (int i = 2; i < comp.length; i++) {
        final int p2 = comp[i];
        if (p0 != p1 && p0 == p2) {
          comp[i - 1] = p1 = p2;
          n++;
        }
        p0 = p1;
        p1 = p2;
      }
    }
    ImageJUtils.log("Re-labelled isolated assignments: %d / %d", n, comp.length);
    return comp;
  }

  /**
   * Class to encapsulate the fitting of the mean square distance data.
   */
  private static class MsdFitter {
    /**
     * Small fraction of MSD for the first jump to set an anchor point close to t0. This fails to
     * handle data where the localisation precision is a large component of the MSD, i.e. MSD(t=0)
     * is non-zero. Setting a value very close to (0,0) can anchor the regression curve incorrectly
     * for non-diffusing data with a flat MSD curve.
     */
    private static final double T1_FRACTION = 1.0 / 50;
    private static final double LOG_T1_FRACTION = -Math.log(50);

    // CHECKSTYLE.OFF: MemberName
    final double deltaT;
    final int window;
    final int wm1;
    final double[] s;

    final LevenbergMarquardtOptimizer optimizer;
    final RealVector observed;
    final RealMatrix weightMatrix;
    final ConvergenceChecker<Evaluation> checker;

    // Linear model for Brownian motion
    final MultivariateJacobianFunction model1;
    final RealVector start1;
    final ParameterValidator paramValidator1;

    // Non-linear model for anomalous diffusion coefficient
    final MultivariateJacobianFunction model2;
    final RealVector start2;
    final ParameterValidator paramValidator2;

    /** Linear regression of: MSD = A t. Slope is 4D. */
    final SimpleRegression reg;
    /**
     * Linear regression of:
     *
     * <pre>
     * log(MSD) = log(A) + alpha * log(t)
     * for MSD = A t^alpha.
     * </pre>
     *
     * <p>Slope is anomalous coefficient alpha. exp(intercept) is 4D.
     */
    final SimpleRegression reg2;
    boolean useT1Fraction;
    final int denominatorDegreesOfFreedom;
    final FDistribution distribution;

    double[] x;
    double[] y;

    boolean positiveSlope;
    Optimum lvmSolution1;
    Optimum lvmSolution2;
    double pValue;
    // CHECKSTYLE.ON: MemberName

    // double alpha;

    /**
     * Create an instance.
     *
     * @param settings the settings
     * @param deltaT the time step of each frame in seconds (delta T)
     */
    MsdFitter(Settings settings, double deltaT) {
      this.deltaT = deltaT;
      window = settings.window;
      wm1 = window - 1;
      // Use for fitting.
      // This uses a weighted linear regression as the MSD(n) is a mean of individual distances
      // using the frame gap n.
      s = new double[wm1];
      optimizer = new LevenbergMarquardtOptimizer();
      observed = new ArrayRealVector(s, false);
      final double[] weight = SimpleArrayUtils.newArray(wm1, wm1, -1.0);
      weightMatrix = new DiagonalMatrix(weight, false);
      checker = (iteration, previous,
          current) -> DoubleEquality.relativeError(previous.getCost(), current.getCost()) < 1e-6;

      // Linear model for Brownian motion
      model1 = new BrownianDiffusionFunction(wm1, deltaT);
      start1 = new ArrayRealVector(2);
      paramValidator1 = point -> {
        // Ensure diffusion coefficient and precision are positive
        final double d = point.getEntry(0);
        final double sigma = point.getEntry(1);
        // Do not use MIN_VALUE here to avoid sub-normal numbers
        if (d < MIN_D) {
          point.setEntry(0, MIN_D);
        }
        if (sigma < MIN_SIGMA) {
          point.setEntry(1, MIN_SIGMA);
        }
        return point;
      };

      // Non-linear model for anomalous diffusion coefficient
      model2 = new FbmDiffusionFunction(wm1, deltaT);
      start2 = new ArrayRealVector(3);
      final double minAlpha = Math.max(settings.minAlpha, MIN_ALPHA);
      final double maxAlpha = settings.maxAlpha;
      paramValidator2 = point -> {
        // Ensure diffusion coefficient and precision are positive.
        paramValidator1.validate(point);
        // Ensure alpha in the specified range. Default is (0, 2].
        final double alpha = point.getEntry(2);
        // Since the computations require (alpha + 1) limit this to the ULP of 1.0
        // so the parameter always has an effect.
        if (alpha < minAlpha) {
          point.setEntry(2, MIN_ALPHA);
        } else if (alpha > maxAlpha) {
          point.setEntry(2, maxAlpha);
        }
        return point;
      };

      // For linear fit estimation and simple alpha coefficient fit estimation
      reg = new SimpleRegression(true);
      reg2 = new SimpleRegression(true);
      // For significance test of the least squares fit.
      // numeratorDegreesOfFreedom = numberOfParameters2 - numberOfParameters1
      // denominatorDegreesOfFreedom = numberOfPoints - numberOfParameters2
      denominatorDegreesOfFreedom = (int) MathUtils.sum(weight) - 3;
      distribution = FDistribution.of(1, denominatorDegreesOfFreedom);
      useT1Fraction = settings.useT1Fraction;
    }

    /**
     * Sets the data.
     *
     * @param x the x
     * @param y the y
     */
    void setData(double[] x, double[] y) {
      this.x = x;
      this.y = y;
    }

    /**
     * Fit the provided mean squared jump distances. This ignores any data provided by the
     * {@link #setData(double[], double[])} method as the MSDs are already computed. The input array
     * must match the expected length defined by the window size minus 1.
     *
     * <p>Fitting exceptions are thrown.
     *
     * @param msds the msds
     */
    void fit(SumOfSquaredDeviations[] msds) {
      ValidationUtils.checkArgument(msds.length == s.length);
      // Initialise the regression used to estimate the slope.
      // This is not a true weighted regression but is only used as a start point so this is OK.
      reg.clear();
      reg2.clear();
      final double[] weights = new double[msds.length];
      double ssy = 0;
      for (int i = 0; i < msds.length; i++) {
        final double t = (i + 1) * deltaT;
        s[i] = msds[i].getMean();
        reg.addData(t, s[i]);
        final double logT = Math.log(t);
        reg2.addData(logT, Math.log(s[i]));
        weights[i] = msds[i].getN();
        ssy += msds[i].getSumOfSquaredDeviations();
      }
      if (useT1Fraction) {
        reg2.addData(Math.log(deltaT) + LOG_T1_FRACTION, Math.log(s[0] * T1_FRACTION));
      }
      // Set the weights and degrees of freedom
      final RealMatrix weightMatrix = new DiagonalMatrix(weights, false);
      final int denominatorDegreesOfFreedom = (int) MathUtils.sum(weights) - 3;
      final FDistribution distribution = FDistribution.of(1, denominatorDegreesOfFreedom);
      fit(weightMatrix, ssy, denominatorDegreesOfFreedom, distribution);
    }

    /**
     * Fit the mean squared jump distances computed from the data points {@code k} to
     * {@code k + window - 1}. The window was defined in the settings passed to the constructor.
     *
     * <p>Fitting exceptions are thrown.
     *
     * @param k the k
     * @param msds the msds (can be null)
     * @param regressionOnly if true only perform the regression fit
     */
    void fit(int k, SumOfSquaredDeviations[] msds, boolean regressionOnly) {
      // First point in window = k
      // Last point in window = k + w - 1 = k + wm1
      final int end = k + wm1;

      // 1. Anomalous exponent.

      // Compute the MSD for all distances from m=0, to m=window-1
      // For the MSD fit compute the gradient and intercept using a linear regression.
      // (This could exploit pre-computation of the regression x components.)

      // Do a regression of log(MSD) vs log(t):
      // MSD = A t^alpha
      // log(MSD) = log(A) + alpha * log(t)

      reg.clear();
      reg2.clear();
      double ssy = 0;
      for (int m = 1; m <= wm1; m++) {
        // Use intermediate points to compute an average
        final SumOfSquaredDeviations msd = new SumOfSquaredDeviations();
        final double t = m * deltaT;
        for (int i = end - m; i >= k; i--) {
          final double d = MathUtils.distance2(x[i], y[i], x[i + m], y[i + m]);
          msd.add(d);
          reg.addData(t, d);
        }
        s[m - 1] = msd.getMean();
        final double logT = Math.log(t);
        reg2.addData(logT, Math.log(msd.getMean()));
        ssy += msd.getSumOfSquaredDeviations();
        if (msds != null) {
          msds[m - 1] = msd;
        }
      }

      if (useT1Fraction) {
        reg2.addData(Math.log(deltaT) + LOG_T1_FRACTION, Math.log(s[0] * T1_FRACTION));
      }

      // System.out.println(Arrays.toString(s));
      // System.out.println(reg2.getSlope());
      // System.out.println(Math.exp(reg2.getIntercept()));

      if (!regressionOnly) {
        fit(weightMatrix, ssy, denominatorDegreesOfFreedom, distribution);
      }
    }

    /**
     * Fit the prepared mean squared jump distances.
     *
     * <p>The total sum of squared deviations of the raw jump distances from the mean is required to
     * compute the correct residual sum of squares for the significance test. The weighted fit (with
     * weights equal to the count of data points for each mean) is therefore a least squared fit of
     * all jump distances.
     *
     * <p>Fitting exceptions are thrown.
     *
     * @param weightMatrix the weight matrix
     * @param ssy the sum-of-squared deviations of all data points y (jump distances) from the mean
     *        squared distance (msd) for each jump interval
     * @param denominatorDegreesOfFreedom the denominator degrees of freedom
     * @param distribution the distribution
     */
    private void fit(RealMatrix weightMatrix, double ssy, double denominatorDegreesOfFreedom,
        FDistribution distribution) {
      final int maxEvaluations = Integer.MAX_VALUE;
      final int maxIterations = 3000;
      final boolean lazyEvaluation = false;

      // Note:
      // If the alpha cannot significantly improve the fit then its value is likely to be
      // bad and the feature is noise. This can be judged using an F-test.
      // The local MSD curve can be highly variable, especially
      // if the molecule has switched from bound to unbound or vice versa. In this case there
      // is no good value for alpha.
      // Setting alpha to 1.0 for insignificant fits can result in no variance in the feature
      // and failure to invert the covariance matrix during Expectation-Maximisation.
      // Thus the user has the option to remove alpha from the features
      // before fitting the Gaussian mixture. The value is still computed for the results.

      // Compute linear fit with the (n-1/3) correction factor:
      // MSD = 4D n - (4D) / 3 + 4 s^2
      final double slope = reg.getSlope();
      positiveSlope = slope > 4 * MIN_D;
      if (positiveSlope) {
        start1.setEntry(0, slope / 4);
        start1.setEntry(1, Math.sqrt(Math.max(0, reg.getIntercept() + slope / 3) / 4));
      } else {
        // Do not allow negative slope. Use a flat line and
        // precision is derived from the mean of the MSD.
        // Using D=0 has no gradient and the fit does not explore the parameter space.
        // Note: If there is no gradient then the alpha value is very likely to
        // be insignificant as the FBM model assumes an upward slope.
        start1.setEntry(0, MIN_D);
        start1.setEntry(1, Math.sqrt((MathUtils.sum(s) / wm1) / 4));
      }
      final LeastSquaresProblem problem1 = LeastSquaresFactory.create(model1, observed, start1,
          weightMatrix, checker, maxEvaluations, maxIterations, lazyEvaluation, paramValidator1);
      lvmSolution1 = optimizer.optimize(problem1);

      // Fit the anomalous exponent alpha. Start from the linear fit result.
      final RealVector fit1 = lvmSolution1.getPoint();
      start2.setEntry(0, fit1.getEntry(0));
      start2.setEntry(1, fit1.getEntry(1));
      start2.setEntry(2, 1.0);
      final LeastSquaresProblem problem2 = LeastSquaresFactory.create(model2, observed, start2,
          weightMatrix, checker, maxEvaluations, maxIterations, lazyEvaluation, paramValidator2);
      lvmSolution2 = optimizer.optimize(problem2);

      // Check for model improvement
      final double rss1 = getResidualSumOfSquares(lvmSolution1);
      final double rss2 = getResidualSumOfSquares(lvmSolution2);

      // Compare significance using the method from RegressionUtils
      // F = ((rss1 - rss2) / (p2 - p1)) / (rss2 / (n - p2))
      final double num = rss1 - rss2;

      // Adjust RSS from the weighted model to the equivalent for a full regression
      // on all the data (not the means).
      // For a collection yi sum of squares to a point x:
      // sum (yi - x)^2 = sum (yi - u)^2 + sum (x - u)^2 * n
      // x = the point
      // u = mean of yi
      // n = number of yi
      // Currently RSS is sum ( (x - u)^2 * n ) where:
      // x is the model value
      // u is the MSD
      // n is the number of MSD observations
      // Thus we adjust using the sum (yi - u)^2 for each MSD.

      final double denom = (ssy + rss2) / denominatorDegreesOfFreedom;
      // Note: If ss2 = 0 then f-statistic is +infinity and the cumulative probability is NaN
      // and the p-value will be accepted
      final double fStatistic = MathUtils.div0(num, denom);
      pValue = 1.0 - distribution.cumulativeProbability(fStatistic);

      // // Debug
      // if (
      // pValue < 0.2
      // // alpha > 0.0 && alpha < 0.2
      // //slope < 0
      // ) {
      // final RealVector p = lvmSolution2.getPoint();
      // final String title = "anomalous exponent";
      // final Plot plot = new Plot(title, "time (s)", "MSD (um^^2^^)");
      // final double[] t = SimpleArrayUtils.newArray(s.length, deltaT, deltaT);
      // plot.addLabel(0, 0, lvmSolution2.getPoint().toString() + " p=" + pValue + ". "
      // + lvmSolution1.getPoint().toString());
      // plot.addPoints(t, s, Plot.CROSS);
      // plot.addPoints(t, model2.value(p).getFirst().toArray(), Plot.LINE);
      // plot.setColor(Color.BLUE);
      // plot.addPoints(t, model1.value(lvmSolution1.getPoint()).getFirst().toArray(),
      // Plot.LINE);
      // plot.setColor(Color.RED);
      // final double[] yy = Arrays.stream(t).map(reg::predict).toArray();
      // plot.addPoints(t, yy, Plot.CIRCLE);
      // ImageJUtils.display(title, plot, ImageJUtils.NO_TO_FRONT);
      // System.out.println(lvmSolution2.getPoint());
      // }
    }

    /**
     * Gets the alpha estimated from a linear regression fit of log(MSD) = log(A) + alpha * log(t).
     *
     * @return the alpha
     */
    double getRegressionAlpha() {
      return reg2.getSlope();
    }

    /**
     * Gets the diffusion coefficient D estimated from a linear regression fit of log(MSD) = log(A)
     * + alpha * log(t).
     *
     * @return the diffusion coefficient D
     */
    double getRegressionD() {
      // log(MSD) = log(A) + alpha * log(t)
      // A = 4D
      return Math.exp(reg2.getIntercept()) * 0.25f;
    }
  }

  /**
   * Define a Brownian motion diffusion function for use in fitting MSD for the case of pure
   * Brownian motion. The function is corrected for static (localisation precision) and dynamic
   * error (localisation error due to diffusion):
   *
   * <pre>
   * MSD = 4Dt * (n - 1/3) + 4s^2
   *     = 4Dtn - 4Dt/3 + 4s^2
   * </pre>
   *
   * <p>D is the diffusion coefficient; t is the acquisition time step; n is the number of frames;
   * and s is the static localisation precision (see Backlund, et al (2015), Physical Review E 91:
   * 062716, eq. 1). Note s is composed of the localisation precision of an unmoving particle (s0)
   * and a contribution due to the effect of a moving particle (see eq. 5).
   *
   * <pre>
   * s^2 = s0^2 + 2Dt / 3p
   * </pre>
   *
   * <p>where p is the number of photons. This term typically results in a small inflation of s and
   * is ignored in fitting.
   *
   * <p>This function assumes {@code n} is an integer in {@code [1, size]}.
   */
  @VisibleForTesting
  static final class BrownianDiffusionFunction implements MultivariateJacobianFunction {
    private static final double THIRD = 1 / 3.0;
    private final int size;
    /** The pre-compute 4t * (n - 1/3). */
    private final double[] scale;

    /**
     * Create an instance.
     *
     * @param size the maximum size of n (inclusive)
     * @param deltaT the time step (delta T)
     */
    BrownianDiffusionFunction(int size, double deltaT) {
      this.size = size;
      // Pre-computation of the scale for n
      scale = new double[size];
      for (int n = 1; n <= size; n++) {
        scale[n - 1] = 4 * deltaT * (n - THIRD);
      }
    }

    @Override
    public Pair<RealVector, RealMatrix> value(RealVector point) {
      final double[] value = new double[size];
      final double[][] jacobian = new double[size][2];
      final double d = point.getEntry(0);
      final double s = point.getEntry(1);
      // MSD = 4Dt * (n - 1/3) + 4s^2
      // dy_dD = 4t * (n - 1/3)
      // dy_ds = 8s
      final double ss4 = 4 * s * s;
      final double s8 = 8 * s;
      // Here the variable n represents (n+1) but is used as an index to a pre-computed scale
      for (int n = 0; n < size; n++) {
        value[n] = d * scale[n] + ss4;
        jacobian[n][0] = scale[n];
        jacobian[n][1] = s8;
      }
      return new Pair<>(new ArrayRealVector(value, false),
          new Array2DRowRealMatrix(jacobian, false));
    }
  }

  /**
   * Define a fractional Brownian motion (FBM) diffusion function for use in fitting MSD for the
   * case of anomalous diffusion. The function is corrected for static (localisation precision) and
   * dynamic error (localisation error due to diffusion):
   *
   * <pre>
   * MSD = [4Dt^a / (a+2)(a+1)] * [(n+1)^(a+2) + (n-1)^(a+2) - 2n^(a+2)]
   *       - [8Dt^a / (a+2)(a+1)] + 4s^2
   * </pre>
   *
   * <p>D is the diffusion coefficient; t is the acquisition time step; n is the number of frames; a
   * is the anomalous coefficient alpha; and s is the static localisation precision (see Backlund,
   * et al (2015), Physical Review E 91: 062716, eq. 8). Note s is composed of the localisation
   * precision of an unmoving particle (s0) and a contribution due to the effect of a moving
   * particle (see eq. 6).
   *
   * <pre>
   * s^2 = s0^2 + 4Dt / (a+2)(a+1)p
   * </pre>
   *
   * <p>where p is the number of photons. This term typically results in a small inflation of s and
   * is ignored in fitting. The factor alpha should be in the range {@code (0, 2]}.
   *
   * <p>This function assumes {@code n} is an integer in {@code [1, size]}.
   */
  @VisibleForTesting
  static final class FbmDiffusionFunction implements MultivariateJacobianFunction {
    private final int size;
    private final double[] log;
    private final double deltaT;
    private final double logt;

    /**
     * Create an instance.
     *
     * @param size the maximum size of n (inclusive)
     * @param deltaT the time step (delta T)
     */
    FbmDiffusionFunction(int size, double deltaT) {
      this.size = size;
      log = new double[size + 2];
      for (int n = 2; n < log.length; n++) {
        log[n] = Math.log(n);
      }
      this.deltaT = deltaT;
      logt = Math.log(deltaT);
    }

    @Override
    public Pair<RealVector, RealMatrix> value(RealVector point) {
      final double[] value = new double[size];
      final double[][] jacobian = new double[size][3];
      final double d = point.getEntry(0);
      final double s = point.getEntry(1);
      final double a = point.getEntry(2);
      // MSD = [4Dt^a / (a+2)(a+1)] * [(n+1)^(a+2) + (n-1)^(a+2) - 2n^(a+2)]
      // - [8Dt^a / (a+2)(a+1)] + 4s^2
      // dy_dD = [4t^a / (a+2)(a+1)] * [(n+1)^(a+2) + (n-1)^(a+2) - 2n^(a+2)] - 8t^a / (a+2)(a+1)
      // dy_ds = 8s
      // Use product rule with the following parts:
      // e(a) = t^a
      // e'(a) = t^a * log(t)
      // f(a) = 4D / (a+2)(a+1)
      // f'(a) = -4D / (a+2)(a+1)(a+2) -4 / (a+2)(a+1)(a+1)
      // g(a) = (n+1)^(a+2) + (n-1)^(a+2) - 2n^(a+2)
      // g'(a) = (n+1)^(a+2)*log(n+1) + (n-1)^(a+2)*log(n-1) - 2n^(a+2)*log(n)
      // y = e(a) * f(a) * g(a) - 2 * e(a) * f(a)
      // dy_da = e(a) * (f'(a) * g(a) + f * g'(a)) + e'(a) * (f(a) * g(a))
      // - 2 * (e(a) * f'(a) + e'(a) + f(a))
      final double ss4 = 4 * s * s;
      final double s8 = 8 * s;
      final double a1 = a + 1;
      final double a2 = a + 2;
      final double a2a1 = a2 * a1;
      final double ea = Math.pow(deltaT, a);
      final double epa = ea * logt;
      // Missing the factor D which is multiplied in later
      final double fa = 4 / a2a1;

      // Pre-compute powers of n^(a+2)
      final double[] na = new double[size + 2];
      for (int n = 1; n < na.length; n++) {
        na[n] = Math.pow(n, a2);
      }

      for (int n = 1; n <= size; n++) {
        final int i = n - 1;
        final double ga = na[n + 1] + na[i] - 2 * na[n];
        jacobian[i][0] = ea * fa * ga - 2 * ea * fa;
        value[i] = d * jacobian[i][0] + ss4;
        jacobian[i][1] = s8;
        final double fpa = -fa / a2 - fa / a1;
        final double gpa = (na[n + 1] * log[n + 1] + na[i] * log[i] - 2 * na[n] * log[n]);
        jacobian[i][2] =
            d * ((ea * (fpa * ga + fa * gpa) + epa * fa * ga) - 2 * (ea * fpa + epa * fa));
      }
      return new Pair<>(new ArrayRealVector(value, false),
          new Array2DRowRealMatrix(jacobian, false));
    }
  }
}
