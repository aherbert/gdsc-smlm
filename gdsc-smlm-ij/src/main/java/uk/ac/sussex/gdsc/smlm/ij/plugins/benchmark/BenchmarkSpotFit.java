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

package uk.ac.sussex.gdsc.smlm.ij.plugins.benchmark;

import com.thoughtworks.xstream.XStreamException;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.plugin.PlugIn;
import ij.text.TextWindow;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import java.awt.Checkbox;
import java.awt.Color;
import java.awt.Rectangle;
import java.awt.TextArea;
import java.awt.TextField;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Function;
import java.util.function.IntConsumer;
import org.apache.commons.lang3.concurrent.ConcurrentRuntimeException;
import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot.HistogramPlotBuilder;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.logging.Ticker;
import uk.ac.sussex.gdsc.core.match.Assignment;
import uk.ac.sussex.gdsc.core.match.AssignmentComparator;
import uk.ac.sussex.gdsc.core.match.BasePoint;
import uk.ac.sussex.gdsc.core.match.Coordinate;
import uk.ac.sussex.gdsc.core.match.FractionClassificationResult;
import uk.ac.sussex.gdsc.core.match.FractionalAssignment;
import uk.ac.sussex.gdsc.core.match.ImmutableFractionalAssignment;
import uk.ac.sussex.gdsc.core.match.PointPair;
import uk.ac.sussex.gdsc.core.utils.Correlator;
import uk.ac.sussex.gdsc.core.utils.FastCorrelator;
import uk.ac.sussex.gdsc.core.utils.FileUtils;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.OpenHashMaps.CustomInt2ObjectOpenHashMap;
import uk.ac.sussex.gdsc.core.utils.RampedScore;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.SortUtils;
import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.XmlUtils;
import uk.ac.sussex.gdsc.core.utils.concurrent.ConcurrencyUtils;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FitSolver;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.NoiseEstimatorMethod;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.PrecisionMethod;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtosHelper;
import uk.ac.sussex.gdsc.smlm.engine.FitConfiguration;
import uk.ac.sussex.gdsc.smlm.engine.FitEngineConfiguration;
import uk.ac.sussex.gdsc.smlm.engine.FitParameters;
import uk.ac.sussex.gdsc.smlm.engine.FitParameters.FitTask;
import uk.ac.sussex.gdsc.smlm.engine.FitWorker;
import uk.ac.sussex.gdsc.smlm.engine.ParameterisedFitJob;
import uk.ac.sussex.gdsc.smlm.filters.MaximaSpotFilter;
import uk.ac.sussex.gdsc.smlm.filters.Spot;
import uk.ac.sussex.gdsc.smlm.fitting.FitResult;
import uk.ac.sussex.gdsc.smlm.fitting.FitStatus;
import uk.ac.sussex.gdsc.smlm.ij.IJImageSource;
import uk.ac.sussex.gdsc.smlm.ij.plugins.HelpUrls;
import uk.ac.sussex.gdsc.smlm.ij.plugins.PeakFit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.PeakFit.FitEngineConfigurationProvider;
import uk.ac.sussex.gdsc.smlm.ij.plugins.PsfCalculator;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsMatchCalculator;
import uk.ac.sussex.gdsc.smlm.ij.plugins.SmlmUsageTracker;
import uk.ac.sussex.gdsc.smlm.ij.plugins.benchmark.BenchmarkSpotFilter.BenchmarkSpotFilterResult;
import uk.ac.sussex.gdsc.smlm.ij.plugins.benchmark.BenchmarkSpotFilter.FilterResult;
import uk.ac.sussex.gdsc.smlm.ij.plugins.benchmark.BenchmarkSpotFilter.ScoredSpot;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.GUIFilterSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.ij.utils.ImageJImageConverter;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResultPoint;
import uk.ac.sussex.gdsc.smlm.results.PeakResults;
import uk.ac.sussex.gdsc.smlm.results.SynchronizedPeakResults;
import uk.ac.sussex.gdsc.smlm.results.count.NullFailCounter;
import uk.ac.sussex.gdsc.smlm.results.filter.BasePreprocessedPeakResult;
import uk.ac.sussex.gdsc.smlm.results.filter.CoordinateStoreFactory;
import uk.ac.sussex.gdsc.smlm.results.filter.DirectFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.EShiftFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.Filter;
import uk.ac.sussex.gdsc.smlm.results.filter.FilterSet;
import uk.ac.sussex.gdsc.smlm.results.filter.FilterValidationFlag;
import uk.ac.sussex.gdsc.smlm.results.filter.FilterXStreamUtils;
import uk.ac.sussex.gdsc.smlm.results.filter.MultiFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.MultiFilter2;
import uk.ac.sussex.gdsc.smlm.results.filter.MultiFilterCrlb;
import uk.ac.sussex.gdsc.smlm.results.filter.MultiPathFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.MultiPathFilter.FractionScoreStore;
import uk.ac.sussex.gdsc.smlm.results.filter.MultiPathFitResult;
import uk.ac.sussex.gdsc.smlm.results.filter.MultiPathFitResults;
import uk.ac.sussex.gdsc.smlm.results.filter.ParameterType;
import uk.ac.sussex.gdsc.smlm.results.filter.PeakFractionalAssignment;
import uk.ac.sussex.gdsc.smlm.results.filter.PrecisionFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.PreprocessedPeakResult;
import uk.ac.sussex.gdsc.smlm.results.filter.ResultAssignment;
import uk.ac.sussex.gdsc.smlm.results.filter.ShiftFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.SignalFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.SnrFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.WidthFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.WidthFilter2;

/**
 * Fits all the candidate spots identified by the benchmark spot filter plugin.
 */
public class BenchmarkSpotFit implements PlugIn, ItemListener {
  // TODO - Add support for optimising z-depth during 3D fitting

  /** The title. */
  static final String TITLE = "Fit Spot Data";

  private static FilterCriteria[] filterCriteria;
  private static final int FILTER_SIGNAL = 0;
  private static final int FILTER_SNR = 1;
  private static final int FILTER_MIN_WIDTH = 2;
  private static final int FILTER_MAX_WIDTH = 3;
  private static final int FILTER_SHIFT = 4;
  private static final int FILTER_ESHIFT = 5;
  private static final int FILTER_PRECISION = 6;
  private static final int FILTER_ITERATIONS = 7;
  private static final int FILTER_EVALUATIONS = 8;

  /** The coordinate cache. This stores the coordinates for a simulation Id. */
  private static AtomicReference<
      Pair<Integer, Int2ObjectOpenHashMap<List<Coordinate>>>> coordinateCache =
          new AtomicReference<>(Pair.of(-1, null));

  private static AtomicReference<TextWindow> summaryTableRef = new AtomicReference<>();

  /**
   * The prefix for the results table header. Contains all the standard header data about the input
   * results data. This is effectively immutable.
   */
  private static String tablePrefix;

  private static AtomicReference<CandidateData> candidateDataCache = new AtomicReference<>();

  /** A reference to the most recent results. */
  private static AtomicReference<BenchmarkSpotFitResult> spotFitResults = new AtomicReference<>();

  private final FitEngineConfiguration config;
  private MultiPathFilter multiFilter;

  private TextArea taFilterXml;
  private TextField textFailLimit;
  private Checkbox cbIncludeNeighbours;
  private TextField textNeighbourHeight;
  private Checkbox cbComputeDoublets;

  private boolean extraOptions;
  /** Flag used when being called by another plugin to suppress dialogs. */
  private boolean silent;
  /** Flag used when being called by another plugin to indicate success. */
  private boolean finished;

  private ImagePlus imp;
  private MemoryPeakResults results;
  private CreateData.SimulationParameters simulationParameters;
  private BenchmarkSpotFilterResult filterResult;

  private MaximaSpotFilter spotFilter;

  /** The distance in pixels. */
  private double distanceInPixels;

  /** The lower distance in pixels. */
  private double lowerDistanceInPixels;


  /** The plugin settings. */
  private final Settings settings;

  // Used to try and guess the range for filtering the results
  private enum LowerLimit {
    ZERO(false), MIN(false), ONE_PERCENT(false), MAX_NEGATIVE_CUMUL_DELTA(
        true), HALF_MAX_JACCARD_VALUE(false, true);

    final boolean requiresDelta;
    final boolean requiresJaccard;

    LowerLimit(boolean requiresDelta) {
      this(requiresDelta, false);
    }

    LowerLimit(boolean requiresDelta, boolean requiresJaccard) {
      this.requiresDelta = requiresDelta;
      this.requiresJaccard = requiresJaccard;
    }

    public boolean requiresDeltaHistogram() {
      return requiresDelta;
    }
  }

  private enum UpperLimit {
    ZERO(false), MAX_POSITIVE_CUMUL_DELTA(true), NINETY_NINE_PERCENT(
        false), NINETY_NINE_NINE_PERCENT(false), MAX_JACCARD2(false, true);

    final boolean requiresDelta;
    final boolean requiresJaccard;

    UpperLimit(boolean requiresDelta) {
      this(requiresDelta, false);
    }

    UpperLimit(boolean requiresDelta, boolean requiresJaccard) {
      this.requiresDelta = requiresDelta;
      this.requiresJaccard = requiresJaccard;
    }

    public boolean requiresDeltaHistogram() {
      return requiresDelta;
    }
  }

  private static class FilterCriteria {
    final ParameterType type;
    final String name;
    final LowerLimit lower;
    final UpperLimit upper;
    final int minBinWidth;
    final boolean restrictRange;
    final boolean requireLabel;

    FilterCriteria(ParameterType type, LowerLimit lower, UpperLimit upper) {
      this(type, type.toString(), lower, upper, 0, true, true);
    }

    FilterCriteria(ParameterType type, String name, LowerLimit lower, UpperLimit upper,
        int minBinWidth, boolean restrictRange, boolean requireLabel) {
      this.type = type;
      this.name = name;
      this.lower = lower;
      this.upper = upper;
      this.minBinWidth = minBinWidth;
      this.restrictRange = restrictRange;
      this.requireLabel = requireLabel;
    }
  }

  /**
   * Store the filter candidates data.
   */
  private static class CandidateData {
    final Int2ObjectOpenHashMap<FilterCandidates> filterCandidates;
    final double fractionPositive;
    final double fractionNegative;
    final int countPositive;
    final int countNegative;

    final int filterId;
    final double fractionPositives;
    final double fractionNegativesAfterAllPositives;
    final int negativesAfterAllPositives;
    final int width;

    CandidateData(Int2ObjectOpenHashMap<FilterCandidates> filterCandidates, int filterId,
        double fractionPositive, double fractionNegative, int countPositive, int countNegative,
        Settings settings, int width) {
      this.filterCandidates = filterCandidates;
      this.filterId = filterId;
      this.fractionPositive = fractionPositive;
      this.fractionNegative = fractionNegative;
      this.countPositive = countPositive;
      this.countNegative = countNegative;
      fractionPositives = settings.fractionPositives;
      fractionNegativesAfterAllPositives = settings.fractionNegativesAfterAllPositives;
      negativesAfterAllPositives = settings.negativesAfterAllPositives;
      this.width = width;
    }

    /**
     * Return true if the settings used in the candidate data are different.
     *
     * @param filterId the filter id
     * @param settings the settings
     * @param width the width
     * @return true if different
     */
    boolean differentSettings(int filterId, Settings settings, int width) {
      return this.filterId != filterId || fractionPositives != settings.fractionPositives
          || fractionNegativesAfterAllPositives != settings.fractionNegativesAfterAllPositives
          || negativesAfterAllPositives != settings.negativesAfterAllPositives
          || this.width != width;
    }
  }

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The default MultiPathFilter. */
    static final MultiPathFilter defaultMultiFilter;
    /** The default parameters. */
    static final double[] defaultParameters;

    static {
      // Add a filter to use for storing the slice results:
      // Use the standard configuration to ensure sensible fits are stored as the current slice
      // results.
      final FitConfiguration tmp = FitConfiguration.create();
      final PrecisionMethod precisionMethod = PrecisionMethod.MORTENSEN_LOCAL_BACKGROUND;
      tmp.setPrecisionMethod(precisionMethod);

      final DirectFilter primaryFilter = tmp.getDefaultSmartFilter();

      // We might as well use the doublet fits given we will compute them.
      final double residualsThreshold = 0.4;
      defaultMultiFilter = new MultiPathFilter(primaryFilter,
          FitWorker.createMinimalFilter(precisionMethod), residualsThreshold);
      defaultParameters = createParameters(createDefaultConfig());
    }

    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    FitEngineConfiguration config;
    MultiPathFilter multiFilter;

    double fractionPositives;
    double fractionNegativesAfterAllPositives;
    int negativesAfterAllPositives;
    double distance;
    double lowerDistance;
    double signalFactor;
    double lowerSignalFactor;
    boolean computeDoublets;
    boolean showFilterScoreHistograms;
    boolean saveFilterRange;
    boolean showCorrelation;
    boolean rankByIntensity;

    Settings() {
      // Set defaults
      config = createDefaultConfig();
      multiFilter = defaultMultiFilter.copy();
      fractionPositives = 100;
      fractionNegativesAfterAllPositives = 50;
      negativesAfterAllPositives = 10;
      distance = 1.5;
      lowerDistance = 1.5;
      signalFactor = 2;
      lowerSignalFactor = 1;
      computeDoublets = true;
      saveFilterRange = true;
    }

    Settings(Settings source) {
      config = source.config.createCopy();
      multiFilter = source.multiFilter.copy();
      fractionPositives = source.fractionPositives;
      fractionNegativesAfterAllPositives = source.fractionNegativesAfterAllPositives;
      negativesAfterAllPositives = source.negativesAfterAllPositives;
      distance = source.distance;
      lowerDistance = source.lowerDistance;
      signalFactor = source.signalFactor;
      lowerSignalFactor = source.lowerSignalFactor;
      computeDoublets = source.computeDoublets;
      showFilterScoreHistograms = source.showFilterScoreHistograms;
      saveFilterRange = source.saveFilterRange;
      showCorrelation = source.showCorrelation;
      rankByIntensity = source.rankByIntensity;
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
    }

    private static FitEngineConfiguration createDefaultConfig() {
      final FitEngineConfiguration config = FitEngineConfiguration.create();
      final FitConfiguration fitConfig = config.getFitConfiguration();
      // Set some default fit settings here ...
      fitConfig.setDisableSimpleFilter(false);
      fitConfig.setMinPhotons(1); // Do not allow negative photons
      fitConfig.setSignalStrength(2); // Lower than the default of 5
      fitConfig.setCoordinateShiftFactor(0);
      fitConfig.setPrecisionThreshold(0);
      fitConfig.setMinWidthFactor(0);
      fitConfig.setMaxWidthFactor(0);

      fitConfig.setBackgroundFitting(true);
      fitConfig.setNoise(0);
      config.setNoiseMethod(NoiseEstimatorMethod.QUICK_RESIDUALS_LEAST_MEAN_OF_SQUARES);

      // Use bounded so that we can fit the neighbours
      config.setIncludeNeighbours(true);
      fitConfig.setFitSolver(FitSolver.LVM_LSE);
      return config;
    }
  }

  //////////////////////
  // TODO
  // Move all static field into a class that can be used as the latest set of results
  // by other plugins.
  // This class will hold a reference to its own results and an atomic reference will
  // contain the latest results.
  // Move all benchmarking code to a plugins package?
  //////////////////////

  /**
   * Contains the results from the latest execution.
   */
  static class BenchmarkSpotFitResult {
    // Used by the Benchmark Spot Fit plugin
    private static AtomicInteger fitResultsId = new AtomicInteger(1);

    /** The simulation id. */
    final int simulationId;

    /** The fit results id. */
    final int id;

    /** The fit results. */
    CustomInt2ObjectOpenHashMap<FilterCandidates> fitResults;

    /** The distance in pixels. */
    double distanceInPixels;

    /** The lower distance in pixels. */
    double lowerDistanceInPixels;

    /** Allow access to the analysis time. */
    StopWatch stopWatch;

    /**
     * The prefix for the results table entry. Contains all the standard data about the input
     * results data. This is constructed using all the settings.
     */
    String resultPrefix;

    /** The minimum for each spot parameter. Use to set limits for spot filters. */
    double[] min;

    /** The maximum for each spot parameter. Use to set limits for spot filters. */
    double[] max;

    /**
     * Instantiates a new benchmark spot fit results.
     *
     * @param simulationId the simulation id
     * @param fitResults the fit results
     */
    BenchmarkSpotFitResult(int simulationId,
        CustomInt2ObjectOpenHashMap<FilterCandidates> fitResults) {
      id = fitResultsId.incrementAndGet();
      this.simulationId = simulationId;
      this.fitResults = fitResults;
    }
  }

  private static class MultiPathPoint extends BasePoint {
    static final int SPOT = -1;
    static final int SINGLE = 0;
    static final int MULTI = 1;
    static final int DOUBLET = 2;
    @SuppressWarnings("unused")
    static final int MULTIDOUBLET = 3;

    final PreprocessedPeakResult result;
    final int id;
    final int type;
    final int index;

    MultiPathPoint(PreprocessedPeakResult result, int id, int type, int index) {
      super(result.getX(), result.getY());
      this.result = result;
      this.id = id;
      this.type = type;
      this.index = index;
    }

    MultiPathPoint(float x, float y, int id, int type, int index) {
      super(x, y);
      this.result = null;
      this.id = id;
      this.type = type;
      this.index = index;
    }

    @Override
    public boolean equals(Object object) {
      // Ignore extra fields
      return super.equals(object);
    }

    @Override
    public int hashCode() {
      // Ignore extra fields
      return super.hashCode();
    }
  }

  /**
   * Store details of spot candidates that match actual spots.
   */
  public abstract static class SpotMatch {
    /**
     * The index for the spot candidate.
     */
    final int index;
    /**
     * The distance to the spot.
     */
    final double distance;
    /**
     * The depth of the actual spot.
     */
    final double zdepth;
    /**
     * The score.
     */
    double score;

    /**
     * Instantiates a new spot match.
     *
     * @param index The index for the spot candidate
     * @param distance The distance to the spot
     * @param zdepth The depth of the actual spot
     */
    public SpotMatch(int index, double distance, double zdepth) {
      this.index = index;
      this.distance = distance;
      this.zdepth = zdepth;
    }

    /**
     * Checks if the spot candidate was successfully fitted.
     *
     * @return True if the spot candidate was successfully fitted.
     */
    public abstract boolean isFitResult();

    /**
     * Return a score for the difference between the fitted and actual signal. Zero is no
     * difference. Negative is the fitted is below the actual. Positive means the fitted is above
     * the actual.
     *
     * @return The factor difference between the successfully fitted signal and the actual signal.
     */
    public abstract double getSignalFactor();

    /**
     * Return a score for the difference between the fitted and actual signal. Zero is no
     * difference. Positive means the fitted is difference from the actual.
     *
     * @return The factor difference between the successfully fitted signal and the actual signal.
     */
    public double getAbsoluteSignalFactor() {
      return Math.abs(getSignalFactor());
    }

    /**
     * Gets the score.
     *
     * @return the score
     */
    public double getScore() {
      return score;
    }
  }

  /**
   * Store details of a fitted spot candidate that matches an actual spot.
   */
  public static class FitMatch extends SpotMatch {
    /** The point. */
    final MultiPathPoint point;

    /** The predicted signal. */
    final double predictedSignal;

    /** The actual signal. */
    final double actualSignal;

    /** The signal factor. */
    final double sf;

    /**
     * Instantiates a new fit match.
     *
     * @param point the point
     * @param distance The distance to the spot
     * @param zdepth The depth of the actual spot
     * @param predictedSignal the predicted signal
     * @param actualSignal the actual signal
     */
    public FitMatch(MultiPathPoint point, double distance, double zdepth, double predictedSignal,
        double actualSignal) {
      super(point.index, distance, zdepth);
      this.point = point;
      this.predictedSignal = predictedSignal;
      this.actualSignal = actualSignal;
      this.sf = BenchmarkSpotFit.computeSignalFactor(predictedSignal, actualSignal);
    }

    @Override
    public boolean isFitResult() {
      return true;
    }

    @Override
    public double getSignalFactor() {
      return sf;
    }
  }

  /**
   * Computes the signal factor.
   *
   * @param predictedSignal the predicted signal
   * @param actualSignal the actual signal
   * @return the signal factor
   */
  static double computeSignalFactor(double predictedSignal, double actualSignal) {
    final double rsf = predictedSignal / actualSignal;
    // The relative signal factor is 1 for a perfect fit, less than 1 for below and above 1 for
    // above. Reset the signal factor from 0.
    // abs(log(s/t)) / log(2): 1x=0, 1.5x= 0.58, 2x=1, 3x=1.58, 4x=2
    return Math.abs(MathUtils.log2(rsf));
  }

  /**
   * Store details of a spot candidate that matches an actual spot.
   */
  public static class CandidateMatch extends SpotMatch {

    /**
     * Instantiates a new candidate match.
     *
     * @param index The index for the spot candidate
     * @param distance The distance to the spot
     * @param z The depth of the actual spot
     */
    public CandidateMatch(int index, double distance, double z) {
      super(index, distance, z);
    }

    @Override
    public boolean isFitResult() {
      return false;
    }

    @Override
    public double getSignalFactor() {
      return 0;
    }
  }

  /**
   * Store the filter candidates.
   */
  static class FilterCandidates implements Cloneable {
    /** Integer counts of positives (matches). */
    final int pos;
    /** Integer counts of negatives. */
    final int neg;
    /** Double sums of the fractions match score. */
    final double np;
    /** Double sums of the fractions match antiscore. */
    final double nn;

    /** The spots. */
    final ScoredSpot[] spots;

    /** The max candidate. */
    final int maxCandidate;
    /** True positives. */
    double tp;
    /** False positives. */
    double fp;
    /** True negatives. */
    double tn;
    /** False negatives. */
    double fn;

    /** The fit result. */
    MultiPathFitResult[] fitResult;

    /** The noise. */
    float noise;

    /**
     * Store the z-position of the actual spots for later analysis. Size is the number of actual
     * spots
     */
    double[] zPosition;

    /**
     * Store details about the actual spots that were matched by spot candidates or fitted spot
     * candidates.
     */
    SpotMatch[] match;

    /**
     * Instantiates a new filter candidates.
     *
     * @param pos the positives
     * @param neg the negatives
     * @param np the sum of the fraction match score
     * @param nn the sum of the fraction match antiscore
     * @param spots the spots
     * @param maxCandidate the max candidate
     */
    FilterCandidates(int pos, int neg, double np, double nn, ScoredSpot[] spots, int maxCandidate) {
      this.pos = pos;
      this.neg = neg;
      this.np = np;
      this.nn = nn;
      this.spots = spots;
      this.maxCandidate = maxCandidate;
    }

    @Override
    public FilterCandidates clone() {
      try {
        return (FilterCandidates) super.clone();
      } catch (final CloneNotSupportedException ex) {
        return null;
      }
    }
  }

  private static class Ranking {
    final double value;
    final int index;

    Ranking(double value, int index) {
      this.value = value;
      this.index = index;
    }
  }

  /**
   * Used to allow multi-threading of the fitting method.
   */
  private class Worker implements Runnable {
    volatile boolean finished;
    final BlockingQueue<Integer> jobs;
    final ImageStack stack;
    final FitWorker fitWorker;
    final Int2ObjectOpenHashMap<List<Coordinate>> actualCoordinates;
    final Int2ObjectOpenHashMap<FilterCandidates> filterCandidates;
    final Int2ObjectOpenHashMap<FilterCandidates> results;
    final Rectangle bounds;
    final MultiPathFilter multiFilter;
    final Ticker ticker;

    float[] data;
    List<PointPair> matches = new ArrayList<>();

    public Worker(BlockingQueue<Integer> jobs, ImageStack stack,
        Int2ObjectOpenHashMap<List<Coordinate>> actualCoordinates,
        Int2ObjectOpenHashMap<FilterCandidates> filterCandidates, PeakResults peakResults,
        Ticker ticker) {
      this.jobs = jobs;
      this.stack = stack;
      this.fitWorker = new FitWorker(config.createCopy(), peakResults, null);

      final int fitting = config.getFittingWidth();
      fitWorker.setSearchParameters((MaximaSpotFilter) spotFilter.copy(), fitting);

      this.actualCoordinates = actualCoordinates;
      this.filterCandidates = filterCandidates;
      this.results = new Int2ObjectOpenHashMap<>();
      bounds = new Rectangle(0, 0, stack.getWidth(), stack.getHeight());
      // Instance copy
      multiFilter = BenchmarkSpotFit.this.multiFilter.copy();
      this.ticker = ticker;
    }

    @Override
    public void run() {
      try {
        for (;;) {
          final Integer job = jobs.take();
          if (job == null || job.intValue() < 0) {
            break;
          }
          if (!finished) {
            // Only run jobs when not finished. This allows the queue to be emptied.
            run(job.intValue());
            ticker.tick();
          }
        }
      } catch (final InterruptedException ex) {
        ConcurrencyUtils.interruptAndThrowUncheckedIf(!finished, ex);
      } finally {
        finished = true;
      }
    }

    private void run(int frame) {
      if (ImageJUtils.isInterrupted()) {
        finished = true;
        return;
      }

      // Extract the data
      data = ImageJImageConverter.getData(stack.getPixels(frame), stack.getWidth(),
          stack.getHeight(), null, data);

      FilterCandidates candidates = filterCandidates.get(frame);
      final int totalCandidates = candidates.spots.length;
      final MultiPathFitResult[] fitResult = new MultiPathFitResult[totalCandidates];

      // Fit the candidates and store the results
      final FitParameters parameters = new FitParameters();
      final Spot[] spots = new Spot[candidates.spots.length];
      for (int i = 0; i < spots.length; i++) {
        spots[i] = candidates.spots[i].spot;
      }
      // Debug candidates...
      // if (frame == 5)
      // System.out.printf("Fit %d [%d,%d = %.1f]\n", i, spots[i].x, spots[i].y,
      // spots[i].intensity);
      parameters.spots = spots;
      parameters.maxCandidate = candidates.maxCandidate;
      parameters.fitTask = FitTask.BENCHMARKING;
      parameters.benchmarkFilter = multiFilter;

      final ParameterisedFitJob job = new ParameterisedFitJob(parameters, frame, data, bounds);
      fitWorker.run(job); // Results will be stored in the fit job

      // TODO - Check how may candidates after the maxCandidate are fit (because a good estimate
      // was stored during multi-fit).
      // This may be taking a long time to perform. Somehow we need to speed up fitting during
      // the iteration.

      // Check why we need to store all the fit results.
      for (int i = 0; i < totalCandidates; i++) {
        fitResult[i] = job.getMultiPathFitResult(i);
      }

      // Compute the matches of the fitted spots to the simulated positions.
      // We will match all fitting results so providing the upper limit for the match score after
      // filtering.
      final Coordinate[] actual = ResultsMatchCalculator.getCoordinates(actualCoordinates, frame);
      final double[] zPosition = new double[actual.length];
      SpotMatch[] match = new SpotMatch[actual.length];
      int matchCount = 0;
      final RampedScore rampedScore =
          RampedScore.of(distanceInPixels, lowerDistanceInPixels, false);
      final RampedScore signalScore = (settings.signalFactor > 0)
          ? RampedScore.of(settings.signalFactor, settings.lowerSignalFactor, false)
          : null;
      if (actual.length > 0) {
        // Build a list of the coordinates z-depth using the PeakResultPoint
        for (int i = 0; i < actual.length; i++) {
          final PeakResultPoint p = (PeakResultPoint) actual[i];
          zPosition[i] = p.getPeakResult().getZPosition();
        }

        // Allow for doublets the predicted array
        // TODO - Check why we need to store all candidates in the predicted list. Those that are
        // not fitted are stored as MultiPathPoint.SPOT. Why do we need these?
        final ArrayList<MultiPathPoint> predicted = new ArrayList<>(totalCandidates * 2);
        matches.clear();

        for (int i = 0; i < totalCandidates; i++) {
          // Use all the results.
          // Store results using the candidate Id. The best match for each Id will be chosen.
          final int size = predicted.size();

          // Allow all new results to be processed as multi-fit can return more than 1 new result.
          add(predicted, fitResult[i].getSingleFitResult(), MultiPathPoint.SINGLE, i);
          add(predicted, fitResult[i].getMultiFitResult(), MultiPathPoint.MULTI, i);
          add(predicted, fitResult[i].getDoubletFitResult(), MultiPathPoint.DOUBLET, i);

          // TODO - why is this not MULTI_DOUBLET?

          add(predicted, fitResult[i].getMultiDoubletFitResult(), MultiPathPoint.MULTI, i);
          if (size == predicted.size()) {
            // Use the candidate position instead
            // TODO - Why store these?
            predicted.add(new MultiPathPoint(spots[i].x + 0.5f, spots[i].y + 0.5f, i,
                MultiPathPoint.SPOT, i));
          }
        }
        // If we made any fits then score them
        if (!predicted.isEmpty()) {
          // Match fit results/candidates with their closest actual spot

          // TODO - Is using the closest match the best way to do this for high density data?
          // Perhaps we should pair up the closest matches using the signal factor as well.

          final double matchDistance = distanceInPixels * distanceInPixels;
          final ArrayList<Assignment> assignments = new ArrayList<>();

          // Match all the fit results to spots. We want to match all fit results to actual spots.
          // All remaining candidate spots can then be matched to any remaining actual spots.
          final ResultAssignment[][] resultAssignments = new ResultAssignment[predicted.size()][];
          final int[] resultAssignmentsSize = new int[resultAssignments.length];
          for (int ii = 0; ii < actual.length; ii++) {
            final float x = actual[ii].getX();
            final float y = actual[ii].getY();
            for (int jj = 0; jj < predicted.size(); jj++) {
              final double d2 = predicted.get(jj).distanceSquared(x, y);
              if (d2 <= matchDistance) {
                // Get the score
                double distance = d2;
                double score = 0;

                if (predicted.get(jj).type != MultiPathPoint.SPOT) {
                  // Use the signal and ramped distance scoring
                  score = rampedScore.score(Math.sqrt(d2));
                  if (signalScore != null) {
                    final PeakResultPoint p3 = (PeakResultPoint) actual[ii];
                    // Assume the simulation is in photons
                    final double sf = computeSignalFactor(predicted.get(jj).result.getSignal(),
                        p3.getPeakResult().getIntensity());
                    score *= signalScore.score(Math.abs(sf));

                    if (score == 0) {
                      // This doesn't match
                      continue;
                    }
                  }

                  // Invert for the ranking (i.e. low is best)
                  distance = 1 - score;

                  // Ensure a perfect match can still be ranked ... closest first
                  if (distance == 0) {
                    distance -= (matchDistance - d2);
                  }

                  // Store the assignments for this result for later filter analysis
                  if (resultAssignments[jj] == null) {
                    resultAssignments[jj] = new ResultAssignment[5];
                  } else if (resultAssignments[jj].length == resultAssignmentsSize[jj]) {
                    resultAssignments[jj] =
                        Arrays.copyOf(resultAssignments[jj], resultAssignmentsSize[jj] * 2);
                  }
                  resultAssignments[jj][resultAssignmentsSize[jj]++] =
                      new ResultAssignment(ii, distance, score);
                } else {
                  // This is not a fit. Ensure that remaining candidate spots are assigned
                  // after any fit results.
                  distance += matchDistance + 1;
                }

                // Store the index of the predicted point in the score
                assignments
                    .add(new ImmutableFractionalAssignment(ii, predicted.get(jj).id, distance, jj));
              }
            }
          }

          // Set assignments
          for (int jj = 0; jj < predicted.size(); jj++) {
            if (predicted.get(jj).type == MultiPathPoint.SPOT) {
              continue;
            }
            if (resultAssignmentsSize[jj] != 0) {
              final BasePreprocessedPeakResult result =
                  (BasePreprocessedPeakResult) predicted.get(jj).result;
              result
                  .setAssignments(Arrays.copyOf(resultAssignments[jj], resultAssignmentsSize[jj]));
            }
          }

          AssignmentComparator.sort(assignments);

          final boolean[] actualAssignment = new boolean[actual.length];
          // We use the candidate Id as the id
          final boolean[] predictedAssignment = new boolean[fitResult.length];

          for (final Assignment assignment : assignments) {
            if (!actualAssignment[assignment.getTargetId()]
                && !predictedAssignment[assignment.getPredictedId()]) {
              actualAssignment[assignment.getTargetId()] = true;
              predictedAssignment[assignment.getPredictedId()] = true;

              final PeakResultPoint p3 = (PeakResultPoint) actual[assignment.getTargetId()];
              final int jj = (int) ((ImmutableFractionalAssignment) assignment).getScore();
              final MultiPathPoint point = predicted.get(jj);

              final double d = point.distanceXy(p3);

              if (point.type != MultiPathPoint.SPOT) {
                // This is a fitted candidate

                final double a = p3.getPeakResult().getIntensity(); // Should be in photons
                final double p = point.result.getSignal();

                match[matchCount++] =
                    new FitMatch(point, d, p3.getPeakResult().getZPosition(), p, a);
              } else {
                // This is a candidate that could not be fitted
                match[matchCount++] =
                    new CandidateMatch(point.index, d, p3.getPeakResult().getZPosition());
              }
            }
          }
        }
      }
      match = Arrays.copyOf(match, matchCount);

      // Mark the results
      double tp = 0;
      double fp = 0;
      double tn = 0;
      double fn = 0;

      for (int i = 0; i < match.length; i++) {
        // Score using just the distance.
        // The signal has been used only to compute the best match when two spots are close.
        final double s = rampedScore.scoreAndFlatten(match[i].distance, 256);
        match[i].score = s;

        if (match[i].isFitResult()) {
          // This is a fitted result so is a positive
          tp += s;
          fp += 1 - s;
        } else {
          // This is an unfitted result that matches so is a negative
          fn += s;
          tn += 1 - s;
        }
      }

      // Store the results using a copy of the original (to preserve the candidates for repeat
      // analysis)
      candidates = candidates.clone();
      candidates.tp = tp;
      candidates.fp = fp;
      candidates.tn = tn;
      candidates.fn = fn;
      candidates.fitResult = fitResult;
      candidates.match = match;
      candidates.zPosition = zPosition;
      candidates.noise = fitWorker.getNoise();
      results.put(frame, candidates);
    }

    private void add(ArrayList<MultiPathPoint> predicted,
        uk.ac.sussex.gdsc.smlm.results.filter.MultiPathFitResult.FitResult fitResult, int type,
        int spotId) {
      if (fitResult == null || fitResult.status != 0 || fitResult.getResults() == null) {
        return;
      }

      for (int i = 0; i < fitResult.getResults().length; i++) {
        if (fitResult.getResults()[i].isNewResult()) {
          predicted.add(new MultiPathPoint(fitResult.getResults()[i],
              fitResult.getResults()[i].getCandidateId(), type, spotId));
        }
      }
    }
  }

  /**
   * Create a new instance.
   */
  public BenchmarkSpotFit() {
    // Get the latest setting (by copy)
    settings = Settings.load();
    config = settings.config;
    multiFilter = settings.multiFilter;
    // Save (by reference)
    settings.save();
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    extraOptions = ImageJUtils.isExtraOptions();

    silent = false;
    finished = false;
    if (!initialise()) {
      return;
    }

    if (!showDialog()) {
      return;
    }

    runFitting();
    finished = true;
  }

  /**
   * Run the analysis non-interactively using the given filter settings.
   *
   * @param filter the filter
   * @param residualsThreshold the residuals threshold
   * @param failuresLimit the failures limit
   * @param duplicateDistance the duplicate distance
   * @param duplicateDistanceAbsolute the duplicate distance absolute
   */
  void run(DirectFilter filter, double residualsThreshold, int failuresLimit,
      double duplicateDistance, boolean duplicateDistanceAbsolute) {
    clearFitResults();

    silent = true;
    finished = false;
    if (!initialise()) {
      return;
    }

    // Reset the filter
    multiFilter = new MultiPathFilter(filter,
        FitWorker.createMinimalFilter(getPrecisionMethod(filter)), residualsThreshold);
    // Update the appropriate fitting settings
    config.setFailuresLimit(failuresLimit);
    config.setDuplicateDistance(duplicateDistance);
    config.setDuplicateDistanceAbsolute(duplicateDistanceAbsolute);

    runFitting();
    finished = true;
  }

  private boolean initialise() {
    simulationParameters = CreateData.getSimulationParameters();
    if (simulationParameters == null) {
      if (!silent) {
        IJ.error(TITLE, "No benchmark spot parameters in memory");
      }
      return false;
    }
    imp = CreateData.getImage();
    if (imp == null) {
      if (!silent) {
        IJ.error(TITLE, "No benchmark image");
      }
      return false;
    }
    results = CreateData.getResults();
    if (results == null) {
      if (!silent) {
        IJ.error(TITLE, "No benchmark results in memory");
      }
      return false;
    }
    filterResult = BenchmarkSpotFilter.getBenchmarkFilterResult();
    if (filterResult == null) {
      if (!silent) {
        IJ.error(TITLE, "No benchmark spot candidates in memory");
      }
      return false;
    }
    if (filterResult.simulationId != simulationParameters.id) {
      if (!silent) {
        IJ.error(TITLE, "Update the benchmark spot candidates for the latest simulation");
      }
      return false;
    }
    // This is required to initialise the FitWorker
    spotFilter = filterResult.spotFilter;
    return true;
  }

  private boolean showDialog() {
    ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("fit-spot-data"));

    ImageJUtils.addMessage(gd,
        "Fit candidate spots in the benchmark image created by " + CreateData.TITLE
            + " plugin\nand identified by the " + BenchmarkSpotFilter.TITLE
            + " plugin.\nPSF width = %s nm (Square pixel adjustment = %s nm)\n \n"
            + "Configure the fitting:",
        MathUtils.rounded(simulationParameters.sd), MathUtils.rounded(getSa()));

    gd.addSlider("Fraction_positives", 50, 100, settings.fractionPositives);
    gd.addSlider("Fraction_negatives_after_positives", 0, 100,
        settings.fractionNegativesAfterAllPositives);
    gd.addSlider("Min_negatives_after_positives", 0, 30, settings.negativesAfterAllPositives);
    gd.addSlider("Match_distance", 0.5, 3.5, settings.distance);
    gd.addSlider("Lower_distance", 0, 3.5, settings.lowerDistance);
    gd.addSlider("Match_signal", 0, 3.5, settings.signalFactor);
    gd.addSlider("Lower_signal", 0, 3.5, settings.lowerSignalFactor);

    final FitEngineConfigurationProvider fitEngineConfigurationProvider =
        new PeakFit.SimpleFitEngineConfigurationProvider(config);

    // Collect options for fitting
    final double sa = getSa();
    final FitConfiguration fitConfig = config.getFitConfiguration();
    fitConfig.setInitialPeakStdDev(MathUtils.round(sa / simulationParameters.pixelPitch));
    PeakFit.addPsfOptions(gd, fitConfig);
    PeakFit.addFittingOptions(gd, fitEngineConfigurationProvider);
    gd.addChoice("Fit_solver", SettingsManager.getFitSolverNames(),
        FitProtosHelper.getName(fitConfig.getFitSolver()));

    gd.addMessage("Multi-path filter (used to pick optimum results during fitting)");

    // Allow loading the best filter for these results
    final boolean benchmarkSettingsCheckbox =
        BenchmarkSpotFitResult.fitResultsId.get() == BenchmarkFilterAnalysis.getLastFittingId();

    Checkbox cbBenchmark = null;
    if (benchmarkSettingsCheckbox) {
      // This should always be an opt-in decision. Otherwise the user cannot use the previous
      // settings.
      cbBenchmark = gd.addAndGetCheckbox("Benchmark_settings", false);
    }

    gd.addTextAreas(XmlUtils.convertQuotes(multiFilter.toXml()), null, 6, 60);

    textFailLimit = gd.addAndGetNumericField("Fail_limit", config.getFailuresLimit(), 0);
    cbIncludeNeighbours = gd.addAndGetCheckbox("Include_neighbours", config.isIncludeNeighbours());
    gd.addAndGetSlider("Neighbour_height", 0.01, 1, config.getNeighbourHeightThreshold());
    textNeighbourHeight = gd.getLastTextField();
    cbComputeDoublets = gd.addAndGetCheckbox("Compute_doublets", settings.computeDoublets);
    PeakFit.addDuplicateDistanceOptions(gd, fitEngineConfigurationProvider);
    gd.addCheckbox("Show_score_histograms", settings.showFilterScoreHistograms);
    gd.addCheckbox("Show_correlation", settings.showCorrelation);
    gd.addCheckbox("Plot_rank_by_intensity", settings.rankByIntensity);
    gd.addCheckbox("Save_filter_range", settings.saveFilterRange);

    if (extraOptions) {
      // No extra options
    }

    // Add a mouse listener to the config file field
    if (cbBenchmark != null && ImageJUtils.isShowGenericDialog()) {
      taFilterXml = gd.getTextArea1();
      cbBenchmark.addItemListener(this);
    }

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.fractionPositives = Math.abs(gd.getNextNumber());
    settings.fractionNegativesAfterAllPositives = Math.abs(gd.getNextNumber());
    settings.negativesAfterAllPositives = (int) Math.abs(gd.getNextNumber());
    settings.distance = Math.abs(gd.getNextNumber());
    settings.lowerDistance = Math.abs(gd.getNextNumber());
    settings.signalFactor = Math.abs(gd.getNextNumber());
    settings.lowerSignalFactor = Math.abs(gd.getNextNumber());

    fitConfig.setPsfType(PeakFit.getPsfTypeValues()[gd.getNextChoiceIndex()]);
    config.setFitting(gd.getNextNumber());
    // Some enum values are not supported
    fitConfig.setFitSolver(SettingsManager.getFitSolverValues()[gd.getNextChoiceIndex()]);

    boolean myUseBenchmarkSettings = false;
    if (benchmarkSettingsCheckbox) {
      // useBenchmarkSettings =
      myUseBenchmarkSettings = gd.getNextBoolean();
    }

    // Read dialog settings
    final String xml = gd.getNextText();
    final int failLimit = (int) gd.getNextNumber();
    final boolean includeNeighbours = gd.getNextBoolean();
    final double neighbourHeightThreshold = gd.getNextNumber();
    final boolean myComputeDoublets = gd.getNextBoolean();
    final double myDuplicateDistance = gd.getNextNumber();

    gd.collectOptions();

    MultiPathFilter myMultiFilter = null;
    if (myUseBenchmarkSettings && !ImageJUtils.isShowGenericDialog()) {
      // Only copy the benchmark settings if not interactive
      final FitEngineConfiguration tmp = FitEngineConfiguration.create();
      final FitConfiguration tmpFitConfig = tmp.getFitConfiguration();
      tmpFitConfig.setComputeResiduals(true); // Collect the residuals threshold
      if (BenchmarkFilterAnalysis.updateConfiguration(tmp, false)) {
        config.setFailuresLimit(tmp.getFailuresLimit());
        config.setIncludeNeighbours(tmp.isIncludeNeighbours());
        config.setNeighbourHeightThreshold(tmp.getNeighbourHeightThreshold());
        settings.computeDoublets = (tmp.getResidualsThreshold() < 1);
        config.setDuplicateDistance(tmp.getDuplicateDistance());
        config.setDuplicateDistanceAbsolute(tmp.getDuplicateDistanceAbsolute());

        final DirectFilter primaryFilter = tmpFitConfig.getSmartFilter();
        final double residualsThreshold = tmp.getResidualsThreshold();
        myMultiFilter = new MultiPathFilter(primaryFilter,
            FitWorker.createMinimalFilter(tmpFitConfig.getFilterPrecisionMethod()),
            residualsThreshold);
      }
    } else {
      myMultiFilter = MultiPathFilter.fromXml(xml);

      config.setFailuresLimit(failLimit);
      config.setIncludeNeighbours(includeNeighbours);
      config.setNeighbourHeightThreshold(neighbourHeightThreshold);
      settings.computeDoublets = myComputeDoublets;
      config.setDuplicateDistance(myDuplicateDistance);
    }

    if (myMultiFilter == null) {
      gd = new ExtendedGenericDialog(TITLE);
      gd.addMessage("The multi-path filter was invalid.\n \nContinue with a default filter?");
      gd.enableYesNoCancel();
      gd.hideCancelButton();
      gd.showDialog();
      if (!gd.wasOKed()) {
        return false;
      }
    } else {
      multiFilter = myMultiFilter;
    }

    if (settings.computeDoublets) {
      config.setResidualsThreshold(0);
      fitConfig.setComputeResiduals(true);
    } else {
      config.setResidualsThreshold(1);
      fitConfig.setComputeResiduals(false);
    }
    settings.showFilterScoreHistograms = gd.getNextBoolean();
    settings.showCorrelation = gd.getNextBoolean();
    settings.rankByIntensity = gd.getNextBoolean();
    settings.saveFilterRange = gd.getNextBoolean();

    // Avoid stupidness, i.e. things that move outside the fit window and are bad widths

    // TODO - Fix this for simple or smart filter...
    fitConfig.setDisableSimpleFilter(false);

    fitConfig.setMinPhotons(15); // Realistically we cannot fit lower than this
    // Disable shift as candidates may be re-mapped to alternative candidates so the initial
    // position is wrong.
    fitConfig.setCoordinateShiftFactor(0);
    fitConfig.setMinWidthFactor(1.0 / 5);
    fitConfig.setMaxWidthFactor(5);
    // Disable the direct filter
    fitConfig.setDirectFilter(null);

    if (extraOptions) {
      // No extra options
    }

    if (gd.invalidNumber()) {
      return false;
    }

    if (settings.lowerDistance > settings.distance) {
      settings.lowerDistance = settings.distance;
    }
    if (settings.lowerSignalFactor > settings.signalFactor) {
      settings.lowerSignalFactor = settings.signalFactor;
    }

    // Distances relative to sa (not s) as this is the same as the BenchmarkSpotFilter plugin
    distanceInPixels = settings.distance * sa / simulationParameters.pixelPitch;
    lowerDistanceInPixels = settings.lowerDistance * sa / simulationParameters.pixelPitch;

    // Copy simulation defaults if a new simulation.
    // if (lastSimulationId.get() != simulationParameters.id) {
    // This is needed to configure the fit solver.
    fitConfig.setNmPerPixel(simulationParameters.pixelPitch);
    fitConfig.setGain(simulationParameters.gain);
    fitConfig.setQuantumEfficiency(simulationParameters.qe);
    fitConfig.setReadNoise(simulationParameters.readNoise);
    fitConfig.setBias(simulationParameters.bias);
    fitConfig.setCameraType(simulationParameters.cameraType);
    fitConfig.setCameraModel(CreateData.getCameraModel(simulationParameters));
    // }
    if (!PeakFit.configurePsfModel(config)) {
      return false;
    }
    return PeakFit.configureFitSolver(config, IJImageSource.getBounds(imp), null,
        (extraOptions) ? PeakFit.FLAG_EXTRA_OPTIONS : 0);
  }

  private BenchmarkSpotFitResult runFitting() {
    // Extract all the results in memory into a list per frame. This can be cached
    boolean refresh = false;
    Pair<Integer, Int2ObjectOpenHashMap<List<Coordinate>>> coords = coordinateCache.get();

    if (coords.getKey() != simulationParameters.id) {
      // Do not get integer coordinates
      // The Coordinate objects will be PeakResultPoint objects that store the original PeakResult
      // from the MemoryPeakResults
      coords =
          Pair.of(simulationParameters.id, ResultsMatchCalculator.getCoordinates(results, false));
      coordinateCache.set(coords);
      refresh = true;
    }

    final Int2ObjectOpenHashMap<List<Coordinate>> actualCoordinates = coords.getValue();

    // Extract all the candidates into a list per frame. This can be cached if the settings have not
    // changed
    final int width = (config.isIncludeNeighbours()) ? config.getFittingWidth() : 0;
    CandidateData candidateData = candidateDataCache.get();
    if (refresh || candidateData == null
        || candidateData.differentSettings(filterResult.id, settings, width)) {
      candidateData = subsetFilterResults(filterResult.filterResults, width);
      candidateDataCache.set(candidateData);
    }

    final StopWatch stopWatch = StopWatch.createStarted();
    final ImageStack stack = imp.getImageStack();

    clearFitResults();

    // Save results to memory
    final MemoryPeakResults peakResults = new MemoryPeakResults();
    peakResults.copySettings(this.results);
    peakResults.setName(TITLE);
    config.configureOutputUnits();
    final FitConfiguration fitConfig = config.getFitConfiguration();
    peakResults.setCalibration(fitConfig.getCalibration());
    MemoryPeakResults.addResults(peakResults);

    // Create a pool of workers
    final int nThreads = Prefs.getThreads();
    final BlockingQueue<Integer> jobs = new ArrayBlockingQueue<>(nThreads * 2);
    final List<Worker> workers = new LinkedList<>();
    final List<Thread> threads = new LinkedList<>();
    final Ticker ticker = ImageJUtils.createTicker(stack.getSize(), nThreads, "Fitting frames ...");
    final PeakResults syncResults = SynchronizedPeakResults.create(peakResults, nThreads);
    for (int i = 0; i < nThreads; i++) {
      final Worker worker = new Worker(jobs, stack, actualCoordinates,
          candidateData.filterCandidates, syncResults, ticker);
      final Thread t = new Thread(worker);
      workers.add(worker);
      threads.add(t);
      t.start();
    }

    // Fit the frames
    final long startTime = System.nanoTime();
    for (int i = 1; i <= stack.getSize(); i++) {
      put(jobs, i);
    }
    // Finish all the worker threads by passing in a null job
    for (int i = 0; i < threads.size(); i++) {
      put(jobs, -1);
    }

    // Wait for all to finish
    for (int i = 0; i < threads.size(); i++) {
      try {
        threads.get(i).join();
      } catch (final InterruptedException ex) {
        Thread.currentThread().interrupt();
        throw new ConcurrentRuntimeException(ex);
      }
    }
    final long runTime = System.nanoTime() - startTime;

    threads.clear();

    ImageJUtils.finished();

    if (ImageJUtils.isInterrupted()) {
      return null;
    }

    stopWatch.stop();
    final String timeString = stopWatch.toString();
    IJ.log("Spot fit time : " + timeString);

    IJ.showStatus("Collecting results ...");

    if (fitConfig.isFitCameraCounts()) {
      // Convert to photons for consistency
      results.convertToPreferredUnits();
    }

    final CustomInt2ObjectOpenHashMap<FilterCandidates> fitResults =
        new CustomInt2ObjectOpenHashMap<>();
    for (final Worker w : workers) {
      fitResults.putAll(w.results);
    }

    // Assign a unique ID to each result
    int count = 0;
    // Materialise into an array since we use it twice
    final FilterCandidates[] candidates =
        fitResults.values().toArray(new FilterCandidates[fitResults.size()]);
    for (final FilterCandidates result : candidates) {
      for (int i = 0; i < result.fitResult.length; i++) {
        final MultiPathFitResult fitResult = result.fitResult[i];
        count += count(fitResult.getSingleFitResult());
        count += count(fitResult.getMultiFitResult());
        count += count(fitResult.getDoubletFitResult());
        count += count(fitResult.getMultiDoubletFitResult());
      }
    }
    final PreprocessedPeakResult[] preprocessedPeakResults = new PreprocessedPeakResult[count];
    count = 0;
    for (final FilterCandidates result : candidates) {
      for (int i = 0; i < result.fitResult.length; i++) {
        final MultiPathFitResult fitResult = result.fitResult[i];
        count = store(fitResult.getSingleFitResult(), count, preprocessedPeakResults);
        count = store(fitResult.getMultiFitResult(), count, preprocessedPeakResults);
        count = store(fitResult.getDoubletFitResult(), count, preprocessedPeakResults);
        count = store(fitResult.getMultiDoubletFitResult(), count, preprocessedPeakResults);
      }
    }

    final BenchmarkSpotFitResult newSpotFitResults =
        new BenchmarkSpotFitResult(simulationParameters.id, fitResults);

    newSpotFitResults.distanceInPixels = distanceInPixels;
    newSpotFitResults.lowerDistanceInPixels = lowerDistanceInPixels;
    newSpotFitResults.stopWatch = stopWatch;

    summariseResults(newSpotFitResults, runTime, preprocessedPeakResults, count, candidateData,
        actualCoordinates);

    IJ.showStatus("");

    spotFitResults.set(newSpotFitResults);

    return newSpotFitResults;
  }

  /**
   * Clear old results to free memory.
   */
  private static void clearFitResults() {
    spotFitResults.set(null);
  }

  private static int
      count(uk.ac.sussex.gdsc.smlm.results.filter.MultiPathFitResult.FitResult fitResult) {
    if (fitResult == null || fitResult.getResults() == null) {
      return 0;
    }
    int count = 0;
    for (int i = 0; i < fitResult.getResults().length; i++) {
      final PreprocessedPeakResult result = fitResult.getResults()[i];
      if (result.isNewResult()) {
        count++;
      }
    }
    return count;
  }

  private static int store(
      uk.ac.sussex.gdsc.smlm.results.filter.MultiPathFitResult.FitResult fitResult, int count,
      PreprocessedPeakResult[] preprocessedPeakResults) {
    if (fitResult == null || fitResult.getResults() == null) {
      return count;
    }
    for (int i = 0; i < fitResult.getResults().length; i++) {
      final BasePreprocessedPeakResult result =
          (BasePreprocessedPeakResult) fitResult.getResults()[i];
      if (result.isNewResult()) {
        result.uniqueId = count++;
        preprocessedPeakResults[result.uniqueId] = result;
      }
    }
    return count;
  }

  /**
   * Extract all the filter candidates in order until the desired number of positives have been
   * reached and the number of negatives matches the configured parameters.
   *
   * @param filterResults the filter results
   * @param fitting the fitting
   * @return The filter candidate data
   */
  private CandidateData subsetFilterResults(CustomInt2ObjectOpenHashMap<FilterResult> filterResults,
      int fitting) {
    // Convert fractions from percent
    final double f1 = Math.min(1, settings.fractionPositives / 100.0);
    final double f2 = settings.fractionNegativesAfterAllPositives / 100.0;

    final int[] counter = new int[2];

    final Int2ObjectOpenHashMap<FilterCandidates> subset = new Int2ObjectOpenHashMap<>();
    final double[] fX = new double[2];
    final int[] nX = new int[2];
    filterResults.forEach((int frame, FilterResult result) -> {
      // Determine the number of positives to find. This score may be fractional.
      fX[0] += result.result.getTruePositives();
      fX[1] += result.result.getFalsePositives();

      // Q. Is r.result.getTruePositives() not the same as the total of r.spots[i].match?
      // A. Not if we used fractional scoring.
      int count = 0;
      for (int i = result.spots.length; i-- > 0;) {
        if (result.spots[i].match) {
          count++;
        }
      }
      nX[0] += count;
      nX[1] += (result.spots.length - count);

      // Make the target use the fractional score
      final double np2 = result.result.getTruePositives() * f1;
      double targetP = np2;

      // Set the target using the closest
      if (f1 < 1) {
        double np = 0;
        double min = result.result.getTruePositives();
        for (final ScoredSpot spot : result.spots) {
          if (spot.match) {
            np += spot.getScore();
            final double d = np2 - np;
            if (d < min) {
              min = d;
              targetP = np;
            } else {
              break;
            }
          }
        }
      }

      // Count the number of positive & negatives
      int pos = 0;
      int neg = 0;
      double np = 0;
      double nn = 0;

      boolean reachedTarget = false;
      int countAfter = 0;

      count = 0;
      for (final ScoredSpot spot : result.spots) {
        count++;
        nn += spot.antiScore();
        if (spot.match) {
          np += spot.getScore();
          pos++;
          if (!reachedTarget) {
            reachedTarget = np >= targetP;
          }
        } else {
          neg++;
          if (reachedTarget) {
            countAfter++;
          }
        }

        // Check if we have reached both the limits
        if (reachedTarget && countAfter >= settings.negativesAfterAllPositives
            && (double) neg / (neg + pos) >= f2) {
          break;
        }
      }

      counter[0] += count;
      counter[1] += result.spots.length;

      // We can use all the candidates but only fit up to count
      subset.put(frame, new FilterCandidates(pos, neg, np, nn, result.spots, count));
    });

    // We now add all the candidates but only fit the first N
    final int target = counter[0];
    final int total = counter[1];
    final int added = total - target;

    if (extraOptions && added > target) {
      ImageJUtils.log("Added %s to %s (total = %d)", TextUtils.pleural(added, "neighbour"),
          TextUtils.pleural(target, "candidate"), total);
    }

    return new CandidateData(subset, filterResult.id, fX[0], fX[1], nX[0], nX[1], settings,
        fitting);
  }

  private static void put(BlockingQueue<Integer> jobs, int slice) {
    try {
      jobs.put(slice);
    } catch (final InterruptedException ex) {
      Thread.currentThread().interrupt();
      throw new ConcurrentRuntimeException("Unexpected interruption", ex);
    }
  }

  /**
   * Create an abstract class to allow a count to be passed to the constructor. The consumer can be
   * coded inline using final object references.
   */
  private abstract static class CustomIntConsumer implements IntConsumer {
    int count;

    CustomIntConsumer(int count) {
      this.count = count;
    }
  }

  private void summariseResults(BenchmarkSpotFitResult spotFitResults, long runTime,
      final PreprocessedPeakResult[] preprocessedPeakResults, int uniqueIdCount,
      CandidateData candidateData, Int2ObjectOpenHashMap<List<Coordinate>> actualCoordinates) {

    // Summarise the fitting results. N fits, N failures.
    // Optimal match statistics if filtering is perfect (since fitting is not perfect).
    final StoredDataStatistics distanceStats = new StoredDataStatistics();
    final StoredDataStatistics depthStats = new StoredDataStatistics();

    // Get stats for all fitted results and those that match
    // Signal, SNR, Width, xShift, yShift, Precision
    createFilterCriteria();
    final StoredDataStatistics[][] stats = new StoredDataStatistics[3][filterCriteria.length];
    for (int i = 0; i < stats.length; i++) {
      for (int j = 0; j < stats[i].length; j++) {
        stats[i][j] = new StoredDataStatistics();
      }
    }

    final double nmPerPixel = simulationParameters.pixelPitch;
    double tp = 0;
    double fp = 0;
    int failCtp = 0;
    int failCfp = 0;
    int ctp = 0;
    int cfp = 0;
    final int[] singleStatus = new int[FitStatus.values().length];
    final int[] multiStatus = new int[singleStatus.length];
    final int[] doubletStatus = new int[singleStatus.length];
    final int[] multiDoubletStatus = new int[singleStatus.length];

    // Easier to materialise the values since we have a lot of non final variables to manipulate
    final CustomInt2ObjectOpenHashMap<FilterCandidates> fitResults = spotFitResults.fitResults;
    final int[] frames = new int[fitResults.size()];
    final FilterCandidates[] candidates = new FilterCandidates[fitResults.size()];
    final int[] counter = new int[1];
    fitResults.forEach((int key, FilterCandidates value) -> {
      frames[counter[0]] = key;
      candidates[counter[0]] = value;
      counter[0]++;
    });

    for (final FilterCandidates result : candidates) {
      // Count the number of fit results that matched (tp) and did not match (fp)
      tp += result.tp;
      fp += result.fp;

      for (int i = 0; i < result.fitResult.length; i++) {
        if (result.spots[i].match) {
          ctp++;
        } else {
          cfp++;
        }
        final MultiPathFitResult fitResult = result.fitResult[i];

        if (singleStatus != null && result.spots[i].match) {
          // Debugging reasons for fit failure
          addStatus(singleStatus, fitResult.getSingleFitResult());
          addStatus(multiStatus, fitResult.getMultiFitResult());
          addStatus(doubletStatus, fitResult.getDoubletFitResult());
          addStatus(multiDoubletStatus, fitResult.getMultiDoubletFitResult());
        }

        if (noMatch(fitResult)) {
          if (result.spots[i].match) {
            failCtp++;
          } else {
            failCfp++;
          }
        }

        // We have multi-path results.
        // We want statistics for:
        // [0] all fitted spots
        // [1] fitted spots that match a result
        // [2] fitted spots that do not match a result
        addToStats(fitResult.getSingleFitResult(), stats);
        addToStats(fitResult.getMultiFitResult(), stats);
        addToStats(fitResult.getDoubletFitResult(), stats);
        addToStats(fitResult.getMultiDoubletFitResult(), stats);
      }

      // Statistics on spots that fit an actual result
      for (int i = 0; i < result.match.length; i++) {
        if (!result.match[i].isFitResult()) {
          // For now just ignore the candidates that matched
          continue;
        }

        final FitMatch fitMatch = (FitMatch) result.match[i];
        distanceStats.add(fitMatch.distance * nmPerPixel);
        depthStats.add(fitMatch.zdepth * nmPerPixel);
      }
    }

    if (tp == 0) {
      IJ.error(TITLE, "No fit results matched the simulation actual results");
      return;
    }

    // Store data for computing correlation
    final double[] i1 = new double[depthStats.getN()];
    final double[] i2 = new double[i1.length];
    final double[] is = new double[i1.length];
    int ci = 0;
    for (final FilterCandidates result : candidates) {
      for (int i = 0; i < result.match.length; i++) {
        if (!result.match[i].isFitResult()) {
          // For now just ignore the candidates that matched
          continue;
        }

        final FitMatch fitMatch = (FitMatch) result.match[i];
        final ScoredSpot spot = result.spots[fitMatch.index];
        i1[ci] = fitMatch.predictedSignal;
        i2[ci] = fitMatch.actualSignal;
        is[ci] = spot.spot.intensity;
        ci++;
      }
    }

    // We want to compute the Jaccard against the spot metric

    // Filter the results using the multi-path filter
    final ArrayList<MultiPathFitResults> multiPathResults = new ArrayList<>(fitResults.size());
    for (int i = 0; i < frames.length; i++) {
      final int frame = frames[i];
      final MultiPathFitResult[] multiPathFitResults = candidates[i].fitResult;
      final int totalCandidates = candidates[i].spots.length;
      final List<Coordinate> list = actualCoordinates.get(frame);
      final int nActual = (list == null) ? 0 : list.size();
      multiPathResults
          .add(new MultiPathFitResults(frame, multiPathFitResults, totalCandidates, nActual));
    }
    // Score the results and count the number returned
    final List<FractionalAssignment[]> assignments = new ArrayList<>();
    final IntOpenHashSet set = new IntOpenHashSet(uniqueIdCount);
    final FractionScoreStore scoreStore = set::add;

    final MultiPathFitResults[] multiResults = multiPathResults.toArray(new MultiPathFitResults[0]);
    // Filter with no filter
    final MultiPathFilter mpf =
        new MultiPathFilter(new SignalFilter(0), null, multiFilter.residualsThreshold);
    mpf.fractionScoreSubset(multiResults, NullFailCounter.INSTANCE, this.results.size(),
        assignments, scoreStore, CoordinateStoreFactory.create(0, 0, imp.getWidth(),
            imp.getHeight(), config.convertUsingHwhMax(config.getDuplicateDistanceParameter())));

    final double[][] matchScores = new double[set.size()][];
    int count = 0;
    for (int i = 0; i < assignments.size(); i++) {
      final FractionalAssignment[] a = assignments.get(i);
      if (a == null) {
        continue;
      }
      for (int j = 0; j < a.length; j++) {
        final PreprocessedPeakResult r = ((PeakFractionalAssignment) a[j]).peakResult;
        set.remove(r.getUniqueId());

        final double precision = Math.sqrt(r.getLocationVariance());
        final double signal = r.getSignal();
        final double snr = r.getSnr();
        final double width = r.getXSdFactor();
        final double xShift = r.getXRelativeShift2();
        final double yShift = r.getYRelativeShift2();
        // Since these two are combined for filtering and the max is what matters.
        final double shift = (xShift > yShift) ? Math.sqrt(xShift) : Math.sqrt(yShift);
        final double eshift = Math.sqrt(xShift + yShift);

        final double[] score = new double[8];
        score[FILTER_SIGNAL] = signal;
        score[FILTER_SNR] = snr;
        score[FILTER_MIN_WIDTH] = width;
        score[FILTER_MAX_WIDTH] = width;
        score[FILTER_SHIFT] = shift;
        score[FILTER_ESHIFT] = eshift;
        score[FILTER_PRECISION] = precision;
        score[FILTER_PRECISION + 1] = a[j].getScore();
        matchScores[count++] = score;
      }
    }
    // Add the rest
    set.forEach(new CustomIntConsumer(count) {
      @Override
      public void accept(int uniqueId) {
        // This should not be null or something has gone wrong
        final PreprocessedPeakResult r = preprocessedPeakResults[uniqueId];
        if (r == null) {
          throw new IllegalArgumentException("Missing result: " + uniqueId);
        }
        final double precision = Math.sqrt(r.getLocationVariance());
        final double signal = r.getSignal();
        final double snr = r.getSnr();
        final double width = r.getXSdFactor();
        final double xShift = r.getXRelativeShift2();
        final double yShift = r.getYRelativeShift2();
        // Since these two are combined for filtering and the max is what matters.
        final double shift = (xShift > yShift) ? Math.sqrt(xShift) : Math.sqrt(yShift);
        final double eshift = Math.sqrt(xShift + yShift);

        final double[] score = new double[8];
        score[FILTER_SIGNAL] = signal;
        score[FILTER_SNR] = snr;
        score[FILTER_MIN_WIDTH] = width;
        score[FILTER_MAX_WIDTH] = width;
        score[FILTER_SHIFT] = shift;
        score[FILTER_ESHIFT] = eshift;
        score[FILTER_PRECISION] = precision;
        matchScores[count++] = score;
      }
    });

    final FitConfiguration fitConfig = config.getFitConfiguration();

    // Debug the reasons the fit failed
    if (singleStatus != null) {
      String name = PeakFit.getSolverName(fitConfig);
      if (fitConfig.getFitSolver() == FitSolver.MLE && fitConfig.isModelCamera()) {
        name += " Camera";
      }
      IJ.log("Failure counts: " + name);
      printFailures("Single", singleStatus);
      printFailures("Multi", multiStatus);
      printFailures("Doublet", doubletStatus);
      printFailures("Multi doublet", multiDoubletStatus);
    }

    final StringBuilder sb = new StringBuilder(300);

    // Add information about the simulation
    final double signal = simulationParameters.averageSignal;
    final int n = results.size();
    sb.append(imp.getStackSize()).append('\t');
    final int w = imp.getWidth();
    final int h = imp.getHeight();
    sb.append(w).append('\t');
    sb.append(h).append('\t');
    sb.append(n).append('\t');
    final double density = ((double) n / imp.getStackSize()) / (w * h)
        / (simulationParameters.pixelPitch * simulationParameters.pixelPitch / 1e6);
    sb.append(MathUtils.rounded(density)).append('\t');
    sb.append(MathUtils.rounded(signal)).append('\t');
    sb.append(MathUtils.rounded(simulationParameters.sd)).append('\t');
    sb.append(MathUtils.rounded(simulationParameters.pixelPitch)).append('\t');
    sb.append(MathUtils.rounded(simulationParameters.depth)).append('\t');
    sb.append(simulationParameters.fixedDepth).append('\t');
    sb.append(MathUtils.rounded(simulationParameters.gain)).append('\t');
    sb.append(MathUtils.rounded(simulationParameters.readNoise)).append('\t');
    sb.append(MathUtils.rounded(simulationParameters.background)).append('\t');
    sb.append(MathUtils.rounded(simulationParameters.noise)).append('\t');

    if (simulationParameters.fullSimulation) {
      // The total signal is spread over frames
    }

    sb.append(MathUtils.rounded(signal / simulationParameters.noise)).append('\t');
    sb.append(MathUtils.rounded(simulationParameters.sd / simulationParameters.pixelPitch))
        .append('\t');

    sb.append(spotFilter.getDescription());

    // nP and nN is the fractional score of the spot candidates
    addCount(sb, (double) candidateData.countPositive + candidateData.countNegative);
    addCount(sb, candidateData.countPositive);
    addCount(sb, candidateData.countNegative);
    addCount(sb, candidateData.fractionPositive);
    addCount(sb, candidateData.fractionNegative);
    String name = PeakFit.getSolverName(fitConfig);
    if (fitConfig.getFitSolver() == FitSolver.MLE && fitConfig.isModelCamera()) {
      name += " Camera";
    }
    add(sb, name);
    add(sb, config.getFitting());

    spotFitResults.resultPrefix = sb.toString();

    // Q. Should I add other fit configuration here?

    // The fraction of positive and negative candidates that were included
    add(sb, (100.0 * ctp) / candidateData.countPositive);
    add(sb, (100.0 * cfp) / candidateData.countNegative);

    // Score the fitting results compared to the original simulation.

    // Score the candidate selection:
    add(sb, ctp + cfp);
    add(sb, ctp);
    add(sb, cfp);
    // TP are all candidates that can be matched to a spot
    // FP are all candidates that cannot be matched to a spot
    // FN = The number of missed spots
    FractionClassificationResult match =
        new FractionClassificationResult(ctp, cfp, 0, simulationParameters.molecules - ctp);
    add(sb, match.getRecall());
    add(sb, match.getPrecision());
    add(sb, match.getF1Score());
    add(sb, match.getJaccard());

    // Score the fitting results:
    add(sb, failCtp);
    add(sb, failCfp);

    // TP are all fit results that can be matched to a spot
    // FP are all fit results that cannot be matched to a spot
    // FN = The number of missed spots
    add(sb, tp);
    add(sb, fp);
    match = new FractionClassificationResult(tp, fp, 0, simulationParameters.molecules - tp);
    add(sb, match.getRecall());
    add(sb, match.getPrecision());
    add(sb, match.getF1Score());
    add(sb, match.getJaccard());

    // Do it again but pretend we can perfectly filter all the false positives
    // add(sb, tp);
    match = new FractionClassificationResult(tp, 0, 0, simulationParameters.molecules - tp);
    // Recall is unchanged
    // Precision will be 100%
    add(sb, match.getF1Score());
    add(sb, match.getJaccard());

    // The mean may be subject to extreme outliers so use the median
    double median = distanceStats.getMedian();
    add(sb, median);

    final WindowOrganiser wo = new WindowOrganiser();

    String label = String.format("Recall = %s. n = %d. Median = %s nm. SD = %s nm",
        MathUtils.rounded(match.getRecall()), distanceStats.getN(), MathUtils.rounded(median),
        MathUtils.rounded(distanceStats.getStandardDeviation()));
    new HistogramPlotBuilder(TITLE, distanceStats, "Match Distance (nm)").setPlotLabel(label)
        .show(wo);

    median = depthStats.getMedian();
    add(sb, median);

    // Sort by spot intensity and produce correlation
    double[] correlation = null;
    double[] rankCorrelation = null;
    double[] rank = null;
    final FastCorrelator fastCorrelator = new FastCorrelator();
    final ArrayList<Ranking> pc1 = new ArrayList<>();
    final ArrayList<Ranking> pc2 = new ArrayList<>();
    ci = 0;
    if (settings.showCorrelation) {
      final int[] indices = SimpleArrayUtils.natural(i1.length);
      SortUtils.sortData(indices, is, settings.rankByIntensity, true);
      correlation = new double[i1.length];
      rankCorrelation = new double[i1.length];
      rank = new double[i1.length];
      for (final int ci2 : indices) {
        fastCorrelator.add(Math.round(i1[ci2]), Math.round(i2[ci2]));
        pc1.add(new Ranking(i1[ci2], ci));
        pc2.add(new Ranking(i2[ci2], ci));
        correlation[ci] = fastCorrelator.getCorrelation();
        rankCorrelation[ci] = Correlator.correlation(rank(pc1), rank(pc2));
        if (settings.rankByIntensity) {
          rank[ci] = is[0] - is[ci];
        } else {
          rank[ci] = ci;
        }
        ci++;
      }
    } else {
      for (int i = 0; i < i1.length; i++) {
        fastCorrelator.add(Math.round(i1[i]), Math.round(i2[i]));
        pc1.add(new Ranking(i1[i], i));
        pc2.add(new Ranking(i2[i], i));
      }
    }

    final double pearsonCorr = fastCorrelator.getCorrelation();
    final double rankedCorr = Correlator.correlation(rank(pc1), rank(pc2));

    // Get the regression
    final SimpleRegression regression = new SimpleRegression(false);
    for (int i = 0; i < pc1.size(); i++) {
      regression.addData(pc1.get(i).value, pc2.get(i).value);
    }
    // final double intercept = regression.getIntercept();
    final double slope = regression.getSlope();

    if (settings.showCorrelation) {
      String title = TITLE + " Intensity";
      Plot plot = new Plot(title, "Candidate", "Spot");
      final double[] limits1 = MathUtils.limits(i1);
      final double[] limits2 = MathUtils.limits(i2);
      plot.setLimits(limits1[0], limits1[1], limits2[0], limits2[1]);
      label = String.format("Correlation=%s; Ranked=%s; Slope=%s", MathUtils.rounded(pearsonCorr),
          MathUtils.rounded(rankedCorr), MathUtils.rounded(slope));
      plot.addLabel(0, 0, label);
      plot.setColor(Color.red);
      plot.addPoints(i1, i2, Plot.DOT);
      if (slope > 1) {
        plot.drawLine(limits1[0], limits1[0] * slope, limits1[1], limits1[1] * slope);
      } else {
        plot.drawLine(limits2[0] / slope, limits2[0], limits2[1] / slope, limits2[1]);
      }
      ImageJUtils.display(title, plot, wo);

      title = TITLE + " Correlation";
      plot = new Plot(title, "Spot Rank", "Correlation");
      final double[] xlimits = MathUtils.limits(rank);
      double[] ylimits = MathUtils.limits(correlation);
      ylimits = MathUtils.limits(ylimits, rankCorrelation);
      plot.setLimits(xlimits[0], xlimits[1], ylimits[0], ylimits[1]);
      plot.setColor(Color.red);
      plot.addPoints(rank, correlation, Plot.LINE);
      plot.setColor(Color.blue);
      plot.addPoints(rank, rankCorrelation, Plot.LINE);
      plot.setColor(Color.black);
      plot.addLabel(0, 0, label);
      ImageJUtils.display(title, plot, wo);
    }

    add(sb, pearsonCorr);
    add(sb, rankedCorr);
    add(sb, slope);

    label = String.format("n = %d. Median = %s nm", depthStats.getN(), MathUtils.rounded(median));
    new HistogramPlotBuilder(TITLE, depthStats, "Match Depth (nm)").setRemoveOutliersOption(1)
        .setPlotLabel(label).show(wo);

    // Plot histograms of the stats on the same window
    final double[] lower = new double[filterCriteria.length];
    final double[] upper = new double[lower.length];
    final double[] min = new double[lower.length];
    final double[] max = new double[lower.length];
    for (int i = 0; i < stats[0].length; i++) {
      final double[] limits = showDoubleHistogram(stats, i, wo, matchScores);
      lower[i] = limits[0];
      upper[i] = limits[1];
      min[i] = limits[2];
      max[i] = limits[3];
    }

    // Reconfigure some of the range limits
    upper[FILTER_SIGNAL] *= 2; // Make this a bit bigger
    upper[FILTER_SNR] *= 2; // Make this a bit bigger
    final double factor = 0.25;
    if (lower[FILTER_MIN_WIDTH] != 0) {
      // (assuming lower is less than 1)
      upper[FILTER_MIN_WIDTH] = 1 - Math.max(0, factor * (1 - lower[FILTER_MIN_WIDTH]));
    }
    if (upper[FILTER_MIN_WIDTH] != 0) {
      // (assuming upper is more than 1)
      lower[FILTER_MAX_WIDTH] = 1 + Math.max(0, factor * (upper[FILTER_MAX_WIDTH] - 1));
    }

    // Round the ranges
    final double[] interval = new double[stats[0].length];
    interval[FILTER_SIGNAL] = SignalFilter.DEFAULT_INCREMENT;
    interval[FILTER_SNR] = SnrFilter.DEFAULT_INCREMENT;
    interval[FILTER_MIN_WIDTH] = WidthFilter2.DEFAULT_MIN_INCREMENT;
    interval[FILTER_MAX_WIDTH] = WidthFilter.DEFAULT_INCREMENT;
    interval[FILTER_SHIFT] = ShiftFilter.DEFAULT_INCREMENT;
    interval[FILTER_ESHIFT] = EShiftFilter.DEFAULT_INCREMENT;
    interval[FILTER_PRECISION] = PrecisionFilter.DEFAULT_INCREMENT;
    interval[FILTER_ITERATIONS] = 0.1;
    interval[FILTER_EVALUATIONS] = 0.1;

    // Create a range increment
    final double[] increment = new double[lower.length];
    for (int i = 0; i < increment.length; i++) {
      lower[i] = MathUtils.floor(lower[i], interval[i]);
      upper[i] = MathUtils.ceil(upper[i], interval[i]);
      final double range = upper[i] - lower[i];
      // Allow clipping if the range is small compared to the min increment
      double multiples = range / interval[i];
      // Use 8 multiples for the equivalent of +/- 4 steps around the centre
      if (multiples < 8) {
        multiples = Math.ceil(multiples);
      } else {
        multiples = 8;
      }
      increment[i] = MathUtils.ceil(range / multiples, interval[i]);

      if (i == FILTER_MIN_WIDTH) {
        // Requires clipping based on the upper limit
        lower[i] = upper[i] - increment[i] * multiples;
      } else {
        upper[i] = lower[i] + increment[i] * multiples;
      }
    }

    for (int i = 0; i < stats[0].length; i++) {
      lower[i] = MathUtils.round(lower[i]);
      upper[i] = MathUtils.round(upper[i]);
      min[i] = MathUtils.round(min[i]);
      max[i] = MathUtils.round(max[i]);
      increment[i] = MathUtils.round(increment[i]);
      sb.append('\t').append(min[i]).append(':').append(lower[i]).append('-').append(upper[i])
          .append(':').append(max[i]);
    }

    // Disable some filters
    increment[FILTER_SIGNAL] = Double.POSITIVE_INFINITY;
    // increment[FILTER_SHIFT] = Double.POSITIVE_INFINITY;
    increment[FILTER_ESHIFT] = Double.POSITIVE_INFINITY;

    wo.tile();

    sb.append('\t').append(TextUtils.nanosToString(runTime));

    createTable().append(sb.toString());

    if (settings.saveFilterRange) {
      GUIFilterSettings filterSettings = SettingsManager.readGuiFilterSettings(0);

      String filename = (silent) ? filterSettings.getFilterSetFilename()
          : ImageJUtils.getFilename("Filter_range_file", filterSettings.getFilterSetFilename());
      if (filename == null) {
        return;
      }
      // Remove extension to store the filename
      filename = FileUtils.replaceExtension(filename, ".xml");
      filterSettings = filterSettings.toBuilder().setFilterSetFilename(filename).build();

      // Create a filter set using the ranges
      final ArrayList<Filter> filters = new ArrayList<>(4);
      // Create the multi-filter using the same precision type as that used during fitting.
      // Currently no support for z-filter as 3D astigmatism fitting is experimental.
      final PrecisionMethod precisionMethod =
          getPrecisionMethod((DirectFilter) multiFilter.getFilter());
      Function<double[], Filter> generator;
      if (precisionMethod == PrecisionMethod.POISSON_CRLB) {
        generator = parameters -> new MultiFilterCrlb(parameters[FILTER_SIGNAL],
            (float) parameters[FILTER_SNR], parameters[FILTER_MIN_WIDTH],
            parameters[FILTER_MAX_WIDTH], parameters[FILTER_SHIFT], parameters[FILTER_ESHIFT],
            parameters[FILTER_PRECISION], 0f, 0f);
      } else if (precisionMethod == PrecisionMethod.MORTENSEN) {
        generator = parameters -> new MultiFilter(parameters[FILTER_SIGNAL],
            (float) parameters[FILTER_SNR], parameters[FILTER_MIN_WIDTH],
            parameters[FILTER_MAX_WIDTH], parameters[FILTER_SHIFT], parameters[FILTER_ESHIFT],
            parameters[FILTER_PRECISION], 0f, 0f);
      } else {
        // Default
        generator = parameters -> new MultiFilter2(parameters[FILTER_SIGNAL],
            (float) parameters[FILTER_SNR], parameters[FILTER_MIN_WIDTH],
            parameters[FILTER_MAX_WIDTH], parameters[FILTER_SHIFT], parameters[FILTER_ESHIFT],
            parameters[FILTER_PRECISION], 0f, 0f);
      }
      filters.add(generator.apply(lower));
      filters.add(generator.apply(upper));
      filters.add(generator.apply(increment));
      if (saveFilters(filename, filters)) {
        SettingsManager.writeSettings(filterSettings);
      }

      // Create a filter set using the min/max and the initial bounds.
      // Set sensible limits
      min[FILTER_SIGNAL] = Math.max(min[FILTER_SIGNAL], 30);
      max[FILTER_SNR] = Math.min(max[FILTER_SNR], 10000);
      max[FILTER_PRECISION] = Math.min(max[FILTER_PRECISION], 100);

      // Make the 4-set filters the same as the 3-set filters.
      filters.clear();
      filters.add(generator.apply(min));
      filters.add(generator.apply(lower));
      filters.add(generator.apply(upper));
      filters.add(generator.apply(max));
      saveFilters(FileUtils.replaceExtension(filename, ".4.xml"), filters);
    }

    spotFitResults.min = min;
    spotFitResults.max = max;
  }

  private static FilterCriteria[] createFilterCriteria() {
    FilterCriteria[] criteria = filterCriteria;
    if (criteria == null) {
      criteria = new FilterCriteria[9];
      int index = 0;
      //@formatter:off
      criteria[index++] = new FilterCriteria(ParameterType.SIGNAL,     LowerLimit.MIN, UpperLimit.MAX_POSITIVE_CUMUL_DELTA);
      criteria[index++] = new FilterCriteria(ParameterType.SNR,        LowerLimit.MIN, UpperLimit.MAX_POSITIVE_CUMUL_DELTA);
      criteria[index++] = new FilterCriteria(ParameterType.MIN_WIDTH,  LowerLimit.MIN, UpperLimit.ZERO);
      criteria[index++] = new FilterCriteria(ParameterType.MAX_WIDTH,  LowerLimit.ZERO,        UpperLimit.NINETY_NINE_PERCENT);
      criteria[index++] = new FilterCriteria(ParameterType.SHIFT,      LowerLimit.MAX_NEGATIVE_CUMUL_DELTA, UpperLimit.NINETY_NINE_PERCENT);
      criteria[index++] = new FilterCriteria(ParameterType.ESHIFT,     LowerLimit.MAX_NEGATIVE_CUMUL_DELTA, UpperLimit.NINETY_NINE_PERCENT);
      // Precision has enough discrimination power to be able to use the Jaccard score
      criteria[index++] = new FilterCriteria(ParameterType.PRECISION,  LowerLimit.HALF_MAX_JACCARD_VALUE, UpperLimit.MAX_JACCARD2);
      // These are not filters but are used for stats analysis
      criteria[index++] = new FilterCriteria(null, "Iterations",  LowerLimit.ONE_PERCENT, UpperLimit.NINETY_NINE_NINE_PERCENT, 1, false, false);
      criteria[index]   = new FilterCriteria(null, "Evaluations", LowerLimit.ONE_PERCENT, UpperLimit.NINETY_NINE_NINE_PERCENT, 1, false, false);
      //@formatter:on

      filterCriteria = criteria;
    }
    return filterCriteria;
  }

  private static boolean saveFilters(String filename, ArrayList<Filter> filters) {
    final ArrayList<FilterSet> filterList = new ArrayList<>(1);
    // Add Range keyword to identify as a range filter set
    filterList.add(new FilterSet("Range", filters));
    try (OutputStream out = new BufferedOutputStream(Files.newOutputStream(Paths.get(filename)))) {
      // Use the instance (not .toXML() method) to allow the exception to be caught
      FilterXStreamUtils.getXStreamInstance().toXML(filterList, out);
      return true;
    } catch (final IOException | XStreamException ex) {
      IJ.log("Unable to save the filter set to file: " + ex.getMessage());
    }
    return false;
  }

  private static void printFailures(String title, int[] status) {
    int total = 0;
    // Count failures
    for (int i = 1; i < status.length; i++) {
      if (status[i] != 0) {
        total += status[i];
      }
    }
    // Print failures
    if (total != 0) {
      for (int i = 1; i < status.length; i++) {
        if (status[i] != 0) {
          ImageJUtils.log("%s %s = %d / %d  (%.2f)", title, FitStatus.values()[i].toString(),
              status[i], total, 100.0 * status[i] / total);
        }
      }
    }
    // Print total
    final int all = total + status[0];
    if (all != 0) {
      ImageJUtils.log("%s %s = %d / %d  (%.2f)", title, "Total", total, all, 100.0 * total / all);
    }
  }

  private static void addStatus(int[] status,
      uk.ac.sussex.gdsc.smlm.results.filter.MultiPathFitResult.FitResult fitResult) {
    if (fitResult != null) {
      status[fitResult.status]++;
    }
  }

  private static void addToStats(
      uk.ac.sussex.gdsc.smlm.results.filter.MultiPathFitResult.FitResult fitResult,
      StoredDataStatistics[][] stats) {
    if (fitResult == null) {
      return;
    }

    final FitResult actualFitResult = (FitResult) fitResult.getData();

    if (fitResult.status != 0) {
      // Add the evaluations for spots that were not OK
      stats[0][FILTER_ITERATIONS].add(actualFitResult.getIterations());
      stats[0][FILTER_EVALUATIONS].add(actualFitResult.getEvaluations());
      return;
    }

    if (fitResult.getResults() == null) {
      return;
    }

    boolean isMatch = false;

    for (int resultIndex = 0; resultIndex < fitResult.getResults().length; resultIndex++) {
      final BasePreprocessedPeakResult result =
          (BasePreprocessedPeakResult) fitResult.getResults()[resultIndex];

      // Q. Only build stats on new results?
      if (!result.isNewResult()) {
        continue;
      }

      // This was fit - Get statistics
      final double precision = Math.sqrt(result.getLocationVariance());

      final double signal = result.getSignal();
      final double snr = result.getSnr();
      final double width = result.getXSdFactor();
      final double xShift = result.getXRelativeShift2();
      final double yShift = result.getYRelativeShift2();
      // Since these two are combined for filtering and the max is what matters.
      final double shift = (xShift > yShift) ? Math.sqrt(xShift) : Math.sqrt(yShift);
      final double eshift = Math.sqrt(xShift + yShift);

      stats[0][FILTER_SIGNAL].add(signal);
      stats[0][FILTER_SNR].add(snr);
      if (width < 1) {
        stats[0][FILTER_MIN_WIDTH].add(width);
      } else {
        stats[0][FILTER_MAX_WIDTH].add(width);
      }
      stats[0][FILTER_SHIFT].add(shift);
      stats[0][FILTER_ESHIFT].add(eshift);
      stats[0][FILTER_PRECISION].add(precision);
      if (resultIndex == 0) {
        stats[0][FILTER_ITERATIONS].add(actualFitResult.getIterations());
        stats[0][FILTER_EVALUATIONS].add(actualFitResult.getEvaluations());
      }

      // Add to the TP or FP stats
      // If it has assignments then it was a match to something
      isMatch |= result.hasAssignments();
      final int index = (result.hasAssignments()) ? 1 : 2;
      stats[index][FILTER_SIGNAL].add(signal);
      stats[index][FILTER_SNR].add(snr);
      if (width < 1) {
        stats[index][FILTER_MIN_WIDTH].add(width);
      } else {
        stats[index][FILTER_MAX_WIDTH].add(width);
      }
      stats[index][FILTER_SHIFT].add(shift);
      stats[index][FILTER_ESHIFT].add(eshift);
      stats[index][FILTER_PRECISION].add(precision);
      if (resultIndex == 0) {
        // Nothing else at current
      }
    }

    final int index = (isMatch) ? 1 : 2;
    stats[index][FILTER_ITERATIONS].add(actualFitResult.getIterations());
    stats[index][FILTER_EVALUATIONS].add(actualFitResult.getEvaluations());
  }

  private static boolean noMatch(MultiPathFitResult fitResult) {
    if (isMatch(fitResult.getSingleFitResult())) {
      return false;
    }
    if (isMatch(fitResult.getMultiFitResult())) {
      return false;
    }
    if (isMatch(fitResult.getDoubletFitResult())) {
      return false;
    }
    return !isMatch(fitResult.getMultiDoubletFitResult());
  }

  private static boolean
      isMatch(uk.ac.sussex.gdsc.smlm.results.filter.MultiPathFitResult.FitResult fitResult) {
    if (fitResult == null) {
      return false;
    }
    if (fitResult.status != 0) {
      return false;
    }
    if (fitResult.getResults() == null) {
      return false;
    }
    for (int resultIndex = 0; resultIndex < fitResult.getResults().length; resultIndex++) {
      final BasePreprocessedPeakResult result =
          (BasePreprocessedPeakResult) fitResult.getResults()[resultIndex];
      if (!result.isNewResult()) {
        continue;
      }
      if (result.hasAssignments()) {
        return true;
      }
    }
    return false;
  }

  private static int[] rank(ArrayList<Ranking> pc) {
    Collections.sort(pc, (r1, r2) -> Double.compare(r1.value, r2.value));
    final int[] ranking = new int[pc.size()];
    int rank = 1;
    for (final Ranking r : pc) {
      ranking[r.index] = rank++;
    }
    return ranking;
  }

  private double[] showDoubleHistogram(StoredDataStatistics[][] stats, final int index,
      WindowOrganiser wo, double[][] matchScores) {
    final String xLabel = filterCriteria[index].name;
    LowerLimit lower = filterCriteria[index].lower;
    UpperLimit upper = filterCriteria[index].upper;

    double[] jaccard = null;
    double[] metric = null;
    double maxJaccard = 0;
    if (index <= FILTER_PRECISION
        && (settings.showFilterScoreHistograms || upper.requiresJaccard || lower.requiresJaccard)) {
      // Jaccard score verses the range of the metric
      for (final double[] d : matchScores) {
        if (!Double.isFinite(d[index])) {
          ImageJUtils.log("Error in fit data [%d]: %s", index, d[index]);
        }
      }

      // Do not use Double.compare(double, double) so we get exceptions in the sort for inf/nan
      Arrays.sort(matchScores, (o1, o2) -> {
        if (o1[index] < o2[index]) {
          return -1;
        }
        if (o1[index] > o2[index]) {
          return 1;
        }
        return 0;
      });

      final int scoreIndex = FILTER_PRECISION + 1;
      final int n = results.size();
      double tp = 0;
      double fp = 0;
      jaccard = new double[matchScores.length + 1];
      metric = new double[jaccard.length];
      for (int k = 0; k < matchScores.length; k++) {
        final double score = matchScores[k][scoreIndex];
        tp += score;
        fp += (1 - score);
        jaccard[k + 1] = tp / (fp + n);
        metric[k + 1] = matchScores[k][index];
      }
      metric[0] = metric[1];
      maxJaccard = MathUtils.max(jaccard);

      if (settings.showFilterScoreHistograms) {
        final String title = TITLE + " Jaccard " + xLabel;
        final Plot plot = new Plot(title, xLabel, "Jaccard");
        plot.addPoints(metric, jaccard, Plot.LINE);
        // Remove outliers
        final double[] limitsx = MathUtils.limits(metric);
        final Percentile p = new Percentile();
        final double l = p.evaluate(metric, 25);
        final double u = p.evaluate(metric, 75);
        final double iqr = 1.5 * (u - l);
        limitsx[1] = Math.min(limitsx[1], u + iqr);
        plot.setLimits(limitsx[0], limitsx[1], 0, MathUtils.max(jaccard));
        ImageJUtils.display(title, plot, wo);
      }
    }

    // [0] is all
    // [1] is matches
    // [2] is no match
    final StoredDataStatistics s1 = stats[0][index];
    final StoredDataStatistics s2 = stats[1][index];
    final StoredDataStatistics s3 = stats[2][index];

    if (s1.getN() == 0) {
      return new double[4];
    }

    final DescriptiveStatistics d = s1.getStatistics();
    double median = 0;
    Plot plot = null;
    String title = null;

    if (settings.showFilterScoreHistograms) {
      median = d.getPercentile(50);
      final String label =
          String.format("n = %d. Median = %s nm", s1.getN(), MathUtils.rounded(median));
      final HistogramPlot histogramPlot = new HistogramPlotBuilder(TITLE, s1, xLabel)
          .setMinBinWidth(filterCriteria[index].minBinWidth)
          .setRemoveOutliersOption((filterCriteria[index].restrictRange) ? 1 : 0)
          .setPlotLabel(label).build();
      final PlotWindow plotWindow = histogramPlot.show(wo);
      if (plotWindow == null) {
        IJ.log("Failed to show the histogram: " + xLabel);
        return new double[4];
      }

      title = plotWindow.getTitle();

      // Reverse engineer the histogram settings
      plot = histogramPlot.getPlot();
      final double[] xvalues = histogramPlot.getPlotXValues();
      final int bins = xvalues.length;
      final double yMin = xvalues[0];
      final double binSize = xvalues[1] - xvalues[0];
      final double yMax = xvalues[0] + (bins - 1) * binSize;

      if (s2.getN() > 0) {
        final double[] values = s2.getValues();
        final double[][] hist = HistogramPlot.calcHistogram(values, yMin, yMax, bins);

        if (hist[0].length > 0) {
          plot.setColor(Color.red);
          plot.addPoints(hist[0], hist[1], Plot.BAR);
          ImageJUtils.display(title, plot);
        }
      }

      if (s3.getN() > 0) {
        final double[] values = s3.getValues();
        final double[][] hist = HistogramPlot.calcHistogram(values, yMin, yMax, bins);

        if (hist[0].length > 0) {
          plot.setColor(Color.blue);
          plot.addPoints(hist[0], hist[1], Plot.BAR);
          ImageJUtils.display(title, plot);
        }
      }
    }

    // Do cumulative histogram
    final double[][] h1 = MathUtils.cumulativeHistogram(s1.getValues(), true);
    final double[][] h2 = MathUtils.cumulativeHistogram(s2.getValues(), true);
    final double[][] h3 = MathUtils.cumulativeHistogram(s3.getValues(), true);

    if (settings.showFilterScoreHistograms) {
      title = TITLE + " Cumul " + xLabel;
      plot = new Plot(title, xLabel, "Frequency");
      // Find limits
      double[] xlimit = MathUtils.limits(h1[0]);
      xlimit = MathUtils.limits(xlimit, h2[0]);
      xlimit = MathUtils.limits(xlimit, h3[0]);
      // Restrict using the inter-quartile range
      if (filterCriteria[index].restrictRange) {
        final double q1 = d.getPercentile(25);
        final double q2 = d.getPercentile(75);
        final double iqr = (q2 - q1) * 2.5;
        xlimit[0] = MathUtils.max(xlimit[0], median - iqr);
        xlimit[1] = MathUtils.min(xlimit[1], median + iqr);
      }
      plot.setLimits(xlimit[0], xlimit[1], 0, 1.05);
      plot.addPoints(h1[0], h1[1], Plot.LINE);
      plot.setColor(Color.red);
      plot.addPoints(h2[0], h2[1], Plot.LINE);
      plot.setColor(Color.blue);
      plot.addPoints(h3[0], h3[1], Plot.LINE);
    }

    // Determine the maximum difference between the TP and FP
    double maxx1 = 0;
    double maxx2 = 0;
    double max1 = 0;
    double max2 = 0;

    // We cannot compute the delta histogram, or use percentiles
    if (s2.getN() == 0) {
      upper = UpperLimit.ZERO;
      lower = LowerLimit.ZERO;
    }

    final boolean requireLabel =
        (settings.showFilterScoreHistograms && filterCriteria[index].requireLabel);
    if (requireLabel || upper.requiresDeltaHistogram() || lower.requiresDeltaHistogram()) {
      if (s2.getN() != 0 && s3.getN() != 0) {
        final LinearInterpolator li = new LinearInterpolator();
        final PolynomialSplineFunction f1 = li.interpolate(h2[0], h2[1]);
        final PolynomialSplineFunction f2 = li.interpolate(h3[0], h3[1]);
        for (final double x : h1[0]) {
          if (x < h2[0][0] || x < h3[0][0]) {
            continue;
          }
          try {
            final double v1 = f1.value(x);
            final double v2 = f2.value(x);
            final double diff = v2 - v1;
            if (diff > 0) {
              if (max1 < diff) {
                max1 = diff;
                maxx1 = x;
              }
            } else if (max2 > diff) {
              max2 = diff;
              maxx2 = x;
            }
          } catch (final OutOfRangeException ex) {
            // Because we reached the end
            break;
          }
        }
      }
    }

    if (plot != null) {
      // We use bins=1 on charts where we do not need a label
      if (requireLabel) {
        final String label = String.format("Max+ %s @ %s, Max- %s @ %s", MathUtils.rounded(max1),
            MathUtils.rounded(maxx1), MathUtils.rounded(max2), MathUtils.rounded(maxx2));
        plot.setColor(Color.black);
        plot.addLabel(0, 0, label);
      }
      ImageJUtils.display(title, plot, wo);
    }

    // Now compute the bounds using the desired limit
    double lowerBound;
    double upperBound;
    switch (lower) {
      case MAX_NEGATIVE_CUMUL_DELTA:
        // Switch to percentiles if we have no delta histogram
        if (maxx2 < 0) {
          lowerBound = maxx2;
          break;
        }
        // fall-through
      case ONE_PERCENT:
        lowerBound = getPercentile(h2, 0.01);
        break;
      case MIN:
        lowerBound = getPercentile(h2, 0.0);
        break;
      case ZERO:
        lowerBound = 0;
        break;
      case HALF_MAX_JACCARD_VALUE:
        lowerBound = getXValue(metric, jaccard, maxJaccard * 0.5);
        break;
      default:
        throw new IllegalStateException("Missing lower limit method");
    }
    switch (upper) {
      case MAX_POSITIVE_CUMUL_DELTA:
        // Switch to percentiles if we have no delta histogram
        if (maxx1 > 0) {
          upperBound = maxx1;
          break;
        }
        // fall-through
      case NINETY_NINE_PERCENT:
        upperBound = getPercentile(h2, 0.99);
        break;
      case NINETY_NINE_NINE_PERCENT:
        upperBound = getPercentile(h2, 0.999);
        break;
      case ZERO:
        upperBound = 0;
        break;
      case MAX_JACCARD2:
        upperBound = getXValue(metric, jaccard, maxJaccard) * 2;
        // System.out.printf("MaxJ = %.4f @ %.3f\n", maxJ, u / 2);
        break;
      default:
        throw new IllegalStateException("Missing upper limit method");
    }
    final double min = getPercentile(h1, 0);
    final double max = getPercentile(h1, 1);
    return new double[] {lowerBound, upperBound, min, max};
  }

  /**
   * Gets the percentile.
   *
   * @param histgram The cumulative histogram
   * @param fraction The fraction
   * @return The value for the given fraction
   */
  private static double getPercentile(double[][] histgram, double fraction) {
    final double[] x = histgram[0];
    final double[] y = histgram[1];
    return getXValue(x, y, fraction);
  }

  /**
   * Gets the value from x corresponding to the value in the y values.
   *
   * @param x the x
   * @param y the y
   * @param yvalue the y-value
   * @return the x-value
   */
  private static double getXValue(double[] x, double[] y, double yvalue) {
    for (int i = 0; i < x.length; i++) {
      if (y[i] >= yvalue) {
        if (i == 0 || y[i] == yvalue) {
          return x[i];
        }
        // Interpolation
        return MathUtils.interpolateX(x[i - 1], y[i - 1], x[i], y[i], yvalue);
      }
    }
    return x[x.length - 1];
  }

  private static void add(StringBuilder sb, String value) {
    sb.append('\t').append(value);
  }

  private static void add(StringBuilder sb, int value) {
    sb.append('\t').append(value);
  }

  private static void add(StringBuilder sb, double value) {
    add(sb, MathUtils.rounded(value));
  }

  private static void addCount(StringBuilder sb, double value) {
    // Check if the double holds an integer count
    if ((int) value == value) {
      sb.append('\t').append((int) value);
      // Otherwise add the counts using at least 2 decimal places
    } else if (value > 100) {
      sb.append('\t').append(IJ.d2s(value));
    } else {
      add(sb, MathUtils.rounded(value));
    }
  }

  private static TextWindow createTable() {
    return ImageJUtils.refresh(summaryTableRef,
        () -> new TextWindow(TITLE, createHeader(), "", 1000, 300));
  }

  private static String createHeader() {
    final StringBuilder sb = new StringBuilder(
        "Frames\tW\tH\tMolecules\tDensity (um^-2)\tN\ts (nm)\ta (nm)\tDepth (nm)\t"
            + "Fixed\tGain\tReadNoise (ADUs)\tB (photons)\tNoise (photons)\tSNR\ts (px)\t");
    sb.append("Filter\t");
    sb.append("Spots\t");
    sb.append("nP\t");
    sb.append("nN\t");
    sb.append("fP\t");
    sb.append("fN\t");
    sb.append("Solver\t");
    sb.append("Fitting");

    tablePrefix = sb.toString();

    sb.append('\t');
    sb.append("% nP\t");
    sb.append("% nN\t");
    sb.append("Total\t");
    sb.append("cTP\t");
    sb.append("cFP\t");

    sb.append("cRecall\t");
    sb.append("cPrecision\t");
    sb.append("cF1\t");
    sb.append("cJaccard\t");

    sb.append("Fail cTP\t");
    sb.append("Fail cFP\t");
    sb.append("TP\t");
    sb.append("FP\t");
    sb.append("Recall\t");
    sb.append("Precision\t");
    sb.append("F1\t");
    sb.append("Jaccard\t");

    sb.append("pF1\t");
    sb.append("pJaccard\t");

    sb.append("Med.Distance (nm)\t");
    sb.append("Med.Depth (nm)\t");
    sb.append("Correlation\t");
    sb.append("Ranked\t");
    sb.append("Slope\t");

    createFilterCriteria();
    for (final FilterCriteria f : filterCriteria) {
      sb.append(f.name).append('\t');
    }

    sb.append("Run time");
    return sb.toString();
  }

  private double getSa() {
    return PsfCalculator.squarePixelAdjustment(simulationParameters.sd,
        simulationParameters.pixelPitch);
  }

  /**
   * Updates the given configuration using the latest settings used in benchmarking.
   *
   * @param targetConfiguration the configuration
   * @return true, if successful
   */
  public static boolean updateConfiguration(FitEngineConfiguration targetConfiguration) {
    final FitConfiguration targetFitConfiguration = targetConfiguration.getFitConfiguration();

    // Q. Why are the settings set to themselves?
    // Removed this for now.
    // targetFitConfiguration.setPsf(targetFitConfiguration.getPsf());
    // targetFitConfiguration.setFitSolverSettings(targetFitConfiguration.getFitSolverSettings());
    // targetFitConfiguration.setFilterSettings(targetFitConfiguration.getFilterSettings());

    final FitEngineConfiguration sourceConfig = Settings.lastSettings.get().config;
    final FitConfiguration sourceFitConfiguration = sourceConfig.getFitConfiguration();
    // Set the fit engine settings manually to avoid merging all child settings
    // i.e. do not do a global update using:
    // targetConfiguration.setFitEngineSettings(sourceConfig.getFitEngineSettings());
    targetFitConfiguration.setPsf(sourceFitConfiguration.getPsf());
    targetFitConfiguration.setFitSolverSettings(sourceFitConfiguration.getFitSolverSettings());
    targetFitConfiguration.setFilterSettings(sourceFitConfiguration.getFilterSettings());
    targetConfiguration.setFitting(sourceConfig.getFitting());
    targetConfiguration.setIncludeNeighbours(sourceConfig.isIncludeNeighbours());
    targetConfiguration.setNeighbourHeightThreshold(sourceConfig.getNeighbourHeightThreshold());
    targetConfiguration.setDuplicateDistance(sourceConfig.getDuplicateDistance());
    targetConfiguration.setDuplicateDistanceAbsolute(sourceConfig.getDuplicateDistanceAbsolute());

    if (getComputeDoublets()) {
      targetConfiguration.setResidualsThreshold(0);
      targetFitConfiguration.setComputeResiduals(true);
    } else {
      targetConfiguration.setResidualsThreshold(1);
      targetFitConfiguration.setComputeResiduals(false);
    }

    // We used simple filtering.
    targetFitConfiguration.setSmartFilter(false);

    return true;
  }

  @Override
  public void itemStateChanged(ItemEvent event) {
    if (event.getSource() instanceof Checkbox) {
      final Checkbox checkbox = (Checkbox) event.getSource();

      int failLimit;
      boolean includeNeighbours;
      double neighbourHeightThrehsold;
      boolean computeDoublets;
      MultiPathFilter myMultiFilter;

      if (checkbox.getState()) {
        final FitEngineConfiguration tmp = FitEngineConfiguration.create();
        final FitConfiguration tmpFitConfig = tmp.getFitConfiguration();
        tmpFitConfig.setComputeResiduals(true); // Collect residuals threshold
        if (BenchmarkFilterAnalysis.updateConfiguration(tmp, false)) {
          failLimit = tmp.getFailuresLimit();
          includeNeighbours = tmp.isIncludeNeighbours();
          neighbourHeightThrehsold = tmp.getNeighbourHeightThreshold();
          computeDoublets = tmp.getResidualsThreshold() < 1;

          final DirectFilter primaryFilter = tmpFitConfig.getSmartFilter();
          final double residualsThreshold = tmp.getResidualsThreshold();
          myMultiFilter = new MultiPathFilter(primaryFilter,
              FitWorker.createMinimalFilter(tmpFitConfig.getFilterPrecisionMethod()),
              residualsThreshold);
        } else {
          IJ.log("Failed to update settings using the filter analysis");
          checkbox.setState(false);
          return;
        }
      } else {
        failLimit = config.getFailuresLimit();
        includeNeighbours = config.isIncludeNeighbours();
        neighbourHeightThrehsold = config.getNeighbourHeightThreshold();
        computeDoublets = BenchmarkSpotFit.this.settings.computeDoublets;
        myMultiFilter = multiFilter;
      }

      // Update the dialog
      taFilterXml.setText(myMultiFilter.toXml());
      textFailLimit.setText("" + failLimit);
      cbIncludeNeighbours.setState(includeNeighbours);
      textNeighbourHeight.setText(MathUtils.rounded(neighbourHeightThrehsold));
      cbComputeDoublets.setState(computeDoublets);
    }
  }

  /**
   * Reset the multi path filter and non-filter parameters. This may have been updated when copying
   * benchmark filter settings or by the user within the dialog.
   *
   * @return true, if a reset was required
   */
  boolean resetMultiPathFilter() {
    if (Settings.defaultMultiFilter.equals(multiFilter)
        && equals(createParameters(config), Settings.defaultParameters)) {
      return false;
    }
    multiFilter = Settings.defaultMultiFilter;
    config.setFailuresLimit((int) Settings.defaultParameters[0]);
    // Note we are not resetting the residuals threshold in the config.
    // Only the threshold in the multi-filter matters.
    config.setDuplicateDistance(Settings.defaultParameters[1]);
    config.setDuplicateDistanceAbsolute(Settings.defaultParameters[2] != 0);
    return true;
  }

  private static boolean equals(double[] currentParameters, double[] previousParameters) {
    for (int i = 0; i < previousParameters.length; i++) {
      if (previousParameters[i] != currentParameters[i]) {
        return false;
      }
    }
    return true;
  }

  private static double[] createParameters(FitEngineConfiguration config) {
    return new double[] {config.getFailuresLimit(), config.getDuplicateDistance(),
        (config.getDuplicateDistanceAbsolute()) ? 1 : 0};
  }

  /**
   * Gets a copy of the fit engine configuration from the last execution of the plugin.
   *
   * @return the fit engine configuration
   */
  static FitEngineConfiguration getFitEngineConfiguration() {
    return Settings.lastSettings.get().config.createCopy();
  }

  /**
   * Gets a copy of the multi filter from the last execution of the plugin.
   *
   * @return the multi filter
   */
  static MultiPathFilter getMultiFilter() {
    return Settings.lastSettings.get().multiFilter.copy();
  }

  /**
   * Gets the signal factor from the last execution of the plugin.
   *
   * @return the signal factor
   */
  static double getSignalFactor() {
    return Settings.lastSettings.get().signalFactor;
  }

  /**
   * Gets the lower signal factor from the last execution of the plugin.
   *
   * @return the lower signal factor
   */
  static double getLowerSignalFactor() {
    return Settings.lastSettings.get().lowerSignalFactor;
  }

  /**
   * Gets the compute doublets from the last execution of the plugin.
   *
   * @return the compute doublets
   */
  static boolean getComputeDoublets() {
    return Settings.lastSettings.get().computeDoublets;
  }

  /**
   * Gets the min value of the most recent fit data for the given parameter name.
   *
   * @param type the type
   * @return the min
   */
  static double getMin(ParameterType type) {
    final BenchmarkSpotFitResult results = spotFitResults.get();
    return results == null ? 0 : getValue(type, results.min, 0);
  }

  /**
   * Gets the max value of the most recent fit data for the given parameter name.
   *
   * @param type the type
   * @return the max
   */
  static double getMax(ParameterType type) {
    final BenchmarkSpotFitResult results = spotFitResults.get();
    return results == null ? Double.MAX_VALUE : getValue(type, results.max, Double.MAX_VALUE);
  }

  private static double getValue(ParameterType type, double[] array, double defaultValue) {
    if (type == null || array == null || filterCriteria == null) {
      return defaultValue;
    }

    // Assume these are roughly the same
    if (type == ParameterType.PRECISION2 || type == ParameterType.PRECISION_CRLB) {
      type = ParameterType.PRECISION;
    }

    for (int j = 0; j < filterCriteria.length; j++) {
      if (filterCriteria[j].type == type) {
        // Some metrics (e.g. SNR) can be infinite
        return MathUtils.clip(-Double.MAX_VALUE, Double.MAX_VALUE, array[j]);
      }
    }

    // All other types will have a default value
    return defaultValue;
  }

  /**
   * Checks if the plugin finished successfully.
   *
   * @return true if finished
   */
  boolean isFinished() {
    return finished;
  }

  /**
   * Gets the prefix for the results table header. Contains all the standard header data about the
   * input results data.
   *
   * @return the table prexix
   */
  static String getTablePrefix() {
    return tablePrefix;
  }

  /**
   * Gets the benchmark spot fit results.
   *
   * @return the benchmark spot fit results
   */
  static BenchmarkSpotFitResult getBenchmarkSpotFitResults() {
    return spotFitResults.get();
  }

  /**
   * Gets the distance in pixels from the last analysis.
   *
   * @return the distance in pixels
   */
  static double getDistanceInPixels() {
    final BenchmarkSpotFitResult results = spotFitResults.get();
    return results == null ? 0 : results.distanceInPixels;
  }

  /**
   * Gets the lower distance in pixels from the last analysis.
   *
   * @return the lower distance in pixels
   */
  static double getLowerDistanceInPixels() {
    final BenchmarkSpotFitResult results = spotFitResults.get();
    return results == null ? 0 : results.lowerDistanceInPixels;
  }

  /**
   * Gets the precision method from the direct filter, or use the default of Mortensen with a local
   * background.
   *
   * @param filter the filter
   * @return the precision method
   */
  private static PrecisionMethod getPrecisionMethod(DirectFilter filter) {
    final int flags = filter.getValidationFlags();
    if (DirectFilter.areSet(flags, FilterValidationFlag.LOCATION_VARIANCE_CRLB)) {
      return PrecisionMethod.POISSON_CRLB;
    }
    if (DirectFilter.areSet(flags, FilterValidationFlag.LOCATION_VARIANCE)) {
      return PrecisionMethod.MORTENSEN;
    }
    // Default
    return PrecisionMethod.MORTENSEN_LOCAL_BACKGROUND;
  }
}
