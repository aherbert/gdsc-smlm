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

package uk.ac.sussex.gdsc.smlm.ij.plugins.benchmark;

import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.procedure.TIntObjectProcedure;
import gnu.trove.set.hash.TIntHashSet;
import ij.IJ;
import ij.ImageListener;
import ij.ImagePlus;
import ij.Macro;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.ImageWindow;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.plugin.PlugIn;
import ij.plugin.frame.Recorder;
import ij.text.TextWindow;
import java.awt.Checkbox;
import java.awt.Color;
import java.awt.Rectangle;
import java.awt.TextArea;
import java.awt.TextField;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Consumer;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.lang3.concurrent.ConcurrentRuntimeException;
import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.rng.UniformRandomProvider;
import uk.ac.sussex.gdsc.core.annotation.Nullable;
import uk.ac.sussex.gdsc.core.ij.BufferedTextWindow;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.logging.Ticker;
import uk.ac.sussex.gdsc.core.logging.TrackProgress;
import uk.ac.sussex.gdsc.core.match.ClassificationResult;
import uk.ac.sussex.gdsc.core.match.Coordinate;
import uk.ac.sussex.gdsc.core.match.FractionClassificationResult;
import uk.ac.sussex.gdsc.core.match.FractionalAssignment;
import uk.ac.sussex.gdsc.core.match.MatchResult;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.FileUtils;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.RampedScore;
import uk.ac.sussex.gdsc.core.utils.SettingsList;
import uk.ac.sussex.gdsc.core.utils.StoredData;
import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.UnicodeReader;
import uk.ac.sussex.gdsc.core.utils.concurrent.ConcurrencyUtils;
import uk.ac.sussex.gdsc.core.utils.rng.UniformRandomProviders;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.PrecisionMethod;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.GUIFilterSettings;
import uk.ac.sussex.gdsc.smlm.data.config.TemplateProtos.TemplateSettings;
import uk.ac.sussex.gdsc.smlm.engine.FitConfiguration;
import uk.ac.sussex.gdsc.smlm.engine.FitEngineConfiguration;
import uk.ac.sussex.gdsc.smlm.filters.MaximaSpotFilter;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.ga.Chromosome;
import uk.ac.sussex.gdsc.smlm.ga.FitnessFunction;
import uk.ac.sussex.gdsc.smlm.ga.Population;
import uk.ac.sussex.gdsc.smlm.ga.RampedSelectionStrategy;
import uk.ac.sussex.gdsc.smlm.ga.Recombiner;
import uk.ac.sussex.gdsc.smlm.ga.SelectionStrategy;
import uk.ac.sussex.gdsc.smlm.ga.SimpleMutator;
import uk.ac.sussex.gdsc.smlm.ga.SimpleRecombiner;
import uk.ac.sussex.gdsc.smlm.ga.SimpleSelectionStrategy;
import uk.ac.sussex.gdsc.smlm.ga.ToleranceChecker;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ConfigurationTemplate;
import uk.ac.sussex.gdsc.smlm.ij.plugins.HelpUrls;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ParameterUtils;
import uk.ac.sussex.gdsc.smlm.ij.plugins.PeakFit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsMatchCalculator;
import uk.ac.sussex.gdsc.smlm.ij.plugins.SmlmUsageTracker;
import uk.ac.sussex.gdsc.smlm.ij.plugins.SpotFinderPreview;
import uk.ac.sussex.gdsc.smlm.ij.plugins.benchmark.BenchmarkSpotFilter.BenchmarkSpotFilterResult;
import uk.ac.sussex.gdsc.smlm.ij.plugins.benchmark.BenchmarkSpotFit.BenchmarkSpotFitResult;
import uk.ac.sussex.gdsc.smlm.ij.plugins.benchmark.BenchmarkSpotFit.FilterCandidates;
import uk.ac.sussex.gdsc.smlm.ij.results.ResultsImageSampler;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.PeakResultGridManager;
import uk.ac.sussex.gdsc.smlm.results.count.ConsecutiveFailCounter;
import uk.ac.sussex.gdsc.smlm.results.count.Counter;
import uk.ac.sussex.gdsc.smlm.results.count.FrameCounter;
import uk.ac.sussex.gdsc.smlm.results.filter.BasePreprocessedPeakResult;
import uk.ac.sussex.gdsc.smlm.results.filter.CoordinateStore;
import uk.ac.sussex.gdsc.smlm.results.filter.CoordinateStoreFactory;
import uk.ac.sussex.gdsc.smlm.results.filter.DirectFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.Filter;
import uk.ac.sussex.gdsc.smlm.results.filter.FilterScore;
import uk.ac.sussex.gdsc.smlm.results.filter.FilterSet;
import uk.ac.sussex.gdsc.smlm.results.filter.FilterType;
import uk.ac.sussex.gdsc.smlm.results.filter.FilterXStreamUtils;
import uk.ac.sussex.gdsc.smlm.results.filter.GridCoordinateStore;
import uk.ac.sussex.gdsc.smlm.results.filter.IDirectFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.MultiPathFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.MultiPathFilter.FractionScoreStore;
import uk.ac.sussex.gdsc.smlm.results.filter.MultiPathFitResult;
import uk.ac.sussex.gdsc.smlm.results.filter.MultiPathFitResults;
import uk.ac.sussex.gdsc.smlm.results.filter.ParameterType;
import uk.ac.sussex.gdsc.smlm.results.filter.PeakFractionalAssignment;
import uk.ac.sussex.gdsc.smlm.results.filter.PreprocessedPeakResult;
import uk.ac.sussex.gdsc.smlm.results.filter.ResultAssignment;
import uk.ac.sussex.gdsc.smlm.results.filter.ResultAssignmentDistanceComparator;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;
import uk.ac.sussex.gdsc.smlm.search.ConvergenceChecker;
import uk.ac.sussex.gdsc.smlm.search.ConvergenceToleranceChecker;
import uk.ac.sussex.gdsc.smlm.search.FixedDimension;
import uk.ac.sussex.gdsc.smlm.search.FullScoreFunction;
import uk.ac.sussex.gdsc.smlm.search.ScoreFunctionHelper;
import uk.ac.sussex.gdsc.smlm.search.SearchDimension;
import uk.ac.sussex.gdsc.smlm.search.SearchResult;
import uk.ac.sussex.gdsc.smlm.search.SearchSpace;

/**
 * Run different filtering methods on a set of benchmark fitting results outputting performance
 * statistics on the success of the filter. The fitting results are generated by the
 * BenchmarkSpotFit plugin.
 *
 * <p>Filtering is done using e.g. SNR threshold, Precision thresholds, etc. The statistics reported
 * are shown in a table, e.g. precision, Jaccard, F-score.
 */
public class BenchmarkFilterAnalysis
    implements PlugIn, FitnessFunction<FilterScore>, TrackProgress, FullScoreFunction<FilterScore> {
  private static final String TITLE = "Benchmark Filter Analysis";
  private static final UniqueIdPeakResult[] EMPTY = new UniqueIdPeakResult[0];
  private static final int FLAG_OPTIMISE_FILTER = 1;
  private static final int FLAG_OPTIMISE_PARAMS = 2;

  /**
   * This can be used during filtering. However the min filter is not used to determine if
   * candidates are valid (that is the primary filter). It is used to store poor estimates during
   * fitting. So we can set it to null.
   */
  private static final DirectFilter defaultMinimalFilter = null;

  private static AtomicReference<TextWindow> resultsWindowRef = new AtomicReference<>();
  private static AtomicReference<TextWindow> summaryWindowRef = new AtomicReference<>();
  private static AtomicReference<TextWindow> sensitivityWindowRef = new AtomicReference<>();
  private static AtomicReference<TextWindow> gaWindowRef = new AtomicReference<>();
  private static AtomicReference<TextWindow> componentAnalysisWindowRef = new AtomicReference<>();

  /**
   * This map stores a flag to indicate if each dimension of a filter set was optimised. It is used
   * to remember settings for the dialog.
   */
  private static ConcurrentHashMap<Integer, boolean[]> searchRangeMap = new ConcurrentHashMap<>();
  /**
   * This map stores the step size for each dimension of a filter set was optimised. It is used to
   * remember settings for the dialog.
   */
  private static ConcurrentHashMap<Integer, double[]> stepSizeMap = new ConcurrentHashMap<>();

  /** A reference to the latest score filter. */
  private static AtomicReference<DirectFilter> scoreFilterRef = new AtomicReference<>();

  /** A reference to the image and its sampler. */
  private static AtomicReference<Pair<Integer, ResultsImageSampler>> samplerRef =
      new AtomicReference<>();

  /** A reference to the filename, modified timestamp and list of filters in a filter file. */
  private static AtomicReference<Triple<String, Long, List<FilterSet>>> lastFilterList =
      new AtomicReference<>(Triple.of("", 0L, null));

  /** The coordinate cache. This stores the coordinates for a simulation Id. */
  private static AtomicReference<
      Pair<Integer, TIntObjectHashMap<UniqueIdPeakResult[]>>> coordinateCache =
          new AtomicReference<>(Pair.of(-1, null));

  /** A reference to the last prepared fit results. */
  private static AtomicReference<FitResultData> fitResultDataCache =
      new AtomicReference<>(new FitResultData());

  /** A reference to the the best filter results from the latest iteration analysis. */
  private static final AtomicReference<Map<String, ComplexFilterScore>> iterBestFilter =
      new AtomicReference<>();

  /** The actual coordinates from the simulation. */
  private TIntObjectHashMap<UniqueIdPeakResult[]> actualCoordinates;

  /** The prepared fit result data. */
  private FitResultData fitResultData;

  /**
   * Reference to the working results.
   */
  private BenchmarkFilterAnalysisResult filterAnalysisResult;

  /**
   * The first results prefix used in results tables. This is set when reading the dialog settings.
   */
  private String resultsPrefix;
  /**
   * The second results prefix used in results tables. This is set when reading the dialog settings.
   */
  private String resultsPrefix2;
  /**
   * The limit for the fail count used in results tables. This is set when reading the dialog
   * settings.
   */
  private String limitFailCount;

  /**
   * A reference to the results window if not in headless mode. This is used with a
   * BufferedTextWindow to output the large results set.
   */
  private TextWindow resultsWindow;
  /** The genetic analysis summary results output. */
  private Consumer<String> gaWindow;

  private double residualsThreshold = 1; // Disabled
  private double minCriteria;
  private boolean invertCriteria;
  private boolean invertScore;

  private final boolean isHeadless;
  private boolean debug;
  private CoordinateStore coordinateStore;

  private CreateData.SimulationParameters simulationParameters;
  private BenchmarkSpotFilterResult filterResult;
  private MemoryPeakResults results;
  private boolean extraOptions;

  /** The latest results from the BenchmarkSpotFit plugin. */
  private BenchmarkSpotFitResult spotFitResults;
  /** Set to true if the BenchmarkSpotFit plugin computed doublets. */
  private boolean computeDoublets;
  /** The signal factor from the BenchmarkSpotFit plugin. */
  private double fitSignalFactor;

  // Used to tile plot windows
  private final WindowOrganiser wo = new WindowOrganiser();

  private StopWatch filterAnalysisStopWatch;
  private StopWatch parameterAnalysisStopWatch;
  private StopWatch analysisStopWatch;
  private StopWatch iterationStopWatch;

  private int[] uniqueIds;
  private int uniqueIdCount;

  // Used to implement the FitnessFunction interface
  private double limit;
  private String gaStatusPrefix = "";
  private Population<FilterScore> gaPopulation;

  // Used to set the strength on a filter
  private double[] strengthLower;
  private double[] strengthUpper;

  // Used for the scoring of filter sets
  private MultiPathFitResults[] gaResultsList;
  private MultiPathFitResults[] gaResultsListToScore;
  private boolean gaSubset;
  private int gaIteration;
  private DirectFilter searchScoreFilter;
  private FilterScoreResult[] gaScoreResults;
  private int gaScoreIndex;

  private SimpleFilterScore filterScoreOptimum;
  private SimpleParameterScore parameterScoreOptimum;

  private Rectangle bounds;
  private double distanceScallingFactor;

  /** The score filter. */
  private DirectFilter scoreFilter;

  /** The plugin settings. */
  private Settings settings;

  /** Flag used to indicate that the template filename has been set during iterative analysis. */
  private boolean saveTemplateIsSet;

  /**
   * Store the filter candidates data.
   */
  private static class FitResultData {
    SettingsList scoreSettings;
    double duplicateDistance;
    boolean duplicateDistanceAbsolute;

    double distanceInPixels;
    double lowerDistanceInPixels;
    double signalFactor;
    double lowerSignalFactor;

    String resultsPrefix3;
    String limitRange;

    int fittingId;
    MultiPathFitResults[] resultsList;
    int matches;
    int fittedResults;
    int totalResults;
    int notDuplicateCount;
    int newResultCount;
    int maxUniqueId;
    int countActual;
    StoredData depthStats;
    StoredData depthFitStats;
    StoredDataStatistics signalFactorStats;
    StoredDataStatistics distanceStats;

    FitResultData() {
      // Identify this as not based on any fitting results
      fittingId = -1;
    }

    FitResultData(int fittingId, SettingsList settingsList, Settings settings) {
      this.fittingId = fittingId;
      this.scoreSettings = settingsList;
      this.duplicateDistance = settings.duplicateDistance;
      this.duplicateDistanceAbsolute = settings.duplicateDistanceAbsolute;
    }

    public boolean differentSettings(Settings settings) {
      return duplicateDistance != settings.duplicateDistance
          || duplicateDistanceAbsolute != settings.duplicateDistanceAbsolute;
    }
  }

  /**
   * Store the results that are the global results of the plugin. If these are updated then it
   * should be done on a copy and the copy saved.
   */
  private static class BenchmarkFilterAnalysisResult {
    static AtomicReference<BenchmarkFilterAnalysisResult> lastResult =
        new AtomicReference<>(new BenchmarkFilterAnalysisResult());

    SettingsList lastAnalyseSettings;
    SettingsList lastAnalyseParametersSettings;
    List<NamedPlot> plots;
    Map<String, ComplexFilterScore> bestFilter;
    List<FilterResult> scores;

    BenchmarkFilterAnalysisResult() {
      plots = Collections.emptyList();
      bestFilter = Collections.emptyMap();
      scores = Collections.emptyList();
    }

    BenchmarkFilterAnalysisResult(BenchmarkFilterAnalysisResult source) {
      // Copy by reference
      lastAnalyseSettings = source.lastAnalyseSettings;
      lastAnalyseParametersSettings = source.lastAnalyseParametersSettings;
      // Copy
      plots = new ArrayList<>(source.plots);
      bestFilter = new LinkedHashMap<>(source.bestFilter);
      scores = new ArrayList<>(source.scores);
    }

    BenchmarkFilterAnalysisResult copy() {
      return new BenchmarkFilterAnalysisResult(this);
    }

    /**
     * Load a copy of the result.
     *
     * @return the settings
     */
    static BenchmarkFilterAnalysisResult load() {
      return lastResult.get().copy();
    }

    /**
     * Save the result.
     */
    void save() {
      lastResult.set(this);
    }
  }

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    static final String KEY_FILTER_FILENAME = "gdsc.filteranalysis.filterfilename";
    static final String KEY_FILTERSET_FILENAME = "gdsc.filteranalysis.filtersetfilename";
    static final String KEY_TEMPLATE_FILENAME = "gdsc.filteranalysis.templatefilename";

    static final String[] COMPONENT_ANALYSIS_OPTIONS =
        {"None", "Best Ranked", "Ranked", "Best All", "All"};
    static final String[] EVOLVE_OPTIONS =
        {"None", "Genetic Algorithm", "Range Search", "Enrichment Search", "Step Search"};
    static final String[] SEARCH_OPTIONS =
        {"Range Search", "Enrichment Search", "Step Search", "Enumerate"};

    static final String[] COLUMNS = {
        // Scores using integer scoring
        "TP", "FP", "FN", "Precision", "Recall", "F1", "Jaccard",
        // Scores using fractional scoring
        "fTP", "fFP", "fFN", "fPrecision", "fRecall", "fF1", "fJaccard",};

    /** The last settings used by the plugin. This should be updated after plugin execution. */
    static final AtomicReference<Settings> lastSettings = new AtomicReference<>(new Settings());

    int failCount;
    int minFailCount;
    int maxFailCount;

    double residualsThreshold;
    double minResidualsThreshold;
    double maxResidualsThreshold;
    double duplicateDistance;
    // This is a flag that is passed around but is only set once, i.e.
    // analysis is done using either absolute or relative distances.
    boolean duplicateDistanceAbsolute;
    double minDuplicateDistance;
    double maxDuplicateDistance;
    boolean reset;
    boolean showResultsTable;
    boolean showSummaryTable;
    boolean clearTables;
    String filterFilename;
    String filterSetFilename;
    String templateFilename;
    int summaryTopN;
    double summaryDepth;
    int plotTopN;
    boolean saveBestFilter;
    boolean saveTemplate;
    boolean calculateSensitivity;
    double delta;
    int criteriaIndex;
    double criteriaLimit;
    int scoreIndex;
    double upperMatchDistance;
    double partialMatchDistance;
    double upperSignalFactor;
    double partialSignalFactor;
    boolean depthRecallAnalysis;
    boolean scoreAnalysis;
    int componentAnalysis;

    int evolve;
    boolean repeatEvolve;
    int rangeSearchWidth;
    double rangeSearchReduce;
    int maxIterations;
    int refinementMode;
    int enrichmentSamples;
    int seedSize;
    double enrichmentFraction;
    double enrichmentPadding;

    int searchParam;
    boolean repeatSearch;

    // Settings for parameter analysis have a 'pa' prefix.
    // They otherwise have the same name as for the filter analysis.
    int paRangeSearchWidth;
    double paRangeSearchReduce;
    int paMaxIterations;
    int paRefinementMode;
    int paEnrichmentSamples;
    int paSeedSize;
    double paEnrichmentFraction;
    double paEnrichmentPadding;
    int paConvergedCount;

    boolean showTP;
    boolean showFP;
    boolean showFN;

    int populationSize;
    int failureLimit;
    double tolerance;
    int convergedCount;
    double crossoverRate;
    double meanChildren;
    double mutationRate;
    double selectionFraction;
    boolean rampedSelection;
    boolean saveOption;
    double iterationScoreTolerance;
    double iterationFilterTolerance;
    boolean iterationCompareResults;
    double iterationCompareDistance;
    int iterationMaxIterations;
    double iterationMinRangeReduction;
    int iterationMinRangeReductionIteration;
    boolean iterationConvergeBeforeRefit;

    // For the template example
    int countNo;
    int countLow;
    int countHigh;

    boolean reUseFilters;
    boolean expandFilters;

    boolean[] showColumns;
    boolean requireIntegerResults;

    int scoreFailCount;
    double scoreResidualsThreshold;
    double scoreDuplicateDistance;

    String resultsTitle;

    Settings() {
      // Set defaults
      failCount = 1;
      maxFailCount = 10;
      residualsThreshold = 0.3;
      minResidualsThreshold = 0.1;
      maxResidualsThreshold = 0.6;
      duplicateDistanceAbsolute = true;
      maxDuplicateDistance = 5;
      reset = true;
      showSummaryTable = true;
      filterFilename = Prefs.get(KEY_FILTER_FILENAME, "");
      filterSetFilename = Prefs.get(KEY_FILTERSET_FILENAME, "");
      templateFilename = Prefs.get(KEY_TEMPLATE_FILENAME, "");
      if (TextUtils.isNullOrEmpty(templateFilename)) {
        final String currentUsersHomeDir = System.getProperty("user.home");
        templateFilename =
            currentUsersHomeDir + File.separator + "gdsc.smlm" + File.separator + "template";
      }
      summaryDepth = 500;
      delta = 0.1;
      // Use the precision as criteria to ensure a set confidence on results labelled as true
      criteriaIndex = COLUMNS.length - 4;
      criteriaLimit = 0.95;
      // Score using Jaccard
      scoreIndex = COLUMNS.length - 1;
      upperMatchDistance = 100;
      partialMatchDistance = 33;
      upperSignalFactor = 100;
      partialSignalFactor = 50;
      depthRecallAnalysis = true;
      scoreAnalysis = true;
      componentAnalysis = 3;
      rangeSearchWidth = 2;
      rangeSearchReduce = 0.3;
      maxIterations = 30;
      refinementMode = SearchSpace.RefinementMode.SINGLE_DIMENSION.ordinal();
      enrichmentSamples = 5000;
      seedSize = 5000;
      enrichmentFraction = 0.2;
      enrichmentPadding = 0.1;
      searchParam = 3;
      paRangeSearchWidth = 2;
      paRangeSearchReduce = 0.3;
      paMaxIterations = 30;
      paRefinementMode = SearchSpace.RefinementMode.MULTI_DIMENSION.ordinal();
      paEnrichmentSamples = 500;
      paSeedSize = 500;
      paEnrichmentFraction = 0.2;
      paEnrichmentPadding = 0.1;
      paConvergedCount = 2;
      populationSize = 5000;
      failureLimit = 5;
      tolerance = 1e-4;
      convergedCount = 2;
      crossoverRate = 1;
      meanChildren = 2;
      mutationRate = 1;
      selectionFraction = 0.2;
      rampedSelection = true;
      iterationScoreTolerance = 1e-4;
      iterationFilterTolerance = 1e-3;
      iterationCompareDistance = 0.1;
      iterationMaxIterations = 10;
      iterationMinRangeReduction = 0.2;
      iterationMinRangeReductionIteration = 5;
      countNo = 2;
      countLow = 4;
      countHigh = 4;
      reUseFilters = true;
      expandFilters = true;
      showColumns = new boolean[COLUMNS.length];
      Arrays.fill(showColumns, true);
      scoreDuplicateDistance = -1;
      resultsTitle = "";
    }

    Settings(Settings source) {
      failCount = source.failCount;
      minFailCount = source.minFailCount;
      maxFailCount = source.maxFailCount;
      residualsThreshold = source.residualsThreshold;
      minResidualsThreshold = source.minResidualsThreshold;
      maxResidualsThreshold = source.maxResidualsThreshold;
      duplicateDistance = source.duplicateDistance;
      duplicateDistanceAbsolute = source.duplicateDistanceAbsolute;
      minDuplicateDistance = source.minDuplicateDistance;
      maxDuplicateDistance = source.maxDuplicateDistance;
      reset = source.reset;
      showResultsTable = source.showResultsTable;
      showSummaryTable = source.showSummaryTable;
      clearTables = source.clearTables;
      filterFilename = source.filterFilename;
      filterSetFilename = source.filterSetFilename;
      templateFilename = source.templateFilename;
      summaryTopN = source.summaryTopN;
      summaryDepth = source.summaryDepth;
      plotTopN = source.plotTopN;
      saveBestFilter = source.saveBestFilter;
      saveTemplate = source.saveTemplate;
      calculateSensitivity = source.calculateSensitivity;
      delta = source.delta;
      criteriaIndex = source.criteriaIndex;
      criteriaLimit = source.criteriaLimit;
      scoreIndex = source.scoreIndex;
      upperMatchDistance = source.upperMatchDistance;
      partialMatchDistance = source.partialMatchDistance;
      upperSignalFactor = source.upperSignalFactor;
      partialSignalFactor = source.partialSignalFactor;
      depthRecallAnalysis = source.depthRecallAnalysis;
      scoreAnalysis = source.scoreAnalysis;
      componentAnalysis = source.componentAnalysis;
      evolve = source.evolve;
      repeatEvolve = source.repeatEvolve;
      rangeSearchWidth = source.rangeSearchWidth;
      rangeSearchReduce = source.rangeSearchReduce;
      maxIterations = source.maxIterations;
      refinementMode = source.refinementMode;
      enrichmentSamples = source.enrichmentSamples;
      seedSize = source.seedSize;
      enrichmentFraction = source.enrichmentFraction;
      enrichmentPadding = source.enrichmentPadding;
      searchParam = source.searchParam;
      repeatSearch = source.repeatSearch;
      paRangeSearchWidth = source.paRangeSearchWidth;
      paRangeSearchReduce = source.paRangeSearchReduce;
      paMaxIterations = source.paMaxIterations;
      paRefinementMode = source.paRefinementMode;
      paEnrichmentSamples = source.paEnrichmentSamples;
      paSeedSize = source.paSeedSize;
      paEnrichmentFraction = source.paEnrichmentFraction;
      paEnrichmentPadding = source.paEnrichmentPadding;
      paConvergedCount = source.paConvergedCount;
      showTP = source.showTP;
      showFP = source.showFP;
      showFN = source.showFN;
      populationSize = source.populationSize;
      failureLimit = source.failureLimit;
      tolerance = source.tolerance;
      convergedCount = source.convergedCount;
      crossoverRate = source.crossoverRate;
      meanChildren = source.meanChildren;
      mutationRate = source.mutationRate;
      selectionFraction = source.selectionFraction;
      rampedSelection = source.rampedSelection;
      saveOption = source.saveOption;
      iterationScoreTolerance = source.iterationScoreTolerance;
      iterationFilterTolerance = source.iterationFilterTolerance;
      iterationCompareResults = source.iterationCompareResults;
      iterationCompareDistance = source.iterationCompareDistance;
      iterationMaxIterations = source.iterationMaxIterations;
      iterationMinRangeReduction = source.iterationMinRangeReduction;
      iterationMinRangeReductionIteration = source.iterationMinRangeReductionIteration;
      iterationConvergeBeforeRefit = source.iterationConvergeBeforeRefit;
      countNo = source.countNo;
      countLow = source.countLow;
      countHigh = source.countHigh;
      reUseFilters = source.reUseFilters;
      expandFilters = source.expandFilters;
      showColumns = source.showColumns.clone();
      requireIntegerResults = source.requireIntegerResults;
      scoreFailCount = source.scoreFailCount;
      scoreResidualsThreshold = source.scoreResidualsThreshold;
      scoreDuplicateDistance = source.scoreDuplicateDistance;
      resultsTitle = source.resultsTitle;
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
  }

  private static class CustomFractionalAssignment extends PeakFractionalAssignment {
    /** The distance to the other true location in nm. */
    public final double distToTarget;
    public final PeakResult peak;

    public CustomFractionalAssignment(int targetId, int predictedId, double distance, double score,
        double distToTarget, PreprocessedPeakResult spot, PeakResult peak) {
      super(targetId, predictedId, distance, score, spot);
      this.distToTarget = distToTarget;
      this.peak = peak;
    }

    public double getSignalFactor() {
      return BenchmarkSpotFit.computeSignalFactor(peakResult.getSignal(), peak.getIntensity());
    }
  }

  private static class CustomResultAssignment extends ResultAssignment {
    /** The distance to the other true location in nm. */
    public final double distToTarget;
    public final PeakResult peak;

    public CustomResultAssignment(int targetId, double distance, double score, double distToTarget,
        PeakResult peak) {
      super(targetId, distance, score);
      this.distToTarget = distToTarget;
      this.peak = peak;
    }

    @Override
    public FractionalAssignment toFractionalAssignment(int predictedId,
        PreprocessedPeakResult spot) {
      return new CustomFractionalAssignment(targetId, predictedId, distance, score, distToTarget,
          spot, peak);
    }
  }

  private static class UniqueIdPeakResult extends PeakResult {
    private static final long serialVersionUID = 20190319L;
    final int id;
    final int uniqueId;

    public UniqueIdPeakResult(int id, int uniqueId, PeakResult result) {
      super(result.getFrame(), result.getOrigX(), result.getOrigY(), result.getOrigValue(),
          result.getError(), result.getNoise(), result.getMeanIntensity(), result.getParameters(),
          null);
      this.id = id;
      this.uniqueId = uniqueId;
    }

    protected UniqueIdPeakResult(UniqueIdPeakResult source) {
      super(source);
      id = source.id;
      uniqueId = source.uniqueId;
    }

    @Override
    public PeakResult copy() {
      return new UniqueIdPeakResult(this);
    }
  }

  private static class Job {
    final int frame;
    final FilterCandidates candidates;

    Job(int frame, FilterCandidates candidates) {
      this.frame = frame;
      this.candidates = candidates;
    }
  }

  /**
   * Used to allow multi-threading of scoring the fit results.
   */
  private class FitResultsWorker implements Runnable {
    volatile boolean finished;
    final BlockingQueue<Job> jobs;
    final List<MultiPathFitResults> results;
    final double matchDistance;
    final RampedScore distanceScore;
    final RampedScore signalScore;
    final AtomicInteger uniqueId;
    int matches;
    int total;
    int included;
    int includedActual;
    int notDuplicateCount;
    int newResultCount;
    StoredData depthStats;
    StoredData depthFitStats;
    StoredDataStatistics signalFactorStats;
    StoredDataStatistics distanceStats;
    private final boolean checkBorder;
    final float border;
    final float xlimit;
    final float ylimit;
    final CoordinateStore coordinateStore;
    final Ticker ticker;
    final TIntObjectHashMap<UniqueIdPeakResult[]> actualCoordinates;

    FitResultsWorker(BlockingQueue<Job> jobs, List<MultiPathFitResults> syncResults,
        double matchDistance, RampedScore distanceScore, RampedScore signalScore,
        AtomicInteger uniqueId, CoordinateStore coordinateStore, Ticker ticker,
        TIntObjectHashMap<UniqueIdPeakResult[]> actualCoordinates) {
      this.jobs = jobs;
      this.results = syncResults;
      this.matchDistance = matchDistance;
      this.distanceScore = distanceScore;
      this.signalScore = signalScore;
      this.uniqueId = uniqueId;
      this.coordinateStore = coordinateStore;
      this.ticker = ticker;
      this.actualCoordinates = actualCoordinates;

      depthStats = new StoredData();
      depthFitStats = new StoredData();
      signalFactorStats = new StoredDataStatistics();
      distanceStats = new StoredDataStatistics();

      checkBorder = (filterResult.analysisBorder != null && filterResult.analysisBorder.x != 0);
      if (checkBorder) {
        final Rectangle lastAnalysisBorder = filterResult.analysisBorder;
        border = lastAnalysisBorder.x;
        xlimit = lastAnalysisBorder.x + lastAnalysisBorder.width;
        ylimit = lastAnalysisBorder.y + lastAnalysisBorder.height;
      } else {
        border = xlimit = ylimit = 0;
      }
    }

    @Override
    public void run() {
      try {
        for (;;) {
          final Job job = jobs.take();
          if (job == null || job.candidates == null) {
            break;
          }
          if (!finished) {
            // Only run jobs when not finished. This allows the queue to be emptied.
            run(job);
          }
        }
      } catch (final InterruptedException ex) {
        ConcurrencyUtils.interruptAndThrowUncheckedIf(!finished, ex);
      } finally {
        finished = true;
      }
    }

    private void run(Job job) {
      if (ImageJUtils.isInterrupted()) {
        finished = true;
        return;
      }

      final int frame = job.frame;
      final FilterCandidates result = job.candidates;

      depthStats.add(result.zPosition);

      UniqueIdPeakResult[] actual = getCoordinates(actualCoordinates, frame);
      final int nActual = actual.length;
      final boolean[] matched = new boolean[nActual];
      actual = filterUsingBorder(actual);
      // We could use distanceInPixels for the resolution. Using a bigger size may allow the
      // different fit locations to be in the same cell and so the grid manager can use it's cache.
      final double resolution = 2 * getDistanceInPixels();
      final PeakResultGridManager resultGrid = new PeakResultGridManager(actual, resolution);

      final MultiPathFitResult[] multiPathFitResults =
          new MultiPathFitResult[result.fitResult.length];
      int size = 0;

      coordinateStore.clear();

      // TODO - support a multi-pass filter.
      // The results are in order they were fit.
      // For a single pass fitter this will be in order of candidate ranking.
      // For a multi pass fitter this will be in order of candidate ranking, then repeat.
      for (int index = 0; index < multiPathFitResults.length; index++) {
        final MultiPathFitResult fitResult = result.fitResult[index].copy(true);

        // Score the results. Do in order of those likely to be in the same position
        // thus the grid manager can cache the neighbours
        boolean fitted = score(1, fitResult.getSingleFitResult(), resultGrid, matched);
        fitted |= score(2, fitResult.getDoubletFitResult(), resultGrid, matched);
        fitted |= score(3, fitResult.getMultiFitResult(), resultGrid, matched);
        fitted |= score(4, fitResult.getMultiDoubletFitResult(), resultGrid, matched);

        coordinateStore.flush();

        // XXX - comment out while debugging
        if (fitted) {
          multiPathFitResults[size++] = fitResult;
        }
      }

      // Count number of results that had a match
      for (int i = 0; i < matched.length; i++) {
        if (matched[i]) {
          matches++;
        }
      }

      included += size;
      includedActual += actual.length;
      total += multiPathFitResults.length;

      results.add(new MultiPathFitResults(frame, Arrays.copyOf(multiPathFitResults, size),
          result.spots.length, nActual));

      ticker.tick();
    }

    /**
     * Filter using border.
     *
     * @param results the results
     * @return the results that are inside the border
     */
    private UniqueIdPeakResult[] filterUsingBorder(UniqueIdPeakResult[] results) {
      // Respect the border used in BenchmarkSpotFilter?
      if (filterResult.analysisBorder == null || filterResult.analysisBorder.x == 0) {
        return results;
      }
      // Just implement a hard border
      // Note: This is mainly relevant when loading in simulated ground truth data.
      // The simulations performed by Create Data already use a border when doing a
      // random distribution.
      final Rectangle lastAnalysisBorder = filterResult.analysisBorder;
      final float border = lastAnalysisBorder.x;
      final float xlimit = lastAnalysisBorder.x + lastAnalysisBorder.width;
      final float ylimit = lastAnalysisBorder.y + lastAnalysisBorder.height;
      final UniqueIdPeakResult[] results2 = new UniqueIdPeakResult[results.length];
      int count = 0;
      for (int i = 0; i < results.length; i++) {
        final UniqueIdPeakResult p = results[i];
        if (outsideBorder(p, border, xlimit, ylimit)) {
          continue;
        }
        results2[count++] = p;
      }
      if (count < results.length) {
        return Arrays.copyOf(results2, count);
      }
      return results;
    }

    /**
     * Score the new results in the fit result.
     *
     * @param set the set
     * @param fitResult the fit result
     * @param resultGrid the result grid
     * @param matched array of actual results that have been matched
     * @return true, if the fit status was ok
     */
    private boolean score(int set,
        final uk.ac.sussex.gdsc.smlm.results.filter.MultiPathFitResult.FitResult fitResult,
        final PeakResultGridManager resultGrid, final boolean[] matched) {
      if (fitResult != null && fitResult.status == 0) {
        // Get the new results
        for (int i = 0; i < fitResult.getResults().length; i++) {
          final BasePreprocessedPeakResult peak =
              (BasePreprocessedPeakResult) fitResult.getResults()[i];
          peak.setAssignments(null);
          peak.setIgnore(false);

          // Note: We cannot remove bad candidates as we do not know what the minimum filter will
          // be.
          // Instead this is done when we create a subset for scoring.

          if (!peak.isNewResult()) {
            continue;
          }

          peak.uniqueId = uniqueId.incrementAndGet();

          if (checkBorder && outsideBorder(peak, border, xlimit, ylimit)) {
            // Leave the spot in as it is used when picking the results.
            // Flag to ignore from scoring.
            peak.setIgnore(true);
            continue;
          }

          // Flag if it is possible to be a duplicate
          final boolean notDuplicate =
              !coordinateStore.contains(peak.getX(), peak.getY(), peak.getZ());
          if (notDuplicate) {
            notDuplicateCount++;
          }
          newResultCount++;
          peak.setNotDuplicate(notDuplicate);
          coordinateStore.addToQueue(peak.getX(), peak.getY(), peak.getZ());

          // Compare to actual results
          // We do this using the PeakResultGridManager to generate a sublist to score against
          final PeakResult[] actual =
              resultGrid.getPeakResultNeighbours((int) peak.getX(), (int) peak.getY());
          if (actual.length == 0) {
            continue;
          }

          final ArrayList<ResultAssignment> assignments = new ArrayList<>(actual.length);
          for (int j = 0; j < actual.length; j++) {
            final double d2 = actual[j].distance2(peak.getX(), peak.getY());
            if (d2 <= matchDistance) {
              double dist = Math.sqrt(d2);

              // Score ...
              double score = distanceScore.score(dist);
              if (score != 0) {
                final double signalFactor = BenchmarkSpotFit.computeSignalFactor(peak.getSignal(),
                    actual[j].getIntensity());
                if (signalScore != null) {
                  score *= signalScore.score(Math.abs(signalFactor));
                  if (score == 0) {
                    continue;
                  }
                }

                final int id = ((UniqueIdPeakResult) actual[j]).id;

                // Invert for the ranking (i.e. low is best)
                double distance = 1 - score;

                // Ensure a perfect match can still be ranked ... closest first
                if (distance == 0) {
                  distance -= (matchDistance - d2);
                }

                // Store distance in nm
                dist *= simulationParameters.pixelPitch;
                assignments.add(new CustomResultAssignment(id, distance, score, dist, actual[j]));

                // Accumulate for each actual result
                if (!matched[id]) {
                  matched[id] = true;
                  // Depth is stored in the error field
                  depthFitStats.add(actual[j].getZPosition());
                }

                // Accumulate for all possible matches
                signalFactorStats.add(signalFactor);
                distanceStats.add(dist);
              }
            }
          }

          // Save
          if (!assignments.isEmpty()) {
            final ResultAssignment[] tmp = assignments.toArray(new ResultAssignment[0]);
            // Sort here to speed up a later sort of merged assignment arrays
            Arrays.sort(tmp, ResultAssignmentDistanceComparator.INSTANCE);
            peak.setAssignments(tmp);
          }
        }

        return (fitResult.getResults().length != 0);
      }
      return false;
    }

    private UniqueIdPeakResult[] getCoordinates(TIntObjectHashMap<UniqueIdPeakResult[]> coords,
        int time) {
      final UniqueIdPeakResult[] tmp = coords.get(time);
      return (tmp == null) ? EMPTY : tmp;
    }
  }

  // Store the best filter scores
  private static class FilterResult {
    final double score;
    final int failCount;
    final double residualsThreshold;
    final double duplicateDistance;
    final boolean duplicateDistanceAbsolute;
    final ComplexFilterScore filterScore;

    public FilterResult(int failCount, double residualsThreshold, double duplicateDistance,
        boolean duplicateDistanceAbsolute, ComplexFilterScore filterScore) {
      this.score = filterScore.score;
      this.failCount = failCount;
      this.residualsThreshold = residualsThreshold;
      this.duplicateDistance = duplicateDistance;
      this.duplicateDistanceAbsolute = duplicateDistanceAbsolute;
      this.filterScore = filterScore;
    }

    DirectFilter getFilter() {
      return filterScore.getFilter();
    }
  }

  private static enum FilterScoreCompararor implements Comparator<ComplexFilterScore> {
    /** An instance of the comparator. */
    INSTANCE;

    @Override
    public int compare(ComplexFilterScore o1, ComplexFilterScore o2) {
      // Sort by size, smallest first
      final int result = o1.size - o2.size;
      if (result != 0) {
        return result;
      }
      return o1.compareTo(o2);
    }
  }

  private class ComplexFilterScore extends SimpleFilterScore {
    static final char WITHIN = '-';
    static final char BELOW = '<';
    static final char FLOOR = 'L';
    static final char ABOVE = '>';
    static final char CEIL = 'U';

    int index;
    final ClassificationResult result2;
    final int size;
    final int[] combinations;
    final boolean[] enable;
    char[] atLimit;
    String algorithm;
    long time;
    String paramAlgorithm;
    long paramTime;

    private ComplexFilterScore(FilterScoreResult result, boolean allSameType, char[] atLimit,
        int index, long time, ClassificationResult result2, int size, int[] combinations,
        boolean[] enable) {
      super(result, allSameType, result.criteria >= minCriteria);
      this.index = index;
      this.time = time;
      this.result2 = result2;
      this.size = size;
      this.combinations = combinations;
      this.enable = enable;
      this.atLimit = atLimit;
    }

    public ComplexFilterScore(FilterScoreResult result, boolean allSameType, int index, long time,
        ClassificationResult result2, int size, int[] combinations, boolean[] enable) {
      this(result, allSameType, null, index, time, result2, size, combinations, enable);
    }

    public ComplexFilterScore(FilterScoreResult result, char[] atLimit, String algorithm, long time,
        String paramAlgorithm, long paramTime) {
      // This may be used in comparisons of different type so set allSameType to false
      this(result, false, atLimit, 0, time, null, 0, null, null);
      this.algorithm = algorithm;
      this.paramAlgorithm = paramAlgorithm;
      this.paramTime = paramTime;
    }

    public DirectFilter getFilter() {
      return result.filter;
    }

    public char atLimit(int index) {
      return atLimit[index];
    }

    public char[] atLimit() {
      return atLimit;
    }

    String getParamAlgorithm() {
      return (paramAlgorithm == null) ? "" : paramAlgorithm;
    }
  }

  private static class NamedPlot {
    String name;
    String xAxisName;
    double[] xValues;
    double[] yValues;
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

    /**
     * Compare the two results.
     *
     * @param r1 the first result
     * @param r2 the second result
     * @return -1, 0 or 1
     */
    static int compare(NamedPlot r1, NamedPlot r2) {
      return Double.compare(r2.score, r1.score);
    }
  }

  /**
   * Allow the genetic algorithm to be stopped using the escape key.
   */
  private class InterruptChecker extends ToleranceChecker<FilterScore> {
    final int convergedCount;
    int count;

    public InterruptChecker(double relative, double absolute, int convergedCount) {
      super(relative, absolute);
      this.convergedCount = convergedCount;
    }

    @Override
    public boolean converged(Chromosome<FilterScore> previous, Chromosome<FilterScore> current) {
      if (super.converged(previous, current)) {
        count++;
      } else {
        count = 0;
      }
      // Allow no convergence except when escape is pressed
      if (convergedCount >= 0 && count > convergedCount) {
        return true;
      }
      if (IJ.escapePressed()) {
        ImageJUtils.log("STOPPED " + gaStatusPrefix);
        IJ.resetEscape(); // Allow the plugin to continue processing
        return true;
      }
      return false;
    }

    @Override
    protected boolean converged(FilterScore previous, FilterScore current) {
      // Check the score only if both have criteria achieved
      if (current.criteriaPassed) {
        if (previous.criteriaPassed) {
          return converged(previous.score, current.score);
        }
        return false;
      }
      if (previous.criteriaPassed) {
        // This should not happen as current should be better than previous
        return false;
      }

      return converged(previous.criteria, current.criteria);
    }
  }

  /**
   * Allow the range search to be stopped using the escape key.
   */
  private class InterruptConvergenceChecker extends ConvergenceToleranceChecker<FilterScore> {
    /**
     * The number of times it must have already converged before convergence is achieved.
     */
    final int convergedCount;
    int count;

    public InterruptConvergenceChecker(double relative, double absolute, int maxIterations) {
      super(relative, absolute, false, false, maxIterations);
      this.convergedCount = 0;
    }

    /**
     * Instantiates a new interrupt convergence checker.
     *
     * @param relative the relative
     * @param absolute the absolute
     * @param maxIterations the max iterations
     * @param convergedCount The number of times it must have already converged before convergence
     *        is achieved. Set to negative to disable point coordinate comparison.
     */
    public InterruptConvergenceChecker(double relative, double absolute, int maxIterations,
        int convergedCount) {
      super(relative, absolute, false, convergedCount >= 0, maxIterations);
      this.convergedCount = Math.max(0, convergedCount);
    }

    /**
     * Instantiates a new interrupt convergence checker.
     *
     * @param relative the relative
     * @param absolute the absolute
     * @param checkScore the check score
     * @param checkSequence the check sequence
     * @param maxIterations the max iterations
     */
    public InterruptConvergenceChecker(double relative, double absolute, boolean checkScore,
        boolean checkSequence, int maxIterations) {
      super(relative, absolute, checkScore, checkSequence, maxIterations);
      this.convergedCount = 0;
    }

    @Override
    protected void noConvergenceCriteria() {
      // Ignore this as we can stop using interrupt
    }

    @Override
    public boolean converged(SearchResult<FilterScore> previous,
        SearchResult<FilterScore> current) {
      if (super.converged(previous, current)) {
        // Max iterations is a hard limit even if we have a converged count configured
        if (maxIterations != 0 && getIterations() >= maxIterations) {
          return true;
        }

        count++;
      } else {
        count = 0;
      }

      // Note: if the converged count was negative then no convergence can be achieved
      // using the point coordinates as this is disabled. So the code only reaches here with
      // a count of zero, effectively skipping this as it is always false.
      if (count > convergedCount) {
        return true;
      }
      // Stop if interrupted
      if (IJ.escapePressed()) {
        ImageJUtils.log("STOPPED " + gaStatusPrefix);
        IJ.resetEscape(); // Allow the plugin to continue processing
        return true;
      }
      return false;
    }
  }

  /**
   * Configure the convergence for iterative optimisation.
   */
  private class IterationConvergenceChecker {
    InterruptChecker scoreChecker;
    InterruptConvergenceChecker filterChecker;
    TIntObjectHashMap<List<Coordinate>> previousResults;
    boolean canContinue = true;

    public IterationConvergenceChecker(FilterScore current) {
      // We have two different relative thresholds so use 2 convergence checkers,
      // one for the score and one for the filter sequence
      scoreChecker = new InterruptChecker(settings.iterationScoreTolerance,
          settings.iterationScoreTolerance * 1e-3, 0);
      filterChecker = new InterruptConvergenceChecker(settings.iterationFilterTolerance,
          settings.iterationFilterTolerance * 1e-3, false, true, settings.iterationMaxIterations);
      if (settings.iterationCompareResults) {
        previousResults = getResults(current);
      }
    }

    private TIntObjectHashMap<List<Coordinate>> getResults(FilterScore current) {
      return ResultsMatchCalculator
          .getCoordinates(createResults(null, (DirectFilter) current.filter, false));
    }

    public boolean converged(String prefix, FilterScore previous, FilterScore current,
        double[] previousParameters, double[] currentParameters) {
      // Stop if interrupted
      if (IJ.escapePressed()) {
        ImageJUtils.log("STOPPED");
        // Do not reset escape: IJ.resetEscape()
        canContinue = false;
        return true;
      }

      // Must converge on the non-filter parameters
      if (!converged(previousParameters, currentParameters)) {
        if (settings.iterationCompareResults) {
          previousResults = getResults(current);
        }
        return false;
      }

      final SearchResult<FilterScore> p =
          new SearchResult<>(previous.filter.getParameters(), previous);
      final SearchResult<FilterScore> c =
          new SearchResult<>(current.filter.getParameters(), current);

      if (filterChecker.converged(p, c)) {
        logConvergence(prefix, "filter parameters");
        return true;
      }
      // Directly call the method with the scores
      if (scoreChecker.converged(p.getScore(), c.getScore())) {
        logConvergence(prefix, "score");
        return true;
      }

      if (settings.iterationCompareResults) {
        final TIntObjectHashMap<List<Coordinate>> currentResults = getResults(current);
        final MatchResult r = ResultsMatchCalculator.compareCoordinates(currentResults,
            previousResults, settings.iterationCompareDistance);
        if (r.getJaccard() == 1) {
          logConvergence(prefix, "results coordinates");
          return true;
        }
        previousResults = currentResults;
      }

      return false;
    }

    private boolean converged(double[] previousParameters, double[] currentParameters) {
      for (int i = 0; i < previousParameters.length; i++) {
        if (previousParameters[i] != currentParameters[i]) {
          return false;
        }
      }
      return true;
    }

    private void logConvergence(String prefix, String component) {
      if (settings.iterationMaxIterations != 0
          && getIterations() >= settings.iterationMaxIterations) {
        component = "iterations";
        canContinue = false;
      }
      ImageJUtils.log(prefix + " converged on " + component);
    }

    public int getIterations() {
      return filterChecker.getIterations();
    }
  }

  private static class ScoreJob {
    final DirectFilter filter;
    final int index;

    ScoreJob(DirectFilter filter, int index) {
      this.filter = filter;
      this.index = index;
    }
  }

  private static class ParameterScoreJob {
    final double[] point;
    final int index;

    ParameterScoreJob(double[] point, int index) {
      this.point = point;
      this.index = index;
    }
  }

  /**
   * Used to allow multi-threading of the scoring the filters.
   */
  private class ScoreWorker implements Runnable {
    volatile boolean finished;
    final BlockingQueue<ScoreJob> jobs;
    final FilterScoreResult[] scoreResults;
    final boolean createTextResult;
    final DirectFilter minFilter;
    final CoordinateStore coordinateStore;
    final Ticker ticker;

    public ScoreWorker(BlockingQueue<ScoreJob> jobs, FilterScoreResult[] scoreResults,
        boolean createTextResult, CoordinateStore coordinateStore, Ticker ticker) {
      this.jobs = jobs;
      this.scoreResults = scoreResults;
      this.createTextResult = createTextResult;
      this.minFilter =
          (defaultMinimalFilter != null) ? (DirectFilter) defaultMinimalFilter.clone() : null;
      this.coordinateStore = coordinateStore;
      this.ticker = ticker;
    }

    @Override
    public void run() {
      try {
        for (;;) {
          final ScoreJob job = jobs.take();
          if (job == null || job.index == -1) {
            break;
          }
          if (!finished) {
            // Only run jobs when not finished. This allows the queue to be emptied.
            run(job);
          }
        }
      } catch (final InterruptedException ex) {
        ConcurrencyUtils.interruptAndThrowUncheckedIf(!finished, ex);
      } finally {
        finished = true;
      }
    }

    private void run(ScoreJob job) {
      if (ImageJUtils.isInterrupted()) {
        finished = true;
        return;
      }
      // Directly write to the result array, this is thread safe
      scoreResults[job.index] =
          scoreFilter(job.filter, minFilter, createTextResult, coordinateStore);
      ticker.tick();
    }
  }

  /**
   * Used to allow multi-threading of the scoring the filters.
   */
  private class ParameterScoreWorker implements Runnable {
    volatile boolean finished;
    final BlockingQueue<ParameterScoreJob> jobs;
    final ParameterScoreResult[] scoreResults;
    final boolean createTextResult;
    final DirectFilter filter;
    final DirectFilter minFilter;
    final GridCoordinateStore gridCoordinateStore;
    final Ticker ticker;

    public ParameterScoreWorker(BlockingQueue<ParameterScoreJob> jobs,
        ParameterScoreResult[] scoreResults, boolean createTextResult, Ticker ticker) {
      this.jobs = jobs;
      this.scoreResults = scoreResults;
      this.createTextResult = createTextResult;
      this.filter = (DirectFilter) searchScoreFilter.clone();
      this.minFilter =
          (defaultMinimalFilter != null) ? (DirectFilter) defaultMinimalFilter.clone() : null;
      getBounds();
      this.gridCoordinateStore =
          new GridCoordinateStore(bounds.x, bounds.y, bounds.width, bounds.height, 0, 0);
      this.ticker = ticker;
    }

    @Override
    public void run() {
      try {
        for (;;) {
          final ParameterScoreJob job = jobs.take();
          if (job == null || job.index == -1) {
            break;
          }
          if (!finished) {
            // Only run jobs when not finished. This allows the queue to be emptied.
            run(job);
          }
        }
      } catch (final InterruptedException ex) {
        ConcurrencyUtils.interruptAndThrowUncheckedIf(!finished, ex);
      } finally {
        finished = true;
      }
    }

    private void run(ParameterScoreJob job) {
      if (ImageJUtils.isInterrupted()) {
        finished = true;
        return;
      }

      final double[] point = job.point;
      final int failCount = (int) Math.round(point[0]);
      final double localResidualsThreshold = point[1];
      final double duplicateDistance = point[2];
      CoordinateStore coordinateStore2;

      // Re-use
      gridCoordinateStore.changeXyResolution(duplicateDistance * distanceScallingFactor);
      coordinateStore2 = gridCoordinateStore;

      // New
      // coordinateStore2 = new GridCoordinateStore(bounds[0], bounds[1], duplicateDistance *
      // distanceScallingFactor);

      // From factory
      // coordinateStore2 = createCoordinateStore(duplicateDistance);

      // Directly write to the result array, this is thread safe
      scoreResults[job.index] = scoreFilter(filter, minFilter, failCount, localResidualsThreshold,
          duplicateDistance, coordinateStore2, createTextResult);

      // Allow debugging the score
      // if (failCount == 3 && DoubleEquality.almostEqualRelativeOrAbsolute(residualsThreshold,
      // 0.35, 0, 0.01) &&
      // DoubleEquality.almostEqualRelativeOrAbsolute(duplicateDistance, 2.5, 0, 0.01))
      // {
      // System.out.printf("%s @ %s : %d %f %f \n", Double.toString(scoreResults[job.index].score),
      // Double.toString(scoreResults[job.index].criteria), failCount, residualsThreshold,
      // duplicateDistance);
      // }

      ticker.tick();
    }
  }

  private class ParameterScoreFunction implements FullScoreFunction<FilterScore> {
    @Override
    public SearchResult<FilterScore> findOptimum(double[][] points) {
      gaIteration++;
      SimpleParameterScore max = parameterScoreOptimum;

      // Sort points to allow the CoordinateStore to be reused with the same duplicate distance
      Arrays.sort(points, (o1, o2) -> Double.compare(o1[2], o2[2]));

      final ParameterScoreResult[] scoreResults = scoreFilters(points, settings.showResultsTable);

      if (scoreResults == null) {
        return null;
      }

      for (int index = 0; index < scoreResults.length; index++) {
        final ParameterScoreResult scoreResult = scoreResults[index];
        final SimpleParameterScore result = new SimpleParameterScore(searchScoreFilter, scoreResult,
            scoreResult.criteria >= minCriteria);
        if (result.compareTo(max) < 0) {
          max = result;
        }
      }

      if (settings.showResultsTable) {
        addToResultsWindow(scoreResults);
      }

      parameterScoreOptimum = max;

      // Add the best filter to the table
      // This filter may not have been part of the scored subset so use the entire results set for
      // reporting
      final double[] parameters = max.result.parameters;
      final int failCount = (int) Math.round(parameters[0]);
      final double localResidualsThreshold = parameters[1];
      final double duplicateDistance = parameters[2];

      final MultiPathFilter multiPathFilter =
          new MultiPathFilter(searchScoreFilter, defaultMinimalFilter, localResidualsThreshold);
      final FractionClassificationResult r = multiPathFilter.fractionScoreSubset(
          gaResultsListToScore, createFailCounter(failCount), fitResultData.countActual,
          null, null, createCoordinateStore(duplicateDistance));

      final StringBuilder text = createResult(searchScoreFilter, r,
          buildResultsPrefix2(failCount, localResidualsThreshold, duplicateDistance));
      add(text, gaIteration);
      gaWindow.accept(text.toString());

      return new SearchResult<>(parameters, max);
    }

    private void addToResultsWindow(ParameterScoreResult[] scoreResults) {
      if (resultsWindow != null) {
        try (BufferedTextWindow tw = new BufferedTextWindow(resultsWindow)) {
          tw.setIncrement(0);
          for (int index = 0; index < scoreResults.length; index++) {
            addToResultsOutput(tw::append, scoreResults[index].text);
          }
        }
      } else {
        for (int index = 0; index < scoreResults.length; index++) {
          addToResultsOutput(IJ::log, scoreResults[index].text);
        }
      }
    }

    @Nullable
    @Override
    public SearchResult<FilterScore>[] score(double[][] points) {
      gaIteration++;
      SimpleParameterScore max = parameterScoreOptimum;

      // Sort points to allow the CoordinateStore to be reused with the same duplicate distance
      Arrays.sort(points, (o1, o2) -> Double.compare(o1[2], o2[2]));

      final ParameterScoreResult[] scoreResults = scoreFilters(points, settings.showResultsTable);

      if (scoreResults == null) {
        return null;
      }

      @SuppressWarnings("unchecked")
      final SearchResult<FilterScore>[] scores = new SearchResult[scoreResults.length];
      for (int index = 0; index < scoreResults.length; index++) {
        final ParameterScoreResult scoreResult = scoreResults[index];
        final SimpleParameterScore result = new SimpleParameterScore(searchScoreFilter, scoreResult,
            scoreResult.criteria >= minCriteria);
        if (result.compareTo(max) < 0) {
          max = result;
        }
        scores[index] = new SearchResult<>(points[index], result);
      }

      if (settings.showResultsTable) {
        addToResultsWindow(scoreResults);
      }

      parameterScoreOptimum = max;

      // Add the best filter to the table
      // This filter may not have been part of the scored subset so use the entire results set for
      // reporting
      final double[] parameters = max.result.parameters;
      final int failCount = (int) Math.round(parameters[0]);
      final double residualsThreshold = parameters[1];
      final double duplicateDistance = parameters[2];
      final MultiPathFilter multiPathFilter =
          new MultiPathFilter(searchScoreFilter, defaultMinimalFilter, residualsThreshold);
      final FractionClassificationResult r = multiPathFilter.fractionScoreSubset(
          gaResultsListToScore, createFailCounter(failCount), fitResultData.countActual,
          null, null, createCoordinateStore(duplicateDistance));

      final StringBuilder text = createResult(searchScoreFilter, r,
          buildResultsPrefix2(failCount, residualsThreshold, duplicateDistance));
      add(text, gaIteration);
      gaWindow.accept(text.toString());

      return scores;
    }

    @Override
    public SearchResult<FilterScore>[] cut(SearchResult<FilterScore>[] scores, int size) {
      return ScoreFunctionHelper.cut(scores, size);
    }
  }

  /**
   * Abstract class to allow the array storage to be reused.
   */
  private abstract class CustomTIntObjectProcedure
      implements TIntObjectProcedure<UniqueIdPeakResult[]> {
    float[] x;
    float[] y;
    float[] x2;
    float[] y2;

    CustomTIntObjectProcedure(float[] x, float[] y, float[] x2, float[] y2) {
      this.x = x;
      this.y = y;
      this.x2 = x2;
      this.y2 = y2;
    }
  }

  /**
   * Instantiates a new benchmark filter analysis.
   */
  public BenchmarkFilterAnalysis() {
    isHeadless = java.awt.GraphicsEnvironment.isHeadless();
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    simulationParameters = CreateData.getSimulationParameters();
    if (simulationParameters == null) {
      IJ.error(TITLE, "No benchmark spot parameters in memory");
      return;
    }
    results = CreateData.getResults();
    if (results == null) {
      IJ.error(TITLE, "No benchmark results in memory");
      return;
    }

    settings = Settings.load();
    filterAnalysisResult = BenchmarkFilterAnalysisResult.load();
    // Save now. Any concurrent thread execution will load a copy of the result.
    filterAnalysisResult.save();

    // Iterative optimisation: [Fit; Optimise] x N
    if ("iterate".equals(arg)) {
      iterate();
      return;
    }

    if (invalidBenchmarkSpotFitResults(false)) {
      return;
    }

    extraOptions = ImageJUtils.isExtraOptions();
    // For now do not provide an option, just debug
    debug = true; // extraOptions

    if (!loadFitResults()) {
      return;
    }

    // Score a given filter
    if ("score".equals(arg)) {
      scoreSingleFilter();
      return;
    }
    // Score a given filter
    if ("parameters".equals(arg)) {
      optimiseParameters();
      return;
    }

    optimiseFilter();
  }

  private boolean invalidBenchmarkSpotFitResults(boolean silent) {
    spotFitResults = BenchmarkSpotFit.getBenchmarkSpotFitResults();
    if (spotFitResults == null || spotFitResults.fitResults == null) {
      if (!silent) {
        IJ.error(TITLE, "No benchmark fitting results in memory");
      }
      return true;
    }
    if (spotFitResults.simulationId != simulationParameters.id) {
      if (!silent) {
        IJ.error(TITLE, "Update the benchmark spot fitting for the latest simulation");
      }
      return true;
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
      return true;
    }

    computeDoublets = BenchmarkSpotFit.getComputeDoublets();
    fitSignalFactor = BenchmarkSpotFit.getSignalFactor();

    return false;
  }

  private boolean loadFitResults() {
    if (readResults() == null) {
      IJ.error(TITLE, "No results could be loaded");
      return false;
    }
    return true;
  }

  private void iterate() {
    // If this is run again immediately then provide options for reporting the results
    final Map<String, ComplexFilterScore> localIterBestFilter = iterBestFilter.get();
    // The result best filter will never be a null reference. If it the same reference then
    // the result is from an iteration analysis.
    if (filterAnalysisResult.bestFilter.equals(localIterBestFilter)) {
      final GenericDialog gd = new GenericDialog(TITLE);
      gd.enableYesNoCancel();
      gd.addMessage("Iteration results are held in memory.\n \nReport these results?");
      gd.addHelp(HelpUrls.getUrl("iterate-filter-analysis"));
      gd.showDialog();
      if (gd.wasCanceled()) {
        return;
      }
      if (gd.wasOKed()) {
        reportIterationResults();
        return;
      }
    }

    // Show this dialog first so we can run fully automated after interactive dialogs
    if (!showIterationDialog()) {
      return;
    }

    // Total the time from the interactive plugins
    long time = 0;

    // Run the benchmark fit once interactively, keep the instance
    final BenchmarkSpotFit fit = new BenchmarkSpotFit();
    // Provide ability to skip this step if the fitting has already been done.
    // It must have been done with the default multi-path filter for consistency
    // when re-optimising with different settings (i.e. do not start from a filter
    // than has been optimised before)
    if (fit.resetMultiPathFilter() || invalidBenchmarkSpotFitResults(true)) {
      fit.run(null);
      if (!fit.isFinished()) {
        // The plugin did not complete
        return;
      }
      resetParametersFromFitting();
    }
    if (invalidBenchmarkSpotFitResults(false)) {
      return;
    }
    if (spotFitResults.stopWatch != null) {
      time += spotFitResults.stopWatch.getTime();
    }

    // Run filter analysis once interactively
    if (!loadFitResults()) {
      return;
    }

    // Collect parameters for optimising the parameters
    if (!showDialog(FLAG_OPTIMISE_FILTER | FLAG_OPTIMISE_PARAMS)) {
      return;
    }

    // Load filters from file
    final List<FilterSet> filterSets = readFilterSets();
    if (filterSets == null || filterSets.isEmpty()) {
      IJ.error(TITLE, "No filters specified");
      return;
    }

    ComplexFilterScore current = analyse(filterSets);
    if (current == null) {
      return;
    }
    time += filterAnalysisStopWatch.getTime();

    current = analyseParameters(current);
    if (current == null) {
      return;
    }
    time += parameterAnalysisStopWatch.getTime();

    // Time the non-interactive plugins as a continuous section
    iterationStopWatch = StopWatch.createStarted();

    // Remove the previous iteration results
    iterBestFilter.set(null);

    ImageJUtils.log(TITLE + " Iterating ...");

    final IterationConvergenceChecker checker = new IterationConvergenceChecker(current);

    // Iterate ...
    boolean outerConverged = false;
    int outerIteration = 1;
    double outerRangeReduction = 1;
    while (!outerConverged) {
      if (settings.iterationConvergeBeforeRefit) {
        // Optional inner loop so that the non-filter and filter parameters converge
        // before a refit
        boolean innerConverged = false;
        double innerRangeReduction = 1;
        // -------
        // TODO: Check this. The code does not do anything as the interpolation using x=0
        // so there is no interpolation.
        // int innerIteration = 0;
        // if (settings.iterationMinRangeReduction < 1) {
        // // Linear interpolate down to the min range reduction
        // innerRangeReduction = Math.max(settings.iterationMinRangeReduction,
        // MathUtils.interpolateY(0, 1, settings.iterationMinRangeReductionIteration,
        // settings.iterationMinRangeReduction, innerIteration++));
        // }
        // -------

        // This would make the range too small...
        // innerRangeReduction *= outerRangeReduction
        while (!innerConverged) {
          final ComplexFilterScore previous = current;
          // Re-use the filters as the user may be loading a custom set.
          current = analyse(filterSets, current, innerRangeReduction);
          if (current == null) {
            return;
          }
          final double[] previousParameters = createParameters();
          current = analyseParameters(current, innerRangeReduction);
          if (current == null) {
            return;
          }
          final double[] currentParameters = createParameters();

          innerConverged =
              checker.converged("Filter", previous, current, previousParameters, currentParameters);
        }
        // Check if we can continue (e.g. not max iterations or escape pressed)
        if (!checker.canContinue) {
          break;
        }
      }

      // Do the fit (using the current optimum filter)
      fit.run(current.result.filter, residualsThreshold, settings.failCount,
          settings.duplicateDistance, settings.duplicateDistanceAbsolute);
      if (invalidBenchmarkSpotFitResults(false)) {
        return;
      }
      if (!loadFitResults()) {
        return;
      }

      // Reduce the range over which the filter parameters are searched. Note that this range
      // is centred around the current optimum.
      if (settings.iterationMinRangeReduction < 1) {
        // Linear interpolate down to the min range reduction
        outerRangeReduction = MathUtils.max(settings.iterationMinRangeReduction,
            MathUtils.interpolateY(0, 1, settings.iterationMinRangeReductionIteration,
                settings.iterationMinRangeReduction, outerIteration++));
      }

      // Optimise the filter again.
      final ComplexFilterScore previous = current;
      // Re-use the filters as the user may be loading a custom set.
      current = analyse(filterSets, current, outerRangeReduction);
      if (current == null) {
        break;
      }
      final double[] previousParameters = createParameters();
      current = analyseParameters(current, outerRangeReduction);
      if (current == null) {
        return;
      }
      final double[] currentParameters = createParameters();

      outerConverged =
          checker.converged("Fit+Filter", previous, current, previousParameters, currentParameters);
    }

    if (current != null) {
      // Set-up the plugin so that it can be run again (in iterative mode)
      // and the results reported for the top filter.
      // If the user runs the non-iterative mode then the results will be lost.
      iterBestFilter.set(filterAnalysisResult.bestFilter);
    }

    time += iterationStopWatch.getTime();
    IJ.log("Iteration analysis time : " + TextUtils.millisToString(time));

    showFinished();
  }

  private void resetParametersFromFitting() {
    final FitEngineConfiguration config = BenchmarkSpotFit.getFitEngineConfiguration();
    settings.failCount = config.getFailuresLimit();
    settings.duplicateDistance = config.getDuplicateDistance();
    settings.duplicateDistanceAbsolute = config.getDuplicateDistanceAbsolute();
    residualsThreshold = settings.residualsThreshold = (BenchmarkSpotFit.getComputeDoublets())
        ? BenchmarkSpotFit.getMultiFilter().residualsThreshold
        : 1;
  }

  private double[] createParameters() {
    // Ignore the duplicate distance absolute as this is just for convergence checking
    return new double[] {settings.failCount, residualsThreshold, settings.duplicateDistance};
  }

  private void reportIterationResults() {
    residualsThreshold = settings.residualsThreshold;
    spotFitResults = BenchmarkSpotFit.getBenchmarkSpotFitResults();
    fitResultData = fitResultDataCache.get();
    if (!showReportDialog()) {
      return;
    }
    reportResults(false);
  }

  /**
   * Score a single multi-path filter and report the results. This is used for testing changes to
   * the filter parameters.
   */
  private void scoreSingleFilter() {
    // Show dialog to allow the user to change the settings.
    if (!showScoreDialog()) {
      return;
    }

    // Set the variables used for scoring. Note: these are used throughout the plugin in various
    // methods so it is easier to just set them here and reset them later than pass the score
    // versions through to every method.
    final double[] stash = createParameters();
    settings.failCount = settings.scoreFailCount;
    residualsThreshold = settings.scoreResidualsThreshold;
    settings.duplicateDistance = settings.scoreDuplicateDistance;

    // Create a dummy result, the filter will be rescored in reportResults(...)
    final FilterScoreResult sr = new FilterScoreResult(0, 0, scoreFilter, "");
    final ComplexFilterScore newFilterScore = new ComplexFilterScore(sr, null, "", 0, "", 0);

    // Report to summary window
    reportResults(true, newFilterScore);

    // Reset the variable used for scoring
    settings.failCount = (int) stash[0];
    residualsThreshold = stash[1];
    settings.duplicateDistance = stash[2];
  }

  /**
   * Creates the local score filter from the last saved filter or creates a new one if missing or
   * requested.
   *
   * @param force set to true to force creation of a new filter
   */
  private void createScoreFilter(boolean force) {
    scoreFilter = scoreFilterRef.get();
    if (scoreFilter == null || force) {
      // Reset to default only on first run
      if (settings.scoreDuplicateDistance == -1) {
        settings.scoreFailCount = settings.failCount;
        settings.scoreDuplicateDistance = settings.duplicateDistance;
        settings.scoreResidualsThreshold = residualsThreshold;
      }

      // Use the best result if we have one
      final FilterResult r = getBestResult(filterAnalysisResult.scores);
      if (r != null) {
        scoreFilter = r.getFilter();
        settings.scoreFailCount = r.failCount;
        settings.scoreDuplicateDistance = r.duplicateDistance;
        settings.scoreResidualsThreshold = r.residualsThreshold;
      } else {
        // Default to the fit config settings
        final FitConfiguration tmp = new FitConfiguration();
        // So we get a MultiFilter2
        tmp.setPrecisionMethod(PrecisionMethod.MORTENSEN_LOCAL_BACKGROUND);
        scoreFilter = tmp.getDefaultSmartFilter();
      }

      // Save for reuse
      scoreFilterRef.set(scoreFilter);
    }
  }

  private void optimiseFilter() {
    if (!showDialog(FLAG_OPTIMISE_FILTER)) {
      return;
    }

    // Load filters from file
    final List<FilterSet> filterSets = readFilterSets();

    if (filterSets == null || filterSets.isEmpty()) {
      IJ.error(TITLE, "No filters specified");
      return;
    }

    analyse(filterSets);

    showFinished();
  }

  private void optimiseParameters() {
    final FilterResult fr = getBestResult(filterAnalysisResult.scores);
    if (fr == null) {
      IJ.error(TITLE, "No filter scores in memory");
      return;
    }

    if (!showDialog(FLAG_OPTIMISE_PARAMS)) {
      return;
    }

    analyseParameters(false, fr.filterScore, 0);

    showFinished();
  }

  private static void showFinished() {
    IJ.showStatus("Finished");
  }

  @Nullable
  @SuppressWarnings("unchecked")
  private List<FilterSet> readFilterSets() {
    if (extraOptions) {
      final MultiPathFilter multiFilter = BenchmarkSpotFit.getMultiFilter();
      if (multiFilter != null) {
        final IDirectFilter f = multiFilter.getFilter();
        if (f instanceof DirectFilter) {
          final GenericDialog gd = new GenericDialog(TITLE);
          gd.addMessage("Use an identical filter to " + BenchmarkSpotFit.TITLE);
          gd.enableYesNoCancel();
          gd.hideCancelButton();
          gd.showDialog();
          if (gd.wasOKed()) {
            final List<FilterSet> filterSets = new ArrayList<>(1);
            final List<Filter> filters = new ArrayList<>(1);
            filters.add((DirectFilter) f);
            final FilterSet filterSet = new FilterSet(filters);
            filterSets.add(filterSet);
            resetParametersFromFitting();
            createResultsPrefix2();
            return filterSets;
          }
        }
      }
    }

    GUIFilterSettings filterSettings = SettingsManager.readGuiFilterSettings(0);

    final String filename =
        ImageJUtils.getFilename("Filter_File", filterSettings.getFilterSetFilename());
    if (filename != null) {
      IJ.showStatus("Reading filters ...");
      filterSettings = filterSettings.toBuilder().setFilterSetFilename(filename).build();

      // Allow the filters to be cached
      final Triple<String, Long, List<FilterSet>> filterCache = lastFilterList.get();
      if (isSameFile(filename, filterCache)) {
        final GenericDialog gd = new GenericDialog(TITLE);
        gd.hideCancelButton();
        gd.addMessage("The same filter file was selected.");
        gd.addCheckbox("Re-use_filters", settings.reUseFilters);
        gd.showDialog();
        if (!gd.wasCanceled()) {
          settings.reUseFilters = gd.getNextBoolean();
          if (settings.reUseFilters) {
            SettingsManager.writeSettings(filterSettings);
            return filterCache.getRight();
          }
        }
      }

      final File file = new File(filename);
      try (BufferedReader input =
          new BufferedReader(new UnicodeReader(new FileInputStream(file), null))) {
        // Use the instance so we can catch the exception
        final Object o = FilterXStreamUtils.getXStreamInstance().fromXML(input);

        if (!(o instanceof List<?>)) {
          IJ.log("No filter sets defined in the specified file: " + filename);
          return null;
        }
        SettingsManager.writeSettings(filterSettings);
        List<FilterSet> filterSets = (List<FilterSet>) o;

        if (containsStandardFilters(filterSets)) {
          IJ.log("Filter sets must contain 'Direct' filters");
          return null;
        }

        // Check they are not empty lists
        final List<FilterSet> filterSets2 = new LinkedList<>();
        for (final FilterSet filterSet : filterSets) {
          if (filterSet.size() != 0) {
            filterSets2.add(filterSet);
          } else {
            IJ.log("Filter set empty: " + filterSet.getName());
          }
        }

        if (filterSets2.isEmpty()) {
          IJ.log("All Filter sets are empty");
          return null;
        }

        // Maintain the same list type
        filterSets.clear();
        filterSets.addAll(filterSets2);

        // Option to enumerate filters
        filterSets = expandFilters(filterSets);

        // Save for re-use
        lastFilterList.set(Triple.of(filename, getLastModified(file), filterSets));

        return filterSets;
      } catch (final Exception ex) {
        IJ.log("Unable to load the filter sets from file: " + ex.getMessage());
      } finally {
        IJ.showStatus("");
      }
    }
    return null;
  }

  private static boolean containsStandardFilters(List<FilterSet> filterSets) {
    for (final FilterSet filterSet : filterSets) {
      for (final Filter f : filterSet.getFilters()) {
        if (f.getFilterType() == FilterType.STANDARD) {
          return true;
        }
      }
    }
    return false;
  }

  /**
   * If filters have been provided in FiltersSets of 3 then expand the filters into a set assuming
   * the three represent min:max:increment.
   *
   * @param filterSets the filter sets
   * @return the expanded filter sets
   */
  private List<FilterSet> expandFilters(List<FilterSet> filterSets) {
    // Note:
    // Do not clear the search range map and step size map when reading a new set of filters.
    // The filters may be the same with slight modifications and so it is useful to keep the last
    // settings.

    final long[] expanded = new long[filterSets.size()];
    final String[] name = new String[expanded.length];
    int count = 0;
    boolean doIt = false;
    for (final FilterSet filterSet : filterSets) {
      if (filterSet.size() == 3 && filterSet.allSameType()) {
        name[count] = filterSet.getName();

        // Check we have min:max:increment by counting the combinations
        final Filter f1 = filterSet.getFilters().get(0);
        final Filter f2 = filterSet.getFilters().get(1);
        final Filter f3 = filterSet.getFilters().get(2);
        final int n = f1.getNumberOfParameters();
        final double[] parameters = new double[n];
        final double[] parameters2 = new double[n];
        final double[] increment = new double[n];
        for (int i = 0; i < n; i++) {
          parameters[i] = f1.getParameterValue(i);
          parameters2[i] = f2.getParameterValue(i);
          increment[i] = f3.getParameterValue(i);
        }
        final long combinations = countCombinations(parameters, parameters2, increment);
        if (combinations > 1) {
          expanded[count] = combinations;
          doIt = true;
        }
      }
      count++;
    }

    if (!doIt) {
      return filterSets;
    }

    final GenericDialog gd = new GenericDialog(TITLE);
    gd.hideCancelButton();
    final StringBuilder sb =
        new StringBuilder("The filter file contains potential triples of min:max:increment.\n \n");
    for (count = 0; count < expanded.length; count++) {
      if (expanded[count] > 0) {
        sb.append("Expand set [").append((count + 1)).append("]");
        if (!TextUtils.isNullOrEmpty(name[count])) {
          sb.append(" ").append(name[count]);
        }
        sb.append(" to ").append(expanded[count]).append(" filters\n");
      }
    }
    gd.addMessage(sb.toString());
    gd.addCheckbox("Expand_filters", settings.expandFilters);
    gd.showDialog();
    if (!gd.wasCanceled()) {
      settings.expandFilters = gd.getNextBoolean();
      if (!settings.expandFilters) {
        return filterSets;
      }
    }

    IJ.showStatus("Expanding filters ...");

    final List<FilterSet> filterList2 = new ArrayList<>(filterSets.size());
    for (final FilterSet filterSet : filterSets) {
      count = filterList2.size();
      if (expanded[count] == 0) {
        filterList2.add(filterSet);
        continue;
      }

      final Filter f1 = filterSet.getFilters().get(0);
      final Filter f2 = filterSet.getFilters().get(1);
      final Filter f3 = filterSet.getFilters().get(2);
      final int n = f1.getNumberOfParameters();
      final double[] parameters = new double[n];
      final double[] parameters2 = new double[n];
      final double[] increment = new double[n];
      for (int i = 0; i < n; i++) {
        parameters[i] = f1.getParameterValue(i);
        parameters2[i] = f2.getParameterValue(i);
        increment[i] = f3.getParameterValue(i);
      }

      final List<Filter> list = expandFilters(f1, parameters, parameters2, increment);

      filterList2.add(new FilterSet(filterSet.getName(), list));
    }

    IJ.showStatus("");

    ImageJUtils.log("Expanded input to %d filters in %s", countFilters(filterList2),
        TextUtils.pleural(filterList2.size(), "set"));
    return filterList2;
  }

  /**
   * Expand filters. Set the increment for any parameter not expanded to zero. Set the input
   * parameters array to the lower bounds. Set the input parameters2 array to the upper bounds.
   *
   * <p>If a parameter is not expanded since the increment is Infinity the parameter is disabled. If
   * it was not expanded for any other reason (increment is zero, NaN values, etc) the weakest
   * parameter of the two input array is used and set as the lower and upper bounds.
   *
   * @param baseFilter the base filter
   * @param parameters the parameters
   * @param parameters2 the parameters 2
   * @param increment the increment
   * @return the list
   */
  private static List<Filter> expandFilters(Filter baseFilter, double[] parameters,
      double[] parameters2, double[] increment) {
    final int n = baseFilter.getNumberOfParameters();
    if (parameters.length < n || parameters2.length < n || increment.length < n) {
      throw new IllegalArgumentException(
          "Input arrays  must be  at least the length of the number of parameters");
    }

    final double[] increment2 = increment.clone();
    final int capacity = (int) countCombinations(parameters, parameters2, increment);

    final ArrayList<Filter> list = new ArrayList<>(capacity);

    // Initialise with a filter set at the minimum for each parameter.
    // Get the weakest parameters for those not expanded.
    Filter f1 = baseFilter.create(parameters);
    final double[] p = parameters2.clone();
    f1.weakestParameters(p);
    for (int i = 0; i < n; i++) {
      if (increment[i] == 0) {
        // Disable if Infinite increment, otherwise use the weakest parameter
        parameters[i] = parameters2[i] =
            (Double.isInfinite(increment2[i])) ? baseFilter.getDisabledParameterValue(i) : p[i];
      }
    }
    f1 = baseFilter.create(parameters);

    list.add(f1);
    for (int i = 0; i < n; i++) {
      if (increment[i] == 0) {
        continue;
      }

      final double min = parameters[i];
      final double max = parameters2[i];
      final double inc = increment[i];
      final double max2 = max + inc;

      // Set the upper bounds for the output and store the expansion params
      final StoredDataStatistics stats = new StoredDataStatistics(10);
      for (double value = min + inc; value < max || value - max < max2 - value; value += inc) {
        parameters2[i] = value;
        stats.add(value);
      }
      final double[] values = stats.getValues();

      final ArrayList<Filter> list2 = new ArrayList<>((values.length + 1) * list.size());
      for (int k = 0; k < list.size(); k++) {
        final Filter f = list.get(k);

        // Copy params of the filter
        for (int j = 0; j < n; j++) {
          p[j] = f.getParameterValue(j);
        }

        for (int l = 0; l < values.length; l++) {
          p[i] = values[l];
          list2.add(f.create(p));
        }
      }
      list.addAll(list2);
    }

    // Sort the filters
    Collections.sort(list);

    return list;
  }

  /**
   * Count combinations. Set the increment for any parameter not expanded to zero.
   *
   * @param parameters the parameters
   * @param parameters2 the parameters 2
   * @param increment the increment
   * @return the combinations
   */
  private static long countCombinations(double[] parameters, double[] parameters2,
      double[] increment) {
    final int n = parameters.length;
    long combinations = 1;
    for (int i = 0; i < n; i++) {
      final double inc = increment[i];
      increment[i] = 0;

      if (!Double.isFinite(inc) || !Double.isFinite(parameters[i])
          || !Double.isFinite(parameters2[i]) || parameters[i] > parameters2[i] || inc < 0) {
        continue;
      }

      increment[i] = inc;

      final double min = parameters[i];
      final double max = parameters2[i];
      final double max2 = max + inc;

      long extra = 1;
      // Check the current value is less than the max or (to avoid small round-off error)
      // that the current value is closer to the max than the next value after the max.
      for (double value = min + inc; value < max || value - max < max2 - value; value += inc) {
        extra++;
      }
      combinations *= extra;
    }
    return combinations;
  }

  /**
   * Checks if is same file.
   *
   * @param filename the filename
   * @param filterCache the filter cache
   * @return true, if is same file
   */
  private static boolean isSameFile(String filename,
      Triple<String, Long, List<FilterSet>> filterCache) {
    if (filterCache.getRight() == null) {
      return false;
    }
    if (filename.equals(filterCache.getLeft())) {
      final File file = new File(filename);
      final long timeStamp = getLastModified(file);
      return timeStamp == filterCache.getMiddle();
    }
    return false;
  }

  /**
   * Gets the last modified time.
   *
   * @param file the file
   * @return the last modified time
   */
  private static long getLastModified(File file) {
    try {
      return file.lastModified();
    } catch (final Exception ex) {
      return 0;
    }
  }

  private MultiPathFitResults[] readResults() {
    // Extract all the results in memory into a list per frame. This can be cached
    boolean update = false;
    Pair<Integer, TIntObjectHashMap<UniqueIdPeakResult[]>> coords = coordinateCache.get();

    if (coords.getKey() != simulationParameters.id) {
      coords = Pair.of(simulationParameters.id, getCoordinates(results));
      coordinateCache.set(coords);
      update = true;
    }

    actualCoordinates = coords.getValue();

    spotFitResults = BenchmarkSpotFit.getBenchmarkSpotFitResults();
    FitResultData localFitResultData = fitResultDataCache.get();

    final SettingsList scoreSettings = new SettingsList(settings.partialMatchDistance,
        settings.upperMatchDistance, settings.partialSignalFactor, settings.upperSignalFactor);
    final boolean equalScoreSettings = scoreSettings.equals(localFitResultData.scoreSettings);

    if (update || localFitResultData.fittingId != spotFitResults.id || !equalScoreSettings
        || localFitResultData.differentSettings(settings)) {
      IJ.showStatus("Reading results ...");

      if (localFitResultData.fittingId < 0) {
        // Copy the settings from the fitter if this is the first run.
        // This just starts the plugin with sensible settings.
        // Q. Should this be per new simulation or fitting result instead?
        final FitEngineConfiguration config = BenchmarkSpotFit.getFitEngineConfiguration();
        settings.failCount = config.getFailuresLimit();
        settings.duplicateDistance = config.getDuplicateDistance();
        settings.duplicateDistanceAbsolute = config.getDuplicateDistanceAbsolute();
        settings.residualsThreshold = (BenchmarkSpotFit.getComputeDoublets())
            ? BenchmarkSpotFit.getMultiFilter().residualsThreshold
            : 1;
      }

      // Only cache results for the same score analysis settings.
      // This functionality is for choosing the optimum filter for the given scoring metric.
      if (!equalScoreSettings) {
        filterAnalysisResult.scores.clear();
      }

      localFitResultData = new FitResultData(spotFitResults.id, scoreSettings, settings);

      // @formatter:off
      // -=-=-=-
      // The scoring is designed to find the best fitter+filter combination for the given spot
      // candidates. The ideal combination would correctly fit+pick all the candidate positions
      // that are close to a localisation.
      //
      // Use the following scoring scheme for all candidates:
      //
      //  Candidates
      // +----------------------------------------+
      // |   Actual matches                       |
      // |  +-----------+                TN       |
      // |  |  FN       |                         |
      // |  |      +----------                    |
      // |  |      | TP |    | Fitted             |
      // |  +-----------+    | spots              |
      // |         |     FP  |                    |
      // |         +---------+                    |
      // +----------------------------------------+
      //
      // Candidates     = All the spot candidates
      // Actual matches = Any spot candidate or fitted spot candidate that matches a localisation
      // Fitted spots   = Any spot candidate that was successfully fitted
      //
      // TP = A spot candidate that was fitted and matches a localisation and is accepted
      // FP = A spot candidate that was fitted but does not match a localisation and is accepted
      // FN = A spot candidate that failed to be fitted but matches a localisation
      //    = A spot candidate that was fitted and matches a localisation and is rejected
      // TN = A spot candidate that failed to be fitted and does not match a localisation
      //    = A spot candidate that was fitted and does not match a localisation and is rejected
      //
      // When fitting only produces one result it is possible to compute the TN score.
      // Since unfitted candidates can only be TN or FN we could accumulate these scores and cache
      // them. This was the old method of benchmarking single spot fitting and allowed more scores
      // to be computed.
      //
      // When fitting produces multiple results then we have to score each fit result against all
      // possible actual results and keep a record of the scores. These can then be assessed when
      // the specific results have been chosen by result filtering.
      //
      // Using a distance ramped scoring function the degree of match can be varied from 0 to 1.
      // Using a signal-factor ramped scoring function the degree of fitted can be varied from 0
      // to 1. When using ramped scoring functions the fractional allocation of scores using the
      // above scheme is performed, i.e. candidates are treated as if they both match and unmatch.
      // This results in an equivalent to multiple analysis using different thresholds and averaging
      // of the scores.
      //
      // The totals TP+FP+TN+FN must equal the number of spot candidates. This allows different
      // fitting methods to be compared since the total number of candidates is the same.
      //
      // Precision = TP / (TP+FP)    : This is always valid as a minimum criteria score
      // Recall    = TP / (TP+FN)    : This is valid between different fitting methods since a
      //                               method that fits more spots will have a potentially lower FN
      // Jaccard   = TP / (TP+FN+FP) : This is valid between fitting methods
      //
      // -=-=-=-
      // As an alternative scoring system, different fitting methods can be compared using the same
      // TP value but calculating FN = localisations - TP and FP as Positives - TP. This creates a
      // score against the original number of simulated molecules using everything that was passed
      // through the filter (Positives). This score is comparable when a different spot candidate
      // filter has been used and the total number of candidates is different, e.g. Mean filtering
      // vs. Gaussian filtering
      // -=-=-=-
      // @formatter:on

      final RampedScore distanceScore =
          new RampedScore(spotFitResults.distanceInPixels * settings.partialMatchDistance / 100.0,
              spotFitResults.distanceInPixels * settings.upperMatchDistance / 100.0);
      localFitResultData.lowerDistanceInPixels = distanceScore.lower;
      localFitResultData.distanceInPixels = distanceScore.upper;
      final double matchDistance = MathUtils.pow2(localFitResultData.distanceInPixels);

      localFitResultData.resultsPrefix3 =
          "\t" + MathUtils.rounded(distanceScore.lower * simulationParameters.pixelPitch) + "\t"
              + MathUtils.rounded(distanceScore.upper * simulationParameters.pixelPitch);
      localFitResultData.limitRange =
          ", d=" + MathUtils.rounded(distanceScore.lower * simulationParameters.pixelPitch) + "-"
              + MathUtils.rounded(distanceScore.upper * simulationParameters.pixelPitch);

      // Signal factor must be greater than 1
      final RampedScore signalScore;
      final double spotSignalFactor = BenchmarkSpotFit.getSignalFactor();
      if (spotSignalFactor > 0 && settings.upperSignalFactor > 0) {
        signalScore = new RampedScore(spotSignalFactor * settings.partialSignalFactor / 100.0,
            spotSignalFactor * settings.upperSignalFactor / 100.0);
        localFitResultData.lowerSignalFactor = signalScore.lower;
        localFitResultData.signalFactor = signalScore.upper;
        localFitResultData.resultsPrefix3 += "\t" + MathUtils.rounded(signalScore.lower) + "\t"
            + MathUtils.rounded(signalScore.upper);
        localFitResultData.limitRange += ", s=" + MathUtils.rounded(signalScore.lower) + "-"
            + MathUtils.rounded(signalScore.upper);
      } else {
        signalScore = null;
        localFitResultData.resultsPrefix3 += "\t0\t0";
        localFitResultData.lowerSignalFactor = localFitResultData.signalFactor = 0;
      }

      // Store all the results
      final ArrayList<MultiPathFitResults> multiPathFitResults =
          new ArrayList<>(spotFitResults.fitResults.size());
      final List<MultiPathFitResults> syncResults =
          Collections.synchronizedList(multiPathFitResults);

      // This could be multi-threaded ...
      final int nThreads = getThreads(spotFitResults.fitResults.size());
      final BlockingQueue<Job> jobs = new ArrayBlockingQueue<>(nThreads * 2);
      final List<FitResultsWorker> workers = new LinkedList<>();
      final List<Thread> threads = new LinkedList<>();
      final AtomicInteger uniqueId = new AtomicInteger();
      final CoordinateStore localCoordinateStore = createCoordinateStore();

      final Ticker ticker =
          ImageJUtils.createTicker(spotFitResults.fitResults.size(), nThreads, null);
      for (int i = 0; i < nThreads; i++) {
        final FitResultsWorker worker =
            new FitResultsWorker(jobs, syncResults, matchDistance, distanceScore, signalScore,
                uniqueId, localCoordinateStore.newInstance(), ticker, actualCoordinates);
        final Thread t = new Thread(worker);
        workers.add(worker);
        threads.add(t);
        t.start();
      }

      spotFitResults.fitResults.forEachEntry((frame, candidates) -> {
        put(jobs, new Job(frame, candidates));
        return true;
      });
      // Finish all the worker threads by passing in a null job
      for (int i = 0; i < threads.size(); i++) {
        put(jobs, new Job(0, null));
      }

      // Wait for all to finish
      for (int i = 0; i < threads.size(); i++) {
        try {
          threads.get(i).join();
          final FitResultsWorker worker = workers.get(i);
          localFitResultData.matches += worker.matches;
          localFitResultData.fittedResults += worker.included;
          localFitResultData.totalResults += worker.total;
          localFitResultData.notDuplicateCount += worker.notDuplicateCount;
          localFitResultData.newResultCount += worker.newResultCount;
          localFitResultData.countActual += worker.includedActual;
          if (i == 0) {
            localFitResultData.depthStats = worker.depthStats;
            localFitResultData.depthFitStats = worker.depthFitStats;
            localFitResultData.signalFactorStats = worker.signalFactorStats;
            localFitResultData.distanceStats = worker.distanceStats;
          } else {
            localFitResultData.depthStats.add(worker.depthStats);
            localFitResultData.depthFitStats.add(worker.depthFitStats);
            localFitResultData.signalFactorStats.add(worker.signalFactorStats);
            localFitResultData.distanceStats.add(worker.distanceStats);
          }
        } catch (final InterruptedException ex) {
          Thread.currentThread().interrupt();
          throw new ConcurrentRuntimeException("Unexpected interrupt", ex);
        }
      }
      threads.clear();
      ImageJUtils.finished();

      localFitResultData.maxUniqueId = uniqueId.get();

      localFitResultData.resultsList = multiPathFitResults.toArray(new MultiPathFitResults[0]);

      Arrays.sort(localFitResultData.resultsList,
          (o1, o2) -> Integer.compare(o1.getFrame(), o2.getFrame()));

      MultiPathFilter.resetValidationFlag(localFitResultData.resultsList);

      fitResultDataCache.set(localFitResultData);
    }
    fitResultData = localFitResultData;

    return localFitResultData.resultsList;
  }

  private static TIntObjectHashMap<UniqueIdPeakResult[]> getCoordinates(MemoryPeakResults results) {
    final TIntObjectHashMap<UniqueIdPeakResult[]> coords = new TIntObjectHashMap<>();
    if (results.size() > 0) {
      // Do not use HashMap directly to build the coords object since there
      // will be many calls to getEntry(). Instead sort the results and use
      // a new list for each time point
      results.sort();

      final Counter uniqueId = new Counter();
      final FrameCounter counter = new FrameCounter();
      final LocalList<PeakResult> tmp = new LocalList<>();
      // Add the results to the lists
      results.forEach((PeakResultProcedure) result -> {
        if (counter.advanceAndReset(result.getFrame()) && !tmp.isEmpty()) {
          coords.put(counter.previousFrame(), tmp.toArray(new UniqueIdPeakResult[0]));
          tmp.clear();
        }
        tmp.add(new UniqueIdPeakResult(tmp.size(), uniqueId.getAndIncrement(), result));
      });

      if (!tmp.isEmpty()) {
        coords.put(counter.currentFrame(), tmp.toArray(new UniqueIdPeakResult[0]));
      }
    }
    return coords;
  }

  private boolean showDialog(int optimiseParameters) {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

    final boolean showOptimiseFilter = (optimiseParameters & FLAG_OPTIMISE_FILTER) != 0;
    final boolean showOptimiseParams = (optimiseParameters & FLAG_OPTIMISE_PARAMS) != 0;

    final String helpKey =
        showOptimiseParams ? "benchmark-filter-parameters" : "benchmark-filter-analysis";

    addSimulationData(gd);

    // TODO - Make minimal filter configurable?

    gd.addSlider("Fail_count", 0, 20, settings.failCount);
    if (showOptimiseParams) {
      gd.addNumericField("Min_fail_count", settings.minFailCount, 0);
      gd.addNumericField("Max_fail_count", settings.maxFailCount, 0);
    }
    if (computeDoublets) {
      gd.addSlider("Residuals_threshold", 0.01, 1, settings.residualsThreshold);
      if (showOptimiseParams) {
        gd.addNumericField("Min_residuals_threshold", settings.minResidualsThreshold, 2);
        gd.addNumericField("Max_residuals_threshold", settings.maxResidualsThreshold, 2);
      }
    }
    final FitEngineConfiguration tmp = new FitEngineConfiguration();
    tmp.setDuplicateDistance(settings.duplicateDistance);
    tmp.setDuplicateDistanceAbsolute(settings.duplicateDistanceAbsolute);
    PeakFit.addDuplicateDistanceOptions(gd, new PeakFit.SimpleFitEngineConfigurationProvider(tmp));
    if (showOptimiseParams) {
      gd.addNumericField("Min_duplicate_distance", settings.minDuplicateDistance, 2);
      gd.addNumericField("Max_duplicate_distance", settings.maxDuplicateDistance, 2);
    }
    gd.addCheckbox("Reset", settings.reset);
    gd.addCheckbox("Show_table", settings.showResultsTable);
    gd.addCheckbox("Show_summary", settings.showSummaryTable);
    gd.addCheckbox("Clear_tables", settings.clearTables);
    gd.addSlider("Summary_top_n", 0, 20, settings.summaryTopN);
    gd.addNumericField("Summary_depth (nm)", settings.summaryDepth, 0);
    gd.addSlider("Plot_top_n", 0, 20, settings.plotTopN);
    gd.addCheckbox("Save_best_filter", settings.saveBestFilter);
    gd.addCheckbox("Save_template", settings.saveTemplate);
    gd.addCheckbox("Calculate_sensitivity", settings.calculateSensitivity);
    gd.addSlider("Delta", 0.01, 1, settings.delta);
    gd.addMessage("Match scoring");
    gd.addChoice("Criteria", Settings.COLUMNS, settings.criteriaIndex);
    gd.addNumericField("Criteria_limit", settings.criteriaLimit, 4);
    gd.addChoice("Score", Settings.COLUMNS, settings.scoreIndex);
    ImageJUtils.addMessage(gd, "Fitting match distance = %s nm; signal factor = %s",
        MathUtils.rounded(spotFitResults.distanceInPixels * simulationParameters.pixelPitch),
        MathUtils.rounded(fitSignalFactor));
    gd.addSlider("Upper_match_distance (%)", 0, 100, settings.upperMatchDistance);
    gd.addSlider("Partial_match_distance (%)", 0, 100, settings.partialMatchDistance);
    gd.addSlider("Upper_signal_factor (%)", 0, 100, settings.upperSignalFactor);
    gd.addSlider("Partial_signal_factor (%)", 0, 100, settings.partialSignalFactor);
    if (!simulationParameters.fixedDepth) {
      gd.addCheckbox("Depth_recall_analysis", settings.depthRecallAnalysis);
    }
    gd.addCheckbox("Score_analysis", settings.scoreAnalysis);
    gd.addChoice("Component_analysis", Settings.COMPONENT_ANALYSIS_OPTIONS,
        settings.componentAnalysis);
    if (showOptimiseFilter) {
      gd.addChoice("Evolve", Settings.EVOLVE_OPTIONS, settings.evolve);
      gd.addCheckbox("Repeat_evolve", settings.repeatEvolve);
    }
    if (showOptimiseParams) {
      gd.addChoice("Search", Settings.SEARCH_OPTIONS, settings.searchParam);
      gd.addCheckbox("Repeat_search", settings.repeatSearch);
    }
    gd.addStringField("Title", settings.resultsTitle, 20);
    final String[] labels = {"Show_TP", "Show_FP", "Show_FN"};
    gd.addCheckboxGroup(1, 3, labels,
        new boolean[] {settings.showTP, settings.showFP, settings.showFN});

    gd.addHelp(HelpUrls.getUrl(helpKey));
    gd.showDialog();

    if (gd.wasCanceled() || !readDialog(gd, optimiseParameters, tmp)) {
      return false;
    }

    if (!selectTableColumns()) {
      return false;
    }

    // We may have to read the results again if the ranking option has changed.
    // Also we must read the results with the maximum duplicate distance we may encounter.
    final double dd = settings.duplicateDistance;
    if (showOptimiseParams) {
      settings.duplicateDistance = settings.maxDuplicateDistance;
    }
    readResults();
    settings.duplicateDistance = dd;

    return true;
  }

  private void addSimulationData(GenericDialog gd) {
    double signal = simulationParameters.minSignal;
    if (simulationParameters.maxSignal > signal) {
      signal += simulationParameters.maxSignal;
      signal *= 0.5;
    }
    final double pSignal =
        CreateData.getPrecisionN(simulationParameters.pixelPitch, simulationParameters.sd, signal,
            MathUtils.pow2(simulationParameters.noise), simulationParameters.isEmCcd());
    final double pLse = Gaussian2DPeakResultHelper.getPrecision(simulationParameters.pixelPitch,
        simulationParameters.sd, signal, simulationParameters.noise,
        simulationParameters.isEmCcd());
    final double pMle = Gaussian2DPeakResultHelper.getMLPrecision(simulationParameters.pixelPitch,
        simulationParameters.sd, signal, simulationParameters.noise,
        simulationParameters.isEmCcd());
    StringBuilder msg = new StringBuilder();
    TextUtils.formatTo(msg,
        "Fit %d/%d results, %d True-Positives, %d unique\nExpected signal = %.3f +/- %.3f\n"
            + "Expected X precision = %.3f (LSE), %.3f (MLE)\nNot duplicates : %d / %d (%.2f%%)",
        fitResultData.fittedResults, fitResultData.totalResults, fitResultData.matches,
        fitResultData.maxUniqueId, signal, pSignal, pLse, pMle, fitResultData.notDuplicateCount,
        fitResultData.newResultCount,
        (100.0 * fitResultData.notDuplicateCount) / fitResultData.newResultCount);

    final FilterResult best = getBestResult(filterAnalysisResult.scores);
    if (best != null) {
      TextUtils.formatTo(msg, "\nCurrent Best=%s, FailCount=%d", MathUtils.rounded(best.score),
          best.failCount);
    }
    gd.addMessage(msg.toString());
  }

  private boolean selectTableColumns() {
    if (settings.showResultsTable || settings.showSummaryTable) {
      final GenericDialog gd = new GenericDialog(TITLE);
      gd.addHelp(HelpUrls.getUrl("benchmark-filter-analysis"));

      gd.addMessage("Select the results:");
      for (int i = 0; i < Settings.COLUMNS.length; i++) {
        gd.addCheckbox(Settings.COLUMNS[i], settings.showColumns[i]);
      }
      gd.showDialog();

      if (gd.wasCanceled()) {
        return false;
      }

      for (int i = 0; i < Settings.COLUMNS.length; i++) {
        settings.showColumns[i] = gd.getNextBoolean();
      }

      settings.requireIntegerResults = false;
      for (int i = 0; i < 7; i++) {
        if (settings.showColumns[i]) {
          settings.requireIntegerResults = true;
          break;
        }
      }
    }
    return true;
  }

  private boolean readDialog(ExtendedGenericDialog gd, int optimiseParameters,
      FitEngineConfiguration tmp) {
    final boolean showOptimiseFilter = (optimiseParameters & FLAG_OPTIMISE_FILTER) != 0;
    final boolean showOptimiseParams = (optimiseParameters & FLAG_OPTIMISE_PARAMS) != 0;

    settings.failCount = (int) Math.abs(gd.getNextNumber());
    if (showOptimiseParams) {
      settings.minFailCount = (int) Math.abs(gd.getNextNumber());
      settings.maxFailCount = (int) Math.abs(gd.getNextNumber());
    }
    if (computeDoublets) {
      // Round to the precision of the min/max
      residualsThreshold =
          settings.residualsThreshold = MathUtils.round(Math.abs(gd.getNextNumber()), 4);
      if (showOptimiseParams) {
        settings.minResidualsThreshold = MathUtils.round(Math.abs(gd.getNextNumber()), 4);
        settings.maxResidualsThreshold = MathUtils.round(Math.abs(gd.getNextNumber()), 4);
      }
    }
    settings.duplicateDistance = MathUtils.round(Math.abs(gd.getNextNumber()), 4);
    if (showOptimiseParams) {
      settings.minDuplicateDistance = MathUtils.round(Math.abs(gd.getNextNumber()), 4);
      settings.maxDuplicateDistance = MathUtils.round(Math.abs(gd.getNextNumber()), 4);
    }
    settings.reset = gd.getNextBoolean();
    settings.showResultsTable = gd.getNextBoolean();
    settings.showSummaryTable = gd.getNextBoolean();
    settings.clearTables = gd.getNextBoolean();
    settings.summaryTopN = (int) Math.abs(gd.getNextNumber());
    settings.summaryDepth = Math.abs(gd.getNextNumber());
    settings.plotTopN = (int) Math.abs(gd.getNextNumber());
    settings.saveBestFilter = gd.getNextBoolean();
    settings.saveTemplate = gd.getNextBoolean();
    settings.calculateSensitivity = gd.getNextBoolean();
    settings.delta = gd.getNextNumber();
    settings.criteriaIndex = gd.getNextChoiceIndex();
    settings.criteriaLimit = gd.getNextNumber();
    settings.scoreIndex = gd.getNextChoiceIndex();
    settings.upperMatchDistance = Math.abs(gd.getNextNumber());
    settings.partialMatchDistance = Math.abs(gd.getNextNumber());
    settings.upperSignalFactor = Math.abs(gd.getNextNumber());
    settings.partialSignalFactor = Math.abs(gd.getNextNumber());
    if (!simulationParameters.fixedDepth) {
      settings.depthRecallAnalysis = gd.getNextBoolean();
    }
    settings.scoreAnalysis = gd.getNextBoolean();
    settings.componentAnalysis = gd.getNextChoiceIndex();
    if (showOptimiseFilter) {
      settings.evolve = gd.getNextChoiceIndex();
      settings.repeatEvolve = gd.getNextBoolean();
    }
    if (showOptimiseParams) {
      settings.searchParam = gd.getNextChoiceIndex();
      settings.repeatSearch = gd.getNextBoolean();
    }
    settings.resultsTitle = gd.getNextString();
    settings.showTP = gd.getNextBoolean();
    settings.showFP = gd.getNextBoolean();
    settings.showFN = gd.getNextBoolean();

    gd.collectOptions();
    settings.duplicateDistanceAbsolute = tmp.getDuplicateDistanceAbsolute();

    resultsPrefix = spotFitResults.resultPrefix + "\t" + settings.resultsTitle + "\t";
    createResultsPrefix2();
    settings.save();

    // Check there is one output
    if (!settings.showResultsTable && !settings.showSummaryTable && !settings.calculateSensitivity
        && settings.plotTopN < 1 && !settings.saveBestFilter) {
      IJ.error(TITLE, "No output selected");
      return false;
    }

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Delta", settings.delta);
      ParameterUtils.isBelow("Delta", settings.delta, 1);
      ParameterUtils.isAboveZero("Upper match distance", settings.upperMatchDistance);
      if (settings.partialMatchDistance > settings.upperMatchDistance) {
        settings.partialMatchDistance = settings.upperMatchDistance;
      }
      if (settings.partialSignalFactor > settings.upperSignalFactor) {
        settings.partialSignalFactor = settings.upperSignalFactor;
      }
      if (showOptimiseParams) {
        ParameterUtils.isEqualOrBelow("Fail count", settings.failCount, settings.maxFailCount);
        ParameterUtils.isEqualOrAbove("Fail count", settings.failCount, settings.minFailCount);
        if (computeDoublets) {
          ParameterUtils.isEqualOrBelow("Residuals threshold", settings.residualsThreshold,
              settings.maxResidualsThreshold);
          ParameterUtils.isEqualOrAbove("Residuals threshold", settings.residualsThreshold,
              settings.minResidualsThreshold);
        }
        ParameterUtils.isEqualOrBelow("Duplicate distance", settings.duplicateDistance,
            settings.maxDuplicateDistance);
        ParameterUtils.isEqualOrAbove("Duplicate distance", settings.duplicateDistance,
            settings.minDuplicateDistance);
      }
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    invertCriteria = requiresInversion(settings.criteriaIndex);
    minCriteria = (invertCriteria) ? -settings.criteriaLimit : settings.criteriaLimit;
    invertScore = requiresInversion(settings.scoreIndex);

    return !gd.invalidNumber();
  }

  private void createResultsPrefix2() {
    createResultsPrefix2(settings.failCount, residualsThreshold, settings.duplicateDistance);
  }

  private void createResultsPrefix2(int failCount, double residualsThreshold,
      double duplicateDistance) {
    resultsPrefix2 = buildResultsPrefix2(failCount, residualsThreshold, duplicateDistance);
    if (!TextUtils.isNullOrEmpty(settings.resultsTitle)) {
      limitFailCount = settings.resultsTitle + ", ";
    } else {
      limitFailCount = "";
    }
    limitFailCount += "f=" + failCount;
    limitFailCount += ", r=" + MathUtils.rounded(residualsThreshold);
  }

  private String buildResultsPrefix2(int failCount, double residualsThreshold,
      double duplicateDistance) {
    return "\t" + failCount + "\t" + MathUtils.rounded(residualsThreshold) + "\t"
        + MathUtils.rounded(duplicateDistance)
        + ((settings.duplicateDistanceAbsolute) ? " (absolute)" : " (relative)");
  }

  private boolean showIterationDialog() {
    final GenericDialog gd = new GenericDialog(TITLE);
    final StringBuilder sb = new StringBuilder();
    sb.append("Iterate ").append(BenchmarkSpotFit.TITLE).append(" & ").append(TITLE).append(".\n");
    sb.append(BenchmarkSpotFit.TITLE)
        .append(" will be run once interactively if results cannot be loaded.\n");
    sb.append(TITLE).append(" will be run once interactively to obtain settings.\n \n");
    sb.append("Configure the convergence criteria for iteration:");
    gd.addMessage(sb.toString());
    gd.addNumericField("Score_Tolerance", settings.iterationScoreTolerance, -1);
    gd.addNumericField("Filter_Tolerance", settings.iterationFilterTolerance, -1);
    gd.addCheckbox("Compare_Results", settings.iterationCompareResults);
    gd.addNumericField("Compare_Distance", settings.iterationCompareDistance, 2);
    gd.addNumericField("Iter_Max_Iterations", settings.iterationMaxIterations, 0);
    gd.addMessage("Configure how the parameter range is updated per iteration:");
    gd.addSlider("Min_range_reduction", 0.05, 1, settings.iterationMinRangeReduction);
    gd.addSlider("Min_range_reduction_iteration", 1, 10,
        settings.iterationMinRangeReductionIteration);
    gd.addCheckbox("Converge_before_refit", settings.iterationConvergeBeforeRefit);

    gd.addHelp(HelpUrls.getUrl("iterate-filter-analysis"));
    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.iterationScoreTolerance = gd.getNextNumber();
    settings.iterationFilterTolerance = gd.getNextNumber();
    settings.iterationCompareResults = gd.getNextBoolean();
    settings.iterationCompareDistance = Math.abs(gd.getNextNumber());
    settings.iterationMaxIterations = (int) gd.getNextNumber();
    settings.iterationMinRangeReduction = Math.abs(gd.getNextNumber());
    settings.iterationMinRangeReductionIteration = (int) Math.abs(gd.getNextNumber());
    settings.iterationConvergeBeforeRefit = gd.getNextBoolean();
    settings.save();

    return true;
  }

  private boolean showReportDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

    addSimulationData(gd);

    gd.addCheckbox("Show_table", settings.showResultsTable);
    gd.addCheckbox("Show_summary", settings.showSummaryTable);
    gd.addCheckbox("Clear_tables", settings.clearTables);
    gd.addSlider("Summary_top_n", 0, 20, settings.summaryTopN);
    gd.addCheckbox("Save_best_filter", settings.saveBestFilter);
    gd.addCheckbox("Save_template", settings.saveTemplate);
    gd.addCheckbox("Calculate_sensitivity", settings.calculateSensitivity);
    gd.addSlider("Delta", 0.01, 1, settings.delta);
    if (!simulationParameters.fixedDepth) {
      gd.addCheckbox("Depth_recall_analysis", settings.depthRecallAnalysis);
    }
    gd.addCheckbox("Score_analysis", settings.scoreAnalysis);
    gd.addChoice("Component_analysis", Settings.COMPONENT_ANALYSIS_OPTIONS,
        settings.componentAnalysis);
    gd.addStringField("Title", settings.resultsTitle, 20);
    final String[] labels = {"Show_TP", "Show_FP", "Show_FN"};
    gd.addCheckboxGroup(1, 3, labels,
        new boolean[] {settings.showTP, settings.showFP, settings.showFN});

    gd.addHelp(HelpUrls.getUrl("benchmark-filter-analysis"));
    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.showResultsTable = gd.getNextBoolean();
    settings.showSummaryTable = gd.getNextBoolean();
    settings.clearTables = gd.getNextBoolean();
    settings.summaryTopN = (int) Math.abs(gd.getNextNumber());
    settings.saveBestFilter = gd.getNextBoolean();
    settings.saveTemplate = gd.getNextBoolean();
    settings.calculateSensitivity = gd.getNextBoolean();
    settings.delta = gd.getNextNumber();
    if (!simulationParameters.fixedDepth) {
      settings.depthRecallAnalysis = gd.getNextBoolean();
    }
    settings.scoreAnalysis = gd.getNextBoolean();
    settings.componentAnalysis = gd.getNextChoiceIndex();
    settings.resultsTitle = gd.getNextString();
    settings.showTP = gd.getNextBoolean();
    settings.showFP = gd.getNextBoolean();
    settings.showFN = gd.getNextBoolean();
    settings.save();

    if (gd.invalidNumber()) {
      return false;
    }

    resultsPrefix = spotFitResults.resultPrefix + "\t" + settings.resultsTitle + "\t";
    createResultsPrefix2();

    // Check there is one output
    if (!settings.showResultsTable && !settings.showSummaryTable && !settings.calculateSensitivity
        && !settings.saveBestFilter && !settings.saveTemplate) {
      IJ.error(TITLE, "No output selected");
      return false;
    }

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Delta", settings.delta);
      ParameterUtils.isBelow("Delta", settings.delta, 1);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return selectTableColumns();
  }

  private boolean showScoreDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

    addSimulationData(gd);

    // Get the last scored filter or default to the best filter
    createScoreFilter(false);

    gd.addSlider("Fail_count", 0, 20, settings.scoreFailCount);
    computeDoublets = BenchmarkSpotFit.getComputeDoublets();
    if (computeDoublets) {
      gd.addSlider("Residuals_threshold", 0.01, 1, settings.scoreResidualsThreshold);
    }
    gd.addNumericField("Duplicate_distance", settings.scoreDuplicateDistance, 2);

    gd.addTextAreas(uk.ac.sussex.gdsc.core.utils.XmlUtils.convertQuotes(scoreFilter.toXml()), null,
        6, 60);

    gd.addCheckbox("Reset_filter", false);
    gd.addCheckbox("Show_summary", settings.showSummaryTable);
    gd.addCheckbox("Clear_tables", settings.clearTables);
    gd.addCheckbox("Save_best_filter", settings.saveBestFilter);
    gd.addCheckbox("Save_template", settings.saveTemplate);
    gd.addCheckbox("Calculate_sensitivity", settings.calculateSensitivity);
    gd.addSlider("Delta", 0.01, 1, settings.delta);
    if (!simulationParameters.fixedDepth) {
      gd.addCheckbox("Depth_recall_analysis", settings.depthRecallAnalysis);
    }
    gd.addCheckbox("Score_analysis", settings.scoreAnalysis);
    gd.addChoice("Component_analysis", Settings.COMPONENT_ANALYSIS_OPTIONS,
        settings.componentAnalysis);
    gd.addStringField("Title", settings.resultsTitle, 20);
    final String[] labels = {"Show_TP", "Show_FP", "Show_FN"};
    gd.addCheckboxGroup(1, 3, labels,
        new boolean[] {settings.showTP, settings.showFP, settings.showFN});

    // Dialog to have a reset checkbox. This reverts back to the default.
    if (ImageJUtils.isShowGenericDialog()) {
      final Checkbox cb = (Checkbox) (gd.getCheckboxes().get(0));
      final Vector<TextField> v = gd.getNumericFields();
      final TextArea ta = gd.getTextArea1();
      cb.addItemListener(event -> {
        if (cb.getState()) {
          createScoreFilter(true);
          int index = 0;
          v.get(index++).setText(Integer.toString(settings.scoreFailCount));
          if (computeDoublets) {
            v.get(index++).setText(Double.toString(settings.scoreResidualsThreshold));
          }
          v.get(index++).setText(Double.toString(settings.scoreDuplicateDistance));
          ta.setText(uk.ac.sussex.gdsc.core.utils.XmlUtils.convertQuotes(scoreFilter.toXml()));
        }
      });
    }

    gd.addHelp(HelpUrls.getUrl("score-filter"));
    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.scoreFailCount = (int) Math.abs(gd.getNextNumber());
    if (computeDoublets) {
      settings.scoreResidualsThreshold = Math.abs(gd.getNextNumber());
    }
    settings.scoreDuplicateDistance = Math.abs(gd.getNextNumber());

    final String xml = gd.getNextText();
    boolean reset = gd.getNextBoolean();
    if (!reset) {
      // Try and parse the filter
      try {
        scoreFilter = (DirectFilter) Filter.fromXml(xml);
      } catch (final Exception ex) {
        // Invalid so reset
        reset = true;
      }
    }

    createScoreFilter(reset);

    settings.showSummaryTable = gd.getNextBoolean();
    settings.clearTables = gd.getNextBoolean();
    settings.saveBestFilter = gd.getNextBoolean();
    settings.saveTemplate = gd.getNextBoolean();
    settings.calculateSensitivity = gd.getNextBoolean();
    settings.delta = gd.getNextNumber();
    if (!simulationParameters.fixedDepth) {
      settings.depthRecallAnalysis = gd.getNextBoolean();
    }
    settings.scoreAnalysis = gd.getNextBoolean();
    settings.componentAnalysis = gd.getNextChoiceIndex();
    settings.resultsTitle = gd.getNextString();
    settings.showTP = gd.getNextBoolean();
    settings.showFP = gd.getNextBoolean();
    settings.showFN = gd.getNextBoolean();
    settings.save();

    if (gd.invalidNumber()) {
      return false;
    }

    resultsPrefix = spotFitResults.resultPrefix + "\t" + settings.resultsTitle + "\t";
    createResultsPrefix2(settings.scoreFailCount, settings.scoreResidualsThreshold,
        settings.scoreDuplicateDistance);

    // Check there is one output
    if (!settings.showSummaryTable && !settings.calculateSensitivity && !settings.saveBestFilter
        && !settings.saveTemplate) {
      IJ.error(TITLE, "No output selected");
      return false;
    }

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Delta", settings.delta);
      ParameterUtils.isBelow("Delta", settings.delta, 1);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return selectTableColumns();
  }

  /**
   * Run different filtering methods on a set of labelled peak results outputting performance
   * statistics on the success of the filter to an ImageJ table.
   *
   * <p>For each filter set a plot is shown of the score verses the filter value, thus filters
   * should be provided in ascending numerical order otherwise they are sorted.
   *
   * @param filterSets the filter sets
   * @param optimum the optimum
   * @param rangeReduction the range reduction
   * @return the best filter
   */
  private ComplexFilterScore analyse(List<FilterSet> filterSets, ComplexFilterScore optimum,
      double rangeReduction) {
    return analyse(filterSets, true, optimum, rangeReduction);
  }

  /**
   * Run different filtering methods on a set of labelled peak results outputting performance
   * statistics on the success of the filter to an ImageJ table.
   *
   * <p>For each filter set a plot is shown of the score verses the filter value, thus filters
   * should be provided in ascending numerical order otherwise they are sorted.
   *
   * @param filterSets the filter sets
   * @return the best filter
   */
  private ComplexFilterScore analyse(List<FilterSet> filterSets) {
    return analyse(filterSets, false, null, 0);
  }

  /**
   * Run different filtering methods on a set of labelled peak results outputting performance
   * statistics on the success of the filter to an ImageJ table.
   *
   * <p>For each filter set a plot is shown of the score verses the filter value, thus filters
   * should be provided in ascending numerical order otherwise they are sorted.
   *
   * @param filterSets the filter sets
   * @param iterative the iterative
   * @param optimum the optimum
   * @param rangeReduction the range reduction
   * @return the best filter
   */
  private ComplexFilterScore analyse(List<FilterSet> filterSets, boolean iterative,
      ComplexFilterScore optimum, double rangeReduction) {
    // Non-zero modes are used for the iterative optimisation which require new results
    boolean newResults = iterative;

    if (optimum != null) {
      // Non-interactive re-run when iterating
      filterAnalysisResult.scores.clear();
      runAnalysis(filterSets, optimum, rangeReduction);

      if (ImageJUtils.isInterrupted()) {
        return null;
      }
    } else {
      // Interactive run, this may be the first run during iterative optimisation
      if (settings.reset) {
        filterAnalysisResult.scores.clear();
      }

      createResultsWindow();

      // Only repeat analysis if necessary
      double evolveSetting = settings.evolve;
      if (settings.evolve == 1) {
        // The delta effects the step size for the Genetic Algorithm
        evolveSetting *= settings.delta;
      }
      final SettingsList settingsList = new SettingsList(filterSets, fitResultData.resultsList,
          settings.failCount, residualsThreshold, settings.duplicateDistance,
          settings.duplicateDistanceAbsolute, settings.plotTopN, settings.summaryDepth,
          settings.criteriaIndex, settings.criteriaLimit, settings.scoreIndex, evolveSetting);

      final boolean equalSettings = settingsList.equals(filterAnalysisResult.lastAnalyseSettings);

      if (!equalSettings || (settings.evolve != 0 && settings.repeatEvolve)) {
        newResults = true;
        filterAnalysisResult.lastAnalyseSettings = settingsList;

        runAnalysis(filterSets);

        if (ImageJUtils.isInterrupted()) {
          return null;
        }
      }
    }

    return reportResults(newResults);
  }

  /**
   * Run the optimum filter on a set of labelled peak results using various parameter settings
   * outputting performance statistics on the success of the filter to an ImageJ table.
   *
   * @param optimum the optimum
   * @return the best filter
   */
  private ComplexFilterScore analyseParameters(ComplexFilterScore optimum) {
    return analyseParameters(false, optimum, 0);
  }

  /**
   * Run the optimum filter on a set of labelled peak results using various parameter settings
   * outputting performance statistics on the success of the filter to an ImageJ table.
   *
   * @param optimum the optimum
   * @param rangeReduction the range reduction
   * @return the best filter
   */
  private ComplexFilterScore analyseParameters(ComplexFilterScore optimum, double rangeReduction) {
    return analyseParameters(true, optimum, rangeReduction);
  }

  /**
   * Run the optimum filter on a set of labelled peak results using various parameter settings
   * outputting performance statistics on the success of the filter to an ImageJ table.
   *
   * @param iterative the iterative
   * @param optimum the optimum
   * @param rangeReduction the range reduction
   * @return the best filter
   */
  private ComplexFilterScore analyseParameters(boolean iterative, ComplexFilterScore optimum,
      double rangeReduction) {
    // Non-zero modes are used for the iterative optimisation which require new results
    boolean newResults = iterative;

    if (!iterative) {
      // Interactive run, this may be the first run during iterative optimisation
      createResultsWindow();

      // Only repeat analysis if necessary
      double min = settings.minResidualsThreshold;
      double max = settings.maxResidualsThreshold;
      if (computeDoublets) {
        min = max = 0;
      }
      final SettingsList settingsList = new SettingsList(optimum, fitResultData.resultsList,
          settings.failCount, settings.minFailCount, settings.maxFailCount, residualsThreshold, min,
          max, settings.duplicateDistance, settings.duplicateDistanceAbsolute,
          settings.minDuplicateDistance, settings.maxDuplicateDistance, settings.summaryDepth,
          settings.criteriaIndex, settings.criteriaLimit, settings.scoreIndex,
          settings.searchParam);

      if (settings.repeatSearch
          || !settingsList.equals(filterAnalysisResult.lastAnalyseParametersSettings)) {
        newResults = true;
        filterAnalysisResult.lastAnalyseParametersSettings = settingsList;
      }
    }

    if (newResults) {
      optimum = runParameterAnalysis(iterative, optimum, rangeReduction);

      if (optimum == null || ImageJUtils.isInterrupted()) {
        return null;
      }
    }

    return reportResults(newResults, optimum);
  }

  private ComplexFilterScore reportResults(boolean newResults) {
    return reportResults(newResults, new ArrayList<>(filterAnalysisResult.bestFilter.values()));
  }

  private ComplexFilterScore reportResults(boolean newResults, ComplexFilterScore optimum) {
    return reportResults(newResults, Arrays.asList(optimum));
  }

  private ComplexFilterScore reportResults(boolean newResults, List<ComplexFilterScore> filters) {
    if (filters.isEmpty()) {
      IJ.log("Warning: No filters pass the criteria");
      return null;
    }

    getCoordinateStore();

    Collections.sort(filters);

    FractionClassificationResult topFilterClassificationResult = null;
    ArrayList<FractionalAssignment[]> topFilterResults = null;
    String topFilterSummary = null;
    if (settings.showSummaryTable || settings.saveTemplate) {
      final Consumer<String> summaryWindow = createSummaryWindow();
      int count = 0;
      final double range = (settings.summaryDepth / simulationParameters.pixelPitch) * 0.5;
      final int[] np = {0};
      fitResultData.depthStats.forEach(depth -> {
        if (Math.abs(depth) < range) {
          np[0]++;
        }
      });
      for (final ComplexFilterScore fs : filters) {
        final ArrayList<FractionalAssignment[]> list =
            new ArrayList<>(fitResultData.resultsList.length);
        final FractionClassificationResult r = scoreFilter(fs.getFilter(), defaultMinimalFilter,
            fitResultData.resultsList, list, coordinateStore);
        final StringBuilder sb = createResult(fs.getFilter(), r);

        if (topFilterResults == null) {
          topFilterResults = list;
          topFilterClassificationResult = r;
        }

        // Show the recall at the specified depth. Sum the distance and signal factor of all scored
        // spots.
        int scored = 0;
        double tp = 0;
        double distance = 0;
        double sf = 0;
        double rmsd = 0;

        final SimpleRegression regression = new SimpleRegression(false);
        for (final FractionalAssignment[] assignments : list) {
          if (assignments == null) {
            continue;
          }
          for (int i = 0; i < assignments.length; i++) {
            final CustomFractionalAssignment c = (CustomFractionalAssignment) assignments[i];
            if (Math.abs(c.peak.getZPosition()) <= range) {
              tp += c.getScore();
            }
            distance += c.distToTarget;
            sf += c.getSignalFactor();
            rmsd += c.distToTarget * c.distToTarget;
            regression.addData(c.peakResult.getSignal(), c.peak.getIntensity());
          }
          scored += assignments.length;
        }
        final double slope = regression.getSlope();

        sb.append('\t');
        sb.append(MathUtils.rounded(tp / np[0])).append('\t');
        sb.append(MathUtils.rounded(distance / scored)).append('\t');
        sb.append(MathUtils.rounded(sf / scored)).append('\t');
        // RMSD to be the root mean square deviation in a single dimension so divide by 2.
        // (This assumes 2D Euclidean distances.)
        sb.append(MathUtils.rounded(Math.sqrt(MathUtils.div0(rmsd / 2, scored)))).append('\t');
        sb.append(MathUtils.rounded(slope)).append('\t');
        if (fs.atLimit() != null) {
          sb.append(fs.atLimit());
        }

        String text = sb.toString();
        if (topFilterSummary == null) {
          topFilterSummary = text;
          if (!settings.showSummaryTable) {
            break;
          }
        }

        if (fs.time != 0) {
          sb.append('\t');
          sb.append(fs.algorithm);
          sb.append('\t');
          sb.append(org.apache.commons.lang3.time.DurationFormatUtils.formatDurationHMS(fs.time));
        } else {
          sb.append("\t\t");
        }

        if (fs.paramTime != 0) {
          sb.append('\t');
          sb.append(fs.getParamAlgorithm());
          sb.append('\t');
          sb.append(
              org.apache.commons.lang3.time.DurationFormatUtils.formatDurationHMS(fs.paramTime));
        } else {
          sb.append("\t\t");
        }
        text = sb.toString();

        summaryWindow.accept(text);
        count++;
        if (settings.summaryTopN > 0 && count >= settings.summaryTopN) {
          break;
        }
      }
      // Add a spacer to the summary table if we have multiple results
      if (count > 1 && settings.showSummaryTable) {
        summaryWindow.accept("");
      }
    }

    final DirectFilter localBestFilter = filters.get(0).getFilter();
    if (settings.saveBestFilter) {
      saveFilter(localBestFilter);
    }

    if (topFilterClassificationResult == null) {
      topFilterResults = new ArrayList<>(fitResultData.resultsList.length);
      scoreFilter(localBestFilter, defaultMinimalFilter, fitResultData.resultsList,
          topFilterResults, coordinateStore);
    }
    if (newResults || filterAnalysisResult.scores.isEmpty()) {
      filterAnalysisResult.scores.add(new FilterResult(settings.failCount, residualsThreshold,
          settings.duplicateDistance, settings.duplicateDistanceAbsolute, filters.get(0)));
    }

    if (settings.saveTemplate) {
      saveTemplate(topFilterSummary);
    }

    showPlots();
    calculateSensitivity();
    topFilterResults = depthAnalysis(topFilterResults, localBestFilter);
    topFilterResults = scoreAnalysis(topFilterResults, localBestFilter);
    componentAnalysis(filters.get(0));
    PreprocessedPeakResult[] filterResults = null;
    if (isShowOverlay()) {
      filterResults = showOverlay(topFilterResults, localBestFilter);
    }
    saveResults(filterResults, localBestFilter);

    wo.tile();

    return filters.get(0);
  }

  private void runAnalysis(List<FilterSet> filterSets) {
    runAnalysis(filterSets, null, 0);
  }

  private void runAnalysis(List<FilterSet> filterSets, ComplexFilterScore optimum,
      double rangeReduction) {
    filterAnalysisResult.plots.clear();
    filterAnalysisResult.bestFilter.clear();

    getCoordinateStore();

    filterAnalysisStopWatch = StopWatch.createStarted();
    IJ.showStatus("Analysing filters ...");
    int setNumber = 0;
    final DirectFilter currentOptimum = (optimum != null) ? optimum.result.filter : null;
    for (final FilterSet filterSet : filterSets) {
      setNumber++;
      if (filterAnalysis(filterSet, setNumber, currentOptimum, rangeReduction) < 0) {
        break;
      }
    }
    filterAnalysisStopWatch.stop();
    ImageJUtils.finished();

    final String timeString = filterAnalysisStopWatch.toString();
    IJ.log("Filter analysis time : " + timeString);
  }

  /**
   * Run the optimum filter on a set of labelled peak results using various parameter settings
   * outputting performance statistics on the success of the filter to an ImageJ table.
   *
   * <p>If a new optimum is found the class level static parameters are updated.
   *
   * @param nonInteractive True if non interactive
   * @param optimum the optimum
   * @param rangeReduction the range reduction
   * @return the best filter
   */
  private ComplexFilterScore runParameterAnalysis(boolean nonInteractive,
      ComplexFilterScore optimum, double rangeReduction) {
    parameterAnalysisStopWatch = StopWatch.createStarted();
    IJ.showStatus("Analysing parameters ...");
    optimum = parameterAnalysis(nonInteractive, optimum, rangeReduction);
    parameterAnalysisStopWatch.stop();
    ImageJUtils.finished();

    final String timeString = parameterAnalysisStopWatch.toString();
    IJ.log("Parameter analysis time : " + timeString);
    return optimum;
  }

  private CoordinateStore getCoordinateStore() {
    if (coordinateStore == null) {
      coordinateStore = createCoordinateStore();
    }
    return coordinateStore;
  }

  private static int countFilters(List<FilterSet> filterSets) {
    int count = 0;
    for (final FilterSet filterSet : filterSets) {
      count += filterSet.size();
    }
    return count;
  }

  private void showPlots() {
    if (filterAnalysisResult.plots.isEmpty()) {
      return;
    }

    // Display the top N plots
    final int[] list = new int[filterAnalysisResult.plots.size()];
    int index = 0;
    for (final NamedPlot p : filterAnalysisResult.plots) {
      final Plot plot = new Plot(p.name, p.xAxisName, Settings.COLUMNS[settings.scoreIndex]);
      plot.setLimits(p.xValues[0], p.xValues[p.xValues.length - 1], 0, 1);
      plot.setColor(Color.RED);
      plot.addPoints(p.xValues, p.yValues, Plot.LINE);
      plot.setColor(Color.BLUE);
      plot.addPoints(p.xValues, p.yValues, Plot.CROSS);
      final PlotWindow plotWindow = ImageJUtils.display(p.name, plot);
      list[index++] = plotWindow.getImagePlus().getID();
    }
    WindowOrganiser.tileWindows(list);
  }

  private void calculateSensitivity() {
    if (!settings.calculateSensitivity) {
      return;
    }
    if (!filterAnalysisResult.bestFilter.isEmpty()) {
      final Consumer<String> output = createSensitivityWindow();

      final Ticker ticker = ImageJUtils.createTicker(filterAnalysisResult.bestFilter.size(), 0,
          "Calculating sensitivity ...");
      for (final String type : filterAnalysisResult.bestFilter.keySet()) {

        final DirectFilter filter = filterAnalysisResult.bestFilter.get(type).getFilter();

        FractionClassificationResult score =
            scoreFilter(filter, defaultMinimalFilter, fitResultData.resultsList, coordinateStore);
        score = getOriginalScore(score);

        final String message = type + "\t\t\t" + MathUtils.rounded(score.getJaccard(), 4) + "\t\t"
            + MathUtils.rounded(score.getPrecision(), 4) + "\t\t"
            + MathUtils.rounded(score.getRecall(), 4);

        output.accept(message);

        // List all the parameters that can be adjusted.
        final int parameters = filter.getNumberOfParameters();
        for (int index = 0; index < parameters; index++) {
          // For each parameter compute as upward + downward delta and get the average gradient
          final DirectFilter higher = (DirectFilter) filter.adjustParameter(index, settings.delta);
          final DirectFilter lower = (DirectFilter) filter.adjustParameter(index, -settings.delta);

          FractionClassificationResult scoreHigher =
              scoreFilter(higher, defaultMinimalFilter, fitResultData.resultsList, coordinateStore);
          scoreHigher = getOriginalScore(scoreHigher);
          FractionClassificationResult scoreLower =
              scoreFilter(lower, defaultMinimalFilter, fitResultData.resultsList, coordinateStore);
          scoreLower = getOriginalScore(scoreLower);

          final StringBuilder sb = new StringBuilder();
          sb.append('\t').append(filter.getParameterName(index)).append('\t');
          sb.append(MathUtils.rounded(filter.getParameterValue(index), 4)).append('\t');

          final double dx1 = higher.getParameterValue(index) - filter.getParameterValue(index);
          final double dx2 = filter.getParameterValue(index) - lower.getParameterValue(index);
          addSensitivityScore(sb, score.getJaccard(), scoreHigher.getJaccard(),
              scoreLower.getJaccard(), dx1, dx2);
          addSensitivityScore(sb, score.getPrecision(), scoreHigher.getPrecision(),
              scoreLower.getPrecision(), dx1, dx2);
          addSensitivityScore(sb, score.getRecall(), scoreHigher.getRecall(),
              scoreLower.getRecall(), dx1, dx2);

          output.accept(sb.toString());
        }

        ticker.tick();
      }

      final String message = "-=-=-=-";
      output.accept(message);

      ImageJUtils.finished();
    }
  }

  private static void addSensitivityScore(StringBuilder sb, double score, double score1,
      double score2, double dx1, double dx2) {
    // Use absolute in case this is not a local maximum. We are mainly interested in how
    // flat the curve is at this point in relation to parameter changes.
    double abs = 0;
    double dydx = 0;
    int count = 0;
    if (dx1 > 0) {
      final double abs1 = Math.abs(score - score1);
      final double dydx1 = abs1 / dx1;
      abs += abs1;
      dydx += dydx1;
      count++;
    }
    if (dx2 > 0) {
      final double abs2 = Math.abs(score - score2);
      final double dydx2 = abs2 / dx2;
      abs += abs2;
      dydx += dydx2;
      count++;
    }

    double relativeSensitivity = 0;
    double sensitivity = 0;
    if (count != 0) {
      relativeSensitivity = abs / count;
      sensitivity = dydx / count;
    }

    sb.append(MathUtils.rounded(relativeSensitivity, 4)).append('\t');
    sb.append(MathUtils.rounded(sensitivity, 4)).append('\t');
  }

  private void createResultsWindow() {
    if (!settings.showResultsTable) {
      return;
    }

    if (isHeadless) {
      IJ.log(createResultsHeader(false));
    } else {
      resultsWindow = ImageJUtils.refresh(resultsWindowRef, () -> {
        final String header = createResultsHeader(false);
        return new TextWindow(TITLE + " Results", header, "", 900, 300);
      });
      if (settings.clearTables) {
        resultsWindow.getTextPanel().clear();
      }
    }
  }

  private Consumer<String> createSummaryWindow() {
    if (!settings.showSummaryTable) {
      return null;
    }

    if (isHeadless) {
      IJ.log(createResultsHeader(true));
      return IJ::log;
    }
    final TextWindow window = ImageJUtils.refresh(summaryWindowRef, () -> {
      final String header = createResultsHeader(true);
      return new TextWindow(TITLE + " Summary", header, "", 900, 300);
    });
    if (settings.clearTables) {
      window.getTextPanel().clear();
    }
    return window::append;
  }

  private void createGaWindow() {
    if (isHeadless) {
      String header = createResultsHeader(false);
      header += "\tIteration";
      IJ.log(header);
      gaWindow = IJ::log;
    } else {
      final TextWindow window = ImageJUtils.refresh(gaWindowRef, () -> {
        String header = createResultsHeader(false);
        header += "\tIteration";
        return new TextWindow(TITLE + " Evolution", header, "", 900, 300);
      });
      if (settings.clearTables) {
        window.getTextPanel().clear();
      }
      gaWindow = window::append;
    }
  }

  private Consumer<String> createComponentAnalysisWindow() {
    if (isHeadless) {
      IJ.log(createComponentAnalysisHeader());
      return IJ::log;
    }
    final TextWindow window = ImageJUtils.refresh(componentAnalysisWindowRef, () -> {
      final String header = createComponentAnalysisHeader();
      return new TextWindow(TITLE + " Component Analysis", header, "", 900, 300);
    });
    if (settings.clearTables) {
      window.getTextPanel().clear();
    }
    return window::append;
  }

  private String createComponentAnalysisHeader() {
    String header = createResultsHeader(false);
    header += "\tSize\tName\tValue\tLimit\t% Criteria\t% Score\tTime\tOverlap P\t"
        + "Overlap R\tOverlap J\tNames";
    return header;
  }

  private static void addToComponentAnalysisWindow(Consumer<String> output,
      ComplexFilterScore filterScore, ComplexFilterScore bestFilterScore, String[] names) {
    final FilterScoreResult result = filterScore.result;
    final StringBuilder sb = new StringBuilder(result.text);
    sb.append('\t').append(filterScore.size);
    final int index = filterScore.index;
    if (index != -1) {
      sb.append('\t').append(names[index]);
      sb.append('\t').append(MathUtils.rounded(filterScore.getFilter().getParameterValue(index)));
      sb.append('\t');
      if (bestFilterScore.atLimit != null) {
        sb.append(bestFilterScore.atLimit(index));
      }
    } else {
      // Broken
      sb.append("\t\t");
    }
    sb.append('\t').append(MathUtils.rounded(100.0 * result.criteria / bestFilterScore.criteria));
    sb.append('\t').append(MathUtils.rounded(100.0 * result.score / bestFilterScore.score));
    sb.append('\t').append(filterScore.time);
    sb.append('\t').append(MathUtils.rounded(filterScore.result2.getPrecision()));
    sb.append('\t').append(MathUtils.rounded(filterScore.result2.getRecall()));
    sb.append('\t').append(MathUtils.rounded(filterScore.result2.getJaccard()));
    sb.append('\t');
    for (int i = 0; i < filterScore.combinations.length; i++) {
      if (i != 0) {
        sb.append(", ");
      }
      sb.append(names[filterScore.combinations[i]]);
    }
    output.accept(sb.toString());
  }

  private String createResultsHeader(boolean summary) {
    final StringBuilder sb = new StringBuilder(BenchmarkSpotFit.getTablePrefix());
    sb.append(
        "\tTitle\tName\tFail\tRes\tDup D\tLower D (nm)\tUpper D (nm)\tLower factor\tUpper factor");

    for (int i = 0; i < Settings.COLUMNS.length; i++) {
      if (settings.showColumns[i]) {
        sb.append('\t').append(Settings.COLUMNS[i]);
      }
    }

    if (summary) {
      sb.append("\tDepth Recall\tDistance\tSignal Factor\tRMSD\tSlope\tAt limit\tEvolve\t"
          + "Time\tSearch\tTime");
    }
    return sb.toString();
  }

  private Consumer<String> createSensitivityWindow() {
    if (isHeadless) {
      IJ.log(createSensitivityHeader());
      return IJ::log;
    }
    return ImageJUtils.refresh(sensitivityWindowRef, () -> new TextWindow(TITLE + " Sensitivity",
        createSensitivityHeader(), "", 900, 300))::append;
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

  private int filterAnalysis(FilterSet filterSet, int setNumber, DirectFilter currentOptimum,
      double rangeReduction) {
    // Check if the filters are the same so allowing optimisation
    final boolean allSameType = filterSet.allSameType();

    this.gaResultsList = fitResultData.resultsList;

    Chromosome<FilterScore> best = null;
    String algorithm = "";

    // All the search algorithms use search dimensions.
    // Create search dimensions if needed (these are used for testing if the optimum is at the
    // limit).
    searchScoreFilter = null;
    strengthLower = null;
    strengthUpper = null;
    FixedDimension[] originalDimensions = null;
    boolean rangeInput = false;
    boolean[] disabled = null;
    double[][] seed = null;
    // This flag is set if the analysis is not interactive. This occurs when running after the
    // first iteration of an iterative analysis.
    boolean nonInteractive = false;
    if (allSameType) {
      // There should always be 1 filter
      searchScoreFilter = (DirectFilter) filterSet.getFilters().get(0);
      final int n = searchScoreFilter.getNumberOfParameters();

      // Option to configure a range
      rangeInput = filterSet.getName().contains("Range");
      final double[] range = new double[n];

      if (rangeInput && filterSet.size() == 4) {
        originalDimensions = new FixedDimension[n];
        // This is used as min/lower/upper/max
        final Filter minF = searchScoreFilter;
        final Filter lowerF = filterSet.getFilters().get(1);
        final Filter upperF = filterSet.getFilters().get(2);
        final Filter maxF = filterSet.getFilters().get(3);
        for (int i = 0; i < n; i++) {
          final double min = minF.getParameterValue(i);
          final double lower = lowerF.getParameterValue(i);
          final double upper = upperF.getParameterValue(i);
          range[i] = upper - lower;
          final double max = maxF.getParameterValue(i);
          final double minIncrement = searchScoreFilter.getParameterIncrement(i);
          try {
            originalDimensions[i] = new FixedDimension(min, max, minIncrement, lower, upper);
          } catch (final IllegalArgumentException ex) {
            ImageJUtils.log(TITLE + " : Unable to configure dimension [%d] %s: " + ex.getMessage(),
                i, searchScoreFilter.getParameterName(i));
            originalDimensions = null;
            rangeInput = false;
            break;
          }
        }
      }

      if (rangeInput && (filterSet.size() == 3 || filterSet.size() == 2)) {
        originalDimensions = new FixedDimension[n];
        // This is used as lower/upper[/increment]
        final Filter lowerF = searchScoreFilter;
        final Filter upperF = filterSet.getFilters().get(1);

        for (int i = 0; i < n; i++) {
          // Do not disable if the increment is not set. This is left to the user to decide
          // which parameters to optimise with the enabled checkboxes in the dialog.

          final double lower = lowerF.getParameterValue(i);
          final double upper = upperF.getParameterValue(i);
          range[i] = upper - lower;
          final ParameterType type = searchScoreFilter.getParameterType(i);
          final double min = BenchmarkSpotFit.getMin(type);
          final double max = BenchmarkSpotFit.getMax(type);
          final double minIncrement = searchScoreFilter.getParameterIncrement(i);
          try {
            originalDimensions[i] = new FixedDimension(min, max, minIncrement, lower, upper);
          } catch (final IllegalArgumentException ex) {
            ImageJUtils.log(TITLE + " : Unable to configure dimension [%d] %s: " + ex.getMessage(),
                i, searchScoreFilter.getParameterName(i));
            originalDimensions = null;
            rangeInput = false;
            break;
          }
        }
      }

      // Get the dimensions from the filters
      if (originalDimensions == null) {
        originalDimensions = new FixedDimension[n];

        // Allow inputing a filter set (e.g. saved from previous optimisation)
        // Find the limits in the current scores
        final double[] lower = searchScoreFilter.getParameters().clone();
        final double[] upper = lower.clone();
        // Allow the SearchSpace algorithms to be seeded with an initial population
        // for the first evaluation of the optimum. This is done when the input filter
        // set is not a range.
        seed = new double[filterSet.size()][];
        int count = 0;
        for (final Filter f : filterSet.getFilters()) {
          final double[] point = f.getParameters();
          seed[count++] = point;
          for (int j = 0; j < lower.length; j++) {
            if (lower[j] > point[j]) {
              lower[j] = point[j];
            }
            if (upper[j] < point[j]) {
              upper[j] = point[j];
            }
          }
        }

        // Get the search dimensions from the data.
        // Min/max must be set using values from BenchmarkSpotFit.
        boolean hasRange = false;
        for (int i = 0; i < n; i++) {
          if (lower[i] == upper[i]) {
            // Not enabled
            originalDimensions[i] = new FixedDimension(lower[i]);
            continue;
          }
          hasRange = true;
          final ParameterType type = searchScoreFilter.getParameterType(i);
          double min = BenchmarkSpotFit.getMin(type);
          double max = BenchmarkSpotFit.getMax(type);
          final double minIncrement = searchScoreFilter.getParameterIncrement(i);
          if (min > lower[i]) {
            min = lower[i];
          }
          if (max < upper[i]) {
            max = upper[i];
          }
          try {
            originalDimensions[i] = new FixedDimension(min, max, minIncrement, lower[i], upper[i]);
          } catch (final IllegalArgumentException ex) {
            ImageJUtils.log(TITLE + " : Unable to configure dimension [%d] %s: " + ex.getMessage(),
                i, searchScoreFilter.getParameterName(i));
            originalDimensions = null;
            break;
          }
        }

        if (!hasRange || originalDimensions == null) {
          // Failed to work out the dimensions. No optimisation will be possible.
          originalDimensions = null;

          // Sort so that the filters are in a nice order for reporting
          filterSet.sort();

          // This will not be used when the dimensions are null
          seed = null;
        }
      }

      if (originalDimensions != null) {
        // Use the current optimum if we are doing a range optimisation
        if (currentOptimum != null && rangeInput
            && currentOptimum.getType().equals(searchScoreFilter.getType())
            && settings.evolve != 0) {
          // Suppress dialogs and use the current settings
          nonInteractive = true;

          final double[] p = currentOptimum.getParameters();

          // If the optimum is that same filter type as the filter set, and we are using
          // a search method then centre the dimensions on the optimum.
          // Note:
          // Enrichment search and GA just use the upper and lower bounds to seed the population
          // These can use FixedDimension and we add the current optimum to the seed.
          // Range search uses SearchDimension and we must centre on the optimum after creation.
          for (int i = 0; i < originalDimensions.length; i++) {
            final double centre = p[i];
            double rangeFactor = 0;
            if (originalDimensions[i].isActive()) {
              // Set the range around the centre.
              // This uses the range for each param when we read the filters.
              rangeFactor = range[i];
              // Optionally reduce the width of the dimensions.
              if (rangeReduction > 0 && rangeReduction < 1) {
                rangeFactor *= rangeReduction;
              }
            }
            final double lower = centre - rangeFactor * 0.5;
            final double upper = centre + rangeFactor * 0.5;
            originalDimensions[i] = originalDimensions[i].create(lower, upper);
          }
        }

        // Store the dimensions so we can do an 'at limit' check
        disabled = new boolean[originalDimensions.length];
        strengthLower = new double[originalDimensions.length];
        strengthUpper = new double[originalDimensions.length];
        for (int i = 0; i < disabled.length; i++) {
          disabled[i] = !originalDimensions[i].isActive();
          strengthLower[i] = originalDimensions[i].lower;
          strengthUpper[i] = originalDimensions[i].upper;
        }
      }
    } else {
      // Sort so that the filters are in a nice order for reporting
      filterSet.sort();
    }

    analysisStopWatch = StopWatch.createStarted();

    // Note:
    // The input filters have yet to be scored.

    if (settings.evolve == 1 && originalDimensions != null) {
      // Collect parameters for the genetic algorithm
      pauseFilterTimer();

      // Remember the step size settings
      double[] stepSize = stepSizeMap.get(setNumber);
      if (stepSize == null || stepSize.length != searchScoreFilter.length()) {
        stepSize = searchScoreFilter.mutationStepRange().clone();
        for (int j = 0; j < stepSize.length; j++) {
          stepSize[j] *= settings.delta;
        }
        // See if the same number of parameters have been optimised in other algorithms
        final boolean[] enabled = searchRangeMap.get(setNumber);
        if (enabled != null && enabled.length == stepSize.length) {
          for (int j = 0; j < stepSize.length; j++) {
            if (!enabled[j]) {
              stepSize[j] *= -1;
            }
          }
        }
      }

      GenericDialog gd = null;
      final int[] indices = searchScoreFilter.getChromosomeParameters();
      boolean runAlgorithm = nonInteractive;
      if (!nonInteractive) {
        // Ask the user for the mutation step parameters.
        gd = new GenericDialog(TITLE);
        final String prefix = setNumber + "_";
        gd.addMessage(
            "Configure the genetic algorithm for [" + setNumber + "] " + filterSet.getName());
        gd.addNumericField(prefix + "Population_size", settings.populationSize, 0);
        gd.addNumericField(prefix + "Failure_limit", settings.failureLimit, 0);
        gd.addNumericField(prefix + "Tolerance", settings.tolerance, -1);
        gd.addNumericField(prefix + "Converged_count", settings.convergedCount, 0);
        gd.addSlider(prefix + "Mutation_rate", 0.05, 1, settings.mutationRate);
        gd.addSlider(prefix + "Crossover_rate", 0.05, 1, settings.crossoverRate);
        gd.addSlider(prefix + "Mean_children", 0.05, 3, settings.meanChildren);
        gd.addSlider(prefix + "Selection_fraction", 0.05, 0.5, settings.selectionFraction);
        gd.addCheckbox(prefix + "Ramped_selection", settings.rampedSelection);
        gd.addCheckbox(prefix + "Save_option", settings.saveOption);

        gd.addMessage("Configure the step size for each parameter");
        for (int j = 0; j < indices.length; j++) {
          // Do not mutate parameters that were not expanded, i.e. the input did not vary them.
          final double step = (originalDimensions[indices[j]].isActive()) ? stepSize[j] : 0;
          gd.addNumericField(getDialogName(prefix, searchScoreFilter, indices[j]), step, 2);
        }

        gd.showDialog();
        runAlgorithm = !gd.wasCanceled();
      }

      if (runAlgorithm) {
        // Used to create random sample
        final FixedDimension[] dimensions =
            Arrays.copyOf(originalDimensions, originalDimensions.length);
        if (!nonInteractive) {
          if (gd == null) {
            throw new IllegalStateException("The dialog has not been shown");
          }
          settings.populationSize = (int) Math.abs(gd.getNextNumber());
          if (settings.populationSize < 10) {
            settings.populationSize = 10;
          }
          settings.failureLimit = (int) Math.abs(gd.getNextNumber());
          settings.tolerance = gd.getNextNumber();
          settings.convergedCount = (int) gd.getNextNumber(); // Allow negatives
          settings.mutationRate = Math.abs(gd.getNextNumber());
          settings.crossoverRate = Math.abs(gd.getNextNumber());
          settings.meanChildren = Math.abs(gd.getNextNumber());
          settings.selectionFraction = Math.abs(gd.getNextNumber());
          settings.rampedSelection = gd.getNextBoolean();
          settings.saveOption = gd.getNextBoolean();

          for (int j = 0; j < indices.length; j++) {
            stepSize[j] = gd.getNextNumber();
          }
          // Store for repeat analysis
          stepSizeMap.put(setNumber, stepSize);
        }

        if (disabled == null) {
          disabled = new boolean[originalDimensions.length];
        }
        for (int j = 0; j < indices.length; j++) {
          // Disable values with a negative step size.
          // A zero step size will keep the parameter but prevent range mutation.
          if (stepSize[j] < 0) {
            dimensions[indices[j]] =
                new FixedDimension(searchScoreFilter.getDisabledParameterValue(indices[j]));
            disabled[indices[j]] = true;
          }
        }

        // Create the genetic algorithm
        final UniformRandomProvider rng = UniformRandomProviders.create();
        final SimpleMutator<FilterScore> mutator = new SimpleMutator<>(rng, settings.mutationRate);
        // Override the settings with the step length, a min of zero and the configured upper
        final double[] upper = searchScoreFilter.upperLimit();
        mutator.overrideChromosomeSettings(stepSize, new double[stepSize.length], upper);
        final Recombiner<FilterScore> recombiner =
            new SimpleRecombiner<>(rng, settings.crossoverRate, settings.meanChildren);
        SelectionStrategy<FilterScore> selectionStrategy;
        // If the initial population is huge ensure that the first selection culls to the correct
        // size
        final int selectionMax = (int) (settings.selectionFraction * settings.populationSize);
        if (settings.rampedSelection) {
          selectionStrategy =
              new RampedSelectionStrategy<>(rng, settings.selectionFraction, selectionMax);
        } else {
          selectionStrategy =
              new SimpleSelectionStrategy<>(rng, settings.selectionFraction, selectionMax);
        }
        final ToleranceChecker<FilterScore> gaChecker = new InterruptChecker(settings.tolerance,
            settings.tolerance * 1e-3, settings.convergedCount);

        // Create new random filters if the population is initially below the population size
        List<Filter> filters = filterSet.getFilters();
        if (filterSet.size() < settings.populationSize) {
          filters = new ArrayList<>(settings.populationSize);
          // Add the existing filters if they are not a range input file
          if (!rangeInput) {
            filters.addAll(filterSet.getFilters());
          }
          // Add current optimum to seed
          if (currentOptimum != null && nonInteractive) {
            filters.add(currentOptimum);
          }
          // The GA does not use the min interval grid so sample without rounding
          final double[][] sample = SearchSpace.sampleWithoutRounding(dimensions,
              settings.populationSize - filters.size(), null);
          filters.addAll(searchSpaceToFilters(sample));
        }
        gaPopulation = new Population<>(filters);
        gaPopulation.setPopulationSize(settings.populationSize);
        gaPopulation.setFailureLimit(settings.failureLimit);
        selectionStrategy.setTracker(this);

        // Evolve
        algorithm = Settings.EVOLVE_OPTIONS[settings.evolve];
        gaStatusPrefix = algorithm + " [" + setNumber + "] " + filterSet.getName() + " ... ";
        gaIteration = 0;
        gaPopulation.setTracker(this);

        createGaWindow();
        resumeFilterTimer();

        best = gaPopulation.evolve(mutator, recombiner, this, selectionStrategy, gaChecker);

        // In case optimisation was stopped
        IJ.resetEscape();

        if (best != null) {
          // The GA may produce coordinates off the min interval grid
          best = enumerateMinInterval(best, stepSize, indices);

          // Now update the filter set for final assessment
          filterSet = new FilterSet(filterSet.getName(),
              populationToFilters(gaPopulation.getIndividuals()));

          // Option to save the filters
          if (settings.saveOption) {
            saveFilterSet(filterSet, setNumber, !nonInteractive);
          }
        }
      } else {
        resumeFilterTimer();
      }
    }

    if ((settings.evolve == 2 || settings.evolve == 4) && originalDimensions != null) {
      // Collect parameters for the range search algorithm
      pauseFilterTimer();

      final boolean isStepSearch = settings.evolve == 4;

      // The step search should use a multi-dimension refinement and no range reduction
      SearchSpace.RefinementMode myRefinementMode = SearchSpace.RefinementMode.MULTI_DIMENSION;

      // Remember the enabled settings
      boolean[] enabled = searchRangeMap.get(setNumber);
      final int n = searchScoreFilter.getNumberOfParameters();
      if (enabled == null || enabled.length != n) {
        enabled = new boolean[n];
        Arrays.fill(enabled, true);
        // See if the same number of parameters have been optimised in other algorithms
        final double[] stepSize = stepSizeMap.get(setNumber);
        if (stepSize != null && enabled.length == stepSize.length) {
          for (int j = 0; j < stepSize.length; j++) {
            if (stepSize[j] < 0) {
              enabled[j] = false;
            }
          }
        }
      }

      GenericDialog gd = null;
      boolean runAlgorithm = nonInteractive;
      if (!nonInteractive) {
        // Ask the user for the search parameters.
        gd = new GenericDialog(TITLE);
        final String prefix = setNumber + "_";
        gd.addMessage("Configure the " + Settings.EVOLVE_OPTIONS[settings.evolve]
            + " algorithm for [" + setNumber + "] " + filterSet.getName());
        gd.addSlider(prefix + "Width", 1, 5, settings.rangeSearchWidth);
        if (!isStepSearch) {
          gd.addCheckbox(prefix + "Save_option", settings.saveOption);
          gd.addNumericField(prefix + "Max_iterations", settings.maxIterations, 0);
          final String[] modes =
              SettingsManager.getNames((Object[]) SearchSpace.RefinementMode.values());
          gd.addSlider(prefix + "Reduce", 0.01, 0.99, settings.rangeSearchReduce);
          gd.addChoice("Refinement", modes, modes[settings.refinementMode]);
        }
        gd.addNumericField(prefix + "Seed_size", settings.seedSize, 0);

        // Add choice of fields to optimise
        for (int i = 0; i < n; i++) {
          gd.addCheckbox(getDialogName(prefix, searchScoreFilter, i), enabled[i]);
        }

        gd.showDialog();
        runAlgorithm = !gd.wasCanceled();
      }

      if (runAlgorithm) {
        final SearchDimension[] dimensions = new SearchDimension[n];
        if (!nonInteractive) {
          if (gd == null) {
            throw new IllegalStateException("The dialog has not been shown");
          }
          settings.rangeSearchWidth = (int) gd.getNextNumber();
          if (!isStepSearch) {
            settings.saveOption = gd.getNextBoolean();
            settings.maxIterations = (int) gd.getNextNumber();
            settings.rangeSearchReduce = gd.getNextNumber();
            settings.refinementMode = gd.getNextChoiceIndex();
          }
          settings.seedSize = (int) gd.getNextNumber();
          for (int i = 0; i < n; i++) {
            enabled[i] = gd.getNextBoolean();
          }

          searchRangeMap.put(setNumber, enabled);
        }

        if (!isStepSearch) {
          myRefinementMode = SearchSpace.RefinementMode.values()[settings.refinementMode];
        }
        for (int i = 0; i < n; i++) {
          if (enabled[i]) {
            try {
              dimensions[i] = originalDimensions[i].create(settings.rangeSearchWidth);
              dimensions[i].setPad(true);
              // Prevent range reduction so that the step search just does a single refinement step
              dimensions[i].setReduceFactor((isStepSearch) ? 1 : settings.rangeSearchReduce);
              // Centre on current optimum
              if (nonInteractive && currentOptimum != null) {
                dimensions[i].setCentre(currentOptimum.getParameterValue(i));
              }
            } catch (final IllegalArgumentException ex) {
              IJ.error(TITLE,
                  String.format("Unable to configure dimension [%d] %s: " + ex.getMessage(), i,
                      searchScoreFilter.getParameterName(i)));
              return -1;
            }
          } else {
            dimensions[i] = new SearchDimension(searchScoreFilter.getDisabledParameterValue(i));
          }
        }
        if (disabled == null) {
          disabled = new boolean[originalDimensions.length];
        }
        for (int i = 0; i < disabled.length; i++) {
          disabled[i] = !dimensions[i].active;
        }

        // Check the number of combinations is OK
        long combinations = SearchSpace.countCombinations(dimensions);
        if (!nonInteractive && combinations > 10000) {
          gd = new GenericDialog(TITLE);
          ImageJUtils.addMessage(gd,
              "%d combinations for the configured dimensions.\n \nClick 'Yes' to optimise.",
              combinations);
          gd.enableYesNoCancel();
          gd.hideCancelButton();
          gd.showDialog();
          if (!gd.wasOKed()) {
            combinations = 0;
          }
        }

        if (combinations == 0) {
          resumeFilterTimer();
        } else {
          algorithm = Settings.EVOLVE_OPTIONS[settings.evolve] + " " + settings.rangeSearchWidth;
          gaStatusPrefix = algorithm + " [" + setNumber + "] " + filterSet.getName() + " ... ";
          gaIteration = 0;
          filterScoreOptimum = null;

          final SearchSpace ss = new SearchSpace();
          ss.setTracker(this);
          if (settings.seedSize > 0) {
            double[][] sample;
            // Add current optimum to seed
            if (nonInteractive && currentOptimum != null) {
              sample = new double[1][];
              sample[0] = currentOptimum.getParameters();
              seed = merge(seed, sample);
            }
            final int size = (seed == null) ? 0 : seed.length;
            // Sample without rounding as the seed will be rounded
            sample = SearchSpace.sampleWithoutRounding(dimensions, settings.seedSize - size, null);
            seed = merge(seed, sample);
          }
          // Note: If we have an optimum and we are not seeding this should not matter as the
          // dimensions have been centred on the current optimum
          ss.seed(seed);
          final ConvergenceChecker<FilterScore> checker =
              new InterruptConvergenceChecker(0, 0, settings.maxIterations);

          createGaWindow();
          resumeFilterTimer();

          final SearchResult<FilterScore> optimum =
              ss.search(dimensions, this, checker, myRefinementMode);

          // In case optimisation was stopped
          IJ.resetEscape();

          if (optimum != null) {
            best = ((SimpleFilterScore) optimum.getScore()).result.filter;

            // Now update the filter set for final assessment
            filterSet = new FilterSet(filterSet.getName(),
                searchSpaceToFilters((DirectFilter) best, ss.getSearchSpace()));

            // Option to save the filters
            if (settings.saveOption) {
              saveFilterSet(filterSet, setNumber, !nonInteractive);
            }
          }
        }
      } else {
        resumeFilterTimer();
      }
    }

    if (settings.evolve == 3 && originalDimensions != null) {
      // Collect parameters for the enrichment search algorithm
      pauseFilterTimer();

      boolean[] enabled = searchRangeMap.get(setNumber);
      final int n = searchScoreFilter.getNumberOfParameters();
      if (enabled == null || enabled.length != n) {
        enabled = new boolean[n];
        Arrays.fill(enabled, true);
        // See if the same number of parameters have been optimised in other algorithms
        final double[] stepSize = stepSizeMap.get(setNumber);
        if (stepSize != null && enabled.length == stepSize.length) {
          for (int j = 0; j < stepSize.length; j++) {
            if (stepSize[j] < 0) {
              enabled[j] = false;
            }
          }
        }
      }

      GenericDialog gd = null;
      boolean runAlgorithm = nonInteractive;
      if (!nonInteractive) {
        // Ask the user for the search parameters.
        gd = new GenericDialog(TITLE);
        final String prefix = setNumber + "_";
        gd.addMessage("Configure the enrichment search algorithm for [" + setNumber + "] "
            + filterSet.getName());
        gd.addCheckbox(prefix + "Save_option", settings.saveOption);
        gd.addNumericField(prefix + "Max_iterations", settings.maxIterations, 0);
        gd.addNumericField(prefix + "Converged_count", settings.convergedCount, 0);
        gd.addNumericField(prefix + "Samples", settings.enrichmentSamples, 0);
        gd.addSlider(prefix + "Fraction", 0.01, 0.99, settings.enrichmentFraction);
        gd.addSlider(prefix + "Padding", 0, 0.99, settings.enrichmentPadding);

        // Add choice of fields to optimise
        for (int i = 0; i < n; i++) {
          gd.addCheckbox(getDialogName(prefix, searchScoreFilter, i), enabled[i]);
        }

        gd.showDialog();
        runAlgorithm = !gd.wasCanceled();
      }

      if (runAlgorithm) {
        final FixedDimension[] dimensions =
            Arrays.copyOf(originalDimensions, originalDimensions.length);
        if (!nonInteractive && gd != null) {
          settings.saveOption = gd.getNextBoolean();
          settings.maxIterations = (int) gd.getNextNumber();
          settings.convergedCount = (int) gd.getNextNumber();
          settings.enrichmentSamples = (int) gd.getNextNumber();
          settings.enrichmentFraction = gd.getNextNumber();
          settings.enrichmentPadding = gd.getNextNumber();
          for (int i = 0; i < n; i++) {
            enabled[i] = gd.getNextBoolean();
          }

          searchRangeMap.put(setNumber, enabled);
        }

        for (int i = 0; i < n; i++) {
          if (!enabled[i]) {
            dimensions[i] = new FixedDimension(searchScoreFilter.getDisabledParameterValue(i));
          }
        }
        if (disabled == null) {
          disabled = new boolean[originalDimensions.length];
        }
        for (int i = 0; i < disabled.length; i++) {
          disabled[i] = !dimensions[i].active;
        }

        algorithm = Settings.EVOLVE_OPTIONS[settings.evolve];
        gaStatusPrefix = algorithm + " [" + setNumber + "] " + filterSet.getName() + " ... ";
        gaIteration = 0;
        filterScoreOptimum = null;

        final SearchSpace ss = new SearchSpace();
        ss.setTracker(this);
        // Add current optimum to seed
        if (nonInteractive && currentOptimum != null) {
          final double[][] sample = new double[1][];
          sample[0] = currentOptimum.getParameters();
          seed = merge(seed, sample);
        }
        ss.seed(seed);
        final ConvergenceChecker<FilterScore> checker =
            new InterruptConvergenceChecker(0, 0, settings.maxIterations, settings.convergedCount);

        createGaWindow();
        resumeFilterTimer();

        final SearchResult<FilterScore> optimum = ss.enrichmentSearch(dimensions, this, checker,
            settings.enrichmentSamples, settings.enrichmentFraction, settings.enrichmentPadding);

        // In case optimisation was stopped
        IJ.resetEscape();

        if (optimum != null) {
          best = ((SimpleFilterScore) optimum.getScore()).result.filter;

          // Now update the filter set for final assessment
          filterSet = new FilterSet(filterSet.getName(),
              searchSpaceToFilters((DirectFilter) best, ss.getSearchSpace()));

          // Option to save the filters
          if (settings.saveOption) {
            saveFilterSet(filterSet, setNumber, !nonInteractive);
          }
        }
      } else {
        resumeFilterTimer();
      }
    }

    IJ.showStatus("Analysing [" + setNumber + "] " + filterSet.getName() + " ...");

    // Do not support plotting if we used optimisation
    final double[] xValues = (best != null || isHeadless || (settings.plotTopN == 0)) ? null
        : new double[filterSet.size()];
    final double[] yValues = (xValues != null) ? new double[xValues.length] : null;
    SimpleFilterScore max = null;

    // Final evaluation does not need to assess all the filters if we have run optimisation.
    // It can just assess the top 1 required for the summary.
    if (best != null) {
      // Only assess the top 1 filter for the summary
      final List<Filter> list = new ArrayList<>();
      list.add((DirectFilter) best);
      filterSet = new FilterSet(filterSet.getName(), list);
    }

    // Score the filters and report the results if configured.

    final FilterScoreResult[] scoreResults =
        scoreFilters(setUncomputedStrength(filterSet), settings.showResultsTable);
    if (scoreResults == null) {
      return -1;
    }

    analysisStopWatch.stop();

    for (int index = 0; index < scoreResults.length; index++) {
      final FilterScoreResult scoreResult = scoreResults[index];

      if (xValues != null && yValues != null) {
        xValues[index] = scoreResult.filter.getNumericalValue();
        yValues[index] = scoreResult.score;
      }

      final SimpleFilterScore result =
          new SimpleFilterScore(scoreResult, allSameType, scoreResult.criteria >= minCriteria);
      if (result.compareTo(max) < 0) {
        max = result;
      }
    }

    if (max == null) {
      return -1;
    }

    if (settings.showResultsTable) {
      addToResultsWindow(scoreResults);
    }

    // Check the top filter against the limits of the original dimensions
    char[] atLimit = null;
    if (allSameType && originalDimensions != null) {
      if (disabled == null) {
        disabled = new boolean[originalDimensions.length];
      }

      final DirectFilter filter = max.result.filter;
      final int[] indices = filter.getChromosomeParameters();
      atLimit = new char[indices.length];
      final StringBuilder sb = new StringBuilder(200);
      for (int j = 0; j < indices.length; j++) {
        atLimit[j] = ComplexFilterScore.WITHIN;
        final int p = indices[j];
        if (disabled[p]) {
          continue;
        }

        final double value = filter.getParameterValue(p);
        final double lowerLimit = originalDimensions[p].getLower();
        final double upperLimit = originalDimensions[p].getUpper();

        final int c1 = Double.compare(value, lowerLimit);
        if (c1 <= 0) {
          atLimit[j] = ComplexFilterScore.FLOOR;
          sb.append(" : ").append(filter.getParameterName(p)).append(' ').append(atLimit[j])
              .append('[').append(MathUtils.rounded(value));
          if (c1 == -1) {
            atLimit[j] = ComplexFilterScore.BELOW;
            sb.append("<").append(MathUtils.rounded(lowerLimit));
          }
          sb.append("]");
        } else {
          final int c2 = Double.compare(value, upperLimit);
          if (c2 >= 0) {
            atLimit[j] = ComplexFilterScore.CEIL;
            sb.append(" : ").append(filter.getParameterName(p)).append(' ').append(atLimit[j])
                .append('[').append(MathUtils.rounded(value));
            if (c2 == 1) {
              atLimit[j] = ComplexFilterScore.ABOVE;
              sb.append(">").append(MathUtils.rounded(upperLimit));
            }
            sb.append("]");
          }
        }
      }

      if (sb.length() > 0) {
        if (max.criteriaPassed) {
          ImageJUtils.log(
              "Warning: Top filter (%s @ %s|%s) [%s] at the limit of the expanded range%s",
              filter.getName(), MathUtils.rounded((invertScore) ? -max.score : max.score),
              MathUtils.rounded((invertCriteria) ? -minCriteria : minCriteria),
              limitFailCount + fitResultData.limitRange, sb.toString());
        } else {
          ImageJUtils.log(
              "Warning: Top filter (%s @ -|%s) [%s] at the limit of the expanded range%s",
              filter.getName(), MathUtils.rounded((invertCriteria) ? -max.criteria : max.criteria),
              limitFailCount + fitResultData.limitRange, sb.toString());
        }
      }
    }

    // Note that max should never be null since this method is not run with an empty filter set

    // We may have no filters that pass the criteria
    final String type = max.result.filter.getType();
    if (!max.criteriaPassed) {
      ImageJUtils.log("Warning: Filter does not pass the criteria: %s : Best = %s using %s", type,
          MathUtils.rounded((invertCriteria) ? -max.criteria : max.criteria),
          max.result.filter.getName());
      return 0;
    }

    // This could be an option?
    final boolean allowDuplicates = true;

    // No requirement to be the same type to store for later analysis.
    // All filter sets should be able to have a best
    // filter irrespective of whether they were the same type or not.
    final ComplexFilterScore newFilterScore =
        new ComplexFilterScore(max.result, atLimit, algorithm, analysisStopWatch.getTime(), "", 0);
    addBestFilter(type, allowDuplicates, newFilterScore);

    // Add spacer at end of each result set
    if (isHeadless) {
      if (settings.showResultsTable && filterSet.size() > 1) {
        IJ.log("");
      }
    } else {
      if (settings.showResultsTable && filterSet.size() > 1) {
        resultsWindow.append("");
      }

      if (settings.plotTopN > 0 && xValues != null) {
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
        if (unique) {
          // Check the values all refer to the same property
          for (final Filter filter : filterSet.getFilters()) {
            if (!xAxisName.equals(filter.getNumericalValueName())) {
              unique = false;
              break;
            }
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
          filterAnalysisResult.plots.add(new NamedPlot(title, xAxisName, xValues, yValues));
        } else {
          plot.updateValues(xAxisName, xValues, yValues);
        }

        if (filterAnalysisResult.plots.size() > settings.plotTopN) {
          Collections.sort(filterAnalysisResult.plots, NamedPlot::compare);
          filterAnalysisResult.plots.remove(filterAnalysisResult.plots.size() - 1);
        }
      }
    }

    return 0;
  }

  /**
   * Adds the filter score for the given type to the best filters.
   *
   * <p>This method is called in filter analysis and parameter analysis to update the best result
   * for the filter type.
   *
   * @param type the type
   * @param allowDuplicates the allow duplicates
   * @param newFilterScore the new filter score
   */
  private void addBestFilter(String type, boolean allowDuplicates,
      ComplexFilterScore newFilterScore) {
    final ComplexFilterScore filterScore = filterAnalysisResult.bestFilter.get(type);
    if (filterScore != null) {
      if (allowDuplicates) {
        // Duplicate type: create a unique key
        // Start at 2 to show it is the second one of the same type
        int count = 2;
        while (filterAnalysisResult.bestFilter.containsKey(type + count)) {
          count++;
        }
        type += count;
        filterAnalysisResult.bestFilter.put(type, newFilterScore);

        // Replace (even if the same so that the latest results settings are stored)
      } else if (newFilterScore.compareTo(filterScore) <= 0) {
        filterAnalysisResult.bestFilter.put(type, newFilterScore);
      }
    } else {
      filterAnalysisResult.bestFilter.put(type, newFilterScore);
    }
  }

  private static double[][] merge(double[][] seed, double[][] sample) {
    if (seed == null) {
      seed = sample;
    } else if (sample != null) {
      // Merge
      final ArrayList<double[]> merged = new ArrayList<>(sample.length + seed.length);
      merged.addAll(Arrays.asList(seed));
      merged.addAll(Arrays.asList(sample));
      seed = merged.toArray(new double[0][]);
    }
    return seed;
  }

  /**
   * Enumerate on the min interval to convert an off grid result to one on the grid.
   *
   * @param best The optimum
   * @param stepSize Array specifying the step size for each of the parameter indices
   * @param indices Array specifying which parameter indices to search
   * @return The optimum on the min interval grid
   */
  private Chromosome<FilterScore> enumerateMinInterval(Chromosome<FilterScore> best,
      double[] stepSize, int[] indices) {
    final boolean[] enabled = new boolean[stepSize.length];
    for (int i = 0; i < indices.length; i++) {
      final int j = indices[i];
      enabled[j] = stepSize[j] > 0;
    }
    return enumerateMinInterval(best, enabled, indices);
  }

  /**
   * Enumerate on the min interval to convert an off grid result to one on the grid.
   *
   * @param best The optimum
   * @param enabled Array specifying which parameters are enabled
   * @param indices Array specifying which parameter indices to search
   * @return The optimum on the min interval grid
   */
  private Chromosome<FilterScore> enumerateMinInterval(Chromosome<FilterScore> best,
      boolean[] enabled, int[] indices) {
    // Enumerate on the min interval to produce the final filter
    searchScoreFilter = (DirectFilter) best;
    filterScoreOptimum = null;
    SearchDimension[] dimensions2 = new SearchDimension[searchScoreFilter.getNumberOfParameters()];
    for (int i = 0; i < indices.length; i++) {
      final int j = indices[i];
      if (enabled[j]) {
        final double minIncrement = searchScoreFilter.getParameterIncrement(j);
        try {
          final double value = searchScoreFilter.getParameterValue(j);
          final double max = MathUtils.ceil(value, minIncrement);
          final double min = MathUtils.floor(value, minIncrement);
          dimensions2[i] = new SearchDimension(min, max, minIncrement, 1);
          dimensions2[i].setCentre(value);
          dimensions2[i].setIncrement(minIncrement);
        } catch (final IllegalArgumentException ex) {
          IJ.error(TITLE, String.format("Unable to configure dimension [%d] %s: %s", j,
              searchScoreFilter.getParameterName(j), ex.getMessage()));
          dimensions2 = null;
          break;
        }
      }
    }
    if (dimensions2 != null) {
      // Add dimensions that have been missed
      for (int i = 0; i < dimensions2.length; i++) {
        if (dimensions2[i] == null) {
          dimensions2[i] = new SearchDimension(searchScoreFilter.getParameterValue(i));
        }
      }

      final SearchSpace ss = new SearchSpace();
      ss.setTracker(this);
      final SearchResult<FilterScore> optimum = ss.findOptimum(dimensions2, this);

      if (optimum != null) {
        best = ((SimpleFilterScore) optimum.getScore()).result.filter;
      }
    }
    return best;
  }

  /**
   * Run the optimum filter on a set of labelled peak results using various parameter settings
   * outputting performance statistics on the success of the filter to an ImageJ table.
   *
   * <p>If a new optimum is found the class level static parameters are updated.
   *
   * @param nonInteractive True if non interactive
   * @param currentOptimum the optimum
   * @param rangeReduction the range reduction
   * @return the best filter
   */
  private ComplexFilterScore parameterAnalysis(boolean nonInteractive,
      ComplexFilterScore currentOptimum, double rangeReduction) {
    this.gaResultsList = fitResultData.resultsList;
    String algorithm = "";

    // All the search algorithms use search dimensions.
    searchScoreFilter = currentOptimum.result.filter;
    final FixedDimension[] originalDimensions = new FixedDimension[3];
    double[] point = createParameters();
    final String[] names = {"Fail count", "Residuals threshold", "Duplicate distance"};

    // Start at -1 so an exception constructing the dimension can use the same index
    // to log the error.
    int index = -1;
    try {
      originalDimensions[++index] =
          new FixedDimension(settings.minFailCount, settings.maxFailCount, 1);
      // TODO - let the min intervals be configured, maybe via extra options
      if (computeDoublets) {
        originalDimensions[++index] = new FixedDimension(settings.minResidualsThreshold,
            settings.maxResidualsThreshold, 0.05);
      } else {
        originalDimensions[++index] = new FixedDimension(1, 1, 0.05);
      }
      originalDimensions[++index] =
          new FixedDimension(settings.minDuplicateDistance, settings.maxDuplicateDistance, 0.5);
    } catch (final IllegalArgumentException ex) {
      ImageJUtils.log(TITLE + " : Unable to configure dimension [%d] %s: " + ex.getMessage(), index,
          names[index]);
      return null;
    }

    // Check for a search
    boolean active = false;
    for (int i = 0; i < originalDimensions.length; i++) {
      if (originalDimensions[i].isActive()) {
        active = true;
        break;
      }
    }
    if (!active) {
      ImageJUtils.log(TITLE + " : No search range");
      return currentOptimum;
    }

    // Optionally use a reduced range (this is used for iteration)
    if (rangeReduction > 0 && rangeReduction < 1) {
      // Suppress dialogs and use the current settings
      nonInteractive = true;

      for (int i = 0; i < originalDimensions.length; i++) {
        final double centre = point[i];
        double rangeFactor = 0;
        if (originalDimensions[i].isActive()) {
          rangeFactor = (originalDimensions[i].max - originalDimensions[i].min) * rangeReduction;
        }
        final double lower = centre - rangeFactor * 0.5;
        final double upper = centre + rangeFactor * 0.5;
        originalDimensions[i] = originalDimensions[i].create(lower, upper);
      }
    }

    analysisStopWatch = StopWatch.createStarted();

    SearchResult<FilterScore> optimum = null; // Store this for later debugging

    if (settings.searchParam == 0 || settings.searchParam == 2) {
      // Collect parameters for the range search algorithm
      pauseParameterTimer();

      final boolean isStepSearch = settings.searchParam == 2;

      // The step search should use a multi-dimension refinement and no range reduction
      SearchSpace.RefinementMode myRefinementMode = SearchSpace.RefinementMode.MULTI_DIMENSION;

      GenericDialog gd = null;
      boolean runAlgorithm = nonInteractive;
      if (!nonInteractive) {
        // Ask the user for the search parameters.
        gd = new GenericDialog(TITLE);
        gd.addMessage("Configure the " + Settings.SEARCH_OPTIONS[settings.searchParam]
            + " algorithm for " + searchScoreFilter.getType());
        gd.addSlider("Width", 1, 5, settings.paRangeSearchWidth);
        if (!isStepSearch) {
          gd.addNumericField("Max_iterations", settings.paMaxIterations, 0);
          final String[] modes =
              SettingsManager.getNames((Object[]) SearchSpace.RefinementMode.values());
          gd.addSlider("Reduce", 0.01, 0.99, settings.paRangeSearchReduce);
          gd.addChoice("Refinement", modes, modes[settings.paRefinementMode]);
        }
        gd.addNumericField("Seed_size", settings.paSeedSize, 0);

        gd.showDialog();
        runAlgorithm = !gd.wasCanceled();
      }

      if (runAlgorithm) {
        final SearchDimension[] dimensions = new SearchDimension[originalDimensions.length];
        if (!nonInteractive && gd != null) {
          settings.paRangeSearchWidth = (int) gd.getNextNumber();
          if (!isStepSearch) {
            settings.paMaxIterations = (int) gd.getNextNumber();
            settings.paRangeSearchReduce = gd.getNextNumber();
            settings.paRefinementMode = gd.getNextChoiceIndex();
          }
          settings.paSeedSize = (int) gd.getNextNumber();
        }

        if (!isStepSearch) {
          myRefinementMode = SearchSpace.RefinementMode.values()[settings.paRefinementMode];
        }
        for (int i = 0; i < dimensions.length; i++) {
          if (originalDimensions[i].isActive()) {
            try {
              dimensions[i] = originalDimensions[i].create(settings.paRangeSearchWidth);
              dimensions[i].setPad(true);
              // Prevent range reduction so that the step search just does a single refinement step
              dimensions[i].setReduceFactor((isStepSearch) ? 1 : settings.paRangeSearchReduce);
              // Centre on current optimum
              dimensions[i].setCentre(point[i]);
            } catch (final IllegalArgumentException ex) {
              IJ.error(TITLE, String.format("Unable to configure dimension [%d] %s: %s", i,
                  names[i], ex.getMessage()));
              return null;
            }
          } else {
            dimensions[i] = new SearchDimension(point[i]);
          }
        }

        // Check the number of combinations is OK
        long combinations = SearchSpace.countCombinations(dimensions);
        if (!nonInteractive && combinations > 2000) {
          gd = new GenericDialog(TITLE);
          ImageJUtils.addMessage(gd,
              "%d combinations for the configured dimensions.\n \nClick 'Yes' to optimise.",
              combinations);
          gd.enableYesNoCancel();
          gd.hideCancelButton();
          gd.showDialog();
          if (!gd.wasOKed()) {
            combinations = 0;
          }
        }

        if (combinations == 0) {
          resumeParameterTimer();
        } else {
          algorithm =
              Settings.SEARCH_OPTIONS[settings.searchParam] + " " + settings.paRangeSearchWidth;
          gaStatusPrefix = algorithm + " " + searchScoreFilter.getName() + " ... ";
          gaIteration = 0;
          parameterScoreOptimum = null;

          final SearchSpace ss = new SearchSpace();
          ss.setTracker(this);
          if (settings.paSeedSize > 0) {
            // Add current optimum to seed
            // Note: If we have an optimum and we are not seeding this should not matter as the
            // dimensions have been centred on the current optimum
            final double[][] seed = new double[1][];
            seed[0] = point;
            // Sample without rounding as the seed will be rounded
            final double[][] sample =
                SearchSpace.sampleWithoutRounding(dimensions, settings.paSeedSize - 1, null);
            ss.seed(merge(seed, sample));
          }
          final ConvergenceChecker<FilterScore> checker =
              new InterruptConvergenceChecker(0, 0, settings.paMaxIterations);

          createGaWindow();
          resumeParameterTimer();

          optimum = ss.search(dimensions, new ParameterScoreFunction(), checker, myRefinementMode);

          // In case optimisation was stopped
          IJ.resetEscape();

          if (optimum != null) {
            // Now update the parameters for final assessment
            point = optimum.getPoint();
          }
        }
      } else {
        resumeParameterTimer();
      }
    }

    if (settings.searchParam == 1) {
      // Collect parameters for the enrichment search algorithm
      pauseParameterTimer();

      GenericDialog gd = null;
      boolean runAlgorithm = nonInteractive;
      if (!nonInteractive) {
        // Ask the user for the search parameters.
        gd = new GenericDialog(TITLE);
        gd.addMessage("Configure the " + Settings.SEARCH_OPTIONS[settings.searchParam]
            + " algorithm for " + searchScoreFilter.getType());
        gd.addNumericField("Max_iterations", settings.paMaxIterations, 0);
        gd.addNumericField("Converged_count", settings.paConvergedCount, 0);
        gd.addNumericField("Samples", settings.paEnrichmentSamples, 0);
        gd.addSlider("Fraction", 0.01, 0.99, settings.paEnrichmentFraction);
        gd.addSlider("Padding", 0, 0.99, settings.paEnrichmentPadding);

        gd.showDialog();
        runAlgorithm = !gd.wasCanceled();
      }

      if (runAlgorithm) {
        final FixedDimension[] dimensions =
            Arrays.copyOf(originalDimensions, originalDimensions.length);
        if (!nonInteractive && gd != null) {
          settings.paMaxIterations = (int) gd.getNextNumber();
          settings.paConvergedCount = (int) gd.getNextNumber();
          settings.paEnrichmentSamples = (int) gd.getNextNumber();
          settings.paEnrichmentFraction = gd.getNextNumber();
          settings.paEnrichmentPadding = gd.getNextNumber();
        }

        algorithm = Settings.SEARCH_OPTIONS[settings.searchParam];
        gaStatusPrefix = algorithm + " " + searchScoreFilter.getName() + " ... ";
        gaIteration = 0;
        parameterScoreOptimum = null;

        final SearchSpace ss = new SearchSpace();
        ss.setTracker(this);
        // Add current optimum to seed
        final double[][] seed = new double[1][];
        seed[0] = point;
        ss.seed(seed);
        final ConvergenceChecker<FilterScore> checker = new InterruptConvergenceChecker(0, 0,
            settings.paMaxIterations, settings.paConvergedCount);

        createGaWindow();
        resumeParameterTimer();

        optimum = ss.enrichmentSearch(dimensions, new ParameterScoreFunction(), checker,
            settings.paEnrichmentSamples, settings.paEnrichmentFraction,
            settings.paEnrichmentPadding);

        // In case optimisation was stopped
        IJ.resetEscape();

        if (optimum != null) {
          point = optimum.getPoint();
        }
      } else {
        resumeParameterTimer();
      }
    }

    if (settings.searchParam == 3) {
      // Collect parameters for the enumeration search algorithm
      pauseParameterTimer();

      final SearchDimension[] dimensions = new SearchDimension[originalDimensions.length];
      for (int i = 0; i < dimensions.length; i++) {
        if (originalDimensions[i].isActive()) {
          try {
            dimensions[i] = originalDimensions[i].create(0);
          } catch (final IllegalArgumentException ex) {
            IJ.error(TITLE, String.format("Unable to configure dimension [%d] %s: %s", i, names[i],
                ex.getMessage()));
            return null;
          }
        } else {
          dimensions[i] = new SearchDimension(point[i]);
        }
      }

      GenericDialog gd = null;
      long combinations = SearchSpace.countCombinations(dimensions);
      if (!nonInteractive && combinations > 2000) {
        gd = new GenericDialog(TITLE);
        ImageJUtils.addMessage(gd,
            "%d combinations for the configured dimensions.\n \nClick 'Yes' to optimise.",
            combinations);
        gd.enableYesNoCancel();
        gd.hideCancelButton();
        gd.showDialog();
        if (!gd.wasOKed()) {
          combinations = 0;
        }
      }

      if (combinations == 0) {
        resumeParameterTimer();
      } else {
        algorithm = Settings.SEARCH_OPTIONS[settings.searchParam];
        gaStatusPrefix = algorithm + " " + searchScoreFilter.getName() + " ... ";
        gaIteration = 0;
        parameterScoreOptimum = null;

        final SearchSpace ss = new SearchSpace();

        ss.setTracker(this);
        createGaWindow();
        resumeParameterTimer();

        optimum = ss.findOptimum(dimensions, new ParameterScoreFunction());

        // In case optimisation was stopped
        IJ.resetEscape();

        if (optimum != null) {
          // Now update the parameters for final assessment
          point = optimum.getPoint();
        }
      }
    }

    IJ.showStatus("Analysing " + searchScoreFilter.getName() + " ...");

    // Update the parameters using the optimum
    settings.failCount = (int) Math.round(point[0]);
    residualsThreshold = settings.residualsThreshold = point[1];
    settings.duplicateDistance = point[2];
    // Refresh the coordinate store
    if (coordinateStore == null
        // Due to the scaling factor the distance may not be exactly the same
        || DoubleEquality.relativeError(settings.duplicateDistance,
            coordinateStore.getXyResolution() / distanceScallingFactor) > 0.01) {
      coordinateStore = createCoordinateStore();
    }
    createResultsPrefix2();

    // (Re) Score the filter.

    // TODO - check this is now OK. Maybe remove the enumeration on the min interval grid

    // If scoring of filter here is different to scoring in the optimisation routine it is probably
    // an ss_filter.clone() issue,
    // i.e. multi-threading use of the filter clone is not working.
    // Or it could be that the optimisation produced params off the min-interval grid
    final FilterScoreResult scoreResult = scoreFilter(searchScoreFilter);
    if (optimum != null && (scoreResult.score != optimum.getScore().score
        && scoreResult.criteria != optimum.getScore().criteria)) {
      final ParameterScoreResult r = scoreFilter((DirectFilter) searchScoreFilter.clone(),
          defaultMinimalFilter, settings.failCount, residualsThreshold, settings.duplicateDistance,
          createCoordinateStore(settings.duplicateDistance), false);

      ImageJUtils.log("Weird re-score of the filter: %f!=%f or %f!=%f (%f:%f)", scoreResult.score,
          optimum.getScore().score, scoreResult.criteria, optimum.getScore().criteria, r.score,
          r.criteria);
    }

    final SimpleFilterScore max =
        new SimpleFilterScore(scoreResult, true, scoreResult.criteria >= minCriteria);

    analysisStopWatch.stop();

    if (settings.showResultsTable) {
      addToResultsWindow(scoreResult.text);
    }

    // Check the top result against the limits of the original dimensions
    final StringBuilder sb = new StringBuilder(200);
    for (int j = 0; j < originalDimensions.length; j++) {
      if (!originalDimensions[j].isActive()) {
        continue;
      }

      final double value = point[j];
      final double lowerLimit = originalDimensions[j].getLower();
      final double upperLimit = originalDimensions[j].getUpper();

      final int c1 = Double.compare(value, lowerLimit);
      if (c1 <= 0) {
        sb.append(" : ").append(names[j]).append(' ').append(ComplexFilterScore.FLOOR).append('[')
            .append(MathUtils.rounded(value));
        if (c1 == -1) {
          sb.append("<").append(MathUtils.rounded(lowerLimit));
        }
        sb.append("]");
      } else {
        final int c2 = Double.compare(value, upperLimit);
        if (c2 >= 0) {
          sb.append(" : ").append(names[j]).append(' ').append(ComplexFilterScore.CEIL).append('[')
              .append(MathUtils.rounded(value));
          if (c2 == 1) {
            sb.append(">").append(MathUtils.rounded(upperLimit));
          }
          sb.append("]");
        }
      }
    }

    if (sb.length() > 0) {
      if (max.criteriaPassed) {
        ImageJUtils.log(
            "Warning: Top filter (%s @ %s|%s) [%s] at the limit of the expanded range%s",
            searchScoreFilter.getName(), MathUtils.rounded((invertScore) ? -max.score : max.score),
            MathUtils.rounded((invertCriteria) ? -minCriteria : minCriteria),
            limitFailCount + fitResultData.limitRange, sb.toString());
      } else {
        ImageJUtils.log("Warning: Top filter (%s @ -|%s) [%s] at the limit of the expanded range%s",
            searchScoreFilter.getName(),
            MathUtils.rounded((invertCriteria) ? -max.criteria : max.criteria),
            limitFailCount + fitResultData.limitRange, sb.toString());
      }
    }

    // We may have no filters that pass the criteria
    final String type = max.result.filter.getType();
    if (!max.criteriaPassed) {
      ImageJUtils.log("Warning: Filter does not pass the criteria: %s : Best = %s using %s", type,
          MathUtils.rounded((invertCriteria) ? -max.criteria : max.criteria),
          max.result.filter.getName());
      return null;
    }

    // Update without duplicates
    final boolean allowDuplicates = false;
    // Re-use the atLimit and algorithm for the input optimum
    final ComplexFilterScore newFilterScore =
        new ComplexFilterScore(max.result, currentOptimum.atLimit, currentOptimum.algorithm,
            currentOptimum.time, algorithm, analysisStopWatch.getTime());
    addBestFilter(type, allowDuplicates, newFilterScore);

    // Add spacer at end of each result set
    if (settings.showResultsTable) {
      addToResultsWindow("");
    }

    if (newFilterScore.compareTo(currentOptimum) <= 0) {
      return newFilterScore;
    }

    // Update the algorithm and time
    currentOptimum.paramAlgorithm = algorithm;
    currentOptimum.paramTime = analysisStopWatch.getTime();
    return currentOptimum;
  }

  private void pauseFilterTimer() {
    filterAnalysisStopWatch.suspend();
    analysisStopWatch.suspend();
    if (iterationStopWatch != null) {
      iterationStopWatch.suspend();
    }
  }

  private void resumeFilterTimer() {
    filterAnalysisStopWatch.resume();
    analysisStopWatch.resume();
    if (iterationStopWatch != null) {
      iterationStopWatch.resume();
    }
  }

  private void pauseParameterTimer() {
    parameterAnalysisStopWatch.suspend();
    analysisStopWatch.suspend();
    if (iterationStopWatch != null) {
      iterationStopWatch.suspend();
    }
  }

  private void resumeParameterTimer() {
    parameterAnalysisStopWatch.resume();
    analysisStopWatch.resume();
    if (iterationStopWatch != null) {
      iterationStopWatch.resume();
    }
  }


  private void addToResultsWindow(final String text) {
    if (text != null) {
      if (resultsWindow != null) {
        resultsWindow.append(text);
      } else {
        IJ.log(text);
      }
    }
  }

  private void addToResultsWindow(FilterScoreResult[] scoreResults) {
    if (resultsWindow != null) {
      try (BufferedTextWindow tw = new BufferedTextWindow(resultsWindow)) {
        tw.setIncrement(0);
        for (int index = 0; index < scoreResults.length; index++) {
          addToResultsOutput(tw::append, scoreResults[index].text);
        }
      }
    } else {
      for (int index = 0; index < scoreResults.length; index++) {
        addToResultsOutput(IJ::log, scoreResults[index].text);
      }
    }
  }

  private static void addToResultsOutput(Consumer<String> output, final String text) {
    if (text != null) {
      output.accept(text);
    }
  }

  /**
   * When the two filters have equal scores, select the filter using the filter parameters.
   *
   * @param filter1 the filter 1
   * @param filter2 the filter 2
   * @return The chosen filter (the one with the strongest parameters)
   */
  public Filter selectFilter(Filter filter1, Filter filter2) {
    return (filter2.weakest(filter1) < 0) ? filter1 : filter2;
  }

  private static String getDialogName(String prefix, Filter filter, int index) {
    final ParameterType type = filter.getParameterType(index);
    final String parameterName = prefix + type.getShortName().replace(" ", "_");
    if (type.getShortName().equals(type.getName())) {
      return parameterName;
    }
    return parameterName + " (" + type.getName() + ")";
  }

  private double getCriteria(FractionClassificationResult result) {
    return getScore(result, settings.criteriaIndex, invertCriteria);
  }

  private double getScore(FractionClassificationResult result) {
    return getScore(result, settings.scoreIndex, invertScore);
  }

  private double getScore(FractionClassificationResult result, final int index,
      final boolean invert) {
    final double score = getScore(result, index);
    return (invert) ? -score : score;
  }

  private double getScore(FractionClassificationResult result, final int index) {
    // This order must match the COLUMNS order
    switch (index) {
      case 0:
        return result.getNumberOfPositives();
      case 1:
        return result.getNumberOfNegatives();
      case 2:
        return (double) fitResultData.countActual - result.getNumberOfPositives();
      case 3:
        return createIntegerResult(result).getPrecision();
      case 4:
        return createIntegerResult(result).getRecall();
      case 5:
        return createIntegerResult(result).getF1Score();
      case 6:
        return createIntegerResult(result).getJaccard();
      case 7:
        return result.getTruePositives();
      case 8:
        return result.getFalsePositives();
      case 9:
        return result.getFalseNegatives();
      case 10:
        return result.getPrecision();
      case 11:
        return result.getRecall();
      case 12:
        return result.getF1Score();
      case 13:
        return result.getJaccard();
      default:
        return 0;
    }
  }

  private static boolean requiresInversion(final int index) {
    switch (index) {
      case 1: // FP
      case 2: // FN
      case 8: // fFP
      case 9: // fFN
        return true;

      default:
        return false;
    }
  }

  private NamedPlot getNamedPlot(String title) {
    for (final NamedPlot p : filterAnalysisResult.plots) {
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

  /**
   * Score the filter using the results list and the configured fail count.
   *
   * @param filter the filter
   * @param minFilter the min filter
   * @param resultsList the results list
   * @param coordinateStore the coordinate store
   * @return The score
   */
  private FractionClassificationResult scoreFilter(DirectFilter filter, DirectFilter minFilter,
      MultiPathFitResults[] resultsList, CoordinateStore coordinateStore) {
    return scoreFilter(filter, minFilter, resultsList, null, coordinateStore);
  }

  private FractionScoreStore scoreStore;

  /**
   * Score the filter using the results list and the configured fail count.
   *
   * @param filter the filter
   * @param minFilter the min filter
   * @param resultsList the results list
   * @param allAssignments all the assignments
   * @param coordinateStore the coordinate store
   * @return The score
   */
  private FractionClassificationResult scoreFilter(DirectFilter filter, DirectFilter minFilter,
      MultiPathFitResults[] resultsList, List<FractionalAssignment[]> allAssignments,
      CoordinateStore coordinateStore) {
    final MultiPathFilter multiPathFilter = createMpf(filter, minFilter);
    // Note: We always use the subset method since fail counts have been accumulated when we read in
    // the results.
    return multiPathFilter.fractionScoreSubset(resultsList,
        createFailCounter(settings.failCount), fitResultData.countActual,
        allAssignments, scoreStore, coordinateStore);
  }

  private FilterScoreResult scoreFilter(DirectFilter filter, DirectFilter minFilter,
      boolean createTextResult, CoordinateStore coordinateStore) {
    final FractionClassificationResult r =
        scoreFilter(filter, minFilter, gaResultsListToScore, coordinateStore);
    final double score = getScore(r);
    final double criteria = getCriteria(r);

    // Show the result if it achieves the criteria limit
    final String text =
        (createTextResult && criteria >= minCriteria) ? createResult(filter, r).toString() : null;

    return new FilterScoreResult(score, criteria, filter, text);
  }

  private FilterScoreResult scoreFilter(DirectFilter filter) {
    final FractionClassificationResult r =
        scoreFilter(filter, defaultMinimalFilter, fitResultData.resultsList, coordinateStore);
    final double score = getScore(r);
    final double criteria = getCriteria(r);
    // Create the score output
    final String text = createResult(filter, r).toString();
    return new FilterScoreResult(score, criteria, filter, text);
  }

  private ParameterScoreResult scoreFilter(DirectFilter filter, DirectFilter minFilter,
      int failCount, double residualsThreshold, double duplicateDistance,
      CoordinateStore coordinateStore, boolean createTextResult) {
    final MultiPathFilter multiPathFilter =
        new MultiPathFilter(filter, minFilter, residualsThreshold);

    final FractionClassificationResult r = multiPathFilter.fractionScoreSubset(gaResultsListToScore,
        createFailCounter(failCount), fitResultData.countActual, null, null,
        coordinateStore);

    final double score = getScore(r);
    final double criteria = getCriteria(r);
    // Create the score output
    final String text = (createTextResult && criteria >= minCriteria)
        ? createResult(filter, r,
            buildResultsPrefix2(failCount, residualsThreshold, duplicateDistance)).toString()
        : null;
    final double[] parameters = new double[] {failCount, residualsThreshold, duplicateDistance};
    return new ParameterScoreResult(score, criteria, parameters, text);
  }

  private MultiPathFilter createMpf(DirectFilter filter, DirectFilter minFilter) {
    return new MultiPathFilter(filter, minFilter, residualsThreshold);
  }

  /**
   * Score the filter using the results list and the configured fail count.
   *
   * @param filter the filter
   * @return The score
   */
  private ArrayList<FractionalAssignment[]> getAssignments(DirectFilter filter) {
    final MultiPathFilter multiPathFilter = createMpf(filter, defaultMinimalFilter);

    final ArrayList<FractionalAssignment[]> allAssignments =
        new ArrayList<>(fitResultData.resultsList.length);
    multiPathFilter.fractionScoreSubset(fitResultData.resultsList,
        createFailCounter(settings.failCount), fitResultData.countActual,
        allAssignments, null, coordinateStore);
    return allAssignments;
  }

  private StringBuilder createResult(DirectFilter filter, FractionClassificationResult result) {
    return createResult(filter, result, resultsPrefix2);
  }

  private StringBuilder createResult(DirectFilter filter, FractionClassificationResult result,
      String resultsPrefix2) {
    final StringBuilder sb = new StringBuilder(resultsPrefix);
    sb.append(filter.getName()).append(resultsPrefix2).append(fitResultData.resultsPrefix3);

    int index = 0;

    // TODO - Fix the scores that we show since we no longer have TN results
    // We could set:
    // TN as any candidate that does not match a true result.
    // FN as any candidate that does match a true result (that is not matched by any fit result)
    // To do this properly would require that we store all the matches of candidates to the data.
    // These can then be totalled up given the candidates that have not been used to create a
    // positive.

    // Integer results
    if (settings.requireIntegerResults) {
      final ClassificationResult r2 = createIntegerResult(result);
      add(sb, r2.getTruePositives(), index++);
      add(sb, r2.getFalsePositives(), index++);
      add(sb, r2.getFalseNegatives(), index++);
      add(sb, r2.getPrecision(), index++);
      add(sb, r2.getRecall(), index++);
      add(sb, r2.getF1Score(), index++);
      add(sb, r2.getJaccard(), index++);
    } else {
      index += 7;
    }

    addCount(sb, result.getTruePositives(), index++);
    addCount(sb, result.getFalsePositives(), index++);
    addCount(sb, result.getFalseNegatives(), index++);
    add(sb, result.getPrecision(), index++);
    add(sb, result.getRecall(), index++);
    add(sb, result.getF1Score(), index++);
    add(sb, result.getJaccard(), index);

    return sb;
  }

  private ClassificationResult createIntegerResult(FractionClassificationResult result) {
    return new ClassificationResult(result.getNumberOfPositives(), result.getNumberOfNegatives(), 0,
        fitResultData.countActual - result.getNumberOfPositives());
  }

  /**
   * Gets the original score.
   *
   * @param result the result
   * @return the original score
   */
  private FractionClassificationResult getOriginalScore(FractionClassificationResult result) {
    // TODO: There was a note to fix this function.
    // Determine if the code below is not adequate.

    // Score the fitting results against the original simulated data:
    // TP are all fit results that can be matched to a spot
    // FP are all fit results that cannot be matched to a spot
    // FN are the number of missed spots
    // Note: We cannot calculate TN since this is the number of fit candidates that are
    // filtered after fitting that do not match a spot or were not fitted.
    final double fp = result.getPositives() - result.getTruePositives();
    final double fn = fitResultData.countActual - result.getTruePositives();
    return new FractionClassificationResult(result.getTruePositives(), fp, 0, fn);
  }

  private static void add(StringBuilder sb, int value) {
    sb.append('\t').append(value);
  }

  private static void add(StringBuilder sb, String value) {
    sb.append('\t').append(value);
  }

  private void add(StringBuilder sb, int value, int index) {
    if (settings.showColumns[index]) {
      add(sb, value);
    }
  }

  private void add(StringBuilder sb, double value, int index) {
    if (settings.showColumns[index]) {
      add(sb, MathUtils.rounded(value));
    }
  }

  private void addCount(StringBuilder sb, double value, int index) {
    if (settings.showColumns[index]) {
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
  }

  private void saveFilter(DirectFilter filter) {
    // Save the filter to file
    final String filename = getFilename("Best_Filter_File", settings.filterFilename);
    if (filename != null) {
      settings.filterFilename = filename;
      Prefs.set(Settings.KEY_FILTER_FILENAME, filename);

      final List<Filter> filters = new ArrayList<>(1);
      filters.add(filter);
      final FilterSet filterSet = new FilterSet(filter.getName(), filters);
      saveFilterSet(filterSet, filename);
    }
  }

  /**
   * Gets the filename (with a .json extension).
   *
   * @param title the title
   * @param filename the filename
   * @return the filename
   */
  static String getFilename(String title, String filename) {
    filename = ImageJUtils.getFilename(title, filename);
    // Use JSON extension
    if (filename != null) {
      filename = FileUtils.replaceExtension(filename, ".json");
    }
    return filename;
  }

  private static void saveFilterSet(FilterSet filterSet, String filename) {
    try (OutputStreamWriter out = new OutputStreamWriter(new FileOutputStream(filename), "UTF-8")) {
      final List<FilterSet> list = new ArrayList<>(1);
      list.add(filterSet);
      // Use the instance so we can catch the exception
      out.write(FilterXStreamUtils.getXStreamInstance().toXML(list));
    } catch (final Exception ex) {
      IJ.log("Unable to save the filter sets to file: " + ex.getMessage());
    }
  }

  /**
   * Save the filter set to a file prompted from the user.
   *
   * <p>If non-interactive then the last filename will be used. This will overwrite existing files
   * if multiple filter sets are used.
   *
   * @param filterSet the filter set
   * @param setNumber the set number
   * @param interactive Set to true to prompt the user
   */
  private void saveFilterSet(FilterSet filterSet, int setNumber, boolean interactive) {
    pauseFilterTimer();

    if (interactive) {
      final String filename = getFilename("Filter_set_" + setNumber, settings.filterSetFilename);
      if (filename != null) {
        settings.filterSetFilename = filename;
        Prefs.set(Settings.KEY_FILTERSET_FILENAME, filename);
        saveFilterSet(filterSet, filename);
      }
    } else {
      // Add support for multiple filter sets
      saveFilterSet(filterSet, settings.filterSetFilename);
    }

    resumeFilterTimer();
  }

  /**
   * Save PeakFit configuration template using the current benchmark settings.
   *
   * @param topFilterSummary the top filter summary
   */
  private void saveTemplate(String topFilterSummary) {
    final FitEngineConfiguration config = new FitEngineConfiguration();
    if (!updateAllConfiguration(config, true)) {
      IJ.log("Unable to create the template configuration");
      return;
    }

    // Remove the PSF width to make the template generic
    config.getFitConfiguration().setInitialPeakStdDev(0);

    // Only get this once when doing iterative analysis
    String filename;
    final boolean localSaveTemplateIsSet = saveTemplateIsSet;
    if (localSaveTemplateIsSet) {
      filename = settings.templateFilename;
    } else {
      filename = getFilename("Template_File", settings.templateFilename);
      saveTemplateIsSet = true;
    }
    if (filename != null) {
      settings.templateFilename = filename;
      Prefs.set(Settings.KEY_TEMPLATE_FILENAME, filename);
      final TemplateSettings.Builder templateSettings = TemplateSettings.newBuilder();
      getNotes(templateSettings, topFilterSummary);
      templateSettings.setFitEngineSettings(config.getFitEngineSettings());
      if (!SettingsManager.toJson(templateSettings.build(), filename,
          SettingsManager.FLAG_SILENT | SettingsManager.FLAG_JSON_WHITESPACE)) {
        IJ.log("Unable to save the template configuration");
        return;
      }

      // The rest of the code below extracts an example image for the template.
      // This need only be performed once as the sample image is the same for all iterations.
      if (localSaveTemplateIsSet) {
        return;
      }

      // Save some random frames from the test image data
      final ImagePlus imp = CreateData.getImage();
      if (imp == null) {
        return;
      }

      // Get the number of frames
      final ResultsImageSampler sampler = getSampler(results, imp);
      if (!sampler.isValid()) {
        return;
      }

      // Iteratively show the example until the user is happy.
      // Yes = OK, No = Repeat, Cancel = Do not save
      final String keyNo = "nNo";
      final String keyLow = "nLower";
      final String keyHigh = "nHigher";
      if (ImageJUtils.isMacro()) {
        // Collect the options if running in a macro
        final String options = Macro.getOptions();
        settings.countNo =
            Integer.parseInt(Macro.getValue(options, keyNo, Integer.toString(settings.countNo)));
        settings.countLow =
            Integer.parseInt(Macro.getValue(options, keyLow, Integer.toString(settings.countLow)));
        settings.countHigh = Integer
            .parseInt(Macro.getValue(options, keyHigh, Integer.toString(settings.countHigh)));
      } else if (settings.countLow + settings.countHigh == 0) {
        settings.countLow = settings.countHigh = 1;
      }

      final ImagePlus[] out = new ImagePlus[1];
      out[0] = sampler.getSample(settings.countNo, settings.countLow, settings.countHigh);

      if (!ImageJUtils.isMacro()) {
        // Show the template results
        final ConfigurationTemplate configTemplate = new ConfigurationTemplate();

        // Interactively show the sample image data
        final boolean[] close = new boolean[1];
        final ImagePlus[] outImp = new ImagePlus[1];
        if (out[0] != null) {
          final WindowOrganiser windowOrganiser = new WindowOrganiser();
          outImp[0] = display(out[0], windowOrganiser);
          if (windowOrganiser.isNotEmpty()) {
            close[0] = true;
            // Zoom a bit
            final ImageWindow iw = outImp[0].getWindow();
            for (int i = 7; i-- > 0 && Math.max(iw.getWidth(), iw.getHeight()) < 512;) {
              iw.getCanvas().zoomIn(0, 0);
            }
          }
          configTemplate.createResults(outImp[0]);
        }

        // TODO - fix this when a second sample is made as the results are not updated.

        final ImageListener listener = new ImageListener() {
          @Override
          public void imageOpened(ImagePlus imp) {
            // Do nothing
          }

          @Override
          public void imageClosed(ImagePlus imp) {
            // Do nothing
          }

          @Override
          public void imageUpdated(ImagePlus imp) {
            if (imp != null && imp == outImp[0]) {
              configTemplate.updateResults(imp.getCurrentSlice());
            }
          }
        };
        ImagePlus.addImageListener(listener);

        // Turn off the recorder when the dialog is showing
        final boolean record = Recorder.record;
        Recorder.record = false;
        final NonBlockingGenericDialog gd = new NonBlockingGenericDialog(TITLE);
        ImageJUtils.addMessage(gd,
            "Showing image data for the template example.\n \nSample Frames:\nEmpty = %d\n"
                + "Lower density = %d\nHigher density = %d\n",
            sampler.getNumberOfEmptySamples(), sampler.getNumberOfLowDensitySamples(),
            sampler.getNumberOfHighDensitySamples());
        gd.addSlider(keyNo, 0, 10, settings.countNo);
        gd.addSlider(keyLow, 0, 10, settings.countLow);
        gd.addSlider(keyHigh, 0, 10, settings.countHigh);
        gd.addDialogListener((genDialog, event) -> {
          // If the event is null then this is the final call when the
          // dialog has been closed. We ignore this to prevent generating a
          // image the user has not seen.
          if (event == null) {
            return true;
          }
          settings.countNo = (int) genDialog.getNextNumber();
          settings.countLow = (int) genDialog.getNextNumber();
          settings.countHigh = (int) genDialog.getNextNumber();
          out[0] = sampler.getSample(settings.countNo, settings.countLow, settings.countHigh);
          if (out[0] != null) {
            final WindowOrganiser windowOrganiser = new WindowOrganiser();
            outImp[0] = display(out[0], windowOrganiser);
            if (windowOrganiser.isNotEmpty()) {
              close[0] = true;
              // Zoom a bit
              final ImageWindow iw = outImp[0].getWindow();
              for (int i = 7; i-- > 0 && Math.max(iw.getWidth(), iw.getHeight()) < 512;) {
                iw.getCanvas().zoomIn(0, 0);
              }
            }
            configTemplate.createResults(outImp[0]);
          }
          return true;
        });
        gd.showDialog();
        if (gd.wasCanceled()) {
          out[0] = null;
          settings.countNo = settings.countLow = settings.countHigh = 0; // For the recorder
        }
        if (close[0]) {
          // Because closing the image sets the stack pixels array to null
          if (out[0] != null) {
            out[0] = out[0].duplicate();
          }
          outImp[0].close();
        }
        configTemplate.closeResults();
        ImagePlus.removeImageListener(listener);

        if (record) {
          Recorder.record = true;
          Recorder.recordOption(keyNo, Integer.toString(settings.countNo));
          Recorder.recordOption(keyLow, Integer.toString(settings.countLow));
          Recorder.recordOption(keyHigh, Integer.toString(settings.countHigh));
        }
      }

      if (out[0] == null) {
        return;
      }

      final ImagePlus example = out[0];
      filename = FileUtils.replaceExtension(filename, ".tif");
      IJ.save(example, filename);
    }
  }

  private static synchronized ResultsImageSampler getSampler(MemoryPeakResults results,
      ImagePlus imp) {
    Pair<Integer, ResultsImageSampler> sampler = samplerRef.get();
    if (sampler == null || imp.getID() != sampler.getKey()
        || sampler.getValue().getResults() != results) {
      final ResultsImageSampler imageSampler =
          new ResultsImageSampler(results, imp.getImageStack(), 32);
      imageSampler.analyse();
      sampler = Pair.of(imp.getID(), imageSampler);
      samplerRef.set(sampler);
    }
    return sampler.getValue();
  }

  private static ImagePlus display(ImagePlus example, WindowOrganiser windowOrganiser) {
    final String title = "Template Example";
    // Update the example image
    example.setTitle(title);

    // Display as a new image. This is so we can close it later.
    return ConfigurationTemplate.displayTemplate(title, example, windowOrganiser);
  }

  private void getNotes(TemplateSettings.Builder templateSettings, String topFilterSummary) {
    templateSettings.addNotes("Benchmark template");
    if (!TextUtils.isNullOrEmpty(settings.resultsTitle)) {
      addField(templateSettings, "Filter Analysis Title", settings.resultsTitle);
    }
    // Add create data settings.
    // Just add the columns and the data from the summary window
    final String header = createResultsHeader(true);
    addField(templateSettings, "Filter Analysis Summary Fields", header);
    addField(templateSettings, "Filter Analysis Summary Values", topFilterSummary);
    // Now pick out key values...
    addKeyFields(templateSettings, header, topFilterSummary, new String[] {"Molecules", "Density",
        "SNR", "s (nm)", "a (nm)", "Lower D", "Upper D", "Lower factor", "Upper factor"});

    // Add any other settings that may be useful in the template
    addField(templateSettings, "Created", getCurrentTimeStamp());
  }

  /**
   * Adds the field to the template settings.
   *
   * @param settings the settings
   * @param field the field
   * @param value the value
   */
  static void addField(TemplateSettings.Builder settings, String field, String value) {
    settings.addNotes(field + ": " + value);
  }

  /**
   * Adds the key fields to the template settings.
   *
   * @param settings the settings
   * @param header the header
   * @param summary the summary
   * @param fields the fields
   */
  static void addKeyFields(TemplateSettings.Builder settings, String header, String summary,
      String[] fields) {
    final String[] labels = header.split("\t");
    final String[] values = summary.split("\t");
    for (final String field : fields) {
      for (int i = 0; i < labels.length; i++) {
        if (labels[i].startsWith(field)) {
          addField(settings, labels[i], values[i]);
          break;
        }
      }
    }
  }

  /**
   * Gets the current time stamp.
   *
   * @return the current time stamp
   */
  static String getCurrentTimeStamp() {
    final SimpleDateFormat sdfDate = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
    final Date now = new Date();
    return sdfDate.format(now);
  }

  /**
   * Depth analysis.
   *
   * @param allAssignments The assignments generated from running the filter (or null)
   * @param filter the filter
   * @return the assignments
   */
  @Nullable
  private ArrayList<FractionalAssignment[]>
      depthAnalysis(ArrayList<FractionalAssignment[]> allAssignments, DirectFilter filter) {
    // TODO : This analysis ignores the partial match distance.
    // Use the score for each result to get a weighted histogram.

    if (!settings.depthRecallAnalysis || simulationParameters.fixedDepth) {
      return null;
    }

    // Build a histogram of the number of spots at different depths
    final double[] depths = fitResultData.depthStats.getValues();
    double[] limits = MathUtils.limits(depths);

    final int bins = HistogramPlot.getBinsSqrtRule(depths.length);
    final double[][] h1 = HistogramPlot.calcHistogram(depths, limits[0], limits[1], bins);
    final double[][] h2 = HistogramPlot.calcHistogram(fitResultData.depthFitStats.getValues(),
        limits[0], limits[1], bins);

    // To get the number of TP at each depth will require that the filter is run
    // manually to get the results that pass.
    if (allAssignments == null) {
      allAssignments = getAssignments(filter);
    }

    double[] depths2 = new double[results.size()];
    int count = 0;
    for (final FractionalAssignment[] assignments : allAssignments) {
      if (assignments == null) {
        continue;
      }
      for (int i = 0; i < assignments.length; i++) {
        final CustomFractionalAssignment c = (CustomFractionalAssignment) assignments[i];
        depths2[count++] = c.peak.getZPosition();
      }
    }
    depths2 = Arrays.copyOf(depths2, count);

    // Build a histogram using the same limits
    final double[][] h3 = HistogramPlot.calcHistogram(depths2, limits[0], limits[1], bins);

    // Convert pixel depth to nm
    for (int i = 0; i < h1[0].length; i++) {
      h1[0][i] *= simulationParameters.pixelPitch;
    }
    limits[0] *= simulationParameters.pixelPitch;
    limits[1] *= simulationParameters.pixelPitch;

    // Produce a histogram of the number of spots at each depth
    final String title1 = TITLE + " Depth Histogram";
    final Plot plot1 = new Plot(title1, "Depth (nm)", "Frequency");
    plot1.setLimits(limits[0], limits[1], 0, MathUtils.max(h1[1]));
    plot1.setColor(Color.black);
    plot1.addPoints(h1[0], h1[1], Plot.BAR);
    plot1.addLabel(0, 0, "Black = Spots; Blue = Fitted; Red = Filtered");
    plot1.setColor(Color.blue);
    plot1.addPoints(h1[0], h2[1], Plot.BAR);
    plot1.setColor(Color.red);
    plot1.addPoints(h1[0], h3[1], Plot.BAR);
    plot1.setColor(Color.magenta);
    ImageJUtils.display(title1, plot1, wo);

    // Interpolate
    final double halfBinWidth = (h1[0][1] - h1[0][0]) * 0.5;
    // Remove final value of the histogram as this is at the upper limit of the range (i.e. count
    // zero)
    h1[0] = Arrays.copyOf(h1[0], h1[0].length - 1);
    h1[1] = Arrays.copyOf(h1[1], h1[0].length);
    h2[1] = Arrays.copyOf(h2[1], h1[0].length);
    h3[1] = Arrays.copyOf(h3[1], h1[0].length);

    // TODO : Fix the smoothing since LOESS sometimes does not work.
    // Perhaps allow configuration of the number of histogram bins and the smoothing bandwidth.

    // Use minimum of 3 points for smoothing
    // Ensure we use at least x% of data
    final double bandwidth = Math.max(3.0 / h1[0].length, 0.15);
    final LoessInterpolator loess = new LoessInterpolator(bandwidth, 1);
    final PolynomialSplineFunction spline1 = loess.interpolate(h1[0], h1[1]);
    final PolynomialSplineFunction spline2 = loess.interpolate(h1[0], h2[1]);
    final PolynomialSplineFunction spline3 = loess.interpolate(h1[0], h3[1]);
    // Use a second interpolator in case the LOESS fails
    final LinearInterpolator lin = new LinearInterpolator();
    final PolynomialSplineFunction spline1b = lin.interpolate(h1[0], h1[1]);
    final PolynomialSplineFunction spline2b = lin.interpolate(h1[0], h2[1]);
    final PolynomialSplineFunction spline3b = lin.interpolate(h1[0], h3[1]);

    // Increase the number of points to show a smooth curve
    final double[] points = new double[bins * 5];
    limits = MathUtils.limits(h1[0]);
    final double interval = (limits[1] - limits[0]) / (points.length - 1);
    final double[] v = new double[points.length];
    final double[] v2 = new double[points.length];
    final double[] v3 = new double[points.length];
    for (int i = 0; i < points.length - 1; i++) {
      points[i] = limits[0] + i * interval;
      v[i] = getSplineValue(spline1, spline1b, points[i]);
      v2[i] = getSplineValue(spline2, spline2b, points[i]);
      v3[i] = getSplineValue(spline3, spline3b, points[i]);
      points[i] += halfBinWidth;
    }
    // Final point on the limit of the spline range
    final int ii = points.length - 1;
    v[ii] = getSplineValue(spline1, spline1b, limits[1]);
    v2[ii] = getSplineValue(spline2, spline2b, limits[1]);
    v3[ii] = getSplineValue(spline3, spline3b, limits[1]);
    points[ii] = limits[1] + halfBinWidth;

    // Calculate recall
    for (int i = 0; i < v.length; i++) {
      v2[i] = v2[i] / v[i];
      v3[i] = v3[i] / v[i];
    }

    final double halfSummaryDepth = settings.summaryDepth * 0.5;

    final String title2 = TITLE + " Depth Histogram (normalised)";
    final Plot plot2 = new Plot(title2, "Depth (nm)", "Recall");
    plot2.setLimits(limits[0] + halfBinWidth, limits[1] + halfBinWidth, 0,
        MathUtils.min(1, MathUtils.max(v2)));
    plot2.setColor(Color.black);
    plot2.addLabel(0, 0, "Blue = Fitted; Red = Filtered");
    plot2.setColor(Color.blue);
    plot2.addPoints(points, v2, Plot.LINE);
    plot2.setColor(Color.red);
    plot2.addPoints(points, v3, Plot.LINE);
    plot2.setColor(Color.magenta);
    if (-halfSummaryDepth - halfBinWidth >= limits[0]) {
      plot2.drawLine(-halfSummaryDepth, 0, -halfSummaryDepth,
          getSplineValue(spline3, spline3b, -halfSummaryDepth - halfBinWidth)
              / getSplineValue(spline1, spline1b, -halfSummaryDepth - halfBinWidth));
    }
    if (halfSummaryDepth - halfBinWidth <= limits[1]) {
      plot2.drawLine(halfSummaryDepth, 0, halfSummaryDepth,
          getSplineValue(spline3, spline3b, halfSummaryDepth - halfBinWidth)
              / getSplineValue(spline1, spline1b, halfSummaryDepth - halfBinWidth));
    }
    ImageJUtils.display(title2, plot2, wo);

    return allAssignments;
  }

  private static double getSplineValue(PolynomialSplineFunction spline,
      PolynomialSplineFunction spline2, double x) {
    double y = spline.value(x);
    if (Double.isNaN(y)) {
      y = spline2.value(x);
    }
    return y;
  }

  /**
   * Score analysis.
   *
   * @param allAssignments The assignments generated from running the filter (or null)
   * @param filter the filter
   * @return the assignments
   */
  @Nullable
  private ArrayList<FractionalAssignment[]>
      scoreAnalysis(ArrayList<FractionalAssignment[]> allAssignments, DirectFilter filter) {
    if (!settings.scoreAnalysis) {
      return null;
    }

    // Build a histogram of the fitted spots that were available to be scored
    final double[] signal = fitResultData.signalFactorStats.getValues();
    final double[] distance = fitResultData.distanceStats.getValues();

    double[] limits1;
    if (fitSignalFactor > 0 && settings.upperSignalFactor > 0) {
      final double range = fitSignalFactor * settings.upperSignalFactor / 100.0;
      limits1 = new double[] {-range, range};
    } else {
      limits1 = MathUtils.limits(signal);
      // Prevent the auto-range being too big
      final double bound = 3;
      if (limits1[0] < -bound) {
        limits1[0] = -bound;
      }
      if (limits1[1] > bound) {
        limits1[1] = bound;
      }
    }

    double[] limits2;
    if (spotFitResults.distanceInPixels > 0 && settings.upperMatchDistance > 0) {
      final double range = simulationParameters.pixelPitch * spotFitResults.distanceInPixels
          * settings.upperMatchDistance / 100.0;
      limits2 = new double[] {0, range};
    } else {
      limits2 = MathUtils.limits(distance);
    }

    final int bins = HistogramPlot.getBinsSqrtRule(signal.length);
    final double[][] h1 = HistogramPlot.calcHistogram(signal, limits1[0], limits1[1], bins);
    final double[][] h2 = HistogramPlot.calcHistogram(distance, limits2[0], limits2[1], bins);

    // Run the filter manually to get the results that pass.
    if (allAssignments == null) {
      allAssignments = getAssignments(filter);
    }

    double[] signal2 = new double[results.size()];
    double[] distance2 = new double[results.size()];
    int count = 0;
    double sumSignal = 0;
    double sumDistance = 0;
    for (final FractionalAssignment[] assignments : allAssignments) {
      if (assignments == null) {
        continue;
      }
      for (int i = 0; i < assignments.length; i++) {
        final CustomFractionalAssignment c = (CustomFractionalAssignment) assignments[i];
        sumDistance += distance2[count] = c.distToTarget;
        sumSignal += signal2[count] = c.getSignalFactor();
        count++;
      }
    }
    signal2 = Arrays.copyOf(signal2, count);
    distance2 = Arrays.copyOf(distance2, count);

    // Build a histogram using the same limits
    final double[][] h1b = HistogramPlot.calcHistogram(signal2, limits1[0], limits1[1], bins);
    final double[][] h2b = HistogramPlot.calcHistogram(distance2, limits2[0], limits2[1], bins);

    // Since the distance and signal factor are computed for all fits (single, multi, doublet)
    // there will be far more of them so we normalise and just plot the histogram profile.
    double s1 = 0;
    double s2 = 0;
    double s1b = 0;
    double s2b = 0;
    for (int i = 0; i < h1b[0].length; i++) {
      s1 += h1[1][i];
      s2 += h2[1][i];
      s1b += h1b[1][i];
      s2b += h2b[1][i];
    }
    for (int i = 0; i < h1b[0].length; i++) {
      h1[1][i] /= s1;
      h2[1][i] /= s2;
      h1b[1][i] /= s1b;
      h2b[1][i] /= s2b;
    }

    // Draw distance histogram first
    final String title2 = TITLE + " Distance Histogram";
    final Plot plot2 = new Plot(title2, "Distance (nm)", "Frequency");
    plot2.setLimits(limits2[0], limits2[1], 0, MathUtils.maxDefault(MathUtils.max(h2[1]), h2b[1]));
    plot2.setColor(Color.black);
    plot2.addLabel(0, 0,
        String.format("Blue = Fitted (%s); Red = Filtered (%s)",
            MathUtils.rounded(fitResultData.distanceStats.getMean()),
            MathUtils.rounded(sumDistance / count)));
    plot2.setColor(Color.blue);
    plot2.addPoints(h2[0], h2[1], Plot.BAR);
    plot2.setColor(Color.red);
    plot2.addPoints(h2b[0], h2b[1], Plot.BAR);
    ImageJUtils.display(title2, plot2, wo);

    // Draw signal factor histogram
    final String title1 = TITLE + " Signal Factor Histogram";
    final Plot plot1 = new Plot(title1, "Signal Factor", "Frequency");
    plot1.setLimits(limits1[0], limits1[1], 0, MathUtils.maxDefault(MathUtils.max(h1[1]), h1b[1]));
    plot1.setColor(Color.black);
    plot1.addLabel(0, 0,
        String.format("Blue = Fitted (%s); Red = Filtered (%s)",
            MathUtils.rounded(fitResultData.signalFactorStats.getMean()),
            MathUtils.rounded(sumSignal / count)));
    plot1.setColor(Color.blue);
    plot1.addPoints(h1[0], h1[1], Plot.BAR);
    plot1.setColor(Color.red);
    plot1.addPoints(h1b[0], h1b[1], Plot.BAR);
    ImageJUtils.display(title1, plot1, wo);

    return allAssignments;
  }

  private void componentAnalysis(ComplexFilterScore bestFilterScore) {
    if (settings.componentAnalysis == 0) {
      return;
    }
    final Consumer<String> output = createComponentAnalysisWindow();

    final DirectFilter localBestFilter = bestFilterScore.getFilter();
    final String[] names = getNames(localBestFilter);

    // Skip disabled parameters
    int paramCount = localBestFilter.getNumberOfParameters();
    final boolean[] enable = new boolean[paramCount];
    final int[] map = new int[paramCount];
    for (int n = 0, i = 0; n < map.length; n++) {
      if (localBestFilter.getParameterValue(n) == localBestFilter.getDisabledParameterValue(n)) {
        enable[n] = true; // Mark to ignore
        paramCount--;
      } else {
        map[i++] = n;
      }
    }

    // Score the best filter just so we have the unique Ids of the results
    scoreComponents(localBestFilter, -1, paramCount, null, null, null, 0);
    final int[] uniqueIds1 = uniqueIds;
    final int uniqueIdCount1 = uniqueIdCount;

    // Limit to 12 params == 4095 combinations (this is the max for two multi filters combined)
    if (settings.componentAnalysis >= 3 && paramCount <= 12) {
      // Enumeration of combinations
      final long count = countComponentCombinations(paramCount);

      // Enumerate all combinations
      final ComplexFilterScore[] localScores = new ComplexFilterScore[(int) count];
      int index = 0;

      for (int k = 1; k <= paramCount; k++) {
        final Iterator<int[]> it = CombinatoricsUtils.combinationsIterator(paramCount, k);
        while (it.hasNext()) {
          final int[] combinations = it.next();
          final boolean[] enable2 = enable.clone();
          for (int i = 0; i < k; i++) {
            combinations[i] = map[combinations[i]];
            enable2[combinations[i]] = true;
          }
          final DirectFilter f = (DirectFilter) localBestFilter.create(enable2);
          localScores[index++] =
              scoreComponents(f, 0, k, combinations, enable2, uniqueIds1, uniqueIdCount1);
        }
      }

      // Report
      Arrays.sort(localScores, FilterScoreCompararor.INSTANCE);

      int lastSize = 0;
      for (int i = 0; i < localScores.length; i++) {
        if (settings.componentAnalysis == 3) {
          if (lastSize == localScores[i].size) {
            // Only add the best result for each size
            continue;
          }
          lastSize = localScores[i].size;
        }

        // Find the last component that was added
        if (localScores[i].size == 1) {
          localScores[i].index = localScores[i].combinations[0];
        } else {
          // For each size k, find the best result with k-1 components and set the add index
          // appropriately
          int add = -1;
          int target = -1;
          for (int l = 0; l < enable.length; l++) {
            if (localScores[i].enable[l]) {
              target++;
            }
          }
          final int size1 = localScores[i].size - 1;
          for (int ii = 0; ii < i; ii++) {
            if (localScores[ii].size < size1) {
              continue;
            }
            if (localScores[ii].size > size1) {
              break; // Broken
            }
            // Count matches. It must be 1 less than the current result
            int matches = 0;
            for (int l = 0; l < enable.length; l++) {
              if (localScores[i].enable[l] && localScores[ii].enable[l]) {
                matches++;
              }
            }
            if (matches == target) {
              // Find the additional parameter added
              for (int l = 0; l < enable.length; l++) {
                if (localScores[i].enable[l]) {
                  if (!localScores[ii].enable[l]) {
                    add = l;
                    break;
                  }
                }
              }
              break;
            }
          }
          localScores[i].index = add;
        }

        addToComponentAnalysisWindow(output, localScores[i], bestFilterScore, names);
      }

      return;
    }

    // Preserve the option to output the best or all results if we fell through from above
    final int myComponentAnalysis =
        (settings.componentAnalysis >= 3) ? settings.componentAnalysis - 2
            : settings.componentAnalysis;

    // Progressively add components until all are the same as the input bestFilter
    int enabled = 0;
    int[] previousCombinations = new int[0];
    for (int ii = 0; ii < paramCount; ii++) {
      // Create a set of filters by enabling each component that is not currently enabled.
      final ComplexFilterScore[] localScores = new ComplexFilterScore[paramCount - enabled];
      final int k = enabled + 1;
      for (int i = 0, j = 0; i < paramCount; i++) {
        final int n = map[i];
        if (enable[n]) {
          continue;
        }
        enable[n] = true;
        final DirectFilter f = (DirectFilter) localBestFilter.create(enable);
        enable[n] = false;
        final int[] combinations = new int[k];
        System.arraycopy(previousCombinations, 0, combinations, 0, previousCombinations.length);
        combinations[k - 1] = n;
        Arrays.sort(combinations);
        localScores[j++] = scoreComponents(f, n, k, combinations, null, uniqueIds1, uniqueIdCount1);
      }

      // Rank them
      Arrays.sort(localScores);
      for (int i = 0; i < localScores.length; i++) {
        addToComponentAnalysisWindow(output, localScores[i], bestFilterScore, names);
        if (myComponentAnalysis == 1) {
          // Only add the best result
          break;
        }
      }

      // Flag the best component added
      enable[localScores[0].index] = true;
      enabled++;
      previousCombinations = localScores[0].combinations;
    }
  }

  private ComplexFilterScore scoreComponents(DirectFilter filter, int index, int size,
      int[] combinations, boolean[] enable, int[] uniqueIds1, int uniqueIdCount1) {
    setupFractionScoreStore();

    // Score them
    long time = System.nanoTime();
    final FilterScoreResult r = scoreFilter(filter);
    time = System.nanoTime() - time;

    endFractionScoreStore();

    ClassificationResult r2 = null;
    if (uniqueIds1 != null) {
      // Build an overlap between the results created by this filter and the best filter.
      // Note that the overlap may be very low given the number of different ways we can generate
      // fit results (i.e. it is not as simple as just single-fitting on each candidate)
      // To do this we assign a unique ID to each possible result.
      // The two sets can be compared to produce a Precision, Recall, Jaccard, etc.

      // PreprocessedPeakResult will need an additional mutable Id field.
      // Sort by Id then iterate through both arrays concurrently counting the matching Ids.
      final int[] uniqueIds2 = uniqueIds;
      final int uniqueIdCount2 = uniqueIdCount;

      int tp = 0;
      int fp = 0;
      int fn = 0;

      // Compare Ids (must be sorted in ascending order)
      int i1 = 0;
      int i2 = 0;
      while (i1 < uniqueIdCount1 && i2 < uniqueIdCount2) {
        final int result = uniqueIds1[i1] - uniqueIds2[i2];
        if (result > 0) {
          i2++;
          fp++;
        } else if (result < 0) {
          i1++;
          fn++;
        }
        i1++;
        i2++;
        tp++;
      }
      // Count the remaining ids
      fn += (uniqueIdCount1 - i1);
      fp += (uniqueIdCount2 - i2);

      r2 = new ClassificationResult(tp, fp, 0, fn);
    }

    // These scores are used when the same filter type so set allSameType to true
    return new ComplexFilterScore(r, true, index, time, r2, size, combinations, enable);
  }

  private static String[] getNames(DirectFilter bestFilter) {

    final int nParams = bestFilter.getNumberOfParameters();
    final String[] names = new String[nParams];

    for (int n = 0; n < nParams; n++) {
      final String name = bestFilter.getParameterName(n);
      String uniqueName = name;
      int count = 1;
      // Avoid duplicates
      while (contains(names, n, uniqueName)) {
        uniqueName = name + (++count);
      }
      names[n] = uniqueName;
    }
    return names;
  }

  private static boolean contains(String[] names, int n, String uniqueName) {
    for (int i = 0; i < n; i++) {
      if (names[i].equals(uniqueName)) {
        return true;
      }
    }
    return false;
  }

  /**
   * Count component combinations.
   *
   * @param n the number of parameters
   * @return the combinations
   */
  private static long countComponentCombinations(int n) {
    return (long) FastMath.pow(2, n) - 1;

    // This returns the same as (2^n)-1:
    // long total = 0
    // for k : 1 <= n
    // total += CombinatoricsUtils.binomialCoefficient(n, k)
    // return total
  }

  private void setupFractionScoreStore() {
    uniqueIds = new int[fitResultData.maxUniqueId];
    uniqueIdCount = 0;
    // Store the Id of each result
    scoreStore = uniqueId -> uniqueIds[uniqueIdCount++] = uniqueId;
  }

  private void endFractionScoreStore() {
    scoreStore = null;
    if (uniqueIdCount != 0) {
      Arrays.sort(uniqueIds, 0, uniqueIdCount);
    }
  }

  /**
   * Sets the strength on all the filters.
   *
   * @param filterSet the filter set
   * @return the filter set
   */
  private FilterSet setStrength(FilterSet filterSet) {
    if (strengthLower != null) {
      for (final Filter f : filterSet.getFilters()) {
        final DirectFilter df = (DirectFilter) f;
        df.setStrength(df.computeStrength(strengthLower, strengthUpper));
      }
    }
    return filterSet;
  }

  /**
   * Sets the strength on all the filters if not computed.
   *
   * @param filterSet the filter set
   * @return the filter set
   */
  private FilterSet setUncomputedStrength(FilterSet filterSet) {
    if (strengthLower != null) {
      for (final Filter f : filterSet.getFilters()) {
        final DirectFilter df = (DirectFilter) f;
        if (Float.isNaN(df.getStrength())) {
          df.setStrength(df.computeStrength(strengthLower, strengthUpper));
        }
      }
    }
    return filterSet;
  }

  @Nullable
  private FilterScoreResult[] scoreFilters(FilterSet filterSet, boolean createTextResult) {
    if (filterSet.size() == 0) {
      return null;
    }

    initialiseScoring(filterSet);

    FilterScoreResult[] scoreResults = new FilterScoreResult[filterSet.size()];

    if (scoreResults.length == 1) {
      // No need to multi-thread this
      scoreResults[0] = scoreFilter((DirectFilter) filterSet.getFilters().get(0),
          defaultMinimalFilter, createTextResult, coordinateStore);
    } else {
      // Multi-thread score all the result
      final int nThreads = getThreads(scoreResults.length);
      final BlockingQueue<ScoreJob> jobs = new ArrayBlockingQueue<>(nThreads * 2);
      final List<Thread> threads = new LinkedList<>();
      final Ticker ticker =
          ImageJUtils.createTicker(scoreResults.length, nThreads, "Scoring Filters");
      for (int i = 0; i < nThreads; i++) {
        final ScoreWorker worker = new ScoreWorker(jobs, scoreResults, createTextResult,
            (coordinateStore == null) ? null : coordinateStore.newInstance(), ticker);
        final Thread t = new Thread(worker);
        threads.add(t);
        t.start();
      }

      int index = 0;
      for (final Filter filter : filterSet.getFilters()) {
        if (IJ.escapePressed()) {
          break;
        }
        put(jobs, new ScoreJob((DirectFilter) filter, index++));
      }
      // Finish all the worker threads by passing in a null job
      for (int i = 0; i < threads.size(); i++) {
        put(jobs, new ScoreJob(null, -1));
      }

      // Wait for all to finish
      for (int i = 0; i < threads.size(); i++) {
        try {
          threads.get(i).join();
        } catch (final InterruptedException ex) {
          Logger.getLogger(BenchmarkFilterAnalysis.class.getName()).log(Level.WARNING,
              "Interrupted!", ex);
          Thread.currentThread().interrupt();
          throw new ConcurrentRuntimeException("Unexpected interruption", ex);
        }
      }
      threads.clear();
      ImageJUtils.finished();

      // In case the threads were interrupted
      if (ImageJUtils.isInterrupted()) {
        scoreResults = null;
      }
    }

    finishScoring();

    return scoreResults;
  }

  /**
   * Score filters.
   *
   * @param points the points (must be sorted by duplicate distance)
   * @param createTextResult set to true to create the text result
   * @return the score results
   */
  @Nullable
  private ParameterScoreResult[] scoreFilters(double[][] points, boolean createTextResult) {
    if (points == null || points.length == 0) {
      return null;
    }

    gaResultsListToScore = gaResultsList;
    gaSubset = false;

    ParameterScoreResult[] scoreResults = new ParameterScoreResult[points.length];

    if (scoreResults.length == 1) {
      // No need to multi-thread this
      final int failCount = (int) Math.round(points[0][0]);
      final double residualsThreshold = points[0][1];
      final double duplicateDistance = points[0][2];
      scoreResults[0] =
          scoreFilter(searchScoreFilter, defaultMinimalFilter, failCount, residualsThreshold,
              duplicateDistance, createCoordinateStore(duplicateDistance), createTextResult);
    } else {
      // Multi-thread score all the result
      final int nThreads = getThreads(scoreResults.length);
      final BlockingQueue<ParameterScoreJob> jobs = new ArrayBlockingQueue<>(nThreads * 2);
      final List<Thread> threads = new LinkedList<>();
      final Ticker ticker =
          ImageJUtils.createTicker(scoreResults.length, nThreads, "Scoring Filters");
      for (int i = 0; i < nThreads; i++) {
        final ParameterScoreWorker worker =
            new ParameterScoreWorker(jobs, scoreResults, createTextResult, ticker);
        final Thread t = new Thread(worker);
        threads.add(t);
        t.start();
      }

      ticker.start();
      for (int i = 0; i < points.length; i++) {
        if (IJ.escapePressed()) {
          break;
        }
        put(jobs, new ParameterScoreJob(points[i], i));
      }
      // Finish all the worker threads by passing in a null job
      for (int i = 0; i < threads.size(); i++) {
        put(jobs, new ParameterScoreJob(null, -1));
      }

      // Wait for all to finish
      for (int i = 0; i < threads.size(); i++) {
        try {
          threads.get(i).join();
        } catch (final InterruptedException ex) {
          Logger.getLogger(BenchmarkFilterAnalysis.class.getName()).log(Level.WARNING,
              "Interrupted!", ex);
          Thread.currentThread().interrupt();
          throw new ConcurrentRuntimeException("Unexpected interruption", ex);
        }
      }
      threads.clear();
      ImageJUtils.finished();

      // In case the threads were interrupted
      if (ImageJUtils.isInterrupted()) {
        scoreResults = null;
      }
    }

    finishScoring();

    return scoreResults;
  }

  private static int getThreads(int length) {
    return Math.max(1, Math.min(Prefs.getThreads(), length));
  }

  private static <T> void put(BlockingQueue<T> jobs, T job) {
    try {
      jobs.put(job);
    } catch (final InterruptedException ex) {
      Logger.getLogger(BenchmarkFilterAnalysis.class.getName()).log(Level.WARNING, "Interrupted!",
          ex);
      Thread.currentThread().interrupt();
      throw new ConcurrentRuntimeException("Unexpected interruption", ex);
    }
  }

  /**
   * Initialise the results list used for scoring the filters. This is shared with the genetic
   * algorithm.
   *
   * @param filterSet the filter set
   */
  private void initialiseScoring(FilterSet filterSet) {
    // Initialise with the candidate true and false negative scores
    gaResultsListToScore = gaResultsList;
    gaSubset = false;

    if (filterSet.size() < 2) {
      return;
    }

    if (filterSet.allSameType() && debug) {
      // Would there be a speed increase if the results list were partitioned into multiple sets
      // by filtering with not just the weakest but a step set of weaker and weaker filters.
      // This could be done using, e.g. precision, to partition the filters into a range.

      // // Plot the cumulative histogram of precision for the filters in the set.
      // try
      // {
      // StoredDataStatistics stats = new StoredDataStatistics(filterSet.size());
      // for (Filter f : filterSet.getFilters())
      // {
      // IMultiFilter f2 = (IMultiFilter) f;
      // stats.add(f2.getPrecision());
      // }
      // double[][] h1 = Maths.cumulativeHistogram(stats.getValues(), false);
      // String title = TITLE + " Cumul Precision";
      // Plot plot = new Plot(title, "Precision", "Frequency");
      // // Find limits
      // double[] xlimit = Maths.limits(h1[0]);
      // plot.setLimits(xlimit[0] - 1, xlimit[1] + 1, 0, Maths.max(h1[1]) * 1.05);
      // plot.addPoints(h1[0], h1[1], Plot.BAR);
      // Utils.display(title, plot);
      // }
      // catch (ClassCastException e)
      // {
      //
      // }

      // Debug the limits of all parameters
      final double[] lower = filterSet.getFilters().get(0).getParameters().clone();
      final double[] upper = lower.clone();
      for (final Filter f : filterSet.getFilters()) {
        final double[] point = f.getParameters();
        for (int j = 0; j < lower.length; j++) {
          if (lower[j] > point[j]) {
            lower[j] = point[j];
          }
          if (upper[j] < point[j]) {
            upper[j] = point[j];
          }
        }
      }
      final StringBuilder sb = new StringBuilder("Scoring (");
      sb.append(filterSet.size()).append("):");
      for (int j = 0; j < lower.length; j++) {
        sb.append(' ').append(MathUtils.rounded(lower[j])).append('-')
            .append(MathUtils.rounded(upper[j]));
      }
      ImageJUtils.log(sb.toString());
    }

    final Filter weakest = filterSet.createWeakestFilter();
    if (weakest != null) {
      gaSubset = true;

      gaResultsListToScore = createMpf((DirectFilter) weakest, defaultMinimalFilter)
          .filterSubset(gaResultsList, createFailCounter(settings.failCount), true);

      // MultiPathFilter.resetValidationFlag(ga_resultsListToScore);

      // ga_resultsListToScore = ga_resultsList;

      // System.out.printf("Weakest %d => %d : %s\n", count(ga_resultsList),
      // count(ga_resultsListToScore),
      // weakest.getName());
    }
  }

  /**
   * Finish scoring and reset the subset.
   */
  private void finishScoring() {
    if (gaSubset) {
      // Reset the validation flag
      MultiPathFilter.resetValidationFlag(gaResultsListToScore);
    }
  }

  @Override
  public void initialise(List<? extends Chromosome<FilterScore>> individuals) {
    gaIteration++;
    gaScoreIndex = 0;
    gaScoreResults =
        scoreFilters(setStrength(new FilterSet(populationToFilters(individuals))), false);
  }

  private static ArrayList<Filter>
      populationToFilters(List<? extends Chromosome<FilterScore>> individuals) {
    final ArrayList<Filter> filters = new ArrayList<>(individuals.size());
    for (final Chromosome<FilterScore> c : individuals) {
      filters.add((DirectFilter) c);
    }
    return filters;
  }

  private List<Filter> searchSpaceToFilters(double[][] searchSpace) {
    return searchSpaceToFilters(null, searchSpace);
  }

  private List<Filter> searchSpaceToFilters(DirectFilter filter, double[][] searchSpace) {
    final int size = ((searchSpace == null) ? 0 : searchSpace.length) + ((filter == null) ? 0 : 1);
    final ArrayList<Filter> filters = new ArrayList<>(size);
    if (filter != null) {
      filters.add(filter);
    }
    if (searchSpace != null) {
      for (int i = 0; i < searchSpace.length; i++) {
        filters.add(searchScoreFilter.create(searchSpace[i]));
      }
    }
    return filters;
  }

  @Override
  public FilterScore fitness(Chromosome<FilterScore> chromosome) {
    // In case the user aborted with Escape
    if (gaScoreResults == null) {
      return null;
    }

    // Assume that fitness will be called in the order of the individuals passed to the initialise
    // function.
    final FilterScoreResult scoreResult = gaScoreResults[gaScoreIndex++];

    // Set this to null and it will be removed at the next population selection
    if (scoreResult.score == 0) {
      return null;
    }

    return new SimpleFilterScore(scoreResult, true, scoreResult.criteria >= minCriteria);
  }

  @Override
  public void shutdown() {
    // Report the score for the best filter
    final List<? extends Chromosome<FilterScore>> individuals = gaPopulation.getIndividuals();

    FilterScore max = null;
    for (final Chromosome<FilterScore> c : individuals) {
      final FilterScore f = c.getFitness();
      if (f != null && f.compareTo(max) < 0) {
        max = f;
      }
    }

    if (max == null) {
      return;
    }

    final DirectFilter filter = (DirectFilter) ((SimpleFilterScore) max).filter;

    // This filter may not have been part of the scored subset so use the entire results set for
    // reporting
    final FractionClassificationResult r =
        scoreFilter(filter, defaultMinimalFilter, gaResultsList, coordinateStore);

    final StringBuilder text = createResult(filter, r);
    add(text, gaIteration);
    gaWindow.accept(text.toString());
  }

  @Nullable
  @Override
  public SearchResult<FilterScore> findOptimum(double[][] points) {
    gaIteration++;
    SimpleFilterScore max = filterScoreOptimum;

    final FilterScoreResult[] scoreResults =
        scoreFilters(setStrength(new FilterSet(searchSpaceToFilters(points))), false);

    if (scoreResults == null) {
      return null;
    }

    for (int index = 0; index < scoreResults.length; index++) {
      final FilterScoreResult scoreResult = scoreResults[index];
      final SimpleFilterScore result =
          new SimpleFilterScore(scoreResult, true, scoreResult.criteria >= minCriteria);
      if (result.compareTo(max) < 0) {
        max = result;
      }
    }

    filterScoreOptimum = max;

    // Add the best filter to the table
    // This filter may not have been part of the scored subset so use the entire results set for
    // reporting
    final DirectFilter filter = max.result.filter;
    final FractionClassificationResult r =
        scoreFilter(filter, defaultMinimalFilter, gaResultsList, coordinateStore);

    final StringBuilder text = createResult(filter, r);
    add(text, gaIteration);
    gaWindow.accept(text.toString());

    return new SearchResult<>(filter.getParameters(), max);
  }

  @Nullable
  @Override
  public SearchResult<FilterScore>[] score(double[][] points) {
    gaIteration++;
    SimpleFilterScore max = filterScoreOptimum;

    final FilterScoreResult[] scoreResults =
        scoreFilters(setStrength(new FilterSet(searchSpaceToFilters(points))), false);

    if (scoreResults == null) {
      return null;
    }

    @SuppressWarnings("unchecked")
    final SearchResult<FilterScore>[] scores = new SearchResult[scoreResults.length];
    for (int index = 0; index < scoreResults.length; index++) {
      final FilterScoreResult scoreResult = scoreResults[index];
      final SimpleFilterScore result =
          new SimpleFilterScore(scoreResult, true, scoreResult.criteria >= minCriteria);
      if (result.compareTo(max) < 0) {
        max = result;
      }
      scores[index] = new SearchResult<>(result.result.filter.getParameters(), result);
    }

    filterScoreOptimum = max;

    // Add the best filter to the table
    // This filter may not have been part of the scored subset so use the entire results set for
    // reporting
    final DirectFilter filter = max.result.filter;
    final FractionClassificationResult r =
        scoreFilter(filter, defaultMinimalFilter, gaResultsList, coordinateStore);

    final StringBuilder text = createResult(filter, r);
    add(text, gaIteration);
    gaWindow.accept(text.toString());

    return scores;
  }

  @Override
  @SuppressWarnings("unchecked")
  public SearchResult<FilterScore>[] cut(SearchResult<FilterScore>[] scores, int size) {
    // Do a full sort and truncation
    // return ScoreFunctionHelper.cut(scores, size);

    // Split the list into those that pass the criteria and those that do not
    final SearchResult<FilterScore>[] passList = new SearchResult[scores.length];
    final SearchResult<FilterScore>[] failList = new SearchResult[scores.length];
    int pass = 0;
    int fail = 0;
    for (int i = 0; i < scores.length; i++) {
      // Ignore this
      if (scores[i].getScore().score == 0) {
        continue;
      }
      if (scores[i].getScore().criteriaPassed) {
        passList[pass++] = scores[i];
      } else {
        failList[fail++] = scores[i];
      }
    }

    // If those that pass are bigger than size then sort that list and return.
    if (pass >= size) {
      return ScoreFunctionHelper.cut(Arrays.copyOf(passList, pass), size);
    }

    // Sort the fail list
    Arrays.sort(failList, 0, fail);

    // Find the top score and put it first.
    if (pass != 0) {
      int best = 0;
      for (int i = 1; i < pass; i++) {
        if (passList[i].compareTo(passList[best]) < 0) {
          best = i;
        }
      }
      final SearchResult<FilterScore> tmp = passList[best];
      passList[best] = passList[0];
      passList[0] = tmp;
    }

    // Combine the lists. Account for removing zero scores (i.e. pass+fail<=size)
    size = Math.min(size, pass + fail);
    System.arraycopy(failList, 0, passList, pass, size - pass);

    return Arrays.copyOf(passList, size);
  }

  @Override
  public void progress(double fraction) {
    if (fraction == 1) {
      // Reset
      limit = 0;
      IJ.showProgress(fraction);
      return;
    }

    // Show only 2% changes
    if (fraction < limit) {
      return;
    }

    limit = fraction + 0.02;
    IJ.showProgress(fraction);
  }

  @Override
  public void progress(long position, long total) {
    progress((double) position / total);
  }

  @Override
  public void incrementProgress(double fraction) {
    // Ignore
  }

  @Override
  public void log(String format, Object... args) {
    // Ignore
  }

  @Override
  public void status(String format, Object... args) {
    IJ.showStatus(gaStatusPrefix + String.format(format, args));
  }

  @Override
  public boolean isEnded() {
    // Ignore
    return false;
  }

  @Override
  public boolean isProgress() {
    return true;
  }

  @Override
  public boolean isLog() {
    return false;
  }

  @Override
  public boolean isStatus() {
    return true;
  }

  /**
   * Updates the given configuration using the latest settings used in benchmarking spot filters,
   * fitting and filtering.
   *
   * @param config the configuration
   * @param useLatest Use the latest best filter. Otherwise use the highest scoring.
   * @return true, if successful
   */
  public static boolean updateAllConfiguration(FitEngineConfiguration config, boolean useLatest) {
    // Do this first as it sets the initial SD
    if (!updateAllConfiguration(config)) {
      return false;
    }
    return updateConfiguration(config, useLatest);
  }

  /**
   * Updates the given configuration using the latest settings used in benchmarking spot filters and
   * fitting.
   *
   * @param config the configuration
   * @return true, if successful
   */
  private static boolean updateAllConfiguration(FitEngineConfiguration config) {
    // Do this first as it sets the initial SD
    if (!BenchmarkSpotFit.updateConfiguration(config)) {
      return false;
    }
    return BenchmarkSpotFilter.updateConfiguration(config);
  }

  /**
   * Updates the given configuration using the latest settings used in benchmarking filtering. The
   * residuals threshold will be copied only if the input FitConfiguration has isComputeResiduals()
   * set to true.
   *
   * <p>This calls {@link FitConfiguration#setDirectFilter(DirectFilter)} and sets the precision
   * method using the method in the direct filter.
   *
   * @param config the configuration
   * @param useLatest Use the latest best filter. Otherwise use the highest scoring.
   * @return true, if successful
   */
  public static boolean updateConfiguration(FitEngineConfiguration config, boolean useLatest) {
    final BenchmarkFilterAnalysisResult lastResult = BenchmarkFilterAnalysisResult.lastResult.get();
    if (lastResult.scores.isEmpty()) {
      return false;
    }

    FilterResult best;
    if (useLatest) {
      best = lastResult.scores.get(lastResult.scores.size() - 1);
    } else {
      best = getBestResult(lastResult.scores);
    }

    // New smart filter support
    final FitConfiguration fitConfig = config.getFitConfiguration();
    fitConfig.setDirectFilter(best.getFilter());

    // Set the precision method using the direct filter
    fitConfig.setPrecisionMethod(fitConfig.getFilterPrecisionMethod());

    if (fitConfig.isComputeResiduals()) {
      config.setResidualsThreshold(best.residualsThreshold);
      fitConfig.setComputeResiduals(true);
    } else {
      config.setResidualsThreshold(1);
      fitConfig.setComputeResiduals(false);
    }

    // Note:
    // We leave the simple filter settings alone. These may be enabled as well, e.g. by the
    // BenchmarkSpotFit plugin

    // We could set the fail count range dynamically using a window around the best filter

    config.setFailuresLimit(best.failCount);

    config.setDuplicateDistance(best.duplicateDistance);
    config.setDuplicateDistanceAbsolute(best.duplicateDistanceAbsolute);

    return true;
  }

  private static FilterResult getBestResult(List<FilterResult> scores) {
    if (!scores.isEmpty()) {
      // Copy so we can sort the list
      final LocalList<FilterResult> scores2 = new LocalList<>(scores);

      // Ideally we pick the top score. However if many filters are close we should
      // pick the filter with the most favourable properties for fitting.

      // Sort by score.
      scores2.sort((o1, o2) -> Double.compare(o2.score, o1.score));

      // For any scores close to the top score sort again using
      // other filter parameters: fail count, residuals threshold, duplicate distance
      // This should favour filters with lower fail count, higher residuals threshold and
      // lower duplicate distance.
      final double threshold = scores2.unsafeGet(0).score * 0.9999;
      final int index = scores2.findIndex(result -> result.score < threshold);
      if (index > 0) {
        scores2.subList(0, index).sort((o1, o2) -> {
          // Lower fail count will be faster during fitting
          int result = Integer.compare(o1.failCount, o2.failCount);
          if (result != 0) {
            return result;
          }
          // Higher residuals threshold is more conservative
          result = Double.compare(o2.residualsThreshold, o1.residualsThreshold);
          if (result != 0) {
            return result;
          }
          // Lower duplicate distance is more conservative
          return Double.compare(o1.duplicateDistance, o2.duplicateDistance);
        });
      }
      return scores2.unsafeGet(0);
    }
    return null;
  }

  /**
   * Checks if any of the settings require an overlay.
   *
   * @return true if an overlay is required
   */
  private boolean isShowOverlay() {
    return (settings.showTP || settings.showFP || settings.showFN);
  }

  /**
   * Show overlay.
   *
   * <ul>
   *
   * <li>Green = TP
   *
   * <li>Red = FP
   *
   * <li>Magenta = FP (Ignored from analysis)
   *
   * <li>Yellow = FN
   *
   * <li>Orange = FN (Outside border)
   *
   * </ul>
   *
   * @param allAssignments The assignments generated from running the filter (or null)
   * @param filter the filter
   * @return The results from running the filter (or null)
   */
  @Nullable
  @SuppressWarnings("null")
  private PreprocessedPeakResult[] showOverlay(ArrayList<FractionalAssignment[]> allAssignments,
      DirectFilter filter) {
    final ImagePlus imp = CreateData.getImage();
    if (imp == null) {
      return null;
    }

    // Run the filter manually to get the results that pass.
    if (allAssignments == null) {
      allAssignments = getAssignments(filter);
    }

    final Overlay o = new Overlay();

    // Do TP
    final TIntHashSet actual = new TIntHashSet();
    final TIntHashSet predicted = new TIntHashSet();
    for (final FractionalAssignment[] assignments : allAssignments) {
      if (assignments == null || assignments.length == 0) {
        continue;
      }
      float[] tx = null;
      float[] ty = null;
      int count = 0;
      if (settings.showTP) {
        tx = new float[assignments.length];
        ty = new float[assignments.length];
      }
      int frame = 0;
      for (int i = 0; i < assignments.length; i++) {
        final CustomFractionalAssignment c = (CustomFractionalAssignment) assignments[i];
        final UniqueIdPeakResult peak = (UniqueIdPeakResult) c.peak;
        final BasePreprocessedPeakResult spot = (BasePreprocessedPeakResult) c.peakResult;
        actual.add(peak.uniqueId);
        predicted.add(spot.getUniqueId());
        frame = spot.getFrame();
        if (settings.showTP) {
          tx[count] = spot.getX();
          ty[count++] = spot.getY();
        }
      }
      if (settings.showTP) {
        SpotFinderPreview.addRoi(frame, o, tx, ty, count, Color.green);
      }
    }

    float[] x = new float[10];
    float[] y = new float[x.length];
    float[] x2 = new float[10];
    float[] y2 = new float[x2.length];

    // Do FP (all remaining results that are not a TP)
    PreprocessedPeakResult[] filterResults = null;
    if (settings.showFP) {
      final MultiPathFilter multiPathFilter = createMpf(filter, defaultMinimalFilter);

      filterResults = filterResults(multiPathFilter);

      int frame = 0;
      int c1 = 0;
      int c2 = 0;
      for (int i = 0; i < filterResults.length; i++) {
        if (frame != filterResults[i].getFrame()) {
          if (c1 != 0) {
            SpotFinderPreview.addRoi(frame, o, x, y, c1, Color.red);
          }
          if (c2 != 0) {
            SpotFinderPreview.addRoi(frame, o, x2, y2, c2, Color.magenta);
          }
          c1 = c2 = 0;
        }
        frame = filterResults[i].getFrame();
        if (predicted.contains(filterResults[i].getUniqueId())) {
          continue;
        }
        if (filterResults[i].ignore()) {
          if (x2.length == c2) {
            x2 = Arrays.copyOf(x2, c2 * 2);
            y2 = Arrays.copyOf(y2, c2 * 2);
          }
          x2[c2] = filterResults[i].getX();
          y2[c2++] = filterResults[i].getY();
        } else {
          if (x.length == c1) {
            x = Arrays.copyOf(x, c1 * 2);
            y = Arrays.copyOf(y, c1 * 2);
          }
          x[c1] = filterResults[i].getX();
          y[c1++] = filterResults[i].getY();
        }
      }
      if (c1 != 0) {
        SpotFinderPreview.addRoi(frame, o, x, y, c1, Color.red);
      }
      if (c2 != 0) {
        SpotFinderPreview.addRoi(frame, o, x2, y2, c2, Color.magenta);
      }
    }

    // Do FN (all remaining peaks that have not been matched)
    if (settings.showFN) {
      final boolean checkBorder =
          (filterResult.analysisBorder != null && filterResult.analysisBorder.x != 0);
      final float border;
      final float xlimit;
      final float ylimit;
      if (checkBorder) {
        final Rectangle lastAnalysisBorder = filterResult.analysisBorder;
        border = lastAnalysisBorder.x;
        xlimit = lastAnalysisBorder.x + lastAnalysisBorder.width;
        ylimit = lastAnalysisBorder.y + lastAnalysisBorder.height;
      } else {
        border = xlimit = ylimit = 0;
      }

      // Add the results to the lists
      actualCoordinates.forEachEntry(new CustomTIntObjectProcedure(x, y, x2, y2) {
        @Override
        public boolean execute(int frame, UniqueIdPeakResult[] results) {
          int c1 = 0;
          int c2 = 0;
          if (x.length <= results.length) {
            x = new float[results.length];
            y = new float[results.length];
          }
          if (x2.length <= results.length) {
            x2 = new float[results.length];
            y2 = new float[results.length];
          }

          for (int i = 0; i < results.length; i++) {
            // Ignore those that were matched by TP
            if (actual.contains(results[i].uniqueId)) {
              continue;
            }

            if (checkBorder && outsideBorder(results[i], border, xlimit, ylimit)) {
              x2[c2] = results[i].getXPosition();
              y2[c2++] = results[i].getYPosition();
            } else {
              x[c1] = results[i].getXPosition();
              y[c1++] = results[i].getYPosition();
            }
          }
          if (c1 != 0) {
            SpotFinderPreview.addRoi(frame, o, x, y, c1, Color.yellow);
          }
          if (c2 != 0) {
            SpotFinderPreview.addRoi(frame, o, x2, y2, c2, Color.orange);
          }
          return true;
        }
      });
    }

    imp.setOverlay(o);

    return filterResults;
  }

  /**
   * Check if the result is outside border.
   *
   * @param result the result
   * @param border the border
   * @param xlimit the xlimit
   * @param ylimit the ylimit
   * @return true, if outside the border (should be ignored)
   */
  private static boolean outsideBorder(PeakResult result, final float border, final float xlimit,
      final float ylimit) {
    // Respect the border used in BenchmarkSpotFilter?
    // Just implement a hard border
    return (result.getXPosition() < border || result.getXPosition() > xlimit
        || result.getYPosition() < border || result.getYPosition() > ylimit);
  }

  /**
   * Check if the result is outside border.
   *
   * @param result the result
   * @param border the border
   * @param xlimit the xlimit
   * @param ylimit the ylimit
   * @return true, if outside the border (should be ignored)
   */
  private static boolean outsideBorder(BasePreprocessedPeakResult result, final float border,
      final float xlimit, final float ylimit) {
    // Respect the border used in BenchmarkSpotFilter?
    // Just implement a hard border
    return (result.getX() < border || result.getX() > xlimit || result.getY() < border
        || result.getY() > ylimit);
  }

  /**
   * Save the results to memory.
   *
   * @param filterResults The results from running the filter (or null)
   * @param filter the filter
   */
  private void saveResults(PreprocessedPeakResult[] filterResults, DirectFilter filter) {
    final MemoryPeakResults newResults = createResults(filterResults, filter, true);
    MemoryPeakResults.addResults(newResults);
  }

  /**
   * Create peak results.
   *
   * @param filterResults The results from running the filter (or null)
   * @param filter the filter
   */
  private MemoryPeakResults createResults(PreprocessedPeakResult[] filterResults,
      DirectFilter filter, boolean withBorder) {
    if (filterResults == null) {
      final MultiPathFilter multiPathFilter = createMpf(filter, defaultMinimalFilter);
      filterResults = filterResults(multiPathFilter);
    }

    final MemoryPeakResults newResults = new MemoryPeakResults();
    newResults.copySettings(this.results);
    newResults.setName(TITLE);

    if (withBorder) {
      // To produce the same results as the PeakFit plugin we must implement the border
      // functionality used in the FitWorker. This respects the border of the spot filter.
      final FitEngineConfiguration config = new FitEngineConfiguration();
      updateAllConfiguration(config);
      final MaximaSpotFilter spotFilter = config.createSpotFilter();
      final int border = spotFilter.getBorder();
      final Rectangle bounds = getBounds();
      final int borderLimitX = bounds.x + bounds.width - border;
      final int borderLimitY = bounds.y + bounds.height - border;

      for (final PreprocessedPeakResult spot : filterResults) {
        if (spot.getX() > border && spot.getX() < borderLimitX && spot.getY() > border
            && spot.getY() < borderLimitY) {
          final double[] p = spot.toGaussian2DParameters();
          final float[] params = new float[p.length];
          for (int j = 0; j < p.length; j++) {
            params[j] = (float) p[j];
          }
          final int frame = spot.getFrame();
          final int origX = (int) p[Gaussian2DFunction.X_POSITION];
          final int origY = (int) p[Gaussian2DFunction.Y_POSITION];

          newResults.add(frame, origX, origY, 0, 0, spot.getNoise(), spot.getMeanSignal(), params,
              null);
        }
      }
    } else {
      for (final PreprocessedPeakResult spot : filterResults) {
        final double[] p = spot.toGaussian2DParameters();
        final float[] params = new float[p.length];
        for (int j = 0; j < p.length; j++) {
          params[j] = (float) p[j];
        }
        final int frame = spot.getFrame();
        final int origX = (int) p[Gaussian2DFunction.X_POSITION];
        final int origY = (int) p[Gaussian2DFunction.Y_POSITION];

        newResults.add(frame, origX, origY, 0, 0, spot.getNoise(), spot.getMeanSignal(), params,
            null);
      }
    }

    return newResults;
  }

  private PreprocessedPeakResult[] filterResults(final MultiPathFilter multiPathFilter) {
    return multiPathFilter.filter(fitResultData.resultsList,
        createFailCounter(settings.failCount), true, coordinateStore);
  }

  private CoordinateStore createCoordinateStore() {
    return createCoordinateStore(settings.duplicateDistance);
  }

  private CoordinateStore createCoordinateStore(double duplicateDistance) {
    getBounds();
    duplicateDistance *= distanceScallingFactor;
    return CoordinateStoreFactory.create(bounds.x, bounds.y, bounds.width, bounds.height,
        duplicateDistance);
  }

  private Rectangle getBounds() {
    if (bounds == null) {
      bounds = createBounds();
    }
    return bounds;
  }

  private synchronized Rectangle createBounds() {
    if (bounds == null) {
      if (settings.duplicateDistanceAbsolute) {
        distanceScallingFactor = 1;
      } else {
        // The duplicate distance is scaled
        distanceScallingFactor = BenchmarkSpotFit.getFitEngineConfiguration().getHwhmMax();
      }

      final ImagePlus imp = CreateData.getImage();
      if (imp != null) {
        return new Rectangle(0, 0, imp.getWidth(), imp.getHeight());
      }
      return results.getBounds(true);
    }
    return bounds;
  }

  /**
   * Creates the fail counter.
   *
   * @param failCount the fail count
   * @return the consecutive fail counter
   */
  private static ConsecutiveFailCounter createFailCounter(final int failCount) {
    // TODO - What about pass rate?
    return ConsecutiveFailCounter.create(failCount);
  }

  /**
   * Gets the distance in pixels from the latest read of the fit results.
   *
   * @return the distance in pixels
   */
  static double getDistanceInPixels() {
    return fitResultDataCache.get().distanceInPixels;
  }

  /**
   * Gets the lower distance in pixels from the latest read of the fit results.
   *
   * @return the lower distance in pixels
   */
  static double getLowerDistanceInPixels() {
    return fitResultDataCache.get().lowerDistanceInPixels;
  }


  /**
   * Gets the signal factor from the latest read of the fit results.
   *
   * @return the signal factor
   */
  static double getSignalFactor() {
    return fitResultDataCache.get().signalFactor;
  }

  /**
   * Gets the lower signal factor from the latest read of the fit results.
   *
   * @return the lower signal factor
   */
  static double getLowerSignalFactor() {
    return fitResultDataCache.get().lowerSignalFactor;
  }

  /**
   * Gets the last fitting id that was analysed.
   *
   * @return the last fitting id
   */
  static int getLastFittingId() {
    return fitResultDataCache.get().fittingId;
  }
}
