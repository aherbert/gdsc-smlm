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

import uk.ac.sussex.gdsc.core.ij.BufferedTextWindow;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot;
import uk.ac.sussex.gdsc.core.ij.ImageJTrackProgress;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.Plot2;
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
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.RampedScore;
import uk.ac.sussex.gdsc.core.utils.SettingsList;
import uk.ac.sussex.gdsc.core.utils.StoredData;
import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.TurboList;
import uk.ac.sussex.gdsc.core.utils.UnicodeReader;
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
import uk.ac.sussex.gdsc.smlm.ij.plugins.BenchmarkSpotFit.FilterCandidates;
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

import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.procedure.TIntObjectProcedure;
import gnu.trove.set.hash.TIntHashSet;

import ij.IJ;
import ij.ImageListener;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Macro;
import ij.Prefs;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.ImageWindow;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.plugin.PlugIn;
import ij.plugin.frame.Recorder;
import ij.text.TextWindow;

import org.apache.commons.lang3.concurrent.ConcurrentRuntimeException;
import org.apache.commons.lang3.time.DurationFormatUtils;
import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.Well44497b;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.util.FastMath;

import java.awt.AWTEvent;
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
import java.io.Serializable;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Consumer;
import java.util.logging.Level;
import java.util.logging.Logger;

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
  private static AtomicReference<TextWindow> resultsWindowRef = new AtomicReference<>();
  private static AtomicReference<TextWindow> summaryWindowRef = new AtomicReference<>();
  private static AtomicReference<TextWindow> sensitivityWindowRef = new AtomicReference<>();
  private static AtomicReference<TextWindow> gaWindowRef = new AtomicReference<>();
  private static AtomicReference<TextWindow> componentAnalysisWindowRef = new AtomicReference<>();
  private static int failCount = 1;
  private static int minFailCount;
  private static int maxFailCount = 10;

  /**
   * A reference to the results window if not in headless mode. This is used with a
   * BufferedTextWindow to output the large results set.
   */
  private TextWindow resultsWindow;
  /** The genetic analysis summary results output. */
  private Consumer<String> gaWindow;

  // This can be used during filtering.
  // However the min filter is not used to determine if candidates are valid (that is the primary
  // filter).
  // It is used to store poor estimates during fitting. So we can set it to null.
  private static final DirectFilter minimalFilter = null;
  private static double sResidualsThreshold = 0.3;
  private static double minResidualsThreshold = 0.1;
  private static double maxResidualsThreshold = 0.6;
  private double residualsThreshold = 1; // Disabled
  private static double duplicateDistance;
  // This is a flag that is passed around but is only set once, i.e.
  // analysis is done using either absolute or relative distances.
  private static boolean duplicateDistanceAbsolute = true;
  private static double minDuplicateDistance;
  private static double maxDuplicateDistance = 5;
  private static boolean reset = true;
  private static boolean showResultsTable;
  private static boolean showSummaryTable = true;
  private static boolean clearTables;
  private static final String KEY_FILTER_FILENAME = "gdsc.filteranalysis.filterfilename";
  private static final String KEY_FILTERSET_FILENAME = "gdsc.filteranalysis.filtersetfilename";
  private static final String KEY_TEMPLATE_FILENAME = "gdsc.filteranalysis.templatefilename";
  private static String filterFilename = Prefs.get(KEY_FILTER_FILENAME, "");
  private static String filterSetFilename = Prefs.get(KEY_FILTERSET_FILENAME, "");
  private static String templateFilename = Prefs.get(KEY_TEMPLATE_FILENAME, "");
  private static int summaryTopN;
  private static double summaryDepth = 500;
  private static int plotTopN;
  private static boolean saveBestFilter;
  private static boolean saveTemplate;
  private static boolean calculateSensitivity;
  private static double delta = 0.1;
  private static int criteriaIndex;
  private static double criteriaLimit = 0.95;
  private double minCriteria;
  private boolean invertCriteria;
  private static int scoreIndex;
  private boolean invertScore;
  private static double upperMatchDistance = 100;
  private static double partialMatchDistance = 33;
  /** The distance in pixels. */
  static double distanceInPixels;
  /** The lower distance in pixels. */
  static double lowerDistanceInPixels;
  private static double upperSignalFactor = 100;
  private static double partialSignalFactor = 50;
  /** The signal factor. */
  static double signalFactor;
  /** The lower signal factor. */
  static double lowerSignalFactor;
  private static boolean depthRecallAnalysis = true;
  private static boolean scoreAnalysis = true;
  private static final String[] COMPONENT_ANALYSIS =
      {"None", "Best Ranked", "Ranked", "Best All", "All"};
  private static int componentAnalysis = 3;

  private static final String[] EVOLVE =
      {"None", "Genetic Algorithm", "Range Search", "Enrichment Search", "Step Search"};
  private static int evolve;
  private static boolean repeatEvolve;
  private static int rangeSearchWidth = 2;
  private static double rangeSearchReduce = 0.3;
  private static int maxIterations = 30;
  private static int refinementMode = SearchSpace.RefinementMode.SINGLE_DIMENSION.ordinal();
  private static int enrichmentSamples = 5000;
  private static int seedSize = 5000;
  private static double enrichmentFraction = 0.2;
  private static double enrichmentPadding = 0.1;

  private static final String[] SEARCH =
      {"Range Search", "Enrichment Search", "Step Search", "Enumerate"};
  private static int searchParam = 3;
  private static boolean repeatSearch;
  private static int pRangeSearchWidth = 2;
  private static double pRangeSearchReduce = 0.3;
  private static int pMaxIterations = 30;
  private static int pRefinementMode = SearchSpace.RefinementMode.MULTI_DIMENSION.ordinal();
  private static int pEnrichmentSamples = 500;
  private static int pSeedSize = 500;
  private static double pEnrichmentFraction = 0.2;
  private static double pEnrichmentPadding = 0.1;
  private static int pConvergedCount = 2;

  private static TIntObjectHashMap<boolean[]> searchRangeMap = new TIntObjectHashMap<>();
  private static TIntObjectHashMap<double[]> stepSizeMap = new TIntObjectHashMap<>();

  private static boolean showTP;
  private static boolean showFP;
  private static boolean showFN;

  private static int populationSize = 5000;
  private static int failureLimit = 5;
  private static double tolerance = 1e-4;
  private static int convergedCount = 2;
  private static double crossoverRate = 1;
  private static double meanChildren = 2;
  private static double mutationRate = 1;
  private static double selectionFraction = 0.2;
  private static boolean rampedSelection = true;
  private static boolean saveOption;
  private static double iterationScoreTolerance = 1e-4;
  private static double iterationFilterTolerance = 1e-3;
  private static boolean iterationCompareResults;
  private static double iterationCompareDistance = 0.1;
  private static int iterationMaxIterations = 10;
  private static double iterationMinRangeReduction = 0.2;
  private static int iterationMinRangeReductionIteration = 5;
  private static boolean iterationConvergeBeforeRefit;

  // For the template example
  private static int nNo = 2;
  private static int nLow = 4;
  private static int nHigh = 4;

  private static String resultsTitle;
  private String resultsPrefix;
  private String resultsPrefix2;
  private String limitFailCount;
  private static String resultsPrefix3;
  private static String limitRange;

  private static ArrayList<NamedPlot> plots = new ArrayList<>();
  private static HashMap<String, ComplexFilterScore> bestFilter = new HashMap<>();
  private static LinkedList<String> bestFilterOrder = new LinkedList<>();

  private static HashMap<String, ComplexFilterScore> iterBestFilter;

  private static boolean reUseFilters = true;
  private static boolean expandFilters;
  private static String oldFilename = "";
  private static long lastModified;
  private static List<FilterSet> filterList;
  /** The last id. */
  static int lastId;
  private static TIntObjectHashMap<UniqueIdPeakResult[]> actualCoordinates;
  private static MultiPathFitResults[] resultsList;
  // private static MultiPathFitResults[] clonedResultsList = null;
  private static int matches;
  private static int fittedResults;
  private static int totalResults;
  private static int notDuplicateCount;
  private static int newResultCount;
  private static int maxUniqueId;
  private static int nActual;
  private static StoredData depthStats;
  private static StoredData depthFitStats;
  private static StoredDataStatistics signalFactorStats;
  private static StoredDataStatistics distanceStats;

  private final boolean isHeadless;
  private boolean debug;
  private CoordinateStore coordinateStore;

  // Used to tile plot windows
  private final WindowOrganiser wo = new WindowOrganiser();

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
      return BenchmarkSpotFit.getSignalFactor(peakResult.getSignal(), peak.getIntensity());
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
   * Used to allow multi-threading of the scoring the fit results.
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

    public FitResultsWorker(BlockingQueue<Job> jobs, List<MultiPathFitResults> syncResults,
        double matchDistance, RampedScore distanceScore, RampedScore signalScore,
        AtomicInteger uniqueId, CoordinateStore coordinateStore) {
      this.jobs = jobs;
      this.results = syncResults;
      this.matchDistance = matchDistance;
      this.distanceScore = distanceScore;
      this.signalScore = signalScore;
      this.uniqueId = uniqueId;
      this.coordinateStore = coordinateStore;

      depthStats = new StoredData();
      depthFitStats = new StoredData();
      signalFactorStats = new StoredDataStatistics();
      distanceStats = new StoredDataStatistics();

      checkBorder = (BenchmarkSpotFilter.lastAnalysisBorder != null
          && BenchmarkSpotFilter.lastAnalysisBorder.x != 0);
      if (checkBorder) {
        final Rectangle lastAnalysisBorder = BenchmarkSpotFilter.lastAnalysisBorder;
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
        while (true) {
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
        System.out.println(ex.toString());
        throw new RuntimeException(ex);
      } finally {
        finished = true;
      }
    }

    private void run(Job job) {
      if (ImageJUtils.isInterrupted()) {
        finished = true;
        return;
      }

      showProgress();

      final int frame = job.frame;
      final FilterCandidates result = job.candidates;

      depthStats.add(result.zPosition);

      UniqueIdPeakResult[] actual = getCoordinates(actualCoordinates, frame);
      final int nActual = actual.length;
      final boolean[] matched = new boolean[nActual];
      actual = filterUsingBorder(actual);
      // We could use distanceInPixels for the resolution. Using a bigger size may allow the
      // different fit locations to be in the same cell and so the grid manager can use it's cache.
      final double resolution = 2 * distanceInPixels;
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
        boolean fitted = score(index, 1, fitResult.getSingleFitResult(), resultGrid, matched);
        fitted |= score(index, 2, fitResult.getDoubletFitResult(), resultGrid, matched);
        fitted |= score(index, 3, fitResult.getMultiFitResult(), resultGrid, matched);
        fitted |= score(index, 4, fitResult.getMultiDoubletFitResult(), resultGrid, matched);

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
    }

    /**
     * Score the new results in the fit result.
     *
     * @param n the n
     * @param set the set
     * @param fitResult the fit result
     * @param resultGrid the result grid
     * @param matched array of actual results that have been matched
     * @return true, if the fit status was ok
     */
    private boolean score(int n, int set,
        final uk.ac.sussex.gdsc.smlm.results.filter.MultiPathFitResult.FitResult fitResult,
        final PeakResultGridManager resultGrid, final boolean[] matched) {
      if (fitResult != null && fitResult.status == 0) {
        // Get the new results
        for (int i = 0; i < fitResult.results.length; i++) {
          final BasePreprocessedPeakResult peak = (BasePreprocessedPeakResult) fitResult.results[i];
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
                final double signalFactor =
                    BenchmarkSpotFit.getSignalFactor(peak.getSignal(), actual[j].getIntensity());
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

        return (fitResult.results.length != 0);
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

  private static ArrayList<FilterResult> scores = new ArrayList<>();

  private static String[] COLUMNS = {
      // Scores using integer scoring
      "TP", "FP", "FN", "Precision", "Recall", "F1", "Jaccard",
      // Scores using fractional scoring
      "fTP", "fFP", "fFN", "fPrecision", "fRecall", "fF1", "fJaccard",};

  private static boolean[] showColumns;
  private boolean requireIntegerResults;

  static {
    showColumns = new boolean[COLUMNS.length];
    Arrays.fill(showColumns, true);
    // showColumns[0] = false; // nP

    // Use the precision as criteria to ensure a set confidence on results labelled as true
    criteriaIndex = COLUMNS.length - 4;
    // Score using Jaccard
    scoreIndex = COLUMNS.length - 1;

    if (TextUtils.isNullOrEmpty(templateFilename)) {
      final String currentUsersHomeDir = System.getProperty("user.home");
      templateFilename =
          currentUsersHomeDir + File.separator + "gdsc.smlm" + File.separator + "template";
    }
  }

  private CreateData.SimulationParameters simulationParameters;
  private MemoryPeakResults results;
  private boolean extraOptions;

  /**
   * Instantiates a new benchmark filter analysis.
   */
  public BenchmarkFilterAnalysis() {
    isHeadless = java.awt.GraphicsEnvironment.isHeadless();
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    simulationParameters = CreateData.simulationParameters;
    if (simulationParameters == null) {
      IJ.error(TITLE, "No benchmark spot parameters in memory");
      return;
    }
    results = CreateData.getResults();
    if (results == null) {
      IJ.error(TITLE, "No benchmark results in memory");
      return;
    }

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
    if (BenchmarkSpotFit.fitResults == null) {
      if (!silent) {
        IJ.error(TITLE, "No benchmark fitting results in memory");
      }
      return true;
    }
    if (BenchmarkSpotFit.lastId != simulationParameters.id) {
      if (!silent) {
        IJ.error(TITLE, "Update the benchmark spot fitting for the latest simulation");
      }
      return true;
    }
    if (BenchmarkSpotFit.lastFilterId != BenchmarkSpotFilter.filterResult.id) {
      if (!silent) {
        IJ.error(TITLE, "Update the benchmark spot fitting for the latest filter");
      }
      return true;
    }
    return false;
  }

  private boolean loadFitResults() {
    resultsList = readResults();
    if (resultsList == null) {
      IJ.error(TITLE, "No results could be loaded");
      return false;
    }
    return true;
  }

  private void iterate() {
    // If this is run again immediately then provide options for reporting the results
    if (iterBestFilter != null && iterBestFilter == bestFilter) {
      final GenericDialog gd = new GenericDialog(TITLE);
      gd.enableYesNoCancel();
      gd.addMessage("Iteration results are held in memory.\n \nReport these results?");
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
    // TODO - collect this in the iteration dialog
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
      if (!fit.finished) {
        // The plugin did not complete
        return;
      }
      resetParametersFromFitting();
    }
    if (invalidBenchmarkSpotFitResults(false)) {
      return;
    }
    if (BenchmarkSpotFit.stopWatch != null) {
      time += BenchmarkSpotFit.stopWatch.getTime();
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
    iterBestFilter = null;

    ImageJUtils.log(TITLE + " Iterating ...");

    final IterationConvergenceChecker checker = new IterationConvergenceChecker(current);

    // Iterate ...
    boolean outerConverged = false;
    int outerIteration = 1;
    double outerRangeReduction = 1;
    while (!outerConverged) {
      if (iterationConvergeBeforeRefit) {
        // Optional inner loop so that the non-filter and filter parameters converge
        // before a refit
        boolean innerConverged = false;
        int innerIteration = 0;
        double innerRangeReduction = 1;
        if (iterationMinRangeReduction < 1) {
          // Linear interpolate down to the min range reduction
          innerRangeReduction = MathUtils.max(iterationMinRangeReduction,
              MathUtils.interpolateY(0, 1, iterationMinRangeReductionIteration,
                  iterationMinRangeReduction, innerIteration++));
        }
        // This would make the range too small...
        // innerRangeReduction *= outerRangeReduction;
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
      fit.run(current.result.filter, residualsThreshold, failCount, duplicateDistance,
          duplicateDistanceAbsolute);
      if (invalidBenchmarkSpotFitResults(false)) {
        return;
      }
      if (!loadFitResults()) {
        return;
      }

      // Reduce the range over which the filter parameters are searched. Note that this range
      // is centred around the current optimum.
      if (iterationMinRangeReduction < 1) {
        // Linear interpolate down to the min range reduction
        outerRangeReduction = MathUtils.max(iterationMinRangeReduction, MathUtils.interpolateY(0, 1,
            iterationMinRangeReductionIteration, iterationMinRangeReduction, outerIteration++));
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
      iterBestFilter = bestFilter;
    }

    time += iterationStopWatch.getTime();
    IJ.log("Iteration analysis time : " + DurationFormatUtils.formatDurationHMS(time));

    IJ.showStatus("Finished");
  }

  private void resetParametersFromFitting() {
    failCount = BenchmarkSpotFit.config.getFailuresLimit();
    duplicateDistance = BenchmarkSpotFit.config.getDuplicateDistance();
    duplicateDistanceAbsolute = BenchmarkSpotFit.config.getDuplicateDistanceAbsolute();
    residualsThreshold = sResidualsThreshold =
        (BenchmarkSpotFit.computeDoublets) ? BenchmarkSpotFit.multiFilter.residualsThreshold : 1;
  }

  private double[] createParameters() {
    // Ignore the duplicate distance absolute as this is just for convergence checking
    return new double[] {failCount, residualsThreshold, duplicateDistance};
  }

  private void reportIterationResults() {
    residualsThreshold = sResidualsThreshold;
    if (!showReportDialog()) {
      return;
    }
    reportResults(false);
  }

  private static DirectFilter scoreFilter;
  private static int scoreFailCount;
  private static double scoreResidualsThreshold;
  private static double scoreDuplicateDistance = -1;

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
    failCount = scoreFailCount;
    residualsThreshold = scoreResidualsThreshold;
    duplicateDistance = scoreDuplicateDistance;

    // Create a dummy result, the filter will be rescored in reportResults(...)
    final FilterScoreResult sr = new FilterScoreResult(0, 0, scoreFilter, "");
    final ComplexFilterScore newFilterScore = new ComplexFilterScore(sr, null, "", 0, "", 0);

    // Report to summary window
    reportResults(true, newFilterScore);

    // Reset the variable used for scoring
    failCount = (int) stash[0];
    residualsThreshold = stash[1];
    duplicateDistance = stash[2];
  }

  private void getScoreFilter() {
    if (scoreFilter == null) {
      // Reset to default only on first run
      if (scoreDuplicateDistance == -1) {
        scoreFailCount = failCount;
        scoreDuplicateDistance = duplicateDistance;
        scoreResidualsThreshold = residualsThreshold;
      }

      // Use the best result if we have one
      final FilterResult r = getBestResult();
      if (r != null) {
        scoreFilter = r.getFilter();
        scoreFailCount = r.failCount;
        scoreResidualsThreshold = r.residualsThreshold;
      } else {
        // Default to the fit config settings
        final FitConfiguration tmp = new FitConfiguration();
        // So we get a MultiFilter2
        tmp.setPrecisionMethod(PrecisionMethod.MORTENSEN_LOCAL_BACKGROUND);
        scoreFilter = tmp.getDefaultSmartFilter();
      }
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

    IJ.showStatus("Finished");
  }

  private void optimiseParameters() {
    final FilterResult fr = getBestResult();
    if (fr == null) {
      IJ.error(TITLE, "No filter scores in memory");
      return;
    }

    if (!showDialog(FLAG_OPTIMISE_PARAMS)) {
      return;
    }

    analyseParameters(false, fr.filterScore, 0);

    IJ.showStatus("Finished");
  }

  @SuppressWarnings("unchecked")
  private List<FilterSet> readFilterSets() {
    if (extraOptions && BenchmarkSpotFit.multiFilter != null) {
      final IDirectFilter f = BenchmarkSpotFit.multiFilter.getFilter();
      if (f instanceof DirectFilter) {
        final GenericDialog gd = new GenericDialog(TITLE);
        gd.addMessage("Use an identical filter to " + BenchmarkSpotFit.TITLE);
        gd.enableYesNoCancel();
        gd.hideCancelButton();
        gd.showDialog();
        if (gd.wasOKed()) {
          setLastFile(null);
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

    GUIFilterSettings filterSettings = SettingsManager.readGuiFilterSettings(0);

    final String filename =
        ImageJUtils.getFilename("Filter_File", filterSettings.getFilterSetFilename());
    if (filename != null) {
      IJ.showStatus("Reading filters ...");
      filterSettings = filterSettings.toBuilder().setFilterSetFilename(filename).build();

      // Allow the filters to be cached
      if (isSameFile(filename)) {
        final GenericDialog gd = new GenericDialog(TITLE);
        gd.hideCancelButton();
        gd.addMessage("The same filter file was selected.");
        gd.addCheckbox("Re-use_filters", reUseFilters);
        gd.showDialog();
        if (!gd.wasCanceled()) {
          reUseFilters = gd.getNextBoolean();
          if (reUseFilters) {
            SettingsManager.writeSettings(filterSettings);
            return filterList;
          }
        }
      }

      setLastFile(null);
      try (BufferedReader input =
          new BufferedReader(new UnicodeReader(new FileInputStream(filename), null))) {
        // Use the instance so we can catch the exception
        final Object o = FilterXStreamUtils.getXStreamInstance().fromXML(input);
        input.close();

        if (!(o instanceof List<?>)) {
          IJ.log("No filter sets defined in the specified file: " + filename);
          return null;
        }
        SettingsManager.writeSettings(filterSettings);
        final List<FilterSet> filterSets = (List<FilterSet>) o;

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

        filterList = filterSets;

        // Option to enumerate filters
        expandFilters();

        setLastFile(filename);
        return filterList;
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
   */
  private static void expandFilters() {
    // Do not clear these when reading a new set of filters.
    // The filters may be the same with slight modifications and so it is useful to keep the last
    // settings.
    // searchRangeMap.clear();
    // stepSizeMap.clear();

    final long[] expanded = new long[filterList.size()];
    final String[] name = new String[expanded.length];
    int count = 0;
    boolean doIt = false;
    for (final FilterSet filterSet : filterList) {
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
      return;
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
    gd.addCheckbox("Expand_filters", expandFilters);
    gd.showDialog();
    if (!gd.wasCanceled()) {
      if (!(expandFilters = gd.getNextBoolean())) {
        return;
      }
    }

    IJ.showStatus("Expanding filters ...");

    final List<FilterSet> filterList2 = new ArrayList<>(filterList.size());
    for (final FilterSet filterSet : filterList) {
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

    filterList = filterList2;
    ImageJUtils.log("Expanded input to %d filters in %s", countFilters(filterList),
        TextUtils.pleural(filterList.size(), "set"));
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

      if (Double.isNaN(inc) || Double.isInfinite(inc)) {
        continue;
      }
      if (Double.isNaN(parameters[i]) || Double.isInfinite(parameters[i])) {
        continue;
      }
      if (Double.isNaN(parameters2[i]) || Double.isInfinite(parameters2[i])) {
        continue;
      }
      if (parameters[i] > parameters2[i]) {
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
   * @return true, if is same file
   */
  public boolean isSameFile(String filename) {
    if (filterList == null) {
      return false;
    }
    if (filename.equals(oldFilename)) {
      try {
        final File f = new File(filename);
        if (lastModified == f.lastModified()) {
          return true;
        }
      } catch (final Exception ex) {
        // Ignore
      }
    }
    return false;
  }

  private static void setLastFile(String filename) {
    oldFilename = filename;
    if (oldFilename != null) {
      try {
        final File f = new File(filename);
        lastModified = f.lastModified();
      } catch (final Exception ex) {
        lastModified = 0;
      }
    }
  }

  private static SettingsList lastReadResultsSettings;
  private static double lastDuplicateDistance = -1;
  private static boolean lastDuplicateDistanceAbsolute;

  private MultiPathFitResults[] readResults() {
    boolean update = resultsList == null; // XXX set to true when debugging
    if (lastId != BenchmarkSpotFit.fitResultsId) {
      if (lastId == 0) {
        // Copy the settings from the fitter if this is the first run
        failCount = BenchmarkSpotFit.config.getFailuresLimit();
        duplicateDistance = BenchmarkSpotFit.config.getDuplicateDistance();
        duplicateDistanceAbsolute = BenchmarkSpotFit.config.getDuplicateDistanceAbsolute();
        sResidualsThreshold =
            (BenchmarkSpotFit.computeDoublets) ? BenchmarkSpotFit.multiFilter.residualsThreshold
                : 1;
      }
      lastId = BenchmarkSpotFit.fitResultsId;
      update = true;
      actualCoordinates = getCoordinates(results);
    }

    final SettingsList settings = new SettingsList(partialMatchDistance, upperMatchDistance,
        partialSignalFactor, upperSignalFactor);
    final boolean equalScoreSettings = settings.equals(lastReadResultsSettings);

    if (update || !equalScoreSettings || lastDuplicateDistance != duplicateDistance
        || lastDuplicateDistanceAbsolute != duplicateDistanceAbsolute) {
      IJ.showStatus("Reading results ...");

      // Only cache results for the same score analysis settings.
      // This functionality is for choosing the optimum filter for the given scoring metric.
      if (!equalScoreSettings) {
        scores.clear();
      }

      lastReadResultsSettings = settings;
      lastDuplicateDistance = duplicateDistance;
      lastDuplicateDistanceAbsolute = duplicateDistanceAbsolute;
      depthStats = null;
      depthFitStats = null;
      signalFactorStats = null;
      distanceStats = null;
      matches = 0;
      fittedResults = 0;
      totalResults = 0;
      notDuplicateCount = 0;
      newResultCount = 0;
      maxUniqueId = 0;
      nActual = 0;

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
          new RampedScore(BenchmarkSpotFit.distanceInPixels * partialMatchDistance / 100.0,
              BenchmarkSpotFit.distanceInPixels * upperMatchDistance / 100.0);
      lowerDistanceInPixels = distanceScore.lower;
      distanceInPixels = distanceScore.upper;
      final double matchDistance = distanceInPixels * distanceInPixels;

      resultsPrefix3 =
          "\t" + MathUtils.rounded(distanceScore.lower * simulationParameters.pixelPitch) + "\t"
              + MathUtils.rounded(distanceScore.upper * simulationParameters.pixelPitch);
      limitRange = ", d=" + MathUtils.rounded(distanceScore.lower * simulationParameters.pixelPitch)
          + "-" + MathUtils.rounded(distanceScore.upper * simulationParameters.pixelPitch);

      // Signal factor must be greater than 1
      final RampedScore signalScore;
      if (BenchmarkSpotFit.signalFactor > 0 && upperSignalFactor > 0) {
        signalScore = new RampedScore(BenchmarkSpotFit.signalFactor * partialSignalFactor / 100.0,
            BenchmarkSpotFit.signalFactor * upperSignalFactor / 100.0);
        lowerSignalFactor = signalScore.lower;
        signalFactor = signalScore.upper;
        resultsPrefix3 += "\t" + MathUtils.rounded(signalScore.lower) + "\t"
            + MathUtils.rounded(signalScore.upper);
        limitRange += ", s=" + MathUtils.rounded(signalScore.lower) + "-"
            + MathUtils.rounded(signalScore.upper);
      } else {
        signalScore = null;
        resultsPrefix3 += "\t0\t0";
        lowerSignalFactor = signalFactor = 0;
      }

      // Store all the results
      final ArrayList<MultiPathFitResults> results =
          new ArrayList<>(BenchmarkSpotFit.fitResults.size());
      final List<MultiPathFitResults> syncResults = Collections.synchronizedList(results);

      // This could be multi-threaded ...
      final int nThreads = getThreads(BenchmarkSpotFit.fitResults.size());
      final BlockingQueue<Job> jobs = new ArrayBlockingQueue<>(nThreads * 2);
      final List<FitResultsWorker> workers = new LinkedList<>();
      final List<Thread> threads = new LinkedList<>();
      final AtomicInteger uniqueId = new AtomicInteger();
      final CoordinateStore coordinateStore = createCoordinateStore();
      for (int i = 0; i < nThreads; i++) {
        final FitResultsWorker worker = new FitResultsWorker(jobs, syncResults, matchDistance,
            distanceScore, signalScore, uniqueId, coordinateStore.newInstance());
        final Thread t = new Thread(worker);
        workers.add(worker);
        threads.add(t);
        t.start();
      }

      totalProgress = BenchmarkSpotFit.fitResults.size();
      stepProgress = ImageJUtils.getProgressInterval(totalProgress);
      progress = 0;
      BenchmarkSpotFit.fitResults.forEachEntry((frame, candidates) -> {
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
          matches += worker.matches;
          fittedResults += worker.included;
          totalResults += worker.total;
          notDuplicateCount += worker.notDuplicateCount;
          newResultCount += worker.newResultCount;
          nActual += worker.includedActual;
          if (i == 0) {
            depthStats = worker.depthStats;
            depthFitStats = worker.depthFitStats;
            signalFactorStats = worker.signalFactorStats;
            distanceStats = worker.distanceStats;
          } else {
            depthStats.add(worker.depthStats);
            depthFitStats.add(worker.depthFitStats);
            signalFactorStats.add(worker.signalFactorStats);
            distanceStats.add(worker.distanceStats);
          }
        } catch (final InterruptedException ex) {
          ex.printStackTrace();
        }
      }
      threads.clear();
      IJ.showProgress(1);
      IJ.showStatus("");

      maxUniqueId = uniqueId.get();

      resultsList = results.toArray(new MultiPathFitResults[results.size()]);

      Arrays.sort(resultsList, new Comparator<MultiPathFitResults>() {
        @Override
        public int compare(MultiPathFitResults o1, MultiPathFitResults o2) {
          return o1.getFrame() - o2.getFrame();
        }
      });
    }

    // In case a previous run was interrupted
    if (resultsList != null) {
      MultiPathFilter.resetValidationFlag(resultsList);
    }

    return resultsList;
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
      final TurboList<PeakResult> tmp = new TurboList<>();
      // Add the results to the lists
      results.forEach((PeakResultProcedure) result -> {
        if (counter.advanceAndReset(result.getFrame()) && !tmp.isEmpty()) {
          coords.put(counter.previousFrame(), tmp.toArray(new UniqueIdPeakResult[tmp.size()]));
          tmp.clear();
        }
        tmp.add(new UniqueIdPeakResult(tmp.size(), uniqueId.getAndIncrement(), result));
      });

      if (!tmp.isEmpty()) {
        coords.put(counter.previousFrame(), tmp.toArray(new UniqueIdPeakResult[tmp.size()]));
      }
    }
    return coords;
  }

  private static final int FLAG_OPTIMISE_FILTER = 1;
  private static final int FLAG_OPTIMISE_PARAMS = 2;

  private boolean showDialog(int optimiseParameters) {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);

    final boolean showOptimiseFilter = (optimiseParameters & FLAG_OPTIMISE_FILTER) != 0;
    final boolean showOptimiseParams = (optimiseParameters & FLAG_OPTIMISE_PARAMS) != 0;

    addSimulationData(gd);

    // TODO - Make minimal filter configurable?

    gd.addSlider("Fail_count", 0, 20, failCount);
    if (showOptimiseParams) {
      gd.addNumericField("Min_fail_count", minFailCount, 0);
      gd.addNumericField("Max_fail_count", maxFailCount, 0);
    }
    if (BenchmarkSpotFit.computeDoublets) {
      gd.addSlider("Residuals_threshold", 0.01, 1, sResidualsThreshold);
      if (showOptimiseParams) {
        gd.addNumericField("Min_residuals_threshold", minResidualsThreshold, 2);
        gd.addNumericField("Max_residuals_threshold", maxResidualsThreshold, 2);
      }
    }
    final FitEngineConfiguration tmp = new FitEngineConfiguration();
    tmp.setDuplicateDistance(duplicateDistance);
    tmp.setDuplicateDistanceAbsolute(duplicateDistanceAbsolute);
    PeakFit.addDuplicateDistanceOptions(gd, new PeakFit.SimpleFitEngineConfigurationProvider(tmp));
    if (showOptimiseParams) {
      gd.addNumericField("Min_duplicate_distance", minDuplicateDistance, 2);
      gd.addNumericField("Max_duplicate_distance", maxDuplicateDistance, 2);
    }
    gd.addCheckbox("Reset", reset);
    gd.addCheckbox("Show_table", showResultsTable);
    gd.addCheckbox("Show_summary", showSummaryTable);
    gd.addCheckbox("Clear_tables", clearTables);
    gd.addSlider("Summary_top_n", 0, 20, summaryTopN);
    gd.addNumericField("Summary_depth (nm)", summaryDepth, 0);
    gd.addSlider("Plot_top_n", 0, 20, plotTopN);
    gd.addCheckbox("Save_best_filter", saveBestFilter);
    gd.addCheckbox("Save_template", saveTemplate);
    gd.addCheckbox("Calculate_sensitivity", calculateSensitivity);
    gd.addSlider("Delta", 0.01, 1, delta);
    gd.addMessage("Match scoring");
    gd.addChoice("Criteria", COLUMNS, COLUMNS[criteriaIndex]);
    gd.addNumericField("Criteria_limit", criteriaLimit, 4);
    gd.addChoice("Score", COLUMNS, COLUMNS[scoreIndex]);
    gd.addMessage(String.format("Fitting match distance = %s nm; signal factor = %s",
        MathUtils.rounded(BenchmarkSpotFit.distanceInPixels * simulationParameters.pixelPitch),
        MathUtils.rounded(BenchmarkSpotFit.signalFactor)));
    gd.addSlider("Upper_match_distance (%)", 0, 100, upperMatchDistance);
    gd.addSlider("Partial_match_distance (%)", 0, 100, partialMatchDistance);
    gd.addSlider("Upper_signal_factor (%)", 0, 100, upperSignalFactor);
    gd.addSlider("Partial_signal_factor (%)", 0, 100, partialSignalFactor);
    if (!simulationParameters.fixedDepth) {
      gd.addCheckbox("Depth_recall_analysis", depthRecallAnalysis);
    }
    gd.addCheckbox("Score_analysis", scoreAnalysis);
    gd.addChoice("Component_analysis", COMPONENT_ANALYSIS, COMPONENT_ANALYSIS[componentAnalysis]);
    if (showOptimiseFilter) {
      gd.addChoice("Evolve", EVOLVE, EVOLVE[evolve]);
      gd.addCheckbox("Repeat_evolve", repeatEvolve);
    }
    if (showOptimiseParams) {
      gd.addChoice("Search", SEARCH, SEARCH[searchParam]);
      gd.addCheckbox("Repeat_search", repeatSearch);
    }
    gd.addStringField("Title", resultsTitle, 20);
    final String[] labels = {"Show_TP", "Show_FP", "Show_FN"};
    gd.addCheckboxGroup(1, 3, labels, new boolean[] {showTP, showFP, showFN});

    gd.showDialog();

    if (gd.wasCanceled() || !readDialog(gd, optimiseParameters, tmp)) {
      return false;
    }

    if (!selectTableColumns()) {
      return false;
    }

    // We may have to read the results again if the ranking option has changed.
    // Also we must read the results with the maximum duplicate distance we may encounter.
    final double dd = duplicateDistance;
    if (showOptimiseParams) {
      duplicateDistance = maxDuplicateDistance;
    }
    readResults();
    duplicateDistance = dd;

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
    String msg = String.format(
        "Fit %d/%d results, %d True-Positives, %d unique\nExpected signal = %.3f +/- %.3f\n"
            + "Expected X precision = %.3f (LSE), %.3f (MLE)\nNot duplicates : %d / %d (%.2f%%)",
        fittedResults, totalResults, matches, maxUniqueId, signal, pSignal, pLse, pMle,
        notDuplicateCount, newResultCount, (100.0 * notDuplicateCount) / newResultCount);

    final FilterResult best = getBestResult();
    if (best != null) {
      msg += String.format("\nCurrent Best=%s, FailCount=%d", MathUtils.rounded(best.score),
          best.failCount);
    }
    gd.addMessage(msg);
  }

  private boolean selectTableColumns() {
    if (showResultsTable || showSummaryTable) {
      final GenericDialog gd = new GenericDialog(TITLE);
      gd.addHelp(About.HELP_URL);

      gd.addMessage("Select the results:");
      for (int i = 0; i < COLUMNS.length; i++) {
        gd.addCheckbox(COLUMNS[i], showColumns[i]);
      }
      gd.showDialog();

      if (gd.wasCanceled()) {
        return false;
      }

      for (int i = 0; i < COLUMNS.length; i++) {
        showColumns[i] = gd.getNextBoolean();
      }

      requireIntegerResults = false;
      for (int i = 0; i < 7; i++) {
        if (showColumns[i]) {
          requireIntegerResults = true;
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

    failCount = (int) Math.abs(gd.getNextNumber());
    if (showOptimiseParams) {
      minFailCount = (int) Math.abs(gd.getNextNumber());
      maxFailCount = (int) Math.abs(gd.getNextNumber());
    }
    if (BenchmarkSpotFit.computeDoublets) {
      // Round to the precision of the min/max
      residualsThreshold = sResidualsThreshold = MathUtils.round(Math.abs(gd.getNextNumber()), 4);
      if (showOptimiseParams) {
        minResidualsThreshold = MathUtils.round(Math.abs(gd.getNextNumber()), 4);
        maxResidualsThreshold = MathUtils.round(Math.abs(gd.getNextNumber()), 4);
      }
    }
    duplicateDistance = MathUtils.round(Math.abs(gd.getNextNumber()), 4);
    if (showOptimiseParams) {
      minDuplicateDistance = MathUtils.round(Math.abs(gd.getNextNumber()), 4);
      maxDuplicateDistance = MathUtils.round(Math.abs(gd.getNextNumber()), 4);
    }
    reset = gd.getNextBoolean();
    showResultsTable = gd.getNextBoolean();
    showSummaryTable = gd.getNextBoolean();
    clearTables = gd.getNextBoolean();
    summaryTopN = (int) Math.abs(gd.getNextNumber());
    summaryDepth = Math.abs(gd.getNextNumber());
    plotTopN = (int) Math.abs(gd.getNextNumber());
    saveBestFilter = gd.getNextBoolean();
    saveTemplate = gd.getNextBoolean();
    calculateSensitivity = gd.getNextBoolean();
    delta = gd.getNextNumber();
    criteriaIndex = gd.getNextChoiceIndex();
    criteriaLimit = gd.getNextNumber();
    scoreIndex = gd.getNextChoiceIndex();
    upperMatchDistance = Math.abs(gd.getNextNumber());
    partialMatchDistance = Math.abs(gd.getNextNumber());
    upperSignalFactor = Math.abs(gd.getNextNumber());
    partialSignalFactor = Math.abs(gd.getNextNumber());
    if (!simulationParameters.fixedDepth) {
      depthRecallAnalysis = gd.getNextBoolean();
    }
    scoreAnalysis = gd.getNextBoolean();
    componentAnalysis = gd.getNextChoiceIndex();
    if (showOptimiseFilter) {
      evolve = gd.getNextChoiceIndex();
      repeatEvolve = gd.getNextBoolean();
    }
    if (showOptimiseParams) {
      searchParam = gd.getNextChoiceIndex();
      repeatSearch = gd.getNextBoolean();
    }
    resultsTitle = gd.getNextString();
    showTP = gd.getNextBoolean();
    showFP = gd.getNextBoolean();
    showFN = gd.getNextBoolean();

    gd.collectOptions();
    duplicateDistanceAbsolute = tmp.getDuplicateDistanceAbsolute();

    resultsPrefix = BenchmarkSpotFit.resultPrefix + "\t" + resultsTitle + "\t";
    createResultsPrefix2();

    // Check there is one output
    if (!showResultsTable && !showSummaryTable && !calculateSensitivity && plotTopN < 1
        && !saveBestFilter) {
      IJ.error(TITLE, "No output selected");
      return false;
    }

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Delta", delta);
      ParameterUtils.isBelow("Delta", delta, 1);
      ParameterUtils.isAboveZero("Upper match distance", upperMatchDistance);
      if (partialMatchDistance > upperMatchDistance) {
        partialMatchDistance = upperMatchDistance;
      }
      if (partialSignalFactor > upperSignalFactor) {
        partialSignalFactor = upperSignalFactor;
      }
      if (showOptimiseParams) {
        ParameterUtils.isEqualOrBelow("Fail count", failCount, maxFailCount);
        ParameterUtils.isEqualOrAbove("Fail count", failCount, minFailCount);
        if (BenchmarkSpotFit.computeDoublets) {
          ParameterUtils.isEqualOrBelow("Residuals threshold", sResidualsThreshold,
              maxResidualsThreshold);
          ParameterUtils.isEqualOrAbove("Residuals threshold", sResidualsThreshold,
              minResidualsThreshold);
        }
        ParameterUtils.isEqualOrBelow("Duplicate distance", duplicateDistance,
            maxDuplicateDistance);
        ParameterUtils.isEqualOrAbove("Duplicate distance", duplicateDistance,
            minDuplicateDistance);
      }
      // Parameters.isEqualOrBelow("Partial match distance", partialMatchDistance,
      // upperMatchDistance);
      // Parameters.isEqualOrBelow("Partial signal factor", partialSignalFactor, upperSignalFactor);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    invertCriteria = requiresInversion(criteriaIndex);
    minCriteria = (invertCriteria) ? -criteriaLimit : criteriaLimit;
    invertScore = requiresInversion(scoreIndex);

    return !gd.invalidNumber();
  }

  private void createResultsPrefix2() {
    createResultsPrefix2(failCount, residualsThreshold, duplicateDistance);
  }

  private void createResultsPrefix2(int failCount, double residualsThreshold,
      double duplicateDistance) {
    resultsPrefix2 = buildResultsPrefix2(failCount, residualsThreshold, duplicateDistance);
    if (!TextUtils.isNullOrEmpty(resultsTitle)) {
      limitFailCount = resultsTitle + ", ";
    } else {
      limitFailCount = "";
    }
    limitFailCount += "f=" + failCount;
    limitFailCount += ", r=" + MathUtils.rounded(residualsThreshold);
  }

  private static String buildResultsPrefix2(int failCount, double residualsThreshold,
      double duplicateDistance) {
    return "\t" + failCount + "\t" + MathUtils.rounded(residualsThreshold) + "\t"
        + MathUtils.rounded(duplicateDistance)
        + ((duplicateDistanceAbsolute) ? " (absolute)" : " (relative)");
  }

  private static boolean showIterationDialog() {
    final GenericDialog gd = new GenericDialog(TITLE);
    final StringBuilder sb = new StringBuilder();
    sb.append("Iterate ").append(BenchmarkSpotFit.TITLE).append(" & ").append(TITLE).append(".\n");
    sb.append(BenchmarkSpotFit.TITLE)
        .append(" will be run once interactively if results cannot be loaded.\n");
    sb.append(TITLE).append(" will be run once interactively to obtain settings.\n \n");
    sb.append("Configure the convergence criteria for iteration:");
    gd.addMessage(sb.toString());
    gd.addNumericField("Score_Tolerance", iterationScoreTolerance, -1);
    gd.addNumericField("Filter_Tolerance", iterationFilterTolerance, -1);
    gd.addCheckbox("Compare_Results", iterationCompareResults);
    gd.addNumericField("Compare_Distance", iterationCompareDistance, 2);
    gd.addNumericField("Iter_Max_Iterations", iterationMaxIterations, 0);
    gd.addMessage("Configure how the parameter range is updated per iteration:");
    gd.addSlider("Min_range_reduction", 0, 1, iterationMinRangeReduction);
    gd.addSlider("Min_range_reduction_iteration", 1, 10, iterationMinRangeReductionIteration);
    gd.addCheckbox("Converge_before_refit", iterationConvergeBeforeRefit);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    iterationScoreTolerance = gd.getNextNumber();
    iterationFilterTolerance = gd.getNextNumber();
    iterationCompareResults = gd.getNextBoolean();
    iterationCompareDistance = Math.abs(gd.getNextNumber());
    iterationMaxIterations = (int) gd.getNextNumber();
    iterationMinRangeReduction = Math.abs(gd.getNextNumber());
    iterationMinRangeReductionIteration = (int) Math.abs(gd.getNextNumber());
    iterationConvergeBeforeRefit = gd.getNextBoolean();

    return true;
  }

  private boolean showReportDialog() {
    final GenericDialog gd = new GenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);

    addSimulationData(gd);

    gd.addCheckbox("Show_table", showResultsTable);
    gd.addCheckbox("Show_summary", showSummaryTable);
    gd.addCheckbox("Clear_tables", clearTables);
    gd.addSlider("Summary_top_n", 0, 20, summaryTopN);
    gd.addCheckbox("Save_best_filter", saveBestFilter);
    gd.addCheckbox("Save_template", saveTemplate);
    gd.addCheckbox("Calculate_sensitivity", calculateSensitivity);
    gd.addSlider("Delta", 0.01, 1, delta);
    if (!simulationParameters.fixedDepth) {
      gd.addCheckbox("Depth_recall_analysis", depthRecallAnalysis);
    }
    gd.addCheckbox("Score_analysis", scoreAnalysis);
    gd.addChoice("Component_analysis", COMPONENT_ANALYSIS, COMPONENT_ANALYSIS[componentAnalysis]);
    gd.addStringField("Title", resultsTitle, 20);
    final String[] labels = {"Show_TP", "Show_FP", "Show_FN"};
    gd.addCheckboxGroup(1, 3, labels, new boolean[] {showTP, showFP, showFN});

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    showResultsTable = gd.getNextBoolean();
    showSummaryTable = gd.getNextBoolean();
    clearTables = gd.getNextBoolean();
    summaryTopN = (int) Math.abs(gd.getNextNumber());
    saveBestFilter = gd.getNextBoolean();
    saveTemplate = gd.getNextBoolean();
    calculateSensitivity = gd.getNextBoolean();
    delta = gd.getNextNumber();
    if (!simulationParameters.fixedDepth) {
      depthRecallAnalysis = gd.getNextBoolean();
    }
    scoreAnalysis = gd.getNextBoolean();
    componentAnalysis = gd.getNextChoiceIndex();
    resultsTitle = gd.getNextString();
    showTP = gd.getNextBoolean();
    showFP = gd.getNextBoolean();
    showFN = gd.getNextBoolean();

    if (gd.invalidNumber()) {
      return false;
    }

    resultsPrefix = BenchmarkSpotFit.resultPrefix + "\t" + resultsTitle + "\t";
    createResultsPrefix2();

    // Check there is one output
    if (!showResultsTable && !showSummaryTable && !calculateSensitivity && !saveBestFilter
        && !saveTemplate) {
      IJ.error(TITLE, "No output selected");
      return false;
    }

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Delta", delta);
      ParameterUtils.isBelow("Delta", delta, 1);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    if (!selectTableColumns()) {
      return false;
    }

    return true;
  }

  private boolean showScoreDialog() {
    final GenericDialog gd = new GenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);

    addSimulationData(gd);

    // Get the last scored filter or default to the best filter
    getScoreFilter();

    gd.addSlider("Fail_count", 0, 20, scoreFailCount);
    if (BenchmarkSpotFit.computeDoublets) {
      gd.addSlider("Residuals_threshold", 0.01, 1, scoreResidualsThreshold);
    }
    gd.addNumericField("Duplicate_distance", scoreDuplicateDistance, 2);

    gd.addTextAreas(uk.ac.sussex.gdsc.core.utils.XmlUtils.convertQuotes(scoreFilter.toXml()), null,
        6, 60);

    gd.addCheckbox("Reset_filter", false);
    // gd.addCheckbox("Show_table", showResultsTable);
    gd.addCheckbox("Show_summary", showSummaryTable);
    gd.addCheckbox("Clear_tables", clearTables);
    // gd.addSlider("Summary_top_n", 0, 20, summaryTopN);
    gd.addCheckbox("Save_best_filter", saveBestFilter);
    gd.addCheckbox("Save_template", saveTemplate);
    gd.addCheckbox("Calculate_sensitivity", calculateSensitivity);
    gd.addSlider("Delta", 0.01, 1, delta);
    if (!simulationParameters.fixedDepth) {
      gd.addCheckbox("Depth_recall_analysis", depthRecallAnalysis);
    }
    gd.addCheckbox("Score_analysis", scoreAnalysis);
    gd.addChoice("Component_analysis", COMPONENT_ANALYSIS, COMPONENT_ANALYSIS[componentAnalysis]);
    gd.addStringField("Title", resultsTitle, 20);
    final String[] labels = {"Show_TP", "Show_FP", "Show_FN"};
    gd.addCheckboxGroup(1, 3, labels, new boolean[] {showTP, showFP, showFN});

    // Dialog to have a reset checkbox. This reverts back to the default.
    if (ImageJUtils.isShowGenericDialog()) {
      final Checkbox cb = (Checkbox) (gd.getCheckboxes().get(0));
      final Vector<TextField> v = gd.getNumericFields();
      final TextArea ta = gd.getTextArea1();
      cb.addItemListener(event -> {
        if (cb.getState()) {
          scoreFilter = null;
          getScoreFilter();
          int index = 0;
          v.get(index++).setText(Integer.toString(scoreFailCount));
          if (BenchmarkSpotFit.computeDoublets) {
            v.get(index++).setText(Double.toString(scoreResidualsThreshold));
          }
          v.get(index++).setText(Double.toString(scoreDuplicateDistance));
          ta.setText(uk.ac.sussex.gdsc.core.utils.XmlUtils.convertQuotes(scoreFilter.toXml()));
        }
      });
    }

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    scoreFailCount = (int) Math.abs(gd.getNextNumber());
    if (BenchmarkSpotFit.computeDoublets) {
      scoreResidualsThreshold = Math.abs(gd.getNextNumber());
    }
    scoreDuplicateDistance = Math.abs(gd.getNextNumber());

    final String xml = gd.getNextText();
    try {
      scoreFilter = (DirectFilter) Filter.fromXml(xml);
    } catch (final Exception ex) {
      scoreFilter = null;
      getScoreFilter();
    }

    final boolean reset = gd.getNextBoolean();
    if (reset) {
      scoreFilter = null;
      getScoreFilter();
    }

    // showResultsTable = gd.getNextBoolean();
    showSummaryTable = gd.getNextBoolean();
    clearTables = gd.getNextBoolean();
    // summaryTopN = (int) Math.abs(gd.getNextNumber());
    saveBestFilter = gd.getNextBoolean();
    saveTemplate = gd.getNextBoolean();
    calculateSensitivity = gd.getNextBoolean();
    delta = gd.getNextNumber();
    if (!simulationParameters.fixedDepth) {
      depthRecallAnalysis = gd.getNextBoolean();
    }
    scoreAnalysis = gd.getNextBoolean();
    componentAnalysis = gd.getNextChoiceIndex();
    resultsTitle = gd.getNextString();
    showTP = gd.getNextBoolean();
    showFP = gd.getNextBoolean();
    showFN = gd.getNextBoolean();

    if (gd.invalidNumber()) {
      return false;
    }

    resultsPrefix = BenchmarkSpotFit.resultPrefix + "\t" + resultsTitle + "\t";
    createResultsPrefix2(scoreFailCount, scoreResidualsThreshold, scoreDuplicateDistance);

    // Check there is one output
    if (!showSummaryTable && !calculateSensitivity && !saveBestFilter && !saveTemplate) {
      IJ.error(TITLE, "No output selected");
      return false;
    }

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Delta", delta);
      ParameterUtils.isBelow("Delta", delta, 1);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    if (!selectTableColumns()) {
      return false;
    }

    return true;
  }

  private static SettingsList lastAnalyseSettings;
  private static SettingsList lastAnalyseParametersSettings;

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
      scores.clear();
      runAnalysis(filterSets, optimum, rangeReduction);

      if (ImageJUtils.isInterrupted()) {
        return null;
      }
    } else {
      // Interactive run, this may be the first run during iterative optimisation
      if (reset) {
        scores.clear();
      }

      createResultsWindow();

      final boolean debugSpeed = false;

      // Only repeat analysis if necessary
      double evolveSetting = evolve;
      if (evolve == 1) {
        // The delta effects the step size for the Genetic Algorithm
        evolveSetting *= delta;
      }
      final SettingsList settings = new SettingsList(filterSets, resultsList, failCount,
          residualsThreshold, duplicateDistance, duplicateDistanceAbsolute, plotTopN, summaryDepth,
          criteriaIndex, criteriaLimit, scoreIndex, evolveSetting);

      final boolean equalSettings = settings.equals(lastAnalyseSettings);

      if (debugSpeed || !equalSettings || (evolve != 0 && repeatEvolve)) {
        newResults = true;
        lastAnalyseSettings = settings;

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
      // if (reset)
      // scores.clear();

      createResultsWindow();

      // Only repeat analysis if necessary
      double min = minResidualsThreshold;
      double max = maxResidualsThreshold;
      if (BenchmarkSpotFit.computeDoublets) {
        min = max = 0;
      }
      final SettingsList settings = new SettingsList(optimum, resultsList, failCount, minFailCount,
          maxFailCount, residualsThreshold, min, max, duplicateDistance, duplicateDistanceAbsolute,
          minDuplicateDistance, maxDuplicateDistance, summaryDepth, criteriaIndex, criteriaLimit,
          scoreIndex, searchParam);

      if (repeatSearch || !settings.equals(lastAnalyseParametersSettings)) {
        newResults = true;
        lastAnalyseParametersSettings = settings;
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
    return reportResults(newResults, new ArrayList<>(bestFilter.values()));
  }

  private ComplexFilterScore reportResults(boolean newResults, ComplexFilterScore optimum) {
    final List<ComplexFilterScore> filters = new ArrayList<>(1);
    filters.add(optimum);
    return reportResults(newResults, filters);
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
    if (showSummaryTable || saveTemplate) {
      final Consumer<String> summaryWindow = createSummaryWindow();
      int count = 0;
      final double range = (summaryDepth / simulationParameters.pixelPitch) * 0.5;
      int np = 0;
      for (final double depth : depthStats) {
        if (Math.abs(depth) < range) {
          np++;
        }
      }
      for (final ComplexFilterScore fs : filters) {
        final ArrayList<FractionalAssignment[]> list = new ArrayList<>(resultsList.length);
        final FractionClassificationResult r =
            scoreFilter(fs.getFilter(), minimalFilter, resultsList, list, coordinateStore);
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
        sb.append(MathUtils.rounded(tp / np)).append('\t');
        sb.append(MathUtils.rounded(distance / scored)).append('\t');
        sb.append(MathUtils.rounded(sf / scored)).append('\t');
        sb.append(MathUtils.rounded(Math.sqrt(rmsd / scored))).append('\t');
        sb.append(MathUtils.rounded(slope)).append('\t');
        if (fs.atLimit() != null) {
          sb.append(fs.atLimit());
        }

        String text = sb.toString();
        if (topFilterSummary == null) {
          topFilterSummary = text;
          if (!showSummaryTable) {
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
        if (summaryTopN > 0 && count >= summaryTopN) {
          break;
        }
      }
      // Add a spacer to the summary table if we have multiple results
      if (count > 1 && showSummaryTable) {
        summaryWindow.accept("");
      }
    }

    final DirectFilter bestFilter = filters.get(0).getFilter();
    if (saveBestFilter) {
      saveFilter(bestFilter);
    }

    if (topFilterClassificationResult == null) {
      topFilterResults = new ArrayList<>(resultsList.length);
      scoreFilter(bestFilter, minimalFilter, resultsList, topFilterResults, coordinateStore);
    }
    if (newResults || scores.isEmpty()) {
      scores.add(new FilterResult(failCount, residualsThreshold, duplicateDistance,
          duplicateDistanceAbsolute, filters.get(0)));
    }

    if (saveTemplate) {
      saveTemplate(topFilterSummary);
    }

    showPlots();
    calculateSensitivity();
    topFilterResults = depthAnalysis(topFilterResults, bestFilter);
    topFilterResults = scoreAnalysis(topFilterResults, bestFilter);
    componentAnalysis(filters.get(0));
    PreprocessedPeakResult[] filterResults = null;
    if (isShowOverlay()) {
      filterResults = showOverlay(topFilterResults, bestFilter);
    }
    saveResults(filterResults, bestFilter);

    wo.tile();

    return filters.get(0);
  }

  private void runAnalysis(List<FilterSet> filterSets) {
    runAnalysis(filterSets, null, 0);
  }

  private void runAnalysis(List<FilterSet> filterSets, ComplexFilterScore optimum,
      double rangeReduction) {
    plots.clear();
    plots.ensureCapacity(plotTopN);
    bestFilter.clear();
    bestFilterOrder.clear();

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
    IJ.showProgress(1);
    IJ.showStatus("");

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
    IJ.showProgress(1);
    IJ.showStatus("");

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

  private static void showPlots() {
    if (plots.isEmpty()) {
      return;
    }

    // Display the top N plots
    final int[] list = new int[plots.size()];
    int index = 0;
    for (final NamedPlot p : plots) {
      final Plot2 plot = new Plot2(p.name, p.xAxisName, COLUMNS[scoreIndex], p.xValues, p.yValues);
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

  private void calculateSensitivity() {
    if (!calculateSensitivity) {
      return;
    }
    if (!bestFilter.isEmpty()) {
      IJ.showStatus("Calculating sensitivity ...");
      final Consumer<String> output = createSensitivityWindow();

      final Ticker ticker =
          Ticker.createStarted(new ImageJTrackProgress(), bestFilter.size(), false);
      for (final String type : bestFilterOrder) {

        final DirectFilter filter = bestFilter.get(type).getFilter();

        FractionClassificationResult score =
            scoreFilter(filter, minimalFilter, resultsList, coordinateStore);
        score = getOriginalScore(score);

        final String message = type + "\t\t\t" + MathUtils.rounded(score.getJaccard(), 4) + "\t\t"
            + MathUtils.rounded(score.getPrecision(), 4) + "\t\t"
            + MathUtils.rounded(score.getRecall(), 4);

        output.accept(message);

        // List all the parameters that can be adjusted.
        final int parameters = filter.getNumberOfParameters();
        for (int index = 0; index < parameters; index++) {
          // For each parameter compute as upward + downward delta and get the average gradient
          final DirectFilter higher = (DirectFilter) filter.adjustParameter(index, delta);
          final DirectFilter lower = (DirectFilter) filter.adjustParameter(index, -delta);

          FractionClassificationResult scoreHigher =
              scoreFilter(higher, minimalFilter, resultsList, coordinateStore);
          scoreHigher = getOriginalScore(scoreHigher);
          FractionClassificationResult scoreLower =
              scoreFilter(lower, minimalFilter, resultsList, coordinateStore);
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

      IJ.showProgress(1);
      IJ.showStatus("");
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
    if (!showResultsTable) {
      return;
    }

    if (isHeadless) {
      IJ.log(createResultsHeader(false));
    } else {
      resultsWindow = ImageJUtils.refresh(resultsWindowRef, () -> {
        final String header = createResultsHeader(false);
        return new TextWindow(TITLE + " Results", header, "", 900, 300);
      });
      if (clearTables) {
        resultsWindow.getTextPanel().clear();
      }
    }
  }

  private Consumer<String> createSummaryWindow() {
    if (!showSummaryTable) {
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
    if (clearTables) {
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
      if (clearTables) {
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
    if (clearTables) {
      window.getTextPanel().clear();
    }
    return window::append;
  }

  private static String createComponentAnalysisHeader() {
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

  private static String createResultsHeader(boolean summary) {
    final StringBuilder sb = new StringBuilder(BenchmarkSpotFit.tablePrefix);
    sb.append(
        "\tTitle\tName\tFail\tRes\tDup D\tLower D (nm)\tUpper D (nm)\tLower factor\tUpper factor");

    for (int i = 0; i < COLUMNS.length; i++) {
      if (showColumns[i]) {
        sb.append('\t').append(COLUMNS[i]);
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

    this.gaResultsList = resultsList;

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
        // final Filter incF = filterSet.getFilters().get(2);

        for (int i = 0; i < n; i++) {
          // Do not disable if the increment is not set. This is left to the user to decide.
          // if (incF.getParameterValue(i) == incF.getDisabledParameterValue(i) ||
          // Double.isInfinite(incF.getParameterValue(i)))
          // {
          // // Not enabled
          // dimensions[i] = new SearchDimension(incF.getDisabledParameterValue(i));
          // continue;
          // }

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
        for (int i = 0; i < n; i++) {
          if (lower[i] == upper[i]) {
            // Not enabled
            originalDimensions[i] = new FixedDimension(lower[i]);
            continue;
          }
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

        if (originalDimensions == null) {
          // Failed to work out the dimensions. No optimisation will be possible.

          // Sort so that the filters are in a nice order for reporting
          filterSet.sort();

          // This will not be used when the dimensions are null
          seed = null;
        }
      }

      if (originalDimensions != null) {
        // Use the current optimum if we are doing a range optimisation
        if (currentOptimum != null && rangeInput
            && currentOptimum.getType().equals(searchScoreFilter.getType()) && evolve != 0) {
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

    if (evolve == 1 && originalDimensions != null) {
      // Collect parameters for the genetic algorithm
      pauseFilterTimer();

      // Remember the step size settings
      double[] stepSize = stepSizeMap.get(setNumber);
      if (stepSize == null || stepSize.length != searchScoreFilter.length()) {
        stepSize = searchScoreFilter.mutationStepRange().clone();
        for (int j = 0; j < stepSize.length; j++) {
          stepSize[j] *= delta;
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
        gd.addNumericField(prefix + "Population_size", populationSize, 0);
        gd.addNumericField(prefix + "Failure_limit", failureLimit, 0);
        gd.addNumericField(prefix + "Tolerance", tolerance, -1);
        gd.addNumericField(prefix + "Converged_count", convergedCount, 0);
        gd.addSlider(prefix + "Mutation_rate", 0.05, 1, mutationRate);
        gd.addSlider(prefix + "Crossover_rate", 0.05, 1, crossoverRate);
        gd.addSlider(prefix + "Mean_children", 0.05, 3, meanChildren);
        gd.addSlider(prefix + "Selection_fraction", 0.05, 0.5, selectionFraction);
        gd.addCheckbox(prefix + "Ramped_selection", rampedSelection);
        gd.addCheckbox(prefix + "Save_option", saveOption);

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
            throw new RuntimeException("The dialog has no been shown");
          }
          populationSize = (int) Math.abs(gd.getNextNumber());
          if (populationSize < 10) {
            populationSize = 10;
          }
          failureLimit = (int) Math.abs(gd.getNextNumber());
          tolerance = gd.getNextNumber();
          convergedCount = (int) gd.getNextNumber(); // Allow negatives
          mutationRate = Math.abs(gd.getNextNumber());
          crossoverRate = Math.abs(gd.getNextNumber());
          meanChildren = Math.abs(gd.getNextNumber());
          selectionFraction = Math.abs(gd.getNextNumber());
          rampedSelection = gd.getNextBoolean();
          saveOption = gd.getNextBoolean();

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

        // // Reset negatives to zero
        // stepSize = stepSize.clone();
        // for (int j = 0; j < stepSize.length; j++)
        // if (stepSize[j] < 0)
        // stepSize[j] = 0;

        // Create the genetic algorithm
        final RandomDataGenerator random = new RandomDataGenerator(new Well44497b());
        final SimpleMutator<FilterScore> mutator = new SimpleMutator<>(random, mutationRate);
        // Override the settings with the step length, a min of zero and the configured upper
        final double[] upper = searchScoreFilter.upperLimit();
        mutator.overrideChromosomeSettings(stepSize, new double[stepSize.length], upper);
        final Recombiner<FilterScore> recombiner =
            new SimpleRecombiner<>(random, crossoverRate, meanChildren);
        SelectionStrategy<FilterScore> selectionStrategy;
        // If the initial population is huge ensure that the first selection culls to the correct
        // size
        final int selectionMax = (int) (selectionFraction * populationSize);
        if (rampedSelection) {
          selectionStrategy =
              new RampedSelectionStrategy<>(random, selectionFraction, selectionMax);
        } else {
          selectionStrategy =
              new SimpleSelectionStrategy<>(random, selectionFraction, selectionMax);
        }
        final ToleranceChecker<FilterScore> ga_checker =
            new InterruptChecker(tolerance, tolerance * 1e-3, convergedCount);

        // Create new random filters if the population is initially below the population size
        List<Filter> filters = filterSet.getFilters();
        if (filterSet.size() < populationSize) {
          filters = new ArrayList<>(populationSize);
          // Add the existing filters if they are not a range input file
          if (!rangeInput) {
            filters.addAll(filterSet.getFilters());
          }
          // Add current optimum to seed
          if (currentOptimum != null && nonInteractive) {
            filters.add(currentOptimum);
          }
          // The GA does not use the min interval grid so sample without rounding
          final double[][] sample =
              SearchSpace.sampleWithoutRounding(dimensions, populationSize - filters.size(), null);
          filters.addAll(searchSpaceToFilters(sample));
        }
        gaPopulation = new Population<>(filters);
        gaPopulation.setPopulationSize(populationSize);
        gaPopulation.setFailureLimit(failureLimit);
        selectionStrategy.setTracker(this);

        // Evolve
        algorithm = EVOLVE[evolve];
        gaStatusPrefix = algorithm + " [" + setNumber + "] " + filterSet.getName() + " ... ";
        gaIteration = 0;
        gaPopulation.setTracker(this);

        createGaWindow();
        resumeFilterTimer();

        best = gaPopulation.evolve(mutator, recombiner, this, selectionStrategy, ga_checker);

        if (best != null) {
          // In case optimisation was stopped
          IJ.resetEscape();

          // The GA may produce coordinates off the min interval grid
          best = enumerateMinInterval(best, stepSize, indices);

          // Now update the filter set for final assessment
          filterSet = new FilterSet(filterSet.getName(),
              populationToFilters(gaPopulation.getIndividuals()));

          // Option to save the filters
          if (saveOption) {
            saveFilterSet(filterSet, setNumber, !nonInteractive);
          }
        }
      } else {
        resumeFilterTimer();
      }
    }

    if ((evolve == 2 || evolve == 4) && originalDimensions != null) {
      // Collect parameters for the range search algorithm
      pauseFilterTimer();

      final boolean isStepSearch = evolve == 4;

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
        gd.addMessage("Configure the " + EVOLVE[evolve] + " algorithm for [" + setNumber + "] "
            + filterSet.getName());
        gd.addSlider(prefix + "Width", 1, 5, rangeSearchWidth);
        if (!isStepSearch) {
          gd.addCheckbox(prefix + "Save_option", saveOption);
          gd.addNumericField(prefix + "Max_iterations", maxIterations, 0);
          final String[] modes =
              SettingsManager.getNames((Object[]) SearchSpace.RefinementMode.values());
          gd.addSlider(prefix + "Reduce", 0.01, 0.99, rangeSearchReduce);
          gd.addChoice("Refinement", modes, modes[refinementMode]);
        }
        gd.addNumericField(prefix + "Seed_size", seedSize, 0);

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
          rangeSearchWidth = (int) gd.getNextNumber();
          if (!isStepSearch) {
            saveOption = gd.getNextBoolean();
            maxIterations = (int) gd.getNextNumber();
            refinementMode = gd.getNextChoiceIndex();
            rangeSearchReduce = gd.getNextNumber();
          }
          seedSize = (int) gd.getNextNumber();
          for (int i = 0; i < n; i++) {
            enabled[i] = gd.getNextBoolean();
          }

          searchRangeMap.put(setNumber, enabled);
        }

        if (!isStepSearch) {
          myRefinementMode = SearchSpace.RefinementMode.values()[refinementMode];
        }
        for (int i = 0; i < n; i++) {
          if (enabled[i]) {
            try {
              dimensions[i] = originalDimensions[i].create(rangeSearchWidth);
              dimensions[i].setPad(true);
              // Prevent range reduction so that the step search just does a single refinement step
              dimensions[i].setReduceFactor((isStepSearch) ? 1 : rangeSearchReduce);
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
          gd.addMessage(String.format(
              "%d combinations for the configured dimensions.\n \nClick 'Yes' to optimise.",
              combinations));
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
          algorithm = EVOLVE[evolve] + " " + rangeSearchWidth;
          gaStatusPrefix = algorithm + " [" + setNumber + "] " + filterSet.getName() + " ... ";
          gaIteration = 0;
          filterScoreOptimum = null;

          final SearchSpace ss = new SearchSpace();
          ss.setTracker(this);
          if (seedSize > 0) {
            double[][] sample;
            // Add current optimum to seed
            if (nonInteractive && currentOptimum != null) {
              sample = new double[1][];
              sample[0] = currentOptimum.getParameters();
              seed = merge(seed, sample);
            }
            final int size = (seed == null) ? 0 : seed.length;
            // Sample without rounding as the seed will be rounded
            sample = SearchSpace.sampleWithoutRounding(dimensions, seedSize - size, null);
            seed = merge(seed, sample);
          }
          // Note: If we have an optimum and we are not seeding this should not matter as the
          // dimensions
          // have been centred on the current optimum
          ss.seed(seed);
          final ConvergenceChecker<FilterScore> checker =
              new InterruptConvergenceChecker(0, 0, maxIterations);

          createGaWindow();
          resumeFilterTimer();

          final SearchResult<FilterScore> optimum =
              ss.search(dimensions, this, checker, myRefinementMode);

          if (optimum != null) {
            // In case optimisation was stopped
            IJ.resetEscape();

            best = ((SimpleFilterScore) optimum.getScore()).result.filter;

            if (seedSize > 0) {
              // Not required as the search now respects the min interval
              // The optimum may be off grid if it was from the seed
              // best = enumerateMinInterval(best, enabled);
            }

            // Now update the filter set for final assessment
            filterSet = new FilterSet(filterSet.getName(),
                searchSpaceToFilters((DirectFilter) best, ss.getSearchSpace()));

            // Option to save the filters
            if (saveOption) {
              saveFilterSet(filterSet, setNumber, !nonInteractive);
            }
          }
        }
      } else {
        resumeFilterTimer();
      }
    }

    if (evolve == 3 && originalDimensions != null) {
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
        gd.addCheckbox(prefix + "Save_option", saveOption);
        gd.addNumericField(prefix + "Max_iterations", maxIterations, 0);
        gd.addNumericField(prefix + "Converged_count", convergedCount, 0);
        gd.addNumericField(prefix + "Samples", enrichmentSamples, 0);
        gd.addSlider(prefix + "Fraction", 0.01, 0.99, enrichmentFraction);
        gd.addSlider(prefix + "Padding", 0, 0.99, enrichmentPadding);

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
          saveOption = gd.getNextBoolean();
          maxIterations = (int) gd.getNextNumber();
          convergedCount = (int) gd.getNextNumber();
          enrichmentSamples = (int) gd.getNextNumber();
          enrichmentFraction = gd.getNextNumber();
          enrichmentPadding = gd.getNextNumber();
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

        algorithm = EVOLVE[evolve];
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
            new InterruptConvergenceChecker(0, 0, maxIterations, convergedCount);

        createGaWindow();
        resumeFilterTimer();

        final SearchResult<FilterScore> optimum = ss.enrichmentSearch(dimensions, this, checker,
            enrichmentSamples, enrichmentFraction, enrichmentPadding);

        if (optimum != null) {
          // In case optimisation was stopped
          IJ.resetEscape();

          best = ((SimpleFilterScore) optimum.getScore()).result.filter;

          // Not required as the search now respects the min interval
          // Enumerate on the min interval to produce the final filter
          // best = enumerateMinInterval(best, enabled);

          // Now update the filter set for final assessment
          filterSet = new FilterSet(filterSet.getName(),
              searchSpaceToFilters((DirectFilter) best, ss.getSearchSpace()));

          // Option to save the filters
          if (saveOption) {
            saveFilterSet(filterSet, setNumber, !nonInteractive);
          }
        }
      } else {
        resumeFilterTimer();
      }
    }

    IJ.showStatus("Analysing [" + setNumber + "] " + filterSet.getName() + " ...");

    // Do not support plotting if we used optimisation
    final double[] xValues =
        (best != null || isHeadless || (plotTopN == 0)) ? null : new double[filterSet.size()];
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
        scoreFilters(setUncomputedStrength(filterSet), showResultsTable);
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

    if (showResultsTable) {
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
              limitFailCount + limitRange, sb.toString());
        } else {
          ImageJUtils.log(
              "Warning: Top filter (%s @ -|%s) [%s] at the limit of the expanded range%s",
              filter.getName(), MathUtils.rounded((invertCriteria) ? -max.criteria : max.criteria),
              limitFailCount + limitRange, sb.toString());
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

    final boolean allowDuplicates = true; // This could be an option?

    // XXX - Commented out the requirement to be the same type to store for later analysis.
    // This may break the code, however I think that all filter sets should be able to have a best
    // filter
    // irrespective of whether they were the same type or not.
    // if (allSameType)
    // {
    final ComplexFilterScore newFilterScore =
        new ComplexFilterScore(max.result, atLimit, algorithm, analysisStopWatch.getTime(), "", 0);
    addBestFilter(type, allowDuplicates, newFilterScore);
    // }

    // Add spacer at end of each result set
    if (isHeadless) {
      if (showResultsTable && filterSet.size() > 1) {
        IJ.log("");
      }
    } else {
      if (showResultsTable && filterSet.size() > 1) {
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
            xValues[ii] = ii + 1;
          }
        }

        final String title = filterSet.getName();

        // Check if a previous filter set had the same name, update if necessary
        NamedPlot plot = getNamedPlot(title);
        if (plot == null) {
          plots.add(new NamedPlot(title, xAxisName, xValues, yValues));
        } else {
          plot.updateValues(xAxisName, xValues, yValues);
        }

        if (plots.size() > plotTopN) {
          Collections.sort(plots);
          plot = plots.remove(plots.size() - 1);
        }
      }
    }

    return 0;
  }

  private static void addBestFilter(String type, boolean allowDuplicates,
      ComplexFilterScore newFilterScore) {
    final ComplexFilterScore filterScore = bestFilter.get(type);
    if (filterScore != null) {
      if (allowDuplicates) {
        // Duplicate type: create a unique key
        // Start at 2 to show it is the second one of the same type
        int count = 2;
        while (bestFilter.containsKey(type + count)) {
          count++;
        }
        type += count;
        bestFilter.put(type, newFilterScore);
        bestFilterOrder.add(type);

        // Replace (even if the same so that the latest results settings are stored)
      } else if (newFilterScore.compareTo(filterScore) <= 0) {
        bestFilter.put(type, newFilterScore);
        // filterScore.update(max.r, atLimit, algorithm, filterSetStopWatch.getTime());
      }
    } else {
      bestFilter.put(type, newFilterScore);
      bestFilterOrder.add(type);
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
      seed = merged.toArray(new double[merged.size()][]);
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
          IJ.error(TITLE, String.format("Unable to configure dimension [%d] %s: " + ex.getMessage(),
              j, searchScoreFilter.getParameterName(j)));
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
    this.gaResultsList = resultsList;
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
      originalDimensions[++index] = new FixedDimension(minFailCount, maxFailCount, 1);
      // TODO - let the min intervals be configured, maybe via extra options
      if (BenchmarkSpotFit.computeDoublets) {
        originalDimensions[++index] =
            new FixedDimension(minResidualsThreshold, maxResidualsThreshold, 0.05);
      } else {
        originalDimensions[++index] = new FixedDimension(1, 1, 0.05);
      }
      originalDimensions[++index] =
          new FixedDimension(minDuplicateDistance, maxDuplicateDistance, 0.5);
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

    if (searchParam == 0 || searchParam == 2) {
      // Collect parameters for the range search algorithm
      pauseParameterTimer();

      final boolean isStepSearch = searchParam == 2;

      // The step search should use a multi-dimension refinement and no range reduction
      SearchSpace.RefinementMode myRefinementMode = SearchSpace.RefinementMode.MULTI_DIMENSION;

      GenericDialog gd = null;
      boolean runAlgorithm = nonInteractive;
      if (!nonInteractive) {
        // Ask the user for the search parameters.
        gd = new GenericDialog(TITLE);
        gd.addMessage("Configure the " + SEARCH[searchParam] + " algorithm for "
            + searchScoreFilter.getType());
        gd.addSlider("Width", 1, 5, pRangeSearchWidth);
        if (!isStepSearch) {
          gd.addNumericField("Max_iterations", pMaxIterations, 0);
          final String[] modes =
              SettingsManager.getNames((Object[]) SearchSpace.RefinementMode.values());
          gd.addSlider("Reduce", 0.01, 0.99, pRangeSearchReduce);
          gd.addChoice("Refinement", modes, modes[pRefinementMode]);
        }
        gd.addNumericField("Seed_size", pSeedSize, 0);

        gd.showDialog();
        runAlgorithm = !gd.wasCanceled();
      }

      if (runAlgorithm) {
        final SearchDimension[] dimensions = new SearchDimension[originalDimensions.length];
        if (!nonInteractive && gd != null) {
          pRangeSearchWidth = (int) gd.getNextNumber();
          if (!isStepSearch) {
            pMaxIterations = (int) gd.getNextNumber();
            pRangeSearchReduce = gd.getNextNumber();
            pRefinementMode = gd.getNextChoiceIndex();
          }
          pSeedSize = (int) gd.getNextNumber();
        }

        if (!isStepSearch) {
          myRefinementMode = SearchSpace.RefinementMode.values()[pRefinementMode];
        }
        for (int i = 0; i < dimensions.length; i++) {
          if (originalDimensions[i].isActive()) {
            try {
              dimensions[i] = originalDimensions[i].create(pRangeSearchWidth);
              dimensions[i].setPad(true);
              // Prevent range reduction so that the step search just does a single refinement step
              dimensions[i].setReduceFactor((isStepSearch) ? 1 : pRangeSearchReduce);
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
        if (!nonInteractive && combinations > 10000) {
          gd = new GenericDialog(TITLE);
          gd.addMessage(String.format(
              "%d combinations for the configured dimensions.\n \nClick 'Yes' to optimise.",
              combinations));
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
          algorithm = SEARCH[searchParam] + " " + pRangeSearchWidth;
          gaStatusPrefix = algorithm + " " + searchScoreFilter.getName() + " ... ";
          gaIteration = 0;
          parameterScoreOptimum = null;

          final SearchSpace ss = new SearchSpace();
          ss.setTracker(this);
          if (pSeedSize > 0) {
            // Add current optimum to seed
            // Note: If we have an optimum and we are not seeding this should not matter as the
            // dimensions
            // have been centred on the current optimum
            final double[][] seed = new double[1][];
            seed[0] = point;
            // Sample without rounding as the seed will be rounded
            final double[][] sample =
                SearchSpace.sampleWithoutRounding(dimensions, pSeedSize - 1, null);
            ss.seed(merge(sample, seed));
          }
          final ConvergenceChecker<FilterScore> checker =
              new InterruptConvergenceChecker(0, 0, pMaxIterations);

          createGaWindow();
          resumeParameterTimer();

          optimum = ss.search(dimensions, new ParameterScoreFunction(), checker, myRefinementMode);

          if (optimum != null) {
            // In case optimisation was stopped
            IJ.resetEscape();

            // Now update the parameters for final assessment
            point = optimum.getPoint();

            // Not required as the seed in now rounded
            // if (pSeedSize > 0)
            // {
            // // The optimum may be off grid if it was from the seed
            // point = enumerateMinInterval(point, names, originalDimensions);
            // }
          }
        }
      } else {
        resumeParameterTimer();
      }
    }

    if (searchParam == 1) {
      // Collect parameters for the enrichment search algorithm
      pauseParameterTimer();

      GenericDialog gd = null;
      boolean runAlgorithm = nonInteractive;
      if (!nonInteractive) {
        // Ask the user for the search parameters.
        gd = new GenericDialog(TITLE);
        gd.addMessage("Configure the " + SEARCH[searchParam] + " algorithm for "
            + searchScoreFilter.getType());
        gd.addNumericField("Max_iterations", pMaxIterations, 0);
        gd.addNumericField("Converged_count", pConvergedCount, 0);
        gd.addNumericField("Samples", pEnrichmentSamples, 0);
        gd.addSlider("Fraction", 0.01, 0.99, pEnrichmentFraction);
        gd.addSlider("Padding", 0, 0.99, pEnrichmentPadding);

        gd.showDialog();
        runAlgorithm = !gd.wasCanceled();
      }

      if (runAlgorithm) {
        final FixedDimension[] dimensions =
            Arrays.copyOf(originalDimensions, originalDimensions.length);
        if (!nonInteractive && gd != null) {
          pMaxIterations = (int) gd.getNextNumber();
          pConvergedCount = (int) gd.getNextNumber();
          pEnrichmentSamples = (int) gd.getNextNumber();
          pEnrichmentFraction = gd.getNextNumber();
          pEnrichmentPadding = gd.getNextNumber();
        }

        algorithm = SEARCH[searchParam];
        gaStatusPrefix = algorithm + " " + searchScoreFilter.getName() + " ... ";
        gaIteration = 0;
        parameterScoreOptimum = null;

        final SearchSpace ss = new SearchSpace();
        ss.setTracker(this);
        // Add current optimum to seed
        final double[][] seed = new double[1][];
        seed[0] = point;
        ss.seed(seed);
        final ConvergenceChecker<FilterScore> checker =
            new InterruptConvergenceChecker(0, 0, pMaxIterations, pConvergedCount);

        createGaWindow();
        resumeParameterTimer();

        optimum = ss.enrichmentSearch(dimensions, new ParameterScoreFunction(), checker,
            pEnrichmentSamples, pEnrichmentFraction, pEnrichmentPadding);

        if (optimum != null) {
          // In case optimisation was stopped
          IJ.resetEscape();

          point = optimum.getPoint();

          // Not required as the search now respects the min interval
          // Enumerate on the min interval to produce the final filter
          // point = enumerateMinInterval(point, names, originalDimensions);
        }
      } else {
        resumeParameterTimer();
      }
    }

    if (searchParam == 3) {
      // Collect parameters for the enumeration search algorithm
      pauseParameterTimer();

      final SearchDimension[] dimensions = new SearchDimension[originalDimensions.length];
      for (int i = 0; i < dimensions.length; i++) {
        if (originalDimensions[i].isActive()) {
          try {
            dimensions[i] = originalDimensions[i].create(0);
          } catch (final IllegalArgumentException ex) {
            IJ.error(TITLE, String
                .format("Unable to configure dimension [%d] %s: " + ex.getMessage(), i, names[i]));
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
        gd.addMessage(String.format(
            "%d combinations for the configured dimensions.\n \nClick 'Yes' to optimise.",
            combinations));
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
        algorithm = SEARCH[searchParam];
        gaStatusPrefix = algorithm + " " + searchScoreFilter.getName() + " ... ";
        gaIteration = 0;
        parameterScoreOptimum = null;

        final SearchSpace ss = new SearchSpace();

        ss.setTracker(this);
        createGaWindow();
        resumeParameterTimer();

        optimum = ss.findOptimum(dimensions, new ParameterScoreFunction());

        if (optimum != null) {
          // In case optimisation was stopped
          IJ.resetEscape();

          // Now update the parameters for final assessment
          point = optimum.getPoint();
        }
      }
    }

    IJ.showStatus("Analysing " + searchScoreFilter.getName() + " ...");

    // Update the parameters using the optimum
    failCount = (int) Math.round(point[0]);
    residualsThreshold = sResidualsThreshold = point[1];
    duplicateDistance = point[2];
    // Refresh the coordinate store
    if (coordinateStore == null
        // Due to the scaling factor the distance may not be exactly the same
        || DoubleEquality.relativeError(duplicateDistance,
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
          minimalFilter, failCount, residualsThreshold, duplicateDistance,
          createCoordinateStore(duplicateDistance), false);

      ImageJUtils.log("Weird re-score of the filter: %f!=%f or %f!=%f (%f:%f)", scoreResult.score,
          optimum.getScore().score, scoreResult.criteria, optimum.getScore().criteria, r.score,
          r.criteria);
    }

    final SimpleFilterScore max =
        new SimpleFilterScore(scoreResult, true, scoreResult.criteria >= minCriteria);

    analysisStopWatch.stop();

    if (showResultsTable) {
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
            limitFailCount + limitRange, sb.toString());
      } else {
        ImageJUtils.log("Warning: Top filter (%s @ -|%s) [%s] at the limit of the expanded range%s",
            searchScoreFilter.getName(),
            MathUtils.rounded((invertCriteria) ? -max.criteria : max.criteria),
            limitFailCount + limitRange, sb.toString());
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
    if (showResultsTable) {
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

  private static StopWatch filterAnalysisStopWatch;
  private static StopWatch parameterAnalysisStopWatch;
  private StopWatch analysisStopWatch;
  private StopWatch iterationStopWatch;

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
    return getScore(result, criteriaIndex, invertCriteria);
  }

  private double getScore(FractionClassificationResult result) {
    return getScore(result, scoreIndex, invertScore);
  }

  private static double getScore(FractionClassificationResult result, final int index,
      final boolean invert) {
    final double score = getScore(result, index);
    return (invert) ? -score : score;
  }

  private static double getScore(FractionClassificationResult result, final int index) {
    // This order must match the COLUMNS order
    switch (index) {
      case 0:
        return result.getNumberOfPositives();
      case 1:
        return result.getNumberOfNegatives();
      case 2:
        return (double) nActual - result.getNumberOfPositives();
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

  private static NamedPlot getNamedPlot(String title) {
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

    // multiPathFilter.setDebugFile("/tmp/fractionScoreSubset.txt");

    // Note: We always use the subset method since fail counts have been accumulated when we read in
    // the results.
    return multiPathFilter.fractionScoreSubset(resultsList,
        ConsecutiveFailCounter.create(failCount), nActual, allAssignments, scoreStore,
        coordinateStore);
  }

  private FilterScoreResult scoreFilter(DirectFilter filter, DirectFilter minFilter,
      boolean createTextResult, CoordinateStore coordinateStore) {
    final FractionClassificationResult r =
        scoreFilter(filter, minFilter, gaResultsListToScore, coordinateStore);

    // // DEBUG - Test if the two methods produce the same results
    // FractionClassificationResult r2 = scoreFilter(filter, minFilter,
    // BenchmarkFilterAnalysis.clonedResultsList);
    // if
    // (!uk.ac.sussex.gdsc.core.utils.DoubleEquality.almostEqualRelativeOrAbsolute(
    // r.getTruePositives(), r2.getTruePositives(), 1e-6, 1e-10) ||
    // !uk.ac.sussex.gdsc.core.utils.DoubleEquality.almostEqualRelativeOrAbsolute(
    // r.getFalsePositives(), r2.getFalsePositives(), 1e-6, 1e-10) ||
    // !uk.ac.sussex.gdsc.core.utils.DoubleEquality.almostEqualRelativeOrAbsolute(r.getFN(),
    // r2.getFN(), 1e-6, 1e-10))
    // {
    // System.out.printf("TP %f != %f, FP %f != %f, FN %f != %f : %s\n", r.getTruePositives(),
    // r2.getTruePositives(), r.getFalsePositives(),
    // r2.getFalsePositives(), r.getFN(), r2.getFN(), filter.getName());
    //
    // // // Debug
    // // MultiPathFilter multiPathFilter = createMPF(filter, minFilter);
    // // multiPathFilter.setDebugFile("/tmp/1.txt");
    // // multiPathFilter.fractionScoreSubset(ga_resultsListToScore, failCount, nActual, null);
    // // multiPathFilter = createMPF(filter, minFilter);
    // // multiPathFilter.setDebugFile("/tmp/2.txt");
    // // multiPathFilter.fractionScoreSubset(BenchmarkFilterAnalysis.clonedResultsList, failCount,
    // // nActual, null);
    // }
    // else
    // {
    // //System.out.println("Matched scores");
    // }

    final double score = getScore(r);
    final double criteria = getCriteria(r);

    // Show the result if it achieves the criteria limit
    final String text =
        (createTextResult && criteria >= minCriteria) ? createResult(filter, r).toString() : null;

    return new FilterScoreResult(score, criteria, filter, text);
  }

  private FilterScoreResult scoreFilter(DirectFilter filter) {
    final FractionClassificationResult r =
        scoreFilter(filter, minimalFilter, resultsList, coordinateStore);
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
        ConsecutiveFailCounter.create(failCount), nActual, null, null, coordinateStore);

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
    final MultiPathFilter multiPathFilter = createMpf(filter, minimalFilter);

    final ArrayList<FractionalAssignment[]> allAssignments = new ArrayList<>(resultsList.length);
    multiPathFilter.fractionScoreSubset(resultsList, ConsecutiveFailCounter.create(failCount),
        nActual, allAssignments, null, coordinateStore);
    return allAssignments;
  }

  private StringBuilder createResult(DirectFilter filter, FractionClassificationResult result) {
    return createResult(filter, result, resultsPrefix2);
  }

  private StringBuilder createResult(DirectFilter filter, FractionClassificationResult result,
      String resultsPrefix2) {
    final StringBuilder sb = new StringBuilder(resultsPrefix);
    sb.append(filter.getName()).append(resultsPrefix2).append(resultsPrefix3);

    int index = 0;

    // TODO - Fix the scores that we show since we no longer have TN results
    // We could set:
    // TN as any candidate that does not match a true result.
    // FN as any candidate that does match a true result (that is not matched by any fit result)
    // To do this properly would require that we store all the matches of candidates to the data.
    // These can then be totalled up given the candidates that have not been used to create a
    // positive.

    // Integer results
    if (requireIntegerResults) {
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

  private static ClassificationResult createIntegerResult(FractionClassificationResult result) {
    return new ClassificationResult(result.getNumberOfPositives(), result.getNumberOfNegatives(), 0,
        nActual - result.getNumberOfPositives());
  }

  /**
   * Gets the original score.
   *
   * @param result the result
   * @return the original score
   */
  private static FractionClassificationResult
      getOriginalScore(FractionClassificationResult result) {
    throw new IllegalStateException("fix this");

    // // Score the fitting results against the original simulated data:
    // // TP are all fit results that can be matched to a spot
    // // FP are all fit results that cannot be matched to a spot
    // // FN are the number of missed spots
    // // Note: We cannot calculate TN since this is the number of fit candidates that are
    // // filtered after fitting that do not match a spot or were not fitted.
    // final double fp = result.getPositives() - result.getTruePositives();
    // final double fn = nActual - result.getTruePositives();
    // return new FractionClassificationResult(result.getTruePositives(), fp, 0, fn);
  }

  private static void add(StringBuilder sb, int value) {
    sb.append('\t').append(value);
  }

  private static void add(StringBuilder sb, String value) {
    sb.append('\t').append(value);
  }

  private static void add(StringBuilder sb, int value, int index) {
    if (showColumns[index]) {
      add(sb, value);
    }
  }

  private static void add(StringBuilder sb, double value, int index) {
    if (showColumns[index]) {
      add(sb, MathUtils.rounded(value));
    }
  }

  private static void addCount(StringBuilder sb, double value, int index) {
    if (showColumns[index]) {
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

  private static void saveFilter(DirectFilter filter) {
    // Save the filter to file
    final String filename = getFilename("Best_Filter_File", filterFilename);
    if (filename != null) {
      filterFilename = filename;
      Prefs.set(KEY_FILTER_FILENAME, filename);

      final List<Filter> filters = new ArrayList<>(1);
      filters.add(filter);
      final FilterSet filterSet = new FilterSet(filter.getName(), filters);
      final List<FilterSet> list = new ArrayList<>(1);
      list.add(filterSet);
      saveFilterSet(filterSet, filename);
    }
  }

  /**
   * Gets the filename (with a .xml extension).
   *
   * @param title the title
   * @param filename the filename
   * @return the filename
   */
  static String getFilename(String title, String filename) {
    filename = ImageJUtils.getFilename(title, filename);
    // Use XML extension
    if (filename != null) {
      filename = FileUtils.replaceExtension(filename, ".xml");
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
      final String filename = getFilename("Filter_set_" + setNumber, filterSetFilename);
      if (filename != null) {
        filterSetFilename = filename;
        Prefs.set(KEY_FILTERSET_FILENAME, filename);
        saveFilterSet(filterSet, filename);
      }
    } else {
      // Add support for multiple filter sets
      saveFilterSet(filterSet, filterSetFilename);
    }

    resumeFilterTimer();
  }

  private static ResultsImageSampler sampler;

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

    String filename = getFilename("Template_File", templateFilename);
    if (filename != null) {
      templateFilename = filename;
      Prefs.set(KEY_TEMPLATE_FILENAME, filename);
      final TemplateSettings.Builder settings = TemplateSettings.newBuilder();
      getNotes(settings, topFilterSummary);
      settings.setFitEngineSettings(config.getFitEngineSettings());
      if (!SettingsManager.toJson(settings.build(), filename, SettingsManager.FLAG_SILENT)) {
        IJ.log("Unable to save the template configuration");
        return;
      }

      // Save some random frames from the test image data
      final ImagePlus imp = CreateData.getImage();
      if (imp == null) {
        return;
      }

      // Get the number of frames
      final ImageStack stack = imp.getImageStack();

      if (sampler == null || sampler.getResults() != results) {
        sampler = new ResultsImageSampler(results, stack, 32);
        sampler.analyse();
      }
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
        nNo = Integer.parseInt(Macro.getValue(options, keyNo, Integer.toString(nNo)));
        nLow = Integer.parseInt(Macro.getValue(options, keyLow, Integer.toString(nLow)));
        nHigh = Integer.parseInt(Macro.getValue(options, keyHigh, Integer.toString(nHigh)));
      } else if (nLow + nHigh == 0) {
        nLow = nHigh = 1;
      }

      final ImagePlus[] out = new ImagePlus[1];
      out[0] = sampler.getSample(nNo, nLow, nHigh);

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

        // For the dialog
        final String msg = String.format(
            "Showing image data for the template example.\n \nSample Frames:\nEmpty = %d\n"
                + "Lower density = %d\nHigher density = %d\n",
            sampler.getNumberOfEmptySamples(), sampler.getNumberOfLowDensitySamples(),
            sampler.getNumberOfHighDensitySamples());

        // Turn off the recorder when the dialog is showing
        final boolean record = Recorder.record;
        Recorder.record = false;
        final NonBlockingGenericDialog gd = new NonBlockingGenericDialog(TITLE);
        gd.addMessage(msg);
        // gd.enableYesNoCancel(" Save ", " Resample ");
        gd.addSlider(keyNo, 0, 10, nNo);
        gd.addSlider(keyLow, 0, 10, nLow);
        gd.addSlider(keyHigh, 0, 10, nHigh);
        gd.addDialogListener(new DialogListener() {
          @Override
          public boolean dialogItemChanged(GenericDialog gd, AWTEvent event) {
            // If the event is null then this is the final call when the
            // dialog has been closed. We ignore this to prevent generating a
            // image the user has not seen.
            if (event == null) {
              return true;
            }
            nNo = (int) gd.getNextNumber();
            nLow = (int) gd.getNextNumber();
            nHigh = (int) gd.getNextNumber();
            out[0] = sampler.getSample(nNo, nLow, nHigh);
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
          }
        });
        gd.showDialog();
        if (gd.wasCanceled()) {
          out[0] = null;
          nNo = nLow = nHigh = 0; // For the recorder
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
          Recorder.recordOption(keyNo, Integer.toString(nNo));
          Recorder.recordOption(keyLow, Integer.toString(nLow));
          Recorder.recordOption(keyHigh, Integer.toString(nHigh));
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

  private static ImagePlus display(ImagePlus example, WindowOrganiser windowOrganiser) {
    final String title = "Template Example";
    // Update the example image
    example.setTitle(title);

    // Display as a new image. This is so we can close it later.
    return ConfigurationTemplate.displayTemplate(title, example, windowOrganiser);
  }

  private static void getNotes(TemplateSettings.Builder settings, String topFilterSummary) {
    settings.addNotes("Benchmark template");
    if (!TextUtils.isNullOrEmpty(resultsTitle)) {
      addField(settings, "Filter Analysis Title", resultsTitle);
    }
    // Add create data settings.
    // Just add the columns and the data from the summary window
    final String header = createResultsHeader(true);
    addField(settings, "Filter Analysis Summary Fields", header);
    addField(settings, "Filter Analysis Summary Values", topFilterSummary);
    // Now pick out key values...
    addKeyFields(settings, header, topFilterSummary, new String[] {"Molecules", "Density", "SNR",
        "s (nm)", "a (nm)", "Lower D", "Upper D", "Lower factor", "Upper factor"});

    // Add any other settings that may be useful in the template
    addField(settings, "Created", getCurrentTimeStamp());
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
    final String strDate = sdfDate.format(now);
    return strDate;
  }

  /**
   * Depth analysis.
   *
   * @param allAssignments The assignments generated from running the filter (or null)
   * @param filter the filter
   * @return the assignments
   */
  private ArrayList<FractionalAssignment[]>
      depthAnalysis(ArrayList<FractionalAssignment[]> allAssignments, DirectFilter filter) {
    // TODO : This analysis ignores the partial match distance.
    // Use the score for each result to get a weighted histogram.

    if (!depthRecallAnalysis || simulationParameters.fixedDepth) {
      return null;
    }

    // Build a histogram of the number of spots at different depths
    final double[] depths = depthStats.getValues();
    double[] limits = MathUtils.limits(depths);

    // final int bins = Math.max(10, nActual / 100);
    // final int bins = HistogramPlot.getBinsSturges(depths.length);
    final int bins = HistogramPlot.getBinsSqrtRule(depths.length);
    final double[][] h1 = HistogramPlot.calcHistogram(depths, limits[0], limits[1], bins);
    final double[][] h2 =
        HistogramPlot.calcHistogram(depthFitStats.getValues(), limits[0], limits[1], bins);

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
    final Plot2 plot1 = new Plot2(title1, "Depth (nm)", "Frequency");
    plot1.setLimits(limits[0], limits[1], 0, MathUtils.max(h1[1]));
    plot1.setColor(Color.black);
    plot1.addPoints(h1[0], h1[1], Plot2.BAR);
    plot1.addLabel(0, 0, "Black = Spots; Blue = Fitted; Red = Filtered");
    plot1.setColor(Color.blue);
    plot1.addPoints(h1[0], h2[1], Plot2.BAR);
    plot1.setColor(Color.red);
    plot1.addPoints(h1[0], h3[1], Plot2.BAR);
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

    final double halfSummaryDepth = summaryDepth * 0.5;

    final String title2 = TITLE + " Depth Histogram (normalised)";
    final Plot2 plot2 = new Plot2(title2, "Depth (nm)", "Recall");
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
  private ArrayList<FractionalAssignment[]>
      scoreAnalysis(ArrayList<FractionalAssignment[]> allAssignments, DirectFilter filter) {
    if (!scoreAnalysis) {
      return null;
    }

    // Build a histogram of the fitted spots that were available to be scored
    final double[] signal = signalFactorStats.getValues();
    final double[] distance = distanceStats.getValues();

    double[] limits1;
    if (BenchmarkSpotFit.signalFactor > 0 && upperSignalFactor > 0) {
      final double range = BenchmarkSpotFit.signalFactor * upperSignalFactor / 100.0;
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
    if (BenchmarkSpotFit.distanceInPixels > 0 && upperMatchDistance > 0) {
      final double range = simulationParameters.pixelPitch * BenchmarkSpotFit.distanceInPixels
          * upperMatchDistance / 100.0;
      limits2 = new double[] {0, range};
    } else {
      limits2 = MathUtils.limits(distance);
    }

    // final int bins = Math.max(10, nActual / 100);
    // final int bins = HistogramPlot.getBinsSturges(signal.length);
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
    final Plot2 plot2 = new Plot2(title2, "Distance (nm)", "Frequency");
    plot2.setLimits(limits2[0], limits2[1], 0, MathUtils.maxDefault(MathUtils.max(h2[1]), h2b[1]));
    plot2.setColor(Color.black);
    plot2.addLabel(0, 0, String.format("Blue = Fitted (%s); Red = Filtered (%s)",
        MathUtils.rounded(distanceStats.getMean()), MathUtils.rounded(sumDistance / count)));
    plot2.setColor(Color.blue);
    plot2.addPoints(h2[0], h2[1], Plot2.BAR);
    plot2.setColor(Color.red);
    plot2.addPoints(h2b[0], h2b[1], Plot2.BAR);
    ImageJUtils.display(title2, plot2, wo);

    // Draw signal factor histogram
    final String title1 = TITLE + " Signal Factor Histogram";
    final Plot2 plot1 = new Plot2(title1, "Signal Factor", "Frequency");
    plot1.setLimits(limits1[0], limits1[1], 0, MathUtils.maxDefault(MathUtils.max(h1[1]), h1b[1]));
    plot1.setColor(Color.black);
    plot1.addLabel(0, 0, String.format("Blue = Fitted (%s); Red = Filtered (%s)",
        MathUtils.rounded(signalFactorStats.getMean()), MathUtils.rounded(sumSignal / count)));
    plot1.setColor(Color.blue);
    plot1.addPoints(h1[0], h1[1], Plot2.BAR);
    plot1.setColor(Color.red);
    plot1.addPoints(h1b[0], h1b[1], Plot2.BAR);
    ImageJUtils.display(title1, plot1, wo);

    return allAssignments;
  }

  private void componentAnalysis(ComplexFilterScore bestFilterScore) {
    if (componentAnalysis == 0) {
      return;
    }
    final Consumer<String> output = createComponentAnalysisWindow();

    final DirectFilter bestFilter = bestFilterScore.getFilter();
    final String[] names = getNames(bestFilter);

    // Skip disabled parameters
    int paramCount = bestFilter.getNumberOfParameters();
    final boolean[] enable = new boolean[paramCount];
    final int[] map = new int[paramCount];
    for (int n = 0, i = 0; n < map.length; n++) {
      if (bestFilter.getParameterValue(n) == bestFilter.getDisabledParameterValue(n)) {
        enable[n] = true; // Mark to ignore
        paramCount--;
      } else {
        map[i++] = n;
      }
    }

    // Score the best filter just so we have the unique Ids of the results
    scoreComponents(bestFilter, -1, paramCount, null, null, null, 0);
    final int[] uniqueIds1 = uniqueIds;
    final int uniqueIdCount1 = uniqueIdCount;

    // Limit to 12 params == 4095 combinations (this is the max for two multi filters combined)
    if (componentAnalysis >= 3 && paramCount <= 12) {
      // Enumeration of combinations
      final long count = countComponentCombinations(paramCount);

      // Enumerate all combinations
      final ComplexFilterScore[] scores = new ComplexFilterScore[(int) count];
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
          final DirectFilter f = (DirectFilter) bestFilter.create(enable2);
          scores[index++] =
              scoreComponents(f, 0, k, combinations, enable2, uniqueIds1, uniqueIdCount1);
        }
      }

      // Report
      Arrays.sort(scores, FilterScoreCompararor.INSTANCE);

      int lastSize = 0;
      for (int i = 0; i < scores.length; i++) {
        if (componentAnalysis == 3) {
          if (lastSize == scores[i].size) {
            // Only add the best result for each size
            continue;
          }
          lastSize = scores[i].size;
        }

        // Find the last component that was added
        if (scores[i].size == 1) {
          scores[i].index = scores[i].combinations[0];
        } else {
          // For each size k, find the best result with k-1 components and set the add index
          // appropriately
          int add = -1;
          int target = -1;
          for (int l = 0; l < enable.length; l++) {
            if (scores[i].enable[l]) {
              target++;
            }
          }
          final int size1 = scores[i].size - 1;
          for (int ii = 0; ii < i; ii++) {
            if (scores[ii].size < size1) {
              continue;
            }
            if (scores[ii].size > size1) {
              break; // Broken
            }
            // Count matches. It must be 1 less than the current result
            int matches = 0;
            for (int l = 0; l < enable.length; l++) {
              if (scores[i].enable[l] && scores[ii].enable[l]) {
                matches++;
              }
            }
            if (matches == target) {
              // Find the additional parameter added
              for (int l = 0; l < enable.length; l++) {
                if (scores[i].enable[l]) {
                  if (scores[ii].enable[l]) {
                    continue;
                  }
                  add = l;
                  break;
                }
              }
              break;
            }
          }
          scores[i].index = add;
        }

        addToComponentAnalysisWindow(output, scores[i], bestFilterScore, names);
      }

      return;
    }

    // Preserve the option to output the best or all results if we fell through from above
    final int myComponentAnalysis =
        (componentAnalysis >= 3) ? componentAnalysis - 2 : componentAnalysis;

    // Progressively add components until all are the same as the input bestFilter
    int enabled = 0;
    int[] previousCombinations = new int[0];
    for (int ii = 0; ii < paramCount; ii++) {
      // Create a set of filters by enabling each component that is not currently enabled.
      final ComplexFilterScore[] scores = new ComplexFilterScore[paramCount - enabled];
      final int k = enabled + 1;
      for (int i = 0, j = 0; i < paramCount; i++) {
        final int n = map[i];
        if (enable[n]) {
          continue;
        }
        enable[n] = true;
        final DirectFilter f = (DirectFilter) bestFilter.create(enable);
        enable[n] = false;
        final int[] combinations = new int[k];
        System.arraycopy(previousCombinations, 0, combinations, 0, previousCombinations.length);
        combinations[k - 1] = n;
        Arrays.sort(combinations);
        scores[j++] = scoreComponents(f, n, k, combinations, null, uniqueIds1, uniqueIdCount1);
      }

      // Rank them
      Arrays.sort(scores);
      for (int i = 0; i < scores.length; i++) {
        addToComponentAnalysisWindow(output, scores[i], bestFilterScore, names);
        if (myComponentAnalysis == 1) {
          // Only add the best result
          break;
        }
      }

      // Flag the best component added
      enable[scores[0].index] = true;
      enabled++;
      previousCombinations = scores[0].combinations;
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

    // This returns the same as (2^n)-1
    // long total = 0;
    // for (int k = 1; k <= n; k++)
    // total += CombinatoricsUtils.binomialCoefficient(n, k);
    // return total;
  }

  private int[] uniqueIds;
  private int uniqueIdCount;

  private void setupFractionScoreStore() {
    uniqueIds = new int[maxUniqueId];
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

  private static class FilterScoreCompararor
      implements Comparator<ComplexFilterScore>, Serializable {
    private static final long serialVersionUID = 1L;
    /** The instance. */
    static final FilterScoreCompararor INSTANCE = new FilterScoreCompararor();

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

  private static class NamedPlot implements Comparable<NamedPlot> {
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

    @Override
    public int compareTo(NamedPlot other) {
      return Double.compare(other.score, score);
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
    TIntObjectHashMap<ArrayList<Coordinate>> previousResults;
    boolean canContinue = true;

    public IterationConvergenceChecker(FilterScore current) {
      // We have two different relative thresholds so use 2 convergence checkers,
      // one for the score and one for the filter sequence
      scoreChecker =
          new InterruptChecker(iterationScoreTolerance, iterationScoreTolerance * 1e-3, 0);
      filterChecker = new InterruptConvergenceChecker(iterationFilterTolerance,
          iterationFilterTolerance * 1e-3, false, true, iterationMaxIterations);
      if (iterationCompareResults) {
        previousResults = getResults(current);
      }
    }

    private TIntObjectHashMap<ArrayList<Coordinate>> getResults(FilterScore current) {
      return ResultsMatchCalculator
          .getCoordinates(createResults(null, (DirectFilter) current.filter, false));
    }

    public boolean converged(String prefix, FilterScore previous, FilterScore current,
        double[] previousParameters, double[] currentParameters) {
      // Stop if interrupted
      if (IJ.escapePressed()) {
        ImageJUtils.log("STOPPED");
        // Do not reset escape
        // IJ.resetEscape();
        canContinue = false;
        return true;
      }

      // Must converge on the non-filter parameters
      if (!converged(previousParameters, currentParameters)) {
        if (iterationCompareResults) {
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

      if (iterationCompareResults) {
        final TIntObjectHashMap<ArrayList<Coordinate>> currentResults = getResults(current);
        final MatchResult r = ResultsMatchCalculator.compareCoordinates(currentResults,
            previousResults, iterationCompareDistance);
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
      if (iterationMaxIterations != 0 && getIterations() >= iterationMaxIterations) {
        component = "iterations";
        canContinue = false;
      }
      ImageJUtils.log(prefix + " converged on " + component);
    }

    public int getIterations() {
      return filterChecker.getIterations();
    }
  }

  // Used to implement the FitnessFunction interface
  private String gaStatusPrefix = "";
  private Population<FilterScore> gaPopulation;

  // Used to set the strength on a filter
  private double[] strengthLower;
  private double[] strengthUpper;

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

  // Used for the scoring of filter sets
  private MultiPathFitResults[] gaResultsList;
  private MultiPathFitResults[] gaResultsListToScore;
  private boolean gaSubset;
  private int gaIteration;
  private DirectFilter searchScoreFilter;
  private FilterScoreResult[] gaScoreResults;
  private int gaScoreIndex;

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

    public ScoreWorker(BlockingQueue<ScoreJob> jobs, FilterScoreResult[] scoreResults,
        boolean createTextResult, CoordinateStore coordinateStore) {
      this.jobs = jobs;
      this.scoreResults = scoreResults;
      this.createTextResult = createTextResult;
      this.minFilter = (minimalFilter != null) ? (DirectFilter) minimalFilter.clone() : null;
      this.coordinateStore = coordinateStore;
    }

    @Override
    public void run() {
      try {
        while (true) {
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
        Logger.getLogger(BenchmarkFilterAnalysis.class.getName()).log(Level.WARNING, "Interrupted!",
            ex);
        Thread.currentThread().interrupt();
      } finally {
        finished = true;
      }
    }

    private void run(ScoreJob job) {
      if (ImageJUtils.isInterrupted()) {
        finished = true;
        return;
      }
      showProgress();
      // Directly write to the result array, this is thread safe
      scoreResults[job.index] =
          scoreFilter(job.filter, minFilter, createTextResult, coordinateStore);
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

    public ParameterScoreWorker(BlockingQueue<ParameterScoreJob> jobs,
        ParameterScoreResult[] scoreResults, boolean createTextResult) {
      this.jobs = jobs;
      this.scoreResults = scoreResults;
      this.createTextResult = createTextResult;
      this.filter = (DirectFilter) searchScoreFilter.clone();
      this.minFilter = (minimalFilter != null) ? (DirectFilter) minimalFilter.clone() : null;
      getBounds();
      this.gridCoordinateStore =
          new GridCoordinateStore(bounds.x, bounds.y, bounds.width, bounds.height, 0, 0);
    }

    @Override
    public void run() {
      try {
        while (true) {
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
        Logger.getLogger(BenchmarkFilterAnalysis.class.getName()).log(Level.WARNING, "Interrupted!",
            ex);
        Thread.currentThread().interrupt();
      } finally {
        finished = true;
      }
    }

    private void run(ParameterScoreJob job) {
      if (ImageJUtils.isInterrupted()) {
        finished = true;
        return;
      }
      showProgress();

      final double[] point = job.point;
      final int failCount = (int) Math.round(point[0]);
      final double residualsThreshold = point[1];
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
      scoreResults[job.index] = scoreFilter(filter, minFilter, failCount, residualsThreshold,
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
    }
  }

  private int progress;
  private int stepProgress;
  private int totalProgress;

  /**
   * Show progress.
   */
  private synchronized void showProgress() {
    if (progress % stepProgress == 0) {
      if (ImageJUtils.showStatus("Scoring Filter: " + progress + " / " + totalProgress)) {
        IJ.showProgress(progress, totalProgress);
      }
    }
    progress++;
  }

  private FilterScoreResult[] scoreFilters(FilterSet filterSet, boolean createTextResult) {
    if (filterSet.size() == 0) {
      return null;
    }

    initialiseScoring(filterSet);

    FilterScoreResult[] scoreResults = new FilterScoreResult[filterSet.size()];

    if (scoreResults.length == 1) {
      // No need to multi-thread this
      scoreResults[0] = scoreFilter((DirectFilter) filterSet.getFilters().get(0), minimalFilter,
          createTextResult, coordinateStore);
    } else {
      // Multi-thread score all the result
      final int nThreads = getThreads(scoreResults.length);
      final BlockingQueue<ScoreJob> jobs = new ArrayBlockingQueue<>(nThreads * 2);
      final List<ScoreWorker> workers = new LinkedList<>();
      final List<Thread> threads = new LinkedList<>();
      for (int i = 0; i < nThreads; i++) {
        final ScoreWorker worker = new ScoreWorker(jobs, scoreResults, createTextResult,
            (coordinateStore == null) ? null : coordinateStore.newInstance());
        final Thread t = new Thread(worker);
        workers.add(worker);
        threads.add(t);
        t.start();
      }

      int index = 0;
      totalProgress = scoreResults.length;
      stepProgress = ImageJUtils.getProgressInterval(totalProgress);
      progress = 0;
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
      IJ.showProgress(1);

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
      scoreResults[0] = scoreFilter(searchScoreFilter, minimalFilter, failCount, residualsThreshold,
          duplicateDistance, createCoordinateStore(duplicateDistance), createTextResult);
    } else {
      // Multi-thread score all the result
      final int nThreads = getThreads(scoreResults.length);
      final BlockingQueue<ParameterScoreJob> jobs = new ArrayBlockingQueue<>(nThreads * 2);
      final List<ParameterScoreWorker> workers = new LinkedList<>();
      final List<Thread> threads = new LinkedList<>();
      for (int i = 0; i < nThreads; i++) {
        final ParameterScoreWorker worker =
            new ParameterScoreWorker(jobs, scoreResults, createTextResult);
        final Thread t = new Thread(worker);
        workers.add(worker);
        threads.add(t);
        t.start();
      }

      totalProgress = scoreResults.length;
      stepProgress = ImageJUtils.getProgressInterval(totalProgress);
      progress = 0;
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
      IJ.showProgress(1);

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
      // Plot2 plot = new Plot2(title, "Precision", "Frequency");
      // // Find limits
      // double[] xlimit = Maths.limits(h1[0]);
      // plot.setLimits(xlimit[0] - 1, xlimit[1] + 1, 0, Maths.max(h1[1]) * 1.05);
      // plot.addPoints(h1[0], h1[1], Plot2.BAR);
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

      gaResultsListToScore = createMpf((DirectFilter) weakest, minimalFilter)
          .filterSubset(gaResultsList, ConsecutiveFailCounter.create(failCount), true);

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
        scoreFilter(filter, minimalFilter, gaResultsList, coordinateStore);

    final StringBuilder text = createResult(filter, r);
    add(text, gaIteration);
    gaWindow.accept(text.toString());
  }

  private SimpleFilterScore filterScoreOptimum;
  private SimpleParameterScore parameterScoreOptimum;

  private class ParameterScoreFunction implements FullScoreFunction<FilterScore> {
    @Override
    public SearchResult<FilterScore> findOptimum(double[][] points) {
      gaIteration++;
      SimpleParameterScore max = parameterScoreOptimum;

      // Sort points to allow the CoordinateStore to be reused with the same duplicate distance
      Arrays.sort(points, (o1, o2) -> {
        return BenchmarkFilterAnalysis.compare(o1[2], o2[2]);
      });

      final ParameterScoreResult[] scoreResults = scoreFilters(points, showResultsTable);

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

      if (showResultsTable) {
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

      // System.out.println(Arrays.toString(parameters));

      final MultiPathFilter multiPathFilter =
          new MultiPathFilter(searchScoreFilter, minimalFilter, residualsThreshold);
      final FractionClassificationResult r = multiPathFilter.fractionScoreSubset(
          gaResultsListToScore, ConsecutiveFailCounter.create(failCount), nActual, null, null,
          createCoordinateStore(duplicateDistance));

      final StringBuilder text = createResult(searchScoreFilter, r,
          buildResultsPrefix2(failCount, residualsThreshold, duplicateDistance));
      add(text, gaIteration);
      gaWindow.accept(text.toString());

      return new SearchResult<>(parameters, max);
    }

    @Override
    public SearchResult<FilterScore>[] score(double[][] points) {
      gaIteration++;
      SimpleParameterScore max = parameterScoreOptimum;

      // Sort points to allow the CoordinateStore to be reused with the same duplicate distance
      Arrays.sort(points, (o1, o2) -> BenchmarkFilterAnalysis.compare(o1[2], o2[2]));

      final ParameterScoreResult[] scoreResults = scoreFilters(points, showResultsTable);

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

      if (showResultsTable) {
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
          new MultiPathFilter(searchScoreFilter, minimalFilter, residualsThreshold);
      final FractionClassificationResult r = multiPathFilter.fractionScoreSubset(
          gaResultsListToScore, ConsecutiveFailCounter.create(failCount), nActual, null, null,
          createCoordinateStore(duplicateDistance));

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

  private static int compare(final double d1, final double d2) {
    if (d1 < d2) {
      return -1;
    }
    if (d1 > d2) {
      return 1;
    }
    return 0;
  }

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
        scoreFilter(filter, minimalFilter, gaResultsList, coordinateStore);

    final StringBuilder text = createResult(filter, r);
    add(text, gaIteration);
    gaWindow.accept(text.toString());

    return new SearchResult<>(filter.getParameters(), max);
  }

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
        scoreFilter(filter, minimalFilter, gaResultsList, coordinateStore);

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

  private double limit;

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
    if (!updateConfiguration(config, useLatest)) {
      return false;
    }
    return true;
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
    if (!BenchmarkSpotFilter.updateConfiguration(config)) {
      return false;
    }
    return true;
  }

  /**
   * Updates the given configuration using the latest settings used in benchmarking filtering. The
   * residuals threshold will be copied only if the input FitConfiguration has isComputeResiduals()
   * set to true.
   *
   * @param config the configuration
   * @param useLatest Use the latest best filter. Otherwise use the highest scoring.
   * @return true, if successful
   */
  public static boolean updateConfiguration(FitEngineConfiguration config, boolean useLatest) {
    if (scores.isEmpty()) {
      return false;
    }

    FilterResult best;
    if (useLatest) {
      best = scores.get(scores.size() - 1);
    } else {
      best = getBestResult();
    }

    // New smart filter support
    final FitConfiguration fitConfig = config.getFitConfiguration();
    fitConfig.setDirectFilter(best.getFilter());

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

  private static FilterResult getBestResult() {
    if (!scores.isEmpty()) {
      // Clone so we can sort the list by fail count
      @SuppressWarnings("unchecked")
      final ArrayList<FilterResult> scores2 = (ArrayList<FilterResult>) scores.clone();
      Collections.sort(scores2, new Comparator<FilterResult>() {
        @Override
        public int compare(FilterResult o1, FilterResult o2) {
          return o1.failCount - o2.failCount;
        }
      });

      int maxi = 0;
      double max = 0;
      for (int i = 0; i < scores2.size(); i++) {
        final double score = scores2.get(i).score;
        if (max >= score) {
          continue;
        }

        // Round this so that small differences are ignored.
        // This should favour filters with lower fail count.
        final double diff = MathUtils.round(score - max, 3);
        if (diff <= 0) {
          continue;
        }

        max = score;
        maxi = i;
      }
      return scores2.get(maxi);
    }
    return null;
  }

  private static boolean isShowOverlay() {
    return (showTP || showFP || showFN);
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
   * Show overlay.
   *
   * @param allAssignments The assignments generated from running the filter (or null)
   * @param filter the filter
   * @return The results from running the filter (or null)
   */
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
    // int tp = 0, fp = 0, fn = 0;
    for (final FractionalAssignment[] assignments : allAssignments) {
      if (assignments == null || assignments.length == 0) {
        continue;
      }
      float[] tx = null;
      float[] ty = null;
      int count = 0;
      // tp += assignments.length;
      if (showTP) {
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
        if (showTP) {
          tx[count] = spot.getX();
          ty[count++] = spot.getY();
        }
      }
      if (showTP) {
        SpotFinderPreview.addRoi(frame, o, tx, ty, count, Color.green);
      }
    }

    float[] x = new float[10];
    float[] y = new float[x.length];
    float[] x2 = new float[10];
    float[] y2 = new float[x2.length];

    // Do FP (all remaining results that are not a TP)
    PreprocessedPeakResult[] filterResults = null;
    if (showFP) {
      final MultiPathFilter multiPathFilter = createMpf(filter, minimalFilter);
      // multiPathFilter.setDebugFile("/tmp/filter.txt");

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
      // fp += c;
      if (c1 != 0) {
        SpotFinderPreview.addRoi(frame, o, x, y, c1, Color.red);
      }
      if (c2 != 0) {
        SpotFinderPreview.addRoi(frame, o, x2, y2, c2, Color.magenta);
      }
    }

    // Do TN (all remaining peaks that have not been matched)
    if (showFN) {
      final boolean checkBorder = (BenchmarkSpotFilter.lastAnalysisBorder != null
          && BenchmarkSpotFilter.lastAnalysisBorder.x != 0);
      final float border;
      final float xlimit;
      final float ylimit;
      if (checkBorder) {
        final Rectangle lastAnalysisBorder = BenchmarkSpotFilter.lastAnalysisBorder;
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
          // fn += c;
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

    // System.out.printf("TP=%d, FP=%d, FN=%d, N=%d (%d)\n", tp, fp, fn, tp + fn, results.size());

    imp.setOverlay(o);

    return filterResults;
  }

  /**
   * Filter using border.
   *
   * @param results the results
   * @return the results that are inside the border
   */
  private static UniqueIdPeakResult[] filterUsingBorder(UniqueIdPeakResult[] results) {
    // Respect the border used in BenchmarkSpotFilter?
    if (BenchmarkSpotFilter.lastAnalysisBorder == null
        || BenchmarkSpotFilter.lastAnalysisBorder.x == 0) {
      return results;
    }
    // Just implement a hard border
    // Note: This is mainly relevant when loading in simulated ground truth data.
    // The simulations performed by Create Data already use a border when doing a
    // random distribution.
    final Rectangle lastAnalysisBorder = BenchmarkSpotFilter.lastAnalysisBorder;
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
      // System.out.printf("Removed %d\n", results.length - j);
      return Arrays.copyOf(results2, count);
    }
    return results;
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
    if (result.getXPosition() < border || result.getXPosition() > xlimit
        || result.getYPosition() < border || result.getYPosition() > ylimit) {
      return true;
    }
    return false;
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
    if (result.getX() < border || result.getX() > xlimit || result.getY() < border
        || result.getY() > ylimit) {
      return true;
    }
    return false;
  }

  /**
   * Save the results to memory.
   *
   * @param filterResults The results from running the filter (or null)
   * @param filter the filter
   */
  private void saveResults(PreprocessedPeakResult[] filterResults, DirectFilter filter) {
    final MemoryPeakResults results = createResults(filterResults, filter, true);
    MemoryPeakResults.addResults(results);
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
      final MultiPathFilter multiPathFilter = createMpf(filter, minimalFilter);
      // multiPathFilter.setDebugFile("/tmp/filter.txt");
      filterResults = filterResults(multiPathFilter);
    }

    final MemoryPeakResults results = new MemoryPeakResults();
    results.copySettings(this.results);
    results.setName(TITLE);

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

          results.add(frame, origX, origY, 0, 0, spot.getNoise(), spot.getMeanSignal(), params,
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

        results.add(frame, origX, origY, 0, 0, spot.getNoise(), spot.getMeanSignal(), params, null);
      }
    }

    return results;
  }

  private PreprocessedPeakResult[] filterResults(final MultiPathFilter multiPathFilter) {
    return multiPathFilter.filter(resultsList, ConsecutiveFailCounter.create(failCount), true,
        coordinateStore);
  }

  private CoordinateStore createCoordinateStore() {
    return createCoordinateStore(duplicateDistance);
  }

  private CoordinateStore createCoordinateStore(double duplicateDistance) {
    getBounds();
    duplicateDistance *= distanceScallingFactor;
    return CoordinateStoreFactory.create(bounds.x, bounds.y, bounds.width, bounds.height,
        duplicateDistance);
  }

  private Rectangle bounds;
  private double distanceScallingFactor;

  private Rectangle getBounds() {
    if (bounds == null) {
      bounds = createBounds();
    }
    return bounds;
  }

  private synchronized Rectangle createBounds() {
    if (bounds == null) {
      if (duplicateDistanceAbsolute) {
        distanceScallingFactor = 1;
      } else {
        // The duplicate distance is scaled
        distanceScallingFactor = BenchmarkSpotFit.config.getHwhmMax();
      }

      final ImagePlus imp = CreateData.getImage();
      if (imp != null) {
        return new Rectangle(0, 0, imp.getWidth(), imp.getHeight());
      }
      return results.getBounds(true);
    }
    return bounds;
  }
}
