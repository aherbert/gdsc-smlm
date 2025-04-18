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

package uk.ac.sussex.gdsc.smlm.ij.plugins.benchmark;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.plugin.PlugIn;
import ij.text.TextWindow;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.IntConsumer;
import java.util.logging.Level;
import org.apache.commons.lang3.concurrent.ConcurrentRuntimeException;
import org.apache.commons.lang3.tuple.Pair;
import uk.ac.sussex.gdsc.core.ij.BufferedTextWindow;
import uk.ac.sussex.gdsc.core.ij.ImageJPluginLoggerHelper;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.logging.Ticker;
import uk.ac.sussex.gdsc.core.match.ClassificationResult;
import uk.ac.sussex.gdsc.core.match.Coordinate;
import uk.ac.sussex.gdsc.core.match.FractionClassificationResult;
import uk.ac.sussex.gdsc.core.threshold.AutoThreshold;
import uk.ac.sussex.gdsc.core.threshold.FloatHistogram;
import uk.ac.sussex.gdsc.core.threshold.Histogram;
import uk.ac.sussex.gdsc.core.utils.ImageExtractor;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.OpenHashMaps.CustomInt2ObjectOpenHashMap;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.concurrent.ConcurrencyUtils;
import uk.ac.sussex.gdsc.smlm.engine.FitEngineConfiguration;
import uk.ac.sussex.gdsc.smlm.engine.FitWorker;
import uk.ac.sussex.gdsc.smlm.filters.Spot;
import uk.ac.sussex.gdsc.smlm.fitting.Gaussian2DFitter;
import uk.ac.sussex.gdsc.smlm.ij.plugins.HelpUrls;
import uk.ac.sussex.gdsc.smlm.ij.plugins.PsfCalculator;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsMatchCalculator;
import uk.ac.sussex.gdsc.smlm.ij.plugins.SmlmUsageTracker;
import uk.ac.sussex.gdsc.smlm.ij.plugins.benchmark.BenchmarkSpotFilter.BenchmarkSpotFilterResult;
import uk.ac.sussex.gdsc.smlm.ij.plugins.benchmark.BenchmarkSpotFilter.FilterResult;
import uk.ac.sussex.gdsc.smlm.ij.plugins.benchmark.BenchmarkSpotFilter.ScoredSpot;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.ij.utils.ImageJImageConverter;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResultPoint;

/**
 * Attempt to classify the spot candidates into those that do match a result and those that do not.
 * This smart ranking can be used to ensure all the good candidates are processed per frame before
 * the fail limit takes effect.
 */
public class BenchmarkSmartSpotRanking implements PlugIn {
  private static final String TITLE = "Smart Spot Ranking";

  private static final byte TP = (byte) 1;
  private static final byte FP = (byte) 2;
  private static final byte TN = (byte) 3;
  private static final byte FN = (byte) 4;

  private static final AtomicReference<TextWindow> SUMMARY_TABLE_REF = new AtomicReference<>();

  /** The reference to the latest fit engine configuration. */
  private static final AtomicReference<FitEngineConfiguration> CONFIG_REF =
      new AtomicReference<>(FitEngineConfiguration.create());

  /** The coordinate cache. This stores the coordinates for a simulation Id. */
  private static final AtomicReference<
      Pair<Integer, Int2ObjectOpenHashMap<List<Coordinate>>>> COORDINATE_CACHE =
          new AtomicReference<>(Pair.of(-1, null));

  private static final AtomicReference<CandidateData> CANDIDATE_DATA_CACHE =
      new AtomicReference<>();

  /** The fit engine configuration. */
  FitEngineConfiguration config;
  /** The threshold methods to compute. */
  AutoThreshold.Method[] methods;
  /** The SNR levels to compute. */
  double[] levels;
  private String[] methodNames;

  private boolean extraOptions;

  private MemoryPeakResults results;
  private CreateData.SimulationParameters simulationParameters;
  private BenchmarkSpotFilterResult filterResult;

  private ImagePlus imp;

  /** The plugin settings. */
  Settings settings;

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

    CandidateData(Int2ObjectOpenHashMap<FilterCandidates> filterCandidates, int filterId,
        double fractionPositive, double fractionNegative, int countPositive, int countNegative,
        Settings settings) {
      this.filterCandidates = filterCandidates;
      this.filterId = filterId;
      this.fractionPositive = fractionPositive;
      this.fractionNegative = fractionNegative;
      this.countPositive = countPositive;
      this.countNegative = countNegative;
      fractionPositives = settings.fractionPositives;
      fractionNegativesAfterAllPositives = settings.fractionNegativesAfterAllPositives;
      negativesAfterAllPositives = settings.negativesAfterAllPositives;
    }

    /**
     * Return true if the settings used in the candidate data are different.
     *
     * @param filterId the filter id
     * @param settings the settings
     * @return true if different
     */
    boolean differentSettings(int filterId, Settings settings) {
      return this.filterId != filterId || fractionPositives != settings.fractionPositives
          || fractionNegativesAfterAllPositives != settings.fractionNegativesAfterAllPositives
          || negativesAfterAllPositives != settings.negativesAfterAllPositives;
    }
  }

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    static final AutoThreshold.Method[] THRESHOLD_METHODS;
    static final double[] SNR_LEVELS;
    static final String[] THRESHOLD_METHOD_NAMES;
    static final String[] SORT_METHODS = {"(None)", "tp", "fp", "tn", "fn", "Precision", "Recall",
        "F0.5", "F1", "F2", "Jaccard", "MCC"};

    static {
      final DoubleArrayList list = new DoubleArrayList(20);
      for (int snr = 20; snr <= 70; snr += 5) {
        list.add(snr);
      }
      SNR_LEVELS = list.toDoubleArray();

      THRESHOLD_METHODS = AutoThreshold.Method.values();
      THRESHOLD_METHOD_NAMES = new String[THRESHOLD_METHODS.length + SNR_LEVELS.length];

      // Enable all methods
      int count = 0;
      while (count < THRESHOLD_METHODS.length) {
        THRESHOLD_METHOD_NAMES[count] = THRESHOLD_METHODS[count].toString();
        count++;
      }
      // Add signal-to-noise threshold methods
      for (final double snr : SNR_LEVELS) {
        THRESHOLD_METHOD_NAMES[count] = "SNR" + snr;
        count++;
      }
    }

    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> INSTANCE = new AtomicReference<>(new Settings());

    boolean[] thresholdMethodOptions;
    double fractionPositives;
    double fractionNegativesAfterAllPositives;
    int negativesAfterAllPositives;
    boolean selectMethods;
    int compactBins;
    int sortIndex;
    boolean useFractionScores;
    boolean showOverlay;

    Settings() {
      // Set defaults
      thresholdMethodOptions = new boolean[THRESHOLD_METHOD_NAMES.length];
      Arrays.fill(thresholdMethodOptions, true);
      // Turn some off
      // These often fail to converge
      thresholdMethodOptions[AutoThreshold.Method.INTERMODES.ordinal()] = false;
      thresholdMethodOptions[AutoThreshold.Method.MINIMUM.ordinal()] = false;
      // These are slow
      thresholdMethodOptions[AutoThreshold.Method.SHANBHAG.ordinal()] = false;
      thresholdMethodOptions[AutoThreshold.Method.RENYI_ENTROPY.ordinal()] = false;
      fractionPositives = 100;
      fractionNegativesAfterAllPositives = 50;
      negativesAfterAllPositives = 10;
      selectMethods = true;
      compactBins = 1024;
      sortIndex = SORT_METHODS.length - 3; // F2 to favour recall
      useFractionScores = true;
    }

    Settings(Settings source) {
      thresholdMethodOptions = source.thresholdMethodOptions.clone();
      fractionPositives = source.fractionPositives;
      fractionNegativesAfterAllPositives = source.fractionNegativesAfterAllPositives;
      negativesAfterAllPositives = source.negativesAfterAllPositives;
      selectMethods = source.selectMethods;
      compactBins = source.compactBins;
      sortIndex = source.sortIndex;
      useFractionScores = source.useFractionScores;
      showOverlay = source.showOverlay;
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
   * The filter candidate data.
   */
  private static class FilterCandidates {
    // Integer counts of positives (matches) and negatives
    final int pos;
    final int neg;
    // Double sums of the fractions match score and antiscore
    final double np;
    final double nn;
    final ScoredSpot[] spots;

    FilterCandidates(int pos, int neg, double np, double nn, ScoredSpot[] spots) {
      this.pos = pos;
      this.neg = neg;
      this.np = np;
      this.nn = nn;
      this.spots = spots;
    }
  }

  /**
   * The rank result data.
   */
  private static class RankResult {
    final float threshold;
    final FractionClassificationResult fresult;
    final ClassificationResult cresult;
    /**
     * Store details about the spots that were accepted.
     */
    final byte[] good;
    final long time;

    RankResult(float threshold, FractionClassificationResult fresult, ClassificationResult cresult,
        byte[] good, long time) {
      this.threshold = threshold;
      this.fresult = fresult;
      this.cresult = cresult;
      this.good = good;
      this.time = time;
    }
  }

  /**
   * The ranked results data.
   */
  private static class RankResults {
    final ScoredSpot[] spots;
    /**
     * Store the z-position of the actual spots for later analysis. Size is the number of actual
     * spots
     */
    final double[] zPosition;

    List<RankResult> results = new ArrayList<>();

    RankResults(ScoredSpot[] spots, double[] zPosition) {
      this.spots = spots;
      this.zPosition = zPosition;
    }
  }

  /**
   * Used to allow multi-threading of the fitting method.
   */
  private class Worker implements Runnable {
    volatile boolean finished;
    final BlockingQueue<Integer> jobs;
    final ImageStack stack;
    final Int2ObjectOpenHashMap<List<Coordinate>> actualCoordinates;
    final Int2ObjectOpenHashMap<FilterCandidates> filterCandidates;
    final Int2ObjectOpenHashMap<RankResults> results;
    final int fitting;
    final boolean requireSnr;
    final Ticker ticker;

    float[] data;
    double[] region;

    Worker(BlockingQueue<Integer> jobs, ImageStack stack,
        Int2ObjectOpenHashMap<List<Coordinate>> actualCoordinates,
        Int2ObjectOpenHashMap<FilterCandidates> filterCandidates, Ticker ticker) {
      this.jobs = jobs;
      this.stack = stack;
      this.actualCoordinates = actualCoordinates;
      this.filterCandidates = filterCandidates;
      this.results = new Int2ObjectOpenHashMap<>();
      this.ticker = ticker;
      fitting = config.getFittingWidth();
      requireSnr = (levels.length > 0);
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
            // Only run if not finished to allow queue to be emptied
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

    @SuppressWarnings("null")
    private void run(int frame) {
      if (ImageJUtils.isInterrupted()) {
        finished = true;
        return;
      }

      final FilterCandidates candidates = filterCandidates.get(frame);

      // Extract the spot intensities
      final ScoredSpot[] spots = candidates.spots;
      final float[] intensity = new float[spots.length];
      for (int i = 0; i < spots.length; i++) {
        final Spot spot = spots[i].spot;
        intensity[i] = spot.intensity;
      }

      final boolean simpleBackground = true;

      // Estimate SNR
      double[] snr = null;
      long timeSnr = 0;
      if (requireSnr) {
        timeSnr = System.nanoTime();
        data = ImageJImageConverter.getData(stack.getPixels(frame), stack.getWidth(),
            stack.getHeight(), null, data);
        final int maxx = stack.getWidth();
        final int maxy = stack.getHeight();
        final ImageExtractor ie = ImageExtractor.wrap(data, maxx, maxy);
        final float noise = FitWorker.estimateNoise(data, maxx, maxy, config.getNoiseMethod());
        snr = new double[spots.length];
        for (int i = 0; i < spots.length; i++) {
          final Spot spot = spots[i].spot;
          final Rectangle regionBounds = ie.getBoxRegionBounds(spot.x, spot.y, fitting);
          region = ie.crop(regionBounds, region);
          double sum = 0;
          final int width = regionBounds.width;
          final int height = regionBounds.height;
          final int size = width * height;
          for (int k = size; k-- > 0;) {
            sum += region[k];
          }

          final double b;
          if (simpleBackground) {
            // TODO - The number of peaks could use the other candidates in the fit region
            b = Gaussian2DFitter.getBackground(region, width, height, 1);
          } else {
            // Use a very wide region to find the local background with the lowest % of the smoothed
            // data
            final Rectangle regionBounds2 = ie.getBoxRegionBounds(spot.x, spot.y, fitting * 2);
            region = ie.crop(regionBounds2, region);
            final int width2 = regionBounds2.width;
            final int height2 = regionBounds2.height;
            int size2 = 0;
            for (int y = 0; y < height2; y++) {
              // Assume neighbour pixels should have equal noise and average them
              for (int x = 1, index = y * width2; x < width2; x += 2, index += 2) {
                region[size2++] = region[index] + region[index + 1];
              }
            }
            Arrays.sort(region, 0, size2);
            double sumB = 0;
            int count = 0;
            for (int k = (int) Math.ceil(size2 * 0.2); k-- > 0;) {
              sumB += region[k];
              count++;
            }
            b = sumB / (count * 2); // Account for averaging
          }

          final double signal = sum - b * size;
          snr[i] = signal / noise;
        }
        timeSnr = System.nanoTime() - timeSnr;
      }

      final Coordinate[] actual = ResultsMatchCalculator.getCoordinates(actualCoordinates, frame);
      final double[] zPosition = new double[actual.length];
      for (int i = 0; i < actual.length; i++) {
        final PeakResultPoint p = (PeakResultPoint) actual[i];
        zPosition[i] = p.getPeakResult().getZPosition();
      }

      final RankResults rankResults = new RankResults(spots, zPosition);
      this.results.put(frame, rankResults);

      long t1 = System.nanoTime();
      final FloatHistogram histogram = FloatHistogram.buildHistogram(intensity, true);
      // Only compact once
      final Histogram histogram2 = histogram.compact(settings.compactBins);
      t1 = System.nanoTime() - t1;

      for (final AutoThreshold.Method m : methods) {
        long t2 = System.nanoTime();
        final float t = histogram2.getAutoThreshold(m);
        t2 = System.nanoTime() - t2;

        final long time = (m == AutoThreshold.Method.NONE) ? 0 : t1 + t2;

        // Score
        final byte[] category = new byte[spots.length];
        double tp = 0;
        double fp = 0;
        double tn = 0;
        double fn = 0;
        int itp = 0;
        int ifp = 0;
        int itn = 0;
        int ifn = 0;
        for (int i = 0; i < spots.length; i++) {
          if (intensity[i] >= t) {
            tp += spots[i].getScore();
            fp += spots[i].antiScore();
            if (spots[i].match) {
              category[i] = TP;
              itp++;
            } else {
              category[i] = FP;
              ifp++;
            }
          } else {
            fn += spots[i].getScore();
            tn += spots[i].antiScore();
            if (spots[i].match) {
              category[i] = FN;
              ifn++;
            } else {
              category[i] = TN;
              itn++;
            }
          }
        }

        // Store the results using a copy of the original (to preserve the candidates for repeat
        // analysis)
        rankResults.results.add(new RankResult(t, new FractionClassificationResult(tp, fp, tn, fn),
            new ClassificationResult(itp, ifp, itn, ifn), category, time));
      }

      for (final double l : levels) {
        // Get the intensity of the lowest spot
        double minIntensity = Double.POSITIVE_INFINITY;

        // Score
        final byte[] category = new byte[spots.length];
        double tp = 0;
        double fp = 0;
        double tn = 0;
        double fn = 0;
        int itp = 0;
        int ifp = 0;
        int itn = 0;
        int ifn = 0;
        for (int i = 0; i < spots.length; i++) {
          if (snr[i] >= l) {
            tp += spots[i].getScore();
            fp += spots[i].antiScore();
            if (spots[i].match) {
              category[i] = TP;
              itp++;
            } else {
              category[i] = FP;
              ifp++;
            }
            if (minIntensity > intensity[i]) {
              minIntensity = intensity[i];
            }
          } else {
            fn += spots[i].getScore();
            tn += spots[i].antiScore();
            if (spots[i].match) {
              category[i] = FN;
              ifn++;
            } else {
              category[i] = TN;
              itn++;
            }
          }
        }

        // Store the results using a copy of the original (to preserve the candidates for repeat
        // analysis)
        rankResults.results.add(
            new RankResult((float) minIntensity, new FractionClassificationResult(tp, fp, tn, fn),
                new ClassificationResult(itp, ifp, itn, ifn), category, timeSnr));
      }
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    extraOptions = ImageJUtils.isExtraOptions();

    simulationParameters = CreateData.getSimulationParameters();
    if (simulationParameters == null) {
      IJ.error(TITLE, "No benchmark spot parameters in memory");
      return;
    }
    imp = CreateData.getImage();
    if (imp == null) {
      IJ.error(TITLE, "No simulation image");
      return;
    }
    results = CreateData.getResults();
    if (results == null) {
      IJ.error(TITLE, "No benchmark results in memory");
      return;
    }
    filterResult = BenchmarkSpotFilter.getBenchmarkFilterResult();
    if (filterResult == null) {
      IJ.error(TITLE, "No benchmark spot candidates in memory");
      return;
    }
    if (filterResult.simulationId != simulationParameters.id) {
      IJ.error(TITLE, "Update the benchmark spot candidates for the latest simulation");
      return;
    }

    if (!showDialog()) {
      return;
    }

    runAnalysis();
  }

  private boolean showDialog() {
    ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("smart-spot-ranking"));

    settings = Settings.load();
    config = CONFIG_REF.get().createCopy();

    ImageJUtils.addMessage(gd,
        "Rank candidate spots in the benchmark image created by " + CreateData.TITLE
            + " plugin\nand identified by the " + BenchmarkSpotFilter.TITLE
            + " plugin.\nPSF width = %s nm (Square pixel adjustment = %s nm)\n \n"
            + "Configure the fitting:",
        MathUtils.rounded(simulationParameters.sd), MathUtils.rounded(getSa()));

    gd.addSlider("Fraction_positives", 50, 100, settings.fractionPositives);
    gd.addSlider("Fraction_negatives_after_positives", 0, 100,
        settings.fractionNegativesAfterAllPositives);
    gd.addSlider("Min_negatives_after_positives", 0, 10, settings.negativesAfterAllPositives);
    gd.addCheckbox("Select_methods", settings.selectMethods);
    gd.addNumericField("Compact_bins", settings.compactBins, 0);
    gd.addChoice("Sort", Settings.SORT_METHODS, settings.sortIndex);
    gd.addCheckbox("Use_fraction_scores", settings.useFractionScores);

    // Collect options for fitting that may effect ranking
    final double sa = getSa();
    gd.addNumericField("Initial_StdDev", sa / simulationParameters.pixelPitch, 3);
    gd.addSlider("Fitting_width", 2, 4.5, config.getFitting());

    // Output options
    gd.addCheckbox("Show_overlay", settings.showOverlay);

    if (extraOptions) {
      gd.addChoice("Noise_method", SettingsManager.getNoiseEstimatorMethodNames(),
          config.getNoiseMethod().getNumber());
    }

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.fractionPositives = Math.abs(gd.getNextNumber());
    settings.fractionNegativesAfterAllPositives = Math.abs(gd.getNextNumber());
    settings.negativesAfterAllPositives = (int) Math.abs(gd.getNextNumber());
    settings.selectMethods = gd.getNextBoolean();
    settings.compactBins = (int) Math.abs(gd.getNextNumber());
    settings.sortIndex = gd.getNextChoiceIndex();
    settings.useFractionScores = gd.getNextBoolean();

    // Collect options for fitting that may effect ranking
    config.getFitConfiguration().setInitialPeakStdDev(gd.getNextNumber());
    config.setFitting(gd.getNextNumber());

    settings.showOverlay = gd.getNextBoolean();

    if (extraOptions) {
      config
          .setNoiseMethod(SettingsManager.getNoiseEstimatorMethodValues()[gd.getNextChoiceIndex()]);
    }

    settings.save();
    CONFIG_REF.set(config);

    if (gd.invalidNumber()) {
      return false;
    }

    methodNames = Settings.THRESHOLD_METHOD_NAMES.clone();
    if (settings.selectMethods) {
      gd = new ExtendedGenericDialog(TITLE);
      gd.addHelp(HelpUrls.getUrl("smart-spot-ranking"));
      for (int i = 0; i < Settings.THRESHOLD_METHOD_NAMES.length; i++) {
        gd.addCheckbox(Settings.THRESHOLD_METHOD_NAMES[i], settings.thresholdMethodOptions[i]);
      }

      gd.showDialog();

      if (gd.wasCanceled()) {
        return false;
      }

      int count = 0;
      int count1 = 0;
      int count2 = 0;
      methods = new AutoThreshold.Method[Settings.THRESHOLD_METHODS.length];
      levels = new double[Settings.SNR_LEVELS.length];
      for (int i = 0, j = 0; i < Settings.THRESHOLD_METHOD_NAMES.length; i++) {
        settings.thresholdMethodOptions[i] = gd.getNextBoolean();
        if (settings.thresholdMethodOptions[i]) {
          methodNames[count++] = Settings.THRESHOLD_METHOD_NAMES[i];
          if (i < Settings.THRESHOLD_METHODS.length) {
            methods[count1++] = Settings.THRESHOLD_METHODS[i];
          } else {
            levels[count2++] = Settings.SNR_LEVELS[j++];
          }
        }
      }

      methodNames = Arrays.copyOf(methodNames, count);
      methods = Arrays.copyOf(methods, count1);
      levels = Arrays.copyOf(levels, count2);
    } else {
      // Do them all
      methods = Settings.THRESHOLD_METHODS.clone();
      levels = Settings.SNR_LEVELS.clone();
    }

    if (methodNames.length == 0) {
      IJ.error(TITLE, "No methods selected");
      return false;
    }

    return true;
  }

  private void runAnalysis() {
    // Extract all the results in memory into a list per frame. This can be cached
    boolean refresh = false;
    final Pair<Integer, Int2ObjectOpenHashMap<List<Coordinate>>> coords = COORDINATE_CACHE.get();

    Int2ObjectOpenHashMap<List<Coordinate>> actualCoordinates;
    if (coords.getKey() != simulationParameters.id) {
      // Do not get integer coordinates
      // The Coordinate objects will be PeakResultPoint objects that store the original PeakResult
      // from the MemoryPeakResults
      actualCoordinates = ResultsMatchCalculator.getCoordinates(results, false);
      COORDINATE_CACHE.set(Pair.of(simulationParameters.id, actualCoordinates));
      refresh = true;
    } else {
      actualCoordinates = coords.getValue();
    }

    // Extract all the candidates into a list per frame. This can be cached if the settings have not
    // changed.
    CandidateData candidateData = CANDIDATE_DATA_CACHE.get();
    if (refresh || candidateData == null
        || candidateData.differentSettings(filterResult.id, settings)) {
      candidateData = subsetFilterResults(filterResult.filterResults);
      CANDIDATE_DATA_CACHE.set(candidateData);
    }

    final Int2ObjectOpenHashMap<FilterCandidates> filterCandidates = candidateData.filterCandidates;
    final ImageStack stack = imp.getImageStack();

    // Create a pool of workers
    final int nThreads = Prefs.getThreads();
    final BlockingQueue<Integer> jobs = new ArrayBlockingQueue<>(nThreads * 2);
    final List<Worker> workers = new LinkedList<>();
    final List<Thread> threads = new LinkedList<>();
    final Ticker ticker = ImageJUtils.createTicker(filterCandidates.size(), nThreads);
    for (int i = 0; i < nThreads; i++) {
      final Worker worker = new Worker(jobs, stack, actualCoordinates, filterCandidates, ticker);
      final Thread t = new Thread(worker);
      workers.add(worker);
      threads.add(t);
      t.start();
    }

    // Process the frames
    filterCandidates.keySet().forEach((IntConsumer) value -> put(jobs, value));
    // Finish all the worker threads by passing in a null job
    for (int i = threads.size(); i-- != 0;) {
      put(jobs, -1);
    }

    // Wait for all to finish
    for (final Thread thread : threads) {
      try {
        thread.join();
      } catch (final InterruptedException ex) {
        Thread.currentThread().interrupt();
        throw new ConcurrentRuntimeException("Unexpected interrupt", ex);
      }
    }
    threads.clear();

    IJ.showProgress(1);

    if (ImageJUtils.isInterrupted()) {
      IJ.showStatus("Aborted");
      return;
    }

    IJ.showStatus("Collecting results ...");

    final CustomInt2ObjectOpenHashMap<RankResults> rankResults =
        new CustomInt2ObjectOpenHashMap<>();
    for (final Worker w : workers) {
      rankResults.putAll(w.results);
    }

    summariseResults(rankResults, candidateData);

    IJ.showStatus("");
  }

  /**
   * Extract all the filter candidates in order until the desired number of positives have been
   * reached and the number of negatives matches the configured parameters.
   *
   * @param filterResults the filter results
   * @return The filter candidate data
   */
  private CandidateData
      subsetFilterResults(CustomInt2ObjectOpenHashMap<FilterResult> filterResults) {
    // Convert fractions from percent
    final double f1 = Math.min(1, settings.fractionPositives / 100.0);
    final double f2 = settings.fractionNegativesAfterAllPositives / 100.0;

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
      int negAfter = 0;

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
            negAfter++;
          }
        }

        if (reachedTarget
            // Check if we have reached both the limits
            && negAfter >= settings.negativesAfterAllPositives
            && (double) neg / (neg + pos) >= f2) {
          break;
        }
      }

      // TODO - This is different from BenchmarkSpotFit where all the candidates are
      // included but only the first N are processed. Should this be changed here too.

      subset.put(frame, new FilterCandidates(pos, neg, np, nn, Arrays.copyOf(result.spots, count)));
    });

    return new CandidateData(subset, filterResult.id, fX[0], fX[1], nX[0], nX[1], settings);
  }

  private static void put(BlockingQueue<Integer> jobs, int index) {
    try {
      jobs.put(index);
    } catch (final InterruptedException ex) {
      ImageJPluginLoggerHelper.getLogger(BenchmarkSmartSpotRanking.class).log(Level.WARNING,
          "Unexpected interruption", ex);
      Thread.currentThread().interrupt();
    }
  }

  /**
   * The scored result data.
   */
  private static class ScoredResult {
    int index;
    double score;
    String result;

    ScoredResult(int index, double score, String result) {
      this.index = index;
      this.score = score;
      this.result = result;
    }

    /**
     * Compare the two results.
     *
     * @param r1 the first result
     * @param r2 the second result
     * @return -1, 0 or 1
     */
    static int compare(ScoredResult r1, ScoredResult r2) {
      return Double.compare(r2.score, r1.score);
    }
  }

  private void summariseResults(CustomInt2ObjectOpenHashMap<RankResults> rankResults,
      CandidateData candidateData) {
    // Summarise the ranking results.
    final StringBuilder sb = new StringBuilder(filterResult.resultPrefix);

    // countPositive and countNegative is the fractional score of the spot candidates
    addCount(sb, (double) candidateData.countPositive + candidateData.countNegative);
    addCount(sb, candidateData.countPositive);
    addCount(sb, candidateData.countNegative);
    addCount(sb, candidateData.fractionPositive);
    addCount(sb, candidateData.fractionNegative);

    final double[] counter1 = new double[2];
    final int[] counter2 = new int[2];
    candidateData.filterCandidates.values().forEach(result -> {
      counter1[0] += result.np;
      counter1[1] += result.nn;
      counter2[0] += result.pos;
      counter2[1] += result.neg;
    });
    final int countTp = counter2[0];
    final int countFp = counter2[1];

    // The fraction of positive and negative candidates that were included
    add(sb, (100.0 * countTp) / candidateData.countPositive);
    add(sb, (100.0 * countFp) / candidateData.countNegative);

    // Add counts of the candidates
    add(sb, countTp + countFp);
    add(sb, countTp);
    add(sb, countFp);

    // Add fractional counts of the candidates
    double tp = counter1[0];
    double fp = counter1[1];
    add(sb, tp + fp);
    add(sb, tp);
    add(sb, fp);

    // Materialise rankResults
    final int[] frames = new int[rankResults.size()];
    final RankResults[] rankResultsArray = new RankResults[rankResults.size()];
    final int[] counter = new int[1];
    rankResults.forEach((int key, RankResults value) -> {
      frames[counter[0]] = key;
      rankResultsArray[counter[0]] = value;
      counter[0]++;
    });

    // Summarise actual and candidate spots per frame
    final Statistics actual = new Statistics();
    final Statistics candidates = new Statistics();
    for (final RankResults rr : rankResultsArray) {
      actual.add(rr.zPosition.length);
      candidates.add(rr.spots.length);
    }
    add(sb, actual.getMean());
    add(sb, actual.getStandardDeviation());
    add(sb, candidates.getMean());
    add(sb, candidates.getStandardDeviation());

    final String resultPrefix = sb.toString();

    // ---
    // TODO: Further explore pre-filtering of spot candidates.

    // Add good label to spot candidates and have the benchmark spot filter respect this before
    // applying the fail count limit.

    // Correlation between intensity and SNR ...

    // SNR is very good at low density
    // SNR fails at high density. The SNR estimate is probably wrong for high intensity spots.

    // Triangle is very good when there are a large number of good spots in a region of the image
    // (e.g. a mask is used).
    // Triangle is poor when there are few good spots in an image.

    // Perhaps we can estimate the density of the spots and choose the correct thresholding method?

    // ---

    // Do a full benchmark through different Spot SNR, image sizes, densities and mask structures
    // and see if there are patterns for a good threshold method.

    // ---

    // Allow using the fitted results from benchmark spot fit. Will it make a difference if we fit
    // the candidates (some will fail if weak).
    // Can this be done by allowing the user to select the input (spot candidates or fitted
    // positions)?

    // Perhaps I need to produce a precision estimate for all simulated spots and then only use
    // those that achieve a certain precision, i.e. are reasonably in focus. Can this be done?
    // Does the image PSF have a width estimate for the entire stack?

    // Perhaps I should filter, fit and then filter all spots using no fail count. These then become
    // the spots to work with for creating a smart fail count filter.

    // ---

    // Pre-compute the results and have optional sort
    final ArrayList<ScoredResult> list = new ArrayList<>(methodNames.length);

    for (int i = 0; i < methodNames.length; i++) {
      tp = 0;
      fp = 0;
      double tn = 0;
      int itp = 0;
      int ifp = 0;
      int itn = 0;
      final Statistics s = new Statistics();
      long time = 0;

      for (final RankResults rr : rankResultsArray) {
        final RankResult r = rr.results.get(i);
        // Some results will not have a threshold
        if (!Float.isInfinite(r.threshold)) {
          s.add(r.threshold);
        }
        time += r.time;
        tp += r.fresult.getTruePositives();
        fp += r.fresult.getFalsePositives();
        tn += r.fresult.getTrueNegatives();
        itp += r.cresult.getTruePositives();
        ifp += r.cresult.getFalsePositives();
        itn += r.cresult.getTrueNegatives();
      }

      sb.setLength(0);
      sb.append(resultPrefix);
      add(sb, methodNames[i]);
      if (methodNames[i].startsWith("SNR")) {
        sb.append('\t');
      } else {
        add(sb, settings.compactBins);
      }
      add(sb, s.getMean());
      add(sb, s.getStandardDeviation());
      add(sb, TextUtils.nanosToString(time));

      // TP are all accepted candidates that can be matched to a spot
      // FP are all accepted candidates that cannot be matched to a spot
      // TN are all accepted candidates that cannot be matched to a spot
      // FN = The number of missed spots

      // Raw counts of match or no-match
      final FractionClassificationResult f1 =
          new FractionClassificationResult(itp, ifp, itn, simulationParameters.molecules - itp);
      final double s1 = addScores(sb, f1);

      // Fractional scoring
      final FractionClassificationResult f2 =
          new FractionClassificationResult(tp, fp, tn, simulationParameters.molecules - tp);
      final double s2 = addScores(sb, f2);

      // Store for sorting
      list.add(new ScoredResult(i, (settings.useFractionScores) ? s2 : s1, sb.toString()));
    }

    if (list.isEmpty()) {
      return;
    }

    Collections.sort(list, ScoredResult::compare);

    final TextWindow summaryTable = createTable();
    if (summaryTable.getTextPanel().getLineCount() > 0) {
      summaryTable.append("");
    }
    try (BufferedTextWindow tw = new BufferedTextWindow(summaryTable)) {
      tw.setIncrement(0);
      for (final ScoredResult r : list) {
        tw.append(r.result);
      }
    }

    if (settings.showOverlay) {
      final int bestMethod = list.get(0).index;
      final Overlay o = new Overlay();
      for (int j = 0; j < rankResultsArray.length; j++) {
        final int frame = frames[j];
        final RankResults rr = rankResultsArray[j];
        final RankResult r = rr.results.get(bestMethod);
        final int[] x1 = new int[r.good.length];
        final int[] y1 = new int[r.good.length];
        int c1 = 0;
        final int[] x2 = new int[r.good.length];
        final int[] y2 = new int[r.good.length];
        int c2 = 0;
        final int[] x3 = new int[r.good.length];
        final int[] y3 = new int[r.good.length];
        int c3 = 0;
        final int[] x4 = new int[r.good.length];
        final int[] y4 = new int[r.good.length];
        int c4 = 0;
        for (int i = 0; i < x1.length; i++) {
          if (r.good[i] == TP) {
            x1[c1] = rr.spots[i].spot.x;
            y1[c1] = rr.spots[i].spot.y;
            c1++;
          } else if (r.good[i] == FP) {
            x2[c2] = rr.spots[i].spot.x;
            y2[c2] = rr.spots[i].spot.y;
            c2++;
          } else if (r.good[i] == TN) {
            x3[c3] = rr.spots[i].spot.x;
            y3[c3] = rr.spots[i].spot.y;
            c3++;
          } else if (r.good[i] == FN) {
            x4[c4] = rr.spots[i].spot.x;
            y4[c4] = rr.spots[i].spot.y;
            c4++;
          }
        }
        addToOverlay(o, frame, x1, y1, c1, Color.green);
        addToOverlay(o, frame, x2, y2, c2, Color.red);
        if (IJ.debugMode) {
          addToOverlay(o, frame, x3, y3, c3, new Color(153, 255, 153)); // light green
        }
        addToOverlay(o, frame, x4, y4, c4, new Color(255, 153, 153)); // light red
      }

      imp.setOverlay(o);
    }
  }

  private static void addToOverlay(Overlay overlay, int frame, int[] x, int[] y, int count,
      Color color) {
    final PointRoi roi = new PointRoi(x, y, count);
    roi.setFillColor(color);
    roi.setStrokeColor(color);
    roi.setPosition(frame);
    roi.setShowLabels(false);
    overlay.add(roi);
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
    return ImageJUtils.refresh(SUMMARY_TABLE_REF,
        () -> new TextWindow(TITLE, createHeader(), "", 1000, 300));
  }

  private static String createHeader() {
    final StringBuilder sb = new StringBuilder(BenchmarkSpotFilter.getTablePrefix());
    sb.append("\tSpots\tcountPositive\tcountNegative\tfP\tfN\t% nP\t% nN\tcTotal\tcP\tcN\t"
        + "cfTotal\tcfP\tcfN\tSpot Av\tSpot SD\tCandidate Av\tCandidate SD\t"
        + "Method\tBins\tT Av\tT SD\tTime\t");

    addScoreColumns(sb, null);
    addScoreColumns(sb, "f ");

    return sb.toString();
  }

  private static void addScoreColumns(StringBuilder sb, String prefix) {
    for (int i = 1; i < Settings.SORT_METHODS.length; i++) {
      addScoreColumn(sb, prefix, Settings.SORT_METHODS[i]);
    }
  }

  private double addScores(StringBuilder sb, FractionClassificationResult result) {
    final double[] scores = new double[Settings.SORT_METHODS.length - 1];
    int index = 0;
    scores[index++] = result.getTruePositives();
    scores[index++] = result.getFalsePositives();
    scores[index++] = result.getTrueNegatives();
    scores[index++] = result.getFalseNegatives();
    scores[index++] = result.getPrecision();
    scores[index++] = result.getRecall();
    scores[index++] = result.getFScore(0.5);
    scores[index++] = result.getF1Score();
    scores[index++] = result.getFScore(2);
    scores[index++] = result.getJaccard();
    scores[index] = result.getMatthewsCorrelationCoefficient();
    for (final double s : scores) {
      add(sb, s);
    }
    return (settings.sortIndex != 0) ? scores[settings.sortIndex - 1] : 0;
  }

  private static void addScoreColumn(StringBuilder sb, String prefix, String name) {
    if (prefix != null) {
      sb.append(prefix);
    }
    sb.append(name).append('\t');
  }

  private double getSa() {
    return PsfCalculator.squarePixelAdjustment(simulationParameters.sd,
        simulationParameters.pixelPitch);
  }
}
