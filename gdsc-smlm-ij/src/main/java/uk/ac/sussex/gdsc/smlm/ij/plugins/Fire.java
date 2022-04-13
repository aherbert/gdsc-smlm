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

package uk.ac.sussex.gdsc.smlm.ij.plugins;

import com.google.common.util.concurrent.AtomicDouble;
import ij.IJ;
import ij.ImagePlus;
import ij.Macro;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.plugin.PlugIn;
import ij.plugin.frame.Recorder;
import ij.process.ImageProcessor;
import ij.process.LUT;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import java.awt.AWTEvent;
import java.awt.Checkbox;
import java.awt.Color;
import java.awt.Rectangle;
import java.awt.Scrollbar;
import java.awt.TextField;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.Rectangle2D;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.SimpleCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.optim.univariate.BracketFinder;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariateOptimizer;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.MathArrays;
import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.core.data.utils.ConversionException;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.NonBlockingExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper.LutColour;
import uk.ac.sussex.gdsc.core.logging.TrackProgress;
import uk.ac.sussex.gdsc.core.logging.TrackProgressAdapter;
import uk.ac.sussex.gdsc.core.utils.DoubleMedianWindow;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.concurrent.ConcurrencyUtils;
import uk.ac.sussex.gdsc.core.utils.rng.RandomUtils;
import uk.ac.sussex.gdsc.core.utils.rng.UniformRandomProviders;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageMode;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageSettings;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageSizeMode;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageType;
import uk.ac.sussex.gdsc.smlm.data.config.UnitHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.function.Erf;
import uk.ac.sussex.gdsc.smlm.ij.frc.Frc;
import uk.ac.sussex.gdsc.smlm.ij.frc.Frc.FourierMethod;
import uk.ac.sussex.gdsc.smlm.ij.frc.Frc.FrcCurve;
import uk.ac.sussex.gdsc.smlm.ij.frc.Frc.FrcCurveResult;
import uk.ac.sussex.gdsc.smlm.ij.frc.Frc.FrcFireResult;
import uk.ac.sussex.gdsc.smlm.ij.frc.Frc.SamplingMethod;
import uk.ac.sussex.gdsc.smlm.ij.frc.Frc.ThresholdMethod;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.results.ImageJImagePeakResults;
import uk.ac.sussex.gdsc.smlm.ij.results.ImagePeakResultsFactory;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.count.Counter;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.PrecisionResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.XyrResultProcedure;
import uk.ac.sussex.gdsc.smlm.utils.StdMath;

/**
 * Computes the Fourier Image Resolution of an image
 *
 * <p>Implements the FIRE (Fourier Image REsolution) method described in:<br> Niewenhuizen, et al
 * (2013). Measuring image resolution in optical nanoscopy. Nature Methods, 10, 557<br>
 * http://www.nature.com/nmeth/journal/v10/n6/full/nmeth.2448.html
 *
 * <p>A second plugin allows estimation of the spurious correlation component contributed by the
 * same molecule being present in both super-resolution images due to splitting of repeat
 * localisations. The correction Q factor is the number of times a molecule is repeat localised
 * (i.e. average blinks per molecule). This code was developed using the Matlab examples provided by
 * Bernd Reiger.
 */
public class Fire implements PlugIn {
  private String pluginTitle = "Fourier Image REsolution (FIRE)";
  private static final String KEY_MEAN = "mean_estimate";
  private static final String KEY_SIGMA = "sigma_estimate";
  private static final String KEY_Q = "q_estimate";

  /**
   * Hold a copy of the ImageJ preferences for the number of threads. This is atomically updated
   * using synchronised methods.
   */
  private static int imagejNThreads = Prefs.getThreads();
  /**
   * Hold a copy of the user choice for the number of threads. This is atomically updated using
   * synchronised methods.
   */
  private static int lastNThreads = imagejNThreads;

  private boolean extraOptions;
  private boolean myUseSignal;
  private Rectangle roiBounds;
  private int roiImageWidth;
  private int roiImageHeight;

  // Stored in initialisation

  /** The results. */
  MemoryPeakResults results;
  private MemoryPeakResults results2;
  private Rectangle2D dataBounds;
  /** The spatial units. */
  String units;
  private double nmPerUnit = 1;

  // Stored in setCorrectionParameters
  private double correctionQValue;
  private double correctionMean;
  private double correctionSigma;

  private TrackProgress progress = new ParallelTrackProgress(1);

  private int numberOfThreads;

  /** The plugin settings. */
  Settings settings;

  private FourierMethod fourierMethod;
  private SamplingMethod samplingMethod;
  /** The threshold method. */
  ThresholdMethod thresholdMethod;
  private PrecisionMethod precisionMethod;

  private PrecisionResultProcedure pp;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    static final String[] SCALE_ITEMS;
    static final int[] SCALE_VALUES = {0, 1, 2, 4, 8, 16, 32, 64, 128};
    static final String[] IMAGE_SIZE_ITEMS;
    static final int[] IMAGE_SIZE_VALUES;

    static {
      SCALE_ITEMS = new String[SCALE_VALUES.length];
      SCALE_ITEMS[0] = "Auto";
      for (int i = 1; i < SCALE_VALUES.length; i++) {
        SCALE_ITEMS[i] = Integer.toString(SCALE_VALUES[i]);
      }

      // Create size for Fourier transforms. Must be power of 2.
      final int[] imageSizeValues = new int[32];
      final String[] imageSizeItems = new String[imageSizeValues.length];
      int size = 512; // Start at a reasonable size. Too small does not work.
      int count = 0;
      while (size <= 16384) {
        // Assumes the scaling works correctly for exact power of 2
        // (no rounding errors creating images 1 pixel too large)
        imageSizeValues[count] = size;
        imageSizeItems[count] = Integer.toString(size);
        size *= 2;
        count++;
      }
      IMAGE_SIZE_VALUES = Arrays.copyOf(imageSizeValues, count);
      IMAGE_SIZE_ITEMS = Arrays.copyOf(imageSizeItems, count);
    }

    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> INSTANCE = new AtomicReference<>(new Settings());

    String inputOption;
    String inputOption2;

    int repeats;
    boolean useSignal;
    int maxPerBin; // 5 in the Niewenhuizen paper
    boolean randomSplit;
    int blockSize;
    int imageScaleIndex;
    int imageSizeIndex;

    // The Q value and the mean and sigma for spurious correlation correction
    boolean spuriousCorrelationCorrection;
    double qvalue;
    double mean;
    double sigma;

    double perimeterSamplingFactor;
    int fourierMethodIndex;
    int samplingMethodIndex;
    int thresholdMethodIndex;
    boolean showFrcCurve;
    boolean showFrcCurveRepeats;
    boolean showFrcTimeEvolution;
    int precisionMethodIndex;
    boolean sampleDecay;
    boolean loessSmoothing;
    boolean fitPrecision;
    double minQ;
    double maxQ;

    boolean chooseRoi;
    String roiImage;

    Settings() {
      // Set defaults
      inputOption = "";
      inputOption2 = "";
      repeats = 1;
      randomSplit = true;
      blockSize = 50;
      imageSizeIndex = Arrays.binarySearch(IMAGE_SIZE_VALUES, 2047);
      perimeterSamplingFactor = 1;
      fourierMethodIndex = FourierMethod.JTRANSFORMS.ordinal();
      samplingMethodIndex = SamplingMethod.RADIAL_SUM.ordinal();
      thresholdMethodIndex = ThresholdMethod.FIXED_1_OVER_7.ordinal();
      showFrcCurve = true;
      precisionMethodIndex = PrecisionMethod.CALCULATE.ordinal();
      minQ = 0.2;
      maxQ = 0.45;
      roiImage = "";
    }

    Settings(Settings source) {
      inputOption = source.inputOption;
      inputOption2 = source.inputOption2;
      repeats = source.repeats;
      useSignal = source.useSignal;
      maxPerBin = source.maxPerBin; // 5 in the Niewenhuizen paper
      randomSplit = source.randomSplit;
      blockSize = source.blockSize;
      imageScaleIndex = source.imageScaleIndex;
      imageSizeIndex = source.imageSizeIndex;
      spuriousCorrelationCorrection = source.spuriousCorrelationCorrection;
      qvalue = source.qvalue;
      mean = source.mean;
      sigma = source.sigma;
      perimeterSamplingFactor = source.perimeterSamplingFactor;
      fourierMethodIndex = source.fourierMethodIndex;
      samplingMethodIndex = source.samplingMethodIndex;
      thresholdMethodIndex = source.thresholdMethodIndex;
      showFrcCurve = source.showFrcCurve;
      showFrcCurveRepeats = source.showFrcCurveRepeats;
      showFrcTimeEvolution = source.showFrcTimeEvolution;
      precisionMethodIndex = source.precisionMethodIndex;
      sampleDecay = source.sampleDecay;
      loessSmoothing = source.loessSmoothing;
      fitPrecision = source.fitPrecision;
      minQ = source.minQ;
      maxQ = source.maxQ;
      chooseRoi = source.chooseRoi;
      roiImage = source.roiImage;
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
   * Specify the method to use to determine the parameters for the distribution of the localisation
   * precision (assumed to be Gaussian).
   */
  private enum PrecisionMethod {
    //@formatter:off
    /**
     * Use a fixed value for the precision distribution mean and standard deviation
     */
    FIXED{ @Override
    public String getName() { return "Fixed"; }},
    /**
     * Use the precision value that is stored in the results, e.g. from loaded data.
     * The values can then be used to fit the entire distribution using a Gaussian or
     * sampled to construct a decay curve from which the parameters are estimated.
     */
    STORED{ @Override
    public String getName() { return "Stored"; }},
    /**
     * Calculate the precision of each localisation using the formula of Mortensen.
     * The values can then be used to fit the entire distribution using a Gaussian or
     * sampled to construct a decay curve from which the parameters are estimated.
     */
    CALCULATE{ @Override
    public String getName() { return "Calculate"; }};
    //@formatter:on

    @Override
    public String toString() {
      return getName();
    }

    /**
     * Gets the name.
     *
     * @return the name
     */
    public abstract String getName();
  }

  /**
   * Store images for FIRE analysis.
   */
  public static class FireImages {
    /** The first super-resolution image. */
    final ImageProcessor ip1;
    /** The second super-resolution image. */
    final ImageProcessor ip2;
    /** The nm per pixel in the super-resolution images. */
    final double nmPerPixel;

    /**
     * Instantiates a new fire images.
     *
     * @param ip1 the first image
     * @param ip2 the second image
     * @param nmPerPixel the nm per pixel
     */
    FireImages(ImageProcessor ip1, ImageProcessor ip2, double nmPerPixel) {
      this.ip1 = ip1;
      this.ip2 = ip2;
      this.nmPerPixel = nmPerPixel;
    }
  }

  /**
   * Contains the Fourier Image REsolution (FIRE) result.
   */
  public static class FireResult {
    /** The fire number (in nm). */
    final double fireNumber;

    /** The correlation at the given resolution. */
    final double correlation;

    /** The FRC curve used to compute the resolution. */
    final FrcCurve frcCurve;

    /** The original correlation curve, i.e. the raw curve before smoothing. */
    final double[] originalCorrelationCurve;

    /**
     * Instantiates a new fire result.
     *
     * @param fireNumber the fire number
     * @param correlation the correlation
     * @param frcCurve the frc curve
     * @param originalCorrelationCurve the original correlation curve
     */
    FireResult(double fireNumber, double correlation, FrcCurve frcCurve,
        double[] originalCorrelationCurve) {
      this.fireNumber = fireNumber;
      this.correlation = correlation;
      this.frcCurve = frcCurve;
      this.originalCorrelationCurve = originalCorrelationCurve;
    }

    /**
     * Gets the nm per pixel for the super-resolution images used to construct the FRC curve.
     *
     * @return the nm per pixel
     */
    double getNmPerPixel() {
      return frcCurve.nmPerPixel;
    }
  }

  /**
   * Worker to compute the FIRE result.
   */
  private class FireWorker implements Runnable {
    final double fourierImageScale;
    final int imageSize;

    String name;
    FireResult result;
    Plot plot;
    /**
     * Flag to denote that an out-of-memory error occurred. This is probably due to using too many
     * threads to compute large Fourier transforms.
     */
    boolean oom;

    FireWorker(int id, double fourierImageScale, int imageSize) {
      this.fourierImageScale = fourierImageScale;
      this.imageSize = imageSize;
      name = results.getName() + " [" + id + "]";
    }

    @Override
    public void run() {
      try {
        result = calculateFireNumber(fourierMethod, samplingMethod, thresholdMethod,
            fourierImageScale, imageSize);
        if (settings.showFrcCurve) {
          plot = createFrcCurve(name, result, thresholdMethod);
          if (settings.showFrcCurveRepeats) {
            // Do this on the thread
            plot.draw();
          }
        }
      } catch (final OutOfMemoryError ex) {
        oom = true;
      }
    }
  }

  /**
   * Dumb implementation of the track progress interface for parallel threads. Uses simple
   * synchronisation to increment total progress.
   */
  private static class ParallelTrackProgress extends TrackProgressAdapter {
    final AtomicDouble done = new AtomicDouble();
    final int total;

    ParallelTrackProgress(int repeats) {
      total = repeats;
    }

    @Override
    public void incrementProgress(double fraction) {
      // Avoid synchronisation for nothing
      if (fraction == 0) {
        return;
      }
      IJ.showProgress(done.addAndGet(fraction) / this.total);
    }
  }

  @Override
  public void run(String arg) {
    extraOptions = ImageJUtils.isExtraOptions();
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    // Require some fit results and selected regions
    final int size = MemoryPeakResults.countMemorySize();
    if (size == 0) {
      IJ.error(pluginTitle, "There are no fitting results in memory");
      return;
    }

    settings = Settings.load();
    settings.save();

    if ("q".equals(arg)) {
      pluginTitle += " Q estimation";
      runQEstimation();
      return;
    }

    IJ.showStatus(pluginTitle + " ...");

    if (!showInputDialog()) {
      return;
    }

    MemoryPeakResults inputResults1 =
        ResultsManager.loadInputResults(settings.inputOption, false, null, null);
    if (MemoryPeakResults.isEmpty(inputResults1)) {
      IJ.error(pluginTitle, "No results could be loaded");
      return;
    }
    inputResults1 = cropToRoi(inputResults1);
    if (inputResults1.size() <= 1) {
      IJ.error(pluginTitle, "No results within the crop region");
      return;
    }

    MemoryPeakResults inputResults2 =
        ResultsManager.loadInputResults(settings.inputOption2, false, null, null);
    if (inputResults2 != null) {
      inputResults2 = cropToRoi(inputResults2);
      if (inputResults2.size() <= 1) {
        IJ.error(pluginTitle, "No results2 within the crop region");
        return;
      }
    }

    initialise(inputResults1, inputResults2);

    if (!showDialog()) {
      return;
    }

    final long start = System.currentTimeMillis();

    // Compute FIRE

    String name = inputResults1.getName();
    final double fourierImageScale = Settings.SCALE_VALUES[settings.imageScaleIndex];
    final int imageSize = Settings.IMAGE_SIZE_VALUES[settings.imageSizeIndex];

    if (this.results2 == null) {
      FireResult result = null;

      final int repeats = (settings.randomSplit) ? Math.max(1, settings.repeats) : 1;
      setProgress(repeats);
      if (repeats == 1) {
        result = calculateFireNumber(fourierMethod, samplingMethod, thresholdMethod,
            fourierImageScale, imageSize);

        if (result != null) {
          logResult(name, result);

          if (settings.showFrcCurve) {
            showFrcCurve(name, result, thresholdMethod);
          }
        }
      } else {
        // Multi-thread this ...
        final int nThreads = MathUtils.min(repeats, getThreads());
        final ExecutorService executor = Executors.newFixedThreadPool(nThreads);
        final LocalList<Future<?>> futures = new LocalList<>(repeats);
        final LocalList<FireWorker> workers = new LocalList<>(repeats);
        IJ.showProgress(0);
        IJ.showStatus(pluginTitle + " computing ...");
        for (int i = 1; i <= repeats; i++) {
          final FireWorker w = new FireWorker(i, fourierImageScale, imageSize);
          workers.add(w);
          futures.add(executor.submit(w));
        }

        // Wait for all to finish
        executor.shutdown();
        ConcurrencyUtils.waitForCompletionUnchecked(futures);
        IJ.showProgress(1);

        // Show a combined FRC curve plot of all the smoothed curves if we have multiples.
        final LUT valuesLut = LutHelper.createLut(LutColour.FIRE_GLOW);
        final LutHelper.DefaultLutMapper mapper = new LutHelper.DefaultLutMapper(0, repeats);
        final FrcCurvePlot curve = new FrcCurvePlot();

        final Statistics stats = new Statistics();
        final WindowOrganiser wo = new WindowOrganiser();
        boolean oom = false;
        for (int i = 0; i < repeats; i++) {
          final FireWorker w = workers.get(i);
          if (w.oom) {
            oom = true;
          }
          if (w.result == null) {
            continue;
          }
          result = w.result;
          if (!Double.isNaN(result.fireNumber)) {
            stats.add(result.fireNumber);
          }

          if (settings.showFrcCurveRepeats) {
            // Output each FRC curve using a suffix.
            logResult(w.name, result);
            wo.add(ImageJUtils.display(w.plot.getTitle(), w.plot));
          }
          if (settings.showFrcCurve) {
            final int index = mapper.map(i + 1);
            curve.add(name, result, thresholdMethod, LutHelper.getColour(valuesLut, index),
                Color.blue, null);
          }
        }

        if (result != null) {
          wo.cascade();
          final double mean = stats.getMean();
          logResult(name, result, mean, stats);
          if (settings.showFrcCurve) {
            curve.addResolution(mean);
            final Plot plot = curve.getPlot();
            ImageJUtils.display(plot.getTitle(), plot);
          }
        }

        if (oom) {
          //@formatter:off
          IJ.error(pluginTitle,
              "ERROR - Parallel computation out-of-memory.\n \n" +
          TextUtils.wrap("The number of results will be reduced. " +
                  "Please reduce the size of the Fourier image " +
                  "or change the number of threads " +
                  "using the extra options (hold down the 'Shift' " +
                  "key when running the plugin).",
                  80));
          //@formatter:on
        }
      }

      // Only do this once
      if (settings.showFrcTimeEvolution && result != null && !Double.isNaN(result.fireNumber)) {
        showFrcTimeEvolution(name, result.fireNumber, thresholdMethod,
            nmPerUnit / result.getNmPerPixel(), imageSize);
      }
    } else {
      name += " vs " + this.results2.getName();

      final FireResult result = calculateFireNumber(fourierMethod, samplingMethod, thresholdMethod,
          fourierImageScale, imageSize);

      if (result != null) {
        logResult(name, result);

        if (settings.showFrcCurve) {
          showFrcCurve(name, result, thresholdMethod);
        }
      }
    }

    IJ.showStatus(pluginTitle + " complete : "
        + TextUtils.millisToString(System.currentTimeMillis() - start));
  }

  private void logResult(String name, FireResult result) {
    IJ.log(String.format("%s : FIRE number = %s %s (Fourier scale = %s)", name,
        MathUtils.rounded(result.fireNumber, 4), units,
        MathUtils.rounded(nmPerUnit / result.getNmPerPixel(), 3)));
    if (Double.isNaN(result.fireNumber)) {
      ImageJUtils.log(
          "%s Warning: NaN result possible if the resolution is below the pixel size of the"
              + " input Fourier image (%s %s).",
          pluginTitle, MathUtils.rounded(result.getNmPerPixel()), units);
    }
  }

  private void logResult(String name, FireResult result, double mean, Statistics stats) {
    IJ.log(String.format("%s : FIRE number = %s +/- %s %s [95%% CI, n=%d] (Fourier scale = %s)",
        name, MathUtils.rounded(mean, 4), MathUtils.rounded(stats.getConfidenceInterval(0.95), 4),
        units, stats.getN(), MathUtils.rounded(nmPerUnit / result.getNmPerPixel(), 3)));
  }

  private MemoryPeakResults cropToRoi(MemoryPeakResults results) {
    if (roiBounds == null) {
      return results;
    }

    // Adjust bounds relative to input results image
    final Rectangle bounds = results.getBounds(true);
    final double xscale = (double) roiImageWidth / bounds.width;
    final double yscale = (double) roiImageHeight / bounds.height;

    final float minX = (float) (bounds.x + roiBounds.x / xscale);
    final float maxX = (float) (minX + roiBounds.width / xscale);
    final float minY = (float) (bounds.y + (roiBounds.y / yscale));
    final float maxY = (float) (minY + roiBounds.height / yscale);

    // Create a new set of results within the bounds
    final MemoryPeakResults newResults = new MemoryPeakResults();
    newResults.begin();
    results.forEach(DistanceUnit.PIXEL, (XyrResultProcedure) (x, y, result) -> {
      if (x < minX || x > maxX || y < minY || y > maxY) {
        return;
      }
      newResults.add(result);
    });
    newResults.end();
    newResults.copySettings(results);
    newResults.setBounds(new Rectangle((int) minX, (int) minY, (int) Math.ceil(maxX - minX),
        (int) Math.ceil(maxY - minY)));
    return newResults;
  }

  private boolean showInputDialog() {
    ExtendedGenericDialog gd = new ExtendedGenericDialog(pluginTitle);
    gd.addMessage("Compute the resolution using Fourier Ring Correlation");
    gd.addHelp(HelpUrls.getUrl("fourier-image-resolution"));

    // Build a list of all images with a region ROI
    final List<String> titles = new LinkedList<>();
    for (final int imageId : ImageJUtils.getIdList()) {
      final ImagePlus imp = WindowManager.getImage(imageId);
      if (imp != null && imp.getRoi() != null && imp.getRoi().isArea()) {
        titles.add(imp.getTitle());
      }
    }

    ResultsManager.addInput(gd, settings.inputOption, InputSource.MEMORY);
    ResultsManager.addInput(gd, "Input2", settings.inputOption2, InputSource.NONE,
        InputSource.MEMORY);

    if (!titles.isEmpty()) {
      gd.addCheckbox((titles.size() == 1) ? "Use_ROI" : "Choose_ROI", settings.chooseRoi);
    }

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.inputOption = ResultsManager.getInputSource(gd);
    settings.inputOption2 = ResultsManager.getInputSource(gd);

    if (!titles.isEmpty()) {
      settings.chooseRoi = gd.getNextBoolean();
    }

    if (!titles.isEmpty() && settings.chooseRoi) {
      if (titles.size() == 1) {
        settings.roiImage = titles.get(0);
        Recorder.recordOption("Image", settings.roiImage);
      } else {
        final String[] items = titles.toArray(new String[0]);
        gd = new ExtendedGenericDialog(pluginTitle);
        gd.addMessage("Select the source image for the ROI");
        gd.addChoice("Image", items, settings.roiImage);
        gd.showDialog();
        if (gd.wasCanceled()) {
          return false;
        }
        settings.roiImage = gd.getNextChoice();
      }
      final ImagePlus imp = WindowManager.getImage(settings.roiImage);

      roiBounds = imp.getRoi().getBounds();
      roiImageWidth = imp.getWidth();
      roiImageHeight = imp.getHeight();
    } else {
      roiBounds = null;
    }

    return true;
  }

  private boolean showDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(pluginTitle);
    gd.addMessage("Compute the resolution using Fourier Ring Correlation");
    gd.addHelp(HelpUrls.getUrl("fourier-image-resolution"));

    final boolean single = results2 == null;

    gd.addMessage("Image construction options:");
    gd.addChoice("Image_scale", Settings.SCALE_ITEMS, settings.imageScaleIndex);
    gd.addChoice("Auto_image_size", Settings.IMAGE_SIZE_ITEMS, settings.imageSizeIndex);
    if (extraOptions) {
      gd.addCheckbox("Use_signal (if present)", settings.useSignal);
    }
    gd.addNumericField("Max_per_bin", settings.maxPerBin, 0);

    gd.addMessage("Fourier options:");
    final String[] fourierMethodNames =
        SettingsManager.getNames((Object[]) Frc.FourierMethod.values());
    gd.addChoice("Fourier_method", fourierMethodNames,
        fourierMethodNames[settings.fourierMethodIndex]);
    final String[] samplingMethodNames =
        SettingsManager.getNames((Object[]) Frc.SamplingMethod.values());
    gd.addChoice("Sampling_method", samplingMethodNames,
        samplingMethodNames[settings.samplingMethodIndex]);
    gd.addSlider("Sampling_factor", 0.2, 4, settings.perimeterSamplingFactor);

    gd.addMessage("FIRE options:");
    final String[] thresholdMethodNames =
        SettingsManager.getNames((Object[]) Frc.ThresholdMethod.values());
    gd.addChoice("Threshold_method", thresholdMethodNames,
        thresholdMethodNames[settings.thresholdMethodIndex]);
    gd.addCheckbox("Show_FRC_curve", settings.showFrcCurve);

    if (single) {
      gd.addMessage("For single datasets:");
      gd.addNumericField("Block_size", settings.blockSize, 0);
      gd.addCheckbox("Random_split", settings.randomSplit);
      gd.addNumericField("Repeats", settings.repeats, 0);
      gd.addCheckbox("Show_FRC_curve_repeats", settings.showFrcCurveRepeats);
      gd.addCheckbox("Show_FRC_time_evolution", settings.showFrcTimeEvolution);
      gd.addCheckbox("Spurious correlation correction", settings.spuriousCorrelationCorrection);
      gd.addNumericField("Q-value", settings.qvalue, 3);
      gd.addNumericField("Precision_Mean", settings.mean, 2, 6, "nm");
      gd.addNumericField("Precision_Sigma", settings.sigma, 2, 6, "nm");
      if (extraOptions) {
        gd.addNumericField("Threads", getLastNThreads(), 0);
      }
    }

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.imageScaleIndex = gd.getNextChoiceIndex();
    settings.imageSizeIndex = gd.getNextChoiceIndex();
    if (extraOptions) {
      myUseSignal = settings.useSignal = gd.getNextBoolean();
    }
    settings.maxPerBin = Math.abs((int) gd.getNextNumber());

    settings.fourierMethodIndex = gd.getNextChoiceIndex();
    fourierMethod = FourierMethod.values()[settings.fourierMethodIndex];
    settings.samplingMethodIndex = gd.getNextChoiceIndex();
    samplingMethod = SamplingMethod.values()[settings.samplingMethodIndex];
    settings.perimeterSamplingFactor = gd.getNextNumber();

    settings.thresholdMethodIndex = gd.getNextChoiceIndex();
    thresholdMethod = Frc.ThresholdMethod.values()[settings.thresholdMethodIndex];
    settings.showFrcCurve = gd.getNextBoolean();

    if (single) {
      settings.blockSize = Math.max(1, (int) gd.getNextNumber());
      settings.randomSplit = gd.getNextBoolean();
      settings.repeats = Math.max(1, (int) gd.getNextNumber());
      settings.showFrcCurveRepeats = gd.getNextBoolean();
      settings.showFrcTimeEvolution = gd.getNextBoolean();
      settings.spuriousCorrelationCorrection = gd.getNextBoolean();
      settings.qvalue = Math.abs(gd.getNextNumber());
      settings.mean = Math.abs(gd.getNextNumber());
      settings.sigma = Math.abs(gd.getNextNumber());
      if (extraOptions) {
        setThreads((int) gd.getNextNumber());
        setLastNThreads(this.numberOfThreads);
      }
    }

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Perimeter sampling factor", settings.perimeterSamplingFactor);
      if (single && settings.spuriousCorrelationCorrection) {
        ParameterUtils.isAboveZero("Q-value", settings.qvalue);
        ParameterUtils.isAboveZero("Precision Mean", settings.mean);
        ParameterUtils.isAboveZero("Precision Sigma", settings.sigma);
        // Set these for use in FIRE computation
        setCorrectionParameters(settings.qvalue, settings.mean, settings.sigma);
      }
    } catch (final IllegalArgumentException ex) {
      IJ.error(pluginTitle, ex.getMessage());
      return false;
    }

    return true;
  }

  /**
   * Initialise this instance with localisation results for the FIRE computation.
   *
   * @param results the results
   * @param results2 the second set of results (can be null)
   */
  public void initialise(MemoryPeakResults results, MemoryPeakResults results2) {
    this.results = verify(results);
    this.results2 = verify(results2);

    if (this.results == null) {
      return;
    }

    nmPerUnit = 1;
    DistanceUnit unit = null;
    units = "unknown";

    final CalibrationReader cal = results.getCalibrationReader();
    if (cal == null) {
      IJ.log(pluginTitle + " Warning: No calibration exists for primary results");
    } else {
      try {
        nmPerUnit = cal.getDistanceConverter(DistanceUnit.NM).convert(1);
        units = UnitHelper.getShortName(DistanceUnit.NM);
        unit = DistanceUnit.NM;
      } catch (final ConversionException ex) {
        IJ.log(pluginTitle + " Warning: Ignoring invalid distance calibration for primary results");
      }
    }

    // Calibration must match between datasets
    if (this.results2 != null) {
      CalibrationReader cal2 = results.getCalibrationReader();
      if (unit == null) {
        if (cal2 != null) {
          IJ.log(pluginTitle
              + " Warning: Ignoring calibration for secondary results since no calibration"
              + " exists for primary results");
        }
      } else {
        // The calibration must match
        try {
          // Try to create a converter and check it is the same conversion
          if (cal2 != null && cal2.getDistanceConverter(DistanceUnit.NM).convert(1) != nmPerUnit) {
            // Set to null to mark invalid
            cal2 = null;
          }
        } catch (final ConversionException ex) {
          // Set to null to mark invalid
          cal2 = null;
        } finally {
          if (cal2 == null) {
            this.results = null;
            IJ.error(pluginTitle,
                "Error: Calibration between the two input datasets does not match");
          }
        }
      }
    }

    // Use the float data bounds. This prevents problems if the data is far from the origin.
    dataBounds = results.getDataBounds(null);

    if (this.results2 != null) {
      final Rectangle2D dataBounds2 = results.getDataBounds(null);
      dataBounds = dataBounds.createUnion(dataBounds2);
    }
  }

  /**
   * Sets the correction parameters for spurious correlation correction. Only relevant for single
   * images.
   *
   * @param qvalue the q value
   * @param mean the mean of the localisation precision
   * @param sigma the standard deviation of the localisation precision
   */
  public void setCorrectionParameters(double qvalue, double mean, double sigma) {
    if (qvalue > 0 && mean > 0 && sigma > 0) {
      correctionQValue = qvalue;
      correctionMean = mean;
      correctionSigma = sigma;
    } else {
      correctionQValue = correctionMean = correctionSigma = 0;
    }
  }

  /**
   * Shallow copy this instance so skipping initialisation. Only variables required for the FIRE
   * calculation are copied.
   *
   * @return the new FIRE instance
   */
  private Fire copy() {
    final Fire f = new Fire();
    f.results = results;
    f.results2 = results2;
    f.nmPerUnit = nmPerUnit;
    f.units = units;
    f.dataBounds = dataBounds;
    f.correctionQValue = correctionQValue;
    f.correctionMean = correctionMean;
    f.correctionSigma = correctionSigma;
    f.settings = settings;
    return f;
  }

  /**
   * Verify the results can be used for FIRE. Results are sorted in time order if the block size is
   * above 1.
   *
   * @param results the results
   * @return the memory peak results
   */
  private MemoryPeakResults verify(MemoryPeakResults results) {
    if (results == null || results.size() <= 1) {
      return null;
    }
    if (settings.blockSize > 1) {
      // Results must be in time order when processing blocks
      results.sort();
    }
    return results;
  }

  /**
   * Creates the images to use for the FIRE calculation. This must be called after
   * {@link #initialise(MemoryPeakResults, MemoryPeakResults)}.
   *
   * @param fourierImageScale the fourier image scale (set to zero to auto compute)
   * @param imageSize the image size
   * @return the fire images
   */
  public FireImages createImages(double fourierImageScale, int imageSize) {
    return createImages(fourierImageScale, imageSize, myUseSignal);
  }

  /**
   * Get the signal from the result.
   */
  private interface SignalProvider {
    /**
     * Gets the signal.
     *
     * @param result the result
     * @return the signal
     */
    float getSignal(PeakResult result);
  }

  /**
   * Creates the images to use for the FIRE calculation. This must be called after
   * {@link #initialise(MemoryPeakResults, MemoryPeakResults)}.
   *
   * @param fourierImageScale the fourier image scale (set to zero to auto compute)
   * @param imageSize the image size
   * @param useSignal Use the localisation signal to weight the intensity. The default uses a value
   *        of 1 per localisation.
   * @return the fire images
   */
  public FireImages createImages(double fourierImageScale, int imageSize, boolean useSignal) {
    if (results == null) {
      return null;
    }

    final boolean fixedSignal = useSignal && results.hasIntensity();
    final SignalProvider signalProvider = fixedSignal ? x -> x.getIntensity() : x -> 1f;

    // Draw images using the existing IJ routines.
    final Rectangle bounds = new Rectangle((int) Math.ceil(dataBounds.getWidth()),
        (int) Math.ceil(dataBounds.getHeight()));

    final ResultsImageSettings.Builder builder =
        ResultsImageSettings.newBuilder().setImageType(ResultsImageType.DRAW_NONE).setWeighted(true)
            .setEqualised(false).setImageMode(ResultsImageMode.IMAGE_ADD);
    if (fourierImageScale > 0) {
      builder.setImageSizeMode(ResultsImageSizeMode.SCALED);
      builder.setScale(fourierImageScale);
    } else {
      builder.setImageSizeMode(ResultsImageSizeMode.IMAGE_SIZE);
      builder.setImageSize(imageSize);
    }

    ImageJImagePeakResults image1 = createPeakResultsImage(bounds, builder, "IP1");
    ImageJImagePeakResults image2 = createPeakResultsImage(bounds, builder, "IP2");

    final float minx = (float) dataBounds.getX();
    final float miny = (float) dataBounds.getY();

    if (this.results2 == null) {
      // Block sampling.
      // Ensure we have at least 2 even sized blocks.
      int blockSize = Math.min(results.size() / 2, Math.max(1, settings.blockSize));
      int nblocks = (int) Math.ceil((double) results.size() / blockSize);
      while (nblocks <= 1 && blockSize > 1) {
        blockSize /= 2;
        nblocks = (int) Math.ceil((double) results.size() / blockSize);
      }
      if (nblocks <= 1) {
        // This should not happen since the results should contain at least 2 localisations
        return null;
      }
      if (blockSize != settings.blockSize) {
        IJ.log(pluginTitle + " Warning: Changed block size to " + blockSize);
      }

      final Counter i = new Counter();
      final Counter block = new Counter();
      final int finalBlockSize = blockSize;
      final PeakResult[][] blocks = new PeakResult[nblocks][blockSize];
      results.forEach((PeakResultProcedure) result -> {
        if (i.getCount() == finalBlockSize) {
          block.increment();
          i.reset();
        }
        blocks[block.getCount()][i.getAndIncrement()] = result;
      });
      // Truncate last block
      blocks[block.getCount()] = Arrays.copyOf(blocks[block.getCount()], i.getCount());

      final int[] indices = SimpleArrayUtils.natural(block.getCount() + 1);
      if (settings.randomSplit) {
        MathArrays.shuffle(indices);
      }

      for (final int index : indices) {
        // Split alternating so just rotate
        final ImageJImagePeakResults image = image1;
        image1 = image2;
        image2 = image;
        for (final PeakResult p : blocks[index]) {
          final float x = p.getXPosition() - minx;
          final float y = p.getYPosition() - miny;
          image.add(x, y, signalProvider.getSignal(p));
        }
      }
    } else {
      // Two image comparison
      final ImageJImagePeakResults i1 = image1;
      results.forEach((PeakResultProcedure) result -> {
        final float x = result.getXPosition() - minx;
        final float y = result.getYPosition() - miny;
        i1.add(x, y, signalProvider.getSignal(result));
      });
      final ImageJImagePeakResults i2 = image2;
      results2.forEach((PeakResultProcedure) result -> {
        final float x = result.getXPosition() - minx;
        final float y = result.getYPosition() - miny;
        i2.add(x, y, signalProvider.getSignal(result));
      });
    }

    image1.end();
    final ImageProcessor ip1 = image1.getImagePlus().getProcessor();

    image2.end();
    final ImageProcessor ip2 = image2.getImagePlus().getProcessor();

    if (settings.maxPerBin > 0 && fixedSignal) {
      // We can eliminate over-sampled pixels
      for (int i = ip1.getPixelCount(); i-- > 0;) {
        if (ip1.getf(i) > settings.maxPerBin) {
          ip1.setf(i, settings.maxPerBin);
        }
        if (ip2.getf(i) > settings.maxPerBin) {
          ip2.setf(i, settings.maxPerBin);
        }
      }
    }

    return new FireImages(ip1, ip2, nmPerUnit / image1.getScale());
  }

  /**
   * Creates the peak results image.
   *
   * @param bounds the bounds
   * @param builder the builder
   * @param title the title
   * @return the peak results image
   */
  private static ImageJImagePeakResults createPeakResultsImage(final Rectangle bounds,
      final ResultsImageSettings.Builder builder, String title) {
    final ImageJImagePeakResults image1 =
        ImagePeakResultsFactory.createPeakResultsImage(builder, title, bounds, 1);
    image1.setDisplayImage(false);
    image1.setUncalibrated(true);
    image1.begin();
    return image1;
  }

  /**
   * Encapsulate plotting the FRC curve to allow multiple curves to be plotted together.
   */
  private class FrcCurvePlot {
    double[] xValues;
    double[] threshold;
    Plot plot;

    void add(String name, FireResult result, ThresholdMethod thresholdMethod, Color colorValues,
        Color colorThreshold, Color colorNoSmooth) {
      final FrcCurve frcCurve = result.frcCurve;

      final double[] yValues = new double[frcCurve.getSize()];

      if (plot == null) {
        final String title = name + " FRC Curve";
        plot = new Plot(title, String.format("Spatial Frequency (%s^-1)", units), "FRC");

        xValues = new double[frcCurve.getSize()];
        final double l = frcCurve.fieldOfView;
        final double conversion = 1.0 / (l * result.getNmPerPixel());
        for (int i = 0; i < xValues.length; i++) {
          xValues[i] = frcCurve.get(i).getRadius() * conversion;
        }

        // The threshold curve is the same
        threshold = Frc.calculateThresholdCurve(frcCurve, thresholdMethod);
        addLine(colorThreshold, threshold);
      }

      for (int i = 0; i < xValues.length; i++) {
        yValues[i] = frcCurve.get(i).getCorrelation();
      }

      addLine(colorValues, yValues);
      addLine(colorNoSmooth, result.originalCorrelationCurve);
    }

    void addResolution(double resolution) {
      // Convert back to nm^-1
      final double x = 1 / resolution;

      // Find the intersection with the threshold line
      for (int i = 1; i < xValues.length; i++) {
        if (x < xValues[i]) {
          double correlation;
          // Interpolate
          final double upper = xValues[i];
          final double lower = xValues[i - 1];
          final double xx = (x - lower) / (upper - lower);
          correlation = threshold[i - 1] + xx * (threshold[i] - threshold[i - 1]);
          addResolution(resolution, correlation);
          return;
        }
      }
    }

    void addResolution(double resolution, double correlation) {
      addResolution(resolution, Double.NaN, correlation);
    }

    void addResolution(double resolution, double originalResolution, double correlation) {
      // Convert back to nm^-1
      final double x = 1 / resolution;
      plot.setColor(Color.MAGENTA);
      plot.drawLine(x, 0, x, correlation);
      plot.setColor(Color.BLACK);
      if (Double.isNaN(originalResolution)) {
        plot.addLabel(0, 0,
            String.format("Resolution = %s %s", MathUtils.rounded(resolution), units));
      } else {
        plot.addLabel(0, 0, String.format("Resolution = %s %s (Original = %s %s)",
            MathUtils.rounded(resolution), units, MathUtils.rounded(originalResolution), units));
      }
    }

    void addLine(Color color, double[] y) {
      if (color == null) {
        return;
      }
      plot.setColor(color);
      plot.addPoints(xValues, y, Plot.LINE);
    }

    Plot getPlot() {
      plot.setLimitsToFit(false);
      // Q. For some reason the limits calculated are ignored,
      // so set them as the defaults.
      // The FRC should not go above 1 so limit Y.
      // Auto-range:
      plot.setLimits(Double.NaN, Double.NaN, Double.NaN, 1.05);
      return plot;
    }
  }

  /**
   * Creates the frc curve.
   *
   * @param name the name
   * @param result the result
   * @param thresholdMethod the threshold method
   * @return the plot
   */
  Plot createFrcCurve(String name, FireResult result, ThresholdMethod thresholdMethod) {
    final FrcCurvePlot curve = new FrcCurvePlot();
    curve.add(name, result, thresholdMethod, Color.red, Color.blue, Color.black);
    curve.addResolution(result.fireNumber, result.correlation);
    return curve.getPlot();
  }

  private PlotWindow showFrcCurve(String name, FireResult result, ThresholdMethod thresholdMethod) {
    return showFrcCurve(name, result, thresholdMethod, 0);
  }

  private PlotWindow showFrcCurve(String name, FireResult result, ThresholdMethod thresholdMethod,
      int flags) {
    final Plot plot = createFrcCurve(name, result, thresholdMethod);
    return ImageJUtils.display(plot.getTitle(), plot, flags);
  }

  private void showFrcTimeEvolution(String name, double fireNumber, ThresholdMethod thresholdMethod,
      double fourierImageScale, int imageSize) {
    IJ.showStatus("Calculating FRC time evolution curve...");

    // Sort by time
    results.sort();

    final int nSteps = 10;
    int maxT = results.getLastFrame();
    if (maxT == 0) {
      maxT = results.size();
    }
    final int step = maxT / nSteps;

    final DoubleArrayList x = new DoubleArrayList();
    final DoubleArrayList y = new DoubleArrayList();

    double yMin = fireNumber;
    double yMax = fireNumber;

    final MemoryPeakResults newResults = new MemoryPeakResults();
    newResults.copySettings(results);
    int index = 0;

    final PeakResult[] list = results.toArray();
    for (int t = step; t <= maxT - step; t += step) {
      while (index < list.length) {
        final PeakResult r = list[index];
        if (r.getFrame() <= t) {
          newResults.add(r);
          index++;
        } else {
          break;
        }
      }

      x.add(t);

      final Fire f = this.copy();
      final FireResult result = f.calculateFireNumber(fourierMethod, samplingMethod,
          thresholdMethod, fourierImageScale, imageSize);
      final double fire = (result == null) ? 0 : result.fireNumber;
      y.add(fire);

      yMin = Math.min(yMin, fire);
      yMax = Math.max(yMax, fire);
    }

    // Add the final fire number
    x.add(maxT);
    y.add(fireNumber);

    final double[] xValues = x.toDoubleArray();
    final double[] yValues = y.toDoubleArray();

    String units = "px";
    if (results.getCalibration() != null) {
      nmPerUnit = results.getNmPerPixel();
      units = "nm";
    }

    final String title = name + " FRC Time Evolution";
    final Plot plot = new Plot(title, "Frames", "Resolution (" + units + ")");
    final double range = Math.max(1, yMax - yMin) * 0.05;
    plot.setLimits(xValues[0], xValues[xValues.length - 1], yMin - range, yMax + range);
    plot.setColor(Color.red);
    plot.addPoints(xValues, yValues, Plot.CONNECTED_CIRCLES);

    ImageJUtils.display(title, plot);
  }

  /**
   * Calculate the Fourier Image REsolution (FIRE) number using the chosen threshold method. Should
   * be called after {@link #initialise(MemoryPeakResults, MemoryPeakResults)}.
   *
   * @param fourierMethod the fourier method
   * @param samplingMethod the sampling method
   * @param thresholdMethod the threshold method
   * @param fourierImageScale The scale to use when reconstructing the super-resolution images (0
   *        for auto)
   * @param imageSize The width of the super resolution images when using auto scale (should be a
   *        power of two minus 1 for optimum memory usage)
   * @return The FIRE number
   */
  public FireResult calculateFireNumber(FourierMethod fourierMethod, SamplingMethod samplingMethod,
      ThresholdMethod thresholdMethod, double fourierImageScale, int imageSize) {
    final FireImages images = createImages(fourierImageScale, imageSize);
    return calculateFireNumber(fourierMethod, samplingMethod, thresholdMethod, images);
  }

  /**
   * Calculate the Fourier Image REsolution (FIRE) number using the chosen threshold method. Should
   * be called after {@link #initialise(MemoryPeakResults, MemoryPeakResults)}.
   *
   * @param fourierMethod the fourier method
   * @param samplingMethod the sampling method
   * @param thresholdMethod the threshold method
   * @param images the images
   * @return The FIRE number
   */
  public FireResult calculateFireNumber(FourierMethod fourierMethod, SamplingMethod samplingMethod,
      ThresholdMethod thresholdMethod, FireImages images) {
    if (images == null) {
      return null;
    }

    final Frc frc = new Frc();
    // Allow a progress tracker to be input.
    // This should be setup for the total number of repeats.
    // If parallelised then do not output the text status messages as they conflict.
    frc.setTrackProgress(progress);
    frc.setFourierMethod(fourierMethod);
    frc.setSamplingMethod(samplingMethod);
    frc.setPerimeterSamplingFactor(settings.perimeterSamplingFactor);
    final FrcCurve frcCurve = frc.calculateFrcCurve(images.ip1, images.ip2, images.nmPerPixel);
    if (frcCurve == null) {
      return null;
    }
    if (correctionQValue > 0) {
      Frc.applyQCorrection(frcCurve, correctionQValue, correctionMean, correctionSigma);
    }
    final double[] originalCorrelationCurve = frcCurve.getCorrelationValues();
    Frc.getSmoothedCurve(frcCurve, true);

    // Resolution in pixels
    final FrcFireResult result = Frc.calculateFire(frcCurve, thresholdMethod);
    if (result == null) {
      return new FireResult(Double.NaN, Double.NaN, frcCurve, originalCorrelationCurve);
    }
    final double fireNumber = result.fireNumber;

    // The FRC paper states that the super-resolution pixel size should be smaller
    // than 1/4 of R (the resolution).
    if (fireNumber > 0 && (images.nmPerPixel > fireNumber / 4)) {
      // Q. Should this be output somewhere else?
      ImageJUtils.log(
          "%s Warning: The super-resolution pixel size (%s) should be smaller than 1/4 of R"
              + " (the resolution %s)",
          pluginTitle, MathUtils.rounded(images.nmPerPixel), MathUtils.rounded(fireNumber));
    }

    return new FireResult(fireNumber, result.correlation, frcCurve, originalCorrelationCurve);
  }

  @SuppressWarnings("null")
  private void runQEstimation() {
    IJ.showStatus(pluginTitle + " ...");

    if (!showQEstimationInputDialog()) {
      return;
    }

    MemoryPeakResults inputResults =
        ResultsManager.loadInputResults(settings.inputOption, false, null, null);
    if (MemoryPeakResults.isEmpty(inputResults)) {
      IJ.error(pluginTitle, "No results could be loaded");
      return;
    }
    if (inputResults.getCalibration() == null) {
      IJ.error(pluginTitle, "The results are not calibrated");
      return;
    }

    inputResults = cropToRoi(inputResults);
    if (inputResults.size() <= 1) {
      IJ.error(pluginTitle, "No results within the crop region");
      return;
    }

    initialise(inputResults, null);

    // We need localisation precision.
    // Build a histogram of the localisation precision.
    // Get the initial mean and SD and plot as a Gaussian.
    final PrecisionHistogram histogram = calculatePrecisionHistogram();
    if (histogram == null) {
      IJ.error(pluginTitle, "No localisation precision available.\n \nPlease choose "
          + PrecisionMethod.FIXED + " and enter a precision mean and SD.");
      return;
    }

    final double fourierImageScale = Settings.SCALE_VALUES[settings.imageScaleIndex];
    final int imageSize = Settings.IMAGE_SIZE_VALUES[settings.imageSizeIndex];

    // Create the image and compute the numerator of FRC.
    // Do not use the signal so results.size() is the number of localisations.
    IJ.showStatus("Computing FRC curve ...");
    final FireImages images = createImages(fourierImageScale, imageSize, false);

    // DEBUGGING - Save the two images to disk. Load the images into the Matlab
    // code that calculates the Q-estimation and make this plugin match the functionality.
    // IJ.save(new ImagePlus("i1", images.ip1), "/scratch/i1.tif");
    // IJ.save(new ImagePlus("i2", images.ip2), "/scratch/i2.tif");

    final Frc frc = new Frc();
    frc.setTrackProgress(progress);
    frc.setFourierMethod(fourierMethod);
    frc.setSamplingMethod(samplingMethod);
    frc.setPerimeterSamplingFactor(settings.perimeterSamplingFactor);
    final FrcCurve frcCurve = frc.calculateFrcCurve(images.ip1, images.ip2, images.nmPerPixel);
    if (frcCurve == null) {
      IJ.error(pluginTitle, "Failed to compute FRC curve");
      return;
    }

    IJ.showStatus("Running Q-estimation ...");

    // Note:
    // The method implemented here is based on Matlab code provided by Bernd Rieger.
    // The idea is to compute the spurious correlation component of the FRC Numerator
    // using an initial estimate of distribution of the localisation precision (assumed
    // to be Gaussian). This component is the contribution of repeat localisations of
    // the same molecule to the numerator and is modelled as an exponential decay
    // (exp_decay). The component is scaled by the Q-value which
    // is the average number of times a molecule is seen in addition to the first time.
    // At large spatial frequencies the scaled component should match the numerator,
    // i.e. at high resolution (low FIRE number) the numerator is made up of repeat
    // localisations of the same molecule and not actual structure in the image.
    // The best fit is where the numerator equals the scaled component, i.e. num / (q*exp_decay) ==
    // 1.
    // The FRC Numerator is plotted and Q can be determined by
    // adjusting Q and the precision mean and SD to maximise the cost function.
    // This can be done interactively by the user with the effect on the FRC curve
    // dynamically updated and displayed.

    // Compute the scaled FRC numerator
    final double qNorm = (1 / frcCurve.mean1 + 1 / frcCurve.mean2);
    final double[] frcnum = new double[frcCurve.getSize()];
    for (int i = 0; i < frcnum.length; i++) {
      final FrcCurveResult r = frcCurve.get(i);
      frcnum[i] = qNorm * r.getNumerator() / r.getNumberOfSamples();
    }

    // Compute the spatial frequency and the region for curve fitting
    final double[] q = Frc.computeQ(frcCurve, false);
    int low = 0;
    int high = q.length;
    while (high > 0 && q[high - 1] > settings.maxQ) {
      high--;
    }
    while (low < q.length && q[low] < settings.minQ) {
      low++;
    }
    // Require we fit at least 10% of the curve
    if (high - low < q.length * 0.1) {
      IJ.error(pluginTitle, "Not enough points for Q estimation");
      return;
    }

    // Obtain initial estimate of Q plateau height and decay.
    // This can be done by fitting the precision histogram and then fixing the mean and sigma.
    // Or it can be done by allowing the precision to be sampled and the mean and sigma
    // become parameters for fitting.

    // Check if we can sample precision values
    final StoredDataStatistics precision = histogram.precision;
    final boolean sampleDecay = precision != null && settings.sampleDecay;

    double[] expDecay;
    if (sampleDecay) {
      // Random sample of precision values from the distribution is used to
      // construct the decay curve
      final int[] sample =
          RandomUtils.sample(10000, precision.getN(), UniformRandomProviders.create());

      final double fourPi2 = 4 * Math.PI * Math.PI;
      final double[] pre = new double[q.length];
      for (int i = 1; i < q.length; i++) {
        pre[i] = -fourPi2 * q[i] * q[i];
      }

      // Sample
      final int n = sample.length;
      final double[] hq = new double[n];
      for (int j = 0; j < n; j++) {
        // Scale to SR pixels
        double s2 = precision.getValue(sample[j]) / images.nmPerPixel;
        s2 *= s2;
        for (int i = 1; i < q.length; i++) {
          hq[i] += StdMath.exp(pre[i] * s2);
        }
      }
      for (int i = 1; i < q.length; i++) {
        hq[i] /= n;
      }

      expDecay = new double[q.length];
      expDecay[0] = 1;
      for (int i = 1; i < q.length; i++) {
        final double sincPiQ = sinc(Math.PI * q[i]);
        expDecay[i] = sincPiQ * sincPiQ * hq[i];
      }
    } else {
      // Note: The sigma mean and std should be in the units of super-resolution
      // pixels so scale to SR pixels
      expDecay = computeExpDecay(histogram.mean / images.nmPerPixel,
          histogram.sigma / images.nmPerPixel, q);
    }

    // Smoothing
    double[] smooth;
    if (settings.loessSmoothing) {
      // Note: This computes the log then smooths it
      final double bandwidth = 0.1;
      final int robustness = 0;
      final double[] l = new double[expDecay.length];
      for (int i = 0; i < l.length; i++) {
        // Original Matlab code computes the log for each array.
        // This is equivalent to a single log on the fraction of the two.
        // Perhaps the two log method is more numerically stable.
        // l[i] = Math.log(Math.abs(frcnum[i])) - Math.log(exp_decay[i]);
        l[i] = Math.log(Math.abs(frcnum[i] / expDecay[i]));
      }
      try {
        final LoessInterpolator loess = new LoessInterpolator(bandwidth, robustness);
        smooth = loess.smooth(q, l);
      } catch (final Exception ex) {
        IJ.error(pluginTitle, "LOESS smoothing failed");
        return;
      }
    } else {
      // Note: This smooths the curve before computing the log

      final double[] norm = new double[expDecay.length];
      for (int i = 0; i < norm.length; i++) {
        norm[i] = frcnum[i] / expDecay[i];
      }
      // Median window of 5 == radius of 2
      final DoubleMedianWindow mw = DoubleMedianWindow.wrap(norm, 2);
      smooth = new double[expDecay.length];
      for (int i = 0; i < norm.length; i++) {
        smooth[i] = Math.log(Math.abs(mw.getMedian()));
        mw.increment();
      }
    }

    // Fit with quadratic to find the initial guess.
    // Note: example Matlab code frc_Qcorrection7.m identifies regions of the
    // smoothed log curve with low derivative and only fits those. The fit is
    // used for the final estimate. Fitting a subset with low derivative is not
    // implemented here since the initial estimate is subsequently optimised
    // to maximise a cost function.
    final Quadratic curve = new Quadratic();
    final SimpleCurveFitter fit = SimpleCurveFitter.create(curve, new double[2]);
    final WeightedObservedPoints points = new WeightedObservedPoints();
    for (int i = low; i < high; i++) {
      points.add(q[i], smooth[i]);
    }
    final double[] estimate = fit.fit(points.toList());
    double qvalue = StdMath.exp(estimate[0]);

    // This could be made an option. Just use for debugging
    final boolean debug = false;
    if (debug) {
      // Plot the initial fit and the fit curve
      final double[] qScaled = Frc.computeQ(frcCurve, true);
      final double[] line = new double[q.length];
      for (int i = 0; i < q.length; i++) {
        line[i] = curve.value(q[i], estimate);
      }
      final String title = pluginTitle + " Initial fit";
      final Plot plot = new Plot(title, "Spatial Frequency (nm^-1)", "FRC Numerator");
      final String label = String.format("Q = %.3f", qvalue);
      plot.addPoints(qScaled, smooth, Plot.LINE);
      plot.setColor(Color.red);
      plot.addPoints(qScaled, line, Plot.LINE);
      plot.setColor(Color.black);
      plot.addLabel(0, 0, label);
      ImageJUtils.display(title, plot, ImageJUtils.NO_TO_FRONT);
    }

    if (settings.fitPrecision) {
      // Q - Should this be optional?
      if (sampleDecay) {
        // If a sample of the precision was used to construct the data for the initial fit
        // then update the estimate using the fit result since it will be a better start point.

        histogram.sigma = precision.getStandardDeviation();
        // Normalise sum-of-squares to the SR pixel size
        final double meanSumOfSquares =
            (precision.getSumOfSquares() / (images.nmPerPixel * images.nmPerPixel))
                / precision.getN();
        histogram.mean =
            images.nmPerPixel * Math.sqrt(meanSumOfSquares - estimate[1] / (4 * Math.PI * Math.PI));
      }

      // Do a multivariate fit ...
      final SimplexOptimizer opt = new SimplexOptimizer(1e-6, 1e-10);
      PointValuePair pair = null;
      final MultiPlateauness f = new MultiPlateauness(frcnum, q, low, high);
      final double[] initial =
          {histogram.mean / images.nmPerPixel, histogram.sigma / images.nmPerPixel, qvalue};
      pair = findMin(pair, opt, f, scale(initial, 0.1));
      pair = findMin(pair, opt, f, scale(initial, 0.5));
      pair = findMin(pair, opt, f, initial);
      pair = findMin(pair, opt, f, scale(initial, 2));
      pair = findMin(pair, opt, f, scale(initial, 10));

      if (pair != null) {
        final double[] point = pair.getPointRef();
        histogram.mean = point[0] * images.nmPerPixel;
        histogram.sigma = point[1] * images.nmPerPixel;
        qvalue = point[2];
      }
    } else {
      // Reset to theoretical curve. This is what will be used to compute the final correction.
      // TODO - check if the Matlab code uses a sampled curve to compute the correction.
      // If so then this should be optional.
      if (sampleDecay) {
        // If a sample of the precision was used to construct the data for the initial fit
        // then update the estimate using the fit result since it will be a better start point.

        if (precisionMethod != PrecisionMethod.FIXED) {
          histogram.sigma = precision.getStandardDeviation();
          // Normalise sum-of-squares to the SR pixel size
          final double meanSumOfSquares =
              (precision.getSumOfSquares() / (images.nmPerPixel * images.nmPerPixel))
                  / precision.getN();
          histogram.mean = images.nmPerPixel
              * Math.sqrt(meanSumOfSquares - estimate[1] / (4 * Math.PI * Math.PI));
        }

        expDecay = computeExpDecay(histogram.mean / images.nmPerPixel,
            histogram.sigma / images.nmPerPixel, q);
      }

      // Estimate spurious component by promoting plateauness.
      // The Matlab code used random initial points for a Simplex optimiser.

      // A Brent line search should be pretty deterministic so do simple repeats.
      // However it will proceed downhill so if the initial point is wrong then
      // it will find a sub-optimal result.
      final UnivariateOptimizer o = new BrentOptimizer(1e-3, 1e-6);
      final Plateauness f = new Plateauness(frcnum, expDecay, low, high);
      UnivariatePointValuePair result = null;
      result = findMin(result, o, f, qvalue, 0.1);
      result = findMin(result, o, f, qvalue, 0.2);
      result = findMin(result, o, f, qvalue, 0.333);
      result = findMin(result, o, f, qvalue, 0.5);

      // Do some Simplex repeats as well
      final SimplexOptimizer opt = new SimplexOptimizer(1e-6, 1e-10);
      result = findMin(result, opt, f, qvalue * 0.1);
      result = findMin(result, opt, f, qvalue * 0.5);
      result = findMin(result, opt, f, qvalue);
      result = findMin(result, opt, f, qvalue * 2);
      result = findMin(result, opt, f, qvalue * 10);

      if (result != null) {
        qvalue = result.getPoint();
      }
    }

    final QPlot qplot = new QPlot(frcCurve, qvalue, low, high);

    // Interactive dialog to estimate Q (blinking events per flourophore) using
    // sliders for the mean and standard deviation of the localisation precision.
    showQEstimationDialog(histogram, qplot, images.nmPerPixel);

    IJ.showStatus(pluginTitle + " complete");
  }

  private static double[] scale(double[] a, double factor) {
    a = a.clone();
    for (int i = 0; i < a.length; i++) {
      a[i] *= factor;
    }
    return a;
  }

  /**
   * Compute the exponential decay.
   *
   * @param mean the mean
   * @param sigma the sigma
   * @param qvalues the qvalues
   * @return the decay
   */
  static double[] computeExpDecay(double mean, double sigma, double[] qvalues) {
    final double[] hq = Frc.computeHq(qvalues, mean, sigma);
    final double[] expDecay = new double[qvalues.length];
    expDecay[0] = 1;
    for (int i = 1; i < qvalues.length; i++) {
      final double sincPiQ = sinc(Math.PI * qvalues[i]);
      expDecay[i] = sincPiQ * sincPiQ * hq[i];
    }
    return expDecay;
  }

  /**
   * Compute the Sinc function.
   *
   * @param value the value
   * @return the sinc value
   */
  static double sinc(double value) {
    return Math.sin(value) / value;
  }

  /**
   * Quadratic function.
   */
  private static class Quadratic implements ParametricUnivariateFunction {
    @Override
    public double value(double x, double... parameters) {
      return parameters[0] + parameters[1] * x * x;
    }

    @Override
    public double[] gradient(double x, double... parameters) {
      return new double[] {1, x * x};
    }
  }

  private static UnivariatePointValuePair findMin(UnivariatePointValuePair current,
      UnivariateOptimizer optimiser, UnivariateFunction func, double qvalue, double factor) {
    try {
      final BracketFinder bracket = new BracketFinder();
      bracket.search(func, GoalType.MINIMIZE, qvalue * factor, qvalue / factor);
      final UnivariatePointValuePair next = optimiser.optimize(GoalType.MINIMIZE, new MaxEval(3000),
          new SearchInterval(bracket.getLo(), bracket.getHi(), bracket.getMid()),
          new UnivariateObjectiveFunction(func));
      if (next == null) {
        return current;
      }
      if (current != null) {
        return (next.getValue() < current.getValue()) ? next : current;
      }
      return next;
    } catch (final Exception ex) {
      return current;
    }
  }

  private static UnivariatePointValuePair findMin(UnivariatePointValuePair current,
      SimplexOptimizer optimiser, MultivariateFunction func, double qvalue) {
    try {
      final NelderMeadSimplex simplex = new NelderMeadSimplex(1);
      final double[] initialSolution = {qvalue};
      final PointValuePair solution =
          optimiser.optimize(new MaxEval(1000), new InitialGuess(initialSolution), simplex,
              new ObjectiveFunction(func), GoalType.MINIMIZE);
      final UnivariatePointValuePair next = (solution == null) ? null
          : new UnivariatePointValuePair(solution.getPointRef()[0], solution.getValue());
      if (next == null) {
        return current;
      }
      if (current != null) {
        return (next.getValue() < current.getValue()) ? next : current;
      }
      return next;
    } catch (final Exception ex) {
      return current;
    }
  }

  private static PointValuePair findMin(PointValuePair current, SimplexOptimizer optimiser,
      MultivariateFunction func, double[] initialSolution) {
    try {
      final NelderMeadSimplex simplex = new NelderMeadSimplex(initialSolution.length);
      final PointValuePair next =
          optimiser.optimize(new MaxEval(1000), new InitialGuess(initialSolution), simplex,
              new ObjectiveFunction(func), GoalType.MINIMIZE);
      if (next == null) {
        return current;
      }
      if (current != null) {
        return (next.getValue() < current.getValue()) ? next : current;
      }
      return next;
    } catch (final Exception ex) {
      return current;
    }
  }

  /**
   * Univariate plateauness function.
   */
  private static class Plateauness implements UnivariateFunction, MultivariateFunction {
    private static final double EPSILON = Math.ulp(1.0);
    static final double FRCNUM_NOISE_VAR = 0.1;
    final double[] pre;
    final double n2;

    /**
     * Create an instance.
     *
     * @param frcnum the scaled FRC numerator
     * @param expDecay the precomputed exponential decay (hq)
     * @param low the lower bound of the array for optimisation
     * @param high the higher bound of the array for optimisation
     */
    Plateauness(double[] frcnum, double[] expDecay, int low, int high) {
      // Precompute
      pre = new double[high - low];
      for (int i = 0; i < pre.length; i++) {
        final int index = i + low;
        pre[i] = frcnum[index] / expDecay[index];
      }
      n2 = FRCNUM_NOISE_VAR * FRCNUM_NOISE_VAR;
    }

    @Override
    public double value(double qvalue) {
      if (qvalue < EPSILON) {
        qvalue = EPSILON;
      }
      double value = 0;
      for (final double p : pre) {
        // Original cost function. Note that each observation has a
        // contribution of 0 to 1.
        final double diff = (p / qvalue) - 1;
        value += 1 - StdMath.exp(-diff * diff / n2);

        // Modified cost function so that the magnitude of difference over or
        // under 1 is penalised the same. This has a problem if FRC numerator
        // is negative. Also the range is unchecked so observation can have
        // unequal contributions.
        // double diff = Math.abs(pre[i]) / qvalue;
        // v += Math.abs(Math.log(diff));
      }
      return value;
    }

    @Override
    public double value(double[] point) {
      return value(point[0]);
    }
  }

  /**
   * Multivariate plateauness function.
   */
  private static class MultiPlateauness implements MultivariateFunction {
    private static final double EPSILON = Math.ulp(1.0);
    static final double FRCNUM_NOISE_VAR = 0.1;
    static final double FOUR_PI2 = 4 * Math.PI * Math.PI;
    final double[] pre;
    final double[] q2;
    final double n2;

    /**
     * Create an instance.
     *
     * @param frcnum the scaled FRC numerator
     * @param qvalues the precomputed exponential decay (hq)
     * @param low the lower bound of the array for optimisation
     * @param high the higher bound of the array for optimisation
     */
    MultiPlateauness(double[] frcnum, double[] qvalues, int low, int high) {
      q2 = new double[qvalues.length];

      // Precompute
      pre = new double[high - low];
      for (int i = 0; i < pre.length; i++) {
        final int index = i + low;
        final double sincPiQ = (index == 0) ? 1 : sinc(Math.PI * qvalues[index]);
        pre[i] = frcnum[index] / (sincPiQ * sincPiQ);
        q2[i] = qvalues[index] * qvalues[index];
      }
      n2 = FRCNUM_NOISE_VAR * FRCNUM_NOISE_VAR;
    }

    @Override
    public double value(double[] point) {
      final double mean = point[0];
      final double sigma = point[1];
      double qvalue = point[2];

      if (qvalue < EPSILON) {
        qvalue = EPSILON;
      }

      // Fast computation of a subset of hq
      final double eightPi2S2 = 2 * FOUR_PI2 * sigma * sigma;
      final double factor = -FOUR_PI2 * mean * mean;

      double value = 0;
      for (int i = 0; i < pre.length; i++) {
        final double d = 1 + eightPi2S2 * q2[i];
        final double hq = StdMath.exp((factor * q2[i]) / d) / Math.sqrt(d);

        // Original cost function. Note that each observation has a
        // contribution of 0 to 1.
        final double diff = (pre[i] / (qvalue * hq)) - 1;
        value += 1 - StdMath.exp(-diff * diff / n2);
      }
      return value;
    }
  }

  private boolean showQEstimationInputDialog() {
    ExtendedGenericDialog gd = new ExtendedGenericDialog(pluginTitle);
    gd.addHelp(HelpUrls.getUrl("fourier-image-resolution"));

    // Build a list of all images with a region ROI
    final List<String> titles = new LinkedList<>();
    for (final int imageId : ImageJUtils.getIdList()) {
      final ImagePlus imp = WindowManager.getImage(imageId);
      if (imp != null && imp.getRoi() != null && imp.getRoi().isArea()) {
        titles.add(imp.getTitle());
      }
    }

    gd.addMessage("Estimate the blinking correction parameter Q for Fourier Ring Correlation");

    ResultsManager.addInput(gd, settings.inputOption, InputSource.MEMORY);
    if (!titles.isEmpty()) {
      gd.addCheckbox((titles.size() == 1) ? "Use_ROI" : "Choose_ROI", settings.chooseRoi);
    }

    gd.addMessage("Image construction options:");
    gd.addChoice("Image_scale", Settings.SCALE_ITEMS, settings.imageScaleIndex);
    gd.addChoice("Auto_image_size", Settings.IMAGE_SIZE_ITEMS, settings.imageSizeIndex);
    gd.addNumericField("Block_size", settings.blockSize, 0);
    gd.addCheckbox("Random_split", settings.randomSplit);
    gd.addNumericField("Max_per_bin", settings.maxPerBin, 0);

    gd.addMessage("Fourier options:");
    final String[] fourierMethodNames =
        SettingsManager.getNames((Object[]) Frc.FourierMethod.values());
    gd.addChoice("Fourier_method", fourierMethodNames,
        fourierMethodNames[settings.fourierMethodIndex]);
    final String[] samplingMethodNames =
        SettingsManager.getNames((Object[]) Frc.SamplingMethod.values());
    gd.addChoice("Sampling_method", samplingMethodNames,
        samplingMethodNames[settings.samplingMethodIndex]);
    gd.addSlider("Sampling_factor", 0.2, 4, settings.perimeterSamplingFactor);

    gd.addMessage("Estimation options:");
    final String[] thresholdMethodNames =
        SettingsManager.getNames((Object[]) Frc.ThresholdMethod.values());
    gd.addChoice("Threshold_method", thresholdMethodNames,
        thresholdMethodNames[settings.thresholdMethodIndex]);
    final String[] precisionMethodNames =
        SettingsManager.getNames((Object[]) PrecisionMethod.values());
    gd.addChoice("Precision_method", precisionMethodNames,
        precisionMethodNames[settings.precisionMethodIndex]);
    gd.addNumericField("Precision_Mean", settings.mean, 2, 6, "nm");
    gd.addNumericField("Precision_Sigma", settings.sigma, 2, 6, "nm");
    gd.addCheckbox("Sample_decay", settings.sampleDecay);
    gd.addCheckbox("LOESS_smoothing", settings.loessSmoothing);
    gd.addCheckbox("Fit_precision", settings.fitPrecision);
    gd.addSlider("MinQ", 0, 0.4, settings.minQ);
    gd.addSlider("MaxQ", 0.1, 0.5, settings.maxQ);

    gd.addHelp(HelpUrls.getUrl("fire-q-estimation"));
    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.inputOption = ResultsManager.getInputSource(gd);
    if (!titles.isEmpty()) {
      settings.chooseRoi = gd.getNextBoolean();
    }

    settings.imageScaleIndex = gd.getNextChoiceIndex();
    settings.imageSizeIndex = gd.getNextChoiceIndex();
    settings.blockSize = Math.max(1, (int) gd.getNextNumber());
    settings.randomSplit = gd.getNextBoolean();
    settings.maxPerBin = Math.abs((int) gd.getNextNumber());

    settings.fourierMethodIndex = gd.getNextChoiceIndex();
    fourierMethod = FourierMethod.values()[settings.fourierMethodIndex];
    settings.samplingMethodIndex = gd.getNextChoiceIndex();
    samplingMethod = SamplingMethod.values()[settings.samplingMethodIndex];
    settings.perimeterSamplingFactor = gd.getNextNumber();

    settings.thresholdMethodIndex = gd.getNextChoiceIndex();
    thresholdMethod = Frc.ThresholdMethod.values()[settings.thresholdMethodIndex];
    settings.precisionMethodIndex = gd.getNextChoiceIndex();
    precisionMethod = PrecisionMethod.values()[settings.precisionMethodIndex];
    settings.mean = Math.abs(gd.getNextNumber());
    settings.sigma = Math.abs(gd.getNextNumber());
    settings.sampleDecay = gd.getNextBoolean();
    settings.loessSmoothing = gd.getNextBoolean();
    settings.fitPrecision = gd.getNextBoolean();
    settings.minQ = MathUtils.clip(0, 0.5, gd.getNextNumber());
    settings.maxQ = MathUtils.clip(0, 0.5, gd.getNextNumber());

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Perimeter sampling factor", settings.perimeterSamplingFactor);
      if (precisionMethod == PrecisionMethod.FIXED) {
        ParameterUtils.isAboveZero("Precision Mean", settings.mean);
        ParameterUtils.isAboveZero("Precision Sigma", settings.sigma);
      }
      ParameterUtils.isAbove("MaxQ", settings.maxQ, settings.minQ);
    } catch (final IllegalArgumentException ex) {
      IJ.error(pluginTitle, ex.getMessage());
      return false;
    }

    if (!titles.isEmpty() && settings.chooseRoi) {
      if (titles.size() == 1) {
        settings.roiImage = titles.get(0);
        Recorder.recordOption("Image", settings.roiImage);
      } else {
        final String[] items = titles.toArray(new String[0]);
        gd = new ExtendedGenericDialog(pluginTitle);
        gd.addMessage("Select the source image for the ROI");
        gd.addChoice("Image", items, settings.roiImage);
        gd.showDialog();
        if (gd.wasCanceled()) {
          return false;
        }
        settings.roiImage = gd.getNextChoice();
      }
      final ImagePlus imp = WindowManager.getImage(settings.roiImage);

      roiBounds = imp.getRoi().getBounds();
      roiImageWidth = imp.getWidth();
      roiImageHeight = imp.getHeight();
    } else {
      roiBounds = null;
    }

    return true;
  }

  /**
   * Represent the Q-plot data.
   */
  private class QPlot {
    final FrcCurve frcCurve;
    final double nmPerPixel;
    final double qnorm;
    final double[] vq;
    final double[] sinc2;
    final double[] qvalues;
    final double[] qscaled;
    final int low;
    final int high;
    String title;
    String title2;
    FireResult originalFireResult;
    double originalFireNumber = Double.NaN;

    // Store the last plotted value
    double mean;
    double sigma;
    double qvalue;

    QPlot(FrcCurve frcCurve, double qvalue, int low, int high) {
      this.nmPerPixel = frcCurve.nmPerPixel;
      this.frcCurve = frcCurve;
      this.qvalue = qvalue;
      this.low = low;
      this.high = high;

      // For normalisation
      qnorm = (1 / frcCurve.mean1 + 1 / frcCurve.mean2);

      // Compute v(q) - The numerator of the FRC divided by the number of pixels
      // in the Fourier circle (2*pi*q*L)
      vq = new double[frcCurve.getSize()];
      for (int i = 0; i < vq.length; i++) {
        vq[i] = qnorm * frcCurve.get(i).getNumerator() / frcCurve.get(i).getNumberOfSamples();
      }

      qvalues = Frc.computeQ(frcCurve, false);
      // For the plot
      qscaled = Frc.computeQ(frcCurve, true);

      // Compute sinc factor
      sinc2 = new double[frcCurve.getSize()];
      sinc2[0] = 1; // By definition
      for (int i = 1; i < sinc2.length; i++) {
        // Note the original equation given in the paper: sinc(pi*q*L)^2 is a typo.
        // Matlab code provided by Bernd Rieger removes L to compute: sinc(pi*q)^2
        // with q == 1/L, 2/L, ... (i.e. no unit conversion to nm). This means that
        // the function will start at 1 and drop off to zero at q=L.

        // sinc(pi*q)^2
        sinc2[i] = sinc(Math.PI * qvalues[i]);
        sinc2[i] *= sinc2[i];
      }

      // For the plot
      title = results.getName() + " FRC Numerator Curve";
      title2 = results.getName() + " FRC Numerator/Correction Ratio";

      // Reset
      Frc.applyQCorrection(frcCurve, 0, 0, 0);
      final FrcCurve smoothedFrcCurve = Frc.getSmoothedCurve(frcCurve, false);
      originalFireResult = new FireResult(0, 0, smoothedFrcCurve, frcCurve.getCorrelationValues());
      final FrcFireResult result = Frc.calculateFire(smoothedFrcCurve, thresholdMethod);
      if (result != null) {
        originalFireNumber = result.fireNumber;
      }
    }

    private double sinc(double x) {
      return Math.sin(x) / x;
    }

    void plot(double mean, double sigma, double qvalue, MyWindowOrganiser windowOrganiser) {
      this.mean = mean;
      this.sigma = sigma;
      this.qvalue = qvalue;

      final double mu = mean / nmPerPixel;
      final double sd = sigma / nmPerPixel;

      final double[] hq = Frc.computeHq(qvalues, mu, sd);
      final double[] correction = new double[hq.length];
      final double[] vqCorr = new double[vq.length];
      final double[] ratio = new double[vq.length];

      for (int i = 0; i < hq.length; i++) {
        // Note: vq already has the qNorm factor applied so we do not
        // divide qvalue by qNorm.
        correction[i] = qvalue * sinc2[i] * hq[i];
        // This is not actually the corrected numerator since it is made absolute
        // vq_corr[i] = Math.abs(vq[i] - correction[i]);
        vqCorr[i] = vq[i] - correction[i];
        ratio[i] = vq[i] / correction[i];
      }

      // Add this to aid is manual adjustment
      final double plateauness = computePlateauness(qvalue, mu, sd);

      Plot plot = new Plot(title, "Spatial Frequency (nm^-1)", "FRC Numerator");

      String label = String.format("Q = %.3f (Precision = %.3f +/- %.3f)", qvalue, mean, sigma);
      plot.setColor(Color.red);
      final double[] plotVq = makeStrictlyPositive(this.vq, Double.POSITIVE_INFINITY);
      plot.addPoints(qscaled, plotVq, Plot.LINE);
      final double min = MathUtils.min(plotVq);
      if (qvalue > 0) {
        label += String.format(". Cost = %.3f", plateauness);
        plot.setColor(Color.darkGray);
        plot.addPoints(qscaled, correction, Plot.DOT);
        plot.setColor(Color.blue);
        plot.addPoints(qscaled, makeStrictlyPositive(vqCorr, min), Plot.LINE);
        plot.setColor(Color.black);
        plot.addLegend("Numerator\nCorrection\nCorrected Numerator", "top-right");
      }
      plot.setColor(Color.magenta);
      plot.drawLine(qscaled[low], min, qscaled[low], plotVq[0]);
      plot.drawLine(qscaled[high], min, qscaled[high], plotVq[0]);
      plot.setColor(Color.black);
      plot.addLabel(0, 0, label);

      plot.setAxisYLog(true);
      windowOrganiser.display(title, plot, ImageJUtils.NO_TO_FRONT);
      plot.setLimitsToFit(true); // For the log scale this seems to only work after drawing

      // Show how the resolution changes

      Frc.applyQCorrection(frcCurve, qvalue, mean, sigma);
      final FrcCurve smoothedFrcCurve = Frc.getSmoothedCurve(frcCurve, false);

      // Resolution in pixels
      final FrcFireResult result = Frc.calculateFire(smoothedFrcCurve, thresholdMethod);
      if (result != null) {
        final double fireNumber = result.fireNumber;

        final FrcCurvePlot curve = new FrcCurvePlot();
        final FireResult fireResult = new FireResult(fireNumber, result.correlation,
            smoothedFrcCurve, frcCurve.getCorrelationValues());
        double orig = Double.NaN;
        if (qvalue > 0) {
          curve.add(results.getName(), originalFireResult, thresholdMethod, Color.orange,
              Color.blue, Color.lightGray);
          orig = originalFireNumber;
        }
        curve.add(results.getName(), fireResult, thresholdMethod, Color.red, Color.blue,
            Color.black);
        curve.addResolution(fireNumber, orig, result.correlation);
        plot = curve.getPlot();

        windowOrganiser.display(plot.getTitle(), plot, ImageJUtils.NO_TO_FRONT);
      }

      // Produce a ratio plot. Plateauness is designed to achieve a value of 1 for this ratio.
      plot = new Plot(title2, "Spatial Frequency (nm^-1)", "FRC Numerator / Spurious component");
      final double xMax = qscaled[qscaled.length - 1];
      if (qvalue > 0) {
        plot.addLabel(0, 0, String.format("Cost = %.3f", plateauness));
        plot.setColor(Color.blue);
        plot.addPoints(qscaled, ratio, Plot.LINE);
      }
      plot.setColor(Color.black);
      plot.drawLine(0, 1, xMax, 1);
      plot.setColor(Color.magenta);
      plot.drawLine(qscaled[low], 0, qscaled[low], 2);
      plot.drawLine(qscaled[high], 0, qscaled[high], 2);
      plot.setLimits(0, xMax, 0, 2);
      windowOrganiser.display(title2, plot, ImageJUtils.NO_TO_FRONT);

      windowOrganiser.tile();
    }

    double computePlateauness(double qvalue, double mu, double sd) {
      final double[] expDecay = computeExpDecay(mu, sd, qvalues);
      final Plateauness p = new Plateauness(vq, expDecay, low, high);
      return p.value(qvalue);
    }

    private double[] makeStrictlyPositive(double[] data, double min) {
      data = data.clone();
      if (min == Double.POSITIVE_INFINITY) {
        // Find min positive value
        for (final double d : data) {
          if (d > 0 && d < min) {
            min = d;
          }
        }
      }
      for (int i = 0; i < data.length; i++) {
        if (data[i] < min) {
          data[i] = min;
        }
      }
      return data;
    }
  }

  /**
   * Represent the precision histogram.
   */
  private static class PrecisionHistogram {
    final float[] x;
    final float[] y;
    final String title;
    final double standardAmplitude;
    final float[] x2;
    final StoredDataStatistics precision;

    /**
     * The mean of the localisation precision distribution (in nm). This value can be updated by the
     * {@link #plot(double, double, MyWindowOrganiser)} method.
     */
    double mean;

    /**
     * The standard deviation of the localisation precision distribution (in nm). This value can be
     * updated by the {@link #plot(double, double, MyWindowOrganiser)} method.
     */
    double sigma;

    PrecisionHistogram(float[][] hist, StoredDataStatistics precision, String title) {
      this.title = title;
      x = hist[0];
      y = hist[1];
      this.precision = precision;

      // Sum the area under the histogram to use for normalisation.
      // Amplitude = volume / (sigma * sqrt(2*pi))
      // Precompute the correct amplitude for a standard width Gaussian
      double dx = (x[1] - x[0]);
      standardAmplitude = precision.getN() * dx / Math.sqrt(2 * Math.PI);

      // Set up for drawing the Gaussian curve
      final double min = x[0];
      final double max = x[x.length - 1];
      final int n = 100;
      dx = (max - min) / n;
      x2 = new float[n + 1];
      for (int i = 0; i <= n; i++) {
        x2[i] = (float) (min + i * dx);
      }
    }

    /**
     * Instantiates a new precision histogram.
     *
     * @param title the title
     */
    public PrecisionHistogram(String title) {
      this.title = title;
      // Set some defaults
      this.mean = 20;
      this.sigma = 2;
      x = y = x2 = null;
      precision = null;
      standardAmplitude = 0;
    }

    void plot(double mean, double sigma, MyWindowOrganiser windowOrganiser) {
      this.mean = mean;
      this.sigma = sigma;
      plot(windowOrganiser);
    }

    void plot(MyWindowOrganiser windowOrganiser) {
      final Plot plot = new Plot(title, "Precision (nm)", "Frequency");
      if (x == null) {
        // There is no base histogram.
        // Just plot a Gaussian +/- 4 SD.
        plot.addLabel(0, 0, String.format("Precision = %.3f +/- %.3f", mean, sigma));
        final double min = Math.max(0, mean - 4 * sigma);
        final double max = mean + 4 * sigma;
        final int n = 100;
        final double dx = (max - min) / n;
        final float[] xdata = new float[n + 1];
        final Gaussian g = new Gaussian(1, mean, sigma);
        final float[] ydata = new float[xdata.length];
        for (int i = 0; i <= n; i++) {
          xdata[i] = (float) (min + i * dx);
          ydata[i] = (float) g.value(xdata[i]);
        }
        plot.setColor(Color.red);
        plot.addPoints(xdata, ydata, Plot.LINE);

        // Always put min = 0 otherwise the plot does not change.
        plot.setLimits(0, max, 0, 1.05);
      } else {
        plot.setColor(Color.black);
        plot.addPoints(x, y, Plot.BAR);
        plot.addLabel(0, 0, String.format("Precision = %.3f +/- %.3f", mean, sigma));
        // Add the Gaussian line
        // Compute the integral of the standard gaussian between the min and max
        final double denom0 = 1.0 / (Math.sqrt(2.0) * sigma);
        final double integral =
            0.5 * Erf.erf((x2[0] - mean) * denom0, (x2[x2.length - 1] - mean) * denom0);
        // Normalise so the integral has the same volume as the histogram
        final Gaussian g = new Gaussian(this.standardAmplitude / (sigma * integral), mean, sigma);
        final float[] ydata = new float[x2.length];
        for (int i = 0; i < ydata.length; i++) {
          ydata[i] = (float) g.value(x2[i]);
        }
        // Normalise
        plot.setColor(Color.red);
        plot.addPoints(x2, ydata, Plot.LINE);
        float max = MathUtils.max(ydata);
        max = MathUtils.maxDefault(max, y);
        final double rangex = 0;
        plot.setLimits(x2[0] - rangex, x2[x2.length - 1] + rangex, 0, max * 1.05);
      }
      windowOrganiser.display(title, plot, ImageJUtils.NO_TO_FRONT);
      windowOrganiser.tile();
    }
  }

  /**
   * Calculate a histogram of the precision. The precision can be either stored in the results or
   * calculated using the Mortensen formula. If the precision method for Q estimation is not fixed
   * then the histogram is fitted with a Gaussian to create an initial estimate.
   *
   * @return The precision histogram
   */
  private PrecisionHistogram calculatePrecisionHistogram() {
    final String title = results.getName() + " Precision Histogram";

    // Check if the results has the precision already or if it can be computed.
    final boolean canUseStored = canUseStoredPrecision(results);
    final boolean canCalculatePrecision = canCalculatePrecision(results);

    // Set the method to compute a histogram. Default to the user selected option.
    PrecisionMethod method = null;
    if ((canUseStored && precisionMethod == PrecisionMethod.STORED)
        || (canCalculatePrecision && precisionMethod == PrecisionMethod.CALCULATE)) {
      method = precisionMethod;
    }

    if (method == null) {
      // We get here if the choice of the user is not available.
      // We only have two choices so if one is available then select it.
      if (canUseStored) {
        method = PrecisionMethod.STORED;
      } else if (canCalculatePrecision) {
        method = PrecisionMethod.CALCULATE;
      }
      // If the user selected a method not available then log a warning
      if (method != null && precisionMethod != PrecisionMethod.FIXED) {
        IJ.log(String.format("%s : Selected precision method '%s' not available, switching to '%s'",
            pluginTitle, precisionMethod, method.getName()));
      }

      if (method == null) {
        // We cannot compute a precision histogram.
        // This does not matter if the user has provide a fixed input.
        if (precisionMethod == PrecisionMethod.FIXED) {
          final PrecisionHistogram histogram = new PrecisionHistogram(title);
          histogram.mean = settings.mean;
          histogram.sigma = settings.sigma;
          return histogram;
        }
        // No precision
        return null;
      }
    }

    // We get here if we can compute precision.
    // Build the histogram
    StoredDataStatistics precision = new StoredDataStatistics(results.size());
    if (method == PrecisionMethod.STORED) {
      final StoredDataStatistics p = precision;
      results.forEach((PeakResultProcedure) result -> p.add(result.getPrecision()));
    } else {
      precision.add(pp.precisions);
    }

    double yMin = Double.NEGATIVE_INFINITY;
    double yMax = 0;

    // Set the min and max y-values using 1.5 x IQR
    final DescriptiveStatistics stats = precision.getStatistics();
    final double lower = stats.getPercentile(25);
    final double upper = stats.getPercentile(75);
    final boolean logFitParameters = IJ.debugMode;
    if (Double.isNaN(lower) || Double.isNaN(upper)) {
      if (logFitParameters) {
        ImageJUtils.log("Error computing IQR: %f - %f", lower, upper);
      }
    } else {
      final double iqr = upper - lower;

      yMin = Math.max(lower - iqr, stats.getMin());
      yMax = Math.min(upper + iqr, stats.getMax());

      if (logFitParameters) {
        ImageJUtils.log("  Data range: %f - %f. Plotting 1.5x IQR: %f - %f", stats.getMin(),
            stats.getMax(), yMin, yMax);
      }
    }

    if (yMin == Double.NEGATIVE_INFINITY) {
      final int n = 5;
      yMin = Math.max(stats.getMin(), stats.getMean() - n * stats.getStandardDeviation());
      yMax = Math.min(stats.getMax(), stats.getMean() + n * stats.getStandardDeviation());

      if (logFitParameters) {
        ImageJUtils.log("  Data range: %f - %f. Plotting mean +/- %dxSD: %f - %f", stats.getMin(),
            stats.getMax(), n, yMin, yMax);
      }
    }

    // Get the data within the range
    final double[] data = precision.getValues();
    precision = new StoredDataStatistics(data.length);
    for (final double d : data) {
      if (d < yMin || d > yMax) {
        continue;
      }
      precision.add(d);
    }

    final int histogramBins = HistogramPlot.getBins(precision, HistogramPlot.BinMethod.SCOTT);
    final float[][] hist =
        HistogramPlot.calcHistogram(precision.getFloatValues(), yMin, yMax, histogramBins);
    final PrecisionHistogram histogram = new PrecisionHistogram(hist, precision, title);

    if (precisionMethod == PrecisionMethod.FIXED) {
      histogram.mean = settings.mean;
      histogram.sigma = settings.sigma;
      return histogram;
    }

    // Fitting of the histogram to produce the initial estimate

    // Extract non-zero data
    float[] x = Arrays.copyOf(hist[0], hist[0].length);
    float[] y = Arrays.copyOf(hist[1], hist[1].length);
    int count = 0;
    for (int i = 0; i < y.length; i++) {
      if (y[i] > 0) {
        x[count] = x[i];
        y[count] = y[i];
        count++;
      }
    }
    x = Arrays.copyOf(x, count);
    y = Arrays.copyOf(y, count);

    // Sense check to fitted data. Get mean and SD of histogram
    final double[] stats2 = HistogramPlot.getHistogramStatistics(x, y);
    if (logFitParameters) {
      ImageJUtils.log("  Initial Statistics: %f +/- %f", stats2[0], stats2[1]);
    }
    histogram.mean = stats2[0];
    histogram.sigma = stats2[1];

    // Standard Gaussian fit
    final double[] parameters = fitGaussian(x, y);
    if (parameters == null) {
      ImageJUtils.log("  Failed to fit initial Gaussian");
      return histogram;
    }
    final double newMean = parameters[1];
    final double error = Math.abs(stats2[0] - newMean) / stats2[1];
    if (error > 3) {
      ImageJUtils.log("  Failed to fit Gaussian: %f standard deviations from histogram mean",
          error);
      return histogram;
    }
    if (newMean < yMin || newMean > yMax) {
      ImageJUtils.log("  Failed to fit Gaussian: %f outside data range %f - %f", newMean, yMin,
          yMax);
      return histogram;
    }

    if (logFitParameters) {
      ImageJUtils.log("  Initial Gaussian: %f @ %f +/- %f", parameters[0], parameters[1],
          parameters[2]);
    }

    histogram.mean = parameters[1];
    histogram.sigma = parameters[2];

    return histogram;
  }

  private boolean canCalculatePrecision(MemoryPeakResults results) {
    try {
      pp = new PrecisionResultProcedure(results);
      pp.getLsePrecision();
    } catch (final DataException ex) {
      return false;
    }

    // Check they are different
    for (int i = 0; i < pp.size(); i++) {
      // Check this is valid
      if (Double.isFinite(pp.precisions[i])) {
        final double p1 = pp.precisions[i];
        for (int j = i + 1; j < pp.size(); j++) {
          if (Double.isFinite(pp.precisions[j]) && pp.precisions[j] != p1) {
            return true;
          }
        }
        // All the results are the same, this is not valid
        break;
      }
    }

    return false;
  }

  private static boolean canUseStoredPrecision(MemoryPeakResults results) {
    return results.hasPrecision();
  }

  /**
   * Fit gaussian.
   *
   * @param x the x
   * @param y the y
   * @return new double[] { norm, mean, sigma }
   */
  private static double[] fitGaussian(float[] x, float[] y) {
    final WeightedObservedPoints obs = new WeightedObservedPoints();
    for (int i = 0; i < x.length; i++) {
      obs.add(x[i], y[i]);
    }

    final Collection<WeightedObservedPoint> observations = obs.toList();
    final GaussianCurveFitter fitter = GaussianCurveFitter.create().withMaxIterations(2000);
    final GaussianCurveFitter.ParameterGuesser guess =
        new GaussianCurveFitter.ParameterGuesser(observations);
    double[] initialGuess = null;
    try {
      initialGuess = guess.guess();
      return fitter.withStartPoint(initialGuess).fit(observations);
    } catch (final Exception ignored) {
      // We are expecting TooManyEvaluationsException.
      // Just in case there is another exception type, or the initial estimate failed
      // we catch all exceptions.
    }
    return initialGuess;
  }

  /**
   * Used to tile the windows from the worker threads on the first plot.
   */
  private static class MyWindowOrganiser {
    final WindowOrganiser windowOrganiser = new WindowOrganiser();
    AtomicInteger size;

    MyWindowOrganiser(int size) {
      this.size = new AtomicInteger(size);
    }

    void display(String title, Plot plot, int flags) {
      if (isTiled()) {
        ImageJUtils.display(title, plot, flags);
      } else {
        // The windows have not yet been tiled so track all new windows
        final WindowOrganiser localWindowOrganiser = new WindowOrganiser();
        ImageJUtils.display(title, plot, flags, localWindowOrganiser);
        if (localWindowOrganiser.isNotEmpty()) {
          synchronized (windowOrganiser) {
            windowOrganiser.add(localWindowOrganiser);
          }
        }
      }
    }

    boolean isTiled() {
      return size.get() == 0;
    }

    void tile() {
      // Get the value and count down to zero
      if (size.getAndUpdate(MyWindowOrganiser::decrementToZero) == 1) {
        // When the value is 1 the next value is 0 so tile
        synchronized (windowOrganiser) {
          windowOrganiser.tile();
        }
      }
    }

    static int decrementToZero(int value) {
      return Math.max(0, value - 1);
    }
  }

  /**
   * Hold the work settings.
   */
  private static class WorkSettings {
    double mean;
    double sigma;
    double qvalue;

    WorkSettings(double mean, double sigma, double qvalue) {
      this.mean = mean;
      this.sigma = sigma;
      this.qvalue = qvalue;
    }
  }

  /**
   * Base worker for interactive analysis.
   */
  private class BaseWorker implements WorkflowWorker<WorkSettings, Object> {
    final MyWindowOrganiser wo;

    BaseWorker(MyWindowOrganiser wo) {
      this.wo = wo;
    }

    @Override
    public boolean equalSettings(WorkSettings current, WorkSettings previous) {
      return (current.mean == previous.mean && current.sigma == previous.sigma);
    }

    @Override
    public boolean equalResults(Object current, Object previous) {
      // We never create any results so ignore this
      return true;
    }

    @Override
    public Pair<WorkSettings, Object> doWork(Pair<WorkSettings, Object> work) {
      return work;
    }
  }

  /**
   * Worker for the histogram.
   */
  private class HistogramWorker extends BaseWorker {
    final PrecisionHistogram histogram;

    HistogramWorker(MyWindowOrganiser wo, PrecisionHistogram histogram) {
      super(wo);
      this.histogram = histogram;
    }

    @Override
    public Pair<WorkSettings, Object> doWork(Pair<WorkSettings, Object> work) {
      // Plot the histogram
      histogram.plot(work.getKey().mean, work.getKey().sigma, wo);
      return work;
    }
  }

  /**
   * Worker for the Q plot.
   */
  private class QPlotWorker extends BaseWorker {
    final QPlot qplot;

    QPlotWorker(MyWindowOrganiser wo, QPlot qplot) {
      super(wo);
      this.qplot = qplot;
    }

    @Override
    public boolean equalSettings(WorkSettings current, WorkSettings previous) {
      if (current.qvalue != previous.qvalue) {
        return false;
      }
      return super.equalSettings(current, previous);
    }

    @Override
    public Pair<WorkSettings, Object> doWork(Pair<WorkSettings, Object> work) {
      // Compute Q and then plot the scaled FRC numerator
      final WorkSettings workSettings = work.getKey();
      qplot.plot(workSettings.mean, workSettings.sigma, workSettings.qvalue, wo);
      return work;
    }
  }

  private boolean showQEstimationDialog(final PrecisionHistogram histogram, final QPlot qplot,
      final double nmPerPixel) {
    // This is used for the initial layout of windows
    final MyWindowOrganiser wo = new MyWindowOrganiser(2);

    // Use a simple workflow
    final Workflow<WorkSettings, Object> workflow = new Workflow<>();

    // Split the work to two children with a dummy initial worker
    final int previous = workflow.add(new BaseWorker(wo));
    workflow.add(new HistogramWorker(wo, histogram), previous);
    workflow.add(new QPlotWorker(wo, qplot), previous);

    workflow.start();

    final String macroOptions = Macro.getOptions();
    if (macroOptions == null) {
      // Draw the plots with the first set of work
      workflow.run(new WorkSettings(histogram.mean, histogram.sigma, qplot.qvalue));

      // Build the dialog
      final NonBlockingExtendedGenericDialog gd = new NonBlockingExtendedGenericDialog(pluginTitle);
      gd.addHelp(HelpUrls.getUrl("fourier-image-resolution"));

      final double mu = histogram.mean / nmPerPixel;
      final double sd = histogram.sigma / nmPerPixel;
      final double plateauness = qplot.computePlateauness(qplot.qvalue, mu, sd);

      ImageJUtils.addMessage(gd,
          "Estimate the blinking correction parameter Q for Fourier Ring Correlation\n \n"
              + "Initial estimate:\nPrecision = %.3f +/- %.3f\nQ = %s\nCost = %.3f",
          histogram.mean, histogram.sigma, MathUtils.rounded(qplot.qvalue), plateauness);

      final double mean10 = histogram.mean * 10;
      final double sd10 = histogram.sigma * 10;
      final double q10 = qplot.qvalue * 10;

      gd.addSlider("Mean (x10)", Math.max(0, mean10 - sd10 * 2), mean10 + sd10 * 2, mean10);
      gd.addSlider("Sigma (x10)", Math.max(0, sd10 / 2), sd10 * 2, sd10);
      gd.addSlider("Q (x10)", 0, Math.max(50, q10 * 2), q10);
      gd.addCheckbox("Reset_all", false);
      gd.addMessage("Double-click a slider to reset");

      gd.addDialogListener(new FireDialogListener(gd, histogram, qplot, workflow));

      // Show this when the workers have finished drawing the plots so it is on top
      try {
        final long timeout = System.currentTimeMillis() + 5000;
        while (!wo.isTiled()) {
          Thread.sleep(50);
          if (System.currentTimeMillis() > timeout) {
            break;
          }
        }
      } catch (final InterruptedException ex) {
        // Propagate the interruption
        Thread.currentThread().interrupt();
      }

      gd.addHelp(HelpUrls.getUrl("fire-q-estimation"));
      gd.showDialog();

      // Finish the worker threads
      final boolean cancelled = gd.wasCanceled();
      workflow.shutdown(cancelled);
      if (cancelled) {
        return false;
      }
    } else {
      // If inside a macro then just get the options and run the work
      final double mean = Double
          .parseDouble(Macro.getValue(macroOptions, KEY_MEAN, Double.toString(histogram.mean)));
      final double sigma = Double
          .parseDouble(Macro.getValue(macroOptions, KEY_SIGMA, Double.toString(histogram.sigma)));
      final double qvalue =
          Double.parseDouble(Macro.getValue(macroOptions, KEY_Q, Double.toString(qplot.qvalue)));
      workflow.run(new WorkSettings(mean, sigma, qvalue));
      workflow.shutdown(false);
    }

    // Store the Q value and the mean and sigma
    settings.qvalue = qplot.qvalue;
    settings.mean = qplot.mean;
    settings.sigma = qplot.sigma;

    // Record the values for Macros since the NonBlockingDialog doesn't
    if (Recorder.record) {
      Recorder.recordOption(KEY_MEAN, Double.toString(settings.mean));
      Recorder.recordOption(KEY_SIGMA, Double.toString(settings.sigma));
      Recorder.recordOption(KEY_Q, Double.toString(settings.qvalue));
    }

    return true;
  }

  /**
   * Listen to dialog changes to queue work for interactive analysis.
   */
  private static class FireDialogListener extends MouseAdapter implements DialogListener {
    /**
     * Delay (in milliseconds) used when entering new values in the dialog before the preview is
     * processed.
     */
    @SuppressWarnings("unused")
    static final long DELAY = 500;

    long time;
    boolean notActive = true;
    volatile int ignore;
    Workflow<WorkSettings, Object> workflow;
    double defaultMean;
    double defaultSigma;
    double defaultQValue;
    String textM;
    String textS;
    String textQ;
    TextField tf1;
    TextField tf2;
    TextField tf3;
    Scrollbar sl1;
    Scrollbar sl2;
    Scrollbar sl3;
    Checkbox cb;
    final boolean isMacro;

    FireDialogListener(ExtendedGenericDialog gd, PrecisionHistogram histogram, QPlot qplot,
        Workflow<WorkSettings, Object> workflow) {
      time = System.currentTimeMillis() + 1000;
      this.workflow = workflow;
      this.defaultMean = histogram.mean;
      this.defaultSigma = histogram.sigma;
      this.defaultQValue = qplot.qvalue;
      isMacro = ImageJUtils.isMacro();
      // For the reset
      tf1 = (TextField) gd.getNumericFields().get(0);
      tf2 = (TextField) gd.getNumericFields().get(1);
      tf3 = (TextField) gd.getNumericFields().get(2);
      cb = (Checkbox) (gd.getCheckboxes().get(0));
      // Sliders
      sl1 = (Scrollbar) gd.getSliders().get(0);
      sl2 = (Scrollbar) gd.getSliders().get(1);
      sl3 = (Scrollbar) gd.getSliders().get(2);
      sl1.addMouseListener(this);
      sl2.addMouseListener(this);
      sl3.addMouseListener(this);
      textM = tf1.getText();
      textS = tf2.getText();
      textQ = tf3.getText();

      // Implement a delay to allow typing.
      // This is also applied to the sliders which we do not want.
      // Ideally we would have no delay for sliders (since they are in the correct place already
      // but a delay for typing in the text field). Unfortunately the AWTEvent raised by ImageJ
      // for the slider is actually from the TextField so we cannot tell the difference.
      // For now just have no delay.
      // if (!isMacro)
      // workflow.startPreview();
    }

    @Override
    public boolean dialogItemChanged(GenericDialog gd, AWTEvent event) {
      // Delay reading the dialog when in interactive mode. This is a workaround for a bug
      // where the dialog has not yet been drawn.
      if ((notActive && !isMacro && System.currentTimeMillis() < time) || (ignore-- > 0)) {
        return true;
      }

      notActive = false;

      double mean = Math.abs(gd.getNextNumber()) / 10;
      double sigma = Math.abs(gd.getNextNumber()) / 10;
      double qvalue = Math.abs(gd.getNextNumber()) / 10;
      final boolean reset = gd.getNextBoolean();

      // Even events from the slider come through as TextEvent from the TextField
      // since ImageJ captures the slider event as just updates the TextField.

      // Allow reset to default
      if (reset) {
        // This does not trigger the event
        cb.setState(false);
        mean = this.defaultMean;
        sigma = this.defaultSigma;
        qvalue = this.defaultQValue;
      }

      final WorkSettings work = new WorkSettings(mean, sigma, qvalue);

      // Offload this work onto a thread that just picks up the most recent dialog input.
      workflow.run(work);

      if (reset) {
        // These trigger dialogItemChanged(...) so do them after we added
        // work to the queue and ignore the events
        ignore = 3;
        tf1.setText(textM);
        tf2.setText(textS);
        tf3.setText(textQ);
      }

      return true;
    }

    @Override
    public void mouseClicked(MouseEvent event) {
      // Reset the slider on double-click
      if ((event.getClickCount() <= 1) || !(event.getSource() instanceof Scrollbar)) {
        return;
      }
      final Scrollbar sl = (Scrollbar) event.getSource();
      if (sl == sl1) {
        tf1.setText(textM);
      }
      if (sl == sl2) {
        tf2.setText(textS);
      }
      if (sl == sl3) {
        tf3.setText(textQ);
      }
    }
  }

  /**
   * Gets the last N threads used in the input dialog.
   *
   * @return the last N threads
   */
  private static synchronized int getLastNThreads() {
    // See if ImageJ preference were updated
    if (imagejNThreads != Prefs.getThreads()) {
      lastNThreads = imagejNThreads = Prefs.getThreads();
    }
    // Otherwise use the last user input
    return lastNThreads;
  }

  /**
   * SGets the last N threads used in the input dialog.
   *
   * @param numberOfThreads the new last N threads
   */
  private static synchronized void setLastNThreads(int numberOfThreads) {
    lastNThreads = numberOfThreads;
  }

  /**
   * Gets the threads to use for multi-threaded computation.
   *
   * @return the threads
   */
  private int getThreads() {
    if (numberOfThreads == 0) {
      numberOfThreads = Prefs.getThreads();
    }
    return numberOfThreads;
  }

  /**
   * Sets the threads to use for multi-threaded computation.
   *
   * @param threads the new threads
   */
  public void setThreads(int threads) {
    this.numberOfThreads = Math.max(1, threads);
  }

  private void setProgress(int repeats) {
    progress = new ParallelTrackProgress(repeats);
  }
}
