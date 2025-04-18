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
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.ImageRoi;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.Roi;
import ij.gui.RoiListener;
import ij.io.FileSaver;
import ij.io.Opener;
import ij.measure.Measurements;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import java.awt.Color;
import java.awt.Label;
import java.awt.Rectangle;
import java.awt.TextField;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collections;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.Future;
import java.util.concurrent.RejectedExecutionException;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Consumer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.swing.SwingUtilities;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.concurrent.ConcurrentRuntimeException;
import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.ContinuousSampler;
import org.apache.commons.rng.sampling.distribution.DiscreteSampler;
import org.apache.commons.rng.sampling.distribution.NormalizedGaussianSampler;
import org.apache.commons.rng.sampling.distribution.SharedStateContinuousSampler;
import org.apache.commons.statistics.inference.KolmogorovSmirnovTest;
import org.apache.commons.statistics.inference.SignificanceResult;
import org.apache.commons.statistics.inference.TTest;
import org.apache.commons.statistics.inference.WilcoxonSignedRankTest;
import uk.ac.sussex.gdsc.core.data.IntegerType;
import uk.ac.sussex.gdsc.core.data.SiPrefix;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot.HistogramPlotBuilder;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.SimpleImageJTrackProgress;
import uk.ac.sussex.gdsc.core.ij.gui.NonBlockingExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.io.CustomTiffEncoder;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.ij.roi.RoiHelper;
import uk.ac.sussex.gdsc.core.logging.Ticker;
import uk.ac.sussex.gdsc.core.logging.TrackProgress;
import uk.ac.sussex.gdsc.core.math.ArrayMoment;
import uk.ac.sussex.gdsc.core.math.IntegerArrayMoment;
import uk.ac.sussex.gdsc.core.math.RollingArrayMoment;
import uk.ac.sussex.gdsc.core.math.SimpleArrayMoment;
import uk.ac.sussex.gdsc.core.utils.DoubleData;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.SoftLock;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.concurrent.ConcurrencyUtils;
import uk.ac.sussex.gdsc.core.utils.rng.Pcg32;
import uk.ac.sussex.gdsc.core.utils.rng.PoissonSamplers;
import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;
import uk.ac.sussex.gdsc.core.utils.rng.UniformRandomProviders;
import uk.ac.sussex.gdsc.smlm.ij.SeriesImageSource;
import uk.ac.sussex.gdsc.smlm.ij.utils.ImageJImageConverter;
import uk.ac.sussex.gdsc.smlm.model.camera.PerPixelCameraModel;

/**
 * Analyse the per pixel offset, variance and gain from a sCMOS camera.
 *
 * <p>See Huang et al (2013) Video-rate nanoscopy using sCMOS camera–specific single-molecule
 * localization algorithms. Nature Methods 10, 653-658 (Supplementary Information).
 */
public class CmosAnalysis implements PlugIn {
  private static final String TITLE = "sCMOS Analysis";
  private static final String REJECT = "reject";
  private static final String ACCEPT = "accept";

  private int numberOfThreads;
  // The simulated offset, variance and gain
  private ImagePlus simulationImp;
  // The measured offset, variance and gain
  private ImageStack measuredStack;
  // The sub-directories containing the sCMOS images
  private LocalList<SubDir> subDirs;

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    private static final String DEFAULT_PHOTONS = "50, 100, 200, 400, 800";
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> INSTANCE = new AtomicReference<>(new Settings());

    String directory;
    String modelDirectory;
    String modelName;
    int originX;
    int originY;
    int width;
    int height;
    boolean rollingAlgorithm;
    boolean reuseProcessedData;
    boolean logFileProgress;
    int histogramBins;
    double offset;
    double variance;
    double gain;
    double gainStdDev;
    int size;
    int frames;
    int imagejNThreads;
    int lastNumberOfThreads;
    String simulationPhotons;

    Settings() {
      // Set defaults
      directory = Prefs.get(PrefsKey.sCMOSAnalysisDirectory, "");
      modelDirectory = Prefs.get(PrefsKey.sCMOSModelDirectory, "");
      modelName = "sCMOS Camera";
      reuseProcessedData = true;
      histogramBins = 200;

      // The simulation can default roughly to the values displayed
      // in the Huang sCMOS paper supplementary figure 1:

      // Offset = Approximately Normal or Poisson. We use Poisson
      // since that is an integer distribution which would be expected
      // for an offset & Poisson approaches the Gaussian at high mean.
      offset = 100;

      // Variance = Exponential (equivalent to chi-squared with k=1, i.e.
      // sum of the squares of 1 normal distribution).
      // We want 99.9% @ 400 ADU based on supplementary figure 1.a/1.b
      // cumul = 1 - e^-lx (l = 1/mean)
      // => e^-lx = 1 - cumul
      // => -lx = log(1-0.999)
      // => l = -log(0.001) / 400 (since x==400)
      // => 1/l = 57.9
      variance = 57.9; // SD = 7.6

      // Gain = Approximately Normal
      gain = 2.2;
      gainStdDev = 0.2;

      size = 512;
      frames = 500;

      imagejNThreads = Prefs.getThreads();
      lastNumberOfThreads = imagejNThreads;
      simulationPhotons = DEFAULT_PHOTONS;
    }

    Settings(Settings source) {
      directory = source.directory;
      modelDirectory = source.modelDirectory;
      modelName = source.modelName;
      originX = source.originX;
      originY = source.originY;
      width = source.width;
      height = source.height;
      rollingAlgorithm = source.rollingAlgorithm;
      reuseProcessedData = source.reuseProcessedData;
      logFileProgress = source.logFileProgress;
      histogramBins = source.histogramBins;
      offset = source.offset;
      variance = source.variance;
      gain = source.gain;
      gainStdDev = source.gainStdDev;
      size = source.size;
      frames = source.frames;
      imagejNThreads = source.imagejNThreads;
      lastNumberOfThreads = source.lastNumberOfThreads;
      simulationPhotons = source.simulationPhotons;
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
      Prefs.set(PrefsKey.sCMOSAnalysisDirectory, directory);
      Prefs.set(PrefsKey.sCMOSModelDirectory, modelDirectory);
    }

    /**
     * Gets the photons from the simulation photons string.
     *
     * @return the photons
     */
    int[] getPhotons() {
      final IntOpenHashSet list = new IntOpenHashSet();
      list.add(0);
      for (final String key : simulationPhotons.split(",")) {
        try {
          final int value = Integer.parseInt(key.trim());
          if (value < 0 || (value * gain) > (1 << 16)) {
            throw new NumberFormatException("Invalid number of photons after gain: " + value);
          }
          list.add(value);
        } catch (final NumberFormatException ex) {
          final NumberFormatException ex2 =
              new NumberFormatException("Invalid number of photons: " + key);
          ex2.initCause(ex);
          throw ex;
        }
      }
      final int[] result = list.toIntArray();
      Arrays.sort(result);
      return result;
    }
  }

  /**
   * Class to compute a sCMOS simulated image.
   */
  private static class SimulationWorker implements Runnable {
    /** The maximum value for a short integer. */
    static final short MAX_SHORT = (short) 0xffff;

    final Ticker ticker;
    final UniformRandomProvider rg;
    final String out;
    final float[] pixelOffset;
    final float[] pixelVariance;
    final float[] pixelGain;
    final int from;
    final int to;
    final int blockSize;
    final int photons;
    final NormalizedGaussianSampler gauss;

    SimulationWorker(Ticker ticker, UniformRandomProvider rng, String out, ImageStack stack,
        int from, int to, int blockSize, int photons) {
      this.ticker = ticker;
      rg = rng;
      pixelOffset = (float[]) stack.getPixels(1);
      pixelVariance = (float[]) stack.getPixels(2);
      pixelGain = (float[]) stack.getPixels(3);
      this.out = out;
      this.from = from;
      this.to = to;
      this.blockSize = blockSize;
      this.photons = photons;
      gauss = SamplerUtils.createNormalizedGaussianSampler(rg);
    }

    @Override
    public void run() {
      // Avoid the status bar talking to the current image
      WindowManager.setTempCurrentImage(null);

      // Convert variance to SD
      final float[] pixelSd = new float[pixelVariance.length];
      for (int i = 0; i < pixelVariance.length; i++) {
        pixelSd[i] = (float) Math.sqrt(pixelVariance[i]);
      }

      final int size = (int) Math.sqrt(pixelVariance.length);

      // Pre-compute a set of Poisson numbers since this is slow
      final DiscreteSampler pd = PoissonSamplers.createPoissonSampler(rg, photons);
      // // For speed we can precompute a set of random numbers to reuse
      // final int[] poisson = new int[pixelVariance.length];
      // for (int i = poisson.length; i-- > 0;) {
      // poisson[i] = pd.sample();
      // }

      // Save image in blocks
      ImageStack stack = new ImageStack(size, size);
      int start = from;
      for (int i = from; i < to; i++) {
        // Create image
        final short[] pixels = new short[pixelOffset.length];
        if (photons == 0) {
          for (int j = 0; j < pixelOffset.length; j++) {
            // Fixed offset per pixel plus a variance
            final double p = pixelOffset[j] + gauss.sample() * pixelSd[j];
            pixels[j] = clip16bit(p);
          }
        } else {
          for (int j = 0; j < pixelOffset.length; j++) {
            // Fixed offset per pixel plus a variance plus a
            // fixed gain multiplied by a Poisson sample of the photons
            final double p =
                pixelOffset[j] + gauss.sample() * pixelSd[j] + (pd.sample() * pixelGain[j]);
            pixels[j] = clip16bit(p);
          }

          // Rotate Poisson numbers.
          // Shuffling what we have is faster than generating new values
          // and we should have enough.
          // RandomUtils.shuffle(poisson, rg);
        }

        // Save image
        stack.addSlice(null, pixels);
        if (stack.getSize() == blockSize) {
          save(stack, start);
          start = i + 1;
          stack = new ImageStack(size, size);
        }

        ticker.tick();
      }
      // This should not happen if we control the to-from range correctly
      if (stack.getSize() != 0) {
        save(stack, start);
      }
    }

    /**
     * Clip to the range for a 16-bit image.
     *
     * @param value the value
     * @return the clipped value
     */
    private static short clip16bit(double value) {
      if (value < 0) {
        return 0;
      }
      final int i = (int) Math.round(value);
      if (i > 0xffff) {
        return MAX_SHORT;
      }
      return (short) i;
    }

    private void save(ImageStack stack, int start) {
      final ImagePlus imp = new ImagePlus("", stack);
      final CustomTiffEncoder file = new CustomTiffEncoder(imp.getFileInfo());
      try (
          OutputStream os =
              Files.newOutputStream(Paths.get(out, String.format("image%06d.tif", start)));
          BufferedOutputStream bos = new BufferedOutputStream(os)) {
        file.write(bos);
      } catch (final IOException ex) {
        throw new ConcurrentRuntimeException("Failed to save image", ex);
      }
    }
  }

  /**
   * Sub-directory containing images with a specified exposure time.
   */
  private static class SubDir {
    int exposureTime;
    File path;
    String name;

    SubDir(int exposureTime, File path, String name) {
      this.exposureTime = exposureTime;
      this.path = path;
      this.name = name;
    }

    /**
     * Compare the two results.
     *
     * @param r1 the first result
     * @param r2 the second result
     * @return -1, 0 or 1
     */
    static int compare(SubDir r1, SubDir r2) {
      return Integer.compare(r1.exposureTime, r2.exposureTime);
    }
  }

  /**
   * Used to allow multi-threading of the scoring the filters.
   */
  private static class ImageWorker implements Runnable {
    /** The signal to stop processing jobs. */
    static final Object STOP_SIGNAL = new Object();
    final Ticker ticker;
    volatile boolean finished;
    final BlockingQueue<Object> jobs;
    final ArrayMoment moment;
    int bitDepth;

    ImageWorker(Ticker ticker, BlockingQueue<Object> jobs, ArrayMoment moment) {
      this.ticker = ticker;
      this.jobs = jobs;
      this.moment = moment.newInstance();
    }

    @Override
    public void run() {
      try {
        for (;;) {
          final Object pixels = jobs.take();
          if (pixels == STOP_SIGNAL) {
            break;
          }
          if (!finished) {
            // Only run jobs when not finished. This allows the queue to be emptied.
            run(pixels);
          }
        }
      } catch (final InterruptedException ex) {
        ConcurrencyUtils.interruptAndThrowUncheckedIf(!finished, ex);
      } finally {
        finished = true;
      }
    }

    private void run(Object pixels) {
      if (ImageJUtils.isInterrupted()) {
        finished = true;
        return;
      }
      if (bitDepth == 0) {
        bitDepth = ImageJUtils.getBitDepth(pixels);
      }
      // Most likely first
      if (bitDepth == 16) {
        moment.addUnsigned((short[]) pixels);
      } else if (bitDepth == 32) {
        moment.add((float[]) pixels);
      } else if (bitDepth == 8) {
        moment.addUnsigned((byte[]) pixels);
      } else {
        throw new IllegalStateException("Unsupported bit depth");
      }
      ticker.tick();
    }
  }

  /**
   * Gets the last number of threads used in the input dialog.
   *
   * @return the last number of threads
   */
  private int getLastNumberOfThreads() {
    // See if ImageJ preference were updated
    if (settings.imagejNThreads != Prefs.getThreads()) {
      settings.lastNumberOfThreads = settings.imagejNThreads = Prefs.getThreads();
    }
    // Otherwise use the last user input
    return settings.lastNumberOfThreads;
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
   * @param numberOfThreads the new threads
   */
  private void setThreads(int numberOfThreads) {
    this.numberOfThreads = Math.max(1, numberOfThreads);
    // Save user input
    settings.lastNumberOfThreads = this.numberOfThreads;
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    final boolean extraOptions = ImageJUtils.isExtraOptions();
    // Avoid the status bar talking to the current image
    WindowManager.setTempCurrentImage(null);

    //@formatter:off
    IJ.log(TextUtils.wrap(
        TITLE + ": Analyse the per-pixel offset, variance and gain of sCMOS images. " +
        "See Huang et al (2013) Video-rate nanoscopy using sCMOS camera–specific " +
        "single-molecule localization algorithms. Nature Methods 10, 653-658 " +
        "(Supplementary Information).",
        80));
    //@formatter:on

    settings = Settings.load();
    final String dir = ImageJUtils.getDirectory(TITLE, settings.directory);
    if (TextUtils.isNullOrEmpty(dir)) {
      return;
    }
    // Reset the crop if the analysis directory changes
    if (!dir.equals(settings.directory)) {
      settings.originX = settings.originY = settings.width = settings.height = 0;
    }
    settings.directory = dir;
    settings.save();

    final boolean simulate = "simulate".equals(arg);
    if (simulate || extraOptions) {
      if (!showSimulateDialog()) {
        return;
      }
      try {
        simulate();
      } catch (final IOException ex) {
        IJ.error(TITLE, "Failed to perform simulation: " + ex.getMessage());
        return;
      }
    }

    if (!showDialog()) {
      return;
    }

    runAnalysis();

    if (simulationImp == null) {
      // Just in case an old simulation is in the directory
      final ImagePlus imp =
          IJ.openImage(new File(settings.directory, "perPixelSimulation.tif").getPath());
      if (imp != null && imp.getStackSize() == 3 && imp.getWidth() == measuredStack.getWidth()
          && imp.getHeight() == measuredStack.getHeight()) {
        simulationImp = imp;
      }
    }

    if (simulationImp != null) {
      computeError();
    }
  }

  private void simulate() throws IOException {
    // Create the offset, variance and gain for each pixel
    final int n = settings.size * settings.size;
    final float[] pixelOffset = new float[n];
    final float[] pixelVariance = new float[n];
    final float[] pixelGain = new float[n];

    IJ.showStatus("Creating random per-pixel readout");
    final long start = System.currentTimeMillis();

    final UniformRandomProvider rg = UniformRandomProviders.create();

    final DiscreteSampler pd = PoissonSamplers.createPoissonSampler(rg, settings.offset);
    final ContinuousSampler ed = SamplerUtils.createExponentialSampler(rg, settings.variance);
    final SharedStateContinuousSampler gauss =
        SamplerUtils.createGaussianSampler(rg, settings.gain, settings.gainStdDev);
    Ticker ticker = ImageJUtils.createTicker(n, 0);
    for (int i = 0; i < n; i++) {
      // Q. Should these be clipped to a sensible range?
      pixelOffset[i] = pd.sample();
      pixelVariance[i] = (float) ed.sample();
      pixelGain[i] = (float) gauss.sample();
      ticker.tick();
    }
    IJ.showProgress(1);

    // Save to the directory as a stack
    final ImageStack simulationStack = new ImageStack(settings.size, settings.size);
    simulationStack.addSlice("Offset", pixelOffset);
    simulationStack.addSlice("Variance", pixelVariance);
    simulationStack.addSlice("Gain", pixelGain);
    simulationImp = new ImagePlus("PerPixel", simulationStack);
    // Only the info property is saved to the TIFF file
    simulationImp.setProperty("Info",
        String.format("Offset=%s; Variance=%s; Gain=%s +/- %s", MathUtils.rounded(settings.offset),
            MathUtils.rounded(settings.variance), MathUtils.rounded(settings.gain),
            MathUtils.rounded(settings.gainStdDev)));
    IJ.save(simulationImp, new File(settings.directory, "perPixelSimulation.tif").getPath());

    // Create thread pool and workers
    final int threadCount = getThreads();
    final ExecutorService executor = Executors.newFixedThreadPool(threadCount);
    final LocalList<Future<?>> futures = new LocalList<>(numberOfThreads);

    // Simulate the exposure input.
    final int[] photons = settings.getPhotons();

    final int blockSize = 10; // For saving stacks
    int numberPerThread = (int) Math.ceil((double) settings.frames / numberOfThreads);
    // Convert to fit the block size
    numberPerThread = (int) Math.ceil((double) numberPerThread / blockSize) * blockSize;
    final Pcg32 rng = Pcg32.xshrs(start);
    // Note the bias is increased by 3-fold so add 2 to the length
    ticker = ImageJUtils.createTicker((long) (photons.length + 2) * settings.frames, threadCount);
    for (final int p : photons) {
      ImageJUtils.showStatus(() -> "Simulating " + TextUtils.pleural(p, "photon"));

      // Create the directory
      final Path out = Paths.get(settings.directory, String.format("photon%03d", p));
      Files.createDirectories(out);

      // Increase frames for bias image
      final int frames = settings.frames * (p == 0 ? 3 : 1);
      for (int from = 0; from < frames;) {
        final int to = Math.min(from + numberPerThread, frames);
        futures.add(executor.submit(new SimulationWorker(ticker, rng.split(), out.toString(),
            simulationStack, from, to, blockSize, p)));
        from = to;
      }

      ConcurrencyUtils.waitForCompletionUnchecked(futures);
      futures.clear();
    }

    final String msg =
        "Simulation time = " + TextUtils.millisToString(System.currentTimeMillis() - start);

    IJ.showStatus(msg);
    ImageJUtils.clearSlowProgress();

    executor.shutdown();

    ImageJUtils.log(msg);
  }

  private boolean showSimulateDialog() {
    final GenericDialog gd = new GenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("scmos-analysis"));

    gd.addMessage("Simulate per-pixel offset, variance and gain of sCMOS images.");

    gd.addNumericField("nThreads", getLastNumberOfThreads(), 0);
    gd.addNumericField("Offset (Poisson)", settings.offset, 3);
    gd.addNumericField("Variance (Exponential)", settings.variance, 3);
    gd.addNumericField("Gain (Gaussian)", settings.gain, 3);
    gd.addNumericField("Gain_SD", settings.gainStdDev, 3);
    gd.addNumericField("Size", settings.size, 0);
    gd.addNumericField("Frames", settings.frames, 0);
    gd.addStringField("Photons (comma-delimited)", settings.simulationPhotons, 20);
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }

    setThreads((int) gd.getNextNumber());
    settings.offset = Math.abs(gd.getNextNumber());
    settings.variance = Math.abs(gd.getNextNumber());
    settings.gain = Math.abs(gd.getNextNumber());
    settings.gainStdDev = Math.abs(gd.getNextNumber());
    settings.size = Math.abs((int) gd.getNextNumber());
    settings.frames = Math.abs((int) gd.getNextNumber());
    settings.simulationPhotons = gd.getNextString();

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Offset", settings.offset);
      ParameterUtils.isAboveZero("Variance", settings.variance);
      ParameterUtils.isAboveZero("Gain", settings.gain);
      ParameterUtils.isAboveZero("Gain SD", settings.gainStdDev);
      ParameterUtils.isAboveZero("Size", settings.size);
      ParameterUtils.isAboveZero("Frames", settings.frames);
    } catch (final IllegalArgumentException ex) {
      ImageJUtils.log(TITLE + ": " + ex.getMessage());
      return false;
    }

    return true;
  }

  private boolean showDialog() {
    // Determine sub-directories to process
    final File dir = new File(settings.directory);
    final File[] dirs = dir.listFiles(File::isDirectory);

    if (ArrayUtils.isEmpty(dirs)) {
      IJ.error(TITLE, "No sub-directories");
      return false;
    }

    // Get only those with numbers at the end.
    // These should correspond to exposure times
    subDirs = new LocalList<>();
    final Pattern p = Pattern.compile("([0-9]+)$");
    for (final File path : dirs) {
      final String name = path.getName();
      final Matcher m = p.matcher(name);
      if (m.find()) {
        final int t = Integer.parseInt(m.group(1));
        subDirs.add(new SubDir(t, path, name));
      }
    }

    if (subDirs.size() <= 1) {
      IJ.error(TITLE, "Not enough sub-directories with exposure time suffix");
      return false;
    }

    Collections.sort(subDirs, SubDir::compare);

    if (subDirs.get(0).exposureTime != 0) {
      IJ.error(TITLE, "No sub-directories with exposure time 0");
      return false;
    }

    for (final SubDir sd : subDirs) {
      ImageJUtils.log("Sub-directory: %s. Exposure time = %d", sd.name, sd.exposureTime);
    }

    final GenericDialog gd = new GenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("scmos-analysis"));

    //@formatter:off
    gd.addMessage("Analyse the per-pixel offset, variance and gain of sCMOS images.\n \n" +
        TextUtils.wrap(
        "See Huang et al (2013) Video-rate nanoscopy using sCMOS camera–specific " +
        "single-molecule localization algorithms. Nature Methods 10, 653-658 " +
        "(Supplementary Information).",
        80));
    //@formatter:on

    gd.addNumericField("nThreads", getLastNumberOfThreads(), 0);
    gd.addMessage(TextUtils.wrap("A rolling algorithm can handle any size of data but is slower. "
        + "Otherwise the camera is assumed to produce a maximum of 16-bit unsigned data.", 80));
    gd.addCheckbox("Rolling_algorithm", settings.rollingAlgorithm);
    gd.addCheckbox("Re-use_processed_data", settings.reuseProcessedData);
    gd.addCheckbox("Log_file_progress", settings.logFileProgress);
    gd.addNumericField("Histogram_bins", settings.histogramBins);
    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    setThreads((int) gd.getNextNumber());
    settings.rollingAlgorithm = gd.getNextBoolean();
    settings.reuseProcessedData = gd.getNextBoolean();
    settings.logFileProgress = gd.getNextBoolean();
    settings.histogramBins = (int) gd.getNextNumber();

    return true;
  }

  private void runAnalysis() {
    final long start = System.currentTimeMillis();

    // Create thread pool and workers. The system is likely to be IO limited
    // so reduce the computation threads to allow the reading thread in the
    // SeriesImageSource to run.
    // If the images are small enough to fit into memory then 3 threads are used,
    // otherwise it is 1.
    final int nThreads = Math.max(1, getThreads() - 3);
    final ExecutorService executor = Executors.newFixedThreadPool(nThreads);
    final LocalList<Future<?>> futures = new LocalList<>(nThreads);
    final LocalList<ImageWorker> workers = new LocalList<>(nThreads);

    final double[][] data = new double[subDirs.size() * 2][];
    double[] pixelOffset = null;
    double[] pixelVariance = null;
    Statistics statsOffset = null;
    Statistics statsVariance = null;

    // For each sub-directory compute the mean and variance
    final int nSubDirs = subDirs.size();
    boolean error = false;
    int width = 0;
    int height = 0;
    final TrackProgress trackProgress =
        settings.logFileProgress ? SimpleImageJTrackProgress.getInstance() : null;
    ImagePlus lastMeanVar = null;

    for (int n = 0; n < nSubDirs; n++) {
      ImageJUtils.showSlowProgress(n, nSubDirs);

      final SubDir sd = subDirs.unsafeGet(n);
      ImageJUtils.showStatus(() -> "Analysing " + sd.name);
      final StopWatch sw = StopWatch.createStarted();

      final File file = new File(settings.directory, "perPixel" + sd.name + ".tif");
      boolean found = false;
      if (settings.reuseProcessedData && file.exists()) {
        final Opener opener = new Opener();
        opener.setSilentMode(true);
        final ImagePlus imp = opener.openImage(file.getPath());
        if (imp != null && imp.getStackSize() == 2 && imp.getBitDepth() == 32) {
          if (n == 0) {
            width = imp.getWidth();
            height = imp.getHeight();
          } else if (width != imp.getWidth() || height != imp.getHeight()) {
            error = true;
            IJ.error(TITLE,
                "Image width/height mismatch in image series: " + file.getPath()
                    + String.format("\n \nExpected %dx%d, Found %dx%d", width, height,
                        imp.getWidth(), imp.getHeight()));
            break;
          }

          final ImageStack stack = imp.getImageStack();
          data[2 * n] = SimpleArrayUtils.toDouble((float[]) stack.getPixels(1));
          data[2 * n + 1] = SimpleArrayUtils.toDouble((float[]) stack.getPixels(2));
          found = true;
          lastMeanVar = imp;
        }
      }

      if (!found) {
        // Open the series
        final SeriesImageSource source = new SeriesImageSource(sd.name, sd.path.getPath());
        source.setTrackProgress(trackProgress);
        if (!source.open()) {
          error = true;
          IJ.error(TITLE, "Failed to open image series: " + sd.path.getPath());
          break;
        }

        if (n == 0) {
          width = source.getWidth();
          height = source.getHeight();
        } else if (width != source.getWidth() || height != source.getHeight()) {
          error = true;
          IJ.error(TITLE,
              "Image width/height mismatch in image series: " + sd.path.getPath()
                  + String.format("\n \nExpected %dx%d, Found %dx%d", width, height,
                      source.getWidth(), source.getHeight()));
          break;
        }

        // So the bar remains at 99% when workers have finished use frames + 1
        final Ticker ticker = ImageJUtils.createTicker(source.getFrames() + 1L, nThreads);

        // Open the first frame to get the bit depth.
        // Assume the first pixels are not empty as the source is open.
        Object pixels = source.nextRaw();
        final int bitDepth = ImageJUtils.getBitDepth(pixels);

        ArrayMoment moment;
        if (settings.rollingAlgorithm) {
          moment = new RollingArrayMoment();
          // We assume 16-bit camera at the maximum
        } else if (bitDepth <= 16
            && IntegerArrayMoment.isValid(IntegerType.UNSIGNED_16, source.getFrames())) {
          moment = new IntegerArrayMoment();
        } else {
          moment = new SimpleArrayMoment();
        }

        final BlockingQueue<Object> jobs = new ArrayBlockingQueue<>(nThreads * 2);
        for (int i = 0; i < nThreads; i++) {
          final ImageWorker worker = new ImageWorker(ticker, jobs, moment);
          workers.add(worker);
          futures.add(executor.submit(worker));
        }

        // Process the raw pixel data
        long lastTime = 0;
        while (pixels != null) {
          final long time = System.currentTimeMillis();
          if (time - lastTime > 150) {
            if (ImageJUtils.isInterrupted()) {
              error = true;
              break;
            }
            lastTime = time;
            IJ.showStatus("Analysing " + sd.name + " Frame " + source.getStartFrameNumber());
          }
          put(jobs, pixels);
          pixels = source.nextRaw();
        }
        source.close();

        if (error) {
          // Kill the workers
          workers.forEach(worker -> worker.finished = true);
          // Clear the queue
          jobs.clear();
          // Signal any waiting workers
          workers.forEach(worker -> jobs.add(ImageWorker.STOP_SIGNAL));
          // Cancel by interruption. We set the finished flag so the ImageWorker should
          // ignore the interrupt.
          futures.forEach(future -> future.cancel(true));
          break;
        }

        // Finish all the worker threads cleanly
        workers.forEach(worker -> jobs.add(ImageWorker.STOP_SIGNAL));

        // Wait for all to finish
        ConcurrencyUtils.waitForCompletionUnchecked(futures);

        // Create the final aggregate statistics
        for (final ImageWorker w : workers) {
          moment.add(w.moment);
        }
        data[2 * n] = moment.getMean();
        data[2 * n + 1] = moment.getVariance();

        // Get the processing speed.
        sw.stop();
        // ticker holds the number of number of frames processed
        final double bits =
            (double) bitDepth * source.getFrames() * source.getWidth() * source.getHeight();
        final double bps = bits / sw.getTime(TimeUnit.SECONDS);
        final SiPrefix prefix = SiPrefix.getSiPrefix(bps);
        ImageJUtils.log("Processed %d frames. Time = %s. Rate = %s %sbits/s", moment.getN(),
            sw.toString(), MathUtils.rounded(prefix.convert(bps)), prefix.getPrefix());

        // Reset
        futures.clear();
        workers.clear();

        final ImageStack stack = new ImageStack(width, height);
        stack.addSlice("Mean", SimpleArrayUtils.toFloat(data[2 * n]));
        stack.addSlice("Variance", SimpleArrayUtils.toFloat(data[2 * n + 1]));
        lastMeanVar = new ImagePlus("PerPixel", stack);
        IJ.save(lastMeanVar, file.getPath());
      }

      final Statistics s = Statistics.create(data[2 * n]);

      if (pixelOffset == null) {
        // Set the offset assuming the first sub-directory is the bias image
        pixelOffset = data[0];
        pixelVariance = data[1];
        statsOffset = s;
        statsVariance = Statistics.create(pixelVariance);
        ImageJUtils.log("%s Offset = %s +/- %s. Variance = %s +/- %s", sd.name,
            MathUtils.rounded(s.getMean()), MathUtils.rounded(s.getStandardDeviation()),
            MathUtils.rounded(statsVariance.getMean()),
            MathUtils.rounded(statsVariance.getStandardDeviation()));
      } else {
        // Compute mean ADU
        final Statistics signal = new Statistics();
        final double[] mean = data[2 * n];
        for (int i = 0; i < pixelOffset.length; i++) {
          signal.add(mean[i] - pixelOffset[i]);
        }
        ImageJUtils.log("%s Mean = %s +/- %s. Signal = %s +/- %s ADU", sd.name,
            MathUtils.rounded(s.getMean()), MathUtils.rounded(s.getStandardDeviation()),
            MathUtils.rounded(signal.getMean()), MathUtils.rounded(signal.getStandardDeviation()));
      }

      IJ.showProgress(1);
    }
    ImageJUtils.clearSlowProgress();

    if (error) {
      executor.shutdownNow();
      IJ.showStatus(TITLE + " cancelled");
      return;
    }

    executor.shutdown();

    if (pixelOffset == null || pixelVariance == null) {
      IJ.showStatus(TITLE + " error: no bias image");
      return;
    }

    // Compute the gain
    ImageJUtils.showStatus("Computing gain");

    final double[] pixelGain = new double[pixelOffset.length];
    final double[] bibiT = new double[pixelGain.length];
    final double[] biaiT = new double[pixelGain.length];

    // Ignore first as this is the 0 exposure image
    for (int n = 1; n < nSubDirs; n++) {
      // Use equation 2.5 from the Huang et al paper.
      final double[] b = data[2 * n];
      final double[] a = data[2 * n + 1];
      for (int i = 0; i < pixelGain.length; i++) {
        final double bi = b[i] - pixelOffset[i];
        final double ai = a[i] - pixelVariance[i];
        bibiT[i] += bi * bi;
        biaiT[i] += bi * ai;
      }
    }
    for (int i = 0; i < pixelGain.length; i++) {
      pixelGain[i] = biaiT[i] / bibiT[i];
    }

    final Statistics statsGain = Statistics.create(pixelGain);
    ImageJUtils.log("Gain Mean = %s +/- %s", MathUtils.rounded(statsGain.getMean()),
        MathUtils.rounded(statsGain.getStandardDeviation()));

    // Histogram of offset, variance and gain
    final int bins = settings.histogramBins;
    final WindowOrganiser wo = new WindowOrganiser();
    showHistogram("Offset (ADU)", pixelOffset, bins, statsOffset, wo);
    showHistogram("Variance (ADU^^2^^)", pixelVariance, bins, statsVariance, wo);
    showHistogram("Gain (ADU/e)", pixelGain, bins, statsGain, wo);
    wo.tile();

    // Save
    final float[] bias = SimpleArrayUtils.toFloat(pixelOffset);
    final float[] variance = SimpleArrayUtils.toFloat(pixelVariance);
    final float[] gain = SimpleArrayUtils.toFloat(pixelGain);
    measuredStack = new ImageStack(width, height);
    measuredStack.addSlice("Offset", bias);
    measuredStack.addSlice("Variance", variance);
    measuredStack.addSlice("Gain", gain);

    final byte[] badPixels = findBadPixels(bias, gain, variance);

    ImageJUtils
        .log("Analysis time = " + TextUtils.millisToString(System.currentTimeMillis() - start));

    final NonBlockingExtendedGenericDialog egd = new NonBlockingExtendedGenericDialog(TITLE);
    egd.addHelp(HelpUrls.getUrl("scmos-analysis"));
    egd.addMessage("Save the sCMOS camera model?");
    Rectangle reset;
    if (badPixels != null) {
      egd.addMessage("Note: The model contains bad pixels and a crop is required.\n"
          + "Use the crop button to select the crop region.");
      reset = RoiHelper.largestBoundingRoi(new ByteProcessor(width, height, badPixels.clone()),
          width / 2, height / 2).getBounds();
    } else {
      reset = new Rectangle(width, height);
    }
    if (TextUtils.isNullOrEmpty(settings.modelDirectory)) {
      settings.modelDirectory = settings.directory;
    }
    egd.addStringField("Model_name", settings.modelName, 30);
    egd.addDirectoryField("Model_directory", settings.modelDirectory);
    // Configure crop
    final int ox = MathUtils.clip(0, width - 1, settings.originX);
    final int oy = MathUtils.clip(0, height - 1, settings.originY);
    final int w = MathUtils.clip(1, width - ox, settings.width <= 0 ? width : settings.width);
    final int h = MathUtils.clip(1, height - oy, settings.height <= 0 ? height : settings.height);
    final Rectangle bounds = reset.intersection(new Rectangle(ox, oy, w, h));
    final TextField textOriginX = egd.addAndGetNumericField("Origin_x", bounds.x, 0);
    final TextField textOriginY = egd.addAndGetNumericField("Origin_y", bounds.y, 0);
    final TextField textWidth = egd.addAndGetNumericField("Width", bounds.width, 0);
    final TextField textHeight = egd.addAndGetNumericField("Height", bounds.height, 0);
    if (ImageJUtils.isShowGenericDialog()) {
      final Consumer<Rectangle> action = r -> {
        textOriginX.setText(Integer.toString(r.x));
        textOriginY.setText(Integer.toString(r.y));
        textWidth.setText(Integer.toString(r.width));
        textHeight.setText(Integer.toString(r.height));
      };
      final ImagePlus cropImp = lastMeanVar;
      // Do not perform the work on the event dispatcher thread
      egd.addAndGetButton("Crop", e -> ForkJoinPool.commonPool().submit(() -> {
        final Rectangle crop = new Rectangle(Integer.parseInt(textOriginX.getText()),
            Integer.parseInt(textOriginY.getText()), Integer.parseInt(textWidth.getText()),
            Integer.parseInt(textHeight.getText()));
        final ImagePlus imp = cropImp.duplicate();
        imp.setTitle(cropImp.getTitle());
        showCropDialog(bias, gain, variance, imp, crop, action, badPixels);
      }));
      egd.addAndGetButton("Reset", e -> action.accept(reset));
    }
    egd.showDialog();
    if (!egd.wasCanceled()) {
      settings.modelName = egd.getNextString();
      settings.modelDirectory = egd.getNextString();
      settings.originX = (int) egd.getNextNumber();
      settings.originY = (int) egd.getNextNumber();
      settings.width = (int) egd.getNextNumber();
      settings.height = (int) egd.getNextNumber();
      settings.save();
      saveCameraModel(width, height, bias, gain, variance,
          new Rectangle(settings.originX, settings.originY, settings.width, settings.height));
    }
    IJ.showStatus(""); // Remove the status from the ij.io.ImageWriter class
  }

  private static void showHistogram(String name, double[] values, int bins, Statistics stats,
      WindowOrganiser wo) {
    showHistogram(TITLE, name, values, bins, stats, wo);
  }

  /**
   * Show a histogram of the values.
   *
   * @param title the title to prepend to the name
   * @param name the name of the value
   * @param values the values
   * @param bins the number of histogram bins
   * @param stats the summary statistics
   * @param wo the window organiser used for the plot window
   * @return the plot
   */
  static Plot showHistogram(String title, String name, double[] values, int bins, Statistics stats,
      WindowOrganiser wo) {
    final DoubleData data = DoubleData.wrap(values);
    final double minWidth = 0;
    final int removeOutliers = 0;
    final int shape = Plot.CIRCLE;
    final String label = String.format("Mean = %s +/- %s", MathUtils.rounded(stats.getMean()),
        MathUtils.rounded(stats.getStandardDeviation()));

    final HistogramPlot histogramPlot = new HistogramPlotBuilder(title, data, name)
        .setMinBinWidth(minWidth).setRemoveOutliersOption(removeOutliers).setNumberOfBins(bins)
        .setPlotShape(shape).setPlotLabel(label).build();

    histogramPlot.show(wo);
    // Redraw using a log scale. This requires a non-zero y-min
    final Plot plot = histogramPlot.getPlot();
    final double[] limits = plot.getLimits();
    plot.setLimits(limits[0], limits[1], 1, limits[3]);
    plot.setAxisYLog(true);
    plot.updateImage();
    return plot;
  }

  private static <T> void put(BlockingQueue<T> jobs, T job) {
    try {
      jobs.put(job);
    } catch (final InterruptedException ex) {
      Thread.currentThread().interrupt();
      throw new ConcurrentRuntimeException("Unexpected interruption", ex);
    }
  }

  private void computeError() {
    // Assume the simulation stack and measured stack are not null.
    ImageJUtils.log("Comparison to simulation: %s (%dx%dpx)", simulationImp.getInfoProperty(),
        simulationImp.getWidth(), simulationImp.getHeight());
    final WindowOrganiser wo = new WindowOrganiser();
    final ImageStack simulationStack = simulationImp.getImageStack();
    for (int slice = 1; slice <= 3; slice++) {
      computeError(slice, simulationStack, wo);
    }
    wo.tile();
  }

  private void computeError(int slice, ImageStack simulationStack, WindowOrganiser wo) {
    final String label = simulationStack.getSliceLabel(slice);
    final float[] e = (float[]) simulationStack.getPixels(slice);
    final float[] o = (float[]) measuredStack.getPixels(slice);

    // Get the mean error
    final double[] error = new double[e.length];
    for (int i = e.length; i-- > 0;) {
      error[i] = (double) o[i] - e[i];
    }

    final Statistics s = new Statistics();
    s.add(error);
    final StringBuilder result = new StringBuilder(128);
    result.append("Error ").append(label).append(" = ").append(MathUtils.rounded(s.getMean()))
        .append(" +/- ").append(MathUtils.rounded(s.getStandardDeviation()));

    // Do statistical tests
    final double[] x = SimpleArrayUtils.toDouble(e);
    final double[] y = SimpleArrayUtils.toDouble(o);

    final PearsonsCorrelation c = new PearsonsCorrelation();
    result.append(" : R=").append(MathUtils.rounded(c.correlation(x, y)));

    // Plot these
    final String title = TITLE + " " + label + " Simulation vs Measured";
    final Plot plot = new Plot(title, "simulated", "measured");
    plot.addPoints(e, o, Plot.DOT);
    plot.addLabel(0, 0, result.toString());
    ImageJUtils.display(title, plot, wo);

    // Histogram the error
    new HistogramPlotBuilder(TITLE + " " + label, DoubleData.wrap(error), "Error")
        .setPlotLabel(result.toString()).show(wo);

    // Kolmogorov–Smirnov test that the distributions are the same
    SignificanceResult r = KolmogorovSmirnovTest.withDefaults().test(x, y);
    result.append(" : Kolmogorov–Smirnov p=").append(MathUtils.rounded(r.getPValue())).append(' ')
        .append((r.reject(0.001) ? REJECT : ACCEPT));

    if (slice == 3) {
      // Paired T-Test compares two related samples to assess whether their
      // population means differ.
      // T-Test is valid when the difference between the means is normally
      // distributed, e.g. gain
      r = TTest.withDefaults().test(x, y);
      result.append(" : Paired T-Test p=").append(MathUtils.rounded(r.getPValue())).append(' ')
          .append((r.reject(0.001) ? REJECT : ACCEPT));
    } else {
      // Wilcoxon Signed Rank test compares two related samples to assess whether their
      // population mean ranks differ
      r = WilcoxonSignedRankTest.withDefaults().test(x, y);
      result.append(" : Wilcoxon Signed Rank p=").append(MathUtils.rounded(r.getPValue()))
          .append(' ').append((r.reject(0.001) ? REJECT : ACCEPT));
    }

    ImageJUtils.log(result.toString());
  }

  /**
   * Save the camera model. The data is validated. If invalid then the raw data files are saved,
   * otherwise the camera model is saved.
   *
   * @param width the width
   * @param height the height
   * @param bias the bias
   * @param gain the gain
   * @param variance the variance
   * @param crop the crop
   */
  private void saveCameraModel(int width, int height, float[] bias, float[] gain, float[] variance,
      Rectangle crop) {
    PerPixelCameraModel cameraModel;
    try {
      // Crop
      cameraModel = new PerPixelCameraModel(crop.x, crop.y, crop.width, crop.height,
          ImageJImageConverter.getData(bias, width, height, crop, null),
          ImageJImageConverter.getData(gain, width, height, crop, null),
          ImageJImageConverter.getData(variance, width, height, crop, null));
    } catch (final IllegalArgumentException ex) {
      ImageJUtils.log("Invalid per-pixel camera model: %s", ex.getMessage());

      // Save raw data
      final ImageStack stack = new ImageStack(width, height);
      stack.addSlice("Bias", bias);
      stack.addSlice("Gain", gain);
      stack.addSlice("Variance", variance);
      final ImagePlus imp = new ImagePlus(settings.modelName, stack);
      final String filename =
          Paths.get(settings.modelDirectory, settings.modelName + ".error.tif").toString();
      ImageJUtils.log("Saving raw data to: %s", filename);
      new FileSaver(imp).saveAsTiffStack(filename);

      // Check for bad data
      int errors = 0;
      for (int y = 0; y < crop.height && errors < 10; y++) {
        for (int x = 0; x < crop.width && errors < 10; x++) {
          final int i = (crop.y + y) * width + crop.x + x;
          if (!Float.isFinite(bias[i])) {
            ImageJUtils.log("Pixel [%d,%d] Bias %s", crop.x + x, crop.y + y, bias[i]);
            errors++;
          }
          if (!(gain[i] <= Float.MAX_VALUE && gain[i] > 0)) {
            ImageJUtils.log("Pixel [%d,%d] Gain %s", crop.x + x, crop.y + y, gain[i]);
            errors++;
          }
          if (!(variance[i] <= Float.MAX_VALUE && variance[i] >= 0)) {
            ImageJUtils.log("Pixel [%d,%d] Variance %s", crop.x + x, crop.y + y, variance[i]);
            errors++;
          }
        }
      }

      IJ.error(TITLE, "Invalid per-pixel camera model.\n \nCheck the ImageJ log for errors.");
      return;
    }
    if (!CameraModelManager.save(cameraModel,
        new File(settings.modelDirectory, settings.modelName).getPath())) {
      IJ.error(TITLE, "Failed to save model to file");
    }
  }

  /**
   * Show a dialog to allow interactive cropping of the camera model.
   *
   * @param bias the bias
   * @param gain the gain
   * @param variance the variance
   * @param imp the last mean-variance image in the image series
   * @param crop the current camera model crop
   * @param action consumer for the crop (when it is updated)
   * @param badPixels the bad pixels (can be null)
   */
  private static void showCropDialog(float[] bias, float[] gain, float[] variance, ImagePlus imp,
      Rectangle crop, Consumer<Rectangle> action, byte[] badPixels) {
    final int width = imp.getWidth();
    final int height = imp.getHeight();
    ByteProcessor bp;
    Roi roi;
    if (badPixels != null) {
      bp = new ByteProcessor(width, height, badPixels.clone());
      // Find largest bounding region
      roi = RoiHelper.largestBoundingRoi(bp, width / 2, height / 2);

      // Draw on image as red overlay
      final ColorProcessor cp = new ColorProcessor(width, height);
      cp.setColor(Color.RED);
      cp.fillOutside(roi);
      final Overlay o = new Overlay();
      final ImageRoi imageRoi = new ImageRoi(0, 0, cp);
      imageRoi.setZeroTransparent(true);
      imageRoi.setOpacity(0.5);
      o.add(imageRoi);
      imp.setOverlay(o);
    } else {
      bp = new ByteProcessor(width, height);
      roi = new Roi(0, 0, width, height);
    }
    final byte[] pixels = (byte[]) bp.getPixels();

    // Use the bad pixel region if the image is not already cropped
    if (crop.width == width && crop.height == height) {
      imp.setRoi(roi);
      action.accept(roi.getBounds());
    } else {
      imp.setRoi(new Roi(crop));
    }
    final Roi original = (Roi) imp.getRoi().clone();

    // Get the intensity limits of the image
    final ImageProcessor ip = imp.getImageStack().getProcessor(1);
    ip.resetMinAndMax();
    // Min of the crop
    ip.setRoi(imp.getRoi());
    final ImageStatistics stats = ImageStatistics.getStatistics(ip, Measurements.MIN_MAX, null);

    final NonBlockingExtendedGenericDialog gd =
        new NonBlockingExtendedGenericDialog(TITLE + " crop");
    gd.setImage(imp);
    gd.addMessage("Crop the camera model.\n \nThe last mean and variance image is shown:\n"
        + imp.getTitle() + "\n \nKnown bad model pixels were used to outline the\n"
        + "largest valid central region on the image.\n \nCrop using an ROI on the image\n"
        + "or use the foreground intensity.");
    final TextField textForeground =
        gd.addAndGetSlider("Foreground", ip.getMin(), ip.getMax(), stats.min);
    gd.addMessage("");
    final Label label = (Label) gd.getMessage();
    gd.addAndGetButton("Reset",
        e -> SwingUtilities.invokeLater(() -> imp.setRoi((Roi) original.clone())));

    final Consumer<Rectangle> labelAction =
        r -> label.setText(String.format("%,d px", r.width * r.height));
    labelAction.accept(imp.getRoi().getBounds());
    final Consumer<Rectangle> combinedAction = action.andThen(labelAction);
    final RoiListener l = (ImagePlus i, int id) -> {
      if (i != null && i.getID() == imp.getID()) {
        combinedAction.accept(i.getRoi().getBounds());
      }
    };
    Roi.addRoiListener(l);

    imp.show();

    // Dynamically update the ROI using the slider.
    // Maintain a queue of the next value to process and use only
    // 1 thread to compute until the value stops changing.
    final SoftLock lock = new SoftLock();
    final AtomicReference<String> nextValue = new AtomicReference<>();
    textForeground.addTextListener(e -> {
      // queue change
      nextValue.set(textForeground.getText());
      if (lock.acquire()) {
        try {
          ForkJoinPool.commonPool().submit(() -> {
            try {
              for (String next = nextValue.getAndSet(null); next != null;
                  next = nextValue.getAndSet(null)) {
                final float v = Float.parseFloat(next);
                // Find the largest bounding region for pixels above the threshold
                final float[] values = (float[]) ip.getPixels();
                for (int i = 0; i < pixels.length; i++) {
                  pixels[i] = (byte) (values[i] >= v ? 0 : 1);
                }
                imp.setRoi(RoiHelper.largestBoundingRoi(bp, width / 2, height / 2));
              }
            } finally {
              lock.release();
            }
          });
        } catch (final RejectedExecutionException ex) {
          // Something is wrong with the ForkJoinPool,
          // just re-throw after releasing the lock.
          lock.release();
          throw ex;
        }
      }
    });
    gd.hideCancelButton();
    gd.setOKLabel("Done");
    gd.showDialog();

    // Clean-up
    Roi.removeRoiListener(l);
    imp.changes = false;
    imp.close();
  }

  /**
   * Create mask of bad pixels. Returns null if there are no bad pixels.
   *
   * @param bias the bias
   * @param gain the gain
   * @param variance the variance
   * @return the mask (or null)
   */
  private static byte[] findBadPixels(float[] bias, float[] gain, float[] variance) {
    final byte[] pixels = new byte[bias.length];
    boolean found = false;
    for (int i = 0; i < pixels.length; i++) {
      if (!Float.isFinite(bias[i]) || !(gain[i] <= Float.MAX_VALUE && gain[i] > 0)
          || !(variance[i] <= Float.MAX_VALUE && variance[i] >= 0)) {
        pixels[i] = 1;
        found = true;
      }
    }
    return found ? pixels : null;
  }
}
