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

import uk.ac.sussex.gdsc.core.data.IntegerType;
import uk.ac.sussex.gdsc.core.data.SiPrefix;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot.HistogramPlotBuilder;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.logging.Ticker;
import uk.ac.sussex.gdsc.core.math.ArrayMoment;
import uk.ac.sussex.gdsc.core.math.IntegerArrayMoment;
import uk.ac.sussex.gdsc.core.math.RollingArrayMoment;
import uk.ac.sussex.gdsc.core.math.SimpleArrayMoment;
import uk.ac.sussex.gdsc.core.utils.DoubleData;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.StoredData;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.TurboList;
import uk.ac.sussex.gdsc.core.utils.concurrent.CloseableBlockingQueue;
import uk.ac.sussex.gdsc.core.utils.concurrent.ConcurrencyUtils;
import uk.ac.sussex.gdsc.core.utils.rng.Pcg32;
import uk.ac.sussex.gdsc.core.utils.rng.PoissonSamplerUtils;
import uk.ac.sussex.gdsc.core.utils.rng.RandomUtils;
import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;
import uk.ac.sussex.gdsc.smlm.ij.SeriesImageSource;
import uk.ac.sussex.gdsc.smlm.ij.settings.Constants;
import uk.ac.sussex.gdsc.smlm.model.camera.PerPixelCameraModel;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.io.FileSaver;
import ij.io.Opener;
import ij.plugin.PlugIn;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.concurrent.ConcurrentRuntimeException;
import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
import org.apache.commons.math3.stat.inference.TestUtils;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.AhrensDieterExponentialSampler;
import org.apache.commons.rng.sampling.distribution.ContinuousSampler;
import org.apache.commons.rng.sampling.distribution.DiscreteSampler;
import org.apache.commons.rng.sampling.distribution.GaussianSampler;
import org.apache.commons.rng.sampling.distribution.NormalizedGaussianSampler;
import org.apache.commons.rng.sampling.distribution.SharedStateContinuousSampler;
import org.apache.commons.rng.simple.RandomSource;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicReference;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Analyse the per pixel offset, variance and gain from a sCMOS camera.
 *
 * <p>See Huang et al (2013) Video-rate nanoscopy using sCMOS camera–specific single-molecule
 * localization algorithms. Nature Methods 10, 653-658 (Supplementary Information).
 */
public class CmosAnalysis implements PlugIn {
  private static final String TITLE = "sCMOS Analysis";

  private int numberOfThreads;
  // The simulated offset, variance and gain
  private ImagePlus simulationImp;
  // The measured offset, variance and gain
  private ImageStack measuredStack;
  // The sub-directories containing the sCMOS images
  private TurboList<SubDir> subDirs;

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    String directory;
    String modelDirectory;
    String modelName;
    boolean rollingAlgorithm;
    boolean reuseProcessedData;
    double offset;
    double variance;
    double gain;
    double gainStdDev;
    int size;
    int frames;
    int imagejNThreads;
    int lastNumberOfThreads;

    Settings() {
      // Set defaults
      directory = Prefs.get(Constants.sCMOSAnalysisDirectory, "");
      reuseProcessedData = true;

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
      frames = 512;

      imagejNThreads = Prefs.getThreads();
      lastNumberOfThreads = imagejNThreads;
    }

    Settings(Settings source) {
      directory = source.directory;
      modelDirectory = source.modelDirectory;
      modelName = source.modelName;
      rollingAlgorithm = source.rollingAlgorithm;
      reuseProcessedData = source.reuseProcessedData;
      offset = source.offset;
      variance = source.variance;
      gain = source.gain;
      gainStdDev = source.gainStdDev;
      size = source.size;
      frames = source.frames;
      imagejNThreads = source.imagejNThreads;
      lastNumberOfThreads = source.lastNumberOfThreads;
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
     * Save the settings. This can be called only once as it saves via a reference.
     */
    void save() {
      lastSettings.set(this);
    }
  }

  private static class SimulationWorker implements Runnable {
    static final short MIN_SHORT = 0;
    // Require cast since this is out-of-range for a signed short
    static final short MAX_SHORT = (short) 65335;

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
      int[] poisson = null;
      if (photons != 0) {
        // For speed we can precompute a set of random numbers to reuse
        final DiscreteSampler pd = PoissonSamplerUtils.createPoissonSampler(rg, photons);
        poisson = new int[pixelVariance.length];
        for (int i = poisson.length; i-- > 0;) {
          poisson[i] = pd.sample();
        }
      }

      // Save image in blocks
      ImageStack stack = new ImageStack(size, size);
      int start = from;
      for (int i = from; i < to; i++) {
        // Create image
        final short[] pixels = new short[pixelOffset.length];
        if (poisson == null) {
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
                pixelOffset[j] + gauss.sample() * pixelSd[j] + (poisson[j] * pixelGain[j]);
            pixels[j] = clip16bit(p);
          }

          // Rotate Poisson numbers.
          // Shuffling what we have is faster than generating new values
          // and we should have enough.
          RandomUtils.shuffle(poisson, rg);
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
      final int i = (int) Math.round(value);
      if (i < 0) {
        return MIN_SHORT;
      }
      if (i > 65335) {
        return MAX_SHORT;
      }
      return (short) i;
    }

    private void save(ImageStack stack, int start) {
      final ImagePlus imp = new ImagePlus("", stack);
      final String path = new File(out, String.format("image%06d.tif", start)).getPath();
      final FileSaver fs = new FileSaver(imp);
      fs.saveAsTiffStack(path);
    }
  }

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
          if (pixels == null) {
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
    settings.directory = dir;
    settings.save();
    Prefs.set(Constants.sCMOSAnalysisDirectory, dir);

    final boolean simulate = "simulate".equals(arg);
    if (simulate || extraOptions) {
      if (!showSimulateDialog()) {
        return;
      }
      try {
        simulate();
      } catch (IOException ex) {
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

    final UniformRandomProvider rg = RandomSource.create(RandomSource.XOR_SHIFT_1024_S);

    final DiscreteSampler pd = PoissonSamplerUtils.createPoissonSampler(rg, settings.offset);
    final ContinuousSampler ed = new AhrensDieterExponentialSampler(rg, settings.variance);
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
    final TurboList<Future<?>> futures = new TurboList<>(numberOfThreads);

    // Simulate the zero exposure input.
    // Simulate 20 - 200 photon images.
    final int[] photons = new int[] {0, 20, 50, 100, 200};

    final int blockSize = 10; // For saving stacks
    int numberPerThread = (int) Math.ceil((double) settings.frames / numberOfThreads);
    // Convert to fit the block size
    numberPerThread = (int) Math.ceil((double) numberPerThread / blockSize) * blockSize;
    final Pcg32 rng = Pcg32.xshrs(start);
    ticker = ImageJUtils.createTicker((long) photons.length * settings.frames, threadCount);
    for (final int p : photons) {
      ImageJUtils.showStatus(() -> "Simulating " + TextUtils.pleural(p, "photon"));

      // Create the directory
      Path out = Paths.get(settings.directory, String.format("photon%03d", p));
      Files.createDirectories(out);

      for (int from = 0; from < settings.frames;) {
        final int to = Math.min(from + numberPerThread, settings.frames);
        futures.add(executor.submit(new SimulationWorker(ticker, rng.split(),
            out.toString(), simulationStack, from, to, blockSize, p)));
        from = to;
      }

      ConcurrencyUtils.waitForCompletionUnchecked(futures);
      futures.clear();
    }

    ImageJUtils.finished();

    executor.shutdown();

    ImageJUtils
        .log("Simulation time = " + TextUtils.millisToString(System.currentTimeMillis() - start));
  }

  private boolean showSimulateDialog() {
    final GenericDialog gd = new GenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);

    gd.addMessage("Simulate per-pixel offset, variance and gain of sCMOS images.");

    gd.addNumericField("nThreads", getLastNumberOfThreads(), 0);
    gd.addNumericField("Offset (Poisson)", settings.offset, 3);
    gd.addNumericField("Variance (Exponential)", settings.variance, 3);
    gd.addNumericField("Gain (Gaussian)", settings.gain, 3);
    gd.addNumericField("Gain_SD", settings.gainStdDev, 3);
    gd.addNumericField("Size", settings.size, 0);
    gd.addNumericField("Frames", settings.frames, 0);
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
    subDirs = new TurboList<>();
    final Pattern p = Pattern.compile("([0-9]+)$");
    for (final File path : dirs) {
      final String name = path.getName();
      final Matcher m = p.matcher(name);
      if (m.find()) {
        final int t = Integer.parseInt(m.group(1));
        subDirs.add(new SubDir(t, path, name));
      }
    }

    if (subDirs.size() < 2) {
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
    gd.addHelp(About.HELP_URL);

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
    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    setThreads((int) gd.getNextNumber());
    settings.rollingAlgorithm = gd.getNextBoolean();

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
    final TurboList<Future<?>> futures = new TurboList<>(nThreads);
    final TurboList<ImageWorker> workers = new TurboList<>(nThreads);

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

    for (int n = 0; n < nSubDirs; n++) {
      ImageJUtils.showSlowProgress(0, nSubDirs);

      final SubDir sd = subDirs.getf(n);
      ImageJUtils.showStatus(() -> "Analysing " + sd.name);
      final StopWatch sw = StopWatch.createStarted();

      // Option to reuse data
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
        }
      }

      if (!found) {
        // Open the series
        final SeriesImageSource source = new SeriesImageSource(sd.name, sd.path.getPath());
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

        final CloseableBlockingQueue<Object> jobs = new CloseableBlockingQueue<>(nThreads * 2);
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
          jobs.close(true);
          workers.stream().forEach(worker -> worker.finished = true);
          // Cancel by interruption. We set the finished flag so the ImageWorker should
          // ignore the interrupt.
          futures.stream().forEach(future -> future.cancel(true));
          break;
        }

        // Finish all the worker threads by passing in a null job
        jobs.close(false);

        // Wait for all to finish
        ConcurrencyUtils.waitForCompletionUnchecked(futures);

        // Create the final aggregate statistics
        for (final ImageWorker w : workers) {
          moment.add(w.moment);
        }
        data[2 * n] = moment.getFirstMoment();
        data[2 * n + 1] = moment.getVariance();

        // Get the processing speed.
        sw.stop();
        // ticker holds the number of number of frames processed
        final double bits =
            (double) bitDepth * ticker.getCurrent() * source.getWidth() * source.getHeight();
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
        IJ.save(new ImagePlus("PerPixel", stack), file.getPath());
      }

      final Statistics s = Statistics.create(data[2 * n]);

      if (pixelOffset != null) {
        // Compute mean ADU
        final Statistics signal = new Statistics();
        final double[] mean = data[2 * n];
        for (int i = 0; i < pixelOffset.length; i++) {
          signal.add(mean[i] - pixelOffset[i]);
        }
        ImageJUtils.log("%s Mean = %s +/- %s. Signal = %s +/- %s ADU", sd.name,
            MathUtils.rounded(s.getMean()), MathUtils.rounded(s.getStandardDeviation()),
            MathUtils.rounded(signal.getMean()), MathUtils.rounded(signal.getStandardDeviation()));
      } else {
        // Set the offset assuming the first sub-directory is the bias image
        pixelOffset = data[0];
        pixelVariance = data[1];
        statsOffset = s;
        statsVariance = Statistics.create(pixelVariance);
        ImageJUtils.log("%s Offset = %s +/- %s. Variance = %s +/- %s", sd.name,
            MathUtils.rounded(s.getMean()), MathUtils.rounded(s.getStandardDeviation()),
            MathUtils.rounded(statsVariance.getMean()),
            MathUtils.rounded(statsVariance.getStandardDeviation()));
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
    final int bins = 2 * HistogramPlot.getBinsSturgesRule(pixelGain.length);
    final WindowOrganiser wo = new WindowOrganiser();
    showHistogram("Offset (ADU)", pixelOffset, bins, statsOffset, wo);
    showHistogram("Variance (ADU^2)", pixelVariance, bins, statsVariance, wo);
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

    final ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE);
    egd.addMessage("Save the sCMOS camera model?");
    if (settings.modelDirectory == null) {
      settings.modelDirectory = settings.directory;
      settings.modelName = "sCMOS Camera";
    }
    egd.addStringField("Model_name", settings.modelName, 30);
    egd.addDirectoryField("Model_directory", settings.modelDirectory);
    egd.showDialog();
    if (!egd.wasCanceled()) {
      settings.modelName = egd.getNextString();
      settings.modelDirectory = egd.getNextString();
      final PerPixelCameraModel cameraModel =
          new PerPixelCameraModel(width, height, bias, gain, variance);
      if (!CameraModelManager.save(cameraModel,
          new File(settings.directory, settings.modelName).getPath())) {
        IJ.error(TITLE, "Failed to save model to file");
      }
    }
    IJ.showStatus(""); // Remove the status from the ij.io.ImageWriter class

    ImageJUtils
        .log("Analysis time = " + TextUtils.millisToString(System.currentTimeMillis() - start));
  }

  private static void showHistogram(String name, double[] values, int bins, Statistics stats,
      WindowOrganiser wo) {
    final DoubleData data = new StoredData(values, false);
    final double minWidth = 0;
    final int removeOutliers = 0;
    final int shape = Plot.CIRCLE;
    final String label = String.format("Mean = %s +/- %s", MathUtils.rounded(stats.getMean()),
        MathUtils.rounded(stats.getStandardDeviation()));

    final HistogramPlot histogramPlot = new HistogramPlotBuilder(TITLE, data, name)
        .setMinBinWidth(minWidth).setRemoveOutliersOption(removeOutliers).setNumberOfBins(bins)
        .setPlotShape(shape).setPlotLabel(label).build();

    histogramPlot.show(wo);
    // Redraw using a log scale. This requires a non-zero y-min
    final Plot plot = histogramPlot.getPlot();
    final double[] limits = plot.getLimits();
    plot.setLimits(limits[0], limits[1], 1, limits[3]);
    plot.setAxisYLog(true);
    plot.updateImage();
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
    ImageJUtils.log("Comparison to simulation: %s", simulationImp.getInfoProperty());
    final ImageStack simulationStack = simulationImp.getImageStack();
    for (int slice = 1; slice <= 3; slice++) {
      computeError(slice, simulationStack);
    }
  }

  private void computeError(int slice, ImageStack simulationStack) {
    final String label = simulationStack.getSliceLabel(slice);
    final float[] e = (float[]) simulationStack.getPixels(slice);
    final float[] o = (float[]) measuredStack.getPixels(slice);

    // Get the mean error
    final Statistics s = new Statistics();
    for (int i = e.length; i-- > 0;) {
      s.add(o[i] - e[i]);
    }

    final StringBuilder result = new StringBuilder("Error ").append(label);
    result.append(" = ").append(MathUtils.rounded(s.getMean()));
    result.append(" +/- ").append(MathUtils.rounded(s.getStandardDeviation()));

    // Do statistical tests
    final double[] x = SimpleArrayUtils.toDouble(e);
    final double[] y = SimpleArrayUtils.toDouble(o);

    final PearsonsCorrelation c = new PearsonsCorrelation();
    result.append(" : R=").append(MathUtils.rounded(c.correlation(x, y)));

    // Mann-Whitney U is valid for any distribution, e.g. variance
    final MannWhitneyUTest test = new MannWhitneyUTest();
    double pvalue = test.mannWhitneyUTest(x, y);
    result.append(" : Mann-Whitney U p=").append(MathUtils.rounded(pvalue)).append(' ')
        .append(((pvalue < 0.05) ? "reject" : "accept"));

    if (slice != 2) {
      // T-Test is valid for approximately Normal distributions, e.g. offset and gain
      pvalue = TestUtils.tTest(x, y);
      result.append(" : T-Test p=").append(MathUtils.rounded(pvalue)).append(' ')
          .append(((pvalue < 0.05) ? "reject" : "accept"));
      pvalue = TestUtils.pairedTTest(x, y);
      result.append(" : Paired T-Test p=").append(MathUtils.rounded(pvalue)).append(' ')
          .append(((pvalue < 0.05) ? "reject" : "accept"));
    }

    ImageJUtils.log(result.toString());
  }
}
