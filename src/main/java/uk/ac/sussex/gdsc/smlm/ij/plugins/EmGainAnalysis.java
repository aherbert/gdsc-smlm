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

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.gui.Roi;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Point;
import java.awt.Rectangle;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicReference;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateFunctionMappingAdapter;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.NormalizedGaussianSampler;
import org.apache.commons.rng.sampling.distribution.PoissonSampler;
import uk.ac.sussex.gdsc.core.data.ComputationException;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.logging.Ticker;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.core.utils.rng.MarsagliaTsangGammaSampler;
import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;
import uk.ac.sussex.gdsc.core.utils.rng.UniformRandomProviders;
import uk.ac.sussex.gdsc.smlm.function.LikelihoodFunction;
import uk.ac.sussex.gdsc.smlm.function.PoissonFunction;
import uk.ac.sussex.gdsc.smlm.function.PoissonGammaFunction;
import uk.ac.sussex.gdsc.smlm.function.PoissonGammaGaussianFunction;
import uk.ac.sussex.gdsc.smlm.function.PoissonGaussianFunction2;
import uk.ac.sussex.gdsc.smlm.math3.optim.nonlinear.scalar.noderiv.CustomPowellOptimizer;
import uk.ac.sussex.gdsc.smlm.utils.Convolution;

/**
 * Analysis a white light image from an EM-CCD camera, construct a histogram of pixel intensity and
 * fit the histogram to obtain the bias, EM-gain, read noise and photons per pixel.
 *
 * <p>See Ulbrich &amp; Isacoff (2007) Subunit counting in membrane-bound proteins. Nature Methods
 * 4, 319-321 (Supplementary Information).
 */
public class EmGainAnalysis implements PlugInFilter {
  private static final String TITLE = "EM-Gain Analysis";
  private static final int FLAGS = DOES_8G | DOES_16 | NO_CHANGES | NO_UNDO;
  private static final int MINIMUM_PIXELS = 1000000;

  private boolean extraOptions;

  private ImagePlus imp;

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    static final String[] APPROXIMATION_TYPE =
        {"PoissonGammaGaussian", "PoissonGamma", "PoissonGaussian", "Poisson"};

    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    double bias;
    double gain;
    double noise;
    boolean settingSimulate;
    boolean showApproximation;
    boolean relativeDelta;
    int approximationType;
    double settingPhotons;
    double settingBias;
    double settingGain;
    double settingNoise;
    double head;
    double tail;
    double settingOffset;
    int simulationSize;
    boolean usePdf;

    Settings() {
      // Set defaults
      bias = 500;
      gain = 40;
      noise = 3;
      settingPhotons = 1;
      settingBias = 500;
      settingGain = 40;
      settingNoise = 3;
      head = 0.01;
      tail = 0.025;
      simulationSize = 20000;
    }

    Settings(Settings source) {
      bias = source.bias;
      gain = source.gain;
      noise = source.noise;
      settingSimulate = source.settingSimulate;
      showApproximation = source.showApproximation;
      relativeDelta = source.relativeDelta;
      approximationType = source.approximationType;
      settingPhotons = source.settingPhotons;
      settingBias = source.settingBias;
      settingGain = source.settingGain;
      settingNoise = source.settingNoise;
      head = source.head;
      tail = source.tail;
      settingOffset = source.settingOffset;
      simulationSize = source.simulationSize;
      usePdf = source.usePdf;
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

  /**
   * Store the probability density function (PDF).
   */
  private static class Pdf {
    /** The observed value x. */
    final double[] x;
    /** The probability. */
    final double[] probability;

    /**
     * Instantiates a new pdf.
     *
     * @param x the observed value x
     * @param probability the probability
     */
    Pdf(double[] x, double[] probability) {
      this.x = x;
      this.probability = probability;
    }
  }

  @Override
  public int setup(String arg, ImagePlus imp) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    extraOptions = ImageJUtils.isExtraOptions();
    settings = Settings.load();

    if ("pmf".equals(arg)) {
      plotPmf();
      return DONE;
    }

    if (imp == null && !extraOptions) {
      IJ.noImage();
      return DONE;
    }
    this.imp = imp;
    return showDialog();
  }

  @Override
  public void run(ImageProcessor ip) {
    // Calculate the histogram
    final int[] histogram =
        (settings.settingSimulate) ? simulateHistogram(settings.usePdf) : buildHistogram(imp);

    // We need > 10^7 pixels from flat white-light images under constant exposure ...
    final int size = getSize(histogram);
    if (imp != null) {
      final Roi roi = imp.getRoi();
      Rectangle bounds;
      if (roi == null) {
        bounds = new Rectangle(0, 0, imp.getWidth(), imp.getHeight());
      } else {
        bounds = roi.getBounds();
      }
      ImageJUtils.log("Analysing %s [x=%d,y=%d,width=%d,height=%d]", imp.getTitle(), bounds.x,
          bounds.y, bounds.width, bounds.height);
    }
    ImageJUtils.log("Histogram contains %d pixels", size);
    if (size < MINIMUM_PIXELS) {
      ImageJUtils.log("WARNING : Recommend at least %d pixels (%sx more)", MINIMUM_PIXELS,
          MathUtils.rounded((double) MINIMUM_PIXELS / size));
    }

    fit(histogram);
  }

  /**
   * Simulate the histogram for fitting. Sample from the fitted PDF or sample from a
   * Poisson-Gamma-Gaussian.
   *
   * @param usePdf the use pdf
   * @return The histogram
   */
  private int[] simulateHistogram(boolean usePdf) {
    IJ.showStatus("Simulating histogram ...");
    final int[] histogram = usePdf ? simulateFromPdf() : simulateFromPoissonGammaGaussian();
    IJ.showStatus("");
    return histogram;
  }

  /**
   * Random sample from the cumulative probability distribution function that is used during
   * fitting.
   *
   * @return The histogram
   */
  private int[] simulateFromPdf() {
    final double step =
        getStepSize(settings.settingPhotons, settings.settingGain, settings.settingNoise);

    final Pdf pdf =
        pdf(0, step, settings.settingPhotons, settings.settingGain, settings.settingNoise);

    // Debug this
    final double[] g = pdf.probability;
    final double[] x = pdf.x;
    final Plot plot = new Plot(TITLE + " PDF", "ADU", "P");
    plot.addPoints(x, Arrays.copyOf(g, g.length), Plot.LINE);
    ImageJUtils.display(TITLE + " PDF", plot);

    // Get cumulative probability
    double sum = 0;
    for (int i = 0; i < g.length; i++) {
      final double p = g[i];
      g[i] += sum;
      sum += p;
    }

    for (int i = 0; i < g.length; i++) {
      g[i] /= sum;
    }
    g[g.length - 1] = 1; // Ensure value of 1 at the end

    // Randomly sample
    final UniformRandomProvider rng = UniformRandomProviders.create();
    final int bias = (int) settings.settingBias;
    final int[] bins = new int[x.length];
    for (int i = 0; i < x.length; i++) {
      bins[i] = bias + (int) x[i];
    }
    final int[] h = new int[bins[bins.length - 1] + 1];
    final int steps = settings.simulationSize;
    final Ticker ticker = ImageJUtils.createTicker(steps, 1);
    for (int n = 0; n < steps; n++) {
      int index = binarySearch(g, rng.nextDouble());
      if (index < 0) {
        index = -(index + 1);
      }
      h[bins[index]]++;
      ticker.tick();
    }
    return h;
  }

  private static int binarySearch(double[] a, double key) {
    int low = 0;
    int high = a.length - 1;

    while (low <= high) {
      final int mid = (low + high) >>> 1;
      final double midVal = a[mid];

      if (midVal < key) {
        low = mid + 1; // Neither val is NaN, thisVal is smaller
      } else if (midVal > key) {
        high = mid - 1; // Neither val is NaN, thisVal is larger
      } else {
        return mid; // Key found
      }
    }
    return -(low + 1); // key not found.
  }

  /**
   * Randomly generate a histogram from poisson-gamma-gaussian samples.
   *
   * @return The histogram
   */
  private int[] simulateFromPoissonGammaGaussian() {
    // Randomly sample
    final UniformRandomProvider rng = UniformRandomProviders.create();

    final PoissonSampler poisson = new PoissonSampler(rng, settings.settingPhotons);
    final MarsagliaTsangGammaSampler gamma =
        new MarsagliaTsangGammaSampler(rng, settings.settingPhotons, settings.settingGain);
    final NormalizedGaussianSampler gauss = SamplerUtils.createNormalizedGaussianSampler(rng);

    final int steps = settings.simulationSize;
    final int[] samples = new int[steps];
    for (int n = 0; n < steps; n++) {
      if (n % 64 == 0) {
        IJ.showProgress(n, steps);
      }

      // Poisson
      double sample = poisson.sample();

      // Gamma
      if (sample > 0) {
        gamma.setAlpha(sample);
        sample = gamma.sample();
      }

      // Gaussian
      sample += settings.settingNoise * gauss.sample();

      // Convert the sample to a count
      samples[n] = (int) Math.round(sample + settings.settingBias);
    }

    final int max = MathUtils.max(samples);
    final int[] histogram = new int[max + 1];
    for (final int s : samples) {
      histogram[s]++;
    }
    return histogram;
  }

  /**
   * Build a histogram using pixels within the image ROI.
   *
   * @param imp The image
   * @return The image histogram
   */
  private static int[] buildHistogram(ImagePlus imp) {
    final ImageStack stack = imp.getImageStack();
    final Roi roi = imp.getRoi();
    final int[] data = getHistogram(stack.getProcessor(1), roi);
    for (int n = 2; n <= stack.getSize(); n++) {
      final int[] tmp = getHistogram(stack.getProcessor(n), roi);
      for (int i = 0; i < tmp.length; i++) {
        data[i] += tmp[i];
      }
    }

    // Avoid super-saturated pixels by using 98% of histogram
    long sum = 0;
    for (int i = 0; i < data.length; i++) {
      sum += data[i];
    }

    final long sum2 = (long) (sum * 0.99);
    sum = 0;
    for (int i = 0; i < data.length; i++) {
      sum += data[i];
      if (sum > sum2) {
        for (; i < data.length; i++) {
          data[i] = 0;
        }
        break;
      }
    }

    return data;
  }

  private static int[] getHistogram(ImageProcessor ip, Roi roi) {
    ip.setRoi(roi);
    return ip.getHistogram();
  }

  /**
   * Fit the EM-gain distribution (Gaussian * Gamma).
   *
   * @param histogram The distribution
   */
  private void fit(int[] histogram) {
    final int[] limits = limits(histogram);
    final double[] x = getX(limits);
    final double[] y = getY(histogram, limits);

    Plot plot = new Plot(TITLE, "ADU", "Frequency");
    double yMax = MathUtils.max(y);
    plot.setLimits(limits[0], limits[1], 0, yMax);
    plot.setColor(Color.black);
    plot.addPoints(x, y, Plot.DOT);
    ImageJUtils.display(TITLE, plot);

    // Estimate remaining parameters.
    // Assuming a gamma_distribution(shape,scale) then mean = shape * scale
    // scale = gain
    // shape = Photons = mean / gain
    double mean = getMean(histogram) - settings.bias;
    // Note: if the bias is too high then the mean will be negative. Just move the bias.
    while (mean < 0) {
      settings.bias -= 1;
      mean += 1;
    }
    double photons = mean / settings.gain;

    if (settings.settingSimulate) {
      ImageJUtils.log("Simulated bias=%d, gain=%s, noise=%s, photons=%s",
          (int) settings.settingBias, MathUtils.rounded(settings.settingGain),
          MathUtils.rounded(settings.settingNoise), MathUtils.rounded(settings.settingPhotons));
    }

    ImageJUtils.log("Estimate bias=%d, gain=%s, noise=%s, photons=%s", (int) settings.bias,
        MathUtils.rounded(settings.gain), MathUtils.rounded(settings.noise),
        MathUtils.rounded(photons));

    final int max = (int) x[x.length - 1];
    double[] pg = pdf(max, photons, settings.gain, settings.noise, (int) settings.bias);

    plot.setColor(Color.blue);
    plot.addPoints(x, pg, Plot.LINE);
    ImageJUtils.display(TITLE, plot);

    // Perform a fit
    final CustomPowellOptimizer o = new CustomPowellOptimizer(1e-6, 1e-16, 1e-6, 1e-16);
    final double[] startPoint =
        new double[] {photons, settings.gain, settings.noise, settings.bias};
    int maxEval = 3000;

    final String[] paramNames = {"Photons", "Gain", "Noise", "Bias"};
    // Set bounds
    final double[] lower = new double[] {0, 0.5 * settings.gain, 0, settings.bias - settings.noise};
    final double[] upper = new double[] {2 * photons, 2 * settings.gain, settings.gain,
        settings.bias + settings.noise};

    // Restart until converged.
    // TODO - Maybe fix this with a better optimiser. This needs to be tested on real data.
    PointValuePair solution = null;
    for (int iter = 0; iter < 3; iter++) {
      IJ.showStatus("Fitting histogram ... Iteration " + iter);

      try {
        // Basic Powell optimiser
        final MultivariateFunction fun = getFunction(limits, y, max, maxEval);
        final PointValuePair optimum =
            o.optimize(new MaxEval(maxEval), new ObjectiveFunction(fun), GoalType.MINIMIZE,
                new InitialGuess((solution == null) ? startPoint : solution.getPointRef()));
        if (solution == null || optimum.getValue() < solution.getValue()) {
          final double[] point = optimum.getPointRef();
          // Check the bounds
          for (int i = 0; i < point.length; i++) {
            if (point[i] < lower[i] || point[i] > upper[i]) {
              throw new ComputationException(
                  String.format("Fit out of of estimated range: %s %f", paramNames[i], point[i]));
            }
          }
          solution = optimum;
        }
      } catch (final Exception ex) {
        IJ.log("Powell error: " + ex.getMessage());
        if (ex instanceof TooManyEvaluationsException) {
          maxEval = (int) (maxEval * 1.5);
        }
      }
      try {
        // Bounded Powell optimiser
        final MultivariateFunction fun = getFunction(limits, y, max, maxEval);
        final MultivariateFunctionMappingAdapter adapter =
            new MultivariateFunctionMappingAdapter(fun, lower, upper);
        PointValuePair optimum = o.optimize(new MaxEval(maxEval), new ObjectiveFunction(adapter),
            GoalType.MINIMIZE, new InitialGuess(adapter
                .boundedToUnbounded((solution == null) ? startPoint : solution.getPointRef())));
        final double[] point = adapter.unboundedToBounded(optimum.getPointRef());
        optimum = new PointValuePair(point, optimum.getValue());

        if (solution == null || optimum.getValue() < solution.getValue()) {
          solution = optimum;
        }
      } catch (final Exception ex) {
        IJ.log("Bounded Powell error: " + ex.getMessage());
        if (ex instanceof TooManyEvaluationsException) {
          maxEval = (int) (maxEval * 1.5);
        }
      }
    }

    ImageJUtils.finished();

    if (solution == null) {
      ImageJUtils.log("Failed to fit the distribution");
      return;
    }

    final double[] point = solution.getPointRef();
    photons = point[0];
    settings.gain = point[1];
    settings.noise = point[2];
    settings.bias = (int) Math.round(point[3]);
    final String label = String.format("Fitted bias=%d, gain=%s, noise=%s, photons=%s",
        (int) settings.bias, MathUtils.rounded(settings.gain), MathUtils.rounded(settings.noise),
        MathUtils.rounded(photons));
    ImageJUtils.log(label);

    if (settings.settingSimulate) {
      ImageJUtils.log("Relative Error bias=%s, gain=%s, noise=%s, photons=%s",
          MathUtils.rounded(relativeError(settings.bias, settings.settingBias)),
          MathUtils.rounded(relativeError(settings.gain, settings.settingGain)),
          MathUtils.rounded(relativeError(settings.noise, settings.settingNoise)),
          MathUtils.rounded(relativeError(photons, settings.settingPhotons)));
    }

    // Show the PoissonGammaGaussian approximation
    double[] approxValues = null;
    if (settings.showApproximation) {
      approxValues = new double[x.length];
      final PoissonGammaGaussianFunction fun =
          new PoissonGammaGaussianFunction(1.0 / settings.gain, settings.noise);
      final double expected = photons * settings.gain;
      for (int i = 0; i < approxValues.length; i++) {
        approxValues[i] = fun.likelihood(x[i] - settings.bias, expected);
      }
      yMax = MathUtils.maxDefault(yMax, approxValues);
    }

    // Replot
    pg = pdf(max, photons, settings.gain, settings.noise, (int) settings.bias);
    plot = new Plot(TITLE, "ADU", "Frequency");
    plot.setLimits(limits[0], limits[1], 0, yMax * 1.05);
    plot.setColor(Color.black);
    plot.addPoints(x, y, Plot.DOT);
    plot.setColor(Color.red);
    plot.addPoints(x, pg, Plot.LINE);

    plot.addLabel(0, 0, label);

    if (settings.showApproximation) {
      plot.setColor(Color.blue);
      plot.addPoints(x, approxValues, Plot.LINE);
    }

    ImageJUtils.display(TITLE, plot);
  }

  private static double relativeError(double v1, double v2) {
    final double delta = v1 - v2; // Math.abs(v1 - v2)
    return delta / v2;
  }

  private static MultivariateFunction getFunction(final int[] limits, final double[] y,
      final int max, final int maxEval) {
    return new MultivariateFunction() {
      int eval;

      @Override
      public double value(double[] point) {
        IJ.showProgress(++eval, maxEval);
        if (ImageJUtils.isInterrupted()) {
          throw new TooManyEvaluationsException(maxEval);
        }
        // Compute the sum of squares between the two functions
        final double photons = point[0];
        final double gain = point[1];
        final double noise = point[2];
        final int bias = (int) Math.round(point[3]);
        final double[] g = pdf(max, photons, gain, noise, bias);
        double ss = 0;
        for (int c = limits[0]; c <= limits[1]; c++) {
          final double d = g[c] - y[c];
          ss += d * d;
        }
        return ss;
      }
    };
  }

  /**
   * Calculate the probability density function for EM-gain. The maximum count to evaluate is
   * calculated dynamically so that the cumulative probability does not change.
   *
   * <p>See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
   *
   * @param step the step between counts to evaluate
   * @param photons The average number of photons per pixel input to the EM-camera
   * @param gain The multiplication factor (gain)
   * @return The PDF
   */
  private static double[] pdfEmGain(final double step, final double photons, final double gain) {
    final StoredDataStatistics stats = new StoredDataStatistics(100);
    stats.add(FastMath.exp(-photons));
    for (int c = 1;; c++) {
      final double g = probabilityEmGain(c * step, photons, gain);
      stats.add(g);
      final double delta = g / stats.getSum();
      if (delta < 1e-5) {
        break;
      }
    }
    return stats.getValues();
  }

  /**
   * Calculate the probability density function for EM-gain.
   *
   * <p>See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
   *
   * @param max The maximum count to evaluate
   * @param step the step between counts to evaluate
   * @param photons The average number of photons per pixel input to the EM-camera
   * @param gain The multiplication factor (gain)
   * @return The PDF
   */
  private static double[] pdfEmGain(final int max, final double step, final double photons,
      final double gain) {
    if (max == 0) {
      return pdfEmGain(step, photons, gain);
    }
    final double[] g = new double[max + 1];
    g[0] = FastMath.exp(-photons);
    for (int c = 1;; c++) {
      final double count = c * step;
      g[c] = probabilityEmGain(count, photons, gain);
      if (g[c] == 0 || count >= max) {
        break;
      }
    }
    return g;
  }

  /**
   * Calculate the probability density function for EM-gain.
   *
   * <p>See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
   *
   * @param count The count to evaluate
   * @param photons The average number of photons per pixel input to the EM-camera
   * @param gain The multiplication factor (gain)
   * @return The PDF
   */
  private static double probabilityEmGain(double count, double photons, double gain) {
    return PoissonGammaFunction.poissonGamma(count, photons, gain);
  }

  /**
   * Calculate the probability density function for EM-gain, convolve with a Gaussian and then add a
   * constant offset.
   *
   * <p>See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI.
   *
   * @param max The maximum count to evaluate
   * @param step the step between counts to evaluate
   * @param photons The average number of photons per pixel input to the EM-camera
   * @param gain The multiplication factor (gain)
   * @param sd The read noise (Gaussian standard deviation)
   * @return The PDF
   */
  private Pdf pdf(final int max, final double step, final double photons, final double gain,
      final double sd) {
    final double[] g = pdfEmGain(max, step, photons, gain);
    double[] gg;

    int zero = 0;

    if (sd > 0) {
      // Convolve with Gaussian kernel up to 4 times the standard deviation
      final int radius = (int) Math.ceil(Math.abs(sd) * 4 / step) + 1;
      final double[] kernel = new double[2 * radius + 1];
      final double norm = -0.5 / (sd * sd);
      for (int i = 0, j = radius, jj = radius; j < kernel.length; i++, j++, jj--) {
        kernel[j] = kernel[jj] = FastMath.exp(norm * MathUtils.pow2(i * step));
      }
      // Normalise
      double sum = 0;
      for (int j = 0; j < kernel.length; j++) {
        sum += kernel[j];
      }
      for (int j = 0; j < kernel.length; j++) {
        kernel[j] /= sum;
      }

      if (extraOptions) {
        // Debug
        String title = "Poisson-Gamma";
        Plot plot = new Plot(title, "x", "y");
        plot.addPoints(SimpleArrayUtils.newArray(g.length, 0, step), g, Plot.LINE);
        ImageJUtils.display(title, plot);

        title = "Gaussian";
        plot = new Plot(title, "x", "y");
        plot.addPoints(SimpleArrayUtils.newArray(kernel.length, radius * -step, step), kernel,
            Plot.LINE);
        ImageJUtils.display(title, plot);
      }

      gg = Convolution.convolveFast(g, kernel);
      // The convolution will have created a larger array so we must adjust the offset for this
      zero = radius;
    } else {
      gg = g;
    }

    final double[] x = new double[gg.length];
    for (int i = 0, j = -zero; i < x.length; i++, j++) {
      x[i] = j * step;
    }

    return new Pdf(x, gg);
  }

  /**
   * Calculate the probability density function for EM-gain, convolve with a Gaussian and then add a
   * constant offset.
   *
   * <p>See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI.
   *
   * @param max The maximum count to evaluate
   * @param photons The average number of photons per pixel input to the EM-camera
   * @param gain The multiplication factor (gain)
   * @param sd The read noise (Gaussian standard deviation)
   * @param c0 The constant offset (bias)
   * @return The PDF
   */
  private static double[] pdf(final int max, final double photons, final double gain,
      final double sd, int c0) {
    final double[] g = pdfEmGain(max, photons, gain);
    double[] gg;

    if (sd > 0) {
      // Convolve with Gaussian kernel up to 4 times the standard deviation
      final int radius = (int) Math.ceil(Math.abs(sd) * 4) + 1;
      final double[] kernel = new double[2 * radius + 1];
      final double norm = -0.5 / (sd * sd);
      for (int i = 0, j = radius, jj = radius; j < kernel.length; i++, j++, jj--) {
        kernel[j] = kernel[jj] = FastMath.exp(norm * i * i);
      }
      // Normalise
      double sum = 0;
      for (int j = 0; j < kernel.length; j++) {
        sum += kernel[j];
      }
      for (int j = 0; j < kernel.length; j++) {
        kernel[j] /= sum;
      }

      gg = Convolution.convolveFast(g, kernel);
      // The convolution will have created a larger array so we must adjust the offset for this
      c0 -= radius;
    } else {
      gg = g;
    }

    // Pad with constant c0
    final double[] g0 = new double[gg.length + c0];
    for (int c = 0; c < gg.length; c++) {
      if (c + c0 >= 0) {
        g0[c + c0] = gg[c];
      }
    }
    return g0;
  }

  private static int[] limits(int[] histogram) {
    int min = 0;
    while (histogram[min] == 0) {
      min++;
    }
    int max = histogram.length - 1;
    while (histogram[max] == 0) {
      max--;
    }
    return new int[] {min, max};
  }

  private static double[] getX(int[] limits) {
    final int min = 0; // limits[0]
    final int range = limits[1] - min + 1;
    final double[] x = new double[range];
    for (int i = 0; i < range; i++) {
      x[i] = (double) min + i;
    }
    return x;
  }

  private static double[] getY(int[] histogram, int[] limits) {
    final int min = 0; // limits[0]
    final int range = limits[1] - min + 1;
    final double[] y = new double[range];
    double sum = 0;
    for (int i = 0; i < range; i++) {
      y[i] = histogram[min + i];
      sum += y[i];
    }
    for (int i = 0; i < range; i++) {
      y[i] /= sum;
    }
    return y;
  }

  private static int getSize(int[] histogram) {
    int size = 0;
    for (int i = 0; i < histogram.length; i++) {
      size += histogram[i];
    }
    return size;
  }

  private static double getMean(int[] histogram) {
    int size = 0;
    double sum = 0;
    for (int i = 0; i < histogram.length; i++) {
      size += histogram[i];
      sum += histogram[i] * i;
    }
    return sum / size;
  }

  private int showDialog() {
    final GenericDialog gd = new GenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("em-gain-analysis"));

    gd.addMessage(
        "Analyse the white-light histogram of an image stack to determine EM-gain parameters.\n \n"
            + "See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321"
            + " (Supplementary Information).");

    if (extraOptions) {
      gd.addCheckbox("Simulate", settings.settingSimulate);
      gd.addNumericField("Bias", settings.settingBias, 0);
      gd.addNumericField("Gain", settings.settingGain, 2);
      gd.addNumericField("Noise", settings.settingNoise, 2);
      gd.addNumericField("Photons", settings.settingPhotons, 2);
      gd.addNumericField("Samples", settings.simulationSize, 0);
      gd.addCheckbox("Sample_PDF", settings.usePdf);
    }

    gd.addNumericField("Bias_estimate", settings.bias, 0);
    gd.addNumericField("Gain_estimate", settings.gain, 2);
    gd.addNumericField("Noise_estimate", settings.noise, 2);
    gd.addCheckbox("Show_approximation", settings.showApproximation);
    gd.showDialog();

    if (gd.wasCanceled()) {
      return DONE;
    }

    settings.save();

    if (extraOptions) {
      settings.settingSimulate = gd.getNextBoolean();
      settings.settingBias = gd.getNextNumber();
      settings.settingGain = gd.getNextNumber();
      settings.settingNoise = Math.abs(gd.getNextNumber());
      settings.settingPhotons = Math.abs(gd.getNextNumber());
      settings.simulationSize = (int) Math.abs(gd.getNextNumber());
      settings.usePdf = gd.getNextBoolean();
      if (gd.invalidNumber() || settings.settingBias < 0 || settings.settingGain < 1
          || settings.settingPhotons == 0 || settings.simulationSize == 0) {
        return DONE;
      }
    }

    settings.bias = gd.getNextNumber();
    settings.gain = gd.getNextNumber();
    settings.noise = Math.abs(gd.getNextNumber());
    settings.showApproximation = gd.getNextBoolean();

    if (gd.invalidNumber() || settings.bias < 0 || settings.gain < 1) {
      return DONE;
    }

    return (settings.settingSimulate) ? NO_IMAGE_REQUIRED : FLAGS;
  }

  @SuppressWarnings("unused")
  private void plotPmf() {
    if (!showPmfDialog()) {
      return;
    }

    final double step =
        getStepSize(settings.settingPhotons, settings.settingGain, settings.settingNoise);

    final Pdf pdf =
        pdf(0, step, settings.settingPhotons, settings.settingGain, settings.settingNoise);
    double[] pmf = pdf.probability;
    double yMax = MathUtils.max(pmf);

    // Get the approximation
    LikelihoodFunction fun;
    switch (settings.approximationType) {
      case 3:
        fun = new PoissonFunction(1.0 / settings.settingGain);
        break;
      case 2:
        // Use adaptive normalisation
        fun = PoissonGaussianFunction2.createWithStandardDeviation(1.0 / settings.settingGain,
            settings.settingNoise);
        break;
      case 1:
        // Create Poisson-Gamma (no Gaussian noise)
        fun = createPoissonGammaGaussianFunction(0);
        break;
      case 0:
      default:
        fun = createPoissonGammaGaussianFunction(settings.settingNoise);
    }
    double expected = settings.settingPhotons;
    if (settings.settingOffset != 0) {
      expected += settings.settingOffset * expected / 100.0;
    }

    // Normalise
    final boolean normalise = false;
    if (normalise) {
      final double sum = MathUtils.sum(pmf);
      for (int i = pmf.length; i-- > 0;) {
        pmf[i] /= sum;
      }
    }

    // Get CDF
    double sum = 0;
    double sum2 = 0;
    double[] x = pdf.x;
    double[] fvalues = new double[x.length];
    double[] cdf1 = new double[pmf.length];
    double[] cdf2 = new double[pmf.length];
    for (int i = 0; i < cdf1.length; i++) {
      sum += pmf[i] * step;
      cdf1[i] = sum;
      fvalues[i] = fun.likelihood(x[i], expected);
      sum2 += fvalues[i] * step;
      cdf2[i] = sum2;
    }

    // Truncate x for plotting
    int max = 0;
    double plimit = 1 - settings.tail;
    while (sum < plimit && max < pmf.length) {
      sum += pmf[max] * step;
      if (sum > 0.5 && pmf[max] == 0) {
        break;
      }
      max++;
    }

    int min = pmf.length;
    sum = 0;
    plimit = 1 - settings.head;
    while (sum < plimit && min > 0) {
      min--;
      sum += pmf[min] * step;
      if (sum > 0.5 && pmf[min] == 0) {
        break;
      }
    }

    pmf = Arrays.copyOfRange(pmf, min, max);
    x = Arrays.copyOfRange(x, min, max);
    fvalues = Arrays.copyOfRange(fvalues, min, max);

    if (settings.showApproximation) {
      yMax = MathUtils.maxDefault(yMax, fvalues);
    }

    final String label =
        String.format("Gain=%s, noise=%s, photons=%s", MathUtils.rounded(settings.settingGain),
            MathUtils.rounded(settings.settingNoise), MathUtils.rounded(settings.settingPhotons));

    final Plot plot = new Plot("PMF", "ADUs", "p");
    plot.setLimits(x[0], x[x.length - 1], 0, yMax);
    plot.setColor(Color.red);
    plot.addPoints(x, pmf, Plot.LINE);
    if (settings.showApproximation) {
      plot.setColor(Color.blue);
      plot.addPoints(x, fvalues, Plot.LINE);
    }

    plot.setColor(Color.magenta);
    plot.drawLine(settings.settingPhotons * settings.settingGain, 0,
        settings.settingPhotons * settings.settingGain, yMax);
    plot.setColor(Color.black);
    plot.addLabel(0, 0, label);
    final PlotWindow win1 = ImageJUtils.display("PMF", plot);

    // Plot the difference between the actual and approximation
    final double[] delta = new double[fvalues.length];
    for (int i = 0; i < fvalues.length; i++) {
      if (pmf[i] == 0 && fvalues[i] == 0) {
        continue;
      }
      if (settings.relativeDelta) {
        delta[i] =
            DoubleEquality.relativeError(fvalues[i], pmf[i]) * Math.signum(fvalues[i] - pmf[i]);
      } else {
        delta[i] = fvalues[i] - pmf[i];
      }
    }

    final Plot plot2 =
        new Plot("PMF delta", "ADUs", (settings.relativeDelta) ? "Relative delta" : "delta");
    final double[] limits = MathUtils.limits(delta);
    plot2.setLimits(x[0], x[x.length - 1], limits[0], limits[1]);
    plot2.setColor(Color.red);
    plot2.addPoints(x, delta, Plot.LINE);
    plot2.setColor(Color.magenta);
    plot2.drawLine(settings.settingPhotons * settings.settingGain, limits[0],
        settings.settingPhotons * settings.settingGain, limits[1]);
    plot2.setColor(Color.black);
    plot2.addLabel(0, 0, label + ((settings.settingOffset == 0) ? ""
        : ", expected = " + MathUtils.rounded(expected / settings.settingGain)));
    final WindowOrganiser wo = new WindowOrganiser();
    final PlotWindow win2 = ImageJUtils.display("PMF delta", plot2, wo);

    if (wo.isNotEmpty()) {
      final Point p2 = win1.getLocation();
      p2.y += win1.getHeight();
      win2.setLocation(p2);
    }

    // Plot the CDF of each distribution.
    // Compute the Kolmogorov distance as the supremum (maximum)
    // difference between the two cumulative probability distributions.
    // https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test
    double kolmogorovDistance = 0;
    double xd = x[0];
    for (int i = 0; i < cdf1.length; i++) {
      final double dist = Math.abs(cdf1[i] - cdf2[i]);
      if (kolmogorovDistance < dist) {
        kolmogorovDistance = dist;
        xd = pdf.x[i];
      }
    }
    cdf1 = Arrays.copyOfRange(cdf1, min, max);
    cdf2 = Arrays.copyOfRange(cdf2, min, max);

    final Plot plot3 = new Plot("CDF", "ADUs", "p");
    yMax = 1.05;
    plot3.setLimits(x[0], x[x.length - 1], 0, yMax);
    plot3.setColor(Color.red);
    plot3.addPoints(x, cdf1, Plot.LINE);
    plot3.setColor(Color.blue);
    plot3.addPoints(x, cdf2, Plot.LINE);

    plot3.setColor(Color.magenta);
    plot3.drawLine(settings.settingPhotons * settings.settingGain, 0,
        settings.settingPhotons * settings.settingGain, yMax);
    plot3.drawDottedLine(xd, 0, xd, yMax, 2);
    plot3.setColor(Color.black);
    plot3.addLabel(0, 0,
        label + ", Kolmogorov distance = " + MathUtils.rounded(kolmogorovDistance) + " @ " + xd);
    plot3.addLegend("CDF\nApprox");
    final int size = wo.size();
    final PlotWindow win3 = ImageJUtils.display("CDF", plot3, wo);

    if (size != wo.size()) {
      final Point p2 = win1.getLocation();
      p2.x += win1.getWidth();
      win3.setLocation(p2);
    }
  }

  private LikelihoodFunction createPoissonGammaGaussianFunction(double noise) {
    final PoissonGammaGaussianFunction fun =
        new PoissonGammaGaussianFunction(1.0 / settings.settingGain, noise);
    fun.setMinimumProbability(0);
    return fun;
  }

  /**
   * Gets the step size.
   *
   * <p>Note: Currently this just returns 1 as it should be a PMF so only uses discrete values.
   *
   * @param photons the photons
   * @param gain the gain
   * @param noise the noise
   * @return the step size
   */
  private static double getStepSize(double photons, double gain, double noise) {
    // Note: This is not valid as the PMF should only accept integer input.
    return 1;

    // // Determine the best step to plot the PMF.
    // // Ensure there are enough points on the chart.
    // double step = 1.0;
    //
    // int poissonGammaWidth = pdfEMGain(step, photons, gain).length;
    //
    // while (step > 0.001)
    // {
    // int gaussianWidth = (int) Math.ceil(Math.abs(noise) * 4 / step);
    // if (gaussianWidth * 2 + poissonGammaWidth / step < 300)
    // step /= 2;
    // else
    // break;
    // }
    // return step;
  }

  private boolean showPmfDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("em-gain-pmf"));

    gd.addMessage("Plot the probability mass function for EM-gain");

    gd.addNumericField("Gain", settings.settingGain, 2, 6, "Count/electron");
    gd.addNumericField("Noise", settings.settingNoise, 2, 6, "Count");
    gd.addNumericField("Photons", settings.settingPhotons, 2);
    gd.addChoice("Approx", Settings.APPROXIMATION_TYPE, settings.approximationType);
    gd.addCheckbox("Show_approximation", settings.showApproximation);
    if (extraOptions) {
      gd.addNumericField("Approximation_offset (%)", settings.settingOffset, 2);
    }
    gd.addNumericField("Remove_head", settings.head, 3);
    gd.addNumericField("Remove_tail", settings.tail, 3);
    gd.addCheckbox("Relative_delta", settings.relativeDelta);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.settingGain = gd.getNextNumber();
    settings.settingNoise = Math.abs(gd.getNextNumber());
    settings.settingPhotons = Math.abs(gd.getNextNumber());
    settings.approximationType = gd.getNextChoiceIndex();
    settings.showApproximation = gd.getNextBoolean();
    if (extraOptions) {
      settings.settingOffset = gd.getNextNumber();
    }
    settings.head = Math.abs(gd.getNextNumber());
    settings.tail = Math.abs(gd.getNextNumber());
    settings.relativeDelta = gd.getNextBoolean();
    settings.save();

    return (!gd.invalidNumber() && settings.settingBias >= 0 && settings.settingGain >= 1
        && settings.settingPhotons != 0 && settings.tail <= 0.5 && settings.head <= 0.5);
  }
}
