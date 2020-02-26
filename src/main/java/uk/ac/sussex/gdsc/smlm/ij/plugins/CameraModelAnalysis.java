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

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.ImageProcessor;
import java.awt.AWTEvent;
import java.awt.Color;
import java.util.Arrays;
import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.exception.MathArithmeticException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.NormalizedGaussianSampler;
import org.apache.commons.rng.sampling.distribution.PoissonSampler;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog.OptionCollectedEvent;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog.OptionCollectedListener;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog.OptionListener;
import uk.ac.sussex.gdsc.core.ij.gui.NonBlockingExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.Plot2;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.math.GeometryUtils;
import uk.ac.sussex.gdsc.core.threshold.IntHistogram;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.rng.MarsagliaTsangGammaSampler;
import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;
import uk.ac.sussex.gdsc.core.utils.rng.UniformRandomProviders;
import uk.ac.sussex.gdsc.smlm.data.NamedObject;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.CameraModelAnalysisSettings;
import uk.ac.sussex.gdsc.smlm.function.InterpolatedPoissonFunction;
import uk.ac.sussex.gdsc.smlm.function.LikelihoodFunction;
import uk.ac.sussex.gdsc.smlm.function.LogFactorial;
import uk.ac.sussex.gdsc.smlm.function.PoissonFunction;
import uk.ac.sussex.gdsc.smlm.function.PoissonGammaFunction;
import uk.ac.sussex.gdsc.smlm.function.PoissonGammaGaussianConvolutionFunction;
import uk.ac.sussex.gdsc.smlm.function.PoissonGammaGaussianFunction;
import uk.ac.sussex.gdsc.smlm.function.PoissonGammaGaussianFunction.ConvolutionMode;
import uk.ac.sussex.gdsc.smlm.function.PoissonGaussianConvolutionFunction;
import uk.ac.sussex.gdsc.smlm.function.PoissonGaussianFunction2;
import uk.ac.sussex.gdsc.smlm.function.PoissonPoissonFunction;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.math3.analysis.integration.CustomSimpsonIntegrator;
import uk.ac.sussex.gdsc.smlm.utils.Convolution;
import uk.ac.sussex.gdsc.smlm.utils.GaussianKernel;

/**
 * Model the on-chip amplification from an EM-CCD camera, CCD or sCMOS camera.
 */
public class CameraModelAnalysis
    implements ExtendedPlugInFilter, DialogListener, OptionCollectedListener {
  // Add mode to compute the distance over a range of photons at a set gain and variance
  // Compute distance to simulation. Compute distance to another distribution.

  private static final String TITLE = "Camera Model Analysis";
  private static final KolmogorovSmirnovTest kolmogorovSmirnovTest = new KolmogorovSmirnovTest();

  /** The relative change in cumulative probability to stop computation. */
  private static final double PROBABILITY_DELTA = 1e-6;
  /** The lower limit on the cumulative probability. */
  private static final double LOWER = 1e-6;
  /** The upper limit on the cumulative probability. */
  private static final double UPPER = 1 - LOWER;

  private static final double SIMPSON_RELATIVE_ACCURACY = 1e-4;
  private static final double SIMPSON_ABSOLUTE_ACCURACY = 1e-8;
  private static final int SIMPSON_MIN_ITERATION_COUNT = 3;

  // Use Simpson's integration with n=4 to get the integral of the probability
  // over the range of each count.
  private static final int SIMPSON_N = 4;
  private static final int SIMPSON_N_2 = SIMPSON_N / 2;
  private static final double SIMPSON_H = 1.0 / SIMPSON_N;

  private static final String[] MODE = {"CCD", "EM-CCD", "sCMOS"};
  private static final int MODE_CCD = 0;
  private static final int MODE_EM_CCD = 1;
  private static final int MODE_SCMOS = 2;

  private static final String[] MODEL = SettingsManager.getNames((Object[]) Model.values());

  private CameraModelAnalysisSettings.Builder settings;

  private boolean extraOptions;
  private boolean dirty = true;
  private CameraModelAnalysisSettings lastSettings;
  private ExtendedGenericDialog gd;
  private IntHistogram lastHistogram;
  private double[][] floatHistogram;
  private CameraModelAnalysisSettings lastSimulationSettings;

  private enum Model implements NamedObject {
    ///////////////
    // CCD / sCMOS
    ///////////////

    //@formatter:off

    POISSON_PMF { @Override
    public String getName() { return "Poisson PMF"; } },
    POISSON_DISRECTE { @Override
    public String getName() { return "Poisson (Discrete)"; } },
    POISSON_CONTINUOUS { @Override
    public String getName() { return "Poisson (Continuous)"; } },
    POISSON_GAUSSIAN_PDF { @Override
    public String getName() { return "Poisson+Gaussian PDF integration"; } },

    // Best for CCD/sCMOS
    POISSON_GAUSSIAN_PMF { @Override
    public String getName() { return "Poisson+Gaussian PMF integration"; } },
    // Saddle-point approximation.
    // Very good. Relatively worse than POISSON_GAUSSIAN_PMF at very low photons.
    POISSON_GAUSSIAN_APPROX { @Override
    public String getName() { return "Poisson+Gaussian approximation"; } },
     // Mixed Poisson distribution (Noise is added as a second Poisson)
    POISSON_POISSON { @Override
    public String getName() { return "Poisson+Poisson"; } },

    ///////////////
    // EMCCD
    ///////////////
    // There is no obvious best for EM-CCD:
    // This requires a full range test to determine the best function for which
    // parameters.

    // Good when no noise.
    // Under-estimates total probability when gain is low (<15) and photons are low (<2).
    POISSON_GAMMA_PMF { @Override
    public String getName() { return "Poisson+Gamma PMF"; } },

    // Good but relatively worse as the read noise increases.
    // Requires full integration when read noise is low (<1) and photons are low (<5).
    POISSON_GAMMA_GAUSSIAN_APPROX { @Override
    public String getName() { return "Poisson+Gamma+Gaussian approximation"; } },
    // Good when read noise is >>1.
    // Requires full integration when read noise is low (<1).
    POISSON_GAMMA_GAUSSIAN_PDF_INTEGRATION { @Override
    public String getName() { return "Poisson+Gamma+Gaussian PDF integration"; } },
    // Good
    // Slow
    POISSON_GAMMA_GAUSSIAN_PMF_INTEGRATION { @Override
    public String getName() { return "Poisson+Gamma+Gaussian PMF integration"; } },
    // Best for EM-CCD.
    // Very robust (computes the full convolution
    // of the Gaussian and the Poisson-Gamma plus the delta function PMF contribution).
    POISSON_GAMMA_GAUSSIAN_SIMPSON_INTEGRATION { @Override
    public String getName() { return "Poisson+Gamma+Gaussian Simpson's integration"; } },
    // Good
    POISSON_GAMMA_GAUSSIAN_LEGENDRE_GAUSS_INTEGRATION { @Override
    public String getName() { return "Poisson+Gamma+Gaussian Legendre-Gauss integration"; } },
    // Independent implementation of POISSON_GAMMA_GAUSSIAN_PDF_INTEGRATION.
    // Good when read noise is >>1.
    // Requires full integration when read noise is low (<1).
    POISSON_GAMMA_GAUSSIAN_PDF_CONVOLUTION { @Override
    public String getName() { return "Poisson+Gamma*Gaussian convolution"; } },
    ;
    //@formatter:on

    @Override
    public String getShortName() {
      return getName();
    }

    public static Model forNumber(int number) {
      final Model[] values = Model.values();
      if (number < 0 || number >= values.length) {
        number = 0;
      }
      return values[number];
    }
  }

  private static interface Round {
    int round(double value);
  }

  private static class RoundDown implements Round {
    private static final Round INSTANCE = new RoundDown();

    @Override
    public int round(double value) {
      return (int) Math.floor(value);
    }
  }

  private static class RoundUpDown implements Round {
    private static final Round INSTANCE = new RoundUpDown();

    @Override
    public int round(double value) {
      return (int) Math.round(value);
    }
  }

  private static Round getRound(CameraModelAnalysisSettings settings) {
    return (settings.getRoundDown()) ? RoundDown.INSTANCE : RoundUpDown.INSTANCE;
  }

  private UniformRandomProvider getRandomGenerator() {
    return UniformRandomProviders.create(settings.getSeed());
  }

  @Override
  public int setup(String arg, ImagePlus imp) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);
    extraOptions = ImageJUtils.isExtraOptions();
    return NO_IMAGE_REQUIRED;
  }

  @Override
  public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr) {
    settings = SettingsManager.readCameraModelAnalysisSettings(0).toBuilder();

    gd = new NonBlockingExtendedGenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);

    gd.addMessage("Simulate on-chip camera applification.");

    gd.addNumericField("Photons", settings.getPhotons(), 2);
    gd.addChoice("Mode", MODE, settings.getMode(), new OptionListener<Integer>() {
      @Override
      public boolean collectOptions(Integer value) {
        settings.setMode(value);
        return collectOptions(false);
      }

      @Override
      public boolean collectOptions() {
        return collectOptions(true);
      }

      private boolean collectOptions(boolean silent) {
        final int mode = settings.getMode();
        final ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE);
        if (mode == MODE_CCD) {
          egd.addNumericField("Gain", settings.getGain(), 2, 6, "Count/electrons");
          egd.addNumericField("Noise", settings.getNoise(), 2, 6, "Count");
        } else if (mode == MODE_EM_CCD) {
          egd.addNumericField("Gain", settings.getEmGain(), 2, 6, "Count/electrons");
          egd.addNumericField("Noise", settings.getEmNoise(), 2, 6, "Count");
          egd.addNumericField("EM_samples", settings.getEmSamples(), 0);
        } else if (mode == MODE_SCMOS) {
          egd.addNumericField("Gain", settings.getCmosGain(), 2, 6, "Count/electrons");
          egd.addNumericField("Noise", settings.getCmosNoise(), 2, 6, "Count");
        } else {
          throw new IllegalStateException();
        }
        egd.setSilent(silent);
        egd.showDialog(true, gd);
        if (egd.wasCanceled()) {
          return false;
        }
        if (mode == MODE_CCD) {
          settings.setGain(egd.getNextNumber());
          settings.setNoise(egd.getNextNumber());
        } else if (mode == MODE_EM_CCD) {
          settings.setEmGain(egd.getNextNumber());
          settings.setEmNoise(egd.getNextNumber());
          settings.setEmSamples(Math.max(1, (int) egd.getNextNumber()));
        } else {
          // MODE_SCMOS
          settings.setCmosGain(egd.getNextNumber());
          settings.setCmosNoise(egd.getNextNumber());
        }
        return true;
      }
    });
    if (extraOptions) {
      gd.addNumericField("Seed", settings.getSeed(), 0);
    }
    gd.addNumericField("Samples", settings.getSamples(), 0);
    gd.addNumericField("Noise_samples", settings.getNoiseSamples(), 0);
    gd.addCheckbox("Round_down", settings.getRoundDown());
    gd.addChoice("Model", MODEL, settings.getModel());
    gd.addCheckbox("Full_integration", settings.getSimpsonIntegration());
    gd.addOptionCollectedListener(this);
    gd.addDialogListener(this);
    gd.addPreviewCheckbox(pfr);
    gd.showDialog();

    SettingsManager.writeSettings(settings);

    if (!gd.wasCanceled() && dirty) {
      execute();
    }

    return DONE;
  }

  @Override
  public boolean dialogItemChanged(GenericDialog gd, AWTEvent event) {
    dirty = true;
    settings.setPhotons(gd.getNextNumber());
    settings.setMode(gd.getNextChoiceIndex());
    if (extraOptions) {
      settings.setSeed((int) gd.getNextNumber());
    }
    settings.setSamples(Math.max(1, (int) gd.getNextNumber()));
    settings.setNoiseSamples(Math.max(1, (int) gd.getNextNumber()));
    settings.setRoundDown(gd.getNextBoolean());
    settings.setModel(gd.getNextChoiceIndex());
    settings.setSimpsonIntegration(gd.getNextBoolean());
    if (gd.getPreviewCheckbox().getState()) {
      return execute();
    }
    return true;
  }

  @Override
  public void optionCollected(OptionCollectedEvent event) {
    if (gd.getPreviewCheckbox().getState()) {
      boolean ok = false;
      try {
        ok = execute();
      } catch (final Exception ex) {
        // Catch as this is run within a AWT dispatch thread
        ImageJUtils.log(TITLE + "Error: " + ex.getMessage());
      } finally {
        if (!ok) {
          gd.getPreviewCheckbox().setState(false);
        }
      }
    }
  }

  @Override
  public void setNPasses(int passes) {
    // Ignore
  }

  @Override
  public void run(ImageProcessor ip) {
    execute();
  }

  /**
   * Execute the analysis.
   */
  private boolean execute() {
    dirty = false;

    final CameraModelAnalysisSettings settings = this.settings.build();
    if (!(getGain(settings) > 0)) {
      ImageJUtils.log(TITLE + "Error: No total gain");
      return false;
    }
    if (!(settings.getPhotons() > 0)) {
      ImageJUtils.log(TITLE + "Error: No photons");
      return false;
    }
    // Avoid repeating the same analysis
    if (settings.equals(lastSettings)) {
      return true;
    }
    lastSettings = settings;

    final IntHistogram h = getHistogram(settings);

    // Build cumulative distribution
    final double[][] cdf1 = cumulativeHistogram(h);
    final double[] x1 = cdf1[0];
    final double[] y1 = cdf1[1];

    // Interpolate to 300 steps faster evaluation?

    // Get likelihood function
    final LikelihoodFunction f = getLikelihoodFunction(settings);

    // Create likelihood cumulative distribution
    final double[][] cdf2 = cumulativeDistribution(settings, cdf1, f);

    // Compute Komolgorov distance
    final double[] distanceAndValue = getDistance(cdf1, cdf2);
    final double distance = distanceAndValue[0];
    final double value = distanceAndValue[1];
    final double area = distanceAndValue[2];

    final double[] x2 = cdf2[0];
    final double[] y2 = cdf2[1];

    // Fill y1
    int offset = 0;
    while (x2[offset] < x1[0]) {
      offset++;
    }
    final double[] y1b = new double[y2.length];
    System.arraycopy(y1, 0, y1b, offset, y1.length);
    Arrays.fill(y1b, offset + y1.length, y2.length, y1[y1.length - 1]);

    // KolmogorovSmirnovTest
    // n is the number of samples used to build the probability distribution.
    final int n = (int) MathUtils.sum(h.histogramCounts);

    // From KolmogorovSmirnovTest.kolmogorovSmirnovTest(RealDistribution distribution, double[]
    // data, boolean exact):
    // Returns the p-value associated with the null hypothesis that data is a sample from
    // distribution.
    // E.g. If p<0.05 then the null hypothesis is rejected and the data do not match the
    // distribution.
    double pvalue = Double.NaN;
    try {
      pvalue = 1d - kolmogorovSmirnovTest.cdf(distance, n);
    } catch (final MathArithmeticException ex) {
      // Cannot be computed to leave at NaN
    }

    // Plot
    final WindowOrganiser wo = new WindowOrganiser();
    String title = TITLE + " CDF";
    Plot2 plot = new Plot2(title, "Count", "CDF");
    final double max = 1.05 * MathUtils.maxDefault(1, y2);
    plot.setLimits(x2[0], x2[x2.length - 1], 0, max);
    plot.setColor(Color.blue);
    plot.addPoints(x2, y1b, Plot2.BAR);
    plot.setColor(Color.red);
    plot.addPoints(x2, y2, Plot2.BAR);
    plot.setColor(Color.magenta);
    plot.drawLine(value, 0, value, max);
    plot.setColor(Color.black);
    plot.addLegend("CDF\nModel");
    plot.addLabel(0, 0, String.format("Distance=%s @ %.0f (Mean=%s) : p=%s",
        MathUtils.rounded(distance), value, MathUtils.rounded(area), MathUtils.rounded(pvalue)));
    ImageJUtils.display(title, plot, ImageJUtils.NO_TO_FRONT, wo);

    // Show the histogram
    title = TITLE + " Histogram";
    plot = new Plot2(title, "Count", "Frequency");
    // Update X1 so that the histogram bars are centred over the x value
    for (int i = x1.length; i-- > 0;) {
      x1[i] -= 0.5;
    }
    plot.setLimits(x1[0] - 0.5, x1[x1.length - 1] + 1.5, 0,
        MathUtils.max(h.histogramCounts) * 1.05);
    plot.setColor(Color.blue);
    plot.addPoints(x1, SimpleArrayUtils.toDouble(h.histogramCounts), Plot2.BAR);

    plot.setColor(Color.red);
    final double[] x = floatHistogram[0].clone();
    final double[] y = floatHistogram[1].clone();
    final double scale = n / (MathUtils.sum(y) * (x[1] - x[0]));
    for (int i = 0; i < y.length; i++) {
      y[i] *= scale;
    }
    plot.addPoints(x, y, Plot.LINE);

    plot.setColor(Color.black);
    plot.addLegend("Sample\nExpected");
    ImageJUtils.display(title, plot, ImageJUtils.NO_TO_FRONT, wo);

    wo.tile();

    return true;
  }

  private IntHistogram getHistogram(CameraModelAnalysisSettings settings) {
    if (lastHistogram == null || newSimulationSettings(settings, lastSimulationSettings)) {
      IJ.showStatus("Simulating histogram ...");
      final StopWatch sw = StopWatch.createStarted();
      lastHistogram = simulateHistogram(settings, getRandomGenerator());
      lastSimulationSettings = settings;
      IJ.showStatus("Simulated in " + sw.toString());

      // Convolve
      floatHistogram = convolveHistogram(settings);
    }
    return lastHistogram;
  }

  private static boolean newSimulationSettings(CameraModelAnalysisSettings s1,
      CameraModelAnalysisSettings s2) {
    // Check those settings for the simulation
    if (s1.getPhotons() != s2.getPhotons() || s1.getMode() != s2.getMode()
        || s1.getSamples() != s2.getSamples() || s1.getNoiseSamples() != s2.getNoiseSamples()
        || s1.getSeed() != s2.getSeed() || s1.getRoundDown() != s2.getRoundDown()) {
      return true;
    }
    if (s1.getMode() == MODE_CCD) {
      return (s1.getGain() != s2.getGain() || s1.getNoise() != s2.getNoise());
    } else if (s1.getMode() == MODE_SCMOS) {
      return (s1.getCmosGain() != s2.getCmosGain() || s1.getCmosNoise() != s2.getCmosNoise());
    }
    // MODE_EM_CCD
    return (s1.getEmGain() != s2.getEmGain() || s1.getEmNoise() != s2.getEmNoise()
        || s1.getEmSamples() != s2.getEmSamples());
  }

  /**
   * Simulate the histogram for fitting.
   *
   * @param settings the settings
   * @param random the random
   * @return The histogram
   */
  private static IntHistogram simulateHistogram(CameraModelAnalysisSettings settings,
      UniformRandomProvider random) {
    switch (settings.getMode()) {
      // EM-CCD
      case 1:
        return simulatePoissonGammaGaussian(settings, random);

      // CCD or sCMOS
      case 2:
      case 0:
        return simulatePoissonGaussian(settings, random);

      default:
        throw new IllegalStateException();
    }
  }

  /**
   * Randomly generate a histogram from poisson-gamma-gaussian samples.
   *
   * @param settings the settings
   * @param random the random
   * @return The histogram
   */
  private static IntHistogram simulatePoissonGammaGaussian(CameraModelAnalysisSettings settings,
      UniformRandomProvider random) {
    final int[] poissonSample = simulatePoisson(settings, random);

    final NormalizedGaussianSampler gauss = SamplerUtils.createNormalizedGaussianSampler(random);

    final double gain = getGain(settings);

    // Note that applying a separate EM-gain and then the camera gain later
    // is the same as applying the total gain in the gamma distribution and no camera gain
    // later, i.e. the Gamma distribution is just squeezed.
    // Thus we sample from a Gamma distribution with a fixed gain (the scale) and the
    // number of photons is the shape.

    final double noise = getReadNoise(settings);
    final int samples = settings.getSamples();
    final int emSamples = settings.getEmSamples();
    final int noiseSamples = (noise > 0) ? settings.getNoiseSamples() : 1;
    final int[] sample = new int[samples * emSamples * noiseSamples];
    int count = 0;
    final Round round = getRound(settings);
    if (gain < 1) {
      throw new IllegalStateException("EM-gain must be >= 1: " + gain);
    }
    final MarsagliaTsangGammaSampler gamma = new MarsagliaTsangGammaSampler(random, 1, gain);
    for (int n = poissonSample.length; n-- > 0;) {
      if (poissonSample[n] != 0) {
        // Gamma
        gamma.setAlpha(poissonSample[n]);

        // Q. The Poisson-Gamma Function does not exactly match this
        // when the mean is around 1. The model from Ulbrich+Isacoff
        // is using a gamma function. This returns a continuous value.
        // Should it be rounded to make it a discrete PMF?
        // Or should the model actually be integrated.
        // Basically the value of the Poisson-Gamma from 0-0.5 should be
        // added to the delta function at c=0 for the count at c=0.
        // This would fix the function.

        // Should we use the Tubb's model which uses:
        // final double shape = count
        // final double scale = gain - 1 + 1 / shape
        // final double electrons = random.nextGamma(shape, scale) - 1
        // final double output = count + electrons

        // The Tubb's model is for additional electrons. So a count of 1
        // can never generate an output of 0. This does not fit the
        // Ulbrich+Isacoff model where c=0 is actually defined, i.e. you can
        // model zero output for EM-gain.

        // Over-sample the Gamma
        for (int k = emSamples; k-- > 0;) {
          // Do rounding to simulate a discrete PMF.
          final double d2 = round.round(gamma.sample());

          // Over-sample the Gaussian
          for (int j = noiseSamples; j-- > 0;) {
            // Convert the sample to a count
            sample[count++] = round.round(d2 + noise * gauss.sample());
          }
        }
      } else {
        // Still over-sample the Gamma even though it was zero
        for (int k = emSamples; k-- > 0;) {
          // Over-sample the Gaussian
          for (int j = noiseSamples; j-- > 0;) {
            // Convert the sample to a count
            sample[count++] = round.round(noise * gauss.sample());
          }
        }
      }
    }

    return createHistogram(sample, count);
  }

  private static int[] simulatePoisson(CameraModelAnalysisSettings settings,
      UniformRandomProvider random) {
    return SamplerUtils.createSamples(settings.getSamples(),
        new PoissonSampler(random, settings.getPhotons()));
  }

  /**
   * Randomly generate a histogram from poisson-gaussian samples.
   *
   * @param settings the settings
   * @param random the random
   * @return The histogram
   */
  private static IntHistogram simulatePoissonGaussian(CameraModelAnalysisSettings settings,
      UniformRandomProvider random) {
    final int[] poissonSample = simulatePoisson(settings, random);

    final NormalizedGaussianSampler gauss = SamplerUtils.createNormalizedGaussianSampler(random);

    final double gain = getGain(settings);
    final double noise = getReadNoise(settings);
    final int samples = settings.getSamples();
    final int noiseSamples = (noise > 0) ? settings.getNoiseSamples() : 1;
    final int[] sample = new int[samples * noiseSamples];
    int count = 0;
    final Round round = getRound(settings);
    for (int n = poissonSample.length; n-- > 0;) {
      // Fixed camera gain
      final double d = poissonSample[n] * gain;

      // Over-sample the Gaussian
      for (int j = noiseSamples; j-- > 0;) {
        // Convert the sample to a count
        sample[count++] = round.round(d + noise * gauss.sample());
      }
    }

    return createHistogram(sample, count);
  }

  private static IntHistogram createHistogram(int[] sample, int count) {
    final int[] limits = MathUtils.limits(sample);
    final int min = limits[0];
    final int max = limits[1];

    final int[] hist = new int[max - min + 1];
    for (int i = count; i-- > 0;) {
      hist[sample[i] - min]++;
    }
    return new IntHistogram(hist, min);
  }

  /**
   * Convolve the histogram. The output is a discrete probability distribution.
   *
   * @param settings the settings
   * @return The histogram
   */
  private static double[][] convolveHistogram(CameraModelAnalysisSettings settings) {

    // Find the range of the Poisson
    final PoissonDistribution poisson = new PoissonDistribution(null, settings.getPhotons(),
        PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);
    final int maxn = poisson.inverseCumulativeProbability(UPPER);

    final double gain = getGain(settings);
    final double noise = getReadNoise(settings);

    final boolean debug = false;

    // Build the Probability Mass/Density Function (PDF) of the distribution:
    // either a Poisson (PMF) or Poisson-Gamma (PDF). The PDF is 0 at all
    // values apart from the step interval.
    // Note: The Poisson-Gamma is computed without the Dirac delta contribution
    // at c=0. This allows correct convolution with the Gaussian of the dirac delta
    // and the rest of the Poisson-Gamma (so matching the simulation).
    final TDoubleArrayList list = new TDoubleArrayList();
    double step;
    String name;

    int upsample = 100;

    // Store the Dirac delta value at c=0. This must be convolved separately.
    double dirac = 0;

    // EM-CCD
    if (settings.getMode() == MODE_EM_CCD) {
      name = "Poisson-Gamma";

      final double m = gain;
      final double p = settings.getPhotons();

      dirac = PoissonGammaFunction.dirac(p);

      // Chose whether to compute a discrete PMF or a PDF using the approximation.
      // Note: The delta function at c=0 is from the PMF of the Poisson. So it is
      // a discrete contribution. This is omitted from the PDF and handled in
      // a separate convolution.
      final boolean discrete = false; // noise != 0;
      if (discrete) {
        // Note: This is obsolete as the Poisson-Gamma function is continuous.
        // Sampling it at integer intervals is not valid, especially for low gain.
        // The Poisson-Gamma PDF should be integrated to form a discrete PMF.

        step = 1.0;

        double upper;
        if (settings.getPhotons() < 20) {
          upper = maxn;
        } else {
          // Approximate reasonable range of Poisson as a Gaussian
          upper = settings.getPhotons() + 3 * Math.sqrt(settings.getPhotons());
        }

        final GammaDistribution gamma = new GammaDistribution(null,
            upper, gain, GammaDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);
        final int maxc = (int) gamma.inverseCumulativeProbability(0.999);

        final int minn = Math.max(1, poisson.inverseCumulativeProbability(LOWER));

        // See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.

        // Note this is not a convolution of a single Gamma distribution since the shape
        // is modified not the count. So it is a convolution of a distribution made with
        // a gamma of fixed count and variable shape.

        // The count=0 is a special case.
        list.add(PoissonGammaFunction.poissonGammaN(0, p, m));

        final long total = (maxn - minn) * (long) maxc;

        if (total < 1000000) {
          // Full computation

          // G(c) = sum n { (1 / n!) p^n e^-p (1 / ((n-1!)m^n)) c^n-1 e^-c/m }
          // Compute as a log
          // - log(n!) + n*log(p)-p -log((n-1)!) - n * log(m) + (n-1) * log(c) -c/m

          // Note: Both methods work

          LogFactorial.increaseTableMaxN(maxn);
          final double[] f = new double[maxn + 1];
          final double logm = Math.log(m);
          final double logp = Math.log(p);
          for (int n = minn; n <= maxn; n++) {
            f[n] = -LogFactorial.logF(n) + n * logp - p - LogFactorial.logF(n - 1) - n * logm;
          }

          // Use Poisson + Gamma distribution
          // double[] pd = new double[maxn + 1];
          // CustomGammaDistribution[] gd = new CustomGammaDistribution[maxn + 1];
          // for (int n = minn; n <= maxn; n++)
          // {
          // pd[n] = poisson.probability(n);
          // gd[n] = new CustomGammaDistribution(null, n, m);
          // }

          // double total = list.getQuick(0);
          // double total2 = total;
          for (int c = 1; c <= maxc; c++) {
            double sum = 0;
            final double c_m = c / m;
            final double logc = Math.log(c);
            for (int n = minn; n <= maxn; n++) {
              sum += FastMath.exp(f[n] + (n - 1) * logc - c_m);
            }
            // sum2 += pd[n] * gd[n].density(c);
            list.add(sum);
            // total += sum;

            // This should match the approximation
            // double approx = PoissonGammaFunction.poissonGamma(c, p, m);
            // total2 += approx;
            // System.out.printf("c=%d sum=%g approx=%g error=%g\n", c, sum2, approx,
            // uk.ac.sussex.gdsc.core.utils.DoubleEquality.relativeError(sum2, approx));
          }

          // System.out.printf("sum=%g approx=%g error=%g\n", total, total2,
          // uk.ac.sussex.gdsc.core.utils.DoubleEquality.relativeError(total, total2));
        } else {
          // Approximate
          for (int c = 1; c <= maxc; c++) {
            list.add(PoissonGammaFunction.poissonGammaN(c, p, m));
          }
        }
      } else {
        // This integrates the PDF using the approximation and up-samples together.
        // Compute the sampling interval.
        step = 1.0 / upsample;
        upsample = 1; // Reset

        // Compute the integral of [-step/2:step/2] for each point.
        // Use trapezoid integration.
        final double step_2 = step / 2;
        double prev = PoissonGammaFunction.poissonGammaN(0, p, m);
        double next = PoissonGammaFunction.poissonGammaN(step_2, p, m);
        list.add((prev + next) * 0.25);
        double max = 0;
        for (int i = 1;; i++) {
          // Each remaining point is modelling a PMF for the range [-step/2:step/2]
          prev = next;
          next = PoissonGammaFunction.poissonGammaN(i * step + step_2, p, m);
          final double pp = (prev + next) * 0.5;
          if (max < pp) {
            max = pp;
          }
          if (pp / max < 1e-5) {
            // Use this if non-zero since it has been calculated
            if (pp != 0) {
              list.add(pp);
            }
            break;
          }
          list.add(pp);
        }
      }

      // Ensure the combined sum of PDF and Dirac is 1
      final double expected = 1 - dirac;
      // Compute the sum using Simpson's rule:
      // Require an odd number to get an even number (n) of sub-intervals:
      if (list.size() % 2 == 0) {
        list.add(0);
      }
      final double[] g = list.toArray();
      // Number of sub intervals
      final int n = g.length - 1;
      final double h = 1; // h = (a-b) / n = sub-interval width
      double sum2 = 0;
      double sum4 = 0;
      for (int j = 1; j <= n / 2 - 1; j++) {
        sum2 += g[2 * j];
      }
      for (int j = 1; j <= n / 2; j++) {
        sum4 += g[2 * j - 1];
      }
      final double sum = (h / 3) * (g[0] + 2 * sum2 + 4 * sum4 + g[n]);
      // Check
      // System.out.printf("Sum=%g Expected=%g\n", sum * step, expected);
      SimpleArrayUtils.multiply(g, expected / sum);
      list.resetQuick();
      list.add(g);
    } else {
      name = "Poisson";
      // Apply fixed gain. Just change the step interval of the PMF.
      step = gain;

      for (int n = 0; n <= maxn; n++) {
        list.add(poisson.probability(n));
      }
      final double p = poisson.probability(list.size());
      if (p != 0) {
        list.add(p);
      }
    }

    // Debug
    if (debug) {
      final String title = name;
      final Plot plot = new Plot(title, "x", "y");
      plot.addPoints(SimpleArrayUtils.newArray(list.size(), 0, step), list.toArray(), Plot.LINE);
      ImageJUtils.display(title, plot);
    }

    double zero = 0;
    double[] pg = list.toArray();

    // Sample Gaussian
    if (noise > 0) {
      step /= upsample;
      pg = list.toArray();

      // Convolve with Gaussian kernel
      final double[] kernel = GaussianKernel.makeGaussianKernel(Math.abs(noise) / step, 6, true);

      if (upsample != 1) {
        // Use scaled convolution. This is faster that zero filling distribution g.
        pg = Convolution.convolve(kernel, pg, upsample);
      } else if (dirac > 0.01) {
        // The Poisson-Gamma may be stepped at low mean causing wrap artifacts in the FFT.
        // This is a problem if most of the probability is in the Dirac.
        // Otherwise it can be ignored and the FFT version is OK.
        pg = Convolution.convolve(kernel, pg);
      } else {
        pg = Convolution.convolveFast(kernel, pg);
      }

      // The convolution will have created a larger array so we must adjust the offset for this
      final int radius = kernel.length / 2;
      zero -= radius * step;

      // Add convolution of the dirac delta function.
      if (dirac != 0) {
        // We only need to convolve the Gaussian at c=0
        for (int i = 0; i < kernel.length; i++) {
          pg[i] += kernel[i] * dirac;
        }
      }

      // Debug
      if (debug) {
        String title = "Gaussian";
        Plot plot = new Plot(title, "x", "y");
        plot.addPoints(SimpleArrayUtils.newArray(kernel.length, -radius * step, step), kernel,
            Plot.LINE);
        ImageJUtils.display(title, plot);

        title = name + "-Gaussian";
        plot = new Plot(title, "x", "y");
        plot.addPoints(SimpleArrayUtils.newArray(pg.length, zero, step), pg, Plot.LINE);
        ImageJUtils.display(title, plot);
      }

      zero = downSampleCdfToPmf(settings, list, step, zero, pg, 1.0);

      pg = list.toArray();
      zero = (int) Math.floor(zero);
      step = 1.0;

      // No convolution means we have the Poisson PMF/Poisson-Gamma PDF already
    } else if (step != 1) {
      // Sample to 1.0 pixel step interval.
      if (settings.getMode() == MODE_EM_CCD) {
        // Poisson-Gamma PDF
        zero = downSampleCdfToPmf(settings, list, step, zero, pg, 1 - dirac);
        pg = list.toArray();
        zero = (int) Math.floor(zero);

        // Add the dirac delta function.
        if (dirac != 0) {
          // Note: zero is the start of the x-axis. This value should be -1.
          assert (int) zero == -1;
          // Use as an offset to find the actual zero.
          pg[-(int) zero] += dirac;
        }
      } else {
        // Poisson PMF

        // Simple non-interpolated expansion.
        // This should be used when there is no Gaussian convolution.
        final double[] pd = pg;
        list.resetQuick();

        // Account for rounding.
        final Round round = getRound(settings);

        final int maxc = round.round(pd.length * step + 1);
        pg = new double[maxc];
        for (int n = pd.length; n-- > 0;) {
          pg[round.round(n * step)] += pd[n];
        }

        if (pg[0] != 0) {
          list.add(0);
          list.add(pg);
          pg = list.toArray();
          zero--;
        }
      }

      step = 1.0;
    } else {
      // Add the dirac delta function.
      list.setQuick(0, list.getQuick(0) + dirac);
    }

    return new double[][] {SimpleArrayUtils.newArray(pg.length, zero, step), pg};
  }

  private static double downSampleCdfToPmf(CameraModelAnalysisSettings settings,
      TDoubleArrayList list, double step, double zero, double[] pg, double sum) {
    // Down-sample to 1.0 pixel step interval.

    // Build cumulative distribution.
    double lowerSum = 0;
    for (int i = 0; i < pg.length; i++) {
      lowerSum += pg[i];
      pg[i] = lowerSum;
    }
    SimpleArrayUtils.multiply(pg, sum / lowerSum);
    pg[pg.length - 1] = sum;

    final double offset = (settings.getRoundDown()) ? 0 : -0.5;

    // Note the interpolation of the CDF is good when the step is much smaller than 1.
    // When the step is above 1 then the gain is likely to be very low and thus
    // unrealistic for modelling. This case is ignored.

    // Subtract the CDF to the upper bounds from the CDF of the lower bound
    // to get the discrete PMF

    // // Pad the CDF to avoid index-out-of bounds during interpolation
    final int padSize = 0;
    // padSize = (int) Math.ceil(1 / step) + 2;
    // list.resetQuick();
    // list.add(new double[padSize]);
    // list.add(g);
    // for (int i = padSize; i-- > 0;)
    // list.add(1);
    // double[] pd = list.toArray();
    final double[] pd = pg;

    list.resetQuick();

    final double[] x = SimpleArrayUtils.newArray(pd.length, zero - padSize * step, step);

    // Q. If the EM-CCD the distribution may have a Dirac delta at c=0 which
    // could break interpolation using a spline?
    final UnivariateInterpolator in =
        // (settings.getMode() == MODE_EM_CCD) ? new LinearInterpolator() :
        new SplineInterpolator();
    final UnivariateFunction f = in.interpolate(x, pd);

    int bound = (int) Math.floor(zero);

    list.add(0);
    zero--;
    double upperSum = 0;
    final double min = x[0];
    final double max = x[x.length - 1];
    while (upperSum < sum) {
      bound++;

      // Find the point at which the CDF should be computed
      lowerSum = upperSum;
      final double point = bound + offset;
      if (point < min) {
        upperSum = 0;
      } else if (point > max) {
        upperSum = sum;
      } else {
        upperSum = f.value(point);
      }
      list.add(upperSum - lowerSum);
    }
    list.add(0);
    return zero;
  }

  private static LikelihoodFunction getLikelihoodFunction(CameraModelAnalysisSettings settings) {
    final double alpha = 1.0 / getGain(settings);
    final double noise = getReadNoise(settings);
    final Model model = Model.forNumber(settings.getModel());
    switch (model) {
      case POISSON_PMF:
        return new PoissonFunction(alpha);
      case POISSON_DISRECTE:
        return new InterpolatedPoissonFunction(alpha, false);
      case POISSON_CONTINUOUS:
        return new InterpolatedPoissonFunction(alpha, true);
      case POISSON_GAUSSIAN_PDF:
      case POISSON_GAUSSIAN_PMF:
        final PoissonGaussianConvolutionFunction f1 =
            PoissonGaussianConvolutionFunction.createWithStandardDeviation(alpha, noise);
        f1.setComputePmf(model == Model.POISSON_GAUSSIAN_PMF);
        return f1;
      case POISSON_GAUSSIAN_APPROX:
        return PoissonGaussianFunction2.createWithStandardDeviation(alpha, noise);
      case POISSON_POISSON:
        return PoissonPoissonFunction.createWithStandardDeviation(alpha, noise);

      case POISSON_GAMMA_GAUSSIAN_PDF_CONVOLUTION:
        return PoissonGammaGaussianConvolutionFunction.createWithStandardDeviation(alpha, noise);

      case POISSON_GAMMA_PMF:
        return PoissonGammaFunction.createWithAlpha(alpha);

      case POISSON_GAMMA_GAUSSIAN_APPROX:
      case POISSON_GAMMA_GAUSSIAN_PDF_INTEGRATION:
      case POISSON_GAMMA_GAUSSIAN_PMF_INTEGRATION:
      case POISSON_GAMMA_GAUSSIAN_SIMPSON_INTEGRATION:
      case POISSON_GAMMA_GAUSSIAN_LEGENDRE_GAUSS_INTEGRATION:
        final PoissonGammaGaussianFunction f2 = new PoissonGammaGaussianFunction(alpha, noise);
        f2.setMinimumProbability(0);
        f2.setConvolutionMode(getConvolutionMode(model));
        // The function should return a PMF/PDF depending on how it is used
        f2.setPmfMode(!settings.getSimpsonIntegration());

        return f2;

      default:
        throw new IllegalStateException();
    }
  }

  private static boolean isPoissonGammaLikelihoodFunction(CameraModelAnalysisSettings settings) {
    final Model model = Model.forNumber(settings.getModel());
    switch (model) {
      case POISSON_GAMMA_GAUSSIAN_PDF_CONVOLUTION:
      case POISSON_GAMMA_PMF:
      case POISSON_GAMMA_GAUSSIAN_APPROX:
      case POISSON_GAMMA_GAUSSIAN_PDF_INTEGRATION:
      case POISSON_GAMMA_GAUSSIAN_PMF_INTEGRATION:
      case POISSON_GAMMA_GAUSSIAN_SIMPSON_INTEGRATION:
      case POISSON_GAMMA_GAUSSIAN_LEGENDRE_GAUSS_INTEGRATION:
        return true;

      default:
        return false;
    }
  }

  private static ConvolutionMode getConvolutionMode(Model model) {
    switch (model) {
      case POISSON_GAMMA_GAUSSIAN_APPROX:
        return ConvolutionMode.APPROXIMATION;
      case POISSON_GAMMA_GAUSSIAN_PDF_INTEGRATION:
        return ConvolutionMode.DISCRETE_PDF;
      case POISSON_GAMMA_GAUSSIAN_PMF_INTEGRATION:
        return ConvolutionMode.DISCRETE_PMF;
      case POISSON_GAMMA_GAUSSIAN_SIMPSON_INTEGRATION:
        return ConvolutionMode.SIMPSON_PDF;
      case POISSON_GAMMA_GAUSSIAN_LEGENDRE_GAUSS_INTEGRATION:
        return ConvolutionMode.LEGENDRE_GAUSS_PDF;

      default:
        throw new IllegalStateException();
    }
  }

  private static double getGain(CameraModelAnalysisSettings settings) {
    switch (settings.getMode()) {
      case MODE_CCD:
        return settings.getGain();
      case MODE_EM_CCD:
        return settings.getEmGain();
      case MODE_SCMOS:
        return settings.getCmosGain();
      default:
        throw new IllegalStateException();
    }
  }

  private static double getReadNoise(CameraModelAnalysisSettings settings) {
    switch (settings.getMode()) {
      case MODE_CCD:
        return settings.getNoise();
      case MODE_EM_CCD:
        return settings.getEmNoise();
      case MODE_SCMOS:
        return settings.getCmosNoise();
      default:
        throw new IllegalStateException();
    }
  }

  private static double[][] cumulativeHistogram(IntHistogram histogram) {
    final int[] h = histogram.histogramCounts;
    final double[] x = new double[h.length];
    final double[] y = new double[x.length];
    double sum = 0;
    for (int i = 0; i < x.length; i++) {
      // The cumulative histogram represents the probability of all values up to this one.
      // However the histogram is discrete so this is the probability of all values up to
      // but not including the next one.

      x[i] = histogram.getValue(i);
      sum += h[i];
      y[i] = sum;
    }
    SimpleArrayUtils.multiply(y, 1.0 / sum);
    y[y.length - 1] = 1; // Ensure total is 1
    return new double[][] {x, y};
  }

  private static class CachingUnivariateFunction implements UnivariateFunction {
    final LikelihoodFunction fun;
    final double theta;
    final TDoubleArrayList list = new TDoubleArrayList();

    public CachingUnivariateFunction(LikelihoodFunction fun, double theta) {
      this.fun = fun;
      this.theta = theta;
    }

    @Override
    public double value(double x) {
      final double v = fun.likelihood(x, theta);
      list.add(x);
      list.add(v);
      return v;
    }

    public void reset() {
      list.resetQuick();
    }
  }

  private static double[][] cumulativeDistribution(CameraModelAnalysisSettings settings,
      double[][] cdf, final LikelihoodFunction fun) {
    // Q. How to match this is the discrete cumulative histogram using the continuous
    // likelihood function:
    // 1. Compute integral up to the value
    // 2. Compute integral up to but not including the next value using trapezoid integration
    // 3. Compute integral up to but not including the next value using flat-top integration
    // Since the function will be used on continuous float data when fitting PSFs the best
    // match for how it will perform in practice is a continuous (trapezoid) integration.
    // The simplest is a flat-top integration.

    // Compute the probability at each value
    final double e = settings.getPhotons();
    double[] x = cdf[0];
    double[] y = new double[x.length];
    for (int i = 0; i < x.length; i++) {
      y[i] = fun.likelihood(x[i], e);
    }
    // Add more until the probability change is marginal
    double sum = MathUtils.sum(y);
    final TDoubleArrayList list = new TDoubleArrayList(y);
    for (int o = (int) x[x.length - 1] + 1;; o++) {
      final double p = fun.likelihood(o, e);
      list.add(p);
      if (p == 0) {
        break;
      }
      sum += p;
      if (p / sum < PROBABILITY_DELTA) {
        break;
      }
    }
    final TDoubleArrayList list2 = new TDoubleArrayList(10);
    for (int o = (int) x[0] - 1;; o--) {
      final double p = fun.likelihood(o, e);
      list2.add(p);
      if (p == 0) {
        break;
      }
      sum += p;
      if (p / sum < PROBABILITY_DELTA) {
        break;
      }
    }
    // Insert at start
    double start = x[0];
    if (!list2.isEmpty()) {
      start -= list2.size();
      list2.reverse();
      list.insert(0, list2.toArray());
    }

    y = list.toArray();
    x = SimpleArrayUtils.newArray(y.length, start, 1.0);

    if (settings.getSimpsonIntegration()) {
      int c0 = -1;
      double dirac = 0;
      int minc = 0;
      int maxc = 0;
      CachingUnivariateFunction uf = null;

      if (settings.getMode() == MODE_EM_CCD && isPoissonGammaLikelihoodFunction(settings)) {
        // A spike is expected at c=0 due to the Dirac delta contribution.
        // This breaks integration, especially when noise < 0.1.
        // Fix by integrating around c=0 fully then integrating the rest.
        c0 = Arrays.binarySearch(x, 0);
        final double noise = getReadNoise(settings);
        final double p = settings.getPhotons();
        if (noise == 0) {
          // Pure Poisson-Gamma. Just subtract the delta, do the simple integration
          // below and add the delta back. Only functions that support noise==0
          // will be allowed so this solution works.
          dirac = PoissonGammaFunction.dirac(p);
          if (c0 != -1) {
            y[c0] -= dirac;
          }
        } else {
          // Fix integration around c=0 using the range of the Gaussian
          minc = (int) Math.max(x[0], Math.floor(-5 * noise));
          maxc = (int) Math.min(x[x.length - 1], Math.ceil(5 * noise));
          uf = new CachingUnivariateFunction(fun, p);
        }
      }

      // Use Simpson's integration with n=4 to get the integral of the probability
      // over the range of each count.

      // Note the Poisson-Gamma function cannot be integrated with the
      // Dirac delta function at c==0

      // Compute the extra function points
      final double[] f = new double[y.length * SIMPSON_N + 1];
      int index = f.length;

      final int mod;
      if (settings.getRoundDown()) {
        // Do this final value outside the loop as y[index/n] does not exist
        mod = 0;
        index--;
        f[index] = fun.likelihood(start + index * SIMPSON_H, e);
      } else {
        // Used when computing for rounding up/down
        mod = SIMPSON_N_2;
        start -= SIMPSON_N_2 * SIMPSON_H;
      }

      while (index-- > 0) {
        if (index % SIMPSON_N == mod) {
          f[index] = y[index / SIMPSON_N];
        } else {
          f[index] = fun.likelihood(start + index * SIMPSON_H, e);
        }
      }

      // Compute indices for the integral
      final TIntArrayList tmp = new TIntArrayList(SIMPSON_N);
      for (int j = 1; j <= SIMPSON_N_2 - 1; j++) {
        tmp.add(2 * j);
      }
      final int[] i2 = tmp.toArray();
      tmp.resetQuick();
      for (int j = 1; j <= SIMPSON_N_2; j++) {
        tmp.add(2 * j - 1);
      }
      final int[] i4 = tmp.toArray();

      // Compute integral
      for (int i = 0; i < y.length; i++) {
        final int a = i * SIMPSON_N;
        final int b = a + SIMPSON_N;
        sum = f[a] + f[b];
        for (int j = i2.length; j-- > 0;) {
          sum += 2 * f[a + i2[j]];
        }
        for (int j = i4.length; j-- > 0;) {
          sum += 4 * f[a + i4[j]];
        }
        sum *= SIMPSON_H / 3;
        // System.out.printf("y[%d] = %f => %f\n", i, y[i], s);
        y[i] = sum;
      }

      // Fix Poisson-Gamma ...
      if (c0 != -1) {
        if (uf != null) {
          // Convolved Poisson-Gamma. Fix in the range of the Gaussian around c=0
          final CustomSimpsonIntegrator in = new CustomSimpsonIntegrator(SIMPSON_RELATIVE_ACCURACY,
              SIMPSON_ABSOLUTE_ACCURACY, SIMPSON_MIN_ITERATION_COUNT,
              CustomSimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT);
          final double lower = (settings.getRoundDown()) ? 0 : -0.5;
          final double upper = lower + 1;
          // Switch from c<=maxc to c<maxc. Avoid double computation at minc==maxc
          if (maxc != minc) {
            maxc++;
          }
          maxc++;
          for (int c = minc, i = Arrays.binarySearch(x, minc); c < maxc; c++, i++) {
            uf.reset();
            try {
              y[i] = in.integrate(2000, uf, c + lower, c + upper);
            } catch (final TooManyEvaluationsException ex) {
              // System.out.printf("Integration failed: c=%g-%g\n", c + lower, c + upper);
              // Q. Is the last sum valid?
              if (in.getLastSum() > 0) {
                y[i] = in.getLastSum();
              } else {
                // Otherwise use all the cached values to compute a sum
                // using the trapezoid rule. This will underestimate the sum.

                // Note: The Simpson integrator will have computed the edge values
                // as the first two values in the cache.
                final double[] g = uf.list.toArray();
                final double dx = (g[3] - g[1]) / in.getN();
                final int total = 1 + 2 * ((int) in.getN());
                sum = 0;
                for (int j = 4; j < total; j += 2) {
                  sum += g[j];
                }
                y[i] = (g[0] + g[2] + 2 * sum) / dx;
              }
            }
          }
        } else {
          // Pure Poisson-Gamma. Just add back the delta.
          y[c0] += dirac;
        }
      }
    }

    // Simple flat-top integration
    sum = 0;
    for (int i = 0; i < y.length; i++) {
      sum += y[i];
      y[i] = sum;
    }

    // Check if sum is approximately 1
    // System.out.printf("gain=%g, sum=%g\n", getGain(settings), sum);

    return new double[][] {x, y};
  }

  /**
   * Compute the Kolmogorov distance as the supremum (maximum) difference between the two cumulative
   * probability distributions. https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test
   *
   * <p>Also compute the mean distance between the two CDFs over the range of CDF 1.
   *
   * @param cdf1 the cdf 1
   * @param cdf2 the cdf 2
   * @return [distance,value,mean distance]
   */
  private static double[] getDistance(double[][] cdf1, double[][] cdf2) {
    // Find the offset
    int offset = 0;
    final double[] x1 = cdf1[0];
    final double[] x2 = cdf2[0];
    final double[] y1 = cdf1[1];
    final double[] y2 = cdf2[1];
    while (x2[offset] < x1[0]) {
      offset++;
    }
    double distance = 0;
    double value = x1[0];
    double area = 0;
    for (int i = 0; i < x1.length; i++) {
      final double d = Math.abs(y1[i] - y2[offset++]);
      if (distance < d) {
        distance = d;
        value = x1[i];
      }

      // Compute area:

      // This assumes both are discrete distributions
      area += d;

      // Note: This assumes both distributions are continuous between the values
      // and computes the actual area, including intersecting lines.
      // if (i != 0)
      // {
      // area += area(y1[i - 1], y1[i], y2[i], y2[i - 1]);
      // }
    }
    return new double[] {distance, value, area / x1.length};
  }

  private static final double[] areaX = {0, 1, 1, 0};

  @SuppressWarnings("unused")
  private static double area(double y1, double y2, double y3, double y4) {
    // Check if they cross
    if (!((y1 > y4 && y2 > y3) || (y1 < y4 && y2 < y3))) {
      final double[] intersection = new double[2];
      if (GeometryUtils.getIntersection(0, y1, 1, y2, 1, y3, 0, y4, intersection)) {
        // Compute area as two triangles
        return GeometryUtils.getArea(0, y1, 0, y4, intersection[0], intersection[1])
            + GeometryUtils.getArea(1, y2, 1, y3, intersection[0], intersection[1]);
      }
    }

    return Math.abs(GeometryUtils.getArea(areaX, new double[] {y1, y2, y3, y4}));
  }
}
