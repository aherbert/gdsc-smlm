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
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.process.LUT;
import ij.text.TextWindow;
import java.awt.Color;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import java.util.stream.IntStream;
import org.apache.commons.rng.sampling.distribution.NormalizedGaussianSampler;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper.LutColour;
import uk.ac.sussex.gdsc.core.logging.Ticker;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.concurrent.ConcurrencyUtils;
import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;
import uk.ac.sussex.gdsc.core.utils.rng.UniformRandomProviders;

/**
 * Move a set of molecules within a depth-of-field (DoF) to generate a distribution of the
 * probability of remaining within the DoF.
 *
 * <p>This empirical distribution can be used to fit correction factors {@code a} and {@code b} used
 * scale the depth of field {@code dz} to {@code dz_corr} for a given diffusion coefficient
 * {@code D}. The scaling factors depend on the DoF {@code dz}, the diffusion time step {@code dt},
 * and the gap size {@code g} allowed in tracks.
 *
 * <pre>
 * dz_corr = dz + a * sqrt(D) + b
 * </pre>
 * 
 * <p>The corrected depth of field can be used to compute the probability that a diffusing molecule
 * remains within the depth-of-field after a given time. This is based on the Spot-On model
 * described in the methods section of the paper:
 *
 * <p>Hansen, A.S., Woringer, M., Grimm, J.B., Lavis, L.D., Tjian, R., and Darzacq, X. (2018) Robust
 * model-based analysis of single-particle tracking experiments with Spot-On. eLife 7, e33125.
 * doi:10.7554/eLife.33125.
 */
public class DiffusionDepthOfField implements PlugIn {
  private static final String TITLE = "Diffusion Depth of Field";

  private static final AtomicReference<TextWindow> TABLE_REF = new AtomicReference<>();

  // private final WindowOrganiser windowOrganiser = new WindowOrganiser();

  /** The plugin settings. */
  private Settings pluginSettings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> INSTANCE = new AtomicReference<>(new Settings());

    double depthOfField;
    double exposureTime;
    int gap;

    int numberOfMolecules;
    int maxT;
    double minD;
    double maxD;
    int sampleD;

    Settings() {
      depthOfField = 750;
      exposureTime = 10;
      gap = 1;
      // Set simulation defaults from the SpotOn paper
      numberOfMolecules = 50000;
      maxT = 15;
      minD = 1;
      maxD = 12;
      sampleD = 13;
    }

    Settings(Settings source) {
      depthOfField = source.depthOfField;
      exposureTime = source.exposureTime;
      gap = source.gap;

      numberOfMolecules = source.numberOfMolecules;
      maxT = source.maxT;
      minD = source.minD;
      maxD = source.maxD;
      sampleD = source.sampleD;
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
     * Save the settings. This can be called only once as it saves via a reference.
     */
    void save() {
      INSTANCE.set(this);
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    pluginSettings = Settings.load();
    pluginSettings.save();

    if (!showDialog()) {
      return;
    }

    IJ.showStatus("Simulating tracks...");
    final double halfDz = pluginSettings.depthOfField / 2000;
    final double dt = pluginSettings.exposureTime / 1000;
    final double g = pluginSettings.gap;
    final int maxT = pluginSettings.maxT;

    // Fit correction coefficient a and b.
    // Compute: dz_corr = dz + a * sqrt(D) + b
    //
    // dz_corr is used to compute the probability of remaining within the depth-of-field
    // P(dt, dz_corr, D) and fit to the empirical distribution from the simulation.

    // Simulate tracks across the depth of field.
    // Each track is scaled using the diffusion coefficient and the molecule tested if
    // it remains in the depth-of-field.

    // Sample diffusion coefficients
    final double[] sampleD = IntStream.rangeClosed(1, pluginSettings.sampleD).mapToDouble(i -> {
      // interpolate D in [min, max] using [1, n] as [0, 1]
      final double p = (i - 1.0) / (pluginSettings.sampleD - 1);
      return (1 - p) * pluginSettings.minD + p * pluginSettings.maxD;
    }).toArray();
    // Create the Gaussian step for each diffusion coefficient
    final double[] step = Arrays.stream(sampleD).map(
        // Std.dev of Gaussian for the step size: sqrt(2D * dt)
        d -> Math.sqrt(2 * d * dt)).toArray();

    final int threadCount = Prefs.getThreads();
    final ExecutorService executor = Executors.newFixedThreadPool(threadCount);
    final List<Future<int[][]>> futures = new LinkedList<>();

    final int total = pluginSettings.numberOfMolecules;
    final Ticker ticker = ImageJUtils.createTicker(total, threadCount);

    // TODO: Add a MSD total for the steps and get the simulated D

    final AtomicInteger position = new AtomicInteger(total);
    for (int n = 0; n < threadCount; n++) {
      final NormalizedGaussianSampler sampler =
          SamplerUtils.createNormalizedGaussianSampler(UniformRandomProviders.create());
      futures.add(executor.submit(() -> {
        // Simulated track z-position from origin (unscaled)
        final double[] zn = new double[maxT];
        // Count of number of molecules inside the depth of field
        // for each diffusion coefficient and time frame
        final int[][] observed = new int[step.length][maxT];
        for (;;) {
          final int p = position.decrementAndGet();
          if (p < 0) {
            break;
          }
          // Simulate the track from the origin
          double z = 0;
          for (int i = 0; i < maxT; i++) {
            zn[i] = z += sampler.sample();
          }
          // Simulate across the depth of field
          z = (p / (total - 1.0)) * halfDz;
          // for each diffusion rate test: |z| < z/2
          for (int j = 0; j < step.length; j++) {
            final double s = step[j];
            final int[] obs = observed[j];
            // Last frame the molecule was inside the DoF
            int last = -1;
            for (int i = 0; i < maxT; i++) {
              if (Math.abs(z + zn[i] * s) < halfDz) {
                obs[i]++;
                last = i;
              } else if (i - last >= g) {
                // Allow tracks with the configured gap size; otherwise molecule has been lost
                break;
              }
            }
          }
          ticker.tick();
        }
        return observed;
      }));
    }

    // Finish processing data
    ConcurrencyUtils.waitForCompletionUncheckedT(futures);

    executor.shutdown();

    // Collate results
    final double[][] probability = new double[step.length][maxT];
    for (Future<int[][]> f : futures) {
      try {
        int[][] observed = f.get();
        for (int j = 0; j < observed.length; j++) {
          for (int i = 0; i < maxT; i++) {
            probability[j][i] += observed[j][i];
          }
        }
      } catch (InterruptedException | ExecutionException e) {
        throw new RuntimeException(e);
      }
    }
    for (int j = 0; j < probability.length; j++) {
      for (int i = 0; i < maxT; i++) {
        probability[j][i] /= pluginSettings.numberOfMolecules;
      }
    }

    ImageJUtils.finished(TITLE + " done");

    // Plot the observed distribution.
    // For each D add a line for all time steps
    final float[] time = SimpleArrayUtils.newArray(maxT, 1.0f, 1.0f);
    SimpleArrayUtils.apply(time, x -> (float) (x * dt));
    final String title = TITLE + " observed";
    final Plot plot = new Plot(title, "Time (seconds)", "Probability");
    LUT lut = LutHelper.createLut(LutColour.RED_BLUE);
    for (int i = 0; i < probability.length; i++) {
      plot.setColor(LutHelper.getColour(lut, i + 1, 1, probability.length));
      plot.addPoints(time, SimpleArrayUtils.toFloat(probability[i]), Plot.LINE);
    }
    plot.setColor(Color.black);
    plot.setLimits(dt * 0.5, dt * maxT + dt * 0.5, 0, 1);

    ImageJUtils.display(title, plot);

    // TODO:
    // Move simulation to a method ???
    // Fit the observed probability
    // Add a table of results
    // Plot the fitted p(remaining) for each diffusion coefficient (use points or change LUTs)
  }

  private boolean showDialog() {
    final GenericDialog gd = new GenericDialog(TITLE);

    gd.addNumericField("Depth_of_field", pluginSettings.depthOfField, 1, 6, "nm");
    gd.addNumericField("Exposure_time", pluginSettings.exposureTime, 1, 6, "ms");
    gd.addSlider("Gap", 1, 5, pluginSettings.gap);

    gd.addNumericField("Number_of_molecules", pluginSettings.numberOfMolecules);
    gd.addNumericField("Max_t", pluginSettings.maxT);
    gd.addNumericField("Min_D", pluginSettings.minD, 3, 6, "um^2/s");
    gd.addNumericField("Max_D", pluginSettings.maxD, 3, 6, "um^2/s");
    gd.addSlider("Sample_D", 2, 15, pluginSettings.sampleD);

    gd.addHelp(HelpUrls.getUrl("diffusion-depth-of-field"));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }

    pluginSettings.depthOfField = gd.getNextNumber();
    pluginSettings.exposureTime = gd.getNextNumber();
    pluginSettings.gap = (int) gd.getNextNumber();

    pluginSettings.numberOfMolecules = (int) gd.getNextNumber();
    pluginSettings.maxT = (int) gd.getNextNumber();
    pluginSettings.minD = gd.getNextNumber();
    pluginSettings.maxD = gd.getNextNumber();
    pluginSettings.sampleD = (int) gd.getNextNumber();

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Depth of Field", pluginSettings.depthOfField);
      ParameterUtils.isAboveZero("Exposure time", pluginSettings.exposureTime);
      ParameterUtils.isAboveZero("Gap", pluginSettings.gap);

      ParameterUtils.isAboveZero("Number of molecules", pluginSettings.numberOfMolecules);
      ParameterUtils.isAboveZero("Max T", pluginSettings.maxT);
      ParameterUtils.isAboveZero("Min D", pluginSettings.minD);
      ParameterUtils.isAboveZero("Max D", pluginSettings.maxD);
      ParameterUtils.isAbove("Sample D", pluginSettings.sampleD, 1);

      ParameterUtils.isEqualOrAbove("Max D", pluginSettings.maxD, pluginSettings.minD);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    if (gd.invalidNumber()) {
      return false;
    }

    return true;
  }
  //
  // private TextWindow createMsdTable(double baseMsd) {
  // return ImageJUtils.refresh(TABLE_REF, () -> {
  // return new TextWindow("MSD Analysis", createHeader(baseMsd), "", 800, 300);
  // });
  // }
  //
  // private String createHeader(double baseMsd) {
  // final double apparentD = baseMsd / 4;
  // final StringBuilder sb = new StringBuilder(256);
//    //@formatter:off
//    sb.append(settings.getDiffusionRate()).append('\t')
//      .append(myPrecision).append('\t')
//      .append(MathUtils.rounded(apparentD)).append('\t')
//      .append(MathUtils.rounded(1.0 / settings.getStepsPerSecond())).append('\t')
//      .append(myAggregateSteps).append('\t');
//    //@formatter:on
  // // Exposure time is the aggregated frame time
  // exposureTime = myAggregateSteps / settings.getStepsPerSecond();
  // sb.append(MathUtils.rounded(exposureTime)).append('\t');
  // prefix = sb.toString();
  // return "D (um^2/s)\tPrecision (nm)\tDsim (um^2/s)\tStep (s)\tResolution\t"
  // + "Frame (s)\tt (s)\tn\tN\tMSD (um^2)\tD (um^2/s)";
  // }
  //
  // private String addResult(int step, double sum, int count) {
  // final StringBuilder sb = new StringBuilder();
  // // Exposure time is the aggregated frame time
  // final double msd = (sum / count) / conversionFactor;
  // // Jump distance separation is the number of steps
  // final double t = step / settings.getStepsPerSecond();
//    //@formatter:off
//    sb.append(prefix)
//      .append(MathUtils.rounded(t)).append('\t')
//      .append(MathUtils.rounded(t / exposureTime)).append('\t')
//      .append(count).append('\t')
//      // Not rounded to preserve precision
//      .append(msd).append('\t')
//      .append(msd / (4 * t));
//    //@formatter:on
  // return sb.toString();
  // }
}
