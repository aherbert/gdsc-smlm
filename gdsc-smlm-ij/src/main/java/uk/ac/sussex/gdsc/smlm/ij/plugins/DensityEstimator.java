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

import ij.IJ;
import ij.Prefs;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import java.awt.Rectangle;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicReference;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot.HistogramPlotBuilder;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.StoredData;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.concurrent.ConcurrencyUtils;
import uk.ac.sussex.gdsc.core.utils.function.IntDoubleConsumer;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.results.LocalDensity;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.count.FrameCounter;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;

/**
 * Estimate the localisation density.
 */
public class DensityEstimator implements PlugIn {
  private static final String TITLE = "Density Estimator";

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> INSTANCE = new AtomicReference<>(new Settings());

    String inputOption;
    int border;
    boolean includeSingles;

    Settings() {
      // Set defaults
      inputOption = "";
      border = 10;
    }

    Settings(Settings source) {
      inputOption = source.inputOption;
      border = source.border;
      includeSingles = source.includeSingles;
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
   * Store the density for a frame.
   */
  private static class FrameDensity implements Runnable, IntDoubleConsumer {
    final int frame;
    int[] x;
    int[] y;
    final int border;
    final boolean includeSingles;
    int[] counts;
    double[] values;
    int size;
    int singles;

    /**
     * Create an instance.
     *
     * @param frame the frame
     * @param x the x
     * @param y the y
     * @param border the border
     * @param includeSingles the include singles
     */
    FrameDensity(int frame, int[] x, int[] y, int border, boolean includeSingles) {
      this.frame = frame;
      this.x = x;
      this.y = y;
      this.border = border;
      this.includeSingles = includeSingles;
    }

    @Override
    public void run() {
      counts = new int[x.length];
      values = new double[x.length];
      LocalDensity.estimate(x, y, border, this);
      // Free memory; prevent running again by forcing a NPE
      x = y = null;
      counts = Arrays.copyOf(counts, size);
      values = Arrays.copyOf(values, size);
    }

    @Override
    public void accept(int t, double value) {
      if (t == 1) {
        singles++;
        if (!includeSingles) {
          return;
        }
      }
      counts[size] = t;
      values[size] = value;
      size++;
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    // Require some fit results and selected regions
    if (MemoryPeakResults.countMemorySize() == 0) {
      IJ.error(TITLE, "There are no fitting results in memory");
      return;
    }

    if (!showDialog()) {
      return;
    }

    // Currently this only supports pixel distance units
    final MemoryPeakResults results =
        ResultsManager.loadInputResults(settings.inputOption, false, DistanceUnit.PIXEL, null);
    if (MemoryPeakResults.isEmpty(results)) {
      IJ.error(TITLE, "No results could be loaded");
      IJ.showStatus("");
      return;
    }

    final long start = System.currentTimeMillis();
    IJ.showStatus("Calculating density ...");

    // Scale to um^2 from px^2
    final double scale = Math.pow(results.getDistanceConverter(DistanceUnit.UM).convertBack(1), 2);

    results.sort();
    final FrameCounter counter = results.newFrameCounter();

    final double localisationsPerFrame =
        (double) results.size() / (results.getLastFrame() - counter.currentFrame() + 1);
    final Rectangle bounds = results.getBounds(true);
    final double globalDensity = localisationsPerFrame / bounds.width / bounds.height;

    final int border = settings.border;
    final boolean includeSingles = settings.includeSingles;
    final int size = 2 * border + 1;
    final double minDensity = Math.pow(size, -2);
    ImageJUtils.log("%s : %s : Global density %s. Minimum density in %dx%d px = %s um^-2", TITLE,
        results.getName(), MathUtils.rounded(globalDensity * scale), size, size,
        MathUtils.rounded(minDensity * scale));

    final IntArrayList x = new IntArrayList();
    final IntArrayList y = new IntArrayList();
    final ExecutorService es = Executors.newFixedThreadPool(Prefs.getThreads());
    final LocalList<FrameDensity> densities = new LocalList<>();
    final LocalList<Future<?>> futures = new LocalList<>();
    results.forEach((PeakResultProcedure) (peak) -> {
      if (counter.advance(peak.getFrame())) {
        final FrameDensity fd = new FrameDensity(peak.getFrame(), x.toIntArray(), y.toIntArray(),
            border, includeSingles);
        densities.add(fd);
        futures.add(es.submit(fd));
        x.clear();
        y.clear();
      }
      x.add((int) peak.getXPosition());
      y.add((int) peak.getYPosition());
    });
    densities.add(new FrameDensity(counter.currentFrame(), x.toIntArray(), y.toIntArray(), border,
        includeSingles));
    futures.add(es.submit(densities.get(densities.size() - 1)));
    es.shutdown();

    // Wait
    ConcurrencyUtils.waitForCompletionUnchecked(futures);

    densities.sort((o1, o2) -> Integer.compare(o1.frame, o2.frame));

    final int total = densities.stream().mapToInt(fd -> fd.counts.length).sum();

    // Plot density
    final Statistics stats = new Statistics();
    final float[] frame = new float[total];
    final float[] density = new float[total];
    densities.stream().forEach(fd -> {
      for (int i = 0; i < fd.counts.length; i++) {
        final double d = (fd.counts[i] / fd.values[i]) * scale;
        frame[stats.getN()] = fd.frame;
        density[stats.getN()] = (float) d;
        stats.add(d);
      }
    });
    final double mean = stats.getMean();
    final double sd = stats.getStandardDeviation();

    final String label =
        String.format("Density = %s +/- %s um^-2", MathUtils.rounded(mean), MathUtils.rounded(sd));

    final Plot plot = new Plot("Frame vs Density", "Frame", "Density (um^-2)");
    plot.addPoints(frame, density, Plot.CIRCLE);
    plot.addLabel(0, 0, label);
    final WindowOrganiser wo = new WindowOrganiser();
    ImageJUtils.display(plot.getTitle(), plot, wo);

    // Histogram density
    new HistogramPlotBuilder("Local", StoredData.create(density), "Density (um^-2)")
        .setPlotLabel(label).show(wo);

    wo.tile();

    // Log the number of singles
    final int singles = densities.stream().mapToInt(fd -> fd.singles).sum();
    ImageJUtils.log("Singles %d / %d (%s%%)", singles, results.size(),
        MathUtils.rounded(100.0 * singles / results.size()));

    IJ.showStatus(
        TITLE + " complete : " + TextUtils.millisToString(System.currentTimeMillis() - start));
  }

  private boolean showDialog() {
    settings = Settings.load();

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("density-estimator"));

    gd.addMessage(TextUtils
        .wrap("Compute local density using NxN squares around localisations. The border size "
            + "should reflect the extent of the PSF in pixels (N=2*border+1). All overlapping "
            + "PSFs are joined to an area and the density of localisations in each area is "
            + "computed.", 80));
    ResultsManager.addInput(gd, settings.inputOption, InputSource.MEMORY);
    gd.addSlider("Border", 0, 20, settings.border);
    gd.addCheckbox("Include_singles", settings.includeSingles);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.inputOption = ResultsManager.getInputSource(gd);
    settings.border = (int) gd.getNextNumber();
    settings.includeSingles = gd.getNextBoolean();

    settings.save();

    return true;
  }
}
