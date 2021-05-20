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

import gnu.trove.list.array.TIntArrayList;
import ij.IJ;
import ij.Prefs;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicReference;
import org.apache.commons.lang3.concurrent.ConcurrentRuntimeException;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot.HistogramPlotBuilder;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.StoredData;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
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
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    String inputOption;
    int border;

    Settings() {
      // Set defaults
      inputOption = "";
      border = 10;
    }

    Settings(Settings source) {
      inputOption = source.inputOption;
      border = source.border;
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
   * Store the density for a frame.
   */
  private static class FrameDensity {
    final int frame;
    final double density;

    /**
     * Create an instance.
     *
     * @param frame the frame
     * @param density the density
     */
    FrameDensity(int frame, double density) {
      this.frame = frame;
      this.density = density;
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

    // Scale to um^2 from px^2
    final double scale = Math.pow(results.getDistanceConverter(DistanceUnit.UM).convert(1), -2);

    final long start = System.currentTimeMillis();
    IJ.showStatus("Calculating density ...");


    final int border = settings.border;
    final int size = 2 * border + 1;
    final double minDensity = Math.pow(size, -2);
    ImageJUtils.log("%s : %s : Minimum density in %dx%d px = %s um^-2", TITLE, results.getName(),
        size, size, MathUtils.rounded(minDensity * scale));

    results.sort();
    final FrameCounter counter = results.newFrameCounter();
    final TIntArrayList x = new TIntArrayList();
    final TIntArrayList y = new TIntArrayList();
    final ExecutorService es = Executors.newFixedThreadPool(Prefs.getThreads());
    final LocalList<Future<FrameDensity>> futures = new LocalList<>();
    results.forEach((PeakResultProcedure) (peak) -> {
      if (counter.advance(peak.getFrame())) {
        final int frame = peak.getFrame();
        final int[] xx = x.toArray();
        final int[] yy = y.toArray();
        futures.add(es.submit(() -> {
          return new FrameDensity(frame, LocalDensity.estimate(xx, yy, border));
        }));
        x.resetQuick();
        y.resetQuick();
      }
      x.add((int) peak.getXPosition());
      y.add((int) peak.getYPosition());
    });
    futures.add(es.submit(() -> {
      return new FrameDensity(counter.currentFrame(),
          LocalDensity.estimate(x.toArray(), y.toArray(), border));
    }));
    es.shutdown();

    // Collect results
    final LocalList<FrameDensity> densities = new LocalList<>(futures.size());
    futures.stream().forEach(f -> {
      FrameDensity fd;
      try {
        fd = f.get();
        if (fd.density != 0) {
          densities.push(fd);
        }
      } catch (InterruptedException | ExecutionException ex) {
        if (ex instanceof InterruptedException) {
          // Restore interrupted state...
          Thread.currentThread().interrupt();
        }
        throw new ConcurrentRuntimeException(ex);
      }
    });

    densities.sort((o1, o2) -> Integer.compare(o1.frame, o2.frame));

    // Plot density
    final Statistics stats = new Statistics();
    final float[] frame = new float[densities.size()];
    final float[] density = new float[densities.size()];
    densities.stream().forEach(fd -> {
      final double d = fd.density * scale;
      frame[stats.getN()] = fd.frame;
      density[stats.getN()] = (float) d;
      stats.add(d);
    });
    final double mean = stats.getMean();
    final double sd = stats.getStandardDeviation();
    final Plot plot = new Plot("Frame vs Density", "Frame", "Density (um^-2)");
    plot.addPoints(frame, density, Plot.CIRCLE);
    plot.addLabel(0, 0,
        String.format("Density = %s +/- %s um^-2", MathUtils.rounded(mean), MathUtils.rounded(sd)));
    final WindowOrganiser wo = new WindowOrganiser();
    ImageJUtils.display(plot.getTitle(), plot, wo);

    // Histogram density
    new HistogramPlotBuilder("Local", StoredData.create(density), "Density (um^-2)")
        .setPlotLabel(String.format("Minimum density in %dx%d px = %s um^-2", size, size,
            MathUtils.rounded(minDensity * scale)))
        .show(wo);

    wo.tile();

    IJ.showStatus(
        TITLE + " complete : " + TextUtils.millisToString(System.currentTimeMillis() - start));
  }

  private boolean showDialog() {
    settings = Settings.load();

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("density-estimator"));

    gd.addMessage(TextUtils
        .wrap("Compute local density using NxN squares around localisations. The border size N "
            + "should reflect the extent of the PSF in pixels. The density is zero unless PSFs "
            + "are overlapping. All overlapping PSFs are joined to an area and the density of "
            + "localisations in each area is computed.", 80));
    ResultsManager.addInput(gd, settings.inputOption, InputSource.MEMORY);
    gd.addSlider("Border (N)", 0, 20, settings.border);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.inputOption = ResultsManager.getInputSource(gd);
    settings.border = (int) gd.getNextNumber();

    settings.save();

    return true;
  }
}
