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

import ij.IJ;
import ij.ImagePlus;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.util.Tools;

import org.apache.commons.math3.util.FastMath;

import java.awt.AWTEvent;
import java.awt.Label;
import java.util.concurrent.atomic.AtomicReference;

/**
 * Filters pixels using the surrounding region.
 */
public class PixelFilter implements ExtendedPlugInFilter, DialogListener {
  private static final String TITLE = "Pixel Filter";
  private static final int FLAGS = DOES_8G | DOES_16 | DOES_32 | PARALLELIZE_STACKS;

  private PlugInFilterRunner pfr;
  private double[] cachedRollingSum;
  private double[] cachedRollingSumSq;
  private boolean preview;
  private Label label;

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    int radius = 1;
    double error = 3;

    Settings() {
      // Set defaults
      radius = 1;
      error = 3;
    }

    Settings(Settings source) {
      radius = source.radius;
      error = source.error;
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

  @Override
  public int setup(String arg, ImagePlus imp) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (imp == null) {
      IJ.noImage();
      return DONE;
    }
    return FLAGS;
  }

  @Override
  public void run(ImageProcessor ip) {
    // Compute rolling sums
    final FloatProcessor fp = ip.toFloat(0, null);
    final float[] data = (float[]) ip.toFloat(0, null).getPixels();
    double[] rollingSum = null;
    double[] rollingSumSq = null;
    if (preview && cachedRollingSum != null) {
      rollingSum = cachedRollingSum;
      rollingSumSq = cachedRollingSumSq;
    }
    if (rollingSum == null || rollingSumSq == null) {
      rollingSum = new double[ip.getPixelCount()];
      rollingSumSq = new double[rollingSum.length];
      calculateRollingSums(fp, rollingSum, rollingSumSq);
    }

    int count = 0;
    final int maxx = ip.getWidth();
    final int maxy = ip.getHeight();
    for (int y = 0, i = 0; y < maxy; y++) {
      for (int x = 0; x < maxx; x++, i++) {
        double sum = 0;
        double sumSquares = 0;

        int minU = x - settings.radius - 1;
        final int maxU = FastMath.min(x + settings.radius, maxx - 1);
        final int maxV = FastMath.min(y + settings.radius, maxy - 1);

        // Compute sum from rolling sum using:
        // sum(u,v) =
        // + s(u+N,v+N)
        // - s(u-N-1,v+N)
        // - s(u+N,v-N-1)
        // + s(u-N-1,v-N-1)
        // Note:
        // s(u,v) = 0 when either u,v < 0
        // s(u,v) = s(umax,v) when u>umax
        // s(u,v) = s(u,vmax) when v>vmax
        // s(u,v) = s(umax,vmax) when u>umax,v>vmax
        // Likewise for ss

        // + s(u+N-1,v+N-1)
        int index = maxV * maxx + maxU;
        sum += rollingSum[index];
        sumSquares += rollingSumSq[index];

        if (minU >= 0) {
          // - s(u-1,v+N-1)
          index = maxV * maxx + minU;
          sum -= rollingSum[index];
          sumSquares -= rollingSumSq[index];
        }
        int minV = y - settings.radius - 1;
        if (minV >= 0) {
          // - s(u+N-1,v-1)
          index = minV * maxx + maxU;
          sum -= rollingSum[index];
          sumSquares -= rollingSumSq[index];

          if (minU >= 0) {
            // + s(u-1,v-1)
            index = minV * maxx + minU;
            sum += rollingSum[index];
            sumSquares += rollingSumSq[index];
          }
        }

        // Reset to bounds to calculate the number of pixels
        if (minU < 0) {
          minU = -1;
        }
        if (minV < 0) {
          minV = -1;
        }

        final int n = (maxU - minU) * (maxV - minV);

        if (n < 2) {
          continue;
        }

        // Get the sum of squared differences
        final double residuals = sumSquares - (sum * sum) / n;

        if (residuals > 0.0) {
          final double stdDev = Math.sqrt(residuals / (n - 1.0));
          final double mean = sum / n;

          if (Math.abs(data[i] - mean) / stdDev > settings.error) {
            ip.setf(i, (float) mean);
            count++;
          }
        }
      }
    }
    if (preview) {
      cachedRollingSum = rollingSum;
      cachedRollingSumSq = rollingSumSq;
      label.setText("Replaced " + count);
    } else if (pfr != null && count > 0) {
      IJ.log(String.format("Slice %d : Replaced %d pixels", pfr.getSliceNumber(), count));
    }
  }

  private static void calculateRollingSums(FloatProcessor ip, double[] sum, double[] sumSq) {
    // Compute the rolling sum and sum of squares
    // s(u,v) = f(u,v) + s(u-1,v) + s(u,v-1) - s(u-1,v-1)
    // ss(u,v) = f(u,v) * f(u,v) + ss(u-1,v) + ss(u,v-1) - ss(u-1,v-1)
    // where s(u,v) = ss(u,v) = 0 when either u,v < 0

    final int maxx = ip.getWidth();
    final int maxy = ip.getHeight();
    final float[] originalData = (float[]) ip.getPixels();
    final double[] data = Tools.toDouble(originalData);

    // First row
    double cs = 0; // Column sum
    double css = 0; // Column sum-squares
    for (int i = 0; i < maxx; i++) {
      cs += data[i];
      css += data[i] * data[i];
      sum[i] = cs;
      sumSq[i] = css;
    }

    // Remaining rows:
    // sum = rolling sum of row + sum of row above
    for (int y = 1; y < maxy; y++) {
      int index = y * maxx;
      cs = 0;
      css = 0;

      // Remaining columns
      for (int x = 0; x < maxx; x++, index++) {
        cs += data[index];
        css += data[index] * data[index];

        sum[index] = sum[index - maxx] + cs;
        sumSq[index] = sumSq[index - maxx] + css;
      }
    }
  }

  @Override
  public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr) {
    this.pfr = pfr;
    preview = true;

    final GenericDialog gd = new GenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);

    gd.addMessage("Replace pixels with mean if they are N StdDevs from the mean");

    settings = Settings.load();
    settings.save();

    gd.addSlider("Radius", 1, 5, settings.radius);
    gd.addSlider("Error (SD units)", 2.5, 7, settings.error);

    gd.addPreviewCheckbox(pfr);
    gd.addDialogListener(this);

    gd.addMessage("");
    label = (Label) gd.getMessage();

    gd.showDialog();

    if (gd.wasCanceled() || !dialogItemChanged(gd, null)) {
      return DONE;
    }

    preview = false;
    cachedRollingSum = cachedRollingSumSq = null;
    label = null;
    return IJ.setupDialog(imp, FLAGS);
  }

  @Override
  public boolean dialogItemChanged(GenericDialog gd, AWTEvent event) {
    label.setText("");
    settings.radius = (int) gd.getNextNumber();
    settings.error = gd.getNextNumber();
    return (!gd.invalidNumber() && settings.radius >= 1 && settings.error >= 0);
  }

  @Override
  public void setNPasses(int passes) {
    // Ignore
  }
}
