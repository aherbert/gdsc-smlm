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
import ij.gui.GenericDialog;
import ij.plugin.ZProjector;
import ij.plugin.filter.GaussianBlur;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import java.util.concurrent.atomic.AtomicReference;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;

/**
 * Produces a background intensity image and a mask from a sample image.
 *
 * <p>The input image should be representative of the super-resolution imaging conditions and so
 * will produce suitable input for the Create Data plugin to create realistic images.
 */
public class ImageBackground implements PlugInFilter {
  private static final String TITLE = "Image Background";

  private static final int FLAGS = DOES_16 | DOES_8G | DOES_32 | NO_CHANGES;
  private ImagePlus imp;

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    float bias;
    double sigma;

    Settings() {
      // Set defaults
      bias = 500;
      sigma = 2;
    }

    Settings(Settings source) {
      bias = source.bias;
      sigma = source.sigma;
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

  @Override
  public int setup(String arg, ImagePlus imp) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (imp == null) {
      IJ.noImage();
      return DONE;
    }

    this.imp = imp;

    return showDialog();
  }

  private int showDialog() {
    final GenericDialog gd = new GenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("image-background"));

    gd.addMessage(
        "Creates a background and mask image from a sample input stack\nusing a median projection");

    settings = Settings.load();
    gd.addNumericField("Bias", settings.bias, 0);
    gd.addSlider("Blur", 0, 20, settings.sigma);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return DONE;
    }

    settings.bias = (float) gd.getNextNumber();
    settings.sigma = gd.getNextNumber();
    settings.save();

    // Check arguments
    try {
      ParameterUtils.isPositive("Bias", settings.bias);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return DONE;
    }

    return FLAGS;
  }

  @Override
  public void run(ImageProcessor ip) {
    final ImageProcessor median = getProjection();

    final ImageProcessor background = applyBlur(median);
    subtractBias(background);

    ImageJUtils.display("Background", background);

    // Q. Is there a better way to do the thresholding for foreground pixels.
    // Ideally we want to outline cell shapes.
    final ImageProcessor mask = median.convertToByte(true);
    mask.autoThreshold();

    ImageJUtils.display("Mask", mask);
  }

  private ImageProcessor getProjection() {
    // Get median intensity projection
    final ZProjector p = new ZProjector(imp);
    p.setMethod(ZProjector.MEDIAN_METHOD);
    p.doProjection();
    return p.getProjection().getProcessor();
  }

  private ImageProcessor applyBlur(ImageProcessor median) {
    ImageProcessor blur = median;
    if (settings.sigma > 0) {
      blur = median.duplicate();
      final GaussianBlur gb = new GaussianBlur();
      gb.blurGaussian(blur, settings.sigma, settings.sigma, 0.0002);
    }
    return blur;
  }

  private void subtractBias(ImageProcessor background) {
    final float[] data = (float[]) background.getPixels();
    for (int i = 0; i < data.length; i++) {
      data[i] = Math.max(0f, data[i] - settings.bias);
    }
    background.resetMinAndMax();
  }
}
