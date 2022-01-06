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
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.awt.AWTEvent;
import java.util.concurrent.atomic.AtomicReference;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.logging.Ticker;
import uk.ac.sussex.gdsc.smlm.filters.FhtFilter;
import uk.ac.sussex.gdsc.smlm.filters.FhtFilter.Operation;
import uk.ac.sussex.gdsc.smlm.filters.KernelFilter;
import uk.ac.sussex.gdsc.smlm.filters.ZeroKernelFilter;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;

/**
 * Convolve an image with a kernel from another image.
 */
public class ImageKernelFilter implements ExtendedPlugInFilter, DialogListener {
  private static final String TITLE = "Image Kernel Filter";
  private static final int FLAGS = DOES_8G | DOES_16 | DOES_32 | KEEP_PREVIEW | PARALLELIZE_STACKS
      | CONVERT_TO_FLOAT | FINAL_PROCESSING;

  // Ensure not null
  private Ticker ticker = Ticker.getDefaultInstance();

  private int lastId;
  private int lastMethod = -1;
  private int lastFilter = -1;
  private boolean lastZero;
  private KernelFilter kf;
  private FhtFilter ff;
  private ImagePlus dataImp;
  private ImagePlus kernelImp;

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    private static final String[] METHODS = {"Spatial domain", "FHT"};
    private static final int METHOD_SPATIAL = 0;
    private static final int METHOD_FHT = 1;
    private static final String[] FILTERS;

    static {
      FILTERS = SettingsManager.getNames((Object[]) Operation.values());
    }

    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    String imageTitle;
    int method;
    int filter;
    int border;
    boolean zero;

    Settings() {
      // Set defaults
      imageTitle = "";
      method = METHOD_FHT;
      filter = Operation.CORRELATION.ordinal();
    }

    Settings(Settings source) {
      imageTitle = source.imageTitle;
      method = source.method;
      filter = source.filter;
      border = source.border;
      zero = source.zero;
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
    if ("final".equals(arg)) {
      imp.getProcessor().resetMinAndMax();
      imp.updateAndDraw();
      return DONE;
    }

    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (imp == null) {
      IJ.noImage();
      return DONE;
    }
    return FLAGS;
  }

  @Override
  public void run(ImageProcessor ip) {
    final float[] data = (float[]) ip.getPixels();
    final int w = ip.getWidth();
    final int h = ip.getHeight();
    if (settings.method == Settings.METHOD_SPATIAL) {
      kf.convolve(data, w, h, settings.border);
    } else {
      // Use a clone for thread safety
      final FhtFilter f = (ticker.getTotal() > 1) ? ff.copy() : ff;
      f.filter(data, w, h, settings.border);
    }
    if (ticker.getTotal() == 1) {
      ip.resetMinAndMax();
    }
    ticker.tick();
  }

  @Override
  public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr) {
    // Get available kernels
    final String[] names = ImageJUtils.getImageList(ImageJUtils.GREY_SCALE | ImageJUtils.SINGLE);
    if (names.length == 0) {
      IJ.error(TITLE, "No suitable kernel images");
      return DONE;
    }

    this.dataImp = imp;

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("image-kernel-filter"));

    gd.addMessage("Convolve an image using another image as the convolution kernel");

    settings = Settings.load();
    gd.addChoice("Kernel_image", names, settings.imageTitle);
    gd.addChoice("Method", Settings.METHODS, settings.method);
    gd.addChoice("Filter", Settings.FILTERS, settings.filter);
    gd.addSlider("Border", 0, 10, settings.border);
    gd.addCheckbox("Zero_outside_image", settings.zero);

    gd.addDialogListener(this);
    gd.addPreviewCheckbox(pfr);

    // Only need do this once. Do it here to save settings from the preview
    settings.save();

    gd.showDialog();

    if (gd.wasCanceled() || !dialogItemChanged(gd, null)) {
      return DONE;
    }

    return IJ.setupDialog(imp, FLAGS);
  }

  @Override
  public boolean dialogItemChanged(GenericDialog gd, AWTEvent event) {
    settings.imageTitle = gd.getNextChoice();
    settings.method = gd.getNextChoiceIndex();
    settings.filter = gd.getNextChoiceIndex();
    settings.border = (int) gd.getNextNumber();
    settings.zero = gd.getNextBoolean();

    kernelImp = WindowManager.getImage(settings.imageTitle);
    return (kernelImp != null);
  }

  @Override
  public void setNPasses(int passes) {
    // Create the kernel from the image
    boolean build = kernelImp.getID() != lastId || settings.method != lastMethod
        || settings.filter != lastFilter;
    build = build || (settings.method == Settings.METHOD_SPATIAL && kf == null);
    build = build || (settings.method == Settings.METHOD_FHT && ff == null);
    if (build) {
      final Operation operation = Operation.forOrdinal(settings.filter);
      FloatProcessor fp = kernelImp.getProcessor().toFloat(0, null);
      if (settings.method == Settings.METHOD_SPATIAL) {
        if (kf == null || kernelImp.getID() != lastId || settings.zero != lastZero) {
          fp = KernelFilter.pad(fp);
          final int kw = fp.getWidth();
          final int kh = fp.getHeight();
          final float[] kernel = (float[]) fp.getPixels();
          kf = (settings.zero) ? new ZeroKernelFilter(kernel, kw, kh)
              : new KernelFilter(kernel, kw, kh);
        }
        switch (operation) {
          case CONVOLUTION:
            kf.setConvolution(true);
            break;
          case CORRELATION:
            kf.setConvolution(false);
            break;
          case DECONVOLUTION:
          default:
            // Spatial filtering does not support anything other than convolution or correlation.
            ImageJUtils.log("Unsupported operation (%s), default to correlation",
                operation.getName());
            kf.setConvolution(false);
            break;
        }
      } else {
        if (ff == null || kernelImp.getID() != lastId) {
          final int kw = fp.getWidth();
          final int kh = fp.getHeight();
          final float[] kernel = (float[]) fp.getPixels();
          ff = new FhtFilter(kernel, kw, kh);
          ff.initialiseKernel(dataImp.getWidth(), dataImp.getHeight());
        }
        ff.setOperation(operation);
      }
      lastId = kernelImp.getID();
      lastMethod = settings.method;
      lastFilter = settings.filter;
      lastZero = settings.zero;
    }

    ticker = ImageJUtils.createTicker(passes, passes);
  }
}
