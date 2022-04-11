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
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Roi;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.awt.AWTEvent;
import java.awt.Rectangle;
import java.util.concurrent.atomic.AtomicReference;
import uk.ac.sussex.gdsc.core.ij.ImageAdapter;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.DataFilterMethod;
import uk.ac.sussex.gdsc.smlm.engine.FitEngineConfiguration;
import uk.ac.sussex.gdsc.smlm.filters.DataProcessor;
import uk.ac.sussex.gdsc.smlm.filters.DifferenceSpotFilter;
import uk.ac.sussex.gdsc.smlm.filters.MaximaSpotFilter;
import uk.ac.sussex.gdsc.smlm.filters.SingleSpotFilter;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;

/**
 * Smooths the selected rectangular ROI using a mean filter.
 */
public class SmoothImage implements ExtendedPlugInFilter, DialogListener {
  private static final String TITLE = "Smooth Image";
  private static final DataFilterMethod[] FILTERS;
  private static final String[] FILTER_NAMES;

  static {
    FILTERS = SettingsManager.getDataFilterMethodValues();
    FILTER_NAMES = SettingsManager.getDataFilterMethodNames();
  }

  private static final int FLAGS =
      DOES_16 | DOES_8G | DOES_32 | PARALLELIZE_STACKS | FINAL_PROCESSING;

  // Track image contrast adjustments when the dialog is visible

  /** The image ID. */
  private int imageId;
  /** Flag to show when the dialog is visible. */
  private boolean dialogVisible;
  /** The display min of the image. */
  private double min;
  /** The display max of the image. */
  private double max;

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> INSTANCE = new AtomicReference<>(new Settings());

    int filter1;
    double smooth1;
    boolean differenceFilter;
    int filter2;
    double smooth2;
    boolean autoAdjustConstrast;
    boolean allowInversion;

    Settings() {
      smooth1 = 1;
      smooth2 = 3;
      autoAdjustConstrast = true;
    }

    Settings(Settings source) {
      filter1 = source.filter1;
      smooth1 = source.smooth1;
      differenceFilter = source.differenceFilter;
      filter2 = source.filter2;
      smooth2 = source.smooth2;
      autoAdjustConstrast = source.autoAdjustConstrast;
      allowInversion = source.allowInversion;
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

  @Override
  public int setup(String arg, ImagePlus imp) {
    if ("final".equals(arg)) {
      imp.updateAndDraw();

      return DONE;
    }

    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (imp == null) {
      IJ.noImage();
      return DONE;
    }

    final Roi roi = imp.getRoi();
    if (roi != null && roi.getType() != Roi.RECTANGLE) {
      IJ.error("Rectangular ROI required");
      return DONE;
    }

    // Support simultaneous contrast adjustment by saving the min/max
    // only when in the preview.
    imageId = imp.getID();
    saveMinMax(imp);
    ImagePlus.addImageListener(new ImageAdapter() {
      @Override
      public void imageUpdated(ImagePlus imp) {
        // Only when the dialog is open
        if (imp.getID() == imageId && dialogVisible) {
          saveMinMax(imp);
        }
      }
    });

    return FLAGS;
  }


  private void saveMinMax(ImagePlus imp) {
    final ImageProcessor ip = imp.getProcessor();
    min = ip.getMin();
    max = ip.getMax();
  }

  @Override
  public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr) {
    // Note: We cannot use a NonBlockingGenericDialog as scrolling through the image
    // throws away the snap shot. The pixel data for the previous slice is then fixed
    // with the preview. So we can only support a single slice.

    final NonBlockingGenericDialog gd = new NonBlockingGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("smooth-image"));

    settings = Settings.load();
    gd.addMessage("Smooth image:");
    gd.addChoice("Spot_filter", FILTER_NAMES, FILTER_NAMES[settings.filter1]);
    gd.addSlider("Smoothing", 0, 4.5, settings.smooth1);
    gd.addCheckbox("Difference_filter", settings.differenceFilter);
    gd.addChoice("Spot_filter2", FILTER_NAMES, FILTER_NAMES[settings.filter2]);
    gd.addSlider("Smoothing2", 1.5, 6, settings.smooth2);
    gd.addCheckbox("Auto_adjust_contrast", settings.autoAdjustConstrast);
    gd.addCheckbox("Allow_inversion", settings.allowInversion);

    gd.addPreviewCheckbox(pfr);
    gd.addDialogListener(this);
    dialogVisible = true;
    gd.showDialog();
    dialogVisible = false;

    if (gd.wasCanceled() || !dialogItemChanged(gd, null)) {
      return DONE;
    }

    return IJ.setupDialog(imp, FLAGS);
  }

  @Override
  public boolean dialogItemChanged(GenericDialog gd, AWTEvent event) {
    settings.filter1 = gd.getNextChoiceIndex();
    settings.smooth1 = gd.getNextNumber();
    settings.differenceFilter = gd.getNextBoolean();
    if (settings.differenceFilter) {
      settings.filter2 = gd.getNextChoiceIndex();
      settings.smooth2 = gd.getNextNumber();
    }
    settings.autoAdjustConstrast = gd.getNextBoolean();
    settings.allowInversion = gd.getNextBoolean();
    settings.save();
    return !gd.invalidNumber();
  }

  @Override
  public void run(ImageProcessor ip) {
    final Rectangle bounds = ip.getRoi();

    // Crop to the ROI
    FloatProcessor fp = ip.crop().toFloat(0, null);

    float[] data = (float[]) fp.getPixels();

    final MaximaSpotFilter filter = createSpotFilter();
    final int width = fp.getWidth();
    final int height = fp.getHeight();
    data = filter.preprocessData(data, width, height);

    fp = new FloatProcessor(width, height, data);
    ip.insert(fp, bounds.x, bounds.y);
    if (settings.autoAdjustConstrast) {
      ip.setMinAndMax(fp.getMin(), fp.getMax());
    } else {
      ip.setMinAndMax(min, max);
    }
  }

  private MaximaSpotFilter createSpotFilter() {
    final int search = 1;
    final int border = 0;
    final DataProcessor processor0 = FitEngineConfiguration.createDataProcessor(border,
        FILTERS[settings.filter1], settings.smooth1);
    if (settings.differenceFilter) {
      final DataProcessor processor1 = FitEngineConfiguration.createDataProcessor(border,
          FILTERS[settings.filter2], settings.smooth2);
      return new DifferenceSpotFilter(search, border, processor0, processor1,
          settings.allowInversion);
    }
    return new SingleSpotFilter(search, border, processor0);
  }

  @Override
  public void setNPasses(int passes) {
    // Nothing to do
  }
}
