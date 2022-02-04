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
import ij.ImageStack;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.ImageProcessor;
import ij.text.TextWindow;
import ij.util.Tools;
import java.awt.AWTEvent;
import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.atomic.AtomicReference;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.NoiseEstimator;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.NoiseEstimatorMethod;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtosHelper;
import uk.ac.sussex.gdsc.smlm.ij.IJImageSource;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.ij.utils.ImageJImageConverter;
import uk.ac.sussex.gdsc.smlm.model.camera.CameraModel;
import uk.ac.sussex.gdsc.smlm.model.camera.CcdCameraModel;
import uk.ac.sussex.gdsc.smlm.model.camera.EmCcdCameraModel;
import uk.ac.sussex.gdsc.smlm.model.camera.NullCameraModel;

/**
 * Contains methods to find the noise in the provided image data.
 */
public class Noise implements ExtendedPlugInFilter, DialogListener {
  private static final String TITLE = "Noise Estimator";
  private static final int FLAGS =
      DOES_8G | DOES_16 | DOES_32 | PARALLELIZE_STACKS | FINAL_PROCESSING | NO_CHANGES;
  private static final String Y_AXIS_COUNT = "Noise (counts)";
  private static final String Y_AXIS_PHOTON = "Noise (photons)";

  private List<double[]> results;
  private PlugInFilterRunner pfr;
  private ImagePlus imp;
  private CalibrationWriter calibration;
  private CameraModel cameraModel;

  private String yAxisTitle;

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    int algorithm;
    int algorithm2;
    int lowestPixelsRange;

    Settings() {
      // Set defaults
      algorithm = NoiseEstimatorMethod.ALL_PIXELS_VALUE;
      algorithm2 = NoiseEstimatorMethod.QUICK_RESIDUALS_LEAST_TRIMMED_OF_SQUARES_VALUE;
      lowestPixelsRange = 6;
    }

    Settings(Settings source) {
      algorithm = source.algorithm;
      algorithm2 = source.algorithm2;
      lowestPixelsRange = source.lowestPixelsRange;
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
    if (arg.equalsIgnoreCase("final")) {
      showResults();
      return DONE;
    }
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (imp == null) {
      IJ.noImage();
      return DONE;
    }
    this.imp = imp;
    results = Collections.synchronizedList(new ArrayList<double[]>(imp.getStackSize()));
    return FLAGS;
  }

  @Override
  public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr) {
    // Select a camera model
    calibration = CalibrationWriter.create(SettingsManager.readCalibration(0));
    ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addMessage("Preprocess camera image and estimate global noise");
    PeakFit.addCameraOptions(gd, calibration);
    gd.showDialog();
    if (gd.wasCanceled()) {
      return DONE;
    }
    calibration.setCameraType(SettingsManager.getCameraTypeValues()[gd.getNextChoiceIndex()]);
    gd.collectOptions();
    SettingsManager.writeSettings(calibration.getCalibration());

    try {
      createCameraModel();
    } catch (final Exception ex) {
      IJ.error(TITLE, ex.getMessage());
      return DONE;
    }

    settings = Settings.load();

    // Do this now to prevent doing it in the preview.
    settings.save();

    // If using a stack, provide a preview graph of the noise for two methods
    if (imp.getStackSize() > 1) {
      this.pfr = pfr;

      drawPlot();

      gd = new ExtendedGenericDialog(TITLE);
      gd.addHelp(HelpUrls.getUrl("noise-estimator"));

      final String[] methodNames = SettingsManager.getNoiseEstimatorMethodNames();

      gd.addChoice("Method1 (blue)", methodNames, methodNames[settings.algorithm]);
      gd.addChoice("Method2 (red)", methodNames, methodNames[settings.algorithm2]);
      gd.addSlider("Lowest_radius", 1, 15, settings.lowestPixelsRange);

      gd.addDialogListener(this);
      gd.addMessage("Click OK to compute noise table using all methods");
      gd.showDialog();

      if (gd.wasCanceled() || !dialogItemChanged(gd, null)) {
        return DONE;
      }
    }

    return IJ.setupDialog(imp, FLAGS);
  }

  private void createCameraModel() {
    yAxisTitle = Y_AXIS_PHOTON;
    switch (calibration.getCameraType()) {
      case CCD:
        cameraModel = new CcdCameraModel(calibration.getBias(), calibration.getCountPerPhoton());
        break;
      case EMCCD:
        cameraModel = new EmCcdCameraModel(calibration.getBias(), calibration.getCountPerPhoton());
        break;
      case SCMOS:
        cameraModel = CameraModelManager.load(calibration.getCameraModelName());
        if (cameraModel == null) {
          throw new IllegalStateException(
              "No camera model for camera type: " + calibration.getCameraType());
        }
        cameraModel =
            PeakFit.cropCameraModel(cameraModel, IJImageSource.getBounds(imp), null, false);
        // Store for next time
        final Rectangle bounds = cameraModel.getBounds();
        final int ox = bounds.x;
        final int oy = bounds.y;
        // Reset origin for filtering
        if (ox != 0 || oy != 0) {
          cameraModel = cameraModel.copy();
          cameraModel.setOrigin(0, 0);
        }
        break;
      case CAMERA_TYPE_NA:
      case UNRECOGNIZED:
      default:
        cameraModel = new NullCameraModel();
        yAxisTitle = Y_AXIS_COUNT;
        break;
    }
  }

  @Override
  public boolean dialogItemChanged(GenericDialog gd, AWTEvent event) {
    settings.algorithm = gd.getNextChoiceIndex();
    settings.algorithm2 = gd.getNextChoiceIndex();
    settings.lowestPixelsRange = (int) gd.getNextNumber();
    if (gd.invalidNumber() || settings.lowestPixelsRange < 1) {
      return false;
    }
    if (gd.isShowing()) {
      drawPlot();
    }
    return true;
  }

  /**
   * Build a plot of the noise estimate from the current frame. Limit the preview to 100 frames.
   */
  private void drawPlot() {
    final NoiseEstimatorMethod[] values = SettingsManager.getNoiseEstimatorMethodValues();
    final NoiseEstimator.Method method1 =
        FitProtosHelper.convertNoiseEstimatorMethod(values[settings.algorithm]);
    final NoiseEstimator.Method method2 =
        FitProtosHelper.convertNoiseEstimatorMethod(values[settings.algorithm2]);
    IJ.showStatus("Estimating noise ...");

    final boolean twoMethods = method1 != method2;
    final boolean preserveResiduals =
        method1.name().contains("Residuals") && method2.name().contains("Residuals") && twoMethods;

    final int current = imp.getCurrentSlice();
    final int stackSize = imp.getStackSize();
    final int preview = 100;
    int start = current;
    int end = current + preview;
    if (end > stackSize) {
      final int shift = end - stackSize;
      start -= shift;
      end = stackSize;
      start = Math.max(1, start);
    }

    final int size = end - start + 1;
    final double[] xValues = new double[size];
    final double[] yValues1 = new double[size];
    final double[] yValues2 = (twoMethods) ? new double[size] : null;

    final ImageStack stack = imp.getImageStack();
    final Rectangle bounds = imp.getProcessor().getRoi();
    float[] buffer = null;
    for (int slice = start, i = 0; slice <= end; slice++, i++) {
      IJ.showProgress(i, size);
      final ImageProcessor ip = stack.getProcessor(slice);
      buffer = ImageJImageConverter.getData(ip.getPixels(), ip.getWidth(), ip.getHeight(), bounds,
          buffer);
      cameraModel.removeBiasAndGain(bounds, buffer);
      final NoiseEstimator ne = NoiseEstimator.wrap(buffer, bounds.width, bounds.height);
      ne.setPreserveResiduals(preserveResiduals);
      ne.setRange(settings.lowestPixelsRange);
      xValues[i] = slice;
      yValues1[i] = ne.getNoise(method1);
      if (yValues2 != null) {
        yValues2[i] = ne.getNoise(method2);
      }
    }
    IJ.showProgress(1);

    IJ.showStatus("Plotting noise ...");

    // Get limits
    final double[] a = Tools.getMinMax(xValues);
    final double[] b1 = Tools.getMinMax(yValues1);
    if (twoMethods) {
      final double[] b2 = Tools.getMinMax(yValues2);
      b1[0] = Math.min(b1[0], b2[0]);
      b1[1] = Math.max(b1[1], b2[1]);
    }

    final String title = imp.getTitle() + " Noise";
    final Plot plot = new Plot(title, "Slice", yAxisTitle);
    double range = b1[1] - b1[0];
    if (range == 0) {
      range = 1;
    }
    plot.setLimits(a[0], a[1], b1[0] - 0.05 * range, b1[1] + 0.05 * range);
    plot.setColor(Color.blue);
    plot.addPoints(xValues, yValues1, Plot.LINE);
    String label = String.format("%s (Blue) = %s", trim(method1.getName()),
        MathUtils.rounded(Statistics.create(yValues1).getMean()));
    if (twoMethods) {
      plot.setColor(Color.red);
      plot.addPoints(xValues, yValues2, Plot.LINE);
      label += String.format(", %s (Red) = %s", trim(method2.getName()),
          MathUtils.rounded(Statistics.create(yValues2).getMean()));
    }
    plot.setColor(Color.black);
    plot.addLabel(0, 0, label);

    ImageJUtils.display(title, plot);

    IJ.showStatus("");
  }

  private static Object trim(String name) {
    return (name.length() > 20) ? name.substring(0, 20) + "..." : name;
  }

  @Override
  public void run(ImageProcessor ip) {
    // Perform all methods and add to the results
    final double[] result = new double[NoiseEstimator.Method.values().length + 1];
    int index = 0;
    result[index++] = (pfr == null) ? 1 : pfr.getSliceNumber();
    final Rectangle bounds = ip.getRoi();
    final float[] buffer =
        ImageJImageConverter.getData(ip.getPixels(), ip.getWidth(), ip.getHeight(), bounds, null);
    cameraModel.removeBiasAndGain(bounds, buffer);
    final NoiseEstimator ne = NoiseEstimator.wrap(buffer, bounds.width, bounds.height);
    ne.setPreserveResiduals(true);
    for (final NoiseEstimator.Method m : NoiseEstimator.Method.values()) {
      result[index++] = ne.getNoise(m);
    }
    results.add(result);
  }

  @Override
  public void setNPasses(int passes) {
    // Do nothing
  }

  private void showResults() {
    // Sort on slice number
    Collections.sort(results, (o1, o2) -> Double.compare(o1[0], o2[0]));

    // Slow when there are lots of results ... Could change the output options in the future
    final TextWindow tw = new TextWindow(imp.getTitle() + " Noise", createHeader(), "", 800, 400);
    for (final double[] result : results) {
      tw.append(createResult(result));
    }
  }

  private static String createHeader() {
    final StringBuilder sb = new StringBuilder("Slice");
    for (final NoiseEstimator.Method m : NoiseEstimator.Method.values()) {
      sb.append('\t').append(m);
    }
    return sb.toString();
  }

  private static String createResult(double[] result) {
    final StringBuilder sb = new StringBuilder("" + (int) result[0]);
    for (int i = 1; i < result.length; i++) {
      sb.append('\t').append(result[i]);
    }
    return sb.toString();
  }
}
