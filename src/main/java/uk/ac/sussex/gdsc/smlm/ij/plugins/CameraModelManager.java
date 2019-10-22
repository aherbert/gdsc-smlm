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

import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.utils.ExtendedStatistics;
import uk.ac.sussex.gdsc.core.utils.FileUtils;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.TurboList;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraModelResource;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraModelSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.CameraModelManagerSettings;
import uk.ac.sussex.gdsc.smlm.ij.IJImageSource;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.model.camera.CameraModel;
import uk.ac.sussex.gdsc.smlm.model.camera.PerPixelCameraModel;
import uk.ac.sussex.gdsc.smlm.results.ImageSource;

import ij.CompositeImage;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import org.apache.commons.lang3.ArrayUtils;

import java.awt.Rectangle;
import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicReference;

/**
 * This plugin handle the save and load of per-pixel camera models.
 */
public class CameraModelManager implements PlugIn {
  private static final String TITLE = "Camera Model Manager";
  private static final String INFO_TAG = "Per-pixel camera model data";

  private static AtomicReference<String> directory = new AtomicReference<>("");
  private static AtomicReference<String> filename = new AtomicReference<>("");

  /** Cache camera models for speed. */
  private static Map<String, PerPixelCameraModel> cameraModels =
      Collections.synchronizedMap(new LinkedHashMap<>());

  //@formatter:off
  private static final String[] OPTIONS = {
      "Load a camera model",
      "Load from directory",
      "Print all model details",
      "View a camera model",
      "Delete a camera model",
      "Filter an image"
  };
  //@formatter:on

  private CameraModelManagerSettings.Builder pluginSettings;

  /**
   * Provide a lazy initialisation holder for the CameraModelSettings.
   */
  private static class CameraModelSettingsHolder {
    private static CameraModelSettings settings;

    static {
      settings = SettingsManager.readCameraModelSettings(SettingsManager.FLAG_SILENT);
    }

    /**
     * Gets the settings held in memory.
     *
     * @return the settings
     */
    static CameraModelSettings getSettings() {
      return settings;
    }

    /**
     * Sets the settings to memory and save them to file.
     *
     * @param settings the settings
     * @return true, if successful writing to file
     */
    static synchronized boolean setSettings(CameraModelSettings settings) {
      CameraModelSettingsHolder.settings = settings;
      return SettingsManager.writeSettings(settings);
    }
  }

  /**
   * Save the camera model. The model will be named in the resources using the filename without the
   * extension or leading path entries.
   *
   * @param cameraModel the camera model
   * @param filename the filename
   * @return true, if successful
   */
  public static boolean save(PerPixelCameraModel cameraModel, String filename) {
    if (cameraModel == null || TextUtils.isNullOrEmpty(filename)) {
      return false;
    }

    // Try to save to file
    final String name = getName(filename);

    final ImageStack stack = new ImageStack(cameraModel.getWidth(), cameraModel.getHeight());
    stack.addSlice("Bias", cameraModel.getBias());
    stack.addSlice("Gain", cameraModel.getGain());
    stack.addSlice("Variance", cameraModel.getVariance());
    final ImagePlus imp = new ImagePlus(name, stack);
    imp.setIgnoreGlobalCalibration(true);
    final Calibration cal = imp.getCalibration();
    cal.xOrigin = cameraModel.getXOrigin();
    cal.yOrigin = cameraModel.getYOrigin();
    imp.setProperty("Info", INFO_TAG);
    // Do this to allow the filename to be something other than .tif
    final boolean ok = new FileSaver(imp).saveAsTiffStack(filename);

    if (ok) {
      saveResource(cameraModel, filename, name);
    }

    return ok;
  }

  private static String getName(String filename) {
    return FileUtils.removeExtension(FileUtils.getName(filename));
  }

  private static void saveResource(PerPixelCameraModel cameraModel, String filename, String name) {
    final CameraModelResource.Builder resource = CameraModelResource.newBuilder();
    resource.setX(cameraModel.getXOrigin());
    resource.setY(cameraModel.getYOrigin());
    resource.setWidth(cameraModel.getWidth());
    resource.setHeight(cameraModel.getHeight());
    resource.setFilename(filename);

    final CameraModelSettings settings = CameraModelSettingsHolder.getSettings();
    CameraModelSettingsHolder
        .setSettings(settings.toBuilder().putCameraModelResources(name, resource.build()).build());

    // Cache this
    cameraModels.put(name, cameraModel);
  }

  /**
   * Load the camera model. Returns null if the named model does not exist. Writes to the ImageJ log
   * if a problems occurred loading the model.
   *
   * @param name the name
   * @return the per pixel camera model (or null)
   */
  public static PerPixelCameraModel load(String name) {
    PerPixelCameraModel model = cameraModels.get(name);
    if (model == null) {
      final CameraModelSettings settings = CameraModelSettingsHolder.getSettings();
      // Try and get the named resource
      final CameraModelResource resource = settings.getCameraModelResourcesMap().get(name);
      if (resource == null) {
        return null;
      }
      model = loadFromFile(name, resource.getFilename());

      // Cache this
      cameraModels.put(name, model);
    }
    return model;
  }

  private static PerPixelCameraModel loadFromFile(String name, String filename) {
    // Try and load the resource
    final ImagePlus imp = IJ.openImage(filename);
    IJ.showStatus(""); // Remove the status from the ij.io.ImageWriter class

    if (imp == null) {
      ImageJUtils.log("Failed to load camera model %s data from file: ", name, filename);
      return null;
    }
    // Check stack size
    final ImageStack stack = imp.getImageStack();
    if (stack.getSize() != 3) {
      ImageJUtils.log("Camera model %s requires 3 image stack from file: %s", name, filename);
      return null;
    }
    // Get the origin
    imp.setIgnoreGlobalCalibration(true);
    final Calibration cal = imp.getCalibration();
    final Rectangle bounds =
        new Rectangle((int) cal.xOrigin, (int) cal.yOrigin, stack.getWidth(), stack.getHeight());
    try {
      final float[] bias = (float[]) stack.getPixels(1);
      final float[] gain = (float[]) stack.getPixels(2);
      final float[] variance = (float[]) stack.getPixels(3);
      return PerPixelCameraModel.create(bounds, bias, gain, variance);
    } catch (final Exception ex) {
      ImageJUtils.log("Failed to load camera model %s from file: %s. %s", name, filename,
          ex.getMessage());
    }
    return null;
  }

  /**
   * List the camera models.
   *
   * @param includeNone Set to true to include an invalid none model string
   * @return the list
   */
  public static String[] listCameraModels(boolean includeNone) {
    final CameraModelSettings settings = CameraModelSettingsHolder.getSettings();
    final List<String> list = createList(includeNone);
    list.addAll(settings.getCameraModelResourcesMap().keySet());
    return list.toArray(new String[0]);
  }

  /**
   * List the camera models with the correct dimensions.
   *
   * @param includeNone Set to true to include an empty string
   * @param width the width
   * @param height the height
   * @return the list
   */
  public static String[] listCameraModels(boolean includeNone, int width, int height) {
    final CameraModelSettings settings = CameraModelSettingsHolder.getSettings();
    final List<String> list = createList(includeNone);
    for (final Map.Entry<String, CameraModelResource> entry : settings.getCameraModelResourcesMap()
        .entrySet()) {
      final CameraModelResource resource = entry.getValue();
      if (resource.getWidth() == width && resource.getHeight() == height) {
        list.add(entry.getKey());
      }
    }
    return list.toArray(new String[0]);
  }

  private static List<String> createList(boolean includeNone) {
    final List<String> list = new TurboList<>();
    if (includeNone) {
      list.add("[None]");
    }
    return list;
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    String[] options = OPTIONS;
    final CameraModelSettings settings = CameraModelSettingsHolder.getSettings();
    if (settings.getCameraModelResourcesCount() == 0) {
      options = Arrays.copyOf(OPTIONS, 2);
    }

    pluginSettings = readCameraModelManagerSettings();

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addChoice("Option", options, pluginSettings.getOption());
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    pluginSettings.setOption(gd.getNextChoiceIndex());

    switch (pluginSettings.getOption()) {
      case 5:
        runFilterImage();
        break;
      case 4:
        runDeleteCameraModel(settings);
        break;
      case 3:
        runViewCameraModel(settings);
        break;
      case 2:
        runPrintCameraModels(settings);
        break;
      case 1:
        runLoadFromDirectory();
        break;
      case 0:
      default:
        runLoadFromFile();
    }

    writeCameraModelManagerSettings(pluginSettings);
  }

  private static synchronized CameraModelManagerSettings.Builder readCameraModelManagerSettings() {
    return SettingsManager.readCameraModelManagerSettings(0).toBuilder();
  }

  private static synchronized void
      writeCameraModelManagerSettings(CameraModelManagerSettings.Builder pluginSettings) {
    SettingsManager.writeSettings(pluginSettings);
  }

  private void runFilterImage() {
    // Select an image
    GenericDialog gd = new GenericDialog(TITLE);
    final String[] list = ImageJUtils.getImageList(ImageJUtils.GREY_SCALE);
    gd.addChoice("Image", list, pluginSettings.getImage());
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    final String image = gd.getNextChoice();
    pluginSettings.setImage(image);
    final ImagePlus imp = WindowManager.getImage(image);
    if (imp == null) {
      IJ.error(TITLE, "Failed to find image: " + image);
      return;
    }

    // Select the model
    gd = new GenericDialog(TITLE);
    final String[] models = listCameraModels(false);
    gd.addChoice("Model", models, pluginSettings.getSelected());
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    final String name = gd.getNextChoice();
    pluginSettings.setSelected(name);

    CameraModel cameraModel = load(name);

    if (cameraModel == null) {
      IJ.error(TITLE, "Failed to find camera data for model: " + name);
      return;
    }

    // Crop the model if appropriate
    try {
      Rectangle bounds;
      try {
        bounds = IJImageSource.getBounds(imp);
      } catch (final IllegalArgumentException ex) {
        bounds = new Rectangle(pluginSettings.getOriginX(), pluginSettings.getOriginY(),
            imp.getWidth(), imp.getHeight());
      }
      cameraModel = PeakFit.cropCameraModel(cameraModel, bounds, null, false);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return;
    }
    final Rectangle bounds = cameraModel.getBounds();
    pluginSettings.setOriginX(bounds.x);
    pluginSettings.setOriginY(bounds.y);

    // Reset origin for fast filtering
    cameraModel = cameraModel.copy();
    cameraModel.setOrigin(0, 0);

    // Filter all the frames
    final ImageSource source = new IJImageSource(imp);
    if (!source.open()) {
      IJ.error(TITLE, "Cannot open image: " + image);
    }
    final ImageStack stack = new ImageStack(imp.getWidth(), imp.getHeight());
    for (float[] data = source.next(); data != null; data = source.next()) {
      cameraModel.removeBiasAndGain(data);
      stack.addSlice(null, data);
    }

    final ImagePlus imp2 = new ImagePlus(imp.getTitle() + " Filtered", stack);
    imp2.copyScale(imp);
    imp2.show();
  }

  private void runDeleteCameraModel(CameraModelSettings settings) {
    final GenericDialog gd = new GenericDialog(TITLE);
    final String[] models = listCameraModels(false);
    gd.addChoice("Model", models, pluginSettings.getSelected());
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    final String name = gd.getNextChoice();
    pluginSettings.setSelected(name);

    final CameraModelResource resource = settings.getCameraModelResourcesMap().get(name);
    if (resource == null) {
      IJ.error(TITLE, "Failed to find camera data for model: " + name);
      return;
    }

    CameraModelSettingsHolder
        .setSettings(settings.toBuilder().removeCameraModelResources(name).build());

    ImageJUtils.log("Deleted camera model: %s\n%s", name, resource);
  }

  private static void runLoadFromDirectory() {
    final ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE);
    egd.addMessage("Load camera models from a directory.");
    egd.addDirectoryField("Directory", directory.get());
    egd.showDialog();
    if (egd.wasCanceled()) {
      return;
    }

    final String dir = egd.getNextString();
    directory.set(dir);

    final File[] fileList = (new File(dir)).listFiles(File::isFile);
    if (!ArrayUtils.isEmpty(fileList)) {
      for (final File file : fileList) {
        loadFromFileAndSaveResource(file.getPath());
      }
    }
  }

  private static void runLoadFromFile() {
    final ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE);
    egd.addMessage("Load a camera model from file.");
    egd.addFilenameField("Filename", filename.get());
    egd.showDialog();
    if (egd.wasCanceled()) {
      return;
    }

    final String file = egd.getNextString();
    filename.set(file);

    loadFromFileAndSaveResource(file);
  }

  private static void loadFromFileAndSaveResource(String filename) {
    final String name = getName(filename);
    final PerPixelCameraModel model = loadFromFile(name, filename);

    if (model != null) {
      saveResource(model, filename, name);
    }
  }

  private void runViewCameraModel(CameraModelSettings settings) {
    final GenericDialog gd = new GenericDialog(TITLE);
    final String[] models = listCameraModels(false);
    gd.addChoice("Model", models, pluginSettings.getSelected());
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    final String name = gd.getNextChoice();
    pluginSettings.setSelected(name);

    // Try and get the named resource
    final CameraModelResource resource = settings.getCameraModelResourcesMap().get(name);
    if (resource == null) {
      IJ.error(TITLE, "Failed to find camera data for model: " + name);
      return;
    }
    // Try and load the resource.
    // Do not use loadFromFile as that validates the model data. We just want
    // to view the raw image.
    ImagePlus imp = IJ.openImage(resource.getFilename());
    IJ.showStatus(""); // Remove the status from the ij.io.ImageWriter class
    if (imp == null) {
      IJ.error(TITLE, "Failed to load camera data for model: " + name);
      return;
    }
    final ImageStack stack = imp.getImageStack();
    if (stack.getSize() != 3) {
      IJ.error(TITLE, "Failed to load camera data for model: " + name
          + ".\nExpecting stack of size 3 but it was " + stack.getSize());
      return;
    }
    ImageJUtils.log("Camera model: %s\n%s", name, resource);
    for (int n = 1; n <= stack.getSize(); n++) {
      logStats(stack.getSliceLabel(n), stack.getProcessor(n));
    }

    // Show normalised variance: var/g^2
    final float[] varG2 = new float[stack.getWidth() * stack.getHeight()];
    try {
      final float[] gain = (float[]) stack.getPixels(2);
      final float[] variance = (float[]) stack.getPixels(3);

      final ExtendedStatistics stats1 = new ExtendedStatistics();
      final ExtendedStatistics stats2 = new ExtendedStatistics();
      for (int i = 0; i < gain.length; i++) {
        final double v1 = variance[i] / MathUtils.pow2(gain[i]);
        varG2[i] = (float) v1;
        stats1.add(Math.sqrt(v1));
        stats2.add(v1);
      }
      logStats("var/g^2", stats2);
      logStats("sqrt(var/g^2)", stats1);
    } catch (final Exception ex) {
      ImageJUtils.log("Failed to load camera model %s from file: %s. %s", name,
          resource.getFilename(), ex.getMessage());
    }
    stack.addSlice("var/g^2", varG2);

    // Create as a hyper stack
    imp = new ImagePlus(name, stack);
    imp.setDimensions(4, 1, 1);
    imp.setOpenAsHyperStack(true);
    final CompositeImage cimp = new CompositeImage(imp, CompositeImage.GRAYSCALE);
    cimp.resetDisplayRanges();
    for (int n = 1; n <= 4; n++) {
      cimp.setSliceWithoutUpdate(n);
      ImageJUtils.autoAdjust(cimp, true);
    }
    cimp.setSliceWithoutUpdate(1);

    imp = WindowManager.getImage(name);
    if (imp != null) {
      imp.setImage(cimp);
    } else {
      cimp.show();
    }
  }

  private static void logStats(String name, ImageProcessor ip) {
    final ExtendedStatistics stats = new ExtendedStatistics();
    if (ip instanceof FloatProcessor) {
      stats.add((float[]) ip.getPixels());
    } else {
      for (int i = ip.getPixelCount(); i-- > 0;) {
        stats.add(ip.getf(i));
      }
    }
    logStats(name, stats);
  }

  private static void logStats(String name, ExtendedStatistics stats) {
    ImageJUtils.log("%s : %s += %s : [%s to %s]", name, MathUtils.rounded(stats.getMean()),
        MathUtils.rounded(stats.getStandardDeviation()), MathUtils.rounded(stats.getMin()),
        MathUtils.rounded(stats.getMax()));
  }

  private static void runPrintCameraModels(CameraModelSettings settings) {
    IJ.log(settings.toString());
  }
}
