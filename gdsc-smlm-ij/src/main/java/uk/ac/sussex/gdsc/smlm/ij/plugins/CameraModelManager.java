/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Package
 *
 * Software for single molecule localisation microscopy (SMLM) in ImageJ
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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

import ij.CompositeImage;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PointRoi;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.Colors;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Rectangle;
import java.awt.TextField;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicReference;
import javax.swing.SwingUtilities;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.statistics.distribution.ContinuousDistribution;
import org.apache.commons.statistics.distribution.ExponentialDistribution;
import org.apache.commons.statistics.distribution.NormalDistribution;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.NonBlockingExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.utils.ExtendedStatistics;
import uk.ac.sussex.gdsc.core.utils.FileUtils;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.OpenHashMaps.CustomInt2IntOpenHashMap;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraModelResource;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraModelSettings;
import uk.ac.sussex.gdsc.smlm.ij.IJImageSource;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.CameraModelManagerSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.CameraModelManagerSettings.Builder;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.ij.utils.ImageJImageConverter;
import uk.ac.sussex.gdsc.smlm.model.camera.CameraModel;
import uk.ac.sussex.gdsc.smlm.model.camera.PerPixelCameraModel;
import uk.ac.sussex.gdsc.smlm.results.ImageSource;

/**
 * This plugin handle the save and load of per-pixel camera models.
 */
public class CameraModelManager implements PlugIn {
  private static final String TITLE = "Camera Model Manager";
  private static final String INFO_TAG = "Per-pixel camera model data";
  private static final Color LOWER_COLOR = new Color(238, 119, 51); // Orange
  private static final Color UPPER_COLOR = new Color(0, 153, 136); // Teal

  private static final AtomicReference<String> DIRECTORY = new AtomicReference<>("");
  private static final AtomicReference<String> FILENAME = new AtomicReference<>("");

  /** Cache camera models for speed. */
  private static final Map<String, PerPixelCameraModel> CAMERA_MODELS = new ConcurrentHashMap<>();

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
    CAMERA_MODELS.put(name, cameraModel);
  }

  /**
   * Load the camera model. Returns null if the named model does not exist. Writes to the ImageJ log
   * if a problems occurred loading the model.
   *
   * @param name the name
   * @return the per pixel camera model (or null)
   */
  public static PerPixelCameraModel load(String name) {
    PerPixelCameraModel model = CAMERA_MODELS.get(name);
    if (model == null) {
      final CameraModelSettings settings = CameraModelSettingsHolder.getSettings();
      // Try and get the named resource
      final CameraModelResource resource = settings.getCameraModelResourcesMap().get(name);
      if (resource == null) {
        return null;
      }
      model = loadFromFile(name, resource.getFilename());

      // Cache this
      CAMERA_MODELS.put(name, model);
    }
    return model;
  }

  private static PerPixelCameraModel loadFromFile(String name, String filename) {
    // Try and load the resource
    final ImagePlus imp = IJ.openImage(filename);
    IJ.showStatus(""); // Remove the status from the ij.io.ImageWriter class

    if (imp == null) {
      ImageJUtils.log("Failed to load camera model %s data from file: %s", name, filename);
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
    final List<String> list = new LocalList<>();
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
    gd.addHelp(HelpUrls.getUrl("camera-model-manager"));
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
        break;
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
    gd.addHelp(HelpUrls.getUrl("camera-model-manager-filter"));
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
    gd.addHelp(HelpUrls.getUrl("camera-model-manager-delete"));
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
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addMessage("Load camera models from a directory.");
    gd.addDirectoryField("Directory", DIRECTORY.get());
    gd.addHelp(HelpUrls.getUrl("camera-model-manager-load-dir"));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }

    final String dir = gd.getNextString();
    DIRECTORY.set(dir);

    final File[] fileList = (new File(dir)).listFiles(File::isFile);
    if (!ArrayUtils.isEmpty(fileList)) {
      for (final File file : fileList) {
        loadFromFileAndSaveResource(file.getPath());
      }
    }
  }

  private static void runLoadFromFile() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addMessage("Load a camera model from file.");
    gd.addFilenameField("Filename", FILENAME.get());
    gd.addHelp(HelpUrls.getUrl("camera-model-manager-load-file"));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }

    final String file = gd.getNextString();
    FILENAME.set(file);

    loadFromFileAndSaveResource(file);
  }

  private static void loadFromFileAndSaveResource(String filename) {
    final String name = getName(filename);
    final PerPixelCameraModel model = loadFromFile(name, filename);

    if (model != null) {
      saveResource(model, filename, name);
      ImageJUtils.log("Load camera model %s data from file: %s", name, filename);
    }
  }

  private void runViewCameraModel(CameraModelSettings settings) {
    final GenericDialog gd = new GenericDialog(TITLE);
    final String[] models = listCameraModels(false);
    gd.addChoice("Model", models, pluginSettings.getSelected());
    gd.addCheckbox("Show_histograms", pluginSettings.getShowHistograms());
    gd.addNumericField("Histogram_bins", pluginSettings.getHistogramBins());
    gd.addCheckbox("Outlier_analysis", pluginSettings.getShowOutliers());
    gd.addHelp(HelpUrls.getUrl("camera-model-manager-view"));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    final String name = gd.getNextChoice();
    pluginSettings.setSelected(name);
    pluginSettings.setShowHistograms(gd.getNextBoolean());
    pluginSettings.setHistogramBins(Math.max((int) gd.getNextNumber(), 0));
    pluginSettings.setShowOutliers(gd.getNextBoolean());

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
    final WindowOrganiser wo = new WindowOrganiser();
    final Plot[] plots = new Plot[3];
    final ExtendedStatistics[] stats = new ExtendedStatistics[3];
    for (int n = 1; n <= stack.getSize(); n++) {
      final ImageProcessor ip = stack.getProcessor(n);
      final ExtendedStatistics s = logStats(stack.getSliceLabel(n), ip);
      stats[n - 1] = s;
      if (pluginSettings.getShowHistograms()) {
        plots[n - 1] = CmosAnalysis.showHistogram(name, stack.getSliceLabel(n),
            ImageJImageConverter.getDoubleData(ip, null),
            // data,
            pluginSettings.getHistogramBins(), s, wo);
      }
    }
    wo.tile();

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
    if (imp == null) {
      cimp.show();
      imp = cimp;
    } else {
      imp.setImage(cimp);
    }

    if (pluginSettings.getShowOutliers()) {
      doOutlierAnalysis(plots, stats, imp, pluginSettings);
    }
  }

  private static ExtendedStatistics logStats(String name, ImageProcessor ip) {
    final ExtendedStatistics stats = new ExtendedStatistics();
    if (ip instanceof FloatProcessor) {
      stats.add((float[]) ip.getPixels());
    } else {
      for (int i = ip.getPixelCount(); i-- > 0;) {
        stats.add(ip.getf(i));
      }
    }
    logStats(name, stats);
    return stats;
  }

  private static void logStats(String name, ExtendedStatistics stats) {
    ImageJUtils.log("%s : %s += %s : [%s to %s]", name, MathUtils.rounded(stats.getMean()),
        MathUtils.rounded(stats.getStandardDeviation()), MathUtils.rounded(stats.getMin()),
        MathUtils.rounded(stats.getMax()));
  }

  /**
   * Do the outlier analysis.
   *
   * @param plots the plots (plots can be null)
   * @param stats the statistics for the camera model data
   * @param imp the image
   * @param pluginSettings the plugin settings
   */
  private static void doOutlierAnalysis(Plot[] plots, ExtendedStatistics[] stats, ImagePlus imp,
      Builder pluginSettings) {
    // Identify outliers.
    // Assume the same distribution as the simulation model (see CMOS analysis plugin)
    // n=1,2 : Bias, Gain ~ normal
    // n=3 : Variance ~ exponential
    final ContinuousDistribution[] dist =
        {NormalDistribution.of(stats[0].getMean(), stats[0].getStandardDeviation()),
            NormalDistribution.of(stats[1].getMean(), stats[1].getStandardDeviation()),
            ExponentialDistribution.of(stats[2].getMean())};

    final boolean nonBlocking = ImageJUtils.isShowGenericDialog();
    final ExtendedGenericDialog gd =
        nonBlocking ? new NonBlockingExtendedGenericDialog("Outlier analysis: " + imp.getTitle())
            : new ExtendedGenericDialog("Outlier analysis: " + imp.getTitle());
    gd.addMessage("Identify outlier pixels");
    gd.addHelp(HelpUrls.getUrl("camera-model-manager-view"));
    final double threshold = 1e-5;
    final ImageStack stack = imp.getImageStack();
    final double[] limits = new double[6];
    final LocalList<TextField> fields = new LocalList<>(6);
    final Color[] colors = {getColor(pluginSettings.getLowerColour(), UPPER_COLOR),
        getColor(pluginSettings.getUpperColour(), LOWER_COLOR)};
    for (int i = 0; i < 3; i++) {
      final double lower = dist[i].inverseCumulativeProbability(threshold);
      final double upper = dist[i].inverseSurvivalProbability(threshold);
      limits[i * 2] = lower;
      limits[i * 2 + 1] = upper;
      if (plots[i] != null) {
        // Create a snapshot to allow redraw of the original plot
        plots[i].savePlotObjects();
      }
      fields.add(bind(gd.addAndGetSlider("Lower_" + stack.getSliceLabel(i + 1), stats[i].getMin(),
          stats[i].getMax(), lower), limits, i * 2, plots[i], colors));
      fields.add(bind(gd.addAndGetSlider("Upper_" + stack.getSliceLabel(i + 1), stats[i].getMin(),
          stats[i].getMax(), upper), limits, i * 2 + 1, plots[i], colors));
    }
    // Bind to colour selection events
    addColorField(gd, "Lower_colour", colors, 0, nonBlocking, plots, limits, pluginSettings);
    addColorField(gd, "Upper_colour", colors, 1, nonBlocking, plots, limits, pluginSettings);
    if (nonBlocking) {
      ((NonBlockingExtendedGenericDialog) gd).setImage(imp);
      // Note: We cannot hide the cancel button as then the help button will close the dialog.
      // This is a feature of the ImageJ GenericDialog actions for the buttons.
      final double[] original = limits.clone();
      gd.addAndGetButton("Reset", e -> SwingUtilities.invokeLater(() -> {
        System.arraycopy(original, 0, limits, 0, original.length);
        // Note: Setting the field will trigger an update to the histogram plot via an event
        for (int i = 0; i < limits.length; i++) {
          fields.unsafeGet(i).setText(Double.toString(limits[i]));
        }
      }));
      gd.addAndGetButton("Overlay", e -> overlayOutliers(imp, limits, colors));
      gd.addAndGetButton("Clear Overlay", e -> imp.setOverlay(null));
      gd.addAndGetButton("Save", e -> saveOutliers(imp, limits, pluginSettings));
      addLimits(plots, limits, colors);
    }
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    // Read the dialog to collect/record possible macro settings.
    for (int i = 0; i < limits.length; i++) {
      limits[i] = gd.getNextNumber();
    }
    colors[0] = getColor(gd.getNextColor(), UPPER_COLOR);
    colors[1] = getColor(gd.getNextColor(), LOWER_COLOR);
    // Hash code is the colour value
    pluginSettings.setLowerColour(colors[0].hashCode());
    pluginSettings.setUpperColour(colors[1].hashCode());
    if (!nonBlocking) {
      // Run once to show the outliers
      addLimits(plots, limits, colors);
      overlayOutliers(imp, limits, colors);
    }
  }

  private static TextField bind(TextField tf, double[] limits, int i, Plot plot, Color[] colors) {
    tf.addTextListener(e -> {
      try {
        limits[i] = Double.parseDouble(tf.getText());
        if (plot != null) {
          SwingUtilities.invokeLater(() -> {
            plot.restorePlotObjects();
            final int j = (i >> 1) << 1;
            addLimits(plot, limits[j], limits[j + 1], colors);
          });
        }
      } catch (final NumberFormatException ignored) {
        // Ignore the update
      }
    });
    return tf;
  }

  private static void addLimits(Plot[] plots, double[] limits, Color[] colors) {
    for (int i = 0; i < 3; i++) {
      if (plots[i] != null) {
        addLimits(plots[i], limits[i * 2], limits[i * 2 + 1], colors);
      }
    }
  }

  private static void addLimits(Plot plot, double lower, double upper, Color[] colors) {
    plot.setColor(colors[0]);
    plot.drawLine(lower, 0, lower, 1e10);
    plot.setColor(colors[1]);
    plot.drawLine(upper, 0, upper, 1e10);
    plot.updateImage();
  }

  private static TextField addColorField(ExtendedGenericDialog gd, String label, Color[] colors,
      int index, boolean nonBlocking, Plot[] plots, double[] limits, Builder pluginSettings) {
    gd.addColorField(label, colors[index]);
    if (nonBlocking) {
      final Color original = colors[index];
      final TextField tf = gd.getLastTextField();
      tf.addTextListener(e -> {
        colors[index] = Colors.decode(tf.getText(), original);
        if (tf.isVisible()) {
          addLimits(plots, limits, colors);
          pluginSettings.setLowerColour(colors[0].hashCode());
          pluginSettings.setUpperColour(colors[1].hashCode());
          writeCameraModelManagerSettings(pluginSettings);
        }
      });
      return tf;
    }
    return null;
  }

  private static Color getColor(int color, Color defaultValue) {
    return color == 0 ? defaultValue : new Color(color);
  }

  private static Color getColor(Color color, Color defaultValue) {
    return color == null ? defaultValue : color;
  }

  private static void overlayOutliers(ImagePlus imp, double[] limits, Color[] colors) {
    final ImageStack stack = imp.getImageStack();
    final int width = stack.getWidth();
    final Overlay o = new Overlay();
    for (int slice = 1; slice <= 3; slice++) {
      final ImageProcessor ip = stack.getProcessor(slice);
      final int n = ip.getPixelCount();
      final float lower = (float) limits[2 * (slice - 1)];
      final float upper = (float) limits[2 * (slice - 1) + 1];
      float[] x = new float[100];
      float[] y = new float[x.length];
      float[] x2 = new float[100];
      float[] y2 = new float[x.length];
      int count = 0;
      int count2 = 0;
      for (int i = 0; i < n; i++) {
        final float f = ip.getf(i);
        if (f < lower) {
          x[count] = i % width;
          y[count] = i / width;
          count++;
          if (count == x.length) {
            x = Arrays.copyOf(x, count * 2);
            y = Arrays.copyOf(y, count * 2);
          }
        }
        if (f > upper) {
          x2[count2] = i % width;
          y2[count2] = i / width;
          count2++;
          if (count2 == x2.length) {
            x2 = Arrays.copyOf(x2, count2 * 2);
            y2 = Arrays.copyOf(y2, count2 * 2);
          }
        }
      }
      ImageJUtils.log("%s outside [%s to %s] = %d + %d = %d", stack.getSliceLabel(slice), lower,
          upper, count, count2, count + count2);
      if (count != 0) {
        final PointRoi roi = new PointRoi(x, y, count);
        roi.setPosition(slice);
        roi.setStrokeColor(colors[0]);
        o.add(roi);
      }
      if (count2 != 0) {
        final PointRoi roi = new PointRoi(x2, y2, count2);
        roi.setPosition(slice);
        roi.setStrokeColor(colors[1]);
        o.add(roi);
      }
    }
    imp.setOverlay(o);
    imp.getWindow().toFront();
  }

  private static void saveOutliers(ImagePlus imp, double[] limits, Builder pluginSettings) {
    final ImageStack stack = imp.getImageStack();
    final int width = stack.getWidth();
    final CustomInt2IntOpenHashMap map = new CustomInt2IntOpenHashMap();
    for (int slice = 1; slice <= 3; slice++) {
      final ImageProcessor ip = stack.getProcessor(slice);
      final int n = ip.getPixelCount();
      final float lower = (float) limits[2 * (slice - 1)];
      final float upper = (float) limits[2 * (slice - 1) + 1];
      final int bit = 1 << (slice - 1);
      for (int i = 0; i < n; i++) {
        final float f = ip.getf(i);
        if (f < lower || f > upper) {
          map.mergeInt(i, bit, (a, b) -> a | b);
        }
      }
    }
    if (map.isEmpty()) {
      return;
    }
    final String filename =
        ImageJUtils.getFilename("Outlier_filename", pluginSettings.getOutlierFilename());
    if (filename == null) {
      return;
    }
    pluginSettings.setOutlierFilename(filename);
    writeCameraModelManagerSettings(pluginSettings);
    final ImageProcessor ip1 = stack.getProcessor(1);
    final ImageProcessor ip2 = stack.getProcessor(2);
    final ImageProcessor ip3 = stack.getProcessor(3);
    try (PrintWriter out = new PrintWriter(Files.newBufferedWriter(Paths.get(filename)))) {
      out.print("x,y,outlier");
      for (int slice = 1; slice <= 3; slice++) {
        out.print(',');
        out.print(stack.getSliceLabel(slice));
      }
      out.println();
      // Unsorted output
      map.forEach((int index, int v) -> {
        out.print(index % width);
        out.print(',');
        out.print(index / width);
        out.print(',');
        out.print(v);
        out.print(',');
        out.print(ip1.getf(index));
        out.print(',');
        out.print(ip2.getf(index));
        out.print(',');
        out.print(ip3.getf(index));
        out.println();
      });
    } catch (final IOException e) {
      throw new UncheckedIOException(e);
    }
  }

  private static void runPrintCameraModels(CameraModelSettings settings) {
    IJ.log(settings.toString());
  }
}
