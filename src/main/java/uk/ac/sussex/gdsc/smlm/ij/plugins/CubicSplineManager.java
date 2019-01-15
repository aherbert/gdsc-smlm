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

import uk.ac.sussex.gdsc.core.data.procedures.FloatStackTrivalueProcedure;
import uk.ac.sussex.gdsc.core.ij.ImageJTrackProgress;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.NonBlockingExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.logging.Ticker;
import uk.ac.sussex.gdsc.core.logging.TrackProgress;
import uk.ac.sussex.gdsc.core.math.interpolation.CustomTricubicFunction;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SoftLock;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.TurboList;
import uk.ac.sussex.gdsc.core.utils.concurrent.ConcurrencyUtils;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.CubicSplineManagerSettings;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.CubicSplineResource;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.CubicSplineSettings;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.ImagePSF;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.ImagePSFOrBuilder;
import uk.ac.sussex.gdsc.smlm.function.StandardValueProcedure;
import uk.ac.sussex.gdsc.smlm.function.cspline.CubicSplineCalculator;
import uk.ac.sussex.gdsc.smlm.function.cspline.CubicSplineData;
import uk.ac.sussex.gdsc.smlm.function.cspline.CubicSplineFunction;
import uk.ac.sussex.gdsc.smlm.function.cspline.SingleCubicSplineFunction;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.ij.utils.IJImageConverter;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.plugin.PlugIn;

import java.awt.AWTEvent;
import java.awt.Label;
import java.awt.TextField;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * This plugin handle the save and load of per-pixel spline models.
 */
public class CubicSplineManager implements PlugIn {
  /**
   * Contains the information used to represent a point spread function using a cubic spline.
   */
  public static class CubicSplinePSF {
    /** The image PSF. */
    ImagePSF imagePSF;

    /** The spline data. */
    CubicSplineData splineData;

    /**
     * Instantiates a new cubic spline PSF.
     *
     * @param imagePSF the image PSF
     * @param splineData the spline data
     * @throws IllegalArgumentException If the centre is not within the function range
     */
    public CubicSplinePSF(ImagePSF imagePSF, CubicSplineData splineData) {
      this.imagePSF = imagePSF;
      this.splineData = splineData;
      //@formatter:off
      if (imagePSF.getXCentre() < 0 ||
          imagePSF.getYCentre() < 0||
          imagePSF.getZCentre() < 0 ||
          imagePSF.getXCentre() > splineData.getMaxX() ||
          imagePSF.getYCentre() > splineData.getMaxY() ||
          imagePSF.getZCentre() > splineData.getMaxZ()
          //@formatter:on
      ) {
        throw new IllegalArgumentException("The centre is not within the function");
      }
    }

    /**
     * Creates the cubic spline function.
     *
     * @param maxy the maxy
     * @param maxx the maxx
     * @param scale the scale
     * @return the cubic spline function
     */
    public CubicSplineFunction createCubicSplineFunction(int maxy, int maxx, int scale) {
      return new SingleCubicSplineFunction(splineData, maxx, maxy, imagePSF.getXCentre(),
          imagePSF.getYCentre(), imagePSF.getZCentre(), scale);
    }

    /**
     * Gets the scale.
     *
     * @param nmPerPixel the nm per pixel
     * @return the scale
     */
    public int getScale(double nmPerPixel) {
      return CubicSplineManager.getScale(imagePSF.getPixelSize(), nmPerPixel);
    }
  }

  private static final String TITLE = "Cubic Spline Manager";
  private static String directory = "";
  private static String filename = "";
  private static String cacheName = "";
  private static CubicSplinePSF cache;

  private static CubicSplineSettings.Builder settings;

  private static CubicSplineSettings.Builder getSettings() {
    return getSettings(0);
  }

  private static CubicSplineSettings.Builder getSettings(int flags) {
    if (settings == null) {
      settings = SettingsManager.readCubicSplineSettings(flags).toBuilder();
    }
    return settings;
  }

  /**
   * Creates the cubic spline.
   *
   * @param imagePSF the image PSF details
   * @param image the image
   * @param singlePrecision Set to true to use single precision (float values) to store the cubic
   *        spline coefficients
   * @return the cubic spline PSF
   */
  public static CubicSplinePSF createCubicSpline(ImagePSFOrBuilder imagePSF, ImageStack image,
      final boolean singlePrecision) {
    final int maxx = image.getWidth();
    final int maxy = image.getHeight();
    final int maxz = image.getSize();
    final float[][] psf = new float[maxz][];
    for (int z = 0; z < maxz; z++) {
      psf[z] = IJImageConverter.getData(image.getPixels(z + 1), null);
    }

    // We reduce by a factor of 3
    final int maxi = (maxx - 1) / 3;
    final int maxj = (maxy - 1) / 3;
    final int maxk = (maxz - 1) / 3;

    final int size = maxi * maxj;
    final CustomTricubicFunction[][] splines = new CustomTricubicFunction[maxk][size];
    final Ticker ticker = Ticker.create(new ImageJTrackProgress(), (long) maxi * maxj * maxk, true);
    ticker.start();
    final ExecutorService threadPool = Executors.newFixedThreadPool(Prefs.getThreads());
    final TurboList<Future<?>> futures = new TurboList<>(maxk);
    // Create all the spline nodes by processing continuous blocks of 4x4x4 from the image stack.
    // Note that the function is enlarge x3 so a 4x4x4 block samples the voxel at [0,1/3,2/3,1]
    // in each dimension. There should be a final pixel on the end of the data for the final
    // spline node along each dimension, i.e. dimension length = n*3 + 1 with n the number of nodes.
    for (int k = 0; k < maxk; k++) {
      final int kk = k;
      futures.add(threadPool.submit(() -> {
        final CubicSplineCalculator calc = new CubicSplineCalculator();
        final double[] value = new double[64];
        final int zz = 3 * kk;
        for (int j = 0, index = 0; j < maxj; j++) {
          // 4x4 block origin in the XY data
          int index0 = 3 * j * maxx;
          for (int i = 0; i < maxi; i++, index++) {
            ticker.tick();
            int count = 0;
            for (int z = 0; z < 4; z++) {
              final float[] data = psf[zz + z];
              for (int y = 0; y < 4; y++) {
                for (int x = 0, ii = index0 + y * maxx; x < 4; x++) {
                  value[count++] = data[ii++];
                }
              }
            }
            splines[kk][index] = CustomTricubicFunction.create(calc.compute(value));
            if (singlePrecision) {
              splines[kk][index] = splines[kk][index].toSinglePrecision();
            }
            index0 += 3;
          }
        }
      }));
    }
    ticker.stop();

    ConcurrencyUtils.waitForCompletionUnchecked(futures);

    threadPool.shutdownNow();

    // Normalise
    double maxSum = 0;
    for (int k = 0; k < maxk; k++) {
      double sum = 0;
      for (int i = 0; i < size; i++) {
        sum += splines[k][i].value000();
      }
      if (maxSum < sum) {
        maxSum = sum;
      }
    }
    if (maxSum == 0) {
      throw new IllegalStateException("The cubic spline has no maximum signal");
    }
    final double scale = 1.0 / maxSum;
    for (int k = 0; k < maxk; k++) {
      for (int i = 0; i < size; i++) {
        splines[k][i].scale(scale);
      }
    }

    // Create on an integer scale
    final CubicSplineData f = new CubicSplineData(maxi, maxj, splines);

    // Create a new info with the PSF details
    final ImagePSF.Builder b = ImagePSF.newBuilder();
    b.setImageCount(imagePSF.getImageCount());
    // Reducing the image has the effect of enlarging the pixel size
    b.setPixelSize(imagePSF.getPixelSize() * 3.0);
    b.setPixelDepth(imagePSF.getPixelDepth() * 3.0);
    // The centre has to be moved as we reduced the image size by 3.
    // In the ImagePSF the XY centre puts 0.5 at the centre of the pixel.
    // The spline puts 0,0 at the centre of each pixel for convenience.
    double cx = maxi / 2.0;
    if (imagePSF.getXCentre() != 0) {
      cx = (imagePSF.getXCentre() - 0.5) / 3;
    }
    double cy = maxj / 2.0;
    if (imagePSF.getYCentre() != 0) {
      cy = (imagePSF.getYCentre() - 0.5) / 3;
    }
    double cz = maxk / 2.0;
    if (imagePSF.getZCentre() != 0) {
      cz = imagePSF.getZCentre() / 3;
    } else if (imagePSF.getCentreImage() != 0) {
      cz = (imagePSF.getCentreImage() - 1) / 3.0;
    }

    b.setXCentre(cx);
    b.setYCentre(cy);
    b.setZCentre(cz);

    return new CubicSplinePSF(b.build(), f);
  }

  /**
   * Save the spline model. The model will be named in the resources using the filename without the
   * extension or leading path entries.
   *
   * @param psfModel the spline model
   * @param filename the filename
   * @return true, if successful
   */
  public static boolean save(CubicSplinePSF psfModel, String filename) {
    if (psfModel == null || TextUtils.isNullOrEmpty(filename)) {
      return false;
    }

    // Try to save to file
    try (FileOutputStream os = new FileOutputStream(filename)) {
      final TrackProgress progress = new ImageJTrackProgress();

      psfModel.imagePSF.writeDelimitedTo(os);
      psfModel.splineData.write(os, progress);

      saveResource(psfModel, filename, getName(filename));

      return true;
    } catch (final Exception ex) {
      ImageJUtils.log("Failed to save spline model to file: %s. %s", filename, ex.getMessage());
    }

    return false;
  }

  private static String getName(String filename) {
    final File file = new File(filename);
    return ImageJUtils.removeExtension(file.getName());
  }

  private static void saveResource(CubicSplinePSF psfModel, String filename, String name) {
    final CubicSplineResource.Builder resource = CubicSplineResource.newBuilder();
    resource.setFilename(filename);
    resource.setSplineScale(psfModel.imagePSF.getPixelSize());

    final CubicSplineSettings.Builder settings = getSettings();
    settings.putCubicSplineResources(name, resource.build());
    SettingsManager.writeSettings(settings.build());

    cacheName = name;
    cache = psfModel;
  }

  /**
   * Load the spline model. Returns null if the named model does not exist. Writes to the ImageJ log
   * if a problems occurred loading the model.
   *
   * @param name the name
   * @return the per pixel spline model (or null)
   */
  public static CubicSplinePSF load(String name) {
    if (cache != null && cacheName.equals(name)) {
      return cache;
    }

    final CubicSplineSettings.Builder settings = getSettings();
    // Try and get the named resource
    final CubicSplineResource resource = settings.getCubicSplineResourcesMap().get(name);
    if (resource == null) {
      return null;
    }
    final CubicSplinePSF f = loadFromFile(name, resource.getFilename());
    if (f != null) {
      cache = f;
      cacheName = name;
    }
    return f;
  }

  private static CubicSplinePSF loadFromFile(String name, String filename) {
    // Try to load from file
    try (InputStream is = new BufferedInputStream(new FileInputStream(filename))) {
      IJ.showStatus("Loading cubic spline: " + name);
      final TrackProgress progress = new ImageJTrackProgress();

      final ImagePSF imagePSF = ImagePSF.parseDelimitedFrom(is);
      final CubicSplineData function = CubicSplineData.read(is, progress);

      return new CubicSplinePSF(imagePSF, function);
    } catch (final Exception ex) {
      ImageJUtils.log("Failed to load spline model %s from file: %s. %s", name, filename,
          ex.getMessage());
    } finally {
      IJ.showStatus("");
    }
    return null;
  }

  /**
   * List the spline models.
   *
   * @param includeNone Set to true to include an invalid none model string
   * @return the list
   */
  public static String[] listCubicSplines(boolean includeNone) {
    final CubicSplineSettings.Builder settings = getSettings();
    final List<String> list = createList(includeNone);
    list.addAll(settings.getCubicSplineResourcesMap().keySet());
    return list.toArray(new String[list.size()]);
  }

  /**
   * List the spline models with the a spline scale (pixel size) that is an integer scale of the
   * given nm/pixel scale.
   *
   * @param includeNone Set to true to include an empty string
   * @param nmPerPixel the nm per pixel
   * @return the list
   */
  public static String[] listCubicSplines(boolean includeNone, double nmPerPixel) {
    final CubicSplineSettings.Builder settings = getSettings();
    final List<String> list = createList(includeNone);
    for (final Map.Entry<String, CubicSplineResource> entry : settings.getCubicSplineResourcesMap()
        .entrySet()) {
      final CubicSplineResource resource = entry.getValue();
      final int scale = getScale(resource.getSplineScale(), nmPerPixel);
      if (scale == 0) {
        continue;
      }
      list.add(entry.getKey());
    }
    return list.toArray(new String[list.size()]);
  }

  private static List<String> createList(boolean includeNone) {
    final List<String> list = new TurboList<>();
    if (includeNone) {
      list.add("[None]");
    }
    return list;
  }

  private static int getScale(double splineSize, double nmPerPixel) {
    if (splineSize > 0) {
      final double factor = nmPerPixel / splineSize;
      final int i = (int) Math.round(factor);
      if (Math.abs(factor - i) < 1e-6) {
        return i;
      }
    }
    return 0;
  }

  //@formatter:off
  private static final String[] OPTIONS = {
      "Print all model details",
      "View a spline model",
      "Load a spline model",
      "Load from directory",
      "Delete a spline model",
      "Render the spline function" };
  //@formatter:on
  private CubicSplineManagerSettings.Builder pluginSettings;

  @Override
  public void run(String arg) {
    SMLMUsageTracker.recordPlugin(this.getClass(), arg);

    final CubicSplineSettings.Builder settings = getSettings(SettingsManager.FLAG_SILENT);
    if (settings.getCubicSplineResourcesCount() == 0) {
      IJ.error(TITLE, "No spline models found");
      return;
    }

    pluginSettings = SettingsManager.readCubicSplineManagerSettings(0).toBuilder();

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addChoice("Option", OPTIONS, pluginSettings.getOption());
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    pluginSettings.setOption(gd.getNextChoiceIndex());

    switch (pluginSettings.getOption()) {
      case 5:
        runRenderCubicSpline();
        break;
      case 4:
        runDeleteCubicSpline();
        break;
      case 3:
        runLoadFromDirectory();
        break;
      case 2:
        runLoadFromFile();
        break;
      case 1:
        runViewCubicSpline();
        break;
      default:
        runPrintCubicSplines();
    }

    SettingsManager.writeSettings(pluginSettings);
  }

  private void runRenderCubicSpline() {
    final String[] modesl = listCubicSplines(false);
    final GenericDialog gd = new GenericDialog(TITLE);
    gd.addChoice("Model", modesl, pluginSettings.getSelected());
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    final String name = gd.getNextChoice();
    pluginSettings.setSelected(name);

    final CubicSplinePSF psfModel = load(name);

    if (psfModel == null) {
      IJ.log("Failed to find spline data for model: " + name);
      return;
    }

    // Interactive render
    new CSplineRenderer(psfModel).run();
  }

  private class CSplineRenderer implements DialogListener {
    CubicSplinePSF psfModel;
    ImagePSF imagePSF;
    CubicSplineData function;
    double padx;
    double pady;
    Label label;

    public CSplineRenderer(CubicSplinePSF psfModel) {
      this.psfModel = psfModel;
      imagePSF = psfModel.imagePSF;
      function = psfModel.splineData;

      // Find the limits of the model. This works if the centre is within the image
      padx = Math.max(imagePSF.getXCentre(), function.getMaxX() - imagePSF.getXCentre());
      pady = Math.max(imagePSF.getYCentre(), function.getMaxY() - imagePSF.getYCentre());
    }

    public void run() {
      final NonBlockingExtendedGenericDialog gd = new NonBlockingExtendedGenericDialog(TITLE);
      gd.addSlider("Scale", 1, 5, pluginSettings.getScale());

      // Find the limits of the PSF
      final int width = (int) Math.ceil(function.getMaxX() * imagePSF.getPixelSize());
      final int height = (int) Math.ceil(function.getMaxY() * imagePSF.getPixelSize());
      final int minZ = (int) Math.floor(-imagePSF.getZCentre() * imagePSF.getPixelDepth());
      final int maxZ =
          (int) Math.ceil((function.getMaxZ() - imagePSF.getZCentre()) * imagePSF.getPixelDepth());

      gd.addSlider("x_shift (nm)", -width, width, pluginSettings.getXShift());
      final TextField tfxshift = gd.getLastTextField();
      gd.addSlider("y_shift (nm)", -height, height, pluginSettings.getYShift());
      final TextField tfyshift = gd.getLastTextField();
      gd.addSlider("z_shift (nm)", minZ, maxZ, pluginSettings.getZShift());
      final TextField tfzshift = gd.getLastTextField();
      gd.addDialogListener(this);
      if (ImageJUtils.isShowGenericDialog()) {
        gd.addAndGetButton("Reset", event -> {
          pluginSettings.setXShift(0);
          pluginSettings.setYShift(0);
          pluginSettings.setZShift(0);
          update();
          // The events triggered by setting these should be ignored now
          tfxshift.setText("0");
          tfyshift.setText("0");
          tfzshift.setText("0");
        });
        gd.addMessage("");
        label = gd.getLastLabel();
        draw();
      }
      gd.showDialog();
    }

    @Override
    public boolean dialogItemChanged(GenericDialog gd, AWTEvent event) {
      pluginSettings.setScale((int) gd.getNextNumber());
      pluginSettings.setXShift(gd.getNextNumber());
      pluginSettings.setYShift(gd.getNextNumber());
      pluginSettings.setZShift(gd.getNextNumber());

      update();
      return true;
    }

    private int scale = -1;
    private double xshift;
    private double yshift;
    private double zshift;
    CubicSplineFunction cspline;
    double[] params = {0, 1, 0, 0, 0}; // Set Intensity == 1

    private void draw() {
      scale = pluginSettings.getScale();
      xshift = pluginSettings.getXShift();
      yshift = pluginSettings.getYShift();
      zshift = pluginSettings.getZShift();

      boolean updateCalibration = false;
      if (cspline == null || cspline.getScale() != scale) {
        updateCalibration = true;
        final int padX = (int) Math.ceil(padx / scale);
        final int padY = (int) Math.ceil(pady / scale);

        // Create a function
        final int rangeX = 1 + 2 * padX;
        final int rangeY = 1 + 2 * padY;
        cspline = psfModel.createCubicSplineFunction(rangeX, rangeY, scale);
        params[PeakResult.X] = padX;
        params[PeakResult.Y] = padY;
      }

      // Render
      final StandardValueProcedure valueProcedure = new StandardValueProcedure();

      // Put the spot in the centre of the image
      final double[] params2 = this.params.clone();

      // Adjust the centre
      final double nmPerPixel = imagePSF.getPixelSize() * scale;
      final double nmPerSlice = imagePSF.getPixelDepth() * scale;

      params2[PeakResult.X] += xshift / nmPerPixel;
      params2[PeakResult.Y] += yshift / nmPerPixel;
      params2[PeakResult.Z] += zshift / nmPerSlice;

      final double[] values = valueProcedure.getValues(cspline, params2);

      final WindowOrganiser wo = new WindowOrganiser();
      final ImagePlus imp = ImageJUtils.display(pluginSettings.getSelected() + " (slice)", values,
          cspline.getMaxX(), cspline.getMaxY(), ImageJUtils.NO_TO_FRONT, wo);
      if (wo.size() != 0) {
        updateCalibration = true;
        imp.getWindow().toFront();
      }

      if (updateCalibration) {
        final Calibration c = imp.getLocalCalibration();
        c.setUnit("nm");
        c.pixelWidth = c.pixelHeight = nmPerPixel;
        c.pixelDepth = nmPerSlice;
      }

      if (label != null) {
        label.setText("Intensity = " + MathUtils.rounded(MathUtils.sum(values)));
      }
    }

    SoftLock lock = new SoftLock();

    private void update() {
      if (lock.acquire()) {
        // Run in a new thread to allow the GUI to continue updating
        new Thread(() -> {
          try {
            // Continue while the parameter is changing
            while (scale != pluginSettings.getScale() || xshift != pluginSettings.getXShift()
                || yshift != pluginSettings.getYShift() || zshift != pluginSettings.getZShift()) {
              draw();
            }
          } finally {
            // Ensure the running flag is reset
            lock.release();
          }
        }).start();
      }
    }
  }

  private void runDeleteCubicSpline() {
    final GenericDialog gd = new GenericDialog(TITLE);
    final String[] models = listCubicSplines(false);
    gd.addChoice("Model", models, pluginSettings.getSelected());
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    final String name = gd.getNextChoice();
    pluginSettings.setSelected(name);

    final CubicSplineResource resource = settings.getCubicSplineResourcesMap().get(name);
    if (resource == null) {
      IJ.log("Failed to find spline data for model: " + name);
      return;
    }

    settings.removeCubicSplineResources(name);
    SettingsManager.writeSettings(settings.build());

    ImageJUtils.log("Deleted spline model: %s\n%s", name, resource);
  }

  private static void runLoadFromDirectory() {
    final ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE);
    egd.addMessage("Load spline models from a directory.");
    egd.addFilenameField("Directory", directory);
    egd.showDialog();
    if (egd.wasCanceled()) {
      return;
    }

    directory = egd.getNextString();

    final File[] fileList = (new File(directory)).listFiles(File::isFile);

    for (final File file : fileList) {
      loadFromFileAndSaveResource(file.getPath());
    }
  }

  private static void runLoadFromFile() {
    final ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE);
    egd.addMessage("Load a spline model from file.");
    egd.addFilenameField("Filename", filename);
    egd.showDialog();
    if (egd.wasCanceled()) {
      return;
    }

    filename = egd.getNextString();

    loadFromFileAndSaveResource(filename);
  }

  private static void loadFromFileAndSaveResource(String filename) {
    final String name = getName(filename);
    final CubicSplinePSF model = loadFromFile(name, filename);

    if (model != null) {
      saveResource(model, filename, name);
    }
  }

  private void runViewCubicSpline() {
    final GenericDialog gd = new GenericDialog(TITLE);
    final String[] models = listCubicSplines(false);
    gd.addChoice("Model", models, pluginSettings.getSelected());
    gd.addSlider("Magnification", 1, 5, pluginSettings.getMagnification());
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    final String name = gd.getNextChoice();
    pluginSettings.setSelected(name);
    final int magnification = (int) gd.getNextNumber();
    pluginSettings.setMagnification(magnification);

    final CubicSplinePSF psfModel = load(name);

    if (psfModel == null) {
      IJ.log("Failed to find spline data for model: " + name);
      return;
    }

    IJ.showStatus("Drawing cubic spline");
    final FloatStackTrivalueProcedure p = new FloatStackTrivalueProcedure();
    psfModel.splineData.sample(magnification, p, new ImageJTrackProgress());

    final ImageStack stack = new ImageStack(p.getXAxis().length, p.getYAxis().length);
    for (final float[] pixels : p.getValue()) {
      stack.addSlice(null, pixels);
    }

    final ImagePlus imp = ImageJUtils.display(name + " (upsampled)", stack);
    final Calibration c = imp.getLocalCalibration();
    c.setUnit("nm");
    c.pixelWidth = c.pixelHeight = psfModel.imagePSF.getPixelSize() * magnification;
    c.pixelDepth = psfModel.imagePSF.getPixelDepth() * magnification;

    final int centre = 1 + (int) Math.round(psfModel.imagePSF.getZCentre() * magnification);
    imp.setSlice(centre);
    imp.resetDisplayRange();
    imp.updateAndDraw();

    IJ.showStatus("");
  }

  private static void runPrintCubicSplines() {
    IJ.log(settings.toString());
  }
}
