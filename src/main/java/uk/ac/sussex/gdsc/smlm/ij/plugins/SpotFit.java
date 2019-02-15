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

import uk.ac.sussex.gdsc.core.annotation.Nullable;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.utils.ImageExtractor;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FitSolver;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.SpotFitSettings;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.data.config.PsfProtosHelper;
import uk.ac.sussex.gdsc.smlm.engine.FitConfiguration;
import uk.ac.sussex.gdsc.smlm.engine.SimplePeakResultValidationData;
import uk.ac.sussex.gdsc.smlm.filters.BlockMeanFilter;
import uk.ac.sussex.gdsc.smlm.fitting.FitResult;
import uk.ac.sussex.gdsc.smlm.fitting.FitStatus;
import uk.ac.sussex.gdsc.smlm.fitting.Gaussian2DFitter;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.utils.ImageConverter;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.gui.Toolbar;
import ij.plugin.PlugIn;
import ij.plugin.tool.PlugInTool;
import ij.process.FloatPolygon;
import ij.process.ImageProcessor;
import ij.text.TextPanel;
import ij.text.TextWindow;

import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicReference;
import java.util.regex.Pattern;

/**
 * Plugin to fit spots interactively.
 */
public class SpotFit implements PlugIn {
  private static final String TITLE = "Spot Fit";

  /**
   * All the work for this plugin is done with the plugin tool. It handles mouse click events from
   * an image.
   */
  private static class SpotFitPluginTool extends PlugInTool {
    private static SpotFitPluginTool toolInstance = new SpotFitPluginTool();

    private static final AtomicReference<TextWindow> resultsWindow = new AtomicReference<>();
    private static final Pattern pattern = Pattern.compile("\t");

    /**
     * The single reference to the settings.
     *
     * <p>This should only be only accessed within a synchronized block.
     */
    private SpotFitSettings settingsInstance;

    /** The active flag. */
    private final AtomicBoolean active = new AtomicBoolean();

    // The following settings are used for fitting. They are all
    // used within code synchronized on the instance.
    private final BlockMeanFilter filter = new BlockMeanFilter();
    private final FitConfiguration config;
    private SimplePeakResultValidationData validationData;
    private final Gaussian2DFitter gf;
    private final double[] lower;
    private final double[] upper;

    private static class ComparisonResult {
      final int channel;
      final double background;
      final double intensity;

      ComparisonResult(int channel, double background, double intensity) {
        this.channel = channel;
        this.background = background;
        this.intensity = intensity;
      }
    }


    private SpotFitPluginTool() {
      settingsInstance = SettingsManager.readSpotFitSettings(0);
      updateActive(settingsInstance);

      config = createFitConfiguration();
      gf = new Gaussian2DFitter(config);

      // Support bounded fit on the coordinates
      final int n = Gaussian2DFunction.PARAMETERS_PER_PEAK + 1;
      lower = new double[n];
      upper = new double[n];
      for (int i = 0; i < n; i++) {
        lower[i] = Double.NEGATIVE_INFINITY;
        upper[i] = Double.POSITIVE_INFINITY;
      }
    }

    /**
     * Gets the spot fit settings. This is synchronized.
     *
     * @return the spot fit settings
     */
    private synchronized SpotFitSettings getSpotFitSettings() {
      return settingsInstance;
    }

    /**
     * Sets the spot fit settings. This is synchronized.
     *
     * @param settings the new spot fit settings
     */
    private synchronized void setSpotFitSettings(SpotFitSettings settings) {
      settingsInstance = settings;
      updateActive(settings);
    }

    private boolean isActive() {
      return active.get();
    }

    /**
     * Update the active flag using the settings.
     *
     * @param settings the settings
     */
    private void updateActive(SpotFitSettings settings) {
      active.set(settings.getFitRadius() > 1);
    }

    @Override
    public String getToolName() {
      return TITLE + " Tool";
    }

    @Override
    public String getToolIcon() {
      // A blue dot with red circle outside (a bull's eye target)
      return "C00f0o4466Cf00O11bb";
    }

    @Override
    public void showOptionsDialog() {
      final GenericDialog gd = new GenericDialog(TITLE + " Tool Options");
      gd.addMessage(TextUtils.wrap("Click on an image and fit a spot in a selected channel. "
          + "The maxima within a search range is used to centre the "
          + "fit window for Gaussian 2D fitting.", 80));
      SpotFitSettings settings = getSpotFitSettings();
      gd.addNumericField("Channel", settings.getChannel(), 0);
      gd.addSlider("Search_range", 1, 10, settings.getSearchRadius());
      gd.addSlider("Fit_radius", 3, 10, settings.getFitRadius());
      gd.addNumericField("SNR_threshold", settings.getSnrThreshold(), 0);
      gd.addCheckbox("Show_fit_ROI", settings.getShowFitRoi());
      gd.addCheckbox("Show_overlay", settings.getShowOverlay());
      gd.addCheckbox("Attach_to_slice", settings.getAttachToSlice());
      gd.addCheckbox("Log_progress", settings.getLogProgress());
      gd.addMessage(TextUtils.wrap(
          "Optionally perform a weighted mean intensity in a second "
              + "channel using the fitted Gaussian to weight the region. Channel 0 is ignored",
          80));
      gd.addNumericField("Comparison_channel", settings.getComparisonChannel(), 0);

      gd.showDialog();
      if (gd.wasCanceled()) {
        return;
      }

      final SpotFitSettings.Builder builder = settings.toBuilder();
      builder.setChannel(Math.max(1, (int) gd.getNextNumber()));
      builder.setSearchRadius((int) gd.getNextNumber());
      builder.setFitRadius((int) gd.getNextNumber());
      builder.setSnrThreshold(gd.getNextNumber());
      builder.setShowFitRoi(gd.getNextBoolean());
      builder.setShowOverlay(gd.getNextBoolean());
      builder.setAttachToSlice(gd.getNextBoolean());
      builder.setLogProgress(gd.getNextBoolean());
      builder.setComparisonChannel(Math.max(0, (int) gd.getNextNumber()));
      settings = builder.build();

      // Only active if the settings are valid
      setSpotFitSettings(settings);
      if (!isActive()) {
        IJ.error(TITLE, "Settings are invalid");
      }

      SettingsManager.writeSettings(settings, 0);
    }

    @Override
    public void mouseClicked(ImagePlus imp, MouseEvent event) {
      if (!isActive()) {
        return;
      }

      final SpotFitSettings settings = getSpotFitSettings();

      if (settings.getChannel() > imp.getNChannels()) {
        // Always warn if the channel is incorrect for the image
        ImageJUtils.log(TITLE + ": Image %s does not contain channel %d", imp.getTitle(),
            settings.getChannel());
        return;
      }

      // Mark this event as handled
      event.consume();

      // TODO - More control over fitting.

      // Ensure rapid mouse click / new options does not break things
      synchronized (this) {
        final ImageCanvas ic = imp.getCanvas();
        int x = ic.offScreenX(event.getX());
        int y = ic.offScreenY(event.getY());

        if (settings.getLogProgress()) {
          ImageJUtils.log("Clicked %d,%d", x, y);
        }

        // Get the data
        final int channel = settings.getChannel();
        final int slice = imp.getZ();
        final int frame = imp.getFrame();

        final int stackIndex = imp.getStackIndex(channel, slice, frame);

        final ImageExtractor ie = ImageExtractor.wrap(null, imp.getWidth(), imp.getHeight());

        if (isRemoveEvent(event)) {
          removeSpots(imp, channel, slice, frame, x, y, ie, settings.getSearchRadius());
          return;
        }

        final ImageStack stack = imp.getImageStack();
        final ImageProcessor ip = stack.getProcessor(stackIndex);

        // Search for the maxima using the search radius
        final int index = findMaxima(ip, ie, x, y, settings.getSearchRadius());

        // Fit the maxima
        x = index % ip.getWidth();
        y = index / ip.getWidth();
        if (settings.getLogProgress()) {
          ImageJUtils.log("Fitting %d,%d", x, y);
        }
        final Rectangle bounds = ie.getBoxRegionBounds(x, y, settings.getFitRadius());
        if (settings.getShowFitRoi()) {
          imp.setRoi(bounds);
        }
        final FitResult fitResult = fitMaxima(ip, bounds, x, y, settings);

        if (settings.getLogProgress()) {
          ImageJUtils.log("Fit estimate = %s", Arrays.toString(fitResult.getInitialParameters()));
        }
        if (settings.getLogProgress()) {
          String msg = "Fit status = " + fitResult.getStatus();
          final Object data = fitResult.getStatusData();
          if (data != null) {
            msg += " : " + SimpleArrayUtils.toString(data);
          }
          ImageJUtils.log(msg);
        }

        if (fitResult.getStatus() != FitStatus.OK) {
          // Q. Do something?
          return;
        }

        // Perform a weighted mean using the fitted Gaussian on the comparison channel
        final ComparisonResult comparisonResult = createComparisonResult(imp, settings, fitResult);

        // Add result
        addResult(imp, channel, slice, frame, bounds, fitResult, comparisonResult);

        if (settings.getShowOverlay()) {
          addOverlay(imp, channel, slice, frame, fitResult, settings.getAttachToSlice());
        }
      }
    }

    private static boolean isRemoveEvent(MouseEvent event) {
      return event.isAltDown() || event.isShiftDown() || event.isControlDown();
    }

    // "data" will not be null as the width and height from the image processor are correct
    @SuppressWarnings("null")
    private int findMaxima(ImageProcessor ip, ImageExtractor ie, int x, int y, int searchRadius) {
      if (searchRadius <= 0) {
        // No search
        return y * ip.getWidth() + x;
      }
      // Get a region
      final Rectangle bounds = ie.getBoxRegionBounds(x, y, searchRadius);
      final float[] data =
          new ImageConverter().getData(ip.getPixels(), ip.getWidth(), ip.getHeight(), bounds, null);
      // Smooth
      filter.blockFilter(data, bounds.width, bounds.height, 1);
      int index = 0;
      // Find maxima
      for (int i = 1; i < data.length; i++) {
        if (data[index] < data[i]) {
          index = i;
        }
      }
      // Convert back
      x = bounds.x + index % bounds.width;
      y = bounds.y + index / bounds.width;
      return ip.getWidth() * y + x;
    }

    private FitResult fitMaxima(ImageProcessor ip, Rectangle bounds, int x, int y,
        SpotFitSettings settings) {
      config.setInitialPeakStdDev(0);

      // Get a region
      final double[] data = new ImageConverter().getDoubleData(ip.getPixels(), ip.getWidth(),
          ip.getHeight(), bounds, null);

      setupPeakFiltering(config, bounds, ip, settings);

      // Find the index of the maxima
      final int ox = x - bounds.x;
      final int oy = y - bounds.y;
      final int index = oy * bounds.width + ox;

      // Limit the range for the XY position
      final double range = Math.max(1, settings.getFitRadius() / 2);
      lower[Gaussian2DFunction.X_POSITION] = ox - range;
      lower[Gaussian2DFunction.Y_POSITION] = oy - range;
      upper[Gaussian2DFunction.X_POSITION] = ox + range;
      upper[Gaussian2DFunction.Y_POSITION] = oy + range;
      gf.setBounds(lower, upper);

      // Leave to the fitter to estimate background, width and height
      final FitResult fitResult = gf.fit(data, bounds.width, bounds.height, new int[] {index});
      if (fitResult.getStatus() == FitStatus.OK) {
        final double[] params = fitResult.getParameters();
        // Add the pixel offset
        params[Gaussian2DFunction.X_POSITION] += bounds.x + 0.5;
        params[Gaussian2DFunction.Y_POSITION] += bounds.y + 0.5;
      }
      return fitResult;
    }

    private static FitConfiguration createFitConfiguration() {
      final FitConfiguration config = new FitConfiguration();
      config.setFitSolver(FitSolver.LVM_LSE);
      config.setPsf(PsfProtosHelper.getDefaultPsf(PSFType.ONE_AXIS_GAUSSIAN_2D));
      config.setInitialPeakStdDev(0);
      config.setComputeDeviations(false);

      // Set-up peak filtering
      config.setDisableSimpleFilter(false);

      config.setBackgroundFitting(true);

      return config;
    }

    protected void setupPeakFiltering(FitConfiguration config, Rectangle bounds, ImageProcessor ip,
        SpotFitSettings settings) {
      config.setCoordinateShift(settings.getSearchRadius());
      config.setSignalStrength(settings.getSnrThreshold());
      config.setMinWidthFactor(0.5);

      validationData = new SimplePeakResultValidationData(
          new GaussianFunctionFactory(config.getFunctionFlags(), config.getAstigmatismZModel()),
          bounds.x, bounds.y, ip.getPixels(), ip.getWidth(), ip.getHeight());
      config.setPeakResultValidationData(validationData);
    }

    /**
     * Create the result window (if it is not available).
     *
     * @return the text window
     */
    private static TextWindow createResultsWindow() {
      TextWindow window = resultsWindow.get();
      if (window == null || !window.isShowing()) {
        window = new TextWindow(TITLE + " Results", createHeader(), "", 700, 300);
        resultsWindow.set(window);
      }
      return window;
    }

    private static String createHeader() {
      final StringBuilder sb = new StringBuilder();
      sb.append("Image\t");
      sb.append("C\t");
      sb.append("Z\t");
      sb.append("T\t");
      sb.append("Region\t");
      sb.append("Background\t");
      sb.append("Intensity\t");
      sb.append("X\t");
      sb.append("Y\t");
      sb.append("S\t");
      sb.append("Noise\t");
      sb.append("SNR\t");
      sb.append("C2\t");
      sb.append("C2 Background\t");
      sb.append("C2 Intensity\t");
      return sb.toString();
    }

    private void addResult(ImagePlus imp, int channel, int slice, int frame, Rectangle bounds,
        FitResult fitResult, ComparisonResult comparisonResult) {
      final StringBuilder sb = new StringBuilder();
      sb.append(imp.getTitle()).append('\t');
      sb.append(channel).append('\t');
      sb.append(slice).append('\t');
      sb.append(frame).append('\t');
      sb.append(bounds.x).append(',');
      sb.append(bounds.y).append(' ');
      sb.append(bounds.width).append('x');
      sb.append(bounds.height);
      final double[] params = fitResult.getParameters();
      sb.append('\t').append(MathUtils.rounded(params[Gaussian2DFunction.BACKGROUND]));
      final double signal = params[Gaussian2DFunction.SIGNAL];
      sb.append('\t').append(MathUtils.rounded(signal));
      sb.append('\t').append(MathUtils.rounded(params[Gaussian2DFunction.X_POSITION]));
      sb.append('\t').append(MathUtils.rounded(params[Gaussian2DFunction.Y_POSITION]));
      final double xsd = params[Gaussian2DFunction.X_SD];
      sb.append('\t').append(MathUtils.rounded(xsd));
      final double noise = validationData.getNoise();
      final double snr = Gaussian2DPeakResultHelper.getMeanSignalUsingP05(signal, xsd, xsd) / noise;
      sb.append('\t').append(MathUtils.rounded(noise));
      sb.append('\t').append(MathUtils.rounded(snr));
      if (comparisonResult != null) {
        sb.append('\t').append(comparisonResult.channel);
        sb.append('\t').append(MathUtils.rounded(comparisonResult.background));
        sb.append('\t').append(MathUtils.rounded(comparisonResult.intensity));
      } else {
        sb.append("\t\t\t");
      }

      createResultsWindow().append(sb.toString());
    }

    private static void addOverlay(ImagePlus imp, int channel, int slice, int frame,
        FitResult fitResult, boolean attachToSlice) {
      final double[] params = fitResult.getParameters();
      Overlay overlay = imp.getOverlay();
      if (overlay == null) {
        overlay = new Overlay();
      }
      final PointRoi roi = new PointRoi(params[Gaussian2DFunction.X_POSITION],
          params[Gaussian2DFunction.Y_POSITION]);
      roi.setPointType(3);
      if (imp.isDisplayedHyperStack()) {
        roi.setPosition(channel, (attachToSlice) ? slice : 0, frame);
      } else if (attachToSlice) {
        roi.setPosition(imp.getCurrentSlice());
      }
      overlay.add(roi);
      imp.setOverlay(overlay);
    }

    private static void removeSpots(ImagePlus imp, int channel, int slice, int frame, int x, int y,
        ImageExtractor ie, int searchRadius) {
      // Get a region to search for spots
      final Rectangle bounds = ie.getBoxRegionBounds(x, y, Math.max(0, searchRadius));
      if (bounds.width == 0 || bounds.height == 0) {
        return;
      }

      final boolean isDisplayedHyperStack = imp.isDisplayedHyperStack();

      // Remove all the overlay components
      Overlay overlay = imp.getOverlay();
      if (overlay != null) {
        final Roi[] rois = overlay.toArray();
        final int size = overlay.size();
        overlay = new Overlay();
        for (int i = 0; i < rois.length; i++) {
          if (rois[i] instanceof PointRoi) {
            final PointRoi roi = (PointRoi) rois[i];
            boolean boundsCheck = true;
            if (isDisplayedHyperStack) {
              // Must be on the same channel/slice/frame
              boundsCheck = roi.getCPosition() == channel && roi.getTPosition() == frame
                  && (roi.getZPosition() == 0 || roi.getZPosition() == slice);
            }
            if (boundsCheck) {
              final FloatPolygon poly = roi.getFloatPolygon();
              if (bounds.contains(poly.xpoints[0], poly.ypoints[0])) {
                continue;
              }
            }
          }
          overlay.add(rois[i]);
        }
        if (overlay.size() != size) {
          if (overlay.size() == 0) {
            imp.setOverlay(null);
          } else {
            imp.setOverlay(overlay);
          }
        }
      }

      final TextWindow window = resultsWindow.get();
      if (window != null && window.isShowing()) {
        final TextPanel tp = window.getTextPanel();
        final String title = imp.getTitle();
        for (int i = 0; i < tp.getLineCount(); i++) {
          final String line = tp.getLine(i);
          // Check the image name
          final int startIndex = line.indexOf('\t');
          if (startIndex == -1 || !title.equals(line.substring(0, startIndex))) {
            continue;
          }

          final String[] fields = pattern.split(line, 0);

          try {
            if (isCorrectSlice(channel, frame, isDisplayedHyperStack, fields)) {
              final float xp = Float.parseFloat(fields[7]);
              final float yp = Float.parseFloat(fields[8]);
              if (bounds.contains(xp, yp)) {
                tp.setSelection(i, i);
                tp.clearSelection();
                // Since i will be incremented for the next line,
                // decrement to check the current line again.
                i--;
              }
            }
          } catch (final NumberFormatException ex) {
            // Ignore
          }
        }
      }
    }

    private static boolean isCorrectSlice(int channel, int frame, boolean isDisplayedHyperStack,
        String[] fields) {
      if (isDisplayedHyperStack) {
        final int c = Integer.parseInt(fields[1]);
        // Ignore z as the user may be on the wrong slice but can still see
        // the overlay if it is not tied to the slice position
        final int t = Integer.parseInt(fields[3]);
        if (c != channel || t != frame) {
          return false;
        }
      }
      return true;
    }

    // "data" will not be null as the width and height from the image processor are correct
    @SuppressWarnings("null")
    private static @Nullable ComparisonResult createComparisonResult(ImagePlus imp,
        SpotFitSettings settings, FitResult fitResult) {
      if (settings.getComparisonChannel() == 0
          || settings.getComparisonChannel() > imp.getNChannels()) {
        // No comparison channel
        return null;
      }

      final int channel = settings.getComparisonChannel();
      final int slice = imp.getZ();
      final int frame = imp.getFrame();

      final int stackIndex = imp.getStackIndex(channel, slice, frame);

      final ImageStack stack = imp.getImageStack();
      final ImageProcessor ip = stack.getProcessor(stackIndex);

      // Extract the region
      final double[] params = fitResult.getParameters();
      final double x = params[Gaussian2DFunction.X_POSITION];
      final double y = params[Gaussian2DFunction.Y_POSITION];
      final double xsd = params[Gaussian2DFunction.X_SD];

      final int cx = (int) x;
      final int cy = (int) y;
      final int width = (int) Math.ceil(3 * xsd);
      final int ox = cx - width;
      final int oy = cy - width;
      // Clip to the image
      final Rectangle bounds = new Rectangle(ox, oy, 2 * width + 1, 2 * width + 1)
          .intersection(new Rectangle(imp.getWidth(), imp.getHeight()));

      final double[] data = new ImageConverter().getDoubleData(ip.getPixels(), ip.getWidth(),
          ip.getHeight(), bounds, null);

      // Find the background using the edge pixels assuming a single peak
      final double background =
          Gaussian2DFitter.getBackground(data, bounds.width, bounds.height, 1);

      // Find the weighted intensity using the 2D Gaussian
      final double[] gaussParams = new double[Gaussian2DFunction.PARAMETERS_PER_PEAK];
      gaussParams[Gaussian2DFunction.SIGNAL] = 1;
      gaussParams[Gaussian2DFunction.X_POSITION] = x - bounds.x;
      gaussParams[Gaussian2DFunction.Y_POSITION] = y - bounds.y;
      gaussParams[Gaussian2DFunction.X_SD] = xsd;
      final double[] weights = GaussianFunctionFactory
          .create2D(1, bounds.width, bounds.height, GaussianFunctionFactory.FIT_ERF_NB_CIRCLE, null)
          .computeValues(gaussParams);

      double intensity = 0;
      // The weights should sum to 1 but the bounds may be clipped so compute
      // the sum for normalisation.
      double sumWeights = 0;
      for (int i = 0; i < data.length; i++) {
        // Clip the data below the background to zero.
        sumWeights += weights[i];
        intensity += weights[i] * Math.max(0, data[i] - background);
      }

      // If there is no Gaussian then return null
      if (sumWeights == 0) {
        return null;
      }
      return new ComparisonResult(channel, background, intensity / sumWeights);
    }
  }

  /**
   * Initialise the spot fit tool. This is to allow support for calling within macro toolsets.
   */
  public static void addPluginTool() {
    // Add the tool
    Toolbar.addPlugInTool(SpotFitPluginTool.toolInstance);
    IJ.showStatus("Added " + TITLE + " Tool");
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    addPluginTool();

    // Fiji restores the toolbar from the last session.
    // Do not show the options if this is happening.
    final ImageJ ij = IJ.getInstance();
    if (ij == null || !ij.isVisible()) {
      return;
    }

    SpotFitPluginTool.toolInstance.showOptionsDialog();
  }
}
