/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Package
 *
 * Software for single molecule localisation microscopy (SMLM) in ImageJ
 * %%
 * Copyright (C) 2011 - 2025 Alex Herbert
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
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.gui.Toolbar;
import ij.plugin.PlugIn;
import ij.plugin.tool.PlugInTool;
import ij.process.ByteProcessor;
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
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import uk.ac.sussex.gdsc.core.annotation.Nullable;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ObjectOutliner;
import uk.ac.sussex.gdsc.core.ij.gui.OffsetPointRoi;
import uk.ac.sussex.gdsc.core.math.RadialStatisticsUtils;
import uk.ac.sussex.gdsc.core.utils.FastCorrelator;
import uk.ac.sussex.gdsc.core.utils.ImageExtractor;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FitSolver;
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
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.SpotFitSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.ij.utils.TextPanelMouseListener;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.utils.ImageConverter;

/**
 * Plugin to fit spots interactively.
 */
public class SpotFit implements PlugIn {
  private static final String TITLE = "Spot Fit";

  /**
   * Custom text window to redraw the fit region when the result is double-clicked in the table.
   */
  private static class CustomTextWindow extends TextWindow {
    private static final long serialVersionUID = 1L;

    boolean draw;

    @SuppressWarnings("unused")
    CustomTextWindow(String title, String headings, int width, int height) {
      super(title, headings, "", width, height);

      // Operates on the provided text panel
      new TextPanelMouseListener(getTextPanel()) {
        @Override
        public void selected(int selectionStart, int selectionEnd) {
          // Only if drawing fit region is enabled
          if (!draw || selectionStart < 0 || selectionStart >= textPanel.getLineCount()) {
            return;
          }
          final String line = this.textPanel.getLine(selectionStart);
          String[] fields = line.split("\t");
          // Require Image, C, Z, T, Region
          if (fields.length < 5) {
            return;
          }
          final ImagePlus imp = WindowManager.getImage(fields[0]);
          if (imp == null) {
            return;
          }
          try {
            final int c = Integer.parseInt(fields[1]);
            final int z = Integer.parseInt(fields[2]);
            final int t = Integer.parseInt(fields[3]);
            // Parse the rectangle data
            fields = fields[4].split("[, x]");
            final int x = Integer.parseInt(fields[0]);
            final int y = Integer.parseInt(fields[1]);
            final int w = Integer.parseInt(fields[2]);

            imp.setPosition(c, z, t);
            imp.setRoi(new Rectangle(x, y, w, w));
          } catch (final NumberFormatException ex) {
            return;
          }
        }
      };
    }
  }

  /**
   * All the work for this plugin is done with the plugin tool. It handles mouse click events from
   * an image.
   */
  private static class SpotFitPluginTool extends PlugInTool {
    static final SpotFitPluginTool INSTANCE = new SpotFitPluginTool();

    private static final AtomicReference<CustomTextWindow> RESULTS_WINDOW = new AtomicReference<>();
    private static final Pattern PATTERN = Pattern.compile("\t");
    private static final AtomicReference<Pair<Integer, Roi>> ANALYSIS_OUTLINE =
        new AtomicReference<>(Pair.of(0, null));

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

    /**
     * A result of comparison between 2 channels.
     */
    private static class ComparisonResult {
      final int channel;
      final double background;
      final double intensity;
      final double m1;
      final double m2;
      final double corr;
      final double[] r1;
      final double[] r2;

      ComparisonResult(int channel, double background, double intensity, double m1, double m2,
          double corr, double[] r1, double[] r2) {
        this.channel = channel;
        this.background = background;
        this.intensity = intensity;
        this.m1 = m1;
        this.m2 = m2;
        this.corr = corr;
        this.r1 = r1;
        this.r2 = r2;
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
      final CustomTextWindow tw = RESULTS_WINDOW.get();
      if (tw != null) {
        tw.draw = settings.getShowFitRoi();
      }
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
          "Optionally perform comparison analysis in a second channel. Channel 0 is ignored.", 80));
      gd.addNumericField("Comparison_channel", settings.getComparisonChannel(), 0);
      gd.addSlider("Analysis_radius", 3, 10, settings.getAnalysisRadius());

      gd.addHelp(HelpUrls.getUrl("spot-fit-tool"));
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
      builder.setAnalysisRadius((int) gd.getNextNumber());
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

        final ImageExtractor ie = ImageExtractor.wrap(null, imp.getWidth(), imp.getHeight());

        if (isRemoveEvent(event)) {
          removeSpots(imp, channel, slice, frame, x, y, ie, settings.getSearchRadius());
          return;
        }

        final int stackIndex = imp.getStackIndex(channel, slice, frame);
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

        // Perform comparison analysis against a second channel
        final ComparisonResult comparisonResult =
            createComparisonResult(ip, imp, settings, fitResult);

        // Add result
        addResult(imp, channel, slice, frame, bounds, fitResult, comparisonResult,
            settings.getShowFitRoi());

        if (settings.getShowOverlay()) {
          addOverlay(imp, channel, slice, frame, fitResult, settings.getAttachToSlice(),
              settings.getAnalysisRadius());
        }
      }
    }

    private static boolean isRemoveEvent(MouseEvent event) {
      return event.isAltDown() || event.isShiftDown() || event.isControlDown();
    }

    // "data" will not be null as the width and height from the image processor are correct
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
      final FitConfiguration config = FitConfiguration.create();
      config.setFitSolver(FitSolver.LVM_LSE);
      config.setPsf(PsfProtosHelper.getDefaultPsf(PSFType.ONE_AXIS_GAUSSIAN_2D));
      config.setInitialPeakStdDev(0);
      config.setComputeDeviations(false);

      // Set-up peak filtering
      config.setDisableSimpleFilter(false);

      config.setBackgroundFitting(true);

      return config;
    }

    void setupPeakFiltering(FitConfiguration config, Rectangle bounds, ImageProcessor ip,
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
     * @param drawSelected Set to true to enable drawing of the selected results
     * @return the text window
     */
    private static CustomTextWindow createResultsWindow(boolean drawSelected) {
      return ImageJUtils.refresh(RESULTS_WINDOW, () -> {
        final CustomTextWindow tw =
            new CustomTextWindow(TITLE + " Results", createHeader(), 700, 300);
        tw.draw = drawSelected;
        return tw;
      });
    }

    private static String createHeader() {
      final StringBuilder sb = new StringBuilder(256);
      //@formatter:off
      final String[] names = {
        "Image",
        "C",
        "Z",
        "T",
        "Region",
        "Background",
        "Intensity",
        "X",
        "Y",
        "S",
        "Mean (PWHM)",
        "Noise",
        "SNR",
        "C2",
        "C2 Background",
        "C2 Intensity",
        "Corr",
        "C1 Mean",
        "C2 Mean",
        "C1 Radial means",
        "C2 Radial means",
      };
      //@formatter:on
      Arrays.stream(names).forEach(s -> sb.append(s).append('\t'));
      sb.setLength(sb.length() - 1);
      return sb.toString();
    }

    private void addResult(ImagePlus imp, int channel, int slice, int frame, Rectangle bounds,
        FitResult fitResult, ComparisonResult comparisonResult, boolean drawSelected) {
      final StringBuilder sb = new StringBuilder(256);
      //@formatter:off
      sb.append(imp.getTitle()).append('\t')
        .append(channel).append('\t')
        .append(slice).append('\t')
        .append(frame).append('\t')
        .append(bounds.x).append(',')
        .append(bounds.y).append(' ')
        .append(bounds.width).append('x')
        .append(bounds.height);
      final double[] params = fitResult.getParameters();
      sb.append('\t').append(MathUtils.rounded(params[Gaussian2DFunction.BACKGROUND]));
      final double signal = params[Gaussian2DFunction.SIGNAL];
      sb.append('\t').append(MathUtils.rounded(signal))
        .append('\t').append(MathUtils.rounded(params[Gaussian2DFunction.X_POSITION]))
        .append('\t').append(MathUtils.rounded(params[Gaussian2DFunction.Y_POSITION]));
      final double xsd = params[Gaussian2DFunction.X_SD];
      sb.append('\t').append(MathUtils.rounded(xsd));
      final double noise = validationData.getNoise();
      final double mean = Gaussian2DPeakResultHelper.getMeanSignalUsingP05(signal, xsd, xsd);
      final double snr = mean / noise;
      sb.append('\t').append(MathUtils.rounded(mean))
        .append('\t').append(MathUtils.rounded(noise))
        .append('\t').append(MathUtils.rounded(snr));
      if (comparisonResult != null) {
        sb.append('\t').append(comparisonResult.channel)
          .append('\t').append(MathUtils.rounded(comparisonResult.background))
          .append('\t').append(MathUtils.rounded(comparisonResult.intensity))
          .append('\t').append(MathUtils.rounded(comparisonResult.corr))
          .append('\t').append(MathUtils.rounded(comparisonResult.m1))
          .append('\t').append(MathUtils.rounded(comparisonResult.m2))
          .append('\t').append(Arrays.stream(comparisonResult.r1)
            .mapToObj(i -> MathUtils.rounded(i)).collect(Collectors.joining(", ")))
          .append('\t').append(Arrays.stream(comparisonResult.r2)
            .mapToObj(i -> MathUtils.rounded(i)).collect(Collectors.joining(", ")));
        //@formatter:on
      } else {
        sb.append("\t\t\t\t\t\t\t\t");
      }

      // TODO Experiments with the radial means plot ...
      // Plot p = new Plot("", "", "");
      // p.addPoints(x, y, shape);

      createResultsWindow(drawSelected).append(sb.toString());
    }

    private static void addOverlay(ImagePlus imp, int channel, int slice, int frame,
        FitResult fitResult, boolean attachToSlice, int analysisRadius) {
      final double[] params = fitResult.getParameters();
      Overlay overlay = imp.getOverlay();
      if (overlay == null) {
        overlay = new Overlay();
      }
      final PointRoi roi = new OffsetPointRoi(params[Gaussian2DFunction.X_POSITION],
          params[Gaussian2DFunction.Y_POSITION]);
      roi.setPointType(3);
      setPosition(imp, channel, slice, frame, attachToSlice, roi);
      overlay.add(roi);

      final Roi analysisRoi = createAnalysisRoi(analysisRadius);
      if (analysisRoi != null) {
        final int ox = (int) params[Gaussian2DFunction.X_POSITION];
        final int oy = (int) params[Gaussian2DFunction.Y_POSITION];
        analysisRoi.setLocation(ox - analysisRadius, oy - analysisRadius);
        setPosition(imp, channel, slice, frame, attachToSlice, analysisRoi);
        overlay.add(analysisRoi);
      }

      imp.setOverlay(overlay);
    }

    private static void setPosition(ImagePlus imp, int channel, int slice, int frame,
        boolean attachToSlice, final Roi roi) {
      if (imp.isDisplayedHyperStack()) {
        roi.setPosition(channel, (attachToSlice) ? slice : 0, frame);
      } else if (attachToSlice) {
        roi.setPosition(imp.getStackIndex(channel, slice, frame));
      }
    }

    private static Roi createAnalysisRoi(int analysisRadius) {
      final Pair<Integer, Roi> outline = ANALYSIS_OUTLINE.get();
      Roi roi;
      if (outline.getKey() != analysisRadius) {
        final int size = 2 * analysisRadius + 1;
        final ByteProcessor bp = new ByteProcessor(size, size);
        // Create offsets
        final double[] d2 = IntStream.rangeClosed(-analysisRadius, analysisRadius)
            .mapToDouble(i -> i * i).toArray();
        // Analysis of pixels within radius r in [n, n+1). Set max n = analysisRadius.
        final double threshold = MathUtils.pow2(analysisRadius + 1);
        for (int y = 0, i = 0; y < size; y++) {
          for (int x = 0; x < size; x++, i++) {
            if (d2[y] + d2[x] < threshold) {
              bp.set(i, 255);
            }
          }
        }
        final ObjectOutliner oo = new ObjectOutliner(bp);
        roi = oo.outline(analysisRadius, analysisRadius);
        ANALYSIS_OUTLINE.set(Pair.of(analysisRadius, roi));
      } else {
        roi = outline.getValue();
      }
      return (Roi) roi.clone();
    }

    private static void removeSpots(ImagePlus imp, int channel, int slice, int frame, int x, int y,
        ImageExtractor ie, int searchRadius) {
      // Get a region to search for spots
      final Rectangle bounds = ie.getBoxRegionBounds(x, y, Math.max(0, searchRadius));
      if (bounds.width == 0 || bounds.height == 0) {
        return;
      }

      final boolean isDisplayedHyperStack = imp.isDisplayedHyperStack();

      final int currentSlice = imp.getStackIndex(channel, slice, frame);

      // Remove all the overlay components
      Overlay overlay = imp.getOverlay();
      if (overlay != null) {
        final Roi[] rois = overlay.toArray();
        final int size = overlay.size();
        overlay = new Overlay();
        for (final Roi nextRoi : rois) {
          if (nextRoi instanceof PointRoi) {
            final PointRoi roi = (PointRoi) nextRoi;
            boolean boundsCheck = true;
            if (isDisplayedHyperStack) {
              // Must be on the same channel/slice/frame
              boundsCheck = roi.getCPosition() == channel && roi.getTPosition() == frame
                  && (roi.getZPosition() == 0 || roi.getZPosition() == slice);
            } else if (roi.getPosition() != 0) {
              // Must be on the same slice
              boundsCheck = roi.getPosition() == currentSlice;
            }
            if (boundsCheck) {
              final FloatPolygon poly = roi.getFloatPolygon();
              if (bounds.contains(poly.xpoints[0], poly.ypoints[0])) {
                continue;
              }
            }
          }
          overlay.add(nextRoi);
        }
        if (overlay.size() != size) {
          if (overlay.size() == 0) {
            imp.setOverlay(null);
          } else {
            imp.setOverlay(overlay);
          }
        }
      }

      final TextWindow window = RESULTS_WINDOW.get();
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

          final String[] fields = PATTERN.split(line, 0);

          try {
            if (isCorrectSlice(channel, slice, frame, isDisplayedHyperStack, fields)) {
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
          } catch (final NumberFormatException ignored) {
            // Ignore
          }
        }
      }
    }

    private static boolean isCorrectSlice(int channel, int slice, int frame,
        boolean isDisplayedHyperStack, String[] fields) {
      // Match channel and frame
      // Ignore z for hyperstacks as the overlay may not be tied to the slice position.
      // For standard stacks the click location should be 1 for 2 of 3 out of CZT,
      // i.e. c*z*t = stack index (1-based), so match all.
      return isMatch(fields, 3, frame) && isMatch(fields, 1, channel)
          && (isDisplayedHyperStack || isMatch(fields, 2, slice));
    }

    private static boolean isMatch(String[] fields, int index, int value) {
      return Integer.parseInt(fields[index]) == value;
    }

    // "data" will not be null as the width and height from the image processor are correct
    private static @Nullable ComparisonResult createComparisonResult(ImageProcessor ip1,
        ImagePlus imp, SpotFitSettings settings, FitResult fitResult) {
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
      final ImageProcessor ip2 = stack.getProcessor(stackIndex);

      // Use the fitted Gaussian to compute a single pixel Gaussian blur on the comparison channel.
      // This creates an approximate value for the other channel at the fit location.

      // Extract the fitted region
      final double[] params = fitResult.getParameters();
      final double x = params[Gaussian2DFunction.X_POSITION];
      final double y = params[Gaussian2DFunction.Y_POSITION];
      final double xsd = params[Gaussian2DFunction.X_SD];

      final int cx = (int) Math.round(x);
      final int cy = (int) Math.round(y);
      Rectangle bounds = getBounds(cx, cy, (int) Math.ceil(3 * xsd), imp);

      final ImageConverter ic = new ImageConverter();
      final double[] data =
          ic.getDoubleData(ip2.getPixels(), ip2.getWidth(), ip2.getHeight(), bounds, null);

      // Find the background using the edge pixels assuming a single peak
      final double background =
          Gaussian2DFitter.getBackground(data, bounds.width, bounds.height, 1);

      final double intensity = getIntensity(cx, cy, xsd, bounds, data, background);

      // Extract the pixels for each region using the analysis radius
      bounds = getBounds(cx, cy, settings.getAnalysisRadius(), imp);

      final float[] data1 =
          ic.getData(ip1.getPixels(), ip1.getWidth(), ip1.getHeight(), bounds, null);
      final float[] data2 =
          ic.getData(ip2.getPixels(), ip2.getWidth(), ip2.getHeight(), bounds, null);

      // Compute a correlation between the pixels
      double corr;
      double m1;
      double m2;
      if (imp.getBitDepth() == 32) {
        // 32-bit image
        final SimpleRegression regression = new SimpleRegression();
        m1 = 0;
        m2 = 0;
        for (int i = 0; i < data1.length; i++) {
          regression.addData(data1[i], data2[i]);
          m1 += data1[i];
          m2 += data2[i];
        }
        corr = regression.getR();
        m1 /= regression.getN();
        m2 /= regression.getN();
      } else {
        // 8/16-bit images
        final FastCorrelator fc = new FastCorrelator();
        for (int i = 0; i < data1.length; i++) {
          fc.add((int) data1[i], (int) data2[i]);
        }
        corr = fc.getCorrelation();
        m1 = (double) fc.getSumX() / fc.getN();
        m2 = (double) fc.getSumY() / fc.getN();
      }

      // If square:
      // Compute radial statistics
      double[] r1;
      double[] r2;
      if (bounds.width == bounds.height) {
        final double[][] stats =
            RadialStatisticsUtils.radialSumAndCount(bounds.width, data1, data2);
        r1 = stats[0];
        r2 = stats[1];
        final double[] c = stats[2];
        // Create a mean for each from the radial statistics or the entire pixel set.
        m1 = m2 = 0;
        double n = 0;
        for (int i = 0; i < r1.length; i++) {
          m1 += r1[i];
          m2 += r2[i];
          n += c[i];
          // Normalise radial sums to radial means
          r1[i] /= c[i];
          r2[i] /= c[i];
        }
        m1 /= n;
        m2 /= n;
      } else {
        r1 = r2 = null;
      }

      return new ComparisonResult(channel, background, intensity, m1, m2, corr, r1, r2);
    }

    /**
     * Gets the bounds around the given centre. The bounds are clipped to the image bounds.
     *
     * @param cx centre x
     * @param cy centre y
     * @param width the width
     * @param imp the imp
     * @return the bounds
     */
    private static Rectangle getBounds(int cx, int cy, int width, ImagePlus imp) {
      final int ox = cx - width;
      final int oy = cy - width;

      // Clip to the image
      return new Rectangle(ox, oy, 2 * width + 1, 2 * width + 1)
          .intersection(new Rectangle(imp.getWidth(), imp.getHeight()));
    }

    private static double getIntensity(double x, double y, double xsd, Rectangle bounds,
        double[] data, double background) {
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

      // Avoid divide by zero error
      return MathUtils.div0(intensity, sumWeights);
    }
  }

  /**
   * Initialise the spot fit tool. This is to allow support for calling within macro toolsets.
   */
  public static void addPluginTool() {
    // Add the tool
    Toolbar.addPlugInTool(SpotFitPluginTool.INSTANCE);
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

    SpotFitPluginTool.INSTANCE.showOptionsDialog();
  }
}
