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
import ij.WindowManager;
import ij.gui.ImageWindow;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.FloatPolygon;
import ij.process.LUT;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Function;
import java.util.function.ToDoubleFunction;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.OffsetLineRoi;
import uk.ac.sussex.gdsc.core.ij.gui.OffsetPointRoi;
import uk.ac.sussex.gdsc.core.ij.gui.OffsetPolygonRoi;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.Trace;
import uk.ac.sussex.gdsc.smlm.results.TraceManager;

/**
 * Compares the coordinates in sets of traced results and computes the match statistics.
 */
public class DrawClusters implements PlugIn {
  private static final String TITLE = "Draw Clusters";
  private static final String[] COLOURS =
      {"None", "ID", "Time", "Size", "Length", "MSD", "Mean/Frame", "Category"};

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> INSTANCE = new AtomicReference<>(new Settings());

    String inputOption;
    String title;
    int imageSize;
    boolean expandToSingles;
    int minSize;
    int maxSize;
    boolean drawLines;
    int colour;
    boolean splineFit;
    boolean useStackPosition;
    int lut;
    boolean invertLut;
    float lineWidth;

    Settings() {
      // Set defaults
      inputOption = "";
      title = "";
      imageSize = 20;
      minSize = 2;
      drawLines = true;
    }

    Settings(Settings source) {
      inputOption = source.inputOption;
      title = source.title;
      imageSize = source.imageSize;
      expandToSingles = source.expandToSingles;
      minSize = source.minSize;
      maxSize = source.maxSize;
      drawLines = source.drawLines;
      colour = source.colour;
      splineFit = source.splineFit;
      useStackPosition = source.useStackPosition;
      lut = source.lut;
      invertLut = source.invertLut;
      lineWidth = source.lineWidth;
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
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "No localisations in memory");
      return;
    }

    if (!showDialog()) {
      return;
    }

    // Load the results
    final MemoryPeakResults results =
        ResultsManager.loadInputResults(settings.inputOption, false, DistanceUnit.PIXEL);
    if (MemoryPeakResults.isEmpty(results)) {
      IJ.error(TITLE, "No results could be loaded");
      return;
    }

    // Get the traces
    final Trace[] traces = TraceManager.convert(results);
    if (traces == null || traces.length == 0) {
      IJ.error(TITLE, "No traces could be loaded");
      return;
    }

    // Filter traces to a min size
    int maxFrame = 0;
    int count = 0;
    final int myMaxSize =
        (settings.maxSize < settings.minSize) ? Integer.MAX_VALUE : settings.maxSize;
    for (final Trace trace : traces) {
      if (settings.expandToSingles) {
        trace.expandToSingles();
      }
      if (trace.size() >= settings.minSize && trace.size() <= myMaxSize) {
        traces[count++] = trace;
        trace.sort();
        final int end = trace.getTail().getFrame();
        maxFrame = maxFrame < end ? end : maxFrame;
      }
    }

    if (count == 0) {
      IJ.error(TITLE, "No traces achieved the size limits");
      return;
    }

    // There are only traces if each trace has only 1 localisation per frame. Otherwise
    // this is some type of clustering result. Only allow drawing lines from traces.
    final boolean myDrawLines = settings.drawLines && myMaxSize > 1 && isTraced(traces);

    final ToDoubleFunction<PeakResult> perLocalisationColour =
        getPerLocalisationColour(settings.colour);
    // Spline-fit is only used the entire trace is added. If split up into per-localisation
    // lines then the spline is not used.
    final boolean doSplineFit = settings.splineFit && perLocalisationColour == null;

    final String msg =
        String.format(TITLE + ": %d / %s (%s)", count, TextUtils.pleural(traces.length, "trace"),
            TextUtils.pleural(results.size(), "localisation"));
    IJ.showStatus(msg);

    final Rectangle bounds = results.getBounds(true);
    ImagePlus imp = WindowManager.getImage(settings.title);
    boolean isUseStackPosition = settings.useStackPosition;
    if (imp == null) {
      // Create a default image using 100 pixels as the longest edge
      final double maxD = (bounds.width > bounds.height) ? bounds.width : bounds.height;
      int width;
      int height;
      if (maxD == 0) {
        // Note that imageSize can be zero (for auto sizing)
        width = height = (settings.imageSize == 0) ? 20 : settings.imageSize;
      } else if (settings.imageSize == 0) {
        // Note that imageSize can be zero (for auto sizing)
        width = bounds.width;
        height = bounds.height;
      } else {
        width = (int) (settings.imageSize * (bounds.width / maxD));
        height = (int) (settings.imageSize * (bounds.height / maxD));
      }
      final ByteProcessor bp = new ByteProcessor(width, height);
      if (isUseStackPosition) {
        final ImageStack stack = new ImageStack(width, height, maxFrame);
        for (int i = 1; i <= maxFrame; i++) {
          stack.setPixels(bp.getPixels(), i); // Do not duplicate pixel as the image is empty
        }
        imp = ImageJUtils.display(TITLE, stack);
      } else {
        imp = ImageJUtils.display(TITLE, bp);
      }

      // Enlarge
      final ImageWindow iw = imp.getWindow();
      for (int i = 9; i-- > 0 && iw.getWidth() < 500 && iw.getHeight() < 500;) {
        iw.getCanvas().zoomIn(imp.getWidth() / 2, imp.getHeight() / 2);
      }

      // Check if the image has enough frames for all the traces
    } else if (maxFrame > imp.getNFrames()) {
      isUseStackPosition = false;
    }

    final float xScale = (float) (imp.getWidth() / bounds.getWidth());
    final float yScale = (float) (imp.getHeight() / bounds.getHeight());

    // Create ROIs and store data to sort them
    final Roi[] rois = new Roi[count];
    final int[][] frames = (isUseStackPosition) ? new int[count][] : null;
    final double[] values = new double[count];
    for (int i = 0; i < count; i++) {
      final Trace trace = traces[i];
      final int npoints = trace.size();
      final float[] xpoints = new float[npoints];
      final float[] ypoints = new float[npoints];
      if (frames != null) {
        frames[i] = new int[npoints];
      }
      for (int k = 0; k < npoints; k++) {
        final PeakResult result = trace.get(k);
        xpoints[k] = (result.getXPosition() - bounds.x) * xScale;
        ypoints[k] = (result.getYPosition() - bounds.y) * yScale;
        if (frames != null) {
          frames[i][k] = result.getFrame();
        }
      }
      Roi roi;
      if (myDrawLines) {
        roi = new OffsetPolygonRoi(xpoints, ypoints, npoints, Roi.POLYLINE);
        if (doSplineFit) {
          ((PolygonRoi) roi).fitSpline();
        }
      } else {
        final OffsetPointRoi r = new OffsetPointRoi(xpoints, ypoints, npoints);
        r.setPointType(2);
        r.setShowLabels(false);
        roi = r;
      }

      rois[i] = roi;
      switch (settings.colour) {
        case 1: // ID (positive only)
          values[i] = Math.max(0, trace.getId());
          break;
        case 2: // Time
          values[i] = trace.getHead().getFrame();
          break;
        case 3: // Cluster size
          values[i] = trace.size();
          break;
        case 4: // Track length
          values[i] = roi.getLength();
          break;
        case 5: // Mean Square Displacement
          values[i] = trace.getMsd();
          break;
        case 6: // Mean / Frame
          values[i] = trace.getMeanDistance();
          break;
        case 7: // Category
          values[i] = maxCategory(trace);
          break;
        case 0: // No sort
        default:
          values[i] = i;
          break;
      }
    }

    // Draw the traces as ROIs on an overlay
    final Overlay o = new Overlay();
    final LUT lut = settings.invertLut ? LutHelper.createLut(settings.lut).createInvertedLut()
        : LutHelper.createLut(settings.lut);
    // Colour assumes the values are a linear scale from [0, max]
    final double max = MathUtils.max(values);
    final double scale = 255.0 / max;
    final Function<PeakResult,
        Color> localisationColour = perLocalisationColour == null ? null
            : result -> LutHelper.getColour(lut,
                (int) Math.round(perLocalisationColour.applyAsDouble(result) * scale));

    final LocalList<Roi> localisations = new LocalList<>(11);
    if (frames == null) {
      // Add the tracks as a single overlay
      for (int i = 0; i < count; i++) {
        final PolygonRoi roi = (PolygonRoi) rois[i];
        // Colour per localisation
        if (localisationColour == null) {
          roi.setStrokeColor(new Color(lut.getRGB((int) Math.round(values[i] * scale))));
          // roi.setStrokeWidth(settings.lineWidth);
          roi.updateWideLine(settings.lineWidth);
          o.add(roi);
        } else {
          toLocalisations(roi, traces[i], localisationColour, myDrawLines, settings.lineWidth,
              localisations);
          localisations.forEach(o::add);
        }
        localisations.clear();
      }
    } else {
      // Add the tracks on the frames containing the results
      final boolean isHyperStack = imp.isDisplayedHyperStack();
      for (int i = 0; i < count; i++) {
        Color c = null;
        final Trace trace = traces[i];
        final PolygonRoi roi = (PolygonRoi) rois[i];
        // Colour per localisation
        if (localisationColour == null) {
          c = LutHelper.getColour(lut, (int) Math.round(values[i] * scale));
          roi.setFillColor(c);
          roi.setStrokeColor(c);
          roi.updateWideLine(settings.lineWidth);
          localisations.push(roi);
        } else {
          toLocalisations(roi, trace, localisationColour, myDrawLines, settings.lineWidth,
              localisations);
        }
        final FloatPolygon fp = roi.getNonSplineFloatPolygon();
        // For each frame in the track, add a point ROI for the current position.
        // If this is a trace then add the entire the ROI(s) track to the frame
        // for convenience. The current position will be highlighted with a circle.
        for (int j = 0; j < frames[i].length; j++) {
          final int frame = frames[i][j];
          // It only makes sense to add all localisations to a single frame if this is
          // a trace. Otherwise add only the relevant point per frame.
          final PointRoi pointRoi = new OffsetPointRoi(fp.xpoints[j], fp.ypoints[j]);
          if (myDrawLines) {
            localisations.forEach(r -> addToOverlay(o, (Roi) r.clone(), isHyperStack, frame));
            pointRoi.setPointType(3);
          } else {
            pointRoi.setPointType(2);
          }
          if (localisationColour != null) {
            c = localisationColour.apply(trace.get(j));
          }
          pointRoi.setFillColor(c);
          pointRoi.setStrokeColor(c);
          addToOverlay(o, pointRoi, isHyperStack, frame);
        }
        localisations.clear();
      }
    }
    imp.setOverlay(o);

    IJ.showStatus(msg);
  }

  /**
   * Checks if is the entire dataset is traced. This requires each trace to have only 1 localisation
   * per frame.
   *
   * @param traces the traces
   * @return true if traced
   */
  private static boolean isTraced(Trace[] traces) {
    final IntOpenHashSet set = new IntOpenHashSet();
    for (final Trace trace : traces) {
      set.clear();
      for (int k = 0; k < trace.size(); k++) {
        if (!set.add(trace.get(k).getFrame())) {
          // Already seen
          return false;
        }
      }
    }
    return true;
  }

  /**
   * Get the maximum category in the trace.
   *
   * @param trace the trace
   * @return the maximum catergory
   */
  private static int maxCategory(Trace trace) {
    int max = 0;
    for (int k = 0; k < trace.size(); k++) {
      final PeakResult result = trace.get(k);
      max = Math.max(max, result.getCategory());
    }
    return max;
  }

  /**
   * Create a function to generate colour is per localisation.
   *
   * @param colour the colour
   * @return the function (or null)
   */
  private static ToDoubleFunction<PeakResult> getPerLocalisationColour(int colour) {
    // Only per category colour is currently supported
    if (colour == 7) {
      return PeakResult::getCategory;
    }
    return null;
  }

  /**
   * Convert the ROI points to a set of per-localisation ROIs. Optionally draw lines; otherwise draw
   * points. Get colour from the trace results using the colour function.
   *
   * @param roi the roi
   * @param trace the trace
   * @param localisationColour the localisation colour
   * @param drawLines the draw lines
   * @param lut the lut
   * @param scale the scale
   * @param lineWidth the line width
   * @param localisations the per-localisation ROIs(output)
   */
  private static void toLocalisations(PolygonRoi roi, Trace trace,
      Function<PeakResult, Color> localisationColour, boolean drawLines, float lineWidth,
      LocalList<Roi> localisations) {
    // Extract all the points from the ROI FloatPolygon
    final FloatPolygon fp = roi.getNonSplineFloatPolygon();
    // If draw traces then draw lines; otherwise draw points.
    // Get colour from the trace results.
    final int npoints = trace.size();
    if (drawLines) {
      for (int k = 1; k < npoints; k++) {
        final PeakResult result = trace.get(k - 1);
        final Color c = localisationColour.apply(result);
        final OffsetLineRoi line =
            new OffsetLineRoi(fp.xpoints[k - 1], fp.ypoints[k - 1], fp.xpoints[k], fp.ypoints[k]);
        line.setStrokeColor(c);
        localisations.add(line);
      }
    } else {
      for (int k = 0; k < npoints; k++) {
        final PeakResult result = trace.get(k);
        final Color c = localisationColour.apply(result);
        final OffsetPointRoi pointRoi = new OffsetPointRoi(fp.xpoints[k], fp.ypoints[k]);
        pointRoi.setPointType(2);
        pointRoi.setStrokeColor(c);
        pointRoi.updateWideLine(lineWidth);
        localisations.add(pointRoi);
      }
    }
  }

  private boolean showDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

    final ArrayList<String> titles = new ArrayList<>(WindowManager.getImageCount());
    titles.add("[None]");
    final int[] idList = WindowManager.getIDList();
    if (idList != null) {
      for (final int id : idList) {
        final ImagePlus imp = WindowManager.getImage(id);
        if (imp != null) {
          titles.add(imp.getTitle());
        }
      }
    }

    settings = Settings.load();
    gd.addMessage("Draw the clusters on an image");
    ResultsManager.addInput(gd, "Input", settings.inputOption, InputSource.MEMORY_CLUSTERED);
    gd.addChoice("Image", titles.toArray(new String[0]), settings.title);
    gd.addNumericField("Image_size", settings.imageSize, 0);
    gd.addCheckbox("Expand_to_singles", settings.expandToSingles);
    gd.addSlider("Min_size", 1, 15, settings.minSize);
    gd.addSlider("Max_size", 0, 20, settings.maxSize);
    gd.addCheckbox("Traces (draw lines)", settings.drawLines);
    gd.addMessage("Spline fit cannot be used when colouring per localisation (e.g. catagory)");
    gd.addCheckbox("Spline_fit (traces only)", settings.splineFit);
    gd.addCheckbox("Use_stack_position", settings.useStackPosition);
    gd.addChoice("Colour", COLOURS, COLOURS[settings.colour]);
    gd.addChoice("LUT", LutHelper.getLutNames(), settings.lut);
    gd.addCheckbox("Invert_LUT", settings.invertLut);
    gd.addSlider("Line_width", 0, 0.5, settings.lineWidth);

    gd.addHelp(HelpUrls.getUrl("draw-clusters"));
    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.inputOption = ResultsManager.getInputSource(gd);
    settings.title = gd.getNextChoice();
    settings.imageSize = (int) Math.abs(gd.getNextNumber());
    settings.expandToSingles = gd.getNextBoolean();
    settings.minSize = (int) Math.abs(gd.getNextNumber());
    settings.maxSize = (int) Math.abs(gd.getNextNumber());
    settings.drawLines = gd.getNextBoolean();
    settings.splineFit = gd.getNextBoolean();
    settings.useStackPosition = gd.getNextBoolean();
    settings.colour = gd.getNextChoiceIndex();
    settings.lut = gd.getNextChoiceIndex();
    settings.invertLut = gd.getNextBoolean();
    settings.lineWidth = (float) gd.getNextNumber();
    settings.save();

    return true;
  }

  private static void addToOverlay(Overlay overlay, Roi roi, boolean isHyperStack, int frame) {
    if (isHyperStack) {
      roi.setPosition(0, 0, frame);
    } else {
      roi.setPosition(frame);
    }
    overlay.add(roi);
  }
}
