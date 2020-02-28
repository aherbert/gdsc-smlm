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
import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicReference;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.SortUtils;
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
  private static final String[] sorts =
      new String[] {"None", "ID", "Time", "Size", "Length", "MSD", "Mean/Frame"};

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    String inputOption;
    String title;
    int imageSize;
    boolean expandToSingles;
    int minSize;
    int maxSize;
    boolean drawLines;
    int sort;
    boolean splineFit;
    boolean useStackPosition;
    int lut;
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
      sort = source.sort;
      splineFit = source.splineFit;
      useStackPosition = source.useStackPosition;
      lut = source.lut;
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
    final boolean myDrawLines = myMaxSize > 1 && settings.drawLines;
    for (int i = 0; i < traces.length; i++) {
      if (settings.expandToSingles) {
        traces[i].expandToSingles();
      }
      if (traces[i].size() >= settings.minSize && traces[i].size() <= myMaxSize) {
        traces[count++] = traces[i];
        traces[i].sort();
        if (maxFrame < traces[i].getTail().getFrame()) {
          maxFrame = traces[i].getTail().getFrame();
        }
      }
    }

    if (count == 0) {
      IJ.error(TITLE, "No traces achieved the size limits");
      return;
    }

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
        width = (int) (settings.imageSize * bounds.width / maxD);
        height = (int) (settings.imageSize * bounds.height / maxD);
      }
      final ByteProcessor bp = new ByteProcessor(width, height);
      if (isUseStackPosition) {
        final ImageStack stack = new ImageStack(width, height, maxFrame);
        for (int i = 1; i <= maxFrame; i++) {
          stack.setPixels(bp.getPixels(), i); // Do not clone as the image is empty
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
    final int[] indices = SimpleArrayUtils.natural(count);
    final double[] values = new double[count];
    for (int i = 0; i < count; i++) {
      final Trace trace = traces[i];
      final int npoints = trace.size();
      final float[] xpoints = new float[npoints];
      final float[] ypoints = new float[npoints];
      int ii = 0;
      if (frames != null) {
        frames[i] = new int[npoints];
      }
      for (int k = 0; k < trace.size(); k++) {
        final PeakResult result = trace.get(k);
        xpoints[ii] = (result.getXPosition() - bounds.x) * xScale;
        ypoints[ii] = (result.getYPosition() - bounds.y) * yScale;
        if (frames != null) {
          frames[i][ii] = result.getFrame();
        }
        ii++;
      }
      Roi roi;
      if (myDrawLines) {
        roi = new PolygonRoi(xpoints, ypoints, npoints, Roi.POLYLINE);
        if (settings.splineFit) {
          ((PolygonRoi) roi).fitSpline();
        }
      } else {
        roi = new PointRoi(xpoints, ypoints, npoints);
        ((PointRoi) roi).setShowLabels(false);
      }

      rois[i] = roi;
      switch (settings.sort) {
        case 1: // Sort by ID
          values[i] = traces[i].getId();
          break;
        case 2: // Sort by time
          values[i] = traces[i].getHead().getFrame();
          break;
        case 3: // Sort by size descending
          values[i] = -traces[i].size();
          break;
        case 4: // Sort by length descending
          values[i] = -roi.getLength();
          break;
        case 5: // Mean Square Displacement
          values[i] = -traces[i].getMsd();
          break;
        case 6: // Mean / Frame
          values[i] = -traces[i].getMeanDistance();
          break;
        case 0: // No sort
        default:
          break;
      }
    }

    if (settings.sort > 0) {
      SortUtils.sortIndices(indices, values, true);
    }

    // Draw the traces as ROIs on an overlay
    final Overlay o = new Overlay();
    final LUT lut = LutHelper.createLut(settings.lut);
    final double scale = 256.0 / count;
    if (frames != null) {
      // Add the tracks on the frames containing the results
      final boolean isHyperStack = imp.isDisplayedHyperStack();
      for (int i = 0; i < count; i++) {
        final int index = indices[i];
        final Color c = LutHelper.getColour(lut, (int) (i * scale));
        final PolygonRoi roi = (PolygonRoi) rois[index];
        roi.setFillColor(c);
        roi.setStrokeColor(c);
        //roi.setStrokeWidth(settings.lineWidth);
        roi.updateWideLine(settings.lineWidth);
        final FloatPolygon fp = roi.getNonSplineFloatPolygon();
        // For each frame in the track, add the ROI track and a point ROI for the current position
        for (int j = 0; j < frames[index].length; j++) {
          addToOverlay(o, (Roi) roi.clone(), isHyperStack, frames[index][j]);
          final PointRoi pointRoi = new PointRoi(fp.xpoints[j], fp.ypoints[j]);
          pointRoi.setPointType(3);
          pointRoi.setFillColor(c);
          pointRoi.setStrokeColor(Color.black);
          addToOverlay(o, pointRoi, isHyperStack, frames[index][j]);
        }
      }
    } else {
      // Add the tracks as a single overlay
      for (int i = 0; i < count; i++) {
        final Roi roi = rois[indices[i]];
        roi.setStrokeColor(new Color(lut.getRGB((int) (i * scale))));
        //roi.setStrokeWidth(settings.lineWidth);
        roi.updateWideLine(settings.lineWidth);
        o.add(roi);
      }
    }
    imp.setOverlay(o);

    IJ.showStatus(msg);
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
    gd.addChoice("Sort", sorts, sorts[settings.sort]);
    gd.addCheckbox("Spline_fit (traces only)", settings.splineFit);
    gd.addCheckbox("Use_stack_position", settings.useStackPosition);
    gd.addChoice("LUT", LutHelper.getLutNames(), settings.lut);
    gd.addSlider("Line_width", 0, 0.5, settings.lineWidth);

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
    settings.sort = gd.getNextChoiceIndex();
    settings.splineFit = gd.getNextBoolean();
    settings.useStackPosition = gd.getNextBoolean();
    settings.lut = gd.getNextChoiceIndex();
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
