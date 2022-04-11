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
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicReference;
import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.logging.Ticker;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.GUIFilterSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GuiSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.procedures.PrecisionResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.StandardResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.WidthResultProcedure;

/**
 * Filters PeakFit results that are stored in memory using various fit criteria.
 */
public class FilterResults implements PlugIn {
  private static final String TITLE = "Filter Results";
  private static AtomicReference<String> inputOptionRef = new AtomicReference<>("");

  private MemoryPeakResults results;

  private GUIFilterSettings.Builder filterSettings =
      GuiSettings.DefaultGUIFilterSettings.INSTANCE.toBuilder();

  // Used to pass data from analyseResults() to checkLimits()
  private float minDrift = Float.MAX_VALUE;
  private float maxDrift;
  private float minSignal = Float.MAX_VALUE;
  private float maxSignal;
  private float minSnr = Float.MAX_VALUE;
  private float maxSnr;
  private double minPrecision = Float.MAX_VALUE;
  private double maxPrecision;
  private double averageWidth;
  private float minWidth = Float.MAX_VALUE;
  private float maxWidth;

  private StandardResultProcedure sp;
  private PrecisionResultProcedure pp;
  private WidthResultProcedure wp;

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "There are no fitting results in memory");
      return;
    }

    String inputOption = inputOptionRef.get();

    // Show a dialog allowing the results set to be filtered
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addMessage("Select a dataset to filter");
    ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    inputOption = ResultsManager.getInputSource(gd);
    inputOptionRef.set(inputOption);

    results = ResultsManager.loadInputResults(inputOption, false, null, null);
    if (MemoryPeakResults.isEmpty(results)) {
      IJ.error(TITLE, "No results could be loaded");
      IJ.showStatus("");
      return;
    }

    if (!analyseResults()) {
      return;
    }

    if (!showDialog()) {
      return;
    }

    filterResults();
  }

  /**
   * Analyse the results and determine the range for each filter.
   */
  private boolean analyseResults() {
    IJ.showStatus("Analysing results ...");

    final ArrayList<String> error = new ArrayList<>();

    try {
      wp = new WidthResultProcedure(results, DistanceUnit.PIXEL);
      wp.getW();
      final float[] limits = MathUtils.limits(wp.wx);
      maxWidth = limits[1];
      minWidth = limits[0];
      averageWidth = MathUtils.sum(wp.wx) / wp.size();
    } catch (final DataException ex) {
      error.add(ex.getMessage());
      wp = null;
      maxWidth = minWidth = 0;
    }

    try {
      pp = new PrecisionResultProcedure(results);
      pp.getPrecision();

      final double[] limits = MathUtils.limits(pp.precisions);
      maxPrecision = limits[1];
      minPrecision = limits[0];
    } catch (final DataException ex) {
      error.add(ex.getMessage());
      pp = null;
      maxPrecision = minPrecision = 0;
    }

    try {
      sp = new StandardResultProcedure(results, DistanceUnit.PIXEL);
      sp.getXyr();

      // Re-use for convenience
      sp.intensity = new float[sp.x.length];
      sp.background = new float[sp.x.length];
      sp.z = new float[sp.x.length];

      for (int i = 0; i < sp.size(); i++) {
        if (i % 64 == 0) {
          IJ.showProgress(i, sp.size());
        }

        final PeakResult result = sp.peakResults[i];

        final float drift = getDrift(result, sp.x[i], sp.y[i]);
        if (maxDrift < drift) {
          maxDrift = drift;
        }
        if (minDrift > drift) {
          minDrift = drift;
        }

        final float signal = result.getIntensity();
        if (maxSignal < signal) {
          maxSignal = signal;
        }
        if (minSignal > signal) {
          minSignal = signal;
        }

        final float snr = getSnr(result);
        if (maxSnr < snr) {
          maxSnr = snr;
        }
        if (minSnr > snr) {
          minSnr = snr;
        }

        // for convenience
        sp.z[i] = drift;
        sp.intensity[i] = signal;
        sp.background[i] = snr;
      }
    } catch (final DataException ex) {
      error.add(ex.getMessage());
      sp = null;
    }

    if (error.size() == 3 || sp == null) {
      final StringBuilder sb = new StringBuilder("Unable to analyse the results:\n");
      for (final String s : error) {
        sb.append(s).append(".\n");
      }
      IJ.error(TITLE, sb.toString());
      return false;
    }

    ImageJUtils.finished();
    return true;
  }

  private static float getDrift(PeakResult result, float x, float y) {
    return Math.max(Math.abs(result.getOrigX() + 0.5f - x), Math.abs(result.getOrigY() + 0.5f - y));
  }

  private static float getSnr(PeakResult result) {
    if (result.getNoise() <= 0) {
      return 0;
    }
    return result.getIntensity() / result.getNoise();
  }

  /**
   * Check that none of the filter values are outside the limits.
   */
  private void checkLimits() {
    if (filterSettings.getMaxDrift() > maxDrift || filterSettings.getMaxDrift() < minDrift) {
      filterSettings.setMaxDrift(maxDrift);
    }

    if (filterSettings.getMinSignal() > maxSignal || filterSettings.getMinSignal() < minSignal) {
      filterSettings.setMinSignal(minSignal);
    }

    if (filterSettings.getMinSnr() > maxSnr || filterSettings.getMinSnr() < minSnr) {
      filterSettings.setMinSnr(minSnr);
    }

    if (filterSettings.getMaxPrecision() > maxPrecision
        || filterSettings.getMaxPrecision() < minPrecision) {
      filterSettings.setMaxPrecision(maxPrecision);
    }

    if (filterSettings.getMinWidth() > maxWidth || filterSettings.getMinWidth() < minWidth) {
      filterSettings.setMinWidth(minWidth);
    }

    if (filterSettings.getMaxWidth() > maxWidth || filterSettings.getMaxWidth() < minWidth) {
      filterSettings.setMaxWidth(maxWidth);
    }

    if (filterSettings.getMinWidth() > filterSettings.getMaxWidth()) {
      final float tmp = filterSettings.getMaxWidth();
      filterSettings.setMaxWidth(filterSettings.getMinWidth());
      filterSettings.setMinWidth(tmp);
    }
  }

  /**
   * Simple interface for filtering the coordinates.
   */
  private interface CoordFilter {
    /**
     * Test if the coordinates are allowed.
     *
     * @param frame the frame
     * @param x the x
     * @param y the y
     * @return true, if allowed
     */
    boolean match(int frame, float x, float y);
  }

  /**
   * Apply the filters to the data.
   */
  private void filterResults() {
    checkLimits();

    final MemoryPeakResults newResults = new MemoryPeakResults();
    newResults.copySettings(results);
    newResults.setName(results.getName() + " Filtered");

    // Initialise the mask
    CoordFilter maskFilter;
    final ImagePlus maskImp = WindowManager.getImage(filterSettings.getMaskTitle());
    if (maskImp != null) {
      final int maxx = maskImp.getWidth();
      final int maxy = maskImp.getHeight();
      final Rectangle bounds = results.getBounds();
      final double ox = bounds.getX();
      final double oy = bounds.getY();
      final double scaleX = bounds.getWidth() / maxx;
      final double scaleY = bounds.getHeight() / maxy;

      // Improve to allow stacks
      if (maskImp.getStackSize() > 1) {
        // 3D filter. The frame is mapped to the stack index.
        final ImageStack stack = maskImp.getImageStack();
        final ImageProcessor ip = stack.getProcessor(1);
        maskFilter = (t, x, y) -> {
          // Check stack index
          if (t >= 1 && t <= stack.size()) {
            int ix = (int) ((x - ox) / scaleX);
            int iy = (int) ((y - oy) / scaleY);
            if (ix >= 0 && ix < maxx && iy >= 0 && iy < maxy) {
              ip.setPixels(stack.getPixels(t));
              return ip.get(iy * maxx + ix) != 0;
            }
          }
          return false;
        };
      } else {
        // 2D filter.
        final ImageProcessor ip = maskImp.getProcessor();
        maskFilter = (t, x, y) -> {
          int ix = (int) ((x - ox) / scaleX);
          int iy = (int) ((y - oy) / scaleY);
          return (ix >= 0 && ix < maxx && iy >= 0 && iy < maxy && ip.get(iy * maxx + ix) != 0);
        };
      }
    } else {
      maskFilter = (t, x, y) -> true;
    }

    // Create the ticker with size+1 as we tick at the start of the loop
    final Ticker ticker = ImageJUtils.createTicker(results.size() + 1, 0);
    for (int i = 0, size = results.size(); i < size; i++) {
      ticker.tick();

      // sp will not be null

      // We stored the drift=z, intensity=signal, background=snr
      if (sp.z[i] > filterSettings.getMaxDrift()) {
        continue;
      }

      if (sp.intensity[i] < filterSettings.getMinSignal()) {
        continue;
      }

      if (sp.background[i] < filterSettings.getMinSnr()) {
        continue;
      }

      final PeakResult peakResult = sp.peakResults[i];

      // Check the coordinates are inside the mask
      if (!maskFilter.match(peakResult.getFrame(), sp.x[i], sp.y[i])) {
        continue;
      }

      if (pp != null && pp.precisions[i] > maxPrecision) {
        continue;
      }

      if (wp != null) {
        final float width = wp.wx[i];
        if (width < filterSettings.getMinWidth() || width > filterSettings.getMaxWidth()) {
          continue;
        }
      }

      // Passed all filters. Add to the results
      newResults.add(peakResult);
    }

    ImageJUtils.finished(TextUtils.pleural(newResults.size(), "Filtered localisation"));
    MemoryPeakResults.addResults(newResults);
  }

  private boolean showDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("filter-results"));

    filterSettings = SettingsManager.readGuiFilterSettings(0).toBuilder();

    checkLimits();

    gd.addSlider("Max_drift", minDrift, maxDrift, filterSettings.getMaxDrift());
    gd.addSlider("Min_Signal", minSignal, maxSignal, filterSettings.getMinSignal());
    gd.addSlider("Min_SNR", minSnr, maxSnr, filterSettings.getMinSnr());
    gd.addSlider("Min_Precision", minPrecision, maxPrecision, filterSettings.getMaxPrecision());

    // TODO - If calibrated present the widths in nm
    gd.addMessage("Average Width = " + IJ.d2s(averageWidth, 3));

    gd.addSlider("Min_Width", minWidth, maxWidth, filterSettings.getMinWidth());
    gd.addSlider("Max_Width", minWidth, maxWidth, filterSettings.getMaxWidth());

    // Get a list of potential mask images
    final String[] items = getImageList();
    gd.addChoice("Mask", items, filterSettings.getMaskTitle());

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    filterSettings.setMaxDrift((float) gd.getNextNumber());
    filterSettings.setMinSignal((float) gd.getNextNumber());
    filterSettings.setMinSnr((float) gd.getNextNumber());
    filterSettings.setMaxPrecision((float) gd.getNextNumber());
    filterSettings.setMinWidth((float) gd.getNextNumber());
    filterSettings.setMaxWidth((float) gd.getNextNumber());
    filterSettings.setMaskTitle(gd.getNextChoice());

    return SettingsManager.writeSettings(filterSettings.build());
  }

  /**
   * Build a list of all the image names.
   *
   * @return The list of images
   */
  public static String[] getImageList() {
    final ArrayList<String> newImageList = new ArrayList<>();
    newImageList.add("[None]");

    for (final int id : ImageJUtils.getIdList()) {
      final ImagePlus imp = WindowManager.getImage(id);
      if (imp == null || !imp.getProcessor().isBinary()) {
        continue;
      }
      newImageList.add(imp.getTitle());
    }

    return newImageList.toArray(new String[0]);
  }
}
