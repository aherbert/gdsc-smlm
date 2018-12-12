/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
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

import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.GUIFilterSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.model.MaskDistribution;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.procedures.PrecisionResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.StandardResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.WidthResultProcedure;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;

import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;

/**
 * Filters PeakFit results that are stored in memory using various fit criteria.
 */
public class FilterResults implements PlugIn {
  private static final String TITLE = "Filter Results";
  private static String inputOption = "";

  private MemoryPeakResults results;

  private GUIFilterSettings.Builder filterSettings =
      GUIProtosHelper.defaultGUIFilterSettings.toBuilder();

  // Used to pass data from analyseResults() to checkLimits()
  private float minDrift = Float.MAX_VALUE;
  private float maxDrift;
  private float minSignal = Float.MAX_VALUE;
  private float maxSignal;
  private float minSNR = Float.MAX_VALUE;
  private float maxSNR;
  private double minPrecision = Float.MAX_VALUE;
  private double maxPrecision;
  private double averageWidth;
  private float minWidth = Float.MAX_VALUE;
  private float maxWidth;

  private StandardResultProcedure sp;
  private PrecisionResultProcedure pp;
  private WidthResultProcedure wp;

  /** {@inheritDoc} */
  @Override
  public void run(String arg) {
    SMLMUsageTracker.recordPlugin(this.getClass(), arg);

    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "There are no fitting results in memory");
      return;
    }

    // Show a dialog allowing the results set to be filtered
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addMessage("Select a dataset to filter");
    ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    inputOption = ResultsManager.getInputSource(gd);
    results = ResultsManager.loadInputResults(inputOption, false, null, null);
    if (results == null || results.size() == 0) {
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
      sp.getXYR();

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
        if (maxSNR < snr) {
          maxSNR = snr;
        }
        if (minSNR > snr) {
          minSNR = snr;
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

    IJ.showProgress(1);
    IJ.showStatus("");
    return true;
  }

  private static float getDrift(PeakResult result, float x, float y) {
    return FastMath.max(Math.abs(result.getOrigX() + 0.5f - x),
        Math.abs(result.getOrigY() + 0.5f - y));
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

    if (filterSettings.getMinSnr() > maxSNR || filterSettings.getMinSnr() < minSNR) {
      filterSettings.setMinSnr(minSNR);
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
   * Apply the filters to the data.
   */
  private void filterResults() {
    checkLimits();

    final MemoryPeakResults newResults = new MemoryPeakResults();
    newResults.copySettings(results);
    newResults.setName(results.getName() + " Filtered");

    // Initialise the mask
    final ByteProcessor mask = getMask(filterSettings.getMaskTitle());
    MaskDistribution maskFilter = null;
    final float centreX = results.getBounds().width / 2.0f;
    final float centreY = results.getBounds().height / 2.0f;
    if (mask != null) {
      final double scaleX = (double) results.getBounds().width / mask.getWidth();
      final double scaleY = (double) results.getBounds().height / mask.getHeight();
      maskFilter = new MaskDistribution((byte[]) mask.getPixels(), mask.getWidth(),
          mask.getHeight(), 0, scaleX, scaleY);
    }

    for (int i = 0, size = results.size(); i < size; i++) {
      if (i % 64 == 0) {
        IJ.showProgress(i, size);
      }

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

      if (maskFilter != null) {
        // Check the coordinates are inside the mask
        final double[] xy = new double[] {sp.x[i] - centreX, sp.y[i] - centreY};
        if (!maskFilter.isWithinXy(xy)) {
          continue;
        }
      }

      if (pp != null) {
        if (pp.precisions[i] > maxPrecision) {
          continue;
        }
      }

      if (wp != null) {
        final float width = wp.wx[i];
        if (width < filterSettings.getMinWidth() || width > filterSettings.getMaxWidth()) {
          continue;
        }
      }

      // Passed all filters. Add to the results
      newResults.add(sp.peakResults[i]);
    }

    IJ.showProgress(1);
    IJ.showStatus(newResults.size() + " Filtered localisations");
    MemoryPeakResults.addResults(newResults);
  }

  private static ByteProcessor getMask(String maskTitle) {
    final ImagePlus imp = WindowManager.getImage(maskTitle);
    if (imp != null) {
      return (ByteProcessor) imp.getProcessor().convertToByte(false);
    }
    return null;
  }

  private boolean showDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);

    filterSettings = SettingsManager.readGUIFilterSettings(0).toBuilder();

    checkLimits();

    gd.addSlider("Max_drift", minDrift, maxDrift, filterSettings.getMaxDrift());
    gd.addSlider("Min_Signal", minSignal, maxSignal, filterSettings.getMinSignal());
    gd.addSlider("Min_SNR", minSNR, maxSNR, filterSettings.getMinSnr());
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
      if (imp == null) {
        continue;
      }
      if (!imp.getProcessor().isBinary()) {
        continue;
      }
      newImageList.add(imp.getTitle());
    }

    return newImageList.toArray(new String[0]);
  }
}
