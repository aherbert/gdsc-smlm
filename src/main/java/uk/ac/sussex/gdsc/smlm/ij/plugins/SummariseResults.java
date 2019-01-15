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

import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot.HistogramPlotBuilder;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.PrecisionMethod;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.SummariseResultsSettings;
import uk.ac.sussex.gdsc.smlm.data.config.UnitHelper;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.count.Counter;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.PrecisionResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.SNRResultProcedure;

import ij.IJ;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.text.TextPanel;
import ij.text.TextWindow;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

/**
 * Produces a summary table of the results that are stored in memory.
 */
public class SummariseResults implements PlugIn, MouseListener {
  private static final String TITLE = "Summarise Results";
  private static final String[] REMOVE_OUTLIERS = {"None", "1.5x IQR", "Top 2%"};

  private static TextWindow summary;
  private int histgramBins;
  private int removeOutliers;

  @Override
  public void run(String arg) {
    SMLMUsageTracker.recordPlugin(this.getClass(), arg);

    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "There are no fitting results in memory");
      clearSummaryTable();
      return;
    }

    createSummaryTable();
    final StringBuilder sb = new StringBuilder();
    int i = 0;
    final int nextFlush = 9;
    for (final MemoryPeakResults result : MemoryPeakResults.getAllResults()) {
      addSummary(sb, result);
      if (++i == nextFlush) {
        summary.append(sb.toString());
        sb.setLength(0);
      }
    }
    summary.append(sb.toString());
    summary.append("");
    summary.toFront();
  }

  /**
   * Remove all entries in the summary table if it showing.
   */
  public static void clearSummaryTable() {
    if (summary != null && summary.isShowing()) {
      summary.getTextPanel().clear();
    }
  }

  private void createSummaryTable() {
    if (summary == null || !summary.isShowing()) {
      final StringBuilder sb = new StringBuilder("Dataset\tN\tFrames\tTime\tMemory\tBounds");
      // Calibration
      sb.append("\tnm/pixel\tms/frame\tCamera\tDUnit\tIUnit\t3D\tPrecision Method");
      for (final String statName : new String[] {"Precision (nm)", "SNR"}) {
        sb.append("\tAv ").append(statName);
        sb.append("\tMedian ").append(statName);
        sb.append("\tMin ").append(statName);
        sb.append("\tMax ").append(statName);
      }
      summary = new TextWindow("Peak Results Summary", sb.toString(), "", 800, 300);
      summary.setVisible(true);

      summary.getTextPanel().addMouseListener(this);
    }
    // This could be optional but at current there is no dialog and it seems unnecessary
    clearSummaryTable();
  }

  private static int NO = -1;
  private static int UNKNOWN;
  private static int YES = 1;
  private int removeNullResults = UNKNOWN;

  private void addSummary(StringBuilder sb, MemoryPeakResults result) {
    final DescriptiveStatistics[] stats = new DescriptiveStatistics[2];
    for (int i = 0; i < stats.length; i++) {
      stats[i] = new DescriptiveStatistics();
    }

    if (result.hasNullResults()) {
      IJ.log("Null results in dataset: " + result.getName());
      if (removeNullResults == UNKNOWN) {
        final GenericDialog gd = new GenericDialog(TITLE);
        gd.addMessage("There are invalid results in memory.\n \nClean these results?");
        gd.enableYesNoCancel();
        gd.hideCancelButton();
        gd.showDialog();
        removeNullResults = (gd.wasOKed()) ? YES : NO;
      }
      if (removeNullResults == NO) {
        result = result.copy();
      }
      result.removeNullResults();
    }

    final CalibrationReader calibration = result.getCalibrationReader();

    PrecisionMethod precisionMethod = PrecisionMethod.PRECISION_METHOD_NA;
    boolean stored = false;
    final int size = result.size();
    if (size > 0) {
      // Precision
      try {
        final PrecisionResultProcedure p = new PrecisionResultProcedure(result);
        // Use stored precision if possible
        stored = result.hasPrecision();
        precisionMethod = p.getPrecision(stored);
        for (final double v : p.precisions) {
          stats[0].addValue(v);
        }
      } catch (final DataException ex) {
        // Ignore
      }

      // SNR
      try {
        final SNRResultProcedure p = new SNRResultProcedure(result);
        p.getSNR();
        for (final double v : p.snr) {
          stats[1].addValue(v);
        }
      } catch (final DataException ex) {
        // Ignore
      }
    }

    sb.append(result.getName());
    int maxT = 0;
    if (result.size() == 0) {
      sb.append("\t0\t0");
    } else {
      sb.append('\t').append(result.size());
      maxT = result.getMaxFrame();
      sb.append('\t').append(maxT);
    }
    if (calibration != null && calibration.hasExposureTime()) {
      sb.append('\t')
          .append(TextUtils.millisToString((long) Math.ceil(maxT * calibration.getExposureTime())));
    } else {
      sb.append("\t-");
    }
    if (size > 0) {
      final boolean includeDeviations = result.hasDeviations();
      final long memorySize = MemoryPeakResults.estimateMemorySize(size, includeDeviations);
      final String memory = TextUtils.bytesToString(memorySize);
      sb.append('\t').append(memory);
    } else {
      sb.append("\t-");
    }
    final Rectangle bounds = result.getBounds(true);
    sb.append(String.format("\t%d,%d,%d,%d", bounds.x, bounds.y, bounds.x + bounds.width,
        bounds.y + bounds.height));
    if (calibration != null) {
      //@formatter:off
      sb.append('\t').append(calibration.hasNmPerPixel() ? MathUtils.rounded(calibration.getNmPerPixel()) : '-');
      sb.append('\t').append(calibration.hasExposureTime() ? MathUtils.rounded(calibration.getExposureTime()) : '-');

      if (calibration.hasCameraType())
      {
        sb.append('\t').append(CalibrationProtosHelper.getName(calibration.getCameraType()));
        if (calibration.isCCDCamera())
        {
          sb.append(" bias=").append(calibration.getBias());
          sb.append(" gain=").append(calibration.getCountPerPhoton());
        }
      } else {
        sb.append("\t-");
      }

      sb.append('\t').append(calibration.hasDistanceUnit() ? UnitHelper.getShortName(calibration.getDistanceUnit()) : '-');
      sb.append('\t').append(calibration.hasIntensityUnit() ? UnitHelper.getShortName(calibration.getIntensityUnit()) : '-');
      //@formatter:on
    } else {
      sb.append("\t\t\t\t\t");
    }

    if (result.is3D()) {
      sb.append("\tY");
    } else {
      sb.append("\tN");
    }

    sb.append("\t").append(FitProtosHelper.getName(precisionMethod));
    if (stored) {
      sb.append(" (Stored)");
    }
    for (int i = 0; i < stats.length; i++) {
      if (Double.isNaN(stats[i].getMean())) {
        sb.append("\t-\t-\t-\t-");
      } else {
        sb.append('\t').append(IJ.d2s(stats[i].getMean(), 3));
        sb.append('\t').append(IJ.d2s(stats[i].getPercentile(50), 3));
        sb.append('\t').append(IJ.d2s(stats[i].getMin(), 3));
        sb.append('\t').append(IJ.d2s(stats[i].getMax(), 3));
      }
    }
    sb.append("\n");
  }

  @Override
  public void mouseClicked(MouseEvent event) {
    if (event.getClickCount() > 1) {
      showStatistics();
      event.consume();
    }
  }

  private void showStatistics() {
    final TextPanel textPanel = summary.getTextPanel();
    final int selectedIndex = textPanel.getSelectionStart();
    if (selectedIndex < 0 || selectedIndex >= textPanel.getLineCount()) {
      return;
    }
    final String line = textPanel.getLine(selectedIndex);
    final int endIndex = line.indexOf('\t');
    if (endIndex == -1) {
      return;
    }
    final String name = line.substring(0, endIndex);
    final MemoryPeakResults result = MemoryPeakResults.getResults(name);
    if (result == null) {
      return;
    }

    // Do this is a thread so the click-event does not block
    new Thread(new Runnable() {
      @Override
      public void run() {
        showStatistics(result);
      }
    }).run();
  }

  private void showStatistics(MemoryPeakResults result) {
    final SummariseResultsSettings.Builder settings =
        SettingsManager.readSummariseResultsSettings(0).toBuilder();
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addMessage("Show histograms of the results properties (if available)");
    gd.addCheckbox("Plot_background", settings.getPlotBackground());
    gd.addCheckbox("Plot_signal", settings.getPlotSignal());
    gd.addCheckbox("Plot_x", settings.getPlotX());
    gd.addCheckbox("Plot_y", settings.getPlotY());
    gd.addCheckbox("Plot_z", settings.getPlotZ());
    gd.addCheckbox("Plot_noise", settings.getPlotNoise());
    gd.addCheckbox("Plot_SNR", settings.getPlotSnr());
    gd.addCheckbox("Plot_precision", settings.getPlotPrecision());
    gd.addNumericField("Histgram_bins", settings.getHistgramBins(), 0);
    gd.addChoice("Remove_outliers", REMOVE_OUTLIERS, settings.getRemoveOutliers());
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    settings.setPlotBackground(gd.getNextBoolean());
    settings.setPlotSignal(gd.getNextBoolean());
    settings.setPlotX(gd.getNextBoolean());
    settings.setPlotY(gd.getNextBoolean());
    settings.setPlotZ(gd.getNextBoolean());
    settings.setPlotNoise(gd.getNextBoolean());
    settings.setPlotSnr(gd.getNextBoolean());
    settings.setPlotPrecision(gd.getNextBoolean());
    histgramBins = Math.max(0, (int) gd.getNextNumber());
    removeOutliers = gd.getNextChoiceIndex();
    settings.setHistgramBins(histgramBins);
    settings.setRemoveOutliers(removeOutliers);
    SettingsManager.writeSettings(settings);

    final WindowOrganiser wo = new WindowOrganiser();

    if (settings.getPlotBackground()) {
      plot(wo, "Background", result, PeakResult.BACKGROUND);
    }
    if (settings.getPlotSignal()) {
      plot(wo, "Signal", result, PeakResult.INTENSITY);
    }
    if (settings.getPlotX()) {
      plot(wo, "X", result, PeakResult.X);
    }
    if (settings.getPlotY()) {
      plot(wo, "Y", result, PeakResult.Y);
    }
    if (settings.getPlotZ()) {
      plot(wo, "Z", result, PeakResult.Z);
    }
    if ((settings.getPlotNoise() || settings.getPlotSnr()) && result.hasNoise()) {
      if (settings.getPlotSnr()) {
        try {
          plot(wo, "SNR", new SNRResultProcedure(result).getSNR());
        } catch (final DataException ex) {
          // Ignore
        }
      }
      if (settings.getPlotNoise()) {
        final Counter counter = new Counter();
        final float[] noise = new float[result.size()];
        result.forEach(new PeakResultProcedure() {
          @Override
          public void execute(PeakResult peakResult) {
            final int i = counter.getAndIncrement();
            noise[i] = peakResult.getNoise();
          }
        });
        plot(wo, "Noise", noise);
      }
    }
    if (settings.getPlotPrecision()) {
      // Precision
      try {
        final PrecisionResultProcedure p = new PrecisionResultProcedure(result);
        // Use stored precision if possible
        final boolean stored = result.hasPrecision();
        final PrecisionMethod precisionMethod = p.getPrecision(stored);
        String name = FitProtosHelper.getName(precisionMethod);
        if (stored) {
          name += " (Stored)";
        }
        plot(wo, "Precision: " + name, StoredDataStatistics.create(p.precisions));
      } catch (final DataException ex) {
        // Ignore
      }
    }

    wo.tile();
  }

  private void plot(WindowOrganiser wo, String title, MemoryPeakResults result, final int index) {
    final StoredDataStatistics data = new StoredDataStatistics(result.size());
    result.forEach(new PeakResultProcedure() {
      @Override
      public void execute(PeakResult peakResult) {
        data.add(peakResult.getParameter(index));
      }
    });
    plot(wo, title, data);
  }

  private void plot(WindowOrganiser wo, String title, float[] data) {
    plot(wo, title, StoredDataStatistics.create(data));
  }

  private void plot(WindowOrganiser wo, String title, StoredDataStatistics data) {
    new HistogramPlotBuilder(TITLE, data, title).setRemoveOutliersOption(removeOutliers)
        .setNumberOfBins(histgramBins).show(wo);
  }

  @Override
  public void mousePressed(MouseEvent event) {
    // Ignore
  }

  @Override
  public void mouseReleased(MouseEvent event) {
    // Ignore
  }

  @Override
  public void mouseEntered(MouseEvent event) {
    // Ignore
  }

  @Override
  public void mouseExited(MouseEvent event) {
    // Ignore
  }
}
