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

import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot.HistogramPlotBuilder;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.NonBlockingExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.SoftLock;
import uk.ac.sussex.gdsc.core.utils.SortUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.StoredData;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.TimeUnit;
import uk.ac.sussex.gdsc.smlm.fitting.JumpDistanceAnalysis;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.procedures.PrecisionResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.sort.IdFramePeakResultComparator;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;

import ij.IJ;
import ij.gui.Plot;
import ij.plugin.PlugIn;

import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import java.awt.Color;
import java.util.Arrays;

/**
 * Analyses the track lengths of traced data.
 */
public class TraceLengthAnalysis implements PlugIn {
  private static final String TITLE = "Trace Length Analysis";
  private static final float WIDTH = 0.3f;

  private static String inputOption = "";
  private static double msdThreshold;
  private static boolean normalise;

  private TypeConverter<DistanceUnit> distanceConverter;
  private TypeConverter<TimeUnit> timeConverter;
  private final SoftLock lock = new SoftLock();
  private double lastMsdThreshold = -1;
  private boolean lastNormalise;
  private int lastIndex;
  private double error;
  private double[] msds; // MSD of trace
  private int[] lengths; // Length of trace
  private int[] ids; // trace id
  private double minX;
  private double maxX;
  private int[] h1;
  private int[] h2;
  private float[] x1;
  private float[] x2;
  private float[] y1;
  private float[] y2;

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
    MemoryPeakResults results = ResultsManager.loadInputResults(inputOption, false, null, null);
    if (results == null || results.size() == 0) {
      IJ.error(TITLE, "No results could be loaded");
      return;
    }

    try {
      distanceConverter = results.getDistanceConverter(DistanceUnit.UM);
      timeConverter = results.getTimeConverter(TimeUnit.SECOND);
    } catch (final Exception ex) {
      IJ.error(TITLE, "Cannot convert units to um or seconds: " + ex.getMessage());
      return;
    }

    // Get the localisation error (4s^2) in raw units^2
    double precision = 0;
    try {
      final PrecisionResultProcedure p = new PrecisionResultProcedure(results);
      p.getPrecision();

      // Precision in nm using the median
      precision = new Percentile().evaluate(p.precisions, 50);
      // Convert from nm to um to raw units
      final double rawPrecision = distanceConverter.convertBack(precision / 1e3);
      // Get the localisation error (4s^2) in units^2
      error = 4 * rawPrecision * rawPrecision;
    } catch (final Exception ex) {
      ImageJUtils.log(TITLE + " - Unable to compute precision: " + ex.getMessage());
    }

    // Analyse the track lengths
    results = results.copy();
    results.sort(IdFramePeakResultComparator.INSTANCE);
    // Ensure the first result triggers an id change
    lastid = results.getFirst().getId() - 1;
    results.forEach(this::processTrackLength);
    store(); // For the final track
    msds = msdList.toArray();
    lengths = lengthList.toArray();
    ids = idList.toArray();
    final int[] limits = MathUtils.limits(lengths);
    minX = limits[0] - 1.0;
    maxX = limits[1] + 1.0;
    h1 = new int[limits[1] + 1];
    h2 = new int[h1.length];
    x1 = createHistogramAxis(h1.length, -WIDTH);
    x2 = createHistogramAxis(h1.length, 0f);
    y1 = new float[x1.length];
    y2 = new float[x1.length];

    // Sort by MSD
    final int[] indices = SimpleArrayUtils.natural(msds.length);
    SortUtils.sortIndices(indices, msds, false);
    final double[] msds2 = msds.clone();
    final int[] lengths2 = lengths.clone();
    final int[] ids2 = ids.clone();
    for (int i = 0; i < indices.length; i++) {
      msds[i] = msds2[indices[i]];
      lengths[i] = lengths2[indices[i]];
      ids[i] = ids2[indices[i]];
    }

    // Interactive analysis
    final NonBlockingExtendedGenericDialog gd = new NonBlockingExtendedGenericDialog(TITLE);
    gd.addMessage(String.format(
        "Split traces into fixed or moving using the track diffusion coefficient (D).\n"
            + "Localistion error has been subtracted from jumps (%s nm).",
        MathUtils.rounded(precision)));
    final Statistics s = Statistics.create(msds);
    final double av = s.getMean();
    final String msg = String.format("Average D per track = %s um^2/s", MathUtils.rounded(av));
    gd.addMessage(msg);
    // Histogram the diffusion coefficients
    final WindowOrganiser wo = new WindowOrganiser();
    final HistogramPlot histogramPlot =
        new HistogramPlotBuilder("Trace diffusion coefficient", new StoredData(msds), "D (um^2/s)")
            .setRemoveOutliersOption(1).setPlotLabel(msg).build();
    histogramPlot.show(wo);
    final double[] xvalues = histogramPlot.getPlotXValues();
    final double min = xvalues[0];
    final double max = xvalues[xvalues.length - 1];
    // see if we can build a nice slider range from the histogram limits
    if (max - min < 5) {
      // Because sliders are used when the range is <5 and floating point
      gd.addSlider("D_threshold", min, max, msdThreshold);
    } else {
      gd.addNumericField("D_threshold", msdThreshold, 2, 6, "um^2/s");
    }
    gd.addCheckbox("Normalise", normalise);
    gd.addDialogListener((gd1, event) -> {
      msdThreshold = gd1.getNextNumber();
      normalise = gd1.getNextBoolean();
      update();
      return true;
    });
    if (ImageJUtils.isShowGenericDialog()) {
      draw(wo);
      wo.tile();
    }
    gd.setOKLabel("Save datasets");
    gd.setCancelLabel("Close");
    gd.showDialog();

    if (gd.wasCanceled()) {
      return;
    }

    // Sort by ID
    final PeakResult[] list = results.toArray();
    Arrays.sort(list, new IdFramePeakResultComparator());

    createResults(results, "Fixed", 0, lastIndex, list);
    createResults(results, "Moving", lastIndex, msds.length, list);
  }

  private void createResults(MemoryPeakResults results, String suffix, int from, int to,
      PeakResult[] list) {
    final MemoryPeakResults out = new MemoryPeakResults();
    out.copySettings(results);
    out.setName(results.getName() + " " + suffix);

    // Sort target ids
    Arrays.sort(ids, from, to);

    for (int i = 0; i < list.length && from < to;) {
      final int nextId = ids[from++];
      // Move forward
      while (i < list.length && list[i].getId() < nextId) {
        i++;
      }
      // Write out
      while (i < list.length && list[i].getId() == nextId) {
        out.add(list[i++]);
      }
    }

    MemoryPeakResults.addResults(out);
  }

  private static boolean showDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addMessage("Analyse the track length of traced data");
    ResultsManager.addInput(gd, "Input", inputOption, InputSource.MEMORY_CLUSTERED);
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }
    inputOption = ResultsManager.getInputSource(gd);
    return true;
  }

  private void draw(WindowOrganiser wo) {
    lastMsdThreshold = msdThreshold;
    lastNormalise = normalise;

    // Find the index in the MSD array
    int index = Arrays.binarySearch(msds, lastMsdThreshold);
    if (index < 0) {
      index = -(index + 1);
    }
    lastIndex = index;

    // Histogram the distributions
    computeHistogram(0, index, lengths, h1);
    computeHistogram(index, lengths.length, lengths, h2);
    final int sum1 = (int) MathUtils.sum(h1);
    final int sum2 = (int) MathUtils.sum(h2);

    final float max1 = createHistogramValues(h1, (lastNormalise) ? sum1 : 1, y1);
    final float max2 = createHistogramValues(h2, (lastNormalise) ? sum2 : 2, y2);

    final String title = "Trace length distribution";
    final Plot plot = new Plot(title, "Length", "Frequency");
    plot.setLimits(minX, maxX, 0, Math.max(max1, max2) * 1.05);
    plot.setColor(Color.red);
    plot.addPoints(x1, y1, Plot.LINE);
    plot.setColor(Color.blue);
    plot.addPoints(x2, y2, Plot.LINE);
    plot.setColor(Color.black);
    final double p = 100.0 * sum1 / (sum1 + sum2);
    plot.addLabel(0, 0, String.format("Fixed (red) = %d (%s%%), Moving (blue) = %d (%s%%)", sum1,
        MathUtils.rounded(p), sum2, MathUtils.rounded(100 - p)));
    ImageJUtils.display(title, plot, ImageJUtils.NO_TO_FRONT, wo);
  }

  private static void computeHistogram(int index, int end, int[] length, int[] histogram) {
    Arrays.fill(histogram, 0);
    while (index < end) {
      histogram[length[index++]]++;
    }
  }

  /**
   * For the provided histogram x-axis bins, produce an x-axis for plotting. This functions doubles
   * up the histogram x-positions to allow plotting a square line profile using the ImageJ plot
   * command.
   *
   * @param length the length
   * @param offset the offset
   * @return the x-axis
   */
  public static float[] createHistogramAxis(int length, float offset) {
    final float[] axis = new float[length * 4];
    int index = 0;
    final float offset2 = offset + WIDTH;
    for (int i = 0; i < length; ++i) {
      axis[index++] = i + offset;
      axis[index++] = i + offset;
      axis[index++] = i + offset2;
      axis[index++] = i + offset2;
    }
    return axis;
  }

  /**
   * For the provided histogram y-axis values, produce a y-axis for plotting. This functions doubles
   * up the histogram values to allow plotting a square line profile using the ImageJ plot command.
   *
   * @param histogramY the histogram Y
   * @param norm the normalisation
   * @param axis the axis
   * @return the float
   */
  public static float createHistogramValues(int[] histogramY, double norm, float[] axis) {
    // Assume axis[0] and axis[axis.length-1] == 0
    int index = 1;
    float max = 0;
    for (int i = 0; i < histogramY.length; ++i) {
      final float v = (float) (histogramY[i] / norm);
      axis[index++] = v;
      axis[index++] = v;
      if (max < v) {
        max = v;
      }
      index += 2;
    }
    return max;
  }

  private void update() {
    if (lock.acquire()) {
      // Run in a new thread to allow the GUI to continue updating
      new Thread(() -> {
        try {
          // Continue while the parameter is changing
          while (lastMsdThreshold != msdThreshold || lastNormalise != normalise) {
            draw(null);
          }
        } finally {
          // Ensure the running flag is reset
          lock.release();
        }
      }).start();
    }
  }

  private int lastid = -1;
  private float lastx;
  private float lasty;
  private int startFrame;
  private int lastFrame;
  private int totalJump;
  private double sumSquared;
  private final TDoubleArrayList msdList = new TDoubleArrayList();
  private final TIntArrayList lengthList = new TIntArrayList();
  private final TIntArrayList idList = new TIntArrayList();

  private void processTrackLength(PeakResult peakResult) {
    final int id = peakResult.getId();
    final int frame = peakResult.getFrame();
    final float x = peakResult.getXPosition();
    final float y = peakResult.getYPosition();
    if (lastid != id) {
      store();
      lastid = id;
      startFrame = frame;
      sumSquared = 0;
      totalJump = 0;
    } else {
      // Compute the jump
      final int jump = frame - lastFrame;
      // Get the raw distance but subtract the expected localisation error
      final double d2 = Math.max(0, MathUtils.distance2(lastx, lasty, x, y) - error);
      // We expect the Mean Squared Distance (MSD) to scale linearly
      // with time so just weight each jump by the time gap.
      // However we apply a correction factor for diffusion with frames.
      sumSquared += JumpDistanceAnalysis.convertObservedToActual(d2, jump);
      totalJump += jump;
    }
    lastFrame = frame;
    lastx = x;
    lasty = y;
  }

  private void store() {
    if (lastid == 0 || totalJump == 0) {
      return;
    }
    final int len = lastFrame - startFrame;
    final double msd = sumSquared / totalJump;
    // 4D = MSD => D = MSD / 4
    // Mean squared distance is in raw units squared.
    // Convert twice since this is a squared distance then from frames to seconds.
    // Use the convertBack() since this is a divide in the final units: um^2/s
    msdList.add(
        timeConverter.convertBack(distanceConverter.convert(distanceConverter.convert(msd))) / 4.0);
    lengthList.add(len);
    idList.add(lastid);
  }
}
