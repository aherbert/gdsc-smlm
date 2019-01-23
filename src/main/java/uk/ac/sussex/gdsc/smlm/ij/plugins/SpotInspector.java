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

import uk.ac.sussex.gdsc.core.ij.HistogramPlot.HistogramPlotBuilder;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.Plot2;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.data.config.PsfHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.engine.FitConfiguration;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.results.ImageJTablePeakResults;
import uk.ac.sussex.gdsc.smlm.results.ImageSource;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.count.Counter;
import uk.ac.sussex.gdsc.smlm.results.procedures.HeightResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.PrecisionResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.StandardResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.WidthResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.XyResultProcedure;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.PointRoi;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.text.TextPanel;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.FastMath;

import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Extract the spots from the original image into a stack, ordering the spots by various rankings.
 */
public class SpotInspector implements PlugIn, MouseListener {
  private static final String TITLE = "Spot Inspector";

  private static String inputOption = "";
  private static String[] SORT_ORDER = new String[] {"SNR", "Precision", "Amplitude", "Signal",
      "Error", "Original Value", "X SD", "Y SD", "Width factor", "Shift"};
  private static int sortOrderIndex = 1;
  private static int radius = 5;
  private static boolean showCalibratedValues = true;
  private static boolean plotScore = true;
  private static boolean plotHistogram = true;
  private static int histogramBins = 100;
  private static boolean removeOutliers = true;

  private MemoryPeakResults results;
  private TextPanel textPanel;
  private List<PeakResultRank> rankedResults;

  private static int currentId;
  private int id;

  private static class PeakResultRank {
    int rank;
    PeakResult peakResult;
    float score;
    float originalScore;

    public PeakResultRank(PeakResult result, float score, float original) {
      peakResult = result;
      this.score = score;
      originalScore = original;
    }

    static int compare(PeakResultRank o1, PeakResultRank o2) {
      // High is better
      return Double.compare(o2.score, o1.score);
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
    results = ResultsManager.loadInputResults(inputOption, false, DistanceUnit.PIXEL, null);
    if (results == null || results.size() == 0) {
      IJ.error(TITLE, "No results could be loaded");
      IJ.showStatus("");
      return;
    }

    // Check if the original image is open
    ImageSource source = results.getSource();
    if (source == null) {
      IJ.error(TITLE, "Unknown original source image");
      return;
    }
    source = source.getOriginal();
    if (!source.open()) {
      IJ.error(TITLE, "Cannot open original source image: " + source.toString());
      return;
    }
    final float stdDevMax = getStandardDeviation(results);
    if (stdDevMax < 0) {
      // TODO - Add dialog to get the initial peak width
      IJ.error(TITLE, "Fitting configuration (for initial peak width) is not available");
      return;
    }

    // Rank spots
    rankedResults = new ArrayList<>(results.size());

    // Data for the sorting
    final PrecisionResultProcedure pp;
    if (sortOrderIndex == 1) {
      pp = new PrecisionResultProcedure(results);
      pp.getPrecision();
    } else {
      pp = null;
    }

    // Build procedures to get:
    // Shift = position in pixels - originXY
    final StandardResultProcedure sp;
    if (sortOrderIndex == 9) {
      sp = new StandardResultProcedure(results, DistanceUnit.PIXEL);
      sp.getXyr();
    } else {
      sp = null;
    }

    // SD = gaussian widths only for Gaussian PSFs
    final WidthResultProcedure wp;
    if (sortOrderIndex >= 6 && sortOrderIndex <= 8) {
      wp = new WidthResultProcedure(results, DistanceUnit.PIXEL);
      wp.getWxWy();
    } else {
      wp = null;
    }

    // Amplitude for Gaussian PSFs
    final HeightResultProcedure hp;
    if (sortOrderIndex == 2) {
      hp = new HeightResultProcedure(results, IntensityUnit.PHOTON);
      hp.getH();
    } else {
      hp = null;
    }

    final Counter c = new Counter();
    results.forEach((PeakResultProcedure) result -> {
      final float[] score = getScore(result, c.getAndIncrement(), pp, sp, wp, hp, stdDevMax);
      rankedResults.add(new PeakResultRank(result, score[0], score[1]));
    });
    Collections.sort(rankedResults, PeakResultRank::compare);

    // Prepare results table
    final ImageJTablePeakResults table = new ImageJTablePeakResults(false, results.getName(), true);
    table.copySettings(results);
    table.setTableTitle(TITLE);
    table.setAddCounter(true);
    table.setShowZ(results.is3D());
    // TODO - Add to settings
    table.setShowFittingData(true);
    table.setShowNoiseData(true);

    if (showCalibratedValues) {
      table.setDistanceUnit(DistanceUnit.NM);
      table.setIntensityUnit(IntensityUnit.PHOTON);
    } else {
      table.setDistanceUnit(DistanceUnit.PIXEL);
      table.setIntensityUnit(IntensityUnit.COUNT);
    }
    table.begin();

    // Add a mouse listener to jump to the frame for the clicked line
    textPanel = table.getResultsWindow().getTextPanel();

    // We must ignore old instances of this class from the mouse listeners
    id = ++currentId;
    textPanel.addMouseListener(this);

    // Add results to the table
    int count = 0;
    for (final PeakResultRank rank : rankedResults) {
      rank.rank = count++;
      table.add(rank.peakResult);
    }
    table.end();

    if (plotScore || plotHistogram) {
      // Get values for the plots
      float[] xValues = null;
      float[] yValues = null;
      double yMin;
      double yMax;

      int spotNumber = 0;
      xValues = new float[rankedResults.size()];
      yValues = new float[xValues.length];
      for (final PeakResultRank rank : rankedResults) {
        xValues[spotNumber] = spotNumber + 1;
        yValues[spotNumber++] = recoverScore(rank.score);
      }

      // Set the min and max y-values using 1.5 x IQR
      final DescriptiveStatistics stats = new DescriptiveStatistics();
      for (final float v : yValues) {
        stats.addValue(v);
      }
      if (removeOutliers) {
        final double lower = stats.getPercentile(25);
        final double upper = stats.getPercentile(75);
        final double iqr = upper - lower;

        yMin = FastMath.max(lower - iqr, stats.getMin());
        yMax = FastMath.min(upper + iqr, stats.getMax());

        IJ.log(String.format("Data range: %f - %f. Plotting 1.5x IQR: %f - %f", stats.getMin(),
            stats.getMax(), yMin, yMax));
      } else {
        yMin = stats.getMin();
        yMax = stats.getMax();

        IJ.log(String.format("Data range: %f - %f", yMin, yMax));
      }

      plotScore(xValues, yValues, yMin, yMax);
      plotHistogram(yValues);
    }

    // Extract spots into a stack
    final int w = source.getWidth();
    final int h = source.getHeight();
    final int size = 2 * radius + 1;
    final ImageStack spots = new ImageStack(size, size, rankedResults.size());

    // To assist the extraction of data from the image source, process them in time order to allow
    // frame caching. Then set the appropriate slice in the result stack
    Collections.sort(rankedResults,
        (o1, o2) -> Integer.compare(o1.peakResult.getFrame(), o2.peakResult.getFrame()));

    for (final PeakResultRank rank : rankedResults) {
      final PeakResult r = rank.peakResult;

      // Extract image
      // Note that the coordinates are relative to the middle of the pixel (0.5 offset)
      // so do not round but simply convert to int
      final int x = (int) (r.getXPosition());
      final int y = (int) (r.getYPosition());

      // Extract a region but crop to the image bounds
      int minX = x - radius;
      int minY = y - radius;
      final int maxX = FastMath.min(x + radius + 1, w);
      final int maxY = FastMath.min(y + radius + 1, h);

      int padX = 0;
      int padY = 0;
      if (minX < 0) {
        padX = -minX;
        minX = 0;
      }
      if (minY < 0) {
        padY = -minY;
        minY = 0;
      }
      final int sizeX = maxX - minX;
      final int sizeY = maxY - minY;

      float[] data = source.get(r.getFrame(), new Rectangle(minX, minY, sizeX, sizeY));
      // Prevent errors with missing data
      if (data == null) {
        data = new float[sizeX * sizeY];
      }
      ImageProcessor spotIp = new FloatProcessor(sizeX, sizeY, data, null);

      // Pad if necessary, i.e. the crop is too small for the stack
      if (padX > 0 || padY > 0 || sizeX < size || sizeY < size) {
        final ImageProcessor spotIp2 = spotIp.createProcessor(size, size);
        spotIp2.insert(spotIp, padX, padY);
        spotIp = spotIp2;
      }
      final int slice = rank.rank + 1;
      spots.setPixels(spotIp.getPixels(), slice);
      spots.setSliceLabel(MathUtils.rounded(rank.originalScore), slice);
    }

    source.close();

    final ImagePlus imp = ImageJUtils.display(TITLE, spots);
    imp.setRoi((PointRoi) null);

    // Make bigger
    for (int i = 10; i-- > 0;) {
      imp.getWindow().getCanvas().zoomIn(imp.getWidth() / 2, imp.getHeight() / 2);
    }
  }

  private static float getStandardDeviation(MemoryPeakResults results2) {
    // Standard deviation is only needed for the width filtering
    if (sortOrderIndex != 8) {
      return 0;
    }
    final PSF psf = results2.getPsf();
    if (!PsfHelper.isGaussian2D(psf)) {
      return -1;
    }
    final FitConfiguration fitConfig = new FitConfiguration();
    fitConfig.setPsf(psf);
    return (float) MathUtils.max(1, fitConfig.getInitialXSd(), fitConfig.getInitialYSd());
  }

  private static void plotScore(float[] xValues, float[] yValues, double yMin, double yMax) {
    if (plotScore) {
      final String title = TITLE + " Score";
      final Plot2 plot = new Plot2(title, "Rank", SORT_ORDER[sortOrderIndex], xValues, yValues);
      plot.setLimits(1, xValues.length, yMin, yMax);

      ImageJUtils.display(title, plot);
    }
  }

  private static void plotHistogram(float[] data) {
    if (plotHistogram) {
      final String title = TITLE + " Histogram";
      new HistogramPlotBuilder(title, StoredDataStatistics.create(data), SORT_ORDER[sortOrderIndex])
          .setRemoveOutliersOption((removeOutliers) ? 1 : 0).setNumberOfBins(histogramBins).show();
    }
  }

  private static float[] getScore(PeakResult result, int index, PrecisionResultProcedure pp,
      StandardResultProcedure sp, WidthResultProcedure wp, HeightResultProcedure hp,
      float stdDevMax) {
    // Return score so high is better
    float score;
    boolean negative = false;
    switch (sortOrderIndex) {
      case 9: // Shift
        // We do not have the original centroid so use the original X/Y
        score = FastMath.max(sp.x[index] - result.getOrigX() + 0.5f,
            sp.y[index] - result.getOrigY() + 0.5f);
        negative = true;
        break;
      case 8: // Width factor
        score = getFactor(FastMath.max(wp.wx[index], wp.wy[index]), stdDevMax);
        negative = true;
        break;
      case 7:
        score = wp.wy[index];
        negative = true;
        break;
      case 6:
        score = wp.wx[index];
        negative = true;
        break;
      case 5: // Original value
        score = result.getOrigValue();
        break;
      case 4: // Error
        score = (float) result.getError();
        negative = true;
        break;
      case 3: // Signal
        score = (result.getIntensity());
        break;
      case 2: // Amplitude
        score = hp.heights[index];
        break;
      case 1: // Precision
        score = (float) pp.precisions[index];
        negative = true;
        break;
      default: // SNR
        score = result.getIntensity() / result.getNoise();
    }
    return new float[] {(negative) ? -score : score, score};
  }

  private static float recoverScore(float score) {
    // Reset the sign of the score
    switch (sortOrderIndex) {
      case 9: // Shift
        return -score;
      case 8: // Width factor
        return -score;
      case 7:
        return -score;
      case 6:
        return -score;
      case 5: // Original value
        return score;
      case 4: // Error
        return -score;
      case 3: // Signal
        return score;
      case 2: // Amplitude
        return score;
      case 1: // Precision
        return -score;
      default: // SNR
        return score;
    }
  }

  /**
   * Get the relative change factor between v1 and v2.
   *
   * @param v1 the first value
   * @param v2 the second value
   * @return the factor
   */
  private static float getFactor(float v1, float v2) {
    if (v1 > v2) {
      return v1 / v2;
    }
    return v2 / v1;
  }

  private static boolean showDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);

    ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);

    gd.addChoice("Ranking", SORT_ORDER, SORT_ORDER[sortOrderIndex]);
    gd.addSlider("Radius", 1, 15, radius);
    gd.addCheckbox("Calibrated_table", showCalibratedValues);
    gd.addCheckbox("Plot_score", plotScore);
    gd.addCheckbox("Plot_histogram", plotHistogram);
    gd.addNumericField("Histogram_bins", histogramBins, 0);
    gd.addCheckbox("Remove_outliers", removeOutliers);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    inputOption = ResultsManager.getInputSource(gd);
    sortOrderIndex = gd.getNextChoiceIndex();
    radius = (int) gd.getNextNumber();
    showCalibratedValues = gd.getNextBoolean();
    plotScore = gd.getNextBoolean();
    plotHistogram = gd.getNextBoolean();
    histogramBins = (int) gd.getNextNumber();
    removeOutliers = gd.getNextBoolean();

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Radius", radius);
      ParameterUtils.isAbove("Histogram bins", histogramBins, 1);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return true;
  }

  @Override
  public void mouseClicked(MouseEvent event) {
    if (id != currentId) {
      return;
    }
    // Show the result that was double clicked in the result table
    if (event.getClickCount() > 1) {
      final int rank = textPanel.getSelectionStart() + 1;

      // Show the spot that was double clicked
      final ImagePlus imp = WindowManager.getImage(TITLE);
      if (imp != null && rank > 0 && rank <= imp.getStackSize()) {
        imp.setSlice(rank);
        if (imp.getWindow() != null) {
          imp.getWindow().toFront();
        }

        // Add an ROI to to the image containing all the spots in that part of the frame
        final PeakResult r = rankedResults.get(rank - 1).peakResult;

        final int x = (int) (r.getXPosition());
        final int y = (int) (r.getYPosition());

        // Find bounds
        final int minX = x - radius;
        final int minY = y - radius;
        final int maxX = x + radius + 1;
        final int maxY = y + radius + 1;

        // Create ROIs
        final ArrayList<float[]> spots = new ArrayList<>();
        results.forEach(DistanceUnit.PIXEL, new XyResultProcedure() {
          @Override
          public void executeXy(float x, float y) {
            if (x > minX && x < maxX && y > minY && y < maxY) {
              // Use only unique points
              final float xPosition = x - minX;
              final float yPosition = y - minY;
              if (!contains(spots, xPosition, yPosition)) {
                spots.add(new float[] {xPosition, yPosition});
              }
            }
          }
        });

        final int points = spots.size();
        final float[] ox = new float[points];
        final float[] oy = new float[points];
        for (int i = 0; i < points; i++) {
          ox[i] = spots.get(i)[0];
          oy[i] = spots.get(i)[1];
        }
        imp.setRoi(new PointRoi(ox, oy, points));
      }
    }
  }

  private static boolean contains(ArrayList<float[]> spots, float xPosition, float yPosition) {
    for (final float[] data : spots) {
      if (data[0] == xPosition && data[1] == yPosition) {
        return true;
      }
    }
    return false;
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
