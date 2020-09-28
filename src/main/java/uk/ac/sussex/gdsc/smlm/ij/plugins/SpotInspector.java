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
import ij.gui.Plot;
import ij.gui.PointRoi;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.text.TextPanel;
import java.awt.Rectangle;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot.HistogramPlotBuilder;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
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
import uk.ac.sussex.gdsc.smlm.results.ImageSource.ReadHint;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.count.Counter;
import uk.ac.sussex.gdsc.smlm.results.procedures.HeightResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.PrecisionResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.StandardResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.WidthResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.XyResultProcedure;

/**
 * Extract the spots from the original image into a stack, ordering the spots by various rankings.
 */
public class SpotInspector implements PlugIn {
  private static final String TITLE = "Spot Inspector";

  private MemoryPeakResults results;
  private TextPanel textPanel;
  private List<PeakResultRank> rankedResults;

  private AtomicInteger currentId = new AtomicInteger();
  private int id;

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    static final String[] SORT_ORDER = new String[] {"SNR", "Precision", "Amplitude", "Signal",
        "Error", "Original Value", "X SD", "Y SD", "Width factor", "Shift"};

    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    String inputOption;
    int sortOrderIndex;
    int radius;
    boolean showCalibratedValues;
    boolean plotScore;
    boolean plotHistogram;
    int histogramBins;
    boolean removeOutliers;

    Settings() {
      // Set defaults
      inputOption = "";
      sortOrderIndex = 1;
      radius = 5;
      showCalibratedValues = true;
      plotScore = true;
      plotHistogram = true;
      histogramBins = 0;
      removeOutliers = true;
    }

    Settings(Settings source) {
      inputOption = source.inputOption;
      sortOrderIndex = source.sortOrderIndex;
      radius = source.radius;
      showCalibratedValues = source.showCalibratedValues;
      plotScore = source.plotScore;
      plotHistogram = source.plotHistogram;
      histogramBins = source.histogramBins;
      removeOutliers = source.removeOutliers;
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
    results =
        ResultsManager.loadInputResults(settings.inputOption, false, DistanceUnit.PIXEL, null);
    if (MemoryPeakResults.isEmpty(results)) {
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
    source.setReadHint(ReadHint.NONSEQUENTIAL);
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
    if (settings.sortOrderIndex == 1) {
      pp = new PrecisionResultProcedure(results);
      pp.getPrecision();
    } else {
      pp = null;
    }

    // Build procedures to get:
    // Shift = position in pixels - originXY
    final StandardResultProcedure sp;
    if (settings.sortOrderIndex == 9) {
      sp = new StandardResultProcedure(results, DistanceUnit.PIXEL);
      sp.getXyr();
    } else {
      sp = null;
    }

    // SD = gaussian widths only for Gaussian PSFs
    final WidthResultProcedure wp;
    if (settings.sortOrderIndex >= 6 && settings.sortOrderIndex <= 8) {
      wp = new WidthResultProcedure(results, DistanceUnit.PIXEL);
      wp.getWxWy();
    } else {
      wp = null;
    }

    // Amplitude for Gaussian PSFs
    final HeightResultProcedure hp;
    if (settings.sortOrderIndex == 2) {
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

    if (settings.showCalibratedValues) {
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
    id = currentId.incrementAndGet();
    textPanel.addMouseListener(new MouseAdapter() {
      @Override
      public void mouseClicked(MouseEvent event) {
        SpotInspector.this.mouseClicked(event);
      }
    });

    // Add results to the table
    int count = 0;
    for (final PeakResultRank rank : rankedResults) {
      rank.rank = count++;
      table.add(rank.peakResult);
    }
    table.end();

    if (settings.plotScore || settings.plotHistogram) {
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
      if (settings.removeOutliers) {
        final double lower = stats.getPercentile(25);
        final double upper = stats.getPercentile(75);
        final double iqr = upper - lower;

        yMin = Math.max(lower - iqr, stats.getMin());
        yMax = Math.min(upper + iqr, stats.getMax());

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
    final int size = 2 * settings.radius + 1;
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
      int minX = x - settings.radius;
      int minY = y - settings.radius;
      final int maxX = Math.min(x + settings.radius + 1, w);
      final int maxY = Math.min(y + settings.radius + 1, h);

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

    // Reset to the rank order
    Collections.sort(rankedResults, PeakResultRank::compare);

    final ImagePlus imp = ImageJUtils.display(TITLE, spots);
    imp.setRoi((PointRoi) null);

    // Make bigger
    for (int i = 10; i-- > 0;) {
      imp.getWindow().getCanvas().zoomIn(imp.getWidth() / 2, imp.getHeight() / 2);
    }
  }

  private float getStandardDeviation(MemoryPeakResults results2) {
    // Standard deviation is only needed for the width filtering
    if (settings.sortOrderIndex != 8) {
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

  private void plotScore(float[] xValues, float[] yValues, double yMin, double yMax) {
    if (settings.plotScore) {
      final String title = TITLE + " Score";
      final Plot plot = new Plot(title, "Rank", Settings.SORT_ORDER[settings.sortOrderIndex]);
      plot.addPoints(xValues, yValues, Plot.LINE);
      plot.setLimits(1, xValues.length, yMin, yMax);

      ImageJUtils.display(title, plot);
    }
  }

  private void plotHistogram(float[] data) {
    if (settings.plotHistogram) {
      final String title = TITLE + " Histogram";
      new HistogramPlotBuilder(title, StoredDataStatistics.create(data),
          Settings.SORT_ORDER[settings.sortOrderIndex])
              .setRemoveOutliersOption((settings.removeOutliers) ? 1 : 0)
              .setNumberOfBins(settings.histogramBins).show();
    }
  }

  private float[] getScore(PeakResult result, int index, PrecisionResultProcedure pp,
      StandardResultProcedure sp, WidthResultProcedure wp, HeightResultProcedure hp,
      float stdDevMax) {
    // Return score so high is better
    float score;
    boolean negative = false;
    switch (settings.sortOrderIndex) {
      case 9: // Shift
        // We do not have the original centroid so use the original X/Y
        score = Math.max(sp.x[index] - result.getOrigX() + 0.5f,
            sp.y[index] - result.getOrigY() + 0.5f);
        negative = true;
        break;
      case 8: // Width factor
        score = getFactor(Math.max(wp.wx[index], wp.wy[index]), stdDevMax);
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

  private float recoverScore(float score) {
    // Reset the sign of the score
    switch (settings.sortOrderIndex) {
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

  private boolean showDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("spot-inspector"));

    settings = Settings.load();

    ResultsManager.addInput(gd, settings.inputOption, InputSource.MEMORY);

    gd.addChoice("Ranking", Settings.SORT_ORDER, settings.sortOrderIndex);
    gd.addSlider("Radius", 1, 15, settings.radius);
    gd.addCheckbox("Calibrated_table", settings.showCalibratedValues);
    gd.addCheckbox("Plot_score", settings.plotScore);
    gd.addCheckbox("Plot_histogram", settings.plotHistogram);
    gd.addNumericField("Histogram_bins", settings.histogramBins, 0);
    gd.addCheckbox("Remove_outliers", settings.removeOutliers);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.inputOption = ResultsManager.getInputSource(gd);
    settings.sortOrderIndex = gd.getNextChoiceIndex();
    settings.radius = (int) gd.getNextNumber();
    settings.showCalibratedValues = gd.getNextBoolean();
    settings.plotScore = gd.getNextBoolean();
    settings.plotHistogram = gd.getNextBoolean();
    settings.histogramBins = (int) gd.getNextNumber();
    settings.removeOutliers = gd.getNextBoolean();
    settings.save();

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Radius", settings.radius);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return true;
  }

  private void mouseClicked(MouseEvent event) {
    if (id != currentId.get()) {
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

        final PeakResult r = rankedResults.get(rank - 1).peakResult;

        final TypeConverter<DistanceUnit> dc = results.getDistanceConverter(DistanceUnit.PIXEL);

        final float rx = dc.convert(r.getXPosition());
        final float ry = dc.convert(r.getYPosition());
        final int x = (int) rx;
        final int y = (int) ry;

        // Find bounds
        final int minX = x - settings.radius;
        final int minY = y - settings.radius;

        // Add an ROI to the image containing the clicked spot / all spots in the region.

        // Require the Shift key to add all spots
        if (!event.isShiftDown()) {
          // Add the single clicked spot
          imp.setRoi(new PointRoi(rx - minX, ry - minY));
          return;
        }

        // Add all the spots
        final int maxX = x + settings.radius + 1;
        final int maxY = y + settings.radius + 1;

        // Create ROIs
        final HashSet<Point2D.Float> spots = new HashSet<>();
        results.forEach(DistanceUnit.PIXEL, (XyResultProcedure) (xp, yp) -> {
          if (xp > minX && xp < maxX && yp > minY && yp < maxY) {
            // Use only unique points
            spots.add(new Point2D.Float(xp - minX, yp - minY));
          }
        });

        final int points = spots.size();
        final float[] ox = new float[points];
        final float[] oy = new float[points];
        final Counter c = new Counter();
        spots.forEach(p -> {
          ox[c.getCount()] = p.x;
          oy[c.getAndIncrement()] = p.y;
        });
        imp.setRoi(new PointRoi(ox, oy, points));
      }
    }
  }
}
