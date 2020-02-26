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

import gnu.trove.list.array.TDoubleArrayList;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.plugin.PlugIn;
import ij.plugin.frame.Recorder;
import java.awt.Color;
import java.awt.Rectangle;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.atomic.AtomicReference;
import org.apache.commons.math3.random.HaltonSequenceGenerator;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.rng.simple.internal.SeedFactory;
import uk.ac.sussex.gdsc.core.clustering.DensityManager;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.Plot2;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageMode;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageType;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.results.ImageJImagePeakResults;
import uk.ac.sussex.gdsc.smlm.ij.results.ImagePeakResultsFactory;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.procedures.StandardResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.XyrResultProcedure;

/**
 * Produces an image on localisation using their density.
 */
public class DensityImage implements PlugIn {
  private static final String TITLE = "Density Image";

  private Rectangle roiBounds;
  private double scaledRoiMinX;
  private double scaledRoiMaxX;
  private double scaledRoiMinY;
  private double scaledRoiMaxY;
  private int roiImageWidth;
  private int roiImageHeight;

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    static final String[] scoreMethods = new String[] {"Density", "Ripley's K", "Ripley's K / Area",
        "Ripley's L", "Ripley's L - r", "Ripley's L / r", "Ripley's (L - r) / r"};

    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    String inputOption;
    float radius;
    boolean chooseRoi;
    String roiImage;
    boolean adjustForBorder;
    int imageScale;
    boolean cumulativeImage;
    boolean useSquareApproximation;
    int resolution;

    int scoreMethodIndex;

    boolean filterLocalisations;
    double filterThreshold;
    boolean computeRipleysPlot;

    double minR;
    double maxR;
    double incrementR;
    boolean confidenceIntervals;

    Settings() {
      // Set defaults
      inputOption = "";
      radius = 1.5f;
      roiImage = "";
      adjustForBorder = true;
      imageScale = 2;
      resolution = 10;
      filterLocalisations = true;
      minR = 0.2;
      maxR = 3;
      incrementR = 0.2;
    }

    Settings(Settings source) {
      inputOption = source.inputOption;
      radius = source.radius;
      chooseRoi = source.chooseRoi;
      roiImage = source.roiImage;
      adjustForBorder = source.adjustForBorder;
      imageScale = source.imageScale;
      cumulativeImage = source.cumulativeImage;
      useSquareApproximation = source.useSquareApproximation;
      resolution = source.resolution;
      scoreMethodIndex = source.scoreMethodIndex;
      filterLocalisations = source.filterLocalisations;
      filterThreshold = source.filterThreshold;
      computeRipleysPlot = source.computeRipleysPlot;
      minR = source.minR;
      maxR = source.maxR;
      incrementR = source.incrementR;
      confidenceIntervals = source.confidenceIntervals;
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

    // Require some fit results and selected regions
    final int size = MemoryPeakResults.countMemorySize();
    if (size == 0) {
      IJ.error(TITLE, "There are no fitting results in memory");
      return;
    }

    if (!showDialog()) {
      return;
    }

    MemoryPeakResults results =
        ResultsManager.loadInputResults(settings.inputOption, false, DistanceUnit.PIXEL, null);
    if (MemoryPeakResults.isEmpty(results)) {
      IJ.error(TITLE, "No results could be loaded");
      IJ.showStatus("");
      return;
    }

    final boolean[] isWithin = new boolean[1];
    results = cropWithBorder(results, isWithin);
    if (results.size() == 0) {
      IJ.error(TITLE, "No results within the crop region");
      IJ.showStatus("");
      return;
    }

    final long start = System.currentTimeMillis();
    IJ.showStatus("Calculating density ...");

    final boolean useAdjustment = settings.adjustForBorder && !isWithin[0];

    final DensityManager dm = createDensityManager(results);
    int[] density = null;
    if (settings.useSquareApproximation) {
      density = dm.calculateSquareDensity(settings.radius, settings.resolution, useAdjustment);
    } else {
      density = dm.calculateDensity(settings.radius, useAdjustment);
    }

    density = cropBorder(results, density);

    // Convert to float
    final ScoreCalculator calc = createCalculator(results);
    final float[] densityScore = calc.calculate(density);

    final int filtered = plotResults(results, densityScore, calc);

    logDensityResults(results, density, settings.radius, filtered);

    if (settings.computeRipleysPlot) {
      computeRipleysPlot(results);
    }

    final double seconds = (System.currentTimeMillis() - start) / 1000.0;
    IJ.showStatus(TITLE + " complete : " + seconds + "s");
  }

  private static DensityManager createDensityManager(MemoryPeakResults results) {
    if (MemoryPeakResults.isEmpty(results)) {
      throw new IllegalArgumentException("Results are null or empty");
    }

    final StandardResultProcedure sp = new StandardResultProcedure(results, DistanceUnit.PIXEL);
    sp.getXy();
    return new DensityManager(sp.x, sp.y, results.getBounds());
  }

  private ScoreCalculator createCalculator(MemoryPeakResults results) {
    switch (settings.scoreMethodIndex) {
      case 1:
        // Ripley's K (Density / av. density)
        return new KScoreCalculator(results, 0);

      case 2:
        // Ripley's K / area
        return new KScoreCalculator(results, 1);

      case 3:
        // Ripley's L
        return new LScoreCalculator(results, 0);

      case 4:
        // Ripley's L - r
        return new LScoreCalculator(results, 1);

      case 5:
        // Ripley's L / r
        return new LScoreCalculator(results, 2);

      case 6:
        // Ripley's (L - r) / r
        return new LScoreCalculator(results, 3);

      case 0:
      default:
        return new DensityScoreCalculator(results);
    }
  }

  private interface ScoreCalculator {
    /**
     * Get the density score for the input density counts.
     *
     * @param density the density
     * @return the float[]
     */
    float[] calculate(int[] density);

    /**
     * Get the score threshold for filtering results using the configured filter threshold.
     *
     * @return the threshold
     */
    float getThreshold();
  }

  private class DensityScoreCalculator implements ScoreCalculator {
    MemoryPeakResults results;

    public DensityScoreCalculator(MemoryPeakResults results) {
      this.results = results;
    }

    @Override
    public float[] calculate(int[] density) {
      final float[] score = new float[density.length];
      for (int i = 0; i < score.length; i++) {
        score[i] = density[i];
      }
      return score;
    }

    protected float getAverageDensity() {
      final Rectangle bounds = results.getBounds();
      final float area = (float) bounds.width * bounds.height;
      return results.size() / area;
    }

    protected float getRegionArea() {
      return settings.radius * settings.radius
          * ((settings.useSquareApproximation) ? 4 : (float) Math.PI);
    }

    @Override
    public float getThreshold() {
      final float expected = getAverageDensity() * getRegionArea();
      return (float) (expected * settings.filterThreshold);
    }
  }

  private class KScoreCalculator extends DensityScoreCalculator {
    int mode;

    public KScoreCalculator(MemoryPeakResults results, int mode) {
      super(results);
      this.mode = mode;
    }

    @Override
    public float[] calculate(int[] density) {
      final float[] score = new float[density.length];
      // K(r)
      float regionDivisor = getAverageDensity();
      if (mode == 1) {
        // K(r) / area
        regionDivisor *= getRegionArea();
      }
      for (int i = 0; i < score.length; i++) {
        score[i] = density[i] / regionDivisor;
      }
      return score;
    }

    @Override
    public float getThreshold() {
      // Note: K(r) ~ Area
      // Since K(r) should be equal to the area to make the filter threshold scale appropriately
      // we adjust the threshold by the area
      if (mode == 0) {
        return (float) settings.filterThreshold * getRegionArea();
      }
      // K(r) / area == 1
      // => no adjustment as this is radius scale independent
      return (float) settings.filterThreshold;
    }
  }

  private class LScoreCalculator extends KScoreCalculator {
    public LScoreCalculator(MemoryPeakResults results, int mode) {
      super(results, mode);
    }

    @Override
    public float[] calculate(int[] density) {
      // Compute a normalised variance stabilised per particle L-score.
      // As in: Scarselli, et al.
      // Cell type-specific β2-adrenergic receptor clusters identified using PALM
      // microscopy are not lipid raft related, but depend on actin cytoskeleton integrity.
      // J Biol Chem. 2012 May 11;287(20):16768-80
      // Note:
      // I have re-arranged the score to be:
      // Li(r) = Math.sqrt((Sample density / Average density) / pi) - r
      // This should be above zero if the density around the spot is higher than the average sample
      // density.

      final float[] score = new float[density.length];
      final float regionDivisor =
          getAverageDensity() * ((settings.useSquareApproximation) ? 4 : (float) Math.PI);
      for (int i = 0; i < score.length; i++) {
        // L(r)
        score[i] = (float) Math.sqrt(density[i] / regionDivisor);
      }
      if (mode == 1 || mode == 3) {
        // L(r) - r
        // (L(r) - r) / r
        for (int i = 0; i < score.length; i++) {
          score[i] -= settings.radius;
        }
      }
      if (mode == 2 || mode == 3) {
        // L(r) / r
        // (L(r) - r) / r
        for (int i = 0; i < score.length; i++) {
          score[i] /= settings.radius;
        }
      }
      return score;
    }

    @Override
    public float getThreshold() {
      // Note:
      // L(r) is proportional to radius
      // K(r) is proportional to area
      // To make the filtered results the same to the K(r) function we could use the
      // sqrt of the filterThreshold

      final double threshold = settings.filterThreshold;

      // Note: L(r) ~ r
      // Since L(r) should be equal to the radius to make the filter threshold scale appropriately
      // we adjust the threshold by the radius
      if (mode == 0) {
        return (float) threshold * settings.radius;
      }
      // L(r) - r == 0
      if (mode == 1) {
        return (float) threshold * settings.radius - settings.radius;
      }
      // L(r) - r / r == 0
      if (mode == 3) {
        return (float) (threshold * settings.radius - settings.radius) / settings.radius;
      }
      // L(r) / r == 1
      // => no adjustment as this is radius scale independent
      return (float) threshold;
    }
  }

  /**
   * Crop the results to the ROI. Add a border using the sampling radius so that counts do not have
   * to be approximated (i.e. all spots around the edge of the ROI will have a complete image to
   * sample from). The results are modified in place.
   *
   * @param results the results
   * @param isWithin Set to true if the added border is within the original bounds (i.e. no
   *        adjustment for missing counts is required)
   * @return the cropped results
   */
  private MemoryPeakResults cropWithBorder(MemoryPeakResults results, boolean[] isWithin) {
    isWithin[0] = false;
    if (roiBounds == null) {
      return results;
    }

    // Adjust bounds relative to input results image:
    // Use the ROI relative to the frame the ROI is drawn on.
    // Map those fractional coordinates back to the original data bounds.
    final Rectangle bounds = results.getBounds();
    final double xscale = (double) roiImageWidth / bounds.width;
    final double yscale = (double) roiImageHeight / bounds.height;

    // Compute relative to the results bounds (if present)
    scaledRoiMinX = bounds.x + roiBounds.x / xscale;
    scaledRoiMaxX = scaledRoiMinX + roiBounds.width / xscale;
    scaledRoiMinY = bounds.y + roiBounds.y / yscale;
    scaledRoiMaxY = scaledRoiMinY + roiBounds.height / yscale;

    // Allow for the border
    final float minX = (int) (scaledRoiMinX - settings.radius);
    final float maxX = (int) Math.ceil(scaledRoiMaxX + settings.radius);
    final float minY = (int) (scaledRoiMinY - settings.radius);
    final float maxY = (int) Math.ceil(scaledRoiMaxY + settings.radius);

    // Create a new set of results within the bounds
    final MemoryPeakResults newResults = new MemoryPeakResults();
    newResults.begin();
    results.forEach(DistanceUnit.PIXEL, (XyrResultProcedure) (x, y, result) -> {
      if (x >= minX && x <= maxX && y >= minY && y <= maxY) {
        newResults.add(result);
      }
    });
    newResults.end();
    newResults.copySettings(results);
    newResults
        .setBounds(new Rectangle((int) minX, (int) minY, (int) (maxX - minX), (int) (maxY - minY)));
    isWithin[0] = minX >= bounds.x && minY >= bounds.y && maxX <= (bounds.x + bounds.width)
        && maxY <= (bounds.y + bounds.height);
    return newResults;
  }

  /**
   * Remove any results which fall in the radius border added around the ROI. Results are modified
   * in place and a new density array is returned.
   *
   * @param results the results
   * @param density the density
   * @return the density array
   */
  private int[] cropBorder(MemoryPeakResults results, int[] density) {
    if (roiBounds == null) {
      return density;
    }
    final float minX = (int) (scaledRoiMinX);
    final float maxX = (int) Math.ceil(scaledRoiMaxX);
    final float minY = (int) (scaledRoiMinY);
    final float maxY = (int) Math.ceil(scaledRoiMaxY);
    // Clone the results then add back those that are within the bounds
    final PeakResult[] peakResults = results.toArray();
    results.begin();
    int count = 0;
    for (int i = 0; i < peakResults.length; i++) {
      final PeakResult peakResult = peakResults[i];
      final float x = peakResult.getXPosition();
      final float y = peakResult.getYPosition();
      if (x < minX || x > maxX || y < minY || y > maxY) {
        continue;
      }
      results.add(peakResult);
      density[count++] = density[i];
    }
    results.end();
    results
        .setBounds(new Rectangle((int) minX, (int) minY, (int) (maxX - minX), (int) (maxY - minY)));
    return Arrays.copyOf(density, count);
  }

  /**
   * Output a log message of the results including the average density for localisations and the
   * expected average.
   *
   * @param results the results
   * @param density the density
   * @param radius the radius
   * @param filtered the filtered
   * @return the summary statistics
   */
  private SummaryStatistics logDensityResults(MemoryPeakResults results, int[] density,
      float radius, int filtered) {
    final float region =
        (float) (radius * radius * ((settings.useSquareApproximation) ? 4 : Math.PI));

    final Rectangle bounds = results.getBounds();
    final float area = (float) bounds.width * bounds.height;
    final float expected = results.size() * region / area;
    final SummaryStatistics summary = new SummaryStatistics();

    for (int i = 0; i < results.size(); i++) {
      summary.addValue(density[i]);
    }

    final DensityManager dm = createDensityManager(results);

    // Compute this using the input density scores since the radius is the same.
    final double l = (settings.useSquareApproximation) ? dm.ripleysLFunction(radius)
        : dm.ripleysLFunction(density, radius);

    String msg = String.format(
        "Density %s : N=%d, %.0fpx : Radius=%s : L(r) - r = %s : E = %s, Obs = %s (%sx)",
        results.getName(), summary.getN(), area, rounded(radius), rounded(l - radius),
        rounded(expected), rounded(summary.getMean()), rounded(summary.getMean() / expected));
    if (settings.filterLocalisations) {
      msg += String.format(" : Filtered=%d (%s%%)", filtered,
          rounded(filtered * 100.0 / density.length));
    }
    IJ.log(msg);

    return summary;
  }

  private static String rounded(double value) {
    return MathUtils.rounded(value, 3);
  }

  /**
   * Draw an image of the density for each localisation. Optionally filter results below a
   * threshold.
   *
   * @param results the results
   * @param density the density
   * @param scoreCalculator the score calculator
   * @return the number of localisations drawn
   */
  private int plotResults(MemoryPeakResults results, float[] density,
      ScoreCalculator scoreCalculator) {
    // Filter results using 5x higher than average density of the sample in a 150nm radius:
    // Annibale, et al (2011). Identification of clustering artifacts in photoactivated localization
    // microscopy.
    // Nature Methods, 8, pp527-528
    MemoryPeakResults newResults = null;
    float densityThreshold = Float.NEGATIVE_INFINITY; // No filtering
    if (settings.filterLocalisations) {
      densityThreshold = scoreCalculator.getThreshold();
      newResults = new MemoryPeakResults();
      newResults.copySettings(results);
      newResults.setName(results.getName() + " Density Filter");
    }

    // Draw an image - Use error so that a floating point value can be used on a single pixel
    final ImageJImagePeakResults image = ImagePeakResultsFactory.createPeakResultsImage(
        ResultsImageType.DRAW_INTENSITY, false, false, results.getName() + " Density",
        results.getBounds(), results.getNmPerPixel(), results.getGain(), settings.imageScale, 0,
        (settings.cumulativeImage) ? ResultsImageMode.IMAGE_ADD : ResultsImageMode.IMAGE_MAX);
    image.setDisplayFlags(image.getDisplayFlags() | ImageJImagePeakResults.DISPLAY_NEGATIVES);
    image.setLutName("grays");
    image.begin();
    final StandardResultProcedure sp = new StandardResultProcedure(newResults, DistanceUnit.PIXEL);
    sp.getXyr();
    for (int i = 0; i < density.length; i++) {
      if (density[i] < densityThreshold) {
        continue;
      }
      image.add(sp.x[i], sp.y[i], density[i]);
      if (newResults != null) {
        newResults.add(sp.peakResults[i]);
      }
    }
    image.end();

    // Add to memory
    if (newResults != null && newResults.size() > 0) {
      MemoryPeakResults.addResults(newResults);
    }

    return image.size();
  }

  private boolean showDialog() {
    settings = Settings.load();

    ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);

    // Build a list of all images with a region ROI
    final List<String> titles = new LinkedList<>();
    for (final int imageId : ImageJUtils.getIdList()) {
      final ImagePlus imp = WindowManager.getImage(imageId);
      if (imp != null && imp.getRoi() != null && imp.getRoi().isArea()) {
        titles.add(imp.getTitle());
      }
    }

    gd.addMessage("Show an image using the localisation density");

    ResultsManager.addInput(gd, settings.inputOption, InputSource.MEMORY);

    gd.addNumericField("Radius", settings.radius, 3);
    if (!titles.isEmpty()) {
      gd.addCheckbox((titles.size() == 1) ? "Use_ROI" : "Choose_ROI", settings.chooseRoi);
    }
    gd.addCheckbox("Adjust_for_border", settings.adjustForBorder);
    gd.addSlider("Image_Scale", 1, 15, settings.imageScale);
    gd.addCheckbox("Cumulative_image", settings.cumulativeImage);

    gd.addCheckbox("Use_square_approx", settings.useSquareApproximation);
    gd.addNumericField("Square_resolution", settings.resolution, 0);
    gd.addChoice("Score", Settings.scoreMethods, settings.scoreMethodIndex);

    gd.addMessage("Filter localisations using the L-score / Relative density.\n"
        + "Filtered results will be added to memory:");
    gd.addCheckbox("Filter_localisations", settings.filterLocalisations);
    gd.addNumericField("Filter_threshold", settings.filterThreshold, 2);

    gd.addCheckbox("Compute_Ripleys_L_plot", settings.computeRipleysPlot);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.inputOption = ResultsManager.getInputSource(gd);

    settings.radius = (float) gd.getNextNumber();
    if (!titles.isEmpty()) {
      settings.chooseRoi = gd.getNextBoolean();
    }
    settings.adjustForBorder = gd.getNextBoolean();
    settings.imageScale = (int) gd.getNextNumber();
    settings.cumulativeImage = gd.getNextBoolean();

    settings.useSquareApproximation = gd.getNextBoolean();
    settings.resolution = (int) gd.getNextNumber();
    settings.scoreMethodIndex = gd.getNextChoiceIndex();

    settings.filterLocalisations = gd.getNextBoolean();
    settings.filterThreshold = gd.getNextNumber();

    settings.computeRipleysPlot = gd.getNextBoolean();

    settings.save();

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Radius", settings.radius);
      ParameterUtils.isAboveZero("Image scale", settings.imageScale);
      ParameterUtils.isAboveZero("Resolution", settings.resolution);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    if (!titles.isEmpty() && settings.chooseRoi) {
      if (titles.size() == 1) {
        settings.roiImage = titles.get(0);
        Recorder.recordOption("Image", settings.roiImage);
      } else {
        final String[] items = titles.toArray(new String[titles.size()]);
        gd = new ExtendedGenericDialog(TITLE);
        gd.addMessage("Select the source image for the ROI");
        gd.addChoice("Image", items, settings.roiImage);
        gd.showDialog();
        if (gd.wasCanceled()) {
          return false;
        }
        settings.roiImage = gd.getNextChoice();
      }
      final ImagePlus imp = WindowManager.getImage(settings.roiImage);

      roiBounds = imp.getRoi().getBounds();
      roiImageWidth = imp.getWidth();
      roiImageHeight = imp.getHeight();
    } else {
      roiBounds = null;
    }

    return true;
  }

  /**
   * Compute the Ripley's L-function for user selected radii and show it on a plot.
   *
   * @param results the results
   */
  private void computeRipleysPlot(MemoryPeakResults results) {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addMessage("Compute Ripley's L(r) - r plot");
    gd.addNumericField("Min_radius", settings.minR, 2);
    gd.addNumericField("Max_radius", settings.maxR, 2);
    gd.addNumericField("Increment", settings.incrementR, 2);
    gd.addCheckbox("Confidence_intervals", settings.confidenceIntervals);

    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    settings.minR = gd.getNextNumber();
    settings.maxR = gd.getNextNumber();
    settings.incrementR = gd.getNextNumber();
    settings.confidenceIntervals = gd.getNextBoolean();

    if (settings.minR > settings.maxR || settings.incrementR < 0 || gd.invalidNumber()) {
      IJ.error(TITLE, "Invalid radius parameters");
      return;
    }

    final DensityManager dm = createDensityManager(results);
    final double[][] values = calculateLScores(dm);

    // 99% confidence intervals
    final int iterations = (settings.confidenceIntervals) ? 99 : 0;
    double[] upper = null;
    double[] lower = null;
    final Rectangle bounds = results.getBounds();
    // Use a uniform distribution for the coordinates
    final HaltonSequenceGenerator dist = new HaltonSequenceGenerator(2);
    dist.skipTo(SeedFactory.createInt());
    for (int i = 0; i < iterations; i++) {
      IJ.showProgress(i, iterations);
      IJ.showStatus(String.format("L-score confidence interval %d / %d", i + 1, iterations));

      // Randomise coordinates
      final float[] x = new float[results.size()];
      final float[] y = new float[x.length];
      for (int j = x.length; j-- > 0;) {
        final double[] d = dist.nextVector();
        x[j] = (float) (d[0] * bounds.width);
        y[j] = (float) (d[1] * bounds.height);
      }
      final double[][] values2 = calculateLScores(new DensityManager(x, y, bounds));
      if (upper == null || lower == null) {
        upper = values2[1];
        lower = upper.clone();
      } else {
        for (int m = upper.length; m-- > 0;) {
          if (upper[m] < values2[1][m]) {
            upper[m] = values2[1][m];
          }
          if (lower[m] > values2[1][m]) {
            lower[m] = values2[1][m];
          }
        }
      }
    }

    final String title = results.getName() + " Ripley's (L(r) - r) / r";
    final Plot2 plot = new Plot2(title, "Radius", "(L(r) - r) / r", values[0], values[1]);
    // Get the limits
    double yMin = min(0, values[1]);
    double yMax = max(0, values[1]);
    if (iterations > 0) {
      yMin = min(yMin, lower);
      yMax = max(yMax, upper);
    }
    plot.setLimits(0, values[0][values[0].length - 1], yMin, yMax);
    if (iterations > 0) {
      plot.setColor(Color.BLUE);
      plot.addPoints(values[0], upper, 1);
      plot.setColor(Color.RED);
      plot.addPoints(values[0], lower, 1);
      plot.setColor(Color.BLACK);
    }
    ImageJUtils.display(title, plot);
  }

  private static double min(double min, double[] data) {
    for (final double d : data) {
      if (min > d) {
        min = d;
      }
    }
    return min;
  }

  private static double max(double max, double[] data) {
    for (final double d : data) {
      if (max < d) {
        max = d;
      }
    }
    return max;
  }

  private double[][] calculateLScores(DensityManager dm) {
    final TDoubleArrayList x = new TDoubleArrayList();
    final TDoubleArrayList y = new TDoubleArrayList();
    x.add(0.0);
    y.add(0.0);

    for (double r = settings.minR; r < settings.maxR; r += settings.incrementR) {
      final double l = dm.ripleysLFunction(r);
      x.add(r);
      final double score = (r > 0) ? (l - r) / r : 0;
      y.add(score);
    }

    final double[][] values = new double[2][];
    values[0] = x.toArray();
    values[1] = y.toArray();
    return values;
  }
}
