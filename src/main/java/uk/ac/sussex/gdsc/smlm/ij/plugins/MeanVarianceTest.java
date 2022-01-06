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
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.plugin.PlugIn;
import ij.text.TextWindow;
import java.awt.Color;
import java.awt.Frame;
import java.awt.Point;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.atomic.AtomicReference;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.apache.commons.math3.util.MathArrays;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot.HistogramPlotBuilder;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.SeriesOpener;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.SortUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;

/**
 * Opens a folder of images and computes a Mean-Variance Test.
 *
 * <p>Each image must contain a 2-slice stack of white light images. The image filename must contain
 * the exposure time separated by whitespace, e.g. 'MVT 30.tif' for 30 milliseconds.
 *
 * <p>Gain calculations for standard CCD and EM-CCD cameras are based on the paper: Hirsch M,
 * Wareham RJ, Martin-Fernandez ML, Hobson MP, Rolfe DJ (2013) A Stochastic Model for Electron
 * Multiplication Charge-Coupled Devices – From Theory to Practice. PLoS ONE 8(1): e53671.
 * doi:10.1371/journal.pone.0053671
 */
public class MeanVarianceTest implements PlugIn {
  private static final String TITLE = "Mean Variance Test";

  private int exposureCounter;
  private boolean singleImage;

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    double cameraGain;
    double bias;
    boolean showTable;
    boolean showCharts;

    Settings() {
      // Set defaults
      bias = 500;
      showTable = true;
      showCharts = true;
    }

    Settings(Settings source) {
      cameraGain = source.cameraGain;
      bias = source.bias;
      showTable = source.showTable;
      showCharts = source.showCharts;
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
     * Save the settings. This can be called only once as it saves via a reference.
     */
    void save() {
      lastSettings.set(this);
    }
  }

  private static class PairSample {
    int slice1;
    int slice2;
    double mean1;
    double mean2;
    double variance;

    public PairSample(int slice1, int slice2, double mean1, double mean2, double variance) {
      this.slice1 = slice1;
      this.slice2 = slice2;
      this.mean1 = mean1;
      this.mean2 = mean2;
      this.variance = variance;
    }

    public double getMean() {
      return (mean1 + mean2) * 0.5;
    }
  }

  private class ImageSample {
    String title;
    float[][] slices;
    double exposure;
    double[] means;
    List<PairSample> samples;

    public ImageSample(ImagePlus imp, double start, double end) {
      // Check stack has two slices
      if (imp.getStackSize() < 2) {
        throw new IllegalArgumentException("Image must have at least 2-slices: " + imp.getTitle());
      }

      // Count all the valid input images
      exposureCounter++;

      // Extract the exposure time
      exposure = -1;
      final String[] tokens = imp.getTitle().split("[ .]");
      for (final String token : tokens) {
        try {
          exposure = Double.parseDouble(token);
          if (exposure >= 0) {
            break;
          }
        } catch (final NumberFormatException ex) {
          // Ignore
        }
      }

      if (exposure < 0) {
        // If no exposure was found: assume exposure 0 for first input image otherwise set an
        // arbitrary exposure
        exposure = (exposureCounter == 1) ? 0 : 9999;
      }

      title = imp.getTitle();

      // Get all the pixels into a float stack.
      // Look for saturated pixels that will invalidate the test.
      final int size = imp.getStackSize();
      slices = new float[size][];
      final float saturated = getSaturation(imp);
      final ImageStack stack = imp.getImageStack();
      final double step = (end - start) / size;
      for (int slice = 1, c = 0; slice <= size; slice++) {
        if (c++ % 16 == 0) {
          IJ.showProgress(start + c * step);
        }
        final float[] thisSlice =
            slices[slice - 1] = (float[]) stack.getProcessor(slice).toFloat(0, null).getPixels();
        checkSaturation(slice, thisSlice, saturated);
      }
    }

    private float getSaturation(ImagePlus imp) {
      switch (imp.getBitDepth()) {
        case 8:
        case 24:
          return 255f;
        case 16:
          return 65535f;
        case 32:
          // float images cannot be saturated
          return Float.NaN;
        default:
          throw new IllegalArgumentException(
              "Cannot determine saturation level for image: " + imp.getTitle());
      }
    }

    private void checkSaturation(int slice, float[] data, float saturated) {
      if (Float.isNaN(saturated)) {
        // float images cannot be saturated
        return;
      }
      for (final float f : data) {
        if (f >= saturated) {
          throw new IllegalArgumentException(
              "Image " + title + " has saturated pixels in slice: " + slice);
        }
      }
    }

    public void compute(boolean consecutive, double start, double end) {
      final int size = slices.length;
      final int nSamples = (consecutive) ? size - 1 : ((size - 1) * size) / 2;
      samples = new ArrayList<>(nSamples);

      // Cache data
      means = new double[size];
      for (int slice1 = 0; slice1 < size; slice1++) {
        means[slice1] = Statistics.create(slices[slice1]).getMean();
      }

      // Compute mean and variance.
      // See http://www.photometrics.com/resources/whitepapers/mean-variance.php
      final double step = (end - start) / nSamples;
      for (int slice1 = 0, c = 0; slice1 < size; slice1++) {
        final float[] data1 = slices[slice1];
        for (int slice2 = slice1 + 1; slice2 < size; slice2++) {
          if (c++ % 16 == 0) {
            IJ.showProgress(start + c * step);
          }
          final float[] data2 = slices[slice2];
          final Statistics s = new Statistics();
          for (int i = 0; i < data1.length; i++) {
            s.add(data1[i] - data2[i]);
          }
          final double variance = s.getVariance() / 2.0;
          samples
              .add(new PairSample(slice1 + 1, slice2 + 1, means[slice1], means[slice2], variance));

          if (consecutive) {
            break;
          }
        }
        slices[slice1] = null; // Allow garbage collection
      }
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    settings = Settings.load();
    settings.save();

    String helpKey = "mean-variance-test";
    if (ImageJUtils.isExtraOptions()) {
      final ImagePlus imp = WindowManager.getCurrentImage();
      if (imp.getStackSize() > 1) {
        final GenericDialog gd = new GenericDialog(TITLE);
        gd.addMessage("Perform single image analysis on the current image?");
        gd.addNumericField("Bias", settings.bias, 0);
        gd.addHelp(HelpUrls.getUrl(helpKey));
        gd.showDialog();
        if (gd.wasCanceled()) {
          return;
        }
        singleImage = true;
        settings.bias = Math.abs(gd.getNextNumber());
      } else {
        IJ.error(TITLE, "Single-image mode requires a stack");
        return;
      }
    }

    List<ImageSample> images;
    String inputDirectory = "";
    if (singleImage) {
      IJ.showStatus("Loading images...");
      images = getImages();
      if (images.size() == 0) {
        IJ.error(TITLE, "Not enough images for analysis");
        return;
      }
    } else {
      inputDirectory = IJ.getDirectory("Select image series ...");
      if (inputDirectory == null) {
        return;
      }

      final SeriesOpener series = new SeriesOpener(inputDirectory);
      series.setVariableSize(true);
      if (series.getNumberOfImages() < 3) {
        IJ.error(TITLE, "Not enough images in the selected directory");
        return;
      }
      if (!IJ.showMessageWithCancel(TITLE, String.format("Analyse %d images, first image:\n%s",
          series.getNumberOfImages(), series.getImageList()[0]))) {
        return;
      }

      IJ.showStatus("Loading images");
      images = getImages(series);

      if (images.size() < 3) {
        IJ.error(TITLE, "Not enough images for analysis");
        return;
      }
      if (images.get(0).exposure != 0) {
        IJ.error(TITLE, "First image in series must have exposure 0 (Bias image)");
        return;
      }
    }

    final boolean emMode = (arg != null && arg.contains("em"));

    GenericDialog gd = new GenericDialog(TITLE);
    gd.addMessage("Set the output options:");
    gd.addCheckbox("Show_table", settings.showTable);
    gd.addCheckbox("Show_charts", settings.showCharts);
    if (emMode) {
      // Ask the user for the camera gain ...
      gd.addMessage("Estimating the EM-gain requires the camera gain without EM readout enabled");
      gd.addNumericField("Camera_gain (Count/e-)", settings.cameraGain, 4);
    }
    if (emMode) {
      helpKey += "-em-ccd";
    }
    gd.addHelp(HelpUrls.getUrl(helpKey));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    settings.showTable = gd.getNextBoolean();
    settings.showCharts = gd.getNextBoolean();
    if (emMode) {
      settings.cameraGain = gd.getNextNumber();
    }

    IJ.showStatus("Computing mean & variance");
    final double nImages = images.size();
    for (int i = 0; i < images.size(); i++) {
      IJ.showStatus(String.format("Computing mean & variance %d/%d", i + 1, images.size()));
      images.get(i).compute(singleImage, i / nImages, (i + 1) / nImages);
    }
    IJ.showProgress(1);

    IJ.showStatus("Computing results");

    // Allow user to input multiple bias images
    int start = 0;
    final Statistics biasStats = new Statistics();
    final Statistics noiseStats = new Statistics();
    final double bias;
    if (singleImage) {
      bias = settings.bias;
    } else {
      while (start < images.size()) {
        final ImageSample sample = images.get(start);
        if (sample.exposure == 0) {
          biasStats.add(sample.means);
          for (final PairSample pair : sample.samples) {
            noiseStats.add(pair.variance);
          }
          start++;
        } else {
          break;
        }
      }
      bias = biasStats.getMean();
    }

    // Get the mean-variance data
    int total = 0;
    for (int i = start; i < images.size(); i++) {
      total += images.get(i).samples.size();
    }

    if (settings.showTable && total > 2000) {
      gd = new GenericDialog(TITLE);
      gd.addMessage(
          "Table output requires " + total + " entries.\n \nYou may want to disable the table.");
      gd.addCheckbox("Show_table", settings.showTable);
      gd.showDialog();
      if (gd.wasCanceled()) {
        return;
      }
      settings.showTable = gd.getNextBoolean();
    }

    final TextWindow results = (settings.showTable) ? createResultsWindow() : null;

    double[] mean = new double[total];
    double[] variance = new double[mean.length];
    final Statistics gainStats = (singleImage) ? new StoredDataStatistics(total) : new Statistics();
    final WeightedObservedPoints obs = new WeightedObservedPoints();
    for (int i = (singleImage) ? 0 : start, j = 0; i < images.size(); i++) {
      final StringBuilder sb = (settings.showTable) ? new StringBuilder() : null;
      final ImageSample sample = images.get(i);
      for (final PairSample pair : sample.samples) {
        if (j % 16 == 0) {
          IJ.showProgress(j, total);
        }

        mean[j] = pair.getMean();
        variance[j] = pair.variance;
        // Gain is in Count / e
        double gain = variance[j] / (mean[j] - bias);
        gainStats.add(gain);
        obs.add(mean[j], variance[j]);

        if (emMode) {
          gain /= (2 * settings.cameraGain);
        }

        if (sb != null) {
          sb.append(sample.title).append('\t');
          sb.append(sample.exposure).append('\t');
          sb.append(pair.slice1).append('\t');
          sb.append(pair.slice2).append('\t');
          sb.append(IJ.d2s(pair.mean1, 2)).append('\t');
          sb.append(IJ.d2s(pair.mean2, 2)).append('\t');
          sb.append(IJ.d2s(mean[j], 2)).append('\t');
          sb.append(IJ.d2s(variance[j], 2)).append('\t');
          sb.append(MathUtils.rounded(gain, 4)).append("\n");
        }
        j++;
      }
      if (results != null && sb != null) {
        results.append(sb.toString());
      }
    }
    IJ.showProgress(1);

    if (singleImage) {
      StoredDataStatistics stats = (StoredDataStatistics) gainStats;
      ImageJUtils.log(TITLE);
      if (emMode) {
        final double[] values = stats.getValues();
        MathArrays.scaleInPlace(0.5, values);
        stats = StoredDataStatistics.create(values);
      }

      if (settings.showCharts) {
        // Plot the gain over time
        final String title = TITLE + " Gain vs Frame";
        final Plot plot = new Plot(title, "Slice", "Gain");
        plot.addPoints(SimpleArrayUtils.newArray(gainStats.getN(), 1, 1.0), stats.getValues(),
            Plot.LINE);
        final PlotWindow pw = ImageJUtils.display(title, plot);

        // Show a histogram
        final String label = String.format("Mean = %s, Median = %s",
            MathUtils.rounded(stats.getMean()), MathUtils.rounded(stats.getMedian()));
        final WindowOrganiser wo = new WindowOrganiser();
        final PlotWindow pw2 = new HistogramPlotBuilder(TITLE, stats, "Gain")
            .setRemoveOutliersOption(1).setPlotLabel(label).show(wo);
        if (wo.isNotEmpty()) {
          final Point point = pw.getLocation();
          point.y += pw.getHeight();
          pw2.setLocation(point);
        }
      }

      ImageJUtils.log("Single-image mode: %s camera", (emMode) ? "EM-CCD" : "Standard");
      final double gain = stats.getMedian();

      if (emMode) {
        final double totalGain = gain;
        final double emGain = totalGain / settings.cameraGain;
        ImageJUtils.log("  Gain = 1 / %s (Count/e-)", MathUtils.rounded(settings.cameraGain, 4));
        ImageJUtils.log("  EM-Gain = %s", MathUtils.rounded(emGain, 4));
        ImageJUtils.log("  Total Gain = %s (Count/e-)", MathUtils.rounded(totalGain, 4));
      } else {
        settings.cameraGain = gain;
        ImageJUtils.log("  Gain = 1 / %s (Count/e-)", MathUtils.rounded(settings.cameraGain, 4));
      }
    } else {
      IJ.showStatus("Computing fit");

      // Sort
      final int[] indices = rank(mean);
      mean = reorder(mean, indices);
      variance = reorder(variance, indices);

      // Compute optimal coefficients.
      final double[] init = {0, 1 / gainStats.getMean()}; // a - b x
      final PolynomialCurveFitter fitter = PolynomialCurveFitter.create(2).withStartPoint(init);
      final double[] best = fitter.fit(obs.toList());

      // Construct the polynomial that best fits the data.
      final PolynomialFunction fitted = new PolynomialFunction(best);

      if (settings.showCharts) {
        // Plot mean verses variance. Gradient is gain in Count/e.
        final String title = TITLE + " results";
        final Plot plot = new Plot(title, "Mean", "Variance");
        final double[] xlimits = MathUtils.limits(mean);
        final double[] ylimits = MathUtils.limits(variance);
        double xrange = (xlimits[1] - xlimits[0]) * 0.05;
        if (xrange == 0) {
          xrange = 0.05;
        }
        double yrange = (ylimits[1] - ylimits[0]) * 0.05;
        if (yrange == 0) {
          yrange = 0.05;
        }
        plot.setLimits(xlimits[0] - xrange, xlimits[1] + xrange, ylimits[0] - yrange,
            ylimits[1] + yrange);
        plot.setColor(Color.blue);
        plot.addPoints(mean, variance, Plot.CROSS);
        plot.setColor(Color.red);
        plot.addPoints(new double[] {mean[0], mean[mean.length - 1]},
            new double[] {fitted.value(mean[0]), fitted.value(mean[mean.length - 1])}, Plot.LINE);
        ImageJUtils.display(title, plot);
      }

      final double avBiasNoise = Math.sqrt(noiseStats.getMean());

      ImageJUtils.log(TITLE);
      ImageJUtils.log("  Directory = %s", inputDirectory);
      ImageJUtils.log("  Bias = %s +/- %s (Count)", MathUtils.rounded(bias, 4),
          MathUtils.rounded(avBiasNoise, 4));
      ImageJUtils.log("  Variance = %s + %s * mean", MathUtils.rounded(best[0], 4),
          MathUtils.rounded(best[1], 4));
      if (emMode) {
        // The gradient is the observed gain of the noise.
        // In an EM-CCD there is a noise factor of 2.
        // Q. Is this true for a correct noise factor calibration:
        // double noiseFactor = (Read Noise EM-CCD) / (Read Noise CCD)

        // Em-gain is the observed gain divided by the noise factor multiplied by camera gain
        final double emGain = best[1] / (2 * settings.cameraGain);

        // Compute total gain
        final double totalGain = emGain * settings.cameraGain;

        final double readNoise = avBiasNoise / settings.cameraGain;
        // Effective noise is standard deviation of the bias image divided by the total gain (in
        // Count/e-)
        final double readNoiseE = avBiasNoise / totalGain;
        ImageJUtils.log("  Read Noise = %s (e-) [%s (Count)]", MathUtils.rounded(readNoise, 4),
            MathUtils.rounded(avBiasNoise, 4));

        ImageJUtils.log("  Gain = 1 / %s (Count/e-)",
            MathUtils.rounded(1 / settings.cameraGain, 4));
        ImageJUtils.log("  EM-Gain = %s", MathUtils.rounded(emGain, 4));
        ImageJUtils.log("  Total Gain = %s (Count/e-)", MathUtils.rounded(totalGain, 4));
        ImageJUtils.log("  Effective Read Noise = %s (e-) (Read Noise/Total Gain)",
            MathUtils.rounded(readNoiseE, 4));
      } else {
        // The gradient is the observed gain of the noise.
        settings.cameraGain = best[1];
        // Noise is standard deviation of the bias image divided by the gain (in Count/e-)
        final double readNoise = avBiasNoise / settings.cameraGain;
        ImageJUtils.log("  Read Noise = %s (e-) [%s (Count)]", MathUtils.rounded(readNoise, 4),
            MathUtils.rounded(avBiasNoise, 4));

        ImageJUtils.log("  Gain = 1 / %s (Count/e-)",
            MathUtils.rounded(1 / settings.cameraGain, 4));
      }
    }
    IJ.showStatus("");
  }

  private static TextWindow createResultsWindow() {
    final Frame f = WindowManager.getFrame(TITLE);
    if (f instanceof TextWindow) {
      return (TextWindow) f;
    }
    return new TextWindow(TITLE,
        "Image\tExposure\tSlice1\tSlice2\tMean1\tMean2\tMean\tVariance\tGain", "", 800, 500);
  }

  private List<ImageSample> getImages(SeriesOpener series) {
    final double nImages = series.getNumberOfImages();
    final List<ImageSample> images = new ArrayList<>((int) nImages);
    ImagePlus imp = series.nextImage();
    int count = 0;
    while (imp != null) {
      try {
        images.add(new ImageSample(imp, count / nImages, (count + 1) / nImages));
      } catch (final IllegalArgumentException ex) {
        ImageJUtils.log(ex.getMessage());
      }
      count++;
      imp.close();
      imp = series.nextImage();
    }
    IJ.showProgress(1);
    // Sort to ensure all 0 exposure images are first, the remaining order is arbitrary
    Collections.sort(images, (o1, o2) -> Double.compare(o1.exposure, o2.exposure));
    return images;
  }

  private List<ImageSample> getImages() {
    final List<ImageSample> images = new ArrayList<>(1);
    final ImagePlus imp = WindowManager.getCurrentImage();
    if (imp != null) {
      try {
        images.add(new ImageSample(imp, 0, 1));
      } catch (final IllegalArgumentException ex) {
        ImageJUtils.log(ex.getMessage());
      }
    }
    IJ.showProgress(1);
    return images;
  }

  /**
   * Returns a sorted list of indices of the specified double array.
   *
   * @param values the values
   * @return the indices
   */
  public static int[] rank(double[] values) {
    final int[] indices = SimpleArrayUtils.natural(values.length);
    SortUtils.sortData(indices, values, false, false);
    return indices;
  }

  private static double[] reorder(double[] data, int[] indices) {
    final double[] array = new double[indices.length];
    for (int i = 0; i < array.length; i++) {
      array[i] = data[indices[i]];
    }
    return array;
  }
}
