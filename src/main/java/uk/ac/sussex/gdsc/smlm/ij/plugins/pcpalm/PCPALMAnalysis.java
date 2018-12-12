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

package uk.ac.sussex.gdsc.smlm.ij.plugins.pcpalm;

import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.Plot2;
import uk.ac.sussex.gdsc.core.ij.process.Fht;
import uk.ac.sussex.gdsc.core.utils.ImageWindow;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.smlm.ij.plugins.About;
import uk.ac.sussex.gdsc.smlm.ij.plugins.Parameters;
import uk.ac.sussex.gdsc.smlm.ij.plugins.SMLMUsageTracker;
import uk.ac.sussex.gdsc.smlm.model.MaskDistribution;

import com.thoughtworks.xstream.XStream;
import com.thoughtworks.xstream.XStreamException;
import com.thoughtworks.xstream.io.xml.DomDriver;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.Roi;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.text.TextWindow;

import org.apache.commons.math3.util.FastMath;
import org.jtransforms.fft.DoubleFFT_2D;
import org.jtransforms.fft.FloatFFT_2D;

import java.awt.Color;
import java.awt.Frame;
import java.awt.Rectangle;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

/**
 * Use the PC-PALM protocol to analyse a set of molecules to produce a correlation curve.
 *
 * <p>Frequency domain analysis see Sengupta, et al (2013). Quantifying spatial resolution in
 * point-localisation superresolution images using pair correlation analysis. Nature Protocols 8,
 * pp345-354.
 *
 * <p>Spatial domain analysis see Puchnar, et al (2013). Counting molecules in single organelles
 * with superresolution microscopy allows tracking of the endosome maturation trajectory. PNAS.
 * doi:10.1073/pnas.1309676110
 */
public class PCPALMAnalysis implements PlugInFilter {
  /** The title. */
  static String TITLE = "PC-PALM Analysis";

  private static String resultsDirectory = "";
  private static double correlationDistance = 800; // nm
  private static double correlationInterval = 20; // nm
  private static boolean binaryImage;

  /** The blinking rate. */
  static double blinkingRate = -1;
  private static double copiedBlinkingRate = -1;
  private static double nmPerPixel = -1;
  private static double copiedNmPerPixel = -1;
  private static boolean showErrorBars;
  private static boolean applyWindow;
  private static boolean showHighResolutionImage;
  private static boolean showCorrelationImages;
  private static boolean useBorder = true;

  // Limits for the molecules constructed from the input ROI
  private double minx;
  private double miny;
  private double maxx;
  private double maxy;

  private double area;
  private double weightedArea;
  private double weightedAreaInPx;
  private int areaInPx;
  private int noOfMolecules;
  private double uniquePoints;

  // Cache the ImageWindow for faster repeat processing
  private static ImageWindow imageWindow = new ImageWindow();

  // Used for the results table
  private static TextWindow resultsTable;

  /** The results. */
  static ArrayList<CorrelationResult> results = new ArrayList<>();

  private boolean spatialDomain;

  /** Area of the region cropped from the PCPALM Molecules list. */
  double croppedArea;

  /** {@inheritDoc} */
  @Override
  public int setup(String arg, ImagePlus imp) {
    SMLMUsageTracker.recordPlugin(this.getClass(), arg);

    if ("save".equalsIgnoreCase(arg)) {
      return saveResults();
    }
    if ("load".equalsIgnoreCase(arg)) {
      return loadResults();
    }

    spatialDomain = "spatial".equalsIgnoreCase(arg);

    if (imp == null || (!spatialDomain && (imp.getRoi() == null || !imp.getRoi().isArea()))) {
      error("Require an input image with an area ROI.\n"
          + "Please create a binary molecule image using " + PCPALMMolecules.TITLE);
      return DONE;
    }
    if (PCPALMMolecules.molecules == null || PCPALMMolecules.molecules.size() < 2) {
      error("Require a set of molecules for analysis.\n"
          + "Please create a binary molecule image using " + PCPALMMolecules.TITLE);
      return DONE;
    }

    if (!showDialog()) {
      return DONE;
    }

    PCPALMMolecules.logSpacer();
    log(TITLE);
    PCPALMMolecules.logSpacer();

    final ArrayList<Molecule> molecules = cropToRoi(imp);
    if (molecules.size() < 2) {
      error("No results within the crop region");
      return DONE;
    }

    log("Using %d molecules", molecules.size());
    final long start = System.currentTimeMillis();

    analyse(molecules);

    final double seconds = (System.currentTimeMillis() - start) / 1000.0;
    final String msg = TITLE + " complete : " + seconds + "s";
    IJ.showStatus(msg);
    log(msg);
    return DONE;
  }

  /**
   * Show a directory selection dialog for the results directory.
   *
   * @return True if a directory was selected
   */
  private static boolean getDirectory() {
    resultsDirectory = ImageJUtils.getDirectory("Results_directory", resultsDirectory);
    return resultsDirectory != null;
  }

  /**
   * Save all the results to a directory.
   *
   * @return {@link PlugInFilter#DONE }
   */
  private static int saveResults() {
    if (results.isEmpty()) {
      error("No results in memory");
    } else if (getDirectory()) {
      final XStream xs = new XStream(new DomDriver());
      XStream.setupDefaultSecurity(xs); // to be removed after 1.5
      xs.allowTypes(new Class[] {CorrelationResult.class});
      for (final CorrelationResult result : results) {
        saveResult(xs, result);
      }
    }
    return DONE;
  }

  private static void saveResult(XStream xs, CorrelationResult result) {
    final String outputFilename = String.format("%s/%s.%d.xml", resultsDirectory,
        (result.spatialDomain) ? "Spatial" : "Frequency", result.id);
    try (FileOutputStream fs = new FileOutputStream(outputFilename)) {
      xs.toXML(result, fs);
    } catch (final XStreamException ex) {
      // ex.printStackTrace();
      IJ.log("Failed to save correlation result to file: " + outputFilename);
    } catch (final Exception ex) {
      IJ.log("Failed to save correlation result to file: " + outputFilename);
    }
  }

  /**
   * Load all the results from a directory. File must have the XML suffix.
   *
   * @return {@link PlugInFilter#DONE }
   */
  private static int loadResults() {
    if (getDirectory()) {
      final File[] fileList = (new File(resultsDirectory)).listFiles(new FilenameFilter() {
        @Override
        public boolean accept(File arg0, String arg1) {
          return arg1.endsWith("xml");
        }
      });
      if (fileList == null) {
        return DONE;
      }

      int count = 0;
      for (int i = 0; i < fileList.length; i++) {
        final XStream xs = new XStream(new DomDriver());
        XStream.setupDefaultSecurity(xs); // to be removed after 1.5
        xs.allowTypes(new Class[] {CorrelationResult.class});
        if (fileList[i].isFile()) {
          if (loadResult(xs, fileList[i].getPath())) {
            count++;
          }
        }
      }
      if (count > 0) {
        Collections.sort(results);
      }
      log("Loaded %d results", count);
    }
    return DONE;
  }

  private static boolean loadResult(XStream xs, String path) {
    try (FileInputStream fs = new FileInputStream(path)) {
      CorrelationResult result = (CorrelationResult) xs.fromXML(fs);
      // Replace a result with the same id
      for (int i = 0; i < results.size(); i++) {
        if (results.get(i).id == result.id) {
          results.set(i, result);
          result = null;
          break;
        }
      }
      // Add to the results if we did not replace any
      if (result != null) {
        results.add(result);
      }
      return true;
    } catch (final ClassCastException ex) {
      // ex.printStackTrace();
      IJ.log("Failed to load correlation result from file: " + path);
    } catch (final XStreamException ex) {
      // ex.printStackTrace();
      IJ.log("Failed to load correlation result from file: " + path);
    } catch (final Exception ex) {
      IJ.log("Failed to load correlation result from file: " + path);
    }
    return false;
  }

  private static void error(String message) {
    log("ERROR : " + message);
    IJ.error(TITLE, message);
  }

  /** {@inheritDoc} */
  @Override
  public void run(ImageProcessor ip) {
    // Nothing to do
  }

  private boolean showDialog() {
    final GenericDialog gd = new GenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);

    if (blinkingRate < 1 || copiedBlinkingRate != PCPALMMolecules.blinkingRate) {
      copiedBlinkingRate = blinkingRate = PCPALMMolecules.blinkingRate;
    }

    if (nmPerPixel < 1 || copiedNmPerPixel != PCPALMMolecules.nmPerPixel) {
      copiedNmPerPixel = nmPerPixel = PCPALMMolecules.nmPerPixel;
    }

    gd.addMessage("Analyse clusters using Pair Correlation.");

    gd.addNumericField("Correlation_distance (nm)", correlationDistance, 0);
    if (!spatialDomain) {
      gd.addMessage("-=- Frequency domain analysis -=-");
      gd.addCheckbox("Binary_image", binaryImage);
      gd.addNumericField("Blinking_rate", blinkingRate, 2);
      gd.addNumericField("nm_per_pixel", nmPerPixel, 2);
      gd.addCheckbox("Show_error_bars", showErrorBars);
      gd.addCheckbox("Apply_window", applyWindow);
      gd.addCheckbox("Show_high_res_image", showHighResolutionImage);
      gd.addCheckbox("Show_correlation_images", showCorrelationImages);
    } else {
      gd.addMessage("-=- Spatial domain analysis -=-");
      gd.addCheckbox("Use_border", useBorder);
      gd.addNumericField("Correlation_interval (nm)", correlationInterval, 0);
    }

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    correlationDistance = gd.getNextNumber();
    if (!spatialDomain) {
      binaryImage = gd.getNextBoolean();
      blinkingRate = gd.getNextNumber();
      nmPerPixel = gd.getNextNumber();
      showErrorBars = gd.getNextBoolean();
      applyWindow = gd.getNextBoolean();
      showHighResolutionImage = gd.getNextBoolean();
      showCorrelationImages = gd.getNextBoolean();
    } else {
      useBorder = gd.getNextBoolean();
      correlationInterval = gd.getNextNumber();
    }

    // Check arguments
    try {
      if (!spatialDomain) {
        Parameters.isAbove("Correlation distance", correlationDistance, 1);
        Parameters.isEqualOrAbove("Blinking_rate", blinkingRate, 1);
        Parameters.isAboveZero("nm per pixel", nmPerPixel);
      } else {
        Parameters.isAboveZero("Correlation interval", correlationInterval);
      }
    } catch (final IllegalArgumentException ex) {
      error(ex.getMessage());
      return false;
    }

    return true;
  }

  /**
   * Extract all the PC-PALM molecules that are within the image ROI region. The coordinates bounds
   * are converted using relative scaling to the limits of the PC-PALM molecules. If a
   * non-rectangular ROI is used then a mask is extracted and used for the crop. If no image is
   * provided then the full set of molecules is returned.
   *
   * <p>Set the area property to the region covered by the molecules.
   *
   * @param imp the image
   * @return the array list
   */
  ArrayList<Molecule> cropToRoi(ImagePlus imp) {
    croppedArea = PCPALMMolecules.area;
    if (PCPALMMolecules.molecules == null || PCPALMMolecules.molecules.isEmpty()) {
      return PCPALMMolecules.molecules;
    }

    final double pcw = PCPALMMolecules.maxx - PCPALMMolecules.minx;
    final double pch = PCPALMMolecules.maxy - PCPALMMolecules.miny;

    final Roi roi = (imp == null) ? null : imp.getRoi();
    if (roi == null || !roi.isArea()) {
      log("Roi = %s nm x %s nm = %s um^2", MathUtils.rounded(pcw, 3), MathUtils.rounded(pch, 3),
          MathUtils.rounded(croppedArea, 3));
      minx = PCPALMMolecules.minx;
      maxx = PCPALMMolecules.maxx;
      miny = PCPALMMolecules.miny;
      maxy = PCPALMMolecules.maxy;
      return PCPALMMolecules.molecules;
    }

    // To avoid null pointer warnings
    if (imp == null) {
      throw new NullPointerException("image is null");
    }

    final int w = imp.getWidth();
    final int h = imp.getHeight();

    final Rectangle bounds = roi.getBounds();

    // Construct relative coordinates
    minx = PCPALMMolecules.minx + pcw * ((double) bounds.x / w);
    miny = PCPALMMolecules.miny + pch * ((double) bounds.y / h);
    maxx = PCPALMMolecules.minx + pcw * ((double) (bounds.x + bounds.width) / w);
    maxy = PCPALMMolecules.miny + pch * ((double) (bounds.y + bounds.height) / h);

    final double roix = maxx - minx;
    final double roiy = maxy - miny;

    final ArrayList<Molecule> newMolecules = new ArrayList<>(PCPALMMolecules.molecules.size() / 2);

    // Support non-rectangular ROIs
    if (roi.getMask() != null) {
      final MaskDistribution dist = createMaskDistribution(roi.getMask(), roix, roiy);

      final double fraction =
          (double) dist.getSize() / (roi.getMask().getWidth() * roi.getMask().getHeight());
      log("Roi is a mask of %d pixels", dist.getSize());
      croppedArea = fraction * roix * roiy / 1e6;
      log("Roi area %s x %s nm x %s nm = %s um^2", MathUtils.rounded(fraction),
          MathUtils.rounded(roix, 3), MathUtils.rounded(roiy, 3),
          MathUtils.rounded(croppedArea, 3));

      final double[] xyz = new double[3];
      // The mask is 0,0 in the centre so offset by this as well
      final double xoffset = minx + roix / 2;
      final double yoffset = miny + roiy / 2;
      for (final Molecule m : PCPALMMolecules.molecules) {
        xyz[0] = m.x - xoffset;
        xyz[1] = m.y - yoffset;
        if (dist.isWithinXy(xyz)) {
          newMolecules.add(m);
        }
      }
    } else {
      croppedArea = roix * roiy / 1e6;
      log("Roi = %s nm x %s nm = %s um^2", MathUtils.rounded(roix, 3), MathUtils.rounded(roiy, 3),
          MathUtils.rounded(croppedArea, 3));

      for (final Molecule m : PCPALMMolecules.molecules) {
        if (m.x < minx || m.x >= maxx || m.y < miny || m.y >= maxy) {
          continue;
        }
        newMolecules.add(m);
      }
    }

    return newMolecules;
  }

  private static MaskDistribution createMaskDistribution(ImageProcessor ip, double roix,
      double roiy) {
    // Calculate the scale of the mask
    final int w = ip.getWidth();
    final int h = ip.getHeight();
    final double scaleX = roix / w;
    final double scaleY = roiy / h;

    // Use an image for the distribution
    final int[] mask = extractMask(ip);
    return new MaskDistribution(mask, w, h, 0, scaleX, scaleY);
  }

  private static int[] extractMask(ImageProcessor ip) {
    final int[] mask = new int[ip.getPixelCount()];
    for (int i = 0; i < mask.length; i++) {
      mask[i] = ip.get(i);
    }
    return mask;
  }

  /**
   * Log a message to the IJ log window.
   *
   * @param format the format
   * @param args the args
   */
  private static void log(String format, Object... args) {
    IJ.log(String.format(format, args));
  }

  /**
   * Perform the PC Analysis.
   *
   * @param molecules the molecules
   */
  private void analyse(ArrayList<Molecule> molecules) {
    // Check if the plots are currently shown
    final String spatialPlotTitle = TITLE + " molecules/um^2";
    final String frequencyDomainTitle = TITLE + " g(r)";

    boolean noPlots;
    String topPlotTitle;

    final long start = System.nanoTime();
    if (spatialDomain) {
      // -----------------
      // Analysis in the spatial domain
      // -----------------

      log("---");
      log("Spatial domain analysis");
      log("Computing density histogram");

      // Compute all-vs-all distance matrix.
      // Create histogram of distances at different radii.
      final int nBins = (int) (correlationDistance / correlationInterval) + 1;
      final double maxDistance2 = correlationDistance * correlationDistance;
      final int[] H = new int[nBins];

      // TODO - Update this using a grid with a resolution of maxDistance to increase speed
      // by only comparing to neighbours within range.

      // An all-vs-all analysis does not account for a border.
      // A simple solution is to only process molecules within the border but compare them
      // to all molecules within the region. Thus every molecule has a complete circle of the max
      // radius around them to use:
      // ----------------------
      // | |
      // | -------------- |
      // | | Within | |
      // | | Border | |
      // | | | |
      // | -------------- |
      // | Region |
      // ----------------------
      // If the fraction of points within the correlation distance of the edge is low then this
      // will not make much difference.

      final double boundaryMinx = (useBorder) ? minx + correlationDistance : minx;
      final double boundaryMaxx = (useBorder) ? maxx - correlationDistance : maxx;
      final double boundaryMiny = (useBorder) ? miny + correlationDistance : miny;
      final double boundaryMaxy = (useBorder) ? maxy - correlationDistance : maxy;

      int N = 0;
      if (boundaryMaxx <= boundaryMinx || boundaryMaxy <= boundaryMiny) {
        log("ERROR: 'Use border' option of %s nm is not possible: Width = %s nm, Height = %s nm",
            MathUtils.rounded(correlationDistance, 4), MathUtils.rounded(maxx - minx, 3),
            MathUtils.rounded(maxy - miny, 3));
        return;
      }

      for (int i = molecules.size(); i-- > 0;) {
        final Molecule m = molecules.get(i);
        // Optionally ignore molecules that are near the edge of the boundary
        if (useBorder && (m.x < boundaryMinx || m.x > boundaryMaxx || m.y < boundaryMiny
            || m.y > boundaryMaxy)) {
          continue;
        }
        N++;

        for (int j = molecules.size(); j-- > 0;) {
          if (i == j) {
            continue;
          }

          final double d = m.distance2(molecules.get(j));
          if (d < maxDistance2) {
            H[(int) (Math.sqrt(d) / correlationInterval)]++;
          }
        }
      }

      double[] r = new double[nBins + 1];
      for (int i = 0; i <= nBins; i++) {
        r[i] = i * correlationInterval;
      }
      double[] pcf = new double[nBins];
      if (N > 0) {
        // Note: Convert nm^2 to um^2
        final double N_pi = N * Math.PI / 1000000.0;
        for (int i = 0; i < nBins; i++) {
          // Pair-correlation is the count at the given distance divided by N and the area at
          // distance ri:
          // H(r_i) / (N x (pi x (r_i+1)^2 - pi x r_i^2))
          pcf[i] = H[i] / (N_pi * (r[i + 1] * r[i + 1] - r[i] * r[i]));
        }
      }

      // The final bin may be empty if the correlation interval was a factor of the correlation
      // distance
      if (pcf[pcf.length - 1] == 0) {
        r = Arrays.copyOf(r, nBins - 1);
        pcf = Arrays.copyOf(pcf, nBins - 1);
      } else {
        r = Arrays.copyOf(r, nBins);
      }

      final double[][] gr = new double[][] {r, pcf, null};

      final CorrelationResult result = new CorrelationResult(results.size() + 1,
          PCPALMMolecules.results.getSource(), boundaryMinx, boundaryMiny, boundaryMaxx,
          boundaryMaxy, N, correlationInterval, 0, false, gr, true);
      results.add(result);

      noPlots = WindowManager.getFrame(spatialPlotTitle) == null;
      topPlotTitle = frequencyDomainTitle;

      plotCorrelation(gr, 0, spatialPlotTitle, "molecules/um^2", true, false);
    } else {
      // -----------------
      // Image correlation in the Frequency Domain
      // -----------------
      log("Frequency domain analysis");

      // Create a binary image for the molecules

      ImageProcessor im = PCPALMMolecules.drawImage(molecules, minx, miny, maxx, maxy, nmPerPixel,
          true, binaryImage);
      log("Image scale = %.2f nm/pixel : %d x %d pixels", nmPerPixel, im.getWidth(),
          im.getHeight());

      // Apply a window function to the image to reduce FFT edge artifacts.
      if (applyWindow) {
        im = applyWindow(im, imageWindow);
      }

      if (showHighResolutionImage) {
        PCPALMMolecules.displayImage(PCPALMMolecules.results.getName() + " "
            + ((binaryImage) ? "Binary" : "Count") + " Image (high res)", im, nmPerPixel);
      }

      // Create weight image (including windowing)
      ImageProcessor w = createWeightImage(im, applyWindow);

      // Store the area of the image in um^2
      weightedAreaInPx = areaInPx = im.getWidth() * im.getHeight();
      if (applyWindow) {
        weightedAreaInPx *= w.getStatistics().mean;
      }
      area = areaInPx * nmPerPixel * nmPerPixel / 1e6;
      weightedArea = weightedAreaInPx * nmPerPixel * nmPerPixel / 1e6;
      noOfMolecules = molecules.size();

      // Pad the images to the largest scale being investigated by the correlation function.
      // Original Sengupta paper uses 800nm for the padding size.
      // Limit to within 80% of the minimum dimension of the image.
      double maxRadius = correlationDistance / nmPerPixel;
      final int imageSize = FastMath.min(im.getWidth(), im.getHeight());
      if (imageSize < 1.25 * maxRadius) {
        maxRadius = imageSize / 1.25;
      }
      final int pad = (int) Math.round(maxRadius);
      log("Analysing up to %.0f nm = %d pixels", maxRadius * nmPerPixel, pad);
      im = padImage(im, pad);
      w = padImage(w, pad);

      // // Used for debugging
      // {
      // ImageProcessor w2 = w.duplicate();
      // w2.setMinAndMax(0, 1);
      // PCPALMMolecules.displayImage(PCPALMMolecules.results.getName() + " Binary Image Mask", w2,
      // nmPerPixel);
      // }

      final double peakDensity = getDensity(im);

      // Create 2D auto-correlation
      double[][] gr;
      try {
        // Use the FFT library as it is multi-threaded. This may not be in the user's path.
        gr = computeAutoCorrelationCurveFFT(im, w, pad, nmPerPixel, peakDensity);
      } catch (final Exception ex) {
        // Default to the ImageJ built-in FHT
        gr = computeAutoCorrelationCurveFHT(im, w, pad, nmPerPixel, peakDensity);
      }
      if (gr == null) {
        return;
      }

      // Add the g(r) curve to the results
      addResult(peakDensity, gr);

      noPlots = WindowManager.getFrame(frequencyDomainTitle) == null;
      topPlotTitle = spatialPlotTitle;

      // Do not plot r=0 value on the curve
      plotCorrelation(gr, 1, frequencyDomainTitle, "g(r)", false, showErrorBars);
    }

    if (noPlots) {
      // Position the plot underneath the other one
      final Frame f1 = WindowManager.getFrame(topPlotTitle);
      if (f1 != null) {
        final String bottomPlotTitle =
            (topPlotTitle.equals(spatialPlotTitle) ? frequencyDomainTitle : spatialPlotTitle);
        final Frame f2 = WindowManager.getFrame(bottomPlotTitle);
        if (f2 != null) {
          f2.setLocation(f2.getLocation().x, f2.getLocation().y + f1.getHeight());
        }
      }
    }

    log("%s domain analysis computed in %s ms", (spatialDomain) ? "Spatial" : "Frequency",
        MathUtils.rounded((System.nanoTime() - start) * 1e-6, 4));
    log("---");
  }

  private static ImageProcessor applyWindow(ImageProcessor im, ImageWindow imageWindow) {
    float[] image = (float[]) im.toFloat(0, null).getPixels();
    image = imageWindow.applySeperable(image, im.getWidth(), im.getHeight(),
        ImageWindow.WindowMethod.TUKEY);
    return new FloatProcessor(im.getWidth(), im.getHeight(), image, null);
  }

  /**
   * Plot the correlation.
   *
   * @param gr the correlation curve
   * @param offset the offset
   * @param plotTitle the plot title
   * @param yAxisTitle the y axis title
   * @param barChart the bar chart
   * @param showErrorBars the show error bars
   * @return the plot
   */
  public static Plot2 plotCorrelation(double[][] gr, int offset, String plotTitle,
      String yAxisTitle, boolean barChart, boolean showErrorBars) {
    final double[] x = new double[gr[1].length - offset];
    final double[] y = new double[x.length];
    System.arraycopy(gr[0], offset, x, 0, x.length);
    System.arraycopy(gr[1], offset, y, 0, y.length);

    final Plot2 plot = new Plot2(plotTitle, "r (nm)", yAxisTitle);
    plot.setLimits(0, x[x.length - 1], MathUtils.min(y) * 0.95, MathUtils.max(y) * 1.05);
    plot.addPoints(x, y, (barChart) ? Plot2.BAR : Plot.LINE);

    ImageJUtils.display(plotTitle, plot);

    if (showErrorBars && !barChart) {
      plot.setColor(Color.magenta);
      for (int i = 0; i < x.length; i++) {
        final double sd = gr[2][i + offset];
        plot.drawLine(x[i], y[i] - sd, x[i], y[i] + sd);
      }
      ImageJUtils.display(plotTitle, plot);
    }

    return plot;
  }

  /**
   * Pad the image by the specified number of pixels.
   *
   * @param im the image
   * @param pad the pad
   * @return the padded image
   */
  private static ImageProcessor padImage(ImageProcessor im, int pad) {
    // int newW = pad * 2 + im.getWidth();
    // int newH = pad * 2 + im.getHeight();
    final int newW = pad + im.getWidth();
    final int newH = pad + im.getHeight();

    final ImageProcessor im2 = im.createProcessor(newW, newH);
    // im2.insert(im, pad, pad);
    im2.insert(im, 0, 0);
    return im2;
  }

  /**
   * Create a weight image of the same size. All pixels corresponding to the original image area are
   * set to 1. A window function is optionally applied.
   *
   * @param im the image
   * @param applyWindow the apply window flag
   * @return The weight image
   */
  private static ImageProcessor createWeightImage(ImageProcessor im, boolean applyWindow) {
    final float[] data = new float[im.getWidth() * im.getHeight()];
    Arrays.fill(data, 1);
    ImageProcessor w = new FloatProcessor(im.getWidth(), im.getHeight(), data, null);
    if (applyWindow) {
      w = applyWindow(w, imageWindow);
    }
    return w;
  }

  /**
   * Compute the auto-correlation curve using FHT (ImageJ built-in). Computes the correlation image
   * and then samples the image at radii up to the specified length to get the average correlation
   * at a given radius.
   *
   * @param im the image
   * @param w the w
   * @param maxRadius the max radius
   * @param nmPerPixel the nm per pixel
   * @param density the density
   * @return { distances[], gr[], gr_se[] }
   */
  private static double[][] computeAutoCorrelationCurveFHT(ImageProcessor im, ImageProcessor w,
      int maxRadius, double nmPerPixel, double density) {
    log("Creating Hartley transforms");
    final Fht fht2Im = padToFHT2(im);
    final Fht fht2W = padToFHT2(w);
    if (fht2Im == null || fht2W == null) {
      error("Unable to perform Hartley transform");
      return null;
    }

    log("Performing correlation");
    final FloatProcessor corrIm = computeAutoCorrelationFHT(fht2Im);
    final FloatProcessor corrW = computeAutoCorrelationFHT(fht2W);

    IJ.showProgress(1);

    final int centre = corrIm.getHeight() / 2;
    final Rectangle crop =
        new Rectangle(centre - maxRadius, centre - maxRadius, maxRadius * 2, maxRadius * 2);
    if (showCorrelationImages) {
      displayCorrelation(corrIm, "Image correlation", crop);
      displayCorrelation(corrW, "Window correlation", crop);
    }

    log("Normalising correlation");
    final FloatProcessor correlation = normaliseCorrelation(corrIm, corrW, density);

    if (showCorrelationImages) {
      displayCorrelation(correlation, "Normalised correlation", crop);
    }

    return computeRadialAverage(maxRadius, nmPerPixel, correlation);
  }

  /**
   * Gets the density of peaks in the image. The density is in squared pixels.
   *
   * @param im the image
   * @return The density (in pixels^-2)
   */
  private double getDensity(ImageProcessor im) {
    // PCPALMMolecules.densityPeaks is in nm^-2
    final double density = PCPALMMolecules.densityPeaks * 1e6;

    // Alternatively use the density in the sample
    final double sampleDensity = noOfMolecules / area;

    // Actually count the density in the image
    uniquePoints = 0;
    for (int i = im.getPixelCount(); i-- > 0;) {
      // The image may not be binary so use the number
      uniquePoints += im.getf(i);
    }
    final double imageDensity = uniquePoints / weightedArea;

    log("  %d molecules plotted as %.1f unique points", noOfMolecules, uniquePoints);
    log("  Total Density = %g um^-2, Sample density = %g (%.2fx), Image density = %g (%.2fx)",
        density, sampleDensity, sampleDensity / density, imageDensity, imageDensity / density);

    // This is the method used by the PC-PALM MATLAB code.
    // The sum of the image divided by the sum of the normalisation window function
    return uniquePoints / weightedAreaInPx;
  }

  /**
   * Compute the auto-correlation curve using FFT. Computes the correlation image and then samples
   * the image at radii up to the specified length to get the average correlation at a given radius.
   *
   * @param im the image
   * @param w the w
   * @param maxRadius the max radius
   * @param nmPerPixel the nm per pixel
   * @param density the density
   * @return { distances[], gr[], gr_se[] }
   */
  private static double[][] computeAutoCorrelationCurveFFT(ImageProcessor im, ImageProcessor w,
      int maxRadius, double nmPerPixel, double density) {
    log("Performing FFT correlation");
    final FloatProcessor corrIm = computeAutoCorrelationFFT(im);
    final FloatProcessor corrW = computeAutoCorrelationFFT(w);
    if (corrIm == null || corrW == null) {
      error("Unable to perform Fourier transform");
      return null;
    }

    final int centre = corrIm.getHeight() / 2;
    final Rectangle crop =
        new Rectangle(centre - maxRadius, centre - maxRadius, maxRadius * 2, maxRadius * 2);
    if (showCorrelationImages) {
      displayCorrelation(corrIm, "Image correlation FFT", crop);
      displayCorrelation(corrW, "Window correlation FFT", crop);
    }

    log("  Normalising correlation");
    final FloatProcessor correlation = normaliseCorrelation(corrIm, corrW, density);

    if (showCorrelationImages) {
      displayCorrelation(correlation, "Normalised correlation FFT", crop);
    }

    return computeRadialAverage(maxRadius, nmPerPixel, correlation);
  }

  /**
   * Compute the radial average correlation function (gr).
   *
   * @param maxRadius the maximum radius to process (in pixels)
   * @param nmPerPixel covert pixel distances to nm
   * @param correlation auto-correlation
   * @return { distances[], gr[], gr_se[] }
   */
  private static double[][] computeRadialAverage(int maxRadius, double nmPerPixel,
      FloatProcessor correlation) {
    // Perform averaging of the correlation function using integer distance bins
    log("  Computing distance vs correlation curve");
    final int centre = correlation.getHeight() / 2;

    // Count the number of pixels at each distance and sum the correlations
    final Statistics[] gr = new Statistics[maxRadius + 1];

    // Cache distances
    final int[] d2 = new int[maxRadius + 1];
    for (int dy = 0; dy <= maxRadius; dy++) {
      gr[dy] = new Statistics();
      d2[dy] = dy * dy;
    }
    final int[][] distance = new int[maxRadius + 1][maxRadius + 1];
    for (int dy = 0; dy <= maxRadius; dy++) {
      for (int dx = dy; dx <= maxRadius; dx++) {
        distance[dy][dx] = distance[dx][dy] = (int) Math.round(Math.sqrt(d2[dx] + d2[dy]));
      }
    }

    final float[] data = (float[]) correlation.getPixels();
    for (int dy = -maxRadius; dy <= maxRadius; dy++) {
      final int absY = Math.abs(dy);
      int index = (centre + dy) * correlation.getWidth() + centre - maxRadius;
      for (int dx = -maxRadius; dx <= maxRadius; dx++, index++) {
        final int d = distance[absY][Math.abs(dx)];
        if (d > maxRadius || d == 0) {
          continue;
        }
        gr[d].add(data[index]);
      }
    }

    // Create the final data: a curve showing distance (in nm) verses the average correlation
    final double[] x = new double[maxRadius + 1];
    final double[] y = new double[maxRadius + 1];
    final double[] sd = new double[maxRadius + 1];
    for (int i = 0; i < x.length; i++) {
      x[i] = i * nmPerPixel;
      y[i] = gr[i].getMean();
      sd[i] = gr[i].getStandardError();
    }
    y[0] = correlation.getf(centre, centre);

    // For debugging
    // double[] H = new double[x.length];
    // for (int i = 0; i < x.length; i++)
    // H[i] = gr[i].getN();
    // Plot2 p = new Plot2("Histogram", "r", "F", x, H);
    // Utils.display("Histogram", p);

    return new double[][] {x, y, sd};
  }

  private static void displayCorrelation(FloatProcessor correlation, String title, Rectangle crop) {
    correlation.setRoi(crop);
    final ImageProcessor ip = correlation.crop();
    ip.resetMinAndMax();
    ImageJUtils.display(title, ip);
  }

  /**
   * Compute auto correlation FHT.
   *
   * @param fftIm in frequency domain
   * @return the auto correlation FHT.
   */
  private static FloatProcessor computeAutoCorrelationFHT(Fht fftIm) {
    final Fht FHT2 = fftIm.conjugateMultiply(fftIm);
    FHT2.inverseTransform();
    FHT2.swapQuadrants();
    return FHT2;
  }

  private static FloatProcessor normaliseCorrelation(FloatProcessor corrIm, FloatProcessor corrW,
      double density) {
    final float[] data = new float[corrIm.getWidth() * corrIm.getHeight()];
    final float[] dataIm = (float[]) corrIm.getPixels();
    final float[] dataW = (float[]) corrW.getPixels();

    // Square the density for normalisation
    density *= density;

    for (int i = 0; i < data.length; i++) {
      data[i] = (float) (dataIm[i] / (density * dataW[i]));
    }
    final FloatProcessor correlation =
        new FloatProcessor(corrIm.getWidth(), corrIm.getHeight(), data, null);
    return correlation;
  }

  /**
   * Pads the image to the next power of two and transforms into the frequency domain.
   *
   * @param ip the image
   * @return An FHT2 image in the frequency domain
   */
  private static Fht padToFHT2(ImageProcessor ip) {
    final FloatProcessor im2 = pad(ip);
    if (im2 == null) {
      return null;
    }
    final Fht FHT2 = new Fht(im2);
    FHT2.transform();
    return FHT2;
  }

  /**
   * Pads the image to the next power of two.
   *
   * @param ip the image
   * @return padded image
   */
  private static FloatProcessor pad(ImageProcessor ip) {
    // Pad to a power of 2
    final int size = FastMath.max(ip.getWidth(), ip.getHeight());
    final int newSize = nextPowerOfTwo(size);
    if (size > newSize) {
      return null; // Error
    }

    final FloatProcessor im2 = new FloatProcessor(newSize, newSize);
    // If the binary processor has a min and max this breaks the conversion to float
    // since values are mapped from 0-255 using the min-max look-up calibration table
    ip.resetMinAndMax();
    im2.insert(ip, 0, 0);
    return im2;
  }

  private static int nextPowerOfTwo(final int size) {
    return MathUtils.nextPow2(size);

    // int newSize = 0;
    // for (int i = 4; i < 15; i++)
    // {
    // newSize = (int) Math.pow(2.0, i);
    // if (size <= newSize)
    // {
    // break;
    // }
    // }
    // return newSize;
  }

  /**
   * Compute the auto-correlation using the JTransforms FFT library.
   *
   * @param ip the image
   * @return the auto correlation.
   */
  private static FloatProcessor computeAutoCorrelationFFT(ImageProcessor ip) {
    final FloatProcessor paddedIp = pad(ip);
    if (paddedIp == null) {
      return null;
    }
    final int size = paddedIp.getWidth();

    final boolean doubleFFT = false;
    final float[] pixels = (float[]) paddedIp.getPixels();

    final float[] correlation = new float[size * size];

    // JTransform library
    // ------------------
    // The data is stored in 1D array in row-major order. Complex number is
    // stored as two float values in sequence: the real and imaginary part

    // Correlation = fft^-1(abs(fft).^2)
    // The absolute value of a complex number z = x + y*i is the value sqrt(x*x+y*y).

    if (doubleFFT) {
      final DoubleFFT_2D fft = new DoubleFFT_2D(size, size);
      final double[] data = new double[size * size * 2];
      for (int i = 0; i < pixels.length; i++) {
        data[i] = pixels[i];
      }
      fft.realForwardFull(data);

      // Re-use data
      for (int i = 0, j = 0; i < data.length; i += 2, j++) {
        data[j] = data[i] * data[i] + data[i + 1] * data[i + 1];
      }
      // Zero fill
      for (int j = correlation.length; j < data.length; j++) {
        data[j] = 0;
      }

      // Re-use the pre-computed object
      // fft = new DoubleFFT_2D(size, size);
      fft.realInverseFull(data, true);

      // Get the real part of the data
      for (int i = 0, j = 0; i < data.length; i += 2, j++) {
        correlation[j] = (float) data[i];
      }
    } else {
      final FloatFFT_2D fft = new FloatFFT_2D(size, size);
      final float[] data = new float[size * size * 2];
      System.arraycopy(pixels, 0, data, 0, pixels.length);
      fft.realForwardFull(data);

      // Re-use data
      for (int i = 0, j = 0; i < data.length; i += 2, j++) {
        data[j] = data[i] * data[i] + data[i + 1] * data[i + 1];
      }
      // Zero fill
      for (int j = correlation.length; j < data.length; j++) {
        data[j] = 0;
      }

      // Re-use the pre-computed object
      // fft = new FloatFFT_2D(size, size);
      fft.realInverseFull(data, true);

      // Get the real part of the data
      for (int i = 0, j = 0; i < data.length; i += 2, j++) {
        correlation[j] = data[i];
      }
    }

    // Swap quadrants
    final FloatProcessor fp = new FloatProcessor(size, size, correlation, null);
    Fht.swapQuadrants(fp);
    return fp;
  }

  private void addResult(double peakDensity, double[][] gr) {
    final int id = results.size() + 1;

    // Convert density from pixel^-2 to um^-2
    peakDensity *= 1e6 / (nmPerPixel * nmPerPixel);

    final CorrelationResult result = new CorrelationResult(id, PCPALMMolecules.results.getSource(),
        minx, miny, maxx, maxy, uniquePoints, nmPerPixel, peakDensity, binaryImage, gr, false);
    results.add(result);

    createResultsTable();

    final double pcw = (PCPALMMolecules.maxx - PCPALMMolecules.minx) / 100.0;
    final double pch = (PCPALMMolecules.maxy - PCPALMMolecules.miny) / 100.0;

    final StringBuilder sb = new StringBuilder();
    sb.append(id).append('\t');
    sb.append(PCPALMMolecules.results.getName()).append('\t');
    sb.append(IJ.d2s(minx)).append('\t');
    sb.append(IJ.d2s((minx) / pcw)).append('\t');
    sb.append(IJ.d2s(miny)).append('\t');
    sb.append(IJ.d2s((miny) / pch)).append('\t');
    sb.append(IJ.d2s(maxx - minx)).append('\t');
    sb.append(IJ.d2s((maxx - minx) / pcw)).append('\t');
    sb.append(IJ.d2s(maxy - miny)).append('\t');
    sb.append(IJ.d2s((maxy - miny) / pch)).append('\t');
    sb.append(MathUtils.rounded(uniquePoints, 4)).append('\t');
    sb.append(MathUtils.rounded(peakDensity, 4)).append('\t');
    sb.append(MathUtils.rounded(nmPerPixel, 4)).append('\t');
    sb.append(binaryImage).append('\t');
    resultsTable.append(sb.toString());
  }

  private static void createResultsTable() {
    if (resultsTable == null || !resultsTable.isVisible()) {
      final StringBuilder sb = new StringBuilder();
      sb.append("ID\t");
      sb.append("Image Source\t");
      sb.append("X\t");
      sb.append("X %\t");
      sb.append("Y\t");
      sb.append("Y %\t");
      sb.append("Width\t");
      sb.append("Width %\t");
      sb.append("Height\t");
      sb.append("Height %\t");
      sb.append("N\t");
      sb.append("PeakDensity (um^-2)\t");
      sb.append("nm/pixel\t");
      sb.append("Binary\t");
      resultsTable = new TextWindow(TITLE, sb.toString(), (String) null, 800, 300);
    }
  }
}
