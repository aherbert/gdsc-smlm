/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Package
 *
 * Software for single molecule localisation microscopy (SMLM) in ImageJ
 * %%
 * Copyright (C) 2011 - 2025 Alex Herbert
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
import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.text.TextWindow;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import java.awt.Color;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicReference;
import java.util.stream.Collectors;
import org.apache.commons.statistics.descriptive.Median;
import org.apache.commons.statistics.inference.KolmogorovSmirnovTest;
import org.apache.commons.statistics.inference.KolmogorovSmirnovTest.TwoResult;
import uk.ac.sussex.gdsc.core.data.utils.ConversionException;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.ConfigurationException;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.count.FrameCounter;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.PrecisionResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.sort.IdFramePeakResultComparator;

/**
 * Compare the jump distances of two datasets of traced data.
 */
public class CompareJumpDistances implements PlugIn {
  private static final String TITLE = "Compare Jump Distances";

  private static final AtomicReference<TextWindow> TABLE_REF = new AtomicReference<>();

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> INSTANCE = new AtomicReference<>(new Settings());

    String inputOption1;
    String inputOption2;
    int frames;
    boolean precisionCorrection;

    Settings() {
      // Set defaults
      inputOption1 = "";
      inputOption2 = "";
      frames = 1;
    }

    Settings(Settings source) {
      inputOption1 = source.inputOption1;
      inputOption2 = source.inputOption2;
      frames = source.frames;
      precisionCorrection = source.precisionCorrection;
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
      return INSTANCE.get().copy();
    }

    /**
     * Save the settings.
     */
    void save() {
      INSTANCE.set(this);
    }
  }

  /**
   * Store the localisation position.
   */
  private static class Position {
    int t;
    float x;
    float y;

    Position(int t, float x, float y) {
      this.t = t;
      this.x = x;
      this.y = y;
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
    MemoryPeakResults results1 =
        ResultsManager.loadInputResults(settings.inputOption1, false, null, null);
    MemoryPeakResults results2 =
        ResultsManager.loadInputResults(settings.inputOption2, false, null, null);
    if (MemoryPeakResults.isEmpty(results1) || MemoryPeakResults.isEmpty(results2)) {
      IJ.error(TITLE, "No results could be loaded");
      return;
    }

    // Results must have the same calibration
    CalibrationReader cal1 = results1.getCalibrationReader();
    CalibrationReader cal2 = results2.getCalibrationReader();
    if (cal1.getNmPerPixel() != cal2.getNmPerPixel()) {
      IJ.error(TITLE, String.format("Distance calibration mismatch: %.3f != %.3f nm/px",
          cal1.getNmPerPixel(), cal2.getNmPerPixel()));
      return;
    }
    if (cal1.getExposureTime() != cal2.getExposureTime()) {
      IJ.error(TITLE, String.format("Exposure time mismatch: %.3f != %.3f ms/frame",
          cal1.getExposureTime(), cal2.getExposureTime()));
      return;
    }

    // Extract the jump distances
    double[] distances1 = getDistances(results1, settings.frames);
    double[] distances2 = getDistances(results2, settings.frames);

    if ((distances1.length & distances2.length) == 0) {
      IJ.error(TITLE, "No distances for time delay: " + settings.frames);
      return;
    }

    TypeConverter<DistanceUnit> distanceConverter;
    try {
      distanceConverter = cal1.getDistanceConverter(DistanceUnit.UM);
    } catch (final ConversionException | ConfigurationException ex) {
      IJ.error(TITLE, "Cannot convert units to um or seconds: " + ex.getMessage());
      return;
    }

    Arrays.sort(distances1);
    Arrays.sort(distances2);

    // Apply precision correction
    double error1 = 0;
    double error2 = 0;
    if (settings.precisionCorrection) {
      // Get the localisation error (4s^2) in raw units^2
      error1 = getLocalisationError(results1, distanceConverter);
      error2 = getLocalisationError(results2, distanceConverter);
      applyCorrection(distances1, error1);
      applyCorrection(distances2, error2);
    }

    // KS test
    TwoResult r = KolmogorovSmirnovTest.withDefaults().test(distances1, distances2);

    // Plot cumulative histogram
    double[][] h1 = MathUtils.cumulativeHistogram(distances1, true);
    double[][] h2 = MathUtils.cumulativeHistogram(distances2, true);
    Plot plot = new Plot(TITLE,
        String.format("Distance (um/%s)", TextUtils.pleural(settings.frames, "frame")),
        "Probability");
    SimpleArrayUtils.apply(h1[0], distanceConverter::convert);
    SimpleArrayUtils.apply(h2[0], distanceConverter::convert);
    plot.setColor(Color.RED);
    plot.addPoints(h1[0], h1[1], Plot.LINE);
    plot.setColor(Color.BLUE);
    plot.addPoints(h2[0], h2[1], Plot.LINE);
    plot.setColor(Color.BLACK);
    plot.addLabel(0, 0, String.format("KS Test: %.4g", r.getPValue()));
    ImageJUtils.display(TITLE, plot);

    // QQ plot

    // Table of results
    // convert (4s^2) to s in nm
    final double scale = distanceConverter.convert(1) * 0.5e3;
    addResult(settings, scale * Math.sqrt(error1), scale * Math.sqrt(error2),
        distances1.length, distances2.length, r);
  }

  private boolean showDialog() {
    settings = Settings.load();
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addMessage("Compare the jump distances of two traced datasets");
    ResultsManager.addInput(gd, "Input_1", settings.inputOption1, InputSource.MEMORY_CLUSTERED);
    ResultsManager.addInput(gd, "Input_2", settings.inputOption2, InputSource.MEMORY_CLUSTERED);
    gd.addSlider("Frames", 1, 10, settings.frames);
    gd.addCheckbox("Precision_correction", settings.precisionCorrection);
    gd.addHelp(HelpUrls.getUrl("compare-jump-distances"));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }
    settings.inputOption1 = ResultsManager.getInputSource(gd);
    settings.inputOption2 = ResultsManager.getInputSource(gd);
    settings.frames = (int) gd.getNextNumber();
    settings.precisionCorrection = gd.getNextBoolean();
    settings.save();
    return true;
  }

  private static double getLocalisationError(MemoryPeakResults results,
      TypeConverter<DistanceUnit> distanceConverter) {
    try {
      final PrecisionResultProcedure p = new PrecisionResultProcedure(results);
      p.getPrecision();

      // Precision in nm using the median
      double precision = Median.withDefaults().evaluate(p.precisions);
      // Convert from nm to um to raw units
      final double rawPrecision = distanceConverter.convertBack(precision * 1e-3);
      // Get the localisation error (4s^2) in units^2
      return 4 * rawPrecision * rawPrecision;
    } catch (final ConversionException | ConfigurationException ex) {
      ImageJUtils.log(TITLE + " - Unable to compute precision: " + ex.getMessage());
    }
    return 0;
  }

  /**
   * Get the jump distances for the specified time delay. Adapted from
   * {@link TrackDiffusionAnalysis} for a single time delay.
   *
   * @param results the results
   * @param t the time delay
   * @return distances
   */
  private static double[] getDistances(MemoryPeakResults results, int t) {
    final DoubleArrayList distances = new DoubleArrayList();
    final LocalList<Position> track = new LocalList<>();

    results.sort(IdFramePeakResultComparator.INSTANCE);
    final FrameCounter id = new FrameCounter(-1);
    results.forEach((PeakResultProcedure) r -> {
      if (id.advance(r.getId())) {
        if (!track.isEmpty()) {
          // Process track
          final int maxStart = track.size() - t;
          for (int i = 0; i < maxStart; i++) {
            final Position origin = track.unsafeGet(i);
            final Position position = track.unsafeGet(i + t);
            if (position.t - origin.t == t) {
              distances.add(MathUtils.distance(origin.x, origin.y, position.x, position.y));
            }
          }
          track.clear();
        }
      }
      track.add(new Position(r.getFrame(), r.getXPosition(), r.getYPosition()));
    });
    // Process final track
    if (!track.isEmpty()) {
      // Process track
      final int maxStart = track.size() - t;
      for (int i = 0; i < maxStart; i++) {
        final Position origin = track.unsafeGet(i);
        final Position position = track.unsafeGet(i + t);
        if (position.t - origin.t == t) {
          distances.add(MathUtils.distance(origin.x, origin.y, position.x, position.y));
        }
      }
    }
    return distances.toDoubleArray();
  }

  private static void applyCorrection(double[] distances, double error) {
    if (error == 0) {
      return;
    }
    int i = 0;
    while (i < distances.length && distances[i] < error) {
      distances[i] = 0;
      i++;
    }
    while (i < distances.length) {
      distances[i] -= error;
      i++;
    }
  }

  private TextWindow createTable() {
    return ImageJUtils.refresh(TABLE_REF, () -> {
      return new TextWindow(TITLE + " Results", createHeader(), "", 1500, 300);
    });
  }

  private String createHeader() {
    return Arrays.stream(new String[] {"Input1", "Input2", "Precision1 (nm)", "Precision2 (nm)",
        "Frames", "N1", "N2", "KS D", "p(D)"}).collect(Collectors.joining("\t"));
  }

  private void addResult(Settings settings, double precision1, double precision2,
      int n1, int n2, TwoResult r) {
    final StringBuilder sb = new StringBuilder(1024);
    //@formatter:off
    sb.append(settings.inputOption1).append('\t')
      .append(settings.inputOption2).append('\t')
      .append(MathUtils.rounded(precision1)).append('\t')
      .append(MathUtils.rounded(precision2)).append('\t')
      .append(settings.frames).append('\t')
      .append(n1).append('\t')
      .append(n2).append('\t')
      .append(r.getStatistic()).append('\t')
      .append(r.getPValue());
    //@formatter:on
    createTable().append(sb.toString());
  }
}
