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
import java.util.Collections;
import java.util.List;
import java.util.concurrent.atomic.AtomicReference;
import java.util.stream.Collectors;
import org.apache.commons.statistics.descriptive.Median;
import org.apache.commons.statistics.inference.KolmogorovSmirnovTest;
import org.apache.commons.statistics.inference.KolmogorovSmirnovTest.TwoResult;
import uk.ac.sussex.gdsc.core.data.utils.ConversionException;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.MultiDialog;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.ConfigurationException;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.MemoryResultsList;
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

    List<String> selected;
    int frames;
    boolean precisionCorrection;

    Settings() {
      // Set defaults
      selected = Collections.emptyList();
      frames = 1;
    }

    Settings(Settings source) {
      selected = source.selected;
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

    final MemoryResultsList items = new MemoryResultsList(MemoryPeakResults::hasId);

    if (items.isEmpty()) {
      IJ.error(TITLE, "No traced localisations in memory");
      return;
    }

    final List<MemoryPeakResults> results = new LocalList<>();

    if (!showDialog() || !showMultiDialog(results, items)) {
      return;
    }

    // Extract the jump distances
    final double[][] distances =
        results.stream().map(r -> getDistances(r, settings.frames)).toArray(double[][]::new);

    for (int i = 0; i < distances.length; i++) {
      if (distances[i].length == 0) {
        IJ.error(TITLE,
            results.get(i).getName() + ": No distances for time delay " + settings.frames);
        return;
      }
    }

    TypeConverter<DistanceUnit> distanceConverter;
    try {
      distanceConverter = results.get(0).getDistanceConverter(DistanceUnit.UM);
    } catch (final ConversionException | ConfigurationException ex) {
      IJ.error(TITLE, "Cannot convert units to um: " + ex.getMessage());
      return;
    }

    for (int i = 0; i < distances.length; i++) {
      Arrays.sort(distances[i]);
    }

    // Apply precision correction
    final double[] error = new double[distances.length];
    if (settings.precisionCorrection) {
      // Get the localisation error (4s^2) in raw units^2
      for (int i = 0; i < distances.length; i++) {
        error[i] = getLocalisationError(results.get(i), distanceConverter);
        applyCorrection(distances[i], error[i]);
      }
    }

    // Table of results convert (4s^2) to s in nm
    final double scale = distanceConverter.convert(1) * 0.5e3;

    // All-vs-all KS test
    for (int i = 0; i < distances.length; i++) {
      for (int j = i + 1; j < distances.length; j++) {
        final TwoResult r = KolmogorovSmirnovTest.withDefaults().test(distances[i], distances[j]);
        addResult(settings, results.get(i).getName(), results.get(j).getName(),
            scale * Math.sqrt(error[i]), scale * Math.sqrt(error[j]), distances[i].length,
            distances[j].length, r);
      }
    }

    if (results.size() > 2) {
      return;
    }

    // Plot cumulative histogram
    final double[][] h1 = MathUtils.cumulativeHistogram(distances[0], true);
    final double[][] h2 = MathUtils.cumulativeHistogram(distances[1], true);
    final Plot plot = new Plot(TITLE,
        String.format("Distance (um/%s)", TextUtils.pleural(settings.frames, "frame")),
        "Probability");
    SimpleArrayUtils.apply(h1[0], distanceConverter::convert);
    SimpleArrayUtils.apply(h2[0], distanceConverter::convert);
    plot.setColor(Color.RED);
    plot.addPoints(h1[0], h1[1], Plot.LINE);
    plot.setColor(Color.BLUE);
    plot.addPoints(h2[0], h2[1], Plot.LINE);
    plot.setColor(Color.BLACK);
    ImageJUtils.display(TITLE, plot);

    // QQ plot

  }

  private boolean showDialog() {
    settings = Settings.load();

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addMessage("Compare the jump distances of traced datasets");
    gd.addSlider("Frames", 1, 10, settings.frames);
    gd.addCheckbox("Precision_correction", settings.precisionCorrection);
    gd.addHelp(HelpUrls.getUrl("compare-jump-distances"));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }
    settings.frames = (int) gd.getNextNumber();
    settings.precisionCorrection = gd.getNextBoolean();
    settings.save();
    return true;
  }

  private boolean showMultiDialog(List<MemoryPeakResults> allResults, MemoryResultsList items) {
    // Show a list box containing all the results. This should remember the last set of chosen
    // items.
    final MultiDialog md = new MultiDialog(TITLE, items);
    md.setDisplayConverter(items.getDisplayConverter());
    md.setSelected(settings.selected);
    md.setHelpUrl(HelpUrls.getUrl("compare-jump-distances"));

    md.showDialog();

    if (md.wasCancelled()) {
      return false;
    }

    final List<String> selected = md.getSelectedResults();
    if (selected.isEmpty()) {
      IJ.error(TITLE, "No results were selected");
      return false;
    }
    settings.selected = selected;

    for (final String name : selected) {
      final MemoryPeakResults r = MemoryPeakResults.getResults(name);
      if (r != null) {
        allResults.add(r);
      }
    }

    // Check calibration exists for the first set of results
    if (allResults.isEmpty()) {
      return false;
    }

    // Check the calibration is the same for the rest
    final CalibrationReader cal = allResults.get(0).getCalibrationReader();
    if (cal == null) {
      IJ.error(TITLE, "Uncalibrated results were selected");
      return false;
    }
    final double nmPerPixel = cal.getNmPerPixel();
    final double exposureTime = cal.getExposureTime();
    final DistanceUnit distanceUnit = cal.getDistanceUnit();
    for (int i = 1; i < allResults.size(); i++) {
      final MemoryPeakResults results = allResults.get(i);

      if (!results.hasCalibration()
          || results.getCalibrationReader().getExposureTime() != exposureTime
          || results.getNmPerPixel() != nmPerPixel || results.getDistanceUnit() != distanceUnit) {
        IJ.error(TITLE,
            "The exposure time, pixel pitch and distance unit must match across all the results");
        return false;
      }
    }

    return true;
  }

  private static double getLocalisationError(MemoryPeakResults results,
      TypeConverter<DistanceUnit> distanceConverter) {
    try {
      final PrecisionResultProcedure p = new PrecisionResultProcedure(results);
      p.getPrecision();

      // Precision in nm using the median
      final double precision = Median.withDefaults().evaluate(p.precisions);
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

    // Note that during processing we cannot use use j = i+t as the
    // track may have frame gaps

    results.sort(IdFramePeakResultComparator.INSTANCE);
    final FrameCounter id = new FrameCounter(-1);
    results.forEach((PeakResultProcedure) r -> {
      if (id.advance(r.getId())) {
        if (!track.isEmpty()) {
          // Process track
          for (int i = 1; i < track.size(); i++) {
            final Position origin = track.unsafeGet(i - 1);
            for (int j = i; j < track.size(); j++) {
              final Position position = track.unsafeGet(j);
              final int gap = position.t - origin.t;
              if (gap >= t) {
                if (gap == t) {
                  distances.add(MathUtils.distance(origin.x, origin.y, position.x, position.y));
                }
                break;
              }
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
      for (int i = 1; i < track.size(); i++) {
        final Position origin = track.unsafeGet(i - 1);
        for (int j = i; j < track.size(); j++) {
          final Position position = track.unsafeGet(j);
          final int gap = position.t - origin.t;
          if (gap >= t) {
            if (gap == t) {
              distances.add(MathUtils.distance(origin.x, origin.y, position.x, position.y));
            }
            break;
          }
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

  private void addResult(Settings settings, String input1, String input2, double precision1,
      double precision2, int n1, int n2, TwoResult r) {
    final StringBuilder sb = new StringBuilder(1024);
    //@formatter:off
    sb.append(input1).append('\t')
      .append(input2).append('\t')
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
