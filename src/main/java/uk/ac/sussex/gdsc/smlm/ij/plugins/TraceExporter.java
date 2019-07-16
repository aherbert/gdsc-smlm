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
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;
import uk.ac.sussex.gdsc.core.utils.rng.SplitMix;
import uk.ac.sussex.gdsc.smlm.data.NamedObject;
import uk.ac.sussex.gdsc.smlm.data.config.UnitConverterUtils;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.TimeUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.MultiDialog.MemoryResultsItems;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.AttributePeakResult;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.XyrResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.sort.IdFramePeakResultComparator;

import gnu.trove.set.hash.TIntHashSet;

import ij.IJ;
import ij.plugin.PlugIn;

import org.apache.commons.rng.sampling.distribution.NormalizedGaussianSampler;
import org.apache.commons.rng.simple.internal.SeedFactory;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Plugin to export traced datasets.
 */
public class TraceExporter implements PlugIn {
  private static final String TITLE = "Trace Exporter";
  private static List<String> selected;
  private static String directory = "";
  private static int minLength = 2;
  private static int maxJump = 1;
  private static double wobble;

  private static String[] formatNames;
  private static int format;

  private enum ExportFormat implements NamedObject {
    SPOT_ON("Spot-On");

    private static final ExportFormat[] values = values();

    private final String name;

    ExportFormat(String name) {
      this.name = name;
    }

    @Override
    public String getName() {
      return name;
    }

    @Override
    public String getShortName() {
      return name;
    }

    /**
     * Get the enum value from the index ordinal.
     *
     * @param index the index
     * @return the export format
     */
    @SuppressWarnings("unused")
    public static ExportFormat fromOrdinal(int index) {
      return values[MathUtils.clip(0, values.length - 1, index)];
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    final MemoryResultsItems items = new MemoryResultsItems(MemoryPeakResults::hasId);

    if (items.size() == 0) {
      IJ.error(TITLE, "No traced localisations in memory");
      return;
    }

    // Get input options
    if (!showDialog()) {
      return;
    }

    final ArrayList<MemoryPeakResults> allResults = new ArrayList<>();

    // Pick multiple input datasets together using a list box.
    if (!showMultiDialog(allResults, items)) {
      return;
    }

    for (final MemoryPeakResults results : allResults) {
      export(results);
    }
  }

  private static boolean showDialog() {
    if (formatNames == null) {
      formatNames = SettingsManager.getNames((Object[]) ExportFormat.values());
    }
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addMessage("Export traces to a directory");
    gd.addDirectoryField("Directory", directory, 30);
    gd.addNumericField("Min_length", minLength, 0);
    gd.addMessage("Specify the maximum jump allowed within a trace.\n"
        + "Traces with larger jumps will be split.");
    gd.addNumericField("Max_jump", maxJump, 0);
    gd.addMessage("Specify localistion precision (wobble) to add");
    gd.addNumericField("Wobble", wobble, 0, 6, "nm");
    gd.addChoice("Format", formatNames, formatNames[format]);
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }
    directory = gd.getNextString();
    minLength = (int) Math.abs(gd.getNextNumber());
    maxJump = (int) Math.abs(gd.getNextNumber());
    wobble = Math.abs(gd.getNextNumber());
    format = gd.getNextChoiceIndex();
    return true;
  }

  private static boolean showMultiDialog(ArrayList<MemoryPeakResults> allResults,
      MemoryResultsItems items) {
    // Show a list box containing all the results. This should remember the last set of chosen
    // items.
    final MultiDialog md = new MultiDialog(TITLE, items);
    md.addSelected(selected);

    md.showDialog();

    if (md.wasCancelled()) {
      return false;
    }

    selected = md.getSelectedResults();
    if (selected.isEmpty()) {
      IJ.error(TITLE, "No results were selected");
      return false;
    }

    for (final String name : selected) {
      final MemoryPeakResults r = MemoryPeakResults.getResults(name);
      if (r != null) {
        allResults.add(r);
      }
    }

    return !allResults.isEmpty();
  }

  private static void export(MemoryPeakResults results) {
    // Copy to allow manipulation
    results = results.copy();

    // Strip results with no trace Id
    results.removeIf(result -> result.getId() <= 0);

    // Sort by ID then time
    results.sort(IdFramePeakResultComparator.INSTANCE);

    // Split traces with big jumps
    results = splitTraces(results);

    // Count each ID and remove short traces
    int id = 0;
    int count = 0;
    final TIntHashSet remove = new TIntHashSet();
    for (int i = 0, size = results.size(); i < size; i++) {
      if (results.get(i).getId() != id) {
        if (count < minLength) {
          remove.add(id);
        }
        count = 0;
        id = results.get(i).getId();
      }
      count++;
    }
    // Final ID
    if (count < minLength) {
      remove.add(id);
    }

    if (!remove.isEmpty()) {
      results.removeIf(result -> remove.contains(result.getId()));
      results.sort(IdFramePeakResultComparator.INSTANCE);
    }

    if (wobble > 0) {
      // Just leave any exceptions to trickle up and kill the plugin
      final TypeConverter<DistanceUnit> c = results.getDistanceConverter(DistanceUnit.NM);
      final double w = c.convertBack(wobble);
      final SplitMix sm = new SplitMix(SeedFactory.createLong());
      final NormalizedGaussianSampler gauss = SamplerUtils.createNormalizedGaussianSampler(sm);
      final boolean is3D = results.is3D();
      results.forEach((PeakResultProcedure) peakResult -> {
        peakResult.setXPosition((float) (peakResult.getXPosition() + w * gauss.sample()));
        peakResult.setYPosition((float) (peakResult.getYPosition() + w * gauss.sample()));
        if (is3D) {
          peakResult.setZPosition((float) (peakResult.getZPosition() + w * gauss.sample()));
        }
      });
    }

    // Only one format at the moment ...
    exportSpotOn(results);
  }

  private static MemoryPeakResults splitTraces(MemoryPeakResults results) {
    if (maxJump < 1) {
      // Disabled
      return results;
    }

    int id = 0;
    final int lastT = 0;
    for (int i = 0, size = results.size(); i < size; i++) {
      final PeakResult r = results.get(i);
      if (r.getId() != id) {
        id = r.getId();
      } else if (r.getFrame() - lastT > maxJump) {
        return doSplit(results);
      }
    }
    return results;
  }

  private static MemoryPeakResults doSplit(MemoryPeakResults results) {
    final MemoryPeakResults results2 = new MemoryPeakResults(results.size());
    results2.copySettings(results);
    int nextId = results.getLast().getId();
    int id = 0;
    int idOut = 0;
    int lastT = 0;
    for (int i = 0, size = results.size(); i < size; i++) {
      final PeakResult r = results.get(i);
      if (r.getId() != id) {
        id = r.getId();
        idOut = id;
      } else if (r.getFrame() - lastT > maxJump) {
        idOut = ++nextId;
      }
      final AttributePeakResult r2 = new AttributePeakResult(r);
      r2.setId(idOut);
      results2.add(r2);
      lastT = r.getEndFrame();
    }
    return results2;
  }

  private static void exportSpotOn(MemoryPeakResults results) {
    // Simple Spot-On CSV file format:
    // https://spoton.berkeley.edu/SPTGUI/docs/latest#input-formats
    // frame, t (seconds), trajectory (trace id), x (um), y (um)

    try (BufferedWriter out =
        Files.newBufferedWriter(Paths.get(directory, results.getName() + ".csv"))) {
      out.write("frame,t,trajectory,x,y");
      out.newLine();

      final TypeConverter<TimeUnit> converter = UnitConverterUtils.createConverter(TimeUnit.FRAME,
          TimeUnit.SECOND, results.getCalibrationReader().getExposureTime());

      @SuppressWarnings("resource")
      final BufferedWriter writer = out;
      results.forEach(DistanceUnit.UM, (XyrResultProcedure) (x, y, result) -> {
        try {
          if (result.hasEndFrame()) {
            final String sId = Integer.toString(result.getId());
            final String sx = Float.toString(x);
            final String sy = Float.toString(y);
            for (int t = result.getFrame(); t <= result.getEndFrame(); t++) {
              writer.write(Integer.toString(t));
              writer.write(",");
              writer.write(Float.toString(converter.convert(t)));
              writer.write(",");
              writer.write(sId);
              writer.write(",");
              writer.write(sx);
              writer.write(",");
              writer.write(sy);
              writer.newLine();
            }
          } else {
            writer.write(Integer.toString(result.getFrame()));
            writer.write(",");
            writer.write(Float.toString(converter.convert(result.getFrame())));
            writer.write(",");
            writer.write(Integer.toString(result.getId()));
            writer.write(",");
            writer.write(Float.toString(x));
            writer.write(",");
            writer.write(Float.toString(y));
            writer.newLine();
          }
        } catch (final IOException ex) {
          handleException(ex);
        }
      });
    } catch (final IOException ex) {
      handleException(ex);
    }
  }

  private static void handleException(Exception ex) {
    Logger.getLogger(TraceExporter.class.getName()).log(Level.SEVERE, "Failed to export", ex);
  }
}
