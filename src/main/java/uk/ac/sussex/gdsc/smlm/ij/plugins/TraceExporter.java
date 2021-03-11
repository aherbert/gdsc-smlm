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

import gnu.trove.set.hash.TIntHashSet;
import ij.IJ;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicReference;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.NormalizedGaussianSampler;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.MultiDialog;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;
import uk.ac.sussex.gdsc.core.utils.rng.UniformRandomProviders;
import uk.ac.sussex.gdsc.smlm.data.NamedObject;
import uk.ac.sussex.gdsc.smlm.data.config.UnitConverterUtils;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.TimeUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.MemoryResultsList;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.AttributePeakResult;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.count.Counter;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.XyrResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.sort.IdFramePeakResultComparator;
import us.hebi.matlab.mat.format.Mat5;
import us.hebi.matlab.mat.types.MatFile;
import us.hebi.matlab.mat.types.Matrix;

/**
 * Plugin to export traced datasets.
 */
public class TraceExporter implements PlugIn {
  private static final String TITLE = "Trace Exporter";
  private static AtomicReference<List<String>> selectedRef = new AtomicReference<>();

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());
    static final String[] formatNames;

    static {
      formatNames = SettingsManager.getNames((Object[]) ExportFormat.values());
      // We do not want capitalisation of 'ana'
      formatNames[ExportFormat.ANA_DNA.ordinal()] = ExportFormat.ANA_DNA.getName();
    }

    String directory;
    int minLength;
    int maxLength;
    int maxJump;
    double wobble;
    int format;
    boolean showTraceLengths;

    Settings() {
      // Set defaults
      directory = "";
      minLength = 2;
      maxJump = 1;
    }

    Settings(Settings source) {
      directory = source.directory;
      minLength = source.minLength;
      maxLength = source.maxLength;
      maxJump = source.maxJump;
      wobble = source.wobble;
      format = source.format;
      showTraceLengths = source.showTraceLengths;
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

  private enum ExportFormat implements NamedObject {
    SPOT_ON("Spot-On"), ANA_DNA("anaDDA");

    private final String name;

    ExportFormat(String name) {
      this.name = name;
    }

    @Override
    public String getName() {
      return name;
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

  private boolean showDialog() {
    settings = Settings.load();
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addMessage("Export traces to a directory");
    gd.addDirectoryField("Directory", settings.directory, 30);
    gd.addNumericField("Min_length", settings.minLength, 0);
    gd.addNumericField("Max_length", settings.maxLength, 0);
    gd.addMessage("Specify the maximum jump allowed within a trace.\n"
        + "Traces with larger jumps will be split.");
    gd.addNumericField("Max_jump", settings.maxJump, 0);
    gd.addMessage("Specify localistion precision (wobble) to add");
    gd.addNumericField("Wobble", settings.wobble, 0, 6, "nm");
    gd.addChoice("Format", Settings.formatNames, settings.format);
    gd.addCheckbox("Histogram_trace_lengths", settings.showTraceLengths);
    gd.addHelp(HelpUrls.getUrl("trace-exporter"));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }
    settings.directory = gd.getNextString();
    settings.minLength = Math.max(0, (int) gd.getNextNumber());
    settings.maxLength = Math.max(0, (int) gd.getNextNumber());
    settings.maxJump = Math.max(0, (int) gd.getNextNumber());
    settings.wobble = Math.max(0, gd.getNextNumber());
    settings.format = gd.getNextChoiceIndex();
    settings.showTraceLengths = gd.getNextBoolean();
    settings.save();
    return true;
  }

  private static boolean showMultiDialog(ArrayList<MemoryPeakResults> allResults,
      MemoryResultsList items) {
    // Show a list box containing all the results. This should remember the last set of chosen
    // items.
    final MultiDialog md = new MultiDialog(TITLE, items);
    md.setDisplayConverter(items.getDisplayConverter());
    md.setSelected(selectedRef.get());

    md.showDialog();

    if (md.wasCancelled()) {
      return false;
    }

    final List<String> selected = md.getSelectedResults();
    if (selected.isEmpty()) {
      IJ.error(TITLE, "No results were selected");
      return false;
    }
    selectedRef.set(selected);

    for (final String name : selected) {
      final MemoryPeakResults r = MemoryPeakResults.getResults(name);
      if (r != null) {
        allResults.add(r);
      }
    }

    return !allResults.isEmpty();
  }

  private void export(MemoryPeakResults results) {
    // Copy to allow manipulation
    results = results.copy();

    // Strip results with no trace Id
    results.removeIf(result -> result.getId() <= 0);

    // Sort by ID then time
    results.sort(IdFramePeakResultComparator.INSTANCE);

    // Split traces with big jumps and long tracks into smaller tracks
    results = splitTraces(results);
    results = splitLongTraces(results);

    // Count each ID and remove short traces
    int id = 0;
    int count = 0;
    int tracks = 0;
    int maxLength = 0;
    final TIntHashSet remove = new TIntHashSet();
    for (int i = 0, size = results.size(); i < size; i++) {
      final PeakResult result = results.get(i);
      if (result.getId() != id) {
        if (count < settings.minLength) {
          remove.add(id);
        } else {
          tracks++;
          maxLength = Math.max(maxLength, count);
        }
        count = 0;
        id = result.getId();
      }
      count += getLength(result);
    }
    // Final ID
    if (count < settings.minLength) {
      remove.add(id);
    } else {
      tracks++;
      maxLength = Math.max(maxLength, count);
    }

    if (!remove.isEmpty()) {
      results.removeIf(result -> remove.contains(result.getId()));
      results.sort(IdFramePeakResultComparator.INSTANCE);
    }

    if (settings.wobble > 0) {
      // Just leave any exceptions to trickle up and kill the plugin
      final TypeConverter<DistanceUnit> c = results.getDistanceConverter(DistanceUnit.NM);
      final double w = c.convertBack(settings.wobble);
      final UniformRandomProvider rng = UniformRandomProviders.create();
      final NormalizedGaussianSampler gauss = SamplerUtils.createNormalizedGaussianSampler(rng);
      final boolean is3D = results.is3D();
      results.forEach((PeakResultProcedure) peakResult -> {
        peakResult.setXPosition((float) (peakResult.getXPosition() + w * gauss.sample()));
        peakResult.setYPosition((float) (peakResult.getYPosition() + w * gauss.sample()));
        if (is3D) {
          peakResult.setZPosition((float) (peakResult.getZPosition() + w * gauss.sample()));
        }
      });
    }

    if (settings.format == 1) {
      exportAnaDda(results);
    } else {
      exportSpotOn(results);
    }
    ImageJUtils.log("Exported %s: %s in %s", results.getName(),
        TextUtils.pleural(results.size(), "localisation"), TextUtils.pleural(tracks, "track"));

    if (settings.showTraceLengths) {
      // We store and index (count-1)
      final int[] h = new int[maxLength];
      id = 0;
      for (int i = 0, size = results.size(); i < size; i++) {
        final PeakResult result = results.get(i);
        if (result.getId() != id) {
          h[count - 1]++;
          count = 0;
          id = result.getId();
        }
        count += getLength(result);
      }
      h[count - 1]++;
      final String title = TITLE + ": " + results.getName();
      final Plot plot = new Plot(title, "Length", "Frequency");
      plot.addPoints(SimpleArrayUtils.newArray(h.length, 1, 1.0f), SimpleArrayUtils.toFloat(h),
          Plot.BAR);
      plot.setLimits(SimpleArrayUtils.findIndex(h, i -> i != 0), maxLength + 1, 0, Double.NaN);
      ImageJUtils.display(title, plot);
    }
  }

  private static int getLength(PeakResult result) {
    return (result.hasEndFrame()) ? (result.getEndFrame() - result.getFrame() + 1) : 1;
  }

  private MemoryPeakResults splitTraces(MemoryPeakResults results) {
    if (settings.maxJump < 1) {
      // Disabled
      return results;
    }

    int id = 0;
    final int lastT = 0;
    for (int i = 0, size = results.size(); i < size; i++) {
      final PeakResult r = results.get(i);
      if (r.getId() != id) {
        id = r.getId();
      } else if (r.getFrame() - lastT > settings.maxJump) {
        return doSplit(results);
      }
    }
    return results;
  }

  private MemoryPeakResults doSplit(MemoryPeakResults results) {
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
      } else if (r.getFrame() - lastT > settings.maxJump) {
        idOut = ++nextId;
      }
      final AttributePeakResult r2 = new AttributePeakResult(r);
      r2.setId(idOut);
      results2.add(r2);
      lastT = r.getEndFrame();
    }
    return results2;
  }

  private MemoryPeakResults splitLongTraces(MemoryPeakResults results) {
    if (settings.maxLength < 1) {
      // Disabled
      return results;
    }

    final MemoryPeakResults results2 = new MemoryPeakResults(results.size());
    results2.copySettings(results);
    int nextId = results.getLast().getId();
    int id = 0;
    int idOut = 0;
    // The length that has been output under the current ID
    int length = 0;
    for (int i = 0, size = results.size(); i < size; i++) {
      final PeakResult r = results.get(i);
      if (r.getId() != id) {
        id = r.getId();
        idOut = id;
        length = 0;
      } else if (length >= settings.maxLength) {
        idOut = ++nextId;
        length = 0;
      }
      final AttributePeakResult r2 = new AttributePeakResult(r);
      r2.setId(idOut);
      results2.add(r2);
      length += getLength(r);
    }
    return results2;
  }

  private void exportSpotOn(MemoryPeakResults results) {
    // Simple Spot-On CSV file format:
    // https://spoton.berkeley.edu/SPTGUI/docs/latest#input-formats
    // frame, t (seconds), trajectory (trace id), x (um), y (um)

    try (BufferedWriter out =
        Files.newBufferedWriter(Paths.get(settings.directory, results.getName() + ".csv"))) {
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
          throw new RuntimeException(ex);
        }
      });
    } catch (final IOException | RuntimeException ex) {
      handleException(ex);
    }
  }

  @SuppressWarnings("resource")
  private void exportAnaDda(MemoryPeakResults results) {
    // anaDDA list of tracked localisations file format:
    // https://github.com/HohlbeinLab/anaDDA
    // Matlab matrix file.
    // 5 columns for n rows of localisations
    // 1. x coordinate (μm)
    // 2. y coordinate (μm)
    // 3. frame number
    // 4. track id
    // 5. frame time (s)

    // Count the number of localisations including start/end frames
    final Counter row = new Counter();
    results.forEach((PeakResultProcedure) result -> {
      row.increment(getLength(result));
    });

    // Create the matrix
    final int rows = row.getCount();
    final Matrix out = Mat5.newMatrix(rows, 5);
    // Set up column offsets
    final int col1 = rows * 1;
    final int col2 = rows * 2;
    final int col3 = rows * 3;
    final int col4 = rows * 4;

    // Frame time in seconds. This is not the frame time point converted to seconds
    // but the exposure duration of the frame.
    final double frameTime = results.getCalibrationReader().getExposureTime() / 1000;

    row.reset();
    results.forEach(DistanceUnit.UM, (XyrResultProcedure) (x, y, result) -> {
      if (result.hasEndFrame()) {
        for (int t = result.getFrame(); t <= result.getEndFrame(); t++) {
          final int index = row.getAndIncrement();
          out.setDouble(index, x);
          out.setDouble(index + col1, y);
          out.setDouble(index + col2, t);
          out.setDouble(index + col3, result.getId());
          out.setDouble(index + col4, frameTime);
        }
      } else {
        // Column major index: row + rows * col
        final int index = row.getAndIncrement();
        out.setDouble(index, x);
        out.setDouble(index + col1, y);
        out.setDouble(index + col2, result.getFrame());
        out.setDouble(index + col3, result.getId());
        out.setDouble(index + col4, frameTime);
      }
    });

    try (MatFile matFile = Mat5.newMatFile()) {
      matFile.addArray("tracks", out);
      Mat5.writeToFile(matFile, Paths.get(settings.directory, results.getName() + ".mat").toFile());
    } catch (final IOException ex) {
      handleException(ex);
    }
  }

  private static void handleException(Exception ex) {
    Logger.getLogger(TraceExporter.class.getName()).log(Level.SEVERE, "Failed to export", ex);
  }
}
