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
import ij.gui.Plot;
import ij.plugin.PlugIn;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
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
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;
import uk.ac.sussex.gdsc.core.utils.rng.UniformRandomProviders;
import uk.ac.sussex.gdsc.smlm.data.NamedObject;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.UnitConverterUtils;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.TimeUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.MemoryResultsList;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.AttributePeakResult;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.count.Counter;
import uk.ac.sussex.gdsc.smlm.results.count.FrameCounter;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.XyrResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.XyzrResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.sort.IdFramePeakResultComparator;
import us.hebi.matlab.mat.format.Mat5;
import us.hebi.matlab.mat.types.Cell;
import us.hebi.matlab.mat.types.MatFile;
import us.hebi.matlab.mat.types.MatlabType;
import us.hebi.matlab.mat.types.Matrix;
import us.hebi.matlab.mat.types.Struct;

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
      // We do not want capitalisation of 'ana' or 'vb'
      formatNames[ExportFormat.ANA_DNA.ordinal()] = ExportFormat.ANA_DNA.getName();
      formatNames[ExportFormat.VB_SPT.ordinal()] = ExportFormat.VB_SPT.getName();
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
    SPOT_ON("Spot-On"), ANA_DNA("anaDDA"), VB_SPT("vbSPT"), NOBIAS("NOBIAS");

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
    final IntOpenHashSet remove = new IntOpenHashSet();
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

    if (settings.format == 3) {
      exportNobias(results);
    } else if (settings.format == 2) {
      exportVbSpt(results);
    } else if (settings.format == 1) {
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

  @SuppressWarnings("resource")
  private void exportVbSpt(MemoryPeakResults results) {
    // vbSPT file format:
    // https://sourceforge.net/projects/vbspt/
    // Matlab matrix file (.mat) containing at least one variable that is a cell
    // array where each element, representing a trajectory, is a matrix
    // where the rows define the coordinates in one, two or three dimensions
    // in subsequent timesteps. The number of dimensions to be used for the
    // analysis will be set by the runinputfile.
    // The units are arbitrary but vbSPT starting estimates must be in the same
    // units. Either nm or μm are recommended.
    // 3 columns for n rows of localisations
    // 1. x coordinate (μm)
    // 2. y coordinate (μm)
    // 3. z coordinate (μm)
    //
    // Note: An extra column is added containing the frame. This allows results to
    // be uniquely identified using frame,x,y,z

    // Count the IDs. Each new result ID will increment the count.
    final FrameCounter idCounter = new FrameCounter(results.getFirst().getId() - 1);
    results.forEach((PeakResultProcedure) result -> {
      if (idCounter.advance(result.getId())) {
        idCounter.increment();
      }
    });

    // Create the cell array as 1xN
    final Cell out = Mat5.newCell(1, idCounter.getCount());

    // This will reset the counter to zero and ensure the current frame does not match
    // in the event of a single track
    idCounter.advanceAndReset(idCounter.currentFrame() + 1);

    final boolean is3d = results.is3D();

    // Write the tracks
    final LocalList<double[]> list = new LocalList<>();
    results.forEach(DistanceUnit.UM, (XyzrResultProcedure) (x, y, z, result) -> {
      if (idCounter.advance(result.getId())) {
        addVbSptTrack(out, idCounter.getCount() - 1, list, is3d);
        idCounter.increment();
        list.clear();
      }
      list.add(new double[] {x, y, z, result.getFrame()});
    });
    addVbSptTrack(out, idCounter.getCount() - 1, list, is3d);

    try (MatFile matFile = Mat5.newMatFile()) {
      matFile.addArray("tracks", out);
      Mat5.writeToFile(matFile, Paths.get(settings.directory, results.getName() + ".mat").toFile());
    } catch (final IOException ex) {
      handleException(ex);
    }
  }

  /**
   * Adds the track to the cell at the given row index.
   *
   * @param cell the output cell
   * @param index the index
   * @param list the list
   * @param is3d true if the data is 3D, otherwise write 2D data
   */
  @SuppressWarnings("resource")
  private static void addVbSptTrack(Cell cell, int index, LocalList<double[]> list, boolean is3d) {
    if (list.isEmpty()) {
      return;
    }
    // Create the matrix
    final int rows = list.size();
    final Matrix m = Mat5.newMatrix(rows, is3d ? 4 : 3);
    // Set up column offsets
    final int col1 = rows * 1;
    final int col2 = rows * 2;
    final int col3 = rows * 3;
    for (int i = 0; i < rows; i++) {
      final double[] xyz = list.unsafeGet(i);
      m.setDouble(i, xyz[0]);
      m.setDouble(i + col1, xyz[1]);
      if (is3d) {
        m.setDouble(i + col2, xyz[2]);
        m.setDouble(i + col3, xyz[3]);
      } else {
        m.setDouble(i + col2, xyz[3]);
      }
    }
    cell.set(index, m);
  }

  @SuppressWarnings("resource")
  private void exportNobias(MemoryPeakResults results) {
    // NOBIAS file format:
    // https://github.com/BiteenMatlab/NOBIAS
    // Matlab matrix file (.mat) input data variable should be a structure variable with
    // at least two fields: ‘obs’ and ‘TrID’. ‘obs’ is a 2 × T matrix of observations where T is the
    // total amount of steps and 2 indicates the dimensionality (2D tracks). Each element of ‘obs’
    // is the step size in that dimension. ‘TrID’ is the trajectory ID number, which denotes the
    // track that step is from. It should be a 1 × T integer vector. To do the motion blur
    // correction, the input data must also have a third field, ‘obs_corr’ (edited), that denotes
    // the correlation steps.

    // Note: From reading NOBIAS_get_simu_data_corrsteps.m
    // the 'obs_corr' is the product of two consecutive steps. The final step is nan.

    // The NOBIAS ‘Params’ structure includes two parameters, ‘Params.frametime’ and
    // ‘Params.pixelsize’, that depend on your experimental settings. Please indicate your imaging
    // frame time in units of seconds and your imaging pixel size in units of μm.

    // Count the IDs. Each new result ID will increment the count.
    final FrameCounter idCounter = new FrameCounter(results.getFirst().getId() - 1);
    results.forEach((PeakResultProcedure) result -> {
      if (idCounter.advance(result.getId())) {
        idCounter.increment();
      }
    });

    // This will reset the counter to zero and ensure the current frame does not match
    // in the event of a single track
    idCounter.advanceAndReset(idCounter.currentFrame() + 1);

    // Collect the jumps for the tracks
    final LocalList<double[]> list = new LocalList<>(results.size());
    final double[] last = new double[3];
    results.forEach(DistanceUnit.PIXEL, (XyrResultProcedure) (x, y, result) -> {
      if (idCounter.advance(result.getId())) {
        last[0] = x;
        last[1] = y;
        // Use the original track ID.
        // NOBIAS is robust to this as it checks for a change in ID using: if (TrID(t)~=TrID(t-1))
        last[2] = result.getId();
        //last[2] = idCounter.incrementAndGet();
      } else {
        last[0] = x - last[0];
        last[1] = y - last[1];
        list.add(last.clone());
        last[0] = x;
        last[1] = y;
      }
    });

    // Create the 'data' cell array as 1x1
    final Struct data = Mat5.newStruct(1, 1);

    final int rows = list.size();

    final Matrix trid = Mat5.newMatrix(rows, 1, MatlabType.Int32);
    for (int i = 0; i < rows; i++) {
      final double[] xyz = list.unsafeGet(i);
      trid.setInt(i, (int) xyz[2]);
    }
    data.set("TrID", trid);

    final Matrix obs = Mat5.newMatrix(2, rows);
    for (int i = 0; i < rows; i++) {
      final double[] xyz = list.unsafeGet(i);
      obs.setDouble(2 * i, xyz[0]);
      obs.setDouble(2 * i + 1, xyz[1]);
    }
    data.set("obs", obs);

    final Matrix obsCorr = Mat5.newMatrix(2, rows);
    if (rows != 0) {
      double[] xyz2 = list.unsafeGet(0);
      for (int i = 1; i < rows; i++) {
        final double[] xyz = list.unsafeGet(i);
        obsCorr.setDouble(2 * i - 2, xyz[0] * xyz2[0]);
        obsCorr.setDouble(2 * i - 1, xyz[1] * xyz2[1]);
        xyz2 = xyz;
      }
      obsCorr.setDouble(2 * rows - 2, Double.NaN);
      obsCorr.setDouble(2 * rows - 1, Double.NaN);
    }
    data.set("obs_corr", obsCorr);

    Struct params = null;
    if (results.hasCalibration()) {
      params = Mat5.newStruct(1, 1);
      final CalibrationReader cr = results.getCalibrationReader();
      params.set("frametime", Mat5.newScalar(cr.getExposureTime() / 1000));
      params.set("pixelsize", Mat5.newScalar(cr.getNmPerPixel() / 1000));
    }

    try (MatFile matFile = Mat5.newMatFile()) {
      matFile.addArray("data", data);
      if (params != null) {
        matFile.addArray("Params", params);
      }
      Mat5.writeToFile(matFile, Paths.get(settings.directory, results.getName() + ".mat").toFile());
    } catch (final IOException ex) {
      handleException(ex);
    }
  }

  private static void handleException(Exception ex) {
    Logger.getLogger(TraceExporter.class.getName()).log(Level.SEVERE, "Failed to export", ex);
  }
}
