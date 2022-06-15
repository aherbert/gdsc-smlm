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
import ij.plugin.PlugIn;
import java.util.concurrent.atomic.AtomicReference;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.SimpleImageJTrackProgress;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.utils.FileUtils;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.TextFilePeakResults;
import uk.ac.sussex.gdsc.smlm.results.Trace;
import uk.ac.sussex.gdsc.smlm.results.TraceManager;

/**
 * Run a tracing algorithm on the peak results to trace neighbours across the frames.
 */
public class NeighbourAnalysis implements PlugIn {
  private static final String TITLE = "Neighbour Analysis";

  private MemoryPeakResults results;

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> INSTANCE = new AtomicReference<>(new Settings());

    String inputOption;
    double distanceThreshold;
    int timeThreshold;

    String filename;

    Settings() {
      // Set defaults
      inputOption = "";
      distanceThreshold = 0.6;
      timeThreshold = 1;
      filename = "";
    }

    Settings(Settings source) {
      inputOption = source.inputOption;
      distanceThreshold = source.distanceThreshold;
      timeThreshold = source.timeThreshold;
      filename = source.filename;
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

    final TraceManager manager = new TraceManager(results);

    // Run the tracing
    manager.setTracker(SimpleImageJTrackProgress.getInstance());
    final Trace[] traces =
        manager.findNeighbours(settings.distanceThreshold, settings.timeThreshold);

    saveTraces(traces);
  }

  private void saveTraces(Trace[] traces) {
    final String filename =
        ImageJUtils.getFilename("Traces_File", settings.filename);
    if (filename != null) {
      // Remove extension and replace with .xls
      settings.filename = FileUtils.replaceExtension(filename, ".xls");

      final boolean showDeviations = results.hasDeviations();
      final TextFilePeakResults traceResults =
          new TextFilePeakResults(settings.filename, showDeviations);
      traceResults.copySettings(results);
      traceResults.begin();
      if (!traceResults.isActive()) {
        IJ.error(TITLE, "Failed to write to file: " + settings.filename);
        return;
      }
      traceResults.addComment(createSettingsComment());
      for (final Trace trace : traces) {
        traceResults.addCluster(trace); // addTrace(...) does a sort on the results
      }
      traceResults.end();
    }
  }

  private String createSettingsComment() {
    return String.format("Neighbour tracing : distance-threshold = %f : time-threshold = %d",
        settings.distanceThreshold, settings.timeThreshold);
  }

  private boolean showDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("neighbour-analysis"));

    settings = Settings.load();
    ResultsManager.addInput(gd, settings.inputOption, InputSource.MEMORY);

    gd.addNumericField("Distance_Threshold (px)", settings.distanceThreshold, 4);
    gd.addNumericField("Time_Threshold (frames)", settings.timeThreshold, 0);

    gd.showDialog();

    if (gd.wasCanceled() || !readDialog(gd)) {
      return false;
    }

    // Load the results
    results = ResultsManager.loadInputResults(settings.inputOption, false, DistanceUnit.PIXEL);
    if (MemoryPeakResults.isEmpty(results)) {
      IJ.error(TITLE, "No results could be loaded");
      IJ.showStatus("");
      return false;
    }

    return true;
  }

  private boolean readDialog(ExtendedGenericDialog gd) {
    settings.inputOption = ResultsManager.getInputSource(gd);
    settings.distanceThreshold = gd.getNextNumber();
    settings.timeThreshold = (int) gd.getNextNumber();

    if (settings.distanceThreshold < 0) {
      settings.distanceThreshold = 0;
    }
    if (settings.timeThreshold < 0) {
      settings.timeThreshold = 0;
    }
    settings.save();
    if (settings.timeThreshold == 0 && settings.distanceThreshold == 0) {
      IJ.error(TITLE, "No thresholds specified");
      return false;
    }
    return !gd.invalidNumber();
  }
}
