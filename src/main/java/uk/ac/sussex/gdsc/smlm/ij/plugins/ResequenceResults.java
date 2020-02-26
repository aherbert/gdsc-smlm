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
import ij.plugin.PlugIn;
import java.util.concurrent.atomic.AtomicReference;
import uk.ac.sussex.gdsc.core.ij.SimpleImageJTrackProgress;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.logging.TrackProgress;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;

/**
 * Updates the frame numbers on results that are stored in memory.
 */
public class ResequenceResults implements PlugIn {
  private static final String TITLE = "Resequence Results";

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    String inputOption;
    int start;
    int block;
    int skip;
    boolean logMapping;

    Settings() {
      // Set defaults
      inputOption = "";
      start = 1;
      block = 1;
    }

    Settings(Settings source) {
      inputOption = source.inputOption;
      start = source.start;
      block = source.block;
      skip = source.skip;
      logMapping = source.logMapping;
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

    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "There are no fitting results in memory");
      return;
    }

    if (!showDialog()) {
      return;
    }

    final MemoryPeakResults results =
        ResultsManager.loadInputResults(settings.inputOption, true, null, null);
    if (MemoryPeakResults.isEmpty(results)) {
      IJ.error(TITLE, "No results could be loaded");
      return;
    }

    if (resequenceResults(results, settings.start, settings.block, settings.skip,
        (settings.logMapping) ? SimpleImageJTrackProgress.getInstance() : null)) {
      IJ.showStatus("Resequenced " + results.getName());
    }
  }

  private boolean showDialog() {
    settings = Settings.load();

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);

    gd.addMessage("Resequence the results in memory (assumed to be continuous from 1).\n"
        + "Describe the regular repeat of the original image:\n"
        + "Start = The first frame that contained the data\n"
        + "Block = The number of continuous frames containing data\n"
        + "Skip = The number of continuous frames to ignore before the next data\n \n"
        + "E.G. 2:9:1 = Data was imaged from frame 2 for 9 frames, 1 frame to ignore,"
        + " then repeat.");

    ResultsManager.addInput(gd, settings.inputOption, InputSource.MEMORY);
    gd.addNumericField("Start", settings.start, 0);
    gd.addNumericField("Block", settings.block, 0);
    gd.addNumericField("Skip", settings.skip, 0);
    gd.addCheckbox("Log_mapping", settings.logMapping);

    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }

    settings.inputOption = ResultsManager.getInputSource(gd);
    settings.start = (int) gd.getNextNumber();
    settings.block = (int) gd.getNextNumber();
    settings.skip = (int) gd.getNextNumber();
    settings.logMapping = gd.getNextBoolean();
    settings.save();

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Start", settings.start);
      ParameterUtils.isAboveZero("Block", settings.block);
      ParameterUtils.isPositive("Skip", settings.skip);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return true;
  }

  private static class ResequencePeakResultProcedure implements PeakResultProcedure {
    final int start;
    final TrackProgress tracker;
    final int block;
    final int skip;

    ResequencePeakResultProcedure(int start, TrackProgress tracker, int block, int skip) {
      this.start = start;
      this.tracker = tracker;
      this.block = block;
      this.skip = skip;
    }

    @Override
    public void execute(PeakResult result) {
      int time = 1; // The current frame in the results
      int mapped = start; // The mapped frame in the results
      int blockSize = 1; // The current block size

      boolean print = (tracker != null);

      if (time != result.getFrame()) {
        // Update the mapped position
        while (time < result.getFrame()) {
          // Move to the next position
          mapped++;

          // Check if this move will make the current block too large
          if (++blockSize > block) {
            // Skip
            mapped += skip;
            blockSize = 1;
          }

          time++;
        }

        time = result.getFrame();
        print = (tracker != null);
      }

      result.setFrame(mapped);

      if (print) {
        tracker.log("Map %d -> %d", time, mapped);
      }
    }
  }

  /**
   * Resequence the results for the original imaging sequence provided. Results are assumed to be
   * continuous from 1.
   *
   * @param results the results
   * @param start The first frame that contained the data
   * @param block The number of continuous frames containing data
   * @param skip The number of continuous frames to ignore before the next data
   * @param tracker Used to report the mapping
   * @return true, if successful
   */
  private static boolean resequenceResults(MemoryPeakResults results, final int start,
      final int block, final int skip, final TrackProgress tracker) {
    if (MemoryPeakResults.isEmpty(results)) {
      return false;
    }

    results.sort();

    // Assume the results start from frame 1 (or above)
    if (results.getFirstFrame() < 1) {
      return false;
    }

    results.forEach(new ResequencePeakResultProcedure(start, tracker, block, skip));

    return true;
  }
}
