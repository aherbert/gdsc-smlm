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
import ij.plugin.PlugIn;
import java.util.concurrent.atomic.AtomicReference;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsSettings;
import uk.ac.sussex.gdsc.smlm.ij.gui.TraceDataTableModel;
import uk.ac.sussex.gdsc.smlm.ij.gui.TraceDataTableModelFrame;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.MemoryResultsList;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;

/**
 * Plugin to display traced datasets.
 */
public class TraceViewer implements PlugIn {
  private static final String TITLE = "Trace Viewer";

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> INSTANCE = new AtomicReference<>(new Settings());

    String input;

    Settings() {
      // Set defaults
      input = "";
    }

    Settings(Settings source) {
      input = source.input;
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

    final MemoryResultsList items = new MemoryResultsList(MemoryPeakResults::hasId);

    if (items.isEmpty()) {
      IJ.error(TITLE, "No traced localisations in memory");
      return;
    }

    // Get input options
    if (!showDialog()) {
      return;
    }

    final MemoryPeakResults results =
        ResultsManager.loadInputResults(settings.input, false, null, null);
    if (results == null) {
      IJ.error(TITLE, "Failed to load results: " + settings.input);
      return;
    }

    // Create a table for the traces
    final TraceDataTableModel model = new TraceDataTableModel(results,
        SettingsManager.readResultsSettings(0).getResultsTableSettings());
    final TraceDataTableModelFrame frame = new TraceDataTableModelFrame(model);
    frame.setTitle(results.getName());
    frame.setVisible(true);
    // Save changes to the interactive table settings
    model.addSettingsUpdatedAction(s -> {
      // Load the latest settings and save the table options
      final ResultsSettings.Builder resultsSettings =
          SettingsManager.readResultsSettings(0).toBuilder();
      resultsSettings.setResultsTableSettings(s);
      SettingsManager.writeSettings(resultsSettings.build());
    });
  }

  private boolean showDialog() {
    settings = Settings.load();
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    ResultsManager.addInput(gd, settings.input, InputSource.MEMORY);
    gd.addHelp(HelpUrls.getUrl("trace-viewer"));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }
    settings.input = ResultsManager.getInputSource(gd);
    settings.save();
    return true;
  }
}
