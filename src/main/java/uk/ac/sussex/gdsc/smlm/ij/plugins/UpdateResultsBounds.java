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

import ij.IJ;
import ij.plugin.PlugIn;
import java.awt.Rectangle;
import java.util.concurrent.atomic.AtomicReference;
import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;

/**
 * Allows results held in memory to be calibrated.
 */
public class UpdateResultsBounds implements PlugIn {
  private static final String TITLE = "Update Results Bounds";

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

    Settings() {
      inputOption = "";
    }

    Settings(Settings source) {
      inputOption = source.inputOption;
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

    if (!showInputDialog()) {
      return;
    }

    final MemoryPeakResults results =
        ResultsManager.loadInputResults(settings.inputOption, true, null, null);
    if (MemoryPeakResults.isEmpty(results)) {
      IJ.error(TITLE, "No results could be loaded");
      return;
    }

    if (!showDialog(results)) {
      return;
    }

    IJ.showStatus("Updated " + results.getName());
  }

  private boolean showInputDialog() {
    final int size = MemoryPeakResults.countMemorySize();
    if (size == 0) {
      IJ.error(TITLE, "There are no fitting results in memory");
      return false;
    }

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);
    gd.addMessage("Select results to update");

    settings = Settings.load();
    ResultsManager.addInput(gd, settings.inputOption, InputSource.MEMORY);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.inputOption = ResultsManager.getInputSource(gd);
    settings.save();

    return true;
  }

  private static boolean showDialog(MemoryPeakResults results) {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);

    // Force computation of bounds
    Rectangle currentBounds = results.getBounds();
    results.setBounds(null);

    Rectangle autoBounds;
    try {
      autoBounds = results.getBounds(true);
    } catch (final DataException ex) {
      IJ.error(TITLE, "No calibration found to convert to pixel units");
      return false;
    } finally {
      // Reset after forcing computation
      if (currentBounds != null) {
        results.setBounds(currentBounds);
      }
    }

    // Re-acquire the bounds (either the existing bounds or the auto-bounds)
    currentBounds = results.getBounds();

    gd.addMessage(TextUtils
        .wrap("Set the bounds of the original source image.\n \nAuto-bounds = " + format(autoBounds)
            + "\n \nThe new bounds will be the union of the auto-bounds and specified bounds "
            + "to ensure all data is within the bounds.", 80));
    gd.addNumericField("Min_x", currentBounds.x, 0, 6, "px");
    gd.addNumericField("Min_y", currentBounds.y, 0, 6, "px");
    gd.addNumericField("Width", currentBounds.width, 0, 6, "px");
    gd.addNumericField("Height", currentBounds.height, 0, 6, "px");

    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }

    final int x = (int) gd.getNextNumber();
    final int y = (int) gd.getNextNumber();
    final int width = (int) gd.getNextNumber();
    final int height = (int) gd.getNextNumber();

    // Check the bounds are not smaller than the auto-bounds
    final Rectangle newBounds = new Rectangle(x, y, width, height).union(autoBounds);
    if (newBounds.isEmpty()) {
      IJ.error(TITLE, "New bounds are not valid: " + format(newBounds));
      return false;
    }

    results.setBounds(newBounds);

    return true;
  }

  private static String format(Rectangle bounds) {
    return "[x=" + bounds.x + ",y=" + bounds.y + ",width=" + bounds.width + ",height="
        + bounds.height + "]";
  }
}
