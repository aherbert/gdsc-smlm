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
import java.util.EnumSet;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.config.UnitHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.TimeUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.SaveLocalisationsSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.SaveLocalisationsSettings.Builder;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;

/**
 * Saves localisations to a user-specified delimited text file.
 */
public class SaveLocalisations implements PlugIn {
  private static final String TITLE = "Save Localisations";

  /**
   * Lazy load the TimeUnit names and values.
   */
  private static class TimeUnitSaveer {
    // Time units for the exposure time cannot be in frames as this makes no sense
    private static final String[] TIME_UNITS;
    private static final TimeUnit[] TIME_UNIT_VALUES;

    static {
      final EnumSet<TimeUnit> set = EnumSet.allOf(TimeUnit.class);
      set.remove(TimeUnit.FRAME);
      set.remove(TimeUnit.UNRECOGNIZED);
      TIME_UNITS = new String[set.size()];
      TIME_UNIT_VALUES = new TimeUnit[set.size()];
      int index = 0;
      for (final TimeUnit t : set) {
        TIME_UNITS[index] =
            SettingsManager.getName(UnitHelper.getName(t), UnitHelper.getShortName(t));
        TIME_UNIT_VALUES[index] = t;
        index++;
      }
    }

    /**
     * Gets the time units.
     *
     * @return the time units
     */
    static String[] getTimeUnits() {
      return TIME_UNITS;
    }

    /**
     * Gets the time unit values.
     *
     * @return the time unit values
     */
    static TimeUnit[] getTimeUnitValues() {
      return TIME_UNIT_VALUES;
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "There are no fitting results in memory");
      return;
    }

    final SaveLocalisationsSettings.Builder settings =
        SettingsManager.readSaveLocalisationsSettings(0).toBuilder();

    // Get results
    if (!showInputDialog(settings)) {
      return;
    }

    SettingsManager.writeSettings(settings.build());

    final MemoryPeakResults results =
        ResultsManager.loadInputResults(settings.getInput(), false, null, null);
    if (MemoryPeakResults.isEmpty(results)) {
      IJ.error(TITLE, "No results could be loaded");
      return;
    }

    // TODO

    // Check calibration to determine if units can be mapped
    // Check PSF to obtain custom fields

    // Custom field for a PeakResult
    // It has the following properties:
    // ID
    // Name (for column headers)
    // Unit (reuse PSFParameterUnit)
    // toString() to show the summary of properties

    // Build a LinkedHashMap<String, Field> from the results
    // Show the available fields in the dialog (in order)
    // Strip non-existent fields from the current output 'format'
    // Add options to the dialog to choose units if the calibration exists.
    // Convert fields to formatters.
    // Save results.

    //
    // final String[] path = ImageJUtils.decodePath(settings.getLocalisationsFilename());
    // final OpenDialog chooser = new OpenDialog("Localisations_File", path[0], path[1]);
    // if (chooser.getFileName() == null) {
    // return;
    // }
    //
    // settings.setLocalisationsFilename(chooser.getDirectory() + chooser.getFileName());

    SettingsManager.writeSettings(settings.build());

    final String msg = "Saved " + TextUtils.pleural(results.size(), "localisation");
    IJ.showStatus(msg);
    // ImageJUtils.log(msg + " to " + settings.getLocalisationsFilename());
  }


  /**
   * Show a dialog to obtain the name of the input localisations.
   *
   * @param settings the settings
   * @return true, if successful
   */
  private static boolean showInputDialog(Builder settings) {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("save-localisations"));

    gd.addMessage("Save localisations in a user-specified text format");

    ResultsManager.addInput(gd, settings.getInput(), InputSource.MEMORY);

    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }

    settings.setInput(ResultsManager.getInputSource(gd));
    return true;
  }
}
