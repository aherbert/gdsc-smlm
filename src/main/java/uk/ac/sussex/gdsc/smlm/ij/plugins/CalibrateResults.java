/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
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

import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;

import ij.IJ;
import ij.plugin.PlugIn;

/**
 * Allows results held in memory to be calibrated.
 */
public class CalibrateResults implements PlugIn {
  private static final String TITLE = "Calibrate Results";

  private static String inputOption = "";
  private static boolean updateAll;

  @Override
  public void run(String arg) {
    SMLMUsageTracker.recordPlugin(this.getClass(), arg);

    if (!showInputDialog()) {
      return;
    }

    final MemoryPeakResults results =
        ResultsManager.loadInputResults(inputOption, false, null, null);
    if (results == null || results.size() == 0) {
      IJ.error(TITLE, "No results could be loaded");
      return;
    }

    if (!showDialog(results)) {
      return;
    }

    IJ.showStatus("Calibrated " + results.getName());
  }

  private static boolean showInputDialog() {
    final int size = MemoryPeakResults.countMemorySize();
    if (size == 0) {
      IJ.error(TITLE, "There are no fitting results in memory");
      return false;
    }

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);
    gd.addMessage("Select results to calibrate");

    ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    inputOption = ResultsManager.getInputSource(gd);

    return true;
  }

  private static boolean showDialog(MemoryPeakResults results) {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);

    final Calibration oldCalibration = results.getCalibration();

    final boolean existingCalibration = oldCalibration != null;
    if (!existingCalibration) {
      gd.addMessage("No calibration found, using defaults");
    }

    final CalibrationWriter cw = results.getCalibrationWriterSafe();

    gd.addStringField("Name", results.getName(),
        Math.max(Math.min(results.getName().length(), 60), 20));
    if (existingCalibration) {
      gd.addCheckbox("Update_all_linked_results", updateAll);
    }
    PeakFit.addCameraOptions(gd, PeakFit.FLAG_QUANTUM_EFFICIENCY, cw);
    gd.addNumericField("Calibration (nm/px)", cw.getNmPerPixel(), 2);
    gd.addNumericField("Exposure_time (ms)", cw.getExposureTime(), 2);

    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }

    final String name = gd.getNextString();
    if (!results.getName().equals(name)) {
      // If the name is changed then remove and add back to memory
      final MemoryPeakResults existingResults = MemoryPeakResults.removeResults(results.getName());
      if (existingResults != null) {
        results = existingResults;
        results.setName(name);
        MemoryPeakResults.addResults(results);
      }
    }

    if (existingCalibration) {
      updateAll = gd.getNextBoolean();
    }

    cw.setCameraType(SettingsManager.getCameraTypeValues()[gd.getNextChoiceIndex()]);
    cw.setNmPerPixel(Math.abs(gd.getNextNumber()));
    cw.setExposureTime(Math.abs(gd.getNextNumber()));

    gd.collectOptions();

    final Calibration newCalibration = cw.getCalibration();
    results.setCalibration(newCalibration);

    if (updateAll) {
      // Calibration is stored as a reference to an immutable object.
      // Update any in memory results with the same object.
      // Note that if any plugins have modified the calibration (rather
      // than just copy it through to a new results set) then they will
      // be missed.
      for (final MemoryPeakResults r : MemoryPeakResults.getAllResults()) {
        if (r.getCalibration() == oldCalibration) {
          r.setCalibration(newCalibration);
        }
      }
    }

    return true;
  }
}
