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
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;

/**
 * Allows results held in memory to be converted to different units.
 */
public class ConvertResults implements PlugIn {
  private static final String TITLE = "Convert Results";

  private static AtomicReference<String> inputOptionRef = new AtomicReference<>("");

  private String inputOption;

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (!showInputDialog()) {
      return;
    }

    final MemoryPeakResults results =
        ResultsManager.loadInputResults(inputOption, false, null, null);
    if (MemoryPeakResults.isEmpty(results)) {
      IJ.error(TITLE, "No results could be loaded");
      return;
    }

    if (!showDialog(results)) {
      return;
    }

    IJ.showStatus("Converted " + results.getName());
  }

  private boolean showInputDialog() {
    final int size = MemoryPeakResults.countMemorySize();
    if (size == 0) {
      IJ.error(TITLE, "There are no fitting results in memory");
      return false;
    }

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("convert-results"));
    gd.addMessage("Select results to convert");

    inputOption = inputOptionRef.get();

    ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    inputOption = ResultsManager.getInputSource(gd);
    inputOptionRef.set(inputOption);

    return true;
  }

  private static boolean showDialog(MemoryPeakResults results) {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addMessage("Convert the current units for the results");
    gd.addHelp(HelpUrls.getUrl("convert-results"));

    final CalibrationReader cr = CalibrationWriter.create(results.getCalibration());

    gd.addChoice("Distance_unit", SettingsManager.getDistanceUnitNames(),
        cr.getDistanceUnitValue());
    gd.addNumericField("Calibration (nm/px)", cr.getNmPerPixel(), 2);
    gd.addChoice("Intensity_unit", SettingsManager.getIntensityUnitNames(),
        cr.getIntensityUnitValue());
    gd.addNumericField("Gain (Count/photon)", cr.getCountPerPhoton(), 2);
    gd.addChoice("Angle_unit", SettingsManager.getAngleUnitNames(), cr.getAngleUnitValue());

    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }

    final CalibrationWriter cw = results.getCalibrationWriterSafe();
    final DistanceUnit distanceUnit =
        SettingsManager.getDistanceUnitValues()[gd.getNextChoiceIndex()];
    cw.setNmPerPixel(Math.abs(gd.getNextNumber()));

    // Don't set the calibration with bad values
    if (distanceUnit.getNumber() > 0 && !(cw.getNmPerPixel() > 0)) {
      IJ.error(TITLE, "Require positive nm/pixel for conversion");
      return false;
    }

    final IntensityUnit intensityUnit =
        SettingsManager.getIntensityUnitValues()[gd.getNextChoiceIndex()];
    cw.setCountPerPhoton(Math.abs(gd.getNextNumber()));
    if (intensityUnit.getNumber() > 0 && !(cw.getCountPerPhoton() > 0)) {
      IJ.error(TITLE, "Require positive Count/photon for conversion");
      return false;
    }

    final Calibration newCalibration = cw.getCalibration();
    results.setCalibration(newCalibration);

    final AngleUnit angleUnit = SettingsManager.getAngleUnitValues()[gd.getNextChoiceIndex()];
    if (!results.convertToUnits(distanceUnit, intensityUnit, angleUnit)) {
      IJ.error(TITLE, "Conversion failed");
      return false;
    }

    return true;
  }
}
