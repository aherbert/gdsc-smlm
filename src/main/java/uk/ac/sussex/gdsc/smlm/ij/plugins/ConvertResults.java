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

import ij.IJ;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.UnitHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;

/**
 * Allows results held in memory to be converted to different units.
 */
public class ConvertResults implements PlugIn
{
    private static final String TITLE = "Convert Results";

    private static String inputOption = "";

    /*
     * (non-)
     *
     * @see ij.plugin.PlugIn#run(java.lang.String)
     */
    @Override
    public void run(String arg)
    {
        SMLMUsageTracker.recordPlugin(this.getClass(), arg);

        if (!showInputDialog())
            return;

        final MemoryPeakResults results = ResultsManager.loadInputResults(inputOption, false, null, null);
        if (results == null || results.size() == 0)
        {
            IJ.error(TITLE, "No results could be loaded");
            return;
        }

        if (!showDialog(results))
            return;

        IJ.showStatus("Converted " + results.getName());
    }

    private static boolean showInputDialog()
    {
        final int size = MemoryPeakResults.countMemorySize();
        if (size == 0)
        {
            IJ.error(TITLE, "There are no fitting results in memory");
            return false;
        }

        final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
        gd.addHelp(About.HELP_URL);
        gd.addMessage("Select results to convert");

        ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);

        gd.showDialog();

        if (gd.wasCanceled())
            return false;

        inputOption = ResultsManager.getInputSource(gd);

        return true;
    }

    private static boolean showDialog(MemoryPeakResults results)
    {
        final GenericDialog gd = new GenericDialog(TITLE);
        gd.addMessage("Convert the current units for the results");
        gd.addHelp(About.HELP_URL);

        final CalibrationReader cr = CalibrationWriter.create(results.getCalibration());

        gd.addChoice("Distance_unit", SettingsManager.getDistanceUnitNames(), UnitHelper.getName(cr.getDistanceUnit()));
        gd.addNumericField("Calibration (nm/px)", cr.getNmPerPixel(), 2);
        gd.addChoice("Intensity_unit", SettingsManager.getIntensityUnitNames(),
                UnitHelper.getName(cr.getIntensityUnit()));
        gd.addNumericField("Gain (Count/photon)", cr.getCountPerPhoton(), 2);
        gd.addChoice("Angle_unit", SettingsManager.getAngleUnitNames(), UnitHelper.getName(cr.getAngleUnit()));

        gd.showDialog();
        if (gd.wasCanceled())
            return false;

        final CalibrationWriter cw = results.getCalibrationWriterSafe();
        final DistanceUnit distanceUnit = SettingsManager.getDistanceUnitValues()[gd.getNextChoiceIndex()];
        cw.setNmPerPixel(Math.abs(gd.getNextNumber()));
        final IntensityUnit intensityUnit = SettingsManager.getIntensityUnitValues()[gd.getNextChoiceIndex()];
        cw.setCountPerPhoton(Math.abs(gd.getNextNumber()));
        final AngleUnit angleUnit = SettingsManager.getAngleUnitValues()[gd.getNextChoiceIndex()];

        // Don't set the calibration with bad values
        if (distanceUnit.getNumber() > 0 && !(cw.getNmPerPixel() > 0))
        {
            IJ.error(TITLE, "Require positive nm/pixel for conversion");
            return false;
        }
        if (intensityUnit.getNumber() > 0 && !(cw.getCountPerPhoton() > 0))
        {
            IJ.error(TITLE, "Require positive Count/photon for conversion");
            return false;
        }

        final Calibration newCalibration = cw.getCalibration();
        results.setCalibration(newCalibration);

        if (!results.convertToUnits(distanceUnit, intensityUnit, angleUnit))
        {
            IJ.error(TITLE, "Conversion failed");
            return false;
        }

        return true;
    }
}
