package gdsc.smlm.ij.plugins;

import gdsc.smlm.data.config.CalibrationProtos.Calibration;
import gdsc.smlm.data.config.CalibrationWriter;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.results.MemoryPeakResults;
import ij.IJ;
import ij.gui.ExtendedGenericDialog;
import ij.plugin.PlugIn;

/**
 * Allows results held in memory to be calibrated.
 */
public class CalibrateResults implements PlugIn
{
	private static final String TITLE = "Calibrate Results";

	private static String inputOption = "";
	private static boolean updateAll = false;

	/*
	 * (non-)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (!showInputDialog())
			return;

		MemoryPeakResults results = ResultsManager.loadInputResults(inputOption, false, null, null);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			return;
		}

		if (!showDialog(results))
			return;

		IJ.showStatus("Calibrated " + results.getName());
	}

	private boolean showInputDialog()
	{
		int size = MemoryPeakResults.countMemorySize();
		if (size == 0)
		{
			IJ.error(TITLE, "There are no fitting results in memory");
			return false;
		}

		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);
		gd.addMessage("Select results to calibrate");

		ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		inputOption = ResultsManager.getInputSource(gd);

		return true;
	}

	private boolean showDialog(MemoryPeakResults results)
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		Calibration oldCalibration = results.getCalibration();

		boolean existingCalibration = oldCalibration != null;
		if (!existingCalibration)
		{
			gd.addMessage("No calibration found, using defaults");
		}

		final CalibrationWriter cw = results.getCalibrationWriterSafe();

		gd.addStringField("Name", results.getName(), Math.max(Math.min(results.getName().length(), 60), 20));
		if (existingCalibration)
			gd.addCheckbox("Update_all_linked_results", updateAll);
		PeakFit.addCameraOptions(gd, PeakFit.FLAG_QUANTUM_EFFICIENCY, cw);
		gd.addNumericField("Calibration (nm/px)", cw.getNmPerPixel(), 2);
		gd.addNumericField("Exposure_time (ms)", cw.getExposureTime(), 2);

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		String name = gd.getNextString();
		if (!results.getName().equals(name))
		{
			// If the name is changed then remove and add back to memory
			results = MemoryPeakResults.removeResults(results.getName());
			if (results != null)
			{
				results.setName(name);
				MemoryPeakResults.addResults(results);
			}
		}

		if (existingCalibration)
		{
			updateAll = gd.getNextBoolean();
		}

		cw.setCameraType(SettingsManager.getCameraTypeValues()[gd.getNextChoiceIndex()]);
		cw.setNmPerPixel(Math.abs(gd.getNextNumber()));
		cw.setExposureTime(Math.abs(gd.getNextNumber()));

		gd.collectOptions();

		Calibration newCalibration = cw.getCalibration();
		results.setCalibration(newCalibration);

		if (updateAll)
		{
			// Calibration is stored as a reference to an immutable object.
			// Update any in memory results with the same object.
			// Note that if any plugins have modified the calibration (rather 
			// than just copy it through to a new results set) then they will
			// be missed.
			for (MemoryPeakResults r : MemoryPeakResults.getAllResults())
			{
				if (r.getCalibration() == oldCalibration)
					r.setCalibration(newCalibration);
			}
		}

		return true;
	}
}
