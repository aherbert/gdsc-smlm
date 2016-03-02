package gdsc.smlm.ij.plugins;

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
import gdsc.smlm.results.Calibration;
import gdsc.smlm.results.MemoryPeakResults;
import ij.IJ;
import ij.gui.GenericDialog;
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
		PluginTracker.recordPlugin(this.getClass(), arg);
		
		if (!showInputDialog())
			return;

		MemoryPeakResults results = ResultsManager.loadInputResults(inputOption, false);
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

		GenericDialog gd = new GenericDialog(TITLE);
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
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		Calibration calibration = results.getCalibration();

		gd.addStringField("Name", results.getName(), Math.max(Math.min(results.getName().length(), 60), 20));
		gd.addCheckbox("Update_all_linked_results", updateAll);
		gd.addNumericField("Calibration (nm/px)", calibration.nmPerPixel, 2);
		gd.addNumericField("Gain (ADU/photon)", calibration.gain, 2);
		gd.addCheckbox("EM-CCD", calibration.emCCD);
		gd.addNumericField("Exposure_time (ms)", calibration.exposureTime, 2);
		gd.addNumericField("Camera_bias (ADUs)", calibration.bias, 2);
		gd.addNumericField("Read_noise (ADUs)", calibration.readNoise, 2);
		
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
		
		updateAll = gd.getNextBoolean();

		// Calibration is stored as a reference. 
		// To avoid changing all datasets with the same calibration we create a copy
		
		if (!updateAll)
			calibration = calibration.clone();
		
		calibration.nmPerPixel = Math.abs(gd.getNextNumber());
		calibration.gain = Math.abs(gd.getNextNumber());
		calibration.emCCD = gd.getNextBoolean();
		calibration.exposureTime = Math.abs(gd.getNextNumber());
		calibration.bias = Math.abs(gd.getNextNumber());
		calibration.readNoise = Math.abs(gd.getNextNumber());
		
		if (!updateAll)
			results.setCalibration(calibration);

		return true;
	}
}
