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
import ij.gui.ExtendedGenericDialog;
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
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);
		
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
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		Calibration calibration = results.getCalibration();
		
		boolean newCalibration = false;
		if (calibration == null)
		{
			newCalibration = true;
			calibration = new Calibration();
			gd.addMessage("No calibration found, using defaults");
		}

		gd.addStringField("Name", results.getName(), Math.max(Math.min(results.getName().length(), 60), 20));
		if (!newCalibration)
			gd.addCheckbox("Update_all_linked_results", updateAll);
		gd.addNumericField("Calibration (nm/px)", calibration.getNmPerPixel(), 2);
		gd.addNumericField("Gain (ADU/photon)", calibration.getGain(), 2);
		gd.addCheckbox("EM-CCD", calibration.isEmCCD());
		gd.addNumericField("Exposure_time (ms)", calibration.getExposureTime(), 2);
		gd.addNumericField("Camera_bias (ADUs)", calibration.getBias(), 2);
		gd.addNumericField("Read_noise (ADUs)", calibration.getReadNoise(), 2);
		gd.addNumericField("Amplification (ADUs/electron)", calibration.getAmplification(), 2);
		
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
		
		if (!newCalibration)
		{
			updateAll = gd.getNextBoolean();
			
			// Calibration is stored as a reference. 
			// To avoid changing all datasets with the same calibration we create a copy
			if (!updateAll)
			{
				newCalibration = true;
				calibration = calibration.clone();				
			}
		}
		
		calibration.setNmPerPixel(Math.abs(gd.getNextNumber()));
		calibration.setGain(Math.abs(gd.getNextNumber()));
		calibration.setEmCCD(gd.getNextBoolean());
		calibration.setExposureTime(Math.abs(gd.getNextNumber()));
		calibration.setBias(Math.abs(gd.getNextNumber()));
		calibration.setReadNoise(Math.abs(gd.getNextNumber()));
		calibration.setAmplification(Math.abs(gd.getNextNumber()));
		
		if (newCalibration)
			results.setCalibration(calibration);

		return true;
	}
}
