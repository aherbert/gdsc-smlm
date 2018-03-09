package gdsc.smlm.ij.plugins;

import gdsc.core.data.DataException;
import gdsc.core.data.utils.TypeConverter;
import gdsc.smlm.data.config.GUIProtos.TranslateResultsSettings;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;

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
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.procedures.PeakResultProcedure;
import ij.IJ;
import ij.gui.ExtendedGenericDialog;
import ij.plugin.PlugIn;

/**
 * Translate fit results. This can be used if the reference frame of the results is incorrect.
 */
public class TranslateResults implements PlugIn
{
	private static final String TITLE = "Translate Results";

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (MemoryPeakResults.isMemoryEmpty())
		{
			IJ.error(TITLE, "There are no fitting results in memory");
			return;
		}
		
		TranslateResultsSettings.Builder settings = SettingsManager.readTranslateResultsSettings(0).toBuilder();

		// Show a dialog allowing the results set to be filtered
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addMessage("Select a dataset to translate");
		ResultsManager.addInput(gd, settings.getInputOption(), InputSource.MEMORY);
		gd.addNumericField("x", settings.getDx(), 3);
		gd.addNumericField("y", settings.getDy(), 3);
		gd.addNumericField("z", settings.getDz(), 3);
		gd.addChoice("Distance_unit", SettingsManager.getDistanceUnitNames(), settings.getDistanceUnitValue());
		gd.showDialog();
		if (gd.wasCanceled())
			return;

		settings.setInputOption(ResultsManager.getInputSource(gd));
		settings.setDx(gd.getNextNumber());
		settings.setDy(gd.getNextNumber());
		settings.setDz(gd.getNextNumber());
		settings.setDistanceUnitValue(gd.getNextChoiceIndex());
		
		SettingsManager.writeSettings(settings);		
		
		MemoryPeakResults results = ResultsManager.loadInputResults(settings.getInputOption(), false, null, null);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			return;
		}

		TypeConverter<DistanceUnit> c;
		try
		{
			c = results.getDistanceConverter(settings.getDistanceUnit());
		}
		catch (DataException e)
		{
			IJ.error(TITLE, "Unit conversion error: " + e.getMessage());
			return;
		}

		final float x = (float) c.convertBack(settings.getDx());
		final float y = (float) c.convertBack(settings.getDy());
		final float z = (float) c.convertBack(settings.getDz());

		// Reset the 2D bounds
		if (x != 0 || y != 0)
			results.setBounds(null);

		results.forEach(new PeakResultProcedure()
		{
			public void execute(PeakResult peakResult)
			{
				// Requires a direct reference!
				float[] params = peakResult.getParameters();
				params[PeakResult.X] += x;
				params[PeakResult.Y] += y;
				params[PeakResult.Z] += z;
			}
		});
	}
}