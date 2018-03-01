package gdsc.smlm.ij.plugins;

import gdsc.core.data.DataException;
import gdsc.core.data.utils.TypeConverter;
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
	private static String inputOption = "";
	private static double dx = 0;
	private static double dy = 0;
	private static double dz = 0;
	private static DistanceUnit distanceUnit = DistanceUnit.PIXEL;

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

		// Show a dialog allowing the results set to be filtered
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addMessage("Select a dataset to translate");
		ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);
		gd.addNumericField("x", dx, 3);
		gd.addNumericField("y", dy, 3);
		gd.addNumericField("z", dz, 3);
		gd.addChoice("Distance_unit", SettingsManager.getDistanceUnitNames(), distanceUnit.getNumber());
		gd.showDialog();
		if (gd.wasCanceled())
			return;

		inputOption = ResultsManager.getInputSource(gd);
		dx = gd.getNextNumber();
		dy = gd.getNextNumber();
		dz = gd.getNextNumber();
		distanceUnit = DistanceUnit.forNumber(gd.getNextChoiceIndex());

		MemoryPeakResults results = ResultsManager.loadInputResults(inputOption, false, null, null);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			return;
		}

		TypeConverter<DistanceUnit> c;
		try
		{
			c = results.getDistanceConverter(distanceUnit);
		}
		catch (DataException e)
		{
			IJ.error(TITLE, "Unit conversion error: " + e.getMessage());
			return;
		}

		final float x = (float) c.convertBack(dx);
		final float y = (float) c.convertBack(dy);
		final float z = (float) c.convertBack(dz);

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