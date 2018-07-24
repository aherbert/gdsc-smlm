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
import ij.plugin.PlugIn;
import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.TranslateResultsSettings;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;

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
	@Override
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (MemoryPeakResults.isMemoryEmpty())
		{
			IJ.error(TITLE, "There are no fitting results in memory");
			return;
		}

		final TranslateResultsSettings.Builder settings = SettingsManager.readTranslateResultsSettings(0).toBuilder();

		// Show a dialog allowing the results set to be filtered
		final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
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

		final MemoryPeakResults results = ResultsManager.loadInputResults(settings.getInputOption(), false, null, null);
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
		catch (final DataException e)
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
			@Override
			public void execute(PeakResult peakResult)
			{
				// Requires a direct reference!
				final float[] params = peakResult.getParameters();
				params[PeakResult.X] += x;
				params[PeakResult.Y] += y;
				params[PeakResult.Z] += z;
			}
		});
	}
}
