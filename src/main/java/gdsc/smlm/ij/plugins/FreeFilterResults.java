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
import gdsc.smlm.ij.settings.FilterSettings;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.filter.AndFilter;
import gdsc.smlm.results.filter.CoordinateFilter;
import gdsc.smlm.results.filter.Filter;
import gdsc.smlm.results.filter.OrFilter;
import gdsc.smlm.results.filter.PrecisionFilter;
import gdsc.smlm.results.filter.PrecisionHysteresisFilter;
import gdsc.smlm.results.filter.SNRFilter;
import gdsc.smlm.results.filter.SNRFilter2;
import gdsc.smlm.results.filter.SNRHysteresisFilter;
import gdsc.smlm.results.filter.SignalFilter;
import gdsc.smlm.results.filter.TraceFilter;
import gdsc.smlm.results.filter.WidthFilter;
import gdsc.smlm.results.filter.WidthFilter2;
import gdsc.smlm.utils.TextUtils;
import gdsc.smlm.utils.XmlUtils;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

import java.awt.Checkbox;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

/**
 * Filters PeakFit results that are stored in memory using the configured filters.
 */
public class FreeFilterResults implements PlugIn, ItemListener
{
	private static final String TITLE = "Free Filter Results";
	private static String inputOption = "";

	private FilterSettings filterSettings;
	private MemoryPeakResults results;

	/*
	 * (non-)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		if (MemoryPeakResults.countMemorySize() == 0)
		{
			// Ask user if they want to show the demo filters
			GenericDialog gd = new GenericDialog(TITLE);
			gd.enableYesNoCancel();
			gd.hideCancelButton();
			gd.addMessage("No results in memory. Show the demo filters?");
			gd.showDialog();
			if (gd.wasOKed())
				logDemoFilters();
			return;
		}

		if (!showDialog())
			return;

		results = ResultsManager.loadInputResults(inputOption, false);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			IJ.showStatus("");
			return;
		}

		// Filter results
		Filter filter = Filter.fromXML(filterSettings.freeFilter);
		if (filter != null)
		{
			MemoryPeakResults newResults = filter.filter(results);
			if (newResults.size() > 0)
			{
				newResults.setName(results.getName() + " Free Filtered");
				MemoryPeakResults.addResults(newResults);
			}
			IJ.showStatus(String.format("Filtered %d results to %d", results.size(), newResults.size()));
		}
		else
		{
			IJ.showStatus("ERROR: Unable to create filter");
		}
	}

	private boolean showDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Select a dataset to filter");
		ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);

		GlobalSettings gs = SettingsManager.loadSettings();
		filterSettings = gs.getFilterSettings();

		gd.addTextAreas(XmlUtils.prettyPrintXml(filterSettings.freeFilter), null, 20, 80);
		gd.addCheckbox("Show_demo_filters", false);

		if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro()))
		{
			Checkbox cb = (Checkbox) gd.getCheckboxes().get(0);
			cb.addItemListener(this);
		}

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		inputOption = ResultsManager.getInputSource(gd);
		filterSettings.freeFilter = gd.getNextText();
		boolean demoFilters = gd.getNextBoolean();

		if (demoFilters)
		{
			logDemoFilters();
			return false;
		}

		return SettingsManager.saveSettings(gs);
	}

	public void itemStateChanged(ItemEvent e)
	{
		// When the checkbox is clicked, output the list of available filters to the ImageJ log

		Checkbox cb = (Checkbox) e.getSource();
		if (cb.getState())
		{
			cb.setState(false);

			logDemoFilters();
		}
	}

	private void logDemoFilters()
	{
		comment(TITLE + " example filters");
		IJ.log("");
		comment("Filters are described using XML");
		comment("Multiple filters can be combined using AND/OR filters");
		IJ.log("");

		comment("Single filters");
		IJ.log("");
		demo(new WidthFilter(2));
		demo(new WidthFilter2(0.7f, 2));
		demo(new SignalFilter(1000));
		demo(new SNRFilter(10));
		demo(new SNRFilter2(10, 0.7f, 2));
		demo(new PrecisionFilter(30));
		demo(new SNRHysteresisFilter(2, 10, 20));
		demo(new PrecisionHysteresisFilter(2, 20, 30));
		demo(new TraceFilter(0.5, 1));
		demo(new CoordinateFilter(15.5f, 234.5f, 80.99f, 133f));

		comment("Combined filters");
		IJ.log("");
		demo(new AndFilter(new SNRFilter(10), new WidthFilter(2)));
		demo(new OrFilter(new SNRFilter(10), new PrecisionFilter(30)));
		demo(new OrFilter(new AndFilter(new SNRFilter(10), new PrecisionFilter(30)), new TraceFilter(0.5, 1)));
	}

	private void demo(Filter filter)
	{
		comment(filter.getClass().getSimpleName() + ": " + filter.getDescription());
		IJ.log(filter.toXML());
		IJ.log("");
	}

	private void comment(String text)
	{
		IJ.log(TextUtils.wrap("<!-- " + text + " -->", 80));
	}
}
