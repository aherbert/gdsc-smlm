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

import gdsc.core.ij.IJTrackProgress;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.core.ij.Utils;
import gdsc.smlm.results.TextFilePeakResults;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.Trace;
import gdsc.smlm.results.TraceManager;
import ij.IJ;
import ij.gui.ExtendedGenericDialog;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

/**
 * Run a tracing algorithm on the peak results to trace neighbours across the frames.
 */
public class NeighbourAnalysis implements PlugIn
{
	private static final String TITLE = "Neighbour Analysis";
	private static String inputOption = "";
	private static double distanceThreshold = 0.6;
	private static int timeThreshold = 1;

	private static String filename = "";

	private MemoryPeakResults results;

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
			IJ.error(TITLE, "No localisations in memory");
			return;
		}

		if (!showDialog())
			return;

		TraceManager manager = new TraceManager(results);

		// Run the tracing
		manager.setTracker(new IJTrackProgress());
		Trace[] traces = manager.findNeighbours(distanceThreshold, timeThreshold);

		saveTraces(traces);
	}

	private void saveTraces(Trace[] traces)
	{
		String[] path = Utils.decodePath(filename);
		OpenDialog chooser = new OpenDialog("Traces_File", path[0], path[1]);
		if (chooser.getFileName() != null)
		{
			filename = chooser.getDirectory() + chooser.getFileName();

			// Remove extension and replace with .xls
			int index = filename.lastIndexOf('.');
			if (index > 0)
			{
				filename = filename.substring(0, index);
			}
			filename += ".xls";

			boolean showDeviations = (!results.getResults().isEmpty() && results.getHead().paramsStdDev != null);
			TextFilePeakResults traceResults = new TextFilePeakResults(filename, showDeviations);
			traceResults.copySettings(results);
			traceResults.begin();
			if (!traceResults.isActive())
			{
				IJ.error(TITLE, "Failed to write to file: " + filename);
				return;
			}
			traceResults.addComment(createSettingsComment());
			for (Trace trace : traces)
				traceResults.addCluster(trace); // addTrace(...) does a sort on the results
			traceResults.end();
		}
	}

	private String createSettingsComment()
	{
		return String.format("Neighbour tracing : distance-threshold = %f : time-threshold = %d", distanceThreshold,
				timeThreshold);
	}

	private boolean showDialog()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);

		gd.addNumericField("Distance_Threshold (px)", distanceThreshold, 4);
		gd.addNumericField("Time_Threshold (frames)", timeThreshold, 0);

		gd.showDialog();

		if (gd.wasCanceled() || !readDialog(gd))
			return false;

		// Load the results
		results = ResultsManager.loadInputResults(inputOption, false);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			IJ.showStatus("");
			return false;
		}

		return true;
	}

	private boolean readDialog(ExtendedGenericDialog gd)
	{
		inputOption = ResultsManager.getInputSource(gd);
		distanceThreshold = gd.getNextNumber();
		timeThreshold = (int) gd.getNextNumber();

		if (distanceThreshold < 0)
			distanceThreshold = 0;
		if (timeThreshold < 0)
			timeThreshold = 0;
		if (timeThreshold == 0 && distanceThreshold == 0)
		{
			IJ.error(TITLE, "No thresholds specified");
			return false;
		}
		return !gd.invalidNumber();
	}
}
