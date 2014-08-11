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
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.Trace;
import gdsc.smlm.results.TraceManager;
import gdsc.smlm.utils.Maths;
import gdsc.smlm.utils.StoredDataStatistics;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.PlugIn;

import java.util.Arrays;

/**
 * Computes a graph of the dark time and estimates the time threshold for the specified point in the
 * cumulative histogram.
 */
public class DarkTimeAnalysis implements PlugIn
{
	private static String TITLE = "Dark-time Analysis";

	private static String inputOption = "";
	private static double msPerFrame = 50;
	private static double searchDistance = 100;
	private static double maxDarkTime = 0;
	private static double percentile = 99;
	private static int nBins = 100;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		// Require some fit results and selected regions
		if (MemoryPeakResults.countMemorySize() == 0)
		{
			IJ.error(TITLE, "There are no fitting results in memory");
			return;
		}

		if (!showDialog())
			return;

		MemoryPeakResults results = ResultsManager.loadInputResults(inputOption, true);
		msPerFrame = results.getCalibration().exposureTime;
		Utils.log("%s: %d localisations", TITLE, results.size());

		if (results.size() == 0)
		{
			IJ.error(TITLE, "No results were loaded");
			return;
		}

		analyse(results);
	}

	private boolean showDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Compute the cumulative dark-time histogram");
		ResultsManager.addInput(gd, inputOption, InputSource.Memory);

		gd.addSlider("Search_distance (nm)", 5, 150, searchDistance);
		gd.addNumericField("Max_dark_time (seconds)", maxDarkTime, 2);
		gd.addSlider("Percentile", 0, 100, percentile);
		gd.addSlider("Histogram_bins", 0, 100, nBins);
		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		inputOption = gd.getNextChoice();
		searchDistance = gd.getNextNumber();
		maxDarkTime = gd.getNextNumber();
		percentile = gd.getNextNumber();
		nBins = (int) Math.abs(gd.getNextNumber());

		// Check arguments
		try
		{
			Parameters.isAboveZero("Search distance", searchDistance);
			Parameters.isPositive("Percentile", percentile);
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		return true;
	}

	private void analyse(MemoryPeakResults results)
	{
		// Find min and max time frames
		results.sort();
		int min = results.getResults().get(0).peak;
		int max = results.getResults().get(results.size() - 1).getEndFrame();

		// Trace results
		TraceManager tm = new TraceManager(results);
		int range = max - min + 1;
		if (maxDarkTime > 0 && results.getCalibration() != null && results.getCalibration().exposureTime > 0)
			range = Math.max(1, (int) Math.round(maxDarkTime * 1000 / results.getCalibration().exposureTime));
		tm.traceMolecules(searchDistance / results.getCalibration().nmPerPixel, range);
		Trace[] traces = tm.getTraces();

		// Build dark-time histogram
		int[] times = new int[range];
		StoredDataStatistics stats = new StoredDataStatistics();
		for (Trace trace : traces)
		{
			if (trace.getNBlinks() > 1)
			{
				for (int t : trace.getOffTimes())
				{
					times[t]++;
				}
				stats.add(trace.getOffTimes());
			}
		}

		plotDarkTimeHistogram(stats);

		// Cumulative histogram
		for (int i = 1; i < times.length; i++)
			times[i] += times[i - 1];
		int total = times[times.length - 1];

		// Plot dark-time up to 100%
		double[] x = new double[range];
		double[] y = new double[range];
		int truncate = 0;
		for (int i = 0; i < x.length; i++)
		{
			x[i] = i * msPerFrame;
			y[i] = (100.0 * times[i]) / total;
			if (times[i] == total) // 100%
			{
				truncate = i + 1;
				break;
			}
		}
		if (truncate > 0)
		{
			x = Arrays.copyOf(x, truncate);
			y = Arrays.copyOf(y, truncate);
		}

		String title = "Cumulative Dark-time";
		Plot plot = new Plot(title, "Time (ms)", "Percentile", x, y);
		Utils.display(title, plot);

		// Report percentile
		for (int i = 0; i < y.length; i++)
		{
			if (y[i] >= percentile)
			{
				Utils.log("Dark-time Percentile %.1f @ %s ms = %s s", percentile, Utils.rounded(x[i]),
						Utils.rounded(x[i] / 1000));
				break;
			}
		}
	}

	private void plotDarkTimeHistogram(StoredDataStatistics stats)
	{
		if (nBins > 0)
		{
			String title = "Dark-time";

			// Ensure the bin width is never less than 1
			float yMax = (int) Math.ceil(Maths.max(stats.getFloatValues()));
			int newBins = (int) (yMax + 1);
			float[][] hist = Utils.calcHistogram(stats.getFloatValues(), 0, yMax, (nBins > newBins) ? newBins : nBins);

			// Create the axes
			float[] xValues = Utils.createHistogramAxis(hist[0]);
			float[] yValues = Utils.createHistogramValues(hist[1]);

			// Convert the X-axis to milliseconds
			for (int i = 0; i < xValues.length; i++)
				xValues[i] *= msPerFrame;

			// Plot
			Plot plot = new Plot("Dark-time", "Time (ms)", "Frequency", xValues, yValues);
			if (xValues.length > 0)
			{
				double xPadding = 0.05 * (xValues[xValues.length - 1] - xValues[0]);
				plot.setLimits(xValues[0] - xPadding, xValues[xValues.length - 1] + xPadding, 0,
						Maths.max(yValues) * 1.05);
			}
			Utils.display(title, plot);
		}
	}
}
