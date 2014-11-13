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

import java.awt.Rectangle;

import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import ij.IJ;
import ij.plugin.PlugIn;
import ij.text.TextWindow;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 * Produces a summary table of the results that are stored in memory.
 */
public class SummariseResults implements PlugIn
{
	private static final String TITLE = "Summarise Results";

	private static TextWindow summary = null;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		if (MemoryPeakResults.countMemorySize() == 0)
		{
			IJ.error(TITLE, "There are no fitting results in memory");
			clearSummaryTable();
			return;
		}

		createSummaryTable();
		for (MemoryPeakResults result : MemoryPeakResults.getAllResults())
		{
			addSummary(result);
		}
		summary.append("");
		summary.toFront();
	}

	/**
	 * Remove all entries in the summary table if it showing
	 */
	public static void clearSummaryTable()
	{
		if (summary != null && summary.isShowing())
			summary.getTextPanel().clear();
	}

	private void createSummaryTable()
	{
		if (summary == null || !summary.isShowing())
		{
			StringBuilder sb = new StringBuilder("Dataset\tN\tFrames\tTime\tMemory\tBounds\tnm/pixel\tGain\tms/frame");
			for (String statName : new String[] { "Precision (nm)", "SNR" })
			{
				sb.append("\tAv ").append(statName);
				sb.append("\tMedian ").append(statName);
				sb.append("\tMin ").append(statName);
				sb.append("\tMax ").append(statName);
			}
			summary = new TextWindow("Peak Results Summary", sb.toString(), "", 800, 300);
			summary.setVisible(true);
		}
		// This could be optional but at current there is no dialog and it seems unnecessary
		clearSummaryTable();
	}

	private void addSummary(MemoryPeakResults result)
	{
		DescriptiveStatistics[] stats = new DescriptiveStatistics[2];
		for (int i = 0; i < stats.length; i++)
		{
			stats[i] = new DescriptiveStatistics();
		}

		// Only process the statistics if we have a noise component
		final int size = result.size();
		if (size > 0 && result.getResults().get(0).noise > 0)
		{

			int ii = 0;
			final double nmPerPixel = result.getNmPerPixel();
			final double gain = result.getGain();
			final boolean emCCD = result.isEMCCD();
			for (PeakResult peakResult : result.getResults())
			{
				if (peakResult == null)
				{
					System.out.printf("Null result in summary @ %d\n", ++ii);
					continue;
				}
				stats[0].addValue(peakResult.getPrecision(nmPerPixel, gain, emCCD));
				stats[1].addValue(peakResult.getSignal() / peakResult.noise);
			}
		}

		StringBuilder sb = new StringBuilder();
		sb.append(result.getName());
		sb.append("\t").append(result.size());
		int maxT = getMaxT(result);
		sb.append("\t").append(maxT);
		final double exposureTime = (result.getCalibration() != null) ? result.getCalibration().exposureTime : 0;
		sb.append("\t").append(Utils.timeToString(maxT * exposureTime));
		if (size > 0)
		{
			boolean includeDeviations = result.getResults().get(0).paramsStdDev != null;
			long memorySize = MemoryPeakResults.estimateMemorySize(size, includeDeviations);
			String memory = MemoryPeakResults.memorySizeString(memorySize);
			sb.append("\t").append(memory);
		}
		else
		{
			sb.append("\t-");
		}
		Rectangle bounds = result.getBounds(true);
		sb.append(String.format("\t%d,%d,%d,%d\t%s\t%s\t%s", bounds.x, bounds.y, bounds.x + bounds.width, bounds.y +
				bounds.height, Utils.rounded(result.getNmPerPixel(), 4), Utils.rounded(result.getGain(), 4),
				Utils.rounded(exposureTime, 4)));
		for (int i = 0; i < stats.length; i++)
		{
			if (Double.isNaN(stats[i].getMean()))
			{
				sb.append("\t-\t-\t-\t-");
			}
			else
			{
				sb.append("\t").append(IJ.d2s(stats[i].getMean(), 3));
				sb.append("\t").append(IJ.d2s(stats[i].getPercentile(50), 3));
				sb.append("\t").append(IJ.d2s(stats[i].getMin(), 3));
				sb.append("\t").append(IJ.d2s(stats[i].getMax(), 3));
			}
		}

		summary.append(sb.toString());
	}

	private int getMaxT(MemoryPeakResults result)
	{
		int maxT = 0;
		for (PeakResult r : result.getResults())
			if (maxT < r.getEndFrame())
				maxT = r.getEndFrame();
		return maxT;
	}
}
