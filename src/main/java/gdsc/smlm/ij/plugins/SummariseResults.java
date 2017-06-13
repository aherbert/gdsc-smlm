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

import gdsc.core.ij.Utils;
import gdsc.smlm.results.Calibration;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import ij.IJ;
import ij.gui.GenericDialog;
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
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (MemoryPeakResults.isMemoryEmpty())
		{
			IJ.error(TITLE, "There are no fitting results in memory");
			clearSummaryTable();
			return;
		}

		createSummaryTable();
		StringBuilder sb = new StringBuilder();
		int i = 0;
		int nextFlush = 9;
		for (MemoryPeakResults result : MemoryPeakResults.getAllResults())
		{
			addSummary(sb, result);
			if (++i == nextFlush)
			{
				summary.append(sb.toString());
				sb.setLength(0);
			}
		}
		summary.append(sb.toString());
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
			StringBuilder sb = new StringBuilder("Dataset\tN\tFrames\tTime\tMemory\tBounds");
			// Calibration
			sb.append("\tnm/pixel\tGain\tms/frame\tCamera\tDUnit\tIUnit");
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

	private static int NO = -1, UNKNOWN = 0, YES = 1;
	private int removeNullResults = UNKNOWN;

	private void addSummary(StringBuilder sb, MemoryPeakResults result)
	{
		DescriptiveStatistics[] stats = new DescriptiveStatistics[2];
		char[] suffix = new char[2];
		for (int i = 0; i < stats.length; i++)
		{
			stats[i] = new DescriptiveStatistics();
			suffix[i] = 0;
		}

		if (result.hasNullResults())
		{
			IJ.log("Null results in dataset: " + result.getName());
			if (removeNullResults == UNKNOWN)
			{
				GenericDialog gd = new GenericDialog(TITLE);
				gd.addMessage("There are invalid results in memory.\n \nClean these results?");
				gd.enableYesNoCancel();
				gd.hideCancelButton();
				gd.showDialog();
				removeNullResults = (gd.wasOKed()) ? YES : NO;
			}
			if (removeNullResults == NO)
				result = result.copy();
			result.removeNullResults();
		}

		final int size = result.size();
		if (size > 0)
		{
			// Precision
			// Check if we can use the stored precision
			if (result.hasStoredPrecision())
			{
				suffix[0] = '*';
				for (PeakResult peakResult : result.getResults())
				{
					stats[0].addValue(peakResult.getPrecision());
				}
			}
			else if (result.isCalibratedForPrecision())
			{
				final double nmPerPixel = result.getNmPerPixel();
				final double gain = result.getGain();
				final boolean emCCD = result.isEMCCD();
				for (PeakResult peakResult : result.getResults())
				{
					stats[0].addValue(peakResult.getPrecision(nmPerPixel, gain, emCCD));
				}
			}

			// SNR requires noise
			if (result.hasNoise())
			{
				for (PeakResult peakResult : result.getResults())
				{
					stats[1].addValue(peakResult.getSignal() / peakResult.noise);
				}
			}
		}

		Calibration calibration = result.getCalibration();

		sb.append(result.getName());
		sb.append('\t').append(result.size());
		int maxT = getMaxT(result);
		sb.append('\t').append(maxT);
		if (calibration != null && calibration.hasExposureTime())
		{
			sb.append('\t').append(Utils.timeToString(maxT * calibration.getExposureTime()));
		}
		else
		{
			sb.append("\t-");
		}
		if (size > 0)
		{
			boolean includeDeviations = result.hasDeviations();
			long memorySize = MemoryPeakResults.estimateMemorySize(size, includeDeviations);
			String memory = MemoryPeakResults.memorySizeString(memorySize);
			sb.append('\t').append(memory);
		}
		else
		{
			sb.append("\t-");
		}
		Rectangle bounds = result.getBounds(true);
		sb.append(
				String.format("\t%d,%d,%d,%d", bounds.x, bounds.y, bounds.x + bounds.width, bounds.y + bounds.height));
		if (calibration != null)
		{
			sb.append('\t').append(calibration.hasNmPerPixel() ? Utils.rounded(calibration.getNmPerPixel()) : '-');
			sb.append('\t').append(calibration.hasGain() ? Utils.rounded(calibration.getGain()) : '-');
			sb.append('\t').append(calibration.hasExposureTime() ? Utils.rounded(calibration.getExposureTime()) : '-');
			sb.append('\t').append(calibration.hasCameraType() ? calibration.getCameraType() : '-');
			sb.append('\t').append(calibration.hasDistanceUnit() ? calibration.getDistanceUnit() : '-');
			sb.append('\t').append(calibration.hasIntensityUnit() ? calibration.getIntensityUnit() : '-');
		}
		else
		{
			sb.append("\t\t\t\t\t\t");
		}
		for (int i = 0; i < stats.length; i++)
		{
			if (Double.isNaN(stats[i].getMean()))
			{
				sb.append("\t-\t-\t-\t-");
			}
			else
			{
				sb.append('\t').append(IJ.d2s(stats[i].getMean(), 3));
				if (suffix[i] != 0)
					sb.append(suffix[i]);
				sb.append('\t').append(IJ.d2s(stats[i].getPercentile(50), 3));
				sb.append('\t').append(IJ.d2s(stats[i].getMin(), 3));
				sb.append('\t').append(IJ.d2s(stats[i].getMax(), 3));
			}
		}
		sb.append("\n");
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
