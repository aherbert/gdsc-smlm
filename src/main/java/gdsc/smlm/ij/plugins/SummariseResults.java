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

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import gdsc.core.data.DataException;
import gdsc.core.ij.Utils;
import gdsc.smlm.data.config.CalibrationProtosHelper;
import gdsc.smlm.data.config.CalibrationReader;
import gdsc.smlm.data.config.FitProtos.PrecisionMethod;
import gdsc.smlm.data.config.FitProtosHelper;
import gdsc.smlm.data.config.UnitHelper;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.procedures.PeakResultProcedure;
import gdsc.smlm.results.procedures.PrecisionResultProcedure;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.text.TextWindow;

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
			sb.append("\tnm/pixel\tms/frame\tCamera\tDUnit\tIUnit\tPrecision Method");
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
		final DescriptiveStatistics[] stats = new DescriptiveStatistics[2];
		for (int i = 0; i < stats.length; i++)
		{
			stats[i] = new DescriptiveStatistics();
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

		CalibrationReader calibration = result.getCalibrationReader();

		PrecisionMethod precisionMethod = PrecisionMethod.PRECISION_METHOD_NA;
		boolean stored = false;
		final int size = result.size();
		if (size > 0)
		{
			// Precision 
			try
			{
				PrecisionResultProcedure pp = new PrecisionResultProcedure(result);
				// Use stored precision if possible
				stored = result.hasPrecision();
				precisionMethod = pp.getPrecision(stored);
				for (double v : pp.precision)
					stats[0].addValue(v);
			}
			catch (DataException e)
			{
				// Ignore
			}

			// SNR requires noise
			if (result.hasNoise())
			{
				result.forEach(new PeakResultProcedure()
				{
					public void execute(PeakResult peakResult)
					{
						stats[1].addValue(peakResult.getSignal() / peakResult.getNoise());
					}
				});
			}
		}

		sb.append(result.getName());
		int maxT = 0;
		if (result.size() == 0)
		{
			sb.append("\t0\t0");
		}
		else
		{
			sb.append('\t').append(result.size());
			maxT = result.getMaxFrame();
			sb.append('\t').append(maxT);
		}
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
			//@formatter:off
			sb.append('\t').append(calibration.hasNmPerPixel() ? Utils.rounded(calibration.getNmPerPixel()) : '-');
			sb.append('\t').append(calibration.hasExposureTime() ? Utils.rounded(calibration.getExposureTime()) : '-');
			
			if (calibration.hasCameraType())
			{
				sb.append('\t').append(CalibrationProtosHelper.getName(calibration.getCameraType()));
				if (calibration.isCCDCamera())
				{
					sb.append(" bias=").append(calibration.getBias());
					sb.append(" gain=").append(calibration.getCountPerPhoton());
				}
			}
			else
			{
				sb.append("\t-");
			}
			
			sb.append('\t').append(calibration.hasDistanceUnit() ? UnitHelper.getShortName(calibration.getDistanceUnit()) : '-');
			sb.append('\t').append(calibration.hasIntensityUnit() ? UnitHelper.getShortName(calibration.getIntensityUnit()) : '-');
			//@formatter:on
		}
		else
		{
			sb.append("\t\t\t\t\t");
		}
		sb.append("\t").append(FitProtosHelper.getName(precisionMethod));
		if (stored)
			sb.append(" (Stored)");
		for (int i = 0; i < stats.length; i++)
		{
			if (Double.isNaN(stats[i].getMean()))
			{
				sb.append("\t-\t-\t-\t-");
			}
			else
			{
				sb.append('\t').append(IJ.d2s(stats[i].getMean(), 3));
				sb.append('\t').append(IJ.d2s(stats[i].getPercentile(50), 3));
				sb.append('\t').append(IJ.d2s(stats[i].getMin(), 3));
				sb.append('\t').append(IJ.d2s(stats[i].getMax(), 3));
			}
		}
		sb.append("\n");
	}
}
