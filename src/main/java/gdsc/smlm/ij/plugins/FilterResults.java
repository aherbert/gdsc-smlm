package gdsc.smlm.ij.plugins;

import java.util.ArrayList;

import org.apache.commons.math3.util.FastMath;

import gdsc.core.data.DataException;
import gdsc.core.utils.Maths;
import gdsc.smlm.data.config.GUIProtos.GUIFilterSettings;
import gdsc.smlm.data.config.GUIProtosHelper;
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
import gdsc.smlm.model.MaskDistribution;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.procedures.PrecisionResultProcedure;
import gdsc.smlm.results.procedures.StandardResultProcedure;
import gdsc.smlm.results.procedures.WidthResultProcedure;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.ExtendedGenericDialog;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;

/**
 * Filters PeakFit results that are stored in memory using various fit criteria.
 */
public class FilterResults implements PlugIn
{
	private static final String TITLE = "Filter Results";
	private static String inputOption = "";

	private MemoryPeakResults results;

	private GUIFilterSettings.Builder filterSettings = GUIProtosHelper.defaultGUIFilterSettings.toBuilder();

	// Used to pass data from analyseResults() to checkLimits()
	private float minDrift = Float.MAX_VALUE;
	private float maxDrift = 0;
	private float minSignal = Float.MAX_VALUE;
	private float maxSignal = 0;
	private float minSNR = Float.MAX_VALUE;
	private float maxSNR = 0;
	private double minPrecision = Float.MAX_VALUE;
	private double maxPrecision = 0;
	private double averageWidth = 0;
	private float minWidth = Float.MAX_VALUE;
	private float maxWidth = 0;

	private StandardResultProcedure sp;
	private PrecisionResultProcedure pp;
	private WidthResultProcedure wp;

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
		gd.addMessage("Select a dataset to filter");
		ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		inputOption = ResultsManager.getInputSource(gd);
		results = ResultsManager.loadInputResults(inputOption, false);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			IJ.showStatus("");
			return;
		}

		if (!analyseResults())
			return;

		if (!showDialog())
			return;

		filterResults();
	}

	/**
	 * Analyse the results and determine the range for each filter
	 */
	private boolean analyseResults()
	{
		IJ.showStatus("Analysing results ...");

		ArrayList<String> error = new ArrayList<String>();

		try
		{
			wp = new WidthResultProcedure(results, DistanceUnit.PIXEL);
			wp.getW();
			float[] limits = Maths.limits(wp.wx);
			maxWidth = limits[1];
			minWidth = limits[0];
			averageWidth = Maths.sum(wp.wx) / wp.size();
		}
		catch (DataException e)
		{
			error.add(e.getMessage());
			wp = null;
			maxWidth = minWidth = 0;
		}

		try
		{
			pp = new PrecisionResultProcedure(results);
			pp.getPrecision();

			double[] limits = Maths.limits(pp.precision);
			maxPrecision = limits[1];
			minPrecision = limits[0];
		}
		catch (DataException e)
		{
			error.add(e.getMessage());
			pp = null;
			maxPrecision = minPrecision = 0;
		}

		try
		{
			sp = new StandardResultProcedure(results, DistanceUnit.PIXEL);
			sp.getXYR();

			// Re-use for convenience
			sp.intensity = new float[sp.x.length];
			sp.background = new float[sp.x.length];
			sp.z = new float[sp.x.length];

			for (int i = 0; i < sp.size(); i++)
			{
				if (i % 64 == 0)
					IJ.showProgress(i, sp.size());

				PeakResult result = sp.peakResults[i];

				final float drift = getDrift(result, sp.x[i], sp.y[i]);
				if (maxDrift < drift)
					maxDrift = drift;
				if (minDrift > drift)
					minDrift = drift;

				final float signal = result.getSignal();
				if (maxSignal < signal)
					maxSignal = signal;
				if (minSignal > signal)
					minSignal = signal;

				final float snr = getSnr(result);
				if (maxSNR < snr)
					maxSNR = snr;
				if (minSNR > snr)
					minSNR = snr;

				// for convenience
				sp.z[i] = drift;
				sp.intensity[i] = signal;
				sp.background[i] = snr;
			}
		}
		catch (DataException e)
		{
			error.add(e.getMessage());
			sp = null;
		}

		if (error.size() == 3 || sp == null)
		{
			StringBuilder sb = new StringBuilder("Unable to analyse the results:\n");
			for (String s : error)
				sb.append(s).append(".\n");
			IJ.error(TITLE, sb.toString());
			return true;
		}

		IJ.showProgress(1);
		IJ.showStatus("");
		return false;
	}

	private float getDrift(PeakResult result, float x, float y)
	{
		return FastMath.max(Math.abs(result.origX - x), Math.abs(result.origY - y));
	}

	private float getSnr(PeakResult result)
	{
		if (result.noise <= 0)
			return 0;
		return result.getSignal() / result.noise;
	}

	/**
	 * Check that none of the filter values are outside the limits
	 */
	private void checkLimits()
	{
		if (filterSettings.getMaxDrift() > maxDrift || filterSettings.getMaxDrift() < minDrift)
			filterSettings.setMaxDrift(maxDrift);

		if (filterSettings.getMinSignal() > maxSignal || filterSettings.getMinSignal() < minSignal)
			filterSettings.setMinSignal(minSignal);

		if (filterSettings.getMinSnr() > maxSNR || filterSettings.getMinSnr() < minSNR)
			filterSettings.setMinSnr(minSNR);

		if (filterSettings.getMaxPrecision() > maxPrecision || filterSettings.getMaxPrecision() < minPrecision)
			filterSettings.setMaxPrecision(maxPrecision);

		if (filterSettings.getMinWidth() > maxWidth || filterSettings.getMinWidth() < minWidth)
			filterSettings.setMinWidth(minWidth);

		if (filterSettings.getMaxWidth() > maxWidth || filterSettings.getMaxWidth() < minWidth)
			filterSettings.setMaxWidth(maxWidth);

		if (filterSettings.getMinWidth() > filterSettings.getMaxWidth())
		{
			float tmp = filterSettings.getMaxWidth();
			filterSettings.setMaxWidth(filterSettings.getMinWidth());
			filterSettings.setMinWidth(tmp);
		}
	}

	/**
	 * Apply the filters to the data
	 */
	private void filterResults()
	{
		checkLimits();

		MemoryPeakResults newResults = new MemoryPeakResults();
		newResults.copySettings(results);
		newResults.setName(results.getName() + " Filtered");

		// Initialise the mask
		ByteProcessor mask = getMask(filterSettings.getMaskTitle());
		MaskDistribution maskFilter = null;
		final float centreX = results.getBounds().width / 2.0f;
		final float centreY = results.getBounds().height / 2.0f;
		if (mask != null)
		{
			double scaleX = (double) results.getBounds().width / mask.getWidth();
			double scaleY = (double) results.getBounds().height / mask.getHeight();
			maskFilter = new MaskDistribution((byte[]) mask.getPixels(), mask.getWidth(), mask.getHeight(), 0, scaleX,
					scaleY);
		}

		for (int i = 0, size = results.size(); i < size; i++)
		{
			if (i % 64 == 0)
				IJ.showProgress(i, size);

			// sp will not be null

			// We stored the drift=z, intensity=signal, background=snr 
			if (sp.z[i] > filterSettings.getMaxDrift())
				continue;

			if (sp.intensity[i] < filterSettings.getMinSignal())
				continue;

			if (sp.background[i] < filterSettings.getMinSnr())
				continue;

			if (maskFilter != null)
			{
				// Check the coordinates are inside the mask
				double[] xy = new double[] { sp.x[i] - centreX, sp.y[i] - centreY };
				if (!maskFilter.isWithinXY(xy))
					continue;
			}

			if (pp != null)
				if (pp.precision[i] > maxPrecision)
					continue;

			if (wp != null)
			{
				final float width = wp.wx[i];
				if (width < filterSettings.getMinWidth() || width > filterSettings.getMaxWidth())
					continue;
			}

			// Passed all filters. Add to the results
			newResults.add(sp.peakResults[i]);
		}

		IJ.showProgress(1);
		IJ.showStatus(newResults.size() + " Filtered localisations");
		MemoryPeakResults.addResults(newResults);
	}

	private ByteProcessor getMask(String maskTitle)
	{
		ImagePlus imp = WindowManager.getImage(maskTitle);
		if (imp != null)
		{
			return (ByteProcessor) imp.getProcessor().convertToByte(false);
		}
		return null;
	}

	private boolean showDialog()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		filterSettings = SettingsManager.readGUIFilterSettings(0).toBuilder();

		checkLimits();

		gd.addSlider("Max_drift", minDrift, maxDrift, filterSettings.getMaxDrift());
		gd.addSlider("Min_Signal", minSignal, maxSignal, filterSettings.getMinSignal());
		gd.addSlider("Min_SNR", minSNR, maxSNR, filterSettings.getMinSnr());
		gd.addSlider("Min_Precision", minPrecision, maxPrecision, filterSettings.getMaxPrecision());

		// TODO - If calibrated present the widths in nm
		gd.addMessage("Average Width = " + IJ.d2s(averageWidth, 3));

		gd.addSlider("Min_Width", minWidth, maxWidth, filterSettings.getMinWidth());
		gd.addSlider("Max_Width", minWidth, maxWidth, filterSettings.getMaxWidth());

		// Get a list of potential mask images
		String[] items = getImageList();
		gd.addChoice("Mask", items, filterSettings.getMaskTitle());

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		filterSettings.setMaxDrift((float) gd.getNextNumber());
		filterSettings.setMinSignal((float) gd.getNextNumber());
		filterSettings.setMinSnr((float) gd.getNextNumber());
		filterSettings.setMaxPrecision((float) gd.getNextNumber());
		filterSettings.setMinWidth((float) gd.getNextNumber());
		filterSettings.setMaxWidth((float) gd.getNextNumber());
		filterSettings.setMaskTitle(gd.getNextChoice());

		return SettingsManager.writeSettings(filterSettings.build());
	}

	/**
	 * Build a list of all the image names.
	 * 
	 * @return The list of images
	 */
	public static String[] getImageList()
	{
		ArrayList<String> newImageList = new ArrayList<String>();
		newImageList.add("[None]");

		for (int id : getIDList())
		{
			ImagePlus imp = WindowManager.getImage(id);
			if (imp == null)
				continue;
			if (!imp.getProcessor().isBinary())
				continue;
			newImageList.add(imp.getTitle());
		}

		return newImageList.toArray(new String[0]);
	}

	public static int[] getIDList()
	{
		int[] list = WindowManager.getIDList();
		return (list != null) ? list : new int[0];
	}
}
