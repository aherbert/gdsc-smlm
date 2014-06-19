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

import gdsc.smlm.fitting.function.Gaussian2DFunction;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.settings.FilterSettings;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.model.MaskDistribution;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;

import java.util.ArrayList;

/**
 * Filters PeakFit results that are stored in memory using various fit criteria.
 */
public class FilterResults implements PlugIn
{
	private static final String TITLE = "Filter Results";
	private static String inputOption = "";

	private MemoryPeakResults results;

	private FilterSettings filterSettings = new FilterSettings();

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
			return;
		}

		// Show a dialog allowing the results set to be filtered
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addMessage("Select a dataset to filter");
		ResultsManager.addInput(gd, inputOption, InputSource.Memory);
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		inputOption = ResultsManager.getInputSource(gd);
		results = ResultsManager.loadInputResults(inputOption, false);

		analyseResults();

		if (!showDialog())
			return;

		filterResults();
	}

	/**
	 * Analyse the results and determine the range for each filter
	 */
	private void analyseResults()
	{
		IJ.showStatus("Analysing results ...");

		int size = results.size();
		int i = 0;
		for (PeakResult result : results.getResults())
		{
			if (i % 64 == 0)
				IJ.showProgress(i, size);

			float drift = getDrift(result);
			if (maxDrift < drift)
				maxDrift = drift;
			if (minDrift > drift)
				minDrift = drift;

			float signal = result.getSignal();
			if (maxSignal < signal)
				maxSignal = signal;
			if (minSignal > signal)
				minSignal = signal;

			float snr = getSNR(signal, result);
			if (maxSNR < snr)
				maxSNR = snr;
			if (minSNR > snr)
				minSNR = snr;

			double precision = getPrecision(result, signal);
			if (maxPrecision < precision)
				maxPrecision = precision;
			if (minPrecision > precision)
				minPrecision = precision;

			float width = getWidth(result);
			averageWidth += width;
			if (maxWidth < width)
				maxWidth = width;
			if (minWidth > width)
				minWidth = width;
		}
		averageWidth /= results.size();

		IJ.showProgress(1);
		IJ.showStatus("");
	}

	private float getDrift(PeakResult result)
	{
		float drift = Math.max(Math.abs(result.origX - result.params[Gaussian2DFunction.X_POSITION]),
				Math.abs(result.origY - result.params[Gaussian2DFunction.Y_POSITION]));
		return drift;
	}

	private float getSNR(float signal, PeakResult result)
	{
		if (result.noise <= 0)
			return 0;
		return signal / result.noise;
	}

	private double getPrecision(PeakResult result, float signal)
	{
		final double s = result.getWidth() * results.getNmPerPixel();
		double precision = PeakResult.getPrecision(results.getNmPerPixel(), s, signal / results.getGain(),
				result.noise / results.getGain());
		return precision;
	}

	private float getWidth(PeakResult result)
	{
		// The X-width should be the largest (major axis)
		// Q. Should a filter be used for the Y-width too?
		return result.params[Gaussian2DFunction.X_WIDTH];
	}

	/**
	 * Check that none of the filter values are outside the limits
	 */
	private void checkLimits()
	{
		if (filterSettings.maxDrift > maxDrift || filterSettings.maxDrift < minDrift)
			filterSettings.maxDrift = maxDrift;

		if (filterSettings.minSignal > maxSignal || filterSettings.minSignal < minSignal)
			filterSettings.minSignal = minSignal;

		if (filterSettings.minSNR > maxSNR || filterSettings.minSNR < minSNR)
			filterSettings.minSNR = minSNR;

		if (filterSettings.maxPrecision > maxPrecision || filterSettings.maxPrecision < minPrecision)
			filterSettings.maxPrecision = maxPrecision;

		if (filterSettings.minWidth > maxWidth || filterSettings.minWidth < minWidth)
			filterSettings.minWidth = minWidth;

		if (filterSettings.maxWidth > maxWidth || filterSettings.maxWidth < minWidth)
			filterSettings.maxWidth = maxWidth;

		if (filterSettings.minWidth > filterSettings.maxWidth)
		{
			float tmp = filterSettings.maxWidth;
			filterSettings.maxWidth = filterSettings.minWidth;
			filterSettings.minWidth = tmp;
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
		ByteProcessor mask = getMask(filterSettings.maskTitle);
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

		int i = 0;
		int size = results.size();
		for (PeakResult result : results.getResults())
		{
			if (i % 64 == 0)
				IJ.showProgress(i, size);

			float drift = getDrift(result);
			if (drift > filterSettings.maxDrift)
				continue;

			float signal = result.getSignal();
			if (signal < filterSettings.minSignal)
				continue;

			float snr = getSNR(signal, result);
			if (snr < filterSettings.minSNR)
				continue;

			double precision = getPrecision(result, signal);
			if (precision > filterSettings.maxPrecision)
				continue;

			float width = getWidth(result);
			if (width < filterSettings.minWidth || width > filterSettings.maxWidth)
				continue;

			if (maskFilter != null)
			{
				// Check the coordinates are inside the mask
				double[] xy = new double[] { result.getXPosition() - centreX, result.getYPosition() - centreY };
				if (!maskFilter.isWithinXY(xy))
					continue;
			}

			// Passed all filters. Add to the results
			newResults.add(result);
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
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		GlobalSettings gs = SettingsManager.loadSettings();
		filterSettings = gs.getFilterSettings();

		checkLimits();

		gd.addSlider("Max_drift", minDrift, maxDrift, filterSettings.maxDrift);
		gd.addSlider("Min_Signal", minSignal, maxSignal, filterSettings.minSignal);
		gd.addSlider("Min_SNR", minSNR, maxSNR, filterSettings.minSNR);
		gd.addSlider("Min_Precision", minPrecision, maxPrecision, filterSettings.maxPrecision);

		// TODO - If calibrated present the widths in nm
		gd.addMessage("Average Width = " + IJ.d2s(averageWidth, 3));

		gd.addSlider("Min_Width", minWidth, maxWidth, filterSettings.minWidth);
		gd.addSlider("Max_Width", minWidth, maxWidth, filterSettings.maxWidth);

		// Get a list of potential mask images
		String[] items = getImageList();
		gd.addChoice("Mask", items, filterSettings.maskTitle);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		filterSettings.maxDrift = (float) gd.getNextNumber();
		filterSettings.minSignal = (float) gd.getNextNumber();
		filterSettings.minSNR = (float) gd.getNextNumber();
		filterSettings.maxPrecision = (float) gd.getNextNumber();
		filterSettings.minWidth = (float) gd.getNextNumber();
		filterSettings.maxWidth = (float) gd.getNextNumber();
		filterSettings.maskTitle = gd.getNextChoice();

		return SettingsManager.saveSettings(gs);
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
