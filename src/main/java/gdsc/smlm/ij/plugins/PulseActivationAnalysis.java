package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.Rectangle;
import java.util.ArrayList;

import gdsc.core.ij.Utils;
import gdsc.core.utils.Settings;
import gdsc.core.utils.TurboList;

/*----------------------------------------------------------------------------- 
 * GDSC Plugins for ImageJ
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.results.IJImagePeakResults;
import gdsc.smlm.ij.results.ImagePeakResultsFactory;
import gdsc.smlm.ij.results.ResultsImage;
import gdsc.smlm.ij.results.ResultsMode;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.ResultsSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResultsList;
import gdsc.smlm.results.Trace;
import gdsc.smlm.results.TraceManager;
import ij.IJ;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

/**
 * Count the number of times a molecules is localised immediately after an activation pulse compared to any other
 * activation.
 */
public class PulseActivationAnalysis implements PlugIn, DialogListener
{
	private static final String TITLE = "Pulse Activation Analysis";

	private static String inputOption = "";
	private static int pulseInterval = 10;
	private static int afterPulseStart = 1;
	private static int afterPulseEnd = 1;
	private static int darkFramesForNewActivation = 1;
	private static double pulseActivationFraction = 0.5;

	private GlobalSettings settings;
	private ResultsSettings resultsSettings;
	private MemoryPeakResults results;
	private Trace[] traces;

	private class PulseActivationResult
	{
		Trace trace;
		@SuppressWarnings("unused")
		int activations;
		@SuppressWarnings("unused")
		int pulseActivation;
		double fraction;

		PulseActivationResult(Trace trace)
		{
			this.trace = trace;
		}

		void setCounts(int activations, int pulseActivation)
		{
			this.activations = activations;
			this.pulseActivation = pulseActivation;
			fraction = (double) pulseActivation / activations;
		}
	}

	private TurboList<PulseActivationResult> activations;

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

		// Load the results
		results = ResultsManager.loadInputResults(inputOption, false);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			return;
		}

		// Get the traces
		traces = TraceManager.convert(results);
		if (traces == null || traces.length == 0)
		{
			IJ.error(TITLE, "No traces could be loaded");
			return;
		}

		// Initialisation
		activations = new TurboList<PulseActivationResult>(traces.length);
		for (Trace trace : traces)
		{
			trace.sort(); // Time-order
			this.activations.add(new PulseActivationResult(trace));
		}

		showAnalysisDialog();
	}

	private boolean showDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);

		gd.addMessage("Count & plot molecules activated after a pulse");
		ResultsManager.addInput(gd, "Input", inputOption, InputSource.MEMORY_CLUSTERED);
		gd.addNumericField("Pulse_interval", pulseInterval, 0);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		inputOption = ResultsManager.getInputSource(gd);
		pulseInterval = (int) Math.abs(gd.getNextNumber());

		if (pulseInterval < 2)
		{
			IJ.error(TITLE, "Pulse interval must be above 1");
			return false;
		}

		return true;
	}

	private boolean showAnalysisDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);

		gd.addMessage("Count & plot molecules activated after a pulse.\nDataset = " + inputOption +
				"\nPulse Interval = " + pulseInterval);
		gd.addNumericField("After_pulse_start", afterPulseStart, 0);
		gd.addNumericField("After_pulse_end", afterPulseEnd, 0);
		gd.addNumericField("Dark_frames_for_new_activation", darkFramesForNewActivation, 0);
		gd.addSlider("Pulse_activation_fraction", 0.05, 1, pulseActivationFraction);

		settings = SettingsManager.loadSettings();
		resultsSettings = settings.getResultsSettings();

		gd.addMessage("--- Image output ---");
		String[] imageNames = SettingsManager.getNames((Object[]) ResultsImage.values());
		gd.addChoice("Image", imageNames, imageNames[resultsSettings.getResultsImage().ordinal()]);
		gd.addCheckbox("Weighted", resultsSettings.weightedImage);
		gd.addCheckbox("Equalised", resultsSettings.equalisedImage);
		gd.addSlider("Image_Precision (nm)", 5, 30, resultsSettings.precision);
		gd.addSlider("Image_Scale", 1, 15, resultsSettings.imageScale);

		gd.addCheckbox("Preview", false);
		gd.addDialogListener(this);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		SettingsManager.saveSettings(settings);

		return true;
	}

	/** The changed flag indicating that work has yet to be done. */
	private boolean changed = true;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.gui.DialogListener#dialogItemChanged(ij.gui.GenericDialog, java.awt.AWTEvent)
	 */
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
	{
		// The event is not null when the GenericDialog's components are updated
		if (e != null)
			changed = true;

		afterPulseStart = (int) gd.getNextNumber();
		afterPulseEnd = (int) gd.getNextNumber();
		darkFramesForNewActivation = Math.max(1, (int) gd.getNextNumber());
		pulseActivationFraction = gd.getNextNumber();

		// Check arguments
		try
		{
			Parameters.isAboveZero("After pulse start", afterPulseStart);
			Parameters.isEqualOrBelow("After pulse start", afterPulseStart, pulseInterval);
			Parameters.isEqualOrAbove("After pulse end", afterPulseEnd, afterPulseStart);
			Parameters.isEqualOrBelow("After pulse end", afterPulseEnd, pulseInterval);
			Parameters.isEqualOrAbove("Pulse activation ratio", pulseActivationFraction, 0);
			Parameters.isEqualOrBelow("Pulse activation ratio", pulseActivationFraction, 1);
		}
		catch (IllegalArgumentException ex)
		{
			IJ.error(TITLE, ex.getMessage());
			return false;
		}

		resultsSettings.setResultsImage(gd.getNextChoiceIndex());
		resultsSettings.weightedImage = gd.getNextBoolean();
		resultsSettings.equalisedImage = gd.getNextBoolean();
		resultsSettings.precision = gd.getNextNumber();
		resultsSettings.imageScale = gd.getNextNumber();
		boolean preview = gd.getNextBoolean();

		if (gd.invalidNumber())
			return false;

		// Only run if there are changes.
		// The event is null when the GenericDialog is OK'd. 
		// In that case we run to process the changes. 
		if (changed && (preview || e == null))
			run();

		return true;
	}

	private void run()
	{
		changed = false;

		getActivationResults();

		// Set-up outputs
		PeakResultsList output = new PeakResultsList();
		output.copySettings(results);
		output.setName(results.getName() + " " + TITLE);

		// Store the set in memory
		MemoryPeakResults memoryResults = new MemoryPeakResults(this.results.size());
		MemoryPeakResults.addResults(memoryResults);
		output.addOutput(memoryResults);

		// Draw the super-resolution image
		Rectangle bounds = results.getBounds(true);
		addImageResults(output, results.getName(), bounds, results.getNmPerPixel(), results.getGain());

		output.begin();

		// Create a results set with only those molecules above the ratio
		int count = 0;
		for (int i = activations.size(); i-- > 0;)
		{
			PulseActivationResult result = activations.getf(i);
			if (result.fraction >= pulseActivationFraction)
			{
				count++;
				output.addAll(activations.get(i).trace.getPoints());
			}
		}

		output.end();

		IJ.showStatus(String.format("%d/%s, %d/%s", count, Utils.pleural(traces.length, "Trace"),
				output.size(), Utils.pleural(results.size(), "Result")));
	}

	private Settings lastSettings;

	private synchronized void getActivationResults()
	{
		// Check if any settings have changed. If not then there is no need to re-compute.
		// Note: pulseInterval is not collected in the second dialog 
		Settings settings = new Settings(afterPulseStart, afterPulseEnd, darkFramesForNewActivation);
		if (settings.equals(lastSettings))
			return;
		lastSettings = settings;

		// Activations are only counted if there are at least 
		// n frames between localisations.
		final int n = darkFramesForNewActivation + 1;

		// For each trace
		for (int i = activations.size(); i-- > 0;)
		{
			PulseActivationResult result = activations.getf(i);
			// Count the number of activations
			ArrayList<PeakResult> points = result.trace.getPoints();
			int nextPulseStart = Integer.MIN_VALUE;
			int activationCount = 0;
			int pulseActivationCount = 0;
			for (int j = 0; j < points.size(); j++)
			{
				PeakResult p = points.get(j);
				// Check if this is an activation
				if (p.getFrame() >= nextPulseStart)
				{
					activationCount++;
					// Classify if within a pulse window
					int mod = p.getFrame() % pulseInterval;
					if (mod >= afterPulseStart && mod <= afterPulseEnd)
					{
						pulseActivationCount++;
					}
				}
				nextPulseStart = p.getEndFrame() + n;
			}
			result.setCounts(activationCount, pulseActivationCount);
		}
	}

	private void addImageResults(PeakResultsList resultsList, String title, Rectangle bounds, double nmPerPixel,
			double gain)
	{
		if (resultsSettings.getResultsImage() != ResultsImage.NONE)
		{
			IJImagePeakResults image = ImagePeakResultsFactory.createPeakResultsImage(resultsSettings.getResultsImage(),
					resultsSettings.weightedImage, resultsSettings.equalisedImage, title, bounds, nmPerPixel, gain,
					resultsSettings.imageScale, resultsSettings.precision, ResultsMode.ADD);
			image.setLiveImage(false);
			resultsList.addOutput(image);
		}
	}
}
