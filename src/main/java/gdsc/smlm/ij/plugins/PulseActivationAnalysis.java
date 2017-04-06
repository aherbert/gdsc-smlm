package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.Rectangle;
import java.util.ArrayList;

import gdsc.core.ij.Utils;
import gdsc.core.utils.Settings;
import gdsc.core.utils.TurboList;
import gdsc.smlm.fitting.linear.EJMLLinearSolver;

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
 * Perform multi-channel super-resolution imaging by means of photo-switchable probes and pulsed light activation.
 * 
 * This plugin is based on the methods described in: Mark Bates, Bo Huang, Graham T. Dempsey, Xiaowei Zhuang (2007).
 * Multicolor Super-Resolution Imaging with Photo-Switchable Fluorescent Probes. Science 317, 1749. DOI:
 * 10.1126/science.1146598.
 */
public class PulseActivationAnalysis implements PlugIn, DialogListener
{
	private String TITLE = " Activation Analysis";

	private static String inputOption = "";
	private static int channels = 1;
	private static final int MAX_CHANNELS = 3;

	private static int pulseInterval = 10;
	private static int afterPulseStart = 1;
	private static int afterPulseEnd = 1;
	private static int darkFramesForNewActivation = 1;
	private static double pulseActivationPercentage = 50;

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
		double percentage;

		PulseActivationResult(Trace trace)
		{
			this.trace = trace;
		}

		void setCounts(int activations, int pulseActivation)
		{
			this.activations = activations;
			this.pulseActivation = pulseActivation;
			percentage = (100.0 * pulseActivation) / activations;
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

		boolean crosstalkMode = "crosstalk".equals(arg);

		if (!showDialog(crosstalkMode))
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

		if (crosstalkMode)
			runCrosstalkAnalysis();
		else
			runPulseAnalysis();
	}

	private boolean showDialog(boolean crosstalkMode)
	{
		TITLE = ((crosstalkMode) ? "Crosstalk" : "Pulse") + TITLE;

		GenericDialog gd = new GenericDialog(TITLE);

		if (crosstalkMode)
			gd.addMessage("Analyse crosstalk activation rate");
		else
			gd.addMessage("Count & plot molecules activated after a pulse");

		ResultsManager.addInput(gd, "Input", inputOption, InputSource.MEMORY_CLUSTERED);
		if (!crosstalkMode)
			gd.addSlider("Channels", 1, MAX_CHANNELS, channels);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		inputOption = ResultsManager.getInputSource(gd);
		if (!crosstalkMode)
		{
			channels = (int) gd.getNextNumber();

			if (channels < 1 || channels > MAX_CHANNELS)
			{
				IJ.error(TITLE, "Channels must be between 1 and " + MAX_CHANNELS);
				return false;
			}
		}

		return true;
	}

	private void runCrosstalkAnalysis()
	{
		// TODO 

		// Determine the cross talk ratio.
		// This is done by imaging only a single photo-switchable probe with the 
		// same activation pulse imaging routine used for multi-colour imaging.
		// Concept:
		// A probe is meant to turn on in a frame following a pulse from a specific wavelength.
		// Multi-wavelengths can be used with probes responding to each wavelength. However
		// each probe may be activated bby the 'wrong' wavelength. This is crosstalk.
		// The idea is to understand how many times the probe will turn on in a 
		// frame following a pulse from the other lasers.

		// To determine the crosstalk ratio we must have a single probe imaged with the full
		// multi-wavelength pulse cycle. We then count how many times a probe activated by 
		// the correct wavelength is activated by the others.
		// This is done by grouping localisation into clusters to represent the same molecule.
		// Crosstalk for each wavelength is then the fraction of times molecules were activated
		// by the 'wrong' wavelength.

		// Require:
		// N-channels
		// The full pulse activation cycle
		// The target wavelength (i.e. the wavelength that is not 'wrong')
		// The clustering distance to group activations into molecules
		// The number of activation to define a molecule

	}

	/**
	 * Unmix the observed local densities into the actual densities for 2-channels.
	 * <p>
	 * Crosstalk from M into N is defined as the number of times the molecule that should be activated by a pulse from
	 * channel N is activated by a pulse from channel M. A value less than 1 is expected (otherwise the fluorophore is
	 * being non-specifically activated by is target channel, N).
	 *
	 * @param D1
	 *            the observed density in channel 1
	 * @param D2
	 *            the observed density in channel 2
	 * @param C21
	 *            the crosstalk from channel 2 into channel 1
	 * @param C12
	 *            the crosstalk from channel 1 into channel 2
	 * @return the actual densities [d1, d2]
	 */
	public static double[] unmix(double D1, double D2, double C21, double C12)
	{
		double d1, d2;
		// Solve the equations:
		// D1 = d1 + C21 * d2
		// D2 = d2 + C12 * d1
		// This is done by direct substitution
		d1 = (D1 - C21 * D2) / (1 - C12 * C21);
		d2 = D2 - C12 * d1;
		return new double[] { d1, d2 };
	}

	/**
	 * Unmix the observed local densities into the actual densities for 3-channels.
	 * <p>
	 * Crosstalk from M into N is defined as the number of times the molecule that should be activated by a pulse from
	 * channel N is activated by a pulse from channel M. A value less than 1 is expected (otherwise the fluorophore is
	 * being non-specifically activated by is target channel, N).
	 *
	 * @param D1
	 *            the observed density in channel 1
	 * @param D2
	 *            the observed density in channel 2
	 * @param D3
	 *            the observed density in channel 3
	 * @param C21
	 *            the crosstalk from channel 2 into channel 1
	 * @param C31
	 *            the crosstalk from channel 3 into channel 1
	 * @param C12
	 *            the crosstalk from channel 1 into channel 2
	 * @param C32
	 *            the crosstalk from channel 3 into channel 2
	 * @param C13
	 *            the crosstalk from channel 1 into channel 3
	 * @param C23
	 *            the crosstalk from channel 2 into channel 3
	 * @return the actual densities [d1, d2, d3]
	 */
	public static double[] unmix(double D1, double D2, double D3, double C21, double C31, double C12, double C32,
			double C13, double C23)
	{
		// Solve the linear equations
		// D1 = d1 + C21 * d2 + C31 * d3
		// D2 = d2 + C12 * d1 + C32 * d3
		// D3 = d3 + C13 * d1 + C23 * d2
		// This is done using matrix decomposition

		EJMLLinearSolver solver = new EJMLLinearSolver();
		// @formatter:off
		double[][] a = { 
				{ 1, C21, C31 },
				{ C12, 1, C32 },
				{ C13, C23, 1 }
		};
		// @formatter:on
		double[] b = { D1, D2, D3 };
		if (!solver.solveLinear(a, b))
		{
			// Unmix failed so reset
			b[0] = D1;
			b[1] = D2;
			b[2] = D3;
		}

		return b;
	}

	private boolean runPulseAnalysis()
	{
		GenericDialog gd = new GenericDialog(TITLE);

		gd.addMessage("Count & plot molecules activated after a pulse.\nDataset = " + inputOption +
				"\nPulse interval = " + pulseInterval);
		gd.addNumericField("After_pulse_start", afterPulseStart, 0);
		gd.addNumericField("After_pulse_end", afterPulseEnd, 0);
		gd.addNumericField("Dark_frames_for_new_activation", darkFramesForNewActivation, 0);
		gd.addSlider("Pulse_activation_percentage", 0, 100, pulseActivationPercentage);

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
		pulseActivationPercentage = gd.getNextNumber();

		// Check arguments
		try
		{
			Parameters.isAboveZero("After pulse start", afterPulseStart);
			Parameters.isEqualOrBelow("After pulse start", afterPulseStart, pulseInterval);
			Parameters.isEqualOrAbove("After pulse end", afterPulseEnd, afterPulseStart);
			Parameters.isEqualOrBelow("After pulse end", afterPulseEnd, pulseInterval);
			Parameters.isEqualOrAbove("Pulse activation percentage", pulseActivationPercentage, 0);
			Parameters.isEqualOrBelow("Pulse activation percentage", pulseActivationPercentage, 100);
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
			if (result.percentage >= pulseActivationPercentage)
			{
				count++;
				output.addAll(activations.get(i).trace.getPoints());
			}
		}

		output.end();

		IJ.showStatus(String.format("%d/%s, %d/%s", count, Utils.pleural(traces.length, "Trace"), output.size(),
				Utils.pleural(results.size(), "Result")));
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
