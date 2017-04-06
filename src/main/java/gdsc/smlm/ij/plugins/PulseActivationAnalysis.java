package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;

import gdsc.core.ij.Utils;
import gdsc.core.utils.Maths;
import gdsc.core.utils.TextUtils;
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
import gdsc.smlm.results.Cluster.CentroidMethod;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResults;
import gdsc.smlm.results.PeakResultsList;
import gdsc.smlm.results.Trace;
import gdsc.smlm.results.TraceManager;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.Plot2;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

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

	/**
	 * Specify the method to use to determine the parameters for the distribution of the localisation precision (assumed
	 * to be Gaussian)
	 */
	private enum CrosstalkCorrection
	{
		//@formatter:off
		NONE{ public String getName() { return "None"; }},
		SUBTRACTION{ public String getName() { return "Subtraction"; }},
		SWITCH{ public String getName() { return "Switch"; }};
		//@formatter:on

		@Override
		public String toString()
		{
			return getName();
		}

		/**
		 * Gets the name.
		 *
		 * @return the name
		 */
		abstract public String getName();
	}

	private static String inputOption = "";
	private static int channels = 1;
	private static final int MAX_CHANNELS = 3;

	private static int repeatInterval = 30;
	private static int[] startFrame = { 1, 11, 21 };
	// Crosstalk 
	private static double[] ct = new double[6];
	private static String[] ctNames = { "21", "31", "12", "32", "13", "23" };
	private static final int C21 = 0;
	private static final int C31 = 1;
	private static final int C12 = 2;
	private static final int C32 = 3;
	private static final int C13 = 4;
	private static final int C23 = 5;
	private static int darkFramesForNewActivation = 1;

	private static int targetChannel = 1;

	private static double densityRadius = 35;
	private static int crosstalkCorrection = CrosstalkCorrection.SUBTRACTION.ordinal();
	private CrosstalkCorrection correction = CrosstalkCorrection.NONE;
	private static double[] subtractionCutoff = { 50, 50, 50 };
	private static boolean nonspecificAssignment = false;
	private static double nonspecificAssignmentCutoff = 50;

	private GlobalSettings settings;
	private ResultsSettings resultsSettings;
	private MemoryPeakResults results;
	private Trace[] traces;

	private class Activation
	{
		final Trace trace;
		final float x, y;
		final int channel;
		int currentChannel;

		Activation(Trace trace, int channel)
		{
			this.trace = trace;
			float[] centroid = trace.getCentroid(CentroidMethod.SIGNAL_WEIGHTED);
			x = centroid[0];
			y = centroid[1];
			this.channel = channel;
			currentChannel = channel;
		}
	}

	private TurboList<Activation> activations;

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

		if (!showPulseCycleDialog())
			return;

		createActivations();

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
		int min = (crosstalkMode) ? 2 : 1;
		gd.addSlider("Channels", min, MAX_CHANNELS, channels);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		inputOption = ResultsManager.getInputSource(gd);
		channels = (int) gd.getNextNumber();
		if (channels < min || channels > MAX_CHANNELS)
		{
			IJ.error(TITLE, "Channels must be between " + min + " and " + MAX_CHANNELS);
			return false;
		}

		return true;
	}

	private boolean showPulseCycleDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);

		gd.addMessage("Specify the pulse cycle");

		gd.addNumericField("Repeat_interval", repeatInterval, 0);
		gd.addNumericField("Dark_frames_for_new_activation", darkFramesForNewActivation, 0);
		for (int c = 1; c <= channels; c++)
			gd.addNumericField("Activation_frame_C" + c, startFrame[c - 1], 0);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		repeatInterval = (int) gd.getNextNumber();
		if (repeatInterval < channels)
		{
			IJ.error(TITLE, "Repeat interval must be greater than the number of channels: " + channels);
			return false;
		}
		darkFramesForNewActivation = Math.max(1, (int) gd.getNextNumber());
		for (int c = 1; c <= channels; c++)
		{
			int frame = (int) gd.getNextNumber();
			if (frame < 1 || frame > repeatInterval)
			{
				IJ.error(TITLE, "Channel " + c + " activation frame must within the repeat interval");
				return false;
			}
			startFrame[c - 1] = frame;
		}

		// Check all start frames are unique
		for (int i = 0; i < channels; i++)
			for (int j = i + 1; j < channels; j++)
				if (startFrame[i] == startFrame[j])
				{
					IJ.error(TITLE, "Start frames must be unique for each channel");
					return false;
				}

		return true;
	}

	/**
	 * Creates the activations. This splits the input traces into continuous chains of localisations. Each chain is an
	 * activation. A new activation is created if there are more than the configured number of dark frames since the
	 * last localisation. The start frame for the activation defines the channel the activation is assigned to (this may
	 * be channel 0 if the start frame is not in a pulse start frame).
	 */
	private void createActivations()
	{
		activations = new TurboList<Activation>(traces.length);

		// Activations are only counted if there are at least 
		// n frames between localisations.
		final int n = darkFramesForNewActivation + 1;

		for (Trace trace : traces)
		{
			trace.sort(); // Time-order			

			ArrayList<PeakResult> points = trace.getPoints();

			// Define the frame for a new activation 
			int nextActivationStartFrame = Integer.MIN_VALUE;
			Trace current = null;
			int channel = 0;
			for (int j = 0; j < points.size(); j++)
			{
				PeakResult p = points.get(j);
				// Check if this is an activation
				if (p.getFrame() >= nextActivationStartFrame)
				{
					if (current != null)
						// Store the last
						activations.add(new Activation(trace, channel));

					// Create a new activation
					current = new Trace(p);
					channel = getChannel(p);
				}
				else
				{
					// This is the same chain of localisations
					current.add(p);
				}
				nextActivationStartFrame = p.getEndFrame() + n;
			}

			if (current != null)
				activations.add(new Activation(trace, channel));
		}
	}

	private int getChannel(PeakResult p)
	{
		// Classify if within a channel activation start frame
		final int mod = p.getFrame() % repeatInterval;
		for (int i = 0; i < channels; i++)
			if (mod == startFrame[i])
				return i + 1;
		return 0;
	}

	private void runCrosstalkAnalysis()
	{
		// Determine the cross talk ratio.
		// This is done by imaging only a single photo-switchable probe with the 
		// same activation pulse imaging routine used for multi-colour imaging.
		// Concept:
		// A probe is meant to turn on in a frame following a pulse from a specific wavelength.
		// Multi-wavelengths can be used with probes responding to each wavelength. However
		// each probe may be activated by the 'wrong' wavelength. This is crosstalk.
		// The idea is to understand how many times the probe will turn on in a 
		// frame following a pulse from the other lasers.

		// To determine the crosstalk ratio we must have a single probe imaged with the full
		// multi-wavelength pulse cycle. We then count how many times a probe activated by 
		// the correct wavelength is activated by the others.
		// Crosstalk for each wavelength is then the fraction of times molecules were activated
		// by the 'wrong' wavelength.

		if (!showCrossTalkAnalysisDialog())
			return;

		// Count the activations per channel
		int[] count = new int[channels + 1]; // for convenience use 1-index channel
		for (int i = activations.size(); i-- > 0;)
		{
			Activation result = activations.getf(i);
			if (result.channel != 0)
				count[result.channel]++;
		}

		double[] crosstalk = new double[count.length];
		long sum = Maths.sum(count);
		for (int c = 1; c <= channels; c++)
			crosstalk[c] = (double) count[c] / sum;

		// Store the cross talk
		int index1, index2 = -1;
		if (channels == 2)
		{
			if (targetChannel == 1)
				index1 = setCrosstalk(C21, crosstalk[2]);
			else
				index1 = setCrosstalk(C12, crosstalk[1]);
		}
		else
		{
			// 3-channel
			if (targetChannel == 1)
			{
				index1 = setCrosstalk(C21, crosstalk[2]);
				index2 = setCrosstalk(C31, crosstalk[3]);
			}
			else if (targetChannel == 2)
			{
				index1 = setCrosstalk(C12, crosstalk[1]);
				index2 = setCrosstalk(C32, crosstalk[3]);
			}
			else
			{
				index1 = setCrosstalk(C13, crosstalk[1]);
				index2 = setCrosstalk(C23, crosstalk[2]);
			}
		}

		// Plot a histogram
		double[] x = Utils.newArray(channels, 0.5, 1);
		double[] y = Arrays.copyOfRange(crosstalk, 1, crosstalk.length);
		Plot2 plot = new Plot2(TITLE, "Channel", "Fraction activations");
		plot.setLimits(0, channels + 1, 0, Maths.max(y) * 1.05);
		plot.setXMinorTicks(false);
		plot.addPoints(x, y, Plot2.BAR);
		String label = String.format("Crosstalk %s = %s", ctNames[index1], Maths.round(ct[index1]));
		if (index2 > -1)
			label += String.format(", %s = %s", ctNames[index2], Maths.round(ct[index2]));
		plot.addLabel(0, 0, label);
		Utils.display(TITLE, plot);
	}

	private int setCrosstalk(int index, double value)
	{
		ct[index] = value;
		return index;
	}

	private boolean showCrossTalkAnalysisDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);

		gd.addMessage(TextUtils.wrap(
				"Crosstalk analysis requires a sample singly labelled with only one photo-switable probe and imaged with the full pulse lifecycle.",
				80));

		String[] ch = new String[channels];
		for (int i = 0; i < ch.length; i++)
			ch[i] = "Channel " + (i + 1);

		gd.addChoice("Target", ch, "Channel " + targetChannel);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		targetChannel = gd.getNextChoiceIndex() + 1;

		return true;
	}

	/**
	 * Unmix the observed local densities into the actual densities for 2-channels.
	 * <p>
	 * Crosstalk from M into N is defined as the number of times the molecule that should be activated by a pulse from
	 * channel N is activated by a pulse from channel M. A value less than 0.5 is expected (otherwise the fluorophore is
	 * not being specifically activated by channel N).
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
		// Solve the equations:
		// D1 = d1 + C21 * d2
		// D2 = d2 + C12 * d1
		// This is done by direct substitution
		double d1 = (D1 - C21 * D2) / (1 - C12 * C21);
		double d2 = D2 - C12 * d1;
		// Assuming D1 and D2 are positive and C12 and C21 are 
		// between 0 and 0.5 then we do not need to check the bounds.
		//d1 = Maths.clip(0, D1, d1);
		//d2 = Maths.clip(0, D2, d2);
		return new double[] { d1, d2 };
	}

	/**
	 * Unmix the observed local densities into the actual densities for 3-channels.
	 * <p>
	 * Crosstalk from M into N is defined as the number of times the molecule that should be activated by a pulse from
	 * channel N is activated by a pulse from channel M. A total value for crosstalk into N is expected to be less than
	 * 2/3 (otherwise the fluorophore is not being non-specifically activated by channel N).
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
				{   1, C21, C31 },
				{ C12,   1, C32 },
				{ C13, C23,   1 }
		};
		// @formatter:on
		double[] b = { D1, D2, D3 };
		if (!solver.solveLinear(a, b))
		{
			// Unmix failed so reset to the observed densities
			b[0] = D1;
			b[1] = D2;
			b[2] = D3;
		}
		else
		{
			// Due to floating-point error in the decomposition we check the bounds
			b[0] = Maths.clip(0, D1, b[0]);
			b[1] = Maths.clip(0, D2, b[1]);
			b[2] = Maths.clip(0, D3, b[2]);
		}

		return b;
	}

	private boolean runPulseAnalysis()
	{
		GenericDialog gd = new GenericDialog(TITLE);

		gd.addMessage("Plot molecules activated after a pulse");
		if (channels > 1)
		{
			if (channels == 2)
			{
				gd.addNumericField("Crosstalk_21", ct[C21], 3);
				gd.addNumericField("Crosstalk_12", ct[C12], 3);
			}
			else
			{
				gd.addNumericField("Crosstalk_21", ct[C21], 3);
				gd.addNumericField("Crosstalk_31", ct[C31], 3);
				gd.addNumericField("Crosstalk_12", ct[C12], 3);
				gd.addNumericField("Crosstalk_32", ct[C32], 3);
				gd.addNumericField("Crosstalk_13", ct[C13], 3);
				gd.addNumericField("Crosstalk_23", ct[C23], 3);
			}

			gd.addNumericField("Local_density_radius", densityRadius, 0, 6, "nm");
			String[] correctionNames = SettingsManager.getNames((Object[]) CrosstalkCorrection.values());
			gd.addChoice("Crosstalk_correction", correctionNames, correctionNames[crosstalkCorrection]);
			for (int c = 1; c <= channels; c++)
				gd.addSlider("Subtraction_cutoff_C" + c + "(%)", 0, 100, subtractionCutoff[c - 1]);
			gd.addCheckbox("Nonspecific_assignment", nonspecificAssignment);
			gd.addSlider("Nonspecific_assignment_cutoff (%)", 0, 100, nonspecificAssignmentCutoff);
		}

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

		// Check arguments
		try
		{
			if (channels > 1)
			{
				if (channels == 2)
				{
					ct[C21] = gd.getNextNumber();
					ct[C12] = gd.getNextNumber();
					validateCrosstalk(C21);
					validateCrosstalk(C12);
				}
				else
				{
					ct[C21] = gd.getNextNumber();
					ct[C31] = gd.getNextNumber();
					ct[C12] = gd.getNextNumber();
					ct[C32] = gd.getNextNumber();
					ct[C13] = gd.getNextNumber();
					ct[C23] = gd.getNextNumber();
					for (int i = 0; i < ct.length; i += 2)
						validateCrosstalk(i, i + 1);
				}

				densityRadius = Math.abs(gd.getNextNumber());
				crosstalkCorrection = gd.getNextChoiceIndex();
				if (crosstalkCorrection >= 0 && crosstalkCorrection < CrosstalkCorrection.values().length)
					correction = CrosstalkCorrection.values()[crosstalkCorrection];
				for (int c = 1; c <= channels; c++)
				{
					subtractionCutoff[c - 1] = (int) gd.getNextNumber();
					validatePercentage("Subtraction_cutoff_C" + c, subtractionCutoff[c - 1]);
				}
				nonspecificAssignment = gd.getNextBoolean();
				nonspecificAssignmentCutoff = gd.getNextNumber();
				validatePercentage("Nonspecific_assignment_cutoff", nonspecificAssignmentCutoff);
			}
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

	private void validateCrosstalk(int index)
	{
		String name = "Crosstalk " + ctNames[index];
		Parameters.isPositive(name, ct[index]);
		Parameters.isBelow(name, ct[index], 0.5);
	}

	private void validateCrosstalk(int index1, int index2)
	{
		validateCrosstalk(index1);
		validateCrosstalk(index2);
		Parameters.isBelow("Crosstalk " + ctNames[index1] + " + " + ctNames[index2], ct[index1] + ct[index2], 0.5);
	}

	private void validatePercentage(String name, double d)
	{
		Parameters.isPositive(name, d);
		Parameters.isEqualOrBelow(name, d, 100);
	}

	private void run()
	{
		changed = false;

		// Assign all activations to a channel.
		// This is only necessary when we have more than 1 channel. If we have 1 channel then 
		// no correction method is specified.
		if (correction != CrosstalkCorrection.NONE)
		{
			// TODO - Build a density manager variant that can put all the activations on a grid
			// It has a method to count the number of activations within a radius from each channel
			// using the currentChannel property.

			// Reset
			for (int i = activations.size(); i-- > 0;)
			{
				Activation result = activations.getf(i);
				result.currentChannel = result.channel;
			}

			int[] newChannel = new int[activations.size()];
			for (int i = activations.size(); i-- > 0;)
			{
				Activation result = activations.getf(i);

				if (result.channel == 0)
					// Ignore non-specific activations
					continue;

				// Find the observed local densities in each channel

				// Compute the true local densities

				// Apply crosstalk correction
				if (correction == CrosstalkCorrection.SUBTRACTION)
				{
					// Compute the probability it is correct

					// Remove it if below the subtraction threshold
				}
				else
				{
					// Switch
					// Compute the probability of each channel and randomly select

				}
			}

			// Update the channel assignment
			for (int i = activations.size(); i-- > 0;)
			{
				activations.getf(i).currentChannel = newChannel[i];
			}

			// Assign non-specific activations
			if (nonspecificAssignment)
			{
				for (int i = activations.size(); i-- > 0;)
				{
					Activation result = activations.getf(i);

					if (result.channel != 0)
						continue;

					// Find the observed local densities in each channel

					// Take the observed as the true local densities (since correction has been applied)

					// Compute the probability of each channel and randomly select

				}
			}
		}

		// Set-up outputs for each channel
		PeakResultsList[] output = new PeakResultsList[channels];
		for (int c = 0; c < channels; c++)
			output[c] = createOutput(c + 1);

		// Create a results set with only those molecules assigned to a channel
		int count = 0;
		for (int i = activations.size(); i-- > 0;)
		{
			Activation result = activations.getf(i);
			if (result.currentChannel == 0)
				continue;
			output[result.currentChannel - 1].addAll(result.trace.getPoints());
		}

		for (int c = 0; c < channels; c++)
			output[c].end();

		// Collate image into a stack
		if (channels > 1 && resultsSettings.getResultsImage() != ResultsImage.NONE)
		{
			ImageStack stack = null; // We do not yet know the size
			for (int c = 1; c <= channels; c++)
			{
				ImageProcessor ip = getImage(output[c - 1]);
				if (stack == null)
					stack = new ImageStack(ip.getWidth(), ip.getHeight());
				stack.addSlice("C" + c, ip);
			}
			String name = results.getName() + " " + TITLE;
			ImagePlus imp = Utils.display(name, stack);
			imp.setDimensions(channels, 1, 1);
		}

		IJ.showStatus(String.format("%d/%s, %d/%s", count, Utils.pleural(traces.length, "Trace"), output[0].size(),
				Utils.pleural(results.size(), "Result")));
	}

	private PeakResultsList createOutput(int c)
	{
		PeakResultsList output = new PeakResultsList();
		output.copySettings(results);
		if (channels > 1)
			output.setName(results.getName() + " " + TITLE + " C" + c);
		else
			output.setName(results.getName() + " " + TITLE);

		// Store the set in memory
		MemoryPeakResults memoryResults = new MemoryPeakResults(this.results.size());
		output.addOutput(memoryResults);
		MemoryPeakResults.addResults(memoryResults);

		// Draw the super-resolution image
		Rectangle bounds = results.getBounds(true);
		addImageResults(output, results.getName(), bounds, results.getNmPerPixel(), results.getGain());

		output.begin();

		return output;
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
			image.setDisplayImage(channels == 1);
			resultsList.addOutput(image);
		}
	}

	private ImageProcessor getImage(PeakResultsList peakResultsList)
	{
		PeakResults[] list = peakResultsList.toArray();
		IJImagePeakResults image = (IJImagePeakResults) list[1];
		return image.getImagePlus().getProcessor();
	}
}
