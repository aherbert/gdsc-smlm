package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.LinearSolver;
import org.ejml.factory.LinearSolverFactory;

import gdsc.core.clustering.DensityCounter;
import gdsc.core.clustering.DensityCounter.Molecule;
import gdsc.core.ij.Utils;
import gdsc.core.utils.Maths;
import gdsc.core.utils.TextUtils;
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
import gdsc.smlm.results.Cluster.CentroidMethod;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResults;
import gdsc.smlm.results.PeakResultsList;
import gdsc.smlm.results.Trace;
import gdsc.smlm.results.TraceManager;
import ij.CompositeImage;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Plot2;
import ij.plugin.PlugIn;
import ij.plugin.frame.Recorder;
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

	/**
	 * Specify the method to use to determine the parameters for the distribution of the localisation precision (assumed
	 * to be Gaussian)
	 */
	private enum NonSpecificAssignment
	{
		//@formatter:off
		NONE{ public String getName() { return "None"; }},
		MOST_LIKELY{ public String getName() { return "Most likely"; }};
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
	private static int crosstalkCorrectionIndex = CrosstalkCorrection.SUBTRACTION.ordinal();
	private static double[] subtractionCutoff = { 50, 50, 50 };
	private static int nonSpecificAssignmentIndex = NonSpecificAssignment.NONE.ordinal();
	private static double nonSpecificAssignmentCutoff = 50;

	private GlobalSettings settings;
	private ResultsSettings resultsSettings;
	private MemoryPeakResults results;
	private Trace[] traces;

	private class Activation implements Molecule
	{
		final Trace trace;
		float x, y;
		final int channelx;
		int currentChannel;

		Activation(Trace trace, int channel)
		{
			this.trace = trace;
			float[] centroid = trace.getCentroid(CentroidMethod.SIGNAL_WEIGHTED);
			x = centroid[0];
			y = centroid[1];
			// zero index
			this.channelx = channel - 1;
			currentChannel = channel;
		}

		boolean hasSpecificChannel()
		{
			return channelx != -1;
		}

		int getChannel()
		{
			return channelx;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see gdsc.core.clustering.DensityCounter.Molecule#getX()
		 */
		public float getX()
		{
			return x;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see gdsc.core.clustering.DensityCounter.Molecule#getY()
		 */
		public float getY()
		{
			return y;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see gdsc.core.clustering.DensityCounter.Molecule#getID()
		 */
		public int getID()
		{
			// Allow the ID to be updated from the original channel by using a current channel field
			return currentChannel;
		}
	}

	private Activation[] specificActivations, nonSpecificActivations;
	private int[] count;

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
		TurboList<Activation> activations = new TurboList<Activation>(traces.length);

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

		save(activations);
	}

	private void save(TurboList<Activation> list)
	{
		// Count the activations per channel
		// Note: Channels are 0-indexed in the activations
		count = new int[channels];
		for (int i = list.size(); i-- > 0;)
		{
			Activation result = list.getf(i);
			if (result.hasSpecificChannel())
				count[result.getChannel()]++;
		}

		// Store specific activations
		int sum = (int) Maths.sum(count);
		specificActivations = new Activation[sum];
		nonSpecificActivations = new Activation[list.size() - sum];
		for (int i = list.size(), c1 = 0, c2 = 0; i-- > 0;)
		{
			Activation result = list.getf(i);
			if (result.hasSpecificChannel())
				specificActivations[c1++] = result;
			else
				nonSpecificActivations[c2++] = result;
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

		double[] crosstalk = new double[count.length];
		long sum = Maths.sum(count);
		for (int c = 0; c < channels; c++)
			crosstalk[c] = (double) count[c] / sum;

		// Store the cross talk
		int index1, index2 = -1;
		if (channels == 2)
		{
			if (targetChannel == 1)
				index1 = setCrosstalk(C21, crosstalk[1]);
			else
				index1 = setCrosstalk(C12, crosstalk[0]);
		}
		else
		{
			// 3-channel
			if (targetChannel == 1)
			{
				index1 = setCrosstalk(C21, crosstalk[1]);
				index2 = setCrosstalk(C31, crosstalk[2]);
			}
			else if (targetChannel == 2)
			{
				index1 = setCrosstalk(C12, crosstalk[0]);
				index2 = setCrosstalk(C32, crosstalk[2]);
			}
			else
			{
				index1 = setCrosstalk(C13, crosstalk[0]);
				index2 = setCrosstalk(C23, crosstalk[1]);
			}
		}

		// Plot a histogram
		double[] x = Utils.newArray(channels, 0.5, 1);
		double[] y = crosstalk;
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
		// Solve the linear equations: A * X = B
		// D1 = d1 + C21 * d2 + C31 * d3
		// D2 = d2 + C12 * d1 + C32 * d3
		// D3 = d3 + C13 * d1 + C23 * d2
		// This is done using matrix decomposition

		// Note: The linear solver uses LU decomposition. 
		// We cannot use a faster method as the matrix A is not symmetric.
		final LinearSolver<DenseMatrix64F> linearSolver = LinearSolverFactory.linear(3);
		final DenseMatrix64F A = new DenseMatrix64F(3, 3);
		A.set(0, 1);
		A.set(1, C21);
		A.set(2, C31);
		A.set(3, C12);
		A.set(4, 1);
		A.set(5, C32);
		A.set(6, C13);
		A.set(7, C23);
		A.set(8, 1);
		final DenseMatrix64F B = new DenseMatrix64F(3, 1);
		B.set(0, D1);
		B.set(1, D2);
		B.set(2, D3);

		if (!linearSolver.setA(A))
		{
			// Failed so reset to the observed densities
			return new double[] { D1, D2, D3 };
		}

		// Input B is not modified to so we can re-use for output X
		linearSolver.solve(B, B);

		final double[] b = B.getData();

		// Due to floating-point error in the decomposition we check the bounds
		b[0] = Maths.clip(0, D1, b[0]);
		b[1] = Maths.clip(0, D2, b[1]);
		b[2] = Maths.clip(0, D3, b[2]);

		return b;
	}

	private boolean runPulseAnalysis()
	{
		NonBlockingGenericDialog gd = new NonBlockingGenericDialog(TITLE);

		gd.addMessage("Plot molecules activated after a pulse");
		String[] correctionNames = null;
		String[] assigmentNames = null;

		if (channels > 1)
		{
			if (channels == 2)
			{
				gd.addNumericField("Crosstalk_21", ct[C21], 3);
				gd.addNumericField("Crosstalk_12", ct[C12], 3);
			}
			else
			{
				for (int i = 0; i < ctNames.length; i++)
					gd.addNumericField("Crosstalk_" + ctNames[i], ct[i], 3);
			}

			gd.addNumericField("Local_density_radius", densityRadius, 0, 6, "nm");
			correctionNames = SettingsManager.getNames((Object[]) CrosstalkCorrection.values());
			gd.addChoice("Crosstalk_correction", correctionNames, correctionNames[crosstalkCorrectionIndex]);
			for (int c = 1; c <= channels; c++)
				gd.addSlider("Subtraction_cutoff_C" + c + "(%)", 0, 100, subtractionCutoff[c - 1]);
			assigmentNames = SettingsManager.getNames((Object[]) NonSpecificAssignment.values());
			gd.addChoice("Nonspecific_assigment", assigmentNames, assigmentNames[nonSpecificAssignmentIndex]);
			gd.addSlider("Nonspecific_assignment_cutoff (%)", 0, 100, nonSpecificAssignmentCutoff);
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
		{
			if (executor != null)
				executor.shutdown();
			return false;
		}

		run();

		if (executor != null)
			executor.shutdown();

		// Record options for a macro since the NonBlockingDialog does not
		if (Recorder.record)
		{
			if (channels > 1)
			{
				if (channels == 2)
				{
					Recorder.recordOption("Crosstalk_21", Double.toString(ct[C21]));
					Recorder.recordOption("Crosstalk_12", Double.toString(ct[C12]));
				}
				else
				{
					for (int i = 0; i < ctNames.length; i++)
						Recorder.recordOption("Crosstalk_" + ctNames[i], Double.toString(ct[i]));
				}

				Recorder.recordOption("Local_density_radius", Double.toString(densityRadius));

				Recorder.recordOption("Crosstalk_correction", correctionNames[crosstalkCorrectionIndex]);
				for (int c = 1; c <= channels; c++)
					Recorder.recordOption("Subtraction_cutoff_C" + c, Double.toString(subtractionCutoff[c - 1]));

				Recorder.recordOption("Nonspecific_assigment", assigmentNames[nonSpecificAssignmentIndex]);
				Recorder.recordOption("Nonspecific_assignment_cutoff (%)",
						Double.toString(nonSpecificAssignmentCutoff));
			}

			Recorder.recordOption("Image", imageNames[resultsSettings.getResultsImage().ordinal()]);
			if (resultsSettings.weightedImage)
				Recorder.recordOption("Weighted");
			if (resultsSettings.equalisedImage)
				Recorder.recordOption("Equalised");
			Recorder.recordOption("Image_Precision", Double.toString(resultsSettings.precision));
			Recorder.recordOption("Image_Scale", Double.toString(resultsSettings.imageScale));
		}

		SettingsManager.saveSettings(settings);

		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.gui.DialogListener#dialogItemChanged(ij.gui.GenericDialog, java.awt.AWTEvent)
	 */
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
	{
		// The event is null when the NonBlockingGenericDialog is first shown
		if (e == null)
		{
			// Do not ignore this if a macro
			if (Utils.isMacro())
				return true;
		}

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
				crosstalkCorrectionIndex = gd.getNextChoiceIndex();
				for (int c = 1; c <= channels; c++)
				{
					subtractionCutoff[c - 1] = (int) gd.getNextNumber();
					validatePercentage("Subtraction_cutoff_C" + c, subtractionCutoff[c - 1]);
				}
				nonSpecificAssignmentIndex = gd.getNextChoiceIndex();
				nonSpecificAssignmentCutoff = gd.getNextNumber();
				validatePercentage("Nonspecific_assignment_cutoff", nonSpecificAssignmentCutoff);
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

		runSettings = new RunSettings();
		if (preview)
		{
			run();
		}

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

	private class RunSettings
	{
		double densityRadius;
		CrosstalkCorrection crosstalkCorrection = CrosstalkCorrection.NONE;
		double[] subtractionCutoff;
		NonSpecificAssignment nonSpecificAssignment = NonSpecificAssignment.NONE;
		double nonSpecificAssignmentCutoff;

		@SuppressWarnings("unused")
		ResultsSettings resultsSettings;

		RunSettings()
		{
			// Copy the current settings required
			if (channels > 1)
			{
				this.densityRadius = PulseActivationAnalysis.densityRadius;
				int crosstalkCorrectionIndex = PulseActivationAnalysis.crosstalkCorrectionIndex;
				if (crosstalkCorrectionIndex >= 0 && crosstalkCorrectionIndex < CrosstalkCorrection.values().length)
					crosstalkCorrection = CrosstalkCorrection.values()[crosstalkCorrectionIndex];
				this.subtractionCutoff = new double[channels];
				for (int i = channels; i-- > 0;)
					// Convert from percentage to a probability
					this.subtractionCutoff[i] = PulseActivationAnalysis.subtractionCutoff[i] / 100.0;
				int nonSpecificAssignmentIndex = PulseActivationAnalysis.nonSpecificAssignmentIndex;
				if (nonSpecificAssignmentIndex >= 0 &&
						nonSpecificAssignmentIndex < NonSpecificAssignment.values().length)
					nonSpecificAssignment = NonSpecificAssignment.values()[nonSpecificAssignmentIndex];
				this.nonSpecificAssignmentCutoff = PulseActivationAnalysis.nonSpecificAssignmentCutoff / 100.0;
			}
			this.resultsSettings = PulseActivationAnalysis.this.resultsSettings.clone();
		}

		public boolean newUnmixSettings(RunSettings lastRunSettings)
		{
			if (lastRunSettings == null)
				return true;
			if (lastRunSettings.densityRadius != densityRadius)
				return true;
			if (lastRunSettings.crosstalkCorrection != crosstalkCorrection)
				return true;
			if (crosstalkCorrection == CrosstalkCorrection.SUBTRACTION)
			{
				for (int i = channels; i-- > 0;)
					if (lastRunSettings.subtractionCutoff[i] != subtractionCutoff[i])
						return true;
			}
			return false;
		}

		public boolean newNonSpecificAssignmentSettings(RunSettings lastRunSettings)
		{
			if (lastRunSettings == null)
				return true;
			if (lastRunSettings.densityRadius != densityRadius)
				return true;
			if (lastRunSettings.nonSpecificAssignment != nonSpecificAssignment)
				return true;
			if (nonSpecificAssignment == NonSpecificAssignment.MOST_LIKELY)
			{
				if (lastRunSettings.nonSpecificAssignmentCutoff != nonSpecificAssignmentCutoff)
					return true;
			}
			return false;
		}
	}

	private DensityCounter dc = null;
	private int[][] density = null;
	private int nThreads;
	private ExecutorService executor = null;
	private TurboList<Future<?>> futures = null;
	private RunSettings lastRunSettings = null, runSettings = null;

	private synchronized void run()
	{
		// This is synchronized since it updates the class results. 
		// Note: We check against the last settings and only repeat what is necessary ...

		if (runSettings == null)
		{
			lastRunSettings = null;
			return;
		}

		// Assign all activations to a channel.
		// This is only necessary when we have more than 1 channel. If we have 1 channel then 
		// no correction method is specified.
		if (runSettings.newUnmixSettings(lastRunSettings))
		{
			// Reset
			for (int i = specificActivations.length; i-- > 0;)
			{
				Activation result = specificActivations[i];
				result.currentChannel = result.channelx;
			}

			if (runSettings.crosstalkCorrection != CrosstalkCorrection.NONE)
			{
				// Use a density counter that can put all the activations on a grid.
				// It has a method to count the number of activations within a radius that 
				// belong to each channel.

				// Add only those with specific activations. Non-specific activations are ignored.
				createDensityCounter((float) runSettings.densityRadius);

				// Do this all together: it uses a faster algorithm and we can cache the results
				if (density == null)
				{
					IJ.showStatus("Computing observed density");
					density = dc.countAll(channels);
				}

				long seed = System.currentTimeMillis();

				// -=-=-=--=-=-
				// Unmix the specific activations to their correct channel.
				// -=-=-=--=-=-
				IJ.showStatus("Unmixing");
				createThreadPool();

				int[] newChannel = new int[specificActivations.length];

				int nPerThread = (int) Math.ceil((double) specificActivations.length / nThreads);
				for (int from = 0; from < specificActivations.length;)
				{
					int to = Math.min(from + nPerThread, specificActivations.length);
					futures.add(executor.submit(new UnmixWorker(runSettings, density, newChannel, from, to, seed++)));
					from = to;
				}
				waitToFinish();

				// Update the channel assignment
				for (int i = specificActivations.length; i-- > 0;)
				{
					specificActivations[i].currentChannel = newChannel[i];
				}
			}
		}

		// -=-=-=--=-=-
		// Assign non-specific activations
		// -=-=-=--=-=-
		if (runSettings.newNonSpecificAssignmentSettings(lastRunSettings))
		{
			// Reset
			for (int i = nonSpecificActivations.length; i-- > 0;)
			{
				Activation result = nonSpecificActivations[i];
				result.currentChannel = result.channelx;
			}

			if (runSettings.nonSpecificAssignment != NonSpecificAssignment.NONE)
			{
				createDensityCounter((float) runSettings.densityRadius);

				long seed = System.currentTimeMillis();

				IJ.showStatus("Non-specific assignment");
				createThreadPool();

				int[] newChannel = new int[nonSpecificActivations.length];

				int nPerThread = (int) Math.ceil((double) nonSpecificActivations.length / nThreads);
				for (int from = 0; from < nonSpecificActivations.length;)
				{
					int to = Math.min(from + nPerThread, nonSpecificActivations.length);
					futures.add(
							executor.submit(new NonSpecificUnmixWorker(runSettings, dc, newChannel, from, to, seed++)));
					from = to;
				}
				waitToFinish();

				// Update the channel assignment
				for (int i = nonSpecificActivations.length; i-- > 0;)
				{
					nonSpecificActivations[i].currentChannel = newChannel[i];
				}
			}
		}

		// Set-up outputs for each channel
		IJ.showStatus("Creating outputs");
		PeakResultsList[] output = new PeakResultsList[channels];
		for (int c = 0; c < channels; c++)
			output[c] = createOutput(c + 1);

		// Create a results set with only those molecules assigned to a channel
		int count = write(output, specificActivations, 0);
		count = write(output, specificActivations, count);

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
				ip.setColorModel(null);
				stack.addSlice("C" + c, ip);
			}

			// Create a composite
			String name = results.getName() + " " + TITLE;
			ImagePlus imp = new ImagePlus(name, stack);
			imp.setDimensions(channels, 1, 1);
			CompositeImage ci = new CompositeImage(imp, IJ.COMPOSITE);

			imp = WindowManager.getImage(name);
			if (imp != null && imp.isComposite())
			{
				imp.setImage(ci);
			}
			else
			{
				ci.show();
			}
		}

		lastRunSettings = runSettings;
		runSettings = null;

		IJ.showStatus(String.format("%d/%s, %d/%s", count, Utils.pleural(traces.length, "Trace"), output[0].size(),
				Utils.pleural(results.size(), "Result")));
	}

	private void createThreadPool()
	{
		if (executor == null)
		{
			nThreads = Prefs.getThreads();
			executor = Executors.newFixedThreadPool(nThreads);
			futures = new TurboList<Future<?>>(nThreads);
		}
	}

	private void createDensityCounter(float densityRadius)
	{
		if (dc == null || dc.radius != densityRadius)
		{
			dc = new DensityCounter(specificActivations, densityRadius, false);
			// Clear cache of density
			density = null;
		}
	}

	private void waitToFinish()
	{
		// Wait for all to finish
		for (int t = futures.size(); t-- > 0;)
		{
			try
			{
				// The future .get() method will block until completed
				futures.get(t).get();
			}
			catch (Exception e)
			{
				// This should not happen. 
				e.printStackTrace();
			}
		}
		futures.clear();
	}

	/**
	 * For processing the unmixing of specific channel activations
	 */
	private class UnmixWorker implements Runnable
	{
		final RunSettings runSettings;
		final int[][] density;
		final int[] newChannel;
		final int from;
		final int to;
		RandomGenerator random;

		public UnmixWorker(RunSettings runSettings, int[][] density, int[] newChannel, int from, int to, long seed)
		{
			this.runSettings = runSettings;
			this.density = density;
			this.newChannel = newChannel;
			this.from = from;
			this.to = to;
			if (runSettings.crosstalkCorrection != CrosstalkCorrection.SUBTRACTION)
				random = new Well19937c(seed);
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Runnable#run()
		 */
		public void run()
		{
			double[] p = new double[channels];
			for (int i = from; i < to; i++)
			{
				// Current channel
				int c = specificActivations[i].channelx;

				// Observed density
				int[] D = density[i];

				// Compute the true local densities
				double[] d;
				if (channels == 2)
				{
					d = unmix(D[0], D[1], ct[C12], ct[C12]);
				}
				else
				{
					d = unmix(D[0], D[1], D[2], ct[C21], ct[C31], ct[C12], ct[C32], ct[C13], ct[C23]);
				}

				// Apply crosstalk correction
				if (runSettings.crosstalkCorrection == CrosstalkCorrection.SUBTRACTION)
				{
					// Compute the probability it is correct:
					double pc = d[c] / D[c];

					// Remove it if below the subtraction threshold
					if (pc < runSettings.subtractionCutoff[c])
						c = -1;
				}
				else
				{
					// Switch.
					// Compute the probability of each channel and perform a weighted random selection
					double sum = 0;
					for (int j = channels; j-- > 0;)
					{
						sum += (p[j] = d[j] / D[j]);
					}
					final double sum2 = sum * random.nextDouble();
					sum = 0;
					for (c = channels; c-- > 0;)
					{
						sum += p[c];
						if (sum >= sum2)
							break;
					}
				}

				newChannel[i] = c;
			}
		}
	}

	/**
	 * For processing the unmixing of specific channel activations
	 */
	private class NonSpecificUnmixWorker implements Runnable
	{
		final RunSettings runSettings;
		final DensityCounter dc;
		final int[] newChannel;
		final int from;
		final int to;
		RandomGenerator random;

		public NonSpecificUnmixWorker(RunSettings runSettings, DensityCounter dc, int[] newChannel, int from, int to,
				long seed)
		{
			this.runSettings = runSettings;
			this.dc = dc;
			this.newChannel = newChannel;
			this.from = from;
			this.to = to;
			random = new Well19937c(seed);
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Runnable#run()
		 */
		public void run()
		{
			// TODO - We could do other non-specific assignments.
			// e.g. Compute probability for each channel and assign
			// using a weighted random selection

			int[] c = new int[channels];
			for (int i = from; i < to; i++)
			{
				// Assigned channel
				c[0] = -1;

				// Assume the observed density is the true local density 
				// (i.e. cross talk correction of specific activations is perfect)
				int[] d = dc.count(nonSpecificActivations[i], channels);

				// Compute the probability of each channel as:
				// p(i) = di / (d1 + d2 + ... + dn)
				double sum = 0;
				for (int j = channels; j-- > 0;)
				{
					sum += d[j];
				}

				// Assign to most probable channel above the cut-off
				double max = runSettings.nonSpecificAssignmentCutoff;
				int size = 0;
				for (int j = channels; j-- > 0;)
				{
					double p = d[j] / sum;
					if (p > max)
					{
						size = 1;
						max = p;
						c[0] = j;
					}
					else if (p == max)
					{
						// Equal so store all for a random pick
						c[size++] = j;
					}
				}

				newChannel[i] = (size > 1) ? c[random.nextInt(size)] : c[0];
			}
		}
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

	private int write(PeakResultsList[] output, Activation[] activations, int count)
	{
		for (int i = activations.length; i-- > 0;)
		{
			Activation result = activations[i];
			if (result.hasSpecificChannel())
			{
				count++;
				output[result.getChannel()].addAll(result.trace.getPoints());
			}
		}
		return count;
	}

	private ImageProcessor getImage(PeakResultsList peakResultsList)
	{
		PeakResults[] list = peakResultsList.toArray();
		IJImagePeakResults image = (IJImagePeakResults) list[1];
		return image.getImagePlus().getProcessor();
	}
}
