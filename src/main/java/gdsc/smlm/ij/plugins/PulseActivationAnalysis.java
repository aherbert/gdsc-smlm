package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.Rectangle;
import java.awt.Window;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.UnitSphereRandomVectorGenerator;
import org.apache.commons.math3.random.Well19937c;

import gdsc.core.clustering.DensityCounter;
import gdsc.core.clustering.DensityCounter.Molecule;
import gdsc.core.ij.Utils;
import gdsc.core.utils.Maths;
import gdsc.core.utils.NotImplementedException;
import gdsc.core.utils.Random;
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
import gdsc.smlm.results.Calibration;
import gdsc.smlm.results.Cluster.CentroidMethod;
import gdsc.smlm.results.IdPeakResult;
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
import ij.process.ImageStatistics;

/**
 * Perform multi-channel super-resolution imaging by means of photo-switchable probes and pulsed light activation.
 * 
 * This plugin is based on the methods described in: Mark Bates, Bo Huang, Graham T. Dempsey, Xiaowei Zhuang (2007).
 * Multicolor Super-Resolution Imaging with Photo-Switchable Fluorescent Probes. Science 317, 1749. DOI:
 * 10.1126/science.1146598.
 */
public class PulseActivationAnalysis implements PlugIn, DialogListener
{
	private String TITLE = "Activation Analysis";

	private enum Correction
	{
		//@formatter:off
		NONE{ public String getName() { return "None"; }},
		SUBTRACTION{ public String getName() { return "Subtraction"; }},
		MOST_LIKELY{ public String getName() { return "Most likely"; }},
		WEIGHTED_RANDOM{ public String getName() { return "Weighted random"; }};
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

	private static Correction[] specificCorrection, nonSpecificCorrection;
	static
	{
		EnumSet<Correction> correction = EnumSet.allOf(Correction.class);
		specificCorrection = correction.toArray(new Correction[correction.size()]);
		correction.remove(Correction.SUBTRACTION);
		nonSpecificCorrection = correction.toArray(new Correction[correction.size()]);
	}

	private enum SimulationDistribution
	{
		//@formatter:off
		POINT{ public String getName() { return "Point"; }},
		LINE{ public String getName() { return "Line"; }},
		CIRCLE{ public String getName() { return "Circle"; }};
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

	private abstract class Shape
	{
		float x, y;

		Shape(float x, float y)
		{
			this.x = x;
			this.y = y;
		}

		boolean canSample()
		{
			return true;
		}

		float[] getPosition()
		{
			return new float[] { x, y };
		}

		abstract float[] sample(RandomGenerator rand);
	}

	private class Point extends Shape
	{
		float[] xy;

		Point(float x, float y)
		{
			super(x, y);
			xy = super.getPosition();
		}

		@Override
		boolean canSample()
		{
			return false;
		}

		@Override
		float[] getPosition()
		{
			return xy;
		}

		@Override
		float[] sample(RandomGenerator rand)
		{
			throw new NotImplementedException();
		}
	}

	private class Line extends Shape
	{
		double radius, length;
		float sina, cosa;

		/**
		 * Instantiates a new line.
		 *
		 * @param x
		 *            the x
		 * @param y
		 *            the y
		 * @param angle
		 *            the angle (in radians)
		 * @param radius
		 *            the radius
		 */
		Line(float x, float y, double angle, double radius)
		{
			super(x, y);
			this.radius = radius;
			length = 2 * radius;
			sina = (float) Math.sin(angle);
			cosa = (float) Math.cos(angle);
		}

		@Override
		float[] sample(RandomGenerator rand)
		{
			float p = (float) (-radius + rand.nextDouble() * length);
			return new float[] { sina * p + x, cosa * p + y };
		}
	}

	private class Circle extends Shape
	{
		double radius;

		Circle(float x, float y, double radius)
		{
			super(x, y);
			this.radius = radius;
		}

		@Override
		float[] sample(RandomGenerator rand)
		{
			double[] v = new UnitSphereRandomVectorGenerator(2, rand).nextVector();
			return new float[] { (float) (v[0] * radius + x), (float) (v[1] * radius + y) };
		}

	}

	private static String inputOption = "";
	// Note: Set defaults to work with the 3-channel simulation
	private static int channels = 3;
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
	private static int specificCorrectionIndex = Correction.SUBTRACTION.ordinal();
	private static double[] specificCorrectionCutoff = { 50, 50, 50 };
	private static int nonSpecificCorrectionIndex = 0;
	private static double nonSpecificCorrectionCutoff = 50;

	// Simulation settings
	private RandomDataGenerator rdg = null;

	private RandomDataGenerator getRandomDataGenerator()
	{
		if (rdg == null)
			rdg = new RandomDataGenerator(new Well19937c());
		return rdg;
	}

	private static int[] sim_nMolecules = { 1000, 1000, 1000 };
	private static SimulationDistribution[] sim_distribution = { SimulationDistribution.CIRCLE,
			SimulationDistribution.LINE, SimulationDistribution.POINT };
	private static double sim_precision[] = { 15, 15, 15 }; // nm
	private static int sim_cycles = 1000;
	private static int sim_size = 256;
	private static double sim_nmPerPixel = 100;
	private static double sim_activationDensity = 0.1; // molecules/micrometer
	private static double sim_nonSpecificFrequency = 0.01;

	private GlobalSettings settings;
	private ResultsSettings resultsSettings;
	private MemoryPeakResults results;
	private Trace[] traces;

	private class Activation implements Molecule
	{
		final Trace trace;
		float x, y;
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

		boolean hasChannel()
		{
			return channel != 0;
		}

		int getChannel()
		{
			return channel - 1;
		}

		boolean hasCurrentChannel()
		{
			return currentChannel != 0;
		}

		int getCurrentChannel()
		{
			return currentChannel - 1;
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
			// Note: the ID must be zero or above
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

		switch (isSimulation())
		{
			case 0:
				// Most common to not run the simulation
				break;
			case -1:
				// Cancelled
				return;
			case 1:
				//OK'd
				runSimulation();
				return;
		}

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

		if (!results.isCalibrated())
		{
			IJ.error(TITLE, "Results must have basic calibration (pixel pitch and gain)");
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
		TITLE = ((crosstalkMode) ? "Crosstalk " : "Pulse ") + TITLE;

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
						activations.add(new Activation(current, channel));

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
				activations.add(new Activation(current, channel));
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
			if (result.hasChannel())
				count[result.getChannel()]++;
		}

		// Store specific activations
		int sum = (int) Maths.sum(count);
		specificActivations = new Activation[sum];
		nonSpecificActivations = new Activation[list.size() - sum];
		for (int i = list.size(), c1 = 0, c2 = 0; i-- > 0;)
		{
			Activation result = list.getf(i);
			if (result.hasChannel())
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

		double[] crosstalk = computeCrosstalk(count, targetChannel - 1);

		// Store the cross talk.
		// Crosstalk from M into N is defined as the number of times the molecule that should be 
		// activated by a pulse from channel M is activated by a pulse from channel N.
		// targetChannel = M
		// activationChannel = N
		int index1, index2 = -1;
		if (channels == 2)
		{
			if (targetChannel == 1)
				index1 = setCrosstalk(C12, crosstalk[1]);
			else
				index1 = setCrosstalk(C21, crosstalk[0]);
		}
		else
		{
			// 3-channel
			if (targetChannel == 1)
			{
				index1 = setCrosstalk(C12, crosstalk[1]);
				index2 = setCrosstalk(C13, crosstalk[2]);
			}
			else if (targetChannel == 2)
			{
				index1 = setCrosstalk(C21, crosstalk[0]);
				index2 = setCrosstalk(C23, crosstalk[2]);
			}
			else
			{
				index1 = setCrosstalk(C31, crosstalk[0]);
				index2 = setCrosstalk(C32, crosstalk[1]);
			}
		}

		// Show fraction activations histogram. So we have to set the sum to 1
		double sum = Maths.sum(crosstalk);
		for (int i = 0; i < crosstalk.length; i++)
			crosstalk[i] /= sum;

		// Plot a histogram
		double[] x = Utils.newArray(channels, 0.5, 1);
		double[] y = crosstalk;
		Plot2 plot = new Plot2(TITLE, "Channel", "Fraction activations");
		plot.setLimits(0, channels + 1, 0, 1);
		plot.setXMinorTicks(false);
		plot.addPoints(x, y, Plot2.BAR);
		String label = String.format("Crosstalk %s = %s", ctNames[index1], Maths.round(ct[index1]));
		if (index2 > -1)
			label += String.format(", %s = %s", ctNames[index2], Maths.round(ct[index2]));
		plot.addLabel(0, 0, label);
		Utils.display(TITLE, plot);
	}

	/**
	 * Compute crosstalk.
	 * <p>
	 * "The crosstalk ratios can be calculated from the ratios of
	 * incorrectly to correctly colored localizations."
	 *
	 * @param count
	 *            the count
	 * @param target
	 *            the target
	 * @return the double[]
	 */
	private static double[] computeCrosstalk(int[] count, int target)
	{
		double[] crosstalk = new double[count.length];
		for (int c = 0; c < count.length; c++)
			crosstalk[c] = (double) count[c] / count[target];
		return crosstalk;
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
				"Crosstalk analysis requires a sample singly labelled with only one photo-switchable probe and imaged with the full pulse lifecycle. The probe should be activated by the pulse in the target channel. Activations from the pulse in other channels is crosstalk.",
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
	 * channel M is activated by a pulse from channel N. A value less than 1 is expected (otherwise the fluorophore is
	 * not being specifically activated by channel M).
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
		// between 0 and 1 then we do not need to check the bounds.
		//d1 = Maths.clip(0, D1, d1);
		//d2 = Maths.clip(0, D2, d2);
		return new double[] { d1, d2 };
	}

	/**
	 * Unmix the observed local densities into the actual densities for 3-channels.
	 * <p>
	 * Crosstalk from M into N is defined as the number of times the molecule that should be activated by a pulse from
	 * channel M is activated by a pulse from channel N. A value less than 1 is expected (otherwise the fluorophore is
	 * not being specifically activated by channel M).
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

		// Use matrix inversion so that: X = A^-1 * B 
		double a = 1;
		double b = C21;
		double c = C31;
		double d = C12;
		double e = 1;
		double f = C32;
		double g = C13;
		double h = C23;
		double i = 1;

		double A = (e * i - f * h);
		double B = -(d * i - f * g);
		double C = (d * h - e * g);

		double det = a * A + b * B + c * C;

		double det_recip = 1.0 / det;

		if (!Maths.isFinite(det_recip))
			// Failed so reset to the observed densities
			return new double[] { D1, D2, D3 };

		double D = -(b * i - c * h);
		double E = (a * i - c * g);
		double F = -(a * h - b * g);
		double G = (b * f - c * e);
		double H = -(a * f - c * d);
		double I = (a * e - b * d);

		a = det_recip * A;
		b = det_recip * D;
		c = det_recip * G;
		d = det_recip * B;
		e = det_recip * E;
		f = det_recip * H;
		g = det_recip * C;
		h = det_recip * F;
		i = det_recip * I;

		final double[] x = new double[3];
		x[0] = a * D1 + b * D2 + c * D3;
		x[1] = d * D1 + e * D2 + f * D3;
		x[2] = g * D1 + h * D2 + i * D3;

		// Use matrix decomposition
		//		// Note: The linear solver uses LU decomposition. 
		//		// We cannot use a faster method as the matrix A is not symmetric.
		//		final LinearSolver<DenseMatrix64F> linearSolver = LinearSolverFactory.linear(3);
		//		final DenseMatrix64F A = new DenseMatrix64F(3, 3);
		//		A.set(0, 1);
		//		A.set(1, C21);
		//		A.set(2, C31);
		//		A.set(3, C12);
		//		A.set(4, 1);
		//		A.set(5, C32);
		//		A.set(6, C13);
		//		A.set(7, C23);
		//		A.set(8, 1);
		//		final DenseMatrix64F B = new DenseMatrix64F(3, 1);
		//		B.set(0, D1);
		//		B.set(1, D2);
		//		B.set(2, D3);
		//
		//		if (!linearSolver.setA(A))
		//		{
		//			// Failed so reset to the observed densities
		//			return new double[] { D1, D2, D3 };
		//		}
		//
		//		// Input B is not modified to so we can re-use for output X
		//		linearSolver.solve(B, B);
		//
		//		final double[] x = B.getData();

		// Due to floating-point error in the decomposition we check the bounds
		x[0] = Maths.clip(0, D1, x[0]);
		x[1] = Maths.clip(0, D2, x[1]);
		x[2] = Maths.clip(0, D3, x[2]);

		return x;
	}

	private class RunSettings
	{
		double densityRadius;
		Correction specificCorrection = Correction.NONE;
		double[] specificCorrectionCutoff;
		Correction nonSpecificCorrection = Correction.NONE;
		double nonSpecificCorrectionCutoff;

		@SuppressWarnings("unused")
		ResultsSettings resultsSettings;

		RunSettings()
		{
			// Copy the current settings required
			if (channels > 1)
			{
				this.densityRadius = PulseActivationAnalysis.densityRadius / results.getNmPerPixel();

				specificCorrection = getCorrection(PulseActivationAnalysis.specificCorrection,
						PulseActivationAnalysis.specificCorrectionIndex);
				this.specificCorrectionCutoff = new double[channels];
				for (int i = channels; i-- > 0;)
					// Convert from percentage to a probability
					this.specificCorrectionCutoff[i] = PulseActivationAnalysis.specificCorrectionCutoff[i] / 100.0;

				nonSpecificCorrection = getCorrection(PulseActivationAnalysis.nonSpecificCorrection,
						PulseActivationAnalysis.nonSpecificCorrectionIndex);
				this.nonSpecificCorrectionCutoff = PulseActivationAnalysis.nonSpecificCorrectionCutoff / 100.0;
			}
			this.resultsSettings = PulseActivationAnalysis.this.resultsSettings.clone();
		}

		Correction getCorrection(Correction[] correction, int index)
		{
			if (index >= 0 && index < correction.length)
				return correction[index];
			return Correction.NONE;
		}

		public boolean newUnmixSettings(RunSettings lastRunSettings)
		{
			if (lastRunSettings == null)
				return true;
			if (lastRunSettings.densityRadius != densityRadius)
				return true;
			if (lastRunSettings.specificCorrection != specificCorrection)
				return true;
			if (specificCorrection != Correction.NONE)
			{
				for (int i = channels; i-- > 0;)
					if (lastRunSettings.specificCorrectionCutoff[i] != specificCorrectionCutoff[i])
						return true;
			}
			return false;
		}

		public boolean newNonSpecificCorrectionSettings(RunSettings lastRunSettings)
		{
			if (lastRunSettings == null)
				return true;
			if (lastRunSettings.densityRadius != densityRadius)
				return true;
			if (lastRunSettings.nonSpecificCorrection != nonSpecificCorrection)
				return true;
			if (nonSpecificCorrection != Correction.NONE)
			{
				if (lastRunSettings.nonSpecificCorrectionCutoff != nonSpecificCorrectionCutoff)
					return true;
			}
			return false;
		}
	}

	// Here we use a simple workflow with only one worker since the results are 
	// written straight back to this class' objects
	private Workflow<RunSettings, Object> workflow;

	private void runPulseAnalysis()
	{
		// Use a simple workflow with one worker
		workflow = new Workflow<RunSettings, Object>();
		workflow.add(new WorkflowWorker<RunSettings, Object>()
		{
			@Override
			public boolean equalSettings(RunSettings current, RunSettings previous)
			{
				return false;
			}

			@Override
			public boolean equalResults(Object current, Object previous)
			{
				return false;
			}

			@Override
			public Object createResults(RunSettings settings, Object results)
			{
				PulseActivationAnalysis.this.run(settings);
				return null;
			}
		});

		workflow.start();

		boolean cancelled = !showPulseAnalysisDialog();

		workflow.shutdown(cancelled);

		if (executor != null)
		{
			// Shutdown the executor used to do the work

			if (cancelled)
				// Stop immediately
				executor.shutdownNow();
			else
				executor.shutdown();
		}
	}

	private boolean showPulseAnalysisDialog()
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
			correctionNames = SettingsManager.getNames((Object[]) specificCorrection);
			gd.addChoice("Crosstalk_correction", correctionNames, correctionNames[specificCorrectionIndex]);
			for (int c = 1; c <= channels; c++)
				gd.addSlider("Crosstalk_correction_cutoff_C" + c + "(%)", 0, 100, specificCorrectionCutoff[c - 1]);
			assigmentNames = SettingsManager.getNames((Object[]) nonSpecificCorrection);
			gd.addChoice("Nonspecific_assigment", assigmentNames, assigmentNames[nonSpecificCorrectionIndex]);
			gd.addSlider("Nonspecific_assignment_cutoff (%)", 0, 100, nonSpecificCorrectionCutoff);
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

		// The dialog was OK'd so run if work was staged in the workflow.
		if (workflow.isStaged())
			workflow.runStaged();

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

				Recorder.recordOption("Crosstalk_correction", correctionNames[specificCorrectionIndex]);
				for (int c = 1; c <= channels; c++)
					Recorder.recordOption("Crosstalk_correction_cutoff_C" + c,
							Double.toString(specificCorrectionCutoff[c - 1]));

				Recorder.recordOption("Nonspecific_assigment", assigmentNames[nonSpecificCorrectionIndex]);
				Recorder.recordOption("Nonspecific_assignment_cutoff (%)",
						Double.toString(nonSpecificCorrectionCutoff));
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
				specificCorrectionIndex = gd.getNextChoiceIndex();
				for (int c = 1; c <= channels; c++)
				{
					specificCorrectionCutoff[c - 1] = (int) gd.getNextNumber();
					validatePercentage("Crosstalk_correction_cutoff_C" + c, specificCorrectionCutoff[c - 1]);
				}
				nonSpecificCorrectionIndex = gd.getNextChoiceIndex();
				nonSpecificCorrectionCutoff = gd.getNextNumber();
				validatePercentage("Nonspecific_assignment_cutoff", nonSpecificCorrectionCutoff);
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

		RunSettings settings = new RunSettings();
		if (preview)
		{
			// Run the settings
			workflow.run(settings);
			workflow.startPreview();
		}
		else
		{
			workflow.stopPreview();
			// Stage the work but do not run
			workflow.stage(settings);
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

	private DensityCounter dc = null;
	private int[][] density = null;
	private int nThreads;
	private ExecutorService executor = null;
	private TurboList<Future<?>> futures = null;
	private RunSettings lastRunSettings = null;

	private synchronized void run(RunSettings runSettings)
	{
		// This is synchronized since it updates the class results. 
		// Note: We check against the last settings and only repeat what is necessary ...

		if (runSettings == null)
		{
			lastRunSettings = null;
			return;
		}

		IJ.showStatus("Analysing ...");

		// Assign all activations to a channel.
		// This is only necessary when we have more than 1 channel. If we have 1 channel then 
		// no correction method is specified.
		boolean changed = false;
		if (runSettings.newUnmixSettings(lastRunSettings))
		{
			changed = true;

			// Reset
			for (int i = specificActivations.length; i-- > 0;)
			{
				Activation result = specificActivations[i];
				result.currentChannel = result.channel;
			}

			if (runSettings.specificCorrection != Correction.NONE)
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
					futures.add(executor
							.submit(new SpecificUnmixWorker(runSettings, density, newChannel, from, to, seed++)));
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
		if (changed || runSettings.newNonSpecificCorrectionSettings(lastRunSettings))
		{
			// Reset
			for (int i = nonSpecificActivations.length; i-- > 0;)
			{
				Activation result = nonSpecificActivations[i];
				result.currentChannel = result.channel;
			}

			if (runSettings.nonSpecificCorrection != Correction.NONE)
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
		count = write(output, nonSpecificActivations, count);

		int size = 0;
		for (int c = 0; c < channels; c++)
		{
			output[c].end();
			size += output[c].size();
		}

		// Collate image into a stack
		if (channels > 1 && resultsSettings.getResultsImage() != ResultsImage.NONE)
		{
			ImageProcessor[] images = new ImageProcessor[channels];
			for (int c = 0; c < channels; c++)
				images[c] = getImage(output[c]);
			displayComposite(images, results.getName() + " " + TITLE);
		}

		lastRunSettings = runSettings;
		runSettings = null;

		IJ.showStatus(String.format("%d/%s, %d/%s", count, Utils.pleural(traces.length, "Trace"), size,
				Utils.pleural(results.size(), "Result")));
	}

	private void displayComposite(ImageProcessor[] images, String name)
	{
		ImageStack stack = null; // We do not yet know the size
		for (int i = 0; i < images.length; i++)
		{
			ImageProcessor ip = images[i];
			if (stack == null)
				stack = new ImageStack(ip.getWidth(), ip.getHeight());
			ip.setColorModel(null);
			stack.addSlice("C" + (i + 1), ip);
		}

		// Create a composite
		ImagePlus imp = new ImagePlus(name, stack);
		imp.setDimensions(images.length, 1, 1);
		CompositeImage ci = new CompositeImage(imp, IJ.COMPOSITE);

		// Make it easier to see
		//ij.plugin.ContrastEnhancerce = new ij.plugin.ContrastEnhancer();
		//double saturated = 0.35;
		//ce.stretchHistogram(ci, saturated);

		autoAdjust(ci, ci.getProcessor());

		imp = WindowManager.getImage(name);
		if (imp != null && imp.isComposite())
		{
			ci.setMode(imp.getCompositeMode());
			imp.setImage(ci);
			imp.getWindow().toFront();
		}
		else
		{
			ci.show();
			imp = ci;
		}

		if (WindowManager.getWindow("Channels") == null)
		{
			IJ.run("Channels Tool...");
			Window w = WindowManager.getWindow("Channels");
			if (w == null)
				return;
			Window w2 = imp.getWindow();
			if (w2 == null)
				return;
			java.awt.Point p = w2.getLocation();
			p.x += w2.getWidth();
			w.setLocation(p);
		}
	}

	/**
	 * Auto adjust. Copied from ij.plgin.frame.ContrastAdjuster
	 * <p>
	 * Although the ContrastAdjuster records its actions as 'run("Enhance Contrast", "saturated=0.35");' it actually
	 * does something else which makes the image easier to see than the afore mentioned command.
	 *
	 * @param imp
	 *            the imp
	 * @param ip
	 *            the ip
	 */
	private void autoAdjust(ImagePlus imp, ImageProcessor ip)
	{
		ij.measure.Calibration cal = imp.getCalibration();
		imp.setCalibration(null);
		ImageStatistics stats = imp.getStatistics(); // get uncalibrated stats
		imp.setCalibration(cal);
		int limit = stats.pixelCount / 10;
		int[] histogram = stats.histogram;
		int autoThreshold = 0;
		if (autoThreshold < 10)
			autoThreshold = 5000;
		else
			autoThreshold /= 2;
		int threshold = stats.pixelCount / autoThreshold;
		int i = -1;
		boolean found = false;
		int count;
		do
		{
			i++;
			count = histogram[i];
			if (count > limit)
				count = 0;
			found = count > threshold;
		} while (!found && i < 255);
		int hmin = i;
		i = 256;
		do
		{
			i--;
			count = histogram[i];
			if (count > limit)
				count = 0;
			found = count > threshold;
		} while (!found && i > 0);
		int hmax = i;
		if (hmax >= hmin)
		{
			double min = stats.histMin + hmin * stats.binSize;
			double max = stats.histMin + hmax * stats.binSize;
			if (min == max)
			{
				min = stats.min;
				max = stats.max;
			}
			imp.setDisplayRange(min, max);
		}
		else
		{
			reset(imp, ip);
			return;
		}
	}

	void reset(ImagePlus imp, ImageProcessor ip)
	{
		int bitDepth = imp.getBitDepth();
		double defaultMin, defaultMax;
		if (bitDepth == 16 || bitDepth == 32)
		{
			imp.resetDisplayRange();
			defaultMin = imp.getDisplayRangeMin();
			defaultMax = imp.getDisplayRangeMax();
		}
		else
		{
			defaultMin = 0;
			defaultMax = 255;
		}
		imp.setDisplayRange(defaultMin, defaultMax);
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

	private abstract class UnmixWorker
	{
		final int[] newChannel;
		final int from;
		final int to;
		RandomGenerator random;
		int[] assignedChannel = new int[channels];
		double[] p = new double[channels];

		public UnmixWorker(int[] newChannel, int from, int to, long seed)
		{
			this.newChannel = newChannel;
			this.from = from;
			this.to = to;
			random = new Well19937c(seed);
		}

		int weightedRandomSelection(double cutoff)
		{
			double sum = 0;
			for (int j = channels; j-- > 0;)
			{
				if (p[j] > cutoff)
					sum += p[j];
				else
					p[j] = 0;
			}
			if (sum == 0)
				return 0;

			final double sum2 = sum * random.nextDouble();
			sum = 0;
			for (int j = channels; j-- > 0;)
			{
				sum += p[j];
				if (sum >= sum2)
				{
					return j + 1;
				}
			}
			// This should not happen
			return 0;
		}

		int mostLikelySelection(double cutoff)
		{
			double max = cutoff;
			int size = 0;
			for (int j = channels; j-- > 0;)
			{
				if (p[j] > max)
				{
					size = 1;
					max = p[j];
					assignedChannel[0] = j;
				}
				else if (p[j] == max)
				{
					// Equal so store all for a random pick
					assignedChannel[size++] = j;
				}
			}

			if (size == 0)
				return 0;

			return (size > 1) ? assignedChannel[random.nextInt(size)] + 1 : assignedChannel[0] + 1;
		}
	}

	/**
	 * For processing the unmixing of specific channel activations
	 */
	private class SpecificUnmixWorker extends UnmixWorker implements Runnable
	{
		final RunSettings runSettings;
		final int[][] density;

		public SpecificUnmixWorker(RunSettings runSettings, int[][] density, int[] newChannel, int from, int to,
				long seed)
		{
			super(newChannel, from, to, seed);
			this.runSettings = runSettings;
			this.density = density;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Runnable#run()
		 */
		public void run()
		{
			for (int i = from; i < to; i++)
			{
				// Current channel (1-indexed)
				int c = specificActivations[i].channel;

				// Observed density
				int[] D = density[i];

				// Compute the true local densities
				double[] d;
				if (channels == 2)
				{
					d = unmix(D[1], D[2], ct[C12], ct[C12]);
				}
				else
				{
					d = unmix(D[1], D[2], D[3], ct[C21], ct[C31], ct[C12], ct[C32], ct[C13], ct[C23]);
				}

				// Apply crosstalk correction
				if (runSettings.specificCorrection == Correction.SUBTRACTION)
				{
					// Compute the probability it is correct:
					// This is a measure of how much crosstalk effected the observed density.
					// (This is taken from Bates et al, 2007)
					double pc = d[c - 1] / D[c];

					// Remove it if below the subtraction threshold
					if (pc < runSettings.specificCorrectionCutoff[c - 1])
						c = 0;
				}
				else
				{
					// Compute the probability of each channel as:
					// p(i) = di / (d1 + d2 + ... + dn)
					// Note this is different from computing the probability of the channel being correct.
					// That probability is an indication of how much crosstalk has effected the observed density.
					// This value is a simple probability using the local density in each channel. 
					double sum = 0;
					for (int j = channels; j-- > 0;)
						sum += d[j];
					// Note that since this is a specific activation we can assume the molecule will be 
					// self-counted within the radius and that d will never be zero in every channel
					for (int j = channels; j-- > 0;)
						p[j] = d[j] / sum;

					if (runSettings.specificCorrection == Correction.WEIGHTED_RANDOM)
						c = weightedRandomSelection(runSettings.specificCorrectionCutoff[c - 1]);
					else //if (runSettings.specificCorrection == CrosstalkCorrection.MOST_LIKELY)
						c = mostLikelySelection(runSettings.specificCorrectionCutoff[c - 1]);
				}

				newChannel[i] = c;
			}
		}
	}

	/**
	 * For processing the unmixing of specific channel activations
	 */
	private class NonSpecificUnmixWorker extends UnmixWorker implements Runnable
	{
		final RunSettings runSettings;
		final DensityCounter dc;

		public NonSpecificUnmixWorker(RunSettings runSettings, DensityCounter dc, int[] newChannel, int from, int to,
				long seed)
		{
			super(newChannel, from, to, seed);
			this.runSettings = runSettings;
			this.dc = dc;
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

			for (int i = from; i < to; i++)
			{
				int c = 0;

				// Assume the observed density is the true local density 
				// (i.e. cross talk correction of specific activations is perfect)
				int[] d = dc.count(nonSpecificActivations[i], channels);

				// Compute the probability of each channel as:
				// p(i) = di / (d1 + d2 + ... + dn)
				double sum = 0;
				for (int j = 1; j <= channels; j++)
				{
					sum += d[j];
				}

				if (sum != 0)
				{
					for (int j = channels; j-- > 0;)
						p[j] = d[j + 1] / sum;

					if (runSettings.nonSpecificCorrection == Correction.WEIGHTED_RANDOM)
						c = weightedRandomSelection(runSettings.nonSpecificCorrectionCutoff);
					else //if (runSettings.nonSpecificCorrection == NonSpecificCorrection.MOST_LIKELY)
						c = mostLikelySelection(runSettings.nonSpecificCorrectionCutoff);
				}

				newChannel[i] = c;
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
			if (result.hasCurrentChannel())
			{
				count++;
				output[result.getCurrentChannel()].addAll(result.trace.getPoints());
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

	private int isSimulation()
	{
		if (Utils.isExtraOptions())
		{
			GenericDialog gd = new GenericDialog(TITLE);
			gd.addMessage("Perform a crosstalk simulation?");
			gd.enableYesNoCancel();
			gd.showDialog();
			if (gd.wasOKed())
				return 1;
			if (gd.wasCanceled())
				return -1;
		}
		return 0;
	}

	private void runSimulation()
	{
		TITLE += " Simulation";

		if (!showSimulationDialog())
			return;

		long start = System.currentTimeMillis();
		RandomDataGenerator rdg = getRandomDataGenerator();

		// Draw the molecule positions
		Utils.showStatus("Simulating molecules ...");
		float[][][] molecules = new float[3][][];
		MemoryPeakResults[] results = new MemoryPeakResults[3];
		Rectangle bounds = new Rectangle(0, 0, sim_size, sim_size);
		for (int c = 0; c < 3; c++)
		{
			molecules[c] = simulateMolecules(rdg, c);

			// Create a dataset to store the activations
			MemoryPeakResults r = new MemoryPeakResults();
			r.setCalibration(new Calibration(sim_nmPerPixel, 1, 100));
			r.setBounds(bounds);
			r.setName(TITLE + " C" + (c + 1));
			results[c] = r;
		}

		// Simulate activation
		Utils.showStatus("Simulating activations ...");
		for (int c = 0; c < 3; c++)
			simulateActivations(rdg, molecules, c, results);

		// Combine
		Utils.showStatus("Producing simulation output ...");
		MemoryPeakResults r = new MemoryPeakResults();
		r.setCalibration(new Calibration(sim_nmPerPixel, 1, 100));
		r.setBounds(new Rectangle(0, 0, sim_size, sim_size));
		r.setName(TITLE);

		ImageProcessor[] images = new ImageProcessor[3];
		for (int c = 0; c < 3; c++)
		{
			ArrayList<PeakResult> list = (ArrayList<PeakResult>) results[c].getResults();
			r.addAllf(list);

			// Draw the unmixed activations
			IJImagePeakResults image = ImagePeakResultsFactory.createPeakResultsImage(ResultsImage.LOCALISATIONS, true,
					true, TITLE, bounds, sim_nmPerPixel, 1, 1024.0 / sim_size, 0, ResultsMode.ADD);
			image.setLiveImage(false);
			image.setDisplayImage(false);
			image.begin();
			image.addAll(list);
			image.end();
			images[c] = image.getImagePlus().getProcessor();
		}
		displayComposite(images, TITLE);

		// Add to memory. Set the composite dataset first.
		MemoryPeakResults.addResults(r);
		for (int c = 0; c < 3; c++)
			MemoryPeakResults.addResults(results[c]);

		// TODO:
		// Show an image of what it looks like with no unmixing, i.e. colours allocated 
		// from the frame

		Utils.showStatus("Simulation complete: " + Utils.timeToString(System.currentTimeMillis() - start));
	}

	private float[][] simulateMolecules(RandomDataGenerator rdg, int c)
	{
		int n = sim_nMolecules[c];
		float[][] molecules = new float[n][];
		if (n == 0)
			return molecules;

		// Draw the shapes
		Shape[] shapes = createShapes(rdg, c);

		// Sample positions from within the shapes
		boolean canSample = shapes[0].canSample();
		RandomGenerator rand = rdg.getRandomGenerator();
		while (n-- > 0)
		{
			float[] coords;
			if (canSample)
			{
				int next = rand.nextInt(shapes.length);
				coords = shapes[next].sample(rand);
			}
			else
			{
				coords = shapes[n % shapes.length].getPosition();
			}

			// Avoid out-of-bounds positions
			if (outOfBounds(coords[0]) || outOfBounds(coords[1]))
				n++;
			else
				molecules[n] = coords;
		}
		return molecules;
	}

	private Shape[] createShapes(RandomDataGenerator rdg, int c)
	{
		RandomGenerator rand = rdg.getRandomGenerator();
		Shape[] shapes;
		double min = sim_size / 20;
		double max = sim_size / 10;
		double range = max - min;
		switch (sim_distribution[c])
		{
			case CIRCLE:
				shapes = new Shape[10];
				for (int i = 0; i < shapes.length; i++)
				{
					float x = nextCoordinate(rand);
					float y = nextCoordinate(rand);
					double radius = rand.nextDouble() * range + min;
					shapes[i] = new Circle(x, y, radius);
				}
				break;

			case LINE:
				shapes = new Shape[10];
				for (int i = 0; i < shapes.length; i++)
				{
					float x = nextCoordinate(rand);
					float y = nextCoordinate(rand);
					double angle = rand.nextDouble() * Math.PI;
					double radius = rand.nextDouble() * range + min;
					shapes[i] = new Line(x, y, angle, radius);
				}

				break;

			case POINT:
			default:
				shapes = new Shape[sim_nMolecules[c]];
				for (int i = 0; i < shapes.length; i++)
				{
					float x = nextCoordinate(rand);
					float y = nextCoordinate(rand);
					shapes[i] = new Point(x, y);
				}
		}
		return shapes;
	}

	private float nextCoordinate(RandomGenerator rand)
	{
		return (float) rand.nextDouble() * sim_size;
	}

	private boolean outOfBounds(float f)
	{
		return f < 0 || f > sim_size;
	}

	private void simulateActivations(RandomDataGenerator rdg, float[][][] molecules, int c, MemoryPeakResults[] results)
	{
		int n = molecules[c].length;
		if (n == 0)
			return;

		// Compute desired number per frame
		double umPerPixel = sim_nmPerPixel / 1000;
		double um2PerPixel = umPerPixel * umPerPixel;
		double area = sim_size * sim_size * um2PerPixel;
		double nPerFrame = area * sim_activationDensity;

		// Compute the activation probability (but set an upper limit so not all are on in every frame)
		double p = Math.min(0.5, nPerFrame / n);

		// Determine the other channels activation probability using crosstalk
		double[] p0 = { p, p, p };
		int index1, index2, c1, c2;
		switch (c)
		{
			case 0:
				index1 = C12;
				index2 = C13;
				c1 = 1;
				c2 = 2;
				break;
			case 1:
				index1 = C21;
				index2 = C23;
				c1 = 0;
				c2 = 2;
				break;
			case 2:
			default:
				index1 = C31;
				index2 = C32;
				c1 = 0;
				c2 = 1;
				break;
		}
		p0[c1] *= ct[index1];
		p0[c2] *= ct[index2];

		// Assume 10 frames after each channel pulse => 30 frames per cycle
		double precision = sim_precision[c] / sim_nmPerPixel;
		int id = c + 1;

		RandomGenerator rand = rdg.getRandomGenerator();
		BinomialDistribution[] bd = new BinomialDistribution[4];
		for (int i = 0; i < 3; i++)
			bd[i] = createBinomialDistribution(rand, n, p0[i]);

		int[] frames = new int[27];
		for (int i = 1, j = 0; i <= 30; i++)
		{
			if (i % 10 == 1)
				// Skip specific activation frames 
				continue;
			frames[j++] = i;
		}
		bd[3] = createBinomialDistribution(rand, n, p * sim_nonSpecificFrequency);

		// Count the actual cross talk
		int[] count = new int[3];

		for (int i = 0, t = 1; i < sim_cycles; i++, t += 30)
		{
			count[0] += simulateActivations(rdg, bd[0], molecules[c], results[c], t, precision, id);
			count[1] += simulateActivations(rdg, bd[1], molecules[c], results[c], t + 10, precision, id);
			count[2] += simulateActivations(rdg, bd[2], molecules[c], results[c], t + 20, precision, id);
			// Add non-specific activations
			if (bd[3] != null)
			{
				for (int t2 : frames)
					simulateActivations(rdg, bd[3], molecules[c], results[c], t2, precision, id);
			}
		}

		// Report simulated cross talk
		double[] crosstalk = computeCrosstalk(count, c);
		Utils.log("Simulated crosstalk C%s  %s=>%s, C%s  %s=>%s", ctNames[index1], Utils.rounded(ct[index1]),
				Utils.rounded(crosstalk[c1]), ctNames[index2], Utils.rounded(ct[index2]), Utils.rounded(crosstalk[c2]));
	}

	private BinomialDistribution createBinomialDistribution(RandomGenerator rand, int n, double p)
	{
		if (p == 0)
			return null;
		return new BinomialDistribution(rand, n, p);
	}

	private int simulateActivations(RandomDataGenerator rdg, BinomialDistribution bd, float[][] molecules,
			MemoryPeakResults results, int t, double precision, int id)
	{
		if (bd == null)
			return 0;
		int n = molecules.length;
		int k = bd.sample();
		// Sample
		RandomGenerator rand = rdg.getRandomGenerator();
		int[] sample = Random.sample(k, n, rand);
		while (k-- > 0)
		{
			float[] xy = molecules[sample[k]];
			float x, y;
			do
			{
				x = (float) (xy[0] + rand.nextGaussian() * precision);
			} while (outOfBounds(x));
			do
			{
				y = (float) (xy[1] + rand.nextGaussian() * precision);
			} while (outOfBounds(y));

			results.add(createResult(t, x, y));
		}
		return sample.length;
	}

	private int id = 0;

	private IdPeakResult createResult(int t, float x, float y)
	{
		// We add them as if tracing is perfect. So each peak result has a new ID.
		// This allows the output of the simulation to be used directly by the pulse analysis code.
		IdPeakResult r = new IdPeakResult(t, x, y, 1, 1, ++id);
		r.noise = 1; // So it appears calibrated
		return r;
	}

	private boolean showSimulationDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);

		SimulationDistribution[] distributionValues = SimulationDistribution.values();
		String[] distribution = SettingsManager.getNames((Object[]) distributionValues);

		// Random crosstalk if not set
		if (Maths.max(ct) == 0)
		{
			RandomDataGenerator rdg = getRandomDataGenerator();
			for (int i = 0; i < ct.length; i++)
				ct[i] = rdg.nextUniform(0.05, 0.15); // Have some crosstalk
		}

		// Three channel
		for (int c = 0; c < 3; c++)
		{
			String ch = "_C" + (c + 1);
			gd.addNumericField("Molcules" + ch, sim_nMolecules[c], 0);
			gd.addChoice("Distribution" + ch, distribution, distribution[sim_distribution[c].ordinal()]);
			gd.addNumericField("Precision_" + ch, sim_precision[c], 3);
			gd.addNumericField("Crosstalk_" + ctNames[2 * c], ct[2 * c], 3);
			gd.addNumericField("Crosstalk_" + ctNames[2 * c + 1], ct[2 * c + 1], 3);
		}
		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		int count = 0;
		for (int c = 0; c < 3; c++)
		{
			sim_nMolecules[c] = (int) Math.abs(gd.getNextNumber());
			if (sim_nMolecules[c] > 0)
				count++;
			sim_distribution[c] = distributionValues[gd.getNextChoiceIndex()];
			sim_precision[c] = Math.abs(gd.getNextNumber());
			ct[2 * c] = Math.abs(gd.getNextNumber());
			ct[2 * c + 1] = Math.abs(gd.getNextNumber());
		}

		if (gd.invalidNumber())
			return false;
		if (count < 2)
		{
			IJ.error(TITLE, "Simulation requires at least 2 channels");
			return false;
		}

		try
		{
			for (int i = 0; i < ct.length; i += 2)
			{
				if (sim_nMolecules[i / 2] > 0)
					validateCrosstalk(i, i + 1);
			}
		}
		catch (IllegalArgumentException ex)
		{
			IJ.error(TITLE, ex.getMessage());
			return false;
		}

		return true;
	}
}
