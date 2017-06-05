package gdsc.smlm.ij.plugins;

import java.awt.Checkbox;
import java.awt.Color;
import java.awt.Component;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Rectangle;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.commons.math3.distribution.CustomGammaDistribution;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.distribution.RealDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.random.EmpiricalDistribution;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.SobolSequenceGenerator;
import org.apache.commons.math3.random.Well44497b;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.util.FastMath;

import com.thoughtworks.xstream.XStream;
import com.thoughtworks.xstream.io.xml.DomDriver;

import gdsc.core.clustering.DensityManager;
import gdsc.core.ij.Utils;
import gdsc.core.threshold.AutoThreshold;
import gdsc.core.utils.Maths;
import gdsc.core.utils.NoiseEstimator;
import gdsc.core.utils.Random;
import gdsc.core.utils.Statistics;
import gdsc.core.utils.StoredDataStatistics;
import gdsc.core.utils.TextUtils;
import gdsc.core.utils.UnicodeReader;

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

import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.engine.FitWorker;
import gdsc.smlm.filters.GaussianFilter;
import gdsc.smlm.fitting.FitConfiguration;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.ij.IJImageSource;
import gdsc.smlm.ij.plugins.LoadLocalisations.LocalisationList;
import gdsc.smlm.ij.settings.Atom;
import gdsc.smlm.ij.settings.Compound;
import gdsc.smlm.ij.settings.CreateDataSettings;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.PSFOffset;
import gdsc.smlm.ij.settings.PSFSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.model.ActivationEnergyImageModel;
import gdsc.smlm.model.AiryPSFModel;
import gdsc.smlm.model.AiryPattern;
import gdsc.smlm.model.CompoundMoleculeModel;
import gdsc.smlm.model.DiffusionType;
import gdsc.smlm.model.FixedLifetimeImageModel;
import gdsc.smlm.model.FluorophoreSequenceModel;
import gdsc.smlm.model.GaussianPSFModel;
import gdsc.smlm.model.GridDistribution;
import gdsc.smlm.model.ImageModel;
import gdsc.smlm.model.ImagePSFModel;
import gdsc.smlm.model.LocalisationModel;
import gdsc.smlm.model.LocalisationModelSet;
import gdsc.smlm.model.MaskDistribution;
import gdsc.smlm.model.MaskDistribution3D;
import gdsc.smlm.model.MoleculeModel;
import gdsc.smlm.model.PSFModel;
import gdsc.smlm.model.RadialFalloffIllumination;
import gdsc.smlm.model.RandomGeneratorFactory;
import gdsc.smlm.model.SpatialDistribution;
import gdsc.smlm.model.SpatialIllumination;
import gdsc.smlm.model.SphericalDistribution;
import gdsc.smlm.model.UniformDistribution;
import gdsc.smlm.model.UniformIllumination;
import gdsc.smlm.results.Calibration;
import gdsc.smlm.results.ExtendedPeakResult;
import gdsc.smlm.results.FilePeakResults;
import gdsc.smlm.results.IdPeakResult;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResultsReader;
import gdsc.smlm.utils.XmlUtils;
import gnu.trove.set.hash.TIntHashSet;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.ExtendedGenericDialog;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import ij.plugin.WindowOrganiser;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import ij.text.TextWindow;

/**
 * Creates data using a simulated PSF
 */
public class CreateData implements PlugIn, ItemListener, RandomGeneratorFactory
{
	public static final String TITLE = "Create Data";
	private static final String CREATE_DATA_IMAGE_TITLE = "Localisation Data";

	private static String[] ILLUMINATION = { "Uniform", "Radial" };
	private static int RADIAL = 1;
	private static String[] DISTRIBUTION = { "Uniform RNG", "Uniform Halton", "Uniform Sobol", "Mask", "Grid" };
	//private static int UNIFORM = 0;
	private static int UNIFORM_HALTON = 1;
	private static int UNIFORM_SOBOL = 2;
	private static int MASK = 3;
	private static int GRID = 4;
	private static String[] CONFINEMENT = { "None", "Mask", "Sphere", "Within Image" };
	private static int CONFINEMENT_MASK = 1;
	private static int CONFINEMENT_SPHERE = 2;
	private static int CONFINEMENT_WITHIN_IMAGE = 3;
	private static String[] PHOTON_DISTRIBUTION = { "Uniform", "Gamma", "Custom", "Fixed", "Correlated" };
	private static int PHOTON_UNIFORM = 0;
	private static int PHOTON_GAMMA = 1;
	private static int PHOTON_CUSTOM = 2;
	private static int PHOTON_FIXED = 3;
	private static int PHOTON_CORRELATED = 4;

	private static String[] PSF_MODELS = new String[] { "2D Gaussian", "Airy", "Image" };

	private static TextWindow summaryTable = null;
	private static int datasetNumber = 0;
	private static double areaInUm = 0;
	private static String header = null;

	private GlobalSettings globalSettings;
	private CreateDataSettings settings;

	private static final String[] NAMES = new String[] { "Samples/Frame", "Signal/Frame", "Signal/Frame (continuous)",
			"Total Signal", "Blinks", "t-On", "t-Off", "Sampled blinks", "Sampled t-On", "Sampled t-Off", "Noise",
			"SNR", "SNR (continuous)", "Density", "Precision", "Precision (in-focus)", "X", "Y", "Z", "Width" };
	private static boolean[] displayHistograms = new boolean[NAMES.length];
	static
	{
		for (int i = 0; i < displayHistograms.length; i++)
			displayHistograms[i] = true;
	}
	private static final int SAMPLES = 0;
	private static final int SIGNAL = 1;
	private static final int SIGNAL_CONTINUOUS = 2;
	private static final int TOTAL_SIGNAL = 3;
	private static final int BLINKS = 4;
	private static final int T_ON = 5;
	private static final int T_OFF = 6;
	private static final int SAMPLED_BLINKS = 7;
	private static final int SAMPLED_T_ON = 8;
	private static final int SAMPLED_T_OFF = 9;
	private static final int NOISE = 10;
	private static final int SNR = 11;
	private static final int SNR_CONTINUOUS = 12;
	private static final int DENSITY = 13;
	private static final int PRECISION = 14;
	private static final int PRECISION_IN_FOCUS = 15;
	private static final int X = 16;
	private static final int Y = 17;
	private static final int Z = 18;
	private static final int WIDTH = 19;

	private static boolean[] integerDisplay;
	static
	{
		integerDisplay = new boolean[NAMES.length];
		integerDisplay[SAMPLES] = true;
		// Signal is an integer but there will be a large range
		integerDisplay[SIGNAL] = false;
		integerDisplay[SIGNAL_CONTINUOUS] = false;
		integerDisplay[TOTAL_SIGNAL] = false;
		integerDisplay[BLINKS] = true;
		integerDisplay[SAMPLED_BLINKS] = true;
		integerDisplay[SAMPLED_T_ON] = false;
		integerDisplay[SAMPLED_T_OFF] = false;
		integerDisplay[DENSITY] = true;
	}
	private static boolean[] alwaysRemoveOutliers;
	static
	{
		alwaysRemoveOutliers = new boolean[NAMES.length];
		alwaysRemoveOutliers[PRECISION] = true;
		alwaysRemoveOutliers[PRECISION_IN_FOCUS] = true;
	}

	private String resultsFileHeader = null;
	private AtomicInteger photonsRemoved;
	private AtomicInteger t1Removed;
	private AtomicInteger tNRemoved;
	private SummaryStatistics photonStats;
	private boolean imagePSF;
	private double hwhm = 0;

	private TIntHashSet movingMolecules;
	private boolean maskListContainsStacks;

	// Created by drawImage(...)
	private MemoryPeakResults results = null;

	// Used by the ImageGenerator to show progress when the thread starts
	private int frame, maxT, totalFrames;

	private XStream xs = null;
	private boolean simpleMode = false;
	private boolean benchmarkMode = false;
	private boolean spotMode = false;
	private boolean trackMode = false;
	private boolean extraOptions = false;

	// Hold private variables for settings that are ignored in simple/benchmark mode 
	private boolean poissonNoise = true;
	private double minPhotons = 0, minSNRt1 = 0, minSNRtN = 0;

	// Store the parameters for the last simulation for spot data
	public static class SimulationParameters
	{
		private static int nextId = 1;

		/**
		 * The parameter set identifier
		 */
		final int id;
		/**
		 * Number of molecules in the simulated image
		 */
		final int molecules;
		/**
		 * True if using a full simulation of fluorophores with a lifetime. False is for single random localisations per
		 * frame.
		 */
		final boolean fullSimulation;
		/**
		 * Gaussian standard deviation
		 */
		final double s;
		/**
		 * Pixel pitch in nm
		 */
		final double a;
		/**
		 * The min number of photons per frame
		 */
		final double minSignal;
		/**
		 * The max number of photons per frame
		 */
		final double maxSignal;
		/**
		 * The average signal per frame
		 */
		final double signalPerFrame;
		/**
		 * The z-position depth
		 */
		final double depth;
		/**
		 * True if the depth is fixed
		 */
		final boolean fixedDepth;
		/**
		 * The camera bias
		 */
		final double bias;
		/**
		 * True if EM-gain was modelled
		 */
		final boolean emCCD;
		/**
		 * Total gain (ADUs/photon)
		 */
		final double gain;
		/**
		 * Amplification gain (ADUs/electron)
		 */
		final double amplification;
		/**
		 * Read noise in ADUs
		 */
		final double readNoise;
		/**
		 * Background
		 */
		final double b;
		/**
		 * Background noise in photons per pixel (used in the precision calculations)
		 */
		final double b2;

		private boolean loaded;

		public SimulationParameters(int molecules, boolean fullSimulation, double s, double a, double minSignal,
				double maxSignal, double signalPerFrame, double depth, boolean fixedDepth, double bias, boolean emCCD,
				double gain, double amplification, double readNoise, double b, double b2)
		{
			id = nextId++;
			this.molecules = molecules;
			this.fullSimulation = fullSimulation;
			this.s = s;
			this.a = a;
			this.minSignal = minSignal;
			this.maxSignal = maxSignal;
			this.signalPerFrame = signalPerFrame;
			this.depth = depth;
			// We must have a fixed depth if the depth is zero
			this.fixedDepth = (depth > 0) ? fixedDepth : true;
			this.bias = bias;
			this.emCCD = emCCD;
			this.gain = gain;
			this.amplification = amplification;
			this.readNoise = readNoise;
			this.b = b;
			this.b2 = b2;
		}

		/**
		 * Checks if is a loaded simulation.
		 *
		 * @return true, if is loaded
		 */
		public boolean isLoaded()
		{
			return loaded;
		}
	}

	// Store the parameters for the last benchmark
	public static class BenchmarkParameters
	{
		private static int nextId = 1;

		/**
		 * The parameter set identifier
		 */
		final int id;
		/**
		 * Number of frames in the simulated image
		 */
		final int frames;
		/**
		 * Gaussian standard deviation
		 */
		final double s;
		/**
		 * Pixel pitch in nm
		 */
		final double a;
		/**
		 * The average number of photons per frame
		 */
		private double signal;
		/**
		 * The x,y,z position of the localisation in each frame
		 */
		final double x, y, z;
		final double bias;
		/**
		 * True if EM-gain was modelled
		 */
		final boolean emCCD;
		/**
		 * Total gain
		 */
		final double gain;
		/**
		 * Amplification gain (ADUs/electron)
		 */
		final double amplification;
		/**
		 * Read noise in ADUs
		 */
		final double readNoise;
		/**
		 * Background
		 */
		private double b;
		/**
		 * Background noise in photons per pixel (used in the precision calculations)
		 */
		final double b2;
		/**
		 * The actual number of simulated ADUs in each frame of the benchmark image. Some frames may be empty (due to
		 * signal filtering or Poisson sampling). The count is after gain has been applied to the photons.
		 */
		final double[] p;
		/**
		 * The actual number of simulated background ADUs in each frame of the benchmark image. The count is after gain
		 * has been applied to the photons.
		 */
		final double[] background;
		/**
		 * The number of frames with a simulated photon count
		 */
		private int molecules;

		final double precisionN, precisionX, precisionXML;

		public BenchmarkParameters(int frames, double s, double a, double signal, double x, double y, double z,
				double bias, boolean emCCD, double gain, double amplification, double readNoise, double b, double b2,
				double precisionN, double precisionX, double precisionXML)
		{
			id = nextId++;
			this.frames = frames;
			this.s = s;
			this.a = a;
			this.signal = signal;
			this.x = x;
			this.y = y;
			this.z = z;
			this.bias = bias;
			this.emCCD = emCCD;
			this.gain = gain;
			this.amplification = amplification;
			this.readNoise = readNoise;
			this.b = b;
			this.b2 = b2;
			this.precisionN = precisionN;
			this.precisionX = precisionX;
			this.precisionXML = precisionXML;
			p = new double[frames];
			background = new double[frames];
		}

		public void setPhotons(MemoryPeakResults results)
		{
			molecules = results.size();
			double sum = 0, sum2 = 0;
			for (PeakResult result : results)
			{
				int i = result.getFrame() - 1;
				if (p[i] != 0)
					throw new RuntimeException("Multiple peaks on the same frame: " + result.getFrame());
				p[i] = result.getSignal();
				background[i] = result.getBackground() - bias;
				sum += p[i];
				sum2 += background[i];
			}
			sum /= molecules;
			sum2 /= molecules;
			Utils.log(
					"Created %d frames, %d molecules. Simulated signal %s : average %s. Simulated background %s : average %s",
					frames, molecules, Utils.rounded(signal), Utils.rounded(sum / gain), Utils.rounded(b),
					Utils.rounded(sum2 / gain));
			// Reset the average signal and background (in photons)
			signal = sum / gain;
			b = sum2 / gain;
		}

		/**
		 * @return The average number of photons per frame
		 */
		public double getSignal()
		{
			return signal;
		}

		/**
		 * @return
		 * 		The number of frames with a simulated photon count
		 */
		public int getMolecules()
		{
			return molecules;
		}

		/**
		 * @return the average number of background photons per frame
		 */
		public double getBackground()
		{
			return b;
		}
	}

	static BenchmarkParameters benchmarkParameters = null;
	static SimulationParameters simulationParameters = null;

	private static String benchmarkFile = "";
	private static String benchmarkImage = "";
	private static boolean benchmarkAuto = true;
	private static int benchmarkImageId = 0;
	private static String benchmarkResultsName = "";

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		extraOptions = Utils.isExtraOptions();
		simpleMode = (arg != null && arg.contains("simple"));
		benchmarkMode = (arg != null && arg.contains("benchmark"));
		spotMode = (arg != null && arg.contains("spot"));
		trackMode = (arg != null && arg.contains("track"));

		if ("load".equals(arg))
		{
			loadBenchmarkData();
			return;
		}

		// Each localisation is a simulated emission of light from a point in space and time
		List<LocalisationModel> localisations = null;

		// Each localisation set is a collection of localisations that represent all localisations
		// with the same ID that are on in the same image time frame (Note: the simulation
		// can create many localisations per fluorophore per time frame which is useful when 
		// modelling moving particles)
		List<LocalisationModelSet> localisationSets = null;

		// Each fluorophore contains the on and off times when light was emitted 
		List<? extends FluorophoreSequenceModel> fluorophores = null;

		if (simpleMode || benchmarkMode || spotMode)
		{
			if (!showSimpleDialog())
				return;
			resetMemory();

			settings.exposureTime = 1000; // 1 second frames
			areaInUm = settings.size * settings.pixelPitch * settings.size * settings.pixelPitch / 1e6;

			// Number of spots per frame
			int n = 0;
			int[] nextN = null;
			SpatialDistribution dist;

			if (benchmarkMode)
			{
				// --------------------
				// BENCHMARK SIMULATION
				// --------------------
				// Draw the same point on the image repeatedly
				n = 1;
				dist = createFixedDistribution();

				reportAndSaveFittingLimits(dist);
			}
			else if (spotMode)
			{
				// ---------------
				// SPOT SIMULATION
				// ---------------
				// The spot simulation draws 0 or 1 random point per frame. 
				// Ensure we have 50% of the frames with a spot.
				nextN = new int[settings.particles * 2];
				Arrays.fill(nextN, 0, settings.particles, 1);
				Random rand = new Random();
				rand.shuffle(nextN);

				// Only put spots in the central part of the image
				double border = settings.size / 4.0;
				dist = createUniformDistribution(border);
			}
			else
			{
				// -----------------
				// SIMPLE SIMULATION
				// -----------------
				// The simple simulation draws n random points per frame to achieve a specified density.
				// No points will appear in multiple frames.
				// Each point has a random number of photons sampled from a range.

				// We can optionally use a mask. Create his first as it updates the areaInUm
				dist = createDistribution();

				// Randomly sample (i.e. not uniform density in all frames)
				if (settings.samplePerFrame)
				{
					final double mean = areaInUm * settings.density;
					System.out.printf("Mean samples = %f\n", mean);
					if (mean < 0.5)
					{
						GenericDialog gd = new GenericDialog(TITLE);
						gd.addMessage("The mean samples per frame is low: " + Utils.rounded(mean) + "\n \nContinue?");
						gd.enableYesNoCancel();
						gd.hideCancelButton();
						gd.showDialog();
						if (!gd.wasOKed())
							return;
					}
					PoissonDistribution poisson = new PoissonDistribution(createRandomGenerator(), mean,
							PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);
					StoredDataStatistics samples = new StoredDataStatistics(settings.particles);
					while (samples.getSum() < settings.particles)
					{
						samples.add(poisson.sample());
					}
					nextN = new int[samples.getN()];
					for (int i = 0; i < nextN.length; i++)
						nextN[i] = (int) samples.getValue(i);
				}
				else
				{
					// Use the density to get the number per frame
					n = (int) FastMath.max(1, Math.round(areaInUm * settings.density));
				}
			}

			RandomGenerator random = null;

			localisations = new ArrayList<LocalisationModel>(settings.particles);
			localisationSets = new ArrayList<LocalisationModelSet>(settings.particles);

			final int minPhotons = (int) settings.photonsPerSecond;
			final int range = (int) settings.photonsPerSecondMaximum - minPhotons + 1;
			if (range > 1)
				random = createRandomGenerator();

			// Add frames at the specified density until the number of particles has been reached
			int id = 0;
			int t = 0;
			while (id < settings.particles)
			{
				// Allow the number per frame to be specified
				if (nextN != null)
				{
					if (t >= nextN.length)
						break;
					n = nextN[t];
				}

				// Simulate random positions in the frame for the specified density
				t++;
				for (int j = 0; j < n; j++)
				{
					final double[] xyz = dist.next();

					// Ignore within border. We do not want to draw things we cannot fit.
					//if (!distBorder.isWithinXY(xyz))
					//	continue;

					// Simulate random photons
					final int intensity = minPhotons + ((random != null) ? random.nextInt(range) : 0);

					LocalisationModel m = new LocalisationModel(id, t, xyz, intensity, LocalisationModel.CONTINUOUS);
					localisations.add(m);

					// Each localisation can be a separate localisation set
					LocalisationModelSet set = new LocalisationModelSet(id, t);
					set.add(m);
					localisationSets.add(set);

					id++;
				}
			}
		}
		else
		{
			// This is used for the track mode as well as the full simulation.

			if (!showDialog())
				return;
			resetMemory();

			areaInUm = settings.size * settings.pixelPitch * settings.size * settings.pixelPitch / 1e6;

			int totalSteps;
			double correlation = 0;
			ImageModel imageModel;

			if (trackMode)
			{
				// ----------------
				// TRACK SIMULATION
				// ----------------
				// In track mode we create fixed lifetime fluorophores that do not overlap in time.
				// This is the simplest simulation to test moving molecules.
				settings.seconds = (int) Math.ceil(settings.particles * (settings.exposureTime + settings.tOn) / 1000);
				totalSteps = 0;

				final double simulationStepsPerFrame = (settings.stepsPerSecond * settings.exposureTime) / 1000.0;
				imageModel = new FixedLifetimeImageModel(settings.stepsPerSecond * settings.tOn / 1000.0,
						simulationStepsPerFrame);
			}
			else
			{
				// ---------------
				// FULL SIMULATION
				// ---------------
				// The full simulation draws n random points in space.
				// The same molecule may appear in multiple frames, move and blink.
				//
				// Points are modelled as fluorophores that must be activated and then will 
				// blink and photo-bleach. The molecules may diffuse and this can be simulated 
				// with many steps per image frame. All steps from a frame are collected
				// into a localisation set which can be drawn on the output image.

				SpatialIllumination activationIllumination = createIllumination(settings.pulseRatio,
						settings.pulseInterval);

				// Generate additional frames so that each frame has the set number of simulation steps
				totalSteps = (int) Math.ceil(settings.seconds * settings.stepsPerSecond);

				// Since we have an exponential decay of activations
				// ensure half of the particles have activated by 30% of the frames.
				double eAct = totalSteps * 0.3 * activationIllumination.getAveragePhotons();

				// Q. Does tOn/tOff change depending on the illumination strength?
				imageModel = new ActivationEnergyImageModel(eAct, activationIllumination,
						settings.stepsPerSecond * settings.tOn / 1000.0,
						settings.stepsPerSecond * settings.tOffShort / 1000.0,
						settings.stepsPerSecond * settings.tOffLong / 1000.0, settings.nBlinksShort,
						settings.nBlinksLong);
				imageModel.setUseGeometricDistribution(settings.nBlinksGeometricDistribution);

				// Only use the correlation if selected for the distribution
				if (PHOTON_DISTRIBUTION[PHOTON_CORRELATED].equals(settings.photonDistribution))
					correlation = settings.correlation;
			}

			imageModel.setRandomGenerator(createRandomGenerator());
			imageModel.setPhotonBudgetPerFrame(true);
			imageModel.setDiffusion2D(settings.diffuse2D);
			imageModel.setRotation2D(settings.rotate2D);

			IJ.showStatus("Creating molecules ...");
			SpatialDistribution distribution = createDistribution();
			List<CompoundMoleculeModel> compounds = createCompoundMolecules();
			if (compounds == null)
				return;
			List<CompoundMoleculeModel> molecules = imageModel.createMolecules(compounds, settings.particles,
					distribution, settings.rotateInitialOrientation);

			// Activate fluorophores
			IJ.showStatus("Creating fluorophores ...");

			// Note: molecules list will be converted to compounds containing fluorophores
			fluorophores = imageModel.createFluorophores(molecules, totalSteps);

			if (fluorophores.isEmpty())
			{
				IJ.error(TITLE, "No fluorophores created");
				return;
			}

			IJ.showStatus("Creating localisations ...");

			// TODO - Output a molecule Id for each fluorophore if using compound molecules. This allows analysis
			// of the ratio of trimers, dimers, monomers, etc that could be detected.

			totalSteps = checkTotalSteps(totalSteps, fluorophores);
			if (totalSteps == 0)
				return;

			imageModel.setPhotonDistribution(createPhotonDistribution());
			imageModel.setConfinementDistribution(createConfinementDistribution());
			// This should be optimised
			imageModel.setConfinementAttempts(10);

			localisations = imageModel.createImage(molecules, settings.fixedFraction, totalSteps,
					(double) settings.photonsPerSecond / settings.stepsPerSecond, correlation,
					settings.rotateDuringSimulation);

			// Re-adjust the fluorophores to the correct time
			if (settings.stepsPerSecond != 1)
			{
				final double scale = 1.0 / settings.stepsPerSecond;
				for (FluorophoreSequenceModel f : fluorophores)
					f.adjustTime(scale);
			}

			// Integrate the frames
			localisationSets = combineSimulationSteps(localisations);

			localisationSets = filterToImageBounds(localisationSets);
		}

		datasetNumber++;

		localisations = drawImage(localisationSets);

		if (localisations == null || localisations.isEmpty())
		{
			IJ.error(TITLE, "No localisations created");
			return;
		}

		fluorophores = removeFilteredFluorophores(fluorophores, localisations);

		double signalPerFrame = showSummary(fluorophores, localisations);

		if (!benchmarkMode)
		{
			boolean fullSimulation = (!(simpleMode || spotMode));
			saveSimulationParameters(localisations.size(), fullSimulation, signalPerFrame);
		}

		IJ.showStatus("Saving data ...");

		//convertRelativeToAbsolute(molecules);
		saveFluorophores(fluorophores);
		saveImageResults(results);
		saveLocalisations(localisations);

		// The settings for the filenames may have changed
		SettingsManager.saveSettings(globalSettings);

		IJ.showStatus("Done");
	}

	private void resetMemory()
	{
		benchmarkParameters = null;
		simulationParameters = null;
		setBenchmarkResults(null, null);
		// Run the garbage collector to free memory
		MemoryPeakResults.runGC();
	}

	/**
	 * Output the theoretical limits for fitting a Gaussian and store the benchmark settings
	 * 
	 * @param dist
	 *            The distribution
	 */
	private void reportAndSaveFittingLimits(SpatialDistribution dist)
	{
		final double totalGain = (settings.getTotalGain() > 0) ? settings.getTotalGain() : 1;

		// Background is in photons
		double backgroundVariance = settings.background;
		// Do not add EM-CCD noise factor. The Mortensen formula also includes this factor 
		// so this is "double-counting" the EM-CCD.  
		//if (settings.getEmGain() > 1)
		//	backgroundVariance *= 2;

		final double backgroundVarianceInADUs = settings.background * totalGain * totalGain *
				((settings.getEmGain() > 1) ? 2 : 1);

		// Read noise is in electrons. Convert to Photons
		double readNoise = settings.readNoise / totalGain;
		if (settings.getCameraGain() != 0)
			readNoise *= settings.getCameraGain();

		final double readVariance = readNoise * readNoise;

		double readVarianceInADUs = settings.readNoise *
				((settings.getCameraGain() != 0) ? settings.getCameraGain() : 1);
		readVarianceInADUs *= readVarianceInADUs;

		// Get the expected value at each pixel in photons. Assuming a Poisson distribution this 
		// is equal to the total variance at the pixel.
		final double b2 = backgroundVariance + readVariance;

		boolean emCCD = settings.getEmGain() > 1;
		double sd = getPsfSD() * settings.pixelPitch;

		// The precision calculation is dependent on the model. The classic Mortensen formula
		// is for a Gaussian Mask Estimator. Use other equation for MLE. The formula provided 
		// for WLSE requires an offset to the background used to stabilise the fitting. This is
		// not implemented (i.e. we used an offset of zero) and in this case the WLSE precision 
		// is the same as MLE with the caveat of numerical instability.

		double lowerP = PeakResult.getPrecisionX(settings.pixelPitch, sd, settings.photonsPerSecondMaximum, b2, emCCD);
		double upperP = PeakResult.getPrecisionX(settings.pixelPitch, sd, settings.photonsPerSecond, b2, emCCD);
		double lowerMLP = PeakResult.getMLPrecisionX(settings.pixelPitch, sd, settings.photonsPerSecondMaximum, b2,
				emCCD);
		double upperMLP = PeakResult.getMLPrecisionX(settings.pixelPitch, sd, settings.photonsPerSecond, b2, emCCD);
		double lowerN = getPrecisionN(settings.pixelPitch, sd, settings.photonsPerSecond, b2, emCCD);
		double upperN = getPrecisionN(settings.pixelPitch, sd, settings.photonsPerSecondMaximum, b2, emCCD);
		//final double b = Math.sqrt(b2);
		Utils.log(TITLE + " Benchmark");
		double[] xyz = dist.next().clone();
		double offset = settings.size * 0.5;
		for (int i = 0; i < 2; i++)
			xyz[i] += offset;
		Utils.log("X = %s nm : %s px", Utils.rounded(xyz[0] * settings.pixelPitch), Utils.rounded(xyz[0], 6));
		Utils.log("Y = %s nm : %s px", Utils.rounded(xyz[1] * settings.pixelPitch), Utils.rounded(xyz[1], 6));
		Utils.log("Width (s) = %s nm : %s px", Utils.rounded(sd), Utils.rounded(sd / settings.pixelPitch));
		final double sa = PSFCalculator.squarePixelAdjustment(sd, settings.pixelPitch);
		Utils.log("Adjusted Width (sa) = %s nm : %s px", Utils.rounded(sa), Utils.rounded(sa / settings.pixelPitch));
		Utils.log("Signal (N) = %s - %s photons : %s - %s ADUs", Utils.rounded(settings.photonsPerSecond),
				Utils.rounded(settings.photonsPerSecondMaximum), Utils.rounded(settings.photonsPerSecond * totalGain),
				Utils.rounded(settings.photonsPerSecondMaximum * totalGain));
		final double noiseInADUs = Math.sqrt(readVarianceInADUs + backgroundVarianceInADUs);
		Utils.log("Pixel noise = %s photons : %s ADUs", Utils.rounded(noiseInADUs / totalGain),
				Utils.rounded(noiseInADUs));
		Utils.log(
				"Expected background variance pre EM-gain (b^2) = %s photons^2 (%s ADUs^2) " +
						"[includes read variance converted to photons]",
				Utils.rounded(b2), Utils.rounded(b2 * totalGain * totalGain));
		Utils.log("Localisation precision (LSE): %s - %s nm : %s - %s px", Utils.rounded(lowerP), Utils.rounded(upperP),
				Utils.rounded(lowerP / settings.pixelPitch), Utils.rounded(upperP / settings.pixelPitch));
		Utils.log("Localisation precision (MLE): %s - %s nm : %s - %s px", Utils.rounded(lowerMLP),
				Utils.rounded(upperMLP), Utils.rounded(lowerMLP / settings.pixelPitch),
				Utils.rounded(upperMLP / settings.pixelPitch));
		Utils.log("Signal precision: %s - %s photons : %s - %s ADUs", Utils.rounded(lowerN), Utils.rounded(upperN),
				Utils.rounded(lowerN * totalGain), Utils.rounded(upperN * totalGain));

		// Store the benchmark settings when not using variable photons
		if (settings.photonsPerSecond == settings.photonsPerSecondMaximum)
		{
			final double amplification = totalGain /
					((settings.getQuantumEfficiency() == 0) ? 1 : settings.getQuantumEfficiency());
			// Store read noise in ADUs
			readNoise = settings.readNoise * ((settings.getCameraGain() > 0) ? settings.getCameraGain() : 1);
			benchmarkParameters = new BenchmarkParameters(settings.particles, sd, settings.pixelPitch,
					settings.photonsPerSecond, xyz[0], xyz[1], xyz[2], settings.bias, emCCD, totalGain, amplification,
					readNoise, settings.background, b2, lowerN, lowerP, lowerMLP);
		}
		else
		{
			Utils.log(
					"Warning: Benchmark settings are only stored in memory when the number of photons is fixed. Min %s != Max %s",
					Utils.rounded(settings.photonsPerSecond), Utils.rounded(settings.photonsPerSecondMaximum));
		}
	}

	/**
	 * Store the simulation settings
	 * 
	 * @param particles
	 */
	private void saveSimulationParameters(int particles, boolean fullSimulation, double signalPerFrame)
	{
		final double totalGain = (settings.getTotalGain() > 0) ? settings.getTotalGain() : 1;

		// Background is in photons
		double backgroundVariance = settings.background;
		// Do not add EM-CCD noise factor. The Mortensen formula also includes this factor 
		// so this is "double-counting" the EM-CCD.  
		//if (settings.getEmGain() > 1)
		//	backgroundVariance *= 2;

		// Read noise is in electrons. Convert to Photons to get contribution to background variance
		double readNoise = settings.readNoise / totalGain;
		if (settings.getCameraGain() != 0)
			readNoise *= settings.getCameraGain();

		final double readVariance = readNoise * readNoise;

		// Get the expected value at each pixel in photons. Assuming a Poisson distribution this 
		// is equal to the total variance at the pixel.
		final double b2 = backgroundVariance + readVariance;

		boolean emCCD = settings.getEmGain() > 1;
		double sd = getPsfSD() * settings.pixelPitch;

		final double amplification = totalGain /
				((settings.getQuantumEfficiency() == 0) ? 1 : settings.getQuantumEfficiency());

		// Store read noise in ADUs
		readNoise = settings.readNoise * ((settings.getCameraGain() > 0) ? settings.getCameraGain() : 1);
		simulationParameters = new SimulationParameters(particles, fullSimulation, sd, settings.pixelPitch,
				settings.photonsPerSecond, settings.photonsPerSecondMaximum, signalPerFrame, settings.depth,
				settings.fixedDepth, settings.bias, emCCD, totalGain, amplification, readNoise, settings.background,
				b2);
	}

	/**
	 * Calculate the signal precision for least squares fitting. Uses the Thompson formula:
	 * (Thompson, et al (2002) Biophysical Journal 82, 2775-2783), equation 19
	 * 
	 * @param a
	 *            The size of the pixels in nm
	 * @param s
	 *            The peak standard deviation in nm
	 * @param N
	 *            The peak signal in photons
	 * @param b2
	 *            The expected number of photons per pixel from a background with spatially constant
	 *            expectation value across the image (Note that this is b^2 not b, which could be the standard deviation
	 *            of the image pixels)
	 * @param emCCD
	 *            True if an emCCD camera
	 * @return The signal precision in photons
	 */
	public static double getPrecisionN(double a, double s, double N, double b2, boolean emCCD)
	{
		// EM-CCD noise factor
		final double F = (emCCD) ? 2 : 1;
		final double a2 = a * a;
		// 4 * pi = 12.56637061

		// Adjustment for square pixels
		//final double sa2 = s * s + a2 / 12.0;

		// Original Thompson formula modified for EM-gain noise factor.

		// TODO - Investigate if this limit is correct

		// My fitters approach this limit when background is 0 photon and EM-gain = 0.
		// The fitters are above this limit when background is >0 photon and EM-gain = 0.

		// The MLE fitter can approach this limit when background is 0 photon and EM-gain = 25.

		return Math.sqrt(F * (N + (12.56637061 * s * s * b2) / a2));
		//return Math.sqrt(F * (N + (12.56637061 * sa2 * b2) / a2));
	}

	/**
	 * Check if the total steps can fit all the fluorophores end times. If not then ask the user if they want to draw
	 * extra
	 * frames. Return the total steps to simulate (either the original steps or a larger number to fit all the data).
	 * 
	 * @param totalSteps
	 * @param fluorophores
	 * @return The new total steps to simulate
	 */
	private int checkTotalSteps(int totalSteps, List<? extends FluorophoreSequenceModel> fluorophores)
	{
		int max = totalSteps;
		for (FluorophoreSequenceModel f : fluorophores)
		{
			if (max < f.getEndTime())
				max = (int) (f.getEndTime() + 1);
		}
		if (max > totalSteps)
		{
			GenericDialog gd = new GenericDialog(TITLE);
			gd.enableYesNoCancel();
			gd.hideCancelButton();
			final double simulationStepsPerFrame = (settings.stepsPerSecond * settings.exposureTime) / 1000.0;
			int newFrames = 1 + (int) (max / simulationStepsPerFrame);

			if (totalSteps != 0)
			{
				int totalFrames = (int) Math.ceil(settings.seconds * 1000 / settings.exposureTime);
				gd.addMessage(String.format(
						"Require %d (%s%%) additional frames to draw all fluorophores.\nDo you want to add extra frames?",
						newFrames - totalFrames, Utils.rounded((100.0 * (newFrames - totalFrames)) / totalFrames, 3)));
			}
			else
			{
				gd.addMessage(String.format("Require %d frames to draw all fluorophores.\nDo you want to proceed?",
						newFrames));
			}
			gd.showDialog();
			if (gd.wasOKed())
				totalSteps = max;
		}
		return totalSteps;
	}

	private SpatialDistribution createDistribution()
	{
		if (settings.distribution.equals(DISTRIBUTION[MASK]))
		{
			ImagePlus imp = WindowManager.getImage(settings.distributionMask);
			if (imp != null)
			{
				return createMaskDistribution(imp, settings.distributionMaskSliceDepth, true);
			}
		}
		else if (settings.distribution.equals(DISTRIBUTION[GRID]))
		{
			return new GridDistribution(settings.size, settings.depth / settings.pixelPitch, settings.cellSize,
					settings.probabilityBinary, settings.minBinaryDistance / settings.pixelPitch,
					settings.maxBinaryDistance / settings.pixelPitch);
		}

		return createUniformDistributionWithPSFWidthBorder();
	}

	private SpatialDistribution createMaskDistribution(ImagePlus imp, double sliceDepth, boolean updateArea)
	{
		// Calculate the scale of the mask
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final double scaleX = (double) settings.size / w;
		final double scaleY = (double) settings.size / h;

		// Use an image for the distribution
		if (imp.getStackSize() > 1)
		{
			ImageStack stack = imp.getImageStack();
			List<int[]> masks = new ArrayList<int[]>(stack.getSize());
			int[] maxMask = new int[w * h];
			for (int slice = 1; slice <= stack.getSize(); slice++)
			{
				int[] mask = extractMask(stack.getProcessor(slice));
				if (updateArea)
				{
					for (int i = 0; i < mask.length; i++)
						if (mask[i] != 0)
							maxMask[i] = 1;
				}
				masks.add(mask);
			}
			if (updateArea)
				updateArea(maxMask, w, h);

			if (sliceDepth == 0)
				// Auto configure to the full depth of the simulation
				sliceDepth = settings.depth / masks.size();

			return new MaskDistribution3D(masks, w, h, sliceDepth / settings.pixelPitch, scaleX, scaleY,
					createRandomGenerator());
		}
		else
		{
			int[] mask = extractMask(imp.getProcessor());
			if (updateArea)
				updateArea(mask, w, h);
			return new MaskDistribution(mask, w, h, settings.depth / settings.pixelPitch, scaleX, scaleY,
					createRandomGenerator());
		}
	}

	private void updateArea(int[] mask, int w, int h)
	{
		// Note: The apparent area will be bigger due to the PSF width blurring the edges.
		// Assume the entire max intensity mask is in focus and the PSF is Gaussian in the focal plane.
		// Blur the mask
		float[] pixels = new float[mask.length];
		for (int i = 0; i < pixels.length; i++)
			pixels[i] = mask[i];

		GaussianFilter blur = new GaussianFilter();
		final double scaleX = (double) settings.size / w;
		final double scaleY = (double) settings.size / h;
		double extra = 1; // Allow extra?
		double sd = getPsfSD() * extra;
		blur.convolve(pixels, w, h, sd / scaleX, sd / scaleY);

		// Count pixels in blurred mask. Ignore those that are very faint (at the edge of the region)
		int c = 0;

		//		// By fraction of max value
		//		float limit = 0.1f;
		//		//float min = (float) (Maths.max(pixels) * limit);
		//		float min = limit;
		//		for (float f : pixels)
		//			if (f > min)
		//				c++;

		//		// Rank in order and get fraction of sum
		//		Arrays.sort(pixels);
		//		double sum = 0;
		//		double stop = Maths.sum(pixels) * 0.95;
		//		float after = Maths.max(pixels) + 1;
		//		for (float f : pixels)
		//		{
		//			sum += f;
		//			if (sum > stop)
		//			{
		//				break;
		//				//after = f;
		//				//stop = Float.POSITIVE_INFINITY;
		//			}
		//			if (f > after)
		//				break;
		//			c++;
		//		}

		// Threshold to a mask
		FloatProcessor fp = new FloatProcessor(w, h, pixels);
		ShortProcessor sp = (ShortProcessor) fp.convertToShort(true);
		final int t = AutoThreshold.getThreshold(AutoThreshold.Method.OTSU, sp.getHistogram());
		//Utils.display("Blurred", fp);
		for (int i = 0; i < mask.length; i++)
			if (sp.get(i) >= t)
				c++;

		// Convert 
		final double scale = ((double) c) / mask.length;
		//System.out.printf("Scale = %f\n", scale);
		areaInUm = scale * settings.size * settings.pixelPitch * settings.size * settings.pixelPitch / 1e6;
	}

	private int[] extractMask(ImageProcessor ip)
	{
		//ip = ip.duplicate();
		//ip.setInterpolationMethod(ImageProcessor.BILINEAR);
		//ip = ip.resize(settings.size, settings.size);
		int[] mask = new int[ip.getPixelCount()];
		for (int i = 0; i < mask.length; i++)
		{
			mask[i] = ip.get(i);
		}
		return mask;
	}

	private UniformDistribution createUniformDistributionWithPSFWidthBorder()
	{
		double border = getHWHM() * 3;
		border = FastMath.min(border, settings.size / 4);
		return createUniformDistribution(border);
	}

	private SpatialDistribution createFixedDistribution()
	{
		SpatialDistribution dist;
		dist = new SpatialDistribution()
		{
			private double[] xyz = new double[] { settings.xPosition / settings.pixelPitch,
					settings.yPosition / settings.pixelPitch, settings.zPosition / settings.pixelPitch };

			public double[] next()
			{
				return xyz;
			}

			public boolean isWithinXY(double[] xyz)
			{
				return true;
			}

			public boolean isWithin(double[] xyz)
			{
				return true;
			}

			public void initialise(double[] xyz)
			{
			}
		};
		return dist;
	}

	/**
	 * Get the PSF half-width at half-maxima
	 * 
	 * @return
	 */
	private double getHWHM()
	{
		if (hwhm == 0)
		{
			if (imagePSF)
			{
				hwhm = getImageHWHM();
			}
			else
			{
				final double sd = (settings.enterWidth) ? settings.psfSD
						: PSFCalculator.calculateStdDev(settings.wavelength, settings.numericalAperture);

				hwhm = Gaussian2DFunction.SD_TO_HWHM_FACTOR * sd / settings.pixelPitch;
			}
		}
		return hwhm;
	}

	/**
	 * Get the PSF standard deviation for a Gaussian using the PSF half-width at half-maxima
	 * 
	 * @return
	 */
	private double getPsfSD()
	{
		return getHWHM() / Gaussian2DFunction.SD_TO_HWHM_FACTOR;
	}

	/**
	 * Get the PSF half-width at half-maxima from the Image PSF
	 * 
	 * @return
	 */
	private double getImageHWHM()
	{
		ImagePlus imp = WindowManager.getImage(settings.psfImageName);
		if (imp == null)
		{
			IJ.error(TITLE, "Unable to create the PSF model from image: " + settings.psfImageName);
			return -1;
		}
		Object o = XmlUtils.fromXML(imp.getProperty("Info").toString());
		if (!(o != null && o instanceof PSFSettings))
		{
			IJ.error(TITLE, "Unknown PSF settings for image: " + imp.getTitle());
			return -1;
		}
		PSFSettings psfSettings = (PSFSettings) o;
		if (psfSettings.fwhm <= 0)
		{
			IJ.error(TITLE, "Unknown PSF FWHM setting for image: " + imp.getTitle());
			return -1;
		}
		if (psfSettings.nmPerPixel <= 0)
		{
			IJ.error(TITLE, "Unknown PSF nm/pixel setting for image: " + imp.getTitle());
			return -1;
		}

		// The width of the PSF is specified in pixels of the PSF image. Convert to the pixels of the 
		// output image
		return 0.5 * psfSettings.fwhm * psfSettings.nmPerPixel / settings.pixelPitch;
	}

	/**
	 * Create distribution within an XY border
	 * 
	 * @param border
	 * @return
	 */
	private UniformDistribution createUniformDistribution(double border)
	{
		double depth = (settings.fixedDepth) ? settings.depth / settings.pixelPitch
				: settings.depth / (2 * settings.pixelPitch);

		// Ensure the focal plane is in the middle of the zDepth
		double[] max = new double[] { settings.size / 2 - border, settings.size / 2 - border, depth };
		double[] min = new double[3];
		for (int i = 0; i < 3; i++)
			min[i] = -max[i];
		if (settings.fixedDepth)
			min[2] = max[2];

		// Try using different distributions:
		final RandomGenerator rand1 = createRandomGenerator();

		if (settings.distribution.equals(DISTRIBUTION[UNIFORM_HALTON]))
		{
			return new UniformDistribution(min, max, rand1.nextInt());
		}

		if (settings.distribution.equals(DISTRIBUTION[UNIFORM_SOBOL]))
		{
			SobolSequenceGenerator rvg = new SobolSequenceGenerator(3);
			rvg.skipTo(rand1.nextInt());
			return new UniformDistribution(min, max, rvg);
		}

		// Create a distribution using random generators for each dimension 
		UniformDistribution distribution = new UniformDistribution(min, max, this);
		return distribution;
	}

	private SpatialDistribution createConfinementDistribution()
	{
		if (settings.diffusionRate <= 0 || settings.fixedFraction >= 1)
			return null;

		if (settings.confinement.equals(CONFINEMENT[CONFINEMENT_MASK]))
		{
			ImagePlus imp = WindowManager.getImage(settings.confinementMask);
			if (imp != null)
			{
				return createMaskDistribution(imp, settings.confinementMaskSliceDepth, false);
			}
		}
		else if (settings.confinement.equals(CONFINEMENT[CONFINEMENT_SPHERE]))
		{
			return new SphericalDistribution(settings.confinementRadius / settings.pixelPitch);
		}
		else if (settings.confinement.equals(CONFINEMENT[CONFINEMENT_WITHIN_IMAGE]))
		{
			//return createUniformDistribution(0);
			return createUniformDistributionWithPSFWidthBorder();
		}

		return null;
	}

	private SpatialIllumination createIllumination(double intensity, int pulseInterval)
	{
		if (settings.illumination.equals(ILLUMINATION[RADIAL]))
		{
			if (pulseInterval > 1)
				return new RadialFalloffIllumination(1, settings.size / 2, intensity, pulseInterval);
			return new RadialFalloffIllumination(intensity, settings.size / 2);
		}
		else
		{
			if (pulseInterval > 1)
				return new UniformIllumination(1, intensity, pulseInterval);
			uniformBackground = true;
			return new UniformIllumination(intensity);
		}
	}

	/**
	 * Filter those not in the distribution
	 * 
	 * @param localisationSets
	 * @return
	 */
	private List<LocalisationModelSet> filterToImageBounds(List<LocalisationModelSet> localisationSets)
	{
		List<LocalisationModelSet> newLocalisations = new ArrayList<LocalisationModelSet>(localisationSets.size());
		SpatialDistribution bounds = createUniformDistribution(0);
		for (LocalisationModelSet s : localisationSets)
		{
			if (bounds.isWithinXY(s.toLocalisation().getCoordinates()))
				newLocalisations.add(s);
		}
		return newLocalisations;
	}

	/**
	 * @return A photon distribution loaded from a file of floating-point values with the specified population mean.
	 */
	private RealDistribution createPhotonDistribution()
	{
		if (PHOTON_DISTRIBUTION[PHOTON_CUSTOM].equals(settings.photonDistribution))
		{
			// Get the distribution file
			String filename = Utils.getFilename("Photon_distribution", settings.photonDistributionFile);
			if (filename != null)
			{
				settings.photonDistributionFile = filename;
				try
				{
					InputStream is = new FileInputStream(new File(settings.photonDistributionFile));
					BufferedReader in = new BufferedReader(new UnicodeReader(is, null));
					StoredDataStatistics stats = new StoredDataStatistics();
					try
					{
						String str = null;
						double val = 0.0d;
						while ((str = in.readLine()) != null)
						{
							val = Double.parseDouble(str);
							stats.add(val);
						}
					}
					finally
					{
						in.close();
					}

					if (stats.getSum() > 0)
					{
						// Update the statistics to the desired mean.
						double scale = (double) settings.photonsPerSecond / stats.getMean();
						double[] values = stats.getValues();
						for (int i = 0; i < values.length; i++)
							values[i] *= scale;

						// TODO - Investigate the limits of this distribution. 
						// How far above and below the input data will values be generated.

						// Create the distribution using the recommended number of bins
						final int binCount = stats.getN() / 10;
						EmpiricalDistribution dist = new EmpiricalDistribution(binCount, createRandomGenerator());
						dist.load(values);
						return dist;
					}
				}
				catch (IOException e)
				{
					// Ignore
				}
				catch (NullArgumentException e)
				{
					// Ignore 
				}
				catch (NumberFormatException e)
				{
					// Ignore
				}
			}
			Utils.log("Failed to load custom photon distribution from file: %s. Default to fixed.",
					settings.photonDistributionFile);
		}
		else if (PHOTON_DISTRIBUTION[PHOTON_UNIFORM].equals(settings.photonDistribution))
		{
			if (settings.photonsPerSecond < settings.photonsPerSecondMaximum)
			{
				UniformRealDistribution dist = new UniformRealDistribution(createRandomGenerator(),
						settings.photonsPerSecond, settings.photonsPerSecondMaximum);
				return dist;
			}
		}
		else if (PHOTON_DISTRIBUTION[PHOTON_GAMMA].equals(settings.photonDistribution))
		{
			final double scaleParameter = settings.photonsPerSecond / settings.photonShape;
			GammaDistribution dist = new GammaDistribution(createRandomGenerator(), settings.photonShape,
					scaleParameter, ExponentialDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);
			return dist;
		}
		else if (PHOTON_DISTRIBUTION[PHOTON_CORRELATED].equals(settings.photonDistribution))
		{
			// No distribution required
			return null;
		}

		settings.photonDistribution = PHOTON_DISTRIBUTION[PHOTON_FIXED];
		return null;
	}

	private List<LocalisationModelSet> combineSimulationSteps(List<LocalisationModel> localisations)
	{
		// Allow fractional integration steps
		final double simulationStepsPerFrame = (settings.stepsPerSecond * settings.exposureTime) / 1000.0;

		List<LocalisationModelSet> newLocalisations = new ArrayList<LocalisationModelSet>(
				(int) (localisations.size() / simulationStepsPerFrame));

		//System.out.printf("combineSimulationSteps @ %f\n", simulationStepsPerFrame);

		final double gain = (settings.getTotalGain() > 0) ? settings.getTotalGain() : 1;
		sortLocalisationsByIdThenTime(localisations);
		int[] idList = getIds(localisations);
		movingMolecules = new TIntHashSet(idList.length);
		int index = 0;
		for (int id : idList)
		{
			int fromIndex = findIndexById(localisations, index, id);
			if (fromIndex > -1)
			{
				int toIndex = findLastIndexById(localisations, fromIndex, id);
				List<LocalisationModel> subset = localisations.subList(fromIndex, toIndex + 1);
				index = toIndex;

				// Store the IDs of any moving molecules
				if (isMoving(subset))
					movingMolecules.add(id);

				// The frames may be longer or shorter than the simulation steps. Allocate the step
				// proportionately to each frame it overlaps:
				//
				// Steps:  |-- 0 --|-- 1 --|-- 2 --|--
				// Frames: |--- 0 ---|--- 1 ---|--- 2 ---|
				//
				//         ^       ^
				//         |       |
				//         |       End frame
				//         |
				//         Start frame

				final double firstFrame = getStartFrame(subset.get(0), simulationStepsPerFrame);
				final double lastFrame = getEndFrame(subset.get(subset.size() - 1), simulationStepsPerFrame);

				// Get the first frame offset and allocate space to store all potential frames  
				final int intFirstFrame = (int) firstFrame;
				final int intLastFrame = (int) Math.ceil(lastFrame);
				LocalisationModelSet[] sets = new LocalisationModelSet[intLastFrame - intFirstFrame + 1];

				// Process each step
				for (LocalisationModel l : subset)
				{
					// Get the fractional start and end frames 
					double startFrame = getStartFrame(l, simulationStepsPerFrame);
					double endFrame = getEndFrame(l, simulationStepsPerFrame);

					// Round down to get the actual frames that are overlapped
					int start = (int) startFrame;
					int end = (int) endFrame;

					// Check if the span covers a fraction of the end frame, otherwise decrement to ignore that frame
					if (end > start && endFrame == end)
					{
						// E.g. convert 
						// Steps:      |-- 0 --|
						// Frames: |- 0 -|- 1 -|- 2 -|
						// to
						// Steps:      |-- 0 --|
						// Frames: |- 0 -|- 1 -|
						end--;
					}

					if (start == end)
					{
						// If the step falls within one frame then add it to the set
						int tIndex = start - intFirstFrame;
						if (sets[tIndex] == null)
							sets[tIndex] = new LocalisationModelSet(id, start);
						sets[tIndex].add(l);
					}
					else
					{
						// Add the localisation to all the frames that the step spans
						final double total = endFrame - startFrame;
						//double t = 0;
						for (int frame = start; frame <= end; frame++)
						{
							// Get the fraction to allocate to this frame
							double fraction;
							int state = (l.isContinuous() ? LocalisationModel.CONTINUOUS : 0);
							if (frame == start)
							{
								state |= LocalisationModel.NEXT | (l.hasPrevious() ? LocalisationModel.PREVIOUS : 0);
								// |-----|====|
								// |     |    ceil(startFrame)
								// |     startFrame
								// start
								if (startFrame == start)
									fraction = 1;
								else
									fraction = (Math.ceil(startFrame) - startFrame);
							}
							else if (frame == end)
							{
								state |= LocalisationModel.PREVIOUS | (l.hasNext() ? LocalisationModel.NEXT : 0);
								// |=====|----|
								// |     |     
								// |     endFrame
								// end
								fraction = (endFrame - end);
							}
							else
							{
								state |= LocalisationModel.CONTINUOUS;
								fraction = 1;
							}

							//t += fraction;

							// Add to the set
							int tIndex = frame - intFirstFrame;
							if (sets[tIndex] == null)
								sets[tIndex] = new LocalisationModelSet(id, frame);
							sets[tIndex].add(getFraction(l, fraction / total, state));
						}
						//if (t < total * 0.98)
						//{
						//	System.out.printf("Total error %g < %g : %f (%d) -> %f (%d)\n", t, total, startFrame,
						//			start, endFrame, end);
						//}
					}
				}

				LocalisationModelSet previous = null;
				for (int i = 0; i < sets.length; i++)
				{
					if (sets[i] != null)
					{
						sets[i].setPrevious(previous);

						// Create a data array and store the current intensity after gain. 
						// This is used later to filter based on SNR
						sets[i].setData(new double[] { 0, 0, 0, 0, sets[i].getIntensity() * gain });

						newLocalisations.add(sets[i]);
					}
					previous = sets[i];
				}
			}
		}
		// Sort by time
		Collections.sort(newLocalisations);
		return newLocalisations;
	}

	/**
	 * Check if any of the coordinates for the subset are different
	 * 
	 * @param subset
	 * @return True if the coordinates move
	 */
	private boolean isMoving(List<LocalisationModel> subset)
	{
		if (subset.size() < 2)
			return false;
		final double[] xyz = subset.get(0).getCoordinates();
		for (int i = 1; i < subset.size(); i++)
		{
			double[] xyz2 = subset.get(i).getCoordinates();
			for (int j = 0; j < 3; j++)
				if (xyz[j] != xyz2[j])
					return true;
		}
		return false;
	}

	/**
	 * Get the simulation frame start point for the localisation
	 * 
	 * @param localisationModel
	 * @param simulationStepsPerFrame
	 * @return
	 */
	private double getStartFrame(LocalisationModel localisationModel, double simulationStepsPerFrame)
	{
		// Time is 1-based (not 0)
		return (localisationModel.getTime() - 1) / simulationStepsPerFrame + 1;
	}

	/**
	 * Get the simulation frame end point for the localisation
	 * 
	 * @param localisationModel
	 * @param simulationStepsPerFrame
	 * @return
	 */
	private double getEndFrame(LocalisationModel localisationModel, double simulationStepsPerFrame)
	{
		// Time is 1-based (not 0)
		return (localisationModel.getTime()) / simulationStepsPerFrame + 1;
	}

	/**
	 * Create a new localisation model with the same id and position but with a fraction of the intensity and the
	 * specified state
	 * 
	 * @param l
	 * @param fraction
	 * @param state
	 * @return
	 */
	private LocalisationModel getFraction(LocalisationModel l, double fraction, int state)
	{
		return new LocalisationModel(l.getId(), l.getTime(), l.getCoordinates(), l.getIntensity() * fraction, state);
	}

	private int[] getIds(List<LocalisationModel> localisations)
	{
		int[] ids = new int[localisations.size()];
		if (localisations.isEmpty())
			return ids;
		int count = 0;
		int id = localisations.get(0).getId();
		ids[count++] = id;
		for (LocalisationModel l : localisations)
		{
			if (id != l.getId())
			{
				id = l.getId();
				ids[count++] = id;
			}
		}
		return Arrays.copyOf(ids, count);
	}

	//StoredDataStatistics rawPhotons = new StoredDataStatistics();
	//StoredDataStatistics drawPhotons = new StoredDataStatistics();

	//	private synchronized void addRaw(double d)
	//	{
	//		//rawPhotons.add(d);
	//	}
	//
	//	private synchronized void addDraw(double d)
	//	{
	//		//drawPhotons.add(d);
	//	}

	/**
	 * Create an image from the localisations using the configured PSF width. Draws a new stack
	 * image.
	 * <p>
	 * Note that the localisations are filtered using the signal. The input list of localisations will be updated.
	 * 
	 * @param localisationSets
	 * @return The localisations
	 */
	private List<LocalisationModel> drawImage(final List<LocalisationModelSet> localisationSets)
	{
		if (localisationSets.isEmpty())
			return null;

		// Create a new list for all localisation that are drawn (i.e. pass the signal filters)
		List<LocalisationModelSet> newLocalisations = Collections
				.synchronizedList(new ArrayList<LocalisationModelSet>(localisationSets.size()));
		photonsRemoved = new AtomicInteger();
		t1Removed = new AtomicInteger();
		tNRemoved = new AtomicInteger();
		photonStats = new SummaryStatistics();

		// Add drawn spots to memory
		results = new MemoryPeakResults();
		Calibration c = new Calibration(settings.pixelPitch, settings.getTotalGain(), settings.exposureTime);
		c.setEmCCD((settings.getEmGain() > 1));
		c.setBias(settings.bias);
		c.setReadNoise(settings.readNoise * ((settings.getCameraGain() > 0) ? settings.getCameraGain() : 1));
		c.setAmplification(settings.getAmplification());
		results.setCalibration(c);
		results.setSortAfterEnd(true);
		results.begin();

		maxT = localisationSets.get(localisationSets.size() - 1).getTime();

		// Display image
		ImageStack stack = new ImageStack(settings.size, settings.size, maxT);

		final double psfSD = getPsfSD();
		if (psfSD <= 0)
			return null;
		ImagePSFModel imagePSFModel = null;

		if (imagePSF)
		{
			// Create one Image PSF model that can be copied
			imagePSFModel = createImagePSF(localisationSets);
			if (imagePSFModel == null)
				return null;
		}

		IJ.showStatus("Drawing image ...");

		// Multi-thread for speed
		// Note that the default Executors.newCachedThreadPool() will continue to make threads if
		// new tasks are added. We need to limit the tasks that can be added using a fixed size
		// blocking queue.
		// http://stackoverflow.com/questions/1800317/impossible-to-make-a-cached-thread-pool-with-a-size-limit
		// ExecutorService threadPool = Executors.newCachedThreadPool();
		ExecutorService threadPool = Executors.newFixedThreadPool(Prefs.getThreads());
		List<Future<?>> futures = new LinkedList<Future<?>>();

		// Count all the frames to process
		frame = 0;
		totalFrames = maxT;

		// Collect statistics on the number of photons actually simulated

		// Process all frames
		int i = 0;
		int lastT = -1;
		for (LocalisationModelSet l : localisationSets)
		{
			if (Utils.isInterrupted())
				break;
			if (l.getTime() != lastT)
			{
				lastT = l.getTime();
				futures.add(threadPool.submit(
						new ImageGenerator(localisationSets, newLocalisations, i, lastT, createPSFModel(imagePSFModel),
								results, stack, poissonNoise, new RandomDataGenerator(createRandomGenerator()))));
			}
			i++;
		}
		// Finish processing data
		Utils.waitForCompletion(futures);
		futures.clear();
		if (Utils.isInterrupted())
		{
			IJ.showProgress(1);
			return null;
		}

		// Do all the frames that had no localisations
		for (int t = 1; t <= maxT; t++)
		{
			if (Utils.isInterrupted())
				break;
			if (stack.getPixels(t) == null)
			{
				futures.add(threadPool.submit(new ImageGenerator(localisationSets, newLocalisations, maxT, t, null,
						results, stack, poissonNoise, new RandomDataGenerator(createRandomGenerator()))));
			}
		}

		// Finish
		Utils.waitForCompletion(futures);
		threadPool.shutdown();
		IJ.showProgress(1);
		if (Utils.isInterrupted())
		{
			return null;
		}
		results.end();

		// Clear memory
		imagePSFModel = null;
		threadPool = null;
		futures.clear();
		futures = null;

		if (photonsRemoved.get() > 0)
			Utils.log("Removed %d localisations with less than %.1f rendered photons", photonsRemoved.get(),
					settings.minPhotons);
		if (t1Removed.get() > 0)
			Utils.log("Removed %d localisations with no neighbours @ SNR %.2f", t1Removed.get(), settings.minSNRt1);
		if (tNRemoved.get() > 0)
			Utils.log("Removed %d localisations with valid neighbours @ SNR %.2f", tNRemoved.get(), settings.minSNRtN);
		if (photonStats.getN() > 0)
			Utils.log("Average photons rendered = %s +/- %s", Utils.rounded(photonStats.getMean()),
					Utils.rounded(photonStats.getStandardDeviation()));

		//System.out.printf("rawPhotons = %f\n", rawPhotons.getMean());
		//System.out.printf("drawPhotons = %f\n", drawPhotons.getMean());
		//Utils.showHistogram("draw photons", drawPhotons, "photons", true, 0, 1000);

		// Update with all those localisation that have been drawn
		localisationSets.clear();
		localisationSets.addAll(newLocalisations);
		newLocalisations = null;

		IJ.showStatus("Displaying image ...");

		ImageStack newStack = stack;

		if (!settings.rawImage)
		{
			// Get the global limits and ensure all values can be represented
			Object[] imageArray = stack.getImageArray();
			float[] limits = Maths.limits((float[]) imageArray[0]);
			for (int j = 1; j < imageArray.length; j++)
				limits = Maths.limits(limits, (float[]) imageArray[j]);
			limits[0] = 0; // Leave bias in place
			// Check if the image will fit in a 16-bit range
			if ((limits[1] - limits[0]) < 65535)
			{
				// Convert to 16-bit
				newStack = new ImageStack(stack.getWidth(), stack.getHeight(), stack.getSize());
				// Account for rounding
				final float min = (float) (limits[0] - 0.5);
				for (int j = 0; j < imageArray.length; j++)
				{
					float[] image = (float[]) imageArray[j];
					short[] pixels = new short[image.length];
					for (int k = 0; k < pixels.length; k++)
					{
						pixels[k] = (short) (image[k] - min);
					}
					newStack.setPixels(pixels, j + 1);
					// Free memory
					imageArray[j] = null;
					// Attempt to stay within memory (check vs 32MB)
					if (MemoryPeakResults.freeMemory() < 33554432L)
						MemoryPeakResults.runGCOnce();
				}
			}
			else
			{
				// Keep as 32-bit but round to whole numbers
				for (int j = 0; j < imageArray.length; j++)
				{
					float[] pixels = (float[]) imageArray[j];
					for (int k = 0; k < pixels.length; k++)
					{
						pixels[k] = Math.round(pixels[k]);
					}
				}
			}
		}

		// Show image
		ImagePlus imp = Utils.display(CREATE_DATA_IMAGE_TITLE, newStack);

		ij.measure.Calibration cal = new ij.measure.Calibration();
		String unit = "nm";
		double unitPerPixel = settings.pixelPitch;
		if (unitPerPixel > 100)
		{
			unit = "um";
			unitPerPixel /= 1000.0;
		}
		cal.setUnit(unit);
		cal.pixelHeight = cal.pixelWidth = unitPerPixel;
		imp.setCalibration(cal);

		imp.setDimensions(1, 1, newStack.getSize());
		imp.resetDisplayRange();
		imp.updateAndDraw();

		saveImage(imp);

		results.setSource(new IJImageSource(imp));
		results.setName(CREATE_DATA_IMAGE_TITLE + " (" + TITLE + ")");
		results.setConfiguration(createConfiguration((float) psfSD));
		results.setBounds(new Rectangle(0, 0, settings.size, settings.size));
		MemoryPeakResults.addResults(results);

		setBenchmarkResults(imp, results);

		if (benchmarkMode && benchmarkParameters != null)
			benchmarkParameters.setPhotons(results);

		List<LocalisationModel> localisations = toLocalisations(localisationSets);

		savePulses(localisations, results, CREATE_DATA_IMAGE_TITLE);

		// Saved the fixed and moving localisations into different datasets
		saveFixedAndMoving(results, CREATE_DATA_IMAGE_TITLE);

		return localisations;
	}

	private synchronized void addPhotons(double p)
	{
		photonStats.addValue(p);
	}

	/**
	 * Create a PSF model from the image that contains all the z-slices needed to draw the given localisations
	 * 
	 * @param localisationSets
	 * @return
	 */
	private ImagePSFModel createImagePSF(List<LocalisationModelSet> localisationSets)
	{
		ImagePlus imp = WindowManager.getImage(settings.psfImageName);
		if (imp == null)
		{
			IJ.error(TITLE, "Unable to create the PSF model from image: " + settings.psfImageName);
			return null;
		}
		try
		{
			Object o = XmlUtils.fromXML(imp.getProperty("Info").toString());
			if (!(o != null && o instanceof PSFSettings))
				throw new RuntimeException("Unknown PSF settings for image: " + imp.getTitle());
			PSFSettings psfSettings = (PSFSettings) o;

			// Check all the settings have values
			if (psfSettings.nmPerPixel <= 0)
				throw new RuntimeException("Missing nmPerPixel calibration settings for image: " + imp.getTitle());
			if (psfSettings.nmPerSlice <= 0)
				throw new RuntimeException("Missing nmPerSlice calibration settings for image: " + imp.getTitle());
			if (psfSettings.zCentre <= 0)
				throw new RuntimeException("Missing zCentre calibration settings for image: " + imp.getTitle());
			if (psfSettings.fwhm <= 0)
				throw new RuntimeException("Missing FWHM calibration settings for image: " + imp.getTitle());

			// To save memory construct the Image PSF using only the slices that are within 
			// the depth of field of the simulation
			double minZ = Double.POSITIVE_INFINITY, maxZ = Double.NEGATIVE_INFINITY;
			for (LocalisationModelSet l : localisationSets)
			{
				for (LocalisationModel m : l.getLocalisations())
				{
					final double z = m.getZ();
					if (minZ > z)
						minZ = z;
					if (maxZ < z)
						maxZ = z;
				}
			}

			int nSlices = imp.getStackSize();
			// z-centre should be an index and not the ImageJ slice number so subtract 1
			int zCentre = psfSettings.zCentre - 1;

			// Calculate the start/end slices to cover the depth of field
			// This logic must match the ImagePSFModel.
			final double unitsPerSlice = psfSettings.nmPerSlice / settings.pixelPitch;
			// We assume the PSF was imaged axially with increasing z-stage position (moving the stage 
			// closer to the objective). Thus higher z-coordinate are for higher slice numbers.
			int lower = (int) Math.round(minZ / unitsPerSlice) + zCentre;
			int upper = (int) Math.round(maxZ / unitsPerSlice) + zCentre;
			upper = (upper < 0) ? 0 : (upper >= nSlices) ? nSlices - 1 : upper;
			lower = (lower < 0) ? 0 : (lower >= nSlices) ? nSlices - 1 : lower;

			// We cannot just extract the correct slices since the 
			// Image PSF requires the z-centre for normalisation
			if (!(lower <= zCentre && upper >= zCentre))
			{
				// Ensure we include the zCentre
				lower = Math.min(lower, zCentre);
				upper = Math.max(upper, zCentre);
			}

			final double noiseFraction = 1e-3;
			final ImagePSFModel model = new ImagePSFModel(extractImageStack(imp, lower, upper), zCentre - lower,
					psfSettings.nmPerPixel / settings.pixelPitch, unitsPerSlice, psfSettings.fwhm, noiseFraction);

			// Add the calibrated centres
			if (psfSettings.offset != null)
			{
				final int sliceOffset = lower + 1;
				for (PSFOffset offset : psfSettings.offset)
				{
					model.setRelativeCentre(offset.slice - sliceOffset, offset.cx, offset.cy);
				}
			}

			// Initialise the HWHM table so that it can be cloned
			model.initialiseHWHM();

			return model;
		}
		catch (Exception e)
		{
			IJ.error(TITLE, "Unable to create the image PSF model:\n" + e.getMessage());
			return null;
		}
	}

	public static float[][] extractImageStack(ImagePlus imp, int start, int end)
	{
		int size = end - start + 1;
		ImageStack stack = imp.getImageStack();
		float[][] image = new float[size][];
		for (int i = 0; i < image.length; i++)
		{
			image[i] = (float[]) stack.getProcessor(i + start + 1).toFloat(0, null).getPixels();
		}
		return image;
	}

	private PSFModel createPSFModel(ImagePSFModel imagePSFModel)
	{
		if (imagePSF)
		{
			PSFModel copy = imagePSFModel.copy();
			copy.setRandomGenerator(createRandomGenerator());
			return copy;
		}
		else if (settings.psfModel.equals(PSF_MODELS[0]))
		{
			// Calibration based on imaging fluorescent beads at 20nm intervals.
			// Set the PSF to 1.5 x FWHM at 450nm
			double sd = getPsfSD();
			return new GaussianPSFModel(createRandomGenerator(), sd, sd, 450.0 / settings.pixelPitch);
		}
		else
		{
			// Airy pattern
			double width = getPsfSD() / PSFCalculator.AIRY_TO_GAUSSIAN;
			AiryPSFModel m = new AiryPSFModel(createRandomGenerator(), width, width, 450.0 / settings.pixelPitch);
			m.setRing(2);
			return m;
		}
	}

	private synchronized void showProgress()
	{
		IJ.showProgress(frame++, totalFrames);
	}

	private class Spot
	{
		final double[] psf;
		final int x0min, x0max, x1min, x1max;
		final int[] samplePositions;

		public Spot(double[] psf, int x0min, int x0max, int x1min, int x1max, int[] samplePositions)
		{
			this.psf = psf;
			this.x0min = x0min;
			this.x0max = x0max;
			this.x1min = x1min;
			this.x1max = x1max;
			this.samplePositions = samplePositions;
		}
	}

	/**
	 * Use a runnable for the image generation to allow multi-threaded operation. Input parameters
	 * that are manipulated should have synchronized methods.
	 */
	private class ImageGenerator implements Runnable
	{
		final List<LocalisationModelSet> localisations;
		final List<LocalisationModelSet> newLocalisations;
		final int startIndex;
		final int t;
		final PSFModel psfModel;
		final MemoryPeakResults results;
		final ImageStack stack;
		final boolean poissonNoise;
		final RandomDataGenerator random;

		public ImageGenerator(final List<LocalisationModelSet> localisationSets,
				List<LocalisationModelSet> newLocalisations, int startIndex, int t, PSFModel psfModel,
				MemoryPeakResults results, ImageStack stack, boolean poissonNoise, RandomDataGenerator random)
		{
			this.localisations = localisationSets;
			this.newLocalisations = newLocalisations;
			this.startIndex = startIndex;
			this.t = t;
			this.psfModel = psfModel;
			this.results = results;
			this.stack = stack;
			this.poissonNoise = poissonNoise;
			this.random = random;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Runnable#run()
		 */
		public void run()
		{
			if (Utils.isInterrupted())
				return;

			final double psfSD = getPsfSD();

			showProgress();

			final boolean checkSNR = minSNRt1 > 0 || minSNRtN > 0;
			final double totalGain = (settings.getTotalGain() > 0) ? settings.getTotalGain() : 1;

			// Adjust XY dimensions since they are centred on zero
			final double xoffset = settings.size * 0.5;

			float[] image = createBackground(random);
			float[] imageCache = Arrays.copyOf(image, image.length);

			// Create read noise now so that we can calculate the true background noise  
			float[] imageReadNoise = new float[image.length];
			if (settings.readNoise > 0)
			{
				// Read noise is in electrons. Apply camera gain to get the noise in ADUs.
				float readNoise = (float) settings.readNoise;
				if (settings.getCameraGain() != 0)
					readNoise *= settings.getCameraGain();

				RandomGenerator r = random.getRandomGenerator();
				for (int i = 0; i < imageReadNoise.length; i++)
					imageReadNoise[i] += readNoise * r.nextGaussian();
			}

			// Extract the localisations and draw if we have a PSF model
			int fromIndex = findIndexByTime(localisations, startIndex, t);
			if (fromIndex > -1 && psfModel != null)
			{
				int toIndex = findLastIndexByTime(localisations, fromIndex, t);
				List<LocalisationModelSet> subset = localisations.subList(fromIndex, toIndex + 1);
				float[] data = new float[settings.size * settings.size];
				for (LocalisationModelSet localisationSet : subset)
				{
					if (Utils.isInterrupted())
						return;

					if (localisationSet.size() == 0)
						continue;

					// Draw each localisation in the set. Store the PSF so we can remove it later
					double totalPhotonsRendered = 0;
					Spot[] spots = new Spot[localisationSet.size()];
					int spotCount = 0;
					for (LocalisationModel localisation : localisationSet.getLocalisations())
					{
						// Adjust to centre of image
						double[] xyz = localisation.getCoordinates();
						xyz[0] += xoffset;
						xyz[1] += xoffset;

						//addRaw(localisation.getIntensity());

						double photonsRendered = 0;
						double intensity = localisation.getIntensity();
						int[] samplePositions = null;
						if (intensity > 0)
						{
							// TODO record all the positions we draw and the number of photons.
							// This can be used later to compute the Fisher information for each spot.
							if (poissonNoise)
							{
								final int samples = (int) random.nextPoisson(localisation.getIntensity());
								intensity = samples;
								photonsRendered = psfModel.sample3D(data, settings.size, settings.size, samples,
										localisation.getX(), localisation.getY(), localisation.getZ());
								samplePositions = psfModel.getSamplePositions();
							}
							else
							{
								intensity = localisation.getIntensity();
								photonsRendered = psfModel.create3D(data, settings.size, settings.size, intensity,
										localisation.getX(), localisation.getY(), localisation.getZ(), false);
							}
						}
						//addDraw(photons);
						if (photonsRendered > 0)
						{
							totalPhotonsRendered += photonsRendered;
							spots[spotCount++] = new Spot(psfModel.getPSF(), psfModel.getX0min(), psfModel.getX0max(),
									psfModel.getX1min(), psfModel.getX1max(), samplePositions);
						}

						// Update the intensity using the gain.
						// Use the sampled intensity and not the photons rendered. This is the intensity that should be
						// fitted by any function irrespective of whether the photons were actually sampled on the image.
						localisation.setIntensity(intensity * totalGain);
					}

					// Skip if nothing has been drawn. Note that if the localisation set is skipped then the 
					// intensity must be set to zero to prevent the SNR checks using the eliminated neighbours.
					if (totalPhotonsRendered == 0)
					{
						photonsRemoved.incrementAndGet();
						localisationSet.setData(new double[5]);
						continue;
					}
					if (totalPhotonsRendered < minPhotons)
					{
						photonsRemoved.incrementAndGet();
						for (int i = 0; i < spotCount; i++)
						{
							erase(data, spots[i]);
						}
						localisationSet.setData(new double[5]);
						continue;
					}

					addPhotons(totalPhotonsRendered);

					LocalisationModel localisation = localisationSet.toLocalisation();

					// Account for gain 
					totalPhotonsRendered *= totalGain;

					// Add to memory. 0.5 is the centre of the pixel so just round down.
					int origX = (int) localisation.getX();
					int origY = (int) localisation.getY();
					float[] params = new float[7];
					// Background and noise should be calculated using the
					// region covered by the PSF.
					double[] localStats = getStatistics(imageCache, imageReadNoise, origX, origY);
					params[Gaussian2DFunction.BACKGROUND] = (float) (localStats[0] * totalGain + settings.bias);
					params[Gaussian2DFunction.X_POSITION] = (float) localisation.getX();
					params[Gaussian2DFunction.Y_POSITION] = (float) localisation.getY();

					// Note: The width estimate does not account for diffusion
					if (psfModel instanceof GaussianPSFModel)
					{
						GaussianPSFModel m = (GaussianPSFModel) psfModel;
						params[Gaussian2DFunction.X_SD] = (float) m.getS0();
						params[Gaussian2DFunction.Y_SD] = (float) m.getS1();
					}
					else if (psfModel instanceof AiryPSFModel)
					{
						AiryPSFModel m = (AiryPSFModel) psfModel;
						params[Gaussian2DFunction.X_SD] = (float) (m.getW0() * AiryPattern.FACTOR);
						params[Gaussian2DFunction.Y_SD] = (float) (m.getW1() * AiryPattern.FACTOR);
					}
					else if (psfModel instanceof ImagePSFModel)
					{
						ImagePSFModel m = (ImagePSFModel) psfModel;
						params[Gaussian2DFunction.X_SD] = (float) (m.getHWHM0() / Gaussian2DFunction.SD_TO_HWHM_FACTOR);
						params[Gaussian2DFunction.Y_SD] = (float) (m.getHWHM1() / Gaussian2DFunction.SD_TO_HWHM_FACTOR);
					}
					else
					{
						params[Gaussian2DFunction.X_SD] = params[Gaussian2DFunction.Y_SD] = (float) psfSD;
					}
					// Use the actual intensity (not the total photons rendered)
					params[Gaussian2DFunction.SIGNAL] = (float) localisation.getIntensity();

					// The variance of the background image is currently in photons^2. Apply gain to convert to ADUs. 
					double backgroundVariance = localStats[1] * totalGain * totalGain;

					// EM-gain noise factor: Adds sqrt(2) to the electrons input to the register.
					// All data 'read' through the EM-register must have this additional noise factor added.
					if (settings.getEmGain() > 1)
					{
						backgroundVariance *= 2; // Note we are using the variance (std.dev^2) so we use the factor 2
					}

					// Get the actual read noise applied to this part of the image
					double readVariance = localStats[3];

					// *-*-*-*-*
					// Note that the noise we are calculating is the noise that would be in the image with no
					// fluorophore present. This is the true background noise and it is the noise that is  
					// estimated by the Peak Fit plugin. This noise therefore IGNORES THE SHOT NOISE of the 
					// fluorophore SIGNAL.
					// *-*-*-*-*

					// Overall noise can be calculated from the ‘root of sum of squares’ equation
					final double totalNoise = Math.sqrt(backgroundVariance + readVariance);

					// Ensure the new data is added before the intensity is updated. This avoids 
					// syncronisation clashes in the getIntensity(...) function.
					// Use the total photons rendered for signal filtering.
					localisationSet.setData(new double[] { localStats[0], totalNoise, params[Gaussian2DFunction.X_SD],
							params[Gaussian2DFunction.Y_SD], totalPhotonsRendered });

					if (checkSNR)
					{
						if (badLocalisation(localisationSet, totalPhotonsRendered, totalNoise))
						{
							for (int i = 0; i < spotCount; i++)
							{
								erase(data, spots[i]);
							}
							localisationSet.setData(new double[5]);
							continue;
						}
					}

					newLocalisations.add(localisationSet);
					// Use extended result to store the ID.
					// Store the z position in the error.
					results.addSync(new IdPeakResult(t, origX, origY, 0, localisation.getZ(), (float) totalNoise,
							params, null, localisationSet.getId()));
				}

				for (int i = 0; i < image.length; i++)
					image[i] += data[i];
			}

			// Quantum efficiency: Model using binomial distribution
			if (settings.getQuantumEfficiency() < 1)
			{
				final double qe = settings.getQuantumEfficiency();
				for (int i = 0; i < image.length; i++)
					image[i] = random.nextBinomial((int) image[i], qe);
			}

			// Apply EM gain and add Gaussian read noise after all the photons have been simulated
			final boolean tubbsModel = false;
			if (settings.getEmGain() > 1) // This could be >=1 but the rest of the code ignores EM-gain if it is <=1
			{
				// See: https://www.andor.com/learning-academy/sensitivity-making-sense-of-sensitivity
				// there is a statistical variation in the overall number of electrons generated from an initial 
				// charge packet by the gain register. This uncertainty is quantified by a parameter called "Noise Factor" 
				// and detailed theoretical and measured analysis has placed this Noise Factor at a value of √2 (or 1.41) 
				// for EMCCD technology.

				// A Stochastic Model for Electron Multiplication Charge-Coupled Devices – From Theory to Practice
				// (PLoS One. 2013; 8(1): e53671)
				// PGN model:
				// - Poisson for photon shot noise
				// - Gamma for EM gain
				// - Normal for read noise
				// EM gain is essentially a repeated loop of input N and get out M where the each N has a probability of 
				// being amplified. This has been modelled as a series of Poisson or Binomial trials and then the curve 
				// fitted.

				if (tubbsModel)
				{
					// Tubbs's model
					// Equation 14: is a gamma distribution for electrons created in the register 
					for (int i = 0; i < image.length; i++)
					{
						if (image[i] <= 0)
							continue;
						final double scale = settings.getEmGain() - 1 + 1 / image[i];
						final double electrons = random.nextGamma(image[i], scale) - 1;
						image[i] += electrons;
					}
				}
				else
				{
					// Standard gamma distribution
					// This is what is modelled in the Poisson-Gamma-Gaussian likelihood model

					// Since the call random.nextGamma(...) creates a Gamma distribution 
					// which pre-calculates factors only using the scale parameter we 
					// create a custom gamma distribution where the shape can be set as a property.
					CustomGammaDistribution dist = new CustomGammaDistribution(random.getRandomGenerator(), 1,
							settings.getEmGain(), GammaDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);

					for (int i = 0; i < image.length; i++)
					{
						if (image[i] <= 0)
						{
							// What about modelling spontaneous electron events?
							continue;
						}
						//image[i] = (float) random.nextGamma(image[i], settings.getEmGain());
						dist.setShapeUnsafe(image[i]);
						image[i] = (float) dist.sample();
					}
				}
			}

			// Apply camera gain. Note that the noise component of the camera gain is the 
			// read noise. Thus the read noise may change for each camera gain.
			if (settings.getCameraGain() > 0)
			{
				for (int i = 0; i < image.length; i++)
					image[i] *= settings.getCameraGain();
			}

			// Apply read noise (in ADUs)
			if (settings.readNoise > 0)
			{
				for (int i = 0; i < image.length; i++)
					image[i] += imageReadNoise[i];
			}

			{
				// TODO Compute the Fisher information for each spot.
				// Iaa = sum(i) (dYi da) * (dYi da) / Yi
				// 1. Draw each spot perfectly on a new image using psfModel.create3D
				// 2. Apply perfect gain to get the correct scale 
				// 3. For each spot
				// - Combine the other spots + background photons to get the background function
				// - Combine with the target spot and apply gain to get Yi	
				// - Render the target spot using an offset +/- delta on X,Y,Signal (a) to get dYi da
				// - Compute Iaa

				// The simulated image is not used?

			}

			// Add the bias
			for (int i = 0; i < image.length; i++)
				image[i] += settings.bias;

			// Send to output
			stack.setPixels(image, t);
		}

		private void erase(float[] data, Spot spot)
		{
			if (spot.samplePositions != null)
			{
				psfModel.eraseSample(data, settings.size, settings.size, spot.samplePositions);
			}
			else
			{
				psfModel.erase(data, settings.size, settings.size, spot.psf, spot.x0min, spot.x0max, spot.x1min,
						spot.x1max);
			}
		}

		/**
		 * Compute the mean and variance for image 1 and image 2 for the region where the spot was inserted.
		 * The region is defined by an area around the origin.
		 * 
		 * @param image1
		 * @param image2
		 * @param origY
		 * @param origX
		 * @return [mean1, variance1, mean2, variance2]
		 */
		private double[] getStatistics(float[] image1, float[] image2, int origX, int origY)
		{
			// When the random sample of the PSF is small the min and max sample positions
			// may not represent the spot region. Ensure we use an area around the origin.
			final int n = 3;
			final int width = FastMath.min(settings.size, 2 * n + 1);
			int x0min = limit(origX - n);
			int x0max = limit(origX + n);
			int x1min = limit(origY - n);
			int x1max = limit(origY + n);

			// Check we have enough pixels, i.e. in case the origin is outside the image
			int x0range = x0max - x0min + 1;
			while (x0range < width)
			{
				x0min = limit(x0min - 1);
				x0max = limit(x0max + 1);
				x0range = x0max - x0min + 1;
			}
			int x1range = x1max - x1min + 1;
			while (x1range < width)
			{
				x1min = limit(x1min - 1);
				x1max = limit(x1max + 1);
				x1range = x1max - x1min + 1;
			}

			Statistics sum = new Statistics();
			Statistics sum2 = new Statistics();
			for (int y = 0; y < x1range; y++)
			{
				// Locate the insert location
				int indexTo = (y + x1min) * settings.size + x0min;
				for (int x = 0; x < x0range; x++)
				{
					sum.add(image1[indexTo]);
					sum2.add(image2[indexTo]);
					indexTo++;
				}
			}
			return new double[] { sum.getMean(), sum.getVariance(), sum2.getMean(), sum2.getVariance() };
		}
	}

	private int limit(int coord)
	{
		if (coord < 0)
			return 0;
		if (coord >= settings.size)
			return settings.size - 1;
		return coord;
	}

	/**
	 * Check if the localisation, or its neighbours, reach the SNR thresholds. The intensity and noise are after EM-gain
	 * has been applied.
	 * 
	 * @param localisationSet
	 * @param intensity
	 * @param noise
	 * @return
	 */
	public boolean badLocalisation(LocalisationModelSet localisationSet, double intensity, double noise)
	{
		// Set the minimum SNR for either a single spot or for a spot next to a brighter neighbour
		double minSNR = settings.minSNRt1;
		AtomicInteger counter = t1Removed;

		if (localisationSet.hasNeighbour())
		{
			double nextIntensity = getIntensity(localisationSet.getNext());
			double previousIntensity = getIntensity(localisationSet.getPrevious());

			// Check if either neighbour is above the t1 threshold
			if ((nextIntensity / noise > settings.minSNRt1) || (previousIntensity / noise > settings.minSNRt1))
			{
				// If neighbours are bright then use a more lenient threshold
				minSNR = settings.minSNRtN;
				counter = tNRemoved;
			}
		}

		if (intensity / noise < minSNR)
		{
			counter.incrementAndGet();
			return true;
		}
		return false;
	}

	private double getIntensity(LocalisationModelSet localisationSet)
	{
		if (localisationSet != null)
			return localisationSet.getData()[4];
		return 0;
	}

	/**
	 * Create a dummy calibration with the initial PSF width. This is used by the SpotInspector
	 * 
	 * @param psfWidth
	 * @return
	 */
	private String createConfiguration(float psfWidth)
	{
		FitConfiguration fitConfig = new FitConfiguration();
		fitConfig.setInitialPeakStdDev0(psfWidth);
		fitConfig.setInitialPeakStdDev1(psfWidth);
		FitEngineConfiguration config = new FitEngineConfiguration(fitConfig);
		return XmlUtils.toXML(config);
	}

	private float[] backgroundPixels = null;
	private boolean uniformBackground = false;

	private float[] createBackground(RandomDataGenerator random)
	{
		float[] pixels2 = null;

		if (settings.background > 0)
		{
			if (random == null)
				random = new RandomDataGenerator(createRandomGenerator());
			createBackgroundPixels();
			pixels2 = Arrays.copyOf(backgroundPixels, backgroundPixels.length);

			// Add Poisson noise.
			if (uniformBackground)
			{
				// We can do N random samples thus ensuring the background average is constant.
				// Note: The number of samples must be Poisson distributed.
				// currently: pixels2[0] = uniform background level
				// => (pixels2[0] * pixels2.length) = amount of photons that fall on the image. 
				final int samples = (int) random.nextPoisson(pixels2[0] * pixels2.length);

				// Only do sampling if the number of samples is valid
				if (samples >= 1)
				{
					pixels2 = new float[pixels2.length];
					final int upper = pixels2.length - 1;
					for (int i = 0; i < samples; i++)
					{
						pixels2[random.nextInt(0, upper)] += 1;
					}
				}
				else
				{
					// If using a uniform illumination then we can use a fixed Poisson distribution
					PoissonDistribution dist = new PoissonDistribution(random.getRandomGenerator(), pixels2[0],
							PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);
					for (int i = 0; i < pixels2.length; i++)
					{
						pixels2[i] = dist.sample();
					}
				}
			}
			else
			{
				for (int i = 0; i < pixels2.length; i++)
				{
					pixels2[i] = random.nextPoisson(pixels2[i]);
				}
			}
		}
		else
		{
			pixels2 = new float[settings.size * settings.size];
		}

		return pixels2;
	}

	synchronized private void createBackgroundPixels()
	{
		// Cache illumination background
		if (backgroundPixels == null)
		{
			backgroundPixels = new float[settings.size * settings.size];

			ImagePlus imp = WindowManager.getImage(settings.backgroundImage);
			if (imp != null)
			{
				// Use an image for the background
				ImageProcessor ip = imp.getProcessor().duplicate().toFloat(0, null);
				ip.setInterpolationMethod(ImageProcessor.BILINEAR);
				ip = ip.resize(settings.size, settings.size);
				float[] data = (float[]) ip.getPixels();
				final double max = FastMath.max(0, Maths.max(data));
				if (max != 0)
				{
					final double scale = settings.background / max;
					for (int i = 0; i < backgroundPixels.length; i++)
					{
						// Ignore pixels below zero
						backgroundPixels[i] = (float) (FastMath.max(0, data[i]) * scale);
					}
					return;
				}
			}

			// Use the illumination (this is the fall-back method if the background image has no
			// maximum)
			SpatialIllumination illumination = createIllumination(settings.background, 0);
			double[] xyz = new double[3];
			for (int y = 0, i = 0; y < settings.size; y++)
			{
				xyz[1] = y - settings.size / 2;
				for (int x = 0, x2 = -settings.size / 2; x < settings.size; x++, x2++, i++)
				{
					xyz[0] = x2;
					backgroundPixels[i] = (float) illumination.getPhotons(xyz);
				}
			}
		}
	}

	/**
	 * Find the first index from the starting index where the localisation matches the time
	 * 
	 * @param localisations
	 * @param fromIndex
	 *            start index
	 * @param t
	 *            time
	 * @return the index (or -1)
	 */
	private int findIndexByTime(List<LocalisationModelSet> localisations, int fromIndex, int t)
	{
		while (fromIndex < localisations.size() && localisations.get(fromIndex).getTime() != t)
			fromIndex++;
		return fromIndex >= localisations.size() ? -1 : fromIndex;
	}

	/**
	 * Find the last index from the starting index where the localisation matches the time
	 * 
	 * @param localisations
	 * @param fromIndex
	 *            start index
	 * @param t
	 *            time
	 * @return the index (or -1)
	 */
	private int findLastIndexByTime(List<LocalisationModelSet> localisations, int fromIndex, int t)
	{
		// Check the start point is valid
		if (localisations.get(fromIndex).getTime() != t)
		{
			fromIndex = findIndexByTime(localisations, 0, t);
			if (fromIndex == -1)
				return fromIndex;
		}
		while (fromIndex < localisations.size() && localisations.get(fromIndex).getTime() == t)
			fromIndex++;
		return fromIndex - 1;
	}

	/**
	 * Find the first index from the starting index where the localisation matches the id
	 * 
	 * @param localisations
	 * @param fromIndex
	 *            start index
	 * @param id
	 * @return the index (or -1)
	 */
	private int findIndexById(List<LocalisationModel> localisations, int fromIndex, int id)
	{
		while (fromIndex < localisations.size() && localisations.get(fromIndex).getId() != id)
			fromIndex++;
		return fromIndex >= localisations.size() ? -1 : fromIndex;
	}

	/**
	 * Find the last index from the starting index where the localisation matches the Id
	 * 
	 * @param localisations
	 * @param fromIndex
	 *            start index
	 * @param id
	 * @return the index (or -1)
	 */
	private int findLastIndexById(List<LocalisationModel> localisations, int fromIndex, int id)
	{
		// Check the start point is valid
		if (localisations.get(fromIndex).getId() != id)
		{
			fromIndex = findIndexById(localisations, 0, id);
			if (fromIndex == -1)
				return fromIndex;
		}
		while (fromIndex < localisations.size() && localisations.get(fromIndex).getId() == id)
			fromIndex++;
		return fromIndex - 1;
	}

	private List<LocalisationModel> toLocalisations(List<LocalisationModelSet> localisationSets)
	{
		ArrayList<LocalisationModel> localisations = new ArrayList<LocalisationModel>(localisationSets.size());
		for (LocalisationModelSet s : localisationSets)
			localisations.add(s.toLocalisation());
		return localisations;
	}

	/**
	 * Remove all fluorophores which were not drawn
	 * 
	 * @param fluorophores
	 * @param localisations
	 * @return
	 */
	private List<? extends FluorophoreSequenceModel> removeFilteredFluorophores(
			List<? extends FluorophoreSequenceModel> fluorophores, List<LocalisationModel> localisations)
	{
		if (fluorophores == null)
			return null;
		// movingMolecules will be created with an initial capacity to hold all the unique IDs
		TIntHashSet idSet = new TIntHashSet((movingMolecules != null) ? movingMolecules.capacity() : 0);
		for (LocalisationModel l : localisations)
			idSet.add(l.getId());
		List<FluorophoreSequenceModel> newFluorophores = new ArrayList<FluorophoreSequenceModel>(idSet.size());
		for (FluorophoreSequenceModel f : fluorophores)
		{
			if (idSet.contains(f.getId()))
				newFluorophores.add(f);
		}
		return newFluorophores;
	}

	private double showSummary(List<? extends FluorophoreSequenceModel> fluorophores,
			List<LocalisationModel> localisations)
	{
		IJ.showStatus("Calculating statistics ...");

		createSummaryTable();

		Statistics[] stats = new Statistics[NAMES.length];
		for (int i = 0; i < stats.length; i++)
		{
			stats[i] = (settings.showHistograms || alwaysRemoveOutliers[i]) ? new StoredDataStatistics()
					: new Statistics();
		}

		// Find the largest timepoint
		ImagePlus outputImp = WindowManager.getImage(benchmarkImageId);
		int nFrames;
		if (outputImp == null)
		{
			sortLocalisationsByTime(localisations);
			nFrames = localisations.get(localisations.size() - 1).getTime();
		}
		else
		{
			nFrames = outputImp.getStackSize();
		}
		int[] countHistogram = new int[nFrames + 1];

		// Use the localisations that were drawn to create the sampled on/off times
		rebuildNeighbours(localisations);

		// Assume that there is at least one localisation
		LocalisationModel first = localisations.get(0);
		int currentId = first.getId(); // The current localisation
		int lastT = first.getTime(); // The last time this localisation was on
		int blinks = 0; // Number of blinks
		int currentT = 0; // On-time of current pulse
		double signal = 0;
		final double centreOffset = settings.size * 0.5;
		// Used to convert the sampled times in frames into seconds
		final double framesPerSecond = 1000.0 / settings.exposureTime;
		final double gain = (settings.getTotalGain() > 0) ? settings.getTotalGain() : 1;
		for (LocalisationModel l : localisations)
		{
			if (l.getData() == null)
				System.out.println("No localisation data. This should not happen!");
			final double noise = (l.getData() != null) ? l.getData()[1] : 1;
			final double intensity = (l.getData() != null) ? l.getData()[4] : l.getIntensity();
			final double intensityInPhotons = intensity / gain;
			// Q. What if the noise is zero, i.e. no background photon / read noise?
			// Just ignore it at current.
			final double snr = intensity / noise;
			stats[SIGNAL].add(intensityInPhotons);
			stats[NOISE].add(noise / gain);
			if (noise != 0)
				stats[SNR].add(snr);
			// Average intensity only from continuous spots.
			// The continuous flag is for spots that have all the simulation steps continuously on.
			// Try using the neighbour pointers instead to get the 'sampled' continuous spots.
			//if (l.isContinuous())
			if (l.getNext() != null && l.getPrevious() != null)
			{
				stats[SIGNAL_CONTINUOUS].add(intensityInPhotons);
				if (noise != 0)
					stats[SNR_CONTINUOUS].add(snr);
			}

			int id = l.getId();
			// Check if this a new fluorophore
			if (currentId != id)
			{
				// Add previous fluorophore
				stats[SAMPLED_BLINKS].add(blinks);
				stats[SAMPLED_T_ON].add(currentT / framesPerSecond);
				stats[TOTAL_SIGNAL].add(signal);

				// Reset
				blinks = 0;
				currentT = 1;
				currentId = id;
				signal = intensityInPhotons;
			}
			else
			{
				signal += intensityInPhotons;
				// Check if the current fluorophore pulse is broken (i.e. a blink)
				if (l.getTime() - 1 > lastT)
				{
					blinks++;
					stats[SAMPLED_T_ON].add(currentT / framesPerSecond);
					currentT = 1;
					stats[SAMPLED_T_OFF].add(((l.getTime() - 1) - lastT) / framesPerSecond);
				}
				else
				{
					// Continuous on-time
					currentT++;
				}
			}

			lastT = l.getTime();
			countHistogram[lastT]++;

			stats[X].add((l.getX() - centreOffset) * settings.pixelPitch);
			stats[Y].add((l.getY() - centreOffset) * settings.pixelPitch);
			stats[Z].add(l.getZ() * settings.pixelPitch);
		}
		// Final fluorophore
		stats[SAMPLED_BLINKS].add(blinks);
		stats[SAMPLED_T_ON].add(currentT / framesPerSecond);
		stats[TOTAL_SIGNAL].add(signal);

		// Samples per frame
		for (int t = 1; t < countHistogram.length; t++)
			stats[SAMPLES].add(countHistogram[t]);

		if (fluorophores != null)
		{
			for (FluorophoreSequenceModel f : fluorophores)
			{
				stats[BLINKS].add(f.getNumberOfBlinks());
				// On-time
				for (double t : f.getOnTimes())
					stats[T_ON].add(t);
				// Off-time
				for (double t : f.getOffTimes())
					stats[T_OFF].add(t);
			}
		}
		else
		{
			// show no blinks
			stats[BLINKS].add(0);
			stats[T_ON].add(1);
			//stats[T_OFF].add(0);
		}

		if (results != null)
		{
			final boolean emCCD = (settings.getEmGain() > 1);
			// Convert depth-of-field to pixels
			final double depth = settings.depthOfField / settings.pixelPitch;
			for (PeakResult r : results.getResults())
			{
				final double precision = r.getPrecision(settings.pixelPitch, gain, emCCD);
				stats[PRECISION].add(precision);
				// The error stores the z-depth in pixels
				if (Math.abs(r.error) < depth)
					stats[PRECISION_IN_FOCUS].add(precision);
				stats[WIDTH].add(r.getSD());
			}
			// Compute density per frame. Multithread for speed
			if (settings.densityRadius > 0)
			{
				IJ.showStatus("Calculating density ...");
				ExecutorService threadPool = Executors.newFixedThreadPool(Prefs.getThreads());
				List<Future<?>> futures = new LinkedList<Future<?>>();
				final ArrayList<float[]> coords = new ArrayList<float[]>();
				int t = results.getHead().getFrame();
				final Statistics densityStats = stats[DENSITY];
				final float radius = (float) (settings.densityRadius * getHWHM());
				final Rectangle bounds = results.getBounds();
				currentIndex = 0;
				finalIndex = results.getTail().getFrame();
				// Store the density for each result.
				int[] allDensity = new int[results.size()];
				int allIndex = 0;
				for (PeakResult r : results.getResults())
				{
					if (t != r.getFrame())
					{
						allIndex += runDensityCalculation(threadPool, futures, coords, densityStats, radius, bounds,
								allDensity, allIndex);
					}
					coords.add(new float[] { r.getXPosition(), r.getYPosition() });
					t = r.getFrame();
				}
				runDensityCalculation(threadPool, futures, coords, densityStats, radius, bounds, allDensity, allIndex);
				Utils.waitForCompletion(futures);
				threadPool.shutdownNow();
				threadPool = null;
				IJ.showProgress(1);

				// Split results into singles (density = 0) and clustered (density > 0)
				MemoryPeakResults singles = copyMemoryPeakResults("No Density");
				MemoryPeakResults clustered = copyMemoryPeakResults("Density");

				int i = 0;
				for (PeakResult r : results.getResults())
				{
					// Store density in the original value field
					r.origValue = allDensity[i];
					if (allDensity[i++] == 0)
						singles.add(r);
					else
						clustered.add(r);
				}
			}
		}

		StringBuilder sb = new StringBuilder();
		sb.append(datasetNumber).append("\t");
		sb.append((fluorophores == null) ? localisations.size() : fluorophores.size()).append("\t");
		sb.append(stats[SAMPLED_BLINKS].getN() + (int) stats[SAMPLED_BLINKS].getSum()).append("\t");
		sb.append(localisations.size()).append("\t");
		sb.append(nFrames).append("\t");
		sb.append(Utils.rounded(areaInUm)).append("\t");
		sb.append(Utils.rounded(localisations.size() / (areaInUm * nFrames), 4)).append("\t");
		sb.append(Utils.rounded(getHWHM(), 4)).append("\t");
		double s = getPsfSD();
		sb.append(Utils.rounded(s, 4)).append("\t");
		s *= settings.pixelPitch;
		final double sa = PSFCalculator.squarePixelAdjustment(s, settings.pixelPitch) / settings.pixelPitch;
		sb.append(Utils.rounded(sa, 4)).append("\t");
		// Width not valid for the Image PSF
		int nStats = (imagePSF) ? stats.length - 1 : stats.length;
		for (int i = 0; i < nStats; i++)
		{
			double centre = (alwaysRemoveOutliers[i])
					? ((StoredDataStatistics) stats[i]).getStatistics().getPercentile(50) : stats[i].getMean();
			sb.append(Utils.rounded(centre, 4)).append("\t");
		}
		if (java.awt.GraphicsEnvironment.isHeadless())
		{
			IJ.log(sb.toString());
			return stats[SIGNAL].getMean();
		}
		else
		{
			summaryTable.append(sb.toString());
		}

		// Show histograms
		if (settings.showHistograms)
		{
			IJ.showStatus("Calculating histograms ...");
			boolean[] chosenHistograms = getChoosenHistograms();

			WindowOrganiser wo = new WindowOrganiser();

			boolean requireRetile = false;
			for (int i = 0; i < NAMES.length; i++)
			{
				if (chosenHistograms[i])
				{
					wo.add(Utils.showHistogram(TITLE, (StoredDataStatistics) stats[i], NAMES[i],
							(integerDisplay[i]) ? 1 : 0, (settings.removeOutliers || alwaysRemoveOutliers[i]) ? 2 : 0,
							settings.histogramBins * ((integerDisplay[i]) ? 100 : 1)));
					requireRetile = requireRetile || Utils.isNewWindow();
				}
			}

			wo.tile();
		}
		IJ.showStatus("");
		return stats[SIGNAL].getMean();
	}

	private int runDensityCalculation(ExecutorService threadPool, List<Future<?>> futures,
			final ArrayList<float[]> coords, final Statistics densityStats, final float radius, final Rectangle bounds,
			final int[] allDensity, final int allIndex)
	{
		final int size = coords.size();
		final float[] xCoords = new float[size];
		final float[] yCoords = new float[size];
		for (int i = 0; i < xCoords.length; i++)
		{
			float[] xy = coords.get(i);
			xCoords[i] = xy[0];
			yCoords[i] = xy[1];
		}
		futures.add(threadPool.submit(new Runnable()
		{
			public void run()
			{
				incrementProgress();
				final DensityManager dm = new DensityManager(xCoords, yCoords, bounds);
				final int[] density = dm.calculateDensity(radius, true);
				addDensity(densityStats, density);

				// Store the density for each result. This does not need to be synchronised 
				// since the indices in different threads are unique.
				for (int i = 0, index = allIndex; i < density.length; i++, index++)
					allDensity[index] = density[i];
			}
		}));
		coords.clear();
		return size;
	}

	private int currentIndex, finalIndex;

	private synchronized void incrementProgress()
	{
		IJ.showProgress(currentIndex, finalIndex);
	}

	private synchronized void addDensity(Statistics stats, int[] density)
	{
		stats.add(density);
	}

	/**
	 * Copy all the settings from the results into a new results set labelled with the name suffix
	 * 
	 * @param nameSuffix
	 * @return The new results set
	 */
	private MemoryPeakResults copyMemoryPeakResults(String nameSuffix)
	{
		MemoryPeakResults newResults = new MemoryPeakResults();
		newResults.copySettings(this.results);
		newResults.setName(newResults.getSource().getName() + " (" + TITLE + " " + nameSuffix + ")");
		newResults.setSortAfterEnd(true);
		newResults.begin();
		MemoryPeakResults.addResults(newResults);
		return newResults;
	}

	private boolean[] getChoosenHistograms()
	{
		if (settings.chooseHistograms)
			return displayHistograms;

		boolean[] all = new boolean[displayHistograms.length];
		for (int i = 0; i < all.length; i++)
			all[i] = true;
		return all;
	}

	private void sortLocalisationsByIdThenTime(List<LocalisationModel> localisations)
	{
		Collections.sort(localisations, new Comparator<LocalisationModel>()
		{
			public int compare(LocalisationModel o1, LocalisationModel o2)
			{
				// Order by ID then time
				if (o1.getId() < o2.getId())
					return -1;
				if (o1.getId() > o2.getId())
					return 1;
				if (o1.getTime() < o2.getTime())
					return -1;
				if (o1.getTime() > o2.getTime())
					return 1;
				return 0;
			}
		});
	}

	private void sortLocalisationsByTime(List<LocalisationModel> localisations)
	{
		Collections.sort(localisations, new Comparator<LocalisationModel>()
		{
			public int compare(LocalisationModel o1, LocalisationModel o2)
			{
				// Order by n time
				if (o1.getTime() < o2.getTime())
					return -1;
				if (o1.getTime() > o2.getTime())
					return 1;
				return 0;
			}
		});
	}

	/**
	 * Sort by id then time, then rebuild the neighbour pointers.
	 * 
	 * @param localisations
	 */
	private void rebuildNeighbours(List<LocalisationModel> localisations)
	{
		sortLocalisationsByIdThenTime(localisations);

		int id = 0, t = 0;
		LocalisationModel previous = null;
		for (LocalisationModel l : localisations)
		{
			if (l.getId() != id)
			{
				// New spot so no previous neighbour
				previous = null;
			}
			else if (l.getTime() > t + 1)
			{
				// Discontinuous time so no previous neighbour
				previous = null;
			}

			l.setPrevious(previous);
			l.setNext(null);

			id = l.getId();
			t = l.getTime();
			previous = l;
		}
	}

	private void createSummaryTable()
	{
		if (java.awt.GraphicsEnvironment.isHeadless())
		{
			if (header == null)
			{
				header = createHeader();
				IJ.log(header);
			}
		}
		else
		{
			if (summaryTable == null || !summaryTable.isVisible())
			{
				summaryTable = new TextWindow("Data Summary", createHeader(), "", 800, 300);
				summaryTable.setVisible(true);
			}
		}
	}

	private String createHeader()
	{
		StringBuilder sb = new StringBuilder(
				"Dataset\tMolecules\tPulses\tLocalisations\tnFrames\tArea (um^2)\tDensity (mol/um^2)\tHWHM\tS\tSa");
		for (int i = 0; i < NAMES.length; i++)
		{
			sb.append("\t").append(NAMES[i]);
			//if (alwaysRemoveOutliers[i])
			//	sb.append("*");
		}
		return sb.toString();
	}

	/**
	 * Save the image to a TIFF file
	 * 
	 * @param imp
	 */
	private void saveImage(ImagePlus imp)
	{
		if (!settings.saveImage)
			return;
		String[] path = Utils.decodePath(settings.imageFilename);
		OpenDialog chooser = new OpenDialog("Image_File", path[0], path[1]);
		if (chooser.getFileName() != null)
		{
			settings.imageFilename = chooser.getDirectory() + chooser.getFileName();
			settings.imageFilename = Utils.replaceExtension(settings.imageFilename, "tiff");

			FileSaver fs = new FileSaver(imp);
			boolean ok;
			if (imp.getStackSize() > 1)
				ok = fs.saveAsTiffStack(settings.imageFilename);
			else
				ok = fs.saveAsTiff(settings.imageFilename);
			// The following call throws a NoSuchMethodError.
			// ok = IJ.saveAsTiff(imp, settings.imageFilename);

			if (!ok)
				IJ.log("Failed to save image to file: " + settings.imageFilename);
		}
	}

	private void saveImageResults(MemoryPeakResults results)
	{
		if (!settings.saveImageResults)
			return;
		String[] path = Utils.decodePath(settings.imageResultsFilename);
		OpenDialog chooser = new OpenDialog("Image_Results_File", path[0], path[1]);
		if (chooser.getFileName() != null)
		{
			settings.imageResultsFilename = chooser.getDirectory() + chooser.getFileName();
			settings.imageResultsFilename = Utils.replaceExtension(settings.imageResultsFilename, "xls");

			FilePeakResults r = new FilePeakResults(settings.imageResultsFilename, false);
			r.copySettings(results);
			r.setPeakIdColumnName("Frame");
			r.begin();
			r.addAll(results.getResults());
			r.end();
		}
	}

	/**
	 * Create a set of results that represent the molecule continuous on-times (pulses)
	 * 
	 * @param localisations
	 * @param results
	 * @param title
	 */
	private void savePulses(List<LocalisationModel> localisations, MemoryPeakResults results, String title)
	{
		sortLocalisationsByIdThenTime(localisations);

		MemoryPeakResults traceResults = copyMemoryPeakResults("Pulses");
		LocalisationModel start = null;
		int currentId = -1;
		int n = 0;
		float[] params = new float[7];
		double noise = 0;
		int lastT = -1;
		for (LocalisationModel localisation : localisations)
		{
			if (currentId != localisation.getId() || lastT + 1 != localisation.getTime())
			{
				if (n > 0)
				{
					params[Gaussian2DFunction.BACKGROUND] /= n;
					params[Gaussian2DFunction.X_POSITION] /= n;
					params[Gaussian2DFunction.Y_POSITION] /= n;
					params[Gaussian2DFunction.X_SD] /= n;
					params[Gaussian2DFunction.Y_SD] /= n;

					ExtendedPeakResult p = new ExtendedPeakResult(start.getTime(), (int) Math.round(start.getX()),
							(int) Math.round(start.getY()), 0, 0, (float) (Math.sqrt(noise)), params, null, lastT,
							currentId);
					// if (p.getPrecision(107, 1) > 2000)
					// {
					// System.out.printf("Weird precision = %g (%d)\n", p.getPrecision(107, 1), n);
					// }
					traceResults.add(p);
				}
				start = localisation;
				currentId = localisation.getId();
				n = 0;
				params = new float[7];
				noise = 0;
			}

			final double[] data = localisation.getData();
			params[Gaussian2DFunction.BACKGROUND] += data[0];
			params[Gaussian2DFunction.X_POSITION] += localisation.getX();
			params[Gaussian2DFunction.Y_POSITION] += localisation.getY();
			params[Gaussian2DFunction.SIGNAL] += localisation.getIntensity();
			noise += data[1] * data[1];
			params[Gaussian2DFunction.X_SD] += data[2];
			params[Gaussian2DFunction.Y_SD] += data[3];
			n++;
			lastT = localisation.getTime();
		}

		// Final pulse
		if (n > 0)
		{
			params[Gaussian2DFunction.BACKGROUND] /= n;
			params[Gaussian2DFunction.X_POSITION] /= n;
			params[Gaussian2DFunction.Y_POSITION] /= n;
			params[Gaussian2DFunction.X_SD] /= n;
			params[Gaussian2DFunction.Y_SD] /= n;

			traceResults.add(new ExtendedPeakResult(start.getTime(), (int) Math.round(start.getX()),
					(int) Math.round(start.getY()), 0, 0, (float) (Math.sqrt(noise)), params, null, lastT, currentId));
		}

		traceResults.end();
	}

	private void saveFixedAndMoving(MemoryPeakResults results, String title)
	{
		if (simpleMode || benchmarkMode || spotMode)
			return;
		if (settings.diffusionRate <= 0 || settings.fixedFraction >= 1)
			return;

		MemoryPeakResults fixedResults = copyMemoryPeakResults("Fixed");
		MemoryPeakResults movingResults = copyMemoryPeakResults("Moving");

		List<PeakResult> peakResults = results.getResults();
		// Sort using the ID
		Collections.sort(peakResults, new Comparator<PeakResult>()
		{
			public int compare(PeakResult o1, PeakResult o2)
			{
				return o1.getId() - o2.getId();
			}
		});

		int currentId = -1;
		MemoryPeakResults currentResults = movingResults;
		for (PeakResult p : peakResults)
		{
			// The ID was stored in the result's parameter standard deviation array
			if (currentId != p.getId())
			{
				currentId = p.getId();
				currentResults = (movingMolecules.contains(currentId)) ? movingResults : fixedResults;
			}
			currentResults.add(p);
		}

		movingResults.end();
		fixedResults.end();

		// Reset the input results
		results.sort();
	}

	/**
	 * Update the fluorophores relative coordinates to absolute
	 * 
	 * @param molecules
	 */
	@SuppressWarnings("unused")
	private void convertRelativeToAbsolute(List<CompoundMoleculeModel> molecules)
	{
		for (CompoundMoleculeModel c : molecules)
		{
			final double[] xyz = c.getCoordinates();
			for (int n = c.getSize(); n-- > 0;)
			{
				MoleculeModel m = c.getMolecule(n);
				double[] xyz2 = m.getCoordinates();
				for (int i = 0; i < 3; i++)
					xyz2[i] += xyz[i];
			}
		}
	}

	/**
	 * Save the fluorophores to a text file
	 * 
	 * @param fluorophores
	 */
	private void saveFluorophores(List<? extends FluorophoreSequenceModel> fluorophores)
	{
		if (!settings.saveFluorophores || fluorophores == null)
			return;

		String[] path = Utils.decodePath(settings.fluorophoresFilename);
		OpenDialog chooser = new OpenDialog("Fluorophores_File", path[0], path[1]);
		if (chooser.getFileName() != null)
		{
			settings.fluorophoresFilename = chooser.getDirectory() + chooser.getFileName();
			settings.fluorophoresFilename = Utils.replaceExtension(settings.fluorophoresFilename, "xls");

			BufferedWriter output = null;
			try
			{
				output = new BufferedWriter(new FileWriter(settings.fluorophoresFilename));
				output.write(createResultsFileHeader());
				output.write("#Id\tn-Blinks\tStart\tStop\t...");
				output.newLine();
				for (int id = 1; id <= fluorophores.size(); id++)
				{
					FluorophoreSequenceModel f = fluorophores.get(id - 1);
					StringBuffer sb = new StringBuffer();
					sb.append(f.getId()).append("\t");
					sb.append(f.getNumberOfBlinks()).append("\t");
					for (double[] burst : f.getBurstSequence())
					{
						sb.append(Utils.rounded(burst[0], 3)).append("\t").append(Utils.rounded(burst[1], 3))
								.append("\t");
					}
					output.write(sb.toString());
					output.newLine();
				}
			}
			catch (Exception e)
			{
				// Q. Add better handling of errors?
				e.printStackTrace();
				IJ.log("Failed to save fluorophores to file: " + settings.fluorophoresFilename);
			}
			finally
			{
				if (output != null)
				{
					try
					{
						output.close();
					}
					catch (IOException e)
					{
						e.printStackTrace();
					}
				}
			}
		}
	}

	/**
	 * Save the localisations to a text file
	 * 
	 * @param localisations
	 */
	private void saveLocalisations(List<LocalisationModel> localisations)
	{
		if (!settings.saveLocalisations)
			return;

		sortLocalisationsByTime(localisations);

		//		Collections.sort(localisations, new Comparator<LocalisationModel>(){
		//
		//			public int compare(LocalisationModel o1, LocalisationModel o2)
		//			{
		//				int cellx1 = (int)(o1.getX() / settings.cellSize);
		//				int cellx2 = (int)(o2.getX() / settings.cellSize);
		//				int result = cellx2 - cellx1;
		//				if (result != 0)
		//					return result;
		//				int celly1 = (int)(o1.getY() / settings.cellSize);
		//				int celly2 = (int)(o2.getY() / settings.cellSize);
		//				result = celly2 - celly1;
		//				if (result != 0)
		//					return result;
		//				return (o1.getZ() == o2.getZ()) ? 0 : (o1.getZ() == 0) ? -1 : 1;
		//			}});

		String[] path = Utils.decodePath(settings.localisationsFilename);
		OpenDialog chooser = new OpenDialog("Localisations_File", path[0], path[1]);
		if (chooser.getFileName() != null)
		{
			settings.localisationsFilename = chooser.getDirectory() + chooser.getFileName();
			settings.localisationsFilename = Utils.replaceExtension(settings.localisationsFilename, "xls");

			BufferedWriter output = null;
			try
			{
				output = new BufferedWriter(new FileWriter(settings.localisationsFilename));
				output.write(createResultsFileHeader());
				output.write("#T\tId\tX\tY\tZ\tIntensity");
				output.newLine();
				for (LocalisationModel l : localisations)
				{
					StringBuffer sb = new StringBuffer();
					sb.append(l.getTime()).append("\t");
					sb.append(l.getId()).append("\t");
					sb.append(IJ.d2s(l.getX(), 6)).append("\t");
					sb.append(IJ.d2s(l.getY(), 6)).append("\t");
					sb.append(IJ.d2s(l.getZ(), 6)).append("\t");
					sb.append(l.getIntensity());
					output.write(sb.toString());
					output.newLine();
				}
			}
			catch (Exception e)
			{
				// Q. Add better handling of errors?
				e.printStackTrace();
				IJ.log("Failed to save localisations to file: " + settings.localisationsFilename);
			}
			finally
			{
				if (output != null)
				{
					try
					{
						output.close();
					}
					catch (IOException e)
					{
						e.printStackTrace();
					}
				}
			}
		}
	}

	private String createResultsFileHeader()
	{
		if (resultsFileHeader == null)
		{
			String[] backgroundImages = createBackgroundImageList();

			StringBuffer sb = new StringBuffer();
			sb.append("# ").append(TITLE).append(" Parameters:\n");
			addHeaderLine(sb, "Pixel_pitch (nm)", settings.pixelPitch);
			addHeaderLine(sb, "Size", settings.size);
			if (!benchmarkMode)
			{
				addHeaderLine(sb, "Depth", settings.depth);
				addHeaderLine(sb, "Fixed depth", settings.fixedDepth);
			}
			if (!(simpleMode || benchmarkMode))
			{
				if (!trackMode)
					addHeaderLine(sb, "Seconds", settings.seconds);
				addHeaderLine(sb, "Exposure_time", settings.exposureTime);
				addHeaderLine(sb, "Steps_per_second", settings.stepsPerSecond);
				if (!trackMode)
				{
					addHeaderLine(sb, "Illumination", settings.illumination);
					addHeaderLine(sb, "Pulse_interval", settings.pulseInterval);
					addHeaderLine(sb, "Pulse_ratio", settings.pulseRatio);
				}
				if (backgroundImages != null)
					addHeaderLine(sb, "Background_image", settings.backgroundImage);
			}
			addHeaderLine(sb, "Background", settings.background);
			addHeaderLine(sb, "EM_gain", settings.getEmGain());
			addHeaderLine(sb, "Camera_gain", settings.getCameraGain());
			addHeaderLine(sb, "Quantum_efficiency", settings.getQuantumEfficiency());
			addHeaderLine(sb, "Read_noise", settings.readNoise);
			addHeaderLine(sb, "Bias", settings.bias);
			addHeaderLine(sb, "PSF_model", settings.psfModel);
			if (imagePSF)
			{
				addHeaderLine(sb, "PSF_image", settings.psfImageName);
			}
			else
			{
				addHeaderLine(sb, "Wavelength (nm)", settings.wavelength);
				addHeaderLine(sb, "Numerical_aperture", settings.numericalAperture);
			}
			if (!(benchmarkMode || spotMode))
			{
				addHeaderLine(sb, "Distribution", settings.distribution);
				if (settings.distribution.equals(DISTRIBUTION[MASK]))
				{
					addHeaderLine(sb, "Distribution_mask", settings.distributionMask);
				}
				else if (settings.distribution.equals(DISTRIBUTION[GRID]))
				{
					addHeaderLine(sb, "Cell_size", settings.cellSize);
					addHeaderLine(sb, "p-binary", settings.probabilityBinary);
					addHeaderLine(sb, "Min_binary_distance (nm)", settings.minBinaryDistance);
					addHeaderLine(sb, "Max_binary_distance (nm)", settings.maxBinaryDistance);
				}
			}
			addHeaderLine(sb, "Particles", settings.particles);
			if (benchmarkMode)
			{
				addHeaderLine(sb, "X_position", settings.xPosition);
				addHeaderLine(sb, "Y_position", settings.yPosition);
				addHeaderLine(sb, "Z_position", settings.zPosition);
				addHeaderLine(sb, "Min_photons", settings.photonsPerSecond);
				addHeaderLine(sb, "Max_photons", settings.photonsPerSecondMaximum);
			}
			else if (simpleMode)
			{
				addHeaderLine(sb, "Density (um^-2)", settings.density);
				addHeaderLine(sb, "Min_photons", settings.photonsPerSecond);
				addHeaderLine(sb, "Max_photons", settings.photonsPerSecondMaximum);
			}
			else if (spotMode)
			{
				addHeaderLine(sb, "Min_photons", settings.photonsPerSecond);
				addHeaderLine(sb, "Max_photons", settings.photonsPerSecondMaximum);
			}
			else
			{
				addHeaderLine(sb, "Diffusion_rate", settings.diffusionRate);
				addHeaderLine(sb, "Diffusion_type", settings.getDiffusionType());
				addHeaderLine(sb, "Fixed_fraction", settings.fixedFraction);
				if (settings.compoundMolecules)
				{
					addHeaderLine(sb, "Compound_molecules", settings.compoundText.replaceAll("\n *", ""));
					addHeaderLine(sb, "Enable_2D_diffusion", settings.diffuse2D);
					addHeaderLine(sb, "Rotate_initial_orientation", settings.rotateInitialOrientation);
					addHeaderLine(sb, "Rotate_during_simulation", settings.rotateDuringSimulation);
					addHeaderLine(sb, "Enable_2D_rotation", settings.rotate2D);
				}
				addHeaderLine(sb, "Confinement", settings.confinement);
				if (settings.confinement.equals(CONFINEMENT[CONFINEMENT_SPHERE]))
				{
					addHeaderLine(sb, "Confinement_radius", settings.confinementRadius);
				}
				else if (settings.confinement.equals(CONFINEMENT[CONFINEMENT_MASK]))
				{
					addHeaderLine(sb, "Confinement_mask", settings.confinementMask);
				}
				addHeaderLine(sb, "Photon", settings.photonsPerSecond);
				addHeaderLine(sb, "Photon_distribution", settings.photonDistribution);
				if (PHOTON_DISTRIBUTION[PHOTON_CUSTOM].equals(settings.photonDistribution))
					addHeaderLine(sb, "Photon_distribution_file", settings.photonDistributionFile);
				else if (PHOTON_DISTRIBUTION[PHOTON_UNIFORM].equals(settings.photonDistribution))
					addHeaderLine(sb, "Photon_max", settings.photonsPerSecondMaximum);
				else if (PHOTON_DISTRIBUTION[PHOTON_GAMMA].equals(settings.photonDistribution))
					addHeaderLine(sb, "Photon_shape", settings.photonShape);
				else if (PHOTON_DISTRIBUTION[PHOTON_CORRELATED].equals(settings.photonDistribution))
					addHeaderLine(sb, "Correlation", settings.correlation);
				addHeaderLine(sb, "On_time", settings.tOn);
				if (!trackMode)
				{
					addHeaderLine(sb, "Off_time_short", settings.tOffShort);
					addHeaderLine(sb, "Off_time_long", settings.tOffLong);
					addHeaderLine(sb, "n_Blinks_short", settings.nBlinksShort);
					addHeaderLine(sb, "n_Blinks_long", settings.nBlinksLong);
					addHeaderLine(sb, "n_Blinks_Geometric", settings.nBlinksGeometricDistribution);
				}
				addHeaderLine(sb, "Min_photons", settings.minPhotons);
				addHeaderLine(sb, "Min_SNR_t1", settings.minSNRt1);
				addHeaderLine(sb, "Min_SNR_tN", settings.minSNRtN);
			}
			resultsFileHeader = sb.toString();
		}
		return resultsFileHeader;
	}

	private void addHeaderLine(StringBuffer sb, String name, Object o)
	{
		sb.append(String.format("# %-20s = %s\n", name, o.toString()));
	}

	/**
	 * Show a dialog allowing the parameters for a simple/benchmark simulation to be performed
	 * 
	 * @return True if the parameters were collected
	 */
	private boolean showSimpleDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);

		globalSettings = SettingsManager.loadSettings();
		settings = globalSettings.getCreateDataSettings();

		// Image size
		gd.addMessage("--- Image Size ---");
		gd.addNumericField("Pixel_pitch (nm)", settings.pixelPitch, 2);
		gd.addNumericField("Size (px)", settings.size, 0);
		if (!benchmarkMode)
		{
			gd.addNumericField("Depth (nm)", settings.depth, 0);
			gd.addCheckbox("Fixed_depth", settings.fixedDepth);
		}

		// Noise model
		gd.addMessage("--- Noise Model ---");
		if (extraOptions)
			gd.addCheckbox("No_poisson_noise", !settings.poissonNoise);
		gd.addNumericField("Background (photons)", settings.background, 2);
		gd.addNumericField("EM_gain", settings.getEmGain(), 2);
		gd.addNumericField("Camera_gain (ADU/e-)", settings.getCameraGain(), 4);
		gd.addNumericField("Quantum_efficiency", settings.getQuantumEfficiency(), 2);
		gd.addNumericField("Read_noise (e-)", settings.readNoise, 2);
		gd.addNumericField("Bias", settings.bias, 0);

		// PSF Model
		List<String> imageNames = addPSFOptions(gd);

		gd.addMessage("--- Fluorophores ---");
		Component splitLabel = gd.getMessage();
		// Do not allow grid or mask distribution
		if (simpleMode)
		{
			// Allow mask but not the grid
			gd.addChoice("Distribution", Arrays.copyOf(DISTRIBUTION, DISTRIBUTION.length - 1), settings.distribution);
			gd.addCheckbox("Sample_per_frame", settings.samplePerFrame);
		}
		gd.addNumericField("Particles", settings.particles, 0);
		if (simpleMode)
			gd.addNumericField("Density (um^-2)", settings.density, 2);
		else if (benchmarkMode)
		{
			gd.addNumericField("X_position (nm)", settings.xPosition, 2);
			gd.addNumericField("Y_position (nm)", settings.yPosition, 2);
			gd.addNumericField("Z_position (nm)", settings.zPosition, 2);
		}
		gd.addNumericField("Min_Photons", settings.photonsPerSecond, 0);
		gd.addNumericField("Max_Photons", settings.photonsPerSecondMaximum, 0);

		gd.addMessage("--- Save options ---");
		gd.addCheckbox("Raw_image", settings.rawImage);
		gd.addCheckbox("Save_image", settings.saveImage);
		gd.addCheckbox("Save_image_results", settings.saveImageResults);
		gd.addCheckbox("Save_localisations", settings.saveLocalisations);

		gd.addMessage("--- Report options ---");
		gd.addCheckbox("Show_histograms", settings.showHistograms);
		gd.addCheckbox("Choose_histograms", settings.chooseHistograms);
		gd.addNumericField("Histogram_bins", settings.histogramBins, 0);
		gd.addCheckbox("Remove_outliers", settings.removeOutliers);
		if (simpleMode)
			gd.addSlider("Density_radius (N x HWHM)", 0, 4.5, settings.densityRadius);
		gd.addNumericField("Depth-of-field (nm)", settings.depthOfField, 0);

		// Split into two columns
		// Re-arrange the standard layout which has a GridBagLayout with 2 columns (label,field)
		// to 4 columns: (label,field) x 2

		if (gd.getLayout() != null)
		{
			GridBagLayout grid = (GridBagLayout) gd.getLayout();

			int xOffset = 0, yOffset = 0;
			int lastY = -1, rowCount = 0;
			for (Component comp : gd.getComponents())
			{
				// Check if this should be the second major column
				if (comp == splitLabel)
				{
					xOffset += 2;
					yOffset -= rowCount;
					rowCount = 0;
				}
				// Reposition the field
				GridBagConstraints c = grid.getConstraints(comp);
				if (lastY != c.gridy)
					rowCount++;
				lastY = c.gridy;
				c.gridx = c.gridx + xOffset;
				c.gridy = c.gridy + yOffset;
				c.insets.left = c.insets.left + 10 * xOffset;
				c.insets.top = 0;
				c.insets.bottom = 0;
				grid.setConstraints(comp, c);
			}

			if (IJ.isLinux())
				gd.setBackground(new Color(238, 238, 238));
		}

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		settings.pixelPitch = Math.abs(gd.getNextNumber());
		settings.size = Math.abs((int) gd.getNextNumber());
		if (!benchmarkMode)
		{
			// Allow negative depth
			settings.depth = gd.getNextNumber();
			settings.fixedDepth = gd.getNextBoolean();
		}

		if (extraOptions)
			poissonNoise = settings.poissonNoise = !gd.getNextBoolean();
		settings.background = Math.abs(gd.getNextNumber());
		settings.setEmGain(Math.abs(gd.getNextNumber()));
		settings.setCameraGain(Math.abs(gd.getNextNumber()));
		settings.setQuantumEfficiency(Math.abs(gd.getNextNumber()));
		settings.readNoise = Math.abs(gd.getNextNumber());
		settings.bias = Math.abs((int) gd.getNextNumber());

		if (!collectPSFOptions(gd, imageNames))
			return false;

		if (simpleMode)
		{
			settings.distribution = gd.getNextChoice();
			settings.samplePerFrame = gd.getNextBoolean();
		}
		settings.particles = Math.abs((int) gd.getNextNumber());
		if (simpleMode)
			settings.density = Math.abs(gd.getNextNumber());
		else if (benchmarkMode)
		{
			settings.xPosition = gd.getNextNumber();
			settings.yPosition = gd.getNextNumber();
			settings.zPosition = gd.getNextNumber();
		}
		settings.photonsPerSecond = Math.abs((int) gd.getNextNumber());
		settings.photonsPerSecondMaximum = Math.abs((int) gd.getNextNumber());

		settings.rawImage = gd.getNextBoolean();
		settings.saveImage = gd.getNextBoolean();
		settings.saveImageResults = gd.getNextBoolean();
		settings.saveLocalisations = gd.getNextBoolean();

		settings.showHistograms = gd.getNextBoolean();
		settings.chooseHistograms = gd.getNextBoolean();
		settings.histogramBins = (int) Math.abs(gd.getNextNumber());
		settings.removeOutliers = gd.getNextBoolean();
		if (simpleMode)
			settings.densityRadius = (float) gd.getNextNumber();
		settings.depthOfField = (float) Math.abs(gd.getNextNumber());

		// Save before validation so that the current values are preserved.
		SettingsManager.saveSettings(globalSettings);

		if (gd.invalidNumber())
			return false;

		// Check arguments
		try
		{
			Parameters.isAboveZero("Pixel Pitch", settings.pixelPitch);
			Parameters.isAboveZero("Size", settings.size);
			if (!benchmarkMode && !settings.fixedDepth)
				Parameters.isPositive("Depth", settings.depth);
			Parameters.isPositive("Background", settings.background);
			Parameters.isPositive("EM gain", settings.getEmGain());
			Parameters.isPositive("Camera gain", settings.getCameraGain());
			Parameters.isPositive("Read noise", settings.readNoise);
			double noiseRange = settings.readNoise * settings.getCameraGain() * 4;
			Parameters.isEqualOrAbove("Bias must prevent clipping the read noise (@ +/- 4 StdDev) so ", settings.bias,
					noiseRange);
			Parameters.isAboveZero("Particles", settings.particles);
			if (simpleMode)
				Parameters.isAboveZero("Density", settings.density);
			Parameters.isAboveZero("Min Photons", settings.photonsPerSecond);
			if (settings.photonsPerSecondMaximum < settings.photonsPerSecond)
				settings.photonsPerSecondMaximum = settings.photonsPerSecond;
			if (!imagePSF)
			{
				Parameters.isAboveZero("Wavelength", settings.wavelength);
				Parameters.isAboveZero("NA", settings.numericalAperture);
				Parameters.isBelow("NA", settings.numericalAperture, 2);
			}
			Parameters.isPositive("Histogram bins", settings.histogramBins);
			if (simpleMode)
				Parameters.isPositive("Density radius", settings.densityRadius);
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		if (settings.distribution.equals(DISTRIBUTION[MASK]))
		{
			String[] maskImages = createDistributionImageList();
			if (maskImages != null)
			{
				gd = new GenericDialog(TITLE);
				gd.addMessage("Select the mask image for the distribution");
				gd.addChoice("Distribution_mask", maskImages, settings.distributionMask);
				if (maskListContainsStacks)
					gd.addNumericField("Distribution_slice_depth (nm)", settings.distributionMaskSliceDepth, 0);
				gd.showDialog();
				if (gd.wasCanceled())
					return false;
				settings.distributionMask = gd.getNextChoice();
				if (maskListContainsStacks)
					settings.distributionMaskSliceDepth = Math.abs(gd.getNextNumber());
			}

			SettingsManager.saveSettings(globalSettings);
		}

		return getHistogramOptions();
	}

	/**
	 * Check if there are any suitable PSF images open. If so add a choice to allow the selection of the Gaussian or
	 * Image PSF model. If no PSF images are open then add options for the wavelength and NA for the simulated
	 * microscope.
	 * 
	 * @param gd
	 * @return
	 */
	private List<String> addPSFOptions(GenericDialog gd)
	{
		gd.addMessage("--- PSF Model ---");
		List<String> imageNames = PSFCombiner.createImageList();
		String[] models = PSF_MODELS;
		if (imageNames.isEmpty())
		{
			imagePSF = false;
			models = Arrays.copyOf(PSF_MODELS, PSF_MODELS.length - 1);
		}
		gd.addChoice("PSF_model", models, settings.psfModel);
		gd.addCheckbox("Enter_width", settings.enterWidth);
		return imageNames;
	}

	/**
	 * If there are any suitable PSF images open then get the selected PSF model and collect the parameters required to
	 * configure it. If no PSF images are open then collect the wavelength and NA for the simulated microscope.
	 * 
	 * @param gd
	 * @return
	 */
	private boolean collectPSFOptions(GenericDialog gd, List<String> imageNames)
	{
		settings.psfModel = gd.getNextChoice();
		settings.enterWidth = gd.getNextBoolean();
		if (!imageNames.isEmpty())
			imagePSF = settings.psfModel.equals(PSF_MODELS[PSF_MODELS.length - 1]);

		// Show a second dialog to get the PSF parameters we need
		GenericDialog gd2 = new GenericDialog(TITLE);
		gd2.addMessage("Configure the " + settings.psfModel + " PSF model");
		if (imagePSF)
		{
			gd2.addChoice("PSF_image", imageNames.toArray(new String[imageNames.size()]), settings.psfImageName);
		}
		else
		{
			if (settings.enterWidth)
			{
				gd2.addNumericField("PSF_SD (nm)", settings.psfSD, 2);
			}
			else
			{
				gd2.addNumericField("Wavelength (nm)", settings.wavelength, 2);
				gd2.addNumericField("Numerical_aperture", settings.numericalAperture, 2);
			}
		}
		gd2.showDialog();
		if (gd2.wasCanceled())
			return false;
		if (imagePSF)
		{
			settings.psfImageName = gd2.getNextChoice();
		}
		else
		{
			if (settings.enterWidth)
			{
				settings.psfSD = Math.abs(gd2.getNextNumber());
			}
			else
			{
				settings.wavelength = Math.abs(gd2.getNextNumber());
				settings.numericalAperture = Math.abs(gd2.getNextNumber());
			}

			try
			{
				if (settings.enterWidth)
				{
					Parameters.isAboveZero("PSF SD", settings.psfSD);
				}
				else
				{
					Parameters.isAboveZero("Wavelength", settings.wavelength);
					Parameters.isAboveZero("NA", settings.numericalAperture);
				}
			}
			catch (IllegalArgumentException e)
			{
				IJ.error(TITLE, e.getMessage());
				return false;
			}
		}

		return true;
	}

	private boolean getHistogramOptions()
	{
		GenericDialog gd;
		if (settings.showHistograms && settings.chooseHistograms)
		{
			gd = new GenericDialog(TITLE);
			gd.addMessage("Select the histograms to display");
			for (int i = 0; i < displayHistograms.length; i++)
				gd.addCheckbox(NAMES[i].replace(' ', '_'), displayHistograms[i]);
			gd.showDialog();
			if (gd.wasCanceled())
				return false;
			for (int i = 0; i < displayHistograms.length; i++)
				displayHistograms[i] = gd.getNextBoolean();
		}
		return true;
	}

	/**
	 * Show a dialog allowing the parameters for a simulation to be performed
	 * 
	 * @return True if the parameters were collected
	 */
	private boolean showDialog()
	{
		// In track mode we do not need a time, illumination model or blinking model.
		// Fixed length tracks will be drawn, non-overlapping in time. This is the simplest
		// simulation for moving molecules

		GenericDialog gd = new GenericDialog(TITLE);

		globalSettings = SettingsManager.loadSettings();
		settings = globalSettings.getCreateDataSettings();

		if (settings.stepsPerSecond < 1)
			settings.stepsPerSecond = 1;

		String[] backgroundImages = createBackgroundImageList();
		gd.addNumericField("Pixel_pitch (nm)", settings.pixelPitch, 2);
		gd.addNumericField("Size (px)", settings.size, 0);
		gd.addNumericField("Depth (nm)", settings.depth, 0);
		gd.addCheckbox("Fixed_depth", settings.fixedDepth);
		if (!trackMode)
			gd.addNumericField("Seconds", settings.seconds, 1);
		gd.addNumericField("Exposure_time (ms)", settings.exposureTime, 1);
		gd.addSlider("Steps_per_second", 1, 15, settings.stepsPerSecond);
		if (!trackMode)
		{
			gd.addChoice("Illumination", ILLUMINATION, settings.illumination);
			gd.addNumericField("Pulse_interval", settings.pulseInterval, 0);
			gd.addNumericField("Pulse_ratio", settings.pulseRatio, 2);
		}
		if (backgroundImages != null)
			gd.addChoice("Background_image", backgroundImages, settings.backgroundImage);

		if (extraOptions)
			gd.addCheckbox("No_poisson_noise", !settings.poissonNoise);
		gd.addNumericField("Background (photons)", settings.background, 2);
		gd.addNumericField("EM_gain", settings.getEmGain(), 2);
		gd.addNumericField("Camera_gain (ADU/e-)", settings.getCameraGain(), 4);
		gd.addNumericField("Quantum_efficiency", settings.getQuantumEfficiency(), 2);
		gd.addNumericField("Read_noise (e-)", settings.readNoise, 2);
		gd.addNumericField("Bias", settings.bias, 0);

		List<String> imageNames = addPSFOptions(gd);

		gd.addMessage("--- Fluorophores ---");
		Component splitLabel = gd.getMessage();
		gd.addChoice("Distribution", DISTRIBUTION, settings.distribution);
		gd.addNumericField("Particles", settings.particles, 0);
		gd.addCheckbox("Compound_molecules", settings.compoundMolecules);
		gd.addNumericField("Diffusion_rate (um^2/sec)", settings.diffusionRate, 2);
		String[] diffusionTypes = SettingsManager.getNames((Object[]) DiffusionType.values());
		gd.addChoice("Diffusion_type", diffusionTypes, diffusionTypes[settings.getDiffusionType().ordinal()]);
		gd.addSlider("Fixed_fraction (%)", 0, 100, settings.fixedFraction * 100);
		gd.addChoice("Confinement", CONFINEMENT, settings.confinement);
		gd.addNumericField("Photons (sec^-1)", settings.photonsPerSecond, 0);
		// We cannot use the correlation moe with fixed life time tracks 
		String[] dist = (trackMode) ? Arrays.copyOf(PHOTON_DISTRIBUTION, PHOTON_DISTRIBUTION.length - 1)
				: PHOTON_DISTRIBUTION;
		gd.addChoice("Photon_distribution", dist, settings.photonDistribution);
		gd.addNumericField("On_time (ms)", settings.tOn, 2);
		if (!trackMode)
		{
			gd.addNumericField("Off_time_short (ms)", settings.tOffShort, 2);
			gd.addNumericField("Off_time_long (ms)", settings.tOffLong, 2);
			gd.addNumericField("n_Blinks_Short", settings.nBlinksShort, 2);
			gd.addNumericField("n_Blinks_Long", settings.nBlinksLong, 2);
			gd.addCheckbox("Use_geometric_distribution", settings.nBlinksGeometricDistribution);
		}

		gd.addMessage("--- Peak filtering ---");
		gd.addSlider("Min_Photons", 0, 50, settings.minPhotons);
		gd.addSlider("Min_SNR_t1", 0, 20, settings.minSNRt1);
		gd.addSlider("Min_SNR_tN", 0, 10, settings.minSNRtN);

		gd.addMessage("--- Save options ---");
		Component splitLabel2 = gd.getMessage();
		gd.addCheckbox("Raw_image", settings.rawImage);
		gd.addCheckbox("Save_image", settings.saveImage);
		gd.addCheckbox("Save_image_results", settings.saveImageResults);
		gd.addCheckbox("Save_fluorophores", settings.saveFluorophores);
		gd.addCheckbox("Save_localisations", settings.saveLocalisations);

		gd.addMessage("--- Report options ---");
		gd.addCheckbox("Show_histograms", settings.showHistograms);
		gd.addCheckbox("Choose_histograms", settings.chooseHistograms);
		gd.addNumericField("Histogram_bins", settings.histogramBins, 0);
		gd.addCheckbox("Remove_outliers", settings.removeOutliers);
		gd.addSlider("Density_radius (N x HWHM)", 0, 4.5, settings.densityRadius);
		gd.addNumericField("Depth-of-field (nm)", settings.depthOfField, 0);

		// Split into two columns
		// Re-arrange the standard layout which has a GridBagLayout with 2 columns (label,field)
		// to 4 columns: (label,field) x 2

		if (gd.getLayout() != null)
		{
			GridBagLayout grid = (GridBagLayout) gd.getLayout();

			int xOffset = 0, yOffset = 0;
			int lastY = -1, rowCount = 0;
			for (Component comp : gd.getComponents())
			{
				// Check if this should be the second major column
				if (comp == splitLabel || comp == splitLabel2)
				{
					xOffset += 2;
					yOffset -= rowCount;
					rowCount = 0;
				}
				// Reposition the field
				GridBagConstraints c = grid.getConstraints(comp);
				if (lastY != c.gridy)
					rowCount++;
				lastY = c.gridy;
				c.gridx = c.gridx + xOffset;
				c.gridy = c.gridy + yOffset;
				c.insets.left = c.insets.left + 10 * xOffset;
				c.insets.top = 0;
				c.insets.bottom = 0;
				grid.setConstraints(comp, c);
			}

			if (IJ.isLinux())
				gd.setBackground(new Color(238, 238, 238));
		}

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		settings.pixelPitch = Math.abs(gd.getNextNumber());
		settings.size = Math.abs((int) gd.getNextNumber());
		settings.depth = Math.abs(gd.getNextNumber());
		settings.fixedDepth = gd.getNextBoolean();
		if (!trackMode)
			settings.seconds = Math.abs(gd.getNextNumber());
		settings.exposureTime = Math.abs(gd.getNextNumber());
		settings.stepsPerSecond = Math.abs(gd.getNextNumber());
		if (!trackMode)
		{
			settings.illumination = gd.getNextChoice();
			settings.pulseInterval = Math.abs((int) gd.getNextNumber());
			settings.pulseRatio = Math.abs(gd.getNextNumber());
		}
		if (backgroundImages != null)
			settings.backgroundImage = gd.getNextChoice();

		if (extraOptions)
			poissonNoise = settings.poissonNoise = !gd.getNextBoolean();
		settings.background = Math.abs(gd.getNextNumber());
		settings.setEmGain(Math.abs(gd.getNextNumber()));
		settings.setCameraGain(Math.abs(gd.getNextNumber()));
		settings.setQuantumEfficiency(Math.abs(gd.getNextNumber()));
		settings.readNoise = Math.abs(gd.getNextNumber());
		settings.bias = Math.abs((int) gd.getNextNumber());

		if (!collectPSFOptions(gd, imageNames))
			return false;

		settings.distribution = gd.getNextChoice();
		settings.particles = Math.abs((int) gd.getNextNumber());
		settings.compoundMolecules = gd.getNextBoolean();
		settings.diffusionRate = Math.abs(gd.getNextNumber());
		settings.setDiffusionType(gd.getNextChoiceIndex());
		settings.fixedFraction = Math.abs(gd.getNextNumber() / 100.0);
		settings.confinement = gd.getNextChoice();
		settings.photonsPerSecond = Math.abs((int) gd.getNextNumber());
		settings.photonDistribution = gd.getNextChoice();
		settings.tOn = Math.abs(gd.getNextNumber());
		if (!trackMode)
		{
			settings.tOffShort = Math.abs(gd.getNextNumber());
			settings.tOffLong = Math.abs(gd.getNextNumber());
			settings.nBlinksShort = Math.abs(gd.getNextNumber());
			settings.nBlinksLong = Math.abs(gd.getNextNumber());
			settings.nBlinksGeometricDistribution = gd.getNextBoolean();
		}

		minPhotons = settings.minPhotons = gd.getNextNumber();
		minSNRt1 = settings.minSNRt1 = gd.getNextNumber();
		minSNRtN = settings.minSNRtN = gd.getNextNumber();

		settings.rawImage = gd.getNextBoolean();
		settings.saveImage = gd.getNextBoolean();
		settings.saveImageResults = gd.getNextBoolean();
		settings.saveFluorophores = gd.getNextBoolean();
		settings.saveLocalisations = gd.getNextBoolean();

		settings.showHistograms = gd.getNextBoolean();
		settings.chooseHistograms = gd.getNextBoolean();
		settings.histogramBins = (int) gd.getNextNumber();
		settings.removeOutliers = gd.getNextBoolean();
		settings.densityRadius = (float) gd.getNextNumber();
		settings.depthOfField = (float) Math.abs(gd.getNextNumber());

		// Ensure tN threshold is more lenient
		if (settings.minSNRt1 < settings.minSNRtN)
		{
			double tmp = settings.minSNRt1;
			settings.minSNRt1 = settings.minSNRtN;
			settings.minSNRtN = tmp;
		}

		// Save before validation so that the current values are preserved.
		SettingsManager.saveSettings(globalSettings);

		// Check arguments
		try
		{
			Parameters.isAboveZero("Pixel Pitch", settings.pixelPitch);
			Parameters.isAboveZero("Size", settings.size);
			if (!settings.fixedDepth)
				Parameters.isPositive("Depth", settings.depth);
			if (!trackMode)
				Parameters.isAboveZero("Seconds", settings.seconds);
			Parameters.isAboveZero("Exposure time", settings.exposureTime);
			Parameters.isAboveZero("Steps per second", settings.stepsPerSecond);
			Parameters.isPositive("Background", settings.background);
			Parameters.isPositive("EM gain", settings.getEmGain());
			Parameters.isPositive("Camera gain", settings.getCameraGain());
			Parameters.isPositive("Read noise", settings.readNoise);
			double noiseRange = settings.readNoise * settings.getCameraGain() * 4;
			Parameters.isEqualOrAbove("Bias must prevent clipping the read noise (@ +/- 4 StdDev) so ", settings.bias,
					noiseRange);
			Parameters.isAboveZero("Particles", settings.particles);
			Parameters.isAboveZero("Photons", settings.photonsPerSecond);
			if (!imagePSF)
			{
				Parameters.isAboveZero("Wavelength", settings.wavelength);
				Parameters.isAboveZero("NA", settings.numericalAperture);
				Parameters.isBelow("NA", settings.numericalAperture, 2);
			}
			Parameters.isPositive("Diffusion rate", settings.diffusionRate);
			Parameters.isPositive("Fixed fraction", settings.fixedFraction);
			Parameters.isPositive("Pulse interval", settings.pulseInterval);
			Parameters.isAboveZero("Pulse ratio", settings.pulseRatio);
			Parameters.isAboveZero("tOn", settings.tOn);
			if (!trackMode)
			{
				Parameters.isAboveZero("tOff Short", settings.tOffShort);
				Parameters.isAboveZero("tOff Long", settings.tOffLong);
				Parameters.isPositive("n-Blinks Short", settings.nBlinksShort);
				Parameters.isPositive("n-Blinks Long", settings.nBlinksLong);
			}
			Parameters.isPositive("Min photons", settings.minPhotons);
			Parameters.isPositive("Min SNR t1", settings.minSNRt1);
			Parameters.isPositive("Min SNR tN", settings.minSNRtN);
			Parameters.isAbove("Histogram bins", settings.histogramBins, 1);
			Parameters.isPositive("Density radius", settings.densityRadius);
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		if (gd.invalidNumber())
			return false;

		if (!getHistogramOptions())
			return false;

		String[] maskImages = null;
		if (settings.distribution.equals(DISTRIBUTION[MASK]))
		{
			maskImages = createDistributionImageList();
			if (maskImages != null)
			{
				gd = new GenericDialog(TITLE);
				gd.addMessage("Select the mask image for the distribution");
				gd.addChoice("Distribution_mask", maskImages, settings.distributionMask);
				if (maskListContainsStacks)
					gd.addNumericField("Distribution_slice_depth (nm)", settings.distributionMaskSliceDepth, 0);
				gd.showDialog();
				if (gd.wasCanceled())
					return false;
				settings.distributionMask = gd.getNextChoice();
				if (maskListContainsStacks)
					settings.distributionMaskSliceDepth = Math.abs(gd.getNextNumber());
			}
		}
		else if (settings.distribution.equals(DISTRIBUTION[GRID]))
		{
			gd = new GenericDialog(TITLE);
			gd.addMessage("Select grid for the distribution");
			gd.addNumericField("Cell_size", settings.cellSize, 0);
			gd.addSlider("p-binary", 0, 1, settings.probabilityBinary);
			gd.addNumericField("Min_binary_distance (nm)", settings.minBinaryDistance, 0);
			gd.addNumericField("Max_binary_distance (nm)", settings.maxBinaryDistance, 0);
			gd.showDialog();
			if (gd.wasCanceled())
				return false;
			settings.cellSize = (int) gd.getNextNumber();
			settings.probabilityBinary = gd.getNextNumber();
			settings.minBinaryDistance = gd.getNextNumber();
			settings.maxBinaryDistance = gd.getNextNumber();

			// Check arguments
			try
			{
				Parameters.isAboveZero("Cell size", settings.cellSize);
				Parameters.isPositive("p-binary", settings.probabilityBinary);
				Parameters.isEqualOrBelow("p-binary", settings.probabilityBinary, 1);
				Parameters.isPositive("Min binary distance", settings.minBinaryDistance);
				Parameters.isPositive("Max binary distance", settings.maxBinaryDistance);
				Parameters.isEqualOrBelow("Min binary distance", settings.minBinaryDistance,
						settings.maxBinaryDistance);
			}
			catch (IllegalArgumentException e)
			{
				IJ.error(TITLE, e.getMessage());
				return false;
			}
		}

		SettingsManager.saveSettings(globalSettings);

		if (settings.diffusionRate > 0 && settings.fixedFraction < 1)
		{
			if (settings.confinement.equals(CONFINEMENT[CONFINEMENT_SPHERE]))
			{
				gd = new GenericDialog(TITLE);
				gd.addMessage("Select the sphere radius for the diffusion confinement");
				gd.addSlider("Confinement_radius (nm)", 0, 2000, settings.confinementRadius);
				gd.showDialog();
				if (gd.wasCanceled())
					return false;
				settings.confinementRadius = gd.getNextNumber();
			}
			else if (settings.confinement.equals(CONFINEMENT[CONFINEMENT_MASK]))
			{
				if (maskImages == null)
					maskImages = createDistributionImageList();
				if (maskImages != null)
				{
					gd = new GenericDialog(TITLE);
					gd.addMessage("Select the mask image for the diffusion confinement");
					gd.addChoice("Confinement_mask", maskImages, settings.confinementMask);
					if (maskListContainsStacks)
						gd.addNumericField("Confinement_slice_depth (nm)", settings.confinementMaskSliceDepth, 0);
					gd.showDialog();
					if (gd.wasCanceled())
						return false;
					settings.confinementMask = gd.getNextChoice();
					if (maskListContainsStacks)
						settings.confinementMaskSliceDepth = Math.abs(gd.getNextNumber());
				}
			}
		}

		SettingsManager.saveSettings(globalSettings);

		if (settings.compoundMolecules)
		{
			// Show a second dialog where the molecule configuration is specified
			gd = new GenericDialog(TITLE);

			gd.addMessage("Specify the compound molecules");
			gd.addTextAreas(settings.compoundText, null, 20, 80);
			gd.addCheckbox("Enable_2D_diffusion", settings.diffuse2D);
			gd.addCheckbox("Rotate_initial_orientation", settings.rotateInitialOrientation);
			gd.addCheckbox("Rotate_during_simulation", settings.rotateDuringSimulation);
			gd.addCheckbox("Enable_2D_rotation", settings.rotate2D);
			gd.addCheckbox("Show_example_compounds", false);

			if (Utils.isShowGenericDialog())
			{
				@SuppressWarnings("rawtypes")
				Vector v = gd.getCheckboxes();
				Checkbox cb = (Checkbox) v.get(v.size() - 1);
				cb.addItemListener(this);
			}

			gd.showDialog();
			if (gd.wasCanceled())
				return false;

			settings.compoundText = gd.getNextText();
			settings.diffuse2D = gd.getNextBoolean();
			settings.rotateInitialOrientation = gd.getNextBoolean();
			settings.rotateDuringSimulation = gd.getNextBoolean();
			settings.rotate2D = gd.getNextBoolean();

			if (gd.getNextBoolean())
			{
				logExampleCompounds();
				return false;
			}
		}

		SettingsManager.saveSettings(globalSettings);

		gd = new GenericDialog(TITLE);
		gd.addMessage("Configure the photon distribution: " + settings.photonDistribution);
		if (PHOTON_DISTRIBUTION[PHOTON_CUSTOM].equals(settings.photonDistribution))
		{
			// Nothing more to be done
			return true;
		}
		else if (PHOTON_DISTRIBUTION[PHOTON_UNIFORM].equals(settings.photonDistribution))
		{
			gd.addNumericField("Max_Photons (sec^-1)", settings.photonsPerSecondMaximum, 0);
		}
		else if (PHOTON_DISTRIBUTION[PHOTON_GAMMA].equals(settings.photonDistribution))
		{
			gd.addNumericField("Photon_shape", settings.photonShape, 2);
		}
		else if (PHOTON_DISTRIBUTION[PHOTON_CORRELATED].equals(settings.photonDistribution))
		{
			gd.addNumericField("Correlation (to total tOn)", settings.correlation, 2);
		}
		else
		{
			// Nothing more to be done
			return true;
		}

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		try
		{
			if (PHOTON_DISTRIBUTION[PHOTON_UNIFORM].equals(settings.photonDistribution))
			{
				settings.photonsPerSecondMaximum = Math.abs((int) gd.getNextNumber());
				if (settings.photonsPerSecondMaximum < settings.photonsPerSecond)
					settings.photonsPerSecondMaximum = settings.photonsPerSecond;
			}
			else if (PHOTON_DISTRIBUTION[PHOTON_GAMMA].equals(settings.photonDistribution))
			{
				settings.photonShape = Math.abs(gd.getNextNumber());
				Parameters.isAbove("Photon shape", settings.photonShape, 0);
			}
			else if (PHOTON_DISTRIBUTION[PHOTON_CORRELATED].equals(settings.photonDistribution))
			{
				settings.correlation = gd.getNextNumber();
				Parameters.isEqualOrBelow("Correlation", settings.correlation, 1);
				Parameters.isEqualOrAbove("Correlation", settings.correlation, -1);
			}
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		SettingsManager.saveSettings(globalSettings);

		return true;
	}

	/**
	 * Build a list of suitable background images. Images must be greyscale.
	 * 
	 * @return
	 */
	private String[] createBackgroundImageList()
	{
		int[] idList = WindowManager.getIDList();
		if (idList != null)
		{
			String[] list = new String[idList.length + 1];
			list[0] = "[None]";
			int count = 1;
			for (int id : idList)
			{
				ImagePlus imp = WindowManager.getImage(id);
				// Image must be square and greyscale
				if (imp != null && imp.getWidth() == imp.getHeight() &&
						(imp.getBitDepth() == 8 || imp.getBitDepth() == 16 || imp.getBitDepth() == 32))
				{
					list[count++] = imp.getTitle();
				}
			}
			if (count == 1)
				return null;
			return Arrays.copyOf(list, count);
		}
		return null;
	}

	/**
	 * Build a list of suitable distribution images. Images must be square.
	 * 
	 * @return
	 */
	private String[] createDistributionImageList()
	{
		maskListContainsStacks = false;
		int[] idList = WindowManager.getIDList();
		if (idList != null)
		{
			String[] list = new String[idList.length + 1];
			list[0] = "[None]";
			int count = 1;
			for (int id : idList)
			{
				ImagePlus imp = WindowManager.getImage(id);
				if (imp != null && imp.getWidth() == imp.getHeight())
				{
					list[count++] = imp.getTitle();
					if (imp.getStackSize() > 1)
						maskListContainsStacks = true;
				}
			}
			if (count == 1)
				return null;
			return Arrays.copyOf(list, count);
		}
		return null;
	}

	public void itemStateChanged(ItemEvent e)
	{
		// When the checkbox is clicked, output example compounds to the ImageJ log
		Checkbox cb = (Checkbox) e.getSource();
		if (cb.getState())
		{
			cb.setState(false);

			logExampleCompounds();
		}
	}

	private void logExampleCompounds()
	{
		comment(TITLE + " example compounds");
		IJ.log("");
		comment("Compounds are described using XML");
		comment("Multiple compounds can be combined using fractional ratios");
		comment("Coordinates are specified in nanometres");
		comment("Coordinates describe the relative positions of molecules in the compound. Compounds will have a randomly assigned XYZ position for their centre-of-mass. Rotation will be about the centre of mass");
		IJ.log("");

		Atom a1 = new Atom(10, 0, 0, 0);
		Atom a2 = new Atom(30, 0, 0, 0);
		Atom a3 = new Atom(20, 1000, 0, 0);
		Compound m1 = new Compound(1, 0, DiffusionType.RANDOM_WALK.toString(), a1);
		Compound m2 = new Compound(1, 1, DiffusionType.GRID_WALK.toString(), a2, a3);

		// Create a hexamer big enough to see with the default pixel pitch
		Atom b1 = new Atom(1, 0, 0, 0);
		Atom b2 = new Atom(1, 1000, 0, 0);
		Atom b3 = new Atom(1, 1500, 866, 0);
		Atom b4 = new Atom(1, 1000, 1732, 0);
		Atom b5 = new Atom(1, 0, 1732, 0);
		Atom b6 = new Atom(1, -500, 866, 0);
		Compound m3 = new Compound(1, 2, DiffusionType.RANDOM_WALK.toString(), b1, b2, b3, b4, b5, b6);

		comment("Single compounds");
		IJ.log("");
		comment("Monomer");
		demo(m1);
		comment("Dimer");
		demo(m2);
		comment("Hexamer");
		demo(m3);

		comment("Combined compounds");
		IJ.log("");
		comment("Two compounds with a ratio of 2:1");
		m1.fraction = 2;
		demo(m1, m2);
	}

	private void demo(Compound... compounds)
	{
		List<Compound> list = new LinkedList<Compound>();
		for (Compound c : compounds)
			list.add(c);

		IJ.log(createXStream().toXML(list));
		IJ.log("");
	}

	private XStream createXStream()
	{
		if (xs == null)
		{
			xs = new XStream(new DomDriver());
			xs.autodetectAnnotations(true);
			xs.alias("Compound", Compound.class);
			xs.alias("Atom", Atom.class);
		}
		return xs;
	}

	private void comment(String text)
	{
		IJ.log(TextUtils.wrap("<!-- " + text + " -->", 80));
	}

	@SuppressWarnings("unchecked")
	private List<CompoundMoleculeModel> createCompoundMolecules()
	{
		// Diffusion rate is um^2/sec. Convert to pixels per simulation frame.
		final double diffusionFactor = (1000000.0 / (settings.pixelPitch * settings.pixelPitch)) /
				settings.stepsPerSecond;

		List<CompoundMoleculeModel> compounds;
		if (settings.compoundMolecules)
		{
			// Try and load the compounds from the XML specification
			try
			{
				Object fromXML = createXStream().fromXML(settings.compoundText);
				List<Compound> rawCompounds = (List<Compound>) fromXML;

				// Convert from the XML serialised objects to the compound model

				compounds = new ArrayList<CompoundMoleculeModel>(rawCompounds.size());
				int id = 1;
				for (Compound c : rawCompounds)
				{
					MoleculeModel[] molecules = new MoleculeModel[c.atoms.length];
					for (int i = 0; i < c.atoms.length; i++)
					{
						Atom a = c.atoms[i];
						molecules[i] = new MoleculeModel(a.mass, a.x, a.y, a.z);
					}
					CompoundMoleculeModel m = new CompoundMoleculeModel(id++, 0, 0, 0, Arrays.asList(molecules));
					m.setFraction(c.fraction);
					m.setDiffusionRate(c.D * diffusionFactor);
					m.setDiffusionType(DiffusionType.fromString(c.diffusionType));
					compounds.add(m);
				}

				// Convert coordinates from nm to pixels
				final double scaleFactor = 1.0 / settings.pixelPitch;
				for (CompoundMoleculeModel c : compounds)
				{
					c.scale(scaleFactor);
				}
			}
			catch (Exception e)
			{
				IJ.error(TITLE, "Unable to create compound molecules");
				return null;
			}
		}
		else
		{
			// Create a simple compound with one molecule at the origin
			compounds = new ArrayList<CompoundMoleculeModel>(1);
			CompoundMoleculeModel m = new CompoundMoleculeModel(1, 0, 0, 0,
					Arrays.asList(new MoleculeModel(0, 0, 0, 0)));
			m.setDiffusionRate(settings.diffusionRate * diffusionFactor);
			m.setDiffusionType(settings.getDiffusionType());
			compounds.add(m);
		}
		return compounds;
	}

	private static int seedAddition = 0;
	private boolean resetSeed = true;

	private enum SeedMode
	{
		//@formatter:off
		DEFAULT{ 
			@Override boolean identicalOffset() { return false; }
			@Override boolean identicalAddition() { return false; } },
		REPRODUCE_EACH_STARTUP{ 
				@Override boolean identicalOffset() { return true; }
				@Override boolean identicalAddition() { return false; } },
		REPRODUCE_EACH_RUN{ 
			@Override boolean identicalOffset() { return true; }
			@Override boolean identicalAddition() { return true; } };
		//@formatter:on

		/**
		 * Set to true when the seed should not have a time dependent offset. This ensure that the plugin can be run at
		 * any time from start-up and the simulation will be reproducible.
		 *
		 * @return true, if the seed offset should be identical
		 */
		abstract boolean identicalOffset();

		/**
		 * Set to true when the seed should have the same addition for each plugin execution. Set to false will ensure
		 * subsequent runs of the plugin produce different simulations.
		 *
		 * @return true, if the seed addition should be reset for each plugin execution
		 */
		abstract boolean identicalAddition();
	}

	private SeedMode seedMode = SeedMode.DEFAULT;

	private long getSeedOffset()
	{
		return (seedMode.identicalOffset()) ? 1 : System.currentTimeMillis() + System.identityHashCode(this);
	}

	private int getSeedAddition()
	{
		if (seedMode.identicalAddition())
		{
			if (resetSeed)
			{
				// Reset only once per plugin execution
				resetSeed = false;
				seedAddition = 0;
			}
		}
		// Increment the seed to ensure that new generators are created at the same system time point
		return seedAddition++;
	}

	/**
	 * Get a random generator. The generators used in the simulation can be adjusted by changing this method.
	 * 
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.RandomGeneratorFactory#createRandomGenerator()
	 */
	public RandomGenerator createRandomGenerator()
	{
		return new Well44497b(getSeedOffset() + getSeedAddition());
	}

	private static void setBenchmarkResults(ImagePlus imp, MemoryPeakResults results)
	{
		if (imp != null)
		{
			benchmarkImageId = imp.getID();
			benchmarkResultsName = results.getName();
			MemoryPeakResults.addResults(results);
		}
		else
		{
			benchmarkImageId = 0;
			benchmarkResultsName = "";
		}
	}

	/**
	 * Gets the benchmark image.
	 *
	 * @return the image
	 */
	public static ImagePlus getImage()
	{
		return WindowManager.getImage(benchmarkImageId);
	}

	/**
	 * Gets the benchmark image Id.
	 *
	 * @return the image Id
	 */
	public static int getImageId()
	{
		return benchmarkImageId;
	}

	/**
	 * Gets the benchmark results.
	 *
	 * @return the results
	 */
	public static MemoryPeakResults getResults()
	{
		return MemoryPeakResults.getResults(benchmarkResultsName);
	}

	/**
	 * Load benchmark data using an open image and a XYZ text file.
	 */
	private void loadBenchmarkData()
	{
		// Note: Do not reset memory until failure. This allows the load method to use the 
		// last simulation parameters to set settings.

		if (!showLoadDialog())
		{
			//resetMemory();
			return;
		}

		// Load the image
		ImagePlus imp = WindowManager.getImage(benchmarkImage);
		if (imp == null)
		{
			IJ.error(TITLE, "No benchmark image: " + benchmarkImage);
			//resetMemory();
			return;
		}

		// Load the results
		MemoryPeakResults results = getSimulationResults();
		if (results == null)
		{
			IJ.error(TITLE, "No benchmark results: " + benchmarkResultsName);
			//resetMemory();
			return;
		}
		results.setName(imp.getTitle() + " (Results)");
		results.setBounds(new Rectangle(0, 0, imp.getWidth(), imp.getHeight()));
		results.setSource(new IJImageSource(imp));

		// Get the calibration
		simulationParameters = showSimulationParametersDialog(imp, results);
		if (simulationParameters != null)
		{
			setBackground(results);
			setNoise(results, imp);
			setBenchmarkResults(imp, results);
			IJ.showStatus("Loaded " + Utils.pleural(results.size(), "result"));
		}
		else
		{
			resetMemory();
		}
	}

	/**
	 * Sets the background in the results if missing.
	 *
	 * @param results
	 *            the results
	 */
	private void setBackground(MemoryPeakResults results)
	{
		// Loaded results do not have a local background.
		for (PeakResult p : results.getResults())
			if (p.getBackground() != 0)
				return;

		// Simple fix is to use the bias plus the global photon background.
		// TODO - Subtract the spots from the local region and compute the true local background.
		// Note this requires knowing the PSF width. If this is a loaded ground truth dataset then 
		// it probably will not have Gaussian widths.
		final float b = (float) (simulationParameters.bias + simulationParameters.gain * simulationParameters.b);
		for (PeakResult p : results.getResults())
			p.params[Gaussian2DFunction.BACKGROUND] = b;
	}

	/**
	 * Sets the noise in the results if missing.
	 *
	 * @param results
	 *            the results
	 */
	private void setNoise(MemoryPeakResults results, ImagePlus imp)
	{
		// Loaded results do not have noise
		for (PeakResult r : results.getResults())
			if (r.noise != 0)
				return;

		// Compute noise per frame
		ImageStack stack = imp.getImageStack();
		final int width = stack.getWidth();
		final int height = stack.getHeight();
		final IJImageSource source = new IJImageSource(imp);
		final float[] noise = new float[source.getFrames() + 1];
		for (int slice = 1; slice < noise.length; slice++)
		{
			stack.getPixels(slice);
			float[] data = source.next();
			// Use the trimmed method as there may be a lot of spots in the frame
			noise[slice] = (float) FitWorker.estimateNoise(data, width, height,
					NoiseEstimator.Method.QUICK_RESIDUALS_LEAST_TRIMMED_OF_SQUARES);
		}

		Statistics stats = new Statistics(Arrays.copyOfRange(noise, 1, noise.length));
		System.out.printf("Noise = %.3f +/- %.3f (%d)\n", stats.getMean(), stats.getStandardDeviation(), stats.getN());

		for (PeakResult p : results.getResults())
		{
			if (p.getFrame() < noise.length)
				p.noise = noise[p.getFrame()];
		}
	}

	private MemoryPeakResults getSimulationResults()
	{
		if (benchmarkAuto)
		{
			// Load directly from a results file. This is mainly to be used to load simulations 
			// saved to memory then saved to file. This is because the z-depth must be in the 
			// error field of the results.
			PeakResultsReader r = new PeakResultsReader(benchmarkFile);
			MemoryPeakResults results = r.getResults();
			if (results != null)
			{
				ResultsManager.checkCalibration(results);
				return results;
			}
		}

		// Load using a universal text file
		LocalisationList localisations = LoadLocalisations.loadLocalisations(benchmarkFile);
		if (localisations.isEmpty())
			return null;

		return localisations.toPeakResults();
	}

	private SimulationParameters showSimulationParametersDialog(ImagePlus imp, MemoryPeakResults results)
	{
		int molecules = results.size();

		// Get the missing parameters from the user
		boolean fullSimulation = false;
		double s = -1;

		// Get these from the data
		double[] signal = getSignal(results);
		double[] limits = Maths.limits(signal);
		double minSignal = limits[0];
		double maxSignal = limits[1];
		double signalPerFrame = Maths.sum(signal) / molecules;

		double[] depths = getDepth(results);
		limits = Maths.limits(depths);
		double depth = Math.max(Math.abs(limits[0]), Math.abs(limits[1]));
		boolean fixedDepth = Double.compare(limits[0], limits[1]) == 0;

		Calibration cal = new Calibration();
		// Get any calibration we have
		if (results.getCalibration() != null)
			cal = results.getCalibration();

		// Get this from the user
		double b = -1;

		// Use last simulation parameters for missing settings.
		// This is good if we are re-running the plugin to load data.
		if (simulationParameters != null && simulationParameters.isLoaded())
		{
			fullSimulation = simulationParameters.fullSimulation;
			s = simulationParameters.s;
			b = simulationParameters.b;
			if (!cal.hasBias())
				cal.setBias(simulationParameters.bias);
			if (!cal.hasGain())
				cal.setGain(simulationParameters.gain);
			if (!cal.hasAmplification())
				cal.setAmplification(simulationParameters.amplification);
			if (!cal.hasReadNoise())
				cal.setReadNoise(simulationParameters.readNoise);
			if (!cal.hasEMCCD())
				cal.setEmCCD(simulationParameters.emCCD);
			if (!cal.hasNmPerPixel())
				cal.setNmPerPixel(simulationParameters.a);
		}

		// Show a dialog to confirm settings
		GenericDialog gd = new GenericDialog(TITLE);
		StringBuilder sb = new StringBuilder();
		sb.append("Results contain ").append(Utils.pleural(molecules, "molecule")).append('\n');
		sb.append("Min signal = ").append(Utils.rounded(minSignal)).append(" ADU\n");
		sb.append("Max signal = ").append(Utils.rounded(maxSignal)).append(" ADU\n");
		sb.append("Av signal = ").append(Utils.rounded(signalPerFrame)).append(" ADU\n");
		if (fixedDepth)
			sb.append("Fixed depth = ").append(Utils.rounded(depth)).append('\n');
		gd.addMessage(sb.toString());

		gd.addCheckbox("Flourophore_simulation", fullSimulation);
		gd.addNumericField("Gaussian_SD", s, 3, 8, "nm");
		gd.addNumericField("Pixel_pitch", cal.getNmPerPixel(), 3, 8, "nm");
		gd.addNumericField("Background", b, 3, 8, "photon");
		gd.addNumericField("Total_gain", cal.getGain(), 3, 8, "ADU/photon");
		gd.addNumericField("Amplification", cal.getAmplification(), 3, 8, "ADU/e-");
		gd.addCheckbox("EM-CCD", cal.isEmCCD());
		gd.addNumericField("Read_noise", cal.getReadNoise(), 3, 8, "ADU");
		gd.addNumericField("Bias", cal.getBias(), 3, 8, "pixel");

		if (!fixedDepth)
		{
			gd.addNumericField("Depth", depth, 3, 8, "pixel");
		}

		gd.showDialog();
		if (gd.wasCanceled())
			return null;

		fullSimulation = gd.getNextBoolean();
		s = gd.getNextNumber();
		cal.setNmPerPixel(gd.getNextNumber());
		b = gd.getNextNumber();
		cal.setGain(gd.getNextNumber());
		cal.setAmplification(gd.getNextNumber());
		cal.setEmCCD(gd.getNextBoolean());
		cal.setReadNoise(gd.getNextNumber());
		cal.setBias(gd.getNextNumber());
		double myDepth = depth;
		if (!fixedDepth)
		{
			myDepth = gd.getNextNumber();
			if (myDepth < depth)
			{
				IJ.error(TITLE, String.format("Input depth is smaller than the depth guessed from the data: %f < %f",
						myDepth, depth));
				return null;
			}
			depth = myDepth;
		}

		// Validate settings
		try
		{
			Parameters.isAboveZero("Gaussian_SD", s);
			Parameters.isAboveZero("Pixel_pitch", cal.getNmPerPixel());
			Parameters.isPositive("Background", b);
			Parameters.isAboveZero("Total_gain", cal.getGain());
			Parameters.isAboveZero("Amplification", cal.getAmplification());
			Parameters.isPositive("Read_noise", cal.getReadNoise());
			Parameters.isPositive("Bias", cal.getBias());
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return null;
		}

		// Store calibration
		results.setCalibration(cal);

		double gain = cal.getGain();
		double a = cal.getNmPerPixel();
		double bias = cal.getBias();
		double readNoise = cal.getReadNoise();
		double amplification = cal.getAmplification();
		boolean emCCD = cal.isEmCCD();
		
		// Convert ADU values to photons
		minSignal /= gain;
		maxSignal /= gain;
		signalPerFrame /= gain;

		// Convert +/- depth to total depth in nm
		depth *= 2 * a;

		// Compute total background variance in photons
		double backgroundVariance = b;
		// Do not add EM-CCD noise factor. The Mortensen formula also includes this factor 
		// so this is "double-counting" the EM-CCD.  
		//if (emCCD)
		//	backgroundVariance *= 2;

		// Read noise is in ADUs. Convert to Photons to get contribution to background variance
		double readNoiseInPhotons = readNoise / gain;

		// Get the expected value at each pixel in photons. Assuming a Poisson distribution this 
		// is equal to the total variance at the pixel.
		double b2 = backgroundVariance + readNoiseInPhotons * readNoiseInPhotons;

		SimulationParameters p = new SimulationParameters(molecules, fullSimulation, s, a, minSignal, maxSignal,
				signalPerFrame, depth, fixedDepth, bias, emCCD, gain, amplification, readNoise, b, b2);
		p.loaded = true;
		return p;
	}

	private static double[] getSignal(MemoryPeakResults results)
	{
		double[] data = new double[results.size()];
		int i = 0;
		for (PeakResult p : results.getResults())
			data[i++] = p.getSignal();
		return data;
	}

	private static double[] getDepth(MemoryPeakResults results)
	{
		double[] data = new double[results.size()];
		int i = 0;
		for (PeakResult p : results.getResults())
			// Results store the z-depth in the error field
			data[i++] = p.error;
		return data;
	}

	/**
	 * Show a dialog allowing the parameters for a benchmark simulation to be loaded
	 * 
	 * @return True if the parameters were collected
	 */
	private boolean showLoadDialog()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

		String[] images = Utils.getImageList(Utils.GREY_SCALE);
		gd.addChoice("Image", images, benchmarkImage);
		gd.addFilenameField("Results_file", benchmarkFile);
		gd.addCheckbox("Preprocessed_results", benchmarkAuto);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		benchmarkImage = gd.getNextChoice();
		benchmarkFile = gd.getNextString();
		benchmarkAuto = gd.getNextBoolean();

		return true;
	}
}
