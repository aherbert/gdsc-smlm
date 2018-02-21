package gdsc.smlm.ij.plugins;

import java.awt.Checkbox;
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
import java.util.Map;
import java.util.Map.Entry;
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

import com.google.protobuf.TextFormat;

import gdsc.core.clustering.DensityManager;
import gdsc.core.data.DataException;
import gdsc.core.data.utils.ConversionException;
import gdsc.core.data.utils.TypeConverter;
import gdsc.core.ij.Utils;
import gdsc.core.threshold.AutoThreshold;
import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.Maths;
import gdsc.core.utils.Random;
import gdsc.core.utils.Statistics;
import gdsc.core.utils.StoredDataStatistics;
import gdsc.core.utils.TextUtils;
import gdsc.core.utils.TurboList;
import gdsc.core.utils.UnicodeReader;
import gdsc.smlm.data.config.CalibrationProtos.CameraType;
import gdsc.smlm.data.config.CalibrationProtosHelper;
import gdsc.smlm.data.config.CalibrationWriter;
import gdsc.smlm.data.config.ConfigurationException;
import gdsc.smlm.data.config.CreateDataSettingsHelper;
import gdsc.smlm.data.config.FitProtos.NoiseEstimatorMethod;
import gdsc.smlm.data.config.GUIProtos.CreateDataSettings;
import gdsc.smlm.data.config.GUIProtos.LoadLocalisationsSettings;
import gdsc.smlm.data.config.MoleculeProtos.Atom;
import gdsc.smlm.data.config.MoleculeProtos.AtomOrBuilder;
import gdsc.smlm.data.config.MoleculeProtos.Mixture;
import gdsc.smlm.data.config.MoleculeProtos.Molecule;
import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.data.config.PSFProtos.AstigmatismModel;
import gdsc.smlm.data.config.PSFProtos.ImagePSF;
import gdsc.smlm.data.config.PSFProtos.Offset;
import gdsc.smlm.data.config.PSFProtos.PSF;
import gdsc.smlm.data.config.PSFProtosHelper;
import gdsc.smlm.data.config.UnitHelper;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import gdsc.smlm.engine.FitWorker;
import gdsc.smlm.filters.GaussianFilter;
import gdsc.smlm.function.gaussian.AstigmatismZModel;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.HoltzerAstigmatismZModel;
import gdsc.smlm.ij.IJImageSource;
import gdsc.smlm.ij.plugins.LoadLocalisations.LocalisationList;
import gdsc.smlm.ij.settings.ImagePSFHelper;
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
import gdsc.smlm.model.camera.CameraModel;
import gdsc.smlm.model.camera.FixedPixelCameraModel;
import gdsc.smlm.results.ExtendedPeakResult;
import gdsc.smlm.results.Gaussian2DPeakResultHelper;
import gdsc.smlm.results.IdPeakResult;
import gdsc.smlm.results.ImmutableMemoryPeakResults;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResults;
import gdsc.smlm.results.PeakResultsReader;
import gdsc.smlm.results.SynchronizedPeakResults;
import gdsc.smlm.results.TextFilePeakResults;
import gdsc.smlm.results.count.FrameCounter;
import gdsc.smlm.results.procedures.PeakResultProcedure;
import gdsc.smlm.results.procedures.PrecisionResultProcedure;
import gdsc.smlm.results.procedures.RawResultProcedure;
import gdsc.smlm.results.procedures.StandardResultProcedure;
import gdsc.smlm.results.procedures.WidthResultProcedure;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.set.hash.TIntHashSet;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.ExtendedGenericDialog;
import ij.gui.ExtendedGenericDialog.OptionListener;
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

	private static String[] PSF_MODELS = new String[] { "2D Gaussian", "Airy", "Image", "Astigmatism" };
	private static final int PSF_MODEL_GAUSSIAN = 0;
	private static final int PSF_MODEL_AIRY = 1;
	private static final int PSF_MODEL_IMAGE = 2;
	private static final int PSF_MODEL_ASTIGMATISM = 3;

	/**
	 * The PSF model type. This is set when validating the PSF settings.
	 */
	private int psfModelType = -1;
	private AstigmatismModel astigmatismModel = null;

	private static TextWindow summaryTable = null;
	private static int datasetNumber = 0;
	private static double areaInUm = 0;
	private static String header = null;

	private CreateDataSettings.Builder settings;

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
	//private boolean imagePSF;
	private double hwhm = 0;

	private TIntHashSet movingMolecules;
	private TIntIntHashMap idToCompound;
	ArrayList<String> compoundNames;
	private boolean maskListContainsStacks;

	// Created by drawImage(...)
	private MemoryPeakResults results = null;

	// Used by the ImageGenerator to show progress when the thread starts
	private int frame, maxT, totalFrames;

	private boolean simpleMode = false;
	private boolean benchmarkMode = false;
	private boolean spotMode = false;
	private boolean trackMode = false;
	private boolean extraOptions = false;

	// Hold private variables for settings that are ignored in simple/benchmark mode 
	private boolean poissonNoise = true;
	private double minPhotons = 0, minSNRt1 = 0, minSNRtN = 0;

	// Store the parameters
	public static class BaseParameters
	{
		private static int nextId = 1;

		/**
		 * The parameter set identifier
		 */
		final int id;
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
		double averageSignal;
		/**
		 * The camera bias
		 */
		final double bias;
		/**
		 * Total gain (ADUs/photon)
		 */
		final double gain;
		/**
		 * Quantum efficiency (electron/photon)
		 */
		final double qe;
		/**
		 * Read noise in ADUs
		 */
		final double readNoise;
		/**
		 * The camera type.
		 */
		CameraType cameraType;
		/**
		 * The camera model name.
		 */
		String cameraModelName;
		/**
		 * The camera bounds.
		 */
		Rectangle cameraBounds;
		/**
		 * Background
		 */
		double b;
		/**
		 * Background noise in photons per pixel (used in the precision calculations)
		 */
		final double noise;

		public BaseParameters(double s, double a, double minSignal, double maxSignal, double averageSignal, double bias,
				double gain, double qe, double readNoise, CameraType cameraType, String cameraModelName,
				Rectangle cameraBounds, double b, double noise)
		{
			id = nextId++;
			this.s = s;
			this.a = a;
			this.minSignal = minSignal;
			this.maxSignal = maxSignal;
			this.averageSignal = averageSignal;
			this.bias = bias;
			this.gain = gain;
			this.qe = qe;
			this.readNoise = readNoise;
			this.cameraType = cameraType;
			this.cameraModelName = cameraModelName;
			this.cameraBounds = cameraBounds;
			this.b = b;
			this.noise = noise;
		}

		/**
		 * Checks if is emccd.
		 *
		 * @return true, if is emccd
		 */
		public boolean isEMCCD()
		{
			return cameraType == CameraType.EMCCD;
		}
	}

	// Store the parameters for the last simulation for spot data
	public static class SimulationParameters extends BaseParameters
	{
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
		 * The z-position depth
		 */
		final double depth;
		/**
		 * True if the depth is fixed
		 */
		final boolean fixedDepth;

		private boolean loaded;

		public SimulationParameters(int molecules, boolean fullSimulation, double s, double a, double minSignal,
				double maxSignal, double averageSignal, double depth, boolean fixedDepth, double bias, double gain,
				double qe, double readNoise, CameraType cameraType, String cameraModelName, Rectangle cameraBounds,
				double b, double noise)
		{
			super(s, a, minSignal, maxSignal, averageSignal, bias, gain, qe, readNoise, cameraType, cameraModelName,
					cameraBounds, b, noise);
			this.molecules = molecules;
			this.fullSimulation = fullSimulation;
			this.depth = depth;
			// We must have a fixed depth if the depth is zero
			this.fixedDepth = (depth > 0) ? fixedDepth : true;
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
	public static class BenchmarkParameters extends BaseParameters
	{
		/**
		 * Number of frames in the simulated image
		 */
		final int frames;
		/**
		 * The x,y,z position of the localisation in each frame
		 */
		final double x, y, z;
		/**
		 * The actual number of simulated photons in each frame of the benchmark image. Some frames may be empty (due to
		 * signal filtering or Poisson sampling).
		 */
		final double[] p;
		/**
		 * The actual number of simulated background photons in each frame of the benchmark image.
		 */
		final double[] background;
		/**
		 * The number of frames with a simulated photon count
		 */
		private int molecules;

		final double precisionN, precisionX, precisionXML;

		public BenchmarkParameters(int frames, double s, double a, double signal, double x, double y, double z,
				double bias, double gain, double qe, double readNoise, CameraType cameraType, String cameraModelName,
				Rectangle cameraBounds, double b, double noise, double precisionN, double precisionX,
				double precisionXML)
		{
			super(s, a, signal, signal, signal, bias, gain, qe, readNoise, cameraType, cameraModelName, cameraBounds, b,
					noise);

			this.frames = frames;
			this.x = x;
			this.y = y;
			this.z = z;
			this.precisionN = precisionN;
			this.precisionX = precisionX;
			this.precisionXML = precisionXML;
			p = new double[frames];
			background = new double[frames];
		}

		public void setPhotons(MemoryPeakResults results)
		{
			molecules = results.size();
			results.forEach(new PeakResultProcedure()
			{
				public void execute(PeakResult result)
				{
					int i = result.getFrame() - 1;
					if (p[i] != 0)
						throw new RuntimeException("Multiple peaks on the same frame: " + result.getFrame());
					p[i] = result.getSignal();
					background[i] = result.getBackground();
				}
			});
			double av = Maths.sum(p) / molecules;
			double av2 = Maths.sum(background) / molecules;
			Utils.log(
					"Created %d frames, %d molecules. Simulated signal %s : average %s. Simulated background %s : average %s",
					frames, molecules, Utils.rounded(averageSignal), Utils.rounded(av / gain), Utils.rounded(b),
					Utils.rounded(av2 / gain));
			// Reset the average signal and background (in photons)
			averageSignal = av / gain;
			b = av2 / gain;
		}

		/**
		 * @return The average number of photons per frame
		 */
		public double getSignal()
		{
			return averageSignal;
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
	private static boolean benchmarkAuto = false;
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

			settings.setExposureTime(1000); // 1 second frames
			areaInUm = settings.getSize() * settings.getPixelPitch() * settings.getSize() * settings.getPixelPitch() /
					1e6;

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
				nextN = new int[settings.getParticles() * 2];
				Arrays.fill(nextN, 0, settings.getParticles(), 1);
				Random rand = new Random();
				rand.shuffle(nextN);

				// Only put spots in the central part of the image
				double border = settings.getSize() / 4.0;
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
				if (settings.getSamplePerFrame())
				{
					final double mean = areaInUm * settings.getDensity();
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
					StoredDataStatistics samples = new StoredDataStatistics(settings.getParticles());
					while (samples.getSum() < settings.getParticles())
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
					n = (int) FastMath.max(1, Math.round(areaInUm * settings.getDensity()));
				}
			}

			RandomGenerator random = null;

			localisations = new ArrayList<LocalisationModel>(settings.getParticles());
			localisationSets = new ArrayList<LocalisationModelSet>(settings.getParticles());

			final int minPhotons = (int) settings.getPhotonsPerSecond();
			final int range = (int) settings.getPhotonsPerSecondMaximum() - minPhotons + 1;
			if (range > 1)
				random = createRandomGenerator();

			// Add frames at the specified density until the number of particles has been reached
			int id = 0;
			int t = 0;
			while (id < settings.getParticles())
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

			areaInUm = settings.getSize() * settings.getPixelPitch() * settings.getSize() * settings.getPixelPitch() /
					1e6;

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
				settings.setSeconds((int) Math
						.ceil(settings.getParticles() * (settings.getExposureTime() + settings.getTOn()) / 1000));
				totalSteps = 0;

				final double simulationStepsPerFrame = (settings.getStepsPerSecond() * settings.getExposureTime()) /
						1000.0;
				imageModel = new FixedLifetimeImageModel(settings.getStepsPerSecond() * settings.getTOn() / 1000.0,
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

				SpatialIllumination activationIllumination = createIllumination(settings.getPulseRatio(),
						settings.getPulseInterval());

				// Generate additional frames so that each frame has the set number of simulation steps
				totalSteps = (int) Math.ceil(settings.getSeconds() * settings.getStepsPerSecond());

				// Since we have an exponential decay of activations
				// ensure half of the particles have activated by 30% of the frames.
				double eAct = totalSteps * 0.3 * activationIllumination.getAveragePhotons();

				// Q. Does tOn/tOff change depending on the illumination strength?
				imageModel = new ActivationEnergyImageModel(eAct, activationIllumination,
						settings.getStepsPerSecond() * settings.getTOn() / 1000.0,
						settings.getStepsPerSecond() * settings.getTOffShort() / 1000.0,
						settings.getStepsPerSecond() * settings.getTOffLong() / 1000.0, settings.getNBlinksShort(),
						settings.getNBlinksLong());
				imageModel.setUseGeometricDistribution(settings.getNBlinksGeometricDistribution());

				// Only use the correlation if selected for the distribution
				if (PHOTON_DISTRIBUTION[PHOTON_CORRELATED].equals(settings.getPhotonDistribution()))
					correlation = settings.getCorrelation();
			}

			imageModel.setRandomGenerator(createRandomGenerator());
			imageModel.setPhotonBudgetPerFrame(true);
			imageModel.setDiffusion2D(settings.getDiffuse2D());
			imageModel.setRotation2D(settings.getRotate2D());

			IJ.showStatus("Creating molecules ...");
			SpatialDistribution distribution = createDistribution();
			List<CompoundMoleculeModel> compounds = createCompoundMolecules();
			if (compounds == null)
				return;
			List<CompoundMoleculeModel> molecules = imageModel.createMolecules(compounds, settings.getParticles(),
					distribution, settings.getRotateInitialOrientation());

			// Activate fluorophores
			IJ.showStatus("Creating fluorophores ...");

			// Note: molecules list will be converted to compounds containing fluorophores
			fluorophores = imageModel.createFluorophores(molecules, totalSteps);

			if (fluorophores.isEmpty())
			{
				IJ.error(TITLE, "No fluorophores created");
				return;
			}

			// Map the fluorophore ID to the compound for mixtures 
			if (compounds.size() > 1)
			{
				idToCompound = new TIntIntHashMap(fluorophores.size());
				for (FluorophoreSequenceModel l : fluorophores)
				{
					idToCompound.put(l.getId(), l.getLabel());
				}
			}

			IJ.showStatus("Creating localisations ...");

			// TODO - Output a molecule Id for each fluorophore if using compound molecules. This allows analysis
			// of the ratio of trimers, dimers, monomers, etc that could be detected.

			totalSteps = checkTotalSteps(totalSteps, fluorophores);
			if (totalSteps == 0)
				return;

			imageModel.setPhotonDistribution(createPhotonDistribution());
			try
			{
				imageModel.setConfinementDistribution(createConfinementDistribution());
			}
			catch (ConfigurationException e)
			{
				// We asked the user if it was OK to continue and they said no
				return;
			}
			// This should be optimised
			imageModel.setConfinementAttempts(10);

			localisations = imageModel.createImage(molecules, settings.getFixedFraction(), totalSteps,
					(double) settings.getPhotonsPerSecond() / settings.getStepsPerSecond(), correlation,
					settings.getRotateDuringSimulation());

			// Re-adjust the fluorophores to the correct time
			if (settings.getStepsPerSecond() != 1)
			{
				final double scale = 1.0 / settings.getStepsPerSecond();
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
		SettingsManager.writeSettings(settings.build());

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
		Utils.log(TITLE + " Benchmark");

		double[] xyz = dist.next().clone();
		double offset = settings.getSize() * 0.5;
		for (int i = 0; i < 2; i++)
			xyz[i] += offset;
		double sd = getPsfSD() * settings.getPixelPitch();

		Utils.log("X = %s nm : %s px", Utils.rounded(xyz[0] * settings.getPixelPitch()), Utils.rounded(xyz[0], 6));
		Utils.log("Y = %s nm : %s px", Utils.rounded(xyz[1] * settings.getPixelPitch()), Utils.rounded(xyz[1], 6));
		Utils.log("Width (s) = %s nm : %s px", Utils.rounded(sd), Utils.rounded(sd / settings.getPixelPitch()));
		final double sa = PSFCalculator.squarePixelAdjustment(sd, settings.getPixelPitch());
		Utils.log("Adjusted Width (sa) = %s nm : %s px", Utils.rounded(sa),
				Utils.rounded(sa / settings.getPixelPitch()));
		Utils.log("Signal (N) = %s - %s photons", Utils.rounded(settings.getPhotonsPerSecond()),
				Utils.rounded(settings.getPhotonsPerSecondMaximum()));

		boolean emCCD;
		double totalGain;
		double noise = getNoiseEstimate();
		double readNoise;

		if (CalibrationProtosHelper.isCCDCameraType(settings.getCameraType()))
		{
			CreateDataSettingsHelper helper = new CreateDataSettingsHelper(settings);
			emCCD = (settings.getCameraType() == CameraType.EMCCD) ? settings.getEmGain() > 1 : false;
			totalGain = helper.getTotalGainSafe();
			// Store read noise in ADUs
			readNoise = settings.getReadNoise() * ((settings.getCameraGain() > 0) ? settings.getCameraGain() : 1);
		}
		else if (settings.getCameraType() == CameraType.SCMOS)
		{
			// Assume sCMOS amplification is like a CCD for the precision computation.
			emCCD = false;
			// Not required for sCMOS
			totalGain = 0;
			readNoise = 0;
		}
		else
		{
			throw new IllegalStateException("Unknown camera type: " + settings.getCameraType());
		}

		// The precision calculation is dependent on the model. The classic Mortensen formula
		// is for a Gaussian Mask Estimator. Use other equation for MLE. The formula provided 
		// for WLSE requires an offset to the background used to stabilise the fitting. This is
		// not implemented (i.e. we used an offset of zero) and in this case the WLSE precision 
		// is the same as MLE with the caveat of numerical instability.

		double lowerP = Gaussian2DPeakResultHelper.getPrecision(settings.getPixelPitch(), sd,
				settings.getPhotonsPerSecondMaximum(), noise, emCCD);
		double upperP = Gaussian2DPeakResultHelper.getPrecision(settings.getPixelPitch(), sd,
				settings.getPhotonsPerSecond(), noise, emCCD);
		double lowerMLP = Gaussian2DPeakResultHelper.getMLPrecision(settings.getPixelPitch(), sd,
				settings.getPhotonsPerSecondMaximum(), noise, emCCD);
		double upperMLP = Gaussian2DPeakResultHelper.getMLPrecision(settings.getPixelPitch(), sd,
				settings.getPhotonsPerSecond(), noise, emCCD);
		double lowerN = getPrecisionN(settings.getPixelPitch(), sd, settings.getPhotonsPerSecond(), Maths.pow2(noise),
				emCCD);
		double upperN = getPrecisionN(settings.getPixelPitch(), sd, settings.getPhotonsPerSecondMaximum(),
				Maths.pow2(noise), emCCD);

		if (settings.getCameraType() == CameraType.SCMOS)
			Utils.log("sCMOS camera background estimate uses an average read noise");
		Utils.log("Effective background noise = %s photons [includes read variance converted to photons]",
				Utils.rounded(noise));
		Utils.log("Localisation precision (LSE): %s - %s nm : %s - %s px", Utils.rounded(lowerP), Utils.rounded(upperP),
				Utils.rounded(lowerP / settings.getPixelPitch()), Utils.rounded(upperP / settings.getPixelPitch()));
		Utils.log("Localisation precision (MLE): %s - %s nm : %s - %s px", Utils.rounded(lowerMLP),
				Utils.rounded(upperMLP), Utils.rounded(lowerMLP / settings.getPixelPitch()),
				Utils.rounded(upperMLP / settings.getPixelPitch()));
		Utils.log("Signal precision: %s - %s photons", Utils.rounded(lowerN), Utils.rounded(upperN));

		// Store the benchmark settings when not using variable photons
		if (settings.getPhotonsPerSecond() == settings.getPhotonsPerSecondMaximum())
		{
			final double qe = getQuantumEfficiency();
			benchmarkParameters = new BenchmarkParameters(settings.getParticles(), sd, settings.getPixelPitch(),
					settings.getPhotonsPerSecond(), xyz[0], xyz[1], xyz[2], settings.getBias(), totalGain, qe,
					readNoise, settings.getCameraType(), settings.getCameraModelName(), cameraModel.getBounds(),
					settings.getBackground(), noise, lowerN, lowerP, lowerMLP);
		}
		else
		{
			Utils.log(
					"Warning: Benchmark settings are only stored in memory when the number of photons is fixed. Min %s != Max %s",
					Utils.rounded(settings.getPhotonsPerSecond()),
					Utils.rounded(settings.getPhotonsPerSecondMaximum()));
		}
	}

	private double getNoiseEstimate()
	{
		if (CalibrationProtosHelper.isCCDCameraType(settings.getCameraType()))
		{
			// Background is in photons
			double backgroundVariance = settings.getBackground();

			// Read noise is in electrons. Convert to Photons
			double readNoise = settings.getReadNoise() / getQuantumEfficiency();

			// In an EM-CCD camera the read noise (in Counts) is swamped by amplification of the signal. 
			// We get the same result by dividing the read noise (in photons) by the EM-gain.
			if (settings.getCameraType() == CameraType.EMCCD && settings.getEmGain() > 1)
			{
				// Add EM-CCD noise factor. The Mortensen formula also includes this factor 
				// so this is "double-counting" the EM-CCD.  
				backgroundVariance *= 2;
				readNoise /= settings.getEmGain();
			}

			// Get the expected value at each pixel in photons. Assuming a Poisson distribution this 
			// is equal to the total variance at the pixel.
			return Math.sqrt(backgroundVariance + Maths.pow2(readNoise));
		}
		else if (settings.getCameraType() == CameraType.SCMOS)
		{
			// Assume sCMOS amplification is like a CCD. We need an average read noise to get an 
			// approximation of the background noise for the precision computation.

			// We get the total background in photons  
			double backgroundVariance = settings.getBackground();

			// Create the camera noise model
			createPerPixelCameraModelData(cameraModel);

			// Get the average read noise. Convert from electrons to photons
			double readNoise = (Maths.sum(this.readNoise) / this.readNoise.length) / getQuantumEfficiency();

			return Math.sqrt(backgroundVariance + Maths.pow2(readNoise));
		}
		throw new IllegalStateException("Unknown camera type: " + settings.getCameraType());
	}

	/**
	 * Store the simulation settings
	 * 
	 * @param particles
	 */
	private void saveSimulationParameters(int particles, boolean fullSimulation, double signalPerFrame)
	{
		double totalGain;
		double noise = getNoiseEstimate();
		double readNoise;

		if (CalibrationProtosHelper.isCCDCameraType(settings.getCameraType()))
		{
			CreateDataSettingsHelper helper = new CreateDataSettingsHelper(settings);
			totalGain = helper.getTotalGainSafe();
			// Store read noise in ADUs
			readNoise = settings.getReadNoise() * ((settings.getCameraGain() > 0) ? settings.getCameraGain() : 1);
		}
		else if (settings.getCameraType() == CameraType.SCMOS)
		{
			// Not required for sCMOS
			totalGain = 0;
			readNoise = 0;
		}
		else
		{
			throw new IllegalStateException("Unknown camera type: " + settings.getCameraType());
		}

		double sd = getPsfSD() * settings.getPixelPitch();

		final double qe = getQuantumEfficiency();

		simulationParameters = new SimulationParameters(particles, fullSimulation, sd, settings.getPixelPitch(),
				settings.getPhotonsPerSecond(), settings.getPhotonsPerSecondMaximum(), signalPerFrame,
				settings.getDepth(), settings.getFixedDepth(), settings.getBias(), totalGain, qe, readNoise,
				settings.getCameraType(), settings.getCameraModelName(), cameraModel.getBounds(),
				settings.getBackground(), noise);
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
			final double simulationStepsPerFrame = (settings.getStepsPerSecond() * settings.getExposureTime()) / 1000.0;
			int newFrames = 1 + (int) (max / simulationStepsPerFrame);

			if (totalSteps != 0)
			{
				int totalFrames = (int) Math.ceil(settings.getSeconds() * 1000 / settings.getExposureTime());
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
		if (settings.getDistribution().equals(DISTRIBUTION[MASK]))
		{
			ImagePlus imp = WindowManager.getImage(settings.getDistributionMask());
			if (imp != null)
			{
				return createMaskDistribution(imp, settings.getDistributionMaskSliceDepth(), true);
			}
		}
		else if (settings.getDistribution().equals(DISTRIBUTION[GRID]))
		{
			return new GridDistribution(settings.getSize(), settings.getDepth() / settings.getPixelPitch(),
					settings.getCellSize(), settings.getProbabilityBinary(),
					settings.getMinBinaryDistance() / settings.getPixelPitch(),
					settings.getMaxBinaryDistance() / settings.getPixelPitch());
		}

		return createUniformDistributionWithPSFWidthBorder();
	}

	private SpatialDistribution createMaskDistribution(ImagePlus imp, double sliceDepth, boolean updateArea)
	{
		// Calculate the scale of the mask
		final int w = imp.getWidth();
		final int h = imp.getHeight();
		final double scaleX = (double) settings.getSize() / w;
		final double scaleY = (double) settings.getSize() / h;

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
				sliceDepth = settings.getDepth() / masks.size();

			return new MaskDistribution3D(masks, w, h, sliceDepth / settings.getPixelPitch(), scaleX, scaleY,
					createRandomGenerator());
		}
		else
		{
			int[] mask = extractMask(imp.getProcessor());
			if (updateArea)
				updateArea(mask, w, h);
			return new MaskDistribution(mask, w, h, settings.getDepth() / settings.getPixelPitch(), scaleX, scaleY,
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
		final double scaleX = (double) settings.getSize() / w;
		final double scaleY = (double) settings.getSize() / h;
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
		areaInUm = scale * settings.getSize() * settings.getPixelPitch() * settings.getSize() *
				settings.getPixelPitch() / 1e6;
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
		border = FastMath.min(border, settings.getSize() / 4);
		return createUniformDistribution(border);
	}

	private SpatialDistribution createFixedDistribution()
	{
		SpatialDistribution dist;
		dist = new SpatialDistribution()
		{
			private double[] xyz = new double[] { settings.getXPosition() / settings.getPixelPitch(),
					settings.getYPosition() / settings.getPixelPitch(),
					settings.getZPosition() / settings.getPixelPitch() };

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
			if (psfModelType == PSF_MODEL_IMAGE)
			{
				hwhm = getImageHWHM();
			}
			else if (psfModelType == PSF_MODEL_ASTIGMATISM)
			{
				hwhm = getAstigmatismHWHM();
			}
			else
			{
				final double sd = (settings.getEnterWidth()) ? settings.getPsfSd()
						: PSFCalculator.calculateStdDev(settings.getWavelength(), settings.getNumericalAperture());

				hwhm = Gaussian2DFunction.SD_TO_HWHM_FACTOR * sd / settings.getPixelPitch();
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
		ImagePlus imp = WindowManager.getImage(settings.getPsfImageName());
		if (imp == null)
		{
			IJ.error(TITLE, "Unable to create the PSF model from image: " + settings.getPsfImageName());
			return -1;
		}
		ImagePSF psfSettings = ImagePSFHelper.fromString(imp.getProperty("Info").toString());
		if (psfSettings == null)
		{
			IJ.error(TITLE, "Unknown PSF settings for image: " + imp.getTitle());
			return -1;
		}
		if (psfSettings.getFwhm() <= 0)
		{
			IJ.error(TITLE, "Unknown PSF FWHM setting for image: " + imp.getTitle());
			return -1;
		}
		if (psfSettings.getPixelSize() <= 0)
		{
			IJ.error(TITLE, "Unknown PSF pixel size setting for image: " + imp.getTitle());
			return -1;
		}

		// The width of the PSF is specified in pixels of the PSF image. Convert to the pixels of the 
		// output image
		return 0.5 * psfSettings.getFwhm() * psfSettings.getPixelSize() / settings.getPixelPitch();
	}

	private double getAstigmatismHWHM()
	{
		AstigmatismModel model = AstigmatismModelManager.getModel(settings.getAstigmatismModel());
		if (model == null)
		{
			IJ.error(TITLE, "Unknown PSF model: " + settings.getAstigmatismModel());
			return -1;
		}
		try
		{
			// Get the width at z=0 in pixels
			model = AstigmatismModelManager.convert(model, model.getZDistanceUnit(), DistanceUnit.PIXEL);
			AstigmatismZModel zModel = AstigmatismModelManager.create(model);
			double sx = zModel.getSx(0);
			double sy = zModel.getSy(0);
			return Gaussian2DPeakResultHelper.getStandardDeviation(sx, sy) * Gaussian2DFunction.SD_TO_HWHM_FACTOR
			// Scale appropriately
					* model.getNmPerPixel() / settings.getPixelPitch();
		}
		catch (ConversionException e)
		{
			IJ.error(TITLE, "Unknown PSF FWHM setting for model: " + settings.getAstigmatismModel());
			return -1;
		}
	}

	/**
	 * Create distribution within an XY border
	 * 
	 * @param border
	 * @return
	 */
	private UniformDistribution createUniformDistribution(double border)
	{
		double depth = (settings.getFixedDepth()) ? settings.getDepth() / settings.getPixelPitch()
				: settings.getDepth() / (2 * settings.getPixelPitch());

		// Ensure the focal plane is in the middle of the zDepth
		double[] max = new double[] { settings.getSize() / 2 - border, settings.getSize() / 2 - border, depth };
		double[] min = new double[3];
		for (int i = 0; i < 3; i++)
			min[i] = -max[i];
		if (settings.getFixedDepth())
			min[2] = max[2];

		// Try using different distributions:
		final RandomGenerator rand1 = createRandomGenerator();

		if (settings.getDistribution().equals(DISTRIBUTION[UNIFORM_HALTON]))
		{
			return new UniformDistribution(min, max, rand1.nextInt());
		}

		if (settings.getDistribution().equals(DISTRIBUTION[UNIFORM_SOBOL]))
		{
			SobolSequenceGenerator rvg = new SobolSequenceGenerator(3);
			rvg.skipTo(rand1.nextInt());
			return new UniformDistribution(min, max, rvg);
		}

		// Create a distribution using random generators for each dimension 
		UniformDistribution distribution = new UniformDistribution(min, max, this);
		return distribution;
	}

	private SpatialDistribution createConfinementDistribution() throws ConfigurationException
	{
		if (settings.getDiffusionRate() <= 0 || settings.getFixedFraction() >= 1)
			return null;

		// Log a warning if the confinement conflicts with the distribution

		if (settings.getConfinement().equals(CONFINEMENT[CONFINEMENT_MASK]))
		{
			// The mask should be the same
			if (!settings.getDistribution().equals(DISTRIBUTION[MASK]))
				checkConfiguration("Simulation uses a mask confinement but no mask distribution");
			else if (!settings.getConfinementMask().equals(settings.getDistributionMask()))
				checkConfiguration(
						"Simulation uses a mask confinement with a different image to the mask distribution");
			else if (settings.getConfinementMaskSliceDepth() != settings.getDistributionMaskSliceDepth())
				checkConfiguration(
						"Simulation uses a mask confinement with a different depth to the mask distribution");

			ImagePlus imp = WindowManager.getImage(settings.getConfinementMask());
			if (imp != null)
			{
				return createMaskDistribution(imp, settings.getConfinementMaskSliceDepth(), false);
			}
		}
		else if (settings.getConfinement().equals(CONFINEMENT[CONFINEMENT_SPHERE]))
		{
			// This may be an error if the distribution is a mask
			if (settings.getDistribution().equals(DISTRIBUTION[MASK]))
				checkConfiguration("Simulation uses a mask confinement but a " + CONFINEMENT[CONFINEMENT_WITHIN_IMAGE] +
						" confinement");

			return new SphericalDistribution(settings.getConfinementRadius() / settings.getPixelPitch());
		}
		else if (settings.getConfinement().equals(CONFINEMENT[CONFINEMENT_WITHIN_IMAGE]))
		{
			// This may be an error if the distribution is a mask
			if (settings.getDistribution().equals(DISTRIBUTION[MASK]))
				checkConfiguration("Simulation uses a mask confinement but a " + CONFINEMENT[CONFINEMENT_WITHIN_IMAGE] +
						" confinement");

			//return createUniformDistribution(0);
			return createUniformDistributionWithPSFWidthBorder();
		}

		if (settings.getDistribution().equals(DISTRIBUTION[MASK]))
			checkConfiguration("Simulation uses a mask confinement but no confinement");
		return null;
	}

	private void checkConfiguration(String message) throws ConfigurationException
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addMessage(TextUtils.wrap("Warning: " + message, 80));
		gd.setOKLabel("Continue");
		gd.showDialog();
		if (gd.wasCanceled())
			throw new ConfigurationException(message);
	}

	private SpatialIllumination createIllumination(double intensity, int pulseInterval)
	{
		if (settings.getIllumination().equals(ILLUMINATION[RADIAL]))
		{
			if (pulseInterval > 1)
				return new RadialFalloffIllumination(1, settings.getSize() / 2, intensity, pulseInterval);
			return new RadialFalloffIllumination(intensity, settings.getSize() / 2);
		}
		else
		{
			uniformBackground = true;
			if (pulseInterval > 1)
				return new UniformIllumination(1, intensity, pulseInterval);
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
		if (PHOTON_DISTRIBUTION[PHOTON_CUSTOM].equals(settings.getPhotonDistribution()))
		{
			// Get the distribution file
			String filename = Utils.getFilename("Photon_distribution", settings.getPhotonDistributionFile());
			if (filename != null)
			{
				settings.setPhotonDistributionFile(filename);
				try
				{
					InputStream is = new FileInputStream(new File(settings.getPhotonDistributionFile()));
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
						double scale = (double) settings.getPhotonsPerSecond() / stats.getMean();
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
					settings.getPhotonDistributionFile());
		}
		else if (PHOTON_DISTRIBUTION[PHOTON_UNIFORM].equals(settings.getPhotonDistribution()))
		{
			if (settings.getPhotonsPerSecond() < settings.getPhotonsPerSecondMaximum())
			{
				UniformRealDistribution dist = new UniformRealDistribution(createRandomGenerator(),
						settings.getPhotonsPerSecond(), settings.getPhotonsPerSecondMaximum());
				return dist;
			}
		}
		else if (PHOTON_DISTRIBUTION[PHOTON_GAMMA].equals(settings.getPhotonDistribution()))
		{
			final double scaleParameter = settings.getPhotonsPerSecond() / settings.getPhotonShape();
			GammaDistribution dist = new GammaDistribution(createRandomGenerator(), settings.getPhotonShape(),
					scaleParameter, ExponentialDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);
			return dist;
		}
		else if (PHOTON_DISTRIBUTION[PHOTON_CORRELATED].equals(settings.getPhotonDistribution()))
		{
			// No distribution required
			return null;
		}

		settings.setPhotonDistribution(PHOTON_DISTRIBUTION[PHOTON_FIXED]);
		return null;
	}

	private List<LocalisationModelSet> combineSimulationSteps(List<LocalisationModel> localisations)
	{
		// Allow fractional integration steps
		final double simulationStepsPerFrame = (settings.getStepsPerSecond() * settings.getExposureTime()) / 1000.0;

		List<LocalisationModelSet> newLocalisations = new ArrayList<LocalisationModelSet>(
				(int) (localisations.size() / simulationStepsPerFrame));

		//System.out.printf("combineSimulationSteps @ %f\n", simulationStepsPerFrame);

		//final double gain = new CreateDataSettingsHelper(settings).getTotalGainSafe();
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
				{
					movingMolecules.add(id);
				}

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

						// Create a data array and store the current intensity. 
						// This is used later to filter based on SNR
						sets[i].setData(new double[] { 0, 0, 0, 0, sets[i].getIntensity() //* gain
						});

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
		if (localisations.isEmpty())
			return new int[0];
		TIntArrayList ids = new TIntArrayList(settings.getParticles());
		// Assume the localisations are sorted by id
		int id = localisations.get(0).getId();
		ids.add(id);
		for (LocalisationModel l : localisations)
		{
			if (id != l.getId())
			{
				id = l.getId();
				ids.add(id);
			}
		}
		return ids.toArray();
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
		CalibrationWriter c = new CalibrationWriter();
		c.setDistanceUnit(DistanceUnit.PIXEL);
		c.setIntensityUnit(IntensityUnit.PHOTON);
		c.setNmPerPixel(settings.getPixelPitch());
		c.setExposureTime(settings.getExposureTime());

		c.setCameraType(settings.getCameraType());
		if (settings.getCameraType() == CameraType.SCMOS)
		{
			c.setCameraModelName(settings.getCameraModelName());
		}
		else
		{
			CreateDataSettingsHelper helper = new CreateDataSettingsHelper(settings);
			c.setCountPerPhoton(helper.getTotalGainSafe());
			c.setBias(settings.getBias());
			c.setReadNoise(settings.getReadNoise() * ((settings.getCameraGain() > 0) ? settings.getCameraGain() : 1));
			c.setQuantumEfficiency(helper.getQuantumEfficiency());
		}

		results.setCalibration(c.getCalibration());
		results.setSortAfterEnd(true);
		results.begin();

		maxT = localisationSets.get(localisationSets.size() - 1).getTime();

		// Display image
		ImageStack stack = new ImageStack(settings.getSize(), settings.getSize(), maxT);

		final double psfSD = getPsfSD();
		if (psfSD <= 0)
			return null;
		PSFModel psfModel = createPSFModel(localisationSets);
		if (psfModel == null)
			return null;

		// Create the camera noise model
		createPerPixelCameraModelData(cameraModel);

		IJ.showStatus("Drawing image ...");

		// Multi-thread for speed
		// Note that the default Executors.newCachedThreadPool() will continue to make threads if
		// new tasks are added. We need to limit the tasks that can be added using a fixed size
		// blocking queue.
		// http://stackoverflow.com/questions/1800317/impossible-to-make-a-cached-thread-pool-with-a-size-limit
		// ExecutorService threadPool = Executors.newCachedThreadPool();
		int threadCount = Prefs.getThreads();
		PeakResults syncResults = SynchronizedPeakResults.create(results, threadCount);
		ExecutorService threadPool = Executors.newFixedThreadPool(threadCount);
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
						new ImageGenerator(localisationSets, newLocalisations, i, lastT, createPSFModel(psfModel),
								syncResults, stack, poissonNoise, new RandomDataGenerator(createRandomGenerator()))));
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
		float[] limits = null;
		for (int t = 1; t <= maxT; t++)
		{
			if (Utils.isInterrupted())
				break;
			Object pixels = stack.getPixels(t);
			if (pixels == null)
			{
				futures.add(threadPool.submit(new ImageGenerator(localisationSets, newLocalisations, maxT, t, null,
						syncResults, stack, poissonNoise, new RandomDataGenerator(createRandomGenerator()))));
			}
			else if (limits == null)
				limits = Maths.limits((float[]) pixels);
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
		psfModel = null;
		threadPool = null;
		futures.clear();
		futures = null;

		if (photonsRemoved.get() > 0)
			Utils.log("Removed %d localisations with less than %.1f rendered photons", photonsRemoved.get(),
					settings.getMinPhotons());
		if (t1Removed.get() > 0)
			Utils.log("Removed %d localisations with no neighbours @ SNR %.2f", t1Removed.get(),
					settings.getMinSnrT1());
		if (tNRemoved.get() > 0)
			Utils.log("Removed %d localisations with valid neighbours @ SNR %.2f", tNRemoved.get(),
					settings.getMinSnrTN());
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

		if (!settings.getRawImage())
		{
			// Get the global limits and ensure all values can be represented
			Object[] imageArray = stack.getImageArray();
			limits = Maths.limits((float[]) imageArray[0]);
			for (int j = 1; j < imageArray.length; j++)
				limits = Maths.limits(limits, (float[]) imageArray[j]);
			//float limits0 = limits[0];
			float limits0 = 0; // Leave bias in place
			// Check if the image will fit in a 16-bit range
			if ((limits[1] - limits0) < 65535)
			{
				// Convert to 16-bit
				newStack = new ImageStack(stack.getWidth(), stack.getHeight(), stack.getSize());
				// Account for rounding
				final float min = (float) (limits0 - 0.5);
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
				for (int k = 2; k-- > 0;)
					limits[k] = (float) Math.floor(limits[k] - min);
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
				for (int k = 2; k-- > 0;)
					limits[k] = Math.round(limits[k]);
			}
		}

		// Show image
		ImagePlus imp = Utils.display(CREATE_DATA_IMAGE_TITLE, newStack);

		ij.measure.Calibration cal = new ij.measure.Calibration();
		String unit = "nm";
		double unitPerPixel = settings.getPixelPitch();
		if (unitPerPixel > 100)
		{
			unit = "um";
			unitPerPixel /= 1000.0;
		}
		cal.setUnit(unit);
		cal.pixelHeight = cal.pixelWidth = unitPerPixel;
		Rectangle bounds = cameraModel.getBounds();
		if (bounds != null)
		{
			cal.xOrigin = -bounds.x;
			cal.yOrigin = -bounds.y;
		}
		imp.setCalibration(cal);

		imp.setDimensions(1, 1, newStack.getSize());
		imp.setDisplayRange(limits[0], limits[1]);
		//imp.resetDisplayRange();
		imp.updateAndDraw();

		saveImage(imp);

		IJImageSource imageSource = new IJImageSource(imp);
		// Shift simulation image source to correct location
		results.setSource(imageSource);
		results.setName(CREATE_DATA_IMAGE_TITLE + " (" + TITLE + ")");
		// Bounds are relative to the image source
		results.setBounds(new Rectangle(settings.getSize(), settings.getSize()));

		PSF psf;
		if (astigmatismModel != null)
		{
			psf = PSFProtosHelper.createPSF(astigmatismModel, DistanceUnit.PIXEL, DistanceUnit.PIXEL);
		}
		else
		{
			PSF.Builder psfBuilder;
			// Set the PSF as a Gaussian using the width at z=0. 
			// In future this could be improved for other PSFs.S
			psfBuilder = PSFProtosHelper.defaultOneAxisGaussian2DPSF.toBuilder();
			psfBuilder.getParametersBuilder(PSFHelper.INDEX_SX).setValue(psfSD);
			psf = psfBuilder.build();
		}
		results.setPSF(psf);
		MemoryPeakResults.addResults(results);

		setBenchmarkResults(imp, results);

		if (benchmarkMode && benchmarkParameters != null)
			benchmarkParameters.setPhotons(results);

		List<LocalisationModel> localisations = toLocalisations(localisationSets);

		savePulses(localisations, results, CREATE_DATA_IMAGE_TITLE);

		// Saved the fixed and moving localisations into different datasets
		saveFixedAndMoving(results, CREATE_DATA_IMAGE_TITLE);

		saveCompoundMolecules(results, CREATE_DATA_IMAGE_TITLE);

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
		ImagePlus imp = WindowManager.getImage(settings.getPsfImageName());
		if (imp == null)
		{
			IJ.error(TITLE, "Unable to create the PSF model from image: " + settings.getPsfImageName());
			return null;
		}
		try
		{
			ImagePSF psfSettings = ImagePSFHelper.fromString(imp.getProperty("Info").toString());
			if (psfSettings == null)
				throw new RuntimeException("Unknown PSF settings for image: " + imp.getTitle());

			// Check all the settings have values
			if (psfSettings.getPixelSize() <= 0)
				throw new RuntimeException("Missing nmPerPixel calibration settings for image: " + imp.getTitle());
			if (psfSettings.getPixelDepth() <= 0)
				throw new RuntimeException("Missing nmPerSlice calibration settings for image: " + imp.getTitle());
			if (psfSettings.getCentreImage() <= 0)
				throw new RuntimeException("Missing zCentre calibration settings for image: " + imp.getTitle());
			if (psfSettings.getFwhm() <= 0)
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
			int zCentre = psfSettings.getCentreImage() - 1;

			// Calculate the start/end slices to cover the depth of field
			// This logic must match the ImagePSFModel.
			final double unitsPerSlice = psfSettings.getPixelDepth() / settings.getPixelPitch();
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
			float[][] image = extractImageStack(imp, lower, upper);
			final ImagePSFModel model = new ImagePSFModel(image, zCentre - lower,
					psfSettings.getPixelSize() / settings.getPixelPitch(), unitsPerSlice, psfSettings.getFwhm(),
					noiseFraction);

			// Add the calibrated centres. The map will not be null
			Map<Integer, Offset> map = psfSettings.getOffsetsMap();
			if (!map.isEmpty())
			{
				final int sliceOffset = lower + 1;
				for (Entry<Integer, Offset> entry : map.entrySet())
				{
					model.setRelativeCentre(entry.getKey() - sliceOffset, entry.getValue().getCx(),
							entry.getValue().getCy());
				}
			}
			else
			{
				// Use the CoM if present
				double cx = psfSettings.getXCentre();
				double cy = psfSettings.getYCentre();
				if (cx != 0 || cy != 0)
				{
					for (int slice = 0; slice < image.length; slice++)
						model.setCentre(slice, cx, cy);
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

	/**
	 * Extract the image stack using a range of stack indices. The index is 0-based.
	 *
	 * @param imp
	 *            the imp
	 * @param start
	 *            the start index
	 * @param end
	 *            the end index
	 * @return the float image stack
	 */
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

	private PSFModel createPSFModel(List<LocalisationModelSet> localisationSets) throws IllegalArgumentException
	{
		if (psfModelType == PSF_MODEL_IMAGE)
		{
			return createImagePSF(localisationSets);
		}

		if (psfModelType == PSF_MODEL_ASTIGMATISM)
		{
			astigmatismModel = AstigmatismModelManager.getModel(settings.getAstigmatismModel());
			if (astigmatismModel == null)
				throw new IllegalArgumentException("Failed to load model: " + settings.getAstigmatismModel());
			// Convert for simulation
			try
			{
				if (DoubleEquality.relativeError(astigmatismModel.getNmPerPixel(), settings.getPixelPitch()) > 1e-6)
				{
					String message = String.format(
							"Astigmatism model '%s' calibration (%s nm) does not match pixel pitch (%s nm)",
							settings.getAstigmatismModel(), Utils.rounded(astigmatismModel.getNmPerPixel()),
							Utils.rounded(settings.getPixelPitch()));
					// Optionally convert
					GenericDialog gd = new GenericDialog(TITLE);
					gd.addMessage(TextUtils.wrap(message + ". Created data is not suitable for fitting.", 80));
					gd.addMessage(TextUtils.wrap("Click OK to continue anyway (i.e. draw the spot using the " +
							"correct nm width on the different sized pixels).", 80));
					gd.showDialog();
					if (gd.wasCanceled())
						throw new IllegalArgumentException(message);
					// Convert to nm
					astigmatismModel = AstigmatismModelManager.convert(astigmatismModel, DistanceUnit.NM,
							DistanceUnit.NM);
					// Reset pixel pitch. This will draw the spot using the correct size on the different size pixels.
					astigmatismModel = astigmatismModel.toBuilder().setNmPerPixel(settings.getPixelPitch()).build();
				}

				// Convert for simulation in pixels
				astigmatismModel = AstigmatismModelManager.convert(astigmatismModel, DistanceUnit.PIXEL,
						DistanceUnit.PIXEL);
				return new GaussianPSFModel(AstigmatismModelManager.create(astigmatismModel));
			}
			catch (ConversionException e)
			{
				// Wrap so this can be caught as the same type
				throw new IllegalArgumentException(e);
			}
		}

		if (psfModelType == PSF_MODEL_GAUSSIAN)
		{
			double sd = getPsfSD();
			double d = settings.getDepthOfFocus() / settings.getPixelPitch();
			double gamma = 0;
			HoltzerAstigmatismZModel zModel = HoltzerAstigmatismZModel.create(sd, sd, gamma, d, 0, 0, 0, 0);
			return new GaussianPSFModel(createRandomGenerator(), zModel);
		}

		// Default to Airy pattern
		double width = getPsfSD() / PSFCalculator.AIRY_TO_GAUSSIAN;
		AiryPSFModel m = new AiryPSFModel(createRandomGenerator(), width, width, 450.0 / settings.getPixelPitch());
		m.setRing(2);
		return m;
	}

	private PSFModel createPSFModel(PSFModel psfModel) throws IllegalArgumentException
	{
		PSFModel copy = psfModel.copy();
		copy.setRandomGenerator(createRandomGenerator());
		return copy;

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

	private CameraModel cameraModel;

	/**
	 * The read noise.
	 * This is stored in electrons even though camera read noise is measured in ADUs. This allows it
	 * to be used to compute the combined background noise (with the background photon shot noise)
	 * for each localisation.
	 */
	private float[] readNoise;

	private void createPerPixelCameraModelData(CameraModel cameraModel)
	{
		if (readNoise != null)
			return;

		// Note: Store read noise in electrons for SNR computation. The gain is later applied back.

		if (cameraModel.isPerPixelModel())
		{
			Rectangle bounds = cameraModel.getBounds();
			readNoise = cameraModel.getVariance(bounds);
			float[] gain = cameraModel.getGain(bounds);
			for (int i = 0; i < readNoise.length; i++)
				readNoise[i] = (float) (Math.sqrt(readNoise[i]) / gain[i]);
		}
		else
		{
			// Use a dummy bounds to find out the fixed variance and gain 
			Rectangle bounds = new Rectangle(1, 1);
			float variance = cameraModel.getVariance(bounds)[0];
			float gain = cameraModel.getGain(bounds)[0];

			// Avoid sqrt on all the same value
			this.readNoise = new float[settings.getSize() * settings.getSize()];
			Arrays.fill(this.readNoise, (float) (Math.sqrt(variance) / gain));
		}

		// Remove if it will have no effect
		readNoise = nullIf(readNoise, 0f);
	}

	private float[] nullIf(float[] data, float f)
	{
		for (int i = 0; i < data.length; i++)
			if (data[i] != f)
				return data;
		return null;
	}

	private CameraModel createCCDCameraModel()
	{
		float bias = settings.getBias();
		float gain = 1f;
		float readNoise = 0;

		if (settings.getCameraGain() != 0)
			gain = (float) settings.getCameraGain();

		if (settings.getReadNoise() > 0)
		{
			readNoise = (float) settings.getReadNoise();
			// Read noise is in electrons. Apply camera gain to get the noise in ADUs.
			if (settings.getCameraGain() != 0)
				readNoise *= settings.getCameraGain();
		}

		return new FixedPixelCameraModel(bias, gain, (float) Maths.pow2(readNoise));
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
		final PeakResults results;
		final ImageStack stack;
		final boolean poissonNoise;
		final RandomDataGenerator random;
		final double emGain, qe;

		public ImageGenerator(final List<LocalisationModelSet> localisationSets,
				List<LocalisationModelSet> newLocalisations, int startIndex, int t, PSFModel psfModel,
				PeakResults results, ImageStack stack, boolean poissonNoise, RandomDataGenerator random)
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
			// This could be >=1 but the rest of the code ignores EM-gain if it is <=1			
			emGain = (settings.getCameraType() == CameraType.EMCCD && settings.getEmGain() > 1) ? settings.getEmGain()
					: 0;
			qe = getQuantumEfficiency();
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

			showProgress();

			final boolean checkSNR = minSNRt1 > 0 || minSNRtN > 0;

			// Adjust XY dimensions since they are centred on zero
			final double xoffset = settings.getSize() * 0.5;

			float[] image = createBackground(random);
			float[] imageCache = image.clone();

			// Create read noise now so that we can calculate the true background noise
			float[] imageReadNoise = new float[image.length];
			if (readNoise != null)
			{
				RandomGenerator r = random.getRandomGenerator();
				for (int i = 0; i < imageReadNoise.length; i++)
					imageReadNoise[i] = (float) (readNoise[i] * r.nextGaussian());
			}

			// Extract the localisations and draw if we have a PSF model
			int fromIndex = findIndexByTime(localisations, startIndex, t);
			if (fromIndex > -1 && psfModel != null)
			{
				int toIndex = findLastIndexByTime(localisations, fromIndex, t);
				List<LocalisationModelSet> subset = localisations.subList(fromIndex, toIndex + 1);
				float[] data = new float[settings.getSize() * settings.getSize()];
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
								final int samples = (int) random.nextPoisson(intensity);
								intensity = samples;
								photonsRendered = psfModel.sample3D(data, settings.getSize(), settings.getSize(),
										samples, localisation.getX(), localisation.getY(), localisation.getZ());
								samplePositions = psfModel.getSamplePositions();
							}
							else
							{
								photonsRendered = psfModel.create3D(data, settings.getSize(), settings.getSize(),
										intensity, localisation.getX(), localisation.getY(), localisation.getZ(),
										false);
							}
						}
						//addDraw(photons);
						if (photonsRendered > 0)
						{
							totalPhotonsRendered += photonsRendered;
							spots[spotCount++] = new Spot(psfModel.getPSF(), psfModel.getX0min(), psfModel.getX0max(),
									psfModel.getX1min(), psfModel.getX1max(), samplePositions);
						}
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

					// Add to memory.
					// Use the actual intensity (not the total photons rendered)
					float intensity = (float) localisation.getIntensity();
					float x = (float) localisation.getX();
					float y = (float) localisation.getY();
					float z = (float) localisation.getZ();
					// 0.5 is the centre of the pixel so just round down.
					int origX = (int) x;
					int origY = (int) y;
					// Background and noise should be calculated using the
					// region covered by the PSF.
					double[] localStats = getStatistics(imageCache, imageReadNoise, origX, origY);
					float background = (float) (localStats[0]);

					// Note: The width estimate does not account for diffusion
					float sx, sy;
					if (psfModel instanceof GaussianPSFModel)
					{
						GaussianPSFModel m = (GaussianPSFModel) psfModel;
						sx = (float) m.getS0();
						sy = (float) m.getS1();
					}
					else if (psfModel instanceof AiryPSFModel)
					{
						AiryPSFModel m = (AiryPSFModel) psfModel;
						sx = (float) (m.getW0() * AiryPattern.FACTOR);
						sy = (float) (m.getW1() * AiryPattern.FACTOR);
					}
					else if (psfModel instanceof ImagePSFModel)
					{
						ImagePSFModel m = (ImagePSFModel) psfModel;
						sx = (float) (m.getHWHM0() / Gaussian2DFunction.SD_TO_HWHM_FACTOR);
						sy = (float) (m.getHWHM1() / Gaussian2DFunction.SD_TO_HWHM_FACTOR);
					}
					else
					{
						throw new IllegalStateException("Unknown PSF model");
					}

					// *-*-*-*-*
					// We want to compute the background noise at the localisation position in photons.
					// noise = shot noise + read noise
					// *-*-*-*-*
					// Note that the noise we are calculating is the noise that would be in the image with no
					// fluorophore present. This is the true background noise and it is the noise that is  
					// estimated by the Peak Fit plugin. This noise therefore IGNORES THE SHOT NOISE of the 
					// fluorophore SIGNAL.
					// *-*-*-*-*

					// (EM-)CCD camera:
					// 

					// The variance of the background image is currently in photons^2 
					double backgroundVariance = localStats[1];

					// Convert to electrons^2
					if (qe < 1)
						backgroundVariance *= Maths.pow2(qe);

					// Get the actual read noise applied to this part of the image. This is in electrons^2.
					double readVariance = localStats[3];

					// EM-gain noise factor: Adds sqrt(2) to the electrons input to the register.
					// All data 'read' through the EM-register must have this additional noise factor added.
					if (emGain != 0)
					{
						backgroundVariance *= 2; // Note we are using the variance (std.dev^2) so we use the factor 2

						// For an EM-CCD camera the EM amplification massively scales the electrophotons
						// and then read noise is added after. This is the same as descaling the read noise by 
						// the EM-gain so that they can be combined.
						// Note: Previously the noise computation was in ADUs:
						// backgroundVariance(photons) * totalGain^2 + readVariance(electrons) * cameraGain^2
						// backgroundVariance(electrons) * (emGain * cameraGain)^2 + readVariance(electrons) * cameraGain^2
						// To convert to electrons we divide by total gain:
						// (backgroundVariance * (emGain * cameraGain)^2 + readVariance * cameraGain^2) / (emGain * cameraGain)^2 
						// backgroundVariance + readVariance * cameraGain^2 / (emGain * cameraGain)^2
						// backgroundVariance + readVariance * cameraGain^2 / (emGain * cameraGain)^2
						// backgroundVariance + readVariance / emGain^2

						readVariance /= Maths.pow2(emGain);
					}

					// Overall noise can be calculated from the ‘root of sum of squares’ equation
					double totalNoise = Math.sqrt(backgroundVariance + readVariance);

					// Convert noise from electrons back to photons for convenience when computing SNR using the signal in photons.
					if (qe < 1)
						totalNoise /= qe;

					// Ensure the new data is added before the intensity is updated. This avoids 
					// syncronisation clashes in the getIntensity(...) function.
					// Use the total photons rendered for signal filtering.
					// [0] = background (photons)
					// [1] = total noise (photons)
					// [2] = Gaussian Sx
					// [3] = Gaussian Sy
					// [4] = total intensity (photons)
					localisationSet.setData(new double[] { localStats[0], totalNoise, sx, sy, totalPhotonsRendered });

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
					// TODO - This is no longer necessary as we now store the z in the results so update plugins to 
					// use the z-position from the results.
					float[] params = Gaussian2DPeakResultHelper.createTwoAxisParams(background, intensity, x, y, z, sx,
							sy);
					results.add(new IdPeakResult(t, origX, origY, 0, localisation.getZ(), (float) totalNoise, params,
							null, localisationSet.getId()));
				}

				for (int i = 0; i < image.length; i++)
					image[i] += data[i];
			}

			// Quantum efficiency: Model using binomial distribution
			if (qe < 1)
			{
				for (int i = 0; i < image.length; i++)
					image[i] = random.nextBinomial((int) image[i], qe);
			}

			// Apply EM gain and add Gaussian read noise after all the photons have been simulated
			if (emGain != 0)
			{
				final boolean tubbsModel = false;
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
						final double scale = emGain - 1 + 1 / image[i];
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
					double shape = 1;
					CustomGammaDistribution dist = new CustomGammaDistribution(random.getRandomGenerator(), shape,
							emGain, GammaDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);

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

			// Apply read noise (in electrons)
			if (readNoise != null)
			{
				for (int i = 0; i < image.length; i++)
					image[i] += imageReadNoise[i];
			}

			// Apply camera gain. Note that the noise component of the camera gain IS the 
			// read noise. Thus the read noise is expected to be different for each camera gain.
			// Also add the bias.
			cameraModel.applyGainAndBias(cameraModel.getBounds(), image);

			{
				// TODO Compute the Fisher information for each spot.
				// This is valid for a Poisson process (see Smith et al, 2010):
				// Iaa = sum(i) (dYi da) * (dYi da) / Yi
				// Q. What about EM-CCD? We actually require the difference in the 
				// log-likelihood function.
				// 1. Draw each spot perfectly on a new image using psfModel.create3D
				// 2. Apply perfect gain to get the correct scale 
				// 3. For each spot
				// - Combine the other spots + background photons to get the background function
				// - Combine with the target spot and apply gain to get Yi	
				// - Render the target spot using an offset +/- delta on X,Y,Signal (a) to get dYi da
				// - Compute Iaa

				// The simulated image is not used?

			}

			// Send to output
			stack.setPixels(image, t);
		}

		private void erase(float[] data, Spot spot)
		{
			if (spot.samplePositions != null)
			{
				psfModel.eraseSample(data, settings.getSize(), settings.getSize(), spot.samplePositions);
			}
			else
			{
				psfModel.erase(data, settings.getSize(), settings.getSize(), spot.psf, spot.x0min, spot.x0max,
						spot.x1min, spot.x1max);
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
			final int width = FastMath.min(settings.getSize(), 2 * n + 1);
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
				int indexTo = (y + x1min) * settings.getSize() + x0min;
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
		if (coord >= settings.getSize())
			return settings.getSize() - 1;
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
		double minSNR = settings.getMinSnrT1();
		AtomicInteger counter = t1Removed;

		if (localisationSet.hasNeighbour())
		{
			double nextIntensity = getIntensity(localisationSet.getNext());
			double previousIntensity = getIntensity(localisationSet.getPrevious());

			// Check if either neighbour is above the t1 threshold
			if (Math.max(nextIntensity, previousIntensity) / noise > settings.getMinSnrT1())
			{
				// If neighbours are bright then use a more lenient threshold
				minSNR = settings.getMinSnrTN();
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

	private float[] backgroundPixels = null;
	private boolean uniformBackground = false;

	private float[] createBackground(RandomDataGenerator random)
	{
		float[] pixels2 = null;

		if (settings.getBackground() > 0)
		{
			if (random == null)
				random = new RandomDataGenerator(createRandomGenerator());
			createBackgroundPixels();
			pixels2 = new float[backgroundPixels.length];

			// Add Poisson noise.
			if (uniformBackground)
			{
				double mean = backgroundPixels[0];

				// Simulate N photons hitting the image. The total photons (N) is 
				// the mean for each pixel multiplied by the number of pixels.
				// Note: The number of samples (N) must be Poisson distributed, i.e. 
				// the total amount of photons per frame is Poisson noise.
				final int samples = (int) random.nextPoisson(mean * backgroundPixels.length);

				final int upper = pixels2.length - 1;
				for (int i = 0; i < samples; i++)
				{
					pixels2[random.nextInt(0, upper)] += 1;
				}

				//// Alternative is to sample each pixel from a Poisson distribution. This is slow
				//PoissonDistribution dist = new PoissonDistribution(random.getRandomGenerator(), backgroundPixels[0],
				//		PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);
				//for (int i = 0; i < pixels2.length; i++)
				//{
				//	pixels2[i] = dist.sample();
				//}
			}
			else
			{
				for (int i = 0; i < pixels2.length; i++)
				{
					pixels2[i] = random.nextPoisson(backgroundPixels[i]);
				}
			}
		}
		else
		{
			pixels2 = new float[settings.getSize() * settings.getSize()];
		}

		return pixels2;
	}

	synchronized private void createBackgroundPixels()
	{
		// Cache illumination background
		if (backgroundPixels == null)
		{
			backgroundPixels = new float[settings.getSize() * settings.getSize()];

			ImagePlus imp = WindowManager.getImage(settings.getBackgroundImage());
			if (imp != null)
			{
				// Use an image for the background
				ImageProcessor ip = imp.getProcessor().duplicate().toFloat(0, null);
				ip.setInterpolationMethod(ImageProcessor.BILINEAR);
				ip = ip.resize(settings.getSize(), settings.getSize());
				float[] data = (float[]) ip.getPixels();
				final double max = Maths.maxDefault(0, data);
				if (max != 0)
				{
					final double scale = settings.getBackground() / max;
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
			SpatialIllumination illumination = createIllumination(settings.getBackground(), 0);
			double[] xyz = new double[3];
			for (int y = 0, i = 0; y < settings.getSize(); y++)
			{
				xyz[1] = y - settings.getSize() / 2;
				for (int x = 0, x2 = -settings.getSize() / 2; x < settings.getSize(); x++, x2++, i++)
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
			stats[i] = (settings.getShowHistograms() || alwaysRemoveOutliers[i]) ? new StoredDataStatistics()
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
		final double centreOffset = settings.getSize() * 0.5;
		// Used to convert the sampled times in frames into seconds
		final double framesPerSecond = 1000.0 / settings.getExposureTime();
		//final double gain = new CreateDataSettingsHelper(settings).getTotalGainSafe();
		for (LocalisationModel l : localisations)
		{
			double[] data = l.getData();
			if (data == null)
				throw new IllegalStateException("No localisation data. This should not happen!");
			final double noise = data[1];
			final double intensityInPhotons = data[4];
			// Q. What if the noise is zero, i.e. no background photon / read noise?
			// Just ignore it at current.
			final double snr = intensityInPhotons / noise;
			stats[SIGNAL].add(intensityInPhotons);
			stats[NOISE].add(noise);
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

			stats[X].add((l.getX() - centreOffset) * settings.getPixelPitch());
			stats[Y].add((l.getY() - centreOffset) * settings.getPixelPitch());
			stats[Z].add(l.getZ() * settings.getPixelPitch());
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
			// Convert depth-of-field to pixels
			final double depth = settings.getDepthOfField() / settings.getPixelPitch();

			try
			{
				// Get widths
				WidthResultProcedure wp = new WidthResultProcedure(results, DistanceUnit.PIXEL);
				wp.getW();
				stats[WIDTH].add(wp.wx);
			}
			catch (DataException e)
			{
				Utils.log("Unable to compute width: " + e.getMessage());
			}

			try
			{
				// Get z depth
				StandardResultProcedure sp = new StandardResultProcedure(results, DistanceUnit.PIXEL);
				sp.getXYZ();

				// Get precision
				PrecisionResultProcedure pp = new PrecisionResultProcedure(results);
				pp.getPrecision();
				stats[PRECISION].add(pp.precision);
				for (int i = 0; i < pp.size(); i++)
				{
					if (Math.abs(sp.z[i]) < depth)
						stats[PRECISION_IN_FOCUS].add(pp.precision[i]);
				}
			}
			catch (DataException e)
			{
				Utils.log("Unable to compute LSE precision: " + e.getMessage());
			}

			// Compute density per frame. Multithread for speed
			if (settings.getDensityRadius() > 0)
			{
				IJ.showStatus("Calculating density ...");
				final ExecutorService threadPool = Executors.newFixedThreadPool(Prefs.getThreads());
				final List<Future<?>> futures = new LinkedList<Future<?>>();
				final TFloatArrayList coordsX = new TFloatArrayList();
				final TFloatArrayList coordsY = new TFloatArrayList();
				final Statistics densityStats = stats[DENSITY];
				final float radius = (float) (settings.getDensityRadius() * getHWHM());
				final Rectangle bounds = results.getBounds();
				currentIndex = 0;
				finalIndex = results.getLastFrame();
				// Store the density for each result.
				final int[] allDensity = new int[results.size()];
				final FrameCounter counter = results.newFrameCounter();
				results.forEach(new PeakResultProcedure()
				{
					public void execute(PeakResult r)
					{
						if (counter.advance(r.getFrame()))
						{
							counter.increment(runDensityCalculation(threadPool, futures, coordsX, coordsY, densityStats,
									radius, bounds, allDensity, counter.getCount()));
						}
						coordsX.add(r.getXPosition());
						coordsY.add(r.getYPosition());
					}
				});
				runDensityCalculation(threadPool, futures, coordsX, coordsY, densityStats, radius, bounds, allDensity,
						counter.getCount());
				Utils.waitForCompletion(futures);
				threadPool.shutdownNow();
				IJ.showProgress(1);

				// Split results into singles (density = 0) and clustered (density > 0)
				final MemoryPeakResults singles = copyMemoryPeakResults("No Density");
				final MemoryPeakResults clustered = copyMemoryPeakResults("Density");

				counter.reset();
				results.forEach(new PeakResultProcedure()
				{
					public void execute(PeakResult r)
					{
						int density = allDensity[counter.getAndIncrement()];
						r.setOrigValue(density);
						if (density == 0)
							singles.add(r);
						else
							clustered.add(r);
					}
				});
			}
		}

		StringBuilder sb = new StringBuilder();
		sb.append(datasetNumber).append('\t');
		if (settings.getCameraType() == CameraType.SCMOS)
		{
			sb.append("sCMOS (").append(settings.getCameraModelName()).append(") ");
			Rectangle bounds = cameraModel.getBounds();
			sb.append(" ").append(bounds.x).append(",").append(bounds.y);
			int size = settings.getSize();
			sb.append(" ").append(size).append("x").append(size);
		}
		else if (CalibrationProtosHelper.isCCDCameraType(settings.getCameraType()))
		{
			sb.append(CalibrationProtosHelper.getName(settings.getCameraType()));
			int size = settings.getSize();
			sb.append(" ").append(size).append("x").append(size);
			if (settings.getCameraType() == CameraType.EMCCD)
				sb.append(" EM=").append(settings.getEmGain());
			sb.append(" CG=").append(settings.getCameraGain());
			sb.append(" RN=").append(settings.getReadNoise());
			sb.append(" B=").append(settings.getBias());
		}
		else
		{
			throw new IllegalStateException();
		}
		sb.append(" QE=").append(settings.getQuantumEfficiency()).append('\t');
		sb.append(settings.getPsfModel());
		if (psfModelType == PSF_MODEL_IMAGE)
		{
			sb.append(" Image").append(settings.getPsfImageName());
		}
		else if (psfModelType == PSF_MODEL_ASTIGMATISM)
		{
			sb.append(" model=").append(settings.getAstigmatismModel());
		}
		else
		{
			sb.append(" DoF=").append(Utils.rounded(settings.getDepthOfFocus()));
			if (settings.getEnterWidth())
			{
				sb.append(" SD=").append(Utils.rounded(settings.getPsfSd()));
			}
			else
			{
				sb.append(" λ=").append(Utils.rounded(settings.getWavelength()));
				sb.append(" NA=").append(Utils.rounded(settings.getNumericalAperture()));
			}
		}
		sb.append('\t');
		sb.append((fluorophores == null) ? localisations.size() : fluorophores.size()).append('\t');
		sb.append(stats[SAMPLED_BLINKS].getN() + (int) stats[SAMPLED_BLINKS].getSum()).append('\t');
		sb.append(localisations.size()).append('\t');
		sb.append(nFrames).append('\t');
		sb.append(Utils.rounded(areaInUm)).append('\t');
		sb.append(Utils.rounded(localisations.size() / (areaInUm * nFrames), 4)).append('\t');
		sb.append(Utils.rounded(getHWHM(), 4)).append('\t');
		double s = getPsfSD();
		sb.append(Utils.rounded(s, 4)).append('\t');
		s *= settings.getPixelPitch();
		final double sa = PSFCalculator.squarePixelAdjustment(s, settings.getPixelPitch()) / settings.getPixelPitch();
		sb.append(Utils.rounded(sa, 4)).append('\t');
		// Width not valid for the Image PSF.
		// Q. Is this true? We can approximate the FHWM for a spot-like image PSF.
		int nStats = (psfModelType == PSF_MODEL_IMAGE) ? stats.length - 1 : stats.length;
		for (int i = 0; i < nStats; i++)
		{
			double centre = (alwaysRemoveOutliers[i])
					? ((StoredDataStatistics) stats[i]).getStatistics().getPercentile(50) : stats[i].getMean();
			sb.append(Utils.rounded(centre, 4)).append('\t');
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
		if (settings.getShowHistograms())
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
							(integerDisplay[i]) ? 1 : 0,
							(settings.getRemoveOutliers() || alwaysRemoveOutliers[i]) ? 2 : 0,
							settings.getHistogramBins() * ((integerDisplay[i]) ? 100 : 1)));
					requireRetile = requireRetile || Utils.isNewWindow();
				}
			}

			wo.tile();
		}
		IJ.showStatus("");
		return stats[SIGNAL].getMean();
	}

	private int runDensityCalculation(ExecutorService threadPool, List<Future<?>> futures,
			final TFloatArrayList coordsX, final TFloatArrayList coordsY, final Statistics densityStats,
			final float radius, final Rectangle bounds, final int[] allDensity, final int allIndex)
	{
		final int size = coordsX.size();
		final float[] xCoords = coordsX.toArray();
		final float[] yCoords = coordsY.toArray();
		coordsX.resetQuick();
		coordsY.resetQuick();
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
		if (settings.getChooseHistograms())
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
				"Dataset\tCamera\tPSF\tMolecules\tPulses\tLocalisations\tnFrames\tArea (um^2)\tDensity (mol/um^2)\tHWHM\tS\tSa");
		for (int i = 0; i < NAMES.length; i++)
		{
			sb.append('\t').append(NAMES[i]);
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
		if (!settings.getSaveImage())
			return;
		String[] path = Utils.decodePath(settings.getImageFilename());
		OpenDialog chooser = new OpenDialog("Image_File", path[0], path[1]);
		if (chooser.getFileName() != null)
		{
			settings.setImageFilename(chooser.getDirectory() + chooser.getFileName());
			settings.setImageFilename(Utils.replaceExtension(settings.getImageFilename(), "tiff"));

			FileSaver fs = new FileSaver(imp);
			boolean ok;
			if (imp.getStackSize() > 1)
				ok = fs.saveAsTiffStack(settings.getImageFilename());
			else
				ok = fs.saveAsTiff(settings.getImageFilename());
			// The following call throws a NoSuchMethodError.
			// ok = IJ.saveAsTiff(imp, settings.imageFilename);

			if (!ok)
				IJ.log("Failed to save image to file: " + settings.getImageFilename());
		}
	}

	private void saveImageResults(MemoryPeakResults results)
	{
		if (!settings.getSaveImageResults())
			return;
		String[] path = Utils.decodePath(settings.getImageResultsFilename());
		OpenDialog chooser = new OpenDialog("Image_Results_File", path[0], path[1]);
		if (chooser.getFileName() != null)
		{
			settings.setImageResultsFilename(chooser.getDirectory() + chooser.getFileName());
			settings.setImageResultsFilename(Utils.replaceExtension(settings.getImageResultsFilename(), "xls"));

			TextFilePeakResults r = new TextFilePeakResults(settings.getImageResultsFilename(), false);
			r.copySettings(results);
			r.begin();
			r.addAll(results.toArray());
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
		float[] params = Gaussian2DPeakResultHelper.createTwoAxisParams(0, 0, 0, 0, 0, 0, 0);
		final int isx = Gaussian2DPeakResultHelper.INDEX_SX;
		final int isy = Gaussian2DPeakResultHelper.INDEX_SY;
		double noise = 0;
		int lastT = -1;
		for (LocalisationModel localisation : localisations)
		{
			if (currentId != localisation.getId() || lastT + 1 != localisation.getTime())
			{
				if (n > 0)
				{
					params[PeakResult.BACKGROUND] /= n;
					params[PeakResult.X] /= n;
					params[PeakResult.Y] /= n;
					params[isx] /= n;
					params[isy] /= n;

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
			params[PeakResult.BACKGROUND] += data[0];
			params[PeakResult.X] += localisation.getX();
			params[PeakResult.Y] += localisation.getY();
			params[PeakResult.INTENSITY] += localisation.getIntensity();
			noise += data[1] * data[1];
			params[isx] += data[2];
			params[isy] += data[3];
			n++;
			lastT = localisation.getTime();
		}

		// Final pulse
		if (n > 0)
		{
			params[PeakResult.BACKGROUND] /= n;
			params[PeakResult.X] /= n;
			params[PeakResult.Y] /= n;
			params[isx] /= n;
			params[isy] /= n;

			traceResults.add(new ExtendedPeakResult(start.getTime(), (int) Math.round(start.getX()),
					(int) Math.round(start.getY()), 0, 0, (float) (Math.sqrt(noise)), params, null, lastT, currentId));
		}

		traceResults.end();
	}

	private void saveFixedAndMoving(MemoryPeakResults results, String title)
	{
		if (simpleMode || benchmarkMode || spotMode)
			return;
		if (settings.getDiffusionRate() <= 0 || settings.getFixedFraction() >= 1)
			return;

		MemoryPeakResults fixedResults = copyMemoryPeakResults("Fixed");
		MemoryPeakResults movingResults = copyMemoryPeakResults("Moving");

		PeakResult[] peakResults = results.toArray();
		// Sort using the ID
		Arrays.sort(peakResults, new Comparator<PeakResult>()
		{
			public int compare(PeakResult o1, PeakResult o2)
			{
				return o1.getId() - o2.getId();
			}
		});

		MemoryPeakResults currentResults = movingResults;
		FrameCounter counter = new FrameCounter(-1);
		for (PeakResult p : peakResults)
		{
			if (counter.advance(p.getId()))
			{
				currentResults = (movingMolecules.contains(p.getId())) ? movingResults : fixedResults;
			}
			currentResults.add(p);
		}

		movingResults.end();
		fixedResults.end();
	}

	private void saveCompoundMolecules(MemoryPeakResults results, String title)
	{
		if (idToCompound == null)
			return;

		MemoryPeakResults[] set = new MemoryPeakResults[compoundNames.size()];
		for (int i = 0; i < set.length; i++)
			set[i] = copyMemoryPeakResults("Compound " + (i + 1) + ", " + compoundNames.get(i));

		PeakResult[] peakResults = results.toArray();
		// Sort using the ID
		Arrays.sort(peakResults, new Comparator<PeakResult>()
		{
			public int compare(PeakResult o1, PeakResult o2)
			{
				return o1.getId() - o2.getId();
			}
		});

		MemoryPeakResults currentResults = null;
		FrameCounter counter = new FrameCounter(-1);
		for (PeakResult p : peakResults)
		{
			if (counter.advance(p.getId()))
			{
				currentResults = set[idToCompound.get(p.getId())];
			}
			currentResults.add(p);
		}

		for (int i = 0; i < set.length; i++)
			set[i].end();
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
		if (!settings.getSaveFluorophores() || fluorophores == null)
			return;

		String[] path = Utils.decodePath(settings.getFluorophoresFilename());
		OpenDialog chooser = new OpenDialog("Fluorophores_File", path[0], path[1]);
		if (chooser.getFileName() != null)
		{
			settings.setFluorophoresFilename(chooser.getDirectory() + chooser.getFileName());
			settings.setFluorophoresFilename(Utils.replaceExtension(settings.getFluorophoresFilename(), "xls"));

			BufferedWriter output = null;
			try
			{
				output = new BufferedWriter(new FileWriter(settings.getFluorophoresFilename()));
				output.write(createResultsFileHeader());
				output.write("#Id\tn-Blinks\tStart\tStop\t...");
				output.newLine();
				for (int id = 1; id <= fluorophores.size(); id++)
				{
					FluorophoreSequenceModel f = fluorophores.get(id - 1);
					StringBuilder sb = new StringBuilder();
					sb.append(f.getId()).append('\t');
					sb.append(f.getNumberOfBlinks()).append('\t');
					for (double[] burst : f.getBurstSequence())
					{
						sb.append(Utils.rounded(burst[0], 3)).append('\t').append(Utils.rounded(burst[1], 3))
								.append('\t');
					}
					output.write(sb.toString());
					output.newLine();
				}
			}
			catch (Exception e)
			{
				// Q. Add better handling of errors?
				e.printStackTrace();
				IJ.log("Failed to save fluorophores to file: " + settings.getFluorophoresFilename());
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
		if (!settings.getSaveLocalisations())
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

		LoadLocalisationsSettings.Builder settings = SettingsManager.readLoadLocalisationsSettings(0).toBuilder();

		String[] path = Utils.decodePath(settings.getLocalisationsFilename());
		OpenDialog chooser = new OpenDialog("Localisations_File", path[0], path[1]);
		if (chooser.getFileName() != null)
		{
			settings.setLocalisationsFilename(
					Utils.replaceExtension(chooser.getDirectory() + chooser.getFileName(), "xls"));
			SettingsManager.writeSettings(settings.build());

			BufferedWriter output = null;
			try
			{
				output = new BufferedWriter(new FileWriter(settings.getLocalisationsFilename()));
				output.write(createResultsFileHeader());
				output.write("#T\tId\tX\tY\tZ\tIntensity");
				output.newLine();
				for (LocalisationModel l : localisations)
				{
					StringBuilder sb = new StringBuilder();
					sb.append(l.getTime()).append('\t');
					sb.append(l.getId()).append('\t');
					sb.append(l.getX()).append('\t');
					sb.append(l.getY()).append('\t');
					sb.append(l.getZ()).append('\t');
					sb.append(l.getIntensity());
					output.write(sb.toString());
					output.newLine();
				}
			}
			catch (Exception e)
			{
				// Q. Add better handling of errors?
				e.printStackTrace();
				IJ.log("Failed to save localisations to file: " + settings.getLocalisationsFilename());
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

			StringBuilder sb = new StringBuilder();
			sb.append("# ").append(TITLE).append(" Parameters:\n");
			addHeaderLine(sb, "Pixel_pitch (nm)", settings.getPixelPitch());
			addHeaderLine(sb, "Size", settings.getSize());
			if (!benchmarkMode)
			{
				addHeaderLine(sb, "Depth", settings.getDepth());
				addHeaderLine(sb, "Fixed depth", settings.getFixedDepth());
			}
			if (!(simpleMode || benchmarkMode))
			{
				if (!trackMode)
					addHeaderLine(sb, "Seconds", settings.getSeconds());
				addHeaderLine(sb, "Exposure_time", settings.getExposureTime());
				addHeaderLine(sb, "Steps_per_second", settings.getStepsPerSecond());
				if (!trackMode)
				{
					addHeaderLine(sb, "Illumination", settings.getIllumination());
					addHeaderLine(sb, "Pulse_interval", settings.getPulseInterval());
					addHeaderLine(sb, "Pulse_ratio", settings.getPulseRatio());
				}
				if (backgroundImages != null)
					addHeaderLine(sb, "Background_image", settings.getBackgroundImage());
			}
			addHeaderLine(sb, "Background", settings.getBackground());
			addCameraOptions(sb);
			addHeaderLine(sb, "PSF_model", settings.getPsfModel());
			if (psfModelType == PSF_MODEL_IMAGE)
			{
				addHeaderLine(sb, "PSF_image", settings.getPsfImageName());
			}
			else if (psfModelType == PSF_MODEL_ASTIGMATISM)
			{
				addHeaderLine(sb, "Astigmatism_model", settings.getAstigmatismModel());
				// Q. Should the actual model be appended?
				//addHeaderLine(sb, "Astigmatism_model parameters", SettingsManager.toJSON(getTheModel()));
			}
			else
			{
				addHeaderLine(sb, "Depth-of-focus (nm)", settings.getDepthOfFocus());
				if (settings.getEnterWidth())
				{
					addHeaderLine(sb, "PSF_SD", settings.getPsfSd());
				}
				else
				{
					addHeaderLine(sb, "Wavelength (nm)", settings.getWavelength());
					addHeaderLine(sb, "Numerical_aperture", settings.getNumericalAperture());
				}
			}
			if (!(benchmarkMode || spotMode))
			{
				addHeaderLine(sb, "Distribution", settings.getDistribution());
				if (settings.getDistribution().equals(DISTRIBUTION[MASK]))
				{
					addHeaderLine(sb, "Distribution_mask", settings.getDistributionMask());
				}
				else if (settings.getDistribution().equals(DISTRIBUTION[GRID]))
				{
					addHeaderLine(sb, "Cell_size", settings.getCellSize());
					addHeaderLine(sb, "p-binary", settings.getProbabilityBinary());
					addHeaderLine(sb, "Min_binary_distance (nm)", settings.getMinBinaryDistance());
					addHeaderLine(sb, "Max_binary_distance (nm)", settings.getMaxBinaryDistance());
				}
			}
			addHeaderLine(sb, "Particles", settings.getParticles());
			if (benchmarkMode)
			{
				addHeaderLine(sb, "X_position", settings.getXPosition());
				addHeaderLine(sb, "Y_position", settings.getYPosition());
				addHeaderLine(sb, "Z_position", settings.getZPosition());
				addHeaderLine(sb, "Min_photons", settings.getPhotonsPerSecond());
				addHeaderLine(sb, "Max_photons", settings.getPhotonsPerSecondMaximum());
			}
			else if (simpleMode)
			{
				addHeaderLine(sb, "Density (um^-2)", settings.getDensity());
				addHeaderLine(sb, "Min_photons", settings.getPhotonsPerSecond());
				addHeaderLine(sb, "Max_photons", settings.getPhotonsPerSecondMaximum());
			}
			else if (spotMode)
			{
				addHeaderLine(sb, "Min_photons", settings.getPhotonsPerSecond());
				addHeaderLine(sb, "Max_photons", settings.getPhotonsPerSecondMaximum());
			}
			else
			{
				addHeaderLine(sb, "Diffusion_rate", settings.getDiffusionRate());
				addHeaderLine(sb, "Diffusion_type", settings.getDiffusionType());
				addHeaderLine(sb, "Fixed_fraction", settings.getFixedFraction());
				if (settings.getCompoundMolecules())
				{
					addHeaderLine(sb, "Compound_molecules", settings.getCompoundText().replaceAll("\n *", ""));
					addHeaderLine(sb, "Enable_2D_diffusion", settings.getDiffuse2D());
					addHeaderLine(sb, "Rotate_initial_orientation", settings.getRotateInitialOrientation());
					addHeaderLine(sb, "Rotate_during_simulation", settings.getRotateDuringSimulation());
					addHeaderLine(sb, "Enable_2D_rotation", settings.getRotate2D());
				}
				addHeaderLine(sb, "Confinement", settings.getConfinement());
				if (settings.getConfinement().equals(CONFINEMENT[CONFINEMENT_SPHERE]))
				{
					addHeaderLine(sb, "Confinement_radius", settings.getConfinementRadius());
				}
				else if (settings.getConfinement().equals(CONFINEMENT[CONFINEMENT_MASK]))
				{
					addHeaderLine(sb, "Confinement_mask", settings.getConfinementMask());
				}
				addHeaderLine(sb, "Photon", settings.getPhotonsPerSecond());
				addHeaderLine(sb, "Photon_distribution", settings.getPhotonDistribution());
				if (PHOTON_DISTRIBUTION[PHOTON_CUSTOM].equals(settings.getPhotonDistribution()))
					addHeaderLine(sb, "Photon_distribution_file", settings.getPhotonDistributionFile());
				else if (PHOTON_DISTRIBUTION[PHOTON_UNIFORM].equals(settings.getPhotonDistribution()))
					addHeaderLine(sb, "Photon_max", settings.getPhotonsPerSecondMaximum());
				else if (PHOTON_DISTRIBUTION[PHOTON_GAMMA].equals(settings.getPhotonDistribution()))
					addHeaderLine(sb, "Photon_shape", settings.getPhotonShape());
				else if (PHOTON_DISTRIBUTION[PHOTON_CORRELATED].equals(settings.getPhotonDistribution()))
					addHeaderLine(sb, "Correlation", settings.getCorrelation());
				addHeaderLine(sb, "On_time", settings.getTOn());
				if (!trackMode)
				{
					addHeaderLine(sb, "Off_time_short", settings.getTOffShort());
					addHeaderLine(sb, "Off_time_long", settings.getTOffLong());
					addHeaderLine(sb, "n_Blinks_short", settings.getNBlinksShort());
					addHeaderLine(sb, "n_Blinks_long", settings.getNBlinksLong());
					addHeaderLine(sb, "n_Blinks_Geometric", settings.getNBlinksGeometricDistribution());
				}
				addHeaderLine(sb, "Min_photons", settings.getMinPhotons());
				addHeaderLine(sb, "Min_SNR_t1", settings.getMinSnrT1());
				addHeaderLine(sb, "Min_SNR_tN", settings.getMinSnrTN());
			}
			resultsFileHeader = sb.toString();
		}
		return resultsFileHeader;
	}

	private void addHeaderLine(StringBuilder sb, String name, Object o)
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
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

		settings = SettingsManager.readCreateDataSettings(0).toBuilder();

		// Image size
		gd.addMessage("--- Image Size ---");
		gd.addNumericField("Pixel_pitch (nm)", settings.getPixelPitch(), 2);
		gd.addNumericField("Size (px)", settings.getSize(), 0);
		if (!benchmarkMode)
		{
			gd.addNumericField("Depth (nm)", settings.getDepth(), 0);
			gd.addCheckbox("Fixed_depth", settings.getFixedDepth());
		}

		// Noise model
		gd.addMessage("--- Noise Model ---");
		if (extraOptions)
			gd.addCheckbox("No_poisson_noise", !settings.getPoissonNoise());
		gd.addNumericField("Background (photons)", settings.getBackground(), 2);

		addCameraOptions(gd);

		addPSFOptions(gd);

		gd.addMessage("--- Fluorophores ---");
		// Do not allow grid or mask distribution
		if (simpleMode)
		{
			// Allow mask but not the grid
			gd.addChoice("Distribution", Arrays.copyOf(DISTRIBUTION, DISTRIBUTION.length - 1),
					settings.getDistribution());
			gd.addCheckbox("Sample_per_frame", settings.getSamplePerFrame());
		}
		gd.addNumericField("Particles", settings.getParticles(), 0);
		if (simpleMode)
			gd.addNumericField("Density (um^-2)", settings.getDensity(), 2);
		else if (benchmarkMode)
		{
			gd.addNumericField("X_position (nm)", settings.getXPosition(), 2);
			gd.addNumericField("Y_position (nm)", settings.getYPosition(), 2);
			gd.addNumericField("Z_position (nm)", settings.getZPosition(), 2);
		}
		gd.addNumericField("Min_Photons", settings.getPhotonsPerSecond(), 0);
		gd.addNumericField("Max_Photons", settings.getPhotonsPerSecondMaximum(), 0);

		gd.addMessage("--- Save options ---");
		gd.addCheckbox("Raw_image", settings.getRawImage());
		gd.addCheckbox("Save_image", settings.getSaveImage());
		gd.addCheckbox("Save_image_results", settings.getSaveImageResults());
		gd.addCheckbox("Save_localisations", settings.getSaveLocalisations());

		gd.addMessage("--- Report options ---");
		gd.addCheckbox("Show_histograms", settings.getShowHistograms());
		gd.addCheckbox("Choose_histograms", settings.getChooseHistograms());
		gd.addNumericField("Histogram_bins", settings.getHistogramBins(), 0);
		gd.addCheckbox("Remove_outliers", settings.getRemoveOutliers());
		if (simpleMode)
			gd.addSlider("Density_radius (N x HWHM)", 0, 4.5, settings.getDensityRadius());
		gd.addNumericField("Depth-of-field (nm)", settings.getDepthOfField(), 0);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		settings.setPixelPitch(Math.abs(gd.getNextNumber()));
		settings.setSize(Math.abs((int) gd.getNextNumber()));
		if (!benchmarkMode)
		{
			// Allow negative depth
			settings.setDepth(gd.getNextNumber());
			settings.setFixedDepth(gd.getNextBoolean());
		}

		if (extraOptions)
		{
			settings.setPoissonNoise(!gd.getNextBoolean());
			poissonNoise = settings.getPoissonNoise();
		}
		settings.setBackground(Math.abs(gd.getNextNumber()));
		settings.setCameraType(SettingsManager.getCameraTypeValues()[gd.getNextChoiceIndex()]);

		settings.setPsfModel(gd.getNextChoice());

		if (simpleMode)
		{
			settings.setDistribution(gd.getNextChoice());
			settings.setSamplePerFrame(gd.getNextBoolean());
		}
		settings.setParticles(Math.abs((int) gd.getNextNumber()));
		if (simpleMode)
			settings.setDensity(Math.abs(gd.getNextNumber()));
		else if (benchmarkMode)
		{
			settings.setXPosition(gd.getNextNumber());
			settings.setYPosition(gd.getNextNumber());
			settings.setZPosition(gd.getNextNumber());
		}
		settings.setPhotonsPerSecond(Math.abs((int) gd.getNextNumber()));
		settings.setPhotonsPerSecondMaximum(Math.abs((int) gd.getNextNumber()));

		settings.setRawImage(gd.getNextBoolean());
		settings.setSaveImage(gd.getNextBoolean());
		settings.setSaveImageResults(gd.getNextBoolean());
		settings.setSaveLocalisations(gd.getNextBoolean());

		settings.setShowHistograms(gd.getNextBoolean());
		settings.setChooseHistograms(gd.getNextBoolean());
		settings.setHistogramBins((int) Math.abs(gd.getNextNumber()));
		settings.setRemoveOutliers(gd.getNextBoolean());
		if (simpleMode)
			settings.setDensityRadius((float) gd.getNextNumber());
		settings.setDepthOfField((float) Math.abs(gd.getNextNumber()));

		gd.collectOptions();

		// Save before validation so that the current values are preserved.
		SettingsManager.writeSettings(settings.build());

		if (gd.invalidNumber())
			return false;

		// Check arguments
		try
		{
			Parameters.isAboveZero("Pixel Pitch", settings.getPixelPitch());
			Parameters.isAboveZero("Size", settings.getSize());
			if (!benchmarkMode && !settings.getFixedDepth())
				Parameters.isPositive("Depth", settings.getDepth());
			Parameters.isPositive("Background", settings.getBackground());
			Parameters.isAboveZero("Particles", settings.getParticles());
			if (simpleMode)
				Parameters.isAboveZero("Density", settings.getDensity());
			Parameters.isAboveZero("Min Photons", settings.getPhotonsPerSecond());
			if (settings.getPhotonsPerSecondMaximum() < settings.getPhotonsPerSecond())
				settings.setPhotonsPerSecondMaximum(settings.getPhotonsPerSecond());
			Parameters.isPositive("Histogram bins", settings.getHistogramBins());
			if (simpleMode)
				Parameters.isPositive("Density radius", settings.getDensityRadius());

			validateCameraOptions();
			validatePSFOptions();
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		if (!benchmarkMode && settings.getDistribution().equals(DISTRIBUTION[MASK]))
		{
			String[] maskImages = createDistributionImageList();
			if (maskImages != null)
			{
				gd = new ExtendedGenericDialog(TITLE);
				gd.addMessage("Select the mask image for the distribution");
				gd.addChoice("Distribution_mask", maskImages, settings.getDistributionMask());
				if (maskListContainsStacks)
					gd.addNumericField("Distribution_slice_depth (nm)", settings.getDistributionMaskSliceDepth(), 0);
				gd.showDialog();
				if (gd.wasCanceled())
					return false;
				settings.setDistributionMask(gd.getNextChoice());
				if (maskListContainsStacks)
					settings.setDistributionMaskSliceDepth(Math.abs(gd.getNextNumber()));
			}

			SettingsManager.writeSettings(settings.build());
		}

		return getHistogramOptions();
	}

	private void addCameraOptions(final ExtendedGenericDialog gd)
	{
		gd.addChoice("Camera_type", SettingsManager.getCameraTypeNames(),
				CalibrationProtosHelper.getName(settings.getCameraType()), new OptionListener<Integer>()
				{
					public boolean collectOptions(Integer field)
					{
						settings.setCameraType(SettingsManager.getCameraTypeValues()[field]);
						return collectOptions(false);
					}

					public boolean collectOptions()
					{
						return collectOptions(true);
					}

					private boolean collectOptions(boolean silent)
					{
						CameraType cameraType = settings.getCameraType();
						boolean isCCD = CalibrationProtosHelper.isCCDCameraType(cameraType);
						ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE, null);
						if (isCCD)
						{
							if (cameraType == CameraType.EMCCD)
								egd.addNumericField("EM_gain", settings.getEmGain(), 2);
							egd.addNumericField("Camera_gain", settings.getCameraGain(), 4, 6, "count/electron");
							egd.addNumericField("Quantum_efficiency", settings.getQuantumEfficiency(), 2, 6,
									"electron/photon");
							egd.addNumericField("Read_noise", settings.getReadNoise(), 2, 6, "electron");
							egd.addNumericField("Bias", settings.getBias(), 0, 6, "count");
						}
						else if (cameraType == CameraType.SCMOS)
						{
							String[] models = CameraModelManager.listCameraModels(true);
							egd.addChoice("Camera_model_name", models, settings.getCameraModelName());
							egd.addNumericField("Quantum_efficiency", settings.getQuantumEfficiency(), 2, 6,
									"electron/photon");
						}
						else
						{
							IJ.error("Unsupported camera type " + CalibrationProtosHelper.getName(cameraType));
							return false;
						}
						egd.setSilent(silent);
						egd.showDialog(true, gd);
						if (egd.wasCanceled())
							return false;
						if (isCCD)
						{
							if (cameraType == CameraType.EMCCD)
								settings.setEmGain(Math.abs(egd.getNextNumber()));
							settings.setCameraGain(Math.abs(egd.getNextNumber()));
							settings.setQuantumEfficiency(Math.abs(egd.getNextNumber()));
							settings.setReadNoise(Math.abs(egd.getNextNumber()));
							settings.setBias(Math.abs((int) egd.getNextNumber()));
						}
						else if (cameraType == CameraType.SCMOS)
						{
							settings.setCameraModelName(egd.getNextChoice());
							settings.setQuantumEfficiency(Math.abs(egd.getNextNumber()));
						}
						return true;
					}
				});
	}

	private void validateCameraOptions()
	{
		CameraType cameraType = settings.getCameraType();
		boolean isCCD = CalibrationProtosHelper.isCCDCameraType(cameraType);
		if (isCCD)
		{
			if (cameraType == CameraType.EMCCD)
				Parameters.isPositive("EM gain", settings.getEmGain());
			Parameters.isPositive("Camera gain", settings.getCameraGain());
			Parameters.isPositive("Read noise", settings.getReadNoise());
			double noiseRange = settings.getReadNoise() * settings.getCameraGain() * 4;
			Parameters.isEqualOrAbove("Bias must prevent clipping the read noise (@ +/- 4 StdDev) so ",
					settings.getBias(), noiseRange);

			cameraModel = createCCDCameraModel();
		}
		else if (cameraType == CameraType.SCMOS)
		{
			// Load the model
			cameraModel = CameraModelManager.load(settings.getCameraModelName());
			if (cameraModel == null)
			{
				throw new IllegalArgumentException("Unknown camera model for name: " + settings.getCameraModelName());
			}

			// Check the width is above the selected size
			Rectangle modelBounds = cameraModel.getBounds();
			int size = settings.getSize();
			if (modelBounds.width < size || modelBounds.height < size)
			{
				throw new IllegalArgumentException(String.format(
						"Camera model bounds [x=%d,y=%d,width=%d,height=%d] are smaller than simulation size [%d]",
						modelBounds.x, modelBounds.y, modelBounds.width, modelBounds.height, size));
			}

			// Ask for a crop
			if (modelBounds.width > size || modelBounds.height > size)
			{
				GenericDialog gd = new GenericDialog(TITLE);
				//@formatter:off
				gd.addMessage(String.format(
						"WARNING:\n \nCamera model bounds\n[x=%d,y=%d,width=%d,height=%d]\nare larger than the simulation size [=%d].\n \nCrop the model?",
						modelBounds.x, modelBounds.y, modelBounds.width, modelBounds.height, size
						));
				//@formatter:on
				gd.addCheckbox("Random_crop", settings.getRandomCrop());
				int upperx = modelBounds.x + modelBounds.width - size;
				int uppery = modelBounds.y + modelBounds.height - size;
				gd.addSlider("Origin_x", modelBounds.x, upperx,
						Maths.clip(modelBounds.x, upperx, settings.getOriginX()));
				gd.addSlider("Origin_y", modelBounds.y, uppery,
						Maths.clip(modelBounds.y, uppery, settings.getOriginY()));
				gd.showDialog();
				if (gd.wasCanceled())
					throw new IllegalArgumentException("Unknown camera model crop");
				settings.setRandomCrop(gd.getNextBoolean());
				settings.setOriginX((int) gd.getNextNumber());
				settings.setOriginY((int) gd.getNextNumber());
				SettingsManager.writeSettings(settings.build());

				int ox, oy;
				if (settings.getRandomCrop())
				{
					RandomDataGenerator rg = new RandomDataGenerator(createRandomGenerator());
					ox = rg.nextInt(modelBounds.x, upperx);
					oy = rg.nextInt(modelBounds.y, uppery);
				}
				else
				{
					ox = settings.getOriginX();
					oy = settings.getOriginY();
				}
				Rectangle bounds = new Rectangle(ox, oy, size, size);
				cameraModel = cameraModel.crop(bounds, false);
				modelBounds = cameraModel.getBounds();
				if (modelBounds.width != size || modelBounds.height != size)
					throw new IllegalArgumentException("Failed to crop camera model to bounds: " + bounds);
			}
		}
		else
		{
			throw new IllegalArgumentException(
					"Unsupported camera type: " + CalibrationProtosHelper.getName(cameraType));
		}
	}

	private void addCameraOptions(StringBuilder sb)
	{
		CameraType cameraType = settings.getCameraType();
		boolean isCCD = CalibrationProtosHelper.isCCDCameraType(cameraType);
		if (isCCD)
		{
			if (cameraType == CameraType.EMCCD)
				addHeaderLine(sb, "EM_gain", settings.getEmGain());
			addHeaderLine(sb, "Camera_gain", settings.getCameraGain());
			addHeaderLine(sb, "Quantum_efficiency", getQuantumEfficiency());
			addHeaderLine(sb, "Read_noise", settings.getReadNoise());
			addHeaderLine(sb, "Bias", settings.getBias());
		}
		else if (cameraType == CameraType.SCMOS)
		{
			addHeaderLine(sb, "Camera_model_name", settings.getCameraModelName());
			addHeaderLine(sb, "Quantum_efficiency", getQuantumEfficiency());
		}
	}

	private double getQuantumEfficiency()
	{
		double qe = settings.getQuantumEfficiency();
		return (qe > 0 && qe < 1) ? qe : 1;
	}

	/**
	 * Check if there are any suitable PSF images open. If so add a choice to allow the selection of the Gaussian or
	 * Image PSF model. If no PSF images are open then add options for the wavelength and NA for the simulated
	 * microscope.
	 * 
	 * @param gd
	 * @return
	 */
	private void addPSFOptions(final ExtendedGenericDialog gd)
	{
		gd.addMessage("--- PSF Model ---");
		List<String> imageNames = PSFCombiner.createImageList();
		TurboList<String> availableModels = new TurboList<String>();
		availableModels.add(PSF_MODELS[PSF_MODEL_GAUSSIAN]);
		availableModels.add(PSF_MODELS[PSF_MODEL_AIRY]);
		final String[] images;
		if (!imageNames.isEmpty())
		{
			availableModels.add(PSF_MODELS[PSF_MODEL_IMAGE]);
			images = imageNames.toArray(new String[imageNames.size()]);
		}
		else
		{
			images = null;
		}
		final String[] astigmatismModels = AstigmatismModelManager.listAstigmatismModels(false, true);
		if (astigmatismModels.length != 0)
			availableModels.add(PSF_MODELS[PSF_MODEL_ASTIGMATISM]);
		final String[] models = availableModels.toArray(new String[availableModels.size()]);
		gd.addChoice("PSF_model", models, settings.getPsfModel(), new OptionListener<Integer>()
		{
			public boolean collectOptions(Integer value)
			{
				settings.setPsfModel(models[value]);
				return collectOptions(false);
			}

			public boolean collectOptions()
			{
				return collectOptions(true);
			}

			private boolean collectOptions(boolean silent)
			{
				ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE, null);
				egd.addMessage("Configure the " + settings.getPsfModel() + " PSF model");

				int type = 0;

				// Get the image
				if (settings.getPsfModel().equals(PSF_MODELS[PSF_MODEL_IMAGE]))
				{
					egd.addChoice("PSF_image", images, settings.getPsfImageName());
				}
				// Get the astigmatism model
				else if (settings.getPsfModel().equals(PSF_MODELS[PSF_MODEL_ASTIGMATISM]))
				{
					type = 1;
					egd.addChoice("Astigmatism_model", astigmatismModels, settings.getAstigmatismModel());
					egd.addMessage(TextUtils.wrap("Note: The pixel size of the astigmatism model should match " +
							"the pixel pitch if fitting of the data is to be performed (i.e. fitting requires the " +
							"astigmatism model to be calibrated to the image). If not then the model will be " +
							"optionally converted before the simulation.", 80));
				}
				// Get the width of the model
				else
				{
					type = 2;
					egd.addNumericField("Depth-of-focus (nm)", settings.getDepthOfFocus(), 2);
					egd.addCheckbox("Enter_width", settings.getEnterWidth());
					egd.addNumericField("PSF_SD (nm)", settings.getPsfSd(), 2);
					egd.addMessage("Or compute from optics:");
					egd.addNumericField("Wavelength (nm)", settings.getWavelength(), 2);
					egd.addNumericField("Numerical_aperture", settings.getNumericalAperture(), 2);
				}
				egd.showDialog(true, gd);
				if (egd.wasCanceled())
					return false;
				if (type == 0)
				{
					settings.setPsfImageName(egd.getNextChoice());
				}
				else if (type == 1)
				{
					settings.setAstigmatismModel(AstigmatismModelManager.removeFormatting(egd.getNextChoice()));
				}
				else
				{
					settings.setDepthOfFocus(egd.getNextNumber());
					settings.setEnterWidth(egd.getNextBoolean());
					settings.setPsfSd(Math.abs(egd.getNextNumber()));
					settings.setWavelength(Math.abs(egd.getNextNumber()));
					settings.setNumericalAperture(Math.abs(egd.getNextNumber()));
				}
				return true;
			}

		});
	}

	private void validatePSFOptions()
	{
		if (settings.getPsfModel().equals(PSF_MODELS[PSF_MODEL_ASTIGMATISM]))
		{
			psfModelType = PSF_MODEL_ASTIGMATISM;
			AstigmatismModel model = AstigmatismModelManager.getModel(settings.getAstigmatismModel());
			if (model == null)
				throw new IllegalArgumentException("Failed to load model: " + settings.getAstigmatismModel());
		}
		else if (settings.getPsfModel().equals(PSF_MODELS[PSF_MODEL_IMAGE]))
		{
			psfModelType = PSF_MODEL_IMAGE;
		}
		else
		{
			psfModelType = (settings.getPsfModel().equals(PSF_MODELS[PSF_MODEL_GAUSSIAN])) ? PSF_MODEL_GAUSSIAN
					: PSF_MODEL_AIRY;
			Parameters.isAboveZero("Depth-of-focus", settings.getDepthOfFocus());
			if (settings.getEnterWidth())
			{
				Parameters.isAboveZero("PSF SD", settings.getPsfSd());
			}
			else
			{
				Parameters.isAboveZero("Wavelength", settings.getWavelength());
				Parameters.isAboveZero("NA", settings.getNumericalAperture());
				Parameters.isBelow("NA", settings.getNumericalAperture(), 2);
			}
		}
	}

	private boolean getHistogramOptions()
	{
		GenericDialog gd;
		if (settings.getShowHistograms() && settings.getChooseHistograms())
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

		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

		settings = SettingsManager.readCreateDataSettings(0).toBuilder();

		if (settings.getStepsPerSecond() < 1)
			settings.setStepsPerSecond(1);

		String[] backgroundImages = createBackgroundImageList();
		gd.addNumericField("Pixel_pitch (nm)", settings.getPixelPitch(), 2);
		gd.addNumericField("Size (px)", settings.getSize(), 0);
		gd.addNumericField("Depth (nm)", settings.getDepth(), 0);
		gd.addCheckbox("Fixed_depth", settings.getFixedDepth());
		if (!trackMode)
			gd.addNumericField("Seconds", settings.getSeconds(), 1);
		gd.addNumericField("Exposure_time (ms)", settings.getExposureTime(), 1);
		gd.addSlider("Steps_per_second", 1, Maths.clip(15, 1000, settings.getStepsPerSecond()),
				settings.getStepsPerSecond());
		if (!trackMode)
		{
			gd.addChoice("Illumination", ILLUMINATION, settings.getIllumination());
			gd.addNumericField("Pulse_interval", settings.getPulseInterval(), 0);
			gd.addNumericField("Pulse_ratio", settings.getPulseRatio(), 2);
		}
		if (backgroundImages != null)
			gd.addChoice("Background_image", backgroundImages, settings.getBackgroundImage());

		if (extraOptions)
			gd.addCheckbox("No_poisson_noise", !settings.getPoissonNoise());
		gd.addNumericField("Background (photons)", settings.getBackground(), 2);

		addCameraOptions(gd);

		addPSFOptions(gd);

		gd.addMessage("--- Fluorophores ---");
		gd.addChoice("Distribution", DISTRIBUTION, settings.getDistribution());
		gd.addNumericField("Particles", settings.getParticles(), 0);
		gd.addCheckbox("Compound_molecules", settings.getCompoundMolecules());
		gd.addNumericField("Diffusion_rate (um^2/sec)", settings.getDiffusionRate(), 2);
		String[] diffusionTypes = SettingsManager.getNames((Object[]) DiffusionType.values());
		gd.addChoice("Diffusion_type", diffusionTypes,
				diffusionTypes[CreateDataSettingsHelper.getDiffusionType(settings.getDiffusionType()).ordinal()]);
		gd.addSlider("Fixed_fraction (%)", 0, 100, settings.getFixedFraction() * 100);
		gd.addChoice("Confinement", CONFINEMENT, settings.getConfinement());
		gd.addNumericField("Photons (sec^-1)", settings.getPhotonsPerSecond(), 0);
		// We cannot use the correlation moe with fixed life time tracks 
		String[] dist = (trackMode) ? Arrays.copyOf(PHOTON_DISTRIBUTION, PHOTON_DISTRIBUTION.length - 1)
				: PHOTON_DISTRIBUTION;
		gd.addChoice("Photon_distribution", dist, settings.getPhotonDistribution());
		gd.addNumericField("On_time (ms)", settings.getTOn(), 2);
		if (!trackMode)
		{
			gd.addNumericField("Off_time_short (ms)", settings.getTOffShort(), 2);
			gd.addNumericField("Off_time_long (ms)", settings.getTOffLong(), 2);
			gd.addNumericField("n_Blinks_Short", settings.getNBlinksShort(), 2);
			gd.addNumericField("n_Blinks_Long", settings.getNBlinksLong(), 2);
			gd.addCheckbox("Use_geometric_distribution", settings.getNBlinksGeometricDistribution());
		}

		gd.addMessage("--- Peak filtering ---");
		gd.addSlider("Min_Photons", 0, 50, settings.getMinPhotons());
		gd.addSlider("Min_SNR_t1", 0, 20, settings.getMinSnrT1());
		gd.addSlider("Min_SNR_tN", 0, 10, settings.getMinSnrTN());

		gd.addMessage("--- Save options ---");
		gd.addCheckbox("Raw_image", settings.getRawImage());
		gd.addCheckbox("Save_image", settings.getSaveImage());
		gd.addCheckbox("Save_image_results", settings.getSaveImageResults());
		gd.addCheckbox("Save_fluorophores", settings.getSaveFluorophores());
		gd.addCheckbox("Save_localisations", settings.getSaveLocalisations());

		gd.addMessage("--- Report options ---");
		gd.addCheckbox("Show_histograms", settings.getShowHistograms());
		gd.addCheckbox("Choose_histograms", settings.getChooseHistograms());
		gd.addNumericField("Histogram_bins", settings.getHistogramBins(), 0);
		gd.addCheckbox("Remove_outliers", settings.getRemoveOutliers());
		gd.addSlider("Density_radius (N x HWHM)", 0, 4.5, settings.getDensityRadius());
		gd.addNumericField("Depth-of-field (nm)", settings.getDepthOfField(), 0);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		settings.setPixelPitch(Math.abs(gd.getNextNumber()));
		settings.setSize(Math.abs((int) gd.getNextNumber()));
		settings.setDepth(Math.abs(gd.getNextNumber()));
		settings.setFixedDepth(gd.getNextBoolean());
		if (!trackMode)
			settings.setSeconds(Math.abs(gd.getNextNumber()));
		settings.setExposureTime(Math.abs(gd.getNextNumber()));
		settings.setStepsPerSecond(Math.abs(gd.getNextNumber()));
		if (!trackMode)
		{
			settings.setIllumination(gd.getNextChoice());
			settings.setPulseInterval(Math.abs((int) gd.getNextNumber()));
			settings.setPulseRatio(Math.abs(gd.getNextNumber()));
		}
		if (backgroundImages != null)
			settings.setBackgroundImage(gd.getNextChoice());

		if (extraOptions)
		{
			settings.setPoissonNoise(!gd.getNextBoolean());
			poissonNoise = settings.getPoissonNoise();
		}
		settings.setBackground(Math.abs(gd.getNextNumber()));
		settings.setCameraType(SettingsManager.getCameraTypeValues()[gd.getNextChoiceIndex()]);

		settings.setPsfModel(gd.getNextChoice());

		settings.setDistribution(gd.getNextChoice());
		settings.setParticles(Math.abs((int) gd.getNextNumber()));
		settings.setCompoundMolecules(gd.getNextBoolean());
		settings.setDiffusionRate(Math.abs(gd.getNextNumber()));
		settings.setDiffusionType(gd.getNextChoiceIndex());
		settings.setFixedFraction(Math.abs(gd.getNextNumber() / 100.0));
		settings.setConfinement(gd.getNextChoice());
		settings.setPhotonsPerSecond(Math.abs((int) gd.getNextNumber()));
		settings.setPhotonDistribution(gd.getNextChoice());
		settings.setTOn(Math.abs(gd.getNextNumber()));
		if (!trackMode)
		{
			settings.setTOffShort(Math.abs(gd.getNextNumber()));
			settings.setTOffLong(Math.abs(gd.getNextNumber()));
			settings.setNBlinksShort(Math.abs(gd.getNextNumber()));
			settings.setNBlinksLong(Math.abs(gd.getNextNumber()));
			settings.setNBlinksGeometricDistribution(gd.getNextBoolean());
		}

		settings.setMinPhotons(gd.getNextNumber());
		settings.setMinSnrT1(gd.getNextNumber());
		settings.setMinSnrTN(gd.getNextNumber());
		minPhotons = settings.getMinPhotons();
		minSNRt1 = settings.getMinSnrT1();
		minSNRtN = settings.getMinSnrTN();

		settings.setRawImage(gd.getNextBoolean());
		settings.setSaveImage(gd.getNextBoolean());
		settings.setSaveImageResults(gd.getNextBoolean());
		settings.setSaveFluorophores(gd.getNextBoolean());
		settings.setSaveLocalisations(gd.getNextBoolean());

		settings.setShowHistograms(gd.getNextBoolean());
		settings.setChooseHistograms(gd.getNextBoolean());
		settings.setHistogramBins((int) gd.getNextNumber());
		settings.setRemoveOutliers(gd.getNextBoolean());
		settings.setDensityRadius((float) gd.getNextNumber());
		settings.setDepthOfField((float) Math.abs(gd.getNextNumber()));

		// Ensure tN threshold is more lenient
		if (settings.getMinSnrT1() < settings.getMinSnrTN())
		{
			double tmp = settings.getMinSnrT1();
			settings.setMinSnrT1(settings.getMinSnrTN());
			settings.setMinSnrTN(tmp);
		}

		gd.collectOptions();

		// Save before validation so that the current values are preserved.
		SettingsManager.writeSettings(settings.build());

		// Check arguments
		try
		{
			Parameters.isAboveZero("Pixel Pitch", settings.getPixelPitch());
			Parameters.isAboveZero("Size", settings.getSize());
			if (!settings.getFixedDepth())
				Parameters.isPositive("Depth", settings.getDepth());
			if (!trackMode)
				Parameters.isAboveZero("Seconds", settings.getSeconds());
			Parameters.isAboveZero("Exposure time", settings.getExposureTime());
			Parameters.isAboveZero("Steps per second", settings.getStepsPerSecond());
			Parameters.isPositive("Background", settings.getBackground());
			Parameters.isAboveZero("Particles", settings.getParticles());
			Parameters.isAboveZero("Photons", settings.getPhotonsPerSecond());
			Parameters.isPositive("Diffusion rate", settings.getDiffusionRate());
			Parameters.isPositive("Fixed fraction", settings.getFixedFraction());
			Parameters.isPositive("Pulse interval", settings.getPulseInterval());
			Parameters.isAboveZero("Pulse ratio", settings.getPulseRatio());
			Parameters.isAboveZero("tOn", settings.getTOn());
			if (!trackMode)
			{
				Parameters.isAboveZero("tOff Short", settings.getTOffShort());
				Parameters.isAboveZero("tOff Long", settings.getTOffLong());
				Parameters.isPositive("n-Blinks Short", settings.getNBlinksShort());
				Parameters.isPositive("n-Blinks Long", settings.getNBlinksLong());
			}
			Parameters.isPositive("Min photons", settings.getMinPhotons());
			Parameters.isPositive("Min SNR t1", settings.getMinSnrT1());
			Parameters.isPositive("Min SNR tN", settings.getMinSnrTN());
			Parameters.isPositive("Histogram bins", settings.getHistogramBins());
			Parameters.isPositive("Density radius", settings.getDensityRadius());

			validateCameraOptions();
			validatePSFOptions();
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
		if (settings.getDistribution().equals(DISTRIBUTION[MASK]))
		{
			maskImages = createDistributionImageList();
			if (maskImages != null)
			{
				gd = new ExtendedGenericDialog(TITLE);
				gd.addMessage("Select the mask image for the distribution");
				gd.addChoice("Distribution_mask", maskImages, settings.getDistributionMask());
				if (maskListContainsStacks)
					gd.addNumericField("Distribution_slice_depth (nm)", settings.getDistributionMaskSliceDepth(), 0);
				gd.showDialog();
				if (gd.wasCanceled())
					return false;
				settings.setDistributionMask(gd.getNextChoice());
				if (maskListContainsStacks)
					settings.setDistributionMaskSliceDepth(Math.abs(gd.getNextNumber()));
			}
		}
		else if (settings.getDistribution().equals(DISTRIBUTION[GRID]))
		{
			gd = new ExtendedGenericDialog(TITLE);
			gd.addMessage("Select grid for the distribution");
			gd.addNumericField("Cell_size", settings.getCellSize(), 0);
			gd.addSlider("p-binary", 0, 1, settings.getProbabilityBinary());
			gd.addNumericField("Min_binary_distance (nm)", settings.getMinBinaryDistance(), 0);
			gd.addNumericField("Max_binary_distance (nm)", settings.getMaxBinaryDistance(), 0);
			gd.showDialog();
			if (gd.wasCanceled())
				return false;
			settings.setCellSize((int) gd.getNextNumber());
			settings.setProbabilityBinary(gd.getNextNumber());
			settings.setMinBinaryDistance(gd.getNextNumber());
			settings.setMaxBinaryDistance(gd.getNextNumber());

			// Check arguments
			try
			{
				Parameters.isAboveZero("Cell size", settings.getCellSize());
				Parameters.isPositive("p-binary", settings.getProbabilityBinary());
				Parameters.isEqualOrBelow("p-binary", settings.getProbabilityBinary(), 1);
				Parameters.isPositive("Min binary distance", settings.getMinBinaryDistance());
				Parameters.isPositive("Max binary distance", settings.getMaxBinaryDistance());
				Parameters.isEqualOrBelow("Min binary distance", settings.getMinBinaryDistance(),
						settings.getMaxBinaryDistance());
			}
			catch (IllegalArgumentException e)
			{
				IJ.error(TITLE, e.getMessage());
				return false;
			}
		}

		SettingsManager.writeSettings(settings.build());

		if (settings.getDiffusionRate() > 0 && settings.getFixedFraction() < 1)
		{
			if (settings.getConfinement().equals(CONFINEMENT[CONFINEMENT_SPHERE]))
			{
				gd = new ExtendedGenericDialog(TITLE);
				gd.addMessage("Select the sphere radius for the diffusion confinement");
				gd.addSlider("Confinement_radius (nm)", 0, 2000, settings.getConfinementRadius());
				gd.showDialog();
				if (gd.wasCanceled())
					return false;
				settings.setConfinementRadius(gd.getNextNumber());
			}
			else if (settings.getConfinement().equals(CONFINEMENT[CONFINEMENT_MASK]))
			{
				if (maskImages == null)
					maskImages = createDistributionImageList();
				if (maskImages != null)
				{
					gd = new ExtendedGenericDialog(TITLE);
					gd.addMessage("Select the mask image for the diffusion confinement");
					gd.addChoice("Confinement_mask", maskImages, settings.getConfinementMask());
					if (maskListContainsStacks)
						gd.addNumericField("Confinement_slice_depth (nm)", settings.getConfinementMaskSliceDepth(), 0);
					gd.showDialog();
					if (gd.wasCanceled())
						return false;
					settings.setConfinementMask(gd.getNextChoice());
					if (maskListContainsStacks)
						settings.setConfinementMaskSliceDepth(Math.abs(gd.getNextNumber()));
				}
			}
		}

		SettingsManager.writeSettings(settings.build());

		if (settings.getCompoundMolecules())
		{
			// Show a second dialog where the molecule configuration is specified
			gd = new ExtendedGenericDialog(TITLE);

			gd.addMessage("Specify the compound molecules");
			gd.addTextAreas(settings.getCompoundText(), null, 20, 80);
			gd.addCheckbox("Enable_2D_diffusion", settings.getDiffuse2D());
			gd.addCheckbox("Rotate_initial_orientation", settings.getRotateInitialOrientation());
			gd.addCheckbox("Rotate_during_simulation", settings.getRotateDuringSimulation());
			gd.addCheckbox("Enable_2D_rotation", settings.getRotate2D());
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

			settings.setCompoundText(gd.getNextText());
			settings.setDiffuse2D(gd.getNextBoolean());
			settings.setRotateInitialOrientation(gd.getNextBoolean());
			settings.setRotateDuringSimulation(gd.getNextBoolean());
			settings.setRotate2D(gd.getNextBoolean());

			if (gd.getNextBoolean())
			{
				logExampleCompounds();
				return false;
			}
		}

		SettingsManager.writeSettings(settings.build());

		gd = new ExtendedGenericDialog(TITLE);
		gd.addMessage("Configure the photon distribution: " + settings.getPhotonDistribution());
		if (PHOTON_DISTRIBUTION[PHOTON_CUSTOM].equals(settings.getPhotonDistribution()))
		{
			// Nothing more to be done
			return true;
		}
		else if (PHOTON_DISTRIBUTION[PHOTON_UNIFORM].equals(settings.getPhotonDistribution()))
		{
			gd.addNumericField("Max_Photons (sec^-1)", settings.getPhotonsPerSecondMaximum(), 0);
		}
		else if (PHOTON_DISTRIBUTION[PHOTON_GAMMA].equals(settings.getPhotonDistribution()))
		{
			gd.addNumericField("Photon_shape", settings.getPhotonShape(), 2);
		}
		else if (PHOTON_DISTRIBUTION[PHOTON_CORRELATED].equals(settings.getPhotonDistribution()))
		{
			gd.addNumericField("Correlation (to total tOn)", settings.getCorrelation(), 2);
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
			if (PHOTON_DISTRIBUTION[PHOTON_UNIFORM].equals(settings.getPhotonDistribution()))
			{
				settings.setPhotonsPerSecondMaximum(Math.abs((int) gd.getNextNumber()));
				if (settings.getPhotonsPerSecondMaximum() < settings.getPhotonsPerSecond())
					settings.setPhotonsPerSecondMaximum(settings.getPhotonsPerSecond());
			}
			else if (PHOTON_DISTRIBUTION[PHOTON_GAMMA].equals(settings.getPhotonDistribution()))
			{
				settings.setPhotonShape(Math.abs(gd.getNextNumber()));
				Parameters.isAbove("Photon shape", settings.getPhotonShape(), 0);
			}
			else if (PHOTON_DISTRIBUTION[PHOTON_CORRELATED].equals(settings.getPhotonDistribution()))
			{
				settings.setCorrelation(gd.getNextNumber());
				Parameters.isEqualOrBelow("Correlation", settings.getCorrelation(), 1);
				Parameters.isEqualOrAbove("Correlation", settings.getCorrelation(), -1);
			}
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		SettingsManager.writeSettings(settings.build());

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
		comment("Compounds are described using parseable text.");
		comment("Missing fields are initialised to the default (0).");
		comment("Multiple compounds can be combined using fractional ratios.");
		comment("Coordinates are specified in nanometres.");
		comment("Coordinates describe the relative positions of atoms in the molecule; the molcule will have a randomly assigned XYZ position for its centre-of-mass. Rotation will be about the centre-of-mass.");
		IJ.log("");

		Molecule.Builder mb = Molecule.newBuilder();
		addAtom(mb, 10, 1, 1, 1);
		mb.setDiffusionRate(0.5);
		mb.setDiffusionType(DiffusionType.RANDOM_WALK.toString());
		Molecule m1 = mb.build();

		mb.clear();
		addAtom(mb, 30, 0, 0, 0);
		addAtom(mb, 20, 1000, 0, 0);
		mb.setDiffusionRate(1);
		mb.setDiffusionType(DiffusionType.GRID_WALK.toString());
		Molecule m2 = mb.build();

		// Create a hexamer big enough to see with the default pixel pitch
		mb.clear();
		addAtom(mb, 1, 0, 0, 0);
		addAtom(mb, 1, 1000, 0, 0);
		addAtom(mb, 1, 1500, 866, 0);
		addAtom(mb, 1, 1000, 1732, 0);
		addAtom(mb, 1, 0, 1732, 0);
		addAtom(mb, 1, -500, 866, 0);
		Molecule m3 = mb.build();

		comment("Single molecules");
		IJ.log("");
		comment("Moving Monomer");
		demo(m1);
		comment("Moving Dimer");
		demo(m2);
		comment("Fixed Hexamer");
		demo(m3);

		comment("Mixtures of molecules");
		IJ.log("");
		comment("Two molecules with a ratio of 2:1");
		m1 = m1.toBuilder().setFraction(2).build();
		m2 = m2.toBuilder().setFraction(1).build();
		demo(m1, m2);
	}

	private void addAtom(Molecule.Builder mb, double mass, double x, double y, double z)
	{
		Atom.Builder ab = mb.addAtomBuilder();
		ab.setMass(mass);
		ab.setX(x);
		ab.setY(y);
		ab.setZ(z);
	}

	private void demo(Molecule... molecules)
	{
		Mixture.Builder builder = Mixture.newBuilder();
		for (Molecule m : molecules)
			builder.addMolecule(m);

		// The toString() method is more verbose than JSON but easier to read
		IJ.log(builder.toString());
		IJ.log("");
	}

	private void comment(String text)
	{
		IJ.log(TextUtils.wrap("# " + text, 80, "\n# ", false));
	}

	private List<CompoundMoleculeModel> createCompoundMolecules()
	{
		// Diffusion rate is um^2/sec. Convert to pixels per simulation frame.
		final double diffusionFactor = (1000000.0 / (settings.getPixelPitch() * settings.getPixelPitch())) /
				settings.getStepsPerSecond();

		List<CompoundMoleculeModel> compounds;
		if (settings.getCompoundMolecules())
		{
			// Try and load the compounds from the text specification
			try
			{
				// Convert from the serialised objects to the compound model
				String text = settings.getCompoundText();
				Mixture.Builder builder = Mixture.newBuilder();
				TextFormat.merge(text, builder);

				compounds = new ArrayList<CompoundMoleculeModel>(builder.getMoleculeCount());
				int id = 1;
				compoundNames = new ArrayList<String>(builder.getMoleculeCount());
				for (Molecule m : builder.getMoleculeList())
				{
					MoleculeModel[] molecules = new MoleculeModel[m.getAtomCount()];
					for (int i = 0; i < molecules.length; i++)
					{
						AtomOrBuilder a = m.getAtomOrBuilder(i);
						molecules[i] = new MoleculeModel(a.getMass(), a.getX(), a.getY(), a.getZ());
					}
					CompoundMoleculeModel cm = new CompoundMoleculeModel(id++, 0, 0, 0, Arrays.asList(molecules));
					cm.setFraction(m.getFraction());
					cm.setDiffusionRate(m.getDiffusionRate() * diffusionFactor);
					cm.setDiffusionType(DiffusionType.fromString(m.getDiffusionType()));
					compounds.add(cm);
					compoundNames.add(String.format("Fraction=%s, D=%s um^2/s", Utils.rounded(cm.getFraction()),
							Utils.rounded(m.getDiffusionRate())));
				}

				// Convert coordinates from nm to pixels
				final double scaleFactor = 1.0 / settings.getPixelPitch();
				for (CompoundMoleculeModel c : compounds)
				{
					c.scale(scaleFactor);
				}
			}
			catch (IOException e)
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
			m.setDiffusionRate(settings.getDiffusionRate() * diffusionFactor);
			m.setDiffusionType(CreateDataSettingsHelper.getDiffusionType(settings.getDiffusionType()));
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
	public static ImmutableMemoryPeakResults getResults()
	{
		MemoryPeakResults r = MemoryPeakResults.getResults(benchmarkResultsName);
		if (r == null)
			return null;
		return new ImmutableMemoryPeakResults(r);
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
		IJImageSource imageSource = new IJImageSource(imp);
		results.setSource(imageSource);

		// Get the calibration
		settings = SettingsManager.readCreateDataSettings(0).toBuilder();

		simulationParameters = showSimulationParametersDialog(imp, results);
		if (simulationParameters != null)
		{
			setBackground(results);
			setNoise(results, imp);
			setBenchmarkResults(imp, results);
			IJ.showStatus("Loaded " + TextUtils.pleural(results.size(), "result"));
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
		if (results.hasBackground())
			return;

		// Simple fix is to use the bias plus the global photon background.
		// TODO - Subtract the spots from the local region and compute the true local background.
		// Note this requires knowing the PSF width. If this is a loaded ground truth dataset then 
		// it probably will not have Gaussian widths.
		final float b = (float) (simulationParameters.bias + simulationParameters.gain * simulationParameters.b);
		results.setZeroBackground(b);
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
		if (results.hasNoise())
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
					NoiseEstimatorMethod.QUICK_RESIDUALS_LEAST_TRIMMED_OF_SQUARES);
		}

		Statistics stats = new Statistics(Arrays.copyOfRange(noise, 1, noise.length));
		System.out.printf("Noise = %.3f +/- %.3f (%d)\n", stats.getMean(), stats.getStandardDeviation(), stats.getN());

		results.forEach(new PeakResultProcedure()
		{
			public void execute(PeakResult p)
			{
				if (p.getFrame() < noise.length)
					p.setNoise(noise[p.getFrame()]);
			}
		});
	}

	private MemoryPeakResults getSimulationResults()
	{
		if (benchmarkAuto)
		{
			// Load directly from a results file. This is mainly to be used to load simulations 
			// saved to memory then saved to file.
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
		if (localisations == null || localisations.isEmpty())
			return null;

		return localisations.toPeakResults();
	}

	private SimulationParameters showSimulationParametersDialog(ImagePlus imp, MemoryPeakResults results)
	{
		int molecules = results.size();

		// Get the missing parameters from the user
		boolean fullSimulation = false;
		double s = -1;

		if (!results.convertToPreferredUnits())
		{
			IJ.error(TITLE,
					String.format("Results should be in the preferred units (%s,%s)",
							UnitHelper.getName(MemoryPeakResults.PREFERRED_DISTANCE_UNIT),
							UnitHelper.getName(MemoryPeakResults.PREFERRED_INTENSITY_UNIT)));
			return null;
		}

		// Get these from the data
		RawResultProcedure sp = new RawResultProcedure(results);
		sp.getBIXYZ();
		float[] signal = sp.intensity;
		float[] limits = Maths.limits(signal);
		double minSignal = limits[0];
		double maxSignal = limits[1];
		double signalPerFrame = Maths.sum(signal) / molecules;

		float[] depths = sp.z;
		limits = Maths.limits(depths);
		float depth = Math.max(Math.abs(limits[0]), Math.abs(limits[1]));
		boolean fixedDepth = Double.compare(limits[0], limits[1]) == 0;

		final CalibrationWriter cal = results.getCalibrationWriter();
		String iUnits = " " + UnitHelper.getName(cal.getIntensityUnit());
		String zUnits = " " + UnitHelper.getName(cal.getDistanceUnit());

		// Get this from the user
		double b = -1;

		// Use last simulation parameters for missing settings.
		// This is good if we are re-running the plugin to load data.
		Rectangle lastCameraBounds = null;
		if (simulationParameters != null && simulationParameters.isLoaded())
		{
			fullSimulation = simulationParameters.fullSimulation;
			s = simulationParameters.s;
			b = simulationParameters.b;
			if (!cal.hasBias())
				cal.setBias(simulationParameters.bias);
			if (!cal.hasCountPerPhoton())
				cal.setCountPerPhoton(simulationParameters.gain);
			if (!cal.hasQuantumEfficiency())
				cal.setQuantumEfficiency(simulationParameters.qe);
			if (!cal.hasReadNoise())
				cal.setReadNoise(simulationParameters.readNoise);
			if (!cal.hasCameraType())
				cal.setCameraType(simulationParameters.cameraType);
			if (!cal.hasNmPerPixel())
				cal.setNmPerPixel(simulationParameters.a);
			if (!cal.hasCameraModelName())
				cal.setCameraModelName(simulationParameters.cameraModelName);
			lastCameraBounds = simulationParameters.cameraBounds;
		}

		// Show a dialog to confirm settings
		final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		StringBuilder sb = new StringBuilder();
		sb.append("Results contain ").append(TextUtils.pleural(molecules, "molecule")).append('\n');
		sb.append("Min signal = ").append(Utils.rounded(minSignal)).append(iUnits).append('\n');
		sb.append("Max signal = ").append(Utils.rounded(maxSignal)).append(iUnits).append('\n');
		sb.append("Av signal = ").append(Utils.rounded(signalPerFrame)).append(iUnits).append('\n');
		if (fixedDepth)
			sb.append("Fixed depth = ").append(Utils.rounded(depth)).append(zUnits).append('\n');
		gd.addMessage(sb.toString());

		gd.addCheckbox("Flourophore_simulation", fullSimulation);
		gd.addNumericField("Gaussian_SD", s, 3, 8, "nm");
		gd.addNumericField("Pixel_pitch", cal.getNmPerPixel(), 3, 8, "nm");
		gd.addNumericField("Background", b, 3, 8, "photon");

		// Camera type does not need the full simulation settings. Plus the units are different
		// so just re-implement.
		gd.addChoice("Camera_type", SettingsManager.getCameraTypeNames(),
				CalibrationProtosHelper.getName(settings.getCameraType()), new OptionListener<Integer>()
				{
					public boolean collectOptions(Integer field)
					{
						settings.setCameraType(SettingsManager.getCameraTypeValues()[field]);
						boolean result = collectOptions(false);
						return result;
					}

					public boolean collectOptions()
					{
						return collectOptions(true);
					}

					private boolean collectOptions(boolean silent)
					{
						CameraType cameraType = settings.getCameraType();
						boolean isCCD = CalibrationProtosHelper.isCCDCameraType(cameraType);
						ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE, null);
						if (isCCD)
						{
							egd.addNumericField("Total_gain", cal.getCountPerPhoton(), 3, 8, "count/photon");
							egd.addNumericField("Quantum_efficiency", cal.getQuantumEfficiency(), 3, 8, "e-/photon");
							egd.addNumericField("Read_noise", cal.getReadNoise(), 3, 8, "count");
							egd.addNumericField("Bias", cal.getBias(), 3, 8, "count");
						}
						else if (cameraType == CameraType.SCMOS)
						{
							String[] models = CameraModelManager.listCameraModels(true);
							egd.addChoice("Camera_model_name", models, cal.getCameraModelName());
							egd.addNumericField("Quantum_efficiency", cal.getQuantumEfficiency(), 2, 6,
									"electron/photon");
						}
						else
						{
							IJ.error("Unsupported camera type " + CalibrationProtosHelper.getName(cameraType));
							return false;
						}
						egd.setSilent(silent);
						egd.showDialog(true, gd);
						if (egd.wasCanceled())
							return false;
						if (isCCD)
						{
							cal.setCountPerPhoton(egd.getNextNumber());
							cal.setQuantumEfficiency(egd.getNextNumber());
							cal.setReadNoise(egd.getNextNumber());
							cal.setBias(egd.getNextNumber());
						}
						else if (cameraType == CameraType.SCMOS)
						{
							cal.setCameraModelName(egd.getNextChoice());
							cal.setQuantumEfficiency(Math.abs(egd.getNextNumber()));
						}
						return true;
					}
				});

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
		settings.setCameraType(SettingsManager.getCameraTypeValues()[gd.getNextChoiceIndex()]);

		float myDepth = depth;
		if (!fixedDepth)
		{
			myDepth = (float) gd.getNextNumber();
			if (myDepth < depth)
			{
				IJ.error(TITLE, String.format("Input depth is smaller than the depth guessed from the data: %f < %f",
						myDepth, depth));
				return null;
			}
			depth = myDepth;
		}

		gd.collectOptions();

		// Validate settings
		Rectangle modelBounds = null;
		try
		{
			Parameters.isAboveZero("Gaussian SD", s);
			Parameters.isAboveZero("Pixel pitch", cal.getNmPerPixel());
			Parameters.isPositive("Background", b);

			Parameters.isAboveZero("Quantum efficiency", cal.getQuantumEfficiency());
			Parameters.isEqualOrBelow("Quantum efficiency", cal.getQuantumEfficiency(), 1);

			if (cal.isCCDCamera())
			{
				Parameters.isAboveZero("Total gain", cal.getCountPerPhoton());
				Parameters.isPositive("Read noise", cal.getReadNoise());
				Parameters.isPositive("Bias", cal.getBias());
			}
			else if (cal.isSCMOS())
			{
				// Load the model
				cameraModel = CameraModelManager.load(cal.getCameraModelName());
				if (cameraModel == null)
				{
					IJ.error(TITLE, "Unknown camera model for name: " + cal.getCameraModelName());
					return null;
				}

				int ox = 0, oy = 0;
				if (lastCameraBounds != null)
				{
					ox = lastCameraBounds.x;
					oy = lastCameraBounds.y;
				}
				cameraModel = PeakFit.cropCameraModel(cameraModel,
						new Rectangle(ox, oy, imp.getWidth(), imp.getHeight()), null, false);
				modelBounds = cameraModel.getBounds();

				IJImageSource imageSource = (IJImageSource) results.getSource();
				imageSource.setOrigin(modelBounds.x, modelBounds.y);

				cal.clearGlobalCameraSettings();
			}
			else
			{
				IJ.error(TITLE, "Unknown camera type: " + cal.getCameraType());
				return null;
			}
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return null;
		}

		// Store calibration
		results.setCalibration(cal.getCalibration());

		double a = cal.getNmPerPixel();
		double bias = cal.getBias();
		double gain = cal.getCountPerPhoton();
		double readNoise = cal.getReadNoise();
		double qe = cal.getQuantumEfficiency();

		// Note: The calibration will throw an exception if the converter cannot be created.
		// This is OK as the data will be invalid.

		// Convert values to photons
		TypeConverter<IntensityUnit> ic = cal.getIntensityConverter(IntensityUnit.PHOTON);
		minSignal = ic.convert(minSignal);
		maxSignal = ic.convert(maxSignal);
		signalPerFrame = ic.convert(signalPerFrame);

		// Convert +/- depth to total depth in nm
		depth = cal.getDistanceConverter(DistanceUnit.NM).convert(depth * 2);

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
				signalPerFrame, depth, fixedDepth, bias, gain, qe, readNoise, cal.getCameraType(),
				cal.getCameraModelName(), modelBounds, b, b2);
		p.loaded = true;
		return p;
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
		gd.addMessage(
				TextUtils.wrap("Specify if the results are preprocessed. This is true only if the simulation was " +
						"previously loaded and then saved to a GDSC SMLM file format from memory. Set to false " +
						"to load using a universal results loader.", 80));
		gd.addCheckbox("Preprocessed_results", benchmarkAuto);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		benchmarkImage = gd.getNextChoice();
		benchmarkFile = gd.getNextString();
		benchmarkAuto = gd.getNextBoolean();

		return true;
	}

	/**
	 * Gets the camera model for processing frames from the simulation image. The model will have bounds that match the
	 * simulation image dimensions.
	 *
	 * @param parameters
	 *            the parameters
	 * @return the camera model
	 * @throws ConfigurationException
	 *             If the model cannot be created
	 */
	public static CameraModel getCameraModel(BaseParameters parameters) throws ConfigurationException
	{
		// Create the camera model
		switch (parameters.cameraType)
		{
			case CCD:
			case EMCCD:
				return new FixedPixelCameraModel(parameters.bias, parameters.gain);

			case SCMOS:
				CameraModel cameraModel = CameraModelManager.load(parameters.cameraModelName);
				if (cameraModel == null)
					throw new ConfigurationException("Unknown camera model for name: " + parameters.cameraModelName);
				try
				{
					cameraModel = cameraModel.crop(parameters.cameraBounds, true);
				}
				catch (IllegalArgumentException e)
				{
					throw new ConfigurationException(e);
				}
				return cameraModel;

			case UNRECOGNIZED:
			case CAMERA_TYPE_NA:
			default:
				throw new ConfigurationException("Unknown camera model");
		}
	}

	public static StringBuilder addCameraDescription(StringBuilder sb, BaseParameters parameters)
			throws ConfigurationException
	{
		if (parameters.cameraType == CameraType.SCMOS)
		{
			sb.append("sCMOS (").append(parameters.cameraModelName).append(") ");
			Rectangle bounds = parameters.cameraBounds;
			sb.append(" ").append(bounds.x).append(",").append(bounds.y);
			sb.append(" ").append(bounds.width).append("x").append(bounds.height);
		}
		else if (CalibrationProtosHelper.isCCDCameraType(parameters.cameraType))
		{
			sb.append(CalibrationProtosHelper.getName(parameters.cameraType));
			sb.append(" G=").append(parameters.gain);
			sb.append(" RN=").append(parameters.readNoise);
			sb.append(" B=").append(parameters.bias);
		}
		else
		{
			throw new IllegalStateException();
		}
		return sb;
	}
}
