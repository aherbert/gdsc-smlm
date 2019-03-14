/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */

package uk.ac.sussex.gdsc.smlm.ij.plugins;

import uk.ac.sussex.gdsc.core.clustering.DensityManager;
import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.core.data.utils.ConversionException;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot.HistogramPlotBuilder;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog.OptionListener;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.threshold.AutoThreshold;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.FileUtils;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.MemoryUtils;
import uk.ac.sussex.gdsc.core.utils.RandomUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.TurboList;
import uk.ac.sussex.gdsc.core.utils.UnicodeReader;
import uk.ac.sussex.gdsc.core.utils.concurrent.ConcurrencyUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraType;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.ConfigurationException;
import uk.ac.sussex.gdsc.smlm.data.config.CreateDataSettingsHelper;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.NoiseEstimatorMethod;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.CreateDataSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.LoadLocalisationsSettings;
import uk.ac.sussex.gdsc.smlm.data.config.MoleculeProtos.Atom;
import uk.ac.sussex.gdsc.smlm.data.config.MoleculeProtos.AtomOrBuilder;
import uk.ac.sussex.gdsc.smlm.data.config.MoleculeProtos.Mixture;
import uk.ac.sussex.gdsc.smlm.data.config.MoleculeProtos.Molecule;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.AstigmatismModel;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.ImagePSF;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.Offset;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.data.config.PsfHelper;
import uk.ac.sussex.gdsc.smlm.data.config.PsfProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.engine.FitWorker;
import uk.ac.sussex.gdsc.smlm.filters.GaussianFilter;
import uk.ac.sussex.gdsc.smlm.fitting.FisherInformationMatrix;
import uk.ac.sussex.gdsc.smlm.fitting.UnivariateLikelihoodFisherInformationCalculator;
import uk.ac.sussex.gdsc.smlm.function.BasePoissonFisherInformation;
import uk.ac.sussex.gdsc.smlm.function.FunctionHelper;
import uk.ac.sussex.gdsc.smlm.function.InterpolatedPoissonFisherInformation;
import uk.ac.sussex.gdsc.smlm.function.PoissonGaussianFisherInformation;
import uk.ac.sussex.gdsc.smlm.function.StandardValueProcedure;
import uk.ac.sussex.gdsc.smlm.function.gaussian.AstigmatismZModel;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.HoltzerAstigmatismZModel;
import uk.ac.sussex.gdsc.smlm.ij.IJImageSource;
import uk.ac.sussex.gdsc.smlm.ij.plugins.LoadLocalisations.LocalisationList;
import uk.ac.sussex.gdsc.smlm.ij.settings.ImagePsfHelper;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.math3.distribution.CustomGammaDistribution;
import uk.ac.sussex.gdsc.smlm.math3.distribution.CustomPoissonDistribution;
import uk.ac.sussex.gdsc.smlm.model.ActivationEnergyImageModel;
import uk.ac.sussex.gdsc.smlm.model.AiryPattern;
import uk.ac.sussex.gdsc.smlm.model.AiryPsfModel;
import uk.ac.sussex.gdsc.smlm.model.CompoundMoleculeModel;
import uk.ac.sussex.gdsc.smlm.model.DiffusionType;
import uk.ac.sussex.gdsc.smlm.model.FixedLifetimeImageModel;
import uk.ac.sussex.gdsc.smlm.model.FluorophoreSequenceModel;
import uk.ac.sussex.gdsc.smlm.model.GaussianPsfModel;
import uk.ac.sussex.gdsc.smlm.model.GridDistribution;
import uk.ac.sussex.gdsc.smlm.model.ImageModel;
import uk.ac.sussex.gdsc.smlm.model.ImagePsfModel;
import uk.ac.sussex.gdsc.smlm.model.LocalisationModel;
import uk.ac.sussex.gdsc.smlm.model.LocalisationModelSet;
import uk.ac.sussex.gdsc.smlm.model.MaskDistribution;
import uk.ac.sussex.gdsc.smlm.model.MaskDistribution3D;
import uk.ac.sussex.gdsc.smlm.model.MoleculeModel;
import uk.ac.sussex.gdsc.smlm.model.PsfModel;
import uk.ac.sussex.gdsc.smlm.model.PsfModelGradient1Function;
import uk.ac.sussex.gdsc.smlm.model.RadialFalloffIllumination;
import uk.ac.sussex.gdsc.smlm.model.SpatialDistribution;
import uk.ac.sussex.gdsc.smlm.model.SpatialIllumination;
import uk.ac.sussex.gdsc.smlm.model.SphericalDistribution;
import uk.ac.sussex.gdsc.smlm.model.UniformDistribution;
import uk.ac.sussex.gdsc.smlm.model.UniformIllumination;
import uk.ac.sussex.gdsc.smlm.model.camera.CameraModel;
import uk.ac.sussex.gdsc.smlm.model.camera.CcdCameraModel;
import uk.ac.sussex.gdsc.smlm.model.camera.EmCcdCameraModel;
import uk.ac.sussex.gdsc.smlm.results.ExtendedPeakResult;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.IdPeakResult;
import uk.ac.sussex.gdsc.smlm.results.ImageSource.ReadHint;
import uk.ac.sussex.gdsc.smlm.results.ImmutableMemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.PeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResultsReader;
import uk.ac.sussex.gdsc.smlm.results.SynchronizedPeakResults;
import uk.ac.sussex.gdsc.smlm.results.TextFilePeakResults;
import uk.ac.sussex.gdsc.smlm.results.count.FrameCounter;
import uk.ac.sussex.gdsc.smlm.results.procedures.BirResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.PrecisionResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.RawResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.StandardResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.WidthResultProcedure;

import com.google.protobuf.TextFormat;

import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.set.hash.TIntHashSet;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import ij.text.TextWindow;

import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.distribution.RealDistribution;
import org.apache.commons.math3.distribution.UniformIntegerDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.random.EmpiricalDistribution;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.SobolSequenceGenerator;
import org.apache.commons.math3.random.Well44497b;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.rng.simple.RandomSource;

import java.awt.Checkbox;
import java.awt.Rectangle;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
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

/**
 * Creates data using a simulated PSF.
 */
public class CreateData implements PlugIn, ItemListener {
  /** The title. */
  static final String TITLE = "Create Data";
  private static final String CREATE_DATA_IMAGE_TITLE = "Localisation Data";

  private static final String[] ILLUMINATION = {"Uniform", "Radial"};
  private static final int RADIAL = 1;
  private static final String[] DISTRIBUTION =
      {"Uniform RNG", "Uniform Halton", "Uniform Sobol", "Mask", "Grid"};
  private static final int UNIFORM_HALTON = 1;
  private static final int UNIFORM_SOBOL = 2;
  private static final int MASK = 3;
  private static final int GRID = 4;
  private static final String[] CONFINEMENT = {"None", "Mask", "Sphere", "Within Image"};
  private static final int CONFINEMENT_MASK = 1;
  private static final int CONFINEMENT_SPHERE = 2;
  private static final int CONFINEMENT_WITHIN_IMAGE = 3;
  private static final String[] PHOTON_DISTRIBUTION =
      {"Uniform", "Gamma", "Custom", "Fixed", "Correlated"};
  private static final int PHOTON_UNIFORM = 0;
  private static final int PHOTON_GAMMA = 1;
  private static final int PHOTON_CUSTOM = 2;
  private static final int PHOTON_FIXED = 3;
  private static final int PHOTON_CORRELATED = 4;

  private static final String[] PSF_MODELS =
      new String[] {"2D Gaussian", "Airy", "Image", "Astigmatism"};
  private static final int PSF_MODEL_GAUSSIAN = 0;
  private static final int PSF_MODEL_AIRY = 1;
  private static final int PSF_MODEL_IMAGE = 2;
  private static final int PSF_MODEL_ASTIGMATISM = 3;

  /**
   * The PSF model type. This is set when validating the PSF settings.
   */
  private int psfModelType = -1;
  private AstigmatismModel astigmatismModel;
  private PsfModel psfModelCache;

  private static TextWindow summaryTable;
  private static AtomicInteger datasetNumber = new AtomicInteger();
  private static double areaInUm;
  private static String header;

  private CreateDataSettings.Builder settings;

  private static final String[] NAMES = new String[] {"Samples/Frame", "Signal/Frame",
      "Signal/Frame (continuous)", "Total Signal", "Blinks", "t-On", "t-Off", "Sampled blinks",
      "Sampled t-On", "Sampled t-Off", "Noise", "SNR", "SNR (continuous)", "Density", "Precision",
      "Precision (in-focus)", "X", "Y", "Z", "Width"};
  private static boolean[] displayHistograms = new boolean[NAMES.length];

  static {
    for (int i = 0; i < displayHistograms.length; i++) {
      displayHistograms[i] = true;
    }
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

  static {
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

  static {
    alwaysRemoveOutliers = new boolean[NAMES.length];
    alwaysRemoveOutliers[PRECISION] = true;
    alwaysRemoveOutliers[PRECISION_IN_FOCUS] = true;
  }

  private String resultsFileHeader;
  private AtomicInteger photonsRemoved;
  private AtomicInteger removedT1;
  private AtomicInteger removedTn;
  private SummaryStatistics photonStats;
  private double hwhm;
  private PSF psf;

  private TIntHashSet movingMolecules;
  private TIntIntHashMap idToCompound;
  private ArrayList<String> compoundNames;
  private boolean maskListContainsStacks;

  // Created by drawImage(...)
  private MemoryPeakResults results;

  // Used by the ImageGenerator to show progress when the thread starts
  private int frame;
  private int maxT;
  private int totalFrames;

  private boolean simpleMode;
  private boolean benchmarkMode;
  private boolean spotMode;
  private boolean trackMode;
  private boolean extraOptions;

  // Hold private variables for settings that are ignored in simple/benchmark mode
  private boolean poissonNoise = true;
  private double minPhotons;
  private double minSnrT1;
  private double minSnrTn;

  // Compute the CRLB for the PSF using the fisher information
  private BasePoissonFisherInformation[] fiFunction;

  /** Store the parameters. */
  public static class BaseParameters {
    private static int nextId = 1;

    /**
     * The parameter set identifier.
     */
    final int id;
    /**
     * Gaussian standard deviation.
     */
    final double sd;
    /**
     * Pixel pitch in nm.
     */
    final double pixelPitch;
    /**
     * The min number of photons per frame.
     */
    final double minSignal;
    /**
     * The max number of photons per frame.
     */
    final double maxSignal;
    /**
     * The average signal per frame.
     */
    double averageSignal;
    /**
     * The camera bias.
     */
    final double bias;
    /**
     * Total gain (ADUs/photon).
     */
    final double gain;
    /**
     * Quantum efficiency (electron/photon).
     */
    final double qe;
    /**
     * Read noise in ADUs.
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
     * Background.
     */
    double background;
    /**
     * Background noise in photons per pixel (used in the precision calculations).
     */
    final double noise;

    /** The best approximation of the PSF used to fit the simulation data. */
    final PSF psf;

    /**
     * Instantiates a new base parameters.
     *
     * @param sd the s
     * @param pixelPitch the pixel pitch
     * @param minSignal the min signal
     * @param maxSignal the max signal
     * @param averageSignal the average signal
     * @param bias the bias
     * @param gain the gain
     * @param qe the qe
     * @param readNoise the read noise
     * @param cameraType the camera type
     * @param cameraModelName the camera model name
     * @param cameraBounds the camera bounds
     * @param background the background
     * @param noise the noise
     * @param psf the psf
     */
    public BaseParameters(double sd, double pixelPitch, double minSignal, double maxSignal,
        double averageSignal, double bias, double gain, double qe, double readNoise,
        CameraType cameraType, String cameraModelName, Rectangle cameraBounds, double background,
        double noise, PSF psf) {
      id = nextId++;
      this.sd = sd;
      this.pixelPitch = pixelPitch;
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
      this.background = background;
      this.noise = noise;
      this.psf = psf;
    }

    /**
     * Checks if is emccd.
     *
     * @return true, if is emccd
     */
    public boolean isEmCcd() {
      return cameraType == CameraType.EMCCD;
    }
  }

  /** Store the parameters for the last simulation for spot data. */
  public static class SimulationParameters extends BaseParameters {
    /**
     * Number of molecules in the simulated image.
     */
    final int molecules;
    /**
     * True if using a full simulation of fluorophores with a lifetime. False is for single random
     * localisations per frame.
     */
    final boolean fullSimulation;
    /**
     * The z-position depth.
     */
    final double depth;
    /**
     * True if the depth is fixed.
     */
    final boolean fixedDepth;

    private boolean loaded;

    /**
     * Instantiates a new simulation parameters.
     *
     * @param molecules the molecules
     * @param fullSimulation the full simulation
     * @param sd the s
     * @param pixelPitch the pixel pitch
     * @param minSignal the min signal
     * @param maxSignal the max signal
     * @param averageSignal the average signal
     * @param depth the depth
     * @param fixedDepth the fixed depth
     * @param bias the bias
     * @param gain the gain
     * @param qe the qe
     * @param readNoise the read noise
     * @param cameraType the camera type
     * @param cameraModelName the camera model name
     * @param cameraBounds the camera bounds
     * @param background the background
     * @param noise the noise
     * @param psf the psf
     */
    public SimulationParameters(int molecules, boolean fullSimulation, double sd, double pixelPitch,
        double minSignal, double maxSignal, double averageSignal, double depth, boolean fixedDepth,
        double bias, double gain, double qe, double readNoise, CameraType cameraType,
        String cameraModelName, Rectangle cameraBounds, double background, double noise, PSF psf) {
      super(sd, pixelPitch, minSignal, maxSignal, averageSignal, bias, gain, qe, readNoise,
          cameraType, cameraModelName, cameraBounds, background, noise, psf);
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
    public boolean isLoaded() {
      return loaded;
    }
  }

  /** Store the parameters for the last benchmark. */
  public static class BenchmarkParameters extends BaseParameters {
    /**
     * Number of frames in the simulated image.
     */
    final int frames;
    /**
     * The x position of the localisation in each frame.
     */
    final double x;
    /**
     * The y position of the localisation in each frame.
     */
    final double y;
    /**
     * The z position of the localisation in each frame.
     */
    final double z;
    /**
     * The actual number of simulated photons in each frame of the benchmark image. Some frames may
     * be empty (due to signal filtering or Poisson sampling).
     */
    final double[] framePhotons;
    /**
     * The actual number of simulated background photons in each frame of the benchmark image.
     */
    final double[] frameBackground;

    /**
     * The Cramer-Rao lower bounds computed from the PSF model and the Fisher information.
     * [Background,Intensity,X,Y,Z]
     */
    final double[] crlb;

    /**
     * The number of frames with a simulated photon count.
     */
    private int molecules;

    /** The precision of the signal (N) in photons. */
    final double precisionN;

    /** The precision of the (x,y) position in nm assuming least-squares fitting. */
    final double precisionX;

    /** The precision of the (x,y) position in nm assuming Maximum Likelihood fitting. */
    final double precisionXml;

    /**
     * Instantiates a new benchmark parameters.
     *
     * @param frames the frames
     * @param sd the s
     * @param pixelPitch the pixel pitch
     * @param signal the signal
     * @param x the x
     * @param y the y
     * @param z the z
     * @param bias the bias
     * @param gain the gain
     * @param qe the quantum efficiency
     * @param readNoise the read noise
     * @param cameraType the camera type
     * @param cameraModelName the camera model name
     * @param cameraBounds the camera bounds
     * @param background the background
     * @param noise the noise
     * @param precisionN The precision of the signal (N) in photons.
     * @param precisionX The precision of the (x,y) position in nm assuming least-squares fitting.
     * @param precisionXml The precision of the (x,y) position in nm assuming Maximum Likelihood
     *        fitting.
     * @param psf the PSF
     * @param crlb the Cramér–Rao Lower Bounds
     */
    public BenchmarkParameters(int frames, double sd, double pixelPitch, double signal, double x,
        double y, double z, double bias, double gain, double qe, double readNoise,
        CameraType cameraType, String cameraModelName, Rectangle cameraBounds, double background,
        double noise, double precisionN, double precisionX, double precisionXml, PSF psf,
        double[] crlb) {
      super(sd, pixelPitch, signal, signal, signal, bias, gain, qe, readNoise, cameraType,
          cameraModelName, cameraBounds, background, noise, psf);

      this.frames = frames;
      this.x = x;
      this.y = y;
      this.z = z;
      this.precisionN = precisionN;
      this.precisionX = precisionX;
      this.precisionXml = precisionXml;
      framePhotons = new double[frames];
      frameBackground = new double[frames];
      this.crlb = crlb;
    }

    /**
     * Sets the photons.
     *
     * @param results the new photons
     */
    public void setPhotons(MemoryPeakResults results) {
      molecules = results.size();
      // Get the signal and background in photons
      results.forEach(IntensityUnit.PHOTON, (BirResultProcedure) (bg, intensity, result) -> {
        final int i = result.getFrame() - 1;
        if (framePhotons[i] != 0) {
          throw new IllegalArgumentException(
              "Multiple peaks on the same frame: " + result.getFrame());
        }
        framePhotons[i] = intensity;
        frameBackground[i] = bg;
      });
      final double av = MathUtils.sum(framePhotons) / molecules;
      final double av2 = MathUtils.sum(frameBackground) / molecules;
      ImageJUtils.log(
          "Created %d frames, %d molecules. Simulated signal %s : average %s. "
              + "Simulated background %s : average %s",
          frames, molecules, MathUtils.rounded(averageSignal), MathUtils.rounded(av),
          MathUtils.rounded(background), MathUtils.rounded(av2));
      // Reset the average signal and background (in photons)
      averageSignal = av;
      background = av2;
    }

    /**
     * Gets the average number of photons per frame (signal).
     *
     * @return The average number of photons per frame.
     */
    public double getSignal() {
      return averageSignal;
    }

    /**
     * Gets the number of molecules. There should be 0 or 1 molecule per frame.
     *
     * @return The number of molecules
     */
    public int getMolecules() {
      return molecules;
    }

    /**
     * Gets the average number of background photons per frame.
     *
     * @return the average number of background photons per frame.
     */
    public double getBackground() {
      return background;
    }
  }

  /** The last benchmark parameters. */
  static BenchmarkParameters benchmarkParameters;

  /** The last simulation parameters. */
  static SimulationParameters simulationParameters;

  private static String benchmarkFile = "";
  private static LoadLocalisationsSettings.Builder loadSettings;
  private static String benchmarkImage = "";
  private static boolean benchmarkAuto;
  private static int benchmarkImageId;
  private static String benchmarkResultsName = "";

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    extraOptions = ImageJUtils.isExtraOptions();
    simpleMode = (arg != null && arg.contains("simple"));
    benchmarkMode = (arg != null && arg.contains("benchmark"));
    spotMode = (arg != null && arg.contains("spot"));
    trackMode = (arg != null && arg.contains("track"));

    if ("load".equals(arg)) {
      loadBenchmarkData();
      return;
    }

    // Each localisation set is a collection of localisations that represent all localisations
    // with the same ID that are on in the same image time frame (Note: the simulation
    // can create many localisations per fluorophore per time frame which is useful when
    // modelling moving particles)
    List<LocalisationModelSet> localisationSets = null;

    // Each fluorophore contains the on and off times when light was emitted
    List<? extends FluorophoreSequenceModel> fluorophores = null;

    if (simpleMode || benchmarkMode || spotMode) {
      if (!showSimpleDialog()) {
        return;
      }
      resetMemory();

      settings.setExposureTime(1000); // 1 second frames
      areaInUm = settings.getSize() * settings.getPixelPitch() * settings.getSize()
          * settings.getPixelPitch() / 1e6;

      // Number of spots per frame
      int count = 0;
      int[] nextN = null;
      SpatialDistribution dist;

      if (benchmarkMode) {
        // --------------------
        // BENCHMARK SIMULATION
        // --------------------
        // Draw the same point on the image repeatedly
        count = 1;
        dist = createFixedDistribution();
        try {
          reportAndSaveFittingLimits(dist);
        } catch (final Exception ex) {
          // This will be from the computation of the CRLB
          IJ.error(TITLE, ex.getMessage());
          return;
        }
      } else if (spotMode) {
        // ---------------
        // SPOT SIMULATION
        // ---------------
        // The spot simulation draws 0 or 1 random point per frame.
        // Ensure we have 50% of the frames with a spot.
        nextN = new int[settings.getParticles() * 2];
        Arrays.fill(nextN, 0, settings.getParticles(), 1);
        RandomUtils.shuffle(nextN, RandomSource.create(RandomSource.SPLIT_MIX_64));

        // Only put spots in the central part of the image
        final double border = settings.getSize() / 4.0;
        dist = createUniformDistribution(border);
      } else {
        // -----------------
        // SIMPLE SIMULATION
        // -----------------
        // The simple simulation draws n random points per frame to achieve a specified density.
        // No points will appear in multiple frames.
        // Each point has a random number of photons sampled from a range.

        // We can optionally use a mask. Create his first as it updates the areaInUm
        dist = createDistribution();

        // Randomly sample (i.e. not uniform density in all frames)
        if (settings.getSamplePerFrame()) {
          final double mean = areaInUm * settings.getDensity();
          ImageJUtils.log("Mean samples = %f", mean);
          if (mean < 0.5) {
            final GenericDialog gd = new GenericDialog(TITLE);
            gd.addMessage(
                "The mean samples per frame is low: " + MathUtils.rounded(mean) + "\n \nContinue?");
            gd.enableYesNoCancel();
            gd.hideCancelButton();
            gd.showDialog();
            if (!gd.wasOKed()) {
              return;
            }
          }
          final PoissonDistribution poisson = new PoissonDistribution(createRandomGenerator(), mean,
              PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);
          final StoredDataStatistics samples = new StoredDataStatistics(settings.getParticles());
          while (samples.getSum() < settings.getParticles()) {
            samples.add(poisson.sample());
          }
          nextN = new int[samples.getN()];
          for (int i = 0; i < nextN.length; i++) {
            nextN[i] = (int) samples.getValue(i);
          }
        } else {
          // Use the density to get the number per frame
          count = (int) FastMath.max(1, Math.round(areaInUm * settings.getDensity()));
        }
      }

      RandomGenerator random = null;

      // List<LocalisationModel> localisations = new ArrayList<>(settings.getParticles());
      localisationSets = new ArrayList<>(settings.getParticles());

      final int minPhotons = (int) settings.getPhotonsPerSecond();
      final int range = (int) settings.getPhotonsPerSecondMaximum() - minPhotons + 1;
      if (range > 1) {
        random = createRandomGenerator();
      }

      // Add frames at the specified density until the number of particles has been reached
      int id = 0;
      int time = 0;
      while (id < settings.getParticles()) {
        // Allow the number per frame to be specified
        if (nextN != null) {
          if (time >= nextN.length) {
            break;
          }
          count = nextN[time];
        }

        // Simulate random positions in the frame for the specified density
        time++;
        for (int j = 0; j < count; j++) {
          final double[] xyz = dist.next();

          // Ignore within border. We do not want to draw things we cannot fit.
          // if (!distBorder.isWithinXy(xyz))
          // continue;

          // Simulate random photons
          final int intensity = minPhotons + ((random != null) ? random.nextInt(range) : 0);

          final LocalisationModel m =
              new LocalisationModel(id, time, xyz, intensity, LocalisationModel.CONTINUOUS);

          // Each localisation can be a separate localisation set
          final LocalisationModelSet set = new LocalisationModelSet(id, time);
          set.add(m);
          localisationSets.add(set);

          id++;
        }
      }
    } else {
      // This is used for the track mode as well as the full simulation.

      if (!showDialog()) {
        return;
      }
      resetMemory();

      areaInUm = settings.getSize() * settings.getPixelPitch() * settings.getSize()
          * settings.getPixelPitch() / 1e6;

      int totalSteps;
      double correlation = 0;
      ImageModel imageModel;

      if (trackMode) {
        // ----------------
        // TRACK SIMULATION
        // ----------------
        // In track mode we create fixed lifetime fluorophores that do not overlap in time.
        // This is the simplest simulation to test moving molecules.
        settings.setSeconds((int) Math.ceil(
            settings.getParticles() * (settings.getExposureTime() + settings.getTOn()) / 1000));
        totalSteps = 0;

        final double simulationStepsPerFrame =
            (settings.getStepsPerSecond() * settings.getExposureTime()) / 1000.0;
        imageModel = new FixedLifetimeImageModel(
            settings.getStepsPerSecond() * settings.getTOn() / 1000.0, simulationStepsPerFrame);
      } else {
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

        final SpatialIllumination activationIllumination =
            createIllumination(settings.getPulseRatio(), settings.getPulseInterval());

        // Generate additional frames so that each frame has the set number of simulation steps
        totalSteps = (int) Math.ceil(settings.getSeconds() * settings.getStepsPerSecond());

        // Since we have an exponential decay of activations
        // ensure half of the particles have activated by 30% of the frames.
        final double eAct = totalSteps * 0.3 * activationIllumination.getAveragePhotons();

        // Q. Does tOn/tOff change depending on the illumination strength?
        imageModel = new ActivationEnergyImageModel(eAct, activationIllumination,
            settings.getStepsPerSecond() * settings.getTOn() / 1000.0,
            settings.getStepsPerSecond() * settings.getTOffShort() / 1000.0,
            settings.getStepsPerSecond() * settings.getTOffLong() / 1000.0,
            settings.getNBlinksShort(), settings.getNBlinksLong());
        imageModel.setUseGeometricDistribution(settings.getNBlinksGeometricDistribution());

        // Only use the correlation if selected for the distribution
        if (PHOTON_DISTRIBUTION[PHOTON_CORRELATED].equals(settings.getPhotonDistribution())) {
          correlation = settings.getCorrelation();
        }
      }

      imageModel.setRandomGenerator(createRandomGenerator());
      imageModel.setPhotonBudgetPerFrame(true);
      imageModel.setDiffusion2D(settings.getDiffuse2D());
      imageModel.setRotation2D(settings.getRotate2D());

      IJ.showStatus("Creating molecules ...");
      final SpatialDistribution distribution = createDistribution();
      final List<CompoundMoleculeModel> compounds = createCompoundMolecules();
      if (compounds == null) {
        return;
      }
      final List<CompoundMoleculeModel> molecules = imageModel.createMolecules(compounds,
          settings.getParticles(), distribution, settings.getRotateInitialOrientation());

      // Activate fluorophores
      IJ.showStatus("Creating fluorophores ...");

      // Note: molecules list will be converted to compounds containing fluorophores
      fluorophores = imageModel.createFluorophores(molecules, totalSteps);

      if (fluorophores.isEmpty()) {
        IJ.error(TITLE, "No fluorophores created");
        return;
      }

      // Map the fluorophore ID to the compound for mixtures
      if (compounds.size() > 1) {
        idToCompound = new TIntIntHashMap(fluorophores.size());
        for (final FluorophoreSequenceModel l : fluorophores) {
          idToCompound.put(l.getId(), l.getLabel());
        }
      }

      IJ.showStatus("Creating localisations ...");

      // TODO - Output a molecule Id for each fluorophore if using compound molecules. This allows
      // analysis
      // of the ratio of trimers, dimers, monomers, etc that could be detected.

      totalSteps = checkTotalSteps(totalSteps, fluorophores);
      if (totalSteps == 0) {
        return;
      }

      imageModel.setPhotonDistribution(createPhotonDistribution());
      try {
        imageModel.setConfinementDistribution(createConfinementDistribution());
      } catch (final ConfigurationException ex) {
        // We asked the user if it was OK to continue and they said no
        return;
      }
      // This should be optimised
      imageModel.setConfinementAttempts(10);

      final List<LocalisationModel> localisations =
          imageModel.createImage(molecules, settings.getFixedFraction(), totalSteps,
              settings.getPhotonsPerSecond() / settings.getStepsPerSecond(), correlation,
              settings.getRotateDuringSimulation());

      // Re-adjust the fluorophores to the correct time
      if (settings.getStepsPerSecond() != 1) {
        final double scale = 1.0 / settings.getStepsPerSecond();
        for (final FluorophoreSequenceModel f : fluorophores) {
          f.adjustTime(scale);
        }
      }

      // Integrate the frames
      localisationSets = combineSimulationSteps(localisations);

      localisationSets = filterToImageBounds(localisationSets);
    }

    datasetNumber.getAndIncrement();

    final List<LocalisationModel> localisations = drawImage(localisationSets);

    if (localisations == null || localisations.isEmpty()) {
      IJ.error(TITLE, "No localisations created");
      return;
    }

    fluorophores = removeFilteredFluorophores(fluorophores, localisations);

    final double signalPerFrame = showSummary(fluorophores, localisations);

    if (!benchmarkMode) {
      final boolean fullSimulation = (!(simpleMode || spotMode));
      saveSimulationParameters(localisations.size(), fullSimulation, signalPerFrame);
    }

    IJ.showStatus("Saving data ...");

    saveFluorophores(fluorophores);
    saveImageResults(results);
    saveLocalisations(localisations);

    // The settings for the filenames may have changed
    SettingsManager.writeSettings(settings.build());

    IJ.showStatus("Done");
  }

  private static void resetMemory() {
    benchmarkParameters = null;
    simulationParameters = null;
    setBenchmarkResults(null, null);
    // Run the garbage collector to free memory
    MemoryUtils.runGarbageCollector();
  }

  /**
   * Output the theoretical limits for fitting a Gaussian and store the benchmark settings.
   *
   * @param dist The distribution
   */
  private void reportAndSaveFittingLimits(SpatialDistribution dist) {
    ImageJUtils.log(TITLE + " Benchmark");

    final double a = settings.getPixelPitch();
    final double[] xyz = dist.next().clone();
    final int size = settings.getSize();
    final double offset = size * 0.5;
    for (int i = 0; i < 2; i++) {
      xyz[i] += offset;
    }

    // Get the width for the z-depth by using the PSF Model
    final PsfModel psf = createPsfModel(xyz);
    psfModelCache = psf;

    double sd0;
    double sd1;
    if (psf instanceof GaussianPsfModel) {
      sd0 = ((GaussianPsfModel) psf).getS0(xyz[2]);
      sd1 = ((GaussianPsfModel) psf).getS1(xyz[2]);
    } else if (psf instanceof AiryPsfModel) {
      psf.create3D((double[]) null, size, size, 1, xyz[0], xyz[1], xyz[2], false);
      sd0 = ((AiryPsfModel) psf).getW0() * AiryPattern.FACTOR;
      sd1 = ((AiryPsfModel) psf).getW1() * AiryPattern.FACTOR;
    } else if (psf instanceof ImagePsfModel) {
      psf.create3D((double[]) null, size, size, 1, xyz[0], xyz[1], xyz[2], false);
      sd0 = ((ImagePsfModel) psf).getHwhm0() / Gaussian2DFunction.SD_TO_HWHM_FACTOR;
      sd1 = ((ImagePsfModel) psf).getHwhm1() / Gaussian2DFunction.SD_TO_HWHM_FACTOR;
    } else {
      throw new IllegalStateException("Unknown PSF: " + psf.getClass().getSimpleName());
    }

    final double sd = Gaussian2DPeakResultHelper.getStandardDeviation(sd0, sd1) * a;

    ImageJUtils.log("X = %s nm : %s px", MathUtils.rounded(xyz[0] * a),
        MathUtils.rounded(xyz[0], 6));
    ImageJUtils.log("Y = %s nm : %s px", MathUtils.rounded(xyz[1] * a),
        MathUtils.rounded(xyz[1], 6));
    ImageJUtils.log("Z = %s nm : %s px", MathUtils.rounded(xyz[2] * a),
        MathUtils.rounded(xyz[2], 6));
    ImageJUtils.log("Width (s) = %s nm : %s px", MathUtils.rounded(sd), MathUtils.rounded(sd / a));
    final double sa = PsfCalculator.squarePixelAdjustment(sd, a);
    ImageJUtils.log("Adjusted Width (sa) = %s nm : %s px", MathUtils.rounded(sa),
        MathUtils.rounded(sa / a));
    ImageJUtils.log("Signal (N) = %s - %s photons",
        MathUtils.rounded(settings.getPhotonsPerSecond()),
        MathUtils.rounded(settings.getPhotonsPerSecondMaximum()));

    boolean emCcd;
    double totalGain;
    final double qe = getQuantumEfficiency();
    final double noise = getNoiseEstimateInPhotoelectrons(qe);
    double readNoise;

    if (CalibrationProtosHelper.isCcdCameraType(settings.getCameraType())) {
      final CreateDataSettingsHelper helper = new CreateDataSettingsHelper(settings);
      emCcd = helper.isEmCcd;
      totalGain = helper.getTotalGainSafe();
      // Store read noise in ADUs
      readNoise =
          settings.getReadNoise() * ((settings.getCameraGain() > 0) ? settings.getCameraGain() : 1);
    } else if (settings.getCameraType() == CameraType.SCMOS) {
      // Assume sCMOS amplification is like a CCD for the precision computation.
      emCcd = false;
      // Not required for sCMOS
      totalGain = 0;
      readNoise = 0;
    } else {
      throw new IllegalStateException("Unknown camera type: " + settings.getCameraType());
    }

    // The precision calculation is dependent on the model. The classic Mortensen formula
    // is for a Gaussian Mask Estimator. Use other equation for MLE. The formula provided
    // for WLSE requires an offset to the background used to stabilise the fitting. This is
    // not implemented (i.e. we used an offset of zero) and in this case the WLSE precision
    // is the same as MLE with the caveat of numerical instability.

    // Apply QE directly to simulated photons to allow computation of bounds
    // with the captured photons...
    final double min = settings.getPhotonsPerSecond() * qe;
    final double max = settings.getPhotonsPerSecondMaximum() * qe;
    final double lowerP = Gaussian2DPeakResultHelper.getPrecision(a, sd, max, noise, emCcd);
    final double upperP = Gaussian2DPeakResultHelper.getPrecision(a, sd, min, noise, emCcd);
    final double lowerMlP = Gaussian2DPeakResultHelper.getMLPrecision(a, sd, max, noise, emCcd);
    final double upperMlP = Gaussian2DPeakResultHelper.getMLPrecision(a, sd, min, noise, emCcd);
    final double lowerN = getPrecisionN(a, sd, min, MathUtils.pow2(noise), emCcd);
    final double upperN = getPrecisionN(a, sd, max, MathUtils.pow2(noise), emCcd);

    if (settings.getCameraType() == CameraType.SCMOS) {
      ImageJUtils.log("sCMOS camera background estimate uses an average read noise");
    }
    ImageJUtils.log("Effective background noise = %s photo-electron "
        + "[includes read noise and background photons]", MathUtils.rounded(noise));
    ImageJUtils.log("Localisation precision (LSE): %s - %s nm : %s - %s px",
        MathUtils.rounded(lowerP), MathUtils.rounded(upperP), MathUtils.rounded(lowerP / a),
        MathUtils.rounded(upperP / a));
    ImageJUtils.log("Localisation precision (MLE): %s - %s nm : %s - %s px",
        MathUtils.rounded(lowerMlP), MathUtils.rounded(upperMlP), MathUtils.rounded(lowerMlP / a),
        MathUtils.rounded(upperMlP / a));
    ImageJUtils.log("Signal precision: %s - %s photo-electrons", MathUtils.rounded(lowerN),
        MathUtils.rounded(upperN));

    // Wrap to a function
    final PsfModelGradient1Function f = new PsfModelGradient1Function(psf, size, size);

    // Set parameters
    final double[] params = new double[5];
    // No background when computing the SNR
    // params[0] = settings.getBackground() * qe;
    params[1] = min;
    System.arraycopy(xyz, 0, params, 2, 3);

    // Compute SNR using mean signal at 50%. Assume the region covers the entire PSF.
    final double[] v = new StandardValueProcedure().getValues(f, params);
    final double u = FunctionHelper.getMeanValue(v, 0.5);
    final double u0 = MathUtils.max(v);

    // Store the benchmark settings when not using variable photons
    if (min == max) {
      ImageJUtils.log("50%% PSF SNR : %s : Peak SNR : %s", MathUtils.rounded(u / noise),
          MathUtils.rounded(u0 / noise));

      // Compute the true CRLB using the fisher information
      createLikelihoodFunction();

      // Compute Fisher information
      final UnivariateLikelihoodFisherInformationCalculator c =
          new UnivariateLikelihoodFisherInformationCalculator(f, fiFunction);

      // Get limits: Include background in the params
      params[0] = settings.getBackground() * qe;
      final FisherInformationMatrix m = c.compute(params);

      // Report and store the limits
      final double[] crlb = m.crlbSqrt();
      if (crlb != null) {
        ImageJUtils.log("Localisation precision (CRLB): B=%s, I=%s photons",
            MathUtils.rounded(crlb[0]), MathUtils.rounded(crlb[1]));
        ImageJUtils.log("Localisation precision (CRLB): X=%s, Y=%s, Z=%s nm : %s,%s,%s px",
            MathUtils.rounded(crlb[2] * a), MathUtils.rounded(crlb[3] * a),
            MathUtils.rounded(crlb[4] * a), MathUtils.rounded(crlb[2]), MathUtils.rounded(crlb[3]),
            MathUtils.rounded(crlb[4]));
      }

      benchmarkParameters =
          new BenchmarkParameters(settings.getParticles(), sd, a, settings.getPhotonsPerSecond(),
              xyz[0], xyz[1], xyz[2], settings.getBias(), totalGain, qe, readNoise,
              settings.getCameraType(), settings.getCameraModelName(), cameraModel.getBounds(),
              settings.getBackground(), noise, lowerN, lowerP, lowerMlP, createPsf(sd / a), crlb);
    } else {
      // SNR will just scale
      final double scale = max / min;
      ImageJUtils.log("50%% PSF SNR : %s - %s : Peak SNR : %s - %s", MathUtils.rounded(u / noise),
          MathUtils.rounded(scale * u / noise), MathUtils.rounded(u0 / noise),
          MathUtils.rounded(scale * u0 / noise));
      ImageJUtils.log(
          "Warning: Benchmark settings are only stored in memory when the number of photons is "
              + "fixed. Min %s != Max %s",
          MathUtils.rounded(settings.getPhotonsPerSecond()),
          MathUtils.rounded(settings.getPhotonsPerSecondMaximum()));
    }
  }

  private double getNoiseEstimateInPhotoelectrons(double qe) {
    if (CalibrationProtosHelper.isCcdCameraType(settings.getCameraType())) {
      // Background is in photons. Convert to electrons
      double backgroundVariance = settings.getBackground() * qe;

      // Read noise is in electrons.
      double readNoise = settings.getReadNoise();

      // In an EM-CCD camera the read noise (in electrons) is swamped by amplification of the
      // signal.
      // We get the same result by dividing the read noise (in electrons) by the EM-gain.
      if (settings.getCameraType() == CameraType.EMCCD && settings.getEmGain() > 1) {
        // Add EM-CCD noise factor. The implementation of the Mortensen formula
        // using the standard deviation expects this factor.
        // See uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper.getPrecision(double,
        // double, double, double, boolean)
        backgroundVariance *= 2;
        readNoise /= settings.getEmGain();
      }

      // Get the expected value at each pixel in electrons. Assuming a Poisson distribution this
      // is equal to the total variance at the pixel.
      return Math.sqrt(backgroundVariance + MathUtils.pow2(readNoise));
    } else if (settings.getCameraType() == CameraType.SCMOS) {
      // Assume sCMOS amplification is like a CCD. We need an average read noise to get an
      // approximation of the background noise for the precision computation.

      // We get the total background in electrons
      final double backgroundVariance = settings.getBackground() * qe;

      // Create the camera noise model
      createPerPixelCameraModelData(cameraModel);

      // Get the average read noise in electrons
      final double readNoise = (MathUtils.sum(this.readNoise) / this.readNoise.length);

      // Combine as above
      return Math.sqrt(backgroundVariance + MathUtils.pow2(readNoise));
    }
    throw new IllegalStateException("Unknown camera type: " + settings.getCameraType());
  }

  /**
   * Store the simulation settings.
   *
   * @param particles the particles
   * @param fullSimulation the full simulation
   * @param signalPerFrame the signal per frame
   */
  private void saveSimulationParameters(int particles, boolean fullSimulation,
      double signalPerFrame) {
    double totalGain;
    final double noise = getNoiseEstimateInPhotoelectrons(getQuantumEfficiency());
    double readNoise;

    if (CalibrationProtosHelper.isCcdCameraType(settings.getCameraType())) {
      final CreateDataSettingsHelper helper = new CreateDataSettingsHelper(settings);
      totalGain = helper.getTotalGainSafe();
      // Store read noise in ADUs
      readNoise = helper.getReadNoiseInCounts();
    } else if (settings.getCameraType() == CameraType.SCMOS) {
      // Not required for sCMOS
      totalGain = 0;
      readNoise = 0;
    } else {
      throw new IllegalStateException("Unknown camera type: " + settings.getCameraType());
    }

    final double sd = getPsfSd() * settings.getPixelPitch();

    final double qe = getQuantumEfficiency();

    simulationParameters = new SimulationParameters(particles, fullSimulation, sd,
        settings.getPixelPitch(), settings.getPhotonsPerSecond(),
        settings.getPhotonsPerSecondMaximum(), signalPerFrame, settings.getDepth(),
        settings.getFixedDepth(), settings.getBias(), totalGain, qe, readNoise,
        settings.getCameraType(), settings.getCameraModelName(), cameraModel.getBounds(),
        settings.getBackground(), noise, createPsf(sd / settings.getPixelPitch()));
  }

  /**
   * Calculate the signal precision for least squares fitting. Uses the Thompson formula: (Thompson,
   * et al (2002) Biophysical Journal 82, 2775-2783), equation 19
   *
   * @param pixelPitch The size of the pixels in nm
   * @param sd The peak standard deviation in nm
   * @param photons The peak signal in photons
   * @param b2 The expected number of photons per pixel from a background with spatially constant
   *        expectation value across the image (Note that this is b^2 not b, which could be the
   *        standard deviation of the image pixels)
   * @param emCcd True if an emCCD camera
   * @return The signal precision in photons
   */
  public static double getPrecisionN(double pixelPitch, double sd, double photons, double b2,
      boolean emCcd) {
    // EM-CCD noise factor
    final double F = (emCcd) ? 2 : 1;
    final double a2 = pixelPitch * pixelPitch;
    // 4 * pi = 12.56637061

    // Adjustment for square pixels
    // final double sa2 = s * s + a2 / 12.0;

    // Original Thompson formula modified for EM-gain noise factor.

    // TODO - Investigate if this limit is correct

    // My fitters approach this limit when background is 0 photon and EM-gain = 0.
    // The fitters are above this limit when background is >0 photon and EM-gain = 0.

    // The MLE fitter can approach this limit when background is 0 photon and EM-gain = 25.

    return Math.sqrt(F * (photons + (12.56637061 * sd * sd * b2) / a2));
    // return Math.sqrt(F * (N + (12.56637061 * sa2 * b2) / a2));
  }

  /**
   * Check if the total steps can fit all the fluorophores end times. If not then ask the user if
   * they want to draw extra frames. Return the total steps to simulate (either the original steps
   * or a larger number to fit all the data).
   *
   * @param totalSteps the total steps
   * @param fluorophores the fluorophores
   * @return The new total steps to simulate
   */
  private int checkTotalSteps(int totalSteps,
      List<? extends FluorophoreSequenceModel> fluorophores) {
    int max = totalSteps;
    for (final FluorophoreSequenceModel f : fluorophores) {
      if (max < f.getEndTime()) {
        max = (int) (f.getEndTime() + 1);
      }
    }
    if (max > totalSteps) {
      final GenericDialog gd = new GenericDialog(TITLE);
      gd.enableYesNoCancel();
      gd.hideCancelButton();
      final double simulationStepsPerFrame =
          (settings.getStepsPerSecond() * settings.getExposureTime()) / 1000.0;
      final int newFrames = 1 + (int) (max / simulationStepsPerFrame);

      if (totalSteps != 0) {
        final int totalFrames =
            (int) Math.ceil(settings.getSeconds() * 1000 / settings.getExposureTime());
        gd.addMessage(String.format(
            "Require %d (%s%%) additional frames to draw all fluorophores.\n"
                + "Do you want to add extra frames?",
            newFrames - totalFrames,
            MathUtils.rounded((100.0 * (newFrames - totalFrames)) / totalFrames, 3)));
      } else {
        gd.addMessage(String.format(
            "Require %d frames to draw all fluorophores.\nDo you want to proceed?", newFrames));
      }
      gd.showDialog();
      if (gd.wasOKed()) {
        totalSteps = max;
      }
    }
    return totalSteps;
  }

  private SpatialDistribution createDistribution() {
    if (settings.getDistribution().equals(DISTRIBUTION[MASK])) {
      final ImagePlus imp = WindowManager.getImage(settings.getDistributionMask());
      if (imp != null) {
        return createMaskDistribution(imp, settings.getDistributionMaskSliceDepth(), true);
      }
    } else if (settings.getDistribution().equals(DISTRIBUTION[GRID])) {
      return new GridDistribution(settings.getSize(),
          settings.getDepth() / settings.getPixelPitch(), settings.getCellSize(),
          settings.getProbabilityBinary(),
          settings.getMinBinaryDistance() / settings.getPixelPitch(),
          settings.getMaxBinaryDistance() / settings.getPixelPitch());
    }

    return createUniformDistributionWithPsfWidthBorder();
  }

  private SpatialDistribution createMaskDistribution(ImagePlus imp, double sliceDepth,
      boolean updateArea) {
    // Calculate the scale of the mask
    final int w = imp.getWidth();
    final int h = imp.getHeight();
    final double scaleX = (double) settings.getSize() / w;
    final double scaleY = (double) settings.getSize() / h;

    // Use an image for the distribution
    if (imp.getStackSize() > 1) {
      final ImageStack stack = imp.getImageStack();
      final List<int[]> masks = new ArrayList<>(stack.getSize());
      final int[] maxMask = new int[w * h];
      for (int slice = 1; slice <= stack.getSize(); slice++) {
        final int[] mask = extractMask(stack.getProcessor(slice));
        if (updateArea) {
          for (int i = 0; i < mask.length; i++) {
            if (mask[i] != 0) {
              maxMask[i] = 1;
            }
          }
        }
        masks.add(mask);
      }
      if (updateArea) {
        updateArea(maxMask, w, h);
      }

      if (sliceDepth == 0) {
        // Auto configure to the full depth of the simulation
        sliceDepth = settings.getDepth() / masks.size();
      }

      return new MaskDistribution3D(masks, w, h, sliceDepth / settings.getPixelPitch(), scaleX,
          scaleY, createRandomGenerator());
    }
    final int[] mask = extractMask(imp.getProcessor());
    if (updateArea) {
      updateArea(mask, w, h);
    }
    return new MaskDistribution(mask, w, h, settings.getDepth() / settings.getPixelPitch(), scaleX,
        scaleY, createRandomGenerator());
  }

  private void updateArea(int[] mask, int width, int height) {
    // Note: The apparent area will be bigger due to the PSF width blurring the edges.
    // Assume the entire max intensity mask is in focus and the PSF is Gaussian in the focal plane.
    // Blur the mask
    final float[] pixels = new float[mask.length];
    for (int i = 0; i < pixels.length; i++) {
      pixels[i] = mask[i];
    }

    final GaussianFilter blur = new GaussianFilter();
    final double scaleX = (double) settings.getSize() / width;
    final double scaleY = (double) settings.getSize() / height;
    final double extra = 1; // Allow extra?
    final double sd = getPsfSd() * extra;
    blur.convolve(pixels, width, height, sd / scaleX, sd / scaleY);

    // Count pixels in blurred mask. Ignore those that are very faint (at the edge of the region)
    int count = 0;

    // // By fraction of max value
    // float limit = 0.1f;
    // //float min = (float) (Maths.max(pixels) * limit);
    // float min = limit;
    // for (float f : pixels)
    // if (f > min)
    // c++;

    // // Rank in order and get fraction of sum
    // Arrays.sort(pixels);
    // double sum = 0;
    // double stop = Maths.sum(pixels) * 0.95;
    // float after = Maths.max(pixels) + 1;
    // for (float f : pixels)
    // {
    // sum += f;
    // if (sum > stop)
    // {
    // break;
    // //after = f;
    // //stop = Float.POSITIVE_INFINITY;
    // }
    // if (f > after)
    // break;
    // c++;
    // }

    // Threshold to a mask
    final FloatProcessor fp = new FloatProcessor(width, height, pixels);
    final ShortProcessor sp = (ShortProcessor) fp.convertToShort(true);
    final int t = AutoThreshold.getThreshold(AutoThreshold.Method.OTSU, sp.getHistogram());
    // Utils.display("Blurred", fp);
    for (int i = 0; i < mask.length; i++) {
      if (sp.get(i) >= t) {
        count++;
      }
    }

    // Convert
    final double scale = ((double) count) / mask.length;
    // System.out.printf("Scale = %f\n", scale);
    areaInUm = scale * settings.getSize() * settings.getPixelPitch() * settings.getSize()
        * settings.getPixelPitch() / 1e6;
  }

  /**
   * Extract mask.
   *
   * @param ip the ip
   * @return the mask
   */
  private static int[] extractMask(ImageProcessor ip) {
    // ip = ip.duplicate();
    // ip.setInterpolationMethod(ImageProcessor.BILINEAR);
    // ip = ip.resize(settings.size, settings.size);
    final int[] mask = new int[ip.getPixelCount()];
    for (int i = 0; i < mask.length; i++) {
      mask[i] = ip.get(i);
    }
    return mask;
  }

  private UniformDistribution createUniformDistributionWithPsfWidthBorder() {
    double border = getHwhm() * 3;
    border = FastMath.min(border, settings.getSize() / 4.0);
    return createUniformDistribution(border);
  }

  private SpatialDistribution createFixedDistribution() {
    SpatialDistribution dist;
    dist = new SpatialDistribution() {
      private final double[] xyz = new double[] {settings.getXPosition() / settings.getPixelPitch(),
          settings.getYPosition() / settings.getPixelPitch(),
          settings.getZPosition() / settings.getPixelPitch()};

      @Override
      public double[] next() {
        return xyz;
      }

      @Override
      public boolean isWithinXy(double[] xyz) {
        return true;
      }

      @Override
      public boolean isWithin(double[] xyz) {
        return true;
      }

      @Override
      public void initialise(double[] xyz) {
        // Do nothing
      }
    };
    return dist;
  }

  /**
   * Get the PSF half-width at half-maxima.
   *
   * @return the PSF half-width at half-maxima
   */
  private double getHwhm() {
    if (hwhm == 0) {
      if (psfModelType == PSF_MODEL_IMAGE) {
        hwhm = getImageHwhm();
      } else if (psfModelType == PSF_MODEL_ASTIGMATISM) {
        hwhm = getAstigmatismHwhm();
      } else {
        final double sd = (settings.getEnterWidth()) ? settings.getPsfSd()
            : PsfCalculator.calculateStdDev(settings.getWavelength(),
                settings.getNumericalAperture());

        hwhm = Gaussian2DFunction.SD_TO_HWHM_FACTOR * sd / settings.getPixelPitch();
      }
    }
    return hwhm;
  }

  private PSF createPsf(double psfSd) {
    if (psf == null) {
      if (psfModelType == PSF_MODEL_ASTIGMATISM) {
        // Note: the astigmatismModel may not yet be created so create if necessary.
        // This is used to store the best PSF to use to fit the data.
        AstigmatismModel astigmatismModel = this.astigmatismModel;
        if (astigmatismModel == null) {
          astigmatismModel = AstigmatismModelManager.getModel(settings.getAstigmatismModel());
        }
        if (astigmatismModel == null) {
          throw new IllegalArgumentException(
              "Failed to load model: " + settings.getAstigmatismModel());
        }

        // Assume conversion for simulation
        try {
          if (DoubleEquality.relativeError(astigmatismModel.getNmPerPixel(),
              settings.getPixelPitch()) > 1e-6) {
            // Convert to nm
            astigmatismModel =
                AstigmatismModelManager.convert(astigmatismModel, DistanceUnit.NM, DistanceUnit.NM);
            // Reset pixel pitch. This will draw the spot using the correct size on the different
            // size pixels.
            astigmatismModel =
                astigmatismModel.toBuilder().setNmPerPixel(settings.getPixelPitch()).build();
          }

          // Convert for simulation in pixels
          astigmatismModel = AstigmatismModelManager.convert(astigmatismModel, DistanceUnit.PIXEL,
              DistanceUnit.PIXEL);
        } catch (final ConversionException ex) {
          // Wrap so this can be caught as the same type
          throw new IllegalArgumentException(ex);
        }

        psf = PsfProtosHelper.createPsf(astigmatismModel, DistanceUnit.PIXEL, DistanceUnit.PIXEL);
        psf = psf.toBuilder().setModelName(settings.getAstigmatismModel()).build();
      } else {
        PSF.Builder psfBuilder;
        // Set the PSF as a Gaussian using the width at z=0.
        // In future this could be improved for other PSFs.
        psfBuilder = PsfProtosHelper.defaultOneAxisGaussian2DPSF.toBuilder();
        psfBuilder.getParametersBuilder(PsfHelper.INDEX_SX).setValue(psfSd);
        psf = psfBuilder.build();
      }
    }
    return psf;
  }

  /**
   * Get the PSF standard deviation for a Gaussian using the PSF half-width at half-maxima.
   *
   * @return the PSF standard deviation for a Gaussian
   */
  private double getPsfSd() {
    return getHwhm() / Gaussian2DFunction.SD_TO_HWHM_FACTOR;
  }

  /**
   * Get the PSF half-width at half-maxima from the Image PSF.
   *
   * @return the PSF half-width at half-maxima
   */
  private double getImageHwhm() {
    final ImagePlus imp = WindowManager.getImage(settings.getPsfImageName());
    if (imp == null) {
      IJ.error(TITLE, "Unable to create the PSF model from image: " + settings.getPsfImageName());
      return -1;
    }
    final ImagePSF psfSettings = ImagePsfHelper.fromString(imp.getProperty("Info").toString());
    if (psfSettings == null) {
      IJ.error(TITLE, "Unknown PSF settings for image: " + imp.getTitle());
      return -1;
    }
    if (psfSettings.getFwhm() <= 0) {
      IJ.error(TITLE, "Unknown PSF FWHM setting for image: " + imp.getTitle());
      return -1;
    }
    if (psfSettings.getPixelSize() <= 0) {
      IJ.error(TITLE, "Unknown PSF pixel size setting for image: " + imp.getTitle());
      return -1;
    }

    // The width of the PSF is specified in pixels of the PSF image. Convert to the pixels of the
    // output image
    return 0.5 * psfSettings.getFwhm() * psfSettings.getPixelSize() / settings.getPixelPitch();
  }

  private double getAstigmatismHwhm() {
    AstigmatismModel model = AstigmatismModelManager.getModel(settings.getAstigmatismModel());
    if (model == null) {
      IJ.error(TITLE, "Unknown PSF model: " + settings.getAstigmatismModel());
      return -1;
    }
    try {
      // Get the width at z=0 in pixels
      model = AstigmatismModelManager.convert(model, model.getZDistanceUnit(), DistanceUnit.PIXEL);
      final AstigmatismZModel zModel = AstigmatismModelManager.create(model);
      final double sx = zModel.getSx(0);
      final double sy = zModel.getSy(0);
      return Gaussian2DPeakResultHelper.getStandardDeviation(sx, sy)
          * Gaussian2DFunction.SD_TO_HWHM_FACTOR
          // Scale appropriately
          * model.getNmPerPixel() / settings.getPixelPitch();
    } catch (final ConversionException ex) {
      IJ.error(TITLE, "Unknown PSF FWHM setting for model: " + settings.getAstigmatismModel());
      return -1;
    }
  }

  /**
   * Create distribution within an XY border.
   *
   * @param border the border
   * @return the uniform distribution
   */
  private UniformDistribution createUniformDistribution(double border) {
    final double depth = (settings.getFixedDepth()) ? settings.getDepth() / settings.getPixelPitch()
        : settings.getDepth() / (2 * settings.getPixelPitch());

    // Ensure the focal plane is in the middle of the zDepth
    final double[] max =
        new double[] {settings.getSize() / 2.0 - border, settings.getSize() / 2.0 - border, depth};
    final double[] min = new double[3];
    for (int i = 0; i < 3; i++) {
      min[i] = -max[i];
    }
    if (settings.getFixedDepth()) {
      min[2] = max[2];
    }

    // Try using different distributions:
    final RandomGenerator rand1 = createRandomGenerator();

    if (settings.getDistribution().equals(DISTRIBUTION[UNIFORM_HALTON])) {
      return new UniformDistribution(min, max, rand1.nextInt());
    }

    if (settings.getDistribution().equals(DISTRIBUTION[UNIFORM_SOBOL])) {
      final SobolSequenceGenerator rvg = new SobolSequenceGenerator(3);
      rvg.skipTo(rand1.nextInt());
      return new UniformDistribution(min, max, rvg);
    }

    // Create a distribution using random generators for each dimension
    return new UniformDistribution(min, max, this::createRandomGenerator);
  }

  private SpatialDistribution createConfinementDistribution() {
    if (settings.getDiffusionRate() <= 0 || settings.getFixedFraction() >= 1) {
      return null;
    }

    // Log a warning if the confinement conflicts with the distribution

    if (settings.getConfinement().equals(CONFINEMENT[CONFINEMENT_MASK])) {
      // The mask should be the same
      if (!settings.getDistribution().equals(DISTRIBUTION[MASK])) {
        checkConfiguration("Simulation uses a mask confinement but no mask distribution");
      } else if (!settings.getConfinementMask().equals(settings.getDistributionMask())) {
        checkConfiguration(
            "Simulation uses a mask confinement with a different image to the mask distribution");
      } else if (settings.getConfinementMaskSliceDepth() != settings
          .getDistributionMaskSliceDepth()) {
        checkConfiguration(
            "Simulation uses a mask confinement with a different depth to the mask distribution");
      }

      final ImagePlus imp = WindowManager.getImage(settings.getConfinementMask());
      if (imp != null) {
        return createMaskDistribution(imp, settings.getConfinementMaskSliceDepth(), false);
      }
    } else if (settings.getConfinement().equals(CONFINEMENT[CONFINEMENT_SPHERE])) {
      // This may be an error if the distribution is a mask
      if (settings.getDistribution().equals(DISTRIBUTION[MASK])) {
        checkConfiguration("Simulation uses a mask confinement but a "
            + CONFINEMENT[CONFINEMENT_WITHIN_IMAGE] + " confinement");
      }

      return new SphericalDistribution(settings.getConfinementRadius() / settings.getPixelPitch());
    } else if (settings.getConfinement().equals(CONFINEMENT[CONFINEMENT_WITHIN_IMAGE])) {
      // This may be an error if the distribution is a mask
      if (settings.getDistribution().equals(DISTRIBUTION[MASK])) {
        checkConfiguration("Simulation uses a mask confinement but a "
            + CONFINEMENT[CONFINEMENT_WITHIN_IMAGE] + " confinement");
      }

      return createUniformDistributionWithPsfWidthBorder();
    }

    if (settings.getDistribution().equals(DISTRIBUTION[MASK])) {
      checkConfiguration("Simulation uses a mask confinement but no confinement");
    }
    return null;
  }

  private static void checkConfiguration(String message) {
    final GenericDialog gd = new GenericDialog(TITLE);
    gd.addMessage(TextUtils.wrap("Warning: " + message, 80));
    gd.setOKLabel("Continue");
    gd.showDialog();
    if (gd.wasCanceled()) {
      throw new ConfigurationException(message);
    }
  }

  private SpatialIllumination createIllumination(double intensity, int pulseInterval) {
    if (settings.getIllumination().equals(ILLUMINATION[RADIAL])) {
      if (pulseInterval > 1) {
        return new RadialFalloffIllumination(1, settings.getSize() / 2.0, intensity, pulseInterval);
      }
      return new RadialFalloffIllumination(intensity, settings.getSize() / 2.0);
    }
    uniformBackground = true;
    if (pulseInterval > 1) {
      return new UniformIllumination(1, intensity, pulseInterval);
    }
    return new UniformIllumination(intensity);
  }

  /**
   * Filter those not in the bounds of the distribution.
   *
   * @param localisationSets the localisation sets
   * @return the list
   */
  private List<LocalisationModelSet>
      filterToImageBounds(List<LocalisationModelSet> localisationSets) {
    final List<LocalisationModelSet> newLocalisations = new ArrayList<>(localisationSets.size());
    final SpatialDistribution bounds = createUniformDistribution(0);
    for (final LocalisationModelSet s : localisationSets) {
      if (bounds.isWithinXy(s.toLocalisation().getCoordinates())) {
        newLocalisations.add(s);
      }
    }
    return newLocalisations;
  }

  /**
   * Creates the photon distribution.
   *
   * @return A photon distribution loaded from a file of floating-point values with the specified
   *         population mean.
   */
  private RealDistribution createPhotonDistribution() {
    if (PHOTON_DISTRIBUTION[PHOTON_CUSTOM].equals(settings.getPhotonDistribution())) {
      // Get the distribution file
      final String filename =
          ImageJUtils.getFilename("Photon_distribution", settings.getPhotonDistributionFile());
      if (filename != null) {
        settings.setPhotonDistributionFile(filename);
        try (BufferedReader in = new BufferedReader(new UnicodeReader(
            new FileInputStream(new File(settings.getPhotonDistributionFile())), null))) {
          final StoredDataStatistics stats = new StoredDataStatistics();
          try {
            String str = null;
            double val = 0.0d;
            while ((str = in.readLine()) != null) {
              val = Double.parseDouble(str);
              stats.add(val);
            }
          } finally {
            in.close();
          }

          if (stats.getSum() > 0) {
            // Update the statistics to the desired mean.
            final double scale = settings.getPhotonsPerSecond() / stats.getMean();
            final double[] values = stats.getValues();
            for (int i = 0; i < values.length; i++) {
              values[i] *= scale;
            }

            // TODO - Investigate the limits of this distribution.
            // How far above and below the input data will values be generated.

            // Create the distribution using the recommended number of bins
            final int binCount = stats.getN() / 10;
            final EmpiricalDistribution dist =
                new EmpiricalDistribution(binCount, createRandomGenerator());
            dist.load(values);
            return dist;
          }
        } catch (final IOException ex) {
          // Ignore
        } catch (final NullArgumentException ex) {
          // Ignore
        } catch (final NumberFormatException ex) {
          // Ignore
        }
      }
      ImageJUtils.log("Failed to load custom photon distribution from file: %s. Default to fixed.",
          settings.getPhotonDistributionFile());
    } else if (PHOTON_DISTRIBUTION[PHOTON_UNIFORM].equals(settings.getPhotonDistribution())) {
      if (settings.getPhotonsPerSecond() < settings.getPhotonsPerSecondMaximum()) {
        final UniformRealDistribution dist = new UniformRealDistribution(createRandomGenerator(),
            settings.getPhotonsPerSecond(), settings.getPhotonsPerSecondMaximum());
        return dist;
      }
    } else if (PHOTON_DISTRIBUTION[PHOTON_GAMMA].equals(settings.getPhotonDistribution())) {
      final double scaleParameter = settings.getPhotonsPerSecond() / settings.getPhotonShape();
      final GammaDistribution dist =
          new GammaDistribution(createRandomGenerator(), settings.getPhotonShape(), scaleParameter,
              ExponentialDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);
      return dist;
    } else if (PHOTON_DISTRIBUTION[PHOTON_CORRELATED].equals(settings.getPhotonDistribution())) {
      // No distribution required
      return null;
    }

    settings.setPhotonDistribution(PHOTON_DISTRIBUTION[PHOTON_FIXED]);
    return null;
  }

  private List<LocalisationModelSet> combineSimulationSteps(List<LocalisationModel> localisations) {
    // Allow fractional integration steps
    final double simulationStepsPerFrame =
        (settings.getStepsPerSecond() * settings.getExposureTime()) / 1000.0;

    final List<LocalisationModelSet> newLocalisations =
        new ArrayList<>((int) (localisations.size() / simulationStepsPerFrame));

    // System.out.printf("combineSimulationSteps @ %f\n", simulationStepsPerFrame);

    // final double gain = new CreateDataSettingsHelper(settings).getTotalGainSafe();
    sortLocalisationsByIdThenTime(localisations);
    final int[] idList = getIds(localisations);
    movingMolecules = new TIntHashSet(idList.length);
    int index = 0;
    for (final int id : idList) {
      final int fromIndex = findIndexById(localisations, index, id);
      if (fromIndex > -1) {
        final int toIndex = findLastIndexById(localisations, fromIndex, id);
        final List<LocalisationModel> subset = localisations.subList(fromIndex, toIndex + 1);
        index = toIndex;

        // Store the IDs of any moving molecules
        if (isMoving(subset)) {
          movingMolecules.add(id);
        }

        // The frames may be longer or shorter than the simulation steps. Allocate the step
        // proportionately to each frame it overlaps:
        //
        // Steps: |-- 0 --|-- 1 --|-- 2 --|--
        // Frames: |--- 0 ---|--- 1 ---|--- 2 ---|
        //
        // ^ ^
        // | |
        // | End frame
        // |
        // Start frame

        final double firstFrame = getStartFrame(subset.get(0), simulationStepsPerFrame);
        final double lastFrame =
            getEndFrame(subset.get(subset.size() - 1), simulationStepsPerFrame);

        // Get the first frame offset and allocate space to store all potential frames
        final int intFirstFrame = (int) firstFrame;
        final int intLastFrame = (int) Math.ceil(lastFrame);
        final LocalisationModelSet[] sets =
            new LocalisationModelSet[intLastFrame - intFirstFrame + 1];

        // Process each step
        for (final LocalisationModel l : subset) {
          // Get the fractional start and end frames
          final double startFrame = getStartFrame(l, simulationStepsPerFrame);
          final double endFrame = getEndFrame(l, simulationStepsPerFrame);

          // Round down to get the actual frames that are overlapped
          final int start = (int) startFrame;
          int end = (int) endFrame;

          // Check if the span covers a fraction of the end frame, otherwise decrement to ignore
          // that frame
          if (end > start && endFrame == end) {
            // E.g. convert
            // Steps: |-- 0 --|
            // Frames: |- 0 -|- 1 -|- 2 -|
            // to
            // Steps: |-- 0 --|
            // Frames: |- 0 -|- 1 -|
            end--;
          }

          if (start == end) {
            // If the step falls within one frame then add it to the set
            final int tIndex = start - intFirstFrame;
            if (sets[tIndex] == null) {
              sets[tIndex] = new LocalisationModelSet(id, start);
            }
            sets[tIndex].add(l);
          } else {
            // Add the localisation to all the frames that the step spans
            final double total = endFrame - startFrame;
            // double t = 0;
            for (int frame = start; frame <= end; frame++) {
              // Get the fraction to allocate to this frame
              double fraction;
              int state = (l.isContinuous() ? LocalisationModel.CONTINUOUS : 0);
              if (frame == start) {
                state |=
                    LocalisationModel.NEXT | (l.hasPrevious() ? LocalisationModel.PREVIOUS : 0);
                // |-----|====|
                // | | ceil(startFrame)
                // | startFrame
                // start
                if (startFrame == start) {
                  fraction = 1;
                } else {
                  fraction = (Math.ceil(startFrame) - startFrame);
                }
              } else if (frame == end) {
                state |= LocalisationModel.PREVIOUS | (l.hasNext() ? LocalisationModel.NEXT : 0);
                // |=====|----|
                // | |
                // | endFrame
                // end
                fraction = (endFrame - end);
              } else {
                state |= LocalisationModel.CONTINUOUS;
                fraction = 1;
              }

              // t += fraction;

              // Add to the set
              final int tIndex = frame - intFirstFrame;
              if (sets[tIndex] == null) {
                sets[tIndex] = new LocalisationModelSet(id, frame);
              }
              sets[tIndex].add(getFraction(l, fraction / total, state));
            }
            // if (t < total * 0.98)
            // {
            // System.out.printf("Total error %g < %g : %f (%d) -> %f (%d)\n", t, total, startFrame,
            // start, endFrame, end);
            // }
          }
        }

        LocalisationModelSet previous = null;
        for (int i = 0; i < sets.length; i++) {
          if (sets[i] != null) {
            sets[i].setPrevious(previous);

            // Create a data array and store the current intensity.
            // This is used later to filter based on SNR
            sets[i].setData(new double[] {0, 0, 0, 0, sets[i].getIntensity() // * gain
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
   * Check if any of the coordinates for the subset are different.
   *
   * @param subset the subset
   * @return True if the coordinates move
   */
  private static boolean isMoving(List<LocalisationModel> subset) {
    if (subset.size() < 2) {
      return false;
    }
    final double[] xyz = subset.get(0).getCoordinates();
    for (int i = 1; i < subset.size(); i++) {
      final double[] xyz2 = subset.get(i).getCoordinates();
      for (int j = 0; j < 3; j++) {
        if (xyz[j] != xyz2[j]) {
          return true;
        }
      }
    }
    return false;
  }

  /**
   * Get the simulation frame start point for the localisation.
   *
   * @param localisationModel the localisation model
   * @param simulationStepsPerFrame the simulation steps per frame
   * @return the start frame
   */
  private static double getStartFrame(LocalisationModel localisationModel,
      double simulationStepsPerFrame) {
    // Time is 1-based (not 0)
    return (localisationModel.getTime() - 1) / simulationStepsPerFrame + 1;
  }

  /**
   * Get the simulation frame end point for the localisation.
   *
   * @param localisationModel the localisation model
   * @param simulationStepsPerFrame the simulation steps per frame
   * @return the end frame
   */
  private static double getEndFrame(LocalisationModel localisationModel,
      double simulationStepsPerFrame) {
    // Time is 1-based (not 0)
    return (localisationModel.getTime()) / simulationStepsPerFrame + 1;
  }

  /**
   * Create a new localisation model with the same id and position but with a fraction of the
   * intensity and the specified state.
   *
   * @param lm the localisation model
   * @param fraction the fraction
   * @param state the state
   * @return new localisation model
   */
  private static LocalisationModel getFraction(LocalisationModel lm, double fraction, int state) {
    return new LocalisationModel(lm.getId(), lm.getTime(), lm.getCoordinates(),
        lm.getIntensity() * fraction, state);
  }

  private int[] getIds(List<LocalisationModel> localisations) {
    if (localisations.isEmpty()) {
      return new int[0];
    }
    final TIntArrayList ids = new TIntArrayList(settings.getParticles());
    // Assume the localisations are sorted by id
    int id = localisations.get(0).getId();
    ids.add(id);
    for (final LocalisationModel l : localisations) {
      if (id != l.getId()) {
        id = l.getId();
        ids.add(id);
      }
    }
    return ids.toArray();
  }

  // StoredDataStatistics rawPhotons = new StoredDataStatistics();
  // StoredDataStatistics drawPhotons = new StoredDataStatistics();

  // private synchronized void addRaw(double d)
  // {
  // //rawPhotons.add(d);
  // }
  //
  // private synchronized void addDraw(double d)
  // {
  // //drawPhotons.add(d);
  // }

  /**
   * Create an image from the localisations using the configured PSF width. Draws a new stack image.
   *
   * <p>Note that the localisations are filtered using the signal. The input list of localisations
   * will be updated.
   *
   * @param localisationSets the localisation sets
   * @return The localisations
   */
  private List<LocalisationModel> drawImage(final List<LocalisationModelSet> localisationSets) {
    if (localisationSets.isEmpty()) {
      return null;
    }

    // Create a new list for all localisation that are drawn (i.e. pass the signal filters)
    List<LocalisationModelSet> newLocalisations =
        Collections.synchronizedList(new ArrayList<LocalisationModelSet>(localisationSets.size()));
    photonsRemoved = new AtomicInteger();
    removedT1 = new AtomicInteger();
    removedTn = new AtomicInteger();
    photonStats = new SummaryStatistics();

    // Add drawn spots to memory
    results = new MemoryPeakResults();
    final CalibrationWriter c = new CalibrationWriter();
    c.setDistanceUnit(DistanceUnit.PIXEL);
    c.setIntensityUnit(IntensityUnit.PHOTON);
    c.setNmPerPixel(settings.getPixelPitch());
    c.setExposureTime(settings.getExposureTime());

    c.setCameraType(settings.getCameraType());
    if (settings.getCameraType() == CameraType.SCMOS) {
      c.setCameraModelName(settings.getCameraModelName());
    } else {
      final CreateDataSettingsHelper helper = new CreateDataSettingsHelper(settings);
      c.setCountPerPhoton(helper.getTotalGainSafe());
      c.setBias(settings.getBias());
      c.setReadNoise(helper.getReadNoiseInCounts());
      c.setQuantumEfficiency(helper.getQuantumEfficiency());
    }

    results.setCalibration(c.getCalibration());
    results.setSortAfterEnd(true);
    results.begin();

    maxT = localisationSets.get(localisationSets.size() - 1).getTime();

    // Display image
    final ImageStack stack = new ImageStack(settings.getSize(), settings.getSize(), maxT);

    final double psfSd = getPsfSd();
    if (psfSd <= 0) {
      return null;
    }
    PsfModel psfModel = createPsfModel(localisationSets);
    if (psfModel == null) {
      return null;
    }

    // Create the camera noise model
    createPerPixelCameraModelData(cameraModel);
    createBackgroundPixels();

    final int threadCount = Prefs.getThreads();

    IJ.showStatus("Drawing image ...");

    // Multi-thread for speed
    final PeakResults syncResults = SynchronizedPeakResults.create(results, threadCount);
    ExecutorService threadPool = Executors.newFixedThreadPool(threadCount);
    List<Future<?>> futures = new LinkedList<>();

    // Count all the frames to process
    frame = 0;
    totalFrames = maxT;

    // Collect statistics on the number of photons actually simulated

    // Process all frames
    int index = 0;
    int lastT = -1;
    for (final LocalisationModelSet l : localisationSets) {
      if (ImageJUtils.isInterrupted()) {
        break;
      }
      if (l.getTime() != lastT) {
        lastT = l.getTime();
        futures.add(threadPool.submit(new ImageGenerator(localisationSets, newLocalisations, index,
            lastT, copyPsfModel(psfModel), syncResults, stack, poissonNoise,
            new RandomDataGenerator(createRandomGenerator()))));
      }
      index++;
    }
    // Finish processing data
    ConcurrencyUtils.waitForCompletionUnchecked(futures);
    futures.clear();
    if (ImageJUtils.isInterrupted()) {
      IJ.showProgress(1);
      return null;
    }

    // Do all the frames that had no localisations
    float[] limits = null;
    for (int t = 1; t <= maxT; t++) {
      if (ImageJUtils.isInterrupted()) {
        break;
      }
      final Object pixels = stack.getPixels(t);
      if (pixels == null) {
        futures.add(threadPool.submit(
            new ImageGenerator(localisationSets, newLocalisations, maxT, t, null, syncResults,
                stack, poissonNoise, new RandomDataGenerator(createRandomGenerator()))));
      } else if (limits == null) {
        limits = MathUtils.limits((float[]) pixels);
      }
    }

    // Finish
    ConcurrencyUtils.waitForCompletionUnchecked(futures);
    threadPool.shutdown();
    IJ.showProgress(1);
    if (ImageJUtils.isInterrupted() || limits == null) {
      return null;
    }
    results.end();

    // Clear memory
    psfModel = null;
    threadPool = null;
    futures.clear();
    futures = null;

    if (photonsRemoved.get() > 0) {
      ImageJUtils.log("Removed %d localisations with less than %.1f rendered photons",
          photonsRemoved.get(), settings.getMinPhotons());
    }
    if (removedT1.get() > 0) {
      ImageJUtils.log("Removed %d localisations with no neighbours @ SNR %.2f", removedT1.get(),
          settings.getMinSnrT1());
    }
    if (removedTn.get() > 0) {
      ImageJUtils.log("Removed %d localisations with valid neighbours @ SNR %.2f", removedTn.get(),
          settings.getMinSnrTN());
    }
    if (photonStats.getN() > 0) {
      ImageJUtils.log("Average photons rendered = %s +/- %s",
          MathUtils.rounded(photonStats.getMean()),
          MathUtils.rounded(photonStats.getStandardDeviation()));
    }

    // System.out.printf("rawPhotons = %f\n", rawPhotons.getMean());
    // System.out.printf("drawPhotons = %f\n", drawPhotons.getMean());
    // new HistogramPlotBuilder("draw photons", drawPhotons, "photons", true, 0, 1000);

    // Update with all those localisation that have been drawn
    localisationSets.clear();
    localisationSets.addAll(newLocalisations);
    newLocalisations = null;

    IJ.showStatus("Displaying image ...");

    ImageStack newStack = stack;

    if (!settings.getRawImage()) {
      // Get the global limits and ensure all values can be represented
      final Object[] imageArray = stack.getImageArray();
      limits = MathUtils.limits((float[]) imageArray[0]);
      for (int j = 1; j < imageArray.length; j++) {
        limits = MathUtils.limits(limits, (float[]) imageArray[j]);
      }
      // float limits0 = limits[0];
      final float limits0 = 0; // Leave bias in place
      // Check if the image will fit in a 16-bit range
      if ((limits[1] - limits0) < 65535) {
        // Convert to 16-bit
        newStack = new ImageStack(stack.getWidth(), stack.getHeight(), stack.getSize());
        // Account for rounding
        final float min = (float) (limits0 - 0.5);
        for (int j = 0; j < imageArray.length; j++) {
          final float[] image = (float[]) imageArray[j];
          final short[] pixels = new short[image.length];
          for (int k = 0; k < pixels.length; k++) {
            pixels[k] = (short) (image[k] - min);
          }
          newStack.setPixels(pixels, j + 1);
          // Free memory
          imageArray[j] = null;
          // Attempt to stay within memory (check vs 32MB)
          if (MemoryUtils.getFreeMemory() < 33554432L) {
            MemoryUtils.runGarbageCollectorOnce();
          }
        }
        for (int k = 2; k-- > 0;) {
          limits[k] = (float) Math.floor(limits[k] - min);
        }
      } else {
        // Keep as 32-bit but round to whole numbers
        for (int j = 0; j < imageArray.length; j++) {
          final float[] pixels = (float[]) imageArray[j];
          for (int k = 0; k < pixels.length; k++) {
            pixels[k] = Math.round(pixels[k]);
          }
        }
        for (int k = 2; k-- > 0;) {
          limits[k] = Math.round(limits[k]);
        }
      }
    }

    // Show image
    final ImagePlus imp = ImageJUtils.display(CREATE_DATA_IMAGE_TITLE, newStack);

    final ij.measure.Calibration cal = new ij.measure.Calibration();
    String unit = "nm";
    double unitPerPixel = settings.getPixelPitch();
    if (unitPerPixel > 100) {
      unit = "um";
      unitPerPixel /= 1000.0;
    }
    cal.setUnit(unit);
    cal.pixelHeight = cal.pixelWidth = unitPerPixel;
    final Rectangle bounds = cameraModel.getBounds();
    if (bounds != null) {
      cal.xOrigin = -bounds.x;
      cal.yOrigin = -bounds.y;
    }
    imp.setCalibration(cal);

    imp.setDimensions(1, 1, newStack.getSize());
    imp.setDisplayRange(limits[0], limits[1]);
    // imp.resetDisplayRange();
    imp.updateAndDraw();

    saveImage(imp);

    final IJImageSource imageSource = new IJImageSource(imp);
    // Shift simulation image source to correct location
    results.setSource(imageSource);
    results.setName(CREATE_DATA_IMAGE_TITLE + " (" + TITLE + ")");
    // Bounds are relative to the image source
    results.setBounds(new Rectangle(settings.getSize(), settings.getSize()));

    results.setPsf(createPsf(psfSd));
    MemoryPeakResults.addResults(results);

    setBenchmarkResults(imp, results);

    if (benchmarkMode && benchmarkParameters != null) {
      benchmarkParameters.setPhotons(results);
    }

    final List<LocalisationModel> localisations = toLocalisations(localisationSets);

    savePulses(localisations, results);

    // Saved the fixed and moving localisations into different datasets
    saveFixedAndMoving(results);

    saveCompoundMolecules(results);

    return localisations;
  }

  private synchronized void addPhotons(double photons) {
    photonStats.addValue(photons);
  }

  /**
   * Create a PSF model from the image that contains all the z-slices needed to draw the given
   * localisations.
   *
   * @param localisationSets the localisation sets
   * @return the image PSF model
   */
  private ImagePsfModel createImagePsf(List<LocalisationModelSet> localisationSets) {
    final ImagePlus imp = WindowManager.getImage(settings.getPsfImageName());
    if (imp == null) {
      IJ.error(TITLE, "Unable to create the PSF model from image: " + settings.getPsfImageName());
      return null;
    }
    try {
      final ImagePSF psfSettings = ImagePsfHelper.fromString(imp.getProperty("Info").toString());
      if (psfSettings == null) {
        throw new IllegalStateException("Unknown PSF settings for image: " + imp.getTitle());
      }

      // Check all the settings have values
      if (psfSettings.getPixelSize() <= 0) {
        throw new IllegalStateException(
            "Missing nmPerPixel calibration settings for image: " + imp.getTitle());
      }
      if (psfSettings.getPixelDepth() <= 0) {
        throw new IllegalStateException(
            "Missing nmPerSlice calibration settings for image: " + imp.getTitle());
      }
      if (psfSettings.getCentreImage() <= 0) {
        throw new IllegalStateException(
            "Missing zCentre calibration settings for image: " + imp.getTitle());
      }
      if (psfSettings.getFwhm() <= 0) {
        throw new IllegalStateException(
            "Missing FWHM calibration settings for image: " + imp.getTitle());
      }

      // To save memory construct the Image PSF using only the slices that are within
      // the depth of field of the simulation
      double minZ = Double.POSITIVE_INFINITY;
      double maxZ = Double.NEGATIVE_INFINITY;
      for (final LocalisationModelSet l : localisationSets) {
        for (final LocalisationModel m : l.getLocalisations()) {
          final double z = m.getZ();
          if (minZ > z) {
            minZ = z;
          }
          if (maxZ < z) {
            maxZ = z;
          }
        }
      }

      final int nSlices = imp.getStackSize();
      // z-centre should be an index and not the ImageJ slice number so subtract 1
      final int zCentre = psfSettings.getCentreImage() - 1;

      // Calculate the start/end slices to cover the depth of field
      // This logic must match the ImagePSFModel.
      final double unitsPerSlice = psfSettings.getPixelDepth() / settings.getPixelPitch();
      // We assume the PSF was imaged axially with increasing z-stage position (moving the stage
      // closer to the objective). Thus higher z-coordinate are for higher slice numbers.
      int lower = (int) Math.round(minZ / unitsPerSlice) + zCentre;
      int upper = (int) Math.round(maxZ / unitsPerSlice) + zCentre;
      // Add extra to the range so that gradients can be computed.
      lower--;
      upper++;
      upper = (upper < 0) ? 0 : (upper >= nSlices) ? nSlices - 1 : upper;
      lower = (lower < 0) ? 0 : (lower >= nSlices) ? nSlices - 1 : lower;

      // We cannot just extract the correct slices since the
      // Image PSF requires the z-centre for normalisation
      if (!(lower <= zCentre && upper >= zCentre)) {
        // Ensure we include the zCentre
        lower = Math.min(lower, zCentre);
        upper = Math.max(upper, zCentre);
      }

      final double noiseFraction = 1e-3;
      final float[][] image = extractImageStack(imp, lower, upper);
      final ImagePsfModel model = new ImagePsfModel(image, zCentre - lower,
          psfSettings.getPixelSize() / settings.getPixelPitch(), unitsPerSlice, noiseFraction);

      // Add the calibrated centres. The map will not be null
      final Map<Integer, Offset> map = psfSettings.getOffsetsMap();
      if (!map.isEmpty()) {
        final int sliceOffset = lower + 1;
        for (final Entry<Integer, Offset> entry : map.entrySet()) {
          model.setRelativeCentre(entry.getKey() - sliceOffset, entry.getValue().getCx(),
              entry.getValue().getCy());
        }
      } else {
        // Use the CoM if present
        final double cx = psfSettings.getXCentre();
        final double cy = psfSettings.getYCentre();
        if (cx != 0 || cy != 0) {
          for (int slice = 0; slice < image.length; slice++) {
            model.setCentre(slice, cx, cy);
          }
        }
      }

      // Initialise the HWHM table so that it can be cloned
      model.initialiseHwhm();

      return model;
    } catch (final Exception ex) {
      IJ.error(TITLE, "Unable to create the image PSF model:\n" + ex.getMessage());
      return null;
    }
  }

  /**
   * Extract the image stack using a range of stack indices. The index is 0-based.
   *
   * @param imp the image
   * @param start the start index
   * @param end the end index
   * @return the float image stack
   */
  public static float[][] extractImageStack(ImagePlus imp, int start, int end) {
    final int size = end - start + 1;
    final ImageStack stack = imp.getImageStack();
    final float[][] image = new float[size][];
    for (int i = 0; i < image.length; i++) {
      image[i] = (float[]) stack.getProcessor(i + start + 1).toFloat(0, null).getPixels();
    }
    return image;
  }

  private PsfModel createPsfModel(List<LocalisationModelSet> localisationSets) {
    // Allow reuse of the cached model
    if (psfModelCache != null) {
      return psfModelCache;
    }

    if (psfModelType == PSF_MODEL_IMAGE) {
      return createImagePsf(localisationSets);
    }

    if (psfModelType == PSF_MODEL_ASTIGMATISM) {
      astigmatismModel = AstigmatismModelManager.getModel(settings.getAstigmatismModel());
      if (astigmatismModel == null) {
        throw new IllegalArgumentException(
            "Failed to load model: " + settings.getAstigmatismModel());
      }
      // Convert for simulation
      try {
        if (DoubleEquality.relativeError(astigmatismModel.getNmPerPixel(),
            settings.getPixelPitch()) > 1e-6) {
          final String message = String.format(
              "Astigmatism model '%s' calibration (%s nm) does not match pixel pitch (%s nm)",
              settings.getAstigmatismModel(), MathUtils.rounded(astigmatismModel.getNmPerPixel()),
              MathUtils.rounded(settings.getPixelPitch()));
          // Optionally convert
          final GenericDialog gd = new GenericDialog(TITLE);
          gd.addMessage(
              TextUtils.wrap(message + ". Created data is not suitable for fitting.", 80));
          gd.addMessage(TextUtils.wrap("Click OK to continue anyway (i.e. draw the spot using the "
              + "correct nm width on the different sized pixels).", 80));
          gd.showDialog();
          if (gd.wasCanceled()) {
            throw new IllegalArgumentException(message);
          }
          // Convert to nm
          astigmatismModel =
              AstigmatismModelManager.convert(astigmatismModel, DistanceUnit.NM, DistanceUnit.NM);
          // Reset pixel pitch. This will draw the spot using the correct size on the different size
          // pixels.
          astigmatismModel =
              astigmatismModel.toBuilder().setNmPerPixel(settings.getPixelPitch()).build();
        }

        // Convert for simulation in pixels
        astigmatismModel = AstigmatismModelManager.convert(astigmatismModel, DistanceUnit.PIXEL,
            DistanceUnit.PIXEL);
        return new GaussianPsfModel(AstigmatismModelManager.create(astigmatismModel));
      } catch (final ConversionException ex) {
        // Wrap so this can be caught as the same type
        throw new IllegalArgumentException(ex);
      }
    }

    final double d = settings.getDepthOfFocus() / settings.getPixelPitch();

    if (psfModelType == PSF_MODEL_GAUSSIAN) {
      final double sd = getPsfSd();
      final double gamma = 0;
      final HoltzerAstigmatismZModel zModel =
          HoltzerAstigmatismZModel.create(sd, sd, gamma, d, 0, 0, 0, 0);
      final GaussianPsfModel m = new GaussianPsfModel(createRandomGenerator(), zModel);
      // m.setRange(10);
      return m;
    }

    // Default to Airy pattern
    final double width = getPsfSd() / PsfCalculator.AIRY_TO_GAUSSIAN;
    final AiryPsfModel m = new AiryPsfModel(createRandomGenerator(), width, width, d);
    m.setRing(2);
    return m;
  }

  private PsfModel createPsfModel(double[] xyz) {
    // Create a set with a single model
    final List<LocalisationModelSet> localisationSets = new TurboList<>(1);
    final LocalisationModelSet set = new LocalisationModelSet(0, 0);
    final LocalisationModel m = new LocalisationModel(0, 0, xyz, 1, 0);
    set.add(m);
    localisationSets.add(set);
    return createPsfModel(localisationSets);
  }

  private PsfModel copyPsfModel(PsfModel psfModel) {
    return psfModel.copy(createRandomGenerator());
  }

  private synchronized void showProgress() {
    IJ.showProgress(frame++, totalFrames);
  }

  private static class Spot {
    final double[] psf;
    final int x0min;
    final int x0max;
    final int x1min;
    final int x1max;
    final int[] samplePositions;

    public Spot(double[] psf, int x0min, int x0max, int x1min, int x1max, int[] samplePositions) {
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
   * The read noise. This is stored in electrons even though camera read noise is measured in ADUs.
   * This allows it to be used to compute the combined background noise (with the background photon
   * shot noise) for each localisation.
   */
  private float[] readNoise;

  private void createPerPixelCameraModelData(CameraModel cameraModel) {
    if (readNoise != null) {
      return;
    }

    // Note: Store read noise in electrons for SNR computation. The gain is later applied back.

    if (cameraModel.isPerPixelModel()) {
      final Rectangle bounds = cameraModel.getBounds();
      readNoise = cameraModel.getVariance(bounds);
      final float[] gain = cameraModel.getGain(bounds);
      for (int i = 0; i < readNoise.length; i++) {
        readNoise[i] = (float) (Math.sqrt(readNoise[i]) / gain[i]);
      }
    } else {
      // Use a dummy coordinate to find out the fixed variance and gain
      final float variance = cameraModel.getVariance(0, 0);
      final float gain = cameraModel.getGain(0, 0);

      // Avoid sqrt on all the same value
      this.readNoise = new float[settings.getSize() * settings.getSize()];
      Arrays.fill(this.readNoise, (float) (Math.sqrt(variance) / gain));
    }

    // Remove if it will have no effect
    readNoise = nullIfZero(readNoise);
  }

  private static float[] nullIfZero(float[] data) {
    for (final float value : data) {
      if (value != 0) {
        return data;
      }
    }
    return null;
  }

  /**
   * Creates the CCD camera model.
   *
   * <p>Note that the model only has camera gain applied thus the normalised variance is the read
   * noise in electrons. This is standard for a CCD model but omits the EM-gain for an EM-CCD model.
   * This model is intended to be used to generate electron noise during the simulation.
   *
   * @return the camera model
   */
  private CameraModel createCcdCameraModel() {
    final float bias = settings.getBias();
    float gain = 1f;

    if (settings.getCameraGain() != 0) {
      gain = (float) settings.getCameraGain();
    }

    final float readNoise = (float) new CreateDataSettingsHelper(settings).getReadNoiseInCounts();

    return new CcdCameraModel(bias, gain, (float) MathUtils.pow2(readNoise));
  }

  /**
   * Use a runnable for the image generation to allow multi-threaded operation. Input parameters
   * that are manipulated should have synchronized methods.
   */
  private class ImageGenerator implements Runnable {
    final List<LocalisationModelSet> localisations;
    final List<LocalisationModelSet> newLocalisations;
    final int startIndex;
    final int time;
    final PsfModel psfModel;
    final PeakResults results;
    final ImageStack stack;
    final boolean poissonNoise;
    final RandomDataGenerator random;
    final double emGain;
    final double qe;

    public ImageGenerator(final List<LocalisationModelSet> localisationSets,
        List<LocalisationModelSet> newLocalisations, int startIndex, int time, PsfModel psfModel,
        PeakResults results, ImageStack stack, boolean poissonNoise, RandomDataGenerator random) {
      this.localisations = localisationSets;
      this.newLocalisations = newLocalisations;
      this.startIndex = startIndex;
      this.time = time;
      this.psfModel = psfModel;
      this.results = results;
      this.stack = stack;
      this.poissonNoise = poissonNoise;
      this.random = random;
      // This could be >=1 but the rest of the code ignores EM-gain if it is <=1
      emGain = (settings.getCameraType() == CameraType.EMCCD && settings.getEmGain() > 1)
          ? settings.getEmGain()
          : 0;
      qe = getQuantumEfficiency();
    }

    @Override
    public void run() {
      if (ImageJUtils.isInterrupted()) {
        return;
      }

      showProgress();

      final boolean checkSnr = minSnrT1 > 0 || minSnrTn > 0;

      // Adjust XY dimensions since they are centred on zero
      final double xoffset = settings.getSize() * 0.5;

      // The simulation applies the following:
      // Poisson - Photons emitted
      // Binomial - Photons to electrons (QE) on the camera chip
      // Gamma* - Optional EM-gain amplification (electrons)
      // Gaussian* - Camera gain amplification with noise (count)
      // *The output is rounded to remain a discrete distribution

      // In order to analytically determine the noise around each localisation
      // the background is simulated separately. This is added to the
      // simulation for the rendered photons.
      // This is possible because the Gamma(shape,scale) distribution can be added:
      // Gamma(a+b,scale) = Gamma(a,scale) + Gamma(b,scale)
      // See: https://en.wikipedia.org/wiki/Gamma_distribution#Summation

      final float[] background = createBackground(random);
      if (settings.getBackground() > 0) {
        amplify(background);
      }

      final float[] image = new float[background.length];

      // Create read noise now so that we can calculate the true background noise.
      // Note: The read noise is in electrons.
      if (readNoise != null) {
        final RandomGenerator r = random.getRandomGenerator();
        for (int i = 0; i < background.length; i++) {
          background[i] += (float) (readNoise[i] * r.nextGaussian());
        }
      }

      // Extract the localisations and draw if we have a PSF model
      final int fromIndex = findIndexByTime(localisations, startIndex, time);
      if (fromIndex > -1 && psfModel != null) {
        final int toIndex = findLastIndexByTime(localisations, fromIndex, time);
        final List<LocalisationModelSet> subset = localisations.subList(fromIndex, toIndex + 1);
        final float[] data = new float[settings.getSize() * settings.getSize()];
        for (final LocalisationModelSet localisationSet : subset) {
          if (ImageJUtils.isInterrupted()) {
            return;
          }

          if (localisationSet.size() == 0) {
            continue;
          }

          // Draw each localisation in the set. Store the PSF so we can remove it later
          double totalPhotonsRendered = 0;
          final Spot[] spots = new Spot[localisationSet.size()];
          int spotCount = 0;
          for (final LocalisationModel localisation : localisationSet.getLocalisations()) {
            // Adjust to centre of image
            final double[] xyz = localisation.getCoordinates();
            xyz[0] += xoffset;
            xyz[1] += xoffset;

            // addRaw(localisation.getIntensity());

            double photonsRendered = 0;
            double intensity = localisation.getIntensity();
            int[] samplePositions = null;
            if (intensity > 0) {
              // TODO record all the positions we draw and the number of photons.
              // This can be used later to compute the Fisher information for each spot.
              if (poissonNoise) {
                final int samples = (int) random.nextPoisson(intensity);
                intensity = samples;
                photonsRendered = psfModel.sample3D(data, settings.getSize(), settings.getSize(),
                    samples, localisation.getX(), localisation.getY(), localisation.getZ());
                samplePositions = psfModel.getSamplePositions();
              } else {
                photonsRendered =
                    psfModel.create3D(data, settings.getSize(), settings.getSize(), intensity,
                        localisation.getX(), localisation.getY(), localisation.getZ(), false);
              }
            }
            // addDraw(photons);
            if (photonsRendered > 0) {
              totalPhotonsRendered += photonsRendered;
              spots[spotCount++] = new Spot(psfModel.getPsf(), psfModel.getX0min(),
                  psfModel.getX0max(), psfModel.getX1min(), psfModel.getX1max(), samplePositions);
            }
          }

          // Skip if nothing has been drawn. Note that if the localisation set is skipped then the
          // intensity must be set to zero to prevent the SNR checks using the eliminated
          // neighbours.
          if (totalPhotonsRendered == 0) {
            photonsRemoved.incrementAndGet();
            localisationSet.setData(new double[5]);
            continue;
          }
          if (totalPhotonsRendered < minPhotons) {
            photonsRemoved.incrementAndGet();
            for (int i = 0; i < spotCount; i++) {
              erase(data, spots[i]);
            }
            localisationSet.setData(new double[5]);
            continue;
          }

          addPhotons(totalPhotonsRendered);

          final LocalisationModel localisation = localisationSet.toLocalisation();

          // Add to memory.
          // Use the actual intensity (not the total photons rendered)
          final float intensity = (float) localisation.getIntensity();
          final float x = (float) localisation.getX();
          final float y = (float) localisation.getY();
          final float z = (float) localisation.getZ();
          // 0.5 is the centre of the pixel so just round down.
          final int origX = (int) x;
          final int origY = (int) y;
          // Background and noise should be calculated using the
          // region covered by the PSF.
          final double[] localStats = getStatistics(backgroundPixels, background, origX, origY);
          final double backgroundInPhotons = localStats[0];

          // Note: The width estimate does not account for diffusion
          float sx;
          float sy;
          if (psfModel instanceof GaussianPsfModel) {
            final GaussianPsfModel m = (GaussianPsfModel) psfModel;
            sx = (float) m.getS0();
            sy = (float) m.getS1();
          } else if (psfModel instanceof AiryPsfModel) {
            final AiryPsfModel m = (AiryPsfModel) psfModel;
            sx = (float) (m.getW0() * AiryPattern.FACTOR);
            sy = (float) (m.getW1() * AiryPattern.FACTOR);
          } else if (psfModel instanceof ImagePsfModel) {
            final ImagePsfModel m = (ImagePsfModel) psfModel;
            sx = (float) (m.getHwhm0() / Gaussian2DFunction.SD_TO_HWHM_FACTOR);
            sy = (float) (m.getHwhm1() / Gaussian2DFunction.SD_TO_HWHM_FACTOR);
          } else {
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

          // The variance of the background image is currently in electrons^2
          final double totalNoise = Math.sqrt(localStats[3]);

          // Convert noise from electrons back to photons for convenience when computing SNR using
          // the signal in photons.
          double totalNoiseInPhotons = totalNoise / qe;

          // Account for EM-gain.
          // This makes the total noise in photons the equivalent of:
          // [noise in counts] / [total gain]
          // [noise in counts] / [EM-gain * camera-gain * qe]
          // [[noise in counts] / camera-gain] / [EM-gain * qe]
          // [noise in electrons] / [EM-gain * qe]
          if (emGain != 0) {
            totalNoiseInPhotons /= emGain;
          }

          // System.out.printf("Noise = %g e-\n", totalNoiseInPhotons);

          // Ensure the new data is added before the intensity is updated. This avoids
          // syncronisation clashes in the getIntensity(...) function.
          // Use the total photons rendered for signal filtering.
          // [0] = background (photons)
          // [1] = total noise (photons)
          // [2] = Gaussian Sx
          // [3] = Gaussian Sy
          // [4] = total intensity (photons)
          localisationSet.setData(new double[] {backgroundInPhotons, totalNoiseInPhotons, sx, sy,
              totalPhotonsRendered});

          if (checkSnr
              && badLocalisation(localisationSet, totalPhotonsRendered, totalNoiseInPhotons)) {
            for (int i = 0; i < spotCount; i++) {
              erase(data, spots[i]);
            }
            localisationSet.setData(new double[5]);
            continue;
          }

          newLocalisations.add(localisationSet);
          // The parameters should have the background and intensity in
          // photons before noise (as a perfect result).
          final float[] params = Gaussian2DPeakResultHelper
              .createTwoAxisParams((float) backgroundInPhotons, intensity, x, y, z, sx, sy);
          final float u =
              (float) Gaussian2DPeakResultHelper.getMeanSignalUsingP05(intensity, sx, sy);
          // Use extended result to store the ID
          results.add(new IdPeakResult(time, origX, origY, 0, localisation.getZ(),
              (float) totalNoiseInPhotons, u, params, null, localisationSet.getId()));
        }

        for (int i = 0; i < image.length; i++) {
          image[i] += data[i];
        }
      }

      amplify(image);

      // Combine with background + read noise (in electrons)
      for (int i = 0; i < image.length; i++) {
        image[i] += background[i];
      }

      // Apply camera gain. Note that the noise component of the camera gain IS the
      // read noise. Thus the read noise is expected to be different for each camera gain.
      // Also add the bias.
      cameraModel.applyGainAndBias(cameraModel.getBounds(), image);

      // Send to output
      stack.setPixels(image, time);
    }

    /**
     * Amplify the photons (assumed to be integer). Apply QE using a Binomial to generate electrons.
     * Optionally apply EM-gain with a Gamma distribution (more electrons).
     *
     * @param image the image
     */
    private void amplify(float[] image) {
      // Quantum efficiency: Model using binomial distribution
      if (qe < 1) {
        for (int i = 0; i < image.length; i++) {
          if (image[i] != 0) {
            image[i] = random.nextBinomial((int) image[i], qe);
          }
        }
      }

      // Apply EM gain
      if (emGain != 0) {
        final boolean tubbsModel = false;
        // See: https://www.andor.com/learning-academy/sensitivity-making-sense-of-sensitivity
        // there is a statistical variation in the overall number of electrons generated from an
        // initial
        // charge packet by the gain register. This uncertainty is quantified by a parameter called
        // "Noise Factor"
        // and detailed theoretical and measured analysis has placed this Noise Factor at a value of
        // √2 (or 1.41)
        // for EMCCD technology.

        // A Stochastic Model for Electron Multiplication Charge-Coupled Devices – From Theory to
        // Practice
        // (PLoS One. 2013; 8(1): e53671)
        // PGN model:
        // - Poisson for photon shot noise
        // - Gamma for EM gain
        // - Normal for read noise
        // EM gain is essentially a repeated loop of input N and get out M where the each N has a
        // probability of
        // being amplified. This has been modelled as a series of Poisson or Binomial trials and
        // then the curve
        // fitted.

        if (tubbsModel) {
          // Tubbs's model
          // Equation 14: is a gamma distribution for electrons created in the register
          for (int i = 0; i < image.length; i++) {
            if (image[i] <= 0) {
              continue;
            }
            final double scale = emGain - 1 + 1 / image[i];
            final double electrons = Math.round(random.nextGamma(image[i], scale) - 1);
            image[i] += electrons;
          }
        } else {
          // Standard gamma distribution
          // This is what is modelled in the Poisson-Gamma-Gaussian likelihood model

          // Since the call random.nextGamma(...) creates a Gamma distribution
          // which pre-calculates factors only using the scale parameter we
          // create a custom gamma distribution where the shape can be set as a property.
          final double shape = 1;
          final CustomGammaDistribution dist =
              new CustomGammaDistribution(random.getRandomGenerator(), shape, emGain,
                  GammaDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);

          for (int i = 0; i < image.length; i++) {
            if (image[i] <= 0) {
              // What about modelling spontaneous electron events?
              continue;
            }
            dist.setShapeUnsafe(image[i]);
            image[i] = Math.round(dist.sample());
          }
        }
      }
    }

    private void erase(float[] data, Spot spot) {
      if (spot.samplePositions != null) {
        psfModel.eraseSample(data, settings.getSize(), settings.getSize(), spot.samplePositions);
      } else {
        psfModel.erase(data, settings.getSize(), settings.getSize(), spot.psf, spot.x0min,
            spot.x0max, spot.x1min, spot.x1max);
      }
    }

    /**
     * Compute the mean and variance for image 1 and image 2 for the region where the spot was
     * inserted. The region is defined by an area around the origin.
     *
     * @param image1 the image 1
     * @param image2 the image 2
     * @param origX the origin X
     * @param origY the origin Y
     * @return [mean1, variance1, mean2, variance2]
     */
    private double[] getStatistics(float[] image1, float[] image2, int origX, int origY) {
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
      while (x0range < width) {
        x0min = limit(x0min - 1);
        x0max = limit(x0max + 1);
        x0range = x0max - x0min + 1;
      }
      int x1range = x1max - x1min + 1;
      while (x1range < width) {
        x1min = limit(x1min - 1);
        x1max = limit(x1max + 1);
        x1range = x1max - x1min + 1;
      }

      final Statistics sum = new Statistics();
      final Statistics sum2 = new Statistics();
      for (int y = 0; y < x1range; y++) {
        // Locate the insert location
        int indexTo = (y + x1min) * settings.getSize() + x0min;
        for (int x = 0; x < x0range; x++) {
          sum.add(image1[indexTo]);
          sum2.add(image2[indexTo]);
          indexTo++;
        }
      }
      return new double[] {sum.getMean(), sum.getVariance(), sum2.getMean(), sum2.getVariance()};
    }
  }

  private int limit(int coord) {
    if (coord < 0) {
      return 0;
    }
    if (coord >= settings.getSize()) {
      return settings.getSize() - 1;
    }
    return coord;
  }

  /**
   * Check if the localisation, or its neighbours, reach the SNR thresholds. The intensity and noise
   * are after EM-gain has been applied.
   *
   * @param localisationSet the localisation set
   * @param intensity the intensity
   * @param noise the noise
   * @return true, if successful
   */
  public boolean badLocalisation(LocalisationModelSet localisationSet, double intensity,
      double noise) {
    // Set the minimum SNR for either a single spot or for a spot next to a brighter neighbour
    double minSnr = settings.getMinSnrT1();
    AtomicInteger counter = removedT1;

    if (localisationSet.hasNeighbour()) {
      final double nextIntensity = getIntensity(localisationSet.getNext());
      final double previousIntensity = getIntensity(localisationSet.getPrevious());

      // Check if either neighbour is above the t1 threshold
      if (Math.max(nextIntensity, previousIntensity) / noise > settings.getMinSnrT1()) {
        // If neighbours are bright then use a more lenient threshold
        minSnr = settings.getMinSnrTN();
        counter = removedTn;
      }
    }

    if (intensity / noise < minSnr) {
      counter.incrementAndGet();
      return true;
    }
    return false;
  }

  private static double getIntensity(LocalisationModelSet localisationSet) {
    if (localisationSet != null) {
      return localisationSet.getData()[4];
    }
    return 0;
  }

  private float[] backgroundPixels;
  private boolean uniformBackground;

  private float[] createBackground(RandomDataGenerator random) {
    float[] pixels2 = null;

    if (settings.getBackground() > 0) {
      if (random == null) {
        random = new RandomDataGenerator(createRandomGenerator());
      }
      pixels2 = new float[backgroundPixels.length];

      // Add Poisson noise.
      if (uniformBackground) {
        final double mean = backgroundPixels[0];

        // Simulate N photons hitting the image. The total photons (N) is
        // the mean for each pixel multiplied by the number of pixels.
        // Note: The number of samples (N) must be Poisson distributed, i.e.
        // the total amount of photons per frame is Poisson noise.
        long samples = random.nextPoisson(mean * backgroundPixels.length);

        final int upper = pixels2.length - 1;
        final UniformIntegerDistribution d =
            new UniformIntegerDistribution(random.getRandomGenerator(), 0, upper);
        while (samples-- > 0) {
          pixels2[d.sample()] += 1;
        }

        //// Alternative is to sample each pixel from a Poisson distribution. This is slow
        // PoissonDistribution dist = new PoissonDistribution(random.getRandomGenerator(),
        //// backgroundPixels[0],
        // PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);
        // for (int i = 0; i < pixels2.length; i++)
        // {
        // pixels2[i] = dist.sample();
        // }
      } else {
        final CustomPoissonDistribution dist =
            new CustomPoissonDistribution(random.getRandomGenerator(), 1,
                PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);
        for (int i = 0; i < pixels2.length; i++) {
          if (backgroundPixels[i] > 0) {
            dist.setMeanUnsafe(backgroundPixels[i]);
            pixels2[i] = dist.sample();
          }
        }
      }
    } else {
      pixels2 = backgroundPixels.clone();
    }

    return pixels2;
  }

  private void createBackgroundPixels() {
    // Cache illumination background
    if (backgroundPixels == null) {
      backgroundPixels = new float[settings.getSize() * settings.getSize()];

      final ImagePlus imp = WindowManager.getImage(settings.getBackgroundImage());
      if (imp != null) {
        // Use an image for the background
        ImageProcessor ip = imp.getProcessor().duplicate().toFloat(0, null);
        ip.setInterpolationMethod(ImageProcessor.BILINEAR);
        ip = ip.resize(settings.getSize(), settings.getSize());
        final float[] data = (float[]) ip.getPixels();
        final double max = MathUtils.maxDefault(0, data);
        if (max != 0) {
          final double scale = settings.getBackground() / max;
          for (int i = 0; i < backgroundPixels.length; i++) {
            // Ignore pixels below zero
            backgroundPixels[i] = (float) (FastMath.max(0, data[i]) * scale);
          }
          return;
        }
      }

      // Use the illumination (this is the fall-back method if the background image has no
      // maximum)
      final SpatialIllumination illumination = createIllumination(settings.getBackground(), 0);
      final double[] xyz = new double[3];
      for (int y = 0, i = 0; y < settings.getSize(); y++) {
        xyz[1] = y - settings.getSize() / 2;
        for (int x = 0, x2 = -settings.getSize() / 2; x < settings.getSize(); x++, x2++, i++) {
          xyz[0] = x2;
          backgroundPixels[i] = (float) illumination.getPhotons(xyz);
        }
      }
    }
  }

  /**
   * Find the first index from the starting index where the localisation matches the time.
   *
   * @param localisations the localisations
   * @param fromIndex start index
   * @param timt time
   * @return the index (or -1)
   */
  private static int findIndexByTime(List<LocalisationModelSet> localisations, int fromIndex,
      int timt) {
    while (fromIndex < localisations.size() && localisations.get(fromIndex).getTime() != timt) {
      fromIndex++;
    }
    return fromIndex >= localisations.size() ? -1 : fromIndex;
  }

  /**
   * Find the last index from the starting index where the localisation matches the time.
   *
   * @param localisations the localisations
   * @param fromIndex start index
   * @param time time
   * @return the index (or -1)
   */
  private static int findLastIndexByTime(List<LocalisationModelSet> localisations, int fromIndex,
      int time) {
    // Check the start point is valid
    if (localisations.get(fromIndex).getTime() != time) {
      fromIndex = findIndexByTime(localisations, 0, time);
      if (fromIndex == -1) {
        return fromIndex;
      }
    }
    while (fromIndex < localisations.size() && localisations.get(fromIndex).getTime() == time) {
      fromIndex++;
    }
    return fromIndex - 1;
  }

  /**
   * Find the first index from the starting index where the localisation matches the id.
   *
   * @param localisations the localisations
   * @param fromIndex start index
   * @param id the id
   * @return the index (or -1)
   */
  private static int findIndexById(List<LocalisationModel> localisations, int fromIndex, int id) {
    while (fromIndex < localisations.size() && localisations.get(fromIndex).getId() != id) {
      fromIndex++;
    }
    return fromIndex >= localisations.size() ? -1 : fromIndex;
  }

  /**
   * Find the last index from the starting index where the localisation matches the Id.
   *
   * @param localisations the localisations
   * @param fromIndex start index
   * @param id the id
   * @return the index (or -1)
   */
  private static int findLastIndexById(List<LocalisationModel> localisations, int fromIndex,
      int id) {
    // Check the start point is valid
    if (localisations.get(fromIndex).getId() != id) {
      fromIndex = findIndexById(localisations, 0, id);
      if (fromIndex == -1) {
        return fromIndex;
      }
    }
    while (fromIndex < localisations.size() && localisations.get(fromIndex).getId() == id) {
      fromIndex++;
    }
    return fromIndex - 1;
  }

  private static List<LocalisationModel>
      toLocalisations(List<LocalisationModelSet> localisationSets) {
    final ArrayList<LocalisationModel> localisations = new ArrayList<>(localisationSets.size());
    for (final LocalisationModelSet s : localisationSets) {
      localisations.add(s.toLocalisation());
    }
    return localisations;
  }

  /**
   * Remove all fluorophores which were not drawn.
   *
   * @param fluorophores the fluorophores
   * @param localisations the localisations
   * @return the filtered list
   */
  private List<? extends FluorophoreSequenceModel> removeFilteredFluorophores(
      List<? extends FluorophoreSequenceModel> fluorophores,
      List<LocalisationModel> localisations) {
    if (fluorophores == null) {
      return null;
    }
    // movingMolecules will be created with an initial capacity to hold all the unique IDs
    final TIntHashSet idSet =
        new TIntHashSet((movingMolecules != null) ? movingMolecules.capacity() : 0);
    for (final LocalisationModel l : localisations) {
      idSet.add(l.getId());
    }
    final List<FluorophoreSequenceModel> newFluorophores = new ArrayList<>(idSet.size());
    for (final FluorophoreSequenceModel f : fluorophores) {
      if (idSet.contains(f.getId())) {
        newFluorophores.add(f);
      }
    }
    return newFluorophores;
  }

  private double showSummary(List<? extends FluorophoreSequenceModel> fluorophores,
      List<LocalisationModel> localisations) {
    IJ.showStatus("Calculating statistics ...");

    createSummaryTable();

    final Statistics[] stats = new Statistics[NAMES.length];
    for (int i = 0; i < stats.length; i++) {
      stats[i] =
          (settings.getShowHistograms() || alwaysRemoveOutliers[i]) ? new StoredDataStatistics()
              : new Statistics();
    }

    // Find the largest timepoint
    final ImagePlus outputImp = WindowManager.getImage(benchmarkImageId);
    int frameCount;
    if (outputImp == null) {
      sortLocalisationsByTime(localisations);
      frameCount = localisations.get(localisations.size() - 1).getTime();
    } else {
      frameCount = outputImp.getStackSize();
    }
    final int[] countHistogram = new int[frameCount + 1];

    // Use the localisations that were drawn to create the sampled on/off times
    rebuildNeighbours(localisations);

    // Assume that there is at least one localisation
    final LocalisationModel first = localisations.get(0);
    int currentId = first.getId(); // The current localisation
    int lastT = first.getTime(); // The last time this localisation was on
    int blinks = 0; // Number of blinks
    int currentT = 0; // On-time of current pulse
    double signal = 0;
    final double centreOffset = settings.getSize() * 0.5;
    // Used to convert the sampled times in frames into seconds
    final double framesPerSecond = 1000.0 / settings.getExposureTime();
    // final double gain = new CreateDataSettingsHelper(settings).getTotalGainSafe();
    for (final LocalisationModel l : localisations) {
      final double[] data = l.getData();
      if (data == null) {
        throw new IllegalStateException("No localisation data. This should not happen!");
      }
      final double noise = data[1];
      final double sx = data[2];
      final double sy = data[3];
      final double intensityInPhotons = data[4];
      // Q. What if the noise is zero, i.e. no background photon / read noise?
      // Just ignore it at current. This is only an approximation to the SNR estimate
      // if this is not a Gaussian spot.
      final double snr =
          Gaussian2DPeakResultHelper.getMeanSignalUsingP05(intensityInPhotons, sx, sy) / noise;
      stats[SIGNAL].add(intensityInPhotons);
      stats[NOISE].add(noise);
      if (noise != 0) {
        stats[SNR].add(snr);
      }
      // Average intensity only from continuous spots.
      // The continuous flag is for spots that have all the simulation steps continuously on.
      // Try using the neighbour pointers instead to get the 'sampled' continuous spots.
      // if (l.isContinuous())
      if (l.getNext() != null && l.getPrevious() != null) {
        stats[SIGNAL_CONTINUOUS].add(intensityInPhotons);
        if (noise != 0) {
          stats[SNR_CONTINUOUS].add(snr);
        }
      }

      final int id = l.getId();
      // Check if this a new fluorophore
      if (currentId != id) {
        // Add previous fluorophore
        stats[SAMPLED_BLINKS].add(blinks);
        stats[SAMPLED_T_ON].add(currentT / framesPerSecond);
        stats[TOTAL_SIGNAL].add(signal);

        // Reset
        blinks = 0;
        currentT = 1;
        currentId = id;
        signal = intensityInPhotons;
      } else {
        signal += intensityInPhotons;
        // Check if the current fluorophore pulse is broken (i.e. a blink)
        if (l.getTime() - 1 > lastT) {
          blinks++;
          stats[SAMPLED_T_ON].add(currentT / framesPerSecond);
          currentT = 1;
          stats[SAMPLED_T_OFF].add(((l.getTime() - 1) - lastT) / framesPerSecond);
        } else {
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
    for (int t = 1; t < countHistogram.length; t++) {
      stats[SAMPLES].add(countHistogram[t]);
    }

    if (fluorophores != null) {
      for (final FluorophoreSequenceModel f : fluorophores) {
        stats[BLINKS].add(f.getNumberOfBlinks());
        // On-time
        for (final double t : f.getOnTimes()) {
          stats[T_ON].add(t);
        }
        // Off-time
        for (final double t : f.getOffTimes()) {
          stats[T_OFF].add(t);
        }
      }
    } else {
      // show no blinks
      stats[BLINKS].add(0);
      stats[T_ON].add(1);
      // stats[T_OFF].add(0);
    }

    if (results != null) {
      // Convert depth-of-field to pixels
      final double depth = settings.getDepthOfField() / settings.getPixelPitch();

      try {
        // Get widths
        final WidthResultProcedure wp = new WidthResultProcedure(results, DistanceUnit.PIXEL);
        wp.getW();
        stats[WIDTH].add(wp.wx);
      } catch (final DataException ex) {
        ImageJUtils.log("Unable to compute width: " + ex.getMessage());
      }

      try {
        // Get z depth
        final StandardResultProcedure sp = new StandardResultProcedure(results, DistanceUnit.PIXEL);
        sp.getXyz();

        // Get precision
        final PrecisionResultProcedure pp = new PrecisionResultProcedure(results);
        pp.getPrecision();
        stats[PRECISION].add(pp.precisions);
        for (int i = 0; i < pp.size(); i++) {
          if (Math.abs(sp.z[i]) < depth) {
            stats[PRECISION_IN_FOCUS].add(pp.precisions[i]);
          }
        }
      } catch (final DataException ex) {
        ImageJUtils.log("Unable to compute LSE precision: " + ex.getMessage());
      }

      // Compute density per frame. Multithread for speed
      if (settings.getDensityRadius() > 0) {
        IJ.showStatus("Calculating density ...");
        final ExecutorService threadPool = Executors.newFixedThreadPool(Prefs.getThreads());
        final List<Future<?>> futures = new LinkedList<>();
        final TFloatArrayList coordsX = new TFloatArrayList();
        final TFloatArrayList coordsY = new TFloatArrayList();
        final Statistics densityStats = stats[DENSITY];
        final float radius = (float) (settings.getDensityRadius() * getHwhm());
        final Rectangle bounds = results.getBounds();
        currentIndex = 0;
        finalIndex = results.getLastFrame();
        // Store the density for each result.
        final int[] allDensity = new int[results.size()];
        final FrameCounter counter = results.newFrameCounter();
        results.forEach((PeakResultProcedure) result -> {
          if (counter.advance(result.getFrame())) {
            counter.increment(runDensityCalculation(threadPool, futures, coordsX, coordsY,
                densityStats, radius, bounds, allDensity, counter.getCount()));
          }
          coordsX.add(result.getXPosition());
          coordsY.add(result.getYPosition());
        });
        runDensityCalculation(threadPool, futures, coordsX, coordsY, densityStats, radius, bounds,
            allDensity, counter.getCount());
        ConcurrencyUtils.waitForCompletionUnchecked(futures);
        threadPool.shutdownNow();
        IJ.showProgress(1);

        // Split results into singles (density = 0) and clustered (density > 0)
        final MemoryPeakResults singles = copyMemoryPeakResults("No Density");
        final MemoryPeakResults clustered = copyMemoryPeakResults("Density");

        counter.reset();
        results.forEach((PeakResultProcedure) result -> {
          final int density = allDensity[counter.getAndIncrement()];
          result.setOrigValue(density);
          if (density == 0) {
            singles.add(result);
          } else {
            clustered.add(result);
          }
        });
      }
    }

    final StringBuilder sb = new StringBuilder();
    sb.append(datasetNumber).append('\t');
    if (settings.getCameraType() == CameraType.SCMOS) {
      sb.append("sCMOS (").append(settings.getCameraModelName()).append(") ");
      final Rectangle bounds = cameraModel.getBounds();
      sb.append(" ").append(bounds.x).append(",").append(bounds.y);
      final int size = settings.getSize();
      sb.append(" ").append(size).append("x").append(size);
    } else if (CalibrationProtosHelper.isCcdCameraType(settings.getCameraType())) {
      sb.append(CalibrationProtosHelper.getName(settings.getCameraType()));
      final int size = settings.getSize();
      sb.append(" ").append(size).append("x").append(size);
      if (settings.getCameraType() == CameraType.EMCCD) {
        sb.append(" EM=").append(settings.getEmGain());
      }
      sb.append(" CG=").append(settings.getCameraGain());
      sb.append(" RN=").append(settings.getReadNoise());
      sb.append(" B=").append(settings.getBias());
    } else {
      throw new IllegalStateException();
    }
    sb.append(" QE=").append(settings.getQuantumEfficiency()).append('\t');
    sb.append(settings.getPsfModel());
    if (psfModelType == PSF_MODEL_IMAGE) {
      sb.append(" Image").append(settings.getPsfImageName());
    } else if (psfModelType == PSF_MODEL_ASTIGMATISM) {
      sb.append(" model=").append(settings.getAstigmatismModel());
    } else {
      sb.append(" DoF=").append(MathUtils.rounded(settings.getDepthOfFocus()));
      if (settings.getEnterWidth()) {
        sb.append(" SD=").append(MathUtils.rounded(settings.getPsfSd()));
      } else {
        sb.append(" λ=").append(MathUtils.rounded(settings.getWavelength()));
        sb.append(" NA=").append(MathUtils.rounded(settings.getNumericalAperture()));
      }
    }
    sb.append('\t');
    sb.append((fluorophores == null) ? localisations.size() : fluorophores.size()).append('\t');
    sb.append(stats[SAMPLED_BLINKS].getN() + (int) stats[SAMPLED_BLINKS].getSum()).append('\t');
    sb.append(localisations.size()).append('\t');
    sb.append(frameCount).append('\t');
    sb.append(MathUtils.rounded(areaInUm)).append('\t');
    sb.append(MathUtils.rounded(localisations.size() / (areaInUm * frameCount), 4)).append('\t');
    sb.append(MathUtils.rounded(getHwhm(), 4)).append('\t');
    double sd = getPsfSd();
    sb.append(MathUtils.rounded(sd, 4)).append('\t');
    sd *= settings.getPixelPitch();
    final double sa = PsfCalculator.squarePixelAdjustment(sd, settings.getPixelPitch())
        / settings.getPixelPitch();
    sb.append(MathUtils.rounded(sa, 4)).append('\t');
    // Width not valid for the Image PSF.
    // Q. Is this true? We can approximate the FHWM for a spot-like image PSF.
    final int nStats = (psfModelType == PSF_MODEL_IMAGE) ? stats.length - 1 : stats.length;
    for (int i = 0; i < nStats; i++) {
      final double centre = (alwaysRemoveOutliers[i])
          ? ((StoredDataStatistics) stats[i]).getStatistics().getPercentile(50)
          : stats[i].getMean();
      sb.append(MathUtils.rounded(centre, 4)).append('\t');
    }
    if (java.awt.GraphicsEnvironment.isHeadless()) {
      IJ.log(sb.toString());
      return stats[SIGNAL].getMean();
    }
    summaryTable.append(sb.toString());

    // Show histograms
    if (settings.getShowHistograms()) {
      IJ.showStatus("Calculating histograms ...");
      final boolean[] chosenHistograms = getChoosenHistograms();

      final WindowOrganiser wo = new WindowOrganiser();

      final HistogramPlotBuilder builder = new HistogramPlotBuilder(TITLE);
      for (int i = 0; i < NAMES.length; i++) {
        if (chosenHistograms[i]) {
          builder.setData((StoredDataStatistics) stats[i]).setName(NAMES[i])
              .setIntegerBins(integerDisplay[i])
              .setRemoveOutliersOption(
                  (settings.getRemoveOutliers() || alwaysRemoveOutliers[i]) ? 2 : 0)
              .setNumberOfBins(settings.getHistogramBins()).show(wo);
        }
      }

      wo.tile();
    }
    IJ.showStatus("");
    return stats[SIGNAL].getMean();
  }

  private int runDensityCalculation(ExecutorService threadPool, List<Future<?>> futures,
      final TFloatArrayList coordsX, final TFloatArrayList coordsY, final Statistics densityStats,
      final float radius, final Rectangle bounds, final int[] allDensity, final int allIndex) {
    final int size = coordsX.size();
    final float[] xCoords = coordsX.toArray();
    final float[] yCoords = coordsY.toArray();
    coordsX.resetQuick();
    coordsY.resetQuick();
    futures.add(threadPool.submit(() -> {
      incrementProgress();
      final DensityManager dm = new DensityManager(xCoords, yCoords, bounds);
      final int[] density = dm.calculateDensity(radius, true);
      addDensity(densityStats, density);

      // Store the density for each result. This does not need to be synchronised
      // since the indices in different threads are unique.
      for (int i = 0, index = allIndex; i < density.length; i++, index++) {
        allDensity[index] = density[i];
      }
    }));
    return size;
  }

  private int currentIndex;
  private int finalIndex;

  private synchronized void incrementProgress() {
    IJ.showProgress(currentIndex, finalIndex);
  }

  private static synchronized void addDensity(Statistics stats, int[] density) {
    stats.add(density);
  }

  /**
   * Copy all the settings from the results into a new results set labelled with the name suffix.
   *
   * @param nameSuffix the name suffix
   * @return The new results set
   */
  private MemoryPeakResults copyMemoryPeakResults(String nameSuffix) {
    final MemoryPeakResults newResults = new MemoryPeakResults();
    newResults.copySettings(this.results);
    newResults.setName(newResults.getSource().getName() + " (" + TITLE + " " + nameSuffix + ")");
    newResults.setSortAfterEnd(true);
    newResults.begin();
    MemoryPeakResults.addResults(newResults);
    return newResults;
  }

  private boolean[] getChoosenHistograms() {
    if (settings.getChooseHistograms()) {
      return displayHistograms;
    }

    final boolean[] all = new boolean[displayHistograms.length];
    for (int i = 0; i < all.length; i++) {
      all[i] = true;
    }
    return all;
  }

  private static void sortLocalisationsByIdThenTime(List<LocalisationModel> localisations) {
    Collections.sort(localisations, new Comparator<LocalisationModel>() {
      @Override
      public int compare(LocalisationModel o1, LocalisationModel o2) {
        // Order by ID then time
        if (o1.getId() < o2.getId()) {
          return -1;
        }
        if (o1.getId() > o2.getId()) {
          return 1;
        }
        if (o1.getTime() < o2.getTime()) {
          return -1;
        }
        if (o1.getTime() > o2.getTime()) {
          return 1;
        }
        return 0;
      }
    });
  }

  private static void sortLocalisationsByTime(List<LocalisationModel> localisations) {
    Collections.sort(localisations, (o1, o2) -> Integer.compare(o1.getTime(), o2.getTime()));
  }

  /**
   * Sort by id then time, then rebuild the neighbour pointers.
   *
   * @param localisations the localisations
   */
  private static void rebuildNeighbours(List<LocalisationModel> localisations) {
    sortLocalisationsByIdThenTime(localisations);

    int id = 0;
    int time = 0;
    LocalisationModel previous = null;
    for (final LocalisationModel l : localisations) {
      if (l.getId() != id || l.getTime() > time + 1) {
        // New spot so no previous neighbour OR
        // Discontinuous time so no previous neighbour
        previous = null;
      }

      l.setPrevious(previous);
      l.setNext(null);

      id = l.getId();
      time = l.getTime();
      previous = l;
    }
  }

  private static void createSummaryTable() {
    if (java.awt.GraphicsEnvironment.isHeadless()) {
      if (header == null) {
        header = createHeader();
        IJ.log(header);
      }
    } else if (summaryTable == null || !summaryTable.isVisible()) {
      summaryTable = new TextWindow("Data Summary", createHeader(), "", 800, 300);
      summaryTable.setVisible(true);
    }
  }

  private static String createHeader() {
    final StringBuilder sb = new StringBuilder(
        "Dataset\tCamera\tPSF\tMolecules\tPulses\tLocalisations\tnFrames\tArea (um^2)\t"
            + "Density (mol/um^2)\tHWHM\tS\tSa");
    for (int i = 0; i < NAMES.length; i++) {
      sb.append('\t').append(NAMES[i]);
    }
    // if (alwaysRemoveOutliers[i])
    // sb.append("*");
    return sb.toString();
  }

  /**
   * Save the image to a TIFF file.
   *
   * @param imp the imp
   */
  private void saveImage(ImagePlus imp) {
    if (!settings.getSaveImage()) {
      return;
    }
    final String[] path = ImageJUtils.decodePath(settings.getImageFilename());
    final OpenDialog chooser = new OpenDialog("Image_File", path[0], path[1]);
    if (chooser.getFileName() != null) {
      settings.setImageFilename(chooser.getDirectory() + chooser.getFileName());
      settings.setImageFilename(FileUtils.replaceExtension(settings.getImageFilename(), "tiff"));

      final FileSaver fs = new FileSaver(imp);
      boolean ok;
      if (imp.getStackSize() > 1) {
        ok = fs.saveAsTiffStack(settings.getImageFilename());
      } else {
        ok = fs.saveAsTiff(settings.getImageFilename());
        // The following call throws a NoSuchMethodError.
        // ok = IJ.saveAsTiff(imp, settings.imageFilename);
      }

      if (!ok) {
        IJ.log("Failed to save image to file: " + settings.getImageFilename());
      }
    }
  }

  private void saveImageResults(MemoryPeakResults results) {
    if (!settings.getSaveImageResults()) {
      return;
    }
    final String[] path = ImageJUtils.decodePath(settings.getImageResultsFilename());
    final OpenDialog chooser = new OpenDialog("Image_Results_File", path[0], path[1]);
    if (chooser.getFileName() != null) {
      settings.setImageResultsFilename(chooser.getDirectory() + chooser.getFileName());
      settings.setImageResultsFilename(
          FileUtils.replaceExtension(settings.getImageResultsFilename(), "xls"));

      final TextFilePeakResults r =
          new TextFilePeakResults(settings.getImageResultsFilename(), false);
      r.copySettings(results);
      r.begin();
      r.addAll(results.toArray());
      r.end();
    }
  }

  /**
   * Create a set of results that represent the molecule continuous on-times (pulses).
   *
   * @param localisations the localisations
   * @param results the results
   */
  @SuppressWarnings("null")
  private void savePulses(List<LocalisationModel> localisations, MemoryPeakResults results) {
    sortLocalisationsByIdThenTime(localisations);

    final MemoryPeakResults traceResults = copyMemoryPeakResults("Pulses");
    LocalisationModel start = null;
    int currentId = -1;
    int count = 0;
    float[] params = Gaussian2DPeakResultHelper.createTwoAxisParams(0, 0, 0, 0, 0, 0, 0);
    final int isx = Gaussian2DPeakResultHelper.INDEX_SX;
    final int isy = Gaussian2DPeakResultHelper.INDEX_SY;
    double noise = 0;
    int lastT = -1;
    for (final LocalisationModel localisation : localisations) {
      if (currentId != localisation.getId() || lastT + 1 != localisation.getTime()) {
        if (count > 0) {
          params[PeakResult.BACKGROUND] /= count;
          params[PeakResult.X] /= count;
          params[PeakResult.Y] /= count;
          params[isx] /= count;
          params[isy] /= count;

          final ExtendedPeakResult p = new ExtendedPeakResult(start.getTime(),
              (int) Math.round(start.getX()), (int) Math.round(start.getY()), 0, 0,
              (float) (Math.sqrt(noise)), 0, params, null, lastT, currentId);
          // if (p.getPrecision(107, 1) > 2000)
          // {
          // System.out.printf("Weird precision = %g (%d)\n", p.getPrecision(107, 1), n);
          // }
          traceResults.add(p);
        }
        start = localisation;
        currentId = localisation.getId();
        count = 0;
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
      count++;
      lastT = localisation.getTime();
    }

    // Final pulse
    if (count > 0) {
      params[PeakResult.BACKGROUND] /= count;
      params[PeakResult.X] /= count;
      params[PeakResult.Y] /= count;
      params[isx] /= count;
      params[isy] /= count;

      traceResults.add(new ExtendedPeakResult(start.getTime(), (int) Math.round(start.getX()),
          (int) Math.round(start.getY()), 0, 0, (float) (Math.sqrt(noise)), 0, params, null, lastT,
          currentId));
    }

    traceResults.end();
  }

  private void saveFixedAndMoving(MemoryPeakResults results) {
    if (simpleMode || benchmarkMode || spotMode) {
      return;
    }
    if (settings.getDiffusionRate() <= 0 || settings.getFixedFraction() >= 1) {
      return;
    }

    final MemoryPeakResults fixedResults = copyMemoryPeakResults("Fixed");
    final MemoryPeakResults movingResults = copyMemoryPeakResults("Moving");

    final PeakResult[] peakResults = results.toArray();
    // Sort using the ID
    Arrays.sort(peakResults, new Comparator<PeakResult>() {
      @Override
      public int compare(PeakResult o1, PeakResult o2) {
        return o1.getId() - o2.getId();
      }
    });

    MemoryPeakResults currentResults = movingResults;
    final FrameCounter counter = new FrameCounter(-1);
    for (final PeakResult p : peakResults) {
      if (counter.advance(p.getId())) {
        currentResults = (movingMolecules.contains(p.getId())) ? movingResults : fixedResults;
      }
      currentResults.add(p);
    }

    movingResults.end();
    fixedResults.end();
  }

  @SuppressWarnings("null")
  private void saveCompoundMolecules(MemoryPeakResults results) {
    if (idToCompound == null) {
      return;
    }

    final MemoryPeakResults[] set = new MemoryPeakResults[compoundNames.size()];
    for (int i = 0; i < set.length; i++) {
      set[i] = copyMemoryPeakResults("Compound " + (i + 1) + ", " + compoundNames.get(i));
    }

    final PeakResult[] peakResults = results.toArray();
    // Sort using the ID
    Arrays.sort(peakResults, new Comparator<PeakResult>() {
      @Override
      public int compare(PeakResult o1, PeakResult o2) {
        return o1.getId() - o2.getId();
      }
    });

    MemoryPeakResults currentResults = null;
    final FrameCounter counter = new FrameCounter(-1);
    for (final PeakResult p : peakResults) {
      if (counter.advance(p.getId())) {
        currentResults = set[idToCompound.get(p.getId())];
      }
      currentResults.add(p);
    }

    for (int i = 0; i < set.length; i++) {
      set[i].end();
    }
  }

  /**
   * Update the fluorophores relative coordinates to absolute.
   *
   * @param molecules the molecules
   */
  @SuppressWarnings("unused")
  private static void convertRelativeToAbsolute(List<CompoundMoleculeModel> molecules) {
    for (final CompoundMoleculeModel c : molecules) {
      final double[] xyz = c.getCoordinates();
      for (int n = c.getSize(); n-- > 0;) {
        final MoleculeModel m = c.getMolecule(n);
        final double[] xyz2 = m.getCoordinates();
        for (int i = 0; i < 3; i++) {
          xyz2[i] += xyz[i];
        }
      }
    }
  }

  /**
   * Save the fluorophores to a text file.
   *
   * @param fluorophores the fluorophores
   */
  private void saveFluorophores(List<? extends FluorophoreSequenceModel> fluorophores) {
    if (!settings.getSaveFluorophores() || fluorophores == null) {
      return;
    }

    final String[] path = ImageJUtils.decodePath(settings.getFluorophoresFilename());
    final OpenDialog chooser = new OpenDialog("Fluorophores_File", path[0], path[1]);
    if (chooser.getFileName() != null) {
      settings.setFluorophoresFilename(chooser.getDirectory() + chooser.getFileName());
      settings.setFluorophoresFilename(
          FileUtils.replaceExtension(settings.getFluorophoresFilename(), "xls"));

      try (BufferedWriter output =
          Files.newBufferedWriter(Paths.get(settings.getFluorophoresFilename()))) {
        output.write(createResultsFileHeader());
        output.write("#Id\tn-Blinks\tStart\tStop\t...");
        output.newLine();
        for (int id = 1; id <= fluorophores.size(); id++) {
          final FluorophoreSequenceModel f = fluorophores.get(id - 1);
          final StringBuilder sb = new StringBuilder();
          sb.append(f.getId()).append('\t');
          sb.append(f.getNumberOfBlinks()).append('\t');
          for (final double[] burst : f.getBurstSequence()) {
            sb.append(MathUtils.rounded(burst[0], 3)).append('\t')
                .append(MathUtils.rounded(burst[1], 3)).append('\t');
          }
          output.write(sb.toString());
          output.newLine();
        }
      } catch (final Exception ex) {
        // Q. Add better handling of errors?
        ex.printStackTrace();
        IJ.log("Failed to save fluorophores to file: " + settings.getFluorophoresFilename());
      }
    }
  }

  /**
   * Save the localisations to a text file.
   *
   * @param localisations the localisations
   */
  private void saveLocalisations(List<LocalisationModel> localisations) {
    if (!settings.getSaveLocalisations()) {
      return;
    }

    sortLocalisationsByTime(localisations);

    // Collections.sort(localisations, new Comparator<LocalisationModel>(){
    //
    // public int compare(LocalisationModel o1, LocalisationModel o2)
    // {
    // int cellx1 = (int)(o1.getX() / settings.cellSize);
    // int cellx2 = (int)(o2.getX() / settings.cellSize);
    // int result = cellx2 - cellx1;
    // if (result != 0)
    // return result;
    // int celly1 = (int)(o1.getY() / settings.cellSize);
    // int celly2 = (int)(o2.getY() / settings.cellSize);
    // result = celly2 - celly1;
    // if (result != 0)
    // return result;
    // return (o1.getZ() == o2.getZ()) ? 0 : (o1.getZ() == 0) ? -1 : 1;
    // }});

    final LoadLocalisationsSettings.Builder settings =
        SettingsManager.readLoadLocalisationsSettings(0).toBuilder();

    final String[] path = ImageJUtils.decodePath(settings.getLocalisationsFilename());
    final OpenDialog chooser = new OpenDialog("Localisations_File", path[0], path[1]);
    if (chooser.getFileName() != null) {
      settings.setLocalisationsFilename(
          FileUtils.replaceExtension(chooser.getDirectory() + chooser.getFileName(), "xls"));
      SettingsManager.writeSettings(settings.build());

      try (BufferedWriter output =
          Files.newBufferedWriter(Paths.get(settings.getLocalisationsFilename()))) {
        output.write(createResultsFileHeader());
        output.write("#T\tId\tX\tY\tZ\tIntensity");
        output.newLine();
        for (final LocalisationModel l : localisations) {
          final StringBuilder sb = new StringBuilder();
          sb.append(l.getTime()).append('\t');
          sb.append(l.getId()).append('\t');
          sb.append(l.getX()).append('\t');
          sb.append(l.getY()).append('\t');
          sb.append(l.getZ()).append('\t');
          sb.append(l.getIntensity());
          output.write(sb.toString());
          output.newLine();
        }
      } catch (final Exception ex) {
        // Q. Add better handling of errors?
        ex.printStackTrace();
        IJ.log("Failed to save localisations to file: " + settings.getLocalisationsFilename());
      }
    }
  }

  private String createResultsFileHeader() {
    if (resultsFileHeader == null) {
      final String[] backgroundImages = createBackgroundImageList();

      final StringBuilder sb = new StringBuilder();
      sb.append("# ").append(TITLE).append(" Parameters:\n");
      addHeaderLine(sb, "Pixel_pitch (nm)", settings.getPixelPitch());
      addHeaderLine(sb, "Size", settings.getSize());
      if (!benchmarkMode) {
        addHeaderLine(sb, "Depth", settings.getDepth());
        addHeaderLine(sb, "Fixed depth", settings.getFixedDepth());
      }
      if (!(simpleMode || benchmarkMode)) {
        if (!trackMode) {
          addHeaderLine(sb, "Seconds", settings.getSeconds());
        }
        addHeaderLine(sb, "Exposure_time", settings.getExposureTime());
        addHeaderLine(sb, "Steps_per_second", settings.getStepsPerSecond());
        if (!trackMode) {
          addHeaderLine(sb, "Illumination", settings.getIllumination());
          addHeaderLine(sb, "Pulse_interval", settings.getPulseInterval());
          addHeaderLine(sb, "Pulse_ratio", settings.getPulseRatio());
        }
        if (backgroundImages != null) {
          addHeaderLine(sb, "Background_image", settings.getBackgroundImage());
        }
      }
      addHeaderLine(sb, "Background", settings.getBackground());
      addCameraOptionsHeader(sb);
      addHeaderLine(sb, "PSF_model", settings.getPsfModel());
      if (psfModelType == PSF_MODEL_IMAGE) {
        addHeaderLine(sb, "PSF_image", settings.getPsfImageName());
      } else if (psfModelType == PSF_MODEL_ASTIGMATISM) {
        addHeaderLine(sb, "Astigmatism_model", settings.getAstigmatismModel());
        // Q. Should the actual model be appended?
        // addHeaderLine(sb, "Astigmatism_model parameters", SettingsManager.toJSON(getTheModel()));
      } else {
        addHeaderLine(sb, "Depth-of-focus (nm)", settings.getDepthOfFocus());
        if (settings.getEnterWidth()) {
          addHeaderLine(sb, "PSF_SD", settings.getPsfSd());
        } else {
          addHeaderLine(sb, "Wavelength (nm)", settings.getWavelength());
          addHeaderLine(sb, "Numerical_aperture", settings.getNumericalAperture());
        }
      }
      if (!(benchmarkMode || spotMode)) {
        addHeaderLine(sb, "Distribution", settings.getDistribution());
        if (settings.getDistribution().equals(DISTRIBUTION[MASK])) {
          addHeaderLine(sb, "Distribution_mask", settings.getDistributionMask());
        } else if (settings.getDistribution().equals(DISTRIBUTION[GRID])) {
          addHeaderLine(sb, "Cell_size", settings.getCellSize());
          addHeaderLine(sb, "p-binary", settings.getProbabilityBinary());
          addHeaderLine(sb, "Min_binary_distance (nm)", settings.getMinBinaryDistance());
          addHeaderLine(sb, "Max_binary_distance (nm)", settings.getMaxBinaryDistance());
        }
      }
      addHeaderLine(sb, "Particles", settings.getParticles());
      if (benchmarkMode) {
        addHeaderLine(sb, "X_position", settings.getXPosition());
        addHeaderLine(sb, "Y_position", settings.getYPosition());
        addHeaderLine(sb, "Z_position", settings.getZPosition());
        addHeaderLine(sb, "Min_photons", settings.getPhotonsPerSecond());
        addHeaderLine(sb, "Max_photons", settings.getPhotonsPerSecondMaximum());
      } else if (simpleMode) {
        addHeaderLine(sb, "Density (um^-2)", settings.getDensity());
        addHeaderLine(sb, "Min_photons", settings.getPhotonsPerSecond());
        addHeaderLine(sb, "Max_photons", settings.getPhotonsPerSecondMaximum());
      } else if (spotMode) {
        addHeaderLine(sb, "Min_photons", settings.getPhotonsPerSecond());
        addHeaderLine(sb, "Max_photons", settings.getPhotonsPerSecondMaximum());
      } else {
        addHeaderLine(sb, "Diffusion_rate", settings.getDiffusionRate());
        addHeaderLine(sb, "Diffusion_type", settings.getDiffusionType());
        addHeaderLine(sb, "Fixed_fraction", settings.getFixedFraction());
        if (settings.getCompoundMolecules()) {
          addHeaderLine(sb, "Compound_molecules",
              settings.getCompoundText().replaceAll("\n *", ""));
          addHeaderLine(sb, "Enable_2D_diffusion", settings.getDiffuse2D());
          addHeaderLine(sb, "Rotate_initial_orientation", settings.getRotateInitialOrientation());
          addHeaderLine(sb, "Rotate_during_simulation", settings.getRotateDuringSimulation());
          addHeaderLine(sb, "Enable_2D_rotation", settings.getRotate2D());
        }
        addHeaderLine(sb, "Confinement", settings.getConfinement());
        if (settings.getConfinement().equals(CONFINEMENT[CONFINEMENT_SPHERE])) {
          addHeaderLine(sb, "Confinement_radius", settings.getConfinementRadius());
        } else if (settings.getConfinement().equals(CONFINEMENT[CONFINEMENT_MASK])) {
          addHeaderLine(sb, "Confinement_mask", settings.getConfinementMask());
        }
        addHeaderLine(sb, "Photon", settings.getPhotonsPerSecond());
        addHeaderLine(sb, "Photon_distribution", settings.getPhotonDistribution());
        if (PHOTON_DISTRIBUTION[PHOTON_CUSTOM].equals(settings.getPhotonDistribution())) {
          addHeaderLine(sb, "Photon_distribution_file", settings.getPhotonDistributionFile());
        } else if (PHOTON_DISTRIBUTION[PHOTON_UNIFORM].equals(settings.getPhotonDistribution())) {
          addHeaderLine(sb, "Photon_max", settings.getPhotonsPerSecondMaximum());
        } else if (PHOTON_DISTRIBUTION[PHOTON_GAMMA].equals(settings.getPhotonDistribution())) {
          addHeaderLine(sb, "Photon_shape", settings.getPhotonShape());
        } else if (PHOTON_DISTRIBUTION[PHOTON_CORRELATED]
            .equals(settings.getPhotonDistribution())) {
          addHeaderLine(sb, "Correlation", settings.getCorrelation());
        }
        addHeaderLine(sb, "On_time", settings.getTOn());
        if (!trackMode) {
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

  private static void addHeaderLine(StringBuilder sb, String name, Object object) {
    sb.append(String.format("# %-20s = %s%n", name, object.toString()));
  }

  /**
   * Show a dialog allowing the parameters for a simple/benchmark simulation to be performed.
   *
   * @return True if the parameters were collected
   */
  private boolean showSimpleDialog() {
    ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

    settings = SettingsManager.readCreateDataSettings(0).toBuilder();

    // Image size
    gd.addMessage("--- Image Size ---");
    gd.addNumericField("Pixel_pitch (nm)", settings.getPixelPitch(), 2);
    gd.addNumericField("Size (px)", settings.getSize(), 0);
    if (!benchmarkMode) {
      gd.addNumericField("Depth (nm)", settings.getDepth(), 0);
      gd.addCheckbox("Fixed_depth", settings.getFixedDepth());
    }

    // Noise model
    gd.addMessage("--- Noise Model ---");
    if (extraOptions) {
      gd.addCheckbox("No_poisson_noise", !settings.getPoissonNoise());
    }
    gd.addNumericField("Background (photons)", settings.getBackground(), 2);

    addCameraOptions(gd);

    addPsfOptions(gd);

    gd.addMessage("--- Fluorophores ---");
    // Do not allow grid or mask distribution
    if (simpleMode) {
      // Allow mask but not the grid
      gd.addChoice("Distribution", Arrays.copyOf(DISTRIBUTION, DISTRIBUTION.length - 1),
          settings.getDistribution());
      gd.addCheckbox("Sample_per_frame", settings.getSamplePerFrame());
    }
    gd.addNumericField("Particles", settings.getParticles(), 0);
    if (simpleMode) {
      gd.addNumericField("Density (um^-2)", settings.getDensity(), 2);
    } else if (benchmarkMode) {
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
    if (simpleMode) {
      gd.addSlider("Density_radius (N x HWHM)", 0, 4.5, settings.getDensityRadius());
    }
    gd.addNumericField("Depth-of-field (nm)", settings.getDepthOfField(), 0);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.setPixelPitch(Math.abs(gd.getNextNumber()));
    settings.setSize(Math.abs((int) gd.getNextNumber()));
    if (!benchmarkMode) {
      // Allow negative depth
      settings.setDepth(gd.getNextNumber());
      settings.setFixedDepth(gd.getNextBoolean());
    }

    if (extraOptions) {
      settings.setPoissonNoise(!gd.getNextBoolean());
      poissonNoise = settings.getPoissonNoise();
    }
    settings.setBackground(Math.abs(gd.getNextNumber()));
    settings.setCameraType(SettingsManager.getCameraTypeValues()[gd.getNextChoiceIndex()]);

    settings.setPsfModel(gd.getNextChoice());

    if (simpleMode) {
      settings.setDistribution(gd.getNextChoice());
      settings.setSamplePerFrame(gd.getNextBoolean());
    }
    settings.setParticles(Math.abs((int) gd.getNextNumber()));
    if (simpleMode) {
      settings.setDensity(Math.abs(gd.getNextNumber()));
    } else if (benchmarkMode) {
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
    if (simpleMode) {
      settings.setDensityRadius((float) gd.getNextNumber());
    }
    settings.setDepthOfField((float) Math.abs(gd.getNextNumber()));

    gd.collectOptions();

    // Save before validation so that the current values are preserved.
    SettingsManager.writeSettings(settings.build());

    if (gd.invalidNumber()) {
      return false;
    }

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Pixel Pitch", settings.getPixelPitch());
      ParameterUtils.isAboveZero("Size", settings.getSize());
      if (!benchmarkMode && !settings.getFixedDepth()) {
        ParameterUtils.isPositive("Depth", settings.getDepth());
      }
      ParameterUtils.isPositive("Background", settings.getBackground());
      ParameterUtils.isAboveZero("Particles", settings.getParticles());
      if (simpleMode) {
        ParameterUtils.isAboveZero("Density", settings.getDensity());
      }
      ParameterUtils.isAboveZero("Min Photons", settings.getPhotonsPerSecond());
      if (settings.getPhotonsPerSecondMaximum() < settings.getPhotonsPerSecond()) {
        settings.setPhotonsPerSecondMaximum(settings.getPhotonsPerSecond());
      }
      ParameterUtils.isPositive("Histogram bins", settings.getHistogramBins());
      if (simpleMode) {
        ParameterUtils.isPositive("Density radius", settings.getDensityRadius());
      }

      validateCameraOptions();
      validatePsfOptions();
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    if (!benchmarkMode && settings.getDistribution().equals(DISTRIBUTION[MASK])) {
      final String[] maskImages = createDistributionImageList();
      if (maskImages != null) {
        gd = new ExtendedGenericDialog(TITLE);
        gd.addMessage("Select the mask image for the distribution");
        gd.addChoice("Distribution_mask", maskImages, settings.getDistributionMask());
        if (maskListContainsStacks) {
          gd.addNumericField("Distribution_slice_depth (nm)",
              settings.getDistributionMaskSliceDepth(), 0);
        }
        gd.showDialog();
        if (gd.wasCanceled()) {
          return false;
        }
        settings.setDistributionMask(gd.getNextChoice());
        if (maskListContainsStacks) {
          settings.setDistributionMaskSliceDepth(Math.abs(gd.getNextNumber()));
        }
      }

      SettingsManager.writeSettings(settings.build());
    }

    return getHistogramOptions();
  }

  private void addCameraOptions(final ExtendedGenericDialog gd) {
    gd.addChoice("Camera_type", SettingsManager.getCameraTypeNames(),
        CalibrationProtosHelper.getName(settings.getCameraType()), new OptionListener<Integer>() {
          @Override
          public boolean collectOptions(Integer field) {
            settings.setCameraType(SettingsManager.getCameraTypeValues()[field]);
            return collectOptions(false);
          }

          @Override
          public boolean collectOptions() {
            return collectOptions(true);
          }

          private boolean collectOptions(boolean silent) {
            final CameraType cameraType = settings.getCameraType();
            final boolean isCcd = CalibrationProtosHelper.isCcdCameraType(cameraType);
            final ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE, null);
            if (isCcd) {
              if (cameraType == CameraType.EMCCD) {
                egd.addNumericField("EM_gain", settings.getEmGain(), 2);
              }
              egd.addNumericField("Camera_gain", settings.getCameraGain(), 4, 6, "count/electron");
              egd.addNumericField("Quantum_efficiency", settings.getQuantumEfficiency(), 2, 6,
                  "electron/photon");
              egd.addNumericField("Read_noise", settings.getReadNoise(), 2, 6, "electron");
              egd.addNumericField("Bias", settings.getBias(), 0, 6, "count");
            } else if (cameraType == CameraType.SCMOS) {
              final String[] models = CameraModelManager.listCameraModels(true);
              egd.addChoice("Camera_model_name", models, settings.getCameraModelName());
              egd.addNumericField("Quantum_efficiency", settings.getQuantumEfficiency(), 2, 6,
                  "electron/photon");
            } else {
              IJ.error("Unsupported camera type " + CalibrationProtosHelper.getName(cameraType));
              return false;
            }
            egd.setSilent(silent);
            egd.showDialog(true, gd);
            if (egd.wasCanceled()) {
              return false;
            }
            if (isCcd) {
              if (cameraType == CameraType.EMCCD) {
                settings.setEmGain(Math.abs(egd.getNextNumber()));
              }
              settings.setCameraGain(Math.abs(egd.getNextNumber()));
              settings.setQuantumEfficiency(Math.abs(egd.getNextNumber()));
              settings.setReadNoise(Math.abs(egd.getNextNumber()));
              settings.setBias(Math.abs((int) egd.getNextNumber()));
            } else if (cameraType == CameraType.SCMOS) {
              settings.setCameraModelName(egd.getNextChoice());
              settings.setQuantumEfficiency(Math.abs(egd.getNextNumber()));
            }
            return true;
          }
        });
  }

  private void validateCameraOptions() {
    final CameraType cameraType = settings.getCameraType();
    final boolean isCcd = CalibrationProtosHelper.isCcdCameraType(cameraType);
    if (isCcd) {
      if (cameraType == CameraType.EMCCD) {
        ParameterUtils.isPositive("EM gain", settings.getEmGain());
      }
      ParameterUtils.isPositive("Camera gain", settings.getCameraGain());
      ParameterUtils.isPositive("Read noise", settings.getReadNoise());
      final double noiseRange = settings.getReadNoise() * settings.getCameraGain() * 4;
      ParameterUtils.isEqualOrAbove(
          "Bias must prevent clipping the read noise (@ +/- 4 StdDev) so ", settings.getBias(),
          noiseRange);

      cameraModel = createCcdCameraModel();
    } else if (cameraType == CameraType.SCMOS) {
      // Load the model
      cameraModel = CameraModelManager.load(settings.getCameraModelName());
      if (cameraModel == null) {
        throw new IllegalArgumentException(
            "Unknown camera model for name: " + settings.getCameraModelName());
      }

      // Check the width is above the selected size
      Rectangle modelBounds = cameraModel.getBounds();
      final int size = settings.getSize();
      if (modelBounds.width < size || modelBounds.height < size) {
        throw new IllegalArgumentException(String.format(
            "Camera model bounds [x=%d,y=%d,width=%d,height=%d] are smaller than "
                + "simulation size [%d]",
            modelBounds.x, modelBounds.y, modelBounds.width, modelBounds.height, size));
      }

      // Ask for a crop
      if (modelBounds.width > size || modelBounds.height > size) {
        final GenericDialog gd = new GenericDialog(TITLE);
        //@formatter:off
        gd.addMessage(String.format(
            "WARNING:\n \nCamera model bounds\n[x=%d,y=%d,width=%d,height=%d]\n"
            + "are larger than the simulation size [=%d].\n \nCrop the model?",
            modelBounds.x, modelBounds.y, modelBounds.width, modelBounds.height, size
            ));
        //@formatter:on
        gd.addCheckbox("Random_crop", settings.getRandomCrop());
        final int upperx = modelBounds.x + modelBounds.width - size;
        final int uppery = modelBounds.y + modelBounds.height - size;
        gd.addSlider("Origin_x", modelBounds.x, upperx,
            MathUtils.clip(modelBounds.x, upperx, settings.getOriginX()));
        gd.addSlider("Origin_y", modelBounds.y, uppery,
            MathUtils.clip(modelBounds.y, uppery, settings.getOriginY()));
        gd.showDialog();
        if (gd.wasCanceled()) {
          throw new IllegalArgumentException("Unknown camera model crop");
        }
        settings.setRandomCrop(gd.getNextBoolean());
        settings.setOriginX((int) gd.getNextNumber());
        settings.setOriginY((int) gd.getNextNumber());
        SettingsManager.writeSettings(settings.build());

        int ox;
        int oy;
        if (settings.getRandomCrop()) {
          final RandomDataGenerator rg = new RandomDataGenerator(createRandomGenerator());
          ox = rg.nextInt(modelBounds.x, upperx);
          oy = rg.nextInt(modelBounds.y, uppery);
        } else {
          ox = settings.getOriginX();
          oy = settings.getOriginY();
        }
        final Rectangle bounds = new Rectangle(ox, oy, size, size);
        cameraModel = cameraModel.crop(bounds, false);
        modelBounds = cameraModel.getBounds();
        if (modelBounds.width != size || modelBounds.height != size) {
          throw new IllegalArgumentException("Failed to crop camera model to bounds: " + bounds);
        }
      }
    } else {
      throw new IllegalArgumentException(
          "Unsupported camera type: " + CalibrationProtosHelper.getName(cameraType));
    }
  }

  private void addCameraOptionsHeader(StringBuilder sb) {
    final CameraType cameraType = settings.getCameraType();
    final boolean isCcd = CalibrationProtosHelper.isCcdCameraType(cameraType);
    if (isCcd) {
      if (cameraType == CameraType.EMCCD) {
        addHeaderLine(sb, "EM_gain", settings.getEmGain());
      }
      addHeaderLine(sb, "Camera_gain", settings.getCameraGain());
      addHeaderLine(sb, "Quantum_efficiency", getQuantumEfficiency());
      addHeaderLine(sb, "Read_noise", settings.getReadNoise());
      addHeaderLine(sb, "Bias", settings.getBias());
    } else if (cameraType == CameraType.SCMOS) {
      addHeaderLine(sb, "Camera_model_name", settings.getCameraModelName());
      addHeaderLine(sb, "Quantum_efficiency", getQuantumEfficiency());
    }
  }

  private double getQuantumEfficiency() {
    final double qe = settings.getQuantumEfficiency();
    return (qe > 0 && qe < 1) ? qe : 1;
  }

  /**
   * Creates the likelihood function. This is used for CRLB computation.
   */
  private void createLikelihoodFunction() {
    final CameraType cameraType = settings.getCameraType();
    final boolean isCcd = CalibrationProtosHelper.isCcdCameraType(cameraType);
    fiFunction = new BasePoissonFisherInformation[settings.getSize() * settings.getSize()];
    if (isCcd) {
      BasePoissonFisherInformation fi;
      final CreateDataSettingsHelper helper = new CreateDataSettingsHelper(settings);
      final double readNoise = helper.getReadNoiseInCounts();

      if (cameraType == CameraType.EMCCD) {
        // We only want the amplification (without QE applied)
        final double amp = helper.getAmplification();
        // This should be interpolated from a stored curve
        final InterpolatedPoissonFisherInformation i = CameraModelFisherInformationAnalysis
            .loadFunction(CameraModelFisherInformationAnalysis.CameraType.EM_CCD, amp, readNoise);
        if (i == null) {
          throw new IllegalStateException(
              "No stored Fisher information for EM-CCD camera with gain " + amp + " and noise "
                  + readNoise);
        }
        fi = i;
      } else {
        // This is fast enough to compute dynamically.
        // Read noise is in electrons so use directly.
        fi = new PoissonGaussianFisherInformation(settings.getReadNoise());
      }
      Arrays.fill(fiFunction, fi);
    } else if (cameraType == CameraType.SCMOS) {
      // Build per-pixel likelihood function.
      // Get the normalised variance per pixel.
      final float[] v = cameraModel.getNormalisedVariance(cameraModel.getBounds());
      // Build the function
      for (int i = 0; i < fiFunction.length; i++) {
        fiFunction[i] = new PoissonGaussianFisherInformation(Math.sqrt(v[i]));
      }
    } else {
      throw new IllegalArgumentException(
          "Unsupported camera type: " + CalibrationProtosHelper.getName(cameraType));
    }
  }

  /**
   * Check if there are any suitable PSF images open. If so add a choice to allow the selection of
   * the Gaussian or Image PSF model. If no PSF images are open then add options for the wavelength
   * and NA for the simulated microscope.
   *
   * @param gd the gd
   */
  private void addPsfOptions(final ExtendedGenericDialog gd) {
    gd.addMessage("--- PSF Model ---");
    final List<String> imageNames = PsfCombiner.createImageList();
    final TurboList<String> availableModels = new TurboList<>();
    availableModels.add(PSF_MODELS[PSF_MODEL_GAUSSIAN]);
    availableModels.add(PSF_MODELS[PSF_MODEL_AIRY]);
    final String[] images;
    if (!imageNames.isEmpty()) {
      availableModels.add(PSF_MODELS[PSF_MODEL_IMAGE]);
      images = imageNames.toArray(new String[imageNames.size()]);
    } else {
      images = null;
    }
    final String[] astigmatismModels = AstigmatismModelManager.listAstigmatismModels(false, true);
    if (astigmatismModels.length != 0) {
      availableModels.add(PSF_MODELS[PSF_MODEL_ASTIGMATISM]);
    }
    final String[] models = availableModels.toArray(new String[availableModels.size()]);
    gd.addChoice("PSF_model", models, settings.getPsfModel(), new OptionListener<Integer>() {
      @Override
      public boolean collectOptions(Integer value) {
        settings.setPsfModel(models[value]);
        return collectOptions(false);
      }

      @Override
      public boolean collectOptions() {
        return collectOptions(true);
      }

      private boolean collectOptions(boolean silent) {
        final ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE, null);
        egd.addMessage("Configure the " + settings.getPsfModel() + " PSF model");

        int type = 0;

        // Get the image
        if (settings.getPsfModel().equals(PSF_MODELS[PSF_MODEL_IMAGE])) {
          egd.addChoice("PSF_image", images, settings.getPsfImageName());
        } else if (settings.getPsfModel().equals(PSF_MODELS[PSF_MODEL_ASTIGMATISM])) {
          type = 1;
          egd.addChoice("Astigmatism_model", astigmatismModels, settings.getAstigmatismModel());
          egd.addMessage(TextUtils.wrap("Note: The pixel size of the astigmatism model should match"
              + " the pixel pitch if fitting of the data is to be performed (i.e. fitting requires"
              + " the astigmatism model to be calibrated to the image). If not then the model will"
              + " be optionally converted before the simulation.", 80));
        } else {
          // Get the width of the model
          type = 2;
          egd.addNumericField("Depth-of-focus (nm)", settings.getDepthOfFocus(), 2);
          egd.addCheckbox("Enter_width", settings.getEnterWidth());
          egd.addNumericField("PSF_SD (nm)", settings.getPsfSd(), 2);
          egd.addMessage("Or compute from optics:");
          egd.addNumericField("Wavelength (nm)", settings.getWavelength(), 2);
          egd.addNumericField("Numerical_aperture", settings.getNumericalAperture(), 2);
        }
        egd.setSilent(silent);
        egd.showDialog(true, gd);
        if (egd.wasCanceled()) {
          return false;
        }
        if (type == 0) {
          settings.setPsfImageName(egd.getNextChoice());
        } else if (type == 1) {
          settings
              .setAstigmatismModel(AstigmatismModelManager.removeFormatting(egd.getNextChoice()));
        } else {
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

  private void validatePsfOptions() {
    if (settings.getPsfModel().equals(PSF_MODELS[PSF_MODEL_ASTIGMATISM])) {
      psfModelType = PSF_MODEL_ASTIGMATISM;
      final AstigmatismModel model =
          AstigmatismModelManager.getModel(settings.getAstigmatismModel());
      if (model == null) {
        throw new IllegalArgumentException(
            "Failed to load model: " + settings.getAstigmatismModel());
      }
    } else if (settings.getPsfModel().equals(PSF_MODELS[PSF_MODEL_IMAGE])) {
      psfModelType = PSF_MODEL_IMAGE;
    } else {
      psfModelType =
          (settings.getPsfModel().equals(PSF_MODELS[PSF_MODEL_GAUSSIAN])) ? PSF_MODEL_GAUSSIAN
              : PSF_MODEL_AIRY;
      ParameterUtils.isAboveZero("Depth-of-focus", settings.getDepthOfFocus());
      if (settings.getEnterWidth()) {
        ParameterUtils.isAboveZero("PSF SD", settings.getPsfSd());
      } else {
        ParameterUtils.isAboveZero("Wavelength", settings.getWavelength());
        ParameterUtils.isAboveZero("NA", settings.getNumericalAperture());
        ParameterUtils.isBelow("NA", settings.getNumericalAperture(), 2);
      }
    }
  }

  private boolean getHistogramOptions() {
    GenericDialog gd;
    if (settings.getShowHistograms() && settings.getChooseHistograms()) {
      gd = new GenericDialog(TITLE);
      gd.addMessage("Select the histograms to display");
      for (int i = 0; i < displayHistograms.length; i++) {
        gd.addCheckbox(NAMES[i].replace(' ', '_'), displayHistograms[i]);
      }
      gd.showDialog();
      if (gd.wasCanceled()) {
        return false;
      }
      for (int i = 0; i < displayHistograms.length; i++) {
        displayHistograms[i] = gd.getNextBoolean();
      }
    }
    return true;
  }

  /**
   * Show a dialog allowing the parameters for a simulation to be performed.
   *
   * @return True if the parameters were collected
   */
  private boolean showDialog() {
    // In track mode we do not need a time, illumination model or blinking model.
    // Fixed length tracks will be drawn, non-overlapping in time. This is the simplest
    // simulation for moving molecules

    ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

    settings = SettingsManager.readCreateDataSettings(0).toBuilder();

    if (settings.getStepsPerSecond() < 1) {
      settings.setStepsPerSecond(1);
    }

    final String[] backgroundImages = createBackgroundImageList();
    gd.addNumericField("Pixel_pitch (nm)", settings.getPixelPitch(), 2);
    gd.addNumericField("Size (px)", settings.getSize(), 0);
    gd.addNumericField("Depth (nm)", settings.getDepth(), 0);
    gd.addCheckbox("Fixed_depth", settings.getFixedDepth());
    if (!trackMode) {
      gd.addNumericField("Seconds", settings.getSeconds(), 1);
    }
    gd.addNumericField("Exposure_time (ms)", settings.getExposureTime(), 1);
    gd.addSlider("Steps_per_second", 1, MathUtils.clip(15, 1000, settings.getStepsPerSecond()),
        settings.getStepsPerSecond());
    if (!trackMode) {
      gd.addChoice("Illumination", ILLUMINATION, settings.getIllumination());
      gd.addNumericField("Pulse_interval", settings.getPulseInterval(), 0);
      gd.addNumericField("Pulse_ratio", settings.getPulseRatio(), 2);
    }
    if (backgroundImages != null) {
      gd.addChoice("Background_image", backgroundImages, settings.getBackgroundImage());
    }

    if (extraOptions) {
      gd.addCheckbox("No_poisson_noise", !settings.getPoissonNoise());
    }
    gd.addNumericField("Background (photons)", settings.getBackground(), 2);

    addCameraOptions(gd);

    addPsfOptions(gd);

    gd.addMessage("--- Fluorophores ---");
    gd.addChoice("Distribution", DISTRIBUTION, settings.getDistribution());
    gd.addNumericField("Particles", settings.getParticles(), 0);
    gd.addCheckbox("Compound_molecules", settings.getCompoundMolecules());
    gd.addNumericField("Diffusion_rate (um^2/sec)", settings.getDiffusionRate(), 2);
    final String[] diffusionTypes = SettingsManager.getNames((Object[]) DiffusionType.values());
    gd.addChoice("Diffusion_type", diffusionTypes, diffusionTypes[CreateDataSettingsHelper
        .getDiffusionType(settings.getDiffusionType()).ordinal()]);
    gd.addSlider("Fixed_fraction (%)", 0, 100, settings.getFixedFraction() * 100);
    gd.addChoice("Confinement", CONFINEMENT, settings.getConfinement());
    gd.addNumericField("Photons (sec^-1)", settings.getPhotonsPerSecond(), 0);
    // We cannot use the correlation moe with fixed life time tracks
    final String[] dist =
        (trackMode) ? Arrays.copyOf(PHOTON_DISTRIBUTION, PHOTON_DISTRIBUTION.length - 1)
            : PHOTON_DISTRIBUTION;
    gd.addChoice("Photon_distribution", dist, settings.getPhotonDistribution());
    gd.addNumericField("On_time (ms)", settings.getTOn(), 2);
    if (!trackMode) {
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

    if (gd.wasCanceled()) {
      return false;
    }

    settings.setPixelPitch(Math.abs(gd.getNextNumber()));
    settings.setSize(Math.abs((int) gd.getNextNumber()));
    settings.setDepth(Math.abs(gd.getNextNumber()));
    settings.setFixedDepth(gd.getNextBoolean());
    if (!trackMode) {
      settings.setSeconds(Math.abs(gd.getNextNumber()));
    }
    settings.setExposureTime(Math.abs(gd.getNextNumber()));
    settings.setStepsPerSecond(Math.abs(gd.getNextNumber()));
    if (!trackMode) {
      settings.setIllumination(gd.getNextChoice());
      settings.setPulseInterval(Math.abs((int) gd.getNextNumber()));
      settings.setPulseRatio(Math.abs(gd.getNextNumber()));
    }
    if (backgroundImages != null) {
      settings.setBackgroundImage(gd.getNextChoice());
    }

    if (extraOptions) {
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
    if (!trackMode) {
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
    minSnrT1 = settings.getMinSnrT1();
    minSnrTn = settings.getMinSnrTN();

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
    if (settings.getMinSnrT1() < settings.getMinSnrTN()) {
      final double tmp = settings.getMinSnrT1();
      settings.setMinSnrT1(settings.getMinSnrTN());
      settings.setMinSnrTN(tmp);
    }

    gd.collectOptions();

    // Save before validation so that the current values are preserved.
    SettingsManager.writeSettings(settings.build());

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Pixel Pitch", settings.getPixelPitch());
      ParameterUtils.isAboveZero("Size", settings.getSize());
      if (!settings.getFixedDepth()) {
        ParameterUtils.isPositive("Depth", settings.getDepth());
      }
      if (!trackMode) {
        ParameterUtils.isAboveZero("Seconds", settings.getSeconds());
      }
      ParameterUtils.isAboveZero("Exposure time", settings.getExposureTime());
      ParameterUtils.isAboveZero("Steps per second", settings.getStepsPerSecond());
      ParameterUtils.isPositive("Background", settings.getBackground());
      ParameterUtils.isAboveZero("Particles", settings.getParticles());
      ParameterUtils.isAboveZero("Photons", settings.getPhotonsPerSecond());
      ParameterUtils.isPositive("Diffusion rate", settings.getDiffusionRate());
      ParameterUtils.isPositive("Fixed fraction", settings.getFixedFraction());
      ParameterUtils.isPositive("Pulse interval", settings.getPulseInterval());
      ParameterUtils.isAboveZero("Pulse ratio", settings.getPulseRatio());
      ParameterUtils.isAboveZero("tOn", settings.getTOn());
      if (!trackMode) {
        ParameterUtils.isAboveZero("tOff Short", settings.getTOffShort());
        ParameterUtils.isAboveZero("tOff Long", settings.getTOffLong());
        ParameterUtils.isPositive("n-Blinks Short", settings.getNBlinksShort());
        ParameterUtils.isPositive("n-Blinks Long", settings.getNBlinksLong());
      }
      ParameterUtils.isPositive("Min photons", settings.getMinPhotons());
      ParameterUtils.isPositive("Min SNR t1", settings.getMinSnrT1());
      ParameterUtils.isPositive("Min SNR tN", settings.getMinSnrTN());
      ParameterUtils.isPositive("Histogram bins", settings.getHistogramBins());
      ParameterUtils.isPositive("Density radius", settings.getDensityRadius());

      validateCameraOptions();
      validatePsfOptions();
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    if (gd.invalidNumber()) {
      return false;
    }

    if (!getHistogramOptions()) {
      return false;
    }

    String[] maskImages = null;
    if (settings.getDistribution().equals(DISTRIBUTION[MASK])) {
      maskImages = createDistributionImageList();
      if (maskImages != null) {
        gd = new ExtendedGenericDialog(TITLE);
        gd.addMessage("Select the mask image for the distribution");
        gd.addChoice("Distribution_mask", maskImages, settings.getDistributionMask());
        if (maskListContainsStacks) {
          gd.addNumericField("Distribution_slice_depth (nm)",
              settings.getDistributionMaskSliceDepth(), 0);
        }
        gd.showDialog();
        if (gd.wasCanceled()) {
          return false;
        }
        settings.setDistributionMask(gd.getNextChoice());
        if (maskListContainsStacks) {
          settings.setDistributionMaskSliceDepth(Math.abs(gd.getNextNumber()));
        }
      }
    } else if (settings.getDistribution().equals(DISTRIBUTION[GRID])) {
      gd = new ExtendedGenericDialog(TITLE);
      gd.addMessage("Select grid for the distribution");
      gd.addNumericField("Cell_size", settings.getCellSize(), 0);
      gd.addSlider("p-binary", 0, 1, settings.getProbabilityBinary());
      gd.addNumericField("Min_binary_distance (nm)", settings.getMinBinaryDistance(), 0);
      gd.addNumericField("Max_binary_distance (nm)", settings.getMaxBinaryDistance(), 0);
      gd.showDialog();
      if (gd.wasCanceled()) {
        return false;
      }
      settings.setCellSize((int) gd.getNextNumber());
      settings.setProbabilityBinary(gd.getNextNumber());
      settings.setMinBinaryDistance(gd.getNextNumber());
      settings.setMaxBinaryDistance(gd.getNextNumber());

      // Check arguments
      try {
        ParameterUtils.isAboveZero("Cell size", settings.getCellSize());
        ParameterUtils.isPositive("p-binary", settings.getProbabilityBinary());
        ParameterUtils.isEqualOrBelow("p-binary", settings.getProbabilityBinary(), 1);
        ParameterUtils.isPositive("Min binary distance", settings.getMinBinaryDistance());
        ParameterUtils.isPositive("Max binary distance", settings.getMaxBinaryDistance());
        ParameterUtils.isEqualOrBelow("Min binary distance", settings.getMinBinaryDistance(),
            settings.getMaxBinaryDistance());
      } catch (final IllegalArgumentException ex) {
        IJ.error(TITLE, ex.getMessage());
        return false;
      }
    }

    SettingsManager.writeSettings(settings.build());

    if (settings.getDiffusionRate() > 0 && settings.getFixedFraction() < 1) {
      if (settings.getConfinement().equals(CONFINEMENT[CONFINEMENT_SPHERE])) {
        gd = new ExtendedGenericDialog(TITLE);
        gd.addMessage("Select the sphere radius for the diffusion confinement");
        gd.addSlider("Confinement_radius (nm)", 0, 2000, settings.getConfinementRadius());
        gd.showDialog();
        if (gd.wasCanceled()) {
          return false;
        }
        settings.setConfinementRadius(gd.getNextNumber());
      } else if (settings.getConfinement().equals(CONFINEMENT[CONFINEMENT_MASK])) {
        if (maskImages == null) {
          maskImages = createDistributionImageList();
        }
        if (maskImages != null) {
          gd = new ExtendedGenericDialog(TITLE);
          gd.addMessage("Select the mask image for the diffusion confinement");
          gd.addChoice("Confinement_mask", maskImages, settings.getConfinementMask());
          if (maskListContainsStacks) {
            gd.addNumericField("Confinement_slice_depth (nm)",
                settings.getConfinementMaskSliceDepth(), 0);
          }
          gd.showDialog();
          if (gd.wasCanceled()) {
            return false;
          }
          settings.setConfinementMask(gd.getNextChoice());
          if (maskListContainsStacks) {
            settings.setConfinementMaskSliceDepth(Math.abs(gd.getNextNumber()));
          }
        }
      }
    }

    SettingsManager.writeSettings(settings.build());

    if (settings.getCompoundMolecules()) {
      // Show a second dialog where the molecule configuration is specified
      gd = new ExtendedGenericDialog(TITLE);

      gd.addMessage("Specify the compound molecules");
      gd.addTextAreas(settings.getCompoundText(), null, 20, 80);
      gd.addCheckbox("Enable_2D_diffusion", settings.getDiffuse2D());
      gd.addCheckbox("Rotate_initial_orientation", settings.getRotateInitialOrientation());
      gd.addCheckbox("Rotate_during_simulation", settings.getRotateDuringSimulation());
      gd.addCheckbox("Enable_2D_rotation", settings.getRotate2D());
      gd.addCheckbox("Show_example_compounds", false);

      if (ImageJUtils.isShowGenericDialog()) {
        @SuppressWarnings("rawtypes")
        final Vector v = gd.getCheckboxes();
        final Checkbox cb = (Checkbox) v.get(v.size() - 1);
        cb.addItemListener(this);
      }

      gd.showDialog();
      if (gd.wasCanceled()) {
        return false;
      }

      settings.setCompoundText(gd.getNextText());
      settings.setDiffuse2D(gd.getNextBoolean());
      settings.setRotateInitialOrientation(gd.getNextBoolean());
      settings.setRotateDuringSimulation(gd.getNextBoolean());
      settings.setRotate2D(gd.getNextBoolean());

      if (gd.getNextBoolean()) {
        logExampleCompounds();
        return false;
      }
    }

    SettingsManager.writeSettings(settings.build());

    gd = new ExtendedGenericDialog(TITLE);
    gd.addMessage("Configure the photon distribution: " + settings.getPhotonDistribution());
    if (PHOTON_DISTRIBUTION[PHOTON_CUSTOM].equals(settings.getPhotonDistribution())) {
      // Nothing more to be done
      return true;
    } else if (PHOTON_DISTRIBUTION[PHOTON_UNIFORM].equals(settings.getPhotonDistribution())) {
      gd.addNumericField("Max_Photons (sec^-1)", settings.getPhotonsPerSecondMaximum(), 0);
    } else if (PHOTON_DISTRIBUTION[PHOTON_GAMMA].equals(settings.getPhotonDistribution())) {
      gd.addNumericField("Photon_shape", settings.getPhotonShape(), 2);
    } else if (PHOTON_DISTRIBUTION[PHOTON_CORRELATED].equals(settings.getPhotonDistribution())) {
      gd.addNumericField("Correlation (to total tOn)", settings.getCorrelation(), 2);
    } else {
      // Nothing more to be done
      return true;
    }

    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }

    try {
      if (PHOTON_DISTRIBUTION[PHOTON_UNIFORM].equals(settings.getPhotonDistribution())) {
        settings.setPhotonsPerSecondMaximum(Math.abs((int) gd.getNextNumber()));
        if (settings.getPhotonsPerSecondMaximum() < settings.getPhotonsPerSecond()) {
          settings.setPhotonsPerSecondMaximum(settings.getPhotonsPerSecond());
        }
      } else if (PHOTON_DISTRIBUTION[PHOTON_GAMMA].equals(settings.getPhotonDistribution())) {
        settings.setPhotonShape(Math.abs(gd.getNextNumber()));
        ParameterUtils.isAbove("Photon shape", settings.getPhotonShape(), 0);
      } else if (PHOTON_DISTRIBUTION[PHOTON_CORRELATED].equals(settings.getPhotonDistribution())) {
        settings.setCorrelation(gd.getNextNumber());
        ParameterUtils.isEqualOrBelow("Correlation", settings.getCorrelation(), 1);
        ParameterUtils.isEqualOrAbove("Correlation", settings.getCorrelation(), -1);
      }
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    SettingsManager.writeSettings(settings.build());

    return true;
  }

  /**
   * Build a list of suitable background images. Images must be greyscale.
   *
   * @return the list of image titles
   */
  private static String[] createBackgroundImageList() {
    final int[] idList = WindowManager.getIDList();
    if (idList != null) {
      final String[] list = new String[idList.length + 1];
      list[0] = "[None]";
      int count = 1;
      for (final int id : idList) {
        final ImagePlus imp = WindowManager.getImage(id);
        // Image must be square and greyscale
        if (imp != null && imp.getWidth() == imp.getHeight()
            && (imp.getBitDepth() == 8 || imp.getBitDepth() == 16 || imp.getBitDepth() == 32)) {
          list[count++] = imp.getTitle();
        }
      }
      if (count == 1) {
        return null;
      }
      return Arrays.copyOf(list, count);
    }
    return null;
  }

  /**
   * Build a list of suitable distribution images. Images must be square.
   *
   * @return the list of image titles
   */
  private String[] createDistributionImageList() {
    maskListContainsStacks = false;
    final int[] idList = WindowManager.getIDList();
    if (idList != null) {
      final String[] list = new String[idList.length + 1];
      list[0] = "[None]";
      int count = 1;
      for (final int id : idList) {
        final ImagePlus imp = WindowManager.getImage(id);
        if (imp != null && imp.getWidth() == imp.getHeight()) {
          list[count++] = imp.getTitle();
          if (imp.getStackSize() > 1) {
            maskListContainsStacks = true;
          }
        }
      }
      if (count == 1) {
        return null;
      }
      return Arrays.copyOf(list, count);
    }
    return null;
  }

  @Override
  public void itemStateChanged(ItemEvent event) {
    // When the checkbox is clicked, output example compounds to the ImageJ log
    final Checkbox cb = (Checkbox) event.getSource();
    if (cb.getState()) {
      cb.setState(false);

      logExampleCompounds();
    }
  }

  private static void logExampleCompounds() {
    comment(TITLE + " example compounds");
    IJ.log("");
    comment("Compounds are described using parseable text.");
    comment("Missing fields are initialised to the default (0).");
    comment("Multiple compounds can be combined using fractional ratios.");
    comment("Coordinates are specified in nanometres.");
    comment("Coordinates describe the relative positions of atoms in the molecule;"
        + " the molcule will have a randomly assigned XYZ position for its centre-of-mass."
        + " Rotation will be about the centre-of-mass.");
    IJ.log("");

    final Molecule.Builder mb = Molecule.newBuilder();
    addAtom(mb, 10, 1, 1, 1);
    mb.setDiffusionRate(0.5);
    mb.setDiffusionType(DiffusionType.RANDOM_WALK.toString());
    final Molecule m1 = mb.build();

    mb.clear();
    addAtom(mb, 30, 0, 0, 0);
    addAtom(mb, 20, 1000, 0, 0);
    mb.setDiffusionRate(1);
    mb.setDiffusionType(DiffusionType.GRID_WALK.toString());
    final Molecule m2 = mb.build();

    // Create a hexamer big enough to see with the default pixel pitch
    mb.clear();
    addAtom(mb, 1, 0, 0, 0);
    addAtom(mb, 1, 1000, 0, 0);
    addAtom(mb, 1, 1500, 866, 0);
    addAtom(mb, 1, 1000, 1732, 0);
    addAtom(mb, 1, 0, 1732, 0);
    addAtom(mb, 1, -500, 866, 0);
    final Molecule m3 = mb.build();

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
    final Molecule m1a = m1.toBuilder().setFraction(2).build();
    final Molecule m2a = m2.toBuilder().setFraction(1).build();
    demo(m1a, m2a);
  }

  private static void addAtom(Molecule.Builder mb, double mass, double x, double y, double z) {
    final Atom.Builder ab = mb.addAtomBuilder();
    ab.setMass(mass);
    ab.setX(x);
    ab.setY(y);
    ab.setZ(z);
  }

  private static void demo(Molecule... molecules) {
    final Mixture.Builder builder = Mixture.newBuilder();
    for (final Molecule m : molecules) {
      builder.addMolecule(m);
    }

    // The toString() method is more verbose than JSON but easier to read
    IJ.log(builder.toString());
    IJ.log("");
  }

  private static void comment(String text) {
    IJ.log(TextUtils.wrap("# " + text, 80, "\n# ", false));
  }

  private List<CompoundMoleculeModel> createCompoundMolecules() {
    // Diffusion rate is um^2/sec. Convert to pixels per simulation frame.
    final double diffusionFactor =
        (1000000.0 / (settings.getPixelPitch() * settings.getPixelPitch()))
            / settings.getStepsPerSecond();

    List<CompoundMoleculeModel> compounds;
    if (settings.getCompoundMolecules()) {
      // Try and load the compounds from the text specification
      try {
        // Convert from the serialised objects to the compound model
        final String text = settings.getCompoundText();
        final Mixture.Builder builder = Mixture.newBuilder();
        TextFormat.merge(text, builder);

        compounds = new ArrayList<>(builder.getMoleculeCount());
        int id = 1;
        compoundNames = new ArrayList<>(builder.getMoleculeCount());
        for (final Molecule m : builder.getMoleculeList()) {
          final MoleculeModel[] molecules = new MoleculeModel[m.getAtomCount()];
          for (int i = 0; i < molecules.length; i++) {
            final AtomOrBuilder a = m.getAtomOrBuilder(i);
            molecules[i] = new MoleculeModel(a.getMass(), a.getX(), a.getY(), a.getZ());
          }
          final CompoundMoleculeModel cm =
              new CompoundMoleculeModel(id++, 0, 0, 0, Arrays.asList(molecules));
          cm.setFraction(m.getFraction());
          cm.setDiffusionRate(m.getDiffusionRate() * diffusionFactor);
          cm.setDiffusionType(DiffusionType.fromString(m.getDiffusionType()));
          compounds.add(cm);
          compoundNames.add(String.format("Fraction=%s, D=%s um^2/s",
              MathUtils.rounded(cm.getFraction()), MathUtils.rounded(m.getDiffusionRate())));
        }

        // Convert coordinates from nm to pixels
        final double scaleFactor = 1.0 / settings.getPixelPitch();
        for (final CompoundMoleculeModel c : compounds) {
          c.scale(scaleFactor);
        }
      } catch (final IOException ex) {
        IJ.error(TITLE, "Unable to create compound molecules");
        return null;
      }
    } else {
      // Create a simple compound with one molecule at the origin
      compounds = new ArrayList<>(1);
      final CompoundMoleculeModel m =
          new CompoundMoleculeModel(1, 0, 0, 0, Arrays.asList(new MoleculeModel(0, 0, 0, 0)));
      m.setDiffusionRate(settings.getDiffusionRate() * diffusionFactor);
      m.setDiffusionType(CreateDataSettingsHelper.getDiffusionType(settings.getDiffusionType()));
      compounds.add(m);
    }
    return compounds;
  }

  private static int seedAddition;
  private boolean resetSeed = true;

  private enum SeedMode {
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
     * Set to true when the seed should not have a time dependent offset. This ensure that the
     * plugin can be run at any time from start-up and the simulation will be reproducible.
     *
     * @return true, if the seed offset should be identical
     */
    abstract boolean identicalOffset();

    /**
     * Set to true when the seed should have the same addition for each plugin execution. Set to
     * false will ensure subsequent runs of the plugin produce different simulations.
     *
     * @return true, if the seed addition should be reset for each plugin execution
     */
    abstract boolean identicalAddition();
  }

  private final SeedMode seedMode = SeedMode.DEFAULT;

  private long getSeedOffset() {
    return (seedMode.identicalOffset()) ? 1
        : System.currentTimeMillis() + System.identityHashCode(this);
  }

  private int getSeedAddition() {
    if (seedMode.identicalAddition() && resetSeed) {
      // Reset only once per plugin execution
      resetSeed = false;
      seedAddition = 0;
    }
    // Increment the seed to ensure that new generators are created at the same system time point
    return seedAddition++;
  }

  /**
   * Creates a random generator.
   *
   * <p>The generators used in the simulation can be adjusted by changing this method.
   */
  private RandomGenerator createRandomGenerator() {
    return new Well44497b(getSeedOffset() + getSeedAddition());
  }

  private static void setBenchmarkResults(ImagePlus imp, MemoryPeakResults results) {
    if (imp != null) {
      benchmarkImageId = imp.getID();
      benchmarkResultsName = results.getName();
      MemoryPeakResults.addResults(results);
    } else {
      benchmarkImageId = 0;
      benchmarkResultsName = "";
    }
  }

  /**
   * Gets the benchmark image.
   *
   * @return the image
   */
  public static ImagePlus getImage() {
    return WindowManager.getImage(benchmarkImageId);
  }

  /**
   * Gets the benchmark image Id.
   *
   * @return the image Id
   */
  public static int getImageId() {
    return benchmarkImageId;
  }

  /**
   * Gets the benchmark results.
   *
   * @return the results
   */
  public static ImmutableMemoryPeakResults getResults() {
    final MemoryPeakResults r = MemoryPeakResults.getResults(benchmarkResultsName);
    if (r == null) {
      return null;
    }
    return new ImmutableMemoryPeakResults(r);
  }

  /**
   * Load benchmark data using an open image and a XYZ text file.
   */
  private void loadBenchmarkData() {
    // Note: Do not reset memory until failure. This allows the load method to use the
    // last simulation parameters to set settings.

    if (!showLoadDialog()) {
      // resetMemory();
      return;
    }

    // Load the image
    final ImagePlus imp = WindowManager.getImage(benchmarkImage);
    if (imp == null) {
      IJ.error(TITLE, "No benchmark image: " + benchmarkImage);
      // resetMemory();
      return;
    }

    // Load the results
    final MemoryPeakResults results = getSimulationResults();
    if (results == null) {
      IJ.error(TITLE, "No benchmark results: " + benchmarkResultsName);
      // resetMemory();
      return;
    }
    results.setName(imp.getTitle() + " (Results)");
    results.setBounds(new Rectangle(0, 0, imp.getWidth(), imp.getHeight()));
    final IJImageSource imageSource = new IJImageSource(imp);
    results.setSource(imageSource);

    // Load the settings as these are used in the dialog
    settings = SettingsManager.readCreateDataSettings(0).toBuilder();

    simulationParameters = showSimulationParametersDialog(imp, results);
    if (simulationParameters != null) {
      setBackground(results);
      setNoise(results, imp);
      setBenchmarkResults(imp, results);
      IJ.showStatus("Loaded " + TextUtils.pleural(results.size(), "result"));
    } else {
      resetMemory();
    }
  }

  /**
   * Sets the background in the results if missing.
   *
   * @param results the results
   */
  private static void setBackground(MemoryPeakResults results) {
    // Loaded results do not have a local background.
    if (results.hasBackground()) {
      return;
    }

    // Simple fix is to use the global photon background.
    // TODO - Subtract the spots from the local region and compute the true local background.
    // Note this requires knowing the PSF width. If this is a loaded ground truth dataset then
    // it probably will not have Gaussian widths.
    results.setZeroBackground(IntensityUnit.PHOTON, (float) simulationParameters.background);
  }

  /**
   * Sets the noise in the results if missing.
   *
   * @param results the results
   * @param imp the imp
   */
  private static void setNoise(MemoryPeakResults results, ImagePlus imp) {
    // Loaded results do not have noise
    if (results.hasNoise()) {
      return;
    }

    IJ.showStatus("Estimating noise ...");

    // Compute noise per frame
    final ImageStack stack = imp.getImageStack();
    final int width = stack.getWidth();
    final int height = stack.getHeight();
    final IJImageSource source = new IJImageSource(imp);
    final float[] noise = new float[source.getFrames() + 1];
    source.setReadHint(ReadHint.SEQUENTIAL);
    source.open();
    for (int slice = 1; slice < noise.length; slice++) {
      final float[] data = source.next();
      // Use the trimmed method as there may be a lot of spots in the frame
      noise[slice] = FitWorker.estimateNoise(data, width, height,
          NoiseEstimatorMethod.QUICK_RESIDUALS_LEAST_TRIMMED_OF_SQUARES);
    }

    // Statistics stats = Statistics.create(Arrays.copyOfRange(noise, 1, noise.length));
    // System.out.printf("Noise = %.3f +/- %.3f (%d)\n", stats.getMean(),
    // stats.getStandardDeviation(), stats.getN());

    // Convert noise units from counts to the result format
    final TypeConverter<IntensityUnit> c = results.getIntensityConverter(IntensityUnit.COUNT);
    for (int i = 0; i < noise.length; i++) {
      noise[i] = c.convertBack(noise[i]);
    }

    results.forEach((PeakResultProcedure) result -> {
      if (result.getFrame() < noise.length) {
        result.setNoise(noise[result.getFrame()]);
      }
    });
  }

  private static MemoryPeakResults getSimulationResults() {
    if (benchmarkAuto) {
      // Load directly from a results file. This is mainly to be used to load simulations
      // saved to memory then saved to file.
      final PeakResultsReader r = new PeakResultsReader(benchmarkFile);
      final MemoryPeakResults results = r.getResults();
      if (results != null) {
        ResultsManager.checkCalibration(results);
        return results;
      }
    }

    // Load using a universal text file
    if (loadSettings == null) {
      loadSettings = SettingsManager.readLoadLocalisationsSettings(0).toBuilder();
    }
    // String tmp = loadSettings.getLocalisationsFilename();
    loadSettings.setLocalisationsFilename(benchmarkFile);
    // TODO - This could be configurable to ignore fields that are not required,
    // e.g. the dataset name
    final LocalisationList localisations = LoadLocalisations.loadLocalisations(loadSettings);
    // loadSettings.setLocalisationsFilename(tmp);
    SettingsManager.writeSettings(loadSettings.build());
    if (localisations == null || localisations.isEmpty()) {
      return null;
    }

    return localisations.toPeakResults("Dummy");
  }

  private SimulationParameters showSimulationParametersDialog(ImagePlus imp,
      MemoryPeakResults results) {
    final int molecules = results.size();

    // Get the missing parameters from the user
    boolean fullSimulation = false;
    double sd = -1;

    if (!results.convertToPreferredUnits()) {
      IJ.error(TITLE,
          String.format("Results should be in the preferred units (%s,%s)",
              UnitHelper.getName(MemoryPeakResults.PREFERRED_DISTANCE_UNIT),
              UnitHelper.getName(MemoryPeakResults.PREFERRED_INTENSITY_UNIT)));
      return null;
    }

    // Get these from the data
    final RawResultProcedure sp = new RawResultProcedure(results);
    sp.getBixyz();
    final float[] signal = sp.intensity;
    float[] limits = MathUtils.limits(signal);
    final double minSignal = limits[0];
    final double maxSignal = limits[1];
    final double signalPerFrame = MathUtils.sum(signal) / molecules;

    final float[] depths = sp.z;
    limits = MathUtils.limits(depths);
    float depth = Math.max(Math.abs(limits[0]), Math.abs(limits[1]));
    final boolean fixedDepth = Double.compare(limits[0], limits[1]) == 0;

    final CalibrationWriter cal = results.getCalibrationWriter();
    final String iUnits = " " + UnitHelper.getName(cal.getIntensityUnit());
    final String zUnits = " " + UnitHelper.getName(cal.getDistanceUnit());

    // Get this from the user
    double background = -1;

    // Use last simulation parameters for missing settings.
    // This is good if we are re-running the plugin to load data.
    Rectangle lastCameraBounds = null;
    if (simulationParameters != null && simulationParameters.isLoaded()) {
      fullSimulation = simulationParameters.fullSimulation;
      sd = simulationParameters.sd;
      background = simulationParameters.background;
      if (!cal.hasBias()) {
        cal.setBias(simulationParameters.bias);
      }
      if (!cal.hasCountPerPhoton()) {
        cal.setCountPerPhoton(simulationParameters.gain);
      }
      if (!cal.hasQuantumEfficiency()) {
        cal.setQuantumEfficiency(simulationParameters.qe);
      }
      if (!cal.hasReadNoise()) {
        cal.setReadNoise(simulationParameters.readNoise);
      }
      if (!cal.hasCameraType()) {
        cal.setCameraType(simulationParameters.cameraType);
      }
      if (!cal.hasNmPerPixel()) {
        cal.setNmPerPixel(simulationParameters.pixelPitch);
      }
      if (!cal.hasCameraModelName()) {
        cal.setCameraModelName(simulationParameters.cameraModelName);
      }
      lastCameraBounds = simulationParameters.cameraBounds;
    }

    // Show a dialog to confirm settings
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    final StringBuilder sb = new StringBuilder();
    sb.append("Results contain ").append(TextUtils.pleural(molecules, "molecule")).append('\n');
    sb.append("Min signal = ").append(MathUtils.rounded(minSignal)).append(iUnits).append('\n');
    sb.append("Max signal = ").append(MathUtils.rounded(maxSignal)).append(iUnits).append('\n');
    sb.append("Av signal = ").append(MathUtils.rounded(signalPerFrame)).append(iUnits).append('\n');
    if (fixedDepth) {
      sb.append("Fixed depth = ").append(MathUtils.rounded(depth)).append(zUnits).append('\n');
    }
    gd.addMessage(sb.toString());

    gd.addCheckbox("Flourophore_simulation", fullSimulation);
    gd.addNumericField("Gaussian_SD", sd, 3, 8, "nm");
    gd.addNumericField("Pixel_pitch", cal.getNmPerPixel(), 3, 8, "nm");
    gd.addNumericField("Background", background, 3, 8, "photon");

    // Camera type does not need the full simulation settings. Plus the units are different
    // so just re-implement.
    gd.addChoice("Camera_type", SettingsManager.getCameraTypeNames(),
        CalibrationProtosHelper.getName(settings.getCameraType()), new OptionListener<Integer>() {
          @Override
          public boolean collectOptions(Integer field) {
            settings.setCameraType(SettingsManager.getCameraTypeValues()[field]);
            return collectOptions(false);
          }

          @Override
          public boolean collectOptions() {
            return collectOptions(true);
          }

          private boolean collectOptions(boolean silent) {
            final CameraType cameraType = settings.getCameraType();
            final boolean isCcd = CalibrationProtosHelper.isCcdCameraType(cameraType);
            final ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE, null);
            if (isCcd) {
              egd.addNumericField("Total_gain", cal.getCountPerPhoton(), 3, 8, "count/photon");
              egd.addNumericField("Quantum_efficiency", cal.getQuantumEfficiency(), 3, 8,
                  "e-/photon");
              egd.addNumericField("Read_noise", cal.getReadNoise(), 3, 8, "count");
              egd.addNumericField("Bias", cal.getBias(), 3, 8, "count");
            } else if (cameraType == CameraType.SCMOS) {
              final String[] models = CameraModelManager.listCameraModels(true);
              egd.addChoice("Camera_model_name", models, cal.getCameraModelName());
              egd.addNumericField("Quantum_efficiency", cal.getQuantumEfficiency(), 2, 6,
                  "electron/photon");
            } else {
              IJ.error("Unsupported camera type " + CalibrationProtosHelper.getName(cameraType));
              return false;
            }
            egd.setSilent(silent);
            egd.showDialog(true, gd);
            if (egd.wasCanceled()) {
              return false;
            }
            if (isCcd) {
              cal.setCountPerPhoton(egd.getNextNumber());
              cal.setQuantumEfficiency(egd.getNextNumber());
              cal.setReadNoise(egd.getNextNumber());
              cal.setBias(egd.getNextNumber());
            } else if (cameraType == CameraType.SCMOS) {
              cal.setCameraModelName(egd.getNextChoice());
              cal.setQuantumEfficiency(Math.abs(egd.getNextNumber()));
            }
            return true;
          }
        });

    if (!fixedDepth) {
      gd.addNumericField("Depth", depth, 3, 8, "pixel");
    }

    gd.showDialog();
    if (gd.wasCanceled()) {
      return null;
    }

    fullSimulation = gd.getNextBoolean();
    sd = gd.getNextNumber();
    cal.setNmPerPixel(gd.getNextNumber());
    background = gd.getNextNumber();
    settings.setCameraType(SettingsManager.getCameraTypeValues()[gd.getNextChoiceIndex()]);

    float myDepth = depth;
    if (!fixedDepth) {
      myDepth = (float) gd.getNextNumber();
      if (myDepth < depth) {
        IJ.error(TITLE,
            String.format("Input depth is smaller than the depth guessed from the data: %f < %f",
                myDepth, depth));
        return null;
      }
      depth = myDepth;
    }

    gd.collectOptions();

    // Validate settings
    Rectangle modelBounds = null;
    try {
      ParameterUtils.isAboveZero("Gaussian SD", sd);
      ParameterUtils.isAboveZero("Pixel pitch", cal.getNmPerPixel());
      ParameterUtils.isPositive("Background", background);

      ParameterUtils.isAboveZero("Quantum efficiency", cal.getQuantumEfficiency());
      ParameterUtils.isEqualOrBelow("Quantum efficiency", cal.getQuantumEfficiency(), 1);

      if (cal.isCcdCamera()) {
        ParameterUtils.isAboveZero("Total gain", cal.getCountPerPhoton());
        ParameterUtils.isPositive("Read noise", cal.getReadNoise());
        ParameterUtils.isPositive("Bias", cal.getBias());
      } else if (cal.isScmos()) {
        // Load the model
        cameraModel = CameraModelManager.load(cal.getCameraModelName());
        if (cameraModel == null) {
          IJ.error(TITLE, "Unknown camera model for name: " + cal.getCameraModelName());
          return null;
        }

        int ox = 0;
        int oy = 0;
        if (lastCameraBounds != null) {
          ox = lastCameraBounds.x;
          oy = lastCameraBounds.y;
        }
        cameraModel = PeakFit.cropCameraModel(cameraModel,
            new Rectangle(ox, oy, imp.getWidth(), imp.getHeight()), null, false);
        modelBounds = cameraModel.getBounds();

        final IJImageSource imageSource = (IJImageSource) results.getSource();
        imageSource.setOrigin(modelBounds.x, modelBounds.y);

        cal.clearGlobalCameraSettings();
      } else {
        IJ.error(TITLE, "Unknown camera type: " + cal.getCameraType());
        return null;
      }
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return null;
    }

    // Store calibration
    results.setCalibration(cal.getCalibration());

    final double a = cal.getNmPerPixel();
    final double bias = cal.getBias();
    final double gain = cal.getCountPerPhoton();
    final double readNoise = cal.getReadNoise();
    final double qe = cal.getQuantumEfficiency();

    // Note: The calibration will throw an exception if the converter cannot be created.
    // This is OK as the data will be invalid.

    // Convert +/- depth to total depth in nm
    depth = cal.getDistanceConverter(DistanceUnit.NM).convert(depth * 2);

    // Compute total background variance in photons
    final double backgroundVariance = background;
    // Do not add EM-CCD noise factor. The Mortensen formula also includes this factor
    // so this is "double-counting" the EM-CCD.
    // if (emCCD)
    // backgroundVariance *= 2;

    // Read noise is in ADUs. Convert to Photons to get contribution to background variance
    final double readNoiseInPhotons = readNoise / gain;

    // Get the expected value at each pixel in photons. Assuming a Poisson distribution this
    // is equal to the total variance at the pixel.
    final double b2 = backgroundVariance + readNoiseInPhotons * readNoiseInPhotons;

    // Convert values to photons
    final TypeConverter<IntensityUnit> ic = cal.getIntensityConverter(IntensityUnit.PHOTON);

    final SimulationParameters p = new SimulationParameters(molecules, fullSimulation, sd, a,
        ic.convert(minSignal), ic.convert(maxSignal), ic.convert(signalPerFrame), depth, fixedDepth,
        bias, gain, qe, readNoise, cal.getCameraType(), cal.getCameraModelName(), modelBounds,
        background, b2, createPsf(sd / a));
    p.loaded = true;
    return p;
  }

  /**
   * Show a dialog allowing the parameters for a benchmark simulation to be loaded.
   *
   * @return True if the parameters were collected
   */
  private static boolean showLoadDialog() {
    final String[] images = ImageJUtils.getImageList(ImageJUtils.GREY_SCALE);
    if (images.length == 0) {
      IJ.error(TITLE, "No greyscale benchmark images");
      return false;
    }

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addChoice("Image", images, benchmarkImage);
    gd.addFilenameField("Results_file", benchmarkFile);
    gd.addMessage(TextUtils
        .wrap("Specify if the results are preprocessed. This is true only if the simulation was "
            + "previously loaded and then saved to a GDSC SMLM file format from memory. Set to "
            + "false to load using a universal results loader.", 80));
    gd.addCheckbox("Preprocessed_results", benchmarkAuto);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    benchmarkImage = gd.getNextChoice();
    benchmarkFile = gd.getNextString();
    benchmarkAuto = gd.getNextBoolean();

    return true;
  }

  /**
   * Gets the camera model for processing frames from the simulation image. The model will have
   * bounds that match the simulation image dimensions.
   *
   * @param parameters the parameters
   * @return the camera model
   * @throws ConfigurationException If the model cannot be created
   */
  public static CameraModel getCameraModel(BaseParameters parameters) {
    // Create the camera model
    switch (parameters.cameraType) {
      case CCD:
        return new CcdCameraModel(parameters.bias, parameters.gain);
      case EMCCD:
        return new EmCcdCameraModel(parameters.bias, parameters.gain);

      case SCMOS:
        CameraModel cameraModel = CameraModelManager.load(parameters.cameraModelName);
        if (cameraModel == null) {
          throw new ConfigurationException(
              "Unknown camera model for name: " + parameters.cameraModelName);
        }
        try {
          cameraModel = cameraModel.crop(parameters.cameraBounds, true);
        } catch (final IllegalArgumentException ex) {
          throw new ConfigurationException(ex);
        }
        return cameraModel;

      case UNRECOGNIZED:
      case CAMERA_TYPE_NA:
      default:
        throw new ConfigurationException("Unknown camera model");
    }
  }

  /**
   * Adds the camera description.
   *
   * @param sb the string builder
   * @param parameters the parameters
   * @return the string builder
   * @throws ConfigurationException the configuration exception
   */
  public static StringBuilder addCameraDescription(StringBuilder sb, BaseParameters parameters) {
    if (parameters.cameraType == CameraType.SCMOS) {
      sb.append("sCMOS (").append(parameters.cameraModelName).append(") ");
      final Rectangle bounds = parameters.cameraBounds;
      sb.append(" ").append(bounds.x).append(",").append(bounds.y);
      sb.append(" ").append(bounds.width).append("x").append(bounds.height);
    } else if (CalibrationProtosHelper.isCcdCameraType(parameters.cameraType)) {
      sb.append(CalibrationProtosHelper.getName(parameters.cameraType));
      sb.append(" G=").append(parameters.gain);
      sb.append(" RN=").append(parameters.readNoise);
      sb.append(" B=").append(parameters.bias);
    } else {
      throw new IllegalStateException();
    }
    return sb;
  }
}
