package gdsc.smlm.ij.plugins.pcpalm;

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

import gdsc.smlm.function.SkewNormalFunction;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.ij.IJTrackProgress;
import gdsc.smlm.ij.plugins.About;
import gdsc.smlm.ij.plugins.Parameters;
import gdsc.smlm.ij.plugins.ResultsManager;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.model.MaskDistribution;
import gdsc.smlm.model.StandardFluorophoreSequenceModel;
import gdsc.smlm.model.UniformDistribution;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.NullSource;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.Trace;
import gdsc.smlm.results.TraceManager;
import gdsc.smlm.results.clustering.Cluster;
import gdsc.smlm.results.clustering.ClusterPoint;
import gdsc.smlm.results.clustering.ClusteringAlgorithm;
import gdsc.smlm.results.clustering.ClusteringEngine;
import gdsc.smlm.utils.Maths;
import gdsc.smlm.utils.Statistics;
import gdsc.smlm.utils.StoredDataStatistics;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.plugin.frame.Recorder;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.MultivariateMatrixFunction;
import org.apache.commons.math3.analysis.MultivariateVectorFunction;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.fitting.GaussianFitter;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.PointVectorValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.optim.nonlinear.vector.ModelFunction;
import org.apache.commons.math3.optim.nonlinear.vector.ModelFunctionJacobian;
import org.apache.commons.math3.optim.nonlinear.vector.MultivariateVectorOptimizer;
import org.apache.commons.math3.optim.nonlinear.vector.Target;
import org.apache.commons.math3.optim.nonlinear.vector.Weight;
import org.apache.commons.math3.optim.nonlinear.vector.jacobian.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.FastMath;

/**
 * Use the PC-PALM protocol to prepare a set of localisations into molecules. This can be used for for clustering
 * analysis.
 * <p>
 * See Sengupta, et al (2013). Quantifying spatial resolution in point-localisation superresolution images using pair
 * correlation analysis. Nature Protocols 8, pp345-354.
 * <p>
 * See also Veatch, et al (2012). Correlation Functions Quantify Super-Resolution Images and Estimate Apparent
 * Clustering Due to Over-Counting. PLoS One 7, Issue 2, e31457
 */
public class PCPALMMolecules implements PlugIn
{
	static String TITLE = "PC-PALM Molecules";
	private static String inputOption = "";
	private static boolean chooseRoi = false;
	private static double nmPerPixelLimit = 0;

	private static String[] RUN_MODE = { "PC-PALM", "Manual Tracing", "In-memory results", "Simulation" };
	private static int runMode = 0;

	// Mode 0: PC-PALM protocol for estimating localisation precision and then tracing molecules
	private static String roiImage = "";
	private static int histogramBins = 50;
	private static String[] singlesMode = new String[] { "Ignore", "Include in molecules histogram",
			"Include in final filtering" };
	private static int singlesModeIndex = 1;
	private static boolean simplexFitting = false;
	private static boolean showHistograms = true;
	private static boolean binaryImage = true;
	static double blinkingRate = 2;
	private static double p = 0.6;
	private static String[] BLINKING_DISTRIBUTION = new String[] { "Poisson", "Geometric", "None", "Binomial" };
	private static int blinkingDistribution = 0;
	private static boolean clearResults = false;

	// Mode 1. Manual tracing of molecules
	private static double dThreshold = 150;
	private static double tThreshold = 1;

	// Mode 2. Direct use of in-memory results
	// - No parameters needed

	// Mode 3. Random simulation of molecules
	private static int nMolecules = 2000;
	private static double simulationSize = 16;
	private static boolean distanceAnalysis = false;

	private static String[] CLUSTER_SIMULATION = new String[] { "None", "Circles", "Non-overlapping circles",
			"Circles Mask" };
	private static int clusterSimulation = 0;
	private static double clusterNumber = 3;
	private static double clusterNumberSD = 0;
	private static double clusterRadius = 50;
	private static boolean showClusterMask = false;

	// Low resolution image construction
	private static int lowResolutionImageSize = 1024;
	private static double roiSizeInUm = 4;
	private static boolean showHighResolutionImage = false;

	private Rectangle roiBounds;
	private int roiImageWidth, roiImageHeight;
	private long start;

	// These package level variables are used by the PCPALMAnalysis plugin.
	static MemoryPeakResults results;
	static double minx, miny, maxx, maxy;
	static double nmPerPixel;
	static ArrayList<Molecule> molecules = null;
	static double sigmaS = 20;
	static double densityPeaks;
	static double densityProtein;
	static double seconds;
	static double area;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		// Require some fit results and selected regions
		boolean resultsAvailable = MemoryPeakResults.countMemorySize() > 0;

		if (!getRunMode(resultsAvailable))
			return;

		if (runMode != 3)
		{
			results = ResultsManager.loadInputResults(inputOption, true);
			if (results == null || results.size() == 0)
			{
				IJ.error(TITLE, "No results could be loaded");
				return;
			}

			if (results.getCalibration() == null)
			{
				IJ.error(TITLE, "Results are not calibrated");
				return;
			}

			// Get the lifetime before cropping as this is the true representation of the number of frames.
			// This may truncate the lifetime if the first/last localisation are not near the end of the
			// acquisition lifetime
			getLifetime();

			results = cropToRoi(results);
			if (results.size() == 0)
			{
				IJ.error(TITLE, "No results within the crop region");
				return;
			}
		}

		// Clear cached results
		molecules = null;

		// Different run-modes for generating the set of molecules for analysis
		switch (runMode)
		{
			case 0:
				runPCPALM();
				break;
			case 1:
				runManualTracing();
				break;
			case 2:
				runInMemoryResults();
				break;
			case 3:
				runSimulation(resultsAvailable);
				area = simulationSize * simulationSize;
				seconds = 100; // Use an arbitrary lifetime
				break;
		}

		if (molecules == null)
			return;
		if (molecules.size() < 2)
		{
			IJ.error(TITLE, "Not enough molecules to construct a binary image");
			return;
		}

		// Generate binary PALM image
		if (!createImage(molecules))
			return;

		// Density is required for the PC analysis
		densityPeaks = calculatePeakDensity();
		if (runMode == 0 || runMode == 3)
		{
			// Blinking rate is mentioned in the PC-PALM protocol and so we include it here.
			// TODO - Add automated estimation of the blinking rate from the data using the method of 
			// Annibale, et al (2011), Quantitative photo activated localization microscopy: unraveling the 
			// effects of photoblinking. PLoS One, 6(7): e22678 (http://dx.doi.org/10.1371%2Fjournal.pone.0022678)
			densityProtein = densityPeaks / blinkingRate;
			log("Peak Density = %s (um^-2). Protein Density = %s (um^-2)", Utils.rounded(densityPeaks * 1e6),
					Utils.rounded(densityProtein * 1e6));
		}
		else
		{
			// No blinking rate for non PC-PALM methods. This can be configured in later plugins if required.
			blinkingRate = 1;
			densityProtein = densityPeaks;
			log("Molecule Density = %s (um^-2)", Utils.rounded(densityPeaks * 1e6));
		}

		log("Results lifetime = %s s", Utils.rounded(seconds));

		// Use a second plugin filter that will work on a region drawn on the binary image
		// and compute the PALM analysis

		double seconds = (System.currentTimeMillis() - start) / 1000.0;
		String msg = TITLE + " complete : " + seconds + "s";
		IJ.showStatus(msg);
		log(msg);
	}

	private boolean getRunMode(boolean resultsAvailable)
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		// Build a list of all images with a region ROI
		List<String> titles = new LinkedList<String>();
		if (WindowManager.getWindowCount() > 0)
		{
			for (int imageID : WindowManager.getIDList())
			{
				ImagePlus imp = WindowManager.getImage(imageID);
				if (imp != null && imp.getRoi() != null && imp.getRoi().isArea())
					titles.add(imp.getTitle());
			}
		}

		if (!resultsAvailable)
		{
			runMode = 3;

			gd.addMessage("Simulate molecules for cluster analysis.\nComputes a binary image from localisation data");

			gd.addNumericField("Molecules", nMolecules, 0);
			gd.addNumericField("Simulation_size (um)", simulationSize, 2);
			gd.addNumericField("Blinking_rate", blinkingRate, 2);
			gd.addChoice("Blinking_distribution", BLINKING_DISTRIBUTION, BLINKING_DISTRIBUTION[blinkingDistribution]);
			gd.addNumericField("Average_precision (nm)", sigmaS, 2);
			gd.addCheckbox("Show_histograms", showHistograms);
			gd.addCheckbox("Distance_analysis", distanceAnalysis);
			gd.addChoice("Cluster_simulation", CLUSTER_SIMULATION, CLUSTER_SIMULATION[clusterSimulation]);
			gd.addNumericField("Cluster_number", clusterNumber, 2);
			gd.addNumericField("Cluster_variation (SD)", clusterNumberSD, 2);
			gd.addNumericField("Cluster_radius", clusterRadius, 2);
			gd.addCheckbox("Show_cluster_mask", showClusterMask);

			Recorder.recordOption("Run_mode", RUN_MODE[runMode]);
		}
		else
		{
			gd.addMessage("Prepare molecules for cluster analysis.\nComputes a binary image from raw localisation data");
			ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);
			if (!titles.isEmpty())
				gd.addCheckbox((titles.size() == 1) ? "Use_ROI" : "Choose_ROI", chooseRoi);
			gd.addChoice("Run_mode", RUN_MODE, RUN_MODE[runMode]);
		}

		gd.addMessage("Select options for low resolution image:");
		gd.addSlider("Image_size (px)", 512, 2048, lowResolutionImageSize);
		gd.addSlider("ROI_size (um)", 1.5, 4, roiSizeInUm);
		gd.addMessage("Select options for high resolution image:");
		gd.addCheckbox("Show_high_res_image", showHighResolutionImage);
		gd.addSlider("nm_per_pixel_limit", 0, 20, nmPerPixelLimit);
		gd.addMessage("Optionally remove all analysis results from memory");
		gd.addCheckbox("Clear_results", clearResults);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		if (!resultsAvailable)
		{
			nMolecules = (int) Math.abs(gd.getNextNumber());
			simulationSize = Math.abs(gd.getNextNumber());
			blinkingRate = Math.abs(gd.getNextNumber());
			blinkingDistribution = gd.getNextChoiceIndex();
			sigmaS = Math.abs(gd.getNextNumber());
			showHistograms = gd.getNextBoolean();
			distanceAnalysis = gd.getNextBoolean();
			clusterSimulation = gd.getNextChoiceIndex();
			clusterNumber = Math.abs(gd.getNextNumber());
			clusterNumberSD = Math.abs(gd.getNextNumber());
			clusterRadius = Math.abs(gd.getNextNumber());
			showClusterMask = gd.getNextBoolean();
		}
		else
		{
			inputOption = ResultsManager.getInputSource(gd);
			if (!titles.isEmpty())
				chooseRoi = gd.getNextBoolean();
			runMode = gd.getNextChoiceIndex();
		}
		lowResolutionImageSize = (int) gd.getNextNumber();
		roiSizeInUm = gd.getNextNumber();
		showHighResolutionImage = gd.getNextBoolean();
		nmPerPixelLimit = Math.abs(gd.getNextNumber());
		clearResults = gd.getNextBoolean();

		// Check arguments
		try
		{
			if (!resultsAvailable)
			{
				Parameters.isAboveZero("Molecules", nMolecules);
				Parameters.isAboveZero("Simulation size", simulationSize);
				Parameters.isEqualOrAbove("Blinking rate", blinkingRate, 1);
				Parameters.isEqualOrAbove("Cluster number", clusterNumber, 1);
			}
			Parameters.isAbove("Image scale", lowResolutionImageSize, 1);
			Parameters.isAboveZero("ROI size", roiSizeInUm);
		}
		catch (IllegalArgumentException ex)
		{
			IJ.error(TITLE, ex.getMessage());
			return false;
		}

		if (!titles.isEmpty() && chooseRoi && resultsAvailable)
		{
			if (titles.size() == 1)
			{
				roiImage = titles.get(0);
				Recorder.recordOption("Image", roiImage);
			}
			else
			{
				String[] items = titles.toArray(new String[titles.size()]);
				gd = new GenericDialog(TITLE);
				gd.addMessage("Select the source image for the ROI");
				gd.addChoice("Image", items, roiImage);
				gd.showDialog();
				if (gd.wasCanceled())
					return false;
				roiImage = gd.getNextChoice();
			}
			ImagePlus imp = WindowManager.getImage(roiImage);

			roiBounds = imp.getRoi().getBounds();
			roiImageWidth = imp.getWidth();
			roiImageHeight = imp.getHeight();
		}
		else
		{
			roiBounds = null;
		}

		if (!resultsAvailable)
		{
			if (!getPValue())
				return false;
		}

		if (clearResults)
		{
			PCPALMAnalysis.results.clear();
			PCPALMFitting.previous_gr = null;
		}

		return true;
	}

	private MemoryPeakResults cropToRoi(MemoryPeakResults results)
	{
		Rectangle bounds = results.getBounds(true);
		area = (bounds.width * bounds.height * results.getNmPerPixel() * results.getNmPerPixel()) / 1e6;
		if (roiBounds == null)
		{
			return results;
		}

		// Adjust bounds relative to input results image
		double xscale = (double) roiImageWidth / bounds.width;
		double yscale = (double) roiImageHeight / bounds.height;
		roiBounds.x /= xscale;
		roiBounds.width /= xscale;
		roiBounds.y /= yscale;
		roiBounds.height /= yscale;

		float minX = (int) (roiBounds.x);
		float maxX = (int) Math.ceil(roiBounds.x + roiBounds.width);
		float minY = (int) (roiBounds.y);
		float maxY = (int) Math.ceil(roiBounds.y + roiBounds.height);

		// Update the area with the cropped region
		area *= (maxX - minX) / bounds.width;
		area *= (maxY - minY) / bounds.height;

		// Create a new set of results within the bounds
		MemoryPeakResults newResults = new MemoryPeakResults();
		newResults.begin();
		for (PeakResult peakResult : results.getResults())
		{
			float x = peakResult.params[Gaussian2DFunction.X_POSITION];
			float y = peakResult.params[Gaussian2DFunction.Y_POSITION];
			if (x < minX || x > maxX || y < minY || y > maxY)
				continue;
			newResults.add(peakResult);
		}
		newResults.end();
		newResults.copySettings(results);
		newResults.setBounds(new Rectangle((int) minX, (int) minY, (int) (maxX - minX), (int) (maxY - minY)));
		return newResults;
	}

	private void runPCPALM()
	{
		if (!showPCPALMDialog())
			return;

		startLog();

		// Follow the PC-PALM protocol
		log("Fitting localisation precision...");
		ArrayList<Molecule> localisations = extractLocalisations(results);
		double sigmaRaw = calculateAveragePrecision(localisations, "Localisations");
		log("%d localisations with an average precision of %.2f", results.size(), sigmaRaw);

		log("Fitting molecule precision...");
		ArrayList<Molecule> singles = new ArrayList<Molecule>();
		molecules = extractMolecules(results, sigmaRaw, singles);
		if (singlesModeIndex == 1)
			molecules.addAll(singles);
		sigmaS = calculateAveragePrecision(molecules, "Molecules");
		log("%d molecules with an average precision of %.2f", molecules.size(), sigmaS);

		// Q. Should this filter the original localisations or just the grouped peaks?
		if (singlesModeIndex == 2)
			molecules.addAll(singles);
		molecules = filterMolecules(molecules, sigmaS);
		log("%d molecules within precision %.2f", molecules.size(), 3 * sigmaS);
	}

	private void startLog()
	{
		logSpacer();
		log(TITLE);
		logSpacer();
		start = System.currentTimeMillis();
	}

	private boolean showPCPALMDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Estimate the average localisation precision by fitting histograms.\nUse the precision to trace localisations into molecule pulses.");

		gd.addNumericField("Histogram_bins", histogramBins, 0);
		gd.addChoice("Singles_mode", singlesMode, singlesMode[singlesModeIndex]);
		gd.addCheckbox("Simplex_fit", simplexFitting);
		gd.addCheckbox("Show_histograms", showHistograms);
		gd.addCheckbox("Binary_image", binaryImage);
		gd.addNumericField("Blinking_rate", blinkingRate, 2);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		histogramBins = (int) gd.getNextNumber();
		singlesModeIndex = gd.getNextChoiceIndex();
		simplexFitting = gd.getNextBoolean();
		showHistograms = gd.getNextBoolean();
		binaryImage = gd.getNextBoolean();
		blinkingRate = gd.getNextNumber();

		// Check arguments
		try
		{
			Parameters.isAbove("Histogram bins", histogramBins, 1);
			Parameters.isEqualOrAbove("Blinking rate", blinkingRate, 1);
		}
		catch (IllegalArgumentException ex)
		{
			IJ.error(TITLE, ex.getMessage());
			return false;
		}

		return true;
	}

	/**
	 * Extract molecules for the PC-PALM analysis.
	 * <p>
	 * Estimate the localisation uncertainty (precision) of each molecule using the formula of Mortensen, et al (2010),
	 * Nature Methods 7, 377-381. Store distance in nm and signal in photons using the calibration
	 * 
	 * @param results
	 * @return
	 */
	public ArrayList<Molecule> extractLocalisations(MemoryPeakResults results)
	{
		ArrayList<Molecule> molecules = new ArrayList<Molecule>(results.size());
		final double nmPerPixel = results.getNmPerPixel();
		final double gain = results.getGain();
		final boolean emCCD = results.isEMCCD();
		for (PeakResult r : results.getResults())
		{
			double p = r.getPrecision(nmPerPixel, gain, emCCD);
			// Remove EMCCD adjustment
			//p /= Math.sqrt(PeakResult.F);
			molecules.add(new Molecule(r.getXPosition() * nmPerPixel, r.getYPosition() * nmPerPixel, p, r.getSignal() /
					gain));
		}
		return molecules;
	}

	/**
	 * Calculate the average precision by fitting a skewed Gaussian to the histogram of the precision distribution.
	 * 
	 * @param molecules
	 * @param subTitle
	 * @return The average precision
	 */
	private double calculateAveragePrecision(ArrayList<Molecule> molecules, String subTitle)
	{
		String title = (showHistograms) ? TITLE + " Histogram " + subTitle : null;
		return calculateAveragePrecision(molecules, title, histogramBins, true, true);
	}

	/**
	 * Calculate the average precision by fitting a skewed Gaussian to the histogram of the precision distribution.
	 * <p>
	 * A simple mean and SD of the histogram is computed. If the mean of the Skewed Gaussian does not fit within 3 SDs
	 * of the simple mean then the simple mean is returned.
	 * 
	 * @param molecules
	 * @param title
	 *            the plot title (null if no plot should be displayed)
	 * @param histogramBins
	 * @param logFitParameters
	 *            Record the fit parameters to the ImageJ log
	 * @param removeOutliers
	 *            The distribution is created using all values within 1.5x the inter-quartile range (IQR) of the data
	 * @return The average precision
	 */
	public double calculateAveragePrecision(ArrayList<Molecule> molecules, String title, int histogramBins,
			boolean logFitParameters, boolean removeOutliers)
	{
		// Plot histogram of the precision
		float[] data = new float[molecules.size()];
		DescriptiveStatistics stats = new DescriptiveStatistics();
		double yMin = Double.NEGATIVE_INFINITY, yMax = 0;
		for (int i = 0; i < data.length; i++)
		{
			data[i] = (float) molecules.get(i).precision;
			stats.addValue(data[i]);
		}

		// Set the min and max y-values using 1.5 x IQR 
		if (removeOutliers)
		{
			double lower = stats.getPercentile(25);
			double upper = stats.getPercentile(75);
			if (Double.isNaN(lower) || Double.isNaN(upper))
			{
				if (logFitParameters)
					Utils.log("Error computing IQR: %f - %f", lower, upper);
			}
			else
			{
				double iqr = upper - lower;

				yMin = FastMath.max(lower - iqr, stats.getMin());
				yMax = FastMath.min(upper + iqr, stats.getMax());

				if (logFitParameters)
					Utils.log("  Data range: %f - %f. Plotting 1.5x IQR: %f - %f", stats.getMin(), stats.getMax(),
							yMin, yMax);
			}
		}

		if (yMin == Double.NEGATIVE_INFINITY)
		{
			yMin = stats.getMin();
			yMax = stats.getMax();

			if (logFitParameters)
				Utils.log("  Data range: %f - %f", yMin, yMax);
		}

		float[][] hist = Utils.calcHistogram(data, yMin, yMax, histogramBins);

		float[] xValues = Utils.createHistogramAxis(hist[0]);
		float[] yValues = Utils.createHistogramValues(hist[1]);

		Plot plot = null;
		if (title != null)
		{
			plot = new Plot(title, "Precision", "Frequency", xValues, yValues);
			Utils.display(title, plot);
		}

		// Extract non-zero data
		float[] x = Arrays.copyOf(hist[0], hist[0].length);
		float[] y = hist[1];
		int count = 0;
		float dx = (x[1] - x[0]) * 0.5f;
		for (int i = 0; i < y.length; i++)
			if (y[i] > 0)
			{
				x[count] = x[i] + dx;
				y[count] = y[i];
				count++;
			}
		x = Arrays.copyOf(x, count);
		y = Arrays.copyOf(y, count);

		// Sense check to fitted data. Get mean and SD of histogram
		double[] stats2 = Utils.getHistogramStatistics(x, y);
		double mean = stats2[0];
		if (logFitParameters)
			log("  Initial Statistics: %f +/- %f", stats2[0], stats2[1]);

		// Standard Gaussian fit
		double[] parameters = fitGaussian(x, y);
		if (parameters == null)
		{
			log("  Failed to fit initial Gaussian");
			return mean;
		}
		double newMean = parameters[1];
		double error = Math.abs(stats2[0] - newMean) / stats2[1];
		if (error > 3)
		{
			log("  Failed to fit Gaussian: %f standard deviations from histogram mean", error);
			return mean;
		}
		if (newMean < yMin || newMean > yMax)
		{
			log("  Failed to fit Gaussian: %f outside data range %f - %f", newMean, yMin, yMax);
			return mean;
		}

		mean = newMean;

		if (logFitParameters)
			log("  Initial Gaussian: %f @ %f +/- %f", parameters[0], parameters[1], parameters[2]);

		double[] initialSolution = new double[] { parameters[0], parameters[1], parameters[2], -1 };

		// Fit to a skewed Gaussian (or appropriate function)
		double[] skewParameters = fitSkewGaussian(x, y, initialSolution);
		if (skewParameters == null)
		{
			log("  Failed to fit Skewed Gaussian");
			return mean;
		}

		SkewNormalFunction sn = new SkewNormalFunction(skewParameters);
		if (logFitParameters)
			log("  Skewed Gaussian: %f @ %f +/- %f (a = %f) => %f +/- %f", skewParameters[0], skewParameters[1],
					skewParameters[2], skewParameters[3], sn.getMean(), Math.sqrt(sn.getVariance()));

		newMean = sn.getMean();
		error = Math.abs(stats2[0] - newMean) / stats2[1];
		if (error > 3)
		{
			log("  Failed to fit Skewed Gaussian: %f standard deviations from histogram mean", error);
			return mean;
		}
		if (newMean < yMin || newMean > yMax)
		{
			log("  Failed to fit Skewed Gaussian: %f outside data range %f - %f", newMean, yMin, yMax);
			return mean;
		}

		// Use original histogram x-axis to maintain all the bins
		if (plot != null)
		{
			x = hist[0];
			for (int i = 0; i < y.length; i++)
				x[i] += dx;
			plot.setColor(Color.red);
			addToPlot(plot, x, skewParameters, Plot.LINE);

			plot.setColor(Color.black);
			Utils.display(title, plot);
		}

		// Return the average precision from the fitted curve
		return newMean;
	}

	private double[] fitGaussian(float[] x, float[] y)
	{
		MyGaussianFitter fitter = new MyGaussianFitter(
				new org.apache.commons.math3.optim.nonlinear.vector.jacobian.LevenbergMarquardtOptimizer(), 2000);
		for (int i = 0; i < x.length; i++)
			fitter.addObservedPoint(x[i], y[i]);

		double[] parameters;
		try
		{
			parameters = fitter.fit();
		}
		catch (TooManyEvaluationsException e)
		{
			// Use the initial estimate
			parameters = fitter.guess;
		}
		catch (Exception e)
		{
			// Just in case there is another exception type, or the initial estimate failed
			return null;
		}
		return parameters;
	}

	private double[] fitSkewGaussian(float[] x, float[] y, double[] initialSolution)
	{
		try
		{
			double[] skewParameters = (simplexFitting) ? optimiseSimplex(x, y, initialSolution) : optimiseLeastSquares(
					x, y, initialSolution);
			return skewParameters;
		}
		catch (TooManyEvaluationsException e)
		{
			return null;
		}
	}

	private double[] optimiseLeastSquares(float[] x, float[] y, double[] initialSolution)
	{
		// Least-squares optimisation using numerical gradients
		final SkewNormalDifferentiableFunction myFunction = new SkewNormalDifferentiableFunction(initialSolution);
		myFunction.addData(x, y);
		LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer();

		PointVectorValuePair optimum = optimizer.optimize(new MaxIter(3000), new MaxEval(Integer.MAX_VALUE),
				new ModelFunctionJacobian(new MultivariateMatrixFunction()
				{
					public double[][] value(double[] point) throws IllegalArgumentException
					{
						return myFunction.jacobian(point);
					}
				}), new ModelFunction(myFunction), new Target(myFunction.calculateTarget()),
				new Weight(myFunction.calculateWeights()), new InitialGuess(initialSolution));

		double[] skewParameters = optimum.getPoint();
		return skewParameters;
	}

	private double[] optimiseSimplex(float[] x, float[] y, double[] initialSolution)
	{
		// Simplex optimisation
		SkewNormalMultivariateFunction sn2 = new SkewNormalMultivariateFunction(initialSolution);
		sn2.addData(x, y);
		NelderMeadSimplex simplex = new NelderMeadSimplex(4);
		SimplexOptimizer opt = new SimplexOptimizer(1e-6, 1e-10);
		PointValuePair solution = opt.optimize(new MaxEval(1000), new InitialGuess(initialSolution), simplex,
				new ObjectiveFunction(sn2), GoalType.MINIMIZE);

		double[] skewParameters2 = solution.getPointRef();
		return skewParameters2;
	}

	/**
	 * Add the skewed gaussian to the histogram plot
	 * 
	 * @param plot
	 * @param x
	 * @param parameters
	 *            Gaussian parameters
	 * @param alpha
	 * @param shape
	 */
	private void addToPlot(Plot plot, float[] x, double[] parameters, int shape)
	{
		SkewNormalFunction sn = new SkewNormalFunction(parameters);
		float[] y = new float[x.length];
		for (int i = 0; i < x.length; i++)
			y[i] = (float) sn.evaluate(x[i]);
		plot.addPoints(x, y, shape);
	}

	/**
	 * Group all localisations in successive frames within 2.5x of the initial precision estimate into a single molecule
	 * 
	 * @param results
	 *            The results
	 * @param sigmaRaw
	 *            The initial precision estimate
	 * @param singles
	 *            a list of the singles (not grouped into molecules)
	 * @return a list of molecules
	 */
	private ArrayList<Molecule> extractMolecules(MemoryPeakResults results, double sigmaRaw, ArrayList<Molecule> singles)
	{
		return traceMolecules(results, sigmaRaw * 2.5, 1, singles);
	}

	/**
	 * Trace localisations
	 * 
	 * @param results
	 *            The results
	 * @param distance
	 *            The distance threshold (nm)
	 * @param time
	 *            The time threshold (frames)
	 * @param singles
	 *            a list of the singles (not grouped into molecules)
	 * @return a list of molecules
	 */
	private ArrayList<Molecule> traceMolecules(MemoryPeakResults results, double distance, int time,
			ArrayList<Molecule> singles)
	{
		TraceManager tm = new TraceManager(results);
		double distanceThreshold = distance / results.getNmPerPixel();
		tm.traceMolecules(distanceThreshold, time);
		Trace[] traces = tm.getTraces();
		ArrayList<Molecule> molecules = new ArrayList<Molecule>(traces.length);
		final double nmPerPixel = results.getNmPerPixel();
		final double gain = results.getGain();
		final boolean emCCD = results.isEMCCD();
		for (Trace t : traces)
		{
			double p = t.getLocalisationPrecision(nmPerPixel, gain, emCCD);
			if (t.size() == 1)
			{
				float[] centroid = t.getCentroid();
				singles.add(new Molecule(centroid[0] * nmPerPixel, centroid[1] * nmPerPixel, p, t.getSignal() / gain));
			}
			else
			{
				float[] centroid = t.getCentroid();
				molecules
						.add(new Molecule(centroid[0] * nmPerPixel, centroid[1] * nmPerPixel, p, t.getSignal() / gain));
			}
		}
		log("  %d localisations traced to %d molecules (%d singles, %d traces) using d=%.2f nm, t=%d frames (%s s)",
				results.size(), molecules.size(), singles.size(), molecules.size() - singles.size(), distance, time,
				Utils.rounded(time * results.getCalibration().exposureTime / 1000.0));
		return molecules;
	}

	/**
	 * Calculate the density of peaks in the original data
	 * 
	 * @return The peak density
	 */
	private double calculatePeakDensity()
	{
		//double pcw = PCPALMMolecules.maxx - PCPALMMolecules.minx;
		//double pch = PCPALMMolecules.maxy - PCPALMMolecules.miny;
		//double area = pcw * pch;
		// Use the area from the source of the molecules
		return molecules.size() / (area * 1E6);
	}

	/**
	 * Return a new list, removing all molecules with a precision over 3x of the precision estimate
	 * 
	 * @param molecules
	 * @param sigmaS
	 *            The precision estimate
	 * @return
	 */
	private ArrayList<Molecule> filterMolecules(ArrayList<Molecule> molecules, double sigmaS)
	{
		ArrayList<Molecule> newMolecules = new ArrayList<Molecule>(molecules.size());
		final double limit = 3 * sigmaS;
		for (Molecule m : molecules)
		{
			if (m.precision <= limit)
				newMolecules.add(m);
		}
		return newMolecules;
	}

	private void runManualTracing()
	{
		if (!showManualTracingDialog())
			return;

		startLog();

		// Convert seconds to frames
		int timeInFrames = FastMath.max(1,
				(int) Math.round(tThreshold * 1000.0 / results.getCalibration().exposureTime));

		ArrayList<Molecule> singles = new ArrayList<Molecule>();
		molecules = traceMolecules(results, dThreshold, timeInFrames, singles);
		molecules.addAll(singles);
	}

	private boolean showManualTracingDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Use distance and time thresholds to trace localisations into molecules.");

		gd.addNumericField("Distance (nm)", dThreshold, 0);
		gd.addNumericField("Time (seconds)", tThreshold, 2);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		dThreshold = Math.abs(gd.getNextNumber());
		tThreshold = Math.abs(gd.getNextNumber());

		// Check arguments
		try
		{
			Parameters.isAboveZero("Distance threshold", dThreshold);
			Parameters.isAboveZero("Time threshold", tThreshold);
		}
		catch (IllegalArgumentException ex)
		{
			IJ.error(TITLE, ex.getMessage());
			return false;
		}

		return true;
	}

	private void runInMemoryResults()
	{
		startLog();
		molecules = extractLocalisations(results);
	}

	private void runSimulation(boolean resultsAvailable)
	{
		if (resultsAvailable && !showSimulationDialog())
			return;

		startLog();

		log("Simulation parameters");
		if (blinkingDistribution == 3)
		{
			log("  - Clusters = %d", nMolecules);
			log("  - Simulation size = %s um", Utils.rounded(simulationSize, 4));
			log("  - Molecules/cluster = %s", Utils.rounded(blinkingRate, 4));
			log("  - Blinking distribution = %s", BLINKING_DISTRIBUTION[blinkingDistribution]);
			log("  - p-Value = %s", Utils.rounded(p, 4));
		}
		else
		{
			log("  - Molecules = %d", nMolecules);
			log("  - Simulation size = %s um", Utils.rounded(simulationSize, 4));
			log("  - Blinking rate = %s", Utils.rounded(blinkingRate, 4));
			log("  - Blinking distribution = %s", BLINKING_DISTRIBUTION[blinkingDistribution]);
		}
		log("  - Average precision = %s nm", Utils.rounded(sigmaS, 4));
		log("  - Clusters simulation = " + CLUSTER_SIMULATION[clusterSimulation]);
		if (clusterSimulation > 0)
		{
			log("  - Cluster number = %s +/- %s", Utils.rounded(clusterNumber, 4), Utils.rounded(clusterNumberSD, 4));
			log("  - Cluster radius = %s nm", Utils.rounded(clusterRadius, 4));
		}

		final double nmPerPixel = 100;
		double width = simulationSize * 1000.0;
		// Allow a border of 3 x sigma for +/- precision
		//if (blinkingRate > 1)
		width -= 3 * sigmaS;
		RandomGenerator randomGenerator = new Well19937c(System.currentTimeMillis() + System.identityHashCode(this));
		RandomDataGenerator dataGenerator = new RandomDataGenerator(randomGenerator);
		UniformDistribution dist = new UniformDistribution(null, new double[] { width, width, 0 },
				randomGenerator.nextInt());

		molecules = new ArrayList<Molecule>(nMolecules);
		// Create some dummy results since the calibration is required for later analysis
		results = new MemoryPeakResults();
		results.setCalibration(new gdsc.smlm.results.Calibration(nmPerPixel, 1, 100));
		results.setSource(new NullSource("Molecule Simulation"));
		results.begin();
		int count = 0;

		// Generate a sequence of coordinates
		ArrayList<double[]> xyz = new ArrayList<double[]>((int) (nMolecules * 1.1));

		Statistics statsRadius = new Statistics();
		Statistics statsSize = new Statistics();
		String maskTitle = TITLE + " Cluster Mask";
		ByteProcessor bp = null;
		double maskScale = 0;

		// TODO - Add a fluctuations model to this.

		if (clusterSimulation > 0)
		{
			// Simulate clusters.

			// Note: In the Veatch et al. paper (Plos 1, e31457) correlation functions are built using circles
			// with small radii of 4-8 Arbitrary Units (AU) or large radii of 10-30 AU. A fluctuations model is
			// created at T = 1.075 Tc. It is not clear exactly how the particles are distributed.
			// It may be that a mask is created first using the model. The particles are placed on the mask using
			// a specified density. This simulation produces a figure to show either a damped cosine function
			// (circles) or an exponential (fluctuations). The number of particles in each circle may be randomly
			// determined just by density. The figure does not discuss the derivation of the cluster size 
			// statistic.
			// 
			// If this plugin simulation is run with a uniform distribution and blinking rate of 1 then the damped
			// cosine function is reproduced. The curve crosses g(r)=1 at a value equivalent to the average
			// distance to the centre-of-mass of each drawn cluster, not the input cluster radius parameter (which 
			// is a hard upper limit on the distance to centre).

			final int maskSize = lowResolutionImageSize;
			int[] mask = null;
			maskScale = width / maskSize; // scale is in nm/pixel

			ArrayList<double[]> clusterCentres = new ArrayList<double[]>();
			int totalSteps = 1 + (int) Math.ceil(nMolecules / clusterNumber);
			if (clusterSimulation == 2 || clusterSimulation == 3)
			{
				// Clusters are non-overlapping circles

				// Ensure the circles do not overlap by using an exclusion mask that accumulates 
				// out-of-bounds pixels by drawing the last cluster (plus some border) on an image. When no
				// more pixels are available then stop generating molecules.
				// This is done by cumulatively filling a mask and using the MaskDistribution to select 
				// a new point. This may be slow but it works.

				// TODO - Allow clusters of different sizes...

				mask = new int[maskSize * maskSize];
				Arrays.fill(mask, 255);
				MaskDistribution maskDistribution = new MaskDistribution(mask, maskSize, maskSize, 0, maskScale,
						maskScale, randomGenerator);
				double[] centre;
				IJ.showStatus("Computing clusters mask");
				int roiRadius = (int) Math.round((clusterRadius * 2) / maskScale);

				if (clusterSimulation == 3)
				{
					// Generate a mask of circles then sample from that.
					// If we want to fill the mask completely then adjust the total steps to be the number of 
					// circles that can fit inside the mask.
					totalSteps = (int) (maskSize * maskSize / (Math.PI * Math.pow(clusterRadius / maskScale, 2)));
				}

				while ((centre = maskDistribution.next()) != null && clusterCentres.size() < totalSteps)
				{
					IJ.showProgress(clusterCentres.size(), totalSteps);
					// The mask returns the coordinates with the centre of the image at 0,0
					centre[0] += width / 2;
					centre[1] += width / 2;
					clusterCentres.add(centre);

					// Fill in the mask around the centre to exclude any more circles that could overlap
					double cx = centre[0] / maskScale;
					double cy = centre[1] / maskScale;
					fillMask(mask, maskSize, (int) cx, (int) cy, roiRadius, 0);
					//log("[%.1f,%.1f] @ [%.1f,%.1f]", centre[0], centre[1], cx, cy);
					//Utils.display("Mask", new ColorProcessor(maskSize, maskSize, mask));
					try
					{
						maskDistribution = new MaskDistribution(mask, maskSize, maskSize, 0, maskScale, maskScale,
								randomGenerator);
					}
					catch (IllegalArgumentException e)
					{
						// This can happen when there are no more non-zero pixels
						log("WARNING: No more room for clusters on the mask area (created %d of estimated %d)",
								clusterCentres.size(), totalSteps);
						break;
					}
				}
				IJ.showProgress(1);
				IJ.showStatus("");
			}
			else
			{
				// Clusters are overlapping circles

				// Pick centres randomly from the distribution 
				while (clusterCentres.size() < totalSteps)
					clusterCentres.add(dist.next());
			}

			if (showClusterMask || clusterSimulation == 3)
			{
				// Show the mask for the clusters
				if (mask == null)
					mask = new int[maskSize * maskSize];
				else
					Arrays.fill(mask, 0);
				int roiRadius = (int) Math.round((clusterRadius) / maskScale);
				for (double[] c : clusterCentres)
				{
					double cx = c[0] / maskScale;
					double cy = c[1] / maskScale;
					fillMask(mask, maskSize, (int) cx, (int) cy, roiRadius, 1);
				}

				if (clusterSimulation == 3)
				{
					// We have the mask. Now pick points at random from the mask.
					MaskDistribution maskDistribution = new MaskDistribution(mask, maskSize, maskSize, 0, maskScale,
							maskScale, randomGenerator);

					// Allocate each molecule position to a parent circle so defining clusters.
					int[][] clusters = new int[clusterCentres.size()][];
					int[] clusterSize = new int[clusters.length];

					for (int i = 0; i < nMolecules; i++)
					{
						double[] centre = maskDistribution.next();
						// The mask returns the coordinates with the centre of the image at 0,0
						centre[0] += width / 2;
						centre[1] += width / 2;
						xyz.add(centre);

						// Output statistics on cluster size and number.
						// TODO - Finding the closest cluster could be done better than an all-vs-all comparison
						double max = distance2(centre, clusterCentres.get(0));
						int cluster = 0;
						for (int j = 1; j < clusterCentres.size(); j++)
						{
							double d2 = distance2(centre, clusterCentres.get(j));
							if (d2 < max)
							{
								max = d2;
								cluster = j;
							}
						}

						// Assign point i to cluster
						centre[2] = cluster;

						if (clusterSize[cluster] == 0)
						{
							clusters[cluster] = new int[10];
						}
						if (clusters[cluster].length <= clusterSize[cluster])
						{
							clusters[cluster] = Arrays
									.copyOf(clusters[cluster], (int) (clusters[cluster].length * 1.5));
						}
						clusters[cluster][clusterSize[cluster]++] = i;
					}

					// Generate real cluster size statistics
					for (int j = 0; j < clusterSize.length; j++)
					{
						final int size = clusterSize[j];
						if (size == 0)
							continue;

						statsSize.add(size);

						if (size == 1)
						{
							statsRadius.add(0);
							continue;
						}

						// Find centre of cluster and add the distance to each point
						double[] com = new double[2];
						for (int n = 0; n < size; n++)
						{
							double[] xy = xyz.get(clusters[j][n]);
							for (int k = 0; k < 2; k++)
								com[k] += xy[k];
						}
						for (int k = 0; k < 2; k++)
							com[k] /= size;
						for (int n = 0; n < size; n++)
						{
							double dx = xyz.get(clusters[j][n])[0] - com[0];
							double dy = xyz.get(clusters[j][n])[1] - com[1];
							statsRadius.add(Math.sqrt(dx * dx + dy * dy));
						}
					}
				}

				if (showClusterMask)
				{
					bp = new ByteProcessor(maskSize, maskSize);
					for (int i = 0; i < mask.length; i++)
						if (mask[i] != 0)
							bp.set(i, 128);
					Utils.display(maskTitle, bp);
				}
			}

			// Use the simulated cluster centres to create clusters of the desired size
			if (clusterSimulation == 1 || clusterSimulation == 2)
			{
				for (double[] clusterCentre : clusterCentres)
				{
					int clusterN = (int) Math.round((clusterNumberSD > 0) ? dataGenerator.nextGaussian(clusterNumber,
							clusterNumberSD) : clusterNumber);
					if (clusterN < 1)
						continue;
					//double[] clusterCentre = dist.next();
					if (clusterN == 1)
					{
						// No need for a cluster around a point
						xyz.add(clusterCentre);
						statsRadius.add(0);
						statsSize.add(1);
					}
					else
					{
						// Generate N random points within a circle of the chosen cluster radius.
						// Locate the centre-of-mass and the average distance to the centre.
						double[] com = new double[3];
						int j = 0;
						while (j < clusterN)
						{
							// Generate a random point within a circle uniformly
							// http://stackoverflow.com/questions/5837572/generate-a-random-point-within-a-circle-uniformly
							double t = 2.0 * Math.PI * randomGenerator.nextDouble();
							double u = randomGenerator.nextDouble() + randomGenerator.nextDouble();
							double r = clusterRadius * ((u > 1) ? 2 - u : u);
							double x = r * Math.cos(t);
							double y = r * Math.sin(t);
							double[] xy = new double[] { clusterCentre[0] + x, clusterCentre[1] + y };
							xyz.add(xy);
							for (int k = 0; k < 2; k++)
								com[k] += xy[k];
							j++;
						}
						// Add the distance of the points from the centre of the cluster.
						// Note this does not account for the movement due to precision.
						statsSize.add(j);
						if (j == 1)
						{
							statsRadius.add(0);
						}
						else
						{
							for (int k = 0; k < 2; k++)
								com[k] /= j;
							while (j > 0)
							{
								double dx = xyz.get(xyz.size() - j)[0] - com[0];
								double dy = xyz.get(xyz.size() - j)[1] - com[1];
								statsRadius.add(Math.sqrt(dx * dx + dy * dy));
								j--;
							}
						}
					}
				}
			}
		}
		else
		{
			// Random distribution
			for (int i = 0; i < nMolecules; i++)
				xyz.add(dist.next());
		}

		// The Gaussian sigma should be applied so the overall distance from the centre
		// ( sqrt(x^2+y^2) ) has a standard deviation of sigmaS?
		final double sigma1D = sigmaS / Math.sqrt(2);

		// Show optional histograms
		StoredDataStatistics intraDistances = null;
		StoredDataStatistics blinks = null;
		if (showHistograms)
		{
			int capacity = (int) (xyz.size() * blinkingRate);
			intraDistances = new StoredDataStatistics(capacity);
			blinks = new StoredDataStatistics(capacity);
		}

		Statistics statsSigma = new Statistics();
		for (int i = 0; i < xyz.size(); i++)
		{
			int nOccurrences = getBlinks(dataGenerator, blinkingRate);
			if (showHistograms)
				blinks.add(nOccurrences);

			final int size = molecules.size();

			// Get coordinates in nm
			final double[] moleculeXyz = xyz.get(i);

			if (bp != null && nOccurrences > 0)
			{
				bp.putPixel((int) Math.round(moleculeXyz[0] / maskScale), (int) Math.round(moleculeXyz[1] / maskScale),
						255);
			}

			while (nOccurrences-- > 0)
			{
				final double[] localisationXy = Arrays.copyOf(moleculeXyz, 2);
				// Add random precision
				if (sigma1D > 0)
				{
					final double dx = dataGenerator.nextGaussian(0, sigma1D);
					final double dy = dataGenerator.nextGaussian(0, sigma1D);
					localisationXy[0] += dx;
					localisationXy[1] += dy;
					if (!dist.isWithinXY(localisationXy))
						continue;
					// Calculate mean-squared displacement
					statsSigma.add(dx * dx + dy * dy);
				}
				final double x = localisationXy[0];
				final double y = localisationXy[1];
				molecules.add(new Molecule(x, y, i, 1));

				// Store in pixels
				float[] params = new float[7];
				params[Gaussian2DFunction.X_POSITION] = (float) (x / nmPerPixel);
				params[Gaussian2DFunction.Y_POSITION] = (float) (y / nmPerPixel);
				results.add(i + 1, (int) x, (int) y, 0, 0, 0, params, null);
			}

			if (molecules.size() > size)
			{
				count++;
				if (showHistograms)
				{
					int newCount = molecules.size() - size;
					if (newCount == 1)
					{
						// No intra-molecule distances
						//intraDistances.add(0);
						continue;
					}

					// Get the distance matrix between these molecules
					double[][] matrix = new double[newCount][newCount];
					for (int ii = size, x = 0; ii < molecules.size(); ii++, x++)
					{
						for (int jj = size + 1, y = 1; jj < molecules.size(); jj++, y++)
						{
							final double d2 = molecules.get(ii).distance2(molecules.get(jj));
							matrix[x][y] = matrix[y][x] = d2;
						}
					}

					// Get the maximum distance for particle linkage clustering of this molecule
					double max = 0;
					for (int x = 0; x < newCount; x++)
					{
						// Compare to all-other molecules and get the minimum distance 
						// needed to join at least one
						double linkDistance = Double.POSITIVE_INFINITY;
						for (int y = 0; y < newCount; y++)
						{
							if (x == y)
								continue;
							if (matrix[x][y] < linkDistance)
								linkDistance = matrix[x][y];
						}
						// Check if this is larger 
						if (max < linkDistance)
							max = linkDistance;
					}
					intraDistances.add(Math.sqrt(max));
				}
			}
		}
		results.end();

		if (bp != null)
			Utils.display(maskTitle, bp);

		// Used for debugging
		//System.out.printf("  * Molecules = %d (%d activated)\n", xyz.size(), count);
		//if (clusterSimulation > 0)
		//	System.out.printf("  * Cluster number = %s +/- %s. Radius = %s +/- %s\n",
		//			Utils.rounded(statsSize.getMean(), 4), Utils.rounded(statsSize.getStandardDeviation(), 4),
		//			Utils.rounded(statsRadius.getMean(), 4), Utils.rounded(statsRadius.getStandardDeviation(), 4));

		log("Simulation results");
		log("  * Molecules = %d (%d activated)", xyz.size(), count);
		log("  * Blinking rate = %s", Utils.rounded((double) molecules.size() / xyz.size(), 4));
		log("  * Precision (Mean-displacement) = %s nm",
				(statsSigma.getN() > 0) ? Utils.rounded(Math.sqrt(statsSigma.getMean()), 4) : "0");
		if (showHistograms)
		{
			if (intraDistances.getN() == 0)
			{
				log("  * Mean Intra-Molecule particle linkage distance = 0 nm");
				log("  * Fraction of inter-molecule particle linkage @ 0 nm = 0 %%");
			}
			else
			{
				plot(blinks, "Blinks/Molecule", true);
				double[][] intraHist = plot(intraDistances, "Intra-molecule particle linkage distance", false);

				// Determine 95th and 99th percentile
				int p99 = intraHist[0].length - 1;
				double limit1 = 0.99 * intraHist[1][p99];
				double limit2 = 0.95 * intraHist[1][p99];
				while (intraHist[1][p99] > limit1 && p99 > 0)
					p99--;
				int p95 = p99;
				while (intraHist[1][p95] > limit2 && p95 > 0)
					p95--;

				log("  * Mean Intra-Molecule particle linkage distance = %s nm (95%% = %s, 99%% = %s, 100%% = %s)",
						Utils.rounded(intraDistances.getMean(), 4), Utils.rounded(intraHist[0][p95], 4),
						Utils.rounded(intraHist[0][p99], 4), Utils.rounded(intraHist[0][intraHist[0].length - 1], 4));

				if (distanceAnalysis)
				{
					performDistanceAnalysis(intraHist, p99);
				}
			}
		}
		if (clusterSimulation > 0)
		{
			log("  * Cluster number = %s +/- %s", Utils.rounded(statsSize.getMean(), 4),
					Utils.rounded(statsSize.getStandardDeviation(), 4));
			log("  * Cluster radius = %s +/- %s nm (mean distance to centre-of-mass)",
					Utils.rounded(statsRadius.getMean(), 4), Utils.rounded(statsRadius.getStandardDeviation(), 4));
		}
	}

	private double[][] plot(StoredDataStatistics stats, String label, boolean integerBins)
	{
		String title = TITLE + " " + label;
		double yMax = Maths.max(stats.getValues());
		Plot plot;
		double[][] hist = null;
		if (integerBins)
		{
			int nBins = (int) (Math.ceil(yMax) + 1);
			float[][] fHist = Utils.calcHistogram(stats.getFloatValues(), 0, yMax, nBins);

			// Create the axes
			float[] xValues = Utils.createHistogramAxis(fHist[0]);
			float[] yValues = Utils.createHistogramValues(fHist[1]);

			//// Copy for return
			//hist = new double[2][fHist[0].length];
			//for (int i=0; i<fHist[0].length; i++)
			//{
			//	hist[0][i] = fHist[0][i];
			//	hist[1][i] = fHist[1][i];
			//}

			// Plot
			yMax = Maths.max(yValues);
			plot = new Plot(title, label, "Frequency", xValues, yValues);
			if (xValues.length > 0)
			{
				double xPadding = 0.05 * (xValues[xValues.length - 1] - xValues[0]);
				plot.setLimits(xValues[0] - xPadding, xValues[xValues.length - 1] + xPadding, 0, yMax * 1.05);
			}
		}
		else
		{
			// Show a cumulative histogram so that the bin size is not relevant
			hist = Maths.cumulativeHistogram(stats.getValues(), false);

			// Create the axes
			double[] xValues = hist[0];
			double[] yValues = hist[1];

			// Plot
			plot = new Plot(title, label, "Frequency", xValues, yValues);
		}

		Utils.display(title, plot);
		return hist;
	}

	private void performDistanceAnalysis(double[][] intraHist, int p99)
	{
		// We want to know the fraction of distances between molecules at the 99th percentile
		// that are intra- rather than inter-molecule.
		// Do single linkage clustering of closest pair at this distance and count the number of 
		// links that are inter and intra.

		// Convert molecules for clustering
		ArrayList<ClusterPoint> points = new ArrayList<ClusterPoint>(molecules.size());
		for (Molecule m : molecules)
			// Precision was used to store the molecule ID
			points.add(ClusterPoint.newClusterPoint((int) m.precision, m.x, m.y, m.photons));
		ClusteringEngine engine = new ClusteringEngine(Prefs.getThreads(), ClusteringAlgorithm.PARTICLE_SINGLE_LINKAGE,
				new IJTrackProgress());
		IJ.showStatus("Clustering to check inter-molecule distances");
		engine.setTrackJoins(true);
		ArrayList<Cluster> clusters = engine.findClusters(points, intraHist[0][p99]);
		IJ.showStatus("");
		if (clusters != null)
		{
			double[] intraIdDistances = engine.getIntraIdDistances();
			double[] interIdDistances = engine.getInterIdDistances();

			int all = interIdDistances.length + intraIdDistances.length;

			log("  * Fraction of inter-molecule particle linkage @ %s nm = %s %%", Utils.rounded(intraHist[0][p99], 4),
					(all > 0) ? Utils.rounded(100.0 * interIdDistances.length / all, 4) : "0");

			// Show a double cumulative histogram plot
			double[][] intraIdHist = Maths.cumulativeHistogram(intraIdDistances, false);
			double[][] interIdHist = Maths.cumulativeHistogram(interIdDistances, false);

			// Plot
			String title = TITLE + " molecule linkage distance";
			Plot plot = new Plot(title, "Distance", "Frequency", intraIdHist[0], intraIdHist[1]);
			double max = (intraIdHist[1].length > 0) ? intraIdHist[1][intraIdHist[1].length - 1] : 0;
			if (interIdHist[1].length > 0)
				max = FastMath.max(max, interIdHist[1][interIdHist[1].length - 1]);
			plot.setLimits(0, intraIdHist[0][intraIdHist[0].length - 1], 0, max);
			plot.setColor(Color.blue);
			plot.addPoints(interIdHist[0], interIdHist[1], Plot.LINE);
			plot.setColor(Color.black);
			Utils.display(title, plot);
		}
		else
		{
			log("Aborted clustering to check inter-molecule distances");
		}
	}

	private double distance2(double[] centre1, double[] centre2)
	{
		final double dx = centre1[0] - centre2[0];
		final double dy = centre1[1] - centre2[1];
		return dx * dx + dy * dy;
	}

	/**
	 * Fill the given mask with the fill value using a circle with the specified centre and radius
	 * 
	 * @param mask
	 * @param maskSize
	 * @param cx
	 * @param cy
	 * @param roiRadius
	 * @param fill
	 */
	private void fillMask(int[] mask, final int maskSize, final int cx, final int cy, final int roiRadius,
			final int fill)
	{
		int minx = cx - roiRadius;
		int maxx = cx + roiRadius;
		int miny = cy - roiRadius;
		int maxy = cy + roiRadius;
		final int r2 = roiRadius * roiRadius;

		// Pre-calculate x range
		if (minx < 0)
			minx = 0;
		if (maxx >= maskSize)
			maxx = maskSize - 1;
		if (minx > maxx)
			return;

		int n = 0;
		int[] dx2 = new int[roiRadius * 2 + 1];
		for (int x = minx; x <= maxx; x++)
		{
			dx2[n++] = (cx - x) * (cx - x);
		}

		if (miny < 0 || maxy >= maskSize)
		{
			for (int y = miny, dy = -roiRadius; y <= maxy; y++, dy++)
			{
				if (y < 0)
					continue;
				if (y >= maskSize)
					break;
				final int limit = r2 - (dy * dy);
				for (int i = y * maskSize + minx, nn = 0; nn < n; i++, nn++)
				{
					if (dx2[nn] <= limit)
						mask[i] = fill;
				}
			}
		}
		else
		{
			for (int y = miny, dy = -roiRadius; y <= maxy; y++, dy++)
			{
				final int limit = r2 - (dy * dy);
				for (int i = y * maskSize + minx, nn = 0; nn < n; i++, nn++)
				{
					if (dx2[nn] <= limit)
						mask[i] = fill;
				}
			}
		}

	}

	private int getBlinks(RandomDataGenerator dataGenerator, double averageBlinks)
	{
		switch (blinkingDistribution)
		{
			case 3:
				// Binomial distribution
				return dataGenerator.nextBinomial((int) Math.round(averageBlinks), p);

			case 2:
				return (int) Math.round(averageBlinks);
			case 1:
				return StandardFluorophoreSequenceModel.getBlinks(true, dataGenerator, averageBlinks);
			default:
				return StandardFluorophoreSequenceModel.getBlinks(false, dataGenerator, averageBlinks);
		}
	}

	private boolean showSimulationDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Simulate a random distribution of molecules.");

		gd.addNumericField("Molecules", nMolecules, 0);
		gd.addNumericField("Simulation_size (um)", simulationSize, 2);
		gd.addNumericField("Blinking_rate", blinkingRate, 2);
		gd.addChoice("Blinking_distribution", BLINKING_DISTRIBUTION, BLINKING_DISTRIBUTION[blinkingDistribution]);
		gd.addNumericField("Average_precision (nm)", sigmaS, 2);
		gd.addCheckbox("Show_histograms", showHistograms);
		gd.addCheckbox("Distance_analysis", distanceAnalysis);

		gd.addChoice("Cluster_simulation", CLUSTER_SIMULATION, CLUSTER_SIMULATION[clusterSimulation]);
		gd.addNumericField("Cluster_number", clusterNumber, 2);
		gd.addNumericField("Cluster_variation (SD)", clusterNumberSD, 2);
		gd.addNumericField("Cluster_radius", clusterRadius, 2);
		gd.addCheckbox("Show_cluster_mask", showClusterMask);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		nMolecules = (int) Math.abs(gd.getNextNumber());
		simulationSize = Math.abs(gd.getNextNumber());
		blinkingRate = Math.abs(gd.getNextNumber());
		blinkingDistribution = gd.getNextChoiceIndex();
		sigmaS = Math.abs(gd.getNextNumber());
		showHistograms = gd.getNextBoolean();
		distanceAnalysis = gd.getNextBoolean();
		clusterSimulation = gd.getNextChoiceIndex();
		clusterNumber = Math.abs(gd.getNextNumber());
		clusterNumberSD = Math.abs(gd.getNextNumber());
		clusterRadius = Math.abs(gd.getNextNumber());
		showClusterMask = gd.getNextBoolean();

		// Check arguments
		try
		{
			Parameters.isAboveZero("Molecules", nMolecules);
			Parameters.isAboveZero("Simulation size", simulationSize);
			Parameters.isEqualOrAbove("Blinking rate", blinkingRate, 1);
			Parameters.isEqualOrAbove("Cluster number", clusterNumber, 1);
		}
		catch (IllegalArgumentException ex)
		{
			IJ.error(TITLE, ex.getMessage());
			return false;
		}

		return getPValue();
	}

	private boolean getPValue()
	{
		if (blinkingDistribution == 3)
		{
			GenericDialog gd = new GenericDialog(TITLE);
			gd.addMessage("Binomial distribution requires a p-value");
			gd.addSlider("p-Value (%)", 0, 100, 100 * p);
			gd.showDialog();
			if (gd.wasCanceled())
				return false;
			p = FastMath.max(FastMath.min(gd.getNextNumber(), 100), 0) / 100;
		}
		return true;
	}

	/**
	 * Get the lifetime of the results using the earliest and latest frames and the calibrated exposure time.
	 */
	private void getLifetime()
	{
		int start;
		int end;
		List<PeakResult> peakResults = results.getResults();
		if (peakResults.isEmpty())
		{
			seconds = 0;
			return;
		}
		start = end = peakResults.get(0).peak;
		for (PeakResult r : peakResults)
		{
			if (start > r.peak)
				start = r.peak;
			if (end < r.getEndFrame())
				end = r.getEndFrame();
		}
		seconds = (end - start + 1) * results.getCalibration().exposureTime / 1000;
	}

	private boolean createImage(ArrayList<Molecule> molecules)
	{
		if (molecules.isEmpty())
			return false;

		// Find the limits of the image
		minx = maxx = molecules.get(0).x;
		miny = maxy = molecules.get(0).y;

		// Compute limits
		for (int i = molecules.size(); i-- > 0;)
		{
			final Molecule m1 = molecules.get(i);
			if (minx > m1.x)
				minx = m1.x;
			else if (maxx < m1.x)
				maxx = m1.x;
			if (miny > m1.y)
				miny = m1.y;
			else if (maxy < m1.y)
				maxy = m1.y;
		}

		//		// Debug all-vs-all comparison
		//		long start = System.nanoTime();
		//		double dMin2 = Double.POSITIVE_INFINITY;
		//		for (int i = molecules.size(); i-- > 0;)
		//		{
		//			final Molecule m1 = molecules.get(i);
		//			for (int j = i; j-- > 0;)
		//			{
		//				final Molecule m2 = molecules.get(j);
		//				if (dMin2 > m1.distance2(m2))
		//					dMin2 = m1.distance2(m2);
		//			}
		//		}
		//		long t1 = System.nanoTime() - start;

		// Assign to a grid
		final int gridSize = 500;
		final double xBinSize = (maxx - minx) / gridSize;
		final double yBinSize = (maxy - miny) / gridSize;
		final int nXBins = 1 + (int) ((maxx - minx) / xBinSize);
		final int nYBins = 1 + (int) ((maxy - miny) / yBinSize);
		Molecule[][] grid = new Molecule[nXBins][nYBins];
		for (Molecule m : molecules)
		{
			final int xBin = (int) ((m.x - minx) / xBinSize);
			final int yBin = (int) ((m.y - miny) / yBinSize);
			// Build a single linked list
			m.next = grid[xBin][yBin];
			grid[xBin][yBin] = m;
		}

		// Find the minimum distance between molecules.
		double dMin = Double.POSITIVE_INFINITY;

		IJ.showStatus("Computing minimum distance ...");
		IJ.showProgress(0);
		Molecule[] neighbours = new Molecule[5];
		for (int yBin = 0, currentIndex = 0, finalIndex = nXBins * nYBins; yBin < nYBins; yBin++)
		{
			for (int xBin = 0; xBin < nXBins; xBin++)
			{
				IJ.showProgress(currentIndex, finalIndex);

				for (Molecule m1 = grid[xBin][yBin]; m1 != null; m1 = m1.next)
				{
					// Build a list of which cells to compare up to a maximum of 4
					//      | 0,0  |  1,0
					// ------------+-----
					// -1,1 | 0,1  |  1,1

					int count = 0;
					neighbours[count++] = m1.next;

					if (yBin < nYBins - 1)
					{
						neighbours[count++] = grid[xBin][yBin + 1];
						if (xBin > 0)
							neighbours[count++] = grid[xBin - 1][yBin + 1];
					}
					if (xBin < nXBins - 1)
					{
						neighbours[count++] = grid[xBin + 1][yBin];
						if (yBin < nYBins - 1)
							neighbours[count++] = grid[xBin + 1][yBin + 1];
					}

					// Compare to neighbours
					while (count-- > 0)
					{
						for (Molecule m2 = neighbours[count]; m2 != null; m2 = m2.next)
						{
							if (dMin > m1.distance2(m2))
								dMin = m1.distance2(m2);
						}
					}
				}
			}
		}
		IJ.showStatus("");
		IJ.showProgress(1);

		nmPerPixel = Math.sqrt(dMin);
		log("Minimum distance between molecules = %g nm", nmPerPixel);
		if (nmPerPixel == 0 && nmPerPixelLimit == 0)
		{
			IJ.error(TITLE,
					"Zero minimum distance between molecules - please enter a nm/pixel limit for image reconstruction");
			return false;
		}

		if (nmPerPixel < nmPerPixelLimit)
		{
			log("Minimum distance adjusted to user defined limit %.2f nm", nmPerPixelLimit);
			nmPerPixel = nmPerPixelLimit;
		}

		// Compute the minimum size we can use and stay within memory.
		// Assume a 4um x 4um analysis section with 800nm feature radius.
		double limit = (4000.0 + 800.0 + 1) / 4096;
		if (nmPerPixel < limit)
		{
			log("Minimum distance adjusted to %.2f nm to fit in memory", limit);
			nmPerPixel = limit;
		}
		log("X-range %.2f - %.2f : Y-range %.2f - %.2f (nm)", minx, maxx, miny, maxy);

		// Construct a binary representation

		String namePrefix = results.getName() + " " + ((binaryImage) ? "Binary" : "Count") + " Image";

		double lowResNmPerPixel;
		double xrange = maxx - minx;
		double yrange = maxy - miny;
		if (xrange > 0 || yrange > 0)
		{
			lowResNmPerPixel = FastMath.max(xrange, yrange) / lowResolutionImageSize;
		}
		else
		{
			// The resolution does not matter
			lowResNmPerPixel = 100;
		}
		ImagePlus imp = displayImage(namePrefix + " (low res)", molecules, minx, miny, maxx, maxy, lowResNmPerPixel,
				false, binaryImage);

		// Add an ROI to allow the user to select regions. PC-PALM recommends 2x2 to 4x4 um^2
		int size = (int) (roiSizeInUm * 1000.0 / lowResNmPerPixel);
		imp.setRoi(new Rectangle(0, 0, size, size));

		if (showHighResolutionImage)
			displayImage(namePrefix + " (high res)", molecules, minx, miny, maxx, maxy, nmPerPixel, false, binaryImage);

		// Store the molecules, the data range and the dMin. 
		// This will be used by a filter plugin that crops sections from the image for PC analysis
		PCPALMMolecules.molecules = molecules;

		return true;
	}

	static ImageProcessor drawImage(ArrayList<Molecule> molecules, double minx, double miny, double maxx, double maxy,
			double nmPerPixel, boolean checkBounds, boolean binary)
	{
		double scalex = maxx - minx;
		double scaley = maxy - miny;
		int width = (int) Math.round(scalex / nmPerPixel) + 1;
		int height = (int) Math.round(scaley / nmPerPixel) + 1;

		// ***
		// The PC-PALM + PLoS One papers describe using a binary image.
		// However both papers provide MatLab code where the number of particles is 
		// calculated using sum(sum(I)). This indicates a non binary image could be input
		// to the routine to calculate the correlation function g(r).
		// ***
		if (binary)
		{
			byte[] data = new byte[width * height];
			for (Molecule m : molecules)
			{
				if (checkBounds)
				{
					if (m.x < minx || m.x >= maxx || m.y < miny || m.y >= maxy)
						continue;
				}

				// Shift to the origin. This makes the image more memory efficient.
				int x = (int) Math.round((m.x - minx) / nmPerPixel);
				int y = (int) Math.round((m.y - miny) / nmPerPixel);
				int index = y * width + x;

				// Construct a binary image
				data[index] = (byte) 1;
			}

			ByteProcessor ip = new ByteProcessor(width, height, data, null);
			ip.setMinAndMax(0, 1);
			return ip;
		}
		else
		{
			short[] data = new short[width * height];
			for (Molecule m : molecules)
			{
				if (checkBounds)
				{
					if (m.x < minx || m.x >= maxx || m.y < miny || m.y >= maxy)
						continue;
				}

				// Shift to the origin. This makes the image more memory efficient.
				int x = (int) Math.round((m.x - minx) / nmPerPixel);
				int y = (int) Math.round((m.y - miny) / nmPerPixel);
				int index = y * width + x;

				// Construct a count image
				data[index]++;
			}

			ShortProcessor ip = new ShortProcessor(width, height, data, null);
			ip.setMinAndMax(0, Maths.max(data));
			return ip;
		}
	}

	static ImagePlus displayImage(String title, ImageProcessor ip, double nmPerPixel)
	{
		ImagePlus imp = Utils.display(title, ip);
		Calibration cal = new Calibration();
		cal.setUnit("um");
		cal.pixelWidth = cal.pixelHeight = nmPerPixel / 1000;
		imp.setCalibration(cal);
		return imp;
	}

	static ImagePlus displayImage(String title, ArrayList<Molecule> molecules, double minx, double miny, double maxx,
			double maxy, double nmPerPixel, boolean checkBounds, boolean binary)
	{
		ImageProcessor ip = drawImage(molecules, minx, miny, maxx, maxy, nmPerPixel, checkBounds, binary);
		return displayImage(title, ip, nmPerPixel);
	}

	/**
	 * Allow optimisation using Apache Commons Math 3 Optimiser
	 */
	public abstract class SkewNormalOptimiserFunction extends SkewNormalFunction
	{
		public SkewNormalOptimiserFunction(double[] parameters)
		{
			super(parameters);
		}

		protected List<Double> x = null;
		protected List<Double> y = null;

		public void addPoint(double x, double y)
		{
			if (this.x == null)
			{
				this.x = new ArrayList<Double>();
				this.y = new ArrayList<Double>();
			}
			this.x.add(x);
			this.y.add(y);
		}

		public void addData(float[] x, float[] y)
		{
			this.x = new ArrayList<Double>();
			this.y = new ArrayList<Double>();
			for (int i = 0; i < x.length; i++)
			{
				this.x.add((double) x[i]);
				this.y.add((double) y[i]);
			}
		}

		public double[] calculateTarget()
		{
			double[] target = new double[y.size()];
			for (int i = 0; i < y.size(); i++)
			{
				target[i] = y.get(i).doubleValue();
			}
			return target;
		}

		public double[] calculateWeights()
		{
			double[] w = new double[y.size()];
			for (int i = 0; i < y.size(); i++)
			{
				w[i] = 1;
			}
			return w;
		}
	}

	/**
	 * Allow optimisation using Apache Commons Math 3 Gradient Optimiser
	 */
	public class SkewNormalDifferentiableFunction extends SkewNormalOptimiserFunction implements
			MultivariateVectorFunction
	{
		// Adapted from http://commons.apache.org/proper/commons-math/userguide/optimization.html
		// Use the deprecated API since the new one is not yet documented.

		public SkewNormalDifferentiableFunction(double[] parameters)
		{
			super(parameters);
		}

		private double[][] jacobian(double[] variables)
		{
			// Compute the gradients using numerical differentiation

			double[][] jacobian = new double[x.size()][4];
			final double delta = 0.001;
			double[][] d = new double[variables.length][variables.length];
			for (int i = 0; i < variables.length; i++)
				d[i][i] = delta * Math.abs(variables[i]); // Should the delta be changed for each parameter ?
			for (int i = 0; i < jacobian.length; ++i)
			{
				double x = this.x.get(i);
				double value = evaluate(x, variables);
				for (int j = 0; j < variables.length; j++)
				{
					double value2 = evaluate(x, variables[0] + d[0][j], variables[1] + d[1][j], variables[2] + d[2][j],
							variables[3] + d[3][j]);
					jacobian[i][j] = (value2 - value) / d[j][j];
				}
			}
			return jacobian;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see org.apache.commons.math3.analysis.MultivariateVectorFunction#value(double[])
		 */
		public double[] value(double[] variables)
		{
			double[] values = new double[x.size()];
			for (int i = 0; i < values.length; i++)
				values[i] = evaluate(x.get(i), variables);
			return values;
		}
	}

	/**
	 * Allow optimisation using Apache Commons Math 3 Simplex
	 */
	public class SkewNormalMultivariateFunction extends SkewNormalOptimiserFunction implements MultivariateFunction
	{
		public SkewNormalMultivariateFunction(double[] parameters)
		{
			super(parameters);
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see org.apache.commons.math3.analysis.MultivariateFunction#value(double[])
		 */
		public double value(double[] point)
		{
			// Objective function is to minimise sum-of-squares
			double ss = 0;
			for (int i = x.size(); i-- > 0;)
			{
				double dx = y.get(i) - evaluate(x.get(i), point);
				ss += dx * dx;
			}
			return ss;
		}
	}

	/**
	 * Apache GaussianFitter extended to allow maxEvaluations to be set
	 * 
	 */
	public class MyGaussianFitter extends GaussianFitter
	{
		public int maxEvaluations = 1000;
		public double[] guess = null;

		public MyGaussianFitter(MultivariateVectorOptimizer optimizer, int maxEvaluations)
		{
			super(optimizer);
			this.maxEvaluations = maxEvaluations;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see org.apache.commons.math3.fitting.GaussianFitter#fit(double[])
		 */
		public double[] fit(double[] initialGuess)
		{
			final Gaussian.Parametric f = new Gaussian.Parametric()
			{
				@Override
				public double value(double x, double... p)
				{
					double v = Double.POSITIVE_INFINITY;
					try
					{
						v = super.value(x, p);
					}
					catch (NotStrictlyPositiveException e)
					{ // NOPMD
					  // Do nothing.
					}
					return v;
				}

				@Override
				public double[] gradient(double x, double... p)
				{
					double[] v = { Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY };
					try
					{
						v = super.gradient(x, p);
					}
					catch (NotStrictlyPositiveException e)
					{ // NOPMD
					  // Do nothing.
					}
					return v;
				}
			};

			return fit(maxEvaluations, f, initialGuess);
		}

		/**
		 * Fits a Gaussian function to the observed points.
		 * 
		 * @return the parameters of the Gaussian function that best fits the
		 *         observed points (in the same order as above).
		 */
		public double[] fit()
		{
			guess = (new ParameterGuesser(getObservations())).guess();
			return fit(guess);
		}
	}

	/**
	 * Log a message to the IJ log window
	 * 
	 * @param format
	 * @param args
	 */
	private static void log(String format, Object... args)
	{
		Utils.log(format, args);
	}

	/**
	 * Output a spacer to the ImageJ log
	 */
	static void logSpacer()
	{
		log("-=-=-=-=-=-=-");
	}
}
