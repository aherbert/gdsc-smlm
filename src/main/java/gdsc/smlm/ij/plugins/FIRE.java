package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.Checkbox;
import java.awt.Color;
import java.awt.Rectangle;
import java.awt.TextField;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.MathArrays;

import gdsc.core.ij.IJTrackProgress;
import gdsc.core.ij.Utils;
import gdsc.core.logging.NullTrackProgress;
import gdsc.core.logging.TrackProgress;
import gdsc.core.utils.Maths;
import gdsc.core.utils.Statistics;
import gdsc.core.utils.StoredDataStatistics;
import gdsc.core.utils.TurboList;

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

import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.ij.frc.FRC;
import gdsc.smlm.ij.frc.FRC.ThresholdMethod;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.results.IJImagePeakResults;
import gdsc.smlm.ij.results.ImagePeakResultsFactory;
import gdsc.smlm.ij.results.ResultsImage;
import gdsc.smlm.ij.results.ResultsMode;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gnu.trove.list.array.TDoubleArrayList;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Plot;
import ij.gui.Plot2;
import ij.gui.PlotWindow;
import ij.plugin.PlugIn;
import ij.plugin.WindowOrganiser;
import ij.plugin.frame.Recorder;
import ij.process.ImageProcessor;
import ij.process.LUT;
import ij.process.LUTHelper;
import ij.process.LUTHelper.LutColour;

/**
 * Computes the Fourier Image Resolution of an image
 * <p>
 * Implements the FIRE (Fourier Image REsolution) method described in:<br>
 * Niewenhuizen, et al (2013). Measuring image resolution in optical nanoscopy. Nature Methods, 10, 557<br>
 * http://www.nature.com/nmeth/journal/v10/n6/full/nmeth.2448.html
 */
public class FIRE implements PlugIn
{
	private static String TITLE = "Fourier Image REsolution (FIRE)";
	private static String inputOption = "";
	private static String inputOption2 = "";

	private static int repeats = 1;
	private static boolean useSignal = false;
	private static int maxPerBin = 0; // 5 in the Niewenhuizen paper
	private static boolean randomSplit = true;
	private static int blockSize = 50;
	private static String[] SCALE_ITEMS;
	private static int[] SCALE_VALUES = new int[] { 0, 1, 2, 4, 8, 16, 32, 64, 128 };
	private static String[] IMAGE_SIZE_ITEMS;
	private static int[] IMAGE_SIZE_VALUES;
	private static int imageScaleIndex = 0;
	private static int imageSizeIndex;

	// The Q value and the mean and sigma for spurious correlation correction
	private static boolean spuriousCorrelationCorrection = false;
	private static double qValue, mean, sigma;

	static
	{
		SCALE_ITEMS = new String[SCALE_VALUES.length];
		SCALE_ITEMS[0] = "Auto";
		for (int i = 1; i < SCALE_VALUES.length; i++)
			SCALE_ITEMS[i] = Integer.toString(SCALE_VALUES[i]);

		// Create size for Fourier transforms. Must be power of 2.
		IMAGE_SIZE_VALUES = new int[32];
		IMAGE_SIZE_ITEMS = new String[IMAGE_SIZE_VALUES.length];
		int size = 512; // Start at a reasonable size. Too small does not work.
		int count = 0;
		while (size <= 16384)
		{
			if (size == 2048)
				imageSizeIndex = count;

			// Image sizes are 1 smaller so that rounding error when scaling does not create an image too large for the power of 2
			IMAGE_SIZE_VALUES[count] = size - 1;
			IMAGE_SIZE_ITEMS[count] = Integer.toString(size);
			size *= 2;
			count++;
		}
		IMAGE_SIZE_VALUES = Arrays.copyOf(IMAGE_SIZE_VALUES, count);
		IMAGE_SIZE_ITEMS = Arrays.copyOf(IMAGE_SIZE_ITEMS, count);
	}

	private static double perimeterSamplingFactor = 1;
	private static boolean useHalfCircle = true;
	private static int thresholdMethodIndex = 0;
	private static boolean showFRCCurve = true;
	private static boolean showFRCCurveRepeats = false;
	private static boolean showFRCTimeEvolution = false;

	private static boolean chooseRoi = false;
	private static String roiImage = "";

	private boolean extraOptions;
	private Rectangle roiBounds;
	private int roiImageWidth, roiImageHeight;

	// Stored in initialisation
	MemoryPeakResults results, results2;
	Rectangle2D dataBounds;
	String units;
	double nmPerPixel = 1;

	// Stored in setCorrectionParameters
	private double correctionQValue, correctionMean, correctionSigma;

	public class FireImages
	{
		final ImageProcessor ip1, ip2;
		final double nmPerPixel;

		FireImages(ImageProcessor ip1, ImageProcessor ip2, double nmPerPixel)
		{
			this.ip1 = ip1;
			this.ip2 = ip2;
			this.nmPerPixel = nmPerPixel;
		}
	}

	public class FireResult
	{
		final double fireNumber;
		final double correlation;
		final double nmPerPixel;
		final double[][] frcCurve;
		final double[][] smoothedFrcCurve;

		FireResult(double fireNumber, double correlation, double nmPerPixel, double[][] frcCurve,
				double[][] smoothedFrcCurve)
		{
			this.fireNumber = fireNumber;
			this.correlation = correlation;
			this.nmPerPixel = nmPerPixel;
			this.frcCurve = frcCurve;
			this.smoothedFrcCurve = smoothedFrcCurve;
		}
	}

	public class FIREWorker implements Runnable
	{
		final int id;
		final ThresholdMethod method;
		final double fourierImageScale;
		final int imageSize;

		String name;
		FireResult result;
		Plot2 plot;

		public FIREWorker(int id, ThresholdMethod method, double fourierImageScale, int imageSize)
		{
			this.id = id;
			this.method = method;
			this.fourierImageScale = fourierImageScale;
			this.imageSize = imageSize;
			name = results.getName() + " [" + id + "]";
		}

		public void run()
		{
			result = calculateFireNumber(method, fourierImageScale, imageSize);
			if (showFRCCurve)
			{
				plot = createFrcCurve(name, result, method);
				if (showFRCCurveRepeats)
					// Do this on the thread
					plot.draw();
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		extraOptions = Utils.isExtraOptions();
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		// Require some fit results and selected regions
		int size = MemoryPeakResults.countMemorySize();
		if (size == 0)
		{
			IJ.error(TITLE, "There are no fitting results in memory");
			return;
		}

		if ("q".equals(arg))
		{
			runQEstimation();
			return;
		}

		if (!showDialog())
			return;

		MemoryPeakResults results = ResultsManager.loadInputResults(inputOption, false);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			return;
		}
		MemoryPeakResults results2 = ResultsManager.loadInputResults(inputOption2, false);

		results = cropToRoi(results);
		if (results.size() < 2)
		{
			IJ.error(TITLE, "No results within the crop region");
			return;
		}
		if (results2 != null)
		{
			results2 = cropToRoi(results2);
			if (results2.size() < 2)
			{
				IJ.error(TITLE, "No results2 within the crop region");
				return;
			}
		}

		long start = System.currentTimeMillis();

		ThresholdMethod method = FRC.ThresholdMethod.values()[thresholdMethodIndex];

		// Compute FIRE
		initialise(results, results2);

		String name = results.getName();
		double fourierImageScale = SCALE_VALUES[imageScaleIndex];
		int imageSize = IMAGE_SIZE_VALUES[imageSizeIndex];

		if (this.results2 != null)
		{
			name += " vs " + results2.getName();

			FireResult result = calculateFireNumber(method, fourierImageScale, imageSize);

			if (result != null)
			{
				IJ.log(String.format("%s : FIRE number = %s %s (Fourier scale = %s)", name,
						Utils.rounded(result.fireNumber, 4), units, Utils.rounded(nmPerPixel / result.nmPerPixel, 3)));

				if (showFRCCurve)
					showFrcCurve(name, result, method);
			}
		}
		else
		{
			FireResult result = null;

			int repeats = (randomSplit) ? Math.max(1, FIRE.repeats) : 1;
			if (repeats == 1)
			{
				result = calculateFireNumber(method, fourierImageScale, imageSize);

				if (result != null)
				{
					IJ.log(String.format("%s : FIRE number = %s %s (Fourier scale = %s)", name,
							Utils.rounded(result.fireNumber, 4), units,
							Utils.rounded(nmPerPixel / result.nmPerPixel, 3)));

					if (showFRCCurve)
						showFrcCurve(name, result, method);
				}
			}
			else
			{
				// Multi-thread this ... 			
				int nThreads = Maths.min(repeats, Prefs.getThreads());
				ExecutorService executor = Executors.newFixedThreadPool(nThreads);
				TurboList<Future<?>> futures = new TurboList<Future<?>>(repeats);
				TurboList<FIREWorker> workers = new TurboList<FIREWorker>(repeats);
				setProgress(repeats);
				IJ.showProgress(0);
				IJ.showStatus(TITLE + " computing ...");
				for (int i = 1; i <= repeats; i++)
				{
					FIREWorker w = new FIREWorker(i, method, fourierImageScale, imageSize);
					workers.add(w);
					futures.add(executor.submit(w));
				}

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
						// Ignore it and allow processing to continue (the number of neighbour samples will just be smaller).  
						e.printStackTrace();
					}
				}
				IJ.showProgress(1);

				executor.shutdown();

				// Show a combined FRC curve plot of all the smoothed curves if we have multiples.
				LUT valuesLUT = LUTHelper.createLUT(LutColour.FIRE_GLOW);
				@SuppressWarnings("unused")
				LUT noSmoothLUT = LUTHelper.createLUT(LutColour.GRAYS).createInvertedLut(); // Black at max value
				LUTHelper.DefaultLUTMapper mapper = new LUTHelper.DefaultLUTMapper(0, repeats);
				FrcCurve curve = new FrcCurve();

				Statistics stats = new Statistics();
				WindowOrganiser wo = new WindowOrganiser();
				for (int i = 0; i < repeats; i++)
				{
					FIREWorker w = workers.get(i);
					if (w.result == null)
						continue;
					result = w.result;
					if (!Double.isNaN(result.fireNumber))
						stats.add(result.fireNumber);

					if (showFRCCurveRepeats)
					{
						// Output each FRC curve using a suffix.
						IJ.log(String.format("%s : FIRE number = %s %s (Fourier scale = %s)", w.name,
								Utils.rounded(result.fireNumber, 4), units,
								Utils.rounded(nmPerPixel / result.nmPerPixel, 3)));
						wo.add(Utils.display(w.plot.getTitle(), w.plot));
					}
					if (showFRCCurve)
					{
						int index = mapper.map(i + 1);
						//@formatter:off
						curve.add(name, result, method, 
								LUTHelper.getColour(valuesLUT, index),
								Color.blue, 
								null //LUTHelper.getColour(noSmoothLUT, index)
								);
						//@formatter:on
					}
				}

				if (result != null)
				{
					wo.cascade();
					double mean = stats.getMean();
					IJ.log(String.format("%s : FIRE number = %s +/- %s %s [95%% C.I.] (Fourier scale = %s)", name,
							Utils.rounded(mean, 4), Utils.rounded(stats.getConfidenceInterval(0.95), 4), units,
							Utils.rounded(nmPerPixel / result.nmPerPixel, 3)));
					if (showFRCCurve)
					{
						curve.addResolution(mean);
						Plot2 plot = curve.getPlot();
						Utils.display(plot.getTitle(), plot);
					}
				}
			}

			// Only do this once
			if (showFRCTimeEvolution && result != null && !Double.isNaN(result.fireNumber))
				showFrcTimeEvolution(name, result.fireNumber, method, nmPerPixel / result.nmPerPixel, imageSize);
		}

		IJ.showStatus(TITLE + " complete : " + Utils.timeToString(System.currentTimeMillis() - start));
	}

	private MemoryPeakResults cropToRoi(MemoryPeakResults results)
	{
		if (roiBounds == null)
			return results;

		// Adjust bounds relative to input results image
		Rectangle2D.Float bounds = results.getDataBounds();
		double xscale = roiImageWidth / bounds.width;
		double yscale = roiImageHeight / bounds.height;

		float minX = (float) (bounds.x + roiBounds.x / xscale);
		float maxX = (float) (minX + roiBounds.width / xscale);
		float minY = (float) (bounds.y + (roiBounds.y / yscale));
		float maxY = (float) (minY + roiBounds.height / yscale);

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

	private boolean showDialog()
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

		gd.addMessage("Compute the resolution using Fourier Ring Correlation");

		ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);
		ResultsManager.addInput(gd, "Input2", inputOption2, InputSource.NONE, InputSource.MEMORY);

		gd.addCheckbox("Use_signal (if present)", useSignal);
		gd.addNumericField("Max_per_bin", maxPerBin, 0);

		gd.addMessage("For single datsets:");
		gd.addNumericField("Block_size", blockSize, 0);
		gd.addCheckbox("Random_split", randomSplit);
		gd.addNumericField("Repeats", repeats, 0);
		gd.addCheckbox("Show_FRC_curve_repeats", showFRCCurveRepeats);
		gd.addCheckbox("Show_FRC_time_evolution", showFRCTimeEvolution);
		if (extraOptions)
		{
			gd.addCheckbox("Spurious correlation correction", spuriousCorrelationCorrection);
			gd.addNumericField("Q-value", qValue, 3);
			gd.addNumericField("Precision_Mean", mean, 2);
			gd.addNumericField("Precision_Sigma", sigma, 2);
		}

		gd.addMessage("Fourier options:");
		gd.addChoice("Fourier_image_scale", SCALE_ITEMS, SCALE_ITEMS[imageScaleIndex]);
		gd.addChoice("Auto_image_scale", IMAGE_SIZE_ITEMS, IMAGE_SIZE_ITEMS[imageSizeIndex]);
		gd.addSlider("Sampling_factor", 0.2, 4, perimeterSamplingFactor);
		gd.addCheckbox("Half_circle", useHalfCircle);
		String[] methodNames = SettingsManager.getNames((Object[]) FRC.ThresholdMethod.values());
		gd.addChoice("Threshold_method", methodNames, methodNames[thresholdMethodIndex]);
		gd.addCheckbox("Show_FRC_curve", showFRCCurve);
		if (!titles.isEmpty())
			gd.addCheckbox((titles.size() == 1) ? "Use_ROI" : "Choose_ROI", chooseRoi);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		inputOption = ResultsManager.getInputSource(gd);
		inputOption2 = ResultsManager.getInputSource(gd);
		useSignal = gd.getNextBoolean();
		maxPerBin = Math.abs((int) gd.getNextNumber());

		blockSize = Math.max(1, (int) gd.getNextNumber());
		randomSplit = gd.getNextBoolean();
		repeats = Math.max(1, (int) gd.getNextNumber());
		showFRCCurveRepeats = gd.getNextBoolean();
		showFRCTimeEvolution = gd.getNextBoolean();
		if (extraOptions)
		{
			spuriousCorrelationCorrection = gd.getNextBoolean();
			qValue = Math.abs(gd.getNextNumber());
			mean = Math.abs(gd.getNextNumber());
			sigma = Math.abs(gd.getNextNumber());
		}
		
		imageScaleIndex = gd.getNextChoiceIndex();
		imageSizeIndex = gd.getNextChoiceIndex();
		perimeterSamplingFactor = gd.getNextNumber();
		useHalfCircle = gd.getNextBoolean();
		thresholdMethodIndex = gd.getNextChoiceIndex();
		showFRCCurve = gd.getNextBoolean();

		// Check arguments
		try
		{
			Parameters.isAboveZero("Perimeter sampling factor", perimeterSamplingFactor);
			if (extraOptions && spuriousCorrelationCorrection)
			{
				Parameters.isAboveZero("Q-value", qValue);
				Parameters.isAboveZero("Precision Mean", mean);
				Parameters.isAboveZero("Precision Sigma", sigma);
				// Set these for use in FIRE computation 
				setCorrectionParameters(qValue, mean, sigma);
			}
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		if (!titles.isEmpty())
			chooseRoi = gd.getNextBoolean();

		if (!titles.isEmpty() && chooseRoi)
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

		return true;
	}

	/**
	 * Initialise this instance with results.
	 *
	 * @param results
	 *            the results
	 * @param results2
	 *            the second set of results (can be null)
	 */
	public void initialise(MemoryPeakResults results, MemoryPeakResults results2)
	{
		this.results = verify(results);
		this.results2 = verify(results2);
		nmPerPixel = 1;
		units = "px";

		if (this.results == null)
			return;

		// Use the float data bounds. This prevents problems if the data is far from the origin.
		dataBounds = results.getDataBounds();

		if (this.results2 != null)
		{
			Rectangle2D dataBounds2 = results.getDataBounds();
			dataBounds = dataBounds.createUnion(dataBounds2);
		}

		if (results.getCalibration() != null)
		{
			// Calibration must match between datasets
			if (this.results2 != null)
			{
				if (results2.getNmPerPixel() != results.getNmPerPixel())
				{
					IJ.log(TITLE +
							" Error: Calibration between the two input datasets does not match, defaulting to pixels");
					return;
				}
			}

			nmPerPixel = results.getNmPerPixel();
			units = "nm";
		}
	}

	/**
	 * Sets the correction parameters for spurious correlation correction. Only relevant for single images.
	 *
	 * @param qValue
	 *            the q value
	 * @param mean
	 *            the mean of the localisation precision
	 * @param sigma
	 *            the standard deviation of the localisation precision
	 */
	public void setCorrectionParameters(double qValue, double mean, double sigma)
	{
		if (qValue > 0 && mean > 0 && sigma > 0)
		{
			correctionQValue = qValue;
			correctionMean = mean;
			correctionSigma = sigma;
		}
		else
		{
			correctionQValue = correctionMean = correctionSigma = 0;
		}
	}

	/**
	 * Copy this instance so skipping initialisation.
	 *
	 * @return the new FIRE instance
	 */
	private FIRE copy()
	{
		FIRE f = new FIRE();
		f.results = results;
		f.results2 = results2;
		f.nmPerPixel = nmPerPixel;
		f.units = units;
		f.dataBounds = dataBounds;
		f.correctionQValue = correctionQValue;
		f.correctionMean = correctionMean;
		f.correctionSigma = correctionSigma;
		return f;
	}

	/**
	 * Verify.
	 *
	 * @param results
	 *            the results
	 * @return the memory peak results
	 */
	private MemoryPeakResults verify(MemoryPeakResults results)
	{
		if (results == null || results.size() < 2)
			return null;
		if (blockSize > 1)
			// Results must be in time order when processing blocks
			results.sort();
		return results;
	}

	public FireImages createImages(double fourierImageScale, int imageSize)
	{
		return createImages(fourierImageScale, imageSize, useSignal);
	}

	public FireImages createImages(double fourierImageScale, int imageSize, boolean useSignal)
	{
		if (results == null)
			return null;

		final boolean hasSignal = useSignal && (results.getHead().getSignal() > 0);

		// Draw images using the existing IJ routines.
		Rectangle bounds = new Rectangle(0, 0, (int) Math.ceil(dataBounds.getWidth()),
				(int) Math.ceil(dataBounds.getHeight()));

		boolean weighted = true;
		boolean equalised = false;
		double imageScale;
		if (fourierImageScale == 0)
		{
			double size = FastMath.max(bounds.width, bounds.height);
			if (size <= 0)
				size = 1;
			imageScale = imageSize / size;
		}
		else
			imageScale = fourierImageScale;

		IJImagePeakResults image1 = ImagePeakResultsFactory.createPeakResultsImage(ResultsImage.NONE, weighted,
				equalised, "IP1", bounds, 1, 1, imageScale, 0, ResultsMode.ADD);
		image1.setDisplayImage(false);
		image1.begin();

		IJImagePeakResults image2 = ImagePeakResultsFactory.createPeakResultsImage(ResultsImage.NONE, weighted,
				equalised, "IP2", bounds, 1, 1, imageScale, 0, ResultsMode.ADD);
		image2.setDisplayImage(false);
		image2.begin();

		float minx = (float) dataBounds.getX();
		float miny = (float) dataBounds.getY();

		if (this.results2 != null)
		{
			// Two image comparison
			for (PeakResult p : results)
			{
				float x = p.getXPosition() - minx;
				float y = p.getYPosition() - miny;
				float v = (hasSignal) ? p.getSignal() : 1f;
				image1.add(x, y, v);
			}
			for (PeakResult p : results2)
			{
				float x = p.getXPosition() - minx;
				float y = p.getYPosition() - miny;
				float v = (hasSignal) ? p.getSignal() : 1f;
				image2.add(x, y, v);
			}
		}
		else
		{
			// Block sampling.
			// Ensure we have at least 2 even sized blocks.
			int blockSize = Math.min(results.size() / 2, Math.max(1, FIRE.blockSize));
			int nBlocks = (int) Math.ceil((double) results.size() / blockSize);
			while (nBlocks <= 1 && blockSize > 1)
			{
				blockSize /= 2;
				nBlocks = (int) Math.ceil((double) results.size() / blockSize);
			}
			if (nBlocks <= 1)
				// This should not happen since the results should contain at least 2 localisations
				return null;
			if (blockSize != FIRE.blockSize)
				IJ.log(TITLE + " Warning: Changed block size to " + blockSize);

			int i = 0;
			int block = 0;
			PeakResult[][] blocks = new PeakResult[nBlocks][blockSize];
			for (PeakResult p : results)
			{
				if (i == blockSize)
				{
					block++;
					i = 0;
				}
				blocks[block][i++] = p;
			}
			// Truncate last block
			blocks[block] = Arrays.copyOf(blocks[block], i);

			final int[] indices = Utils.newArray(nBlocks, 0, 1);
			if (randomSplit)
				MathArrays.shuffle(indices);

			for (int index : indices)
			{
				// Split alternating so just rotate
				IJImagePeakResults image = image1;
				image1 = image2;
				image2 = image;
				for (PeakResult p : blocks[index])
				{
					float x = p.getXPosition() - minx;
					float y = p.getYPosition() - miny;
					float v = (hasSignal) ? p.getSignal() : 1f;
					image.add(x, y, v);
				}
			}
		}

		image1.end();
		ImageProcessor ip1 = image1.getImagePlus().getProcessor();

		image2.end();
		ImageProcessor ip2 = image2.getImagePlus().getProcessor();

		if (!hasSignal && maxPerBin > 0)
		{
			// We can eliminate over-sampled pixels
			for (int i = ip1.getPixelCount(); i-- > 0;)
			{
				if (ip1.getf(i) > maxPerBin)
					ip1.setf(i, maxPerBin);
				if (ip2.getf(i) > maxPerBin)
					ip2.setf(i, maxPerBin);
			}
		}

		return new FireImages(ip1, ip2, nmPerPixel / imageScale);
	}

	/**
	 * Encapsulate plotting the FRC curve to allow multiple curves to be plotted together
	 */
	private class FrcCurve
	{
		double[] xValues = null;
		double[] threshold = null;
		Plot2 plot = null;

		void add(String name, FireResult result, ThresholdMethod method, Color colorValues, Color colorThreshold,
				Color colorNoSmooth)
		{
			double[][] frcCurve = result.smoothedFrcCurve;
			double[][] frcNoSmooth = result.frcCurve;

			double[] yValues = new double[frcCurve.length];
			double[] yValuesNotSmooth = new double[frcCurve.length];

			if (plot == null)
			{
				String title = name + " FRC Curve";
				plot = new Plot2(title, String.format("Spatial Frequency (%s^-1)", units), "FRC");

				xValues = new double[frcCurve.length];
				final double L = FRC.computeL(frcCurve);
				final double conversion = 1.0 / (L * result.nmPerPixel);
				for (int i = 0; i < xValues.length; i++)
				{
					final double radius = frcCurve[i][0];
					xValues[i] = radius * conversion;
				}

				// The threshold curve is the same
				threshold = FRC.calculateThresholdCurve(frcCurve, method);
				add(colorThreshold, threshold);
			}

			for (int i = 0; i < xValues.length; i++)
			{
				yValues[i] = frcCurve[i][1];
				if (frcNoSmooth != null)
					yValuesNotSmooth[i] = frcNoSmooth[i][1];
			}

			add(colorValues, yValues);
			if (frcNoSmooth != null)
				add(colorNoSmooth, yValuesNotSmooth);
		}

		public void addResolution(double resolution)
		{
			// Convert back to nm^-1
			double x = 1 / resolution;

			// Find the intersection with the threshold line
			for (int i = 1; i < xValues.length; i++)
			{
				if (x < xValues[i])
				{
					double correlation;
					// Interpolate
					double upper = xValues[i], lower = xValues[i - 1];
					double xx = (x - lower) / (upper - lower);
					correlation = threshold[i - 1] + xx * (threshold[i] - threshold[i - 1]);
					addResolution(resolution, correlation);
					return;
				}
			}
		}

		public void addResolution(double resolution, double correlation)
		{
			// Convert back to nm^-1
			double x = 1 / resolution;
			plot.setColor(Color.MAGENTA);
			plot.drawLine(x, 0, x, correlation);
		}

		private void add(Color color, double[] y)
		{
			if (color == null)
				return;
			plot.setColor(color);
			plot.addPoints(xValues, y, Plot2.LINE);
		}

		Plot2 getPlot()
		{
			plot.setLimitsToFit(false);
			// Q. For some reason the limits calculated are ignored,
			// so set them as the defaults.
			double[] limits = plot.getCurrentMinAndMax();
			// The FRC should not go above 1 so limit Y
			plot.setLimits(limits[0], limits[1], limits[2], Math.min(1.05, limits[3]));
			return plot;
		}
	}

	private Plot2 createFrcCurve(String name, FireResult result, ThresholdMethod method)
	{
		FrcCurve curve = new FrcCurve();
		curve.add(name, result, method, Color.red, Color.blue, Color.black);
		curve.addResolution(result.fireNumber, result.correlation);
		return curve.getPlot();
	}

	private void showFrcCurve(String name, FireResult result, ThresholdMethod method)
	{
		Plot2 plot = createFrcCurve(name, result, method);
		Utils.display(plot.getTitle(), plot);
	}

	private void showFrcTimeEvolution(String name, double fireNumber, ThresholdMethod method, double fourierImageScale,
			int imageSize)
	{
		IJ.showStatus("Calculating FRC time evolution curve...");

		List<PeakResult> list = results.getResults();

		int nSteps = 10;
		int maxT = list.get(list.size() - 1).peak;
		if (maxT == 0)
			maxT = list.size();
		int step = maxT / nSteps;

		TDoubleArrayList x = new TDoubleArrayList();
		TDoubleArrayList y = new TDoubleArrayList();

		double yMin = fireNumber;
		double yMax = fireNumber;

		MemoryPeakResults newResults = new MemoryPeakResults();
		newResults.copySettings(results);
		int i = 0;

		for (int t = step; t <= maxT - step; t += step)
		{
			while (i < list.size())
			{
				if (list.get(i).peak <= t)
				{
					newResults.add(list.get(i));
					i++;
				}
				else
					break;
			}

			x.add((double) t);

			FIRE f = this.copy();
			FireResult result = f.calculateFireNumber(method, fourierImageScale, imageSize);
			double fire = (result == null) ? 0 : result.fireNumber;
			y.add(fire);

			yMin = FastMath.min(yMin, fire);
			yMax = FastMath.max(yMax, fire);
		}

		// Add the final fire number
		x.add((double) maxT);
		y.add(fireNumber);

		double[] xValues = x.toArray();
		double[] yValues = y.toArray();

		String units = "px";
		if (results.getCalibration() != null)
		{
			nmPerPixel = results.getNmPerPixel();
			units = "nm";
		}

		String title = name + " FRC Time Evolution";
		Plot2 plot = new Plot2(title, "Frames", "Resolution (" + units + ")", (float[]) null, (float[]) null);
		double range = Math.max(1, yMax - yMin) * 0.05;
		plot.setLimits(xValues[0], xValues[xValues.length - 1], yMin - range, yMax + range);
		plot.setColor(Color.red);
		plot.addPoints(xValues, yValues, Plot.CONNECTED_CIRCLES);

		Utils.display(title, plot);
	}

	/**
	 * Calculate the Fourier Image REsolution (FIRE) number using the chosen threshold method. Should be called after
	 * {@link #initialise(MemoryPeakResults)}
	 * 
	 * @param method
	 * @param method
	 * @param fourierImageScale
	 *            The scale to use when reconstructing the super-resolution images (0 for auto)
	 * @param imageSize
	 *            The width of the super resolution images when using auto scale (should be a power of two minus 1 for
	 *            optimum memory usage)
	 * @return The FIRE number
	 */
	public FireResult calculateFireNumber(ThresholdMethod method, double fourierImageScale, int imageSize)
	{
		FireImages images = createImages(fourierImageScale, imageSize);
		return calculateFireNumber(method, images);
	}

	private TrackProgress progress = new IJTrackProgress();

	private void setProgress(int repeats)
	{
		if (repeats > 1)
			progress = new ParallelTrackProgress(repeats);
		else
			progress = new IJTrackProgress();
	}

	/**
	 * Dumb implementation of the track progress interface for parallel threads. Used simple synchronisation to
	 * increment total progress.
	 */
	private static class ParallelTrackProgress extends NullTrackProgress
	{
		double done = 0;
		final int total;

		ParallelTrackProgress(int repeats)
		{
			total = repeats;
		}

		@Override
		public void incrementProgress(double fraction)
		{
			// Avoid synchronisation for nothing
			if (fraction == 0)
				return;
			double done = add(fraction);
			IJ.showProgress(done / this.total);
		}

		synchronized double add(double d)
		{
			done += d;
			return done;
		}
	}

	/**
	 * Calculate the Fourier Image REsolution (FIRE) number using the chosen threshold method. Should be called after
	 * {@link #initialise(MemoryPeakResults)}
	 *
	 * @param method
	 *            the method
	 * @param ip1
	 *            the ip 1
	 * @param ip2
	 *            the ip 2
	 * @return The FIRE number
	 */
	public FireResult calculateFireNumber(ThresholdMethod method, FireImages images)
	{
		if (images == null)
			return null;

		FRC frc = new FRC();
		// Allow a progress tracker to be input.
		// This should be setup for the total number of repeats. 
		// If parallelised then do not output the text status messages as they conflict. 
		frc.progress = progress;
		frc.perimeterSamplingFactor = perimeterSamplingFactor;
		frc.useHalfCircle = useHalfCircle;
		double[][] frcCurve = frc.calculateFrcCurve(images.ip1, images.ip2);
		if (correctionQValue > 0)
			FRC.applyQCorrection(frcCurve, images.nmPerPixel, correctionQValue, correctionMean, correctionSigma);
		double[][] smoothedFrcCurve = frc.getSmoothedCurve(frcCurve);

		// Resolution in pixels
		double[] result = frc.calculateFire(smoothedFrcCurve, method);
		if (result == null)
			return null;
		double fireNumber = result[0];
		double correlation = result[1];

		// The FRC paper states that the super-resolution pixel size should be smaller
		// than 1/4 of R (the resolution).
		boolean pixelsTooBig = (4 > fireNumber);

		// The FIRE number will be returned in pixels relative to the input images. 
		// However these were generated using an image scale so adjust for this.
		fireNumber *= images.nmPerPixel;

		if (pixelsTooBig)
		{
			// Q. Should this be output somewhere else?
			Utils.log(
					"%s Warning: The super-resolution pixel size (%s) should be smaller than 1/4 of R (the resolution %s)",
					TITLE, Utils.rounded(images.nmPerPixel), Utils.rounded(fireNumber));
		}

		return new FireResult(fireNumber, correlation, images.nmPerPixel, frcCurve, smoothedFrcCurve);
	}

	private void runQEstimation()
	{
		if (!showQEstimationInputDialog())
			return;

		MemoryPeakResults results = ResultsManager.loadInputResults(inputOption, false);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			return;
		}
		if (results.getCalibration() == null)
		{
			IJ.error(TITLE, "The results are not calibrated");
			return;
		}

		results = cropToRoi(results);
		if (results.size() < 2)
		{
			IJ.error(TITLE, "No results within the crop region");
			return;
		}

		initialise(results, null);

		//String name = results.getName();
		double fourierImageScale = SCALE_VALUES[imageScaleIndex];
		int imageSize = IMAGE_SIZE_VALUES[imageSizeIndex];

		// Create the image and compute the numerator of FRC. 
		// Do not use the signal so results.size() is the number of localisations.
		FireImages images = createImages(fourierImageScale, imageSize, false);
		FRC frc = new FRC();
		frc.progress = progress;
		frc.perimeterSamplingFactor = perimeterSamplingFactor;
		frc.useHalfCircle = useHalfCircle;
		double[][] frcCurve = frc.calculateFrcCurve(images.ip1, images.ip2);

		QPlot qplot = new QPlot(images.nmPerPixel, results.size(), frcCurve);

		// Build a histogram of the localisation precision.
		// Get the initial mean and SD and plot as a Gaussian.
		PrecisionHistogram histogram = calculatePrecisionHistogram();

		// Interactive dialog to estimate Q (blinking events per flourophore) using 
		// sliders for the mean and standard deviation of the localisation precision.
		showQEstimationDialog(histogram, qplot);
	}

	private boolean showQEstimationInputDialog()
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

		gd.addMessage("Estimate the blinking correction parameter Q for Fourier Ring Correlation");

		ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);

		//gd.addCheckbox("Use_signal (if present)", useSignal);
		gd.addNumericField("Max_per_bin", maxPerBin, 0);
		gd.addNumericField("Block_size", blockSize, 0);
		gd.addCheckbox("Random_split", randomSplit);
		gd.addMessage("Fourier options:");
		gd.addChoice("Fourier_image_scale", SCALE_ITEMS, SCALE_ITEMS[imageScaleIndex]);
		gd.addChoice("Auto_image_scale", IMAGE_SIZE_ITEMS, IMAGE_SIZE_ITEMS[imageSizeIndex]);
		gd.addSlider("Sampling_factor", 0.2, 4, perimeterSamplingFactor);
		gd.addCheckbox("Half_circle", useHalfCircle);
		if (!titles.isEmpty())
			gd.addCheckbox((titles.size() == 1) ? "Use_ROI" : "Choose_ROI", chooseRoi);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		inputOption = ResultsManager.getInputSource(gd);
		//useSignal = gd.getNextBoolean();
		maxPerBin = Math.abs((int) gd.getNextNumber());
		blockSize = Math.max(1, (int) gd.getNextNumber());
		randomSplit = gd.getNextBoolean();
		imageScaleIndex = gd.getNextChoiceIndex();
		imageSizeIndex = gd.getNextChoiceIndex();
		perimeterSamplingFactor = gd.getNextNumber();
		useHalfCircle = gd.getNextBoolean();

		// Check arguments
		try
		{
			Parameters.isAboveZero("Perimeter sampling factor", perimeterSamplingFactor);
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		if (!titles.isEmpty())
			chooseRoi = gd.getNextBoolean();

		if (!titles.isEmpty() && chooseRoi)
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

		return true;
	}

	public class QPlot
	{
		final double N;
		double[] vq, sinc, q;
		String title;

		// Store the last computed value
		double mean, sigma, qValue;

		QPlot(double nmPerPixel, double N, double[][] frcCurve)
		{
			this.N = N;

			// Compute v(q) - The numerator of the FRC divided by the number of pixels 
			// in the Fourier circle (2*pi*q*L)
			vq = new double[frcCurve.length];
			for (int i = 0; i < vq.length; i++)
			{
				// Use the actual number of samples and not 2*pi*q*L so we normalise 
				// to the correct scale 
				double d = frcCurve[i][2];
				//double d = 2.0 * Math.PI * i;
				vq[i] = frcCurve[i][3] / d;
			}

			q = FRC.computeQ(frcCurve, nmPerPixel);

			final double L = FRC.computeL(frcCurve);

			// Compute sinc factor
			sinc = new double[frcCurve.length];
			sinc[0] = 1; // By definition
			for (int i = 1; i < sinc.length; i++)

			{
				double d;
				// This should be sinc(pi*q*L)^2
				// where pi*q*L is half the number of samples in the Fourier circle.
				//d = Math.PI * i;

				// This does not seem correct. We are essentially computing Sine at intervals of pi
				// which all compute to 0. Try removing the L factor => sinc(pi*q)^2 with 
				// q == 1/L, 2/L, ... (i.e. no unit conversion to nm). This means that the function 
				// will start at 1 and drop off to zero at L.

				// sinc(pi*q)^2
				d = Math.PI * i / L;

				// Use q in the correct units - WRONG curve shape
				//d = Math.PI * q[i] * L;				

				sinc[i] = sinc(d);
				sinc[i] *= sinc[i];

				// TODO - figure out what this should be...
				// Ignore this term for now
				//sinc[i] = 1;
			}

			// For the plot
			title = results.getName() + " FRC Numerator Curve";
		}

		private double sinc(double x)
		{
			return FastMath.sin(x) / x;
		}

		final double LOG_10 = Math.log(10);

		double[] getLog(double[] hq)
		{
			double[] l = new double[hq.length];
			for (int q = 0; q < hq.length; q++)
				// Do this using log10
				l[q] = Math.log(Math.abs(vq[q] / hq[q] / sinc[q])) / LOG_10;
			return l;
		}

		double[] smooth(double[] l)
		{
			double bandwidth = 0.1;
			int robustness = 0;
			try
			{
				LoessInterpolator loess = new LoessInterpolator(bandwidth, robustness);
				return loess.smooth(q, l);
			}
			catch (Exception e)
			{
				e.printStackTrace();
			}
			return l;
		}

		PlotWindow plot(double mean, double sigma)
		{
			this.mean = mean;
			this.sigma = sigma;

			double[] hq = FRC.computeHq(q, mean, sigma);
			double[] l = getLog(hq);
			// Avoid bad value at zero
			l[0] = l[1];
			double[] sl = smooth(l);

			// log(NQ/4) = min of the curve => Q = 4*exp(min) / N
			// where N == Number of localisations
			double min = Maths.min(sl);
			// Use Math.pow since we are using log10. 
			qValue = 4 * Math.pow(10, min) / N;

			Plot2 plot = new Plot2(title, "Spatial Frequency (nm^-1)", "Log10 Scaled FRC Numerator");
			plot.setColor(Color.black);
			plot.addPoints(q, l, Plot.LINE);
			plot.addLabel(0, 0, String.format("Q = %.3f (Precision = %.3f +/- %.3f)", qValue, mean, sigma));
			plot.setColor(Color.red);
			plot.addPoints(q, sl, Plot.LINE);
			plot.setColor(Color.blue);
			plot.drawLine(0, min, q[q.length - 1], min);

			return Utils.display(title, plot);
		}
	}

	public class PrecisionHistogram
	{
		final float[] x, y;
		final String title;
		final double standardAmplitude;
		final float[] x2;

		double mean;
		double sigma;

		PrecisionHistogram(float[][] hist, int nPoints, String title)
		{
			this.title = title;
			x = Utils.createHistogramAxis(hist[0]);
			y = Utils.createHistogramValues(hist[1]);

			// Sum the area under the histogram to use for normalisation.
			// Amplitude = volume / (sigma * sqrt(2*pi)) 
			// Precompute the correct amplitude for a standard width Gaussian
			double dx = (hist[0][1] - hist[0][0]);
			standardAmplitude = nPoints * dx / Math.sqrt(2 * Math.PI);

			// Set up for drawing the Gaussian curve
			double min = x[0];
			double max = x[x.length - 1];
			int n = 100;
			dx = (max - min) / n;
			x2 = new float[n + 1];
			for (int i = 0; i <= n; i++)
				x2[i] = (float) (min + i * dx);
		}

		public PrecisionHistogram(String title)
		{
			this.title = title;
			// Set some defaults
			this.mean = 20;
			this.sigma = 2;
			x = y = x2 = null;
			standardAmplitude = 0;
		}

		PlotWindow plot(double mean, double sigma)
		{
			this.mean = mean;
			this.sigma = sigma;
			return plot();
		}

		PlotWindow plot()
		{
			Plot2 plot = new Plot2(title, "Precision (nm)", "Frequency");
			if (x != null)
			{
				plot.setColor(Color.black);
				plot.addPoints(x, y, Plot.LINE);
				plot.addLabel(0, 0, String.format("Precision = %.3f +/- %.3f", mean, sigma));
				// Add the Gaussian line
				// Compute the intergal of the standard gaussian between the min and max
				final double denom0 = 1.0 / (Math.sqrt(2.0) * sigma);
				double integral = 0.5 * Erf.erf((x2[0] - mean) * denom0, (x2[x2.length - 1] - mean) * denom0);
				// Normalise so the integral has the same volume as the histogram
				Gaussian g = new Gaussian(this.standardAmplitude / (sigma * integral), mean, sigma);
				float[] y2 = new float[x2.length];
				for (int i = 0; i < y2.length; i++)
				{
					y2[i] = (float) g.value(x2[i]);
				}
				// Normalise
				plot.setColor(Color.red);
				plot.addPoints(x2, y2, Plot.LINE);
				float max = Maths.max(y2);
				max = Maths.maxDefault(max, y);
				double rangex = 0; //(x2[x2.length - 1] - x2[0]) * 0.025;
				plot.setLimits(x2[0] - rangex, x2[x2.length - 1] + rangex, 0, max * 1.05);
			}
			else
			{
				// There is no base histogram.
				// Just plot a Gaussian +/- 4 SD.
				plot.addLabel(0, 0, String.format("Precision = %.3f +/- %.3f", mean, sigma));
				double min = Math.max(0, mean - 4 * sigma);
				double max = mean + 4 * sigma;
				int n = 100;
				double dx = (max - min) / n;
				float[] x2 = new float[n + 1];
				Gaussian g = new Gaussian(1, mean, sigma);
				float[] y2 = new float[x2.length];
				for (int i = 0; i <= n; i++)
				{
					x2[i] = (float) (min + i * dx);
					y2[i] = (float) g.value(x2[i]);
				}
				plot.setColor(Color.red);
				plot.addPoints(x2, y2, Plot.LINE);

				// Always put min = 0 otherwise the plot does not change.
				plot.setLimits(0, max, 0, 1.05);
			}
			return Utils.display(title, plot);
		}
	}

	/**
	 * Calculate the average precision by fitting a Gaussian to the histogram of the precision distribution.
	 * 
	 * @return The precision histogram
	 */
	public PrecisionHistogram calculatePrecisionHistogram()
	{
		boolean logFitParameters = false;
		String title = results.getName() + " Precision Histogram";

		// Check we can compute the precision for the results. We require that the widths and signal be valid and 
		// different for at least some of the localisations
		if (invalid(results))
		{
			return new PrecisionHistogram(title);
		}

		final double nmPerPixel = results.getNmPerPixel();
		final double gain = results.getGain();
		final boolean emCCD = results.isEMCCD();
		StoredDataStatistics precision = new StoredDataStatistics(results.size());
		for (PeakResult r : results.getResults())
		{
			precision.add(r.getPrecision(nmPerPixel, gain, emCCD));
		}
		//System.out.printf("Raw p = %f\n", precision.getMean());

		double yMin = Double.NEGATIVE_INFINITY, yMax = 0;

		// Set the min and max y-values using 1.5 x IQR 
		DescriptiveStatistics stats = precision.getStatistics();
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
				Utils.log("  Data range: %f - %f. Plotting 1.5x IQR: %f - %f", stats.getMin(), stats.getMax(), yMin,
						yMax);
		}

		if (yMin == Double.NEGATIVE_INFINITY)
		{
			int n = 5;
			yMin = Math.max(stats.getMin(), stats.getMean() - n * stats.getStandardDeviation());
			yMax = Math.min(stats.getMax(), stats.getMean() + n * stats.getStandardDeviation());

			if (logFitParameters)
				Utils.log("  Data range: %f - %f. Plotting mean +/- %dxSD: %f - %f", stats.getMin(), stats.getMax(), n,
						yMin, yMax);
		}

		// Get the data within the range
		double[] data = precision.getValues();
		precision = new StoredDataStatistics(data.length);
		for (double d : data)
		{
			if (d < yMin || d > yMax)
				continue;
			precision.add(d);
		}

		int histogramBins = Utils.getBins(precision, Utils.BinMethod.SCOTT);
		float[][] hist = Utils.calcHistogram(precision.getFloatValues(), yMin, yMax, histogramBins);

		PrecisionHistogram histogram = new PrecisionHistogram(hist, precision.getN(), title);

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
		if (logFitParameters)
			Utils.log("  Initial Statistics: %f +/- %f", stats2[0], stats2[1]);
		histogram.mean = stats2[0];
		histogram.sigma = stats2[1];

		// Standard Gaussian fit
		double[] parameters = fitGaussian(x, y);
		if (parameters == null)
		{
			Utils.log("  Failed to fit initial Gaussian");
			return histogram;
		}
		double newMean = parameters[1];
		double error = Math.abs(stats2[0] - newMean) / stats2[1];
		if (error > 3)
		{
			Utils.log("  Failed to fit Gaussian: %f standard deviations from histogram mean", error);
			return histogram;
		}
		if (newMean < yMin || newMean > yMax)
		{
			Utils.log("  Failed to fit Gaussian: %f outside data range %f - %f", newMean, yMin, yMax);
			return histogram;
		}

		if (logFitParameters)
			Utils.log("  Initial Gaussian: %f @ %f +/- %f", parameters[0], parameters[1], parameters[2]);

		histogram.mean = parameters[1];
		histogram.sigma = parameters[2];

		return histogram;
	}

	private boolean invalid(MemoryPeakResults results)
	{
		// Check all have a width and signal
		PeakResult[] data = results.toArray();
		for (int i = 0; i < data.length; i++)
		{
			PeakResult p = data[i];
			if (p.getSD() <= 0 || p.getSignal() <= 0)
				return true;
		}

		// Check for variable width that is not 1 and a variable signal
		for (int i = 0; i < data.length; i++)
		{
			PeakResult p = data[i];
			// Check this is valid
			if (p.getSD() != 1)
			{
				// Check the rest for a different value
				float w1 = p.getSD();
				float s1 = p.getSignal();
				for (int j = i + 1; j < data.length; j++)
				{
					PeakResult p2 = data[j];
					if (p2.getSD() != 1 && p2.getSD() != w1 && p.getSignal() != s1)
						return false;
				}
				// All the results are the same, this is not valid
				break;
			}
		}

		return true;
	}

	/**
	 * Fit gaussian.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @return new double[] { norm, mean, sigma }
	 */
	private double[] fitGaussian(float[] x, float[] y)
	{
		WeightedObservedPoints obs = new WeightedObservedPoints();
		for (int i = 0; i < x.length; i++)
			obs.add(x[i], y[i]);

		Collection<WeightedObservedPoint> observations = obs.toList();
		GaussianCurveFitter fitter = GaussianCurveFitter.create().withMaxIterations(2000);
		GaussianCurveFitter.ParameterGuesser guess = new GaussianCurveFitter.ParameterGuesser(observations);
		double[] initialGuess = null;
		try
		{
			initialGuess = guess.guess();
			return fitter.withStartPoint(initialGuess).fit(observations);
		}
		catch (TooManyEvaluationsException e)
		{
			// Use the initial estimate
			return initialGuess;
		}
		catch (Exception e)
		{
			// Just in case there is another exception type, or the initial estimate failed
			return null;
		}
	}

	private boolean showQEstimationDialog(final PrecisionHistogram histogram, final QPlot qplot)
	{
		WindowOrganiser wo = new WindowOrganiser();
		wo.add(histogram.plot());
		wo.add(qplot.plot(histogram.mean, histogram.sigma));
		wo.tile();

		NonBlockingGenericDialog gd = new NonBlockingGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Estimate the blinking correction parameter Q for Fourier Ring Correlation\n \n" +
				String.format("Precision estimate = %.3f +/- %.3f", histogram.mean, histogram.sigma));

		double mean10 = histogram.mean * 10;
		double sd10 = histogram.sigma * 10;
		gd.addSlider("Mean (x10)", Math.max(0, mean10 - sd10 * 2), mean10 + sd10 * 2, mean10);
		gd.addSlider("SD (x10)", Math.max(0, sd10 / 2), sd10 * 2, sd10);
		gd.addCheckbox("reset", false);

		// - create a synchronised work queue and pass it to the dialog listener.
		// - create a worker thread to take the work from the queue and do the computation.
		ArrayList<WorkStack> stacks = new ArrayList<WorkStack>();
		ArrayList<Worker> workers = new ArrayList<Worker>();
		add(stacks, workers, new Worker()
		{
			@Override
			void createResult(Work work)
			{
				histogram.plot(work.mean, work.sigma);
			}
		});
		add(stacks, workers, new Worker()
		{
			@Override
			void createResult(Work work)
			{
				qplot.plot(work.mean, work.sigma);
			}
		});

		gd.addDialogListener(new FIREDialogListener(gd, histogram, stacks));
		gd.showDialog();

		// Finish the worker threads
		for (WorkStack stack : stacks)
			stack.close();

		if (gd.wasCanceled())
			return false;

		// Store the Q value and the mean and sigma
		qValue = qplot.qValue;
		mean = qplot.mean;
		sigma = qplot.sigma;

		return true;
	}

	private void add(ArrayList<WorkStack> stacks, ArrayList<Worker> workers, Worker worker)
	{
		WorkStack stack = new WorkStack();
		worker.inbox = stack;
		stacks.add(stack);
		Thread t = new Thread(worker);
		t.setDaemon(true);
		t.start();
	}

	private class Work
	{
		long time = 0;
		double mean, sigma;

		Work(double mean, double sigma)
		{
			this.mean = mean;
			this.sigma = sigma;
		}
	}

	/**
	 * Allow work to be added to a FIFO stack in a synchronised manner
	 * 
	 * @author Alex Herbert
	 */
	private class WorkStack
	{
		// We only support a stack size of 1
		private Work work = null;

		synchronized void addWork(Work work)
		{
			this.work = work;
			this.notify();
		}

		synchronized void close()
		{
			this.work = null;
			this.notify();
		}

		synchronized Work getWork()
		{
			Work work = this.work;
			this.work = null;
			return work;
		}

		boolean isEmpty()
		{
			return work == null;
		}
	}

	private abstract class Worker implements Runnable
	{
		private boolean running = true;
		private Work lastWork = null;
		private WorkStack inbox;

		public void run()
		{
			while (running)
			{
				try
				{
					Work work = null;
					synchronized (inbox)
					{
						if (inbox.isEmpty())
						{
							debug("Inbox empty, waiting ...");
							inbox.wait();
						}
						work = inbox.getWork();
						if (work != null)
							debug(" Found work");
					}
					if (work == null)
					{
						debug(" No work, stopping");
						break;
					}

					// Delay processing the work. Allows the work to be updated before we process it.
					if (work.time != 0)
					{
						debug(" Checking delay");
						long time = work.time;
						while (System.currentTimeMillis() < time)
						{
							debug(" Delaying");
							Thread.sleep(50);
							// Assume new work can be added to the inbox. Here we are peaking at the inbox
							// so we do not take ownership with synchronized
							if (inbox.work != null)
								time = inbox.work.time;
						}
						// If we intend to modify the inbox then we should take ownership
						synchronized (inbox)
						{
							if (!inbox.isEmpty())
							{
								work = inbox.getWork();
								debug(" Found updated work");
							}
						}
					}

					if (!equals(work, lastWork))
						createResult(work);
					lastWork = work;
				}
				catch (InterruptedException e)
				{
					debug(" Interrupted, stopping");
					break;
				}
			}
		}

		private void debug(String msg)
		{
			boolean debug = false;
			if (debug)
				System.out.println(this.getClass().getSimpleName() + msg);
		}

		boolean equals(Work work, Work lastWork)
		{
			if (lastWork == null)
				return false;

			if (work.mean != lastWork.mean)
				return false;
			if (work.sigma != lastWork.sigma)
				return false;

			return true;
		}

		abstract void createResult(Work work);
	}

	private class FIREDialogListener implements DialogListener
	{
		/**
		 * Delay (in milliseconds) used when entering new values in the dialog before the preview is processed
		 */
		static final long DELAY = 500;
		
		long time;
		boolean notActive = true;
		volatile int ignore = 0;
		ArrayList<WorkStack> stacks;
		double defaultMean, defaultSigma;
		String m, s;
		TextField tf1, tf2;
		Checkbox cb;
		final boolean isMacro;

		FIREDialogListener(GenericDialog gd, PrecisionHistogram histogram, ArrayList<WorkStack> stacks)
		{
			time = System.currentTimeMillis() + 1000;
			this.stacks = stacks;
			this.defaultMean = histogram.mean;
			this.defaultSigma = histogram.sigma;
			isMacro = Utils.isMacro();
			// For the reset
			tf1 = (TextField) gd.getNumericFields().get(0);
			tf2 = (TextField) gd.getNumericFields().get(1);
			cb = (Checkbox) (gd.getCheckboxes().get(0));
			m = tf1.getText();
			s = tf2.getText();
		}

		public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
		{
			if (notActive && System.currentTimeMillis() < time)
				return true;
			if (ignore-- > 0)
			{
				//System.out.println("ignored");
				return true;
			}

			notActive = false;

			double mean = Math.abs(gd.getNextNumber()) / 10;
			double sigma = Math.abs(gd.getNextNumber()) / 10;
			boolean reset = gd.getNextBoolean();

			System.out.printf("Event: %s, %f, %f\n", e, mean, sigma);

			// Allow reset to default
			if (reset)
			{
				// This does not trigger the event
				cb.setState(false);
				mean = this.defaultMean;
				sigma = this.defaultSigma;
			}

			Work work = new Work(mean, sigma);
			// Implement a delay to allow typing.
			// This is also applied to the sliders which we do not want. 
			// Ideally we would have no delay for sliders (since they are in the correct place already
			// but a delay for typing the in the text field). Unfortunately the AWTEvent raised by ImageJ
			// for the slider is actually from the TextField so we cannot tell the difference.
			// For now just have no delay.
			//work.time = (isMacro) ? 0 : System.currentTimeMillis() + DELAY;

			// Offload this work onto a thread that just picks up the most recent dialog input.
			for (WorkStack stack : stacks)
				stack.addWork(work);

			if (reset)
			{
				// These trigger dialogItemChanged(...) so do them after we added 
				// work to the queue and ignore the events
				ignore = 2;
				tf1.setText(m);
				tf2.setText(s);
			}

			return true;
		}
	}
}
