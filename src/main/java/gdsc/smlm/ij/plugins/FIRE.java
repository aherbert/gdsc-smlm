package gdsc.smlm.ij.plugins;

import java.awt.Color;
import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.MathArrays;

import gdsc.core.ij.IJTrackProgress;
import gdsc.core.ij.Utils;
import gdsc.core.logging.NullTrackProgress;
import gdsc.core.logging.TrackProgress;
import gdsc.core.utils.Maths;
import gdsc.core.utils.Statistics;
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
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.Plot2;
import ij.plugin.PlugIn;
import ij.plugin.WindowOrganiser;
import ij.plugin.frame.Recorder;
import ij.process.ImageProcessor;
import ij.process.LUT;
import ij.process.LUTHelper;
import ij.process.LUTHelper.LutColour;

/**
 * Computes the Fourier Image Resolution of an image
 */
public class FIRE implements PlugIn
{
	private static String TITLE = "Fourier Image REsolution (FIRE)";
	private static String inputOption = "";
	private static String inputOption2 = "";

	private static int repeats = 1;
	private static boolean useSignal = true;
	private static int maxPerBin = 5;
	private static boolean randomSplit = true;
	private static int blockSize = 50;
	private static String[] SCALE_ITEMS;
	private static int[] SCALE_VALUES = new int[] { 0, 1, 2, 4, 8, 16, 32, 64, 128 };
	private static String[] IMAGE_SIZE_ITEMS;
	private static int[] IMAGE_SIZE_VALUES;
	private static int imageScaleIndex = 0;
	private static int imageSizeIndex;

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

	private Rectangle roiBounds;
	private int roiImageWidth, roiImageHeight;

	// Stored in initialisation
	MemoryPeakResults results, results2;
	Rectangle2D dataBounds;
	String units;
	double nmPerPixel = 1;

	public class FireImages
	{
		final ImageProcessor ip1, ip2;
		final double imageScale;

		FireImages(ImageProcessor ip1, ImageProcessor ip2, double imageScale)
		{
			this.ip1 = ip1;
			this.ip2 = ip2;
			this.imageScale = imageScale;
		}
	}

	public class FireResult
	{
		final double fireNumber;
		final double imageScale;
		final double[][] frcCurve;
		final double[][] smoothedFrcCurve;

		FireResult(double fireNumber, double imageScale, double[][] frcCurve, double[][] smoothedFrcCurve)
		{
			this.fireNumber = fireNumber;
			this.imageScale = imageScale;
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
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		// Require some fit results and selected regions
		int size = MemoryPeakResults.countMemorySize();
		if (size == 0)
		{
			IJ.error(TITLE, "There are no fitting results in memory");
			return;
		}

		if (!showDialog())
			return;

		MemoryPeakResults results = ResultsManager.loadInputResults(inputOption, false);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			IJ.showStatus("");
			return;
		}
		MemoryPeakResults results2 = ResultsManager.loadInputResults(inputOption2, false);

		results = cropToRoi(results);
		if (results.size() == 0)
		{
			IJ.error(TITLE, "No results within the crop region");
			IJ.showStatus("");
			return;
		}
		if (results2 != null)
		{
			results2 = cropToRoi(results2);
			if (results2.size() == 0)
			{
				IJ.error(TITLE, "No results2 within the crop region");
				IJ.showStatus("");
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
						Utils.rounded(result.fireNumber, 4), units, Utils.rounded(result.imageScale, 3)));

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
							Utils.rounded(result.fireNumber, 4), units, Utils.rounded(result.imageScale, 3)));

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
								Utils.rounded(result.fireNumber, 4), units, Utils.rounded(result.imageScale, 3)));
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
					IJ.log(String.format("%s : FIRE number = %s +/- %s %s (Fourier scale = %s)", name,
							Utils.rounded(stats.getMean(), 4), Utils.rounded(stats.getStandardDeviation(), 4), units,
							Utils.rounded(result.imageScale, 3)));
					if (showFRCCurve)
					{
						Plot2 plot = curve.getPlot();
						Utils.display(plot.getTitle(), plot);
					}
				}
			}

			// Only do this once
			if (showFRCTimeEvolution && result != null && !Double.isNaN(result.fireNumber))
				showFrcTimeEvolution(name, result.fireNumber, method, result.imageScale, imageSize);
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
		if (results == null)
			return null;

		final boolean hasSignal = useSignal && (results.getResults().get(0).getSignal() > 0);

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
			// TODO - The coordinates should be adjusted if the max-min can fit inside a smaller power of 2
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

			for (int b = 0; b < nBlocks; b++)
			{
				// Split alternating
				IJImagePeakResults image = (b % 2 == 0) ? image1 : image2;
				for (PeakResult p : blocks[indices[b]])
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

		return new FireImages(ip1, ip2, imageScale);
	}

	/**
	 * Encapsulate plotting the FRC curve to allow multiple curves to be plotted together
	 */
	private class FrcCurve
	{
		double[] xValues = null;
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
				double nmPerPixel = 1;
				String units = "px";
				if (results.getCalibration() != null)
				{
					nmPerPixel = results.getNmPerPixel();
					units = "nm";
				}
				String title = name + " FRC Curve";
				plot = new Plot2(title, String.format("Spatial Frequency (%s^-1)", units), "FRC");

				xValues = new double[frcCurve.length];
				// Since the Fourier calculation only uses half of the image (from centre to the edge) 
				// we must double the curve length to get the original maximum image width. In addition
				// the computation was up to the edge-1 pixels so add back a pixel to the curve length.
				double frcCurveLength = (frcCurve[(frcCurve.length - 1)][0] + 1) * 2.0;
				double conversion = result.imageScale / (frcCurveLength * nmPerPixel);
				for (int i = 0; i < xValues.length; i++)
				{
					final double radius = frcCurve[i][0];
					xValues[i] = radius * conversion;
				}

				// The threshold curve is the same
				double[] threshold = FRC.calculateThresholdCurve(frcCurve, method);
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

		void add(Color color, double[] y)
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
			plot.setLimits(limits[0], limits[1], limits[2], limits[3]);
			return plot;
		}
	}

	private Plot2 createFrcCurve(String name, FireResult result, ThresholdMethod method)
	{
		FrcCurve curve = new FrcCurve();
		curve.add(name, result, method, Color.red, Color.blue, Color.black);
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

			FIRE f = new FIRE();
			f.initialise(newResults, null);
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
		double imageScale = images.imageScale;
		double[][] frcCurve = frc.calculateFrcCurve(images.ip1, images.ip2);
		double[][] smoothedFrcCurve = frc.getSmoothedCurve(frcCurve);
		// The FIRE number will be returned in pixels relative to the input images. 
		// However these were generated using an image scale so adjust for this.

		double fireNumber = nmPerPixel * frc.calculateFireNumber(smoothedFrcCurve, method) / imageScale;

		return new FireResult(fireNumber, imageScale, frcCurve, smoothedFrcCurve);
	}
}
