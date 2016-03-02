package gdsc.smlm.ij.plugins;

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
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Plot2;
import ij.plugin.PlugIn;
import ij.plugin.frame.Recorder;
import ij.process.ImageProcessor;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

/**
 * Computes the Fourier Image Resolution of an image
 */
public class FIRE implements PlugIn
{
	private static String TITLE = "Fourier Image REsolution (FIRE)";
	private static String inputOption = "";

	private static boolean randomSplit = true;
	private static String[] SCALE_ITEMS = new String[] { "Auto", "1", "2", "4", "8", "16", "32", "64", "128" };
	private static int[] SCALE_VALUES = new int[] { 0, 1, 2, 4, 8, 16, 32, 64, 128 };
	private static String[] IMAGE_SIZE_ITEMS = new String[] { "2048", "4096", "8192", "16384" };
	private static int[] IMAGE_SIZE_VALUES = new int[] { 2047, 4095, 8191, 16383 };
	private static int imageScaleIndex = 0;
	private static int imageSizeIndex = 0;
	private static double perimeterSamplingFactor = 1;
	private static boolean useHalfCircle = true;
	private static int thresholdMethodIndex = 0;
	private static boolean showFRCCurve = true;
	private static boolean showFRCTimeEvolution = false;

	private static boolean chooseRoi = false;
	private static String roiImage = "";

	private Rectangle roiBounds;
	private int roiImageWidth, roiImageHeight;

	MemoryPeakResults results;
	ImageProcessor ip1, ip2;
	float imageScale;
	double[][] frcCurve;
	double[][] smoothedFrcCurve;
	String units;
	double nmPerPixel = 1;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		PluginTracker.recordPlugin(this.getClass(), arg);
		
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

		results = cropToRoi(results);
		if (results.size() == 0)
		{
			IJ.error(TITLE, "No results within the crop region");
			IJ.showStatus("");
			return;
		}

		long start = System.currentTimeMillis();

		ThresholdMethod method = FRC.ThresholdMethod.values()[thresholdMethodIndex];

		// Compute FIRE
		initialise(results);
		double fire = calculateFireNumber(method, SCALE_VALUES[imageScaleIndex], IMAGE_SIZE_VALUES[imageSizeIndex]);

		IJ.log(String.format("%s : FIRE number = %s %s (Fourier scale = %s)", results.getName(),
				Utils.rounded(fire, 4), units, Utils.rounded(imageScale, 3)));

		String name = results.getName();
		if (showFRCCurve)
			showFrcCurve(name, frcCurve, smoothedFrcCurve, method);

		if (showFRCTimeEvolution && !Double.isNaN(fire))
			showFrcTimeEvolution(name, fire, method);

		double seconds = (System.currentTimeMillis() - start) / 1000.0;
		IJ.showStatus(TITLE + " complete : " + seconds + "s");
	}

	private MemoryPeakResults cropToRoi(MemoryPeakResults results)
	{
		if (roiBounds == null)
			return results;

		// Adjust bounds relative to input results image
		double xscale = roiImageWidth / results.getBounds().width;
		double yscale = roiImageHeight / results.getBounds().height;
		roiBounds.x /= xscale;
		roiBounds.width /= xscale;
		roiBounds.y /= yscale;
		roiBounds.height /= yscale;

		float minX = (int) (roiBounds.x);
		float maxX = (int) Math.ceil(roiBounds.x + roiBounds.width);
		float minY = (int) (roiBounds.y);
		float maxY = (int) Math.ceil(roiBounds.y + roiBounds.height);

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

		gd.addCheckbox("Random_split", randomSplit);
		gd.addChoice("Fourier_image_scale", SCALE_ITEMS, SCALE_ITEMS[imageScaleIndex]);
		gd.addChoice("Auto_image_scale", IMAGE_SIZE_ITEMS, IMAGE_SIZE_ITEMS[imageSizeIndex]);
		gd.addSlider("Sampling_factor", 0.2, 4, perimeterSamplingFactor);
		gd.addCheckbox("Half_circle", useHalfCircle);
		String[] methodNames = SettingsManager.getNames((Object[]) FRC.ThresholdMethod.values());
		gd.addChoice("Threshold_method", methodNames, methodNames[thresholdMethodIndex]);
		gd.addCheckbox("Show_FRC_curve", showFRCCurve);
		gd.addCheckbox("Show_FRC_time_evolution", showFRCTimeEvolution);
		if (!titles.isEmpty())
			gd.addCheckbox((titles.size() == 1) ? "Use_ROI" : "Choose_ROI", chooseRoi);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		inputOption = ResultsManager.getInputSource(gd);
		randomSplit = gd.getNextBoolean();
		imageScaleIndex = gd.getNextChoiceIndex();
		imageSizeIndex = gd.getNextChoiceIndex();
		perimeterSamplingFactor = gd.getNextNumber();
		useHalfCircle = gd.getNextBoolean();
		thresholdMethodIndex = gd.getNextChoiceIndex();
		showFRCCurve = gd.getNextBoolean();
		showFRCTimeEvolution = gd.getNextBoolean();

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

	public void initialise(MemoryPeakResults results)
	{
		this.results = results;
		ip1 = ip2 = null;
		nmPerPixel = 1;
		units = "px";
		if (results.getCalibration() != null)
		{
			nmPerPixel = results.getNmPerPixel();
			units = "nm";
		}
	}

	public ImageProcessor[] createImages(float fourierImageScale, int imageSize)
	{
		ip1 = ip2 = null;
		if (results == null || results.size() == 0)
			return null;

		ResultsImage imageType = (results.getResults().get(0).getSignal() > 0) ? ResultsImage.SIGNAL_INTENSITY
				: ResultsImage.LOCALISATIONS;

		// Draw images using the existing IJ routines.
		// TODO - This could be speeded up using a simple 2D-histogram
		Rectangle bounds = results.getBounds(true);
		boolean weighted = true;
		boolean equalised = false;
		if (fourierImageScale == 0)
			this.imageScale = imageSize / (float) FastMath.max(bounds.x + bounds.width, bounds.y + bounds.height);
		else
			// TODO - The coordinates should be adjusted if the max-min can fit inside a smaller power of 2
			this.imageScale = fourierImageScale;

		IJImagePeakResults image1 = ImagePeakResultsFactory.createPeakResultsImage(imageType, weighted, equalised,
				"IP1", bounds, results.getNmPerPixel(), results.getGain(), imageScale, 0, ResultsMode.ADD);
		image1.setDisplayImage(false);
		image1.begin();

		IJImagePeakResults image2 = ImagePeakResultsFactory.createPeakResultsImage(imageType, weighted, equalised,
				"IP2", bounds, results.getNmPerPixel(), results.getGain(), imageScale, 0, ResultsMode.ADD);
		image2.setDisplayImage(false);
		image2.begin();

		List<PeakResult> list = results.getResults();
		if (randomSplit)
		{
			@SuppressWarnings("unchecked")
			ArrayList<PeakResult> dest = (ArrayList<PeakResult>) ((ArrayList<PeakResult>) list).clone();
			Collections.shuffle(dest);
			list = dest;
		}

		// Split alternating
		int i = 0;
		for (PeakResult p : list)
		{
			if (i++ % 2 == 0)
				image1.add(p.peak, p.origX, p.origY, p.origValue, p.error, p.noise, p.params, p.paramsStdDev);
			else
				image2.add(p.peak, p.origX, p.origY, p.origValue, p.error, p.noise, p.params, p.paramsStdDev);
		}

		image1.end();
		ip1 = image1.getImagePlus().getProcessor();

		image2.end();
		ip2 = image2.getImagePlus().getProcessor();

		return new ImageProcessor[] { ip1, ip2 };
	}

	private void showFrcCurve(String name, double[][] frcIn, double[][] frcNoSmooth, ThresholdMethod method)
	{
		double[][] frcCurve = frcIn;

		FRC frc = new FRC();
		double[] threshold = frc.calculateThresholdCurve(frcCurve, method);

		double[] xValues = new double[frcCurve.length];
		double[] yValues = new double[frcCurve.length];
		double[] yValuesNotSmooth = new double[frcCurve.length];

		double nmPerPixel = 1;
		String units = "px";
		if (results.getCalibration() != null)
		{
			nmPerPixel = results.getNmPerPixel();
			units = "nm";
		}

		double yMin = frcCurve[0][1];
		double yMax = frcCurve[0][1];
		// Since the Fourier calculation only uses half of the image (from centre to the edge) 
		// we must double the curve length to get the original maximum image width. In addition
		// the computation was up to the edge-1 pixels so add back a pixel to the curve length.
		double frcCurveLength = (frcCurve[(frcCurve.length - 1)][0] + 1) * 2.0;
		double conversion = imageScale / (frcCurveLength * nmPerPixel);
		for (int i = 0; i < xValues.length; i++)
		{
			final double radius = frcCurve[i][0];
			xValues[i] = radius * conversion;
			yValues[i] = frcCurve[i][1];

			yMin = FastMath.min(yMin, yValues[i]);
			yMax = FastMath.max(yMax, yValues[i]);
			if (frcNoSmooth != null)
				yValuesNotSmooth[i] = frcNoSmooth[i][1];
		}

		String title = name + " FRC Curve";
		double[] dummy = null;
		Plot2 plot = new Plot2(title, String.format("Spatial Frequency (%s^-1)", units), "FRC", dummy, dummy);
		plot.setLimits(0, xValues[xValues.length - 1], yMin, yMax);

		plot.setColor(Color.BLACK);
		plot.addPoints(xValues, yValues, Plot2.LINE);

		plot.setColor(Color.BLUE);
		plot.addPoints(xValues, threshold, Plot2.LINE);

		if (frcNoSmooth != null)
		{
			plot.setColor(Color.RED);
			plot.addPoints(xValues, yValuesNotSmooth, Plot2.LINE);
		}

		Utils.display(title, plot);
	}

	private void showFrcTimeEvolution(String name, double fireNumber, ThresholdMethod method)
	{
		IJ.showStatus("Calculating FRC time evolution curve...");

		results.sort();
		List<PeakResult> list = results.getResults();

		int nSteps = 10;
		int maxT = list.get(list.size() - 1).peak;
		if (maxT == 0)
			maxT = list.size();
		int step = maxT / nSteps;

		ArrayList<Double> x = new ArrayList<Double>();
		ArrayList<Double> y = new ArrayList<Double>();

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
			f.initialise(newResults);
			double fire = f.calculateFireNumber(method, imageScale, ip1.getWidth() - 1);
			y.add(fire);

			yMin = FastMath.min(yMin, fire);
			yMax = FastMath.max(yMax, fire);
		}

		// Add the final fire number
		x.add((double) maxT);
		y.add(fireNumber);

		double[] xValues = new double[x.size()];
		double[] yValues = new double[xValues.length];
		for (int ii = 0; ii < xValues.length; ii++)
		{
			xValues[ii] = x.get(ii);
			yValues[ii] = y.get(ii);
		}

		String title = name + " FRC Time Evolution";
		Plot2 plot = new Plot2(title, "Frames", "Resolution", (float[]) null, (float[]) null);
		plot.setLimits(xValues[0], xValues[xValues.length - 1], yMin, yMax);
		plot.addPoints(xValues, yValues, 1);

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
	public double calculateFireNumber(ThresholdMethod method, float fourierImageScale, int imageSize)
	{
		if (ip1 == null || ip2 == null)
		{
			if (createImages(fourierImageScale, imageSize) == null)
				return 0;
		}
		FRC frc = new FRC();
		frc.perimeterSamplingFactor = perimeterSamplingFactor;
		frc.useHalfCircle = useHalfCircle;
		frcCurve = frc.calculateFrcCurve(ip1, ip2);
		smoothedFrcCurve = frc.getSmoothedCurve(frcCurve);
		// The FIRE number will be returned in pixels relative to the input images. 
		// However these were generated using an image scale so adjust for this.

		return nmPerPixel * frc.calculateFireNumber(smoothedFrcCurve, method) / imageScale;
	}
}
