package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.Color;
import java.awt.Frame;
import java.awt.Point;
import java.awt.Rectangle;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.FastMath;

import gdsc.core.data.utils.Rounder;
import gdsc.core.data.utils.RounderFactory;
import gdsc.core.ij.AlignImagesFFT;
import gdsc.core.ij.AlignImagesFFT.SubPixelMethod;
import gdsc.core.ij.AlignImagesFFT.WindowMethod;
import gdsc.core.ij.IJTrackProgress;
import gdsc.core.ij.Utils;
import gdsc.core.match.BasePoint;
import gdsc.core.math.interpolation.CubicSplinePosition;
import gdsc.core.math.interpolation.CustomTricubicFunction;
import gdsc.core.math.interpolation.CustomTricubicInterpolatingFunction;
import gdsc.core.math.interpolation.CustomTricubicInterpolator;
import gdsc.core.math.interpolation.IndexedCubicSplinePosition;
import gdsc.core.utils.DoubleData;
import gdsc.core.utils.ImageExtractor;
import gdsc.core.utils.ImageWindow;
import gdsc.core.utils.Maths;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.core.utils.Sort;
import gdsc.core.utils.Statistics;
import gdsc.core.utils.StoredData;
import gdsc.core.utils.StoredDataStatistics;
import gdsc.core.utils.TurboList;
import gdsc.smlm.data.config.FitProtos.FitSolver;
import gdsc.smlm.data.config.FitProtosHelper;
import gdsc.smlm.data.config.PSFProtos.ImagePSF;
import gdsc.smlm.data.config.PSFProtos.PSF;
import gdsc.smlm.data.config.PSFProtos.PSFParameter;
import gdsc.smlm.data.config.PSFProtos.PSFParameterUnit;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import gdsc.smlm.engine.FitConfiguration;

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

import gdsc.smlm.engine.FitEngine;
import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.engine.FitParameters;
import gdsc.smlm.engine.FitQueue;
import gdsc.smlm.engine.ParameterisedFitJob;
import gdsc.smlm.ij.settings.ImagePSFHelper;
import gdsc.smlm.ij.utils.ImageConverter;
import gdsc.smlm.results.Counter;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.SynchronizedPeakResults;
import gdsc.smlm.results.procedures.HeightResultProcedure;
import gdsc.smlm.results.procedures.PeakResultProcedure;
import gdsc.smlm.results.procedures.WidthResultProcedure;
import gnu.trove.list.array.TDoubleArrayList;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.ExtendedGenericDialog;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.ImageRoi;
import ij.gui.Line;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.Plot2;
import ij.gui.PlotWindow;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.io.FileInfo;
import ij.measure.Calibration;
import ij.plugin.filter.PlugInFilter;
import ij.process.Blitter;
import ij.process.ByteProcessor;
import ij.process.FloatPolygon;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.LUTHelper;
import ij.process.LUTHelper.LutColour;

/**
 * Produces an average PSF image using selected diffraction limited spots from a sample image.
 * <p>
 * The input image must be a z-stack of diffraction limited spots for example quantum dots or fluorescent beads. Spots
 * will be used only when there are no spots within a specified distance to ensure a clean signal is extracted.
 */
public class PSFCreator implements PlugInFilter
{
	private final static String TITLE = "PSF Creator";
	private final static String TITLE_AMPLITUDE = "Spot Amplitude";
	private final static String TITLE_PSF_PARAMETERS = "Spot PSF";
	private final static String TITLE_INTENSITY = "Spot Intensity";
	private final static String TITLE_WINDOW = "PSF Window Function";
	private final static String TITLE_BACKGROUND = "PSF Cumulative Intensity";

	private final static String[] MODE = { "Projection", "Gaussian Fitting" };
	private static int mode = 0;

	private static double nmPerSlice = 0;
	private static double radius = 10;
	private static double amplitudeFraction = 0.2;
	private static int startBackgroundFrames = 5;
	private static int endBackgroundFrames = 5;
	private static int magnification = 10;
	private static double smoothing = 0.25;
	private static boolean centreEachSlice = false;
	private static double comCutOff = 5e-2;
	private static boolean interactiveMode = false;
	private static int interpolationMethod = ImageProcessor.BICUBIC;

	private int flags = DOES_16 | DOES_8G | DOES_32 | NO_CHANGES;
	private ImagePlus imp, psfImp;

	private static double nmPerPixel;
	private static int projectionMagnification = 2;
	private static int maxIterations = 20;
	private static int psfMagnification = 4;
	private static double backgroundCutoff = 2;
	private float background = 0;
	private static double windowAlpha = 0.1;
	private static Rounder rounder = RounderFactory.create(4);

	private FitEngineConfiguration config = null;
	private FitConfiguration fitConfig;
	private int boxRadius;
	private static Point yesNoPosition = null;

	private ExecutorService threadPool = null;
	private double progress = 0;

	// Private variables that are used during background threaded plotting of the cumulative signal 
	private ImageStack psf = null;
	private ImagePlus[] psfOut = null;
	private int zCentre = 0;
	private double psfWidth = 0;
	private double psfNmPerPixel = 0;

	// Amplitude plot
	private double[] z = null;
	private double[] a;
	private double[] smoothAz;
	private double[] smoothA;

	// PSF plot
	private double[] xCoord = null;
	private double[] yCoord;
	private double[] sd;
	private double[] newZ;
	private double[] smoothX;
	private double[] smoothY;
	private double[] smoothSd;

	// % PSF Signal plot
	private double[] signalZ = null;
	private double[] signal = null;
	private String signalTitle = null;
	private double[] signalLimits = null;

	// Cumulative signal plot
	private int[] indexLookup = null;
	private double[] distances = null;
	private double maxCumulativeSignal = 1;
	private int slice = 0;
	private double distanceThreshold = 0;
	private boolean normalise = false;
	private boolean resetScale = true;

	private boolean plotLock1 = false;
	private boolean plotLock2 = false;
	private boolean plotLock3 = false;
	private boolean plotLock4 = false;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	public int setup(String arg, ImagePlus imp)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (imp == null)
		{
			IJ.noImage();
			return DONE;
		}

		Roi roi = imp.getRoi();
		if (roi == null || roi.getType() != Roi.POINT)
		{
			IJ.error("Point ROI required");
			return DONE;
		}

		this.imp = imp;

		return showDialog();
	}

	private int showDialog()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		guessScale();

		gd.addMessage("Produces an average PSF using selected diffraction limited spots.");

		gd.addChoice("Mode", MODE, mode);
		//		gd.addNumericField("nm_per_slice", nmPerSlice, 0);
		gd.addSlider("Radius", 3, 20, radius);
		//		gd.addSlider("Amplitude_fraction", 0.01, 0.5, amplitudeFraction);
		//		gd.addSlider("Start_background_frames", 1, 20, startBackgroundFrames);
		//		gd.addSlider("End_background_frames", 1, 20, endBackgroundFrames);
		//		gd.addSlider("Magnification", 5, 15, magnification);
		//		gd.addSlider("Smoothing", 0.25, 0.5, smoothing);
		//		gd.addCheckbox("Centre_each_slice", centreEachSlice);
		//		gd.addNumericField("CoM_cut_off", comCutOff, -2);
		gd.addCheckbox("Interactive_mode", interactiveMode);
		//		String[] methods = ImageProcessor.getInterpolationMethods();
		//		gd.addChoice("Interpolation", methods, methods[interpolationMethod]);

		gd.showDialog();

		if (gd.wasCanceled())
			return DONE;

		mode = gd.getNextChoiceIndex();
		//		nmPerSlice = gd.getNextNumber();
		radius = gd.getNextNumber();
		//		amplitudeFraction = gd.getNextNumber();
		//		startBackgroundFrames = (int) gd.getNextNumber();
		//		endBackgroundFrames = (int) gd.getNextNumber();
		//		magnification = (int) gd.getNextNumber();
		//		smoothing = gd.getNextNumber();
		//		centreEachSlice = gd.getNextBoolean();
		//		comCutOff = Maths.max(0, gd.getNextNumber());
		interactiveMode = gd.getNextBoolean();
		//		interpolationMethod = gd.getNextChoiceIndex();

		// Check arguments
		try
		{
			//			Parameters.isPositive("nm/slice", nmPerSlice);
			Parameters.isAbove("Radius", radius, 2);
			//			Parameters.isAbove("Amplitude fraction", amplitudeFraction, 0.01);
			//			Parameters.isBelow("Amplitude fraction", amplitudeFraction, 0.9);
			//			Parameters.isPositive("Start background frames", startBackgroundFrames);
			//			Parameters.isPositive("End background frames", endBackgroundFrames);
			//			Parameters.isAbove("Total background frames", startBackgroundFrames + endBackgroundFrames, 1);
			//			Parameters.isAbove("Magnification", magnification, 1);
			//			Parameters.isAbove("Smoothing", smoothing, 0);
			//			Parameters.isBelow("Smoothing", smoothing, 1);
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return DONE;
		}

		return flags;
	}

	private void guessScale()
	{
		if (nmPerPixel == 0 || nmPerSlice == 0)
		{
			Calibration c = imp.getCalibration();
			nmPerPixel = guessScale(c.getXUnit(), c.pixelWidth);
			nmPerSlice = guessScale(c.getZUnit(), c.pixelDepth);
		}
	}

	private double guessScale(String unit, double units)
	{
		unit = unit.toLowerCase();
		if (unit.equals("nm") || unit.startsWith("nanomet"))
			return units;
		if (unit.equals("\u00B5m") || // Sanitised version of um
				unit.startsWith("micron"))
			return units * 1000;
		return 0;
	}

	private boolean showFittingDialog()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Use PSF fitting to create a combined PSF");

		gd.addNumericField("nm_per_slice", nmPerSlice, 0);
		gd.addSlider("Amplitude_fraction", 0.01, 0.5, amplitudeFraction);
		gd.addSlider("Start_background_frames", 1, 20, startBackgroundFrames);
		gd.addSlider("End_background_frames", 1, 20, endBackgroundFrames);
		gd.addSlider("Magnification", 5, 15, magnification);
		gd.addSlider("Smoothing", 0.25, 0.5, smoothing);
		gd.addCheckbox("Centre_each_slice", centreEachSlice);
		gd.addNumericField("CoM_cut_off", comCutOff, -2);
		String[] methods = ImageProcessor.getInterpolationMethods();
		gd.addChoice("Interpolation", methods, methods[interpolationMethod]);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		nmPerSlice = gd.getNextNumber();
		amplitudeFraction = gd.getNextNumber();
		startBackgroundFrames = (int) gd.getNextNumber();
		endBackgroundFrames = (int) gd.getNextNumber();
		magnification = (int) gd.getNextNumber();
		smoothing = gd.getNextNumber();
		centreEachSlice = gd.getNextBoolean();
		comCutOff = Maths.max(0, gd.getNextNumber());
		interpolationMethod = gd.getNextChoiceIndex();

		// Check arguments
		try
		{
			Parameters.isPositive("nm/slice", nmPerSlice);
			Parameters.isAbove("Amplitude fraction", amplitudeFraction, 0.01);
			Parameters.isBelow("Amplitude fraction", amplitudeFraction, 0.9);
			Parameters.isPositive("Start background frames", startBackgroundFrames);
			Parameters.isPositive("End background frames", endBackgroundFrames);
			Parameters.isAbove("Total background frames", startBackgroundFrames + endBackgroundFrames, 1);
			Parameters.isAbove("Magnification", magnification, 1);
			Parameters.isAbove("Smoothing", smoothing, 0);
			Parameters.isBelow("Smoothing", smoothing, 1);
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
	 */
	public void run(ImageProcessor ip)
	{
		if (mode == 0)
			runUsingProjections();
		else
			runUsingFitting();

		if (threadPool != null)
		{
			threadPool.shutdownNow();
			threadPool = null;
		}
	}

	private void runUsingFitting()
	{
		if (!showFittingDialog())
			return;
		if (!loadConfiguration())
			return;

		BasePoint[] spots = getSpots(0);
		if (spots.length == 0)
		{
			IJ.error(TITLE, "No spots without neighbours within " + (boxRadius * 2) + "px");
			return;
		}

		ImageStack stack = getImageStack();
		final int width = imp.getWidth();
		final int height = imp.getHeight();
		final int currentSlice = imp.getSlice();

		// Adjust settings for a single maxima
		config.setIncludeNeighbours(false);

		ArrayList<double[]> centres = new ArrayList<double[]>(spots.length);
		int iterations = 1;
		LoessInterpolator loess = new LoessInterpolator(smoothing, iterations);

		// TODO - The fitting routine may not produce many points. In this instance the LOESS interpolator
		// fails to smooth the data very well. A higher bandwidth helps this but perhaps 
		// try a different smoothing method.

		// For each spot
		Utils.log(TITLE + ": " + imp.getTitle());
		Utils.log("Finding spot locations...");
		Utils.log("  %d spot%s without neighbours within %dpx", spots.length, ((spots.length == 1) ? "" : "s"),
				(boxRadius * 2));
		StoredDataStatistics averageSd = new StoredDataStatistics();
		StoredDataStatistics averageA = new StoredDataStatistics();
		Statistics averageRange = new Statistics();
		MemoryPeakResults allResults = new MemoryPeakResults();
		allResults.setCalibration(fitConfig.getCalibration());
		allResults.setPSF(fitConfig.getPSF());
		allResults.setName(TITLE);
		allResults.setBounds(new Rectangle(0, 0, width, height));
		MemoryPeakResults.addResults(allResults);
		for (int n = 1; n <= spots.length; n++)
		{
			BasePoint spot = spots[n - 1];
			final int x = (int) spot.getX();
			final int y = (int) spot.getY();

			MemoryPeakResults results = fitSpot(stack, width, height, x, y);
			allResults.add(results);

			if (results.size() < 5)
			{
				Utils.log("  Spot %d: Not enough fit results %d", n, results.size());
				continue;
			}

			// Get the results for the spot centre and width
			final double[] z = new double[results.size()];
			final double[] xCoord = new double[z.length];
			final double[] yCoord = new double[z.length];
			final double[] sd;
			final double[] a;
			final Counter counter = new Counter();

			// We have fit the results so they will be in the preferred units 
			results.forEach(new PeakResultProcedure()
			{
				public void execute(PeakResult peak)
				{
					int i = counter.getAndIncrement();
					z[i] = peak.getFrame();
					xCoord[i] = peak.getXPosition() - x;
					yCoord[i] = peak.getYPosition() - y;
				}
			});

			WidthResultProcedure wp = new WidthResultProcedure(results, DistanceUnit.PIXEL);
			wp.getW();
			sd = SimpleArrayUtils.toDouble(wp.wx);

			HeightResultProcedure hp = new HeightResultProcedure(results, IntensityUnit.COUNT);
			hp.getH();
			a = SimpleArrayUtils.toDouble(hp.h);

			// Smooth the amplitude plot
			double[] smoothA = loess.smooth(z, a);

			// Find the maximum amplitude
			int maximumIndex = findMaximumIndex(smoothA);

			// Find the range at a fraction of the max. This is smoothed to find the X/Y centre
			int start = 0, stop = smoothA.length - 1;
			double limit = smoothA[maximumIndex] * amplitudeFraction;
			for (int j = 0; j < smoothA.length; j++)
			{
				if (smoothA[j] > limit)
				{
					start = j;
					break;
				}
			}
			for (int j = smoothA.length; j-- > 0;)
			{
				if (smoothA[j] > limit)
				{
					stop = j;
					break;
				}
			}
			averageRange.add(stop - start + 1);

			// Extract xy centre coords and smooth
			double[] smoothX = new double[stop - start + 1];
			double[] smoothY = new double[smoothX.length];
			double[] smoothSd = new double[smoothX.length];
			double[] newZ = new double[smoothX.length];
			for (int j = start, k = 0; j <= stop; j++, k++)
			{
				smoothX[k] = xCoord[j];
				smoothY[k] = yCoord[j];
				smoothSd[k] = sd[j];
				newZ[k] = z[j];
			}
			smoothX = loess.smooth(newZ, smoothX);
			smoothY = loess.smooth(newZ, smoothY);
			smoothSd = loess.smooth(newZ, smoothSd);

			// Since the amplitude is not very consistent move from this peak to the 
			// lowest width which is the in-focus spot.
			maximumIndex = findMinimumIndex(smoothSd, maximumIndex - start);

			// Find the centre at the amplitude peak
			double cx = smoothX[maximumIndex] + x;
			double cy = smoothY[maximumIndex] + y;
			int cz = (int) newZ[maximumIndex];
			double csd = smoothSd[maximumIndex];
			double ca = smoothA[maximumIndex + start];

			// The average should weight the SD using the signal for each spot
			averageSd.add(smoothSd[maximumIndex]);
			averageA.add(ca);

			if (ignoreSpot(n, z, a, smoothA, xCoord, yCoord, sd, newZ, smoothX, smoothY, smoothSd, cx, cy, cz, csd))
			{
				Utils.log("  Spot %d was ignored", n);
				continue;
			}

			// Store result - it may have been moved interactively
			maximumIndex += this.slice - cz;
			cz = (int) newZ[maximumIndex];
			csd = smoothSd[maximumIndex];
			ca = smoothA[maximumIndex + start];
			Utils.log("  Spot %d => x=%.2f, y=%.2f, z=%d, sd=%.2f, A=%.2f\n", n, cx, cy, cz, csd, ca);
			centres.add(new double[] { cx, cy, cz, csd, n });
		}

		if (interactiveMode)
		{
			imp.setSlice(currentSlice);
			imp.setOverlay(null);

			// Hide the amplitude and spot plots
			Utils.hide(TITLE_AMPLITUDE);
			Utils.hide(TITLE_PSF_PARAMETERS);
		}

		if (centres.isEmpty())
		{
			String msg = "No suitable spots could be identified";
			Utils.log(msg);
			IJ.error(TITLE, msg);
			return;
		}

		// Find the limits of the z-centre
		int minz = (int) centres.get(0)[2];
		int maxz = minz;
		for (double[] centre : centres)
		{
			if (minz > centre[2])
				minz = (int) centre[2];
			else if (maxz < centre[2])
				maxz = (int) centre[2];
		}

		IJ.showStatus("Creating PSF image");

		// Create a stack that can hold all the data.
		ImageStack psf = createStack(stack, minz, maxz, magnification);

		// For each spot
		Statistics stats = new Statistics();
		boolean ok = true;
		for (int i = 0; ok && i < centres.size(); i++)
		{
			double progress = (double) i / centres.size();
			final double increment = 1.0 / (stack.getSize() * centres.size());
			IJ.showProgress(progress);
			double[] centre = centres.get(i);

			// Extract the spot
			float[][] spot = new float[stack.getSize()][];
			Rectangle regionBounds = null;
			for (int slice = 1; slice <= stack.getSize(); slice++)
			{
				ImageExtractor ie = new ImageExtractor((float[]) stack.getPixels(slice), width, height);
				if (regionBounds == null)
					regionBounds = ie.getBoxRegionBounds((int) centre[0], (int) centre[1], boxRadius);
				spot[slice - 1] = ie.crop(regionBounds);
			}

			int n = (int) centre[4];
			final float b = getBackground(n, spot);
			if (!subtractBackgroundAndWindow(spot, b, regionBounds.width, regionBounds.height, centre, loess))
			{
				Utils.log("  Spot %d was ignored", n);
				continue;
			}

			stats.add(b);

			// Adjust the centre using the crop
			centre[0] -= regionBounds.x;
			centre[1] -= regionBounds.y;

			// This takes a long time so this should track progress
			ok = addToPSF(maxz, magnification, psf, centre, spot, regionBounds, progress, increment, centreEachSlice);
		}

		if (interactiveMode)
		{
			Utils.hide(TITLE_INTENSITY);
		}

		IJ.showProgress(1);

		if (!ok || stats.getN() == 0)
			return;

		final double avSd = getAverage(averageSd, averageA, 2);
		Utils.log("  Average background = %.2f, Av. SD = %s px", stats.getMean(), Utils.rounded(avSd, 4));

		normalise(psf, maxz, avSd * magnification, false);
		IJ.showProgress(1);

		psfImp = Utils.display("PSF", psf);
		psfImp.setSlice(maxz);
		psfImp.resetDisplayRange();
		psfImp.updateAndDraw();

		double[][] fitCom = new double[2][psf.getSize()];
		Arrays.fill(fitCom[0], Double.NaN);
		Arrays.fill(fitCom[1], Double.NaN);
		double fittedSd = fitPSF(psf, loess, maxz, averageRange.getMean(), fitCom);

		// Compute the drift in the PSF:
		// - Use fitted centre if available; otherwise find CoM for each frame
		// - express relative to the average centre

		double[][] com = calculateCentreOfMass(psf, fitCom, nmPerPixel / magnification);
		double[] slice = SimpleArrayUtils.newArray(psf.getSize(), 1, 1.0);
		String title = TITLE + " CoM Drift";
		Plot2 plot = new Plot2(title, "Slice", "Drift (nm)");
		plot.addLabel(0, 0, "Red = X; Blue = Y");
		//double[] limitsX = Maths.limits(com[0]);
		//double[] limitsY = Maths.limits(com[1]);
		double[] limitsX = getLimits(com[0]);
		double[] limitsY = getLimits(com[1]);
		plot.setLimits(1, psf.getSize(), Math.min(limitsX[0], limitsY[0]), Math.max(limitsX[1], limitsY[1]));
		plot.setColor(Color.red);
		plot.addPoints(slice, com[0], Plot.DOT);
		plot.addPoints(slice, loess.smooth(slice, com[0]), Plot.LINE);
		plot.setColor(Color.blue);
		plot.addPoints(slice, com[1], Plot.DOT);
		plot.addPoints(slice, loess.smooth(slice, com[1]), Plot.LINE);
		Utils.display(title, plot);

		// TODO - Redraw the PSF with drift correction applied. 
		// This means that the final image should have no drift.
		// This is relevant when combining PSF images. It doesn't matter too much for simulations 
		// unless the drift is large.

		// Add Image properties containing the PSF details
		final double fwhm = getFWHM(psf, maxz);
		psfImp.setProperty("Info", ImagePSFHelper.toString(
				ImagePSFHelper.create(maxz, nmPerPixel / magnification, nmPerSlice, stats.getN(), fwhm, createNote())));

		Utils.log("%s : z-centre = %d, nm/Pixel = %s, nm/Slice = %s, %d images, PSF SD = %s nm, FWHM = %s px\n",
				psfImp.getTitle(), maxz, Utils.rounded(nmPerPixel / magnification, 3), Utils.rounded(nmPerSlice, 3),
				stats.getN(), Utils.rounded(fittedSd * nmPerPixel, 4), Utils.rounded(fwhm));

		createInteractivePlots(psf, maxz, nmPerPixel / magnification, fittedSd * nmPerPixel);

		IJ.showStatus("");
	}

	/**
	 * Get the limits of the array ignoring outliers more than 1.5x the inter quartile range
	 * 
	 * @param data
	 * @return
	 */
	private double[] getLimits(double[] data)
	{
		double[] limits = Maths.limits(data);
		DescriptiveStatistics stats = new DescriptiveStatistics(data);
		double lower = stats.getPercentile(25);
		double upper = stats.getPercentile(75);
		double iqr = (upper - lower) * 2;
		limits[0] = FastMath.max(lower - iqr, limits[0]);
		limits[1] = FastMath.min(upper + iqr, limits[1]);
		return limits;
	}

	private double getAverage(StoredDataStatistics averageSd, StoredDataStatistics averageA, int averageMethod)
	{
		if (averageMethod == 0)
			return averageSd.getMean();
		double[] sd = averageSd.getValues();
		double[] w = averageA.getValues();
		double sum = 0, sumW = 0;

		if (averageMethod == 1)
		{
			// Weighted average using Amplitude
		}
		else if (averageMethod == 2)
		{
			// Weighted average using signal
			for (int i = 0; i < sd.length; i++)
			{
				w[i] *= sd[i] * sd[i];
			}
		}

		for (int i = 0; i < sd.length; i++)
		{
			sum += sd[i] * w[i];
			sumW += w[i];
		}

		return sum / sumW;
	}

	private MemoryPeakResults fitSpot(ImageStack stack, final int width, final int height, final int x, final int y)
	{
		Rectangle regionBounds = null;

		// Create a fit engine
		MemoryPeakResults results = new MemoryPeakResults();
		results.setCalibration(fitConfig.getCalibration());
		results.setPSF(fitConfig.getPSF());
		results.setSortAfterEnd(true);
		results.begin();
		int threadCount = Prefs.getThreads();
		FitEngine engine = new FitEngine(config, SynchronizedPeakResults.create(results, threadCount), threadCount,
				FitQueue.BLOCKING);

		List<ParameterisedFitJob> jobItems = new ArrayList<ParameterisedFitJob>(stack.getSize());

		for (int slice = 1; slice <= stack.getSize(); slice++)
		{
			// Extract the region from each frame
			ImageExtractor ie = new ImageExtractor((float[]) stack.getPixels(slice), width, height);
			if (regionBounds == null)
				regionBounds = ie.getBoxRegionBounds(x, y, boxRadius);
			float[] region = ie.crop(regionBounds);

			// Fit only a spot in the centre
			FitParameters params = new FitParameters();
			params.maxIndices = new int[] { boxRadius * regionBounds.width + boxRadius };
			ParameterisedFitJob job = new ParameterisedFitJob(slice, params, slice, region, regionBounds);
			jobItems.add(job);
			engine.run(job);
		}

		engine.end(false);
		results.end();
		return results;
	}

	private int findMaximumIndex(double[] data)
	{
		double max = data[0];
		int pos = 0;
		for (int j = 0; j < data.length; j++)
		{
			if (max < data[j])
			{
				max = data[j];
				pos = j;
			}
		}
		return pos;
	}

	private int findMinimumIndex(double[] data, int initialGuess)
	{
		double min = data[initialGuess];
		int pos = initialGuess;
		// Move only downhill from the initial guess.
		for (int j = initialGuess + 1; j < data.length; j++)
		{
			if (min > data[j])
			{
				min = data[j];
				pos = j;
			}
			else
			{
				break;
			}
		}
		for (int j = initialGuess; j-- > 0;)
		{
			if (min > data[j])
			{
				min = data[j];
				pos = j;
			}
			else
			{
				break;
			}
		}
		return pos;
	}

	private boolean ignoreSpot(int n, final double[] z, final double[] a, final double[] smoothA, final double[] xCoord,
			final double[] yCoord, final double[] sd, final double[] newZ, final double[] smoothX,
			final double[] smoothY, double[] smoothSd, final double cx, final double cy, final int cz, double csd)
	{
		this.slice = cz;
		// Allow an interactive mode that shows the plots and allows the user to Yes/No
		// the addition of the data
		if (interactiveMode)
		{
			zCentre = cz;
			psfWidth = csd * nmPerPixel;

			// Store the data for replotting
			this.z = z;
			this.a = a;
			this.smoothAz = z;
			this.smoothA = smoothA;
			this.xCoord = xCoord;
			this.yCoord = yCoord;
			this.sd = sd;
			this.newZ = newZ;
			this.smoothX = smoothX;
			this.smoothY = smoothY;
			this.smoothSd = smoothSd;

			showPlots(z, a, z, smoothA, xCoord, yCoord, sd, newZ, smoothX, smoothY, smoothSd, cz);

			// Draw the region on the input image as an overlay
			imp.setSlice(cz);
			imp.setOverlay(
					new Roi((int) (cx - boxRadius), (int) (cy - boxRadius), 2 * boxRadius + 1, 2 * boxRadius + 1),
					Color.GREEN, 1, null);

			// Ask if the spot should be included
			GenericDialog gd = new GenericDialog(TITLE);
			gd.enableYesNoCancel();
			gd.hideCancelButton();
			gd.addMessage(String.format(
					"Add spot %d to the PSF?\n \nEstimated centre using min PSF width:\n \nx = %.2f\ny = %.2f\nz = %d\nsd = %.2f\n",
					n, cx, cy, cz, csd));
			gd.addSlider("Slice", z[0], z[z.length - 1], slice);
			if (yesNoPosition != null)
			{
				gd.centerDialog(false);
				gd.setLocation(yesNoPosition);
			}
			gd.addDialogListener(new SimpleInteractivePlotListener());
			gd.showDialog();

			yesNoPosition = gd.getLocation();
			return !gd.wasOKed();
		}
		return false;
	}

	private class SimpleInteractivePlotListener implements DialogListener
	{
		public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
		{
			slice = (int) gd.getNextNumber();
			drawPlots(false);
			return true;
		}
	}

	private void showPlots(final double[] z, final double[] a, final double[] smoothAz, final double[] smoothA,
			final double[] xCoord, final double[] yCoord, final double[] sd, final double[] newZ,
			final double[] smoothX, final double[] smoothY, double[] smoothSd, final int cz)
	{
		PlotWindow amplitudeWindow = null;

		// Draw a plot of the amplitude
		if (a != null)
		{
			Plot2 plot = new Plot2(TITLE_AMPLITUDE, "z", "Amplitude", smoothAz, smoothA);
			double[] limits2 = Maths.limits(Maths.limits(a), smoothA);
			plot.setLimits(z[0], z[z.length - 1], limits2[0], limits2[1]);
			plot.addPoints(z, a, Plot2.CIRCLE);

			// Add a line for the z-centre
			plot.setColor(Color.GREEN);
			plot.addPoints(new double[] { cz, cz }, limits2, Plot2.LINE);
			plot.setColor(Color.BLACK);

			double amplitude = Double.NaN;
			for (int i = 0; i < smoothAz.length; i++)
			{
				if (smoothAz[i] == cz)
				{
					amplitude = smoothA[i];
					break;
				}
			}
			double maxAmplitude = Double.NaN;
			for (int i = 0; i < smoothAz.length; i++)
			{
				if (smoothAz[i] == zCentre)
				{
					maxAmplitude = smoothA[i];
					break;
				}
			}
			plot.addLabel(0, 0, String.format("Amplitude = %s (%sx). z = %s nm", Utils.rounded(amplitude),
					Utils.rounded(amplitude / maxAmplitude), Utils.rounded((slice - zCentre) * nmPerSlice)));

			amplitudeWindow = Utils.display(TITLE_AMPLITUDE, plot);
		}

		// Show plot of width, X centre, Y centre
		if (xCoord != null)
		{
			Plot2 plot = new Plot2(TITLE_PSF_PARAMETERS, "z", "px", newZ, smoothSd);
			// Get the limits
			double[] sd2 = invert(sd);
			double[] limits = Maths.limits(Maths.limits(Maths.limits(Maths.limits(xCoord), yCoord), sd), sd2);
			plot.setLimits(z[0], z[z.length - 1], limits[0], limits[1]);
			plot.addPoints(newZ, invert(smoothSd), Plot2.LINE);
			plot.addPoints(z, sd, Plot2.DOT);
			plot.addPoints(z, sd2, Plot2.DOT);
			plot.setColor(Color.BLUE);
			plot.addPoints(z, xCoord, Plot2.DOT);
			plot.addPoints(newZ, smoothX, Plot2.LINE);
			plot.setColor(Color.RED);
			plot.addPoints(z, yCoord, Plot2.DOT);
			plot.addPoints(newZ, smoothY, Plot2.LINE);

			// Add a line for the z-centre
			plot.setColor(Color.GREEN);
			plot.addPoints(new double[] { cz, cz }, limits, Plot2.LINE);
			plot.setColor(Color.BLACK);

			double width = Double.NaN;
			for (int i = 0; i < smoothSd.length; i++)
			{
				if (newZ[i] == cz)
				{
					width = smoothSd[i];
					break;
				}
			}
			plot.addLabel(0, 0, String.format("Width = %s nm (%sx). z = %s nm", Utils.rounded(width * nmPerPixel),
					Utils.rounded(width * nmPerPixel / psfWidth), Utils.rounded((slice - zCentre) * nmPerSlice)));

			// Check if the window will need to be aligned
			boolean alignWindows = (WindowManager.getFrame(TITLE_PSF_PARAMETERS) == null);

			PlotWindow psfWindow = Utils.display(TITLE_PSF_PARAMETERS, plot);

			if (alignWindows && psfWindow != null && amplitudeWindow != null)
			{
				// Put the two plots tiled together so both are visible
				Point l = psfWindow.getLocation();
				l.x = amplitudeWindow.getLocation().x;
				l.y = amplitudeWindow.getLocation().y + amplitudeWindow.getHeight();
				psfWindow.setLocation(l);
			}
		}
	}

	private double[] invert(final double[] data)
	{
		double[] data2 = new double[data.length];
		for (int i = 0; i < data.length; i++)
			data2[i] = -data[i];
		return data2;
	}

	private ImageStack createStack(ImageStack stack, int minz, int maxz, final int magnification)
	{
		// Pad box radius with an extra pixel border to allow offset insertion
		final int w = ((2 * boxRadius + 1) + 2) * magnification;
		final int d = maxz - minz + stack.getSize();
		ImageStack psf = new ImageStack(w, w, d);
		for (int i = 1; i <= d; i++)
			psf.setPixels(new float[w * w], i);
		return psf;
	}

	private float getBackground(int n, float[][] spot)
	{
		// Get the average value of the first and last n frames
		Statistics first = new Statistics();
		Statistics last = new Statistics();
		for (int i = 0; i < startBackgroundFrames; i++)
		{
			first.add(spot[i]);
		}
		for (int i = 0, j = spot.length - 1; i < endBackgroundFrames; i++, j--)
		{
			last.add(spot[j]);
		}
		float av = (float) ((first.getSum() + last.getSum()) / (first.getN() + last.getN()));
		Utils.log("  Spot %d Background: First %d = %.2f, Last %d = %.2f, av = %.2f", n, startBackgroundFrames,
				first.getMean(), endBackgroundFrames, last.getMean(), av);
		return av;
	}

	@SuppressWarnings("unused")
	private float getBackground(final double fraction, DoubleData all)
	{
		double[] allValues = all.values();
		Arrays.sort(allValues);
		int fractionIndex = (int) (allValues.length * fraction);
		double sum = 0;
		for (int i = 0; i <= fractionIndex; i++)
		{
			sum += allValues[i];
		}
		final float min = (float) (sum / (fractionIndex + 1));
		return min;
	}

	private boolean[] dmap = null;
	private int lastWidth = 0;
	private int lastHeight = 0;
	private int minx, maxx, miny, maxy;

	/**
	 * Subtract the background from the spot, compute the intensity within half the box region distance from the centre
	 * and smooth the intensity profile. In interactive mode the user must choose to accept the profile or reject.
	 * If accepted the smoothed profile is user to normalise the image and then the image is rolled off to zero
	 * using a Tukey window function.
	 * 
	 * @param spot
	 * @param background
	 *            The minimum level, all below this is background and set to zero
	 * @param spotWidth
	 * @param spotHeight
	 * @param n
	 *            The spot number
	 * @param loess
	 *            The smoothing interpolator
	 * @return True if accepted
	 */
	private boolean subtractBackgroundAndWindow(float[][] spot, final float background, final int spotWidth,
			final int spotHeight, double[] centre, LoessInterpolator loess)
	{
		//ImageWindow imageWindow = new ImageWindow();
		for (int i = 0; i < spot.length; i++)
		{
			for (int j = 0; j < spot[i].length; j++)
				spot[i][j] = FastMath.max(spot[i][j] - background, 0);
		}

		// Create a distance map from the centre
		if (lastWidth != spotWidth || lastHeight != spotHeight)
		{
			final double cx = spotWidth * 0.5;
			final double cy = spotHeight * 0.5;
			minx = FastMath.max(0, (int) (cx - boxRadius * 0.5));
			maxx = FastMath.min(spotWidth, (int) Math.ceil(cx + boxRadius * 0.5));
			miny = FastMath.max(0, (int) (cy - boxRadius * 0.5));
			maxy = FastMath.min(spotHeight, (int) Math.ceil(cy + boxRadius * 0.5));

			// Precompute square distances
			double[] dx2 = new double[maxx - minx + 1];
			for (int x = minx, i = 0; x < maxx; x++, i++)
			{
				// Use pixel centres with 0.5 offset
				final double dx = x + 0.5 - cx;
				dx2[i] = dx * dx;
			}
			dmap = new boolean[dx2.length * (maxy - miny + 1)];
			final double d2 = boxRadius * boxRadius / 4;
			for (int y = miny, j = 0; y < maxy; y++)
			{
				final double dy = (y + 0.5 - cy);
				final double dy2 = dy * dy;
				final double limit = d2 - dy2;
				for (int x = minx, i = 0; x < maxx; x++, i++, j++)
				{
					dmap[j] = (dx2[i] < limit);
				}
			}
			lastWidth = spotWidth;
			lastHeight = spotHeight;
		}

		// Calculate the intensity profile within half the box radius from the centre
		double[] xValues = new double[spot.length];
		double[] yValues = new double[spot.length];
		for (int i = 0; i < spot.length; i++)
		{
			xValues[i] = i + 1;
			double sum = 0;
			for (int y = miny, j = 0; y < maxy; y++)
			{
				int index = y * spotWidth + minx;
				for (int x = minx; x < maxx; x++, index++, j++)
					if (dmap[j])
						sum += spot[i][index];
			}
			yValues[i] = sum;
		}

		double[] newY = loess.smooth(xValues, yValues);
		// It can happen that the LOESS creates values below zero (e.g. when the curve
		// falls towards zero at the ends)
		for (int i = 0; i < newY.length; i++)
			if (newY[i] < 0)
				newY[i] = yValues[i];

		if (interactiveMode)
		{
			Utils.hide(TITLE_AMPLITUDE);
			Utils.hide(TITLE_PSF_PARAMETERS);

			final int n = (int) centre[4];

			String title = TITLE_INTENSITY;
			Plot plot = new Plot(title, "Slice", "Sum", xValues, yValues);
			plot.setColor(Color.red);
			plot.addPoints(xValues, newY, Plot.LINE);
			plot.setColor(Color.green);
			double[] limits = Maths.limits(yValues);
			plot.drawLine(centre[2], limits[0], centre[2], limits[1]);
			plot.setColor(Color.black);
			plot.addLabel(0, 0, "Spot " + n);
			Utils.display(title, plot);

			GenericDialog gd = new GenericDialog(TITLE);
			gd.enableYesNoCancel();
			gd.hideCancelButton();
			gd.addMessage(String.format(
					"Add spot %d to the PSF?\n(The intensity profile is the sum within half the box region)", n));
			if (yesNoPosition != null)
			{
				gd.centerDialog(false);
				gd.setLocation(yesNoPosition);
			}
			gd.showDialog();

			yesNoPosition = gd.getLocation();
			if (!gd.wasOKed())
				return false;
		}

		for (int i = 0; i < spot.length; i++)
		{
			// Normalise
			final float scale = (float) (newY[i] / yValues[i]);
			for (int j = 0; j < spot[i].length; j++)
				spot[i][j] *= scale;

			// Use a Tukey window to roll-off the image edges
			//spot[i] = imageWindow.applySeperable(spot[i], spotWidth, spotHeight, ImageWindow.WindowFunction.Tukey);
			spot[i] = ImageWindow.applyWindow(spot[i], spotWidth, spotHeight, ImageWindow.WindowFunction.TUKEY);
		}

		return true;
	}

	private boolean addToPSF(int maxz, final int magnification, ImageStack psf, final double[] centre,
			final float[][] spot, final Rectangle regionBounds, double progress, final double increment,
			final boolean centreEachSlice)
	{
		// Calculate insert point in enlargement
		final int z = (int) centre[2];
		int insertZ = maxz - z + 1;

		// Enlargement size
		final int w = regionBounds.width, h = regionBounds.height;
		final int dstWidth = w * magnification;
		final int dstHeight = h * magnification;

		// Multi-thread for speed
		if (threadPool == null)
			threadPool = Executors.newFixedThreadPool(Prefs.getThreads());

		List<Future<?>> futures = new LinkedList<Future<?>>();

		for (int i = 0; i < spot.length; i++)
		{
			//final int slice = i + 1;
			final ImageProcessor ip = psf.getProcessor(insertZ++);
			final float[] originalSpotData = spot[i];

			futures.add(threadPool.submit(new Runnable()
			{
				public void run()
				{
					if (Utils.isInterrupted())
						return;

					incrementProgress(increment);

					double insertX, insertY;

					// Enlarge
					FloatProcessor fp = new FloatProcessor(w, h, originalSpotData, null);
					fp.setInterpolationMethod(interpolationMethod);
					fp = (FloatProcessor) fp.resize(dstWidth, dstHeight);

					// In the case of Bicubic interpolation check for negative values
					if (interpolationMethod == ImageProcessor.BICUBIC)
					{
						float[] pixels = (float[]) fp.getPixels();
						for (int i = 0; i < pixels.length; i++)
							if (pixels[i] < 0)
								pixels[i] = 0;
					}

					// Do all CoM calculations here since we use an interpolation
					// when resizing and the CoM will move.
					if (centreEachSlice)
					{
						final double[] com = calculateCenterOfMass(fp);
						//System.out.printf("CoM %d : %f %f vs %f %f\n", slice, com[0], com[1],
						//		centre[0] * magnification, centre[1] * magnification);

						// Get the insert position by subtracting the centre-of-mass of the enlarged image from the 
						// image centre + allow for a border of 1 pixel * magnification
						insertX = magnification + dstWidth * 0.5 - com[0];
						insertY = magnification + dstHeight * 0.5 - com[1];
						//Utils.log("Insert point = %.2f,%.2f => %.2f,%.2f\n", dstWidth * 0.5 - cx, dstHeight * 0.5 - cy,
						//		insertX, insertY);
					}
					else
					{
						// Get the insert position from the stack centre using enlargement
						insertX = getInsert(centre[0], (int) centre[0], magnification);
						insertY = getInsert(centre[1], (int) centre[1], magnification);
						//Utils.log("Insert point = %.2f,%.2f => %.2f,%.2f\n", centre[0] - (int) centre[0], centre[1] - (int) centre[1], insertX, insertY);
					}

					// Copy the processor using a weighted image
					final int lowerX = (int) insertX;
					final int lowerY = (int) insertY;

					final double wx2 = insertX - lowerX;
					final double wx1 = 1 - wx2;
					final double wy2 = insertY - lowerY;
					final double wy1 = 1 - wy2;

					// Add to the combined PSF using the correct offset and the weighting
					copyBits(ip, fp, lowerX, lowerY, wx1 * wy1);
					copyBits(ip, fp, lowerX + 1, lowerY, wx2 * wy1);
					copyBits(ip, fp, lowerX, lowerY + 1, wx1 * wy2);
					copyBits(ip, fp, lowerX + 1, lowerY + 1, wx2 * wy2);

					//// Check CoM is correct. This is never perfect since the bilinear weighting 
					//// interpolates the data and shifts the CoM.
					//ImageProcessor ip2 = ip.createProcessor(ip.getWidth(), ip.getHeight());
					//copyBits(ip2, fp, lowerX, lowerY, wx1 * wy1);
					//copyBits(ip2, fp, lowerX + 1, lowerY, wx2 * wy1);
					//copyBits(ip2, fp, lowerX, lowerY + 1, wx1 * wy2);
					//copyBits(ip2, fp, lowerX + 1, lowerY + 1, wx2 * wy2);
					//
					//double[] com = getCoM((FloatProcessor) ip2);
					//System.out.printf("Inserted CoM %d : %f %f\n", slice, com[0], com[1]);
				}
			}));

			if (Utils.isInterrupted())
				break;
		}

		Utils.waitForCompletion(futures);

		return !Utils.isInterrupted();
	}

	private static double[] calculateCenterOfMass(FloatProcessor fp)
	{
		final int h = fp.getHeight();
		final int w = fp.getWidth();
		float[] data = (float[]) fp.getPixels();
		final double threshold = Maths.max(data) * comCutOff;
		double sx = 0, sy = 0, s = 0;
		for (int y = 0, i = 0; y < h; y++)
			for (int x = 0; x < w; x++, i++)
			{
				final float v = data[i];
				if (v >= threshold)
				{
					sx += x * v;
					sy += y * v;
					s += v;
				}
			}
		// Allow for centre of pixel to be at 0.5
		return new double[] { 0.5 + sx / s, 0.5 + sy / s };
	}

	/**
	 * Calculate the insertion position so that the spot is added at exactly the centre of the PSF
	 * 
	 * @param coord
	 *            The coordinate
	 * @param iCoord
	 *            The coordinate rounded down to an integer
	 * @param magnification
	 *            The magnification
	 * @return The insert position
	 */
	private final double getInsert(final double coord, final int iCoord, final int magnification)
	{
		// Note that a perfect alignment to the centre of a pixel would be 0.5,0.5.
		// Insert should align the image into the middle:
		// Offset in pixel       Insert
		// 0.0               =>  +0.5
		// 0.1               =>  +0.4
		// 0.2               =>  +0.3
		// 0.3               =>  +0.2
		// 0.4               =>  +0.1
		// 0.5               =>  +0.0
		// 0.6               =>  -0.1
		// 0.7               =>  -0.2
		// 0.8               =>  -0.3
		// 0.9               =>  -0.4
		// 1.0               =>  -0.5

		// Off set should range from 0 to 1
		final double offset = (coord - iCoord);
		// Insert point is in the opposite direction from the offset (range from -0.5 to 0.5)
		final double insert = -1 * (offset - 0.5);
		//return magnification + (int) Math.round(insert * magnification);
		return magnification + (insert * magnification);
	}

	private synchronized void incrementProgress(final double increment)
	{
		progress += increment;
		IJ.showProgress(progress);
	}

	private void copyBits(ImageProcessor ip, FloatProcessor fp, int lowerX, int lowerY, double weight)
	{
		if (weight > 0)
		{
			fp = (FloatProcessor) fp.duplicate();
			fp.multiply(weight);
			ip.copyBits(fp, lowerX, lowerY, Blitter.ADD);
		}
	}

	/**
	 * Normalise the PSF using a given denominator
	 * 
	 * @param psf
	 * @param n
	 *            The denominator
	 */
	public static void normaliseUsingSpots(ImageStack psf, int n)
	{
		if (psf == null || psf.getSize() == 0)
			return;
		if (!(psf.getPixels(1) instanceof float[]))
			return;
		for (int i = 0; i < psf.getSize(); i++)
		{
			float[] data = (float[]) psf.getPixels(i + 1);
			for (int j = 0; j < data.length; j++)
				data[j] /= n;
		}
	}

	/**
	 * Normalise the PSF so the sum of the specified frame foreground pixels is 1.
	 * <p>
	 * Assumes the PSF can be approximated by a Gaussian in the central frame. All pixels within 3 sigma of the centre
	 * are foreground pixels.
	 * 
	 * @param psf
	 * @param n
	 *            The frame number
	 * @param sigma
	 *            the Gaussian standard deviation (in pixels)
	 * @param subtractBackground
	 *            Normalise so everything below the background is zero
	 */
	public static void normalise(ImageStack psf, int n, double sigma, boolean subtractBackground)
	{
		if (psf == null || psf.getSize() == 0)
			return;
		if (!(psf.getPixels(1) instanceof float[]))
			return;
		final double cx = psf.getWidth() * 0.5;

		// Get the sum of the foreground pixels
		float[] data = (float[]) psf.getPixels(n);
		double foregroundSum = 0;
		int foregroundN = 0;
		final int min = FastMath.max(0, (int) (cx - 3 * sigma));
		final int max = FastMath.min(psf.getWidth() - 1, (int) Math.ceil(cx + 3 * sigma));

		// Precompute square distances within 3 sigma of the centre
		final double r2 = 3 * sigma * 3 * sigma;
		double[] d2 = new double[max - min + 1];
		for (int x = min, i = 0; x <= max; x++, i++)
			// Use pixel centres with 0.5 offset
			d2[i] = (x + 0.5 - cx) * (x + 0.5 - cx);

		for (int y = min, i = 0; y <= max; y++, i++)
		{
			int index = y * psf.getWidth() + min;
			final double limit = r2 - d2[i];
			for (int x = min, j = 0; x <= max; x++, index++, j++)
			{
				// Check if the pixel is within 3 sigma of the centre
				if (d2[j] < limit)
				{
					foregroundSum += data[index];
					foregroundN++;
				}
			}
		}

		if (subtractBackground)
		{
			// Normalise so everything below the background is zero

			// Get the average background
			final double backgroundSum = Maths.sum(data) - foregroundSum;
			final double background = backgroundSum / (data.length - foregroundN);

			// Subtract the background from the foreground sum
			final double newForegroundSum = foregroundSum - background * foregroundN;

			for (int i = 0; i < psf.getSize(); i++)
			{
				data = (float[]) psf.getPixels(i + 1);
				for (int j = 0; j < data.length; j++)
				{
					data[j] = (float) (Math.max(0, data[j] - background) / newForegroundSum);
				}
			}
		}
		else
		{
			for (int i = 0; i < psf.getSize(); i++)
			{
				data = (float[]) psf.getPixels(i + 1);
				for (int j = 0; j < data.length; j++)
				{
					// Normalise so the foreground is 1
					data[j] = (float) (data[j] / foregroundSum);
				}
			}
		}
	}

	/**
	 * Normalise the PSF so the sum of specified frame is 1.
	 * 
	 * @param psf
	 * @param n
	 *            The frame number
	 */
	public static void normalise(ImageStack psf, int n)
	{
		if (psf == null || psf.getSize() == 0)
			return;
		if (!(psf.getPixels(1) instanceof float[]))
			return;
		double sum = Maths.sum((float[]) psf.getPixels(n));
		for (int i = 0; i < psf.getSize(); i++)
		{
			float[] data = (float[]) psf.getPixels(i + 1);
			for (int j = 0; j < data.length; j++)
				data[j] /= sum;
		}
	}

	/**
	 * Calculate the centre of mass and express it relative to the average centre
	 * 
	 * @param psf
	 * @param fitCom
	 * @param nmPerPixel
	 * @return The centre of mass
	 */
	private double[][] calculateCentreOfMass(ImageStack psf, double[][] fitCom, double nmPerPixel)
	{
		final int size = psf.getSize();
		double[][] com = new double[2][size];
		final double offset = psf.getWidth() / 2.0;
		for (int i = 0; i < size; i++)
		{
			final double[] com2 = calculateCenterOfMass((FloatProcessor) psf.getProcessor(i + 1));
			com[0][i] = com2[0] - offset;
			com[1][i] = com2[1] - offset;
			//if (!Double.isNaN(fitCom[0][i]))
			//{
			//	// Interlacing the fit centre of mass is not consistent. There appears to be a large discrepancy
			//	// between the pixel centre-of-mass and the fit CoM. A small test shows correlation of
			//	// 0.11 and 0.066. Spearman's rank is 0.16. Basically it messes the data and effects smoothing.
			//	//System.out.printf("CoM = [ %f , %f ] == [ %f , %f ]\n", comX, comY, fitCom[0][i], fitCom[1][i]);
			//	//com[0][i] = fitCom[0][i];
			//	//com[1][i] = fitCom[1][i];
			//}
		}

		// Smooth the curve ...
		//		LoessInterpolator loess = new LoessInterpolator(smoothing, 1);
		//		double[] slice = SimpleArrayUtils.newArray(psf.getSize(), 1, 1.0);
		//		com[0] = loess.smooth(slice, com[0]);
		//		com[1] = loess.smooth(slice, com[1]);

		// Express relative to the average centre
		final double avX = new Statistics(com[0]).getMean();
		final double avY = new Statistics(com[1]).getMean();
		for (int i = 0; i < size; i++)
		{
			com[0][i] = (com[0][i] - avX) * nmPerPixel;
			com[1][i] = (com[1][i] - avY) * nmPerPixel;
		}

		return com;
	}

	private boolean loadConfiguration()
	{
		Configuration c = new Configuration();
		// TODO: We could have a different fit configuration just for the PSF Creator.
		// This would allow it to be saved and not effect PeakFit settings.
		boolean save = true;
		if (!c.showDialog(save))
		{
			IJ.error(TITLE, "No fit configuration loaded");
			return false;
		}

		config = c.getFitEngineConfiguration();
		config.configureOutputUnits();
		config.setResidualsThreshold(1);
		fitConfig = config.getFitConfiguration();
		nmPerPixel = fitConfig.getCalibrationWriter().getNmPerPixel();
		if (radius < 5 * FastMath.max(fitConfig.getInitialXSD(), fitConfig.getInitialYSD()))
		{
			radius = 5 * FastMath.max(fitConfig.getInitialXSD(), fitConfig.getInitialYSD());
			Utils.log("Radius is less than 5 * PSF standard deviation, increasing to %s", Utils.rounded(radius));
		}
		boxRadius = (int) Math.ceil(radius);
		return true;
	}

	/**
	 * @return Extract all the ROI points that are not within twice the box radius of any other spot
	 */
	private BasePoint[] getSpots(float offset)
	{
		float z = imp.getStackSize() / 2;
		Roi roi = imp.getRoi();
		if (roi != null && roi.getType() == Roi.POINT)
		{
			FloatPolygon p = roi.getFloatPolygon();
			int n = p.npoints;

			if (offset != 0)
			{
				// Check if already float coordinates
				if (!SimpleArrayUtils.isInteger(p.xpoints) || !SimpleArrayUtils.isInteger(p.ypoints))
					offset = 0;
			}

			BasePoint[] roiPoints = new BasePoint[n];
			for (int i = 0; i < n; i++)
			{
				roiPoints[i] = new BasePoint(p.xpoints[i] + offset, p.ypoints[i] + offset, z);
			}

			// All vs all distance matrix
			double[][] d = new double[n][n];
			for (int i = 0; i < n; i++)
				for (int j = i + 1; j < n; j++)
					d[i][j] = d[j][i] = roiPoints[i].distanceXY2(roiPoints[j]);

			// Spots must be twice as far apart to have no overlap of the extracted box region
			double d2 = boxRadius * boxRadius * 4;
			int ok = 0;
			for (int i = 0; i < n; i++)
			{
				if (noNeighbours(d, n, i, d2))
					roiPoints[ok++] = roiPoints[i];
			}

			return Arrays.copyOf(roiPoints, ok);
		}
		return new BasePoint[0];
	}

	/**
	 * Check spots box region do not overlap.
	 *
	 * @param centres
	 *            the centres
	 * @return true, if successful
	 */
	private boolean checkSpots(BasePoint[] centres)
	{
		int n = centres.length;

		// Spots must be twice as far apart to have no overlap of the extracted box region
		double d2 = boxRadius * boxRadius * 4;

		for (int i = 0; i < n; i++)
			for (int j = i + 1; j < n; j++)
				if (centres[i].distanceXY2(centres[j]) < d2)
					return false;

		return true;
	}

	/**
	 * Check the spot is not within the given squared distance from any other spot
	 * 
	 * @param d
	 *            The distance matrix
	 * @param n
	 *            The number of spots
	 * @param i
	 *            The spot
	 * @param d2
	 *            The squared distance
	 * @return True if there are no neighbours
	 */
	private boolean noNeighbours(final double[][] d, final int n, final int i, final double d2)
	{
		for (int j = 0; j < n; j++)
		{
			if (i != j && d[i][j] < d2)
			{
				return false;
			}
		}
		return true;
	}

	/**
	 * @return The input image as a 32-bit (float) image stack
	 */
	private ImageStack getImageStack()
	{
		final int width = imp.getWidth();
		final int height = imp.getHeight();
		ImageStack stack = imp.getImageStack();
		ImageStack newStack = new ImageStack(width, height, stack.getSize());
		for (int slice = 1; slice <= stack.getSize(); slice++)
		{
			newStack.setPixels(ImageConverter.getData(stack.getPixels(slice), width, height, null, null), slice);
		}
		return newStack;
	}

	/**
	 * Fit the new PSF image and show a graph of the amplitude/width
	 * 
	 * @param psfStack
	 * @param loess
	 * @param averageRange
	 * @param fitCom
	 * @return The width of the PSF in the z-centre
	 */
	private double fitPSF(ImageStack psfStack, LoessInterpolator loess, int cz, double averageRange,
			final double[][] fitCom)
	{
		IJ.showStatus("Fitting final PSF");

		// Note: Fitting the final PSF does not really work using MLE. This is because the noise model
		// is not appropriate for a normalised PSF. 
		if (fitConfig.getFitSolver() != FitSolver.LVM_LSE)
		{
			Utils.log("  " + FitProtosHelper.getName(fitConfig.getFitSolver()) +
					" is not appropriate for final PSF fitting.");
			Utils.log("  Switching to Least Square Estimation");
			fitConfig.setFitSolver(FitSolver.LVM_LSE);
			if (interactiveMode)
			{
				// This assumes the LVM does not need the calibration
				PeakFit.configureFitSolver(config, 0, 0, 0);
			}
		}

		// Update the box radius since this is used in the fitSpot method.
		boxRadius = psfStack.getWidth() / 2;
		final int x = boxRadius;
		final int y = boxRadius;
		FitConfiguration fitConfig = config.getFitConfiguration();
		final double shift = fitConfig.getCoordinateShiftFactor();

		// Scale the PSF
		PSF.Builder psf = fitConfig.getPSF().toBuilder();
		for (int i = 0; i < psf.getParametersCount(); i++)
		{
			PSFParameter param = psf.getParameters(i);
			if (param.getUnit() == PSFParameterUnit.DISTANCE)
			{
				PSFParameter.Builder b = param.toBuilder();
				b.setValue(b.getValue() * magnification);
				psf.setParameters(i, b);
			}
		}
		fitConfig.setPSF(psf.build());

		// Need to be updated after the widths have been set
		fitConfig.setCoordinateShiftFactor(shift);
		fitConfig.setBackgroundFitting(false);
		// Since the PSF will be normalised
		fitConfig.setMinPhotons(0);
		fitConfig.setBias(0);
		fitConfig.setGain(1);
		// No complex filtering so we get a fit. It should be easy to fit anyway.
		fitConfig.setPrecisionThreshold(0);
		fitConfig.setDirectFilter(null);
		//fitConfig.setDisableSimpleFilter(true);
		//fitConfig.setLog(new gdsc.core.ij.IJLogger());

		MemoryPeakResults results = fitSpot(psfStack, psfStack.getWidth(), psfStack.getHeight(), x, y);

		if (results.size() < 5)
		{
			Utils.log("  Final PSF: Not enough fit results %d", results.size());
			return 0;
		}

		// Get the results for the spot centre and width
		final double[] z = new double[results.size()];
		final double[] xCoord = new double[z.length];
		final double[] yCoord = new double[z.length];
		final double[] sd = new double[z.length];
		final double[] a = new double[z.length];

		// Set limits for the fit
		final float maxWidth = (float) (FastMath.max(fitConfig.getInitialXSD(), fitConfig.getInitialYSD()) *
				magnification * 4);
		final float maxSignal = 2; // PSF is normalised to 1  

		final WidthResultProcedure wp = new WidthResultProcedure(results, DistanceUnit.PIXEL);
		wp.getWxWy();

		final HeightResultProcedure hp = new HeightResultProcedure(results, IntensityUnit.COUNT);
		hp.getH();

		final Counter counter = new Counter();
		final Counter counterOK = new Counter();

		// We have fit the results so they will be in the preferred units 
		results.forEach(new PeakResultProcedure()
		{

			public void execute(PeakResult peak)
			{
				int i = counter.getAndIncrement();

				// Remove bad fits where the width/signal is above the expected
				final float w = FastMath.max(wp.wx[i], wp.wy[i]);
				if (peak.getSignal() > maxSignal || w > maxWidth)
					return;

				i = counterOK.getAndIncrement();
				z[i] = peak.getFrame();
				fitCom[0][peak.getFrame() - 1] = xCoord[i] = peak.getXPosition() - x;
				fitCom[1][peak.getFrame() - 1] = yCoord[i] = peak.getYPosition() - y;
				sd[i] = w;
				a[i] = hp.h[i];
			}
		});

		// Truncate
		double[] z2 = Arrays.copyOf(z, counter.getCount());
		double[] xCoord2 = Arrays.copyOf(xCoord, z2.length);
		double[] yCoord2 = Arrays.copyOf(yCoord, z2.length);
		double[] sd2 = Arrays.copyOf(sd, z2.length);
		double[] a2 = Arrays.copyOf(a, z2.length);

		// Extract the average smoothed range from the individual fits
		int r = (int) Math.ceil(averageRange / 2);
		int start = 0, stop = z2.length - 1;
		for (int j = 0; j < z2.length; j++)
		{
			if (z2[j] > cz - r)
			{
				start = j;
				break;
			}
		}
		for (int j = z2.length; j-- > 0;)
		{
			if (z2[j] < cz + r)
			{
				stop = j;
				break;
			}
		}

		// Extract xy centre coords and smooth
		double[] smoothX = new double[stop - start + 1];
		double[] smoothY = new double[smoothX.length];
		double[] smoothSd = new double[smoothX.length];
		double[] smoothA = new double[smoothX.length];
		double[] newZ = new double[smoothX.length];
		int smoothCzIndex = 0;
		for (int j = start, k = 0; j <= stop; j++, k++)
		{
			smoothX[k] = xCoord2[j];
			smoothY[k] = yCoord2[j];
			smoothSd[k] = sd2[j];
			smoothA[k] = a2[j];
			newZ[k] = z2[j];
			if (newZ[k] == cz)
				smoothCzIndex = k;
		}
		smoothX = loess.smooth(newZ, smoothX);
		smoothY = loess.smooth(newZ, smoothY);
		smoothSd = loess.smooth(newZ, smoothSd);
		smoothA = loess.smooth(newZ, smoothA);

		// Update the widths and positions using the magnification
		final double scale = 1.0 / magnification;
		for (int j = 0; j < xCoord2.length; j++)
		{
			xCoord2[j] *= scale;
			yCoord2[j] *= scale;
			sd2[j] *= scale;
		}
		for (int j = 0; j < smoothX.length; j++)
		{
			smoothX[j] *= scale;
			smoothY[j] *= scale;
			smoothSd[j] *= scale;
		}

		showPlots(z2, a2, newZ, smoothA, xCoord2, yCoord2, sd2, newZ, smoothX, smoothY, smoothSd, cz);

		// Store the data for replotting
		this.z = z2;
		this.a = a2;
		this.smoothAz = newZ;
		this.smoothA = smoothA;
		this.xCoord = xCoord2;
		this.yCoord = yCoord2;
		this.sd = sd2;
		this.newZ = newZ;
		this.smoothX = smoothX;
		this.smoothY = smoothY;
		this.smoothSd = smoothSd;

		//maximumIndex = findMinimumIndex(smoothSd, maximumIndex - start);
		return smoothSd[smoothCzIndex];
	}

	private PlotWindow getPlot(String title)
	{
		Frame f = WindowManager.getFrame(TITLE_AMPLITUDE);
		if (f != null && f instanceof PlotWindow)
			return (PlotWindow) f;
		return null;
	}

	private synchronized boolean aquirePlotLock1()
	{
		if (plotLock1)
			return false;
		return plotLock1 = true;
	}

	private synchronized boolean aquirePlotLock2()
	{
		if (plotLock2)
			return false;
		return plotLock2 = true;
	}

	private synchronized boolean aquirePlotLock3()
	{
		if (plotLock3)
			return false;
		return plotLock3 = true;
	}

	private synchronized boolean aquirePlotLock4()
	{
		if (plotLock4)
			return false;
		return plotLock4 = true;
	}

	private void createInteractivePlots(ImageStack psf, int zCentre, double nmPerPixel, double psfWidth)
	{
		this.psf = psf;
		this.zCentre = zCentre;
		this.psfNmPerPixel = nmPerPixel;
		this.psfWidth = psfWidth;

		this.slice = zCentre;
		this.distanceThreshold = psfWidth * 3;

		GenericDialog gd = new GenericDialog(TITLE);
		gd.addMessage("Plot the cumulative signal verses distance from the PSF centre.\n \nZ-centre = " + zCentre +
				"\nPSF width = " + Utils.rounded(psfWidth) + " nm");
		gd.addSlider("Slice", 1, psf.getSize(), slice);
		final double maxDistance = (psf.getWidth() / 1.414213562) * nmPerPixel;
		gd.addSlider("Distance", 0, maxDistance, distanceThreshold);
		gd.addCheckbox("Normalise", normalise);
		gd.addDialogListener(new InteractivePlotListener());
		if (!IJ.isMacro())
			drawPlots(true);
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		drawPlots(true);
	}

	private class InteractivePlotListener implements DialogListener
	{
		public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
		{
			slice = (int) gd.getNextNumber();

			double myDistanceThreshold = gd.getNextNumber();
			resetScale = resetScale || (myDistanceThreshold != distanceThreshold);
			distanceThreshold = myDistanceThreshold;

			boolean myNormalise = gd.getNextBoolean();
			resetScale = resetScale || (myNormalise != normalise);
			normalise = myNormalise;

			drawPlots(true);
			return true;
		}
	}

	private void drawPlots(boolean doSignalPlots)
	{
		updateAmplitudePlot();
		updatePSFPlot();
		if (doSignalPlots)
		{
			updateSignalAtSpecifiedSDPlot();
			updateCumulativeSignalPlot();
		}
	}

	private void updateAmplitudePlot()
	{
		if (aquirePlotLock1())
		{
			// Run in a new thread to allow the GUI to continue updating
			new Thread(new Runnable()
			{
				public void run()
				{
					try
					{
						// Continue while the parameter is changing
						boolean parametersChanged = true;
						while (parametersChanged)
						{
							// Store the parameters to be processed
							int mySlice = slice;

							// Do something with parameters
							showPlots(z, a, smoothAz, smoothA, null, null, null, null, null, null, null, mySlice);

							// Check if the parameters have changed again
							parametersChanged = (mySlice != slice);
						}
					}
					finally
					{
						// Ensure the running flag is reset
						plotLock1 = false;
					}
				}
			}).start();
		}
	}

	private void updatePSFPlot()
	{
		if (aquirePlotLock2())
		{
			// Run in a new thread to allow the GUI to continue updating
			new Thread(new Runnable()
			{
				public void run()
				{
					try
					{
						// Continue while the parameter is changing
						boolean parametersChanged = true;
						while (parametersChanged)
						{
							// Store the parameters to be processed
							int mySlice = slice;

							// Do something with parameters
							showPlots(z, null, null, null, xCoord, yCoord, sd, newZ, smoothX, smoothY, smoothSd,
									mySlice);

							// Check if the parameters have changed again
							parametersChanged = (mySlice != slice);
						}
					}
					finally
					{
						// Ensure the running flag is reset
						plotLock2 = false;
					}
				}
			}).start();
		}
	}

	private void updateSignalAtSpecifiedSDPlot()
	{
		if (aquirePlotLock3())
		{
			// Run in a new thread to allow the GUI to continue updating
			new Thread(new Runnable()
			{
				public void run()
				{
					try
					{
						// Continue while the parameter is changing
						boolean parametersChanged = true;
						while (parametersChanged)
						{
							// Store the parameters to be processed
							int mySlice = slice;

							// Do something with parameters
							plotSignalAtSpecifiedSD(psf, psfWidth / psfNmPerPixel, 3, mySlice);

							// Check if the parameters have changed again
							parametersChanged = (mySlice != slice);
						}
					}
					finally
					{
						// Ensure the running flag is reset
						plotLock3 = false;
					}
				}
			}).start();
		}
	}

	/**
	 * Show a plot of the amount of signal within N x SD for each z position. This indicates
	 * how much the PSF has spread from the original Gaussian shape.
	 * 
	 * @param psf
	 *            The PSF
	 * @param fittedSd
	 *            The width of the PSF (in pixels)
	 * @param factor
	 *            The factor to use
	 * @param slice
	 *            The slice used to create the label
	 */
	private void plotSignalAtSpecifiedSD(ImageStack psf, double fittedSd, double factor, int slice)
	{
		if (signalZ == null)
		{
			// Get the bounds
			int radius = (int) Math.round(fittedSd * factor);
			int min = FastMath.max(0, psf.getWidth() / 2 - radius);
			int max = FastMath.min(psf.getWidth() - 1, psf.getWidth() / 2 + radius);

			// Create a circle mask of the PSF projection
			ByteProcessor circle = new ByteProcessor(max - min + 1, max - min + 1);
			circle.setColor(255);
			circle.fillOval(0, 0, circle.getWidth(), circle.getHeight());
			final byte[] mask = (byte[]) circle.getPixels();

			// Sum the pixels within the mask for each slice
			signalZ = new double[psf.getSize()];
			signal = new double[psf.getSize()];
			for (int i = 0; i < psf.getSize(); i++)
			{
				double sum = 0;
				float[] data = (float[]) psf.getProcessor(i + 1).getPixels();
				for (int y = min, ii = 0; y <= max; y++)
				{
					int index = y * psf.getWidth() + min;
					for (int x = min; x <= max; x++, ii++, index++)
					{
						if (mask[ii] != 0 && data[index] > 0)
							sum += data[index];
					}
				}
				double total = 0;
				for (float f : data)
					if (f > 0)
						total += f;
				signalZ[i] = i + 1;
				signal[i] = 100 * sum / total;
			}

			signalTitle = String.format("%% PSF signal at %s x SD", Utils.rounded(factor, 3));
			signalLimits = Maths.limits(signal);
		}

		// Plot the sum
		boolean alignWindows = (WindowManager.getFrame(signalTitle) == null);

		final double total = signal[slice - 1];
		Plot2 plot = new Plot2(signalTitle, "z", "Signal", signalZ, signal);
		plot.addLabel(0, 0, String.format("Total = %s. z = %s nm", Utils.rounded(total),
				Utils.rounded((slice - zCentre) * nmPerSlice)));
		plot.setColor(Color.green);
		plot.drawLine(slice, signalLimits[0], slice, signalLimits[1]);
		plot.setColor(Color.blue);
		PlotWindow plotWindow = Utils.display(signalTitle, plot);

		if (alignWindows && plotWindow != null)
		{
			if (alignWindows && plotWindow != null)
			{
				PlotWindow otherWindow = getPlot(TITLE_AMPLITUDE);
				if (otherWindow != null)
				{
					// Put the two plots tiled together so both are visible
					Point l = plotWindow.getLocation();
					l.x = otherWindow.getLocation().x + otherWindow.getWidth();
					l.y = otherWindow.getLocation().y;
					plotWindow.setLocation(l);
				}
			}
		}
	}

	private void updateCumulativeSignalPlot()
	{
		if (aquirePlotLock4())
		{
			// Run in a new thread to allow the GUI to continue updating
			new Thread(new Runnable()
			{
				public void run()
				{
					try
					{
						// Continue while the parameter is changing
						boolean parametersChanged = true;
						while (parametersChanged)
						{
							// Store the parameters to be processed
							int mySlice = slice;
							boolean myResetScale = resetScale;
							double myDistanceThreshold = distanceThreshold;
							boolean myNormalise = normalise;

							resetScale = false;

							// Do something with parameters
							plotCumulativeSignal(mySlice, myNormalise, myResetScale, myDistanceThreshold);

							// Check if the parameters have changed again
							parametersChanged = (mySlice != slice || resetScale || myNormalise != normalise ||
									myDistanceThreshold != distanceThreshold);
						}
					}
					finally
					{
						// Ensure the running flag is reset
						plotLock4 = false;
					}
				}
			}).start();
		}
	}

	/**
	 * Show a plot of the cumulative signal vs distance from the centre
	 * 
	 * @param z
	 *            The slice to plot
	 * @param normalise
	 *            normalise the sum to 1
	 * @param resetScale
	 *            Reset the y-axis maximum
	 * @param distanceThreshold
	 *            The distance threshold for the cumulative total shown in the plot label
	 */
	private void plotCumulativeSignal(int z, boolean normalise, boolean resetScale, double distanceThreshold)
	{
		float[] data = (float[]) psf.getProcessor(z).getPixels();
		final int size = psf.getWidth();

		if (indexLookup == null || indexLookup.length != data.length)
		{
			// Precompute square distances
			double[] d2 = new double[size];
			for (int y = 0, y2 = -size / 2; y < size; y++, y2++)
				d2[y] = y2 * y2;

			// Precompute distances
			double[] d = new double[data.length];
			for (int y = 0, i = 0; y < size; y++)
			{
				for (int x = 0; x < size; x++, i++)
				{
					d[i] = Math.sqrt(d2[y] + d2[x]);
				}
			}

			// Sort
			int[] indices = SimpleArrayUtils.newArray(d.length, 0, 1);
			Sort.sortAscending(indices, d, true);

			// Store a unique cumulative index for each distance
			double lastD = d[0];
			int lastI = 0;
			int counter = 0;
			StoredData distance = new StoredData();
			indexLookup = new int[indices.length];
			for (int i = 0; i < indices.length; i++)
			{
				if (lastD != d[i])
				{
					distance.add(lastD * psfNmPerPixel);
					for (int j = lastI; j < i; j++)
					{
						indexLookup[indices[j]] = counter;
					}
					lastD = d[i];
					lastI = i;
					counter++;
				}
			}
			// Do the final distance
			distance.add(lastD * psfNmPerPixel);
			for (int j = lastI; j < indices.length; j++)
			{
				indexLookup[indices[j]] = counter;
			}
			counter++;

			distances = distance.getValues();
		}

		// Get the signal at each distance
		double[] signal = new double[distances.length];
		for (int i = 0; i < data.length; i++)
		{
			if (data[i] > 0)
				signal[indexLookup[i]] += data[i];
		}

		// Get the cumulative signal
		for (int i = 1; i < signal.length; i++)
			signal[i] += signal[i - 1];

		// Get the total up to the distance threshold
		double sum = 0;
		for (int i = 0; i < signal.length; i++)
		{
			if (distances[i] > distanceThreshold)
				break;
			sum = signal[i];
		}

		if (normalise && distanceThreshold > 0)
		{
			for (int i = 0; i < signal.length; i++)
				signal[i] /= sum;
		}

		if (resetScale)
			maxCumulativeSignal = 0;

		maxCumulativeSignal = Maths.maxDefault(maxCumulativeSignal, signal);

		String title = "Cumulative Signal";

		boolean alignWindows = (WindowManager.getFrame(title) == null);

		Plot2 plot = new Plot2(title, "Distance (nm)", "Signal", distances, signal);
		plot.setLimits(0, distances[distances.length - 1], 0, maxCumulativeSignal);
		plot.addLabel(0, 0, String.format("Total = %s (@ %s nm). z = %s nm", Utils.rounded(sum),
				Utils.rounded(distanceThreshold), Utils.rounded((z - zCentre) * nmPerSlice)));
		plot.setColor(Color.green);
		plot.drawLine(distanceThreshold, 0, distanceThreshold, maxCumulativeSignal);
		plot.setColor(Color.blue);
		PlotWindow plotWindow = Utils.display(title, plot);

		if (alignWindows && plotWindow != null)
		{
			PlotWindow otherWindow = getPlot(TITLE_PSF_PARAMETERS);
			if (otherWindow != null)
			{
				// Put the two plots tiled together so both are visible
				Point l = plotWindow.getLocation();
				l.x = otherWindow.getLocation().x + otherWindow.getWidth();
				l.y = otherWindow.getLocation().y + otherWindow.getHeight();
				plotWindow.setLocation(l);
			}
		}

		// Update the PSF to the correct slice
		if (psfImp != null)
			psfImp.setSlice(z);
	}

	private double getFWHM(ImageStack psf, int maxz)
	{
		// Extract the line profile through the stack
		int size = psf.getWidth();
		int cx = size / 2;
		// Even PSFs have the middle in the centre of two pixels
		int cx2 = (size % 2 == 0) ? cx - 1 : cx;

		double[] p0 = new double[size];
		double[] p1 = new double[size];
		double[] p2 = new double[size];
		double[] p3 = new double[size];
		double[] p4 = new double[size];
		ImageProcessor ip = psf.getProcessor(maxz);
		for (int i = 0, j = size - 1; i < size; i++, j--)
		{
			p0[i] = i;
			p1[i] = (ip.getf(i, cx) + ip.getf(i, cx2)) / 2.0;
			p2[i] = (ip.getf(cx, i) + ip.getf(cx2, i)) / 2.0;
			p3[i] = ip.getf(i, i);
			p4[i] = ip.getf(i, j);
		}

		// Find the FWHM for each line profile.
		// Diagonals need to be scaled to the appropriate distance.
		return (getFWHM(p0, p1) + getFWHM(p0, p2) + Math.sqrt(2) * getFWHM(p0, p3) + Math.sqrt(2) * getFWHM(p0, p4)) /
				4.0;
	}

	private double getFWHM(double[] x, double[] y)
	{
		// Find half max of original data
		double max = 0;
		int position = 0;
		for (int i = 0; i < y.length; i++)
		{
			if (max < y[i])
			{
				max = y[i];
				position = i;
			}
		}

		if (max == 0)
			return y.length;

		// Store half-max
		max *= 0.5;

		// The PSF profile should be a relatively straight line at half-max 
		// so no smoothing. (Note attempts to use a LOESS smoothing function failed,
		// possibly due to the small values in the y-data array)

		// Find points defining the half-max
		double p1 = 0, p2 = y.length;

		for (int i = position; i < y.length; i++)
		{
			if (y[i] < max)
			{
				// Interpolate:
				p2 = i - (max - y[i]) / (y[i - 1] - y[i]);
				break;
			}
		}
		for (int i = position; i-- > 0;)
		{
			if (y[i] < max)
			{
				// Interpolate:
				p1 = i + (max - y[i]) / (y[i + 1] - y[i]);
				break;
			}
		}

		return p2 - p1;
	}

	private HashMap<String, String> createNote()
	{
		HashMap<String, String> note = new HashMap<String, String>();
		note.put("Created", new SimpleDateFormat("d-MMM-yyyy HH:mm").format(new Date()));
		FileInfo info = imp.getOriginalFileInfo();
		if (info != null)
		{
			note.put("File", info.fileName);
			note.put("Dir", info.directory);
		}
		else
		{
			note.put("Title", imp.getTitle());
		}
		return note;
	}

	private void runUsingProjections()
	{
		if (!showProjectionDialog())
			return;

		boxRadius = (int) Math.ceil(radius);

		// Find the selected PSF spots x,y,z centre
		// We offset the centre to the middle of pixel.
		BasePoint[] centres = getSpots(0.5f);
		if (centres.length == 0)
		{
			IJ.error(TITLE, "No spots without neighbours within " + (boxRadius * 2) + "px");
			return;
		}

		// Extract the image data for processing as float
		float[][] image = CreateData.extractImageStack(imp, 0, imp.getStackSize() - 1);

		// Multi-thread for speed
		if (threadPool == null)
			threadPool = Executors.newFixedThreadPool(Prefs.getThreads());

		// Relocate the initial centres using a CoM
		centres = relocateUsingCentreOfMass(image, centres);

		// Extract each PSF into a scaled PSF
		ExtractedPSF[] psfs = extractPSFs(image, centres);

		// Iterate until centres have converged
		boolean converged = false;
		for (int iter = 0; !converged && iter < maxIterations; iter++)
		{
			// Combine all PSFs
			ExtractedPSF combined = combine(psfs);
			combined.createProjections();

			// Align each to the combined projection
			float[][] translation = align(combined, psfs);

			// Find the new centre using the old centre plus the alignment shift
			for (int j = 0; j < psfs.length; j++)
			{
				centres[j] = psfs[j].updateCentre(translation[j]);
				// Update to get the correct scale
				translation[j][0] = centres[j].getX() - psfs[j].centre.getX();
				translation[j][1] = centres[j].getY() - psfs[j].centre.getY();
				translation[j][2] = centres[j].getZ() - psfs[j].centre.getZ();
				Utils.log("[%d] Centre %d : Shift X = %s : Shift Y = %s : Shift Z = %s", iter, j + 1,
						rounder.toString(translation[j][0]), rounder.toString(translation[j][1]),
						rounder.toString(translation[j][2]));
			}

			if (interactiveMode)
			{
				combined.show("PSF");

				// Ask about each centre in turn.
				// Update Point ROI using float coordinates and set image slice to 
				// correct z-centre.
				//imp.saveRoi();
				imp.killRoi();
				ImageCanvas ic = imp.getCanvas();
				//ic.setMagnification(16);
				int reject = 0;
				Point location = null;
				float box = boxRadius + 0.5f;
				int n = imp.getStackSize();
				for (int j = 0; j < centres.length; j++)
				{
					psfs[j].show("Spot PSF");

					Overlay o = new Overlay();
					o.add(createRoi(psfs[j].centre.getX(), psfs[j].centre.getY(), Color.RED));
					float cx = centres[j].getX();
					float cy = centres[j].getY();
					o.add(createRoi(cx, cy, Color.GREEN));
					Roi roi = new Roi(cx - box, cy - box, 2 * box, 2 * box);
					o.add(roi);
					// The centre is relative to the combined PSF so use as an offset
					imp.setSlice(Maths.clip(1, n, centres[j].getZint() + 1 + n / 2));
					Rectangle r = ic.getSrcRect();
					int x = centres[j].getXint();
					int y = centres[j].getYint();
					if (!r.contains(x, y))
					{
						r.x = x - r.width / 2;
						r.y = y - r.height / 2;
						ic.setSourceRect(r);
					}
					imp.setOverlay(o);
					imp.updateAndDraw();
					NonBlockingGenericDialog gd = new NonBlockingGenericDialog(TITLE);
					gd.addMessage(String.format("Shift X = %s\nShift Y = %s\nShift Z = %s",
							rounder.toString(translation[j][0]), rounder.toString(translation[j][1]),
							rounder.toString(translation[j][2])));
					gd.enableYesNoCancel("Accept", "Reject");
					if (location != null)
						gd.setLocation(location.x, location.y);
					gd.showDialog();
					if (gd.wasCanceled())
					{
						imp.restoreRoi();
						imp.setOverlay(null);
						return;
					}
					if (!gd.wasOKed())
					{
						reject++;
						centres[j] = psfs[j].centre;
						Arrays.fill(translation[j], 0f); // For RMSD computation
					}
					location = gd.getLocation();
				}
				imp.restoreRoi();
				imp.setOverlay(null);
				if (reject == psfs.length)
				{
					IJ.error(TITLE, "No centre translations were accepted");
					return;
				}
			}
			else
			{
				if (!checkSpots(centres))
				{
					IJ.error(TITLE, "New spots centers create overlap within " + (boxRadius * 2) + "px");
					return;
				}
			}

			// Find the change in centres
			double[] rmsd = new double[2];
			for (int j = 0; j < psfs.length; j++)
			{
				rmsd[0] += Maths.pow2(translation[j][0]) + Maths.pow2(translation[j][1]);
				rmsd[1] += Maths.pow2(translation[j][2]);
			}
			for (int j = 0; j < 2; j++)
				rmsd[j] = Math.sqrt(rmsd[j] / psfs.length);

			double[] shift = combined.getCentreOfMassShift();
			double shiftd = Math.sqrt(shift[0] * shift[0] + shift[1] * shift[1]);

			Utils.log("[%d] RMSD XY = %s : RMSD Z = %s : Combined CoM shift = %s,%s (%s)", iter,
					rounder.toString(rmsd[0]), rounder.toString(rmsd[1]), rounder.toString(shift[0]),
					rounder.toString(shift[1]), rounder.toString(shiftd));

			if (interactiveMode)
			{
				// Ask if OK to continue?
				GenericDialog gd = new GenericDialog(TITLE);
				gd.addMessage(String.format("RMSD XY = %s\nRMSD Z = %s\nCombined CoM shift = %s,%s (%s)",
						rounder.toString(rmsd[0]), rounder.toString(rmsd[1]), rounder.toString(shift[0]),
						rounder.toString(shift[1]), rounder.toString(shiftd)));
				gd.enableYesNoCancel("Converged", "Continue");
				gd.showDialog();
				if (gd.wasCanceled())
					return;
				converged = gd.wasOKed();
			}
			else
			{
				// Sensible convergence on minimal shift
				converged = rmsd[0] < 0.01 && rmsd[1] < 0.05 && shiftd < 0.001;
			}

			// Update the centres using the centre-of-mass of the combined PSF
			centres = updateUsingCentreOfMassShift(shift, shiftd, combined, centres);

			// Extract each PSF into a scaled PSF
			psfs = extractPSFs(image, centres);
		}

		// Update ROI
		float[] ox = new float[centres.length];
		float[] oy = new float[centres.length];
		for (int i = 0; i < centres.length; i++)
		{
			ox[i] = centres[i].getX();
			oy[i] = centres[i].getY();
		}
		imp.setRoi(new PointRoi(ox, oy));

		// Combine all
		ExtractedPSF combined = combine(psfs);

		// Enlarge the combined PSF
		combined = combined.enlarge(psfMagnification);

		combined.createProjections();
		psfOut = combined.show("PSF");

		double[] com = combined.getCentreOfMass();
		zCentre = Maths.clip(1, combined.psf.length, (int) Math.round(1 + com[2]));

		// Show a dialog to collect processing options and 
		// find background interactively
		NonBlockingGenericDialog gd = new NonBlockingGenericDialog(TITLE);
		gd.addMessage("Configure the output PSF.\nZ-centre = " + zCentre);
		gd.addSlider("Background_slice", 1, combined.psf.length, zCentre);
		gd.addSlider("Background_cutoff (%)", 1, 50, backgroundCutoff);
		gd.addSlider("Window", 0, 0.95, windowAlpha);
		if (Utils.isShowGenericDialog())
		{
			drawPSFPlots();
			gd.addDialogListener(new InteractivePSFListener());
		}
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		drawPSFPlots();

		// Remove interactive guides		
		psfOut[0].killRoi();
		psfOut[0].setOverlay(null);
		psfOut[1].killRoi();
		psfOut[2].killRoi();

		// When click ok the background is subtracted from the PSF
		// All pixels below the background are set to zero
		// Apply a Tukey window to roll-off to zero at the outer pixels

		double[] wx = ImageWindow.tukey(combined.size, windowAlpha);
		double[] wz = ImageWindow.tukey(combined.psf.length, windowAlpha);

		// Normalise so the z centre is 1.
		float[][] psf = combined.psf;

		int cz = zCentre - 1;
		float[] data = psf[cz];
		// Apply background and window with no normalisation
		normalise(data, cz, wx, wz, 1.0);
		// Copmute normalisation from z-centre and apply
		double norm = 1.0 / Maths.sum(data);
		for (int i = 0; i < data.length; i++)
			data[i] *= norm;
		// Normalise the rest
		for (int z = 0; z < psf.length; z++)
		{
			if (z != cz)
				normalise(psf[z], z, wx, wz, norm);
		}

		// Create a new extracted PSF and show
		int magnification = combined.magnification;
		combined = new ExtractedPSF(psf, combined.size, combined.centre, magnification);
		combined.createProjections();
		psfOut = combined.show("PSF");
		psfImp = psfOut[0];
		com = combined.getCentreOfMass();

		// Add image info
		int imageCount = centres.length;
		zCentre = Maths.clip(1, combined.psf.length, (int) Math.round(1 + com[2]));
		ImagePSF.Builder imagePsf = ImagePSFHelper
				.create(zCentre, nmPerPixel / magnification, nmPerSlice, imageCount, 0, createNote()).toBuilder();
		// Add the CoM
		imagePsf.setXCentre(com[0]);
		imagePsf.setYCentre(com[1]);
		imagePsf.setZCentre(com[2]);
		// This is a bit redundant ...
		//		Offset.Builder offsetBuilder = Offset.newBuilder();
		//		offsetBuilder.setCx(com[0]);
		//		offsetBuilder.setCy(com[1]);
		//		Offset offset = offsetBuilder.build();
		//		for (int z = 1; z <= psf.length; z++)
		//			imagePsf.putOffsets(z, offset);
		psfImp.setProperty("Info", ImagePSFHelper.toString(imagePsf));

		psfImp.setRoi(new PointRoi(com[0], com[1]));
		psfImp.setSlice(zCentre);
		psfImp.resetDisplayRange();
		psfImp.updateAndDraw();

		Utils.log("Final Centre-of-mass = %s,%s,%s\n", rounder.toString(com[0]), rounder.toString(com[1]),
				rounder.toString(com[2]));
		Utils.log("%s : z-centre = %d, nm/Pixel = %s, nm/Slice = %s, %d images\n", psfImp.getTitle(), zCentre,
				Utils.rounded(nmPerPixel / magnification, 3), Utils.rounded(nmPerSlice, 3), imageCount);
	}

	/**
	 * Normalise the slice from the PSF using the XY and Z weighting window. The background is subtracted from the data,
	 * the window applied and then the result is normalised.
	 *
	 * @param data
	 *            the data
	 * @param z
	 *            the z position of the PSF data
	 * @param wx
	 *            the XY weighting
	 * @param wz
	 *            the Z weighting
	 * @param norm
	 *            the normalisation factor
	 */
	private void normalise(float[] data, int z, double[] wx, double[] wz, double norm)
	{
		// Weight by the stack position 
		double[] w;
		if (wz[z] == 1)
		{
			w = wx;
		}
		else
		{
			w = wx.clone();
			for (int i = 0; i < w.length; i++)
				w[i] *= wz[z];
		}

		final int size = wx.length;
		for (int y = 0, i = 0; y < size; y++)
		{
			double weight = w[y] * norm;
			for (int x = 0; x < size; x++, i++)
			{
				// Subtract background
				float f = (data[i] - background);
				if (f <= 0)
					data[i] = 0;
				else
					// Window function and normalise
					data[i] = (float) (f * w[x] * weight);
			}
		}
	}

	private class InteractivePSFListener implements DialogListener
	{
		public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
		{
			zCentre = (int) gd.getNextNumber();
			backgroundCutoff = gd.getNextNumber();
			windowAlpha = gd.getNextNumber();

			drawPSFPlots();
			return true;
		}
	}

	private void drawPSFPlots()
	{
		updateWindowPlot();
		updateBackgroundPlot();
		updatePSF();
	}

	double plotWindowAlpha = -1;

	private void updateWindowPlot()
	{
		if (aquirePlotLock1())
		{
			// Run in a new thread to allow the GUI to continue updating
			new Thread(new Runnable()
			{
				public void run()
				{
					try
					{
						// Continue while the parameter is changing
						while (plotWindowAlpha != windowAlpha)
						{
							// Store the parameters to be processed
							plotWindowAlpha = windowAlpha;

							int size = psfOut[0].getWidth();
							double[] x = SimpleArrayUtils.newArray(size, 0, 1.0);
							double[] w = ImageWindow.tukey(size, plotWindowAlpha);

							Plot2 plot = new Plot2(TITLE_WINDOW, "x", "Weight", x, w);
							plot.setLimits(0, size - 1, 0, 1.05);
							Utils.display(TITLE_WINDOW, plot);

							// Find the region where windowing will start
							int i = 0;
							while (i < x.length && w[i] < 1)
								i++;

							Roi roi = null;
							if (i > 0 && i != x.length)
							{
								int width = x.length - 2 * i;
								roi = new Roi(i, i, width, width);
							}
							psfOut[0].setRoi(roi);
						}
					}
					finally
					{
						// Ensure the running flag is reset
						plotLock1 = false;
					}
				}
			}).start();
		}
	}

	private int plotZCentre = -1;
	private double plotBackgroundCutoff = -1;

	private void updateBackgroundPlot()
	{
		if (aquirePlotLock2())
		{
			// Run in a new thread to allow the GUI to continue updating
			new Thread(new Runnable()
			{
				public void run()
				{
					try
					{
						// Continue while the parameter is changing
						while (plotZCentre != zCentre || plotBackgroundCutoff != backgroundCutoff)
						{
							// Store the parameters to be processed
							plotZCentre = zCentre;
							plotBackgroundCutoff = backgroundCutoff;

							// Show a plot of the background using the cumulative histogram of the
							// pixel values
							float[] pixels = (float[]) psfOut[0].getImageStack().getPixels(plotZCentre);
							double[] values = SimpleArrayUtils.toDouble(pixels);
							double p = plotBackgroundCutoff / 100;
							double[][] h = Maths.cumulativeHistogram(values, true);
							double[] x = h[0];
							double[] y = h[1];
							int i = Arrays.binarySearch(y, p);
							if (i < 0)
							{
								int after = Maths.clip(0, x.length - 1, -(i + 1));
								int before = Maths.clip(0, x.length - 1, after - 1);
								background = (float) Maths.interpolateX(x[before], y[before], x[after], y[after], p);
							}
							else
							{
								background = (float) values[i];
							}

							final byte ON = (byte) 255;
							ByteProcessor bp = new ByteProcessor(psfOut[0].getWidth(), psfOut[0].getHeight());
							bp.setLut(LUTHelper.createLUT(LutColour.BLUE, true));
							byte[] bpixels = (byte[]) bp.getPixels();
							for (int j = 0; j < pixels.length; j++)
								if (pixels[j] < background)
									bpixels[j] = ON;
							ImageRoi roi = new ImageRoi(0, 0, bp);
							roi.setZeroTransparent(true);
							Overlay o = new Overlay(roi);
							psfOut[0].setOverlay(o);

							Plot2 plot = new Plot2(TITLE_BACKGROUND, "Intensity", "Fraction", x, y);
							plot.setLimits(x[0], x[x.length - 1], 0, 1.05);
							plot.addLabel(0, 0, String.format("Background = %s @ %s %%", Utils.rounded(background),
									plotBackgroundCutoff));
							plot.setColor(Color.BLUE);
							plot.drawLine(background, 0, background, 1.05);
							plot.setColor(Color.BLACK);
							if (x[0] > 0)
								plot.setLogScaleX();
							Utils.display(TITLE_BACKGROUND, plot);
						}
					}
					finally
					{
						// Ensure the running flag is reset
						plotLock2 = false;
					}
				}
			}).start();
		}
	}

	private int psfZCentre = -1;

	private void updatePSF()
	{
		if (aquirePlotLock3())
		{
			// Run in a new thread to allow the GUI to continue updating
			new Thread(new Runnable()
			{
				public void run()
				{
					try
					{
						// Continue while the parameter is changing
						while (psfZCentre != zCentre)
						{
							// Store the parameters to be processed
							psfZCentre = zCentre;

							// Select the z-centre
							psfOut[0].setSlice(psfZCentre);
							psfOut[0].resetDisplayRange();
							psfOut[0].updateAndDraw();

							// Mark projections
							// X-projection
							psfOut[1].setRoi(new Line(psfZCentre, 0, psfZCentre, psfOut[1].getHeight()));
							// Y-projection
							psfOut[2].setRoi(new Line(0, psfZCentre, psfOut[2].getWidth(), psfZCentre));
						}
					}
					finally
					{
						// Ensure the running flag is reset
						plotLock3 = false;
					}
				}
			}).start();
		}
	}

	private BasePoint[] relocateUsingCentreOfMass(float[][] image, BasePoint[] centres)
	{
		int w = imp.getWidth();
		int h = imp.getHeight();

		// Just for getting the bounds
		ImageExtractor ie = new ImageExtractor(image[0], w, h);

		// This can be reused as a buffer
		float[][] psf = new float[image.length][];

		double mean = 0;

		for (int i = 0; i < centres.length; i++)
		{
			// Extract stack
			int x = centres[i].getXint();
			int y = centres[i].getYint();
			Rectangle bounds = ie.getBoxRegionBounds(x, y, boxRadius);
			for (int z = 0; z < image.length; z++)
				psf[z] = ImageConverter.getData(image[z], w, h, bounds, psf[z]);
			Projection p = new Projection(psf, bounds.width, bounds.height);
			double[] com = p.getCentreOfMass();
			float dx = (float) (com[0] + bounds.x - centres[i].getX());
			float dy = (float) (com[1] + bounds.y - centres[i].getY());
			float dz = (float) (com[2] - centres[i].getZ());
			Utils.log("Centre %d : %s,%s,%s updated to CoM by %s,%s,%s", i + 1, rounder.toString(centres[i].getX()),
					rounder.toString(centres[i].getY()), rounder.toString(centres[i].getZ()), rounder.toString(dx),
					rounder.toString(dy), rounder.toString(dz));
			centres[i] = centres[i].shift(dx, dy, dz);
			mean += centres[i].getZ();
		}

		// z-centres should be relative to the combined stack, not absolute
		mean /= centres.length;
		for (int i = 0; i < centres.length; i++)
			centres[i] = new BasePoint(centres[i].getX(), centres[i].getY(), (float) (centres[i].getZ() - mean));

		return centres;
	}

	private Roi createRoi(float x, float y, Color color)
	{
		Roi roi = new PointRoi(x, y);
		roi.setStrokeColor(color);
		roi.setFillColor(color);
		return roi;
	}

	private ExtractedPSF combine(ExtractedPSF[] psfs)
	{
		// All PSFs have the same size.
		// Just find the range of the relative centres;
		int min = psfs[0].relativeCentre;
		int max = min;
		for (int i = 1; i < psfs.length; i++)
		{
			min = Math.min(min, psfs[i].relativeCentre);
			max = Math.max(max, psfs[i].relativeCentre);
		}
		int range = max - min;
		int totalDepth = psfs[0].psf.length + range;
		float[][] combined = new float[totalDepth][psfs[0].psf[0].length];
		int size = psfs[0].size;
		BasePoint centre = new BasePoint(size / 2f, size / 2f, (min + max) / 2f);
		int[] count = new int[totalDepth];
		for (int i = 0; i < psfs.length; i++)
		{
			// Note: If a stack relative centre is below the centre of the stack then
			// it should be inserted later.
			int offset = max - psfs[i].relativeCentre;
			for (int j = 0; j < psfs[i].psf.length; j++)
			{
				float[] from = psfs[i].psf[j];
				float[] to = combined[j + offset];
				count[j + offset]++;
				for (int k = 0; k < to.length; k++)
					to[k] += from[k];
			}
		}
		// Q. Should the normalisation be done?
		for (int j = 0; j < combined.length; j++)
		{
			float[] to = combined[j];
			int c = count[j];
			for (int k = 0; k < to.length; k++)
				to[k] /= c;
		}
		return new ExtractedPSF(combined, size, centre, psfs[0].magnification);
	}

	private boolean showProjectionDialog()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Use X,Y,Z-projection alignment to create a combined PSF");

		gd.addNumericField("nm_per_pixel", nmPerPixel, 2);
		gd.addNumericField("nm_per_slice", nmPerSlice, 0);
		gd.addSlider("Projection_magnification", 1, 8, projectionMagnification);
		gd.addSlider("Max_iterations", 1, 20, maxIterations);
		gd.addSlider("PSF_magnification", 1, 8, psfMagnification);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		nmPerPixel = gd.getNextNumber();
		nmPerSlice = gd.getNextNumber();
		projectionMagnification = (int) gd.getNextNumber();
		maxIterations = (int) gd.getNextNumber();
		psfMagnification = (int) gd.getNextNumber();

		// Check arguments
		try
		{
			Parameters.isPositive("nm/pixel", nmPerPixel);
			Parameters.isPositive("nm/slice", nmPerSlice);
			Parameters.isEqualOrAbove("Projection magnification", projectionMagnification, 1);
			Parameters.isEqualOrAbove("Max iterations", maxIterations, 1);
			Parameters.isEqualOrAbove("PSF magnification", psfMagnification, 1);
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		return true;
	}

	private ExtractedPSF[] extractPSFs(final float[][] image, final BasePoint[] centres)
	{
		List<Future<?>> futures = new TurboList<Future<?>>(centres.length);

		final int w = imp.getWidth();
		final int h = imp.getHeight();

		final ExtractedPSF[] psfs = new ExtractedPSF[centres.length];
		for (int i = 0; i < centres.length; i++)
		{
			final int index = i;
			futures.add(threadPool.submit(new Runnable()
			{
				public void run()
				{
					psfs[index] = new ExtractedPSF(image, w, h, centres[index], boxRadius, projectionMagnification);
					// Do this here within the thread
					psfs[index].createProjections();

					//psfs[index].show(TITLE + index);
				}
			}));
		}

		Utils.waitForCompletion(futures);

		return psfs;
	}

	private static class Projection
	{
		static final int X = 0;
		static final int Y = 1;
		static final int Z = 2;

		int x, y, z;
		float[] xp, yp, zp;

		Projection(float[][] psf, int x, int y)
		{
			// Maximum project each PSF: X, Y, Z projections
			this.x = x;
			this.y = y;
			z = psf.length;
			xp = new float[y * z];
			yp = new float[x * z];
			zp = new float[x * y];
			float[] stripx = new float[x];
			float[] stripy = new float[y];
			for (int zz = 0; zz < z; zz++)
			{
				float[] data = psf[zz];

				// Z-projection
				if (zz == 0)
				{
					System.arraycopy(data, 0, zp, 0, zp.length);
				}
				else
				{
					for (int i = 0; i < zp.length; i++)
						if (zp[i] < data[i])
							zp[i] = data[i];
				}

				// X projection
				for (int yy = 0; yy < y; yy++)
				{
					System.arraycopy(data, yy * x, stripx, 0, x);
					xp[yy * z + zz] = Maths.max(stripx);
				}

				// Y projection
				for (int xx = 0; xx < x; xx++)
				{
					for (int yy = 0; yy < y; yy++)
						stripy[yy] = data[yy * x + xx];
					yp[zz * x + xx] = Maths.max(stripy);
				}
			}
		}

		int[] getDimensions()
		{
			return new int[] { x, y, z };
		}

		public double[] getCentreOfMass()
		{
			// Start in the centre of the image
			double[] com = new double[] { x / 2.0, y / 2.0, z / 2.0 };
			for (int i = 0; i < 3; i++)
			{
				FloatProcessor fp = getProjection(i);
				double[] c = getCentreOfProjection(fp);
				// Get the shift relative to the centre
				double dx = c[0] - fp.getWidth() * 0.5;
				double dy = c[1] - fp.getHeight() * 0.5;
				//System.out.printf("[%d]  [%d] %.2f = %.2f, [%d] %.2f = %.2f\n", i,
				//		Projection.getXDimension(i), c[0], dx,
				//		Projection.getYDimension(i), c[1], dy);
				// Add the shift to the current centre
				//com[Projection.getXDimension(i)] -= Projection.getXShiftDirection(i) * dx / 2;
				//com[Projection.getYDimension(i)] -= Projection.getYShiftDirection(i) * dy / 2;
				com[Projection.getXDimension(i)] += dx / 2;
				com[Projection.getYDimension(i)] += dy / 2;
			}
			return com;
		}

		public double[] getCentreOfProjection(FloatProcessor fp)
		{
			return centreOfMass((float[]) fp.getPixels(), fp.getWidth(), fp.getHeight());
		}

		private double[] centreOfMass(float[] data, int w, int h)
		{
			double cx = 0;
			double cy = 0;
			double sum = 0;
			for (int v = 0, j = 0; v < h; v++)
			{
				double sumU = 0;
				for (int u = 0; u < w; u++)
				{
					float f = data[j++];
					sumU += f;
					cx += f * u;
				}
				sum += sumU;
				cy += sumU * v;
			}
			// Find centre with 0.5 as the centre of the pixel
			cx = 0.5 + cx / sum;
			cy = 0.5 + cy / sum;
			return new double[] { cx, cy };
		}

		FloatProcessor getProjection(int i)
		{
			switch (i)
			{
				case X:
					return new FloatProcessor(z, y, xp);
				case Y:
					return new FloatProcessor(x, z, yp);
				case Z:
					return new FloatProcessor(x, y, zp);
			}
			return null;
		}

		static int getXDimension(int i)
		{
			switch (i)
			{
				case X:
					return Z;
				case Y:
					return X;
				case Z:
					return X;
			}
			return -1;
		}

		static int getYDimension(int i)
		{
			switch (i)
			{
				case X:
					return Y;
				case Y:
					return Z;
				case Z:
					return Y;
			}
			return -1;
		}

		//		static int getXShiftDirection(int i)
		//		{
		//			switch (i)
		//			{
		//				case X:
		//					return -1;
		//				case Y:
		//					return -1;
		//				case Z:
		//					return -1;
		//			}
		//			return 0;
		//		}
		//
		//		static int getYShiftDirection(int i)
		//		{
		//			switch (i)
		//			{
		//				case X:
		//					return -1;
		//				case Y:
		//					return -1;
		//				case Z:
		//					return -1;
		//			}
		//			return 0;
		//		}
	}

	private static class ExtractedPSF
	{
		BasePoint centre;
		/**
		 * The relative centre. This is just used as a relative position when combining PSFs.
		 */
		int relativeCentre;
		float[][] psf;
		int size;
		Projection projection;
		final int magnification;

		ExtractedPSF(float[][] psf, int size, BasePoint centre, int magnification)
		{
			this.centre = centre;
			this.magnification = magnification;
			relativeCentre = -1; // Not used
			this.psf = psf;
			this.size = size;
		}

		ExtractedPSF(float[][] image, int w, int h, BasePoint centre, int boxRadius, int magnification)
		{
			this.centre = centre;
			this.magnification = magnification;

			// Build using tri-cubic interpolation.

			// Create the ranges we want to interpolate. 
			// Ensure the size is odd so there is a definite centre pixel.
			size = 2 * boxRadius * magnification + 1;
			double[] y0 = new double[size];
			double[] x0 = new double[size];
			double pc = 1.0 / magnification; // Pixel centre in scaled image
			for (int i = 0, j = -size / 2; i < size; i++, j++)
			{
				double delta = pc * j;
				x0[i] = centre.getX() + delta;
				y0[i] = centre.getY() + delta;
			}

			// Extract the data for interpolation. 
			// Account for the centre of the pixel being 0.5.
			int lx = (int) Math.floor(x0[0] - 0.5);
			int ly = (int) Math.floor(y0[0] - 0.5);
			// To interpolate up to we need to have the value after.
			int ux = (int) Math.ceil(x0[size - 1] + 0.5);
			int uy = (int) Math.ceil(y0[size - 1] + 0.5);

			int rangex = ux - lx + 1;
			int rangey = uy - ly + 1;

			// Build an interpolating function
			// We pad with an extra pixel (or a duplicate pixel if at the bounds)
			double[] xval = new double[rangex];
			double[] yval = new double[rangey];
			double[] zval = new double[image.length + 2];
			double[][][] fval = new double[rangex][rangey][zval.length];
			int[] xi = new int[xval.length];
			int[] yi = new int[yval.length];
			for (int i = 0; i < xval.length; i++)
			{
				xval[i] = lx + i + 0.5;
				// Keep within data bounds for the indices
				xi[i] = Maths.clip(0, w - 1, (int) xval[i]);
			}
			for (int i = 0; i < yval.length; i++)
			{
				yval[i] = ly + i + 0.5;
				// Keep within data bounds for the indices
				yi[i] = Maths.clip(0, h - 1, (int) yval[i]);
			}

			for (int z = 0; z < zval.length; z++)
			{
				zval[z] = z - 1;

				// Keep within data bounds for the indices
				float[] data = image[Maths.clip(0, image.length - 1, z - 1)];

				for (int y = 0; y < yval.length; y++)
				{
					final int index = w * yi[y];
					for (int x = 0; x < xval.length; x++)
					{
						fval[x][y][z] = data[index + xi[x]];
					}
				}
			}

			CustomTricubicInterpolatingFunction f = new CustomTricubicInterpolator().interpolate(xval, yval, zval,
					fval);

			// Interpolate

			// Zoom the z-range too
			TDoubleArrayList list = new TDoubleArrayList(image.length * magnification);
			// The centre is a shift relative to the centre of the combined PSF in the 
			// original scale
			int cx = (int) Math.floor(centre.getZint());
			// The centre is relative to the combined stack so for interpolation
			// just get an offset from the stack centre
			double cz = (image.length - 1) / 2.0 + centre.getZ() - cx;
			// Set the relative centre after scaling
			relativeCentre = cx * magnification;
			list.add(cz);
			for (int i = 1;; i++)
			{
				double z = cz + i * pc;
				// Check if possible to interpolate
				if (z >= image.length)
					break;
				list.add(z);
			}
			for (int i = 1;; i++)
			{
				double z = cz - i * pc;
				// Check if possible to interpolate
				if (z < 0)
					break;
				list.add(z);
			}
			double[] z0 = list.toArray();
			Arrays.sort(z0);
			//relativeCentre = Arrays.binarySearch(z0, cz);

			psf = new float[z0.length][size * size];

			// Pre-compute spline positions
			IndexedCubicSplinePosition[] sx = new IndexedCubicSplinePosition[size];
			IndexedCubicSplinePosition[] sy = new IndexedCubicSplinePosition[size];
			for (int i = 0; i < size; i++)
			{
				sx[i] = f.getXSplinePosition(x0[i]);
				sy[i] = f.getYSplinePosition(y0[i]);
			}

			for (int z = 0; z < z0.length; z++)
			{
				IndexedCubicSplinePosition sz = f.getZSplinePosition(z0[z]);
				float[] data = psf[z];
				for (int y = 0, i = 0; y < size; y++)
				{
					for (int x = 0; x < size; x++, i++)
					{
						data[i] = (float) f.value(sx[x], sy[y], sz);
					}
				}
			}
		}

		void createProjections()
		{
			projection = new Projection(psf, size, size);
		}

		ImagePlus[] show(String title)
		{
			ImagePlus[] out = new ImagePlus[4];
			ImageStack stack = new ImageStack(size, size);
			for (float[] pixels : psf)
				stack.addSlice(null, pixels);
			ImagePlus imp = Utils.display(title, stack);
			// The centre is relative to the combined PSF in the original scale so use as an offset
			int n = psf.length;
			imp.setSlice(Maths.clip(1, n, centre.getZint() * magnification + 1 + n / 2));
			setCalibration(imp, 2);

			out[0] = imp;

			// Show the projections
			if (projection == null)
				return out;

			out[1] = setCalibration(Utils.display(title + " X-projection (ZY)", getProjection(0)), 0);
			out[2] = setCalibration(Utils.display(title + " Y-projection (XZ)", getProjection(1)), 1);
			out[3] = setCalibration(Utils.display(title + " Z-projection (XY)", getProjection(2)), 2);
			return out;
		}

		ImagePlus setCalibration(ImagePlus imp, int dimension)
		{
			imp.setCalibration(getCalibration(dimension));
			imp.resetDisplayRange();
			imp.updateAndDraw();
			return imp;
		}

		Calibration getCalibration(int dimension)
		{
			Calibration c = new Calibration();
			c.setUnit("nm");
			switch (dimension)
			{
				case Projection.X:
					c.pixelWidth = nmPerSlice / magnification;
					c.pixelHeight = nmPerPixel / magnification;
					break;
				case Projection.Y:
					c.pixelWidth = nmPerPixel / magnification;
					c.pixelHeight = nmPerSlice / magnification;
					break;
				case Projection.Z:
					c.pixelWidth = nmPerPixel / magnification;
					c.pixelHeight = nmPerPixel / magnification;
					c.pixelDepth = nmPerSlice / magnification;
					break;
			}
			return c;
		}

		FloatProcessor getProjection(int i)
		{
			return projection.getProjection(i);
		}

		/**
		 * Create a new centre using the shift computed from the projection.
		 *
		 * @param translation
		 *            the translation
		 * @return the new base point
		 */
		BasePoint updateCentre(float[] translation)
		{
			return new BasePoint(
					// Centre in X,Y refer to the position extracted from the image
					centre.getX() + translation[0] / magnification, centre.getY() + translation[1] / magnification,
					// The z-position is used as relative to the combined PSF 
					translation[2] / magnification);
		}

		/**
		 * Compute the centre of mass of each projection and then combine
		 *
		 * @return the centre of mass
		 */
		public double[] getCentreOfMass()
		{
			return projection.getCentreOfMass();
		}

		/**
		 * Compute the centre of mass of each projection and then the shift of the CoM from the centre of the projection
		 * image.
		 *
		 * @return the centre of mass shift
		 */
		public double[] getCentreOfMassShift()
		{
			double[] shift = projection.getCentreOfMass();
			// Turn into a shift relative to the centre
			int[] d = projection.getDimensions();
			for (int i = 0; i < 3; i++)
			{
				shift[i] -= d[i] / 2.0;
				// Account for magnification
				shift[i] /= magnification;
			}
			return shift;
		}

		public ExtractedPSF enlarge(int n)
		{
			IJ.showStatus("Enlarging ... creating spline");

			// The scales are actually arbitrary
			// We can enlarge by interpolation between the start and end
			// by evenly sampling each spline node
			double[] xval = SimpleArrayUtils.newArray(size, 0, 1.0);
			double[] yval = xval;
			double[] zval = SimpleArrayUtils.newArray(psf.length, 0, 1.0);
			double[][][] fval = new double[size][size][psf.length];
			for (int z = 0; z < psf.length; z++)
			{
				float[] data = psf[z];
				for (int y = 0, i = 0; y < size; y++)
				{
					for (int x = 0; x < size; x++, i++)
					{
						fval[x][y][z] = data[i];
					}
				}
			}

			CustomTricubicInterpolatingFunction f = new CustomTricubicInterpolator().interpolate(xval, yval, zval, fval,
					new IJTrackProgress());

			// Interpolate
			int maxx = (size - 1) * n;
			int maxy = maxx;
			int maxz = (psf.length - 1) * n;
			double step = 1.0 / n;
			float[][] psf2 = new float[maxz + 1][(maxx + 1) * (maxy + 1)];

			// Pre-compute interpolation tables
			// Use an extra one to have the final x=1 interpolation point.
			int n1 = n + 1;
			double[][] tables = new double[Maths.pow3(n1)][];
			CubicSplinePosition[] s = new CubicSplinePosition[n1];
			for (int x = 0; x < n; x++)
				s[x] = new CubicSplinePosition(x * step);
			// Final interpolation point
			s[n] = new CubicSplinePosition(1);
			for (int z = 0, i = 0; z < n1; z++)
			{
				CubicSplinePosition sz = s[z];
				for (int y = 0; y < n1; y++)
				{
					CubicSplinePosition sy = s[y];
					for (int x = 0; x < n1; x++, i++)
					{
						tables[i] = CustomTricubicFunction.computePowerTable(s[x], sy, sz);
					}
				}
			}

			//			// TODO - remove the old interpolation code
			//			// Pre-compute spline positions
			//			IndexedCubicSplinePosition[] sx = new IndexedCubicSplinePosition[maxx + 1];
			//			IndexedCubicSplinePosition[] sy = sx;
			//			double maxX = f.getMaxX();
			//			for (int i = 0; i < sx.length; i++)
			//			{
			//				sx[i] = f.getXSplinePosition(Math.min(i * step, maxX));
			//			}

			IJ.showStatus("Enlarging ... interpolating");
			//			double maxZ = f.getMaxZ();
			for (int z = 0; z <= maxz; z++)
			{
				IJ.showProgress(z, psf2.length);
				float[] data = psf2[z];

				int zposition = z / n;
				int ztable = z % n;
				if (z == maxz)
				{
					// Final interpolation point
					zposition--;
					ztable = n;
				}

				//				IndexedCubicSplinePosition sz = f.getZSplinePosition(Math.min(z * step, maxZ));
				for (int y = 0, i = 0; y <= maxy; y++)
				{
					int yposition = y / n;
					int ytable = y % n;
					if (y == maxy)
					{
						// Final interpolation point
						yposition--;
						ytable = n;
					}
					int j = n1 * (ytable + n1 * ztable);

					for (int x = 0; x <= maxx; x++, i++)
					{
						int xposition = x / n;
						int xtable = x % n;
						if (x == maxx)
						{
							// Final interpolation point
							xposition--;
							xtable = n;
						}

						//						float f1 = (float) f.value(sx[x], sy[y], sz);
						float f2 = (float) f.value(xposition, yposition, zposition, tables[j + xtable]);

						//						double e = DoubleEquality.relativeError(f1, f2);
						//						if (e > 1e-6)
						//							System.out.printf("%f vs %f  = %f\n", f1, f2, e);
						data[i] = f2;
					}
				}
			}

			IJ.showProgress(1);
			IJ.showStatus("");

			BasePoint newCentre = new BasePoint((maxx + 1) / 2.0f, (maxy + 1) / 2.0f, (maxz + 1) / 2.0f);
			return new ExtractedPSF(psf2, (maxx + 1), newCentre, magnification * n);
		}
	}

	private float[][] align(ExtractedPSF combined, final ExtractedPSF[] psfs)
	{
		int n = psfs.length * 3;
		List<Future<?>> futures = new TurboList<Future<?>>(n);

		final AlignImagesFFT[] align = new AlignImagesFFT[3];
		final Rectangle[] bounds = new Rectangle[3];
		for (int i = 0; i < 3; i++)
		{
			align[i] = new AlignImagesFFT();
			FloatProcessor fp1 = combined.getProjection(i);
			FloatProcessor fp2 = psfs[0].getProjection(i);
			align[i].init(fp1, WindowMethod.TUKEY, false);
			bounds[i] = AlignImagesFFT.createHalfMaxBounds(fp1.getWidth(), fp1.getHeight(), fp2.getWidth(),
					fp2.getHeight());
		}

		final float[][] results = new float[psfs.length][3];

		for (int j = 0; j < psfs.length; j++)
		{
			final int jj = j;
			for (int i = 0; i < 3; i++)
			{
				final int ii = i;
				futures.add(threadPool.submit(new Runnable()
				{
					public void run()
					{
						ExtractedPSF psf = psfs[jj];
						double[] result = align[ii].align(psf.getProjection(ii), WindowMethod.TUKEY, bounds[ii],
								SubPixelMethod.CUBIC);
						// We just average the shift from each projection. There should be
						// two shifts for each dimension
						results[jj][Projection.getXDimension(ii)] -= result[0] / 2;
						results[jj][Projection.getYDimension(ii)] -= result[1] / 2;
						//psfs[index].show(TITLE + index);
					}
				}));
			}
		}

		Utils.waitForCompletion(futures);

		return results;
	}

	private BasePoint[] updateUsingCentreOfMassShift(double[] shift, double shiftd, ExtractedPSF combined,
			BasePoint[] centres)
	{
		float dx = (float) shift[0];
		float dy = (float) shift[1];
		// Ignore this. We just want to keep the centres relative to the combined stack.
		// The actual z-centre of the combined stack does not matter.
		float dz = 0; //(float) shift[2]; 
		Utils.log("Combined PSF has CoM shift %s,%s (%s)", rounder.toString(shift[0]), rounder.toString(shift[1]),
				rounder.toString(shiftd));
		for (int i = 0; i < centres.length; i++)
		{
			centres[i] = centres[i].shift(dx, dy, dz);
		}
		return centres;
	}
}
