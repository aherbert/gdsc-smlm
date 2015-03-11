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

import gdsc.smlm.engine.FitEngine;
import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.engine.FitParameters;
import gdsc.smlm.engine.FitQueue;
import gdsc.smlm.engine.ParameterisedFitJob;
import gdsc.smlm.fitting.FitConfiguration;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.PSFSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.ImageConverter;
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.match.BasePoint;
import gdsc.smlm.utils.ImageExtractor;
import gdsc.smlm.utils.ImageWindow;
import gdsc.smlm.utils.Maths;
import gdsc.smlm.utils.Sort;
import gdsc.smlm.utils.Statistics;
import gdsc.smlm.utils.StoredDataStatistics;
import gdsc.smlm.utils.XmlUtils;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.Plot2;
import ij.gui.PlotWindow;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.filter.PlugInFilter;
import ij.process.Blitter;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.AWTEvent;
import java.awt.Checkbox;
import java.awt.Color;
import java.awt.Frame;
import java.awt.Point;
import java.awt.Polygon;
import java.awt.Rectangle;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.FastMath;

/**
 * Produces an average PSF image using selected diffraction limited spots from a sample image.
 * <p>
 * The input image must be a z-stack of diffraction limited spots for example quantum dots or fluorescent beads. Spots
 * will be used only when there are no spots within a specified distance to ensure a clean signal is extracted.
 */
public class PSFCreator implements PlugInFilter, ItemListener, DialogListener
{
	private final static String TITLE = "PSF Creator";
	private final static String TITLE_AMPLITUDE = "Spot Amplitude";
	private final static String TITLE_PSF_PARAMETERS = "Spot PSF";

	private static double nmPerSlice = 20;
	private static double radius = 10;
	private static double amplitudeFraction = 0.2;
	private static int startBackgroundFrames = 5;
	private static int endBackgroundFrames = 5;
	private static int magnification = 10;
	private static double smoothing = 0.25;
	private static boolean interactiveMode = false;
	private static int interpolationMethod = ImageProcessor.BICUBIC;

	private int flags = DOES_16 | DOES_8G | DOES_32 | NO_CHANGES;
	private ImagePlus imp, psfImp;
	private double nmPerPixel;
	private FitEngineConfiguration config = null;
	private FitConfiguration fitConfig;
	private int boxRadius;
	private static Point yesNoPosition = null;

	private ExecutorService threadPool = null;
	private double progress = 0;

	// Private variables that are used during background threaded plotting of the cumulative signal 
	private ImageStack psf = null;
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
	private double maxy = 1;
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
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Produces an average PSF using selected diffraction limited spots.\nUses the current fit configuration to fit spots.");

		gd.addCheckbox("Update_Fit_Configuration", false);
		gd.addNumericField("nm_per_slice", nmPerSlice, 0);
		gd.addSlider("Radius", 3, 20, radius);
		gd.addSlider("Amplitude_fraction", 0.01, 0.5, amplitudeFraction);
		gd.addSlider("Start_background_frames", 1, 20, startBackgroundFrames);
		gd.addSlider("End_background_frames", 1, 20, endBackgroundFrames);
		gd.addSlider("Magnification", 5, 15, magnification);
		gd.addSlider("Smoothing", 0.25, 0.5, smoothing);
		gd.addCheckbox("Interactive_mode", interactiveMode);
		String[] methods = ImageProcessor.getInterpolationMethods();
		gd.addChoice("Interpolation", methods, methods[interpolationMethod]);

		((Checkbox) gd.getCheckboxes().get(0)).addItemListener(this);

		gd.showDialog();

		if (gd.wasCanceled())
			return DONE;

		gd.getNextBoolean();
		nmPerSlice = gd.getNextNumber();
		radius = gd.getNextNumber();
		amplitudeFraction = gd.getNextNumber();
		startBackgroundFrames = (int) gd.getNextNumber();
		endBackgroundFrames = (int) gd.getNextNumber();
		magnification = (int) gd.getNextNumber();
		smoothing = gd.getNextNumber();
		interactiveMode = gd.getNextBoolean();
		interpolationMethod = gd.getNextChoiceIndex();

		// Check arguments
		try
		{
			Parameters.isPositive("nm/slice", nmPerSlice);
			Parameters.isAbove("Radius", radius, 2);
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
			return DONE;
		}

		return flags;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
	 */
	public void run(ImageProcessor ip)
	{
		loadConfiguration();
		BasePoint[] spots = getSpots();
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
		fitConfig.setDuplicateDistance(0);

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
		for (int n = 1; n <= spots.length; n++)
		{
			BasePoint spot = spots[n - 1];
			final int x = (int) spot.getX();
			final int y = (int) spot.getY();

			MemoryPeakResults results = fitSpot(stack, width, height, x, y);

			if (results.size() < 5)
			{
				Utils.log("  Spot %d: Not enough fit results %d", n, results.size());
				continue;
			}

			// Get the results for the spot centre and width
			double[] z = new double[results.size()];
			double[] xCoord = new double[z.length];
			double[] yCoord = new double[z.length];
			double[] sd = new double[z.length];
			double[] a = new double[z.length];
			int i = 0;
			for (PeakResult peak : results.getResults())
			{
				z[i] = peak.peak;
				xCoord[i] = peak.getXPosition() - x;
				yCoord[i] = peak.getYPosition() - y;
				sd[i] = FastMath.max(peak.getXSD(), peak.getYSD());
				a[i] = peak.getAmplitude();
				i++;
			}

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

			// Store result
			Utils.log("  Spot %d => x=%.2f, y=%.2f, z=%d, sd=%.2f, A=%.2f\n", n, cx, cy, cz, csd, ca);
			centres.add(new double[] { cx, cy, cz, csd });
		}

		if (interactiveMode)
		{
			imp.setSlice(currentSlice);
			imp.setOverlay(null);
		}

		if (centres.isEmpty())
		{
			String msg = "No suitable spots could be identified centres";
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

			float b = getBackground(spot);
			stats.add(b);
			subtractBackgroundAndWindow(spot, b, regionBounds.width, regionBounds.height);

			// This takes a long time so this should track progress
			ok = addToPSF(maxz, magnification, psf, centre, spot, regionBounds, progress, increment);
		}

		IJ.showProgress(1);
		if (threadPool != null)
		{
			threadPool.shutdownNow();
			threadPool = null;
		}

		if (!ok)
			return;

		final double avSd = getAverage(averageSd, averageA, 2);
		Utils.log("  Average background = %.2f, Av. SD = %s px", stats.getMean(), Utils.rounded(avSd, 4));

		normalise(psf, maxz, avSd * magnification);
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
		double[] slice = Utils.newArray(psf.getSize(), 1, 1.0);
		String title = TITLE + " CoM Drift";
		Plot2 plot = new Plot2(title, "Slice", "Drift (nm)");
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
		psfImp.setProperty("Info",
				XmlUtils.toXML(new PSFSettings(maxz, nmPerPixel / magnification, nmPerSlice, centres.size(), fwhm)));

		Utils.log("%s : z-centre = %d, nm/Pixel = %s, nm/Slice = %s, %d images, PSF SD = %s nm, FWHM = %s px\n",
				psfImp.getTitle(), maxz, Utils.rounded(nmPerPixel / magnification, 3), Utils.rounded(nmPerSlice, 3),
				centres.size(), Utils.rounded(fittedSd * nmPerPixel, 4), Utils.rounded(fwhm));

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
		double iqr = upper - lower;
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
		results.setSortAfterEnd(true);
		results.begin();
		FitEngine engine = new FitEngine(config, results, Prefs.getThreads(), FitQueue.BLOCKING);

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

	private boolean ignoreSpot(int n, final double[] z, final double[] a, final double[] smoothA,
			final double[] xCoord, final double[] yCoord, final double[] sd, final double[] newZ,
			final double[] smoothX, final double[] smoothY, double[] smoothSd, final double cx, final double cy,
			final int cz, double csd)
	{
		// Allow an interactive mode that shows the plots and allows the user to Yes/No
		// the addition of the data
		if (interactiveMode)
		{
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
			gd.addMessage(String.format("Add spot %d to the PSF?\n \nx = %.2f\ny = %.2f\nz = %d\nsd = %.2f\n", n, cx,
					cy, cz, csd));
			if (yesNoPosition != null)
			{
				gd.centerDialog(false);
				gd.setLocation(yesNoPosition);
			}

			gd.showDialog();

			yesNoPosition = gd.getLocation();
			return !gd.wasOKed();
		}
		return false;
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
			plot.addLabel(
					0,
					0,
					String.format("Amplitude = %s (%sx). z = %s nm", Utils.rounded(amplitude),
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

	private float getBackground(float[][] spot)
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
		//Utils.log("First %d = %.2f, Last %d = %.2f, av = %.2f", backgroundFrames, first.getMean(), backgroundFrames,
		//		last.getMean(), av);
		return av;
	}

	@SuppressWarnings("unused")
	private float getBackground(final double fraction, StoredDataStatistics all)
	{
		double[] allValues = all.getValues();
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

	private void subtractBackgroundAndWindow(float[][] spot, final float min, final int spotWidth, final int spotHeight)
	{
		//ImageWindow imageWindow = new ImageWindow();
		for (int i = 0; i < spot.length; i++)
		{
			for (int j = 0; j < spot[i].length; j++)
				spot[i][j] = FastMath.max(spot[i][j] - min, 0);

			// Use a Tukey window to roll-off the image edges
			//spot[i] = imageWindow.applySeperable(spot[i], spotWidth, spotHeight, ImageWindow.WindowFunction.Tukey);
			spot[i] = ImageWindow.applyWindow(spot[i], spotWidth, spotHeight, ImageWindow.WindowFunction.TUKEY);
		}
	}

	private boolean addToPSF(int maxz, final int magnification, ImageStack psf, double[] centre, final float[][] spot,
			final Rectangle regionBounds, double progress, final double increment)
	{
		// Calculate insert point in enlargement
		final int x = (int) centre[0];
		final int y = (int) centre[1];
		final int z = (int) centre[2];

		final double insertX = getInsert(centre[0], x, magnification);
		final double insertY = getInsert(centre[1], y, magnification);

		int insertZ = maxz - z + 1;
		//Utils.log("Insert point = %.2f,%.2f => %.2f,%.2f\n", centre[0] - x, centre[1] - y, insertX, insertY);

		// Copy the processor using a weighted image
		final int lowerX = (int) insertX;
		final int lowerY = (int) insertY;

		final double wx2 = insertX - lowerX;
		final double wx1 = 1 - wx2;
		final double wy2 = insertY - lowerY;
		final double wy1 = 1 - wy2;

		// Enlargement size
		final int dstWidth = regionBounds.width * magnification;
		final int dstHeight = regionBounds.height * magnification;

		// Multi-thread for speed
		if (threadPool == null)
			threadPool = Executors.newFixedThreadPool(Prefs.getThreads());

		List<Future<?>> futures = new LinkedList<Future<?>>();

		for (int i = 0; i < spot.length; i++)
		{
			final ImageProcessor ip = psf.getProcessor(insertZ++);
			final float[] spotData = spot[i];

			futures.add(threadPool.submit(new Runnable()
			{
				public void run()
				{
					if (Utils.isInterrupted())
						return;

					incrementProgress(increment);

					// Enlarge
					FloatProcessor fp = new FloatProcessor(regionBounds.width, regionBounds.height, spotData, null);
					fp.setInterpolationMethod(interpolationMethod);
					fp = (FloatProcessor) fp.resize(dstWidth, dstHeight);

					// Add to the combined PSF using the correct offset
					//psf.getProcessor(insertZ++).copyBits(fp, (int) Math.round(insertX), (int) Math.round(insertY), Blitter.ADD);

					// Add to the combined PSF using the correct offset
					copyBits(ip, fp, lowerX, lowerY, wx1 * wy1);
					copyBits(ip, fp, lowerX + 1, lowerY, wx2 * wy1);
					copyBits(ip, fp, lowerX, lowerY + 1, wx1 * wy2);
					copyBits(ip, fp, lowerX + 1, lowerY + 1, wx2 * wy2);
				}
			}));

			if (Utils.isInterrupted())
				break;
		}

		Utils.waitForCompletion(futures);

		return !Utils.isInterrupted();
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
	 * Normalise the PSF so the height of the specified frame foreground pixels is 1.
	 * <p>
	 * Assumes the PSF can be approximated by a Gaussian in the central frame. All pixels within 3 sigma of the centre
	 * are foreground pixels.
	 * 
	 * @param psf
	 * @param n
	 *            The frame number
	 * @param sigma
	 *            the Gaussian standard deviation (in pixels)
	 */
	public static void normalise(ImageStack psf, int n, double sigma)
	{
		if (psf == null || psf.getSize() == 0)
			return;
		if (!(psf.getPixels(1) instanceof float[]))
			return;
		double cx = psf.getWidth() * 0.5;

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
			d2[i] = (x - cx) * (x - cx);

		for (int y = min, i = 0; y <= max; y++, i++)
		{
			int index = y * psf.getWidth() + min;
			for (int x = min, j = 0; x <= max; x++, index++, j++)
			{
				// Check if the pixel is within 3 sigma of the centre
				if (d2[i] + d2[j] < r2)
				{
					foregroundSum += data[index];
					foregroundN++;
				}
			}
		}

		// Get the average background
		double backgroundSum = new Statistics(data).getSum() - foregroundSum;
		double background = backgroundSum / (data.length - foregroundN);

		// Subtract the background from the foreground sum
		foregroundSum -= background * foregroundN;

		for (int i = 0; i < psf.getSize(); i++)
		{
			data = (float[]) psf.getPixels(i + 1);
			for (int j = 0; j < data.length; j++)
				data[j] = (float) ((data[j] - background) / foregroundSum);
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
		double sum = new Statistics((float[]) psf.getPixels(n)).getSum();
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
		final int offset = psf.getWidth() / 2;
		for (int i = 0; i < size; i++)
		{
			float[] data = (float[]) psf.getPixels(i + 1);
			double sumX = 0, sumY = 0, sum = 0;
			for (int y = 0, j = 0; y < psf.getHeight(); y++)
			{
				for (int x = 0; x < psf.getWidth(); x++, j++)
				{
					sum += data[j];
					sumX += x * data[j];
					sumY += y * data[j];
				}
			}
			double comX = sumX / sum - offset;
			double comY = sumY / sum - offset;
			com[0][i] = comX;
			com[1][i] = comY;
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
		//		double[] slice = Utils.newArray(psf.getSize(), 1, 1.0);
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

	private void loadConfiguration()
	{
		String filename = SettingsManager.getSettingsFilename();
		GlobalSettings settings = SettingsManager.loadSettings(filename);
		nmPerPixel = settings.getCalibration().nmPerPixel;
		config = settings.getFitEngineConfiguration();
		fitConfig = config.getFitConfiguration();
		if (radius < 5 * FastMath.max(fitConfig.getInitialPeakStdDev0(), fitConfig.getInitialPeakStdDev1()))
		{
			radius = 5 * FastMath.max(fitConfig.getInitialPeakStdDev0(), fitConfig.getInitialPeakStdDev1());
			Utils.log("Radius is less than 5 * PSF standard deviation, increasing to %s", Utils.rounded(radius));
		}
		boxRadius = (int) Math.ceil(radius);
	}

	/**
	 * @return Extract all the ROI points that are not within twice the box radius of any other spot
	 */
	private BasePoint[] getSpots()
	{
		Roi roi = imp.getRoi();
		if (roi != null && roi.getType() == Roi.POINT)
		{
			Polygon p = ((PolygonRoi) roi).getNonSplineCoordinates();
			int n = p.npoints;
			Rectangle bounds = roi.getBounds();

			BasePoint[] roiPoints = new BasePoint[n];
			for (int i = 0; i < n; i++)
			{
				roiPoints[i] = new BasePoint(bounds.x + p.xpoints[i], bounds.y + p.ypoints[i], 0);
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

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.awt.event.ItemListener#itemStateChanged(java.awt.event.ItemEvent)
	 */
	public void itemStateChanged(ItemEvent e)
	{
		// Run the fit configuration plugin to update the settings.
		if (e.getSource() instanceof Checkbox)
		{
			((Checkbox) e.getSource()).setState(false);
			IJ.run("Fit Configuration");
		}
	}

	/**
	 * Fit the new PSF image and show a graph of the amplitude/width
	 * 
	 * @param psf
	 * @param loess
	 * @param averageRange
	 * @param fitCom
	 * @return The width of the PSF in the z-centre
	 */
	private double fitPSF(ImageStack psf, LoessInterpolator loess, int cz, double averageRange, double[][] fitCom)
	{
		IJ.showStatus("Fitting final PSF");
		// Update the box radius since this is used in the fitSpot method.
		boxRadius = psf.getWidth() / 2;
		int x = boxRadius, y = boxRadius;
		FitConfiguration fitConfig = config.getFitConfiguration();
		final double shift = fitConfig.getCoordinateShiftFactor();
		fitConfig.setInitialPeakStdDev0(fitConfig.getInitialPeakStdDev0() * magnification);
		fitConfig.setInitialPeakStdDev1(fitConfig.getInitialPeakStdDev1() * magnification);
		// Need to be updated after the widths have been set
		fitConfig.setCoordinateShiftFactor(shift);
		fitConfig.setBackgroundFitting(false);
		//fitConfig.setLog(new IJLogger());

		MemoryPeakResults results = fitSpot(psf, psf.getWidth(), psf.getHeight(), x, y);

		if (results.size() < 5)
		{
			Utils.log("  Final PSF: Not enough fit results %d", results.size());
			return 0;
		}

		// Get the results for the spot centre and width
		double[] z = new double[results.size()];
		double[] xCoord = new double[z.length];
		double[] yCoord = new double[z.length];
		double[] sd = new double[z.length];
		double[] a = new double[z.length];
		int i = 0;

		// Set limits for the fit
		final float maxWidth = (float) (FastMath.max(fitConfig.getInitialPeakStdDev0(),
				fitConfig.getInitialPeakStdDev1()) *
				magnification * 4);
		final float maxSignal = 2; // PSF is normalised to 1  

		for (PeakResult peak : results.getResults())
		{
			// Remove bad fits where the width/signal is above the expected
			final float w = FastMath.max(peak.getXSD(), peak.getYSD());
			if (peak.getSignal() > maxSignal || w > maxWidth)
				continue;

			z[i] = peak.peak;
			fitCom[0][peak.peak - 1] = xCoord[i] = peak.getXPosition() - x;
			fitCom[1][peak.peak - 1] = yCoord[i] = peak.getYPosition() - y;
			sd[i] = w;
			a[i] = peak.getAmplitude();
			i++;
		}

		// Truncate
		z = Arrays.copyOf(z, i);
		xCoord = Arrays.copyOf(xCoord, i);
		yCoord = Arrays.copyOf(yCoord, i);
		sd = Arrays.copyOf(sd, i);
		a = Arrays.copyOf(a, i);

		// Extract the average smoothed range from the individual fits
		int r = (int) Math.ceil(averageRange / 2);
		int start = 0, stop = z.length - 1;
		for (int j = 0; j < z.length; j++)
		{
			if (z[j] > cz - r)
			{
				start = j;
				break;
			}
		}
		for (int j = z.length; j-- > 0;)
		{
			if (z[j] < cz + r)
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
			smoothX[k] = xCoord[j];
			smoothY[k] = yCoord[j];
			smoothSd[k] = sd[j];
			smoothA[k] = a[j];
			newZ[k] = z[j];
			if (newZ[k] == cz)
				smoothCzIndex = k;
		}
		smoothX = loess.smooth(newZ, smoothX);
		smoothY = loess.smooth(newZ, smoothY);
		smoothSd = loess.smooth(newZ, smoothSd);
		smoothA = loess.smooth(newZ, smoothA);

		// Update the widths and positions using the magnification
		final double scale = 1.0 / magnification;
		for (int j = 0; j < xCoord.length; j++)
		{
			xCoord[j] *= scale;
			yCoord[j] *= scale;
			sd[j] *= scale;
		}
		for (int j = 0; j < smoothX.length; j++)
		{
			smoothX[j] *= scale;
			smoothY[j] *= scale;
			smoothSd[j] *= scale;
		}

		showPlots(z, a, newZ, smoothA, xCoord, yCoord, sd, newZ, smoothX, smoothY, smoothSd, cz);

		// Store the data for replotting
		this.z = z;
		this.a = a;
		this.smoothAz = newZ;
		this.smoothA = smoothA;
		this.xCoord = xCoord;
		this.yCoord = yCoord;
		this.sd = sd;
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
		gd.addDialogListener(this);
		if (!IJ.isMacro())
			drawPlots();
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		drawPlots();
	}

	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
	{
		slice = (int) gd.getNextNumber();

		double myDistanceThreshold = gd.getNextNumber();
		resetScale = resetScale || (myDistanceThreshold != distanceThreshold);
		distanceThreshold = myDistanceThreshold;

		boolean myNormalise = gd.getNextBoolean();
		resetScale = resetScale || (myNormalise != normalise);
		normalise = myNormalise;

		drawPlots();
		return true;
	}

	private void drawPlots()
	{
		updateAmplitudePlot();
		updatePSFPlot();
		updateSignalAtSpecifiedSDPlot();
		updateCumulativeSignalPlot();
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
		plot.addLabel(
				0,
				0,
				String.format("Total = %s. z = %s nm", Utils.rounded(total),
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
							parametersChanged = (mySlice != slice || resetScale || myNormalise != normalise || myDistanceThreshold != distanceThreshold);
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
			int[] indices = Utils.newArray(d.length, 0, 1);
			Sort.sort(indices, d, true);

			// The sort is made in descending order so invert
			Sort.reverse(indices);
			Sort.reverse(d);

			// Store a unique cumulative index for each distance
			double lastD = d[0];
			int lastI = 0;
			int counter = 0;
			StoredDataStatistics distance = new StoredDataStatistics();
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

		if (normalise)
		{
			for (int i = 0; i < signal.length; i++)
				signal[i] /= sum;
		}

		if (resetScale)
			maxy = 0;

		maxy = Maths.maxDefault(maxy, signal);

		String title = "Cumulative signal";

		boolean alignWindows = (WindowManager.getFrame(title) == null);

		Plot2 plot = new Plot2(title, "Distance (nm)", "Signal", distances, signal);
		plot.setLimits(0, distances[distances.length - 1], 0, maxy);
		plot.addLabel(0, 0, String.format("Total = %s (@ %s nm). z = %s nm", Utils.rounded(sum),
				Utils.rounded(distanceThreshold), Utils.rounded((z - zCentre) * nmPerSlice)));
		plot.setColor(Color.green);
		plot.drawLine(distanceThreshold, 0, distanceThreshold, maxy);
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
		return (getFWHM(p0, p1) + getFWHM(p0, p2) + Math.sqrt(2) * getFWHM(p0, p3) + Math.sqrt(2) * getFWHM(p0, p4)) / 4.0;
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
}
