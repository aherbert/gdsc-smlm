package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.Checkbox;
import java.awt.Color;
import java.awt.Frame;
import java.awt.Label;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
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

import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.FastMath;

import gdsc.core.data.FloatStackTrivalueProvider;
import gdsc.core.data.procedures.FloatStackTrivalueProcedure;
import gdsc.core.data.utils.Rounder;
import gdsc.core.data.utils.RounderFactory;
import gdsc.core.data.utils.TypeConverter;
import gdsc.core.ij.AlignImagesFFT;
import gdsc.core.ij.AlignImagesFFT.SubPixelMethod;
import gdsc.core.ij.AlignImagesFFT.WindowMethod;
import gdsc.core.ij.IJTrackProgress;
import gdsc.core.ij.Utils;
import gdsc.core.match.BasePoint;
import gdsc.core.math.interpolation.CubicSplinePosition;
import gdsc.core.math.interpolation.CustomTricubicFunction;
import gdsc.core.math.interpolation.CustomTricubicInterpolatingFunction;
import gdsc.core.math.interpolation.CustomTricubicInterpolatingFunction.Size;
import gdsc.core.math.interpolation.CustomTricubicInterpolator;
import gdsc.core.utils.DoubleData;
import gdsc.core.utils.ImageExtractor;
import gdsc.core.utils.ImageWindow;
import gdsc.core.utils.Maths;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.core.utils.Sort;
import gdsc.core.utils.Statistics;
import gdsc.core.utils.StoredData;
import gdsc.core.utils.StoredDataStatistics;
import gdsc.core.utils.TextUtils;
import gdsc.core.utils.TurboList;
import gdsc.smlm.data.config.CalibrationReader;
import gdsc.smlm.data.config.CalibrationWriter;
import gdsc.smlm.data.config.FitProtos.FitSolver;
import gdsc.smlm.data.config.FitProtosHelper;
import gdsc.smlm.data.config.GUIProtos.PSFCreatorSettings;
import gdsc.smlm.data.config.GUIProtosHelper;
import gdsc.smlm.data.config.PSFProtos.ImagePSF;
import gdsc.smlm.data.config.PSFProtos.PSF;
import gdsc.smlm.data.config.PSFProtos.PSFParameter;
import gdsc.smlm.data.config.PSFProtos.PSFParameterUnit;
import gdsc.smlm.data.config.UnitConverterFactory;
import gdsc.smlm.data.config.UnitProtos.AngleUnit;
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
import gdsc.smlm.filters.BlockMeanFilter;
import gdsc.smlm.function.Erf;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.ij.plugins.CubicSplineManager.CubicSplinePSF;
import gdsc.smlm.ij.settings.ImagePSFHelper;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.Image2DAligner;
import gdsc.smlm.ij.utils.Image3DAligner;
import gdsc.smlm.ij.utils.ImageConverter;
import gdsc.smlm.model.camera.CameraModel;
import gdsc.smlm.model.camera.FixedPixelCameraModel;
import gdsc.smlm.results.Counter;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.SynchronizedPeakResults;
import gdsc.smlm.results.procedures.HeightResultProcedure;
import gdsc.smlm.results.procedures.PeakResultProcedure;
import gdsc.smlm.results.procedures.WidthResultProcedure;
import gdsc.smlm.utils.Tensor2D;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.ExtendedGenericDialog;
import ij.gui.ExtendedGenericDialog.OptionCollectedEvent;
import ij.gui.ExtendedGenericDialog.OptionCollectedListener;
import ij.gui.ExtendedGenericDialog.OptionListener;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.Line;
import ij.gui.NonBlockingExtendedGenericDialog;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.Plot2;
import ij.gui.PlotWindow;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.io.FileInfo;
import ij.measure.Calibration;
import ij.plugin.WindowOrganiser;
import ij.plugin.filter.PlugInFilter;
import ij.process.Blitter;
import ij.process.ByteProcessor;
import ij.process.FloatPolygon;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

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
	private final static String TITLE_BACKGROUND = "PSF Background";
	private final static String TITLE_FOREGROUND = "PSF Foreground";
	private final static String TITLE_SIGNAL = "PSF Signal";
	private final static String TITLE_HWHM = "PSF HWHM";
	private final static String TITLE_ANGLE = "PSF Angle";
	//private final static String TITLE_CONTRAST = "PSF Contrast";
	private final static String TITLE_PSF = "PSF";
	private final static String TITLE_SPOT_PSF = "Spot PSF";

	private final static String[] MODE = { "Stack Alignment", "Gaussian Fitting" };
	private final static int MODE_ALIGNMENT = 0;
	private final static int MODE_FITTING = 1;
	private final static String[] ALIGNMENT_MODE = { "2D Projections", "3D" };
	private final static int ALIGNMENT_MODE_2D = 0;
	private static String[] PSF_TYPE = { "Spot", "Astigmatism", "Double Helix" };
	@SuppressWarnings("unused")
	private static final int PSF_TYPE_SPOT = 0;
	private static final int PSF_TYPE_ASTIGMATISM = 1;
	private static final int PSF_TYPE_DH = 2;
	private static final String[] OUTPUT_TYPE = { "CSpline", "Image PSF" };
	private static final int OUTPUT_TYPE_CSPLINE = 0;
	private static final int OUTPUT_TYPE_IMAGE_PSF = 1;

	private PSFCreatorSettings.Builder settings;

	private int flags = DOES_16 | DOES_8G | DOES_32 | NO_CHANGES;
	private ImagePlus imp, psfImp;

	private static Rounder rounder = RounderFactory.create(4);
	private PSFCentreSelector zSelector;

	private FitEngineConfiguration config = null;
	private FitConfiguration fitConfig;
	private double nmPerPixel;
	private int boxRadius, zRadius;
	private static Point yesNoPosition = null;

	private ExecutorService threadPool = null;
	private double progress = 0;

	// Private variables that are used during background threaded plotting of the cumulative signal 
	private ImageStack psf = null;
	private ImagePlus[] psfOut = null;
	private int zCentre = 0;
	private double psfWidth = 0;

	// Cache settings for convenience
	private double psfNmPerPixel = 0;
	private boolean checkAlignments = false;

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
		NonBlockingExtendedGenericDialog gd = new NonBlockingExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		settings = SettingsManager.readPSFCreatorSettings(0).toBuilder();

		guessScale();

		gd.addMessage("Produces an average PSF using selected diffraction limited spots.");

		gd.addChoice("Mode", MODE, settings.getMode());
		gd.addSlider("Radius (px)", 3, Maths.max(10, imp.getWidth(), imp.getHeight()) / 2, settings.getRadius());
		gd.addCheckbox("Interactive_mode", settings.getInteractiveMode());

		InteractiveInputListener l = new InteractiveInputListener();
		gd.addDialogListener(l);

		gd.showDialog();

		// Clear the bounding box
		if (plotRadius != -1)
			imp.setOverlay(null);

		SettingsManager.writeSettings(settings);

		if (gd.wasCanceled())
			return DONE;

		// Check arguments
		try
		{
			Parameters.isAbove("Radius", settings.getRadius(), 2);
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
		CalibrationWriter cw = CalibrationWriter.create(settings.getCalibration());
		// It does not matter if we already have settings, try and update them anyway
		//if (cw.getNmPerPixel() == 0 || settings.getNmPerSlice() == 0)
		{
			Calibration c = imp.getCalibration();
			double r = guessScale(c.getXUnit(), c.pixelWidth);
			if (r != 0)
			{
				cw.setNmPerPixel(r);
				settings.setCalibration(cw.getBuilder());
			}
			r = guessScale(c.getZUnit(), c.pixelDepth);
			if (r != 0)
				settings.setNmPerSlice(r);
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

	private class InteractiveInputListener implements DialogListener
	{
		final boolean draw;;

		InteractiveInputListener()
		{
			draw = Utils.isShowGenericDialog();
			if (draw)
			{
				imp.setSlice(imp.getStackSize() / 2);
				drawBoundingBox();
			}
		}

		public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
		{
			settings.setMode(gd.getNextChoiceIndex());
			settings.setRadius(gd.getNextNumber());
			settings.setInteractiveMode(gd.getNextBoolean());

			if (draw)
				drawBoundingBox();

			return settings.getRadius() >= 2;
		}
	}

	double plotRadius = -1;

	private void drawBoundingBox()
	{
		if (aquirePlotLock1())
		{
			// Get the spots here as the user may want to interactively pick new ones 
			final BasePoint[] points = getSpots();
			final Rectangle[] bounds = new Rectangle[points.length];
			final Roi[] rois = new Roi[points.length];

			// Run in a new thread to allow the GUI to continue updating
			new Thread(new Runnable()
			{
				public void run()
				{
					try
					{
						// Continue while the parameter is changing
						while (plotRadius != settings.getRadius())
						{
							// Store the parameters to be processed
							plotRadius = settings.getRadius();
							int boxRadius = (int) Math.ceil(plotRadius);
							int w = 2 * boxRadius + 1;

							for (int i = 0; i < points.length; i++)
							{
								BasePoint p = points[i];
								int cx = p.getXint();
								int cy = p.getYint();
								Roi r = new Roi(cx - plotRadius, cy - plotRadius, w, w);
								bounds[i] = r.getBounds();
								rois[i] = r;
								// Check for overlap
								for (int j = i; j-- > 0;)
								{
									if (bounds[i].intersects(bounds[j]))
									{
										rois[j].setStrokeColor(Color.RED);
										r.setStrokeColor(Color.RED);
									}
								}
							}

							Overlay o = new Overlay();
							for (Roi roi : rois)
								o.add(roi);
							imp.setOverlay(o);
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

	private boolean showFittingDialog()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Use PSF fitting to create a combined PSF");

		gd.addNumericField("nm_per_slice", settings.getNmPerSlice(), 0);
		gd.addSlider("Amplitude_fraction", 0.01, 0.5, settings.getAmplitudeFraction());
		gd.addSlider("Start_background_frames", 1, 20, settings.getStartBackgroundFrames());
		gd.addSlider("End_background_frames", 1, 20, settings.getEndBackgroundFrames());
		gd.addSlider("Magnification", 5, 15, settings.getMagnification());
		gd.addSlider("Smoothing", 0.25, 0.5, settings.getSmoothing());
		gd.addCheckbox("Centre_each_slice", settings.getCentreEachSlice());
		gd.addNumericField("CoM_cut_off", settings.getComCutOff(), -2);
		String[] methods = ImageProcessor.getInterpolationMethods();
		gd.addChoice("Interpolation", methods, methods[settings.getInterpolationMethod()]);

		gd.showDialog();

		SettingsManager.writeSettings(settings);

		if (gd.wasCanceled())
			return false;

		settings.setNmPerSlice(gd.getNextNumber());
		settings.setAmplitudeFraction(gd.getNextNumber());
		settings.setStartBackgroundFrames((int) gd.getNextNumber());
		settings.setEndBackgroundFrames((int) gd.getNextNumber());
		settings.setMagnification((int) gd.getNextNumber());
		settings.setSmoothing(gd.getNextNumber());
		settings.setCentreEachSlice(gd.getNextBoolean());
		settings.setComCutOff(Maths.max(0, gd.getNextNumber()));
		settings.setInterpolationMethod(gd.getNextChoiceIndex());

		// Check arguments
		try
		{
			Parameters.isPositive("nm/slice", settings.getNmPerSlice());
			Parameters.isAbove("Amplitude fraction", settings.getAmplitudeFraction(), 0.01);
			Parameters.isBelow("Amplitude fraction", settings.getAmplitudeFraction(), 0.9);
			Parameters.isPositive("Start background frames", settings.getStartBackgroundFrames());
			Parameters.isPositive("End background frames", settings.getEndBackgroundFrames());
			Parameters.isAbove("Total background frames",
					settings.getStartBackgroundFrames() + settings.getEndBackgroundFrames(), 1);
			Parameters.isAbove("Magnification", settings.getMagnification(), 1);
			Parameters.isAbove("Smoothing", settings.getSmoothing(), 0);
			Parameters.isBelow("Smoothing", settings.getSmoothing(), 1);
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
		if (settings.getMode() == MODE_FITTING)
		{
			runUsingFitting();
		}
		else
		{
			try
			{
				runUsingAlignment();
			}
			catch (OutOfMemoryError e)
			{
				MemoryPeakResults.runGC();
				IJ.log(ExceptionUtils.getStackTrace(e));
				IJ.showMessage(TITLE, TextUtils.wrap("Out-of-memory. You may be using too many spots or " +
						"too large a PSF projection. The default projection size is 2.", 80));
			}
		}

		SettingsManager.writeSettings(settings);

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

		BasePoint[] spots = getSpots(0, true);
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
		LoessInterpolator loess = new LoessInterpolator(settings.getSmoothing(), iterations);

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
			double limit = smoothA[maximumIndex] * settings.getAmplitudeFraction();
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

		if (settings.getInteractiveMode())
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
		ImageStack psf = createStack(stack, minz, maxz, settings.getMagnification());

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
			ok = addToPSF(maxz, settings.getMagnification(), psf, centre, spot, regionBounds, progress, increment,
					settings.getCentreEachSlice());
		}

		if (settings.getInteractiveMode())
		{
			Utils.hide(TITLE_INTENSITY);
		}

		IJ.showProgress(1);

		if (!ok || stats.getN() == 0)
			return;

		final double avSd = getAverage(averageSd, averageA, 2);
		Utils.log("  Average background = %.2f, Av. SD = %s px", stats.getMean(), Utils.rounded(avSd, 4));

		normalise(psf, maxz, avSd * settings.getMagnification(), false);
		IJ.showProgress(1);

		psfImp = Utils.display(TITLE_PSF, psf);
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

		double[][] com = calculateCentreOfMass(psf, fitCom, nmPerPixel / settings.getMagnification());
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
		psfImp.setProperty("Info", ImagePSFHelper.toString(ImagePSFHelper.create(maxz,
				nmPerPixel / settings.getMagnification(), settings.getNmPerSlice(), stats.getN(), fwhm, createNote())));

		Utils.log("%s : z-centre = %d, nm/Pixel = %s, nm/Slice = %s, %d images, PSF SD = %s nm, FWHM = %s px\n",
				psfImp.getTitle(), maxz, Utils.rounded(nmPerPixel / settings.getMagnification(), 3),
				Utils.rounded(settings.getNmPerSlice(), 3), stats.getN(), Utils.rounded(fittedSd * nmPerPixel, 4),
				Utils.rounded(fwhm));

		createInteractivePlots(psf, maxz, nmPerPixel / settings.getMagnification(), fittedSd * nmPerPixel);

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
		if (settings.getInteractiveMode())
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
			plot.addLabel(0, 0,
					String.format("Amplitude = %s (%sx). z = %s nm", Utils.rounded(amplitude),
							Utils.rounded(amplitude / maxAmplitude),
							Utils.rounded((slice - zCentre) * settings.getNmPerSlice())));

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
			plot.addLabel(0, 0,
					String.format("Width = %s nm (%sx). z = %s nm", Utils.rounded(width * nmPerPixel),
							Utils.rounded(width * nmPerPixel / psfWidth),
							Utils.rounded((slice - zCentre) * settings.getNmPerSlice())));

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
		for (int i = 0; i < settings.getStartBackgroundFrames(); i++)
		{
			first.add(spot[i]);
		}
		for (int i = 0, j = spot.length - 1; i < settings.getEndBackgroundFrames(); i++, j--)
		{
			last.add(spot[j]);
		}
		float av = (float) ((first.getSum() + last.getSum()) / (first.getN() + last.getN()));
		Utils.log("  Spot %d Background: First %d = %.2f, Last %d = %.2f, av = %.2f", n,
				settings.getStartBackgroundFrames(), first.getMean(), settings.getEndBackgroundFrames(), last.getMean(),
				av);
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

		if (settings.getInteractiveMode())
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
					fp.setInterpolationMethod(settings.getInterpolationMethod());
					fp = (FloatProcessor) fp.resize(dstWidth, dstHeight);

					// In the case of Bicubic interpolation check for negative values
					if (settings.getInterpolationMethod() == ImageProcessor.BICUBIC)
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

	private double[] calculateCenterOfMass(FloatProcessor fp)
	{
		final int h = fp.getHeight();
		final int w = fp.getWidth();
		float[] data = (float[]) fp.getPixels();
		final double threshold = Maths.max(data) * settings.getComCutOff();
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
		// We have a different fit configuration just for the PSF Creator.
		// This allows it to be saved and not effect PeakFit settings.
		FitEngineConfiguration config = new FitEngineConfiguration(settings.getFitEngineSettings(),
				settings.getCalibration(), settings.getPsf());
		boolean save = false;
		if (!c.showDialog(config, save))
		{
			IJ.error(TITLE, "No fit configuration loaded");
			return false;
		}

		config = c.getFitEngineConfiguration();
		config.configureOutputUnits();
		config.setResidualsThreshold(1);
		fitConfig = config.getFitConfiguration();
		nmPerPixel = fitConfig.getCalibrationWriter().getNmPerPixel();
		if (settings.getRadius() < 5 * FastMath.max(fitConfig.getInitialXSD(), fitConfig.getInitialYSD()))
		{
			settings.setRadius(5 * FastMath.max(fitConfig.getInitialXSD(), fitConfig.getInitialYSD()));
			Utils.log("Radius is less than 5 * PSF standard deviation, increasing to %s",
					Utils.rounded(settings.getRadius()));
		}
		boxRadius = (int) Math.ceil(settings.getRadius());

		settings.setFitEngineSettings(config.getFitEngineSettings());
		settings.setCalibration(fitConfig.getCalibration());
		settings.setPsf(fitConfig.getPSF());
		SettingsManager.writeSettings(settings);

		return true;
	}

	/**
	 * @return Extract all the ROI points
	 */
	private BasePoint[] getSpots()
	{
		float z = imp.getStackSize() / 2;
		Roi roi = imp.getRoi();
		if (roi != null && roi.getType() == Roi.POINT)
		{
			FloatPolygon p = roi.getFloatPolygon();
			int n = p.npoints;

			float offset = 0.5f;

			// Check if already float coordinates
			if (!SimpleArrayUtils.isInteger(p.xpoints) || !SimpleArrayUtils.isInteger(p.ypoints))
				offset = 0;

			BasePoint[] roiPoints = new BasePoint[n];
			for (int i = 0; i < n; i++)
			{
				roiPoints[i] = new BasePoint(p.xpoints[i] + offset, p.ypoints[i] + offset, z);
			}
			return roiPoints;
		}
		return new BasePoint[0];
	}

	/**
	 * Extract all the ROI points and optionally exclude those that have a box region overlapping with any other spot.
	 *
	 * @param offset
	 *            the offset
	 * @return the spots
	 */
	private BasePoint[] getSpots(float offset, boolean checkOverlap)
	{
		float z = imp.getStackSize() / 2;
		//float z = (imp.getStackSize() - 1) / 2.0f; // Interpolate between slices
		Roi roi = imp.getRoi();
		if (roi != null && roi.getType() == Roi.POINT)
		{
			FloatPolygon p = roi.getFloatPolygon();
			int n = p.npoints;
			if (n == 0)
				return new BasePoint[0];

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

			return (checkOverlap) ? checkSpotOverlap(roiPoints) : roiPoints;
		}
		return new BasePoint[0];
	}

	/**
	 * Get all the ROI points that have a box region not overlapping with any other spot.
	 *
	 * @param offset
	 *            the offset
	 * @return the spots
	 */
	private BasePoint[] checkSpotOverlap(BasePoint[] roiPoints)
	{
		return getNonBadSpots(roiPoints, findSpotOverlap(roiPoints, null));
	}

	/**
	 * Find all the ROI points that have a box region overlapping with any other spot.
	 *
	 * @param offset
	 *            the offset
	 * @return the overlap array
	 */
	private boolean[] findSpotOverlap(BasePoint[] roiPoints, boolean[] excluded)
	{
		int n = roiPoints.length;
		boolean[] bad = new boolean[n];
		if (n == 1)
			return bad;

		// Check overlap of box regions
		int w = imp.getWidth();
		int h = imp.getHeight();
		ImageExtractor ie = new ImageExtractor(null, w, h);
		Rectangle[] regions = new Rectangle[n];
		// Check size if not fitting
		int size = (settings.getMode() != MODE_FITTING) ? 2 * boxRadius + 1 : Integer.MAX_VALUE;
		for (int i = 0; i < n; i++)
		{
			if (excluded != null && excluded[i])
				continue;
			Rectangle r = ie.getBoxRegionBounds(roiPoints[i].getXint(), roiPoints[i].getYint(), boxRadius);
			regions[i] = r;
			if (r.width < size || r.height < size)
			{
				Utils.log("Warning: Spot %d region extends beyond the image, border pixels will be duplicated", i + 1);
			}
		}

		// Add support for 3D overlap analysis. Only do this if a zRadius has been specified.
		boolean is3D = (settings.getMode() == MODE_ALIGNMENT && zRadius > 0);

		for (int i = 0; i < n; i++)
		{
			if (excluded != null && excluded[i])
				continue;
			if (bad[i]) // Already found to overlap
				continue;
			// Check intersect with others
			for (int j = i; ++j < n;)
			{
				if (excluded != null && excluded[j])
					continue;
				boolean overlap = regions[i].intersects(regions[j]);
				if (overlap && is3D)
				{
					// Reset (assume non-overlapping)
					overlap = false;

					// Check for 3D overlap:
					// iiiiiiiiiiiiiiiii
					//       jjjjjjjjjjjjjjjjjjj
					int mini = roiPoints[i].getZint() - zRadius;
					int maxi = roiPoints[i].getZint() + zRadius;
					int minj = roiPoints[j].getZint() - zRadius;
					int maxj = roiPoints[j].getZint() + zRadius;
					if (mini <= minj)
					{
						overlap = (maxi >= minj);
					}
					else // (minj< mini)
					{
						overlap = (maxj >= mini);
					}
				}
				if (overlap)
				{
					Utils.log("Warning: Spot %d region overlaps with spot %d, ignoring both", i + 1, j + 1);
					bad[i] = bad[j] = true;
					break;
				}
			}
		}

		return bad;
	}

	private BasePoint[] getNonBadSpots(BasePoint[] roiPoints, boolean[] bad)
	{
		int ok = 0;
		for (int i = 0, n = bad.length; i < n; i++)
		{
			if (bad[i])
				continue;
			roiPoints[ok++] = roiPoints[i];
		}
		return Arrays.copyOf(roiPoints, ok);
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
			if (settings.getInteractiveMode())
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
				b.setValue(b.getValue() * settings.getMagnification());
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
				settings.getMagnification() * 4);
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
		final double scale = 1.0 / settings.getMagnification();
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
				Utils.rounded((slice - zCentre) * settings.getNmPerSlice())));
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
				Utils.rounded(distanceThreshold), Utils.rounded((z - zCentre) * settings.getNmPerSlice())));
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

	private class Smoother
	{
		LoessInterpolator loess;
		double[] xval, yval;

		double[] dsmooth;
		float[] fsmooth;

		public Smoother()
		{
			if (settings.getSmoothing() != 0)
				loess = new LoessInterpolator(settings.getSmoothing(), 1);
		}

		public Smoother smooth(float[] data)
		{
			dsmooth = null;
			fsmooth = null;
			if (loess == null)
			{
				fsmooth = data;
				return this;
			}
			if (xval == null || xval.length != data.length)
				xval = SimpleArrayUtils.newArray(data.length, 0, 1.0);
			// Convert the data for smoothing
			if (yval == null || yval.length != data.length)
				yval = new double[data.length];
			for (int i = 0; i < data.length; i++)
				yval[i] = data[i];
			dsmooth = loess.smooth(xval, yval);
			return this;
		}

		public Smoother smooth(double[] data)
		{
			dsmooth = null;
			fsmooth = null;
			if (loess == null)
			{
				dsmooth = data;
				return this;
			}
			if (xval == null || xval.length != data.length)
				xval = SimpleArrayUtils.newArray(data.length, 0, 1.0);
			dsmooth = loess.smooth(xval, data);
			return this;
		}

		public float[] getFSmooth()
		{
			// Either fsmooth or dsmooth will not be null
			if (fsmooth == null)
				fsmooth = SimpleArrayUtils.toFloat(dsmooth);
			return fsmooth;
		}

		public double[] getDSmooth()
		{
			// Either fsmooth or dsmooth will not be null
			if (dsmooth == null)
				dsmooth = SimpleArrayUtils.toDouble(fsmooth);
			return dsmooth;
		}
	}

	private void runUsingAlignment()
	{
		if (!showAlignmentDialog())
			return;

		boxRadius = (int) Math.ceil(settings.getRadius());

		CalibrationReader calibration = new CalibrationReader(settings.getCalibration());

		// Limit this
		settings.setAnalysisWindow(Math.min(settings.getAnalysisWindow(), boxRadius / 2));

		// Find the selected PSF spots x,y,z centre
		// We offset the centre to the middle of pixel.
		BasePoint[] centres = getSpots(0.5f, false);
		if (centres.length == 0)
		{
			IJ.error(TITLE, "No PSFs");
			return;
		}

		CameraModel cameraModel = null;
		if (calibration.isSCMOS())
		{
			cameraModel = CameraModelManager.load(calibration.getCameraModelName());
			if (cameraModel == null)
			{
				IJ.error(TITLE, "No camera model");
				return;
			}
			cameraModel = PeakFit.cropCameraModel(cameraModel, imp.getWidth(), imp.getHeight(), 0, 0, true);
		}
		else
		{
			cameraModel = new FixedPixelCameraModel(calibration.getBias(), 1);
		}

		// Extract the image data for processing as float
		float[][] image = CreateData.extractImageStack(imp, 0, imp.getStackSize() - 1);

		for (float[] data : image)
			cameraModel.removeBiasAndGain(data);

		zSelector = new PSFCentreSelector();

		// Relocate the initial centres
		Utils.showStatus("Relocating initial centres");
		centres = relocateCentres(image, centres);
		if (centres == null)
			return;

		zRadius = (int) Math.ceil(settings.getAlignmentZRadius());

		// Check the region overlap in 3D and exclude overlapping PSFs
		boolean[] bad = findSpotOverlap(centres, null);
		centres = getNonBadSpots(centres, bad);
		if (centres.length == 0)
		{
			IJ.error(TITLE, "No PSFs without neighbours within the box region");
			return;
		}

		// Multi-thread for speed
		if (threadPool == null)
			threadPool = Executors.newFixedThreadPool(Prefs.getThreads());

		// Extract each PSF into a scaled PSF
		Utils.showStatus(String.format("[%d] Extracting PSFs", 0));
		ExtractedPSF[] psfs = extractPSFs(image, centres);

		Point location = null;

		// Iterate until centres have converged
		boolean converged = false;
		for (int iter = 0; !converged && iter < settings.getMaxIterations(); iter++)
		{
			// Combine all PSFs
			Utils.showStatus(String.format("[%d] Aligning PSFs", iter + 1));
			ExtractedPSF combined = combine(psfs);
			combined.createProjections();

			// Get the current combined z-centre. 
			// This is used to get the centre of mass for repositioning.
			// It also effects the alignment so do it for the first iteration.
			zSelector.setPSF(combined);
			if (iter == 0)
			{
				// TODO - check if the z-centre should be guessed here. 
				// We assume that the combined PSF may be easier to guess if the initial 
				// guess for each individual PSF was OK. It may not be necessary since all 
				// the PSFs are combined around their z-centres. Once alignment has 
				// started we skip this step.
				zSelector.analyse();
				zSelector.guessZCentre();
			}
			if (settings.getInteractiveMode())
			{
				if (iter != 0)
					zSelector.analyse();
				//zSelector.guessZCentre();
				double dz = zSelector.run("Update combined PSF z-centre", true, false, false, null);
				if (dz < 0)
					return;
			}

			// Align each to the combined PSF
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

			final boolean[] excluded = new boolean[psfs.length];
			if (checkAlignments)
			{
				combined.show(TITLE_PSF);

				// Ask about each centre in turn.
				// Update Point ROI using float coordinates and set image slice to 
				// correct z-centre.
				//imp.saveRoi();
				imp.killRoi();
				ImageCanvas ic = imp.getCanvas();
				//ic.setMagnification(16);
				int reject = 0;
				float box = boxRadius + 0.5f;
				int n = imp.getStackSize();
				for (int j = 0; j < centres.length; j++)
				{
					psfs[j].show(TITLE_SPOT_PSF);

					Overlay o = new Overlay();
					o.add(createRoi(psfs[j].centre.getX(), psfs[j].centre.getY(), Color.RED));
					float cx = centres[j].getX();
					float cy = centres[j].getY();
					o.add(createRoi(cx, cy, Color.GREEN));
					Roi roi = new Roi(cx - box, cy - box, 2 * box, 2 * box);
					o.add(roi);
					// The centre is absolute within the original stack
					imp.setSlice(Maths.clip(1, n, (int) Math.round(centres[j].getZ())));
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
					NonBlockingExtendedGenericDialog gd = new NonBlockingExtendedGenericDialog(TITLE);
					gd.addMessage(String.format("Shift X = %s\nShift Y = %s\nShift Z = %s",
							rounder.toString(translation[j][0]), rounder.toString(translation[j][1]),
							rounder.toString(translation[j][2])));
					final int spotIndex = j;
					gd.addAndGetButton("Exclude spot", new ActionListener()
					{
						public void actionPerformed(ActionEvent e)
						{
							if (excluded[spotIndex])
							{
								Utils.log("Included spot %d", spotIndex + 1);
								excluded[spotIndex] = false;
							}
							else
							{
								Utils.log("Excluded spot %d", spotIndex + 1);
								excluded[spotIndex] = true;
							}
						}
					});
					gd.enableYesNoCancel("Accept", "Reject");
					if (location != null)
						gd.setLocation(location.x, location.y);
					gd.showDialog();
					if (gd.wasCanceled())
					{
						resetImp();
						return;
					}
					boolean failed = excluded[spotIndex] || !gd.wasOKed();
					if (failed)
					{
						reject++;
						centres[j] = psfs[j].centre;
						Arrays.fill(translation[j], 0f); // For RMSD computation
					}
					location = gd.getLocation();
				}
				resetImp();
				if (reject == psfs.length)
				{
					IJ.error(TITLE, "No PSF translations were accepted");
					return;
				}
			}

			bad = findSpotOverlap(centres, excluded);
			int badCount = count(bad);
			int excludedCount = count(excluded);
			int ok = bad.length - badCount - excludedCount;
			if (ok < bad.length)
			{
				if (badCount != 0 && settings.getInteractiveMode())
				{
					ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
					gd.addMessage("Warning: Regions now overlap!");
					gd.addMessage("OK = " + TextUtils.pleural(ok, "PSF"));
					gd.addMessage("Overlapping = " + TextUtils.pleural(badCount, "PSF"));
					//gd.addMessage("Excluded = " + TextUtils.pleural(excludedCount, "PSF"));
					gd.enableYesNoCancel("Exclude", "Include");
					if (location != null)
						gd.setLocation(location.x, location.y);
					gd.showDialog();
					if (gd.wasCanceled())
					{
						resetImp();
						return;
					}
					if (!gd.wasOKed())
					{
						Arrays.fill(bad, false); // allow bad spots
						ok = bad.length;
					}
					location = gd.getLocation();
				}

				if (ok == 0)
				{
					IJ.error(TITLE, "No PSFs remaining");
					resetImp();
					return;
				}
			}

			// Merge bad and excluded to get new centres
			for (int i = 0; i < bad.length; i++)
				if (excluded[i])
					bad[i] = true;
			ok = bad.length - count(bad);

			BasePoint[] newCentres = getNonBadSpots(centres, bad);

			// Find the change in centres
			double[] rmsd = new double[2];
			for (int j = 0; j < psfs.length; j++)
			{
				if (bad[j])
					continue;
				rmsd[0] += Maths.pow2(translation[j][0]) + Maths.pow2(translation[j][1]);
				rmsd[1] += Maths.pow2(translation[j][2]);
			}
			for (int j = 0; j < 2; j++)
				rmsd[j] = Math.sqrt(rmsd[j] / ok);

			Utils.showStatus(String.format("[%d] Checking combined PSF", iter + 1));

			// Compute CoM shift using the current z-centre and z-window
			double[] shift = combined.getCentreOfMassXYShift(zSelector.getCentreSlice());
			double shiftd = Math.sqrt(shift[0] * shift[0] + shift[1] * shift[1]);

			Utils.log("[%d] RMSD XY = %s : RMSD Z = %s : Combined CoM shift = %s,%s (%s)", iter,
					rounder.toString(rmsd[0]), rounder.toString(rmsd[1]), rounder.toString(shift[0]),
					rounder.toString(shift[1]), rounder.toString(shiftd));

			if (settings.getInteractiveMode())
			{
				// Ask if OK to continue?
				GenericDialog gd = new GenericDialog(TITLE);
				gd.addMessage(String.format("RMSD XY = %s\nRMSD Z = %s\nCombined CoM shift = %s,%s (%s)",
						rounder.toString(rmsd[0]), rounder.toString(rmsd[1]), rounder.toString(shift[0]),
						rounder.toString(shift[1]), rounder.toString(shiftd)));
				// Check if we can do more iterations
				if (iter + 1 < settings.getMaxIterations())
					gd.enableYesNoCancel("Continue", "Converged");
				else
					gd.setOKLabel("Converged");
				gd.showDialog();
				if (gd.wasCanceled())
					return;
				converged = !gd.wasOKed();
			}
			else
			{
				// Sensible convergence on minimal shift
				converged = rmsd[0] < 0.01 && rmsd[1] < 0.05 && shiftd < 0.001;
			}

			// For the next round we move to the non-overlapping spots
			centres = newCentres;

			// Update the centres using the centre-of-mass of the combined PSF
			centres = updateUsingCentreOfMassXYShift(shift, shiftd, combined, centres);

			// Extract each PSF into a scaled PSF
			Utils.showStatus(String.format("[%d] Extracting PSFs", iter + 1));
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

		// Show an interactive dialog for cropping the PSF and choosing the 
		// final output
		PSFOutputSelector cropSelector = new PSFOutputSelector(combined);
		combined = cropSelector.run();
		if (combined == null)
			return;

		// For an image PSF we can just enlarge the PSF and window.
		// For a CSpline then we already have the 3D cubic spline function.
		// However we want to post-process the function to allow windowing and 
		// normalisation. So we enlarge by 3 in each dimension.
		// The CSpline can be created by solving the coefficients for the 
		// 4x4x4 (64) sampled points on each node. 

		int magnification;
		if (settings.getOutputType() == OUTPUT_TYPE_IMAGE_PSF)
		{
			magnification = settings.getPsfMagnification();
		}
		else
		{
			magnification = 3;
		}

		// Enlarge the combined PSF for final processing
		ExtractedPSF finalPSF = combined.enlarge(magnification, threadPool);

		// Show a dialog to collect final z-centre interactively
		Utils.showStatus("Analysing PSF");
		zSelector.setPSF(finalPSF);
		zSelector.analyse();
		//zSelector.guessZCentre(); // No need to guess the centre
		double dz = zSelector.run("Finalise PSF", true, true, true, null);
		if (dz < 0)
			return;

		zCentre = zSelector.getCentreSlice();

		if (settings.getCropToZCentre())
		{
			finalPSF = finalPSF.cropToZCentre(zCentre);
			zCentre = finalPSF.stackZCentre + 1; // Back to 1-based index
		}

		// When click ok the background is subtracted from the PSF
		// All pixels below the background are set to zero
		// Apply a Tukey window to roll-off to zero at the outer pixels

		Utils.showStatus("Windowing PSF");
		double[] wx = ImageWindow.tukeyEdge(finalPSF.maxx, settings.getWindow());
		double[] wz = ImageWindow.tukeyEdge(finalPSF.psf.length, settings.getWindow());

		// Normalisation so the max intensity frame is one
		float[][] psf = finalPSF.psf;
		int maxz = psf.length;
		double[] sum = new double[maxz];
		for (int z = 0; z < maxz; z++)
		{
			sum[z] = applyWindow(psf[z], z, wx, wz, zSelector.background);
		}

		// Smooth the intensity
		Utils.showStatus("Normalising PSF");
		Smoother smoother = zSelector.ssmoother;
		double[] ssum = smoother.smooth(sum).getDSmooth();

		// Compute normalisation and apply.
		SimpleArrayUtils.multiply(ssum, 1.0 / Maths.max(ssum));
		for (int z = 0; z < psf.length; z++)
		{
			if (sum[z] != 0)
				SimpleArrayUtils.multiply(psf[z], ssum[z] / sum[z]);
			sum[z] = Maths.sum(psf[z]);
		}

		// Show the final intensity profile
		double[] slice = SimpleArrayUtils.newArray(maxz, 1, 1.0);
		Plot plot = new Plot(TITLE_SIGNAL, "Slice", "Signal");
		double[] range = Maths.limits(sum);
		plot.setLimits(1, maxz, range[0], range[1]);
		plot.setColor(Color.black);
		plot.addPoints(slice, sum, Plot.LINE);
		Utils.display(TITLE_SIGNAL, plot);

		// Create a new extracted PSF and show
		Utils.showStatus("Displaying PSF");
		magnification = finalPSF.magnification;
		finalPSF = new ExtractedPSF(psf, finalPSF.maxx, finalPSF.centre, magnification);
		finalPSF.createProjections();
		psfOut = finalPSF.show(TITLE_PSF, zCentre);
		psfImp = psfOut[0];

		// Add image info
		int imageCount = centres.length;
		ImagePSF.Builder imagePsf = ImagePSFHelper.create(zCentre, nmPerPixel / magnification,
				settings.getNmPerSlice() / magnification, imageCount, 0, createNote()).toBuilder();

		// Add the CoM
		// Find the XY centre around the z centre 
		double[] com = getCentreOfMassXY(finalPSF.psf, finalPSF.maxx, finalPSF.maxy, zCentre - 1,
				settings.getComWindow(), getCoMXYBorder(finalPSF.maxx, finalPSF.maxy));

		imagePsf.setXCentre(com[0]);
		imagePsf.setYCentre(com[1]);
		imagePsf.setZCentre(zCentre - 1);
		psfImp.setProperty("Info", ImagePSFHelper.toString(imagePsf));

		psfImp.setRoi(new PointRoi(com[0], com[1]));
		psfImp.setSlice(zCentre);
		psfImp.resetDisplayRange();
		psfImp.updateAndDraw();

		Utils.log("Final Centre-of-mass = %s,%s\n", rounder.toString(com[0]), rounder.toString(com[1]));
		Utils.log("%s : z-centre = %d, nm/Pixel = %s, nm/Slice = %s, %d images\n", psfImp.getTitle(), zCentre,
				Utils.rounded(imagePsf.getPixelSize(), 3), Utils.rounded(imagePsf.getPixelDepth(), 3), imageCount);

		if (settings.getOutputType() == OUTPUT_TYPE_CSPLINE)
		{
			// Ask this again as it is important
			//if (TextUtils.isNullOrEmpty(settings.getSplineFilename()))
			//{
			final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
			gd.addFilenameField("Spline_filename", settings.getSplineFilename());
			gd.showDialog(true);
			if (gd.wasCanceled())
				return;
			settings.setSplineFilename(gd.getNextString());
			//}
			if (!TextUtils.isNullOrEmpty(settings.getSplineFilename()))
			{
				// Save the result ...
				IJ.showStatus("Creating cubic spline");
				CubicSplinePSF cubicSplinePSF = CubicSplineManager.createCubicSpline(imagePsf, psfImp.getImageStack(),
						settings.getSinglePrecision());

				IJ.showStatus("Saving cubic spline");
				CubicSplineManager.save(cubicSplinePSF, settings.getSplineFilename());

				IJ.showStatus("Spline saved to " + settings.getSplineFilename());
				return; // To leave the status message
			}
		}
		IJ.showStatus("");
	}

	private void resetImp()
	{
		imp.restoreRoi();
		imp.setOverlay(null);
		imp.getWindow().toFront();
	}

	private int count(boolean[] flags)
	{
		int c = 0;
		for (int i = 0; i < flags.length; i++)
			if (flags[i])
				c++;
		return c;
	}

	private int getCoMXYBorder(int maxx, int maxy)
	{
		int w = Math.min(maxx, maxy);
		return Maths.clip(0, w / 2 - 1, (int) Math.round(w * settings.getComBorder()));
	}

	/**
	 * Extract the stack for each centre and try and guess the z-centre based on the type of PSF. Relocate the XY centre
	 * using the centre-of-mass around the pixels close to the z-centre slice.
	 *
	 * @param image
	 *            the image
	 * @param centres
	 *            the centres
	 * @return the new centres
	 */
	private BasePoint[] relocateCentres(float[][] image, BasePoint[] centres)
	{
		int w = imp.getWidth();
		int h = imp.getHeight();

		// Just for getting the bounds
		ImageExtractor ie = new ImageExtractor(image[0], w, h);

		// This can be reused as a buffer
		float[][] psf = new float[image.length][];

		//double mean = 0;

		for (int i = 0; i < centres.length; i++)
		{
			// Extract stack
			int x = centres[i].getXint();
			int y = centres[i].getYint();
			Rectangle bounds = ie.getBoxRegionBounds(x, y, boxRadius);
			for (int z = 0; z < image.length; z++)
			{
				psf[z] = ImageConverter.getData(image[z], w, h, bounds, psf[z]);
			}

			if (settings.getInteractiveMode())
			{
				Overlay o = new Overlay();
				Roi roi = new Roi(bounds);
				o.add(roi);
				imp.setOverlay(o);
				imp.updateAndDraw();
			}

			// Create a PSF and select the z-centre
			zSelector.setPSF(new ExtractedPSF(psf, bounds.width, bounds.height));
			zSelector.analyse();
			zSelector.guessZCentre();

			if (settings.getInteractiveMode())
			{
				// Ask user for z-centre confirmation
				double dz = zSelector.run("Confirm PSF Z-centre", true, false, false, Integer.toString(i + 1));
				if (dz == -1)
				{
					resetImp();
					return null;
				}
				if (dz == -2)
				{
					// Exclude this PSF
					centres[i] = null;
					continue;
				}
			}

			zCentre = zSelector.getCentreSlice() - 1;

			// Subtract background
			final float adjust = -zSelector.background;
			for (int z = 0; z < image.length; z++)
			{
				SimpleArrayUtils.add(psf[z], adjust);
			}

			// Update centre
			double[] com = getCentreOfMassXY(psf, bounds.width, bounds.height, zCentre, settings.getComWindow(),
					getCoMXYBorder(bounds.width, bounds.height));
			float dx = (float) (com[0] + bounds.x - centres[i].getX());
			float dy = (float) (com[1] + bounds.y - centres[i].getY());
			float dz = (float) (zSelector.zCentre - centres[i].getZ());
			Utils.log("Centre %d : %s,%s,%s updated by %s,%s,%s", i + 1, rounder.toString(centres[i].getX()),
					rounder.toString(centres[i].getY()), rounder.toString(centres[i].getZ()), rounder.toString(dx),
					rounder.toString(dy), rounder.toString(dz));
			centres[i] = centres[i].shift(dx, dy, dz);
			//mean += centres[i].getZ();
		}

		if (settings.getInteractiveMode())
		{
			imp.setOverlay(null);

			// Check if any centres were excluded
			int size = 0;
			for (int i = 0; i < centres.length; i++)
			{
				if (centres[i] == null)
					continue;
				centres[size++] = centres[i];
			}
			if (size == 0)
			{
				resetImp();
				IJ.error(TITLE, "No remaining PSF centres");
				return null;
			}
			if (size < centres.length)
				centres = Arrays.copyOf(centres, size);
		}

		//// z-centres should be relative to the combined stack, not absolute
		//mean /= centres.length;
		//for (int i = 0; i < centres.length; i++)
		//	centres[i] = new BasePoint(centres[i].getX(), centres[i].getY(), (float) (centres[i].getZ() - mean));

		return centres;
	}

	private static final TypeConverter<AngleUnit> angleConverter = UnitConverterFactory
			.createConverter(AngleUnit.RADIAN, AngleUnit.DEGREE);
	private static final double NO_ANGLE = -360.0;

	private class PSFCentreSelector implements DialogListener
	{
		// The concept of HWHM only applies to a PSF that is a peaked maxima.
		// This may not be true for an image. To approximate this we assume that
		// the peak is Gaussian and find the sum which equals the integral of 
		// a Gaussian at HWHM = SD * 1.17741 (i.e. Gaussian2DFunction.SD_TO_FWHM_FACTOR)
		final double integral = Erf.erf(Gaussian2DFunction.SD_TO_HWHM_FACTOR / Math.sqrt(2));

		ExtractedPSF psf;
		double zCentre;
		Point location = null;

		Smoother bsmoother = new Smoother();
		Smoother fsmoother = new Smoother();
		Smoother ssmoother = new Smoother();
		Smoother w0smoother = new Smoother();
		Smoother w1smoother = new Smoother();
		Smoother asmoother = new Smoother();
		float[][] limits;

		float[] bdata;
		float[] fdata;
		float[] sdata;
		double[] w0, w1, w01;
		double[] adata, asdata;
		int bIndex, fIndex, sIndex, wIndex, aIndex;
		float background;
		private Label backgroundLabel = null;
		boolean hasId, plotBackground, plotEdgeWindow, cropOption;

		// We choose the angle with the first spot
		double targetAngle = NO_ANGLE;

		public PSFCentreSelector()
		{
		}

		public int getCentreSlice()
		{
			return Maths.clip(1, psf.psf.length, (int) Math.round(1 + zCentre));
		}

		public int getAnalysisWindow()
		{
			return (int) Math.round(settings.getAnalysisWindow() * psf.magnification);
		}

		public void setPSF(ExtractedPSF psf)
		{
			this.psf = psf;
			zCentre = psf.stackZCentre;

			// Reset
			bIndex = -1;
			background = Float.NaN;
			sdata = null;
		}

		public void analyse()
		{
			// These are window dependent
			int window = getAnalysisWindow();
			limits = getLimits(psf.psf, psf.maxx, psf.maxy, window);

			bdata = bsmoother.smooth(limits[0]).getFSmooth();
			fdata = fsmoother.smooth(limits[1]).getFSmooth();

			bIndex = SimpleArrayUtils.findMinIndex(bdata);
			fIndex = SimpleArrayUtils.findMaxIndex(fdata);
			background = bdata[bIndex];

			if (sdata == null)
			{
				int maxx = psf.maxx;
				int maxy = psf.maxy;
				int maxz = psf.psf.length;

				// These can be computed once per PSF
				sdata = new float[maxz];
				for (int z = 0; z < maxz; z++)
					sdata[z] = (float) Maths.sum(psf.psf[z]);
				ssmoother.smooth(sdata);
				sIndex = SimpleArrayUtils.findMaxIndex(ssmoother.getDSmooth());

				if (settings.getPsfType() == PSF_TYPE_DH)
				{
					// DoubleHelix - get rotation of moment of inertia.
					// Since we are interested in seeing the change in angle we track that
					// to avoid problems with wrapping the angle in the circle.
					adata = new double[maxz];
					double lastA = 180;
					double origin = 0;
					for (int z = 0; z < maxz; z++)
					{
						Tensor2D t = new Tensor2D(psf.psf[z], maxx, maxy);
						double[][] v = t.getEigenVectors();
						// Q. Which vector to use, small or big Eigen value?
						// Since they are orthogonal it does not matter.
						// Use arc tan to get the result in the domain -pi to pi
						double a = angleConverter.convert(Math.atan2(v[1][1], v[1][0]));
						double d = a - lastA;
						d += (d > 180.0) ? -360.0 : (d < -180.0) ? 360.0 : 0;

						// We can adjust to the opposite direction to get closer
						if (d > 90.0)
						{
							d -= 180.0;
							a -= 180.0;
						}
						else if (d < -90.0)
						{
							d += 180.0;
							a += 180.0;
						}

						lastA = a;
						adata[z] = d;

						if (z == 0)
						{
							origin = lastA;
							adata[z] = 0;
						}
					}
					// Reconstruct the angle from the deltas. 
					// The angle now may be outside the domain -180 to 180
					adata[0] = origin;
					for (int z = 1; z < maxz; z++)
					{
						adata[z] += adata[z - 1];
					}
					asdata = asmoother.smooth(adata).getDSmooth().clone();
					// Convert smoothed to the domain 0 to 180 for finding the centre
					// (We do not care about the direction so we discard the full 0 - 360 domain) 
					for (int z = 0; z < maxz; z++)
					{
						asdata[z] = asdata[z] % 180;
						if (asdata[z] < 0)
							asdata[z] += 180;
					}
					if (targetAngle == NO_ANGLE)
					{
						// Closest to the index which is our best guess at the z-centre
						// Use total signal for now. Assumes the helix loses light as it 
						// moves out of focus.
						double[] data = asdata;
						int best = data.length;
						for (int angle = 0; angle < 180; angle += 15)
						{
							int i = findIndex(data, sIndex, angle);
							int d = Math.abs(i - sIndex);
							if (d < best)
							{
								aIndex = i;
								best = d;
								targetAngle = angle;
							}
						}
					}
					else
					{
						aIndex = findIndex(asdata, sIndex, targetAngle);
					}
				}
				else
				{
					// TODO - check if this is working
					// The PSF should be background subtracted.
					// Perhaps just get the profile, find the centre
					// and the max and then move outward until at half-max
					// If there is more than one peak then the width
					// is invalid.

					// Spot - get FWHM of the X and Y sum projections
					w0 = new double[maxz];
					w1 = new double[maxz];
					double[] xp = new double[maxx];
					double[] yp = new double[maxy];
					for (int z = 0; z < maxz; z++)
					{
						float[] data = psf.psf[z];
						// Ensure the sum is positive
						float min = Maths.min(data);
						double min_by_x = min * maxx;
						double min_by_y = min * maxy;
						// rolling sums for each column/row
						double sr = 0, sc = 0;
						for (int x = 0; x < maxx; x++)
						{
							for (int y = 0, i = x; y < maxy; y++, i += maxx)
								sc += data[i];
							sc -= min_by_y;
							xp[x] = sc;
						}
						for (int y = 0, i = 0; y < maxy; y++)
						{
							for (int x = 0; x < maxx; x++, i++)
								sr += data[i];
							sr -= min_by_x;
							yp[y] = sr;
						}
						// Find centre
						w0[z] = hwhm(yp);
						w1[z] = hwhm(xp);
					}
					w0smoother.smooth(w0);
					w1smoother.smooth(w1);

					// Combine
					w01 = new double[maxz];
					for (int z = 0; z < maxz; z++)
					{
						w01[z] = Math.sqrt(w0[z] * w1[z]);
					}
					w01 = new Smoother().smooth(w01).getDSmooth();

					if (settings.getPsfType() == PSF_TYPE_ASTIGMATISM)
					{
						// Find half-way between min of each
						wIndex = (SimpleArrayUtils.findMinIndex(w0smoother.getDSmooth()) +
								SimpleArrayUtils.findMinIndex(w1smoother.getDSmooth())) / 2;
					}
					else
					{
						// Find min combined
						wIndex = SimpleArrayUtils.findMinIndex(w01);
					}
				}
			}
		}

		private int findIndex(double[] data, int start, double target)
		{
			// Slide along looking for minimim to the target
			int min = 0;
			int minD = data.length;

			// Sliding window of 3
			double d0;
			double d1 = Math.abs(data[0] - target);
			double d2 = Math.abs(data[1] - target);
			for (int i = 2; i < data.length; i++)
			{
				d0 = d1;
				d1 = d2;
				d2 = Math.abs(data[i] - target);
				// Check for a minimum
				if (d1 < d0 && d1 < d2)
				{
					int d = Math.abs(start - (i - 1));
					if (d < minD)
					{
						min = i - 1;
						minD = d;
					}
				}
			}
			return min;
		}

		private double hwhm(double[] xp)
		{
			int upper = xp.length - 1;
			int lx = Arrays.binarySearch(xp, xp[upper] / 2);
			if (lx < 0)
				lx = -(lx + 1);
			int ux = lx + 1;
			if (ux > upper)
				ux = upper;
			double target = integral * xp[upper];
			double s = xp[ux] - xp[lx];
			double lastS = s;
			// Work outwards from the centre
			while (s < target)
			{
				lastS = s;
				lx--;
				if (lx < 0)
					lx = 0;
				ux++;
				if (ux > upper)
					ux = upper;
				s = xp[ux] - xp[lx];
			}
			if (lastS != s)
			{
				double fraction = (target - s) / (lastS - s);
				return ((ux - lx) / 2.0 - fraction);
			}
			return (ux - lx) / 2.0;
		}

		public void guessZCentre()
		{
			if (settings.getPsfType() == PSF_TYPE_DH)
			{
				// DoubleHelix
				// Use foreground
				//zCentre = fIndex;
				// Use angle
				zCentre = aIndex;
			}
			else if (settings.getPsfType() == PSF_TYPE_ASTIGMATISM)
			{
				// Use closest to min width
				zCentre = wIndex;

				//// Use intersection between widths ...
				//// Use the smoothed data
				//double[] w0 = w0smoother.getDSmooth();
				//double[] w1 = w1smoother.getDSmooth();
				//
				//double mind = w0.length;
				//
				//for (int i = 1; i < w0.length; i++)
				//{
				//	// http://en.wikipedia.org/wiki/Line-line_intersection
				//	//
				//	//     x1,y1            x4,y4      
				//	//         **        ++ 
				//	//           **    ++
				//	//             **++ P(x,y)
				//	//            ++ **
				//	//          ++     **
				//	//        ++         **
				//	//    x3,y3            ** 
				//	//                       x2,y2
				//
				//	final double y1 = w0[i - 1];
				//	final double y2 = w0[i];
				//	final double y3 = w1[i - 1];
				//	final double y4 = w1[i];
				//
				//	// Check if they cross
				//	if (!((y3 >= y1 && y4 < y2) || (y1 >= y3 && y2 < y4)))
				//		continue;
				//
				//	final double x1 = i - 1;
				//	final double x2 = i;
				//	final double x3 = x1;
				//	final double x4 = x2;
				//
				//	final double x1_x2 = -1.0; //x1 - x2;
				//	final double x3_x4 = -1.0; //x3 - x4;
				//	final double y1_y2 = y1 - y2;
				//	final double y3_y4 = y3 - y4;
				//
				//	// Check if lines are parallel
				//	if (x1_x2 * y3_y4 - y1_y2 * x3_x4 == 0)
				//	{
				//		if (y1 == y3)
				//		{
				//			double d = Math.abs(x1 - wIndex);
				//			if (mind > d)
				//			{
				//				mind = d;
				//				zCentre = x1;
				//			}
				//		}
				//	}
				//	else
				//	{
				//		// Find intersection
				//		double px = ((x1 * y2 - y1 * x2) * x3_x4 - x1_x2 * (x3 * y4 - y3 * x4)) /
				//				(x1_x2 * y3_y4 - y1_y2 * x3_x4);
				//
				//		// Check if the intersection is within the two points
				//		// Q. Is this necessary given the intersection check above?
				//		if (px >= x1 && px < x2)
				//		{
				//			double d = Math.abs(px - wIndex);
				//			if (mind > d)
				//			{
				//				mind = d;
				//				zCentre = px;
				//			}
				//		}
				//	}
				//}
				//
				//zCentre = (int) Math.round(zCentre);
			}
			else
			{
				// Use foreground
				zCentre = fIndex;
				// Use min width
				// Note: The spot can move to a double ring Airy pattern which has no 
				// width so this doesn't work if the depth-of-field is high
				//zCentre = wIndex;
			}
		}

		public double run(String title, boolean plotBackground, boolean plotEdgeWindow, boolean cropOption, String id)
		{
			hasId = !TextUtils.isNullOrEmpty(id);
			this.plotBackground = plotBackground;
			this.plotEdgeWindow = plotEdgeWindow;
			this.cropOption = cropOption;

			int slice = getCentreSlice();

			psf.createProjections();
			psfOut = psf.show(TITLE_PSF, slice);

			// Show a dialog to collect processing options and 
			// find background interactively
			NonBlockingExtendedGenericDialog gd = new NonBlockingExtendedGenericDialog(TITLE);
			gd.addMessage(title);
			gd.addMessage("Current Z-centre = " + Utils.rounded(1 + zCentre));
			String label = "z-centre";
			if (hasId)
				label += "_" + id;
			gd.addSlider(label, 1, psf.psf.length, 1 + zCentre);
			final TextField tf = gd.getLastTextField();
			final String defaultZ = tf.getText();
			gd.addAndGetButton("Reset", new ActionListener()
			{
				public void actionPerformed(ActionEvent e)
				{
					tf.setText(defaultZ);
				}
			});
			if (hasId)
			{
				gd.addSlider("Z_radius (px)", 0, imp.getStackSize(), settings.getAlignmentZRadius());
			}
			gd.addSlider("CoM_z_window", 0, 8, settings.getComWindow());
			gd.addSlider("CoM_border", 0, 0.5, settings.getComBorder());
			if (plotBackground)
			{
				String label2 = "Analysis_window";
				if (psf.magnification != 1)
					label2 += "_mag_x" + psf.magnification;
				gd.addSliderIncludeDefault(label2, 0, Math.max(3, psf.maxx / 4),
						settings.getAnalysisWindow() * psf.magnification);
			}
			if (plotEdgeWindow)
			{
				gd.addSlider("Edge_window", 0, psf.maxx / 2, settings.getWindow());
			}
			if (cropOption)
			{
				gd.addCheckbox("Crop_to_z-centre", settings.getCropToZCentre());
			}
			gd.addDialogListener(this);
			if (Utils.isShowGenericDialog())
			{
				gd.addMessage("");
				backgroundLabel = gd.getLastLabel();
				drawPSFPlots();
			}
			gd.setLocation(location);
			if (hasId)
			{
				// This is when first selecting spots
				gd.enableYesNoCancel("Include", "Exclude");
			}
			gd.showDialog(true);

			// Remove interactive guides
			if (Utils.isShowGenericDialog())
			{
				for (int i = 0; i < 4; i++)
					psfOut[i].killRoi();
				psfOut[0].setOverlay(null);
				if (zRadius != 0)
				{
					psfOut[1].setOverlay(null);
					psfOut[2].setOverlay(null);
				}
			}

			removeCentreOnPlots();

			location = gd.getLocation();

			if (gd.wasCanceled())
				return -1;
			if (!gd.wasOKed())
				return -2;

			return zCentre;
		}

		private void drawPSFPlots()
		{
			WindowOrganiser wo = new WindowOrganiser();
			if (plotEdgeWindow)
				drawEdgeWindowPlot(wo);
			if (plotBackground)
				drawIntensityPlot(true, wo);
			wo.tile();
			drawPSFCentre();
			drawCentreOnPlots();
			drawCoMBorder();
		}

		public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
		{
			zCentre = gd.getNextNumber() - 1;
			if (hasId)
				settings.setAlignmentZRadius(gd.getNextNumber());
			settings.setComWindow((int) gd.getNextNumber());
			settings.setComBorder(gd.getNextNumber());
			if (plotBackground)
				settings.setAnalysisWindow(gd.getNextNumber() / psf.magnification);
			if (plotEdgeWindow)
				settings.setWindow((int) gd.getNextNumber());
			if (cropOption)
				settings.setCropToZCentre(gd.getNextBoolean());

			updatePSFPlots();
			return true;
		}

		private void updatePSFPlots()
		{
			if (plotEdgeWindow)
				updateEdgeWindowPlot();
			if (plotBackground)
				updateIntensityPlot();
			updatePSFCentre();
			updateCoMBorder();
		}

		int plotWindow = -1;

		private void drawEdgeWindowPlot(WindowOrganiser wo)
		{
			plotWindow = settings.getWindow();

			int size = psf.maxx;
			double[] x = SimpleArrayUtils.newArray(size, 0, 1.0);
			// Convert to an alpha
			double[] w = ImageWindow.tukeyEdge(size, plotWindow);

			Plot2 plot = new Plot2(TITLE_WINDOW, "x", "Weight", x, w);
			plot.setLimits(0, size - 1, 0, 1.05);
			PlotWindow pw = Utils.display(TITLE_WINDOW, plot);
			if (wo != null && Utils.isNewWindow())
				wo.add(pw);

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

		private void updateEdgeWindowPlot()
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
							while (plotWindow != settings.getWindow())
							{
								drawEdgeWindowPlot(null);
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

		private int plotBackgroundWindow = -1;
		private PlotWindow pwBackground, pwForeground, pwSignal, pwWidth, pwAngle;
		private float[] rangeB, rangeF, rangeS, rangeW, rangeA;

		private void drawIntensityPlot(boolean newData, WindowOrganiser wo)
		{
			plotBackgroundWindow = getAnalysisWindow();

			int length = psf.psf.length;
			double[] slice = SimpleArrayUtils.newArray(length, 1, 1.0);

			Plot plot = new Plot(TITLE_BACKGROUND, "Slice", "Background");
			rangeB = Maths.limits(limits[0]);
			plot.setLimits(1, length, rangeB[0], rangeB[1]);
			plot.setColor(Color.blue);
			plot.addPoints(slice, bsmoother.getDSmooth(), Plot.LINE);
			plot.setColor(Color.black);
			plot.addPoints(slice, SimpleArrayUtils.toDouble(limits[0]), Plot.LINE);
			String msgB = "Background = " + Utils.rounded(background);
			plot.addLabel(0, 0, msgB);
			pwBackground = Utils.display(TITLE_BACKGROUND, plot);
			if (wo != null && Utils.isNewWindow())
				wo.add(pwBackground);

			plot = new Plot(TITLE_FOREGROUND, "Slice", "Foreground");
			rangeF = Maths.limits(limits[1]);
			plot.setLimits(1, length, rangeF[0], rangeF[1]);
			plot.setColor(Color.blue);
			plot.addPoints(slice, fsmoother.getDSmooth(), Plot.LINE);
			plot.setColor(Color.black);
			plot.addPoints(slice, SimpleArrayUtils.toDouble(limits[1]), Plot.LINE);
			String msgF = "Foreground = " + Utils.rounded(fdata[fIndex]);
			plot.addLabel(0, 0, msgF);
			pwForeground = Utils.display(TITLE_FOREGROUND, plot);
			if (wo != null && Utils.isNewWindow())
				wo.add(pwForeground);

			if (newData)
			{
				plot = new Plot(TITLE_SIGNAL, "Slice", "Signal");
				rangeS = Maths.limits(sdata);
				plot.setLimits(1, length, rangeS[0], rangeS[1]);
				plot.setColor(Color.blue);
				plot.addPoints(slice, ssmoother.getDSmooth(), Plot.LINE);
				plot.setColor(Color.black);
				plot.addPoints(slice, SimpleArrayUtils.toDouble(sdata), Plot.LINE);
				String msgS = "Signal = " + Utils.rounded(sdata[sIndex]);
				plot.addLabel(0, 0, msgS);
				pwSignal = Utils.display(TITLE_SIGNAL, plot);
				if (wo != null && Utils.isNewWindow())
					wo.add(pwSignal);

				if (settings.getPsfType() == PSF_TYPE_DH)
				{
					// Double-Helix
					plot = new Plot(TITLE_ANGLE, "Slice", "Angle");
					double[] range = Maths.limits(adata);
					rangeA = SimpleArrayUtils.toFloat(range);
					plot.setLimits(1, length, range[0], range[1]);
					plot.setColor(Color.blue);
					plot.addPoints(slice, asmoother.getDSmooth(), Plot.LINE);
					plot.setColor(Color.black);
					plot.addPoints(slice, adata, Plot.LINE);
					String msgA = "Angle = " + Utils.rounded(asdata[aIndex]);
					plot.addLabel(0, 0, msgA);
					pwAngle = Utils.display(TITLE_ANGLE, plot);
					if (wo != null && Utils.isNewWindow())
						wo.add(pwAngle);
				}
				else
				{
					// Spot
					plot = new Plot(TITLE_HWHM, "Slice", "Width");
					double max = Maths.maxDefault(Maths.max(w0), w1);
					rangeW = new float[] { 0, (float) max };
					plot.setLimits(1, length, rangeW[0], rangeW[1]);
					plot.setColor(Color.blue);
					plot.addPoints(slice, w0smoother.getDSmooth(), Plot.LINE);
					plot.setColor(Color.blue.darker());
					plot.addPoints(slice, w0, Plot.LINE);
					plot.setColor(Color.red);
					plot.addPoints(slice, w1smoother.getDSmooth(), Plot.LINE);
					plot.setColor(Color.red.darker());
					plot.addPoints(slice, w1, Plot.LINE);
					plot.setColor(Color.magenta);
					plot.addPoints(slice, w01, Plot.LINE);
					plot.setColor(Color.black);
					String msgW = "Width = " + Utils.rounded(w01[wIndex]);
					plot.addLabel(0, 0, msgW);
					pwWidth = Utils.display(TITLE_HWHM, plot);
					if (wo != null && Utils.isNewWindow())
						wo.add(pwWidth);
				}
			}

			drawCentreOnPlots();

			// Draw a region on the image which contains the background.
			// This is not really useful as it is not linked to the current frame.
			//int cx = bIndex / psf.maxx;
			//int cy = bIndex % psf.maxx;
			//int y = cx - plotBackgroundWindow;
			//int x = cy - plotBackgroundWindow;
			//int w = 2 * plotBackgroundWindow + 1;
			//psfOut[0].setOverlay(new Roi(x, y, w, w), Color.YELLOW, 1, null);

			if (backgroundLabel != null)
			{
				backgroundLabel.setText(msgB);
			}
		}

		private void drawCentreOnPlots()
		{
			drawCentreOnPlot(pwBackground, rangeB);
			drawCentreOnPlot(pwForeground, rangeF);
			drawCentreOnPlot(pwSignal, rangeS);
			drawCentreOnPlot(pwWidth, rangeW);
			drawCentreOnPlot(pwAngle, rangeA);
		}

		private void drawCentreOnPlot(PlotWindow pw, float[] rangeB)
		{
			if (pw != null)
			{
				Plot plot = pw.getPlot();
				int slice = getCentreSlice();

				double x = plot.scaleXtoPxl(slice);
				double min = plot.scaleYtoPxl(rangeB[0]);
				double max = plot.scaleYtoPxl(rangeB[1]);

				pw.getImagePlus().setRoi(new Line(x, min, x, max));

				if (zRadius != 0)
				{
					double lx = plot.scaleXtoPxl(slice - zRadius);
					double ux = plot.scaleXtoPxl(slice + zRadius);
					Overlay o = new Overlay();
					Line l = new Line(lx, min, lx, max);
					l.setStrokeColor(Color.red);
					o.add(l);
					l = new Line(ux, min, ux, max);
					l.setStrokeColor(Color.red);
					o.add(l);
					pw.getImagePlus().setOverlay(o);
				}
				else
				{
					if (pw.getImagePlus().getOverlay() != null)
					{
						pw.getImagePlus().setOverlay(null);
					}
				}
			}
		}

		private void removeCentreOnPlots()
		{
			removeCentreOnPlot(pwBackground);
			removeCentreOnPlot(pwForeground);
			removeCentreOnPlot(pwSignal);
			removeCentreOnPlot(pwWidth);
			removeCentreOnPlot(pwAngle);
		}

		private void removeCentreOnPlot(PlotWindow pw)
		{
			if (pw != null)
				pw.getImagePlus().killRoi();
		}

		private void updateIntensityPlot()
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
							while (plotBackgroundWindow != getAnalysisWindow())
							{
								analyse();
								drawIntensityPlot(false, null);
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

		private int psfZCentre = -1, zRadius = 0;

		private int getZRadius()
		{
			if (hasId)
			{
				return (int) settings.getAlignmentZRadius();
				//return (int) Math.max(settings.getRadius(), settings.getAlignmentZRadius());
			}
			return 0;
		}

		private void drawPSFCentre()
		{
			psfZCentre = getCentreSlice();
			zRadius = getZRadius();

			// Select the z-centre
			psfOut[0].setSlice(psfZCentre);
			psfOut[0].resetDisplayRange();
			psfOut[0].updateAndDraw();

			// If selecting the original PSFs we can scroll the original image
			if (hasId)
			{
				imp.setSlice(psfZCentre);
				imp.setDisplayRange(psfOut[0].getDisplayRangeMin(), psfOut[0].getDisplayRangeMax());
				imp.updateAndDraw();
			}

			// Show border
			int border = getCoMXYBorder(psf.maxx, psf.maxy);
			psfOut[0].setRoi(border, border, psf.maxx - 2 * border, psf.maxy - 2 * border);

			// Mark XY projections
			// X-projection (dimension X is Z)
			psfOut[1].setRoi(new Line(psfZCentre, 0, psfZCentre, psfOut[1].getHeight()));
			// Y-projection (dimension Y is Z)
			psfOut[2].setRoi(new Line(0, psfZCentre, psfOut[2].getWidth(), psfZCentre));

			if (zRadius != 0)
			{
				int lz = psfZCentre - zRadius;
				int uz = psfZCentre + zRadius;

				// X-projection (dimension X is Z)
				psfOut[1].setOverlay(new Roi(lz, 0, uz - lz, psfOut[1].getHeight()), Color.red, 1, null);
				// Y-projection (dimension Y is Z)
				psfOut[2].setOverlay(new Roi(0, lz, psfOut[1].getWidth(), uz - lz), Color.red, 1, null);
			}
			else
			{
				if (psfOut[1].getOverlay() != null)
				{
					psfOut[1].setOverlay(null);
					psfOut[2].setOverlay(null);
				}
			}
		}

		private void updatePSFCentre()
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
							while (psfZCentre != getCentreSlice() || zRadius != getZRadius())
							{
								drawPSFCentre();
								drawCentreOnPlots();
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

		private int psfXYBorder = -1;

		private void drawCoMBorder()
		{
			// Show border
			psfXYBorder = getCoMXYBorder(psf.maxx, psf.maxy);
			psfOut[0].setRoi(psfXYBorder, psfXYBorder, psf.maxx - 2 * psfXYBorder, psf.maxy - 2 * psfXYBorder);
			psfOut[Projection.Z + 1].setRoi(psfOut[0].getRoi());
		}

		private void updateCoMBorder()
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
							while (psfXYBorder != getCoMXYBorder(psf.maxx, psf.maxy))
							{
								drawCoMBorder();
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
	}

	private class PSFOutputSelector implements DialogListener
	{
		ExtractedPSF psf;
		Label label1;
		int slice;

		public PSFOutputSelector(ExtractedPSF psf)
		{
			this.psf = psf;
		}

		public ExtractedPSF run()
		{
			// Show a dialog to collect processing options and 
			// find background interactively
			final NonBlockingExtendedGenericDialog gd = new NonBlockingExtendedGenericDialog(TITLE);
			gd.addMessage("Crop the final PSF");
			int maxz = psf.psf.length;
			gd.addSlider("Slice", 1, maxz, maxz / 2);
			// Take away 1 from the limits to avoid having dimensions 1 on an axis
			gd.addSlider("Crop_border", 0, Math.min(psf.maxx, psf.maxy) / 2 - 1, settings.getCropBorder());
			gd.addSlider("Crop_start", 0, maxz / 2 - 1, settings.getCropStart());
			gd.addSlider("Crop_end", 0, maxz / 2 - 1, settings.getCropEnd());
			gd.addMessage("Click ... to configure the output options");
			gd.addChoice("Output_type", OUTPUT_TYPE, settings.getOutputType(), new OptionListener<Integer>()
			{
				public boolean collectOptions(Integer value)
				{
					settings.setOutputType(value);
					boolean result = collectOptions(false);
					return result;
				}

				public boolean collectOptions()
				{
					return collectOptions(true);
				}

				private boolean collectOptions(boolean silent)
				{
					int outputType = settings.getOutputType();
					ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE);
					if (outputType == OUTPUT_TYPE_IMAGE_PSF)
					{
						// Image PSF
						egd.addSlider("PSF_magnification", 1, 8, settings.getPsfMagnification());
					}
					else
					{
						// CSpline
						egd.addCheckbox("Single_precision", settings.getSinglePrecision());
						egd.addFilenameField("Spline_filename", settings.getSplineFilename());
					}
					egd.setSilent(silent);
					egd.showDialog(true, gd);
					if (egd.wasCanceled())
						return false;
					if (outputType == OUTPUT_TYPE_IMAGE_PSF)
					{
						// Image PSF
						settings.setPsfMagnification((int) egd.getNextNumber());
					}
					else
					{
						// CSpline
						settings.setSinglePrecision(egd.getNextBoolean());
						settings.setSplineFilename(egd.getNextString());
					}
					return true;
				}
			});

			gd.addDialogListener(this);
			if (Utils.isShowGenericDialog())
			{
				gd.addOptionCollectedListener(new OptionCollectedListener()
				{
					public void optionCollected(OptionCollectedEvent e)
					{
						update();
					}
				});
				gd.addMessage("");
				label1 = gd.getLastLabel();

				// Get X and Y projections and use those to show crop start and end
				// with a line ROI.
				// TODO - Maybe show an intensity profile too ...

				psf.createProjections();
				psfOut = psf.show(TITLE_PSF);
				draw();
			}
			gd.showDialog(true);

			// Remove interactive guides
			if (Utils.isShowGenericDialog())
			{
				for (int i = 0; i < 4; i++)
					psfOut[i].killRoi();
				//psfOut[0].setOverlay(null);
				psfOut[1].setOverlay(null);
				psfOut[2].setOverlay(null);
			}

			if (gd.wasCanceled())
				return null;

			gd.collectOptions();

			int start = settings.getCropStart();
			int end = maxz - settings.getCropEnd();
			int size = end - start;
			if (size == maxz && settings.getCropBorder() == 0)
				return psf;

			int maxx = psf.maxx;
			int maxy = psf.maxy;

			Rectangle bounds = new Rectangle(settings.getCropBorder(), settings.getCropBorder(),
					maxx - 2 * settings.getCropBorder(), maxy - 2 * settings.getCropBorder());
			if (size < 2 || bounds.width < 2 || bounds.height < 2)
			{
				IJ.error(TITLE, "Invalid crop");
				return null;
			}

			float[][] psf2 = new float[size][];
			for (int z = start, i = 0; z < end; z++, i++)
			{
				ImageExtractor ie = new ImageExtractor(psf.psf[z], maxx, maxy);
				psf2[i] = ie.crop(bounds);
			}

			return new ExtractedPSF(psf2, bounds.width, bounds.height);
		}

		public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
		{
			slice = (int) gd.getNextNumber();
			settings.setCropBorder((int) gd.getNextNumber());
			settings.setCropStart((int) gd.getNextNumber());
			settings.setCropEnd((int) gd.getNextNumber());
			settings.setOutputType(gd.getNextChoiceIndex());

			update();
			return true;
		}

		private void draw()
		{
			drawLabel();
			drawCentre();
		}

		private void update()
		{
			updateLabel();
			updateCentre();
		}

		int cropBorder = -1;
		int cropStart = -1;
		int cropEnd = -1;
		int outputType = -1;
		int psfMagnification;
		boolean singlePrecision = false;

		private void drawLabel()
		{
			cropBorder = settings.getCropBorder();
			cropStart = settings.getCropStart();
			cropEnd = settings.getCropEnd();
			outputType = settings.getOutputType();
			psfMagnification = settings.getPsfMagnification();
			singlePrecision = settings.getSinglePrecision();

			int[] dimensions = psf.getDimensions();
			// Set limits
			int[] min = new int[] { cropBorder, cropBorder, cropStart };
			int[] max = new int[] { dimensions[0] - cropBorder, dimensions[1] - cropBorder, dimensions[2] - cropEnd };
			dimensions[0] -= 2 * cropBorder;
			dimensions[1] -= 2 * cropBorder;
			dimensions[2] -= (cropStart + cropEnd);

			Size size = CustomTricubicInterpolatingFunction.estimateSize(dimensions);
			if (outputType == OUTPUT_TYPE_IMAGE_PSF)
			{
				Size next = size.enlarge(psfMagnification);
				long future = next.getTotalFunctionPoints() * 4; // 4 bytes for a float
				String futureS = Utils.rounded((double) (future / 1048576));
				label1.setText(String.format("Size for Image PSF = %s MB", futureS));
			}
			else
			{
				long current = size.getMemoryFootprint(singlePrecision);
				String currentS = Utils.rounded((double) (current / 1048576));
				label1.setText(String.format("Size of Cubic Spline = %s MB", currentS));
			}

			// Draw ROI in the image
			psfOut[0].setRoi(min[0], min[1], max[0] - min[0], max[1] - min[1]);
			for (int i = 0; i < 3; i++)
			{
				int xd = Projection.getXDimension(i);
				int yd = Projection.getYDimension(i);
				psfOut[i + 1].setRoi(min[xd], min[yd], max[xd] - min[xd], max[yd] - min[yd]);
			}
		}

		private void updateLabel()
		{
			// Get the dimensions after the crop
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
							//@formatter:off
							while (cropBorder != settings.getCropBorder() || 
									cropStart != settings.getCropStart() ||
									cropEnd != settings.getCropEnd() ||
									outputType != settings.getOutputType() ||
									psfMagnification != settings.getPsfMagnification() ||
									singlePrecision != settings.getSinglePrecision())
							//@formatter:on
							{
								drawLabel();
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

		private int centre = -1;

		private void drawCentre()
		{
			centre = slice;

			if (psfOut[0].getSlice() != centre)
			{
				psfOut[0].setSlice(centre);
				psfOut[0].resetDisplayRange();
				psfOut[0].updateAndDraw();
			}

			// Mark projections using an overlay
			// X-projection
			psfOut[1].setOverlay(new Line(centre, 0, centre, psfOut[1].getHeight()), Color.blue, 1, Color.blue);
			// Y-projection
			psfOut[2].setOverlay(new Line(0, centre, psfOut[2].getWidth(), centre), Color.blue, 1, Color.blue);
		}

		private void updateCentre()
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
							while (centre != slice)
							{
								drawCentre();
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
		ExtractedPSF first = psfs[0];

		// PSFs can have different stack sizes. XY size is the same.
		// Find the biggest insert before and after the centre to find
		// the combined stack size.

		int before = first.stackZCentre;
		int after = first.psf.length - first.stackZCentre;
		for (int i = 1; i < psfs.length; i++)
		{
			before = Math.max(before, psfs[i].stackZCentre);
			after = Math.max(after, psfs[i].psf.length - psfs[i].stackZCentre);
		}

		int totalDepth = before + after;
		float[][] psf = new float[totalDepth][first.psf[0].length];
		int size = first.maxx;
		int[] count = new int[totalDepth];
		for (int i = 0; i < psfs.length; i++)
		{
			int offset = before - psfs[i].stackZCentre;
			for (int j = 0; j < psfs[i].psf.length; j++)
			{
				float[] from = psfs[i].psf[j];
				float[] to = psf[j + offset];
				count[j + offset]++;
				for (int k = 0; k < to.length; k++)
					to[k] += from[k];
			}
		}
		// Q. Should the normalisation be done?
		for (int j = 0; j < psf.length; j++)
		{
			float[] to = psf[j];
			int c = count[j];
			if (c != 0)
			{
				for (int k = 0; k < to.length; k++)
					to[k] /= c;
			}
		}
		BasePoint centre = new BasePoint(size / 2f, size / 2f, before);
		ExtractedPSF combined = new ExtractedPSF(psf, size, centre, first.magnification);
		combined.stackZCentre = before;
		return combined;
	}

	private boolean showAlignmentDialog()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Use XYZ stack alignment to create a combined PSF");

		CalibrationWriter cw = CalibrationWriter.create(settings.getCalibration());

		// These are ignored for reset
		gd.addChoice("Alignment_mode", ALIGNMENT_MODE, settings.getAlignmentMode());
		gd.addSlider("Z_radius (px)", 0, imp.getStackSize(), settings.getAlignmentZRadius());
		gd.addChoice("PSF_type", PSF_TYPE, settings.getPsfType());
		gd.addNumericField("nm_per_pixel", cw.getNmPerPixel(), 2, 6, "nm");
		gd.addNumericField("nm_per_slice", settings.getNmPerSlice(), 0, 6, "nm");
		PeakFit.addCameraOptions(gd, PeakFit.FLAG_NO_GAIN, cw);

		// For reset
		final TurboList<TextField> tf = new TurboList<TextField>();
		final TurboList<Checkbox> cb = new TurboList<Checkbox>();
		tf.add(gd.addAndGetSlider("Analysis_window", 0, 8, settings.getAnalysisWindow()));
		tf.add(gd.addAndGetSlider("Smoothing", 0.1, 0.5, settings.getSmoothing()));
		tf.add(gd.addAndGetSlider("CoM_z_window", 0, 8, settings.getComWindow()));
		tf.add(gd.addAndGetSlider("CoM_border", 0, 0.5, settings.getComBorder()));
		tf.add(gd.addAndGetSlider("Alignment_magnification", 1, 8, settings.getAlignmentMagnification()));
		cb.add(gd.addAndGetCheckbox("Smooth_stack_signal", settings.getSmoothStackSignal()));
		tf.add(gd.addAndGetSlider("Max_iterations", 1, 20, settings.getMaxIterations()));
		if (settings.getInteractiveMode())
			cb.add(gd.addAndGetCheckbox("Check_alignments", settings.getCheckAlignments()));

		if (Utils.isShowGenericDialog())
		{
			gd.addAndGetButton("Reset", new ActionListener()
			{
				public void actionPerformed(ActionEvent e)
				{
					PSFCreatorSettings defaults = GUIProtosHelper.defaultPSFCreatorSettings;
					int t = 0, c = 0;
					tf.get(t++).setText(Double.toString(defaults.getAnalysisWindow()));
					tf.get(t++).setText(Double.toString(defaults.getSmoothing()));
					tf.get(t++).setText(Integer.toString(defaults.getComWindow()));
					tf.get(t++).setText(Double.toString(defaults.getComBorder()));
					tf.get(t++).setText(Integer.toString(defaults.getAlignmentMagnification()));
					cb.get(c++).setState(defaults.getSmoothStackSignal());
					tf.get(t++).setText(Integer.toString(defaults.getMaxIterations()));
					if (PSFCreator.this.settings.getInteractiveMode())
						cb.get(c++).setState(defaults.getCheckAlignments());

					// Reset later options too
					PSFCreator.this.settings.setPsfMagnification(defaults.getPsfMagnification());
					PSFCreator.this.settings.setWindow(defaults.getWindow());
					PSFCreator.this.settings.setSmoothStackSignal(defaults.getSmoothStackSignal());
					PSFCreator.this.settings.setComBorder(defaults.getComBorder());
					PSFCreator.this.settings.setCropBorder(0);
					PSFCreator.this.settings.setCropStart(0);
					PSFCreator.this.settings.setCropEnd(0);
				}
			});
		}

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		settings.setAlignmentMode(gd.getNextChoiceIndex());
		settings.setAlignmentZRadius((int) gd.getNextNumber());
		settings.setPsfType(gd.getNextChoiceIndex());
		cw.setNmPerPixel(gd.getNextNumber());
		settings.setNmPerSlice(gd.getNextNumber());
		cw.setCameraType(SettingsManager.getCameraTypeValues()[gd.getNextChoiceIndex()]);
		settings.setAnalysisWindow((int) gd.getNextNumber());
		settings.setSmoothing(gd.getNextNumber());
		settings.setComWindow((int) gd.getNextNumber());
		settings.setComBorder(gd.getNextNumber());
		settings.setAlignmentMagnification((int) gd.getNextNumber());
		settings.setSmoothStackSignal(gd.getNextBoolean());
		settings.setMaxIterations((int) gd.getNextNumber());
		if (settings.getInteractiveMode())
		{
			checkAlignments = gd.getNextBoolean();
			settings.setCheckAlignments(checkAlignments);
		}
		//settings.setPsfMagnification((int) gd.getNextNumber());

		gd.collectOptions();

		nmPerPixel = cw.getNmPerPixel();

		settings.setCalibration(cw.getBuilder());
		SettingsManager.writeSettings(settings);

		// Check arguments
		try
		{
			Parameters.isPositive("nm/pixel", nmPerPixel);
			Parameters.isPositive("nm/slice", settings.getNmPerSlice());
			// Since we do a local background estimation for each extracted PSF then we
			// do not need the bias for non sCMOS cameras.
			//if (!cw.isSCMOS())
			//	Parameters.isAboveZero("Bias", cw.getBias());
			Parameters.isEqualOrAbove("Projection magnification", settings.getAlignmentMagnification(), 1);
			Parameters.isEqualOrAbove("Max iterations", settings.getMaxIterations(), 1);
			Parameters.isEqualOrAbove("PSF magnification", settings.getPsfMagnification(), 1);
			Parameters.isAbove("Smoothing", settings.getSmoothing(), 0);
			Parameters.isBelow("Smoothing", settings.getSmoothing(), 1);
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
					psfs[index] = new ExtractedPSF(image, w, h, centres[index], boxRadius, zRadius,
							settings.getAlignmentMagnification());
					// Do this here within the thread
					psfs[index].createProjections();

					//psfs[index].show(TITLE + index);
				}
			}));
		}

		Utils.waitForCompletion(futures);

		return psfs;
	}

	private static float[][] getLimits(float[][] psf, int x, int y, int n)
	{
		BlockMeanFilter filter = new BlockMeanFilter();
		float[][] limits = new float[2][psf.length];
		for (int zz = 0; zz < psf.length; zz++)
		{
			float[] data = psf[zz];
			if (n > 0)
			{
				data = data.clone();
				filter.rollingBlockFilterInternal(data, x, y, n);
			}
			float[] l = findLimits(data, x, y, n);
			limits[0][zz] = l[0];
			limits[1][zz] = l[1];
			//float[] l2 = Maths.limits(data);
			//System.out.printf("%f - %f vs %f - %f\n", l[0], l[1], l2[0], l2[1]);
		}
		return limits;
	}

	private static float getBackground(float[][] psf, int x, int y, int n)
	{
		float[][] limits = getLimits(psf, x, y, n);
		return Maths.min(limits[0]);
		//return Maths.min(new Smoother().smooth(limits[0]).getFSmooth());
	}

	/**
	 * Find min/max value. This must put the entire region within the image
	 * 
	 * @param data
	 *            the data
	 * @param maxx
	 *            the maxx
	 * @param maxy
	 *            the maxy
	 * @param n
	 *            the block size of the region
	 * @return the limits
	 */
	private static float[] findLimits(float[] data, int maxx, int maxy, int n)
	{
		int[] l = findLimitsIndex(data, maxx, maxy, n);
		return new float[] { data[l[0]], data[l[1]] };
	}

	/**
	 * Find min/max index. This must put the entire region within the image
	 * 
	 * @param data
	 *            the data
	 * @param maxx
	 *            the maxx
	 * @param maxy
	 *            the maxy
	 * @param n
	 *            the block size of the region
	 * @return the min/max index
	 */
	private static int[] findLimitsIndex(float[] data, int maxx, int maxy, int n)
	{
		int min = n * maxx + n;
		int max = min;
		for (int y = n; y < maxy - n; y++)
		{
			for (int x = n, i = y * maxx + n; x < maxx - n; x++, i++)
			{
				if (data[i] < data[min])
					min = i;
				else if (data[i] > data[max])
					max = i;
			}
		}
		return new int[] { min, max };
	}

	/**
	 * Gets the constrast as the ratio between max and min value.
	 *
	 * @param min
	 *            the min value
	 * @param max
	 *            the max value
	 * @return the constrast
	 */
	@SuppressWarnings("unused")
	private static double[] getConstrast(float[] min, float[] max)
	{
		double[] c = new double[min.length];
		for (int i = 0; i < min.length; i++)
			c[i] = (double) max[i] / min[i];
		return c;
	}

	private static double[] getCentreOfMass(float[][] psf, int w, int h)
	{
		double cx = 0;
		double cy = 0;
		double cz = 0;
		double sumXYZ = 0;
		for (int z = 0; z < psf.length; z++)
		{
			float[] data = psf[z];
			double sumXY = 0;
			for (int y = 0, j = 0; y < h; y++)
			{
				double sumX = 0;
				for (int x = 0; x < w; x++)
				{
					float f = data[j++];
					sumX += f;
					cx += f * x;
				}
				sumXY += sumX;
				cy += sumX * y;
			}
			cz += sumXY * z;
			sumXYZ += sumXY;
		}
		// Find centre with 0.5 as the centre of the pixel
		cx = 0.5 + cx / sumXYZ;
		cy = 0.5 + cy / sumXYZ;
		cz = 0.5 + cz / sumXYZ;
		return new double[] { cx, cy, cz };
	}

	private static double[] getCentreOfMassXY(float[][] psf, int w, int h, int zCentre, int n, int border)
	{
		double cx = 0;
		double cy = 0;
		double sumXYZ = 0;
		for (int z = zCentre - n; z <= zCentre + n; z++)
		{
			// Check bounds
			if (z < 0 || z >= psf.length)
			{
				// Ignore. This will break the com but we only use it for the XY position
				continue;
			}

			float[] data = psf[z];

			double sumXY = 0;
			for (int y = border; y < h - border; y++)
			{
				double sumX = 0;
				for (int x = border, j = y * w + border; x < w - border; x++)
				{
					float f = data[j++];
					sumX += f;
					cx += f * x;
				}
				sumXY += sumX;
				cy += sumX * y;
			}
			sumXYZ += sumXY;
		}
		// Find centre with 0.5 as the centre of the pixel
		cx = 0.5 + cx / sumXYZ;
		cy = 0.5 + cy / sumXYZ;
		return new double[] { cx, cy };
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

		@SuppressWarnings("unused")
		int[] getDimensions()
		{
			return new int[] { x, y, z };
		}

		/**
		 * Gets the centre of mass of each projection and then combine.
		 *
		 * @return the centre of mass
		 */
		public double[] getCentreOfMass()
		{
			// Start in the centre of the image
			double[] com = new double[] { x / 2.0, y / 2.0, z / 2.0 };
			for (int i = 0; i < 3; i++)
			{
				// This should be croped around z-centre but the method isn't used so ignore this
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
			float[] data = (float[]) fp.getPixels();
			double[] com = centreOfMass(data, fp.getWidth(), fp.getHeight());

			//			// Compare will background subtracted centre
			//			float min = Maths.min(data);
			//			float[] data2 = new float[data.length];
			//			for (int i = 0; i < data.length; i++)
			//				data2[i] = data[i] - min;
			//			double[] com2 = centreOfMass(data2, fp.getWidth(), fp.getHeight());
			//			
			//			System.out.printf("COM shift %f,%f\n", com[0]-com2[0], com[1]-com2[1]);

			return com;
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

		FloatProcessor getProjection(int i, int minz, int maxz)
		{
			FloatProcessor fp = getProjection(i);
			if (i < Z)
			{
				int range = maxz - minz + 1;
				switch (i)
				{
					case X:
						fp.setRoi(minz, 0, range, y);
						break;
					case Y:
						fp.setRoi(0, minz, x, range);
					default:
						break;
				}
				fp = (FloatProcessor) fp.crop();
			}
			return fp;
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

	private class ExtractedPSF
	{
		BasePoint centre;
		/**
		 * The centre of the stack. Used to crop the image around the centre for alignment.
		 * (Note that the centre in XY is the middle pixel).
		 * <p>
		 * Also used as a relative position when combining PSFs.
		 */
		int stackZCentre;
		final float[][] psf;
		final int maxx, maxy;
		float background;
		Projection projection;
		final int magnification;

		ExtractedPSF(float[][] psf, int size, BasePoint centre, int magnification)
		{
			this.centre = centre;
			this.magnification = magnification;
			stackZCentre = -1; // Not used
			this.psf = psf;
			maxx = size;
			maxy = size;
		}

		/**
		 * Crop to Z centre.
		 *
		 * @param zCentre
		 *            the z centre (1-based index)
		 */
		public ExtractedPSF cropToZCentre(int zCentre)
		{
			if (zCentre < 1 || zCentre > psf.length)
				throw new IllegalStateException("Cannot crop outside the PSF: " + zCentre);
			zCentre--; // Make 0-based index
			// Number of slices before and after the centre to make it even
			int d = Math.min(zCentre, psf.length - zCentre - 1);
			int from = zCentre - d;
			int to = zCentre + d + 1;
			float[][] psf = Arrays.copyOfRange(this.psf, from, to);
			int zShift = -from;
			ExtractedPSF p = new ExtractedPSF(psf, maxx, centre.shift(0, 0, -zShift), magnification);
			p.stackZCentre = zCentre + zShift;
			return p;
		}

		ExtractedPSF(float[][] psf, int maxx, int maxy)
		{
			this.centre = new BasePoint(0, 0, 0);
			this.magnification = 1;
			stackZCentre = -1; // Not used
			this.psf = psf;
			this.maxx = maxx;
			this.maxy = maxy;
		}

		ExtractedPSF(float[][] image, int w, int h, BasePoint centre, int boxRadius, int zRadius, int magnification)
		{
			this.centre = centre;
			this.magnification = magnification;

			// Build using tricubic interpolation.

			// Account for the centre of the XY pixel being 0.5 by convention (e.g. floating point ROIs).
			// The centre of the Z pixel is 0 (to support a user selected integer z-centre)
			double cx = centre.getX() - 0.5;
			double cy = centre.getY() - 0.5;
			double cz = centre.getZ();
			int icz = (int) cz;

			// We want to sample NxNxN points per voxel (N=magnification).
			// Create the sample points
			double pc = 1.0 / magnification; // Pixel centre in scaled image
			CubicSplinePosition[] sx = createCubicSplinePosition(cx, pc);
			CubicSplinePosition[] sy = createCubicSplinePosition(cy, pc);
			CubicSplinePosition[] sz = createCubicSplinePosition(cz, pc);

			// Create the interpolation tables
			double[][] tables = new double[magnification * magnification * magnification][];
			for (int z = 0, i = 0; z < magnification; z++)
				for (int y = 0; y < magnification; y++)
					for (int x = 0; x < magnification; x++)
					{
						tables[i++] = CustomTricubicFunction.computePowerTable(sx[x], sy[y], sz[z]);
					}

			// Find where centre is within the sample points
			int ix = findCentre(sx, cx);
			int iy = findCentre(sy, cy);
			int iz = findCentre(sz, cz);

			// Find the first and last voxels where we can interpolate. 
			// Correct interpolation needs values at this voxel and the next but we pad by duplication
			// if the range overlaps the XY edge. This should not happen if the user selects good XY centres.
			int lx = (int) cx - boxRadius;
			int ux = (int) cx + boxRadius;
			int ly = (int) cy - boxRadius;
			int uy = (int) cy + boxRadius;
			// Note: we use the full range of the stack.
			int lz, uz;
			if (zRadius == 0)
			{
				lz = 0;
				uz = image.length - 1;
			}
			else
			{
				lz = Math.max(0, (int) cz - zRadius);
				uz = Math.min(image.length - 1, (int) cz + zRadius);
			}

			// Extract the data for interpolation. 
			// For correct interpolation we need an extra +1 after the upper limit on each axis.
			// Note: Cubic interpolation actually requires -1 and +2 around interpolation point 
			// 0 but for simplicity we allow the gradient to evaluate to zero at the bounds thus
			// requiring a range of 0 to 1 for each point. The bounds should not be important
			// if the PSF reduces to nothing at the edge.
			int rangex = ux - lx + 2;
			int rangey = uy - ly + 2;
			// Clip z as we do not pad the stack. 
			// This allows +1 after only if uz is clipped by the z-radius.
			int rangez = Math.min(image.length, uz - lz + 2);

			// Build an interpolating function
			// We pad with duplicate pixels if at the bounds
			int[] xi = new int[rangex];
			int[] yi = new int[rangey];
			for (int i = 0; i < xi.length; i++)
			{
				// Keep within data bounds for the indices
				xi[i] = Maths.clip(0, w - 1, lx + i);
			}
			for (int i = 0; i < yi.length; i++)
			{
				// Keep within data bounds for the indices
				yi[i] = Maths.clip(0, h - 1, ly + i);
			}

			// Extract to a stack to extract the background
			float[][] fval = new float[rangez][rangex * rangey];

			for (int z = 0; z < rangez; z++)
			{
				// Note: rangez is already clipped to the image size
				float[] data = image[lz + z];
				float[] slice = fval[z];

				for (int y = 0, i = 0; y < yi.length; y++)
				{
					final int index = w * yi[y];
					for (int x = 0; x < xi.length; x++, i++)
					{
						slice[i] = data[index + xi[x]];
					}
				}
			}

			int window = (int) Math.round(settings.getAnalysisWindow());
			background = getBackground(fval, rangex, rangey, window);
			double[] sum = new double[rangez];
			for (int z = 0; z < rangez; z++)
			{
				float[] slice = fval[z];
				double s = 0;
				for (int i = slice.length; i-- > 0;)
				{
					slice[i] -= background;
					s += slice[i];
				}
				sum[z] = s;
			}

			// Normalisation so the max intensity frame is one.
			// This allows adding all PSFs fairly.
			// This must be done after the background subtraction.
			// Note that this will only effect the interpolation in Z.
			// It will not effect how the function interpolates in XY.
			//
			// If smoothing is off then the enlarged PSF has a much noiser signal
			if (settings.getSmoothStackSignal())
			{
				// Smooth the intensity
				double[] ssum = new Smoother().smooth(sum).getDSmooth();

				// Compute normalisation and apply.
				SimpleArrayUtils.multiply(ssum, 1.0 / Maths.max(ssum));
				for (int z = 0; z < rangez; z++)
				{
					if (sum[z] != 0)
						SimpleArrayUtils.multiply(fval[z], ssum[z] / sum[z]);
				}
			}
			else
			{
				final double norm = 1.0 / Maths.max(sum);
				for (int z = 0; z < rangez; z++)
				{
					SimpleArrayUtils.multiply(fval[z], norm);
				}
			}

			FloatStackTrivalueProvider value = new FloatStackTrivalueProvider(fval, rangex, rangey);

			// The range x and y is controled to put the PSF in the middle of the slice
			maxx = maxy = 2 * boxRadius * magnification + 1;
			// The output z-stack interpolates all points for each voxel so the middle may not be in the centre of the stack
			int maxz = (rangez - 1) * magnification;
			psf = new float[maxz][maxx * maxy];

			// Compute the centre of the stack: (icz - lz) == number of slices before the centre slice
			stackZCentre = (icz - lz) * magnification + iz;

			// Interpolate: NxNxN points per voxel
			final int n = magnification;

			// Visit each voxel that can be interpolated.
			rangez--;
			rangey--;
			rangex--;

			for (int z = 0, pz = 0; z < rangez; z++, pz += n)
			{
				// We use pointers to the position in the interpolation tables so that the initial edge is
				// treated differently. This makes the X/Y dimension the same for all PSFs
				for (int y = 0, yy = iy, py = 0; y < rangey; y++, py += (n - yy), yy = 0)
				{
					for (int x = 0, xx = ix, px = 0; x < rangex; x++, px += (n - xx), xx = 0)
					{
						// Build the interpolator
						CustomTricubicFunction f = CustomTricubicInterpolator.create(value, x, y, z);

						// Sample NxNxN. 
						// The initial edge is handled by the position indices (pz,py,px).
						// The final edge is handled by the bounds of the PSF (psf.length, maxy, maxx).
						for (int zzz = 0, ppz = pz; zzz < n && ppz < psf.length; zzz++, ppz++)
						{
							float[] data = psf[ppz];
							for (int yyy = yy, ppy = py; yyy < n && ppy < maxy; yyy++, ppy++)
							{
								for (int xxx = xx, ppx = px; xxx < n && ppx < maxx; xxx++, ppx++)
								{
									double[] table = tables[xxx + n * (yyy + n * zzz)];
									data[maxx * ppy + ppx] = (float) f.value(table);
									if (Float.isNaN(data[maxx * ppy + ppx]))
									{
										System.out.printf("%d,%d,%d\n", x, y, z);
									}
								}
							}
						}
					}
				}
			}
		}

		/**
		 * Creates the cubic spline position sampling 'magnification' points within the range [0-1] from the given
		 * centre.
		 *
		 * @param centre
		 *            the centre
		 * @param pc
		 *            the pixel change (1/magnification)
		 * @return the cubic spline position
		 */
		private CubicSplinePosition[] createCubicSplinePosition(double centre, double pc)
		{
			CubicSplinePosition[] c = new CubicSplinePosition[magnification];
			// Find the first position in the range [0-1]
			centre -= Math.floor(centre);
			int j = 0;
			while (centre + j * pc >= 0)
				j--;
			for (int i = 0; i < c.length; i++)
			{
				j++;
				c[i] = new CubicSplinePosition(centre + j * pc);
			}
			return c;
		}

		private int findCentre(CubicSplinePosition[] sx, double centre)
		{
			centre -= Math.floor(centre);
			for (int i = 0; i < sx.length; i++)
				if (sx[i].getX() == centre)
					return i;
			throw new IllegalStateException();
		}

		void createProjections()
		{
			if (projection == null)
				projection = new Projection(psf, maxx, maxy);
		}

		ImagePlus[] show(String title)
		{
			int slice = (stackZCentre != -1) ? stackZCentre :
			// The centre is in the original scale so magnify
					centre.getZint() * magnification + 1;
			return show(title, slice);
		}

		ImagePlus[] show(String title, int slice)
		{
			ImagePlus[] out = new ImagePlus[4];
			ImageStack stack = new ImageStack(maxx, maxy);
			for (float[] pixels : psf)
				stack.addSlice(null, pixels);
			ImagePlus imp = Utils.display(title, stack);
			int n = psf.length;
			imp.setSlice(Maths.clip(1, n, slice));
			setCalibration(imp, 2);

			out[0] = imp;

			// Show the projections
			if (projection == null)
				return out;

			out[1] = setCalibration(Utils.display(title + " X-projection (ZY)", getProjection(0, false)), 0);
			out[2] = setCalibration(Utils.display(title + " Y-projection (XZ)", getProjection(1, false)), 1);
			out[3] = setCalibration(Utils.display(title + " Z-projection (XY)", getProjection(2, false)), 2);
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
					c.pixelWidth = settings.getNmPerSlice() / magnification;
					c.pixelHeight = nmPerPixel / magnification;
					break;
				case Projection.Y:
					c.pixelWidth = nmPerPixel / magnification;
					c.pixelHeight = settings.getNmPerSlice() / magnification;
					break;
				case Projection.Z:
					c.pixelWidth = nmPerPixel / magnification;
					c.pixelHeight = nmPerPixel / magnification;
					c.pixelDepth = settings.getNmPerSlice() / magnification;
					break;
			}
			return c;
		}

		FloatProcessor getProjection(int i, boolean crop)
		{
			if (crop)
			{
				int pad = getZPadding();
				return projection.getProjection(i, stackZCentre - pad, stackZCentre + pad);
			}
			return projection.getProjection(i);
		}

		public ImageStack getImageStack(boolean crop)
		{
			int min, max;
			if (crop)
			{
				int pad = getZPadding();
				min = stackZCentre - pad;
				max = stackZCentre + pad;
			}
			else
			{
				min = 0;
				max = psf.length - 1;
			}
			ImageStack stack = new ImageStack(maxx, maxy, max - min + 1);
			for (int i = min; i <= max; i++)
			{
				stack.setPixels(psf[i], i + 1);
			}
			return stack;
		}

		/**
		 * Gets the z padding to have an equal number of slices before and after the current stack z centre.
		 *
		 * @return the z padding
		 */
		int getZPadding()
		{
			return Math.min(stackZCentre, psf.length - stackZCentre - 1);
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
					// Centre in X,Y,Z refer to the position extracted from the image
					centre.getX() + translation[0] / magnification, centre.getY() + translation[1] / magnification,
					centre.getZ() + translation[2] / magnification);
		}

		/**
		 * Compute the centre of mass.
		 *
		 * @param useProjection
		 *            the use projection
		 * @return the centre of mass
		 */
		@SuppressWarnings("unused")
		public double[] getCentreOfMass(boolean useProjection)
		{
			if (useProjection)
				return projection.getCentreOfMass();
			return PSFCreator.getCentreOfMass(psf, maxx, maxy);
		}

		/**
		 * Compute the centre of mass and then the shift of the CoM from the centre of the
		 * image.
		 *
		 * @return the centre of mass shift
		 */
		public double[] getCentreOfMassXYShift(int zCentre)
		{
			double[] shift = PSFCreator.getCentreOfMassXY(psf, maxx, maxy, zCentre, settings.getComWindow(),
					getCoMXYBorder(maxx, maxy));
			// Turn into a shift relative to the centre
			int[] d = new int[] { maxx, maxy };
			for (int i = 0; i < d.length; i++)
			{
				shift[i] -= d[i] / 2.0;
				// Account for magnification
				shift[i] /= magnification;
			}
			return shift;
		}

		public ExtractedPSF enlarge(int n, ExecutorService threadPool)
		{
			if (n <= 1)
				return this;

			IJTrackProgress progress = new IJTrackProgress();
			FloatStackTrivalueProcedure p = new FloatStackTrivalueProcedure();
			FloatStackTrivalueProvider fval = new FloatStackTrivalueProvider(psf, maxx, maxy);

			// We can enlarge by interpolation between the start and end
			// by evenly sampling each spline node.
			// Do this dynamically. It is slower than creating the entire interpolating function
			// but uses much less memory.
			Utils.showStatus("Enlarging ... interpolating");

			//@formatter:off
			new CustomTricubicInterpolator.Builder()
					.setExecutorService(threadPool)
					.setProgress(progress)
					.build()
					.sample(fval, n, p);
			//@formatter:on

			IJ.showStatus("");

			// Enlarge the centre.  
			BasePoint newCentre = new BasePoint(
					// XY pixels are centred at 0.5.
					n * (centre.getX() - 0.5f) + 0.5f, n * (centre.getY() - 0.5f) + 0.5f,
					// Z stack at 0
					n * centre.getZ());
			ExtractedPSF psf = new ExtractedPSF(p.value, p.x.length, newCentre, n);
			// We will have enlarged all the slices before the centre n times
			if (stackZCentre != -1)
				psf.stackZCentre = (this.stackZCentre - 1) * n + 1;
			return psf;
		}

		public int[] getDimensions()
		{
			return new int[] { maxx, maxy, psf.length };
		}
	}

	/**
	 * Align the PSFs with the combined PSF.
	 *
	 * @param combined
	 *            the combined
	 * @param psfs
	 *            the psfs
	 * @return The XYZ translations for each PSF
	 */
	private float[][] align(ExtractedPSF combined, final ExtractedPSF[] psfs)
	{
		if (settings.getAlignmentMode() == ALIGNMENT_MODE_2D)
			return align2D(combined, psfs);
		return align3D(combined, psfs);
	}

	/**
	 * Align the PSFs with the combined PSF using the Image2DAligner class to align the 2D max intensity projections.
	 * The final alignment shift is the average of the shift from two projection alignments for each dimension.
	 *
	 * @param combined
	 *            the combined
	 * @param psfs
	 *            the psfs
	 * @return The XYZ translations for each PSF
	 */
	private float[][] align2D(ExtractedPSF combined, final ExtractedPSF[] psfs)
	{
		// Note: For alignment we crop the X/Y projections around the current z-centre
		// so the middle of the 2D image is the middle of the projection.

		int n = psfs.length * 3;
		List<Future<?>> futures = new TurboList<Future<?>>(n);

		final Image2DAligner[] align = new Image2DAligner[3];
		for (int i = 0; i < 3; i++)
		{
			align[i] = new Image2DAligner();
			FloatProcessor fp1 = combined.getProjection(i, true);
			align[i].setReference(fp1); // No need to set the bounds as the PSF will be smaller
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
						double[] result = align[ii].copy().align(psf.getProjection(ii, true), 10);
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

	/**
	 * Align the PSFs with the combined PSF using the Image2DAligner class to align the 2D max intensity projections.
	 * The final alignment shift is the average of the shift from two projection alignments for each dimension.
	 *
	 * @param combined
	 *            the combined
	 * @param psfs
	 *            the psfs
	 * @return The XYZ translations for each PSF
	 */
	private float[][] align3D(ExtractedPSF combined, final ExtractedPSF[] psfs)
	{
		// Note: For alignment we extract each PSF around the current z-centre
		// so the middle of the stack is the middle of the PSF.

		List<Future<?>> futures = new TurboList<Future<?>>(psfs.length);

		final Image3DAligner align = new Image3DAligner();
		align.setReference(combined.getImageStack(true));

		final float[][] results = new float[psfs.length][3];

		for (int j = 0; j < psfs.length; j++)
		{
			final int jj = j;
			futures.add(threadPool.submit(new Runnable()
			{
				public void run()
				{
					ExtractedPSF psf = psfs[jj];
					double[] result = align.copy().align(psf.getImageStack(true), 10);
					for (int i = 0; i < 3; i++)
						results[jj][i] = (float) result[i];
				}
			}));
		}

		Utils.waitForCompletion(futures);

		return results;
	}

	/**
	 * Align the PSFs with the combined PSF using the gdsc.core.ij.AlignImagesFFT class
	 *
	 * @param combined
	 *            the combined
	 * @param psfs
	 *            the psfs
	 * @return The XYZ translations for each PSF
	 */
	@SuppressWarnings("unused")
	private float[][] align2(ExtractedPSF combined, final ExtractedPSF[] psfs)
	{
		int n = psfs.length * 3;
		List<Future<?>> futures = new TurboList<Future<?>>(n);

		final AlignImagesFFT[] align = new AlignImagesFFT[3];
		final Rectangle[] bounds = new Rectangle[3];
		for (int i = 0; i < 3; i++)
		{
			align[i] = new AlignImagesFFT();
			FloatProcessor fp1 = combined.getProjection(i, true);
			FloatProcessor fp2 = psfs[0].getProjection(i, true);
			align[i].init(fp1, WindowMethod.TUKEY, true);
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
						double[] result = align[ii].align(psf.getProjection(ii, true), WindowMethod.TUKEY, bounds[ii],
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

	private BasePoint[] updateUsingCentreOfMassXYShift(double[] shift, double shiftd, ExtractedPSF combined,
			BasePoint[] centres)
	{
		float dx = (float) shift[0];
		float dy = (float) shift[1];
		float dz = 0;
		Utils.log("Combined PSF has CoM shift %s,%s (%s)", rounder.toString(shift[0]), rounder.toString(shift[1]),
				rounder.toString(shiftd));
		for (int i = 0; i < centres.length; i++)
		{
			centres[i] = centres[i].shift(dx, dy, dz);
		}
		return centres;
	}

	/**
	 * Apply the XY and Z weighting window. The background is subtracted from the data and the window applied.
	 *
	 * @param data
	 *            the data
	 * @param z
	 *            the z position of the PSF data
	 * @param wx
	 *            the XY weighting
	 * @param wz
	 *            the Z weighting
	 * @param background
	 *            the background
	 * @return the sum of the image
	 */
	private static double applyWindow(float[] data, int z, double[] wx, double[] wz, float background)
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
		double sum = 0;
		for (int y = 0, i = 0; y < size; y++)
		{
			double weight = w[y];
			for (int x = 0; x < size; x++, i++)
			{
				// Subtract background
				float f = (data[i] - background);
				if (f <= 0)
					data[i] = 0;
				else
					// Window function and normalise
					sum += data[i] = (float) (f * w[x] * weight);
			}
		}
		return sum;
	}
}
