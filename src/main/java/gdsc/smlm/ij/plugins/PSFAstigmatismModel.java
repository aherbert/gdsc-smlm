package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.Color;
import java.awt.Rectangle;
import java.util.Arrays;

import org.apache.commons.math3.analysis.MultivariateMatrixFunction;
import org.apache.commons.math3.analysis.MultivariateVectorFunction;
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer.Optimum;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.util.Precision;

import gdsc.core.data.utils.Rounder;
import gdsc.core.data.utils.RounderFactory;
import gdsc.core.ij.IJLogger;
import gdsc.core.ij.IJTrackProgress;
import gdsc.core.ij.Utils;
import gdsc.core.logging.Ticker;
import gdsc.core.utils.Maths;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.core.utils.TextUtils;
import gdsc.smlm.data.config.CalibrationWriter;
import gdsc.smlm.data.config.FitProtos.DataFilterMethod;
import gdsc.smlm.data.config.FitProtos.DataFilterType;
import gdsc.smlm.data.config.GUIProtos.PSFAstigmatismModelSettings;
import gdsc.smlm.data.config.PSFProtos.PSFType;
import gdsc.smlm.data.config.PSFProtosHelper;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.engine.FitConfiguration;
import gdsc.smlm.engine.FitEngine;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2018 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.engine.FitParameters;
import gdsc.smlm.engine.FitQueue;
import gdsc.smlm.engine.ParameterisedFitJob;
import gdsc.smlm.ij.IJImageSource;
import gdsc.smlm.ij.plugins.PeakFit.FitEngineConfigurationProvider;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.SynchronizedPeakResults;
import gdsc.smlm.results.count.Counter;
import gdsc.smlm.results.procedures.PeakResultProcedure;
import gdsc.smlm.results.procedures.WidthResultProcedure;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.DialogListener;
import ij.gui.ExtendedGenericDialog;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingExtendedGenericDialog;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.WindowOrganiser;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatPolygon;
import ij.process.ImageProcessor;

/**
 * Produces a 2D Gaussian astigmatism model for a 2D astigmatic PSF.
 * <p>
 * The input images must be a z-stack of a PSF.
 */
public class PSFAstigmatismModel implements PlugInFilter
{
	private final static String TITLE = "PSF Astigmatism Model";

	private final static int FLAGS = DOES_16 | DOES_8G | DOES_32 | STACK_REQUIRED | NO_CHANGES;
	private PSFAstigmatismModelSettings.Builder settings;
	private ImagePlus imp;
	private FitEngineConfiguration config;
	private FitConfiguration fitConfig;
	private int cx, cy;
	private MemoryPeakResults results;
	double[] z, x, y, I, sx, sy;
	PlotWindow iPlot, xyPlot, sPlot;
	int minz, maxz;
	double[] fitZ, fitSx, fitSy;
	double[] parameters;

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

		settings = SettingsManager.readPSFAstigmatismModelSettings(0).toBuilder();

		guessScale();

		gd.addMessage("Use Gaussian 2D PSF fitting to create an astigmatism z-model");

		gd.addNumericField("nm_per_slice", settings.getNmPerSlice(), 0);

		gd.showDialog();

		SettingsManager.writeSettings(settings);

		if (gd.wasCanceled())
			return DONE;

		settings.setNmPerSlice(gd.getNextNumber());

		// Check arguments
		try
		{
			Parameters.isPositive("nm/slice", settings.getNmPerSlice());
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return DONE;
		}

		return FLAGS;
	}

	private void guessScale()
	{
		CalibrationWriter cw = CalibrationWriter.create(settings.getCalibration());
		// It does not matter if we already have settings, try and update them anyway
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

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
	 */
	public void run(ImageProcessor ip)
	{
		if (!loadConfiguration())
			return;

		if (!findFitRegion())
			return;

		if (!fitRegion())
			return;

		if (!plotData())
			return;

		if (!fitData())
			return;

	}

	private boolean loadConfiguration()
	{
		// We have a different fit configuration just for the PSF Creator.
		// This allows it to be saved and not effect PeakFit settings.
		config = new FitEngineConfiguration(settings.getFitEngineSettings(), settings.getCalibration(),
				settings.getPsf());
		if (!showConfigurationDialog())
		{
			IJ.error(TITLE, "No fit configuration loaded");
			return false;
		}

		SettingsManager.writeSettings(settings);

		if (fitConfig.getPSFType() != PSFType.TWO_AXIS_GAUSSIAN_2D)
		{
			IJ.error(TITLE, "PSF must be " + PSFProtosHelper.getName(PSFType.TWO_AXIS_GAUSSIAN_2D));
			return false;
		}

		// Simple data filter. This is just used to get the initial estimate of amplitude.
		config.setDataFilterType(DataFilterType.SINGLE);
		config.setDataFilter(DataFilterMethod.GAUSSIAN, 1, 0);

		config.setIncludeNeighbours(false);
		config.configureOutputUnits();
		config.setResidualsThreshold(1);
		config.setDuplicateDistance(0);

		settings.setFitEngineSettings(config.getFitEngineSettings());
		settings.setCalibration(fitConfig.getCalibration());
		settings.setPsf(fitConfig.getPSF());
		SettingsManager.writeSettings(settings);

		return true;
	}

	public boolean showConfigurationDialog()
	{
		fitConfig = config.getFitConfiguration();
		CalibrationWriter calibration = fitConfig.getCalibrationWriter();

		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);
		gd.addMessage("Configuration settings for the single-molecule localisation microscopy plugins");

		PeakFit.addCameraOptions(gd, calibration);
		gd.addNumericField("Calibration (nm/px)", calibration.getNmPerPixel(), 2);
		//gd.addNumericField("Exposure_time (ms)", calibration.getExposureTime(), 2);

		PeakFit.addPSFOptions(gd, fitConfig);

		FitEngineConfigurationProvider provider = new PeakFit.SimpleFitEngineConfigurationProvider(config);
		PeakFit.addFittingOptions(gd, provider);
		gd.addChoice("Fit_solver", SettingsManager.getFitSolverNames(), fitConfig.getFitSolver().ordinal());
		gd.addCheckbox("Log_fit_progress", settings.getLogFitProgress());

		gd.addMessage("--- Peak filtering ---\nDiscard fits that shift; are too low; or expand/contract");

		gd.addCheckbox("Smart_filter", fitConfig.isSmartFilter());
		gd.addCheckbox("Disable_simple_filter", fitConfig.isDisableSimpleFilter());
		gd.addSlider("Shift_factor", 0.01, 2, fitConfig.getCoordinateShiftFactor());
		gd.addNumericField("Signal_strength", fitConfig.getSignalStrength(), 2);
		gd.addNumericField("Min_photons", fitConfig.getMinPhotons(), 0);
		gd.addSlider("Min_width_factor", 0, 0.99, fitConfig.getMinWidthFactor());
		// Fitting may need to be extra wide
		double w = fitConfig.getMaxWidthFactor();
		gd.addSlider("Width_factor", 1.01, Math.max(10, w), w);
		PeakFit.addPrecisionOptions(gd, new PeakFit.SimpleFitConfigurationProvider(fitConfig));

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		// In case a template update the calibration
		calibration = fitConfig.getCalibrationWriter();

		calibration.setCameraType(SettingsManager.getCameraTypeValues()[gd.getNextChoiceIndex()]);
		calibration.setNmPerPixel(gd.getNextNumber());
		calibration.setExposureTime(100); // Arbitrary
		fitConfig.setPSFType(PeakFit.getPSFTypeValues()[gd.getNextChoiceIndex()]);
		config.setFitting(gd.getNextNumber());
		fitConfig.setFitSolver(gd.getNextChoiceIndex());
		settings.setLogFitProgress(gd.getNextBoolean());
		fitConfig.setSmartFilter(gd.getNextBoolean());
		fitConfig.setDisableSimpleFilter(gd.getNextBoolean());
		fitConfig.setCoordinateShiftFactor(gd.getNextNumber());
		fitConfig.setSignalStrength(gd.getNextNumber());
		fitConfig.setMinPhotons(gd.getNextNumber());
		fitConfig.setMinWidthFactor(gd.getNextNumber());
		fitConfig.setWidthFactor(gd.getNextNumber());
		fitConfig.setPrecisionThreshold(gd.getNextNumber());

		gd.collectOptions();

		// Check arguments
		try
		{
			Parameters.isAboveZero("nm per pixel", calibration.getNmPerPixel());
			Parameters.isAboveZero("Initial SD0", fitConfig.getInitialXSD());
			if (fitConfig.getPSF().getParametersCount() > 1)
			{
				Parameters.isAboveZero("Initial SD1", fitConfig.getInitialYSD());
			}
			Parameters.isAboveZero("Fitting_width", config.getFitting());

			if (!fitConfig.isSmartFilter())
			{
				Parameters.isPositive("Coordinate Shift factor", fitConfig.getCoordinateShiftFactor());
				Parameters.isPositive("Signal strength", fitConfig.getSignalStrength());
				Parameters.isPositive("Min photons", fitConfig.getMinPhotons());
				Parameters.isPositive("Min width factor", fitConfig.getMinWidthFactor());
				Parameters.isPositive("Width factor", fitConfig.getMaxWidthFactor());
				Parameters.isPositive("Precision threshold", fitConfig.getPrecisionThreshold());
			}
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		if (gd.invalidNumber())
			return false;

		int flags = PeakFit.FLAG_NO_SAVE;
		if (!PeakFit.configureSmartFilter(config, flags))
			return false;
		PeakFit.configureFitSolver(config, null, null, flags);

		return true;
	}

	private boolean findFitRegion()
	{
		// Get the centre
		Roi roi = imp.getRoi();
		if (roi != null && roi.getType() == Roi.POINT)
		{
			FloatPolygon p = roi.getFloatPolygon();
			int n = p.npoints;
			if (n != 1)
			{
				IJ.error(TITLE, "Require a single point ROI");
				return false;
			}
			cx = (int) p.xpoints[0];
			cy = (int) p.ypoints[0];
			return true;
		}
		IJ.error(TITLE, "Require a single point ROI");
		return false;
	}

	private boolean fitRegion()
	{
		int radius = config.getFittingWidth();

		//IJ.log(config.getFitEngineSettings().toString());

		if (settings.getLogFitProgress())
			fitConfig.setLog(new IJLogger());

		// Create a fit engine
		results = new MemoryPeakResults();
		results.setCalibration(fitConfig.getCalibration());
		results.setPSF(fitConfig.getPSF());
		results.setSortAfterEnd(true);
		results.begin();
		int threadCount = Prefs.getThreads();
		FitEngine engine = new FitEngine(config, SynchronizedPeakResults.create(results, threadCount), threadCount,
				FitQueue.BLOCKING);

		//List<ParameterisedFitJob> jobItems = new ArrayList<ParameterisedFitJob>(stack.getSize());

		IJImageSource source = new IJImageSource(imp);
		source.open();

		Rectangle r1 = new Rectangle(cx - radius, cy - radius, 2 * radius + 1, 2 * radius + 1);
		Rectangle regionBounds = r1.intersection(new Rectangle(source.getWidth(), source.getHeight()));
		// Fit only a spot in the centre
		int x = cx - regionBounds.x;
		int y = cy - regionBounds.y;
		int[] maxIndices = new int[] { y * regionBounds.width + x };

		Ticker ticker = Ticker.createStarted(new IJTrackProgress(), source.getFrames(), threadCount > 1);
		IJ.showStatus("Fitting ...");

		boolean shutdown = false;
		while (!shutdown)
		{
			// Extract the region from each frame
			float[] region = source.next(regionBounds);
			if (region == null)
				break;

			FitParameters params = new FitParameters();
			params.maxIndices = maxIndices.clone();
			int slice = (int) ticker.getCurrent();
			ParameterisedFitJob job = new ParameterisedFitJob(slice, params, slice, region, regionBounds);
			//jobItems.add(job);
			engine.run(job);

			ticker.tick();
			shutdown = IJ.escapePressed();
		}

		if (shutdown)
			IJ.showStatus("Cancelled");
		engine.end(shutdown);
		results.end();

		IJ.showProgress(1);

		if (!shutdown)
			Utils.log("Fit %d/%s", results.size(), TextUtils.pleural(source.getFrames(), "spot"));

		return !shutdown;
	}

	private boolean plotData()
	{
		if (results.size() <= imp.getStackSize() / 2)
		{
			IJ.error(TITLE, "Not enough fit results " + results.size());
			return false;
		}

		final double umPerSlice = settings.getNmPerSlice() / 1000.0;
		//final double nmPerPixel = results.getNmPerPixel();

		z = new double[results.size()];
		x = new double[z.length];
		y = new double[z.length];
		I = new double[z.length];
		final Counter counter = new Counter();

		// We have fit the results so they will be in the preferred units 
		results.forEach(new PeakResultProcedure()
		{
			public void execute(PeakResult peak)
			{
				int i = counter.getAndIncrement();
				z[i] = peak.getFrame() * umPerSlice;
				x[i] = (peak.getXPosition() - cx);
				y[i] = (peak.getYPosition() - cy);
				I[i] = peak.getSignal();
			}
		});

		WidthResultProcedure wp = new WidthResultProcedure(results, DistanceUnit.PIXEL);
		wp.getWxWy();
		sx = SimpleArrayUtils.toDouble(wp.wx);
		sy = SimpleArrayUtils.toDouble(wp.wy);

		WindowOrganiser wo = new WindowOrganiser();

		iPlot = plot(wo, z, "Intensity (photon)", I, "Intensity", null, null);
		xyPlot = plot(wo, z, "Position (px)", x, "X", y, "Y");
		sPlot = plot(wo, z, "Width (px)", sx, "Sx", sy, "Sy");

		wo.tile();

		return true;
	}

	private PlotWindow plot(WindowOrganiser wo, double[] z, String yTitle, double[] y1, String y1Title, double[] y2,
			String y2Title)
	{
		String title = TITLE + " " + yTitle;
		Plot plot = new Plot(title, "Z (μm)", yTitle);
		double[] limits = Maths.limits(y1);
		if (y2 != null)
			limits = Maths.limits(limits, y2);
		double rangex = (z[z.length - 1] - z[0]) * 0.05;
		double rangey = (limits[1] - limits[0]) * 0.05;
		plot.setLimits(z[0] - rangex, z[z.length - 1] + rangex, limits[0] - rangey, limits[1] + rangey);
		plot.setColor(Color.RED);
		plot.addPoints(z, y1, Plot.CIRCLE);
		if (y2 != null)
		{
			plot.setColor(Color.BLUE);
			plot.addPoints(z, y2, Plot.CIRCLE);
			plot.addLegend(y1Title + "\n" + y2Title);
		}
		return Utils.display(title, plot, 0, wo);
	}

	private boolean fitData()
	{
		if (!getRange())
			return false;

		if (!doCurveFit())
			return false;

		return true;
	}

	private boolean getRange()
	{
		minz = 0;
		maxz = z.length - 1;

		sPlot.getImagePlus().killRoi();

		NonBlockingExtendedGenericDialog gd = new NonBlockingExtendedGenericDialog(TITLE);
		gd.addMessage("Select z-range for curve fit");
		gd.addSlider("Min_z", minz, maxz, minz);
		gd.addSlider("Max_z", minz, maxz, maxz);
		gd.addMessage("Select smoothing before curve parameter estimation");
		gd.addSlider("Smoothing", 0.05, 0.5, settings.getSmoothing());
		gd.addDialogListener(new ZDialogListener());
		gd.showDialog();

		// Save settings
		SettingsManager.writeSettings(settings);

		if (gd.wasCanceled())
			return false;

		// Ensure there are enough points to fit
		if (maxz - minz < 10)
		{
			IJ.error(TITLE, "Not enough points for a curve fit");
			return false;
		}

		// Extract data for fit
		fitZ = Arrays.copyOfRange(z, minz, maxz + 1);
		fitSx = Arrays.copyOfRange(sx, minz, maxz + 1);
		fitSy = Arrays.copyOfRange(sy, minz, maxz + 1);

		return true;
	}

	private class ZDialogListener implements DialogListener
	{
		boolean showRoi = Utils.isShowGenericDialog();

		public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
		{
			int oldMinz = minz;
			int oldMaxz = maxz;
			minz = (int) gd.getNextNumber();
			maxz = (int) gd.getNextNumber();
			settings.setSmoothing(gd.getNextNumber());
			if (showRoi && (oldMinz != minz || oldMaxz != maxz))
			{
				Plot plot = sPlot.getPlot();
				int x1 = (int) plot.scaleXtoPxl(z[minz]);
				int x2 = (int) plot.scaleXtoPxl(z[maxz]);
				double[] limits = plot.getLimits();
				int y1 = (int) plot.scaleYtoPxl(limits[3]);
				int y2 = (int) plot.scaleYtoPxl(limits[2]);
				sPlot.getImagePlus().setRoi(new Roi(x1, y1, x2 - x1, y2 - y1));
			}
			return maxz > minz;
		}
	}

	private boolean doCurveFit()
	{
		// Estimate:
		// Focal plane = where width is at a minimum
		// s0x/s0y = the min width of x/y
		// gamma = Half the distance between the focal planes
		// z0 = half way between the two focal planes
		// d = depth of focus

		double[] smoothSx = fitSx;
		double[] smoothSy = fitSy;
		if (settings.getSmoothing() > 0)
		{
			LoessInterpolator loess = new LoessInterpolator(settings.getSmoothing(), 0);
			smoothSx = loess.smooth(fitZ, fitSx);
			smoothSy = loess.smooth(fitZ, fitSy);

			Plot plot = sPlot.getPlot();
			plot.setColor(Color.RED);
			plot.addPoints(fitZ, smoothSx, Plot.LINE);
			plot.setColor(Color.BLUE);
			plot.addPoints(fitZ, smoothSy, Plot.LINE);
			plot.setColor(Color.BLACK);
			plot.updateImage();
		}

		int focalPlaneXindex = SimpleArrayUtils.findMinIndex(smoothSx);
		int focalPlaneYindex = SimpleArrayUtils.findMinIndex(smoothSy);
		double s0x = smoothSx[focalPlaneXindex];
		double s0y = smoothSy[focalPlaneYindex];
		double focalPlaneX = fitZ[focalPlaneXindex];
		double focalPlaneY = fitZ[focalPlaneYindex];
		double gamma = Math.abs(focalPlaneY - focalPlaneX) / 2;
		double z0 = (focalPlaneX + focalPlaneY) / 2;
		double d = (estimateD(focalPlaneXindex, fitZ, smoothSx) + estimateD(focalPlaneYindex, fitZ, smoothSy)) / 2;

		// Start with Ax, Bx, Ay, By as zero.
		double Ax = 0, Bx = 0, Ay = 0, By = 0;

		// Equations assume that x direction is focused above (positive).
		if (focalPlaneXindex < focalPlaneYindex)
			gamma = -gamma;

		// Use an LVM fitter with numerical gradients.
		final double initialStepBoundFactor = 100;
		final double costRelativeTolerance = 1e-10;
		final double parRelativeTolerance = 1e-10;
		final double orthoTolerance = 1e-10;
		final double threshold = Precision.SAFE_MIN;

		// We optimise against both sx and sy as a combined y-value.
		final double[] y = new double[fitZ.length * 2];
		System.arraycopy(fitSx, 0, y, 0, fitSx.length);
		System.arraycopy(fitSy, 0, y, fitSx.length, fitSy.length);

		LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer(initialStepBoundFactor,
				costRelativeTolerance, parRelativeTolerance, orthoTolerance, threshold);

		parameters = new double[9];
		parameters[P_GAMMA] = gamma;
		parameters[P_Z0] = z0;
		parameters[P_D] = d;
		parameters[P_S0X] = s0x;
		parameters[P_AX] = Ax;
		parameters[P_BX] = Bx;
		parameters[P_S0Y] = s0y;
		parameters[P_AY] = Ay;
		parameters[P_BY] = By;

		//plotFit(initialSolution);
		record("Initial", parameters);

		//@formatter:off
		LeastSquaresBuilder builder = new LeastSquaresBuilder()
				.maxEvaluations(Integer.MAX_VALUE)
				.maxIterations(3000)
				.start(parameters)
				.target(y);
		//@formatter:on

		AstigmatismVectorFunction vf = new AstigmatismVectorFunction();
		builder.model(vf, new AstigmatismMatrixFunction());

		LeastSquaresProblem problem = builder.build();

		try
		{
			Optimum optimum = optimizer.optimize(problem);

			parameters = optimum.getPoint().toArray();

			record("Final", parameters);
			plotFit(parameters);
		}
		catch (Exception e)
		{
			IJ.error(TITLE, "Failed to fit curve: " + e.getMessage());
			return false;
		}

		return true;
	}

	/**
	 * Get depth of focus as twice the min width.
	 *
	 * @param min
	 *            the index of the min value of the width
	 * @param z
	 *            the z
	 * @param sx
	 *            the width
	 * @return the estimated depth of focus
	 */
	private double estimateD(int min, double[] z, double[] sx)
	{
		double w = sx[min] * 2; // Twice the min width
		int lower = min;
		while (lower > 0 && sx[lower] < w)
			lower--;
		int upper = min;
		while (upper < sx.length - 1 && sx[upper] < w)
			upper++;
		return z[upper] - z[lower];
	}

	private static final int P_GAMMA = 0;
	private static final int P_Z0 = 1;
	private static final int P_D = 2;
	private static final int P_S0X = 3;
	private static final int P_AX = 4;
	private static final int P_BX = 5;
	private static final int P_S0Y = 6;
	private static final int P_AY = 7;
	private static final int P_BY = 8;

	/**
	 * Gets the standard deviation for the z-depth.
	 *
	 * @param s0
	 *            the width in the focal plane
	 * @param z
	 *            the z
	 * @param one_d2
	 *            one over the depth of focus squared (1/d^2)
	 * @param A
	 *            Empirical constant A for the astigmatism of the PSF
	 * @param B
	 *            Empirical constant B for the astigmatism of the PSF
	 * @return the standard deviation
	 */
	public static double getS(double s0, double z, double one_d2, double A, double B)
	{
		final double z2 = z * z;
		final double z3 = z2 * z;
		final double z4 = z2 * z2;
		// Eq. 17a
		return s0 * Math.sqrt(1 + one_d2 * (z2 + A * z3 + B * z4));
	}

	private class AstigmatismVectorFunction implements MultivariateVectorFunction
	{
		public double[] value(double[] p) throws IllegalArgumentException
		{
			double one_d2 = 1.0 / Maths.pow2(p[P_D]);

			double[] value = new double[fitZ.length * 2];
			double z, z2, z3, z4;

			for (int i = 0, j = fitZ.length; i < fitZ.length; i++, j++)
			{
				// X : z -> z-gamma
				z = fitZ[i] - p[P_Z0] - p[P_GAMMA];
				z2 = z * z;
				z3 = z2 * z;
				z4 = z2 * z2;
				value[i] = p[P_S0X] * Math.sqrt(1 + one_d2 * (z2 + p[P_AX] * z3 + p[P_BX] * z4));
				// Y : z -> z+gamma
				z = fitZ[i] - p[P_Z0] + p[P_GAMMA];
				z2 = z * z;
				z3 = z2 * z;
				z4 = z2 * z2;
				value[j] = p[P_S0Y] * Math.sqrt(1 + one_d2 * (z2 + p[P_AY] * z3 + p[P_BY] * z4));
			}

			return value;
		}
	}

	private class AstigmatismMatrixFunction implements MultivariateMatrixFunction
	{
		public double[][] value(double[] p) throws IllegalArgumentException
		{
			double[] pu = p.clone();
			double[] pl = p.clone();

			// Numerical gradients
			double delta = 1e-6;
			double twoDelta = 2 * delta;

			for (int i = 0; i < p.length; i++)
			{
				pu[i] += delta;
				pl[i] -= delta;
			}

			double one_d2 = 1.0 / Maths.pow2(p[P_D]);
			pu[P_D] = 1.0 / Maths.pow2(pu[P_D]);
			pl[P_D] = 1.0 / Maths.pow2(pl[P_D]);

			double[][] value = new double[fitZ.length * 2][p.length];
			double z, z2, z3, z4, v1, v2;

			// X : z -> z-gamma
			for (int i = 0; i < fitZ.length; i++)
			{
				// X : z -> z-gamma
				z = fitZ[i] - p[P_Z0] - pu[P_GAMMA];
				z2 = z * z;
				z3 = z2 * z;
				z4 = z2 * z2;
				v1 = p[P_S0X] * Math.sqrt(1 + one_d2 * (z2 + p[P_AX] * z3 + p[P_BX] * z4));
				z = fitZ[i] - p[P_Z0] - pl[P_GAMMA];
				z2 = z * z;
				z3 = z2 * z;
				z4 = z2 * z2;
				v2 = p[P_S0X] * Math.sqrt(1 + one_d2 * (z2 + p[P_AX] * z3 + p[P_BX] * z4));

				value[i][P_GAMMA] = (v1 - v2) / twoDelta;

				// Since we use the same delta
				value[i][P_Z0] = value[i][P_GAMMA];

				// All other gradients have the same z
				z = fitZ[i] - p[P_Z0] - p[P_GAMMA];
				z2 = z * z;
				z3 = z2 * z;
				z4 = z2 * z2;

				//@formatter:off
				value[i][P_D] = p[P_S0X] * (
						Math.sqrt(1 + pu[P_D] * (z2 + p[P_AX] * z3 + p[P_BX] * z4))-
						Math.sqrt(1 + pl[P_D] * (z2 + p[P_AX] * z3 + p[P_BX] * z4))) / twoDelta;
				// Analytical gradient
				value[i][P_S0X] = Math.sqrt(1 + one_d2 * (z2 + p[P_AX] * z3 + p[P_BX] * z4)); 
				value[i][P_AX] = p[P_S0X] * (
						Math.sqrt(1 + one_d2 * (z2 + pu[P_AX] * z3 + p[P_BX] * z4))-
						Math.sqrt(1 + one_d2 * (z2 + pl[P_AX] * z3 + p[P_BX] * z4))) / twoDelta;
				value[i][P_BX] = p[P_S0X] * (
						Math.sqrt(1 + one_d2 * (z2 + p[P_AX] * z3 + pu[P_BX] * z4))-
						Math.sqrt(1 + one_d2 * (z2 + p[P_AX] * z3 + pl[P_BX] * z4))) / twoDelta;
				//@formatter:on
			}

			// Y : z -> z+gamma
			for (int i = 0, j = fitZ.length; i < fitZ.length; i++, j++)
			{
				// Y : z -> z+gamma
				z = fitZ[i] - p[P_Z0] + pu[P_GAMMA];
				z2 = z * z;
				z3 = z2 * z;
				z4 = z2 * z2;
				v1 = p[P_S0Y] * Math.sqrt(1 + one_d2 * (z2 + p[P_AY] * z3 + p[P_BY] * z4));
				z = fitZ[i] - p[P_Z0] + pl[P_GAMMA];
				z2 = z * z;
				z3 = z2 * z;
				z4 = z2 * z2;
				v2 = p[P_S0Y] * Math.sqrt(1 + one_d2 * (z2 + p[P_AY] * z3 + p[P_BY] * z4));

				value[j][P_GAMMA] = (v1 - v2) / twoDelta;

				// Since we use the same delta
				value[j][P_Z0] = value[j][P_GAMMA];

				// All other gradients have the same z
				z = fitZ[i] - p[P_Z0] + p[P_GAMMA];
				z2 = z * z;
				z3 = z2 * z;
				z4 = z2 * z2;

				//@formatter:off
				value[j][P_D] = p[P_S0Y] * (
						Math.sqrt(1 + pu[P_D] * (z2 + p[P_AY] * z3 + p[P_BY] * z4))-
						Math.sqrt(1 + pl[P_D] * (z2 + p[P_AY] * z3 + p[P_BY] * z4))) / twoDelta;
				// Analytical gradient
				value[j][P_S0Y] = Math.sqrt(1 + one_d2 * (z2 + p[P_AY] * z3 + p[P_BY] * z4)); 
				value[j][P_AY] = p[P_S0Y] * (
						Math.sqrt(1 + one_d2 * (z2 + pu[P_AY] * z3 + p[P_BY] * z4))-
						Math.sqrt(1 + one_d2 * (z2 + pl[P_AY] * z3 + p[P_BY] * z4))) / twoDelta;
				value[j][P_BY] = p[P_S0Y] * (
						Math.sqrt(1 + one_d2 * (z2 + p[P_AY] * z3 + pu[P_BY] * z4))-
						Math.sqrt(1 + one_d2 * (z2 + p[P_AY] * z3 + pl[P_BY] * z4))) / twoDelta;
				//@formatter:on
			}

			return value;
		}
	}

	private void record(String name, double[] parameters)
	{
		StringBuilder sb = new StringBuilder(name);
		Rounder rounder = RounderFactory.create(4);
		sb.append(": ").append("gamma=").append(rounder.round(parameters[P_GAMMA]));
		sb.append("; ").append("z0=").append(rounder.round(parameters[P_Z0]));
		sb.append("; ").append("d=").append(rounder.round(parameters[P_D]));
		sb.append("; ").append("s0x=").append(rounder.round(parameters[P_S0X]));
		sb.append("; ").append("Ax=").append(rounder.round(parameters[P_AX]));
		sb.append("; ").append("Bx=").append(rounder.round(parameters[P_BX]));
		sb.append("; ").append("s0y=").append(rounder.round(parameters[P_S0Y]));
		sb.append("; ").append("Ay=").append(rounder.round(parameters[P_AY]));
		sb.append("; ").append("By=").append(rounder.round(parameters[P_BY]));
		IJ.log(sb.toString());
	}

	private void plotFit(final double[] parameters)
	{
		//System.out.println(Arrays.toString(parameters));

		double gamma = parameters[P_GAMMA];
		double z0 = parameters[P_Z0];
		double d = parameters[P_D];
		double s0x = parameters[P_S0X];
		double Ax = parameters[P_AX];
		double Bx = parameters[P_BX];
		double s0y = parameters[P_S0Y];
		double Ay = parameters[P_AY];
		double By = parameters[P_BY];

		// Draw across the entire data range
		double one_d2 = 1.0 / Maths.pow2(d);

		// Update plot
		double[] sx1 = new double[z.length];
		double[] sy1 = new double[z.length];
		for (int i = 0; i < z.length; i++)
		{
			sx1[i] = getS(s0x, z[i] - z0 - gamma, one_d2, Ax, Bx);
			sy1[i] = getS(s0y, z[i] - z0 + gamma, one_d2, Ay, By);
		}

		// Just redraw the plot
		sPlot = plot(null, z, "Width (px)", sx, "Sx", sy, "Sy");
		Plot plot = sPlot.getPlot();
		plot.setColor(Color.RED);
		plot.addPoints(z, sx1, Plot.LINE);
		plot.setColor(Color.BLUE);
		plot.addPoints(z, sy1, Plot.LINE);

		//double[] y = new AstigmatismVectorFunction().value(parameters);		
		//plot.setColor(Color.MAGENTA);
		//plot.addPoints(fitZ, Arrays.copyOf(y, fitZ.length), Plot.BOX);
		//plot.setColor(Color.YELLOW);
		//plot.addPoints(fitZ, Arrays.copyOfRange(y, fitZ.length, y.length), Plot.BOX);

		plot.setColor(Color.BLACK);
		plot.updateImage();
	}
}
