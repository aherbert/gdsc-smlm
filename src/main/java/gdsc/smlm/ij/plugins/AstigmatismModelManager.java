package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.Color;
import java.awt.Rectangle;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.analysis.MultivariateMatrixFunction;
import org.apache.commons.math3.analysis.MultivariateVectorFunction;
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer.Optimum;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.util.Precision;

import gdsc.core.data.utils.ConversionException;
import gdsc.core.data.utils.Rounder;
import gdsc.core.data.utils.RounderFactory;
import gdsc.core.ij.IJLogger;
import gdsc.core.ij.IJTrackProgress;
import gdsc.core.ij.Utils;
import gdsc.core.logging.Ticker;
import gdsc.core.utils.Maths;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.core.utils.SimpleLock;
import gdsc.core.utils.TextUtils;
import gdsc.core.utils.TurboList;
import gdsc.smlm.data.config.CalibrationWriter;
import gdsc.smlm.data.config.FitProtos.DataFilterMethod;
import gdsc.smlm.data.config.FitProtos.DataFilterType;
import gdsc.smlm.data.config.GUIProtos.AstigmatismModelManagerSettings;
import gdsc.smlm.data.config.PSFProtos.AstigmatismModel;
import gdsc.smlm.data.config.PSFProtos.AstigmatismModelSettings;
import gdsc.smlm.data.config.PSFProtos.PSFType;
import gdsc.smlm.data.config.PSFProtosHelper;
import gdsc.smlm.data.config.UnitHelper;
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
import gdsc.smlm.function.gaussian.AstigmatismZModel;
import gdsc.smlm.function.gaussian.HoltzerAstigmatismZModel;
import gdsc.smlm.ij.IJImageSource;
import gdsc.smlm.ij.plugins.PeakFit.FitEngineConfigurationProvider;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.model.GaussianPSFModel;
import gdsc.smlm.results.Gaussian2DPeakResultHelper;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.SynchronizedPeakResults;
import gdsc.smlm.results.count.Counter;
import gdsc.smlm.results.procedures.PeakResultProcedure;
import gdsc.smlm.results.procedures.WidthResultProcedure;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.ExtendedGenericDialog;
import ij.gui.GenericDialog;
import ij.gui.Line;
import ij.gui.NonBlockingExtendedGenericDialog;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.plugin.WindowOrganiser;
import ij.process.FloatPolygon;
import ij.process.FloatProcessor;
import ij.text.TextWindow;

/**
 * Produces a 2D Gaussian astigmatism model for a 2D astigmatic PSF.
 * <p>
 * The input images must be a z-stack of a PSF.
 */
public class AstigmatismModelManager implements PlugIn
{
	private final static String TITLE = "Astigmatism Model Manager";

	private static AstigmatismModelSettings.Builder settings = null;
	private static TextWindow resultsWindow = null;

	private AstigmatismModelManagerSettings.Builder pluginSettings;
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

	private static AstigmatismModelSettings.Builder getSettings()
	{
		return getSettings(0);
	}

	private static AstigmatismModelSettings.Builder getSettings(int flags)
	{
		if (settings == null)
			settings = SettingsManager.readAstigmatismModelSettings(flags).toBuilder();
		return settings;
	}

	/**
	 * List the astigmatism models.
	 *
	 * @param includeNone
	 *            Set to true to include an invalid none model string
	 * @return the list
	 */
	public static String[] listAstigmatismModels(boolean includeNone)
	{
		AstigmatismModelSettings.Builder settings = getSettings();
		List<String> list = createList(includeNone);
		list.addAll(settings.getAstigmatismModelResourcesMap().keySet());
		return list.toArray(new String[list.size()]);
	}

	private static List<String> createList(boolean includeNone)
	{
		List<String> list = new TurboList<String>();
		if (includeNone)
			list.add("[None]");
		return list;
	}

	/**
	 * List the astigmatism models with the correct pixel scale.
	 *
	 * @param includeNone
	 *            Set to true to include an empty string
	 * @param nmPerPixel
	 *            the nm per pixel
	 * @return the list
	 */
	public static String[] listAstigmatismModels(boolean includeNone, double nmPerPixel)
	{
		return listAstigmatismModels(includeNone, nmPerPixel, 0);
	}

	/**
	 * List the astigmatism models with the correct pixel scale within an error margin.
	 *
	 * @param includeNone
	 *            Set to true to include an empty string
	 * @param nmPerPixel
	 *            the nm per pixel
	 * @param error
	 *            the error margin
	 * @return the list
	 */
	public static String[] listAstigmatismModels(boolean includeNone, double nmPerPixel, double error)
	{
		AstigmatismModelSettings.Builder settings = getSettings();
		List<String> list = createList(includeNone);
		error = Math.abs(error);
		double low = nmPerPixel - error;
		double high = nmPerPixel + error;
		for (Map.Entry<String, AstigmatismModel> entry : settings.getAstigmatismModelResourcesMap().entrySet())
		{
			AstigmatismModel resource = entry.getValue();
			if (resource.getNmPerPixel() >= low && resource.getNmPerPixel() <= high)
				list.add(entry.getKey());
		}
		return list.toArray(new String[list.size()]);
	}

	/**
	 * List the astigmatism models with their pixel scale.
	 *
	 * @param includeNone
	 *            Set to true to include an empty string
	 * @param withNmPerPixel
	 *            Append the nm per pixel to the model names
	 * @return the list
	 */
	public static String[] listAstigmatismModels(boolean includeNone, boolean withNmPerPixel)
	{
		AstigmatismModelSettings.Builder settings = getSettings();
		List<String> list = createList(includeNone);
		for (Map.Entry<String, AstigmatismModel> entry : settings.getAstigmatismModelResourcesMap().entrySet())
		{
			AstigmatismModel resource = entry.getValue();
			if (withNmPerPixel)
				list.add(String.format("%s [%s nm]", entry.getKey(), Utils.rounded(resource.getNmPerPixel())));
			else
				list.add(entry.getKey());
		}
		return list.toArray(new String[list.size()]);
	}

	/**
	 * Remove the extra information added to a name for use in dialogs
	 * 
	 * @param name
	 *            the formatted name
	 * @return The name
	 */
	public static String removeFormatting(String name)
	{
		int index = name.lastIndexOf('[');
		if (index > 0)
			name = name.substring(0, index - 1);
		return name;
	}

	/**
	 * Gets the model.
	 *
	 * @param name
	 *            the name
	 * @return the model (or null)
	 */
	public static AstigmatismModel getModel(String name)
	{
		return getSettings().getAstigmatismModelResourcesMap().get(name);
	}

	/**
	 * Creates the working astigmatism model from the model settings.
	 *
	 * @param model
	 *            the model settings
	 * @return the astigmatism Z model
	 */
	public static AstigmatismZModel create(AstigmatismModel model)
	{
		return HoltzerAstigmatismZModel.create(model.getS0X(), model.getS0Y(), model.getGamma(), model.getD(),
				model.getAx(), model.getBx(), model.getAy(), model.getBy());
	}

	//@formatter:off
	private static String[] OPTIONS = { 
			"Create model",
			"Import model",
			// All option below require models
			"View model",
			"Delete model",
			"Invert model",
			"Export model",
			};
	//@formatter:on
	private static String[] OPTIONS2;
	static
	{
		OPTIONS2 = Arrays.copyOf(OPTIONS, 1);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		String[] options = OPTIONS;
		AstigmatismModelSettings.Builder settings = getSettings(SettingsManager.FLAG_SILENT);
		if (settings.getAstigmatismModelResourcesCount() == 0)
		{
			options = OPTIONS2;
		}

		pluginSettings = SettingsManager.readAstigmatismModelManagerSettings(0).toBuilder();

		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addChoice("Option", options, pluginSettings.getOption());
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		pluginSettings.setOption(gd.getNextChoiceIndex());

		switch (pluginSettings.getOption())
		{
			case 5:
				exportModel();
				break;
			case 4:
				invertModel();
				break;
			case 3:
				deleteModel();
				break;
			case 2:
				viewModel();
				break;
			case 1:
				importModel();
				break;
			default:
				createModel();
		}

		SettingsManager.writeSettings(pluginSettings);
	}

	private void createModel()
	{
		if (!getImage())
			return;

		if (!showFitDialog())
			return;

		if (!loadConfiguration())
			return;

		// TODO - Can this plugin be updated to handle multiple spots on the same image?
		// This would produce an astigmatism curve for each spot. Then optimise the parameters
		// by having a z0 for each spot.
		// - Store the width curve (z, sx, sy) for each spot.
		// - Pick the range (fitZ, fitSx, fitSy) for each spot.
		// - Guess parameters for each and then average them.
		// - Optimisation function to handle the different z-range for each spot.

		if (!findFitRegion())
			return;

		if (!fitRegion())
			return;

		if (!plotData())
			return;

		if (!fitData())
			return;

		saveModel();
	}

	private boolean getImage()
	{
		// Select an image
		GenericDialog gd = new GenericDialog(TITLE);
		String[] list = getImageList();
		if (list.length == 0)
		{
			IJ.error("No suitable images");
			return false;
		}
		gd.addChoice("Image", list, pluginSettings.getImage());
		gd.showDialog();
		if (gd.wasCanceled())
			return false;
		String image = gd.getNextChoice();
		pluginSettings.setImage(image);
		imp = WindowManager.getImage(image);
		if (imp == null)
		{
			IJ.error(TITLE, "Failed to find image: " + image);
			return false;
		}
		Roi roi = imp.getRoi();
		if (roi == null || roi.getType() != Roi.POINT)
		{
			IJ.error("Point ROI required");
			return false;
		}
		return true;
	}

	private static String[] getImageList()
	{
		TurboList<String> newImageList = new TurboList<String>();

		for (int id : Utils.getIDList())
		{
			ImagePlus imp = WindowManager.getImage(id);
			if (imp == null)
				continue;
			if (imp.getNDimensions() != 3)
				continue;
			if (imp.getBitDepth() == 24)
				continue;
			Roi roi = imp.getRoi();
			if (roi == null || roi.getType() != Roi.POINT)
				continue;
			newImageList.add(imp.getTitle());
		}

		return newImageList.toArray(new String[0]);
	}

	private boolean showFitDialog()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		pluginSettings = SettingsManager.readAstigmatismModelManagerSettings(0).toBuilder();

		guessScale();

		gd.addMessage("Use Gaussian 2D PSF fitting to create an astigmatism z-model");

		gd.addNumericField("nm_per_slice", pluginSettings.getNmPerSlice(), 0);

		gd.showDialog();

		SettingsManager.writeSettings(pluginSettings);

		if (gd.wasCanceled())
			return false;

		pluginSettings.setNmPerSlice(gd.getNextNumber());

		// Check arguments
		try
		{
			Parameters.isPositive("nm/slice", pluginSettings.getNmPerSlice());
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		return true;
	}

	private void guessScale()
	{
		CalibrationWriter cw = CalibrationWriter.create(pluginSettings.getCalibration());
		// It does not matter if we already have settings, try and update them anyway
		Calibration c = imp.getCalibration();
		double r = guessScale(c.getXUnit(), c.pixelWidth);
		if (r != 0)
		{
			cw.setNmPerPixel(r);
			pluginSettings.setCalibration(cw.getBuilder());
		}
		r = guessScale(c.getZUnit(), c.pixelDepth);
		if (r != 0)
			pluginSettings.setNmPerSlice(r);
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

	private boolean loadConfiguration()
	{
		// We have a different fit configuration just for the PSF Creator.
		// This allows it to be saved and not effect PeakFit settings.
		config = new FitEngineConfiguration(pluginSettings.getFitEngineSettings(), pluginSettings.getCalibration(),
				pluginSettings.getPsf());
		if (!showConfigurationDialog())
		{
			IJ.error(TITLE, "No fit configuration loaded");
			return false;
		}

		SettingsManager.writeSettings(pluginSettings);

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

		pluginSettings.setFitEngineSettings(config.getFitEngineSettings());
		pluginSettings.setCalibration(fitConfig.getCalibration());
		pluginSettings.setPsf(fitConfig.getPSF());
		SettingsManager.writeSettings(pluginSettings);

		return true;
	}

	public boolean showConfigurationDialog()
	{
		fitConfig = config.getFitConfiguration();
		CalibrationWriter calibration = fitConfig.getCalibrationWriter();

		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);
		gd.addMessage("Configuration settings for the single-molecule localisation microscopy plugins");

		PeakFit.addCameraOptions(gd, fitConfig);
		gd.addNumericField("Calibration (nm/px)", calibration.getNmPerPixel(), 2);
		//gd.addNumericField("Exposure_time (ms)", calibration.getExposureTime(), 2);

		PeakFit.addPSFOptions(gd, fitConfig);

		FitEngineConfigurationProvider provider = new PeakFit.SimpleFitEngineConfigurationProvider(config);
		PeakFit.addFittingOptions(gd, provider);
		gd.addChoice("Fit_solver", SettingsManager.getFitSolverNames(), fitConfig.getFitSolver().ordinal());
		gd.addCheckbox("Log_fit_progress", pluginSettings.getLogFitProgress());

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
		fitConfig.setCalibration(calibration.getCalibration());
		fitConfig.setPSFType(PeakFit.getPSFTypeValues()[gd.getNextChoiceIndex()]);
		config.setFitting(gd.getNextNumber());
		fitConfig.setFitSolver(gd.getNextChoiceIndex());
		pluginSettings.setLogFitProgress(gd.getNextBoolean());
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

		if (pluginSettings.getLogFitProgress())
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

		final double umPerSlice = pluginSettings.getNmPerSlice() / 1000.0;
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
			plot.setColor(Color.BLACK);
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

		NonBlockingExtendedGenericDialog gd = new NonBlockingExtendedGenericDialog(TITLE);
		gd.addMessage("Select z-range for curve fit.\nChoose a region with a smooth width curve and low XY drift.");
		gd.addSlider("Min_z", minz, maxz, minz);
		gd.addSlider("Max_z", minz, maxz, maxz);
		gd.addMessage("Curve parameter estimation");
		gd.addSlider("Smoothing", 0.05, 0.5, pluginSettings.getSmoothing());
		gd.addCheckbox("Show_estimated_curve", pluginSettings.getShowEstimatedCurve());
		gd.addMessage("Fit options");
		gd.addCheckbox("Weighted_fit", pluginSettings.getWeightedFit());
		gd.addDialogListener(new ZDialogListener());
		gd.showDialog();

		// Save settings
		SettingsManager.writeSettings(pluginSettings);

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

		public ZDialogListener()
		{
			sPlot.getImagePlus().killRoi();
			xyPlot.getImagePlus().killRoi();
		}

		public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
		{
			int oldMinz = minz;
			int oldMaxz = maxz;
			minz = (int) gd.getNextNumber();
			maxz = (int) gd.getNextNumber();
			pluginSettings.setSmoothing(gd.getNextNumber());
			pluginSettings.setShowEstimatedCurve(gd.getNextBoolean());
			pluginSettings.setWeightedFit(gd.getNextBoolean());
			if (showRoi && (oldMinz != minz || oldMaxz != maxz))
			{
				addRoi(sPlot);
				addRoi(xyPlot);
			}
			return maxz > minz;
		}

		private void addRoi(PlotWindow pw)
		{
			Plot plot = pw.getPlot();
			int x1 = (int) plot.scaleXtoPxl(z[minz]);
			int x2 = (int) plot.scaleXtoPxl(z[maxz]);
			double[] limits = plot.getLimits();
			int y1 = (int) plot.scaleYtoPxl(limits[3]);
			int y2 = (int) plot.scaleYtoPxl(limits[2]);
			pw.getImagePlus().setRoi(new Roi(x1, y1, x2 - x1, y2 - y1));
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
		if (pluginSettings.getSmoothing() > 0)
		{
			LoessInterpolator loess = new LoessInterpolator(pluginSettings.getSmoothing(), 0);
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
		// If this is not the case we can invert the gamma parameter.
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

		record("Initial", parameters);
		if (pluginSettings.getShowEstimatedCurve())
		{
			plotFit(parameters);
			IJ.showMessage(TITLE, "Showing the estimated curve parameters.\nClick OK to continue.");
		}

		//@formatter:off
		LeastSquaresBuilder builder = new LeastSquaresBuilder()
				.maxEvaluations(Integer.MAX_VALUE)
				.maxIterations(3000)
				.start(parameters)
				.target(y);
		//@formatter:on

		if (pluginSettings.getWeightedFit())
			builder.weight(new DiagonalMatrix(getWeights(smoothSx, smoothSy)));

		AstigmatismVectorFunction vf = new AstigmatismVectorFunction();
		builder.model(vf, new AstigmatismMatrixFunction());

		LeastSquaresProblem problem = builder.build();

		try
		{
			Optimum optimum = optimizer.optimize(problem);

			parameters = optimum.getPoint().toArray();

			record("Final", parameters);
			plotFit(parameters);

			saveResult(optimum);
		}
		catch (Exception e)
		{
			IJ.error(TITLE, "Failed to fit curve: " + e.getMessage());
			return false;
		}

		return true;
	}

	/**
	 * Get depth of focus as the point where width = min width * sqrt(2).
	 *
	 * @param min
	 *            the index of the min value of the width
	 * @param z
	 *            the z
	 * @param sx
	 *            the width
	 * @return the estimated depth of focus
	 */
	private static double estimateD(int min, double[] z, double[] sx)
	{
		// w = w0 * sqrt(1 + z^2/d^2)
		// if z==d then w = w0 * sqrt(2) 

		double w = sx[min] * 1.414213562; // sqrt(2) the min width
		int lower = min;
		while (lower > 0 && sx[lower] < w)
			lower--;
		int upper = min;
		while (upper < sx.length - 1 && sx[upper] < w)
			upper++;
		return (z[upper] - z[lower]) / 2; // Since we searched both directions
	}

	private static double[] getWeights(double[]... y)
	{
		int n = 0;
		for (int i = 0; i < y.length; i++)
			n += y[i].length;
		double[] w = new double[n];
		for (int i = 0, k = 0; i < y.length; i++)
			for (int j = 0; j < y[i].length; j++)
				w[k++] = 1.0 / y[i][j];
		return w;
	}

	private static final int P_GAMMA = 0;
	private static final int P_D = 1;
	private static final int P_S0X = 2;
	private static final int P_AX = 3;
	private static final int P_BX = 4;
	private static final int P_S0Y = 5;
	private static final int P_AY = 6;
	private static final int P_BY = 7;
	private static final int P_Z0 = 8;

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
		sb.append("; ").append("d=").append(rounder.round(parameters[P_D]));
		sb.append("; ").append("s0x=").append(rounder.round(parameters[P_S0X]));
		sb.append("; ").append("Ax=").append(rounder.round(parameters[P_AX]));
		sb.append("; ").append("Bx=").append(rounder.round(parameters[P_BX]));
		sb.append("; ").append("s0y=").append(rounder.round(parameters[P_S0Y]));
		sb.append("; ").append("Ay=").append(rounder.round(parameters[P_AY]));
		sb.append("; ").append("By=").append(rounder.round(parameters[P_BY]));
		sb.append("; ").append("z0=").append(rounder.round(parameters[P_Z0]));
		IJ.log(sb.toString());
	}

	private void plotFit(final double[] parameters)
	{
		//System.out.println(Arrays.toString(parameters));

		double gamma = parameters[P_GAMMA];
		double d = parameters[P_D];
		double s0x = parameters[P_S0X];
		double Ax = parameters[P_AX];
		double Bx = parameters[P_BX];
		double s0y = parameters[P_S0Y];
		double Ay = parameters[P_AY];
		double By = parameters[P_BY];
		double z0 = parameters[P_Z0];

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

	private void saveResult(Optimum optimum)
	{
		createResultWindow();
		StringBuilder sb = new StringBuilder();
		Rounder rounder = RounderFactory.create(4);
		sb.append(fitZ.length * 2);
		sb.append('\t').append(pluginSettings.getWeightedFit());
		sb.append('\t').append(Utils.rounded(optimum.getRMS(), 6));
		sb.append('\t').append(optimum.getIterations());
		sb.append('\t').append(optimum.getEvaluations());
		sb.append('\t').append(rounder.round(parameters[P_GAMMA]));
		sb.append('\t').append(rounder.round(parameters[P_D]));
		sb.append('\t').append(rounder.round(parameters[P_S0X]));
		sb.append('\t').append(rounder.round(parameters[P_AX]));
		sb.append('\t').append(rounder.round(parameters[P_BX]));
		sb.append('\t').append(rounder.round(parameters[P_S0Y]));
		sb.append('\t').append(rounder.round(parameters[P_AY]));
		sb.append('\t').append(rounder.round(parameters[P_BY]));
		sb.append('\t').append(rounder.round(parameters[P_Z0]));
		resultsWindow.append(sb.toString());
	}

	private void createResultWindow()
	{
		if (resultsWindow == null || !resultsWindow.isShowing())
			resultsWindow = new TextWindow(TITLE,
					"N\tWeighted\tRMS\tIter\tEval\tgamma\td\ts0x\tAx\tBx\ts0y\tAy\tBy\tz0", "", 1000, 300);
	}

	private boolean saveModel()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addCheckbox("Save_model", pluginSettings.getSaveModel());
		gd.addStringField("Model_name", pluginSettings.getModelName());
		gd.addCheckbox("Save_fit_width", pluginSettings.getSaveFitWidth());
		//gd.setCancelLabel(" No ");
		gd.showDialog();
		if (gd.wasCanceled())
			return false;
		pluginSettings.setSaveModel(gd.getNextBoolean());
		String name = gd.getNextString();
		pluginSettings.setSaveFitWidth(gd.getNextBoolean());

		if (pluginSettings.getSaveFitWidth())
		{
			// Save the widths in the fit configuration
			fitConfig.setInitialPeakStdDev0(parameters[P_S0X]);
			fitConfig.setInitialPeakStdDev1(parameters[P_S0Y]);
			pluginSettings.setPsf(fitConfig.getPSF());
			SettingsManager.writeSettings(pluginSettings);
		}

		if (pluginSettings.getSaveModel())
		{
			AstigmatismModel.Builder model = AstigmatismModel.newBuilder();
			model.setGamma(parameters[P_GAMMA]);
			model.setD(parameters[P_D]);
			model.setS0X(parameters[P_S0X]);
			model.setAx(parameters[P_AX]);
			model.setBx(parameters[P_BX]);
			model.setS0Y(parameters[P_S0Y]);
			model.setAy(parameters[P_AY]);
			model.setBy(parameters[P_BY]);
			model.setZDistanceUnit(DistanceUnit.UM);
			model.setSDistanceUnit(DistanceUnit.PIXEL);
			model.setNmPerPixel(fitConfig.getCalibrationReader().getNmPerPixel());
			return save(name, model);
		}
		return true;
	}

	private boolean save(String name, AstigmatismModel.Builder model)
	{
		pluginSettings.setModelName(name);

		// Check existing names
		AstigmatismModelSettings.Builder settings = getSettings();
		Map<String, AstigmatismModel> map = settings.getAstigmatismModelResourcesMap();
		if (map.containsKey(name))
		{
			name = suggest(map, name);
			ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
			gd.addMessage("Model name " + pluginSettings.getModelName() + " already exists.\n \nSuggest renaming to:");
			gd.addStringField("Model_name", name);
			gd.enableYesNoCancel("Rename", "Overwrite");
			gd.showDialog(true);
			if (gd.wasCanceled())
				return false;
			if (gd.wasOKed())
				// Rename
				pluginSettings.setModelName(name);
		}

		// Save the model
		settings.putAstigmatismModelResources(pluginSettings.getModelName(), model.build());
		if (!SettingsManager.writeSettings(settings.build()))
		{
			IJ.error(TITLE, "Failed to save the model");
			return false;
		}
		return true;
	}

	private String suggest(Map<String, AstigmatismModel> map, String name)
	{
		name += '_';
		for (int i = 2; i > 0; i++)
		{
			String name2 = name + i;
			if (!map.containsKey(name2))
				return name2;
		}
		return ""; // This happens if there are a lot of models
	}

	private void importModel()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addStringField("Model_name", pluginSettings.getModelName());
		gd.addFilenameField("Filename", pluginSettings.getFilename());
		//gd.setCancelLabel(" No ");
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		String name = gd.getNextString();
		pluginSettings.setModelName(name);
		pluginSettings.setFilename(gd.getNextString());

		AstigmatismModel.Builder model = AstigmatismModel.newBuilder();
		if (!SettingsManager.fromJSON(new File(pluginSettings.getFilename()), model, 0))
			return;

		save(name, model);
	}

	private void exportModel()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		String[] MODELS = listAstigmatismModels(false);
		gd.addChoice("Model", MODELS, pluginSettings.getSelected());
		gd.addFilenameField("Filename", pluginSettings.getFilename());
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		String name = gd.getNextChoice();
		pluginSettings.setSelected(name);
		pluginSettings.setFilename(gd.getNextString());

		// Try and get the named resource
		AstigmatismModel model = settings.getAstigmatismModelResourcesMap().get(name);
		if (model == null)
		{
			IJ.error(TITLE, "Failed to find astigmatism model: " + name);
			return;
		}

		if (!SettingsManager.toJSON(model, new File(pluginSettings.getFilename()), 0))
		{
			IJ.error(TITLE, "Failed to export astigmatism model: " + name);
			return;
		}

		Utils.log("Exported astigmatism model: %s to %s", name, pluginSettings.getFilename());
	}

	private void viewModel()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		String[] MODELS = listAstigmatismModels(false);
		gd.addChoice("Model", MODELS, pluginSettings.getSelected());
		gd.addChoice("z_distance_unit", SettingsManager.getDistanceUnitNames(), pluginSettings.getZDistanceUnitValue());
		gd.addChoice("s_distance_unit", SettingsManager.getDistanceUnitNames(), pluginSettings.getSDistanceUnitValue());
		gd.addCheckbox("Show_depth_of_focus", pluginSettings.getShowDepthOfFocus());
		gd.addCheckbox("Show_combined_width", pluginSettings.getShowCombinedWidth());
		gd.addCheckbox("Show_PSF", pluginSettings.getShowPsf());
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		String name = gd.getNextChoice();
		pluginSettings.setSelected(name);
		pluginSettings.setZDistanceUnitValue(gd.getNextChoiceIndex());
		pluginSettings.setSDistanceUnitValue(gd.getNextChoiceIndex());
		pluginSettings.setShowDepthOfFocus(gd.getNextBoolean());
		pluginSettings.setShowCombinedWidth(gd.getNextBoolean());
		pluginSettings.setShowPsf(gd.getNextBoolean());

		// Try and get the named resource
		AstigmatismModel model = settings.getAstigmatismModelResourcesMap().get(name);
		if (model == null)
		{
			IJ.error(TITLE, "Failed to find astigmatism model: " + name);
			return;
		}

		try
		{
			model = convert(model, pluginSettings.getZDistanceUnit(), pluginSettings.getSDistanceUnit());
		}
		catch (ConversionException e)
		{
			Utils.log("Bad conversion (%s), defaulting to native model units", e.getMessage());
		}

		Utils.log("Astigmatism model: %s\n%s", name, model);

		// Plot the curve. Do this so we encompass twice the depth-of-field.
		double gamma = model.getGamma();
		double d = model.getD();
		double s0x = model.getS0X();
		double Ax = model.getAx();
		double Bx = model.getBx();
		double s0y = model.getS0Y();
		double Ay = model.getAy();
		double By = model.getBy();

		double range = Math.abs(gamma) + 1.5 * d;
		int n = 200;
		double step = range / n;
		double[] z = new double[2 * n + 1];
		double[] sx = new double[z.length];
		double[] sy = new double[z.length];
		// Use the same class that is used during fitting
		HoltzerAstigmatismZModel m = HoltzerAstigmatismZModel.create(s0x, s0y, gamma, d, Ax, Bx, Ay, By);
		for (int i = 0; i < z.length; i++)
		{
			double zz = -range + i * step;
			z[i] = zz;
			sx[i] = m.getSx(zz);
			sy[i] = m.getSy(zz);
		}

		String title = TITLE + " Width Curve";
		Plot plot = new Plot(title, "Z (" + UnitHelper.getShortName(model.getZDistanceUnit()) + ")",
				"Width (" + UnitHelper.getShortName(model.getSDistanceUnit()) + ")");
		double[] limits = Maths.limits(sx);
		limits = Maths.limits(limits, sy);
		double rangex = (z[z.length - 1] - z[0]) * 0.05;
		double rangey = (limits[1] - limits[0]) * 0.05;
		double miny = limits[0] - rangey;
		double maxy = limits[1] + rangey;
		plot.setLimits(z[0] - rangex, z[z.length - 1] + rangex, miny, maxy);
		plot.setColor(Color.RED);
		plot.addPoints(z, sx, Plot.LINE);
		plot.setColor(Color.BLUE);
		plot.addPoints(z, sy, Plot.LINE);
		plot.setColor(Color.YELLOW);

		if (pluginSettings.getShowDepthOfFocus())
		{
			double z0x = gamma, z0y = -gamma;
			plot.setColor(Color.RED.darker());
			plot.drawDottedLine(z0x - d, miny, z0x - d, maxy, 4);
			plot.drawDottedLine(z0x + d, miny, z0x + d, maxy, 4);
			plot.setColor(Color.BLUE.darker());
			plot.drawDottedLine(z0y - d, miny, z0y - d, maxy, 4);
			plot.drawDottedLine(z0y + d, miny, z0y + d, maxy, 4);
		}
		String legend = "Sx\nSy";
		if (pluginSettings.getShowCombinedWidth())
		{
			double [] s = new double[z.length];
			for (int i = 0; i < z.length; i++)
			{
				s[i] = Gaussian2DPeakResultHelper.getStandardDeviation(sx[i], sy[i]);
			}
			plot.setColor(Color.GREEN);
			plot.addPoints(z, s, Plot.LINE);
			legend += "\tS";
		}

		plot.setColor(Color.BLACK);
		plot.addLegend(legend);
		plot.addLabel(0, 0, String.format("Model = %s (%s nm/pixel)", name, Utils.rounded(model.getNmPerPixel())));
		Utils.display(title, plot);

		if (!pluginSettings.getShowPsf())
			return;

		// Get pixel range using 3x[max SD]
		int width = 1 + 2 * ((int) Math.ceil(limits[1] * 3));
		new ModelRenderer(name, model, m, range, width, plot).run();
	}

	private class ModelRenderer implements DialogListener
	{
		String name;
		AstigmatismModel model;
		double range;
		double z = 0;
		int width;
		double cx;
		GaussianPSFModel psf;
		Plot plot;

		public ModelRenderer(String name, AstigmatismModel model, HoltzerAstigmatismZModel m, double range, int width,
				Plot plot)
		{
			this.name = name;
			this.model = model;
			this.range = range;
			this.width = width;
			cx = width * 0.5;
			psf = new GaussianPSFModel(m);
			this.plot = plot;
		}

		public void run()
		{
			final NonBlockingExtendedGenericDialog gd = new NonBlockingExtendedGenericDialog(TITLE);
			gd.addMessage("Model = " + name);
			gd.addSlider("Z (" + UnitHelper.getShortName(model.getZDistanceUnit()) + ")", -range, range, 0);
			final TextField tfz = gd.getLastTextField();
			gd.addCheckbox("Calibrated_image", pluginSettings.getCalibratedImage());
			gd.addDialogListener(this);
			if (Utils.isShowGenericDialog())
			{
				gd.hideCancelButton();
				gd.setOKLabel(" Close ");
				gd.addAndGetButton("Reset", new ActionListener()
				{
					public void actionPerformed(ActionEvent e)
					{
						update();
						// The events triggered by setting these should be ignored now
						tfz.setText("0");
					}
				});
				draw();
			}
			gd.showDialog();

			SettingsManager.writeSettings(pluginSettings);
		}

		public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
		{
			z = gd.getNextNumber();
			pluginSettings.setCalibratedImage(gd.getNextBoolean());
			if (gd.invalidNumber())
				return false;
			update();
			return true;
		}

		private double _z = -1;
		private int _calibratedImage = -1;

		private void draw()
		{
			_z = z;
			int calibratedImage = getCalibratedImage();
			float[] data = new float[width * width];
			psf.create3D(data, width, width, 1, cx, cx, _z, false);
			ImagePlus imp = Utils.display(TITLE + " PSF", new FloatProcessor(width, width, data));
			if (_calibratedImage != calibratedImage)
			{
				if (calibratedImage == 1)
				{
					Calibration cal = new Calibration();
					cal.setXUnit("um");
					cal.pixelWidth = cal.pixelHeight = model.getNmPerPixel() / 1000;
					imp.setCalibration(cal);
				}
				else
				{
					imp.setCalibration(null);
				}
				_calibratedImage = calibratedImage;
			}
			imp.resetDisplayRange();

			// Show the z position using an overlay
			double x = plot.scaleXtoPxl(_z);
			double[] limits = plot.getLimits();
			int y1 = (int) plot.scaleYtoPxl(limits[3]);
			int y2 = (int) plot.scaleYtoPxl(limits[2]);
			plot.getImagePlus().setRoi(new Line(x, y1, x, y2));
		}

		SimpleLock lock = new SimpleLock();

		private void update()
		{
			if (lock.acquire())
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
							while (
									_z != z ||
									_calibratedImage != getCalibratedImage()
									)
							//@formatter:on
							{
								draw();
							}
						}
						finally
						{
							// Ensure the running flag is reset
							lock.release();
						}
					}
				}).start();
			}
		}

		private int getCalibratedImage()
		{
			return pluginSettings.getCalibratedImage() ? 1 : 0;
		}
	}

	/**
	 * Convert the model to the given units.
	 *
	 * @param model
	 *            the model
	 * @param zDistanceUnit
	 *            the desired input z distance unit
	 * @param sDistanceUnit
	 *            the desired output s distance unit
	 * @return the astigmatism model
	 * @throws ConversionException
	 *             if the units cannot be converted
	 */
	public static AstigmatismModel convert(AstigmatismModel model, DistanceUnit zDistanceUnit,
			DistanceUnit sDistanceUnit) throws ConversionException
	{
		return PSFProtosHelper.convert(model, zDistanceUnit, sDistanceUnit);
	}

	private void deleteModel()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		String[] MODELS = listAstigmatismModels(false);
		gd.addChoice("Model", MODELS, pluginSettings.getSelected());
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		String name = gd.getNextChoice();
		pluginSettings.setSelected(name);

		AstigmatismModel model = settings.getAstigmatismModelResourcesMap().get(name);
		if (model == null)
		{
			IJ.error(TITLE, "Failed to find astigmatism model: " + name);
			return;
		}

		settings.removeAstigmatismModelResources(name);
		SettingsManager.writeSettings(settings.build());

		Utils.log("Deleted astigmatism model: %s", name);
	}

	private void invertModel()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		String[] MODELS = listAstigmatismModels(false);
		gd.addMessage("Invert the z-orientation of a model.\n \n" + TextUtils.wrap(
				//@formatter:off
				 "Note that a positive gamma puts the focal plane for the X-dimension " +
				 "above the z-centre (positive Z) and the focal "+
				 "plane for the Y-dimension below the z-centre (negative Z). If gamma " +
				 "is negative then the orientation of the focal "+ 
				 "planes of X and Y are reversed."
				//@formatter:on
				, 80));
		gd.addChoice("Model", MODELS, pluginSettings.getSelected());
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		String name = gd.getNextChoice();
		pluginSettings.setSelected(name);

		AstigmatismModel model = settings.getAstigmatismModelResourcesMap().get(name);
		if (model == null)
		{
			IJ.error(TITLE, "Failed to find astigmatism model: " + name);
			return;
		}

		AstigmatismModel.Builder builder = model.toBuilder();
		// Invert the gamma
		builder.setGamma(-model.getGamma());
		// Invert the constants of z^3 as these have an asymmetric effect on the curve 
		builder.setAx(-model.getAx());
		builder.setAy(-model.getAy());
		// The constants of z^4 have a symmetric effect on the curve
		//builder.setBx(-model.getBx());
		//builder.setBy(-model.getBy());
		settings.putAstigmatismModelResources(name, builder.build());
		SettingsManager.writeSettings(settings.build());

		Utils.log("Inverted astigmatism model: %s", name);
	}
}
