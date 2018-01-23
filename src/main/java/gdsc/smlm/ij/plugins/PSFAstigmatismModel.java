package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.Color;
import java.awt.Label;
import java.awt.Rectangle;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.util.FastMath;

import gdsc.core.ij.IJLogger;
import gdsc.core.ij.IJTrackProgress;
import gdsc.core.ij.Utils;
import gdsc.core.logging.Ticker;
import gdsc.core.utils.ImageExtractor;
import gdsc.core.utils.Maths;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.core.utils.Statistics;
import gdsc.core.utils.TextUtils;
import gdsc.core.utils.TurboList;
import gdsc.smlm.data.config.FitProtos.FitEngineSettings;
import gdsc.smlm.data.config.GUIProtos.PSFAstigmatismModelSettings;
import gdsc.smlm.data.config.GUIProtos.PSFCreatorSettings;
import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.data.config.CalibrationWriter;
import gdsc.smlm.data.config.FitProtosHelper;
import gdsc.smlm.data.config.PSFProtos.ImagePSF;
import gdsc.smlm.data.config.PSFProtos.Offset;
import gdsc.smlm.data.config.PSFProtos.PSFType;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.data.config.PSFProtosHelper;
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
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.FunctionSolver;
import gdsc.smlm.fitting.Gaussian2DFitter;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.ij.IJImageSource;
import gdsc.smlm.ij.settings.ImagePSFHelper;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.model.ImagePSFModel;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.SynchronizedPeakResults;
import gdsc.smlm.results.count.Counter;
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
import ij.gui.Line;
import ij.gui.NonBlockingExtendedGenericDialog;
import ij.gui.Plot;
import ij.gui.Plot2;
import ij.gui.PlotWindow;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
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

		// Fit data
		// fit a z0 point subtracted from the current z

	}

	private boolean loadConfiguration()
	{
		Configuration c = new Configuration();
		// We have a different fit configuration just for the PSF Creator.
		// This allows it to be saved and not effect PeakFit settings.
		config = new FitEngineConfiguration(settings.getFitEngineSettings(), settings.getCalibration(),
				settings.getPsf());
		boolean save = false;
		if (!c.showDialog(config, save))
		{
			IJ.error(TITLE, "No fit configuration loaded");
			return false;
		}

		SettingsManager.writeSettings(settings);

		config = c.getFitEngineConfiguration();
		fitConfig = config.getFitConfiguration();
		if (fitConfig.getPSFType() != PSFType.TWO_AXIS_GAUSSIAN_2D)
		{
			IJ.error(TITLE, "PSF must be " + PSFProtosHelper.getName(PSFType.TWO_AXIS_GAUSSIAN_2D));
			return false;
		}

		config.setIncludeNeighbours(false);
		config.configureOutputUnits();
		config.setResidualsThreshold(1);

		settings.setFitEngineSettings(config.getFitEngineSettings());
		settings.setCalibration(fitConfig.getCalibration());
		settings.setPsf(fitConfig.getPSF());
		SettingsManager.writeSettings(settings);

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

		IJ.log(config.getFitEngineSettings().toString());
		
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
			IJ.log("Fit " + TextUtils.pleural(results.size(), "spot"));

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

		final double[] z = new double[results.size()];
		final double[] x = new double[z.length];
		final double[] y = new double[z.length];
		final double[] I = new double[z.length];
		final double[] sx, sy;
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

		plot(wo, z, "Intensity", x, "Intensity", null, null);
		plot(wo, z, "Position (px)", x, "X", y, "Y");
		plot(wo, z, "Width (px)", sx, "Sx", sy, "Sy");

		wo.tile();

		return true;
	}

	private Plot plot(WindowOrganiser wo, double[] z, String yTitle, double[] y1, String y1Title, double[] y2,
			String y2Title)
	{
		String title = TITLE + " " + yTitle;
		Plot plot = new Plot(title, "Z (μm)", yTitle);
		double[] limits = Maths.limits(y1);
		if (y2 != null)
			limits = Maths.limits(limits, y2);
		plot.setLimits(z[0], z[z.length - 1], limits[0], limits[1]);
		plot.setColor(Color.RED);
		plot.addPoints(z, y1, Plot.CIRCLE);
		if (y2 != null)
		{
			plot.setColor(Color.BLUE);
			plot.addPoints(z, y2, Plot.CIRCLE);
			plot.addLegend(y1Title + "\n" + y2Title);
		}
		Utils.display(title, plot, 0, wo);
		return plot;
	}
}
