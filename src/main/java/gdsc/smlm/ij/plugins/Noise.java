package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

import gdsc.core.ij.Utils;
import gdsc.core.utils.NoiseEstimator;
import gdsc.core.utils.Statistics;
import gdsc.smlm.data.config.CalibrationWriter;
import gdsc.smlm.data.config.FitProtos.NoiseEstimatorMethod;
import gdsc.smlm.data.config.FitProtosHelper;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.ImageConverter;
import gdsc.smlm.model.camera.CameraModel;
import gdsc.smlm.model.camera.FixedPixelCameraModel;
import gdsc.smlm.model.camera.NullCameraModel;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.DialogListener;
import ij.gui.ExtendedGenericDialog;
import ij.gui.GenericDialog;
import ij.gui.Plot2;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.ImageProcessor;
import ij.text.TextWindow;
import ij.util.Tools;

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

/**
 * Contains methods to find the noise in the provided image data.
 */
public class Noise implements ExtendedPlugInFilter, DialogListener
{
	private static final String TITLE = "Noise Estimator";
	private List<double[]> results;
	private final int FLAGS = DOES_8G | DOES_16 | DOES_32 | PARALLELIZE_STACKS | FINAL_PROCESSING | NO_CHANGES;
	private PlugInFilterRunner pfr;
	private ImagePlus imp;
	private static int algorithm = NoiseEstimatorMethod.ALL_PIXELS_VALUE;
	private static int algorithm2 = NoiseEstimatorMethod.QUICK_RESIDUALS_LEAST_TRIMMED_OF_SQUARES_VALUE;
	private static int lowestPixelsRange = 6;
	private static int ox = 0, oy = 0;
	private CalibrationWriter calibration;
	private CameraModel cameraModel;

	private static final String Y_AXIS_COUNT = "Noise (counts)";
	private static final String Y_AXIS_PHOTON = "Noise (photons)";
	private String yAxisTitle;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	public int setup(String arg, ImagePlus imp)
	{
		if (arg.equalsIgnoreCase("final"))
		{
			showResults();
			return DONE;
		}
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (imp == null)
		{
			IJ.noImage();
			return DONE;
		}
		this.imp = imp;
		results = Collections.synchronizedList(new ArrayList<double[]>(imp.getStackSize()));
		return FLAGS;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.ExtendedPlugInFilter#showDialog(ij.ImagePlus, java.lang.String,
	 * ij.plugin.filter.PlugInFilterRunner)
	 */
	public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr)
	{
		// Select a camera model
		calibration = CalibrationWriter.create(SettingsManager.readCalibration(0));
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addMessage("Preprocess camera image and estimate global noise");
		PeakFit.addCameraOptions(gd, calibration);
		gd.showDialog();
		if (gd.wasCanceled())
			return DONE;
		calibration.setCameraType(SettingsManager.getCameraTypeValues()[gd.getNextChoiceIndex()]);
		gd.collectOptions();
		SettingsManager.writeSettings(calibration.getCalibration());

		try
		{
			createCameraModel();
		}
		catch (Exception e)
		{
			IJ.error(TITLE, e.getMessage());
			return DONE;
		}

		// If using a stack, provide a preview graph of the noise for two methods
		if (imp.getStackSize() > 1)
		{
			this.pfr = pfr;

			drawPlot();

			gd = new ExtendedGenericDialog(TITLE);
			gd.addHelp(About.HELP_URL);

			String[] methodNames = SettingsManager.getNoiseEstimatorMethodNames();

			gd.addChoice("Method1 (blue)", methodNames, methodNames[algorithm]);
			gd.addChoice("Method2 (red)", methodNames, methodNames[algorithm2]);
			gd.addSlider("Lowest_radius", 1, 15, lowestPixelsRange);

			//gd.addPreviewCheckbox(pfr);
			gd.addDialogListener(this);
			gd.addMessage("Click OK to compute noise table using all methods");
			gd.showDialog();

			if (gd.wasCanceled() || !dialogItemChanged(gd, null))
				return DONE;
		}

		return IJ.setupDialog(imp, FLAGS);
	}

	private void createCameraModel()
	{
		yAxisTitle = Y_AXIS_PHOTON;
		switch (calibration.getCameraType())
		{
			case CCD:
			case EMCCD:
				cameraModel = new FixedPixelCameraModel(calibration.getBias(), calibration.getCountPerPhoton());
				break;
			case SCMOS:
				cameraModel = CameraModelManager.load(calibration.getCameraModelName());
				if (cameraModel == null)
				{
					throw new IllegalStateException("No camera model for camera type: " + calibration.getCameraType());
				}
				cameraModel = PeakFit.cropCameraModel(cameraModel, imp.getWidth(), imp.getHeight(), ox, oy, false);
				// Store for next time
				Rectangle bounds = cameraModel.getBounds();
				ox = bounds.x;
				oy = bounds.y;
				// Reset origin for filtering
				if (ox != 0 || oy != 0)
				{
					cameraModel = cameraModel.copy();
					cameraModel.setOrigin(0, 0);
				}
				break;
			case CAMERA_TYPE_NA:
			case UNRECOGNIZED:
			default:
				cameraModel = new NullCameraModel();
				yAxisTitle = Y_AXIS_COUNT;
				break;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.gui.DialogListener#dialogItemChanged(ij.gui.GenericDialog, java.awt.AWTEvent)
	 */
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
	{
		algorithm = gd.getNextChoiceIndex();
		algorithm2 = gd.getNextChoiceIndex();
		lowestPixelsRange = (int) gd.getNextNumber();
		if (gd.invalidNumber() || lowestPixelsRange < 1)
			return false;
		if (gd.isShowing())
			drawPlot();
		return true;
	}

	/**
	 * Build a plot of the noise estimate from the current frame.
	 * Limit the preview to 100 frames.
	 */
	private void drawPlot()
	{
		NoiseEstimatorMethod[] values = SettingsManager.getNoiseEstimatorMethodValues();
		NoiseEstimator.Method method1 = FitProtosHelper.convertNoiseEstimatorMethod(values[algorithm]);
		NoiseEstimator.Method method2 = FitProtosHelper.convertNoiseEstimatorMethod(values[algorithm2]);
		IJ.showStatus("Estimating noise ...");

		boolean twoMethods = method1 != method2;
		boolean preserveResiduals = method1.name().contains("Residuals") && method2.name().contains("Residuals") &&
				twoMethods;

		int current = imp.getCurrentSlice();
		int stackSize = imp.getStackSize();
		int preview = 100;
		int start = current;
		int end = current + preview;
		if (end > stackSize)
		{
			int shift = end - stackSize;
			start -= shift;
			end = stackSize;
			start = Math.max(1, start);
		}

		int size = end - start + 1;
		double[] xValues = new double[size];
		double[] yValues1 = new double[size];
		double[] yValues2 = (twoMethods) ? new double[size] : null;

		ImageStack stack = imp.getImageStack();
		Rectangle bounds = imp.getProcessor().getRoi();
		float[] buffer = null;
		for (int slice = start, i = 0; slice <= end; slice++, i++)
		{
			IJ.showProgress(i, size);
			final ImageProcessor ip = stack.getProcessor(slice);
			buffer = ImageConverter.getData(ip.getPixels(), ip.getWidth(), ip.getHeight(), bounds, buffer);
			cameraModel.removeBiasAndGain(bounds, buffer);
			final NoiseEstimator ne = new NoiseEstimator(buffer, bounds.width, bounds.height);
			ne.preserveResiduals = preserveResiduals;
			ne.setRange(lowestPixelsRange);
			xValues[i] = slice;
			yValues1[i] = ne.getNoise(method1);
			if (twoMethods)
				yValues2[i] = ne.getNoise(method2);
		}
		IJ.showProgress(1);

		IJ.showStatus("Plotting noise ...");

		// Get limits
		double[] a = Tools.getMinMax(xValues);
		double[] b1 = Tools.getMinMax(yValues1);
		if (twoMethods)
		{
			double[] b2 = Tools.getMinMax(yValues2);
			b1[0] = FastMath.min(b1[0], b2[0]);
			b1[1] = FastMath.max(b1[1], b2[1]);
		}

		String title = imp.getTitle() + " Noise";
		Plot2 plot = new Plot2(title, "Slice", yAxisTitle);
		double range = b1[1] - b1[0];
		if (range == 0)
			range = 1;
		plot.setLimits(a[0], a[1], b1[0] - 0.05 * range, b1[1] + 0.05 * range);
		plot.setColor(Color.blue);
		plot.addPoints(xValues, yValues1, Plot2.LINE);
		//plot.draw();
		String label = String.format("%s (Blue) = %s", trim(method1.getName()),
				Utils.rounded(new Statistics(yValues1).getMean()));
		if (twoMethods)
		{
			plot.setColor(Color.red);
			plot.addPoints(xValues, yValues2, Plot2.LINE);
			label += String.format(", %s (Red) = %s", trim(method2.getName()),
					Utils.rounded(new Statistics(yValues2).getMean()));
		}
		plot.setColor(Color.black);
		plot.addLabel(0, 0, label);

		Utils.display(title, plot);

		IJ.showStatus("");
	}

	private Object trim(String name)
	{
		return (name.length() > 20) ? name.substring(0, 20) + "..." : name;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
	 */
	public void run(ImageProcessor ip)
	{
		// Perform all methods and add to the results
		double[] result = new double[NoiseEstimator.Method.values().length + 1];
		int i = 0;
		result[i++] = (pfr == null) ? 1 : pfr.getSliceNumber();
		Rectangle bounds = ip.getRoi();
		float[] buffer = ImageConverter.getData(ip.getPixels(), ip.getWidth(), ip.getHeight(), bounds, null);
		cameraModel.removeBiasAndGain(bounds, buffer);
		NoiseEstimator ne = new NoiseEstimator(buffer, bounds.width, bounds.height);
		ne.preserveResiduals = true;
		for (NoiseEstimator.Method m : NoiseEstimator.Method.values())
		{
			result[i++] = ne.getNoise(m);
		}
		results.add(result);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.ExtendedPlugInFilter#setNPasses(int)
	 */
	public void setNPasses(int nPasses)
	{
		// Do nothing
	}

	private void showResults()
	{
		Collections.sort(results, new Comparator<double[]>()
		{
			public int compare(double[] o1, double[] o2)
			{
				// Sort on slice number
				return (o1[0] < o2[0]) ? -1 : 1;
			}
		});

		// Slow when there are lots of results ... Could change the output options in the future 
		TextWindow tw = new TextWindow(imp.getTitle() + " Noise", createHeader(), "", 800, 400);
		for (double[] result : results)
		{
			tw.append(createResult(result));
		}
	}

	private String createHeader()
	{
		StringBuilder sb = new StringBuilder("Slice");
		for (NoiseEstimator.Method m : NoiseEstimator.Method.values())
		{
			sb.append('\t').append(m);
		}
		return sb.toString();
	}

	private String createResult(double[] result)
	{
		StringBuilder sb = new StringBuilder("" + (int) result[0]);
		for (int i = 1; i < result.length; i++)
		{
			sb.append('\t').append(result[i]);
		}
		return sb.toString();
	}
}
