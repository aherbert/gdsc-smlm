package gdsc.smlm.ij.plugins;

import java.awt.Color;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.math3.util.FastMath;

import gdsc.core.ij.IJTrackProgress;
import gdsc.core.ij.Utils;
import gdsc.core.logging.Ticker;
import gdsc.core.utils.Maths;
import gdsc.core.utils.TextUtils;
import gdsc.core.utils.TurboList;
import gdsc.smlm.data.config.GUIProtos.CameraModelFisherInformationAnalysisSettings;
import gdsc.smlm.function.BasePoissonFisherInformation;
import gdsc.smlm.function.HalfPoissonFisherInformation;
import gdsc.smlm.function.InterpolatedPoissonFisherInformation;
import gdsc.smlm.function.PoissonGammaGaussianFisherInformation;
import gdsc.smlm.function.PoissonGaussianApproximationFisherInformation;
import gdsc.smlm.function.PoissonGaussianFisherInformation;
import gdsc.smlm.function.RealPoissonGammaGaussianFisherInformation;
import gdsc.smlm.function.RealPoissonGaussianFisherInformation;
import gdsc.smlm.ij.settings.SettingsManager;
import gnu.trove.list.array.TDoubleArrayList;
import ij.IJ;
import ij.Prefs;
import ij.gui.ExtendedGenericDialog;
import ij.gui.NonBlockingExtendedGenericDialog;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.plugin.WindowOrganiser;

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

/**
 * Model the Fisher information from an EM-CCD camera, CCD or sCMOS camera.
 */
public class CameraModelFisherInformationAnalysis implements PlugIn
{
	// TODO 
	// Options to show the computed convolution across a range of means.

	private static final String TITLE = "Camera Model Fisher Information Analysis";

	private static final String[] POINT_OPTION = { "None", "X", "Circle", "Box", "Cross" };

	private CameraModelFisherInformationAnalysisSettings.Builder settings;

	private ExecutorService es = null;

	//private boolean extraOptions;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);
		//extraOptions = Utils.isExtraOptions();

		if (!showDialog())
			return;

		analyse();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.ExtendedPlugInFilter#showDialog(ij.ImagePlus, java.lang.String,
	 * ij.plugin.filter.PlugInFilterRunner)
	 */
	private boolean showDialog()
	{
		settings = SettingsManager.readCameraModelFisherInformationAnalysisSettings(0).toBuilder();

		NonBlockingExtendedGenericDialog gd = new NonBlockingExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		//@formatter:off
		gd.addMessage(TextUtils.wrap(
				"Compute Fisher information for a CCD/EM-CCD camera model. " +
				"Configure the range of photons using a log10 scale.", 80));
		//@formatter:on

		gd.addSlider("Min_exponent", -50, 4, settings.getMinExponent());
		gd.addSlider("Max_exponent", -10, 4, settings.getMaxExponent());
		gd.addSlider("Sub_divisions", 0, 10, settings.getSubDivisions());
		gd.addNumericField("CCD_gain", settings.getCcdGain(), 2);
		gd.addNumericField("CCD_noise", settings.getCcdNoise(), 2);
		gd.addNumericField("EM-CCD_gain", settings.getEmCcdGain(), 2);
		gd.addNumericField("EM-CCD_noise", settings.getEmCcdNoise(), 2);
		gd.addChoice("Plot_point", POINT_OPTION, settings.getPointOption());

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		settings.setMinExponent((int) gd.getNextNumber());
		settings.setMaxExponent((int) gd.getNextNumber());
		settings.setSubDivisions((int) gd.getNextNumber());
		settings.setCcdGain(gd.getNextNumber());
		settings.setCcdNoise(gd.getNextNumber());
		settings.setEmCcdGain(gd.getNextNumber());
		settings.setEmCcdNoise(gd.getNextNumber());
		settings.setPointOption(gd.getNextChoiceIndex());

		SettingsManager.writeSettings(settings);

		if (settings.getMinExponent() > settings.getMaxExponent())
		{
			IJ.error(TITLE, "Min exponent must be less or equal to max exponent");
			return false;
		}

		return true;
	}

	private void analyse()
	{
		PoissonGaussianFisherInformation pg = createPoissonGaussianFisherInformation(false);
		if (pg == null)
			return;
		PoissonGaussianApproximationFisherInformation pga = createPoissonGaussianApproximationFisherInformation();

		PoissonGammaGaussianFisherInformation pgg = createPoissonGammaGaussianFisherInformation();
		if (pgg == null)
			return;

		double[] exp = createExponents();
		if (exp == null)
			return;

		double[] photons = new double[exp.length];
		for (int i = 0; i < photons.length; i++)
			photons[i] = FastMath.pow(10, exp[i]);

		double[] pgFI = getFisherInformation(photons, pg, true);
		double[] pgaFI = getFisherInformation(photons, pga, false);
		double[] pggFI = getFisherInformation(photons, pgg, true);

		// ==============
		{
			// Debug PGG
			double t = 10;
			pgg.getFisherInformation(t);
			double[][] data1 = pgg.getFisherInformationFunction(false);
			double[][] data2 = pgg.getFisherInformationFunction(true);

			double[] fif = data1[1];
			int max = 0;
			for (int j = 1; j < fif.length; j++)
				if (fif[max] < fif[j])
					max = j;
			System.out.printf("PGG(p=%g) max=%g\n", t, data1[0][max]);

			String title = TITLE + " " + t;
			Plot plot = new Plot(title, "Count", "FI function");
			double yMax = Maths.max(data1[1]);
			yMax = Maths.maxDefault(yMax, data2[1]);
			plot.setLimits(data2[0][0], data2[0][data2[0].length-1], 0, yMax);
			plot.setColor(Color.red);
			plot.addPoints(data1[0], data1[1], Plot.LINE);
			plot.setColor(Color.blue);
			plot.addPoints(data2[0], data2[1], Plot.LINE);
			Utils.display(title, plot);
		}

		// ==============

		// Compute relative to the Poisson Fisher information
		double[] rpgFI = getAlpha(pgFI, photons);
		//double[] rpgaFI = getAlpha(pgaFI, photons);
		double[] rpggFI = getAlpha(pggFI, photons);

		Color color1 = Color.BLUE;
		Color color2 = Color.GREEN;
		Color color3 = Color.RED;

		WindowOrganiser wo = new WindowOrganiser();

		// Test if we can use ImageJ support for a X log scale
		boolean logScaleX = ((float) photons[0] != 0);
		double[] x = (logScaleX) ? photons : exp;
		String xTitle = (logScaleX) ? "photons" : "log10(photons)";

		// Get interpolation for alpha
		double[] logU = exp.clone();
		double scale = Math.log(10);
		for (int i = 0; i < logU.length; i++)
			logU[i] *= scale;
		BasePoissonFisherInformation ipg = new InterpolatedPoissonFisherInformation(logU, rpgFI, true, pga);
		BasePoissonFisherInformation ipga = pga; // No need for interpolation
		// TODO - change this if the PGG FI is not fixed alpha at low mean.
		BasePoissonFisherInformation ipgg = new InterpolatedPoissonFisherInformation(logU, rpggFI, false,
				new HalfPoissonFisherInformation());

		// Interpolate with 300 points for smooth curve
		int n = 300;
		double[] iexp = new double[n + 1];
		double[] iphotons = new double[iexp.length];
		double h = (exp[exp.length - 1] - exp[0]) / n;
		for (int i = 0; i <= n; i++)
		{
			iexp[i] = exp[0] + i * h;
			iphotons[i] = FastMath.pow(10, iexp[i]);
		}
		double[] ix = (logScaleX) ? iphotons : iexp;
		double[] irpgFI = getAlpha(ipg, iphotons);
		double[] irpgaFI = getAlpha(ipga, iphotons);
		double[] irpggFI = getAlpha(ipgg, iphotons);

		int pointShape = getPointShape();

		String title = "Relative Fisher Information";
		Plot plot = new Plot(title, xTitle, "Noise coefficient (alpha)");
		plot.setLimits(x[0], x[x.length - 1], 0, 1);
		if (logScaleX)
			plot.setLogScaleX();
		plot.setColor(color1);
		plot.addPoints(ix, irpgFI, Plot.LINE);
		plot.setColor(color2);
		plot.addPoints(ix, irpgaFI, Plot.LINE);
		plot.setColor(color3);
		plot.addPoints(ix, irpggFI, Plot.LINE);
		plot.setColor(Color.BLACK);
		plot.addLegend("CCD\nCCD approx\nEM CCD");
		// Option to show nodes
		if (pointShape != -1)
		{
			plot.setColor(color1);
			plot.addPoints(x, rpgFI, pointShape);
			plot.setColor(color3);
			plot.addPoints(x, rpggFI, pointShape);
			plot.setColor(Color.BLACK);
		}
		Utils.display(title, plot, 0, wo);

		// The approximation should not produce an infinite computation
		double[] limits = new double[] { pgaFI[pgaFI.length - 1], pgaFI[pgaFI.length - 1] };
		limits = limits(limits, pgFI);
		limits = limits(limits, pgaFI);
		limits = limits(limits, pggFI);

		// Check if we can use ImageJ support for a Y log scale
		boolean logScaleY = ((float) limits[1] <= Float.MAX_VALUE);
		if (!logScaleY)
		{
			for (int i = 0; i < pgFI.length; i++)
			{
				pgFI[i] = Math.log10(pgFI[i]);
				pgaFI[i] = Math.log10(pgaFI[i]);
				pggFI[i] = Math.log10(pggFI[i]);
			}
			limits[0] = Math.log10(limits[0]);
			limits[1] = Math.log10(limits[1]);
		}

		String yTitle = (logScaleY) ? "Fisher Information" : "log10(Fisher Information)";

		title = "Fisher Information";
		plot = new Plot(title, xTitle, yTitle);

		plot.setLimits(x[0], x[x.length - 1], limits[0], limits[1]);
		if (logScaleX)
			plot.setLogScaleX();
		if (logScaleY)
			plot.setLogScaleY();
		plot.setColor(color1);
		plot.addPoints(x, pgFI, Plot.LINE);
		plot.setColor(color2);
		plot.addPoints(x, pgaFI, Plot.LINE);
		plot.setColor(color3);
		plot.addPoints(x, pggFI, Plot.LINE);
		plot.setColor(Color.BLACK);
		plot.addLegend("CCD\nCCD approx\nEM CCD");
		//// Option to show nodes
		// This gets messy as the lines are straight
		//if (pointShape != -1)
		//{
		//	plot.setColor(color1);
		//	plot.addPoints(x, pgFI, pointShape);
		//	plot.setColor(color3);
		//	plot.addPoints(x, pggFI, pointShape);
		//	plot.setColor(Color.BLACK);
		//}
		Utils.display(title, plot, 0, wo);

		wo.tile();
	}

	private double[] limits(double[] limits, double[] f)
	{
		double min = limits[0];
		double max = limits[1];
		for (int i = 0; i < f.length; i++)
		{
			double d = f[i];
			// Find limits of numbers that can be logged
			if (d <= 0 || d == Double.POSITIVE_INFINITY)
				continue;
			if (min > d)
				min = d;
			else if (max < d)
				max = d;
		}
		limits[0] = min;
		limits[1] = max;
		return limits;
	}

	private PoissonGaussianFisherInformation createPoissonGaussianFisherInformation(boolean discrete)
	{
		double s = settings.getCcdNoise() / settings.getCcdGain();
		if (s <= 0)
		{
			IJ.error(TITLE, "CCD noise must be positive");
			return null;
		}

		int sampling = PoissonGaussianFisherInformation.DEFAULT_SAMPLING;
		//sampling <<= 2;
		PoissonGaussianFisherInformation fi = new RealPoissonGaussianFisherInformation(s, sampling);
		//fi.setCumulativeProbability(1 - 1e-12);
		fi.setMinRange(5);
		return fi;
	}

	private PoissonGaussianApproximationFisherInformation createPoissonGaussianApproximationFisherInformation()
	{
		double s = settings.getCcdNoise() / settings.getCcdGain();
		if (s <= 0)
		{
			IJ.error(TITLE, "CCD noise must be positive");
			return null;
		}
		return new PoissonGaussianApproximationFisherInformation(s);
	}

	private PoissonGammaGaussianFisherInformation createPoissonGammaGaussianFisherInformation()
	{
		double s = settings.getEmCcdNoise();
		double m = settings.getEmCcdGain();
		if (s <= 0)
		{
			IJ.error(TITLE, "EM CCD noise must be positive");
			return null;
		}
		if (m <= 0)
		{
			IJ.error(TITLE, "EM CCD gain must be positive");
			return null;
		}

		int sampling = PoissonGammaGaussianFisherInformation.DEFAULT_SAMPLING;
		sampling <<= 2;
		PoissonGammaGaussianFisherInformation fi = new RealPoissonGammaGaussianFisherInformation(m, s, sampling);
		//fi.setCumulativeProbability(1 - 1e-8);
		fi.setMinRange(10);
		//fi.setMaxRange(10);
		fi.setMeanThreshold(200);
		return fi;
	}

	private double[] createExponents()
	{
		int n = 1 + Math.max(0, settings.getSubDivisions());
		double h = 1.0 / n;
		double minExp = settings.getMinExponent();
		double maxExp = settings.getMaxExponent();
		double size = (maxExp - minExp) * n + 1;
		if (size > 100)
		{
			ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
			gd.addMessage("Number of exponents is " + Math.ceil(size) + ". OK to continue?");
			gd.showDialog();
			if (gd.wasCanceled())
				return null;
		}
		TDoubleArrayList list = new TDoubleArrayList();
		for (int i = 0;; i++)
		{
			double e = minExp + i * h;
			list.add(e);
			if (e >= settings.getMaxExponent())
				break;
		}
		return list.toArray();
	}

	private double[] getFisherInformation(final double[] photons, final BasePoissonFisherInformation fi,
			boolean multithread)
	{
		final double[] f = new double[photons.length];
		if (multithread)
		{
			Utils.showStatus("Computing " + fi.getClass().getSimpleName());
			int nThreads = Prefs.getThreads();
			if (es == null)
				es = Executors.newFixedThreadPool(nThreads);
			final Ticker ticker = Ticker.createStarted(new IJTrackProgress(), f.length, nThreads != 1);
			int nPerThread = (int) Math.ceil((double) f.length / nThreads);
			TurboList<Future<?>> futures = new TurboList<Future<?>>(nThreads);
			for (int i = 0; i < f.length; i += nPerThread)
			{
				final int start = i;
				final int end = Math.min(f.length, i + nPerThread);
				futures.add(es.submit(new Runnable()
				{
					public void run()
					{
						BasePoissonFisherInformation fi2 = fi.clone();
						for (int ii = start; ii < end; ii++)
						{
							f[ii] = fi2.getFisherInformation(photons[ii]);
							ticker.tick();
						}
					}
				}));
			}
			Utils.waitForCompletion(futures);
			ticker.stop();
			Utils.showStatus("");
		}
		else
		{
			// Simple single threaded method.
			for (int i = 0; i < f.length; i++)
			{
				f[i] = fi.getFisherInformation(photons[i]);
			}
		}
		return f;
	}

	private double[] getAlpha(double[] I, double[] photons)
	{
		double[] rI = new double[photons.length];
		for (int i = 0; i < photons.length; i++)
		{
			rI[i] = I[i] * photons[i];
		}
		return rI;
	}

	private double[] getAlpha(BasePoissonFisherInformation fi, double[] photons)
	{
		// Do not multi-thread as this is an interpolated/fast function
		double[] rI = new double[photons.length];
		for (int i = 0; i < photons.length; i++)
		{
			rI[i] = fi.getAlpha(photons[i]);
		}
		return rI;
	}

	private int getPointShape()
	{
		switch (settings.getPointOption())
		{
			case 1:
				return Plot.X;
			case 2:
				return Plot.CIRCLE;
			case 3:
				return Plot.BOX;
			case 4:
				return Plot.CROSS;
		}
		return -1;
	}
}
