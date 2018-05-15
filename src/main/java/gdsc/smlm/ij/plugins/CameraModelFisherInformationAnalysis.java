package gdsc.smlm.ij.plugins;

import java.awt.Color;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import gdsc.core.ij.Utils;
import gdsc.core.utils.Maths;
import gdsc.core.utils.TextUtils;
import gdsc.core.utils.TurboList;
import gdsc.smlm.data.config.GUIProtos.CameraModelFisherInformationAnalysisSettings;
import gdsc.smlm.function.BasePoissonFisherInformation;
import gdsc.smlm.function.DiscretePoissonGaussianFisherInformation;
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

		gd.addSlider("Min_exponent", -10, 4, settings.getMinExponent());
		gd.addSlider("Max_exponent", -10, 4, settings.getMaxExponent());
		gd.addSlider("Sub_divisions", 0, 10, settings.getSubDivisions());
		gd.addNumericField("CCD_gain", settings.getCcdGain(), 2);
		gd.addNumericField("CCD_noise", settings.getCcdNoise(), 2);
		gd.addNumericField("EM-CCD_gain", settings.getEmCcdGain(), 2);
		gd.addNumericField("EM-CCD_noise", settings.getEmCcdNoise(), 2);

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
			photons[i] = Math.pow(10, exp[i]);

		double[] pgFI = getFisherInformation(photons, pg, true);
		double[] pgaFI = getFisherInformation(photons, pga, false);
		double[] pggFI = getFisherInformation(photons, pgg, true);

		// Compute relative to the Poisson Fisher information
		double[] rpgFI = getAlpha(pgFI, photons);
		double[] rpgaFI = getAlpha(pgaFI, photons);
		double[] rpggFI = getAlpha(pggFI, photons);

		Color color1 = Color.BLUE;
		Color color2 = Color.GREEN;
		Color color3 = Color.RED;

		WindowOrganiser wo = new WindowOrganiser();

		String title = "Relative Fisher Information";
		Plot plot = new Plot(title, "photons", "Noise coefficient (alpha)");
		plot.setLimits(photons[0], photons[photons.length - 1], 0, 1);
		plot.setLogScaleX();
		plot.setColor(color1);
		plot.addPoints(photons, rpgFI, Plot.LINE);
		plot.setColor(color2);
		plot.addPoints(photons, rpgaFI, Plot.LINE);
		plot.setColor(color3);
		plot.addPoints(photons, rpggFI, Plot.LINE);
		plot.setColor(Color.BLACK);
		plot.addLegend("CCD\nCCD approx\nEM CCD");
		Utils.display(title, plot, 0, wo);

		title = "Fisher Information";
		plot = new Plot(title, "photons", "Fisher Information");
		double[] limits = Maths.limits(pgFI);
		limits = Maths.limits(limits, pgaFI);
		limits = Maths.limits(limits, pggFI);
		plot.setLimits(photons[0], photons[photons.length - 1], limits[0], limits[1]);
		plot.setLogScaleX();
		plot.setLogScaleY();
		plot.setColor(color1);
		plot.addPoints(photons, pgFI, Plot.LINE);
		plot.setColor(color2);
		plot.addPoints(photons, pgaFI, Plot.LINE);
		plot.setColor(color3);
		plot.addPoints(photons, pggFI, Plot.LINE);
		plot.setColor(Color.BLACK);
		plot.addLegend("CCD\nCCD approx\nEM CCD");
		Utils.display(title, plot, 0, wo);

		wo.tile();
	}

	private PoissonGaussianFisherInformation createPoissonGaussianFisherInformation(boolean discrete)
	{
		double s = settings.getCcdNoise() / settings.getCcdGain();
		if (s <= 0)
		{
			IJ.error(TITLE, "CCD noise must be positive");
			return null;
		}
		if (discrete)
			return new DiscretePoissonGaussianFisherInformation(s);
		return new RealPoissonGaussianFisherInformation(s);
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
		return new RealPoissonGammaGaussianFisherInformation(m, s);
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
			int nThreads = Prefs.getThreads();
			if (es == null)
				es = Executors.newFixedThreadPool(nThreads);
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
							f[ii] = fi2.getFisherInformation(photons[ii]);
					}
				}));
			}
			Utils.waitForCompletion(futures);
		}
		else
		{
			// Simple single threaded method.
			for (int i = 0; i < f.length; i++)
			{
				f[i] = fi.getFisherInformation(photons[i]);

				//if (fi instanceof PoissonGammaGaussianFisherInformation)
				//{
				//	PoissonGammaGaussianFisherInformation pgg = (PoissonGammaGaussianFisherInformation) fi;
				//	double[][] data = pgg.getFisherInformationFunction(false);
				//
				//	double[] fif = data[1];
				//	int max = 0;
				//	for (int j = 1; j < fif.length; j++)
				//		if (fif[max] < fif[j])
				//			max = j;
				//	System.out.printf("PGG(p=%g) max=%g\n", photons[i], data[0][max]);
				//
				//	// Debugging
				//	if (photons[i] > 500.999 && photons[i] < 100.001)
				//	{
				//		//PoissonGammaGaussianFisherInformation pgg = (PoissonGammaGaussianFisherInformation)fi;
				//		//double[][] data = pgg.getFisherInformationFunction(false);
				//		String title = TITLE + " " + photons[i];
				//		Plot plot = new Plot(title, "Count", "FI function", data[0], data[1]);
				//		Utils.display(title, plot);
				//	}
				//}
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
}
