package gdsc.smlm.ij.plugins;

import java.awt.Color;

import gdsc.core.ij.Utils;
import gdsc.core.utils.Maths;
import gdsc.core.utils.TextUtils;
import gdsc.smlm.data.config.GUIProtos.CameraModelFisherInformationAnalysisSettings;
import gdsc.smlm.function.FisherInformation;
import gdsc.smlm.function.PoissonGaussianApproximationFisherInformation;
import gdsc.smlm.function.PoissonGaussianFisherInformation;
import gdsc.smlm.function.RealPoissonGaussianFisherInformation;
import gdsc.smlm.ij.settings.SettingsManager;
import gnu.trove.list.array.TDoubleArrayList;
import ij.IJ;
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
	// Let the function store the last computed convolution.
	// Create a plugin to compute the alpha parameter across a range of means.
	// Draw the P-G convolution and the computed function.
	// Is it smooth enough for a Simpson sum?
	// Maybe the simpson sum is not working at low mean as the convolution/function
	// is not smooth.

	private static final String TITLE = "Camera Model Fisher Information Analysis";

	private CameraModelFisherInformationAnalysisSettings.Builder settings;

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

		double[] exp = createExponents();
		if (exp == null)
			return;

		double[] photons = new double[exp.length];
		for (int i = 0; i < photons.length; i++)
			photons[i] = Math.pow(10, exp[i]);

		double[] pgFI = getFisherInformation(photons, pg);
		double[] pgaFI = getFisherInformation(photons, pga);

		// Compute relative to the Poisson Fisher information
		double[] rpgFI = getAlpha(pgFI, photons);
		double[] rpgaFI = getAlpha(pgaFI, photons);

		Color color1 = Color.BLUE;
		Color color2 = Color.GREEN;
		//Color color3 = Color.RED;

		WindowOrganiser wo = new WindowOrganiser();

		String title = "Relative Fisher Information";
		Plot plot = new Plot(title, "photons", "Noise coefficient (alpha)");
		plot.setLimits(photons[0], photons[photons.length - 1], 0, 1);
		plot.setLogScaleX();
		plot.setColor(color1);
		plot.addPoints(photons, rpgFI, Plot.LINE);
		plot.setColor(color2);
		plot.addPoints(photons, rpgaFI, Plot.LINE);
		plot.setColor(Color.BLACK);
		plot.addLegend("CCD\nCCD approx");
		Utils.display(title, plot, 0, wo);

		title = "Fisher Information";
		plot = new Plot(title, "photons", "Fisher Information");
		double[] limits = Maths.limits(pgFI);
		limits = Maths.limits(limits, pgaFI);
		plot.setLimits(photons[0], photons[photons.length - 1], limits[0], limits[1]);
		plot.setLogScaleX();
		plot.setLogScaleY();
		plot.setColor(color1);
		plot.addPoints(photons, pgFI, Plot.LINE);
		plot.setColor(color2);
		plot.addPoints(photons, pgaFI, Plot.LINE);
		plot.setColor(Color.BLACK);
		plot.addLegend("CCD\nCCD approx");
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

	private double[] getFisherInformation(double[] photons, FisherInformation fi)
	{
		// Simple single threaded method.
		// TODO - multithread for speed.
		double[] f = new double[photons.length];
		for (int i = 0; i < f.length; i++)
			f[i] = fi.getFisherInformation(photons[i]);
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
