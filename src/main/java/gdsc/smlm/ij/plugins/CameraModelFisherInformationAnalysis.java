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
import gdsc.smlm.data.NamedObject;
import gdsc.smlm.data.config.GUIProtos.CameraModelFisherInformationAnalysisSettings;
import gdsc.smlm.function.BasePoissonFisherInformation;
import gdsc.smlm.function.HalfPoissonFisherInformation;
import gdsc.smlm.function.InterpolatedPoissonFisherInformation;
import gdsc.smlm.function.PoissonFisherInformation;
import gdsc.smlm.function.PoissonGammaGaussianFisherInformation;
import gdsc.smlm.function.PoissonGaussianApproximationFisherInformation;
import gdsc.smlm.function.PoissonGaussianFisherInformation;
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

	//@formatter:off
	private enum CameraType implements NamedObject
	{
		POISSON { public String getName() { return "Poisson"; } },
		CCD { public String getName() { return "CCD"; } 
			  public boolean isFast() { return false; }
			  public boolean isLowerFixedI() { return true; } },
		CCD_APPROXIMATION { public String getName() { return "CCD Approximation"; }
							public String getShortName() { return "CCD Approx"; }},
		EM_CCD { public String getName() { return "EM-CCD"; } 
		         public boolean isFast() { return false; } };

		public abstract String getName();
		
		/**
		 * Checks if is fast to compute.
		 *
		 * @return true, if is fast
		 */
		public boolean isFast() { return false; }
		
		/**
		 * Checks if is fixed Fisher information at the lower bound.
		 *
		 * @return true, if is lower fixed I
		 */
		public boolean isLowerFixedI() { return false; }

		public String getShortName()
		{
			return getName();
		}
		
		public static CameraType forNumber(int number)
		{
			CameraType[] values = CameraType.values();
			if (number >= 0 && number < values.length)
				return values[number];
			return values[0];
		}
	}
	//@formatter:on

	private static String[] CAMERA_TYPES = SettingsManager.getNames((Object[]) CameraType.values());

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
		gd.addChoice("Camera_1_type", CAMERA_TYPES, settings.getCamera1Type());
		gd.addNumericField("Camera_1_gain", settings.getCamera1Gain(), 2);
		gd.addNumericField("Camera_1_noise", settings.getCamera1Noise(), 2);
		gd.addChoice("Camera_2_type", CAMERA_TYPES, settings.getCamera2Type());
		gd.addNumericField("Camera_2_gain", settings.getCamera2Gain(), 2);
		gd.addNumericField("Camera_2_noise", settings.getCamera2Noise(), 2);
		gd.addChoice("Plot_point", POINT_OPTION, settings.getPointOption());

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		settings.setMinExponent((int) gd.getNextNumber());
		settings.setMaxExponent((int) gd.getNextNumber());
		settings.setSubDivisions((int) gd.getNextNumber());
		settings.setCamera1Type(gd.getNextChoiceIndex());
		settings.setCamera1Gain(gd.getNextNumber());
		settings.setCamera1Noise(gd.getNextNumber());
		settings.setCamera2Type(gd.getNextChoiceIndex());
		settings.setCamera2Gain(gd.getNextNumber());
		settings.setCamera2Noise(gd.getNextNumber());
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
		CameraType type1 = CameraType.forNumber(settings.getCamera1Type());
		CameraType type2 = CameraType.forNumber(settings.getCamera2Type());

		BasePoissonFisherInformation f1 = createPoissonFisherInformation(type1, settings.getCamera1Gain(),
				settings.getCamera1Noise());
		if (f1 == null)
			return;
		BasePoissonFisherInformation f2 = createPoissonFisherInformation(type2, settings.getCamera2Gain(),
				settings.getCamera2Noise());
		if (f2 == null)
			return;

		double[] exp = createExponents();
		if (exp == null)
			return;

		double[] photons = new double[exp.length];
		for (int i = 0; i < photons.length; i++)
			photons[i] = FastMath.pow(10, exp[i]);

		double[] fi1 = getFisherInformation(photons, f1, type1);
		double[] fi2 = getFisherInformation(photons, f2, type2);

		// Compute relative to the Poisson Fisher information
		double[] alpha1 = getAlpha(fi1, photons);
		double[] alpha2 = getAlpha(fi2, photons);

		// ==============
		// Debug PGG
		boolean debug = false;
		// ==============
		if (debug && f2 instanceof PoissonGammaGaussianFisherInformation)
		{
			PoissonGammaGaussianFisherInformation pgg = (PoissonGammaGaussianFisherInformation) f2;
			double t = 200;
			//for (int i=0; i<rpggFI.length; i++)
			//{
			//	if (photons[i] > 100 && rpggFI[i] < 0.3)
			//		t = photons[i];
			//}

			double fi = pgg.getFisherInformation(t);
			double alpha = fi * t;
			double[][] data1 = pgg.getFisherInformationFunction(false);
			double[][] data2 = pgg.getFisherInformationFunction(true);

			double[] fif = data1[1];
			int max = 0;
			for (int j = 1; j < fif.length; j++)
				if (fif[max] < fif[j])
					max = j;
			System.out.printf("PGG(p=%g) max=%g\n", t, data1[0][max]);

			String title = TITLE + " photons=" + Utils.rounded(t) + " alpha=" + Utils.rounded(alpha);
			Plot plot = new Plot(title, "Count", "FI function");
			double yMax = Maths.max(data1[1]);
			yMax = Maths.maxDefault(yMax, data2[1]);
			plot.setLimits(data2[0][0], data2[0][data2[0].length - 1], 0, yMax);
			plot.setColor(Color.red);
			plot.addPoints(data1[0], data1[1], Plot.LINE);
			plot.setColor(Color.blue);
			plot.addPoints(data2[0], data2[1], Plot.LINE);
			Utils.display(title, plot);
		}
		// ==============

		Color color1 = Color.BLUE;
		Color color2 = Color.RED;

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
		BasePoissonFisherInformation if1 = getInterpolatedPoissonFisherInformation(type1, logU, alpha1, f1);
		BasePoissonFisherInformation if2 = getInterpolatedPoissonFisherInformation(type2, logU, alpha2, f2);

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
		double[] ialpha1 = getAlpha(if1, iphotons);
		double[] ialpha2 = getAlpha(if2, iphotons);

		int pointShape = getPointShape();

		String name1 = getName(type1, settings.getCamera1Gain(), settings.getCamera1Noise());
		String name2 = getName(type2, settings.getCamera2Gain(), settings.getCamera2Noise());
		String legend = name1 + "\n" + name2;

		String title = "Relative Fisher Information";
		Plot plot = new Plot(title, xTitle, "Noise coefficient (alpha)");
		plot.setLimits(x[0], x[x.length - 1], -0.05, 1.05);
		if (logScaleX)
			plot.setLogScaleX();
		plot.setColor(color1);
		plot.addPoints(ix, ialpha1, Plot.LINE);
		plot.setColor(color2);
		plot.addPoints(ix, ialpha2, Plot.LINE);
		plot.setColor(Color.BLACK);
		plot.addLegend(legend);
		// Option to show nodes
		if (pointShape != -1)
		{
			plot.setColor(color1);
			plot.addPoints(x, alpha1, pointShape);
			plot.setColor(color2);
			plot.addPoints(x, alpha2, pointShape);
			plot.setColor(Color.BLACK);
		}
		Utils.display(title, plot, 0, wo);

		// The approximation should not produce an infinite computation
		double[] limits = new double[] { fi2[fi2.length - 1], fi2[fi2.length - 1] };
		limits = limits(limits, fi1);
		limits = limits(limits, fi2);

		// Check if we can use ImageJ support for a Y log scale
		boolean logScaleY = ((float) limits[1] <= Float.MAX_VALUE);
		if (!logScaleY)
		{
			for (int i = 0; i < fi1.length; i++)
			{
				fi1[i] = Math.log10(fi1[i]);
				fi2[i] = Math.log10(fi2[i]);
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
		plot.addPoints(x, fi1, Plot.LINE);
		plot.setColor(color2);
		plot.addPoints(x, fi2, Plot.LINE);
		plot.setColor(Color.BLACK);
		plot.addLegend(legend);
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

	private BasePoissonFisherInformation createPoissonFisherInformation(CameraType type, double gain, double noise)
	{
		switch (type)
		{
			case CCD:
				return createPoissonGaussianFisherInformation(noise / gain);
			case CCD_APPROXIMATION:
				return createPoissonGaussianApproximationFisherInformation(noise / gain);
			case EM_CCD:
				return createPoissonGammaGaussianFisherInformation(gain, noise);
			case POISSON:
				return new PoissonFisherInformation();
			default:
				throw new IllegalStateException("Unknown camera type: " + type);
		}
	}

	private String getName(CameraType type, double gain, double noise)
	{
		String name = type.getShortName();
		switch (type)
		{
			case CCD:
			case CCD_APPROXIMATION:
			case EM_CCD:
				name += String.format(" g=%.1f,n=%.1f", gain, noise);
			default:
				break;
		}
		return name;
	}

	private PoissonGaussianFisherInformation createPoissonGaussianFisherInformation(double s)
	{
		if (s < 0)
		{
			IJ.error(TITLE, "CCD noise must be positive");
			return null;
		}

		int sampling = PoissonGaussianFisherInformation.DEFAULT_SAMPLING;
		//sampling <<= 2;
		PoissonGaussianFisherInformation fi = new PoissonGaussianFisherInformation(s, sampling);
		//fi.setCumulativeProbability(1 - 1e-12);
		//fi.setMinRange(5);
		fi.setMeanThreshold(1000);
		return fi;
	}

	private PoissonGaussianApproximationFisherInformation createPoissonGaussianApproximationFisherInformation(double s)
	{
		if (s <= 0)
		{
			IJ.error(TITLE, "CCD noise must be positive");
			return null;
		}
		return new PoissonGaussianApproximationFisherInformation(s);
	}

	private PoissonGammaGaussianFisherInformation createPoissonGammaGaussianFisherInformation(double m, double s)
	{
		if (s < 0)
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
		//sampling = 2;
		PoissonGammaGaussianFisherInformation fi = new PoissonGammaGaussianFisherInformation(m, s, sampling);
		//fi.setRelativeProbabilityThreshold(1e-6);
		//fi.setMinRange(10);
		//fi.setMaxIterations(3);
		//fi.setMaxRange(10);
		//fi.setUse38(false);
		fi.setMeanThreshold(Double.MAX_VALUE);
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
			CameraType type)
	{
		final double[] f = new double[photons.length];
		if (!type.isFast())
		{
			IJ.showStatus("Computing " + type.getName());
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
			IJ.showStatus("");
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

	private BasePoissonFisherInformation getInterpolatedPoissonFisherInformation(CameraType type, double[] logU,
			double[] alpha1, BasePoissonFisherInformation f)
	{
		if (type.isFast())
			return f;
		// This should not matter as the interpolation is only used between the ends of the range
		BasePoissonFisherInformation upperf = f;
		if (type == CameraType.EM_CCD)
			upperf = new HalfPoissonFisherInformation();
		return new InterpolatedPoissonFisherInformation(logU, alpha1, type.isLowerFixedI(), upperf);
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
}
