/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package uk.ac.sussex.gdsc.smlm.ij.plugins;

import java.awt.Color;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.math3.util.FastMath;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import ij.IJ;
import ij.Prefs;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import uk.ac.sussex.gdsc.core.ij.IJTrackProgress;
import uk.ac.sussex.gdsc.core.ij.Utils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.NonBlockingExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.logging.Ticker;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.Maths;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.Sort;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.TurboList;
import uk.ac.sussex.gdsc.smlm.data.NamedObject;
import uk.ac.sussex.gdsc.smlm.data.config.FisherProtos.AlphaSample;
import uk.ac.sussex.gdsc.smlm.data.config.FisherProtos.PoissonFisherInformationCache;
import uk.ac.sussex.gdsc.smlm.data.config.FisherProtos.PoissonFisherInformationData;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.CameraModelFisherInformationAnalysisSettings;
import uk.ac.sussex.gdsc.smlm.function.BasePoissonFisherInformation;
import uk.ac.sussex.gdsc.smlm.function.HalfPoissonFisherInformation;
import uk.ac.sussex.gdsc.smlm.function.InterpolatedPoissonFisherInformation;
import uk.ac.sussex.gdsc.smlm.function.PoissonFisherInformation;
import uk.ac.sussex.gdsc.smlm.function.PoissonGammaGaussianFisherInformation;
import uk.ac.sussex.gdsc.smlm.function.PoissonGaussianApproximationFisherInformation;
import uk.ac.sussex.gdsc.smlm.function.PoissonGaussianFisherInformation;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;

/**
 * Model the Fisher information from an EM-CCD camera, CCD or sCMOS camera.
 */
public class CameraModelFisherInformationAnalysis implements PlugIn
{
	// TODO
	// Options to show the computed convolution across a range of means.

	private static final String TITLE = "Camera Model Fisher Information Analysis";

	private static final String[] POINT_OPTION = { "None", "X", "Circle", "Box", "Cross" };

	/**
	 * The camera type.
	 */
	//@formatter:off
	public enum CameraType implements NamedObject
	{
		/** A camera with pure Poisson shot noise */
		POISSON { @Override public String getName() { return "Poisson"; } },
		/** CCD camera has Poisson shot noise and Gaussian read noise */
		CCD { @Override	public String getName() { return "CCD"; }
			  @Override	public boolean isFast() { return false; }
			  @Override	public boolean isLowerFixedI() { return true; } },
		/**
		 * CCD approximation has Poisson shot noise and uses a Poisson distribution to simulate Gaussian noise (s).
		 * The distribution is Approx(u,s) = Poisson(u) + Poisson(s) = Poisson(u+s)
		 */
		CCD_APPROXIMATION { @Override public String getName() { return "CCD Approximation"; }
							@Override public String getShortName() { return "CCD Approx"; }},
		/**
		 * EM-CCD camera has Poisson shot noise, a Gamma distribution model for EM-gain
		 * and Gaussian read noise
		 */
		EM_CCD { @Override public String getName() { return "EM-CCD"; }
		         @Override public boolean isFast() { return false; } };

		@Override
		public abstract String getName();

		/**
		 * Checks if is fast to compute.
		 *
		 * @return true, if is fast
		 */
		public boolean isFast() { return true; }

		/**
		 * Checks if is fixed Fisher information at the lower bound.
		 *
		 * @return true, if is lower fixed I
		 */
		public boolean isLowerFixedI() { return false; }

		@Override
		public String getShortName()
		{
			return getName();
		}

		/**
		 * Get the camera type for the number.
		 *
		 * @param number
		 *            the number
		 * @return the camera type
		 */
		public static CameraType forNumber(int number)
		{
			final CameraType[] values = CameraType.values();
			if (number >= 0 && number < values.length)
				return values[number];
			return values[0];
		}
	}
	//@formatter:on
	private static CameraType[] cameraTypeValues = CameraType.values();

	private static String[] CAMERA_TYPES = SettingsManager.getNames((Object[]) cameraTypeValues);

	/**
	 * Class for hashing the Fisher information settings
	 */
	private static class FIKey
	{
		// This is not part of the hashcode.
		// It is used to ensure only the latest computations are saved to disk.
		final long age = System.currentTimeMillis();

		final int type;
		final double gain;
		final double noise;

		public FIKey(int type, double gain, double noise)
		{
			this.type = type;
			this.gain = gain;
			this.noise = noise;
		}

		public FIKey(PoissonFisherInformationData d)
		{
			this(d.getType(), d.getGain(), d.getNoise());
		}

		public FIKey(CameraType type, double gain, double noise)
		{
			this(type.ordinal(), gain, noise);
		}

		@Override
		public int hashCode()
		{
			int hash = type;
			hash = 31 * hash + Double.hashCode(gain);
			hash = 31 * hash + Double.hashCode(noise);
			return hash;
		}

		@Override
		public boolean equals(Object o)
		{
			if (this == o)
				return true;
			// null check
			if (o == null)
				return false;
			// type check and cast
			if (getClass() != o.getClass())
				return false;
			final FIKey key = (FIKey) o;
			// field comparison
			return type == key.type && gain == key.gain && noise == key.noise;
		}

		public CameraType getType()
		{
			return cameraTypeValues[type];
		}

		@Override
		public String toString()
		{
			final CameraType type = getType();
			String name = type.getShortName();
			switch (type)
			{
				case CCD:
				case CCD_APPROXIMATION:
				case EM_CCD:
					name += String.format(" g=%s,n=%s", gain, noise);
				default:
					break;
			}
			return name;
		}
	}

	private static HashMap<FIKey, PoissonFisherInformationData> cache = new HashMap<>();

	static
	{
		final PoissonFisherInformationCache cacheData = new SettingsManager.ConfigurationReader<>(
				PoissonFisherInformationCache.getDefaultInstance()).read();
		if (cacheData != null)
			for (final PoissonFisherInformationData data : cacheData.getDataList())
				cache.put(new FIKey(data), data);
	}

	/**
	 * Save the data to the cache.
	 *
	 * @param key
	 *            the key
	 * @param log10photons
	 *            the log 10 photons
	 * @param alpha
	 *            the alpha
	 */
	private static void save(FIKey key, double[] log10photons, double[] alpha)
	{
		final CameraType type = key.getType();
		if (type.isFast())
			return;
		PoissonFisherInformationData data = cache.get(key);
		if (data != null)
		{
			// This should only be called if new values have been computed.
			// so assume we must merge the lists.
			// Note: The lists must be sorted.
			final AlphaSample[] list1 = data.getAlphaSampleList().toArray(new AlphaSample[0]);
			//AlphaSample[] list2 = data.getAlphaSampleList().toArray(new AlphaSample[0]);
			final TurboList<AlphaSample> list = new TurboList<>(list1.length + log10photons.length);
			int i = 0, j = 0;
			final AlphaSample.Builder a2 = AlphaSample.newBuilder();
			while (i < list1.length && j < log10photons.length)
			{
				final AlphaSample a1 = list1[i];
				final double mean = log10photons[j];
				if (a1.getLog10Mean() == mean)
				{
					list.add(a1);
					i++;
					j++;
				}
				else if (a1.getLog10Mean() < mean)
				{
					list.add(a1);
					i++;
				}
				else //if (a1.getLog10Mean() > mean)
				{
					a2.setLog10Mean(mean);
					a2.setAlpha(alpha[j++]);
					list.add(a2.build());
				}
			}
			while (i < list1.length)
				list.add(list1[i++]);
			while (j < log10photons.length)
			{
				a2.setLog10Mean(log10photons[j]);
				a2.setAlpha(alpha[j++]);
				list.add(a2.build());
			}
			final PoissonFisherInformationData.Builder b = data.toBuilder();
			b.clearAlphaSample();
			b.addAllAlphaSample(list);
			data = b.build();

			//System.out.println(data);
		}
		else
		{
			final PoissonFisherInformationData.Builder b = PoissonFisherInformationData.newBuilder();
			final AlphaSample.Builder sample = AlphaSample.newBuilder();
			b.setType(key.type);
			b.setGain(key.gain);
			b.setNoise(key.noise);
			for (int i = 0; i < log10photons.length; i++)
			{
				sample.setLog10Mean(log10photons[i]);
				sample.setAlpha(alpha[i]);
				b.addAlphaSample(sample);
			}
			data = b.build();
		}
		cache.put(key, data);

		// Also save EM-CCD data to file
		if (type == CameraType.EM_CCD)
		{
			final int t = type.ordinal();
			// Get the EM-CCD keys
			final TurboList<FIKey> list = new TurboList<>(cache.size());
			for (final FIKey k : cache.keySet())
				if (k.type == t)
					list.add(k);
			final int MAX = 10;
			if (list.size() > MAX)
				// Sort by age
				list.sort(new Comparator<FIKey>()
				{
					@Override
					public int compare(FIKey o1, FIKey o2)
					{
						// Youngest first
						return -Long.compare(o1.age, o2.age);
					}
				});
			final PoissonFisherInformationCache.Builder cacheData = PoissonFisherInformationCache.newBuilder();
			for (int i = Math.min(list.size(), MAX); i-- > 0;)
				cacheData.addData(cache.get(list.getf(i)));
			SettingsManager.writeSettings(cacheData);
		}
	}

	private static PoissonFisherInformationData load(FIKey key)
	{
		return cache.get(key);
	}

	/**
	 * Load the Poisson Fisher information data for the camera type from the cache. The gain and noise must match within
	 * an
	 * error of 1e-3.
	 *
	 * @param type
	 *            the type
	 * @param gain
	 *            the gain
	 * @param noise
	 *            the noise
	 * @return the poisson fisher information data (or null)
	 */
	public static PoissonFisherInformationData loadData(CameraType type, double gain, double noise)
	{
		return loadData(type, gain, noise, 1e-3);
	}

	/**
	 * Load the Poisson Fisher information data for the camera type from the cache.
	 *
	 * @param type
	 *            the type
	 * @param gain
	 *            the gain
	 * @param noise
	 *            the noise
	 * @param relativeError
	 *            the relative error (used to compare the gain and noise)
	 * @return the poisson fisher information data (or null)
	 */
	public static PoissonFisherInformationData loadData(CameraType type, double gain, double noise,
			double relativeError)
	{
		final FIKey key = new FIKey(type, gain, noise);
		PoissonFisherInformationData data = load(key);
		if (data == null && relativeError > 0 && relativeError < 0.1)
		{
			// Fuzzy matching
			double error = 1;
			for (final PoissonFisherInformationData d : cache.values())
			{
				final double e1 = DoubleEquality.relativeError(gain, d.getGain());
				if (e1 < relativeError)
				{
					final double e2 = DoubleEquality.relativeError(noise, d.getNoise());
					if (e2 < relativeError)
					{
						// Combined error - Euclidean distance
						final double e = e1 * e1 + e2 * e2;
						if (error > e)
						{
							error = e;
							data = d;
						}
					}
				}
			}
		}
		return data;
	}

	/**
	 * Load the Poisson Fisher information for the camera type from the cache. The gain and noise must match within an
	 * error of 1e-3.
	 *
	 * @param type
	 *            the type
	 * @param gain
	 *            the gain
	 * @param noise
	 *            the noise
	 * @return the poisson fisher information (or null)
	 */
	public static InterpolatedPoissonFisherInformation loadFunction(CameraType type, double gain, double noise)
	{
		return loadFunction(type, gain, noise, 1e-3);
	}

	/**
	 * Load the Poisson Fisher information for the camera type from the cache.
	 *
	 * @param type
	 *            the type
	 * @param gain
	 *            the gain
	 * @param noise
	 *            the noise
	 * @param relativeError
	 *            the relative error (used to compare the gain and noise)
	 * @return the poisson fisher information (or null)
	 */
	public static InterpolatedPoissonFisherInformation loadFunction(CameraType type, double gain, double noise,
			double relativeError)
	{
		if (type == null || type.isFast())
			return null;
		final FIKey key = new FIKey(type, gain, noise);
		PoissonFisherInformationData data = load(key);
		if (data == null && relativeError > 0 && relativeError < 0.1)
		{
			// Fuzzy matching
			double error = 1;
			for (final PoissonFisherInformationData d : cache.values())
			{
				final double e1 = DoubleEquality.relativeError(gain, d.getGain());
				if (e1 < relativeError)
				{
					final double e2 = DoubleEquality.relativeError(noise, d.getNoise());
					if (e2 < relativeError)
					{
						// Combined error - Euclidean distance
						final double e = e1 * e1 + e2 * e2;
						if (error > e)
						{
							error = e;
							data = d;
						}
					}
				}
			}
		}
		if (data == null)
			return null;
		// Dump the samples. Convert to base e.
		final double scale = Math.log(10);
		final TDoubleArrayList meanList = new TDoubleArrayList(data.getAlphaSampleCount());
		final TDoubleArrayList alphalist = new TDoubleArrayList(data.getAlphaSampleCount());
		for (final AlphaSample sample : data.getAlphaSampleList())
		{
			meanList.add(sample.getLog10Mean() * scale);
			alphalist.add(sample.getAlpha());
		}
		final double[] means = meanList.toArray();
		final double[] alphas = alphalist.toArray();
		Sort.sortArrays(alphas, means, true);
		final BasePoissonFisherInformation upperf = (type == CameraType.EM_CCD) ? new HalfPoissonFisherInformation()
				: createPoissonGaussianApproximationFisherInformation(noise / gain);
		return new InterpolatedPoissonFisherInformation(means, alphas, type.isLowerFixedI(), upperf);
	}

	private CameraModelFisherInformationAnalysisSettings.Builder settings;

	private ExecutorService es = null;

	/*
	 * (non-Javadoc)
	 *
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	@Override
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (Utils.isExtraOptions() && !cache.isEmpty())
		{
			plotFromCache();
			return;
		}

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

		final NonBlockingExtendedGenericDialog gd = new NonBlockingExtendedGenericDialog(TITLE);
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

	// Set the debug flag
	private final static boolean debug = System.getProperty("gdsc.smlm.debug") != null;

	private void analyse()
	{
		final CameraType type1 = CameraType.forNumber(settings.getCamera1Type());
		final CameraType type2 = CameraType.forNumber(settings.getCamera2Type());

		final FIKey key1 = new FIKey(type1, settings.getCamera1Gain(), settings.getCamera1Noise());
		final FIKey key2 = new FIKey(type2, settings.getCamera2Gain(), settings.getCamera2Noise());

		final BasePoissonFisherInformation f1 = createPoissonFisherInformation(key1);
		if (f1 == null)
			return;
		final BasePoissonFisherInformation f2 = createPoissonFisherInformation(key2);
		if (f2 == null)
			return;

		final double[] exp = createExponents();
		if (exp == null)
			return;

		final double[] photons = new double[exp.length];
		for (int i = 0; i < photons.length; i++)
			photons[i] = FastMath.pow(10, exp[i]);

		// Get the alpha. This may be from the cache
		final double[] alpha1 = getAlpha(photons, exp, f1, key1);
		final double[] alpha2 = getAlpha(photons, exp, f2, key2);

		// Compute the Poisson Fisher information
		final double[] fi1 = getFisherInformation(alpha1, photons);
		final double[] fi2 = getFisherInformation(alpha2, photons);

		// ==============
		// Debug PGG
		// ==============
		if (debug && f2 instanceof PoissonGammaGaussianFisherInformation)
		{
			final PoissonGammaGaussianFisherInformation pgg = (PoissonGammaGaussianFisherInformation) f2;
			final double t = 200;
			//for (int i=0; i<rpggFI.length; i++)
			//{
			//	if (photons[i] > 100 && rpggFI[i] < 0.3)
			//		t = photons[i];
			//}

			final double fi = pgg.getFisherInformation(t);
			final double alpha = fi * t;
			final double[][] data1 = pgg.getFisherInformationFunction(false);
			final double[][] data2 = pgg.getFisherInformationFunction(true);

			final double[] fif = data1[1];
			int max = 0;
			for (int j = 1; j < fif.length; j++)
				if (fif[max] < fif[j])
					max = j;
			System.out.printf("PGG(p=%g) max=%g\n", t, data1[0][max]);

			final String title = TITLE + " photons=" + Utils.rounded(t) + " alpha=" + Utils.rounded(alpha);
			final Plot plot = new Plot(title, "Count", "FI function");
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

		final Color color1 = Color.BLUE;
		final Color color2 = Color.RED;

		final WindowOrganiser wo = new WindowOrganiser();

		// Test if we can use ImageJ support for a X log scale
		final boolean logScaleX = ((float) photons[0] != 0);
		final double[] x = (logScaleX) ? photons : exp;
		final String xTitle = (logScaleX) ? "photons" : "log10(photons)";

		// Get interpolation for alpha. Convert to base e.
		final double[] logU = exp.clone();
		final double scale = Math.log(10);
		for (int i = 0; i < logU.length; i++)
			logU[i] *= scale;
		final BasePoissonFisherInformation if1 = getInterpolatedPoissonFisherInformation(type1, logU, alpha1, f1);
		final BasePoissonFisherInformation if2 = getInterpolatedPoissonFisherInformation(type2, logU, alpha2, f2);

		// Interpolate with 5 points per sample for smooth curve
		final int n = 5 * exp.length;
		final double[] iexp = new double[n + 1];
		final double[] iphotons = new double[iexp.length];
		final double h = (exp[exp.length - 1] - exp[0]) / n;
		for (int i = 0; i <= n; i++)
		{
			iexp[i] = exp[0] + i * h;
			iphotons[i] = FastMath.pow(10, iexp[i]);
		}
		final double[] ix = (logScaleX) ? iphotons : iexp;
		final double[] ialpha1 = getAlpha(if1, iphotons);
		final double[] ialpha2 = getAlpha(if2, iphotons);

		final int pointShape = getPointShape(settings.getPointOption());

		final String name1 = getName(key1);
		final String name2 = getName(key2);
		final String legend = name1 + "\n" + name2;

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
		final boolean logScaleY = ((float) limits[1] <= Float.MAX_VALUE);
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

		final String yTitle = (logScaleY) ? "Fisher Information" : "log10(Fisher Information)";

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

	private static BasePoissonFisherInformation createPoissonFisherInformation(FIKey key)
	{
		switch (key.getType())
		{
			case CCD:
				return createPoissonGaussianFisherInformation(key.noise / key.gain);
			case CCD_APPROXIMATION:
				return createPoissonGaussianApproximationFisherInformation(key.noise / key.gain);
			case EM_CCD:
				return createPoissonGammaGaussianFisherInformation(key.gain, key.noise);
			case POISSON:
				return new PoissonFisherInformation();
			default:
				throw new IllegalStateException("Unknown camera type: " + key.getType());
		}
	}

	private static String getName(FIKey key)
	{
		final CameraType type = key.getType();
		String name = type.getShortName();
		switch (type)
		{
			case CCD:
			case CCD_APPROXIMATION:
			case EM_CCD:
				name += String.format(" g=%.1f,n=%.1f", key.gain, key.noise);
			default:
				break;
		}
		return name;
	}

	private static PoissonGaussianFisherInformation createPoissonGaussianFisherInformation(double s)
	{
		if (s < 0)
		{
			IJ.error(TITLE, "CCD noise must be positive");
			return null;
		}

		final int sampling = PoissonGaussianFisherInformation.DEFAULT_SAMPLING;
		//sampling <<= 2;
		final PoissonGaussianFisherInformation fi = new PoissonGaussianFisherInformation(s, sampling);
		//fi.setCumulativeProbability(1 - 1e-12);
		//fi.setMinRange(5);
		fi.setMeanThreshold(1000);
		return fi;
	}

	private static PoissonGaussianApproximationFisherInformation createPoissonGaussianApproximationFisherInformation(
			double s)
	{
		if (s <= 0)
		{
			IJ.error(TITLE, "CCD noise must be positive");
			return null;
		}
		return new PoissonGaussianApproximationFisherInformation(s);
	}

	private static PoissonGammaGaussianFisherInformation createPoissonGammaGaussianFisherInformation(double m, double s)
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

		final int sampling = PoissonGammaGaussianFisherInformation.DEFAULT_SAMPLING;
		//sampling = 2;
		final PoissonGammaGaussianFisherInformation fi = new PoissonGammaGaussianFisherInformation(m, s, sampling);
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
		final int n = 1 + Math.max(0, settings.getSubDivisions());
		final double h = 1.0 / n;
		final double minExp = settings.getMinExponent();
		final double maxExp = settings.getMaxExponent();
		final double size = (maxExp - minExp) * n + 1;
		if (size > 100)
		{
			final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
			gd.addMessage("Number of exponents is " + Math.ceil(size) + ". OK to continue?");
			gd.showDialog();
			if (gd.wasCanceled())
				return null;
		}
		final TDoubleArrayList list = new TDoubleArrayList();
		for (int i = 0;; i++)
		{
			final double e = minExp + i * h;
			list.add(e);
			if (e >= settings.getMaxExponent())
				break;
		}
		return list.toArray();
	}

	private double[] getAlpha(final double[] photons, double[] exp, final BasePoissonFisherInformation fi, FIKey key)
	{
		final CameraType type = key.getType();
		final double[] alpha = new double[photons.length];
		if (!type.isFast())
		{
			final int[] index;

			// Try and load from the cache
			final PoissonFisherInformationData data = load(key);
			if (data != null)
			{
				// Dump the samples
				final TDoubleArrayList meanList = new TDoubleArrayList(data.getAlphaSampleCount());
				final TDoubleArrayList alphalist = new TDoubleArrayList(data.getAlphaSampleCount());
				for (final AlphaSample sample : data.getAlphaSampleList())
				{
					meanList.add(sample.getLog10Mean());
					alphalist.add(sample.getAlpha());
				}
				final double[] exp2 = meanList.toArray();
				final double[] alphas = alphalist.toArray();
				Sort.sortArrays(alphas, exp2, true);

				// Find any exponent not in the array
				final TIntArrayList list = new TIntArrayList(exp.length);
				for (int i = 0; i < exp.length; i++)
				{
					// Assume exp2 is sorted
					final int j = Arrays.binarySearch(exp2, exp[i]);
					if (j < 0)
						list.add(i); // Add to indices to compute
					else
						alpha[i] = alphas[j]; // Get alpha
				}
				index = list.toArray();
			}
			else
				// Compute all
				index = SimpleArrayUtils.newArray(alpha.length, 0, 1);

			if (index.length > 0)
			{
				IJ.showStatus("Computing " + getName(key));
				final int nThreads = Prefs.getThreads();
				if (es == null)
					es = Executors.newFixedThreadPool(nThreads);
				final Ticker ticker = Ticker.createStarted(new IJTrackProgress(), index.length, nThreads != 1);
				final int nPerThread = (int) Math.ceil((double) index.length / nThreads);
				final TurboList<Future<?>> futures = new TurboList<>(nThreads);
				for (int i = 0; i < index.length; i += nPerThread)
				{
					final int start = i;
					final int end = Math.min(index.length, i + nPerThread);
					futures.add(es.submit(new Runnable()
					{
						@Override
						public void run()
						{
							final BasePoissonFisherInformation fi2 = fi.clone();
							for (int ii = start; ii < end; ii++)
							{
								final int j = index[ii];
								alpha[j] = fi2.getAlpha(photons[j]);
								ticker.tick();
							}
						}
					}));
				}
				Utils.waitForCompletion(futures);
				ticker.stop();
				IJ.showStatus("");

				save(key, exp, alpha);
			}
		}
		else
			// Simple single threaded method.
			for (int i = 0; i < alpha.length; i++)
				alpha[i] = fi.getAlpha(photons[i]);
		return alpha;
	}

	private static double[] getFisherInformation(double[] alpha, double[] photons)
	{
		final double[] I = new double[photons.length];
		for (int i = 0; i < photons.length; i++)
			I[i] = alpha[i] / photons[i];
		return I;
	}

	private static double[] getAlpha(BasePoissonFisherInformation fi, double[] photons)
	{
		// Do not multi-thread as this is an interpolated/fast function
		final double[] rI = new double[photons.length];
		for (int i = 0; i < photons.length; i++)
			rI[i] = fi.getAlpha(photons[i]);
		return rI;
	}

	private static BasePoissonFisherInformation getInterpolatedPoissonFisherInformation(CameraType type, double[] logU,
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

	private static int getPointShape(int pointOption)
	{
		switch (pointOption)
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

	private static double[] limits(double[] limits, double[] f)
	{
		double min = limits[0];
		double max = limits[1];
		for (int i = 0; i < f.length; i++)
		{
			final double d = f[i];
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

	private static String cachePlot = "";
	private static int pointOption = 0;

	private static void plotFromCache()
	{
		// Build a list of curve stored in the cache
		final String[] names = new String[cache.size()];
		final PoissonFisherInformationData[] datas = new PoissonFisherInformationData[cache.size()];
		int c = 0;
		for (final Entry<FIKey, PoissonFisherInformationData> e : cache.entrySet())
		{
			names[c] = e.getKey().toString();
			datas[c] = e.getValue();
			c++;
		}

		final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addChoice("Fisher_information", names, cachePlot);
		gd.addChoice("Plot_point", POINT_OPTION, pointOption);
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		final int index = gd.getNextChoiceIndex();
		pointOption = gd.getNextChoiceIndex();
		cachePlot = names[index];

		final PoissonFisherInformationData data = datas[index];
		final FIKey key = new FIKey(data);

		c = 0;
		final double[] exp = new double[data.getAlphaSampleCount()];
		final double[] alpha = new double[data.getAlphaSampleCount()];
		for (final AlphaSample s : data.getAlphaSampleList())
		{
			exp[c] = s.getLog10Mean();
			alpha[c] = s.getAlpha();
			c++;
		}
		// Just in case
		//Sort.sortArrays(alpha, exp, true);

		// Test if we can use ImageJ support for a X log scale
		final boolean logScaleX = ((float) FastMath.pow(10, exp[0]) != 0);
		double[] x = exp;
		String xTitle = "log10(photons)";
		if (logScaleX)
		{
			final double[] photons = new double[exp.length];
			for (int i = 0; i < photons.length; i++)
				photons[i] = FastMath.pow(10, exp[i]);
			x = photons;
			xTitle = "photons";
		}

		// Get interpolation for alpha. Convert to base e.
		final double[] logU = exp.clone();
		final double scale = Math.log(10);
		for (int i = 0; i < logU.length; i++)
			logU[i] *= scale;
		final BasePoissonFisherInformation if1 = getInterpolatedPoissonFisherInformation(key.getType(), logU, alpha,
				null);

		// Interpolate with 5 points per sample for smooth curve
		final int n = 5;
		final TDoubleArrayList iexp = new TDoubleArrayList();
		final TDoubleArrayList iphotons = new TDoubleArrayList();
		for (int i = 1; i < exp.length; i++)
		{
			final int i_1 = i - 1;
			final double h = (exp[i] - exp[i_1]) / n;
			for (int j = 0; j < n; j++)
			{
				final double e = exp[i_1] + j * h;
				iexp.add(e);
				iphotons.add(FastMath.pow(10, e));
			}
		}
		iexp.add(exp[exp.length - 1]);
		iphotons.add(FastMath.pow(10, exp[exp.length - 1]));
		final double[] photons = iphotons.toArray();
		final double[] ix = (logScaleX) ? photons : iexp.toArray();
		final double[] ialpha1 = getAlpha(if1, photons);

		final int pointShape = getPointShape(pointOption);

		final String title = "Cached Relative Fisher Information";
		final Plot plot = new Plot(title, xTitle, "Noise coefficient (alpha)");
		plot.setLimits(x[0], x[x.length - 1], -0.05, 1.05);
		if (logScaleX)
			plot.setLogScaleX();
		plot.setColor(Color.blue);
		plot.addPoints(ix, ialpha1, Plot.LINE);
		// Option to show nodes
		if (pointShape != -1)
			plot.addPoints(x, alpha, pointShape);
		plot.setColor(Color.BLACK);
		plot.addLabel(0, 0, cachePlot);
		Utils.display(title, plot);
	}
}
