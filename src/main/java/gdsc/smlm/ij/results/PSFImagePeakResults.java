package gdsc.smlm.ij.results;

import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.results.PeakResult;

import java.awt.Rectangle;
import java.util.Collection;

import org.apache.commons.math3.util.FastMath;

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
 * Draws the fit results using the Gaussian PSF to an ImageJ image
 */
public class PSFImagePeakResults extends IJImagePeakResults
{
	private boolean fixedWidth = false;
	private float psfWidth = 0f;
	private boolean calculatedPrecision = false;

	private double nmPerPixel = 100.0;
	private double gain = 1;
	private boolean emCCD = true;

	// Multiplication factors and variables for plotting the fixed Gaussian
	private double[] fixedParams = null;

	/**
	 * @param title
	 *            Title of the image (appended with a suffix)
	 * @param bounds
	 *            Define the bounding rectangle of the image coordinates. Any results outside this will not be
	 *            displayed.
	 * @param scale
	 */
	public PSFImagePeakResults(String title, Rectangle bounds, float scale)
	{
		super(title, bounds, scale);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.IJImagePeakResults#checkDisplayFlags()
	 */
	@Override
	protected void preBegin()
	{
		// Flags should be OK

		// Cache the nmPerPixel and gain
		if (calibration != null)
		{
			nmPerPixel = calibration.nmPerPixel;
			gain = calibration.gain;
			emCCD = calibration.emCCD;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.results.IJImagePeakResults#add(int, float, float, float)
	 */
	@Override
	public void add(int peak, float x, float y, float v)
	{
		throw new RuntimeException(
				"This method is not supported. Some PSF images require the PSF parameters for amplitude and angle.");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.results.IJImagePeakResults#add(int, int, int, float, double, float, float[], float[])
	 */
	public void add(int peak, int origX, int origY, float origValue, double error, float noise, float[] params,
			float[] paramsDev)
	{
		if (!imageActive)
			return;
		addPeak(peak, origX, origY, origValue, error, noise, params, paramsDev);
		updateImage();
	}

	private void addPeak(int peak, int origX, int origY, float origValue, double chiSquared, float noise,
			float[] params, float[] paramsDev)
	{
		float x = (params[3] - bounds.x) * scale;
		float y = (params[4] - bounds.y) * scale;

		// Check bounds
		if (x < 0 || x > imageWidth || y < 0 || y > imageHeight)
			return;

		checkAndUpdateToFrame(peak);

		// Initialise for a free Gaussian function:
		//   f(x,y) = A exp(-(a(x-x0)(x-x0) + 2b(x-x0)(y-y0) + c(y-y0)(y-y0)))
		// See: http://en.wikipedia.org/wiki/Gaussian_function#Two-dimensional_Gaussian_function

		final float amplitude = ((displayFlags & DISPLAY_SIGNAL) != 0) ? PeakResult.getAmplitude(params) : 1;

		final double[] psfParams;
		if (fixedWidth)
		{
			psfParams = fixedParams;
		}
		else
		{
			// Precalculate multiplication factors
			final double t, sx, sy;
			if (calculatedPrecision && nmPerPixel > 0)
			{
				t = 0.0;
				final double N = params[Gaussian2DFunction.SIGNAL] / gain;
				final double s = (params[Gaussian2DFunction.X_SD] + params[Gaussian2DFunction.Y_SD]) * 0.5 * nmPerPixel;
				final double precision = PeakResult.getPrecision(nmPerPixel, s, N, noise / gain, emCCD);
				sx = sy = (precision / nmPerPixel);
			}
			else
			{
				t = params[Gaussian2DFunction.ANGLE];
				sx = params[Gaussian2DFunction.X_SD];
				sy = params[Gaussian2DFunction.Y_SD];
			}
			psfParams = setPSFParameters(t, sx, sy, new double[5]);
		}

		final double a = psfParams[0];
		final double b = psfParams[1];
		final double c = psfParams[2];
		final double width = psfParams[3];
		final double height = psfParams[4];

		// Use 0.5 offset to centre the value in the middle of each pixel 
		x -= 0.5 / scale;
		y -= 0.5 / scale;

		int xmin = (int) Math.floor(x - width * scale);
		int xmax = (int) Math.ceil(x + width * scale);
		int ymin = (int) Math.floor(y - height * scale);
		int ymax = (int) Math.ceil(y + height * scale);

		// Clip range
		xmin = FastMath.max(xmin, 0);
		xmax = (int) FastMath.min(xmax, xlimit);
		ymin = FastMath.max(ymin, 0);
		ymax = (int) FastMath.min(ymax, ylimit);

		// Compute Gaussian PSF
		final int[] index = new int[(xmax - xmin + 1) * (ymax - ymin + 1)];
		final float[] value = new float[index.length];
		int i = 0;
		for (int y0 = ymin; y0 <= ymax; y0++)
			for (int x0 = xmin; x0 <= xmax; x0++)
			{
				int ii = y0 * imageWidth + x0;
				index[i] = ii;
				final float dx = (x0 - x) / scale;
				final float dy = (y0 - y) / scale;
				value[i] = (float) (amplitude * FastMath.exp(a * dx * dx + b * dx * dy + c * dy * dy));
				i++;
			}

		// Now add the values to the configured indices
		synchronized (data)
		{
			size++;
			while (i-- > 0)
			{
				data[index[i]] += value[i];
			}
		}
	}

	public double[] setPSFParameters(double t, double sx, double sy, double[] params)
	{
		double a, b, c;
		double height, width;

		if (t == 0)
		{
			// sin(0) == 0
			// cos(0) == 1
			a = (-1.0 / (2.0 * sx * sx));
			b = 0.0;
			c = (-1.0 / (2.0 * sy * sy));

			// Calculate the range for the PSF as 3 sigma.
			width = 3.0 * sx;
			height = 3.0 * sy;
		}
		else
		{
			a = -(Math.cos(t) * Math.cos(t) / (2.0 * sx * sx) + Math.sin(t) * Math.sin(t) / (2.0 * sy * sy));
			b = -(-Math.sin(2.0 * t) / (2.0 * sx * sx) + Math.sin(2.0 * t) / (2.0 * sy * sy));
			c = -(Math.sin(t) * Math.sin(t) / (2.0 * sx * sx) + Math.cos(t) * Math.cos(t) / (2.0 * sy * sy));

			// Note that the Gaussian2DFitter returns the angle of the major axis (sx) relative to the x-axis.
			// The angle is in the range -pi/2 to pi/2

			// The width and height for the range to be plotted can be derived from the general parametric 
			// form of the ellipse.
			// See: http://en.wikipedia.org/wiki/Ellipse#General_parametric_form

			// Ensure the angle is the correct range (0 to pi)
			if (t < 0)
				t += Math.PI;

			final double phi = t; // Angle between x-axis and major axis of ellipse
			final double t1 = -t; // Angle around the ellipse from 0 to 2pi for the x-axis
			final double t2 = t1 + Math.PI / 2; // Angle around the ellipse from 0 to 2pi for the y-axis

			// Calculate the size of the ellipse at 3 sigma
			width = Math.abs(3.0 * (sx * Math.cos(t1) * Math.cos(phi) - sy * Math.sin(t1) * Math.sin(phi)));
			height = Math.abs(3.0 * (sx * Math.cos(t2) * Math.sin(phi) + sy * Math.sin(t2) * Math.cos(phi)));
		}

		params[0] = a;
		params[1] = b;
		params[2] = c;
		params[3] = width;
		params[4] = height;
		return params;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#addAll(java.util.Collection)
	 */
	public void addAll(Collection<PeakResult> results)
	{
		if (!imageActive)
			return;

		// TODO - Make this more efficient. It could use worker threads to increase speed.
		int i = 0;
		for (PeakResult result : results)
		{
			addPeak(result.peak, result.origX, result.origY, result.origValue, result.error, result.noise,
					result.params, result.paramsStdDev);
			if (++i % 64 == 0)
			{
				updateImage();
				if (!imageActive)
					return;
			}
		}
		updateImage();
	}

	/**
	 * @return the width for a fixed-width Gaussian
	 */
	public float getWidth()
	{
		return psfWidth;
	}

	/**
	 * Set the width of a fixed-width PSF. This is before scale adjustment so is provided in terms of the original
	 * fitted data.
	 * 
	 * @param width
	 *            the width to set for a fixed-width Gaussian
	 */
	public void setWidth(float width)
	{
		this.psfWidth = width;

		if (width > 0)
		{
			fixedWidth = true;
			fixedParams = setPSFParameters(0.0, width, width, new double[5]);
		}
		else
		{
			fixedWidth = false;
			fixedParams = null;
		}
	}

	/**
	 * @return if true plot the width of the PSF using the calculated precision
	 */
	public boolean isCalculatedPrecision()
	{
		return calculatedPrecision;
	}

	/**
	 * Set to true to plot the width of the PSF using the calculated precision
	 * 
	 * @param calculatedPrecision
	 */
	public void setCalculatedPrecision(boolean calculatedPrecision)
	{
		this.calculatedPrecision = calculatedPrecision;
	}
}
