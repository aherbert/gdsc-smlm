package gdsc.smlm.results;

import gdsc.smlm.fitting.function.Gaussian2DFunction;

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
 * Specifies a peak fitting result
 */
public class PeakResult implements Comparable<PeakResult>
{
	/**
	 * Adjustment factor for the precision calculation for an EMCCD camera.
	 * <p>
	 * The Mortensen paper uses a factor of 2 for EMCCD cameras. This value can be calibrated experimentally. The value
	 * for the camera (Photometrics Evolve2) at the GDSC is 1.83.
	 */
	public static final double F = 1.83;

	public int peak;
	public int origX;
	public int origY;
	public float origValue;
	public double error;
	public float noise;
	public float[] params;
	public float[] paramsStdDev;

	public PeakResult(int peak, int origX, int origY, float origValue, double error, float noise, float[] params,
			float[] paramsStdDev)
	{
		this.peak = peak;
		this.origX = origX;
		this.origY = origY;
		this.origValue = origValue;
		this.error = error;
		this.noise = noise;
		this.params = params;
		this.paramsStdDev = paramsStdDev;
	}

	/**
	 * Simple constructor to create a result with location, width and strength
	 * 
	 * @param x
	 * @param y
	 * @param sd
	 * @param signal
	 */
	public PeakResult(float x, float y, float sd, int signal)
	{
		// TODO Auto-generated constructor stub
		origX = (int) Math.round(x);
		origY = (int) Math.round(y);
		params = new float[7];
		params[Gaussian2DFunction.X_POSITION] = x;
		params[Gaussian2DFunction.Y_POSITION] = y;
		params[Gaussian2DFunction.X_SD] = params[Gaussian2DFunction.Y_SD] = sd;
		params[Gaussian2DFunction.AMPLITUDE] = (float) (signal / (Math.PI * 2 * sd * sd));
	}

	/**
	 * Get the signal strength (i.e. the volume under the Gaussian peak, amplitude * 2 * pi * sx * sy)
	 * 
	 * @return The signal of the first peak
	 */
	public float getSignal()
	{
		return getSignal(params);
	}

	/**
	 * Get the signal strength for the nth peak (i.e. the volume under the Gaussian peak, amplitude * 2 * pi * sx * sy)
	 * 
	 * @param peakId
	 *            The peak number
	 * @return The signal of the nth peak
	 */
	public float getSignal(int peakId)
	{
		return getSignal(peakId, params);
	}

	/**
	 * Get the signal strength (i.e. the volume under the Gaussian peak, amplitude * 2 * pi * sx * sy)
	 * 
	 * @param params
	 *            The peak parameters
	 * @return The signal of the first peak
	 */
	public static float getSignal(float[] params)
	{
		return (float) (params[Gaussian2DFunction.AMPLITUDE] * 2 * Math.PI * params[Gaussian2DFunction.X_SD] * params[Gaussian2DFunction.Y_SD]);
	}

	/**
	 * Get the signal strength for the nth peak (i.e. the volume under the Gaussian peak, amplitude * 2 * pi * sx * sy)
	 * 
	 * @param peakId
	 *            The peak number
	 * @param params
	 *            The peak parameters
	 * @return The signal of the nth peak
	 */
	public static float getSignal(int peakId, float[] params)
	{
		if (peakId * 6 + 6 >= params.length)
			return 0;
		return (float) (params[peakId * 6 + Gaussian2DFunction.AMPLITUDE] * 2 * Math.PI *
				params[peakId * 6 + Gaussian2DFunction.X_SD] * params[peakId * 6 + Gaussian2DFunction.Y_SD]);
	}

	/**
	 * Calculate the localisation precision. Uses the Mortensen method for an EMCCD camera
	 * 
	 * @param a
	 *            The size of the pixels in nm
	 * @return The location precision in nm of the first peak
	 */
	public double getPrecision(double a, double gain)
	{
		// Get peak standard deviation in nm. Just use the average of the X & Y.
		final double s = a * getSD();
		final double N = getSignal();
		return getPrecision(a, s, N / gain, noise / gain);
	}

	/**
	 * Calculate the localisation precision. Uses the Mortensen method for an EMCCD camera
	 * 
	 * @param a
	 *            The size of the pixels in nm
	 * @param s
	 *            The peak standard deviation in nm
	 * @param N
	 *            The peak signal in photons
	 * @param b
	 *            The background noise in photons
	 * @return The location precision in nm in each dimension (X/Y)
	 */
	public static double getPrecision(double a, double s, double N, double b)
	{
		final double a2 = a * a;
		final double sa2 = s * s + a2 / 12.0; // Adjustment for square pixels
		// 16 / 9 = 1.7777777778
		// 8 * pi = 25.13274123
		return Math.sqrt(F * (sa2 / N) * (1.7777777778 + (25.13274123 * sa2 * b * b) / (N * a2)));
	}

	public int compareTo(PeakResult o)
	{
		// Sort by peak number: Ascending
		if (peak == o.peak)
		{
			// Sort by peak height: Descending
			if (params[Gaussian2DFunction.AMPLITUDE] > o.params[Gaussian2DFunction.AMPLITUDE])
				return -1;
			if (params[Gaussian2DFunction.AMPLITUDE] < o.params[Gaussian2DFunction.AMPLITUDE])
				return 1;
			return 0;
		}
		return peak - o.peak;
	}

	/**
	 * @return The average peak standard deviation in the X and Y dimension
	 */
	public float getSD()
	{
		return (Math.abs(params[Gaussian2DFunction.X_SD]) + Math.abs(params[Gaussian2DFunction.Y_SD])) * 0.5f;
	}

	/**
	 * @return The background for the first peak
	 */
	public float getBackground()
	{
		return params[Gaussian2DFunction.BACKGROUND];
	}

	/**
	 * @return The amplitude for the first peak
	 */
	public float getAmplitude()
	{
		return params[Gaussian2DFunction.AMPLITUDE];
	}

	/**
	 * @return The angle for the first peak
	 */
	public float getAngle()
	{
		return params[Gaussian2DFunction.ANGLE];
	}

	/**
	 * @return The x position for the first peak
	 */
	public float getXPosition()
	{
		return params[Gaussian2DFunction.X_POSITION];
	}

	/**
	 * @return The y position for the first peak
	 */
	public float getYPosition()
	{
		return params[Gaussian2DFunction.Y_POSITION];
	}

	/**
	 * @return The x-dimension standard deviation for the first peak
	 */
	public float getXSD()
	{
		return params[Gaussian2DFunction.X_SD];
	}

	/**
	 * @return The y-dimension standard deviation for the first peak
	 */
	public float getYSD()
	{
		return params[Gaussian2DFunction.Y_SD];
	}

	/**
	 * @return The last time frame that this result corresponds to
	 */
	public int getEndFrame()
	{
		return peak;
	}

	/**
	 * @return The results identifier
	 */
	public int getId()
	{
		return 0;
	}
}
