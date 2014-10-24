package gdsc.smlm.utils;

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
 * Apply a window function to reduce edge artifacts
 */
public class ImageWindow
{
	public enum WindowFunction {
		HANNING, COSINE, TUKEY
	}

	// Allow cached window weights
	private double[] wx = null;
	private double[] wy = null;
	private WindowFunction windowFunction = null;

	/**
	 * Apply a window function to reduce edge artifacts.
	 * <p>
	 * Applied as two 1-dimensional window functions. Faster than the nonseparable form but has
	 * direction dependent corners.
	 * <p>
	 * Instance method allows caching the weight matrices. 
	 * 
	 * @param image
	 * @param maxx
	 * @param maxy
	 * @param windowFunction
	 * @return
	 */
	public float[] applySeperable(float[] image, final int maxx,
			final int maxy, WindowFunction windowFunction)
	{
		if (this.windowFunction != windowFunction || wx == null
				|| wx.length != maxx || wy == null || wy.length != maxy)
		{
			switch (windowFunction)
			{
			case HANNING:
				wx = hanning(maxx);
				wy = hanning(maxy);
				break;
			case COSINE:
				wx = cosine(maxx);
				wy = cosine(maxy);
				break;
			case TUKEY:
			default:
				wx = tukey(maxx, ALPHA);
				wy = tukey(maxy, ALPHA);
				break;
			}
		}

		float[] data = new float[image.length];

		for (int y = 0, i = 0; y < maxy; y++)
		{
			for (int x = 0; x < maxx; x++, i++)
			{
				data[i] = (float) (image[i] * wx[x] * wy[y]);
			}
		}

		return data;
	}

	/**
	 * Apply a window function to reduce edge artifacts.
	 * <p>
	 * Applied as two 1-dimensional window functions. Faster than the nonseparable form but has
	 * direction dependent corners.
	 * 
	 * @param image
	 * @param maxx
	 * @param maxy
	 * @param windowFunction
	 * @return
	 */
	public static float[] applyWindowSeparable(float[] image, final int maxx,
			final int maxy, WindowFunction windowFunction)
	{
		double[] wx = null;
		double[] wy = null;

		switch (windowFunction)
		{
		case HANNING:
			wx = hanning(maxx);
			wy = hanning(maxy);
			break;
		case COSINE:
			wx = cosine(maxx);
			wy = cosine(maxy);
			break;
		case TUKEY:
			wx = tukey(maxx, ALPHA);
			wy = tukey(maxy, ALPHA);
			break;
		}

		if (wx == null)
			return image;

		float[] data = new float[image.length];

		for (int y = 0, i = 0; y < maxy; y++)
		{
			for (int x = 0; x < maxx; x++, i++)
			{
				data[i] = (float) (image[i] * wx[x] * wy[y]);
			}
		}

		return data;
	}

	/**
	 * Apply a window function to reduce edge artifacts
	 * <p>
	 * Applied as a nonseparable form.
	 * 
	 * @param image
	 * @param maxx
	 * @param maxy
	 * @param windowFunction
	 * @return
	 */
	public static float[] applyWindow(float[] image, final int maxx,
			final int maxy, WindowFunction windowFunction)
	{
		WindowMethod wf = null;
		switch (windowFunction)
		{
		case HANNING:
			wf = instance.new Hanning();
			break;
		case COSINE:
			wf = instance.new Cosine();
			break;
		case TUKEY:
			wf = instance.new Tukey(ALPHA);
		}

		if (wf == null)
			return image;

		float[] data = new float[image.length];

		double cx = maxx * 0.5;
		double cy = maxy * 0.5;
		double maxDistance = Math.sqrt(maxx * maxx + maxy * maxy);

		for (int y = 0, i = 0; y < maxy; y++)
		{
			for (int x = 0; x < maxx; x++, i++)
			{
				double distance = Math.sqrt((x - cx) * (x - cx) + (y - cy)
						* (y - cy));
				double w = wf.weight(0.5 - (distance / maxDistance));
				data[i] = (float) (image[i] * w);
			}
		}

		return data;
	}

	private static ImageWindow instance = new ImageWindow();
	private static double ALPHA = 0.5;

	private interface WindowMethod
	{
		/**
		 * Return the weight for the window at a fraction of the distance from the edge of the
		 * window.
		 * 
		 * @param fractionDistance
		 *            (range 0-1)
		 * @return
		 */
		double weight(double fractionDistance);
	}

	private class Hanning implements WindowMethod
	{
		public double weight(double fractionDistance)
		{
			return 0.5 * (1 - Math.cos(Math.PI * 2 * fractionDistance));
		}
	}

	private class Cosine implements WindowMethod
	{
		public double weight(double fractionDistance)
		{
			return Math.sin(Math.PI * fractionDistance);
		}
	}

	private class Tukey implements WindowMethod
	{
		final double alpha;

		public Tukey(double alpha)
		{
			this.alpha = alpha;
		}

		public double weight(double fractionDistance)
		{
			if (fractionDistance < alpha / 2)
				return 0.5 * (1 + Math.cos(Math.PI
						* (2 * fractionDistance / alpha - 1)));
			if (fractionDistance > 1 - alpha / 2)
				return 0.5 * (1 + Math.cos(Math.PI
						* (2 * fractionDistance / alpha - 2 / alpha + 1)));
			return 1;
		}
	}

	private static double[] createWindow(WindowMethod wf, int N)
	{
		double N_1 = N - 1;
		double[] w = new double[N];
		for (int i = 0; i < N; i++)
			w[i] = wf.weight(i / N_1);
		return w;
	}

	private static double[] hanning(int N)
	{
		return createWindow(instance.new Hanning(), N);
	}

	private static double[] cosine(int N)
	{
		return createWindow(instance.new Cosine(), N);
	}

	private static double[] tukey(int N, double alpha)
	{
		return createWindow(instance.new Tukey(alpha), N);
	}
}
