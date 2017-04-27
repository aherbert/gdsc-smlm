package gdsc.smlm.function.gaussian;

import org.apache.commons.math3.distribution.NormalDistribution;

import gdsc.core.utils.Sort;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Compute the overlap between 2D Gaussian functions.
 * <p>
 * Given an input 2D Gaussian a region is created that covers a range of the function (relative to the SD). The square
 * region is masked using the expected sum of the function within the range. The overlap of other functions within this
 * region can be computed.
 */
public class GaussianOverlapAnalysis
{
	private static final double[] p;
	static
	{
		p = new double[20];
		final NormalDistribution d = new NormalDistribution();
		for (int x = 1; x < p.length; x++)
			p[x] = d.probability(-x, x);
	}

	private final int flags;
	private final double[] params0;
	private final double range;

	private final int maxx, maxy, size;
	private final double centrex, centrey;

	private double[] data = null, overlap = null;
	private boolean[] mask = null;

	/**
	 * Create a new overlap analysis object
	 * 
	 * @param flags
	 *            The flags describing the Gaussian2DFunction function (see GaussianFunctionFactory)
	 * @param params
	 *            The parameters for the Gaussian (assumes a single peak)
	 * @param range
	 *            The range over which to compute the function (factor of the standard deviation)
	 */
	public GaussianOverlapAnalysis(int flags, double[] params, double range)
	{
		this.flags = flags;
		this.params0 = params.clone();
		this.range = range;

		maxx = 2 * ((int) Math.ceil(params[Gaussian2DFunction.X_SD] * range)) + 1;
		maxy = 2 * ((int) Math.ceil(params[Gaussian2DFunction.Y_SD] * range)) + 1;
		size = maxx * maxy;
		// We will sample the Gaussian at integer intervals, i.e. on a pixel grid.
		// Pixels centres should be at 0.5,0.5. So if we want to draw a Gauss 
		// centred in the middle of a pixel we need to adjust each centre 
		centrex = maxx * 0.5 - 0.5;
		centrey = maxy * 0.5 - 0.5;
	}

	/**
	 * Add the Gaussian function data to the overlap region. This is the region that contains the input function within
	 * the range defined in the constructor.
	 * 
	 * @param flags
	 *            The flags describing the Gaussian2DFunction function (see GaussianFunctionFactory)
	 * @param params
	 *            The parameters for the Gaussian (can be multiple peaks)
	 */
	public void add(int flags, double[] params)
	{
		add(flags, params, false);
	}

	/**
	 * Add the Gaussian function data to the overlap region. This is the region that contains the input function within
	 * the range defined in the constructor.
	 * <p>
	 * The square region is masked using the expected sum of the function within the range. The overlap of other
	 * functions within this masked region can be computed, or within the square region.
	 *
	 * @param flags
	 *            The flags describing the Gaussian2DFunction function (see GaussianFunctionFactory)
	 * @param params
	 *            The parameters for the Gaussian (can be multiple peaks)
	 * @param withinMask
	 *            Set to true to only compute the overlap within the mask. This effects the computation of the weighted
	 *            background (see {@link #getWeightedbackground()}.
	 */
	public void add(int flags, double[] params, boolean withinMask)
	{
		// Note: When computing the overlap no adjustment is made for sampling on a pixel grid.
		// This is OK since the method will be used to evaluate the overlap between Gaussians that have
		// been fit using the functions.

		if (data == null)
		{
			// Initialise the input function
			data = new double[size];
			overlap = new double[size];
			Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, maxx, maxy, this.flags);

			// Note that the position should be in the centre of the sample region
			final double cx = params0[Gaussian2DFunction.X_POSITION];
			final double cy = params0[Gaussian2DFunction.Y_POSITION];
			params0[Gaussian2DFunction.X_POSITION] = centrex;
			params0[Gaussian2DFunction.Y_POSITION] = centrey;

			f.initialise(params0);
			final int[] indices = new int[size];
			for (int k = 0; k < size; k++)
			{
				final double v = f.eval(k);
				data[k] = v;
				indices[k] = k;
			}
			// Reset
			params0[Gaussian2DFunction.X_POSITION] = cx;
			params0[Gaussian2DFunction.Y_POSITION] = cy;

			// Compute the expected sum in the range 
			final double expected = getArea(range) * params0[Gaussian2DFunction.SIGNAL];

			Sort.sort(indices, data);
			double sum = 0, last = 0;
			boolean useMask = false;
			mask = new boolean[data.length];
			for (int i = 0; i < data.length; i++)
			{
				final double v = data[indices[i]];
				// Note: We track the value since the Gaussian is symmetric and we want to include
				// all pixels with the same value
				final double newSum = sum + v;
				if (newSum >= expected && last != v)
				{
					// This is a new value that takes us over the expected signal
					useMask = true;
					break;
				}
				sum = newSum;
				mask[indices[i]] = true;
				last = v;
			}
			if (!useMask)
				mask = null;
		}

		// Add the function to the overlap
		final int nPeaks = params.length / 6;
		Gaussian2DFunction f = GaussianFunctionFactory.create2D(nPeaks, maxx, maxy, flags);
		params = params.clone();
		for (int n = 0; n < nPeaks; n++)
		{
			params[n * 6 + Gaussian2DFunction.X_POSITION] += centrex - params0[Gaussian2DFunction.X_POSITION];
			params[n * 6 + Gaussian2DFunction.Y_POSITION] += centrey - params0[Gaussian2DFunction.Y_POSITION];
		}
		f.initialise(params);
		if (mask == null || !withinMask)
		{
			for (int k = 0; k < size; k++)
			{
				overlap[k] += f.eval(k);
			}
		}
		else
		{
			for (int k = 0; k < size; k++)
			{
				if (mask[k])
				{
					overlap[k] += f.eval(k);
				}
			}
		}

		//		// Debug
		//		double[] combined = data.clone();
		//		for (int k = 0; k < size; k++)
		//			combined[k] += overlap[k];
		//		gdsc.core.ij.Utils.display("Spot i", data, maxx, maxx);
		//		gdsc.core.ij.Utils.display("Overlap", overlap, maxx, maxx);
		//		gdsc.core.ij.Utils.display("Combined", combined, maxx, maxx);
		//		System.out.printf("Signal %.2f, sd %.2f, overlap = %s\n", params0[Gaussian2DFunction.SIGNAL],
		//				params0[Gaussian2DFunction.X_SD], Arrays.toString(params));
	}

	/**
	 * Get the overlap data.
	 * <p>
	 * Computes the sum of the central function, and the sum of the overlap
	 * 
	 * @return The data [sum f1, sum overlap]
	 */
	public double[] getOverlapData()
	{
		double[] result = new double[2];
		if (overlap != null)
		{
			double sumF = 0;
			double sumO = 0;
			if (mask == null)
			{
				for (int k = 0; k < size; k++)
				{
					sumF += data[k];
					sumO += overlap[k];
				}
			}
			else
			{
				for (int k = 0; k < size; k++)
				{
					if (mask[k])
					{
						sumF += data[k];
						sumO += overlap[k];
					}
				}
			}

			result[0] = sumF;
			result[1] = sumO;
		}
		return result;
	}

	/**
	 * Get the weighted background
	 * <p>
	 * Computes the convolution of the central function and the overlap for the central pixel of the region. This is an
	 * estimate of the background contributed to the region by overlapping functions.
	 * <p>
	 * The result of this function is effected by how the overlap was computed, either within the mask or within the
	 * entire square region (see {@link #add(int, double[], boolean)})
	 * 
	 * @return The weighted background
	 */
	public double getWeightedbackground()
	{
		if (overlap != null)
		{
			double sum = 0;
			double sumw = 0;
			final double norm = 1.0 / params0[Gaussian2DFunction.SIGNAL];
			//double[] combined = new double[size];
			for (int k = 0; k < size; k++)
			{
				final double v = data[k] * norm;
				sum += overlap[k] * v;
				sumw += v;
			}
			//gdsc.core.ij.Utils.display("O Spot i", data, maxx, maxx);
			//gdsc.core.ij.Utils.display("O Spot overlap", overlap, maxx, maxx);
			return sum / sumw;
		}
		return 0;
	}

	/**
	 * Get the probability of a standard Gaussian within the given range x, i.e. P(-x < X <= x)
	 * 
	 * @param x
	 * @return The probability (0-1)
	 */
	public static double getArea(double x)
	{
		if ((int) x == x && x < p.length)
			return p[(int) x];
		final NormalDistribution d = new NormalDistribution();
		return d.probability(-x, x);
	}
}
