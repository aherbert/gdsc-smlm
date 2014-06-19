package gdsc.smlm.model;

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

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;

/**
 * Contains methods for generating models of a Point Spread Function
 */
public abstract class PSFModel
{
	private RandomDataGenerator rand;
	private double[] psf;
	private int x0min;
	private int x1min;
	private int x0max;
	private int x1max;

	public PSFModel()
	{
		rand = new RandomDataGenerator();
	}

	public PSFModel(RandomGenerator randomGenerator)
	{
		rand = new RandomDataGenerator(randomGenerator);
	}

	public PSFModel(RandomDataGenerator randomDataGenerator)
	{
		rand = randomDataGenerator;
	}

	/**
	 * Construct a PSF function on the provided data.
	 * 
	 * @param data
	 *            The data (can be null)
	 * @param width
	 *            The data width
	 * @param height
	 *            The data height
	 * @param sum
	 *            The integral
	 * @param x0
	 *            The centre in dimension 0
	 * @param x1
	 *            The centre in dimension 1
	 * @param x2
	 *            The centre in dimension 2
	 * @param poissonNoise
	 *            Add Poisson noise
	 * @return The total sum added to the image (useful when poissonNoise is added)
	 */
	public abstract double create3D(float[] data, final int width, final int height, final double sum, double x0,
			double x1, double x2, boolean poissonNoise);

	/**
	 * Construct a PSF function on the provided data.
	 * 
	 * @param data
	 *            The data (can be null)
	 * @param width
	 *            The data width
	 * @param height
	 *            The data height
	 * @param sum
	 *            The integral
	 * @param x0
	 *            The centre in dimension 0
	 * @param x1
	 *            The centre in dimension 1
	 * @param x2
	 *            The centre in dimension 2
	 * @param poissonNoise
	 *            Add Poisson noise
	 * @return The total sum added to the image (useful when poissonNoise is added)
	 */
	public abstract double create3D(double[] data, final int width, final int height, final double sum, double x0,
			double x1, double x2, boolean poissonNoise);

	/**
	 * Construct a PSF function on the provided data.
	 * 
	 * @param data
	 *            The data (can be null)
	 * @param width
	 *            The data width
	 * @param height
	 *            The data height
	 * @param sum
	 *            The integral
	 * @param x0
	 *            The centre in dimension 0
	 * @param x1
	 *            The centre in dimension 1
	 * @param poissonNoise
	 *            Add Poisson noise
	 * @return The total sum added to the image (useful when poissonNoise is added)
	 */
	public double create2D(float[] data, final int width, final int height, final double sum, double x0, double x1,
			boolean poissonNoise)
	{
		return create3D(data, width, height, sum, x0, x1, 0, poissonNoise);
	}

	/**
	 * Construct a PSF function on the provided data.
	 * 
	 * @param data
	 *            The data (can be null)
	 * @param width
	 *            The data width
	 * @param height
	 *            The data height
	 * @param sum
	 *            The integral
	 * @param x0
	 *            The centre in dimension 0
	 * @param x1
	 *            The centre in dimension 1
	 * @param poissonNoise
	 *            Add Poisson noise
	 * @return The total sum added to the image (useful when poissonNoise is added)
	 */
	public double create2D(double[] data, final int width, final int height, final double sum, double x0, double x1,
			boolean poissonNoise)
	{
		return create3D(data, width, height, sum, x0, x1, 0, poissonNoise);
	}

	/**
	 * @return The last drawn PSF
	 */
	public double[] getPSF()
	{
		return psf;
	}

	/**
	 * @return The minimum position in dimension 0 for the last drawn PSF
	 */
	public int getX0min()
	{
		return x0min;
	}

	/**
	 * @return The maximum position in dimension 0 for the last drawn PSF
	 */
	public int getX0max()
	{
		return x0max;
	}

	/**
	 * @return The minimum position in dimension 1 for the last drawn PSF
	 */
	public int getX1min()
	{
		return x1min;
	}

	/**
	 * @return The maximum position in dimension 1 for the last drawn PSF
	 */
	public int getX1max()
	{
		return x1max;
	}

	/**
	 * Insert the psf into the data
	 * 
	 * @param data
	 *            The input data (width*height)
	 * @param x0min
	 *            The minimum position to insert in dimension 0
	 * @param x1min
	 *            The minimum position to insert in dimension 1
	 * @param x0max
	 *            The maximum position to insert in dimension 0
	 * @param x1max
	 *            The maximum position to insert in dimension 1
	 * @param width
	 *            The width of the input data
	 * @param psf
	 *            The PSF data
	 * @param poissonNoise
	 * @return
	 */
	protected double insert(float[] data, int x0min, int x1min, int x0max, int x1max, int width, double[] psf,
			boolean poissonNoise)
	{
		final int x0range = x0max - x0min;
		final int x1range = x1max - x1min;

		if (x0range < 1 || x1range < 1)
		{
			this.psf = null;
			this.x0min = 0;
			this.x0max = 0;
			this.x1min = 0;
			this.x1max = 0;
			return 0;
		}

		this.psf = psf;
		this.x0min = x0min;
		this.x0max = x0max;
		this.x1min = x1min;
		this.x1max = x1max;

		if (poissonNoise)
		{
			for (int i = 0; i < psf.length; i++)
				if (psf[i] > 0)
					psf[i] = rand.nextPoisson(psf[i]);
		}

		// Insert the function into the input data
		for (int y = 0; y < x1range; y++)
		{
			// Locate the insert location
			int indexTo = (y + x1min) * width + x0min;
			int indexFrom = y * x0range;
			for (int x = 0; x < x0range; x++)
			{
				data[indexTo++] += psf[indexFrom++];
			}
		}

		double total = 0;
		for (double d : psf)
			total += d;
		return total;
	}

	/**
	 * Insert the psf into the data
	 * 
	 * @param data
	 *            The input data (width*height)
	 * @param x0min
	 *            The minimum position to insert in dimension 0
	 * @param x1min
	 *            The minimum position to insert in dimension 1
	 * @param x0max
	 *            The maximum position to insert in dimension 0
	 * @param x1max
	 *            The maximum position to insert in dimension 1
	 * @param width
	 *            The width of the input data
	 * @param psf
	 *            The PSF data
	 * @param poissonNoise
	 * @return
	 */
	protected double insert(double[] data, int x0min, int x1min, int x0max, int x1max, int width, double[] psf,
			boolean poissonNoise)
	{
		final int x0range = x0max - x0min;
		final int x1range = x1max - x1min;

		if (x0range < 1 || x1range < 1)
		{
			this.psf = null;
			this.x0min = 0;
			this.x0max = 0;
			this.x1min = 0;
			this.x1max = 0;
			return 0;
		}

		this.psf = psf;
		this.x0min = x0min;
		this.x0max = x0max;
		this.x1min = x1min;
		this.x1max = x1max;

		if (poissonNoise)
		{
			for (int i = 0; i < psf.length; i++)
				if (psf[i] > 0)
					psf[i] = rand.nextPoisson(psf[i]);
		}

		// Insert the function into the input data
		for (int y = 0; y < x1range; y++)
		{
			// Locate the insert location
			int indexTo = (y + x1min) * width + x0min;
			int indexFrom = y * x0range;
			for (int x = 0; x < x0range; x++)
			{
				data[indexTo++] += psf[indexFrom++];
			}
		}

		double total = 0;
		for (double d : psf)
			total += d;
		return total;
	}

	/**
	 * Remove the last added PSF from the data. This can be invoked after any call to draw a
	 * PSF into an input data array.
	 * 
	 * @param data
	 * @param width
	 * @param height
	 */
	public void erase(float[] data, int width, int height)
	{
		erase(data, width, height, psf, x0min, x0max, x1min, x1max);
	}

	/**
	 * Remove the last added PSF from the data. This can be invoked after any call to draw a
	 * PSF into an input data array.
	 * 
	 * @param data
	 * @param width
	 * @param height
	 */
	public void erase(double[] data, int width, int height)
	{
		erase(data, width, height, psf, x0min, x0max, x1min, x1max);
	}

	/**
	 * Remove the PSF from the data. Can be invoked using a saved copy of the PSF previously drawn by the model obtained
	 * from the appropriate get() methods.
	 * 
	 * @param data
	 * @param width
	 * @param height
	 * @param psf
	 * @param x0min
	 * @param x0max
	 * @param x1min
	 * @param x1max
	 */
	public void erase(float[] data, int width, int height, double[] psf, int x0min, int x0max, int x1min, int x1max)
	{
		if (psf == null)
			return;

		final int x0range = x0max - x0min;
		final int x1range = x1max - x1min;

		// min should always be less than max
		if (x0range < 1)
			throw new IllegalArgumentException("Dimension 0 range not within data bounds");
		if (x1range < 1)
			throw new IllegalArgumentException("Dimension 1 range not within data bounds");

		// Remove from the input data
		for (int y = 0; y < x1range; y++)
		{
			// Locate the insert location
			int indexTo = (y + x1min) * width + x0min;
			int indexFrom = y * x0range;
			for (int x = 0; x < x0range; x++)
			{
				data[indexTo++] -= psf[indexFrom++];
			}
		}
	}

	/**
	 * Remove the PSF from the data. Can be invoked using a saved copy of the PSF previously drawn by the model obtained
	 * from the appropriate get() methods.
	 * 
	 * @param data
	 * @param width
	 * @param height
	 * @param psf
	 * @param x0min
	 * @param x0max
	 * @param x1min
	 * @param x1max
	 */
	public void erase(double[] data, int width, int height, double[] psf, int x0min, int x0max, int x1min, int x1max)
	{
		if (psf == null)
			return;

		final int x0range = x0max - x0min;
		final int x1range = x1max - x1min;

		// min should always be less than max
		if (x0range < 1)
			throw new IllegalArgumentException("Dimension 0 range not within data bounds");
		if (x1range < 1)
			throw new IllegalArgumentException("Dimension 1 range not within data bounds");

		// Remove from the input data
		for (int y = 0; y < x1range; y++)
		{
			// Locate the insert location
			int indexTo = (y + x1min) * width + x0min;
			int indexFrom = y * x0range;
			for (int x = 0; x < x0range; x++)
			{
				data[indexTo++] -= psf[indexFrom++];
			}
		}
	}
}
