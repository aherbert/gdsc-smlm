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
import org.apache.commons.math3.special.Erf;

/**
 * Contains methods for generating models of a Point Spread Function using a Gaussian approximation
 */
public class GaussianPSFModel extends PSFModel
{
	private double zeroS0, zeroS1;
	private double s0;
	private double s1;
	private double zDepth = 0;

	/**
	 * @param s0
	 *            The Gaussian standard deviation dimension 0
	 * @param s1
	 *            The Gaussian standard deviation dimension 1
	 */
	public GaussianPSFModel(double s0, double s1)
	{
		super();
		this.zeroS0 = s0;
		this.zeroS1 = s1;
	}

	/**
	 * @param randomGenerator
	 * @param s0
	 *            The Gaussian standard deviation dimension 0
	 * @param s1
	 *            The Gaussian standard deviation dimension 1
	 */
	public GaussianPSFModel(RandomGenerator randomGenerator, double s0, double s1)
	{
		super(randomGenerator);
		this.zeroS0 = s0;
		this.zeroS1 = s1;
	}

	/**
	 * @param randomGenerator
	 * @param s0
	 *            The Gaussian standard deviation dimension 0
	 * @param s1
	 *            The Gaussian standard deviation dimension 1
	 * @param zDepth
	 *            the Z-depth where the 3D PSF is 1.5x the width (1.5 x FWHM)
	 */
	public GaussianPSFModel(RandomGenerator randomGenerator, double s0, double s1, double zDepth)
	{
		super(randomGenerator);
		this.zeroS0 = s0;
		this.zeroS1 = s1;
		setzDepth(zDepth);
	}

	/**
	 * @param randomDataGenerator
	 * @param s0
	 *            The Gaussian standard deviation dimension 0
	 * @param s1
	 *            The Gaussian standard deviation dimension 1
	 */
	public GaussianPSFModel(RandomDataGenerator randomDataGenerator, double s0, double s1)
	{
		super(randomDataGenerator);
		this.zeroS0 = s0;
		this.zeroS1 = s1;
	}

	/**
	 * @param randomDataGenerator
	 * @param s0
	 *            The Gaussian standard deviation dimension 0
	 * @param s1
	 *            The Gaussian standard deviation dimension 1
	 * @param zDepth
	 *            the Z-depth where the 3D PSF is twice the width (2 x FWHM)
	 */
	public GaussianPSFModel(RandomDataGenerator randomDataGenerator, double s0, double s1, double zDepth)
	{
		super(randomDataGenerator);
		this.zeroS0 = s0;
		this.zeroS1 = s1;
		setzDepth(zDepth);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.PSFModel#create3D(float[], int, int, double, double, double, double, boolean)
	 */
	public double create3D(float[] data, final int width, final int height, final double sum, double x0, double x1,
			double x2, boolean poissonNoise)
	{
		if (sum == 0)
			return 0;
		final double scale = createWidthScale(x2);
		try
		{
			final double d = gaussian2D(data, width, height, sum, x0, x1, scale * zeroS0, scale * zeroS1, poissonNoise);
			//			if (d == 0)
			//			{
			//				System.out.printf("No data inserted: %f @ %f %f %f (%f x %f)\n", sum, x0, x1, x2, scale * zeroS0,
			//						scale * zeroS1);
			//			}
			return d;
		}
		catch (IllegalArgumentException e)
		{
			//System.out.println(e.getMessage());
			return 0;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.PSFModel#create3D(double[], int, int, double, double, double, double, boolean)
	 */
	public double create3D(double[] data, final int width, final int height, final double sum, double x0, double x1,
			double x2, boolean poissonNoise)
	{
		if (sum == 0)
			return 0;
		final double scale = createWidthScale(x2);
		try
		{
			return gaussian2D(data, width, height, sum, x0, x1, scale * zeroS0, scale * zeroS1, poissonNoise);
		}
		catch (IllegalArgumentException e)
		{
			//System.out.println(e.getMessage());
			return 0;
		}
	}

	/**
	 * Generate a scale so that at the configured zDepth the scale is 1.5.
	 * 
	 * @param z
	 * @return The scale
	 */
	private double createWidthScale(double z)
	{
		if (zDepth == 0) // Not 3D data
			return 1;

		// PSF fitting on data from the GDSC microscope show that the PSF width spread can be modelled
		// by a simple quadratic up to 1.5 times the width:
		//   width = 1 + z^2 / 2
		//         = 1.5 @ z=1

		z /= zDepth; // Scale so z=1 at the configured z-depth
		return 1.0 + z * z * 0.5;
	}

	/**
	 * Construct a Gaussian 2D function on the provided data. Only evaluates the function within +/- 5 standard
	 * deviations in each direction from the centre (allows populating large images).
	 * <p>
	 * Builds the pixel approximation using the Gaussian error function as described in Smith et al, (2010). Fast,
	 * single-molecule localisation that achieves theoretically minimum uncertainty. Nature Methods 7, 373-375
	 * (supplementary note).
	 * 
	 * @param data
	 *            The data (can be null)
	 * @param width
	 *            The data width
	 * @param height
	 *            The data height
	 * @param sum
	 *            The Gaussian integral
	 * @param x0
	 *            The Gaussian centre in dimension 0
	 * @param x1
	 *            The Gaussian centre in dimension 1
	 * @param s0
	 *            The Gaussian standard deviation dimension 0
	 * @param s1
	 *            The Gaussian standard deviation dimension 1
	 * @param poissonNoise
	 *            Add Poisson noise
	 * @return The total sum added to the image (useful when poissonNoise is added)
	 */
	public double gaussian2D(float[] data, final int width, final int height, final double sum, double x0, double x1,
			double s0, double s1, boolean poissonNoise)
	{
		if (sum == 0)
			return 0;
		// Parameter check
		if (width < 1)
			throw new IllegalArgumentException("Width cannot be less than 1");
		if (height < 1)
			throw new IllegalArgumentException("Height cannot be less than 1");
		if (data == null)
			data = new float[width * height];
		else if (data.length < width * height)
			throw new IllegalArgumentException("Data length cannot be smaller than width * height");

		s0 = Math.abs(s0);
		s1 = Math.abs(s1);

		// Evaluate the Gaussian error function on a pixel grid covering +/- 5 SD
		final int x0min = clip((int) (x0 - 5 * s0), width);
		final int x1min = clip((int) (x1 - 5 * s1), height);
		final int x0max = clip((int) Math.ceil(x0 + 5 * s0), width);
		final int x1max = clip((int) Math.ceil(x1 + 5 * s1), height);

		final int x0range = x0max - x0min;
		final int x1range = x1max - x1min;

		// min should always be less than max
		if (x0range < 1)
			throw new IllegalArgumentException("Gaussian dimension 0 range not within data bounds");
		if (x1range < 1)
			throw new IllegalArgumentException("Gaussian dimension 1 range not within data bounds");

		// Shift centre to origin and compute gaussian
		double[] gauss = gaussian2D(x0range, x1range, sum, x0 - x0min, x1 - x1min, s0, s1);

		return insert(data, x0min, x1min, x0max, x1max, width, gauss, poissonNoise);
	}

	/**
	 * Construct a Gaussian 2D function on the provided data. Only evaluates the function within +/- 5 standard
	 * deviations in each direction from the centre (allows populating large images).
	 * <p>
	 * Builds the pixel approximation using the Gaussian error function as described in Smith et al, (2010). Fast,
	 * single-molecule localisation that achieves theoretically minimum uncertainty. Nature Methods 7, 373-375
	 * (supplementary note).
	 * 
	 * @param data
	 *            The data (can be null)
	 * @param width
	 *            The data width
	 * @param height
	 *            The data height
	 * @param sum
	 *            The Gaussian integral
	 * @param x0
	 *            The Gaussian centre in dimension 0
	 * @param x1
	 *            The Gaussian centre in dimension 1
	 * @param s0
	 *            The Gaussian standard deviation dimension 0
	 * @param s1
	 *            The Gaussian standard deviation dimension 1
	 * @param poissonNoise
	 *            Add Poisson noise
	 * @return The total sum added to the image (useful when poissonNoise is added)
	 */
	public double gaussian2D(double[] data, final int width, final int height, final double sum, double x0, double x1,
			double s0, double s1, boolean poissonNoise)
	{
		if (sum == 0)
			return 0;
		// Parameter check
		if (width < 1)
			throw new IllegalArgumentException("Width cannot be less than 1");
		if (height < 1)
			throw new IllegalArgumentException("Height cannot be less than 1");
		if (data == null)
			data = new double[width * height];
		else if (data.length < width * height)
			throw new IllegalArgumentException("Data length cannot be smaller than width * height");

		s0 = Math.abs(s0);
		s1 = Math.abs(s1);

		// Evaluate the Gaussian error function on a pixel grid covering +/- 5 SD
		final int x0min = clip((int) (x0 - 5 * s0), width);
		final int x1min = clip((int) (x1 - 5 * s1), height);
		final int x0max = clip((int) Math.ceil(x0 + 5 * s0), width);
		final int x1max = clip((int) Math.ceil(x1 + 5 * s1), height);

		final int x0range = x0max - x0min;
		final int x1range = x1max - x1min;

		// min should always be less than max
		if (x0range < 1)
			throw new IllegalArgumentException("Gaussian dimension 0 range not within data bounds");
		if (x1range < 1)
			throw new IllegalArgumentException("Gaussian dimension 1 range not within data bounds");

		// Shift centre to origin and compute gaussian
		double[] gauss = gaussian2D(x0range, x1range, sum, x0 - x0min, x1 - x1min, s0, s1);

		return insert(data, x0min, x1min, x0max, x1max, width, gauss, poissonNoise);
	}

	/**
	 * Construct a Gaussian 2D function based at the origin using the specified range in each dimension.
	 * <p>
	 * Builds the pixel approximation using the Gaussian error function as described in Smith et al, (2010). Fast,
	 * single-molecule localisation that achieves theoretically minimum uncertainty. Nature Methods 7, 373-375
	 * (supplementary note).
	 * 
	 * @param x0range
	 *            The maximum range in dimension 0 (width)
	 * @param x1range
	 *            The maximum range in dimension 1 (height)
	 * @param sum
	 *            The Gaussian integral
	 * @param x0
	 *            The Gaussian centre in dimension 0
	 * @param x1
	 *            The Gaussian centre in dimension 1
	 * @param s0
	 *            The Gaussian standard deviation dimension 0
	 * @param s1
	 *            The Gaussian standard deviation dimension 1
	 * @return The data (packed in yx order, length = x0range * x1range)
	 */
	public double[] gaussian2D(int x0range, int x1range, double sum, double x0, double x1, double s0, double s1)
	{
		s0 = Math.abs(s0);
		s1 = Math.abs(s1);

		this.s0 = s0;
		this.s1 = s1;

		// Compute Gaussian error function grid up to and including the final grid position
		double[] erf0 = new double[x0range + 1];
		double[] erf1 = new double[x1range + 1];

		final double denom0 = 1.0 / (Math.sqrt(2.0) * s0);
		final double denom1 = 1.0 / (Math.sqrt(2.0) * s1);
		//final double denom0 = 1.0 / (s0);
		//final double denom1 = 1.0 / (s1);
		for (int x = 0; x <= x0range; x++)
		{
			erf0[x] = 0.5 * Erf.erf((x - x0) * denom0);
		}
		for (int y = 0; y <= x1range; y++)
		{
			erf1[y] = 0.5 * Erf.erf((y - x1) * denom1);
		}

		// Pre-compute deltaE0
		double[] deltaE0 = new double[x0range];
		for (int x = 0; x < x0range; x++)
			deltaE0[x] = erf0[x + 1] - erf0[x];

		// Compute Gaussian using the difference of the Gaussian error function
		double[] data = new double[x0range * x1range];
		for (int y = 0, i = 0; y < x1range; y++)
		{
			// Include the sum into the first deltaE to get the Gaussian integral
			final double deltaE1 = sum * (erf1[y + 1] - erf1[y]);

			for (int x = 0; x < x0range; x++, i++)
			{
				data[i] = deltaE0[x] * deltaE1;

				//// Validate using numerical integration
				//double sum2 = 0;
				//for (int ii = 0; ii < 100; ii++)
				//{
				//	double xx = x + ii / 100.0;
				//	double dx = (xx - x0) * (xx - x0) * denom0 * denom0;
				//	for (int jj = 0; jj < 100; jj++)
				//	{
				//		double yy = y + jj / 100.0;
				//		sum2 += FastMath.exp(-(dx + (yy - x1) * (yy - x1) * denom1 * denom1));
				//	}
				//}
				//sum2 *= sum / 10000 / (Math.PI * 2 * s0 * s1);
				//System.out.printf("sum=%g, sum2=%g\n", data[i], sum2);
			}
		}

		return data;
	}

	private int clip(int x, int max)
	{
		if (x < 0)
			x = 0;
		if (x > max)
			x = max;
		return x;
	}

	/**
	 * @return the Z-depth where the 3D PSF is 1.5x the width (1.5 x FWHM)
	 */
	public double getzDepth()
	{
		return zDepth;
	}

	/**
	 * @param zDepth
	 *            the Z-depth where the 3D PSF is 1.5x the width (1.5 x FWHM)
	 */
	public void setzDepth(double zDepth)
	{
		this.zDepth = Math.abs(zDepth);
	}

	/**
	 * @return The standard deviation dimension 0 for the last drawn Gaussian
	 */
	public double getS0()
	{
		return s0;
	}

	/**
	 * @return The standard deviation dimension 1 for the last drawn Gaussian
	 */
	public double getS1()
	{
		return s1;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.PSFModel#getFwhm()
	 */
	public double getFwhm()
	{
		return (s0 + s1) * Math.sqrt(2.0 * Math.log(2.0));
	}
}
