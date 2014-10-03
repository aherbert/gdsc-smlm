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
 * Contains methods for generating models of a Point Spread Function using a Airy pattern
 */
public class AiryPSFModel extends PSFModel
{
	private double zeroW0, zeroW1;
	private double w0;
	private double w1;
	private double zDepth = 0;
	private int ring = 2;

	/**
	 * The zeros of J1(x) corresponding to the rings of the Airy pattern
	 */
	public static double[] RINGS = { 0, 3.8317, 7.0156, 10.1735, 13.3237, 16.4706 };

	/**
	 * The Airy power corresponding to the rings of the Airy pattern
	 */
	public static double[] POWER;
	static
	{
		POWER = new double[RINGS.length];
		for (int i = 1; i < POWER.length; i++)
			POWER[i] = AiryPattern.power(RINGS[i]);
	}

	/**
	 * @param w0
	 *            The Airy width for dimension 0
	 * @param w1
	 *            The Airy width for dimension 1
	 */
	public AiryPSFModel(double w0, double w1)
	{
		super();
		this.zeroW0 = w0;
		this.zeroW1 = w1;
	}

	/**
	 * @param randomGenerator
	 * @param w0
	 *            The Airy width for dimension 0
	 * @param w1
	 *            The Airy width for dimension 1
	 */
	public AiryPSFModel(RandomGenerator randomGenerator, double w0, double w1)
	{
		super(randomGenerator);
		this.zeroW0 = w0;
		this.zeroW1 = w1;
	}

	/**
	 * @param randomGenerator
	 * @param w0
	 *            The Airy width for dimension 0
	 * @param w1
	 *            The Airy width for dimension 1
	 * @param zDepth
	 *            the Z-depth where the 3D PSF is 1.5x the width (1.5 x FWHM)
	 */
	public AiryPSFModel(RandomGenerator randomGenerator, double w0, double w1, double zDepth)
	{
		super(randomGenerator);
		this.zeroW0 = w0;
		this.zeroW1 = w1;
		setzDepth(zDepth);
	}

	/**
	 * @param randomDataGenerator
	 * @param w0
	 *            The Airy width for dimension 0
	 * @param w1
	 *            The Airy width for dimension 1
	 */
	public AiryPSFModel(RandomDataGenerator randomDataGenerator, double w0, double w1)
	{
		super(randomDataGenerator);
		this.zeroW0 = w0;
		this.zeroW1 = w1;
	}

	/**
	 * @param randomDataGenerator
	 * @param w0
	 *            The Airy width for dimension 0
	 * @param w1
	 *            The Airy width for dimension 1
	 * @param zDepth
	 *            the Z-depth where the 3D PSF is twice the width (2 x FWHM)
	 */
	public AiryPSFModel(RandomDataGenerator randomDataGenerator, double w0, double w1, double zDepth)
	{
		super(randomDataGenerator);
		this.zeroW0 = w0;
		this.zeroW1 = w1;
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
			final double d = airy2D(data, width, height, sum, x0, x1, scale * zeroW0, scale * zeroW1, poissonNoise);
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
			return airy2D(data, width, height, sum, x0, x1, scale * zeroW0, scale * zeroW1, poissonNoise);
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
	 * Construct a Airy pattern on the provided data. Only evaluates the function up to the configured dark ring.
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
	 * @param w0
	 *            The Airy width for dimension 0
	 * @param w1
	 *            The Airy width for dimension 1
	 * @param poissonNoise
	 *            Add Poisson noise
	 * @return The total sum added to the image (useful when poissonNoise is added)
	 */
	public double airy2D(float[] data, final int width, final int height, final double sum, double x0, double x1,
			double w0, double w1, boolean poissonNoise)
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

		w0 = Math.abs(w0);
		w1 = Math.abs(w1);

		// The second zero (dark ring of an Airy pattern is at 7.0156 of the width
		final int x0min = clip((int) (x0 - RINGS[ring] * w0), width);
		final int x1min = clip((int) (x1 - RINGS[ring] * w1), height);
		final int x0max = clip((int) Math.ceil(x0 + RINGS[ring] * w0), width);
		final int x1max = clip((int) Math.ceil(x1 + RINGS[ring] * w1), height);

		final int x0range = x0max - x0min;
		final int x1range = x1max - x1min;

		// min should always be less than max
		if (x0range < 1)
			throw new IllegalArgumentException("Dimension 0 range not within data bounds");
		if (x1range < 1)
			throw new IllegalArgumentException("Dimension 1 range not within data bounds");

		// Shift centre to origin and compute gaussian
		double[] gauss = airy2D(x0range, x1range, sum, x0 - x0min, x1 - x1min, w0, w1);

		return insert(data, x0min, x1min, x0max, x1max, width, gauss, poissonNoise);
	}

	/**
	 * Construct a Airy pattern on the provided data. Only evaluates the function up to the configured dark ring.
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
	 * @param w0
	 *            The Airy width for dimension 0
	 * @param w1
	 *            The Airy width for dimension 1
	 * @param poissonNoise
	 *            Add Poisson noise
	 * @return The total sum added to the image (useful when poissonNoise is added)
	 */
	public double airy2D(double[] data, final int width, final int height, final double sum, double x0, double x1,
			double w0, double w1, boolean poissonNoise)
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

		w0 = Math.abs(w0);
		w1 = Math.abs(w1);

		// The second zero (dark ring of an Airy pattern is at 7.0156 of the width
		final int x0min = clip((int) (x0 - RINGS[ring] * w0), width);
		final int x1min = clip((int) (x1 - RINGS[ring] * w1), height);
		final int x0max = clip((int) Math.ceil(x0 + RINGS[ring] * w0), width);
		final int x1max = clip((int) Math.ceil(x1 + RINGS[ring] * w1), height);

		final int x0range = x0max - x0min;
		final int x1range = x1max - x1min;

		// min should always be less than max
		if (x0range < 1)
			throw new IllegalArgumentException("Dimension 0 range not within data bounds");
		if (x1range < 1)
			throw new IllegalArgumentException("Dimension 1 range not within data bounds");

		// Shift centre to origin and compute gaussian
		double[] gauss = airy2D(x0range, x1range, sum, x0 - x0min, x1 - x1min, w0, w1);

		return insert(data, x0min, x1min, x0max, x1max, width, gauss, poissonNoise);
	}

	/**
	 * Construct a Airy pattern on the provided data. Only evaluates the function up to the configured dark ring.
	 * <p>
	 * Builds the pixel approximation using the cumulative power function:<br/>
	 * P(x) = 1 - (J0(x))^2- (J1(x))^2 <br/>
	 * Where JN are Bessel functions and x is the distance from the centre
	 * 
	 * @param x0range
	 *            The maximum range in dimension 0 (width)
	 * @param x1range
	 *            The maximum range in dimension 1 (height)
	 * @param sum
	 *            The integral
	 * @param x0
	 *            The centre in dimension 0
	 * @param x1
	 *            The centre in dimension 1
	 * @param w0
	 *            The Airy width for dimension 0
	 * @param w1
	 *            The Airy width for dimension 1
	 * @return The data (packed in yx order, length = x0range * x1range)
	 */
	public double[] airy2D(int x0range, int x1range, double sum, double x0, double x1, double w0, double w1)
	{
		w0 = Math.abs(w0);
		w1 = Math.abs(w1);

		this.w0 = w0;
		this.w1 = w1;

		// Limit to second dark ring
		final double limit = RINGS[ring] * RINGS[ring];
		double[] data = new double[x0range * x1range];

		// Store if the Airy pattern has been clipped
		boolean clipped = (x0 - RINGS[ring] * w0 < 0) || (x1 - RINGS[ring] * w1 < 0) ||
				(x0 + RINGS[ring] * w0 > x0range) || (x1 + RINGS[ring] * w0 > x1range);

		// Single point approximation, offset by 0.5 pixels
		x0 -= 0.5;
		x1 -= 0.5;
		int n = 0;
		double integral = 0;
		for (int y = 0, i = 0; y < x1range; y++)
		{
			double d1 = (y - x1) / w1;
			d1 *= d1;

			for (int x = 0; x < x0range; x++, i++)
			{
				double d0 = (x - x0) / w0;
				d0 *= d0;
				final double distance2 = d0 + d1;
				if (distance2 < limit)
				{
					n++;
					final double r = Math.sqrt(distance2);
					final double a = AiryPattern.intensity(r);
					data[i] = a;
					integral += a;
				}
			}
		}

		// We must normalise the integral we calculated to the correct power of the Airy pattern.
		//System.out.printf("Norm = %g (%g) = %gx\n", n / 12.0, integral / POWER[2], (n / 12.0) / (integral / POWER[2]));
		if (clipped)
		{
			// Use the normalising approximation for square pixels to calculate the integral.
			// TODO - check where this square pixel approximation is valid. It is used in the Mortensen
			// formula for precision.
			sum *= 1 / (n / 12.0);
		}
		else
		{
			// Note: This will be incorrect if the Airy pattern was clipped at the edge (i.e. does not represent
			// the full area of the calculated Airy rings)
			sum *= POWER[ring] / integral;
		}
		for (int i = 0; i < data.length; i++)
			data[i] *= sum;

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
	 * @return The width in dimension 0 for the last drawn Airy pattern
	 */
	public double getW0()
	{
		return w0;
	}

	/**
	 * @return The width in dimension 0 for the last drawn Airy pattern
	 */
	public double getW1()
	{
		return w1;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.PSFModel#getFwhm()
	 */
	public double getFwhm()
	{
		// Use the Gaussian approximation to the Airy pattern to get the FWHM
		// Gaussian SD = 1.323 * Airy width
		// Gaussian FWHM = SD * 2.35 = SD * (2 * Math.sqrt(2.0 * Math.log(2.0)))
		return AiryPattern.FACTOR * (w0 + w1) * Math.sqrt(2.0 * Math.log(2.0));
	}

	/**
	 * @return the ring limit for the calculated Airy pattern
	 */
	public int getRing()
	{
		return ring;
	}

	/**
	 * Set the limit of the Airy pattern, defined by the dark rings where the pattern is zero. Allowed values are 1-5.
	 * 
	 * @param ring
	 *            the ring limit for the calculated Airy pattern
	 */
	public void setRing(int ring)
	{
		if (ring < RINGS.length && ring > 1)
			this.ring = ring;
	}
}
