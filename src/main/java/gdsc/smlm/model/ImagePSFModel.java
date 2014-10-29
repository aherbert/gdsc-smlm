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
import org.apache.commons.math3.util.FastMath;

/**
 * Generates a Point Spread Function using an image constructed from diffraction limited spots imaged axially through
 * the plane of focus.
 * <p>
 * The input image must be square. The X/Y centre is the middle of the square. The stack can have any number of slices.
 * The z-centre must be identified. Any pixels below zero will be set to zero.
 * <p>
 * The model can be used to draw a PSF for a point on an image. If the z coordinate is positive then the PSF image from
 * a negative index (below the z-centre) is used. If the z-coordinate is negative then the PSF image from a positive
 * index (above the z-centre) is used. I.e. the input stack is assumed to be imaged axially with increasing z-stage
 * position moving the stage closer to the objective.
 */
public class ImagePSFModel extends PSFModel
{
	//private double[][] inputImage;
	private double[][] sumImage;
	private int psfWidth, xyCentre, zCentre;
	private double halfPsfWidthInPixels;
	private double unitsPerPixel;
	private double unitsPerSlice;
	private double fwhm;

	/**
	 * @param image
	 *            The image consisting of a stack of square pixel buffers. The buffers are stored in YX order.
	 * @param zCentre
	 *            The centre of the PSF image
	 * @param unitsPerPixel
	 *            The distance between adjacent X/Y pixels
	 * @param unitsPerSlice
	 *            The distance between adjacent Z pixels
	 * @param fwhm
	 *            The full-width at half-maximum for the z-centre
	 */
	public ImagePSFModel(float[][] image, int zCentre, double unitsPerPixel, double unitsPerSlice, double fwhm)
	{
		super();
		init(image, zCentre, unitsPerPixel, unitsPerSlice, fwhm);
	}

	/**
	 * Private constructor used in the {@link #copy()} method
	 */
	private ImagePSFModel()
	{
		super();
	}

	private void init(float[][] image, int zCentre, double unitsPerPixel, double unitsPerSlice, double fwhm)
	{
		if (image == null || image.length == 0)
			throw new IllegalArgumentException("Image cannot be null/empty");
		for (int i = 0; i < image.length; i++)
			if (image[i] == null)
				throw new IllegalArgumentException("Image contains null plane");
		if (zCentre < 0 || zCentre >= image.length)
			throw new IllegalArgumentException("z-centre is not within the bounds of the image stack");
		final int size = image[0].length;
		double edge = Math.sqrt(size);
		if (edge != (int) edge)
			throw new IllegalArgumentException("Image planes are not square");
		psfWidth = (int) edge;
		//if (psfWidth % 2 != 1)
		//	throw new IllegalArgumentException("Image edge length is not an odd number");
		xyCentre = psfWidth / 2;
		for (int i = 1; i < image.length; i++)
			if (image[i].length != size)
				throw new IllegalArgumentException("Image planes are not the same size");
		this.zCentre = zCentre;
		if (unitsPerPixel <= 0 || unitsPerPixel > 1)
			throw new IllegalArgumentException("Units per pixel must be between 0 and 1");
		if (unitsPerSlice <= 0 || unitsPerSlice > 1)
			throw new IllegalArgumentException("Units per slice must be between 0 and 1");
		this.unitsPerPixel = unitsPerPixel;
		this.unitsPerSlice = unitsPerSlice;

		halfPsfWidthInPixels = 0.5 * psfWidth * unitsPerPixel;

		// Create the image.
		this.sumImage = duplicate(image);
		// Normalise so that the highest intensity frame sums to 1.
		normalise(this.sumImage);

		// Used for debugging
		//inputImage = new double[sumImage.length][];
		//for (int i = 0; i < sumImage.length; i++)
		//	inputImage[i] = Arrays.copyOf(sumImage[i], sumImage[i].length);

		// Then create a rolling sum table.
		for (int i = 0; i < sumImage.length; i++)
			calculateRollingSums(sumImage[i]);
	}

	private double[][] duplicate(float[][] image)
	{
		final int size = image[0].length;
		double[][] duplicate = new double[image.length][size];
		for (int i = 0; i < image.length; i++)
			for (int j = 0; j < size; j++)
				duplicate[i][j] = image[i][j];
		return duplicate;
	}

	/**
	 * Normalise the image so that the z-centre frame has a sum of 1. Ensure no pixels are below zero.
	 * 
	 * @param image
	 */
	private void normalise(double[][] image)
	{
		if (image == null || image.length == 0)
			return;

		// Reset negative pixels
		double[] data = image[zCentre];
		for (int j = 0; j < data.length; j++)
			if (data[j] < 0)
				data[j] = 0;

		double max = sum(data);
		if (max <= 0)
			return;

		for (int i = 0; i < image.length; i++)
		{
			data = image[i];
			for (int j = 0; j < data.length; j++)
			{
				if (data[j] < 0)
					data[j] = 0;
				else
					data[j] /= max;
			}
		}
	}

	private double sum(double[] data)
	{
		double sum = 0;
		for (double f : data)
			sum += f;
		return sum;
	}

	private void calculateRollingSums(double[] s)
	{
		// Compute the rolling sum and sum of squares
		// s(u,v) = f(u,v) + s(u-1,v) + s(u,v-1) - s(u-1,v-1) 
		// where s(u,v) = 0 when either u,v < 0

		final int maxx = psfWidth;
		final int maxy = psfWidth;

		// First row
		double cs = 0; // Column sum
		for (int i = 0; i < maxx; i++)
		{
			cs += s[i];
			s[i] = cs;
		}

		// Remaining rows:
		// sum = rolling sum of row + sum of row above
		for (int y = 1, i = maxx; y < maxy; y++)
		{
			cs = 0;

			// Remaining columns
			for (int x = 0; x < maxx; x++, i++)
			{
				cs += s[i];
				s[i] = (s[i - maxx] + cs);
			}
		}
	}

	/**
	 * @param randomGenerator
	 * @param image
	 *            The image consisting of a stack of square pixel buffers. The buffers are stored in YX order.
	 * @param zCentre
	 *            The centre of the PSF image
	 * @param unitsPerPixel
	 *            The distance between adjacent X/Y pixels
	 * @param unitsPerSlice
	 *            The distance between adjacent Z pixels
	 * @param fwhm
	 *            The full-width at half-maximum for the z-centre
	 */
	public ImagePSFModel(RandomGenerator randomGenerator, float[][] image, int zCentre, double unitsPerPixel,
			double unitsPerSlice, double fwhm)
	{
		super(randomGenerator);
		init(image, zCentre, unitsPerPixel, unitsPerSlice, fwhm);
	}

	/**
	 * @param randomDataGenerator
	 * @param image
	 *            The image consisting of a stack of square pixel buffers. The buffers are stored in YX order.
	 * @param zCentre
	 *            The centre of the PSF image
	 * @param unitsPerPixel
	 *            The distance between adjacent X/Y pixels
	 * @param unitsPerSlice
	 *            The distance between adjacent Z pixels
	 * @param fwhm
	 *            The full-width at half-maximum for the z-centre
	 */
	public ImagePSFModel(RandomDataGenerator randomDataGenerator, float[][] image, int zCentre, double unitsPerPixel,
			double unitsPerSlice, double fwhm)
	{
		super(randomDataGenerator);
		init(image, zCentre, unitsPerPixel, unitsPerSlice, fwhm);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.PSFModel#create3D(float[], int, int, double, double, double, double, boolean)
	 */
	public double create3D(float[] data, final int width, final int height, final double sum, double x0, double x1,
			double x2, boolean poissonNoise)
	{
		try
		{
			return drawPSF(data, width, height, sum, x0, x1, x2, poissonNoise);
		}
		catch (IllegalArgumentException e)
		{
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
		try
		{
			return drawPSF(data, width, height, sum, x0, x1, x2, poissonNoise);
		}
		catch (IllegalArgumentException e)
		{
			return 0;
		}
	}

	public double drawPSF(float[] data, final int width, final int height, final double sum, double x0, double x1,
			double x2, boolean poissonNoise)
	{
		// Parameter check
		if (width < 1)
			throw new IllegalArgumentException("Width cannot be less than 1");
		if (height < 1)
			throw new IllegalArgumentException("Height cannot be less than 1");
		if (data == null)
			data = new float[width * height];
		else if (data.length < width * height)
			throw new IllegalArgumentException("Data length cannot be smaller than width * height");

		// Evaluate the PSF over the full range
		final int x0min = clip((int) (x0 - halfPsfWidthInPixels), width);
		final int x1min = clip((int) (x1 - halfPsfWidthInPixels), height);
		final int x0max = clip((int) Math.ceil(x0 + halfPsfWidthInPixels), width);
		final int x1max = clip((int) Math.ceil(x1 + halfPsfWidthInPixels), height);

		final int x0range = x0max - x0min;
		final int x1range = x1max - x1min;

		// min should always be less than max
		if (x0range < 1)
			throw new IllegalArgumentException("Dimension 0 range not within data bounds");
		if (x1range < 1)
			throw new IllegalArgumentException("Dimension 1 range not within data bounds");

		// Shift centre to origin and draw the PSF
		double[] psf = drawPSF(x0range, x1range, sum, x0 - x0min, x1 - x1min, x2);

		return insert(data, x0min, x1min, x0max, x1max, width, psf, poissonNoise);
	}

	public double drawPSF(double[] data, final int width, final int height, final double sum, double x0, double x1,
			double x2, boolean poissonNoise)
	{
		// Parameter check
		if (width < 1)
			throw new IllegalArgumentException("Width cannot be less than 1");
		if (height < 1)
			throw new IllegalArgumentException("Height cannot be less than 1");
		if (data == null)
			data = new double[width * height];
		else if (data.length < width * height)
			throw new IllegalArgumentException("Data length cannot be smaller than width * height");

		// Evaluate the PSF over the full range
		final int x0min = clip((int) (x0 - halfPsfWidthInPixels), width);
		final int x1min = clip((int) (x1 - halfPsfWidthInPixels), height);
		final int x0max = clip((int) Math.ceil(x0 + halfPsfWidthInPixels), width);
		final int x1max = clip((int) Math.ceil(x1 + halfPsfWidthInPixels), height);

		final int x0range = x0max - x0min;
		final int x1range = x1max - x1min;

		// min should always be less than max
		if (x0range < 1)
			throw new IllegalArgumentException("Dimension 0 range not within data bounds");
		if (x1range < 1)
			throw new IllegalArgumentException("Dimension 1 range not within data bounds");

		// Shift centre to origin and draw the PSF
		double[] psf = drawPSF(x0range, x1range, sum, x0 - x0min, x1 - x1min, x2);

		return insert(data, x0min, x1min, x0max, x1max, width, psf, poissonNoise);
	}

	/**
	 * Construct a PSF function based at the origin using the specified range in each dimension.
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
	 * @param x2
	 *            The centre in dimension 2
	 * @return The data (packed in yx order, length = x0range * x1range)
	 */
	public double[] drawPSF(int x0range, int x1range, double sum, double x0, double x1, double x2)
	{
		double[] data = new double[x0range * x1range];

		// Determine the slice of the PSF.
		// We assume the PSF was imaged axially with increasing z-stage position (moving the stage 
		// closer to the objective). Thus we invert the z-coordinate to find the appropriate slice.
		final int slice = (int) Math.round(-x2 / unitsPerSlice) + zCentre;
		if (slice < 0 || slice >= sumImage.length)
			return data;
		final double[] sumPsf = sumImage[slice];

		// Used for debugging
		//final double[] psf = inputImage[slice];
		//int ok=0, error=0;

		// Determine the insert offset
		final int xOffset = (int) Math.round(x0 / unitsPerPixel) - xyCentre;
		final int yOffset = (int) Math.round(x1 / unitsPerPixel) - xyCentre;

		// Determine PSF blocks
		int[] u = createLookup(x0range, xOffset);
		int[] v = createLookup(x1range, yOffset);

		// Data will be the sum of the input pixels from (u,v) to (u+1,v+1)

		// Check data can be inserted from the PSF
		if (u[0] > psfWidth - 1 || v[0] > psfWidth - 1 || u[x0range] < 0 || v[x1range] < 0)
			return data;

		for (int y = 0; y < x1range; y++)
		{
			if (v[y] > psfWidth - 1)
				break;
			if (v[y + 1] < 0)
				continue;
			final int lowerV = v[y];
			final int upperV = FastMath.min(v[y + 1], psfWidth - 1);
			for (int x = 0, i = y * x0range; x < x0range; x++, i++)
			{
				if (u[x] > psfWidth - 1)
					break;
				if (u[x + 1] < 0)
					continue;
				data[i] = sum(sumPsf, u[x], lowerV, u[x + 1], upperV);

				//// Check using the original input image that the sum is correct
				//double s = 0;
				//for (int vv = lowerV + 1; vv <= upperV; vv++)
				//{
				//	if (vv > psfWidth-1)
				//		break;
				//	if (vv < 0)
				//		continue;
				//	for (int uu = u[x] + 1, index = vv * psfWidth + u[x] + 1; uu <= u[x + 1]; uu++, index++)
				//	{
				//		if (uu > psfWidth-1)
				//			break;
				//		if (uu < 0)
				//			continue;
				//		s += psf[index];
				//	}
				//}
				//if (s / data[i] > 1.1 || s / data[i] < 0.9)
				//{
				//	error++;
				//	System.out.printf("Mistake %f != %f\n", s, data[i]);
				//}
				//else
				//	ok++;
			}
		}

		// The PSF is normalised so the brightest plane is 1. Just multiply by the sum to create the integral.
		for (int i = 0; i < data.length; i++)
			data[i] *= sum;

		//System.out.printf("%% OK = %f (Evaluated %f) : Total = %f (%f)\n", 
		//		(100.0*ok) / (ok + error), (100.* (ok + error)) / (x0range*x1range), sum(data), sum * sum(psf));

		return data;
	}

	private int[] createLookup(int range, int offset)
	{
		int[] pixel = new int[range + 1];
		for (int i = 0; i < pixel.length; i++)
		{
			pixel[i] = (int) Math.round(i / unitsPerPixel) - offset - 1;
		}
		return pixel;
	}

	private double sum(double[] s, int lowerU, int lowerV, int upperU, int upperV)
	{
		// Compute sum from rolling sum using:
		// sum = 
		// + s(upperU,upperV) 
		// - s(lowerU,upperV)
		// - s(upperU,lowerV)
		// + s(lowerU,lowerV)
		// Note: 
		// s(u,v) = 0 when either u,v < 0
		// s(u,v) = s(umax,v) when u>umax
		// s(u,v) = s(u,vmax) when v>vmax
		// s(u,v) = s(umax,vmax) when u>umax,v>vmax

		//		if (lowerU > psfWidth - 1 || lowerV > psfWidth - 1)
		//			return 0;
		//		if (upperU < 0 || upperV < 0)
		//			return 0;

		upperU = FastMath.min(upperU, psfWidth - 1);
		//upperV = FastMath.min(upperV, psfWidth - 1);

		int index = upperV * psfWidth + upperU;
		double sum = s[index];

		if (lowerU >= 0)
		{
			index = upperV * psfWidth + lowerU;
			sum -= s[index];
		}
		if (lowerV >= 0)
		{
			index = lowerV * psfWidth + upperU;
			sum -= s[index];

			if (lowerU >= 0)
			{
				// + s(u-1,v-1)
				index = lowerV * psfWidth + lowerU;
				sum += s[index];
			}
		}
		return sum;
	}

	private int clip(int x, int max)
	{
		if (x < 0)
			x = 0;
		if (x > max)
			x = max;
		return x;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.PSFModel#getFwhm()
	 */
	public double getFwhm()
	{
		return fwhm;
	}

	/**
	 * Produce a shallow copy of this object. This shares the pre-computed PSF image data but will allow
	 * the copy to store its own version of the most recently created PSF. The copy has a new random data generator.
	 * 
	 * @return A shallow copy of this object
	 */
	public ImagePSFModel copy()
	{
		ImagePSFModel model = new ImagePSFModel();
		model.sumImage = sumImage;
		//model.inputImage = inputImage;
		model.psfWidth = psfWidth;
		model.xyCentre = xyCentre;
		model.zCentre = zCentre;
		model.halfPsfWidthInPixels = halfPsfWidthInPixels;
		model.unitsPerPixel = unitsPerPixel;
		model.unitsPerSlice = unitsPerSlice;
		model.fwhm = fwhm;
		return model;
	}
}
