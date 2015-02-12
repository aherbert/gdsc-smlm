package gdsc.smlm.filters;

import gdsc.smlm.utils.StoredDataStatistics;

import java.awt.Rectangle;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2015 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Computes the mean using a circular mask.
 * <p>
 * Adapted from ij.plugin.filter.RankFilters
 */
public class CircularMeanFilter implements Cloneable
{
	private int[] kernel = null;
	private double lastRadius = 0;

	/**
	 * Compute the mean.
	 * Pixels within border regions (defined by 1 x radius) are unchanged.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param radius
	 *            The circle radius
	 */
	public void convolveInternal(float[] data, final int maxx, final int maxy, final double radius)
	{
		final int border = getBorder(radius);
		Rectangle roi = new Rectangle(border, border, maxx - 2 * border, maxy - 2 * border);
		if (roi.width < 1 || roi.height < 1)
			return;
		rank(data, roi, maxx, maxy, radius);
	}

	/**
	 * Get the border that will be ignored for the specified radius
	 * 
	 * @param radius
	 *            the radius
	 * @return The border
	 */
	public static int getBorder(double radius)
	{
		return (int) Math.ceil(radius);
	}

	/**
	 * Compute the mean.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param radius
	 *            The circle radius
	 */
	public void convolve(float[] data, final int maxx, final int maxy, final double radius)
	{
		Rectangle roi = new Rectangle(maxx, maxy);
		rank(data, roi, maxx, maxy, radius);
	}

	/**
	 * 
	 * @param radius
	 *            The kernel radius
	 */
	private void rank(float[] data, Rectangle roi, int width, int height, double radius)
	{
		int[] lineRadii = makeLineRadii(radius);

		int kHeight = kHeight(lineRadii);
		int kRadius = kRadius(lineRadii);
		final int cacheWidth = roi.width + 2 * kRadius;
		final int cacheHeight = kHeight;
		// 'cache' is the input buffer. Each line y in the image is mapped onto cache line y%cacheHeight
		final float[] cache = new float[cacheWidth * cacheHeight];

		doFiltering(data, roi, width, height, lineRadii, cache, cacheWidth, cacheHeight);
	}

	// Filter a grayscale image or one channel of an RGB image using one thread
	//
	// Data handling: The area needed for processing a line is written into the array 'cache'.
	// This is a stripe of sufficient width for all threads to have each thread processing one
	// line, and some extra space if one thread is finished to start the next line.
	// This array is padded at the edges of the image so that a surrounding with radius kRadius
	// for each pixel processed is within 'cache'. Out-of-image
	// pixels are set to the value of the nearest edge pixel. When adding a new line, the lines in
	// 'cache' are not shifted but rather the smaller array with the start and end pointers of the
	// kernel area is modified to point at the addresses for the next line.
	//
	// Algorithm: For mean and variance, except for very small radius, usually do not calculate the
	// sum over all pixels. This sum is calculated for the first pixel of every line only. For the
	// following pixels, add the new values and subtract those that are not in the sum any more.
	private void doFiltering(float[] pixels, Rectangle roi, int width, int height, int[] lineRadii, float[] cache,
			int cacheWidth, int cacheHeight)
	{
		int kHeight = kHeight(lineRadii);
		int kRadius = kRadius(lineRadii);
		int kNPoints = kNPoints(lineRadii);

		int xmin = roi.x - kRadius;
		int xmax = roi.x + roi.width + kRadius;
		int[] cachePointers = makeCachePointers(lineRadii, cacheWidth);

		int padLeft = xmin < 0 ? -xmin : 0;
		int padRight = xmax > width ? xmax - width : 0;
		int xminInside = xmin > 0 ? xmin : 0;
		int xmaxInside = xmax < width ? xmax : width;
		int widthInside = xmaxInside - xminInside;

		double[] sums = new double[2];

		boolean smallKernel = kRadius < 2;

		float maxValue = Float.NaN;
		float[] values = pixels;

		int previousY = kHeight / 2 - cacheHeight;

		for (int y = roi.y; y < roi.y + roi.height; y++)
		{
			for (int i = 0; i < cachePointers.length; i++)
				//shift kernel pointers to new line
				cachePointers[i] = (cachePointers[i] + cacheWidth * (y - previousY)) % cache.length;
			previousY = y;

			int yStartReading = y == roi.y ? Math.max(roi.y - kHeight / 2, 0) : y + kHeight / 2;
			for (int yNew = yStartReading; yNew <= y + kHeight / 2; yNew++)
			{ //only 1 line except at start
				readLineToCacheOrPad(pixels, width, height, roi.y, xminInside, widthInside, cache, cacheWidth,
						cacheHeight, padLeft, padRight, kHeight, yNew);
			}

			int cacheLineP = cacheWidth * (y % cacheHeight) + kRadius; //points to pixel (roi.x, y)
			filterLine(values, width, cache, cachePointers, kNPoints, cacheLineP, roi, y, // F I L T E R
					sums, maxValue, smallKernel);
		}
	}

	private void filterLine(float[] values, int width, float[] cache, int[] cachePointers, int kNPoints,
			int cacheLineP, Rectangle roi, int y, double[] sums, float maxValue, boolean smallKernel)
	{
		int valuesP = roi.x + y * width;

		// NOTE:
		// The incremental algorithm does not work.
		// The full calculation is always true in the original source code.
		boolean fullCalculation = true;// smallKernel; //for small kernel, always use the full area, not incremental algorithm

		for (int x = 0; x < roi.width; x++, valuesP++)
		{ // x is with respect to roi.x
			if (fullCalculation)
			{
				getAreaSums(cache, x, cachePointers, sums);
			}
			else
			{
				addSideSums(cache, x, cachePointers, sums);
				if (Double.isNaN(sums[0])) //avoid perpetuating NaNs into remaining line
					fullCalculation = true;
			}
			values[valuesP] = (float) (sums[0] / kNPoints);
		} // for x
	}

	/**
	 * Read a line into the cache (including padding in x).
	 * If y>=height, instead of reading new data, it duplicates the line y=height-1.
	 * If y==0, it also creates the data for y<0, as far as necessary, thus filling the cache with
	 * more than one line (padding by duplicating the y=0 row).
	 */
	private static void readLineToCacheOrPad(Object pixels, int width, int height, int roiY, int xminInside,
			int widthInside, float[] cache, int cacheWidth, int cacheHeight, int padLeft, int padRight, int kHeight,
			int y)
	{
		int lineInCache = y % cacheHeight;
		if (y < height)
		{
			readLineToCache(pixels, y * width, xminInside, widthInside, cache, lineInCache * cacheWidth, padLeft,
					padRight);
			if (y == 0)
				for (int prevY = roiY - kHeight / 2; prevY < 0; prevY++)
				{ //for y<0, pad with y=0 border pixels 
					int prevLineInCache = cacheHeight + prevY;
					System.arraycopy(cache, 0, cache, prevLineInCache * cacheWidth, cacheWidth);
				}
		}
		else
			System.arraycopy(cache, cacheWidth * ((height - 1) % cacheHeight), cache, lineInCache * cacheWidth,
					cacheWidth);
	}

	/** Read a line into the cache (includes conversion to flaot). Pad with edge pixels in x if necessary */
	private static void readLineToCache(Object pixels, int pixelLineP, int xminInside, int widthInside, float[] cache,
			int cacheLineP, int padLeft, int padRight)
	{
		System.arraycopy(pixels, pixelLineP + xminInside, cache, cacheLineP + padLeft, widthInside);

		for (int cp = cacheLineP; cp < cacheLineP + padLeft; cp++)
			cache[cp] = cache[cacheLineP + padLeft];
		for (int cp = cacheLineP + padLeft + widthInside; cp < cacheLineP + padLeft + widthInside + padRight; cp++)
			cache[cp] = cache[cacheLineP + padLeft + widthInside - 1];
	}

	/**
	 * Get sum of values and values squared within the kernel area.
	 * x between 0 and cacheWidth-1
	 * Output is written to array sums[0] = sum
	 */
	private static void getAreaSums(float[] cache, int xCache0, int[] kernel, double[] sums)
	{
		double sum = 0;
		for (int kk = 0; kk < kernel.length; kk++)
		{ // y within the cache stripe (we have 2 kernel pointers per cache line)
			for (int p = kernel[kk++] + xCache0; p <= kernel[kk] + xCache0; p++)
			{
				float v = cache[p];
				sum += v;
			}
		}
		sums[0] = sum;
		return;
	}

	/**
	 * Add all values and values squared at the right border inside minus at the left border outside the kernel area.
	 * Output is added or subtracted to/from array sums[0] += sum when at
	 * the right border, minus when at the left border
	 */
	private static void addSideSums(float[] cache, int xCache0, int[] kernel, double[] sums)
	{
		double sum = 0;
		for (int kk = 0; kk < kernel.length; /* k++;k++ below */)
		{
			float v = cache[kernel[kk++] + (xCache0 - 1)];
			sum -= v;
			v = cache[kernel[kk++] + xCache0];
			sum += v;
		}
		sums[0] += sum;
		return;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public Object clone()
	{
		try
		{
			CircularMeanFilter o = (CircularMeanFilter) super.clone();
			return o;
		}
		catch (CloneNotSupportedException e)
		{
			// Ignore
		}
		return null;
	}

	/**
	 * Create a circular kernel (structuring element) of a given radius.
	 *
	 * @param radius
	 *            Radius = 0.5 includes the 4 neighbors of the pixel in the center,
	 *            radius = 1 corresponds to a 3x3 kernel size.
	 * @return:
	 *          The output is an array that gives the length of each line of the structuring element
	 *          (kernel) to the left (negative) and to the right (positive):
	 *          [0] left in line 0, [1] right in line 0,
	 *          [2] left in line 2, ...
	 *          The maximum (absolute) value should be kernelRadius.
	 *          Array elements at the end:
	 *          length-2: nPoints, number of pixels in the kernel area
	 *          length-1: kernelRadius in x direction (kernel width is 2*kernelRadius+1)
	 *          Kernel height can be calculated as (array length - 1)/2 (odd number);
	 *          Kernel radius in y direction is kernel height/2 (truncating integer division).
	 *          Note that kernel width and height are the same for the circular kernels used here,
	 *          but treated separately for the case of future extensions with non-circular kernels.
	 */
	private int[] makeLineRadii(double radius)
	{
		// Cache the kernel
		if (kernel != null && lastRadius == radius)
			return kernel;
		lastRadius = radius;

		//System.out.println(java.util.Arrays.toString(getRadii(10, 0.1)));

		if (radius >= 1.5 && radius < 1.75) //this code creates the same sizes as the previous RankFilters
			radius = 1.75;
		else if (radius >= 2.5 && radius < 2.85)
			radius = 2.85;
		int r2 = (int) (radius * radius) + 1;
		int kRadius = (int) (Math.sqrt(r2 + 1e-10));
		int kHeight = 2 * kRadius + 1;
		kernel = new int[2 * kHeight + 2];
		kernel[2 * kRadius] = -kRadius;
		kernel[2 * kRadius + 1] = kRadius;
		int nPoints = 2 * kRadius + 1;
		for (int y = 1; y <= kRadius; y++)
		{ //lines above and below center together
			int dx = (int) (Math.sqrt(r2 - y * y + 1e-10));
			kernel[2 * (kRadius - y)] = -dx;
			kernel[2 * (kRadius - y) + 1] = dx;
			kernel[2 * (kRadius + y)] = -dx;
			kernel[2 * (kRadius + y) + 1] = dx;
			nPoints += 4 * dx + 2; //2*dx+1 for each line, above&below
		}
		kernel[kernel.length - 2] = nPoints;
		kernel[kernel.length - 1] = kRadius;
		//for (int i=0; i<kHeight;i++)IJ.log(i+": "+kernel[2*i]+"-"+kernel[2*i+1]);
		return kernel;
	}

	/**
	 * Get the radii where different smoothing will be applied
	 * 
	 * @param max
	 *            The maximum radii to include
	 * @param increment
	 *            The increment to use between radii
	 * @return The radiis where a different smoothing will be applied
	 */
	public static double[] getRadii(double max, double increment)
	{
		StoredDataStatistics radii = new StoredDataStatistics();
		double lastN = 0;
		if (increment < 0)
			increment = 0.1;
		for (double r = 0.5; r <= max; r += increment)
		{
			int nPoints = getNPoints(r);
			if (nPoints > lastN)
				radii.add(r);
			lastN = nPoints;
		}
		return radii.getValues();
	}

	/**
	 * Count the number of points in the circle mask for the given radius
	 * 
	 * @param radius
	 * @return The number of points
	 */
	public static int getNPoints(double radius)
	{
		if (radius >= 1.5 && radius < 1.75) //this code creates the same sizes as the previous RankFilters
			radius = 1.75;
		else if (radius >= 2.5 && radius < 2.85)
			radius = 2.85;
		int r2 = (int) (radius * radius) + 1;
		int kRadius = (int) (Math.sqrt(r2 + 1e-10));
		int nPoints = 2 * kRadius + 1;
		for (int y = 1; y <= kRadius; y++)
		{
			int dx = (int) (Math.sqrt(r2 - y * y + 1e-10));
			nPoints += 4 * dx + 2;
		}
		return nPoints;
	}

	/**
	 * Get the diameter of the pixel region used
	 * 
	 * @param radius
	 * @return The diameter
	 */
	public static int getDiameter(double radius)
	{
		int r2 = (int) (radius * radius) + 1;
		int kRadius = (int) (Math.sqrt(r2 + 1e-10));
		int kHeight = 2 * kRadius + 1;
		return kHeight;
	}

	/**
	 * Count the number of points in the circle mask for the given radius. Then convert it into an approximate radius
	 * using sqrt(nPoints/pi)
	 * 
	 * @param radius
	 * @return The diameter
	 */
	public static double getPixelRadius(double radius)
	{
		return Math.sqrt(getNPoints(radius) / Math.PI);
	}

	//kernel height
	private int kHeight(int[] lineRadii)
	{
		return (lineRadii.length - 2) / 2;
	}

	//kernel radius in x direction. width is 2+kRadius+1
	private int kRadius(int[] lineRadii)
	{
		return lineRadii[lineRadii.length - 1];
	}

	//number of points in kernal area
	private int kNPoints(int[] lineRadii)
	{
		return lineRadii[lineRadii.length - 2];
	}

	//cache pointers for a given kernel
	private int[] makeCachePointers(int[] lineRadii, int cacheWidth)
	{
		int kRadius = kRadius(lineRadii);
		int kHeight = kHeight(lineRadii);
		int[] cachePointers = new int[2 * kHeight];
		for (int i = 0; i < kHeight; i++)
		{
			cachePointers[2 * i] = i * cacheWidth + kRadius + lineRadii[2 * i];
			cachePointers[2 * i + 1] = i * cacheWidth + kRadius + lineRadii[2 * i + 1];
		}
		return cachePointers;
	}
}