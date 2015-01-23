package gdsc.smlm.filters;

import gdsc.smlm.ij.utils.Utils;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

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
 * Identifies candidate spots (local maxima) in an image. The image is smoothed with an average box filter. A simple
 * weighted approximation is used for less than 1 pixel smoothing.
 */
public class AverageSpotFilter extends MaximaSpotFilter
{
	private final double smooth;
	private AverageFilter filter = new AverageFilter();

	/**
	 * Constructor
	 * 
	 * @param search
	 *            The search width for non-maximum suppression
	 * @param border
	 *            The border to ignore for maxima
	 * @param smooth
	 *            The smoothing width to apply to the data
	 * @throws IllegalArgumentException
	 *             if smooth is below zero
	 */
	public AverageSpotFilter(int search, int border, double smooth)
	{
		super(search, border);
		if (smooth < 0)
			throw new IllegalArgumentException("The smoothing width must be positive");
		this.smooth = convert(smooth);
	}

	/**
	 * Convert the smoothing parameter to the value which is used for the AverageFilter.
	 * We only use int smoothing above 1, all values below 1 use the input value. Values below zero are set to zero.
	 * 
	 * @see AverageFilter
	 * 
	 * @param smooth
	 * @return The adjusted value
	 */
	public static double convert(double smooth)
	{
		if (smooth < 0)
			return 0;
		if (smooth < 1)
			return smooth;
		return (int) smooth;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.SpotFilter#isAbsoluteIntensity()
	 */
	public boolean isAbsoluteIntensity()
	{
		return true;
	}

	@Override
	public float[] preprocessData(float[] data, int width, int height)
	{
		if (smooth > 0)
			data = applySmoothing(data, width, height, smooth);
		return data;
	}

	/**
	 * Apply a box smoothing filter to the data
	 * 
	 * @param data
	 * @param width
	 * @param height
	 * @param smooth
	 * @return
	 */
	protected float[] applySmoothing(float[] data, int width, int height, double smooth)
	{
		// Smoothing destructively modifies the data so create a copy
		float[] smoothData = Arrays.copyOf(data, width * height);

		if (smooth < 1)
		{
			// Use a weighted smoothing for less than 1 pixel box size
			if (smooth <= getBorder())
				filter.blockAverage3x3Internal(smoothData, width, height, (float) smooth);
			else
				filter.blockAverage3x3(smoothData, width, height, (float) smooth);
		}
		else
		{
			// Check upper limits are safe
			int tmpSmooth = FastMath.min((int) smooth, FastMath.min(width, height) / 2);

			if (tmpSmooth <= getBorder())
			{
				filter.rollingBlockAverageInternal(smoothData, width, height, tmpSmooth);
			}
			else
			{
				filter.rollingBlockAverage(smoothData, width, height, tmpSmooth);
			}
		}
		return smoothData;
	}

	/**
	 * @return the smoothing width
	 */
	public double getSmooth()
	{
		return smooth;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public Object clone()
	{
		AverageSpotFilter f = (AverageSpotFilter) super.clone();
		// Ensure the object is duplicated and not passed by reference.
		f.filter = (AverageFilter) filter.clone();
		return f;
	}

	@Override
	public String getName()
	{
		return "Average Filter";
	}

	@Override
	public List<String> getParameters()
	{
		List<String> list = super.getParameters();
		list.add("smooth = " + Utils.rounded(smooth));
		return list;
	}
}