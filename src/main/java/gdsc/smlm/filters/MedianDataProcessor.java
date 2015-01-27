package gdsc.smlm.filters;

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
 * Identifies candidate spots (local maxima) in an image. The image is smoothed with an median box filter.
 */
public class MedianDataProcessor extends DataProcessor
{
	private final int smooth;
	private MedianFilter filter = new MedianFilter();

	/**
	 * Constructor
	 * 
	 * @param border
	 *            The border to ignore for maxima
	 * @param smooth
	 *            The smoothing width to apply to the data
	 * @throws IllegalArgumentException
	 *             if smooth is below zero
	 */
	public MedianDataProcessor(int border, double smooth)
	{
		super(border);
		this.smooth = convert(smooth);
	}

	/**
	 * Convert the smoothing parameter to the value which is used for the MedianFilter.
	 * We only use int smoothing. Values below zero are set to zero.
	 * 
	 * @see MedianFilter
	 * 
	 * @param smooth
	 * @return The adjusted value
	 */
	public static int convert(double smooth)
	{
		if (smooth < 0)
			return 0;
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

	/**
	 * @param data
	 * @param width
	 * @param height
	 * @return
	 */
	@Override
	public float[] process(float[] data, int width, int height)
	{
		float[] smoothData = data;
		if (smooth > 0)
		{
			// Smoothing destructively modifies the data so create a copy
			smoothData = Arrays.copyOf(data, width * height);

			// Check upper limits are safe
			final int tmpSmooth = FastMath.min((int) smooth, FastMath.min(width, height) / 2);

			// JUnit speed tests show that the rolling median is faster on windows of n<=3
			if (tmpSmooth <= 3)
			{
				if (tmpSmooth <= getBorder())
				{
					filter.rollingMedianInternal(smoothData, width, height, tmpSmooth);
				}
				else
				{
					filter.rollingMedian(smoothData, width, height, tmpSmooth);
				}
			}
			else
			{
				if (tmpSmooth <= getBorder())
				{
					filter.blockMedianInternal(smoothData, width, height, tmpSmooth);
				}
				else
				{
					filter.blockMedian(smoothData, width, height, tmpSmooth);
				}
			}
		}
		return smoothData;
	}

	/**
	 * @return the smoothing width
	 */
	public int getSmooth()
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
		MedianDataProcessor f = (MedianDataProcessor) super.clone();
		// Ensure the object is duplicated and not passed by reference.
		f.filter = (MedianFilter) filter.clone();
		return f;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.DataProcessor#getName()
	 */
	@Override
	public String getName()
	{
		return "Median Filter";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.DataProcessor#getParameters()
	 */
	@Override
	public List<String> getParameters()
	{
		List<String> list = super.getParameters();
		list.add("smooth = " + smooth);
		return list;
	}
}