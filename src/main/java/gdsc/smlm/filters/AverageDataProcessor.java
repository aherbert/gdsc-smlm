package gdsc.smlm.filters;

import gdsc.core.ij.Utils;

import java.util.Arrays;
import java.util.List;

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
 * Identifies candidate spots (local maxima) in an image. The image is smoothed with an average box filter. When the box
 * size is large then the smoothing switches to using an interpolation between two block sizes. This is an approximation
 * due to incorrect weighting of the corners. Note that at large sizes the corners are a fraction of the total edge
 * pixels so the difference is minor.
 * 
 * @see gdsc.smlm.filters.AreaAverageFilter
 */
public class AverageDataProcessor extends DataProcessor
{
	/**
	 * Define the default smoothing size at which the smoothing switches to using an interpolation between two block
	 * sizes. Below this limit the smoothing uses an exact mean filter.
	 */
	public static final int AREA_FILTER_LIMIT = 3;

	private final float smooth;
	private final int iSmooth;
	private AverageFilter filter = null;
	private AreaAverageFilter areaFilter = null;

	/**
	 * Constructor
	 * 
	 * @param border
	 *            The border to ignore for maxima
	 * @param smooth
	 *            The smoothing width to apply to the data
	 */
	public AverageDataProcessor(int border, double smooth)
	{
		this(border, smooth, AREA_FILTER_LIMIT);
	}

	/**
	 * Constructor
	 * 
	 * @param border
	 *            The border to ignore for maxima
	 * @param smooth
	 *            The smoothing width to apply to the data
	 * @param areaFilterLimit
	 *            The limit to switch to simple area interpolation
	 */
	public AverageDataProcessor(int border, double smooth, int areaFilterLimit)
	{
		super(border);
		this.smooth = (float) convert(smooth);
		// Store the smoothing value as an integer
		iSmooth = ((int) smooth == smooth) ? (int) smooth : 0;

		// Only create the area filter if we need it
		if (iSmooth > 0)
		{
			filter = new AverageFilter();
		}
		else
		{
			if (smooth > areaFilterLimit)
				areaFilter = new AreaAverageFilter();
			else
				filter = new AverageFilter();
		}
	}

	/**
	 * Convert the smoothing parameter to the value which is used for the AreaAverageFilter.
	 * Values below zero are set to zero.
	 * 
	 * @param smooth
	 * @return The adjusted value
	 */
	public static double convert(double smooth)
	{
		if (smooth < 0)
			return 0;
		return smooth;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.DataProcessor#isWeighted()
	 */
	@Override
	public boolean isWeighted()
	{
		return false;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.DataProcessor#setWeights(float[], int, int)
	 */
	@Override
	public void setWeights(float[] weights, int width, int height)
	{

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.DataProcessor#hasWeights()
	 */
	@Override
	public boolean hasWeights()
	{
		return false;
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

			if (iSmooth > 1)
			{
				// Integer smoothing is faster using a rolling block algorithm
				if (smooth <= getBorder())
				{
					filter.rollingBlockAverageInternal(smoothData, width, height, iSmooth);
				}
				else
				{
					filter.rollingBlockAverage(smoothData, width, height, iSmooth);
				}
			}
			else
			{
				// Float smoothing must use the striped block algorithm or the area average filter.
				// The area average filter is faster when above a 7x7 block size.
				// At this point the difference between the filters is small (the area average filter
				// is biased to the corners) so switch to the faster filter.

				if (areaFilter != null)
				{
					if (smooth <= getBorder())
					{
						areaFilter.areaAverageUsingSumsInternal(smoothData, width, height, smooth);
					}
					else
					{
						areaFilter.areaAverageUsingSums(smoothData, width, height, smooth);
					}
				}
				else
				{
					if (smooth <= getBorder())
					{
						filter.stripedBlockAverageInternal(smoothData, width, height, smooth);
					}
					else
					{
						filter.stripedBlockAverage(smoothData, width, height, smooth);
					}
				}
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
	public AverageDataProcessor clone()
	{
		AverageDataProcessor f = (AverageDataProcessor) super.clone();
		// Ensure the object is duplicated and not passed by reference.
		if (filter != null)
			f.filter = filter.clone();
		else
			f.areaFilter = areaFilter.clone();
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
		return "Average";
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
		list.add("smooth = " + Utils.rounded(smooth));
		return list;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.DataProcessor#getSpread()
	 */
	@Override
	public double getSpread()
	{
		return 2 * smooth + 1;
	}
}