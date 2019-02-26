package gdsc.smlm.filters;

import uk.ac.sussex.gdsc.core.ij.ImageJUtils; import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils; import uk.ac.sussex.gdsc.core.utils.TextUtils; import uk.ac.sussex.gdsc.core.utils.MathUtils;

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
 * Identifies candidate spots (local maxima) in an image. The image is smoothed with an average box filter.
 */
public class BlockAverageDataProcessor extends DataProcessor
{
	private final double smooth;
	private AverageFilter filter = new AverageFilter();

	/**
	 * Constructor
	 * 
	 * @param border
	 *            The border to ignore for maxima
	 * @param smooth
	 *            The smoothing width to apply to the data
	 */
	public BlockAverageDataProcessor(int border, double smooth)
	{
		super(border);
		this.smooth = convert(smooth);
	}

	/**
	 * Convert the smoothing parameter to the value which is used for the AverageFilter.
	 * We only use int smoothing. Values below zero are set to zero.
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
		return (int) smooth;
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
	public BlockAverageDataProcessor clone()
	{
		BlockAverageDataProcessor f = (BlockAverageDataProcessor) super.clone();
		// Ensure the object is duplicated and not passed by reference.
		f.filter = filter.clone();
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
		return "Block Average";
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
		list.add("smooth = " + MathUtils.rounded(smooth));
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