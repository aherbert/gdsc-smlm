/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package gdsc.smlm.filters;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

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
	 * @param smooth
	 *            the smoothing parameter
	 * @return The adjusted value
	 * @see MedianFilter
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
		// Do nothing
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

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.DataProcessor#process(float[], int, int)
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
			final int tmpSmooth = FastMath.min(smooth, FastMath.min(width, height) / 2);

			// JUnit speed tests show that the rolling median is not faster.
			// It used to be faster on windows less than 3.

			//if (tmpSmooth <= 3)
			//{
			//	if (tmpSmooth <= getBorder())
			//	{
			//		filter.rollingMedianInternal(smoothData, width, height, tmpSmooth);
			//	}
			//	else
			//	{
			//		filter.rollingMedian(smoothData, width, height, tmpSmooth);
			//	}
			//}
			//else
			//{
			if (tmpSmooth <= getBorder())
			{
				filter.blockMedianInternal(smoothData, width, height, tmpSmooth);
			}
			else
			{
				filter.blockMedian(smoothData, width, height, tmpSmooth);
			}
			//}
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
	@Override
	public MedianDataProcessor clone()
	{
		MedianDataProcessor f = (MedianDataProcessor) super.clone();
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
		return "Median";
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
