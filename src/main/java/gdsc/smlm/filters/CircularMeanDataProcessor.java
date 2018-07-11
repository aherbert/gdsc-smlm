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

import gdsc.core.ij.Utils;

/**
 * Identifies candidate spots (local maxima) in an image. The image is smoothed with a circular mean filter.
 */
public class CircularMeanDataProcessor extends DataProcessor
{
	private final double radius;
	private CircularMeanFilter filter;

	/**
	 * Constructor
	 *
	 * @param border
	 *            The border to ignore for maxima
	 * @param smooth
	 *            The distance into neighbouring pixels to extend
	 */
	public CircularMeanDataProcessor(int border, double smooth)
	{
		super(border);
		this.radius = getSigma(smooth);
		filter = new CircularMeanFilter();
	}

	/**
	 * Get the radius for the desired smoothing distance.
	 *
	 * @param smooth
	 *            the smoothing distance
	 * @return the radius for the desired smoothing distance.
	 */
	public static double getSigma(double smooth)
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
		return true;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.filters.DataProcessor#setWeights(float[], int, int)
	 */
	@Override
	public void setWeights(float[] weights, int width, int height)
	{
		filter.setWeights(weights, width, height);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.filters.DataProcessor#hasWeights()
	 */
	@Override
	public boolean hasWeights()
	{
		return filter.hasWeights();
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
		if (radius > 0)
		{
			// Smoothing destructively modifies the data so create a copy
			smoothData = Arrays.copyOf(data, width * height);
			if (CircularFilter.getBorder(radius) <= getBorder())
				filter.convolveInternal(smoothData, width, height, radius);
			else
				filter.convolve(smoothData, width, height, radius);
		}
		return smoothData;
	}

	/**
	 * Gets the radius.
	 *
	 * @return the smoothing radius
	 */
	public double getRadius()
	{
		return radius;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see java.lang.Object#clone()
	 */
	@Override
	public CircularMeanDataProcessor clone()
	{
		CircularMeanDataProcessor f = (CircularMeanDataProcessor) super.clone();
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
		return "Circular Mean";
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
		list.add("radius = " + Utils.rounded(radius));
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
		return CircularFilter.getPixelRadius(radius) * 2;
	}
}
