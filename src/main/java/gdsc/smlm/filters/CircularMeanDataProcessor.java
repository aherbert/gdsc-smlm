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
	 * @return the radius for the desired smoothing distance.
	 */
	public static double getSigma(double smooth)
	{
		if (smooth < 0)
			return 0;
		return smooth;
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
		if (radius > 0)
		{
			// Smoothing destructively modifies the data so create a copy
			smoothData = Arrays.copyOf(data, width * height);
			if (CircularMeanFilter.getBorder(radius) <= getBorder())
				filter.convolveInternal(smoothData, width, height, radius);
			else
				filter.convolve(smoothData, width, height, radius);
		}
		return smoothData;
	}

	/**
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
		return CircularMeanFilter.getPixelRadius(radius) * 2;
	}
}