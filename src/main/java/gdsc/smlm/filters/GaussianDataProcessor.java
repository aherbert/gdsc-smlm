package gdsc.smlm.filters;

import gdsc.smlm.ij.utils.Utils;

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
 * Identifies candidate spots (local maxima) in an image. The image is smoothed with a Gaussian filter.
 */
public class GaussianDataProcessor extends DataProcessor
{
	private final double sigma;
	private GaussianFilter filter;

	/**
	 * Constructor
	 * 
	 * @param border
	 *            The border to ignore for maxima
	 * @param smooth
	 *            The distance into neighbouring pixels to extend. The resulting standard deviation can be found using
	 *            {@link #getSigma()}
	 * @throws IllegalArgumentException
	 *             if smooth is below zero
	 */
	public GaussianDataProcessor(int border, double smooth)
	{
		super(border);
		// Make the 3 SD reach the desired smoothing distance. We add 0.5 pixels to always cover the target pixel. 		
		this.sigma = (smooth + 0.5) / 3;
		filter = new GaussianFilter(0.02);
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
		if (sigma > 0)
		{
			// Smoothing destructively modifies the data so create a copy
			smoothData = Arrays.copyOf(data, width * height);
			if (GaussianFilter.getBorder(sigma) <= getBorder())
				filter.convolveInternal(smoothData, width, height, sigma);
			else
				filter.convolve(smoothData, width, height, sigma);
		}
		return smoothData;
	}

	/**
	 * @return the Gaussian standard deviation
	 */
	public double getSigma()
	{
		return sigma;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public Object clone()
	{
		GaussianDataProcessor f = (GaussianDataProcessor) super.clone();
		// Ensure the object is duplicated and not passed by reference.
		f.filter = (GaussianFilter) filter.clone();
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
		return "Gaussian Filter";
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
		list.add("sigma = " + Utils.rounded(sigma));
		list.add("width = " + Utils.rounded(filter.getHalfWidth(sigma)));
		return list;
	}
}