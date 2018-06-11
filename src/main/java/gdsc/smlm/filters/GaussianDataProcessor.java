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

import gdsc.core.ij.Utils;

import java.util.Arrays;
import java.util.List;


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
	 */
	public GaussianDataProcessor(int border, double smooth)
	{
		super(border);
		this.sigma = getSigma(smooth);
		filter = new GaussianFilter(0.02);
	}

	/**
	 * Get the Gaussian standard deviation for the desired smoothing distance.
	 * 
	 * @return the Gaussian standard deviation for the desired smoothing distance.
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
	public GaussianDataProcessor clone()
	{
		GaussianDataProcessor f = (GaussianDataProcessor) super.clone();
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
		return "Gaussian";
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

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.DataProcessor#getSpread()
	 */
	@Override
	public double getSpread()
	{
		return 6 * sigma;
	}
}
