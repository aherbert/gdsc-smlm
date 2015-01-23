package gdsc.smlm.filters;

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
 * Identifies candidate spots (local maxima) in an image. The image is smoothed twice with an average box filter and the
 * result of the large smoothing subtracted from the smaller smoothing. A simple weighted approximation is used for less
 * than 1 pixel smoothing.
 * <p>
 * Difference-of-smoothing (Top-Hat Box filter) See: E. J. Breen, G. H. Joss, and K. L. Williams, “Locating objects of
 * interest within biological images: The top hat box ﬁlter,” J. Comput.-Assist. Microsc., vol. 3, no. 2, pp. 97–102,
 * 1991. 
 */
public class TopHatBoxSpotFilter extends AverageSpotFilter
{
	private final double smooth2;

	/**
	 * Constructor
	 * 
	 * @param search
	 *            The search width for non-maximum suppression
	 * @param border
	 *            The border to ignore for maxima
	 * @param smooth
	 *            The first smoothing width to apply to the data
	 * @param smooth2
	 *            The second smoothing width to apply to the data
	 * @throws IllegalArgumentException
	 *             if smooth is below zero, if smooth2 is smaller than or equal to smooth
	 */
	public TopHatBoxSpotFilter(int search, int border, double smooth, double smooth2)
	{
		super(search, border, smooth);
		smooth2 = convert(smooth2);
		if (smooth2 <= smooth)
			throw new IllegalArgumentException("The second smoothing width must be larger than the first");
		this.smooth2 = smooth2;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.SpotFilter#isAbsoluteIntensity()
	 */
	public boolean isAbsoluteIntensity()
	{
		return false;
	}

	@Override
	public float[] preprocessData(float[] data, int width, int height)
	{
		final float[] smoothData;
		if (getSmooth() > 0)
		{
			// Single smoothing (Box filter)
			smoothData = applySmoothing(data, width, height, getSmooth());
		}
		else
		{
			// No smoothing so copy the data 
			smoothData = Arrays.copyOf(data, width * height);
		}

		// The constructor check ensured that smooth2 is bigger than smooth 
		// so we subtract the second smoothed image
		final float[] smooth2Data = applySmoothing(data, width, height, smooth2);
		for (int i = 0; i < smoothData.length; i++)
		{
			smoothData[i] -= smooth2Data[i];
		}
		return smoothData;
	}

	/**
	 * @return the second smoothing width
	 */
	public double getSmooth2()
	{
		return smooth2;
	}

	@Override
	public String getName()
	{
		return "Top-Hat Filter";
	}
	
	@Override
	public List<String> getParameters()
	{
		List<String> list = super.getParameters();
		list.add("smooth2 = " + smooth2);
		return list;
	}
}