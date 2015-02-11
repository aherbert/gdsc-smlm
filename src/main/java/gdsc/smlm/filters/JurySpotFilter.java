package gdsc.smlm.filters;

import gdsc.smlm.utils.NotImplementedException;

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
 * Identifies candidate spots (local maxima) in an image. The image is pre-processed with a collection of filters and
 * the combined height used to identify candidates.
 */
public final class JurySpotFilter extends MaximaSpotFilter
{
	private DataProcessor[] processors;

	/**
	 * Constructor
	 * 
	 * @param search
	 *            The search width for non-maximum suppression
	 * @param border
	 *            The border to ignore for maxima
	 * @param processor
	 *            The data processor
	 * @throws IllegalArgumentException
	 *             if processor is null
	 */
	public JurySpotFilter(int search, int border, DataProcessor... processors)
	{
		super(search, border);
		if (processors == null || processors.length == 0)
			throw new IllegalArgumentException("No processors");
		for (int i = 0; i < processors.length; i++)
			if (processors[i] == null)
				throw new IllegalArgumentException("Processor " + (i + 1) + " is null");
		this.processors = processors;
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

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.MaximaSpotFilter#find(float[], int, int)
	 */
	@Override
	protected Spot[] find(float[] data, int width, int height)
	{
		// Run all the processors and store the total maxima intensity at each index
		final float[] intensity = new float[width * height];
		final float[] sum = new float[intensity.length];
		for (int i = 0; i < processors.length; i++)
		{
			final float[] data2 = processors[i].process(data, width, height);
			for (int j = 0; j < data2.length; j++)
				sum[j] += data2[j];
			final int[] maxIndices = getMaxima(data2, width, height);
			for (final int index : maxIndices)
			{
				intensity[index] += data2[index];
			}
		}

		// Note: A simple jury using the any non-zero point in the maxima intensity will work if the background 
		// noise is uniform. If the noise is sloped then larger filters may result in the centre of the maxima 
		// moving and the jury will return two smaller candidates in adjacent pixels. To mitigate this we get the 
		// maxima again. 
		final int[] maxIndices = getMaxima(intensity, width, height);
		if (maxIndices.length == 0)
			return null;

		// Normalise the intensity across all processors
		final float divisor = (float) (1.0 / processors.length);

		final Spot[] spots = new Spot[maxIndices.length];
		for (int n = 0; n < maxIndices.length; n++)
		{
			final int y = maxIndices[n] / width;
			final int x = maxIndices[n] % width;
			spots[n] = new Spot(x, y, sum[maxIndices[n]] * divisor);
		}
		return spots;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.MaximaSpotFilter#preprocessData(float[], int, int)
	 */
	@Override
	public float[] preprocessData(float[] data, int width, int height)
	{
		// Not needed as we override the filter method directly
		throw new NotImplementedException();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public Object clone()
	{
		JurySpotFilter f = (JurySpotFilter) super.clone();
		// Ensure the object is duplicated and not passed by reference.
		f.processors = processors.clone();
		for (int i = 0; i < processors.length; i++)
			f.processors[i] = (DataProcessor) processors[i].clone();
		return f;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.SpotFilter#getName()
	 */
	@Override
	public String getName()
	{
		return "Jury Filter";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.MaximaSpotFilter#getParameters()
	 */
	@Override
	public List<String> getParameters()
	{
		List<String> list = super.getParameters();
		for (int i = 0; i < processors.length; i++)
			list.add("Filter = " + processors[i].getDescription());
		return list;
	}
}