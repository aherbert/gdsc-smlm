package gdsc.smlm.filters;

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
 * Identifies candidate spots (local maxima) in an image. The image is pre-processed with two filters and the second
 * subtracted from the first.
 */
public class DoublePassSpotFilter extends MaximaSpotFilter
{
	private DataProcessor processor1;
	private DataProcessor processor2;

	/**
	 * Constructor
	 * 
	 * @param search
	 *            The search width for non-maximum suppression
	 * @param border
	 *            The border to ignore for maxima
	 * @param processor1
	 *            The first data processor
	 * @param processor2
	 *            The second data processor
	 * @throws IllegalArgumentException
	 *             if either processor is null
	 */
	public DoublePassSpotFilter(int search, int border, DataProcessor processor1, DataProcessor processor2)
	{
		super(search, border);
		if (processor1 == null)
			throw new IllegalArgumentException("Processor 1 is null");
		if (processor2 == null)
			throw new IllegalArgumentException("Processor 2 is null");
		this.processor1 = processor1;
		this.processor2 = processor2;
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
	 * @see gdsc.smlm.filters.MaximaSpotFilter#preprocessData(float[], int, int)
	 */
	@Override
	public float[] preprocessData(float[] data, int width, int height)
	{
		final float[] data1 = processor1.process(data, width, height);
		final float[] data2 = processor2.process(data, width, height);
		for (int i = 0; i < data1.length; i++)
		{
			data1[i] -= data2[i];
		}
		return data1;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public Object clone()
	{
		DoublePassSpotFilter f = (DoublePassSpotFilter) super.clone();
		// Ensure the object is duplicated and not passed by reference.
		f.processor1 = (DataProcessor) processor1.clone();
		f.processor2 = (DataProcessor) processor2.clone();
		return f;
	}

	/* (non-Javadoc)
	 * @see gdsc.smlm.filters.SpotFilter#getName()
	 */
	@Override
	public String getName()
	{
		return "Double-Pass Filter";
	}

	/* (non-Javadoc)
	 * @see gdsc.smlm.filters.MaximaSpotFilter#getParameters()
	 */
	@Override
	public List<String> getParameters()
	{
		List<String> list = super.getParameters();
		list.add("Filter 1 = " + processor1.getDescription());
		list.add("Filter 2 = " + processor2.getDescription());
		return list;
	}
}