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
public class DifferenceSpotFilter extends MaximaSpotFilter
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
	public DifferenceSpotFilter(int search, int border, DataProcessor processor1, DataProcessor processor2)
	{
		super(search, border);
		if (processor1 == null)
			throw new IllegalArgumentException("Processor 1 is null");
		if (processor2 == null)
			throw new IllegalArgumentException("Processor 2 is null");
		// TODO : This is a simple protection from invalid difference-of-smoothing. It could be improved.
		if (processor2.getSpread() < processor1.getSpread())
			throw new IllegalArgumentException("Processor 2 acts on a smaller spread of data than processor 1");
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
		return false;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.SpotFilter#isWeighted()
	 */
	@Override
	public boolean isWeighted()
	{
		return processor1.isWeighted() || processor2.isWeighted();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.SpotFilter#setWeights(float[], int, int)
	 */
	@Override
	public void setWeights(float[] weights, int width, int height)
	{
		processor1.setWeights(weights, width, height);
		processor2.setWeights(weights, width, height);
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
	public DifferenceSpotFilter clone()
	{
		DifferenceSpotFilter f = (DifferenceSpotFilter) super.clone();
		// Ensure the object is duplicated and not passed by reference.
		f.processor1 = processor1.clone();
		f.processor2 = processor2.clone();
		return f;
	}

	/* (non-Javadoc)
	 * @see gdsc.smlm.filters.SpotFilter#getName()
	 */
	@Override
	public String getName()
	{
		return "Difference";
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

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.SpotFilter#getSpread()
	 */
	@Override
	public double getSpread()
	{
		return Math.max(processor1.getSpread(), processor2.getSpread());
	}
}