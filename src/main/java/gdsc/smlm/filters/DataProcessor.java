package gdsc.smlm.filters;

import java.util.ArrayList;
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
 * Define a class to pre-process the data, ignoring a defined border.
 */
public abstract class DataProcessor implements Cloneable
{
	private final int border;

	/**
	 * @param border
	 *            The border that can be ignored
	 */
	public DataProcessor(int border)
	{
		this.border = border;
	}

	/**
	 * Checks if the data processor is weighted, i.e. supports {@link #setWeights(float[], int, int)}.
	 *
	 * @return true, if is weighted
	 */
	public abstract boolean isWeighted();
	
	/**
	 * Sets the weights of the data. This should be called before {@link #process(float[], int, int)} is called with
	 * data samples.
	 * <p>
	 * Calling this in advance allows efficient caching of pre-computed weightings.
	 *
	 * @param weights
	 *            the weights of the data (can be null)
	 * @param width
	 *            The width of the data
	 * @param height
	 *            The height of the data
	 */
	public abstract void setWeights(final float[] weights, final int width, final int height);

	/**
	 * Process the data
	 * 
	 * @param data
	 *            The data
	 * @param width
	 *            The width of the data
	 * @param height
	 *            The height of the data
	 * @return The new data
	 */
	protected abstract float[] process(final float[] data, final int width, final int height);

	/**
	 * @return the border
	 */
	public int getBorder()
	{
		return border;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public DataProcessor clone()
	{
		try
		{
			return (DataProcessor) super.clone();
		}
		catch (CloneNotSupportedException e)
		{
			return null;
		}
	}

	/**
	 * @return A description of the processor and parameters
	 */
	public String getDescription()
	{
		return getName() + ": " + Arrays.toString(getParameters().toArray());
	}

	/**
	 * @return The name of the filter
	 */
	public abstract String getName();

	/**
	 * @return The parameters of the filter
	 */
	public List<String> getParameters()
	{
		ArrayList<String> list = new ArrayList<String>();
		list.add("border = " + border);
		return list;
	}

	/**
	 * Get the width spread of data used to process each position
	 * 
	 * @return The spread
	 */
	public abstract double getSpread();
}