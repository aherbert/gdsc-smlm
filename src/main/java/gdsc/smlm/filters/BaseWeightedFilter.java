package gdsc.smlm.filters;

// TODO: Auto-generated Javadoc
/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Contains common functionality for weighted filters.
 */
public abstract class BaseWeightedFilter extends BaseFilter
{

	/** The weights. */
	protected float[] weights;

	/** The height. */
	protected int width, height;

	/**
	 * Sets the weights of the data. This should be called before filtering data samples.
	 *
	 * @param weights
	 *            the weights of the data (can be null)
	 * @param width
	 *            The width of the data
	 * @param height
	 *            The height of the data
	 */
	public void setWeights(final float[] weights, final int width, final int height)
	{
		this.weights = weights;
		this.width = width;
		this.height = height;
		newWeights();
	}

	/**
	 * Checks for weights.
	 *
	 * @return true, if successful
	 */
	public boolean hasWeights()
	{
		return weights != null;
	}

	/**
	 * Signal that new weight parameters have been set. Sub-classes can re-initialise for the new weights.
	 */
	protected abstract void newWeights();

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public BaseWeightedFilter clone()
	{
		return (BaseWeightedFilter) super.clone();
	}
}