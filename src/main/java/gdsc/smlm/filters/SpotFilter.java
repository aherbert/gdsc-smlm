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
 * Identifies candidate spots (local maxima) in an image.
 */
public abstract class SpotFilter implements Cloneable
{
	/**
	 * Find the candidate spots in the data
	 * 
	 * @param data
	 *            The data
	 * @param width
	 *            The width of the data
	 * @param height
	 *            The height of the data
	 * @return The candidate spots
	 */
	protected abstract Spot[] find(final float[] data, final int width, final int height);

	/**
	 * List the candidate spots in the data. The list will be in the order the candidates are found.
	 * 
	 * @param data
	 *            The data
	 * @param width
	 *            The width of the data
	 * @param height
	 *            The height of the data
	 * @return The candidate spots (may be an empty array but will not be null)
	 */
	public Spot[] list(float[] data, int width, int height)
	{
		Spot[] spots = find(data, width, height);
		return (spots == null) ? new Spot[0] : spots;
	}

	/**
	 * List and then rank the candidate spots in the data. The list will be in the order defined by sorting the
	 * candidates.
	 * 
	 * @param data
	 *            The data
	 * @param width
	 *            The width of the data
	 * @param height
	 *            The height of the data
	 * @return The candidate spots (may be an empty array but will not be null)
	 */
	public Spot[] rank(float[] data, int width, int height)
	{
		Spot[] spots = find(data, width, height);
		if (spots == null)
			return new Spot[0];
		Arrays.sort(spots);
		return spots;
	}

	/**
	 * Return true if the intensity value of the candidate spots is absolute, i.e. is the height of the candidate
	 * maximum using the original data scale. Otherwise the intensity is relative, for example this could be relative to
	 * the local background.
	 * 
	 * @return True if the intensity value of the candidate spots is absolute
	 */
	public abstract boolean isAbsoluteIntensity();

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public Object clone()
	{
		try
		{
			return (SpotFilter) super.clone();
		}
		catch (CloneNotSupportedException e)
		{
			return null;
		}
	}

	/**
	 * @return A description of the filter and parameters
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
	public abstract List<String> getParameters();

	/**
	 * Get the width spread of data used to process each position
	 * 
	 * @return The spread
	 */
	public abstract double getSpread();
}