package gdsc.smlm.filters;

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
 * Identify a candidate spot (local maximum).
 */
public class Spot implements Comparable<Spot>, Cloneable
{
	final public int x, y;
	final public float intensity;

	public Spot(int x, int y, float intensity)
	{
		this.x = x;
		this.y = y;
		this.intensity = intensity;
	}

	/**
	 * Get the distance between the two spots
	 * 
	 * @param o
	 *            the other spot
	 * @return The distance
	 */
	public double distance(Spot o)
	{
		return Math.sqrt(distance2(o));
	}

	/**
	 * Get the squared distance between the two spots
	 * 
	 * @param o
	 *            the other spot
	 * @return The squared distance
	 */
	public double distance2(Spot o)
	{
		final int dx = x - o.x;
		final int dy = y - o.y;
		return dx * dx + dy * dy;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	@Override
	public int compareTo(Spot o)
	{
		if (intensity > o.intensity)
			return -1;
		if (intensity < o.intensity)
			return 1;
		return 0;
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public Object clone()
	{
		try
		{
			return (Spot) super.clone();
		}
		catch (CloneNotSupportedException e)
		{
			// Ignore
		}
		return null;
	}
}