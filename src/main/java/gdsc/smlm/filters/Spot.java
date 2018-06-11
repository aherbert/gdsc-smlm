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

/**
 * Identify a candidate spot (local maximum).
 */
public class Spot implements Comparable<Spot>, Cloneable
{
	final public int x, y;
	final public float intensity;
	private float score;

	/**
	 * Constructor that sets the score (for sorting) equal to the intensity
	 * 
	 * @param x
	 *            The x-coordinate
	 * @param y
	 *            The y-coordinate
	 * @param intensity
	 *            The intensity of the spot
	 */
	public Spot(int x, int y, float intensity)
	{
		this.x = x;
		this.y = y;
		this.intensity = score = intensity;
	}

	/**
	 * @param x
	 *            The x-coordinate
	 * @param y
	 *            The y-coordinate
	 * @param intensity
	 *            The intensity of the spot
	 * @param score
	 *            The score used for sorting
	 */
	public Spot(int x, int y, float intensity, float score)
	{
		this.x = x;
		this.y = y;
		this.intensity = intensity;
		this.score = score;
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
		if (score > o.score)
			return -1;
		if (score < o.score)
			return 1;
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	public Spot clone()
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

	/**
	 * @return the score used for sorting the spots
	 */
	public float getScore()
	{
		return score;
	}

	/**
	 * @param score
	 *            the score to set for sorting the spots
	 */
	public void setScore(float score)
	{
		this.score = score;
	}
}
