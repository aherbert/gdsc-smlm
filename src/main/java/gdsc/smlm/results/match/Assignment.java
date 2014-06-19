package gdsc.smlm.results.match;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Stores an assignment between two identified points and the distance between them
 */
public class Assignment implements Comparable<Assignment>
{
	private int targetId;
	private int predictedId;
	private double distance;

	public Assignment(int targetId, int predictedId, double distance)
	{
		this.targetId = targetId;
		this.predictedId = predictedId;
		this.distance = distance;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo(Assignment o)
	{
		if (this.distance < o.distance)
			return -1;
		if (this.distance > o.distance)
			return 1;
		return 0;
	}

	/**
	 * @return the targetId
	 */
	public int getTargetId()
	{
		return targetId;
	}

	/**
	 * @return the predictedId
	 */
	public int getPredictedId()
	{
		return predictedId;
	}

	/**
	 * @return the distance
	 */
	public double getDistance()
	{
		return distance;
	}
}