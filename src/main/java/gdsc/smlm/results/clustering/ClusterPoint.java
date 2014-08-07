package gdsc.smlm.results.clustering;

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
 * Used to store information about a point in the clustering analysis
 */
public class ClusterPoint
{
	public double x, y, weight;
	public int id, t;

	// Used to construct a single linked list of points
	public ClusterPoint next = null;

	/**
	 * Create a cluster point
	 * 
	 * @param id
	 * @param x
	 * @param y
	 * @return The cluster point
	 */
	public static ClusterPoint newClusterPoint(int id, double x, double y)
	{
		return new ClusterPoint(id, x, y);
	}

	/**
	 * Create a cluster point with time information
	 * 
	 * @param id
	 * @param x
	 * @param y
	 * @param t
	 * @return The cluster point
	 */
	public static ClusterPoint newTimeClusterPoint(int id, double x, double y, int t)
	{
		return new ClusterPoint(id, x, y, t);
	}

	/**
	 * Create a cluster point with weight information
	 * 
	 * @param id
	 * @param x
	 * @param y
	 * @param weight
	 * @return The cluster point
	 */
	public static ClusterPoint newClusterPoint(int id, double x, double y, double weight)
	{
		return new ClusterPoint(id, x, y, weight);
	}

	/**
	 * Create a cluster point with weight and time information
	 * 
	 * @param id
	 * @param x
	 * @param y
	 * @param weight
	 * @param t
	 * @return The cluster point
	 */
	public static ClusterPoint newTimeClusterPoint(int id, double x, double y, double weight, int t)
	{
		return new ClusterPoint(id, x, y, weight, t);
	}

	protected ClusterPoint(int id, double x, double y)
	{
		this.id = id;
		this.x = x;
		this.y = y;
		weight = 1;
		t = 0;
	}

	protected ClusterPoint(int id, double x, double y, int t)
	{
		this.id = id;
		this.x = x;
		this.y = y;
		weight = 1;
		this.t = t;
	}

	protected ClusterPoint(int id, double x, double y, double weight)
	{
		this.id = id;
		this.x = x;
		this.y = y;
		this.weight = weight;
		t = 0;
	}

	protected ClusterPoint(int id, double x, double y, double weight, int t)
	{
		this.id = id;
		this.x = x;
		this.y = y;
		this.weight = weight;
		this.t = t;
	}

	public double distance(ClusterPoint other)
	{
		final double dx = x - other.x;
		final double dy = y - other.y;
		return Math.sqrt(dx * dx + dy * dy);
	}

	public double distance2(ClusterPoint other)
	{
		final double dx = x - other.x;
		final double dy = y - other.y;
		return dx * dx + dy * dy;
	}
}