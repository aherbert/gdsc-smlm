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
	public int id, start, end;

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
	 * @param start
	 * @param end
	 * @return The cluster point
	 */
	public static ClusterPoint newTimeClusterPoint(int id, double x, double y, int start, int end)
	{
		return new ClusterPoint(id, x, y, start, end);
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
	 * @param start
	 * @param end
	 * @return The cluster point
	 */
	public static ClusterPoint newTimeClusterPoint(int id, double x, double y, double weight, int start, int end)
	{
		return new ClusterPoint(id, x, y, weight, start, end);
	}

	protected ClusterPoint(int id, double x, double y)
	{
		this.id = id;
		this.x = x;
		this.y = y;
		weight = 1;
		start = end = 0;
	}

	protected ClusterPoint(int id, double x, double y, int start, int end)
	{
		this.id = id;
		this.x = x;
		this.y = y;
		weight = 1;
		this.start = start;
		this.end = end;
	}

	protected ClusterPoint(int id, double x, double y, double weight)
	{
		this.id = id;
		this.x = x;
		this.y = y;
		this.weight = weight;
		start = end = 0;
	}

	protected ClusterPoint(int id, double x, double y, double weight, int start, int end)
	{
		this.id = id;
		this.x = x;
		this.y = y;
		this.weight = weight;
		this.start = start;
		this.end = end;
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

	/**
	 * Get the time gap between the two points. If the points overlap then return 0.
	 * 
	 * @param other
	 * @return the time gap
	 */
	public int gap(ClusterPoint other)
	{
		// Overlap:
		// S-----------E
		//         S---------E
		//
		// Gap:
		// S-----------E
		//                  S---------E
		return Math.max(0, Math.max(start, other.start) - Math.min(end, other.end));
	}
}