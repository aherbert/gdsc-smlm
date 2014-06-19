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
 * Used to store all the information about a cluster in the clustering analysis
 */
public class Cluster implements Comparable<Cluster>
{
	public double x, y, sumx, sumy;
	public int n;

	// Used to construct a single linked list of clusters
	public Cluster next = null;

	// Used to store potential clustering links
	public Cluster closest = null;
	public double d2;
	
	// Used to indicate this cluster has a neighbour
	public int neighbour = 0;

	// Used to construct a single linked list of cluster points
	public ClusterPoint head = null;
	public int xBin, yBin;

	public Cluster(ClusterPoint point)
	{
		point.next = null;
		head = point;
		this.x = sumx = point.x;
		this.y = sumy = point.y;
		n = 1;
	}

	public double distance(Cluster other)
	{
		final double dx = x - other.x;
		final double dy = y - other.y;
		return Math.sqrt(dx * dx + dy * dy);
	}

	public double distance2(Cluster other)
	{
		final double dx = x - other.x;
		final double dy = y - other.y;
		return dx * dx + dy * dy;
	}

	public void add(Cluster other)
	{
		// Do not check if the other cluster is null or has no points

		// Add to this list
		// Find the tail of the shortest list
		ClusterPoint big, small;
		if (n < other.n)
		{
			small = head;
			big = other.head;
		}
		else
		{
			small = other.head;
			big = head;
		}

		ClusterPoint tail = small;
		while (tail.next != null)
			tail = tail.next;

		// Join the small list to the long list 
		tail.next = big;
		head = small;

		// Find the new centroid
		sumx += other.sumx;
		sumy += other.sumy;
		n += other.n;
		x = sumx / n;
		y = sumy / n;

		// Free the other cluster
		other.clear();
	}

	protected void clear()
	{
		head = null;
		closest = null;
		n = 0;
		x = y = sumx = sumy = d2 = 0;
	}

	/**
	 * Link the two clusters as potential merge candidates only if the squared distance is smaller than the other
	 * clusters current closest
	 * 
	 * @param other
	 * @param d2
	 */
	public void link(Cluster other, double d2)
	{
		// Check if the other cluster has a closer candidate
		if (other.closest != null && other.d2 < d2)
			return;

		other.closest = this;
		other.d2 = d2;

		this.closest = other;
		this.d2 = d2;
	}

	/**
	 * @return True if the closest cluster links back to this cluster
	 */
	public boolean validLink()
	{
		// Check if the other cluster has an updated link to another candidate
		if (closest != null)
		{
			// Valid if the other cluster links back to this cluster
			return closest.closest == this;
		}
		return false;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo(Cluster o)
	{
		// Sort by size then centroid x then y. 
		// The sort is arbitrary but allows comparison of two lists after sorting.
		if (n < o.n)
			return -1;
		if (n > o.n)
			return 1;
		if (x < o.x)
			return -1;
		if (x > o.x)
			return 1;
		if (y < o.y)
			return -1;
		if (y > o.y)
			return 1;
		return 0;
	}
}