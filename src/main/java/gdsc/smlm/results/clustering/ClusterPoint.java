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
	public double x, y;
	public int id, t;
	
	// Used to construct a single linked list of points
	public ClusterPoint next = null;

	public ClusterPoint(int id, double x, double y)
	{
		this.id = id;
		this.x = x;
		this.y = y;
		t = 0;
	}

	public ClusterPoint(int id, double x, double y, int t)
	{
		this.id = id;
		this.x = x;
		this.y = y;
		this.t = t;
	}
}