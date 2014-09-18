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
 * Extends the {@link gdsc.smlm.results.clustering.ClusterPoint } class to hold additional information for use in
 * clustering
 */
public class ExtendedClusterPoint extends ClusterPoint
{
	public ExtendedClusterPoint nextE = null;
	public boolean inCluster = false;

	public ExtendedClusterPoint(int id, double x, double y, int t, ClusterPoint other)
	{
		super(id, x, y, t, t);
		super.next = other;
	}
}