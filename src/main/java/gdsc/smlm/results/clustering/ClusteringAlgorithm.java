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
 * Define the clustering algorithm
 */
public enum ClusteringAlgorithm
{
	/**
	 * Hierarchical clustering by joining the closest pair of clusters iteratively
	 */
	Closest,
	/**
	 * Join the current set of closest pairs in a greedy algorithm. This method computes the pairwise distances and
	 * joins the closest pairs without updating the centroid of each cluster, and the distances, after every join
	 * (centroids and distances are updated after each pass over the data). This can lead to errors over true
	 * hierarchical clustering where centroid are computed after each link step. For example if A joins B and C joins D
	 * in a single step but the new centroid of AB is closer to C than D.
	 */
	PairwiseWithoutNeighbours,
	/**
	 * A variant of PairwiseWithoutNeighbours is to join the closest pairs only if the number of neighbours for each is
	 * 1. In the event that no pairs has only a single neighbour then only the closest pair is joined.
	 */
	Pairwise,
	/**
	 * Hierarchical clustering by joining the closest pair of clusters iteratively. Clusters are compared using time and
	 * distance thresholds with priority on the closest time gap (within the distance threshold). 
	 */
	ClosestDistancePriority,
	/**
	 * Hierarchical clustering by joining the closest pair of clusters iteratively. Clusters are compared using time and
	 * distance thresholds with priority on the closest distance gap (within the time threshold). 
	 */
	ClosestTimePriority
}