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
	 * Joins the closest pair of particles, one of which must not be in a cluster. Clusters are not joined and can
	 * only grow when particles are added.
	 */
	ParticleSingleLinkage,
	/**
	 * Hierarchical centroid-linkage clustering by joining the closest pair of clusters iteratively
	 */
	Closest,
	/**
	 * Hierarchical centroid-linkage clustering by joining the closest pair of any single particle and another single or
	 * cluster. Clusters are not joined and can only grow when particles are added.
	 */
	ClosestParticle,
	/**
	 * Join the current set of closest pairs in a greedy algorithm. This method computes the pairwise distances and
	 * joins the closest pairs without updating the centroid of each cluster, and the distances, after every join
	 * (centroids and distances are updated after each pass over the data). This can lead to errors over true
	 * hierarchical centroid-linkage clustering where centroid are computed after each link step. For example if A joins
	 * B and C joins D in a single step but the new centroid of AB is closer to C than D.
	 */
	Pairwise,
	/**
	 * A variant of Pairwise is to join the closest pairs only if the number of neighbours for each is
	 * 1. In the event that no pairs has only a single neighbour then only the closest pair is joined.
	 * <p>
	 * In dense images this will return the same results as the Closest algorithm but will be much slower. It may be
	 * faster for sparse density due to the greedy nature of the algorithm.
	 */
	PairwiseWithoutNeighbours,
	/**
	 * Hierarchical centroid-linkage clustering by joining the closest pair of clusters iteratively. Clusters are
	 * compared using time and distance thresholds with priority on the closest time gap (within the distance
	 * threshold).
	 */
	ClosestDistancePriority,
	/**
	 * Hierarchical centroid-linkage clustering by joining the closest pair of clusters iteratively. Clusters are
	 * compared using time and distance thresholds with priority on the closest distance gap (within the time
	 * threshold).
	 */
	ClosestTimePriority,
	/**
	 * Hierarchical centroid-linkage clustering by joining the closest pair of any single particle and another single or
	 * cluster. Clusters are not joined and can only grow when particles are added.
	 * <p>
	 * Clusters are compared using time and distance thresholds with priority on the closest time gap (within the
	 * distance threshold).
	 */
	ClosestParticleDistancePriority,
	/**
	 * Hierarchical centroid-linkage clustering by joining the closest pair of any single particle and another single or
	 * cluster. Clusters are not joined and can only grow when particles are added.
	 * <p>
	 * Clusters are compared using time and distance thresholds with priority on the closest distance gap (within the
	 * time threshold).
	 */
	ClosestParticleTimePriority
}