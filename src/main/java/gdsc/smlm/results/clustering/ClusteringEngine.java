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

import gdsc.smlm.results.DensityManager;
import gdsc.smlm.results.NullTrackProgress;
import gdsc.smlm.results.TrackProgress;

import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Find clusters of points using a partial centroid-linkage hierarchical clustering algorithm.
 * <p>
 * Points are added to the nearest cluster if they are below the distance threshold to the cluster centroid. The cluster
 * centroid is updated. All points above the cluster distance threshold remain as single molecules.
 */
public class ClusteringEngine
{
	private ClusteringAlgorithm clusteringAlgorithm = ClusteringAlgorithm.Pairwise;
	private TrackProgress tracker;
	private int pulseInterval = 0;
	private int nextClusterId;

	public ClusteringEngine()
	{
		setTracker(null);
	}

	/**
	 * Find the clusters of points within the specified radius
	 * 
	 * @param points
	 * @param radius
	 * @return the clusters
	 */
	public ArrayList<Cluster> findClusters(List<ClusterPoint> points, double radius)
	{
		return findClusters(points, radius, 0);
	}

	/**
	 * Find the clusters of points within the specified radius and time separation
	 * 
	 * @param points
	 * @param radius
	 * @param time
	 * @return the clusters
	 */
	public ArrayList<Cluster> findClusters(List<ClusterPoint> points, double radius, int time)
	{
		if (clusteringAlgorithm == ClusteringAlgorithm.ParticleLinkage)
			return runParticleLinkage(points, radius);

		// Get the density around each point. Points with no density cannot be clustered.
		// Ensure that we only ignore points that could never be within radius of the centroid 
		// of any other pair:
		//
		// Set point C on the circle drawn with radius R around both points A and B. Distance A-B > R.
		// When the angle ACB is <90 then the line C-B intersects A's circle. If ACB = 90 then the 
		// line C-B cannot intersect A's circle. Distance AB = sqrt(2R^2) = sqrt(2) * R
		//
		//        -------   ------   
		//       /       \ /       \
		//      /         C         \
		//     /         / \         \
		//     |     A   | |   B     |
		//     \         \ /         /
		//      \         X         /
		//       \       / \       /
		//        -------   -------

		int[] density = calculateDensity(points, 1.4142 * radius);

		// Extract initial cluster points using molecules with a density above 1
		// (All other points cannot be clustered at this radius)
		ArrayList<Cluster> candidates = new ArrayList<Cluster>(density.length);
		ArrayList<Cluster> singles = new ArrayList<Cluster>(density.length);
		int i = 0;
		for (ClusterPoint p : points)
		{
			final Cluster c = new Cluster(p);
			if (density[i] > 0)
				candidates.add(c);
			else
				singles.add(c);
			i++;
		}

		if (candidates.isEmpty())
			return singles;

		// Check for time information if required
		switch (clusteringAlgorithm)
		{
			case ClosestDistancePriority:
			case ClosestTimePriority:
				if (noTimeInformation(candidates))
				{
					tracker.log("No time information among candidates");
					clusteringAlgorithm = ClusteringAlgorithm.Closest;
				}

				// All other methods do not use time information	
			default:
				break;
		}

		tracker.log("Starting clustering : %d singles, %d cluster candidates", singles.size(), candidates.size());
		tracker.log("Algorithm = %s", clusteringAlgorithm.toString());

		// Find the bounds of the candidates
		double minx = candidates.get(0).x;
		double miny = candidates.get(0).y;
		double maxx = minx, maxy = miny;
		for (Cluster c : candidates)
		{
			if (minx > c.x)
				minx = c.x;
			else if (maxx < c.x)
				maxx = c.x;
			if (miny > c.y)
				miny = c.y;
			else if (maxy < c.y)
				maxy = c.y;
		}

		// Assign to a grid
		final int maxBins = 500;
		// If tracking potential neighbours then the cells must be larger to cover the increased distance
		final double cellSize = (clusteringAlgorithm == ClusteringAlgorithm.Pairwise) ? radius : radius * 1.4142;
		final double xBinWidth = Math.max(cellSize, (maxx - minx) / maxBins);
		final double yBinWidth = Math.max(cellSize, (maxy - miny) / maxBins);
		final int nXBins = 1 + (int) ((maxx - minx) / xBinWidth);
		final int nYBins = 1 + (int) ((maxy - miny) / yBinWidth);
		Cluster[][] grid = new Cluster[nXBins][nYBins];
		for (Cluster c : candidates)
		{
			final int xBin = (int) ((c.x - minx) / xBinWidth);
			final int yBin = (int) ((c.y - miny) / yBinWidth);
			// Build a single linked list
			c.xBin = xBin;
			c.yBin = yBin;
			c.next = grid[xBin][yBin];
			grid[xBin][yBin] = c;
		}

		final double r2 = radius * radius;

		tracker.log("Clustering " + clusteringAlgorithm.toString() + " ...");
		ArrayList<Cluster> clusters;
		switch (clusteringAlgorithm)
		{
			case Pairwise:
				clusters = runPairwise(grid, nXBins, nYBins, r2, minx, miny, xBinWidth, yBinWidth, candidates, singles);
				break;

			case PairwiseWithoutNeighbours:
				clusters = runPairwiseWithoutNeighbours(grid, nXBins, nYBins, r2, minx, miny, xBinWidth, yBinWidth,
						candidates, singles);
				break;

			case ClosestTimePriority:
				clusters = runClosestTimePriority(grid, nXBins, nYBins, r2, time, minx, miny, xBinWidth, yBinWidth,
						candidates, singles);
				break;

			case ClosestDistancePriority:
				clusters = runClosestDistancePriority(grid, nXBins, nYBins, r2, time, minx, miny, xBinWidth, yBinWidth,
						candidates, singles);
				break;

			default:
				clusters = runClosest(grid, nXBins, nYBins, r2, minx, miny, xBinWidth, yBinWidth, candidates, singles);
		}

		tracker.progress(1);

		tracker.log("Found %d clusters", (clusters == null) ? 0 : clusters.size());

		return clusters;
	}

	/**
	 * Join the closest unlinked particle to its neighbour particle/cluster
	 * 
	 * @param points
	 * @param radius
	 * @return The clusters
	 */
	private ArrayList<Cluster> runParticleLinkage(List<ClusterPoint> points, double radius)
	{
		int[] density = calculateDensity(points, radius);

		// Extract initial cluster points using molecules with a density above 1
		// (All other points cannot be clustered at this radius)
		ArrayList<ExtendedClusterPoint> candidates = new ArrayList<ExtendedClusterPoint>(density.length);
		ArrayList<Cluster> singles = new ArrayList<Cluster>(density.length);
		int i = 0, id = 0;
		for (ClusterPoint p : points)
		{
			if (density[i] > 0)
				// Store the point using the next pointer of a new point which will be used for clustering
				candidates.add(new ExtendedClusterPoint(id++, p.x, p.y, 0, p));
			else
				singles.add(new Cluster(p));
			i++;
		}

		if (candidates.isEmpty())
			return singles;

		tracker.log("Starting clustering : %d singles, %d cluster candidates", singles.size(), candidates.size());
		tracker.log("Algorithm = %s", clusteringAlgorithm.toString());

		// Find the bounds of the candidates
		double minx = candidates.get(0).x;
		double miny = candidates.get(0).y;
		double maxx = minx, maxy = miny;
		for (ExtendedClusterPoint c : candidates)
		{
			if (minx > c.x)
				minx = c.x;
			else if (maxx < c.x)
				maxx = c.x;
			if (miny > c.y)
				miny = c.y;
			else if (maxy < c.y)
				maxy = c.y;
		}

		// Assign to a grid
		final int maxBins = 500;
		final double cellSize = radius * 1.01; // Add an error margin
		final double xBinWidth = Math.max(cellSize, (maxx - minx) / maxBins);
		final double yBinWidth = Math.max(cellSize, (maxy - miny) / maxBins);
		final int nXBins = 1 + (int) ((maxx - minx) / xBinWidth);
		final int nYBins = 1 + (int) ((maxy - miny) / yBinWidth);
		ExtendedClusterPoint[][] grid = new ExtendedClusterPoint[nXBins][nYBins];
		for (ExtendedClusterPoint c : candidates)
		{
			final int xBin = (int) ((c.x - minx) / xBinWidth);
			final int yBin = (int) ((c.y - miny) / yBinWidth);
			// Build a single linked list
			c.nextE = grid[xBin][yBin];
			grid[xBin][yBin] = c;
		}

		final double r2 = radius * radius;

		tracker.log("Clustering " + clusteringAlgorithm.toString() + " ...");

		ArrayList<Cluster> clusters = runParticleLinkage(grid, nXBins, nYBins, r2, minx, miny, xBinWidth, yBinWidth,
				candidates, singles);

		tracker.progress(1);

		tracker.log("Found %d clusters", (clusters == null) ? 0 : clusters.size());

		return clusters;
	}

	private ArrayList<Cluster> runParticleLinkage(ExtendedClusterPoint[][] grid, int nXBins, int nYBins, double r2,
			double minx, double miny, double xBinWidth, double yBinWidth, ArrayList<ExtendedClusterPoint> candidates,
			ArrayList<Cluster> singles)
	{
		int N = candidates.size();
		int candidatesProcessed = 0;
		nextClusterId = 0; // Incremented within joinClosestParticle(...)

		// Used to store the cluster for each candidate
		int[] clusterId = new int[candidates.size()];

		int nProcessed = 0;
		while ((nProcessed = joinClosestParticle(grid, nXBins, nYBins, r2, minx, miny, xBinWidth, yBinWidth, clusterId)) > 0)
		{
			if (tracker.stop())
				return null;
			candidatesProcessed += nProcessed;
			tracker.progress(candidatesProcessed, N);
		}

		tracker.log("Processed %d / %d", candidatesProcessed, N);
		tracker.log("%d candidates linked into %d clusters", candidates.size(), nextClusterId);

		// Create clusters from the original cluster points using the assignments
		Cluster[] clusters = new Cluster[nextClusterId];
		int failed = singles.size();
		for (int i = 0; i < clusterId.length; i++)
		{
			ClusterPoint originalPoint = candidates.get(i).next;
			if (clusterId[i] == 0)
			{
				//throw new RuntimeException("Failed to assign a cluster to a candidate particle");
				//tracker.log("Failed to assign a cluster to a candidate particle: " + i);
				singles.add(new Cluster(originalPoint));

				// Check is a neighbour
				boolean neighbour = false;
				for (int j = 0; j < clusterId.length; j++)
				{
					if (i == j)
						continue;
					ClusterPoint p = candidates.get(j);
					if (originalPoint.distance2(p) < r2)
					{
						neighbour = true;
						break;
					}
				}
				if (neighbour)
					tracker.log("Failed to assign a cluster to a candidate particle: " + i);
			}
			else
			{
				// The next points were used to store the original cluster points
				int clusterIndex = clusterId[i] - 1;
				if (clusters[clusterIndex] == null)
					clusters[clusterIndex] = new Cluster(originalPoint);
				else
					clusters[clusterIndex].add(originalPoint);
			}
		}
		failed = singles.size() - failed;
		tracker.log("Failed to assign %d candidates", failed);

		for (int i = 0; i < clusters.length; i++)
		{
			if (clusters[i] != null)
				singles.add(clusters[i]);
			else
				; //tracker.log("Empty cluster ID %d", i);
		}
		return singles;
	}

	/**
	 * Search for the closest pair of particles, one of which is not in a cluster, below the squared radius distance and
	 * join them
	 * 
	 * @param grid
	 * @param nXBins
	 * @param nYBins
	 * @param r2
	 *            The squared radius distance
	 * @param yBinWidth
	 * @param xBinWidth
	 * @param miny
	 * @param minx
	 * @param clusterId
	 * @return The number of points assigned to clusters (either 0, 1, or 2)
	 */
	private int joinClosestParticle(ExtendedClusterPoint[][] grid, final int nXBins, final int nYBins, final double r2,
			double minx, double miny, double xBinWidth, double yBinWidth, int[] clusterId)
	{
		double min = r2;
		ExtendedClusterPoint pair1 = null, pair2 = null;
		for (int yBin = 0; yBin < nYBins; yBin++)
		{
			for (int xBin = 0; xBin < nXBins; xBin++)
			{
				for (ExtendedClusterPoint c1 = grid[xBin][yBin]; c1 != null; c1 = c1.nextE)
				{
					final boolean cluster1 = c1.inCluster;

					// Build a list of which cells to compare up to a maximum of 4
					//      | 0,0  |  1,0
					// ------------+-----
					// -1,1 | 0,1  |  1,1

					// Compare to neighbours and find the closest.
					// Use either the radius threshold or the current closest distance
					// which may have been set by an earlier comparison.
					ExtendedClusterPoint other = null;

					for (ExtendedClusterPoint c2 = c1.nextE; c2 != null; c2 = c2.nextE)
					{
						// Ignore comparing points that are both in a cluster
						if (cluster1 && c2.inCluster)
							continue;

						final double d2 = c1.distance2(c2);
						if (d2 < min)
						{
							min = d2;
							other = c2;
						}
					}

					if (yBin < nYBins - 1)
					{
						for (ExtendedClusterPoint c2 = grid[xBin][yBin + 1]; c2 != null; c2 = c2.nextE)
						{
							// Ignore comparing points that are both in a cluster
							if (cluster1 && c2.inCluster)
								continue;
							final double d2 = c1.distance2(c2);
							if (d2 < min)
							{
								min = d2;
								other = c2;
							}
						}
						if (xBin > 0)
						{
							for (ExtendedClusterPoint c2 = grid[xBin - 1][yBin + 1]; c2 != null; c2 = c2.nextE)
							{
								// Ignore comparing points that are both in a cluster
								if (cluster1 && c2.inCluster)
									continue;
								final double d2 = c1.distance2(c2);
								if (d2 < min)
								{
									min = d2;
									other = c2;
								}
							}
						}
					}
					if (xBin < nXBins - 1)
					{
						for (ExtendedClusterPoint c2 = grid[xBin + 1][yBin]; c2 != null; c2 = c2.nextE)
						{
							// Ignore comparing points that are both in a cluster
							if (cluster1 && c2.inCluster)
								continue;
							final double d2 = c1.distance2(c2);
							if (d2 < min)
							{
								min = d2;
								other = c2;
							}
						}
						if (yBin < nYBins - 1)
						{
							for (ExtendedClusterPoint c2 = grid[xBin + 1][yBin + 1]; c2 != null; c2 = c2.nextE)
							{
								// Ignore comparing points that are both in a cluster
								if (cluster1 && c2.inCluster)
									continue;
								final double d2 = c1.distance2(c2);
								if (d2 < min)
								{
									min = d2;
									other = c2;
								}
							}
						}
					}

					// Store the details of the closest pair
					if (other != null)
					{
						pair1 = c1;
						pair2 = other;
					}
				}
			}
		}

		// Assign the closest pair.
		if (pair1 != null)
		{
			int nProcessed = 1;

			// Create a new cluster if necessary
			if (clusterId[pair2.id] == 0)
			{
				nProcessed = 2;
				clusterId[pair2.id] = ++nextClusterId;
				pair2.inCluster = true;
			}

			clusterId[pair1.id] = clusterId[pair2.id];
			pair1.inCluster = true;

			return nProcessed;
		}

		// This should not happen if the density manager counted neighbours correctly
		// (i.e. all points should have at least one neighbour)		
		return 0;
	}

	/**
	 * Count the number of points around each point.
	 * 
	 * @param points
	 * @param radius
	 * @return The density count
	 */
	private int[] calculateDensity(List<ClusterPoint> points, double radius)
	{
		float[] xcoord = new float[points.size()];
		float[] ycoord = new float[points.size()];
		int i = 0;
		for (ClusterPoint p : points)
		{
			xcoord[i] = (float) p.x;
			ycoord[i] = (float) p.y;
			i++;
		}
		// The bounds are not used in the density calculation when not adjusting for borders
		Rectangle bounds = new Rectangle();
		DensityManager dm = new DensityManager(xcoord, ycoord, bounds);
		return dm.calculateDensity((float) radius, false);
	}

	/**
	 * Check there are at least two different time points in the data
	 * 
	 * @param candidates
	 * @return true if there are no different time points
	 */
	private boolean noTimeInformation(ArrayList<Cluster> candidates)
	{
		final int firstT = candidates.get(0).head.t;
		for (Cluster c : candidates)
			if (firstT != c.head.t)
				return false;
		return true;
	}

	/**
	 * Sweep the all-vs-all clusters and make potential links between clusters.
	 * If a link can be made to a closer cluster then break the link and rejoin.
	 * Then join all the links into clusters.
	 * 
	 * @param grid
	 * @param nXBins
	 * @param nYBins
	 * @param r2
	 * @param minx
	 * @param miny
	 * @param xBinWidth
	 * @param yBinWidth
	 * @param candidates
	 * @param singles
	 * @return The clusters
	 */
	private ArrayList<Cluster> runPairwise(Cluster[][] grid, final int nXBins, final int nYBins, final double r2,
			final double minx, final double miny, final double xBinWidth, final double yBinWidth,
			ArrayList<Cluster> candidates, ArrayList<Cluster> singles)
	{
		while (findLinks(grid, nXBins, nYBins, r2))
		{
			if (tracker.stop())
				return null;

			joinLinks(grid, nXBins, nYBins, candidates);

			// Reassign the grid
			for (Cluster c : candidates)
			{
				final int xBin = (int) ((c.x - minx) / xBinWidth);
				final int yBin = (int) ((c.y - miny) / yBinWidth);
				// Build a single linked list
				c.next = grid[xBin][yBin];
				grid[xBin][yBin] = c;
			}
		}

		candidates.addAll(singles);
		return candidates;
	}

	/**
	 * Search for potential links between clusters that are below the squared radius distance. Store if the clusters
	 * have any neighbours within 2*r^2.
	 * 
	 * @param grid
	 * @param nXBins
	 * @param nYBins
	 * @param r2
	 *            The squared radius distance
	 * @return True if any links were made
	 */
	private boolean findLinks(Cluster[][] grid, final int nXBins, final int nYBins, final double r2)
	{
		Cluster[] neighbours = new Cluster[5];
		boolean linked = false;
		for (int yBin = 0; yBin < nYBins; yBin++)
		{
			for (int xBin = 0; xBin < nXBins; xBin++)
			{
				for (Cluster c1 = grid[xBin][yBin]; c1 != null; c1 = c1.next)
				{
					// Build a list of which cells to compare up to a maximum of 5
					//      | 0,0  |  1,0
					// ------------+-----
					// -1,1 | 0,1  |  1,1

					int count = 0;
					neighbours[count++] = c1.next;

					if (yBin < nYBins - 1)
					{
						neighbours[count++] = grid[xBin][yBin + 1];
						if (xBin > 0)
							neighbours[count++] = grid[xBin - 1][yBin + 1];
					}
					if (xBin < nXBins - 1)
					{
						neighbours[count++] = grid[xBin + 1][yBin];
						if (yBin < nYBins - 1)
							neighbours[count++] = grid[xBin + 1][yBin + 1];
					}

					// Compare to neighbours and find the closest.
					// Use either the radius threshold or the current closest distance
					// which may have been set by an earlier comparison.
					double min = (c1.closest == null) ? r2 : c1.d2;
					Cluster other = null;
					while (count-- > 0)
					{
						for (Cluster c2 = neighbours[count]; c2 != null; c2 = c2.next)
						{
							final double d2 = c1.distance2(c2);
							if (d2 < min)
							{
								min = d2;
								other = c2;
							}
						}
					}

					if (other != null)
					{
						// Store the potential link between the two clusters
						c1.link(other, min);
						linked = true;
					}
				}
			}
		}
		return linked;
	}

	/**
	 * Join valid links between clusters. Resets the link candidates.
	 * 
	 * @param grid
	 * @param nXBins
	 * @param nYBins
	 * @param candidates
	 *            Re-populate will all the remaining clusters
	 */
	private void joinLinks(Cluster[][] grid, int nXBins, int nYBins, ArrayList<Cluster> candidates)
	{
		candidates.clear();

		for (int yBin = 0; yBin < nYBins; yBin++)
		{
			for (int xBin = 0; xBin < nXBins; xBin++)
			{
				for (Cluster c1 = grid[xBin][yBin]; c1 != null; c1 = c1.next)
				{
					if (c1.validLink())
					{
						c1.add(c1.closest);
					}
					// Reset the link candidates
					c1.closest = null;

					// Store all remaining clusters
					if (c1.n != 0)
					{
						candidates.add(c1);
					}
				}

				// Reset the grid
				grid[xBin][yBin] = null;
			}
		}
	}

	/**
	 * Sweep the all-vs-all clusters and make potential links between clusters.
	 * If a link can be made to a closer cluster then break the link and rejoin.
	 * Then join all the links into clusters only if the pair has no other neighbours. Default to joining the closest
	 * pair.
	 * 
	 * @param grid
	 * @param nXBins
	 * @param nYBins
	 * @param r2
	 * @param minx
	 * @param miny
	 * @param xBinWidth
	 * @param yBinWidth
	 * @param candidates
	 * @param singles
	 * @return
	 */
	private ArrayList<Cluster> runPairwiseWithoutNeighbours(Cluster[][] grid, final int nXBins, final int nYBins,
			final double r2, final double minx, final double miny, final double xBinWidth, final double yBinWidth,
			ArrayList<Cluster> candidates, ArrayList<Cluster> singles)
	{
		// Store if the clusters have any neighbours within sqrt(2)*r and remove them 
		// from the next loop.
		int N = candidates.size();
		ArrayList<Cluster> joined = new ArrayList<Cluster>();
		while (findLinksAndCountNeighbours(grid, nXBins, nYBins, r2, singles))
		{
			if (tracker.stop())
				return null;

			int joins = joinLinks(grid, nXBins, nYBins, r2, candidates, joined);
			if (joins == 0)
				break; // This should not happen

			tracker.progress(N - candidates.size(), N);

			// TODO - determine at what point it is faster to reassign the grid
			if (joins < candidates.size() / 5)
			{
				// Reassigning the whole grid is a costly step when the number of joins is small.
				// In that case check the clusters that were updated and reassign them to a new
				// grid position if necessary.
				for (Cluster c : candidates)
				{
					if (c.neighbour != 0)
					{
						c.neighbour = 0;
						final int xBin = (int) ((c.x - minx) / xBinWidth);
						final int yBin = (int) ((c.y - miny) / yBinWidth);

						// Check the current grid position.
						if (xBin != c.xBin || yBin != c.yBin)
						{
							remove(grid, c);

							c.xBin = xBin;
							c.yBin = yBin;
							c.next = grid[xBin][yBin];
							grid[xBin][yBin] = c;
						}
					}
				}

				// We must remove the joined clusters from the grid
				for (Cluster c : joined)
				{
					remove(grid, c);
				}
			}
			else
			{
				// Reassign the grid.
				for (int xBin = 0; xBin < nXBins; xBin++)
					for (int yBin = 0; yBin < nYBins; yBin++)
						grid[xBin][yBin] = null;
				for (Cluster c : candidates)
				{
					// Only candidates that have been flagged as a join may have changed 
					// their grid position
					if (c.neighbour != 0)
					{
						c.xBin = (int) ((c.x - minx) / xBinWidth);
						c.yBin = (int) ((c.y - miny) / yBinWidth);
					}
					// Build a single linked list
					c.next = grid[c.xBin][c.yBin];
					grid[c.xBin][c.yBin] = c;
					c.neighbour = 0;
				}
			}
		}

		return combine(singles, grid, nXBins, nYBins);
	}

	/**
	 * Search for potential links between clusters that are below the squared radius distance. Store if the clusters
	 * have any neighbours within 2*r^2.
	 * 
	 * @param grid
	 * @param nXBins
	 * @param nYBins
	 * @param r2
	 *            The squared radius distance
	 * @param singles
	 *            Add remaining clusters that have no neighbours
	 * @return True if any links were made
	 */
	private boolean findLinksAndCountNeighbours(Cluster[][] grid, final int nXBins, final int nYBins, final double r2,
			ArrayList<Cluster> singles)
	{
		Cluster[] neighbours = new Cluster[5];
		boolean linked = false;
		final double neighbourDistance = 2 * r2;
		for (int yBin = 0; yBin < nYBins; yBin++)
		{
			for (int xBin = 0; xBin < nXBins; xBin++)
			{
				Cluster previous = null;
				for (Cluster c1 = grid[xBin][yBin]; c1 != null; c1 = c1.next)
				{
					// Build a list of which cells to compare up to a maximum of 5
					//      | 0,0  |  1,0
					// ------------+-----
					// -1,1 | 0,1  |  1,1

					int count = 0;
					neighbours[count++] = c1.next;

					if (yBin < nYBins - 1)
					{
						neighbours[count++] = grid[xBin][yBin + 1];
						if (xBin > 0)
							neighbours[count++] = grid[xBin - 1][yBin + 1];
					}
					if (xBin < nXBins - 1)
					{
						neighbours[count++] = grid[xBin + 1][yBin];
						if (yBin < nYBins - 1)
							neighbours[count++] = grid[xBin + 1][yBin + 1];
					}

					// Compare to neighbours and find the closest.
					// Use either the radius threshold or the current closest distance
					// which may have been set by an earlier comparison.
					double min = (c1.closest == null) ? r2 : c1.d2;
					Cluster other = null;
					while (count-- > 0)
					{
						for (Cluster c2 = neighbours[count]; c2 != null; c2 = c2.next)
						{
							final double d2 = c1.distance2(c2);
							if (d2 < neighbourDistance)
							{
								c1.neighbour++;
								c2.neighbour++;
								if (d2 < min)
								{
									min = d2;
									other = c2;
								}
							}
						}
					}

					if (other != null)
					{
						// Store the potential link between the two clusters
						c1.link(other, min);
						linked = true;
					}

					// Check for singles
					if (c1.neighbour == 0)
					{
						// Add singles to the singles list and remove from the grid
						singles.add(c1);
						if (previous == null)
							grid[xBin][yBin] = c1.next;
						else
							previous.next = c1.next;
					}
					else
					{
						previous = c1;
					}
				}
			}
		}
		return linked;
	}

	/**
	 * Join valid links between clusters. Resets the link candidates. Use the neighbour count property to flag if the
	 * candidate was joined to another cluster.
	 * 
	 * @param grid
	 * @param nXBins
	 * @param nYBins
	 * @param r2
	 * @param candidates
	 *            Re-populate will all the remaining clusters
	 * @return The number of joins that were made
	 */
	private int joinLinks(Cluster[][] grid, int nXBins, int nYBins, double r2, ArrayList<Cluster> candidates,
			ArrayList<Cluster> joined)
	{
		candidates.clear();
		joined.clear();

		double min = r2;
		Cluster cluster1 = null, cluster2 = null;
		for (int yBin = 0; yBin < nYBins; yBin++)
		{
			for (int xBin = 0; xBin < nXBins; xBin++)
			{
				for (Cluster c1 = grid[xBin][yBin]; c1 != null; c1 = c1.next)
				{
					int joinFlag = 0;
					if (c1.validLink())
					{
						// Check if each cluster only has 1 neighbour
						if (c1.neighbour == 1 && c1.closest.neighbour == 1)
						{
							//System.out.printf("Joining pairs with no neighbours @ %f\n", Math.sqrt(c1.d2));
							c1.add(c1.closest);
							joined.add(c1.closest);
							joinFlag = 1;
						}
						else if (c1.d2 < min)
						{
							// Otherwise find the closest pair in case no joins are made
							min = c1.d2;
							cluster1 = c1;
							cluster2 = c1.closest;
						}
					}
					// Reset the link candidates
					c1.closest = null;

					// Store all remaining clusters
					if (c1.n != 0)
					{
						candidates.add(c1);
						c1.neighbour = joinFlag;
					}
				}
			}
		}

		//// If no joins were made then join the closest pair
		//if (joined.isEmpty() && cluster1 != null)

		// Always join the closest pair if it has not been already
		if (cluster1 != null && cluster1.neighbour == 0)
		{
			//System.out.printf("Joining closest pair @ %f\n", Math.sqrt(cluster1.d2));
			cluster1.add(cluster2);
			// Remove cluster 2 from the candidates
			candidates.remove(cluster2);
			joined.add(cluster2);
			cluster1.neighbour = 1;
		}

		return joined.size();
	}

	/**
	 * The process should iterate finding the closest nodes, joining them and repeating.
	 * The iterative process of joining the closest pair will be slow. Hopefully it will be manageable.
	 * 
	 * @param grid
	 * @param nXBins
	 * @param nYBins
	 * @param r2
	 * @param minx
	 * @param miny
	 * @param xBinWidth
	 * @param yBinWidth
	 * @param candidates
	 * @param singles
	 * @return
	 */
	private ArrayList<Cluster> runClosest(Cluster[][] grid, int nXBins, int nYBins, double r2, double minx,
			double miny, double xBinWidth, double yBinWidth, ArrayList<Cluster> candidates, ArrayList<Cluster> singles)
	{
		int N = candidates.size();
		int candidatesProcessed = 0;
		int s = singles.size();
		final boolean trackProgress = (tracker.getClass() != NullTrackProgress.class);
		while (joinClosest(grid, nXBins, nYBins, r2, minx, miny, xBinWidth, yBinWidth, singles))
		{
			if (tracker.stop())
				return null;

			// The number of candidates that have been processed is incremented by the number of singles
			if (trackProgress)
			{
				candidatesProcessed += singles.size() - s;
				s = singles.size();
				tracker.progress(candidatesProcessed++, N);
			}
		}

		return combine(singles, grid, nXBins, nYBins);
	}

	/**
	 * Search for the closest pair of clusters that are below the squared radius distance and join them
	 * 
	 * @param grid
	 * @param nXBins
	 * @param nYBins
	 * @param r2
	 *            The squared radius distance
	 * @param yBinWidth
	 * @param xBinWidth
	 * @param miny
	 * @param minx
	 * @param singles
	 *            Add remaining clusters that have no neighbours
	 * @return True if a join was made
	 */
	private boolean joinClosest(Cluster[][] grid, final int nXBins, final int nYBins, final double r2, double minx,
			double miny, double xBinWidth, double yBinWidth, ArrayList<Cluster> singles)
	{
		double min = r2;
		Cluster pair1 = null, pair2 = null;
		final double neighbourDistance = 2 * r2;
		for (int yBin = 0; yBin < nYBins; yBin++)
		{
			for (int xBin = 0; xBin < nXBins; xBin++)
			{
				Cluster previous = null;
				for (Cluster c1 = grid[xBin][yBin]; c1 != null; c1 = c1.next)
				{
					// Build a list of which cells to compare up to a maximum of 4
					//      | 0,0  |  1,0
					// ------------+-----
					// -1,1 | 0,1  |  1,1

					// Compare to neighbours and find the closest.
					// Use either the radius threshold or the current closest distance
					// which may have been set by an earlier comparison.
					Cluster other = null;

					for (Cluster c2 = c1.next; c2 != null; c2 = c2.next)
					{
						final double d2 = c1.distance2(c2);
						if (d2 < min)
						{
							min = d2;
							other = c2;
							c1.neighbour++;
							c2.neighbour++;
						}
						else if (d2 < neighbourDistance)
						{
							c1.neighbour++;
							c2.neighbour++;
						}
					}

					if (yBin < nYBins - 1)
					{
						for (Cluster c2 = grid[xBin][yBin + 1]; c2 != null; c2 = c2.next)
						{
							final double d2 = c1.distance2(c2);
							if (d2 < min)
							{
								min = d2;
								other = c2;
								c1.neighbour++;
								c2.neighbour++;
							}
							else if (d2 < neighbourDistance)
							{
								c1.neighbour++;
								c2.neighbour++;
							}
						}
						if (xBin > 0)
						{
							for (Cluster c2 = grid[xBin - 1][yBin + 1]; c2 != null; c2 = c2.next)
							{
								final double d2 = c1.distance2(c2);
								if (d2 < min)
								{
									min = d2;
									other = c2;
									c1.neighbour++;
									c2.neighbour++;
								}
								else if (d2 < neighbourDistance)
								{
									c1.neighbour++;
									c2.neighbour++;
								}
							}
						}
					}
					if (xBin < nXBins - 1)
					{
						for (Cluster c2 = grid[xBin + 1][yBin]; c2 != null; c2 = c2.next)
						{
							final double d2 = c1.distance2(c2);
							if (d2 < min)
							{
								min = d2;
								other = c2;
								c1.neighbour++;
								c2.neighbour++;
							}
							else if (d2 < neighbourDistance)
							{
								c1.neighbour++;
								c2.neighbour++;
							}
						}
						if (yBin < nYBins - 1)
						{
							for (Cluster c2 = grid[xBin + 1][yBin + 1]; c2 != null; c2 = c2.next)
							{
								final double d2 = c1.distance2(c2);
								if (d2 < min)
								{
									min = d2;
									other = c2;
									c1.neighbour++;
									c2.neighbour++;
								}
								else if (d2 < neighbourDistance)
								{
									c1.neighbour++;
									c2.neighbour++;
								}
							}
						}
					}

					// Store the details of the closest pair
					if (other != null)
					{
						pair1 = c1;
						pair2 = other;
					}

					// Check for singles
					if (c1.neighbour == 0)
					{
						// Add singles to the singles list and remove from the grid
						singles.add(c1);
						if (previous == null)
							grid[xBin][yBin] = c1.next;
						else
							previous.next = c1.next;
					}
					else
					{
						previous = c1;
						c1.neighbour = 0;
					}
				}
			}
		}

		// Join the closest pair
		if (pair1 != null)
		{
			pair2.add(pair1);

			remove(grid, pair1);

			// Reassign the updated grid position
			final int xBin = (int) ((pair2.x - minx) / xBinWidth);
			final int yBin = (int) ((pair2.y - miny) / yBinWidth);

			if (xBin != pair2.xBin || yBin != pair2.yBin)
			{
				remove(grid, pair2);

				// Build a single linked list
				pair2.xBin = xBin;
				pair2.yBin = yBin;
				pair2.next = grid[xBin][yBin];
				grid[xBin][yBin] = pair2;
			}

			return true;
		}

		return false;
	}

	/**
	 * Remove cluster from the grid by sweeping the linked list grid position
	 * 
	 * @param grid
	 * @param cluster
	 */
	private void remove(Cluster[][] grid, Cluster cluster)
	{
		Cluster previous = null;
		for (Cluster c1 = grid[cluster.xBin][cluster.yBin]; c1 != null; c1 = c1.next)
		{
			if (c1 == cluster)
			{
				if (previous == null)
					grid[cluster.xBin][cluster.yBin] = c1.next;
				else
					previous.next = c1.next;
				return;
			}
			else
				previous = c1;
		}
	}

	/**
	 * Add the clusters in the grid to the existing singles
	 * 
	 * @param singles
	 * @param grid
	 * @param nXBins
	 * @param nYBins
	 * @return The combined clusters
	 */
	private ArrayList<Cluster> combine(ArrayList<Cluster> singles, Cluster[][] grid, int nXBins, int nYBins)
	{
		for (int xBin = 0; xBin < nXBins; xBin++)
		{
			for (int yBin = 0; yBin < nYBins; yBin++)
			{
				for (Cluster c1 = grid[xBin][yBin]; c1 != null; c1 = c1.next)
				{
					//if (c1.n > 0)
					singles.add(c1);
				}
			}
		}
		return singles;
	}

	private ArrayList<Cluster> runClosestTimePriority(Cluster[][] grid, int nXBins, int nYBins, double r2, int time,
			double minx, double miny, double xBinWidth, double yBinWidth, ArrayList<Cluster> candidates,
			ArrayList<Cluster> singles)
	{
		int N = candidates.size();
		int candidatesProcessed = 0;
		int s = singles.size();
		final boolean trackProgress = (tracker.getClass() != NullTrackProgress.class);
		TimeCluster[][] newGrid = convertGrid(grid, nXBins, nYBins);
		while (joinClosestTimePriority(newGrid, nXBins, nYBins, r2, time, minx, miny, xBinWidth, yBinWidth, singles))
		{
			if (tracker.stop())
				return null;

			// The number of candidates that have been processed is incremented by the number of singles
			if (trackProgress)
			{
				candidatesProcessed += singles.size() - s;
				s = singles.size();
				tracker.progress(candidatesProcessed++, N);
			}
		}

		return combine(singles, newGrid, nXBins, nYBins);
	}

	private TimeCluster[][] convertGrid(Cluster[][] grid, int nXBins, int nYBins)
	{
		TimeCluster[][] newGrid = new TimeCluster[nXBins][nYBins];
		for (int yBin = 0; yBin < nYBins; yBin++)
		{
			for (int xBin = 0; xBin < nXBins; xBin++)
			{
				for (Cluster c1 = grid[xBin][yBin]; c1 != null; c1 = c1.next)
				{
					// Build a single linked list
					TimeCluster c = new TimeCluster(c1.head);
					c.pulse = getPulse(c.start);
					c.xBin = xBin;
					c.yBin = yBin;
					c.next = newGrid[xBin][yBin];
					newGrid[xBin][yBin] = c;
				}
			}
		}
		return newGrid;
	}

	/**
	 * Search for the closest pair of clusters that are below the squared radius distance and join them. Join closest in
	 * time and then distance.
	 * 
	 * @param grid
	 * @param nXBins
	 * @param nYBins
	 * @param r2
	 *            The squared radius distance
	 * @param time
	 *            The time threshold
	 * @param yBinWidth
	 * @param xBinWidth
	 * @param miny
	 * @param minx
	 * @param singles
	 *            Add remaining clusters that have no neighbours
	 * @return True if a join was made
	 */
	private boolean joinClosestTimePriority(TimeCluster[][] grid, final int nXBins, final int nYBins, final double r2,
			int time, double minx, double miny, double xBinWidth, double yBinWidth, ArrayList<Cluster> singles)
	{
		double minD = Double.POSITIVE_INFINITY;
		int minT = Integer.MAX_VALUE;
		TimeCluster pair1 = null, pair2 = null;
		final double neighbourDistance = 2 * r2;
		TimeCluster[] neighbourCells = new TimeCluster[5];
		final boolean checkPulseInterval = pulseInterval > 0;
		for (int yBin = 0; yBin < nYBins; yBin++)
		{
			for (int xBin = 0; xBin < nXBins; xBin++)
			{
				TimeCluster previous = null;
				for (TimeCluster c1 = grid[xBin][yBin]; c1 != null; c1 = (TimeCluster) c1.next)
				{
					// Build a list of which cells to compare up to a maximum of 4
					//      | 0,0  |  1,0
					// ------------+-----
					// -1,1 | 0,1  |  1,1

					// Compare to neighbours and find the closest.
					// Use either the radius threshold or the current closest distance
					// which may have been set by an earlier comparison.
					TimeCluster other = null;
					int cells = 1;
					neighbourCells[0] = (TimeCluster) c1.next;
					if (yBin < nYBins - 1)
					{
						neighbourCells[cells++] = grid[xBin][yBin + 1];
						if (xBin > 0)
							neighbourCells[cells++] = grid[xBin - 1][yBin + 1];
					}
					if (xBin < nXBins - 1)
					{
						neighbourCells[cells++] = grid[xBin + 1][yBin];
						if (yBin < nYBins - 1)
							neighbourCells[cells++] = grid[xBin + 1][yBin + 1];
					}

					for (int c = 0; c < cells; c++)
					{
						for (TimeCluster c2 = neighbourCells[c]; c2 != null; c2 = (TimeCluster) c2.next)
						{
							if (checkPulseInterval && c1.pulse != c2.pulse)
								continue;

							final int gap = c1.gap(c2);
							if (gap <= time)
							{
								final double d2 = c1.distance2(c2);

								if (d2 < neighbourDistance)
								{
									// Check if the two clusters can be merged
									if (!c1.validUnion(c2))
										continue;

									// Count any possible neighbour
									c1.neighbour++;
									c2.neighbour++;
								}
								if (d2 <= r2)
								{
									// This is within the time and distance thresholds.									
									// Find closest pair with time priority
									if ((gap < minT) || (gap <= minT && d2 < minD))
									{
										minD = d2;
										minT = gap;
										other = c2;
									}
								}
							}
						}
					}

					// Store the details of the closest pair
					if (other != null)
					{
						pair1 = c1;
						pair2 = other;
					}

					// Check for singles
					if (c1.neighbour == 0)
					{
						// Add singles to the singles list and remove from the grid
						singles.add(c1);
						if (previous == null)
							grid[xBin][yBin] = (TimeCluster) c1.next;
						else
							previous.next = c1.next;
					}
					else
					{
						previous = c1;
						c1.neighbour = 0;
					}
				}
			}
		}

		// Join the closest pair
		if (pair1 != null)
		{
			pair2.add(pair1);

			remove(grid, pair1);

			// Reassign the updated grid position
			final int xBin = (int) ((pair2.x - minx) / xBinWidth);
			final int yBin = (int) ((pair2.y - miny) / yBinWidth);

			if (xBin != pair2.xBin || yBin != pair2.yBin)
			{
				remove(grid, pair2);

				// Build a single linked list
				pair2.xBin = xBin;
				pair2.yBin = yBin;
				pair2.next = grid[xBin][yBin];
				grid[xBin][yBin] = pair2;
			}

			return true;
		}

		return false;
	}

	private ArrayList<Cluster> runClosestDistancePriority(Cluster[][] grid, int nXBins, int nYBins, double r2,
			int time, double minx, double miny, double xBinWidth, double yBinWidth, ArrayList<Cluster> candidates,
			ArrayList<Cluster> singles)
	{
		int N = candidates.size();
		int candidatesProcessed = 0;
		int s = singles.size();
		final boolean trackProgress = (tracker.getClass() != NullTrackProgress.class);
		TimeCluster[][] newGrid = convertGrid(grid, nXBins, nYBins);
		while (joinClosestDistancePriority(newGrid, nXBins, nYBins, r2, time, minx, miny, xBinWidth, yBinWidth, singles))
		{
			if (tracker.stop())
				return null;

			// The number of candidates that have been processed is incremented by the number of singles
			if (trackProgress)
			{
				candidatesProcessed += singles.size() - s;
				s = singles.size();
				tracker.progress(candidatesProcessed++, N);
			}
		}

		return combine(singles, newGrid, nXBins, nYBins);
	}

	/**
	 * Search for the closest pair of clusters that are below the squared radius distance and join them. Join closest in
	 * distance and then time.
	 * 
	 * @param grid
	 * @param nXBins
	 * @param nYBins
	 * @param r2
	 *            The squared radius distance
	 * @param time
	 *            The time threshold
	 * @param yBinWidth
	 * @param xBinWidth
	 * @param miny
	 * @param minx
	 * @param singles
	 *            Add remaining clusters that have no neighbours
	 * @return True if a join was made
	 */
	private boolean joinClosestDistancePriority(TimeCluster[][] grid, final int nXBins, final int nYBins,
			final double r2, int time, double minx, double miny, double xBinWidth, double yBinWidth,
			ArrayList<Cluster> singles)
	{
		double minD = Double.POSITIVE_INFINITY;
		int minT = Integer.MAX_VALUE;
		TimeCluster pair1 = null, pair2 = null;
		final double neighbourDistance = 2 * r2;
		TimeCluster[] neighbourCells = new TimeCluster[5];
		final boolean checkPulseInterval = pulseInterval > 0;
		for (int yBin = 0; yBin < nYBins; yBin++)
		{
			for (int xBin = 0; xBin < nXBins; xBin++)
			{
				TimeCluster previous = null;
				for (TimeCluster c1 = grid[xBin][yBin]; c1 != null; c1 = (TimeCluster) c1.next)
				{
					// Build a list of which cells to compare up to a maximum of 4
					//      | 0,0  |  1,0
					// ------------+-----
					// -1,1 | 0,1  |  1,1

					// Compare to neighbours and find the closest.
					// Use either the radius threshold or the current closest distance
					// which may have been set by an earlier comparison.
					TimeCluster other = null;
					int cells = 1;
					neighbourCells[0] = (TimeCluster) c1.next;
					if (yBin < nYBins - 1)
					{
						neighbourCells[cells++] = grid[xBin][yBin + 1];
						if (xBin > 0)
							neighbourCells[cells++] = grid[xBin - 1][yBin + 1];
					}
					if (xBin < nXBins - 1)
					{
						neighbourCells[cells++] = grid[xBin + 1][yBin];
						if (yBin < nYBins - 1)
							neighbourCells[cells++] = grid[xBin + 1][yBin + 1];
					}

					for (int c = 0; c < cells; c++)
					{
						for (TimeCluster c2 = neighbourCells[c]; c2 != null; c2 = (TimeCluster) c2.next)
						{
							if (checkPulseInterval && c1.pulse != c2.pulse)
								continue;

							final int gap = c1.gap(c2);
							if (gap <= time)
							{
								final double d2 = c1.distance2(c2);
								if (d2 < neighbourDistance)
								{
									// Check if the two clusters can be merged
									if (!c1.validUnion(c2))
										continue;

									// Count any possible neighbour
									c1.neighbour++;
									c2.neighbour++;
								}
								if (d2 <= r2)
								{
									// This is within the time and distance thresholds.									
									// Find closest pair with distance priority
									if ((d2 < minD) || (d2 <= minD && gap < minT))
									{
										minD = d2;
										minT = gap;
										other = c2;
									}
								}
							}
						}
					}

					// Store the details of the closest pair
					if (other != null)
					{
						pair1 = c1;
						pair2 = other;
					}

					// Check for singles
					if (c1.neighbour == 0)
					{
						// Add singles to the singles list and remove from the grid
						singles.add(c1);
						if (previous == null)
							grid[xBin][yBin] = (TimeCluster) c1.next;
						else
							previous.next = c1.next;
					}
					else
					{
						previous = c1;
						c1.neighbour = 0;
					}
				}
			}
		}

		// Join the closest pair
		if (pair1 != null)
		{
			pair2.add(pair1);

			remove(grid, pair1);

			// Reassign the updated grid position
			final int xBin = (int) ((pair2.x - minx) / xBinWidth);
			final int yBin = (int) ((pair2.y - miny) / yBinWidth);

			if (xBin != pair2.xBin || yBin != pair2.yBin)
			{
				remove(grid, pair2);

				// Build a single linked list
				pair2.xBin = xBin;
				pair2.yBin = yBin;
				pair2.next = grid[xBin][yBin];
				grid[xBin][yBin] = pair2;
			}

			return true;
		}

		return false;
	}

	/**
	 * @return the clustering algorithm
	 */
	public ClusteringAlgorithm getClusteringAlgorithm()
	{
		return clusteringAlgorithm;
	}

	/**
	 * @param clusteringAlgorithm
	 *            The algorithm
	 */
	public void setClusteringAlgorithm(ClusteringAlgorithm clusteringAlgorithm)
	{
		this.clusteringAlgorithm = clusteringAlgorithm;
	}

	/**
	 * @return the tracker
	 */
	public TrackProgress getTracker()
	{
		return tracker;
	}

	/**
	 * @param tracker
	 *            the tracker to set
	 */
	public void setTracker(TrackProgress tracker)
	{
		if (tracker == null)
			tracker = new NullTrackProgress();
		this.tracker = tracker;
	}

	/**
	 * @return the pulse interval
	 */
	public int getPulseInterval()
	{
		return pulseInterval;
	}

	/**
	 * Set a pulse interval. Clusters will only be created by joining localisations within each pulse. Pulses are
	 * assumed to start at t=1.
	 * <p>
	 * This only applies to the algorithms that use time and distance thresholds.
	 * 
	 * @param pulseInterval
	 *            the pulse interval
	 */
	public void setPulseInterval(int pulseInterval)
	{
		this.pulseInterval = Math.max(0, pulseInterval);
	}

	/**
	 * Get the pulse for the specified time. Assumes pulses start at t=1. Returns zero if no pulse interval is defined.
	 * 
	 * @param time
	 * @return the pulse
	 */
	public int getPulse(int time)
	{
		if (pulseInterval == 0)
			return 0;
		return ((time - 1) / pulseInterval);
	}
}
