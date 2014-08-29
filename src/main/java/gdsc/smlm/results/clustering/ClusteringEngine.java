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

import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.results.DensityManager;
import gdsc.smlm.results.NullTrackProgress;
import gdsc.smlm.results.TrackProgress;

import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * Find clusters of points using a clustering algorithm.
 */
public class ClusteringEngine
{
	/**
	 * Used for multi-threaded clustering to store the closest pair in a region of the search space
	 */
	private class ClosestPair
	{
		double distance;
		int time;
		Object point1, point2;

		public ClosestPair(double distance, Object point1, Object point2)
		{
			this.distance = distance;
			this.point1 = point1;
			this.point2 = point2;
		}

		public ClosestPair(double distance, int time, Object point1, Object point2)
		{
			this.distance = distance;
			this.time = time;
			this.point1 = point1;
			this.point2 = point2;
		}

		public ClosestPair()
		{
		}
	}

	/**
	 * Use a runnable to allow multi-threaded operation. Input parameters
	 * that are manipulated should have synchronized methods.
	 */
	private class ParticleLinkageWorker implements Runnable
	{
		ClosestPair pair;
		ExtendedClusterPoint[][] grid;
		int nXBins;
		int nYBins;
		double r2;
		double minx;
		double miny;
		int startXBin;
		int endXBin;
		int startYBin;
		int endYBin;

		public ParticleLinkageWorker(ClosestPair pair, ExtendedClusterPoint[][] grid, int nXBins, int nYBins,
				double r2, double minx, double miny, int startXBin, int endXBin, int startYBin, int endYBin)
		{
			this.pair = pair;
			this.grid = grid;
			this.nXBins = nXBins;
			this.nYBins = nYBins;
			this.r2 = r2;
			this.minx = minx;
			this.miny = miny;
			this.startXBin = startXBin;
			this.endXBin = endXBin;
			this.startYBin = startYBin;
			this.endYBin = endYBin;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Runnable#run()
		 */
		public void run()
		{
			ClosestPair result = findClosestParticle(grid, nXBins, nYBins, r2, minx, miny, startXBin, endXBin,
					startYBin, endYBin);
			if (result != null)
			{
				pair.distance = result.distance;
				pair.point1 = result.point1;
				pair.point2 = result.point2;
			}
		}
	}

	/**
	 * Use a runnable to allow multi-threaded operation. Input parameters
	 * that are manipulated should have synchronized methods.
	 */
	private class ClosestWorker implements Runnable
	{
		ClosestPair pair;
		Cluster[][] grid;
		int nXBins;
		int nYBins;
		double r2;
		int startXBin;
		int endXBin;
		int startYBin;
		int endYBin;
		boolean single;

		public ClosestWorker(ClosestPair pair, Cluster[][] grid, int nXBins, int nYBins, double r2, int startXBin,
				int endXBin, int startYBin, int endYBin, boolean single)
		{
			this.pair = pair;
			this.grid = grid;
			this.nXBins = nXBins;
			this.nYBins = nYBins;
			this.r2 = r2;
			this.startXBin = startXBin;
			this.endXBin = endXBin;
			this.startYBin = startYBin;
			this.endYBin = endYBin;
			this.single = single;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Runnable#run()
		 */
		public void run()
		{
			ClosestPair result;
			if (single)
				result = findClosestParticle(grid, nXBins, nYBins, r2, startXBin, endXBin, startYBin, endYBin);
			else
				result = findClosest(grid, nXBins, nYBins, r2, startXBin, endXBin, startYBin, endYBin);

			if (result != null)
			{
				pair.distance = result.distance;
				pair.point1 = result.point1;
				pair.point2 = result.point2;
			}
		}
	}

	/**
	 * Use a runnable to allow multi-threaded operation. Input parameters
	 * that are manipulated should have synchronized methods.
	 */
	private class ClosestPriorityWorker implements Runnable
	{
		boolean timePriority;
		ClosestPair pair;
		TimeCluster[][] grid;
		int nXBins;
		int nYBins;
		double r2;
		int time;
		int startXBin;
		int endXBin;
		int startYBin;
		int endYBin;
		boolean single;

		public ClosestPriorityWorker(boolean timePriority, ClosestPair pair, TimeCluster[][] grid, int nXBins,
				int nYBins, double r2, int time, int startXBin, int endXBin, int startYBin, int endYBin, boolean single)
		{
			this.timePriority = timePriority;
			this.pair = pair;
			this.grid = grid;
			this.nXBins = nXBins;
			this.nYBins = nYBins;
			this.r2 = r2;
			this.time = time;
			this.startXBin = startXBin;
			this.endXBin = endXBin;
			this.startYBin = startYBin;
			this.endYBin = endYBin;
			this.single = single;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Runnable#run()
		 */
		public void run()
		{
			ClosestPair result = null;
			if (timePriority)
			{
				result = (single) ? findClosestParticleTimePriority(grid, nXBins, nYBins, r2, time, startXBin, endXBin,
						startYBin, endYBin) : findClosestTimePriority(grid, nXBins, nYBins, r2, time, startXBin,
						endXBin, startYBin, endYBin);
			}
			else
			{
				result = (single) ? findClosestParticleDistancePriority(grid, nXBins, nYBins, r2, time, startXBin,
						endXBin, startYBin, endYBin) : findClosestDistancePriority(grid, nXBins, nYBins, r2, time,
						startXBin, endXBin, startYBin, endYBin);
			}

			if (result != null)
			{
				pair.distance = result.distance;
				pair.time = result.time;
				pair.point1 = result.point1;
				pair.point2 = result.point2;
			}
		}
	}

	/**
	 * Use a runnable to allow multi-threaded operation. Input parameters
	 * that are manipulated should have synchronized methods.
	 */
	private class FindLinksWorker implements Runnable
	{
		boolean links = false;
		Cluster[][] grid;
		int nXBins;
		int nYBins;
		double r2;
		int startXBin;
		int endXBin;
		int startYBin;
		int endYBin;

		public FindLinksWorker(Cluster[][] grid, int nXBins, int nYBins, double r2, int startXBin, int endXBin,
				int startYBin, int endYBin)
		{
			this.grid = grid;
			this.nXBins = nXBins;
			this.nYBins = nYBins;
			this.r2 = r2;
			this.startXBin = startXBin;
			this.endXBin = endXBin;
			this.startYBin = startYBin;
			this.endYBin = endYBin;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Runnable#run()
		 */
		public void run()
		{
			links = findLinksAndCountNeighbours(grid, nXBins, nYBins, r2, startXBin, endXBin, startYBin, endYBin);
		}
	}

	private ExecutorService threadPool = null;
	private int xBlock, yBlock;

	private ClusteringAlgorithm clusteringAlgorithm = ClusteringAlgorithm.Pairwise;
	private TrackProgress tracker;
	private int pulseInterval = 0;
	private boolean trackJoins = false;
	private int threadCount = 1;
	private double[] intraIdDistances = null;
	private double[] interIdDistances = null;
	private int intraIdCount, interIdCount;
	private int nextClusterId;
	private boolean useRange = false;

	public ClusteringEngine()
	{
		setTracker(null);
	}

	public ClusteringEngine(int threadCount)
	{
		this.threadCount = threadCount;
		setTracker(null);
	}

	public ClusteringEngine(int threadCount, ClusteringAlgorithm custeringAlgorithm)
	{
		this.threadCount = threadCount;
		this.clusteringAlgorithm = custeringAlgorithm;
		setTracker(null);
	}

	public ClusteringEngine(int threadCount, ClusteringAlgorithm custeringAlgorithm, TrackProgress tracker)
	{
		this.threadCount = threadCount;
		this.clusteringAlgorithm = custeringAlgorithm;
		setTracker(tracker);
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
		if (clusteringAlgorithm == ClusteringAlgorithm.ParticleSingleLinkage)
			return runParticleSingleLinkage(points, radius);

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
		boolean single = false;
		switch (clusteringAlgorithm)
		{
			case Pairwise:
				clusters = runPairwise(grid, nXBins, nYBins, r2, minx, miny, xBinWidth, yBinWidth, candidates, singles);
				break;

			case PairwiseWithoutNeighbours:
				clusters = runPairwiseWithoutNeighbours(grid, nXBins, nYBins, r2, minx, miny, xBinWidth, yBinWidth,
						candidates, singles);
				break;

			case ClosestParticleTimePriority:
				single = true;

			case ClosestTimePriority:
				clusters = runClosestTimePriority(grid, nXBins, nYBins, r2, time, minx, miny, xBinWidth, yBinWidth,
						candidates, singles, single);
				break;

			case ClosestParticleDistancePriority:
				single = true;

			case ClosestDistancePriority:
				clusters = runClosestDistancePriority(grid, nXBins, nYBins, r2, time, minx, miny, xBinWidth, yBinWidth,
						candidates, singles, single);
				break;

			case ClosestParticle:
				single = true;

			case Closest:
			default:
				clusters = runClosest(grid, nXBins, nYBins, r2, minx, miny, xBinWidth, yBinWidth, candidates, singles,
						single);
		}

		tracker.progress(1);
		shutdownMultithreading();

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
	private ArrayList<Cluster> runParticleSingleLinkage(List<ClusterPoint> points, double radius)
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

		ArrayList<Cluster> clusters = runParticleSingleLinkage(grid, nXBins, nYBins, r2, minx, miny, candidates,
				singles);

		tracker.progress(1);
		shutdownMultithreading();

		tracker.log("Found %d clusters", (clusters == null) ? 0 : clusters.size());

		return clusters;
	}

	private ArrayList<Cluster> runParticleSingleLinkage(ExtendedClusterPoint[][] grid, int nXBins, int nYBins,
			double r2, double minx, double miny, ArrayList<ExtendedClusterPoint> candidates, ArrayList<Cluster> singles)
	{
		int N = candidates.size();
		int candidatesProcessed = 0;
		nextClusterId = 0; // Incremented within joinClosestParticle(...)

		// Used to store the cluster for each candidate
		int[] clusterId = new int[candidates.size()];

		int nProcessed = 0;

		if (trackJoins)
		{
			interIdDistances = new double[candidates.size()];
			intraIdDistances = new double[candidates.size()];
			interIdCount = intraIdCount = 0;
		}

		initialiseMultithreading(nXBins, nYBins);
		while ((nProcessed = joinClosestParticle(grid, nXBins, nYBins, r2, minx, miny, clusterId)) > 0)
		{
			if (tracker.isEnded())
				return null;
			candidatesProcessed += nProcessed;
			tracker.progress(candidatesProcessed, N);
		}

		if (trackJoins)
		{
			interIdDistances = Arrays.copyOf(interIdDistances, interIdCount);
			intraIdDistances = Arrays.copyOf(intraIdDistances, intraIdCount);
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

				//// Check is a neighbour
				//boolean neighbour = false;
				//for (int j = 0; j < clusterId.length; j++)
				//{
				//	if (i == j)
				//		continue;
				//	ClusterPoint p = candidates.get(j);
				//	if (originalPoint.distance2(p) < r2)
				//	{
				//		neighbour = true;
				//		break;
				//	}
				//}
				//if (neighbour)
				//	tracker.log("Failed to assign a cluster to a candidate particle: " + i);
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

	private void initialiseMultithreading(int nXBins, int nYBins)
	{
		final int MIN_BLOCK_SIZE = 3;

		// Do not use threads if the number of bins is small
		if (nXBins < MIN_BLOCK_SIZE && nYBins < MIN_BLOCK_SIZE)
			return;

		if (threadCount > 1)
		{
			// Set up for multi-threading
			threadPool = Executors.newFixedThreadPool(threadCount);

			// Ensure a minimum block size to avoid wasting time.
			xBlock = Math.max(nXBins / threadCount, MIN_BLOCK_SIZE);
			yBlock = Math.max(nYBins / threadCount, MIN_BLOCK_SIZE);

			//System.out.printf("Block size %d x %d = %d\n", xBlock, yBlock, countBlocks(nXBins, nYBins));
			// Increment the block size until the number of blocks to process is just above the thread count.
			// This reduces thread overhead but still processes across many threads.
			int counter = 0;
			while (countBlocks(nXBins, nYBins) > 2 * threadCount)
			{
				if (counter++ % 2 == 0)
					xBlock++;
				else
					yBlock++;
			}
			//System.out.printf("Block size %d x %d = %d\n", xBlock, yBlock, countBlocks(nXBins, nYBins));
		}
	}

	private void shutdownMultithreading()
	{
		if (threadPool != null)
		{
			threadPool.shutdownNow();
			threadPool = null;
		}
	}

	private int countBlocks(int nXBins, int nYBins)
	{
		int nBlocks = 0;
		for (int startYBin = 0; startYBin < nYBins; startYBin += yBlock)
		{
			for (int startXBin = 0; startXBin < nXBins; startXBin += xBlock)
			{
				nBlocks++;
			}
		}
		return nBlocks;
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
	 * @param miny
	 * @param minx
	 * @param clusterId
	 * @return The number of points assigned to clusters (either 0, 1, or 2)
	 */
	private int joinClosestParticle(ExtendedClusterPoint[][] grid, final int nXBins, final int nYBins, final double r2,
			double minx, double miny, int[] clusterId)
	{
		ClosestPair closest = null;

		// Blocks must be overlapping by one bin to calculate all distances. There is no point multi-threading
		// if there are not enough blocks since each overlap is double processed.		

		if (threadPool == null)
		{
			closest = findClosestParticle(grid, nXBins, nYBins, r2, minx, miny, 0, nXBins, 0, nYBins);
		}
		else
		{
			// Use threads to find the closest pairs in blocks
			List<Future<?>> futures = new LinkedList<Future<?>>();
			List<ClosestPair> results = new LinkedList<ClosestPair>();

			for (int startYBin = 0; startYBin < nYBins; startYBin += yBlock)
			{
				int endYBin = Math.min(nYBins, startYBin + yBlock);
				for (int startXBin = 0; startXBin < nXBins; startXBin += xBlock)
				{
					int endXBin = Math.min(nXBins, startXBin + xBlock);
					//System.out.printf("Block [%d-%d, %d-%d]\n", startXBin, endXBin, startYBin, endYBin);

					ClosestPair pair = new ClosestPair();
					results.add(pair);
					futures.add(threadPool.submit(new ParticleLinkageWorker(pair, grid, nXBins, nYBins, r2, minx, miny,
							startXBin, endXBin, startYBin, endYBin)));
				}
			}

			// Finish processing data
			Utils.waitForCompletion(futures);
			futures.clear();

			// Find the closest pair from all the results
			for (ClosestPair result : results)
			{
				if (result.point1 != null)
				{
					if (closest == null || result.distance < closest.distance)
						closest = result;
				}
			}
		}

		// Assign the closest pair.
		if (closest != null)
		{
			ExtendedClusterPoint pair1 = (ExtendedClusterPoint) closest.point1;
			ExtendedClusterPoint pair2 = (ExtendedClusterPoint) closest.point2;

			int nProcessed = 1;

			if (pair1.inCluster && pair2.inCluster)
			{
				// Error
				throw new RuntimeException("Linkage between two particles already in a cluster");
			}
			else if (pair1.inCluster)
			{
				clusterId[pair2.id] = clusterId[pair1.id];
				pair2.inCluster = true;
			}
			else if (pair2.inCluster)
			{
				clusterId[pair1.id] = clusterId[pair2.id];
				pair1.inCluster = true;
			}
			// Create a new cluster if necessary
			else
			{
				nProcessed = 2;
				clusterId[pair1.id] = clusterId[pair2.id] = ++nextClusterId;
				pair1.inCluster = pair2.inCluster = true;
			}

			if (trackJoins)
			{
				if (pair1.next.id == pair2.next.id)
					intraIdDistances[intraIdCount++] = Math.sqrt(closest.distance);
				else
					interIdDistances[interIdCount++] = Math.sqrt(closest.distance);
			}

			return nProcessed;
		}

		// This should not happen if the density manager counted neighbours correctly
		// (i.e. all points should have at least one neighbour)		
		return 0;
	}

	/**
	 * Search for the closest pair of particles, one of which is not in a cluster, below the squared radius
	 * distance
	 * 
	 * @param grid
	 * @param nXBins
	 * @param nYBins
	 * @param r2
	 *            The squared radius distance
	 * @param miny
	 * @param minx
	 * @param startXBin
	 * @param endXBin
	 * @param startYBin
	 * @param endYBin
	 * @return The pair of closest points (or null)
	 */
	private ClosestPair findClosestParticle(ExtendedClusterPoint[][] grid, final int nXBins, final int nYBins,
			final double r2, double minx, double miny, int startXBin, int endXBin, int startYBin, int endYBin)
	{
		double min = r2;
		ExtendedClusterPoint pair1 = null, pair2 = null;
		for (int yBin = startYBin; yBin < endYBin; yBin++)
		{
			for (int xBin = startXBin; xBin < endXBin; xBin++)
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
		if (pair1 != null)
			return new ClosestPair(min, pair1, pair2);
		return null;
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
	 * Check there are at least two different time points in the data. Also check if any candidates have a different
	 * start and end time.
	 * 
	 * @param candidates
	 * @return true if there are no different time points
	 */
	private boolean noTimeInformation(ArrayList<Cluster> candidates)
	{
		useRange = checkForTimeRange(candidates);
		final int firstT = candidates.get(0).head.start;
		if (useRange)
		{
			final int lastT = candidates.get(0).head.end;
			for (Cluster c : candidates)
			{
				if (firstT != c.head.start)
					return false;
				if (lastT != c.head.end)
					return false;
			}
		}
		else
		{
			for (Cluster c : candidates)
			{
				if (firstT != c.head.start)
					return false;
			}
		}

		return true;
	}

	/**
	 * Check if any of the candidates have a different start and end time
	 * 
	 * @param candidates
	 * @return
	 */
	private boolean checkForTimeRange(ArrayList<Cluster> candidates)
	{
		for (Cluster c : candidates)
			if (c.head.start != c.head.end)
				return true;
		return false;
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
			if (tracker.isEnded())
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
		// The find method is multi-threaded
		// The remaining join and update of the grid is single threaded
		initialiseMultithreading(nXBins, nYBins);

		// Store if the clusters have any neighbours within sqrt(2)*r and remove them 
		// from the next loop.
		int N = candidates.size();
		ArrayList<Cluster> joined = new ArrayList<Cluster>();
		while (findLinksAndCountNeighbours(grid, nXBins, nYBins, r2, singles))
		{
			if (tracker.isEnded())
				return null;

			int joins = joinLinks(grid, nXBins, nYBins, r2, candidates, joined, singles);
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
		if (threadPool == null)
		{
			return findLinksAndCountNeighbours(grid, nXBins, nYBins, r2, 0, nXBins, 0, nYBins);
		}
		else
		{
			// Use threads to find the closest pairs in blocks
			List<Future<?>> futures = new LinkedList<Future<?>>();
			List<FindLinksWorker> results = new LinkedList<FindLinksWorker>();

			for (int startYBin = 0; startYBin < nYBins; startYBin += yBlock)
			{
				int endYBin = Math.min(nYBins, startYBin + yBlock);
				for (int startXBin = 0; startXBin < nXBins; startXBin += xBlock)
				{
					int endXBin = Math.min(nXBins, startXBin + xBlock);

					FindLinksWorker worker = new FindLinksWorker(grid, nXBins, nYBins, r2, startXBin, endXBin,
							startYBin, endYBin);
					results.add(worker);
					futures.add(threadPool.submit(worker));
				}
			}

			// Finish processing data
			Utils.waitForCompletion(futures);
			futures.clear();

			for (FindLinksWorker worker : results)
			{
				if (worker.links)
					return true;
			}
			return false;
		}
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
	private boolean findLinksAndCountNeighbours(Cluster[][] grid, final int nXBins, final int nYBins, final double r2,
			int startXBin, int endXBin, int startYBin, int endYBin)
	{
		Cluster[] neighbours = new Cluster[5];
		boolean linked = false;
		final double neighbourDistance = 2 * r2;

		for (int yBin = startYBin; yBin < endYBin; yBin++)
		{
			for (int xBin = startXBin; xBin < endXBin; xBin++)
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
							if (d2 < neighbourDistance)
							{
								c1.incrementNeighbour();
								c2.incrementNeighbour();
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
						c1.linkSynchronized(other, min);
						linked = true;
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
	 * @param singles
	 *            Add any clusters with no neighbours
	 * @return The number of joins that were made
	 */
	private int joinLinks(Cluster[][] grid, int nXBins, int nYBins, double r2, ArrayList<Cluster> candidates,
			ArrayList<Cluster> joined, ArrayList<Cluster> singles)
	{
		candidates.clear();
		joined.clear();

		double min = r2;
		Cluster cluster1 = null, cluster2 = null;
		for (int yBin = 0; yBin < nYBins; yBin++)
		{
			for (int xBin = 0; xBin < nXBins; xBin++)
			{
				Cluster previous = null;
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

						// Store all remaining clusters
						if (c1.n != 0)
						{
							candidates.add(c1);
							c1.neighbour = joinFlag;
						}
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
	 * @param single
	 *            True if only singles can be joined to another cluster
	 * @return
	 */
	private ArrayList<Cluster> runClosest(Cluster[][] grid, int nXBins, int nYBins, double r2, double minx,
			double miny, double xBinWidth, double yBinWidth, ArrayList<Cluster> candidates, ArrayList<Cluster> singles,
			boolean single)
	{
		int N = candidates.size();
		int candidatesProcessed = 0;
		int s = singles.size();
		final boolean trackProgress = (tracker.getClass() != NullTrackProgress.class);
		initialiseMultithreading(nXBins, nYBins);
		while (joinClosest(grid, nXBins, nYBins, r2, minx, miny, xBinWidth, yBinWidth, single))
		{
			if (tracker.isEnded())
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
	 * @param single
	 *            True if only singles can be joined to another cluster
	 * @return True if a join was made
	 */
	private boolean joinClosest(Cluster[][] grid, final int nXBins, final int nYBins, final double r2, double minx,
			double miny, double xBinWidth, double yBinWidth, boolean single)
	{
		ClosestPair closest = null;

		if (threadPool == null)
		{
			closest = (single) ? findClosestParticle(grid, nXBins, nYBins, r2, 0, nXBins, 0, nYBins) : findClosest(
					grid, nXBins, nYBins, r2, 0, nXBins, 0, nYBins);
		}
		else
		{
			// Use threads to find the closest pairs in blocks
			List<Future<?>> futures = new LinkedList<Future<?>>();
			List<ClosestPair> results = new LinkedList<ClosestPair>();

			for (int startYBin = 0; startYBin < nYBins; startYBin += yBlock)
			{
				int endYBin = Math.min(nYBins, startYBin + yBlock);
				for (int startXBin = 0; startXBin < nXBins; startXBin += xBlock)
				{
					int endXBin = Math.min(nXBins, startXBin + xBlock);
					//System.out.printf("Block [%d-%d, %d-%d]\n", startXBin, endXBin, startYBin, endYBin);

					ClosestPair pair = new ClosestPair();
					results.add(pair);
					futures.add(threadPool.submit(new ClosestWorker(pair, grid, nXBins, nYBins, r2, startXBin, endXBin,
							startYBin, endYBin, single)));
				}
			}

			// Finish processing data
			Utils.waitForCompletion(futures);

			// Find the closest pair from all the results
			for (ClosestPair result : results)
			{
				if (result.point1 != null)
				{
					if (closest == null || result.distance < closest.distance)
						closest = result;
				}
			}
		}

		// Join the closest pair
		if (closest != null)
		{
			Cluster pair1 = (Cluster) closest.point1;
			Cluster pair2 = (Cluster) closest.point2;

			if (single)
			{
				// Check
				if (pair1.n > 1 && pair2.n > 1)
				{
					throw new RuntimeException(
							"Linkage between two clusters (not a single particle and a single/cluster)");
				}

				// Add the single to the cluster
				if (pair2.n < pair1.n)
				{
					Cluster tmp = pair1;
					pair1 = pair2;
					pair2 = tmp;
				}
			}

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
	 * Search for the closest pair of clusters that are below the squared radius distance
	 * 
	 * @param grid
	 * @param nXBins
	 * @param nYBins
	 * @param r2
	 *            The squared radius distance
	 * @param startXBin
	 * @param endXBin
	 * @param startYBin
	 * @param endYBin
	 * @return The closest pair
	 */
	private ClosestPair findClosest(Cluster[][] grid, final int nXBins, final int nYBins, final double r2,
			int startXBin, int endXBin, int startYBin, int endYBin)
	{
		double min = r2;
		Cluster pair1 = null, pair2 = null;
		for (int yBin = startYBin; yBin < endYBin; yBin++)
		{
			for (int xBin = startXBin; xBin < endXBin; xBin++)
			{
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
		if (pair1 != null)
			return new ClosestPair(min, pair1, pair2);
		return null;
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

	/**
	 * Search for the closest pair of a single particle and any existing single/cluster that are below the squared
	 * radius distance
	 * 
	 * @param grid
	 * @param nXBins
	 * @param nYBins
	 * @param r2
	 *            The squared radius distance
	 * @param startXBin
	 * @param endXBin
	 * @param startYBin
	 * @param endYBin
	 * @return The closest pair
	 */
	private ClosestPair findClosestParticle(Cluster[][] grid, final int nXBins, final int nYBins, final double r2,
			int startXBin, int endXBin, int startYBin, int endYBin)
	{
		double min = r2;
		Cluster pair1 = null, pair2 = null;
		for (int yBin = startYBin; yBin < endYBin; yBin++)
		{
			for (int xBin = startXBin; xBin < endXBin; xBin++)
			{
				for (Cluster c1 = grid[xBin][yBin]; c1 != null; c1 = c1.next)
				{
					final boolean cluster1 = c1.n > 1;

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
						if (cluster1 && c2.n > 1)
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
						for (Cluster c2 = grid[xBin][yBin + 1]; c2 != null; c2 = c2.next)
						{
							if (cluster1 && c2.n > 1)
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
							for (Cluster c2 = grid[xBin - 1][yBin + 1]; c2 != null; c2 = c2.next)
							{
								if (cluster1 && c2.n > 1)
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
						for (Cluster c2 = grid[xBin + 1][yBin]; c2 != null; c2 = c2.next)
						{
							if (cluster1 && c2.n > 1)
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
							for (Cluster c2 = grid[xBin + 1][yBin + 1]; c2 != null; c2 = c2.next)
							{
								if (cluster1 && c2.n > 1)
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
		if (pair1 != null)
			return new ClosestPair(min, pair1, pair2);
		return null;
	}

	/**
	 * The process should iterate finding the closest nodes, joining them and repeating. Join closest in
	 * time and then distance.
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
	 * @param single
	 *            True if only singles can be joined to another cluster
	 * @return
	 */
	private ArrayList<Cluster> runClosestTimePriority(Cluster[][] grid, int nXBins, int nYBins, double r2, int time,
			double minx, double miny, double xBinWidth, double yBinWidth, ArrayList<Cluster> candidates,
			ArrayList<Cluster> singles, boolean single)
	{
		int N = candidates.size();
		int candidatesProcessed = 0;
		int s = singles.size();
		final boolean trackProgress = (tracker.getClass() != NullTrackProgress.class);
		TimeCluster[][] newGrid = convertGrid(grid, nXBins, nYBins);
		initialiseMultithreading(nXBins, nYBins);
		while (joinClosestTimePriority(newGrid, nXBins, nYBins, r2, time, minx, miny, xBinWidth, yBinWidth, singles,
				single))
		{
			if (tracker.isEnded())
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
	 * @param single
	 *            True if only singles can be joined to another cluster
	 * @return True if a join was made
	 */
	private boolean joinClosestTimePriority(TimeCluster[][] grid, final int nXBins, final int nYBins, final double r2,
			final int time, double minx, double miny, double xBinWidth, double yBinWidth, ArrayList<Cluster> singles,
			boolean single)
	{
		ClosestPair closest = null;

		if (threadPool == null)
		{
			closest = (single) ? findClosestParticleTimePriority(grid, nXBins, nYBins, r2, time, 0, nXBins, 0, nYBins)
					: findClosestTimePriority(grid, nXBins, nYBins, r2, time, 0, nXBins, 0, nYBins);
		}
		else
		{
			// Use threads to find the closest pairs in blocks
			List<Future<?>> futures = new LinkedList<Future<?>>();
			List<ClosestPair> results = new LinkedList<ClosestPair>();

			for (int startYBin = 0; startYBin < nYBins; startYBin += yBlock)
			{
				int endYBin = Math.min(nYBins, startYBin + yBlock);
				for (int startXBin = 0; startXBin < nXBins; startXBin += xBlock)
				{
					int endXBin = Math.min(nXBins, startXBin + xBlock);
					//System.out.printf("Block [%d-%d, %d-%d]\n", startXBin, endXBin, startYBin, endYBin);

					ClosestPair pair = new ClosestPair();
					results.add(pair);
					futures.add(threadPool.submit(new ClosestPriorityWorker(true, pair, grid, nXBins, nYBins, r2, time,
							startXBin, endXBin, startYBin, endYBin, single)));
				}
			}

			// Finish processing data
			Utils.waitForCompletion(futures);
			futures.clear();

			// Find the closest pair from all the results
			for (ClosestPair result : results)
			{
				if (result.point1 != null)
				{
					if (closest == null || (result.time < closest.time) ||
							(result.time <= closest.time && result.distance < closest.distance))
						closest = result;
				}
			}
		}

		// Join the closest pair
		if (closest != null)
		{
			TimeCluster pair1 = (TimeCluster) closest.point1;
			TimeCluster pair2 = (TimeCluster) closest.point2;

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
	 * Search for the closest pair of clusters that are below the squared radius distance. Find the closest in
	 * time and then distance.
	 * 
	 * @param grid
	 * @param nXBins
	 * @param nYBins
	 * @param r2
	 *            The squared radius distance
	 * @param time
	 * @param startXBin
	 * @param endXBin
	 * @param startYBin
	 * @param endYBin
	 * @return True if a join was made
	 */
	private ClosestPair findClosestTimePriority(TimeCluster[][] grid, final int nXBins, final int nYBins,
			final double r2, final int time, int startXBin, int endXBin, int startYBin, int endYBin)
	{
		double minD = Double.POSITIVE_INFINITY;
		int minT = Integer.MAX_VALUE;
		TimeCluster pair1 = null, pair2 = null;
		TimeCluster[] neighbourCells = new TimeCluster[5];
		final boolean checkPulseInterval = pulseInterval > 0;
		for (int yBin = startYBin; yBin < endYBin; yBin++)
		{
			for (int xBin = startXBin; xBin < endXBin; xBin++)
			{
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

								if (d2 <= r2)
								{
									// Check if the two clusters can be merged
									if (gap == 0 && validUnion(c1, c2))
										continue;

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
				}
			}
		}
		if (pair1 != null)
			return new ClosestPair(minD, minT, pair1, pair2);
		return null;
	}

	/**
	 * Check if the union of two clusters is valid
	 * 
	 * @param c1
	 * @param c2
	 * @return
	 */
	private boolean validUnion(TimeCluster c1, TimeCluster c2)
	{
		if (useRange)
			return c1.validUnionRange(c2);
		else
			return c1.validUnion(c2);
	}

	/**
	 * Search for the closest pair of a single particle and any existing single/cluster that are below the squared
	 * radius distance. Find the closest in time and then distance.
	 * 
	 * @param grid
	 * @param nXBins
	 * @param nYBins
	 * @param r2
	 *            The squared radius distance
	 * @param time
	 * @param startXBin
	 * @param endXBin
	 * @param startYBin
	 * @param endYBin
	 * @return True if a join was made
	 */
	private ClosestPair findClosestParticleTimePriority(TimeCluster[][] grid, final int nXBins, final int nYBins,
			final double r2, final int time, int startXBin, int endXBin, int startYBin, int endYBin)
	{
		double minD = Double.POSITIVE_INFINITY;
		int minT = Integer.MAX_VALUE;
		TimeCluster pair1 = null, pair2 = null;
		TimeCluster[] neighbourCells = new TimeCluster[5];
		final boolean checkPulseInterval = pulseInterval > 0;
		for (int yBin = startYBin; yBin < endYBin; yBin++)
		{
			for (int xBin = startXBin; xBin < endXBin; xBin++)
			{
				for (TimeCluster c1 = grid[xBin][yBin]; c1 != null; c1 = (TimeCluster) c1.next)
				{
					final boolean cluster1 = c1.n > 1;

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
							if (cluster1 && c2.n > 1)
								continue;
							if (checkPulseInterval && c1.pulse != c2.pulse)
								continue;

							final int gap = c1.gap(c2);
							if (gap <= time)
							{
								final double d2 = c1.distance2(c2);

								if (d2 <= r2)
								{
									// Check if the two clusters can be merged
									if (gap == 0 && validUnion(c1, c2))
										continue;

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
				}
			}
		}
		if (pair1 != null)
			return new ClosestPair(minD, minT, pair1, pair2);
		return null;
	}

	/**
	 * The process should iterate finding the closest nodes, joining them and repeating. Join closest in
	 * distance and then time.
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
	 * @param single
	 *            True if only singles can be joined to another cluster
	 * @return
	 */
	private ArrayList<Cluster> runClosestDistancePriority(Cluster[][] grid, int nXBins, int nYBins, double r2,
			int time, double minx, double miny, double xBinWidth, double yBinWidth, ArrayList<Cluster> candidates,
			ArrayList<Cluster> singles, boolean single)
	{
		int N = candidates.size();
		int candidatesProcessed = 0;
		int s = singles.size();
		final boolean trackProgress = (tracker.getClass() != NullTrackProgress.class);
		TimeCluster[][] newGrid = convertGrid(grid, nXBins, nYBins);
		initialiseMultithreading(nXBins, nYBins);
		while (joinClosestDistancePriority(newGrid, nXBins, nYBins, r2, time, minx, miny, xBinWidth, yBinWidth,
				singles, single))
		{
			if (tracker.isEnded())
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
	 * @param single
	 *            True if only singles can be joined to another cluster
	 * @return True if a join was made
	 */
	private boolean joinClosestDistancePriority(TimeCluster[][] grid, final int nXBins, final int nYBins,
			final double r2, int time, double minx, double miny, double xBinWidth, double yBinWidth,
			ArrayList<Cluster> singles, boolean single)
	{
		ClosestPair closest = null;

		if (threadPool == null)
		{
			closest = (single) ? findClosestParticleDistancePriority(grid, nXBins, nYBins, r2, time, 0, nXBins, 0,
					nYBins) : findClosestDistancePriority(grid, nXBins, nYBins, r2, time, 0, nXBins, 0, nYBins);
		}
		else
		{
			// Use threads to find the closest pairs in blocks
			List<Future<?>> futures = new LinkedList<Future<?>>();
			List<ClosestPair> results = new LinkedList<ClosestPair>();

			for (int startYBin = 0; startYBin < nYBins; startYBin += yBlock)
			{
				int endYBin = Math.min(nYBins, startYBin + yBlock);
				for (int startXBin = 0; startXBin < nXBins; startXBin += xBlock)
				{
					int endXBin = Math.min(nXBins, startXBin + xBlock);

					ClosestPair pair = new ClosestPair();
					results.add(pair);
					futures.add(threadPool.submit(new ClosestPriorityWorker(false, pair, grid, nXBins, nYBins, r2,
							time, startXBin, endXBin, startYBin, endYBin, single)));
				}
			}

			// Finish processing data
			Utils.waitForCompletion(futures);
			futures.clear();

			// Find the closest pair from all the results
			for (ClosestPair result : results)
			{
				if (result.point1 != null)
				{
					if (closest == null || (result.distance < closest.distance) ||
							(result.distance <= closest.distance && result.time < closest.time))
						closest = result;
				}
			}
		}

		// Join the closest pair
		if (closest != null)
		{
			TimeCluster pair1 = (TimeCluster) closest.point1;
			TimeCluster pair2 = (TimeCluster) closest.point2;

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
	 * Search for the closest pair of clusters that are below the squared radius distance. Find the closest in
	 * distance and then time.
	 * 
	 * @param grid
	 * @param nXBins
	 * @param nYBins
	 * @param r2
	 *            The squared radius distance
	 * @param time
	 * @param startXBin
	 * @param endXBin
	 * @param startYBin
	 * @param endYBin
	 * @return True if a join was made
	 */
	private ClosestPair findClosestDistancePriority(TimeCluster[][] grid, final int nXBins, final int nYBins,
			final double r2, final int time, int startXBin, int endXBin, int startYBin, int endYBin)
	{
		double minD = Double.POSITIVE_INFINITY;
		int minT = Integer.MAX_VALUE;
		TimeCluster pair1 = null, pair2 = null;
		TimeCluster[] neighbourCells = new TimeCluster[5];
		final boolean checkPulseInterval = pulseInterval > 0;
		for (int yBin = startYBin; yBin < endYBin; yBin++)
		{
			for (int xBin = startXBin; xBin < endXBin; xBin++)
			{
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
								if (d2 <= r2)
								{
									// Check if the two clusters can be merged
									if (gap == 0 && validUnion(c1, c2))
										continue;

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
				}
			}
		}
		if (pair1 != null)
			return new ClosestPair(minD, minT, pair1, pair2);
		return null;
	}

	/**
	 * Search for the closest pair of a single particle and any existing single/cluster that are below the squared
	 * radius distance. Find the closest in distance and then time.
	 * 
	 * @param grid
	 * @param nXBins
	 * @param nYBins
	 * @param r2
	 *            The squared radius distance
	 * @param time
	 * @param startXBin
	 * @param endXBin
	 * @param startYBin
	 * @param endYBin
	 * @return True if a join was made
	 */
	private ClosestPair findClosestParticleDistancePriority(TimeCluster[][] grid, final int nXBins, final int nYBins,
			final double r2, final int time, int startXBin, int endXBin, int startYBin, int endYBin)
	{
		double minD = Double.POSITIVE_INFINITY;
		int minT = Integer.MAX_VALUE;
		TimeCluster pair1 = null, pair2 = null;
		TimeCluster[] neighbourCells = new TimeCluster[5];
		final boolean checkPulseInterval = pulseInterval > 0;
		for (int yBin = startYBin; yBin < endYBin; yBin++)
		{
			for (int xBin = startXBin; xBin < endXBin; xBin++)
			{
				for (TimeCluster c1 = grid[xBin][yBin]; c1 != null; c1 = (TimeCluster) c1.next)
				{
					final boolean cluster1 = c1.n > 1;

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
							if (cluster1 && c2.n > 1)
								continue;
							if (checkPulseInterval && c1.pulse != c2.pulse)
								continue;

							final int gap = c1.gap(c2);
							if (gap <= time)
							{
								final double d2 = c1.distance2(c2);
								if (d2 <= r2)
								{
									// Check if the two clusters can be merged
									if (gap == 0 && validUnion(c1, c2))
										continue;

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
				}
			}
		}
		if (pair1 != null)
			return new ClosestPair(minD, minT, pair1, pair2);
		return null;
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

	/**
	 * @return true if recording the distances between particles that were joined.
	 */
	public boolean isTrackJoins()
	{
		return trackJoins;
	}

	/**
	 * Set to true to record the distances between particles that were joined. Only applies to the particle linkage
	 * algorithm. The distances can be retrieved after the {@link #findClusters(List, double)} method has been called.
	 * 
	 * @param trackJoins
	 */
	public void setTrackJoins(boolean trackJoins)
	{
		this.trackJoins = trackJoins;
	}

	/**
	 * @return the intra-Id distances from joins in the particle linkage algorithm
	 */
	public double[] getIntraIdDistances()
	{
		return intraIdDistances;
	}

	/**
	 * @return the inter-Id distances from joins in the particle linkage algorithm
	 */
	public double[] getInterIdDistances()
	{
		return interIdDistances;
	}

	/**
	 * @return the thread count for multi-thread compatible algorithms
	 */
	public int getThreadCount()
	{
		return threadCount;
	}

	/**
	 * @param threadCount
	 *            the thread count for multi-thread compatible algorithms
	 */
	public void setThreadCount(int threadCount)
	{
		this.threadCount = threadCount;
	}
}
