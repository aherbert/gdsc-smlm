package gdsc.smlm.search;

import gdsc.core.logging.TrackProgress;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Search a range of parameter space using a window divided into increments.
 */
public class SearchSpace
{
	private final SearchDimension[] dimensions;

	private int iteration = 0;
	private double[][] dimensionValues;
	private double[][] searchSpace;
	// This introduces a dependency on another gdsc.smlm package
	private TrackProgress tracker = null;

	/**
	 * Instantiates a new search space.
	 *
	 * @param dimensions
	 *            the dimensions
	 */
	public SearchSpace(SearchDimension[] dimensions)
	{
		if (dimensions == null || dimensions.length == 0)
			throw new IllegalArgumentException("Dimensions must not be empty");
		this.dimensions = dimensions;
	}

	/**
	 * Search the configured search space until convergence of the optimum.
	 * <p>
	 * At each iteration the search will enumerate all points in the configured search space and find the optimum. If at
	 * the bounds or the range in any dimension then the range is re-centred and the process repeated. If not at the
	 * bounds then the range is re-centred and the width of the range reduced.
	 * <p>
	 * The process iterates until the range cannot be reduced in size, or convergence is reached.
	 *
	 * @param scoreFunction
	 *            the score function
	 * @param checker
	 *            the checker
	 * @return The optimum (or null)
	 */
	public double[] search(ScoreFunction scoreFunction, ConvergenceChecker checker)
	{
		if (scoreFunction == null)
			throw new IllegalArgumentException("Score function is null");

		// Find the best individual
		double[] current = findOptimum(scoreFunction);
		double[] previous;

		boolean converged = false;
		while (!converged)
		{
			previous = current;

			if (!updateSearchSpace(current))
				break;

			// Find the optimum and check convergence
			current = findOptimum(scoreFunction);
			if (checker != null)
				converged = checker.converged(previous, current);
		}
		if (tracker != null)
			tracker.status("Converged [%d]", iteration);
		return current;
	}

	private boolean createSearchSpace()
	{
		start("Create Search Space");

		// Get the values
		int combinations = 1;
		dimensionValues = new double[dimensions.length][];
		double[] v = new double[dimensions.length];
		for (int i = 0; i < dimensions.length; i++)
		{
			dimensionValues[i] = dimensions[i].values();
			combinations *= dimensionValues[i].length;
			v[i] = dimensionValues[i][0];
		}

		// This will be a list of points enumerating the entire range
		// of the dimensions
		searchSpace = new double[combinations][];

		// Start with the min value in each dimension
		searchSpace[0] = v;
		int n = 1;
		try
		{
			// Enumerate each dimension
			for (int i = 0; i < dimensions.length; i++)
			{
				// The number of current points
				final int size = n;
				
				// Values to iterate over for this dimension
				final double[] v1 = dimensionValues[i];

				// For all the current points
				for (int j = 0; j < size; j++)
				{
					// The point values
					final double[] v2 = searchSpace[j];
					
					// We started with the min value for the dimension 
					// so go from index 1 upwards
					for (int k = 1; k < v1.length; k++)
					{
						// Create a new point with an updated value for this dimension
						final double[] v3 = v2.clone();
						v3[k] = v1[k];
						searchSpace[n++] = v3;
					}
				}
			}
		}
		catch (ArrayIndexOutOfBoundsException e)
		{
			// Return false
			e.printStackTrace();
			n = -1;
		}

		end();

		return (n == combinations);
	}

	private double[] findOptimum(ScoreFunction scoreFunction)
	{
		if (!createSearchSpace())
			return null;

		start("Find Optimum");
		double[] optimum = scoreFunction.findOptimum(searchSpace);
		end();

		return optimum;
	}

	private boolean updateSearchSpace(double[] optimum)
	{
		if (optimum == null)
			return false;

		start("Update search space");
		
		// Check if at the bounds of the dimension values
		
		// Move to the centre using the current optimum
		
		// If at bounds then check if the bounds have changed due to the move
		
		// If changed then we stick to the current range
		
		// If not changed then reduce the range
		// First check if the range can be reduced. If not then return false.

		end();

		return true;
	}

	/**
	 * @return the tracker
	 */
	public TrackProgress getTracker()
	{
		return tracker;
	}

	/**
	 * Set a tracker to allow the progress to be followed
	 * 
	 * @param tracker
	 *            the tracker to set
	 */
	public void setTracker(TrackProgress tracker)
	{
		this.tracker = tracker;
	}

	/**
	 * Get the iteration. The iteration is increased each time the population grows as part of the [grow, evaluate,
	 * select] cycle.
	 * 
	 * @return the iteration
	 */
	public int getIteration()
	{
		return iteration;
	}

	private void start(String stage)
	{
		if (tracker != null)
		{
			tracker.status(stage + " [%d]", iteration);
			tracker.progress(0);
		}
	}

	private void end()
	{
		if (tracker != null)
			tracker.progress(1);
	}

	/**
	 * Gets the current search space.
	 *
	 * @return the search space
	 */
	public double[][] getSearchSpace()
	{
		return searchSpace;
	}
}
