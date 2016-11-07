package gdsc.smlm.search;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

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
	private SearchDimension[] dimensions;
	private int iteration = 0;
	private double[][] searchSpace;

	private double[][] dimensionValues;
	private double[][] scoredSearchSpace;
	private ArrayList<String> scoredSearchSpaceHash = new ArrayList<String>();
	private HashSet<String> coveredSpace = new HashSet<String>();
	private StringBuilder sb = new StringBuilder();

	// This introduces a dependency on another gdsc.smlm package
	private TrackProgress tracker = null;

	/**
	 * Clone the search dimensions
	 * 
	 * @param dimensions
	 *            the search dimensions
	 * @return the clone
	 */
	public static SearchDimension[] clone(SearchDimension[] dimensions)
	{
		if (dimensions == null)
			return null;
		SearchDimension[] dimensions2 = new SearchDimension[dimensions.length];
		for (int i = 0; i < dimensions.length; i++)
			dimensions2[i] = dimensions[i].clone();
		return dimensions2;
	}

	/**
	 * Search the configured search space until convergence of the optimum.
	 * <p>
	 * At each iteration the search will enumerate all points in the configured search space and find the optimum. If at
	 * the bounds or the range in any dimension then the range is re-centred and the process repeated. During repeats
	 * only those points that have not yet been evaluated will be passed to the score function. If not at the bounds
	 * then the range is re-centred and the width of the range reduced.
	 * <p>
	 * The process iterates until the range cannot be reduced in size, or convergence is reached. The input dimensions
	 * are modified during the search. Use the clone(SearchDimension[]) method to create a copy.
	 *
	 * @param dimensions
	 *            the dimensions
	 * @param scoreFunction
	 *            the score function
	 * @param checker
	 *            the checker
	 * @return The optimum (or null)
	 */
	public <T extends Comparable<T>> ScoreResult<T> search(SearchDimension[] dimensions, ScoreFunction<T> scoreFunction,
			ConvergenceChecker checker)
	{
		if (dimensions == null || dimensions.length == 0)
			throw new IllegalArgumentException("Dimensions must not be empty");
		if (scoreFunction == null)
			throw new IllegalArgumentException("Score function is null");
		this.dimensions = dimensions;

		// Find the best individual
		ScoreResult<T> current = findOptimum(scoreFunction, null);
		ScoreResult<T> previous = null;

		boolean converged = false;
		while (!converged)
		{
			previous = current;

			if (!updateSearchSpace(current))
				break;

			// Find the optimum and check convergence
			current = findOptimum(scoreFunction, current);
			if (current == null)
				break;
			if (checker != null)
				converged = checker.converged(previous, current);
		}
		if (tracker != null)
			tracker.status("Converged [%d]", iteration);

		// Free memory
		dimensionValues = null;
		scoredSearchSpace = null;
		scoredSearchSpaceHash.clear();
		coveredSpace.clear();

		return current;
	}

	private <T extends Comparable<T>> ScoreResult<T> findOptimum(ScoreFunction<T> scoreFunction, ScoreResult<T> current)
	{
		if (!createSearchSpace())
			return null;

		scoredSearchSpace = searchSpace;
		scoredSearchSpaceHash.clear();

		if (!coveredSpace.isEmpty())
		{
			// Check we do not recompute scores
			scoredSearchSpace = new double[searchSpace.length][];
			scoredSearchSpaceHash.ensureCapacity(searchSpace.length);
			int size = 0;
			for (int i = 0; i < searchSpace.length; i++)
			{
				final String hash = generateHashString(searchSpace[i]);
				if (!coveredSpace.contains(hash))
				{
					scoredSearchSpace[size++] = searchSpace[i];
					scoredSearchSpaceHash.add(hash);
				}
			}
			if (size == 0)
				// We have scored everything already so return the current best
				return current;

			scoredSearchSpace = Arrays.copyOf(scoredSearchSpace, size);
		}

		start("Find Optimum");
		ScoreResult<T> optimum = scoreFunction.findOptimum(scoredSearchSpace);
		end();

		if (current != null)
		{
			if (current.compareTo(optimum) < 0)
				optimum = current;
		}

		return optimum;
	}

	private String generateHashString(double[] values)
	{
		for (int i = 0; i < values.length; i++)
			sb.append(values[i]);
		final String hash = sb.toString();
		sb.setLength(0);
		return hash;
	}

	private boolean createSearchSpace()
	{
		start("Create Search Space");

		// Get the values
		dimensionValues = new double[dimensions.length][];
		for (int i = 0; i < dimensions.length; i++)
		{
			dimensionValues[i] = dimensions[i].values();
		}
		searchSpace = createSearchSpace(dimensions, dimensionValues);

		end();
		return searchSpace != null;
	}

	public static double[][] createSearchSpace(SearchDimension[] dimensions)
	{
		// Get the values
		double[][] dimensionValues = new double[dimensions.length][];
		for (int i = 0; i < dimensions.length; i++)
		{
			dimensionValues[i] = dimensions[i].values();
		}
		return createSearchSpace(dimensions, dimensionValues);
	}

	private static double[][] createSearchSpace(SearchDimension[] dimensions, double[][] dimensionValues)
	{
		// Get the values
		int combinations = 1;
		double[] v = new double[dimensions.length];
		for (int i = 0; i < dimensions.length; i++)
		{
			combinations *= dimensionValues[i].length;
			v[i] = dimensionValues[i][0];
		}

		// This will be a list of points enumerating the entire range
		// of the dimensions
		double[][] searchSpace = new double[combinations][];

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

		return (n == combinations) ? searchSpace : null;
	}

	private boolean updateSearchSpace(ScoreResult<?> current)
	{
		if (current == null)
			return false;

		start("Update search space");

		// Process each dimension
		final double[] p = current.point;
		boolean changed = false;
		for (int i = 0; i < dimensions.length; i++)
		{
			// Check if at the bounds of the dimension values
			final double[] values = dimensionValues[i];
			final boolean atBounds = (p[i] == values[0] || p[i] == values[values.length - 1]);

			// Move to the centre using the current optimum
			dimensions[i].setCentre(p[i]);

			// If at bounds then check if the bounds have changed due to the move
			if (atBounds)
			{
				if (changed(values, dimensions[i].values()))
					changed = true;
			}
		}

		if (changed)
		{
			// If changed then we stick to the current range.
			// Store all the scored search space so we do not recompute it.
			if (scoredSearchSpaceHash.isEmpty())
			{
				// Compute the hash for each item scored
				for (int i = 0; i < scoredSearchSpace.length; i++)
					coveredSpace.add(generateHashString(scoredSearchSpace[i]));
			}
			else
			{
				// Hash was already computed
				for (int i = 0; i < scoredSearchSpaceHash.size(); i++)
					coveredSpace.add(scoredSearchSpaceHash.get(i));

				scoredSearchSpaceHash.clear();
			}
		}
		else
		{
			// If not changed then reduce the range.

			// Clear the memory of the space that has been searched 
			// (as the range increment will be updated so the values may not overlap).
			coveredSpace.clear();

			// First check if the range can be reduced. 
			// If not then return false as nothing can be changed.
			changed = false;
			for (int i = 0; i < dimensions.length; i++)
			{
				if (!dimensions[i].isAtMinIncrement())
				{
					changed = true;
					break;
				}
			}
			if (changed)
			{
				// Reduce the range
				for (int i = 0; i < dimensions.length; i++)
					dimensions[i].reduce();
			}
		}

		end();

		return changed;
	}

	private static boolean changed(double[] v1, double[] v2)
	{
		if (v1.length != v2.length)
			return true;
		for (int i = 0; i < v1.length; i++)
			if (v1[i] != v2[i])
				return true;
		return false;
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
	 * Gets the current search space. This is the space that was most recently evaluated.
	 *
	 * @return the search space
	 */
	public double[][] getSearchSpace()
	{
		return searchSpace;
	}
}
