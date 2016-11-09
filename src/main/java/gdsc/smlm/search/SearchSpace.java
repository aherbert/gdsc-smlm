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

	private double[][] scoredSearchSpace;
	private ArrayList<String> scoredSearchSpaceHash = new ArrayList<String>();
	private HashSet<String> coveredSpace = new HashSet<String>();
	private StringBuilder sb = new StringBuilder();

	private static final byte ENUMERATE = (byte) 0;
	private static final byte REFINE = (byte) 1;
	// The current search mode
	private byte searchMode = ENUMERATE;

	// This introduces a dependency on another gdsc.smlm package
	private TrackProgress tracker = null;

	private boolean refinement = true;

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
	 * the bounds of the range in any dimension then the range is re-centred and the process repeated. During repeats
	 * only those points that have not yet been evaluated will be passed to the score function. If not at the bounds
	 * then the range is re-centred and the width of the range reduced.
	 * <p>
	 * The process iterates until the range cannot be reduced in size, or convergence is reached. The input dimensions
	 * are modified during the search. Use the clone(SearchDimension[]) method to create a copy.
	 * 
	 * @param <T>
	 *            the type of comparable score
	 * @param dimensions
	 *            the dimensions
	 * @param scoreFunction
	 *            the score function
	 * @param checker
	 *            the checker
	 * @return The optimum (or null)
	 */
	public <T extends Comparable<T>> SearchResult<T> search(SearchDimension[] dimensions,
			ScoreFunction<T> scoreFunction, ConvergenceChecker<T> checker)
	{
		if (dimensions == null || dimensions.length == 0)
			throw new IllegalArgumentException("Dimensions must not be empty");
		if (scoreFunction == null)
			throw new IllegalArgumentException("Score function is null");
		this.dimensions = dimensions;

		// Start with enumeration of the space
		searchMode = ENUMERATE;

		// Find the best individual
		SearchResult<T> current = findOptimum(scoreFunction, null);
		SearchResult<T> previous = null;

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
		scoredSearchSpace = null;
		scoredSearchSpaceHash.clear();
		coveredSpace.clear();

		return current;
	}

	/**
	 * Find the optimum. Create the search space using the current dimensions. Score any new point that has not
	 * previously been scored. Compare the result with the current optimum and return the best.
	 *
	 * @param <T>
	 *            the type of comparable score
	 * @param scoreFunction
	 *            the score function
	 * @param current
	 *            the current optimum
	 * @return the new optimum
	 */
	private <T extends Comparable<T>> SearchResult<T> findOptimum(ScoreFunction<T> scoreFunction,
			SearchResult<T> current)
	{
		if (!createSearchSpace(current))
			return null;

		start("Find Optimum");

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
			{
				end();
				// We have scored everything already so return the current best
				return current;
			}

			//System.out.printf("Scoring %d / %d\n", size, searchSpace.length);
			scoredSearchSpace = Arrays.copyOf(scoredSearchSpace, size);
		}

		SearchResult<T> optimum = scoreFunction.findOptimum(scoredSearchSpace);
		if (optimum != null)
		{
			//System.out.printf("Optimum = %s\n", optimum.score);
			//if (current != null)
			//	System.out.printf("Current = %s\n", current.score);
			// Replace if better
			if (optimum.compareTo(current) < 0)
			{
				current = optimum;
				//if (searchMode == REFINE)
				//	System.out.printf("Refine improved = %s\n", current.score);
			}
		}

		end();
		return current;
	}

	/**
	 * Generate hash string.
	 *
	 * @param values
	 *            the values
	 * @return the string
	 */
	private String generateHashString(double[] values)
	{
		for (int i = 0; i < values.length; i++)
			sb.append(values[i]);
		final String hash = sb.toString();
		sb.setLength(0);
		return hash;
	}

	/**
	 * Creates the search space.
	 *
	 * @param <T>
	 *            the type of comparable score
	 * @param current
	 *            the current
	 * @return true, if successful
	 */
	private <T extends Comparable<T>> boolean createSearchSpace(SearchResult<T> current)
	{
		start("Create Search Space");
		if (searchMode == REFINE && current != null)
			searchSpace = createRefineSpace(dimensions, current.point);
		else
			// Enumerate
			searchSpace = createSearchSpace(dimensions);
		end();
		return searchSpace != null;
	}

	public static double[][] createRefineSpace(SearchDimension[] dimensions, double[] point)
	{
		final double[][] dimensionValues = new double[dimensions.length][];
		for (int i = 0; i < dimensions.length; i++)
		{
			double[] values = dimensions[i].values();
			// Find the range values either side of the point value
			double v = point[i];
			double min = v;
			double max = v;
			for (int j = 0; j < values.length; j++)
			{
				if (values[j] < v)
					min = values[j];
				else if (values[j] > v)
				{
					max = values[j];
					break;
				}
			}
			// Create a sequence from min to max.
			// Add option to configure the number of steps?
			dimensionValues[i] = dimensions[i].enumerate(min, max, 100);
		}
		return createRefineSpace(dimensionValues, point);
	}

	/**
	 * Creates the refine space.
	 *
	 * @param dimensionValues
	 *            the dimension values
	 * @param point
	 *            the point
	 * @return the double[][]
	 */
	private static double[][] createRefineSpace(double[][] dimensionValues, double[] point)
	{
		// Get the values
		int combinations = 0;
		for (int i = 0; i < dimensionValues.length; i++)
		{
			combinations += dimensionValues[i].length;
		}

		// This will be a list of points enumerating the entire range
		// of the dimensions
		final double[][] searchSpace = new double[combinations][];

		// Start with the min value in each dimension
		int n = 0;
		// Enumerate each dimension
		for (int i = 0; i < dimensionValues.length; i++)
		{
			// Values to iterate over for this dimension
			final double[] v1 = dimensionValues[i];

			for (int k = 0; k < v1.length; k++)
			{
				// Create a new point with an updated value for this dimension
				final double[] v3 = point.clone();
				v3[i] = v1[k];
				searchSpace[n++] = v3;
			}
		}

		return searchSpace;
	}

	/**
	 * Creates the search space.
	 *
	 * @param dimensions
	 *            the dimensions
	 * @return the double[][]
	 */
	public static double[][] createSearchSpace(SearchDimension[] dimensions)
	{
		// Get the values
		final double[][] dimensionValues = new double[dimensions.length][];
		for (int i = 0; i < dimensions.length; i++)
		{
			dimensionValues[i] = dimensions[i].values();
			//System.out.printf(" [%d] %.3f-%.3f", i, dimensions[i].getLower(), dimensions[i].getUpper());
		}
		//System.out.println();
		return createSearchSpace(dimensionValues);
	}

	/**
	 * Creates the search space.
	 *
	 * @param dimensionValues
	 *            the dimension values
	 * @return the double[][]
	 */
	private static double[][] createSearchSpace(double[][] dimensionValues)
	{
		// Get the values
		int combinations = 1;
		final double[] v = new double[dimensionValues.length];
		for (int i = 0; i < dimensionValues.length; i++)
		{
			combinations *= dimensionValues[i].length;
			v[i] = dimensionValues[i][0];
		}

		// This will be a list of points enumerating the entire range
		// of the dimensions
		final double[][] searchSpace = new double[combinations][];

		// Start with the min value in each dimension
		searchSpace[0] = v;
		int n = 1;
		try
		{
			// Enumerate each dimension
			for (int i = 0; i < dimensionValues.length; i++)
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
						v3[i] = v1[k];
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

	/**
	 * Update search space. Re-centre the dimension on the current optimum. If the current optimum is at the bounds of
	 * the dimensions then check if the bounds change when re-centring. If the bound change then the search must
	 * continue with the current dimension ranges. If the bounds do not change then the dimension ranges can be reduced.
	 *
	 * @param current
	 *            the current optimum
	 * @return true, if successful
	 */
	private boolean updateSearchSpace(SearchResult<?> current)
	{
		if (current == null)
			return false;

		start("Update search space");
		boolean changed = false;

		if (searchMode == REFINE)
		{
			// During refinement we will not repeat the same point if the optimum moves
			coveredSpace.clear();

			// Move to the centre using the current optimum
			final double[] p = current.point;
			for (int i = 0; i < dimensions.length; i++)
			{
				if (p[i] != dimensions[i].getCentre())
				{
					changed = true;
					dimensions[i].setCentre(p[i]);
				}
			}

			if (!changed)
			{
				// Switch back to enumeration of a smaller range
				searchMode = ENUMERATE;
				changed = reduceRange();
			}
		}
		else
		{
			// searchMode == ENUMERATION

			// Process each dimension
			final double[] p = current.point;
			for (int i = 0; i < dimensions.length; i++)
			{
				// Check if at the bounds of the dimension values
				final boolean atBounds = dimensions[i].isAtBounds(p[i]);
				final double[] values = (atBounds) ? dimensions[i].values() : null;

				// Move to the centre using the current optimum
				dimensions[i].setCentre(p[i]);

				// If at bounds then check if the bounds have changed due to the move
				if (atBounds)
				{
					if (changed(values, dimensions[i].values()))
						changed = true;
				}
			}

			// Q. Always reduce the range ... ?
			//changed = false;

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
				// Clear the memory of the space that has been searched 
				// (as the search space is about to be altered so the values may not overlap).
				coveredSpace.clear();

				// Optionally switch to refinement
				if (refinement)
				{
					searchMode = REFINE;
					changed = true;
				}
				else
				{

					changed = reduceRange();
				}
			}
		}

		end();

		return changed;
	}

	/**
	 * Reduce range.
	 *
	 * @return true, if successful
	 */
	private boolean reduceRange()
	{
		// First check if the range can be reduced. 
		// If not then return false as nothing can be changed.
		boolean reduced = false;
		for (int i = 0; i < dimensions.length; i++)
		{
			if (!dimensions[i].isAtMinIncrement())
			{
				reduced = true;
				break;
			}
		}
		if (reduced)
		{
			// Reduce the range
			for (int i = 0; i < dimensions.length; i++)
				dimensions[i].reduce();
		}
		return reduced;
	}

	/**
	 * Check if the two arrays are different.
	 *
	 * @param v1
	 *            the v 1
	 * @param v2
	 *            the v 2
	 * @return true, if different
	 */
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

	/**
	 * Send a start message to the tracker
	 * 
	 * @param stage
	 *            The stage that has started
	 */
	private void start(String stage)
	{
		if (tracker != null)
		{
			tracker.status(stage + " [%d]", iteration);
			tracker.progress(0);
		}
	}

	/**
	 * Send an end signal to the tracker
	 */
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

	/**
	 * Checks if refinement within the search space is enabled.
	 *
	 * @return true, if is refinement
	 */
	public boolean isRefinement()
	{
		return refinement;
	}

	/**
	 * Set to true to allow refinement within the current search space.
	 *
	 * @param refinement
	 *            the new refinement value
	 */
	public void setRefinement(boolean refinement)
	{
		this.refinement = refinement;
	}
}
