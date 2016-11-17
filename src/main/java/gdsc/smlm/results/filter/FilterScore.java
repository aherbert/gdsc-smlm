package gdsc.smlm.results.filter;

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
 * Store the filter score used in benchmarking
 */
public class FilterScore implements Comparable<FilterScore>
{
	final public Filter filter;
	final public double score, criteria;
	final public boolean criteriaPassed;
	final public boolean allSameType;

	public FilterScore(Filter filter, double score, double criteria, boolean allSameType, boolean criteriaPassed)
	{
		this.filter = filter;
		this.score = score;
		this.criteria = criteria;
		this.allSameType = allSameType;
		this.criteriaPassed = criteriaPassed;
	}

	public int compareTo(FilterScore that)
	{
		if (that == null)
			return -1;

		//		if (this.criteriaPassed && !that.criteriaPassed)
		//			return -1;
		//		if (that.criteriaPassed && !this.criteriaPassed)
		//			return 1;

		if (this.criteriaPassed)
		{
			// Must pass criteria first
			if (!that.criteriaPassed)
				return -1;

			// Sort by the score
			if (this.score > that.score)
				return -1;
			if (this.score < that.score)
				return 1;
			if (this.criteria > that.criteria)
				return -1;
			if (this.criteria < that.criteria)
				return 1;
			// If the same type then compare the parameters 
			if (allSameType)
			{
				return compareParameters(that);
			}
			else if (this.filter.getType().equals(that.filter.getType()))
			{
				return compareParameters(that);
			}
			return 0;
		}
		else
		{
			// Must pass criteria first
			if (that.criteriaPassed)
				return 1;

			// Sort by how close we are to passing the criteria
			if (this.criteria > that.criteria)
				return -1;
			if (this.criteria < that.criteria)
				return 1;
			if (this.score > that.score)
				return -1;
			if (this.score < that.score)
				return 1;
			// If the same type then compare the parameters 
			if (allSameType)
			{
				return compareParameters(that);
			}
			else if (this.filter.getType().equals(that.filter.getType()))
			{
				return compareParameters(that);
			}
			return 0;
		}
	}

	private int compareParameters(FilterScore that)
	{
		return compare(this.filter, that.filter);
	}

	/**
	 * Compare the two filters. Get the filter with the strongest parameters.
	 *
	 * @param f1
	 *            filter 1
	 * @param f2
	 *            filter 2
	 * @return the comparison score
	 */
	public static int compare(Filter f1, Filter f2)
	{
		// Get the filter with the weakest params
		//return f1.weakestUnsafe(f2);

		// Get the filter with the strongest params
		return f2.weakestUnsafe(f1);
	}

	@Override
	public String toString()
	{
		// Add the score
		return String.format("%s : %.3f (%.3f)", filter.getName(), score, criteria);
	}
}