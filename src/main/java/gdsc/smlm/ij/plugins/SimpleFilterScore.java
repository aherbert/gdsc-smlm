package gdsc.smlm.ij.plugins;

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
public class SimpleFilterScore implements Comparable<SimpleFilterScore>
{
	ScoreResult r;
	double score, criteria;
	boolean criteriaPassed;
	boolean allSameType;

	public SimpleFilterScore(ScoreResult r, boolean allSameType, boolean criteriaPassed)
	{
		this.r = r;
		this.score = r.score;
		this.criteria = r.criteria;
		this.criteriaPassed = criteriaPassed;
	}

	public int compareTo(SimpleFilterScore that)
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
			if (allSameType)
			{
				// If equal criteria then if the same type get the filter with the strongest params
				return this.r.filter.weakest(that.r.filter);
			}
			else if (this.r.filter.getType().equals(that.r.filter.getType()))
			{
				return this.r.filter.weakest(that.r.filter);
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
			if (allSameType)
			{
				// If equal criteria then if the same type get the filter with the strongest params
				return this.r.filter.weakest(that.r.filter);
			}
			else if (this.r.filter.getType().equals(that.r.filter.getType()))
			{
				return this.r.filter.weakest(that.r.filter);
			}
			return 0;
		}
	}
	
	@Override
	public String toString()
	{
		// Add the score
		return String.format("%s : %.3f (%.3f)", r.filter.getName(), score, criteria);
	}
}