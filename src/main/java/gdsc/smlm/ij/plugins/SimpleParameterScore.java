package gdsc.smlm.ij.plugins;

import gdsc.smlm.results.filter.DirectFilter;
import gdsc.smlm.results.filter.FilterScore;

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
public class SimpleParameterScore extends FilterScore
{
	final ParameterScoreResult r;

	public SimpleParameterScore(DirectFilter filter, ParameterScoreResult r, boolean criteriaPassed)
	{
		super(filter, r.score, r.criteria, true, criteriaPassed);
		this.r = r;
	}

	@Override
	protected int compareParameters(FilterScore that)
	{
		// Compare the parameters and return the strongest, those most likely to restrict the output
		
		// 0 = failCount
		// 1 = residudalsThreshold
		// 2 = duplicateDistance
		double[] p1 = this.r.parameters;
		double[] p2 = ((SimpleParameterScore)that).r.parameters;
		
		// Lowest fail count
		if (p1[0] < p2[0])
			return -1;
		if (p1[0] > p2[0])
			return 1;
		
		// Lowest duplicate distance 
		if (p1[2] < p2[2])
			return -1;
		if (p1[2] > p2[2])
			return 1;

		// Highest residuals threshold 
		if (p1[2] > p2[2])
			return -1;
		if (p1[2] < p2[2])
			return 1;
		
		return 0;
	}
}