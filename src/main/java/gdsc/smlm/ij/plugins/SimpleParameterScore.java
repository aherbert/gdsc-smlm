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
		// TODO - Compare the parameters and return the strongest
		return 0;
	}
}