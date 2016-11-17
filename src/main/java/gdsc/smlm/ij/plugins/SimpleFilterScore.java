package gdsc.smlm.ij.plugins;

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
public class SimpleFilterScore extends FilterScore
{
	final ScoreResult r;

	public SimpleFilterScore(ScoreResult r, boolean allSameType, boolean criteriaPassed)
	{
		super(r.filter, r.score, r.criteria, allSameType, criteriaPassed);
		this.r = r;
	}
}