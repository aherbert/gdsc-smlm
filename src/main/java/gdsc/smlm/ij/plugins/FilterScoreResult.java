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

import gdsc.smlm.results.filter.DirectFilter;

/**
 * Store the score from analysis of the direct filter
 */
public class FilterScoreResult
{
	final double score, criteria;
	final DirectFilter filter;
	final String text;

	public FilterScoreResult(double score, double criteria, DirectFilter filter, String text)
	{
		this.score = score;
		this.criteria = criteria;
		this.filter = filter;
		this.text = text;
	}
}