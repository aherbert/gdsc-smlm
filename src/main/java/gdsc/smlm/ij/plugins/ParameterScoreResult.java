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
 * Store the score from analysis of the non-filter parameters during direct filter analysis
 */
public class ParameterScoreResult
{
	final double score, criteria;
	final double[] parameters;
	final String text;

	public ParameterScoreResult(double score, double criteria, double[] parameters, String text)
	{
		this.score = score;
		this.criteria = criteria;
		this.parameters=parameters;
		this.text = text;
	}
}