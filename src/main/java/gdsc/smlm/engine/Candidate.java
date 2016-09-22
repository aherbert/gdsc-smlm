package gdsc.smlm.engine;

import gdsc.smlm.filters.Spot;

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
 * Extends the Spot class with fields used during fitting
 */
public class Candidate extends Spot
{
	/**
	 * The index
	 */
	final public int index;
	/**
	 * Flag to indicate if the candidate has been fit
	 */
	public boolean fit = false;

	/**
	 * Instantiates a new candidate.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param intensity
	 *            the intensity
	 * @param score
	 *            the score
	 * @param index
	 *            the index
	 */
	public Candidate(int x, int y, float intensity, float score, int index)
	{
		super(x, y, intensity, score);
		this.index = index;
	}

	/**
	 * Instantiates a new candidate.
	 *
	 * @param spot the spot
	 * @param index the index
	 */
	public Candidate(Spot spot, int index)
	{
		super(spot.x, spot.y, spot.intensity, spot.getScore());
		this.index = index;
	}
}