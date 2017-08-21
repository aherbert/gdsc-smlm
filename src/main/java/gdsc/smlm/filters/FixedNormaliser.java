package gdsc.smlm.filters;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Defines normalisation of an area around a point using a fixed normalisation
 */
public class FixedNormaliser implements Normaliser
{
	/** The normalisation. */
	public final float normalisation;

	/**
	 * Instantiates a new fixed normaliser.
	 *
	 * @param normalisation the normalisation
	 */
	public FixedNormaliser(float normalisation)
	{
		this.normalisation = normalisation;
	}

	/**
	 * Normalise the sum.
	 *
	 * @param sum
	 *            the sum
	 * @param index
	 *            the index
	 * @return the normalised value
	 */
	public float normalise(double sum, int index)
	{
		return (float) (sum / normalisation);
	}

	/**
	 * Normalise the sum.
	 *
	 * @param sum
	 *            the sum
	 * @param index
	 *            the index
	 * @return the normalised value
	 */
	public float normalise(float sum, int index)
	{
		return sum / normalisation;
	}
}