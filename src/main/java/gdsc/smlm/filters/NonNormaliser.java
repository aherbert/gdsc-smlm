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
 * Defines no normalisation of an area around a point
 */
public class NonNormaliser implements Normaliser
{
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
		return (float) sum;
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
		return sum;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.Normaliser#normalise(float[], int)
	 */
	public void normalise(float[] data, int size)
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.Normaliser#normalise(float[], float[], int)
	 */
	public void normalise(float[] data, float[] out, int size)
	{
		System.arraycopy(data, 0, out, 0, size);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.Normaliser#normalise(float[], int, int, int)
	 */
	public void normalise(float[] data, int maxx, int maxy, int border)
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.Normaliser#normalise(float[], float[], int, int, int)
	 */
	public void normalise(float[] data, float[] out, int maxx, int maxy, int border)
	{
		int xlimit = maxx - border;
		int ylimit = maxy - border;
		for (int y = border; y < ylimit; y++)
		{
			for (int x = border, i = y * maxx + border; x < xlimit; x++, i++)
			{
				out[i] = data[i];
			}
		}
	}
}