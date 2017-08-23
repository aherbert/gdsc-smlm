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
 * Defines normalisation of an area around a point
 */
public interface Normaliser
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
	public float normalise(double sum, int index);

	/**
	 * Normalise the sum.
	 *
	 * @param sum
	 *            the sum
	 * @param index
	 *            the index
	 * @return the normalised value
	 */
	public float normalise(float sum, int index);

	/**
	 * Normalise the sum of the data (in-place) for all indices.
	 *
	 * @param data
	 *            the data
	 * @param size
	 *            the size of the data
	 */
	public void normalise(float[] data, int size);

	/**
	 * Normalise the sum of the data for all indices to the out data.
	 *
	 * @param data
	 *            the data
	 * @param out
	 *            the out data
	 * @param size
	 *            the size of the data
	 */
	public void normalise(float[] data, float[] out, int size);

	/**
	 * Normalise the sum of the data (in-place) ignoring the border.
	 *
	 * @param data
	 *            the data
	 * @param maxx
	 *            the maxx
	 * @param maxy
	 *            the maxy
	 * @param border
	 *            the border
	 */
	public void normalise(float[] data, int maxx, int maxy, int border);

	/**
	 * Normalise the sum of the data ignoring the border to the out data.
	 *
	 * @param data
	 *            the data
	 * @param maxx
	 *            the maxx
	 * @param maxy
	 *            the maxy
	 * @param border
	 *            the border
	 */
	public void normalise(float[] data, float[] out, int maxx, int maxy, int border);

	/**
	 * Normalise the sum of the data for all indices to the out data.
	 *
	 * @param data
	 *            the data
	 * @param out
	 *            the out data
	 * @param size
	 *            the size of the data
	 */
	public void normalise(double[] data, float[] out, int size);

	/**
	 * Normalise the sum of the data ignoring the border to the out data.
	 *
	 * @param data
	 *            the data
	 * @param maxx
	 *            the maxx
	 * @param maxy
	 *            the maxy
	 * @param border
	 *            the border
	 */
	public void normalise(double[] data, float[] out, int maxx, int maxy, int border);
}