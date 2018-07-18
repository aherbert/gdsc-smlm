/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package uk.ac.sussex.gdsc.smlm.filters;

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
	 * @param out
	 *            the out data
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
	 * @param out
	 *            the out data
	 * @param maxx
	 *            the maxx
	 * @param maxy
	 *            the maxy
	 * @param border
	 *            the border
	 */
	public void normalise(double[] data, float[] out, int maxx, int maxy, int border);
}
