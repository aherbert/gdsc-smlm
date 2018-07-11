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
package gdsc.smlm.filters;

/**
 * Defines normalisation of an area around a point using a per-pixel normalisation
 */
public class PerPixelNormaliser implements Normaliser
{
	/** The normalisation. */
	public final float[] normalisation;

	/**
	 * Instantiates a new fixed normaliser.
	 *
	 * @param normalisation
	 *            the normalisation
	 */
	public PerPixelNormaliser(float[] normalisation)
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
	@Override
	public float normalise(double sum, int index)
	{
		return (float) (sum / normalisation[index]);
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
	@Override
	public float normalise(float sum, int index)
	{
		return sum / normalisation[index];
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.filters.Normaliser#normalise(float[], int)
	 */
	@Override
	public void normalise(float[] data, int size)
	{
		for (int i = 0; i < size; i++)
			data[i] /= normalisation[i];
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.filters.Normaliser#normalise(float[], float[], int)
	 */
	@Override
	public void normalise(float[] data, float[] out, int size)
	{
		for (int i = 0; i < size; i++)
			out[i] = data[i] / normalisation[i];
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.filters.Normaliser#normalise(float[], int, int, int)
	 */
	@Override
	public void normalise(float[] data, int maxx, int maxy, int border)
	{
		int xlimit = maxx - border;
		int ylimit = maxy - border;
		for (int y = border; y < ylimit; y++)
		{
			for (int x = border, i = y * maxx + border; x < xlimit; x++, i++)
			{
				data[i] /= normalisation[i];
			}
		}
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.filters.Normaliser#normalise(float[], float[], int, int, int)
	 */
	@Override
	public void normalise(float[] data, float[] out, int maxx, int maxy, int border)
	{
		int xlimit = maxx - border;
		int ylimit = maxy - border;
		for (int y = border; y < ylimit; y++)
		{
			for (int x = border, i = y * maxx + border; x < xlimit; x++, i++)
			{
				out[i] = data[i] / normalisation[i];
			}
		}
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.filters.Normaliser#normalise(double[], float[], int)
	 */
	@Override
	public void normalise(double[] data, float[] out, int size)
	{
		for (int i = 0; i < size; i++)
			out[i] = (float) (data[i] / normalisation[i]);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.filters.Normaliser#normalise(double[], float[], int, int, int)
	 */
	@Override
	public void normalise(double[] data, float[] out, int maxx, int maxy, int border)
	{
		int xlimit = maxx - border;
		int ylimit = maxy - border;
		for (int y = border; y < ylimit; y++)
		{
			for (int x = border, i = y * maxx + border; x < xlimit; x++, i++)
			{
				out[i] = (float) (data[i] / normalisation[i]);
			}
		}
	}
}
