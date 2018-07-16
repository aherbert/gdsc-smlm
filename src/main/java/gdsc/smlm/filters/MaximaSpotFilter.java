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

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

import gdsc.core.filters.NonMaximumSuppression;

/**
 * Identifies candidate spots (local maxima) in an image using non-maximum suppression.
 */
public abstract class MaximaSpotFilter extends SpotFilter
{
	private final int search;
	private final int border;
	private NonMaximumSuppression nms;
	private float[] data2 = null;

	/**
	 * Create the spot filter
	 *
	 * @param search
	 *            The search width for non-maximum suppression
	 * @param border
	 *            The border to ignore for maxima
	 * @throws IllegalArgumentException
	 *             if search is below 1 or border is below zero
	 */
	public MaximaSpotFilter(int search, int border)
	{
		if (search < 1)
			throw new IllegalArgumentException("Search width must be 1 or above");
		if (border < 0)
			throw new IllegalArgumentException("Border must be 0 or above");
		this.search = search;
		this.border = border;
		nms = new NonMaximumSuppression();
		// Do a neighbour check when using a low block size
		nms.setNeighbourCheck(search < 3);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.filters.SpotFilter#find(float[], int, int)
	 */
	@Override
	protected Spot[] find(final float[] data, final int width, final int height)
	{
		data2 = preprocessData(data, width, height);

		//gdsc.core.ij.Utils.display("Spot Filter", new FloatProcessor(width, height, data2));

		final int[] maxIndices = getMaxima(data2, width, height);
		if (maxIndices.length == 0)
			return null;

		final Spot[] spots = new Spot[maxIndices.length];
		for (int n = 0; n < maxIndices.length; n++)
		{
			final int y = maxIndices[n] / width;
			final int x = maxIndices[n] % width;
			final float intensity = data2[maxIndices[n]];
			spots[n] = new Spot(x, y, intensity);
		}
		return spots;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.filters.SpotFilter#getPreprocessedData()
	 */
	@Override
	public float[] getPreprocessedData()
	{
		return data2;
	}

	/**
	 * Pre-process the data before finding local maxima.
	 *
	 * @param data
	 *            the data
	 * @param width
	 *            the width
	 * @param height
	 *            the height
	 * @return The pre-processed data
	 */
	public abstract float[] preprocessData(final float[] data, final int width, final int height);

	/**
	 * Find the indices of the maxima using the currently configured parameters
	 * <p>
	 * Data must be arranged in yx block order, i.e. height rows of width.
	 *
	 * @param data
	 *            the data
	 * @param width
	 *            the width
	 * @param height
	 *            the height
	 * @return Indices of the maxima
	 */
	protected int[] getMaxima(float[] data, int width, int height)
	{
		// Check upper limits are safe
		final int n = FastMath.min(search, FastMath.min(width, height));
		final int border = FastMath.min(this.border, FastMath.min(width, height) / 2);
		return nms.blockFindInternal(data, width, height, n, border);
	}

	/**
	 * @return the search width for maxima (maximum must be the highest point in a 2n+1 region)
	 */
	public int getSearch()
	{
		return search;
	}

	/**
	 * @return the border at the edge to ignore for maxima
	 */
	public int getBorder()
	{
		return border;
	}

	@Override
	public List<String> getParameters()
	{
		ArrayList<String> list = new ArrayList<>();
		list.add("search = " + search);
		list.add("border = " + border);
		return list;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see java.lang.Object#clone()
	 */
	@Override
	public MaximaSpotFilter clone()
	{
		MaximaSpotFilter f = (MaximaSpotFilter) super.clone();
		// Ensure the object is duplicated and not passed by reference.
		f.nms = nms.clone();
		f.data2 = null;
		return f;
	}
}
