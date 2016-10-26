package gdsc.smlm.filters;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2015 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Identifies candidate spots (local maxima) in an image using non-maximum suppression.
 */
public abstract class MaximaSpotFilter extends SpotFilter
{
	private final int search;
	private final int border;
	private NonMaximumSuppression nms;
	private DataProcessor scoreDataProcessor = null;
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
	protected Spot[] find(final float[] data, final int width, final int height)
	{
		data2 = preprocessData(data, width, height);

		//gdsc.core.ij.Utils.display("Spot Filter", new FloatProcessor(width, height, data2));

		final int[] maxIndices = getMaxima(data2, width, height);
		if (maxIndices.length == 0)
			return null;

		// Use a data processor to generate a ranking score
		final float[] scoreData;
		if (scoreDataProcessor != null)
		{
			scoreData = data.clone();
			scoreDataProcessor.process(scoreData, width, height);
		}
		else
		{
			scoreData = data2;
		}

		final Spot[] spots = new Spot[maxIndices.length];
		for (int n = 0; n < maxIndices.length; n++)
		{
			final int y = maxIndices[n] / width;
			final int x = maxIndices[n] % width;
			final float intensity = data2[maxIndices[n]];
			final float score = scoreData[maxIndices[n]];
			spots[n] = new Spot(x, y, intensity, score);
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
	 * Pre-process the data before finding local maxima
	 * 
	 * @param data
	 * @param width
	 * @param height
	 * @return The pre-processed data
	 */
	public abstract float[] preprocessData(final float[] data, final int width, final int height);

	/**
	 * Find the indices of the maxima using the currently configured parameters
	 * <p>
	 * Data must be arranged in yx block order, i.e. height rows of width.
	 * 
	 * @param data
	 * @param width
	 * @param height
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
		ArrayList<String> list = new ArrayList<String>();
		list.add("search = " + search);
		list.add("border = " + border);
		if (scoreDataProcessor != null)
			list.add("score = " + scoreDataProcessor.getDescription());
		return list;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public MaximaSpotFilter clone()
	{
		MaximaSpotFilter f = (MaximaSpotFilter) super.clone();
		// Ensure the object is duplicated and not passed by reference.
		f.nms = nms.clone();
		f.data2 = null;
		if (scoreDataProcessor != null)
			f.scoreDataProcessor = scoreDataProcessor.clone();
		return f;
	}

	/**
	 * @return the score data processor
	 */
	public DataProcessor getScoreDataProcessor()
	{
		return scoreDataProcessor;
	}

	/**
	 * Set a data processor to use to generate a score for each spot that is used to rank them. If null then the score
	 * is the same as the intensity score of the spot.
	 * <p>
	 * This allows the spots to be ranked independently of the intensity estimate that will be used for fitting.
	 * 
	 * @param scoreDataProcessor
	 *            the score data processor
	 */
	public void setScoreDataProcessor(DataProcessor scoreDataProcessor)
	{
		this.scoreDataProcessor = scoreDataProcessor;
	}
}