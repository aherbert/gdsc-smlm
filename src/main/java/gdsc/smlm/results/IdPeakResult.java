package gdsc.smlm.results;

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
 * Specifies a peak fitting result with an Id
 */
public class IdPeakResult extends PeakResult
{
	private final int id;

	/**
	 * Instantiates a new id peak result.
	 *
	 * @param frame
	 *            the frame
	 * @param origX
	 *            the orig X
	 * @param origY
	 *            the orig Y
	 * @param origValue
	 *            the orig value
	 * @param error
	 *            the error
	 * @param noise
	 *            the noise
	 * @param params
	 *            the params
	 * @param paramsStdDev
	 *            the params std dev
	 * @param id
	 *            the id
	 */
	public IdPeakResult(int frame, int origX, int origY, float origValue, double error, float noise,
			float[] params, float[] paramsStdDev, int id)
	{
		super(frame, origX, origY, origValue, error, noise, params, paramsStdDev);
		this.id = id;
	}

	/**
	 * Instantiates a new id peak result.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param sd
	 *            the sd
	 * @param signal
	 *            the signal
	 * @param id
	 *            the id
	 */
	public IdPeakResult(float x, float y, float sd, float signal, int id)
	{
		super(x, y, sd, signal);
		this.id = id;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResult#hasId()
	 */
	@Override
	public boolean hasId()
	{
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResult#getId()
	 */
	@Override
	public int getId()
	{
		return id;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResult#compareTo(gdsc.smlm.results.PeakResult)
	 */
	public int compareTo(PeakResult o)
	{
		// Sort by peak number: Ascending
		if (getFrame() == o.getFrame())
		{
			return getId() - o.getId();
		}
		return getFrame() - o.getFrame();
	}
}
