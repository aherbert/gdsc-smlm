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
	 * Instantiates a new peak result.
	 *
	 * @param frame
	 *            the frame
	 * @param origX
	 *            the original X position
	 * @param origY
	 *            the original Y position
	 * @param origValue
	 *            the original value
	 * @param error
	 *            the error
	 * @param noise
	 *            the noise
	 * @param meanIntensity
	 *            the mean intensity
	 * @param params
	 *            the params (must not be null and must have at least {@value #STANDARD_PARAMETERS} parameters)
	 * @param paramsStdDev
	 *            the params standard deviations (if not null must match the length of the {@link #params} array)
	 * @param id
	 *            the id
	 * @throws IllegalArgumentException
	 *             the illegal argument exception if the parameters are invalid
	 */
	public IdPeakResult(int frame, int origX, int origY, float origValue, double error, float noise,
			float meanIntensity, float[] params, float[] paramsStdDev, int id) throws IllegalArgumentException
	{
		super(frame, origX, origY, origValue, error, noise, meanIntensity, params, paramsStdDev);
		this.id = id;
	}

	/**
	 * Instantiates a new id peak result.
	 *
	 * @param frame
	 *            the frame
	 * @param x
	 *            the x position
	 * @param y
	 *            the y position
	 * @param intensity
	 *            the intensity
	 * @param id
	 *            the id
	 */
	public IdPeakResult(int frame, float x, float y, float intensity, int id)
	{
		super(frame, x, y, intensity);
		this.id = id;
	}

	/**
	 * Instantiates a new id peak result.
	 *
	 * @param x
	 *            the x position
	 * @param y
	 *            the y position
	 * @param intensity
	 *            the intensity
	 * @param id
	 *            the id
	 */
	public IdPeakResult(float x, float y, float intensity, int id)
	{
		super(x, y, intensity);
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
		if (getFrame() < o.getFrame())
			return -1;
		if (getFrame() > o.getFrame())
			return 1;
		// Sort by ID height: Ascending
		if (getId() > o.getId())
			return -1;
		if (getId() < o.getId())
			return 1;
		return 0;
	}
}
