package gdsc.smlm.results;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Specifies a peak fitting result that spans multiple frames
 */
public class ExtendedPeakResult extends IdPeakResult
{
	private int endFrame;

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
	 * @param params
	 *            the params (must not be null and must have at least {@value #STANDARD_PARAMETERS} parameters)
	 * @param paramsStdDev
	 *            the params standard deviations (if not null must match the length of the {@link #params} array)
	 * @throws IllegalArgumentException
	 *             the illegal argument exception if the parameters are invalid
	 */
	public ExtendedPeakResult(int startFrame, int origX, int origY, float origValue, double error, float noise,
			float[] params, float[] paramsStdDev, int endFrame, int id) throws IllegalArgumentException
	{
		super(startFrame, origX, origY, origValue, error, noise, params, paramsStdDev, id);
		setEndFrame(endFrame);
	}

	/**
	 * Simple constructor to create a result with frame, location, width, strength, and id.
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
	public ExtendedPeakResult(int frame, float x, float y, float intensity, int id)
	{
		super(frame, x, y, intensity, id);
	}

	/**
	 * Simple constructor to create a result with location, width, strength, and id.
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
	public ExtendedPeakResult(float x, float y, float intensity, int id)
	{
		super(x, y, intensity, id);
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResult#hasEndFrame()
	 */
	@Override
	public boolean hasEndFrame()
	{
		return endFrame > getFrame();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResult#getEndFrame()
	 */
	@Override
	public int getEndFrame()
	{
		return endFrame;
	}

	public void setEndFrame(int endFrame)
	{
		// Ensure that the end frame is valid
		this.endFrame = Math.max(endFrame, super.getFrame());
	}

	@Override
	public void setFrame(int frame)
	{
		// Set the new start frame
		super.setFrame(frame);
		// Validate the current end frame
		setEndFrame(endFrame);
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
			// Sort by peak end number: Ascending
			if (endFrame == o.getEndFrame())
			{
				// Sort by peak height: Descending
				if (getSignal() > o.getSignal())
					return -1;
				if (getSignal() < o.getSignal())
					return 1;
				// Finally by Id
				return getId() - o.getId();
			}
			return endFrame - o.getEndFrame();
		}
		return getFrame() - o.getFrame();
	}
}
