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
public class ExtendedPeakResult extends PeakResult
{
	private int endFrame, id;

	public ExtendedPeakResult(int startFrame, int origX, int origY, float origValue, double error, float noise,
			float[] params, float[] paramsStdDev, int endFrame, int id)
	{
		super(startFrame, origX, origY, origValue, error, noise, params, paramsStdDev);
		this.endFrame = endFrame;
		this.id = id;
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
	
	/* (non-Javadoc)
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
		if (peak == o.peak)
		{
			// Sort by peak end number: Descending
			if (endFrame == o.getEndFrame())
			{
				// Sort by peak height: Descending
				if (params[1] > o.params[1])
					return -1;
				if (params[1] < o.params[1])
					return 1;
				// Finally by Id
				return getId() - o.getId();
			}
			return endFrame - o.getEndFrame();
		}
		return peak - o.peak;
	}
}
