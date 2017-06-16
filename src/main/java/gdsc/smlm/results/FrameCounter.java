package gdsc.smlm.results;

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
 * Class to count per frame.
 */
public class FrameCounter extends Counter
{
	private int current, previous;

	/**
	 * Instantiates a new frame counter with the current value of -1.
	 */
	public FrameCounter()
	{
		this(-1);
	}

	/**
	 * Instantiates a new frame counter with the given frame.
	 *
	 * @param frame
	 *            the current frame
	 */
	public FrameCounter(int frame)
	{
		current = previous = frame;
	}

	/**
	 * Advance the frame. If this is a new frame then the count is reset.
	 *
	 * @param frame
	 *            the frame
	 * @return true if the frame has changed
	 */
	public boolean advanceAndReset(int frame)
	{
		if (current != frame)
		{
			previous = current;
			current = frame;
			reset();
			return true;
		}
		return false;
	}
	
	/**
	 * Advance the frame.
	 *
	 * @param frame
	 *            the frame
	 * @return true if the frame has changed
	 */
	public boolean advance(int frame)
	{
		if (current != frame)
		{
			previous = current;
			current = frame;
			return true;
		}
		return false;
	}

	/**
	 * Get the current frame
	 *
	 * @return the current frame
	 */
	public int currentFrame()
	{
		return current;
	}

	/**
	 * Gets the previous frame.
	 * <p>
	 * Note: If the counter has never been advanced to a new frame then the previous will be the same
	 * as the current.
	 *
	 * @return the previous frame
	 */
	public int previousFrame()
	{
		return previous;
	}
}