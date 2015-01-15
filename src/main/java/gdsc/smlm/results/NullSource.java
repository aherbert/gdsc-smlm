package gdsc.smlm.results;

import java.awt.Rectangle;

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
 * Represent a named results source. Does not support data provision.
 */
public class NullSource extends ImageSource
{
	/**
	 * Create a new image source
	 * 
	 * @param name
	 */
	public NullSource(String name)
	{
		super(name);
	}

	@Override
	protected boolean openSource()
	{
		return false;
	}

	@Override
	public void close()
	{
		// Nothing to do
	}

	@Override
	protected float[] nextFrame(Rectangle bounds)
	{
		return null;
	}

	@Override
	protected float[] getFrame(int frame, Rectangle bounds)
	{
		return null;
	}

	@Override
	public boolean isValid(int frame)
	{
		return false;
	}
}
