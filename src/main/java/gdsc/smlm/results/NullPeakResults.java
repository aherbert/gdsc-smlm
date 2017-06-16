package gdsc.smlm.results;

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
 * Does nothing for any of the PeakResults methods
 */
public class NullPeakResults extends AbstractPeakResults implements ThreadSafePeakResults
{
	public void begin()
	{
	}

	public void add(int peak, int origX, int origY, float origValue, double error, float noise, float[] params,
			float[] paramsStdDev)
	{
	}
	
	public void add(PeakResult result)
	{
	}

	public void addAll(PeakResult[] results)
	{
	}

	public int size()
	{
		return 0;
	}

	public void end()
	{
	}

	public boolean isActive()
	{
		return true;
	}
}
