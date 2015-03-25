package gdsc.smlm.results;

import java.util.Collection;

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
 * Stores peak results in memory. The PeakResults interface add methods are thread safe.
 */
public class NullPeakResults extends AbstractPeakResults
{
	@Override
	public void begin()
	{
	}

	@Override
	public void add(int peak, int origX, int origY, float origValue, double chiSquared, float noise, float[] params,
			float[] paramsStdDev)
	{
	}

	@Override
	public void addAll(Collection<PeakResult> results)
	{
	}

	@Override
	public int size()
	{
		return 0;
	}

	@Override
	public void end()
	{
	}

	@Override
	public boolean isActive()
	{
		return true;
	}
}
