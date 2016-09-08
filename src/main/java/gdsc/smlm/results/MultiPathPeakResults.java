package gdsc.smlm.results;

import java.util.Collection;

import gdsc.core.utils.NotImplementedException;

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
 * Provides storage for MultiPathPeakResults.
 * <p>
 * Does nothing for any of the inherited abstract methods for the PeakResults. The base class is extended to support the
 * storage of the configuration and calibration.
 */
public class MultiPathPeakResults extends AbstractPeakResults
{
	// TODO - Store the results of fitting multiple pathways.
	// Then build a MultiPathFilter that can score them. 
	// This can be passed in to the FitWorker when computing the benchmark fit results.
	
	
	/**
	 * This is not implemented
	 * 
	 * @throws gdsc.core.utils.NotImplementedException
	 * @see gdsc.smlm.results.AbstractPeakResults#begin()
	 */
	@Override
	public void begin()
	{
		throw new NotImplementedException();
	}

	/**
	 * This is not implemented
	 * 
	 * @throws gdsc.core.utils.NotImplementedException
	 * @see gdsc.smlm.results.AbstractPeakResults#add(int, int, int, float, double, float, float[], float[])
	 */
	@Override
	public void add(int peak, int origX, int origY, float origValue, double chiSquared, float noise, float[] params,
			float[] paramsStdDev)
	{
	}

	/**
	 * This is not implemented
	 * 
	 * @throws gdsc.core.utils.NotImplementedException
	 * @see gdsc.smlm.results.AbstractPeakResults#addAll(java.util.Collection)
	 */
	@Override
	public void addAll(Collection<PeakResult> results)
	{
	}

	@Override
	public int size()
	{
		return 0;
	}

	/**
	 * This is not implemented
	 * 
	 * @throws gdsc.core.utils.NotImplementedException
	 * @see gdsc.smlm.results.AbstractPeakResults#end()
	 */
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
