package gdsc.smlm.results.procedures;

import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResultValue;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2018 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Class to get min/max value in a set of results
 */
public class MinMaxResultProcedure implements PeakResultProcedure
{
	private float min, max;
	private final PeakResultValue value;

	/**
	 * Instantiates a new min max result procedure.
	 *
	 * @param results
	 *            the results
	 * @param value
	 *            the value
	 * @throws IllegalStateException
	 *             If the results are empty
	 */
	public MinMaxResultProcedure(MemoryPeakResults results, PeakResultValue value) throws IllegalStateException
	{
		this.value = value;
		min = max = value.getValue(results.getFirst());
		results.forEach(this);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.procedures.PeakResultProcedure#execute(gdsc.smlm.results.PeakResult)
	 */
	public void execute(PeakResult peakResult)
	{
		float v = value.getValue(peakResult);
		if (min > v)
			min = v;
		else if (max < v)
			max = v;
	}

	/**
	 * Gets the minimum.
	 *
	 * @return the minimum
	 */
	public float getMinimum()
	{
		return min;
	}

	/**
	 * Gets the maximum.
	 *
	 * @return the maximum
	 */
	public float getMaximum()
	{
		return max;
	}
}
