package gdsc.smlm.results.predicates;

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
 * Accept results with a value within a min/max range
 */
public class MinMaxPeakResultPredicate implements PeakResultPredicate
{
	/** The min of the value range. */
	public final float min;

	/** The max of the value range. */
	public final float max;

	/** The value. */
	public final PeakResultValue value;

	/**
	 * Instantiates a new min max peak result predicate.
	 *
	 * @param min
	 *            the min of the value range
	 * @param max
	 *            the max of the value range
	 * @param value
	 *            the value
	 * @throws IllegalArgumentException
	 *             If the min/max are NaN. If min is greater than max. If value is null.
	 */
	public MinMaxPeakResultPredicate(float min, float max, PeakResultValue value) throws IllegalArgumentException
	{
		if (Float.isNaN(min) || Float.isNaN(max) || min > max)
			throw new IllegalArgumentException("Min/Max range is invalid");
		if (value == null)
			throw new IllegalArgumentException("No value function");
		this.min = min;
		this.max = max;
		this.value = value;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.procedures.PeakResultPredicate#test(gdsc.smlm.results.PeakResult)
	 */
	public boolean test(PeakResult t)
	{
		final float v = value.getValue(t);
		return v >= min && v <= max;
	}
}