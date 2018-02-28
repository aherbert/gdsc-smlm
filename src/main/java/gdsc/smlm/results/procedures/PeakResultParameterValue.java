package gdsc.smlm.results.procedures;

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
 * Gets a parameter value from a results.
 */
public class PeakResultParameterValue implements PeakResultValue
{
	/** The parameter index. */
	public final int index;

	/**
	 * Instantiates a new peak result parameter value.
	 *
	 * @param index
	 *            the index
	 */
	public PeakResultParameterValue(int index)
	{
		this.index = index;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.procedures.PeakResultValue#getValue(gdsc.smlm.results.PeakResult)
	 */
	public float getValue(PeakResult result)
	{
		return result.getParameter(index);
	}
}