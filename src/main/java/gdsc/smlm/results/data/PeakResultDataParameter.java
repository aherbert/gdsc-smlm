package gdsc.smlm.results.data;

import gdsc.smlm.results.PeakResult;

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
 * Gets a parameter data value from a result.
 */
public class PeakResultDataParameter extends PeakResultDataFloat
{
	/** The parameter index. */
	public final int index;

	/**
	 * Instantiates a new peak result parameter value.
	 *
	 * @param index
	 *            the index
	 */
	public PeakResultDataParameter(int index)
	{
		this.index = index;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultData#getValue(gdsc.smlm.results.PeakResult)
	 */
	public Float getValue(PeakResult result)
	{
		return result.getParameter(index);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultData#getValueName()
	 */
	public String getValueName()
	{
		return PeakResult.getParameterName(index);
	}
}