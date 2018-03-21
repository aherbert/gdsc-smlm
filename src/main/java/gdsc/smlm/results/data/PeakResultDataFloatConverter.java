package gdsc.smlm.results.data;

import gdsc.core.data.utils.Converter;

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
 * Gets a float data value from a result.
 */
public abstract class PeakResultDataFloatConverter extends PeakResultDataFloat
{
	final Converter converter;

	/**
	 * Instantiates a new peak result data float converter.
	 *
	 * @param converter
	 *            the converter
	 */
	public PeakResultDataFloatConverter(Converter converter)
	{
		this.converter = converter;
	}
}