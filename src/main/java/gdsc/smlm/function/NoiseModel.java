package gdsc.smlm.function;

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
 * Defines the expected variance of a function value (i.e. the noise)
 */
public interface NoiseModel
{
	/**
	 * Calculate the expected variance of a function value, i.e. the noise.
	 * 
	 * @param value
	 *            The model value
	 * @return The expected variance of the value
	 */
	double variance(final double value);
}
