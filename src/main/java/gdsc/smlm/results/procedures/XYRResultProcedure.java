package gdsc.smlm.results.procedures;

import gdsc.smlm.results.PeakResult;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Interface for accessing the results
 */
public interface XYRResultProcedure
{
	/**
	 * Executes this procedure.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param result
	 *            the result
	 */
	void executeXYR(float x, float y, PeakResult result);
}
