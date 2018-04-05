package gdsc.smlm.results.procedures;

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
 * Interface for accessing the results
 */
public interface BIRResultProcedure
{
	/**
	 * Executes this procedure.
	 *
	 * @param background
	 *            the background
	 * @param intensity
	 *            the intensity
	 * @param result
	 *            the result
	 */
	void executeBIR(float background, float intensity, PeakResult result);
}
