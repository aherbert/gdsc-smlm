package gdsc.smlm.results;

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
 * Provides access methods for a read-only view of peak results.
 */
public interface PeakResultView
{
	/**
	 * Gets the results by frame.
	 *
	 * @param frame
	 *            the frame
	 * @return the results
	 */
	public PeakResult[] getResultsByFrame(int frame);

	/**
	 * Gets the results by id.
	 *
	 * @param id
	 *            the id
	 * @return the results
	 */
	public PeakResult[] getResultsById(int id);
}
