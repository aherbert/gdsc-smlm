package gdsc.smlm.results.filter;

import gdsc.core.match.FractionalAssignment;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Specifies a classified peak fitting result for use in results scoring
 */
public interface ClassifiedPeakResult extends PreprocessedPeakResult
{
	/**
	 * Get the identifier
	 * 
	 * @return The identifier
	 */
	int getID();

	/**
	 * Get the assignments between this result and the true data.
	 * <p>
	 * The assignments should all have the same predicted id as the value returned from {@link #getID()}
	 * 
	 * @return The assignments
	 */
	FractionalAssignment[] getAssignments();
}
