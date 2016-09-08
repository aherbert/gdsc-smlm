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
	 * Get the assignments between this result and the true data.
	 * <p>
	 * The assignments should all have the same predicted Id. The actual Id should be a value starting from 0 and
	 * incrementing for each actual result in the frame that is scored.
	 * 
	 * @param predictedId
	 *            The predicted Id
	 * @return The assignments
	 */
	FractionalAssignment[] getAssignments(int predictedId);
}
