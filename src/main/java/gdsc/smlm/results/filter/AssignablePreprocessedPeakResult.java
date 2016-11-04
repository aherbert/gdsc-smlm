package gdsc.smlm.results.filter;

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
 * Specifies a peak fitting result for use in filtering. Any result implementing this interface can be directly filtered
 * without requiring the filter to be initialised with calibration data. This result can be assigned matches to actual data.
 */
public interface AssignablePreprocessedPeakResult extends PreprocessedPeakResult
{
	/**
	 * Set the assignments between this result and the true data.
	 * 
	 * @param assignments
	 *            The assignments
	 */
	void setAssignments(ResultAssignment[] assignments);
	
	/**
	 * Sets the ignore flag. If true then the result should be ignored from the total counts when scoring.
	 *
	 * @param ignore the new ignore flag
	 */
	void setIgnore(boolean ignore);	
}
