package gdsc.smlm.fitting;

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
 * Define the status of a fit result
 */
public enum FitStatus
{
	//@formatter:off
	OK{ public String getName() { return "OK"; }},
	SINGULAR_NON_LINEAR_MODEL{ public String getName() { return "Singular non-linear model"; }},
	SINGULAR_NON_LINEAR_SOLUTION{ public String getName() { return "Singular non-linear solution"; }},
	INVALID_GRADIENTS{ public String getName() { return "Invalid gradients"; }},
	FAILED_TO_CONVERGE{ public String getName() { return "Failed to converge"; }},
	TOO_MANY_ITERATIONS{ public String getName() { return "Too many iterations"; }},
	TOO_MANY_EVALUATIONS{ public String getName() { return "Too many evaluations"; }},
	INVALID_LIKELIHOOD{ public String getName() { return "Invalid likelihood"; }},
	BAD_PARAMETERS{ public String getName() { return "Bad parameters"; }},
	FAILED_TO_ESTIMATE_WIDTH{ public String getName() { return "Failed to estimate width"; }},
	COORDINATES_MOVED{ public String getName() { return "Coordinates moved"; }},
	OUTSIDE_FIT_REGION{ public String getName() { return "Outside fit region"; }},
	INSUFFICIENT_SIGNAL{ public String getName() { return "Insufficient signal"; }},
	WIDTH_DIVERGED{ public String getName() { return "Width diverged"; }},
	INSUFFICIENT_PRECISION{ public String getName() { return "Insufficient precision"; }},
	NEIGHBOUR_OVERLAP{ public String getName() { return "Neighbour overlap"; }},
	FAILED_SMART_FILTER{ public String getName() { return "Failed smart filter"; }},
	DRIFT_TO_ANOTHER_RESULT{ public String getName() { return "Drift to another result"; }},
	FAILED_VALIDATION{ public String getName() { return "Failed validation"; }},
	NO_MODEL_IMPROVEMENT{ public String getName() { return "No model improvement"; }},
	UNKNOWN{ public String getName() { return "Unknown"; }};
	//@formatter:on

	@Override
	public String toString()
	{
		return getName();
	}

	/**
	 * Gets the name.
	 *
	 * @return the name
	 */
	abstract public String getName();
}