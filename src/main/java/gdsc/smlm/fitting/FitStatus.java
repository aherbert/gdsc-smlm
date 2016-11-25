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
	/**
	 * 
	 */
	OK("OK"),
	/**
	 * 
	 */
	SINGULAR_NON_LINEAR_MODEL("Singular non-linear model"),
	/**
	 * 
	 */
	SINGULAR_NON_LINEAR_SOLUTION("Singular non-linear solution"),
	/**
	 * 
	 */
	INVALID_GRADIENTS_IN_NON_LINEAR_MODEL("Invalid gradients in non-linear model"),
	/**
	 * 
	 */
	FAILED_TO_CONVERGE("Failed to converge"),
	/**
	 * 
	 */
	TOO_MANY_ITERATIONS("Too many iterations"),
	/**
	 * 
	 */
	TOO_MANY_EVALUATIONS("Too many evaluations"),
	/**
	 * 
	 */
	INVALID_LIKELIHOOD("Invalid likelihood"),
	/**
	 * 
	 */
	BAD_PARAMETERS("Bad parameters"),
	/**
	 * 
	 */
	FAILED_TO_ESTIMATE_WIDTH("Failed to estimate width"),
	/**
	 * 
	 */
	COORDINATES_MOVED("Coordinates moved"),
	/**
	 * 
	 */
	OUTSIDE_FIT_REGION("Outside fit region"),
	/**
	 * 
	 */
	INSUFFICIENT_SIGNAL("Insufficient signal"),
	/**
	 * 
	 */
	WIDTH_DIVERGED("Width diverged"),
	/**
	 * 
	 */
	INSUFFICIENT_PRECISION("Insufficient precision"),
	/**
	 * 
	 */
	NEIGHBOUR_OVERLAP("Neighbour overlap"),
	/**
	 * 
	 */
	FAILED_SMART_FILTER("Failed smart filter"),
	/**
	 * 
	 */
	DRIFT_TO_ANOTHER_RESULT("Drift to another result"),
	/**
	 * 
	 */
	FAILED_VALIDATION("Failed validation"),
	/**
	 * 
	 */
	NO_MODEL_IMPROVEMENT("No model improvement"),
	/**
	 * 
	 */
	UNKNOWN("Unknown");

	private String name;

	private FitStatus(String name)
	{
		this.name = name;
	}

	@Override
	public String toString()
	{
		return name;
	}
}