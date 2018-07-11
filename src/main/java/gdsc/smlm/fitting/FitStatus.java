/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package gdsc.smlm.fitting;

/**
 * Define the status of a fit result
 */
public enum FitStatus
{
	//@formatter:off
	OK{ @Override
	public String getName() { return "OK"; }},
	SINGULAR_NON_LINEAR_MODEL{ @Override
	public String getName() { return "Singular non-linear model"; }},
	SINGULAR_NON_LINEAR_SOLUTION{ @Override
	public String getName() { return "Singular non-linear solution"; }},
	INVALID_GRADIENTS{ @Override
	public String getName() { return "Invalid gradients"; }},
	FAILED_TO_CONVERGE{ @Override
	public String getName() { return "Failed to converge"; }},
	TOO_MANY_ITERATIONS{ @Override
	public String getName() { return "Too many iterations"; }},
	TOO_MANY_EVALUATIONS{ @Override
	public String getName() { return "Too many evaluations"; }},
	INVALID_LIKELIHOOD{ @Override
	public String getName() { return "Invalid likelihood"; }},
	BAD_PARAMETERS{ @Override
	public String getName() { return "Bad parameters"; }},
	FAILED_TO_ESTIMATE_WIDTH{ @Override
	public String getName() { return "Failed to estimate width"; }},
	COORDINATES_MOVED{ @Override
	public String getName() { return "Coordinates moved"; }},
	OUTSIDE_FIT_REGION{ @Override
	public String getName() { return "Outside fit region"; }},
	INSUFFICIENT_SIGNAL{ @Override
	public String getName() { return "Insufficient signal"; }},
	INSUFFICIENT_SNR{ @Override
	public String getName() { return "Insufficient SNR"; }},
	WIDTH_DIVERGED{ @Override
	public String getName() { return "Width diverged"; }},
	Z_MOVED{ @Override
	public String getName() { return "Z-coordinate moved"; }},
	INSUFFICIENT_PRECISION{ @Override
	public String getName() { return "Insufficient precision"; }},
	NEIGHBOUR_OVERLAP{ @Override
	public String getName() { return "Neighbour overlap"; }},
	FAILED_SMART_FILTER{ @Override
	public String getName() { return "Failed smart filter"; }},
	DRIFT_TO_ANOTHER_RESULT{ @Override
	public String getName() { return "Drift to another result"; }},
	FAILED_VALIDATION{ @Override
	public String getName() { return "Failed validation"; }},
	NO_MODEL_IMPROVEMENT{ @Override
	public String getName() { return "No model improvement"; }},
	LINE_SEARCH_ERROR{ @Override
	public String getName() { return "Line search error"; }},
	UNKNOWN{ @Override
	public String getName() { return "Unknown"; }};
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
