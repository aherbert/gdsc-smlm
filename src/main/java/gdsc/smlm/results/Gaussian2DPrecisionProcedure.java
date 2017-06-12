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
 * Interface for accessing the results of Gaussian 2D fitting required for computing the precision
 */
public interface Gaussian2DPrecisionProcedure
{
	/**
	 * Executes this procedure.
	 *
	 * @param background
	 *            the background (in photons)
	 * @param intensity
	 *            the intensity (in photons)
	 * @param s
	 *            the Gaussian standard deviation (in nm)
	 */
	void execute(float background, float intensity, float s);
}
