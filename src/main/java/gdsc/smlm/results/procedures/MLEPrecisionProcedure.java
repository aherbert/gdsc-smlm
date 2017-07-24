package gdsc.smlm.results.procedures;

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
 * Interface for accessing the localisation precision of Gaussian 2D fitting computed using the Mortensen formula for
 * Maximum Likelihood Estimation using a local noise estimate.
 * <p>
 * See Mortensen, et al (2010) Nature Methods 7, 377-383, equation 6.
 */
public interface MLEPrecisionProcedure
{
	/**
	 * Executes this procedure.
	 *
	 * @param precision
	 *            the precision
	 */
	void executeMLEPrecision(double precision);
}
