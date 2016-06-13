package gdsc.smlm.function;

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

import org.apache.commons.math3.distribution.PoissonDistribution;

/**
 * Implements the probability density function for a Poisson distribution.
 * <p>
 * This is a simple implementation of the LikelihoodFunction interface.
 */
public class PoissonFunction implements LikelihoodFunction
{
	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.LikelihoodFunction#likelihood(double, double)
	 */
	public double likelihood(double o, double e)
	{

		return new PoissonDistribution(e).probability((int) Math.round(o));
	}
}