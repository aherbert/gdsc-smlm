package gdsc.smlm.function;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2018 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Calculate half the Fisher information for a Poisson distribution. 
 * This is suitable for modelling an EM-CCD camera (Poisson-Gamma-Gaussian distribution) when the mean is high.
 * <p>
 * <a href="https://en.wikipedia.org/wiki/Poisson_distribution">https://en.wikipedia.org/wiki/Poisson_distribution</a>
 */
public class HalfPoissonFisherInformation extends BasePoissonFisherInformation
{
	/*
	 * {@inheritDoc}
	 * <p>
	 * The input parameter refers to the mean of the Poisson distribution. The Fisher information is 1/(2*mean).
	 * 
	 * @see gdsc.smlm.function.FisherInformation#getFisherInformation(double)
	 */
	public double getFisherInformation(double t)
	{
		if (t <= 0)
			throw new IllegalArgumentException("Poisson mean must be positive");
		return 0.5 / t;
	}
	
	@Override
	public double getAlpha(double t)
	{
		return 0.5;
	}
	
	@Override
	protected void postClone()
	{
		// Nothing to do.
	}
}