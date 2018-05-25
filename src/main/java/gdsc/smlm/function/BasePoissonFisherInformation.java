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
 * Base class for any probability distribution that can calculate the Fisher information for a Poisson distributed mean.
 * <p>
 * <a href="https://en.wikipedia.org/wiki/Poisson_distribution">https://en.wikipedia.org/wiki/Poisson_distribution</a>
 */
public abstract class BasePoissonFisherInformation implements FisherInformation, Cloneable
{
	/**
	 * The lowest value for the mean that can be computed. This is the lowest value where the reciprocal is not infinity
	 */
	public static final double MIN_MEAN = Double.longBitsToDouble(0x4000000000001L);

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.FisherInformation#isValid(double)
	 */
	public boolean isValid(double t)
	{
		return t >= MIN_MEAN;
	}

	/**
	 * Gets the alpha scale of the Poisson Fisher information.
	 * This is the Fisher information relative to the Fisher information of a pure Poisson distribution:
	 * 
	 * <pre>
	 * alpha = FI / (Poisson FI) = FI * t
	 * </pre>
	 *
	 * @param t
	 *            parameter θ of a distribution that models X
	 * @return the alpha scale
	 * @throws IllegalArgumentException
	 *             if the parameter is not in the valid range
	 */
	public abstract double getAlpha(double t);

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	public BasePoissonFisherInformation clone()
	{
		try
		{
			BasePoissonFisherInformation fi = (BasePoissonFisherInformation) super.clone();
			fi.postClone();
			return fi;
		}
		catch (CloneNotSupportedException e)
		{
			// Should not happen
			return null;
		}
	}

	/**
	 * Run any actions after clone, for example creating new instance fields if they cannot be shared.
	 */
	protected abstract void postClone();
}