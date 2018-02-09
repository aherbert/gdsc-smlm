package gdsc.smlm.results.filter;

import gdsc.smlm.data.NamedObject;

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
 * Define the type of localisation precision value.
 */
public enum PrecisionType implements NamedObject
{
	/** No precision is available. */
	NONE("None"),
	/**
	 * The precision is estimated.
	 * <p>
	 * For example this can be done using the Mortensen formula for Gaussian 2D fitting. Note that the
	 * estimate is only correct for fitting to a Gaussian 2D peak and so is an approximation to fitting a spot Point
	 * Spread Function (PSF). If using the Mortensen formula the estimate is expected to use the noise from the fit
	 * region (e.g. the region standard deviation) as background noise. This distinguishes it from
	 * {@link #ESTIMATE_USING_LOCAL_BACKGROUND}.
	 * <p>
	 * See Mortensen, et al (2010) Nature Methods 7, 377-383.
	 * <p>
	 * If not using Gaussian 2D fitting then the Mortensen formula is not applicable and another method must be
	 * employed.
	 */
	ESTIMATE("Estimate"),
	/**
	 * The precision is estimated. This is a special case of {@link #ESTIMATE} for use when using the Mortensen formula
	 * to estimate the precision with the Poisson noise from the fit region (i.e. the fitted background level) as
	 * background noise.
	 * <p>
	 * For example this can be done using the Mortensen formula for Gaussian 2D fitting. Note that the
	 * estimate is only correct for fitting to a Gaussian 2D peak and so is an approximation to fitting a spot Point
	 * Spread Function (PSF). If using the Mortensen formula the estimate is expected to use the Poisson noise from the
	 * fit region (i.e. the fitted background level) as background noise.
	 * <p>
	 * See Mortensen, et al (2010) Nature Methods 7, 377-383.
	 * <p>
	 * If not using Gaussian 2D fitting then the Mortensen formula is not applicable and another method must be
	 * employed.
	 */
	ESTIMATE_USING_LOCAL_BACKGROUND("Estimate using local background", "Estimate b2"),
	/**
	 * The precision is computed using the lower bound on the variance of the estimators of a fitting parameter. This
	 * can be done by computing the inverse of the Fisher information matrix. It assumes that the estimator is 100%
	 * efficient, i.e. achieves the theoretical bounds.
	 * <p>
	 * See <A href="https://en.wikipedia.org/wiki/Cram%C3%A9r%E2%80%93Rao_bound">
	 * https://en.wikipedia.org/wiki/Cram%C3%A9r%E2%80%93Rao_bound</a>
	 */
	CRLB("Cramér-Rao lower bound", "CRLB");

	/** The name. */
	private final String name;

	/** The short name. */
	private final String shortName;

	/**
	 * Instantiates a new precision type.
	 *
	 * @param name
	 *            the name
	 */
	private PrecisionType(String name)
	{
		this(name, name);
	}

	/**
	 * Instantiates a new precision type.
	 *
	 * @param name
	 *            the name
	 * @param sname
	 *            the sname
	 */
	private PrecisionType(String name, String sname)
	{
		this.name = name;
		this.shortName = sname;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Enum#toString()
	 */
	@Override
	public String toString()
	{
		return shortName;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.data.NamedObject#getName()
	 */
	public String getName()
	{
		return name;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.data.NamedObject#getShortName()
	 */
	public String getShortName()
	{
		return shortName;
	}
}