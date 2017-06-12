package gdsc.smlm.data.units;

import java.math.BigDecimal;
import java.math.MathContext;

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
 * Class for rounding
 */
public class MathContextRounder implements Rounder
{
	final MathContext mathContext;

	/**
	 * Instantiates a new math context rounder.
	 *
	 * @param mathContext
	 *            the math context
	 * @throws IllegalArgumentException
	 *             if the mathContext is null
	 */
	public MathContextRounder(MathContext mathContext)
	{
		if (mathContext == null)
			throw new IllegalArgumentException("MathContext must not be null");
		this.mathContext = mathContext;
	}

	/**
	 * Instantiates a new math context rounder.
	 *
	 * @param precision
	 *            The non-negative {@code int} precision setting.
	 * @throws IllegalArgumentException
	 *             if the {@code precision} parameter is less
	 *             than zero.
	 */
	public MathContextRounder(int precision)
	{
		mathContext = new MathContext(precision);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.data.Rounder#round(double)
	 */
	public double round(double value)
	{
		if (Math.abs(value) <= Double.MAX_VALUE)
			return new BigDecimal(value).round(mathContext).doubleValue();
		return value; // NaN or infinite
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.data.Rounder#toString(double)
	 */
	public String toString(double value)
	{
		if (Math.abs(value) <= Double.MAX_VALUE)
			return new BigDecimal(value).round(mathContext).toString();
		if (value == Double.POSITIVE_INFINITY)
			return "Infinity";
		if (value == Double.NEGATIVE_INFINITY)
			return "-Infinity";
		return "NaN";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.data.Rounder#round(float)
	 */
	public float round(float value)
	{
		if (Math.abs(value) <= Float.MAX_VALUE)
			return new BigDecimal(value).round(mathContext).floatValue();
		return value; // NaN or infinite
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.data.Rounder#toString(float)
	 */
	public String toString(float value)
	{
		if (Math.abs(value) <= Float.MAX_VALUE)
			return new BigDecimal(value).round(mathContext).toString();
		if (value == Float.POSITIVE_INFINITY)
			return "Infinity";
		if (value == Float.NEGATIVE_INFINITY)
			return "-Infinity";
		return "NaN";
	}
}
