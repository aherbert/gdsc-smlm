package gdsc.smlm.data.units;

import gdsc.core.utils.Maths;

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
 * Perform conversion by addition then multiplication
 */
public class AddMultiplyUnitConverter<T extends Unit> extends MultiplyUnitConverter<T>
{
	private final double addition;

	/**
	 * Instantiates a new add then multiplication unit converter.
	 *
	 * @param from
	 *            unit to convert from
	 * @param to
	 *            unit to convert to
	 * @param addition
	 *            the value to add before multiplication
	 * @param multiplication
	 *            the multiplication
	 * @throws ConversionException
	 *             If the input units are null
	 * @throws ConversionException
	 *             If the multiplication is not finite
	 * @throws ConversionException
	 *             If the addition is not finite
	 */
	public AddMultiplyUnitConverter(T from, T to, double addition, double multiplication)
	{
		super(from, to, multiplication);
		if (!Maths.isFinite(addition))
			throw new ConversionException("addition must be finite");
		this.addition = addition;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.units.UnitConverter#convert(double)
	 */
	public double convert(double value)
	{
		return (value + addition) * multiplication;
	}
}