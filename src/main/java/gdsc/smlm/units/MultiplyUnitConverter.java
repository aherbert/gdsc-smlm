package gdsc.smlm.units;

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
 * Perform conversion by multiplication
 */
public class MultiplyUnitConverter<T extends Unit> extends AbstractUnitConverter<T>
{
	protected final double multiplication;

	/**
	 * Instantiates a new multiplication unit converter.
	 *
	 * @param from
	 *            unit to convert from
	 * @param to
	 *            unit to convert to
	 * @param multiplication
	 *            the multiplication
	 * @throws UnitConversionException
	 *             If the input units are null
	 * @throws UnitConversionException
	 *             If the multiplication is not finite
	 */
	public MultiplyUnitConverter(T from, T to, double multiplication)
	{
		super(from, to);
		if (!Maths.isFinite(multiplication))
			throw new UnitConversionException("multiplication must be finite");
		this.multiplication = multiplication;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.units.UnitConverter#convert(double)
	 */
	public double convert(double value)
	{
		return value * multiplication;
	}
}