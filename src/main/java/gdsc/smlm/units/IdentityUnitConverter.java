package gdsc.smlm.units;

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
 * Perform no conversion
 */
public class IdentityUnitConverter<T extends Unit> extends AbstractUnitConverter<T>
{
	/**
	 * Instantiates a new identity unit converter.
	 *
	 * @param units
	 *            the units (can be null)
	 */
	public IdentityUnitConverter(T units)
	{
		super(units, units, true);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.units.UnitConverter#convert(double)
	 */
	public double convert(double value)
	{
		return value;
	}
}