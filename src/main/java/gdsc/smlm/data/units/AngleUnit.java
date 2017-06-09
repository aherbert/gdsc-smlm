/*
 * 
 */
package gdsc.smlm.data.units;

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
 * Unit for measuring angle
 */
public enum AngleUnit implements Unit
{
	/** Radian units */
	RADIAN
	{
		public String getName()
		{
			return "radian";
		}

		public String getShortName()
		{
			return "rad";
		}

		UnitConverter<AngleUnit> buildConverter(AngleUnit to) throws UnitConversionException
		{
			if (to == AngleUnit.DEGREE)
				return new MultiplyUnitConverter<AngleUnit>(this, to, 180.0 / Math.PI);
			throw new UnitConversionException(this + " to " + to);
		}
	},

	/** Degree units */
	DEGREE
	{
		public String getName()
		{
			return "degree";
		}

		public String getShortName()
		{
			return "°";
		}

		UnitConverter<AngleUnit> buildConverter(AngleUnit to) throws UnitConversionException
		{
			if (to == AngleUnit.RADIAN)
				return new MultiplyUnitConverter<AngleUnit>(this, to, Math.PI / 180.0);
			throw new UnitConversionException(this + " to " + to);
		}
	},;

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Enum#toString()
	 */
	public String toString()
	{
		return getName();
	}

	/**
	 * Creates the converter.
	 *
	 * @param to
	 *            the to
	 * @param nmPerPixel
	 *            the nm per pixel
	 * @return the unit converter
	 * @throws UnitConversionException
	 *             if a converter cannot be created
	 */
	public UnitConverter<AngleUnit> createConverter(AngleUnit to) throws UnitConversionException
	{
		if (this == to)
			return new IdentityUnitConverter<AngleUnit>(this);

		return buildConverter(to);
	}

	/**
	 * Build the converter for the unit.
	 *
	 * @param to
	 *            the to
	 * @return the unit converter
	 * @throws UnitConversionException
	 *             if a converter cannot be created
	 */
	abstract UnitConverter<AngleUnit> buildConverter(AngleUnit to) throws UnitConversionException;
}