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
 * Unit for measuring distance
 */
public enum DistanceUnit implements Unit
{
	/** Pixel units */
	PIXEL
	{
		public String getName()
		{
			return "pixel";
		}

		public String getShortName()
		{
			return "px";
		}
		
		UnitConverter<DistanceUnit> buildConverter(DistanceUnit to, double nmPerPixel) throws UnitConversionException
		{
			if (to == DistanceUnit.NM)
				return new MultiplyUnitConverter<DistanceUnit>(this, to, nmPerPixel);
			if (to == DistanceUnit.UM)
				return new MultiplyUnitConverter<DistanceUnit>(this, to, nmPerPixel / 1e3);
			throw new UnitConversionException(this + " to " + to);
		}
	},

	/** Micrometer units */
	UM
	{
		public String getName()
		{
			return "micrometer";
		}

		public String getShortName()
		{
			return "um";
		}

		UnitConverter<DistanceUnit> buildConverter(DistanceUnit to, double nmPerPixel) throws UnitConversionException
		{
			if (to == DistanceUnit.PIXEL)
				return new MultiplyUnitConverter<DistanceUnit>(this, to, 1e3 / nmPerPixel);
			if (to == DistanceUnit.NM)
				return new MultiplyUnitConverter<DistanceUnit>(this, to, 1e3);
			throw new UnitConversionException(this + " to " + to);
		}
	},

	/** Nanometer units */
	NM
	{
		public String getName()
		{
			return "nanometer";
		}

		public String getShortName()
		{
			return "nm";
		}
		
		UnitConverter<DistanceUnit> buildConverter(DistanceUnit to, double nmPerPixel) throws UnitConversionException
		{
			if (to == DistanceUnit.PIXEL)
				return new MultiplyUnitConverter<DistanceUnit>(this, to, 1.0 / nmPerPixel);
			if (to == DistanceUnit.UM)
				return new MultiplyUnitConverter<DistanceUnit>(this, to, 1e-3);
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
	public UnitConverter<DistanceUnit> createConverter(DistanceUnit to, double nmPerPixel)
			throws UnitConversionException
	{
		if (this == to)
			return new IdentityUnitConverter<DistanceUnit>(this);

		if (!(nmPerPixel > 0 && nmPerPixel <= java.lang.Double.MAX_VALUE))
			throw new UnitConversionException("nm/pixel must be positive");

		return buildConverter(to, nmPerPixel);
	}

	/**
	 * Build the converter for the unit.
	 *
	 * @param to
	 *            the to
	 * @param nmPerPixel
	 *            the nm per pixel
	 * @return the unit converter
	 * @throws UnitConversionException
	 *             if a converter cannot be created
	 */
	abstract UnitConverter<DistanceUnit> buildConverter(DistanceUnit to, double nmPerPixel)
			throws UnitConversionException;

	/**
	 * Creates the converter.
	 *
	 * @param to
	 *            the to
	 * @return the unit converter
	 * @throws UnitConversionException
	 *             if a converter cannot be created
	 */
	public UnitConverter<DistanceUnit> createConverter(DistanceUnit to) throws UnitConversionException
	{
		if (this == to)
			return new IdentityUnitConverter<DistanceUnit>(this);
		if (to == DistanceUnit.PIXEL)
			throw new UnitConversionException(this + " to " + to + " requires nm/pixel");
		return buildConverter(to, 1.0);
	}
}