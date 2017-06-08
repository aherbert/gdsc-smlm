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
 * Unit for measuring intensity
 */
public enum IntensityUnit implements Unit
{
	/** Camera count units */
	COUNT
	{
		public String getName()
		{
			return "count";
		}

		UnitConverter<IntensityUnit> buildConverter(IntensityUnit to, double offset, double countPerPhoton)
				throws UnitConversionException
		{
			if (to == IntensityUnit.PHOTON)
				return new AddMultiplyUnitConverter<IntensityUnit>(this, to, -offset, 1.0 / countPerPhoton);
			throw new UnitConversionException(this + " to " + to);
		}

		UnitConverter<IntensityUnit> buildConverter(IntensityUnit to, double countPerPhoton)
				throws UnitConversionException
		{
			if (to == IntensityUnit.PHOTON)
				return new MultiplyUnitConverter<IntensityUnit>(this, to, 1.0 / countPerPhoton);
			throw new UnitConversionException(this + " to " + to);
		}
	},

	/** Photon units */
	PHOTON
	{
		public String getName()
		{
			return "photon";
		}

		UnitConverter<IntensityUnit> buildConverter(IntensityUnit to, double offset, double countPerPhoton)
				throws UnitConversionException
		{
			if (to == IntensityUnit.COUNT)
				return new MultiplyAddUnitConverter<IntensityUnit>(this, to, countPerPhoton, offset);
			throw new UnitConversionException(this + " to " + to);
		}

		UnitConverter<IntensityUnit> buildConverter(IntensityUnit to, double countPerPhoton)
				throws UnitConversionException
		{
			if (to == IntensityUnit.COUNT)
				return new MultiplyUnitConverter<IntensityUnit>(this, to, countPerPhoton);
			throw new UnitConversionException(this + " to " + to);
		}
	},;

	/**
	 * Gets the name.
	 *
	 * @return the name
	 */
	public abstract String getName();

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
	 * @param offset
	 *            the camera offset (in count units)
	 * @param countPerPhoton
	 *            the count per photon
	 * @return the unit converter
	 * @throws UnitConversionException
	 *             if a converter cannot be created
	 */
	public UnitConverter<IntensityUnit> createConverter(IntensityUnit to, double offset, double countPerPhoton)
			throws UnitConversionException
	{
		if (this == to)
			return new IdentityUnitConverter<IntensityUnit>(this);

		if (!Maths.isFinite(offset))
			throw new UnitConversionException("offset must be finite");
		if (!(countPerPhoton > 0 && countPerPhoton <= java.lang.Double.MAX_VALUE))
			throw new UnitConversionException("count/photon must be positive");

		return buildConverter(to, offset, countPerPhoton);
	}

	/**
	 * Build the converter for the unit.
	 *
	 * @param to
	 *            the to
	 * @param offset
	 *            the camera offset (in count units)
	 * @param countPerPhoton
	 *            the count per photon
	 * @return the unit converter
	 * @throws UnitConversionException
	 *             if a converter cannot be created
	 */
	abstract UnitConverter<IntensityUnit> buildConverter(IntensityUnit to, double offset, double countPerPhoton)
			throws UnitConversionException;

	/**
	 * Creates the converter.
	 *
	 * @param to
	 *            the to
	 * @param countPerPhoton
	 *            the count per photon
	 * @return the unit converter
	 * @throws UnitConversionException
	 *             if a converter cannot be created
	 */
	public UnitConverter<IntensityUnit> createConverter(IntensityUnit to, double countPerPhoton)
			throws UnitConversionException
	{
		if (this == to)
			return new IdentityUnitConverter<IntensityUnit>(this);

		if (!(countPerPhoton > 0 && countPerPhoton <= java.lang.Double.MAX_VALUE))
			throw new UnitConversionException("count/photon must be positive");

		return buildConverter(to, countPerPhoton);
	}

	/**
	 * Build the converter for the unit.
	 *
	 * @param to
	 *            the to
	 * @param countPerPhoton
	 *            the count per photon
	 * @return the unit converter
	 * @throws UnitConversionException
	 *             if a converter cannot be created
	 */
	abstract UnitConverter<IntensityUnit> buildConverter(IntensityUnit to, double countPerPhoton)
			throws UnitConversionException;
}