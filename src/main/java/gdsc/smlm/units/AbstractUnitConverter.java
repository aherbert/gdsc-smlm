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
 * Base class for converters
 */
public abstract class AbstractUnitConverter<T extends Unit> implements UnitConverter<T>
{
	private final T from, to;

	/**
	 * Instantiates a new abstract unit converter.
	 *
	 * @param from
	 *            unit to convert from
	 * @param to
	 *            unit to convert to
	 * @throws UnitConversionException
	 *             If the input units are null
	 */
	public AbstractUnitConverter(T from, T to)
	{
		if (from == null)
			throw new UnitConversionException("From unit is null");
		if (to == null)
			throw new UnitConversionException("To unit is null");
		this.from = from;
		this.to = to;
	}

	/**
	 * Instantiates a new abstract unit converter.
	 *
	 * @param from
	 *            unit to convert from
	 * @param to
	 *            unit to convert to
	 * @param suppressExceptions
	 *            the suppress exceptions flag
	 * @throws UnitConversionException
	 *             If the input units are null (and exception are not suppressed)
	 */
	AbstractUnitConverter(T from, T to, boolean suppressExceptions)
	{
		if (from == null && !suppressExceptions)
			throw new UnitConversionException("From unit is null");
		if (to == null && !suppressExceptions)
			throw new UnitConversionException("To unit is null");
		this.from = from;
		this.to = to;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.units.UnitConverter#from()
	 */
	public T from()
	{
		return from;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.units.UnitConverter#to()
	 */
	public T to()
	{
		return to;
	}
}