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
 * Exception to throw if conversion is not possible
 */
public class UnitConversionException extends RuntimeException
{
	private static final long serialVersionUID = 2470815639465684383L;

	public UnitConversionException()
	{
		super();
	}

	public UnitConversionException(String message)
	{
		super(message);
	}

	public UnitConversionException(String message, Throwable cause)
	{
		super(message, cause);
	}

	public UnitConversionException(Throwable cause)
	{
		super(cause);
	}
}
