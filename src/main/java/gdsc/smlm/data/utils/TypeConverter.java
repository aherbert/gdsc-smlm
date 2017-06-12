package gdsc.smlm.data.utils;

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
 * Define conversion of a type
 */
public interface TypeConverter<T> extends Converter
{
	/**
	 * Specify the source unit to be converted from
	 *
	 * @return the source unit
	 */
	public T from();
	
	/**
	 * Specify the destination unit to be converted to
	 *
	 * @return the destination unit
	 */
	public T to();
}