package gdsc.smlm.filters;

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
 * Contains common functionality for filters.
 */
public abstract class BaseFilter implements Cloneable
{
	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public BaseFilter clone()
	{
		try
		{
			return (BaseFilter) super.clone();
		}
		catch (CloneNotSupportedException e)
		{
			// Ignore
		}
		return null;
	}
}