package gdsc.smlm.utils.logging;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Logs messages to nowhere
 */
public class NullLogger implements Logger
{
	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.utils.logging.Logger#info(java.lang.String)
	 */
	public void info(String message)
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.utils.logging.Logger#info(java.lang.String, java.lang.Object[])
	 */
	public void info(String format, Object... args)
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.utils.logging.Logger#debug(java.lang.String)
	 */
	public void debug(String message)
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.utils.logging.Logger#debug(java.lang.String, java.lang.Object[])
	 */
	public void debug(String format, Object... args)
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.utils.logging.Logger#error(java.lang.String)
	 */
	public void error(String message)
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.utils.logging.Logger#error(java.lang.String, java.lang.Object[])
	 */
	public void error(String format, Object... args)
	{
	}
}
