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
 * Logs messages to the Java console 
 */
public class ConsoleLogger implements Logger
{
	/**
	 * Log the message to the Java console. Appends a new line to the message.
	 * @param message
	 */
	public void info(String message)
	{
		System.out.println(message);
	}

	/**
	 * Log the arguments using the given format to the Java console. Appends a new line to the message.
	 * @param format
	 * @param args
	 */
	public void info(String format, Object... args)
	{
		System.out.printf(format, args);
		System.out.printf("\n");
	}

	@Override
	public void debug(String message)
	{
		info(message);
	}

	@Override
	public void debug(String format, Object... args)
	{
		info(format, args);
	}

	@Override
	public void error(String message)
	{
		info(message);
	}

	@Override
	public void error(String format, Object... args)
	{
		info(format, args);
	}
}
