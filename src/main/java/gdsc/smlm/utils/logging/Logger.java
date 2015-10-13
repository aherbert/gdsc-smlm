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
 * Defines a simple interface for logging information
 */
public interface Logger
{
	/**
	 * Log the message
	 * @param message
	 */
	void info(String message);
	
	/**
	 * Log the arguments using the given format
	 * @param format
	 * @param args
	 */
	void info(String format, Object... args);
	
	/**
	 * Log the message
	 * @param message
	 */
	void debug(String message);
	
	/**
	 * Log the arguments using the given format
	 * @param format
	 * @param args
	 */
	void debug(String format, Object... args);
	
	/**
	 * Log the message
	 * @param message
	 */
	void error(String message);
	
	/**
	 * Log the arguments using the given format
	 * @param format
	 * @param args
	 */
	void error(String format, Object... args);
}
