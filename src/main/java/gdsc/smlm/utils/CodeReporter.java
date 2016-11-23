package gdsc.smlm.utils;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Allow the position of executing code to be reported
 */
public class CodeReporter
{
	/**
	 * Print the classname, methodname and line number from the throwable
	 *
	 * @param t
	 *            the throwable
	 */
	public static void debug(Throwable t)
	{
		StackTraceElement e = t.getStackTrace()[0];
		System.err.printf("%s:%s:%d\n", e.getClassName(), e.getMethodName(), e.getLineNumber());
	}
}
