package gdsc.smlm.function;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2018 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Factory class for creating a fast log instance
 */
public class FastLogFactory
{
	private static FastLog fastLog = null;
	
	/**
	 * Gets the global fast log instance.
	 *
	 * @return the fast log instance
	 */
	public static FastLog getFastLog()
	{
		if (fastLog == null)
			fastLog = new TurboLog();
		return fastLog;
	}
	
	/**
	 * Gets an instance of FastLog that uses Math.log, i.e. this is not fast.
	 *
	 * @return the log instance
	 */
	public static FastLog getLog()
	{
		return NonFastLog.INSTANCE;
	}
}
