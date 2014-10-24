package gdsc.smlm.engine;

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
 * Define the type of queue used within the fit engine
 */
public enum FitQueue
{
	/**
	 * Block additions if there is a backlog
	 */
	BLOCKING("Blocking"),
	/**
	 * Allow all additions if there is a backlog
	 */
	NON_BLOCKING("Non-blocking"),
	/**
	 * Ignore additions if there is a backlog
	 */
	IGNORE("Ignore");

	private String name;

	private FitQueue(String name)
	{
		this.name = name;
	}

	@Override
	public String toString()
	{
		return name;
	}
}