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
	//@formatter:off
	/**
	 * Block additions if there is a backlog
	 */
	BLOCKING{ public String getName() { return "Blocking"; }},
	/**
	 * Allow all additions if there is a backlog
	 */
	NON_BLOCKING{ public String getName() { return "Non-blocking"; }},
	/**
	 * Ignore additions if there is a backlog
	 */
	IGNORE{ public String getName() { return "Ignore"; }};
	//@formatter:on

	@Override
	public String toString()
	{
		return getName();
	}

	/**
	 * Gets the name.
	 *
	 * @return the name
	 */
	abstract public String getName();
}