package gdsc.smlm.results;

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
 * Specify the file format when reading results from file
 */
public enum FileFormat
{
	//@formatter:off
	SMLM_TEXT{ public String getName() { return "SMLM Text"; }}, 
	SMLM_BINARY{ public String getName() { return "SMLM Binary"; }},
	RAPID_STORM{ public String getName() { return "RapidSTORM"; }}, 
	NSTORM{ public String getName() { return "NSTORM"; }},
	SMLM_TABLE{ public String getName() { return "SMLM Table"; }}, 
	MALK{ public String getName() { return "MALK"; }}, 
	UNKNOWN{ public String getName() { return "Unknown"; }};
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