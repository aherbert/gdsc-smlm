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
	SMLM_TEXT("SMLM Text"), SMLM_BINARY("SMLM Binary"), RAPID_STORM("RapidSTORM"), NSTORM("NSTORM"), SMLM_TABLE(
			"SMLM Table"), UNKNOWN("Unknown");

	private String name;

	private FileFormat(String name)
	{
		this.name = name;
	}

	@Override
	public String toString()
	{
		return name;
	}
}