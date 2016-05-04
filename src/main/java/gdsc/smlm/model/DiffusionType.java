package gdsc.smlm.model;

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
 * Define the diffusion type
 */
public enum DiffusionType
{
	/**
	 * A random walk
	 */
	RANDOM_WALK("Random Walk"),
	/**
	 * A grid walk using defined step sizes in each dimension
	 */
	GRID_WALK("Grid Walk"),
	/**
	 * A random walk along a linear axis
	 */
	LINEAR_WALK("Linear Walk");

	private String name;

	private DiffusionType(String name)
	{
		this.name = name;
	}

	@Override
	public String toString()
	{
		return name;
	}

	/**
	 * Get the diffusion type from a given string. Returns null if the text is not a valid type.
	 * 
	 * @param text The text
	 * @return The diffusion type (or null)
	 */
	public static DiffusionType fromString(String text)
	{
		if (text != null)
		{
			text = text.trim();
			for (DiffusionType type : DiffusionType.values())
			{
				if (text.equalsIgnoreCase(type.name))
				{
					return type;
				}
			}
		}
		return null;
	}
}