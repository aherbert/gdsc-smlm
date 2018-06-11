/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 * 
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package gdsc.smlm.model;


/**
 * Define the diffusion type
 */
public enum DiffusionType
{
	//@formatter:off
	/**
	 * A random walk
	 */
	RANDOM_WALK{ public String getName() { return "Random Walk"; }},
	/**
	 * A grid walk using defined step sizes in each dimension
	 */
	GRID_WALK{ public String getName() { return "Grid Walk"; }},
	/**
	 * A random walk along a linear axis
	 */
	LINEAR_WALK{ public String getName() { return "Linear Walk"; }};
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
				if (text.equalsIgnoreCase(type.getName()))
				{
					return type;
				}
			}
		}
		return null;
	}
}
