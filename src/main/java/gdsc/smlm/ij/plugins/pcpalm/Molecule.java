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
package gdsc.smlm.ij.plugins.pcpalm;

/**
 * Used to store all the information required for the PC-PALM analysis
 */
public class Molecule
{
	public double x, y, precision, photons;
	
	// Used to construct a single linked list of molecules
	public Molecule next = null;

	public Molecule(double x, double y, double precision, double photons)
	{
		this.x = x;
		this.y = y;
		this.precision = precision;
		this.photons = photons;
	}

	public double distance(Molecule other)
	{
		final double dx = x - other.x;
		final double dy = y - other.y;
		return Math.sqrt(dx * dx + dy * dy);
	}

	public double distance2(Molecule other)
	{
		final double dx = x - other.x;
		final double dy = y - other.y;
		return dx * dx + dy * dy;
	}
}
