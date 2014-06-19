package gdsc.smlm.ij.settings;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;

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
 * Specify a set of coordinates with mass. 
 * <p>
 * Used for XStream serialisation of the Create Data compounds.
 */
public class Atom
{
	@XStreamAsAttribute
	public double mass;
	@XStreamAsAttribute
	public double x;
	@XStreamAsAttribute
	public double y;
	@XStreamAsAttribute
	public double z;

	public Atom(double mass, double x, double y, double z)
	{
		this.mass = mass;
		this.x = x;
		this.y = y;
		this.z = z;
	}
}