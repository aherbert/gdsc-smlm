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
public class Compound
{
	@XStreamAsAttribute
	public double fraction;
	@XStreamAsAttribute
	public double D;
	@XStreamAsAttribute
	public String diffusionType;
	public Atom[] atoms;

	public Compound(double fraction, double D, String diffusionType, Atom... atoms)
	{
		this.fraction = fraction;
		this.D = D;
		this.diffusionType = diffusionType;
		this.atoms = atoms;
	}
}