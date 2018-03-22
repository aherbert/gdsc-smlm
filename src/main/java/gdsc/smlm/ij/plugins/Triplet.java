package gdsc.smlm.ij.plugins;

/*----------------------------------------------------------------------------- 
 * GDSC Plugins for ImageJ
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
 * A generic triplet.
 * 
 * @author Alex Herbert
 */
public class Triplet<A, B, C>
{
	public final A a;
	public final B b;
	public final C c;

	public Triplet(A a, B b, C c)
	{
		this.a = a;
		this.b = b;
		this.c = c;
	}
}
