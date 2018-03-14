package gdsc.smlm.ij.ij3d;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
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

import ij3d.Content;

/**
 * Extend the Content class to use CustomContentInstant
 */
public class CustomContent extends Content
{
	public CustomContent(String name)
	{
		super(name);
		// Replace the default from the super constructor
		final CustomContentInstant ci = new CustomContentInstant(name + "_#0");
		addInstant(ci);
	}
}
