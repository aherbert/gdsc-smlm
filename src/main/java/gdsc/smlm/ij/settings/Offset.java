package gdsc.smlm.ij.settings;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2015 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Contain the settings for the centre offset of a PSF
 */
public class Offset
{
	private final double cx;
	private final double cy;

	public Offset(double cx, double cy)
	{
		this.cx = cx;
		this.cy = cy;
	}

	public double getCx()
	{
		return cx;
	}

	public double getCy()
	{
		return cy;
	}
}
