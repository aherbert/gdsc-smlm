package gdsc.smlm.ij.plugins;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2014 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import uk.ac.sussex.gdsc.core.match.BasePoint;
import gdsc.smlm.results.PeakResult;

/**
 * Stores details of a simulated localisation. Contains details of the amount of signal that occurs due to overlap
 * with neighbour PSFs.
 */
public class PSFSpot extends BasePoint
{
	public final int t;
	public final PeakResult peakResult;
	/**
	 * The amount of total intensity contributed within the region of this spot from overlapping PSFs, i.e. how much
	 * more signal is in the area of this spot due to other PSFs.
	 */
	public float intensityOffset = 0;
	/**
	 * The amount of total background contributed within the region of this spot from overlapping PSFs, i.e. how much
	 * higher is this spot due to other PSFs.
	 */
	public float backgroundOffset = 0;

	public PSFSpot(int t, float x, float y, PeakResult peakResult)
	{
		super(x, y);
		this.t = t;
		this.peakResult = peakResult;
	}

	public int getTime()
	{
		return t;
	}
}