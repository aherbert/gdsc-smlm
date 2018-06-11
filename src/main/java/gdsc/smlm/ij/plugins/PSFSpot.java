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
package gdsc.smlm.ij.plugins;

import gdsc.core.match.BasePoint;
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

	/** The amplitude. */
	public double amplitude = 0;

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
