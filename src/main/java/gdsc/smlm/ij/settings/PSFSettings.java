package gdsc.smlm.ij.settings;

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
 * Contain the settings for the PSF Creator plugin
 */
public class PSFSettings
{
	public int zCentre = 1;
	public double nmPerPixel = 100;
	public double nmPerSlice = 20;
	public int nImages = 1;

	public PSFSettings(int zCentre, double nmPerPixel, double nmPerSlice, int nImages)
	{
		this.zCentre = zCentre;
		this.nmPerPixel = nmPerPixel;
		this.nmPerSlice = nmPerSlice;
		this.nImages = nImages;
	}
}
