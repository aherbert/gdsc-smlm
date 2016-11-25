package gdsc.smlm.results.filter;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Stores a set of results within a grid arrangement at a given resolution. Allows checking for duplicates.
 */
public class CoordinateStoreFactory
{
	/**
	 * Creates the coordinate store
	 *
	 * @param maxx
	 *            the max x coordinate value
	 * @param maxy
	 *            the max y coordinate value
	 * @param resolution
	 *            the resolution
	 * @return the coordinate store
	 */
	public static CoordinateStore create(int maxx, int maxy, double resolution)
	{
		if (resolution <= 0)
			return new NullCoordinateStore();
		
		// This should be faster (for additions and block lookup) as it has a fixed block resolution of 1.
		// However it may be slower if the distance is much lower than 1 and there are many points close 
		// to the resolution distance as it will have to compute the distance for each. As a compromise
		// we only use it when the resolution is above the min block size of the default store.
		if (resolution >= GridCoordinateStore.MINIMUM_BLOCK_SIZE && resolution <= 1)
			return new GridCoordinateStore1(maxx, maxy, resolution);
		
		return new GridCoordinateStore(maxx, maxy, resolution);
	}
}