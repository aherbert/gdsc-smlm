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
public interface CoordinateStore
{
	/**
	 * Gets the resolution of the store.
	 *
	 * @return the resolution
	 */
	public double getResolution();

	/**
	 * Queue a coordinate to the store.
	 * <p>
	 * It is not added to the store until flush is called.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 */
	public void addToQueue(double x, double y);

	/**
	 * Flush the queue to the store
	 */
	public void flush();

	/**
	 * Add a coordinate to the store. Assumes that the coordinates are within the size of the grid otherwise they will be ignored.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 */
	public void add(double x, double y);

	/**
	 * Clear to the store.
	 */
	public void clear();

	/**
	 * Check if the store contains the coordinates within the configured resolution.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @return true, if the store contains another coordinate closer than the resolution
	 */
	public boolean contains(double x, double y);

	/**
	 * Find the closest coordinate within the configured resolution.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @return the coordinate closer than the resolution (or null)
	 */
	public double[] find(double x, double y);
	
	/**
	 * Create a new instance.
	 *
	 * @return the new coordinate store
	 */
	public CoordinateStore newInstance();
	
	/**
	 * Resize to the given dimensions. If these match the existing dimensions the current store is returned. Otherwise a
	 * new store is returned.
	 *
	 * @param maxx
	 *            the max x coordinate value
	 * @param maxy
	 *            the max y coordinate value
	 * @return the coordinate store
	 */
	public CoordinateStore resize(int maxx, int maxy);
}