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
 * Stores a set of results within a grid arrangement at a given resolution. Allows checking for duplicates. The XY and Z
 * resolution can be different.
 */
public interface CoordinateStore
{
	/**
	 * Gets the XY resolution of the store. If negative then nothing is stored.
	 *
	 * @return the XY resolution
	 */
	public double getXYResolution();

	/**
	 * Gets the Z resolution of the store. If negative then this is ignored and the store behaves as if processing 2D
	 * coordinates.
	 *
	 * @return the Z resolution
	 */
	public double getZResolution();

	/**
	 * Queue a coordinate to the store.
	 * <p>
	 * It is not added to the store until flush is called.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param z
	 *            the z
	 */
	public void addToQueue(double x, double y, double z);

	/**
	 * Flush the queue to the store
	 */
	public void flush();

	/**
	 * Add a coordinate to the store. Assumes that the coordinates are within the size of the grid otherwise they will
	 * be ignored.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param z
	 *            the z
	 */
	public void add(double x, double y, double z);

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
	 * @param z
	 *            the z
	 * @return true, if the store contains another coordinate closer than the resolution
	 */
	public boolean contains(double x, double y, double z);

	/**
	 * Find the closest coordinate within the configured resolution.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param z
	 *            the z
	 * @return the coordinate closer than the resolution (or null)
	 */
	public double[] find(double x, double y, double z);

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
	 * @param minx
	 *            the min x coordinate value
	 * @param miny
	 *            the min y coordinate value
	 * @param width
	 *            the width
	 * @param height
	 *            the height
	 * @return the coordinate store
	 */
	public CoordinateStore resize(int minx, int miny, int width, int height);
}