package gdsc.smlm.results.match;

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
 * Stores a 2D/3D point
 */
public interface Coordinate
{
	/**
	 * @return The X-coordinate
	 */
	public float getX();
	/**
	 * @return The Y-coordinate
	 */
	public float getY();
	/**
	 * @return The Z-coordinate
	 */
	public float getZ();

	/**
	 * Calculate the XY distance to the given coordinates
	 */
	public double distance(float x, float y);

	/**
	 * Calculate the XYZ distance to the given coordinates
	 */
	public double distance(float x, float y, float z);

	/**
	 * Calculate the XY squared distance to the given coordinates
	 */
	public double distance2(float x, float y);
	
	/**
	 * Calculate the XYZ squared distance to the given coordinates
	 */
	public double distance2(float x, float y, float z);
	
	/**
	 * Calculate the XY distance to the given coordinate
	 */
	public double distanceXY(Coordinate other);
	
	/**
	 * Calculate the XY squared distance to the given coordinate
	 */
	public double distanceXY2(Coordinate other);
	
	/**
	 * Calculate the XYZ distance to the given coordinate
	 */
	public double distanceXYZ(Coordinate other);
	
	/**
	 * Calculate the XYZ squared distance to the given coordinate
	 */
	public double distanceXYZ2(Coordinate other);
}