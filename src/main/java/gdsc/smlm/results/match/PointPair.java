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
 * Class to store a pair of coordinates.
 */
public class PointPair
{
	private Coordinate point1;
	private Coordinate point2;

	/**
	 * @param point1
	 * @param point2
	 */
	public PointPair(Coordinate point1, Coordinate point2)
	{
		this.point1 = point1;
		this.point2 = point2;
	}

	/**
	 * @return the point1
	 */
	public Coordinate getPoint1()
	{
		return point1;
	}

	/**
	 * @return the point2
	 */
	public Coordinate getPoint2()
	{
		return point2;
	}
	
	/**
	 * @return the distance (or -1 if either point is null)
	 */
	public double getXYZDistance()
	{
		if (point1 == null || point2 == null)
			return -1;
		
		return point1.distanceXYZ(point2);
	}
	
	/**
	 * @return the squared distance (or -1 if either point is null)
	 */
	public double getXYZDistance2()
	{
		if (point1 == null || point2 == null)
			return -1;
		
		return point1.distanceXYZ2(point2);
	}
	
	/**
	 * @return the XY distance (or -1 if either point is null)
	 */
	public double getXYDistance()
	{
		if (point1 == null || point2 == null)
			return -1;
		
		return point1.distanceXY(point2);
	}
	
	/**
	 * @return the squared XY distance (or -1 if either point is null)
	 */
	public double getXYDistance2()
	{
		if (point1 == null || point2 == null)
			return -1;
		
		return point1.distanceXY2(point2);
	}
}
