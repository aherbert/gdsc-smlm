package gdsc.smlm.search;

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
 * Null implementation of the Dimension interface that does not perform rounding
 */
class NonRoundingDimension implements Dimension
{
	public double getLower()
	{
		return 0;
	}

	public double getUpper()
	{
		return 0;
	}

	public double getCentre()
	{
		return 0;
	}

	public double getMin()
	{
		return 0;
	}

	public double getMax()
	{
		return 0;
	}

	public boolean isActive()
	{
		return true;
	}

	public boolean isAtBounds(double v)
	{
		return false;
	}

	public Dimension create(double lower, double upper)
	{
		return null;
	}

	/**
	 * Does not round the number
	 * 
	 * @see gdsc.smlm.search.Dimension#round(double)
	 */
	public double round(double value)
	{
		return value;
	}

	public boolean canRound()
	{
		return true;
	}		
}