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
package gdsc.smlm.search;

/**
 * Null implementation of the Dimension interface that does not perform rounding
 */
class NonRoundingDimension implements Dimension
{
	@Override
	public double getLower()
	{
		return 0;
	}

	@Override
	public double getUpper()
	{
		return 0;
	}

	@Override
	public double getCentre()
	{
		return 0;
	}

	@Override
	public double getMin()
	{
		return 0;
	}

	@Override
	public double getMax()
	{
		return 0;
	}

	@Override
	public boolean isActive()
	{
		return true;
	}

	@Override
	public boolean isAtBounds(double v)
	{
		return false;
	}

	@Override
	public Dimension create(double lower, double upper)
	{
		return null;
	}

	/**
	 * Does not round the number
	 * 
	 * @see gdsc.smlm.search.Dimension#round(double)
	 */
	@Override
	public double round(double value)
	{
		return value;
	}

	@Override
	public boolean canRound()
	{
		return true;
	}
}
