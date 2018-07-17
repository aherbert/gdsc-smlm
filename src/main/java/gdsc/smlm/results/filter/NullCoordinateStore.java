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
package gdsc.smlm.results.filter;

/**
 * Null implementation of the CoordinateStore interface.
 */
public class NullCoordinateStore implements CoordinateStore
{
	/** An instance to ignore calls to the CoordinateStore interface */
	public static final NullCoordinateStore INSTANCE = new NullCoordinateStore();

	/**
	 * Instantiates a new null coordinate store.
	 */
	private NullCoordinateStore()
	{
	}

	/**
	 * Creates an instance if the argument is null, else return the argument.
	 *
	 * @param coordinateStore
	 *            the coordinate store (may be null)
	 * @return the coordinate store (not null)
	 */
	public static CoordinateStore replaceIfNull(CoordinateStore coordinateStore)
	{
		return (coordinateStore == null) ? INSTANCE : coordinateStore;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.CoordinateStore#getResolution()
	 */
	@Override
	public double getXYResolution()
	{
		return 0;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.CoordinateStore#getZResolution()
	 */
	@Override
	public double getZResolution()
	{
		return 0;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.CoordinateStore#addToQueue(double, double, double)
	 */
	@Override
	public void addToQueue(double x, double y, double z)
	{
		// Do nothing
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.CoordinateStore#flush()
	 */
	@Override
	public void flush()
	{
		// Do nothing
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.CoordinateStore#add(double, double, double)
	 */
	@Override
	public void add(double x, double y, double z)
	{
		// Do nothing
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.CoordinateStore#clear()
	 */
	@Override
	public void clear()
	{
		// Do nothing
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.CoordinateStore#contains(double, double, double)
	 */
	@Override
	public boolean contains(double x, double y, double z)
	{
		return false;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.CoordinateStore#find(double, double, double)
	 */
	@Override
	public double[] find(double x, double y, double z)
	{
		return null;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.CoordinateStore#newInstance()
	 */
	@Override
	public CoordinateStore newInstance()
	{
		return this;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.CoordinateStore#resize(int, int, int, int)
	 */
	@Override
	public CoordinateStore resize(int minx, int miny, int maxx, int maxy)
	{
		return this;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.CoordinateStore#getMinX()
	 */
	@Override
	public int getMinX()
	{
		return 0;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.CoordinateStore#getMinY()
	 */
	@Override
	public int getMinY()
	{
		return 0;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.CoordinateStore#getWidth()
	 */
	@Override
	public int getWidth()
	{
		return 0;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.CoordinateStore#getHeight()
	 */
	@Override
	public int getHeight()
	{
		return 0;
	}
}
