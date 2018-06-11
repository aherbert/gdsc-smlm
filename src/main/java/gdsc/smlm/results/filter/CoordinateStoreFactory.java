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
 * Creates a coordinate store
 */
public class CoordinateStoreFactory
{
	/**
	 * Creates the coordinate store.
	 *
	 * @param minx
	 *            the min x coordinate value
	 * @param miny
	 *            the min y coordinate value
	 * @param width
	 *            the width
	 * @param height
	 *            the height
	 * @param xyResolution
	 *            the xy resolution (if negative then nothing is stored)
	 * @return the coordinate store
	 * @deprecated The z resolution should be specified using {@link #create(int, int, double, double)}
	 */
	@Deprecated
	public static CoordinateStore create(int minx, int miny, int width, int height, double xyResolution)
	{
		return create(minx, miny, width, height, xyResolution, -1);
	}

	/**
	 * Creates the coordinate store.
	 *
	 * @param minx
	 *            the min x coordinate value
	 * @param miny
	 *            the min y coordinate value
	 * @param width
	 *            the width
	 * @param height
	 *            the height
	 * @param xyResolution
	 *            the xy resolution (if negative then nothing is stored)
	 * @param zResolution
	 *            the z resolution (if negative then this is ignored and the store behaves as if processing 2D
	 *            coordinates)
	 * @return the coordinate store
	 */
	public static CoordinateStore create(int minx, int miny, int width, int height, double xyResolution,
			double zResolution)
	{
		if (xyResolution < 0)
			return NullCoordinateStore.INSTANCE;

		// This should be faster (for additions and block lookup) as it has a fixed block resolution of 1.
		// However it may be slower if the distance is much lower than 1 and there are many points close 
		// to the resolution distance as it will have to compute the distance for each. As a compromise
		// we only use it when the resolution is above the min block size of the default store.
		if (xyResolution >= GridCoordinateStore.MINIMUM_BLOCK_SIZE && xyResolution <= 1)
			return new GridCoordinateStore1(minx, miny, width, height, xyResolution, zResolution);

		return new GridCoordinateStore(minx, miny, width, height, xyResolution, zResolution);
	}
}
