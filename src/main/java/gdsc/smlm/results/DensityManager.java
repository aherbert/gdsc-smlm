package gdsc.smlm.results;

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

import gdsc.smlm.fitting.function.Gaussian2DFunction;

import java.awt.Rectangle;

import org.apache.commons.math3.util.FastMath;

/**
 * Calculate the density of localisations around a given position using a square block of specified width
 */
public class DensityManager
{
	public class Molecule
	{
		int id;
		public float x, y;

		// Used to construct a single linked list of molecules
		public Molecule next = null;

		public Molecule(int id, float x, float y, Molecule next)
		{
			this.id = id;
			this.x = x;
			this.y = y;
			this.next = next;
		}

		public float distance(Molecule other)
		{
			final float dx = x - other.x;
			final float dy = y - other.y;
			return (float) Math.sqrt(dx * dx + dy * dy);
		}

		public float distance2(Molecule other)
		{
			final float dx = x - other.x;
			final float dy = y - other.y;
			return dx * dx + dy * dy;
		}
	}

	private TrackProgress tracker = null;
	private float[] xcoord, ycoord;
	private float minXCoord, minYCoord, maxXCoord, maxYCoord;
	private int area;

	/**
	 * @param results
	 * @throws IllegalArgumentException
	 *             if results are null or empty
	 */
	public DensityManager(final MemoryPeakResults results)
	{
		initialise(results);
	}

	/**
	 * Input arrays are modified
	 * 
	 * @param xcoord
	 * @param ycoord
	 * @param bounds
	 * @throws IllegalArgumentException
	 *             if results are null or empty
	 */
	public DensityManager(float[] xcoord, float[] ycoord, Rectangle bounds)
	{
		initialise(xcoord, ycoord, bounds);
	}

	private void initialise(final MemoryPeakResults results)
	{
		if (results == null || results.size() == 0)
			throw new IllegalArgumentException("Results are null or empty");

		xcoord = new float[results.size()];
		ycoord = new float[xcoord.length];
		for (int i = 0; i < xcoord.length; i++)
		{
			PeakResult result = results.getResults().get(i);
			xcoord[i] = result.params[Gaussian2DFunction.X_POSITION];
			ycoord[i] = result.params[Gaussian2DFunction.Y_POSITION];
		}

		initialise(xcoord, ycoord, results.getBounds());
	}

	private void initialise(float[] xcoord, float[] ycoord, Rectangle bounds)
	{
		if (xcoord == null || ycoord == null || xcoord.length == 0 || xcoord.length != ycoord.length)
			throw new IllegalArgumentException("Results are null or empty or mismatched in length");

		this.xcoord = xcoord;
		this.ycoord = ycoord;

		// Assign localisations & get min bounds
		minXCoord = Float.POSITIVE_INFINITY;
		minYCoord = Float.POSITIVE_INFINITY;
		for (int i = 0; i < xcoord.length; i++)
		{
			if (minXCoord > xcoord[i])
				minXCoord = xcoord[i];
			if (minYCoord > ycoord[i])
				minYCoord = ycoord[i];
		}

		// Round down and shift to origin
		minXCoord = (int) Math.floor(minXCoord);
		minYCoord = (int) Math.floor(minYCoord);
		// Get max bounds
		maxXCoord = 0;
		maxYCoord = 0;
		for (int i = 0; i < xcoord.length; i++)
		{
			xcoord[i] -= minXCoord;
			ycoord[i] -= minYCoord;
			if (maxXCoord < xcoord[i])
				maxXCoord = xcoord[i];
			if (maxYCoord < ycoord[i])
				maxYCoord = ycoord[i];
		}

		// Store the area of the input results
		area = bounds.width * bounds.height;
	}

	/**
	 * Calculate the density for the results.
	 * <p>
	 * A square block is used around each result of the specified radius. The results are assigned to a grid using a
	 * cell size of radius / resolution. The totals of each cell are then counted for the range +/- radius around each
	 * result.
	 * <p>
	 * If the block overlaps the border of the image the density will suffer from under-counting. The value can be
	 * optionally scaled using the fraction of the overlap area.
	 * <p>
	 * Note that the score is the number of molecules surrounding the given molecule, so the molecule itself is not
	 * counted.
	 * 
	 * @param radius
	 * @param resolution
	 * @param adjustForBorder
	 * @return
	 */
	public int[] calculateSquareDensity(float radius, int resolution, boolean adjustForBorder)
	{
		if (radius < 0)
			throw new IllegalArgumentException("Radius must be positive");
		if (resolution < 1)
			throw new IllegalArgumentException("Resolution must be positive");

		final float cellSize = radius / resolution;

		int maxx = (int) (maxXCoord / cellSize) + 1;
		int maxy = (int) (maxYCoord / cellSize) + 1;

		// Allocate counts to the cells
		int[] data = new int[maxx * maxy];
		for (int i = 0; i < xcoord.length; i++)
		{
			int x = (int) (xcoord[i] / cellSize);
			int y = (int) (ycoord[i] / cellSize);
			data[y * maxx + x]++;
		}

		// Create rolling sum table. Re-use the storage
		// First row
		int cs_ = 0; // Column sum
		for (int i = 0; i < maxx; i++)
		{
			cs_ += data[i];
			data[i] = cs_;
		}

		// Remaining rows:
		// sum = rolling sum of row + sum of row above
		for (int y = 1; y < maxy; y++)
		{
			int i = y * maxx;
			cs_ = 0;

			// Remaining columns
			for (int x = 0; x < maxx; x++, i++)
			{
				cs_ += data[i];
				data[i] = data[i - maxx] + cs_;
			}
		}

		// Used for debugging
		//		FileWriter out = null;
		//		try
		//		{
		//			out = new FileWriter("/tmp/check.txt");
		//		}
		//		catch (IOException e)
		//		{
		//			// Ignore
		//		}

		// For each localisation, compute the sum of counts within a square box radius
		final float area = 4 * resolution * resolution;
		int[] density = new int[xcoord.length];
		for (int i = 0; i < xcoord.length; i++)
		{
			int u = (int) (xcoord[i] / cellSize);
			int v = (int) (ycoord[i] / cellSize);

			// Note: Subtract 1 to discount the current localisation. Should this be done?
			int sum = -1;

			// Get the bounds
			int minU = u - resolution - 1;
			int maxU = FastMath.min(u + resolution, maxx - 1);
			int minV = v - resolution - 1;
			int maxV = FastMath.min(v + resolution, maxy - 1);

			// Compute sum from rolling sum using:
			// sum(u,v) = 
			// + s(maxU,maxV) 
			// - s(minU,maxV)
			// - s(maxU,minV)
			// + s(minU,minV)
			// Note: 
			// s(u,v) = 0 when either u,v < 0
			// s(u,v) = s(umax,v) when u>umax
			// s(u,v) = s(u,vmax) when v>vmax
			// s(u,v) = s(umax,vmax) when u>umax,v>vmax

			// + s(maxU,maxV) 
			int index = maxV * maxx + maxU;
			sum += data[index];

			boolean clipped = false;
			if (minU >= 0)
			{
				// - s(minU,maxV)
				index = maxV * maxx + minU;
				sum -= data[index];
			}
			else
			{
				clipped = true;
				minU = -1;
			}
			if (minV >= 0)
			{
				// - s(maxU,minV)
				index = minV * maxx + maxU;
				sum -= data[index];

				if (minU >= 0)
				{
					// + s(minU,minV)
					index = minV * maxx + minU;
					sum += data[index];
				}
			}
			else
			{
				clipped = true;
				minV = -1;
			}

			// Adjust for area
			if (adjustForBorder && clipped)
			{
				sum *= area / ((maxU - minU - 1) * (maxV - minV - 1));
			}

			density[i] = sum;

			//			// Check
			//			if (out != null)
			//			{
			//				int sum2 = 0;
			//				float xlower = xcoord[i] - radius;
			//				float xupper = xcoord[i] + radius;
			//				float ylower = ycoord[i] - radius;
			//				float yupper = ycoord[i] + radius;
			//				for (int j = 0; j < xcoord.length; j++)
			//				{
			//					if (j == i)
			//						continue;
			//					if (xcoord[j] < xlower || xcoord[j] > xupper)
			//						continue;
			//					if (ycoord[j] < ylower || ycoord[j] > yupper)
			//						continue;
			//					sum2++;
			//				}
			//
			//				try
			//				{
			//					out.write(String.format("%d %d\n", sum, sum2));
			//				}
			//				catch (IOException e)
			//				{
			//					// Just shutdown
			//					try
			//					{
			//						out.close();
			//					}
			//					catch (IOException e1)
			//					{
			//						// Ignore
			//					}
			//					out = null;
			//				}
			//			}
		}

		//		if (out != null)
		//		{
		//			try
		//			{
		//				out.close();
		//			}
		//			catch (IOException e)
		//			{
		//				// Ignore
		//			}
		//		}

		return density;
	}

	/**
	 * @return the tracker
	 */
	public TrackProgress getTracker()
	{
		return tracker;
	}

	/**
	 * @param tracker
	 *            the tracker to set
	 */
	public void setTracker(TrackProgress tracker)
	{
		this.tracker = tracker;
	}

	/**
	 * Calculate the density for the results.
	 * <p>
	 * A circle is used around each result of the specified radius and the number of neighbours counted for each result.
	 * <p>
	 * If the block overlaps the border of the image the density will suffer from under-counting. The value can be
	 * optionally scaled using the fraction of the overlap area.
	 * <p>
	 * Note that the score is the number of molecules surrounding the given molecule, so the molecule itself is not
	 * counted.
	 * 
	 * @param radius
	 * @param adjustForBorder
	 * @return
	 */
	public int[] calculateDensity(float radius, boolean adjustForBorder)
	{
		if (radius < 0)
			throw new IllegalArgumentException("Radius must be positive");

		// For each localisation, compute the sum of counts within a circle radius
		// TODO - Determine the optimum parameters to switch to using the grid method.
		int[] density = (xcoord.length < 200) ? calculateDensityTriangle(radius) : calculateDensityGrid(radius);

		// Adjust for area
		if (adjustForBorder)
		{
			// Boundary
			final float upperX = maxXCoord - radius;
			final float upperY = maxYCoord - radius;

			for (int i = 0; i < xcoord.length; i++)
			{
				int sum = density[i];
				final float x = xcoord[i];
				final float y = ycoord[i];

				// Calculate the area of the circle that has been missed
				// http://stackoverflow.com/questions/622287/area-of-intersection-between-circle-and-rectangle
				// Assume: Circle centre will be within the rectangle

				//   S1       S2       S3
				//
				//        |        |
				//    A1  |________|   A3      SA
				//        /        \
				//       /|   A2   |\        
				// -----/-|--------|-\-----
				//     |  |        |  |    
				//     |B1|   B2   |B3|        SB
				//     |  |        |  |
				// -----\-|--------|-/-----
				//       \|   C2   |/   C3     SC
				//   C1   \________/
				//        |        |

				// Note: A1,A3,C1,C3 are inside the circle
				// S1 = Slice 1, SA = Slice A, etc

				// Calculate if the upper/lower boundary of the rectangle slices the circle
				// -- Calculate the slice area using the formula for a segment
				// -- Check if the second boundary is slices the circle (i.e. a vertex is inside the circle)
				// ---- Calculate the corner section area to subtract from the overlapping slices
				// Missed = S1 + S3 + SA + SC - A1 - A3 - C1 - C3
				double S1 = 0, S3 = 0, SA = 0, SC = 0, A1 = 0, A3 = 0, C1 = 0, C3 = 0;

				// Note all coords are shifted the origin so simply compare the radius and the 
				// max bounds minus the radius

				if (x < radius)
				{
					S1 = getSegmentArea(radius, radius - x);
					if (y < radius)
					{
						A1 = getCornerArea(radius, x, y);
					}
					if (y > upperY)
					{
						C1 = getCornerArea(radius, x, maxYCoord - y);
					}
				}
				if (x > upperX)
				{
					float dx = maxXCoord - x;
					S1 = getSegmentArea(radius, radius - dx);
					if (y < radius)
					{
						A3 = getCornerArea(radius, dx, y);
					}
					if (y > upperY)
					{
						C3 = getCornerArea(radius, dx, maxYCoord - y);
					}
				}
				if (y < radius)
				{
					SA = getSegmentArea(radius, radius - y);
				}
				if (y > upperY)
				{
					float dy = maxYCoord - y;
					SC = getSegmentArea(radius, radius - dy);
				}

				double missed = S1 + S3 + SA + SC - A1 - A3 - C1 - C3;
				if (missed > 0)
				{
					double adjustment = area / (area - missed);
					//					if (missed > area)
					//					{
					//						System.out.printf("Ooops %f > %f\n", missed, area);
					//					}
					//					else
					//					{
					//						System.out.printf("increase %f @ %f %f\n", adjustment, x, y);
					//					}
					sum *= adjustment;
				}

				density[i] = sum;
			}
		}

		return density;
	}

	/**
	 * Calculate the density for the results using an all-vs-all analysis.
	 * <p>
	 * A circle is used around each result of the specified radius and the number of neighbours counted for each result.
	 * <p>
	 * If the block overlaps the border of the image the density will suffer from under-counting.
	 * <p>
	 * Note that the score is the number of molecules surrounding the given molecule, so the molecule itself is not
	 * counted.
	 * 
	 * @param radius
	 * @return
	 */
	public int[] calculateDensity(float radius)
	{
		final float r2 = radius * radius;
		int[] density = new int[xcoord.length];
		for (int i = 0; i < xcoord.length; i++)
		{
			int sum = density[i];
			final float x = xcoord[i];
			final float y = ycoord[i];
			for (int j = 0; j < xcoord.length; j++)
			{
				if (i == j)
					continue;
				final float dx = x - xcoord[j];
				final float dy = y - ycoord[j];
				if (dx * dx + dy * dy < r2)
				{
					sum++;
				}
			}

			density[i] = sum;
		}
		return density;
	}

	/**
	 * Calculate the density for the results using an all-vs-all analysis in the lower triangle of comparisons.
	 * <p>
	 * A circle is used around each result of the specified radius and the number of neighbours counted for each result.
	 * <p>
	 * If the block overlaps the border of the image the density will suffer from under-counting.
	 * <p>
	 * Note that the score is the number of molecules surrounding the given molecule, so the molecule itself is not
	 * counted.
	 * 
	 * @param radius
	 * @return
	 */
	public int[] calculateDensityTriangle(float radius)
	{
		final float r2 = radius * radius;
		int[] density = new int[xcoord.length];
		for (int i = 0; i < xcoord.length; i++)
		{
			int sum = density[i];
			final float x = xcoord[i];
			final float y = ycoord[i];
			for (int j = i + 1; j < xcoord.length; j++)
			{
				final float dx = x - xcoord[j];
				final float dy = y - ycoord[j];
				if (dx * dx + dy * dy < r2)
				{
					sum++;
					density[j]++;
				}
			}

			density[i] = sum;
		}
		return density;
	}

	/**
	 * Calculate the density for the results using a nearest neighbour cell grid analysis.
	 * <p>
	 * A circle is used around each result of the specified radius and the number of neighbours counted for each result.
	 * <p>
	 * If the block overlaps the border of the image the density will suffer from under-counting.
	 * <p>
	 * Note that the score is the number of molecules surrounding the given molecule, so the molecule itself is not
	 * counted.
	 * 
	 * @param radius
	 * @param adjustForBorder
	 * @return
	 */
	public int[] calculateDensityGrid(float radius)
	{
		int[] density = new int[xcoord.length];

		float minx = xcoord[0];
		float miny = ycoord[0];
		float maxx = minx, maxy = miny;
		for (int i = 0; i < xcoord.length; i++)
		{
			final float x = xcoord[i];
			final float y = ycoord[i];
			if (minx > x)
				minx = x;
			else if (maxx < x)
				maxx = x;
			if (miny > y)
				miny = y;
			else if (maxy < y)
				maxy = y;
		}

		// Assign to a grid
		final float binWidth = radius * 1.01f;
		final int nXBins = 1 + (int) ((maxx - minx) / binWidth);
		final int nYBins = 1 + (int) ((maxy - miny) / binWidth);
		Molecule[][] grid = new Molecule[nXBins][nYBins];
		for (int i = 0; i < xcoord.length; i++)
		{
			final float x = xcoord[i];
			final float y = ycoord[i];
			final int xBin = (int) ((x - minx) / binWidth);
			final int yBin = (int) ((y - miny) / binWidth);
			// Build a single linked list
			grid[xBin][yBin] = new Molecule(i, x, y, grid[xBin][yBin]);
		}

		Molecule[] neighbours = new Molecule[5];
		final float radius2 = radius * radius;
		for (int yBin = 0; yBin < nYBins; yBin++)
		{
			for (int xBin = 0; xBin < nXBins; xBin++)
			{
				for (Molecule m1 = grid[xBin][yBin]; m1 != null; m1 = m1.next)
				{
					// Build a list of which cells to compare up to a maximum of 4
					//      | 0,0  |  1,0
					// ------------+-----
					// -1,1 | 0,1  |  1,1

					int count = 0;
					neighbours[count++] = m1.next;

					if (yBin < nYBins - 1)
					{
						neighbours[count++] = grid[xBin][yBin + 1];
						if (xBin > 0)
							neighbours[count++] = grid[xBin - 1][yBin + 1];
					}
					if (xBin < nXBins - 1)
					{
						neighbours[count++] = grid[xBin + 1][yBin];
						if (yBin < nYBins - 1)
							neighbours[count++] = grid[xBin + 1][yBin + 1];
					}

					// Compare to neighbours
					while (count-- > 0)
					{
						for (Molecule m2 = neighbours[count]; m2 != null; m2 = m2.next)
						{
							if (m1.distance2(m2) < radius2)
							{
								density[m1.id]++;
								density[m2.id]++;
							}
						}
					}
				}
			}
		}

		return density;
	}

	/**
	 * Calculate the area of circular segment, a portion of a disk whose upper boundary is a (circular) arc and whose
	 * lower boundary is a chord making a central angle of theta radians.
	 * <p>
	 * See http://mathworld.wolfram.com/CircularSegment.html
	 * 
	 * @param R
	 *            the radius of the circle
	 * @param h
	 *            the height of the arced portion
	 * @return The area
	 */
	private double getSegmentArea(double R, double h)
	{
		return R * R * Math.acos((R - h) / R) - (R - h) * Math.sqrt(2 * R * h - h * h);
	}

	/**
	 * Get the area taken by a corner of a rectangle within a circle of radius R
	 * 
	 * @param R
	 *            the radius of the circle
	 * @param x
	 *            The corner X position
	 * @param y
	 *            The corner Y position
	 * @return The area
	 */
	private double getCornerArea(double R, double x, double y)
	{
		// 1 vertex is inside the circle: The sum of the areas of a circular segment and a triangle.

		//                            (x,y)
		//	    XXXXX                   XXXXXXXXX p2
		//     X     X       Triangle ->X     _-X
		//    X       X                 X   _-  X 
		//    X    +--X--+              X _-   X <- Circular segment 
		//     X   | X   |              X-  XXX 
		//	    XXXXX    |              XXXX
		//	       |     |             p1

		// Assume: circle at origin, x & y are positive, x^2 + y^2 < radius^2

		// Get the point p1 & p2
		double x1 = x;
		double y1 = otherSide(x, R);
		double x2 = otherSide(y, R);
		double y2 = y;

		// Calculate half the length of the chord cutting the circle between p1 & p2
		final double dx = x2 - x1;
		final double dy = y2 - y1;
		double halfChord = Math.sqrt(dx * dx + dy * dy);

		// Calculate the height of the arced portion
		double h = R - otherSide(halfChord, R);

		// Get the area as the circular segment plus the triangle
		return getSegmentArea(R, h) + 0.5 * dx * dy;
	}

	/**
	 * Returns a from a right angle triangle where a^2 + b^2 = c^2
	 * 
	 * @param b
	 * @param c
	 * @return a
	 */
	private double otherSide(double b, double c)
	{
		return Math.sqrt(c * c - b * b);
	}

	/**
	 * Compute Ripley's K-function.
	 * <p>
	 * See http://en.wikipedia.org/wiki/Spatial_descriptive_statistics#Ripley.27s_K_and_L_functions
	 * 
	 * @param radius
	 *            The radius
	 * @return The K-function score
	 */
	public double ripleysKFunction(double radius)
	{
		if (radius < 0)
			throw new IllegalArgumentException("Radius must be positive");

		// Count the number of points within the distance 
		int sum = calculateSumGrid((float) radius);

		// Normalise
		double scale = area / ((double) xcoord.length * (double) xcoord.length);
		double k = sum * scale;

		return k;
	}

	/**
	 * Calculate the number of pairs within the given radius.
	 * <p>
	 * The sum is over i<n, j<n, i!=j
	 * 
	 * @param radius
	 * @return
	 */
	public int calculateSum(float radius)
	{
		final float r2 = (float) (radius * radius);
		int sum = 0;
		for (int i = 0; i < xcoord.length; i++)
		{
			final float x = xcoord[i];
			final float y = ycoord[i];
			for (int j = i + 1; j < xcoord.length; j++)
			{
				final float dx = x - xcoord[j];
				final float dy = y - ycoord[j];
				if (dx * dx + dy * dy < r2)
				{
					sum++;
				}
			}
		}

		// Note that the sum should be computed over:
		//   i<n, j<n, i!=j
		// Thus it should be doubled to account for j iterating from zero.
		sum *= 2;
		return sum;
	}

	/**
	 * Calculate the number of pairs within the given radius using a nearest neighbour cell grid analysis.
	 * <p>
	 * The sum is over i<n, j<n, i!=j
	 * 
	 * @param radius
	 * @return
	 */
	public int calculateSumGrid(float radius)
	{
		int sum = 0;

		float minx = xcoord[0];
		float miny = ycoord[0];
		float maxx = minx, maxy = miny;
		for (int i = 0; i < xcoord.length; i++)
		{
			final float x = xcoord[i];
			final float y = ycoord[i];
			if (minx > x)
				minx = x;
			else if (maxx < x)
				maxx = x;
			if (miny > y)
				miny = y;
			else if (maxy < y)
				maxy = y;
		}

		// Assign to a grid
		final float binWidth = radius * 1.01f;
		final int nXBins = 1 + (int) ((maxx - minx) / binWidth);
		final int nYBins = 1 + (int) ((maxy - miny) / binWidth);
		Molecule[][] grid = new Molecule[nXBins][nYBins];
		for (int i = 0; i < xcoord.length; i++)
		{
			final float x = xcoord[i];
			final float y = ycoord[i];
			final int xBin = (int) ((x - minx) / binWidth);
			final int yBin = (int) ((y - miny) / binWidth);
			// Build a single linked list
			grid[xBin][yBin] = new Molecule(i, x, y, grid[xBin][yBin]);
		}

		Molecule[] neighbours = new Molecule[5];
		final float radius2 = radius * radius;
		for (int yBin = 0; yBin < nYBins; yBin++)
		{
			for (int xBin = 0; xBin < nXBins; xBin++)
			{
				for (Molecule m1 = grid[xBin][yBin]; m1 != null; m1 = m1.next)
				{
					// Build a list of which cells to compare up to a maximum of 4
					//      | 0,0  |  1,0
					// ------------+-----
					// -1,1 | 0,1  |  1,1

					int count = 0;
					neighbours[count++] = m1.next;

					if (yBin < nYBins - 1)
					{
						neighbours[count++] = grid[xBin][yBin + 1];
						if (xBin > 0)
							neighbours[count++] = grid[xBin - 1][yBin + 1];
					}
					if (xBin < nXBins - 1)
					{
						neighbours[count++] = grid[xBin + 1][yBin];
						if (yBin < nYBins - 1)
							neighbours[count++] = grid[xBin + 1][yBin + 1];
					}

					// Compare to neighbours
					while (count-- > 0)
					{
						for (Molecule m2 = neighbours[count]; m2 != null; m2 = m2.next)
						{
							if (m1.distance2(m2) < radius2)
							{
								sum++;
							}
						}
					}
				}
			}
		}

		return sum * 2;
	}

	/**
	 * Compute Ripley's L-function.
	 * <p>
	 * See http://en.wikipedia.org/wiki/Spatial_descriptive_statistics#Ripley.27s_K_and_L_functions
	 * 
	 * @param radius
	 *            The radius
	 * @return The L-function score
	 */
	public double ripleysLFunction(double radius)
	{
		double k = ripleysKFunction(radius);
		return Math.sqrt(k / Math.PI);
	}

	/**
	 * Compute Ripley's K-function.
	 * <p>
	 * See http://en.wikipedia.org/wiki/Spatial_descriptive_statistics#Ripley.27s_K_and_L_functions
	 * 
	 * @param density
	 *            The density score for each particle
	 * @param radius
	 *            The radius at which the density was computed
	 * @return The K-function score
	 * @see {@link #calculateDensity(float, boolean)}
	 * @see {@link #calculateSquareDensity(float, int, boolean)}
	 */
	public double ripleysKFunction(int[] density, double radius)
	{
		if (radius < 0)
			throw new IllegalArgumentException("Radius must be positive");
		if (density.length != xcoord.length)
			throw new IllegalArgumentException("Input density array must match the number of coordinates");

		// Count the number of points within the distance 
		int sum = 0;
		for (int d : density)
		{
			sum += d;
		}

		// Normalise
		double scale = area / ((double) density.length * (double) density.length);
		double k = sum * scale;

		return k;
	}

	/**
	 * Compute Ripley's L-function.
	 * <p>
	 * See http://en.wikipedia.org/wiki/Spatial_descriptive_statistics#Ripley.27s_K_and_L_functions
	 * 
	 * @param density
	 *            The density score for each particle
	 * @param radius
	 *            The radius at which the density was computed
	 * @return The K-function score
	 * @see {@link #calculateDensity(float, boolean)}
	 * @see {@link #calculateSquareDensity(float, int, boolean)}
	 */
	public double ripleysLFunction(int[] density, double radius)
	{
		double k = ripleysKFunction(density, radius);
		return Math.sqrt(k / Math.PI);
	}
}
