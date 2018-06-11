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
package gdsc.smlm.engine;


/**
 * Performs quadrant analysis on fit residuals to look for asymmetry
 */
public class QuadrantAnalysis
{
	// Make these public for simplicity

	public double ABCD;
	public double A, B, C, D;
	public double ABCD2;
	public double A2, B2, C2, D2;

	public double AC;
	public double BD;
	public double score1;

	public double AC2;
	public double BD2;
	public double score2;

	public int[] vector;
	public double score;

	public double x1;
	public double y1;
	public double x2;
	public double y2;

	public int xi1;
	public int yi1;
	public int xi2;
	public int yi2;

	/**
	 * Perform quadrant analysis as per rapidSTORM
	 * <p>
	 * When two fluorophores emit close to each other, typically the nonlinear fit will result in a suspected
	 * fluorophore position midway between the two fluorophores and with a high amplitude. In this case, the fit
	 * results show a characteristic handle structure: The two true fluorophore emissions leave slightly positive
	 * residues, while there are negative residues on an axis perpendicular to the one connecting the fluorophores.
	 * <p>
	 * This condition is detected well by quadrant-differential residue analysis: The residue matrix is divided into
	 * quadrants, with the pixels above both diagonals forming the upper quadrant, the pixels above the main and
	 * below the off diagonal forming the right quadrants and so on. Pixels right on the diagonals are ignored.
	 * Then, the opposing quadrants are summed, and these sums substracted from another, resulting in two quadrant
	 * differences: upper and lower minus right and left and right and left minus upper and lower. This process is
	 * repeated for the quadrants defined by the central row and the central column.
	 * <p>
	 * The maximum sum obtained in this way divided by the sum of the absolute quadrant contributions is an
	 * observable correlating highly with the true presence of double emitters. Also, the quadrants containing the
	 * positive contribution in the highest sum indicate along which axis the double emission happened.
	 *
	 * @param residuals
	 *            the residuals
	 * @param width
	 *            the width
	 * @param height
	 *            the height
	 * @param cx
	 *            the centre in x
	 * @param cy
	 *            the centre in y
	 * @return true, if successful
	 */
	public boolean quadrantAnalysis(final double[] residuals, final int width, final int height, final int cx,
			final int cy)
	{
		vector = null;

		if (cx < 0 || cx >= width || cy < 0 || cy >= height) // out of bounds
			return false;

		// Compute quadrants

		// X quadrant:
		// .AAA.   
		// D.A.B   
		// DD.BB
		// D.C.B
		// .CCC.
		ABCD = 0;
		A = 0;
		B = 0;
		C = 0;
		D = 0;
		for (int y = cy, x1 = cx, x2 = cx; y < height; y++, x1--, x2++)
		{
			for (int x = 0, index = y * width; x < width; x++, index++)
			{
				ABCD += Math.abs(residuals[index]);
				if (x < x1)
					D += residuals[index];
				else if (x < x2 && x > x1)
					C += residuals[index];
				else if (x > x2)
					B += residuals[index];
				else
					ABCD -= Math.abs(residuals[index]);
			}
		}
		for (int y = cy - 1, x1 = cx - 1, x2 = cx + 1; y >= 0; y--, x1--, x2++)
		{
			for (int x = 0, index = y * width; x < width; x++, index++)
			{
				ABCD += Math.abs(residuals[index]);
				if (x < x1)
					D += residuals[index];
				else if (x < x2 && x > x1)
					A += residuals[index];
				else if (x > x2)
					B += residuals[index];
				else
					ABCD -= Math.abs(residuals[index]);
			}
		}

		// Similar for + quadrants:
		// AA.BB
		// AA.BB
		// .....
		// DD.CC
		// DD.CC
		ABCD2 = 0;
		A2 = 0;
		B2 = 0;
		C2 = 0;
		D2 = 0;
		for (int y = cy + 1; y < height; y++)
		{
			for (int x = 0, index = y * width; x < width; x++, index++)
			{
				ABCD2 += Math.abs(residuals[index]);
				if (x < cx)
					D2 += residuals[index];
				else if (x > cx)
					C2 += residuals[index];
			}
		}
		for (int y = cy - 1; y >= 0; y--)
		{
			for (int x = 0, index = y * width; x < width; x++, index++)
			{
				ABCD2 += Math.abs(residuals[index]);
				if (x < cx)
					A2 += residuals[index];
				else if (x > cx)
					B2 += residuals[index];
			}
		}

		// X quadrant:
		// .AAA.   
		// D.A.B   
		// DD.BB
		// D.C.B
		// .CCC.
		AC = A + C;
		BD = B + D;
		score1 = Math.abs(AC - BD) / ABCD;

		// + quadrant:
		// AA.BB
		// AA.BB
		// .....
		// DD.CC
		// DD.CC
		AC2 = A2 + C2;
		BD2 = B2 + D2;
		score2 = Math.abs(AC2 - BD2) / ABCD2;

		if (score1 > score2)
		{
			vector = (AC > BD) ? new int[] { 0, 1 } : new int[] { 1, 0 };
			score = score1;
		}
		else
		{
			vector = (AC2 > BD2) ? new int[] { 1, 1 } : new int[] { 1, -1 };
			score = score2;
		}

		return true;
	}

	/**
	 * Locate the 2 new centres by moving out into the quadrant defined by the computed vector by the defined shift
	 * <p>
	 * Requires a valid call to {@link #quadrantAnalysis(double[], int, int, int, int)} to create the vector
	 *
	 * @param width
	 *            the width
	 * @param height
	 *            the height
	 * @param cx
	 *            the centre in x
	 * @param cy
	 *            the centre in y
	 * @param shiftx
	 *            the shiftx
	 * @param shifty
	 *            the shifty
	 * @return true, if successful
	 */
	public boolean computeDoubletCentres(final int width, final int height, final int cx, final int cy, double shiftx,
			double shifty)
	{
		if (vector == null)
			return false;
		if (cx < 0 || cx >= width || cy < 0 || cy >= height) // out of bounds
			return false;

		// Pick double coords. The input centres must be shifted to the centre of the pixel by 0.5.
		x1 = cx + 0.5 + vector[0] * shiftx;
		y1 = cy + 0.5 + vector[1] * shifty;
		x2 = cx + 0.5 - vector[0] * shiftx;
		y2 = cy + 0.5 - vector[1] * shifty;

		// Check bounds
		if (x1 < 0)
			x1 = 0;
		else if (x1 >= width)
			x1 = width - 0.01;
		if (y1 < 0)
			y1 = 0;
		else if (y1 >= height)
			y1 = height - 0.01;
		if (x2 < 0)
			x2 = 0;
		else if (x2 >= width)
			x2 = width - 0.01;
		if (y2 < 0)
			y2 = 0;
		else if (y2 >= height)
			y2 = height - 0.01;

		xi1 = (int) (x1);
		yi1 = (int) (y1);
		xi2 = (int) (x2);
		yi2 = (int) (y2);

		// Check the two points are not the same pixel.
		// This is an edge case where the shift is small.
		if (xi1 == xi2 && yi1 == yi2)
		{
			// This can only happen when the shift is zero after rounding. 
			// If they are the same then the value (xi1,yi1) should be cx,cy
			// and we can move along the vector.
			x1 = cx + 0.5 + vector[0];
			y1 = cy + 0.5 + vector[1];
			x2 = cx + 0.5 - vector[0];
			y2 = cy + 0.5 - vector[1];

			// Check bounds again
			if (x1 < 0)
				x1 = 0;
			else if (x1 >= width)
				x1 = width - 0.01;
			if (y1 < 0)
				y1 = 0;
			else if (y1 >= height)
				y1 = height - 0.01;
			if (x2 < 0)
				x2 = 0;
			else if (x2 >= width)
				x2 = width - 0.01;
			if (y2 < 0)
				y2 = 0;
			else if (y2 >= height)
				y2 = height - 0.01;
			
			xi1 = (int) (x1);
			yi1 = (int) (y1);
			xi2 = (int) (x2);
			yi2 = (int) (y2);
			
			return (xi1 != xi2 || yi1 != yi2);
		}

		return true;
	}

	/**
	 * Gets the angle between two vectors
	 *
	 * @param a
	 *            the a
	 * @param b
	 *            the b
	 * @return the angle (in radians)
	 */
	public static double getAngle(int[] a, double[] b)
	{
		double d1 = a[0] * a[0] + a[1] * a[1];
		double d2 = b[0] * b[0] + b[1] * b[1];
		if (d1 > 0.0 && d2 > 0.0)
		{
			d1 = Math.sqrt(d1);
			d2 = Math.sqrt(d2);
			final double sum = a[0] * b[0] + a[1] * b[1];
			final double cosang = sum / (d1 * d2);

			if (cosang > 1.0)
			{
				return 0;
			}
			else if (cosang < -1.0)
			{
				return Math.PI;
			}

			return Math.acos(cosang);
		}
		return 999;
	}

	/**
	 * Gets the angle between two vectors
	 *
	 * @param a
	 *            the a
	 * @param b
	 *            the b
	 * @return the angle (in radians)
	 */
	public static double getAngle(double[] a, double[] b)
	{
		double d1 = a[0] * a[0] + a[1] * a[1];
		double d2 = b[0] * b[0] + b[1] * b[1];
		if (d1 > 0.0 && d2 > 0.0)
		{
			d1 = Math.sqrt(d1);
			d2 = Math.sqrt(d2);
			final double sum = a[0] * b[0] + a[1] * b[1];
			final double cosang = sum / (d1 * d2);

			if (cosang > 1.0)
			{
				return 0;
			}
			else if (cosang < -1.0)
			{
				return Math.PI;
			}

			return Math.acos(cosang);
		}
		return 999;
	}
}
