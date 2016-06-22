package gdsc.smlm.engine;

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
 * Performs quadrant analysis on fit residuals to look for asymmetry
 */
public class QuadrantAnalysis
{
	double ABCD = 0;
	double A = 0, B = 0, C = 0, D = 0;
	double ABCD2 = 0;
	double A2 = 0, B2 = 0, C2 = 0, D2 = 0;

	double AC = A + C;
	double BD = B + D;
	double score1;

	double AC2;
	double BD2;
	double score2;

	int[] vector;
	double score;

	int x1;
	int y1;
	int x2;
	int y2;

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
	 * Locate the 2 new centres by moving out into the quadrant defined by the vector
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

		x1 = (int) Math.round(cx + vector[0] * shiftx);
		y1 = (int) Math.round(cy + vector[1] * shifty);
		x2 = (int) Math.round(cx - vector[0] * shiftx);
		y2 = (int) Math.round(cy - vector[1] * shifty);

		// Check bounds
		if (x1 < 0)
			x1 = 0;
		else if (x1 >= width)
			x1 = width - 1;
		if (y1 < 0)
			y1 = 0;
		else if (y1 >= height)
			y1 = height - 1;
		if (x2 < 0)
			x2 = 0;
		else if (x2 >= width)
			x2 = width - 1;
		if (y2 < 0)
			y2 = 0;
		else if (y2 >= height)
			y2 = height - 1;

		// Check the two points are not the same. 
		if (x1 == x2 && y1 == y2)
		{
			// This can only happen when the shift is zero due to the round() function 
			// moving to the next integer. If they are the same then the value should be cx,cy
			// and we can move along the vector.
			x1 += vector[0];
			y1 += vector[1];
			x2 -= vector[0];
			y2 -= vector[1];

			// Check bounds
			if (x1 < 0)
				x1 = 0;
			else if (x1 >= width)
				x1 = width - 1;
			if (y1 < 0)
				y1 = 0;
			else if (y1 >= height)
				y1 = height - 1;
			if (x2 < 0)
				x2 = 0;
			else if (x2 >= width)
				x2 = width - 1;
			if (y2 < 0)
				y2 = 0;
			else if (y2 >= height)
				y2 = height - 1;

			return (x1 != x2 || y1 != y2);
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
