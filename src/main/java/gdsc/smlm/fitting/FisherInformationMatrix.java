package gdsc.smlm.fitting;

import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.Maths;
import gdsc.smlm.fitting.linear.EJMLLinearSolver;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Container for the Fisher information, a symmetric positive definite matrix containing the amount of information that
 * an observable random variable X carries about an unknown parameter θ of a distribution that models X.
 */
public class FisherInformationMatrix
{
	private static final byte YES = 1;
	private static final byte UNKNOWN = 0;
	private static final byte NO = -1;

	private final double[][] m;
	private final int n;
	private double[] crlb = null;
	private byte inverted = UNKNOWN;
	private DoubleEquality equal = null;

	/**
	 * Instantiates a new fisher information matrix.
	 *
	 * @param m
	 *            the fisher information matrix
	 */
	public FisherInformationMatrix(double[][] m)
	{
		this.m = m;
		n = m.length;
	}

	private void invert()
	{
		if (inverted != UNKNOWN)
			return;

		if (n == 0)
		{
			// Nothing to do
			crlb = new double[0];
			inverted = YES;
			return;
		}

		inverted = NO;

		//		if (n < 5)
		//		{
		//			// We could solve this without linear algebra...
		//			// Look for zero diagonal entries since the EJMLLinearSolver will handle those 
		//			boolean ok = true;
		//			for (int i = 0; i < n; i++)
		//			{
		//				if (m[i][i] == 0)
		//				{
		//					ok = false;
		//					break;
		//				}
		//			}
		//			if (ok)
		//			{
		//				// Most likely to be used with larger matrices
		//				if (n == 4)
		//				{
		//					crlb = computeCRLB4(m);
		//					if (crlb != null)
		//						inverted = YES;
		//					return;
		//				}
		//				if (n == 3)
		//				{
		//					crlb = computeCRLB3(m);
		//					if (crlb != null)
		//						inverted = YES;
		//					return;
		//				}
		//				if (n == 2)
		//				{
		//					crlb = computeCRLB2(m);
		//					if (crlb != null)
		//						inverted = YES;
		//					return;
		//				}
		//				if (n == 1)
		//				{
		//					if (m[0][0] != 0)
		//					{
		//						crlb = new double[] { 1.0 / m[0][0] };
		//						inverted = YES;
		//					}
		//					return;
		//				}
		//			}
		//		}

		// Matrix inversion
		EJMLLinearSolver solver = new EJMLLinearSolver(equal);
		// TODO - Use a DenseMatrix64F
		//double[] crlb = solver.invertSymmPosDefDiagonal(new DenseMatrix64F(data));
		double[] crlb = solver.invertSymmPosDefDiagonal(getMatrix());
		if (crlb != null)
		{
			// Check all diagonal values are zero or above
			for (int i = n; i-- > 0;)
			{
				if (crlb[i] < 0)
				{
					// A small error is OK
					if (crlb[i] > -1e-2 || (equal != null && equal.almostEqualComplement(crlb[i], 0)))
					{
						crlb[i] = 0;
						continue;
					}
					return;
				}
			}

			inverted = YES;
			this.crlb = crlb;
		}
	}

	/**
	 * Compute CRLB of the Fisher Information matrix by direct matrix inversion of a 2x2 matrix.
	 *
	 * @param m
	 *            the matrix
	 * @return the CRLB (or null if the determinant is zero)
	 */
	public static double[] computeCRLB2(double[][] m)
	{
		// https://en.wikipedia.org/wiki/Invertible_matrix#Inversion_of_2_.C3.97_2_matrices		
		double a = m[0][0];
		double b = m[0][1];
		double c = m[1][0];
		double d = m[1][1];

		double det = a * d - b * c;

		double det_recip = 1.0 / det;

		if (!Maths.isFinite(det_recip))
			return null;
		//@formatter:off
		return new double[] {
			det_recip * d,
			det_recip * a
		};
		//@formatter:on
	}

	/**
	 * Compute CRLB of the Fisher Information matrix by direct matrix inversion of a 3x3 matrix.
	 *
	 * @param m
	 *            the matrix
	 * @return the CRLB (or null if the determinant is zero)
	 */
	public static double[] computeCRLB3(double[][] m)
	{
		// https://en.wikipedia.org/wiki/Invertible_matrix#Inversion_of_3_.C3.97_3_matrices		
		double a = m[0][0];
		double b = m[0][1];
		double c = m[0][2];
		double d = m[1][0];
		double e = m[1][1];
		double f = m[1][2];
		double g = m[2][0];
		double h = m[2][1];
		double i = m[2][2];

		double A = (e * i - f * h);
		double B = -(d * i - f * g);
		double C = (d * h - e * g);

		double det = a * A + b * B + c * C;

		double det_recip = 1.0 / det;

		if (!Maths.isFinite(det_recip))
			return null;
		//@formatter:off
		return new double[] {
			det_recip * A,
			det_recip * (a * i - c * g),
			det_recip * (a * e - b * d)
		};
		//@formatter:on
	}

	/**
	 * Compute CRLB of the Fisher Information matrix by direct matrix inversion of a 4x4 matrix.
	 *
	 * @param m
	 *            the matrix
	 * @return the CRLB (or null if the determinant is zero)
	 */
	public static double[] computeCRLB4(double[][] m)
	{
		// Adapted from GraspJ:
		// https://github.com/isman7/graspj/blob/master/graspj/src/main/java/eu/brede/graspj/opencl/src/functions/CRLB_from_fisher_matrix_4x4.cl

		double a0 = m[0][0] * m[1][1] - m[0][1] * m[1][0];
		double a1 = m[0][0] * m[1][2] - m[0][2] * m[1][0];
		double a2 = m[0][0] * m[1][3] - m[0][3] * m[1][0];
		double a3 = m[0][1] * m[1][2] - m[0][2] * m[1][1];
		double a4 = m[0][1] * m[1][3] - m[0][3] * m[1][1];
		double a5 = m[0][2] * m[1][3] - m[0][3] * m[1][2];
		double b0 = m[2][0] * m[3][1] - m[2][1] * m[3][0];
		double b1 = m[2][0] * m[3][2] - m[2][2] * m[3][0];
		double b2 = m[2][0] * m[3][3] - m[2][3] * m[3][0];
		double b3 = m[2][1] * m[3][2] - m[2][2] * m[3][1];
		double b4 = m[2][1] * m[3][3] - m[2][3] * m[3][1];
		double b5 = m[2][2] * m[3][3] - m[2][3] * m[3][2];

		double det = a0 * b5 - a1 * b4 + a2 * b3 + a3 * b2 - a4 * b1 + a5 * b0;

		double det_recip = 1.0 / det;

		if (!Maths.isFinite(det_recip))
			return null;
		//@formatter:off
		return new double[] {
			det_recip * (m[1][1] * b5 - m[1][2] * b4 + m[1][3] * b3),
			det_recip * (m[0][0] * b5 - m[0][2] * b2 + m[0][3] * b1),
			det_recip * (m[3][0] * a4 - m[3][1] * a2 + m[3][3] * a0),
			det_recip * (m[2][0] * a3 - m[2][1] * a1 + m[2][2] * a0)
		};
		//@formatter:on
	}

	/**
	 * Compute the Cramer-Rao Lower Bound (CRLB) variance for fitted variables using the central diagonal of the
	 * inverted Fisher information matrix.
	 * <p>
	 * The information matrix is inverted and the square root of the central diagonal returned.
	 * 
	 * @return CRLB (or null if inversion failed)
	 */
	public double[] crlb()
	{
		return crlb(false);
	}

	/**
	 * Compute the Cramer-Rao Lower Bound (CRLB) variance for fitted variables using the central diagonal of the
	 * inverted Fisher information matrix.
	 * <p>
	 * The information matrix is inverted and the square root of the central diagonal returned. If the inversion fails
	 * then the routine optionally returns the square root of the reciprocal of the diagonal element to find a (possibly
	 * loose) lower bound.
	 *
	 * @param allowReciprocal
	 *            the allow reciprocal flag
	 * @return CRLB (or null if inversion failed and the reciprocal is not used)
	 */
	public double[] crlb(boolean allowReciprocal)
	{
		invert();

		if (inverted == YES)
		{
			return crlb;
		}

		if (allowReciprocal)
		{
			return crlbReciprocal();
		}

		return null;
	}

	/**
	 * Compute the Cramer-Rao Lower Bound (CRLB) variance for fitted variables using the reciprocal of the central
	 * diagonal of the Fisher information matrix.
	 * 
	 * The information matrix is NOT inverted. The square root of the reciprocal of the central diagonal returned for a
	 * (possibly loose) lower bound.
	 *
	 * @return CRLB (or null if inversion failed)
	 */
	public double[] crlbReciprocal()
	{
		final double[] crlb = new double[n];
		for (int i = 0; i < n; i++)
			crlb[i] = reciprocal(m[i][i]);
		return crlb;
	}

	/**
	 * Compute the Cramer-Rao Lower Bound (CRLB) for fitted variables using the central diagonal of the inverted
	 * Fisher information matrix.
	 * <p>
	 * The information matrix is inverted and the square root of the central diagonal returned.
	 * 
	 * @return CRLB (or null if inversion failed)
	 */
	public double[] crlbSqrt()
	{
		return crlbSqrt(false);
	}

	/**
	 * Compute the Cramer-Rao Lower Bound (CRLB) for fitted variables using the central diagonal of the inverted
	 * Fisher information matrix.
	 * <p>
	 * The information matrix is inverted and the square root of the central diagonal returned. If the inversion fails
	 * then the routine optionally returns the square root of the reciprocal of the diagonal element to find a (possibly
	 * loose) lower bound.
	 *
	 * @param allowReciprocal
	 *            the allow reciprocal flag
	 * @return CRLB (or null if inversion failed and the reciprocal is not used)
	 */
	public double[] crlbSqrt(boolean allowReciprocal)
	{
		invert();

		if (inverted == YES)
		{
			final double[] crlb = new double[n];
			for (int i = 0; i < n; i++)
				crlb[i] = Math.sqrt(this.crlb[i]);
			return crlb;
		}

		if (allowReciprocal)
		{
			return crlbReciprocalSqrt();
		}

		return null;
	}

	/**
	 * Compute the Cramer-Rao Lower Bound (CRLB) for fitted variables using the reciprocal of the central diagonal of
	 * the Fisher information matrix.
	 * 
	 * The information matrix is NOT inverted. Uses the square root of the reciprocal of the central diagonal returned
	 * for a
	 * (possibly loose) lower bound.
	 *
	 * @return CRLB (or null if inversion failed)
	 */
	public double[] crlbReciprocalSqrt()
	{
		final double[] crlb = new double[n];
		for (int i = 0; i < n; i++)
			crlb[i] = reciprocalSqrt(m[i][i]);
		return crlb;
	}

	/**
	 * Compute the Cramer-Rao Lower Bound (CRLB) variance for fitted variables using the reciprocal of the central
	 * diagonal of the Fisher information matrix.
	 * 
	 * The information matrix is NOT inverted. Uses the reciprocal of the central diagonal returned for a (possibly
	 * loose) lower bound.
	 *
	 * @param m
	 *            the fisher information matrix
	 * @return CRLB
	 */
	public static double[] crlbReciprocal(double[][] m)
	{
		int n = m.length;
		final double[] crlb = new double[n];
		for (int i = 0; i < n; i++)
			crlb[i] = reciprocal(m[i][i]);
		return crlb;
	}

	/**
	 * Compute the Cramer-Rao Lower Bound (CRLB) for fitted variables using the reciprocal of the central
	 * diagonal of the Fisher information matrix.
	 * 
	 * The information matrix is NOT inverted. Uses the square root of the reciprocal of the central diagonal returned
	 * for a (possibly loose) lower bound.
	 *
	 * @param m
	 *            the fisher information matrix
	 * @return CRLB
	 */
	public static double[] crlbReciprocalSqrt(double[][] m)
	{
		int n = m.length;
		final double[] crlb = new double[n];
		for (int i = 0; i < n; i++)
			crlb[i] = reciprocalSqrt(m[i][i]);
		return crlb;
	}

	/**
	 * Return the reciprocal of the input. If the number is not strictly positive then zero is
	 * returned. Note that zero can only be returned if there was no Fisher information. This is done to match the
	 * return value from matrix inversion when there is no Fisher information for a parameter i within the matrix. In
	 * that case the zero column and row is removed from the matrix before inversion and the inverted matrix contains
	 * zeros.
	 * <p>
	 * The reciprocal of the diagonal element of the Fisher information matrix is a (possibly loose) lower bound.
	 *
	 * @param d
	 *            the input value
	 * @return the reciprocal of the square root of the input value
	 */
	public static double reciprocal(double d)
	{
		return (d > 0) ? 1.0 / d : 0;
	}

	/**
	 * Return the reciprocal of the square root of the input. If the number is not strictly positive then zero is
	 * returned. Note that zero can only be returned if there was no Fisher information. This is done to match the
	 * return value from matrix inversion when there is no Fisher information for a parameter i within the matrix. In
	 * that case the zero column and row is removed from the matrix before inversion and the inverted matrix contains
	 * zeros.
	 * <p>
	 * The square root of the reciprocal of the diagonal element of the Fisher information matrix is a (possibly loose)
	 * lower bound.
	 *
	 * @param d
	 *            the input value
	 * @return the reciprocal of the square root of the input value
	 */
	public static double reciprocalSqrt(double d)
	{
		return (d > 0) ? 1.0 / Math.sqrt(d) : 0;
	}

	/**
	 * @param equal
	 *            the equality class to compare that the solution inversion is correct (A*A_inv = I)
	 */
	public void setEqual(DoubleEquality equal)
	{
		this.equal = equal;
	}

	/**
	 * @return the equality class to compare that the solution inversion is correct (A*A_inv = I)
	 */
	public DoubleEquality getEqual()
	{
		return equal;
	}

	/**
	 * Gets a copy of the matrix.
	 *
	 * @return the matrix
	 */
	public double[][] getMatrix()
	{
		double[][] m2 = new double[n][];
		for (int i = n; i-- > 0;)
			m2[i] = m[i].clone();
		return m2;
	}
}