package gdsc.smlm.fitting;

import gdsc.core.utils.DoubleEquality;
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
	private double[][] m_inv = null;
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

		inverted = NO;

		m_inv = new double[n][];
		for (int i = n; i-- > 0;)
			m_inv[i] = m[i].clone();
		EJMLLinearSolver solver = new EJMLLinearSolver();
		solver.setEqual(equal);
		if (solver.invertSymmPosDef(m_inv))
		{
			// Check all diagonal values are zero or above
			for (int i = n; i-- > 0;)
				if (m_inv[i][i] < 0)
					return;

			inverted = YES;
		}
	}

	/**
	 * Compute the Cramer-Rao Lower Bound (CRLB) for fitted variables using the central diagonal of the inverted
	 * Fisher information matrix.
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
	public double[] crlb(boolean allowReciprocal)
	{
		invert();

		if (inverted == YES)
		{
			final double[] crlb = new double[n];
			for (int i = n; i-- > 0;)
				crlb[i] = Math.sqrt(m_inv[i][i]);
			return crlb;
		}

		if (allowReciprocal)
		{
			return crlbReciprocal();
		}

		return null;
	}

	/**
	 * Compute the Cramer-Rao Lower Bound (CRLB) for fitted variables using the reciprocal of the central diagonal of
	 * the Fisher information matrix.
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
			crlb[i] = reciprocalSqrt(m[i][i]);
		return crlb;
	}

	/**
	 * Compute the Cramer-Rao Lower Bound (CRLB) for fitted variables using the reciprocal of the central diagonal of
	 * the Fisher information matrix.
	 * 
	 * The information matrix is NOT inverted. The square root of the reciprocal of the central diagonal returned for a
	 * (possibly loose) lower bound.
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
			crlb[i] = reciprocalSqrt(m[i][i]);
		return crlb;
	}

	/**
	 * Return the reciprocal of the square root of the input. If the number is not strictly positive then zero is
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
}