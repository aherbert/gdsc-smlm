package gdsc.smlm.fitting;

import org.ejml.data.DenseMatrix64F;

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
	/** The default inversion tolerance. */
	public static final double DEFAULT_INVERSION_TOLERANCE = 1e-2;

	private static final byte YES = 1;
	private static final byte UNKNOWN = 0;
	private static final byte NO = -1;

	private final DenseMatrix64F m;
	private double[] crlb = null;
	private byte inverted = UNKNOWN;
	private double inversionTolerance = 0;
	private boolean preserve = false;

	/**
	 * Instantiates a new fisher information matrix.
	 *
	 * @param m
	 *            the fisher information matrix
	 * @param inversionTolerance
	 *            the inversion tolerance
	 */
	public FisherInformationMatrix(double[][] m, double inversionTolerance)
	{
		this(EJMLLinearSolver.toA(m), inversionTolerance);
	}

	/**
	 * Instantiates a new fisher information matrix.
	 *
	 * @param m
	 *            the fisher information matrix
	 * @param n
	 *            the number of columns/rows
	 * @param inversionTolerance
	 *            the inversion tolerance
	 */
	public FisherInformationMatrix(double[] m, int n, double inversionTolerance)
	{
		this(EJMLLinearSolver.toA(m, n), inversionTolerance);
	}

	/**
	 * Instantiates a new fisher information matrix.
	 *
	 * @param m
	 *            the fisher information matrix
	 * @param inversionTolerance
	 *            the inversion tolerance
	 */
	public FisherInformationMatrix(DenseMatrix64F m, double inversionTolerance)
	{
		this.m = m;
		setInversionTolerance(inversionTolerance);
	}

	/**
	 * Instantiates a new fisher information matrix.
	 *
	 * @param m
	 *            the fisher information matrix
	 */
	public FisherInformationMatrix(double[][] m)
	{
		this(EJMLLinearSolver.toA(m));
	}

	/**
	 * Instantiates a new fisher information matrix.
	 *
	 * @param m
	 *            the fisher information matrix
	 * @param n
	 *            the number of columns/rows
	 */
	public FisherInformationMatrix(double[] m, int n)
	{
		this(EJMLLinearSolver.toA(m, n));
	}

	/**
	 * Instantiates a new fisher information matrix.
	 *
	 * @param m
	 *            the fisher information matrix
	 * @param n
	 *            the number of columns/rows
	 */
	public FisherInformationMatrix(DenseMatrix64F m)
	{
		this(m, DEFAULT_INVERSION_TOLERANCE);
	}

	private void invert()
	{
		if (inverted != UNKNOWN)
			return;

		if (m.numCols == 0)
		{
			// Nothing to do
			crlb = new double[0];
			inverted = YES;
			return;
		}

		inverted = NO;

		// Matrix inversion
		EJMLLinearSolver solver = EJMLLinearSolver.createForInversion(inversionTolerance);
		DenseMatrix64F a = (preserve && solver.invertModifiesA()) ? m.copy() : m;
		double[] crlb = solver.invertDiagonal(a);
		if (crlb != null)
		{
			// Check all diagonal values are zero or above
			if (inversionTolerance > 0)
			{
				// Already checked so just ignore values just below zero
				for (int i = m.numCols; i-- > 0;)
					if (crlb[i] < 0)
						crlb[i] = 0;
			}
			else
			{
				// A small error is OK
				for (int i = m.numCols; i-- > 0;)
				{
					if (crlb[i] < 0)
					{
						if (crlb[i] > -DEFAULT_INVERSION_TOLERANCE)
						{
							crlb[i] = 0;
							continue;
						}
						return;
					}
				}
			}

			// Check all diagonal values are zero or above

			inverted = YES;
			this.crlb = crlb;
		}
	}

	/**
	 * Compute the Cramer-Rao Lower Bound (CRLB) variance for fitted variables using the central diagonal of the
	 * inverted Fisher information matrix.
	 * <p>
	 * The information matrix is inverted and the square root of the central diagonal returned.
	 * <p>
	 * Warning: The matrix may be destroyed unless the preserve flag is set.
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
	 * <p>
	 * Warning: The matrix may be destroyed unless the preserve flag is set, or the reciprocal flag is true.
	 *
	 * @param allowReciprocal
	 *            the allow reciprocal flag
	 * @return CRLB (or null if inversion failed and the reciprocal is not used)
	 */
	public double[] crlb(boolean allowReciprocal)
	{
		// We must preserve the matrix if we are computing the reciprocal upon failure
		preserve |= allowReciprocal;
		
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
		final double[] crlb = new double[m.numCols];
		for (int i = 0, j = 0, n = m.numCols; i < n; i++, j += n + 1)
			crlb[i] = reciprocal(m.data[j]);
		return crlb;
	}

	/**
	 * Compute the Cramer-Rao Lower Bound (CRLB) for fitted variables using the central diagonal of the inverted
	 * Fisher information matrix.
	 * <p>
	 * The information matrix is inverted and the square root of the central diagonal returned.
	 * <p>
	 * Warning: The matrix may be destroyed unless the preserve flag is set.
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
	 * <p>
	 * Warning: The matrix may be destroyed unless the preserve flag is set, or the reciprocal flag is true.
	 *
	 * @param allowReciprocal
	 *            the allow reciprocal flag
	 * @return CRLB (or null if inversion failed and the reciprocal is not used)
	 */
	public double[] crlbSqrt(boolean allowReciprocal)
	{
		// We must preserve the matrix if we are computing the reciprocal upon failure
		preserve |= allowReciprocal;
		
		invert();

		if (inverted == YES)
		{
			final double[] crlb = new double[this.crlb.length];
			for (int i = crlb.length; i-- > 0;)
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
	 * for a (possibly loose) lower bound.
	 *
	 * @return CRLB (or null if inversion failed)
	 */
	public double[] crlbReciprocalSqrt()
	{
		final double[] crlb = new double[m.numCols];
		for (int i = 0, j = 0, n = m.numCols; i < n; i++, j += n + 1)
			crlb[i] = reciprocalSqrt(m.data[j]);
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
	 * Gets the inversion tolerance. Inversions are checked by ensuring that the product matches the identity matrix: A
	 * * A^-1 = I. Elements must be within the tolerance or else the inversion is rejected. Set to zero to disable.
	 *
	 * @return the inversion tolerance
	 */
	public double getInversionTolerance()
	{
		return inversionTolerance;
	}

	/**
	 * Sets the inversion tolerance. Inversions are checked by ensuring that the product matches the identity matrix: A
	 * * A^-1 = I. Elements must be within the tolerance or else the inversion is rejected. Set to zero to disable.
	 *
	 * @param inversionTolerance
	 *            the new inversion tolerance
	 */
	public void setInversionTolerance(double inversionTolerance)
	{
		this.inversionTolerance = inversionTolerance;
	}

	/**
	 * Gets a copy of the matrix.
	 *
	 * @return the matrix
	 */
	public double[][] getSquareMatrix()
	{
		return EJMLLinearSolver.toSquareData(m);
	}

	/**
	 * Gets a reference to the matrix.
	 *
	 * @return the matrix
	 */
	public DenseMatrix64F getMatrix()
	{
		return m;
	}

	/**
	 * Checks if the matrix will be preserved following inversion.
	 *
	 * @return true, if is preserve
	 */
	public boolean isPreserve()
	{
		return preserve;
	}

	/**
	 * Sets if the matrix will be preserved following inversion.
	 *
	 * @param preserve
	 *            the new preserve
	 */
	public void setPreserve(boolean preserve)
	{
		this.preserve = preserve;
	}
}