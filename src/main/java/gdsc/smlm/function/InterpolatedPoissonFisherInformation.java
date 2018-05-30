package gdsc.smlm.function;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NonMonotonicSequenceException;
import org.apache.commons.math3.exception.NumberIsTooSmallException;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2018 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Calculate the Fisher information for a Poisson-distributed random variable using an interpolation
 * of the alpha scale parameter. The alpha scale parameter is the ratio between the Fisher information
 * for a Poisson distribution and the Fisher information of another Poisson-based distribution,
 * e.g. a Poisson-Gaussian convolution.
 */
public class InterpolatedPoissonFisherInformation extends BasePoissonFisherInformation
{
	/** The minimum of the interpolation range (log scale). */
	public final double min;

	/** The maximum of the interpolation range (log scale). */
	public final double max;

	/** The mean at the minimum of the interpolation range. */
	public final double uMin;

	/** The mean at the maximum of the interpolation range. */
	public final double uMax;

	/** The alpha at the minimum of the interpolation range. */
	private final double alphaMin;

	/** The Fisher information at the minimum of the interpolation range. */
	private final double iMin;

	/** The alpha at the maximum of the interpolation range. */
	private final double alphaMax;

	/** Flag indicating if the Fisher information or the alpha is fixed at the minimum of the interpolation range. */
	private final boolean lowerFixedI;

	/** The function to compute the Fisher information above the maximum of the interpolation range. */
	private BasePoissonFisherInformation upperFI;

	/** The function to interpolate alpha in the range min-max. */
	private final PolynomialSplineFunction alphaF;

	/** The fast log function. */
	private FastLog fastLog;

	/**
	 * Instantiates a new interpolated poisson fisher information.
	 * The input series of means must have at least 3 points and be in increasing order.
	 *
	 * @param logU
	 *            the log of the mean
	 * @param alpha
	 *            the alpha for each Poisson mean
	 * @throws DimensionMismatchException
	 *             if {@code x} and {@code y}
	 *             have different sizes.
	 * @throws NumberIsTooSmallException
	 *             if the size of {@code x} is smaller
	 *             than 3.
	 * @throws NonMonotonicSequenceException
	 *             if {@code x} is not sorted in
	 *             strict increasing order.
	 * @throws IllegalArgumentException
	 *             the illegal argument exception
	 */
	public InterpolatedPoissonFisherInformation(double[] logU, double[] alpha)
			throws DimensionMismatchException, NumberIsTooSmallException, NonMonotonicSequenceException
	{
		this(logU, alpha, true, null);
	}

	/**
	 * Instantiates a new interpolated poisson fisher information.
	 * The input series of means must have at least 3 points and be in increasing order.
	 *
	 * @param logU
	 *            the log of the mean
	 * @param alpha
	 *            the alpha for each Poisson mean
	 * @param lowerFixedI
	 *            Flag indicating if the Fisher information or the alpha is fixed at the minimum of the interpolation
	 *            range.
	 * @param upperFI
	 *            The function to compute the Fisher information above the maximum of the interpolation range. If null
	 *            then the alpha is considered fixed.
	 * @throws DimensionMismatchException
	 *             if {@code x} and {@code y}
	 *             have different sizes.
	 * @throws NumberIsTooSmallException
	 *             if the size of {@code x} is smaller
	 *             than 3.
	 * @throws NonMonotonicSequenceException
	 *             if {@code x} is not sorted in
	 *             strict increasing order.
	 * @throws IllegalArgumentException
	 *             the illegal argument exception
	 */
	public InterpolatedPoissonFisherInformation(double[] logU, double[] alpha, boolean lowerFixedI,
			BasePoissonFisherInformation upperFI)
			throws DimensionMismatchException, NumberIsTooSmallException, NonMonotonicSequenceException
	{
		SplineInterpolator si = new SplineInterpolator();
		alphaF = si.interpolate(logU, alpha);

		this.lowerFixedI = lowerFixedI;
		this.upperFI = upperFI;

		min = logU[0];
		alphaMin = alpha[0];

		int n_1 = logU.length - 1;
		max = logU[n_1];
		alphaMax = alpha[n_1];

		// Store the ends in non-log format
		uMin = Math.exp(min);
		uMax = Math.exp(max);

		iMin = alphaMin / uMin;

		fastLog = FastLogFactory.getFastLog();
	}

	/**
	 * {@inheritDoc}
	 * <p>
	 * Gets the approximate Poisson-Gaussian Fisher information.
	 * Approximate the Poisson as a Gaussian (u=t, var=t) and convolve with a Gaussian (u=0,var=s*s).
	 * Gaussian-Gaussian convolution: var1 * var2 => var = var1+var2.
	 * The Fisher information of Gaussian mean is 1/variance.
	 * The Poisson-Gaussian Fisher information is therefore 1 / (t + s*s).
	 * 
	 * @see gdsc.smlm.function.FisherInformation#getFisherInformation(double)
	 */
	public double getFisherInformation(double t) throws IllegalArgumentException
	{
		if (t <= 0)
			throw new IllegalArgumentException("Poisson mean must be positive");
		// Poisson fisher information
		double I = 1.0 / t;
		// The Fisher information is returned using the alpha multiplied by the 
		// Poisson Fisher information.
		if (I != Double.POSITIVE_INFINITY)
			I *= getAlpha(t);
		return I;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.BasePoissonFisherInformation#getAlpha(double)
	 */
	public double getAlpha(double t)
	{
		if (t <= 0)
			throw new IllegalArgumentException("Poisson mean must be positive");

		// At the minimum the Fisher information uses a fixed alpha.
		if (t <= uMin)
			return getAlphaMin(t);

		// At the maximum the Fisher information can use a fixed alpha or a approximation
		// function.
		if (t >= uMax)
			return getAlphaMax(t);

		// Within the range the poisson mean is converted to a log scale and alpha is
		// interpolated. Use a fast log for this as the precision is not critical due 
		// to the assumed error in the interpolation. 
		// At this point t is known to be in the bound >0, but it may be NaN or infinity
		// so allow the checks (i.e. don't use fastLogD(double)).
		final double x = fastLog.logD(t);

		// Check again as fast log may not be precise. 
		// This avoids an out-of-range exception in the interpolating function. 
		if (x <= min)
			return getAlphaMin(t);
		if (x >= max)
			return getAlphaMax(t);

		return alphaF.value(x);
	}

	private double getAlphaMin(double t)
	{
		// alpha = t * I
		return (lowerFixedI) ? t * iMin : alphaMin;
	}

	private double getAlphaMax(double t)
	{
		return (upperFI != null) ? upperFI.getAlpha(t) : alphaMax;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.BasePoissonFisherInformation#postClone()
	 */
	@Override
	protected void postClone()
	{
		// Ensure the function instance is cloned
		if (upperFI != null)
			upperFI = upperFI.clone();
	}
}