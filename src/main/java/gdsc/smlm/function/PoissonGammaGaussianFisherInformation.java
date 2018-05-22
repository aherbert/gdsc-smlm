package gdsc.smlm.function;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optim.AbstractConvergenceChecker;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.util.FastMath;

import gdsc.core.math.NumberUtils;
import gdsc.core.utils.Maths;
import gdsc.core.utils.NotImplementedException;
import gdsc.smlm.utils.Convolution;
import gdsc.smlm.utils.GaussianKernel;
import gnu.trove.list.array.TDoubleArrayList;

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
 * Calculate the Fisher information for a Poisson-Gamma-Gaussian distribution.
 * <p>
 * Uses a modified form of the equation of Chao, et al (2013) Nature Methods, 10, 335-338, SI Eq S8.
 * <p>
 * Performs a convolution with a finite Gaussian kernel.
 */
public class PoissonGammaGaussianFisherInformation extends BasePoissonFisherInformation
{
	/** The default minimum range for the Gaussian kernel (in units of SD). */
	public static final int DEFAULT_MIN_RANGE = 6;

	/** The default maximum range for the Gaussian kernel (in units of SD). */
	public static final int DEFAULT_MAX_RANGE = 38;

	/** The maximum range for the Gaussian kernel (in units of SD). */
	public static final int MAX_RANGE = 38;

	/** The log2 of the maximum scale for the Gaussian kernel (in units of SD). */
	public static final int LOG_2_MAX_SCALE = 16;

	/** The maximum scale for the Gaussian kernel (in units of SD). */
	public static final int MAX_SCALE = 1 << LOG_2_MAX_SCALE;

	/** The default threshold for the relative probability. */
	public static final double DEFAULT_RELATIVE_PROBABILITY_THRESHOLD = 1e-5;

	/** The default threshold for the cumulative probability. */
	public static final double DEFAULT_CUMULATIVE_PROBABILITY = 1 - 1e-6;

	/**
	 * The default sampling of the Gaussian kernel. The kernel will be sampled at s/sampling,
	 * i.e. this is the number of samples to take per standard deviation unit.
	 */
	public static final int DEFAULT_SAMPLING = 4;

	/** The default relative accuracy for convergence. */
	public static final double DEFAULT_RELATIVE_ACCURACY = 1e-6;

	/**
	 * The default number of maximum iterations.
	 * The Gaussian sampling will be doubled during each iteration.
	 */
	public static final int DEFAULT_MAX_ITERATIONS = 4;

	/**
	 * The lowest value for the mean that can be computed. This is the lowest value where the reciprocal is not infinity
	 */
	public static final double MIN_MEAN = Double.longBitsToDouble(0x4000000000001L);

	/** The gain multiplication factor. */
	public final double m;

	/** The standard deviation of the Gaussian. */
	public final double s;

	/** The minimum range of the Gaussian kernel (in SD units). */
	private int minRange = DEFAULT_MIN_RANGE;

	/** The maximum range of the Gaussian kernel (in SD units). */
	private int maxRange = DEFAULT_MAX_RANGE;

	/** The default scale of the kernel. */
	private final int defaultScale;

	/** Working space to store the probabilities. */
	private TDoubleArrayList listP = new TDoubleArrayList();
	/** Working space to store the gradient. */
	private TDoubleArrayList listA = new TDoubleArrayList();

	/**
	 * The relative probability threshold of the Poisson-Gamma distribution that is used. Any probability less than this
	 * is ignored. This prevents using values that do not contribute to the sum.
	 */
	private double relativeProbabilityThreshold = DEFAULT_RELATIVE_PROBABILITY_THRESHOLD;

	/** The upper mean threshold for the switch to half the Poisson Fisher information. */
	private double upperMeanThreshold = 200;

	/** The cumulative probability of the partial gradient of the Poisson-Gamma distribution that is used. */
	private double cumulativeProbability = DEFAULT_CUMULATIVE_PROBABILITY;

	/** Set to true to use Simpson's 3/8 rule for cubic interpolation of the integral. */
	private boolean use38 = true;

	/** Flag to indicate no possible convolution with the gaussian. */
	private final boolean noGaussian;

	/** The gaussian kernel. */
	private GaussianKernel gaussianKernel;

	/** The relative accuracy for convergence. */
	private double relativeAccuracy = DEFAULT_RELATIVE_ACCURACY;

	/** The max iterations. */
	private int maxIterations = DEFAULT_MAX_ITERATIONS;

	/** The step size (h) between values of the last integrated function. */
	private double h;

	/** The probability component of the last integrated function. */
	private double[] P;

	/** The gradient component of the last integrated function. */
	private double[] A;

	/** The start offset of the last integrated function. */
	private double offset;

	/** The Poisson mean of the last integrated function. */
	private double lastT;

	/**
	 * Instantiates a new poisson gamma gaussian fisher information.
	 *
	 * @param m
	 *            the gain multiplication factor
	 * @param s
	 *            the standard deviation of the Gaussian
	 * @throws IllegalArgumentException
	 *             If the standard deviation is not strictly positive
	 */
	public PoissonGammaGaussianFisherInformation(double m, double s) throws IllegalArgumentException
	{
		this(m, s, DEFAULT_SAMPLING);
	}

	/**
	 * Instantiates a new poisson gamma gaussian fisher information.
	 *
	 * @param m
	 *            the gain multiplication factor
	 * @param s
	 *            the standard deviation of the Gaussian
	 * @param sampling
	 *            The number of Gaussian samples to take per standard deviation
	 * @throws IllegalArgumentException
	 *             If the gain is not strictly positive
	 * @throws IllegalArgumentException
	 *             If the standard deviation is below 1
	 * @throws IllegalArgumentException
	 *             If the sampling is below 1
	 * @throws IllegalArgumentException
	 *             If the maximum kernel size after scaling is too large
	 */
	public PoissonGammaGaussianFisherInformation(double m, double s, double sampling) throws IllegalArgumentException
	{
		if (!(m > 0 && m <= Double.MAX_VALUE))
			throw new IllegalArgumentException("Gain multiplication factor must be strictly positive");
		if (!(s >= 0 && s <= Double.MAX_VALUE))
			throw new IllegalArgumentException("Gaussian standard deviation must be positive");
		if (!(sampling >= 1 && sampling <= Double.MAX_VALUE))
			throw new IllegalArgumentException("Gaussian sampling must at least 1");

		this.m = m;
		this.s = s;
		noGaussian = (s * MAX_RANGE < 1);
		if (noGaussian)
		{
			// This is OK. Just return the information for the Poisson-Gamma.
			defaultScale = 0;
		}
		else
		{
			// This is set to work for reasonable values of the Gaussian kernel and sampling
			// e.g. s [0.5:20], sampling from [1:8].

			defaultScale = getPow2Scale(sampling / s);

			// Don't support excess scaling caused by small kernels
			if (defaultScale * s * MAX_RANGE > 1000000000)
			{
				throw new IllegalArgumentException(
						"Maximum Gaussian kernel size too large: " + defaultScale * s * MAX_RANGE);
			}

			gaussianKernel = new GaussianKernel(s);
		}
	}

	/**
	 * Gets a value using the next integer power of 2. This is limited to a size of 2^16.
	 *
	 * @param s
	 *            the value
	 * @return the next power of 2
	 */
	protected static int getPow2Scale(double s)
	{
		double scale = Math.ceil(s);
		if (scale > MAX_SCALE)
			return MAX_SCALE;
		return Maths.nextPow2((int) scale);
	}

	/**
	 * {@inheritDoc}
	 * <p>
	 * The input parameter refers to the mean of the Poisson distribution.
	 * <p>
	 * The Fisher information is computed using the equation of Chao, et al (2013) Nature Methods, 10, 335-338, SI Eq
	 * S8. Note that that equation computes the noise coefficient relative to a Poisson, this computes the Fisher
	 * information. To get the noise coefficient multiply by the input parameter.
	 * 
	 * @see gdsc.smlm.function.FisherInformation#getFisherInformation(double)
	 */
	public double getFisherInformation(double t) throws IllegalArgumentException
	{
		final double I = getPoissonGammaGaussianI(t);

		// Check limits.

		// It should be worse than the Poisson Fisher information (upper limit).
		// Note a low Fisher information is worse as this is the amount of information
		// carried about the parameter.
		final double upper = 1.0 / t; // PoissonFisherInformation.getPoissonI(t);
		// When the mean is high then the lower bound should be half the 
		// Poisson Fisher information. Otherwise check against 0.
		final double lower = (t > 1) ? 0.5 * upper : 0;
		return Maths.clip(lower, upper, I);
	}

	@Override
	public double getAlpha(double t)
	{
		// Don't check against MIN_MEAN as this is done later if it is necessary
		// to compute the Fisher information.

		if (t > upperMeanThreshold)
		{
			// Use an approximation as half the Poisson Fisher information
			return 0.5;
		}

		return t * getFisherInformation(t);
	}

	/**
	 * Gets the Poisson-Gamma-Gaussian Fisher information.
	 * <p>
	 * The input parameter refers to the mean of the Poisson distribution.
	 * <p>
	 * The Fisher information is computed using the equation of Chao, et al (2013) Nature Methods, 10, 335-338, SI Eq
	 * S8. Note that that equation computes the noise coefficient relative to a Poisson, this computes the Fisher
	 * information. To get the noise coefficient multiply by the input parameter.
	 * <p>
	 * Note: This uses a convolution of an infinite integral over a finite range. It may under-estimate the information
	 * when the mean is large.
	 *
	 * @param t
	 *            the Poisson mean
	 * @return the Poisson Gaussian Fisher information
	 * @throws IllegalArgumentException
	 *             the illegal argument exception
	 */
	public double getPoissonGammaGaussianI(double t) throws IllegalArgumentException
	{
		if (t < MIN_MEAN)
		{
			throw new IllegalArgumentException("Poisson mean must be positive");
		}

		lastT = t;

		if (t > upperMeanThreshold)
		{
			// Use an approximation as half the Poisson Fisher information when the mean is large
			return 1.0 / (2 * t);
		}

		final double dirac = PoissonGammaFunction.dirac(t);

		// Special case where the start of the range will have no probability.
		// This occurs when mean > limit of exp(-p) => 746.
		if (dirac == 0)
		{
			// Use half the Poisson fisher information
			return 1.0 / (2 * t);
		}

		// This computes the convolution of a Poisson-Gamma PDF and a Gaussian PDF.
		// The value of this is p(z).

		// The Poisson-Gamma-Gaussian must be differentiated to compute the Fisher information:
		// Expected [ (d ln(p(z)) dv)^2 ]
		// = integral [ (1/p(z) .  d p(z) dv)^2 p(z) dz ]
		// = integral [  1/p(z) . (d p(z) dv)^2 dz ]

		// Gaussian standard deviation = s

		// Note: (fg)' => fg' + f'g
		// e^-v => -e^-v
		// (e^-v * g(v))' => e^-v * g'(v) - e^-v * g(v)
		// So for a probability distribution with a e^-v component:
		// p(z|v)  = e^-v * g(v)
		// p'(z|v) = e^-v * g'(v) - p(z|v)

		// -=-=-=-

		// This Fisher information is based on Chao, et al (2013) Nature Methods, 10, 335-338, SI.

		// Chao et al, Eq 4:
		// (Poisson Binomial Gaussian convolution)

		// Gradient of the Poisson Binomial is:

		// This component can be seen in Chao et al, Eq S8:

		// Substitute the Poisson Binomial Gaussian with Poisson Gamma Gaussian:

		// Gradient of the Poisson Gamma is:

		// Fisher information is:

		// -=-=-=-

		// We need the convolution of:
		// - the Poisson-Gamma with the Gaussian (P)
		// - the Poisson-Gamma partial gradient component with the Gaussian (A)
		// This is combined to a function A^2/P which is integrated.
		// Note that A & P should sum to 1 over the full range used for integration.
		// In practice a relative height threshold of the function can be used
		// to determine the range where it is significant.

		// Note that when the Poisson mean is low then the contribution from the 
		// Dirac delta function at 0 will be large. As the mean increases 
		// the Dirac will have little significance, the function is very smooth
		// and may have a large range. Thus convolution with an extended Gaussian kernel 
		// may not be necessary over the full range.
		// So the integral can be computed as [integral range1] + [integral range2]

		// The function we are integrating is:  
		// E = integral [ (1/p(z) . d p(z) dv)^2 p(z) dz ]
		// Integration is done using all values of this function within a reasonable
		// relative size to the maximum.

		// Determine the extent of the function: A^2/P
		double[] maximum = findMaximum(t, 1e-4);
		double[] upper = findUpperLimit(t, maximum, relativeProbabilityThreshold);
		double upperLimit = upper[0];

		// Note: This method works by using the Poisson-Gamma PDF as a substitute
		// for a Poisson-Binomial PMF. For simplicity the PDF is integrated over 
		// integer intervals to create a discrete PMF. This is a fair approximation
		// of how an EM-CCD camera would work as input electrons should produce
		// discrete output electrons.

		if (noGaussian)
		{
			// Special case. Compute E [ A^2/P ] using only the Poisson-Gamma PMF. 
			computePMF(dirac, t, (int) Math.ceil(upperLimit));

			// No interpolated integration necessary as this is discrete
			double sum = 0;
			double[] p = listP.toArray();
			double[] a = listA.toArray();
			for (int i = 0; i < p.length; i++)
			{
				sum += getF(p[i], a[i]);
			}
			return sum;
		}

		// The exponent provides a rough idea of the size of the mean (i.e. log2(mean))
		int exp = NumberUtils.getSignedExponent(t);
		// As the mean reduces the Poisson distribution is more skewed 
		// and the extent of the kernel must change. Just increase the range
		// for the kernel for each power of 2 the number is below 1.
		int range1 = minRange;
		/// XXX - check this is needed if iteration is working.
		for (int e = exp; range1 < maxRange && e <= 0; e++, range1++)
			;

		// Extend the upper limit using the configured kernel range
		upperLimit += s * range1;

		// TODO:
		// Set range1 using the range of the Gaussian. 
		// This is the range around c=0 where the Dirac has a large effect.
		// Set range2 using the rest of the range.
		// Iteratively integrate range 1 up to a maximum kernel size.
		// Iteratively Integrate range 2.
		// Combine the sums.

		// Range 1:
		// The range around the Dirac delta function is computed using a full convolution
		// with the Gaussian. This is iterated until convergence using more Gaussian samples.

		// Range 2: 
		// The remaining range of the Poisson-Gamma is integrated using an integrator until the
		// the Gaussian kernel is non-trivial.

		// For integration to work efficiently the Gaussian scaling should be in powers of 2. 
		// Adjust the ranges appropriately so that the ranges start and end on powers of 1.
		int iEndRange1 = Maths.nextPow2((int) Math.ceil(s * range1));
		int iEndRange2 = iEndRange1;
		if (upperLimit > iEndRange1)
			iEndRange2 += Maths.nextPow2((int) Math.ceil(upperLimit - iEndRange1));

		// Compute the Poisson-Gamma over the range. This is a PDF but for 
		// simplicity it is integrated to a discrete PMF.
		computePMF(dirac, t, iEndRange2);

		int scale = defaultScale;

		// Testing - do a single convolution

		// Convolve with the Gaussian kernel
		double[] p = listP.toArray();
		double[] a = listA.toArray();
		double[] g = gaussianKernel.getGaussianKernel(scale, range1, false);
		h = 1.0 / scale;

		convolve(p, a, g, true, scale);

		System.out.printf("t=%g 0=%d scale=%d kernel=%d After convolution SumP=%s SumA=%s\n", t, 
				iEndRange2, scale, g.length,
				Maths.sum(P), Maths.sum(A));

		int maxi = P.length;

		// Note: The step size for Simpson integration should be 1/scale.
		// Because the kernel is normalised each doubling in scale above 1
		// results in the Gaussian kernel at each point halving. 
		// The result is that the P and A values halve. A^2/P will then halve. 
		// So it should be increased by scale. 
		// This makes h effectively 1.
		
		// We assume that the function values at the end are zero and so do not 
		// include them in the sum. Just alternate totals.
		double range1sum;
		double sum = 0;
		if (use38)
		{
			// Simpson's 3/8 rule based on cubic interpolation has a lower error.
			// This computes the sum as:
			// 3h/8 * [ f(x0) + 3f(x1) + 3f(x2) + 2f(x3) + 3f(x4) + 3f(x5) + 2f(x6) + ... + f(xn) ]
			// The final step before maxi must be sum 2. 
			final int end = (maxi - 1) % 3;
			double sum3 = 0, sum2 = 0;
			for (int i = 0; i < maxi; i++)
			{
				final double f = getF(P[i], A[i]);
				sum += f;
				if (i % 3 == end)
					sum2 += f;
				else
					sum3 += f;
			}

			// Assume all values before the function are zero.
			range1sum = (
			//h * 
			3.0 / 8) * (sum3 * 3 + sum2 * 2);
		}
		else
		{
			// Simpson's rule.
			// This computes the sum as:
			// h/3 * [ f(x0) + 4f(x1) + 2f(x2) + 4f(x3) + 2f(x4) ... + 4f(xn-1) + f(xn) ]
			// The number of steps must be modulo 2. Assume all values before this are zero.
			// The final step before maxi must be sum 4. 
			final int end = (maxi - 1) % 2;
			double sum4 = 0, sum2 = 0;
			for (int i = 0; i < maxi; i++)
			{
				final double f = getF(P[i], A[i]);
				sum += f;
				if (i % 2 == end)
					sum4 += f;
				else
					sum2 += f;
			}

			// Assume all values before the function are zero.
			range1sum = (
			//h * 
			1.0 / 3) * (sum4 * 4 + sum2 * 2);
		}

		return range1sum - 1;
	}

	private void computePMF(double dirac, double t, int endX)
	{
		double G;
		double[] dG_dp = new double[1];

		listP.resetQuick();
		listA.resetQuick();

		// Initialise at x=0
		G = PoissonGammaFunction.poissonGammaPartial(0, t, m, dG_dp);

		// At x=0 integrate from 0 to 0.5 without the Dirac 
		// (this is a smooth function)
		G -= dirac;
		G = add(t, G, dG_dp, 0, 0.25, 0.5);
		listP.setQuick(0, listP.getQuick(0) + dirac);

		// For the rest of the range integrate from x-0.5 to x+0.5
		for (int x = 1; x <= endX; x++)
		{
			G = add(t, G, dG_dp, x - 0.5, x, x + 0.5);
			if (G == 0)
			{
				listP.removeAt(x);
				listA.removeAt(x);
				break;
			}
		}
		if (listP.size() != endX)
		{
			// Zero fill for the convolution?
		}
		System.out.printf("t=%g 0=%d Dirac=%s SumP=%s SumA=%s\n", t, endX, dirac, listP.sum(), listA.sum());
	}

	/**
	 * Adds the integral of the Poisson-Gamma PDF to the working lists for P and A. It is assumed that the input is the
	 * value and gradient of the lower bound.
	 * When finished return the value and gradient of the upper bound.
	 *
	 * @param G
	 *            the previous computed point value
	 * @param dG_dp
	 *            the previous computed point gradient
	 * @param lower
	 *            the lower bound
	 * @param mid
	 *            the mid point
	 * @param upper
	 *            the upper bound
	 * @return the double
	 */
	private double add(double t, double G, double[] dG_dp, double lower, double mid, double upper)
	{
		// TODO - does this need refinement?

		// Do a single Simpson sum (f0 + 4 * fxi + fn).
		double sG = G;
		double sA = dG_dp[0];
		G = PoissonGammaFunction.poissonGammaPartial(mid, t, m, dG_dp);
		sG += 4 * G;
		sA += 4 * dG_dp[0];
		G = PoissonGammaFunction.poissonGammaPartial(upper, t, m, dG_dp);
		double h_3 = (upper - lower) / 6;
		listP.add((sG + G) * h_3);
		listA.add((sA + dG_dp[0]) * h_3);
		return G;
	}

	private void convolve(double[] p, double[] a, double[] g, boolean requireEdge, int scale)
	{
		// Convolution in the frequency domain may create negatives
		// due to edge wrap artifacts. This is classically reduced by using a 
		// edge window function. However this is undesirable as it 
		// will truncate the dirac delta function and high contribution
		// at c=0 when the Poisson mean is low.

		// When the entire range of the convolution is required use the spatial domain.
		// When only part of the range is required (since the range start has been 
		// computed, and the range end will asymptote to zero) use the FFT

//		double[][] result = (requireEdge) ? Convolution.convolve(g, p, a) : Convolution.convolveFast(g, p, a);
//		P = result[0];
//		A = result[1];
		
		// XXX - fix for a twin scaled convolution
		P = Convolution.convolve(g, p, scale);
		A = Convolution.convolve(g, a, scale);
		offset = g.length / 2 * -h;

//		if (!requireEdge && Convolution.isFFT(g.length, p.length))
//		{
//			// Just remove values below the minimum expected.
//			// If the Gaussian kernel is sufficiently large then 
//			// edge artifacts should be limited to the region where
//			// the actual values are very small.
//			double minP = Maths.min(p) * g[0];
//			double minA = Maths.min(a) * g[0];
//			for (int i = 0; i < P.length; i++)
//			{
//				if (P[i] < minP)
//					P[i] = 0;
//				if (A[i] < minA)
//					A[i] = 0;
//			}
//		}
	}

	private static double getF(double P, double A)
	{
		//		if (NumberUtils.isSubNormal(P))
		//			System.out.printf("Sub-normal: P=%s A=%s\n", P, A);
		if (P == 0)
			return 0;
		final double result = Maths.pow2(A) / P;
		//if (result > 1)
		//	System.out.printf("Strange: P=%s A=%s result=%s\n", P, A, result);
		return result;
	}

	/**
	 * Gets the fisher information function. This is the function that is integrated to get the Fisher information. This
	 * can be called after a called to {@link #getFisherInformation(double)}. It will return the last value of the
	 * function that was computed (before/after convolution with the Gaussian kernel).
	 * <p>
	 * The function should be smooth with zero at the upper end. The lower end may not be smooth due to the Dirac delta
	 * function at c=0.
	 *
	 * @param withGaussian
	 *            Set to true to return the function after convolution with the gaussian
	 * @return the fisher information function
	 */
	public double[][] getFisherInformationFunction(boolean withGaussian)
	{
		if (listP.isEmpty())
			return null;

		double[] A, P;
		double offset;
		if (withGaussian && !noGaussian)
		{
			P = this.P;
			A = this.A;
			offset = this.offset;
		}
		else
		{
			P = listP.toArray();
			A = listA.toArray();
			offset = 0;
		}
		double[] f = new double[A.length];
		double[] x = new double[A.length];
		for (int i = 0; i < A.length; i++)
		{
			x[i] = i * h + offset;
			f[i] = getF(P[i], A[i]);
		}
		return new double[][] { x, f };
	}

	/**
	 * Gets the min range of the Gaussian kernel in SD units.
	 * The Gaussian kernel is scaled dynamically depending on the shape of the Poisson distribution.
	 *
	 * @return the min range
	 */
	public int getMinRange()
	{
		return minRange;
	}

	/**
	 * Sets the min range of the Gaussian kernel in SD units.
	 * The Gaussian kernel is scaled dynamically depending on the shape of the Poisson distribution.
	 *
	 * @param minRange
	 *            the new min range
	 */
	public void setMinRange(int minRange)
	{
		this.minRange = checkRange(minRange);
	}

	private int checkRange(int range)
	{
		// Gaussian = Math.exp(-0.5 * x^2)
		// FastMath.exp(-746) == 0
		// => range for the Gaussian is sqrt(2*746) = 38.6
		return Maths.clip(1, MAX_RANGE, range);
	}

	/**
	 * Gets the max range of the Gaussian kernel in SD units.
	 * The Gaussian kernel is scaled dynamically depending on the shape of the Poisson distribution.
	 *
	 * @return the max range
	 */
	public int getMaxRange()
	{
		return maxRange;
	}

	/**
	 * Sets the max range of the Gaussian kernel in SD units.
	 * The Gaussian kernel is scaled dynamically depending on the shape of the Poisson distribution.
	 *
	 * @param maxRange
	 *            the new max range
	 */
	public void setMaxRange(int maxRange)
	{
		this.maxRange = checkRange(maxRange);
	}

	/**
	 * Gets the mean threshold for the switch to half the Poisson Fisher information.
	 * <p>
	 * When the mean is high then there is no Dirac delta contribution to the probability function. This occurs when
	 * exp(-mean) is zero, which is approximately 746. In practice there is not much difference when the mean is above
	 * 100.
	 *
	 * @return the mean threshold
	 */
	public double getMeanThreshold()
	{
		return upperMeanThreshold;
	}

	/**
	 * Sets the mean threshold for the switch to half the Poisson Fisher information.
	 * <p>
	 * When the mean is high then there is no Dirac delta contribution to the probability function. This occurs when
	 * exp(-mean) is zero, which is approximately 746. In practice there is not much difference when the mean is above
	 * 100.
	 *
	 * @param meanThreshold
	 *            the new mean threshold
	 */
	public void setMeanThreshold(double meanThreshold)
	{
		this.upperMeanThreshold = meanThreshold;
	}

	/**
	 * Gets the cumulative probability of the partial gradient of the Poisson-Gamma distribution that is used.
	 *
	 * @return the cumulative probability
	 */
	public double getCumulativeProbability()
	{
		return cumulativeProbability;
	}

	/**
	 * Sets the cumulative probability of the partial gradient of the Poisson-Gamma distribution that is used.
	 *
	 * @param cumulativeProbability
	 *            the new cumulative probability
	 */
	public void setCumulativeProbability(double cumulativeProbability)
	{
		if (!(cumulativeProbability > 0 && cumulativeProbability <= 1))
			throw new IllegalArgumentException("P must be in the range 0-1");
		this.cumulativeProbability = cumulativeProbability;
	}

	/**
	 * Gets the relative probability threshold of the Poisson-Gamma distribution that
	 * is used. Any probability less than this is ignored. This prevents using values
	 * that do not contribute to the sum.
	 *
	 * @return the relative probability threshold
	 */
	public double getRelativeProbabilityThreshold()
	{
		return relativeProbabilityThreshold;
	}

	/**
	 * Sets the relative probability threshold of the Poisson-Gamma distribution that
	 * is used. Any probability less than this is ignored. This prevents using values
	 * that do not contribute to the sum.
	 * <p>
	 * This must be in the range 0 to 0.5. Larger values would prevent using any
	 * slope of the Poisson-Gamma distribution. A value of zero requires all probability values to be used and may
	 * result in long computation times.
	 *
	 * @param relativeProbabilityThreshold
	 *            the new relative probability threshold
	 */
	public void setRelativeProbabilityThreshold(double relativeProbabilityThreshold)
	{
		if (!(relativeProbabilityThreshold >= 0 && relativeProbabilityThreshold <= 0.5))
			throw new IllegalArgumentException("P must be in the range 0-0.5");
		this.relativeProbabilityThreshold = relativeProbabilityThreshold;
	}

	/**
	 * If true, use Simpson's 3/8 rule for cubic interpolation of the integral. False uses Simpson's rule for
	 * quadratic interpolation.
	 *
	 * @return the use 38
	 */
	public boolean getUse38()
	{
		return use38;
	}

	/**
	 * Set to true to use Simpson's 3/8 rule for cubic interpolation of the integral. False uses Simpson's rule for
	 * quadratic interpolation.
	 *
	 * @param use38
	 *            the new use 38
	 */
	public void setUse38(boolean use38)
	{
		this.use38 = use38;
	}

	/**
	 * Gets the relative accuracy for convergence during iteration.
	 *
	 * @return the relative accuracy
	 */
	public double getRelativeAccuracy()
	{
		return relativeAccuracy;
	}

	/**
	 * Sets the relative accuracy for convergence during iteration.
	 * <p>
	 * Set to below zero to prevent convergence check. This results in a fixed number of iterations.
	 *
	 * @param relativeAccuracy
	 *            the new relative accuracy
	 */
	public void setRelativeAccuracy(double relativeAccuracy)
	{
		this.relativeAccuracy = relativeAccuracy;
	}

	/**
	 * Gets the max iterations for iteration.
	 *
	 * @return the max iterations
	 */
	public int getMaxIterations()
	{
		return maxIterations;
	}

	/**
	 * Sets the max iterations for iteration. Set to zero to prevent iteration.
	 *
	 * @param maxIterations
	 *            the new max iterations
	 */
	public void setMaxIterations(int maxIterations)
	{
		this.maxIterations = maxIterations;
	}

	@Override
	protected void postClone()
	{
		listP = new TDoubleArrayList();
		listA = new TDoubleArrayList();
		if (gaussianKernel != null)
			gaussianKernel = gaussianKernel.clone();
	}

	private static class QuickConvergenceChecker extends AbstractConvergenceChecker<UnivariatePointValuePair>
	{
		int fixedIterations;
		UnivariatePointValuePair previous, current;

		public QuickConvergenceChecker(int fixedIterations)
		{
			super(0, 0);
			this.fixedIterations = fixedIterations;
		}

		@Override
		public boolean converged(int iteration, UnivariatePointValuePair previous, UnivariatePointValuePair current)
		{
			this.previous = previous;
			this.current = current;
			return iteration > fixedIterations;
		}

		UnivariatePointValuePair getBest()
		{
			if (previous == null || previous.getValue() < current.getValue())
				return current;
			return previous;
		}
	}

	/**
	 * Find the maximum of the integrated function A^2/P. P is the Poisson-Gamma convolution, A is the partial gradient.
	 * <p>
	 * When both A and P are convolved with a Gaussian kernel, the integral of this function - 1
	 * is the Fisher information.
	 * <p>
	 * This method is used to determine the maximum of the function.
	 * The integral of A should equal 1. The range can be determined using a fraction of the maximum, or integrating
	 * until the sum is 1.
	 * 
	 * @param t
	 *            the Poisson mean
	 * @param rel
	 *            Relative threshold
	 * @return [max,max value,evaluations,upper,upper value,evaluations]
	 */
	public double[] findMaximum(final double t, double rel)
	{
		if (t < MIN_MEAN)
		{
			throw new IllegalArgumentException("Poisson mean must be positive");
		}

		final double dirac = PoissonGammaFunction.dirac(t);

		final UnivariateFunction f = new UnivariateFunction()
		{
			double[] dG_dp = new double[1];

			public double value(double x)
			{
				// Return F without the dirac. This is so that the maximum of the function
				// is known for relative height analysis.
				double G;
				if (x == 0)
				{
					dG_dp[0] = dirac / m;
					G = dG_dp[0] * t;

					//System.out.printf("m=%g s=%g u=%g max=%s %s\n", m, s, t, 0, getF(G, dG_dp[0]));
					//G = PoissonGammaFunction.poissonGammaPartial(x, t, m, dG_dp);
					//g-=dirac;
				}
				else
					G = PoissonGammaFunction.poissonGammaPartial(x, t, m, dG_dp);
				return getF(G, dG_dp[0]);
			}
		};

		// Determine the extent of the function: A^2/P
		// Without convolution this will have a maximum that is above, and 
		// converges to [Poisson mean * amplification].
		double mean = t * m;

		int iterations = 10;
		QuickConvergenceChecker checker = new QuickConvergenceChecker(iterations);
		BrentOptimizer opt = new BrentOptimizer(rel, Double.MIN_VALUE, checker);
		UnivariatePointValuePair pair;
		try
		{
			double start = (t < 1) ? 0 : mean;
			pair = opt.optimize(new SearchInterval(0, mean * 1.5, start), GoalType.MAXIMIZE, new MaxEval(50000),
					new UnivariateObjectiveFunction(f));
		}
		catch (TooManyEvaluationsException e)
		{
			// This should not occur with the quick checker fixed iterations
			// - just get the current best
			pair = checker.getBest();
		}

		double[] max = new double[3];
		max[0] = pair.getPoint();
		max[1] = pair.getValue();
		max[2] = opt.getEvaluations();
		return max;
	}

	/**
	 * Minimum relative tolerance.
	 */
	private static final double MIN_RELATIVE_TOLERANCE = 2 * FastMath.ulp(1d);

	/**
	 * Find the upper limit of the integrated function A^2/P. P is the Poisson-Gamma convolution, A is the partial
	 * gradient.
	 * <p>
	 * When both A and P are convolved with a Gaussian kernel, the integral of this function - 1
	 * is the Fisher information.
	 * <p>
	 * This method is used to determine the upper limit of the function using a binary search.
	 *
	 * @param t
	 *            the Poisson mean
	 * @param max
	 *            the max of the function (returned from {@link #findMaximum(double, double)})
	 * @param rel
	 *            Relative threshold
	 * @return [upper,upper value,evaluations]
	 */
	public double[] findUpperLimit(final double t, double[] max, double rel)
	{
		if (rel < MIN_RELATIVE_TOLERANCE)
		{
			throw new IllegalArgumentException("Relative tolerance too small: " + rel);
		}

		final UnivariateFunction f = new UnivariateFunction()
		{
			double[] dG_dp = new double[1];

			public double value(double x)
			{
				double G = PoissonGammaFunction.poissonGammaPartial(x, t, m, dG_dp);
				return getF(G, dG_dp[0]);
			}
		};

		int eval = 0;

		// Increase from the max until the tolerance is achieved.

		// Use the mean to get a rough initial step size
		double mean = t * m;
		double step = Math.max(mean, 1);

		double upper = max[0];
		double g = max[1];
		double threshold = g * rel;

		double lower = upper;
		while (g > threshold)
		{
			lower = upper;
			upper += step;
			step *= 2;
			eval++;
			g = f.value(upper);
		}

		// Binary search the bracket between lower and upper
		while (lower + 1 < upper)
		{
			double mid = (lower + upper) * 0.5;
			eval++;
			double midg = f.value(mid);
			if (midg > threshold)
			{
				lower = mid;
			}
			else
			{
				upper = mid;
				g = midg;
			}
		}

		return new double[] { upper, g, eval };
	}
}