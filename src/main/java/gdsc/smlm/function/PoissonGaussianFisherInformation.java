package gdsc.smlm.function;

import java.util.Arrays;

import org.apache.commons.math3.distribution.CustomPoissonDistribution;

import gdsc.core.utils.Maths;
import gdsc.smlm.utils.Convolution;
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
 * Calculate the Fisher information for a Poisson-Gaussian distribution.
 * <p>
 * Uses the equation of Chao, et al (2013) Nature Methods, 10, 335-338, SI Eq S7.
 * <p>
 * Performs a convolution with a finite Gaussian kernel. The Gaussian is constructed using a range of the standard
 * deviation (s) and sampled at least every s/2.
 * <p>
 * An optimisation is used to avoid computation on tiny Gaussian kernels (i.e. too small to be computed)
 * This will occur when the Gaussian standard deviation is less than 0.02. The result is no convolution and the result
 * computes Poisson Fisher information.
 * <p>
 * An optimisation is used when the mean of the Poisson is above a threshold. In this case the Poisson can be
 * approximated as a Gaussian and the Fisher information is returned for the Gaussian-Gaussian convolution.
 */
public class PoissonGaussianFisherInformation implements FisherInformation
{
	public static final double DEFAULT_CUMULATIVE_PROBABILITY = 1 - 1e-10;

	/** Store the limit of the Poisson distribution for small mean for the default cumulative probability. */
	private static final int[] defaultLimits;
	/** Store the limit of the Poisson distribution for tiny mean for the default cumulative probability. */
	private static final int[] defaultTinyLimits;

	static
	{
		CustomPoissonDistribution pd = new CustomPoissonDistribution(null, 1);
		defaultLimits = new int[101];
		for (int i = 1; i < defaultLimits.length; i++)
		{
			defaultLimits[i] = computeLimit(pd, i, DEFAULT_CUMULATIVE_PROBABILITY);
			//System.out.printf("[%d] = %d  scale=%d\n", i, defaultLimits[i], getScale(Math.sqrt(i)));
		}

		// Use exponent of the mean down to -20
		defaultTinyLimits = new int[21];
		for (int i = 1; i < defaultTinyLimits.length; i++)
		{
			defaultTinyLimits[i] = computeTinyLimit(pd, -i, DEFAULT_CUMULATIVE_PROBABILITY);
			//System.out.printf("[%d] = %d : p(0) = %g : cumul = %s : next = %g\n", i, defaultTinyLimits[i],
			//		pd.probability(0), pd.cumulativeProbability(defaultTinyLimits[i]),
			//		pd.probability(defaultTinyLimits[i] + 1));
		}
	}

	/**
	 * Compute the limit of a usable probability above 0.
	 * <p>
	 * Find the point where probability will return above 0.
	 * Using FastMath.exp this is -746. However sub-normal output occurs at -709.
	 * This is a good limit for computation.
	 *
	 * @param mean
	 *            the mean
	 * @return the limit
	 */
	public static int computeLimit(double mean)
	{
		CustomPoissonDistribution pd = new CustomPoissonDistribution(null, mean);
		int x = (int) mean;
		while (pd.logProbability(x + 1) > -709)
			x++;
		return x;
	}

	/**
	 * Compute the limit of a usable probability above 0.
	 *
	 * @param pd
	 *            the pd
	 * @param mean
	 *            the mean
	 * @param cumulativeProbability
	 *            the cumulative probability
	 * @return the limit
	 */
	private static int computeLimit(CustomPoissonDistribution pd, double mean, double cumulativeProbability)
	{
		pd.setMeanUnsafe(mean);
		return pd.inverseCumulativeProbability(cumulativeProbability);
	}

	/**
	 * Compute the limit of a usable probability above 0.
	 *
	 * @param pd
	 *            the pd
	 * @param mean
	 *            the mean
	 * @param cumulativeProbability
	 *            the cumulative probability
	 * @return the limit
	 */
	private static int computeTinyLimit(CustomPoissonDistribution pd, int exp, double cumulativeProbability)
	{
		// Fill all bits of the mantissa
		long bits = 0xfffffffffffffL;
		double mean = Double.longBitsToDouble(bits | (long) (exp + 1023) << 52);
		pd.setMeanUnsafe(mean);
		return pd.inverseCumulativeProbability(cumulativeProbability);
	}

	/** The standard deviation of the Gaussian. */
	public final double s;

	/** The range of the Gaussian kernel (in SD units). */
	public final double range;

	/** The default scale for the kernel. */
	private final int defaultScale;

	/** The poisson distribution used to generate the Poisson probabilities. */
	private CustomPoissonDistribution pd = new CustomPoissonDistribution(null, 1);

	/**
	 * The Gaussian convolution kernels for different scaling. The scale is 2^index, e.g. 1, 2, 4, 8, 16, 32, 64, 128.
	 */
	private final double[][] kernel;

	/** Working space to store the Poisson probabilities. */
	private TDoubleArrayList list = new TDoubleArrayList();

	/** The mean threshold for the switch to a Gaussian-Gaussian convolution. */
	private double meanThreshold = 100;

	/** The cumulative probability of the Poisson distribution that is used. */
	private double cumulativeProbability = DEFAULT_CUMULATIVE_PROBABILITY;

	/** Store the limit of the Poisson distribution for small mean for the cumulative probability. */
	private int[] limits = defaultLimits;

	/** Store the limit of the Poisson distribution for tiny mean for the cumulative probability. */
	private int[] tinyLimits = defaultTinyLimits;

	/**
	 * Instantiates a new poisson gaussian fisher information.
	 *
	 * @param s
	 *            the standard deviation of the Gaussian
	 * @throws IllegalArgumentException
	 *             If the standard deviation is not strictly positive
	 */
	public PoissonGaussianFisherInformation(double s) throws IllegalArgumentException
	{
		this(s, 5);
	}

	/**
	 * Instantiates a new poisson gaussian fisher information.
	 *
	 * @param s
	 *            the standard deviation of the Gaussian
	 * @param range
	 *            the range of the Gaussian kernel (in SD units). This is clipped to the range
	 *            1-38 to provide a meaningful convolution.
	 * @throws IllegalArgumentException
	 *             If the standard deviation is not strictly positive
	 */
	public PoissonGaussianFisherInformation(double s, double range) throws IllegalArgumentException
	{
		if (!(s > 0 && s <= Double.MAX_VALUE))
			throw new IllegalArgumentException("Gaussian variance must be strictly positive");

		// Gaussian = Math.exp(-0.5 * x^2)
		// FastMath.exp(-746) == 0
		// => range for the Gaussian is sqrt(2*746) = 38.6

		if (Double.isNaN(range))
			throw new IllegalArgumentException("Gaussian range must not be NaN");
		range = Maths.clip(1, 38, range);

		this.s = s;

		// Check if the Gaussian standard deviation is above the threshold for computation.
		// Also check if the gaussian filter will touch more than one Poisson value.
		// Otherwise convolution is not possible.
		// The limit s ==0.02 is based on the scale being 2/s = 100. Do not support scaling 
		// greater than this. It is unlikely anyway.
		if (s >= 0.02 && s * range >= 1)
		{
			this.range = range;

			// Determine how much to up-sample so that the convolution with the Gaussian
			// uses multiple values of the Gaussian.
			defaultScale = getScale(s);

			// Store the Gaussian kernels for convolution:
			// 1, 2, 4, 8, 16, 32, 64, 128
			kernel = new double[8][];
		}
		else
		{
			this.range = 0;
			defaultScale = 0;
			kernel = null;
		}
	}

	private static int getScale(double s)
	{
		double scale = Math.ceil(2 / s);
		if (scale > 128)
			return 128;
		return Maths.nextPow2((int) scale);
	}

	/**
	 * {@inheritDoc}
	 * <p>
	 * The input parameter refers to the mean of the Poisson distribution.
	 * <p>
	 * The Fisher information is computed using the equation of Chao, et al (2013) Nature Methods, 10, 335-338, SI Eq
	 * S7. Note that that equation computes the noise coefficient relative to a Poisson, this computes the Fisher
	 * information. To get the noise coefficient multiply by the input parameter.
	 * 
	 * @see gdsc.smlm.function.FisherInformation#getFisherInformation(double)
	 */
	public double getFisherInformation(double t) throws IllegalArgumentException
	{
		if (t <= 0)
		{
			//throw new IllegalArgumentException("Poisson mean must be positive");

			// No Poisson. Return the Fisher information for a Gaussian
			return getGaussianI();
		}

		if (t > meanThreshold)
			return getPoissonGaussianI(t);

		// Get the Fisher information for a Poisson. 
		final double pI = getPoissonI(t);

		if (kernel == null)
			// No Gaussian convolution
			return pI;

		// This computes the convolution of a Poisson PMF and a Gaussian PDF.
		// The value of this is p(z).

		// The Poisson-Gaussian must be differentiated to compute the Fisher information:
		// Expected [ (d ln(p(z)) dv)^2 ]
		// = integral [ (1/p(z) . d p(z) dv)^2 p(z) dz ]
		// = integral [ (1/p(z) . (d p(z) dv)^2 dz ]

		// Gaussian standard deviation = s

		// Chao et al, S5:
		// p(z) = 1/sqrt(2pi)s sum_j=0:Inf  e^-v . v^j / j! . e^-1/2((z-j)/s)^2

		// This is the sum over j of the probability of Poisson(j) * probability of Gaussian(z-j) 

		// Note: (fg)' => f'g + fg'
		// e^-v => -e^-v
		// v^j => j.v^(j-1)
		// e^-v v^j / j! => e^-v v^(j-1) / (j-1)! - e^-v v^j / j!

		// d p(z) dv = 1/sqrt(2pi)s sum_j=1:Inf  e^-v . v^(j-1) / (j-1)! . e^-1/2((z-j)/s)^2 - 
		//                          sum_j=0:Inf  e^-v . v^j / j! . e^-1/2((z-j)/s)^2
		// Note: j=0 differentiates to -e^v since v^j / j! = 1. This removes j=0 from the first sum
		// but not the second.
		// Over the sum the second term adds up to become p(z) so:
		// d p(z) dv = (1/sqrt(2pi)s sum_j=1:Inf  e^-v . v^(j-1) / (j-1)! . e^-1/2((z-j)/s)^2 ) - p(z) 

		// Set the first term to A, the second to P:
		// d p(z) dv = A - P

		// E = integral [ (1/p(z) . d p(z) dv)^2 p(z) dz ]
		//   = integral [ (1/P    . (A - P))^2 * P ]
		//   = integral [ (1/P^2  . (A^2 - 2AP + p^2) * P ]
		//   = integral [ (A^2/P^2 - 2A/P + 1) * P ]
		//   = integral [  A^2/P   - 2A   + P ]
		//   = integral [A^2/P] - integral [2A] + integral [P]

		// Note that the integral of P==1.
		// Since the integral of A is just P offset by j-1, integral of A==1

		// E = integral [A^2/P] - 1

		// P(z) = 1/sqrt(2pi)s sum_j=0:Inf  e^-v . v^j / j!         . e^-1/2((z-j)/s)^2
		// A(z) = 1/sqrt(2pi)s sum_j=1:Inf  e^-v . v^(j-1) / (j-1)! . e^-1/2((z-j)/s)^2
		// A(z) = 1/sqrt(2pi)s sum_j=0:Inf  e^-v . v^j / j!         . e^-1/2((z-(j+1))/s)^2
		// A(z) = P(z-1)

		// We need the convolution of the Poisson with the Gaussian 

		// The Poisson is a PMF. The Gaussian is a PDF.
		// The expected value is integrated over all real z, -Inf:Inf, 
		// (the full range of the Gaussian).

		// Sample the values of the full range and compute a sum using Simpson integration.

		// Build the Poisson distribution.
		pd.setMeanUnsafe(t);

		// Only use part of the cumulative distribution.
		// Start at zero for a 1-tailed truncation of the cumulative distribution.
		// This code will be for low mean values where the contribution at zero is 
		// significant.
		int minx = 0;
		// Find the limit. These can be cached (or may be the defaults).
		// TODO - Determine if the Poisson can be truncated. We may have to use more of 
		// the values (for example those returned by computeLimit(...).
		int maxx;
		if (t < 1)
		{
			int exp = -FastLog.getSignedExponent(t);
			if (exp >= tinyLimits.length)
				exp = tinyLimits.length - 1;
			if (tinyLimits[exp] == 0)
				tinyLimits[exp] = computeTinyLimit(pd, -exp, cumulativeProbability);
			maxx = tinyLimits[exp];
		}
		else
		{
			int x = (int) Math.ceil(t);
			if (x < limits.length)
			{
				if (limits[x] == 0)
					limits[x] = computeLimit(pd, x, cumulativeProbability);
				maxx = limits[x];
			}
			else
			{
				// For large mean the distribution will be far from zero.
				// In this case use a 2-tailed limit.
				double lower = (1 - cumulativeProbability) / 2;
				minx = computeLimit(pd, x, lower);
				maxx = computeLimit(pd, x, 1 - lower);
			}
		}

		list.resetQuick();
		for (int x = minx; x < maxx; x++)
		{
			double pp = pd.probability(x);
			list.add(pp);
		}
		// Final value may be zero
		{
			double pp = pd.probability(maxx);
			if (pp != 0)
				list.add(pp);
		}

		if (list.size() < 2)
		{
			// Extreme case where there is no Poisson for convolution.
			// Assume a Gaussian distribution. Return the Fisher information
			// for the Gaussian with mean 0. This will happen when the cumulative 
			// probability has been altered from the default.
			return getGaussianI();
		}

		// Unscaled Poisson
		double[] p = list.toArray();

		// Choose the kernel. A small mean requires more Gaussian samples.
		// Note the default scale is the minimum required to sample at 0.5 SD units.
		// Find the same for the Poisson using its variance.
		// mean 4 => scale = 1
		// mean <4 => scale = 2
		// mean <1 => scale >= 4
		// This may have to be changed.
		int scale = Math.max(defaultScale, getScale(t / 2));

		//scale = 1;
		//scale = defaultScale;

		// Get the Gaussian kernel
		int index = log2(scale);
		if (kernel[index] == null)
			kernel[index] = Convolution.makeGaussianKernel(s * scale, range, true);
		double[] g = kernel[index];

		// Up-sample the Poisson
		if (scale != 1)
		{
			list.resetQuick();
			double[] pad = new double[scale - 1];
			list.add(p[0]);
			for (int i = 1; i < p.length; i++)
			{
				list.add(pad);
				list.add(p[i]);
			}
			p = list.toArray();
		}

		// Convolve with the Gaussian kernel
		double[] pg = Convolution.convolveFast(p, g);

		// In order for A(z) = P(z-1) to work sum A(z) must be 1
		double sum = 0;
		for (int i = 0; i < pg.length; i++)
			sum += pg[i];
		//System.out.printf("Normalisation = %s\n", sum);
		for (int i = 0; i < pg.length; i++)
			pg[i] /= sum;

		// Integrate function:
		// E = integral [A^2/P] - 1
		// P(z) = Poisson-Gaussian convolution
		// A(z) = P(z-1)

		// The offset for P(z-1) is the scale. 
		// When P(z-1) does not exist assume it is zero. Therefore only
		// integrate over the range scale:P.length
		int mini = scale;
		int maxi = pg.length;
		int offseti = -scale;
		
		// Compute the sum using Simpson's integration. 
		// h = interval = (b-a)/n
		// The integral range is:
		// f(x0) = 0 where x0 = -1;
		// f(xn) = 0 where xn = length;
		// n = length-1 + 2 = length+1
		// a = -1
		// b = length
		// h = (length - -1) / (length+1) = 1
		double h = 1;

		// We assume that the function values at the end are zero and so do not 
		// include them in the sum. Just alternate totals.
		int integrationMode = 1;
		if (integrationMode == 0)
		{
			// Trapezoid sum. This will under-estimate the sum.
			sum = 0;
			for (int i = mini; i < maxi; i++)
			{
				if (pg[i] == 0)
				{
					// No probability so the function is zero.
					// This is the equivalent of only computing the Fisher information over 
					// the valid range of z.
					continue;
				}
				final double f = Maths.pow2(pg[i + offseti]) / pg[i];
				sum += f;
			}
		}
		else if (integrationMode == 1)
		{
			// This computes the sum as:
			// 3h/8 * [ f(x0) + 3f(x1) + 3f(x2) + 2f(x3) + 3f(x4) + 3f(x5) + 2f(x6) + ... + f(xn) ]
			double sum3 = 0, sum2 = 0;
			for (int i = mini; i < maxi; i++)
			{
				if (pg[i] == 0)
				{
					// No probability so the function is zero.
					// This is the equivalent of only computing the Fisher information over 
					// the valid range of z.
					continue;
				}
				final double f = Maths.pow2(pg[i + offseti]) / pg[i];
				if (i % 3 == 2)
					sum2 += f;
				else
					sum3 += f;
			}

			sum = (3 * h / 8) * (sum3 * 3 + sum2 * 2);
		}
		else
		{
			// This computes the sum as:
			// h/3 * [ f(x0) + 4f(x1) + 2f(x2) + 4f(x3) + 2f(x4) ... + 4f(xn-1) + f(xn) ]
			double sum4 = 0, sum2 = 0;
			for (int i = mini; i < maxi; i++)
			{
				if (pg[i] == 0)
				{
					// No probability so the function is zero.
					// This is the equivalent of only computing the Fisher information over 
					// the valid range of z.
					continue;
				}
				final double f = Maths.pow2(pg[i + offseti]) / pg[i];
				if (i % 2 == 0)
					sum4 += f;
				else
					sum2 += f;
			}

			sum = (h / 3) * (sum4 * 4 + sum2 * 2);
		}

		// Subtract the final 1 
		sum -= 1;

		System.out.printf("t=%g x=%d-%d scale=%d   sum=%s  pI = %g   pgI = %g\n", t, minx, maxx, scale, sum, pI,
				1 / (t + s * s));

		// Check limits.

		// It should be worse than the Poisson Fisher information (upper limit) but 
		// better than the Poisson-Gaussian Fisher information (lower limit).
		// Note a low Fisher information is worse as this is the amount of information
		// carried about the parameter.

		return Maths.clip(pI, getPoissonGaussianI(t), sum);
		//return Maths.max(getPoissonGaussianI(t), sum);
	}

	/**
	 * Gets the poisson Fisher information.
	 *
	 * @param t
	 *            the poisson mean
	 * @return the poisson Fisher information
	 */
	private double getPoissonI(double t)
	{
		return 1.0 / t;
	}

	/**
	 * Gets the poisson-gaussian Fisher information.
	 * Approximate the Poisson as a Gaussian with u=t and var=t.
	 * Gaussian-Gaussian convolution: sa * sb => sc = sqrt(sa^2+sb^2).
	 * Fisher information of Gaussian mean is 1/variance.
	 *
	 * @param t
	 *            the poisson mean
	 * @return the poisson gaussian Fisher information
	 */
	private double getPoissonGaussianI(double t)
	{
		return 1.0 / (t + s * s);
	}

	/**
	 * Gets the gaussian Fisher information for mean 0.
	 * Fisher information of Gaussian mean is 1/variance.
	 *
	 * @return the gaussian Fisher information
	 */
	private double getGaussianI()
	{
		return 1.0 / (s * s);
	}

	private int log2(int scale)
	{
		int bits = 30;
		while ((scale & (1 << bits)) == 0)
			bits--;
		return bits;
	}

	/**
	 * Gets the mean threshold for the switch to a Gaussian-Gaussian convolution.
	 *
	 * @return the mean threshold
	 */
	public double getMeanThreshold()
	{
		return meanThreshold;
	}

	/**
	 * Sets the mean threshold for the switch to a Gaussian-Gaussian convolution.
	 *
	 * @param meanThreshold
	 *            the new mean threshold
	 */
	public void setMeanThreshold(double meanThreshold)
	{
		this.meanThreshold = meanThreshold;
	}

	/**
	 * Gets the cumulative probability of the Poisson distribution that is used.
	 *
	 * @return the cumulative probability
	 */
	public double getCumulativeProbability()
	{
		return cumulativeProbability;
	}

	/**
	 * Sets the cumulative probability of the Poisson distribution that is used.
	 *
	 * @param cumulativeProbability
	 *            the new cumulative probability
	 */
	public void setCumulativeProbability(double cumulativeProbability)
	{
		if (!(cumulativeProbability > 0 && cumulativeProbability <= 1))
			throw new IllegalArgumentException("P must be in the range 0-1");
		if (this.cumulativeProbability != cumulativeProbability)
		{
			this.cumulativeProbability = cumulativeProbability;
			if (limits == defaultLimits)
			{
				limits = new int[defaultLimits.length];
				tinyLimits = new int[tinyLimits.length];
			}
			else
			{
				// Reset
				Arrays.fill(limits, 0);
				Arrays.fill(tinyLimits, 0);
			}
		}
	}
}