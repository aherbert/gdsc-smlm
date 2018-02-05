package gdsc.smlm.function;

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;

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
 * Computes likelihood values for a Poisson function
 */
public class PoissonCalculator
{
	/** Avoid repeated computation of log of 2 PI */
	private static final double HALF_LOG_2_PI = 0.5 * Math.log(2.0 * FastMath.PI);

	/** The value of x where the instance method computes x! using an approximation. */
	public static final double APPROXIMATION_X = 1.5;

	// For computation of Stirling series
	private static final double ONE_OVER_12 = 1.0 / 12.0;
	private static final double ONE_OVER_360 = 1.0 / 360;
	private static final double ONE_OVER_1260 = 1.0 / 1260;
	private static final double ONE_OVER_1680 = 1.0 / 1680;

	private double mll = Double.NaN, sumLogXFactorial;
	private final double[] x;

	/**
	 * Instantiates a new poisson calculator. This pre-computes factors for the log-likelihood and log-likelihood ratio.
	 * It should be used for repeat calls to determine the log-likelihood (ratio) when the mean value is different but
	 * the data x is the same.
	 *
	 * @param x
	 *            the values x (must be positive)
	 * @throws IllegalArgumentException
	 *             if the values are not positive
	 */
	public PoissonCalculator(double[] x)
	{
		this.x = x;
	}

	/**
	 * Compute the maximum log-likelihood and the log(x!) term prefactors.
	 */
	private void computePrefactors()
	{
		mll = 0;

		// The maximum log-likelihood (mll) is:
		// mll = (x==0) ? 0 : x * Math.log(x) - x - logFactorial(x)

		// Note that the logFactorial can be approximated as using Stirling's approximation:
		// https://en.m.wikipedia.org/wiki/Stirling%27s_approximation
		// log(x!) = x * log(x) - x + O(ln(x))
		// O(ln(x)) is a final term that can be computed using a series expansion.

		// This makes:
		// mll = x * log(x) - x - x * log(x) + x - O(ln(x))
		// mll = -O(ln(x))

		// The series can be written as:
		// O(ln(x)) = 0.5 * log(2*pi*x) + 
		//            + 1/12x 
		//            - 1/360x^3 
		//            + 1/1260x^5 
		//            - 1/1680x^7 
		//            + ...
		// The error in truncating the series is always of the opposite sign and at most 
		// the same magnitude as the first omitted term.

		for (int i = 0, n = x.length; i < n; i++)
		{
			// Note: When x==0:
			// log(n!) = 0
			// mll = 0
			if (x[i] > 0)
			{
				final double logx = Math.log(x[i]);
				if (x[i] <= 1)
				{
					// When x is below 1 then the factorial can be omitted, i.e. log(x!) = 0
					// Compute the actual maximum log likelihood
					mll += x[i] * logx - x[i];
				}
				else if (x[i] <= APPROXIMATION_X)
				{
					// At low values of log(n!) we use the gamma function as the relative error of the 
					// approximation is high. 
					// Note that the logGamma function uses only 1 Math.log() call when the input is 
					// below 2.5. Above that it uses 2 calls. So the cost for this accuracy is an extra 
					// Math.log() call.
					//final double logXFactorial = Gamma.logGamma(x[i] + 1);
					// Note Math.log1p is faster than FastMath.log1p.
					final double logXFactorial = -Math.log1p(Gamma.invGamma1pm1(x[i]));
					sumLogXFactorial += logXFactorial;
					mll += x[i] * logx - x[i] - logXFactorial;
				}
				else
				{
					// Approximate log(n!) using Stirling's approximation using the first 3 terms.
					// This will have a maximum relative error of approximately 2.87e-4 
					final double ologx = HALF_LOG_2_PI + 0.5 * logx + ONE_OVER_12 / x[i] - ONE_OVER_360 / pow3(x[i]);
					sumLogXFactorial += x[i] * logx - x[i] + ologx;
					mll -= ologx;
				}
			}
			else if (x[i] != 0)
			{
				throw new IllegalArgumentException("Input values x must be positive");
			}
		}
	}

	/**
	 * Raise the value to the power 3
	 *
	 * @param value
	 *            the value
	 * @return value^3
	 */
	private static double pow3(double value)
	{
		return value * value * value;
	}

	/**
	 * Compute log(x!) using the Sterling series approximation with the first n terms of the following series.
	 *
	 * <pre>
	 * ln(x!) = x * log(x) - x
	 *  + 0.5 * log(2*pi*x)   // term 1
	 *  + 1/12x               // term 2 ...
	 *  - 1/360x^3
	 *  + 1/1260x^5
	 *  - 1/1680x^7
	 * </pre>
	 * 
	 * @param x
	 *            the value x
	 * @param n
	 *            the number of terms
	 * @return x!
	 */
	public static double logFactorialApproximation(double x, int n)
	{
		if (x <= 1)
			return 0;
		final double logx = Math.log(x);
		double value = x * logx - x;
		if (n <= 0)
			return value;
		// First term
		// 0.5 * log(2*pi*n) = 0.5 * log(2*pi) + 0.5 * log(n)
		value += HALF_LOG_2_PI + 0.5 * logx;
		if (n <= 1)
			return value;
		// Second term
		value += ONE_OVER_12 / x;
		if (n <= 2)
			return value;
		// Third term
		final double x2 = x * x;
		x *= x2;
		value -= ONE_OVER_360 / x;
		if (n <= 3)
			return value;
		// Fourth term
		x *= x2;
		value += ONE_OVER_1260 / x;
		if (n <= 4)
			return value;
		// Fifth term
		x *= x2;
		value -= ONE_OVER_1680 / x;
		return value;
	}

	/**
	 * Get the pseudo Poisson log likelihood of value x given the mean. The mean must be strictly positive.
	 * <p>
	 * The pseudo log-likelihood is equivalent to the log-likelihood without subtracting the log(x!) term. It can be
	 * converted to the log-likelihood by subtracting {@link #getLogXFactorialTerm()}.
	 * <p>
	 * This term is suitable for use in maximum likelihood routines.
	 * 
	 * <pre>
	 * pseudo ll = x * log(u) - u
	 * </pre>
	 *
	 * @param u
	 *            the mean
	 * @return the pseudo log likelihood
	 */
	public double pseudoLogLikelihood(double[] u)
	{
		double ll = 0.0;
		for (int i = u.length; i-- > 0;)
		{
			if (x[i] == 0)
				ll -= u[i];
			else
				ll += x[i] * Math.log(u[i]) - u[i];
		}
		return ll;
	}

	/**
	 * Gets the log X factorial term to convert the pseudo log-likelihood to the log-likelihood.
	 *
	 * <pre>
	 * ll = pseudo ll - log(x!)
	 * </pre>
	 * 
	 * @return the log X factorial term
	 */
	public double getLogXFactorialTerm()
	{
		if (Double.isNaN(mll))
			computePrefactors();
		return sumLogXFactorial;
	}

	/**
	 * Get the Poisson log likelihood of value x given the mean. The mean must be strictly positive.
	 *
	 * @param u
	 *            the mean
	 * @return the log likelihood
	 */
	public double logLikelihood(double[] u)
	{
		return pseudoLogLikelihood(u) - getLogXFactorialTerm();
	}

	/**
	 * Gets the values x
	 *
	 * @return the values x
	 */
	public double[] getX()
	{
		return x.clone();
	}

	/**
	 * Get the Poisson maximum log likelihood of values x.
	 *
	 * @return the maximum log likelihood
	 */
	public double getMaximumLogLikelihood()
	{
		if (Double.isNaN(mll))
			computePrefactors();
		return mll;
	}

	/**
	 * Get the Poisson log likelihood ratio of the log likelihood. Note that the input must not be the pseudo
	 * log-likelihood
	 *
	 * @param logLikelihood
	 *            the log likelihood
	 * @return the log likelihood ratio
	 */
	public double getLogLikelihoodRatio(double logLikelihood)
	{
		if (Double.isNaN(mll))
			computePrefactors();
		// The log likelihood should be below the maximum log likelihood
		return (logLikelihood > mll) ? 0 : -2.0 * (logLikelihood - mll);
	}

	/**
	 * Get the Poisson log likelihood of value x given the mean. The mean must be strictly positive. x must be positive.
	 *
	 * @param u
	 *            the mean
	 * @param x
	 *            the x
	 * @return the log likelihood
	 */
	public static double logLikelihood(double u, double x)
	{
		if (x == 0)
			return -u;
		return x * Math.log(u) - u - logFactorial(x);
	}

	/**
	 * Get the Poisson log likelihood of value x given the mean. The mean must be strictly positive. x must be positive.
	 *
	 * @param u
	 *            the mean
	 * @param x
	 *            the x
	 * @return the log likelihood
	 */
	public static double logLikelihood(double[] u, double[] x)
	{
		double ll = 0.0;
		for (int i = u.length; i-- > 0;)
		{
			if (x[i] == 0)
				ll -= u[i];
			else
				ll += x[i] * Math.log(u[i]) - u[i] - logFactorial(x[i]);
		}
		return ll;
	}

	/**
	 * Compute log(x!) using logGamma(x+1)
	 *
	 * @param x
	 *            the value x
	 * @param n
	 *            the number of terms
	 * @return x!
	 */
	public static double logFactorial(double x)
	{
		if (x > 1)
		{
			//return Gamma.logGamma(1 + x);
			return logGamma(1 + x);
		}
		return 0.0;
	}

	/**
	 * Copied from Apache Commons FastMath. Removed support for NaN and x < 0.5.
	 * <p>
	 * Returns the value of log&nbsp;&Gamma;(x) for x&nbsp;&gt;&nbsp;0.
	 * </p>
	 * <p>
	 * For x &le; 8, the implementation is based on the double precision
	 * implementation in the <em>NSWC Library of Mathematics Subroutines</em>,
	 * {@code DGAMLN}. For x &gt; 8, the implementation is based on
	 * </p>
	 * <ul>
	 * <li><a href="http://mathworld.wolfram.com/GammaFunction.html">Gamma
	 * Function</a>, equation (28).</li>
	 * <li><a href="http://mathworld.wolfram.com/LanczosApproximation.html">
	 * Lanczos Approximation</a>, equations (1) through (5).</li>
	 * <li><a href="http://my.fit.edu/~gabdo/gamma.txt">Paul Godfrey, A note on
	 * the computation of the convergent Lanczos complex Gamma
	 * approximation</a></li>
	 * </ul>
	 *
	 * @param x
	 *            Argument.
	 * @return the value of {@code log(Gamma(x))}, {@code Double.NaN} if
	 *         {@code x <= 0.0}.
	 */
	private static double logGamma(double x)
	{
		//		if (Double.isNaN(x) || (x <= 0.0))
		//		{
		//			return Double.NaN;
		//		}
		//		else if (x < 0.5)
		//		{
		//			return Gamma.logGamma1p(x) - Math.log(x);
		//		}
		//		else 
		if (x <= 2.5)
		{
			return Gamma.logGamma1p((x - 0.5) - 0.5);
		}
		else if (x <= 8.0)
		{
			final int n = (int) FastMath.floor(x - 1.5);
			double prod = 1.0;
			for (int i = 1; i <= n; i++)
			{
				prod *= x - i;
			}
			return Gamma.logGamma1p(x - (n + 1)) + Math.log(prod);
		}
		else
		{
			double sum = Gamma.lanczos(x);
			double tmp = x + Gamma.LANCZOS_G + .5;
			return ((x + .5) * Math.log(tmp)) - tmp + HALF_LOG_2_PI + Math.log(sum / x);
		}
	}

	/**
	 * Get the Poisson log likelihood of value x given the mean. The mean must be strictly positive. x must be positive.
	 *
	 * @param u
	 *            the mean
	 * @param x
	 *            the x
	 * @return the log likelihood
	 */
	public static double fastLogLikelihood(double u, double x)
	{
		// The log-likelihood (ll) is:
		// ll = (x==0) ? -u : x * Math.log(u) - u - logFactorial(x)
		if (x > 0)
			return fastLogLikelihoodX(u, x);
		else
			return -u;
	}

	/**
	 * Get the Poisson log likelihood of value x given the mean. The mean and x must be strictly positive.
	 * <p>
	 * Computation is done using an approximation to x! when x is above {@link #APPROXIMATION_X}. The number of calls to
	 * Math.log() is 2 for all x over 1.
	 *
	 * @param u
	 *            the mean
	 * @param x
	 *            the x
	 * @return the log likelihood
	 */
	private static double fastLogLikelihoodX(double u, double x)
	{
		// The log-likelihood (ll) is:
		// ll = x * Math.log(u) - u - logFactorial(x)

		// Note that the logFactorial can be approximated as using Stirling's approximation:
		// https://en.m.wikipedia.org/wiki/Stirling%27s_approximation
		// log(x!) = x * log(x) - x + O(ln(x))
		// O(ln(x)) is a final term that can be computed using a series expansion.

		// This makes:
		// ll = x * log(u) - u - x * log(x) + x - O(ln(x))
		
		// Note: This can be rearranged:
		// ll = x * (log(u) - log(x)) - u + x - O(ln(x))
		// ll = x * log(u/x) - u + x - O(ln(x))
		// However the log(x) is needed in the O(ln(x)) computation. 

		// Use the Stirling approximation when appropriate
		if (x <= 1)
		{
			// When x is below 1 then the factorial can be omitted, i.e. log(x!) = 0
			// Compute the actual log likelihood
			return x * Math.log(u) - u;
		}
		else if (x <= APPROXIMATION_X)
		{
			// At low values of log(n!) we use the gamma function as the relative error of the 
			// approximation is high. 
			// Note that the logGamma function uses only 1 Math.log() call when the input is 
			// below 2.5. Above that it uses 2 calls so we switch to the approximation.
			//return x * Math.log(u) - u - Gamma.logGamma(x[i] + 1);
			// Note Math.log1p is faster than FastMath.log1p.
			return x * Math.log(u) - u + Math.log1p(Gamma.invGamma1pm1(x));
		}
		else
		{
			// Approximate log(n!) using Stirling's approximation using the first 3 terms.
			// This will have a maximum relative error of approximately 6.7e-5
			// ll = x * log(u) - u - x * log(x) + x - O(ln(x))
			// O(ln(x)) = 0.5 * log(2*pi) + 0.5 * log(x) + 1/12x - 1/360x^3
			final double logx = Math.log(x);
			return x * Math.log(u) - u - HALF_LOG_2_PI - (x + 0.5) * logx + x - ONE_OVER_12 / x +
					ONE_OVER_360 / pow3(x);
		}
	}

	/**
	 * Get the Poisson log likelihood of value x given the mean. The mean must be strictly positive. x must be positive.
	 * <p>
	 * Computation is done using an approximation to x! when x is above {@link #APPROXIMATION_X}. The number of calls to
	 * Math.log() is 2 for all x over 1.
	 *
	 * @param u
	 *            the mean
	 * @param x
	 *            the x
	 * @return the log likelihood
	 */
	public static double fastLogLikelihood(double[] u, double[] x)
	{
		double ll = 0.0;
		for (int i = u.length; i-- > 0;)
		{
			if (x[i] == 0)
				ll -= u[i];
			else
				ll += fastLogLikelihoodX(u[i], x[i]);
		}
		return ll;
	}

	/**
	 * Get the Poisson likelihood of value x given the mean. The mean must be strictly positive. x must be positive.
	 *
	 * @param u
	 *            the mean
	 * @param x
	 *            the x
	 * @return the likelihood
	 */
	public static double likelihood(double u, double x)
	{
		// This has a smaller range before computation fails:
		//return Math.pow(u, x) * FastMath.exp(-u) / factorial(x);
		return FastMath.exp(logLikelihood(u, x));
	}

	/**
	 * Get the Poisson likelihood of value x given the mean. The mean must be strictly positive. x must be positive.
	 *
	 * @param u
	 *            the mean
	 * @param x
	 *            the x
	 * @return the likelihood
	 */
	public static double likelihood(double[] u, double[] x)
	{
		return FastMath.exp(logLikelihood(u, x));
	}

	/**
	 * Get the Poisson likelihood of value x given the mean. The mean must be strictly positive. x must be positive.
	 *
	 * @param u
	 *            the mean
	 * @param x
	 *            the x
	 * @return the likelihood
	 */
	public static double fastLikelihood(double u, double x)
	{
		return FastMath.exp(fastLogLikelihood(u, x));
	}

	/**
	 * Get the Poisson likelihood of value x given the mean. The mean must be strictly positive. x must be positive.
	 * <p>
	 * Computation is done using an approximation to x! when x is above {@link #APPROXIMATION_X}. The number of calls to
	 * Math.log() is 2 for all x over 1.
	 *
	 * @param u
	 *            the mean
	 * @param x
	 *            the x
	 * @return the likelihood
	 */
	public static double fastLikelihood(double[] u, double[] x)
	{
		return FastMath.exp(fastLogLikelihood(u, x));
	}

	/**
	 * Get the Poisson maximum log likelihood of value x given the mean is value x. x must be positive.
	 * <p>
	 * Computation is done using an approximation to x! when x is above {@link #APPROXIMATION_X}. The number of calls to
	 * Math.log() is 2 for all x over 1.
	 *
	 * @param x
	 *            the x
	 * @return the maximum log likelihood
	 */
	public static double maximumLogLikelihood(double x)
	{
		return (x > 0.0) ? logLikelihood(x, x) : 0.0;
	}

	/**
	 * Get the Poisson maximum log likelihood of value x given the mean is value x. x must be positive.
	 *
	 * @param x
	 *            the x
	 * @return the maximum log likelihood
	 */
	public static double maximumLogLikelihood(double[] x)
	{
		double ll = 0.0;
		for (int i = x.length; i-- > 0;)
			ll += maximumLogLikelihood(x[i]);
		return ll;
	}

	/**
	 * Get the Poisson maximum log likelihood of value x given the mean is value x. x must be positive.
	 * <p>
	 * Computation is done using an approximation to x! when x is above {@link #APPROXIMATION_X}. The number of calls to
	 * Math.log() is 2 for all x over 1.
	 *
	 * @param x
	 *            the x
	 * @return the maximum log likelihood
	 */
	public static double fastMaximumLogLikelihood(double x)
	{
		return (x > 0.0) ? fastLogLikelihoodX(x, x) : 0.0;
	}

	/**
	 * Get the Poisson maximum log likelihood of value x given the mean is value x. x must be positive.
	 * <p>
	 * Computation is done using an approximation to x! when x is above {@link #APPROXIMATION_X}. The number of calls to
	 * Math.log() is 2 for all x over 1.
	 *
	 * @param x
	 *            the x
	 * @return the maximum log likelihood
	 */
	public static double fastMaximumLogLikelihood(double[] x)
	{
		double ll = 0.0;
		for (int i = x.length; i-- > 0;)
			if (x[i] > 0)
				ll += fastLogLikelihoodX(x[i], x[i]);
		return ll;
	}

	/**
	 * Get the Poisson maximum likelihood of value x given the mean is value x. x must be positive.
	 *
	 * @param x
	 *            the x
	 * @return the maximum likelihood
	 */
	public static double maximumLikelihood(double x)
	{
		return (x > 0.0) ? likelihood(x, x) : 1;
	}

	/**
	 * Get the Poisson maximum likelihood of value x given the mean is value x. x must be positive.
	 *
	 * @param x
	 *            the x
	 * @return the maximum likelihood
	 */
	public static double maximumLikelihood(double[] x)
	{
		return FastMath.exp(maximumLogLikelihood(x));
	}

	/**
	 * Get the Poisson maximum likelihood of value x given the mean is value x. x must be positive.
	 *
	 * @param x
	 *            the x
	 * @return the maximum likelihood
	 */
	public static double fastMaximumLikelihood(double x)
	{
		return (x > 0.0) ? FastMath.exp(fastLogLikelihoodX(x, x)) : 1;
	}

	/**
	 * Get the Poisson maximum likelihood of value x given the mean is value x. x must be positive.
	 *
	 * @param x
	 *            the x
	 * @return the maximum likelihood
	 */
	public static double fastMaximumLikelihood(double[] x)
	{
		return FastMath.exp(fastMaximumLogLikelihood(x));
	}

	/**
	 * Get the Poisson log likelihood ratio of value x given the mean. The mean must be strictly positive. x must be
	 * positive.
	 *
	 * @param u
	 *            the mean
	 * @param x
	 *            the x
	 * @return the log likelihood ratio
	 */
	public static double logLikelihoodRatio(double[] u, double[] x)
	{
		// From https://en.wikipedia.org/wiki/Likelihood-ratio_test#Use:
		// LLR = -2 * [ ln(likelihood for alternative model) - ln(likelihood for null model)]
		// The model with more parameters (here alternative) will always fit at least as well—
		// i.e., have the same or greater log-likelihood—than the model with fewer parameters 
		// (here null)

		double ll = 0.0;
		for (int i = u.length; i-- > 0;)
		{
			//ll += logLikelihood(u[i], x[i]) - maximumLogLikelihood(x[i]);

			if (x[i] > 0.0)
			{
				//ll += (x[i] * Math.log(u[i]) - u[i]) - (x[i] * Math.log(x[i]) - x[i]);
				//ll += x[i] * Math.log(u[i]) - u[i] - x[i] * Math.log(x[i]) + x[i];
				//ll += x[i] * (Math.log(u[i]) - Math.log(x[i])) - u[i] + x[i];
				ll += x[i] * Math.log(u[i] / x[i]) - u[i] + x[i];
			}
			else
			{
				ll -= u[i];
			}
		}
		return -2.0 * ll;
	}
}