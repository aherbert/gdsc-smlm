/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.function;

import uk.ac.sussex.gdsc.smlm.utils.StdMath;

/**
 * Computes likelihood values for a Poisson function.
 */
public class PoissonCalculator {
  /** Avoid repeated computation of log of 2 PI. */
  private static final double HALF_LOG_2_PI = 0.9189385332046727417;

  /** The value of x where the instance method computes x! using an approximation. */
  public static final double APPROXIMATION_X = 1.5;

  // For computation of Stirling series

  /** The constant 1/12. */
  protected static final double ONE_OVER_12 = 1.0 / 12.0;
  /** The constant 1/360. */
  protected static final double ONE_OVER_360 = 1.0 / 360;
  /** The constant 1/1260. */
  protected static final double ONE_OVER_1260 = 1.0 / 1260;
  /** The constant 1/1680. */
  protected static final double ONE_OVER_1680 = 1.0 / 1680;

  private boolean uninitialsed = true;

  /** The maximum log likelihood. */
  protected double mll;
  /**
   * The sum of the log X factorial term to convert the pseudo log-likelihood to the log-likelihood.
   */
  protected double sumLogXFactorial;
  /** The Poisson data values x (must be positive). */
  protected final double[] x;

  /**
   * Instantiates a new poisson calculator. This pre-computes factors for the log-likelihood and
   * log-likelihood ratio. It should be used for repeat calls to determine the log-likelihood
   * (ratio) when the mean value is different but the data x is the same.
   *
   * @param x the values x (must be positive)
   * @throws IllegalArgumentException if the values are not positive
   */
  public PoissonCalculator(double[] x) {
    this.x = x;
  }

  /**
   * Compute the maximum log-likelihood and the log(x!) term prefactors.
   */
  protected void computePrefactors() {
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
    // + 1/12x
    // - 1/360x^3
    // + 1/1260x^5
    // - 1/1680x^7
    // + ...
    // The error in truncating the series is always of the opposite sign and at most
    // the same magnitude as the first omitted term.

    for (int i = 0, n = x.length; i < n; i++) {
      // Note: When x==0:
      // log(n!) = 0
      // mll = 0
      if (x[i] > 0.0) {
        final double logx = Math.log(x[i]);
        if (x[i] <= APPROXIMATION_X) {
          // At low values of log(n!) we use the gamma function as the relative error of the
          // approximation is high.
          final double logXFactorial = LogFactorial.value(x[i]);
          sumLogXFactorial += logXFactorial;
          mll += x[i] * logx - x[i] - logXFactorial;
        } else {
          // Approximate log(n!) using Stirling's approximation using the first 3 terms.
          // This will have a maximum relative error of approximately 2.87e-4
          final double ologx =
              HALF_LOG_2_PI + 0.5 * logx + ONE_OVER_12 / x[i] - ONE_OVER_360 / pow3(x[i]);
          sumLogXFactorial += x[i] * logx - x[i] + ologx;
          mll -= ologx;
        }
      } else if (x[i] != 0) {
        throw new IllegalArgumentException("Input values x must be positive");
      }
    }
  }

  /**
   * Raise the value to the power 3.
   *
   * @param value the value
   * @return value^3
   */
  protected static double pow3(double value) {
    return value * value * value;
  }

  /**
   * Compute log(x!) using the Sterling series approximation with the first n terms of the following
   * series.
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
   * @param x the value x
   * @param n the number of terms
   * @return x!
   */
  static double logFactorialApproximation(double x, int n) {
    final double logx = Math.log(x);
    double value = x * logx - x;
    if (n <= 0) {
      return value;
    }
    // First term
    // 0.5 * log(2*pi*n) = 0.5 * log(2*pi) + 0.5 * log(n)
    value += HALF_LOG_2_PI + 0.5 * logx;
    if (n <= 1) {
      return value;
    }
    // Second term
    value += ONE_OVER_12 / x;
    if (n <= 2) {
      return value;
    }
    // Third term
    final double x2 = x * x;
    x *= x2;
    value -= ONE_OVER_360 / x;
    if (n <= 3) {
      return value;
    }
    // Fourth term
    x *= x2;
    value += ONE_OVER_1260 / x;
    if (n <= 4) {
      return value;
    }
    // Fifth term
    x *= x2;
    value -= ONE_OVER_1680 / x;
    return value;
  }

  /**
   * Get the pseudo Poisson log likelihood of value x given the mean. The mean must be strictly
   * positive.
   *
   * <p>The pseudo log-likelihood is equivalent to the log-likelihood without subtracting the
   * log(x!) term. It can be converted to the log-likelihood by subtracting
   * {@link #getLogXFactorialTerm()}.
   *
   * <p>This term is suitable for use in maximum likelihood routines.
   *
   * <pre>
   * pseudo ll = x * log(mean) - mean
   * </pre>
   *
   * @param mean the mean
   * @return the pseudo log likelihood
   */
  public double pseudoLogLikelihood(double[] mean) {
    double ll = 0.0;
    for (int i = mean.length; i-- > 0;) {
      if (x[i] == 0.0) {
        ll -= mean[i];
      } else {
        ll += x[i] * Math.log(mean[i]) - mean[i];
      }
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
  public double getLogXFactorialTerm() {
    initialise();
    return sumLogXFactorial;
  }

  private void initialise() {
    if (uninitialsed) {
      uninitialsed = false;
      computePrefactors();
    }
  }

  /**
   * Gets the values x.
   *
   * @return the values x
   */
  public double[] getX() {
    return x.clone();
  }

  /**
   * Get the Poisson maximum log likelihood of values x.
   *
   * @return the maximum log likelihood
   */
  public double getMaximumLogLikelihood() {
    initialise();
    return mll;
  }

  /**
   * Get the Poisson log likelihood ratio of the log likelihood. Note that the input must not be the
   * pseudo log-likelihood
   *
   * @param logLikelihood the log likelihood
   * @return the log likelihood ratio
   */
  public double getLogLikelihoodRatio(double logLikelihood) {
    initialise();
    // The log likelihood should be below the maximum log likelihood
    return (logLikelihood > mll) ? 0 : -2.0 * (logLikelihood - mll);
  }

  /**
   * Get the Poisson log likelihood of value x given the mean. The mean must be strictly positive.
   *
   * @param mean the mean
   * @return the log likelihood
   */
  public double logLikelihood(double[] mean) {
    return pseudoLogLikelihood(mean) - getLogXFactorialTerm();
  }

  /**
   * Get the Poisson log likelihood of value x given the mean. The mean must be strictly positive. x
   * must be positive.
   *
   * @param mean the mean
   * @param x the x
   * @return the log likelihood
   */
  public static double logLikelihood(double mean, double x) {
    if (x == 0.0) {
      return -mean;
    }
    return x * Math.log(mean) - mean - LogFactorial.value(x);
  }

  /**
   * Get the Poisson log likelihood of value x given the mean. The mean must be strictly positive. x
   * must be positive.
   *
   * @param mean the mean
   * @param x the x
   * @return the log likelihood
   */
  public static double logLikelihood(double[] mean, double[] x) {
    double ll = 0.0;
    for (int i = mean.length; i-- > 0;) {
      if (x[i] == 0.0) {
        ll -= mean[i];
      } else {
        ll += x[i] * Math.log(mean[i]) - mean[i] - LogFactorial.value(x[i]);
      }
    }
    return ll;
  }


  /**
   * Get the Poisson log likelihood of value x given the mean. The mean must be strictly positive. x
   * must be positive.
   *
   * @param mean the mean
   * @param x the x
   * @return the log likelihood
   */
  public static double fastLogLikelihood(double mean, double x) {
    // The log-likelihood (ll) is:
    // ll = (x==0) ? -mean : x * Math.log(mean) - mean - logFactorial(x)
    if (x > 0) {
      return fastLogLikelihoodX(mean, x);
    }
    return -mean;
  }

  /**
   * Get the Poisson log likelihood of value x given the mean. The mean must be strictly positive. x
   * must be positive.
   *
   * <p>Computation is done using an approximation to x! when x is above {@link #APPROXIMATION_X}.
   *
   * @param mean the mean
   * @param x the x
   * @return the log likelihood
   */
  public static double fastLogLikelihood(double[] mean, double[] x) {
    double ll = 0.0;
    for (int i = mean.length; i-- > 0;) {
      if (x[i] == 0.0) {
        ll -= mean[i];
      } else {
        ll += fastLogLikelihoodX(mean[i], x[i]);
      }
    }
    return ll;
  }

  /**
   * Get the Poisson log likelihood of value x given the mean. The mean must be strictly positive. x
   * must be positive.
   *
   * @param mean the mean
   * @param x the x
   * @param fastLog the fast log function
   * @return the log likelihood
   */
  public static double fastLogLikelihood(double mean, double x, FastLog fastLog) {
    // The log-likelihood (ll) is:
    // ll = (x==0) ? -mean : x * Math.log(mean) - mean - logFactorial(x)
    if (x > 0) {
      return fastLogLikelihoodX(mean, x, fastLog);
    }
    return -mean;
  }

  /**
   * Get the Poisson log likelihood of value x given the mean. The mean must be strictly positive. x
   * must be positive.
   *
   * <p>Computation is done using an approximation to x! when x is above {@link #APPROXIMATION_X}.
   *
   * @param mean the mean
   * @param x the x
   * @param fastLog the fast log function
   * @return the log likelihood
   */
  public static double fastLogLikelihood(double[] mean, double[] x, FastLog fastLog) {
    double ll = 0.0;
    for (int i = mean.length; i-- > 0;) {
      if (x[i] == 0.0) {
        ll -= mean[i];
      } else {
        ll += fastLogLikelihoodX(mean[i], x[i], fastLog);
      }
    }
    return ll;
  }

  /**
   * Get the Poisson log likelihood of value x given the mean. The mean and x must be strictly
   * positive.
   *
   * <p>Computation is done using an approximation to x! when x is above {@link #APPROXIMATION_X}.
   *
   * @param mean the mean
   * @param x the x
   * @return the log likelihood
   */
  private static double fastLogLikelihoodX(double mean, double x) {
    // The log-likelihood (ll) is:
    // ll = x * Math.log(mean) - mean - logFactorial(x)

    // Note that the logFactorial can be approximated as using Stirling's approximation:
    // https://en.m.wikipedia.org/wiki/Stirling%27s_approximation
    // log(x!) = x * log(x) - x + O(ln(x))
    // O(ln(x)) is a final term that can be computed using a series expansion.

    // This makes:
    // ll = x * log(mean) - mean - x * log(x) + x - O(ln(x))

    // Note: This can be rearranged:
    // ll = x * (log(mean) - log(x)) - mean + x - O(ln(x))
    // ll = x * log(mean/x) - mean + x - O(ln(x))
    // However the log(x) is needed in the O(ln(x)) computation.

    // Use the Stirling approximation when appropriate
    if (x <= APPROXIMATION_X) {
      // At low values of log(n!) we use the gamma function as the relative error of the
      // approximation is high.
      // Note that the LogFactorial function uses only 1 Math.log() call when the input is
      // below 2.5. Above that it uses 2 calls so we switch to the approximation.
      return x * Math.log(mean) - mean - LogFactorial.value(x);
    }
    // Approximate log(n!) using Stirling's approximation using the first 3 terms.
    // This will have a maximum relative error of approximately 6.7e-5
    // ll = x * log(mean) - mean - x * log(x) + x - O(ln(x))
    // O(ln(x)) = 0.5 * log(2*pi) + 0.5 * log(x) + 1/12x - 1/360x^3
    return x * Math.log(mean) - mean - HALF_LOG_2_PI - (x + 0.5) * Math.log(x) + x - ONE_OVER_12 / x
        + ONE_OVER_360 / pow3(x);
  }

  /**
   * Get the Poisson log likelihood of value x given the mean. The mean and x must be strictly
   * positive.
   *
   * <p>Computation is done using an approximation to x! when x is above {@link #APPROXIMATION_X}.
   *
   * @param mean the mean
   * @param x the x
   * @param fastLog the fast log function
   * @return the log likelihood
   */
  private static double fastLogLikelihoodX(double mean, double x, FastLog fastLog) {
    // The log-likelihood (ll) is:
    // ll = x * Math.log(mean) - mean - logFactorial(x)

    // Note that the logFactorial can be approximated as using Stirling's approximation:
    // https://en.m.wikipedia.org/wiki/Stirling%27s_approximation
    // log(x!) = x * log(x) - x + O(ln(x))
    // O(ln(x)) is a final term that can be computed using a series expansion.

    // This makes:
    // ll = x * log(mean) - mean - x * log(x) + x - O(ln(x))

    // Note: This can be rearranged:
    // ll = x * (log(mean) - log(x)) - mean + x - O(ln(x))
    // ll = x * log(mean/x) - mean + x - O(ln(x))
    // However the log(x) is needed in the O(ln(x)) computation.

    // Use the Stirling approximation when appropriate
    if (x <= APPROXIMATION_X) {
      // At low values of log(n!) we use the gamma function as the relative error of the
      // approximation is high.
      // Note that the LogFactorial function uses only 1 Math.log() call when the input is
      // below 2.5. Above that it uses 2 calls so we switch to the approximation.
      return x * fastLog.log(mean) - mean - LogFactorial.value(x);
    }
    // Approximate log(n!) using Stirling's approximation using the first 3 terms.
    // This will have a maximum relative error of approximately 6.7e-5
    // ll = x * log(mean) - mean - x * log(x) + x - O(ln(x))
    // O(ln(x)) = 0.5 * log(2*pi) + 0.5 * log(x) + 1/12x - 1/360x^3
    return x * fastLog.log(mean) - mean - HALF_LOG_2_PI - (x + 0.5) * fastLog.log(x) + x
        - ONE_OVER_12 / x + ONE_OVER_360 / pow3(x);
  }

  /**
   * Get the Poisson likelihood of value x given the mean. The mean must be strictly positive. x
   * must be positive.
   *
   * @param mean the mean
   * @param x the x
   * @return the likelihood
   */
  public static double likelihood(double mean, double x) {
    // This has a smaller range before computation fails:
    // return Math.pow(mean, x) * StdMath.exp(-mean) / factorial(x)
    return StdMath.exp(logLikelihood(mean, x));
  }

  /**
   * Get the Poisson likelihood of value x given the mean. The mean must be strictly positive. x
   * must be positive.
   *
   * @param mean the mean
   * @param x the x
   * @return the likelihood
   */
  public static double likelihood(double[] mean, double[] x) {
    return StdMath.exp(logLikelihood(mean, x));
  }

  /**
   * Get the Poisson likelihood of value x given the mean. The mean must be strictly positive. x
   * must be positive.
   *
   * @param mean the mean
   * @param x the x
   * @return the likelihood
   */
  public static double fastLikelihood(double mean, double x) {
    return StdMath.exp(fastLogLikelihood(mean, x));
  }

  /**
   * Get the Poisson likelihood of value x given the mean. The mean must be strictly positive. x
   * must be positive.
   *
   * <p>Computation is done using an approximation to x! when x is above {@link #APPROXIMATION_X}.
   *
   * @param mean the mean
   * @param x the x
   * @return the likelihood
   */
  public static double fastLikelihood(double[] mean, double[] x) {
    return StdMath.exp(fastLogLikelihood(mean, x));
  }

  /**
   * Get the Poisson likelihood of value x given the mean. The mean must be strictly positive. x
   * must be positive.
   *
   * @param mean the mean
   * @param x the x
   * @param fastLog the fast log function
   * @return the likelihood
   */
  public static double fastLikelihood(double mean, double x, FastLog fastLog) {
    return StdMath.exp(fastLogLikelihood(mean, x, fastLog));
  }

  /**
   * Get the Poisson likelihood of value x given the mean. The mean must be strictly positive. x
   * must be positive.
   *
   * <p>Computation is done using an approximation to x! when x is above {@link #APPROXIMATION_X}.
   *
   * @param mean the mean
   * @param x the x
   * @param fastLog the fast log function
   * @return the likelihood
   */
  public static double fastLikelihood(double[] mean, double[] x, FastLog fastLog) {
    return StdMath.exp(fastLogLikelihood(mean, x, fastLog));
  }

  /**
   * Get the Poisson maximum log likelihood of value x given the mean is value x. x must be
   * positive.
   *
   * <p>Computation is done using an approximation to x! when x is above {@link #APPROXIMATION_X}.
   *
   * @param x the x
   * @return the maximum log likelihood
   */
  public static double maximumLogLikelihood(double x) {
    return (x > 0.0) ? logLikelihood(x, x) : 0.0;
  }

  /**
   * Get the Poisson maximum log likelihood of value x given the mean is value x. x must be
   * positive.
   *
   * @param x the x
   * @return the maximum log likelihood
   */
  public static double maximumLogLikelihood(double[] x) {
    double ll = 0.0;
    for (int i = x.length; i-- > 0;) {
      ll += maximumLogLikelihood(x[i]);
    }
    return ll;
  }

  /**
   * Get the Poisson maximum log likelihood of value x given the mean is value x. x must be
   * positive.
   *
   * <p>Computation is done using an approximation to x! when x is above {@link #APPROXIMATION_X}.
   *
   * @param x the x
   * @return the maximum log likelihood
   */
  public static double fastMaximumLogLikelihood(double x) {
    return (x > 0.0) ? fastLogLikelihoodX(x, x) : 0.0;
  }

  /**
   * Get the Poisson maximum log likelihood of value x given the mean is value x. x must be
   * positive.
   *
   * <p>Computation is done using an approximation to x! when x is above {@link #APPROXIMATION_X}.
   *
   * @param x the x
   * @return the maximum log likelihood
   */
  public static double fastMaximumLogLikelihood(double[] x) {
    double ll = 0.0;
    for (int i = x.length; i-- > 0;) {
      if (x[i] > 0.0) {
        ll += fastLogLikelihoodX(x[i], x[i]);
      }
    }
    return ll;
  }

  /**
   * Get the Poisson maximum log likelihood of value x given the mean is value x. x must be
   * positive.
   *
   * <p>Computation is done using an approximation to x! when x is above {@link #APPROXIMATION_X}.
   *
   * @param x the x
   * @param fastLog the fast log function
   * @return the maximum log likelihood
   */
  public static double fastMaximumLogLikelihood(double x, FastLog fastLog) {
    return (x > 0.0) ? fastLogLikelihoodX(x, x, fastLog) : 0.0;
  }

  /**
   * Get the Poisson maximum log likelihood of value x given the mean is value x. x must be
   * positive.
   *
   * <p>Computation is done using an approximation to x! when x is above {@link #APPROXIMATION_X}.
   *
   * @param x the x
   * @param fastLog the fast log function
   * @return the maximum log likelihood
   */
  public static double fastMaximumLogLikelihood(double[] x, FastLog fastLog) {
    double ll = 0.0;
    for (int i = x.length; i-- > 0;) {
      if (x[i] > 0.0) {
        ll += fastLogLikelihoodX(x[i], x[i], fastLog);
      }
    }
    return ll;
  }

  /**
   * Get the Poisson maximum likelihood of value x given the mean is value x. x must be positive.
   *
   * @param x the x
   * @return the maximum likelihood
   */
  public static double maximumLikelihood(double x) {
    return (x > 0.0) ? likelihood(x, x) : 1;
  }

  /**
   * Get the Poisson maximum likelihood of value x given the mean is value x. x must be positive.
   *
   * @param x the x
   * @return the maximum likelihood
   */
  public static double maximumLikelihood(double[] x) {
    return StdMath.exp(maximumLogLikelihood(x));
  }

  /**
   * Get the Poisson maximum likelihood of value x given the mean is value x. x must be positive.
   *
   * @param x the x
   * @return the maximum likelihood
   */
  public static double fastMaximumLikelihood(double x) {
    return (x > 0.0) ? StdMath.exp(fastLogLikelihoodX(x, x)) : 1;
  }

  /**
   * Get the Poisson maximum likelihood of value x given the mean is value x. x must be positive.
   *
   * @param x the x
   * @return the maximum likelihood
   */
  public static double fastMaximumLikelihood(double[] x) {
    return StdMath.exp(fastMaximumLogLikelihood(x));
  }

  /**
   * Get the Poisson maximum likelihood of value x given the mean is value x. x must be positive.
   *
   * @param x the x
   * @param fastLog the fast log function
   * @return the maximum likelihood
   */
  public static double fastMaximumLikelihood(double x, FastLog fastLog) {
    return (x > 0.0) ? StdMath.exp(fastLogLikelihoodX(x, x, fastLog)) : 1;
  }

  /**
   * Get the Poisson maximum likelihood of value x given the mean is value x. x must be positive.
   *
   * @param x the x
   * @param fastLog the fast log function
   * @return the maximum likelihood
   */
  public static double fastMaximumLikelihood(double[] x, FastLog fastLog) {
    return StdMath.exp(fastMaximumLogLikelihood(x, fastLog));
  }

  /**
   * Get the Poisson log likelihood ratio of value x given the mean. The mean must be strictly
   * positive. x must be positive.
   *
   * @param mean the mean
   * @param x the x
   * @return the log likelihood ratio
   */
  public static double logLikelihoodRatio(double[] mean, double[] x) {
    // From https://en.wikipedia.org/wiki/Likelihood-ratio_test#Use:
    // LLR = -2 * [ ln(likelihood for alternative model) - ln(likelihood for null model)]
    // The model with more parameters (here alternative) will always fit at least as well—
    // i.e., have the same or greater log-likelihood—than the model with fewer parameters
    // (here null)

    double ll = 0.0;
    for (int i = mean.length; i-- > 0;) {
      if (x[i] > 0.0) {
        // ll += (x[i] * Math.log(mean[i]) - mean[i]) - (x[i] * Math.log(x[i]) - x[i])
        // ll += x[i] * Math.log(mean[i]) - mean[i] - x[i] * Math.log(x[i]) + x[i]
        // ll += x[i] * (Math.log(mean[i]) - Math.log(x[i])) - mean[i] + x[i]
        ll += x[i] * Math.log(mean[i] / x[i]) - mean[i] + x[i];
      } else {
        ll -= mean[i];
      }
    }
    return -2.0 * ll;
  }

  /**
   * Get the Poisson log likelihood ratio of value x given the mean. The mean must be strictly
   * positive. x must be positive.
   *
   * @param mean the mean
   * @param x the x
   * @param fastLog the fast log function
   * @return the log likelihood ratio
   */
  public static double logLikelihoodRatio(double[] mean, double[] x, FastLog fastLog) {
    double ll = 0.0;
    for (int i = mean.length; i-- > 0;) {
      if (x[i] > 0.0) {
        ll += x[i] * fastLog.log(mean[i] / x[i]) - mean[i] + x[i];
      } else {
        ll -= mean[i];
      }
    }
    return -2.0 * ll;
  }

  /**
   * Get the Poisson log likelihood ratio of value x given the mean. The mean must be strictly
   * positive. x must be positive.
   *
   * @param mean the mean
   * @param x the x
   * @return the log likelihood ratio
   */
  public static double logLikelihoodRatio(double mean, double x) {
    if (x > 0.0) {
      return -2.0 * x * Math.log(mean / x) - mean + x;
    }
    return -2.0 * mean;
  }

  /**
   * Get the Poisson log likelihood ratio of value x given the mean. The mean must be strictly
   * positive. x must be positive.
   *
   * @param mean the mean
   * @param x the x
   * @param fastLog the fast log function
   * @return the log likelihood ratio
   */
  public static double logLikelihoodRatio(double mean, double x, FastLog fastLog) {
    if (x > 0.0) {
      return -2.0 * x * fastLog.log(mean / x) - mean + x;
    }
    return -2.0 * mean;
  }
}
