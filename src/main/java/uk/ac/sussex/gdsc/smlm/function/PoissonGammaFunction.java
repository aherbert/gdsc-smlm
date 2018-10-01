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
package uk.ac.sussex.gdsc.smlm.function;

import org.apache.commons.math3.util.FastMath;

/**
 * Implements the probability density function for a Poisson-Gamma Mixture.
 * <p>
 * The implementation uses the Poisson-Gamma mixture described from Ulbrich &amp; Isacoff (2007). Nature Methods 4,
 * 319-321,
 * SI equation 3:<br>
 * <blockquote>
 * {@code Gp,m(0) = e^-p}<br>
 * {@code Gp,m(c|c>0) = sqrt(p/(c*m)) * e^(-c/m -p) * I1(2*sqrt(c*p/m))}<br>
 * </blockquote>
 * Where:<br>
 * c = the observed value at the pixel <br>
 * p = the function value (expected number of photons) <br>
 * m = the gain of the pixel <br>
 * I1 = Modified Bessel function of the first kind <br>
 * <p>
 * The likelihood function is designed to model on-chip amplification of a EMCCD/CCD/sCMOS camera which captures a
 * Poisson process of emitted light, converted to electrons on the camera chip, amplified by a gain and then read with
 * Gaussian noise.
 */
public class PoissonGammaFunction implements LikelihoodFunction, LogLikelihoodFunction, GradientLikelihoodFunction
{
    /**
     * The on-chip gain multiplication factor
     */
    final double m;

    /**
     * Instantiates a new poisson gamma function.
     *
     * @param m
     *            The on-chip gain multiplication factor
     * @throws IllegalArgumentException
     *             if the gain is zero or below
     */
    public PoissonGammaFunction(double m)
    {
        if (!(m > 0))
            throw new IllegalArgumentException("Gain must be strictly positive");
        this.m = m;
    }

    /**
     * Creates the with standard deviation.
     *
     * @param alpha
     *            The inverse of the on-chip gain multiplication factor
     * @return the poisson gamma function
     * @throws IllegalArgumentException
     *             if the gain is zero or below
     */
    public static PoissonGammaFunction createWithAlpha(final double alpha)
    {
        return new PoissonGammaFunction(1.0 / alpha);
    }

    private static final double twoPi = 2 * Math.PI;

    /**
     * Calculate the probability density function for a Poisson-Gamma distribution model of EM-gain.
     * <p>
     * See Ulbrich &amp; Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
     *
     * @param c
     *            The count to evaluate
     * @param p
     *            The average number of photons per pixel input to the EM-camera (must be positive)
     * @param m
     *            The multiplication factor (gain)
     * @return The probability
     */
    public static double poissonGamma(double c, double p, double m)
    {
        // Any observed count above zero
        if (c > 0.0)
        {
            // The observed count converted to photons
            final double c_m = c / m;

            // The current implementation of Bessel.II(x) is Infinity at x==710.
            // Also exp(-c/m -p) will be sub-normal at < -709.
            // Switch to an approximation.
            final double x = 2 * Math.sqrt(p * c_m);
            final double _c_m_p = -c_m - p;
            if (x > 709 || _c_m_p < -709)
                //return FastMath.exp(0.5 * Math.log(p / (c * m)) + _c_m_p + x - 0.5 * Math.log(twoPi * x));
                return (x / (2 * c)) * FastMath.exp(_c_m_p + x - 0.5 * Math.log(twoPi * x));
            //return Math.sqrt(p / (c * m)) * FastMath.exp(_c_m_p) * Bessel.I1(x);
            return (x / (2 * c)) * FastMath.exp(_c_m_p) * Bessel.I1(x);
        }
        else if (c == 0.0)
            return FastMath.exp(-p) * (1 + p / m);
        else
            return 0;
    }

    /**
     * Calculate the probability density function for a Poisson-Gamma distribution model of EM-gain for observed Poisson
     * counts. This avoids the computation of the Dirac delta function at c=0.
     * <p>
     * This method is suitable for use in integration routines.
     * <p>
     * If c==0 then the true probability is obtained by adding Math.exp(-p).
     * <p>
     * See Ulbrich &amp; Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
     *
     * @param c
     *            The count to evaluate
     * @param p
     *            The average number of photons per pixel input to the EM-camera (must be positive)
     * @param m
     *            The multiplication factor (gain)
     * @return The probability function for observed Poisson counts
     * @see #poissonGamma(double, double, double)
     * @see #dirac(double)
     */
    public static double poissonGammaN(double c, double p, double m)
    {
        // As above with no Dirac delta function at c=0

        if (c > 0.0)
        {
            final double c_m = c / m;
            final double x = 2 * Math.sqrt(p * c_m);
            final double _c_m_p = -c_m - p;
            if (x > 709 || _c_m_p < -709)
                //return FastMath.exp(0.5 * Math.log(p / (c * m)) + _c_m_p + x - 0.5 * Math.log(twoPi * x));
                return (x / (2 * c)) * FastMath.exp(_c_m_p + x - 0.5 * Math.log(twoPi * x));
            //return Math.sqrt(p / (c * m)) * FastMath.exp(_c_m_p) * Bessel.I1(x);
            return (x / (2 * c)) * FastMath.exp(_c_m_p) * Bessel.I1(x);
        }
        else if (c == 0.0)
            // No Dirac delta function
            return FastMath.exp(-p) * p / m;
        else
            return 0;
    }

    /**
     * Calculate the probability density function for a Poisson-Gamma distribution model of EM-gain for no observed
     * Poisson counts. This is the Dirac delta function at c=0.
     * <p>
     * See Ulbrich &amp; Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
     *
     * @param p
     *            The average number of photons per pixel input to the EM-camera (must be positive)
     * @return The probability function for observed Poisson counts
     * @see #poissonGamma(double, double, double)
     */
    public static double dirac(double p)
    {
        return FastMath.exp(-p);
    }

    /**
     * Calculate the probability density function for a Poisson-Gamma distribution model of EM-gain for no observed
     * Poisson counts. This is the Dirac delta function at c=0.
     * <p>
     * See Ulbrich &amp; Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
     *
     * @param p
     *            The average number of photons per pixel input to the EM-camera (must be positive)
     * @param dG_dp
     *            the gradient of the function G(c) with respect to parameter p
     * @return The probability function for observed Poisson counts
     * @see #poissonGamma(double, double, double)
     */
    public static double dirac(double p, double[] dG_dp)
    {
        final double exp_p = FastMath.exp(-p);
        dG_dp[0] = -exp_p;
        return exp_p;
    }

    /**
     * Calculate the probability density function for a Poisson-Gamma distribution model of EM-gain.
     * <p>
     * See Ulbrich &amp; Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
     * <p>
     * Note: This implementation will underestimate the cumulative probability ({@code sum<1}) when the mean is close to
     * 1 and the gain is low ({@code <10}).
     *
     * @param c
     *            The count to evaluate
     * @param p
     *            The average number of photons per pixel input to the EM-camera (must be positive)
     * @param m
     *            The multiplication factor (gain)
     * @param dG_dp
     *            the gradient of the function G(c) with respect to parameter p
     * @return The probability
     */
    public static double poissonGamma(double c, double p, double m, double[] dG_dp)
    {
        // Any observed count above zero
        if (c > 0.0)
        {
            // The observed count converted to photons
            final double c_m = c / m;
            final double cp_m = p * c_m;

            // The current implementation of Bessel.II(x) is Infinity at x==710
            // due to the use of Math.exp(x). Switch to an approximation.
            final double x = 2 * Math.sqrt(cp_m);
            final double _c_m_p = -c_m - p;
            if (x > 709 || _c_m_p < -709)
            {
                // Approximate Bessel function i0(x)/i1(x) when using large x:
                // In(x) ~ exp(x)/sqrt(2*pi*x)
                final double exp_transform = FastMath.exp(_c_m_p + x - 0.5 * Math.log(twoPi * x));
                final double G = (x / (2 * c)) * exp_transform;
                dG_dp[0] = exp_transform / m - G;
                return G;
            }

            // G(c) = e^-p . e^-c/m . sum n=1 to inf { 1/(n!(n-1)!) . p^n c^(n-1) / m^n }
            // dG(c)/dp = e^-p . e^-c/m . sum n=1 to inf { 1/(n!(n-1)!) . n * p^(n-1) c^(n-1) / m^n } - e^-p . e^-c/m . sum n=1 to inf { 1/(n!(n-1)!) . p^n c^(n-1) / m^n }
            // dG(c)/dp = e^-p . e^-c/m . 1/m . sum n=1 to inf { 1/((n-1)!(n-1)!) . p^(n-1) c^(n-1) / m^(n-1) } - G(c)
            // dG(c)/dp = e^-p . e^-c/m . 1/m . sum n=0 to inf { 1/(n!^2) . (pc/m)^n } - G(c)

            // Bessel I0 = sum n=0 to inf { 1/(n!^2) . ((x/2)^2)^n }
            // x = 2 * sqrt(cp/m)

            // dG(c)/dp = e^-p . e^-c/m . 1/m . I0(2*sqrt(cp/m)) - G(c)
            // dG(c)/dp = e^(-c/m -p) . I0(2*sqrt(cp/m))/m - G(c)

            final double exp_c_m_p = FastMath.exp(_c_m_p);
            //double G = Math.sqrt(p / (c * m)) * exp_c_m_p * Bessel.I1(x);
            final double G = (x / (2 * c)) * exp_c_m_p * Bessel.I1(x);
            dG_dp[0] = exp_c_m_p * Bessel.I0(x) / m - G;
            return G;
        }
        else if (c == 0.0)
        {
            // f(p) = exp(-p) * (1 + p / m)
            // df/dp = (-exp(-p) * (1 + p / m)) + (exp(-p) / m)
            final double exp_p = FastMath.exp(-p);
            final double G = exp_p * (1 + p / m);
            dG_dp[0] = exp_p / m - G;
            return G;
        }
        else
        {
            dG_dp[0] = 0;
            return 0;
        }
    }

    /**
     * Calculate the probability density function for a Poisson-Gamma distribution model of EM-gain for observed Poisson
     * counts. This avoids the computation of the Dirac delta function at c=0.
     * <p>
     * This method is suitable for use in integration routines.
     * <p>
     * If c==0 then the true probability is obtained by adding Math.exp(-p).
     * <p>
     * See Ulbrich &amp; Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
     *
     * @param c
     *            The count to evaluate
     * @param p
     *            The average number of photons per pixel input to the EM-camera (must be positive)
     * @param m
     *            The multiplication factor (gain)
     * @param dG_dp
     *            the gradient of the function G(c) with respect to parameter p
     * @return The probability function for observed Poisson counts
     * @see #poissonGamma(double, double, double)
     * @see #dirac(double)
     */
    public static double poissonGammaN(double c, double p, double m, double[] dG_dp)
    {
        // As above with no Dirac delta function at c=0

        if (c > 0.0)
        {
            final double c_m = c / m;
            final double cp_m = p * c_m;
            final double x = 2 * Math.sqrt(cp_m);
            final double _c_m_p = -c_m - p;
            if (x > 709 || _c_m_p < -709)
            {
                final double exp_transform = FastMath.exp(_c_m_p + x - 0.5 * Math.log(twoPi * x));
                final double G = (x / (2 * c)) * exp_transform;
                dG_dp[0] = exp_transform / m - G;
                return G;
            }
            final double exp_c_m_p = FastMath.exp(_c_m_p);
            final double G = (x / (2 * c)) * exp_c_m_p * Bessel.I1(x);
            dG_dp[0] = exp_c_m_p * Bessel.I0(x) / m - G;
            return G;
        }
        else if (c == 0.0)
        {
            // No Dirac delta function
            final double exp_p_m = FastMath.exp(-p) / m;
            final double G = exp_p_m * p;
            dG_dp[0] = exp_p_m - G;
            return G;
        }
        else
        {
            dG_dp[0] = 0;
            return 0;
        }
    }

    /**
     * Calculate the probability density function for a Poisson-Gamma distribution model of EM-gain.
     * <p>
     * See Ulbrich &amp; Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
     * <p>
     * This is a special version which computes only part of the gradient.
     * The partial gradient is equal to the actual gradient plus the value of the function.
     *
     * @param c
     *            The count to evaluate
     * @param p
     *            The average number of photons per pixel input to the EM-camera (must be positive)
     * @param m
     *            The multiplication factor (gain)
     * @param dG_dp
     *            the partial gradient of the function G(c) with respect to parameter p
     * @return The probability
     */
    static double poissonGammaPartial(double c, double p, double m, double[] dG_dp)
    {
        // As above but do not subtract the function value G from the gradient.
        if (c > 0.0)
        {
            final double c_m = c / m;
            final double cp_m = p * c_m;
            final double x = 2 * Math.sqrt(cp_m);
            final double _c_m_p = -c_m - p;
            if (x > 709 || _c_m_p < -709)
            {
                final double exp_transform = FastMath.exp(_c_m_p + x - 0.5 * Math.log(twoPi * x));
                final double G = (x / (2 * c)) * exp_transform;
                dG_dp[0] = exp_transform / m;
                return G;
            }
            final double exp_c_m_p = FastMath.exp(_c_m_p);
            //double G = Math.sqrt(p / (c * m)) * exp_c_m_p * Bessel.I1(x);
            final double G = (x / (2 * c)) * exp_c_m_p * Bessel.I1(x);
            dG_dp[0] = exp_c_m_p * Bessel.I0(x) / m;
            return G;
        }
        else if (c == 0.0)
        {
            final double exp_p = FastMath.exp(-p);
            final double G = exp_p * (1 + p / m);
            dG_dp[0] = exp_p / m;
            return G;
        }
        else
        {
            dG_dp[0] = 0;
            return 0;
        }
    }

    /**
     * Calculate the an unscaled probability density function for a Poisson-Gamma distribution model of EM-gain.
     * <p>
     * See Ulbrich &amp; Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
     * <p>
     * This is unscaled as the factor exp^p has been removed. This stabilises computation for large p.
     * <p>
     * This is a special version which computes only part of the gradient.
     * The partial gradient is equal to the actual gradient plus the value of the function.
     *
     * @param c
     *            The count to evaluate
     * @param p
     *            The average number of photons per pixel input to the EM-camera (must be positive)
     * @param m
     *            The multiplication factor (gain)
     * @param dG_dp
     *            the partial gradient of the function G(c) with respect to parameter p
     * @return The unscaled probability
     */
    static double unscaledPoissonGammaPartial(double c, double p, double m, double[] dG_dp)
    {
        // As above but:
        // - do not multiply by exp^-p
        // - do not subtract the function value G from the gradient.
        if (c > 0.0)
        {
            final double c_m = c / m;
            final double cp_m = p * c_m;
            final double x = 2 * Math.sqrt(cp_m);
            if (x > 709 || -c_m < -709)
            {
                final double exp_transform = FastMath.exp(-c_m + x - 0.5 * Math.log(twoPi * x));
                final double G = (x / (2 * c)) * exp_transform;
                dG_dp[0] = exp_transform / m;
                return G;
            }
            final double exp_c_m = FastMath.exp(-c_m);
            final double G = (x / (2 * c)) * exp_c_m * Bessel.I1(x);
            dG_dp[0] = exp_c_m * Bessel.I0(x) / m;
            return G;
        }
        else if (c == 0.0)
        {
            final double G = 1 + p / m;
            dG_dp[0] = 1 / m;
            return G;
        }
        else
        {
            dG_dp[0] = 0;
            return 0;
        }
    }

    /**
     * Calculate the log probability density function for a Poisson-Gamma distribution model of EM-gain.
     * <p>
     * See
     * . Nature Methods 4, 319-321, SI equation 3.
     *
     * @param c
     *            The count to evaluate
     * @param p
     *            The average number of photons per pixel input to the EM-camera (must be positive)
     * @param m
     *            The multiplication factor (gain)
     * @return The log probability
     */
    public static double logPoissonGamma(double c, double p, double m)
    {
        // As above without final exp
        if (c > 0.0)
        {
            final double c_m = c / m;
            final double cp_m = p * c_m;
            final double x = 2 * Math.sqrt(cp_m);
            if (x > 709)
                return 0.5 * Math.log(p / (c * m)) - c_m - p + x - 0.5 * Math.log(twoPi * x);
            return 0.5 * Math.log(p / (c * m)) - c_m - p + Math.log(Bessel.I1(x));
        }
        else if (c == 0.0)
            // log (FastMath.exp(-p) * (1 + p / m))
            return -p + Math.log(1 + p / m);
        else
            return Double.NEGATIVE_INFINITY;
    }

    /** {@inheritDoc} */
    @Override
    public double likelihood(final double o, final double e)
    {
        return poissonGamma(o, e, m);
    }

    /** {@inheritDoc} */
    @Override
    public double logLikelihood(double o, double e)
    {
        return logPoissonGamma(o, e, m);
    }

    /** {@inheritDoc} */
    @Override
    public double likelihood(double o, double t, double[] dp_dt)
    {
        return poissonGamma(o, t, m, dp_dt);
    }
}
