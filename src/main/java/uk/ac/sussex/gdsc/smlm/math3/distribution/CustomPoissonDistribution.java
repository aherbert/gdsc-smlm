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
/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package uk.ac.sussex.gdsc.smlm.math3.distribution;

import org.apache.commons.math3.distribution.AbstractIntegerDistribution;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.MathUtils;

/**
 * Implementation of the Poisson distribution.
 * <p>
 * Copy of the org.apache.commons.math3.distribution.PoissonDistribution but modified to allow the Poisson mean to be
 * set using a property. The Normal distribution used for the approximation is now created only when necessary.
 *
 * @see <a href="http://en.wikipedia.org/wiki/Poisson_distribution">Poisson distribution (Wikipedia)</a>
 * @see <a href="http://mathworld.wolfram.com/CustomPoissonDistribution.html">Poisson distribution (MathWorld)</a>
 */
public class CustomPoissonDistribution extends AbstractIntegerDistribution
{
    /**
     * Default maximum number of iterations for cumulative probability calculations.
     *
     * @since 2.1
     */
    public static final int DEFAULT_MAX_ITERATIONS = 10000000;
    /**
     * Default convergence criterion.
     *
     * @since 2.1
     */
    public static final double DEFAULT_EPSILON = 1e-12;
    /** Serializable version identifier. Different from PoissonDistribution */
    private static final long serialVersionUID = 2117648478379468511L;
    /**
     * Distribution used to compute normal approximation.
     * This is dynamically created when needed.
     */
    private NormalDistribution normal;
    /** Distribution needed for the {@link #sample()} method. */
    private final ExponentialDistribution exponential;
    /** Mean of the distribution. */
    private double mean;
    /**
     * Maximum number of iterations for cumulative probability. Cumulative
     * probabilities are estimated using either Lanczos series approximation
     * of {@link Gamma#regularizedGammaP(double, double, double, int)}
     * or continued fraction approximation of
     * {@link Gamma#regularizedGammaQ(double, double, double, int)}.
     */
    private final int maxIterations;

    /** Convergence criterion for cumulative probability. */
    private final double epsilon;

    /** The uninitialised flag. Indicates that the factors have not been computed for random sampling */
    private boolean uninitialised;

    // Pre-computed for sampling the mean
    private double lambda;
    private double lambdaFractional;
    private double logLambda;
    private double logLambdaFactorial;
    private long y2;
    private double delta;
    private double halfDelta;
    private double twolpd;
    private double p1;
    private double p2;
    private double c1;

    /**
     * Creates a new Poisson distribution with specified mean.
     * <p>
     * <b>Note:</b> this constructor will implicitly create an instance of
     * {@link Well19937c} as random generator to be used for sampling only (see
     * {@link #sample()} and {@link #sample(int)}). In case no sampling is
     * needed for the created distribution, it is advised to pass {@code null}
     * as random generator via the appropriate constructors to avoid the
     * additional initialisation overhead.
     *
     * @param p
     *            the Poisson mean
     * @throws NotStrictlyPositiveException
     *             if {@code p <= 0}.
     */
    public CustomPoissonDistribution(double p) throws NotStrictlyPositiveException
    {
        this(p, DEFAULT_EPSILON, DEFAULT_MAX_ITERATIONS);
    }

    /**
     * Creates a new Poisson distribution with specified mean, convergence
     * criterion and maximum number of iterations.
     * <p>
     * <b>Note:</b> this constructor will implicitly create an instance of
     * {@link Well19937c} as random generator to be used for sampling only (see
     * {@link #sample()} and {@link #sample(int)}). In case no sampling is
     * needed for the created distribution, it is advised to pass {@code null}
     * as random generator via the appropriate constructors to avoid the
     * additional initialisation overhead.
     *
     * @param p
     *            Poisson mean.
     * @param epsilon
     *            Convergence criterion for cumulative probabilities.
     * @param maxIterations
     *            the maximum number of iterations for cumulative
     *            probabilities.
     * @throws NotStrictlyPositiveException
     *             if {@code p <= 0}.
     * @since 2.1
     */
    public CustomPoissonDistribution(double p, double epsilon, int maxIterations) throws NotStrictlyPositiveException
    {
        this(new Well19937c(), p, epsilon, maxIterations);
    }

    /**
     * Creates a new Poisson distribution with specified mean, convergence
     * criterion and maximum number of iterations.
     *
     * @param rng
     *            Random number generator.
     * @param p
     *            Poisson mean.
     * @throws NotStrictlyPositiveException
     *             if {@code p <= 0}.
     */
    public CustomPoissonDistribution(RandomGenerator rng, double p) throws NotStrictlyPositiveException
    {
        this(rng, p, DEFAULT_EPSILON, DEFAULT_MAX_ITERATIONS);
    }

    /**
     * Creates a new Poisson distribution with specified mean, convergence
     * criterion and maximum number of iterations.
     *
     * @param rng
     *            Random number generator.
     * @param p
     *            Poisson mean.
     * @param epsilon
     *            Convergence criterion for cumulative probabilities.
     * @param maxIterations
     *            the maximum number of iterations for cumulative
     *            probabilities.
     * @throws NotStrictlyPositiveException
     *             if {@code p <= 0}.
     * @since 3.1
     */
    public CustomPoissonDistribution(RandomGenerator rng, double p, double epsilon, int maxIterations)
            throws NotStrictlyPositiveException
    {
        super(rng);

        setMean(p);
        this.epsilon = epsilon;
        this.maxIterations = maxIterations;

        // Use the same RNG instance as the parent class.
        exponential = new ExponentialDistribution(rng, 1, ExponentialDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);
    }

    /**
     * Creates a new Poisson distribution with the specified mean and
     * convergence criterion.
     *
     * @param p
     *            Poisson mean.
     * @param epsilon
     *            Convergence criterion for cumulative probabilities.
     * @throws NotStrictlyPositiveException
     *             if {@code p <= 0}.
     * @since 2.1
     */
    public CustomPoissonDistribution(double p, double epsilon) throws NotStrictlyPositiveException
    {
        this(p, epsilon, DEFAULT_MAX_ITERATIONS);
    }

    /**
     * Creates a new Poisson distribution with the specified mean and maximum
     * number of iterations.
     *
     * @param p
     *            Poisson mean.
     * @param maxIterations
     *            Maximum number of iterations for cumulative
     *            probabilities.
     * @since 2.1
     */
    public CustomPoissonDistribution(double p, int maxIterations)
    {
        this(p, DEFAULT_EPSILON, maxIterations);
    }

    /**
     * Get the mean for the distribution.
     *
     * @return the mean for the distribution.
     */
    public double getMean()
    {
        return mean;
    }

    /**
     * Sets the mean.
     *
     * @param p
     *            Poisson mean.
     * @throws NotStrictlyPositiveException
     *             if {@code p <= 0}.
     */
    public void setMean(double p)
    {
        if (p <= 0)
            throw new NotStrictlyPositiveException(LocalizedFormats.MEAN, p);
        setMeanUnsafe(p);
    }

    /**
     * Sets the mean.
     * <p>
     * Does not throw an exception if mean is not strictly positive
     *
     * @param p
     *            Poisson mean.
     */
    public void setMeanUnsafe(double p)
    {
        mean = p;
        normal = null;
        uninitialised = true;
    }

    /** {@inheritDoc} */
    @Override
    public double probability(int x)
    {
        final double logProbability = logProbability(x);
        return logProbability == Double.NEGATIVE_INFINITY ? 0 : FastMath.exp(logProbability);
    }

    /** {@inheritDoc} */
    @Override
    public double logProbability(int x)
    {
        double ret;
        if (x < 0 || x == Integer.MAX_VALUE)
            ret = Double.NEGATIVE_INFINITY;
        else if (x == 0)
            ret = -mean;
        else
            ret = -SaddlePointExpansionCopy.getStirlingError(x) - SaddlePointExpansionCopy.getDeviancePart(x, mean) -
                    0.5 * FastMath.log(MathUtils.TWO_PI) - 0.5 * FastMath.log(x);
        return ret;
    }

    /** {@inheritDoc} */
    @Override
    public double cumulativeProbability(int x)
    {
        if (x < 0)
            return 0;
        if (x == Integer.MAX_VALUE)
            return 1;
        return Gamma.regularizedGammaQ((double) x + 1, mean, epsilon, maxIterations);
    }

    /**
     * Calculates the Poisson distribution function using a normal
     * approximation. The {@code N(mean, sqrt(mean))} distribution is used
     * to approximate the Poisson distribution. The computation uses
     * "half-correction" (evaluating the normal distribution function at
     * {@code x + 0.5}).
     *
     * @param x
     *            Upper bound, inclusive.
     * @return the distribution function value calculated using a normal
     *         approximation.
     */
    public double normalApproximateProbability(int x)
    {
        // Create dynamically
        if (normal == null)
            normal = new NormalDistribution(random, mean, FastMath.sqrt(mean),
                    NormalDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);
        // calculate the probability using half-correction
        return normal.cumulativeProbability(x + 0.5);
    }

    /**
     * {@inheritDoc}
     *
     * For mean parameter {@code p}, the mean is {@code p}.
     */
    @Override
    public double getNumericalMean()
    {
        return getMean();
    }

    /**
     * {@inheritDoc}
     *
     * For mean parameter {@code p}, the variance is {@code p}.
     */
    @Override
    public double getNumericalVariance()
    {
        return getMean();
    }

    /**
     * {@inheritDoc}
     *
     * The lower bound of the support is always 0 no matter the mean parameter.
     *
     * @return lower bound of the support (always 0)
     */
    @Override
    public int getSupportLowerBound()
    {
        return 0;
    }

    /**
     * {@inheritDoc}
     *
     * The upper bound of the support is positive infinity,
     * regardless of the parameter values. There is no integer infinity,
     * so this method returns {@code Integer.MAX_VALUE}.
     *
     * @return upper bound of the support (always {@code Integer.MAX_VALUE} for
     *         positive infinity)
     */
    @Override
    public int getSupportUpperBound()
    {
        return Integer.MAX_VALUE;
    }

    /**
     * {@inheritDoc}
     *
     * The support of this distribution is connected.
     *
     * @return {@code true}
     */
    @Override
    public boolean isSupportConnected()
    {
        return true;
    }

    /**
     * {@inheritDoc}
     * <p>
     * <strong>Algorithm Description</strong>:
     * </p>
     * <ul>
     * <li>For small means, uses simulation of a Poisson process
     * using Uniform deviates, as described
     * <a href="http://mathaa.epfl.ch/cours/PMMI2001/interactive/rng7.htm"> here</a>.
     * The Poisson process (and hence value returned) is bounded by 1000 * mean.
     * </li>
     * <li>For large means, uses the rejection algorithm described in
     * <blockquote>
     * Devroye, Luc. (1981).<i>The Computer Generation of Poisson Random Variables</i><br>
     * <strong>Computing</strong> vol. 26 pp. 197-207.<br>
     * </blockquote>
     * </li>
     * </ul>
     *
     * @return a random value.
     * @since 2.2
     */
    @Override
    public int sample()
    {
        return (int) FastMath.min(nextPoisson(mean), Integer.MAX_VALUE);
    }

    /**
     * @param meanPoisson
     *            Mean of the Poisson distribution.
     * @return the next sample.
     */
    private long nextPoisson(double meanPoisson)
    {
        final double pivot = 40.0d;
        if (meanPoisson < pivot)
        {
            final double p = FastMath.exp(-meanPoisson);
            long n = 0;
            double r = 1.0d;
            double rnd = 1.0d;

            while (n < 1000 * meanPoisson)
            {
                rnd = random.nextDouble();
                r *= rnd;
                if (r >= p)
                    n++;
                else
                    return n;
            }
            return n;
        }

        computeFactors();

        y2 = lambdaFractional < Double.MIN_VALUE ? 0 : nextPoisson(lambdaFractional);

        double x = 0;
        double y = 0;
        double v = 0;
        int a = 0;
        double t = 0;
        double qr = 0;
        double qa = 0;
        for (;;)
        {
            final double u = random.nextDouble();
            if (u <= p1)
            {
                final double n = random.nextGaussian();
                x = n * FastMath.sqrt(lambda + halfDelta) - 0.5d;
                if (x > delta || x < -lambda)
                    continue;
                y = x < 0 ? FastMath.floor(x) : FastMath.ceil(x);
                final double e = exponential.sample();
                v = -e - (n * n / 2) + c1;
            }
            else
            {
                if (u > p1 + p2)
                {
                    y = lambda;
                    break;
                }

                x = delta + (twolpd / delta) * exponential.sample();
                y = FastMath.ceil(x);
                v = -exponential.sample() - delta * (x + 1) / twolpd;
            }
            a = x < 0 ? 1 : 0;
            t = y * (y + 1) / (2 * lambda);
            if (v < -t && a == 0)
            {
                y = lambda + y;
                break;
            }
            qr = t * ((2 * y + 1) / (6 * lambda) - 1);
            qa = qr - (t * t) / (3 * (lambda + a * (y + 1)));
            if (v < qa)
            {
                y = lambda + y;
                break;
            }
            if (v > qr)
                continue;
            if (v < y * logLambda - CombinatoricsUtils.factorialLog((int) (y + lambda)) + logLambdaFactorial)
            {
                y = lambda + y;
                break;
            }
        }
        return y2 + (long) y;
    }

    private void computeFactors()
    {
        if (uninitialised)
        {
            lambda = FastMath.floor(mean);
            lambdaFractional = mean - lambda;
            logLambda = FastMath.log(lambda);
            logLambdaFactorial = CombinatoricsUtils.factorialLog((int) lambda);
            delta = FastMath.sqrt(lambda * FastMath.log(32 * lambda / FastMath.PI + 1));
            halfDelta = delta / 2;
            twolpd = 2 * lambda + delta;
            final double a1 = FastMath.sqrt(FastMath.PI * twolpd) * FastMath.exp(1 / (8 * lambda));
            final double a2 = (twolpd / delta) * FastMath.exp(-delta * (1 + delta) / twolpd);
            final double aSum = a1 + a2 + 1;
            p1 = a1 / aSum;
            p2 = a2 / aSum;
            c1 = 1 / (8 * lambda);
            uninitialised = false;
        }
    }
}
