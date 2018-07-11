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
package org.apache.commons.math3.distribution;

import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;

/**
 * Implementation of the Gamma distribution.
 * <p>
 * Copy of the org.apache.commons.math3.distribution.GammaDistribution but modified to allow the shape parameter to be
 * set using a property.
 *
 * @see <a href="http://en.wikipedia.org/wiki/Gamma_distribution">Gamma distribution (Wikipedia)</a>
 * @see <a href="http://mathworld.wolfram.com/GammaDistribution.html">Gamma distribution (MathWorld)</a>
 */
public class CustomGammaDistribution extends AbstractRealDistribution
{
	/**
	 * Default inverse cumulative probability accuracy.
	 *
	 * @since 2.1
	 */
	public static final double DEFAULT_INVERSE_ABSOLUTE_ACCURACY = 1e-9;
	/** Serializable version identifier. Different from GammaDistribution. */
	private static final long serialVersionUID = 6520774824964116838L;
	/** The shape parameter. */
	private double shape;
	/** The scale parameter. */
	private double scale;

	/** The uninitialised flag. Indicates that the shape factors have not been computed for the current shape. */
	private boolean uninitialised;

	/**
	 * The constant value of {@code shape + g + 0.5}, where {@code g} is the
	 * Lanczos constant {@link Gamma#LANCZOS_G}.
	 */
	private double shiftedShape;
	/**
	 * The constant value of {@code shape / scale * sqrt(e / (2 * pi * (shape + g + 0.5))) / L(shape)},
	 * where {@code L(shape)} is the Lanczos approximation returned by {@link Gamma#lanczos(double)}. This prefactor is
	 * used in {@link #density(double)}, when no overflow occurs with the natural
	 * calculation.
	 */
	private double densityPrefactor1;
	/**
	 * The constant value of {@code shape * sqrt(e / (2 * pi * (shape + g + 0.5))) / L(shape)},
	 * where {@code L(shape)} is the Lanczos approximation returned by {@link Gamma#lanczos(double)}. This prefactor is
	 * used in {@link #density(double)}, when overflow occurs with the natural
	 * calculation.
	 */
	private double densityPrefactor2;
	/**
	 * Lower bound on {@code y = x / scale} for the selection of the computation
	 * method in {@link #density(double)}. For {@code y <= minY}, the natural
	 * calculation overflows.
	 */
	private double minY;
	/**
	 * Upper bound on {@code log(y)} ({@code y = x / scale}) for the selection
	 * of the computation method in {@link #density(double)}. For {@code log(y) >= maxLogY}, the natural calculation
	 * overflows.
	 */
	private double maxLogY;
	/** Inverse cumulative probability accuracy. */
	private final double solverAbsoluteAccuracy;

	/**
	 * Creates a new gamma distribution with specified values of the shape and
	 * scale parameters.
	 *
	 * @param shape
	 *            the shape parameter
	 * @param scale
	 *            the scale parameter
	 * @throws NotStrictlyPositiveException
	 *             if {@code shape <= 0} or {@code scale <= 0}.
	 */
	public CustomGammaDistribution(double shape, double scale) throws NotStrictlyPositiveException
	{
		this(shape, scale, DEFAULT_INVERSE_ABSOLUTE_ACCURACY);
	}

	/**
	 * Creates a new gamma distribution with specified values of the shape and
	 * scale parameters.
	 *
	 * @param shape
	 *            the shape parameter
	 * @param scale
	 *            the scale parameter
	 * @param inverseCumAccuracy
	 *            the maximum absolute error in inverse
	 *            cumulative probability estimates (defaults to {@link #DEFAULT_INVERSE_ABSOLUTE_ACCURACY}).
	 * @throws NotStrictlyPositiveException
	 *             if {@code shape <= 0} or {@code scale <= 0}.
	 * @since 2.1
	 */
	public CustomGammaDistribution(double shape, double scale, double inverseCumAccuracy)
			throws NotStrictlyPositiveException
	{
		this(new Well19937c(), shape, scale, inverseCumAccuracy);
	}

	/**
	 * Creates a Gamma distribution.
	 *
	 * @param rng
	 *            Random number generator.
	 * @param shape
	 *            the shape parameter
	 * @param scale
	 *            the scale parameter
	 * @throws NotStrictlyPositiveException
	 *             if {@code shape <= 0} or {@code scale <= 0}.
	 */
	public CustomGammaDistribution(RandomGenerator rng, double shape, double scale) throws NotStrictlyPositiveException
	{
		this(rng, shape, scale, DEFAULT_INVERSE_ABSOLUTE_ACCURACY);
	}

	/**
	 * Creates a Gamma distribution.
	 *
	 * @param rng
	 *            Random number generator.
	 * @param shape
	 *            the shape parameter
	 * @param scale
	 *            the scale parameter
	 * @param inverseCumAccuracy
	 *            the maximum absolute error in inverse
	 *            cumulative probability estimates (defaults to {@link #DEFAULT_INVERSE_ABSOLUTE_ACCURACY}).
	 * @throws NotStrictlyPositiveException
	 *             if {@code shape <= 0} or {@code scale <= 0}.
	 * @since 3.1
	 */
	public CustomGammaDistribution(RandomGenerator rng, double shape, double scale, double inverseCumAccuracy)
			throws NotStrictlyPositiveException
	{
		super(rng);

		setShape(shape);
		setScale(scale);
		this.solverAbsoluteAccuracy = inverseCumAccuracy;
	}

	/**
	 * Returns the shape parameter of {@code this} distribution.
	 *
	 * @return the shape parameter
	 * @deprecated as of version 3.1, {@link #getShape()} should be preferred.
	 *             This method will be removed in version 4.0.
	 */
	@Deprecated
	public double getAlpha()
	{
		return shape;
	}

	/**
	 * Returns the shape parameter of {@code this} distribution.
	 *
	 * @return the shape parameter
	 * @since 3.1
	 */
	public double getShape()
	{
		return shape;
	}

	/**
	 * Set the shape parameter
	 *
	 * @param shape
	 * @throws NotStrictlyPositiveException
	 *             if {@code shape <= 0}
	 */
	public void setShape(double shape)
	{
		if (shape <= 0)
		{
			throw new NotStrictlyPositiveException(LocalizedFormats.SHAPE, shape);
		}
		setShapeUnsafe(shape);
	}

	/**
	 * Set the shape parameter
	 * <p>
	 * Does not throw an exception if shape is not strictly positive
	 *
	 * @param shape
	 */
	public void setShapeUnsafe(double shape)
	{
		this.shape = shape;
		uninitialised = true;
	}

	/**
	 * Returns the scale parameter of {@code this} distribution.
	 *
	 * @return the scale parameter
	 * @deprecated as of version 3.1, {@link #getScale()} should be preferred.
	 *             This method will be removed in version 4.0.
	 */
	@Deprecated
	public double getBeta()
	{
		return scale;
	}

	/**
	 * Returns the scale parameter of {@code this} distribution.
	 *
	 * @return the scale parameter
	 * @since 3.1
	 */
	public double getScale()
	{
		return scale;
	}

	/**
	 * Set the scale parameter
	 *
	 * @param scale
	 * @throws NotStrictlyPositiveException
	 *             if {@code scale <= 0}
	 */
	public void setScale(double scale)
	{
		if (scale <= 0)
		{
			throw new NotStrictlyPositiveException(LocalizedFormats.SCALE, scale);
		}
		setScaleUnsafe(scale);
	}

	/**
	 * Set the scale parameter
	 * <p>
	 * Does not throw an exception if scale is not strictly positive
	 *
	 * @param scale
	 */
	public void setScaleUnsafe(double scale)
	{
		this.scale = scale;
		uninitialised = true;
	}

	/** {@inheritDoc} */
	@Override
	public double density(double x)
	{
		/*
		 * The present method must return the value of
		 *
		 * 1 x a - x
		 * ---------- (-) exp(---)
		 * x Gamma(a) b b
		 *
		 * where a is the shape parameter, and b the scale parameter.
		 * Substituting the Lanczos approximation of Gamma(a) leads to the
		 * following expression of the density
		 *
		 * a e 1 y a
		 * - sqrt(------------------) ---- (-----------) exp(a - y + g),
		 * x 2 pi (a + g + 0.5) L(a) a + g + 0.5
		 *
		 * where y = x / b. The above formula is the "natural" computation, which
		 * is implemented when no overflow is likely to occur. If overflow occurs
		 * with the natural computation, the following identity is used. It is
		 * based on the BOOST library
		 * http://www.boost.org/doc/libs/1_35_0/libs/math/doc/sf_and_dist/html/math_toolkit/special/sf_gamma/igamma.html
		 * Formula (15) needs adaptations, which are detailed below.
		 *
		 * y a
		 * (-----------) exp(a - y + g)
		 * a + g + 0.5
		 * y - a - g - 0.5 y (g + 0.5)
		 * = exp(a log1pm(---------------) - ----------- + g),
		 * a + g + 0.5 a + g + 0.5
		 *
		 * where log1pm(z) = log(1 + z) - z. Therefore, the value to be
		 * returned is
		 *
		 * a e 1
		 * - sqrt(------------------) ----
		 * x 2 pi (a + g + 0.5) L(a)
		 * y - a - g - 0.5 y (g + 0.5)
		 * * exp(a log1pm(---------------) - ----------- + g).
		 * a + g + 0.5 a + g + 0.5
		 */
		if (x < 0)
		{
			return 0;
		}
		computeFactors();
		final double y = x / scale;
		if ((y <= minY) || (FastMath.log(y) >= maxLogY))
		{
			/*
			 * Overflow.
			 */
			final double aux1 = (y - shiftedShape) / shiftedShape;
			final double aux2 = shape * (FastMath.log1p(aux1) - aux1);
			final double aux3 = -y * (Gamma.LANCZOS_G + 0.5) / shiftedShape + Gamma.LANCZOS_G + aux2;
			return densityPrefactor2 / x * FastMath.exp(aux3);
		}
		/*
		 * Natural calculation.
		 */
		return densityPrefactor1 * FastMath.exp(-y) * FastMath.pow(y, shape - 1);
	}

	private void computeFactors()
	{
		if (uninitialised)
		{
			this.shiftedShape = shape + Gamma.LANCZOS_G + 0.5;
			final double aux = FastMath.E / (2.0 * FastMath.PI * shiftedShape);
			this.densityPrefactor2 = shape * FastMath.sqrt(aux) / Gamma.lanczos(shape);
			this.densityPrefactor1 = this.densityPrefactor2 / scale * FastMath.pow(shiftedShape, -shape) *
					FastMath.exp(shape + Gamma.LANCZOS_G);
			this.minY = shape + Gamma.LANCZOS_G - FastMath.log(Double.MAX_VALUE);
			this.maxLogY = FastMath.log(Double.MAX_VALUE) / (shape - 1.0);
			uninitialised = false;
		}
	}

	/**
	 * {@inheritDoc}
	 *
	 * The implementation of this method is based on:
	 * <ul>
	 * <li>
	 * <a href="http://mathworld.wolfram.com/Chi-SquaredDistribution.html"> Chi-Squared Distribution</a>, equation
	 * (9).</li>
	 * <li>Casella, G., & Berger, R. (1990). <i>Statistical Inference</i>. Belmont, CA: Duxbury Press.</li>
	 * </ul>
	 */
	@Override
	public double cumulativeProbability(double x)
	{
		double ret;

		if (x <= 0)
		{
			ret = 0;
		}
		else
		{
			ret = Gamma.regularizedGammaP(shape, x / scale);
		}

		return ret;
	}

	/** {@inheritDoc} */
	@Override
	protected double getSolverAbsoluteAccuracy()
	{
		return solverAbsoluteAccuracy;
	}

	/**
	 * {@inheritDoc}
	 *
	 * For shape parameter {@code alpha} and scale parameter {@code beta}, the
	 * mean is {@code alpha * beta}.
	 */
	@Override
	public double getNumericalMean()
	{
		return shape * scale;
	}

	/**
	 * {@inheritDoc}
	 *
	 * For shape parameter {@code alpha} and scale parameter {@code beta}, the
	 * variance is {@code alpha * beta^2}.
	 *
	 * @return {@inheritDoc}
	 */
	@Override
	public double getNumericalVariance()
	{
		return shape * scale * scale;
	}

	/**
	 * {@inheritDoc}
	 *
	 * The lower bound of the support is always 0 no matter the parameters.
	 *
	 * @return lower bound of the support (always 0)
	 */
	@Override
	public double getSupportLowerBound()
	{
		return 0;
	}

	/**
	 * {@inheritDoc}
	 *
	 * The upper bound of the support is always positive infinity
	 * no matter the parameters.
	 *
	 * @return upper bound of the support (always Double.POSITIVE_INFINITY)
	 */
	@Override
	public double getSupportUpperBound()
	{
		return Double.POSITIVE_INFINITY;
	}

	/** {@inheritDoc} */
	@Override
	public boolean isSupportLowerBoundInclusive()
	{
		return true;
	}

	/** {@inheritDoc} */
	@Override
	public boolean isSupportUpperBoundInclusive()
	{
		return false;
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
	 * <p>
	 * This implementation uses the following algorithms:
	 * </p>
	 *
	 * <p>
	 * For 0 < shape < 1: <br/>
	 * Ahrens, J. H. and Dieter, U., <i>Computer methods for sampling from gamma, beta, Poisson and binomial
	 * distributions.</i> Computing, 12, 223-246, 1974.
	 * </p>
	 *
	 * <p>
	 * For shape >= 1: <br/>
	 * Marsaglia and Tsang, <i>A Simple Method for Generating Gamma Variables.</i> ACM Transactions on Mathematical
	 * Software, Volume 26 Issue 3, September, 2000.
	 * </p>
	 *
	 * @return random value sampled from the Gamma(shape, scale) distribution
	 */
	@Override
	public double sample()
	{
		if (shape < 1)
		{
			// [1]: p. 228, Algorithm GS

			while (true)
			{
				// Step 1:
				final double u = random.nextDouble();
				final double bGS = 1 + shape / FastMath.E;
				final double p = bGS * u;

				if (p <= 1)
				{
					// Step 2:

					final double x = FastMath.pow(p, 1 / shape);
					final double u2 = random.nextDouble();

					if (u2 > FastMath.exp(-x))
					{
						// Reject
						continue;
					}
					else
					{
						return scale * x;
					}
				}
				else
				{
					// Step 3:

					final double x = -1 * FastMath.log((bGS - p) / shape);
					final double u2 = random.nextDouble();

					if (u2 > FastMath.pow(x, shape - 1))
					{
						// Reject
						continue;
					}
					else
					{
						return scale * x;
					}
				}
			}
		}

		// Now shape >= 1

		final double d = shape - 0.333333333333333333;
		final double c = 1 / (3 * FastMath.sqrt(d));

		while (true)
		{
			final double x = random.nextGaussian();
			final double v = (1 + c * x) * (1 + c * x) * (1 + c * x);

			if (v <= 0)
			{
				continue;
			}

			final double x2 = x * x;
			final double u = random.nextDouble();

			// Squeeze
			if (u < 1 - 0.0331 * x2 * x2)
			{
				return scale * d * v;
			}

			if (FastMath.log(u) < 0.5 * x2 + d * (1 - v + FastMath.log(v)))
			{
				return scale * d * v;
			}
		}
	}
}
