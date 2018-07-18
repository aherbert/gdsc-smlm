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

//import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.util.FastMath;

import uk.ac.sussex.gdsc.core.utils.Maths;

/**
 * Calculate the value of the Skew Normal distribution
 * <p>
 * @see <a href="http://en.wikipedia.org/wiki/Skew_normal_distribution">http://en.wikipedia.org/wiki/Skew_normal_distribution</a>
 */
public class SkewNormalFunction
{
	/** The amplitude. */
	protected double amplitude;

	/** The location. */
	protected double location;

	/** The scale. */
	protected double scale;

	/** The alpha. */
	protected double alpha;

	/**
	 * Create a function with the given parameters
	 *
	 * @param parameters
	 *            [amplitude, location, scale, alpha]
	 */
	public SkewNormalFunction(double[] parameters)
	{
		setParameters(parameters);
	}

	/**
	 * Set the parameters
	 *
	 * @param parameters
	 *            [amplitude, location, scale, alpha]
	 */
	public void setParameters(double[] parameters)
	{
		amplitude = parameters[0];
		location = parameters[1];
		scale = parameters[2];
		alpha = parameters[3];
	}

	/**
	 * Gets the mean.
	 *
	 * @return The mean of the function
	 */
	public double getMean()
	{
		return location + scale * getSigma() * Math.sqrt(2.0 / Math.PI);
	}

	/**
	 * Gets the variance.
	 *
	 * @return The variance of the function
	 */
	public double getVariance()
	{
		final double sigma = getSigma();
		return scale * scale * (1.0 - (2.0 * sigma * sigma / Math.PI));
	}

	/**
	 * Gets the skewness.
	 *
	 * @return the skewness
	 */
	public double getSkewness()
	{
		final double sigma = getSigma();
		return ((4.0 - Math.PI) / 2.0) *
				(Maths.pow3((sigma * Math.sqrt(2.0 / Math.PI))) / Math.pow(1 - 2 * sigma * sigma / Math.PI, 1.5));
	}

	/**
	 * Gets the sigma.
	 *
	 * @return the sigma
	 */
	private double getSigma()
	{
		return alpha / Math.sqrt(1 + alpha * alpha);
	}

	/**
	 * Evaluates the skewed Gaussian at the given point * @param x.
	 *
	 * @param x
	 *            the x
	 * @return the value
	 */
	public double evaluate(double x)
	{
		return evaluate(x, amplitude, location, scale, alpha);
	}

	/**
	 * Evaluates the skewed Gaussian at the given point.
	 *
	 * @param x
	 *            the x
	 * @param parameters
	 *            [amplitude, location, scale, alpha]
	 * @return the value
	 */
	public static double evaluate(double x, double[] parameters)
	{
		return evaluate(x, parameters[0], parameters[1], parameters[2], parameters[3]);
	}

	/**
	 * Evaluates the skewed Gaussian at the given point.
	 *
	 * @param x
	 *            the x
	 * @param amplitude
	 *            the amplitude
	 * @param location
	 *            the location
	 * @param scale
	 *            the scale
	 * @param alpha
	 *            the alpha
	 * @return the value
	 */
	public static double evaluate(double x, double amplitude, double location, double scale, double alpha)
	{
		//return (2 * amplitude / scale) * normal((x - location) / scale) * cumul(alpha * (x - location) / scale);

		// Do not normalise the area under the graph to 1. This allows the amplitude to be correctly modelled.
		return (2 * amplitude) * normal((x - location) / scale) * cumul(alpha * (x - location) / scale);
	}

	/**
	 * Probability density function of the Gaussian.
	 *
	 * @param x
	 *            the x
	 * @return the probability density
	 */
	private static double normal(double x)
	{
		// 1/sqrt(2*pi) = 0.39894228
		//return 0.39894228 * FastMath.exp(-0.5 * x*x);

		// Do not normalise the area under the graph to 1. This allows the amplitude to be correctly modelled.
		return FastMath.exp(-0.5 * x * x);
	}

	/**
	 * Cumulative distribution function of the gaussian.
	 *
	 * @param x
	 *            the x
	 * @return the cumulative
	 */
	private static double cumul(double x)
	{
		// 1/sqrt(2) = 0.707106781
		return 0.5 * (1 + Erf.erf(x * 0.707106781));
	}
}
