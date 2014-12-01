package gdsc.smlm.function;

import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.util.FastMath;

/**
 * Calculate the value of the Skew Normal distribution
 * <p>
 * See http://en.wikipedia.org/wiki/Skew_normal_distribution
 */
public class SkewNormalFunction
{
	protected double amplitude, location, scale, alpha;

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
	 * @return The mean of the function
	 */
	public double getMean()
	{
		return location + scale * getSigma() * Math.sqrt(2.0 / Math.PI);
	}

	/**
	 * @return The variance of the function
	 */
	public double getVariance()
	{
		final double sigma = getSigma();
		return scale * scale * (1.0 - (2.0 * sigma * sigma / Math.PI));
	}

	public double getSkewness()
	{
		final double sigma = getSigma();
		return ((4.0 - Math.PI) / 2.0) *
				(Math.pow((sigma * Math.sqrt(2.0 / Math.PI)), 3) / Math.pow(1 - 2 * sigma * sigma / Math.PI, 1.5));
	}

	private double getSigma()
	{
		return alpha / Math.sqrt(1 + alpha * alpha);
	}

	/**
	 * Evaluates the skewed Gaussian at the given point * @param x
	 * 
	 * @return
	 */
	public double evaluate(double x)
	{
		return evaluate(x, amplitude, location, scale, alpha);
	}

	/**
	 * Evaluates the skewed Gaussian at the given point
	 * 
	 * @param x
	 * @param parameters
	 *            [amplitude, location, scale, alpha]
	 * @return
	 */
	public static double evaluate(double x, double[] parameters)
	{
		return evaluate(x, parameters[0], parameters[1], parameters[2], parameters[3]);
	}

	/**
	 * Evaluates the skewed Gaussian at the given point
	 * 
	 * @param x
	 * @param amplitude
	 * @param location
	 * @param scale
	 * @param alpha
	 * @return
	 */
	public static double evaluate(double x, double amplitude, double location, double scale, double alpha)
	{
		//return (2 * amplitude / scale) * normal((x - location) / scale) * cumul(alpha * (x - location) / scale);

		// Do not normalise the area under the graph to 1. This allows the amplitude to be correctly modelled.
		return (2 * amplitude) * normal((x - location) / scale) * cumul(alpha * (x - location) / scale);
	}

	/**
	 * Probability density function of the Gaussian
	 * 
	 * @param x
	 * @return
	 */
	private static double normal(double x)
	{
		// 1/sqrt(2*pi) = 0.39894228 
		//return 0.39894228 * FastMath.exp(-0.5 * x*x);

		// Do not normalise the area under the graph to 1. This allows the amplitude to be correctly modelled.
		return FastMath.exp(-0.5 * x * x);
	}

	/**
	 * Cumulative distribution function of the gaussian
	 * 
	 * @param x
	 * @return
	 */
	private static double cumul(double x)
	{
		// 1/sqrt(2) = 0.707106781 
		return 0.5 * (1 + Erf.erf(x * 0.707106781));
	}
}