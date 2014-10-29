package gdsc.smlm.utils;

import java.util.Arrays;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.util.FastMath;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Contains methods to find the noise in the provided image data.
 * <p>
 * Certain noise estimation routines have been copied from the estimation routines of ND-Safir (N-dimensional noise
 * reduction software): <br/>
 * http://raweb.inria.fr/rapportsactivite/RA2011/serpico/uid21.html
 */
public class NoiseEstimator
{
	public enum Method
	{
		/**
		 * Use all pixels
		 */
		ALL_PIXELS("All pixels"),
		/**
		 * Use a range around the lowest pixel in the image
		 */
		LOWEST_PIXELS("Lowest pixels"),
		/**
		 * Use the psuedo-residuals and calculate the least median of squares
		 */
		RESIDUALS_LEAST_MEDIAN_OF_SQUARES("Residuals least-median-of-squares"),
		/**
		 * Use the psuedo-residuals and calculate the least trimmed of squares
		 */
		RESIDUALS_LEAST_TRIMMED_OF_SQUARES("Residuals least-trimmed-of-squares"),
		/**
		 * Use the psuedo-residuals and calculate the least mean of squares
		 */
		RESIDUALS_LEAST_MEAN_OF_SQUARES("Residuals least-mean-of-squares"),
		/**
		 * Use the psuedo-residuals ignoring image border and calculate the least median of squares
		 */
		QUICK_RESIDUALS_LEAST_MEDIAN_OF_SQUARES("Quick residuals least-median-of-squares"),
		/**
		 * Use the psuedo-residuals ignoring image border and calculate the least trimmed of squares
		 */
		QUICK_RESIDUALS_LEAST_TRIMMED_OF_SQUARES("Quick residuals least-trimmed-of-squares"),
		/**
		 * Use the psuedo-residuals ignoring image border and calculate the least mean of squares
		 */
		QUICK_RESIDUALS_LEAST_MEAN_OF_SQUARES("Quick residuals least-mean-of-squares");

		private String name;

		private Method(String name)
		{
			this.name = name;
		}

		@Override
		public String toString()
		{
			return name;
		}
	}

	private float[] data;
	private float[] residuals = null;
	private float[] quickResiduals = null;
	private int maxx;
	private int maxy;

	private int range = 6;
	/**
	 * Set this to true if multiple calls will be made to {@link #getNoise(Method)} using methods that modify the
	 * residuals (LeastMedian or LeastTrimmed). If false these methods destroy the residuals which then have to be
	 * recomputed.
	 */
	public boolean preserveResiduals = false;

	/**
	 * @param data
	 * @param maxx
	 * @param maxy
	 */
	public NoiseEstimator(float[] data, int maxx, int maxy)
	{
		if (maxx < 1 || maxy < 1)
			throw new IllegalArgumentException("X/Y dimensions must be larger than 0");
		if (data == null || data.length < maxx * maxy)
			throw new IllegalArgumentException("Data must be at least as large as the given dimensions");
		this.data = data;
		this.maxx = maxx;
		this.maxy = maxy;
	}

	/**
	 * Estimates the noise using random pixels from the image.
	 * 
	 * @param method
	 */
	public double getNoise(Method method)
	{
		Estimator ne;

		switch (method)
		{
			case QUICK_RESIDUALS_LEAST_TRIMMED_OF_SQUARES:
				ne = new ResidualsLeastTrimmedSquareEstimator(true);
				break;

			case QUICK_RESIDUALS_LEAST_MEDIAN_OF_SQUARES:
				ne = new ResidualsLeastMedianSquareEstimator(true);
				break;

			case QUICK_RESIDUALS_LEAST_MEAN_OF_SQUARES:
				ne = new ResidualsLeastMeanSquareEstimator(true);
				break;

			case RESIDUALS_LEAST_TRIMMED_OF_SQUARES:
				ne = new ResidualsLeastTrimmedSquareEstimator(false);
				break;

			case RESIDUALS_LEAST_MEDIAN_OF_SQUARES:
				ne = new ResidualsLeastMedianSquareEstimator(false);
				break;

			case RESIDUALS_LEAST_MEAN_OF_SQUARES:
				ne = new ResidualsLeastMeanSquareEstimator(false);
				break;

			case LOWEST_PIXELS:
				ne = new MinEstimator(range);
				break;

			default:
				ne = new AllEstimator();
		}

		return ne.getNoise();
	}

	/**
	 * Provide the base implementation for all noise estimators
	 */
	private abstract class Estimator
	{
		abstract double getNoise();
	}

	/**
	 * Estimate the noise using standard deviation of all pixels in an image
	 */
	private class AllEstimator extends Estimator
	{
		@Override
		double getNoise()
		{
			SummaryStatistics stats = new SummaryStatistics();
			for (int i = maxx * maxy; i-- > 0;)
				stats.addValue(data[i]);
			return stats.getStandardDeviation();
		}
	}

	/**
	 * Estimate noise using region around lowest pixel in image
	 */
	private class MinEstimator extends Estimator
	{
		final int range;

		public MinEstimator(int range)
		{
			this.range = range;
		}

		@Override
		double getNoise()
		{
			// Get the image minimum
			float min = Float.POSITIVE_INFINITY;
			int index = 0;
			for (int i = maxx * maxy; i-- > 0;)
			{
				if (min > data[i])
				{
					min = data[i];
					index = i;
				}
			}

			int x = index % maxx;
			int y = index / maxx;
			int ys = FastMath.max(y - range, 0);
			int ye = FastMath.min(y + range, maxy - 1);
			int xs = FastMath.max(x - range, 0);
			int xe = FastMath.min(x + range, maxx - 1);

			SummaryStatistics stats = new SummaryStatistics();
			for (int y2 = ys; y2 <= ye; y2++)
			{
				for (int x2 = xs, i = ys * maxx + xs; x2 <= xe; x2++, i++)
				{
					stats.addValue(data[i]);
				}
			}
			return stats.getStandardDeviation();
		}
	}

	private class ResidualsLeastMedianSquareEstimator extends Estimator
	{
		public boolean quick = false;

		public ResidualsLeastMedianSquareEstimator(boolean quick)
		{
			this.quick = quick;
		}

		@Override
		double getNoise()
		{
			float[] buf = (quick) ? getQuickPseudoResiduals() : getPseudoResiduals();
			int n = buf.length;
			if (n < 2)
				return 0;
			if (preserveResiduals)
				buf = Arrays.copyOf(buf, buf.length);
			Arrays.sort(buf);
			float med_i = buf[(int) (.5 * (float) n)];
			for (int j = 0; j < n; j++)
				buf[j] = Math.abs(buf[j] - med_i);
			Arrays.sort(buf);
			double sig = 1.4828 * buf[(int) (.5 * (float) n)];
			if (!preserveResiduals)
			{
				// Residuals have been destroyed
				if (quick)
					quickResiduals = null;
				else
					residuals = null;
			}
			return Math.abs(sig);
		}
	}

	private class ResidualsLeastTrimmedSquareEstimator extends Estimator
	{
		public boolean quick = false;

		public ResidualsLeastTrimmedSquareEstimator(boolean quick)
		{
			this.quick = quick;
		}

		@Override
		double getNoise()
		{
			float[] buf = (quick) ? getQuickPseudoResiduals() : getPseudoResiduals();
			int n = buf.length;
			if (n < 2)
				return 0;
			if (preserveResiduals)
				buf = Arrays.copyOf(buf, buf.length);
			for (int k = 0; k < n; k++)
				buf[k] = buf[k] * buf[k];
			Arrays.sort(buf);
			double a = 0;
			for (int j = 0; j < (int) (.5 * n); j++)
				a += buf[j];
			double sig = 2.6477 * Math.sqrt(a / (int) (.5 * n));
			if (!preserveResiduals)
			{
				// Residuals have been destroyed
				if (quick)
					quickResiduals = null;
				else
					residuals = null;
			}
			return Math.abs(sig);
		}
	}

	private class ResidualsLeastMeanSquareEstimator extends Estimator
	{
		public boolean quick = false;

		public ResidualsLeastMeanSquareEstimator(boolean quick)
		{
			this.quick = quick;
		}

		@Override
		double getNoise()
		{
			float[] buf = (quick) ? getQuickPseudoResiduals() : getPseudoResiduals();
			if (buf.length < 2)
				return 0;
			double a = 0, b = 0;
			for (int i = 0; i < buf.length; i++)
			{
				a += buf[i];
				b += buf[i] * buf[i];
			}
			a /= (double) buf.length;
			b /= (double) buf.length;
			b -= a * a;
			return (b > 0) ? Math.sqrt(b) : 0;
		}
	}

	/**
	 * Compute the pseudo-residuals of the input data.
	 * <p>
	 * The pseudo residual \f$ R(x,y) \f$ of the image \f$ I(x,y) \f$ are defined by \f$ R(x,y) = 4 * I(x,y) - (I(x+1,y)
	 * + I(x-1,y) + I(x,y+1) + I(x,y-1))\f$ and normalized so that \f$ \mathbb{E}[R(x,y)^2] = \mathbb{E}[I(x,y)^2].
	 * 
	 * @return The pseudo residuals
	 */
	public float[] getPseudoResiduals()
	{
		if (residuals == null)
		{
			residuals = new float[maxx * maxy];

			for (int y = 0, index = 0; y < maxy; y++)
				for (int x = 0; x < maxx; x++, index++)
				{
					double t2 = 0;
					if (x == 0)
						t2 += data[index + 1];
					else
						t2 += data[index - 1];
					if (x == maxx - 1)
						t2 += data[index - 1];
					else
						t2 += data[index + 1];
					if (y == 0)
						t2 += data[index + maxx];
					else
						t2 += data[index - maxx];
					if (y == maxy - 1)
						t2 += data[index - maxx];
					else
						t2 += data[index + maxx];

					// 0.223606798 = 1 / sqrt(20)
					residuals[index] = (float) (0.223606798 * (4. * (double) data[index] - t2));
				}
		}

		return residuals;
	}

	/**
	 * Compute the pseudo-residuals of the input data. Ignore the image border so output will be size = (maxx - 2) *
	 * (maxy - 2)
	 * 
	 * @return The pseudo residuals
	 */
	public float[] getQuickPseudoResiduals()
	{
		if (quickResiduals == null)
		{
			if (maxx < 3 || maxy < 3)
			{
				quickResiduals = new float[0];
				return quickResiduals;
			}

			quickResiduals = new float[(maxx - 2) * (maxy - 2)];

			for (int y = 1, i = 0; y < maxy - 1; y++)
				for (int x = 1, index = y * maxx + 1; x < maxx - 1; x++, index++, i++)
				{
					double t2 = data[index - 1] + data[index + 1] + data[index - maxx] + data[index + maxx];
					// 0.223606798 = 1 / sqrt(20)
					quickResiduals[i] = (float) (0.223606798 * (4. * (double) data[index] - t2));
				}
		}

		return quickResiduals;
	}

	/**
	 * @param range
	 *            the range for the search around the local minimum. Must be at least 1.
	 */
	public void setRange(int range)
	{
		if (range < 1)
			range = 1;
		this.range = range;
	}

	/**
	 * @return the range
	 */
	public int getRange()
	{
		return range;
	}
}
