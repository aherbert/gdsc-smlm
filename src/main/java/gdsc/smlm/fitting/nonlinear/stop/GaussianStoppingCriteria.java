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
package gdsc.smlm.fitting.nonlinear.stop;

import gdsc.smlm.fitting.nonlinear.StoppingCriteria;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;

import java.util.Arrays;

/**
 * Defines the stopping criteria for the {@link gdsc.smlm.fitting.nonlinear.NonLinearFit } class.
 * <p>
 * Stop when successive iterations with a reduced error move the fitted X,Y coordinates by less than a specified
 * distance (delta).
 * <p>
 * The criteria also ensure that signal, coordinates and peak-widths are held positive, otherwise fitting is stopped.
 */
public class GaussianStoppingCriteria extends StoppingCriteria
{
	private double delta = 0.01;

	protected int peaks;
	protected Gaussian2DFunction func;

	private double minimumSignal = Float.NEGATIVE_INFINITY;
	private double[] minimumPosition = null;
	private double[] maximumPosition = null;
	private double[] minimumSD = null;
	private double[] maximumSD = null;

	/**
	 * @param func
	 *            The Gaussian function
	 */
	public GaussianStoppingCriteria(Gaussian2DFunction func)
	{
		this.func = func;
		this.peaks = func.getNPeaks();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.stoppingCriteria#evaluate(double, double, double[])
	 */
	@Override
	public void evaluate(double oldError, double newError, double[] a)
	{
		StringBuilder sb = logParameters(oldError, newError, a);

		if (newError > oldError)
		{
			// Fit is worse
			increment(a, false);
		}
		else
		{
			// Fit is improved - Check if the movement is negligible
			if (noCoordinateChange(a))
			{
				//				// Check if all params are within 2sf
				//				FloatEquality eq = new FloatEquality(3, 1e-10f);
				//				for (int i = 0; i < a.length; i++)
				//				{
				//					if (!eq.almostEqualComplement(bestA[i], a[i]))
				//						System.out.printf("Stopping when still moving: %f => %f (%g)\n%s\n%s\n",
				//								bestA[i], a[i], FloatEquality.relativeError(bestA[i], a[i]),
				//								Arrays.toString(bestA), Arrays.toString(a));
				//				}

				areAchieved = true;
				notSatisfied = false;
			}

			// Check the parameters are still valid
			if (invalidCoordinates(a))
			{
				notSatisfied = false;
				if (log != null)
				{
					sb.append(" Bad Coords: ").append(Arrays.toString(a));
				}
			}

			increment(a, true);
		}

		if (log != null)
		{
			sb.append(" Continue=").append(notSatisfied).append("\n");
			log.info(sb.toString());
		}
	}

	/**
	 * Creates a string representation of the peak parameters if logging
	 * 
	 * @param oldError
	 * @param newError
	 * @param a
	 *            The parameters
	 * @return The string
	 */
	protected StringBuilder logParameters(double oldError, double newError, double[] a)
	{
		if (log != null)
		{
			StringBuilder sb = new StringBuilder();
			sb.append("iter = ").append(getIteration() + 1).append(", error = ").append(oldError).append(" -> ")
					.append(newError);
			if (newError <= oldError)
			{
				for (int i = 0; i < peaks; i++)
				{
					sb.append(", Peak").append(i + 1).append("=[");
					for (int j = 0, k = i * Gaussian2DFunction.PARAMETERS_PER_PEAK +
							Gaussian2DFunction.X_POSITION; j < 2; j++, k++)
					{
						if (j > 0)
							sb.append(",");
						sb.append(a[k] - bestA[k]);
					}
					sb.append("]");
				}
			}
			return sb;
		}
		return null;
	}

	protected boolean noCoordinateChange(double[] a)
	{
		for (int i = 0; i < peaks; i++)
		{
			for (int j = 0, k = i * Gaussian2DFunction.PARAMETERS_PER_PEAK +
					Gaussian2DFunction.X_POSITION; j < 2; j++, k++)
			{
				// Check if the coordinates have moved less than the delta limit
				if (Math.abs(bestA[k] - a[k]) > delta)
					return false;
			}
		}
		return true;
	}

	private boolean invalidCoordinates(double[] a)
	{
		for (int i = 0; i < peaks; i++)
		{
			if (a[i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.SIGNAL] < minimumSignal)
				return true;

			if (isBelow(minimumPosition, a, i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_POSITION))
				return true;
			if (isAbove(maximumPosition, a, i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_POSITION))
				return true;

			if (func.evaluatesSD0())
			{
				if (isBelow(minimumSD, a, i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_SD))
					return true;
				if (isAbove(maximumSD, a, i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_SD))
					return true;
			}

			if (func.evaluatesSD1())
			{
				if (isBelow(minimumSD, a, i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.Y_SD))
					return true;
				if (isAbove(maximumSD, a, i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.Y_SD))
					return true;
			}
		}
		return false;
	}

	/**
	 * Check if any dimension is below the threshold
	 * 
	 * @param threshold
	 * @param params
	 * @param paramIndex
	 * @return
	 */
	private boolean isBelow(double[] threshold, double[] params, int paramIndex)
	{
		if (threshold != null)
		{
			for (int j = 0; j < 2; j++, paramIndex++)
			{
				if (params[paramIndex] < threshold[j])
					return true;
			}
		}
		return false;
	}

	/**
	 * Check if any dimension is above the threshold
	 * 
	 * @param threshold
	 * @param params
	 * @param paramIndex
	 * @return
	 */
	private boolean isAbove(double[] threshold, double[] params, int paramIndex)
	{
		if (threshold != null)
		{
			for (int j = 0; j < 2; j++, paramIndex++)
			{
				if (params[paramIndex] > threshold[j])
					return true;
			}
		}
		return false;
	}

	/**
	 * Set the change in error that defines a negligible amount
	 * 
	 * @param delta
	 *            the delta to set
	 */
	public void setDelta(double delta)
	{
		this.delta = delta;
	}

	/**
	 * @return the delta
	 */
	public double getDelta()
	{
		return delta;
	}

	/**
	 * @param minimumSignal
	 *            the minimum signal
	 */
	public void setMinimumSignal(double minimumSignal)
	{
		if (func.evaluatesSignal())
			this.minimumSignal = minimumSignal;
	}

	/**
	 * @return the minimum signal
	 */
	public double getMinimumSignal()
	{
		return minimumSignal;
	}

	/**
	 * @param minimumPosition
	 *            the minimum position for each dimension
	 */
	public void setMinimumPosition(double[] minimumPosition)
	{
		if (func.evaluatesPosition())
			this.minimumPosition = checkArray(minimumPosition);
	}

	/**
	 * @return the minimum position for each dimension
	 */
	public double[] getMinimumPosition()
	{
		return minimumPosition;
	}

	/**
	 * @param maximumPosition
	 *            the maximum position for each dimension
	 */
	public void setMaximumPosition(double[] maximumPosition)
	{
		if (func.evaluatesPosition())
			this.maximumPosition = checkArray(maximumPosition);
	}

	/**
	 * @return the maximum position for each dimension
	 */
	public double[] getMaximumPosition()
	{
		return maximumPosition;
	}

	/**
	 * @param minimumSD
	 *            the minimum SD for each dimension
	 */
	public void setMinimumSD(double[] minimumSD)
	{
		if (func.evaluatesSD0())
			this.minimumSD = checkArray(minimumSD);
	}

	/**
	 * @return the minimum SD for each dimension
	 */
	public double[] getMinimumSD()
	{
		return minimumSD;
	}

	/**
	 * @param maximumSD
	 *            the maximum SD for each dimension
	 */
	public void setMaximumSD(double[] maximumSD)
	{
		if (func.evaluatesSD0())
			this.maximumSD = checkArray(maximumSD);
	}

	/**
	 * @return the maximum SD for each dimension
	 */
	public double[] getMaximumSD()
	{
		return maximumSD;
	}

	private double[] checkArray(double[] array)
	{
		return (array == null || array.length != 2) ? null : array;
	}
}
