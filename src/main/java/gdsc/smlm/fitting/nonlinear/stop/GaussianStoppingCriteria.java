package gdsc.smlm.fitting.nonlinear.stop;

import java.util.Arrays;

import gdsc.smlm.fitting.function.Gaussian2DFunction;
import gdsc.smlm.fitting.function.GaussianFunction;
import gdsc.smlm.fitting.logging.Logger;
import gdsc.smlm.fitting.nonlinear.StoppingCriteria;

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
 * Defines the stopping criteria for the {@link gdsc.smlm.fitting.nonlinear.NonLinearFit } class.
 * <p>
 * Stop when successive iterations with a reduced error move the fitted X,Y coordinates by less than a specified
 * distance (delta).
 * <p>
 * The criteria also ensure that amplitude, coordinates and peak-widths are held positive, otherwise fitting is stopped.
 */
public class GaussianStoppingCriteria extends StoppingCriteria
{
	private double delta = 0.01;

	protected int peaks;
	protected int dimensions;
	protected GaussianFunction func;

	private float minimumAmplitude = Float.NEGATIVE_INFINITY;
	private float[] minimumPosition = null;
	private float[] maximumPosition = null;
	private float[] minimumWidth = null;
	private float[] maximumWidth = null;

	/**
	 * @param func
	 *            The Gaussian function
	 */
	public GaussianStoppingCriteria(GaussianFunction func)
	{
		this.func = func;
		this.peaks = func.getNPeaks();
		this.dimensions = func.getNDimensions();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.stoppingCriteria#evaluate(double, double, float[])
	 */
	@Override
	public void evaluate(double oldError, double newError, float[] a)
	{
		StringBuffer sb = logParameters(oldError, newError, a);

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
	protected StringBuffer logParameters(double oldError, double newError, float[] a)
	{
		if (log != null)
		{
			StringBuffer sb = new StringBuffer();
			sb.append("iter = ").append(getIteration() + 1).append(", error = ").append(oldError).append(" -> ")
					.append(newError);
			if (newError <= oldError)
			{
				for (int i = 0; i < peaks; i++)
				{
					sb.append(", Peak").append(i + 1).append("=[");
					for (int j = 0, k = i * 6 + Gaussian2DFunction.X_POSITION; j < dimensions; j++, k++)
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

	protected boolean noCoordinateChange(float[] a)
	{
		for (int i = 0; i < peaks; i++)
		{
			for (int j = 0, k = i * 6 + Gaussian2DFunction.X_POSITION; j < dimensions; j++, k++)
			{
				// Check if the coordinates have moved less than the delta limit
				if (Math.abs(bestA[k] - a[k]) > delta)
					return false;
			}
		}
		return true;
	}

	private boolean invalidCoordinates(float[] a)
	{
		for (int i = 0; i < peaks; i++)
		{
			if (a[i * 6 + Gaussian2DFunction.AMPLITUDE] < minimumAmplitude)
				return true;

			if (isBelow(minimumPosition, a, i * 6 + Gaussian2DFunction.X_POSITION))
				return true;
			if (isAbove(maximumPosition, a, i * 6 + Gaussian2DFunction.X_POSITION))
				return true;

			if (func.evaluatesWidth0())
			{
				if (isBelow(minimumWidth, a, i * 6 + Gaussian2DFunction.X_WIDTH))
					return true;
				if (isAbove(maximumWidth, a, i * 6 + Gaussian2DFunction.X_WIDTH))
					return true;
			}

			if (func.evaluatesWidth1())
			{
				if (isBelow(minimumWidth, a, i * 6 + Gaussian2DFunction.Y_WIDTH))
					return true;
				if (isAbove(maximumWidth, a, i * 6 + Gaussian2DFunction.Y_WIDTH))
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
	private boolean isBelow(float[] threshold, float[] params, int paramIndex)
	{
		if (threshold != null)
		{
			for (int j = 0; j < dimensions; j++, paramIndex++)
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
	private boolean isAbove(float[] threshold, float[] params, int paramIndex)
	{
		if (threshold != null)
		{
			for (int j = 0; j < dimensions; j++, paramIndex++)
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
	 * @param minimumAmplitude
	 *            the minimum amplitude
	 */
	public void setMinimumAmplitude(float minimumAmplitude)
	{
		if (func.evaluatesAmplitude())
			this.minimumAmplitude = minimumAmplitude;
	}

	/**
	 * @return the minimum amplitude
	 */
	public float getMinimumAmplitude()
	{
		return minimumAmplitude;
	}

	/**
	 * @param minimumPosition
	 *            the minimum position for each dimension
	 */
	public void setMinimumPosition(float[] minimumPosition)
	{
		if (func.evaluatesPosition())
			this.minimumPosition = checkArray(minimumPosition);
	}

	/**
	 * @return the minimum position for each dimension
	 */
	public float[] getMinimumPosition()
	{
		return minimumPosition;
	}

	/**
	 * @param maximumPosition
	 *            the maximum position for each dimension
	 */
	public void setMaximumPosition(float[] maximumPosition)
	{
		if (func.evaluatesPosition())
			this.maximumPosition = checkArray(maximumPosition);
	}

	/**
	 * @return the maximum position for each dimension
	 */
	public float[] getMaximumPosition()
	{
		return maximumPosition;
	}

	/**
	 * @param minimumWidth
	 *            the minimum width for each dimension
	 */
	public void setMinimumWidth(float[] minimumWidth)
	{
		if (func.evaluatesWidth0())
			this.minimumWidth = checkArray(minimumWidth);
	}

	/**
	 * @return the minimum width for each dimension
	 */
	public float[] getMinimumWidth()
	{
		return minimumWidth;
	}

	/**
	 * @param maximumWidth
	 *            the maximum width for each dimension
	 */
	public void setMaximumWidth(float[] maximumWidth)
	{
		if (func.evaluatesWidth0())
			this.maximumWidth = checkArray(maximumWidth);
	}

	/**
	 * @return the maximum width for each dimension
	 */
	public float[] getMaximumWidth()
	{
		return maximumWidth;
	}

	private float[] checkArray(float[] array)
	{
		return (array == null || array.length != func.getNDimensions()) ? null : array;
	}
}
