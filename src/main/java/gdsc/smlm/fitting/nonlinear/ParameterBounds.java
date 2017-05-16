package gdsc.smlm.fitting.nonlinear;

import java.util.Arrays;

import gdsc.smlm.function.GradientFunction;

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
 * 
 * This is an adaption of the C-code contained in the CcpNmr Analysis Program:
 *   CCPN website (http://www.ccpn.ac.uk/). 
 * The CCPN code was based on Numerical Recipes. 
 *---------------------------------------------------------------------------*/

/**
 * Allow restricting parameter steps.
 * <p>
 * Support bounded parameters using a hard-stop limit.
 * <p>
 * Support parameter clamping to prevent large parameter shifts. Optionally update the clamping when the search
 * direction changes.
 */
public class ParameterBounds
{
	private GradientFunction f;
	private int[] gradientIndices;
	private boolean isLower = false, isUpper = false;
	private double[] lower, upper;
	private boolean isClamped = false;
	private double[] clampInitial, clamp;
	private int[] dir;
	private boolean dynamicClamp = false;

	public ParameterBounds(GradientFunction f)
	{
		setGradientFunction(f);
	}

	/**
	 * Apply the step to the parameters a to generate new parameters. The step may be clamped and the new parameters may
	 * be bounded depending on settings.
	 * <p>
	 * Note that the step is the length defined by the gradient indices of the function. The length of a and newA are
	 * the full set of parameters to initialise the function. newA should be a clone of a such that the only difference
	 * is the step for each of the gradient indices.
	 *
	 * @param a
	 *            the a
	 * @param step
	 *            the step for each of the gradient indices
	 * @param newA
	 *            the new A
	 */
	public void applyBounds(double[] a, double[] step, double[] newA)
	{
		if (isClamped)
		{
			for (int j = gradientIndices.length; j-- > 0;)
			{
				if (clampInitial[j] == 0)
				{
					// Use the update parameter directly
					newA[gradientIndices[j]] = a[gradientIndices[j]] + step[j];
				}
				else
				{
					// This parameter is clamped
					newA[gradientIndices[j]] = a[gradientIndices[j]] + step[j] / clamp(step[j], j);
				}
			}
			applyBounds(newA);
		}
		else
		{
			for (int j = gradientIndices.length; j-- > 0;)
			{
				// Use the update parameter directly
				newA[gradientIndices[j]] = a[gradientIndices[j]] + step[j];
			}
			applyBounds(newA);
		}
	}

	/**
	 * Check the point falls within the configured bounds truncating if necessary.
	 *
	 * @param point
	 *            the point
	 */
	private void applyBounds(double[] point)
	{
		if (isUpper)
		{
			for (int j = gradientIndices.length; j-- > 0;)
				if (point[gradientIndices[j]] > upper[j])
				{
					point[gradientIndices[j]] = upper[j];
				}
		}
		if (isLower)
		{
			for (int j = gradientIndices.length; j-- > 0;)
				if (point[gradientIndices[j]] < lower[j])
				{
					point[gradientIndices[j]] = lower[j];
				}
		}
	}

	/**
	 * Produce the clamping value.
	 * <p>
	 * See Stetson PB (1987) DAOPHOT: A compute program for crowded-field stellar photometry. Publ Astrom Soc Pac
	 * 99:191-222. pp207-208
	 *
	 * @param u
	 *            the update parameter
	 * @param k
	 *            the parameter index
	 * @return the clamping value
	 */
	private double clamp(double u, int k)
	{
		if (u == 0)
			// Nothing to clamp
			return 1.0;

		double ck = clamp[k];
		if (dynamicClamp)
		{
			// If the sign has changed then reduce the clamp factor
			final int newDir = (u > 0) ? 1 : -1;

			// This addition overcomes the issue when the direction vector is new (i.e. zero filled)
			if (newDir + dir[k] == 0)
			{
				// Note: By reducing the size of the clamping factor we are restricting the movement
				ck *= 0.5;
			}

			// Note: We do not update the clamp[k] array yet as the move may be rejected. 
		}

		// Denominator for clamping function
		return 1 + (Math.abs(u) / ck);
	}

	/**
	 * Called when a step was accepted.
	 *
	 * @param a
	 *            the a
	 * @param newA
	 *            the new A
	 */
	public void accepted(double[] a, double[] newA)
	{
		if (isClamped && dynamicClamp)
		{
			// Get the direction and update the clamp parameter if the direction has changed
			for (int k = gradientIndices.length; k-- > 0;)
			{
				if (clamp[k] != 0)
				{
					final double u = newA[gradientIndices[k]] - a[gradientIndices[k]];
					if (u == 0)
						continue;
					final int newDir = (u > 0) ? 1 : -1;
					// This addition overcomes the issue when the direction vector is new (i.e. zero filled)
					if (newDir + dir[k] == 0)
					{
						// Note: By reducing the size of the clamping factor we are restricting the movement
						clamp[k] *= 0.5;
					}
					dir[k] = newDir;
				}
			}
		}
	}

	/**
	 * Initialise for clamping.
	 */
	public void initialise()
	{
		// Initialise for clamping
		if (isClamped)
		{
			// Prevent the clamping value being destroyed by dynamic updates
			if (dynamicClamp)
			{
				final int m = gradientIndices.length;
				clamp = Arrays.copyOf(clampInitial, m);
				// Reset the direction 
				if (dir == null || dir.length < m)
				{
					dir = new int[m];
				}
				else
				{
					for (int i = m; i-- > 0;)
						dir[i] = 0;
				}
			}
			else
			{
				clamp = clampInitial;
			}
		}
	}

	/**
	 * Set the bounds for each of the parameters. If a subset of the parameters are fitted then the bounds can be
	 * ignored for the fixed parameters.
	 * <p>
	 * The bounds can be used to set the expected range for a parameter.
	 *
	 * @param lowerBounds
	 *            the lower bounds
	 * @param upperBounds
	 *            the upper bounds
	 */
	public void setBounds(double[] lowerBounds, double[] upperBounds)
	{
		// Extract the bounds for the parameters we are fitting
		if (lowerBounds == null)
		{
			lower = null;
			isLower = false;
		}
		else
		{
			final int[] indices = f.gradientIndices();
			lower = new double[indices.length];
			for (int i = 0; i < indices.length; i++)
			{
				lower[i] = lowerBounds[indices[i]];
			}
			isLower = checkArray(lower, Double.NEGATIVE_INFINITY);
		}
		if (upperBounds == null)
		{
			upper = null;
			isUpper = false;
		}
		else
		{
			final int[] indices = f.gradientIndices();
			upper = new double[indices.length];
			for (int i = 0; i < indices.length; i++)
			{
				upper[i] = upperBounds[indices[i]];
			}
			isUpper = checkArray(upper, Double.POSITIVE_INFINITY);
		}
		// Check that the upper bound is above the lower bound
		if (isUpper && isLower)
		{
			for (int i = 0; i < lower.length; i++)
				if (lower[i] > upper[i])
					throw new IllegalArgumentException(
							"Lower bound is above upper bound: " + lower[i] + " > " + upper[i]);
		}
	}

	/**
	 * Check if the array contains anything other than value.
	 *
	 * @param array
	 *            the array
	 * @param value
	 *            the value
	 * @return True if the array has another value
	 */
	private static boolean checkArray(double[] array, double value)
	{
		for (int i = 0; i < array.length; i++)
			if (value != array[i])
				return true;
		return false;
	}

	/**
	 * Determine if the current solution (a) is at the the bounds
	 *
	 * @param a
	 *            the current parameters
	 * @return true, if the point is at the bounds
	 */
	public boolean atBounds(double[] a)
	{
		if (isUpper)
		{
			for (int i = 0; i < gradientIndices.length; i++)
				if (a[gradientIndices[i]] == upper[i])
				{
					return true;
				}
		}
		if (isLower)
		{
			for (int i = 0; i < gradientIndices.length; i++)
				if (a[gradientIndices[i]] == lower[i])
				{
					return true;
				}
		}
		return false;
	}

	/**
	 * Sets the parameter specific clamp values. This is the maximum permissible update to the parameter.
	 * <p>
	 * See Stetson PB (1987) DAOPHOT: A compute program for crowded-field stellar photometry. Publ Astrom Soc Pac
	 * 99:191-222.
	 * <p>
	 * Warning: If the function is changed then the clamp values may require updating. However setting a new function
	 * does not set the clamp values to null to allow caching when the clamp values are unchanged.
	 *
	 * @param clampValues
	 *            the new clamp values
	 */
	public void setClampValues(double[] clampValues)
	{
		// Extract the bounds for the parameters we are fitting
		if (clampValues == null)
		{
			clampInitial = null;
			isClamped = false;
		}
		else
		{
			clampInitial = new double[gradientIndices.length];
			for (int i = 0; i < gradientIndices.length; i++)
			{
				final double v = clampValues[gradientIndices[i]];
				if (Double.isNaN(v) || Double.isInfinite(v))
					continue;
				clampInitial[i] = Math.abs(v);
			}
			isClamped = checkArray(clampInitial, 0);
		}
	}

	/**
	 * Checks if is clamped.
	 *
	 * @return true, if is clamped
	 */
	public boolean isClamped()
	{
		return isClamped;
	}

	/**
	 * Checks if is dynamic clamping. The clamping factor will be reduced by a factor of 2 when the direction changes.
	 * <p>
	 * Note: Dynamic clamping only applies if {@link #isClamped()} is true, i.e. clamp values have been set.
	 *
	 * @return true, if is dynamic clamping
	 */
	public boolean isDynamicClamp()
	{
		return dynamicClamp;
	}

	/**
	 * Set to true to reduce the clamp factor by a factor of when the direction changes.
	 * <p>
	 * Note: Dynamic clamping only applies if {@link #isClamped()} is true, i.e. clamp values have been set.
	 *
	 * @param dynamicClamp
	 *            the new dynamic clamp
	 */
	public void setDynamicClamp(boolean dynamicClamp)
	{
		this.dynamicClamp = dynamicClamp;
	}

	/**
	 * Warning: If the function is changed then the clamp values may require updating. However setting a new function
	 * does not set the clamp values to null to allow caching when the clamp values are unchanged, e.g. evaluation of a
	 * different function in the same parameter space.
	 * <p>
	 * Setting a new function removes the current bounds.
	 */
	public void setGradientFunction(GradientFunction f)
	{
		if (f == null)
			throw new NullPointerException();

		// Do not do this to allow caching
		//setClampValues(null);

		setBounds(null, null);
		this.f = f;
		gradientIndices = f.gradientIndices();
	}
}
