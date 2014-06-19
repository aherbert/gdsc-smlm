package gdsc.smlm.fitting.nonlinear.stop;

import gdsc.smlm.fitting.function.Gaussian2DFunction;
import gdsc.smlm.fitting.function.GaussianFunction;
import gdsc.smlm.fitting.utils.FloatEquality;

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
 * Stop when successive iterations with a reduced error move the fitted coordinates by less than a specified
 * distance.
 * <p>
 * The criteria also ensure that amplitude, coordinates and peak-widths are held positive, otherwise fitting is stopped.
 */
public class ParameterStoppingCriteria extends GaussianStoppingCriteria
{
	private int significantDigits = 3;
	private float angleLimit = 1e-3f;

	private FloatEquality eq;

	/**
	 * @param func
	 *            The Gaussian function
	 */
	public ParameterStoppingCriteria(GaussianFunction func)
	{
		super(func);
		eq = new FloatEquality(significantDigits, 1e-16f);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.stop.GaussianStoppingCriteria#logParameters(double, double, float[])
	 */
	@Override
	protected StringBuffer logParameters(double oldError, double newError, float[] a)
	{
		if (log != null)
		{
			StringBuffer sb = new StringBuffer();
			sb.append("iter = ").append(getIteration() + 1).append(", error = ").append(oldError).append(" -> ")
					.append(newError);
			if (newError <= oldError)
			{
				if (func.evaluatesBackground())
				{
					sb.append(", Back=[");
					sb.append(FloatEquality.relativeError(bestA[0], a[0]));
					sb.append("]");
				}

				for (int i = 0; i < peaks; i++)
				{
					sb.append(", Peak").append(i + 1).append("=[");
					sb.append(FloatEquality.relativeError(bestA[i * 6 + Gaussian2DFunction.AMPLITUDE], a[i * 6 +
							Gaussian2DFunction.AMPLITUDE]));
					sb.append(",");

					if (func.evaluatesAngle())
					{
						float x = bestA[i * 6 + Gaussian2DFunction.ANGLE];
						float y = a[i * 6 + Gaussian2DFunction.ANGLE];
						sb.append(relativeAngle(x, y));
					}
					else
						sb.append(0);

					for (int j = 0, k = i * 6 + Gaussian2DFunction.X_POSITION; j < 2 * dimensions; j++, k++)
					{
						sb.append(",");
						sb.append(FloatEquality.relativeError(bestA[k], a[k]));
					}
					sb.append("]");
				}
			}
			return sb;
		}
		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.stop.GaussianStoppingCriteria#noCoordinateChange(float[])
	 */
	@Override
	protected boolean noCoordinateChange(float[] a)
	{
		// Old code does not correctly compute difference in angles. This is ignored for now.
		//return eq.almostEqualComplement(bestA, a);

		if (func.evaluatesBackground())
		{
			if (!eq.almostEqualComplement(bestA[Gaussian2DFunction.BACKGROUND], a[Gaussian2DFunction.BACKGROUND]))
				return false;
		}

		for (int i = 0; i < peaks; i++)
		{
			if (!eq.almostEqualComplement(bestA[i * 6 + Gaussian2DFunction.AMPLITUDE], a[i * 6 +
					Gaussian2DFunction.AMPLITUDE]))
				return false;

			// Calculate the smallest angle between the two angles. This should be in the range 0 - 90 degrees.
			// Use this to compare if the angle has changed significantly relative to the maximum it could change.
			if (func.evaluatesAngle())
			{
				float x = bestA[i * 6 + Gaussian2DFunction.ANGLE];
				float y = a[i * 6 + Gaussian2DFunction.ANGLE];
				if (relativeAngle(x, y) > angleLimit)
					return false;
			}

			for (int j = 0, k = i * 6 + Gaussian2DFunction.X_POSITION; j < 2 * dimensions; j++, k++)
			{
				if (!eq.almostEqualComplement(bestA[k], a[k]))
					return false;
			}
		}

		return true;
	}

	private float relativeAngle(float x, float y)
	{
		float angle = (float) Math.atan2(Math.sin(x - y), Math.cos(x - y));
		final float halfPi = (float) Math.PI / 2;
		return Math.abs(angle / halfPi);
	}

	/**
	 * Set the change in parameters that defines a negligible amount
	 * 
	 * @param significantDigits
	 *            the significantDigits to set
	 */
	public void setSignificantDigits(int significantDigits)
	{
		this.significantDigits = significantDigits;
		eq.setSignificantDigits(significantDigits);
		angleLimit = (float) (1.0 / Math.pow(10, significantDigits - 1));
	}

	/**
	 * @return the significantDigits
	 */
	public int getSignificantDigits()
	{
		return significantDigits;
	}
}
