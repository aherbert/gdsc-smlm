package gdsc.smlm.fitting;

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
 *---------------------------------------------------------------------------*/

import gdsc.smlm.function.gaussian.Gaussian2DFunction;

/**
 * Fits a 2-dimensional Gaussian function for the specified peak. Can optionally fit an elliptical Gaussian function.
 * <p>
 * Performs fitting using the configured algorithm.
 * <p>
 * This is based on the Gaussian2DFitter class but does not support estimating the width of each Gaussian.
 * Widths must be provided in the fit configuration. Settings from the fit configuration are cached and thus updates to
 * the configuration after construction may be ignored.
 */
public class FastGaussian2DFitter extends Gaussian2DFitter
{
	// Cache the fitting defaults
	private final boolean isWidth1Fitting, isAngleFitting;
	private final double angle, sx, sy;

	/**
	 * Constructor
	 * 
	 * @param fitConfiguration
	 * @throws IllegalArgumentException
	 *             If the configuration is missing information, e.g. initial widths
	 */
	public FastGaussian2DFitter(FitConfiguration fitConfiguration)
	{
		super(fitConfiguration);

		// Cache the estimate for the Gaussian
		if (fitConfiguration.getInitialPeakStdDev0() > 0)
			sx = fitConfiguration.getInitialPeakStdDev0();
		else
			throw new IllegalArgumentException("No initial width0 estimate");

		isWidth1Fitting = fitConfiguration.isWidth1Fitting();
		if (isWidth1Fitting)
		{
			if (fitConfiguration.getInitialPeakStdDev1() > 0)
				sy = fitConfiguration.getInitialPeakStdDev1();
			else
				throw new IllegalArgumentException("No initial width1 estimate");
		}
		else
		{
			sy = sx;
		}

		isAngleFitting = fitConfiguration.isAngleFitting();
		if (isAngleFitting)
		{
			if (fitConfiguration.getInitialAngle() >= -Math.PI && fitConfiguration.getInitialAngle() <= -Math.PI)
				angle = fitConfiguration.getInitialAngle();
			else
				throw new IllegalArgumentException("No initial angle estimate");
		}
		else
		{
			angle = 0;
		}
	}

	protected boolean checkParameters(final int maxx, final int maxy, final int npeaks, double[] params,
			final boolean[] amplitudeEstimate, final int ySize, final double[] y, final int paramsPerPeak,
			double background, double[] initialParams)
	{
		final int[] dim = new int[] { maxx, maxy };
		final int[] position = new int[2];
		for (int i = 0, j = 0; i < npeaks; i++, j += paramsPerPeak)
		{
			// Get the parameters
			double signal = params[j + Gaussian2DFunction.SIGNAL];
			double xpos = params[j + Gaussian2DFunction.X_POSITION];
			double ypos = params[j + Gaussian2DFunction.Y_POSITION];
			double sx = params[j + Gaussian2DFunction.X_SD];
			double sy = params[j + Gaussian2DFunction.Y_SD];
			double angle = params[j + Gaussian2DFunction.ANGLE];

			// ----
			// Check all input parameters and uses the default values if necessary
			// ----

			// Set-up for estimating peak width at half maximum 
			position[0] = (int) Math.round(xpos);
			position[1] = (int) Math.round(ypos);

			if (sx == 0)
			{
				sx = this.sx;
			}

			if (isWidth1Fitting)
			{
				if (sy == 0)
				{
					sy = this.sy;
				}
			}
			else
			{
				sy = sx;
			}

			// Guess the initial angle if input angle is out-of-bounds
			if (isAngleFitting)
			{
				if (angle == 0)
				{
					angle = this.angle;
				}
			}

			// If the position is on the integer grid then use a centre-of-mass approximation
			if (npeaks == 1 && xpos == position[0] && ypos == position[1])
			{
				// Estimate using centre of mass around peak index 
				// Use 2 * SD estimate to calculate the range around the index that should be considered.
				// SD = (sx+sy)/2 => Range = sx+sy
				final int range = Math.max(1, (int) Math.ceil(sx + sy));
				final double[] com = findCentreOfMass(y, dim, range, position);
				xpos = (double) com[0];
				ypos = (double) com[1];
			}

			// Convert amplitudes to signal
			if (amplitudeEstimate[i])
				signal *= 6.283185307 * sx * sy; // 2 * Math.PI * sx * sy

			// Set all the parameters
			params[j + Gaussian2DFunction.SIGNAL] = signal;
			params[j + Gaussian2DFunction.X_POSITION] = xpos;
			params[j + Gaussian2DFunction.Y_POSITION] = ypos;
			params[j + Gaussian2DFunction.X_SD] = sx;
			params[j + Gaussian2DFunction.Y_SD] = sy;
			params[j + Gaussian2DFunction.ANGLE] = angle;
		}

		return true;
	}
}