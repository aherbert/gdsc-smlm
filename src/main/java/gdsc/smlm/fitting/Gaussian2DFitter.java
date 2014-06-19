package gdsc.smlm.fitting;

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
 * 
 * This is an adaption of the C-code contained in the CcpNmr Analysis Program:
 * CCPN website (http://www.ccpn.ac.uk/)
 *---------------------------------------------------------------------------*/

import gdsc.smlm.fitting.function.Gaussian2DFunction;

import java.util.Arrays;

/**
 * Fits a 2-dimensional Gaussian function for the specified peak. Can optionally fit an elliptical Gaussian function.
 * <p>
 * Performs fitting using the Levenberg-Marquardt algorithm.
 */
public class Gaussian2DFitter
{
	private FitConfiguration fitConfiguration;
	private FunctionSolver solver;

	// The last successful fit. Used to compute the residuals.
	private float[] y_fit = null;
	// Allow calculation of residuals to be turned off (overwrite constructor fit configuration)
	private boolean computeResiduals = true;

	/**
	 * Constructor
	 * 
	 * @param fitConfiguration
	 */
	public Gaussian2DFitter(FitConfiguration fitConfiguration)
	{
		if (fitConfiguration == null)
		{
			throw new NullPointerException("No fit configuration");
		}
		this.fitConfiguration = fitConfiguration;
		computeResiduals = fitConfiguration.isComputeResiduals();
	}

	static float half_max_position(float[] data, int index, int[] point, int[] dim, int dimension, int[] cumul_region,
			int dirn, float background)
	{
		int i, i_start, i_end, i_step;
		float v = data[index];
		float v_half = 0.5f * (v + background);
		float v_prev = v, v_this;
		int jump;

		if (dirn == 1)
		{
			i_start = point[dimension] + 1;
			i_end = dim[dimension];
			i_step = 1;
		}
		else
		{
			i_start = point[dimension] - 1;
			i_end = -1;
			i_step = -1;
		}

		jump = i_step * cumul_region[dimension];

		for (i = i_start; i != i_end; i += i_step)
		{
			index += jump;
			v_this = data[index];

			if (v_this < v_half)
				return i - i_step * (v_half - v_this) / (v_prev - v_this);

			v_prev = v_this;
		}

		// Not reached the half-max point. Return the dimension limit.
		if (dirn == 1)
			return dim[dimension];
		else
			return 0f;
	}

	public static float half_max_linewidth(float[] data, int index, int[] point, int[] dim, int dimension,
			int[] cumul_region, float background)
	{
		float linewidth, a, b;

		a = half_max_position(data, index, point, dim, dimension, cumul_region, 1, background);
		b = half_max_position(data, index, point, dim, dimension, cumul_region, -1, background);

		linewidth = a - b;

		return linewidth;
	}

	/**
	 * Accepts a single array containing 2-dimensional data and a list of the
	 * peaks to fit. Data should be packed in descending dimension order,
	 * e.g. Y,X : Index for [y,z] = MaxX*y + x.
	 * <p>
	 * Performs fitting using the specified method with a Levenberg-Marquardt algorithm.
	 * <p>
	 * Adapted from the CCPN fit_peaks routine for Python.
	 * 
	 * @param data
	 *            The data to fit
	 * @param maxx
	 *            The data size in the x dimension
	 * @param maxy
	 *            The data size in the y dimension
	 * @param peaks
	 *            The index of the peaks
	 * @return The fit result
	 */
	public FitResult fit(final float[] data, final int maxx, final int maxy, final int[] peaks)
	{
		return fit(data, maxx, maxy, peaks, null);
	}

	/**
	 * Accepts a single array containing 2-dimensional data and a list of the
	 * peaks to fit. Data should be packed in descending dimension order,
	 * e.g. Y,X : Index for [y,z] = MaxX*y + x.
	 * <p>
	 * Performs fitting using the specified method with a Levenberg-Marquardt algorithm.
	 * 
	 * @param data
	 *            The data to fit
	 * @param maxx
	 *            The data size in the x dimension
	 * @param maxy
	 *            The data size in the y dimension
	 * @param peaks
	 *            The index of the peaks (must be within the data bounds if the heights are null)
	 * @param heights
	 *            An initial estimate of the peak heights (can be null)
	 * @return The fit result
	 */
	public FitResult fit(final float[] data, final int maxx, final int maxy, final int[] peaks, float[] heights)
	{
		int npeaks = peaks.length;

		final int paramsPerPeak = 6;

		float[] params = new float[1 + paramsPerPeak * npeaks];

		// Get peak heights
		if (heights == null || heights.length != peaks.length)
		{
			heights = new float[peaks.length];
			for (int i = 0; i < peaks.length; i++)
				heights[i] = data[peaks[i]];
		}

		float background = getBackground(data, maxx, maxy, heights);

		// Set the initial parameters
		params[0] = background;
		for (int i = 0, j = 0; i < peaks.length; i++, j += paramsPerPeak)
		{
			int index = peaks[i];
			params[j + Gaussian2DFunction.AMPLITUDE] = heights[i] - background;
			params[j + Gaussian2DFunction.X_POSITION] = index % maxx;
			params[j + Gaussian2DFunction.Y_POSITION] = index / maxx;
		}

		return fit(data, maxx, maxy, npeaks, params);
	}

	/**
	 * Guess the background from the data given the estimated peak heights.
	 * <p>
	 * For a single peak the method assumes the peak is in the centre. In this case the average edge value of the data
	 * is used but if the calculated background is above the peak height the minimum value in the data is used.
	 * <p>
	 * If multiple peaks heights are provided then always use the minimum value in the data since it cannot be assumed
	 * that all peaks are away from the edge of the data.
	 * 
	 * @param data
	 * @param maxx
	 * @param maxy
	 * @param heights
	 * @return The background estimate
	 */
	public static float getBackground(final float[] data, final int maxx, final int maxy, final float[] heights)
	{
		// TODO - What is the best method for setting the background?
		// 1. Min in data
		// 2. Average of border?
		// 3. Evaluate total volume under the initial Gaussian params and subtract that from the sum of the image?
		//    (note initial gaussian requires guess of the amplitude)

		// -----
		// Noted that if the peak height is negative then fitting becomes unstable.
		// This is likely when fitting multiple peaks since the initial edge guess 
		// for the background may be wrong.
		// Use the edge value for single peaks but the minimum value in the data if fitting
		// multiple peaks.
		// -----

		float background = 0;

		if (heights.length == 1)
		{
			// Set background using the average value of the edge in the data
			for (int xi = 0; xi < maxx; xi++)
				background += data[xi] + data[maxx * (maxy - 1) + xi];
			for (int yi = 0; yi < maxy; yi++)
				background += data[maxx * yi] + data[maxx * yi + (maxx - 1)];
			background /= 2 * (maxx + maxy);

			// Avoid negative peak height
			if (background > heights[0])
				background = 0;
		}

		if (background == 0)
		{
			// Set background using the minimum value in the data
			background = data[0];
			for (int i = maxx * maxy; --i > 0;)
				background = Math.min(background, data[i]);
		}

		return background;
	}

	/**
	 * Accepts a single array containing 2-dimensional data and a list of the
	 * peaks to fit. Data should be packed in descending dimension order,
	 * e.g. Y,X : Index for [y,z] = MaxX*y + x.
	 * <p>
	 * Performs fitting using the specified method with a Levenberg-Marquardt algorithm.
	 * <p>
	 * Note that if the background parameter is zero it will be assumed that the input amplitude is the total height of
	 * the peak. A new background will be estimated and the heights lowered by this estimated background to create the
	 * amplitude estimate.
	 * <p>
	 * If a peak location is outside the region bounds and has no input width parameters set or from the fit
	 * configuration then fitting will fail (this is because they cannot be estimated).
	 * 
	 * @param data
	 *            The data to fit
	 * @param maxx
	 *            The data size in the x dimension
	 * @param maxy
	 *            The data size in the y dimension
	 * @param npeaks
	 *            The number of peaks
	 * @param params
	 *            The parameters of the peaks. Must have the Amplitude,Xpos,Ypos set. Other
	 *            parameters that are zero will be estimated.
	 * @return The fit result
	 */
	public FitResult fit(final float[] data, final int maxx, final int maxy, final int npeaks, final float[] params)
	{
		FitResult fitResult = null;
		int[] dim = new int[] { maxx, maxy };

		// Working variables
		int[] cumul_region = new int[] { 1, maxx, maxx * maxy };
		int[] position = new int[2];

		// Fitting variables
		float[] y = data; // Value at index
		y_fit = (computeResiduals) ? new float[cumul_region[2]] : null; // Predicted points
		solver = null;
		float[] params_dev = null; // standard deviations for parameters for the fitting function
		double[] error = { 0 }; // The fit Chi-squared value
		int ySize = cumul_region[2];

		final int paramsPerPeak = 6;

		float background = params[0];
		if (background == 0)
		{
			// Extract the heights
			float[] heights = new float[npeaks];
			for (int i = 0, j = 0; i < npeaks; i++, j += paramsPerPeak)
			{
				heights[i] = params[j + Gaussian2DFunction.AMPLITUDE];
			}
			// Get background
			background = getBackground(data, maxx, maxy, heights);
			params[0] = background;
			// Lower heights to get amplitude
			for (int i = 0, j = 0; i < npeaks; i++, j += paramsPerPeak)
			{
				params[j + Gaussian2DFunction.AMPLITUDE] -= background;
			}
		}

		float[] initialParams = Arrays.copyOf(params, params.length);

		int zeroHeight = 0;
		int parameter = 1;

		for (int i = 0, j = 0; i < npeaks; i++, j += paramsPerPeak)
		{
			// Get the parameters
			float height = params[j + Gaussian2DFunction.AMPLITUDE];
			float angle = params[j + Gaussian2DFunction.ANGLE];
			float xpos = params[j + Gaussian2DFunction.X_POSITION];
			float ypos = params[j + Gaussian2DFunction.Y_POSITION];
			float xWidth = params[j + Gaussian2DFunction.X_WIDTH];
			float yWidth = params[j + Gaussian2DFunction.Y_WIDTH];

			// ----
			// Check all input parameters and estimate them if necessary
			// ----

			if (height <= 0)
			{
				// This should not happen if fitting true peaks with an initial height estimate
				zeroHeight++;
				// Set the height to 1 in the event there are multiple peaks to allow fitting
				if (npeaks > 1)
					height = 1;
				//new ij.ImagePlus("zero", new ij.process.FloatProcessor(maxx, maxy, y, null)).show();
			}

			// Set-up for estimating peak width at half maximum 
			position[0] = Math.round(xpos);
			position[1] = Math.round(ypos);
			int index = position[1] * maxx + position[0];

			if (xWidth == 0)
			{
				if (fitConfiguration.getInitialPeakWidth0() > 0)
				{
					xWidth = fitConfiguration.getInitialPeakWidth0();
				}
				else
				{
					// Fail if the width cannot be estimated due to out of bounds
					if (position[0] < 0 || position[0] > maxx || position[1] < 0 || position[1] > maxy)
						return new FitResult(FitStatus.BAD_PARAMETERS, 0, 0, initialParams, null, null, npeaks, 0, null);

					xWidth = half_max_linewidth(y, index, position, dim, 0, cumul_region, background);
				}
			}

			if (yWidth == 0)
			{
				if (fitConfiguration.isWidth1Fitting())
				{
					if (fitConfiguration.getInitialPeakWidth1() > 0)
					{
						yWidth = fitConfiguration.getInitialPeakWidth1();
					}
					else
					{
						// Fail if the width cannot be estimated
						if (position[0] < 0 || position[0] > maxx || position[1] < 0 || position[1] > maxy)
							return new FitResult(FitStatus.BAD_PARAMETERS, 0, 0, initialParams, null, null, npeaks, 0,
									null);

						yWidth = half_max_linewidth(y, index, position, dim, 1, cumul_region, background);
					}
				}
				else
				{
					yWidth = xWidth;
				}
			}

			// Guess the initial angle if input angle is out-of-bounds
			if (angle == 0)
			{
				if (fitConfiguration.isAngleFitting() && fitConfiguration.getInitialAngle() < -Math.PI)
				{
					if (xWidth == yWidth)
						// There is no angle gradient information if the widths are equal. Zero and it will be ignored
						angle = 0;
					else
						// The fit output angle is relative to the major axis.
						// Initialise using the largest width as the major axis (x coordinate):
						angle = (float) ((xWidth > yWidth) ? Math.atan2(yWidth, xWidth) : Math.atan2(xWidth, yWidth));
				}
				else
				{
					angle = fitConfiguration.getInitialAngle();
				}
			}

			// If the position is on the integer grid then use a centre-of-mass approximation
			if (xpos == position[0] && ypos == position[1])
			{
				// Estimate using centre of mass around peak index 
				// Use 0.5 of the width estimate to calculate the range around the index that should be considered.
				// Note: round(x) = (int)(x+0.5)
				int range = (int) ((xWidth + yWidth + 2) / 4);
				double[] com = findCentreOfMass(y, dim, range, position);
				xpos = (float) com[0];
				ypos = (float) com[1];
			}

			// Set all the parameters
			params[parameter++] = height;
			params[parameter++] = angle;
			params[parameter++] = xpos;
			params[parameter++] = ypos;
			params[parameter++] = xWidth;
			params[parameter++] = yWidth;
		}

		// Re-copy the parameters now they have all been set
		initialParams = Arrays.copyOf(params, params.length);

		// Check there are peaks to fit
		if (zeroHeight == npeaks)
			return new FitResult(FitStatus.BAD_PARAMETERS, 0, 0, initialParams, null, null, npeaks, 0, null);

		// -----------------------
		// Use alternative fitters
		// -----------------------

		fitConfiguration.initialise(npeaks, maxx, initialParams);
		solver = fitConfiguration.getFunctionSolver();
		if (fitConfiguration.isComputeDeviations())
			params_dev = new float[params.length];
		final double noise = 0; // //fitConfiguration.getNoise()
		FitStatus result = solver.fit(ySize, y, y_fit, params, params_dev, error, noise);

		// -----------------------

		if (result == FitStatus.OK)
		{
			// Re-assemble all the parameters
			if (!fitConfiguration.isWidth1Fitting() && fitConfiguration.isWidth0Fitting())
			{
				// Ensure Y width is updated with the fitted X width
				for (int i = Gaussian2DFunction.X_WIDTH; i < params.length; i += 6)
				{
					params[i + 1] = params[i];
					if (params_dev != null)
						params_dev[i + 1] = params_dev[i];
				}
			}
			if (fitConfiguration.isAngleFitting())
			{
				// Ensure the angle is within the correct bounds
				for (int i = Gaussian2DFunction.ANGLE; i < params.length; i += 6)
				{
					correctAngle(i, params, params_dev);
				}
			}
			// Ensure widths are positive
			for (int i = params.length - 1; i > 0; i -= paramsPerPeak)
			{
				params[i] = Math.abs(params[i]);
				params[i - 1] = Math.abs(params[i - 1]);
			}

			// TODO - Compare the Chi-squared value for the curve against the standard deviation of 
			// the data, i.e. how a flat plane would fit the data. If no better then mark as a bad fit.

			// Filter peaks only if single peak fitting
			Object statusData = null;
			if (fitConfiguration.isFitValidation())
			{
				result = fitConfiguration.validateFit(npeaks, initialParams, params);
				statusData = fitConfiguration.getValidationData();
			}

			if (computeResiduals)
			{
				for (int i = 0; i < y_fit.length; i++)
					y_fit[i] = y[i] - y_fit[i];
			}

			fitResult = new FitResult(result, Math.max(ySize - solver.getNumberOfFittedParameters(), 0), error[0],
					initialParams, params, params_dev, npeaks, solver.getNumberOfFittedParameters(), statusData);
		}
		else
		{
			y_fit = null;
			fitResult = new FitResult(result, 0, 0, initialParams, null, null, npeaks,
					solver.getNumberOfFittedParameters(), null);
		}

		return fitResult;
	}

	/**
	 * Finds the centre of the image using the centre of mass within the given range of the specified centre-of-mass.
	 */
	private double[] findCentreOfMass(final float[] subImage, final int[] dimensions, final int range,
			final int[] centre)
	{
		int[] min = new int[2];
		int[] max = new int[2];
		for (int i = 2; i-- > 0;)
		{
			min[i] = centre[i] - range;
			max[i] = centre[i] + range;
			if (min[i] < 0)
				min[i] = 0;
			if (max[i] >= dimensions[i] - 1)
				max[i] = dimensions[i] - 1;
		}

		double[] newCom = new double[2];
		double sum = 0;
		for (int y = min[1]; y <= max[1]; y++)
		{
			int index = dimensions[0] * y + min[0];
			for (int x = min[0]; x <= max[0]; x++, index++)
			{
				float value = subImage[index];
				sum += value;
				newCom[0] += x * value;
				newCom[1] += y * value;
			}
		}

		for (int i = 2; i-- > 0;)
		{
			newCom[i] /= sum;
		}

		return newCom;
	}

	/**
	 * Swap the axes so that the major axis is the X axis.
	 * Correct the fit angle to lie within the -90-90 degree domain from the major-axis
	 * 
	 * @param i
	 *            The angle position within the parameter array
	 * @param params
	 * @param params_dev
	 */
	private void correctAngle(final int i, final float[] params, final float[] params_dev)
	{
		double angle = params[i];

		final double twicePI = 2 * Math.PI;
		double fixed = (angle + Math.PI) % twicePI;
		if (fixed < 0)
		{
			fixed += twicePI;
		}
		angle = fixed - Math.PI;

		//		// Angle should now be in -180 - 180 degrees domain
		//		if (angle < -Math.PI || angle > Math.PI)
		//		{
		//			System.out.printf("angle error %g != %g\n", angle, Math.asin(Math.sin(params[i])));
		//		}

		// Commented out as this interferes with the PSF Estimator
		float xWidth = params[i + 3];
		float yWidth = params[i + 4];
		// The fit will compute the angle from the major axis. 
		// Standardise so it is always from the X-axis
		if (yWidth > xWidth)
		{
			swap(i + 3, params);
			if (params_dev != null)
				swap(i + 3, params_dev);

			// Rotate 90 degrees
			angle += Math.PI / 2.0;
			// Check domain
			if (angle > Math.PI)
			{
				angle -= 2 * Math.PI;
			}
		}

		// Return in 0 - 180 degrees domain since the Gaussian has 2-fold symmetry,
		// i.e. angle -10 == 170
		params[i] = (float) ((angle < 0) ? angle + Math.PI : angle);

		// Return in -90 - 90 degrees domain since 0 should be no angle
		params[i] -= Math.PI / 2;
	}

	private void swap(final int i, final float[] params)
	{
		float tmp = params[i];
		params[i] = params[i + 1];
		params[i + 1] = tmp;
	}

	/**
	 * Convert the Full-Width at Half-Maximum to the Standard Deviation
	 * 
	 * @param fwhm
	 * @return sd
	 */
	public static float fwhm2sd(float fwhm)
	{
		return (float) (fwhm / (2 * Math.sqrt(2 * Math.log(2))));
	}

	/**
	 * Convert the Standard Deviation to the Full-Width at Half-Maximum
	 * 
	 * @param sd
	 * @return fwhm
	 */
	public static float sd2fwhm(final float sd)
	{
		return (float) (sd * 2 * Math.sqrt(2 * Math.log(2)));
	}

	/**
	 * @return the residuals from the last successful fit. If fitting failed then this is null.
	 */
	public float[] getResiduals()
	{
		return y_fit;
	}

	/**
	 * @return the computeResiduals
	 */
	public boolean isComputeResiduals()
	{
		return computeResiduals;
	}

	/**
	 * @param computeResiduals
	 *            Set to true to compute the residuals
	 */
	public void setComputeResiduals(final boolean computeResiduals)
	{
		this.computeResiduals = computeResiduals;
	}

	/**
	 * Return true if the last call to a fit(...) method created a function solver. This allows the properties to be
	 * accessed for the last fit. Otherwise the properties will return zero.
	 * 
	 * @return True if the last call to a fit(...) method created a function solver
	 */
	public boolean solvedLastFit()
	{
		return (solver != null);
	}

	/**
	 * @return The number of iterations used in the last fit
	 */
	public int getIterations()
	{
		return (solver != null) ? solver.getIterations() : 0;
	}

	/**
	 * @return (solver != null) ? the totalSumOfSquares for the last fit
	 */
	public double getTotalSumOfSquares()
	{
		return (solver != null) ? solver.getTotalSumOfSquares() : 0;
	}

	/**
	 * @return (solver != null) ? the finalResidualSumOfSquares for the last fit
	 */
	public double getFinalResidualSumOfSquares()
	{
		return (solver != null) ? solver.getResidualSumOfSquares() : 0;
	}

	/**
	 * @return (solver != null) ? the numberOfFittedParameters for the last fit
	 */
	public int getNumberOfFittedParameters()
	{
		return (solver != null) ? solver.getNumberOfFittedParameters() : 0;
	}

	/**
	 * @return (solver != null) ? the numberOfFittedPoints for the last fit
	 */
	public int getNumberOfFittedPoints()
	{
		return (solver != null) ? solver.getNumberOfFittedPoints() : 0;
	}
}