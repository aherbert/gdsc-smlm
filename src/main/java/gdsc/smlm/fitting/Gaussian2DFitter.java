package gdsc.smlm.fitting;

import java.util.Arrays;

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
 * 
 * This is an adaption of the C-code contained in the CcpNmr Analysis Program:
 * CCPN website (http://www.ccpn.ac.uk/)
 *---------------------------------------------------------------------------*/

import gdsc.smlm.function.gaussian.Gaussian2DFunction;

/**
 * Fits a 2-dimensional Gaussian function for the specified peak. Can optionally fit an elliptical Gaussian function.
 * <p>
 * Performs fitting using the configured algorithm.
 */
public class Gaussian2DFitter
{
	private FitConfiguration fitConfiguration;
	private FunctionSolver solver;

	// The last successful fit. Used to compute the residuals.
	private double[] y_fit = null;
	// Allow calculation of residuals to be turned off (overwrite constructor fit configuration)
	private boolean computeResiduals = true;

	private int maximumWidthFactor = 2;

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

	static double half_max_position(double[] data, int index, int[] point, int[] dim, int dimension, int[] cumul_region,
			int dirn, double background)
	{
		int i, i_start, i_end, i_step;
		double v = data[index];
		double v_half = 0.5f * (v + background);
		double v_prev = v, v_this;
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

	public static double half_max_linewidth(double[] data, int index, int[] point, int[] dim, int dimension,
			int[] cumul_region, double background)
	{
		double linewidth, a, b;

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
	public FitResult fit(final double[] data, final int maxx, final int maxy, final int[] peaks)
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
	public FitResult fit(final double[] data, final int maxx, final int maxy, final int[] peaks, double[] heights)
	{
		int npeaks = peaks.length;

		final int paramsPerPeak = 6;

		double[] params = new double[1 + paramsPerPeak * npeaks];

		// Get peak heights (if multiple peaks)
		boolean amplitudeEstimate = true;
		if (npeaks > 1 && (heights == null || heights.length != peaks.length))
		{
			heights = new double[peaks.length];
			for (int i = 0; i < peaks.length; i++)
				heights[i] = data[peaks[i]];
		}

		final double background = getBackground(data, maxx, maxy, npeaks);

		// Set the initial parameters
		params[0] = background;

		if (npeaks == 1)
		{
			double sum = 0;
			final int size = maxx * maxy;
			for (int i = size; i-- > 0;)
				sum += data[i];
			params[Gaussian2DFunction.SIGNAL] = sum - background * size;
			params[Gaussian2DFunction.X_POSITION] = peaks[0] % maxx;
			params[Gaussian2DFunction.Y_POSITION] = peaks[0] / maxx;
			amplitudeEstimate = false;
		}
		else
		{
			for (int i = 0, j = 0; i < peaks.length; i++, j += paramsPerPeak)
			{
				int index = peaks[i];
				params[j + Gaussian2DFunction.SIGNAL] = heights[i] - background;
				params[j + Gaussian2DFunction.X_POSITION] = index % maxx;
				params[j + Gaussian2DFunction.Y_POSITION] = index / maxx;
			}
		}

		// We have estimated the background already
		return fit(data, maxx, maxy, npeaks, params, amplitudeEstimate, true);
	}

	/**
	 * Guess the background from the data given the estimated peak heights.
	 * <p>
	 * For a single peak the method assumes the peak is in the centre. In this case the average edge value of the data
	 * is used.
	 * <p>
	 * If multiple peaks heights are provided then always use the minimum value in the data since it cannot be assumed
	 * that all peaks are away from the edge of the data.
	 * 
	 * @param data
	 * @param maxx
	 * @param maxy
	 * @param npeaks
	 * @return The background estimate
	 */
	public static double getBackground(final double[] data, final int maxx, final int maxy, final int npeaks)
	{
		// TODO - What is the best method for setting the background?
		// 1. Min in data
		// 2. Average of border?
		// 3. Evaluate total volume under the initial Gaussian params and subtract that from the sum of the image?
		//    (note initial gaussian requires guess of the amplitude which needs background (or bias))

		// -----
		// Noted that if the peak height is negative then fitting becomes unstable.
		// This is likely when fitting multiple peaks since the initial edge guess 
		// for the background may be wrong.
		// Use the edge value for single peaks but the minimum value in the data if fitting
		// multiple peaks.
		// -----

		double background = 0;

		//Utils.display("Spot", data, maxx, maxy);

		if (npeaks == 1)
		{
			// Set background using the average value of the edge in the data
			final int s2 = getIndex(0, maxy - 1, maxx);
			for (int xi = 0, xi2 = s2; xi < maxx; xi++, xi2++)
				background += data[xi] + data[xi2];
			for (int yi = maxx, yi2 = getIndex(maxx - 1, 1, maxx); yi < s2; yi += maxx, yi2 += maxx)
				background += data[yi] + data[yi2];
			background /= 2 * (maxx + maxy - 2);
		}
		else
		{
			// Set background using the minimum value in the data
			background = data[0];
			for (int i = maxx * maxy; --i > 0;)
				if (background > data[i])
					background = data[i];
		}

		return background;
	}

	private static int getIndex(int x, int y, int maxx)
	{
		return y * maxx + x;
	}

	/**
	 * Accepts a single array containing 2-dimensional data and a list of the
	 * peaks to fit. Data should be packed in descending dimension order,
	 * e.g. Y,X : Index for [y,z] = MaxX*y + x.
	 * <p>
	 * Performs fitting using the specified method with a Levenberg-Marquardt algorithm.
	 * <p>
	 * The input parameter can estimate the signal (the total volume of the Gaussian) or the amplitude (the height of
	 * the Gaussian). The signal = amplitude * 2 * pi * sd0 * sd1. The amplitude is the recommended method to estimate
	 * parameters for multiple peaks. The signal can be estimated for a single peak by summing all the pixels (minus the
	 * background).
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
	 *            The parameters of the peaks. Must have the Signal/Amplitude,Xpos,Ypos set. Other
	 *            parameters that are zero will be estimated.
	 * @param amplitudeEstimate
	 *            Set to true if the parameters have amplitude estimated in the
	 *            {@link gdsc.smlm.function.Gaussian2DFunction.SIGNAL} field. The
	 *            default is signal.
	 * @return The fit result
	 */
	public FitResult fit(final double[] data, final int maxx, final int maxy, final int npeaks, final double[] params,
			final boolean amplitudeEstimate)
	{
		return fit(data, maxx, maxy, npeaks, params, amplitudeEstimate, false);
	}

	/**
	 * Accepts a single array containing 2-dimensional data and a list of the
	 * peaks to fit. Data should be packed in descending dimension order,
	 * e.g. Y,X : Index for [y,z] = MaxX*y + x.
	 * <p>
	 * Performs fitting using the specified method with a Levenberg-Marquardt algorithm.
	 * <p>
	 * The input parameter can estimate the signal (the total volume of the Gaussian) or the amplitude (the height of
	 * the Gaussian). The signal = amplitude * 2 * pi * sd0 * sd1. The amplitude is the recommended method to estimate
	 * parameters for multiple peaks. The signal can be estimated for a single peak by summing all the pixels (minus the
	 * background).
	 * <p>
	 * Note that if the background parameter is zero it will be assumed that the input signal/amplitude is the total
	 * volume/height of the peak. A new background will be estimated and the volume/heights lowered by this estimated
	 * background to create the new parameters.
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
	 *            The parameters of the peaks. Must have the Signal/Amplitude,Xpos,Ypos set. Other
	 *            parameters that are zero will be estimated. This can optionally be ignored for the background
	 *            parameter which is valid if zero.
	 * @param zeroBackground
	 *            Set to true if a zero value for the background parameter is the estimate
	 * @param amplitudeEstimate
	 *            Set to true if the parameters have amplitude estimated in the
	 *            {@link gdsc.smlm.function.Gaussian2DFunction.SIGNAL} field. The
	 *            default is signal.
	 * @return The fit result
	 */
	public FitResult fit(final double[] data, final int maxx, final int maxy, final int npeaks, final double[] params,
			final boolean amplitudeEstimate, final boolean zeroBackground)
	{
		FitResult fitResult = null;
		final int[] dim = new int[] { maxx, maxy };

		// Working variables
		final int[] cumul_region = new int[] { 1, maxx, maxx * maxy };
		final int[] position = new int[2];

		// Fitting variables
		final double[] y = data; // Value at index
		y_fit = (computeResiduals) ? new double[cumul_region[2]] : null; // Predicted points
		solver = null;
		double[] params_dev = null; // standard deviations for parameters for the fitting function
		double[] error = { 0 }; // The fit Chi-squared value
		final int ySize = cumul_region[2];

		final int paramsPerPeak = 6;

		double background = params[0];
		if (background == 0 && !zeroBackground)
		{
			// Get background
			background = getBackground(data, maxx, maxy, npeaks);
			params[0] = background;
			if (amplitudeEstimate)
			{
				// For a single peak, check the height is above background
				if (npeaks == 1 && params[Gaussian2DFunction.SIGNAL] < background)
				{
					// Set the background to the min value in the data using the multiple peak option
					background = getBackground(data, maxx, maxy, 2);

					// Check if still below background
					if (params[Gaussian2DFunction.SIGNAL] < background)
					{
						// Set the height to the max value in the data
						double yMax = y[0];
						for (int i = 1; i < ySize; i++)
							if (yMax < y[i])
								yMax = y[i];
						params[Gaussian2DFunction.SIGNAL] = (float) yMax;
					}
				}

				// Lower heights to get amplitude
				for (int j = Gaussian2DFunction.SIGNAL; j < params.length; j += paramsPerPeak)
				{
					params[j] -= background;
				}
			}
		}

		double[] initialParams = Arrays.copyOf(params, params.length);

		// Check all the heights are valid first
		int zeroHeight = 0;
		for (int j = Gaussian2DFunction.SIGNAL; j < params.length; j += paramsPerPeak)
		{
			if (params[j] <= 0)
				zeroHeight++;
		}

		// Check there are peaks to fit
		if (zeroHeight == npeaks)
			return new FitResult(FitStatus.BAD_PARAMETERS, 0, 0, initialParams, null, null, npeaks, 0, null);

		// Set all zero height peaks to a fraction of the maximum to allow fitting
		if (zeroHeight > 0)
		{
			double max = 0;
			for (int j = Gaussian2DFunction.SIGNAL; j < params.length; j += paramsPerPeak)
			{
				if (max < params[j])
					max = params[j];
			}
			max *= 0.1; // Use fraction of the max peak
			for (int j = Gaussian2DFunction.SIGNAL; j < params.length; j += paramsPerPeak)
			{
				if (params[j] <= 0)
					params[j] = max;
			}
		}

		int parameter = 1;
		for (int i = 0, j = 0; i < npeaks; i++, j += paramsPerPeak)
		{
			// Get the parameters
			double signal = params[j + Gaussian2DFunction.SIGNAL];
			double angle = params[j + Gaussian2DFunction.ANGLE];
			double xpos = params[j + Gaussian2DFunction.X_POSITION];
			double ypos = params[j + Gaussian2DFunction.Y_POSITION];
			double sx = params[j + Gaussian2DFunction.X_SD];
			double sy = params[j + Gaussian2DFunction.Y_SD];

			// ----
			// Check all input parameters and estimate them if necessary
			// ----

			// Set-up for estimating peak width at half maximum 
			position[0] = (int) Math.round(xpos);
			position[1] = (int) Math.round(ypos);
			int index = position[1] * maxx + position[0];

			if (sx == 0)
			{
				if (fitConfiguration.getInitialPeakStdDev0() > 0)
				{
					sx = fitConfiguration.getInitialPeakStdDev0();
				}
				else
				{
					// Fail if the width cannot be estimated due to out of bounds
					if (position[0] < 0 || position[0] > maxx || position[1] < 0 || position[1] > maxy)
						return new FitResult(FitStatus.BAD_PARAMETERS, 0, 0, initialParams, null, null, npeaks, 0,
								null);

					sx = fwhm2sd(half_max_linewidth(y, index, position, dim, 0, cumul_region, background));
				}
			}

			if (sy == 0)
			{
				if (fitConfiguration.isWidth1Fitting())
				{
					if (fitConfiguration.getInitialPeakStdDev1() > 0)
					{
						sy = fitConfiguration.getInitialPeakStdDev1();
					}
					else
					{
						// Fail if the width cannot be estimated
						if (position[0] < 0 || position[0] > maxx || position[1] < 0 || position[1] > maxy)
							return new FitResult(FitStatus.BAD_PARAMETERS, 0, 0, initialParams, null, null, npeaks, 0,
									null);

						sy = fwhm2sd(half_max_linewidth(y, index, position, dim, 1, cumul_region, background));
					}
				}
				else
				{
					sy = sx;
				}
			}

			// Guess the initial angle if input angle is out-of-bounds
			if (angle == 0)
			{
				if (fitConfiguration.isAngleFitting() && fitConfiguration.getInitialAngle() >= -Math.PI &&
						fitConfiguration.getInitialAngle() <= -Math.PI)
				{
					if (sx != sy)
					{
						// There is no angle gradient information if the widths are equal. Zero and it will be ignored
						angle = fitConfiguration.getInitialAngle();
					}
				}
			}

			// If the position is on the integer grid then use a centre-of-mass approximation
			if (xpos == position[0] && ypos == position[1])
			{
				// Estimate using centre of mass around peak index 
				// Use 2 * SD estimate to calculate the range around the index that should be considered.
				// SD = (sx+sy)/2 => Range = sx+sy
				final int range = (int) Math.ceil(sx + sy + 0.5);
				final double[] com = findCentreOfMass(y, dim, range, position);
				xpos = (double) com[0];
				ypos = (double) com[1];
			}

			// Convert amplitudes to signal
			if (amplitudeEstimate)
				signal *= 2 * Math.PI * sx * sy;

			// Set all the parameters
			params[parameter++] = signal;
			params[parameter++] = angle;
			params[parameter++] = xpos;
			params[parameter++] = ypos;
			params[parameter++] = sx;
			params[parameter++] = sy;
		}

		// Re-copy the parameters now they have all been set
		initialParams = Arrays.copyOf(params, params.length);

		// -----------------------
		// Use alternative fitters
		// -----------------------

		fitConfiguration.initialise(npeaks, maxx, initialParams);
		solver = fitConfiguration.getFunctionSolver();
		if (fitConfiguration.isComputeDeviations())
			params_dev = new double[params.length];

		// Subtract the bias
		double bias = 0;
		if (fitConfiguration.isRemoveBiasBeforeFitting())
		{
			// Some methods can fit negative data, e.g. PoissonGaussian or PoissonGammaGaussian.
			if (background < fitConfiguration.getBias())
			{
				// Debugging: remove this
				//System.out.printf("Background %f < Bias %f\n", background, fitConfiguration.getBias());
			}

			// No negative data
			//bias = FastMath.min(background, fitConfiguration.getBias());

			// Subtract the full bias. Leave it to the solver to handle negative data.
			bias = fitConfiguration.getBias();

			params[0] -= bias;
			for (int i = 0; i < ySize; i++)
				y[i] -= bias;
		}

		// Bounds are more restrictive than constraints
		if (solver.isBounded())
		{
			setBounds(maxx, maxy, npeaks, params, y, ySize, paramsPerPeak);
		}
		if (solver.isConstrained())
		{
			setConstraints(maxx, maxy, npeaks, params, y, ySize, paramsPerPeak);
		}

		final double noise = 0; // //fitConfiguration.getNoise()
		FitStatus result = solver.fit(ySize, y, y_fit, params, params_dev, error, noise);

		// -----------------------

		if (result == FitStatus.OK)
		{
			// Add the bias back to the background
			params[0] += bias;

			// Re-assemble all the parameters
			if (!fitConfiguration.isWidth1Fitting() && fitConfiguration.isWidth0Fitting())
			{
				// Ensure Y width is updated with the fitted X width
				for (int i = Gaussian2DFunction.X_SD; i < params.length; i += 6)
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

			fitResult = new FitResult(result, FastMath.max(ySize - solver.getNumberOfFittedParameters(), 0), error[0],
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
	private double[] findCentreOfMass(final double[] subImage, final int[] dimensions, final int range,
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
				double value = subImage[index];
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
	 * Sets the bounds for the fitted parameters
	 * 
	 * @param maxx
	 *            The x range of the data
	 * @param maxy
	 *            The y range of the data
	 * @param npeaks
	 *            The number of peaks
	 * @param params
	 *            The estimated parameters
	 * @param y
	 *            The data
	 * @param ySize
	 *            The size of the data
	 * @param paramsPerPeak
	 *            The number of parameters per peak
	 */
	private void setBounds(final int maxx, final int maxy, final int npeaks, final double[] params, final double[] y,
			final int ySize, final int paramsPerPeak)
	{
		// Create appropriate bounds for the parameters
		double[] lower = new double[params.length];
		double[] upper = new double[lower.length];
		double yMax = y[0];
		double yMin = y[0];
		for (int i = 1; i < ySize; i++)
		{
			if (yMax < y[i])
				yMax = y[i];
			else if (yMin > y[i])
				yMin = y[i];
		}
		if (fitConfiguration.isBackgroundFitting())
		{
			if (yMax > params[0])
				upper[0] = yMax;
			else
				upper[0] = params[0] + (params[0] - yMax);

			if (yMin < params[0])
				lower[0] = yMin;
			else
				lower[0] = params[0] - (yMin - params[0]);
		}

		final double wf = (fitConfiguration.getWidthFactor() > 1 &&
				fitConfiguration.getWidthFactor() < maximumWidthFactor) ? fitConfiguration.getWidthFactor()
						: maximumWidthFactor;

		// TODO - Check if the signal bounds are appropriate
		if (npeaks == 1)
		{
			// Allow the signal to explain all the data. This assumes the data window entirely covers the spot. 
			double sum = 0;
			for (int i = 1; i < ySize; i++)
				sum += y[i];
			// Increase sum by 2 to allow for error
			upper[Gaussian2DFunction.SIGNAL] = 2 * sum - yMin * ySize;
		}
		else
		{
			final double height = yMax - yMin;
			// Signal = height * 2 * pi * sd0 * sd1
			// Allow a maximum using the width factor that defines the bounds on the width.
			// Increase the height by 2 to allow for error.
			final double factor = 2 * height * 2 * Math.PI * wf * wf;
			//System.out.printf("%f or %f\n", upper[Gaussian2DFunction.SIGNAL], factor * params[Gaussian2DFunction.X_SD] *
			//		params[Gaussian2DFunction.Y_SD]);
			for (int i = 0, j = 0; i < npeaks; i++, j += paramsPerPeak)
			{
				upper[j + Gaussian2DFunction.SIGNAL] = factor * params[j + Gaussian2DFunction.X_SD] *
						params[j + Gaussian2DFunction.Y_SD];
			}
		}

		for (int i = 0, j = 0; i < npeaks; i++, j += paramsPerPeak)
		{
			if (params[j + Gaussian2DFunction.SIGNAL] < lower[j + Gaussian2DFunction.SIGNAL])
				lower[j + Gaussian2DFunction.SIGNAL] = params[j + Gaussian2DFunction.SIGNAL] -
						(lower[j + Gaussian2DFunction.SIGNAL] - params[j + Gaussian2DFunction.SIGNAL]);
			if (params[j + Gaussian2DFunction.SIGNAL] > upper[j + Gaussian2DFunction.SIGNAL])
				upper[j + Gaussian2DFunction.SIGNAL] = params[j + Gaussian2DFunction.SIGNAL] +
						(params[j + Gaussian2DFunction.SIGNAL] - upper[j + Gaussian2DFunction.SIGNAL]);

			// All functions evaluate the x and y position.
			// Lower bounds on these will be zero when the array is initialised.
			// We may have an estimate outside the bounds (if including neighbours).
			upper[j + Gaussian2DFunction.X_POSITION] = Math.max(maxx,
					params[j + Gaussian2DFunction.X_POSITION] + params[j + Gaussian2DFunction.X_SD]);
			upper[j + Gaussian2DFunction.Y_POSITION] = Math.max(maxy,
					params[j + Gaussian2DFunction.Y_POSITION] + params[j + Gaussian2DFunction.Y_SD]);

			lower[j + Gaussian2DFunction.X_POSITION] = Math.min(0,
					params[j + Gaussian2DFunction.X_POSITION] - params[j + Gaussian2DFunction.X_SD]);
			lower[j + Gaussian2DFunction.Y_POSITION] = Math.min(0,
					params[j + Gaussian2DFunction.Y_POSITION] - params[j + Gaussian2DFunction.Y_SD]);

			if (fitConfiguration.isAngleFitting())
			{
				lower[j + Gaussian2DFunction.ANGLE] = -Math.PI;
				upper[j + Gaussian2DFunction.ANGLE] = Math.PI;
			}
			if (fitConfiguration.isWidth0Fitting())
			{
				lower[j + Gaussian2DFunction.X_SD] = params[j + Gaussian2DFunction.X_SD] / wf;
				upper[j + Gaussian2DFunction.X_SD] = params[j + Gaussian2DFunction.X_SD] * wf;
			}
			if (fitConfiguration.isWidth1Fitting())
			{
				lower[j + Gaussian2DFunction.Y_SD] = params[j + Gaussian2DFunction.Y_SD] / wf;
				upper[j + Gaussian2DFunction.Y_SD] = params[j + Gaussian2DFunction.Y_SD] * wf;
			}
		}

		// TODO - Comment this out
		// Debug check
		for (int i = 0; i < params.length; i++)
		{
			if (params[i] < lower[i])
			{
				System.out.printf("Param %d too low %f < %f\n", i, params[i], lower[i]);
				lower[i] = params[i] - (lower[i] - params[i]);
			}
			if (params[i] > upper[i])
			{
				System.out.printf("Param %d too high %f > %f\n", i, params[i], upper[i]);
				upper[i] = params[i] + (params[i] - upper[i]);
			}
		}

		solver.setBounds(lower, upper);
	}

	/**
	 * Sets the constraints for the fitted parameters. This functions set the lower bounds of the signal to zero and
	 * background to zero (or negative if the background estimate is < 0).
	 * 
	 * @param maxx
	 *            The x range of the data
	 * @param maxy
	 *            The y range of the data
	 * @param npeaks
	 *            The number of peaks
	 * @param params
	 *            The estimated parameters
	 * @param y
	 *            The data
	 * @param ySize
	 *            The size of the data
	 * @param paramsPerPeak
	 *            The number of parameters per peak
	 */
	private void setConstraints(final int maxx, final int maxy, final int npeaks, final double[] params,
			final double[] y, final int ySize, final int paramsPerPeak)
	{
		// Create appropriate bounds for the parameters
		double[] lower = new double[params.length];
		double[] upper = new double[lower.length];
		Arrays.fill(lower, Float.NEGATIVE_INFINITY);
		Arrays.fill(upper, Float.POSITIVE_INFINITY);

		// If the bias is subtracted then we may have negative data and a background estimate that is negative
		if (params[0] < 0)
		{
			double yMin = 0;
			for (int i = 0; i < ySize; i++)
			{
				if (yMin > y[i])
					yMin = y[i];
			}
			if (yMin < params[0])
				lower[0] = yMin;
			else
				lower[0] = params[0] - (yMin - params[0]);
		}

		for (int i = 0, j = 0; i < npeaks; i++, j += paramsPerPeak)
		{
			lower[j + Gaussian2DFunction.SIGNAL] = 0;
		}
		solver.setConstraints(lower, upper);
	}

	/**
	 * Swap the axes so that the major axis is the X axis.
	 * Correct the fit angle to lie within the 0-180 degree domain from the major-axis.
	 * 
	 * @param i
	 *            The angle position within the parameter array
	 * @param params
	 * @param params_dev
	 */
	private void correctAngle(final int i, final double[] params, final double[] params_dev)
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
		double xWidth = params[i + Gaussian2DFunction.X_SD - Gaussian2DFunction.ANGLE];
		double yWidth = params[i + Gaussian2DFunction.Y_SD - Gaussian2DFunction.ANGLE];
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
		params[i] = (angle < 0) ? angle + Math.PI : angle;
	}

	private void swap(final int i, final double[] params)
	{
		double tmp = params[i];
		params[i] = params[i + 1];
		params[i + 1] = tmp;
	}

	/**
	 * Convert the Full-Width at Half-Maximum to the Standard Deviation
	 * 
	 * @param fwhm
	 * @return sd
	 */
	public static double fwhm2sd(double fwhm)
	{
		return (double) (fwhm / (2 * Math.sqrt(2 * Math.log(2))));
	}

	/**
	 * Convert the Standard Deviation to the Full-Width at Half-Maximum
	 * 
	 * @param sd
	 * @return fwhm
	 */
	public static double sd2fwhm(final double sd)
	{
		return (double) (sd * 2 * Math.sqrt(2 * Math.log(2)));
	}

	/**
	 * @return the residuals from the last successful fit. If fitting failed then this is null.
	 */
	public double[] getResiduals()
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
	 * @return the totalSumOfSquares for the last fit
	 */
	public double getTotalSumOfSquares()
	{
		return (solver != null) ? solver.getTotalSumOfSquares() : 0;
	}

	/**
	 * @return the finalResidualSumOfSquares for the last fit
	 */
	public double getFinalResidualSumOfSquares()
	{
		return (solver != null) ? solver.getResidualSumOfSquares() : 0;
	}

	/**
	 * @return the numberOfFittedParameters for the last fit
	 */
	public int getNumberOfFittedParameters()
	{
		return (solver != null) ? solver.getNumberOfFittedParameters() : 0;
	}

	/**
	 * @return the numberOfFittedPoints for the last fit
	 */
	public int getNumberOfFittedPoints()
	{
		return (solver != null) ? solver.getNumberOfFittedPoints() : 0;
	}

	/**
	 * @return the maximum width factor to use to set the bounds for width parameter fitting
	 */
	public int getMaximumWidthFactor()
	{
		return maximumWidthFactor;
	}

	/**
	 * Set the maximum width factor to use to set the bounds for width parameter fitting. If the fit configuration has a
	 * smaller width factor then that will be used instead.
	 * <p>
	 * The bounds are set using the initial estimate w in the range w/f to w*f.
	 * 
	 * @param maximumWidthFactor
	 *            the maximum width factor to use to set the bounds for width parameter fitting (must be above 1)
	 */
	public void setMaximumWidthFactor(int maximumWidthFactor)
	{
		if (maximumWidthFactor > 1)
			this.maximumWidthFactor = maximumWidthFactor;
	}

	/**
	 * @return the optimised function value for the last fit
	 */
	public double getValue()
	{
		return (solver != null) ? solver.getValue() : 0;
	}
}