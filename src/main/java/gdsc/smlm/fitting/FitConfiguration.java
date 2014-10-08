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
import gdsc.smlm.fitting.function.GaussianFunction;
import gdsc.smlm.fitting.function.GaussianFunctionFactory;
import gdsc.smlm.fitting.function.NoiseModel;
import gdsc.smlm.fitting.logging.Logger;
import gdsc.smlm.fitting.nonlinear.ApacheLVMFitter;
import gdsc.smlm.fitting.nonlinear.MaximumLikelihoodFitter;
import gdsc.smlm.fitting.nonlinear.NonLinearFit;
import gdsc.smlm.fitting.nonlinear.StoppingCriteria;
import gdsc.smlm.fitting.nonlinear.stop.ErrorStoppingCriteria;
import gdsc.smlm.fitting.nonlinear.stop.GaussianStoppingCriteria;
import gdsc.smlm.fitting.nonlinear.stop.ParameterStoppingCriteria;
import gdsc.smlm.results.PeakResult;

/**
 * Specifies the fitting configuration for Gaussian fitting
 */
public class FitConfiguration implements Cloneable
{
	private FitCriteria fitCriteria = FitCriteria.LEAST_SQUARED_ERROR;
	private Logger log = null;
	private double delta = 0.0001;
	private float initialAngle = 0; // Radians
	private float initialSD0 = 1;
	private float initialSD1 = 1;
	private boolean computeDeviations = false;
	private FitSolver fitSolver = FitSolver.LVM;
	private int minIterations = 0;
	private int maxIterations = 20;
	private int significantDigits = 5;
	private FitFunction fitFunction;
	private int flags;
	private boolean backgroundFitting = true;
	private float coordinateShift = 1;
	private float signalThreshold = 0;
	private float signalStrength = 30;
	private float precisionThreshold = 0;
	private double nmPerPixel = 0;
	private float gain = 0;
	private float noise = 0;
	private float widthFactor = 2;
	private boolean fitValidation = true;
	private double lambda = 10;
	private boolean computeResiduals = true;
	private float duplicateDistance = 0.5f;

	private StoppingCriteria stoppingCriteria = null;
	private GaussianFunction gaussianFunction = null;
	private NoiseModel noiseModel = null;

	/**
	 * Default constructor
	 */
	public FitConfiguration()
	{
		setFitFunction(FitFunction.CIRCULAR);
	}

	/**
	 * Creates the appropriate stopping criteria and Gaussian function for the configuration
	 * 
	 * @param npeaks
	 *            The number of peaks to fit
	 * @param maxx
	 *            The width of the XY data
	 * @param params
	 *            The Gaussian parameters
	 */
	public void initialise(int npeaks, int maxx, float[] params)
	{
		gaussianFunction = createGaussianFunction(npeaks, maxx, params);
		stoppingCriteria = createStoppingCriteria(gaussianFunction, params);
	}

	/**
	 * Creates the appropriate stopping criteria for the configuration
	 * 
	 * @param func
	 *            The Gaussian function
	 * @param params
	 *            The Gaussian parameters
	 * @return The stopping criteria
	 */
	public StoppingCriteria createStoppingCriteria(GaussianFunction func, float[] params)
	{
		StoppingCriteria stoppingCriteria;
		if (fitCriteria == FitCriteria.PARAMETERS)
		{
			ParameterStoppingCriteria sc = new ParameterStoppingCriteria(func);
			sc.setSignificantDigits(significantDigits);
			addPeakRestrictions(func, sc, params);
			stoppingCriteria = sc;
		}
		else if (fitCriteria == FitCriteria.COORDINATES)
		{
			GaussianStoppingCriteria sc = new GaussianStoppingCriteria(func);
			sc.setDelta(delta);
			addPeakRestrictions(func, sc, params);
			stoppingCriteria = sc;
		}
		else
		{
			ErrorStoppingCriteria sc = new ErrorStoppingCriteria();
			sc.setSignificantDigits(significantDigits);
			sc.setAvoidPlateau(fitCriteria == FitCriteria.LEAST_SQUARED_PLUS);
			stoppingCriteria = sc;
		}
		// Removed to reduce verbosity of logging output
		//stoppingCriteria.setLog(log);
		stoppingCriteria.setMinimumIterations(minIterations);
		stoppingCriteria.setMaximumIterations(maxIterations);
		return stoppingCriteria;
	}

	/**
	 * Add restrictions on peak during fitting
	 * 
	 * @param func
	 * @param sc
	 */
	protected void addPeakRestrictions(GaussianFunction func, GaussianStoppingCriteria sc, float[] params)
	{
		// TODO - Check if it is worth using this to stop fitting early or simply do it at the end.
		sc.setMinimumAmplitude(0);

		//		sc.setMinimumPosition(new float[] { 2, 2 });
		//		sc.setMaximumPosition(new float[] { func.getDimensions()[0], func.getDimensions()[0] });
		//
		//		sc.setMinimumWidth(new float[] { 2, 2 });
		//		sc.setMaximumWidth(new float[] { 12, 12 });
	}

	/**
	 * Creates the appropriate 2D Gaussian function for the configuration
	 * 
	 * @param npeaks
	 *            The number of peaks to fit
	 * @param maxx
	 *            The width of the XY data
	 * @param params
	 *            The Gaussian parameters
	 * @return The function
	 */
	public GaussianFunction createGaussianFunction(int npeaks, int maxx, float[] params)
	{
		int flags = this.flags;

		if (isBackgroundFitting())
			flags |= GaussianFunctionFactory.FIT_BACKGROUND;

		GaussianFunction f = GaussianFunctionFactory.create2D(npeaks, maxx, flags);
		//f.initialise(params);
		return f;
	}

	/**
	 * @param fitCriteria
	 *            the fit criteria to set
	 */
	public void setFitCriteria(int fitCriteria)
	{
		if (fitCriteria >= 0 && fitCriteria < FitCriteria.values().length)
		{
			setFitCriteria(FitCriteria.values()[fitCriteria]);
		}
	}

	/**
	 * @param fitCriteria
	 *            the fit criteria to set
	 */
	public void setFitCriteria(FitCriteria fitCriteria)
	{
		this.fitCriteria = fitCriteria;
	}

	/**
	 * @return the fit criteria
	 */
	public FitCriteria getFitCriteria()
	{
		return fitCriteria;
	}

	/**
	 * @param log
	 *            the log to set. Used to output fit evaluations for each iteration
	 */
	public void setLog(Logger log)
	{
		this.log = log;
	}

	/**
	 * @return the log
	 */
	public Logger getLog()
	{
		return log;
	}

	/**
	 * @param delta
	 *            the delta to set for the fitting criteria
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
	 * @param initialAngle
	 *            the initialAngle to set (in radians)
	 */
	public void setInitialAngle(float initialAngle)
	{
		this.initialAngle = initialAngle;
	}

	/**
	 * @param initialAngle
	 *            the initialAngle to set (in degrees)
	 */
	public void setInitialAngleD(float initialAngle)
	{
		this.initialAngle = (float) (initialAngle * Math.PI / 180.0);
	}

	/**
	 * @return the initialAngle (in radians)
	 */
	public float getInitialAngle()
	{
		return initialAngle;
	}

	/**
	 * @return the initialAngle (in degrees)
	 */
	public float getInitialAngleD()
	{
		return (float) (initialAngle * 180.0 / Math.PI);
	}

	/**
	 * @param initialPeakStdDev
	 *            An estimate for the peak standard deviation used to initialise the fit for all dimensions
	 */
	public void setInitialPeakStdDev(float initialPeakStdDev)
	{
		setInitialPeakStdDev0(initialPeakStdDev);
		setInitialPeakStdDev1(initialPeakStdDev);
	}

	/**
	 * @param initialPeakStdDev0
	 *            An estimate for the peak standard deviation used to initialise the fit for dimension 0
	 */
	public void setInitialPeakStdDev0(float initialPeakStdDev0)
	{
		this.initialSD0 = initialPeakStdDev0;
	}

	/**
	 * @return An estimate for the peak standard deviation used to initialise the fit for dimension 0
	 */
	public float getInitialPeakStdDev0()
	{
		return initialSD0;
	}

	/**
	 * @param initialPeakStdDev1
	 *            An estimate for the peak standard deviation used to initialise the fit for dimension 1
	 */
	public void setInitialPeakStdDev1(float initialPeakStdDev1)
	{
		this.initialSD1 = initialPeakStdDev1;
	}

	/**
	 * @return An estimate for the peak standard deviation used to initialise the fit for dimension 1
	 */
	public float getInitialPeakStdDev1()
	{
		return initialSD1;
	}

	/**
	 * @param computeDeviations
	 *            True if computing the parameter deviations
	 */
	public void setComputeDeviations(boolean computeDeviations)
	{
		this.computeDeviations = computeDeviations;
	}

	/**
	 * @return True if computing the parameter deviations
	 */
	public boolean isComputeDeviations()
	{
		return computeDeviations;
	}

	/**
	 * @return the fit solver used to fit the point spread function (PSF)
	 */
	public FitSolver getFitSolver()
	{
		return fitSolver;
	}

	/**
	 * @param fitSolver
	 *            the fit solver to use to fit the point spread function (PSF)
	 */
	public void setFitSolver(FitSolver fitSolver)
	{
		this.fitSolver = fitSolver;
	}

	/**
	 * @param fitSolver
	 *            the fit solver to use to fit the point spread function (PSF)
	 */
	public void setFitSolver(int fitSolver)
	{
		if (fitSolver >= 0 && fitSolver < FitSolver.values().length)
		{
			setFitSolver(FitSolver.values()[fitSolver]);
		}
	}

	/**
	 * @param minIterations
	 *            the minIterations to set
	 */
	public void setMinIterations(int minIterations)
	{
		this.minIterations = minIterations;
	}

	/**
	 * @return the minIterations
	 */
	public int getMinIterations()
	{
		return minIterations;
	}

	/**
	 * @param maxIterations
	 *            the maxIterations to set
	 */
	public void setMaxIterations(int maxIterations)
	{
		this.maxIterations = maxIterations;
	}

	/**
	 * @return the maxIterations
	 */
	public int getMaxIterations()
	{
		return maxIterations;
	}

	/**
	 * @param significantDigits
	 *            the significant digits for computing if the parameters have changed
	 * @see gdsc.smlm.fitting.FitCriteria
	 */
	public void setSignificantDigits(int significantDigits)
	{
		this.significantDigits = significantDigits;
	}

	/**
	 * @return the significant digits for computing if the parameters have changed
	 */
	public int getSignificantDigits()
	{
		return significantDigits;
	}

	/**
	 * @param backgroundFitting
	 *            True if fitting the peak widths
	 */
	public void setBackgroundFitting(boolean backgroundFitting)
	{
		this.backgroundFitting = backgroundFitting;
	}

	/**
	 * @return True if fitting the background
	 */
	public boolean isBackgroundFitting()
	{
		return backgroundFitting;
	}

	/**
	 * @return True if fitting an elliptical peak (with an angle parameter)
	 */
	public boolean isAngleFitting()
	{
		return (flags & GaussianFunctionFactory.FIT_ANGLE) != 0;
	}

	/**
	 * @return True if fitting the peak width in dimension 0
	 */
	public boolean isWidth0Fitting()
	{
		return (flags & GaussianFunctionFactory.FIT_X_WIDTH) != 0;
	}

	/**
	 * @return True if fitting the peak width in dimension 1
	 */
	public boolean isWidth1Fitting()
	{
		return (flags & GaussianFunctionFactory.FIT_Y_WIDTH) != 0;
	}

	/**
	 * @param fitFunction
	 *            the fitting function to use
	 */
	public void setFitFunction(int fitFunction)
	{
		if (fitFunction >= 0 && fitFunction < FitFunction.values().length)
		{
			setFitFunction(FitFunction.values()[fitFunction]);
		}
	}

	/**
	 * @param fitFunction
	 *            the fitting function to use
	 */
	public void setFitFunction(FitFunction fitFunction)
	{
		this.fitFunction = fitFunction;
		switch (fitFunction)
		{
			case CIRCULAR:
				flags = GaussianFunctionFactory.FIT_NB_CIRCLE;
				break;
			case FIXED:
				flags = GaussianFunctionFactory.FIT_NB_FIXED;
				break;
			case FREE_CIRCULAR:
				flags = GaussianFunctionFactory.FIT_NB_FREE_CIRCLE;
				break;
			case FREE:
				flags = GaussianFunctionFactory.FIT_NB_ELLIPTICAL;
				break;
			default:
				throw new RuntimeException("Unknown fitting function");
		}
	}

	/**
	 * @return the fitting function
	 */
	public FitFunction getFitFunction()
	{
		return fitFunction;
	}

	/**
	 * @param fitValidation
	 *            True if fit should be validated with {@link #validatePeak(int, float[], float[])}
	 */
	public void setFitValidation(boolean fitValidation)
	{
		this.fitValidation = fitValidation;
	}

	/**
	 * @return the fitValidation
	 */
	public boolean isFitValidation()
	{
		return fitValidation;
	}

	/**
	 * @param coordinateShift
	 *            The maximum absolute coordinate shift for a good fit
	 */
	public void setCoordinateShift(float coordinateShift)
	{
		this.coordinateShift = coordinateShift;
	}

	/**
	 * @return the coordinateShift
	 */
	public float getCoordinateShift()
	{
		return coordinateShift;
	}

	/**
	 * @param shiftFactor
	 *            The maximum absolute coordinate shift for a good fit, relative to the largest peak width
	 */
	public void setCoordinateShiftFactor(float shiftFactor)
	{
		if (shiftFactor > 0)
		{
			float widthMax = (initialSD0 > 0) ? initialSD0 : 1;
			if (initialSD1 > 0)
				widthMax = Math.max(initialSD1, widthMax);
			setCoordinateShift(shiftFactor * widthMax / 2);
		}
		else
		{
			setCoordinateShift(Float.POSITIVE_INFINITY);
		}
	}

	/**
	 * @return the coordinateShift relative the the largest peak width
	 */
	public float getCoordinateShiftFactor()
	{
		float widthMax = (initialSD0 > 0) ? initialSD0 : Gaussian2DFitter.sd2fwhm(1);
		if (initialSD1 > 0)
			widthMax = Math.max(initialSD1, widthMax);
		return coordinateShift * 2 / widthMax;
	}

	/**
	 * @param signalStrength
	 *            The signal strength. Used to determine the signal strength for a good fit (signalThreshold = noise x
	 *            signalStrength)
	 */
	public void setSignalStrength(float signalStrength)
	{
		this.signalStrength = signalStrength;
		setSignalThreshold();
	}

	/**
	 * @return the signal strength
	 */
	public float getSignalStrength()
	{
		return signalStrength;
	}

	/**
	 * @return the precision threshold. Used to determine if the peak is a good fit. Requires that the image is
	 *         calibrated
	 */
	public float getPrecisionThreshold()
	{
		return precisionThreshold;
	}

	/**
	 * @param precisionThreshold
	 *            the precisionThreshold to set
	 */
	public void setPrecisionThreshold(float precisionThreshold)
	{
		this.precisionThreshold = precisionThreshold;
	}

	/**
	 * @param noise
	 *            The image noise. Used to determine the signal strength for a good fit (signalThreshold = noise x
	 *            signalStrength)
	 */
	public void setNoise(float noise)
	{
		this.noise = noise;
		setSignalThreshold();
	}

	/**
	 * @return the image noise
	 */
	public float getNoise()
	{
		return noise;
	}

	private void setSignalThreshold()
	{
		signalThreshold = noise * signalStrength;
	}

	/**
	 * @param widthFactor
	 *            The factor difference allowed between widths for a good fit
	 */
	public void setWidthFactor(float widthFactor)
	{
		if (widthFactor > 0)
		{
			this.widthFactor = widthFactor;
		}
		else
		{
			this.widthFactor = Float.POSITIVE_INFINITY;
		}
	}

	/**
	 * @return the widthFactor
	 */
	public float getWidthFactor()
	{
		return widthFactor;
	}

	/**
	 * @return the stoppingCriteria
	 */
	public StoppingCriteria getStoppingCriteria()
	{
		return stoppingCriteria;
	}

	/**
	 * @return the gaussianFunction
	 */
	public GaussianFunction getGaussianFunction()
	{
		return gaussianFunction;
	}

	/**
	 * @return the result
	 */
	public FitStatus getValidationResult()
	{
		return result;
	}

	/**
	 * @return Data associated with the validation result
	 */
	public Object getValidationData()
	{
		return statusData;
	}

	private FitStatus setValidationResult(FitStatus newResult, Object data)
	{
		result = newResult;
		statusData = data;
		return result;
	}

	private FitStatus result;
	private Object statusData;

	/**
	 * Check peaks to see if the fit was sensible
	 * 
	 * @param nPeaks
	 *            The number of peaks
	 * @param initialParams
	 *            The initial peak parameters
	 * @param peakParams
	 *            The fitted peak parameters
	 * @return True if the fit fails the criteria
	 */
	public FitStatus validateFit(int nPeaks, float[] initialParams, float[] peakParams)
	{
		for (int n = 0; n < nPeaks; n++)
		{
			validatePeak(n, initialParams, peakParams);
			if (result != FitStatus.OK)
				break;
		}

		return result;
	}

	/**
	 * Check peak to see if the fit was sensible. Assumes a single peak.
	 * 
	 * @param initialParams
	 *            The initial peak parameters
	 * @param peakParams
	 *            The fitted peak parameters
	 * @return True if the fit fails the criteria
	 */
	public FitStatus validateFit(float[] initialParams, float[] peakParams)
	{
		return validatePeak(0, initialParams, peakParams);
	}

	/**
	 * Check peak to see if the fit was sensible
	 * 
	 * @param n
	 *            The peak number
	 * @param initialParams
	 *            The initial peak parameters
	 * @param peakParams
	 *            The fitted peak parameters
	 * @return True if the fit fails the criteria
	 */
	public FitStatus validatePeak(int n, float[] initialParams, float[] peakParams)
	{
		final int offset = n * 6;
		// Check spot movement
		float xShift = peakParams[Gaussian2DFunction.X_POSITION + offset] -
				initialParams[Gaussian2DFunction.X_POSITION + offset];
		float yShift = peakParams[Gaussian2DFunction.Y_POSITION + offset] -
				initialParams[Gaussian2DFunction.Y_POSITION + offset];
		if (Math.abs(xShift) > coordinateShift || Math.abs(yShift) > coordinateShift)
		{
			if (log != null)
			{
				log.info("Bad peak %d: Fitted coordinates moved (x=%g,y=%g)", n, xShift, yShift);
			}
			return setValidationResult(FitStatus.COORDINATES_MOVED, new float[] { xShift, yShift });
		}

		// Check signal threshold
		float signal = 0;
		if (signalStrength != 0)
		{
			// Signal = Amplitude * 2 * pi * sx * sy
			final float factor = (float) (Math.PI * 2);
			signal = factor * peakParams[Gaussian2DFunction.AMPLITUDE + offset] *
					peakParams[Gaussian2DFunction.X_SD + offset] * peakParams[Gaussian2DFunction.Y_SD + offset];

			// Compare the signal to the desired signal strength
			if (signal < signalThreshold)
			{
				if (log != null)
				{
					log.info("Bad peak %d: Insufficient signal (%gx)\n", n, signal / noise);
				}
				//System.out.printf("Bad peak %d: Insufficient signal (%gx)\n", n, signal / noise);
				return setValidationResult(FitStatus.INSUFFICIENT_SIGNAL, signal / noise);
			}
		}

		// Check widths
		float xFactor = getFactor(peakParams[Gaussian2DFunction.X_SD + offset], initialParams[Gaussian2DFunction.X_SD +
				offset]);
		float yFactor = getFactor(peakParams[Gaussian2DFunction.Y_SD + offset], initialParams[Gaussian2DFunction.Y_SD +
				offset]);
		if (xFactor > widthFactor || yFactor > widthFactor)
		{
			if (log != null)
			{
				log.info(
						"Bad peak %d: Fitted width diverged (x=%gx,y=%gx)\n",
						n,
						(peakParams[Gaussian2DFunction.X_SD + offset] > initialParams[Gaussian2DFunction.X_SD + offset]) ? xFactor
								: -xFactor,
						(peakParams[Gaussian2DFunction.Y_SD + offset] > initialParams[Gaussian2DFunction.Y_SD + offset]) ? yFactor
								: -yFactor);
			}
			return setValidationResult(FitStatus.WIDTH_DIVERGED, new float[] { xFactor, yFactor });
		}

		// Check precision
		if (precisionThreshold > 0 && nmPerPixel > 0 && gain > 0)
		{
			if (signal == 0)
			{
				final float factor = (float) (Math.PI * 2);
				signal = factor * peakParams[Gaussian2DFunction.AMPLITUDE + offset] *
						peakParams[Gaussian2DFunction.X_SD + offset] * peakParams[Gaussian2DFunction.Y_SD + offset];
			}
			float sd = (peakParams[Gaussian2DFunction.X_SD + offset] + peakParams[Gaussian2DFunction.Y_SD + offset]) * 0.5f;
			double p = PeakResult.getPrecision(nmPerPixel, nmPerPixel * sd, signal / gain, noise / gain);
			if (p > precisionThreshold)
			{
				if (log != null)
				{
					log.info("Bad peak %d: Insufficient precision (%gx)\n", n, p);
				}
				return setValidationResult(FitStatus.INSUFFICIENT_PRECISION, p);
			}
		}

		return setValidationResult(FitStatus.OK, null);
	}

	/**
	 * Get the relative change factor between f and g
	 * 
	 * @param f
	 * @param g
	 * @return
	 */
	private static float getFactor(float f, float g)
	{
		if (f > g)
			return f / g;
		return g / f;
	}

	/**
	 * @param lambda
	 *            the lambda to start the Levenberg-Marquardt fitting process
	 */
	public void setLambda(double lambda)
	{
		this.lambda = lambda;
	}

	/**
	 * @return the lambda
	 */
	public double getLambda()
	{
		return lambda;
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
	public void setComputeResiduals(boolean computeResiduals)
	{
		this.computeResiduals = computeResiduals;
	}

	/**
	 * @param duplicateDistance
	 *            The distance within which spots are considered duplicates
	 */
	public void setDuplicateDistance(final float duplicateDistance)
	{
		this.duplicateDistance = duplicateDistance;
	}

	/**
	 * @return The distance within which spots are considered duplicates
	 */
	public float getDuplicateDistance()
	{
		return duplicateDistance;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public Object clone()
	{
		try
		{
			FitConfiguration f = (FitConfiguration) super.clone();
			// Ensure the object reference is passed through.
			// All other fields are primitives and so should be copied by Object.clone().
			f.log = log;
			return f;
		}
		catch (CloneNotSupportedException e)
		{
			// Ignore
		}
		return null;
	}

	/**
	 * Ensure that the internal state of the object is initialised. This is used after deserialisation since some state
	 * is not saved but restored from other property values.
	 */
	public void initialiseState()
	{
		if (fitSolver == null)
			fitSolver = FitSolver.LVM;
		if (fitFunction == null)
			fitFunction = FitFunction.CIRCULAR;
		if (fitCriteria == null)
			fitCriteria = FitCriteria.LEAST_SQUARED_ERROR;
		if (initialSD0 == 0)
			initialSD0 = 1;
		if (initialSD1 == 0)
			initialSD1 = 1;
		setNoise(noise);
		setFitFunction(fitFunction);
	}

	/**
	 * @return the nmPerPixel
	 */
	public double getNmPerPixel()
	{
		return nmPerPixel;
	}

	/**
	 * @param nmPerPixel
	 *            the nm per pixel scale to use when evaluating a fitted peak's localisation precision
	 */
	public void setNmPerPixel(double nmPerPixel)
	{
		this.nmPerPixel = nmPerPixel;
	}

	/**
	 * @return the gain
	 */
	public float getGain()
	{
		return gain;
	}

	/**
	 * @param gain
	 *            the gain to use when evaluating a fitted peak's localisation precision
	 */
	public void setGain(float gain)
	{
		this.gain = gain;
	}

	/**
	 * @return the noise model
	 */
	public NoiseModel getNoiseModel()
	{
		return noiseModel;
	}

	/**
	 * Set the noise model for the fitted function
	 * 
	 * @param noiseModel
	 *            the noise model
	 */
	public void setNoiseModel(NoiseModel noiseModel)
	{
		this.noiseModel = noiseModel;
	}

	/**
	 * @return The function solver for the current configuration
	 */
	public FunctionSolver getFunctionSolver()
	{
		GaussianFunction gf = getGaussianFunction();
		switch (fitSolver)
		{
			case MLE:
				MaximumLikelihoodFitter fitter = new MaximumLikelihoodFitter(gf);
				fitter.setMaxEvaluations(maxIterations);
				// TODO - Configure stopping criteria ...
				return fitter;
			
			case APACHE_LVM:
				if (gf instanceof Gaussian2DFunction)
				{
					ApacheLVMFitter apacheNLinFit = new ApacheLVMFitter((Gaussian2DFunction) gf);
					apacheNLinFit.setMaxEvaluations(maxIterations);
					// TODO - Configure stopping criteria ...
					return apacheNLinFit;
				}
				// else fall through to default fitter

			case LVM_WEIGHTED:
				gf.setNoiseModel(getNoiseModel());

			case LVM:
			default:
				StoppingCriteria sc = getStoppingCriteria();
				NonLinearFit nlinfit = new NonLinearFit();
				nlinfit.setStoppingCriteria(sc);
				nlinfit.setNonLinearFunction(gf);
				nlinfit.setInitialLambda(getLambda());
				return nlinfit;
		}
	}
}