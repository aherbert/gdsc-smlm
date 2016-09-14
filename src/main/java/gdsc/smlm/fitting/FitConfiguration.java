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
 *---------------------------------------------------------------------------*/

import gdsc.smlm.fitting.nonlinear.ApacheLVMFitter;
import gdsc.smlm.fitting.nonlinear.MaximumLikelihoodFitter;
import gdsc.smlm.fitting.nonlinear.MaximumLikelihoodFitter.SearchMethod;
import gdsc.smlm.fitting.nonlinear.NonLinearFit;
import gdsc.smlm.fitting.nonlinear.StoppingCriteria;
import gdsc.smlm.fitting.nonlinear.stop.ErrorStoppingCriteria;
import gdsc.smlm.fitting.nonlinear.stop.GaussianStoppingCriteria;
import gdsc.smlm.fitting.nonlinear.stop.ParameterStoppingCriteria;
import gdsc.smlm.function.NoiseModel;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.filter.PreprocessedPeakResult;
import gdsc.core.logging.Logger;
import gdsc.core.match.FractionalAssignment;
import gdsc.core.utils.Maths;

/**
 * Specifies the fitting configuration for Gaussian fitting
 */
public class FitConfiguration implements Cloneable
{
	private FitCriteria fitCriteria = FitCriteria.LEAST_SQUARED_ERROR;
	private Logger log = null;
	private double delta = 0.0001;
	private double initialAngle = 0; // Radians
	private double initialSD0 = 1;
	private double initialSD1 = 1;
	private boolean computeDeviations = false;
	private FitSolver fitSolver = FitSolver.LVM;
	private int minIterations = 0;
	private int maxIterations = 20;
	private int significantDigits = 5;
	private FitFunction fitFunction;
	private int flags;
	private boolean backgroundFitting = true;
	private boolean notSignalFitting = false;
	private double coordinateShift = 1;
	private double shiftFactor = 1;
	private int fitRegion = 0;
	private double coordinateOffset = 0.5;
	private double signalThreshold = 0;
	private double signalStrength = 0;
	private double minPhotons = 30;
	private double precisionThreshold = 1600; // == 40nm * 40nm;
	private boolean precisionUsingBackground = false;
	private double nmPerPixel = 0;
	private double gain = 0;
	private boolean emCCD = true;
	private boolean modelCamera = false;
	private double noise = 0;
	private double minWidthFactor = 0.5;
	private double widthFactor = 2;
	private boolean fitValidation = true;
	private double lambda = 10;
	private boolean computeResiduals = true;
	private double duplicateDistance = 0.5f;
	private double bias = 0;
	private double readNoise = 0;
	private double amplification = 0;
	private int maxFunctionEvaluations = 2000;
	private SearchMethod searchMethod = SearchMethod.POWELL_BOUNDED; // Best for noisy data since gradients are unstable
	private boolean gradientLineMinimisation = false;
	private double relativeThreshold = 1e-6;
	private double absoluteThreshold = 1e-16;

	private StoppingCriteria stoppingCriteria = null;
	private GaussianFunction gaussianFunction = null;
	private NoiseModel noiseModel = null;
	private FunctionSolver functionSolver = null;

	private double[] peakShiftFactors = null;
	private DynamicPeakResult dynamicPeakResult = new DynamicPeakResult();

	/**
	 * Sets the peak shift factors. This is used to adjust the coordinate shift limit for individual peaks in the
	 * {@link #validatePeak(int, double[], double[])} function. For example if fitting multiple peaks and the validation
	 * of each position should be adjusted based on the uncertainty of the estimated parameters.
	 *
	 * @param peakShiftFactors
	 *            the new peak shift factors
	 */
	public void setPeakShiftFactors(double[] peakShiftFactors)
	{
		this.peakShiftFactors = peakShiftFactors;
	}

	/**
	 * Gets the peak shift factors
	 *
	 * @return a clone copy of the peak shift factors
	 */
	public double[] getPeakShiftFactors()
	{
		return (this.peakShiftFactors == null) ? null : peakShiftFactors.clone();
	}

	/**
	 * Gets the peak shift factor for the specified peak.
	 *
	 * @param peak
	 *            the peak
	 * @return the peak shift factor (1 if not explicitly set using {@link #setPeakShiftFactors(double[])}
	 */
	public double getPeakShiftFactor(int peak)
	{
		if (peakShiftFactors == null || peakShiftFactors.length <= peak)
			return 1;
		return peakShiftFactors[peak];
	}

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
	public void initialise(int npeaks, int maxx, double[] params)
	{
		// Check if the Gaussian function is invalid
		if (gaussianFunction != null &&
				(gaussianFunction.getNPeaks() != npeaks || gaussianFunction.getDimensions()[0] != maxx))
		{
			invalidateGaussianFunction();
		}
		if (gaussianFunction == null)
		{
			invalidateFunctionSolver();
			gaussianFunction = createGaussianFunction(npeaks, maxx, params);
		}
		if (stoppingCriteria == null)
		{
			invalidateFunctionSolver();
			stoppingCriteria = createStoppingCriteria(gaussianFunction, params);
		}
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
	public StoppingCriteria createStoppingCriteria(GaussianFunction func, double[] params)
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
	protected void addPeakRestrictions(GaussianFunction func, GaussianStoppingCriteria sc, double[] params)
	{
		// TODO - Check if it is worth using this to stop fitting early or simply do it at the end.
		sc.setMinimumSignal(0);
		// We could also set min and max positions and widths
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
	public GaussianFunction createGaussianFunction(int npeaks, int maxx, double[] params)
	{
		final int flags = getFunctionFlags();

		GaussianFunction f = GaussianFunctionFactory.create2D(npeaks, maxx, flags);
		//f.initialise(params);
		return f;
	}

	/**
	 * Gets the function flags used for the GaussianFunctionFactory.
	 *
	 * @return the function flags
	 */
	public int getFunctionFlags()
	{
		int flags = this.flags;

		if (isBackgroundFitting())
			flags |= GaussianFunctionFactory.FIT_BACKGROUND;
		if (isNotSignalFitting())
			// Remove signal fitting (on by default)
			flags &= ~GaussianFunctionFactory.FIT_SIGNAL;
		return flags;
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
		invalidateStoppingCriteria();
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
		invalidateStoppingCriteria();
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
	public void setInitialAngle(double initialAngle)
	{
		this.initialAngle = initialAngle;
	}

	/**
	 * @param initialAngle
	 *            the initialAngle to set (in degrees)
	 */
	public void setInitialAngleD(double initialAngle)
	{
		this.initialAngle = (double) (initialAngle * Math.PI / 180.0);
	}

	/**
	 * @return the initialAngle (in radians)
	 */
	public double getInitialAngle()
	{
		return initialAngle;
	}

	/**
	 * @return the initialAngle (in degrees)
	 */
	public double getInitialAngleD()
	{
		return (double) (initialAngle * 180.0 / Math.PI);
	}

	/**
	 * @param initialPeakStdDev
	 *            An estimate for the peak standard deviation used to initialise the fit for all dimensions
	 */
	public void setInitialPeakStdDev(double initialPeakStdDev)
	{
		setInitialPeakStdDev0(initialPeakStdDev);
		setInitialPeakStdDev1(initialPeakStdDev);
	}

	/**
	 * Set an estimate for the peak standard deviation used to initialise the fit for dimension 0
	 * <p>
	 * Setting this will update the value in {@link #getCoordinateShift()}
	 * 
	 * @param initialPeakStdDev0
	 *            An estimate for the peak standard deviation used to initialise the fit for dimension 0
	 */
	public void setInitialPeakStdDev0(double initialPeakStdDev0)
	{
		this.initialSD0 = initialPeakStdDev0;
		updateCoordinateShift();
	}

	/**
	 * @return An estimate for the peak standard deviation used to initialise the fit for dimension 0
	 */
	public double getInitialPeakStdDev0()
	{
		return initialSD0;
	}

	/**
	 * Set an estimate for the peak standard deviation used to initialise the fit for dimension 1
	 * <p>
	 * Setting this will update the value in {@link #getCoordinateShift()}
	 * 
	 * @param initialPeakStdDev1
	 *            An estimate for the peak standard deviation used to initialise the fit for dimension 1
	 */
	public void setInitialPeakStdDev1(double initialPeakStdDev1)
	{
		this.initialSD1 = initialPeakStdDev1;
		updateCoordinateShift();
	}

	/**
	 * @return An estimate for the peak standard deviation used to initialise the fit for dimension 1
	 */
	public double getInitialPeakStdDev1()
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
		invalidateStoppingCriteria();
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
		invalidateStoppingCriteria();
		invalidateFunctionSolver();
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
		invalidateStoppingCriteria();
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
	 *            True if fitting the background
	 */
	public void setBackgroundFitting(boolean backgroundFitting)
	{
		invalidateGaussianFunction();
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
	 * Use this to turn off fitting of the signal. This should be used with caution. The setting only applies to fixed
	 * width fitting and can be used to benchmark position accuracy when fitting signals of known strength.
	 * 
	 * @param noSignalFitting
	 *            True if not fitting the signal
	 */
	public void setNotSignalFitting(boolean noSignalFitting)
	{
		invalidateGaussianFunction();
		this.notSignalFitting = noSignalFitting;
	}

	/**
	 * @return True if not fitting the signal. The setting only applies to fixed width fitting
	 */
	public boolean isNotSignalFitting()
	{
		return notSignalFitting;
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
		invalidateGaussianFunction();
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
	 *            True if fit should be validated with {@link #validatePeak(int, double[], double[])}
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
	 * Set the maximum absolute coordinate shift for a good fit. This is also set when calling
	 * {@link #setCoordinateShiftFactor(double)} or any of the standard deviations, e.g.
	 * {@link #setInitialPeakStdDev(double)}. If these are set then the coordinate shift will change.
	 * 
	 * @param coordinateShift
	 *            The maximum absolute coordinate shift for a good fit
	 */
	public void setCoordinateShift(double coordinateShift)
	{
		this.coordinateShift = coordinateShift;
	}

	/**
	 * @return The maximum absolute coordinate shift for a good fit
	 */
	public double getCoordinateShift()
	{
		return coordinateShift;
	}

	/**
	 * Set the maximum absolute coordinate shift for a good fit, relative to the largest peak width.
	 * Set to zero to disable.
	 * <p>
	 * Setting this will update the value in {@link #getCoordinateShift()}
	 * 
	 * @param shiftFactor
	 *            The maximum absolute coordinate shift for a good fit, relative to the largest peak width
	 */
	public void setCoordinateShiftFactor(double shiftFactor)
	{
		this.shiftFactor = (shiftFactor > 0) ? shiftFactor : 0;
		// Reset to zero for disabled values 
		this.shiftFactor = getCoordinateShiftFactor();
		updateCoordinateShift();
	}

	private void updateCoordinateShift()
	{
		if (shiftFactor > 0)
		{
			final double widthMax = Maths.max(initialSD0, initialSD1);
			if (widthMax > 0)
			{
				setCoordinateShift(shiftFactor * widthMax);
				return;
			}
		}
		setCoordinateShift(Double.POSITIVE_INFINITY);
	}

	/**
	 * @return the coordinateShift relative the the largest peak width
	 */
	public double getCoordinateShiftFactor()
	{
		if (shiftFactor == Double.POSITIVE_INFINITY)
			return 0;
		return shiftFactor;
	}

	/**
	 * Gets the maximum distance a peak will be allowed to shift for the specified peak.
	 * <p>
	 * This is constructed using the coordinate shift factor
	 *
	 * @param peak
	 *            the peak
	 * @return the peak shift
	 */
	public double getMaxShift(int peak)
	{
		return coordinateShift * getPeakShiftFactor(peak);
	}

	/**
	 * @return the size of the fit region used for validation
	 */
	public int getFitRegion()
	{
		return fitRegion;
	}

	/**
	 * Set the size of the fit region (N). Any coordinate outside the region will fail fit validation (see
	 * {@link #validatePeak(int, double[], double[])}). Set to zero to disable.
	 * <p>
	 * Note: it is assumed that the coordinates of the peak are relative to the fit region of size NxN. Coordinates are
	 * offset by the amount defined by {@link #setCoordinateOffset(double)}.
	 * 
	 * @param fitRegion
	 *            the size of the fit region
	 */
	public void setFitRegion(int fitRegion)
	{
		this.fitRegion = Math.max(0, fitRegion);
	}

	/**
	 * @return the coordinate offset when validating the coordinates are within the fit window
	 */
	public double getCoordinateOffset()
	{
		return coordinateOffset;
	}

	/**
	 * @param coordinateOffset
	 *            the coordinate offset when validating the coordinates are within the fit window
	 */
	public void setCoordinateOffset(double coordinateOffset)
	{
		this.coordinateOffset = coordinateOffset;
	}

	/**
	 * @param signalStrength
	 *            The signal strength. Used to determine the signal strength for a good fit (signalThreshold =
	 *            max(gain x minPhotons, noise x signalStrength).
	 */
	public void setSignalStrength(double signalStrength)
	{
		this.signalStrength = signalStrength;
		setSignalThreshold();
	}

	/**
	 * @return the signal strength
	 */
	public double getSignalStrength()
	{
		return signalStrength;
	}

	/**
	 * @return The minimum number of photons
	 */
	public double getMinPhotons()
	{
		return minPhotons;
	}

	/**
	 * @param minPhotons
	 *            The minimum number of photons. Used to determine the signal strength for a good fit (signalThreshold =
	 *            max(gain x minPhotons, noise x signalStrength).
	 */
	public void setMinPhotons(double minPhotons)
	{
		this.minPhotons = minPhotons;
		setSignalThreshold();
	}

	/**
	 * @return the precision threshold. Used to determine if the peak is a good fit. Requires that the image is
	 *         calibrated
	 */
	public double getPrecisionThreshold()
	{
		return (precisionThreshold > 0) ? Math.sqrt(precisionThreshold) : 0;
	}

	/**
	 * @param precisionThreshold
	 *            the precisionThreshold to set
	 */
	public void setPrecisionThreshold(double precisionThreshold)
	{
		// Store the squared threshold
		this.precisionThreshold = precisionThreshold * precisionThreshold;
	}

	/**
	 * @return True if calculating the precision using the fitted background
	 */
	public boolean isPrecisionUsingBackground()
	{
		return precisionUsingBackground;
	}

	/**
	 * Set to true to calculate the precision using the fitted background. Set to false to use the configured noise
	 * (which may not be reflective of the noise at the fit location). Using false will be consistent with results
	 * analysis performed using a global estimated noise.
	 * 
	 * @param precisionUsingBackground
	 *            True if calculating the precision using the fitted background
	 */
	public void setPrecisionUsingBackground(boolean precisionUsingBackground)
	{
		this.precisionUsingBackground = precisionUsingBackground;
	}

	/**
	 * @param noise
	 *            The image noise. Used to determine the signal strength for a good fit (signalThreshold =
	 *            max(gain x minPhotons, noise x signalStrength).
	 */
	public void setNoise(double noise)
	{
		this.noise = noise;
		setSignalThreshold();
	}

	/**
	 * @return the image noise
	 */
	public double getNoise()
	{
		return noise;
	}

	private void setSignalThreshold()
	{
		signalThreshold = Math.max(noise * signalStrength, minPhotons * gain);
	}

	/**
	 * @param widthFactor
	 *            The factor difference allowed between widths for a good fit
	 */
	public void setWidthFactor(double widthFactor)
	{
		if (widthFactor > 1)
		{
			this.widthFactor = widthFactor;
		}
		else
		{
			this.widthFactor = Double.POSITIVE_INFINITY;
		}
	}

	/**
	 * @return the widthFactor (or zero if not configured)
	 */
	public double getWidthFactor()
	{
		return (widthFactor == Double.POSITIVE_INFINITY) ? 0 : widthFactor;
	}

	/**
	 * @param minWidthFactor
	 *            The minimum factor difference allowed between widths for a good fit
	 */
	public void setMinWidthFactor(double minWidthFactor)
	{
		if (minWidthFactor < 1 && minWidthFactor > 0)
		{
			this.minWidthFactor = minWidthFactor;
		}
		else
		{
			this.minWidthFactor = 0;
		}
	}

	/**
	 * @return the minWidthFactor
	 */
	public double getMinWidthFactor()
	{
		return minWidthFactor;
	}

	/**
	 * @return the stoppingCriteria
	 */
	public StoppingCriteria getStoppingCriteria()
	{
		return stoppingCriteria;
	}

	/**
	 * Call this when a property changes that will change the stopping criteria
	 */
	private void invalidateStoppingCriteria()
	{
		stoppingCriteria = null;
	}

	/**
	 * @return the gaussianFunction
	 */
	public GaussianFunction getGaussianFunction()
	{
		return gaussianFunction;
	}

	/**
	 * Call this when a property changes that will change the Gaussian function
	 */
	private void invalidateGaussianFunction()
	{
		gaussianFunction = null;
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
	 * @param params
	 *            The fitted peak parameters
	 * @return True if the fit fails the criteria
	 */
	public FitStatus validateFit(int nPeaks, double[] initialParams, double[] params)
	{
		for (int n = 0; n < nPeaks; n++)
		{
			validatePeak(n, initialParams, params);
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
	 * @param params
	 *            The fitted peak parameters
	 * @return True if the fit fails the criteria
	 */
	public FitStatus validateFit(double[] initialParams, double[] params)
	{
		return validatePeak(0, initialParams, params);
	}

	/**
	 * Check peak to see if the fit was sensible
	 * 
	 * @param n
	 *            The peak number
	 * @param initialParams
	 *            The initial peak parameters
	 * @param params
	 *            The fitted peak parameters
	 * @return True if the fit fails the criteria
	 */
	public FitStatus validatePeak(int n, double[] initialParams, double[] params)
	{
		final int offset = n * 6;
		// Check spot movement
		final double xShift = params[Gaussian2DFunction.X_POSITION + offset] -
				initialParams[Gaussian2DFunction.X_POSITION + offset];
		final double yShift = params[Gaussian2DFunction.Y_POSITION + offset] -
				initialParams[Gaussian2DFunction.Y_POSITION + offset];
		final double maxShift = getMaxShift(n);
		if (Math.abs(xShift) > maxShift || Math.abs(yShift) > maxShift)
		{
			if (log != null)
			{
				log.info("Bad peak %d: Fitted coordinates moved (x=%g,y=%g) > %g", n, xShift, yShift, maxShift);
			}
			return setValidationResult(FitStatus.COORDINATES_MOVED, new double[] { xShift, yShift });
		}

		// Check if outside the fit window.
		// TODO - Make this configurable per peak. At the moment we only use this in BenchmarkSpotFit where 
		// additional peaks will be neighbours. In the future we may want to control this better.
		if (fitRegion != 0 && n == 0)
		{
			final double x = params[Gaussian2DFunction.X_POSITION + offset] + coordinateOffset;
			final double y = params[Gaussian2DFunction.Y_POSITION + offset] + coordinateOffset;
			if (x <= 0 || x >= fitRegion || y <= 0 || y >= fitRegion)
			{
				if (log != null)
				{
					log.info("Bad peak %d: Coordinates outside fit region (x=%g,y=%g) <> %d", n, x, y, fitRegion);
				}
				return setValidationResult(FitStatus.OUTSIDE_FIT_REGION, new double[] { x, y, fitRegion });
			}
		}

		// Check signal threshold
		final double signal = params[Gaussian2DFunction.SIGNAL + offset];
		// Compare the signal to the desired signal strength
		if (signal < signalThreshold)
		{
			if (log != null)
			{
				log.info("Bad peak %d: Insufficient signal %g (SNR=%g)\n", n, signal / ((gain > 0) ? gain : 1),
						signal / noise);
			}
			//System.out.printf("Bad peak %d: Insufficient signal (%gx)\n", n, signal / noise);
			return setValidationResult(FitStatus.INSUFFICIENT_SIGNAL, signal);
		}

		// Check widths
		if (isWidth0Fitting())
		{
			boolean badWidth = false;
			double xFactor = 0, yFactor = 0;

			xFactor = params[Gaussian2DFunction.X_SD + offset] / initialParams[Gaussian2DFunction.X_SD + offset];
			badWidth = (xFactor > widthFactor || xFactor < minWidthFactor);

			// Always do this (even if badWidth=true) since we need the factor for the return value
			if (isWidth1Fitting())
			{
				yFactor = params[Gaussian2DFunction.Y_SD + offset] / initialParams[Gaussian2DFunction.Y_SD + offset];
				badWidth = (yFactor > widthFactor || yFactor < minWidthFactor);
			}
			else
			{
				yFactor = xFactor;
			}

			if (badWidth)
			{
				if (log != null)
				{
					log.info("Bad peak %d: Fitted width diverged (x=%gx,y=%gx)\n", n, xFactor, yFactor);
				}
				return setValidationResult(FitStatus.WIDTH_DIVERGED, new double[] { xFactor, yFactor });
			}
		}

		// Check precision
		if (precisionThreshold > 0 && nmPerPixel > 0 && gain > 0)
		{
			final double sd = (params[Gaussian2DFunction.X_SD + offset] + params[Gaussian2DFunction.Y_SD + offset]) *
					0.5;
			final double variance = getVariance(params[Gaussian2DFunction.BACKGROUND], signal, sd);

			if (variance > precisionThreshold)
			{
				final double precision = Math.sqrt(variance);
				if (log != null)
				{
					log.info("Bad peak %d: Insufficient precision (%gx)\n", n, precision);
				}
				return setValidationResult(FitStatus.INSUFFICIENT_PRECISION, precision);
			}
		}

		return setValidationResult(FitStatus.OK, null);
	}

	/**
	 * Get the localisation variance for fitting a spot with the specified parameters given the configuration (fit
	 * solver, precision using background, gain, nm per pixel)
	 * 
	 * @param background
	 *            The background
	 * @param signal
	 *            The signal (in ADUs)
	 * @param sd
	 *            The spot standard deviation
	 * @return The localisation variance
	 */
	public double getVariance(double background, final double signal, final double sd)
	{
		double variance = 0;
		// We can calculate the precision using the estimated noise for the image or using the expected number
		// of background photons at the location.
		if (precisionUsingBackground)
		{
			// Check using the formula which uses the estimated background.
			// This allows for better filtering when the background is variable, e.g. when imaging cells.
			if (fitSolver == FitSolver.MLE)
			{
				try
				{
					// This may be slow due to the integration required within the formula.
					variance = PeakResult.getMLVarianceX(nmPerPixel, nmPerPixel * sd, signal / gain,
							Math.max(0, background - bias) / gain, emCCD);
				}
				catch (Exception e)
				{
					// Catch all exceptions. They are likely to be a TooManyIterationsException and other
					// problems with the integration
					variance = PeakResult.getVarianceX(nmPerPixel, nmPerPixel * sd, signal / gain,
							Math.max(0, background - bias) / gain, emCCD);
				}
			}
			else
			{
				variance = PeakResult.getVarianceX(nmPerPixel, nmPerPixel * sd, signal / gain,
						Math.max(0, background - bias) / gain, emCCD);
			}
		}
		else
		{
			if (fitSolver == FitSolver.MLE)
			{
				try
				{
					// This may be slow due to the integration required within the formula.
					variance = PeakResult.getMLVariance(nmPerPixel, nmPerPixel * sd, signal / gain, noise / gain,
							emCCD);
				}
				catch (Exception e)
				{
					// Catch all exceptions. They are likely to be a TooManyIterationsException and other
					// problems with the integration
					variance = PeakResult.getVariance(nmPerPixel, nmPerPixel * sd, signal / gain, noise / gain, emCCD);
				}
			}
			else
			{
				variance = PeakResult.getVariance(nmPerPixel, nmPerPixel * sd, signal / gain, noise / gain, emCCD);
			}
		}
		return variance;
	}

	/**
	 * An object that can return the results in a formatted state for the multi-path filter
	 */
	private class DynamicPeakResult implements PreprocessedPeakResult
	{
		int offset;
		double[] initialParams;
		double[] params;

		/**
		 * @param n
		 *            The peak number
		 * @param initialParams
		 *            The initial peak parameters
		 * @param params
		 *            The fitted peak parameters
		 */
		void setParameters(int n, double[] initialParams, double[] params)
		{
			offset = n * 6;

		}

		public int getFrame()
		{
			// Not implemented
			return 0;
		}

		public float getSignal()
		{
			return (float) params[Gaussian2DFunction.SIGNAL + offset];
		}

		public float getPhotons()
		{
			return (float) (params[Gaussian2DFunction.SIGNAL + offset] / FitConfiguration.this.gain);
		}

		public float getSNR()
		{
			return (float) (params[Gaussian2DFunction.SIGNAL + offset] / FitConfiguration.this.noise);
		}

		public float getNoise()
		{
			return (float) FitConfiguration.this.noise;
		}

		public double getLocationVariance()
		{
			final double sd = (params[Gaussian2DFunction.X_SD + offset] + params[Gaussian2DFunction.Y_SD + offset]) *
					0.5;
			return FitConfiguration.this.getVariance(params[Gaussian2DFunction.BACKGROUND],
					params[Gaussian2DFunction.SIGNAL + offset], sd);
		}

		public float getSD()
		{
			return (float) ((params[Gaussian2DFunction.X_SD + offset] + params[Gaussian2DFunction.Y_SD + offset]) *
					0.5);
		}

		public float getBackground()
		{
			return (float) params[Gaussian2DFunction.BACKGROUND];
		}

		public float getAmplitude()
		{
			return (float) (params[Gaussian2DFunction.SIGNAL + offset] /
					(2 * Math.PI * params[Gaussian2DFunction.X_SD + offset] +
							params[Gaussian2DFunction.Y_SD + offset]));
		}

		public float getAngle()
		{
			return (float) params[Gaussian2DFunction.ANGLE + offset];
		}

		public float getX()
		{
			return (float) params[Gaussian2DFunction.X_POSITION + offset];
		}

		public float getY()
		{
			return (float) params[Gaussian2DFunction.Y_POSITION + offset];
		}

		public float getXRelativeShift2()
		{
			final double d = (params[Gaussian2DFunction.X_POSITION + offset] -
					initialParams[Gaussian2DFunction.X_POSITION + offset]) / initialParams[Gaussian2DFunction.X_SD + offset];
			return (float) (d * d);
		}

		public float getYRelativeShift2()
		{
			final double d = (params[Gaussian2DFunction.Y_POSITION + offset] -
					initialParams[Gaussian2DFunction.Y_POSITION + offset]) / initialParams[Gaussian2DFunction.Y_SD + offset];
			return (float) (d * d);
		}

		public float getXSD()
		{
			return (float) params[Gaussian2DFunction.X_SD + offset];
		}

		public float getYSD()
		{
			return (float) params[Gaussian2DFunction.Y_SD + offset];
		}

		public float getXSDFactor()
		{
			return (float) (params[Gaussian2DFunction.X_SD + offset] / initialParams[Gaussian2DFunction.X_SD + offset]);
		}

		public float getYSDFactor()
		{
			return (float) (params[Gaussian2DFunction.Y_SD + offset] / initialParams[Gaussian2DFunction.Y_SD + offset]);
		}

		public boolean isExistingResult()
		{
			return false;
		}

		public boolean isNewResult()
		{
			return false;
		}

		public FractionalAssignment[] getAssignments(int predictedId)
		{
			return null;
		}
	}

	/**
	 * Create an object that can return the results in a formatted state for the multi-path filter
	 * 
	 * @param n
	 *            The peak number
	 * @param initialParams
	 *            The initial peak parameters
	 * @param params
	 *            The fitted peak parameters
	 * @return A preprocessed peak result
	 */
	public PreprocessedPeakResult createPreprocessedPeakResult(int n, double[] initialParams, double[] params)
	{
		dynamicPeakResult.setParameters(n, initialParams, params);
		return dynamicPeakResult;
	}

	/**
	 * @param lambda
	 *            the lambda to start the Levenberg-Marquardt fitting process
	 */
	public void setLambda(double lambda)
	{
		invalidateFunctionSolver();
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
	public void setDuplicateDistance(final double duplicateDistance)
	{
		this.duplicateDistance = duplicateDistance;
	}

	/**
	 * @return The distance within which spots are considered duplicates
	 */
	public double getDuplicateDistance()
	{
		return duplicateDistance;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public FitConfiguration clone()
	{
		try
		{
			FitConfiguration f = (FitConfiguration) super.clone();
			// Ensure the object reference is passed through.
			// All other fields are primitives and so should be copied by Object.clone().
			f.log = log;
			// Reset instance specific objects
			f.stoppingCriteria = null;
			f.gaussianFunction = null;
			f.noiseModel = null;
			f.functionSolver = null;
			f.setValidationResult(null, null);
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
		if (searchMethod == null)
			searchMethod = SearchMethod.POWELL;
		if (maxFunctionEvaluations == 0)
			maxFunctionEvaluations = 1000;
		if (initialSD0 == 0)
			initialSD0 = 1;
		if (initialSD1 == 0)
			initialSD1 = 1;
		if (shiftFactor == 0)
			setCoordinateShiftFactor(1);
		setNoise(noise);
		setFitFunction(fitFunction);
		invalidateFunctionSolver();
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
	public double getGain()
	{
		return gain;
	}

	/**
	 * @param gain
	 *            the gain to use when evaluating a fitted peak's localisation precision. Also used to determine the
	 *            signal threshold (signalThreshold = max(gain x minPhotons, noise x signalStrength)
	 */
	public void setGain(double gain)
	{
		invalidateFunctionSolver();
		this.gain = gain;
		setSignalThreshold();
	}

	/**
	 * @return True if using an EM-CCD camera
	 */
	public boolean isEmCCD()
	{
		return emCCD;
	}

	/**
	 * Specify if an EM-CCD camera is used. This is relevant when validating results using the localisation precision.
	 * 
	 * @param emCCD
	 *            Set to true if using an EM-CCD camera
	 */
	public void setEmCCD(boolean emCCD)
	{
		this.emCCD = emCCD;
	}

	/**
	 * @return True if modelling the camera noise during maximum likelihood fitting
	 */
	public boolean isModelCamera()
	{
		return modelCamera;
	}

	/**
	 * Specify if the camera noise should be modelled during maximum likelihood fitting. If true then the read noise
	 * must be set. If the EmCCD property is true then the gain must also be set.
	 * 
	 * @param modelCamera
	 *            Set to true to model the camera
	 */
	public void setModelCamera(boolean modelCamera)
	{
		invalidateFunctionSolver();
		this.modelCamera = modelCamera;
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
		invalidateFunctionSolver();
		this.noiseModel = noiseModel;
	}

	/**
	 * @return the camera bias (used for maximum likelihood estimation)
	 */
	public double getBias()
	{
		return bias;
	}

	/**
	 * @param bias
	 *            the camera bias (used for maximum likelihood estimation to evaluate the correct value of the
	 *            observed count)
	 */
	public void setBias(double bias)
	{
		this.bias = bias;
	}

	/**
	 * @return the camera read noise (used for maximum likelihood estimation)
	 */
	public double getReadNoise()
	{
		return readNoise;
	}

	/**
	 * @param readNoise
	 *            the camera read noise (used for maximum likelihood estimation)
	 */
	public void setReadNoise(double readNoise)
	{
		invalidateFunctionSolver();
		this.readNoise = readNoise;
	}

	/**
	 * @return The amplification [ADUs/electron] (used for maximum likelihood estimation)
	 */
	public double getAmplification()
	{
		return amplification;
	}

	/**
	 * @param amplification
	 *            The amplification [ADUs/electron] (used for maximum likelihood estimation)
	 */
	public void setAmplification(double amplification)
	{
		invalidateFunctionSolver();
		this.amplification = amplification;
	}

	/**
	 * @return Set to true if the bias should be removed from the data before fitting, e.g. for maximum likelihood
	 *         estimation.
	 */
	public boolean isRemoveBiasBeforeFitting()
	{
		return fitSolver == FitSolver.MLE;
	}

	/**
	 * @return the maximum number of function evaluations for the Maximum Likelihood Estimator
	 */
	public int getMaxFunctionEvaluations()
	{
		return maxFunctionEvaluations;
	}

	/**
	 * @param maxFunctionEvaluations
	 *            the maximum number of function evaluations for the Maximum Likelihood Estimator
	 */
	public void setMaxFunctionEvaluations(int maxFunctionEvaluations)
	{
		invalidateFunctionSolver();
		this.maxFunctionEvaluations = maxFunctionEvaluations;
	}

	/**
	 * @return the search for the Maximum Likelihood Estimator
	 */
	public SearchMethod getSearchMethod()
	{
		return searchMethod;
	}

	/**
	 * @param searchMethod
	 *            the search for the Maximum Likelihood Estimator
	 */
	public void setSearchMethod(int searchMethod)
	{
		if (searchMethod >= 0 && searchMethod < SearchMethod.values().length)
		{
			setSearchMethod(SearchMethod.values()[searchMethod]);
		}
	}

	/**
	 * @param searchMethod
	 *            the search for the Maximum Likelihood Estimator
	 */
	public void setSearchMethod(SearchMethod searchMethod)
	{
		invalidateFunctionSolver();
		this.searchMethod = searchMethod;
	}

	/**
	 * This setting applies to the conjugate gradient method of the Maximum Likelihood Estimator
	 * 
	 * @return the gradientLineMinimisation True if using the gradient for line minimisation
	 */
	public boolean isGradientLineMinimisation()
	{
		return gradientLineMinimisation;
	}

	/**
	 * This setting applies to the conjugate gradient method of the Maximum Likelihood Estimator
	 * 
	 * @param gradientLineMinimisation
	 *            Set to true to use the gradient for line minimisation
	 */
	public void setGradientLineMinimisation(boolean gradientLineMinimisation)
	{
		invalidateFunctionSolver();
		this.gradientLineMinimisation = gradientLineMinimisation;
	}

	/**
	 * @return the relative threshold for convergence in the Maximum Likelihood Estimator
	 */
	public double getRelativeThreshold()
	{
		return relativeThreshold;
	}

	/**
	 * @param relativeThreshold
	 *            the relative threshold for convergence in the Maximum Likelihood Estimator
	 */
	public void setRelativeThreshold(double relativeThreshold)
	{
		invalidateFunctionSolver();
		this.relativeThreshold = relativeThreshold;
	}

	/**
	 * @return the absolute threshold for convergence in the Maximum Likelihood Estimator
	 */
	public double getAbsoluteThreshold()
	{
		return absoluteThreshold;
	}

	/**
	 * @param absoluteThreshold
	 *            the absolute threshold for convergence in the Maximum Likelihood Estimator
	 */
	public void setAbsoluteThreshold(double absoluteThreshold)
	{
		invalidateFunctionSolver();
		this.absoluteThreshold = absoluteThreshold;
	}

	/**
	 * @return The function solver for the current configuration
	 */
	public FunctionSolver getFunctionSolver()
	{
		if (functionSolver == null)
			functionSolver = createFunctionSolver();
		return functionSolver;
	}

	/**
	 * Call this when a property changes that will change the function solver
	 */
	private void invalidateFunctionSolver()
	{
		functionSolver = null;
	}

	private FunctionSolver createFunctionSolver()
	{
		switch (fitSolver)
		{
			case MLE:
				// Only the Poisson likelihood function supports gradients
				if (searchMethod.usesGradients() && modelCamera)
				{
					throw new IllegalArgumentException(String.format(
							"The derivative based search method '%s' can only be used with the " +
									"'%s' likelihood function, i.e. no model camera noise",
							searchMethod, MaximumLikelihoodFitter.LikelihoodFunction.POISSON));
				}

				MaximumLikelihoodFitter fitter = new MaximumLikelihoodFitter(gaussianFunction);
				fitter.setRelativeThreshold(relativeThreshold);
				fitter.setAbsoluteThreshold(absoluteThreshold);
				fitter.setMaxEvaluations(maxFunctionEvaluations);
				fitter.setMaxIterations(maxIterations);
				fitter.setSearchMethod(searchMethod);
				fitter.setGradientLineMinimisation(gradientLineMinimisation);

				// Specify the likelihood function to use
				if (modelCamera)
				{
					// Set the camera read noise
					fitter.setSigma(readNoise);

					if (emCCD)
					{
						// EMCCD = Poisson+Gamma+Gaussian
						fitter.setLikelihoodFunction(MaximumLikelihoodFitter.LikelihoodFunction.POISSON_GAMMA_GAUSSIAN);
					}
					else
					{
						// CCD = Poisson+Gaussian
						fitter.setLikelihoodFunction(MaximumLikelihoodFitter.LikelihoodFunction.POISSON_GAUSSIAN);
					}
				}
				else
				{
					fitter.setLikelihoodFunction(MaximumLikelihoodFitter.LikelihoodFunction.POISSON);
				}

				// All models use the amplification gain (i.e. how many ADUs/electron)
				fitter.setAlpha(1.0 / amplification);

				// TODO - Configure stopping criteria ...

				return fitter;

			case LVM_QUASI_NEWTON:
				if (gaussianFunction instanceof Gaussian2DFunction)
				{
					ApacheLVMFitter apacheNLinFit = new ApacheLVMFitter((Gaussian2DFunction) gaussianFunction);
					apacheNLinFit.setMaxEvaluations(maxIterations);
					// TODO - Configure stopping criteria ...
					return apacheNLinFit;
				}
				// else fall through to default fitter

			case LVM_WEIGHTED:
				// Do not set the noise model on the quasi-newton fall-through case
				if (fitSolver != FitSolver.LVM_QUASI_NEWTON)
					gaussianFunction.setNoiseModel(getNoiseModel());

			case LVM:
			default:
				NonLinearFit nlinfit = new NonLinearFit(gaussianFunction, getStoppingCriteria());
				nlinfit.setInitialLambda(getLambda());
				return nlinfit;
		}
	}
}