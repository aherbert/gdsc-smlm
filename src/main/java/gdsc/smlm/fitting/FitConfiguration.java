package gdsc.smlm.fitting;

import gdsc.core.logging.Logger;
import gdsc.core.match.FractionalAssignment;
import gdsc.core.utils.Maths;
import gdsc.core.utils.NotImplementedException;

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
import gdsc.smlm.fitting.nonlinear.BaseFunctionSolver;
import gdsc.smlm.fitting.nonlinear.BoundedNonLinearFit;
import gdsc.smlm.fitting.nonlinear.MaximumLikelihoodFitter;
import gdsc.smlm.fitting.nonlinear.MaximumLikelihoodFitter.SearchMethod;
import gdsc.smlm.fitting.nonlinear.NonLinearFit;
import gdsc.smlm.fitting.nonlinear.ParameterBounds;
import gdsc.smlm.fitting.nonlinear.StoppingCriteria;
import gdsc.smlm.fitting.nonlinear.stop.ErrorStoppingCriteria;
import gdsc.smlm.fitting.nonlinear.stop.GaussianStoppingCriteria;
import gdsc.smlm.fitting.nonlinear.stop.ParameterStoppingCriteria;
import gdsc.smlm.function.NoiseModel;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.filter.BasePreprocessedPeakResult;
import gdsc.smlm.results.filter.BasePreprocessedPeakResult.ResultType;
import gdsc.smlm.results.filter.DirectFilter;
import gdsc.smlm.results.filter.FilterType;
import gdsc.smlm.results.filter.IDirectFilter;
import gdsc.smlm.results.filter.MultiFilter;
import gdsc.smlm.results.filter.MultiFilter2;
import gdsc.smlm.results.filter.PreprocessedPeakResult;

/**
 * Specifies the fitting configuration for Gaussian fitting
 */
public class FitConfiguration implements Cloneable, IDirectFilter
{
	private FitCriteria fitCriteria;
	private Logger log = null;
	private double delta = 0.0001;
	private double initialAngle = 0; // Radians
	private double initialSD0 = 1;
	private double initialSD1 = 1;
	private boolean computeDeviations = false;
	private FitSolver fitSolver;
	private int minIterations = 0;
	private int maxIterations = 20;
	private int significantDigits = 5;
	private FitFunction fitFunction;
	private int flags;
	private boolean backgroundFitting = true;
	private boolean notSignalFitting = false;
	private double coordinateShift = 1;
	private double shiftFactor = 1;
	private int fitRegionWidth = 0, fitRegionHeight = 0;
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
	private double lambda = 10;
	private boolean computeResiduals = true;
	private double duplicateDistance = 0.5;
	private double bias = 0;
	private double readNoise = 0;
	private double amplification = 0;
	private int maxFunctionEvaluations;
	private SearchMethod searchMethod;
	private boolean gradientLineMinimisation = false;
	private double relativeThreshold = 1e-6;
	private double absoluteThreshold = 1e-16;

	// Options for clamping
	private boolean useClamping = false;
	private boolean useDynamicClamping = false;
	private double[] clampValues;
	private int nClampPeaks;
	private ParameterBounds bounds = null;

	private static double[] defaultClampValues;
	static
	{
		defaultClampValues = new double[7];
		// Taken from the 3D-DAO-STORM paper:
		// (Babcock et al. 2012) A high-density 3D localization algorithm for stochastic optical 
		// reconstruction microscopy. Optical Nanoscopy. 2012 1:6
		// DOI: 10.1186/2192-2853-1-6
		// Page 3
		// Note: It is not clear if the background/signal are in ADUs or photons. I assume photons.
		defaultClampValues[Gaussian2DFunction.BACKGROUND] = 100;
		defaultClampValues[Gaussian2DFunction.SIGNAL] = 1000;
		defaultClampValues[Gaussian2DFunction.SHAPE] = Math.PI;
		defaultClampValues[Gaussian2DFunction.X_POSITION] = 1;
		defaultClampValues[Gaussian2DFunction.Y_POSITION] = 1;
		defaultClampValues[Gaussian2DFunction.X_SD] = 3;
		defaultClampValues[Gaussian2DFunction.Y_SD] = 3;
	}

	private StoppingCriteria stoppingCriteria = null;
	private Gaussian2DFunction gaussianFunction = null;
	private NoiseModel noiseModel = null;
	private BaseFunctionSolver functionSolver = null;

	private DynamicPeakResult dynamicPeakResult = new DynamicPeakResult();

	// Flag to indicate simple filtering is enabled
	//private boolean simpleFilter = true;

	// Support using a smart filter and disabling the simple filtering
	private boolean disableSimpleFilter = false;
	private boolean smartFilter = false;
	private String smartFilterXML = "";
	private DirectFilter directFilter = null;
	private int filterResult = 0;
	private boolean widthEnabled;
	private float offset;

	/**
	 * Default constructor
	 */
	public FitConfiguration()
	{
		initialiseState();
	}

	/**
	 * Creates the appropriate stopping criteria and Gaussian function for the configuration
	 * 
	 * @param npeaks
	 *            The number of peaks to fit
	 * @param maxx
	 *            The width of the XY data
	 * @param maxx
	 *            The height of the XY data
	 * @param params
	 *            The Gaussian parameters
	 */
	public void initialise(int npeaks, int maxx, int maxy, double[] params)
	{
		// Check if the Gaussian function is invalid
		if (gaussianFunction != null && (gaussianFunction.getNPeaks() != npeaks || gaussianFunction.getMaxX() != maxx ||
				gaussianFunction.getMaxY() != maxy))
		{
			// The gaussian function cannot be reused. 
			// Do not call invalidate as it also invalidates the solver and we can re-use that.
			//invalidateGaussianFunction();
			gaussianFunction = null;
		}
		if (gaussianFunction == null)
		{
			// TODO : See if this works ...
			// We can update the function solver with the new function so do not invalidate the solver
			//invalidateFunctionSolver();
			gaussianFunction = createGaussianFunction(npeaks, maxx, maxy, params);
		}
		if (stoppingCriteria == null)
		{
			// Requires a new function solver
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
	public StoppingCriteria createStoppingCriteria(Gaussian2DFunction func, double[] params)
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
	protected void addPeakRestrictions(Gaussian2DFunction func, GaussianStoppingCriteria sc, double[] params)
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
	 * @param maxx
	 *            The height of the XY data
	 * @param params
	 *            The Gaussian parameters
	 * @return The function
	 */
	public Gaussian2DFunction createGaussianFunction(int npeaks, int maxx, int maxy, double[] params)
	{
		final int flags = getFunctionFlags();

		Gaussian2DFunction f = GaussianFunctionFactory.create2D(npeaks, maxx, maxy, flags, null);
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

		if (!isBackgroundFitting())
		{
			// Remove background fitting (on by default)
			flags &= ~GaussianFunctionFactory.FIT_BACKGROUND;
		}
		if (isNotSignalFitting())
		{
			// Remove signal fitting (on by default)
			flags &= ~GaussianFunctionFactory.FIT_SIGNAL;
			flags |= GaussianFunctionFactory.FIT_SIMPLE;
		}
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
				flags = GaussianFunctionFactory.FIT_CIRCLE;
				break;
			case FIXED:
				flags = GaussianFunctionFactory.FIT_FIXED;
				break;
			case FREE_CIRCULAR:
				flags = GaussianFunctionFactory.FIT_FREE_CIRCLE;
				break;
			case FREE:
				flags = GaussianFunctionFactory.FIT_ELLIPTICAL;
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

	//	/**
	//	 * @param fitValidation
	//	 *            True if fit should be validated with {@link #validatePeak(int, double[], double[])}
	//	 */
	//	public void setFitValidation(boolean fitValidation)
	//	{
	//	}

	/**
	 * @return True if fit should be validated with {@link #validatePeak(int, double[], double[])}
	 */
	public boolean isFitValidation()
	{
		return isDirectFilter() || isRegionValidation() || !isDisableSimpleFilter();
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
	 * Set the size of the fit region. Any coordinate outside the region will fail fit validation (see
	 * {@link #validatePeak(int, double[], double[])}). Set to zero to disable.
	 * <p>
	 * Note: it is assumed that the coordinates of the peak are relative to the fit region of size NxN. Coordinates are
	 * offset by the amount defined by {@link #setCoordinateOffset(double)}.
	 *
	 * @param fitRegionWidth
	 *            the fit region width
	 * @param fitRegionHeight
	 *            the fit region height
	 * @param coordinateOffset
	 *            the coordinate offset when validating the coordinates are within the fit window
	 */
	public void setFitRegion(int fitRegionWidth, int fitRegionHeight, double coordinateOffset)
	{
		this.fitRegionWidth = Math.max(0, fitRegionWidth);
		this.fitRegionHeight = Math.max(0, fitRegionHeight);
		this.coordinateOffset = coordinateOffset;
	}

	/**
	 * @return the width of the fit region used for validation
	 */
	public int getFitRegionWidth()
	{
		return fitRegionWidth;
	}

	private boolean isRegionValidation()
	{
		return fitRegionWidth != 0;
	}

	/**
	 * @return the height of the fit region used for validation
	 */
	public int getFitRegionHeight()
	{
		return fitRegionHeight;
	}

	/**
	 * @return the coordinate offset when validating the coordinates are within the fit window
	 */
	public double getCoordinateOffset()
	{
		return coordinateOffset;
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
	public Gaussian2DFunction getGaussianFunction()
	{
		return gaussianFunction;
	}

	/**
	 * Call this when a property changes that will change the Gaussian function
	 */
	private void invalidateGaussianFunction()
	{
		gaussianFunction = null;
		invalidateFunctionSolver();
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
		if (isDirectFilter())
		{
			// Always specify a new result and we have no local background or offset
			PreprocessedPeakResult peak = createPreprocessedPeakResult(0, n, initialParams, params, 0, ResultType.NEW,
					0, 0, false);
			if (directFilter.accept(peak))
				return setValidationResult(FitStatus.OK, null);
			if (log != null)
			{
				log.info("Bad peak %d: %s", peak.getId(),
						DirectFilter.getStatusMessage(peak, directFilter.getResult()));
			}
			if (DirectFilter.anySet(directFilter.getResult(), V_X_SD_FACTOR | V_Y_SD_FACTOR))
			{
				return setValidationResult(FitStatus.WIDTH_DIVERGED, null);
			}
			// At the moment we do not get any other validation data
			return setValidationResult(FitStatus.FAILED_SMART_FILTER, null);
		}

		// Check if outside the fit window.
		// TODO - Make this configurable per peak. At the moment we only use this in BenchmarkSpotFit where 
		// additional peaks will be neighbours. In the future we may want to control this better.
		if (isRegionValidation())
		{
			final int offset = n * 6;
			final double x = params[Gaussian2DFunction.X_POSITION + offset] + coordinateOffset;
			final double y = params[Gaussian2DFunction.Y_POSITION + offset] + coordinateOffset;
			if (x <= 0 || x >= fitRegionWidth || y <= 0 || y >= fitRegionHeight)
			{
				if (log != null)
				{
					log.info("Bad peak %d: Coordinates outside fit region (x=%g,y=%g) <> %d,%d", n, x, y,
							fitRegionWidth, fitRegionHeight);
				}
				return setValidationResult(FitStatus.OUTSIDE_FIT_REGION,
						new double[] { x, y, fitRegionWidth, fitRegionHeight });
			}
		}

		if (isDisableSimpleFilter())
			return setValidationResult(FitStatus.OK, null);

		final int offset = n * 6;
		// Check spot movement
		final double xShift = params[Gaussian2DFunction.X_POSITION + offset] -
				initialParams[Gaussian2DFunction.X_POSITION + offset];
		final double yShift = params[Gaussian2DFunction.Y_POSITION + offset] -
				initialParams[Gaussian2DFunction.Y_POSITION + offset];
		final double maxShift = coordinateShift;
		if (Math.abs(xShift) > maxShift || Math.abs(yShift) > maxShift)
		{
			if (log != null)
			{
				log.info("Bad peak %d: Fitted coordinates moved (x=%g,y=%g) > %g", n, xShift, yShift, maxShift);
			}
			return setValidationResult(FitStatus.COORDINATES_MOVED, new double[] { xShift, yShift });
		}

		// Check signal threshold
		final double signal = params[Gaussian2DFunction.SIGNAL + offset] * gain;
		// Compare the signal to the desired signal strength
		if (signal < signalThreshold)
		{
			if (log != null)
			{
				log.info("Bad peak %d: Insufficient signal %g (SNR=%g)\n", n, signal, signal / noise);
			}
			//if (params.length == 7) // Single peak
			//	System.out.printf("Bad peak %d: Insufficient signal (%gx)\n", n, signal / noise);
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
			final double variance = getVariance(params[Gaussian2DFunction.BACKGROUND],
					params[Gaussian2DFunction.SIGNAL + offset], sd, this.precisionUsingBackground);

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
	 * solver, precision using background, gain, nm per pixel).
	 * <p>
	 * We can calculate the precision using the estimated noise for the image or using the expected number
	 * of background photons at the location.
	 *
	 * @param localBackground
	 *            The background (in photons)
	 * @param signal
	 *            The signal (in photons)
	 * @param sd
	 *            The spot standard deviation
	 * @param precisionUsingBackground
	 *            calculate the precision using expected number of background photons at the location
	 * @return The localisation variance
	 */
	public double getVariance(double localBackground, final double signal, final double sd,
			boolean precisionUsingBackground)
	{
		double variance = 0;
		if (precisionUsingBackground)
		{
			// Check using the formula which uses the estimated background.
			// This allows for better filtering when the background is variable, e.g. when imaging cells.
			if (fitSolver == FitSolver.MLE)
			{
				try
				{
					// This may be slow due to the integration required within the formula.
					variance = PeakResult.getMLVarianceX(nmPerPixel, nmPerPixel * sd, signal,
							Math.max(0, localBackground), emCCD);
				}
				catch (Exception e)
				{
					// Catch all exceptions. They are likely to be a TooManyIterationsException and other
					// problems with the integration
					variance = PeakResult.getVarianceX(nmPerPixel, nmPerPixel * sd, signal,
							Math.max(0, localBackground), emCCD);
				}
			}
			else
			{
				variance = PeakResult.getVarianceX(nmPerPixel, nmPerPixel * sd, signal, Math.max(0, localBackground),
						emCCD);
			}
		}
		else
		{
			if (fitSolver == FitSolver.MLE)
			{
				try
				{
					// This may be slow due to the integration required within the formula.
					variance = PeakResult.getMLVariance(nmPerPixel, nmPerPixel * sd, signal, noise, emCCD);
				}
				catch (Exception e)
				{
					// Catch all exceptions. They are likely to be a TooManyIterationsException and other
					// problems with the integration
					variance = PeakResult.getVariance(nmPerPixel, nmPerPixel * sd, signal, noise, emCCD);
				}
			}
			else
			{
				variance = PeakResult.getVariance(nmPerPixel, nmPerPixel * sd, signal, noise, emCCD);
			}
		}
		return variance;
	}

	/**
	 * An object that can return the results in a formatted state for the multi-path filter
	 */
	private class DynamicPeakResult implements PreprocessedPeakResult
	{
		int id, candidateId;
		int offset;
		double[] initialParams;
		double[] params;
		double localBackground;
		boolean existingResult;
		boolean newResult;
		float offsetx;
		float offsety;
		double var, var2;

		DynamicPeakResult(int candidateId, int n, double[] initialParams, double[] params, double localBackground,
				ResultType resultType, float offsetx, float offsety)
		{
			setParameters(candidateId, n, initialParams, params, localBackground, resultType, offsetx, offsety);
		}

		DynamicPeakResult()
		{
			var = var2 = -1;
		}

		/**
		 * Sets the parameters.
		 *
		 * @param candidateId
		 *            the candidate id
		 * @param n
		 *            The peak number
		 * @param initialParams
		 *            The initial peak parameters
		 * @param params
		 *            The fitted peak parameters
		 * @param localBackground
		 *            the local background
		 * @param resultType
		 *            the result type
		 * @param offsetx
		 *            the offsetx
		 * @param offsety
		 *            the offsety
		 */
		void setParameters(int candidateId, int n, double[] initialParams, double[] params, double localBackground,
				ResultType resultType, float offsetx, float offsety)
		{
			this.id = n;
			this.candidateId = candidateId;
			offset = n * 6;
			this.initialParams = initialParams;
			this.params = params;
			this.localBackground = localBackground;
			this.existingResult = resultType == ResultType.EXISTING;
			this.newResult = resultType == ResultType.NEW;
			this.offsetx = offsetx;
			this.offsety = offsety;
			var = var2 = -1;
		}

		public int getFrame()
		{
			// Not implemented
			return 0;
		}

		public int getUniqueId()
		{
			// In the future we may want to use this so throw an exception so we notice
			throw new NotImplementedException("Unique Id not available");
		}

		public int getId()
		{
			return id;
		}

		public int getCandidateId()
		{
			return candidateId;
		}

		public float getSignal()
		{
			return (float) (params[Gaussian2DFunction.SIGNAL + offset] * FitConfiguration.this.gain);
		}

		public float getPhotons()
		{
			return (float) (params[Gaussian2DFunction.SIGNAL + offset]);
		}

		public float getSNR()
		{
			return (float) (params[Gaussian2DFunction.SIGNAL + offset] / getNoise());
		}

		public float getNoise()
		{
			// Comment this out to use the configured local background.
			// If uncommented then the background will be either the local background or the fitted background.
			final double localBackground = getLocalBackground();

			return (float) ((localBackground > 0) ? PeakResult.localBackgroundToNoise(localBackground, gain, emCCD)
					: FitConfiguration.this.noise);
		}

		private double getLocalBackground()
		{
			return (localBackground > 0) ? localBackground : params[Gaussian2DFunction.BACKGROUND];
		}

		public double getLocationVariance()
		{
			// We do not use the local background so set as zero
			if (var == -1)
				var = FitConfiguration.this.getVariance(0, params[Gaussian2DFunction.SIGNAL + offset], getSD(), false);
			return var;
		}

		public double getLocationVariance2()
		{
			if (var2 == -1)
				var2 = FitConfiguration.this.getVariance(getLocalBackground(),
						params[Gaussian2DFunction.SIGNAL + offset], getSD(), true);
			return var2;
		}

		public float getSD()
		{
			return (float) PeakResult.getSD(params[Gaussian2DFunction.X_SD + offset],
					params[Gaussian2DFunction.Y_SD + offset]);
		}

		public float getBackground()
		{
			return (float) params[Gaussian2DFunction.BACKGROUND];
		}

		public float getAmplitude()
		{
			return (float) (params[Gaussian2DFunction.SIGNAL + offset] / (2 * Math.PI *
					params[Gaussian2DFunction.X_SD + offset] * params[Gaussian2DFunction.Y_SD + offset]));
		}

		public float getAngle()
		{
			return (float) params[Gaussian2DFunction.SHAPE + offset];
		}

		public float getX()
		{
			return (float) params[Gaussian2DFunction.X_POSITION + offset] + offsetx;
		}

		public float getY()
		{
			return (float) params[Gaussian2DFunction.Y_POSITION + offset] + offsety;
		}

		public float getXRelativeShift2()
		{
			final double d = (params[Gaussian2DFunction.X_POSITION + offset] -
					initialParams[Gaussian2DFunction.X_POSITION + offset]) /
					initialParams[Gaussian2DFunction.X_SD + offset];
			return (float) (d * d);
		}

		public float getYRelativeShift2()
		{
			final double d = (params[Gaussian2DFunction.Y_POSITION + offset] -
					initialParams[Gaussian2DFunction.Y_POSITION + offset]) /
					initialParams[Gaussian2DFunction.Y_SD + offset];
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
			return existingResult;
		}

		public boolean isNewResult()
		{
			return newResult;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see gdsc.smlm.results.filter.PreprocessedPeakResult#getAssignments(int)
		 */
		public FractionalAssignment[] getAssignments(int predictedId)
		{
			return null;
		}

		public double[] toGaussian2DParameters()
		{
			final double[] p = new double[7];
			p[Gaussian2DFunction.BACKGROUND] = params[Gaussian2DFunction.BACKGROUND];
			System.arraycopy(params, 1 + offset, p, 1, 6);
			p[Gaussian2DFunction.X_POSITION] += offsetx;
			p[Gaussian2DFunction.Y_POSITION] += offsety;
			return p;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see gdsc.smlm.results.filter.PreprocessedPeakResult#setValidationResult(int)
		 */
		public void setValidationResult(int result)
		{
			throw new NotImplementedException("The validation result should not be set on a dynamic result");
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see gdsc.smlm.results.filter.PreprocessedPeakResult#getValidationResult()
		 */
		public int getValidationResult()
		{
			throw new NotImplementedException("The validation result should not be set on a dynamic result");
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see gdsc.smlm.results.filter.PreprocessedPeakResult#ignore()
		 */
		public boolean ignore()
		{
			return false;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see gdsc.smlm.results.filter.PreprocessedPeakResult#isNotDuplicate()
		 */
		public boolean isNotDuplicate()
		{
			return false;
		}
	}

	/**
	 * Create a dynamic object that can return the results in a formatted state for the multi-path filter.
	 * <p>
	 * The result is dynamic in that it computes the values just-in-time using the input array data.
	 * <p>
	 * The result can be a recycled object that is associated with this fit configuration, or a new object. If using the
	 * recycled object then a second call to this method will replace the array data on all references to the object. If
	 * using a new object then this method can be called again with new data and the old reference is still valid.
	 * <p>
	 * Note: All returned objects will be linked with this fit configuration. Thus changing properties such as the gain,
	 * noise or settings for computing the variance will result in changes to the values returned by the
	 * PreprocessedPeakResult.
	 * <p>
	 * Note: XY position may be wrong if the input parameters have not been updated with an offset from fitting a
	 * sub-region.
	 *
	 * @param candidateId
	 *            the candidate id
	 * @param n
	 *            The peak number
	 * @param initialParams
	 *            The initial peak parameters
	 * @param params
	 *            The fitted peak parameters
	 * @param localBackground
	 *            the local background
	 * @param resultType
	 *            the result type
	 * @param offsetx
	 *            the offsetx to adjust the x-position
	 * @param offsety
	 *            the offsety to adjust the y-position
	 * @return A preprocessed peak result
	 */
	public PreprocessedPeakResult createDynamicPreprocessedPeakResult(int candidateId, int n, double[] initialParams,
			double[] params, double localBackground, ResultType resultType, float offsetx, float offsety)
	{
		return createPreprocessedPeakResult(candidateId, n, initialParams, params, localBackground, resultType, offsetx,
				offsety, true);
	}

	/**
	 * Create a dynamic object that can return the results in a formatted state for the multi-path filter.
	 * <p>
	 * The result is dynamic in that it computes the values just-in-time using the input array data.
	 * <p>
	 * The result can be a recycled object that is associated with this fit configuration, or a new object. If using the
	 * recycled object then a second call to this method will replace the array data on all references to the object. If
	 * using a new object then this method can be called again with new data and the old reference is still valid.
	 * <p>
	 * Note: All returned objects will be linked with this fit configuration. Thus changing properties such as the gain,
	 * noise or settings for computing the variance will result in changes to the values returned by the
	 * PreprocessedPeakResult.
	 * <p>
	 * Note: XY position may be wrong if the input parameters have not been updated with an offset from fitting a
	 * sub-region.
	 *
	 * @param candidateId
	 *            the candidate id
	 * @param n
	 *            The peak number
	 * @param initialParams
	 *            The initial peak parameters
	 * @param params
	 *            The fitted peak parameters
	 * @param localBackground
	 *            the local background
	 * @param resultType
	 *            the result type
	 * @param offsetx
	 *            the offsetx to adjust the x-position
	 * @param offsety
	 *            the offsety to adjust the y-position
	 * @param newObject
	 *            Set to true to create a new object, the default uses the object associated with this fit configuration
	 * @return A preprocessed peak result
	 */
	private PreprocessedPeakResult createPreprocessedPeakResult(int candidateId, int n, double[] initialParams,
			double[] params, double localBackground, ResultType resultType, float offsetx, float offsety,
			boolean newObject)
	{
		if (newObject)
			return new DynamicPeakResult(candidateId, n, initialParams, params, localBackground, resultType, offsetx,
					offsety);

		dynamicPeakResult.setParameters(candidateId, n, initialParams, params, localBackground, resultType, offsetx,
				offsety);
		return dynamicPeakResult;
	}

	/**
	 * Create an object that can return the results in a formatted state for the multi-path filter.
	 * <p>
	 * The result is fixed in that it computes the values on construction using the input array data.
	 * <p>
	 * If a local background is provided then it is used instead of the fitted background. The local background can be
	 * computed if a multi-peak fit has been performed since the background will be the global background, The local
	 * background for a peak will be the global background plus the contribution of all the other peaks in the local
	 * region around the peak of interest.
	 * <p>
	 * The local background will be used to estimate the noise in the local region (as photon shot noise) if it is above
	 * the bias.
	 *
	 * @param frame
	 *            the frame
	 * @param candidateId
	 *            the candidate id
	 * @param n
	 *            The peak number
	 * @param initialParameters
	 *            the initial parameters
	 * @param parameters
	 *            the parameters
	 * @param localBackground
	 *            the local background (set to negative to use the fitted background instead)
	 * @param resultType
	 *            the result type
	 * @param offsetx
	 *            the offsetx to adjust the x-position
	 * @param offsety
	 *            the offsety to adjust the y-position
	 * @return A preprocessed peak result
	 */
	public BasePreprocessedPeakResult createPreprocessedPeakResult(int frame, int candidateId, int n,
			double[] initialParameters, double[] parameters, double localBackground, ResultType resultType,
			float offsetx, float offsety)
	{
		final int offset = n * 6;
		final double signal = parameters[offset + Gaussian2DFunction.SIGNAL] * gain;
		final double photons = parameters[offset + Gaussian2DFunction.SIGNAL];
		final double b = (localBackground > 0) ? localBackground : parameters[Gaussian2DFunction.BACKGROUND];
		final double angle = parameters[offset + Gaussian2DFunction.SHAPE];
		final double x = parameters[offset + Gaussian2DFunction.X_POSITION] + offsetx;
		final double y = parameters[offset + Gaussian2DFunction.Y_POSITION] + offsety;
		final double x0 = initialParameters[offset + Gaussian2DFunction.X_POSITION] + offsetx;
		final double y0 = initialParameters[offset + Gaussian2DFunction.Y_POSITION] + offsety;
		final double xsd = parameters[offset + Gaussian2DFunction.X_SD];
		final double ysd = parameters[offset + Gaussian2DFunction.Y_SD];
		final double xsd0 = initialParameters[offset + Gaussian2DFunction.X_SD];
		final double ysd0 = initialParameters[offset + Gaussian2DFunction.Y_SD];
		final double variance = getVariance(0, signal, PeakResult.getSD(xsd, ysd), false);
		final double variance2 = getVariance(b, signal, PeakResult.getSD(xsd, ysd), true);
		// Q. Should noise be the local background or the estimate from the whole image?

		// This uses the local background if specified or the estimate from the whole image 
		//final double noise = (localBackground > bias)
		//		? PeakResult.localBackgroundToNoise(localBackground - bias, this.gain, this.emCCD) : this.noise;

		// This uses the local fitted background to estimate the noise
		final double noise = (b > 0) ? PeakResult.localBackgroundToNoise(b, 1.0, this.emCCD) : this.noise;
		return new BasePreprocessedPeakResult(frame, n, candidateId, signal, photons, noise, b, angle, x, y, x0, y0,
				xsd, ysd, xsd0, ysd0, variance, variance2, resultType);
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
			f.dynamicPeakResult = new DynamicPeakResult();
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
		// Initialise objects that cannot be null
		if (fitCriteria == null)
			fitCriteria = FitCriteria.LEAST_SQUARED_ERROR;
		if (fitSolver == null)
			fitSolver = FitSolver.LVM;
		if (fitFunction == null)
			fitFunction = FitFunction.CIRCULAR;
		if (searchMethod == null)
			searchMethod = SearchMethod.POWELL;

		if (maxFunctionEvaluations == 0)
			maxFunctionEvaluations = 2000;
		if (initialSD0 == 0)
			initialSD0 = 1;
		if (initialSD1 == 0)
			initialSD1 = 1;
		if (shiftFactor == 0)
			setCoordinateShiftFactor(1);
		if (dynamicPeakResult == null)
			dynamicPeakResult = new DynamicPeakResult();

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
	 * @return the gain (or 1 if the gain is invalid)
	 */
	public double getGainSafe()
	{
		return (gain <= 0) ? 1 : gain;
	}

	/**
	 * @param gain
	 *            the gain to use when evaluating a fitted peak's localisation precision. Also used to determine the
	 *            signal threshold (signalThreshold = max(gain x minPhotons, noise x signalStrength)
	 */
	public void setGain(double gain)
	{
		invalidateFunctionSolver();
		this.gain = Math.abs(gain);
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
		invalidateFunctionSolver();
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
		invalidateGaussianFunction();
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
	 * @return Set to true if the bias should be removed from the data before fitting
	 * @deprecated The bias is always removed before fitting
	 */
	@Deprecated
	public boolean isRemoveBiasBeforeFitting()
	{
		return true;
	}

	/**
	 * @return Set to true if the gain should be removed from the data before fitting.
	 * @deprecated The gain should always be removed before fitting as the system should produce estimates in
	 *             photo-electrons. Only the legacy MLE solvers that explicitly model the camera noise require the data
	 *             and estimate to be in ADUs.
	 */
	@Deprecated
	public boolean isApplyGainBeforeFitting()
	{
		// Only the legacy MLE solvers that explicitly model the camera noise require the data
		// and estimate to be in ADUs.
		if (fitSolver == FitSolver.MLE)
			return false;
		return getGain() != 0;
	}

	/**
	 * @return Set to true if fitting requires the camera counts, i.e. amplification is explicitly modelled during
	 *         fitting.
	 */
	public boolean isFitCameraCounts()
	{
		// Only the legacy MLE solvers that explicitly model the camera noise require the data
		// and estimate to be in ADUs.
		if (fitSolver == FitSolver.MLE)
			return true;
		return false;
	}

	/**
	 * The function solver requires a strictly positive function.
	 *
	 * @return true if requires a strictly positive function.
	 */
	public boolean requireStrictlyPositiveFunction()
	{
		// Only the LSE variants can fit negatives. The MLE variants all require a positive function. 
		switch (fitSolver)
		{
			case BOUNDED_LVM:
			case BOUNDED_LVM_WEIGHTED:
			case LVM:
			case LVM_QUASI_NEWTON:
			case LVM_WEIGHTED:
				return false;

			case LVM_MLE:
			case MLE:
				return true;

			default:
				throw new NotImplementedException("Unknown strictly positive requirement: " + fitSolver.getName());
		}
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
		if (functionSolver == null || gaussianFunction == null)
		{
			// The function solver was invalidated so create a new one
			functionSolver = createFunctionSolver();
		}
		else
		{
			// Update with the solver with the latest function.
			functionSolver.setGradientFunction(gaussianFunction);

			// Note: We must carefully update anything that depends on the function.
			if (bounds != null)
			{
				// We have to update the clamping.
				// Note this code is only executed if the clamp settings have not changed
				// (since changes to those settings invalidate the solver) and the 
				// function settings have not changed (since that invalidates the function
				// and the solver).
				// All that is different is the number of peaks in the function. 
				if (gaussianFunction.getNPeaks() > nClampPeaks)
					setClampValues(bounds);
			}
		}
		return functionSolver;
	}

	/**
	 * Call this when a property changes that will change the function solver
	 */
	private void invalidateFunctionSolver()
	{
		functionSolver = null;
		bounds = null;
	}

	private BaseFunctionSolver createFunctionSolver()
	{
		if (gaussianFunction == null)
		{
			// Other code may want to call getFunctionSolver() to see if exceptions are thrown
			// so create a dummy function so we can return a function solver.
			gaussianFunction = createGaussianFunction(1, 1, 1, null);
		}

		// Remove noise model
		gaussianFunction.setNoiseModel(null);

		NonLinearFit nlinfit;

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
				if (amplification <= 0)
				{
					throw new IllegalArgumentException("The amplification is required for the " + fitSolver.getName());
				}

				fitter.setAlpha(1.0 / amplification);

				// TODO - Configure better stopping criteria ...

				return fitter;

			case BOUNDED_LVM_WEIGHTED:
				gaussianFunction.setNoiseModel(getNoiseModel());
			case BOUNDED_LVM:
			case LVM_MLE:

				if (gain <= 0)
				{
					throw new IllegalArgumentException("The gain is required for the " + fitSolver.getName());
				}

				if (useClamping)
				{
					bounds = new ParameterBounds(gaussianFunction);
					setClampValues(bounds);
				}
				BoundedNonLinearFit bnlinfit = new BoundedNonLinearFit(gaussianFunction, getStoppingCriteria(), bounds);
				nlinfit = bnlinfit;

				break;

			case LVM_QUASI_NEWTON:
				if (gain <= 0)
				{
					throw new IllegalArgumentException("The gain is required for the " + fitSolver.getName());
				}
				// This only works with a Gaussian2DFunction
				if (gaussianFunction instanceof Gaussian2DFunction)
				{
					ApacheLVMFitter apacheNLinFit = new ApacheLVMFitter((Gaussian2DFunction) gaussianFunction);
					apacheNLinFit.setMaxEvaluations(maxIterations);
					// TODO - Configure stopping criteria ...
					return apacheNLinFit;
				}
				// else fall through to default LVM fitter

			case LVM_WEIGHTED:
			case LVM:
			default:
				if (gain <= 0)
				{
					throw new IllegalArgumentException("The gain is required for the " + fitSolver.getName());
				}
				// Only set the weighting function if necessary
				if (fitSolver == FitSolver.LVM_WEIGHTED)
					gaussianFunction.setNoiseModel(getNoiseModel());
				nlinfit = new NonLinearFit(gaussianFunction, getStoppingCriteria());
		}

		nlinfit.setInitialLambda(getLambda());
		if (fitSolver == FitSolver.LVM_MLE)
		{
			nlinfit.setMLE(true);
		}
		return nlinfit;
	}

	private void setClampValues(ParameterBounds bounds)
	{
		double[] clamp = getClampValues();
		// Note: The units are photons. This is OK as all solvers except the legacy MLE fit in photons.
		nClampPeaks = gaussianFunction.getNPeaks();
		double[] clampValues = new double[1 + 6 * nClampPeaks];
		clampValues[Gaussian2DFunction.BACKGROUND] = clamp[Gaussian2DFunction.BACKGROUND];
		for (int i = 0; i < nClampPeaks; i++)
		{
			for (int j = 1; j <= 6; j++)
				clampValues[i + j] = clamp[j];
		}
		bounds.setClampValues(clampValues);
		bounds.setDynamicClamp(useDynamicClamping);
	}

	/**
	 * @return True if simple filtering is disabled
	 */
	public boolean isDisableSimpleFilter()
	{
		return disableSimpleFilter;
	}

	/**
	 * @param Set
	 *            to true to diable simple filtering during validation
	 */
	public void setDisableSimpleFilter(boolean disableSimpleFilter)
	{
		this.disableSimpleFilter = disableSimpleFilter;
	}

	/**
	 * @return True if filtering should use the configured smart filter
	 */
	public boolean isSmartFilter()
	{
		return smartFilter;
	}

	/**
	 * @param smartFilter
	 *            True if filtering should use the configured smart filter
	 */
	public void setSmartFilter(boolean smartFilter)
	{
		this.smartFilter = smartFilter;
	}

	/**
	 * Checks if smart filter is enabled and a valid filter is present.
	 *
	 * @return true, if is direct filtering is enabled.
	 */
	public boolean isDirectFilter()
	{
		return smartFilter && directFilter != null;
	}

	/**
	 * @return the smart filter XML
	 */
	public String getSmartFilterXML()
	{
		return smartFilterXML;
	}

	/**
	 * This returns the representation of this object as a smart filter. This ignores any current smart filter and
	 * only uses the standard filtering settings.
	 * 
	 * @return the smart filter if using this object as a smart filter.
	 */
	public DirectFilter getDefaultSmartFilter()
	{
		double signal = getMinPhotons();
		float snr = (float) getSignalStrength();
		double minWidth = getMinWidthFactor();
		double maxWidth = getWidthFactor();
		double shift = getCoordinateShiftFactor();
		double eshift = 0;
		double precision = getPrecisionThreshold();

		DirectFilter f = (isPrecisionUsingBackground())
				? new MultiFilter2(signal, snr, minWidth, maxWidth, shift, eshift, precision)
				: new MultiFilter(signal, snr, minWidth, maxWidth, shift, eshift, precision);
		return f;
	}

	/**
	 * This returns the XML representation of this object as a smart fitler. This ignores any current smart filter and
	 * only uses the standard filtering settings.
	 * 
	 * @return the smart filter XML if using this object as a smart filter.
	 */
	public String getDefaultSmartFilterXML()
	{
		return getDefaultSmartFilter().toXML();
	}

	/**
	 * @param smartFilterXML
	 *            the smart filter XML to set
	 */
	public void setSmartFilterXML(String smartFilterXML)
	{
		this.smartFilterXML = smartFilterXML;
	}

	/**
	 * Sets the direct filter. This changes the smart filter flag to true and updates the smart filter XML.
	 *
	 * @param directFilter
	 *            the new direct filter
	 */
	public void setDirectFilter(DirectFilter directFilter)
	{
		this.directFilter = directFilter;
		this.smartFilter = directFilter != null;
		if (smartFilter)
		{
			smartFilterXML = directFilter.toXML();
		}
		else
		{
			smartFilterXML = "";
		}
	}

	/**
	 * Gets the smart filter name, if a smart filter exists
	 *
	 * @return the smart filter name
	 */
	public String getSmartFilterName()
	{
		if (directFilter != null)
		{
			return directFilter.getName();
		}
		return "";
	}

	/**
	 * Gets the smart filter, if a smart filter exists. A clone is returned.
	 *
	 * @return the smart filter (or null)
	 */
	public DirectFilter getSmartFilter()
	{
		if (directFilter != null)
		{
			return (DirectFilter) directFilter.clone();
		}
		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IDirectFilter#setup()
	 */
	public void setup()
	{
		setup(0);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IDirectFilter#setup(int)
	 */
	public void setup(int flags)
	{
		if (directFilter != null)
		{
			//if (flags == 0)
			//	directFilter.setup();
			//else
			directFilter.setup(flags);
		}
		else
		{
			widthEnabled = !DirectFilter.areSet(flags, DirectFilter.NO_WIDTH);
			offset = (float) ((shiftFactor > 0) ? shiftFactor * shiftFactor : Float.POSITIVE_INFINITY);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IDirectFilter#accept(gdsc.smlm.results.filter.PreprocessedPeakResult)
	 */
	public boolean accept(PreprocessedPeakResult peak)
	{
		return (filterResult = validate(peak)) == 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IDirectFilter#validate(gdsc.smlm.results.filter.PreprocessedPeakResult)
	 */
	public int validate(PreprocessedPeakResult peak)
	{
		final int flags = doValidate(peak);
		if (log == null)
			return flags;
		// Log the error
		if (flags != 0)
			log.info("Bad peak %d (%.1f,%.1f) [%d]: %s", peak.getCandidateId(), peak.getX(), peak.getY(), peak.getId(),
					DirectFilter.getStatusMessage(peak, flags));
		return flags;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IDirectFilter#validate(gdsc.smlm.results.filter.PreprocessedPeakResult)
	 */
	public int doValidate(PreprocessedPeakResult peak)
	{
		if (directFilter != null)
			return directFilter.validate(peak);

		// Do filtering 
		if (peak.getPhotons() < minPhotons)
			return V_PHOTONS;
		if (peak.getSNR() < this.signalStrength)
			return V_SNR;
		if (widthEnabled)
		{
			if (peak.getXSDFactor() > widthFactor || peak.getXSDFactor() < minWidthFactor)
				return V_X_SD_FACTOR;
		}
		if (peak.getXRelativeShift2() > offset)
			return V_X_RELATIVE_SHIFT;
		if (peak.getYRelativeShift2() > offset)
			return V_Y_RELATIVE_SHIFT;
		// Do not support Euclidian shift
		//if (peak.getXRelativeShift2() + peak.getYRelativeShift2() > offset)
		//	return V_X_RELATIVE_SHIFT | V_Y_RELATIVE_SHIFT;
		final double p = (precisionUsingBackground) ? peak.getLocationVariance2() : peak.getLocationVariance();
		if (p > precisionThreshold)
			return V_LOCATION_VARIANCE;
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IDirectFilter#getFilterType()
	 */
	public FilterType getFilterType()
	{
		return FilterType.DIRECT;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IDirectFilter#getResult()
	 */
	public int getResult()
	{
		return filterResult;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IDirectFilter#copy()
	 */
	public IDirectFilter copy()
	{
		return (IDirectFilter) clone();
	}

	/**
	 * @return Set to true to clamp the parameter update to a maximum value
	 */
	public boolean isUseClamping()
	{
		return useClamping;
	}

	/**
	 * Set to true to clamp the parameter update to a maximum value
	 * 
	 * @param useClamping
	 *            Set to true to clamp the parameter update to a maximum value
	 */
	public void setUseClamping(boolean useClamping)
	{
		invalidateFunctionSolver();
		this.useClamping = useClamping;
	}

	/**
	 * @return Set to true to update the clamp values when the parameter update direction changes
	 */
	public boolean isUseDynamicClamping()
	{
		return useDynamicClamping;
	}

	/**
	 * Set to true to update the clamp values when the parameter update direction changes
	 * 
	 * @param useDynamicClamping
	 *            Set to true to update the clamp values when the parameter update direction changes
	 */
	public void setUseDynamicClamping(boolean useDynamicClamping)
	{
		invalidateFunctionSolver();
		this.useDynamicClamping = useDynamicClamping;
	}

	/**
	 * @return the clampValues
	 */
	private double[] getClampValues()
	{
		if (clampValues == null)
			clampValues = defaultClampValues.clone();
		return clampValues;
	}

	/**
	 * @return The clamp value for the background
	 */
	public double getClampBackground()
	{
		return getClampValues()[Gaussian2DFunction.BACKGROUND];
	}

	/**
	 * Sets the clamp value for the background
	 *
	 * @param value
	 *            the new clamp value
	 */
	public void setClampBackground(double value)
	{
		updateClampValue(Gaussian2DFunction.BACKGROUND, value);
	}

	/**
	 * @return The clamp value for the signal
	 */
	public double getClampSignal()
	{
		return getClampValues()[Gaussian2DFunction.SIGNAL];
	}

	/**
	 * Sets the clamp value for the signal
	 *
	 * @param value
	 *            the new clamp value
	 */
	public void setClampSignal(double value)
	{
		updateClampValue(Gaussian2DFunction.SIGNAL, value);
	}

	/**
	 * @return The clamp value for the angle
	 */
	public double getClampAngle()
	{
		return getClampValues()[Gaussian2DFunction.SHAPE];
	}

	/**
	 * Sets the clamp value for the angle
	 *
	 * @param value
	 *            the new clamp value
	 */
	public void setClampAngle(double value)
	{
		updateClampValue(Gaussian2DFunction.SHAPE, value);
	}

	/**
	 * @return The clamp value for the x position
	 */
	public double getClampX()
	{
		return getClampValues()[Gaussian2DFunction.X_POSITION];
	}

	/**
	 * Sets the clamp value for the x position
	 *
	 * @param value
	 *            the new clamp value
	 */
	public void setClampX(double value)
	{
		updateClampValue(Gaussian2DFunction.X_POSITION, value);
	}

	/**
	 * @return The clamp value for the y position
	 */
	public double getClampY()
	{
		return getClampValues()[Gaussian2DFunction.Y_POSITION];
	}

	/**
	 * Sets the clamp value for the y position
	 *
	 * @param value
	 *            the new clamp value
	 */
	public void setClampY(double value)
	{
		updateClampValue(Gaussian2DFunction.Y_POSITION, value);
	}

	/**
	 * @return The clamp value for the x sd
	 */
	public double getClampXSD()
	{
		return getClampValues()[Gaussian2DFunction.X_SD];
	}

	/**
	 * Sets the clamp value for the x sd
	 *
	 * @param value
	 *            the new clamp value
	 */
	public void setClampXSD(double value)
	{
		updateClampValue(Gaussian2DFunction.X_SD, value);
	}

	/**
	 * @return The clamp value for the y sd
	 */
	public double getClampYSD()
	{
		return getClampValues()[Gaussian2DFunction.Y_SD];
	}

	/**
	 * Sets the clamp value for the y sd
	 *
	 * @param value
	 *            the new clamp value
	 */
	public void setClampYSD(double value)
	{
		updateClampValue(Gaussian2DFunction.Y_SD, value);
	}

	private void updateClampValue(int index, double value)
	{
		invalidateFunctionSolver();
		getClampValues()[index] = value;
	}
}