package gdsc.smlm.engine;

import gdsc.core.ij.Utils;
import gdsc.core.logging.Logger;
import gdsc.core.match.FractionalAssignment;
import gdsc.core.utils.Maths;
import gdsc.core.utils.NotImplementedException;
import gdsc.smlm.data.config.CalibrationProtos.Calibration;
import gdsc.smlm.data.config.CalibrationProtos.CameraType;
import gdsc.smlm.data.config.CalibrationProtosHelper;
import gdsc.smlm.data.config.CalibrationWriter;
import gdsc.smlm.data.config.FitProtos.FilterSettings;
import gdsc.smlm.data.config.FitProtos.FitSettings;
import gdsc.smlm.data.config.FitProtos.FitSolver;
import gdsc.smlm.data.config.FitProtos.FitSolverSettings;
import gdsc.smlm.data.config.FitProtos.FitSolverSettingsOrBuilder;
import gdsc.smlm.data.config.FitProtos.SearchMethod;
import gdsc.smlm.data.config.FitProtosHelper;
import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.data.config.PSFProtos.PSF;
import gdsc.smlm.data.config.PSFProtos.PSFParameter;
import gdsc.smlm.data.config.PSFProtos.PSFParameterUnit;
import gdsc.smlm.data.config.PSFProtos.PSFType;
import gdsc.smlm.data.config.PSFProtosHelper;
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.FunctionSolver;
import gdsc.smlm.fitting.Gaussian2DFitConfiguration;
import gdsc.smlm.fitting.nonlinear.BacktrackingFastMLESteppingFunctionSolver;
import gdsc.smlm.fitting.nonlinear.BaseFunctionSolver;
import gdsc.smlm.fitting.nonlinear.FastMLESteppingFunctionSolver;
import gdsc.smlm.fitting.nonlinear.LSELVMSteppingFunctionSolver;
import gdsc.smlm.fitting.nonlinear.LVMSteppingFunctionSolver;
import gdsc.smlm.fitting.nonlinear.MLELVMSteppingFunctionSolver;
import gdsc.smlm.fitting.nonlinear.MaximumLikelihoodFitter;
import gdsc.smlm.fitting.nonlinear.ParameterBounds;
import gdsc.smlm.fitting.nonlinear.SteppingFunctionSolver;
import gdsc.smlm.fitting.nonlinear.ToleranceChecker;
import gdsc.smlm.fitting.nonlinear.WLSELVMSteppingFunctionSolver;
import gdsc.smlm.function.Gradient2Function;
import gdsc.smlm.function.GradientFunction;
import gdsc.smlm.function.PrecomputedFunctionFactory;
import gdsc.smlm.function.gaussian.AstigmatismZModel;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.results.Gaussian2DPeakResultHelper;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResultHelper;
import gdsc.smlm.results.filter.BasePreprocessedPeakResult;
import gdsc.smlm.results.filter.BasePreprocessedPeakResult.ResultType;
import gdsc.smlm.results.filter.DirectFilter;
import gdsc.smlm.results.filter.FilterType;
import gdsc.smlm.results.filter.IDirectFilter;
import gdsc.smlm.results.filter.MultiFilter;
import gdsc.smlm.results.filter.MultiFilter2;
import gdsc.smlm.results.filter.PreprocessedPeakResult;

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
 * Specifies the fitting configuration
 */
public class FitConfiguration implements Cloneable, IDirectFilter, Gaussian2DFitConfiguration
{
	private FitSettings.Builder fitSettings;

	// Extract the settings for convenience
	private CalibrationWriter calibration;
	PSF.Builder psf;
	private FilterSettings.Builder filterSettings;
	private FitSolverSettings.Builder fitSolverSettings;

	private Logger log = null;
	private boolean computeDeviations = false;
	private int flags;
	private AstigmatismZModel astigmatismZModel = null;
	private double coordinateShift = 1;
	private int fitRegionWidth = 0, fitRegionHeight = 0;
	private double coordinateOffset = 0.5;
	private double signalThreshold = 0;
	private double precisionThreshold = 0;
	private boolean isTwoAxisGaussian2D;
	private double nmPerPixel = 0;
	private double gain = 0;
	private double signalToPhotons, signalToCount;
	private boolean emCCD = false;
	private double noise = 0;
	private double minWidthFactor = 0.5;
	private double widthFactor = 2;
	private boolean computeResiduals = true;

	// Options for clamping
	private double[] clampValues;
	private int nClampPeaks;
	private ParameterBounds bounds = null;

	private ToleranceChecker toleranceChecker = null;
	private Gaussian2DFunction gaussianFunction = null;
	private BaseFunctionSolver functionSolver = null;

	private DynamicPeakResult dynamicPeakResult = new DynamicPeakResult();

	// Support using a smart filter and disabling the simple filtering
	private DirectFilter directFilter = null;
	private int filterResult = 0;
	private boolean widthEnabled;
	private float offset;
	private double varianceThreshold;

	private double[] precomputedFunctionValues = null;

	/**
	 * Instantiates a new fit configuration.
	 */
	public FitConfiguration()
	{
		this(FitProtosHelper.defaultFitSettings, CalibrationProtosHelper.defaultCalibration,
				PSFProtosHelper.defaultOneAxisGaussian2DPSF);
	}

	/**
	 * Instantiates a new fit configuration.
	 *
	 * @param fitSettings
	 *            the fit settings
	 * @param calibration
	 *            the calibration
	 * @param psf
	 *            the psf
	 */
	public FitConfiguration(FitSettings fitSettings, Calibration calibration, PSF psf)
	{
		if (fitSettings == null)
			throw new IllegalArgumentException("FitSettings is null");
		if (calibration == null)
			throw new IllegalArgumentException("Calibration is null");
		if (psf == null)
			throw new IllegalArgumentException("PSF is null");
		init(fitSettings.toBuilder(), calibration.toBuilder(), psf.toBuilder());
	}

	/**
	 * Instantiates a new fit configuration.
	 *
	 * @param fitSettings
	 *            the fit settings
	 * @param calibration
	 *            the calibration
	 * @param psf
	 *            the psf
	 */
	public FitConfiguration(FitSettings.Builder fitSettings, Calibration.Builder calibration, PSF.Builder psf)
	{
		if (fitSettings == null)
			throw new IllegalArgumentException("FitSettings is null");
		if (calibration == null)
			throw new IllegalArgumentException("Calibration is null");
		if (psf == null)
			throw new IllegalArgumentException("PSF is null");
		init(fitSettings, calibration, psf);
	}

	/**
	 * Instantiates a new fit configuration. Does not check for null objects.
	 *
	 * @param fitSettings
	 *            the fit settings
	 * @param calibration
	 *            the calibration
	 * @param psf
	 *            the psf
	 * @param dummy
	 *            the dummy parameter
	 */
	FitConfiguration(FitSettings.Builder fitSettings, Calibration.Builder calibration, PSF.Builder psf, boolean dummy)
	{
		init(fitSettings, calibration, psf);
	}

	/**
	 * Initialise the instance.
	 *
	 * @param fitSettings
	 *            the fit settings
	 * @param calibration
	 *            the calibration
	 * @param psf
	 *            the psf
	 */
	private void init(FitSettings.Builder fitSettings, Calibration.Builder calibration, PSF.Builder psf)
	{
		this.fitSettings = fitSettings;

		// Extract for convenience
		this.calibration = new CalibrationWriter(calibration);
		this.psf = psf;
		fitSolverSettings = fitSettings.getFitSolverSettingsBuilder();
		filterSettings = fitSettings.getFilterSettingsBuilder();

		initialiseState();
	}

	/**
	 * Update the fit settings.
	 *
	 * @param fitSettings
	 *            the fit settings
	 */
	void updateFitSettings(FitSettings.Builder fitSettings)
	{
		fitSolverSettings = fitSettings.getFitSolverSettingsBuilder();
		filterSettings = fitSettings.getFilterSettingsBuilder();
		updateFitSolverSettings();
		updateFilterSettings();
	}

	/**
	 * Ensure that the internal state of the object is initialised. This is used after deserialisation since some state
	 * is not saved but restored from other property values.
	 */
	public void initialiseState()
	{
		// Create the state using the settings
		if (dynamicPeakResult == null)
			dynamicPeakResult = new DynamicPeakResult();
		updateCalibration();
		updatePSF();
		updateFitSolverSettings();
		updateFilterSettings();
	}

	/**
	 * Gets the fit settings.
	 *
	 * @return the fit settings
	 */
	public FitSettings getFitSettings()
	{
		return fitSettings.build();
	}

	/**
	 * Merge fit settings.
	 *
	 * @param fitSettings
	 *            the fit settings
	 */
	public void mergeFitSettings(FitSettings fitSettings)
	{
		this.fitSettings.mergeFrom(fitSettings);
		initialiseState();
	}

	/**
	 * Sets the fit settings.
	 *
	 * @param fitSettings
	 *            the new fit settings
	 */
	public void setFitSettings(FitSettings fitSettings)
	{
		this.fitSettings.clear().mergeFrom(fitSettings);
		initialiseState();
	}

	/**
	 * Gets the calibration.
	 *
	 * @return the calibration
	 */
	public Calibration getCalibration()
	{
		return calibration.getCalibration();
	}

	/**
	 * Gets the calibration writer.
	 *
	 * @return the calibration writer
	 */
	public CalibrationWriter getCalibrationWriter()
	{
		return calibration;
	}

	/**
	 * Merge the calibration.
	 *
	 * @param calibration
	 *            the new calibration
	 */
	public void mergeCalibration(Calibration calibration)
	{
		this.calibration.mergeCalibration(calibration);
		updateCalibration();
	}

	/**
	 * Sets the calibration.
	 *
	 * @param calibration
	 *            the new calibration
	 */
	public void setCalibration(Calibration calibration)
	{
		this.calibration.setCalibration(calibration);
		updateCalibration();
	}

	private void updateCalibration()
	{
		// This uses the camera calibration
		invalidateFunctionSolver();

		nmPerPixel = calibration.getNmPerPixel();
		gain = calibration.getCountPerPhoton();
		emCCD = calibration.isEMCCD();

		if (isFitCameraCounts())
		{
			signalToPhotons = 1.0 / gain;
			signalToCount = 1;
		}
		else
		{
			signalToPhotons = 1;
			signalToCount = gain;
		}

		updateSignalThreshold();
		updatePrecisionThreshold();
	}

	/**
	 * Gets the PSF.
	 *
	 * @return the PSF
	 */
	public PSF getPSF()
	{
		return psf.build();
	}

	/**
	 * Merge the PSF.
	 *
	 * @param psf
	 *            the new PSF
	 */
	public void mergePSF(PSF psf)
	{
		this.psf.mergeFrom(psf);
		updatePSF();
	}

	/**
	 * Sets the psf.
	 *
	 * @param psf
	 *            the new psf
	 */
	public void setPSF(PSF psf)
	{
		this.psf.clear().mergeFrom(psf);
		updatePSF();
	}

	private void updatePSF()
	{
		invalidateGaussianFunction();

		int nParams;
		PSFType psfType = psf.getPsfType();
		switch (psfType)
		{
			case ASTIGMATIC_GAUSSIAN_2D:
				flags = GaussianFunctionFactory.FIT_ERF_ASTIGMATISM;
				nParams = 2;
				break;
			case ONE_AXIS_GAUSSIAN_2D:
				if (isFixedPSF())
					flags = GaussianFunctionFactory.FIT_ERF_FIXED;
				else
					flags = GaussianFunctionFactory.FIT_ERF_CIRCLE;
				nParams = 1;
				break;
			case TWO_AXIS_AND_THETA_GAUSSIAN_2D:
				flags = GaussianFunctionFactory.FIT_ELLIPTICAL;
				nParams = 3;
				break;
			case TWO_AXIS_GAUSSIAN_2D:
				flags = GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE;
				nParams = 2;
				break;
			default:
				throw new IllegalArgumentException("FitSettings must be a Gaussian 2D PSF");
		}
		isTwoAxisGaussian2D = PSFHelper.isTwoAxisGaussian2D(psfType);

		while (psf.getParametersCount() > nParams)
			psf.removeParameters(psf.getParametersCount() - 1);

		// Ensure we have enough parameters
		if (psf.getParametersCount() == 0)
		{
			// Create a dummy Sx
			PSFParameter.Builder p = psf.addParametersBuilder(PSFHelper.INDEX_SX);
			p.setName("Sx");
			p.setValue(1);
			p.setUnit(PSFParameterUnit.DISTANCE);
		}
		if (psf.getParametersCount() == 1 && nParams > 1)
		{
			// Duplicate the Sx to Sy
			PSFParameter.Builder p = psf.addParametersBuilder(PSFHelper.INDEX_SY);
			p.setName("Sy");
			p.setValue(psf.getParameters(PSFHelper.INDEX_SX).getValue());
			p.setUnit(PSFParameterUnit.DISTANCE);
		}
		if (psf.getParametersCount() == 2 && nParams > 2)
		{
			// Create a dummy angle
			PSFParameter.Builder p = psf.addParametersBuilder(PSFHelper.INDEX_THETA);
			p.setName("Angle");
			p.setUnit(PSFParameterUnit.ANGLE);
		}

		updateCoordinateShift();
	}

	/**
	 * Gets the FitSolverSettings.
	 *
	 * @return the FitSolverSettings
	 */
	public FitSolverSettings getFitSolverSettings()
	{
		return fitSolverSettings.build();
	}

	/**
	 * Merge the FitSolverSettings.
	 *
	 * @param fitSolverSettings
	 *            the new FitSolverSettings
	 */
	public void mergeFitSolverSettings(FitSolverSettings fitSolverSettings)
	{
		this.fitSolverSettings.mergeFrom(fitSolverSettings);
		updateFitSolverSettings();
	}

	/**
	 * Sets the fit solver settings.
	 *
	 * @param fitSolverSettings
	 *            the new fit solver settings
	 */
	public void setFitSolverSettings(FitSolverSettings fitSolverSettings)
	{
		fitSettings.setFitSolverSettings(fitSolverSettings);
		this.fitSolverSettings = fitSettings.getFitSolverSettingsBuilder();
		//this.fitSolverSettings.clear().mergeFrom(fitSolverSettings);
		updateFitSolverSettings();
	}

	private void updateFitSolverSettings()
	{
		invalidateGaussianFunction();
		invalidateFunctionSolver();
		invalidateToleranceChecker();
		invalidateClampValues();
	}

	/**
	 * Gets the FilterSettings.
	 *
	 * @return the FilterSettings
	 */
	public FilterSettings getFilterSettings()
	{
		return filterSettings.build();
	}

	/**
	 * Merge the FilterSettings.
	 *
	 * @param filterSettings
	 *            the new FilterSettings
	 */
	public void mergeFilterSettings(FilterSettings filterSettings)
	{
		this.filterSettings.clear().mergeFrom(filterSettings);
		updateFilterSettings();
	}

	/**
	 * Sets the filter settings.
	 *
	 * @param filterSettings
	 *            the new filter settings
	 */
	public void setFilterSettings(FilterSettings filterSettings)
	{
		fitSettings.setFilterSettings(filterSettings);
		this.filterSettings = fitSettings.getFilterSettingsBuilder();
		//this.filterSettings.mergeFrom(filterSettings);
		updateFilterSettings();
	}

	private void updateFilterSettings()
	{
		updateSignalThreshold();
		updatePrecisionThreshold();
		updateCoordinateShift();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public FitConfiguration clone()
	{
		return new FitConfiguration(getFitSettings(), getCalibration(), getPSF());

		//		// This is not a complete duplicate. The settings builder objects with the 
		//		// underlying configuration will be the same between all instances. 
		//		try
		//		{
		//			FitConfiguration f = (FitConfiguration) super.clone();
		//			// Reset instance specific objects
		//			f.toleranceChecker = null;
		//			f.gaussianFunction = null;
		//			f.functionSolver = null;
		//			f.setValidationResult(null, null);
		//			f.dynamicPeakResult = new DynamicPeakResult();
		//			return f;
		//		}
		//		catch (CloneNotSupportedException e)
		//		{
		//			// Ignore
		//		}
		//		return null;
	}

	/**
	 * Creates the appropriate stopping criteria and Gaussian function for the configuration.
	 *
	 * @param npeaks
	 *            The number of peaks to fit
	 * @param maxx
	 *            The height of the XY data
	 * @param maxy
	 *            the maxy
	 * @param params
	 *            The Gaussian parameters
	 */
	public void initialise(int npeaks, int maxx, int maxy, double[] params)
	{
		{
			// XXX: For debugging thread safety require new objects for each fit
			//invalidateGaussianFunction();
			//invalidateToleranceChecker();
		}

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
			gaussianFunction = createGaussianFunction(npeaks, maxx, maxy);
		}
		if (toleranceChecker == null)
		{
			// Requires a new function solver
			invalidateFunctionSolver();
			toleranceChecker = createToleranceChecker(gaussianFunction, params);
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
	public ToleranceChecker createToleranceChecker(Gaussian2DFunction func, double[] params)
	{
		FitSolverSettingsOrBuilder s = getFitSolverSettings();
		double relativeValue = s.getRelativeThreshold();
		double absoluteValue = s.getAbsoluteThreshold();
		double relativeParameters = s.getParameterRelativeThreshold();
		double absoluteParameters = s.getParameterAbsoluteThreshold();
		int maxIterations = s.getMaxIterations();
		ToleranceChecker toleranceChecker = new ToleranceChecker(relativeValue, absoluteValue, relativeParameters,
				absoluteParameters, maxIterations);
		return toleranceChecker;
	}

	/**
	 * Creates the appropriate 2D Gaussian function for the configuration.
	 *
	 * @param npeaks
	 *            The number of peaks to fit
	 * @param maxx
	 *            The width of the XY data
	 * @param maxy
	 *            The height of the XY data
	 * @return The function
	 */
	public Gaussian2DFunction createGaussianFunction(int npeaks, int maxx, int maxy)
	{
		return GaussianFunctionFactory.create2D(npeaks, maxx, maxy, getFunctionFlags(), getAstigmatismZModel());
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
	 * Gets the astigmatism Z model.
	 *
	 * @return the astigmatism Z model
	 */
	public AstigmatismZModel getAstigmatismZModel()
	{
		if (astigmatismZModel == null)
		{
			// TODO - support this within the configuration proto object. 
			// This could be added to the PSF proto.
			astigmatismZModel = null;
		}
		return astigmatismZModel;
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
	 * @param initialAngle
	 *            the initialAngle to set
	 */
	public void setInitialAngle(double initialAngle)
	{
		// Provide backward compatibility
		if (psf.getPsfType() != PSFType.TWO_AXIS_AND_THETA_GAUSSIAN_2D)
			throw new IllegalStateException("Not a 2 axis and theta Gaussian 2D PSF");
		psf.getParametersBuilder(PSFHelper.INDEX_THETA).setValue(initialAngle);
	}

	/**
	 * @return the initialAngle
	 */
	public double getInitialAngle()
	{
		// Provide backward compatibility
		if (psf.getPsfType() != PSFType.TWO_AXIS_AND_THETA_GAUSSIAN_2D)
			throw new IllegalStateException("Not a 2 axis and theta Gaussian 2D PSF");
		return psf.getParameters(PSFHelper.INDEX_THETA).getValue();
	}

	/**
	 * Gets the PSF type.
	 *
	 * @return the PSF type
	 */
	public PSFType getPSFType()
	{
		return psf.getPsfType();
	}

	/**
	 * Sets the PSF type.
	 *
	 * @param psfType
	 *            the new PSF type
	 */
	public void setPSFType(PSFType psfType)
	{
		psf.setPsfType(psfType);
		updatePSF();
	}

	/**
	 * @param initialPeakStdDev
	 *            An estimate for the peak standard deviation used to initialise the fit for all dimensions
	 */
	public void setInitialPeakStdDev(double initialPeakStdDev)
	{
		psf.getParametersBuilder(PSFHelper.INDEX_SX).setValue(initialPeakStdDev);
		if (isTwoAxisGaussian2D)
			psf.getParametersBuilder(PSFHelper.INDEX_SY).setValue(initialPeakStdDev);
		updateCoordinateShift();
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
		psf.getParametersBuilder(PSFHelper.INDEX_SX).setValue(initialPeakStdDev0);
		updateCoordinateShift();
	}

	/**
	 * @return An estimate for the combined peak standard deviation
	 */
	public double getInitialPeakStdDev()
	{
		if (isTwoAxisGaussian2D)
			return Gaussian2DPeakResultHelper.getStandardDeviation(getInitialXSD(), getInitialYSD());
		return getInitialXSD();
	}

	/**
	 * @return An estimate for the peak standard deviation used to initialise the fit for dimension 0
	 */
	public double getInitialXSD()
	{
		return psf.getParameters(PSFHelper.INDEX_SX).getValue();
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
		if (!isTwoAxisGaussian2D)
			throw new IllegalStateException("Not a 2 axis Gaussian 2D PSF");
		psf.getParametersBuilder(PSFHelper.INDEX_SY).setValue(initialPeakStdDev1);
		updateCoordinateShift();
	}

	/**
	 * @return An estimate for the peak standard deviation used to initialise the fit for dimension 1
	 */
	public double getInitialYSD()
	{
		if (isTwoAxisGaussian2D)
			return psf.getParameters(PSFHelper.INDEX_SY).getValue();
		return getInitialXSD();
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
		return fitSolverSettings.getFitSolver();
	}

	/**
	 * @param fitSolver
	 *            the fit solver to use to fit the point spread function (PSF)
	 */
	public void setFitSolver(FitSolver fitSolver)
	{
		fitSolverSettings.setFitSolver(fitSolver);
	}

	/**
	 * @param fitSolver
	 *            the fit solver to use to fit the point spread function (PSF)
	 */
	public void setFitSolver(int fitSolver)
	{
		FitSolver f = FitSolver.forNumber(fitSolver);
		if (f != null)
		{
			setFitSolver(f);
		}
	}

	/**
	 * Sets the fixed iterations flag.
	 *
	 * @param fixedIterations
	 *            the fixedIterations to set
	 */
	public void setFixedIterations(boolean fixedIterations)
	{
		invalidateToleranceChecker();
		fitSolverSettings.setFixedIterations(fixedIterations);
	}

	/**
	 * @return the fixedIterations flag
	 */
	public boolean isFixedIterations()
	{
		return fitSolverSettings.getFixedIterations();
	}

	/**
	 * @param maxIterations
	 *            the maxIterations to set
	 */
	public void setMaxIterations(int maxIterations)
	{
		invalidateToleranceChecker();
		fitSolverSettings.setMaxIterations(Math.max(0, maxIterations));
	}

	/**
	 * @return the maxIterations
	 */
	public int getMaxIterations()
	{
		return Math.max(0, fitSolverSettings.getMaxIterations());
	}

	/**
	 * @param backgroundFitting
	 *            True if fitting the background
	 */
	public void setBackgroundFitting(boolean backgroundFitting)
	{
		invalidateGaussianFunction();
		fitSolverSettings.setDisableBackgroundFitting(!backgroundFitting);
	}

	/**
	 * @return True if fitting the background
	 */
	public boolean isBackgroundFitting()
	{
		return !fitSolverSettings.getDisableBackgroundFitting();
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
		fitSolverSettings.setDisableSignalFitting(noSignalFitting);
	}

	/**
	 * @return True if not fitting the signal. The setting only applies to fixed width fitting
	 */
	public boolean isNotSignalFitting()
	{
		return fitSolverSettings.getDisableBackgroundFitting();
	}

	/**
	 * @return True if fitting an elliptical peak (with an angle parameter)
	 */
	public boolean isAngleFitting()
	{
		return (flags & GaussianFunctionFactory.FIT_ANGLE) != 0;
	}

	/**
	 * @return True if fitting the z-position
	 */
	public boolean isZFitting()
	{
		return (flags & GaussianFunctionFactory.FIT_Z) != 0;
	}

	/**
	 * @return True if fitting the peak width in dimension 0
	 */
	public boolean isXSDFitting()
	{
		return (flags & GaussianFunctionFactory.FIT_X_WIDTH) != 0;
	}

	/**
	 * @return True if fitting the peak width in dimension 1
	 */
	public boolean isYSDFitting()
	{
		return (flags & GaussianFunctionFactory.FIT_Y_WIDTH) != 0;
	}

	/**
	 * Set to true to fix the PSF using the initial parameters. This is only supported for a one-axis Gaussian 2D PSF.
	 *
	 * @param fixed
	 *            the new fixed PSF
	 */
	public void setFixedPSF(boolean fixed)
	{
		fitSolverSettings.setFixedPsf(fixed);
		updatePSF();
	}

	/**
	 * Set to true to fix the PSF using the initial parameters. This is only supported for a one-axis Gaussian 2D PSF.
	 *
	 * @return the fixed PSF
	 */
	public boolean isFixedPSF()
	{
		return fitSolverSettings.getFixedPsf();
	}

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
		filterSettings.setShiftFactor(shiftFactor);
		updateCoordinateShift();
	}

	private void updateCoordinateShift()
	{
		double shiftFactor = getCoordinateShiftFactor();
		if (shiftFactor > 0)
		{
			double widthMax = getInitialXSD();
			if (isTwoAxisGaussian2D)
				widthMax = Math.max(widthMax, getInitialYSD());
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
		return filterSettings.getShiftFactor();
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
	 * Set the signal strength used to determine the signal strength for a good fit (signalThreshold =
	 * max(minSignal, noise x signalStrength).
	 * <p>
	 * Note that minSignal is created appropriately from minPhotons using the
	 * type of fitter, see {@link #isFitCameraCounts()}.
	 * 
	 * @param signalStrength
	 *            The signal strength
	 */
	public void setSignalStrength(double signalStrength)
	{
		filterSettings.setSignalStrength(signalStrength);
		updateSignalThreshold();
	}

	/**
	 * @return the signal strength
	 */
	public double getSignalStrength()
	{
		return filterSettings.getSignalStrength();
	}

	/**
	 * @return The minimum number of photons
	 */
	public double getMinPhotons()
	{
		return filterSettings.getMinPhotons();
	}

	/**
	 * Set the minimum photons used to determine the signal strength for a good fit (signalThreshold =
	 * max(minSignal, noise x signalStrength).
	 * <p>
	 * Note that minSignal is created appropriately from minPhotons using the
	 * type of fitter, see {@link #isFitCameraCounts()}.
	 * 
	 * @param minPhotons
	 *            The minimum number of photons
	 */
	public void setMinPhotons(double minPhotons)
	{
		filterSettings.setMinPhotons(minPhotons);
		updateSignalThreshold();
	}

	/**
	 * @return the precision threshold. Used to determine if the peak is a good fit. Requires that the image is
	 *         calibrated
	 */
	public double getPrecisionThreshold()
	{
		return filterSettings.getPrecisionThreshold();
	}

	/**
	 * @param precisionThreshold
	 *            the precisionThreshold to set
	 */
	public void setPrecisionThreshold(double precisionThreshold)
	{
		if (precisionThreshold > 0)
			filterSettings.setPrecisionThreshold(precisionThreshold);
		else
			filterSettings.clearPrecisionThreshold();
		updatePrecisionThreshold();
	}

	private void updatePrecisionThreshold()
	{
		// Store the squared threshold for speed
		if (nmPerPixel > 0 && gain > 0 && calibration.isCCDCamera())
			this.precisionThreshold = Maths.pow2(getPrecisionThreshold());
		else
			this.precisionThreshold = 0;
	}

	/**
	 * @return True if calculating the precision using the fitted background
	 */
	public boolean isPrecisionUsingBackground()
	{
		return filterSettings.getPrecisionUsingBackground();
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
		filterSettings.setPrecisionUsingBackground(precisionUsingBackground);
	}

	/**
	 * Set the image noise used to determine the signal strength for a good fit (signalThreshold =
	 * max(minSignal, noise x signalStrength).
	 * <p>
	 * Note that minSignal is created appropriately from minPhotons using the
	 * type of fitter, see {@link #isFitCameraCounts()}.
	 * 
	 * @param noise
	 *            The image noise.
	 */
	public void setNoise(double noise)
	{
		this.noise = noise;
		updateSignalThreshold();
	}

	/**
	 * @return the image noise
	 */
	public double getNoise()
	{
		return noise;
	}

	private void updateSignalThreshold()
	{
		double minSignal = (isFitCameraCounts()) ? getMinPhotons() * gain : getMinPhotons();
		signalThreshold = Math.max(noise * getSignalStrength(), minSignal);
	}

	/**
	 * @param widthFactor
	 *            The factor difference allowed between widths for a good fit
	 */
	public void setWidthFactor(double widthFactor)
	{
		filterSettings.setMaxWidthFactor(widthFactor);
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
	public double getMaxWidthFactor()
	{
		double widthFactor = filterSettings.getMaxWidthFactor();
		return (widthFactor > 1) ? widthFactor : 0;
	}

	/**
	 * @param minWidthFactor
	 *            The minimum factor difference allowed between widths for a good fit
	 */
	public void setMinWidthFactor(double minWidthFactor)
	{
		filterSettings.setMinWidthFactor(minWidthFactor);
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
		double minWidthFactor = filterSettings.getMinWidthFactor();
		return (minWidthFactor < 1 && minWidthFactor > 0) ? minWidthFactor : 0;
	}

	/**
	 * @param lambda
	 *            the lambda to start the Levenberg-Marquardt fitting process
	 */
	public void setLambda(double lambda)
	{
		invalidateFunctionSolver();
		fitSolverSettings.setLambda(lambda);
	}

	/**
	 * @return the lambda
	 */
	public double getLambda()
	{
		return fitSolverSettings.getLambda();
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
	 * @param nmPerPixel
	 *            the nm per pixel scale to use when evaluating a fitted peak's localisation precision
	 */
	public void setNmPerPixel(double nmPerPixel)
	{
		calibration.setNmPerPixel(nmPerPixel);
		updateCalibration();
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
	 *            the camera gain to use when evaluating a fitted peak's localisation precision.
	 */
	public void setGain(double gain)
	{
		invalidateFunctionSolver();
		calibration.setCountPerPhoton(gain);
		updateCalibration();
		//updateSignalThreshold();
	}

	/**
	 * Specify the camera type used.
	 * <p>
	 * Specifying a CCD camera is relevant when validating results using the localisation precision.
	 *
	 * @param cameraType
	 *            the new camera type
	 */
	public void setCameraType(CameraType cameraType)
	{
		invalidateFunctionSolver();
		calibration.setCameraType(cameraType);
		updateCalibration();
	}

	/**
	 * @return True if modelling the camera noise during maximum likelihood fitting
	 */
	public boolean isModelCamera()
	{
		return fitSolverSettings.getModelCamera();
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
		fitSolverSettings.setModelCamera(modelCamera);
	}

	/**
	 * @param bias
	 *            the camera bias (used for maximum likelihood estimation to evaluate the correct value of the
	 *            observed count)
	 */
	public void setBias(double bias)
	{
		calibration.setBias(bias);
	}

	/**
	 * @param readNoise
	 *            the camera read noise (used for maximum likelihood estimation)
	 */
	public void setReadNoise(double readNoise)
	{
		invalidateFunctionSolver();
		calibration.setReadNoise(readNoise);
	}

	/**
	 * @param amplification
	 *            The amplification [ADUs/electron] (used for maximum likelihood estimation)
	 */
	public void setAmplification(double amplification)
	{
		invalidateFunctionSolver();
		calibration.setCountPerElectron(amplification);
	}

	/**
	 * @return the maximum number of function evaluations for the Maximum Likelihood Estimator
	 */
	public int getMaxFunctionEvaluations()
	{
		return fitSolverSettings.getMaxFunctionEvaluations();
	}

	/**
	 * @param maxFunctionEvaluations
	 *            the maximum number of function evaluations for the Maximum Likelihood Estimator
	 */
	public void setMaxFunctionEvaluations(int maxFunctionEvaluations)
	{
		invalidateFunctionSolver();
		fitSolverSettings.setMaxFunctionEvaluations(maxFunctionEvaluations);
	}

	/**
	 * @return the search for the Maximum Likelihood Estimator
	 */
	public SearchMethod getSearchMethod()
	{
		return fitSolverSettings.getSearchMethod();
	}

	/**
	 * @param searchMethod
	 *            the search for the Maximum Likelihood Estimator
	 */
	public void setSearchMethod(int searchMethod)
	{
		SearchMethod sm = SearchMethod.forNumber(searchMethod);
		if (sm != null)
		{
			setSearchMethod(sm);
		}
	}

	/**
	 * @param searchMethod
	 *            the search for the Maximum Likelihood Estimator
	 */
	public void setSearchMethod(SearchMethod searchMethod)
	{
		invalidateFunctionSolver();
		fitSolverSettings.getSearchMethod();
	}

	/**
	 * This setting applies to the conjugate gradient method of the Maximum Likelihood Estimator
	 * 
	 * @return the gradientLineMinimisation True if using the gradient for line minimisation
	 */
	public boolean isGradientLineMinimisation()
	{
		return fitSolverSettings.getGradientLineMinimisation();
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
		fitSolverSettings.setGradientLineMinimisation(gradientLineMinimisation);
	}

	/**
	 * @return the relative threshold for convergence
	 */
	public double getRelativeThreshold()
	{
		return fitSolverSettings.getRelativeThreshold();
	}

	/**
	 * @param relativeThreshold
	 *            the relative threshold for convergence
	 */
	public void setRelativeThreshold(double relativeThreshold)
	{
		invalidateToleranceChecker();
		fitSolverSettings.setRelativeThreshold(relativeThreshold);
	}

	/**
	 * @return the absolute threshold for convergence
	 */
	public double getAbsoluteThreshold()
	{
		return fitSolverSettings.getAbsoluteThreshold();
	}

	/**
	 * @param absoluteThreshold
	 *            the absolute threshold for convergence
	 */
	public void setAbsoluteThreshold(double absoluteThreshold)
	{
		invalidateToleranceChecker();
		fitSolverSettings.setAbsoluteThreshold(absoluteThreshold);
	}

	/**
	 * @retParameterurn the parameter relative threshold for convergence
	 */
	public double getParameterRelativeThreshold()
	{
		return fitSolverSettings.getParameterRelativeThreshold();
	}

	/**
	 * @param relativeThreshold
	 *            the parameter relative threshold for convergence
	 */
	public void setParameterRelativeThreshold(double relativeThreshold)
	{
		invalidateToleranceChecker();
		fitSolverSettings.setParameterRelativeThreshold(relativeThreshold);
	}

	/**
	 * @retParameterurn the parameter absolute threshold for convergence
	 */
	public double getParameterAbsoluteThreshold()
	{
		return fitSolverSettings.getParameterAbsoluteThreshold();
	}

	/**
	 * @param absoluteThreshold
	 *            the parameter absolute threshold for convergence
	 */
	public void setParameterAbsoluteThreshold(double absoluteThreshold)
	{
		invalidateToleranceChecker();
		fitSolverSettings.setParameterAbsoluteThreshold(absoluteThreshold);
	}

	/**
	 * @return the toleranceChecker
	 */
	public ToleranceChecker getToleranceChecker()
	{
		if (toleranceChecker == null)
		{
			// This can be set to negative for fixed iterations
			int maxIterations = getMaxIterations();
			if (isFixedIterations())
				maxIterations = -maxIterations;

			boolean minimiseValue = isMinimiseValue();
			toleranceChecker = new ToleranceChecker(minimiseValue, getRelativeThreshold(), getAbsoluteThreshold(),
					getParameterRelativeThreshold(), getParameterAbsoluteThreshold(), maxIterations);
		}
		return toleranceChecker;
	}

	private boolean isMinimiseValue()
	{
		switch (getFitSolver())
		{
			case BACKTRACKING_FAST_MLE:
			case FAST_MLE:
				// Maximum likelihood
				return false;

			case LVM_LSE:
			case LVM_WLSE:
				// Minimise sum-of-squares
				return true;

			case LVM_MLE:
				// Minimises the log-likelihood ratio 
				return true;

			case MLE:
				// The legacy MLE actual minimises the negative likelihood.
				// This should not matter anyway since the tolerance checker is not used
				// for this fitter.
				return true;

			default:
				throw new IllegalStateException("Unrecognised fit solver: " + getFitSolver());
		}
	}

	/**
	 * Call this when a property changes that will change the stopping criteria
	 */
	private void invalidateToleranceChecker()
	{
		toleranceChecker = null;
		invalidateFunctionSolver();
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
			final int offset = n * Gaussian2DFunction.PARAMETERS_PER_PEAK;
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

		final int offset = n * Gaussian2DFunction.PARAMETERS_PER_PEAK;
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

		// Check signal threshold. 
		// The threshold should be set in the same units as those used during fitting. 
		final double signal = params[Gaussian2DFunction.SIGNAL + offset];
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
		if (isXSDFitting())
		{
			boolean badWidth = false;
			double xFactor = 0, yFactor = 0;

			xFactor = params[Gaussian2DFunction.X_SD + offset] / initialParams[Gaussian2DFunction.X_SD + offset];
			badWidth = (xFactor > widthFactor || xFactor < minWidthFactor);

			// Always do this (even if badWidth=true) since we need the factor for the return value
			if (isYSDFitting())
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
		if (precisionThreshold > 0)
		{
			final double sd = (isTwoAxisGaussian2D)
					? Gaussian2DPeakResultHelper.getStandardDeviation(params[Gaussian2DFunction.X_SD + offset],
							params[Gaussian2DFunction.Y_SD + offset])
					: params[Gaussian2DFunction.X_SD + offset];
			final double variance = getVariance(params[Gaussian2DFunction.BACKGROUND],
					params[Gaussian2DFunction.SIGNAL + offset] * signalToPhotons, sd, isPrecisionUsingBackground());

			if (variance > precisionThreshold)
			{
				if (log != null)
				{
					final double precision = Math.sqrt(variance);
					log.info("Bad peak %d: Insufficient precision (%gx)\n", n, precision);
				}
				return setValidationResult(FitStatus.INSUFFICIENT_PRECISION, variance);
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
			if (isModelCameraMLE())
			{
				try
				{
					// This may be slow due to the integration required within the formula.
					variance = Gaussian2DPeakResultHelper.getMLVarianceX(nmPerPixel, nmPerPixel * sd, signal,
							Math.max(0, localBackground), emCCD);
				}
				catch (Exception e)
				{
					// Catch all exceptions. They are likely to be a TooManyIterationsException and other
					// problems with the integration
					variance = Gaussian2DPeakResultHelper.getVarianceX(nmPerPixel, nmPerPixel * sd, signal,
							Math.max(0, localBackground), emCCD);
				}
			}
			else
			{
				variance = Gaussian2DPeakResultHelper.getVarianceX(nmPerPixel, nmPerPixel * sd, signal,
						Math.max(0, localBackground), emCCD);
			}
		}
		else
		{
			if (isModelCameraMLE())
			{
				try
				{
					// This may be slow due to the integration required within the formula.
					variance = Gaussian2DPeakResultHelper.getMLVariance(nmPerPixel, nmPerPixel * sd, signal, noise,
							emCCD);
				}
				catch (Exception e)
				{
					// Catch all exceptions. They are likely to be a TooManyIterationsException and other
					// problems with the integration
					variance = Gaussian2DPeakResultHelper.getVariance(nmPerPixel, nmPerPixel * sd, signal, noise,
							emCCD);
				}
			}
			else
			{
				variance = Gaussian2DPeakResultHelper.getVariance(nmPerPixel, nmPerPixel * sd, signal, noise, emCCD);
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
			offset = n * Gaussian2DFunction.PARAMETERS_PER_PEAK;
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
			return (float) (params[Gaussian2DFunction.SIGNAL + offset] * signalToCount);
		}

		public float getPhotons()
		{
			return (float) (params[Gaussian2DFunction.SIGNAL + offset] * signalToPhotons);
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

			return (float) ((localBackground > 0)
					? (isFitCameraCounts()) ? PeakResultHelper.localBackgroundToNoise(localBackground, gain, emCCD)
							: PeakResultHelper.localBackgroundToNoise(localBackground, emCCD)
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
				var = FitConfiguration.this.getVariance(0, params[Gaussian2DFunction.SIGNAL + offset] * signalToPhotons,
						getSD(), false);
			return var;
		}

		public double getLocationVariance2()
		{
			if (var2 == -1)
				var2 = FitConfiguration.this.getVariance(getLocalBackground(),
						params[Gaussian2DFunction.SIGNAL + offset] * signalToPhotons, getSD(), true);
			return var2;
		}

		public float getSD()
		{
			if (isTwoAxisGaussian2D)
				return (float) Gaussian2DPeakResultHelper.getStandardDeviation(params[Gaussian2DFunction.X_SD + offset],
						params[Gaussian2DFunction.Y_SD + offset]);
			return (float) params[Gaussian2DFunction.X_SD + offset];
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
			return (float) params[Gaussian2DFunction.ANGLE + offset];
		}

		public float getX()
		{
			return (float) params[Gaussian2DFunction.X_POSITION + offset] + offsetx;
		}

		public float getY()
		{
			return (float) params[Gaussian2DFunction.Y_POSITION + offset] + offsety;
		}

		public float getZ()
		{
			return (float) params[Gaussian2DFunction.Z_POSITION + offset];
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
			final double[] p = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
			p[Gaussian2DFunction.BACKGROUND] = params[Gaussian2DFunction.BACKGROUND];
			System.arraycopy(params, 1 + offset, p, 1, Gaussian2DFunction.PARAMETERS_PER_PEAK);
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
		final int offset = n * Gaussian2DFunction.PARAMETERS_PER_PEAK;
		final double signal = parameters[offset + Gaussian2DFunction.SIGNAL] * signalToCount;
		final double photons = parameters[offset + Gaussian2DFunction.SIGNAL] * signalToPhotons;
		final double b = signalToPhotons *
				((localBackground > 0) ? localBackground : parameters[Gaussian2DFunction.BACKGROUND]);
		final double angle = parameters[offset + Gaussian2DFunction.ANGLE];
		final double x = parameters[offset + Gaussian2DFunction.X_POSITION] + offsetx;
		final double y = parameters[offset + Gaussian2DFunction.Y_POSITION] + offsety;
		final double z = parameters[offset + Gaussian2DFunction.Z_POSITION];
		final double x0 = initialParameters[offset + Gaussian2DFunction.X_POSITION] + offsetx;
		final double y0 = initialParameters[offset + Gaussian2DFunction.Y_POSITION] + offsety;
		final double xsd = parameters[offset + Gaussian2DFunction.X_SD];
		final double ysd = parameters[offset + Gaussian2DFunction.Y_SD];
		final double xsd0 = initialParameters[offset + Gaussian2DFunction.X_SD];
		final double ysd0 = initialParameters[offset + Gaussian2DFunction.Y_SD];
		final double sd = (isTwoAxisGaussian2D) ? Gaussian2DPeakResultHelper.getStandardDeviation(xsd, ysd) : xsd;
		final double variance = getVariance(0, signal, sd, false);
		final double variance2 = getVariance(b, signal, sd, true);
		// Q. Should noise be the local background or the estimate from the whole image?

		// This uses the local background if specified or the estimate from the whole image 
		//final double noise = (localBackground > bias)
		//		? PeakResult.localBackgroundToNoise(localBackground - bias, this.gain, this.emCCD) : this.noise;

		// This uses the local fitted background to estimate the noise
		final double noise = (b > 0) ? PeakResultHelper.localBackgroundToNoise(b, 1.0, this.emCCD) : this.noise;
		return new BasePreprocessedPeakResult(frame, n, candidateId, signal, photons, noise, b, angle, x, y, z, x0, y0,
				xsd, ysd, xsd0, ysd0, variance, variance2, resultType);
	}

	/**
	 * @return Set to true if fitting requires the camera counts, i.e. amplification is explicitly modelled during
	 *         fitting.
	 */
	public boolean isFitCameraCounts()
	{
		// Only the legacy MLE solvers that explicitly model the camera noise require the data
		// and estimate to be in ADUs.
		if (getFitSolver() == FitSolver.MLE)
			return true;
		return false;
	}

	/**
	 * Checks if is a full maximum likelihood estimator (MLE) modelling the CCD camera.
	 *
	 * @return true, if is a MLE modelling the camera
	 */
	public boolean isModelCameraMLE()
	{
		return (isModelCamera() && fitSolverSettings.getFitSolverValue() == FitSolver.MLE.getNumber());
	}

	/**
	 * The function solver requires a strictly positive function.
	 *
	 * @return true if requires a strictly positive function.
	 */
	public boolean requireStrictlyPositiveFunction()
	{
		// Only the LSE variants can fit negatives. The MLE variants all require a positive function. 
		switch (getFitSolver())
		{
			case LVM_LSE:
			case LVM_WLSE:
				return false;

			case LVM_MLE:
			case MLE:
			case FAST_MLE:
			case BACKTRACKING_FAST_MLE:
				return true;

			default:
				throw new NotImplementedException("Unknown strictly positive requirement: " + getFitSolver());
		}
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

		if (precomputedFunctionValues != null)
		{
			functionSolver.setGradientFunction((GradientFunction) PrecomputedFunctionFactory
					.wrapFunction(gaussianFunction, precomputedFunctionValues));
			precomputedFunctionValues = null;
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
			gaussianFunction = createGaussianFunction(1, 1, 1);
		}

		// Only support CCD/EM-CCD at the moment
		if (!calibration.isCCDCamera())
		{
			throw new IllegalArgumentException("CCD/EM-CCD camera is required for fit solver: " + getFitSolver());
		}

		if (getFitSolver() == FitSolver.MLE)
		{
			// This requires the gain
			checkCalibration();

			MaximumLikelihoodFitter.SearchMethod searchMethod = convertSearchMethod();

			// Only the Poisson likelihood function supports gradients
			if (searchMethod.usesGradients() && isModelCamera())
			{
				throw new IllegalArgumentException(String.format(
						"The derivative based search method '%s' can only be used with the " +
								"'%s' likelihood function, i.e. no model camera noise",
						searchMethod, MaximumLikelihoodFitter.LikelihoodFunction.POISSON));
			}

			MaximumLikelihoodFitter fitter = new MaximumLikelihoodFitter(gaussianFunction);
			fitter.setRelativeThreshold(getRelativeThreshold());
			fitter.setAbsoluteThreshold(getAbsoluteThreshold());
			fitter.setMaxEvaluations(getMaxFunctionEvaluations());
			fitter.setMaxIterations(getMaxIterations());
			fitter.setSearchMethod(searchMethod);
			fitter.setGradientLineMinimisation(isGradientLineMinimisation());

			// Specify the likelihood function to use
			if (isModelCamera())
			{
				// Set the camera read noise.
				// Do not check if this is set as 0 is a valid option.
				fitter.setSigma(calibration.getReadNoise());

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
			if (!calibration.hasCountPerElectron())
			{
				throw new IllegalArgumentException(
						"The amplification is required for the fit solver: " + getFitSolver());
			}

			fitter.setAlpha(1.0 / calibration.getCountPerElectron());

			// TODO - Configure better stopping criteria ...

			return fitter;
		}

		// All the remaining solvers are based on the stepping function solver
		ToleranceChecker tc = getToleranceChecker();
		ParameterBounds bounds = new ParameterBounds(gaussianFunction);
		if (isUseClamping())
		{
			setClampValues(bounds);
		}

		SteppingFunctionSolver solver;

		switch (getFitSolver())
		{
			case LVM_LSE:
				solver = new LSELVMSteppingFunctionSolver(gaussianFunction, tc, bounds);
				break;

			case LVM_MLE:
				checkCalibration();
				solver = new MLELVMSteppingFunctionSolver(gaussianFunction, tc, bounds);
				break;

			case LVM_WLSE:
				checkCalibration();
				solver = new WLSELVMSteppingFunctionSolver(gaussianFunction, tc, bounds);
				break;

			case FAST_MLE:
				checkCalibration();
				// This may throw a class cast exception is the function does not support
				// the Gradient2Function interface
				solver = new FastMLESteppingFunctionSolver((Gradient2Function) gaussianFunction, tc, bounds);
				break;

			case BACKTRACKING_FAST_MLE:
				checkCalibration();
				solver = new BacktrackingFastMLESteppingFunctionSolver((Gradient2Function) gaussianFunction, tc,
						bounds);
				break;

			default:
				throw new IllegalArgumentException("Unknown fit solver: " + getFitSolver());
		}

		if (solver instanceof LVMSteppingFunctionSolver)
		{
			((LVMSteppingFunctionSolver) solver).setInitialLambda(getLambda());
		}
		return solver;
	}

	private gdsc.smlm.fitting.nonlinear.MaximumLikelihoodFitter.SearchMethod convertSearchMethod()
	{
		return FitProtosHelper.convertSearchMethod(getSearchMethod());
	}

	private void checkCalibration()
	{
		if (gain <= 0)
		{
			throw new IllegalArgumentException("The gain is required for fit solver: " + getFitSolver());
		}
	}

	private void setClampValues(ParameterBounds bounds)
	{
		double[] clamp = getClampValues();
		// Note: The units are photons. This is OK as all solvers except the legacy MLE fit in photons.
		nClampPeaks = gaussianFunction.getNPeaks();
		double[] clampValues = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK * nClampPeaks];
		clampValues[Gaussian2DFunction.BACKGROUND] = clamp[Gaussian2DFunction.BACKGROUND];
		for (int i = 0; i < nClampPeaks; i++)
		{
			for (int j = 1; j <= Gaussian2DFunction.PARAMETERS_PER_PEAK; j++)
				clampValues[i + j] = clamp[j];
		}
		bounds.setClampValues(clampValues);
		bounds.setDynamicClamp(isUseDynamicClamping());
	}

	/**
	 * @return True if simple filtering is disabled
	 */
	public boolean isDisableSimpleFilter()
	{
		return filterSettings.getDisableSimpleFilter();
	}

	/**
	 * @param Set
	 *            to true to diable simple filtering during validation
	 */
	public void setDisableSimpleFilter(boolean disableSimpleFilter)
	{
		filterSettings.setDisableSimpleFilter(disableSimpleFilter);
	}

	/**
	 * @return True if filtering should use the configured smart filter
	 */
	public boolean isSmartFilter()
	{
		return filterSettings.getSmartFilter();
	}

	/**
	 * @param smartFilter
	 *            True if filtering should use the configured smart filter
	 */
	public void setSmartFilter(boolean smartFilter)
	{
		filterSettings.setSmartFilter(smartFilter);
	}

	/**
	 * Checks if smart filter is enabled and a valid filter is present.
	 *
	 * @return true, if is direct filtering is enabled.
	 */
	public boolean isDirectFilter()
	{
		return isSmartFilter() && directFilter != null;
	}

	/**
	 * @return the smart filter string
	 */
	public String getSmartFilterString()
	{
		String s = filterSettings.getSmartFilterString();
		return (Utils.isNullOrEmpty(s)) ? "" : s;
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
		double maxWidth = getMaxWidthFactor();
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
	 * Sets the direct filter. This changes the smart filter flag to true and updates the smart filter string.
	 *
	 * @param directFilter
	 *            the new direct filter
	 */
	public void setDirectFilter(DirectFilter directFilter)
	{
		this.directFilter = directFilter;
		if (directFilter != null)
		{
			setSmartFilter(true);
			filterSettings.setSmartFilterString(directFilter.toXML());
		}
		else
		{
			setSmartFilter(false);
			filterSettings.clearSmartFilterString();
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
			double shiftFactor = getCoordinateShiftFactor();
			offset = (float) ((shiftFactor > 0) ? shiftFactor * shiftFactor : Float.POSITIVE_INFINITY);
			varianceThreshold = (precisionThreshold > 0) ? precisionThreshold : Double.POSITIVE_INFINITY;
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
		if (peak.getPhotons() < getMinPhotons())
			return V_PHOTONS;
		if (peak.getSNR() < getSignalStrength())
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
		final double p = (isPrecisionUsingBackground()) ? peak.getLocationVariance2() : peak.getLocationVariance();
		if (p > varianceThreshold)
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
		return fitSolverSettings.getUseClamping();
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
		fitSolverSettings.setUseClamping(useClamping);
	}

	/**
	 * @return Set to true to update the clamp values when the parameter update direction changes
	 */
	public boolean isUseDynamicClamping()
	{
		return fitSolverSettings.getUseDynamicClamping();
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
		fitSolverSettings.setUseDynamicClamping(useDynamicClamping);
	}

	/**
	 * @return the clampValues
	 */
	private double[] getClampValues()
	{
		if (clampValues == null)
		{
			int n = PeakResult.STANDARD_PARAMETERS + 3;
			if (fitSolverSettings.getClampValuesCount() != n)
				throw new IllegalStateException("Require clamp values for all the Gaussian 2D parameters");

			clampValues = new double[n];
			for (int i = 0; i < n; i++)
				clampValues[i] = fitSolverSettings.getClampValues(i);
		}
		return clampValues;
	}

	/**
	 * Invalidate clamp values.
	 */
	private void invalidateClampValues()
	{
		clampValues = null;
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
		updateClampValues(Gaussian2DFunction.BACKGROUND, value);
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
		updateClampValues(Gaussian2DFunction.SIGNAL, value);
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
		updateClampValues(Gaussian2DFunction.X_POSITION, value);
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
		updateClampValues(Gaussian2DFunction.Y_POSITION, value);
	}

	/**
	 * @return The clamp value for the z position
	 */
	public double getClampZ()
	{
		return getClampValues()[Gaussian2DFunction.Z_POSITION];
	}

	/**
	 * Sets the clamp value for the z position
	 *
	 * @param value
	 *            the new clamp value
	 */
	public void setClampZ(double value)
	{
		updateClampValues(Gaussian2DFunction.Z_POSITION, value);
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
		updateClampValues(Gaussian2DFunction.X_SD, value);
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
		updateClampValues(Gaussian2DFunction.Y_SD, value);
	}

	/**
	 * @return The clamp value for the angle
	 */
	public double getClampAngle()
	{
		return getClampValues()[Gaussian2DFunction.ANGLE];
	}

	/**
	 * Sets the clamp value for the angle
	 *
	 * @param value
	 *            the new clamp value
	 */
	public void setClampAngle(double value)
	{
		updateClampValues(Gaussian2DFunction.ANGLE, value);
	}

	private void updateClampValues(int index, double value)
	{
		invalidateFunctionSolver();
		getClampValues()[index] = value;
	}

	/**
	 * Gets the precomputed function values.
	 *
	 * @return the precomputed function values
	 */
	public double[] getPrecomputedFunctionValues()
	{
		return precomputedFunctionValues;
	}

	/**
	 * Sets the precomputed function values. This is combined with the configured Gaussian function add passed to the
	 * function
	 * solver returned from {@link #getFunctionSolver()}. The precomputed function values are then reset to null so it
	 * must be
	 * set each time the function solver is accessed.
	 *
	 * @param precomputedFunctionValues
	 *            the new precomputed function values
	 */
	public void setPrecomputedFunctionValues(double[] precomputedFunctionValues)
	{
		this.precomputedFunctionValues = precomputedFunctionValues;
	}
}