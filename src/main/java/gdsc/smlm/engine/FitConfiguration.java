package gdsc.smlm.engine;

import gdsc.core.data.utils.ConversionException;
import gdsc.core.logging.Logger;
import gdsc.core.match.FractionalAssignment;
import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.Maths;
import gdsc.core.utils.NotImplementedException;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.core.utils.TextUtils;
import gdsc.smlm.data.config.CalibrationProtos.Calibration;
import gdsc.smlm.data.config.CalibrationProtos.CameraType;
import gdsc.smlm.data.config.CalibrationProtosHelper;
import gdsc.smlm.data.config.CalibrationReader;
import gdsc.smlm.data.config.CalibrationWriter;
import gdsc.smlm.data.config.ConfigurationException;
import gdsc.smlm.data.config.FitProtos.FilterSettings;
import gdsc.smlm.data.config.FitProtos.FitSettings;
import gdsc.smlm.data.config.FitProtos.FitSolver;
import gdsc.smlm.data.config.FitProtos.FitSolverSettings;
import gdsc.smlm.data.config.FitProtos.LineSearchMethod;
import gdsc.smlm.data.config.FitProtos.PrecisionMethod;
import gdsc.smlm.data.config.FitProtos.SearchMethod;
import gdsc.smlm.data.config.FitProtosHelper;
import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.data.config.PSFProtos.AstigmatismModel;
import gdsc.smlm.data.config.PSFProtos.PSF;
import gdsc.smlm.data.config.PSFProtos.PSFParameter;
import gdsc.smlm.data.config.PSFProtos.PSFParameterUnit;
import gdsc.smlm.data.config.PSFProtos.PSFType;
import gdsc.smlm.data.config.PSFProtosHelper;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.FunctionSolver;
import gdsc.smlm.fitting.Gaussian2DFitConfiguration;
import gdsc.smlm.fitting.MLEScaledFunctionSolver;
import gdsc.smlm.fitting.ScaledFunctionSolver;
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
import gdsc.smlm.function.gaussian.HoltzerAstigmatismZModel;
import gdsc.smlm.model.camera.CameraModel;
import gdsc.smlm.model.camera.FixedPixelCameraModel;
import gdsc.smlm.model.camera.NullCameraModel;
import gdsc.smlm.results.Gaussian2DPeakResultHelper;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResultHelper;
import gdsc.smlm.results.filter.BasePreprocessedPeakResult;
import gdsc.smlm.results.filter.BasePreprocessedPeakResult.ResultType;
import gdsc.smlm.results.filter.DirectFilter;
import gdsc.smlm.results.filter.Filter;
import gdsc.smlm.results.filter.FilterSetupData;
import gdsc.smlm.results.filter.FilterType;
import gdsc.smlm.results.filter.IDirectFilter;
import gdsc.smlm.results.filter.MultiFilter;
import gdsc.smlm.results.filter.MultiFilter2;
import gdsc.smlm.results.filter.MultiFilterCRLB;
import gdsc.smlm.results.filter.PreprocessedPeakResult;
import gdsc.smlm.results.filter.ShiftFilterSetupData;

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
	private double signalToPhotons;
	private boolean emCCD = false;
	private boolean isMLE = false;
	private double noise = 0;
	private double minWidthFactor = 0.5;
	/**
	 * The width factor. This is the squared width if a 2 axis PSF
	 */
	private double widthFactor = 2;
	private boolean computeResiduals = true;
	private boolean zEnabled = false;

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
	private int filterSetupFlags;
	private FilterSetupData[] filterSetupData;

	private double[] precomputedFunctionValues = null, observationWeights = null;
	private CameraModel cameraModel = null;

	private BaseVarianceSelector varianceSelector = new BaseVarianceSelector();

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
		updatePSF(true);
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
	 * Gets a reference to the current calibration writer.
	 *
	 * @return the calibration writer
	 */
	CalibrationWriter getCalibrationWriterReference()
	{
		return calibration;
	}

	/**
	 * Gets a new calibration writer. Any changes to the writer must be saved using
	 * {@link #setCalibration(Calibration)}.
	 *
	 * @return the calibration writer
	 */
	public CalibrationWriter getCalibrationWriter()
	{
		return new CalibrationWriter(calibration.getCalibration());
	}

	/**
	 * Gets a new calibration reader.
	 *
	 * @return the calibration reader
	 */
	public CalibrationReader getCalibrationReader()
	{
		return new CalibrationReader(calibration.getCalibrationOrBuilder());
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
		invalidateCameraModel();
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
		invalidateCameraModel();
	}

	/**
	 * Update calibration used for signal and precision filtering (i.e. nmPerPixel, gain, emCCD (camera type)).
	 */
	private void updateCalibration()
	{
		// This uses the camera calibration
		invalidateFunctionSolver();

		nmPerPixel = calibration.getNmPerPixel();
		gain = calibration.getCountPerPhoton();
		emCCD = calibration.isEMCCD();

		if (isRawFit())
		{
			// No camera calibration so assume raw data is in photons 
			gain = 1;
		}

		if (isFitCameraCounts())
		{
			signalToPhotons = 1.0 / gain;
		}
		else
		{
			signalToPhotons = 1;
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
		updatePSF(true);
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
		updatePSF(true);
	}

	private void updatePSF(boolean resetAstigmatismModel)
	{
		invalidateGaussianFunction();

		// Reset the astigmatism model. It will be dynamically created from the PSF settings.
		if (resetAstigmatismModel)
			astigmatismZModel = null;

		int nParams;
		PSFType psfType = psf.getPsfType();
		switch (psfType)
		{
			case ASTIGMATIC_GAUSSIAN_2D:
				flags = GaussianFunctionFactory.FIT_ERF_ASTIGMATISM;
				//nParams = 2; // Store Sx and Sy
				nParams = 8; // The PSF stores the full astigmatism model
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
				throw new IllegalStateException("FitSettings must be a Gaussian 2D PSF");
		}
		isTwoAxisGaussian2D = PSFHelper.isTwoAxisGaussian2D(psfType);

		boolean changed = psf.getParametersCount() > nParams;
		if (changed)
		{
			while (psf.getParametersCount() > nParams)
				psf.removeParameters(psf.getParametersCount() - 1);

			// Updated names when changing from two-axis to one axis
			if (psf.getParametersCount() == 1)
			{
				psf.getParametersBuilder(PSFHelper.INDEX_SX).setName(
						PSFProtosHelper.defaultOneAxisGaussian2DPSF.getParameters(PSFHelper.INDEX_SX).getName());
			}
		}

		// Ensure we have enough parameters
		if (psf.getParametersCount() == 0 && nParams > 0)
		{
			// Create a dummy Sx
			PSFParameter.Builder p = psf.addParametersBuilder();
			p.setName(PSFProtosHelper.defaultOneAxisGaussian2DPSF.getParameters(PSFHelper.INDEX_SX).getName());
			p.setValue(1);
			p.setUnit(PSFParameterUnit.DISTANCE);
		}
		if (psf.getParametersCount() == 1 && nParams > 1)
		{
			// Rename S to Sx
			psf.getParametersBuilder(PSFHelper.INDEX_SX)
					.setName(PSFProtosHelper.defaultTwoAxisGaussian2DPSF.getParameters(PSFHelper.INDEX_SX).getName());

			// Duplicate the Sx to Sy
			PSFParameter.Builder p = psf.addParametersBuilder();
			p.setName(PSFProtosHelper.defaultTwoAxisGaussian2DPSF.getParameters(PSFHelper.INDEX_SY).getName());
			p.setValue(psf.getParameters(PSFHelper.INDEX_SX).getValue());
			p.setUnit(PSFParameterUnit.DISTANCE);
		}
		if (psf.getParametersCount() == 2 && nParams > 2)
		{
			// Create a dummy angle
			PSFParameter.Builder p = psf.addParametersBuilder();
			p.setName(
					PSFProtosHelper.defaultTwoAxisAndThetaGaussian2DPSF.getParameters(PSFHelper.INDEX_THETA).getName());
			p.setUnit(PSFParameterUnit.ANGLE);
		}

		// These depend on the 2-axis Gaussian flag
		updateWidthThreshold();
		updateMinWidthThreshold();

		// This depends on the width. 
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
		updateWidthThreshold();
		updateMinWidthThreshold();
		updateZFilter();

		// Recreate the smart filter
		if (!filterSettings.getSmartFilter())
			return;
		String xml = filterSettings.getSmartFilterString();
		if (TextUtils.isNullOrEmpty(xml))
			return;
		Filter f = DirectFilter.fromXML(xml);
		if (f == null || !(f instanceof DirectFilter))
		{
			// Throw to ensure the filter is OK
			throw new IllegalStateException("Unrecognised smart filter: " + xml);
			// or 
			//setDirectFilter(null);
		}

		// This updates the SmartFilter flag and the SmartFilterString.
		// Just set the filter directly
		//setDirectFilter((DirectFilter) f);		
		this.directFilter = (DirectFilter) f;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public FitConfiguration clone()
	{
		// Make a new initialised instance. This will have new settings builder objects.
		return new FitConfiguration(getFitSettings(), getCalibration(), getPSF()).copySettings(this);

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
	 * Copy settings from the other configuration. This copies all the instance fields not stored in settings objects
	 * that can be shared between instances. It is used in the {@link #clone()} method after a new instance has been
	 * created with the current settings objects.
	 *
	 * @param other
	 *            the other configuration
	 * @return the fit configuration
	 */
	FitConfiguration copySettings(FitConfiguration other)
	{
		// Set all the properties that are not updated by a change in the settings
		log = other.log;
		computeDeviations = other.computeDeviations;
		astigmatismZModel = other.astigmatismZModel;
		fitRegionWidth = other.fitRegionWidth;
		fitRegionHeight = other.fitRegionHeight;
		coordinateOffset = other.coordinateOffset;
		noise = other.noise;
		computeResiduals = other.computeResiduals;

		// Support cloning the initialised state from IDirectFilter.setup(...)
		directFilter = other.getSmartFilter(); // This is a clone
		widthEnabled = other.widthEnabled;
		offset = other.offset;
		varianceThreshold = other.varianceThreshold;
		filterSetupFlags = other.filterSetupFlags;
		filterSetupData = other.filterSetupData;

		cameraModel = other.cameraModel;
		varianceSelector = other.varianceSelector;

		return this;
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
		}
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
	 * Gets the PSF type value.
	 *
	 * @return the PSF type value
	 */
	public int getPSFTypeValue()
	{
		return psf.getPsfTypeValue();
	}

	/**
	 * Sets the PSF type.
	 * <p>
	 * If the type is astigmatism and the astigmatism model cannot be constructed from the current PSF parameters then
	 * the result filtering state may be incorrect. It is safer to call {@link #setAstigmatismModel(AstigmatismModel)}
	 * which also updates the PSF type to astigmatism.
	 *
	 * @param psfType
	 *            the new PSF type
	 */
	public void setPSFType(PSFType psfType)
	{
		psf.setPsfType(psfType);
		updatePSF(true);
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
		astigmatismZModel = null;
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
		astigmatismZModel = null;
		updateCoordinateShift();
	}

	/**
	 * Gets the initial peak std dev.
	 *
	 * @return An estimate for the combined peak standard deviation
	 * @throws ConfigurationException
	 *             if the PSF type is astigmatism and the model cannot be constructed
	 */
	public double getInitialPeakStdDev() throws ConfigurationException
	{
		if (isTwoAxisGaussian2D)
			return Gaussian2DPeakResultHelper.getStandardDeviation(getInitialXSD(), getInitialYSD());
		return getInitialXSD();
	}

	/**
	 * Gets the initial XSD.
	 *
	 * @return An estimate for the peak standard deviation used to initialise the fit for dimension 0
	 * @throws ConfigurationException
	 *             if the PSF type is astigmatism and the model cannot be constructed
	 */
	public double getInitialXSD() throws ConfigurationException
	{
		if (getAstigmatismZModel() != null)
			return astigmatismZModel.getSx(0);
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
		astigmatismZModel = null;
		updateCoordinateShift();
	}

	/**
	 * Gets the initial YSD.
	 *
	 * @return An estimate for the peak standard deviation used to initialise the fit for dimension 1
	 * @throws ConfigurationException
	 *             if the PSF type is astigmatism and the model cannot be constructed
	 */
	public double getInitialYSD() throws ConfigurationException
	{
		if (getAstigmatismZModel() != null)
			return astigmatismZModel.getSy(0);
		if (isTwoAxisGaussian2D)
			return psf.getParameters(PSFHelper.INDEX_SY).getValue();
		return getInitialXSD();
	}

	/**
	 * Sets to true to compute the deviations.
	 *
	 * @param computeDeviations
	 *            True if computing the parameter deviations
	 */
	public void setComputeDeviations(boolean computeDeviations)
	{
		this.computeDeviations = computeDeviations;
	}

	/**
	 * {@inheritDoc}
	 * <p>
	 * Note: This is also true if validation is active and the precision method requires computation of the deviations
	 * (see {@link #isFilterRequiresDeviations()}).
	 * 
	 * @see gdsc.smlm.fitting.Gaussian2DFitConfiguration#isComputeDeviations()
	 */
	public boolean isComputeDeviations()
	{
		if (computeDeviations)
			return true;
		return isFilterRequiresDeviations();
	}

	/**
	 * Get the value of the compute deviations flag. This may be false but {@link #isComputeDeviations()} can still
	 * return true.
	 *
	 * @return the value of the compute deviations flag
	 */
	public boolean getComputeDeviationsFlag()
	{
		return computeDeviations;
	}

	/**
	 * Checks if the current filter settings require deviations.
	 *
	 * @return true, if filtering requires deviations
	 */
	public boolean isFilterRequiresDeviations()
	{
		if (isDirectFilter() && directFilter.requiresParameterDeviations())
			return true;
		if (precisionThreshold > 0 && getPrecisionMethodValue() == PrecisionMethod.POISSON_CRLB_VALUE)
			return true;
		return false;
	}

	/**
	 * @return the fit solver used to fit the point spread function (PSF)
	 */
	public FitSolver getFitSolver()
	{
		return fitSolverSettings.getFitSolver();
	}

	/**
	 * @return the fit solver used to fit the point spread function (PSF)
	 */
	public int getFitSolverValue()
	{
		return fitSolverSettings.getFitSolverValue();
	}

	/**
	 * @param fitSolver
	 *            the fit solver to use to fit the point spread function (PSF)
	 */
	public void setFitSolver(FitSolver fitSolver)
	{
		invalidateFunctionSolver();
		updatePrecisionThreshold();
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
		updatePSF(true);
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
			// It may throw if a model cannot be created for the astigmatism type.
			double widthMax;
			try
			{
				widthMax = getWidthMax();
			}
			catch (ConfigurationException e)
			{
				if (getPSFTypeValue() != PSFType.ASTIGMATIC_GAUSSIAN_2D_VALUE)
					throw e;
				// This is OK as the full astigmatism model may not have been set yet.
				// Just use a dummy value of 1
				widthMax = 1;
			}
			if (widthMax > 0)
			{
				setCoordinateShift(shiftFactor * widthMax);
				return;
			}
		}
		setCoordinateShift(Double.POSITIVE_INFINITY);
	}

	/**
	 * Gets the maximum of the initial X and Y widths.
	 *
	 * @return the width max
	 * @throws ConfigurationException
	 *             if the PSF type is astigmatism and the model cannot be constructed
	 */
	public double getWidthMax() throws ConfigurationException
	{
		double widthMax = getInitialXSD();
		if (isTwoAxisGaussian2D)
			widthMax = Math.max(widthMax, getInitialYSD());
		return widthMax;
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
		// Note: Store the squared threshold for speed.
		precisionThreshold = 0;
		switch (getPrecisionMethodValue())
		{
			case PrecisionMethod.MORTENSEN_VALUE:
			case PrecisionMethod.MORTENSEN_LOCAL_BACKGROUND_VALUE:
				// XXX - Determine if precision filtering for SCMOS is valid.
				// For now we leave this in but it may have to be changed to have a precision
				// computed during the fit which is stored for validation.
				if (nmPerPixel > 0 && gain > 0 &&
						//calibration.isCCDCamera()
						(calibration.isCCDCamera() || calibration.isSCMOS()))
					this.precisionThreshold = Maths.pow2(getPrecisionThreshold());
				break;
			case PrecisionMethod.POISSON_CRLB_VALUE:
				this.precisionThreshold = Maths.pow2(getPrecisionThreshold());
				break;
			default:
				break;
		}
		isMLE = isMLE();
	}

	/**
	 * @return True if calculating the precision using the fitted background
	 */
	public boolean isPrecisionUsingBackground()
	{
		return getPrecisionMethodValue() == PrecisionMethod.MORTENSEN_LOCAL_BACKGROUND_VALUE;
	}

	/**
	 * Gets the precision method used to calculate the precision.
	 *
	 * @return the precision method
	 */
	public PrecisionMethod getPrecisionMethod()
	{
		return filterSettings.getPrecisionMethod();
	}

	/**
	 * Gets the precision method used to calculate the precision.
	 *
	 * @return the precision method
	 */
	public int getPrecisionMethodValue()
	{
		return filterSettings.getPrecisionMethodValue();
	}

	/**
	 * Sets the precision method used to calculate the precision.
	 *
	 * @param precisionMethod
	 *            the new precision method
	 */
	public void setPrecisionMethod(int precisionMethod)
	{
		PrecisionMethod pm = PrecisionMethod.forNumber(precisionMethod);
		if (pm != null)
		{
			setPrecisionMethod(pm);
		}
	}

	/**
	 * Sets the precision method used to calculate the precision.
	 *
	 * @param precisionMethod
	 *            the new precision method
	 */
	public void setPrecisionMethod(PrecisionMethod precisionMethod)
	{
		filterSettings.setPrecisionMethodValue(precisionMethod.getNumber());
		updatePrecisionThreshold();
	}

	// Enable caching the location variance selection
	private class BaseVarianceSelector
	{
		double getLocationVariance(PreprocessedPeakResult peak)
		{
			return 0;
		}
	}

	private class VarianceSelector extends BaseVarianceSelector
	{
		double getLocationVariance(PreprocessedPeakResult peak)
		{
			return peak.getLocationVariance();
		}
	}

	private class VarianceSelector2 extends BaseVarianceSelector
	{
		double getLocationVariance(PreprocessedPeakResult peak)
		{
			return peak.getLocationVariance2();
		}
	}

	private class VarianceSelectorCRLB extends BaseVarianceSelector
	{
		double getLocationVariance(PreprocessedPeakResult peak)
		{
			return peak.getLocationVarianceCRLB();
		}
	}

	/**
	 * Gets the precision method that will be used to produce the precision value for filtering.
	 * <p>
	 * This checks first the direct filter and then the current value for the precision method.
	 * <p>
	 * The value is written into the current calibration to allow the calibration to be used in saved results.
	 *
	 * @return the filter precision method
	 */
	public PrecisionMethod getFilterPrecisionMethod()
	{
		PrecisionMethod m = getFilterPrecisionMethodInternal();
		switch (m)
		{
			case MORTENSEN:
				varianceSelector = new VarianceSelector();
				break;
			case MORTENSEN_LOCAL_BACKGROUND:
				varianceSelector = new VarianceSelector2();
				break;
			case POISSON_CRLB:
				varianceSelector = new VarianceSelectorCRLB();
				break;
			default:
				m = PrecisionMethod.PRECISION_METHOD_NA;
				varianceSelector = new BaseVarianceSelector();
				break;
		}
		calibration.setPrecisionMethod(m);
		return m;
	}

	private PrecisionMethod getFilterPrecisionMethodInternal()
	{
		if (isDirectFilter())
		{
			int flags = directFilter.getValidationFlags();
			if (DirectFilter.areSet(flags, IDirectFilter.V_LOCATION_VARIANCE_CRLB))
				return PrecisionMethod.POISSON_CRLB;
			if (DirectFilter.areSet(flags, IDirectFilter.V_LOCATION_VARIANCE2))
				return PrecisionMethod.MORTENSEN_LOCAL_BACKGROUND;
			if (DirectFilter.areSet(flags, IDirectFilter.V_LOCATION_VARIANCE))
				return PrecisionMethod.MORTENSEN;
		}
		if (precisionThreshold == 0)
			// No simple precision filter
			return PrecisionMethod.PRECISION_METHOD_NA;
		switch (getPrecisionMethodValue())
		{
			case PrecisionMethod.MORTENSEN_VALUE:
				return PrecisionMethod.MORTENSEN;
			case PrecisionMethod.MORTENSEN_LOCAL_BACKGROUND_VALUE:
				return PrecisionMethod.MORTENSEN_LOCAL_BACKGROUND;
			case PrecisionMethod.POISSON_CRLB_VALUE:
				return PrecisionMethod.POISSON_CRLB;
			default:
				break;
		}
		return PrecisionMethod.PRECISION_METHOD_NA;
	}

	/**
	 * Gets the location variance of the result that was used for result filtering. This will use the method matching
	 * that from the last call to {@link #getFilterPrecisionMethod()}. This allows the state to be cached for efficient
	 * selection of the location variance.
	 *
	 * @param peak
	 *            the peak
	 * @return the variance
	 */
	public double getLocationVariance(PreprocessedPeakResult peak)
	{
		return varianceSelector.getLocationVariance(peak);
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
		updateWidthThreshold();
	}

	private void updateWidthThreshold()
	{
		double w = filterSettings.getMaxWidthFactor();
		if (w > 1)
		{
			this.widthFactor = (isTwoAxisGaussian2D) ? w * w : w;
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
		updateMinWidthThreshold();
	}

	private void updateMinWidthThreshold()
	{
		double w = filterSettings.getMinWidthFactor();
		if (w < 1 && w > 0)
		{
			this.minWidthFactor = (isTwoAxisGaussian2D) ? w * w : w;
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
	 * Sets the min Z. If both min and max z are zero then depth filtering is disabled.
	 *
	 * @param minZ
	 *            The minimum z depth
	 */
	public void setMinZ(double minZ)
	{
		filterSettings.setMinZ(minZ);
		updateZFilter();
	}

	/**
	 * Gets the min Z.
	 *
	 * @return the min Z
	 */
	public double getMinZ()
	{
		return filterSettings.getMinZ();
	}

	/**
	 * Sets the max Z. If both min and max z are zero then depth filtering is disabled.
	 *
	 * @param maxZ
	 *            The maximum z depth
	 */
	public void setMaxZ(double maxZ)
	{
		filterSettings.setMaxZ(maxZ);
		updateZFilter();
	}

	/**
	 * Gets the max Z.
	 *
	 * @return the max Z
	 */
	public double getMaxZ()
	{
		return filterSettings.getMaxZ();
	}

	private void updateZFilter()
	{
		// Check if the PSF is 3D but avoid exceptions from creating the z-model
		zEnabled = getPSFTypeValue() == PSFType.ASTIGMATIC_GAUSSIAN_2D_VALUE
		//is3D() // This can throw
				&& (getMaxZ() != 0 || getMinZ() != 0);
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
		invalidateCameraModel();
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
		invalidateCameraModel();
		calibration.setCameraType(cameraType);
		updateCalibration();
	}

	/**
	 * Gets the camera type.
	 *
	 * @return the camera type
	 */
	public CameraType getCameraType()
	{
		return calibration.getCameraType();
	}

	/**
	 * Gets the camera type.
	 *
	 * @return the camera type
	 */
	public int getCameraTypeValue()
	{
		return calibration.getCameraTypeValue();
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
		invalidateCameraModel();
		calibration.setBias(bias);
	}

	/**
	 * @param readNoise
	 *            the camera read noise (used for maximum likelihood estimation)
	 */
	public void setReadNoise(double readNoise)
	{
		invalidateFunctionSolver();
		invalidateCameraModel();
		calibration.setReadNoise(readNoise);
	}

	/**
	 * Sets the quantum efficiency.
	 *
	 * @param quantumEfficiency
	 *            the new quantum efficiency [electron/photon] (used for maximum likelihood estimation)
	 */
	public void setQuantumEfficiency(double quantumEfficiency)
	{
		invalidateFunctionSolver();
		calibration.setQuantumEfficiency(quantumEfficiency);
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
	 * @return the search for the Maximum Likelihood Estimator
	 */
	public int getSearchMethodValue()
	{
		return fitSolverSettings.getSearchMethodValue();
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
		fitSolverSettings.setSearchMethodValue(searchMethod.getNumber());
	}

	/**
	 * @return the line search for the Fast MLE
	 */
	public LineSearchMethod getLineSearchMethod()
	{
		return fitSolverSettings.getLineSearchMethod();
	}

	/**
	 * @return the line search for the Fast MLE
	 */
	public int getLineSearchMethodValue()
	{
		return fitSolverSettings.getLineSearchMethodValue();
	}

	/**
	 * @param lineSearchMethod
	 *            the line search for the Fast MLE
	 */
	public void setLineSearchMethod(int lineSearchMethod)
	{
		LineSearchMethod sm = LineSearchMethod.forNumber(lineSearchMethod);
		if (sm != null)
		{
			setLineSearchMethod(sm);
		}
	}

	/**
	 * @param lineSearchMethod
	 *            the line search for the Fast MLE
	 */
	public void setLineSearchMethod(LineSearchMethod lineSearchMethod)
	{
		invalidateFunctionSolver();
		fitSolverSettings.setLineSearchMethodValue(lineSearchMethod.getNumber());
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
			// If the value is zero then ignore this for convergence. 
			toleranceChecker = new ToleranceChecker(minimiseValue, getIfStrictlyPositive(getRelativeThreshold()),
					getIfStrictlyPositive(getAbsoluteThreshold()),
					getIfStrictlyPositive(getParameterRelativeThreshold()),
					getIfStrictlyPositive(getParameterAbsoluteThreshold()), maxIterations);
		}
		return toleranceChecker;
	}

	private boolean isMinimiseValue()
	{
		switch (getFitSolverValue())
		{
			case FitSolver.BACKTRACKING_FAST_MLE_VALUE:
			case FitSolver.FAST_MLE_VALUE:
				// Maximum likelihood
				return false;

			case FitSolver.LVM_LSE_VALUE:
			case FitSolver.LVM_WLSE_VALUE:
				// Minimise sum-of-squares
				return true;

			case FitSolver.LVM_MLE_VALUE:
				// Minimises the log-likelihood ratio 
				return true;

			case FitSolver.MLE_VALUE:
				// The legacy MLE actual minimises the negative likelihood.
				// This should not matter anyway since the tolerance checker is not used
				// for this fitter.
				return true;

			default:
				throw new IllegalStateException("Unrecognised fit solver: " + getFitSolver());
		}
	}

	/**
	 * Gets the value if positive, otherwise return -1.
	 *
	 * @param value
	 *            the value
	 * @return the value if positive, or -1
	 */
	private double getIfStrictlyPositive(double value)
	{
		return (value > 0) ? value : -1.0;
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

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.Gaussian2DFitConfiguration#validateFit(int, double[], double[], double[])
	 */
	public FitStatus validateFit(int nPeaks, double[] initialParams, double[] params, double[] paramDevs)
	{
		for (int n = 0; n < nPeaks; n++)
		{
			validatePeak(n, initialParams, params, paramDevs);
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
	 * @param paramDevs
	 *            the fitted peak parameter variances (can be null)
	 * @return True if the fit fails the criteria
	 */
	public FitStatus validateFit(double[] initialParams, double[] params, double[] paramDevs)
	{
		return validatePeak(0, initialParams, params, paramDevs);
	}

	/**
	 * Check peak to see if the fit was sensible.
	 *
	 * @param n
	 *            The peak number
	 * @param initialParams
	 *            The initial peak parameters
	 * @param params
	 *            The fitted peak parameters
	 * @param paramDevs
	 *            the fitted peak parameter variances (can be null)
	 * @return True if the fit fails the criteria
	 */
	public FitStatus validatePeak(int n, double[] initialParams, double[] params, double[] paramDevs)
	{
		if (isDirectFilter())
		{
			// Always specify a new result and we have no local background or offset
			PreprocessedPeakResult peak = createPreprocessedPeakResult(0, n, initialParams, params, paramDevs, 0,
					ResultType.NEW, 0, 0, false);
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

		if (zEnabled)
		{
			double z = params[Gaussian2DFunction.Z_POSITION + offset];
			if (z < getMinZ() || z > getMaxZ())
				return setValidationResult(FitStatus.Z_MOVED, z);
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

		double xsd, ysd;
		// Map the width parameters using the z-model if present
		if (getAstigmatismZModel() != null)
		{
			double z = params[Gaussian2DFunction.Z_POSITION + offset];
			xsd = astigmatismZModel.getSx(z);
			ysd = astigmatismZModel.getSy(z);

			// Check widths. This may be the only filter used even if z-fitting 
			// (i.e. a z-depth filter is not used)
			double xFactor = xsd / initialParams[Gaussian2DFunction.X_SD + offset];
			double yFactor = ysd / initialParams[Gaussian2DFunction.Y_SD + offset];
			final double s2 = xFactor * yFactor;

			if (s2 > widthFactor || s2 < minWidthFactor)
			{
				if (log != null)
				{
					log.info("Bad peak %d: Fitted width diverged (x=%gx,y=%gx)\n", n, xFactor, yFactor);
				}
				return setValidationResult(FitStatus.WIDTH_DIVERGED, new double[] { xFactor, yFactor });
			}
		}
		else
		{
			xsd = params[Gaussian2DFunction.X_SD + offset];
			ysd = params[Gaussian2DFunction.Y_SD + offset];
		}

		// Check widths
		if (isXSDFitting())
		{
			boolean badWidth = false;
			double xFactor = 0, yFactor = 0;

			xFactor = xsd / initialParams[Gaussian2DFunction.X_SD + offset];
			if (isTwoAxisGaussian2D)
			{
				yFactor = ysd / initialParams[Gaussian2DFunction.Y_SD + offset];
				final double s2 = xFactor * yFactor;
				badWidth = (s2 > widthFactor || s2 < minWidthFactor);
			}
			else
			{
				badWidth = (xFactor > widthFactor || xFactor < minWidthFactor);
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

		// Check precision. This is above zero if a threshold is present.
		if (precisionThreshold > 0)
		{
			final double variance;
			switch (getPrecisionMethodValue())
			{
				case PrecisionMethod.MORTENSEN_VALUE:
				case PrecisionMethod.MORTENSEN_LOCAL_BACKGROUND_VALUE:
					final double sd = (isTwoAxisGaussian2D) ? Gaussian2DPeakResultHelper.getStandardDeviation(xsd, ysd)
							: xsd;
					variance = getVariance(params[Gaussian2DFunction.BACKGROUND],
							params[Gaussian2DFunction.SIGNAL + offset] * signalToPhotons, sd,
							isPrecisionUsingBackground());
					break;
				case PrecisionMethod.POISSON_CRLB_VALUE:
					variance = getVariance(paramDevs, n);
					break;
				default:
					// This should not happen
					throw new IllegalStateException("Unknown precision method: " + getPrecisionMethod());
			}

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
			if (isMLE)
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
			if (isMLE)
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

	private final static int Y_OFFSET = Gaussian2DFunction.Y_POSITION - Gaussian2DFunction.X_POSITION;

	/**
	 * Gets the variance. This is computed using the mean of the variance for the X and Y parameters.
	 *
	 * @param paramsDev
	 *            the parameter variances
	 * @param n
	 *            the peak number
	 * @return the variance (or zero if there are no deviations)
	 */
	public double getVariance(double[] paramsDev, int n)
	{
		if (paramsDev != null)
		{
			final int offset = n * Gaussian2DFunction.PARAMETERS_PER_PEAK;
			// Scale to nm
			return nmPerPixel * nmPerPixel * (paramsDev[offset + Gaussian2DFunction.X_POSITION] +
					paramsDev[offset + Gaussian2DFunction.Y_POSITION]) / 2.0;
		}
		return 0;
	}

	/**
	 * Update the estimated parameter variance. This method should be run on the deviations produced by fitting.
	 * <p>
	 * If not an EM-CCD then the method does nothing.
	 * <p>
	 * If an EM-CCD the method scales the variance for XY position parameters by a factor of 2. This
	 * allows the computation of the localisation variance using the parameter deviations to match the estimation
	 * formulas of Mortensen. See Mortensen, et al (2010) Nature Methods 7, 377-383, SI 4.3 for assumptions and proof
	 * using MLE. Note that in the EM-CCD case the deviations for the other parameters will be incorrect as they are
	 * computed assuming a Poisson process.
	 *
	 * @param paramsDev
	 *            the params dev
	 */
	public void updateVariance(double[] paramsDev)
	{
		if (emCCD && paramsDev != null)
		{
			for (int i = Gaussian2DFunction.X_POSITION; i < paramsDev.length; i += Gaussian2DFunction.PARAMETERS_PER_PEAK)
			{
				paramsDev[i] *= 2;
				paramsDev[i + Y_OFFSET] *= 2;
			}
		}
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
		double[] paramsDev;
		double xsd, ysd;
		double localBackground;
		boolean existingResult;
		boolean newResult;
		float offsetx;
		float offsety;
		double var, var2, varCRLB;

		DynamicPeakResult(int candidateId, int n, double[] initialParams, double[] params, double[] paramsDev,
				double localBackground, ResultType resultType, float offsetx, float offsety)
		{
			setParameters(candidateId, n, initialParams, params, paramsDev, localBackground, resultType, offsetx,
					offsety);
		}

		DynamicPeakResult()
		{
			var = var2 = varCRLB = -1;
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
		 * @param paramsDev
		 *            the parameter variances (can be null)
		 * @param localBackground
		 *            the local background
		 * @param resultType
		 *            the result type
		 * @param offsetx
		 *            the offsetx
		 * @param offsety
		 *            the offsety
		 */
		void setParameters(int candidateId, int n, double[] initialParams, double[] params, double[] paramsDev,
				double localBackground, ResultType resultType, float offsetx, float offsety)
		{
			this.id = n;
			this.candidateId = candidateId;
			offset = n * Gaussian2DFunction.PARAMETERS_PER_PEAK;
			this.initialParams = initialParams;
			this.params = params;
			this.paramsDev = paramsDev;
			this.localBackground = localBackground;
			this.existingResult = resultType == ResultType.EXISTING;
			this.newResult = resultType == ResultType.NEW;
			this.offsetx = offsetx;
			this.offsety = offsety;
			var = var2 = varCRLB = -1;
			// Map the width parameters using the z-model
			if (getAstigmatismZModel() != null)
			{
				double z = getZ();
				xsd = astigmatismZModel.getSx(z);
				ysd = astigmatismZModel.getSy(z);
			}
			else
			{
				xsd = params[Gaussian2DFunction.X_SD + offset];
				ysd = params[Gaussian2DFunction.Y_SD + offset];
			}
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
			return (float) (params[Gaussian2DFunction.SIGNAL + offset] * signalToPhotons);
		}

		public float getSNR()
		{
			return (float) (params[Gaussian2DFunction.SIGNAL + offset] / getNoise());
		}

		public float getNoise()
		{
			// Comment this out to use the configured local background, i.e. this.localBackground.
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
				var2 = FitConfiguration.this.getVariance(getLocalBackground() * signalToPhotons,
						params[Gaussian2DFunction.SIGNAL + offset] * signalToPhotons, getSD(), true);
			return var2;
		}

		public double getLocationVarianceCRLB()
		{
			if (varCRLB == -1)
				varCRLB = FitConfiguration.this.getVariance(paramsDev, id);
			return varCRLB;
		}

		public float getSD()
		{
			if (isTwoAxisGaussian2D)
				return (float) Gaussian2DPeakResultHelper.getStandardDeviation(xsd, ysd);
			return (float) xsd;
		}

		public float getBackground()
		{
			return (float) params[Gaussian2DFunction.BACKGROUND];
		}

		public float getAmplitude()
		{
			return (float) (params[Gaussian2DFunction.SIGNAL + offset] / (2 * Math.PI * xsd * ysd));
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
			return (float) xsd;
		}

		public float getYSD()
		{
			return (float) ysd;
		}

		public float getXSDFactor()
		{
			return (float) (xsd / initialParams[Gaussian2DFunction.X_SD + offset]);
		}

		public float getYSDFactor()
		{
			return (float) (ysd / initialParams[Gaussian2DFunction.Y_SD + offset]);
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
			// In case of astigmatism
			p[Gaussian2DFunction.X_SD] = xsd;
			p[Gaussian2DFunction.Y_SD] = ysd;
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
	 * @param paramVariances
	 *            the parameter variances (can be null)
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
			double[] params, double[] paramVariances, double localBackground, ResultType resultType, float offsetx,
			float offsety)
	{
		return createPreprocessedPeakResult(candidateId, n, initialParams, params, paramVariances, localBackground,
				resultType, offsetx, offsety, true);
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
	 * @param paramVariances
	 *            the parameter variances (can be null)
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
			double[] params, double[] paramVariances, double localBackground, ResultType resultType, float offsetx,
			float offsety, boolean newObject)
	{
		if (newObject)
			return new DynamicPeakResult(candidateId, n, initialParams, params, paramVariances, localBackground,
					resultType, offsetx, offsety);

		dynamicPeakResult.setParameters(candidateId, n, initialParams, params, paramVariances, localBackground,
				resultType, offsetx, offsety);
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
	 * @param paramVariances
	 *            the parameter variances (can be null)
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
			double[] initialParameters, double[] parameters, double[] paramVariances, double localBackground,
			ResultType resultType, float offsetx, float offsety)
	{
		final int offset = n * Gaussian2DFunction.PARAMETERS_PER_PEAK;
		final double signal = parameters[offset + Gaussian2DFunction.SIGNAL] * signalToPhotons;
		final double b = signalToPhotons *
				((localBackground > 0) ? localBackground : parameters[Gaussian2DFunction.BACKGROUND]);
		final double angle = parameters[offset + Gaussian2DFunction.ANGLE];
		final double x = parameters[offset + Gaussian2DFunction.X_POSITION] + offsetx;
		final double y = parameters[offset + Gaussian2DFunction.Y_POSITION] + offsety;
		final double z = parameters[offset + Gaussian2DFunction.Z_POSITION];
		final double x0 = initialParameters[offset + Gaussian2DFunction.X_POSITION] + offsetx;
		final double y0 = initialParameters[offset + Gaussian2DFunction.Y_POSITION] + offsety;
		double xsd, ysd;
		// Map the width parameters using the z-model
		if (getAstigmatismZModel() != null)
		{
			xsd = astigmatismZModel.getSx(z);
			ysd = astigmatismZModel.getSy(z);
		}
		else
		{
			xsd = parameters[offset + Gaussian2DFunction.X_SD];
			ysd = parameters[offset + Gaussian2DFunction.Y_SD];
		}
		final double xsd0 = initialParameters[offset + Gaussian2DFunction.X_SD];
		final double ysd0 = initialParameters[offset + Gaussian2DFunction.Y_SD];
		final double sd = (isTwoAxisGaussian2D) ? Gaussian2DPeakResultHelper.getStandardDeviation(xsd, ysd) : xsd;
		final double variance = getVariance(0, signal, sd, false);
		final double variance2 = getVariance(b, signal, sd, true);
		final double varianceCRLB = getVariance(paramVariances, n);

		// Q. Should noise be the local background or the estimate from the whole image?

		// This uses the local background if specified or the estimate from the whole image 
		//final double noise = (localBackground > bias)
		//		? PeakResult.localBackgroundToNoise(localBackground - bias, this.gain, this.emCCD) : this.noise;

		// This uses the local fitted background to estimate the noise
		final double noise = (b > 0) ? PeakResultHelper.localBackgroundToNoise(b, 1.0, this.emCCD) : this.noise;
		return new BasePreprocessedPeakResult(frame, n, candidateId, signal, noise, b, angle, x, y, z, x0, y0, xsd, ysd,
				xsd0, ysd0, variance, variance2, varianceCRLB, resultType);
	}

	// The model now entirely defines the width so no unmapping is necessary
	//	/**
	//	 * Unmap the width parameters using the Z model. This assumes the parameters are for a single peak.
	//	 * <p>
	//	 * Note that this is unnecessary if the original widths are known (i.e. at z=0) since they should be identical to
	//	 * the current widths unmapped using the current z.
	//	 *
	//	 * @param params
	//	 *            the params
	//	 */
	//	public void unmapZModel(double[] params)
	//	{
	//		if (getAstigmatismZModel() != null)
	//		{
	//			final double z = params[Gaussian2DFunction.Z_POSITION];
	//			params[Gaussian2DFunction.X_SD] /= astigmatismZModel.getSx(z);
	//			params[Gaussian2DFunction.Y_SD] /= astigmatismZModel.getSy(z);
	//		}
	//	}

	/**
	 * @return Set to true if fitting requires the camera counts, i.e. amplification is explicitly modelled during
	 *         fitting.
	 */
	public boolean isFitCameraCounts()
	{
		// Only the legacy MLE solvers that explicitly model the camera noise require the data
		// and estimate to be in ADUs. This is also true if there is no camera calibration.
		if (fitSolverSettings.getFitSolverValue() == FitSolver.MLE_VALUE ||
				calibration.getCameraTypeValue() == CameraType.CAMERA_TYPE_NA_VALUE)
			return true;
		return false;
	}

	/**
	 * Check if performing raw data fitting. In this mode then no camera model is available to map the data counts to
	 * photons. The gain is assumed to be 1.
	 * 
	 * @return Set to true if raw fitting
	 */
	public boolean isRawFit()
	{
		// This is only true if there is no camera calibration.
		// Otherwise the camera details can be used to convert input data in camera counts
		// to photo-electrons and all fitting is done assuming the signal is photo-electrons.  
		if (calibration.getCameraTypeValue() == CameraType.CAMERA_TYPE_NA_VALUE)
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
		return (isModelCamera() && fitSolverSettings.getFitSolverValue() == FitSolver.MLE_VALUE);
	}

	/**
	 * Checks if is using a maximum likelihood estimator (MLE).
	 *
	 * @return true, if is a MLE
	 */
	public boolean isMLE()
	{
		switch (getFitSolverValue())
		{
			case FitSolver.LVM_MLE_VALUE:
			case FitSolver.MLE_VALUE:
			case FitSolver.FAST_MLE_VALUE:
			case FitSolver.BACKTRACKING_FAST_MLE_VALUE:
				return true;
			default:
				return false;
		}
	}

	/**
	 * The function solver requires a strictly positive function.
	 *
	 * @return true if requires a strictly positive function.
	 */
	public boolean requireStrictlyPositiveFunction()
	{
		// Only the LSE variants can fit negatives. The MLE variants all require a positive function. 
		switch (getFitSolverValue())
		{
			case FitSolver.LVM_LSE_VALUE:
			case FitSolver.LVM_WLSE_VALUE:
				return false;

			case FitSolver.LVM_MLE_VALUE:
			case FitSolver.MLE_VALUE:
			case FitSolver.FAST_MLE_VALUE:
			case FitSolver.BACKTRACKING_FAST_MLE_VALUE:
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

		// Special case where the data and estimate are in counts but
		// the function solver requires the function to output an expected 
		// number of photons.
		boolean doScaling = getFitSolverValue() == FitSolver.MLE_VALUE;
		assert !doScaling || functionSolver instanceof MaximumLikelihoodFitter;

		if (precomputedFunctionValues != null)
		{
			if (doScaling)
			{
				// Scale the pre-computed function. It is assumes that this was
				// set using output fitted parameters (which would be unscaled).
				// This may be reused by the calling code so don't destroy it.
				double[] f = new double[precomputedFunctionValues.length];
				for (int i = precomputedFunctionValues.length; i-- > 0;)
					f[i] = precomputedFunctionValues[i] * signalToPhotons;
				precomputedFunctionValues = f;
			}

			functionSolver.setGradientFunction((GradientFunction) PrecomputedFunctionFactory
					.wrapFunction(gaussianFunction, precomputedFunctionValues));
			precomputedFunctionValues = null;
		}
		if (functionSolver.isWeighted())
		{
			functionSolver.setWeights(observationWeights);
		}

		if (doScaling)
		{
			// Scale the background and signal estimates
			int[] indices = new int[1 + gaussianFunction.getNPeaks()];
			indices[0] = Gaussian2DFunction.BACKGROUND;
			for (int i = 1; i < indices.length; i++)
				indices[i] = (i - 1) * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.SIGNAL;

			// signalToPhotons should be 1.0/gain
			return new MLEScaledFunctionSolver((MaximumLikelihoodFitter) functionSolver, signalToPhotons, indices);
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

		if (getFitSolverValue() == FitSolver.MLE_VALUE)
		{
			// Only support CCD/EM-CCD at the moment
			if (!calibration.isCCDCamera())
			{
				throw new IllegalStateException("CCD/EM-CCD camera is required for fit solver: " + getFitSolver());
			}

			// This requires the gain
			if (gain <= 0)
			{
				throw new IllegalStateException("The gain is required for fit solver: " + getFitSolver());
			}

			MaximumLikelihoodFitter.SearchMethod searchMethod = convertSearchMethod();

			// Only the Poisson likelihood function supports gradients
			if (searchMethod.usesGradients() && isModelCamera())
			{
				throw new IllegalStateException(String.format(
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
				throw new IllegalStateException("The amplification is required for the fit solver: " + getFitSolver());
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

		switch (getFitSolverValue())
		{
			case FitSolver.LVM_LSE_VALUE:
				solver = new LSELVMSteppingFunctionSolver(gaussianFunction, tc, bounds);
				break;

			case FitSolver.LVM_MLE_VALUE:
				checkCameraCalibration();
				solver = new MLELVMSteppingFunctionSolver(gaussianFunction, tc, bounds);
				break;

			case FitSolver.LVM_WLSE_VALUE:
				checkCameraCalibration();
				solver = new WLSELVMSteppingFunctionSolver(gaussianFunction, tc, bounds);
				break;

			case FitSolver.FAST_MLE_VALUE:
				checkCameraCalibration();
				// This may throw a class cast exception if the function does not support
				// the Gradient2Function interface
				solver = new FastMLESteppingFunctionSolver((Gradient2Function) gaussianFunction, tc, bounds);
				break;

			case FitSolver.BACKTRACKING_FAST_MLE_VALUE:
				checkCameraCalibration();
				solver = new BacktrackingFastMLESteppingFunctionSolver((Gradient2Function) gaussianFunction, tc,
						bounds);
				break;

			default:
				throw new IllegalStateException("Unknown fit solver: " + getFitSolver());
		}

		if (solver instanceof LVMSteppingFunctionSolver)
		{
			((LVMSteppingFunctionSolver) solver).setInitialLambda(getLambda());
		}
		else if (solver instanceof FastMLESteppingFunctionSolver)
		{
			((FastMLESteppingFunctionSolver) solver).setLineSearchMethod(convertLineSearchMethod());
		}
		return solver;
	}

	private gdsc.smlm.fitting.nonlinear.MaximumLikelihoodFitter.SearchMethod convertSearchMethod()
	{
		return FitProtosHelper.convertSearchMethod(getSearchMethod());
	}

	private gdsc.smlm.fitting.nonlinear.FastMLESteppingFunctionSolver.LineSearchMethod convertLineSearchMethod()
	{
		return FitProtosHelper.convertLineSearchMethod(getLineSearchMethod());
	}

	private void checkCameraCalibration()
	{
		if (!calibration.hasCameraCalibration())
			throw new IllegalStateException("The camera calibration is required for fit solver: " + getFitSolver());

		switch (getCameraTypeValue())
		{
			// CCD/EMCCD requires gain and bias (but bias could be zero)
			case CameraType.CCD_VALUE:
			case CameraType.EMCCD_VALUE:

				// sCMOS requires per-pixel bias, gain and read noise (var/gain^2)
			case CameraType.SCMOS_VALUE:

				// Handle the camera checks within getCameraModel(). This throws if the 
				// camera model is invalid
				getCameraModel();
				break;

			case CameraType.CAMERA_TYPE_NA_VALUE:
			default:
				throw new IllegalStateException(
						"Unrecognised camera type for for fit solver: " + getFitSolver() + ": " + getCameraType());
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
		return (TextUtils.isNullOrEmpty(s)) ? "" : s;
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
		float minZ = (float) getMinZ();
		float maxZ = (float) getMaxZ();

		switch (getPrecisionMethodValue())
		{
			case PrecisionMethod.MORTENSEN_VALUE:
				return new MultiFilter(signal, snr, minWidth, maxWidth, shift, eshift, precision, minZ, maxZ);
			case PrecisionMethod.MORTENSEN_LOCAL_BACKGROUND_VALUE:
				return new MultiFilter2(signal, snr, minWidth, maxWidth, shift, eshift, precision, minZ, maxZ);
			case PrecisionMethod.POISSON_CRLB_VALUE:
				return new MultiFilterCRLB(signal, snr, minWidth, maxWidth, shift, eshift, precision, minZ, maxZ);
			default:
				throw new IllegalStateException("Unknown precision method: " + getPrecisionMethod());
		}
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
		// Handle switching to a 2 axis width filter
		if (isTwoAxisGaussian2D)
			flags |= IDirectFilter.XY_WIDTH;
		filterSetupFlags = flags;
		this.filterSetupData = null;

		if (directFilter != null)
		{
			directFilter.setup(flags);
		}
		else
		{
			widthEnabled = !DirectFilter.areSet(flags, DirectFilter.NO_WIDTH);
			if (DirectFilter.areSet(flags, DirectFilter.NO_SHIFT))
			{
				offset = Float.POSITIVE_INFINITY;
			}
			else
			{
				double shiftFactor = getCoordinateShiftFactor();
				offset = (float) ((shiftFactor > 0) ? shiftFactor * shiftFactor : Float.POSITIVE_INFINITY);
			}
			varianceThreshold = (precisionThreshold > 0) ? precisionThreshold : Double.POSITIVE_INFINITY;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IDirectFilter#setup(gdsc.smlm.results.filter.FilterSetupData[])
	 */
	public void setup(int flags, FilterSetupData... filterSetupData)
	{
		// Handle switching to a 2 axis width filter
		if (isTwoAxisGaussian2D)
			flags |= IDirectFilter.XY_WIDTH;
		filterSetupFlags = flags;
		this.filterSetupData = filterSetupData;

		if (directFilter != null)
		{
			directFilter.setup(flags, filterSetupData);
		}
		else
		{
			// Note: These variables should be copied in copySettings(...)
			widthEnabled = !DirectFilter.areSet(flags, DirectFilter.NO_WIDTH);
			if (DirectFilter.areSet(flags, DirectFilter.NO_SHIFT))
			{
				offset = Float.POSITIVE_INFINITY;
			}
			else
			{
				double shiftFactor = getCoordinateShiftFactor();
				for (int i = filterSetupData.length; i-- > 0;)
				{
					if (filterSetupData[i] instanceof ShiftFilterSetupData)
					{
						double shift = ((ShiftFilterSetupData) filterSetupData[i]).shift;
						if (shift > 0)
						{
							double widthMax = getWidthMax();
							if (widthMax > 0)
							{
								shiftFactor = shift / widthMax;
							}
						}
						break;
					}
				}
				offset = (float) ((shiftFactor > 0) ? shiftFactor * shiftFactor : Float.POSITIVE_INFINITY);
			}
			varianceThreshold = (precisionThreshold > 0) ? precisionThreshold : Double.POSITIVE_INFINITY;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IDirectFilter#getFilterSetupFlags()
	 */
	public int getFilterSetupFlags() throws IllegalStateException
	{
		// Cached for speed
		return filterSetupFlags;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IDirectFilter#getFilterSetupData()
	 */
	public FilterSetupData[] getFilterSetupData() throws IllegalStateException
	{
		// Cached for speed
		return filterSetupData;
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
	 * @see gdsc.smlm.results.filter.IDirectFilter#getValidationFlags()
	 */
	public int getValidationFlags()
	{
		if (directFilter != null)
			return directFilter.getValidationFlags();
		// Q. Is this necessary? See doValidate below ...
		if (isDisableSimpleFilter())
			return 0;
		// These could be conditional on the filter settings. For now just set them all.
		int flags = V_PHOTONS | V_SNR | V_X_RELATIVE_SHIFT | V_Y_RELATIVE_SHIFT;
		if (widthEnabled)
		{
			flags |= V_X_SD_FACTOR;
			if (isTwoAxisGaussian2D)
				flags |= V_Y_SD_FACTOR;
		}
		switch (getPrecisionMethodValue())
		{
			case PrecisionMethod.MORTENSEN_VALUE:
				flags |= V_LOCATION_VARIANCE;
				break;
			case PrecisionMethod.MORTENSEN_LOCAL_BACKGROUND_VALUE:
				flags |= V_LOCATION_VARIANCE2;
				break;
			case PrecisionMethod.POISSON_CRLB_VALUE:
				flags |= V_LOCATION_VARIANCE_CRLB;
				break;
		}
		return flags;
	}

	/**
	 * Do the validation. Either use the direct filter or perform the same function using the filter configuration.
	 *
	 * @param peak
	 *            the peak
	 * @return the result
	 * @see gdsc.smlm.results.filter.IDirectFilter#validate(gdsc.smlm.results.filter.PreprocessedPeakResult)
	 */
	public int doValidate(PreprocessedPeakResult peak)
	{
		if (directFilter != null)
			return directFilter.validate(peak);

		// Q. Is this necessary? Simple filtering is to support turning off 
		// filtering in the validatePeak(...) method. Set a debug point to check
		// if this is used.
		if (isDisableSimpleFilter())
			return 0;

		// Do filtering 
		if (peak.getSignal() < getMinPhotons())
			return V_PHOTONS;
		if (peak.getSNR() < getSignalStrength())
			return V_SNR;
		if (widthEnabled)
		{
			// Handle switching to a 2 axis width filter. 
			// Note this will ignore any flags passed to the IDirectFilter#setup methods
			// to control the filter.
			if (isTwoAxisGaussian2D)
			{
				final float s2 = peak.getXSDFactor() * peak.getYSDFactor();
				if (s2 > widthFactor || s2 < minWidthFactor)
					return V_X_SD_FACTOR | V_Y_SD_FACTOR;
			}
			else
			{
				final double s = peak.getXSDFactor();
				if (s > widthFactor || s < minWidthFactor)
					return V_X_SD_FACTOR;
			}
		}
		if (peak.getXRelativeShift2() > offset)
			return V_X_RELATIVE_SHIFT;
		if (peak.getYRelativeShift2() > offset)
			return V_Y_RELATIVE_SHIFT;
		// Do not support Euclidian shift
		//if (peak.getXRelativeShift2() + peak.getYRelativeShift2() > offset)
		//	return V_X_RELATIVE_SHIFT | V_Y_RELATIVE_SHIFT;

		if (zEnabled)
		{
			double z = peak.getZ();
			if (z < getMinZ() || z > getMaxZ())
				return V_Z;
		}

		switch (getPrecisionMethodValue())
		{
			case PrecisionMethod.MORTENSEN_VALUE:
				if (peak.getLocationVariance() > varianceThreshold)
					return V_LOCATION_VARIANCE;
				break;
			case PrecisionMethod.MORTENSEN_LOCAL_BACKGROUND_VALUE:
				if (peak.getLocationVariance2() > varianceThreshold)
					return V_LOCATION_VARIANCE2;
				break;
			case PrecisionMethod.POISSON_CRLB_VALUE:
				if (peak.getLocationVarianceCRLB() > varianceThreshold)
					return V_LOCATION_VARIANCE_CRLB;
				break;
			default:
				throw new IllegalStateException("Unknown precision method: " + getPrecisionMethod());
		}

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
	 * Sets the precomputed function values. This is combined with the configured Gaussian function and passed to the
	 * function solver returned from {@link #getFunctionSolver()}. The precomputed function values are then reset to
	 * null so it must be set each time the function solver is accessed.
	 *
	 * @param precomputedFunctionValues
	 *            the new precomputed function values
	 */
	public void setPrecomputedFunctionValues(double[] precomputedFunctionValues)
	{
		this.precomputedFunctionValues = precomputedFunctionValues;
	}

	/**
	 * Sets the observation weights. These are passed to the function solver if it supports weights. These must be set
	 * each time the function solver from {@link #getFunctionSolver()} will be used on new data.
	 *
	 * @param observationWeights
	 *            the new observation weights
	 */
	public void setObservationWeights(double[] observationWeights)
	{
		this.observationWeights = observationWeights;
	}

	/**
	 * Gets the observation weights.
	 *
	 * @return the observation weights
	 */
	public double[] getObservationWeights()
	{
		return observationWeights;
	}

	/**
	 * Call this when a property changes that will change the camera model, e.g. bias, gain, read noise, camera type.
	 */
	private void invalidateCameraModel()
	{
		setCameraModel(null);
	}

	/**
	 * Sets the camera model. This must be set if a sCMOS camera type is used.
	 *
	 * @param cameraModel
	 *            the new camera model
	 */
	public void setCameraModel(CameraModel cameraModel)
	{
		invalidateFunctionSolver();
		// Use this to set the bias and gain
		if (cameraModel != null && cameraModel.isPerPixelModel())
		{
			calibration.clearGlobalCameraSettings();

			// Trigger an update to the calibration used for validation.
			updateCalibration();
		}
		this.cameraModel = cameraModel;
	}

	/**
	 * Return true if the camera type requires a per-pixel camera model
	 *
	 * @return true, if successful
	 */
	public boolean isPerPixelCameraType()
	{
		switch (getCameraTypeValue())
		{
			case CameraType.CAMERA_TYPE_NA_VALUE:
			case CameraType.CCD_VALUE:
			case CameraType.EMCCD_VALUE:
				return false;

			case CameraType.SCMOS_VALUE:
				return true;

			default:
				throw new IllegalStateException("Unknown camera type: " + getCameraType());
		}
	}

	/**
	 * Gets the camera model.
	 *
	 * @return the camera model
	 */
	public CameraModel getCameraModel()
	{
		if (cameraModel == null)
		{
			switch (getCameraTypeValue())
			{
				case CameraType.CAMERA_TYPE_NA_VALUE:
					// We can support this by doing nothing to pixels values
					cameraModel = new NullCameraModel();
					break;

				case CameraType.CCD_VALUE:
				case CameraType.EMCCD_VALUE:
					float bias = (float) calibration.getBias();
					float gain = (float) calibration.getCountPerPhoton();
					float variance = (float) Maths.pow2(calibration.getReadNoise());
					// This will throw an exception if the calibration is invalid
					cameraModel = new FixedPixelCameraModel(bias, gain, variance);
					break;

				case CameraType.SCMOS_VALUE:
					// Impossible to create a per-pixel model with the current configuration.
					// The model must be passed in.

				default:
					throw new IllegalStateException("No camera model for camera type: " + getCameraType());
			}
		}
		return cameraModel;
	}

	/**
	 * Sets the camera model name. This should contain all the information required to load the camera model, e.g. in
	 * the case of a per-pixel camera model for sCMOS cameras.
	 * <p>
	 * This settings is saved to the underlying configuration. If a camera model is used (e.g. for sCMOS camera) then
	 * {@link #setCameraModel(CameraModel)} should be called after setting the new camera model name.
	 *
	 * @param cameraModelName
	 *            the new camera model name
	 */
	public void setCameraModelName(String cameraModelName)
	{
		if (cameraModelName == null || !cameraModelName.equals(getCameraModelName()))
			invalidateCameraModel();
		calibration.setCameraModelName(cameraModelName);
	}

	/**
	 * Gets the camera model name.
	 *
	 * @return the camera model name
	 */
	public String getCameraModelName()
	{
		return calibration.getCameraModelName();
	}

	/**
	 * Sets the PSF model name. This should contain all the information required to load the PSF model, e.g. in
	 * the case of an astigmatic Gaussian 2D PSF.
	 * <p>
	 * This settings is saved to the underlying configuration. If a PSF model is used (e.g. for an astigmatic Gaussian
	 * 2D PSF) then {@link #setAstigmatismModel(AstigmatismModel)} should be called after setting the new PSF model
	 * name.
	 *
	 * @param psfModelName
	 *            the new PSF model name
	 */
	public void setPSFModelName(String psfModelName)
	{
		psf.setModelName(psfModelName);
	}

	/**
	 * Gets the PSF model name.
	 *
	 * @return the PSF model name
	 */
	public String getPSFModelName()
	{
		return psf.getModelName();
	}

	/**
	 * Sets the astigmatism model. This has the effect of changing the PSF to an astigmatic Gaussian 2D function.
	 *
	 * @param model
	 *            the new astigmatism model
	 * @throws ConfigurationException
	 *             if the model pixel pitch does not match the calibration
	 * @throws ConversionException
	 *             if the model cannot be converted to pixel units
	 */
	public void setAstigmatismModel(AstigmatismModel model) throws ConfigurationException, ConversionException
	{
		// Check the calibration
		if (DoubleEquality.relativeError(model.getNmPerPixel(), calibration.getNmPerPixel()) > 1e-3)
		{
			throw new ConfigurationException(
					String.format("The astigmatism model pixel pitch (%s) does not match the calibration (%s)",
							model.getNmPerPixel(), calibration.getNmPerPixel()));
		}

		// Convert to pixels
		model = PSFProtosHelper.convert(model, DistanceUnit.PIXEL, DistanceUnit.PIXEL);

		// Create the working model
		astigmatismZModel = HoltzerAstigmatismZModel.create(model.getS0X(), model.getS0Y(), model.getGamma(),
				model.getD(), model.getAx(), model.getBx(), model.getAy(), model.getBy());

		// Store the parameters in the PSF
		String modelName = getPSFModelName();
		psf.clear().mergeFrom(PSFProtosHelper.createPSF(model, DistanceUnit.PIXEL, DistanceUnit.PIXEL));
		psf.setModelName(modelName);
		updatePSF(false);
	}

	/**
	 * Gets the astigmatism Z model. This is only valid if the PSF type is an astigmatic Gaussian 2D and the parameters
	 * are correctly configured.
	 *
	 * @return the astigmatism Z model (or null)
	 * @throws ConfigurationException
	 *             if the model cannot be created from the parameters
	 */
	public AstigmatismZModel getAstigmatismZModel() throws ConfigurationException
	{
		if (getPSFTypeValue() == PSFType.ASTIGMATIC_GAUSSIAN_2D_VALUE)
		{
			if (astigmatismZModel == null)
			{
				// Use the helper to convert the PSF parameters back to a model
				AstigmatismModel model = PSFProtosHelper.createModel(getPSF(), DistanceUnit.PIXEL, DistanceUnit.PIXEL,
						calibration.getNmPerPixel());
				astigmatismZModel = HoltzerAstigmatismZModel.create(model.getS0X(), model.getS0Y(), model.getGamma(),
						model.getD(), model.getAx(), model.getBx(), model.getAy(), model.getBy());
			}
			return astigmatismZModel;
		}
		return null;
	}

	/**
	 * Checks if the PSF is 3D. Currently only an Astigmatic 2D Gaussian is supported. A check is made for a valid
	 * astigmatism z-model.
	 *
	 * @return true, if is 3D
	 * @throws ConfigurationException
	 *             if the 3D model cannot be created
	 */
	public boolean is3D() throws ConfigurationException
	{
		return (getPSFTypeValue() == PSFType.ASTIGMATIC_GAUSSIAN_2D_VALUE) ? getAstigmatismZModel() != null : false;
	}
}