/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */

package uk.ac.sussex.gdsc.smlm.engine;

import uk.ac.sussex.gdsc.core.data.NotImplementedException;
import uk.ac.sussex.gdsc.core.data.utils.ConversionException;
import uk.ac.sussex.gdsc.core.logging.LoggerUtils;
import uk.ac.sussex.gdsc.core.match.FractionalAssignment;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraType;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.ConfigurationException;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FilterSettings;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FitSettings;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FitSolver;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FitSolverSettings;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.LineSearchMethod;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.PrecisionMethod;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.SearchMethod;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.AstigmatismModel;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFParameter;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFParameterUnit;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.data.config.PsfHelper;
import uk.ac.sussex.gdsc.smlm.data.config.PsfProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.fitting.FitStatus;
import uk.ac.sussex.gdsc.smlm.fitting.FunctionSolver;
import uk.ac.sussex.gdsc.smlm.fitting.Gaussian2DFitConfiguration;
import uk.ac.sussex.gdsc.smlm.fitting.MleScaledFunctionSolver;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.BacktrackingFastMleSteppingFunctionSolver;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.BaseFunctionSolver;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.FastMleSteppingFunctionSolver;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.LseLvmSteppingFunctionSolver;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.LvmSteppingFunctionSolver;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.MaximumLikelihoodFitter;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.MleLvmSteppingFunctionSolver;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.ParameterBounds;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.SteppingFunctionSolver;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.ToleranceChecker;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.WLseLvmSteppingFunctionSolver;
import uk.ac.sussex.gdsc.smlm.function.Gradient2Function;
import uk.ac.sussex.gdsc.smlm.function.GradientFunction;
import uk.ac.sussex.gdsc.smlm.function.OffsetFunctionFactory;
import uk.ac.sussex.gdsc.smlm.function.gaussian.AstigmatismZModel;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.smlm.function.gaussian.HoltzerAstigmatismZModel;
import uk.ac.sussex.gdsc.smlm.model.camera.CameraModel;
import uk.ac.sussex.gdsc.smlm.model.camera.CcdCameraModel;
import uk.ac.sussex.gdsc.smlm.model.camera.EmCcdCameraModel;
import uk.ac.sussex.gdsc.smlm.model.camera.NullCameraModel;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.filter.BasePreprocessedPeakResult;
import uk.ac.sussex.gdsc.smlm.results.filter.BasePreprocessedPeakResult.ResultType;
import uk.ac.sussex.gdsc.smlm.results.filter.DirectFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.Filter;
import uk.ac.sussex.gdsc.smlm.results.filter.FilterSetupData;
import uk.ac.sussex.gdsc.smlm.results.filter.FilterType;
import uk.ac.sussex.gdsc.smlm.results.filter.IDirectFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.MultiFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.MultiFilter2;
import uk.ac.sussex.gdsc.smlm.results.filter.MultiFilterCrlb;
import uk.ac.sussex.gdsc.smlm.results.filter.PreprocessedPeakResult;
import uk.ac.sussex.gdsc.smlm.results.filter.ShiftFilterSetupData;

import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Specifies the fitting configuration.
 */
public class FitConfiguration implements Cloneable, IDirectFilter, Gaussian2DFitConfiguration {
  private FitSettings.Builder fitSettings;

  // Extract the settings for convenience
  private CalibrationWriter calibration;
  private PSF.Builder psf;
  private FilterSettings.Builder filterSettings;
  private FitSolverSettings.Builder fitSolverSettings;

  private Logger log;

  private boolean computeDeviations;
  private int flags;
  private AstigmatismZModel astigmatismZModel;
  private double coordinateShift = 1;
  private int fitRegionWidth;
  private int fitRegionHeight;
  private double coordinateOffset = 0.5;
  private double minSignal;
  private double precisionThreshold;
  private boolean isTwoAxisGaussian2D;
  private double nmPerPixel;
  private double gain;
  private double signalToPhotons;
  private boolean emCcd;
  private boolean isMle;
  private double noise;
  private double minWidthFactor = 0.5;
  /**
   * The width factor. This is the squared width if a 2 axis PSF
   */
  private double widthFactor = 2;
  private boolean computeResiduals = true;
  private boolean zEnabled;

  // Options for clamping
  private double[] clampValues;
  /** The number of peaks for the clamp values. */
  private int clampPeakCount;
  private ParameterBounds bounds;

  private ToleranceChecker toleranceChecker;
  private Gaussian2DFunction gaussianFunction;
  private BaseFunctionSolver functionSolver;

  private DynamicPeakResult dynamicPeakResult = new DynamicPeakResult();

  // Support using a smart filter and disabling the simple filtering
  private DirectFilter directFilter;
  private int filterResult;
  private boolean widthEnabled;
  private float shiftOffset;
  private double varianceThreshold;
  private int filterSetupFlags;
  private FilterSetupData[] filterSetupData;

  private double[] precomputedFunctionValues;
  private double[] observationWeights;
  private CameraModel cameraModel;

  private BaseVarianceSelector varianceSelector = new BaseVarianceSelector();

  /**
   * Instantiates a new fit configuration.
   */
  public FitConfiguration() {
    this(FitProtosHelper.defaultFitSettings, CalibrationProtosHelper.defaultCalibration,
        PsfProtosHelper.defaultOneAxisGaussian2DPSF);
  }

  /**
   * Instantiates a new fit configuration.
   *
   * @param fitSettings the fit settings
   * @param calibration the calibration
   * @param psf the psf
   */
  public FitConfiguration(FitSettings fitSettings, Calibration calibration, PSF psf) {
    if (fitSettings == null) {
      throw new IllegalArgumentException("FitSettings is null");
    }
    if (calibration == null) {
      throw new IllegalArgumentException("Calibration is null");
    }
    if (psf == null) {
      throw new IllegalArgumentException("PSF is null");
    }
    init(fitSettings.toBuilder(), calibration.toBuilder(), psf.toBuilder());
  }

  /**
   * Instantiates a new fit configuration.
   *
   * @param fitSettings the fit settings
   * @param calibration the calibration
   * @param psf the psf
   */
  public FitConfiguration(FitSettings.Builder fitSettings, Calibration.Builder calibration,
      PSF.Builder psf) {
    if (fitSettings == null) {
      throw new IllegalArgumentException("FitSettings is null");
    }
    if (calibration == null) {
      throw new IllegalArgumentException("Calibration is null");
    }
    if (psf == null) {
      throw new IllegalArgumentException("PSF is null");
    }
    init(fitSettings, calibration, psf);
  }

  /**
   * Instantiates a new fit configuration. Does not check for null objects.
   *
   * @param fitSettings the fit settings
   * @param calibration the calibration
   * @param psf the psf
   * @param dummy the dummy parameter
   */
  FitConfiguration(FitSettings.Builder fitSettings, Calibration.Builder calibration,
      PSF.Builder psf, boolean dummy) {
    init(fitSettings, calibration, psf);
  }

  /**
   * Initialise the instance.
   *
   * @param fitSettings the fit settings
   * @param calibration the calibration
   * @param psf the psf
   */
  private void init(FitSettings.Builder fitSettings, Calibration.Builder calibration,
      PSF.Builder psf) {
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
   * @param fitSettings the fit settings
   */
  void updateFitSettings(FitSettings.Builder fitSettings) {
    fitSolverSettings = fitSettings.getFitSolverSettingsBuilder();
    filterSettings = fitSettings.getFilterSettingsBuilder();
    updateFitSolverSettings();
    updateFilterSettings();
  }

  /**
   * Ensure that the internal state of the object is initialised. This is used after deserialisation
   * since some state is not saved but restored from other property values.
   */
  public void initialiseState() {
    // Create the state using the settings
    if (dynamicPeakResult == null) {
      dynamicPeakResult = new DynamicPeakResult();
    }
    updateCalibration();
    updatePsf(true);
    updateFitSolverSettings();
    updateFilterSettings();
  }

  /**
   * Gets the fit settings.
   *
   * @return the fit settings
   */
  public FitSettings getFitSettings() {
    return fitSettings.build();
  }

  /**
   * Merge fit settings.
   *
   * @param fitSettings the fit settings
   */
  public void mergeFitSettings(FitSettings fitSettings) {
    this.fitSettings.mergeFrom(fitSettings);
    initialiseState();
  }

  /**
   * Sets the fit settings.
   *
   * @param fitSettings the new fit settings
   */
  public void setFitSettings(FitSettings fitSettings) {
    this.fitSettings.clear().mergeFrom(fitSettings);
    initialiseState();
  }

  /**
   * Gets the calibration.
   *
   * @return the calibration
   */
  public Calibration getCalibration() {
    return calibration.getCalibration();
  }

  /**
   * Gets a reference to the current calibration writer.
   *
   * @return the calibration writer
   */
  CalibrationWriter getCalibrationWriterReference() {
    return calibration;
  }

  /**
   * Gets a new calibration writer. Any changes to the writer must be saved using
   * {@link #setCalibration(uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration)}.
   *
   * @return the calibration writer
   */
  public CalibrationWriter getCalibrationWriter() {
    return new CalibrationWriter(calibration.getCalibration());
  }

  /**
   * Gets a new calibration reader.
   *
   * @return the calibration reader
   */
  public CalibrationReader getCalibrationReader() {
    return new CalibrationReader(calibration.getCalibrationOrBuilder());
  }

  /**
   * Merge the calibration.
   *
   * @param calibration the new calibration
   */
  public void mergeCalibration(Calibration calibration) {
    this.calibration.mergeCalibration(calibration);
    updateCalibration();
    invalidateCameraModel();
  }

  /**
   * Sets the calibration.
   *
   * @param calibration the new calibration
   */
  public void setCalibration(Calibration calibration) {
    this.calibration.setCalibration(calibration);
    updateCalibration();
    invalidateCameraModel();
  }

  /**
   * Update calibration used for signal and precision filtering (i.e. nmPerPixel, gain, emCCD
   * (camera type)).
   */
  private void updateCalibration() {
    // This uses the camera calibration
    invalidateFunctionSolver();

    nmPerPixel = calibration.getNmPerPixel();
    gain = calibration.getCountPerPhoton();
    emCcd = calibration.isEmCcd();

    if (isRawFit()) {
      // No camera calibration so assume raw data is in photons
      gain = 1;
    }

    if (isFitCameraCounts()) {
      signalToPhotons = 1.0 / gain;
    } else {
      signalToPhotons = 1;
    }

    updateMinSignal();
    updatePrecisionThreshold();
  }

  /**
   * Gets the PSF.
   *
   * @return the PSF
   */
  public PSF getPsf() {
    return psf.build();
  }

  /**
   * Merge the PSF.
   *
   * @param psf the new PSF
   */
  public void mergePsf(PSF psf) {
    this.psf.mergeFrom(psf);
    updatePsf(true);
  }

  /**
   * Sets the psf.
   *
   * @param psf the new psf
   */
  public void setPsf(PSF psf) {
    this.psf.clear().mergeFrom(psf);
    updatePsf(true);
  }

  private void updatePsf(boolean resetAstigmatismModel) {
    invalidateGaussianFunction();

    // Reset the astigmatism model. It will be dynamically created from the PSF settings.
    if (resetAstigmatismModel) {
      astigmatismZModel = null;
    }

    int paramCount;
    final PSFType psfType = psf.getPsfType();
    switch (psfType) {
      case ASTIGMATIC_GAUSSIAN_2D:
        flags = GaussianFunctionFactory.FIT_ERF_ASTIGMATISM;
        // nParams = 2; // Store Sx and Sy
        paramCount = 8; // The PSF stores the full astigmatism model
        break;
      case ONE_AXIS_GAUSSIAN_2D:
        if (isFixedPsf()) {
          flags = GaussianFunctionFactory.FIT_ERF_FIXED;
        } else {
          flags = GaussianFunctionFactory.FIT_ERF_CIRCLE;
        }
        paramCount = 1;
        break;
      case TWO_AXIS_AND_THETA_GAUSSIAN_2D:
        flags = GaussianFunctionFactory.FIT_ELLIPTICAL;
        paramCount = 3;
        break;
      case TWO_AXIS_GAUSSIAN_2D:
        flags = GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE;
        paramCount = 2;
        break;
      default:
        throw new IllegalStateException("FitSettings must be a Gaussian 2D PSF");
    }
    isTwoAxisGaussian2D = PsfHelper.isTwoAxisGaussian2D(psfType);

    final boolean changed = psf.getParametersCount() > paramCount;
    if (changed) {
      while (psf.getParametersCount() > paramCount) {
        psf.removeParameters(psf.getParametersCount() - 1);
      }

      // Updated names when changing from two-axis to one axis
      if (psf.getParametersCount() == 1) {
        psf.getParametersBuilder(PsfHelper.INDEX_SX)
            .setName(PsfProtosHelper.defaultOneAxisGaussian2DPSF.getParameters(PsfHelper.INDEX_SX)
                .getName());
      }
    }

    // Ensure we have enough parameters
    if (psf.getParametersCount() == 0 && paramCount > 0) {
      // Create a dummy Sx
      final PSFParameter.Builder p = psf.addParametersBuilder();
      p.setName(
          PsfProtosHelper.defaultOneAxisGaussian2DPSF.getParameters(PsfHelper.INDEX_SX).getName());
      p.setValue(1);
      p.setUnit(PSFParameterUnit.DISTANCE);
    }
    if (psf.getParametersCount() == 1 && paramCount > 1) {
      // Rename S to Sx
      psf.getParametersBuilder(PsfHelper.INDEX_SX).setName(
          PsfProtosHelper.defaultTwoAxisGaussian2DPSF.getParameters(PsfHelper.INDEX_SX).getName());

      // Duplicate the Sx to Sy
      final PSFParameter.Builder p = psf.addParametersBuilder();
      p.setName(
          PsfProtosHelper.defaultTwoAxisGaussian2DPSF.getParameters(PsfHelper.INDEX_SY).getName());
      p.setValue(psf.getParameters(PsfHelper.INDEX_SX).getValue());
      p.setUnit(PSFParameterUnit.DISTANCE);
    }
    if (psf.getParametersCount() == 2 && paramCount > 2) {
      // Create a dummy angle
      final PSFParameter.Builder p = psf.addParametersBuilder();
      p.setName(PsfProtosHelper.defaultTwoAxisAndThetaGaussian2DPSF
          .getParameters(PsfHelper.INDEX_THETA).getName());
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
  public FitSolverSettings getFitSolverSettings() {
    return fitSolverSettings.build();
  }

  /**
   * Merge the FitSolverSettings.
   *
   * @param fitSolverSettings the new FitSolverSettings
   */
  public void mergeFitSolverSettings(FitSolverSettings fitSolverSettings) {
    this.fitSolverSettings.mergeFrom(fitSolverSettings);
    updateFitSolverSettings();
  }

  /**
   * Sets the fit solver settings.
   *
   * @param fitSolverSettings the new fit solver settings
   */
  public void setFitSolverSettings(FitSolverSettings fitSolverSettings) {
    fitSettings.setFitSolverSettings(fitSolverSettings);
    this.fitSolverSettings = fitSettings.getFitSolverSettingsBuilder();
    // this.fitSolverSettings.clear().mergeFrom(fitSolverSettings);
    updateFitSolverSettings();
  }

  private void updateFitSolverSettings() {
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
  public FilterSettings getFilterSettings() {
    return filterSettings.build();
  }

  /**
   * Merge the FilterSettings.
   *
   * @param filterSettings the new FilterSettings
   */
  public void mergeFilterSettings(FilterSettings filterSettings) {
    this.filterSettings.clear().mergeFrom(filterSettings);
    updateFilterSettings();
  }

  /**
   * Sets the filter settings.
   *
   * @param filterSettings the new filter settings
   */
  public void setFilterSettings(FilterSettings filterSettings) {
    fitSettings.setFilterSettings(filterSettings);
    this.filterSettings = fitSettings.getFilterSettingsBuilder();
    // this.filterSettings.mergeFrom(filterSettings);
    updateFilterSettings();
  }

  private void updateFilterSettings() {
    updateMinSignal();
    updatePrecisionThreshold();
    updateCoordinateShift();
    updateWidthThreshold();
    updateMinWidthThreshold();
    updateZFilter();

    // Recreate the smart filter
    if (!filterSettings.getSmartFilter()) {
      return;
    }
    final String xml = filterSettings.getSmartFilterString();
    if (TextUtils.isNullOrEmpty(xml)) {
      return;
    }
    final Filter f = Filter.fromXml(xml);
    if (f == null || !(f instanceof DirectFilter)) {
      // Throw to ensure the filter is OK
      throw new IllegalStateException("Unrecognised smart filter: " + xml);
      // or
      // setDirectFilter(null);
    }

    // This updates the SmartFilter flag and the SmartFilterString.
    // Just set the filter directly
    // setDirectFilter((DirectFilter) f);
    this.directFilter = (DirectFilter) f;
  }

  @Override
  public FitConfiguration clone() {
    // Make a new initialised instance. This will have new settings builder objects.
    return new FitConfiguration(getFitSettings(), getCalibration(), getPsf()).copySettings(this);

    // // This is not a complete duplicate. The settings builder objects with the
    // // underlying configuration will be the same between all instances.
    // try
    // {
    // FitConfiguration f = (FitConfiguration) super.clone();
    // // Reset instance specific objects
    // f.toleranceChecker = null;
    // f.gaussianFunction = null;
    // f.functionSolver = null;
    // f.setValidationResult(null, null);
    // f.dynamicPeakResult = new DynamicPeakResult();
    // return f;
    // }
    // catch (CloneNotSupportedException e)
    // {
    // // Ignore
    // }
    // return null;
  }

  /**
   * Copy settings from the other configuration. This copies all the instance fields not stored in
   * settings objects that can be shared between instances. It is used in the {@link #clone()}
   * method after a new instance has been created with the current settings objects.
   *
   * @param other the other configuration
   * @return the fit configuration
   */
  FitConfiguration copySettings(FitConfiguration other) {
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
    shiftOffset = other.shiftOffset;
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
   * @param npeaks The number of peaks to fit
   * @param maxx The height of the XY data
   * @param maxy the maxy
   * @param params The Gaussian parameters
   */
  @Override
  public void initialise(int npeaks, int maxx, int maxy, double[] params) {
    {
      // XXX: For debugging thread safety require new objects for each fit
      // invalidateGaussianFunction();
      // invalidateToleranceChecker();
    }

    // Check if the Gaussian function is invalid
    if (gaussianFunction != null && (gaussianFunction.getNPeaks() != npeaks
        || gaussianFunction.getMaxX() != maxx || gaussianFunction.getMaxY() != maxy)) {
      // The gaussian function cannot be reused.
      // Do not call invalidate as it also invalidates the solver and we can re-use that.
      // invalidateGaussianFunction();
      gaussianFunction = null;
    }
    if (gaussianFunction == null) {
      // TODO : See if this works ...
      // We can update the function solver with the new function so do not invalidate the solver
      // invalidateFunctionSolver();
      gaussianFunction = createGaussianFunction(npeaks, maxx, maxy);
    }
    if (toleranceChecker == null) {
      // Requires a new function solver
      invalidateFunctionSolver();
    }
  }

  /**
   * Creates the appropriate 2D Gaussian function for the configuration.
   *
   * @param npeaks The number of peaks to fit
   * @param maxx The width of the XY data
   * @param maxy The height of the XY data
   * @return The function
   */
  public Gaussian2DFunction createGaussianFunction(int npeaks, int maxx, int maxy) {
    return GaussianFunctionFactory.create2D(npeaks, maxx, maxy, getFunctionFlags(),
        getAstigmatismZModel());
  }

  /**
   * Gets the function flags used for the GaussianFunctionFactory.
   *
   * @return the function flags
   */
  public int getFunctionFlags() {
    int flags = this.flags;
    if (!isBackgroundFitting()) {
      // Remove background fitting (on by default)
      flags &= ~GaussianFunctionFactory.FIT_BACKGROUND;
    }
    if (isNotSignalFitting()) {
      // Remove signal fitting (on by default)
      flags &= ~GaussianFunctionFactory.FIT_SIGNAL;
      flags |= GaussianFunctionFactory.FIT_SIMPLE;
    }
    return flags;
  }

  /**
   * Sets the log.
   *
   * <p>Used to output fit evaluations for each iteration.
   *
   * @param log the log to set
   */
  public void setLog(Logger log) {
    this.log = log;
  }

  /**
   * Gets the log.
   *
   * @return the log.
   */
  public Logger getLog() {
    return log;
  }

  /**
   * Sets the initial angle.
   *
   * @param initialAngle the new initial angle
   */
  public void setInitialAngle(double initialAngle) {
    // Provide backward compatibility
    if (psf.getPsfType() != PSFType.TWO_AXIS_AND_THETA_GAUSSIAN_2D) {
      throw new IllegalStateException("Not a 2 axis and theta Gaussian 2D PSF");
    }
    psf.getParametersBuilder(PsfHelper.INDEX_THETA).setValue(initialAngle);
  }

  /**
   * Gets the initial angle.
   *
   * @return the initialAngle.
   */
  @Override
  public double getInitialAngle() {
    // Provide backward compatibility
    if (psf.getPsfType() != PSFType.TWO_AXIS_AND_THETA_GAUSSIAN_2D) {
      throw new IllegalStateException("Not a 2 axis and theta Gaussian 2D PSF");
    }
    return psf.getParameters(PsfHelper.INDEX_THETA).getValue();
  }

  /**
   * Gets the PSF type.
   *
   * @return the PSF type
   */
  public PSFType getPsfType() {
    return psf.getPsfType();
  }

  /**
   * Gets the PSF type value.
   *
   * @return the PSF type value
   */
  public int getPsfTypeValue() {
    return psf.getPsfTypeValue();
  }

  /**
   * Sets the PSF type.
   *
   * <p>If the type is astigmatism and the astigmatism model cannot be constructed from the current
   * PSF parameters then the result filtering state may be incorrect. It is safer to call
   * {@link #setAstigmatismModel(uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.AstigmatismModel)}
   * which also updates the PSF type to astigmatism.
   *
   * @param psfType the new PSF type
   */
  public void setPsfType(PSFType psfType) {
    psf.setPsfType(psfType);
    updatePsf(true);
  }

  /**
   * Sets the initial peak standard deviation.
   *
   * @param initialPeakStdDev An estimate for the peak standard deviation used to initialise the fit
   *        for all dimensions
   */
  public void setInitialPeakStdDev(double initialPeakStdDev) {
    psf.getParametersBuilder(PsfHelper.INDEX_SX).setValue(initialPeakStdDev);
    if (isTwoAxisGaussian2D) {
      psf.getParametersBuilder(PsfHelper.INDEX_SY).setValue(initialPeakStdDev);
    }
    astigmatismZModel = null;
    updateCoordinateShift();
  }

  /**
   * Set an estimate for the peak standard deviation used to initialise the fit for dimension 0.
   *
   * <p>Setting this will update the value in {@link #getCoordinateShift()}.
   *
   * @param initialPeakStdDev0 An estimate for the peak standard deviation used to initialise the
   *        fit for dimension 0
   */
  public void setInitialPeakStdDev0(double initialPeakStdDev0) {
    psf.getParametersBuilder(PsfHelper.INDEX_SX).setValue(initialPeakStdDev0);
    astigmatismZModel = null;
    updateCoordinateShift();
  }

  /**
   * Gets the initial peak standard deviation.
   *
   * @return An estimate for the combined peak standard deviation
   * @throws ConfigurationException if the PSF type is astigmatism and the model cannot be
   *         constructed
   */
  public double getInitialPeakStdDev() {
    if (isTwoAxisGaussian2D) {
      return Gaussian2DPeakResultHelper.getStandardDeviation(getInitialXSd(), getInitialYSd());
    }
    return getInitialXSd();
  }

  /**
   * Gets the initial XSD.
   *
   * @return An estimate for the peak standard deviation used to initialise the fit for dimension 0
   * @throws ConfigurationException if the PSF type is astigmatism and the model cannot be
   *         constructed
   */
  @Override
  public double getInitialXSd() {
    if (getAstigmatismZModel() != null) {
      return astigmatismZModel.getSx(0);
    }
    return psf.getParameters(PsfHelper.INDEX_SX).getValue();
  }

  /**
   * Set an estimate for the peak standard deviation used to initialise the fit for dimension 1.
   *
   * <p>Setting this will update the value in {@link #getCoordinateShift()}.
   *
   * @param initialPeakStdDev1 An estimate for the peak standard deviation used to initialise the
   *        fit for dimension 1
   */
  public void setInitialPeakStdDev1(double initialPeakStdDev1) {
    if (!isTwoAxisGaussian2D) {
      throw new IllegalStateException("Not a 2 axis Gaussian 2D PSF");
    }
    psf.getParametersBuilder(PsfHelper.INDEX_SY).setValue(initialPeakStdDev1);
    astigmatismZModel = null;
    updateCoordinateShift();
  }

  /**
   * Gets the initial YSD.
   *
   * @return An estimate for the peak standard deviation used to initialise the fit for dimension 1
   * @throws ConfigurationException if the PSF type is astigmatism and the model cannot be
   *         constructed
   */
  @Override
  public double getInitialYSd() {
    if (getAstigmatismZModel() != null) {
      return astigmatismZModel.getSy(0);
    }
    if (isTwoAxisGaussian2D) {
      return psf.getParameters(PsfHelper.INDEX_SY).getValue();
    }
    return getInitialXSd();
  }

  /**
   * Sets to true to compute the deviations.
   *
   * @param computeDeviations True if computing the parameter deviations
   */
  public void setComputeDeviations(boolean computeDeviations) {
    this.computeDeviations = computeDeviations;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Note: This is also true if validation is active and the precision method requires
   * computation of the deviations (see {@link #isFilterRequiresDeviations()}).
   */
  @Override
  public boolean isComputeDeviations() {
    if (computeDeviations) {
      return true;
    }
    return isFilterRequiresDeviations();
  }

  /**
   * Get the value of the compute deviations flag. This may be false but
   * {@link #isComputeDeviations()} can still return true.
   *
   * @return the value of the compute deviations flag
   */
  public boolean getComputeDeviationsFlag() {
    return computeDeviations;
  }

  /**
   * Checks if the current filter settings require deviations.
   *
   * @return true, if filtering requires deviations
   */
  public boolean isFilterRequiresDeviations() {
    if (isDirectFilter() && directFilter.requiresParameterDeviations()) {
      return true;
    }
    return (precisionThreshold > 0
        && getPrecisionMethodValue() == PrecisionMethod.POISSON_CRLB_VALUE);
  }

  /**
   * Gets the fit solver used to fit the point spread function (PSF).
   *
   * @return the fit solver
   */
  public FitSolver getFitSolver() {
    return fitSolverSettings.getFitSolver();
  }

  /**
   * Gets the fit solver enum value.
   *
   * @return the fit solver enum value
   */
  public int getFitSolverValue() {
    return fitSolverSettings.getFitSolverValue();
  }

  /**
   * Sets the fit solver to use to fit the point spread function (PSF).
   *
   * @param fitSolver the fit solver
   */
  public void setFitSolver(FitSolver fitSolver) {
    invalidateFunctionSolver();
    updatePrecisionThreshold();
    fitSolverSettings.setFitSolver(fitSolver);
  }

  /**
   * Sets the fit solver to use to fit the point spread function (PSF).
   *
   * @param fitSolver the fit solver
   */
  public void setFitSolver(int fitSolver) {
    final FitSolver f = FitSolver.forNumber(fitSolver);
    if (f != null) {
      setFitSolver(f);
    }
  }

  /**
   * Sets the fixed iterations flag.
   *
   * @param fixedIterations the fixedIterations to set
   */
  public void setFixedIterations(boolean fixedIterations) {
    invalidateToleranceChecker();
    fitSolverSettings.setFixedIterations(fixedIterations);
  }

  /**
   * Checks if using a fixed number of iterations.
   *
   * @return true, using a fixed number of iterations
   */
  public boolean isFixedIterations() {
    return fitSolverSettings.getFixedIterations();
  }

  /**
   * Sets the max iterations.
   *
   * @param maxIterations the new max iterations
   */
  public void setMaxIterations(int maxIterations) {
    invalidateToleranceChecker();
    fitSolverSettings.setMaxIterations(Math.max(0, maxIterations));
  }

  /**
   * Gets the max iterations.
   *
   * @return the maxIterations.
   */
  public int getMaxIterations() {
    return Math.max(0, fitSolverSettings.getMaxIterations());
  }

  /**
   * Sets the background fitting option.
   *
   * @param backgroundFitting True if fitting the background
   */
  public void setBackgroundFitting(boolean backgroundFitting) {
    invalidateGaussianFunction();
    fitSolverSettings.setDisableBackgroundFitting(!backgroundFitting);
  }

  /**
   * Checks if fitting the background.
   *
   * @return True if fitting the background.
   */
  @Override
  public boolean isBackgroundFitting() {
    return !fitSolverSettings.getDisableBackgroundFitting();
  }

  /**
   * Use this to turn off fitting of the signal. This should be used with caution. The setting only
   * applies to fixed width fitting and can be used to benchmark position accuracy when fitting
   * signals of known strength.
   *
   * @param noSignalFitting True if not fitting the signal
   */
  public void setNotSignalFitting(boolean noSignalFitting) {
    invalidateGaussianFunction();
    fitSolverSettings.setDisableSignalFitting(noSignalFitting);
  }

  /**
   * Checks if not fitting the signal.
   *
   * <p>The setting only applies to fixed width fitting.
   *
   * @return True if not fitting the signal
   */
  public boolean isNotSignalFitting() {
    return fitSolverSettings.getDisableBackgroundFitting();
  }

  /**
   * Checks if fitting an elliptical peak (with an angle parameter).
   *
   * @return True if fitting an elliptical peak (with an angle parameter).
   */
  @Override
  public boolean isAngleFitting() {
    return (flags & GaussianFunctionFactory.FIT_ANGLE) != 0;
  }

  /**
   * Checks if fitting the z-position.
   *
   * @return True if fitting the z-position.
   */
  @Override
  public boolean isZFitting() {
    return (flags & GaussianFunctionFactory.FIT_Z) != 0;
  }

  /**
   * Checks if fitting the peak width in dimension 0.
   *
   * @return True if fitting the peak width in dimension 0.
   */
  @Override
  public boolean isXSdFitting() {
    return (flags & GaussianFunctionFactory.FIT_X_WIDTH) != 0;
  }

  /**
   * Checks if fitting the peak width in dimension 1.
   *
   * @return True if fitting the peak width in dimension 1.
   */
  @Override
  public boolean isYSdFitting() {
    return (flags & GaussianFunctionFactory.FIT_Y_WIDTH) != 0;
  }

  /**
   * Set to true to fix the PSF using the initial parameters. This is only supported for a one-axis
   * Gaussian 2D PSF.
   *
   * @param fixed the new fixed PSF
   */
  public void setFixedPsf(boolean fixed) {
    fitSolverSettings.setFixedPsf(fixed);
    updatePsf(true);
  }

  /**
   * Set to true to fix the PSF using the initial parameters. This is only supported for a one-axis
   * Gaussian 2D PSF.
   *
   * @return the fixed PSF
   */
  public boolean isFixedPsf() {
    return fitSolverSettings.getFixedPsf();
  }

  @Override
  public boolean isFitValidation() {
    return isDirectFilter() || isRegionValidation() || !isDisableSimpleFilter();
  }

  /**
   * Set the maximum absolute coordinate shift for a good fit. This is also set when calling
   * {@link #setCoordinateShiftFactor(double)} or any of the standard deviations, e.g.
   * {@link #setInitialPeakStdDev(double)}. If these are set then the coordinate shift will change.
   *
   * @param coordinateShift The maximum absolute coordinate shift for a good fit
   */
  public void setCoordinateShift(double coordinateShift) {
    this.coordinateShift = coordinateShift;
  }

  /**
   * Gets the maximum absolute coordinate shift for a good fit.
   *
   * @return The maximum absolute coordinate shift for a good fit.
   */
  public double getCoordinateShift() {
    return coordinateShift;
  }

  /**
   * Set the maximum absolute coordinate shift for a good fit, relative to the largest peak width.
   * Set to zero to disable.
   *
   * <p>Setting this will update the value in {@link #getCoordinateShift()}
   *
   * @param shiftFactor The maximum absolute coordinate shift for a good fit, relative to the
   *        largest peak width
   */
  public void setCoordinateShiftFactor(double shiftFactor) {
    filterSettings.setShiftFactor(shiftFactor);
    updateCoordinateShift();
  }

  private void updateCoordinateShift() {
    final double shiftFactor = getCoordinateShiftFactor();
    if (shiftFactor > 0) {
      // It may throw if a model cannot be created for the astigmatism type.
      double widthMax;
      try {
        widthMax = getWidthMax();
      } catch (final ConfigurationException ex) {
        if (getPsfTypeValue() != PSFType.ASTIGMATIC_GAUSSIAN_2D_VALUE) {
          throw ex;
        }
        // This is OK as the full astigmatism model may not have been set yet.
        // Just use a dummy value of 1
        widthMax = 1;
      }
      if (widthMax > 0) {
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
   * @throws ConfigurationException if the PSF type is astigmatism and the model cannot be
   *         constructed
   */
  public double getWidthMax() {
    double widthMax = getInitialXSd();
    if (isTwoAxisGaussian2D) {
      widthMax = Math.max(widthMax, getInitialYSd());
    }
    return widthMax;
  }

  /**
   * Gets the coordinate shift relative the the largest peak width.
   *
   * @return the coordinate shift relative the the largest peak width.
   */
  public double getCoordinateShiftFactor() {
    return filterSettings.getShiftFactor();
  }

  /**
   * Set the size of the fit region. Any coordinate outside the region will fail fit validation (see
   * {@link #validatePeak(int, double[], double[], double[])}). Set to zero to disable.
   *
   * <p>Note: it is assumed that the coordinates of the peak are relative to the fit region of size
   * NxN. Coordinates are offset by the amount defined by the coordinate offset, e.g. 0.5 pixels.
   *
   * @param fitRegionWidth the fit region width
   * @param fitRegionHeight the fit region height
   * @param coordinateOffset the coordinate offset when validating the coordinates are within the
   *        fit window
   */
  public void setFitRegion(int fitRegionWidth, int fitRegionHeight, double coordinateOffset) {
    this.fitRegionWidth = Math.max(0, fitRegionWidth);
    this.fitRegionHeight = Math.max(0, fitRegionHeight);
    this.coordinateOffset = coordinateOffset;
  }

  /**
   * Gets the width of the fit region used for validation.
   *
   * @return the width of the fit region used for validation.
   */
  public int getFitRegionWidth() {
    return fitRegionWidth;
  }

  private boolean isRegionValidation() {
    return fitRegionWidth != 0;
  }

  /**
   * Gets the height of the fit region used for validation.
   *
   * @return the height of the fit region used for validation.
   */
  public int getFitRegionHeight() {
    return fitRegionHeight;
  }

  /**
   * Gets the coordinate offset when validating the coordinates are within the fit window.
   *
   * @return the coordinate offset when validating the coordinates are within the fit window.
   */
  public double getCoordinateOffset() {
    return coordinateOffset;
  }

  /**
   * Set the signal strength (Signal-to-Noise Ratio, SNR) for a good fit.
   *
   * @param signalStrength The signal strength
   */
  public void setSignalStrength(double signalStrength) {
    filterSettings.setSignalStrength(signalStrength);
  }

  /**
   * Gets the signal strength.
   *
   * @return the signal strength.
   */
  public double getSignalStrength() {
    return filterSettings.getSignalStrength();
  }

  /**
   * Gets the minimum number of photons.
   *
   * @return The minimum number of photons.
   */
  public double getMinPhotons() {
    return filterSettings.getMinPhotons();
  }

  /**
   * Set the minimum photons used to determine the signal strength for a good fit (signalThreshold =
   * max(minSignal, noise x signalStrength).
   *
   * <p>Note that minSignal is created appropriately from minPhotons using the type of fitter, see
   * {@link #isFitCameraCounts()}.
   *
   * @param minPhotons The minimum number of photons
   */
  public void setMinPhotons(double minPhotons) {
    filterSettings.setMinPhotons(minPhotons);
    updateMinSignal();
  }

  private void updateMinSignal() {
    minSignal = (isFitCameraCounts()) ? getMinPhotons() * gain : getMinPhotons();
  }

  /**
   * Gets the precision threshold used to determine if the peak is a good fit. Requires that the
   * image is calibrated
   *
   * @return the precision threshold
   */
  public double getPrecisionThreshold() {
    return filterSettings.getPrecisionThreshold();
  }

  /**
   * Sets the precision threshold.
   *
   * @param precisionThreshold the precisionThreshold to set
   */
  public void setPrecisionThreshold(double precisionThreshold) {
    if (precisionThreshold > 0) {
      filterSettings.setPrecisionThreshold(precisionThreshold);
    } else {
      filterSettings.clearPrecisionThreshold();
    }
    updatePrecisionThreshold();
  }

  private void updatePrecisionThreshold() {
    // Note: Store the squared threshold for speed.
    precisionThreshold = 0;
    switch (getPrecisionMethodValue()) {
      case PrecisionMethod.MORTENSEN_VALUE:
      case PrecisionMethod.MORTENSEN_LOCAL_BACKGROUND_VALUE:
        // XXX - Determine if precision filtering for SCMOS is valid.
        // For now we leave this in but it may have to be changed to have a precision
        // computed during the fit which is stored for validation.
        if (nmPerPixel > 0 && gain > 0 && (calibration.isCcdCamera() || calibration.isScmos())) {
          this.precisionThreshold = MathUtils.pow2(getPrecisionThreshold());
        }
        break;
      case PrecisionMethod.POISSON_CRLB_VALUE:
        this.precisionThreshold = MathUtils.pow2(getPrecisionThreshold());
        break;
      default:
        break;
    }
    isMle = isMle();
  }

  /**
   * Checks if calculating the precision using the fitted background.
   *
   * @return True if calculating the precision using the fitted background.
   */
  public boolean isPrecisionUsingBackground() {
    return getPrecisionMethodValue() == PrecisionMethod.MORTENSEN_LOCAL_BACKGROUND_VALUE;
  }

  /**
   * Gets the precision method used to calculate the precision.
   *
   * @return the precision method
   */
  public PrecisionMethod getPrecisionMethod() {
    return filterSettings.getPrecisionMethod();
  }

  /**
   * Gets the precision method used to calculate the precision.
   *
   * @return the precision method
   */
  public int getPrecisionMethodValue() {
    return filterSettings.getPrecisionMethodValue();
  }

  /**
   * Sets the precision method used to calculate the precision.
   *
   * @param precisionMethod the new precision method
   */
  public void setPrecisionMethod(int precisionMethod) {
    final PrecisionMethod pm = PrecisionMethod.forNumber(precisionMethod);
    if (pm != null) {
      setPrecisionMethod(pm);
    }
  }

  /**
   * Sets the precision method used to calculate the precision.
   *
   * @param precisionMethod the new precision method
   */
  public void setPrecisionMethod(PrecisionMethod precisionMethod) {
    filterSettings.setPrecisionMethodValue(precisionMethod.getNumber());
    updatePrecisionThreshold();
  }

  // Enable caching the location variance selection
  private class BaseVarianceSelector {
    /**
     * Gets the location variance.
     *
     * @param peak the peak
     * @return the location variance
     */
    double getLocationVariance(PreprocessedPeakResult peak) {
      return 0;
    }
  }

  private class VarianceSelector extends BaseVarianceSelector {
    @Override
    double getLocationVariance(PreprocessedPeakResult peak) {
      return peak.getLocationVariance();
    }
  }

  private class VarianceSelector2 extends BaseVarianceSelector {
    @Override
    double getLocationVariance(PreprocessedPeakResult peak) {
      return peak.getLocationVariance2();
    }
  }

  private class VarianceSelectorCrlb extends BaseVarianceSelector {
    @Override
    double getLocationVariance(PreprocessedPeakResult peak) {
      return peak.getLocationVarianceCrlb();
    }
  }

  /**
   * Gets the precision method that will be used to produce the precision value for filtering.
   *
   * <p>This checks first the direct filter and then the current value for the precision method.
   *
   * <p>The value is written into the current calibration to allow the calibration to be used in
   * saved results.
   *
   * @return the filter precision method
   */
  public PrecisionMethod getFilterPrecisionMethod() {
    PrecisionMethod method = getFilterPrecisionMethodInternal();
    switch (method) {
      case MORTENSEN:
        varianceSelector = new VarianceSelector();
        break;
      case MORTENSEN_LOCAL_BACKGROUND:
        varianceSelector = new VarianceSelector2();
        break;
      case POISSON_CRLB:
        varianceSelector = new VarianceSelectorCrlb();
        break;
      default:
        method = PrecisionMethod.PRECISION_METHOD_NA;
        varianceSelector = new BaseVarianceSelector();
        break;
    }
    calibration.setPrecisionMethod(method);
    return method;
  }

  private PrecisionMethod getFilterPrecisionMethodInternal() {
    if (isDirectFilter()) {
      final int flags = directFilter.getValidationFlags();
      if (DirectFilter.areSet(flags, IDirectFilter.V_LOCATION_VARIANCE_CRLB)) {
        return PrecisionMethod.POISSON_CRLB;
      }
      if (DirectFilter.areSet(flags, IDirectFilter.V_LOCATION_VARIANCE2)) {
        return PrecisionMethod.MORTENSEN_LOCAL_BACKGROUND;
      }
      if (DirectFilter.areSet(flags, IDirectFilter.V_LOCATION_VARIANCE)) {
        return PrecisionMethod.MORTENSEN;
      }
    }
    if (precisionThreshold == 0) {
      // No simple precision filter
      return PrecisionMethod.PRECISION_METHOD_NA;
    }
    switch (getPrecisionMethodValue()) {
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
   * Gets the location variance of the result that was used for result filtering. This will use the
   * method matching that from the last call to {@link #getFilterPrecisionMethod()}. This allows the
   * state to be cached for efficient selection of the location variance.
   *
   * @param peak the peak
   * @return the variance
   */
  public double getLocationVariance(PreprocessedPeakResult peak) {
    return varianceSelector.getLocationVariance(peak);
  }

  /**
   * Set the image noise used to determine the Signal-to-Noise Ratio (SNR) for a good fit.
   *
   * @param noise The image noise.
   */
  public void setNoise(double noise) {
    this.noise = noise;
  }

  /**
   * Gets the image noise.
   *
   * @return the image noise.
   */
  public double getNoise() {
    return noise;
  }

  /**
   * Sets the maximum factor difference allowed between widths for a good fit.
   *
   * @param maxWidthFactor The maximum factor difference allowed between widths for a good fit
   */
  public void setMaxWidthFactor(double maxWidthFactor) {
    filterSettings.setMaxWidthFactor(maxWidthFactor);
    updateWidthThreshold();
  }

  private void updateWidthThreshold() {
    final double w = filterSettings.getMaxWidthFactor();
    if (w > 1) {
      this.widthFactor = (isTwoAxisGaussian2D) ? w * w : w;
    } else {
      this.widthFactor = Double.POSITIVE_INFINITY;
    }
  }

  /**
   * Gets the maximum factor difference allowed between widths for a good fit.
   *
   * <p>Returns zero if not configured.
   *
   * @return the maximum factor difference allowed between widths for a good fit
   */
  @Override
  public double getMaxWidthFactor() {
    final double widthFactor = filterSettings.getMaxWidthFactor();
    return (widthFactor > 1) ? widthFactor : 0;
  }

  /**
   * Sets the minimum factor difference allowed between widths for a good fit.
   *
   * @param minWidthFactor The minimum factor difference allowed between widths for a good fit
   */
  public void setMinWidthFactor(double minWidthFactor) {
    filterSettings.setMinWidthFactor(minWidthFactor);
    updateMinWidthThreshold();
  }

  private void updateMinWidthThreshold() {
    final double w = filterSettings.getMinWidthFactor();
    if (w < 1 && w > 0) {
      this.minWidthFactor = (isTwoAxisGaussian2D) ? w * w : w;
    } else {
      this.minWidthFactor = 0;
    }
  }

  /**
   * Gets the minimum factor difference allowed between widths for a good fit.
   *
   * <p>Returns zero if not configured.
   *
   * @return the minimum factor difference allowed between widths for a good fit
   */
  @Override
  public double getMinWidthFactor() {
    final double minWidthFactor = filterSettings.getMinWidthFactor();
    return (minWidthFactor < 1 && minWidthFactor > 0) ? minWidthFactor : 0;
  }

  /**
   * Sets the min Z. If both min and max z are zero then depth filtering is disabled.
   *
   * @param minZ The minimum z depth
   */
  public void setMinZ(double minZ) {
    filterSettings.setMinZ(minZ);
    updateZFilter();
  }

  /**
   * Gets the min Z.
   *
   * @return the min Z
   */
  public double getMinZ() {
    return filterSettings.getMinZ();
  }

  /**
   * Sets the max Z. If both min and max z are zero then depth filtering is disabled.
   *
   * @param maxZ The maximum z depth
   */
  public void setMaxZ(double maxZ) {
    filterSettings.setMaxZ(maxZ);
    updateZFilter();
  }

  /**
   * Gets the max Z.
   *
   * @return the max Z
   */
  public double getMaxZ() {
    return filterSettings.getMaxZ();
  }

  private void updateZFilter() {
    // Check if the PSF is 3D but avoid exceptions from creating the z-model
    zEnabled = getPsfTypeValue() == PSFType.ASTIGMATIC_GAUSSIAN_2D_VALUE
        // is3D() // This can throw
        && (getMaxZ() != 0 || getMinZ() != 0) && (getMinZ() <= getMaxZ());
  }

  /**
   * Sets the lambda.
   *
   * @param lambda the lambda to start the Levenberg-Marquardt fitting process
   */
  public void setLambda(double lambda) {
    invalidateFunctionSolver();
    fitSolverSettings.setLambda(lambda);
  }

  /**
   * Gets the lambda.
   *
   * @return the lambda
   */
  public double getLambda() {
    return fitSolverSettings.getLambda();
  }

  /**
   * Checks if computing the residuals.
   *
   * @return the computeResiduals
   */
  @Override
  public boolean isComputeResiduals() {
    return computeResiduals;
  }

  /**
   * Set to true to compute the residuals.
   *
   * @param computeResiduals Set to true to compute the residuals
   */
  public void setComputeResiduals(boolean computeResiduals) {
    this.computeResiduals = computeResiduals;
  }

  /**
   * Sets the nm per pixel scale to use when evaluating a fitted peak's localisation precision.
   *
   * @param nmPerPixel the nm per pixel scale to use when evaluating a fitted peak's localisation
   *        precision
   */
  public void setNmPerPixel(double nmPerPixel) {
    calibration.setNmPerPixel(nmPerPixel);
    updateCalibration();
  }

  /**
   * Gets the gain (or 1 if the gain is invalid).
   *
   * @return the gain (or 1 if the gain is invalid)
   */
  public double getGainSafe() {
    return (gain <= 0) ? 1 : gain;
  }

  /**
   * Sets the camera gain to use when evaluating a fitted peak's localisation precision.
   *
   * @param gain the camera gain to use when evaluating a fitted peak's localisation precision.
   */
  public void setGain(double gain) {
    invalidateFunctionSolver();
    invalidateCameraModel();
    calibration.setCountPerPhoton(gain);
    updateCalibration();
    // updateSignalThreshold();
  }

  /**
   * Specify the camera type used.
   *
   * <p>Specifying a CCD camera is relevant when validating results using the localisation
   * precision.
   *
   * @param cameraType the new camera type
   */
  public void setCameraType(CameraType cameraType) {
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
  public CameraType getCameraType() {
    return calibration.getCameraType();
  }

  /**
   * Gets the camera type.
   *
   * @return the camera type
   */
  public int getCameraTypeValue() {
    return calibration.getCameraTypeValue();
  }

  /**
   * Checks if modelling the camera noise during maximum likelihood fitting.
   *
   * @return True if modelling the camera noise during maximum likelihood fitting
   */
  public boolean isModelCamera() {
    return fitSolverSettings.getModelCamera();
  }

  /**
   * Specify if the camera noise should be modelled during maximum likelihood fitting. If true then
   * the read noise must be set. If the EmCCD property is true then the gain must also be set.
   *
   * @param modelCamera Set to true to model the camera
   */
  public void setModelCamera(boolean modelCamera) {
    invalidateFunctionSolver();
    fitSolverSettings.setModelCamera(modelCamera);
  }

  /**
   * Sets the bias (used for maximum likelihood estimation to evaluate the correct value of the
   * observed count).
   *
   * @param bias the camera bias
   */
  public void setBias(double bias) {
    invalidateCameraModel();
    calibration.setBias(bias);
  }

  /**
   * Sets the camera read noise (used for maximum likelihood estimation).
   *
   * @param readNoise the camera read noise (used for maximum likelihood estimation)
   */
  public void setReadNoise(double readNoise) {
    invalidateFunctionSolver();
    invalidateCameraModel();
    calibration.setReadNoise(readNoise);
  }

  /**
   * Sets the quantum efficiency.
   *
   * @param quantumEfficiency the new quantum efficiency [electron/photon] (used for maximum
   *        likelihood estimation)
   */
  public void setQuantumEfficiency(double quantumEfficiency) {
    invalidateFunctionSolver();
    calibration.setQuantumEfficiency(quantumEfficiency);
  }

  /**
   * Gets maximum number of function evaluations for the Maximum Likelihood Estimator.
   *
   * @return the maximum number of function evaluations for the Maximum Likelihood Estimator
   */
  public int getMaxFunctionEvaluations() {
    return fitSolverSettings.getMaxFunctionEvaluations();
  }

  /**
   * Sets the maximum number of function evaluations for the Maximum Likelihood Estimator.
   *
   * @param maxFunctionEvaluations the maximum number of function evaluations for the Maximum
   *        Likelihood Estimator
   */
  public void setMaxFunctionEvaluations(int maxFunctionEvaluations) {
    invalidateFunctionSolver();
    fitSolverSettings.setMaxFunctionEvaluations(maxFunctionEvaluations);
  }

  /**
   * Gets the search method for the Maximum Likelihood Estimator.
   *
   * @return the search for the Maximum Likelihood Estimator
   */
  public SearchMethod getSearchMethod() {
    return fitSolverSettings.getSearchMethod();
  }

  /**
   * Gets the search method value for the Maximum Likelihood Estimator.
   *
   * @return the search for the Maximum Likelihood Estimator
   */
  public int getSearchMethodValue() {
    return fitSolverSettings.getSearchMethodValue();
  }

  /**
   * Sets the search method for the Maximum Likelihood Estimator.
   *
   * @param searchMethod the search for the Maximum Likelihood Estimator
   */
  public void setSearchMethod(int searchMethod) {
    final SearchMethod sm = SearchMethod.forNumber(searchMethod);
    if (sm != null) {
      setSearchMethod(sm);
    }
  }

  /**
   * Sets the search method for the Maximum Likelihood Estimator.
   *
   * @param searchMethod the search method for the Maximum Likelihood Estimator
   */
  public void setSearchMethod(SearchMethod searchMethod) {
    invalidateFunctionSolver();
    fitSolverSettings.setSearchMethodValue(searchMethod.getNumber());
  }

  /**
   * Gets the line search method for the Fast MLE.
   *
   * @return the line search for the Fast MLE
   */
  public LineSearchMethod getLineSearchMethod() {
    return fitSolverSettings.getLineSearchMethod();
  }

  /**
   * Gets the line search method value for the Fast MLE.
   *
   * @return the line search for the Fast MLE
   */
  public int getLineSearchMethodValue() {
    return fitSolverSettings.getLineSearchMethodValue();
  }

  /**
   * Sets the line search method for the Fast MLE.
   *
   * @param lineSearchMethod the line search for the Fast MLE
   */
  public void setLineSearchMethod(int lineSearchMethod) {
    final LineSearchMethod sm = LineSearchMethod.forNumber(lineSearchMethod);
    if (sm != null) {
      setLineSearchMethod(sm);
    }
  }

  /**
   * Sets the line search method for the Fast MLE.
   *
   * @param lineSearchMethod the line search for the Fast MLE
   */
  public void setLineSearchMethod(LineSearchMethod lineSearchMethod) {
    invalidateFunctionSolver();
    fitSolverSettings.setLineSearchMethodValue(lineSearchMethod.getNumber());
  }

  /**
   * This setting applies to the conjugate gradient method of the Maximum Likelihood Estimator.
   *
   * @return the gradientLineMinimisation True if using the gradient for line minimisation
   */
  public boolean isGradientLineMinimisation() {
    return fitSolverSettings.getGradientLineMinimisation();
  }

  /**
   * This setting applies to the conjugate gradient method of the Maximum Likelihood Estimator.
   *
   * @param gradientLineMinimisation Set to true to use the gradient for line minimisation
   */
  public void setGradientLineMinimisation(boolean gradientLineMinimisation) {
    invalidateFunctionSolver();
    fitSolverSettings.setGradientLineMinimisation(gradientLineMinimisation);
  }

  /**
   * Gets the relative threshold for convergence.
   *
   * @return the relative threshold for convergence
   */
  public double getRelativeThreshold() {
    return fitSolverSettings.getRelativeThreshold();
  }

  /**
   * Sets the relative threshold for convergence.
   *
   * @param relativeThreshold the relative threshold for convergence
   */
  public void setRelativeThreshold(double relativeThreshold) {
    invalidateToleranceChecker();
    fitSolverSettings.setRelativeThreshold(relativeThreshold);
  }

  /**
   * Gets the absolute threshold for convergence.
   *
   * @return the absolute threshold for convergence
   */
  public double getAbsoluteThreshold() {
    return fitSolverSettings.getAbsoluteThreshold();
  }

  /**
   * Sets the absolute threshold for convergence.
   *
   * @param absoluteThreshold the absolute threshold for convergence
   */
  public void setAbsoluteThreshold(double absoluteThreshold) {
    invalidateToleranceChecker();
    fitSolverSettings.setAbsoluteThreshold(absoluteThreshold);
  }

  /**
   * Gets the parameter relative threshold for convergence.
   *
   * @return the parameter relative threshold for convergence
   */
  public double getParameterRelativeThreshold() {
    return fitSolverSettings.getParameterRelativeThreshold();
  }

  /**
   * Sets the parameter relative threshold for convergence.
   *
   * @param relativeThreshold the parameter relative threshold for convergence
   */
  public void setParameterRelativeThreshold(double relativeThreshold) {
    invalidateToleranceChecker();
    fitSolverSettings.setParameterRelativeThreshold(relativeThreshold);
  }

  /**
   * Gets the parameter absolute threshold for convergence.
   *
   * @return the parameter absolute threshold for convergence
   */
  public double getParameterAbsoluteThreshold() {
    return fitSolverSettings.getParameterAbsoluteThreshold();
  }

  /**
   * Sets the parameter absolute threshold for convergence.
   *
   * @param absoluteThreshold the parameter absolute threshold for convergence
   */
  public void setParameterAbsoluteThreshold(double absoluteThreshold) {
    invalidateToleranceChecker();
    fitSolverSettings.setParameterAbsoluteThreshold(absoluteThreshold);
  }

  /**
   * Gets the tolerance checker.
   *
   * @return the toleranceChecker
   */
  public ToleranceChecker getToleranceChecker() {
    if (toleranceChecker == null) {
      // This can be set to negative for fixed iterations
      int maxIterations = getMaxIterations();
      if (isFixedIterations()) {
        maxIterations = -maxIterations;
      }

      final boolean minimiseValue = isMinimiseValue();
      // If the value is zero then ignore this for convergence.
      toleranceChecker =
          new ToleranceChecker(minimiseValue, getIfStrictlyPositive(getRelativeThreshold()),
              getIfStrictlyPositive(getAbsoluteThreshold()),
              getIfStrictlyPositive(getParameterRelativeThreshold()),
              getIfStrictlyPositive(getParameterAbsoluteThreshold()), maxIterations);
    }
    return toleranceChecker;
  }

  private boolean isMinimiseValue() {
    switch (getFitSolverValue()) {
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
   * @param value the value
   * @return the value if positive, or -1
   */
  private static double getIfStrictlyPositive(double value) {
    return (value > 0) ? value : -1.0;
  }

  /**
   * Call this when a property changes that will change the stopping criteria.
   */
  private void invalidateToleranceChecker() {
    toleranceChecker = null;
    invalidateFunctionSolver();
  }

  /**
   * Gets the gaussian function.
   *
   * @return the gaussianFunction.
   */
  public Gaussian2DFunction getGaussianFunction() {
    return gaussianFunction;
  }

  /**
   * Call this when a property changes that will change the Gaussian function.
   */
  private void invalidateGaussianFunction() {
    gaussianFunction = null;
    invalidateFunctionSolver();
  }

  /**
   * Interface used to pass additional data for validation of fitting results.
   */
  public interface PeakResultValidationData {
    /**
     * Sets the result.
     *
     * @param peakNumber The peak number
     * @param initialParams The initial peak parameters
     * @param params The fitted peak parameters
     * @param paramDevs the fitted peak parameter variances (can be null)
     */
    void setResult(int peakNumber, double[] initialParams, double[] params, double[] paramDevs);

    /**
     * Gets the local background.
     *
     * @return the local background
     */
    double getLocalBackground();

    /**
     * Gets the noise.
     *
     * @return the noise
     */
    double getNoise();
  }

  private FitStatus result;
  private Object statusData;
  private PeakResultValidationData peakResultValidationData;

  /**
   * Gets the validation result.
   *
   * @return the result.
   */
  public FitStatus getValidationResult() {
    return result;
  }

  /**
   * Gets the validation data associated with the validation result.
   *
   * @return Data associated with the validation result.
   */
  @Override
  public Object getValidationData() {
    return statusData;
  }

  private FitStatus setValidationResult(FitStatus newResult, Object data) {
    result = newResult;
    statusData = data;
    return result;
  }

  /**
   * Sets the peak result validation data. This is used to obtain extra information about each peak
   * during calls to {@link #validatePeak(int, double[], double[], double[])}.
   *
   * <p>The object is discarded after a call to
   * {@link #validateFit(int, double[], double[], double[])} or
   * {@link #validateFit(double[], double[], double[])}.
   *
   * @param peakResultValidationData the new peak result validation data
   */
  public void setPeakResultValidationData(PeakResultValidationData peakResultValidationData) {
    this.peakResultValidationData = peakResultValidationData;
  }

  @Override
  public FitStatus validateFit(int peaks, double[] initialParams, double[] params,
      double[] paramDevs) {
    for (int n = 0; n < peaks; n++) {
      validatePeak(n, initialParams, params, paramDevs);
      if (result != FitStatus.OK) {
        break;
      }
    }
    peakResultValidationData = null;
    return result;
  }

  /**
   * Check peak to see if the fit was sensible. Assumes a single peak.
   *
   * @param initialParams The initial peak parameters
   * @param params The fitted peak parameters
   * @param paramDevs the fitted peak parameter variances (can be null)
   * @return True if the fit fails the criteria
   */
  public FitStatus validateFit(double[] initialParams, double[] params, double[] paramDevs) {
    validatePeak(0, initialParams, params, paramDevs);
    peakResultValidationData = null;
    return result;
  }

  /**
   * Check peak to see if the fit was sensible.
   *
   * @param n The peak number
   * @param initialParams The initial peak parameters
   * @param params The fitted peak parameters
   * @param paramDevs the fitted peak parameter variances (can be null)
   * @return True if the fit fails the criteria
   */
  public FitStatus validatePeak(int n, double[] initialParams, double[] params,
      double[] paramDevs) {
    // This requires local background and noise so that validation works the same
    // way for simple filteroing as it would for a PreprocessedPeakResult.

    // TODO - Update this so it can use a callback function to get
    // local background / noise.
    // These can be cached by the FitWorker so it doesn't have to compute them again.

    if (isDirectFilter()) {
      // Always specify a new result and we have no local background or offset.
      // Use the global noise.
      if (peakResultValidationData != null) {
        peakResultValidationData.setResult(n, initialParams, params, paramDevs);
      }
      final PreprocessedPeakResult peak = createPreprocessedPeakResult(0, n, initialParams, params,
          paramDevs, peakResultValidationData, ResultType.NEW, 0, 0, false);
      if (directFilter.accept(peak)) {
        return setValidationResult(FitStatus.OK, null);
      }
      if (log != null) {
        log.info(() -> String.format("Bad peak %d: %s", peak.getId(),
            DirectFilter.getStatusMessage(peak, directFilter.getResult())));
      }
      if (DirectFilter.anySet(directFilter.getResult(), V_X_SD_FACTOR | V_Y_SD_FACTOR)) {
        return setValidationResult(FitStatus.WIDTH_DIVERGED, null);
      }
      // At the moment we do not get any other validation data
      return setValidationResult(FitStatus.FAILED_SMART_FILTER, null);
    }

    // Check if outside the fit window.
    // TODO - Make this configurable per peak. At the moment we only use this in BenchmarkSpotFit
    // where
    // additional peaks will be neighbours. In the future we may want to control this better.
    if (isRegionValidation()) {
      final int offset = n * Gaussian2DFunction.PARAMETERS_PER_PEAK;
      final double x = params[Gaussian2DFunction.X_POSITION + offset] + coordinateOffset;
      final double y = params[Gaussian2DFunction.Y_POSITION + offset] + coordinateOffset;
      if (x <= 0 || x >= fitRegionWidth || y <= 0 || y >= fitRegionHeight) {
        if (log != null) {
          log.info(() -> String.format(
              "Bad peak %d: Coordinates outside fit region (x=%g,y=%g) <> %d,%d", n, x, y,
              fitRegionWidth, fitRegionHeight));
        }
        return setValidationResult(FitStatus.OUTSIDE_FIT_REGION,
            new double[] {x, y, fitRegionWidth, fitRegionHeight});
      }
    }

    if (isDisableSimpleFilter()) {
      return setValidationResult(FitStatus.OK, null);
    }

    final int offset = n * Gaussian2DFunction.PARAMETERS_PER_PEAK;
    // Check spot movement
    final double xShift = params[Gaussian2DFunction.X_POSITION + offset]
        - initialParams[Gaussian2DFunction.X_POSITION + offset];
    final double yShift = params[Gaussian2DFunction.Y_POSITION + offset]
        - initialParams[Gaussian2DFunction.Y_POSITION + offset];
    final double maxShift = coordinateShift;
    if (Math.abs(xShift) > maxShift || Math.abs(yShift) > maxShift) {
      if (log != null) {
        log.info(() -> String.format("Bad peak %d: Fitted coordinates moved (x=%g,y=%g) > %g", n,
            xShift, yShift, maxShift));
      }
      return setValidationResult(FitStatus.COORDINATES_MOVED, new double[] {xShift, yShift});
    }

    if (zEnabled) {
      final double z = params[Gaussian2DFunction.Z_POSITION + offset];
      if (z < getMinZ() || z > getMaxZ()) {
        return setValidationResult(FitStatus.Z_MOVED, z);
      }
    }

    // Check signal threshold.
    // The threshold should be set in the same units as those used during fitting.
    final double signal = params[Gaussian2DFunction.SIGNAL + offset];
    // Compare the signal to the desired signal strength
    if (signal < minSignal) {
      if (log != null) {
        log.info(() -> String.format("Bad peak %d: Insufficient signal %g", n, signal));
      }
      // if (params.length == 7) // Single peak
      // System.out.printf("Bad peak %d: Insufficient signal (%g)\n", n, signal);
      return setValidationResult(FitStatus.INSUFFICIENT_SIGNAL, signal);
    }

    double xsd;
    double ysd;
    // Map the width parameters using the z-model if present
    if (getAstigmatismZModel() != null) {
      final double z = params[Gaussian2DFunction.Z_POSITION + offset];
      xsd = astigmatismZModel.getSx(z);
      ysd = astigmatismZModel.getSy(z);

      // Check widths. This may be the only filter used even if z-fitting
      // (i.e. a z-depth filter is not used)
      final double xFactor = xsd / initialParams[Gaussian2DFunction.X_SD + offset];
      final double yFactor = ysd / initialParams[Gaussian2DFunction.Y_SD + offset];
      final double s2 = xFactor * yFactor;

      if (s2 > widthFactor || s2 < minWidthFactor) {
        if (log != null) {
          log.info(() -> String.format("Bad peak %d: Fitted width diverged (x=%gx,y=%gx)", n,
              xFactor, yFactor));
        }
        return setValidationResult(FitStatus.WIDTH_DIVERGED, new double[] {xFactor, yFactor});
      }
    } else {
      xsd = params[Gaussian2DFunction.X_SD + offset];
      ysd = params[Gaussian2DFunction.Y_SD + offset];
    }

    double noise = this.noise;
    if (peakResultValidationData != null) {
      peakResultValidationData.setResult(n, initialParams, params, paramDevs);
      noise = peakResultValidationData.getNoise();
    }

    // Check SNR threshold. Assume the noise has been set in the same units as the fit result.
    final double snr = Gaussian2DPeakResultHelper.getMeanSignalUsingP05(signal, xsd, ysd) / noise;

    // Compare the signal to the desired signal strength
    if (snr < getSignalStrength()) {
      if (log != null) {
        log.info(() -> String.format("Bad peak %d: Insufficient SNR %g", n, snr));
      }
      // if (params.length == 7) // Single peak
      // System.out.printf("Bad peak %d: Insufficient SNR (%g)\n", n, snr);
      return setValidationResult(FitStatus.INSUFFICIENT_SNR, snr);
    }

    // Check widths
    if (isXSdFitting()) {
      boolean badWidth = false;
      double xfactor = 0;
      double yfactor = 0;

      xfactor = xsd / initialParams[Gaussian2DFunction.X_SD + offset];
      if (isTwoAxisGaussian2D) {
        yfactor = ysd / initialParams[Gaussian2DFunction.Y_SD + offset];
        final double s2 = xfactor * yfactor;
        badWidth = (s2 > widthFactor || s2 < minWidthFactor);
      } else {
        badWidth = (xfactor > widthFactor || xfactor < minWidthFactor);
        yfactor = xfactor;
      }

      if (badWidth) {
        if (log != null) {
          final double localXFactor = xfactor;
          final double localYFactor = yfactor;
          log.info(() -> String.format("Bad peak %d: Fitted width diverged (x=%gx,y=%gx)", n,
              localXFactor, localYFactor));
        }
        return setValidationResult(FitStatus.WIDTH_DIVERGED, new double[] {xfactor, yfactor});
      }
    }

    // Check precision. This is above zero if a threshold is present.
    if (precisionThreshold > 0) {
      final double variance;
      switch (getPrecisionMethodValue()) {
        case PrecisionMethod.MORTENSEN_VALUE:
        case PrecisionMethod.MORTENSEN_LOCAL_BACKGROUND_VALUE:
          final double sd =
              (isTwoAxisGaussian2D) ? Gaussian2DPeakResultHelper.getStandardDeviation(xsd, ysd)
                  : xsd;
          final double localBackground =
              isPrecisionUsingBackground() && peakResultValidationData != null
                  ? peakResultValidationData.getLocalBackground()
                  : params[Gaussian2DFunction.BACKGROUND];
          variance = getVariance(localBackground,
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

      if (variance > precisionThreshold) {
        if (log != null) {
          final double precision = Math.sqrt(variance);
          log.info(
              () -> String.format("Bad peak %d: Insufficient precision (%gx)", n, precision));
        }
        return setValidationResult(FitStatus.INSUFFICIENT_PRECISION, variance);
      }
    }

    return setValidationResult(FitStatus.OK, null);
  }

  /**
   * Get the localisation variance for fitting a spot with the specified parameters given the
   * configuration (fit solver, precision using background, gain, nm per pixel).
   *
   * <p>We can calculate the precision using the estimated noise for the image or using the expected
   * number of background photons at the location.
   *
   * @param localBackground The background (in photons)
   * @param signal The signal (in photons)
   * @param sd The spot standard deviation
   * @param precisionUsingBackground calculate the precision using expected number of background
   *        photons at the location
   * @return The localisation variance
   */
  public double getVariance(double localBackground, final double signal, final double sd,
      boolean precisionUsingBackground) {
    double variance = 0;
    if (precisionUsingBackground) {
      // Check using the formula which uses the estimated background.
      // This allows for better filtering when the background is variable, e.g. when imaging cells.
      if (isMle) {
        try {
          // This may be slow due to the integration required within the formula.
          variance = Gaussian2DPeakResultHelper.getMLVarianceX(nmPerPixel, nmPerPixel * sd, signal,
              Math.max(0, localBackground), emCcd);
        } catch (final Exception ex) {
          // Catch all exceptions. They are likely to be a TooManyIterationsException and other
          // problems with the integration
          variance = Gaussian2DPeakResultHelper.getVarianceX(nmPerPixel, nmPerPixel * sd, signal,
              Math.max(0, localBackground), emCcd);
        }
      } else {
        variance = Gaussian2DPeakResultHelper.getVarianceX(nmPerPixel, nmPerPixel * sd, signal,
            Math.max(0, localBackground), emCcd);
      }
    } else if (isMle) {
      try {
        // This may be slow due to the integration required within the formula.
        variance = Gaussian2DPeakResultHelper.getMLVariance(nmPerPixel, nmPerPixel * sd, signal,
            noise, emCcd);
      } catch (final Exception ex) {
        // Catch all exceptions. They are likely to be a TooManyIterationsException and other
        // problems with the integration
        variance = Gaussian2DPeakResultHelper.getVariance(nmPerPixel, nmPerPixel * sd, signal,
            noise, emCcd);
      }
    } else {
      variance =
          Gaussian2DPeakResultHelper.getVariance(nmPerPixel, nmPerPixel * sd, signal, noise, emCcd);
    }
    return variance;
  }

  /**
   * Gets the variance. This is computed using the mean of the variance for the X and Y parameters.
   *
   * @param paramsDev the parameter variances
   * @param n the peak number
   * @return the variance (or zero if there are no deviations)
   */
  public double getVariance(double[] paramsDev, int n) {
    if (paramsDev != null) {
      final int offset = n * Gaussian2DFunction.PARAMETERS_PER_PEAK;
      // Scale to nm
      return nmPerPixel * nmPerPixel * (paramsDev[offset + Gaussian2DFunction.X_POSITION]
          + paramsDev[offset + Gaussian2DFunction.Y_POSITION]) / 2.0;
    }
    return 0;
  }

  /**
   * Update the estimated parameter variance. This method should be run on the deviations produced
   * by fitting.
   *
   * <p>If not an EM-CCD then the method does nothing.
   *
   * <p>If an EM-CCD the method scales the variance for all parameters by a factor of 2. This allows
   * the computation of the localisation variance using the parameter deviations to match the
   * estimation formulas of Mortensen. See Mortensen, et al (2010) Nature Methods 7, 377-383, SI 4.3
   * for assumptions and proof using MLE.
   *
   * @param paramsDev the params dev
   */
  public void updateVariance(double[] paramsDev) {
    if (emCcd && paramsDev != null) {
      // Scale all
      for (int i = 0; i < paramsDev.length; i++) {
        paramsDev[i] *= 2;
      }
    }
  }

  /**
   * An object that can return the results in a formatted state for the multi-path filter.
   */
  private class DynamicPeakResult implements PreprocessedPeakResult {
    int id;
    int candidateId;
    int offset;
    double[] initialParams;
    double[] params;
    double[] paramsDev;
    double xsd;
    double ysd;
    PeakResultValidationData peakResultValidationData;
    boolean existingResult;
    boolean newResult;
    float offsetx;
    float offsety;
    double var;
    double var2;
    double varCrlb;

    DynamicPeakResult(int candidateId, int n, double[] initialParams, double[] params,
        double[] paramsDev, PeakResultValidationData peakResultValidationData,
        ResultType resultType, float offsetx, float offsety) {
      setParameters(candidateId, n, initialParams, params, paramsDev, peakResultValidationData,
          resultType, offsetx, offsety);
    }

    DynamicPeakResult() {
      var = var2 = varCrlb = -1;
    }

    /**
     * Sets the parameters.
     *
     * @param candidateId the candidate id
     * @param n The peak number
     * @param initialParams The initial peak parameters
     * @param params The fitted peak parameters
     * @param paramsDev the parameter variances (can be null)
     * @param peakResultValidationData the peak result validation data
     * @param resultType the result type
     * @param offsetx the offsetx
     * @param offsety the offsety
     */
    void setParameters(int candidateId, int n, double[] initialParams, double[] params,
        double[] paramsDev, PeakResultValidationData peakResultValidationData,
        ResultType resultType, float offsetx, float offsety) {
      this.id = n;
      this.candidateId = candidateId;
      offset = n * Gaussian2DFunction.PARAMETERS_PER_PEAK;
      this.initialParams = initialParams;
      this.params = params;
      this.paramsDev = paramsDev;
      this.peakResultValidationData = peakResultValidationData;
      this.existingResult = resultType == ResultType.EXISTING;
      this.newResult = resultType == ResultType.NEW;
      this.offsetx = offsetx;
      this.offsety = offsety;
      var = var2 = varCrlb = -1;
      // Map the width parameters using the z-model
      if (getAstigmatismZModel() != null) {
        final double z = getZ();
        xsd = astigmatismZModel.getSx(z);
        ysd = astigmatismZModel.getSy(z);
      } else {
        xsd = params[Gaussian2DFunction.X_SD + offset];
        ysd = params[Gaussian2DFunction.Y_SD + offset];
      }
    }

    @Override
    public int getFrame() {
      // Not implemented
      return 0;
    }

    @Override
    public int getUniqueId() {
      // In the future we may want to use this so throw an exception so we notice
      throw new NotImplementedException("Unique Id not available");
    }

    @Override
    public int getId() {
      return id;
    }

    @Override
    public int getCandidateId() {
      return candidateId;
    }

    @Override
    public float getSignal() {
      return (float) (params[Gaussian2DFunction.SIGNAL + offset] * signalToPhotons);
    }

    @Override
    public float getMeanSignal() {
      return (float) Gaussian2DPeakResultHelper
          .getMeanSignalUsingP05(params[Gaussian2DFunction.SIGNAL + offset], xsd, ysd);
    }

    @Override
    public float getNoise() {
      if (peakResultValidationData != null) {
        // Because this is dynamic the id must be reset.
        // The other arguments should be identical to what is already stored.
        peakResultValidationData.setResult(id, initialParams, params, paramsDev);
        return (float) peakResultValidationData.getNoise();
      }

      return (float) FitConfiguration.this.noise;

      // This code is obsolete now that a noise estimate is used

      //// Note: This gets the noise estimate using the local background as photon
      //// shot noise. The noise is doubled if an EM-CCD.
      //// If the local background is zero then the global noise estimate is used.
      //
      //// TODO: Change so that the noise uses a region of the data around the spot centre,
      //// e.g. 5x5 pixels.
      //
      //// The local background is important for precision using the local background.
      //// it is not used for anything else.
      //
      //// Comment this out to use the configured local background, i.e. this.localBackground.
      //// If uncommented then the background will be either the local background or the fitted
      //// background.
      // final double localBackground = getLocalBackground();
      //
      // return (float) ((localBackground > 0)
      // ? (isFitCameraCounts()) ? PeakResultHelper.localBackgroundToNoise(localBackground, gain,
      //// emCCD)
      // : PeakResultHelper.localBackgroundToNoise(localBackground, emCCD)
      // : FitConfiguration.this.noise);
    }

    private double getLocalBackground() {
      if (peakResultValidationData != null) {
        // Because this is dynamic the id must be reset.
        // The other arguments should be identical to what is already stored.
        peakResultValidationData.setResult(id, initialParams, params, paramsDev);
        final double localBackground = peakResultValidationData.getLocalBackground();
        return (localBackground > 0) ? localBackground : params[Gaussian2DFunction.BACKGROUND];
      }
      return params[Gaussian2DFunction.BACKGROUND];
    }

    @Override
    public double getLocationVariance() {
      // We do not use the local background so set as zero
      if (var == -1) {
        var = FitConfiguration.this.getVariance(0,
            params[Gaussian2DFunction.SIGNAL + offset] * signalToPhotons, getSd(), false);
      }
      return var;
    }

    @Override
    public double getLocationVariance2() {
      if (var2 == -1) {
        var2 = FitConfiguration.this.getVariance(getLocalBackground() * signalToPhotons,
            params[Gaussian2DFunction.SIGNAL + offset] * signalToPhotons, getSd(), true);
      }
      return var2;
    }

    @Override
    public double getLocationVarianceCrlb() {
      if (varCrlb == -1) {
        varCrlb = FitConfiguration.this.getVariance(paramsDev, id);
      }
      return varCrlb;
    }

    @Override
    public float getSd() {
      if (isTwoAxisGaussian2D) {
        return (float) Gaussian2DPeakResultHelper.getStandardDeviation(xsd, ysd);
      }
      return (float) xsd;
    }

    @Override
    public float getBackground() {
      return (float) params[Gaussian2DFunction.BACKGROUND];
    }

    @Override
    public float getAmplitude() {
      return (float) (params[Gaussian2DFunction.SIGNAL + offset] / (2 * Math.PI * xsd * ysd));
    }

    @Override
    public float getAngle() {
      return (float) params[Gaussian2DFunction.ANGLE + offset];
    }

    @Override
    public float getX() {
      return (float) params[Gaussian2DFunction.X_POSITION + offset] + offsetx;
    }

    @Override
    public float getY() {
      return (float) params[Gaussian2DFunction.Y_POSITION + offset] + offsety;
    }

    @Override
    public float getZ() {
      return (float) params[Gaussian2DFunction.Z_POSITION + offset];
    }

    @Override
    public float getXRelativeShift2() {
      final double d = (params[Gaussian2DFunction.X_POSITION + offset]
          - initialParams[Gaussian2DFunction.X_POSITION + offset])
          / initialParams[Gaussian2DFunction.X_SD + offset];
      return (float) (d * d);
    }

    @Override
    public float getYRelativeShift2() {
      final double d = (params[Gaussian2DFunction.Y_POSITION + offset]
          - initialParams[Gaussian2DFunction.Y_POSITION + offset])
          / initialParams[Gaussian2DFunction.Y_SD + offset];
      return (float) (d * d);
    }

    @Override
    public float getXSd() {
      return (float) xsd;
    }

    @Override
    public float getYSd() {
      return (float) ysd;
    }

    @Override
    public float getXSdFactor() {
      return (float) (xsd / initialParams[Gaussian2DFunction.X_SD + offset]);
    }

    @Override
    public float getYSdFactor() {
      return (float) (ysd / initialParams[Gaussian2DFunction.Y_SD + offset]);
    }

    @Override
    public boolean isExistingResult() {
      return existingResult;
    }

    @Override
    public boolean isNewResult() {
      return newResult;
    }

    @Override
    public FractionalAssignment[] getAssignments(int predictedId) {
      return null;
    }

    @Override
    public double[] toGaussian2DParameters() {
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

    @Override
    public void setValidationResult(int result) {
      throw new NotImplementedException(
          "The validation result should not be set on a dynamic result");
    }

    @Override
    public int getValidationResult() {
      throw new NotImplementedException(
          "The validation result should not be set on a dynamic result");
    }

    @Override
    public boolean ignore() {
      return false;
    }

    @Override
    public boolean isNotDuplicate() {
      return false;
    }
  }

  /**
   * Create a dynamic object that can return the results in a formatted state for the multi-path
   * filter.
   *
   * <p>The result is dynamic in that it computes the values just-in-time using the input array
   * data.
   *
   * <p>The result can be a recycled object that is associated with this fit configuration, or a new
   * object. If using the recycled object then a second call to this method will replace the array
   * data on all references to the object. If using a new object then this method can be called
   * again with new data and the old reference is still valid.
   *
   * <p>Note: All returned objects will be linked with this fit configuration. Thus changing
   * properties such as the gain, noise or settings for computing the variance will result in
   * changes to the values returned by the PreprocessedPeakResult.
   *
   * <p>Note: XY position may be wrong if the input parameters have not been updated with an offset
   * from fitting a sub-region.
   *
   * @param candidateId the candidate id
   * @param n The peak number
   * @param initialParams The initial peak parameters
   * @param params The fitted peak parameters
   * @param paramVariances the parameter variances (can be null)
   * @param peakResultValidationData the peak result validation data
   * @param resultType the result type
   * @param offsetx the offsetx to adjust the x-position
   * @param offsety the offsety to adjust the y-position
   * @return A preprocessed peak result
   */
  public PreprocessedPeakResult createDynamicPreprocessedPeakResult(int candidateId, int n,
      double[] initialParams, double[] params, double[] paramVariances,
      PeakResultValidationData peakResultValidationData, ResultType resultType, float offsetx,
      float offsety) {
    return createPreprocessedPeakResult(candidateId, n, initialParams, params, paramVariances,
        peakResultValidationData, resultType, offsetx, offsety, true);
  }

  /**
   * Create a dynamic object that can return the results in a formatted state for the multi-path
   * filter.
   *
   * <p>The result is dynamic in that it computes the values just-in-time using the input array
   * data.
   *
   * <p>The result can be a recycled object that is associated with this fit configuration, or a new
   * object. If using the recycled object then a second call to this method will replace the array
   * data on all references to the object. If using a new object then this method can be called
   * again with new data and the old reference is still valid.
   *
   * <p>Note: All returned objects will be linked with this fit configuration. Thus changing
   * properties such as the gain, noise or settings for computing the variance will result in
   * changes to the values returned by the PreprocessedPeakResult.
   *
   * <p>Note: XY position may be wrong if the input parameters have not been updated with an offset
   * from fitting a sub-region.
   *
   * @param candidateId the candidate id
   * @param n The peak number
   * @param initialParams The initial peak parameters
   * @param params The fitted peak parameters
   * @param paramVariances the parameter variances (can be null)
   * @param peakResultValidationData the peak result validation data
   * @param resultType the result type
   * @param offsetx the offsetx to adjust the x-position
   * @param offsety the offsety to adjust the y-position
   * @param newObject Set to true to create a new object, the default uses the object associated
   *        with this fit configuration
   * @return A preprocessed peak result
   */
  private PreprocessedPeakResult createPreprocessedPeakResult(int candidateId, int n,
      double[] initialParams, double[] params, double[] paramVariances,
      PeakResultValidationData peakResultValidationData, ResultType resultType, float offsetx,
      float offsety, boolean newObject) {
    if (newObject) {
      return new DynamicPeakResult(candidateId, n, initialParams, params, paramVariances,
          peakResultValidationData, resultType, offsetx, offsety);
    }

    dynamicPeakResult.setParameters(candidateId, n, initialParams, params, paramVariances,
        peakResultValidationData, resultType, offsetx, offsety);
    return dynamicPeakResult;
  }

  /**
   * Create an object that can return the results in a formatted state for the multi-path filter.
   *
   * <p>The result is fixed in that it computes the values on construction using the input array
   * data.
   *
   * <p>If a local background is provided then it is used instead of the fitted background. The
   * local background can be computed if a multi-peak fit has been performed since the background
   * will be the global background, The local background for a peak will be the global background
   * plus the contribution of all the other peaks in the local region around the peak of interest.
   *
   * <p>The local background will be used to estimate the noise in the local region (as photon shot
   * noise) if it is above the bias.
   *
   * @param frame the frame
   * @param candidateId the candidate id
   * @param n The peak number
   * @param initialParameters the initial parameters
   * @param parameters the parameters
   * @param paramVariances the parameter variances (can be null)
   * @param peakResultValidationData the peak result validation data
   * @param resultType the result type
   * @param offsetx the offsetx
   * @param offsety the offsety
   * @return A preprocessed peak result
   */
  public BasePreprocessedPeakResult createPreprocessedPeakResult(int frame, int candidateId, int n,
      double[] initialParameters, double[] parameters, double[] paramVariances,
      PeakResultValidationData peakResultValidationData, ResultType resultType, float offsetx,
      float offsety) {
    peakResultValidationData.setResult(n, initialParameters, parameters, paramVariances);

    final int offset = n * Gaussian2DFunction.PARAMETERS_PER_PEAK;
    final double signal = parameters[offset + Gaussian2DFunction.SIGNAL] * signalToPhotons;
    final double b = signalToPhotons * peakResultValidationData.getLocalBackground();
    final double angle = parameters[offset + Gaussian2DFunction.ANGLE];
    final double x = parameters[offset + Gaussian2DFunction.X_POSITION] + offsetx;
    final double y = parameters[offset + Gaussian2DFunction.Y_POSITION] + offsety;
    final double z = parameters[offset + Gaussian2DFunction.Z_POSITION];
    final double x0 = initialParameters[offset + Gaussian2DFunction.X_POSITION] + offsetx;
    final double y0 = initialParameters[offset + Gaussian2DFunction.Y_POSITION] + offsety;
    double xsd;
    double ysd;
    // Map the width parameters using the z-model
    if (getAstigmatismZModel() != null) {
      xsd = astigmatismZModel.getSx(z);
      ysd = astigmatismZModel.getSy(z);
    } else {
      xsd = parameters[offset + Gaussian2DFunction.X_SD];
      ysd = parameters[offset + Gaussian2DFunction.Y_SD];
    }
    final double xsd0 = initialParameters[offset + Gaussian2DFunction.X_SD];
    final double ysd0 = initialParameters[offset + Gaussian2DFunction.Y_SD];
    final double sd =
        (isTwoAxisGaussian2D) ? Gaussian2DPeakResultHelper.getStandardDeviation(xsd, ysd) : xsd;
    final double variance = getVariance(0, signal, sd, false);
    final double variance2 = getVariance(b, signal, sd, true);
    final double varianceCrlb = getVariance(paramVariances, n);

    final double noise = peakResultValidationData.getNoise();
    final double u = Gaussian2DPeakResultHelper.getMeanSignalUsingP05(signal, xsd, ysd);
    return new BasePreprocessedPeakResult(frame, n, candidateId, signal, u, noise, b, angle, x, y,
        z, x0, y0, xsd, ysd, xsd0, ysd0, variance, variance2, varianceCrlb, resultType);
  }

  /**
   * Checks if fitting requires the camera counts.
   *
   * @return Set to true if fitting requires the camera counts, i.e. amplification is explicitly
   *         modelled during fitting.
   */
  public boolean isFitCameraCounts() {
    // Only the legacy MLE solvers that explicitly model the camera noise require the data
    // and estimate to be in ADUs. This is also true if there is no camera calibration.
    return (fitSolverSettings.getFitSolverValue() == FitSolver.MLE_VALUE
        || calibration.getCameraTypeValue() == CameraType.CAMERA_TYPE_NA_VALUE);
  }

  /**
   * Check if performing raw data fitting. In this mode then no camera model is available to map the
   * data counts to photons. The gain is assumed to be 1.
   *
   * @return Set to true if raw fitting
   */
  public boolean isRawFit() {
    // This is only true if there is no camera calibration.
    // Otherwise the camera details can be used to convert input data in camera counts
    // to photo-electrons and all fitting is done assuming the signal is photo-electrons.
    return (calibration.getCameraTypeValue() == CameraType.CAMERA_TYPE_NA_VALUE);
  }

  /**
   * Checks if is a full maximum likelihood estimator (MLE) modelling the CCD camera.
   *
   * @return true, if is a MLE modelling the camera
   */
  public boolean isModelCameraMle() {
    return (isModelCamera() && fitSolverSettings.getFitSolverValue() == FitSolver.MLE_VALUE);
  }

  /**
   * Checks if is using a maximum likelihood estimator (MLE).
   *
   * @return true, if is a MLE
   */
  public boolean isMle() {
    switch (getFitSolverValue()) {
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
  public boolean requireStrictlyPositiveFunction() {
    // Only the LSE variants can fit negatives. The MLE variants all require a positive function.
    switch (getFitSolverValue()) {
      case FitSolver.LVM_LSE_VALUE:
      case FitSolver.LVM_WLSE_VALUE:
        return false;

      case FitSolver.LVM_MLE_VALUE:
      case FitSolver.MLE_VALUE:
      case FitSolver.FAST_MLE_VALUE:
      case FitSolver.BACKTRACKING_FAST_MLE_VALUE:
        return true;

      default:
        throw new NotImplementedException(
            "Unknown strictly positive requirement: " + getFitSolver());
    }
  }

  /**
   * Gets the function solver for the current configuration.
   *
   * @return The function solver for the current configuration.
   */
  @Override
  public FunctionSolver getFunctionSolver() {
    if (functionSolver == null || gaussianFunction == null) {
      // The function solver was invalidated so create a new one
      functionSolver = createFunctionSolver();
    } else {
      // Update with the solver with the latest function.
      functionSolver.setGradientFunction(gaussianFunction);

      // Note: We must carefully update anything that depends on the function.
      if (bounds != null) {
        // We have to update the clamping.
        // Note this code is only executed if the clamp settings have not changed
        // (since changes to those settings invalidate the solver) and the
        // function settings have not changed (since that invalidates the function
        // and the solver).
        // All that is different is the number of peaks in the function.
        if (gaussianFunction.getNPeaks() > clampPeakCount) {
          setClampValues(bounds);
        }
      }
    }

    // Special case where the data and estimate are in counts but
    // the function solver requires the function to output an expected
    // number of photons.
    final boolean doScaling = getFitSolverValue() == FitSolver.MLE_VALUE;
    assert !doScaling || functionSolver instanceof MaximumLikelihoodFitter;

    if (precomputedFunctionValues != null) {
      if (doScaling) {
        // Scale the pre-computed function. It is assumes that this was
        // set using output fitted parameters (which would be unscaled).
        // This may be reused by the calling code so don't destroy it.
        final double[] f = new double[precomputedFunctionValues.length];
        for (int i = precomputedFunctionValues.length; i-- > 0;) {
          f[i] = precomputedFunctionValues[i] * signalToPhotons;
        }
        precomputedFunctionValues = f;
      }

      functionSolver.setGradientFunction((GradientFunction) OffsetFunctionFactory
          .wrapFunction(gaussianFunction, precomputedFunctionValues));
      precomputedFunctionValues = null;
    }
    if (functionSolver.isWeighted()) {
      functionSolver.setWeights(observationWeights);
    }

    if (doScaling) {
      // Scale the background and signal estimates
      final int[] indices = new int[1 + gaussianFunction.getNPeaks()];
      indices[0] = Gaussian2DFunction.BACKGROUND;
      for (int i = 1; i < indices.length; i++) {
        indices[i] = (i - 1) * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.SIGNAL;
      }

      // signalToPhotons should be 1.0/gain
      return new MleScaledFunctionSolver((MaximumLikelihoodFitter) functionSolver, signalToPhotons,
          indices);
    }

    return functionSolver;
  }

  /**
   * Call this when a property changes that will change the function solver.
   */
  private void invalidateFunctionSolver() {
    functionSolver = null;
    bounds = null;
  }

  private BaseFunctionSolver createFunctionSolver() {
    if (gaussianFunction == null) {
      // Other code may want to call getFunctionSolver() to see if exceptions are thrown
      // so create a dummy function so we can return a function solver.
      gaussianFunction = createGaussianFunction(1, 1, 1);
    }

    if (getFitSolverValue() == FitSolver.MLE_VALUE) {
      // Only support CCD/EM-CCD at the moment
      if (!calibration.isCcdCamera()) {
        throw new IllegalStateException(
            "CCD/EM-CCD camera is required for fit solver: " + getFitSolver());
      }

      // This requires the gain
      if (gain <= 0) {
        throw new IllegalStateException("The gain is required for fit solver: " + getFitSolver());
      }

      final MaximumLikelihoodFitter.SearchMethod searchMethod = convertSearchMethod();

      // Only the Poisson likelihood function supports gradients
      if (searchMethod.usesGradients() && isModelCamera()) {
        throw new IllegalStateException(String.format(
            "The derivative based search method '%s' can only be used with the "
                + "'%s' likelihood function, i.e. no model camera noise",
            searchMethod, MaximumLikelihoodFitter.LikelihoodFunction.POISSON));
      }

      final MaximumLikelihoodFitter fitter = new MaximumLikelihoodFitter(gaussianFunction);
      fitter.setRelativeThreshold(getRelativeThreshold());
      fitter.setAbsoluteThreshold(getAbsoluteThreshold());
      fitter.setMaxEvaluations(getMaxFunctionEvaluations());
      fitter.setMaxIterations(getMaxIterations());
      fitter.setSearchMethod(searchMethod);
      fitter.setGradientLineMinimisation(isGradientLineMinimisation());

      // Specify the likelihood function to use
      if (isModelCamera()) {
        // Set the camera read noise.
        // Do not check if this is set as 0 is a valid option.
        fitter.setSigma(calibration.getReadNoise());

        if (emCcd) {
          // EMCCD = Poisson+Gamma+Gaussian
          fitter.setLikelihoodFunction(
              MaximumLikelihoodFitter.LikelihoodFunction.POISSON_GAMMA_GAUSSIAN);
        } else {
          // CCD = Poisson+Gaussian
          fitter.setLikelihoodFunction(MaximumLikelihoodFitter.LikelihoodFunction.POISSON_GAUSSIAN);
        }
      } else {
        fitter.setLikelihoodFunction(MaximumLikelihoodFitter.LikelihoodFunction.POISSON);
      }

      // All models use the amplification gain (i.e. how many ADUs/electron)
      if (!calibration.hasCountPerElectron()) {
        throw new IllegalStateException(
            "The amplification is required for the fit solver: " + getFitSolver());
      }

      fitter.setAlpha(1.0 / calibration.getCountPerElectron());

      // TODO - Configure better stopping criteria ...

      return fitter;
    }

    // All the remaining solvers are based on the stepping function solver
    final ToleranceChecker tc = getToleranceChecker();
    final ParameterBounds bounds = new ParameterBounds(gaussianFunction);
    if (isUseClamping()) {
      setClampValues(bounds);
    }

    SteppingFunctionSolver solver;

    switch (getFitSolverValue()) {
      case FitSolver.LVM_LSE_VALUE:
        solver = new LseLvmSteppingFunctionSolver(gaussianFunction, tc, bounds);
        break;

      case FitSolver.LVM_MLE_VALUE:
        checkCameraCalibration();
        solver = new MleLvmSteppingFunctionSolver(gaussianFunction, tc, bounds);
        break;

      case FitSolver.LVM_WLSE_VALUE:
        checkCameraCalibration();
        solver = new WLseLvmSteppingFunctionSolver(gaussianFunction, tc, bounds);
        break;

      case FitSolver.FAST_MLE_VALUE:
        checkCameraCalibration();
        // This may throw a class cast exception if the function does not support
        // the Gradient2Function interface
        solver =
            new FastMleSteppingFunctionSolver((Gradient2Function) gaussianFunction, tc, bounds);
        break;

      case FitSolver.BACKTRACKING_FAST_MLE_VALUE:
        checkCameraCalibration();
        solver = new BacktrackingFastMleSteppingFunctionSolver((Gradient2Function) gaussianFunction,
            tc, bounds);
        break;

      default:
        throw new IllegalStateException("Unknown fit solver: " + getFitSolver());
    }

    if (solver instanceof LvmSteppingFunctionSolver) {
      ((LvmSteppingFunctionSolver) solver).setInitialLambda(getLambda());
    } else if (solver instanceof FastMleSteppingFunctionSolver) {
      ((FastMleSteppingFunctionSolver) solver).setLineSearchMethod(convertLineSearchMethod());
    }
    return solver;
  }

  private MaximumLikelihoodFitter.SearchMethod convertSearchMethod() {
    return FitProtosHelper.convertSearchMethod(getSearchMethod());
  }

  private FastMleSteppingFunctionSolver.LineSearchMethod convertLineSearchMethod() {
    return FitProtosHelper.convertLineSearchMethod(getLineSearchMethod());
  }

  private void checkCameraCalibration() {
    if (!calibration.hasCameraCalibration()) {
      throw new IllegalStateException(
          "The camera calibration is required for fit solver: " + getFitSolver());
    }

    switch (getCameraTypeValue()) {
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
        throw new IllegalStateException("Unrecognised camera type for for fit solver: "
            + getFitSolver() + ": " + getCameraType());
    }
  }

  private void setClampValues(ParameterBounds bounds) {
    final double[] clamp = getClampValues();
    // Note: The units are photons. This is OK as all solvers except the legacy MLE fit in photons.
    clampPeakCount = gaussianFunction.getNPeaks();
    final double[] clampValues =
        new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK * clampPeakCount];
    clampValues[Gaussian2DFunction.BACKGROUND] = clamp[Gaussian2DFunction.BACKGROUND];
    for (int i = 0; i < clampPeakCount; i++) {
      for (int j = 1; j <= Gaussian2DFunction.PARAMETERS_PER_PEAK; j++) {
        clampValues[i + j] = clamp[j];
      }
    }
    bounds.setClampValues(clampValues);
    bounds.setDynamicClamp(isUseDynamicClamping());
  }

  /**
   * Checks if simple filtering is disabled.
   *
   * @return True if simple filtering is disabled
   */
  public boolean isDisableSimpleFilter() {
    return filterSettings.getDisableSimpleFilter();
  }

  /**
   * Sets to true to disable simple filtering during validation.
   *
   * @param disableSimpleFilter Set to true to disable simple filtering during validation
   */
  public void setDisableSimpleFilter(boolean disableSimpleFilter) {
    filterSettings.setDisableSimpleFilter(disableSimpleFilter);
  }

  /**
   * Checks if filtering should use the configured smart filter.
   *
   * @return True if filtering should use the configured smart filter
   */
  public boolean isSmartFilter() {
    return filterSettings.getSmartFilter();
  }

  /**
   * Sets if filtering should use the configured smart filter.
   *
   * @param smartFilter True if filtering should use the configured smart filter
   */
  public void setSmartFilter(boolean smartFilter) {
    filterSettings.setSmartFilter(smartFilter);
  }

  /**
   * Checks if smart filter is enabled and a valid filter is present.
   *
   * @return true, if is direct filtering is enabled.
   */
  public boolean isDirectFilter() {
    return isSmartFilter() && directFilter != null;
  }

  /**
   * Gets the smart filter string.
   *
   * @return the smart filter string.
   */
  public String getSmartFilterString() {
    final String s = filterSettings.getSmartFilterString();
    return (TextUtils.isNullOrEmpty(s)) ? "" : s;
  }

  /**
   * This returns the representation of this object as a smart filter. This ignores any current
   * smart filter and only uses the standard filtering settings.
   *
   * @return the smart filter if using this object as a smart filter.
   */
  public DirectFilter getDefaultSmartFilter() {
    final double signal = getMinPhotons();
    final float snr = (float) getSignalStrength();
    final double minWidth = getMinWidthFactor();
    final double maxWidth = getMaxWidthFactor();
    final double shift = getCoordinateShiftFactor();
    final double eshift = 0;
    final double precision = getPrecisionThreshold();
    final float minZ = (float) getMinZ();
    final float maxZ = (float) getMaxZ();

    switch (getPrecisionMethodValue()) {
      case PrecisionMethod.MORTENSEN_VALUE:
        return new MultiFilter(signal, snr, minWidth, maxWidth, shift, eshift, precision, minZ,
            maxZ);
      case PrecisionMethod.MORTENSEN_LOCAL_BACKGROUND_VALUE:
        return new MultiFilter2(signal, snr, minWidth, maxWidth, shift, eshift, precision, minZ,
            maxZ);
      case PrecisionMethod.POISSON_CRLB_VALUE:
        return new MultiFilterCrlb(signal, snr, minWidth, maxWidth, shift, eshift, precision, minZ,
            maxZ);
      default:
        throw new IllegalStateException("Unknown precision method: " + getPrecisionMethod());
    }
  }

  /**
   * This returns the XML representation of this object as a smart fitler. This ignores any current
   * smart filter and only uses the standard filtering settings.
   *
   * @return the smart filter XML if using this object as a smart filter.
   */
  public String getDefaultSmartFilterXml() {
    return getDefaultSmartFilter().toXml();
  }

  /**
   * Sets the direct filter. This changes the smart filter flag to true and updates the smart filter
   * string.
   *
   * @param directFilter the new direct filter
   */
  public void setDirectFilter(DirectFilter directFilter) {
    this.directFilter = directFilter;
    if (directFilter != null) {
      setSmartFilter(true);
      filterSettings.setSmartFilterString(directFilter.toXml());
    } else {
      setSmartFilter(false);
      filterSettings.clearSmartFilterString();
    }
  }

  /**
   * Gets the smart filter name, if a smart filter exists.
   *
   * @return the smart filter name
   */
  public String getSmartFilterName() {
    if (directFilter != null) {
      return directFilter.getName();
    }
    return "";
  }

  /**
   * Gets the smart filter, if a smart filter exists. A clone is returned.
   *
   * @return the smart filter (or null)
   */
  public DirectFilter getSmartFilter() {
    if (directFilter != null) {
      return (DirectFilter) directFilter.clone();
    }
    return null;
  }

  @Override
  public void setup() {
    setup(0);
  }

  @Override
  public void setup(int flags) {
    filterSetupFlags = createFlags(flags);
    this.filterSetupData = null;

    if (directFilter != null) {
      directFilter.setup(flags);
    } else {
      widthEnabled = !DirectFilter.areSet(flags, IDirectFilter.NO_WIDTH);
      if (DirectFilter.areSet(flags, IDirectFilter.NO_SHIFT)) {
        shiftOffset = Float.POSITIVE_INFINITY;
      } else {
        final double shiftFactor = getCoordinateShiftFactor();
        shiftOffset =
            (float) ((shiftFactor > 0) ? shiftFactor * shiftFactor : Float.POSITIVE_INFINITY);
      }
      varianceThreshold = (precisionThreshold > 0) ? precisionThreshold : Double.POSITIVE_INFINITY;
    }
  }

  @Override
  public void setup(int flags, FilterSetupData... filterSetupData) {
    filterSetupFlags = createFlags(flags);
    this.filterSetupData = filterSetupData;

    if (directFilter != null) {
      directFilter.setup(flags, filterSetupData);
    } else {
      // Note: These variables should be copied in copySettings(...)
      widthEnabled = !DirectFilter.areSet(flags, IDirectFilter.NO_WIDTH);
      if (DirectFilter.areSet(flags, IDirectFilter.NO_SHIFT)) {
        shiftOffset = Float.POSITIVE_INFINITY;
      } else {
        double shiftFactor = getCoordinateShiftFactor();
        for (int i = filterSetupData.length; i-- > 0;) {
          if (filterSetupData[i] instanceof ShiftFilterSetupData) {
            final double shift = ((ShiftFilterSetupData) filterSetupData[i]).shift;
            if (shift > 0) {
              final double widthMax = getWidthMax();
              if (widthMax > 0) {
                shiftFactor = shift / widthMax;
              }
            }
            break;
          }
        }
        shiftOffset =
            (float) ((shiftFactor > 0) ? shiftFactor * shiftFactor : Float.POSITIVE_INFINITY);
      }
      varianceThreshold = (precisionThreshold > 0) ? precisionThreshold : Double.POSITIVE_INFINITY;
    }
  }

  private int createFlags(int flags) {
    // Handle switching to a 2 axis width filter
    if (isTwoAxisGaussian2D) {
      flags |= IDirectFilter.XY_WIDTH;
    }
    // Ignore z if not 3D
    if (!zEnabled) {
      flags |= IDirectFilter.NO_Z;
    }
    return flags;
  }

  @Override
  public int getFilterSetupFlags() {
    // Cached for speed
    return filterSetupFlags;
  }

  @Override
  public FilterSetupData[] getFilterSetupData() {
    // Cached for speed
    return filterSetupData;
  }

  @Override
  public boolean accept(PreprocessedPeakResult peak) {
    return (filterResult = validate(peak)) == 0;
  }

  @Override
  public int validate(PreprocessedPeakResult peak) {
    final int validateResult = doValidate(peak);
    if (log != null && validateResult != 0) {
      LoggerUtils.log(log, Level.INFO, "Bad peak %d (%.1f,%.1f) [%d]: %s", peak.getCandidateId(),
          peak.getX(), peak.getY(), peak.getId(),
          DirectFilter.getStatusMessage(peak, validateResult));
    }
    return validateResult;
  }

  @Override
  public int getValidationFlags() {
    if (directFilter != null) {
      return directFilter.getValidationFlags();
    }
    // Q. Is this necessary? See doValidate below ...
    if (isDisableSimpleFilter()) {
      return 0;
    }
    // These could be conditional on the filter settings. For now just set them all.
    int validationFlags = V_PHOTONS | V_SNR | V_X_RELATIVE_SHIFT | V_Y_RELATIVE_SHIFT;
    if (widthEnabled) {
      validationFlags |= V_X_SD_FACTOR;
      if (isTwoAxisGaussian2D) {
        validationFlags |= V_Y_SD_FACTOR;
      }
    }
    switch (getPrecisionMethodValue()) {
      case PrecisionMethod.MORTENSEN_VALUE:
        validationFlags |= V_LOCATION_VARIANCE;
        break;
      case PrecisionMethod.MORTENSEN_LOCAL_BACKGROUND_VALUE:
        validationFlags |= V_LOCATION_VARIANCE2;
        break;
      case PrecisionMethod.POISSON_CRLB_VALUE:
        validationFlags |= V_LOCATION_VARIANCE_CRLB;
        break;
      default:
        break;
    }
    return validationFlags;
  }

  /**
   * Do the validation. Either use the direct filter or perform the same function using the filter
   * configuration.
   *
   * @param peak the peak
   * @return the result
   * @see IDirectFilter#validate(PreprocessedPeakResult)
   */
  public int doValidate(PreprocessedPeakResult peak) {
    if (directFilter != null) {
      return directFilter.validate(peak);
    }

    // Q. Is this necessary? Simple filtering is to support turning off
    // filtering in the validatePeak(...) method. Set a debug point to check
    // if this is used.
    if (isDisableSimpleFilter()) {
      return 0;
    }

    // Do filtering
    if (peak.getSignal() < getMinPhotons()) {
      return V_PHOTONS;
    }
    if (peak.getSnr() < getSignalStrength()) {
      return V_SNR;
    }
    if (widthEnabled) {
      // Handle switching to a 2 axis width filter.
      // Note this will ignore any flags passed to the IDirectFilter#setup methods
      // to control the filter.
      if (isTwoAxisGaussian2D) {
        final float s2 = peak.getXSdFactor() * peak.getYSdFactor();
        if (s2 > widthFactor || s2 < minWidthFactor) {
          return V_X_SD_FACTOR | V_Y_SD_FACTOR;
        }
      } else {
        final double s = peak.getXSdFactor();
        if (s > widthFactor || s < minWidthFactor) {
          return V_X_SD_FACTOR;
        }
      }
    }
    if (peak.getXRelativeShift2() > shiftOffset) {
      return V_X_RELATIVE_SHIFT;
    }
    if (peak.getYRelativeShift2() > shiftOffset) {
      return V_Y_RELATIVE_SHIFT;
    }
    // Note: Do not support Euclidian shift

    if (zEnabled) {
      final double z = peak.getZ();
      if (z < getMinZ() || z > getMaxZ()) {
        return V_Z;
      }
    }

    // Note: The variance threshold is set to infinity if the variance filter is disabled
    switch (getPrecisionMethodValue()) {
      case PrecisionMethod.MORTENSEN_VALUE:
        if (peak.getLocationVariance() > varianceThreshold) {
          return V_LOCATION_VARIANCE;
        }
        break;
      case PrecisionMethod.MORTENSEN_LOCAL_BACKGROUND_VALUE:
        if (peak.getLocationVariance2() > varianceThreshold) {
          return V_LOCATION_VARIANCE2;
        }
        break;
      case PrecisionMethod.POISSON_CRLB_VALUE:
        if (peak.getLocationVarianceCrlb() > varianceThreshold) {
          return V_LOCATION_VARIANCE_CRLB;
        }
        break;
      default:
        throw new IllegalStateException("Unknown precision method: " + getPrecisionMethod());
    }

    return 0;
  }

  @Override
  public FilterType getFilterType() {
    return FilterType.DIRECT;
  }

  @Override
  public int getResult() {
    return filterResult;
  }

  @Override
  public IDirectFilter copy() {
    return clone();
  }

  /**
   * Checks if is use clamping the parameter update to a maximum value.
   *
   * @return Set to true to clamp the parameter update to a maximum value.
   */
  public boolean isUseClamping() {
    return fitSolverSettings.getUseClamping();
  }

  /**
   * Set to true to clamp the parameter update to a maximum value.
   *
   * @param useClamping Set to true to clamp the parameter update to a maximum value
   */
  public void setUseClamping(boolean useClamping) {
    invalidateFunctionSolver();
    fitSolverSettings.setUseClamping(useClamping);
  }

  /**
   * Checks if is using dynamic clamping. If true the clamp values are updated when the parameter
   * update direction changes.
   *
   * @return Set to true to update the clamp values when the parameter update direction changes.
   */
  public boolean isUseDynamicClamping() {
    return fitSolverSettings.getUseDynamicClamping();
  }

  /**
   * Set to true to update the clamp values when the parameter update direction changes.
   *
   * @param useDynamicClamping Set to true to update the clamp values when the parameter update
   *        direction changes
   */
  public void setUseDynamicClamping(boolean useDynamicClamping) {
    invalidateFunctionSolver();
    fitSolverSettings.setUseDynamicClamping(useDynamicClamping);
  }

  /**
   * Gets the clamp values.
   *
   * @return the clampValues.
   */
  private double[] getClampValues() {
    if (clampValues == null) {
      final int n = PeakResult.STANDARD_PARAMETERS + 3;
      if (fitSolverSettings.getClampValuesCount() != n) {
        throw new IllegalStateException("Require clamp values for all the Gaussian 2D parameters");
      }

      clampValues = new double[n];
      for (int i = 0; i < n; i++) {
        clampValues[i] = fitSolverSettings.getClampValues(i);
      }
    }
    return clampValues;
  }

  /**
   * Invalidate clamp values.
   */
  private void invalidateClampValues() {
    clampValues = null;
  }

  /**
   * Gets the clamp value for the background.
   *
   * @return The clamp value for the background.
   */
  public double getClampBackground() {
    return getClampValues()[Gaussian2DFunction.BACKGROUND];
  }

  /**
   * Sets the clamp value for the background.
   *
   * @param value the new clamp value
   */
  public void setClampBackground(double value) {
    updateClampValues(Gaussian2DFunction.BACKGROUND, value);
  }

  /**
   * Gets the clamp value for the signal.
   *
   * @return The clamp value for the signal.
   */
  public double getClampSignal() {
    return getClampValues()[Gaussian2DFunction.SIGNAL];
  }

  /**
   * Sets the clamp value for the signal.
   *
   * @param value the new clamp value
   */
  public void setClampSignal(double value) {
    updateClampValues(Gaussian2DFunction.SIGNAL, value);
  }

  /**
   * Gets the clamp value for the x position.
   *
   * @return The clamp value for the x position.
   */
  public double getClampX() {
    return getClampValues()[Gaussian2DFunction.X_POSITION];
  }

  /**
   * Sets the clamp value for the x position.
   *
   * @param value the new clamp value
   */
  public void setClampX(double value) {
    updateClampValues(Gaussian2DFunction.X_POSITION, value);
  }

  /**
   * Gets the clamp value for the y position.
   *
   * @return The clamp value for the y position.
   */
  public double getClampY() {
    return getClampValues()[Gaussian2DFunction.Y_POSITION];
  }

  /**
   * Sets the clamp value for the y position.
   *
   * @param value the new clamp value
   */
  public void setClampY(double value) {
    updateClampValues(Gaussian2DFunction.Y_POSITION, value);
  }

  /**
   * Gets the clamp value for the z position.
   *
   * @return The clamp value for the z position.
   */
  public double getClampZ() {
    return getClampValues()[Gaussian2DFunction.Z_POSITION];
  }

  /**
   * Sets the clamp value for the z position.
   *
   * @param value the new clamp value
   */
  public void setClampZ(double value) {
    updateClampValues(Gaussian2DFunction.Z_POSITION, value);
  }

  /**
   * Gets the value for the x standard deviation.
   *
   * @return The clamp value for the x sd.
   */
  public double getClampXSd() {
    return getClampValues()[Gaussian2DFunction.X_SD];
  }

  /**
   * Sets the clamp value for the x sd.
   *
   * @param value the new clamp value
   */
  public void setClampXSd(double value) {
    updateClampValues(Gaussian2DFunction.X_SD, value);
  }

  /**
   * Gets the value for the y standard deviation.
   *
   * @return The clamp value for the y sd.
   */
  public double getClampYSd() {
    return getClampValues()[Gaussian2DFunction.Y_SD];
  }

  /**
   * Sets the clamp value for the y sd.
   *
   * @param value the new clamp value
   */
  public void setClampYSd(double value) {
    updateClampValues(Gaussian2DFunction.Y_SD, value);
  }

  /**
   * Gets the value for the angle.
   *
   * @return The clamp value for the angle.
   */
  public double getClampAngle() {
    return getClampValues()[Gaussian2DFunction.ANGLE];
  }

  /**
   * Sets the clamp value for the angle.
   *
   * @param value the new clamp value
   */
  public void setClampAngle(double value) {
    updateClampValues(Gaussian2DFunction.ANGLE, value);
  }

  private void updateClampValues(int index, double value) {
    invalidateFunctionSolver();
    getClampValues()[index] = value;
  }

  /**
   * Gets the precomputed function values.
   *
   * @return the precomputed function values
   */
  public double[] getPrecomputedFunctionValues() {
    return precomputedFunctionValues;
  }

  /**
   * Sets the precomputed function values. This is combined with the configured Gaussian function
   * and passed to the function solver returned from {@link #getFunctionSolver()}. The precomputed
   * function values are then reset to null so it must be set each time the function solver is
   * accessed.
   *
   * @param precomputedFunctionValues the new precomputed function values
   */
  public void setPrecomputedFunctionValues(double[] precomputedFunctionValues) {
    this.precomputedFunctionValues = precomputedFunctionValues;
  }

  /**
   * Sets the observation weights. These are passed to the function solver if it supports weights.
   * These must be set each time the function solver from {@link #getFunctionSolver()} will be used
   * on new data.
   *
   * @param observationWeights the new observation weights
   */
  public void setObservationWeights(double[] observationWeights) {
    this.observationWeights = observationWeights;
  }

  /**
   * Gets the observation weights.
   *
   * @return the observation weights
   */
  public double[] getObservationWeights() {
    return observationWeights;
  }

  /**
   * Call this when a property changes that will change the camera model, e.g. bias, gain, read
   * noise, camera type.
   */
  private void invalidateCameraModel() {
    setCameraModel(null);
  }

  /**
   * Sets the camera model. This must be set if a sCMOS camera type is used.
   *
   * @param cameraModel the new camera model
   */
  public void setCameraModel(CameraModel cameraModel) {
    invalidateFunctionSolver();
    // Use this to set the bias and gain
    if (cameraModel != null && cameraModel.isPerPixelModel()) {
      calibration.clearGlobalCameraSettings();

      // Trigger an update to the calibration used for validation.
      updateCalibration();
    }
    this.cameraModel = cameraModel;
  }

  /**
   * Return true if the camera type requires a per-pixel camera model.
   *
   * @return true, if successful
   */
  public boolean isPerPixelCameraType() {
    switch (getCameraTypeValue()) {
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
   * @throws IllegalStateException if no camera model exists for the camera type
   */
  public CameraModel getCameraModel() {
    if (cameraModel == null) {
      final int value = getCameraTypeValue();
      switch (value) {
        case CameraType.CAMERA_TYPE_NA_VALUE:
          // We can support this by doing nothing to pixels values
          cameraModel = new NullCameraModel();
          break;

        case CameraType.CCD_VALUE:
        case CameraType.EMCCD_VALUE:
          final double bias = calibration.getBias();
          final double gain = calibration.getCountPerPhoton();
          final double variance = MathUtils.pow2(calibration.getReadNoise());
          // This will throw an exception if the calibration is invalid.
          cameraModel =
              (value == CameraType.EMCCD_VALUE) ? new EmCcdCameraModel(bias, gain, variance)
                  : new CcdCameraModel(bias, gain, variance);
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
   * Checks for a valid camera model type. This does not validate if a camera model is present, only
   * that an attempt can be made to create one.
   *
   * @return true, if successful
   */
  public boolean hasValidCameraModelType() {
    switch (getCameraTypeValue()) {
      case CameraType.CCD_VALUE:
      case CameraType.EMCCD_VALUE:
      case CameraType.SCMOS_VALUE:
        return true;

      default:
        return false;
    }
  }

  /**
   * Sets the camera model name. This should contain all the information required to load the camera
   * model, e.g. in the case of a per-pixel camera model for sCMOS cameras.
   *
   * <p>This settings is saved to the underlying configuration. If a camera model is used (e.g. for
   * sCMOS camera) then {@link #setCameraModel(CameraModel)} should be called after setting the new
   * camera model name.
   *
   * @param cameraModelName the new camera model name
   */
  public void setCameraModelName(String cameraModelName) {
    if (cameraModelName == null || !cameraModelName.equals(getCameraModelName())) {
      invalidateCameraModel();
    }
    calibration.setCameraModelName(cameraModelName);
  }

  /**
   * Gets the camera model name.
   *
   * @return the camera model name
   */
  public String getCameraModelName() {
    return calibration.getCameraModelName();
  }

  /**
   * Sets the PSF model name. This should contain all the information required to load the PSF
   * model, e.g. in the case of an astigmatic Gaussian 2D PSF.
   *
   * <p>This settings is saved to the underlying configuration. If a PSF model is used (e.g. for an
   * astigmatic Gaussian 2D PSF) then
   * {@link #setAstigmatismModel(uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.AstigmatismModel)}
   * should be called after setting the new PSF model name.
   *
   * @param psfModelName the new PSF model name
   */
  public void setPsfModelName(String psfModelName) {
    psf.setModelName(psfModelName);
  }

  /**
   * Gets the PSF model name.
   *
   * @return the PSF model name
   */
  public String getPsfModelName() {
    return psf.getModelName();
  }

  /**
   * Sets the astigmatism model. This has the effect of changing the PSF to an astigmatic Gaussian
   * 2D function.
   *
   * @param model the new astigmatism model
   * @throws ConfigurationException if the model pixel pitch does not match the calibration
   * @throws ConversionException if the model cannot be converted to pixel units
   */
  public void setAstigmatismModel(AstigmatismModel model) {
    // Check the calibration
    if (DoubleEquality.relativeError(model.getNmPerPixel(), calibration.getNmPerPixel()) > 1e-3) {
      throw new ConfigurationException(String.format(
          "The astigmatism model pixel pitch (%s) does not match the calibration (%s)",
          model.getNmPerPixel(), calibration.getNmPerPixel()));
    }

    // Convert to pixels
    model = PsfProtosHelper.convert(model, DistanceUnit.PIXEL, DistanceUnit.PIXEL);

    // Create the working model
    astigmatismZModel = HoltzerAstigmatismZModel.create(model.getS0X(), model.getS0Y(),
        model.getGamma(), model.getD(), model.getAx(), model.getBx(), model.getAy(), model.getBy());

    // Store the parameters in the PSF
    final String modelName = getPsfModelName();
    psf.clear().mergeFrom(PsfProtosHelper.createPsf(model, DistanceUnit.PIXEL, DistanceUnit.PIXEL));
    psf.setModelName(modelName);
    updatePsf(false);
  }

  /**
   * Gets the astigmatism Z model. This is only valid if the PSF type is an astigmatic Gaussian 2D
   * and the parameters are correctly configured.
   *
   * @return the astigmatism Z model (or null)
   * @throws ConfigurationException if the model cannot be created from the parameters
   */
  public AstigmatismZModel getAstigmatismZModel() {
    if (getPsfTypeValue() == PSFType.ASTIGMATIC_GAUSSIAN_2D_VALUE) {
      if (astigmatismZModel == null) {
        // Use the helper to convert the PSF parameters back to a model
        final AstigmatismModel model = PsfProtosHelper.createModel(getPsf(), DistanceUnit.PIXEL,
            DistanceUnit.PIXEL, calibration.getNmPerPixel());
        astigmatismZModel =
            HoltzerAstigmatismZModel.create(model.getS0X(), model.getS0Y(), model.getGamma(),
                model.getD(), model.getAx(), model.getBx(), model.getAy(), model.getBy());
      }
      return astigmatismZModel;
    }
    return null;
  }

  /**
   * Checks if the PSF is 3D. Currently only an Astigmatic 2D Gaussian is supported. A check is made
   * for a valid astigmatism z-model.
   *
   * @return true, if is 3D
   * @throws ConfigurationException if the 3D model cannot be created
   */
  public boolean is3D() {
    return (getPsfTypeValue() == PSFType.ASTIGMATIC_GAUSSIAN_2D_VALUE
        && getAstigmatismZModel() != null);
  }

  /**
   * Gets the Gaussian 2D x-width and y-width for the PSF parameters.
   *
   * @return the Gaussian 2D x-width and y-width for the PSF parameters.
   * @throws ConfigurationException if the psf is null, or not a Gaussian 2D function
   */
  double[] getGaussian2DWxWy() {
    return PsfHelper.getGaussian2DWxWy(psf);
  }
}
