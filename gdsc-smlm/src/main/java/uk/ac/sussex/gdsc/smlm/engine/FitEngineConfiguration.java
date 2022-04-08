/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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

import java.util.Objects;
import uk.ac.sussex.gdsc.core.data.NotImplementedException;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.DataFilter;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.DataFilterMethod;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.DataFilterSettings;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.DataFilterType;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FitEngineSettings;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.NoiseEstimatorMethod;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.RelativeParameter;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.data.config.PsfProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.filters.AverageDataProcessor;
import uk.ac.sussex.gdsc.smlm.filters.BlockAverageDataProcessor;
import uk.ac.sussex.gdsc.smlm.filters.CircularMeanDataProcessor;
import uk.ac.sussex.gdsc.smlm.filters.DataProcessor;
import uk.ac.sussex.gdsc.smlm.filters.DifferenceSpotFilter;
import uk.ac.sussex.gdsc.smlm.filters.GaussianDataProcessor;
import uk.ac.sussex.gdsc.smlm.filters.JurySpotFilter;
import uk.ac.sussex.gdsc.smlm.filters.MaximaSpotFilter;
import uk.ac.sussex.gdsc.smlm.filters.MedianDataProcessor;
import uk.ac.sussex.gdsc.smlm.filters.SingleSpotFilter;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.model.camera.CameraModel;
import uk.ac.sussex.gdsc.smlm.results.count.CombinedAndFailCounter;
import uk.ac.sussex.gdsc.smlm.results.count.ConsecutiveFailCounter;
import uk.ac.sussex.gdsc.smlm.results.count.FailCounter;
import uk.ac.sussex.gdsc.smlm.results.count.NullFailCounter;
import uk.ac.sussex.gdsc.smlm.results.count.PassRateFailCounter;

/**
 * Specifies the configuration for the fit engine.
 */
public final class FitEngineConfiguration {
  private FitEngineSettings.Builder fitEngineSettings;
  private FitConfiguration fitConfiguration;

  /**
   * Creates a new fit engine configuration.
   *
   * @return the fit engine configuration
   */
  public static FitEngineConfiguration create() {
    return create(FitProtosHelper.defaultFitEngineSettings,
        CalibrationProtosHelper.defaultCalibration, PsfProtosHelper.defaultOneAxisGaussian2DPSF);
  }

  /**
   * Creates a new fit engine configuration.
   *
   * @param fitEngineSettings the fit engine settings
   * @param calibration the calibration
   * @param psf the psf
   * @return the fit engine configuration
   */
  public static FitEngineConfiguration create(FitEngineSettings fitEngineSettings,
      Calibration calibration, PSF psf) {
    Objects.requireNonNull(fitEngineSettings, "fitEngineSettings");
    final FitConfiguration fitConfiguration =
        FitConfiguration.create(fitEngineSettings.getFitSettings(), calibration, psf);
    return new FitEngineConfiguration(fitEngineSettings.toBuilder(), fitConfiguration);
  }

  /**
   * Creates a new fit engine configuration.
   *
   * @param fitEngineSettings the fit engine settings
   * @param calibration the calibration
   * @param psf the psf
   * @return the fit engine configuration
   */
  public static FitEngineConfiguration create(FitEngineSettings.Builder fitEngineSettings,
      Calibration.Builder calibration, PSF.Builder psf) {
    Objects.requireNonNull(fitEngineSettings, "fitEngineSettings");
    final FitConfiguration fitConfiguration =
        FitConfiguration.create(fitEngineSettings.getFitSettingsBuilder(), calibration, psf);
    return new FitEngineConfiguration(fitEngineSettings, fitConfiguration);
  }

  /**
   * Creates a new fit engine configuration.
   *
   * @param fitEngineSettings the fit engine settings
   * @param fitConfiguration the fit configuration
   */
  private FitEngineConfiguration(FitEngineSettings.Builder fitEngineSettings,
      FitConfiguration fitConfiguration) {
    this.fitEngineSettings = fitEngineSettings;
    this.fitConfiguration = fitConfiguration;
  }

  /**
   * Gets the fit engine settings.
   *
   * @return the fit engine settings
   */
  public FitEngineSettings getFitEngineSettings() {
    return fitEngineSettings.build();
  }

  /**
   * Merge fit engine settings.
   *
   * @param fitEngineSettings the fit engine settings
   */
  public void mergeFitEngineSettings(FitEngineSettings fitEngineSettings) {
    this.fitEngineSettings.mergeFrom(fitEngineSettings);
  }

  /**
   * Sets the fit engine settings. After calling this then the object returned from
   * {@link #getFitConfiguration()} must be refreshed.
   *
   * @param fitEngineSettings the new fit engine settings
   */
  public void setFitEngineSettings(FitEngineSettings fitEngineSettings) {
    this.fitEngineSettings.clear().mergeFrom(fitEngineSettings);
    fitConfiguration.updateFitSettings(this.fitEngineSettings.getFitSettingsBuilder());
  }

  /**
   * Gets the fit configuration.
   *
   * @return the fit configuration
   */
  public FitConfiguration getFitConfiguration() {
    return fitConfiguration;
  }

  /**
   * Gets the search parameter.
   *
   * @return the size of the region to search for local maxima. The actual window is calculated
   *         dynamically in conjunction with the peak widths using {@link #getHwhmMax()}.
   */
  public RelativeParameter getSearchParameter() {
    return fitEngineSettings.getSearch();
  }

  /**
   * Gets the size of the region to search for local maxima.
   *
   * @return the size of the region to search for local maxima. The actual window is calculated
   *         dynamically in conjunction with the peak widths using {@link #getHwhmMax()}.
   */
  public double getSearch() {
    return fitEngineSettings.getSearch().getValue();
  }

  /**
   * Sets the size of the region to search for local maxima.
   *
   * @param search the size of the region to search for local maxima. The actual window is
   *        calculated dynamically in conjunction with the peak widths using {@link #getHwhmMax()}.
   */
  public void setSearch(double search) {
    fitEngineSettings.getSearchBuilder().setValue(search);
  }

  /**
   * Sets the size of the region to search for local maxima.
   *
   * @param search the size of the region to search for local maxima. The actual window is
   *        calculated dynamically in conjunction with the peak widths using {@link #getHwhmMax()}.
   * @param absolute the absolute flag
   */
  public void setSearch(double search, boolean absolute) {
    fitEngineSettings.getSearchBuilder().setValue(search).setAbsolute(absolute);
  }

  /**
   * Gets if the Search parameter is absolute.
   *
   * @return True if the Search parameter is absolute
   */
  public boolean getSearchAbsolute() {
    return fitEngineSettings.getSearch().getAbsolute();
  }

  /**
   * Sets if the Search parameter is absolute.
   *
   * @param absolute True if the Search parameter is absolute
   */
  public void setSearchAbsolute(boolean absolute) {
    fitEngineSettings.getSearchBuilder().setAbsolute(absolute);
  }

  /**
   * Gets the size of the border region to ignore.
   *
   * @return the border the size of the border region to ignore. The actual window is calculated
   *         dynamically in conjunction with the peak widths using {@link #getHwhmMax()}.
   */
  public RelativeParameter getBorderParameter() {
    return fitEngineSettings.getBorder();
  }

  /**
   * Gets the size of the border region to ignore.
   *
   * @return the border the size of the border region to ignore. The actual window is calculated
   *         dynamically in conjunction with the peak widths using {@link #getHwhmMax()}.
   */
  public double getBorder() {
    return fitEngineSettings.getBorder().getValue();
  }

  /**
   * Sets the size of the border region to ignore.
   *
   * @param border the size of the border region to ignore. The actual window is calculated
   *        dynamically in conjunction with the peak widths using {@link #getHwhmMax()}.
   */
  public void setBorder(double border) {
    fitEngineSettings.getBorderBuilder().setValue(border);
  }

  /**
   * Sets the size of the border region to ignore.
   *
   * @param border the size of the border region to ignore. The actual window is calculated
   *        dynamically in conjunction with the peak widths using {@link #getHwhmMax()}.
   * @param absolute the absolute flag
   */
  public void setBorder(double border, boolean absolute) {
    fitEngineSettings.getBorderBuilder().setValue(border).setAbsolute(absolute);
  }

  /**
   * Gets if the Border parameter is absolute.
   *
   * @return True if the Border parameter is absolute
   */
  public boolean getBorderAbsolute() {
    return fitEngineSettings.getBorder().getAbsolute();
  }

  /**
   * Sets if the Border parameter is absolute.
   *
   * @param absolute True if the Border parameter is absolute
   */
  public void setBorderAbsolute(boolean absolute) {
    fitEngineSettings.getBorderBuilder().setAbsolute(absolute);
  }

  /**
   * Gets the size of the window used for fitting around a maxima.
   *
   * @return the size of the window used for fitting around a maxima. The actual window is
   *         calculated dynamically in conjunction with the peak widths using {@link #getSdMax()}.
   */
  public RelativeParameter getFittingParameter() {
    return fitEngineSettings.getFitting();
  }

  /**
   * Gets the size of the window used for fitting around a maxima.
   *
   * @return the size of the window used for fitting around a maxima. The actual window is
   *         calculated dynamically in conjunction with the peak widths using {@link #getSdMax()}.
   */
  public double getFitting() {
    return fitEngineSettings.getFitting().getValue();
  }

  /**
   * Sets the size of the window used for fitting around a maxima.
   *
   * @param fitting the size of the window used for fitting around a maxima. The actual window is
   *        calculated dynamically in conjunction with the peak widths using {@link #getSdMax()}.
   */
  public void setFitting(double fitting) {
    fitEngineSettings.getFittingBuilder().setValue(fitting);
  }

  /**
   * Sets the size of the window used for fitting around a maxima.
   *
   * @param fitting the size of the window used for fitting around a maxima. The actual window is
   *        calculated dynamically in conjunction with the peak widths using {@link #getSdMax()}.
   * @param absolute the absolute flag
   */
  public void setFitting(double fitting, boolean absolute) {
    fitEngineSettings.getFittingBuilder().setValue(fitting).setAbsolute(absolute);
  }

  /**
   * Gets if the Fitting parameter is absolute.
   *
   * @return True if the Fitting parameter is absolute
   */
  public boolean getFittingAbsolute() {
    return fitEngineSettings.getFitting().getAbsolute();
  }

  /**
   * Sets if the Fitting parameter is absolute.
   *
   * @param absolute True if the Fitting parameter is absolute
   */
  public void setFittingAbsolute(boolean absolute) {
    fitEngineSettings.getFittingBuilder().setAbsolute(absolute);
  }

  /**
   * Set the failures limit. When failures exceeds the failures limit then stop fitting.
   *
   * @return the failuresLimit
   */
  public int getFailuresLimit() {
    return fitEngineSettings.getFailuresLimit();
  }

  /**
   * Set the failures limit. When failures exceeds the failures limit then stop fitting. i.e.
   * failures=0 will stop on the first failure, failures=1 will stop on the second consecutive
   * failure. If negative then this is disabled and all candidates will be processed.
   *
   * @param failuresLimit the number of consecutive failures that stops the fitting process on the
   *        frame
   */
  public void setFailuresLimit(int failuresLimit) {
    if (failuresLimit != getFailuresLimit()) {
      failCounter = null;
    }
    fitEngineSettings.setFailuresLimit(failuresLimit);
  }

  /**
   * Gets the pass rate. If the fraction of accepted fits falls below this threshold then stop
   * fitting of the remaining candidates.
   *
   * @return the pass rate
   */
  public double getPassRate() {
    return fitEngineSettings.getPassRate();
  }

  /**
   * Sets the pass rate (range 0-1) to continue fitting. If the fraction of accepted fits falls
   * below this threshold then stop fitting of the remaining candidates. Set to zero to disable.
   *
   * @param passRate the new pass rate
   */
  public void setPassRate(double passRate) {
    if (passRate != getPassRate()) {
      failCounter = null;
    }
    fitEngineSettings.setPassRate(passRate);
  }

  /**
   * Reset the fail counter. This disables stopping criteria so that all candidates will be fit.
   */
  public void resetFailCounter() {
    failCounter = null;
    fitEngineSettings.setFailuresLimit(-1);
    fitEngineSettings.setPassRate(0);
  }

  private FailCounter failCounter;

  /**
   * Gets the fail counter.
   *
   * @return the fail counter
   */
  public FailCounter getFailCounter() {
    // TODO - Make this more configurable. This could be done by allowing the fail counter
    // to be serialised to a string and this can be passed in. Or the fail counter could be
    // manually passed in. Either would allow more complex fail counters to be used.
    if (failCounter == null) {
      final int failuresLimit = getFailuresLimit();
      final FailCounter f1 =
          // Negative will disable
          (failuresLimit >= 0) ? ConsecutiveFailCounter.create(failuresLimit) : null;
      final double passRate = getPassRate();
      // TODO - the allowed counts could be an input
      final FailCounter f2 = (passRate > 0) ? PassRateFailCounter.create(5, passRate) : null;

      // All fail counters must pass to continue fitting
      failCounter = CombinedAndFailCounter.join(f1, f2);
      if (failCounter == null) {
        failCounter = NullFailCounter.INSTANCE;
      }
    }
    return failCounter.newCounter();
  }

  /**
   * Checks if is include neighbours.
   *
   * @return true, if is include neighbours
   */
  public boolean isIncludeNeighbours() {
    return fitEngineSettings.getIncludeNeighbours();
  }

  /**
   * Include neighbour maxima in the fitting.
   *
   * <p>Use this option when the fitting search region is large relative to the smoothing, thus
   * other peaks may be within the region used for fitting.
   *
   * @param includeNeighbours the new include neighbours
   */
  public void setIncludeNeighbours(boolean includeNeighbours) {
    fitEngineSettings.setIncludeNeighbours(includeNeighbours);
  }

  /**
   * Gets the neighbour height threshold.
   *
   * @return the neighbour height threshold
   */
  public double getNeighbourHeightThreshold() {
    return fitEngineSettings.getNeighbourHeightThreshold();
  }

  /**
   * Sets the neighbour height threshold.
   *
   * @param neighbourHeightThreshold Set the height threshold that determines if a neighbour peak
   *        should be fitted (specified as a fraction of the central peak relative to the
   *        background)
   */
  public void setNeighbourHeightThreshold(double neighbourHeightThreshold) {
    fitEngineSettings.setNeighbourHeightThreshold(neighbourHeightThreshold);
  }

  /**
   * Gets the residuals threshold.
   *
   * @return the residuals threshold
   */
  public double getResidualsThreshold() {
    return fitEngineSettings.getResidualsThreshold();
  }

  /**
   * Sets the residuals threshold.
   *
   * @param residualsThreshold Set the threshold for the residuals analysis that determines if a
   *        two-kernel model should be fitted
   */
  public void setResidualsThreshold(double residualsThreshold) {
    fitEngineSettings.setResidualsThreshold(residualsThreshold);
  }

  /**
   * Gets the noise method used to estimate the image noise.
   *
   * @return the method used to estimate the image noise
   */
  public NoiseEstimatorMethod getNoiseMethod() {
    return fitEngineSettings.getNoiseMethod();
  }

  /**
   * Sets the noise method used to estimate the image noise.
   *
   * @param noiseMethod Set the method used to estimate the image noise
   */
  public void setNoiseMethod(NoiseEstimatorMethod noiseMethod) {
    fitEngineSettings.setNoiseMethod(noiseMethod);
  }

  /**
   * Sets the noise method used to estimate the image noise.
   *
   * @param noiseMethod Set the method used to estimate the image noise
   */
  public void setNoiseMethod(int noiseMethod) {
    fitEngineSettings.setNoiseMethodValue(noiseMethod);
  }

  /**
   * Sets the duplicate distance within which spots are considered duplicates.
   *
   * @param duplicateDistance The distance within which spots are considered duplicates
   */
  public void setDuplicateDistance(final double duplicateDistance) {
    fitEngineSettings.getDuplicateDistanceBuilder().setValue(duplicateDistance);
  }

  /**
   * Sets if the duplicate distance is absolute.
   *
   * @param absolute True if the duplicate distance is absolute
   */
  public void setDuplicateDistanceAbsolute(boolean absolute) {
    fitEngineSettings.getDuplicateDistanceBuilder().setAbsolute(absolute);
  }

  /**
   * Gets the duplicate distance parameter.
   *
   * @return The distance within which spots are considered duplicates
   */
  public RelativeParameter getDuplicateDistanceParameter() {
    return fitEngineSettings.getDuplicateDistance();
  }

  /**
   * Gets the duplicate distance.
   *
   * @return The distance within which spots are considered duplicates
   */
  public double getDuplicateDistance() {
    return fitEngineSettings.getDuplicateDistance().getValue();
  }

  /**
   * Gets if the duplicate distance is absolute.
   *
   * @return True if the duplicate distance is absolute
   */
  public boolean getDuplicateDistanceAbsolute() {
    return fitEngineSettings.getDuplicateDistance().getAbsolute();
  }

  /**
   * Creates a copy of the configuration.
   *
   * @return the copy
   */
  public FitEngineConfiguration createCopy() {
    // Manual copy using the current proto objects
    final FitEngineConfiguration clone = create(getFitEngineSettings(),
        getFitConfiguration().getCalibration(), getFitConfiguration().getPsf());
    // Copy anything else not in a proto object
    clone.getFitConfiguration().copySettings(getFitConfiguration());
    return clone;
  }

  /**
   * Gets the type of filter to apply to the data before identifying local maxima.
   *
   * @return the type of filter to apply to the data before identifying local maxima
   */
  public DataFilterType getDataFilterType() {
    return fitEngineSettings.getDataFilterSettings().getDataFilterType();
  }

  /**
   * Sets the type of filter to apply to the data before identifying local maxima.
   *
   * @param dataFilterType the type of filter to apply to the data before identifying local maxima
   */
  public void setDataFilterType(int dataFilterType) {
    final DataFilterType t = DataFilterType.forNumber(dataFilterType);
    if (t != null) {
      setDataFilterType(t);
    }
  }

  /**
   * Sets the type of filter to apply to the data before identifying local maxima.
   *
   * @param dataFilterType the type of filter to apply to the data before identifying local maxima
   */
  public void setDataFilterType(DataFilterType dataFilterType) {
    final DataFilterSettings.Builder b = fitEngineSettings.getDataFilterSettingsBuilder();
    b.setDataFilterType(dataFilterType);

    // Resize the filter arrays to discard unused filters
    switch (dataFilterType) {
      case JURY:
        return;

      case DIFFERENCE:
        truncateFilters(b, 2);
        break;

      case SINGLE:
      default:
        truncateFilters(b, 1);
    }
  }

  private static void truncateFilters(DataFilterSettings.Builder builder, int size) {
    while (builder.getDataFiltersCount() > size) {
      builder.removeDataFilters(builder.getDataFiltersCount() - 1);
    }
  }

  /**
   * Gets the data filter count.
   *
   * @return the data filter count
   */
  public int getDataFiltersCount() {
    return this.fitEngineSettings.getDataFilterSettings().getDataFiltersCount();
  }

  /**
   * Gets the data filter method.
   *
   * @param n The filter number
   * @return the filter to apply to the data before identifying local maxima
   */
  public DataFilterMethod getDataFilterMethod(int n) {
    if (n < getDataFiltersCount()) {
      return this.fitEngineSettings.getDataFilterSettings().getDataFilters(n).getDataFilterMethod();
    }
    throw new IndexOutOfBoundsException(n + " >= " + getDataFiltersCount());
  }

  /**
   * Gets the data filter method.
   *
   * @param n The filter number
   * @param defaultValue the default value
   * @return the filter to apply to the data before identifying local maxima
   */
  public DataFilterMethod getDataFilterMethod(int n, DataFilterMethod defaultValue) {
    if (n < getDataFiltersCount()) {
      return this.fitEngineSettings.getDataFilterSettings().getDataFilters(n).getDataFilterMethod();
    }
    return defaultValue;
  }

  /**
   * Gets the data filter parameter.
   *
   * @param n The filter number
   * @return the filter parameter
   */
  public RelativeParameter getDataFilterParameter(int n) {
    if (n < getDataFiltersCount()) {
      return this.fitEngineSettings.getDataFilterSettings().getDataFilters(n).getParameters(0);
    }
    throw new IndexOutOfBoundsException(n + " >= " + getDataFiltersCount());
  }

  /**
   * Gets if the smoothing parameter is absolute.
   *
   * @param n The filter number
   * @return if the smoothing parameter is absolute
   */
  public boolean getDataFilterParameterAbsolute(int n) {
    if (n < getDataFiltersCount()) {
      return this.fitEngineSettings.getDataFilterSettings().getDataFilters(n).getParameters(0)
          .getAbsolute();
    }
    throw new IndexOutOfBoundsException(n + " >= " + getDataFiltersCount());
  }

  /**
   * Gets if the smoothing parameter is absolute.
   *
   * @param n The filter number
   * @param defaultValue the default value
   * @return if the smoothing parameter is absolute
   */
  public boolean getDataFilterParameterAbsolute(int n, boolean defaultValue) {
    if (n < getDataFiltersCount()) {
      return this.fitEngineSettings.getDataFilterSettings().getDataFilters(n).getParameters(0)
          .getAbsolute();
    }
    return defaultValue;
  }

  /**
   * Gets the data filter parameter value.
   *
   * @param n The filter number
   * @return the smoothing window size
   */
  public double getDataFilterParameterValue(int n) {
    if (n < getDataFiltersCount()) {
      return this.fitEngineSettings.getDataFilterSettings().getDataFilters(n).getParameters(0)
          .getValue();
    }
    throw new IndexOutOfBoundsException(n + " >= " + getDataFiltersCount());
  }

  /**
   * Gets smoothing window size.
   *
   * @param n The filter number
   * @param defaultValue the default value
   * @return the smoothing window size
   */
  public double getDataFilterParameterValue(int n, double defaultValue) {
    if (n < getDataFiltersCount()) {
      return this.fitEngineSettings.getDataFilterSettings().getDataFilters(n).getParameters(0)
          .getValue();
    }
    return defaultValue;
  }

  /**
   * Sets the data filter.
   *
   * @param dataFilterMethod the filter to apply to the data before identifying local maxima
   * @param smooth the size of the smoothing window. The actual window is calculated dynamically in
   *        conjunction with the peak widths.
   * @param absolute the absolute
   * @param n The filter number
   */
  public void setDataFilter(DataFilterMethod dataFilterMethod, double smooth, boolean absolute,
      int n) {
    final DataFilterSettings.Builder b = fitEngineSettings.getDataFilterSettingsBuilder();
    final DataFilter.Builder b2 =
        (b.getDataFiltersCount() == n) ? b.addDataFiltersBuilder() : b.getDataFiltersBuilder(n);
    b2.setDataFilterMethod(dataFilterMethod);
    b2.clearParameters();
    final RelativeParameter.Builder b3 = b2.addParametersBuilder();
    b3.setValue(smooth);
    b3.setAbsolute(absolute);
  }

  /**
   * Sets the data filter without the absolute flag. If no parameter exists then the absolute flag
   * will have the default value (false).
   *
   * @param dataFilterMethod the filter to apply to the data before identifying local maxima
   * @param smooth the size of the smoothing window. The actual window is calculated dynamically in
   *        conjunction with the peak widths.
   * @param n The filter number
   */
  public void setDataFilter(DataFilterMethod dataFilterMethod, double smooth, int n) {
    final DataFilterSettings.Builder b = fitEngineSettings.getDataFilterSettingsBuilder();
    final DataFilter.Builder b2 =
        (b.getDataFiltersCount() == n) ? b.addDataFiltersBuilder() : b.getDataFiltersBuilder(n);
    b2.setDataFilterMethod(dataFilterMethod);
    final RelativeParameter.Builder b3 =
        (b2.getParametersCount() == 0) ? b2.addParametersBuilder() : b2.getParametersBuilder(0);
    b3.setValue(smooth);
  }

  /**
   * Sets the data filter.
   *
   * @param dataFilter the data filter
   * @param smooth the size of the smoothing window. The actual window is calculated dynamically in
   *        conjunction with the peak widths using {@link #getHwhmMin()}.
   * @param absolute Set to true to use the absolute value. Otherwise compute relative to
   *        {@link #getHwhmMin()}.
   * @param n The filter number
   */
  public void setDataFilter(int dataFilter, double smooth, boolean absolute, int n) {
    final DataFilterMethod m = DataFilterMethod.forNumber(dataFilter);
    if (m != null) {
      setDataFilter(m, smooth, absolute, n);
    }
  }

  /**
   * Sets the data filter without the absolute flag. If no parameter exists then the absolute flag
   * will have the default value (false).
   *
   * @param dataFilter the data filter
   * @param smooth the size of the smoothing window. The actual window is calculated dynamically in
   *        conjunction with the peak widths using {@link #getHwhmMin()}.
   * @param n The filter number
   */
  public void setDataFilter(int dataFilter, double smooth, int n) {
    final DataFilterMethod m = DataFilterMethod.forNumber(dataFilter);
    if (m != null) {
      setDataFilter(m, smooth, n);
    }
  }

  /**
   * Set the number of filters to use. Call this method when all the filters have been set to clear
   * any other stored filters from memory.
   *
   * @param n The number of filters
   */
  public void setNumberOfFilters(int n) {
    if (n < 1) {
      n = 1;
    }
    final DataFilterSettings.Builder b = fitEngineSettings.getDataFilterSettingsBuilder();
    // Truncate the filter
    truncateFilters(b, n);
  }

  /**
   * Get the fitting width. If the fitting parameter is relative this is calculated using the
   * maximum peak standard deviation multiplied by the fitting parameter, otherwise the absolute
   * value is used. The value is then rounded up to the next integer and a minimum of 2 is returned.
   *
   * @return The fitting width
   */
  public int getFittingWidth() {
    // Region for peak fitting
    int fitting = (int) Math.ceil(convertUsingSdMax(getFittingParameter()));
    if (fitting < 2) {
      fitting = 2;
    }
    return fitting;
  }

  /**
   * Create the spot filter for identifying candidate maxima. The actual border, search width and
   * smoothing parameters can be configured relative to the configured standard deviations or left
   * absolute. The standard deviation is used to determine the Half-Width at Half-Maximum (HWHM) for
   * each dimension and the parameters set as follows.
   *
   * <pre>
   * int search = (int) Math.ceil(getSearch() * hwhmMax);
   * int border = (int) Math.floor(getBorder() * hwhmMax);
   * // For each filter
   * double smooth = getSmooth(i) * hwhmMin;
   * </pre>
   *
   * @return the maxima spot filter
   */
  public MaximaSpotFilter createSpotFilter() {
    // Respect the absolute parameter absolute flag for all the distances.

    // Get the half-width at half maximum
    final double hwhmMin = getHwhmMin();
    final double hwhmMax = getHwhmMax();

    // Note: rounding to 2 decimal places is a simple method for removing small errors
    // in floating point precision from creating an incorrect integer

    // Region for maxima finding
    int search = (int) Math.ceil(convert(getSearchParameter(), hwhmMax, 2));
    if (search < 1) {
      search = 1;
    }

    // Border where peaks are ignored
    int border = (int) Math.floor(convert(getBorderParameter(), hwhmMax, 2));
    if (border < 0) {
      border = 0;
    }

    final DataProcessor processor0 = createDataProcessor(border, 0, hwhmMin);
    final DataFilterSettings f = fitEngineSettings.getDataFilterSettings();
    final int filterCount = f.getDataFiltersCount();

    final MaximaSpotFilter spotFilter;
    if (f.getDataFilterType() == DataFilterType.JURY && filterCount > 1) {
      final DataProcessor[] processors = new DataProcessor[filterCount];
      processors[0] = processor0;
      for (int i = 1; i < filterCount; i++) {
        processors[i] = createDataProcessor(border, i, hwhmMin);
      }
      spotFilter = new JurySpotFilter(search, border, processors);
    } else if (f.getDataFilterType() == DataFilterType.DIFFERENCE && filterCount > 1) {
      final DataProcessor processor1 = createDataProcessor(border, 1, hwhmMin);
      spotFilter = new DifferenceSpotFilter(search, border, processor0, processor1);
    } else {
      spotFilter = new SingleSpotFilter(search, border, processor0);
    }

    if (getFitConfiguration().isPerPixelCameraType()) {
      if (!spotFilter.isWeighted()) {
        throw new IllegalStateException(
            "Camera type requires a weighted spot filter: " + fitConfiguration.getCameraType());
      }
      final CameraModel model = fitConfiguration.getCameraModel();
      if (model == null || !model.isPerPixelModel()) {
        throw new IllegalStateException("Weighted spot filter requires a per-pixel camera model");
      }
    }

    return spotFilter;
  }

  /**
   * Gets the minimum SD using the initial standard deviations.
   *
   * <p>Note: This uses initial peak SD0 and optionally initial peak SD1 if width fitting is enabled
   * for the second dimension.
   *
   * @return the SD min
   */
  public double getSdMin() {
    final double[] w = fitConfiguration.getGaussian2DWxWy();

    final double initialPeakStdDev0 = w[0];
    final double initialPeakStdDev1 = w[1];

    // Use 1 if zero to get at least a single pixel width
    double widthMin = (initialPeakStdDev0 > 0) ? initialPeakStdDev0 : 1;

    // Only use the second width if this is part of the function.
    // This should be taken care of within the PSF helper.
    if (initialPeakStdDev1 > 0) {
      widthMin = Math.min(initialPeakStdDev1, widthMin);
    }

    return widthMin;
  }

  /**
   * Gets the minimum HWHM using the initial standard deviations.
   *
   * <p>Note: This uses initial peak SD0 and optionally initial peak SD1 if width fitting is enabled
   * for the second dimension.
   *
   * @return the HWHM min
   */
  public double getHwhmMin() {
    // Get the half-width at half maximim
    return Gaussian2DFunction.SD_TO_HWHM_FACTOR * getSdMin();
  }

  /**
   * Gets the maximum SD using the initial standard deviations.
   *
   * <p>Note: This uses initial peak SD0 and optionally initial peak SD1 if width fitting is enabled
   * for the second dimension.
   *
   * @return the SD max
   */
  public double getSdMax() {
    final double[] w = fitConfiguration.getGaussian2DWxWy();

    final double initialPeakStdDev0 = w[0];
    final double initialPeakStdDev1 = w[1];

    // Use 1 if zero to get at least a single pixel width
    double widthMax = (initialPeakStdDev0 > 0) ? initialPeakStdDev0 : 1;

    // Only use the second width if this is part of the function.
    // This should be taken care of within the PSF helper.
    if (initialPeakStdDev1 > 0) {
      widthMax = Math.max(initialPeakStdDev1, widthMax);
    }

    return widthMax;
  }

  /**
   * Gets the maximum HWHM using the initial standard deviations.
   *
   * <p>Note: This uses initial peak SD0 and optionally initial peak SD1 if width fitting is enabled
   * for the second dimension.
   *
   * @return the HWHM max
   */
  public double getHwhmMax() {
    // Get the half-width at half maximim
    return Gaussian2DFunction.SD_TO_HWHM_FACTOR * getSdMax();
  }

  /**
   * Convert the relative parameter using the scale. If the parameter is absolute then return the
   * unscaled value.
   *
   * @param parameter the parameter
   * @param scale the scale
   * @param decimalPlaces the decimal places to round the scaled number (set to negative to ignore)
   * @return the double
   */
  public static double convert(RelativeParameter parameter, double scale, int decimalPlaces) {
    return (parameter.getAbsolute()) ? parameter.getValue()
        : round(parameter.getValue() * scale, decimalPlaces);
  }

  /**
   * Round to the given number of decimal places.
   *
   * @param value the value
   * @param decimalPlaces the decimal places
   * @return the double
   */
  public static double round(double value, int decimalPlaces) {
    return (decimalPlaces >= 0) ? MathUtils.roundUsingDecimalPlaces(value, decimalPlaces) : value;
  }

  /**
   * Convert using half-width half-max as the scale.
   *
   * @param parameter the parameter
   * @return the converted value
   */
  public double convertUsingHwhMax(RelativeParameter parameter) {
    return (parameter.getAbsolute()) ? parameter.getValue() : parameter.getValue() * getHwhmMax();
  }

  /**
   * Convert using half-width half-min as the scale.
   *
   * @param parameter the parameter
   * @return the converted value
   */
  public double convertUsingHwhMin(RelativeParameter parameter) {
    return (parameter.getAbsolute()) ? parameter.getValue() : parameter.getValue() * getHwhmMin();
  }

  /**
   * Convert using SD max as the scale.
   *
   * @param parameter the parameter
   * @return the converted value
   */
  public double convertUsingSdMax(RelativeParameter parameter) {
    return (parameter.getAbsolute()) ? parameter.getValue() : parameter.getValue() * getSdMax();
  }

  /**
   * Gets the number of filters for the configured filter type.
   *
   * @return the number of filters
   */
  public int getNumberOfFilters() {
    final DataFilterSettings f = fitEngineSettings.getDataFilterSettings();
    final int filterCount = f.getDataFiltersCount();
    if (f.getDataFilterType() == DataFilterType.JURY && filterCount > 1) {
      return filterCount;
    } else if (f.getDataFilterType() == DataFilterType.DIFFERENCE && filterCount > 1) {
      return 2;
    } else {
      // Single filter
      return 1;
    }
  }

  private static double getSmoothingWindow(DataFilter filter, double hwhmMin) {
    final RelativeParameter rp = filter.getParameters(0);

    // Q. Why is this rounded. Is it just to make a nicer number?
    return convert(rp, hwhmMin, 2);
  }

  private DataProcessor createDataProcessor(int border, int n, double hwhm) {
    if (n < this.fitEngineSettings.getDataFilterSettings().getDataFiltersCount()) {
      final DataFilter f = fitEngineSettings.getDataFilterSettings().getDataFilters(n);
      return createDataProcessor(border, f.getDataFilterMethod(), getSmoothingWindow(f, hwhm));
    }
    return null;
  }

  /**
   * Create a data processor for the spot filter.
   *
   * @param border the border
   * @param dataFilter the data filter
   * @param parameter the parameter
   * @return the data processor
   */
  public static DataProcessor createDataProcessor(int border, DataFilterMethod dataFilter,
      double parameter) {
    switch (dataFilter) {
      case MEAN:
        return new AverageDataProcessor(border, parameter);

      case BLOCK_MEAN:
        return new BlockAverageDataProcessor(border, parameter);

      case CIRCULAR_MEAN:
        return new CircularMeanDataProcessor(border, parameter);

      case MEDIAN:
        return new MedianDataProcessor(border, parameter);

      case GAUSSIAN:
        return new GaussianDataProcessor(border, parameter);

      default:
        throw new NotImplementedException("Not yet implemented: " + dataFilter.toString());
    }
  }

  /**
   * Copy data filter.
   *
   * @param config the config
   */
  public void copyDataFilter(FitEngineConfiguration config) {
    fitEngineSettings.setDataFilterSettings(config.fitEngineSettings.getDataFilterSettings());
  }

  /**
   * Configure the output units from fitting using the current calibration and fit solver settings.
   *
   * <p>This method should be called before the calibration is passed to any object that will handle
   * the fitting output.
   *
   * <p>It will update the calibration units and the precision method to match that used by
   * precision filter.
   */
  public void configureOutputUnits() {
    final FitConfiguration fitConfig = getFitConfiguration();
    // If there is no calibration then the writer will just have the defaults
    final CalibrationWriter calibration = fitConfig.getCalibrationWriterReference();

    // Fitting is always done pixels and radians
    calibration.setDistanceUnit(DistanceUnit.PIXEL);
    calibration.setAngleUnit(AngleUnit.RADIAN);

    // Most fitters fit in photons unless we have no calibration.
    IntensityUnit intensityUnit = IntensityUnit.PHOTON;

    if (fitConfig.isFitCameraCounts()) {
      intensityUnit = IntensityUnit.COUNT;
    }

    calibration.setIntensityUnit(intensityUnit);

    // This initialises the calibration precision method
    fitConfig.getFilterPrecisionMethod();
  }
}
