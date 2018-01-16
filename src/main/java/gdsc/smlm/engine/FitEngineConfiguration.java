package gdsc.smlm.engine;

import org.apache.commons.math3.util.FastMath;

import gdsc.core.utils.Maths;
import gdsc.smlm.data.config.CalibrationProtos.Calibration;
import gdsc.smlm.data.config.CalibrationProtosHelper;
import gdsc.smlm.data.config.CalibrationWriter;
import gdsc.smlm.data.config.FitProtos.DataFilter;
import gdsc.smlm.data.config.FitProtos.DataFilterMethod;
import gdsc.smlm.data.config.FitProtos.DataFilterSettings;
import gdsc.smlm.data.config.FitProtos.DataFilterType;
import gdsc.smlm.data.config.FitProtos.FitEngineSettings;
import gdsc.smlm.data.config.FitProtos.NoiseEstimatorMethod;
import gdsc.smlm.data.config.FitProtos.RelativeParameter;
import gdsc.smlm.data.config.FitProtosHelper;
import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.data.config.PSFProtos.PSF;
import gdsc.smlm.data.config.PSFProtosHelper;
import gdsc.smlm.data.config.UnitProtos.AngleUnit;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import gdsc.smlm.filters.AverageDataProcessor;
import gdsc.smlm.filters.BlockAverageDataProcessor;
import gdsc.smlm.filters.CircularMeanDataProcessor;
import gdsc.smlm.filters.DataProcessor;
import gdsc.smlm.filters.DifferenceSpotFilter;
import gdsc.smlm.filters.GaussianDataProcessor;
import gdsc.smlm.filters.JurySpotFilter;
import gdsc.smlm.filters.MaximaSpotFilter;
import gdsc.smlm.filters.MedianDataProcessor;
import gdsc.smlm.filters.SingleSpotFilter;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.model.camera.CameraModel;
import gdsc.smlm.results.count.CombinedAndFailCounter;
import gdsc.smlm.results.count.ConsecutiveFailCounter;
import gdsc.smlm.results.count.FailCounter;
import gdsc.smlm.results.count.NullFailCounter;
import gdsc.smlm.results.count.PassRateFailCounter;

/**
 * Specifies the configuration for the fit engine
 */
public class FitEngineConfiguration implements Cloneable
{
	private FitEngineSettings.Builder fitEngineSettings;
	private FitConfiguration fitConfiguration = null;

	/**
	 * Instantiates a new fit engine configuration.
	 */
	public FitEngineConfiguration()
	{
		this(FitProtosHelper.defaultFitEngineSettings, CalibrationProtosHelper.defaultCalibration,
				PSFProtosHelper.defaultOneAxisGaussian2DPSF);
	}

	/**
	 * Instantiates a new fit engine configuration.
	 *
	 * @param fitEngineSettings
	 *            the fit engine settings
	 * @param calibration
	 *            the calibration
	 * @param psf
	 *            the psf
	 */
	public FitEngineConfiguration(FitEngineSettings fitEngineSettings, Calibration calibration, PSF psf)
	{
		if (fitEngineSettings == null)
			throw new IllegalArgumentException("FitEngineSettings is null");
		if (calibration == null)
			throw new IllegalArgumentException("Calibration is null");
		if (psf == null)
			throw new IllegalArgumentException("PSF is null");
		init(fitEngineSettings.toBuilder(), calibration.toBuilder(), psf.toBuilder());
	}

	/**
	 * Instantiates a new fit engine configuration.
	 *
	 * @param fitEngineSettings
	 *            the fit engine settings
	 * @param calibration
	 *            the calibration
	 * @param psf
	 *            the psf
	 */
	public FitEngineConfiguration(FitEngineSettings.Builder fitEngineSettings, Calibration.Builder calibration,
			PSF.Builder psf)
	{
		if (fitEngineSettings == null)
			throw new IllegalArgumentException("FitEngineSettings is null");
		if (calibration == null)
			throw new IllegalArgumentException("Calibration is null");
		if (psf == null)
			throw new IllegalArgumentException("PSF is null");
		init(fitEngineSettings, calibration, psf);
	}

	private void init(FitEngineSettings.Builder fitEngineSettings, Calibration.Builder calibration, PSF.Builder psf)
	{
		this.fitEngineSettings = fitEngineSettings;
		fitConfiguration = new FitConfiguration(fitEngineSettings.getFitSettingsBuilder(), calibration, psf, false);
	}

	/**
	 * Gets the fit engine settings.
	 *
	 * @return the fit engine settings
	 */
	public FitEngineSettings getFitEngineSettings()
	{
		return fitEngineSettings.build();
	}

	/**
	 * Merge fit engine settings.
	 *
	 * @param fitEngineSettings
	 *            the fit engine settings
	 */
	public void mergeFitEngineSettings(FitEngineSettings fitEngineSettings)
	{
		this.fitEngineSettings.mergeFrom(fitEngineSettings);
	}

	/**
	 * Sets the fit engine settings. After calling this then the object returned from {@link #getFitConfiguration()}
	 * must be refreshed.
	 *
	 * @param fitEngineSettings
	 *            the new fit engine settings
	 */
	public void setFitEngineSettings(FitEngineSettings fitEngineSettings)
	{
		this.fitEngineSettings.clear().mergeFrom(fitEngineSettings);
		fitConfiguration.updateFitSettings(this.fitEngineSettings.getFitSettingsBuilder());
	}

	/**
	 * Gets the fit configuration.
	 *
	 * @return the fit configuration
	 */
	public FitConfiguration getFitConfiguration()
	{
		return fitConfiguration;
	}

	/**
	 * @return the size of the region to search for local maxima. The actual window is calculated dynamically
	 *         in conjunction with the peak widths using {@link #getHWHMMax()}.
	 */
	public RelativeParameter getSearchParameter()
	{
		return fitEngineSettings.getSearch();
	}

	/**
	 * @return the size of the region to search for local maxima. The actual window is calculated dynamically
	 *         in conjunction with the peak widths using {@link #getHWHMMax()}.
	 */
	public double getSearch()
	{
		return fitEngineSettings.getSearch().getValue();
	}

	/**
	 * @param search
	 *            the size of the region to search for local maxima. The actual window is calculated dynamically
	 *            in conjunction with the peak widths using {@link #getHWHMMax()}.
	 */
	public void setSearch(double search)
	{
		fitEngineSettings.getSearchBuilder().setValue(search);
	}

	/**
	 * Sets the search.
	 *
	 * @param search
	 *            the size of the region to search for local maxima. The actual window is calculated dynamically
	 *            in conjunction with the peak widths using {@link #getHWHMMax()}.
	 * @param absolute
	 *            the absolute flag
	 */
	public void setSearch(double search, boolean absolute)
	{
		fitEngineSettings.getSearchBuilder().setValue(search).setAbsolute(absolute);
	}

	/**
	 * @return True if the Search parameter is absolute
	 */
	public boolean getSearchAbsolute()
	{
		return fitEngineSettings.getSearch().getAbsolute();
	}

	/**
	 * @param absolute
	 *            True if the Search parameter is absolute
	 */
	public void setSearchAbsolute(boolean absolute)
	{
		fitEngineSettings.getSearchBuilder().setAbsolute(absolute);
	}

	/**
	 * @return the border
	 *         the size of the border region to ignore. The actual window is calculated dynamically
	 *         in conjunction with the peak widths using {@link #getHWHMMax()}.
	 */
	public RelativeParameter getBorderParameter()
	{
		return fitEngineSettings.getBorder();
	}

	/**
	 * @return the border
	 *         the size of the border region to ignore. The actual window is calculated dynamically
	 *         in conjunction with the peak widths using {@link #getHWHMMax()}.
	 */
	public double getBorder()
	{
		return fitEngineSettings.getBorder().getValue();
	}

	/**
	 * @param border
	 *            the size of the border region to ignore. The actual window is calculated dynamically
	 *            in conjunction with the peak widths using {@link #getHWHMMax()}.
	 */
	public void setBorder(double border)
	{
		fitEngineSettings.getBorderBuilder().setValue(border);
	}

	/**
	 * Sets the border.
	 *
	 * @param border
	 *            the size of the border region to ignore. The actual window is calculated dynamically
	 *            in conjunction with the peak widths using {@link #getHWHMMax()}.
	 * @param absolute
	 *            the absolute flag
	 */
	public void setBorder(double border, boolean absolute)
	{
		fitEngineSettings.getBorderBuilder().setValue(border).setAbsolute(absolute);
	}

	/**
	 * @return True if the Border parameter is absolute
	 */
	public boolean getBorderAbsolute()
	{
		return fitEngineSettings.getBorder().getAbsolute();
	}

	/**
	 * @param absolute
	 *            True if the Border parameter is absolute
	 */
	public void setBorderAbsolute(boolean absolute)
	{
		fitEngineSettings.getBorderBuilder().setAbsolute(absolute);
	}

	/**
	 * @return the size of the window used for fitting around a maxima. The actual window is calculated dynamically
	 *         in conjunction with the peak widths using {@link #getSDMax()}.
	 */
	public RelativeParameter getFittingParameter()
	{
		return fitEngineSettings.getFitting();
	}

	/**
	 * @return the size of the window used for fitting around a maxima. The actual window is calculated dynamically
	 *         in conjunction with the peak widths using {@link #getSDMax()}.
	 */
	public double getFitting()
	{
		return fitEngineSettings.getFitting().getValue();
	}

	/**
	 * @param fitting
	 *            the size of the window used for fitting around a maxima. The actual window is calculated dynamically
	 *            in conjunction with the peak widths using {@link #getSDMax()}.
	 */
	public void setFitting(double fitting)
	{
		fitEngineSettings.getFittingBuilder().setValue(fitting);
	}

	/**
	 * Sets the fitting.
	 *
	 * @param fitting
	 *            the size of the window used for fitting around a maxima. The actual window is calculated dynamically
	 *            in conjunction with the peak widths using {@link #getSDMax()}.
	 * @param absolute
	 *            the absolute flag
	 */
	public void setFitting(double fitting, boolean absolute)
	{
		fitEngineSettings.getFittingBuilder().setValue(fitting).setAbsolute(absolute);
	}

	/**
	 * @return True if the Fitting parameter is absolute
	 */
	public boolean getFittingAbsolute()
	{
		return fitEngineSettings.getFitting().getAbsolute();
	}

	/**
	 * @param absolute
	 *            True if the Fitting parameter is absolute
	 */
	public void setFittingAbsolute(boolean absolute)
	{
		fitEngineSettings.getFittingBuilder().setAbsolute(absolute);
	}

	/**
	 * Set the failures limit. When failures exceeds the failures limit then stop fitting.
	 * 
	 * @return the failuresLimit
	 */
	public int getFailuresLimit()
	{
		return fitEngineSettings.getFailuresLimit();
	}

	/**
	 * Set the failures limit. When failures exceeds the failures limit then stop fitting.
	 * i.e. failures=0 will stop on the first failure, failures=1 will stop on the second consecutive failure.
	 * If negative then this is disabled and all candidates will be processed.
	 * 
	 * @param failuresLimit
	 *            the number of consecutive failures that stops the fitting process on the frame
	 */
	public void setFailuresLimit(int failuresLimit)
	{
		if (failuresLimit != getFailuresLimit())
			failCounter = null;
		fitEngineSettings.setFailuresLimit(failuresLimit);
	}

	/**
	 * Gets the pass rate. If the fraction of accepted fits falls below this threshold then stop fitting of the
	 * remaining candidates.
	 *
	 * @return the pass rate
	 */
	public double getPassRate()
	{
		return fitEngineSettings.getPassRate();
	}

	/**
	 * Sets the pass rate (range 0-1) to continue fitting. If the fraction of accepted fits falls below this threshold
	 * then stop fitting of the remaining candidates. Set to zero to disable.
	 *
	 * @param passRate
	 *            the new pass rate
	 */
	public void setPassRate(double passRate)
	{
		if (passRate != getPassRate())
			failCounter = null;
		fitEngineSettings.setPassRate(passRate);
	}

	/**
	 * Reset the fail counter. This disables stopping criteria so that all candidates will be fit.
	 */
	public void resetFailCounter()
	{
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
	public FailCounter getFailCounter()
	{
		// TODO - Make this more configurable. This could be done by allowing the fail counter
		// to be serialised to a string and this can be passed in. Or the fail counter could be
		// manually passed in. Either would allow more complex fail counters to be used.
		if (failCounter == null)
		{
			int failuresLimit = getFailuresLimit();
			FailCounter f1 = (failuresLimit >= 1) ? ConsecutiveFailCounter.create(failuresLimit) : null;
			double passRate = getPassRate();
			// TODO - the allowed counts could be an input 
			FailCounter f2 = (passRate > 0) ? PassRateFailCounter.create(5, passRate) : null;

			// All fail counters must pass to continue fitting
			failCounter = CombinedAndFailCounter.join(f1, f2);
			if (failCounter == null)
				failCounter = NullFailCounter.INSTANCE;
		}
		return failCounter.newCounter();
	}

	/**
	 * @return the includeNeighbours
	 */
	public boolean isIncludeNeighbours()
	{
		return fitEngineSettings.getIncludeNeighbours();
	}

	/**
	 * Include neighbour maxima in the fitting.
	 * <p>
	 * Use this option when the fitting search region is large relative the the smoothing, thus other peaks may be
	 * within the region used for fitting.
	 * 
	 * @param includeNeighbours
	 */
	public void setIncludeNeighbours(boolean includeNeighbours)
	{
		fitEngineSettings.setIncludeNeighbours(includeNeighbours);
	}

	/**
	 * @return the neighbourHeightThreshold
	 */
	public double getNeighbourHeightThreshold()
	{
		return fitEngineSettings.getNeighbourHeightThreshold();
	}

	/**
	 * @param neighbourHeightThreshold
	 *            Set the height threshold that determines if a neighbour peak should be fitted (specified as a fraction
	 *            of the central peak relative to the background)
	 */
	public void setNeighbourHeightThreshold(double neighbourHeightThreshold)
	{
		fitEngineSettings.setNeighbourHeightThreshold(neighbourHeightThreshold);
	}

	/**
	 * @return the residuals threshold
	 */
	public double getResidualsThreshold()
	{
		return fitEngineSettings.getResidualsThreshold();
	}

	/**
	 * @param residualsThreshold
	 *            Set the threshold for the residuals analysis that determines if a two-kernel model should be fitted
	 */
	public void setResidualsThreshold(double residualsThreshold)
	{
		fitEngineSettings.setResidualsThreshold(residualsThreshold);
	}

	/**
	 * @return the method used to estimate the image noise
	 */
	public NoiseEstimatorMethod getNoiseMethod()
	{
		return fitEngineSettings.getNoiseMethod();
	}

	/**
	 * @param noiseMethod
	 *            Set the method used to estimate the image noise
	 */
	public void setNoiseMethod(NoiseEstimatorMethod noiseMethod)
	{
		fitEngineSettings.setNoiseMethod(noiseMethod);
	}

	/**
	 * @param noiseMethod
	 *            Set the method used to estimate the image noise
	 */
	public void setNoiseMethod(int noiseMethod)
	{
		fitEngineSettings.setNoiseMethodValue(noiseMethod);
	}

	/**
	 * @param duplicateDistance
	 *            The distance within which spots are considered duplicates
	 */
	public void setDuplicateDistance(final double duplicateDistance)
	{
		fitEngineSettings.getDuplicateDistanceBuilder().setValue(duplicateDistance);
	}

	/**
	 * Sets the duplicate distance absolute flas.
	 *
	 * @param absolute
	 *            True if the duplicate distance is absolute
	 */
	public void setDuplicateDistanceAbsolute(boolean absolute)
	{
		fitEngineSettings.getDuplicateDistanceBuilder().setAbsolute(absolute);
	}

	/**
	 * @return The distance within which spots are considered duplicates
	 */
	public RelativeParameter getDuplicateDistanceParameter()
	{
		return fitEngineSettings.getDuplicateDistance();
	}

	/**
	 * @return The distance within which spots are considered duplicates
	 */
	public double getDuplicateDistance()
	{
		return fitEngineSettings.getDuplicateDistance().getValue();
	}

	/**
	 * @return True if the duplicate distance is absolute
	 */
	public boolean getDuplicateDistanceAbsolute()
	{
		return fitEngineSettings.getDuplicateDistance().getAbsolute();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public FitEngineConfiguration clone()
	{
		// Manual copy using the current proto objects
		FitEngineConfiguration clone = new FitEngineConfiguration(getFitEngineSettings(),
				getFitConfiguration().getCalibration(), getFitConfiguration().getPSF());
		// Copy anything else not in a proto object
		clone.getFitConfiguration().copySettings(getFitConfiguration());
		return clone;

		//		// This is not a complete duplicate. The settings builder objects with the 
		//		// underlying configuration will be the same between all instances.
		//		try
		//		{
		//			FitEngineConfiguration f = (FitEngineConfiguration) super.clone();
		//			// Ensure the object is duplicated and not passed by reference.
		//			if (fitConfiguration != null)
		//				f.fitConfiguration = fitConfiguration.clone();
		//			return f;
		//		}
		//		catch (CloneNotSupportedException e)
		//		{
		//			// Ignore
		//		}
		//		return null;
	}

	/**
	 * @return the type of filter to apply to the data before identifying local maxima
	 */
	public DataFilterType getDataFilterType()
	{
		return fitEngineSettings.getDataFilterSettings().getDataFilterType();
	}

	/**
	 * @param DataFilterType
	 *            the type of filter to apply to the data before identifying local maxima
	 */
	public void setDataFilterType(DataFilterType dataFilterType)
	{
		DataFilterSettings.Builder b = fitEngineSettings.getDataFilterSettingsBuilder();
		b.setDataFilterType(dataFilterType);

		// Resize the filter arrays to discard unused filters
		final int n;
		switch (dataFilterType)
		{
			case JURY:
				return;

			case DIFFERENCE:
				n = 2;
				break;

			case SINGLE:
			default:
				n = 1;
		}
		truncateFilters(b, n);
	}

	private void truncateFilters(DataFilterSettings.Builder b, int n)
	{
		while (b.getDataFiltersCount() > n)
			b.removeDataFilters(b.getDataFiltersCount() - 1);
	}

	/**
	 * @param DataFilterType
	 *            the type of filter to apply to the data before identifying local maxima
	 */
	public void setDataFilterType(int dataFilterType)
	{
		DataFilterType t = DataFilterType.forNumber(dataFilterType);
		if (t != null)
		{
			setDataFilterType(t);
		}
	}

	//	/**
	//	 * @param n
	//	 *            The filter number
	//	 * @return the filter to apply to the data before identifying local maxima
	//	 */
	//	public DataFilter getDataFilter(int n)
	//	{
	//		if (n < this.fitEngineSettings.getDataFilterSettings().getDataFiltersCount())
	//			return this.fitEngineSettings.getDataFilterSettings().getDataFilter(n);
	//		return null;
	//	}

	/**
	 * Gets the data filter count.
	 *
	 * @return the data filter count
	 */
	public int getDataFiltersCount()
	{
		return this.fitEngineSettings.getDataFilterSettings().getDataFiltersCount();
	}

	/**
	 * @param n
	 *            The filter number
	 * @return the filter to apply to the data before identifying local maxima
	 */
	public DataFilterMethod getDataFilterMethod(int n)
	{
		if (n < getDataFiltersCount())
			return this.fitEngineSettings.getDataFilterSettings().getDataFilters(n).getDataFilterMethod();
		throw new IndexOutOfBoundsException(n + " >= " + getDataFiltersCount());
	}

	/**
	 * Gets the data filter method.
	 *
	 * @param n
	 *            The filter number
	 * @param defaultValue
	 *            the default value
	 * @return the filter to apply to the data before identifying local maxima
	 */
	public DataFilterMethod getDataFilterMethod(int n, DataFilterMethod defaultValue)
	{
		if (n < getDataFiltersCount())
			return this.fitEngineSettings.getDataFilterSettings().getDataFilters(n).getDataFilterMethod();
		return defaultValue;
	}

	/**
	 * @param n
	 *            The filter number
	 * @return the filter parameter
	 */
	public RelativeParameter getDataFilterParameter(int n)
	{
		if (n < getDataFiltersCount())
			return this.fitEngineSettings.getDataFilterSettings().getDataFilters(n).getParameters(0);
		throw new IndexOutOfBoundsException(n + " >= " + getDataFiltersCount());
	}

	/**
	 * @param n
	 *            The filter number
	 * @return if the smoothing parameter is absolute
	 */
	public boolean getDataFilterParameterAbsolute(int n)
	{
		if (n < getDataFiltersCount())
			return this.fitEngineSettings.getDataFilterSettings().getDataFilters(n).getParameters(0).getAbsolute();
		throw new IndexOutOfBoundsException(n + " >= " + getDataFiltersCount());
	}

	/**
	 * Gets the data filter absolute.
	 *
	 * @param n
	 *            The filter number
	 * @param defaultValue
	 *            the default value
	 * @return if the smoothing parameter is absolute
	 */
	public boolean getDataFilterParameterAbsolute(int n, boolean defaultValue)
	{
		if (n < getDataFiltersCount())
			return this.fitEngineSettings.getDataFilterSettings().getDataFilters(n).getParameters(0).getAbsolute();
		return defaultValue;
	}

	/**
	 * @param n
	 *            The filter number
	 * @return the smoothing window size
	 */
	public double getDataFilterParameterValue(int n)
	{
		if (n < getDataFiltersCount())
			return this.fitEngineSettings.getDataFilterSettings().getDataFilters(n).getParameters(0).getValue();
		throw new IndexOutOfBoundsException(n + " >= " + getDataFiltersCount());
	}

	/**
	 * Gets the smooth.
	 *
	 * @param n
	 *            The filter number
	 * @param defaultValue
	 *            the default value
	 * @return the smoothing window size
	 */
	public double getDataFilterParameterValue(int n, double defaultValue)
	{
		if (n < getDataFiltersCount())
			return this.fitEngineSettings.getDataFilterSettings().getDataFilters(n).getParameters(0).getValue();
		return defaultValue;
	}

	/**
	 * Sets the data filter.
	 *
	 * @param dataFilterMethod
	 *            the filter to apply to the data before identifying local maxima
	 * @param smooth
	 *            the size of the smoothing window. The actual window is calculated dynamically in conjunction with the
	 *            peak widths.
	 * @param absolute
	 *            the absolute
	 * @param n
	 *            The filter number
	 */
	public void setDataFilter(DataFilterMethod dataFilterMethod, double smooth, boolean absolute, int n)
	{
		DataFilterSettings.Builder b = fitEngineSettings.getDataFilterSettingsBuilder();
		//truncateFilters(b, n + 1);
		DataFilter.Builder b2 = (b.getDataFiltersCount() == n) ? b.addDataFiltersBuilder() : b.getDataFiltersBuilder(n);
		b2.setDataFilterMethod(dataFilterMethod);
		b2.clearParameters();
		RelativeParameter.Builder b3 = b2.addParametersBuilder();
		b3.setValue(smooth);
		b3.setAbsolute(absolute);
	}

	/**
	 * Sets the data filter without the absolute flag. If no parameter exists then the absolute flag will have the
	 * default value (false).
	 *
	 * @param dataFilterMethod
	 *            the filter to apply to the data before identifying local maxima
	 * @param smooth
	 *            the size of the smoothing window. The actual window is calculated dynamically in conjunction with the
	 *            peak widths.
	 * @param absolute
	 *            the absolute
	 * @param n
	 *            The filter number
	 */
	public void setDataFilter(DataFilterMethod dataFilterMethod, double smooth, int n)
	{
		DataFilterSettings.Builder b = fitEngineSettings.getDataFilterSettingsBuilder();
		//truncateFilters(b, n + 1);
		DataFilter.Builder b2 = (b.getDataFiltersCount() == n) ? b.addDataFiltersBuilder() : b.getDataFiltersBuilder(n);
		b2.setDataFilterMethod(dataFilterMethod);
		RelativeParameter.Builder b3 = (b2.getParametersCount() == 0) ? b2.addParametersBuilder()
				: b2.getParametersBuilder(0);
		b3.setValue(smooth);
	}

	/**
	 * Sets the data filter
	 *
	 * @param dataFilter
	 *            the data filter
	 * @param smooth
	 *            the size of the smoothing window. The actual window is calculated dynamically in conjunction with the
	 *            peak widths using {@link #getHWHMMin()}.
	 * @param absolute
	 *            Set to true to use the absolute value. Otherwise compute relative to {@link #getHWHMMin()}.
	 * @param n
	 *            The filter number
	 */
	public void setDataFilter(int dataFilter, double smooth, boolean absolute, int n)
	{
		DataFilterMethod m = DataFilterMethod.forNumber(dataFilter);
		if (m != null)
		{
			setDataFilter(m, smooth, absolute, n);
		}
	}

	/**
	 * Sets the data filter without the absolute flag. If no parameter exists then the absolute flag will have the
	 * default value (false).
	 *
	 * @param dataFilter
	 *            the data filter
	 * @param smooth
	 *            the size of the smoothing window. The actual window is calculated dynamically in conjunction with the
	 *            peak widths using {@link #getHWHMMin()}.
	 * @param n
	 *            The filter number
	 */
	public void setDataFilter(int dataFilter, double smooth, int n)
	{
		DataFilterMethod m = DataFilterMethod.forNumber(dataFilter);
		if (m != null)
		{
			setDataFilter(m, smooth, n);
		}
	}

	/**
	 * Set the number of filters to use. Call this method when all the filters have been set to clear any other stored
	 * filters from memory. This allows the user to call {@link #setDataFilter(DataFilter, double, int)} 3 times when
	 * then configuration has more than 3 filters already stored.
	 * 
	 * @param n
	 *            The number of filters
	 */
	public void setNumberOfFilters(int n)
	{
		if (n < 1)
			n = 1;
		DataFilterSettings.Builder b = fitEngineSettings.getDataFilterSettingsBuilder();
		// Truncate the filter
		truncateFilters(b, n);
	}

	/**
	 * Get the fitting width. If the fitting parameter is relative this is calculated using the maximum peak standard
	 * deviation multiplied by the fitting parameter, otherwise the absolute value is used. The value is then rounded up
	 * to the next integer and a minimum of 2 is returned.
	 * 
	 * @return The fitting width
	 */
	public int getFittingWidth()
	{
		// Region for peak fitting
		int fitting = (int) Math.ceil(convertUsingSDMax(getFittingParameter()));
		if (fitting < 2)
			fitting = 2;
		return fitting;
	}

	/**
	 * Create the spot filter for identifying candidate maxima. The actual border, search width and smoothing parameters
	 * can be configured relative to the configured standard deviations or left absolute. The standard deviation is used
	 * to determine the Half-Width at Half-Maximum (HWHM) for each dimension and the parameters set as follows.
	 * 
	 * <pre>
	 * 
	 * int search = (int) Math.ceil(getSearch() * hwhmMax);
	 * int border = (int) Math.floor(getBorder() * hwhmMax);
	 * // For each filter
	 * double smooth = getSmooth(i) * hwhmMin;
	 * </pre>
	 *
	 * @return the maxima spot filter
	 */
	public MaximaSpotFilter createSpotFilter()
	{
		// Respect the absolute parameter absolute flag for all the distances.

		// Get the half-width at half maximum
		double hwhmMin = getHWHMMin();
		double hwhmMax = getHWHMMax();

		// Note: rounding to 2 decimal places is a simple method for removing small errors
		// in floating point precision from creating an incorrect integer  

		// Region for maxima finding
		int search = (int) Math.ceil(convert(getSearchParameter(), hwhmMax, 2));
		if (search < 1)
			search = 1;

		// Border where peaks are ignored
		int border = (int) Math.floor(convert(getBorderParameter(), hwhmMax, 2));
		if (border < 0)
			border = 0;

		DataProcessor processor0 = createDataProcessor(border, 0, hwhmMin);
		DataFilterSettings f = fitEngineSettings.getDataFilterSettings();
		final int nFilters = f.getDataFiltersCount();

		final MaximaSpotFilter spotFilter;
		switch (f.getDataFilterType())
		{
			case JURY:
				if (nFilters > 1)
				{
					DataProcessor[] processors = new DataProcessor[nFilters];
					processors[0] = processor0;
					for (int i = 1; i < nFilters; i++)
						processors[i] = createDataProcessor(border, i, hwhmMin);
					spotFilter = new JurySpotFilter(search, border, processors);
					break;
				}

			case DIFFERENCE:
				if (nFilters > 1)
				{
					DataProcessor processor1 = createDataProcessor(border, 1, hwhmMin);
					spotFilter = new DifferenceSpotFilter(search, border, processor0, processor1);
					break;
				}

			case SINGLE:
			default:
				spotFilter = new SingleSpotFilter(search, border, processor0);
		}

		if (getFitConfiguration().isPerPixelCameraType())
		{
			if (!spotFilter.isWeighted())
				throw new IllegalStateException(
						"Camera type requires a weighted spot filter: " + fitConfiguration.getCameraType());
			CameraModel model = fitConfiguration.getCameraModel();
			if (model == null || !model.isPerPixelModel())
				throw new IllegalStateException("Weighted spot filter requires a per-pixel camera model");
		}

		return spotFilter;
	}

	/**
	 * Gets the minimum SD using the initial standard deviations
	 * <p>
	 * Note: This uses initial peak SD0 and optionally initial peak SD1 if width fitting is enabled for the second
	 * dimension.
	 *
	 * @return the SD min
	 */
	public double getSDMin()
	{
		double[] w = PSFHelper.getGaussian2DWxWy(fitConfiguration.psf);

		final double initialPeakStdDev0 = w[0];
		final double initialPeakStdDev1 = w[1];

		// Use 1 if zero to get at least a single pixel width
		double widthMin = (initialPeakStdDev0 > 0) ? initialPeakStdDev0 : 1;

		// Only use the second width if this is part of the function. 
		// This should be taken care of within the PSF helper.
		if (initialPeakStdDev1 > 0)
			widthMin = FastMath.min(initialPeakStdDev1, widthMin);

		return widthMin;
	}

	/**
	 * Gets the minimum HWHM using the initial standard deviations
	 * <p>
	 * Note: This uses initial peak SD0 and optionally initial peak SD1 if width fitting is enabled for the second
	 * dimension.
	 *
	 * @return the HWHM min
	 */
	public double getHWHMMin()
	{
		// Get the half-width at half maximim
		return Gaussian2DFunction.SD_TO_HWHM_FACTOR * getSDMin();
	}

	/**
	 * Gets the maximum SD using the initial standard deviations
	 * <p>
	 * Note: This uses initial peak SD0 and optionally initial peak SD1 if width fitting is enabled for the second
	 * dimension.
	 *
	 * @return the SD max
	 */
	public double getSDMax()
	{
		double[] w = PSFHelper.getGaussian2DWxWy(fitConfiguration.psf);

		final double initialPeakStdDev0 = w[0];
		final double initialPeakStdDev1 = w[1];

		// Use 1 if zero to get at least a single pixel width
		double widthMax = (initialPeakStdDev0 > 0) ? initialPeakStdDev0 : 1;

		// Only use the second width if this is part of the function.
		// This should be taken care of within the PSF helper.
		if (initialPeakStdDev1 > 0)
			widthMax = FastMath.max(initialPeakStdDev1, widthMax);

		return widthMax;
	}

	/**
	 * Gets the maximum HWHM using the initial standard deviations
	 * <p>
	 * Note: This uses initial peak SD0 and optionally initial peak SD1 if width fitting is enabled for the second
	 * dimension.
	 *
	 * @return the HWHM max
	 */
	public double getHWHMMax()
	{
		// Get the half-width at half maximim
		return Gaussian2DFunction.SD_TO_HWHM_FACTOR * getSDMax();
	}

	/**
	 * Convert the relative parameter using the scale. If the parameter is absolute then return the unscaled value.
	 *
	 * @param p
	 *            the p
	 * @param scale
	 *            the scale
	 * @param decimalPlaces
	 *            the decimal places to round the scaled number (set to negative to ignore)
	 * @return the double
	 */
	public static double convert(RelativeParameter p, double scale, int decimalPlaces)
	{
		return (p.getAbsolute()) ? p.getValue() : round(p.getValue() * scale, decimalPlaces);
	}

	/**
	 * Round to the given number of decimal places.
	 *
	 * @param value
	 *            the value
	 * @param decimalPlaces
	 *            the decimal places
	 * @return the double
	 */
	public static double round(double value, int decimalPlaces)
	{
		return (decimalPlaces >= 0) ? Maths.roundUsingDecimalPlaces(value, decimalPlaces) : value;
	}

	/**
	 * Convert using half-width half-max as the scale.
	 *
	 * @param p
	 *            the p
	 * @return the converted value
	 */
	public double convertUsingHWHMax(RelativeParameter p)
	{
		return (p.getAbsolute()) ? p.getValue() : p.getValue() * getHWHMMax();
	}

	/**
	 * Convert using half-width half-max as the scale.
	 *
	 * @param p
	 *            the p
	 * @return the converted value
	 */
	public double convertUsingHWHMin(RelativeParameter p)
	{
		return (p.getAbsolute()) ? p.getValue() : p.getValue() * getHWHMMin();
	}

	/**
	 * Convert using SD max as the scale.
	 *
	 * @param p
	 *            the p
	 * @return the converted value
	 */
	public double convertUsingSDMax(RelativeParameter p)
	{
		return (p.getAbsolute()) ? p.getValue() : p.getValue() * getSDMax();
	}

	/**
	 * Gets the number of filters for the configured filter type.
	 *
	 * @return the number of filters
	 */
	public int getNumberOfFilters()
	{
		DataFilterSettings f = fitEngineSettings.getDataFilterSettings();
		final int nFilters = f.getDataFiltersCount();
		switch (f.getDataFilterType())
		{
			case JURY:
				if (nFilters > 1)
				{
					return nFilters;
				}

			case DIFFERENCE:
				if (nFilters > 1)
				{
					return 2;
				}

			case SINGLE:
			default:
				return 1;
		}
	}

	private double getSmoothingWindow(DataFilter f, double hwhmMin)
	{
		//return BlockAverageDataProcessor.convert(smoothingParameter * hwhmMin);
		RelativeParameter rp = f.getParameters(0);

		// Q. Why is this rounded. Is it just to make a nicer number?
		return convert(rp, hwhmMin, 2);
	}

	private DataProcessor createDataProcessor(int border, int n, double hwhm)
	{
		if (n < this.fitEngineSettings.getDataFilterSettings().getDataFiltersCount())
		{
			DataFilter f = fitEngineSettings.getDataFilterSettings().getDataFilters(n);
			return createDataProcessor(border, f.getDataFilterMethod(), getSmoothingWindow(f, hwhm));
		}
		return null;
	}

	/**
	 * Create a data processor for the spot filter
	 * 
	 * @param border
	 * @param dataFilter
	 * @param parameter
	 * @return the data processor
	 */
	public static DataProcessor createDataProcessor(int border, DataFilterMethod dataFilter, double parameter)
	{
		switch (dataFilter)
		{
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
				throw new RuntimeException("Not yet implemented: " + dataFilter.toString());
		}
	}

	/**
	 * Copy data filter.
	 *
	 * @param config
	 *            the config
	 */
	public void copyDataFilter(FitEngineConfiguration config)
	{
		fitEngineSettings.setDataFilterSettings(config.fitEngineSettings.getDataFilterSettings());
	}

	/**
	 * Configure the output units from fitting using the current calibration and fit solver settings.
	 * <p>
	 * This method should be called before the calibration is passed to any object that will handle the fitting output.
	 */
	public void configureOutputUnits()
	{
		FitConfiguration fitConfig = getFitConfiguration();
		// If there is no calibration then the writer will just have the defaults
		CalibrationWriter calibration = fitConfig.getCalibrationWriter();

		// Fitting is always done pixels and radians
		calibration.setDistanceUnit(DistanceUnit.PIXEL);
		calibration.setAngleUnit(AngleUnit.RADIAN);

		// Most fitters fit in photons unless we have no calibration.
		IntensityUnit intensityUnit = IntensityUnit.PHOTON;

		if (//calibration.getCountPerPhoton() == 0 || 
			fitConfig.isFitCameraCounts())
			intensityUnit = IntensityUnit.COUNT;

		calibration.setIntensityUnit(intensityUnit);
		
		// This initialises the calibration precision method
		fitConfig.getFilterPrecisionMethod();
	}
}