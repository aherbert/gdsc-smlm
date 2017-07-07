package gdsc.smlm.engine;

import org.apache.commons.math3.util.FastMath;

import gdsc.core.utils.Maths;
import gdsc.smlm.data.config.CalibrationConfig.Calibration;
import gdsc.smlm.data.config.CalibrationConfigHelper;
import gdsc.smlm.data.config.CalibrationWriter;
import gdsc.smlm.data.config.FitConfig.DataFilter;
import gdsc.smlm.data.config.FitConfig.DataFilterMethod;
import gdsc.smlm.data.config.FitConfig.DataFilterSettings;
import gdsc.smlm.data.config.FitConfig.DataFilterType;
import gdsc.smlm.data.config.FitConfig.FitEngineSettings;
import gdsc.smlm.data.config.FitConfig.NoiseEstimatorMethod;
import gdsc.smlm.data.config.FitConfig.RelativeParameter;
import gdsc.smlm.data.config.FitConfigHelper;
import gdsc.smlm.data.config.PSFConfig.PSF;
import gdsc.smlm.data.config.PSFConfigHelper;
import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.data.config.UnitConfig.AngleUnit;
import gdsc.smlm.data.config.UnitConfig.DistanceUnit;
import gdsc.smlm.data.config.UnitConfig.IntensityUnit;
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
		this(FitConfigHelper.defaultFitEngineSettings, CalibrationConfigHelper.defaultCalibration,
				PSFConfigHelper.defaultOneAxisGaussian2DPSF);
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
	 * @return the size of the region to search for local maxima
	 */
	public double getSearch()
	{
		return fitEngineSettings.getSearch().getValue();
	}

	/**
	 * @param search
	 *            the size of the region to search for local maxima. The actual window is calculated dynamically
	 *            in conjunction with the peak widths.
	 */
	public void setSearch(double search)
	{
		fitEngineSettings.getSearchBuilder().setValue(search);
	}

	/**
	 * @return the border
	 *         the size of the border region to ignore. The actual window is calculated dynamically
	 *         in conjunction with the peak widths.
	 */
	public double getBorder()
	{
		return fitEngineSettings.getBorder().getValue();
	}

	/**
	 * @param border
	 *            the size of the border region to ignore
	 */
	public void setBorder(double border)
	{
		fitEngineSettings.getBorderBuilder().setValue(border);
	}

	/**
	 * @return the fitting window size
	 */
	public double getFitting()
	{
		return fitEngineSettings.getFitting().getValue();
	}

	/**
	 * @param fitting
	 *            the size of the window used for fitting around a maxima. The actual window is calculated dynamically
	 *            in conjunction with the peak widths.
	 */
	public void setFitting(double fitting)
	{
		fitEngineSettings.getFittingBuilder().setValue(fitting);
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
		fitEngineSettings.setFailuresLimit(failuresLimit);
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
	 * @return The distance within which spots are considered duplicates
	 */
	public double getDuplicateDistance()
	{
		return fitEngineSettings.getDuplicateDistance().getValue();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public FitEngineConfiguration clone()
	{
		// This is not a complete duplicate. The settings builder objects with the 
		// underlying configuration will be the same between all instances. 
		try
		{
			FitEngineConfiguration f = (FitEngineConfiguration) super.clone();
			// Ensure the object is duplicated and not passed by reference.
			if (fitConfiguration != null)
				f.fitConfiguration = fitConfiguration.clone();
			return f;
		}
		catch (CloneNotSupportedException e)
		{
			// Ignore
		}
		return null;

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
		while (b.getDataFilterCount() > n)
			b.removeDataFilter(b.getDataFilterCount() - 1);
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
	//		if (n < this.fitEngineSettings.getDataFilterSettings().getDataFilterCount())
	//			return this.fitEngineSettings.getDataFilterSettings().getDataFilter(n);
	//		return null;
	//	}

	/**
	 * Gets the data filter count.
	 *
	 * @return the data filter count
	 */
	public int getDataFilterCount()
	{
		return this.fitEngineSettings.getDataFilterSettings().getDataFilterCount();
	}

	/**
	 * @param n
	 *            The filter number
	 * @return the filter to apply to the data before identifying local maxima
	 */
	public DataFilterMethod getDataFilterMethod(int n)
	{
		if (n < getDataFilterCount())
			return this.fitEngineSettings.getDataFilterSettings().getDataFilter(n).getDataFilterMethod();
		throw new IndexOutOfBoundsException(n + " >= " + getDataFilterCount());
	}

	/**
	 * @param n
	 *            The filter number
	 * @return if the smoothing parameter is absolute
	 */
	public boolean getDataFilterAbsolute(int n)
	{
		if (n < getDataFilterCount())
			return this.fitEngineSettings.getDataFilterSettings().getDataFilter(n).getParameter(0).getAbsolute();
		throw new IndexOutOfBoundsException(n + " >= " + getDataFilterCount());
	}

	/**
	 * @param n
	 *            The filter number
	 * @return the smoothing window size
	 */
	public double getSmooth(int n)
	{
		if (n < getDataFilterCount())
			return this.fitEngineSettings.getDataFilterSettings().getDataFilter(n).getParameter(0).getValue();
		throw new IndexOutOfBoundsException(n + " >= " + getDataFilterCount());
	}

	/**
	 * @param n
	 *            The filter number
	 * @return the absolute flag for the smoothing window size
	 */
	public boolean getAbsolute(int n)
	{
		if (n < getDataFilterCount())
			return this.fitEngineSettings.getDataFilterSettings().getDataFilter(n).getParameter(0).getAbsolute();
		throw new IndexOutOfBoundsException(n + " >= " + getDataFilterCount());
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
		truncateFilters(b, n + 1);
		DataFilter.Builder b2 = (b.getDataFilterCount() == n) ? b.addDataFilterBuilder() : b.getDataFilterBuilder(n);
		b2.setDataFilterMethod(dataFilterMethod);
		b2.clearParameter();
		RelativeParameter.Builder b3 = b2.addParameterBuilder();
		b3.setValue(smooth);
		b3.setAbsolute(absolute);
	}

	/**
	 * Sets the data filter.
	 *
	 * @param dataFilter
	 *            the data filter
	 * @param smooth
	 *            the size of the smoothing window. The actual window is calculated dynamically in conjunction with the
	 *            peak widths.
	 * @param absolute
	 *            the absolute
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
	 * Get the relative fitting width. This is calculated using the maximum peak standard deviation multiplied by the
	 * fitting parameter.
	 * 
	 * @return The fitting width
	 */
	public int getRelativeFitting()
	{
		double[] w = PSFHelper.getGaussian2DWxWy(getFitConfiguration().getPSF());

		final double initialPeakStdDev0 = w[0];
		final double initialPeakStdDev1 = w[1];

		double widthMax = (initialPeakStdDev0 > 0) ? initialPeakStdDev0 : 1;
		if (initialPeakStdDev1 > 0)
			widthMax = FastMath.max(initialPeakStdDev1, widthMax);

		// Region for peak fitting
		int fitting = (int) Math.ceil(getFitting() * widthMax);
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
	 * 
	 * </pre>
	 *
	 * @param relative
	 *            True if the parameters should be made relative to the configured standard deviations
	 * @return the maxima spot filter
	 * @deprecated The relative flag will be removed and the configured parameters must be set as relative/absolute
	 */
	@Deprecated
	public MaximaSpotFilter createSpotFilter(boolean relative)
	{
		final double hwhmMin, hwhmMax;

		if (relative)
		{
			// Get the half-width at half maximim
			hwhmMin = getHWHMMin();
			hwhmMax = getHWHMMax();
		}
		else
		{
			hwhmMin = hwhmMax = 1;
		}

		// Region for maxima finding
		int search = (int) Math.ceil(Maths.round(getSearch() * hwhmMax, 0.01));
		if (search < 1)
			search = 1;

		// Border where peaks are ignored
		int border = (int) Math.floor(Maths.round(getBorder() * hwhmMax, 0.01));
		if (border < 0)
			border = 0;

		DataProcessor processor0 = createDataProcessor(border, 0, hwhmMin);
		DataFilterSettings f = fitEngineSettings.getDataFilterSettings();
		final int nFilters = f.getDataFilterCount();

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

		// Note: It is possible to configure the score data processor here. However small tests 
		// show this often reduces performance and the additional parameters make it harder to 
		// configure. It is a subject for future work.

		return spotFilter;
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
		double[] w = PSFHelper.getGaussian2DWxWy(fitConfiguration.psf);

		final double initialPeakStdDev0 = w[0];
		final double initialPeakStdDev1 = w[1];

		// Use 1 if zero to get at least a single pixel width
		double widthMin = (initialPeakStdDev0 > 0) ? initialPeakStdDev0 : 1;

		// Only use the second width if this is part of the function. 
		// This should be taken care of within the PSF helper.
		if (initialPeakStdDev1 > 0)
			widthMin = FastMath.min(initialPeakStdDev1, widthMin);

		// Get the half-width at half maximim
		return Gaussian2DFunction.SD_TO_HWHM_FACTOR * widthMin;
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
		double[] w = PSFHelper.getGaussian2DWxWy(fitConfiguration.psf);

		final double initialPeakStdDev0 = w[0];
		final double initialPeakStdDev1 = w[1];

		// Use 1 if zero to get at least a single pixel width
		double widthMax = (initialPeakStdDev0 > 0) ? initialPeakStdDev0 : 1;

		// Only use the second width if this is part of the function.
		// This should be taken care of within the PSF helper.
		if (initialPeakStdDev1 > 0)
			widthMax = FastMath.max(initialPeakStdDev1, widthMax);

		// Get the half-width at half maximim
		return Gaussian2DFunction.SD_TO_HWHM_FACTOR * widthMax;
	}

	/**
	 * Gets the number of filters for the configured filter type.
	 *
	 * @return the number of filters
	 */
	public int getNumberOfFilters()
	{
		DataFilterSettings f = fitEngineSettings.getDataFilterSettings();
		final int nFilters = f.getDataFilterCount();
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
		RelativeParameter rp = f.getParameter(0);
		double p = rp.getValue();
		if (!rp.getAbsolute())
			p *= hwhmMin;
		return Maths.round(p, 0.01);
	}

	private DataProcessor createDataProcessor(int border, int n, double hwhm)
	{
		if (n < this.fitEngineSettings.getDataFilterSettings().getDataFilterCount())
		{
			DataFilter f = fitEngineSettings.getDataFilterSettings().getDataFilter(n);
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
		getFitConfiguration();

		// If there is no calibration then the writer will just have the defaults
		CalibrationWriter calibration = fitConfiguration.getCalibrationWriter();

		// Fitting is always done pixels and radians
		calibration.setDistanceUnit(DistanceUnit.PIXEL);
		calibration.setAngleUnit(AngleUnit.RADIAN);

		// Most fitters fit in photons unless we have no calibration.
		IntensityUnit intensityUnit = IntensityUnit.PHOTON;

		if (calibration.getGain() == 0)
			intensityUnit = IntensityUnit.COUNT;

		calibration.setIntensityUnit(intensityUnit);
	}
}