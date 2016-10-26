package gdsc.smlm.engine;

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

import java.util.Arrays;

import org.apache.commons.math3.util.FastMath;

import gdsc.core.utils.Maths;
import gdsc.core.utils.NoiseEstimator;
import gdsc.core.utils.NoiseEstimator.Method;
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
import gdsc.smlm.fitting.FitConfiguration;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;

/**
 * Specifies the configuration for the fit engine
 */
public class FitEngineConfiguration implements Cloneable
{
	private FitConfiguration fitConfiguration;

	// Analysis* shows the best area-under-precision-recall curve (AUC) using a mean filter or
	// a Gaussian filter with ~1.2 SD smoothing. The Gaussian filter is more robust to width mismatch but
	// the mean filter will be faster as it uses a smaller block size. The Gaussian filter has higher 
	// recall but lower precision as it identifies more spots due to the shape of the smoothing filter.
	// The overall AUC is very similar.
	//
	// Note: Setting the parameter at a higher level allows the filter to work on out-of-focus spots which
	// will have a wider PSF.
	//
	// *Analysis was performed on simulated data using a Image PSF with spots of 20-100 photons at a 
	// depth of up to 1380nm (the PSF limit).

	private double search = 1;
	private double border = 1;
	private double fitting = 3;
	private int failuresLimit = 3;
	private boolean includeNeighbours = true;
	private double neighbourHeightThreshold = 0.3;
	private double residualsThreshold = 1;
	private NoiseEstimator.Method noiseMethod = Method.QUICK_RESIDUALS_LEAST_MEAN_OF_SQUARES;
	private DataFilterType dataFilterType = DataFilterType.SINGLE;
	private double[] smooth = new double[] { 1.2 };
	private DataFilter[] dataFilter = new DataFilter[] { DataFilter.MEAN };

	/**
	 * Constructor
	 * 
	 * @param fitConfiguration
	 */
	public FitEngineConfiguration(FitConfiguration fitConfiguration)
	{
		this.fitConfiguration = fitConfiguration;
	}

	/**
	 * @return the size of the region to search for local maxima
	 */
	public double getSearch()
	{
		return search;
	}

	/**
	 * @param search
	 *            the size of the region to search for local maxima. The actual window is calculated dynamically
	 *            in conjunction with the peak widths.
	 */
	public void setSearch(double search)
	{
		this.search = search;
	}

	/**
	 * @return the border
	 *         the size of the border region to ignore. The actual window is calculated dynamically
	 *         in conjunction with the peak widths.
	 */
	public double getBorder()
	{
		return border;
	}

	/**
	 * @param border
	 *            the size of the border region to ignore
	 */
	public void setBorder(double border)
	{
		this.border = border;
	}

	/**
	 * @return the fitting window size
	 */
	public double getFitting()
	{
		return fitting;
	}

	/**
	 * @param fitting
	 *            the size of the window used for fitting around a maxima. The actual window is calculated dynamically
	 *            in conjunction with the peak widths.
	 */
	public void setFitting(double fitting)
	{
		this.fitting = fitting;
	}

	/**
	 * Set the failures limit. When failures exceeds the failures limit then stop fitting.
	 * 
	 * @return the failuresLimit
	 */
	public int getFailuresLimit()
	{
		return failuresLimit;
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
		this.failuresLimit = failuresLimit;
	}

	/**
	 * @return the fitConfiguration
	 */
	public FitConfiguration getFitConfiguration()
	{
		return fitConfiguration;
	}

	/**
	 * @return the includeNeighbours
	 */
	public boolean isIncludeNeighbours()
	{
		return includeNeighbours;
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
		this.includeNeighbours = includeNeighbours;
	}

	/**
	 * @return the neighbourHeightThreshold
	 */
	public double getNeighbourHeightThreshold()
	{
		return neighbourHeightThreshold;
	}

	/**
	 * @param neighbourHeightThreshold
	 *            Set the height threshold that determines if a neighbour peak should be fitted (specified as a fraction
	 *            of the central peak relative to the background)
	 */
	public void setNeighbourHeightThreshold(double neighbourHeightThreshold)
	{
		this.neighbourHeightThreshold = neighbourHeightThreshold;
	}

	/**
	 * @return the residuals threshold
	 */
	public double getResidualsThreshold()
	{
		return residualsThreshold;
	}

	/**
	 * @param residualsThreshold
	 *            Set the threshold for the residuals analysis that determines if a two-kernel model should be fitted
	 */
	public void setResidualsThreshold(double residualsThreshold)
	{
		this.residualsThreshold = residualsThreshold;
	}

	/**
	 * @return the method used to estimate the image noise
	 */
	public NoiseEstimator.Method getNoiseMethod()
	{
		return noiseMethod;
	}

	/**
	 * @param noiseMethod
	 *            Set the method used to estimate the image noise
	 */
	public void setNoiseMethod(NoiseEstimator.Method noiseMethod)
	{
		this.noiseMethod = noiseMethod;
	}

	/**
	 * @param noiseMethod
	 *            Set the method used to estimate the image noise
	 */
	public void setNoiseMethod(int noiseMethod)
	{
		if (noiseMethod >= 0 && noiseMethod < NoiseEstimator.Method.values().length)
		{
			setNoiseMethod(NoiseEstimator.Method.values()[noiseMethod]);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public FitEngineConfiguration clone()
	{
		try
		{
			FitEngineConfiguration f = (FitEngineConfiguration) super.clone();
			// Ensure the object is duplicated and not passed by reference.
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
	 * Ensure that the internal state of the object is initialised. This is used after deserialisation since some state
	 * is not saved but restored from other property values.
	 */
	public void initialiseState()
	{
		if (fitConfiguration == null)
			fitConfiguration = new FitConfiguration();
		fitConfiguration.initialiseState();
		if (noiseMethod == null)
			noiseMethod = Method.QUICK_RESIDUALS_LEAST_MEAN_OF_SQUARES;
		if (dataFilter == null || smooth == null)
			setDataFilter(DataFilter.MEAN, 1.3, 0);
		// Do this last as it resizes the dataFilter and smooth arrays
		if (dataFilterType == null)
			setDataFilterType(DataFilterType.SINGLE);
	}

	/**
	 * @return the type of filter to apply to the data before identifying local maxima
	 */
	public DataFilterType getDataFilterType()
	{
		return dataFilterType;
	}

	/**
	 * @param DataFilterType
	 *            the type of filter to apply to the data before identifying local maxima
	 */
	public void setDataFilterType(DataFilterType dataFilterType)
	{
		this.dataFilterType = dataFilterType;

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
		resizeFilters(n);
	}

	private void resizeFilters(int n)
	{
		if (this.dataFilter == null || this.dataFilter.length < n)
		{
			this.dataFilter = Arrays.copyOf(this.dataFilter, n);
			this.smooth = Arrays.copyOf(this.smooth, n);
		}
	}

	/**
	 * @param DataFilterType
	 *            the type of filter to apply to the data before identifying local maxima
	 */
	public void setDataFilterType(int dataFilterType)
	{
		if (dataFilterType >= 0 && dataFilterType < DataFilterType.values().length)
		{
			setDataFilterType(DataFilterType.values()[dataFilterType]);
		}
	}

	/**
	 * @param n
	 *            The filter number
	 * @return the filter to apply to the data before identifying local maxima
	 */
	public DataFilter getDataFilter(int n)
	{
		if (n < this.dataFilter.length)
			return dataFilter[n];
		return DataFilter.MEAN;
	}

	/**
	 * @param n
	 *            The filter number
	 * @return the smoothing window size
	 */
	public double getSmooth(int n)
	{
		if (n < this.smooth.length)
			return smooth[n];
		return 0;
	}

	/**
	 * @param DataFilter
	 *            the filter to apply to the data before identifying local maxima
	 * @param smooth
	 *            the size of the smoothing window. The actual window is calculated dynamically in conjunction with the
	 *            peak widths.
	 * @param n
	 *            The filter number
	 */
	public void setDataFilter(DataFilter dataFilter, double smooth, int n)
	{
		resizeFilters(n + 1);
		this.dataFilter[n] = dataFilter;
		this.smooth[n] = smooth;
	}

	/**
	 * @param DataFilter
	 *            the filter to apply to the data before identifying local maxima
	 * @param smooth
	 *            the size of the smoothing window. The actual window is calculated dynamically in conjunction with the
	 *            peak widths.
	 * @param n
	 *            The filter number
	 */
	public void setDataFilter(int dataFilter, double smooth, int n)
	{
		if (dataFilter >= 0 && dataFilter < DataFilter.values().length)
		{
			setDataFilter(DataFilter.values()[dataFilter], smooth, n);
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
		// Truncate the filter
		if (this.dataFilter != null)
		{
			this.dataFilter = Arrays.copyOf(this.dataFilter, n);
			this.smooth = Arrays.copyOf(this.smooth, n);
		}
	}

	/**
	 * Get the relative fitting width. This is calculated using the maximum peak standard deviation multiplied by the
	 * fitting parameter.
	 * 
	 * @return The fitting width
	 */
	public int getRelativeFitting()
	{
		final double initialPeakStdDev0 = fitConfiguration.getInitialPeakStdDev0();
		final double initialPeakStdDev1 = fitConfiguration.getInitialPeakStdDev1();

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
	 * @return
	 */
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
		final int nFilters = Math.min(dataFilter.length, smooth.length);

		final MaximaSpotFilter spotFilter;
		switch (dataFilterType)
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
		final FitConfiguration fitConfiguration = getFitConfiguration();
		final double initialPeakStdDev0 = fitConfiguration.getInitialPeakStdDev0();
		// Use 1 if zero to get at least a single pixel width
		double widthMin = (initialPeakStdDev0 > 0) ? initialPeakStdDev0 : 1;

		// Only use the second width if this is part of the function
		if (fitConfiguration.isWidth1Fitting())
		{
			final double initialPeakStdDev1 = fitConfiguration.getInitialPeakStdDev1();
			if (initialPeakStdDev1 > 0)
				widthMin = FastMath.min(initialPeakStdDev1, widthMin);
		}

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
		final FitConfiguration fitConfiguration = getFitConfiguration();
		final double initialPeakStdDev0 = fitConfiguration.getInitialPeakStdDev0();
		// Use 1 if zero to get at least a single pixel width
		double widthMax = (initialPeakStdDev0 > 0) ? initialPeakStdDev0 : 1;

		// Only use the second width if this is part of the function
		if (fitConfiguration.isWidth1Fitting())
		{
			final double initialPeakStdDev1 = fitConfiguration.getInitialPeakStdDev1();
			if (initialPeakStdDev1 > 0)
				widthMax = FastMath.max(initialPeakStdDev1, widthMax);
		}

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
		final int nFilters = Math.min(dataFilter.length, smooth.length);
		switch (dataFilterType)
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

	private double getSmoothingWindow(double smoothingParameter, double hwhmMin)
	{
		//return BlockAverageDataProcessor.convert(smoothingParameter * hwhmMin);
		return Maths.round(smoothingParameter * hwhmMin, 0.01);
	}

	private DataProcessor createDataProcessor(int border, int n, double hwhm)
	{
		if (n < dataFilter.length && n < smooth.length)
			return createDataProcessor(border, getDataFilter(n), getSmoothingWindow(getSmooth(n), hwhm));
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
	public static DataProcessor createDataProcessor(int border, DataFilter dataFilter, double parameter)
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

	public void copyDataFilter(FitEngineConfiguration config)
	{
		setDataFilterType(config.getDataFilterType());
		final int nFilters = config.getNumberOfFilters();
		for (int n = 0; n < nFilters; n++)
		{
			setDataFilter(config.getDataFilter(n), config.getSmooth(n), n);
		}
	}
}