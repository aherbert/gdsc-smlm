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

import gdsc.smlm.engine.FitParameters.FitTask;
import gdsc.smlm.engine.filter.DistanceResultFilter;
import gdsc.smlm.engine.filter.ResultFilter;
import gdsc.smlm.filters.AverageFilter;
import gdsc.smlm.filters.NonMaximumSuppression;
import gdsc.smlm.fitting.FitConfiguration;
import gdsc.smlm.fitting.FitResult;
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.Gaussian2DFitter;
import gdsc.smlm.fitting.function.Gaussian2DFunction;
import gdsc.smlm.fitting.function.GaussianFunction;
import gdsc.smlm.fitting.logging.Logger;
import gdsc.smlm.results.ExtendedPeakResult;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResults;
import gdsc.smlm.utils.ImageExtractor;
import gdsc.smlm.utils.NoiseEstimator;
import gdsc.smlm.utils.Sort;

import java.awt.Rectangle;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.concurrent.BlockingQueue;

/**
 * Fits local maxima using a 2D Gaussian.
 */
public class FitWorker implements Runnable
{
	/**
	 * Testings on a idealised dataset of simulated data show that multiple peaks increases the iterations
	 * but the increase asymptotes. Initial rate is 2-fold for 7x7 region, decreasing to 1.5-fold for 17x17 region.
	 * Best solution is to add a set of iterations for each additional peak.
	 */
	private static final int ITERATION_INCREASE_FOR_MULTIPLE_PEAKS = 1; // 0 for no effect
	/**
	 * Testings on a idealised dataset of simulated data show that fitting the doublets increases the iterations by
	 * approx 3.5-fold for 7x7 region, decreasing to 2.5-fold for 17x17 region.
	 * An additional doubling of iterations were used when the fit of the doublet resulted in one peak being eliminated
	 * for moving.
	 */
	private static final int ITERATION_INCREASE_FOR_DOUBLETS = 4; // 1 for no effect

	private Logger logger = null;
	private long time = 0;

	private double smooth = 0;
	private double smooth2 = 0;
	private int border = 0;
	private int search;
	private int fitting = 1;

	// Used for fitting
	private FitEngineConfiguration config;
	private FitConfiguration fitConfig;

	private PeakResults results;
	private BlockingQueue<FitJob> jobs;
	private AverageFilter filter = new AverageFilter();
	private Gaussian2DFitter gf;

	// Used for fitting methods
	private LinkedList<PeakResult> sliceResults;
	private double fittedBackground;
	private int slice;
	private int endT;
	private Rectangle bounds;
	private int borderLimitX, borderLimitY;
	private float[] data;
	private boolean differenceOfSmoothing;
	private float noise, background;
	private final boolean calculateNoise;
	// The coordinate of the maxima	
	private int[] xcoord = null, ycoord = null;
	// Contains the index in the list of maxima for any neighbours
	private int neighbourCount = 0;
	private int[] neighbourIndices = null;
	// Contains the index in the list of fitted results for any neighbours 
	private int fittedNeighbourCount = 0;
	private int[] fittedNeighbourIndices = null;
	private float duplicateDistance2;

	private volatile boolean finished = false;

	static int ID = 0;
	int id;

	public FitWorker(FitEngineConfiguration config, PeakResults results, BlockingQueue<FitJob> jobs)
	{
		this.config = config;
		this.fitConfig = config.getFitConfiguration();
		this.results = results;
		this.jobs = jobs;
		this.logger = fitConfig.getLog();
		gf = new Gaussian2DFitter(fitConfig);
		duplicateDistance2 = fitConfig.getDuplicateDistance() * fitConfig.getDuplicateDistance();
		calculateNoise = config.getFitConfiguration().getNoise() <= 0;
		if (!calculateNoise)
		{
			noise = config.getFitConfiguration().getNoise();
		}

		id = ID;
		ID += 1;
	}

	/**
	 * Set the parameters for smoothing the image, searching for maxima and fitting maxima
	 * 
	 * @param smooth
	 *            The first smoothing parameter
	 * @param smooth2
	 *            The second smoothing parameter (used for difference smoothing)
	 * @param border
	 *            The border of the image to ignore for fitting
	 * @param search
	 *            The block size to be used for searching for maxima
	 * @param fitting
	 *            The block size to be used for fitting
	 */
	void setSearchParameters(double smooth, double smooth2, int border, int search, int fitting)
	{
		this.smooth = smooth;
		this.smooth2 = smooth2;
		this.border = border;
		this.search = search;
		this.fitting = fitting;
	}

	/**
	 * Locate all the peaks in the image specified by the fit job
	 * 
	 * @param job
	 *            The fit job
	 */
	public void run(FitJob job)
	{
		long start = System.nanoTime();
		job.start();
		this.slice = job.slice;

		// Used for debugging
		//if (logger == null) logger = new gdsc.fitting.logging.ConsoleLogger();

		// Crop to the ROI
		bounds = job.bounds;
		int width = bounds.width;
		int height = bounds.height;
		borderLimitX = width - border;
		borderLimitY = height - border;
		data = job.data;

		FitParameters params = job.getFitParameters();
		this.endT = (params != null) ? params.endT : -1;

		// TODO: Test if using a true Gaussian filter would be faster.
		// Smoothing:
		// Either use a single smoothing pass or do a difference-of-smoothing
		float[] smoothData = null;
		differenceOfSmoothing = false;
		if (smooth > 0)
		{
			// Single smoothing (Box filter)
			smoothData = applySmoothing(data, width, height, smooth);
			//new ij.ImagePlus("smooth", new ij.process.FloatProcessor(width, height, Arrays.copyOf(smoothData,
			//		smoothData.length), null)).show();
		}
		if (smooth2 > smooth && smooth2 > 0)
		{
			if (smoothData == null)
				smoothData = Arrays.copyOf(data, width * height);

			// Difference-of-smoothing (Top-Hat Box filter)
			// See: E. J. Breen, G. H. Joss, and K. L. Williams, “Locating objects of interest
			// within biological images: The top hat box ﬁlter,” J. Comput.-Assist.
			// Microsc., vol. 3, no. 2, pp. 97–102, 1991.
			float[] smooth2Data = applySmoothing(data, width, height, smooth2);
			//new ij.ImagePlus("smooth2", new ij.process.FloatProcessor(width, height, smooth2Data, null)).show();
			for (int i = 0; i < smoothData.length; i++)
			{
				smoothData[i] -= smooth2Data[i];
			}
			differenceOfSmoothing = true;
			//new ij.ImagePlus("smoothDiff", new ij.process.FloatProcessor(width, height, smoothData, null)).show();
		}
		if (smoothData == null)
			smoothData = data;

		int[] maxIndices;
		if (params != null && params.maxIndices != null)
		{
			maxIndices = params.maxIndices;
		}
		else
		{
			// Non-maximum suppression
			maxIndices = getMaxima(smoothData, width, height);
		}

		if (logger != null)
			logger.info("%d: Slice %d: %d candidates", id, slice, maxIndices.length);

		sliceResults = new LinkedList<PeakResult>();
		job.setResults(sliceResults);
		job.setIndices(maxIndices);

		if (maxIndices.length == 0)
		{
			finish(job, start);
			return;
		}

		// Sort the maxima
		if (maxIndices.length > 1)
			// sortIgnoreEqual is numerically unstable due to the implementation of the standard libraries
			Sort.sort(maxIndices, smoothData);

		fittedBackground = 0;

		allocateArrays(maxIndices.length);

		// Compute the x/y locations of each maxima
		for (int n = 0; n < maxIndices.length; n++)
		{
			ycoord[n] = maxIndices[n] / width;
			xcoord[n] = maxIndices[n] % width;
		}

		// Always get the noise and store it with the results.
		boolean updatedNoise = false;
		if (params != null && !Float.isNaN(params.noise))
		{
			noise = params.noise;
			updatedNoise = true;
		}
		else if (calculateNoise)
		{
			noise = estimateNoise(data, width, height, config.getNoiseMethod());
			updatedNoise = true;
		}
		if (updatedNoise)
			fitConfig.setNoise(noise);

		// Get the background
		if (params != null && !Float.isNaN(params.background))
		{
			background = params.background;
		}
		else
		{
			background = estimateBackground(width, height);
		}

		//System.out.printf("Slice %d : Noise = %g\n", slice, noise);
		if (logger != null)
			logger.info("Slice %d: Noise = %f : Background = %f", slice, noise, background);

		ResultFilter filter = null;

		if (params != null && params.fitTask == FitTask.MAXIMA_IDENITIFICATION)
		{
			final float sd0 = fitConfig.getInitialPeakStdDev0();
			final float sd1 = fitConfig.getInitialPeakStdDev1();
			ImageExtractor ie = new ImageExtractor(data, width, height);
			float[] region = null;
			for (int n = 0; n < maxIndices.length; n++)
			{
				// Find the background using the perimeter of the data.
				// TODO - Perhaps the Gaussian Fitter should be used to produce the initial estimates but no actual fit done.
				// This would produce coords using the centre-of-mass.
				final int x = xcoord[n];
				final int y = ycoord[n];
				Rectangle regionBounds = ie.getBoxRegionBounds(x, y, fitting);
				region = ie.crop(regionBounds, region);
				float b = Gaussian2DFitter.getBackground(region, regionBounds.width, regionBounds.height,
						new float[] { data[maxIndices[n]] });

				// Offset the coords to the centre of the pixel. Note the bounds will be added later.
				// Subtract the background to get the amplitude estimate.
				float amplitude = smoothData[maxIndices[n]] - ((differenceOfSmoothing) ? 0 : b);
				float[] peakParams = new float[] { b, amplitude, 0, xcoord[n] + 0.5f, ycoord[n] + 0.5f, sd0, sd1 };
				sliceResults.add(createResult(bounds.x + xcoord[n], bounds.y + ycoord[n], data[maxIndices[n]], 0,
						noise, peakParams, null));
			}
		}
		else
		{
			// Perform the Gaussian fit
			int success = 0;

			// Track failures
			int[] failed = new int[maxIndices.length];
			int failedCount = 0;

			// Allow the results to be filtered for certain peaks
			if (params != null && params.filter != null)
			{
				// TODO - Check if this works.
				filter = new DistanceResultFilter(params.filter, params.distanceThreshold, maxIndices.length);
				//filter = new OptimumDistanceResultFilter(params.filter, params.distanceThreshold, maxIndices.length);
			}

			// Extract each peak and fit individually until N consecutive failures
			ImageExtractor ie = new ImageExtractor(data, width, height);
			float[] region = null;
			for (int n = 0, failures = 0; n < maxIndices.length && !finished; n++)
			{
				failures++;

				final int x = xcoord[n];
				final int y = ycoord[n];

				Rectangle regionBounds = ie.getBoxRegionBounds(x, y, fitting);
				region = ie.crop(regionBounds, region);

				if (logger != null)
					logger.info("Fitting %d:%d [%d,%d]", slice, n + 1, x + bounds.x, y + bounds.y);

				FitResult fitResult = fit(gf, region, regionBounds, maxIndices, n, smoothData);
				job.setFitResult(n, fitResult);

				// Check fit result
				if (fitResult.getStatus() == FitStatus.OK)
				{
					float[] peakParams = fitResult.getParameters();
					convertParameters(peakParams);

					float[] peakParamsDev = null;
					if (fitConfig.isComputeDeviations())
					{
						peakParamsDev = fitResult.getParameterStdDev();
						convertParameters(peakParamsDev);
					}

					int npeaks = addResults(sliceResults, n, x, y, bounds, regionBounds, peakParams, peakParamsDev,
							data[maxIndices[n]], fitResult.getError(), noise);

					// Add to filtered output results
					if (filter != null && npeaks != 0)
					{
						// Check all result peaks for the distance to the filter positions
						Iterator<PeakResult> iter = sliceResults.descendingIterator();
						PeakResult[] results = new PeakResult[npeaks];
						while (iter.hasNext() && npeaks-- > 0)
						{
							results[npeaks] = iter.next();
						}
						filter.filter(fitResult, maxIndices[n], results);
					}

					// Q. Should this be a failure if the npeaks == 0 (i.e. outside the border; was a duplicate)?
					// At the moment just let the fitting continue until real failures occur.
					failures = 0;
					success++;
				}
				else
				{
					// Add to filtered output results
					if (filter != null)
					{
						// Check the start position for the distance to the filter positions
						filter.filter(fitResult, maxIndices[n], x + 0.5f, y + 0.5f);
					}

					// Add the failed jobs to a list.
					failed[failedCount++] = n;
				}

				if (config.getFailuresLimit() != 0 && failures >= config.getFailuresLimit())
					break;
			}

			if (logger != null)
				logger.info("Slice %d: %d / %d", slice, success, maxIndices.length);
		}

		float offsetx = bounds.x;
		float offsety = bounds.y;

		if (params != null && params.getOffset() != null)
		{
			offsetx += params.getOffset()[0];
			offsety += params.getOffset()[1];
		}

		// Add the ROI bounds to the fitted peaks
		for (PeakResult result : sliceResults)
		{
			result.params[Gaussian2DFunction.X_POSITION] += offsetx;
			result.params[Gaussian2DFunction.Y_POSITION] += offsety;
		}

		// Check if only selected peaks should be added to the results
		if (filter != null)
		{
			filter.finalise();
			job.setResults(filter.getResults());
			job.setIndices(filter.getMaxIndices());

			for (int i = 0; i < filter.getFilteredCount(); i++)
				job.setFitResult(i, filter.getFitResults()[i]);

			this.results.addAll(filter.getResults());
		}
		else
		{
			this.results.addAll(sliceResults);
		}

		finish(job, start);
	}

	/**
	 * Estimate the noise in the data
	 * 
	 * @param data
	 * @param width
	 * @param height
	 * @param method
	 * @return The noise
	 */
	public static float estimateNoise(float[] data, int width, int height, NoiseEstimator.Method method)
	{
		NoiseEstimator ne = new NoiseEstimator(data, width, height);
		return (float) ne.getNoise(method);
	}

	private void finish(FitJob job, long start)
	{
		time += System.nanoTime() - start;
		job.finished();
	}

	private float[] applySmoothing(float[] data, int width, int height, double smooth)
	{
		// Smoothing destructively modifies the data so create a copy
		float[] smoothData = Arrays.copyOf(data, width * height);

		if (smooth < 1)
		{
			// Do a 3x3 Gaussian approximation for less than 1 pixel smoothing
			//   1 2 1
			//   2 4 2
			//   1 2 1
			if (smooth <= border)
				filter.blockGaussian3x3Internal(smoothData, width, height);
			else
				filter.blockGaussian3x3(smoothData, width, height);
		}
		else
		{
			// Check upper limits are safe
			int tmpSmooth = Math.min((int) smooth, Math.min(width, height) / 2);

			if (tmpSmooth <= border)
			{
				filter.rollingBlockAverageInternal(smoothData, width, height, tmpSmooth);
			}
			else
			{
				filter.rollingBlockAverage(smoothData, width, height, tmpSmooth);
			}
		}
		return smoothData;
	}

	private void allocateArrays(int length)
	{
		xcoord = allocateArray(xcoord, length);
		ycoord = allocateArray(ycoord, length);
		neighbourIndices = allocateArray(neighbourIndices, length);
		// Allocate enough room for all fits to be doublets 
		fittedNeighbourIndices = allocateArray(fittedNeighbourIndices, length *
				((config.getResidualsThreshold() < 1) ? 2 : 1));
	}

	private int[] allocateArray(int[] array, int length)
	{
		if (array == null || array.length < length)
			array = new int[length];
		return array;
	}

	private int addResults(LinkedList<PeakResult> peakResults, int n, int x, int y, Rectangle bounds,
			Rectangle regionBounds, float[] peakParams, float[] peakParamsDev, float value, double error, float noise)
	{
		x += bounds.x;
		y += bounds.y;

		int npeaks = peakParams.length / 6;
		int count = 0;

		// Note that during processing the data is assumed to refer to the top-left
		// corner of the pixel. The coordinates should be represented in the middle of the pixel 
		// so add a 0.5 shift to the coordinates.

		if (npeaks == 1)
		{
			peakParams[Gaussian2DFunction.X_POSITION] += 0.5f;
			peakParams[Gaussian2DFunction.Y_POSITION] += 0.5f;
			if (addSingleResult(peakResults, x, y, regionBounds, peakParams, peakParamsDev, value, error, noise))
				count++;
		}
		else
		{
			final int currentResultsSize = peakResults.size();
			for (int i = 0; i < npeaks; i++)
			{
				float[] params = new float[] { peakParams[0], peakParams[i * 6 + 1], peakParams[i * 6 + 2],
						peakParams[i * 6 + 3] + 0.5f, peakParams[i * 6 + 4] + 0.5f, peakParams[i * 6 + 5],
						peakParams[i * 6 + 6] };
				float[] paramsStdDev = null;
				if (peakParamsDev != null)
				{
					paramsStdDev = new float[] { peakParamsDev[0], peakParamsDev[i * 6 + 1], peakParamsDev[i * 6 + 2],
							peakParamsDev[i * 6 + 3], peakParamsDev[i * 6 + 4], peakParamsDev[i * 6 + 5],
							peakParamsDev[i * 6 + 6] };
				}
				if (addSingleResult(peakResults, currentResultsSize, x, y, regionBounds, params, paramsStdDev, value,
						error, noise))
					count++;
			}
		}
		return count;
	}

	/**
	 * Add the result to the list
	 * 
	 * @param peakResults
	 * @param x
	 * @param y
	 * @param regionBounds
	 * @param peakParams
	 * @param peakParamsDev
	 * @param value
	 * @param error
	 * @param noise
	 * @return
	 */
	private boolean addSingleResult(LinkedList<PeakResult> peakResults, int x, int y, Rectangle regionBounds,
			float[] peakParams, float[] peakParamsDev, float value, double error, float noise)
	{
		peakParams[Gaussian2DFunction.X_POSITION] += regionBounds.x;
		peakParams[Gaussian2DFunction.Y_POSITION] += regionBounds.y;

		// Check if the position is inside the border
		if (insideBorder(peakParams[Gaussian2DFunction.X_POSITION], peakParams[Gaussian2DFunction.Y_POSITION]))
		{
			// Check for duplicates
			if (duplicateDistance2 > 0)
			{
				for (PeakResult result : peakResults)
				{
					if (distance2(result.params, peakParams) < duplicateDistance2)
					{
						//System.out.printf("Duplicate  [%d] %.2f,%.2f\n", slice,
						//		peakParams[Gaussian2DFunction.X_POSITION],
						//		peakParams[Gaussian2DFunction.Y_POSITION]);
						return false;
					}
				}
			}

			peakResults.add(createResult(x, y, value, error, noise, peakParams, peakParamsDev));
			fittedBackground += peakParams[Gaussian2DFunction.BACKGROUND];
			return true;
		}
		else if (logger != null)
		{
			logger.info("Ignoring peak within image border @ %.2f,%.2f", peakParams[Gaussian2DFunction.X_POSITION],
					peakParams[Gaussian2DFunction.Y_POSITION]);
		}
		return false;
	}

	/**
	 * Add the result to the list. Only check for duplicates up to the currentResultsSize
	 * 
	 * @param peakResults
	 * @param currentResultsSize
	 * @param x
	 * @param y
	 * @param regionBounds
	 * @param peakParams
	 * @param peakParamsDev
	 * @param value
	 * @param error
	 * @param noise
	 * @return
	 */
	private boolean addSingleResult(LinkedList<PeakResult> peakResults, final int currentResultsSize, int x, int y,
			Rectangle regionBounds, float[] peakParams, float[] peakParamsDev, float value, double error, float noise)
	{
		peakParams[Gaussian2DFunction.X_POSITION] += regionBounds.x;
		peakParams[Gaussian2DFunction.Y_POSITION] += regionBounds.y;

		// Check if the position is inside the border
		if (insideBorder(peakParams[Gaussian2DFunction.X_POSITION], peakParams[Gaussian2DFunction.Y_POSITION]))
		{
			// Check for duplicates
			if (duplicateDistance2 > 0)
			{
				int i = 0;
				for (PeakResult result : peakResults)
				{
					// Only check up to the current results size
					if (i++ == currentResultsSize)
						break;
					if (distance2(result.params, peakParams) < duplicateDistance2)
					{
						//System.out.printf("Duplicate  [%d] %.2f,%.2f\n", slice,
						//		peakParams[Gaussian2DFunction.X_POSITION],
						//		peakParams[Gaussian2DFunction.Y_POSITION]);
						return false;
					}
				}
			}

			peakResults.add(createResult(x, y, value, error, noise, peakParams, peakParamsDev));
			fittedBackground += peakParams[Gaussian2DFunction.BACKGROUND];
			return true;
		}
		else if (logger != null)
		{
			logger.info("Ignoring peak within image border @ %.2f,%.2f", peakParams[Gaussian2DFunction.X_POSITION],
					peakParams[Gaussian2DFunction.Y_POSITION]);
		}
		return false;
	}

	private PeakResult createResult(int origX, int origY, float origValue, double error, float noise, float[] params,
			float[] paramsStdDev)
	{
		if (endT >= 0 && slice != endT)
		{
			return new ExtendedPeakResult(slice, origX, origY, origValue, error, noise, params, paramsStdDev, endT, 0);
		}
		return new PeakResult(slice, origX, origY, origValue, error, noise, params, paramsStdDev);
	}

	/**
	 * @param x
	 * @param y
	 * @return True if the fitted position iis inside the border
	 */
	private boolean insideBorder(float x, float y)
	{
		return (x > border && x < borderLimitX && y > border && y < borderLimitY);
	}

	/**
	 * @param params1
	 * @param params2
	 * @return The squared distance between the two points
	 */
	private float distance2(float[] params1, float[] params2)
	{
		final float dx = params1[Gaussian2DFunction.X_POSITION] - params2[Gaussian2DFunction.X_POSITION];
		final float dy = params1[Gaussian2DFunction.Y_POSITION] - params2[Gaussian2DFunction.Y_POSITION];
		return dx * dx + dy * dy;
	}

	/**
	 * Find the indices of the maxima using the currently configured parameters
	 * <p>
	 * Data must be arranged in yx block order, i.e. height rows of width.
	 * 
	 * @param data
	 * @param width
	 * @param height
	 * @return Indices of the maxima
	 */
	private int[] getMaxima(float[] data, int width, int height)
	{
		// Find maxima
		NonMaximumSuppression nms = new NonMaximumSuppression();

		int[] maxIndices;

		// Do a neighbour check when using a low block size
		nms.setNeighbourCheck(search < 3);

		// Check upper limits are safe
		int n = Math.min(search, Math.min(width, height));
		int border = Math.min(this.border, Math.min(width, height) / 2);

		maxIndices = nms.blockFindInternal(data, width, height, n, border);

		return maxIndices;
	}

	private float[] convertParameters(float[] params)
	{
		// Convert radians to degrees (if elliptical fitting)
		if (fitConfig.isAngleFitting())
		{
			for (int i = 6; i < params.length; i += 6)
			{
				params[i - 4] *= 180.0 / Math.PI;
			}
		}
		return params;
	}

	/**
	 * Fits a 2D Gaussian to the given data. Fits all the specified peaks.
	 * <p>
	 * Data must be arranged in yx block order, i.e. height rows of width.
	 * <p>
	 * Can perform differential residue analysis on fit results to detect doublets. These are re-fitted and accepted if
	 * the residuals are significantly improved.
	 * 
	 * @param region
	 * @param w
	 * @param height
	 * @param maxIndices
	 *            All the indices of the maxima
	 * @param n
	 *            The maxima
	 * @param smoothData
	 *            The smoothed image data
	 * @return Array containing the fitted curve data: The first value is the Background. The remaining values are
	 *         Amplitude, Angle, PosX, PosY, StdDevX, StdDevY for each fitted peak.
	 *         <p>
	 *         Null if no fit is possible.
	 */
	private FitResult fit(Gaussian2DFitter gf, float[] region, Rectangle regionBounds, int[] maxIndices, int n,
			float[] smoothData)
	{
		int width = regionBounds.width;
		int height = regionBounds.height;
		int x = xcoord[n];
		int y = ycoord[n];

		// Analyse neighbours and include them in the fit if they are within a set height of this peak.
		int neighbours = (config.isIncludeNeighbours()) ? findNeighbours(regionBounds, n, x, y, smoothData, maxIndices)
				: 0;

		// Estimate background
		float background = 0;
		if (fittedNeighbourCount > 0)
		{
			// Use the average previously fitted background

			// Add the details of the already fitted peaks
			for (int i = 0; i < fittedNeighbourCount; i++)
			{
				PeakResult result = sliceResults.get(fittedNeighbourIndices[i]);
				background += result.params[Gaussian2DFunction.BACKGROUND];
			}

			background /= fittedNeighbourCount;
		}
		else
		{
			if (!sliceResults.isEmpty())
			{
				// Use the average background from all results
				background = (float) (fittedBackground / sliceResults.size());
			}
			else
			{
				// Initial guess using image mean
				background = this.background;
			}
		}

		FitResult fitResult = null;
		if (neighbours > 0)
		{
			boolean subtractFittedPeaks = false;

			// Multiple-fit ...
			if (logger != null)
				logger.info("Slice %d: Multiple-fit (%d peaks : neighbours [%d + %d])", slice, neighbours + 1,
						neighbourCount, fittedNeighbourCount);

			neighbours = (subtractFittedPeaks) ? neighbourCount : neighbourCount + fittedNeighbourCount;

			// Create the parameters for the fit
			final int parametersPerPeak = 6;
			int npeaks = 1 + neighbours;
			float[] params = new float[1 + npeaks * parametersPerPeak];

			// Note: If difference-of-smoothing is performed the heights have background subtracted so 
			// it must be added back 

			// The main peak
			params[Gaussian2DFunction.AMPLITUDE] = smoothData[maxIndices[n]] +
					((differenceOfSmoothing) ? background : 0);
			params[Gaussian2DFunction.X_POSITION] = x - regionBounds.x;
			params[Gaussian2DFunction.Y_POSITION] = y - regionBounds.y;

			// The neighbours
			for (int i = 0, j = parametersPerPeak; i < neighbourCount; i++, j += parametersPerPeak)
			{
				int n2 = neighbourIndices[i];
				params[j + Gaussian2DFunction.AMPLITUDE] = smoothData[maxIndices[n2]] +
						((differenceOfSmoothing) ? background : 0);
				params[j + Gaussian2DFunction.X_POSITION] = xcoord[n2] - regionBounds.x;
				params[j + Gaussian2DFunction.Y_POSITION] = ycoord[n2] - regionBounds.y;
			}

			if (!differenceOfSmoothing)
			{
				// In the case of uneven illumination the peaks may be below the background.
				// Check the heights are positive. Otherwise set it to zero and it will be re-estimated
				// in the fitter.
				for (int i = 0, j = 0; i <= neighbourCount; i++, j += parametersPerPeak)
				{
					if (params[j + Gaussian2DFunction.AMPLITUDE] < background)
					{
						background = 0;
						break;
					}
				}
			}

			// The fitted neighbours
			// TODO - Test which is better: (1) subtracting fitted peaks; or (2) including them
			// Initial tests show that Chi-squared is much lower when including them in the fit.

			if (fittedNeighbourCount > 0)
			{
				if (subtractFittedPeaks)
				{
					// Subtract the already fitted peaks from the data. 
					// This will speed up evaluation of the fitting function.
					region = Arrays.copyOf(region, width * height);

					float[] funcParams = new float[1 + parametersPerPeak * fittedNeighbourCount];
					for (int i = 0, j = 0; i < fittedNeighbourCount; i++, j += parametersPerPeak)
					{
						PeakResult result = sliceResults.get(fittedNeighbourIndices[i]);
						// Copy Amplitude,Angle,Xpos,Ypos,Xwidth,Ywidth
						for (int k = 1; k <= 6; k++)
							funcParams[j + k] = result.params[k];
						// Adjust position relative to extracted region
						funcParams[j + Gaussian2DFunction.X_POSITION] -= regionBounds.x;
						funcParams[j + Gaussian2DFunction.Y_POSITION] -= regionBounds.y;
					}

					GaussianFunction func = fitConfig.createGaussianFunction(fittedNeighbourCount, regionBounds.width,
							funcParams);
					func.initialise(funcParams);

					// Subtract fitted peaks
					float[] dyda = new float[funcParams.length];
					for (int i = 0; i < region.length; i++)
					{
						region[i] -= func.eval(i, dyda);
					}
				}
				else
				{
					// Add the details of the already fitted peaks
					for (int i = 0, j = (1 + neighbourCount) * parametersPerPeak; i < fittedNeighbourCount; i++, j += parametersPerPeak)
					{
						PeakResult result = sliceResults.get(fittedNeighbourIndices[i]);
						// Copy Amplitude,Angle,Xpos,Ypos,Xwidth,Ywidth
						for (int k = 1; k <= parametersPerPeak; k++)
							params[j + k] = result.params[k];

						// Add background to amplitude (required since the background will be re-estimated)
						params[j + Gaussian2DFunction.AMPLITUDE] += result.params[Gaussian2DFunction.BACKGROUND];
						// Adjust position relative to extracted region
						params[j + Gaussian2DFunction.X_POSITION] -= regionBounds.x;
						params[j + Gaussian2DFunction.Y_POSITION] -= regionBounds.y;
					}
				}
			}

			// Subtract the background from all estimated peak amplitudes
			if (background != 0)
			{
				params[0] = background;
				for (int j = Gaussian2DFunction.AMPLITUDE; j < params.length; j += parametersPerPeak)
					params[j] -= background;
			}

			// Increase the iterations for a multiple fit.

			int maxIterations = fitConfig.getMaxIterations();
			fitConfig.setMaxIterations(maxIterations + maxIterations * (npeaks - 1) *
					ITERATION_INCREASE_FOR_MULTIPLE_PEAKS);

			fitResult = gf.fit(region, width, height, npeaks, params);

			printFitResults(fitResult, region, width, height, npeaks, 0, gf.getIterations());

			fitConfig.setMaxIterations(maxIterations);

			if (fitResult.getStatus() == FitStatus.OK)
			{
				// Extract the first fitted peak
				final int degreesOfFreedom = fitResult.getDegreesOfFreedom();
				double error = fitResult.getError();
				final float[] initialParameters = truncate(fitResult.getInitialParameters());
				final float[] parameters = truncate(fitResult.getParameters());
				final float[] parametersDev = (fitResult.getParameterStdDev() != null) ? truncate(fitResult
						.getParameterStdDev()) : null;
				final int nPeaks = 1;
				final int nFittedParameters = fitResult.getNumberOfFittedParameters() % 6 +
						fitResult.getNumberOfFittedParameters() / npeaks;

				double r2 = 1 - (gf.getFinalResidualSumOfSquares() / gf.getTotalSumOfSquares());
				error = r2;

				fitResult = new FitResult(FitStatus.OK, degreesOfFreedom, error, initialParameters, parameters,
						parametersDev, nPeaks, nFittedParameters, null);
			}
			else
			{
				if (logger != null)
					logger.info("Multiple-fit failed, resorting to single-fit");
				fitResult = null;
			}
		}

		if (fitResult == null)
		{
			// Single fit
			//int newIndex = (y - regionBounds.y) * regionBounds.width + x - regionBounds.x;
			//int[] peaks = new int[] { newIndex };
			//float[] heights = new float[] { smoothData[maxIndices[n]] + ((differenceOfSmoothing) ? 0 : background) };
			//fitResult = gf.fit(region, width, height, peaks, heights);

			float amplitude = smoothData[maxIndices[n]] - ((differenceOfSmoothing) ? 0 : background);
			if (amplitude < 0)
			{
				amplitude += background;
				background = 0;
			}
			float[] params = new float[] { background, amplitude, 0, x - regionBounds.x, y - regionBounds.y, 0, 0 };
			fitResult = gf.fit(region, width, height, 1, params);
			printFitResults(fitResult, region, width, height, 1, -neighbours, gf.getIterations());

			final int degreesOfFreedom = fitResult.getDegreesOfFreedom();
			double error = fitResult.getError();
			final float[] initialParameters = fitResult.getInitialParameters();
			final float[] parameters = fitResult.getParameters();
			final float[] parametersDev = fitResult.getParameterStdDev();
			final int nPeaks = 1;
			final int nFittedParameters = fitResult.getNumberOfFittedParameters();

			double r2 = 1 - (gf.getFinalResidualSumOfSquares() / gf.getTotalSumOfSquares());
			error = r2;

			fitResult = new FitResult(fitResult.getStatus(), degreesOfFreedom, error, initialParameters, parameters,
					parametersDev, nPeaks, nFittedParameters, fitResult.getStatusData());

			// Only compute quadrant improvement if fitted as a single peak. If fitting multiple peaks
			// then we would expect the quadrant to be skewed by the neighbours.
			// TODO - Subtract the fitted peaks surrounding the central peak and then perform quadrant analysis
			// Irrespective of the neighbours.
			if (//neighbours == 0 && 
			canPerformQuadrantAnalysis(fitResult, width, height))
			{
				FitResult newFitResult = quadrantAnalysis(fitResult, region, regionBounds, smoothData);
				if (newFitResult != null)
					fitResult = newFitResult;
			}
		}

		if (logger != null)
		{
			switch (fitResult.getStatus())
			{
				case OK:
					// Show the shift, signal and width spread
					float[] initialParams = fitResult.getInitialParameters();
					float[] params = fitResult.getParameters();
					for (int i = 0, j = 0; i < fitResult.getNumberOfPeaks(); i++, j += 6)
					{
						final float factor = (float) (2 * Math.PI);
						float signal = factor * params[j + Gaussian2DFunction.AMPLITUDE] *
								params[j + Gaussian2DFunction.X_SD] * params[j + Gaussian2DFunction.Y_SD];

						logger.info(
								"Fit OK [%d]. Shift = %.3f,%.3f : SNR = %.2f : Width = %.2f,%.2f",
								i,
								params[j + Gaussian2DFunction.X_POSITION] -
										initialParams[j + Gaussian2DFunction.X_POSITION],
								params[j + Gaussian2DFunction.Y_POSITION] -
										initialParams[j + Gaussian2DFunction.Y_POSITION],
								signal / noise,
								getFactor(params[j + Gaussian2DFunction.X_SD], initialParams[j +
										Gaussian2DFunction.X_SD]),
								getFactor(params[j + Gaussian2DFunction.Y_SD], initialParams[j +
										Gaussian2DFunction.Y_SD]));
					}
					break;

				case BAD_PARAMETERS:
					float[] p = fitResult.getInitialParameters();
					logger.info("Bad parameters: %f,%f,%f,%f,%f,%f,%f", p[0], p[1], p[2], p[3], p[4], p[5], p[6]);
					break;

				default:
					logger.info(fitResult.getStatus().toString());
					break;
			}
		}

		return fitResult;
	}

	private void printFitResults(FitResult fitResult, float[] region, int width, int height, int npeaks, int doublet,
			int iterations)
	{
		// **********
		// Comment out for production code
		// **********

		//		if (fitResult.getStatus() == FitStatus.BAD_PARAMETERS)
		//			return;
		//
		//		int length = width * height;
		//		//double adjustedR2 = getAdjustedCoefficientOfDetermination(gf.getFinalResdiualSumOfSquares(),
		//		//		gf.getTotalSumOfSquares(), length, fitResult.getInitialParameters().length);
		//		double ic = getInformationCriterion(gf.getFinalResdiualSumOfSquares(), length,
		//				fitResult.getInitialParameters().length);
		//
		//		System.out.printf("[%dx%d] p %d i %d d %d : SS %f : %f : %f (%s). %f\n", width, height, npeaks, iterations,
		//				doublet, gf.getTotalSumOfSquares(), gf.getInitialResdiualSumOfSquares(),
		//				gf.getFinalResdiualSumOfSquares(), fitResult.getStatus().toString(), ic);
	}

	private void printFitResults(FitResult fitResult, float[] region, int width, int height, int npeaks, int doublet,
			int iterations, double AICc1, double AICc2)
	{
		// **********
		// Comment out for production code
		// **********

		//		if (fitResult.getStatus() == FitStatus.BAD_PARAMETERS)
		//			return;
		//
		//		System.out.printf("[%dx%d] p %d i %d d %d : SS %f : %f : %f (%s). %f -> %f\n", width, height, npeaks,
		//				iterations, doublet, gf.getTotalSumOfSquares(), gf.getInitialResdiualSumOfSquares(),
		//				gf.getFinalResdiualSumOfSquares(), fitResult.getStatus().toString(), AICc1, AICc2);
	}

	/**
	 * Get the relative change factor between f and g
	 * 
	 * @param f
	 * @param g
	 * @return The relative change factor (negative if g is bigger than f)
	 */
	private static float getFactor(float f, float g)
	{
		if (f > g)
			return f / g;
		return -g / f;
	}

	private float[] truncate(float[] array)
	{
		float[] newArray = new float[7];
		for (int i = 0; i < newArray.length; i++)
			newArray[i] = array[i];
		return newArray;
	}

	/**
	 * Search for any peak within a set height of the specified peak that is within the search region bounds
	 * 
	 * @param regionBounds
	 * @param n
	 * @param x
	 * @param y
	 * @param smoothData
	 * @param maxIndices
	 * @return The number of neighbours
	 */
	private int findNeighbours(Rectangle regionBounds, int n, int x, int y, float[] smoothData, int[] maxIndices)
	{
		int xmin = regionBounds.x;
		int xmax = xmin + regionBounds.width - 1;
		int ymin = regionBounds.y;
		int ymax = ymin + regionBounds.height - 1;

		float heightThreshold;
		if (differenceOfSmoothing)
		{
			// No background when performing difference-of-smoothing
			heightThreshold = (float) (smoothData[maxIndices[n]] * config.getNeighbourHeightThreshold());
		}
		else
		{
			if (smoothData[maxIndices[n]] < background)
				heightThreshold = smoothData[maxIndices[n]];
			else
				heightThreshold = (float) ((smoothData[maxIndices[n]] - background) *
						config.getNeighbourHeightThreshold() + background);
		}

		// TODO - Speed this up by storing the maxima in a grid and only compare the neighbouring grid cells

		// Check all maxima that are lower than this
		neighbourCount = 0;
		for (int i = n + 1; i < maxIndices.length; i++)
		{
			if (canIgnore(xcoord[i], ycoord[i], xmin, xmax, ymin, ymax, smoothData[maxIndices[i]], heightThreshold))
				continue;
			neighbourIndices[neighbourCount++] = i;
		}

		// Check all existing maxima. 
		// Since these will be higher than the current peak it is prudent to extend the range that should be considered.
		// Use the configured peak standard deviation.

		fittedNeighbourCount = 0;
		if (!sliceResults.isEmpty())
		{
			float range = Math.max(fitConfig.getInitialPeakStdDev0(), fitConfig.getInitialPeakStdDev1());
			float xmin2 = xmin - range;
			float xmax2 = xmax + range;
			float ymin2 = ymin - range;
			float ymax2 = ymax + range;
			for (int i = 0; i < sliceResults.size(); i++)
			{
				PeakResult result = sliceResults.get(i);
				// No height threshold check as this is a validated peak
				if (canIgnore(result.params[Gaussian2DFunction.X_POSITION],
						result.params[Gaussian2DFunction.Y_POSITION], xmin2, xmax2, ymin2, ymax2))
					continue;
				fittedNeighbourIndices[fittedNeighbourCount++] = i;
			}
		}

		return neighbourCount + fittedNeighbourCount;
	}

	private boolean canIgnore(int x, int y, int xmin, int xmax, int ymin, int ymax, float height, float heightThreshold)
	{
		return (x < xmin || x > xmax || y < ymin || y > ymax || height < heightThreshold);
	}

	private boolean canIgnore(float x, float y, float xmin, float xmax, float ymin, float ymax)
	{
		return (x < xmin || x > xmax || y < ymin || y > ymax);
	}

	private boolean canPerformQuadrantAnalysis(FitResult fitResult, final int width, final int height)
	{
		// Only do this for OK/BAD_FIT results 
		if (fitConfig.isComputeResiduals() && config.getResidualsThreshold() < 1)
		{
			if (fitResult.getStatus() == FitStatus.OK)
				return true;
			// Check why it was a bad fit. If it due to width divergence then 
			// check the width is reasonable given the size of the fitted region.
			if (fitResult.getStatus() == FitStatus.WIDTH_DIVERGED)
			{
				float[] params = fitResult.getParameters();
				float regionSize = Math.max(width, height) / 2;
				//int tmpSmooth = (int) Math.max(smooth, 1);
				//float regionSize = 2 * tmpSmooth + 1;
				if ((params[Gaussian2DFunction.X_SD] > 0 && params[Gaussian2DFunction.X_SD] < regionSize) ||
						(params[Gaussian2DFunction.Y_SD] > 0 && params[Gaussian2DFunction.Y_SD] < regionSize))
					return true;
			}
		}
		return false;
	}

	private FitResult quadrantAnalysis(FitResult fitResult, float[] region, Rectangle regionBounds, float[] smoothData)
	{
		// Perform quadrant analysis as per rapidSTORM:
		/*
		 * When two fluorophores emit close to each other, typically the nonlinear fit will result in a suspected
		 * fluorophore position midway between the two fluorophores and with a high amplitude. In this case, the fit
		 * results show a characteristic handle structure: The two true fluorophore emissions leave slightly positive
		 * resiudes, while there are negative residues on an axis perpendicular to the one connecting the fluorophores.
		 * 
		 * This condition is detected well by quadrant-differential residue analysis: The residue matrix is divided into
		 * quadrants, with the pixels above both diagonals forming the upper quadrant, the pixels above the main and
		 * below the off diagonal forming the right quadrants and so on. Pixels right on the diagonals are ignored.
		 * Then, the opposing quadrants are summed, and these sums substracted from another, resulting in two quadrant
		 * differences: upper and lower minus right and left and right and left minus upper and lower. This process is
		 * repeated for the quadrants defined by the central row and the central column.
		 * 
		 * The maximum sum obtained in this way divided by the sum of the absolute quadrant contributions is an
		 * observable correlating highly with the true presence of double emitters. Also, the quadrants containing the
		 * positive contribution in the highest sum indicate along which axis the double emission happened.
		 */

		final int width = regionBounds.width;
		final int height = regionBounds.height;
		float[] params = fitResult.getParameters();
		int cx = (int) Math.round(params[Gaussian2DFunction.X_POSITION]);
		int cy = (int) Math.round(params[Gaussian2DFunction.Y_POSITION]);
		if (cx < 0 || cx >= width || cy < 0 || cy >= height) // Bad fits may be out of bounds
			return null;

		float[] residuals = gf.getResiduals();

		// Compute quadrants

		// X quadrant:
		// .AAA.   
		// D.A.B   
		// DD.BB
		// D.C.B
		// .CCC.
		double ABCD = 0;
		double A = 0, B = 0, C = 0, D = 0;
		for (int y = cy, x1 = cx, x2 = cx; y < height; y++, x1--, x2++)
		{
			for (int x = 0, index = y * width; x < width; x++, index++)
			{
				ABCD += Math.abs(residuals[index]);
				if (x < x1)
					D += residuals[index];
				else if (x < x2 && x > x1)
					C += residuals[index];
				else if (x > x2)
					B += residuals[index];
				else
					ABCD -= Math.abs(residuals[index]);
			}
		}
		for (int y = cy - 1, x1 = cx - 1, x2 = cx + 1; y >= 0; y--, x1--, x2++)
		{
			for (int x = 0, index = y * width; x < width; x++, index++)
			{
				ABCD += Math.abs(residuals[index]);
				if (x < x1)
					D += residuals[index];
				else if (x < x2 && x > x1)
					A += residuals[index];
				else if (x > x2)
					B += residuals[index];
				else
					ABCD -= Math.abs(residuals[index]);
			}
		}

		// Similar for + quadrants:
		// AA.BB
		// AA.BB
		// .....
		// DD.CC
		// DD.CC
		double ABCD2 = 0;
		double A2 = 0, B2 = 0, C2 = 0, D2 = 0;
		for (int y = cy + 1; y < height; y++)
		{
			for (int x = 0, index = y * width; x < width; x++, index++)
			{
				ABCD2 += Math.abs(residuals[index]);
				if (x < cx)
					D2 += residuals[index];
				else if (x > cx)
					C2 += residuals[index];
			}
		}
		for (int y = cy - 1; y >= 0; y--)
		{
			for (int x = 0, index = y * width; x < width; x++, index++)
			{
				ABCD2 += Math.abs(residuals[index]);
				if (x < cx)
					A2 += residuals[index];
				else if (x > cx)
					B2 += residuals[index];
			}
		}

		double AC = A + C;
		double BD = B + D;
		double score1 = Math.abs(AC - BD) / ABCD;

		double AC2 = A2 + C2;
		double BD2 = B2 + D2;
		double score2 = Math.abs(AC2 - BD2) / ABCD2;

		int[] vector;
		if (score1 > score2)
		{
			vector = (AC > BD) ? new int[] { 0, 1 } : new int[] { 1, 0 };
		}
		else
		{
			vector = (AC2 > BD2) ? new int[] { 1, 1 } : new int[] { 1, -1 };
		}
		double score = Math.max(score1, score2);

		if (logger != null)
			logger.info("Residue analysis = %f (%d,%d)", score, vector[0], vector[1]);

		// If differential residue analysis indicates a doublet then re-fit as two spots.
		if (score < config.getResidualsThreshold())
			return null;

		if (logger != null)
			logger.info("Computing 2-kernel model");

		// TODO - Locate the 2 new centres by moving out into the quadrant defined by the vector
		// and finding the maxima on the original image data.

		// Guess maxima using the fitted width as a single peak		
		int x1 = Math.round(cx + vector[0] * params[Gaussian2DFunction.X_SD]);
		int y1 = Math.round(cy + vector[1] * params[Gaussian2DFunction.Y_SD]);
		int x2 = Math.round(cx - vector[0] * params[Gaussian2DFunction.X_SD]);
		int y2 = Math.round(cy - vector[1] * params[Gaussian2DFunction.Y_SD]);

		// Check bounds
		if (x1 < 0)
			x1 = 0;
		else if (x1 >= regionBounds.width)
			x1 = regionBounds.width - 1;
		if (y1 < 0)
			y1 = 0;
		else if (y1 >= regionBounds.height)
			y1 = regionBounds.height - 1;

		// Check regionBounds
		if (x2 < 0)
			x2 = 0;
		else if (x2 >= regionBounds.width)
			x2 = regionBounds.width - 1;
		if (y2 < 0)
			y2 = 0;
		else if (y2 >= regionBounds.height)
			y2 = regionBounds.height - 1;

		// Check the two points are not the same. 
		if (x1 == x2 && y1 == y2)
		{
			//System.out.println("matching points");
			// This can only happen when the fitted with is zero due to the round() function 
			// moving to the next integer. If they are the same then the value should be cx,cy
			// and we can move along the vector.
			x1 += vector[0];
			y1 += vector[1];
			x2 -= vector[0];
			y2 -= vector[1];

			// Check regionBounds
			if (x1 < 0)
				x1 = 0;
			else if (x1 >= regionBounds.width)
				x1 = regionBounds.width - 1;
			if (y1 < 0)
				y1 = 0;
			else if (y1 >= regionBounds.height)
				y1 = regionBounds.height - 1;

			// Check regionBounds
			if (x2 < 0)
				x2 = 0;
			else if (x2 >= regionBounds.width)
				x2 = regionBounds.width - 1;
			if (y2 < 0)
				y2 = 0;
			else if (y2 >= regionBounds.height)
				y2 = regionBounds.height - 1;
		}

		// -*-*-
		// Allow the routine to estimate the background
		// -*-*-

		//int index1 = y1 * width + x1;
		//int originalIndex1 = (y1 + regionBounds.y) * bounds.width + x1 + regionBounds.x;
		//int index2 = y2 * width + x2;
		//int originalIndex2 = (y2 + regionBounds.y) * bounds.width + x2 + regionBounds.x;

		//int[] peaks = new int[] { index1, index2 };
		//float[] heights = new float[] { smoothData[originalIndex1], smoothData[originalIndex2] };
		// -*-*-

		// -+-+-
		// Estimate params using the single fitted peak
		// -+-+-
		float[] doubletParams = new float[1 + 2 * 6];

		doubletParams[Gaussian2DFunction.BACKGROUND] = params[Gaussian2DFunction.BACKGROUND];
		doubletParams[Gaussian2DFunction.AMPLITUDE] = params[Gaussian2DFunction.AMPLITUDE] * 0.5f;
		doubletParams[Gaussian2DFunction.X_POSITION] = x1;
		doubletParams[Gaussian2DFunction.Y_POSITION] = y1;
		doubletParams[6 + Gaussian2DFunction.AMPLITUDE] = params[Gaussian2DFunction.AMPLITUDE] * 0.5f;
		doubletParams[6 + Gaussian2DFunction.X_POSITION] = x2;
		doubletParams[6 + Gaussian2DFunction.Y_POSITION] = y2;
		// -+-+-

		// Store the residual sum-of-squares from the previous fit
		double singleSumOfSquares = gf.getFinalResidualSumOfSquares();

		// Disable checking of position movement since we guessed the 2-peak location.
		// Increase the iterations level then reset afterwards.

		// TODO - Should width and signal validation be disabled too?
		float shift = fitConfig.getCoordinateShift();
		int maxIterations = fitConfig.getMaxIterations();

		fitConfig.setCoordinateShift(Math.min(width, height));
		fitConfig.setMaxIterations(maxIterations * ITERATION_INCREASE_FOR_DOUBLETS);

		//FitResult newFitResult = gf.fit(region, width, height, peaks, heights);
		FitResult newFitResult = gf.fit(region, width, height, 2, doubletParams);

		fitConfig.setCoordinateShift(shift);
		fitConfig.setMaxIterations(maxIterations);

		// Allow bad fits since these are due to peak width divergence or low signal for all peaks.
		if (newFitResult.getStatus() == FitStatus.OK) // || newFitResult.getResult() == Result.BAD_FIT)
		{
			// Check if the residuals are better using the adjusted coefficient of determination.
			// See http://www.mathworks.co.uk/help/matlab/data_analysis/linear-regression.html
			// Adjusted Coefficient of determination is not good for non-linear models. Use the 
			// bias corrected Akaike Information Criterion (AICc):
			// http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2892436/

			double doubleSumOfSquares = gf.getFinalResidualSumOfSquares();

			int length = width * height;
			//double SStotal = gf.getTotalSumOfSquares();
			//double adjustedR2 = getAdjustedCoefficientOfDetermination(singleSumOfSquares, SStotal, length,
			//		fitResult.getInitialParameters().length);
			//double newAdjustedR2 = getAdjustedCoefficientOfDetermination(doubleSumOfSquares, SStotal, length,
			//		newFitResult.getInitialParameters().length);

			double ic1 = getInformationCriterion(singleSumOfSquares, length, fitResult.getInitialParameters().length);
			double ic2 = getInformationCriterion(doubleSumOfSquares, length, newFitResult.getInitialParameters().length);

			if (logger != null)
				logger.info("Model improvement - Sum-of-squares (AIC) : %f (%f) => %f (%f) : %f", singleSumOfSquares,
						ic1, doubleSumOfSquares, ic2, ic1 - ic2);

			// Check if the predictive power of the model is better with two peaks:
			// AIC should be lower
			// Adjusted R2 should be higher
			if (ic2 > ic1)
			//if (newAdjustedR2 < adjustedR2)
			{
				printFitResults(newFitResult, region, width, height, 2, 1, gf.getIterations(), ic1, ic2);
				return null;
			}

			// Check if either coordinate is outside the fitted region.
			float[] newParams = newFitResult.getParameters();
			for (int n = 0; n < 2; n++)
			{
				// Note that during processing the data is assumed to refer to the top-left
				// corner of the pixel. The coordinates should be represented in the middle of the pixel 
				// so add a 0.5 shift to the coordinates.
				final float xpos = newParams[Gaussian2DFunction.X_POSITION + n * 6] + 0.5f;
				final float ypos = newParams[Gaussian2DFunction.Y_POSITION + n * 6] + 0.5f;
				if (xpos < 0 || xpos > regionBounds.width || ypos < 0 || ypos > regionBounds.height)
				{
					if (logger != null)
						logger.info("Fitted coordinates outside the fitted region (x %g <> %d || y %g <> %d)", xpos,
								regionBounds.width, ypos, regionBounds.height);
					printFitResults(newFitResult, region, width, height, 2, 1, gf.getIterations(), ic1, ic2);
					return null;
				}
			}

			// Check the distance of the peaks to the centre of the region. It is possible that a second peak
			// at the edge of the region has been fitted (note that no coordinate shift check was performed).
			if (shift != 0)
			{
				// Allow extra large shifts since this is a doublet and the centre will be in the middle of the pair
				if (fitConfig.isWidth0Fitting())
				{
					// Add the fitted standard deviation to the allowed shift
					shift += Math.max(params[Gaussian2DFunction.X_SD], params[Gaussian2DFunction.Y_SD]);
				}
				else
				{
					// Quadrant analysis only happens when there are no neighbours or the neighbour fit failed.
					if (config.isIncludeNeighbours())
						// Assume that neighbours are insignificant and allow the shift to span half of the 
						// fitted window.
						shift = 0.5f * Math.max(regionBounds.width, regionBounds.height);
					else
						// Expand the allowed shift by the configured SD fit to allow close by peaks to be
						// included. Duplicate filtering will eliminate fits onto close by neighbours.
						shift += 2 * Math.max(fitConfig.getInitialPeakStdDev0(), fitConfig.getInitialPeakStdDev1());
				}

				int[] position = new int[2];
				int nPeaks = 0;
				for (int n = 0; n < 2; n++)
				{
					float xShift = newParams[Gaussian2DFunction.X_POSITION + n * 6] -
							params[Gaussian2DFunction.X_POSITION];
					float yShift = newParams[Gaussian2DFunction.Y_POSITION + n * 6] -
							params[Gaussian2DFunction.Y_POSITION];
					if (Math.abs(xShift) > shift || Math.abs(yShift) > shift)
					{
						if (logger != null)
							logger.info("Bad peak %d: Fitted coordinates moved (x=%g,y=%g)", n, xShift, yShift);
					}
					else
					{
						position[nPeaks++] = n;
					}
				}

				printFitResults(newFitResult, region, width, height, 2, 3 + nPeaks, gf.getIterations(), ic1, ic2);
				if (nPeaks == 0)
				{
					return null;
				}

				// Copy the OK peaks into a new result
				float[] newInitialParams = newFitResult.getInitialParameters();
				float[] newParamStdDev = newFitResult.getParameterStdDev();

				float[] okInitialParams = new float[1 + nPeaks * 6];
				float[] okParams = new float[1 + nPeaks * 6];
				float[] okParamStdDev = new float[1 + nPeaks * 6];

				okInitialParams[0] = newInitialParams[0];
				okParams[0] = params[0];
				if (newParamStdDev != null)
					okParamStdDev[0] = newParamStdDev[0];

				int destPos = 1;
				for (int i = 0; i < nPeaks; i++)
				{
					int srcPos = position[i] * 6 + 1;
					System.arraycopy(newInitialParams, srcPos, okInitialParams, destPos, 6);
					System.arraycopy(newParams, srcPos, okParams, destPos, 6);
					if (newParamStdDev != null)
						System.arraycopy(newParamStdDev, srcPos, okParamStdDev, destPos, 6);
					destPos += 6;
				}

				final int nFittedParameters = newFitResult.getNumberOfFittedParameters() % 6 + nPeaks *
						newFitResult.getNumberOfFittedParameters() / 2;

				double error = newFitResult.getError();
				double r2 = 1 - (gf.getFinalResidualSumOfSquares() / gf.getTotalSumOfSquares());
				error = r2;

				return new FitResult(newFitResult.getStatus(), newFitResult.getDegreesOfFreedom(), error,
						okInitialParams, okParams, okParamStdDev, nPeaks, nFittedParameters,
						newFitResult.getStatusData());
			}
			else
			{
				printFitResults(newFitResult, region, width, height, 2, 2, gf.getIterations(), ic1, ic2);
			}

			// TODO - Should a check still be made for duplicate peaks? 
			// It could be that splitting the single into a double peak has resulted in a fit of a 
			// previously fitted peak.  

			return newFitResult;
		}

		return null;
	}

	/**
	 * @param sumOfSquaredResiduals
	 *            the sum of squared residuals from the nonlinear least-squares fit
	 * @param n
	 *            The number of data points
	 * @param p
	 *            The number of fitted parameters
	 * @return The Information Criterion
	 */
	private double getInformationCriterion(double sumOfSquaredResiduals, int n, int p)
	{
		final double logLikelihood = 0.5 * (-n * (Math.log(2 * Math.PI) + 1 - Math.log(n) + Math
				.log(sumOfSquaredResiduals)));

		// Note: The true bias corrected AIC is derived from the 2nd, 3rd and 4th derivatives of the 
		// negative log-likelihood function. This is complex and so is not implemented.
		// See: 
		// http://www.math.sci.hiroshima-u.ac.jp/stat/TR/TR11/TR11-06.pdf
		// http://www.sciencedirect.com/science/article/pii/S0024379512000821#

		// This paper explains that the AIC or BIC are much better than the Adjusted coefficient of determination
		// for model selection:
		// http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2892436/

		//double aic = 2.0 * p - 2.0 * logLikelihood;

		// The Bias Corrected Akaike Information Criterion (AICc)
		// http://en.wikipedia.org/wiki/Akaike_information_criterion#AICc
		// Assumes a univariate linear model. We have a multivariate model (x,y): does this matter?
		//aic = aic + (2.0 * p * (p + 1)) / (n - p - 1);

		// Optimised 
		final double aic = 2.0 * (p - logLikelihood) + (2.0 * p * (p + 1)) / (n - p - 1);

		// Bayesian Information Criterion (BIC), which gives a higher penalty on the number of parameters
		// http://en.wikipedia.org/wiki/Bayesian_information_criterion
		//final double bic = p * Math.log(n) - 2.0 * logLikelihood;

		return aic;
	}

	/**
	 * @param SSresid
	 *            Sum of squared residuals from the model
	 * @param SStotal
	 *            SStotal is the sum of the squared differences from the mean of the dependent variable (total sum of
	 *            squares)
	 * @param n
	 *            Number of observations
	 * @param d
	 *            Number of parameters in the model
	 * @return
	 */
	@SuppressWarnings("unused")
	private double getAdjustedCoefficientOfDetermination(double SSresid, double SStotal, int n, int d)
	{
		return 1 - (SSresid / SStotal) * ((n - 1) / (n - d - 1));
	}

	/**
	 * Get an estimate of the background level using the mean of image.
	 * 
	 * @param width
	 * @param height
	 * @return
	 */
	private float estimateBackground(int width, int height)
	{
		// Compute average of the entire image
		double sum = 0;
		for (int i = width * height; i-- > 0;)
			sum += data[i];
		return (float) sum / (width * height);
	}

	/**
	 * Get an estimate of the background level using the fitted peaks. If no fits available then estimate background
	 * using mean of image.
	 * 
	 * @param peakResults
	 * @param width
	 * @param height
	 * @return
	 */
	@SuppressWarnings("unused")
	private float estimateBackground(LinkedList<PeakResult> peakResults, int width, int height)
	{
		if (peakResults.size() > 0)
		{
			// Compute average background of the fitted peaks
			double sum = 0;
			for (PeakResult result : peakResults)
				sum += result.params[Gaussian2DFunction.BACKGROUND];
			float av = (float) sum / peakResults.size();
			if (logger != null)
				logger.info("Average background %f", av);
			return av;
		}
		else
		{
			// Compute average of the entire image
			double sum = 0;
			for (int i = width * height; i-- > 0;)
				sum += data[i];
			float av = (float) sum / (width * height);
			if (logger != null)
				logger.info("Image background %f", av);
			return av;
		}
	}

	/**
	 * Identify failed peaks that seem quite high, e.g. if above the background + 3X the noise.
	 * <p>
	 * Updates the input failed array to contain the candidates.
	 * 
	 * @param failed
	 * @param failedCount
	 * @param background
	 * @param noise
	 * @param maxIndices
	 * @param smoothData
	 * @return The number of re-fit candidates
	 */
	@SuppressWarnings("unused")
	private int identifyRefitCandidates(int[] failed, int failedCount, float background, float noise, int[] maxIndices,
			float[] smoothData)
	{
		int candidates = 0;
		float threshold = background + 3 * noise;
		for (int i = 0; i < failedCount; i++)
		{
			if (smoothData[maxIndices[failed[i]]] > threshold)
				failed[candidates++] = i;
		}
		return candidates;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Runnable#run()
	 */
	public void run()
	{
		try
		{
			while (!finished)
			{
				FitJob job = jobs.take();
				if (job == null || job.data == null || finished)
					break;
				run(job);
			}
		}
		catch (InterruptedException e)
		{
			System.out.println(e.toString());
			throw new RuntimeException(e);
		}
		finally
		{
			finished = true;
			//notifyAll();
		}
	}

	/**
	 * @return the total time used for fitting
	 */
	public long getTime()
	{
		return time;
	}

	/**
	 * Signal that the worker should end
	 */
	public void finish()
	{
		finished = true;
	}

	/**
	 * @return True if the worker has finished
	 */
	public boolean isFinished()
	{
		return finished;
	}
}
