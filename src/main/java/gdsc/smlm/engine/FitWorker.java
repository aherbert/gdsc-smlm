package gdsc.smlm.engine;

import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.concurrent.BlockingQueue;

import org.apache.commons.math3.util.FastMath;

import gdsc.core.ij.Utils;
import gdsc.core.logging.Logger;
import gdsc.core.utils.ImageExtractor;
import gdsc.core.utils.Maths;
import gdsc.core.utils.NoiseEstimator;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
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
import gdsc.smlm.filters.MaximaSpotFilter;
import gdsc.smlm.filters.Spot;
import gdsc.smlm.fitting.FitConfiguration;
import gdsc.smlm.fitting.FitResult;
import gdsc.smlm.fitting.FitSolver;
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.Gaussian2DFitter;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunction;
import gdsc.smlm.function.gaussian.GaussianOverlapAnalysis;
import gdsc.smlm.results.ExtendedPeakResult;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResults;
import gdsc.smlm.results.filter.BasePreprocessedPeakResult;
import gdsc.smlm.results.filter.MultiFilter;
import gdsc.smlm.results.filter.MultiPathFilter;
import gdsc.smlm.results.filter.MultiPathFitResult;
import gdsc.smlm.results.filter.PreprocessedPeakResult;

/**
 * Fits local maxima using a 2D Gaussian.
 */
public class FitWorker implements Runnable
{
	/**
	 * The number of additional iterations to use for multiple peaks
	 * <p>
	 * Testings on a idealised dataset of simulated data show that multiple peaks increases the iterations but the
	 * increase asymptotes. Initial rate is 2-fold for 7x7 region, decreasing to 1.5-fold for 17x17 region. Best
	 * solution is to add a set of iterations for each additional peak.
	 */
	public static final int ITERATION_INCREASE_FOR_MULTIPLE_PEAKS = 1; // 0 for no effect
	/**
	 * The number of additional iterations to use for doublets
	 * <p>
	 * Testings on a idealised dataset of simulated data show that fitting the doublets increases the iterations by
	 * approx 3.5-fold for 7x7 region, decreasing to 2.5-fold for 17x17 region. An additional doubling of iterations
	 * were used when the fit of the doublet resulted in one peak being eliminated for moving.
	 */
	public static final int ITERATION_INCREASE_FOR_DOUBLETS = 4; // 1 for no effect

	/** The number of additional evaluations to use for multiple peaks */
	public static final int EVALUATION_INCREASE_FOR_MULTIPLE_PEAKS = 1; // 0 for no effect

	/** The number of additional evaluations to use for doublets */
	public static final int EVALUATION_INCREASE_FOR_DOUBLETS = 4; // 1 for no effect

	private Logger logger = null, logger2 = null;
	private FitTypeCounter counter = null;
	private long time = 0;

	private MaximaSpotFilter spotFilter = null;
	private int fitting = 1;

	// Used for fitting
	private FitEngineConfiguration config;
	private FitConfiguration fitConfig;

	private PeakResults results;
	private BlockingQueue<FitJob> jobs;
	private Gaussian2DFitter gf;

	// Used for fitting methods
	private ArrayList<PeakResult> sliceResults;
	private boolean useFittedBackground = false;
	private double fittedBackground;
	private int slice;
	private int endT;
	private Rectangle bounds;
	private int border, borderLimitX, borderLimitY;
	private float[] data;
	private boolean relativeIntensity;
	private float noise, background;
	private final boolean calculateNoise;
	private boolean estimateSignal;
	private ResultGridManager gridManager = null;
	private Spot[] spotNeighbours = null;
	private PeakResult[] peakNeighbours = null;
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

	/**
	 * Store an estimate for a spot candidate. This may be aquired during multi fitting of neighbours.
	 */
	private class Estimate
	{
		final double d2;
		final double[] params;

		public Estimate(double d2, double[] params)
		{
			this.d2 = d2;
			this.params = params;
		}
	}

	private Estimate[] estimates = null;

	public FitWorker(FitEngineConfiguration config, PeakResults results, BlockingQueue<FitJob> jobs)
	{
		this.config = config;
		this.fitConfig = config.getFitConfiguration();
		this.results = results;
		this.jobs = jobs;
		this.logger = fitConfig.getLog();
		gf = new Gaussian2DFitter(fitConfig);
		duplicateDistance2 = (float) (fitConfig.getDuplicateDistance() * fitConfig.getDuplicateDistance());
		calculateNoise = config.getFitConfiguration().getNoise() <= 0;
		if (!calculateNoise)
		{
			noise = (float) config.getFitConfiguration().getNoise();
		}

		id = ID;
		ID += 1;
	}

	/**
	 * Set the parameters for smoothing the image, searching for maxima and fitting maxima. This should be called before
	 * the {@link #run()} method to configure the fitting.
	 * 
	 * @param spotFilter
	 *            The spot filter for identifying fitting candidates
	 * @param fitting
	 *            The block size to be used for fitting
	 */
	public void setSearchParameters(MaximaSpotFilter spotFilter, int fitting)
	{
		this.spotFilter = spotFilter;
		this.border = spotFilter.getBorder();
		this.fitting = fitting;
		this.relativeIntensity = !spotFilter.isAbsoluteIntensity();
		// We can estimate the signal for a single peak when the fitting window covers enough of the Gaussian
		estimateSignal = 2.5 * config.getHWHMMax() / Gaussian2DFunction.SD_TO_HWHM_FACTOR < fitting;
	}

	/**
	 * Locate all the peaks in the image specified by the fit job
	 * <p>
	 * WARNING: The FitWorker fits a sub-region of the data for each maxima. It then updates the FitResult parameters
	 * with an offset reflecting the position. The initialParameters are not updated with this offset unless configured.
	 * 
	 * @param job
	 *            The fit job
	 */
	public void run(FitJob job)
	{
		final long start = System.nanoTime();
		job.start();
		this.slice = job.slice;

		// Used for debugging
		//if (logger == null) logger = new gdsc.fitting.logging.ConsoleLogger();

		// Crop to the ROI
		bounds = job.bounds;
		final int width = bounds.width;
		final int height = bounds.height;
		borderLimitX = width - border;
		borderLimitY = height - border;
		data = job.data;

		FitParameters params = job.getFitParameters();
		this.endT = (params != null) ? params.endT : -1;

		Spot[] spots = indentifySpots(job, width, height, params);

		if (spots.length == 0)
		{
			finish(job, start);
			return;
		}

		fittedBackground = 0;

		// Always get the noise and store it with the results.
		if (params != null && !Float.isNaN(params.noise))
		{
			noise = params.noise;
			fitConfig.setNoise(noise);
		}
		else if (calculateNoise)
		{
			noise = estimateNoise(data, width, height, config.getNoiseMethod());
			fitConfig.setNoise(noise);
		}

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

		ResultFilter resultFilter = null;
		final ImageExtractor ie = new ImageExtractor(data, width, height);
		double[] region = null;

		if (params != null && params.fitTask == FitTask.MAXIMA_IDENITIFICATION)
		{
			final float sd0 = (float) fitConfig.getInitialPeakStdDev0();
			final float sd1 = (float) fitConfig.getInitialPeakStdDev1();
			for (int n = 0; n < spots.length; n++)
			{
				// Find the background using the perimeter of the data.
				// TODO - Perhaps the Gaussian Fitter should be used to produce the initial estimates but no actual fit done.
				// This would produce coords using the centre-of-mass.
				final int x = spots[n].x;
				final int y = spots[n].y;
				Rectangle regionBounds = ie.getBoxRegionBounds(x, y, fitting);
				region = ie.crop(regionBounds, region);
				final float b = (float) Gaussian2DFitter.getBackground(region, regionBounds.width, regionBounds.height,
						1);

				// Offset the coords to the centre of the pixel. Note the bounds will be added later.
				// Subtract the background to get the amplitude estimate then convert to signal.
				final float amplitude = spots[n].intensity - ((relativeIntensity) ? 0 : b);
				final float signal = (float) (amplitude * 2.0 * Math.PI * sd0 * sd1);
				final float[] peakParams = new float[] { b, signal, 0, x + 0.5f, y + 0.5f, sd0, sd1 };
				final int index = y * width + x;
				sliceResults.add(createResult(bounds.x + x, bounds.y + y, data[index], 0, noise, peakParams, null));
			}
		}
		else if (params != null && params.fitTask == FitTask.BENCHMARKING)
		{
			initialiseFitting(spots);

			// Fit all spot candidates using all fitting pathways

			// Since later fitting results depend on earlier fitting results we must choose which of
			// the pathways produces the best result. This is usually done by fit validation. However
			// fit validation may be severely crippled in order to generate fit results for all pathways.
			// In this case we can provide a MultiPathFilter to decide what results we should keep as we go...

			MultiPathFilter filter = params.benchmarkFilter;
			duplicateDistance2 = (float) (params.duplicateDistance * params.duplicateDistance);

			if (filter == null)
			{
				// Create a default filter using the standard FitConfiguration to ensure sensible fits
				// are stored as the current slice results.
				final FitConfiguration tmp = new FitConfiguration();
				final FitEngineConfiguration fec = new FitEngineConfiguration(tmp);
				filter = new MultiPathFilter(
						new MultiFilter(tmp.getMinPhotons(), (float) tmp.getSignalStrength(), tmp.getMinWidthFactor(),
								tmp.getWidthFactor(), tmp.getCoordinateShiftFactor(), 0, tmp.getPrecisionThreshold()),
						(tmp.isComputeResiduals()) ? fec.getResidualsThreshold() : 1);
			}

			filter.setup();

			// Extract each peak and fit individually
			for (int n = 0; n < spots.length; n++)
			{
				final int x = spots[n].x;
				final int y = spots[n].y;

				final Rectangle regionBounds = ie.getBoxRegionBounds(x, y, fitting);
				region = ie.crop(regionBounds, region);

				// Fit all paths. The method uses the filter to store the best results for the slice.
				MultiPathFitResult fitResult = benchmarkFit(gf, region, regionBounds, spots, n, filter);
				job.setMultiPathFitResult(n, fitResult);
			}
		}
		else
		{
			initialiseFitting(spots);

			// Perform the Gaussian fit
			int success = 0;

			// Allow the results to be filtered for certain peaks
			if (params != null && params.filter != null)
			{
				// TODO - Check if this works.
				resultFilter = new DistanceResultFilter(params.filter, params.distanceThreshold, spots.length);
				//filter = new OptimumDistanceResultFilter(params.filter, params.distanceThreshold, maxIndices.length);
			}

			// Extract each peak and fit individually until N consecutive failures
			for (int n = 0, failures = 0; n < spots.length && !finished; n++)
			{
				failures++;

				final int x = spots[n].x;
				final int y = spots[n].y;
				final int index = y * width + x;

				final Rectangle regionBounds = ie.getBoxRegionBounds(x, y, fitting);
				region = ie.crop(regionBounds, region);

				if (logger != null)
					logger.info("Fitting %d:%d [%d,%d]", slice, n + 1, x + bounds.x, y + bounds.y);

				FitResult fitResult = fit(gf, region, regionBounds, spots, n);
				job.setFitResult(n, fitResult);

				// Debugging
				if (logger2 != null)
				{
					double[] peakParams = fitResult.getParameters();
					if (peakParams != null)
					{
						peakParams = Arrays.copyOf(peakParams, peakParams.length);
						int npeaks = peakParams.length / 6;
						for (int i = 0; i < npeaks; i++)
						{
							peakParams[i * 6 + Gaussian2DFunction.X_POSITION] += 0.5 + regionBounds.x;
							peakParams[i * 6 + Gaussian2DFunction.Y_POSITION] += 0.5 + regionBounds.y;
						}
					}
					String msg = String.format("%d:%d [%d,%d] %s (%s) = %s\n", slice, n + 1, x + bounds.x, y + bounds.y,
							fitResult.getStatus(), fitResult.getStatusData(), Arrays.toString(peakParams));
					logger2.debug(msg);
				}

				// Q. Should we add the regionBounds offset here to the params and initialParams in the FitResult?

				// Check fit result
				if (fitResult.getStatus() == FitStatus.OK)
				{
					final double[] peakParams = fitResult.getParameters();
					convertParameters(peakParams);

					double[] peakParamsDev = null;
					if (fitConfig.isComputeDeviations())
					{
						peakParamsDev = fitResult.getParameterStdDev();
						if (peakParamsDev != null)
							convertParameters(peakParamsDev);
					}

					int npeaks = addResults(x, y, bounds, regionBounds, peakParams, peakParamsDev, data[index],
							fitResult.getError(), noise);

					// Add to filtered output results
					if (resultFilter != null && npeaks != 0)
					{
						// Check all result peaks for the distance to the filter positions
						PeakResult[] results = new PeakResult[npeaks];
						for (int i = sliceResults.size(); npeaks-- > 0;)
						{
							results[npeaks] = sliceResults.get(--i);
						}
						resultFilter.filter(fitResult, index, results);
					}

					// Q. Should this be a failure if the npeaks == 0 (i.e. outside the border; was a duplicate)
					// or if the filter removed a peak?
					// At the moment just let the fitting continue until real failures occur.
					failures = 0;
					success++;
				}
				else
				{
					// Add to filtered output results
					if (resultFilter != null)
					{
						// Check the start position for the distance to the filter positions
						resultFilter.filter(fitResult, index, x + 0.5f, y + 0.5f);
					}
				}

				if (config.getFailuresLimit() >= 0 && failures > config.getFailuresLimit())
					break;
			}

			if (logger != null)
				logger.info("Slice %d: %d / %d", slice, success, spots.length);
		}

		// Add the ROI bounds to the fitted peaks
		float offsetx = bounds.x;
		float offsety = bounds.y;

		if (params != null && params.getOffset() != null)
		{
			offsetx += params.getOffset()[0];
			offsety += params.getOffset()[1];
		}

		for (PeakResult result : sliceResults)
		{
			result.params[Gaussian2DFunction.X_POSITION] += offsetx;
			result.params[Gaussian2DFunction.Y_POSITION] += offsety;
		}

		// Check if only selected peaks should be added to the results
		if (resultFilter != null)
		{
			resultFilter.finalise();
			job.setResults(resultFilter.getResults());
			job.setIndices(resultFilter.getMaxIndices());

			for (int i = 0; i < resultFilter.getFilteredCount(); i++)
			{
				job.setFitResult(i, resultFilter.getFitResults()[i]);
			}

			this.results.addAll(resultFilter.getResults());
		}
		else
		{
			this.results.addAll(sliceResults);
		}

		finish(job, start);
	}

	private Spot[] indentifySpots(FitJob job, int width, int height, FitParameters params)
	{
		Spot[] spots = null;
		int[] maxIndices = null;

		// Only sub-classes may require the indices
		boolean requireIndices = (job.getClass() != FitJob.class);

		if (params != null)
		{
			if (params.spots != null)
			{
				spots = params.spots;
				// Get the indices
				maxIndices = new int[spots.length];
				for (int n = 0; n < maxIndices.length; n++)
				{
					maxIndices[n] = spots[n].y * width + spots[n].x;
				}
			}
			else if (params.maxIndices != null)
			{
				// Extract the desired spots
				maxIndices = params.maxIndices;
				float[] data2 = spotFilter.preprocessData(data, width, height);
				spots = new Spot[maxIndices.length];
				for (int n = 0; n < maxIndices.length; n++)
				{
					final int y = maxIndices[n] / width;
					final int x = maxIndices[n] % width;
					final float intensity = data2[maxIndices[n]];
					spots[n] = new Spot(x, y, intensity);
				}
				// Sort the maxima
				Arrays.sort(spots);
			}
		}

		if (spots == null)
		{
			// Run the filter to get the spot
			spots = spotFilter.rank(data, width, height);
			// Extract the indices
			if (requireIndices)
			{
				maxIndices = new int[spots.length];
				for (int n = 0; n < maxIndices.length; n++)
				{
					maxIndices[n] = spots[n].y * width + spots[n].x;
				}
			}
		}

		if (logger != null)
			logger.info("%d: Slice %d: %d candidates", id, slice, spots.length);

		sliceResults = new ArrayList<PeakResult>(spots.length);
		if (requireIndices)
		{
			job.setResults(sliceResults);
			job.setIndices(maxIndices);
		}
		return spots;
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

	private void finish(FitJob job, final long start)
	{
		time += System.nanoTime() - start;
		job.finished();
	}

	private void initialiseFitting(Spot[] spots)
	{
		neighbourIndices = allocateArray(neighbourIndices, spots.length);
		// Allocate enough room for all fits to be doublets 
		fittedNeighbourIndices = allocateArray(fittedNeighbourIndices,
				spots.length * ((config.getResidualsThreshold() < 1) ? 2 : 1));

		clearEstimates(spots.length);

		gridManager = new ResultGridManager(bounds.width, bounds.height, 2 * fitting + 1);
		for (int i = 0; i < spots.length; i++)
		{
			// Store ID in the score field
			spots[i].setScore(i);
			gridManager.addToGrid(spots[i]);
		}
	}

	private int[] allocateArray(int[] array, int length)
	{
		if (array == null || array.length < length)
			array = new int[length];
		return array;
	}

	private void clearEstimates(int length)
	{
		if (estimates == null || estimates.length < length)
		{
			estimates = new Estimate[length];
		}
		else
		{
			for (int i = 0; i < length; i++)
				estimates[i] = null;
		}
	}

	private int addResults(int x, int y, Rectangle bounds, Rectangle regionBounds, double[] peakParams,
			double[] peakParamsDev, float value, double error, float noise)
	{
		// Neighbours to check for duplicates
		final PeakResult[] neighbours = (duplicateDistance2 > 0) ? gridManager.getPeakResultNeighbours(x, y) : null;

		x += bounds.x;
		y += bounds.y;

		// Note the bounds will be added to the params at the end of processing

		final int npeaks = peakParams.length / 6;
		int count = 0;

		// Note that during processing the data is assumed to refer to the top-left
		// corner of the pixel. The coordinates should be represented in the middle of the pixel 
		// so add a 0.5 shift to the coordinates.

		// WARNING: We do not update the initialParameters in the FitResult with the same offset.

		for (int i = 0; i < npeaks; i++)
		{
			peakParams[i * 6 + Gaussian2DFunction.X_POSITION] += 0.5 + regionBounds.x;
			peakParams[i * 6 + Gaussian2DFunction.Y_POSITION] += 0.5 + regionBounds.y;
		}

		if (npeaks == 1)
		{
			if (addSingleResult(neighbours, x, y, Utils.toFloat(peakParams), Utils.toFloat(peakParamsDev), value, error,
					noise))
				count++;
		}
		else
		{
			for (int i = 0; i < npeaks; i++)
			{
				final float[] params = new float[] { (float) peakParams[0], (float) peakParams[i * 6 + 1],
						(float) peakParams[i * 6 + 2], (float) peakParams[i * 6 + 3], (float) peakParams[i * 6 + 4],
						(float) peakParams[i * 6 + 5], (float) peakParams[i * 6 + 6] };
				float[] paramsStdDev = null;
				if (peakParamsDev != null)
				{
					paramsStdDev = new float[] { (float) peakParamsDev[0], (float) peakParamsDev[i * 6 + 1],
							(float) peakParamsDev[i * 6 + 2], (float) peakParamsDev[i * 6 + 3],
							(float) peakParamsDev[i * 6 + 4], (float) peakParamsDev[i * 6 + 5],
							(float) peakParamsDev[i * 6 + 6] };
				}
				if (addSingleResult(neighbours, x, y, params, paramsStdDev, value, error, noise))
					count++;
			}
		}
		return count;
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
	private boolean addSingleResult(final PeakResult[] neighbours, int x, int y, float[] peakParams,
			float[] peakParamsDev, float value, double error, float noise)
	{
		// Check if the position is inside the border
		if (insideBorder(peakParams[Gaussian2DFunction.X_POSITION], peakParams[Gaussian2DFunction.Y_POSITION]))
		{
			// Check for duplicates
			if (neighbours != null)
			{
				for (int i = 0; i < neighbours.length; i++)
				{
					if (distance2(neighbours[i].params, peakParams) < duplicateDistance2)
					{
						//System.out.printf("Duplicate  [%d] %.2f,%.2f\n", slice,
						//		peakParams[Gaussian2DFunction.X_POSITION],
						//		peakParams[Gaussian2DFunction.Y_POSITION]);
						return false;
					}
				}
			}

			final PeakResult peakResult = createResult(x, y, value, error, noise, peakParams, peakParamsDev);
			gridManager.addToGrid(peakResult);
			sliceResults.add(peakResult);
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
	 * @return True if the fitted position is inside the border
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

	private double[] convertParameters(double[] params)
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
	 * Provide functionality to fit spots in a region using different methods
	 */
	private class SpotFitter
	{
		final Gaussian2DFitter gf;
		final double[] region;
		final Rectangle regionBounds;
		final Spot[] spots;
		final int n;
		final int width;
		final int height;
		final int x;
		final int y;
		final int parametersPerPeak = 6;

		int neighbours;
		float singleBackground = Float.NaN, multiBackground = Float.NaN;

		FitResult resultMulti = null;
		FitResult resultSingle = null;
		double[] singleRegion = null;
		FitResult resultDoublet = null;
		boolean computedDoublet = false;
		boolean computedDoubletQA = false;
		boolean storeEstimates = true;

		public SpotFitter(Gaussian2DFitter gf, double[] region, Rectangle regionBounds, Spot[] spots, int n)
		{
			this.gf = gf;
			this.region = region;
			this.regionBounds = regionBounds;
			this.spots = spots;
			this.n = n;

			// Initialise
			width = regionBounds.width;
			height = regionBounds.height;
			x = spots[n].x;
			y = spots[n].y;

			// Analyse neighbours and include them in the fit if they are within a set height of this peak.
			resetNeighbours();
			neighbours = (config.isIncludeNeighbours()) ? findNeighbours(regionBounds, n, x, y, spots) : 0;
		}

		private float getMultiFittingBackground()
		{
			if (Float.isNaN(multiBackground))
			{
				multiBackground = 0;
				if (fittedNeighbourCount > 0)
				{
					// Use the average previously fitted background

					// Add the details of the already fitted peaks
					for (int i = 0; i < fittedNeighbourCount; i++)
					{
						multiBackground += sliceResults.get(fittedNeighbourIndices[i]).getBackground();
					}

					multiBackground /= fittedNeighbourCount;
				}
				else
				{
					multiBackground = this.getSingleFittingBackground();
				}
			}
			return multiBackground;
		}

		private float getSingleFittingBackground()
		{
			if (Float.isNaN(singleBackground))
			{
				singleBackground = FitWorker.this.getSingleFittingBackground();
			}
			return singleBackground;
		}

		public FitResult getResultMulti()
		{
			if (neighbours == 0)
				return null;
			if (resultMulti != null)
				return resultMulti;

			// Estimate background.
			// Note that using the background from previous fit results leads to an 
			// inconsistency in the results when the fitting parameters are changed which may be unexpected, 
			// e.g. altering the max iterations.
			float background = getMultiFittingBackground();

			// TODO

			// If we have fitted neighbours:
			// subtract them from the region and then try a single fit. It should work if
			// something is there. This can be used as an initial estimate for one of the
			// peaks in the multiple fit (i.e. the closest one)

			// If the fits fails then we can guess that the region has no good peaks and 
			// it is not worth doing a multiple fit.			

			int subtractFittedPeaks = 0;
			boolean[] subtract = null;

			// Determine if any fitted neighbours are outside the region
			// Examples show that fitting spots with a centre 
			// outside the region can result in large drift.
			if (fittedNeighbourCount > 0)
			{
				subtract = new boolean[fittedNeighbourCount];

				// TODO - optionally subtract all fitted peaks including those inside the fit region.

				// The fitted result will be relative to (0,0) in the fit data and already 
				// have an offset applied so that 0.5 is the centre of a pixel. We can test 
				// the coordinates exactly against the fit frame.
				final float xmin = regionBounds.x;
				final float xmax = xmin + regionBounds.width;
				final float ymin = regionBounds.y;
				final float ymax = ymin + regionBounds.height;

				for (int i = 0; i < fittedNeighbourCount; i++)
				{
					final PeakResult result = sliceResults.get(fittedNeighbourIndices[i]);

					// Subtract peaks from the data if they are outside the fit region
					if (result.getXPosition() < xmin || result.getXPosition() > xmax || result.getYPosition() < ymin ||
							result.getYPosition() > ymax)
					{
						subtract[i] = true;
						subtractFittedPeaks++;
					}
				}
			}

			// Multiple-fit ...
			if (logger != null)
				logger.info("Slice %d: Multiple-fit (%d peaks : neighbours [%d + %d - %d])", slice, neighbours + 1,
						neighbourCount, fittedNeighbourCount, subtractFittedPeaks);

			neighbours = neighbourCount + fittedNeighbourCount - subtractFittedPeaks;

			// Create the parameters for the fit
			int npeaks = 1 + neighbours;
			double[] params = new double[1 + npeaks * parametersPerPeak];

			// Note: If difference-of-smoothing is performed the heights have background subtracted so 
			// it must be added back 

			// The main peak
			params[Gaussian2DFunction.SIGNAL] = spots[n].intensity + ((relativeIntensity) ? background : 0);
			params[Gaussian2DFunction.X_POSITION] = x - regionBounds.x;
			params[Gaussian2DFunction.Y_POSITION] = y - regionBounds.y;

			// The neighbours
			for (int i = 0, j = parametersPerPeak; i < neighbourCount; i++, j += parametersPerPeak)
			{
				final int n2 = neighbourIndices[i];
				final double[] estimatedParams = getEstimate(n2);
				if (estimatedParams != null)
				{
					// Re-use previous good multi-fit results to estimate the peak params...
					// Convert signal into amplitude
					params[j + Gaussian2DFunction.SIGNAL] = estimatedParams[Gaussian2DFunction.SIGNAL] / (2 * Math.PI *
							estimatedParams[Gaussian2DFunction.X_SD] * estimatedParams[Gaussian2DFunction.Y_SD]);
					params[j + Gaussian2DFunction.X_POSITION] = estimatedParams[Gaussian2DFunction.X_SD] -
							regionBounds.x;
					params[j + Gaussian2DFunction.Y_POSITION] = estimatedParams[Gaussian2DFunction.Y_SD] -
							regionBounds.y;
					params[j + Gaussian2DFunction.ANGLE] = estimatedParams[Gaussian2DFunction.ANGLE];
					params[j + Gaussian2DFunction.X_SD] = estimatedParams[Gaussian2DFunction.X_SD];
					params[j + Gaussian2DFunction.Y_SD] = estimatedParams[Gaussian2DFunction.Y_SD];
				}
				else
				{
					params[j + Gaussian2DFunction.SIGNAL] = spots[n2].intensity +
							((relativeIntensity) ? background : 0);
					params[j + Gaussian2DFunction.X_POSITION] = spots[n2].x - regionBounds.x;
					params[j + Gaussian2DFunction.Y_POSITION] = spots[n2].y - regionBounds.y;
				}
			}

			if (!relativeIntensity)
			{
				// In the case of uneven illumination the peaks may be below the background.
				// Check the heights are positive. Otherwise set it to zero and it will be re-estimated
				// in the fitter.
				for (int i = 0, j = 0; i <= neighbourCount; i++, j += parametersPerPeak)
				{
					if (params[j + Gaussian2DFunction.SIGNAL] < background)
					{
						background = 0;
						break;
					}
				}
			}

			// The fitted neighbours
			// TODO - Test which is better: (1) subtracting fitted peaks; or (2) including them
			// Initial tests show that Chi-squared is much lower when including them in the fit.

			double[] region = this.region;
			if (fittedNeighbourCount > 0)
			{
				// The fitted result will be relative to (0,0) and have a 0.5 pixel offset applied
				final double xOffset = regionBounds.x + 0.5;
				final double yOffset = regionBounds.y + 0.5;

				if (subtractFittedPeaks > 0)
				{
					// Subtract the already fitted peaks from the data. 
					// This will speed up evaluation of the fitting function.
					region = Arrays.copyOf(region, width * height);
					//Utils.display("Region", region, width, height);

					final double[] funcParams = new double[1 + parametersPerPeak * subtractFittedPeaks];
					for (int i = 0, j = 0; i < fittedNeighbourCount; i++)
					{
						if (!subtract[i])
							continue;
						PeakResult result = sliceResults.get(fittedNeighbourIndices[i]);
						// Copy Signal,Angle,Xpos,Ypos,Xwidth,Ywidth
						for (int k = 1; k <= parametersPerPeak; k++)
							funcParams[j + k] = result.params[k];
						// Adjust position relative to extracted region
						funcParams[j + Gaussian2DFunction.X_POSITION] -= xOffset;
						funcParams[j + Gaussian2DFunction.Y_POSITION] -= yOffset;
						j += parametersPerPeak;
					}

					GaussianFunction func = fitConfig.createGaussianFunction(subtractFittedPeaks, regionBounds.width,
							funcParams);
					func.initialise(funcParams);

					// Subtract fitted peaks
					for (int i = 0; i < region.length; i++)
					{
						region[i] -= func.eval(i);
					}

					//Utils.display("Region2", region, width, height);
				}

				// Add the details of the already fitted peaks
				for (int i = 0, j = (1 + neighbourCount) * parametersPerPeak; i < fittedNeighbourCount; i++)
				{
					if (subtract[i])
						continue;
					PeakResult result = sliceResults.get(fittedNeighbourIndices[i]);

					// Get the Amplitude (since the params array stores the signal)
					params[j + Gaussian2DFunction.SIGNAL] = result.getAmplitude();
					// Copy Angle,Xpos,Ypos,Xwidth,Ywidth
					for (int k = 2; k <= parametersPerPeak; k++)
						params[j + k] = result.params[k];

					// Add background to amplitude (required since the background will be re-estimated)
					params[j + Gaussian2DFunction.SIGNAL] += result.params[Gaussian2DFunction.BACKGROUND];
					// Adjust position relative to extracted region
					params[j + Gaussian2DFunction.X_POSITION] -= xOffset;
					params[j + Gaussian2DFunction.Y_POSITION] -= yOffset;
					j += parametersPerPeak;
				}
			}

			// Subtract the background from all estimated peak amplitudes
			if (background != 0)
			{
				params[Gaussian2DFunction.BACKGROUND] = background;
				for (int j = Gaussian2DFunction.SIGNAL; j < params.length; j += parametersPerPeak)
					params[j] -= background;
			}

			// Note that if the input XY positions are on the integer grid then the fitter will estimate
			// the position using a CoM estimate. To avoid this adjust the centres to be off-grid.
			for (int i = 0; i < npeaks; i++)
			{
				if ((int) params[i * parametersPerPeak + Gaussian2DFunction.X_POSITION] == params[i *
						parametersPerPeak + Gaussian2DFunction.X_POSITION])
					params[i * parametersPerPeak + Gaussian2DFunction.X_POSITION] += 0.001;
				if ((int) params[i * parametersPerPeak + Gaussian2DFunction.Y_POSITION] == params[i *
						parametersPerPeak + Gaussian2DFunction.Y_POSITION])
					params[i * parametersPerPeak + Gaussian2DFunction.Y_POSITION] += 0.001;
			}

			// -=-=-=-

			// Note: Some of the unfitted neighbours may be bad candidates.
			// Only validate the fitted neighbours and the current candidate.
			// The unfitted neighbours are allowed to fail.

			// Turn off validation of peaks
			final boolean fitValidation = fitConfig.isFitValidation();
			fitConfig.setFitValidation(false);

			// Increase the iterations for a multiple fit.
			final int maxIterations = fitConfig.getMaxIterations();
			final int maxEvaluations = fitConfig.getMaxFunctionEvaluations();
			fitConfig.setMaxIterations(
					maxIterations + maxIterations * (npeaks - 1) * ITERATION_INCREASE_FOR_MULTIPLE_PEAKS);
			fitConfig.setMaxIterations(
					maxEvaluations + maxEvaluations * (npeaks - 1) * EVALUATION_INCREASE_FOR_MULTIPLE_PEAKS);

			resultMulti = gf.fit(region, width, height, npeaks, params, true);

			printFitResults(resultMulti, region, width, height, npeaks, 0, gf.getIterations());

			// Restore
			fitConfig.setFitValidation(fitValidation);
			fitConfig.setMaxIterations(maxIterations);
			fitConfig.setMaxFunctionEvaluations(maxEvaluations);

			if (resultMulti.getStatus() != FitStatus.OK)
				// Probably failed to converge 
				return resultMulti;

			// Validate only the peaks we are interested in...
			if (fitValidation)
			{
				final double[] initialParams = resultMulti.getInitialParameters();
				final double[] fittedParams = resultMulti.getParameters();

				// The main candidate must be valid
				if (fitConfig.validatePeak(0, initialParams, fittedParams) != FitStatus.OK)
				{
					resultMulti.setStatus(fitConfig.getValidationResult(), fitConfig.getValidationData());
					return resultMulti;
				}

				// Already fitted peaks must be valid
				for (int n = neighbourCount + 1; n < npeaks; n++)
				{
					if (fitConfig.validatePeak(n, initialParams, fittedParams) != FitStatus.OK)
					{
						// Current this is a fail.
						// TODO - repeat fitting but subtract all the fitted peaks
						resultMulti.setStatus(fitConfig.getValidationResult(), fitConfig.getValidationData());
						return resultMulti;
					}
				}

				// The neighbours can be validated but are allowed to fail.
				GaussianOverlapAnalysis overlap = null;
				int flags = 0;
				for (int n = 1; n <= neighbourCount; n++)
				{
					if (fitConfig.validatePeak(n, initialParams, fittedParams) == FitStatus.OK)
					{
						// This can be stored as an estimate for the candidate
						if (storeEstimates)
							storeEstimate(neighbourIndices[n - 1], extractParams(fittedParams, n));
					}
					else
					{
						// The candidate failed
						if (fitConfig.getValidationResult() == FitStatus.COORDINATES_MOVED)
						{
							// If valid but has drifted then check for location within the image.
							final double shift = fitConfig.getCoordinateShift();
							fitConfig.setCoordinateShift(Double.MAX_VALUE);

							if (fitConfig.validatePeak(n, initialParams, fittedParams) == FitStatus.OK)
							{
								// This peak is OK, but has moved.
								// (TODO - Align all drifted spots with the neighbours as they may swap around)
								// Check for overlap with the main candidate
								if (overlap == null)
								{
									flags = fitConfig.getFunctionFlags();
									overlap = new GaussianOverlapAnalysis(flags, extractParams(fittedParams, n), 2);
								}
								overlap.add(flags, fittedParams, true);
							}

							fitConfig.setCoordinateShift(shift);
						}
					}
				}
				if (overlap != null)
				{
					double[] overlapData = overlap.getOverlapData();
					// Ensure the sum of the fitted function is greater than the sum of the 
					// overlap from neighbours
					if (overlapData[0] < overlapData[1])
					{
						resultMulti.setStatus(FitStatus.NEIGHBOUR_OVERLAP, overlapData);
						return resultMulti;
					}
				}
			}

			updateError(resultMulti);

			return resultMulti;
		}

		private void storeEstimate(int i, double[] params)
		{
			//if (!storeEstimates)
			//	return;

			FitWorker.this.storeEstimate(spots[i], i, params, regionBounds.x, regionBounds.x);
		}

		private double[] getEstimate(int i)
		{
			return (estimates[i] == null) ? null : estimates[i].params;
		}

		/**
		 * Update error to the coefficient of determination (so that the error describes how much of the data is
		 * encapsulated in the fit)
		 *
		 * @param result
		 *            the result multi
		 */
		private void updateError(FitResult result)
		{
			final double r2 = 1 - (gf.getFinalResidualSumOfSquares() / gf.getTotalSumOfSquares());
			result.setError(r2);
		}

		@SuppressWarnings("unused")
		public FitResult getResultMultiDoublet()
		{
			// TODO - Add this method
			// Fit a multi-fit. If all peaks are valid then subtract them from the image, apart from the central candidate,
			// then perform residuals analysis. Fit as a doublet where both spots must be closer to the central candidate
			// fitted position than any other candidate or existing fit result (including those from the multi-fit)

			// This can use new dedicated code since:
			// 1. we need to compute the IC correctly given the number of fitted parameters in the multi-fit
			// 2. the output doublet should be compared to all the results of the multi-fit (plus any fitted neighbours) 

			return null;
		}

		public FitResult getResultSingle()
		{
			if (resultSingle != null)
				return resultSingle;

			background = getMultiFittingBackground();
			double[] region = this.region;

			// Subtract all fitted neighbours from the region
			if (fittedNeighbourCount != 0)
			{
				// The fitted result will be relative to (0,0) and have a 0.5 pixel offset applied
				final double xOffset = regionBounds.x + 0.5;
				final double yOffset = regionBounds.y + 0.5;

				region = Arrays.copyOf(region, width * height);
				//Utils.display("Region", region, width, height);

				final double[] funcParams = new double[1 + parametersPerPeak * fittedNeighbourCount];
				for (int i = 0, j = 0; i < fittedNeighbourCount; i++)
				{
					PeakResult result = sliceResults.get(fittedNeighbourIndices[i]);
					// Copy Signal,Angle,Xpos,Ypos,Xwidth,Ywidth
					for (int k = 1; k <= parametersPerPeak; k++)
						funcParams[j + k] = result.params[k];
					// Adjust position relative to extracted region
					funcParams[j + Gaussian2DFunction.X_POSITION] -= xOffset;
					funcParams[j + Gaussian2DFunction.Y_POSITION] -= yOffset;
					j += parametersPerPeak;
				}

				GaussianFunction func = fitConfig.createGaussianFunction(fittedNeighbourCount, regionBounds.width,
						funcParams);
				func.initialise(funcParams);

				for (int i = 0; i < region.length; i++)
				{
					region[i] -= func.eval(i);
				}
			}

			// Store this as it is used in doublet fitting
			this.singleRegion = region;

			// TODO - See if this changes fitting accuracy/speed
			// Estimate the signal here instead of using the amplitude.
			// Do this when the fitting window covers enough of the Gaussian (e.g. 2.5xSD).
			boolean amplitudeEstimate = false;
			float signal = 0;
			if (estimateSignal)
			{
				double sum = 0;
				final int size = width * height;
				for (int i = size; i-- > 0;)
					sum += region[i];
				signal = (float) (sum - background * size);
			}
			if (signal <= 0)
			{
				amplitudeEstimate = true;
				signal = spots[n].intensity - ((relativeIntensity) ? 0 : background);
				if (signal < 0)
				{
					signal += background;
					background = 0;
				}
			}
			final double[] params = new double[] { background, signal, 0, x - regionBounds.x, y - regionBounds.y, 0,
					0 };

			// If there were unfitted neighbours then use off grid pixels to prevent re-estimate of CoM
			// since the CoM estimate will be skewed by the neighbours
			if (neighbourCount != 0)
			{
				params[Gaussian2DFunction.X_POSITION] += 0.001;
				params[Gaussian2DFunction.Y_POSITION] += 0.001;
			}

			resultSingle = gf.fit(region, width, height, 1, params, amplitudeEstimate);

			printFitResults(resultSingle, region, width, height, 1, -neighbours, gf.getIterations());

			updateError(resultSingle);

			return resultSingle;
		}

		public FitResult getResultDoublet(boolean force, MultiPathFilter filter)
		{
			if (computedDoublet)
				return resultDoublet;
			computedDoublet = true;

			// Ensure we have a single result
			getResultSingle();

			// Use the region from the single fit which had fitted peaks subtracted
			final double[] region = singleRegion;

			// Only compute quadrant improvement if fitted as a single peak. If fitting multiple peaks
			// then we would expect the quadrant to be skewed by the neighbours.
			if (force || canPerformQuadrantAnalysis(resultSingle, width, height))
			{
				computedDoubletQA = true;
				resultDoublet = quadrantAnalysis(spots, n, resultSingle, region, regionBounds, filter);
				if (resultDoublet != null)
				{
					updateError(resultDoublet);
				}
			}

			return resultDoublet;
		}
	}

	private void storeEstimate(Spot spot, int i, double[] params, double offsetx, double offsety)
	{
		// Add region offset
		params[Gaussian2DFunction.X_POSITION] += offsetx;
		params[Gaussian2DFunction.Y_POSITION] += offsety;
		// Compute distance to spot
		final double dx = spot.x - params[Gaussian2DFunction.X_POSITION];
		final double dy = spot.y - params[Gaussian2DFunction.Y_POSITION];
		final double d2 = dx * dx + dy * dy;
		if (estimates[i] == null || d2 < estimates[i].d2)
			estimates[i] = new Estimate(d2, params);
	}

	/**
	 * Extract parameters for the specified peak.
	 *
	 * @param params
	 *            the params
	 * @param n
	 *            the peak
	 * @return the extracted params
	 */
	private static double[] extractParams(double[] params, int n)
	{
		double[] newParams = new double[7];
		newParams[Gaussian2DFunction.BACKGROUND] = params[Gaussian2DFunction.BACKGROUND];
		for (int i = 0, j = n * 6 + Gaussian2DFunction.SIGNAL; i < 6; i++, j++)
		{
			newParams[i + Gaussian2DFunction.SIGNAL] = params[j];
		}
		return newParams;
	}

	/**
	 * Fits a 2D Gaussian to the given data. Fits all the specified peaks.
	 * <p>
	 * Data must be arranged in yx block order, i.e. height rows of width.
	 * <p>
	 * Can perform differential residue analysis on fit results to detect doublets. These are re-fitted and accepted if
	 * the residuals are significantly improved.
	 * 
	 * @param gf
	 * @param region
	 * @param regionBounds
	 * @param spots
	 *            All the maxima
	 * @param n
	 *            The maxima to fit
	 * @return Array containing the fitted curve data: The first value is the Background. The remaining values are
	 *         Amplitude, Angle, PosX, PosY, StdDevX, StdDevY for each fitted peak.
	 *         <p>
	 *         Null if no fit is possible.
	 */
	private FitResult fit(Gaussian2DFitter gf, double[] region, Rectangle regionBounds, Spot[] spots, int n)
	{
		// When multi-fit fails (of the candidate) then ensure that a single fit is attempted on the
		// neighbour-fit subtracted data + residuals analysis. This is done since
		// fitting multiple initial start points may be prone to error in estimates.

		// This makes the system a greedy single-pass algorithm that should 
		// fit an addition spot each time, or accumulate the fail count.

		// It should still work well on low density data and can gracefully fail HD data.

		// This is used to track the decision tree
		final FitType fitType = new FitType();

		SpotFitter spotFitter = new SpotFitter(gf, region, regionBounds, spots, n);

		// Attempt multi-fit if settings allow
		FitResult fitResult = spotFitter.getResultMulti();
		if (fitResult != null)
		{
			fitType.setNeighbours(true);
			if (fitResult.getStatus() == FitStatus.OK)
			{
				fitType.setNeighboursOK(true);

				// Extract the first fitted peak
				final int degreesOfFreedom = fitResult.getDegreesOfFreedom();
				double error = fitResult.getError();
				final double[] initialParameters = truncate(fitResult.getInitialParameters());
				final double[] parameters = truncate(fitResult.getParameters());
				final double[] parametersDev = (fitResult.getParameterStdDev() != null)
						? truncate(fitResult.getParameterStdDev()) : null;
				final int nPeaks = fitResult.getParameters().length / 6;
				final int offset = (fitConfig.isBackgroundFitting()) ? 1 : 0;
				final int nFittedParameters = offset + (fitResult.getNumberOfFittedParameters() - offset) / nPeaks;

				double r2 = 1 - (gf.getFinalResidualSumOfSquares() / gf.getTotalSumOfSquares());
				error = r2;

				fitResult = new FitResult(FitStatus.OK, degreesOfFreedom, error, initialParameters, parameters,
						parametersDev, 1, nFittedParameters, null, fitResult.getIterations(),
						fitResult.getEvaluations());

				// Note that if we have neighbours it is very possible we will have doublets
				// due to high density data. It makes sense to attempt a doublet fit here.

				// TODO - Doublet analysis here if all neighbour peaks were valid
			}
			else
			{
				// Fall back to single-fit
				fitResult = null;
			}
		}

		if (fitResult == null)
		{
			// Multi-fit failed, do a single fit
			fitResult = spotFitter.getResultSingle();
			if (fitResult.getStatus() == FitStatus.OK)
			{
				// Attempt doublet fit if settings allow
				FitResult newFitResult = spotFitter.getResultDoublet(false, null);
				if (spotFitter.computedDoubletQA)
				{
					fitType.setDoublet(true);
					// The doublet result is null if nothing was valid
					if (newFitResult != null && newFitResult.getStatus() == FitStatus.OK)
					{
						fitType.setDoubletOK(true);
						// The doublet already has the valid peak data extracted
						fitResult = newFitResult;
					}
				}
			}
		}

		if (logger != null)
		{
			switch (fitResult.getStatus())
			{
				case OK:
					// Show the shift, signal and width spread
					final double[] initialParams = fitResult.getInitialParameters();
					final double[] params = fitResult.getParameters();
					for (int i = 0, j = 0; i < fitResult.getNumberOfPeaks(); i++, j += 6)
					{
						logger.info("Fit OK [%d]. Shift = %.3f,%.3f : SNR = %.2f : Width = %.2f,%.2f", i,
								params[j + Gaussian2DFunction.X_POSITION] -
										initialParams[j + Gaussian2DFunction.X_POSITION],
								params[j + Gaussian2DFunction.Y_POSITION] -
										initialParams[j + Gaussian2DFunction.Y_POSITION],
								params[j + Gaussian2DFunction.SIGNAL] / noise,
								getFactor(params[j + Gaussian2DFunction.X_SD],
										initialParams[j + Gaussian2DFunction.X_SD]),
								getFactor(params[j + Gaussian2DFunction.Y_SD],
										initialParams[j + Gaussian2DFunction.Y_SD]));
					}
					break;

				case BAD_PARAMETERS:
					final double[] p = fitResult.getInitialParameters();
					logger.info("Bad parameters: %f,%f,%f,%f,%f,%f,%f", p[0], p[1], p[2], p[3], p[4], p[5], p[6]);
					break;

				default:
					logger.info(fitResult.getStatus().toString());
					break;
			}
		}

		if (fitResult.getStatus() == FitStatus.OK)
			fitType.setOK(true);
		add(fitType);
		return fitResult;
	}

	/**
	 * Fits a 2D Gaussian to the given data. Fits all possible paths and stores the multi-path result.
	 * <p>
	 * The best result is chosen using the provided MultiPathFilter and stored in the sliceResults array.
	 * 
	 * @param gf
	 * @param region
	 * @param regionBounds
	 * @param spots
	 *            All the maxima
	 * @param n
	 *            The maxima to fit
	 * @param filter
	 *            The filter for choosing the best result to add to the slice results
	 * @return The full multi-path result
	 */
	private MultiPathFitResult benchmarkFit(Gaussian2DFitter gf, double[] region, Rectangle regionBounds, Spot[] spots,
			int n, MultiPathFilter filter)
	{
		final SpotFitter spotFitter = new SpotFitter(gf, region, regionBounds, spots, n);
		// Disable estimate within the multi-fit as the filtering of peaks is probably disabled. 
		// We will do the estimates using the input filter.
		spotFitter.storeEstimates = false;

		final MultiPathFitResult result = new MultiPathFitResult();
		result.frame = slice;
		result.width = regionBounds.width;
		result.height = regionBounds.height;

		// Offsets to convert fit coordinates to the data frame
		final float offsetx = bounds.x + regionBounds.x + 0.5f;
		final float offsety = bounds.y + regionBounds.y + 0.5f;

		// Attempt multi-fit if settings allow
		FitResult fitResult = spotFitter.getResultMulti();
		if (fitResult != null)
		{
			result.multiFitResultStatus = fitResult.getStatus();
			if (fitResult.getStatus() == FitStatus.OK)
			{
				// Extract the peaks
				convertParameters(fitResult.getParameters());
				BasePreprocessedPeakResult[] multiFitResult = new BasePreprocessedPeakResult[fitResult
						.getNumberOfPeaks()];
				result.multiFitResult = multiFitResult;

				// Assume the first fitted peak is the new result
				multiFitResult[0] = createNewPreprocessedPeakResult(0, fitResult.getParameters(),
						fitResult.getInitialParameters(), offsetx, offsety);
				// Validate to check if we can use the candidates as estimates  
				boolean ok = filter.accept(multiFitResult[0]);

				// Add the existing results
				for (int i = neighbourCount + 1; i < multiFitResult.length; i++)
				{
					multiFitResult[i] = createExistingPreprocessedPeakResult(i, fitResult.getParameters(),
							fitResult.getInitialParameters(), offsetx, offsety);
					ok = ok && filter.accept(multiFitResult[i]);
				}

				// Then add the neighbours, storing estimates if they are good fits.
				// Note that the PreprocessedPeakResult will have corrected coordinates so
				// we subtract them and make relative to the fit region bounds.
				final float correctionx = regionBounds.x - offsetx;
				final float correctiony = regionBounds.y - offsety;
				for (int i = 1; i <= neighbourCount; i++)
				{
					multiFitResult[i] = createCandidatePreprocessedPeakResult(i, fitResult.getParameters(),
							fitResult.getInitialParameters(), offsetx, offsety);
					if (ok && filter.accept(multiFitResult[i]))
					{
						// If the multi-fit results are good then we can store the estimates using 
						// candidate results that pass the filter.
						final int ii = neighbourIndices[i - 1];
						storeEstimate(spots[ii], ii, multiFitResult[i].toGaussian2DParameters(), correctionx,
								correctiony);
					}
				}

				// Note that if we have neighbours it is very possible we will have doublets
				// due to high density data. It makes sense to attempt a doublet fit here.

				// TODO - Doublet analysis here if all neighbour peaks were valid
			}
		}

		// Do a single fit
		fitResult = spotFitter.getResultSingle();
		result.singleFitResultStatus = fitResult.getStatus();
		if (fitResult.getStatus() == FitStatus.OK)
		{
			// Extract the peaks
			convertParameters(fitResult.getParameters());
			result.singleFitResult = new BasePreprocessedPeakResult[1];
			result.singleFitResult[0] = createNewPreprocessedPeakResult(0, fitResult.getParameters(),
					fitResult.getInitialParameters(), offsetx, offsety);

			// Doublet fit
			// Force this if the residuals threshold is configured. The MultiPathFilter can
			// use the residuals QA score in its decision path.
			if (config.getResidualsThreshold() < 1)
			{
				fitResult = spotFitter.getResultDoublet(true, filter);

				// Store residuals analysis data
				result.singleQAScore = lastQAscore;

				// The doublet result is null if nothing was valid
				if (fitResult != null)
				{
					result.doubletFitResultStatus = fitResult.getStatus();
					if (fitResult.getStatus() == FitStatus.OK)

					{
						// Extract the peaks
						convertParameters(fitResult.getParameters());
						result.doubletFitResult = new BasePreprocessedPeakResult[fitResult.getNumberOfPeaks()];
						for (int i = 0; i < result.doubletFitResult.length; i++)
							result.doubletFitResult[i] = createNewPreprocessedPeakResult(i, fitResult.getParameters(),
									fitResult.getInitialParameters(), offsetx, offsety);
					}
				}
			}
		}

		// Pick the best result
		final PreprocessedPeakResult[] results = filter.accept(result);
		if (results != null)
		{
			// Add to the current slice results. Note that the PreprocessedPeakResult will have the bounds added already.
			final double[] params = new double[1 + results.length * 6];
			params[Gaussian2DFunction.BACKGROUND] = results[0].getBackground();
			for (int i = 0; i < results.length; i++)
			{
				final double[] p = ((BasePreprocessedPeakResult) results[i]).toGaussian2DParameters();
				// Convert back to the initial positions
				p[Gaussian2DFunction.X_POSITION] -= offsetx;
				p[Gaussian2DFunction.Y_POSITION] -= offsety;
				System.arraycopy(p, 1, params, 1 + i * 6, 6);
			}
			addResults(spots[n].x, spots[n].y, bounds, regionBounds, params, null, 0, 0, noise);
		}

		return result;
	}

	private BasePreprocessedPeakResult createNewPreprocessedPeakResult(int i, double[] parameters,
			double[] initialParameters, float offsetx, float offsety)
	{
		return createPreprocessedPeakResult(i, parameters, initialParameters, BasePreprocessedPeakResult.ResultType.NEW,
				offsetx, offsety);
	}

	private BasePreprocessedPeakResult createCandidatePreprocessedPeakResult(int i, double[] parameters,
			double[] initialParameters, float offsetx, float offsety)
	{
		return createPreprocessedPeakResult(i, parameters, initialParameters,
				BasePreprocessedPeakResult.ResultType.CANDIDATE, offsetx, offsety);
	}

	private BasePreprocessedPeakResult createExistingPreprocessedPeakResult(int i, double[] parameters,
			double[] initialParameters, float offsetx, float offsety)
	{
		return createPreprocessedPeakResult(i, parameters, initialParameters,
				BasePreprocessedPeakResult.ResultType.EXISTING, offsetx, offsety);
	}

	private BasePreprocessedPeakResult createPreprocessedPeakResult(int i, double[] parameters,
			double[] initialParameters, BasePreprocessedPeakResult.ResultType resultType, float offsetx, float offsety)
	{
		final int offset = i * 6;
		int frame = slice;
		float signal = (float) parameters[offset + Gaussian2DFunction.SIGNAL];
		float photons = (float) (parameters[offset + Gaussian2DFunction.SIGNAL] / fitConfig.getGain());
		float b = (float) parameters[Gaussian2DFunction.BACKGROUND];
		float angle = (float) parameters[offset + Gaussian2DFunction.ANGLE];
		float x = (float) parameters[offset + Gaussian2DFunction.X_POSITION] + offsetx;
		float y = (float) parameters[offset + Gaussian2DFunction.Y_POSITION] + offsety;
		float x0 = (float) initialParameters[offset + Gaussian2DFunction.X_POSITION] + offsetx;
		float y0 = (float) initialParameters[offset + Gaussian2DFunction.Y_POSITION] + offsety;
		float xsd = (float) parameters[offset + Gaussian2DFunction.X_SD];
		float ysd = (float) parameters[offset + Gaussian2DFunction.Y_SD];
		float xsd0 = (float) initialParameters[offset + Gaussian2DFunction.X_SD];
		float ysd0 = (float) initialParameters[offset + Gaussian2DFunction.Y_SD];
		double variance = fitConfig.getVariance(b, signal, (xsd + ysd) * 0.5);
		return new BasePreprocessedPeakResult(frame, signal, photons, i, b, angle, x, y, x0, y0, xsd, ysd, xsd0, ysd0,
				variance, resultType);
	}

	/**
	 * @return The background estimate when fitting a single peak
	 */
	private float getSingleFittingBackground()
	{
		float background;
		if (useFittedBackground && !sliceResults.isEmpty())
		{
			// Use the average background from all results
			background = (float) (fittedBackground / sliceResults.size());
		}
		else
		{
			// Initial guess using image mean
			background = this.background;
		}
		return background;
	}

	private void printFitResults(FitResult fitResult, double[] region, int width, int height, int npeaks, int doublet,
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
		//		double ic = getInformationCriterion(gf.getFinalResidualSumOfSquares(), length,
		//				fitResult.getInitialParameters().length);
		//
		//		System.out.printf("[%dx%d] p %d i %d d %d : SS %f : %f (%s). %f\n", width, height, npeaks, iterations,
		//				doublet, gf.getTotalSumOfSquares(), 
		//				gf.getFinalResidualSumOfSquares(), fitResult.getStatus().toString(), ic);
	}

	private void printFitResults(FitResult fitResult, double[] region, int width, int height, int npeaks, int doublet,
			int iterations, double ICc1, double ICc2)
	{
		// **********
		// Comment out for production code
		// **********

		//		if (fitResult.getStatus() == FitStatus.BAD_PARAMETERS)
		//			return;
		//
		//		System.out.printf("[%dx%d] p %d i %d d %d : SS %f : %f : %f (%s). %f -> %f\n", width, height, npeaks,
		//				iterations, doublet, gf.getTotalSumOfSquares(), gf.getInitialResdiualSumOfSquares(),
		//				gf.getFinalResdiualSumOfSquares(), fitResult.getStatus().toString(), ICc1, ICc2);
	}

	/**
	 * Get the relative change factor between f and g
	 * 
	 * @param f
	 * @param g
	 * @return The relative change factor (negative if g is bigger than f)
	 */
	private static double getFactor(double f, double g)
	{
		if (f > g)
			return f / g;
		return -g / f;
	}

	private static double[] truncate(double[] array)
	{
		double[] newArray = new double[7];
		System.arraycopy(array, 0, newArray, 0, 7);
		return newArray;
	}

	private void resetNeighbours()
	{
		neighbourCount = 0;
		fittedNeighbourCount = 0;
		spotNeighbours = null;
		peakNeighbours = null;
	}

	private Spot[] findSpotNeighbours(Spot[] spots, int n)
	{
		if (spotNeighbours == null)
		{
			// Using the neighbour grid 
			spotNeighbours = gridManager.getSpotNeighbours(spots[n].x, spots[n].y, n + 1);
			// Ensure they are sorted in ID order
			Arrays.sort(spotNeighbours);
			Collections.reverse(Arrays.asList(spotNeighbours));
		}
		return spotNeighbours;
	}

	private PeakResult[] findPeakNeighbours(Spot[] spots, int n)
	{
		if (peakNeighbours == null)
		{
			// Using the neighbour grid 
			peakNeighbours = gridManager.getPeakResultNeighbours(spots[n].x, spots[n].y);
		}
		return peakNeighbours;
	}

	/**
	 * Search for any peak within a set height of the specified peak that is within the search region bounds
	 * 
	 * @param regionBounds
	 * @param n
	 * @param x
	 * @param y
	 * @param spots
	 * @return The number of neighbours
	 */
	private int findNeighbours(Rectangle regionBounds, int n, int x, int y, Spot[] spots)
	{
		int xmin = regionBounds.x;
		int xmax = xmin + regionBounds.width - 1;
		int ymin = regionBounds.y;
		int ymax = ymin + regionBounds.height - 1;

		float heightThreshold;
		if (relativeIntensity)
		{
			// No background when spot filter has relative intensity
			heightThreshold = (float) (spots[n].intensity * config.getNeighbourHeightThreshold());
		}
		else
		{
			if (spots[n].intensity < background)
				heightThreshold = spots[n].intensity;
			else
				heightThreshold = (float) ((spots[n].intensity - background) * config.getNeighbourHeightThreshold() +
						background);
		}

		// Check all maxima that are lower than this
		neighbourCount = 0;

		// Using the neighbour grid 
		final Spot[] spotNeighbours = findSpotNeighbours(spots, n);
		for (int i = 0; i < spotNeighbours.length; i++)
		{
			final int id = (int) spotNeighbours[i].getScore();
			if (canIgnore(spots[id].x, spots[id].y, xmin, xmax, ymin, ymax, spots[id].intensity, heightThreshold))
				continue;
			neighbourIndices[neighbourCount++] = id;
		}

		// Processing all lower spots.
		// XXX - This is the old code and is here as a check. Remove this.
		int c = 0;
		for (int i = n + 1; i < spots.length; i++)
		{
			if (canIgnore(spots[i].x, spots[i].y, xmin, xmax, ymin, ymax, spots[i].intensity, heightThreshold))
				continue;
			//neighbourIndices[c++] = i;
			if (neighbourIndices[c++] != i)
				throw new RuntimeException("invalid grid neighbours");
		}

		// Check all existing maxima. 

		fittedNeighbourCount = 0;
		if (!sliceResults.isEmpty())
		{
			// Since these will be higher than the current peak it is prudent to extend the range that should be considered.
			// Use 2x the configured peak standard deviation.
			//final float range = 2f *
			//		(float) FastMath.max(fitConfig.getInitialPeakStdDev0(), fitConfig.getInitialPeakStdDev1());
			//final float xmin2 = regionBounds.x - range;
			//final float xmax2 = regionBounds.x + regionBounds.width + range;
			//final float ymin2 = regionBounds.y - range;
			//final float ymax2 = regionBounds.y + regionBounds.height + range;
			//for (int i = 0; i < sliceResults.size(); i++)
			//{
			//	final PeakResult result = sliceResults.get(i);
			//	// No height threshold check as this is a validated peak
			//	if (canIgnore(result.getXPosition(), result.getYPosition(), xmin2, xmax2, ymin2, ymax2))
			//		continue;
			//	fittedNeighbourIndices[fittedNeighbourCount++] = i;
			//}

			// Note: A smarter filter would be to compute the bounding rectangle of each fitted result and see if it 
			// overlaps the target region. This would involve overlap analysis
			final double x0min = regionBounds.x;
			final double y0min = regionBounds.y;
			final double x0max = regionBounds.x + regionBounds.width;
			final double y0max = regionBounds.y + regionBounds.height;

			final PeakResult[] peakNeighbours = findPeakNeighbours(spots, n);

			for (int i = 0; i < peakNeighbours.length; i++)
			{
				final PeakResult result = peakNeighbours[i];
				// No height threshold check as this is a validated peak
				final double xw = 2 * result.getXSD();
				final double yw = 2 * result.getYSD();
				if (intersects(x0min, y0min, x0max, y0max, result.getXPosition() - xw, result.getYPosition() - yw,
						result.getXPosition() + xw, result.getYPosition() + yw))
					continue;
				fittedNeighbourIndices[fittedNeighbourCount++] = i;
			}
		}

		return neighbourCount + fittedNeighbourCount;
	}

	/**
	 * Copied from java.awt.geom.Rectangle2D and modified assuming width and height is non-zero
	 * 
	 * @param r1
	 * @param r2
	 * @return true if they intersect
	 */
	public boolean intersects(double x0min, double y0min, double x0max, double y0max, double x1min, double y1min,
			double x1max, double y1max)
	{
		return (x1max > x0min && y1max > y0min && x1min < x0max && y1min < y0max);
	}

	private boolean canIgnore(int x, int y, int xmin, int xmax, int ymin, int ymax, float height, float heightThreshold)
	{
		return (x < xmin || x > xmax || y < ymin || y > ymax || height < heightThreshold);
	}

	@SuppressWarnings("unused")
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
			// Check why it was a bad fit. 

			// If it due to width divergence then 
			// check the width is reasonable given the size of the fitted region.
			if (fitResult.getStatus() == FitStatus.WIDTH_DIVERGED)
			{
				final double[] params = fitResult.getParameters();
				final double regionSize = FastMath.max(width, height) * 0.5;
				//int tmpSmooth = (int) FastMath.max(smooth, 1);
				//float regionSize = 2 * tmpSmooth + 1;
				if ((params[Gaussian2DFunction.X_SD] > 0 && params[Gaussian2DFunction.X_SD] < regionSize) ||
						(params[Gaussian2DFunction.Y_SD] > 0 && params[Gaussian2DFunction.Y_SD] < regionSize))
					return true;
			}

			// If moved then it could be a close neighbour ...
			if (fitResult.getStatus() == FitStatus.COORDINATES_MOVED)
			{

			}
		}
		return false;
	}

	private double lastQAscore;

	private FitResult quadrantAnalysis(Spot[] spots, int candidate, FitResult fitResult, double[] region,
			Rectangle regionBounds, MultiPathFilter filter)
	{
		// Perform quadrant analysis as per rapidSTORM:

		final int width = regionBounds.width;
		final int height = regionBounds.height;
		final double[] params = fitResult.getParameters();
		// Use rounding since the fit coords are not yet offset by 0.5 pixel to centre them
		final int cx = (int) Math.round(params[Gaussian2DFunction.X_POSITION]);
		final int cy = (int) Math.round(params[Gaussian2DFunction.Y_POSITION]);
		final double[] residuals = gf.getResiduals();

		lastQAscore = 2;

		QuadrantAnalysis qa = new QuadrantAnalysis();
		if (!qa.quadrantAnalysis(residuals, width, height, cx, cy))
			return null;

		lastQAscore = qa.score;

		if (logger != null)
			logger.info("Residue analysis = %f (%d,%d)", qa.score, qa.vector[0], qa.vector[1]);

		// If differential residue analysis indicates a doublet then re-fit as two spots.
		if (qa.score < config.getResidualsThreshold())
			return null;

		if (logger != null)
			logger.info("Computing 2-kernel model");

		if (!qa.computeDoubletCentres(width, height, cx, cy, params[Gaussian2DFunction.X_SD],
				params[Gaussian2DFunction.Y_SD]))
			return null;

		// TODO - Locate the 2 new centres by moving out into the quadrant defined by the vector
		// and finding the maxima on the original image data.

		// -+-+-
		// Estimate params using the single fitted peak
		// -+-+-
		final double[] doubletParams = new double[1 + 2 * 6];

		// Note: Quadrant analysis sets the positions using 0.5,0.5 as the centre of the pixel.
		// The fitting does not, so subtract 0.5.

		doubletParams[Gaussian2DFunction.BACKGROUND] = params[Gaussian2DFunction.BACKGROUND];
		doubletParams[Gaussian2DFunction.SIGNAL] = params[Gaussian2DFunction.SIGNAL] * 0.5;
		doubletParams[Gaussian2DFunction.X_POSITION] = (float) (qa.x1 - 0.5);
		doubletParams[Gaussian2DFunction.Y_POSITION] = (float) (qa.y1 - 0.5);
		doubletParams[6 + Gaussian2DFunction.SIGNAL] = params[Gaussian2DFunction.SIGNAL] * 0.5;
		doubletParams[6 + Gaussian2DFunction.X_POSITION] = (float) (qa.x2 - 0.5);
		doubletParams[6 + Gaussian2DFunction.Y_POSITION] = (float) (qa.y2 - 0.5);
		// -+-+-

		// Store the residual sum-of-squares from the previous fit
		final double singleSumOfSquares = gf.getFinalResidualSumOfSquares();
		final double singleValue = gf.getValue();

		// Disable checking of position movement since we guessed the 2-peak location.
		// Increase the iterations level then reset afterwards.

		// TODO - Should width and signal validation be disabled too?
		double shift = fitConfig.getCoordinateShift();
		final int maxIterations = fitConfig.getMaxIterations();
		final int maxEvaluations = fitConfig.getMaxFunctionEvaluations();

		fitConfig.setCoordinateShift(FastMath.min(width, height));
		fitConfig.setMaxIterations(maxIterations * ITERATION_INCREASE_FOR_DOUBLETS);
		fitConfig.setMaxFunctionEvaluations(maxEvaluations * FitWorker.EVALUATION_INCREASE_FOR_DOUBLETS);

		//FitResult newFitResult = gf.fit(region, width, height, peaks, heights);
		gf.setComputeResiduals(false);
		final FitResult newFitResult = gf.fit(region, width, height, 2, doubletParams, false);
		gf.setComputeResiduals(true);

		fitConfig.setCoordinateShift(shift);
		fitConfig.setMaxIterations(maxIterations);
		fitConfig.setMaxFunctionEvaluations(maxEvaluations);

		// Q. Allow bad fits since these are due to peak width divergence or low signal for all peaks?
		if (newFitResult.getStatus() == FitStatus.OK) // || newFitResult.getResult() == Result.BAD_FIT)
		{
			// Adjusted Coefficient of determination is not good for non-linear models. Use the 
			// Bayesian Information Criterion (BIC):

			final double doubleSumOfSquares = gf.getFinalResidualSumOfSquares();

			final int length = width * height;
			//double SStotal = gf.getTotalSumOfSquares();
			//double adjustedR2 = getAdjustedCoefficientOfDetermination(singleSumOfSquares, SStotal, length,
			//		fitResult.getInitialParameters().length);
			//double newAdjustedR2 = getAdjustedCoefficientOfDetermination(doubleSumOfSquares, SStotal, length,
			//		newFitResult.getInitialParameters().length);

			final double ic1, ic2;

			// Get the likelihood for the fit.
			if (fitConfig.getFitSolver() == FitSolver.MLE && fitConfig.isModelCamera())
			{
				// This is computed directly by the maximum likelihood estimator. 
				// The MLE is only good if we are modelling the camera noise. 
				// The MLE put out by the Poisson model is not better than using the IC from the fit residuals.
				final double doubleValue = gf.getValue();
				ic1 = Maths.getBayesianInformationCriterion(singleValue, length,
						fitResult.getNumberOfFittedParameters());
				ic2 = Maths.getBayesianInformationCriterion(doubleValue, length,
						newFitResult.getNumberOfFittedParameters());
				if (logger != null)
					logger.info("Model improvement - Sum-of-squares, MLE (IC) : %f, %f (%f) => %f, %f (%f) : %f",
							singleSumOfSquares, singleValue, ic1, doubleSumOfSquares, doubleValue, ic2, ic1 - ic2);
			}
			else
			{
				// If using the least squares estimator then we can get the log likelihood from an approximation
				ic1 = Maths.getBayesianInformationCriterionFromResiduals(singleSumOfSquares, length,
						fitResult.getNumberOfFittedParameters());
				ic2 = Maths.getBayesianInformationCriterionFromResiduals(doubleSumOfSquares, length,
						newFitResult.getNumberOfFittedParameters());
				if (logger != null)
					logger.info("Model improvement - Sum-of-squares (IC) : %f (%f) => %f (%f) : %f", singleSumOfSquares,
							ic1, doubleSumOfSquares, ic2, ic1 - ic2);
			}

			if (logger2 != null)
			{
				double[] peakParams = newFitResult.getParameters();
				if (peakParams != null)
				{
					peakParams = Arrays.copyOf(peakParams, peakParams.length);
					int npeaks = peakParams.length / 6;
					for (int i = 0; i < npeaks; i++)
					{
						peakParams[i * 6 + Gaussian2DFunction.X_POSITION] += 0.5 + regionBounds.x;
						peakParams[i * 6 + Gaussian2DFunction.Y_POSITION] += 0.5 + regionBounds.y;
					}
				}
				String msg = String.format("Doublet %d [%d,%d] %s (%s) [%f -> %f] SS [%f -> %f] IC [%f -> %f] = %s\n",
						slice, cx + bounds.x + regionBounds.x, cy + bounds.y + regionBounds.y, newFitResult.getStatus(),
						newFitResult.getStatusData(), singleValue, gf.getValue(), singleSumOfSquares,
						doubleSumOfSquares, ic1, ic2, Arrays.toString(peakParams));
				logger2.debug(msg);
			}

			// Check if the predictive power of the model is better with two peaks:
			// IC should be lower
			// Adjusted R2 should be higher
			if (ic2 > ic1)
			//if (newAdjustedR2 < adjustedR2)
			{
				printFitResults(newFitResult, region, width, height, 2, 1, gf.getIterations(), ic1, ic2);
				return null;
			}

			// TODO - get the distance of each new centre from the original centre
			// If the shift is too far (e.g. half the distance to the edge), the centre must be in the correct
			// quadrant. Then check if there is an unfit candidate spot closer than the current candidate.
			// This represents drift out to fit another spot that will be fit later.

			// Check if either coordinate is outside the fitted region.
			final double[] newParams = newFitResult.getParameters();
			for (int n = 0; n < 2; n++)
			{
				// Note that during processing the data is assumed to refer to the top-left
				// corner of the pixel. The coordinates should be represented in the middle of the pixel 
				// so add a 0.5 shift to the coordinates.
				final double xpos = newParams[Gaussian2DFunction.X_POSITION + n * 6] + 0.5;
				final double ypos = newParams[Gaussian2DFunction.Y_POSITION + n * 6] + 0.5;

				// Set the bounds using the expected HWHM
				final double hwhm = Maths.max(fitConfig.getInitialPeakStdDev0(), fitConfig.getInitialPeakStdDev1(), 1) *
						GaussianFunction.SD_TO_HWHM_FACTOR;

				if (xpos < -hwhm || xpos > regionBounds.width + hwhm || ypos < -hwhm ||
						ypos > regionBounds.height + hwhm)
				{
					if (logger != null)
						logger.info("Fitted coordinates too far outside the fitted region (x %g || y %g) in %dx%d",
								xpos, ypos, regionBounds.width, regionBounds.height);
					printFitResults(newFitResult, region, width, height, 2, 1, gf.getIterations(), ic1, ic2);
					return null;
				}
			}

			// Check the distance of the peaks to the centre of the region. It is possible that a second peak
			// at the edge of the region has been fitted (note that no coordinate shift check was performed).

			if (shift == 0 || shift == Double.POSITIVE_INFINITY)
			{
				// Allow the shift to span half of the fitted window.
				shift = 0.5 * FastMath.min(regionBounds.width, regionBounds.height);
			}

			// Set an upper limit on the shift that is not too far outside the fit window
			final double maxShiftX, maxShiftY;
			final double factor = Gaussian2DFunction.SD_TO_HWHM_FACTOR;
			if (fitConfig.isWidth0Fitting())
			{
				// Add the fitted standard deviation to the allowed shift
				maxShiftX = regionBounds.width * 0.5 + factor * params[Gaussian2DFunction.X_SD];
				maxShiftY = regionBounds.height * 0.5 + factor * params[Gaussian2DFunction.Y_SD];
			}
			else
			{
				// Add the configured standard deviation to the allowed shift
				maxShiftX = regionBounds.width * 0.5 + factor * fitConfig.getInitialPeakStdDev0();
				maxShiftY = regionBounds.height * 0.5 + factor * fitConfig.getInitialPeakStdDev1();
			}

			int[] position = new int[2];
			int nPeaks = 0;
			NEXT_PEAK: for (int n = 0; n < 2; n++)
			{
				// Exit this loop if we found an error and we are not logging
				//if (logger == null && n != nPeaks)
				//	break;

				final double xShift = newParams[Gaussian2DFunction.X_POSITION + n * 6] -
						params[Gaussian2DFunction.X_POSITION];
				final double yShift = newParams[Gaussian2DFunction.Y_POSITION + n * 6] -
						params[Gaussian2DFunction.Y_POSITION];
				if (Math.abs(xShift) > maxShiftX || Math.abs(yShift) > maxShiftY)
				{
					if (logger != null)
					{
						logger.info("Bad peak %d: Fitted coordinates moved outside fit region (x=%g,y=%g)", n, xShift,
								yShift);
					}
					continue;
				}
				if (Math.abs(xShift) > shift || Math.abs(yShift) > shift)
				{
					// Allow large shifts if they are along the vector
					final double a = QuadrantAnalysis.getAngle(qa.vector, new double[] { xShift, yShift });
					// Check the domain is OK (the angle is in radians). 
					// Allow up to a 45 degree difference to show the shift is along the vector
					if (a > 0.785398 && a < 2.356194)
					{
						if (logger != null)
						{
							logger.info("Bad peak %d: Fitted coordinates moved into wrong quadrant (x=%g,y=%g,a=%f)", n,
									xShift, yShift, a * 57.29578);
						}
						continue;
					}
				}

				// Spots from the doublet must be closer to the single fit than any other spots.
				// This is true for candidates we have yet to fit or already fitted candidates.

				// Distance to current candidate fitted as a single
				final double d2 = xShift * xShift + yShift * yShift;

				// Check if there are any candidates closer than the current candidate with a 
				// fit window that contains this spot.
				// This represents drift out to fit another spot that will be fit later.

				// Position - do not add 0.5 pixel offset to allow distance comparison to integer candidate positions.
				float fcx2 = (float) (regionBounds.x + newParams[Gaussian2DFunction.X_POSITION + n * 6]);
				float fcy2 = (float) (regionBounds.y + newParams[Gaussian2DFunction.Y_POSITION + n * 6]);
				// Add offset to know the pixel this fit in within
				final int cx2 = (int) (fcx2 + 0.5);
				final int cy2 = (int) (fcy2 + 0.5);
				final int xmin = cx2 - fitting;
				final int xmax = cx2 + fitting;
				final int ymin = cy2 - fitting;
				final int ymax = cy2 + fitting;
				final Spot[] spotNeighbours = findSpotNeighbours(spots, candidate);
				for (int j = 0; j < spotNeighbours.length; j++)
				{
					final int i = (int) spotNeighbours[j].getScore();
					//for (int i = candidate + 1; i < spots.length; i++)
					//{
					if (
					//i == candidate || 
					spots[i].x < xmin || spots[i].x > xmax || spots[i].y < ymin || spots[i].y > ymax)
						continue;
					if (d2 > distance2(fcx2, fcy2, spots[i]))
					{
						if (logger != null)
						{
							logger.info(
									"Bad peak %d: Fitted coordinates moved closer to another candidate (%d,%d : x=%.1f,y=%.1f : %d,%d)",
									n, spots[candidate].x, spots[candidate].y, fcx2 + 0.5f, fcy2 + 0.5f, spots[i].x,
									spots[i].y);
						}

						// Store the estimate. This has passed filtering in the FitConfig object.
						// It must optionally pass an additional filter 
						if (filter != null)
						{
							if (filter.accept(fitConfig.createPreprocessedPeakResult(n,
									newFitResult.getInitialParameters(), newParams)))
								;
							storeEstimate(spots[i], i, extractParams(newParams, n), regionBounds.x, regionBounds.x);
						}
						else
							storeEstimate(spots[i], i, extractParams(newParams, n), regionBounds.x, regionBounds.x);

						// There is another candidate to be fit later that is closer
						continue NEXT_PEAK;
					}
				}

				// Note: We cannot ignore already fitted spots as they may be on the edge of the fit window and 
				// so the result does not quite enter the duplicate distance.
				fcx2 += 0.5f;
				fcy2 += 0.5f;
				final float fxmin = fcx2 - fitting;
				final float fxmax = fcx2 + fitting;
				final float fymin = fcy2 - fitting;
				final float fymax = fcy2 + fitting;
				final PeakResult[] peakNeighbours = findPeakNeighbours(spots, candidate);
				for (PeakResult result : peakNeighbours) //sliceResults)
				{
					if (result.getXPosition() < fxmin || result.getXPosition() > fxmax ||
							result.getYPosition() < fymin || result.getYPosition() > fymax)
						continue;
					if (d2 > distance2(fcx2, fcy2, result.params))
					{
						if (logger != null)
						{
							logger.info(
									"Bad peak %d: Fitted coordinates moved closer to another result (%d,%d : x=%.1f,y=%.1f : %.1f,%.1f)",
									n, spots[candidate].x, spots[candidate].y, fcx2, fcy2, result.getXPosition(),
									result.getYPosition());
						}
						// There is another fitted result that is closer.
						// This indicates that the previously fitted result should be repeated. Currently
						// this is ignored and the current fit discarded. 
						continue NEXT_PEAK;
					}
				}

				position[nPeaks++] = n;
			}

			// Debug print here so we can see the number of valid peaks in the doublet
			printFitResults(newFitResult, region, width, height, 2, 3 + nPeaks, gf.getIterations(), ic1, ic2);

			if (nPeaks == 0)
			{
				return null;
			}

			// This code ensures the doublet is tightly located around the original single fit centre
			//if (nPeaks != 2)
			//{
			//	return null;
			//}

			// The code below allows any valid result from the doublet ...

			// Copy the OK peaks into a new result
			final double[] newInitialParams = newFitResult.getInitialParameters();
			final double[] newParamStdDev = newFitResult.getParameterStdDev();

			final double[] okInitialParams = new double[1 + nPeaks * 6];
			final double[] okParams = new double[1 + nPeaks * 6];
			final double[] okParamStdDev = new double[1 + nPeaks * 6];

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

			final int offset = (fitConfig.isBackgroundFitting()) ? 1 : 0;
			final int nFittedParameters = offset + (newFitResult.getNumberOfFittedParameters() - offset) / nPeaks;

			double error = newFitResult.getError();
			final double r2 = 1 - (gf.getFinalResidualSumOfSquares() / gf.getTotalSumOfSquares());
			error = r2;

			return new FitResult(newFitResult.getStatus(), newFitResult.getDegreesOfFreedom(), error, okInitialParams,
					okParams, okParamStdDev, nPeaks, nFittedParameters, newFitResult.getStatusData(),
					fitResult.getIterations(), fitResult.getEvaluations());
		}
		else
		{
			if (logger != null)
				logger.info("Unable to fit 2-kernel model : %s", newFitResult.getStatus());

			if (logger2 != null)
			{
				double[] peakParams = newFitResult.getParameters();
				if (peakParams != null)
				{
					peakParams = Arrays.copyOf(peakParams, peakParams.length);
					int npeaks = peakParams.length / 6;
					for (int i = 0; i < npeaks; i++)
					{
						peakParams[i * 6 + Gaussian2DFunction.X_POSITION] += 0.5 + regionBounds.x;
						peakParams[i * 6 + Gaussian2DFunction.Y_POSITION] += 0.5 + regionBounds.y;
					}
				}
				String msg = String.format("Doublet %d [%d,%d] %s (%s) = %s\n", slice, cx + bounds.x + regionBounds.x,
						cy + bounds.y + regionBounds.y, newFitResult.getStatus(), newFitResult.getStatusData(),
						Arrays.toString(peakParams));
				logger2.debug(msg);
			}

		}

		return null;
	}

	private float distance2(float cx, float cy, Spot spot)
	{
		final float dx = cx - spot.x;
		final float dy = cy - spot.y;
		return dx * dx + dy * dy;
	}

	private float distance2(float cx, float cy, float[] params)
	{
		final float dx = cx - params[Gaussian2DFunction.X_POSITION];
		final float dy = cy - params[Gaussian2DFunction.Y_POSITION];
		return dx * dx + dy * dy;
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

	/**
	 * @return Use the average fitted background as the background estimate for new fits
	 */
	public boolean isUseFittedBackground()
	{
		return useFittedBackground;
	}

	/**
	 * Set to true to use the average fitted background as the background estimate for new fits. The default is to use
	 * the image average as the background for all fits.
	 * 
	 * @param Use
	 *            the average fitted background as the background estimate for new fits
	 */
	public void setUseFittedBackground(boolean useFittedBackground)
	{
		this.useFittedBackground = useFittedBackground;
	}

	/**
	 * Sets the 2nd logger instance. This can be used for capturing debugging information.
	 *
	 * @param logger
	 *            the new logger
	 */
	public void setLogger2(Logger logger)
	{
		this.logger2 = logger;
	}

	/**
	 * Set the counter. This can be used to count the type of fitting process that was performed.
	 * 
	 * @param counter
	 *            The counter
	 */
	public void setCounter(FitTypeCounter counter)
	{
		this.counter = counter;
	}

	private void add(FitType fitType)
	{
		if (counter != null)
			counter.add(fitType);
	}
}
