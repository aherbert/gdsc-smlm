package gdsc.smlm.engine;

import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.concurrent.BlockingQueue;

import org.apache.commons.math3.util.FastMath;

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
import gdsc.smlm.engine.ResultGridManager.CandidateList;
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
import gdsc.smlm.results.filter.BasePreprocessedPeakResult.ResultType;
import gdsc.smlm.results.filter.IMultiPathFitResults;
import gdsc.smlm.results.filter.MultiFilter2;
import gdsc.smlm.results.filter.MultiPathFilter;
import gdsc.smlm.results.filter.MultiPathFilter.SelectedResult;
import gdsc.smlm.results.filter.MultiPathFilter.SelectedResultStore;
import gdsc.smlm.results.filter.MultiPathFitResult;
import gdsc.smlm.results.filter.PreprocessedPeakResult;

/**
 * Fits local maxima using a 2D Gaussian.
 */
public class FitWorker implements Runnable, IMultiPathFitResults, SelectedResultStore
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
	private Rectangle dataBounds;
	//private Rectangle regionBounds;
	private int border, borderLimitX, borderLimitY;
	private FitJob job;
	private float[] data;
	private boolean relativeIntensity;
	private float noise, background;
	private final boolean calculateNoise;
	private boolean estimateSignal;
	private ResultGridManager gridManager = null;
	private PeakResult[] peakNeighbours = null;
	// Contains the index in the list of maxima for any neighbours
	private int neighbourCount = 0;
	private int[] neighbourIndices = null;
	// Contains the index in the list of fitted results for any neighbours 
	private int fittedNeighbourCount = 0;
	private int[] fittedNeighbourIndices = null;
	private float duplicateDistance2;

	private volatile boolean finished = false;

	static int WORKER_ID = 0;
	int workerId;

	private static byte FILTER_RANK_MINIMAL = (byte) 0;
	private static byte FILTER_RANK_PRIMARY = (byte) 1;

	/**
	 * Store an estimate for a spot candidate. This may be aquired during multi fitting of neighbours.
	 */
	private class Estimate
	{
		final double[] params;
		final byte filterRank;
		final double d2;
		final double precision;

		public Estimate(double[] params, byte filterRank, double d2, double precision)
		{
			this.params = params;
			this.filterRank = filterRank;
			this.d2 = d2;
			this.precision = precision;
		}

		boolean isWeaker(byte filterRank, double d2, double precision)
		{
			if (this.filterRank < filterRank)
				return true;
			// The spot must be close to the estimate.
			// Note that if fitting uses a bounded fitter the estimates
			// should always be within 2 (1 pixel max in each dimension).
			// This check ensure that unbounded fitters store close estimates.
			if (d2 > 2)
			{
				// This is not very close so make sure the closest estimate is stored
				return (this.d2 > d2);
			}
			// If it is close enough then we use to fit precision.			
			if (this.precision > precision)
				return true;
			return false;
		}
	}

	private Estimate[] estimates = null;

	private Candidate[] candidates = null;
	private CandidateList neighbours = null;

	/**
	 * Instantiates a new fit worker.
	 * <p>
	 * Note that if the fit configuration has fit validation enabled then the initial fit results will be validated
	 * using
	 * only the basic filtering setting of the fit configuration. The use of the smart filter will be disabled. Once all
	 * results have passed the basic validation the results are then filtered again using the IDirectFilter
	 * implementation
	 * of the fit configuration. This will use a configured smart filter if present.
	 *
	 * @param config
	 *            the configuration
	 * @param results
	 *            the results
	 * @param jobs
	 *            the jobs
	 */
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

		// Disable the use of the direct filter within the FitConfiguration validate method. 
		// This allows validate() to be used for basic filtering of all fit results (e.g. using the fit region bounds).
		// The validation of each result will be performed by the FitConfiguration implementation
		// of the IDirectFilter interface. This may involve the DirectFilter object.
		fitConfig.setSmartFilter(false);

		workerId = WORKER_ID++;
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
		this.job = job;
		this.slice = job.slice;

		// Used for debugging
		//if (logger == null) logger = new gdsc.fitting.logging.ConsoleLogger();

		// Crop to the ROI
		dataBounds = job.bounds;
		final int width = dataBounds.width;
		final int height = dataBounds.height;
		borderLimitX = width - border;
		borderLimitY = height - border;
		data = job.data;

		FitParameters params = job.getFitParameters();
		this.endT = (params != null) ? params.endT : -1;

		candidates = indentifySpots(job, width, height, params);

		if (candidates.length == 0)
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

		final ImageExtractor ie = new ImageExtractor(data, width, height);
		double[] region = null;

		if (params != null && params.fitTask == FitTask.MAXIMA_IDENITIFICATION)
		{
			final float sd0 = (float) fitConfig.getInitialPeakStdDev0();
			final float sd1 = (float) fitConfig.getInitialPeakStdDev1();
			for (int n = 0; n < candidates.length; n++)
			{
				// Find the background using the perimeter of the data.
				// TODO - Perhaps the Gaussian Fitter should be used to produce the initial estimates but no actual fit done.
				// This would produce coords using the centre-of-mass.
				final int x = candidates[n].x;
				final int y = candidates[n].y;
				final Rectangle regionBounds = ie.getBoxRegionBounds(x, y, fitting);
				region = ie.crop(regionBounds, region);
				final float b = (float) Gaussian2DFitter.getBackground(region, regionBounds.width, regionBounds.height,
						1);

				// Offset the coords to the centre of the pixel. Note the bounds will be added later.
				// Subtract the background to get the amplitude estimate then convert to signal.
				final float amplitude = candidates[n].intensity - ((relativeIntensity) ? 0 : b);
				final float signal = (float) (amplitude * 2.0 * Math.PI * sd0 * sd1);
				final float[] peakParams = new float[] { b, signal, 0, x + 0.5f, y + 0.5f, sd0, sd1 };
				final int index = y * width + x;
				sliceResults
						.add(createResult(dataBounds.x + x, dataBounds.y + y, data[index], 0, noise, peakParams, null));
			}
		}
		else if (params != null && params.fitTask == FitTask.BENCHMARKING)
		{
			initialiseFitting();

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
				// Note the current fit configuration for benchmarking may have minimal filtering settings
				// so we do not use that object.
				final FitConfiguration tmp = new FitConfiguration();
				final double residualsThreshold = 0.4;
				filter = new MultiPathFilter(tmp, createMinimalFilter(), residualsThreshold);
			}

			filter.setup();

			dynamicMultiPathFitResult = new DynamicMultiPathFitResult(ie, false);

			// Fit all candidates
			for (int n = 0; n < candidates.length; n++)
			{
				// Fit all paths. 
				// The benchmarkFit method uses the filter to store the best results for the slice.
				job.setMultiPathFitResult(n, benchmarkFit(n, filter));
			}
		}
		else
		{
			initialiseFitting();

			// Perform the Gaussian fit

			// Allow the results to be filtered for certain peaks
			if (params != null && params.filter != null)
			{
				resultFilter = new DistanceResultFilter(params.filter, params.distanceThreshold, candidates.length);
				//filter = new OptimumDistanceResultFilter(params.filter, params.distanceThreshold, maxIndices.length);
			}

			// Use the SpotFitter is used to create a dynamic MultiPathFitResult object.
			// This is then passed to a multi-path filter. Thus the same fitting decision process 
			// is used when benchmarking and when running on actual data.

			// Note: The SpotFitter labels each PreprocessedFitResult using the offset in the FitResult object.
			// The initial params and deviations can then be extracted for the results that pass the filter.

			MultiPathFilter filter = new MultiPathFilter(fitConfig, createMinimalFilter(),
					config.getResidualsThreshold());
			IMultiPathFitResults multiPathResults = this;
			SelectedResultStore store = this;
			dynamicMultiPathFitResult = new DynamicMultiPathFitResult(ie, true);

			filter.select(multiPathResults, config.getFailuresLimit(), true, store);

			if (logger != null)
				logger.info("Slice %d: %d / %d", slice, success, candidates.length);

			// Result filter post-processing
			if (resultFilter != null)
			{
				resultFilter.finalise();
				job.setResults(resultFilter.getResults());
				job.setIndices(resultFilter.getMaxIndices());

				for (int i = 0; i < resultFilter.getFilteredCount(); i++)
				{
					job.setFitResult(i, resultFilter.getFitResults()[i]);
				}

				sliceResults.clear();
				sliceResults.addAll(resultFilter.getResults());
			}
		}

		// Add the ROI bounds to the fitted peaks
		float offsetx = dataBounds.x;
		float offsety = dataBounds.y;

		if (params != null && params.getOffset() != null)
		{
			offsetx += params.getOffset()[0];
			offsety += params.getOffset()[1];
		}

		for (int i = 0; i < sliceResults.size(); i++)
		{
			final PeakResult result = sliceResults.get(i);
			result.params[Gaussian2DFunction.X_POSITION] += offsetx;
			result.params[Gaussian2DFunction.Y_POSITION] += offsety;
		}

		this.results.addAll(sliceResults);

		finish(job, start);
	}

	private Candidate[] indentifySpots(FitJob job, int width, int height, FitParameters params)
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
			logger.info("%d: Slice %d: %d candidates", workerId, slice, spots.length);

		sliceResults = new ArrayList<PeakResult>(spots.length);
		if (requireIndices)
		{
			job.setResults(sliceResults);
			job.setIndices(maxIndices);
		}

		final Candidate[] candidates = new Candidate[spots.length];
		for (int i = 0; i < candidates.length; i++)
			candidates[i] = new Candidate(spots[i], i);
		return candidates;
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

	private void initialiseFitting()
	{
		neighbourIndices = allocateArray(neighbourIndices, candidates.length);
		// Allocate enough room for all fits to be doublets 
		fittedNeighbourIndices = allocateArray(fittedNeighbourIndices,
				candidates.length * ((config.getResidualsThreshold() < 1) ? 2 : 1));

		clearEstimates(candidates.length);

		gridManager = new ResultGridManager(dataBounds.width, dataBounds.height, 2 * fitting + 1);
		for (int i = 0; i < candidates.length; i++)
		{
			gridManager.putOnGrid(candidates[i]);
		}
	}

	// Used to add results to the grid for the current fit position.
	// This prevents filtering duplicates within the current fit results, 
	// only with all that has been fit before.
	private PeakResult[] queue = new PeakResult[5];
	private int queueSize = 0;

	private void queueToGrid(PeakResult result)
	{
		if (queueSize == queue.length)
			queue = Arrays.copyOf(queue, queueSize * 2);
		queue[queueSize++] = result;
	}

	private void flushToGrid()
	{
		for (int i = 0; i < queueSize; i++)
			gridManager.putOnGrid(queue[i]);
		queueSize = 0;
		clearGridCache();
	}

	private void clearGridCache()
	{
		gridManager.clearCache();
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

	//	private int addResults(int n, Rectangle regionBounds, double[] peakParams, double[] peakParamsDev, float value,
	//			double error, float noise)
	//	{
	//		final int npeaks = peakParams.length / 6;
	//		int count = 0;
	//
	//		// Note that during processing the data is assumed to refer to the top-left
	//		// corner of the pixel. The coordinates should be represented in the middle of the pixel 
	//		// so add a 0.5 shift to the coordinates.
	//
	//		// WARNING: We do not update the initialParameters in the FitResult with the same offset.
	//
	//		for (int i = 0; i < npeaks; i++)
	//		{
	//			peakParams[i * 6 + Gaussian2DFunction.X_POSITION] += 0.5 + regionBounds.x;
	//			peakParams[i * 6 + Gaussian2DFunction.Y_POSITION] += 0.5 + regionBounds.y;
	//		}
	//
	//		if (npeaks == 1)
	//		{
	//			if (addSingleResult(neighbours, x, y, Utils.toFloat(peakParams), Utils.toFloat(peakParamsDev), value, error,
	//					noise))
	//				count++;
	//		}
	//		else
	//		{
	//			for (int i = 0; i < npeaks; i++)
	//			{
	//				final float[] params = new float[] { (float) peakParams[0], (float) peakParams[i * 6 + 1],
	//						(float) peakParams[i * 6 + 2], (float) peakParams[i * 6 + 3], (float) peakParams[i * 6 + 4],
	//						(float) peakParams[i * 6 + 5], (float) peakParams[i * 6 + 6] };
	//				float[] paramsStdDev = null;
	//				if (peakParamsDev != null)
	//				{
	//					paramsStdDev = new float[] { (float) peakParamsDev[0], (float) peakParamsDev[i * 6 + 1],
	//							(float) peakParamsDev[i * 6 + 2], (float) peakParamsDev[i * 6 + 3],
	//							(float) peakParamsDev[i * 6 + 4], (float) peakParamsDev[i * 6 + 5],
	//							(float) peakParamsDev[i * 6 + 6] };
	//				}
	//				if (addSingleResult(neighbours, x, y, params, paramsStdDev, value, error, noise))
	//					count++;
	//			}
	//		}
	//
	//		return count;
	//	}

	/**
	 * Add the result to the list. Only check for duplicates in the current results grid.
	 *
	 * @param candidateId
	 *            the candidate id
	 * @param peakParams
	 *            the peak params
	 * @param peakParamsDev
	 *            the peak params dev
	 * @param error
	 *            the error
	 * @param noise
	 *            the noise
	 * @return true, if successful
	 */
	private boolean addSingleResult(int candidateId, float[] peakParams, float[] peakParamsDev, double error)
	{
		int x = candidates[candidateId].x;
		int y = candidates[candidateId].y;
		float value = data[y * dataBounds.width + x];

		if (duplicateDistance2 > 0)
		{
			// Check for duplicates
			final PeakResult[] neighbours = gridManager.getPeakResultNeighbours(x, y);
			for (int i = 0; i < neighbours.length; i++)
			{
				if (distance2(neighbours[i].params, peakParams) < duplicateDistance2)
				{
					if (logger != null)
						logger.info("[%d] Ignoring duplicate peak @ %.2f,%.2f", slice,
								peakParams[Gaussian2DFunction.X_POSITION], peakParams[Gaussian2DFunction.Y_POSITION]);
					return false;
				}
			}
		}

		// Update to the global bounds.
		// (Note the global bounds will be added to the params at the end of processing the frame
		// so we leave those untouched)
		x += dataBounds.x;
		y += dataBounds.y;

		// This was fit OK so add it to the grid of results (so we do not fit it again)
		final PeakResult peakResult = createResult(x, y, value, error, noise, peakParams, peakParamsDev);
		queueToGrid(peakResult);
		candidates[candidateId].fit = true;

		// Check if the position is inside the border tolerance
		if (insideBorder(peakParams[Gaussian2DFunction.X_POSITION], peakParams[Gaussian2DFunction.Y_POSITION]))
		{
			sliceResults.add(peakResult);
			fittedBackground += peakParams[Gaussian2DFunction.BACKGROUND];
		}
		else if (logger != null)
		{
			logger.info("[%d] Ignoring peak within image border @ %.2f,%.2f", slice,
					peakParams[Gaussian2DFunction.X_POSITION], peakParams[Gaussian2DFunction.Y_POSITION]);
		}
		return true;
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

	/**
	 * Provide the ability to convert raw fitted results into PreprocessedPeakResult for validation
	 */
	private abstract class ResultFactory
	{
		final float offsetx, offsety;

		ResultFactory(float offsetx, float offsety)
		{
			this.offsetx = offsetx;
			this.offsety = offsety;
		}

		abstract PreprocessedPeakResult createPreprocessedPeakResult(int candidateId, int n, double[] initialParams,
				double[] params, double localBackground, ResultType resultType);
	}

	/**
	 * Provide dynamic PreprocessedPeakResult. This is basically a wrapper around the result arrays that provides
	 * properties on-the-fly.
	 */
	private class DynamicResultFactory extends ResultFactory
	{
		DynamicResultFactory(float offsetx, float offsety)
		{
			super(offsetx, offsety);
		}

		PreprocessedPeakResult createPreprocessedPeakResult(int candidateId, int n, double[] initialParams,
				double[] params, double localBackground, ResultType resultType)
		{
			return fitConfig.createDynamicPreprocessedPeakResult(candidateId, n, initialParams, params, localBackground,
					resultType, offsetx, offsety);
		}
	}

	/**
	 * Provide a materialised PreprocessedPeakResult as a new object with all properties computed.
	 */
	private class FixedResultFactory extends ResultFactory
	{
		FixedResultFactory(float offsetx, float offsety)
		{
			super(offsetx, offsety);
		}

		PreprocessedPeakResult createPreprocessedPeakResult(int candidateId, int n, double[] initialParams,
				double[] params, double localBackground, ResultType resultType)
		{
			return fitConfig.createPreprocessedPeakResult(slice, candidateId, n, initialParams, params, localBackground,
					resultType, offsetx, offsety);
		}
	}

	/**
	 * Provide functionality to fit spots in a region using different methods. Decisions about what to accept are not
	 * performed. The fit results are just converted to PreprocessedPeakResult objects for validation.
	 */
	private class SpotFitter
	{
		final Gaussian2DFitter gf;
		final ResultFactory resultFactory;
		final double[] region;
		final Rectangle regionBounds;
		final int candidateId;
		final int width;
		final int height;
		final int parametersPerPeak = 6;

		int neighbours;
		float singleBackground = Float.NaN, multiBackground = Float.NaN;

		MultiPathFitResult.FitResult resultMulti = null;
		MultiPathFitResult.FitResult resultSingle = null;
		double[] singleRegion = null;
		MultiPathFitResult.FitResult resultDoublet = null;
		boolean computedDoublet = false;
		QuadrantAnalysis singleQA = null;

		public SpotFitter(Gaussian2DFitter gf, ResultFactory resultFactory, double[] region, Rectangle regionBounds,
				int n)
		{
			this.gf = gf;
			this.resultFactory = resultFactory;
			this.region = region;
			this.regionBounds = regionBounds;
			this.candidateId = n;

			// Initialise
			width = regionBounds.width;
			height = regionBounds.height;

			// Analyse neighbours and include them in the fit if they are within a set height of this peak.
			resetNeighbours();
			neighbours = (config.isIncludeNeighbours()) ? findNeighbours(regionBounds, n) : 0;
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

		public MultiPathFitResult.FitResult getResultMulti()
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

			// Support bounds on the known fitted peaks
			double[] lower = new double[params.length];
			double[] upper = new double[params.length];
			for (int i = 0; i < lower.length; i++)
			{
				lower[i] = Double.NEGATIVE_INFINITY;
				upper[i] = Double.POSITIVE_INFINITY;
			}

			// Note: If difference-of-smoothing is performed the heights have background subtracted so 
			// it must be added back 

			// The main peak
			params[Gaussian2DFunction.SIGNAL] = candidates[candidateId].intensity +
					((relativeIntensity) ? background : 0);
			params[Gaussian2DFunction.X_POSITION] = candidates[candidateId].x - regionBounds.x;
			params[Gaussian2DFunction.Y_POSITION] = candidates[candidateId].y - regionBounds.y;

			// The neighbours
			for (int i = 0, j = parametersPerPeak; i < neighbourCount; i++, j += parametersPerPeak)
			{
				final int n2 = neighbourIndices[i];
				final double candidateX = candidates[n2].x - regionBounds.x;
				final double candidateY = candidates[n2].y - regionBounds.y;
				
				final double[] estimatedParams = getEstimate(n2);
				if (estimatedParams != null)
				{
					// Re-use previous good multi-fit results to estimate the peak params...
					// Convert signal into amplitude
					params[j + Gaussian2DFunction.SIGNAL] = estimatedParams[Gaussian2DFunction.SIGNAL] / (2 * Math.PI *
							estimatedParams[Gaussian2DFunction.X_SD] * estimatedParams[Gaussian2DFunction.Y_SD]);
					params[j + Gaussian2DFunction.X_POSITION] = estimatedParams[Gaussian2DFunction.X_POSITION] -
							regionBounds.x;
					params[j + Gaussian2DFunction.Y_POSITION] = estimatedParams[Gaussian2DFunction.Y_POSITION] -
							regionBounds.y;
					params[j + Gaussian2DFunction.ANGLE] = estimatedParams[Gaussian2DFunction.ANGLE];
					params[j + Gaussian2DFunction.X_SD] = estimatedParams[Gaussian2DFunction.X_SD];
					params[j + Gaussian2DFunction.Y_SD] = estimatedParams[Gaussian2DFunction.Y_SD];
				}
				else
				{
					params[j + Gaussian2DFunction.SIGNAL] = candidates[n2].intensity +
							((relativeIntensity) ? background : 0);
					params[j + Gaussian2DFunction.X_POSITION] = candidateX;
					params[j + Gaussian2DFunction.Y_POSITION] = candidateY;
				}
				
				// Constrain the location using the candidate position.
				// Do not use the current estimate as this will create drift over time if the estimate is updated.
				lower[j + Gaussian2DFunction.X_POSITION] = candidateX - 1;
				upper[j + Gaussian2DFunction.X_POSITION] = candidateX + 1;
				lower[j + Gaussian2DFunction.Y_POSITION] = candidateY - 1;
				upper[j + Gaussian2DFunction.Y_POSITION] = candidateY + 1;
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

					final GaussianFunction func = fitConfig.createGaussianFunction(subtractFittedPeaks,
							regionBounds.width, funcParams);
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
					final PeakResult result = sliceResults.get(fittedNeighbourIndices[i]);

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
					
					// Add support for constraining the known fit results using bounded coordinates.
					// Currently we just constrain the location.
					lower[j + Gaussian2DFunction.X_POSITION] = params[j + Gaussian2DFunction.X_POSITION] - 1;
					upper[j + Gaussian2DFunction.X_POSITION] = params[j + Gaussian2DFunction.X_POSITION] + 1;
					lower[j + Gaussian2DFunction.Y_POSITION] = params[j + Gaussian2DFunction.Y_POSITION] - 1;
					upper[j + Gaussian2DFunction.Y_POSITION] = params[j + Gaussian2DFunction.Y_POSITION] + 1;					
					
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

			gf.setBounds(lower, upper);
			final FitResult fitResult = gf.fit(region, width, height, npeaks, params, true);
			gf.setBounds(null, null);			

			// Restore
			fitConfig.setFitValidation(fitValidation);
			fitConfig.setMaxIterations(maxIterations);
			fitConfig.setMaxFunctionEvaluations(maxEvaluations);

			updateError(fitResult);

			// Create the results
			PreprocessedPeakResult[] results = null;
			if (fitResult.getStatus() == FitStatus.OK)
			{
				final double[] fitParams = fitResult.getParameters();
				final double[] initialParams = fitResult.getInitialParameters();
				convertParameters(fitParams);

				results = new PreprocessedPeakResult[npeaks];

				// We must compute a local background for all the spots
				final int flags = fitConfig.getFunctionFlags();

				results[0] = resultFactory.createPreprocessedPeakResult(this.candidateId, 0, initialParams, fitParams,
						getLocalBackground(0, npeaks, fitParams, flags), ResultType.NEW);

				// Already fitted peaks
				for (int n = neighbourCount + 1; n < npeaks; n++)
				{
					results[n] = resultFactory.createPreprocessedPeakResult(this.candidateId, n, initialParams,
							fitParams, getLocalBackground(n, npeaks, fitParams, flags), ResultType.EXISTING);
				}

				// Neighbours
				for (int n = 1; n <= neighbourCount; n++)
				{
					results[n] = resultFactory.createPreprocessedPeakResult(this.candidateId, n, initialParams,
							fitParams, getLocalBackground(n, npeaks, fitParams, flags), ResultType.CANDIDATE);
				}
			}

			return resultMulti = createResult(fitResult, results);
		}

		private double getLocalBackground(int n, int npeaks, double[] params, final int flags)
		{
			GaussianOverlapAnalysis overlap = new GaussianOverlapAnalysis(flags, extractSpotParams(params, n), 2);
			overlap.add(flags, extractOtherParams(params, n, npeaks), true);
			double[] overlapData = overlap.getOverlapData();
			return overlapData[1] + params[0];
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
		 *            the result
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

		public MultiPathFitResult.FitResult getResultSingle()
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
				signal = candidates[candidateId].intensity - ((relativeIntensity) ? 0 : background);
				if (signal < 0)
				{
					signal += background;
					background = 0;
				}
			}
			final double[] params = new double[] { background, signal, 0, candidates[candidateId].x - regionBounds.x,
					candidates[candidateId].y - regionBounds.y, 0, 0 };

			// If there were unfitted neighbours then use off grid pixels to prevent re-estimate of CoM
			// since the CoM estimate will be skewed by the neighbours
			if (neighbourCount != 0)
			{
				params[Gaussian2DFunction.X_POSITION] += 0.001;
				params[Gaussian2DFunction.Y_POSITION] += 0.001;
			}

			final FitResult fitResult = gf.fit(region, width, height, 1, params, amplitudeEstimate);

			updateError(fitResult);

			// Create the results
			PreprocessedPeakResult[] results = null;
			if (fitResult.getStatus() == FitStatus.OK)
			{
				results = new PreprocessedPeakResult[1];
				final double[] fitParams = fitResult.getParameters();
				final double[] initialParams = fitResult.getInitialParameters();
				convertParameters(params);
				results[0] = resultFactory.createPreprocessedPeakResult(this.candidateId, 0, initialParams, fitParams,
						0, ResultType.NEW);
			}

			resultSingle = createResult(fitResult, results);

			return resultSingle;
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

		private MultiPathFitResult.FitResult createResult(FitResult fitResult, PreprocessedPeakResult[] results)
		{
			MultiPathFitResult.FitResult mfitResult = new MultiPathFitResult.FitResult(fitResult.getStatus().ordinal(),
					fitResult);
			mfitResult.results = results;
			return mfitResult;
		}

		public double getSingleQAScore()
		{
			if (singleQA != null)
				return singleQA.score;

			// Ensure we have a single result
			getResultSingle();

			// Note this assumes that this method will be called after a single fit and that the 
			// residuals were computed.
			final double[] residuals = gf.getResiduals();
			singleQA = computeQA((FitResult) resultSingle.data, regionBounds, residuals);

			return singleQA.score;
		}

		public MultiPathFitResult.FitResult getResultSingleDoublet(double residualsThreshold)
		{
			if (computedDoublet)
				return resultDoublet;

			if (residualsThreshold >= 1 || residualsThreshold < 0)
				return null;

			if (getSingleQAScore() < residualsThreshold)
				return null;

			computedDoublet = true;

			final CandidateList neighbours = findNeighbours(candidates[candidateId]);
			final PeakResult[] peakNeighbours = findPeakNeighbours(candidates[candidateId]);

			// Use the region from the single fit which had fitted peaks subtracted
			final double[] region = singleRegion;

			resultDoublet = fitAsDoublet((FitResult) resultSingle.data, region, regionBounds, residualsThreshold,
					neighbours, peakNeighbours, singleQA);

			return resultDoublet;
		}

		/**
		 * Perform quadrant analysis on the residuals
		 *
		 * @param fitResult
		 *            the fit result
		 * @param regionBounds
		 *            the region bounds
		 * @param residuals
		 *            the residuals
		 * @return the multi path fit result. fit result
		 */
		private QuadrantAnalysis computeQA(FitResult fitResult, Rectangle regionBounds, final double[] residuals)
		{
			if (residuals == null)
			{
				QuadrantAnalysis qa = new QuadrantAnalysis();
				qa.score = -1; // Set so that any residuals threshold will be ignored.
				return qa;
			}

			final double[] params = fitResult.getParameters();
			final int width = regionBounds.width;
			final int height = regionBounds.height;

			// Use rounding since the fit coords are not yet offset by 0.5 pixel to centre them
			final int cx = (int) Math.round(params[Gaussian2DFunction.X_POSITION]);
			final int cy = (int) Math.round(params[Gaussian2DFunction.Y_POSITION]);

			QuadrantAnalysis qa = new QuadrantAnalysis();
			qa.quadrantAnalysis(residuals, width, height, cx, cy);

			if (logger != null)
				logger.info("Residue analysis = %f (%d,%d)", qa.score, qa.vector[0], qa.vector[1]);

			return qa;
		}

		/**
		 * Fit the single spot location as two spots (a doublet).
		 * <p>
		 * Perform quadrant analysis as per rapidSTORM. If the residuals of the the fit are skewed around the single fit
		 * result then an attempt is made to fit two spots (a doublet)
		 *
		 * @param fitResult
		 *            the fit result
		 * @param region
		 *            the region
		 * @param regionBounds
		 *            the region bounds
		 * @param residuals
		 *            the residuals
		 * @param residualsThreshold
		 *            the residuals threshold
		 * @param neighbours
		 *            the neighbours
		 * @param peakNeighbours
		 *            the peak neighbours
		 * @return the multi path fit result. fit result
		 */
		private MultiPathFitResult.FitResult fitAsDoublet(FitResult fitResult, double[] region, Rectangle regionBounds,
				double residualsThreshold, CandidateList neighbours, PeakResult[] peakNeighbours, QuadrantAnalysis qa)
		{
			final int width = regionBounds.width;
			final int height = regionBounds.height;
			final double[] params = fitResult.getParameters();
			// Use rounding since the fit coords are not yet offset by 0.5 pixel to centre them
			final int cx = (int) Math.round(params[Gaussian2DFunction.X_POSITION]);
			final int cy = (int) Math.round(params[Gaussian2DFunction.Y_POSITION]);

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

			// If validation is on in the fit configuration,
			// disable checking of position movement since we guessed the 2-peak location.
			// Also increase the iterations level then reset afterwards.

			// TODO - Should width and signal validation be disabled too?
			double shift = fitConfig.getCoordinateShift();
			final int maxIterations = fitConfig.getMaxIterations();
			final int maxEvaluations = fitConfig.getMaxFunctionEvaluations();

			fitConfig.setCoordinateShift(FastMath.min(width, height));
			fitConfig.setMaxIterations(maxIterations * ITERATION_INCREASE_FOR_DOUBLETS);
			fitConfig.setMaxFunctionEvaluations(maxEvaluations * FitWorker.EVALUATION_INCREASE_FOR_DOUBLETS);

			// We assume that residuals calculation is on but just in case something else turned it off we get the state.
			final boolean isComputeResiduals = gf.isComputeResiduals();
			gf.setComputeResiduals(false);
			final FitResult newFitResult = gf.fit(region, width, height, 2, doubletParams, false);
			gf.setComputeResiduals(isComputeResiduals);

			fitConfig.setCoordinateShift(shift);
			fitConfig.setMaxIterations(maxIterations);
			fitConfig.setMaxFunctionEvaluations(maxEvaluations);

			if (newFitResult.getStatus() == FitStatus.OK)
			{
				// Adjusted Coefficient of determination is not good for non-linear models. Use the 
				// Bayesian Information Criterion (BIC):

				final double doubleSumOfSquares = gf.getFinalResidualSumOfSquares();

				final int length = width * height;
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
						logger.info("Model improvement - Sum-of-squares (IC) : %f (%f) => %f (%f) : %f",
								singleSumOfSquares, ic1, doubleSumOfSquares, ic2, ic1 - ic2);
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
					String msg = String.format(
							"Doublet %d [%d,%d] %s (%s) [%f -> %f] SS [%f -> %f] IC [%f -> %f] = %s\n", slice,
							cx + dataBounds.x + regionBounds.x, cy + dataBounds.y + regionBounds.y,
							newFitResult.getStatus(), newFitResult.getStatusData(), singleValue, gf.getValue(),
							singleSumOfSquares, doubleSumOfSquares, ic1, ic2, Arrays.toString(peakParams));
					logger2.debug(msg);
				}

				// Check if the predictive power of the model is better with two peaks:
				// IC should be lower
				if (ic2 > ic1)
				{
					return null;
				}

				// Validation of fit. For each spot:
				// 1. Check the spot is inside the region
				// 2. Check the distance of each new centre from the original centre.
				//    If the shift is too far (e.g. half the distance to the edge), the centre 
				//    must be in the correct quadrant. 
				// 3. Check we are not closer to a fitted spot. This has already had a chance at 
				//    fitting a doublet so is ignored..
				// 4. Check if there is an unfit candidate spot closer than the current candidate.
				//    This represents drift out to fit another spot that will be fit later.

				final double[] newParams = newFitResult.getParameters();

				// Allow the shift to span half of the fitted window.
				shift = 0.5 * FastMath.min(regionBounds.width, regionBounds.height);

				final int[] position = new int[2];
				final int[] candidateIndex = new int[2];
				int nPeaks = 0;
				NEXT_PEAK: for (int n = 0; n < 2; n++)
				{
					// 1. Check the spot is inside the region

					// Note that during processing the data is assumed to refer to the top-left
					// corner of the pixel. The coordinates should be represented in the middle of the pixel 
					// so add a 0.5 shift to the coordinates.
					final double xpos = newParams[Gaussian2DFunction.X_POSITION + n * 6] + 0.5;
					final double ypos = newParams[Gaussian2DFunction.Y_POSITION + n * 6] + 0.5;

					if (xpos < 0 || xpos > width || ypos < 0 || ypos > height)
					{
						if (logger != null)
							logger.info("Fitted coordinates too far outside the fitted region (x %g || y %g) in %dx%d",
									xpos, ypos, width, height);
						continue;
					}

					// 2. Check the distance of each new centre from the original centre.
					//    If the shift is too far (e.g. half the distance to the edge), the centre 
					//    must be in the correct quadrant.

					final double xShift = newParams[Gaussian2DFunction.X_POSITION + n * 6] -
							params[Gaussian2DFunction.X_POSITION];
					final double yShift = newParams[Gaussian2DFunction.Y_POSITION + n * 6] -
							params[Gaussian2DFunction.Y_POSITION];
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
								logger.info(
										"Bad peak %d: Fitted coordinates moved into wrong quadrant (x=%g,y=%g,a=%f)", n,
										xShift, yShift, a * 57.29578);
							}
							continue;
						}
					}

					// Spots from the doublet must be closer to the single fit than any other spots.
					// This is true for candidates we have yet to fit or already fitted candidates.

					// Distance to current candidate fitted as a single
					final double distanceToSingleFit2 = xShift * xShift + yShift * yShift;

					// 3. Check we are not closer to a fitted spot. This has already had a chance at 
					//    fitting a doublet so is ignored..

					if (peakNeighbours.length != 0)
					{
						// Coords for comparison to the real positions
						final float fcx2 = (float) (regionBounds.x + newParams[Gaussian2DFunction.X_POSITION + n * 6] +
								0.5);
						final float fcy2 = (float) (regionBounds.y + newParams[Gaussian2DFunction.Y_POSITION + n * 6] +
								0.5);
						final float d2 = (float) distanceToSingleFit2;
						for (int i = 0; i < peakNeighbours.length; i++)
						{
							if (d2 > distance2(fcx2, fcy2, peakNeighbours[i].params))
							{
								if (logger != null)
								{
									logger.info(
											"Bad peak %d: Fitted coordinates moved closer to another result (%d,%d : x=%.1f,y=%.1f : %.1f,%.1f)",
											n, candidates[candidateId].x, candidates[candidateId].y, fcx2, fcy2,
											peakNeighbours[i].getXPosition(), peakNeighbours[i].getYPosition());
								}
								// There is another fitted result that is closer.
								// Note: The fit region is not centred on the other spot so this fit will probably
								// be worse and is discarded (not compared to the existing fit to get the best one). 
								continue NEXT_PEAK;
							}
						}
					}

					// 4. Check if there is an unfit candidate spot closer than the current candidate.
					//    This represents drift out to fit another unfitted spot.

					int otherId = candidateId;
					if (neighbours.size != 0)
					{
						// Position - do not add 0.5 pixel offset to allow distance comparison to integer candidate positions.
						float fcx2 = (float) (regionBounds.x + newParams[Gaussian2DFunction.X_POSITION + n * 6]);
						float fcy2 = (float) (regionBounds.y + newParams[Gaussian2DFunction.Y_POSITION + n * 6]);
						double mind2 = distanceToSingleFit2;
						for (int j = 0; j < neighbours.size; j++)
						{
							final int id = neighbours.list[j].index;
							if (isFit(id))
								// This will be in the already fitted results instead so ignore...
								continue;
							final double d2 = distance2(fcx2, fcy2, candidates[id]);
							if (mind2 > d2)
							{
								mind2 = d2;
								otherId = id;
							}
						}
						if (otherId != candidateId)
						{
							if (logger != null)
							{
								logger.info(
										"Bad peak %d: Fitted coordinates moved closer to another candidate (%d,%d : x=%.1f,y=%.1f : %d,%d)",
										n, candidates[candidateId].x, candidates[candidateId].y, fcx2 + 0.5f,
										fcy2 + 0.5f, candidates[otherId].x, candidates[otherId].y);
							}

							// There is another candidate to be fit later that is closer.
							// This may be used as an estimate so we return it as such (i.e we do not ignore it)
						}
					}

					candidateIndex[nPeaks] = otherId;
					position[nPeaks++] = n;
				}

				if (nPeaks == 0)
				{
					return null;
				}

				// Return results for validation
				final double[] fitParams = fitResult.getParameters();
				final double[] initialParams = fitResult.getInitialParameters();
				convertParameters(params);

				final PreprocessedPeakResult[] results = new PreprocessedPeakResult[nPeaks];
				for (int i = 0; i < nPeaks; i++)
				{
					// If it is this candidate, or an earlier one that was not fit then this is a new result.
					// Otherwise it is a candidate we will process later
					final ResultType resultType = (candidateIndex[i] <= candidateId) ? ResultType.NEW
							: ResultType.CANDIDATE;
					results[i] = resultFactory.createPreprocessedPeakResult(candidateIndex[i], position[i],
							initialParams, fitParams, 0, resultType);
				}

				return createResult(fitResult, results);
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
					String msg = String.format("Doublet %d [%d,%d] %s (%s) = %s\n", slice,
							cx + dataBounds.x + regionBounds.x, cy + dataBounds.y + regionBounds.y,
							newFitResult.getStatus(), newFitResult.getStatusData(), Arrays.toString(peakParams));
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
	}

	private double estimateOffsetx, estimateOffsety;

	private void storeEstimate(int i, PreprocessedPeakResult peak, byte filterRank)
	{
		double[] params = peak.toGaussian2DParameters();
		double precision = (fitConfig.isPrecisionUsingBackground()) ? peak.getLocationVariance2()
				: peak.getLocationVariance();
		storeEstimate(i, params, precision, filterRank);
	}

	private void storeEstimate(int i, double[] params, double precision, byte filterRank)
	{
		// Add region offset
		params[Gaussian2DFunction.X_POSITION] += estimateOffsetx;
		params[Gaussian2DFunction.Y_POSITION] += estimateOffsety;

		// Compute distance to spot
		final double dx = candidates[i].x - params[Gaussian2DFunction.X_POSITION];
		final double dy = candidates[i].y - params[Gaussian2DFunction.Y_POSITION];
		final double d2 = dx * dx + dy * dy;

		if (estimates[i] == null || estimates[i].isWeaker(filterRank, d2, precision))
			estimates[i] = new Estimate(params, filterRank, d2, precision);
	}

	/**
	 * Extract parameters for the specified peak. The background is ignored.
	 *
	 * @param params
	 *            the params
	 * @param n
	 *            the peak
	 * @return the extracted params
	 */
	private static double[] extractSpotParams(double[] params, int n)
	{
		final double[] newParams = new double[7];
		System.arraycopy(params, n * 6 + 1, newParams, 1, 6);
		return newParams;
	}

	/**
	 * Extract parameters other than the specified peak. The background is ignored.
	 *
	 * @param params
	 *            the params
	 * @param n
	 *            the peak
	 * @param nPeaks
	 *            the n peaks
	 * @return the extracted params
	 */
	private static double[] extractOtherParams(double[] params, int n, int nPeaks)
	{
		final double[] newParams = new double[params.length - 6];
		if (n > 0)
		{
			System.arraycopy(params, 1, newParams, 1, n * 6);
		}
		final int left = nPeaks - (n + 1);
		if (left > 0)
		{
			System.arraycopy(params, (n + 1) * 6 + 1, newParams, n * 6 + 1, left * 6);
		}
		return newParams;
	}

	/**
	 * Fits a 2D Gaussian to the given data. Fits all possible paths and stores the multi-path result.
	 * <p>
	 * The best result is chosen using the provided MultiPathFilter and stored in the sliceResults array.
	 *
	 * @param n
	 *            The candidate to fit
	 * @param filter
	 *            The filter for choosing the best result to add to the slice results
	 * @return The full multi-path result
	 */
	private MultiPathFitResult benchmarkFit(int n, MultiPathFilter filter)
	{
		clearGridCache();

		dynamicMultiPathFitResult.reset(n);

		MultiPathFitResult.FitResult result;

		result = dynamicMultiPathFitResult.getMultiFitResult();
		if (result != null && result.getStatus() == 0)
		{
			// Note that if we have neighbours it is very possible we will have doublets
			// due to high density data. It makes sense to attempt a doublet fit here.

			// TODO - Doublet analysis here if all neighbour peaks were valid using the input filter.
		}

		// Do a single fit
		result = dynamicMultiPathFitResult.getSingleFitResult();
		if (result != null && result.getStatus() == 0)
		{
			dynamicMultiPathFitResult.getSingleQAScore();
			dynamicMultiPathFitResult.getDoubletFitResult();
		}

		// Pick the best result.
		// Store the new results and use the passed candidates as estimates.
		// This is done by passing this class in as the SelectedResultStore.
		// Results passing the primary filter as added to the current slice results.
		// Results passing the minimal filter are used as estimates.

		final MultiPathFitResult multiPathFitResult = dynamicMultiPathFitResult.copy();
		final SelectedResult selectedResult = filter.select(multiPathFitResult, true, this);

		if (selectedResult != null)
		{
			add(selectedResult);
		}

		flushToGrid();

		return multiPathFitResult;
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

	private void resetNeighbours()
	{
		neighbourCount = 0;
		fittedNeighbourCount = 0;
		neighbours = null;
		peakNeighbours = null;
	}

	private CandidateList findNeighbours(Candidate candidate)
	{
		if (neighbours == null)
		{
			// Using the neighbour grid
			neighbours = gridManager.getNeighbours(candidate);
			neighbours.sort();
		}
		return neighbours;
	}

	private PeakResult[] findPeakNeighbours(Candidate candidate)
	{
		if (peakNeighbours == null)
		{
			// Using the neighbour grid 
			peakNeighbours = gridManager.getPeakResultNeighbours(candidate.x, candidate.y);
		}
		return peakNeighbours;
	}

	/**
	 * Search for any neighbours within a set height of the specified peak that is within the search region bounds.
	 *
	 * @param regionBounds
	 *            the region bounds
	 * @param n
	 *            the candidate index
	 * @return The number of neighbours
	 */
	private int findNeighbours(Rectangle regionBounds, int n)
	{
		int xmin = regionBounds.x;
		int xmax = xmin + regionBounds.width - 1;
		int ymin = regionBounds.y;
		int ymax = ymin + regionBounds.height - 1;

		final Candidate spot = candidates[n];

		final float heightThreshold;
		if (relativeIntensity)
		{
			// No background when spot filter has relative intensity
			heightThreshold = (float) (spot.intensity * config.getNeighbourHeightThreshold());
		}
		else
		{
			if (spot.intensity < background)
				heightThreshold = spot.intensity;
			else
				heightThreshold = (float) ((spot.intensity - background) * config.getNeighbourHeightThreshold() +
						background);
		}

		// Check all maxima that are lower than this
		neighbourCount = 0;

		// Using the neighbour grid.
		// Note this will also include all higher intensity spots that failed to be fit. 
		// These may still have estimates.
		final CandidateList neighbours = findNeighbours(spot);
		final Candidate[] candidates = neighbours.list;
		for (int i = 0; i < neighbours.size; i++)
		{
			final int id = candidates[i].index;
			if (isFit(id))
				continue;
			if (canIgnore(candidates[i].x, candidates[i].y, xmin, xmax, ymin, ymax, candidates[i].intensity,
					heightThreshold))
				continue;
			neighbourIndices[neighbourCount++] = id;
		}

		int c = 0;
		// Processing all lower spots.
		//for (int i = n + 1; i < candidates.length; i++)
		// Processing all spots.
		for (int i = 0; i < candidates.length; i++)
		{
			if (i == n || isFit(i))
				continue;
			if (canIgnore(candidates[i].x, candidates[i].y, xmin, xmax, ymin, ymax, candidates[i].intensity,
					heightThreshold))
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

			final PeakResult[] peakNeighbours = findPeakNeighbours(spot);

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

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IMultiPathFitResults#getFrame()
	 */
	public int getFrame()
	{
		return slice;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IMultiPathFitResults#getNumberOfResults()
	 */
	public int getNumberOfResults()
	{
		return getTotalCandidates();
	}

	/**
	 * Use for result filtering
	 */
	private ResultFilter resultFilter = null;

	/**
	 * Provide the multi-path fit results dynamically
	 */
	private class DynamicMultiPathFitResult extends MultiPathFitResult
	{
		/** The constant for no Quadrant Analysis score */
		private static final double NO_QA_SCORE = 2;

		final ImageExtractor ie;
		boolean dynamic;
		Rectangle regionBounds;
		double[] region;
		SpotFitter spotFitter;
		FitType fitType;

		public DynamicMultiPathFitResult(ImageExtractor ie, boolean dynamic)
		{
			this.frame = FitWorker.this.slice;
			this.width = FitWorker.this.dataBounds.width;
			this.height = FitWorker.this.dataBounds.height;
			this.ie = ie;
			this.dynamic = dynamic;
		}

		public void reset(int candidateId)
		{
			this.candidateId = candidateId;
			fitType = new FitType();

			// Reset results
			this.setSingleQAScore(NO_QA_SCORE);
			this.setMultiFitResult(null);
			this.setSingleFitResult(null);
			this.setDoubletFitResult(null);

			// Set fitting region
			regionBounds = ie.getBoxRegionBounds(candidates[candidateId].x, candidates[candidateId].y, fitting);
			region = ie.crop(regionBounds, region);

			// Offsets to convert fit coordinates to the data frame
			final float offsetx = dataBounds.x + regionBounds.x + 0.5f;
			final float offsety = dataBounds.y + regionBounds.y + 0.5f;

			// Note that the PreprocessedPeakResult will have corrected coordinates so
			// we subtract them and make relative to the fit region bounds when storing 
			// estimates.
			estimateOffsetx = regionBounds.x - offsetx; // == -bounds.x - 0.5f
			estimateOffsety = regionBounds.y - offsety; // == -bounds.y - 0.5f

			final ResultFactory factory = (dynamic) ? new DynamicResultFactory(offsetx, offsety)
					: new FixedResultFactory(offsetx, offsety);
			spotFitter = new SpotFitter(gf, factory, region, regionBounds, candidateId);
		}

		@Override
		public FitResult getMultiFitResult()
		{
			FitResult result = super.getMultiFitResult();
			if (result == null)
			{
				result = spotFitter.getResultMulti();
				setMultiFitResult(result);
				fitType.setNeighbours(true);
			}
			return result;
		}

		public FitResult getSuperMultiFitResult()
		{
			// Pass through the reference to the result
			return super.getMultiFitResult();
		}

		@Override
		public FitResult getSingleFitResult()
		{
			FitResult result = super.getSingleFitResult();
			if (result == null)
			{
				result = spotFitter.getResultSingle();
				setSingleFitResult(result);
			}
			return result;
		}

		@Override
		public double getSingleQAScore()
		{
			double score = super.getSingleQAScore();
			if (score == NO_QA_SCORE)
			{
				score = spotFitter.getSingleQAScore();
				this.setSingleQAScore(score);
			}
			return score;
		}

		@Override
		public FitResult getDoubletFitResult()
		{
			FitResult result = super.getDoubletFitResult();
			if (result == null)
			{
				result = spotFitter.getResultSingleDoublet(config.getResidualsThreshold());
				setDoubletFitResult(result);
				fitType.setDoublet(true);
			}
			return result;
		}

		public FitResult getSuperDoubletFitResult()
		{
			// Pass through the reference to the result
			return super.getDoubletFitResult();
		}
	}

	private DynamicMultiPathFitResult dynamicMultiPathFitResult;

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IMultiPathFitResults#getResult(int)
	 */
	public MultiPathFitResult getResult(int index)
	{
		dynamicMultiPathFitResult.reset(index);
		return dynamicMultiPathFitResult;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IMultiPathFitResults#getTotalCandidates()
	 */
	public int getTotalCandidates()
	{
		return candidates.length;
	}

	/**
	 * Count the number of successful fits
	 */
	private int success = 0;

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.MultiPathFilter.SelectedResultStore#add(gdsc.smlm.results.filter.MultiPathFilter.
	 * SelectedResult)
	 */
	public void add(SelectedResult selectedResult)
	{
		// Add to the slice results.
		final PreprocessedPeakResult[] results = selectedResult.results;
		if (results == null)
			return;

		final int currrentSize = sliceResults.size();
		final int candidateId = dynamicMultiPathFitResult.candidateId;

		final FitResult fitResult = (FitResult) selectedResult.fitResult.data;

		// The background for each result was the local background. We want the fitted global background
		final float background = (float) fitResult.getParameters()[0];
		final double[] dev = fitResult.getParameterStdDev();

		for (int i = 0; i < results.length; i++)
		{
			if (results[i].isExistingResult())
				continue;
			if (results[i].isNewResult())
			{
				final double[] p = ((BasePreprocessedPeakResult) results[i]).toGaussian2DParameters();
				// Convert back so that 0,0 is the top left of the data bounds
				p[Gaussian2DFunction.X_POSITION] -= dataBounds.x;
				p[Gaussian2DFunction.Y_POSITION] -= dataBounds.y;

				final float[] params = new float[7];
				params[Gaussian2DFunction.BACKGROUND] = background;
				for (int j = 1; j < 7; j++)
					params[j] = (float) p[j];

				final float[] paramsDev;
				if (dev == null)
				{
					paramsDev = null;
				}
				else
				{
					paramsDev = new float[7];
					paramsDev[Gaussian2DFunction.BACKGROUND] = (float) dev[Gaussian2DFunction.BACKGROUND];
					final int offset = results[i].getId() * 6;
					for (int j = 1; j < 7; j++)
						paramsDev[j] = (float) dev[offset + j];
				}

				addSingleResult(results[i].getCandidateId(), params, paramsDev, fitResult.getError());
			}
			else
			{
				// This is a candidate that passed validation. Store the estimate as passing the primary filter.
				storeEstimate(results[i].getCandidateId(), results[i], FILTER_RANK_PRIMARY);
			}
		}

		job.setFitResult(candidateId, fitResult);

		// Reporting
		if (this.counter != null)
		{
			FitType fitType = dynamicMultiPathFitResult.fitType;
			if (selectedResult.fitResult.getStatus() == 0)
			{
				fitType.setOK(true);
				if (dynamicMultiPathFitResult.getSuperMultiFitResult() == selectedResult.fitResult)
					fitType.setNeighboursOK(true);
				else if (dynamicMultiPathFitResult.getSuperDoubletFitResult() == selectedResult.fitResult)
					fitType.setDoubletOK(true);
			}
			add(fitType);
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
								params[j + Gaussian2DFunction.X_SD] / initialParams[j + Gaussian2DFunction.X_SD],
								params[j + Gaussian2DFunction.Y_SD] / initialParams[j + Gaussian2DFunction.Y_SD]);
					}
					break;

				case BAD_PARAMETERS:
					logger.info("Bad parameters: %s", Arrays.toString(fitResult.getInitialParameters()));
					break;

				default:
					logger.info(fitResult.getStatus().toString());
					break;
			}
		}

		// Debugging
		if (logger2 != null)
		{
			double[] peakParams = fitResult.getParameters();
			final Rectangle regionBounds = dynamicMultiPathFitResult.regionBounds;
			if (peakParams != null)
			{
				// Parameters are the raw values from fitting the region. Convert for logging.
				peakParams = Arrays.copyOf(peakParams, peakParams.length);
				int npeaks = peakParams.length / 6;
				for (int i = 0; i < npeaks; i++)
				{
					peakParams[i * 6 + Gaussian2DFunction.X_POSITION] += 0.5 + dataBounds.x + regionBounds.x;
					peakParams[i * 6 + Gaussian2DFunction.Y_POSITION] += 0.5 + dataBounds.y + regionBounds.y;
					peakParams[i * 6 + Gaussian2DFunction.ANGLE] *= 180.0 / Math.PI;
				}
			}
			final int x = candidates[candidateId].x;
			final int y = candidates[candidateId].y;
			logger2.debug("%d:%d [%d,%d] %s (%s) = %s\n", slice, candidateId, x + dataBounds.x, y + dataBounds.y,
					fitResult.getStatus(), fitResult.getStatusData(), Arrays.toString(peakParams));
		}

		// Check if there were any new results
		int npeaks = sliceResults.size() - currrentSize;
		if (npeaks != 0)
		{
			success++;

			// Support for post-processing filter 
			if (resultFilter != null)
			{
				// Check all result peaks for the distance to the filter positions
				PeakResult[] peakResults = new PeakResult[npeaks];
				for (int i = sliceResults.size(); npeaks-- > 0;)
				{
					peakResults[npeaks] = sliceResults.get(--i);
				}
				resultFilter.filter(fitResult, candidateId, peakResults);
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.MultiPathFilter.SelectedResultStore#isFit(int)
	 */
	public boolean isFit(int candidateId)
	{
		// Return if we already have a fit result for this candidate
		return candidates[candidateId].fit;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.MultiPathFilter.SelectedResultStore#isValid(int)
	 */
	public boolean isValid(int candidateId)
	{
		// If we have an estimate then this is a valid candidate for fitting
		return estimates[candidateId] != null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.MultiPathFilter.SelectedResultStore#pass(gdsc.smlm.results.filter.
	 * PreprocessedPeakResult)
	 */
	public void pass(PreprocessedPeakResult result)
	{
		// Ignore results that pass the primary filter. We will deal with these in add(SelectedResult) 
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.MultiPathFilter.SelectedResultStore#passMin(gdsc.smlm.results.filter.
	 * PreprocessedPeakResult)
	 */
	public void passMin(PreprocessedPeakResult result)
	{
		// This is a candidate that passed validation. Store the estimate as passing the minimal filter.
		storeEstimate(result.getCandidateId(), result, FILTER_RANK_MINIMAL);
	}

	/**
	 * Create a minimum filter to use for storing estimates
	 * 
	 * @return The minimal filter
	 */
	public static MultiFilter2 createMinimalFilter()
	{
		double signal = 30;
		float snr = 20;
		double minWidth = 0.5;
		double maxWidth = 4;
		double shift = 2;
		double eshift = 0;
		double precision = 60;
		return new MultiFilter2(signal, snr, minWidth, maxWidth, shift, eshift, precision);
	}
}
