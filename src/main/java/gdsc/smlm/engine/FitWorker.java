package gdsc.smlm.engine;

import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;
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
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.results.ExtendedPeakResult;
import gdsc.smlm.results.IdPeakResult;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResults;
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
	private CoordinateConverter cc;
	//private Rectangle regionBounds;
	private int border, borderLimitX, borderLimitY;
	private FitJob job;
	private boolean benchmarking;
	private float[] data;
	private DataEstimator dataEstimator = null;
	//private float[] filteredData;
	private boolean relativeIntensity;
	private float noise, background;
	private final boolean calculateNoise;
	private boolean estimateSignal;
	private ResultGridManager gridManager = null;
	private PeakResult[] peakNeighbours = null;
	// Contains the index in the list of maxima for any neighbours
	private int candidateNeighbourCount = 0;
	private Candidate[] candidateNeighbours = null;
	// Contains the index in the list of fitted results for any neighbours 
	private int fittedNeighbourCount = 0;
	private PeakResult[] fittedNeighbours = null;
	private float duplicateDistance2;

	private volatile boolean finished = false;

	static int WORKER_ID = 0;
	int workerId;

	private static byte FILTER_RANK_MINIMAL = (byte) 0;
	private static byte FILTER_RANK_PRIMARY = (byte) 1;

	/**
	 * Encapsulate all conversion of coordinates between the frame of data (data bounds) and the sub-section currently
	 * used in fitting (region bounds) and the global coordinate system.
	 */
	private class CoordinateConverter
	{
		/** The data bounds. */
		final Rectangle dataBounds;

		/** The region bounds. */
		Rectangle regionBounds;

		CoordinateConverter(Rectangle dataBounds)
		{
			this.dataBounds = dataBounds;
		}

		void setRegionBounds(Rectangle regionBounds)
		{
			this.regionBounds = regionBounds;
		}

		/**
		 * Convert from the data bounds to the global bounds.
		 *
		 * @param x
		 *            the x coordinate
		 * @return the x coordinate
		 */
		public int fromDataToGlobalX(int x)
		{
			return x + dataBounds.x;
		}

		/**
		 * Convert from the data bounds to the global bounds.
		 *
		 * @param y
		 *            the y coordinate
		 * @return the y coordinate
		 */
		public int fromDataToGlobalY(int y)
		{
			return y + dataBounds.y;
		}

		/**
		 * Convert from the region bounds to the global bounds.
		 *
		 * @param x
		 *            the x coordinate
		 * @return the x coordinate
		 */
		public int fromRegionToGlobalX(int x)
		{
			return x + dataBounds.x + regionBounds.x;
		}

		/**
		 * Convert from the region bounds to the global bounds.
		 *
		 * @param y
		 *            the y coordinate
		 * @return the y coordinate
		 */
		public int fromRegionToGlobalY(int y)
		{
			return y + dataBounds.y + regionBounds.y;
		}

		/**
		 * Conversion from the raw fit coordindates to the global bounds.
		 *
		 * @return the x offset
		 */
		public double fromFitRegionToGlobalX()
		{
			return 0.5 + dataBounds.x + regionBounds.x;
		}

		/**
		 * Conversion from the raw fit coordindates to the global bounds.
		 *
		 * @return the y offset
		 */
		public double fromFitRegionToGlobalY()
		{
			return 0.5 + dataBounds.y + regionBounds.y;
		}
	}

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

	private Estimate[] estimates = new Estimate[0], estimates2 = null;
	private boolean[] isValid = null;

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
		benchmarking = false;
		this.slice = job.slice;

		// Used for debugging
		//if (logger == null) logger = new gdsc.fitting.logging.ConsoleLogger();

		// Crop to the ROI
		cc = new CoordinateConverter(job.bounds);
		final int width = cc.dataBounds.width;
		final int height = cc.dataBounds.height;
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

		// TODO - Better estimate of the background and the noise. Using all the image pixels
		// results in an estimate that is too high when there are many spots in the image.
		// Create a method that thresholds the image and finds the mean/sd of the thresholded image.

		// Note: Other code calls the static estimateNoise method.
		// So add a private instance method to estimate the noise and background using a static helper
		// class. This can also be called from the static estiamteNoise method.

		// Always get the noise and store it with the results.
		if (params != null && !Float.isNaN(params.noise))
		{
			noise = params.noise;
			fitConfig.setNoise(noise);
		}
		else if (calculateNoise)
		{
			noise = estimateNoise(width, height);
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
				sliceResults.add(createResult(cc.fromDataToGlobalX(x), cc.fromDataToGlobalY(y), data[index], 0, noise,
						peakParams, null, n));
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

			// The SpotFitter is used to create a dynamic MultiPathFitResult object.
			// This is then passed to a multi-path filter. Thus the same fitting decision process 
			// is used when benchmarking and when running on actual data.

			// Note: The SpotFitter labels each PreprocessedFitResult using the offset in the FitResult object.
			// The initial params and deviations can then be extracted for the results that pass the filter.

			MultiPathFilter filter;
			IMultiPathFitResults multiPathResults = this;
			SelectedResultStore store = this;

			if (params != null && params.fitTask == FitTask.BENCHMARKING)
			{
				// Run filtering as normal. However in the event that a candidate is missed or some
				// results are not generated we must generate them. This is done in the complete(int)
				// method if we set the benchmarking flag.
				benchmarking = true;

				// Filter using the benchmark filter
				filter = params.benchmarkFilter;
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
			}
			else
			{
				// Filter using the configuration
				filter = new MultiPathFilter(fitConfig, createMinimalFilter(), config.getResidualsThreshold());
			}

			// If we are benchmarking then do not generate results dynamically since we will store all 
			// results in the fit job
			dynamicMultiPathFitResult = new DynamicMultiPathFitResult(ie, !benchmarking);
			//dynamicMultiPathFitResult = new DynamicMultiPathFitResult(ie, false);

			// Debug where the fit config may be different between benchmarking and fitting
			if (slice == -1)
			{
				SettingsManager.saveFitEngineConfiguration(config, String.format("/tmp/config.%b.xml", benchmarking));
				Utils.write(String.format("/tmp/filter.%b.xml", benchmarking), filter.toXML());
				//filter.setDebugFile(String.format("/tmp/fitWorker.%b.txt", benchmarking));
				StringBuilder sb = new StringBuilder();
				sb.append((benchmarking) ? ((gdsc.smlm.results.filter.Filter) filter.getFilter()).toXML()
						: fitConfig.getSmartFilter().toXML()).append("\n");
				sb.append(((gdsc.smlm.results.filter.Filter) filter.getMinimalFilter()).toXML()).append("\n");
				sb.append(filter.residualsThreshold).append("\n");
				sb.append(config.getFailuresLimit()).append("\n");
				if (spotFilter != null)
					sb.append(spotFilter.getDescription()).append("\n");
				for (int i = 0; i < candidates.length; i++)
				{
					sb.append(String.format("Fit %d [%d,%d = %.1f]\n", i, candidates[i].x, candidates[i].y,
							candidates[i].intensity));
				}
				Utils.write(String.format("/tmp/candidates.%b.xml", benchmarking), sb.toString());
			}

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
		final float offsetx = cc.dataBounds.x;
		final float offsety = cc.dataBounds.y;

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
			//filteredData = spotFilter.getPreprocessedData();
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

	private void finish(FitJob job, final long start)
	{
		time += System.nanoTime() - start;
		job.finished();
	}

	private void initialiseFitting()
	{
		candidateNeighbours = allocateArray(candidateNeighbours, candidates.length);
		// Allocate enough room for all fits to be doublets 
		fittedNeighbours = allocateArray(fittedNeighbours, candidates.length * 2);

		clearEstimates(candidates.length);

		gridManager = new ResultGridManager(cc.dataBounds.width, cc.dataBounds.height, 2 * fitting + 1);
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

	/**
	 * Queue to grid.
	 * <p>
	 * Used to add results to the grid for the current fit position. This prevents filtering duplicates within the
	 * current fit results, only with all that has been fit before.
	 *
	 * @param result
	 *            the result
	 */
	private void queueToGrid(PeakResult result)
	{
		if (queueSize == queue.length)
			queue = Arrays.copyOf(queue, queueSize * 2);
		queue[queueSize++] = result;
	}

	/**
	 * Flush the results for the current position to the grid.
	 */
	private void flushToGrid()
	{
		for (int i = 0; i < queueSize; i++)
			gridManager.putOnGrid(queue[i]);
		queueSize = 0;
	}

	/**
	 * Clear the grid cache of the local neighbourhood.
	 */
	private void clearGridCache()
	{
		gridManager.clearCache();
	}

	private Candidate[] allocateArray(Candidate[] array, int length)
	{
		if (array == null || array.length < length)
			array = new Candidate[length];
		return array;
	}

	private PeakResult[] allocateArray(PeakResult[] array, int length)
	{
		if (array == null || array.length < length)
			array = new PeakResult[length];
		return array;
	}

	private void clearEstimates(int length)
	{
		if (estimates.length < length)
		{
			estimates = new Estimate[length];
			estimates2 = new Estimate[length];
			isValid = new boolean[length];
		}
		else
		{
			for (int i = 0; i < length; i++)
			{
				estimates[i] = null;
				estimates2[i] = null;
				isValid[i] = false;
			}
		}
	}

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
	 * @return true, if successful
	 */
	private boolean addSingleResult(int candidateId, float[] peakParams, float[] peakParamsDev, double error)
	{
		int x = candidates[candidateId].x;
		int y = candidates[candidateId].y;
		float value = data[y * cc.dataBounds.width + x];

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
		x = cc.fromDataToGlobalX(x);
		y = cc.fromDataToGlobalY(y);

		// This was fit OK so add it to the grid of results (so we do not fit it again)
		final PeakResult peakResult = createResult(x, y, value, error, noise, peakParams, peakParamsDev, candidateId);
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
			float[] paramsStdDev, int id)
	{
		if (endT >= 0 && slice != endT)
		{
			return new ExtendedPeakResult(slice, origX, origY, origValue, error, noise, params, paramsStdDev, endT, id);
		}
		return new IdPeakResult(slice, origX, origY, origValue, error, noise, params, paramsStdDev, id);
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

		PreprocessedPeakResult createPreprocessedPeakResult(int candidateId, int n, double[] initialParams,
				double[] params, double localBackground, ResultType resultType)
		{
			// Update the initial params since we may have used an estimate
			// This will ensure that the width factor is computed correctly.

			// Q. Should this be ignored for existing results? They have already passed validation. 
			// So we do not have to be as strict on their width and could just use the drift from 
			// the the initial estimate.
			// For now do a full validation since multi-fit results are only accepted if existing
			// results are still valid.

			final int offset = n * 6;
			initialParams[Gaussian2DFunction.X_SD + offset] = fitConfig.getInitialPeakStdDev0();
			initialParams[Gaussian2DFunction.Y_SD + offset] = fitConfig.getInitialPeakStdDev1();
			return createResult(candidateId, n, initialParams, params, localBackground, resultType);
		}

		abstract PreprocessedPeakResult createResult(int candidateId, int n, double[] initialParams, double[] params,
				double localBackground, ResultType resultType);
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

		PreprocessedPeakResult createResult(int candidateId, int n, double[] initialParams, double[] params,
				double localBackground, ResultType resultType)
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

		PreprocessedPeakResult createResult(int candidateId, int n, double[] initialParams, double[] params,
				double localBackground, ResultType resultType)
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
		double singleBackground = Double.NaN, multiBackground = Double.NaN;

		MultiPathFitResult.FitResult resultMulti = null;
		boolean computedMulti = false;
		double[] multiResiduals = null;
		MultiPathFitResult.FitResult resultMultiDoublet = null;
		boolean computedMultiDoublet = false;
		QuadrantAnalysis multiQA = null;

		MultiPathFitResult.FitResult resultSingle = null;
		double[] singleRegion = null;
		double[] singleResiduals = null;
		double singleSumOfSquares = 0;
		double singleValue = 0;
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

			fitConfig.setFitRegion(width, height, 0.5);

			// Analyse neighbours and include them in the fit if they are within a set height of this peak.
			resetNeighbours();
			neighbours = findNeighbours(regionBounds, n);

			if (benchmarking)
			{
				// When benchmarking we may compute additional results after the main filtering routine
				// has been run. In this case we must pre-compute some values with the current results.

				// Cache the background estimate
				getSingleFittingBackground();
			}
		}

		private double getMultiFittingBackground()
		{
			if (Double.isNaN(multiBackground))
			{
				multiBackground = 0;
				if (fittedNeighbourCount > 0)
				{
					// Use the average previously fitted background

					// Add the details of the already fitted peaks
					for (int i = 0; i < fittedNeighbourCount; i++)
					{
						multiBackground += fittedNeighbours[i].getBackground();
					}

					multiBackground /= fittedNeighbourCount;
				}
				else
				{
					multiBackground = this.getSingleFittingBackground();
				}
				multiBackground = limitBackground(multiBackground);
			}
			return multiBackground;
		}

		private double getSingleFittingBackground()
		{
			if (Double.isNaN(singleBackground))
			{
				singleBackground = limitBackground(FitWorker.this.getSingleFittingBackground());
			}
			return singleBackground;
		}

		private double getDefaultBackground(double[] region, int width, int height)
		{
			// Use the minimum in the data.
			// This is what is done in the in the fitter if the background is zero.
			return limitBackground(Gaussian2DFitter.getBackground(region, width, height, 2));
		}

		private double limitBackground(double b)
		{
			// Ensure we do not get a negative background
			final double limit = (fitConfig.isRemoveBiasBeforeFitting()) ? fitConfig.getBias() : 0;
			if (b < limit)
				b = limit;
			return b;
		}

		private double getMax(double[] region, int width, int height)
		{
			double max = region[0];
			for (int i = width * height; --i > 0;)
				if (max < region[i])
					max = region[i];
			return max;
		}

		public MultiPathFitResult.FitResult getResultMulti()
		{
			if (computedMulti)
				return resultMulti;

			computedMulti = true;

			// Do not do a multi-fit if the configuration is not set to include neighbours 
			if (neighbours == 0 || !config.isIncludeNeighbours())
				return null;

			// Estimate background.
			// Note that using the background from previous fit results leads to an 
			// inconsistency in the results when the fitting parameters are changed which may be unexpected, 
			// e.g. altering the max iterations.
			double background = getMultiFittingBackground();

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
					final PeakResult result = fittedNeighbours[i];

					// Subtract peaks from the data if they are outside the fit region
					if (result.getXPosition() < xmin || result.getXPosition() > xmax || result.getYPosition() < ymin ||
							result.getYPosition() > ymax)
					{
						subtract[i] = true;
						subtractFittedPeaks++;
					}
				}
			}

			neighbours = candidateNeighbourCount + fittedNeighbourCount - subtractFittedPeaks;
			if (neighbours == 0)
			{
				// There are no neighbours after subtraction.
				// This will be the same result as the single result so return.
				return null;
			}

			// Create the parameters for the fit
			final int npeaks = 1 + neighbours;

			// Multiple-fit ...
			if (logger != null)
				logger.info("Slice %d: Multiple-fit (%d peaks : neighbours [%d + %d - %d])", slice, neighbours + 1,
						candidateNeighbourCount, fittedNeighbourCount, subtractFittedPeaks);

			double[] params = new double[1 + npeaks * parametersPerPeak];
			params[Gaussian2DFunction.BACKGROUND] = background;

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

			final boolean[] amplitudeEstimate = new boolean[npeaks];

			// The main peak. We use a close estimate if we have one.
			amplitudeEstimate[0] = getEstimate(candidates[candidateId], params, 0, true);

			// The neighbours
			for (int i = 0, j = parametersPerPeak; i < candidateNeighbourCount; i++, j += parametersPerPeak)
			{
				final Candidate candidateNeighbour = candidateNeighbours[i];
				amplitudeEstimate[i + 1] = getEstimate(candidateNeighbour, params, j, true);

				// Constrain the location using the candidate position.
				// Do not use the current estimate as this will create drift over time if the estimate is updated.
				final double candidateX = candidateNeighbour.x - regionBounds.x;
				final double candidateY = candidateNeighbour.y - regionBounds.y;
				lower[j + Gaussian2DFunction.X_POSITION] = candidateX - 1;
				upper[j + Gaussian2DFunction.X_POSITION] = candidateX + 1;
				lower[j + Gaussian2DFunction.Y_POSITION] = candidateY - 1;
				upper[j + Gaussian2DFunction.Y_POSITION] = candidateY + 1;
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
						PeakResult result = fittedNeighbours[i];
						// Copy Signal,Angle,Xpos,Ypos,Xwidth,Ywidth
						for (int k = 1; k <= parametersPerPeak; k++)
							funcParams[j + k] = result.params[k];
						// Adjust position relative to extracted region
						funcParams[j + Gaussian2DFunction.X_POSITION] -= xOffset;
						funcParams[j + Gaussian2DFunction.Y_POSITION] -= yOffset;
						j += parametersPerPeak;
					}

					final GaussianFunction func = fitConfig.createGaussianFunction(subtractFittedPeaks, width,
							funcParams);
					func.initialise(funcParams);

					// Subtract fitted peaks
					//double[] f = new double[region.length];
					for (int i = 0; i < region.length; i++)
					{
						//f[i] = func.eval(i);
						region[i] -= func.eval(i);
					}

					//gdsc.core.ij.Utils.display("Function2", f, width, height);
					//gdsc.core.ij.Utils.display("Region2", region, width, height);
				}

				// Add the details of the already fitted peaks
				for (int i = 0, j = (1 + candidateNeighbourCount) * parametersPerPeak; i < fittedNeighbourCount; i++)
				{
					if (subtract[i])
						continue;
					final PeakResult result = fittedNeighbours[i];

					// Copy Signal,Angle,Xpos,Ypos,Xwidth,Ywidth
					for (int k = 1; k <= parametersPerPeak; k++)
						params[j + k] = result.params[k];

					// Adjust position relative to extracted region
					params[j + Gaussian2DFunction.X_POSITION] -= xOffset;
					params[j + Gaussian2DFunction.Y_POSITION] -= yOffset;

					// Add support for constraining the known fit results using bounded coordinates.
					// Currently we just constrain the location.
					lower[j + Gaussian2DFunction.X_POSITION] = params[j + Gaussian2DFunction.X_POSITION] - 0.5;
					upper[j + Gaussian2DFunction.X_POSITION] = params[j + Gaussian2DFunction.X_POSITION] + 0.5;
					lower[j + Gaussian2DFunction.Y_POSITION] = params[j + Gaussian2DFunction.Y_POSITION] - 0.5;
					upper[j + Gaussian2DFunction.Y_POSITION] = params[j + Gaussian2DFunction.Y_POSITION] + 0.5;

					j += parametersPerPeak;
				}
			}

			// XXX Debugging the bad parameters 
			//double bbefore = background;
			//double[] before = params.clone();

			// In the case of a bad background estimate (e.g. uneven illumination) the peaks may 
			// be below the background.
			// Check the heights are positive.

			if (containsAmplitudeEstimates(amplitudeEstimate))
			{
				// Find the min amplitude
				double minSignal = Double.POSITIVE_INFINITY;
				for (int j = Gaussian2DFunction.SIGNAL, i = 0; j < params.length; j += parametersPerPeak, i++)
					if (amplitudeEstimate[i] && minSignal > params[j])
						minSignal = params[j];

				if (minSignal <= 0)
				{
					// Reset to the minimum value in the data.
					// This is what is done in the in the fitter if the background is zero.
					params[Gaussian2DFunction.BACKGROUND] = getDefaultBackground(region, width, height);
					final double backgroundChange = background - params[Gaussian2DFunction.BACKGROUND];

					// Make the peaks higher by the change in background
					for (int j = Gaussian2DFunction.SIGNAL, i = 0; j < params.length; j += parametersPerPeak, i++)
						if (amplitudeEstimate[i])
							params[j] += backgroundChange;

					if ((minSignal + backgroundChange) <= 0)
					{
						// This is probably extremely rare and the result of a poor candidate estimate.
						// Set a small height based on the data range
						final double defaultHeight = Math.max(1,
								0.1 * (getMax(region, width, height) - params[Gaussian2DFunction.BACKGROUND]));
						for (int j = Gaussian2DFunction.SIGNAL, i = 0; j < params.length; j += parametersPerPeak, i++)
							if (amplitudeEstimate[i] && params[j] <= 0)
								params[j] = defaultHeight;
					}
				}
			}

			// Note that if the input XY positions are on the integer grid then the fitter will estimate
			// the position using a CoM estimate. To avoid this adjust the centres to be off-grid.
			for (int i = Gaussian2DFunction.X_POSITION; i < params.length; i += parametersPerPeak)
			{
				for (int j = 0; j < 2; j++)
					if ((int) params[i + j] == params[i + j])
						// If at the limit then decrement instead
						params[i + j] += (params[i + j] == upper[i + j]) ? -0.001 : 0.001;
			}

			// -=-=-=-

			// Note: Some of the unfitted neighbours may be bad candidates.
			// Only validate the fitted neighbours and the current candidate.
			// The unfitted neighbours are allowed to fail.

			// Turn off validation of peaks
			final boolean smartFilter = fitConfig.isSmartFilter();
			final boolean disableSimpleFilter = fitConfig.isDisableSimpleFilter();
			fitConfig.setSmartFilter(false);
			fitConfig.setDisableSimpleFilter(true);
			fitConfig.setFitRegion(0, 0, 0);

			// Increase the iterations for a multiple fit.
			final int maxIterations = fitConfig.getMaxIterations();
			final int maxEvaluations = fitConfig.getMaxFunctionEvaluations();
			fitConfig.setMaxIterations(
					maxIterations + maxIterations * (npeaks - 1) * ITERATION_INCREASE_FOR_MULTIPLE_PEAKS);
			fitConfig.setMaxFunctionEvaluations(
					maxEvaluations + maxEvaluations * (npeaks - 1) * EVALUATION_INCREASE_FOR_MULTIPLE_PEAKS);

			gf.setBounds(lower, upper);
			final FitResult fitResult = gf.fit(region, width, height, npeaks, params, amplitudeEstimate,
					background == 0);
			multiResiduals = gf.getResiduals();
			gf.setBounds(null, null);

			//			if (fitResult.getStatus() == FitStatus.BAD_PARAMETERS)
			//			{
			//				int x = candidates[candidateId].x;
			//				int y = candidates[candidateId].y;
			//				int index = (y-regionBounds.y) * width + (x-regionBounds.x);
			//				System.out.printf("Bad : [%d,%d] %d,%d %.1f (%.1f) B=%.1f (%.1f) : %s\n", slice, candidateId,
			//						x, y, candidates[candidateId].intensity, region[index],
			//						bbefore, background, Arrays.toString(before));
			//				if (filteredData != null)
			//					gdsc.core.ij.Utils.display("Filtered", new FloatProcessor(dataBounds.width, dataBounds.height, filteredData));
			//			}

			// Restore
			fitConfig.setSmartFilter(smartFilter);
			fitConfig.setDisableSimpleFilter(disableSimpleFilter);
			fitConfig.setFitRegion(width, height, 0.5);
			fitConfig.setMaxIterations(maxIterations);
			fitConfig.setMaxFunctionEvaluations(maxEvaluations);

			updateError(fitResult);

			// Ensure the initial parameters are at the candidate position since we may have used an estimate.
			// This will ensure that drift is computed correctly.
			final double[] fitParams = fitResult.getParameters();
			final double[] initialParams = fitResult.getInitialParameters();

			initialParams[Gaussian2DFunction.X_POSITION] = candidates[candidateId].x - regionBounds.x;
			initialParams[Gaussian2DFunction.Y_POSITION] = candidates[candidateId].y - regionBounds.y;
			initialParams[Gaussian2DFunction.X_SD] = fitConfig.getInitialPeakStdDev0();
			initialParams[Gaussian2DFunction.Y_SD] = fitConfig.getInitialPeakStdDev1();

			// Perform validation of the candidate and existing peaks (other candidates are allowed to fail)
			if (fitResult.getStatus() == FitStatus.OK)
			{
				// The candidate peak
				if (fitConfig.validatePeak(0, initialParams, fitParams) != FitStatus.OK)
					return resultMulti = createResult(fitResult, null, fitConfig.getValidationResult());

				// Existing peaks
				for (int n = candidateNeighbourCount + 1; n < npeaks; n++)
				{
					if (fitConfig.validatePeak(n, initialParams, fitParams) != FitStatus.OK)
						return resultMulti = createResult(fitResult, null, fitConfig.getValidationResult());
				}
			}

			// Create the results
			PreprocessedPeakResult[] results = null;
			if (fitResult.getStatus() == FitStatus.OK)
			{
				// The primary candidate is not bounded. Check it has not drifted close to 
				// a neighbour. 

				// 3. Check we are not closer to a fitted spot. This has already had a chance at 
				//    fitting a doublet so is ignored.
				// 4. Check if there is an unfit candidate spot closer than the current candidate.
				//    This represents drift out to fit another spot that will be fit later.

				int otherId = candidateId;
				ResultType resultType = ResultType.NEW;

				final double xShift = fitParams[Gaussian2DFunction.X_POSITION] -
						initialParams[Gaussian2DFunction.X_POSITION];
				final double yShift = fitParams[Gaussian2DFunction.Y_POSITION] -
						initialParams[Gaussian2DFunction.Y_POSITION];

				// We must be closer to the current candidate than any other spots.
				// This is true for candidates we have yet to fit or already fitted candidates.

				// Distance to current candidate fitted as a single
				final double distanceToSingleFit2 = xShift * xShift + yShift * yShift;

				// 3. Check we are not closer to a fitted spot. This has already had a chance at 
				//    fitting a doublet so is ignored..

				final PeakResult[] peakNeighbours = findPeakNeighbours(candidates[candidateId]);
				if (peakNeighbours.length != 0)
				{
					// Coords for comparison to the real positions
					final float fcx2 = (float) (regionBounds.x + fitParams[Gaussian2DFunction.X_POSITION] + 0.5);
					final float fcy2 = (float) (regionBounds.y + fitParams[Gaussian2DFunction.Y_POSITION] + 0.5);
					final float d2 = (float) distanceToSingleFit2;
					int ii = -1;
					for (int i = 0; i < peakNeighbours.length; i++)
					{
						if (d2 > distance2(fcx2, fcy2, peakNeighbours[i].params))
						{
							// There is another fitted result that is closer.
							// Note: The fit region is not centred on the other spot so this fit will probably
							// be worse and is discarded (not compared to the existing fit to get the best one). 

							ii = i;
							otherId = peakNeighbours[i].getId();
						}
					}
					if (otherId != candidateId)
					{
						if (logger != null)
						{
							logger.info(
									"Bad peak: Fitted coordinates moved closer to another result (%d,%d : x=%.1f,y=%.1f : %.1f,%.1f)",
									candidates[candidateId].x, candidates[candidateId].y, fcx2, fcy2,
									peakNeighbours[ii].getXPosition(), peakNeighbours[ii].getYPosition());
						}
						//System.out.printf("Multi drift to another result: [%d,%d] %d\n", slice, candidateId, otherId);
						resultType = ResultType.EXISTING;

						// Update the initial parameters to the position of the existing result so 
						// that drift is correct for filtering
						initialParams[Gaussian2DFunction.X_POSITION] = peakNeighbours[ii].getXPosition() -
								cc.fromFitRegionToGlobalX();
						initialParams[Gaussian2DFunction.Y_POSITION] = peakNeighbours[ii].getYPosition() -
								cc.fromFitRegionToGlobalY();
					}
				}

				// 4. Check if there is an unfit candidate spot closer than the current candidate.
				//    This represents drift out to fit another unfitted spot.

				if (otherId != candidateId)
				{
					final CandidateList neighbours = findNeighbours(candidates[candidateId]);
					if (neighbours.size != 0)
					{
						// Position - do not add 0.5 pixel offset to allow distance comparison to integer candidate positions.
						float fcx2 = (float) (regionBounds.x + fitParams[Gaussian2DFunction.X_POSITION]);
						float fcy2 = (float) (regionBounds.y + fitParams[Gaussian2DFunction.Y_POSITION]);
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
										"Bad peak: Fitted coordinates moved closer to another candidate (%d,%d : x=%.1f,y=%.1f : %d,%d)",
										candidates[candidateId].x, candidates[candidateId].y, fcx2 + 0.5f, fcy2 + 0.5f,
										candidates[otherId].x, candidates[otherId].y);
							}
							//System.out.printf("Multi drift to another candidate: [%d,%d] %d\n", slice, candidateId, otherId);

							// There is another candidate to be fit later that is closer.
							// This may be used as an estimate so we return it as such (i.e we do not ignore it)
							//otherId = candidateId;
							if (otherId > candidateId)
								resultType = ResultType.CANDIDATE;

							// Update the initial parameters to the position of the candidate so 
							// that drift is correct for filtering
							initialParams[Gaussian2DFunction.X_POSITION] = candidates[otherId].x - regionBounds.x;
							initialParams[Gaussian2DFunction.Y_POSITION] = candidates[otherId].y - regionBounds.y;
						}
					}
				}

				convertParameters(fitParams);

				results = new PreprocessedPeakResult[npeaks];

				// We must compute a local background for all the spots
				final int flags = fitConfig.getFunctionFlags();

				// Note: This could be the current candidate or drift to another candidate
				results[0] = resultFactory.createPreprocessedPeakResult(otherId, 0, initialParams, fitParams,
						getLocalBackground(0, npeaks, fitParams, flags), resultType);

				// Neighbours
				int n = 1;
				for (int i = 0; i < candidateNeighbourCount; i++)
				{
					final Candidate candidateNeighbour = candidateNeighbours[i];
					results[n] = resultFactory.createPreprocessedPeakResult(candidateNeighbour.index, n, initialParams,
							fitParams, getLocalBackground(n, npeaks, fitParams, flags), ResultType.CANDIDATE);
					n++;
				}

				// Already fitted peaks
				for (int i = 0; i < fittedNeighbourCount; i++)
				{
					if (subtract[i])
						continue;
					final PeakResult result = fittedNeighbours[i];
					results[n] = resultFactory.createPreprocessedPeakResult(result.getId(), n, initialParams, fitParams,
							getLocalBackground(n, npeaks, fitParams, flags), ResultType.EXISTING);
					n++;
				}
			}

			return resultMulti = createResult(fitResult, results);
		}

		private boolean containsAmplitudeEstimates(boolean[] amplitudeEstimate)
		{
			for (int i = 0; i < amplitudeEstimate.length; i++)
				if (amplitudeEstimate[i])
					return true;
			return false;
		}

		private double getLocalBackground(int n, int npeaks, double[] params, final int flags)
		{
			GaussianOverlapAnalysis overlap = new GaussianOverlapAnalysis(flags, extractSpotParams(params, n), 2);
			overlap.add(flags, extractOtherParams(params, n, npeaks), true);
			double[] overlapData = overlap.getOverlapData();
			return overlapData[1] + params[0];
		}

		private boolean getEstimate(Candidate candidate, double[] params, int j, boolean close)
		{
			final double[] estimatedParams = getEstimate(candidate.index, close);
			if (estimatedParams != null)
			{
				// Re-use previous good multi-fit results to estimate the peak params...
				params[j + Gaussian2DFunction.SIGNAL] = estimatedParams[Gaussian2DFunction.SIGNAL];
				params[j + Gaussian2DFunction.X_POSITION] = estimatedParams[Gaussian2DFunction.X_POSITION] -
						regionBounds.x;
				params[j + Gaussian2DFunction.Y_POSITION] = estimatedParams[Gaussian2DFunction.Y_POSITION] -
						regionBounds.y;
				params[j + Gaussian2DFunction.ANGLE] = estimatedParams[Gaussian2DFunction.ANGLE];
				params[j + Gaussian2DFunction.X_SD] = estimatedParams[Gaussian2DFunction.X_SD];
				params[j + Gaussian2DFunction.Y_SD] = estimatedParams[Gaussian2DFunction.Y_SD];
				return false;
			}
			else
			{
				// Amplitude estimate
				params[j + Gaussian2DFunction.SIGNAL] = candidate.intensity -
						((relativeIntensity) ? 0 : params[Gaussian2DFunction.BACKGROUND]);
				params[j + Gaussian2DFunction.X_POSITION] = candidate.x - regionBounds.x;
				params[j + Gaussian2DFunction.Y_POSITION] = candidate.y - regionBounds.y;
				return true;
			}
		}

		/**
		 * Gets the estimate. Note estimates are classed as close (within 1 pixel) of the candidate position, or not. A
		 * candidate may have either or both types of estimate. The close estimate is used in preference to the other.
		 *
		 * @param i
		 *            the candidate index
		 * @param close
		 *            True if only considering the close estimate
		 * @return the estimate
		 */
		private double[] getEstimate(int i, boolean close)
		{
			// Check the close estimate
			if (estimates[i] != null)
				return estimates[i].params;

			// Only return the second estimate if we do not require the close estimate
			return (close || estimates2[i] == null) ? null : estimates2[i].params;
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

		public double getMultiQAScore()
		{
			if (multiQA != null)
				return multiQA.score;

			// Ensure we have a multi result
			getResultMulti();

			// Note this assumes that this method will be called after a multi fit and that the 
			// residuals were computed.
			final double[] residuals = multiResiduals;
			multiQA = computeQA((residuals == null) ? null : (FitResult) resultMulti.data, regionBounds, residuals);

			return multiQA.score;
		}

		public MultiPathFitResult.FitResult getResultMultiDoublet(double residualsThreshold)
		{
			// Fit a multi-fit. If all peaks are valid then subtract them from the image, apart from the central candidate,
			// then perform residuals analysis and fit as a doublet. Validate against an updated neighbourhood using the 
			// multi-fit results

			if (computedMultiDoublet)
				return resultMultiDoublet;

			if (residualsThreshold >= 1 || residualsThreshold < 0)
				return null;

			if (getMultiQAScore() < residualsThreshold)
				return null;

			computedMultiDoublet = true;

			if (resultMulti.status != 0)
				return null;

			// Note this assumes that this method will be called after a multi fit and that the 
			// residuals were computed.
			final double[] residuals = multiResiduals;
			if (residuals == null)
				return null;

			// Ideally all multi-results must be valid to proceed with doublet fitting. 
			// We do not perform validation of the results. So we assume that the results have 
			// been checked and are valid and continue.

			PreprocessedPeakResult[] fitResults = resultMulti.results;

			// Get the background for the multi-fit result
			final FitResult multiFitResult = (FitResult) resultMulti.data;
			final double[] fittedParams = multiFitResult.getParameters();
			final float background = (float) fittedParams[0];

			// Get the neighbours 
			final CandidateList neighbours = findNeighbours(candidates[candidateId]).copy();

			// Exclude the fitted candidate neighbours from the candidate neighbours
			int size = 0;
			NEXT_NEIGHBOUR: for (int i = 0; i < neighbours.size; i++)
			{
				final int otherId = neighbours.list[i].index;
				for (int j = 0; j < fitResults.length; j++)
				{
					if (fitResults[j].getCandidateId() == otherId)
					{
						continue NEXT_NEIGHBOUR;
					}
				}
				neighbours.list[size++] = neighbours.list[i];
			}
			neighbours.size = size;

			// Get the fitted neighbours
			final PeakResult[] peakNeighbours = findPeakNeighbours(candidates[candidateId]);
			size = peakNeighbours.length;
			PeakResult[] peakNeighbours2 = Arrays.copyOf(peakNeighbours, size + fitResults.length);

			// Update with the fitted results from the multi fit
			NEXT_RESULT: for (int j = 0; j < fitResults.length; j++)
			{
				final int otherId = fitResults[j].getCandidateId();
				if (otherId == candidateId)
					// Ignore this as it is the current candidate
					continue;
				// Check if this is already a fitted neighbour
				for (int i = 0; i < peakNeighbours.length; i++)
				{
					if (otherId == peakNeighbours[i].getId())
					{
						// Options:
						// 1, Change the coordinates to those of the multi-fit
						// 2. Add the new multi-fit result and keep the old result (leaving 2 positions for the neighbour)
						// 3. Leave the position as that originally fitted.

						// Choose option 3 for simplicity. Note that the original fitted coordinates
						// should be more accurate as the fit region was centred around the spot. Also
						// note that due to bounding the fit coordinates of any existing spot will be 
						// limited to a single pixel shift in XY.
						continue NEXT_RESULT;
					}
				}
				// Create a new fitted neighbour
				// (Use similar logic to when we create the actual results in #add(SelectedResult))
				final double[] p = fitResults[j].toGaussian2DParameters();
				final float[] params = new float[p.length];
				params[Gaussian2DFunction.BACKGROUND] = background;
				for (int i = 1; i < params.length; i++)
					params[i] = (float) p[i];
				// Store slice results relative to the data frame (not the global bounds)
				// Convert back so that 0,0 is the top left of the data bounds
				params[Gaussian2DFunction.X_POSITION] -= cc.dataBounds.x;
				params[Gaussian2DFunction.Y_POSITION] -= cc.dataBounds.y;
				peakNeighbours2[size++] = FitWorker.this.createResult(0, 0, 0, 0, 0, params, null, otherId);
			}

			peakNeighbours2 = Arrays.copyOf(peakNeighbours2, size);

			// Add the single fitted peak to the residuals. This is the data to fit as a doublet.
			final double[] region = residuals.clone();
			GaussianFunction func = fitConfig.createGaussianFunction(1, width, fittedParams);
			func.initialise(fittedParams);
			for (int i = 0; i < region.length; i++)
			{
				region[i] += func.eval(i);
			}

			// Build a fit result object for fitting the single candidate to the data.
			// This will be a single evaluation of the function against the data.
			int nFittedParameters = func.gradientIndices().length;
			int degreesOfFreedom = FastMath.max(region.length - nFittedParameters, 0);
			double error = 0;
			double[] initialParameters = null;
			double[] parameters = Arrays.copyOf(fittedParams, 7);
			double[] parametersDev = null;
			int nPeaks = 1;
			Object data = null;
			int iterations = 0;
			int evaluations = 0;
			FitResult fitResult = new FitResult(FitStatus.OK, degreesOfFreedom, error, initialParameters, parameters,
					parametersDev, nPeaks, nFittedParameters, data, iterations, evaluations);

			// Evaluate the multi fit as if fitted as a single peak. 
			// The resulting sum-of-squares and function value are used in the doublet fit
			try
			{
				gf.setComputeResiduals(false);
				if (!gf.evaluate(region, regionBounds.width, regionBounds.height, 1, parameters))
					return null;
			}
			finally
			{
				gf.setComputeResiduals(true);
			}

			double singleSumOfSquares = gf.getFinalResidualSumOfSquares();
			double singleValue = gf.getValue();

			//          // Debugging:
			//			// The evaluate computes the residuals. These should be similar to the original residuals
			//			double[] residuals2 = gf.getResiduals();
			//			for (int i = 0; i < region.length; i++)
			//				if (!gdsc.core.utils.DoubleEquality.almostEqualRelativeOrAbsolute(residuals[i], residuals2[i], 1e-5,
			//						1e-6))
			//				{
			//					System.out.printf("Residuals error [%d] %f != %f\n", i, residuals[i], residuals2[i]);
			//					gdsc.core.ij.Utils.display("Residuals1", residuals, width, height);
			//					gdsc.core.ij.Utils.display("Residuals2", residuals2, width, height);
			//					gdsc.core.ij.Utils.display("Region", region, width, height);
			//					break;
			//				}

			resultMultiDoublet = fitAsDoublet(fitResult, region, regionBounds, residualsThreshold, neighbours,
					peakNeighbours2, multiQA, singleSumOfSquares, singleValue);

			//			if (resultMultiDoublet != null && resultMultiDoublet.status == FitStatus.BAD_PARAMETERS.ordinal())
			//			{
			//				System.out.println("Bad params: " + Arrays.toString(parameters));
			//				//gdsc.core.ij.Utils.display("Region", region, width, height);
			//				//gdsc.core.ij.Utils.display("Residuals1", residuals, width, height);
			//			}

			if (resultMultiDoublet != null && resultMultiDoublet.status == 0 && resultMultiDoublet.results != null)
			{
				// The code below builds a combined result for the multi-fit and the primary candidate
				// fitted as a doublet. However this means that validation will be repeated on the spots that have 
				// been validated before. This is a small overhead but allows using the multi-doublet result to 
				// return all the fit results combined.

				// Fitting 2 spots is better than 1 and does not clash with any of the other multi-fit results.
				// Build a combined result with the doublet and the other multi-fit results
				FitResult doubletFitResult = (gdsc.smlm.fitting.FitResult) resultMultiDoublet.getData();

				double[] initialParams1 = doubletFitResult.getInitialParameters();
				double[] params1 = doubletFitResult.getParameters();
				double[] paramsDev1 = doubletFitResult.getParameterStdDev();

				double[] initialParams2 = multiFitResult.getInitialParameters();
				double[] params2 = multiFitResult.getParameters();
				double[] paramsDev2 = multiFitResult.getParameterStdDev();

				// Create the initial parameters by adding an extra spot
				double[] initialParams = new double[initialParams2.length + 6];
				double[] params = new double[initialParams.length];
				double[] paramsDev = (paramsDev2 == null) ? null : new double[initialParams.length];

				System.arraycopy(initialParams1, 0, initialParams, 0, initialParams1.length);
				System.arraycopy(initialParams2, 7, initialParams, initialParams1.length, initialParams2.length - 7);

				System.arraycopy(params1, 0, params, 0, params1.length);
				System.arraycopy(params2, 7, params, params1.length, params2.length - 7);

				if (paramsDev != null)
				{
					System.arraycopy(paramsDev1, 0, paramsDev, 0, paramsDev1.length);
					System.arraycopy(paramsDev2, 7, paramsDev, paramsDev1.length, paramsDev2.length - 7);
				}

				// Create all the output results
				PreprocessedPeakResult[] results = new PreprocessedPeakResult[resultMultiDoublet.results.length +
						resultMulti.results.length - 1];

				final int npeaks = multiFitResult.getNumberOfPeaks() + 1;
				// We must compute a local background for all the spots
				final int flags = fitConfig.getFunctionFlags();

				int n = 0;
				for (int i = 0; i < resultMultiDoublet.results.length; i++)
				{
					PreprocessedPeakResult r = resultMultiDoublet.results[i];
					results[n] = resultFactory.createPreprocessedPeakResult(r.getCandidateId(), r.getId(),
							initialParams, params, getLocalBackground(n, npeaks, params, flags), (r.isExistingResult())
									? ResultType.EXISTING : (r.isNewResult()) ? ResultType.NEW : ResultType.CANDIDATE);
					n++;
				}
				// Ignore the first result (this was removed and fit as the doublet)
				for (int i = 1; i < resultMulti.results.length; i++)
				{
					PreprocessedPeakResult r = resultMulti.results[i];
					// Increment the ID by one since the position in the parameters array is moved to 
					// accommodate 2 preceding peaks and not 1 
					results[n] = resultFactory.createPreprocessedPeakResult(r.getCandidateId(), r.getId() + 1,
							initialParams, params, getLocalBackground(n, npeaks, params, flags), (r.isExistingResult())
									? ResultType.EXISTING : (r.isNewResult()) ? ResultType.NEW : ResultType.CANDIDATE);
					n++;
				}

				final int adjust = (fitConfig.isBackgroundFitting()) ? 1 : 0;
				nFittedParameters = npeaks *
						((multiFitResult.getNumberOfFittedParameters() - adjust) / multiFitResult.getNumberOfPeaks()) +
						adjust;
				degreesOfFreedom = FastMath.max(region.length - nFittedParameters, 0);
				// Use error from the doublet fit
				error = doubletFitResult.getError();
				data = null;
				iterations = doubletFitResult.getIterations();
				evaluations = doubletFitResult.getEvaluations();

				FitResult mdoubletFitResult = new FitResult(FitStatus.OK, degreesOfFreedom, error, initialParams,
						params, paramsDev, npeaks, nFittedParameters, data, iterations, evaluations);

				resultMultiDoublet = createResult(mdoubletFitResult, results);
			}
			else
			{
				// Nothing to return. Do not set to null to allow reporting of the errors
				//resultMultiDoublet = null;
			}

			return resultMultiDoublet;
		}

		public MultiPathFitResult.FitResult getResultSingle()
		{
			if (resultSingle != null)
				return resultSingle;

			double background = getMultiFittingBackground();
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
					final PeakResult result = fittedNeighbours[i];
					// Copy Signal,Angle,Xpos,Ypos,Xwidth,Ywidth
					for (int k = 1; k <= parametersPerPeak; k++)
						funcParams[j + k] = result.params[k];
					// Adjust position relative to extracted region
					funcParams[j + Gaussian2DFunction.X_POSITION] -= xOffset;
					funcParams[j + Gaussian2DFunction.Y_POSITION] -= yOffset;
					j += parametersPerPeak;
				}

				GaussianFunction func = fitConfig.createGaussianFunction(fittedNeighbourCount, width, funcParams);
				func.initialise(funcParams);

				for (int i = 0; i < region.length; i++)
				{
					region[i] -= func.eval(i);
				}
			}

			// Store this as it is used in doublet fitting
			this.singleRegion = region;

			final double[] params = new double[] { background, 0, 0, 0, 0, 0, 0 };

			final boolean[] amplitudeEstimate = new boolean[1];

			// Re-use an estimate if we have it. Note that this may be quite far from the candidate.
			amplitudeEstimate[0] = getEstimate(candidates[candidateId], params, 0, false);
			boolean usingEstimate = getEstimate(candidates[candidateId].index, false) != null;
			if (!usingEstimate)
			{
				// If we have no estimate the default will be an amplitude estimate.
				// We can estimate the signal here instead of using the amplitude.
				// Do this when the fitting window covers enough of the Gaussian (e.g. 2.5xSD).
				float signal = 0;
				if (estimateSignal)
				{
					double sum = 0;
					final int size = width * height;
					for (int i = size; i-- > 0;)
						sum += region[i];
					signal = (float) (sum - background * size);
				}
				if (signal > 0)
				{
					params[Gaussian2DFunction.SIGNAL] = signal;
				}
				else
				{
					// Resort to default amplitude estimate
					amplitudeEstimate[0] = true;

					if (params[Gaussian2DFunction.SIGNAL] <= 0)
					{
						// Reset to the minimum value in the data.
						params[Gaussian2DFunction.BACKGROUND] = getDefaultBackground(region, width, height);
						final double backgroundChange = background - params[Gaussian2DFunction.BACKGROUND];

						params[Gaussian2DFunction.SIGNAL] += backgroundChange;

						if (params[Gaussian2DFunction.SIGNAL] <= 0)
						{
							// This is probably extremely rare and the result of a poor candidate estimate.
							// Set a small height based on the data range
							final double defaultHeight = Math.max(1,
									0.1 * (getMax(region, width, height) - params[Gaussian2DFunction.BACKGROUND]));
							params[Gaussian2DFunction.SIGNAL] = defaultHeight;
						}
					}
				}
			}

			// If there were unfitted neighbours or we have an estimate 
			// then use off grid pixels to prevent re-estimate of CoM
			// (since the CoM estimate will be skewed by the neighbours, or is not needed)
			if (candidateNeighbourCount != 0 || usingEstimate)
			{
				if ((int) params[Gaussian2DFunction.X_POSITION] == params[Gaussian2DFunction.X_POSITION])
					params[Gaussian2DFunction.X_POSITION] += 0.001;
				if ((int) params[Gaussian2DFunction.Y_POSITION] == params[Gaussian2DFunction.Y_POSITION])
					params[Gaussian2DFunction.Y_POSITION] += 0.001;
			}

			final FitResult fitResult = gf.fit(region, width, height, 1, params, amplitudeEstimate,
					params[Gaussian2DFunction.BACKGROUND] == 0);
			singleResiduals = gf.getResiduals();
			singleSumOfSquares = gf.getFinalResidualSumOfSquares();
			singleValue = gf.getValue();
			updateError(fitResult);

			// Ensure the initial parameters are at the candidate position since we may have used an estimate.
			// This will ensure that drift is computed correctly.
			final double[] initialParams = fitResult.getInitialParameters();
			initialParams[Gaussian2DFunction.X_POSITION] = candidates[candidateId].x - regionBounds.x;
			initialParams[Gaussian2DFunction.Y_POSITION] = candidates[candidateId].y - regionBounds.y;

			// Create the results
			PreprocessedPeakResult[] results = null;
			if (fitResult.getStatus() == FitStatus.OK)
			{
				// The primary candidate is not bounded. Check it has not drifted close to 
				// a neighbour. 

				// 3. Check we are not closer to a fitted spot. This has already had a chance at 
				//    fitting a doublet so is ignored.
				// 4. Check if there is an unfit candidate spot closer than the current candidate.
				//    This represents drift out to fit another spot that will be fit later.

				final double[] fitParams = fitResult.getParameters();

				int otherId = candidateId;
				ResultType resultType = ResultType.NEW;

				final double xShift = fitParams[Gaussian2DFunction.X_POSITION] -
						initialParams[Gaussian2DFunction.X_POSITION];
				final double yShift = fitParams[Gaussian2DFunction.Y_POSITION] -
						initialParams[Gaussian2DFunction.Y_POSITION];

				// We must be closer to the current candidate than any other spots.
				// This is true for candidates we have yet to fit or already fitted candidates.

				// Distance to current candidate fitted as a single
				final double distanceToSingleFit2 = xShift * xShift + yShift * yShift;

				// 3. Check we are not closer to a fitted spot. This has already had a chance at 
				//    fitting a doublet so is ignored..

				final PeakResult[] peakNeighbours = findPeakNeighbours(candidates[candidateId]);
				if (peakNeighbours.length != 0)
				{
					// Coords for comparison to the real positions
					final float fcx2 = (float) (regionBounds.x + fitParams[Gaussian2DFunction.X_POSITION] + 0.5);
					final float fcy2 = (float) (regionBounds.y + fitParams[Gaussian2DFunction.Y_POSITION] + 0.5);
					final float d2 = (float) distanceToSingleFit2;
					int ii = -1;
					for (int i = 0; i < peakNeighbours.length; i++)
					{
						if (d2 > distance2(fcx2, fcy2, peakNeighbours[i].params))
						{
							// There is another fitted result that is closer.
							// Note: The fit region is not centred on the other spot so this fit will probably
							// be worse and is discarded (not compared to the existing fit to get the best one). 

							ii = i;
							otherId = peakNeighbours[i].getId();
						}
					}
					if (otherId != candidateId)
					{
						if (logger != null)
						{
							logger.info(
									"Bad peak: Fitted coordinates moved closer to another result (%d,%d : x=%.1f,y=%.1f : %.1f,%.1f)",
									candidates[candidateId].x, candidates[candidateId].y, fcx2, fcy2,
									peakNeighbours[ii].getXPosition(), peakNeighbours[ii].getYPosition());
						}
						//System.out.printf("Single drift to another result: [%d,%d] %d\n", slice, candidateId, otherId);						
						resultType = ResultType.EXISTING;

						// Update the initial parameters to the position of the existing result so 
						// that drift is correct for filtering
						initialParams[Gaussian2DFunction.X_POSITION] = peakNeighbours[ii].getXPosition() -
								cc.fromFitRegionToGlobalX();
						initialParams[Gaussian2DFunction.Y_POSITION] = peakNeighbours[ii].getYPosition() -
								cc.fromFitRegionToGlobalY();
					}
				}

				// 4. Check if there is an unfit candidate spot closer than the current candidate.
				//    This represents drift out to fit another unfitted spot.

				if (otherId != candidateId)
				{
					final CandidateList neighbours = findNeighbours(candidates[candidateId]);
					if (neighbours.size != 0)
					{
						// Position - do not add 0.5 pixel offset to allow distance comparison to integer candidate positions.
						float fcx2 = (float) (regionBounds.x + fitParams[Gaussian2DFunction.X_POSITION]);
						float fcy2 = (float) (regionBounds.y + fitParams[Gaussian2DFunction.Y_POSITION]);
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
										"Bad peak: Fitted coordinates moved closer to another candidate (%d,%d : x=%.1f,y=%.1f : %d,%d)",
										candidates[candidateId].x, candidates[candidateId].y, fcx2 + 0.5f, fcy2 + 0.5f,
										candidates[otherId].x, candidates[otherId].y);
							}
							//System.out.printf("Single drift to another candidate: [%d,%d] %d\n", slice, candidateId, otherId);

							// There is another candidate to be fit later that is closer.
							// This may be used as an estimate so we return it as such (i.e we do not ignore it)
							//otherId = candidateId;
							if (otherId > candidateId)
								resultType = ResultType.CANDIDATE;

							// Update the initial parameters to the position of the candidate so 
							// that drift is correct for filtering
							initialParams[Gaussian2DFunction.X_POSITION] = candidates[otherId].x - regionBounds.x;
							initialParams[Gaussian2DFunction.Y_POSITION] = candidates[otherId].y - regionBounds.y;
						}
					}
				}

				convertParameters(fitParams);

				results = new PreprocessedPeakResult[1];

				results[0] = resultFactory.createPreprocessedPeakResult(otherId, 0, initialParams, fitParams, 0,
						resultType);
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
			return createResult(fitResult, results, fitResult.getStatus());
		}

		private MultiPathFitResult.FitResult createResult(FitResult fitResult, PreprocessedPeakResult[] results,
				FitStatus fitStatus)
		{
			MultiPathFitResult.FitResult mfitResult = new MultiPathFitResult.FitResult(fitStatus.ordinal(), fitResult);
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
			final double[] residuals = singleResiduals;
			singleQA = computeQA((residuals == null) ? null : (FitResult) resultSingle.data, regionBounds, residuals);

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
					neighbours, peakNeighbours, singleQA, singleSumOfSquares, singleValue);

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

			// Q. The primary candidate may have drifted. Should we check it is reasonably centred in the region?. 

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
		 * @param residualsThreshold
		 *            the residuals threshold
		 * @param neighbours
		 *            the neighbours
		 * @param peakNeighbours
		 *            the peak neighbours
		 * @param qa
		 *            the qa object that performed quadrant analysis
		 * @return the multi path fit result. fit result
		 */
		private MultiPathFitResult.FitResult fitAsDoublet(FitResult fitResult, double[] region, Rectangle regionBounds,
				double residualsThreshold, CandidateList neighbours, PeakResult[] peakNeighbours, QuadrantAnalysis qa,
				double singleSumOfSquares, double singelValue)
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

			// If validation is on in the fit configuration:
			// - Disable checking of position movement since we guessed the 2-peak location.
			// - Disable checking within the fit region as we do that per peak 
			//   (which is better than failing if either peak is outside the region)
			// - Increase the iterations level then reset afterwards.

			// TODO - Should width and signal validation be disabled too?
			final double shift = fitConfig.getCoordinateShift();
			final int maxIterations = fitConfig.getMaxIterations();
			final int maxEvaluations = fitConfig.getMaxFunctionEvaluations();

			fitConfig.setCoordinateShift(FastMath.min(width, height));
			fitConfig.setFitRegion(0, 0, 0);
			fitConfig.setMaxIterations(maxIterations * ITERATION_INCREASE_FOR_DOUBLETS);
			fitConfig.setMaxFunctionEvaluations(maxEvaluations * FitWorker.EVALUATION_INCREASE_FOR_DOUBLETS);

			// We assume that residuals calculation is on but just in case something else turned it off we get the state.
			final boolean isComputeResiduals = gf.isComputeResiduals();
			gf.setComputeResiduals(false);
			final boolean[] amplitudeEstimate = new boolean[2];
			final FitResult newFitResult = gf.fit(region, width, height, 2, doubletParams, amplitudeEstimate,
					doubletParams[Gaussian2DFunction.BACKGROUND] == 0);
			gf.setComputeResiduals(isComputeResiduals);

			fitConfig.setCoordinateShift(shift);
			fitConfig.setFitRegion(width, height, 0.5);
			fitConfig.setMaxIterations(maxIterations);
			fitConfig.setMaxFunctionEvaluations(maxEvaluations);

			updateError(newFitResult);

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
							cc.fromRegionToGlobalX(cx), cc.fromRegionToGlobalY(cy), newFitResult.getStatus(),
							newFitResult.getStatusData(), singleValue, gf.getValue(), singleSumOfSquares,
							doubleSumOfSquares, ic1, ic2, Arrays.toString(peakParams));
					logger2.debug(msg);
				}

				// Check if the predictive power of the model is better with two peaks:
				// IC should be lower
				if (ic2 > ic1)
				{
					return createResult(newFitResult, null, FitStatus.NO_MODEL_IMPROVEMENT);
				}

				// Validation of fit. For each spot:
				// 1. Check the spot is inside the region
				// 2. Check the distance of each new centre from the original centre.
				//    If the shift is too far (e.g. half the distance to the edge), the centre 
				//    must be in the correct quadrant. 
				// 3. Check we are not closer to a fitted spot. This has already had a chance at 
				//    fitting a doublet so is ignored.
				// 4. Check if there is an unfit candidate spot closer than the current candidate.
				//    This represents drift out to fit another spot that will be fit later.

				final double[] fitParams = newFitResult.getParameters();
				final double[] initialParams = newFitResult.getInitialParameters();

				// Allow the shift to span half of the fitted window.
				final double halfWindow = 0.5 * FastMath.min(regionBounds.width, regionBounds.height);

				final int[] position = new int[2];
				final int[] candidateIndex = new int[2];
				int nPeaks = 0;
				NEXT_PEAK: for (int n = 0; n < 2; n++)
				{
					final int offset = n * 6;

					// Ensure the initial parameters are at the candidate position since we may have used an estimate.
					// This will ensure that drift is computed correctly.
					initialParams[Gaussian2DFunction.X_POSITION + offset] = candidates[candidateId].x - regionBounds.x;
					initialParams[Gaussian2DFunction.Y_POSITION + offset] = candidates[candidateId].y - regionBounds.y;

					// 1. Check the spot is inside the region

					// Note that during processing the data is assumed to refer to the top-left
					// corner of the pixel. The coordinates should be represented in the middle of the pixel 
					// so add a 0.5 shift to the coordinates.
					final double xpos = fitParams[Gaussian2DFunction.X_POSITION + offset] + 0.5;
					final double ypos = fitParams[Gaussian2DFunction.Y_POSITION + offset] + 0.5;

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

					final double xShift = fitParams[Gaussian2DFunction.X_POSITION + offset] -
							params[Gaussian2DFunction.X_POSITION];
					final double yShift = fitParams[Gaussian2DFunction.Y_POSITION + offset] -
							params[Gaussian2DFunction.Y_POSITION];
					if (Math.abs(xShift) > halfWindow || Math.abs(yShift) > halfWindow)
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
						final float fcx2 = (float) (regionBounds.x + fitParams[Gaussian2DFunction.X_POSITION + offset] +
								0.5);
						final float fcy2 = (float) (regionBounds.y + fitParams[Gaussian2DFunction.Y_POSITION + offset] +
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
								// Q. Should this be returned for validation? 
								// A. Currently the MultiPathFilter accepts any new result and all existing results must pass.
								// So we could return this as an existing result which would make validation tougher.
								// However the existing result should have been subtracted from the input data so it will not 
								// be a full peak making validation incorrect. So at the moment we ignore this result.
								// Note that any new result will still have to be valid.
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
						float fcx2 = (float) (regionBounds.x + fitParams[Gaussian2DFunction.X_POSITION + offset]);
						float fcy2 = (float) (regionBounds.y + fitParams[Gaussian2DFunction.Y_POSITION + offset]);
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

							// Update the initial parameters to the position of the candidate so 
							// that drift is correct for filtering
							initialParams[Gaussian2DFunction.X_POSITION + offset] = candidates[otherId].x -
									regionBounds.x;
							initialParams[Gaussian2DFunction.Y_POSITION + offset] = candidates[otherId].y -
									regionBounds.y;
						}
					}

					candidateIndex[nPeaks] = otherId;
					position[nPeaks++] = n;
				}

				if (nPeaks == 0)
				{
					return createResult(newFitResult, null, FitStatus.FAILED_VALIDATION);
				}

				// Return results for validation
				convertParameters(fitParams);

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

				return createResult(newFitResult, results);
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
					String msg = String.format("Doublet %d [%d,%d] %s (%s) = %s\n", slice, cc.fromRegionToGlobalX(cx),
							cc.fromRegionToGlobalY(cy), newFitResult.getStatus(), newFitResult.getStatusData(),
							Arrays.toString(peakParams));
					logger2.debug(msg);
				}

				return createResult(newFitResult, null);
				//return null;
			}
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

		final Estimate[] estimates;

		// dx and dy should be <=1 pixel when a candidate is being fit since we use bounds.
		// They can be larger if we drifted close to another candidate (e.g. during doublet fitting)
		// or if this is the result of fitting the current candidate (which is not bounded).		
		if (dx < -1 || dx > 1 || dy < -1 || dy > 1)
		{
			//if (dynamicMultiPathFitResult.candidateId != i)
			//	System.out.printf("Drift error: [%d,%d]  %d  %.1f %.1f\n", slice, dynamicMultiPathFitResult.candidateId,
			//			i, dx, dy);

			// Ignore this as it is not a good estimate
			if (d2 > 2)
				return;

			// Store as a non-local estimate
			estimates = this.estimates2;
		}
		else
		{
			// Store as a close estimate
			estimates = this.estimates;
		}

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
			if (fitConfig.getBias() != 0 && this.noise != 0)
			{
				// Initial guess using the noise (assuming all noise is from Poisson background).
				// EMCCD will have increase noise by a factor of sqrt(2)
				background = (float) (fitConfig.getBias() + this.noise / ((fitConfig.isEmCCD()) ? 1.414213562 : 1));

				background = (float) (fitConfig.getBias() +
						PeakResult.noiseToLocalBackground(noise, fitConfig.getGain(), fitConfig.isEmCCD()));
			}
			else
			{
				// Initial guess using image mean
				background = this.background;
			}

		}
		return background;
	}

	private void resetNeighbours()
	{
		clearGridCache();
		candidateNeighbourCount = 0;
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
		candidateNeighbourCount = 0;

		// Using the neighbour grid.
		// Note this will also include all higher intensity spots that failed to be fit. 
		// These may still have estimates.
		final CandidateList neighbours = findNeighbours(spot);
		final Candidate[] candidates = neighbours.list;
		for (int i = 0; i < neighbours.size; i++)
		{
			final Candidate neighbour = candidates[i];
			if (isFit(neighbour.index))
				continue;
			if (canIgnore(neighbour.x, neighbour.y, xmin, xmax, ymin, ymax, neighbour.intensity, heightThreshold))
				continue;
			candidateNeighbours[candidateNeighbourCount++] = neighbour;
		}

		// XXX Debugging
		//		int c = 0;
		//		// Processing all lower spots.
		//		//for (int i = n + 1; i < candidates.length; i++)
		//		// Processing all spots.
		//		for (int i = 0; i < this.candidates.length; i++)
		//		{
		//			if (i == n || isFit(i))
		//				continue;
		//			if (canIgnore(this.candidates[i].x, this.candidates[i].y, xmin, xmax, ymin, ymax,
		//					this.candidates[i].intensity, heightThreshold))
		//				continue;
		//			//neighbourIndices[c++] = i;
		//			if (neighbourIndices[c++] != i)
		//				throw new RuntimeException("invalid grid neighbours");
		//		}

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
				final PeakResult neighbour = peakNeighbours[i];
				// No height threshold check as this is a validated peak
				final double xw = 2 * neighbour.getXSD();
				final double yw = 2 * neighbour.getYSD();
				if (intersects(x0min, y0min, x0max, y0max, neighbour.getXPosition() - xw, neighbour.getYPosition() - yw,
						neighbour.getXPosition() + xw, neighbour.getYPosition() + yw))
					fittedNeighbours[fittedNeighbourCount++] = neighbour;
			}
		}

		return candidateNeighbourCount + fittedNeighbourCount;
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
		createDataEstimator(width, height);
		return dataEstimator.getBackground();
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

	private float estimateNoise(int width, int height)
	{
		createDataEstimator(width, height);

		if (dataEstimator.isBackgroundRegion())
			return dataEstimator.getNoise();
		return dataEstimator.getNoise(config.getNoiseMethod());
	}

	/**
	 * Estimate the noise in the data.
	 *
	 * @param data
	 *            the data
	 * @param width
	 *            the width
	 * @param height
	 *            the height
	 * @param method
	 *            the method
	 * @return The noise
	 */
	public static float estimateNoise(float[] data, int width, int height, NoiseEstimator.Method method)
	{
		// Do the same logic as the non-static method 
		DataEstimator dataEstimator = newDataEstimator(data, width, height);

		if (dataEstimator.isBackgroundRegion())
			return dataEstimator.getNoise();
		return dataEstimator.getNoise(method);
	}

	private void createDataEstimator(int width, int height)
	{
		if (dataEstimator == null)
			dataEstimator = newDataEstimator(data, width, height);
	}

	private static DataEstimator newDataEstimator(float[] data, int width, int height)
	{
		// TODO - add options to control the thresholding method and the background fraction

		return new DataEstimator(data, width, height);
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
		private static final double NO_QA_SCORE = -1;

		final ImageExtractor ie;
		boolean dynamic;
		Rectangle regionBounds;
		double[] region;
		SpotFitter spotFitter;
		FitType fitType;

		public DynamicMultiPathFitResult(ImageExtractor ie, boolean dynamic)
		{
			this.frame = FitWorker.this.slice;
			this.width = cc.dataBounds.width;
			this.height = cc.dataBounds.height;
			this.ie = ie;
			this.dynamic = dynamic;
		}

		public void reset(int candidateId)
		{
			this.candidateId = candidateId;
			fitType = new FitType();

			// Reset results
			this.setMultiQAScore(NO_QA_SCORE);
			this.setSingleQAScore(NO_QA_SCORE);
			this.setMultiFitResult(null);
			this.setMultiDoubletFitResult(null);
			this.setSingleFitResult(null);
			this.setDoubletFitResult(null);

			// Set fitting region
			regionBounds = ie.getBoxRegionBounds(candidates[candidateId].x, candidates[candidateId].y, fitting);
			region = ie.crop(regionBounds, region);

			cc.setRegionBounds(regionBounds);

			// Offsets to convert fit coordinates to the global reference frame
			final float offsetx = cc.dataBounds.x + regionBounds.x + 0.5f;
			final float offsety = cc.dataBounds.y + regionBounds.y + 0.5f;

			// Note that the PreprocessedPeakResult will have coordinates 
			// in the global reference frame. We store estimates relative to 
			// the data bounds without the pixel offset making them suitable 
			// for initialising fitting.
			estimateOffsetx = -cc.dataBounds.x - 0.5f;
			estimateOffsety = -cc.dataBounds.y - 0.5f;

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
				if (result != null)
					fitType.setMulti(true);
			}
			return result;
		}

		public FitResult getSuperMultiFitResult()
		{
			// Pass through the reference to the result
			return super.getMultiFitResult();
		}

		@Override
		public double getMultiQAScore()
		{
			double score = super.getMultiQAScore();
			if (score == NO_QA_SCORE)
			{
				score = spotFitter.getMultiQAScore();
				this.setMultiQAScore(score);
			}
			return score;
		}

		@Override
		public FitResult getMultiDoubletFitResult()
		{
			FitResult result = super.getMultiDoubletFitResult();
			if (result == null)
			{
				result = spotFitter.getResultMultiDoublet(config.getResidualsThreshold());
				setMultiDoubletFitResult(result);
				fitType.setMultiDoublet(spotFitter.computedMultiDoublet);
			}
			return result;
		}

		public FitResult getSuperMultiDoubletFitResult()
		{
			// Pass through the reference to the result
			return super.getMultiDoubletFitResult();
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
				fitType.setDoublet(spotFitter.computedDoublet);
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
	 * @see gdsc.smlm.results.filter.IMultiPathFitResults#complete(int)
	 */
	public void complete(int index)
	{
		if (benchmarking)
		{
			// When benchmarking we must generate all the results possible
			// and store them in the job.
			// We do not assess the results and we do not store estimates.
			// This means that fitting results for the candidates that are not
			// processed by the main routine may not be representative
			// of fitting using a higher fail count or different residuals/neighbour
			// thresholds.

			// Calling the spot fitter with zero residuals will force the result to be computed if possible			
			dynamicMultiPathFitResult.spotFitter.getResultMultiDoublet(0);
			dynamicMultiPathFitResult.spotFitter.getResultSingleDoublet(0);

			// Now update the result if they were previously null
			dynamicMultiPathFitResult.getMultiFitResult();
			dynamicMultiPathFitResult.getMultiQAScore();
			dynamicMultiPathFitResult.getMultiDoubletFitResult();
			dynamicMultiPathFitResult.getSingleFitResult();
			dynamicMultiPathFitResult.getSingleQAScore();
			dynamicMultiPathFitResult.getDoubletFitResult();

			job.setMultiPathFitResult(index, dynamicMultiPathFitResult.copy(false));
		}

		// Send the actual results to the neighbour grid
		flushToGrid();
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
		// TODO - Print the current state of the dynamicMultiPathFitResult to file.
		// This will allow debugging what is different between the benchmark fit and the PeakFit.
		// Output:
		// slice
		// candidate Id
		// Initial and final params for each fit result.
		// Details of the selected result.
		// Then try to figure out why the benchmark fit deviates from PeakFit.

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

		if (queueSize != 0)
			throw new RuntimeException("There are results queued already!");

		for (int i = 0; i < results.length; i++)
		{
			if (results[i].isExistingResult())
				continue;
			if (results[i].isNewResult())
			{
				final double[] p = results[i].toGaussian2DParameters();

				// Store slice results relative to the data frame (not the global bounds)
				// Convert back so that 0,0 is the top left of the data bounds
				p[Gaussian2DFunction.X_POSITION] -= cc.dataBounds.x;
				p[Gaussian2DFunction.Y_POSITION] -= cc.dataBounds.y;

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

				if (logger != null)
				{
					// Show the shift, signal and width spread
					PreprocessedPeakResult peak = results[i];
					logger.info("Fit OK %d (%.1f,%.1f) [%d]: Shift = %.3f,%.3f : SNR = %.2f : Width = %.2f,%.2f",
							peak.getCandidateId(), peak.getX(), peak.getY(), peak.getId(),
							Math.sqrt(peak.getXRelativeShift2()), Math.sqrt(peak.getYRelativeShift2()), peak.getSNR(),
							peak.getXSDFactor(), peak.getYSDFactor());
				}
			}
			else
			{
				// This is a candidate that passed validation. Store the estimate as passing the primary filter.

				// We now do this is the pass() method.
				//storeEstimate(results[i].getCandidateId(), results[i], FILTER_RANK_PRIMARY);
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
					fitType.setMultiOK(true);
				else if (dynamicMultiPathFitResult.getSuperMultiDoubletFitResult() == selectedResult.fitResult)
					fitType.setMultiDoubletOK(true);
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
					// We log good results in the loop above. 
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
			if (peakParams != null)
			{
				// Parameters are the raw values from fitting the region. Convert for logging.
				peakParams = Arrays.copyOf(peakParams, peakParams.length);
				int npeaks = peakParams.length / 6;
				for (int i = 0; i < npeaks; i++)
				{
					peakParams[i * 6 + Gaussian2DFunction.X_POSITION] += cc.fromFitRegionToGlobalX();
					peakParams[i * 6 + Gaussian2DFunction.Y_POSITION] += cc.fromFitRegionToGlobalY();
					peakParams[i * 6 + Gaussian2DFunction.ANGLE] *= 180.0 / Math.PI;
				}
			}
			final int x = candidates[candidateId].x;
			final int y = candidates[candidateId].y;
			logger2.debug("%d:%d [%d,%d] %s (%s) = %s\n", slice, candidateId, cc.fromDataToGlobalX(x),
					cc.fromDataToGlobalY(y), fitResult.getStatus(), fitResult.getStatusData(),
					Arrays.toString(peakParams));
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
		// If we have an estimate then this is a valid candidate for fitting.
		// Q. Should we attempt fitting is we have only passed the min filter?
		//return (estimates[candidateId] != null && estimates[candidateId].filterRank == FILTER_RANK_PRIMARY) ||
		//		(estimates2[candidateId] != null && estimates2[candidateId].filterRank == FILTER_RANK_PRIMARY);
		return isValid[candidateId];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.MultiPathFilter.SelectedResultStore#pass(gdsc.smlm.results.filter.
	 * PreprocessedPeakResult)
	 */
	public void pass(PreprocessedPeakResult result)
	{
		// Do not ignore these. They may be from a fit result that is eventually not selected so we cannot 
		// wait until the add(...) method is called with the selected result.
		storeEstimate(result.getCandidateId(), result, FILTER_RANK_PRIMARY);

		// We must implement the same logic as the default SimpleSelectedResultStore which visits every
		// candidate that has passed the main filter 
		isValid[result.getCandidateId()] = true;
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

	/**
	 * Gets the noise estimate for the last processed job.
	 *
	 * @return the noise estimate
	 */
	public float getNoise()
	{
		return noise;
	}
}
