package gdsc.smlm.engine;

import java.awt.Rectangle;
import java.io.File;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.concurrent.BlockingQueue;

import org.apache.commons.math3.util.FastMath;

import gdsc.core.logging.Logger;
import gdsc.core.utils.ImageExtractor;
import gdsc.core.utils.Maths;
import gdsc.core.utils.NoiseEstimator;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.core.utils.Statistics;
import gdsc.core.utils.TextUtils;
import gdsc.core.utils.TurboList;
import gdsc.smlm.data.config.CalibrationReader;
import gdsc.smlm.data.config.ConfigurationException;
import gdsc.smlm.data.config.FitProtos.NoiseEstimatorMethod;
import gdsc.smlm.data.config.FitProtosHelper;
import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.data.config.PSFProtos.PSF;
import gdsc.smlm.data.config.PSFProtos.PSFType;

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
import gdsc.smlm.filters.BlockAverageDataProcessor;
import gdsc.smlm.filters.MaximaSpotFilter;
import gdsc.smlm.filters.Spot;
import gdsc.smlm.fitting.FastGaussian2DFitter;
import gdsc.smlm.fitting.FitResult;
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.FunctionSolver;
import gdsc.smlm.fitting.FunctionSolverType;
import gdsc.smlm.fitting.Gaussian2DFitter;
import gdsc.smlm.fitting.LSEFunctionSolver;
import gdsc.smlm.fitting.MLEFunctionSolver;
import gdsc.smlm.fitting.WLSEFunctionSolver;
import gdsc.smlm.function.StandardValueProcedure;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.function.gaussian.GaussianOverlapAnalysis;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.model.camera.CameraModel;
import gdsc.smlm.results.AttributePeakResult;
import gdsc.smlm.results.ExtendedPeakResult;
import gdsc.smlm.results.Gaussian2DPeakResultHelper;
import gdsc.smlm.results.IdPeakResult;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResultHelper;
import gdsc.smlm.results.PeakResults;
import gdsc.smlm.results.count.FailCounter;
import gdsc.smlm.results.filter.BasePreprocessedPeakResult.ResultType;
import gdsc.smlm.results.filter.CoordinateStore;
import gdsc.smlm.results.filter.CoordinateStoreFactory;
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
	private Rectangle lastBounds = null;
	private int fitting = 1;

	// Used for fitting
	private FitEngineConfiguration config;
	private FitConfiguration fitConfig;

	private PeakResults results;
	private PSFType psfType;
	private BlockingQueue<FitJob> jobs;
	private Gaussian2DFitter gf;
	private final double xsd, ysd;

	// Used for fitting methods
	private TurboList<PeakResult> sliceResults;
	private boolean useFittedBackground = false;
	private Statistics fittedBackground;
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
	private float noise;
	private final boolean calculateNoise;
	private boolean estimateSignal;
	private CandidateGridManager gridManager = null;
	// Contains the index in the list of maxima for any neighbours
	private int candidateNeighbourCount = 0;
	private Candidate[] candidateNeighbours = null;
	// Contains the index in the list of fitted results for any neighbours 
	private int fittedNeighbourCount = 0;
	// The fitted neighbours use the same parameters and result 
	// as output from fitting the function. They should be converted to PeakResults
	// at the end of fitting. This allows using different representations of the PSF.
	private Candidate[] fittedNeighbours = null;
	private CoordinateStore coordinateStore;

	private volatile boolean finished = false;

	static int WORKER_ID = 0;
	int workerId;

	private static byte FILTER_RANK_MINIMAL = (byte) 0;
	private static byte FILTER_RANK_PRIMARY = (byte) 1;

	private final boolean isFitCameraCounts;
	private final CameraModel cameraModel;

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

	/**
	 * Allow recording the pass/fail events sent to the FailCounter from the MultiPathFilter
	 */
	private class RecordingFailCounter implements FailCounter
	{
		final boolean[] pass;
		final FailCounter failCounter;

		RecordingFailCounter(boolean[] pass, FailCounter failCounter)
		{
			this.pass = pass;
			this.failCounter = failCounter;
		}

		public String getDescription()
		{
			return failCounter.getDescription();
		}

		public void pass()
		{
			// We record that this candidate generated new fit results
			pass[dynamicMultiPathFitResult.candidateId] = true;
			failCounter.pass();
		}

		public void pass(int n)
		{
			throw new IllegalStateException("Cannot record multiple passes");
		}

		public void fail()
		{
			failCounter.fail();
		}

		public void fail(int n)
		{
			throw new IllegalStateException("Cannot record multiple fails");
		}

		public boolean isOK()
		{
			return failCounter.isOK();
		}

		public FailCounter newCounter()
		{
			throw new IllegalStateException("Cannot record to a new instance");
		}

		public void reset()
		{
			failCounter.reset();
		}
	}

	private Estimate[] estimates = new Estimate[0], estimates2 = null;
	private boolean[] isValid = null;

	private CandidateList candidates = null;
	private CandidateList allNeighbours = null;
	private CandidateList allFittedNeighbours = null;

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
	 * @throws ConfigurationException
	 *             if the configuration is invalid
	 */
	public FitWorker(FitEngineConfiguration config, PeakResults results, BlockingQueue<FitJob> jobs)
			throws ConfigurationException
	{
		this.config = config;
		this.fitConfig = config.getFitConfiguration();

		// The fitting method is current tied to a Gaussian 2D function
		PSF psf = fitConfig.getPSF();
		if (!PSFHelper.isGaussian2D(psf))
			throw new ConfigurationException("Gaussian 2D PSF required");
		psfType = psf.getPsfType();

		this.results = results;
		this.jobs = jobs;
		this.logger = fitConfig.getLog();
		gf = new FastGaussian2DFitter(fitConfig);
		// Cache for convenience
		xsd = fitConfig.getInitialXSD();
		ysd = fitConfig.getInitialYSD();

		// Used for duplicate checking
		coordinateStore = CoordinateStoreFactory.create(0, 0, 0, 0,
				config.convertUsingHWHMax(config.getDuplicateDistanceParameter()));
		calculateNoise = fitConfig.getNoise() <= 0;
		if (!calculateNoise)
		{
			noise = (float) fitConfig.getNoise();
		}

		// Disable the use of the direct filter within the FitConfiguration validate method. 
		// This allows validate() to be used for basic filtering of all fit results (e.g. using the fit region bounds).
		// The validation of each result will be performed by the FitConfiguration implementation
		// of the IDirectFilter interface. This may involve the DirectFilter object.
		fitConfig.setSmartFilter(false);

		workerId = WORKER_ID++;

		// Store this flag so we know how to process the data
		isFitCameraCounts = fitConfig.isFitCameraCounts();
		cameraModel = fitConfig.getCameraModel();
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
		lastBounds = null;
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
		dataEstimator = null; // This is tied to the input data 

		// 06-Jun-2017
		// The data model was changed to store the signal in photons.
		// This allows support for per-pixel bias and gain (sCMOS cameras).

		// Remove the bias and gain. This is done for all solvers except:
		// - the legacy MLE solvers which model camera amplification
		// - the basic LVM solver without a camera calibration

		// Note: Assume that the camera model has been correctly initialised to be 
		// relative to the global origin.
		if (isFitCameraCounts)
		{
			cameraModel.removeBias(cc.dataBounds, data);
		}
		else
		{
			cameraModel.removeBiasAndGain(cc.dataBounds, data);
		}

		FitParameters params = job.getFitParameters();
		this.endT = (params != null) ? params.endT : -1;

		candidates = indentifySpots(job, width, height, params);

		if (candidates.getSize() == 0)
		{
			finish(job, start);
			return;
		}

		fittedBackground = new Statistics();

		// TODO - Better estimate of the background and the noise. Using all the image pixels
		// results in an estimate that is too high when there are many spots in the image.
		// Create a method that thresholds the image and finds the mean/sd of the thresholded image.

		// Note: Other code calls the static estimateNoise method.
		// So add a private instance method to estimate the noise and background using a static helper
		// class. This can also be called from the static estimateNoise method.

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

		//System.out.printf("Slice %d : Noise = %g\n", slice, noise);
		if (logger != null)
			logger.info("Slice %d: Noise = %f", slice, noise);

		final ImageExtractor ie = new ImageExtractor(data, width, height);
		double[] region = null;

		final float offsetx = cc.dataBounds.x;
		final float offsety = cc.dataBounds.y;

		if (params != null && params.fitTask == FitTask.MAXIMA_IDENITIFICATION)
		{
			final float sd0 = (float) xsd;
			final float sd1 = (float) ysd;
			for (int n = 0; n < candidates.getSize(); n++)
			{
				// Find the background using the perimeter of the data.
				// TODO - Perhaps the Gaussian Fitter should be used to produce the initial estimates but no actual fit done.
				// This would produce coords using the centre-of-mass.
				final Candidate candidate = candidates.get(n);
				int x = candidate.x;
				int y = candidate.y;
				final Rectangle regionBounds = ie.getBoxRegionBounds(x, y, fitting);
				region = ie.crop(regionBounds, region);
				final float b = (float) FastGaussian2DFitter.getBackground(region, regionBounds.width,
						regionBounds.height, 1);

				// Offset the coords to the centre of the pixel. Note the bounds will be added later.
				// Subtract the background to get the amplitude estimate then convert to signal.
				final float amplitude = candidate.intensity - ((relativeIntensity) ? 0 : b);
				final float signal = (float) (amplitude * 2.0 * Math.PI * sd0 * sd1);
				final int index = y * width + x;

				x += offsetx;
				y += offsety;
				final float[] peakParams = new float[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
				peakParams[Gaussian2DFunction.BACKGROUND] = b;
				peakParams[Gaussian2DFunction.SIGNAL] = signal;
				peakParams[Gaussian2DFunction.X_POSITION] = x + 0.5f;
				peakParams[Gaussian2DFunction.Y_POSITION] = y + 0.5f;
				//peakParams[Gaussian2DFunction.Z_POSITION] = 0;
				peakParams[Gaussian2DFunction.X_SD] = sd0;
				peakParams[Gaussian2DFunction.Y_SD] = sd1;
				//peakParams[Gaussian2DFunction.ANGLE] = 0;
				sliceResults.add(createResult(x, y, data[index], 0, noise, peakParams, null, n, 0));
			}
		}
		else
		{
			initialiseFitting();

			// Smooth the data to provide initial background estimates
			final BlockAverageDataProcessor processor = new BlockAverageDataProcessor(1, 1);
			final float[] smoothedData = processor.process(data, width, height);
			final ImageExtractor ie2 = new ImageExtractor(smoothedData, width, height);

			// Perform the Gaussian fit

			// The SpotFitter is used to create a dynamic MultiPathFitResult object.
			// This is then passed to a multi-path filter. Thus the same fitting decision process 
			// is used when benchmarking and when running on actual data.

			// Note: The SpotFitter labels each PreprocessedFitResult using the offset in the FitResult object.
			// The initial params and deviations can then be extracted for the results that pass the filter.

			MultiPathFilter filter;
			IMultiPathFitResults multiPathResults = this;
			SelectedResultStore store = this;
			coordinateStore = coordinateStore.resize(cc.dataBounds.x, cc.dataBounds.y, width, height);

			// TODO - Test if duplicate distance is now obsolete ...

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
			dynamicMultiPathFitResult = new DynamicMultiPathFitResult(ie, ie2, !benchmarking);
			//dynamicMultiPathFitResult = new DynamicMultiPathFitResult(ie, false);

			// Debug where the fit config may be different between benchmarking and fitting
			if (slice == -1)
			{
				fitConfig.initialise(1, 1, 1, null);

				SettingsManager.toJSON(config.getFitEngineSettings(),
						new File(String.format("/tmp/config.%d.txt", slice)), SettingsManager.FLAG_JSON_WHITESPACE);
				TextUtils.write(String.format("/tmp/filter.%d.xml", slice), filter.toXML());
				//filter.setDebugFile(String.format("/tmp/fitWorker.%b.txt", benchmarking));
				StringBuilder sb = new StringBuilder();
				sb.append((benchmarking) ? ((gdsc.smlm.results.filter.Filter) filter.getFilter()).toXML()
						: fitConfig.getSmartFilterString()).append("\n");
				sb.append(((gdsc.smlm.results.filter.Filter) filter.getMinimalFilter()).toXML()).append("\n");
				sb.append(filter.residualsThreshold).append("\n");
				sb.append(config.getFailuresLimit()).append("\n");
				sb.append(config.getDuplicateDistance()).append(":");
				sb.append(config.getDuplicateDistanceAbsolute()).append("\n");
				if (spotFilter != null)
					sb.append(spotFilter.getDescription()).append("\n");
				sb.append("MaxCandidate = ").append(candidates.getSize()).append("\n");
				for (int i = 0, len = candidates.getLength(); i < len; i++)
				{
					sb.append(String.format("Fit %d [%d,%d = %.1f]\n", i, candidates.get(i).x, candidates.get(i).y,
							candidates.get(i).intensity));
				}
				TextUtils.write(String.format("/tmp/candidates.%d.xml", slice), sb.toString());
			}

			FailCounter failCounter = config.getFailCounter();
			if (!benchmarking && params != null && params.pass != null)
			{
				// We want to store the pass/fail for consecutive candidates
				params.pass = new boolean[candidates.getLength()];
				failCounter = new RecordingFailCounter(params.pass, failCounter);
				filter.select(multiPathResults, failCounter, true, store, coordinateStore);
			}
			else
			{
				filter.select(multiPathResults, failCounter, true, store, coordinateStore);
			}

			// Note: We go deeper into the candidate list than max candidate
			// for any candidate where we have a good fit result as an estimate.
			// Q. Should this only be for benchmarking?

			//if (benchmarking)
			//	System.out.printf("Slice %d: %d + %d\n", slice, dynamicMultiPathFitResult.extra, candidates.getSize());

			// Create the slice results
			CandidateList fitted = gridManager.getFittedCandidates();
			sliceResults.ensureCapacity(fitted.getSize());
			for (int i = 0; i < fitted.getSize(); i++)
			{
				if (fitted.get(i).fit)
					sliceResults.addf(createResult(offsetx, offsety, fitted.get(i)));
			}

			if (logger != null)
				logger.info("Slice %d: %d / %d = %s", slice, success, candidates.getSize(),
						TextUtils.pleural(fitted.getSize(), "result"));
		}

		this.results.addAll(sliceResults);

		finish(job, start);
	}

	private CandidateList indentifySpots(FitJob job, int width, int height, FitParameters params)
	{
		Spot[] spots = null;
		int maxCandidate = 0;
		int[] maxIndices = null;

		// Only sub-classes may require the indices
		boolean requireIndices = (job.getClass() != FitJob.class);

		if (params != null)
		{
			maxCandidate = params.maxCandidate;
			if (params.spots != null)
			{
				spots = params.spots;
				if (maxCandidate <= 0 || maxCandidate > spots.length)
					maxCandidate = spots.length;
				// Get the indices for all candidates, even above the max candidate
				//maxIndices = new int[maxCandidate];
				//for (int n = 0; n < maxCandidate; n++)
				//{
				//	maxIndices[n] = spots[n].y * width + spots[n].x;
				//}
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
				if (maxCandidate <= 0 || maxCandidate > maxIndices.length)
					maxCandidate = maxIndices.length;
				else
					maxIndices = Arrays.copyOf(maxIndices, maxCandidate);
				float[] data2 = initialiseSpotFilter().preprocessData(data, width, height);
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
			spots = initialiseSpotFilter().rank(data, width, height);
			maxCandidate = spots.length;
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
			logger.info("%d: Slice %d: %d candidates", workerId, slice, maxCandidate);

		sliceResults = new TurboList<PeakResult>(maxCandidate);
		if (requireIndices)
		{
			job.setResults(sliceResults);
			job.setIndices(maxIndices);
		}

		final Candidate[] list = new Candidate[spots.length];
		for (int i = 0; i < spots.length; i++)
			list[i] = new Candidate(spots[i], i);
		return new CandidateList(maxCandidate, list);
	}

	private MaximaSpotFilter initialiseSpotFilter()
	{
		if (cameraModel.isPerPixelModel() && spotFilter.isWeighted())
		{
			// Use a per-pixel variance for weighting. 
			// Only get this if the bounds have changed to enable efficient caching.
			if (!cc.dataBounds.equals(lastBounds))
			{
				lastBounds = cc.dataBounds;
				float[] w = cameraModel.getWeights(lastBounds);
				spotFilter.setWeights(w, lastBounds.width, lastBounds.height);
			}
		}
		return spotFilter;
	}

	private void finish(FitJob job, final long start)
	{
		time += System.nanoTime() - start;
		job.finished();
	}

	private void initialiseFitting()
	{
		// Note that a ParameterisedFitJob can pass in a maxCandidate and the list is longer 
		// than candidates.getSize(). We must allocate for the maximum candidate Id (even if not processed).
		final int length = candidates.getLength();

		candidateNeighbours = allocateArray(candidateNeighbours, length);
		// Allocate enough room for all fits to be doublets 
		fittedNeighbours = allocateArray(fittedNeighbours, length * 2);

		clearEstimates(length);

		final int width = cc.dataBounds.width;
		final int height = cc.dataBounds.height;
		gridManager = new CandidateGridManager(width, height, 2 * fitting + 1);
		for (int i = 0; i < length; i++)
		{
			gridManager.putCandidateOnGrid(candidates.get(i));
		}
	}

	// Used to add fitted results to the grid for the current fit position.
	// This prevents filtering duplicates within the current fit results, 
	// only with all that has been fit before.
	private Candidate[] queue = new Candidate[5];
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
	private void queueToGrid(Candidate result)
	{
		if (queueSize == queue.length)
			queue = Arrays.copyOf(queue, queueSize * 2);
		queue[queueSize++] = result;
	}

	/**
	 * Flush the results for the current position to the grid.
	 *
	 * @return true, if successful
	 */
	private boolean flushToGrid()
	{
		if (queueSize == 0)
			return false;
		for (int i = 0; i < queueSize; i++)
			gridManager.putFittedOnGrid(queue[i]);
		queueSize = 0;
		return true;
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
	 * @param peakParamDevs
	 *            the peak params dev
	 * @param error
	 *            the error
	 * @param noise
	 *            the noise
	 * @param locationVariance
	 *            the location variance (in nm)
	 * @return true, if successful
	 */
	private boolean addSingleResult(int candidateId, float[] peakParams, float[] peakParamDevs, double error,
			float noise, double locationVariance)
	{
		Candidate c = candidates.get(candidateId);

		// Check if inside the allowed border
		boolean inside = insideBorder(peakParams[Gaussian2DFunction.X_POSITION],
				peakParams[Gaussian2DFunction.Y_POSITION]);

		// Add it to the grid of results (so we do not fit it again)
		int x = (int) peakParams[Gaussian2DFunction.X_POSITION];
		int y = (int) peakParams[Gaussian2DFunction.Y_POSITION];
		Candidate fitted = c.createFitted(x, y, candidateId, peakParams, peakParamDevs, error, noise, inside);
		if (locationVariance > 0)
			fitted.precision = Math.sqrt(locationVariance);
		queueToGrid(fitted);
		c.fit = true;

		// Check if the position is inside the border tolerance
		if (inside)
		{
			fittedBackground.add(peakParams[Gaussian2DFunction.BACKGROUND]);
		}
		else if (logger != null)
		{
			logger.info("[%d] Ignoring peak within image border @ %.2f,%.2f", slice,
					peakParams[Gaussian2DFunction.X_POSITION], peakParams[Gaussian2DFunction.Y_POSITION]);
		}
		return true;
	}

	private PeakResult createResult(float offsetx, float offsety, Candidate fitted)
	{
		int candidateId = fitted.index;
		Candidate c = candidates.get(candidateId);
		int x = c.x;
		int y = c.y;
		float value = data[y * cc.dataBounds.width + x];

		// Update to the global bounds.
		x += offsetx;
		y += offsety;

		float[] params = fitted.params;
		params[Gaussian2DFunction.X_POSITION] += offsetx;
		params[Gaussian2DFunction.Y_POSITION] += offsety;

		return createResult(x, y, value, fitted.error, fitted.noise, params, fitted.paramDevs, candidateId,
				c.precision);
	}

	private PeakResult createResult(int origX, int origY, float origValue, double error, float noise, float[] params,
			float[] paramDevs, int id, double precision)
	{
		// Convert to a variable PSF parameter PeakResult
		params = Gaussian2DPeakResultHelper.createParams(psfType, params);
		if (paramDevs != null)
		{
			paramDevs = Gaussian2DPeakResultHelper.createParams(psfType, paramDevs);
			// Convert variances to standard deviations
			for (int i = 0; i < paramDevs.length; i++)
				paramDevs[i] = (float) Math.sqrt(paramDevs[i]);
		}

		if (precision > 0)
		{
			AttributePeakResult r = new AttributePeakResult(slice, origX, origY, origValue, error, noise, params,
					paramDevs);
			r.setId(id);
			r.setPrecision(precision);
			if (endT >= 0 && slice != endT)
				r.setEndFrame(endT);
			return r;
		}

		if (endT >= 0 && slice != endT)
		{
			return new ExtendedPeakResult(slice, origX, origY, origValue, error, noise, params, paramDevs, endT, id);
		}
		return new IdPeakResult(slice, origX, origY, origValue, error, noise, params, paramDevs, id);
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
	@SuppressWarnings("unused")
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
				double[] params, double[] paramVariances, double localBackground, ResultType resultType)
		{
			//if (dynamicMultiPathFitResult.candidateId < candidateId && resultType == ResultType.NEW)
			//	System.out.println("WTF");

			// Update the initial params since we may have used an estimate
			// This will ensure that the width factor is computed correctly.

			// Q. Should this be ignored for existing results? They have already passed validation. 
			// So we do not have to be as strict on their width and could just use the drift from 
			// the the initial estimate.
			// For now do a full validation since multi-fit results are only accepted if existing
			// results are still valid.

			final int offset = n * Gaussian2DFunction.PARAMETERS_PER_PEAK;
			initialParams[Gaussian2DFunction.X_SD + offset] = xsd;
			initialParams[Gaussian2DFunction.Y_SD + offset] = ysd;
			return createResult(candidateId, n, initialParams, params, paramVariances, localBackground, resultType);
		}

		abstract PreprocessedPeakResult createResult(int candidateId, int n, double[] initialParams, double[] params,
				double[] paramVariances, double localBackground, ResultType resultType);
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
				double[] paramVariances, double localBackground, ResultType resultType)
		{
			return fitConfig.createDynamicPreprocessedPeakResult(candidateId, n, initialParams, params, paramVariances,
					localBackground, resultType, offsetx, offsety);
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
				double[] paramVariances,

				double localBackground, ResultType resultType)
		{
			return fitConfig.createPreprocessedPeakResult(slice, candidateId, n, initialParams, params, paramVariances,
					localBackground, resultType, offsetx, offsety);
		}
	}

	/**
	 * Provide functionality to fit spots in a region using different methods. Decisions about what to accept are not
	 * performed. The fit results are just converted to PreprocessedPeakResult objects for validation.
	 */
	private class CandidateSpotFitter
	{
		// TODO: When using an astigmatism z-model it is possible to fit two spots that are colocated.
		// So the colocation checks to existing peaks should be refined to allow both spots if they 
		// have suitably different z-depths (i.e. X/Y widths). 

		final Gaussian2DFitter gf;
		final ResultFactory resultFactory;
		final double[] region, region2; //, var_g2;
		final Rectangle regionBounds;
		final int candidateId;
		final int width;
		final int height;
		final int parametersPerPeak = Gaussian2DFunction.PARAMETERS_PER_PEAK;

		int neighbours;
		double singleBackground = Double.NaN, multiBackground = Double.NaN;

		// Flag each neighbour peak that is pre-computed for the multi-fit. 
		// These are fitted peaks outside the fit region. 
		private boolean[] precomputed = null;
		private int precomputedFittedNeighbourCount = -1;

		double[] precomputedFunctionParamsMulti = null;
		double[] precomputedFittedNeighboursMulti = null;
		MultiPathFitResult.FitResult resultMulti = null;
		boolean computedMulti = false;
		double[] residualsMulti = null;
		double valueMulti = 0;
		MultiPathFitResult.FitResult resultDoubletMulti = null;
		boolean computedDoubletMulti = false;
		QuadrantAnalysis qaMulti = null;

		double[] functionParamsSingle = null;
		double[] precomputedFittedNeighboursSingle = null;
		MultiPathFitResult.FitResult resultSingle = null;
		double[] residualsSingle = null;
		double valueSingle = 0;
		double localBackgroundSingle = 0;
		MultiPathFitResult.FitResult resultDoubletSingle = null;
		boolean computedDoubletSingle = false;
		QuadrantAnalysis qaSingle = null;

		public CandidateSpotFitter(Gaussian2DFitter gf, ResultFactory resultFactory, double[] region, double[] region2,
				double[] var_g2, Rectangle regionBounds, int candidateId)
		{
			this.gf = gf;
			this.resultFactory = resultFactory;
			this.region = region;
			this.region2 = region2;
			//this.var_g2 = var_g2;
			this.regionBounds = regionBounds;
			this.candidateId = candidateId;

			// Initialise
			width = regionBounds.width;
			height = regionBounds.height;

			fitConfig.setFitRegion(width, height, 0.5);
			// The variance is always needed in each fit of the same data 
			fitConfig.setObservationWeights(var_g2);

			// Analyse neighbours and include them in the fit if they are within a set height of this peak.
			resetNeighbours();
			neighbours = findNeighbours(regionBounds, candidateId, (float) getFittingBackgroundSingle());

			if (benchmarking)
			{
				// When benchmarking we may compute additional results after the main filtering routine
				// has been run. In this case we must pre-compute some values with the current results.
			}
		}

		@SuppressWarnings("unused")
		private double getFittingBackgroundMulti()
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
						multiBackground += fittedNeighbours[i].params[Gaussian2DFunction.BACKGROUND];
					}

					multiBackground /= fittedNeighbourCount;
					multiBackground = limitBackground(multiBackground);
				}
				else
				{
					multiBackground = this.getFittingBackgroundSingle();
				}
			}
			return multiBackground;
		}

		private double getFittingBackgroundSingle()
		{
			if (Double.isNaN(singleBackground))
			{
				// Note:
				// Looking at various images simulated with low, medium and high density spots
				// it is clear that a global estimate is not appropriate.
				//singleBackground = limitBackground(FitWorker.this.getSingleFittingBackground());

				// Use the min in the data
				//singleBackground = getDefaultBackground(region, width, height);

				// Use the min in smoothed data. This avoids noise
				singleBackground = getDefaultBackground(region2, width, height);
			}
			return singleBackground;
		}

		private double getDefaultBackground(double[] region, int width, int height)
		{
			// Use the minimum in the data.
			// This is what is done in the fitter if the background is zero.
			return limitBackground(FastGaussian2DFitter.getBackground(region, width, height, 2));
		}

		private double limitBackground(double b)
		{
			// Ensure we do not get a negative background
			return (b < 0) ? 0 : b;
		}

		private double getMax(double[] region, int width, int height)
		{
			double max = region[0];
			for (int i = width * height; --i > 0;)
				if (max < region[i])
					max = region[i];
			return max;
		}

		private int getPrecomputedNeighbourCount()
		{
			if (precomputedFittedNeighbourCount == -1)
			{
				precomputedFittedNeighbourCount = 0;
				if (fittedNeighbourCount > 0)
				{
					precomputed = new boolean[fittedNeighbourCount];

					// The fitted result will be relative to (0,0) in the fit data and already 
					// have an offset applied so that 0.5 is the centre of a pixel. (Note: the 
					// parameters were created from a PreprocessedPeakResult generated by the 
					// PreprocessedPeakResults factory.)
					// We can test the coordinates exactly against the fit frame.
					final float xmin = regionBounds.x;
					final float xmax = xmin + regionBounds.width;
					final float ymin = regionBounds.y;
					final float ymax = ymin + regionBounds.height;

					for (int i = 0; i < fittedNeighbourCount; i++)
					{
						final float[] params = fittedNeighbours[i].params;
						final float x = params[Gaussian2DFunction.X_POSITION];
						final float y = params[Gaussian2DFunction.Y_POSITION];

						// Pre-compute peaks if they are outside the fit region
						if (x < xmin || x > xmax || y < ymin || y > ymax)
						{
							precomputed[i] = true;
							precomputedFittedNeighbourCount++;
						}
					}
				}
			}
			return precomputedFittedNeighbourCount;
		}

		private double[] getPrecomputedFittedNeighbours()
		{
			if (precomputedFittedNeighboursMulti == null && precomputedFittedNeighbourCount != 0)
			{
				// The fitted result will be relative to (0,0) and have a 0.5 pixel offset applied
				final double xOffset = regionBounds.x + 0.5;
				final double yOffset = regionBounds.y + 0.5;

				// Pre-compute the already fitted peaks. 

				precomputedFunctionParamsMulti = new double[1 + parametersPerPeak * precomputedFittedNeighbourCount];
				for (int i = 0, j = 0; i < fittedNeighbourCount; i++)
				{
					if (!precomputed[i])
						continue;
					getFittedParams(i, precomputedFunctionParamsMulti, j);
					// Adjust position relative to extracted region
					precomputedFunctionParamsMulti[j + Gaussian2DFunction.X_POSITION] -= xOffset;
					precomputedFunctionParamsMulti[j + Gaussian2DFunction.Y_POSITION] -= yOffset;
					j += parametersPerPeak;
				}

				final Gaussian2DFunction func = fitConfig.createGaussianFunction(precomputedFittedNeighbourCount, width,
						height);
				precomputedFittedNeighboursMulti = new StandardValueProcedure().getValues(func,
						precomputedFunctionParamsMulti);

				//gdsc.core.ij.Utils.display("precomputedFunctionValues", precomputedFittedNeighboursMulti, width, height);
			}
			return precomputedFittedNeighboursMulti;
		}

		public MultiPathFitResult.FitResult getResultMulti()
		{
			if (computedMulti)
				return resultMulti;

			computedMulti = true;

			// Do not do a multi-fit if the configuration is not set to include neighbours 
			if (neighbours == 0 || !config.isIncludeNeighbours())
				return null;

			// -=-=-=-
			// TODO
			//
			// If we have fitted neighbours:
			// Precompute them and then try a single fit. It should work if
			// something is there. This can be used as an initial estimate for one of the
			// peaks in the multiple fit (i.e. the closest one)
			//
			// If the fits fails then we can guess that the region has no good peaks and 
			// it is not worth doing a multiple fit.			
			// -=-=-=-

			// Flag each peak that is precomputed
			getPrecomputedNeighbourCount();

			neighbours = candidateNeighbourCount + fittedNeighbourCount - precomputedFittedNeighbourCount;
			if (neighbours == 0)
			{
				// There are no neighbours after precomputation.
				// This will be the same result as the single result so return.
				return null;
			}

			// Background of fitted peaks within the region
			double background = 0;
			int backgroundCount = 0;

			for (int i = 0; i < fittedNeighbourCount; i++)
			{
				if (!precomputed[i])
				{
					// Used to estimate the background.
					// Q. Should this only use those within the region?
					background += fittedNeighbours[i].params[Gaussian2DFunction.BACKGROUND];
					backgroundCount++;
				}
			}

			// Create the parameters for the fit
			final int npeaks = 1 + neighbours;

			// Multiple-fit ...
			if (logger != null)
				logger.info("Slice %d: Multiple-fit (%d peaks : neighbours [%d + %d - %d])", slice, npeaks,
						candidateNeighbourCount, fittedNeighbourCount, precomputedFittedNeighbourCount);

			double[] params = new double[1 + npeaks * parametersPerPeak];

			// Estimate background.
			// Use the fitted peaks within the region or fall-back to the estimate for 
			// a single peak (i.e. low point of the region).
			params[Gaussian2DFunction.BACKGROUND] = (backgroundCount == 0) ? getFittingBackgroundSingle()
					: background / backgroundCount;
			// Store for debugging 
			background = params[Gaussian2DFunction.BACKGROUND];

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
			amplitudeEstimate[0] = getEstimate(candidates.get(candidateId), params, 0, true);

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
			// TODO - Test which is better: (1) precomputing fitted peaks; or (2) including them.
			// Initial tests show that Chi-squared is much lower when including them in the fit.

			if (fittedNeighbourCount > 0)
			{
				// The fitted result will be relative to (0,0) and have a 0.5 pixel offset applied
				final double xOffset = regionBounds.x + 0.5;
				final double yOffset = regionBounds.y + 0.5;

				getPrecomputedFittedNeighbours();

				// Add the details of the already fitted peaks
				for (int i = 0, j = (1 + candidateNeighbourCount) * parametersPerPeak; i < fittedNeighbourCount; i++)
				{
					if (precomputed[i])
						continue;
					getFittedParams(i, params, j);

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
			//double bbefore = params[Gaussian2DFunction.BACKGROUND];
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

				// Note: Amplitude estimates are amplitude above the background so we compare to zero
				if (minSignal <= 0)
				{
					// Reset background to the minimum value in the data.
					final double oldBackground = params[Gaussian2DFunction.BACKGROUND];
					params[Gaussian2DFunction.BACKGROUND] = getDefaultBackground(region, width, height);
					final double backgroundChange = oldBackground - params[Gaussian2DFunction.BACKGROUND];

					// Make the amplitude estimates higher by the change in background
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
			fitConfig.setPrecomputedFunctionValues(precomputedFittedNeighboursMulti);
			final FitResult fitResult = gf.fit(region, width, height, npeaks, params, amplitudeEstimate,
					params[Gaussian2DFunction.BACKGROUND] == 0);
			fitConfig.setPrecomputedFunctionValues(null);
			valueMulti = getFitValue();
			gf.setBounds(null, null);

			//			if (fitResult.getStatus() == FitStatus.BAD_PARAMETERS)
			//			{
			//				int x = candidates.get(candidateId).x;
			//				int y = candidates.get(candidateId).y;
			//				int index = (y-regionBounds.y) * width + (x-regionBounds.x);
			//				System.out.printf("Bad : [%d,%d] %d,%d %.1f (%.1f) B=%.1f (%.1f) : %s\n", slice, candidateId,
			//						x, y, candidates.get(candidateId).intensity, region[index],
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
			// Note: We ignore those parameters from peaks that were pre-computed in precomputedFittedNeighboursMulti.
			// These are outside the fit region and so should not usually overlap enough to effect the computation
			// of the fit deviations.
			final double[] fitParamStdDevs = fitResult.getParameterDeviations();

			initialParams[Gaussian2DFunction.X_POSITION] = candidates.get(candidateId).x - regionBounds.x;
			initialParams[Gaussian2DFunction.Y_POSITION] = candidates.get(candidateId).y - regionBounds.y;
			initialParams[Gaussian2DFunction.X_SD] = xsd;
			initialParams[Gaussian2DFunction.Y_SD] = ysd;

			// Perform validation of the candidate and existing peaks (other candidates are allowed to fail)
			if (fitResult.getStatus() == FitStatus.OK)
			{
				// The candidate peak
				if (fitConfig.validatePeak(0, initialParams, fitParams, fitParamStdDevs) != FitStatus.OK)
					return resultMulti = createResult(fitResult, null, fitConfig.getValidationResult());

				// Existing peaks
				for (int n = candidateNeighbourCount + 1; n < npeaks; n++)
				{
					if (fitConfig.validatePeak(n, initialParams, fitParams, fitParamStdDevs) != FitStatus.OK)
						return resultMulti = createResult(fitResult, null, fitConfig.getValidationResult());
				}
			}

			// Create the results
			PreprocessedPeakResult[] results = null;
			if (fitResult.getStatus() == FitStatus.OK)
			{
				residualsMulti = gf.getResiduals();

				//				// Debug background estimates
				//				double base = 1; //params[0] - fitConfig.getBias();
				//				System.out.printf("[%d] %d %.1f : %.1f  %.2f  %.1f  %.2f  %.1f  %.2f  %.1f  %.2f  %.1f  %.2f\n", slice,
				//						candidateId, params[0], background, (background - params[0]) / base,
				//						getMultiFittingBackground(), (getMultiFittingBackground() - params[0]) / base,
				//						getSingleFittingBackground(), (getSingleFittingBackground() - params[0]) / base,
				//						getDefaultBackground(this.region, width, height),
				//						(getDefaultBackground(this.region, width, height) - params[0]) / base,
				//						getDefaultBackground(region, width, height),
				//						(getDefaultBackground(region, width, height) - params[0]) / base);

				// Debug estimates verses fit.
				// Distinguish those we have estimated using the amplitudeEstimate array.
				// Trivial analysis shows that estimates have a lower relative error than a default initial guess.

				//				// Primary candidate
				//				int offset = 0;
				//				System.out.printf("[%d] %d MC %b %.2f %.2f %.2f %.2f %.2f %.2f\n", slice, candidateId, !amplitudeEstimate[0], 
				//						DoubleEquality.relativeError(initialParams[offset + 1], params[offset + 1]), 
				//						DoubleEquality.relativeError(initialParams[offset + 2], params[offset + 2]), 
				//						DoubleEquality.relativeError(initialParams[offset + 3], params[offset + 3]), 
				//						DoubleEquality.relativeError(initialParams[offset + 4], params[offset + 4]), 
				//						DoubleEquality.relativeError(initialParams[offset + 5], params[offset + 5]), 
				//						DoubleEquality.relativeError(initialParams[offset + 6], params[offset + 6])); 
				//				
				//				// Neighbours
				//				int nn = 1;
				//				for (int i = 0; i < candidateNeighbourCount; i++)
				//				{
				//					offset += Gaussian2DFunction.NUMBER_PER_PEAK;
				//					System.out.printf("[%d] %d N %b %.2f %.2f %.2f %.2f %.2f %.2f\n", slice, candidateId, !amplitudeEstimate[nn], 
				//							DoubleEquality.relativeError(initialParams[offset + 1], params[offset + 1]), 
				//							DoubleEquality.relativeError(initialParams[offset + 2], params[offset + 2]), 
				//							DoubleEquality.relativeError(initialParams[offset + 3], params[offset + 3]), 
				//							DoubleEquality.relativeError(initialParams[offset + 4], params[offset + 4]), 
				//							DoubleEquality.relativeError(initialParams[offset + 5], params[offset + 5]), 
				//							DoubleEquality.relativeError(initialParams[offset + 6], params[offset + 6])); 
				//					nn++;
				//				}
				//
				//				// Already fitted peaks
				//				for (int i = 0; i < fittedNeighbourCount; i++)
				//				{
				//					if (subtract[i])
				//						continue;
				//					offset += Gaussian2DFunction.NUMBER_PER_PEAK;
				//					System.out.printf("[%d] %d F %b %.2f %.2f %.2f %.2f %.2f %.2f\n", slice, candidateId, !amplitudeEstimate[nn], 
				//							DoubleEquality.relativeError(initialParams[offset + 1], params[offset + 1]), 
				//							DoubleEquality.relativeError(initialParams[offset + 2], params[offset + 2]), 
				//							DoubleEquality.relativeError(initialParams[offset + 3], params[offset + 3]), 
				//							DoubleEquality.relativeError(initialParams[offset + 4], params[offset + 4]), 
				//							DoubleEquality.relativeError(initialParams[offset + 5], params[offset + 5]), 
				//							DoubleEquality.relativeError(initialParams[offset + 6], params[offset + 6]));
				//					nn++;
				//				}

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

				final CandidateList peakNeighbours = findPeakNeighbours(candidates.get(candidateId));
				if (peakNeighbours.getSize() != 0)
				{
					// Coords for comparison to the real positions
					final float fcx2 = (float) (regionBounds.x + fitParams[Gaussian2DFunction.X_POSITION] + 0.5);
					final float fcy2 = (float) (regionBounds.y + fitParams[Gaussian2DFunction.Y_POSITION] + 0.5);
					float mind2 = (float) distanceToSingleFit2;
					int ii = -1;
					for (int i = 0; i < peakNeighbours.getSize(); i++)
					{
						final float d2 = distance2(fcx2, fcy2,
								peakNeighbours.get(i).params[Gaussian2DFunction.X_POSITION],
								peakNeighbours.get(i).params[Gaussian2DFunction.Y_POSITION]);
						if (mind2 > d2)
						{
							// There is another fitted result that is closer.
							// Note: The fit region is not centred on the other spot so this fit will probably
							// be worse and is discarded (not compared to the existing fit to get the best one). 

							mind2 = d2;
							ii = i;
							otherId = peakNeighbours.get(i).index;
						}
					}
					if (otherId != candidateId)
					{
						if (logger != null)
						{
							logger.info(
									"Bad peak: Fitted coordinates moved closer to another result (%d,%d : x=%.1f,y=%.1f : %.1f,%.1f)",
									candidates.get(candidateId).x, candidates.get(candidateId).y, fcx2, fcy2,
									peakNeighbours.get(ii).params[Gaussian2DFunction.X_POSITION],
									peakNeighbours.get(ii).params[Gaussian2DFunction.Y_POSITION]);
						}
						//System.out.printf("Multi drift to another result: [%d,%d] %d\n", slice, candidateId, otherId);
						resultType = ResultType.EXISTING;

						// Update the initial parameters to the position of the existing result so 
						// that drift is correct for filtering
						initialParams[Gaussian2DFunction.X_POSITION] = peakNeighbours
								.get(ii).params[Gaussian2DFunction.X_POSITION] - cc.fromFitRegionToGlobalX();
						initialParams[Gaussian2DFunction.Y_POSITION] = peakNeighbours
								.get(ii).params[Gaussian2DFunction.Y_POSITION] - cc.fromFitRegionToGlobalY();
					}
				}

				// 4. Check if there is an unfit candidate spot closer than the current candidate.
				//    This represents drift out to fit another unfitted spot.

				if (otherId != candidateId)
				{
					final CandidateList neighbours = findNeighbours(candidates.get(candidateId));
					if (neighbours.getSize() != 0)
					{
						// Position - do not add 0.5 pixel offset to allow distance comparison to integer candidate positions.
						float fcx2 = (float) (regionBounds.x + fitParams[Gaussian2DFunction.X_POSITION]);
						float fcy2 = (float) (regionBounds.y + fitParams[Gaussian2DFunction.Y_POSITION]);
						double mind2 = distanceToSingleFit2;
						for (int j = 0; j < neighbours.getSize(); j++)
						{
							final int id = neighbours.get(j).index;
							if (isFit(id))
								// This will be in the already fitted results instead so ignore...
								continue;
							final double d2 = distance2(fcx2, fcy2, candidates.get(id));
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
										candidates.get(candidateId).x, candidates.get(candidateId).y, fcx2 + 0.5f,
										fcy2 + 0.5f, candidates.get(otherId).x, candidates.get(otherId).y);
							}
							//System.out.printf("Multi drift to another candidate: [%d,%d] %d\n", slice, candidateId, otherId);

							// There is another candidate to be fit later that is closer.
							// This may be used as an estimate so we return it as such (i.e we do not ignore it)
							//otherId = candidateId;
							if (otherId > candidateId)
								resultType = ResultType.CANDIDATE;

							// Update the initial parameters to the position of the candidate so 
							// that drift is correct for filtering
							initialParams[Gaussian2DFunction.X_POSITION] = candidates.get(otherId).x - regionBounds.x;
							initialParams[Gaussian2DFunction.Y_POSITION] = candidates.get(otherId).y - regionBounds.y;
						}
					}
				}

				results = new PreprocessedPeakResult[npeaks];

				// We must compute a local background for all the spots
				final double[] frozenParams = fitParams.clone();
				final int flags = GaussianFunctionFactory.freeze(fitConfig.getFunctionFlags(),
						fitConfig.getAstigmatismZModel(), frozenParams);

				// Note: This could be the current candidate or drift to another candidate
				results[0] = resultFactory.createPreprocessedPeakResult(otherId, 0, initialParams, fitParams,
						fitParamStdDevs, getLocalBackground(0, npeaks, frozenParams, flags), resultType);

				// Neighbours
				int n = 1;
				for (int i = 0; i < candidateNeighbourCount; i++)
				{
					final Candidate candidateNeighbour = candidateNeighbours[i];
					results[n] = resultFactory.createPreprocessedPeakResult(candidateNeighbour.index, n, initialParams,
							fitParams, fitParamStdDevs, getLocalBackground(n, npeaks, frozenParams, flags),
							ResultType.CANDIDATE);
					n++;
				}

				// Already fitted peaks
				for (int i = 0; i < fittedNeighbourCount; i++)
				{
					if (precomputed[i])
						continue;
					int candidateId = fittedNeighbours[i].index;
					results[n] = resultFactory.createPreprocessedPeakResult(candidateId, n, initialParams, fitParams,
							fitParamStdDevs, getLocalBackground(n, npeaks, frozenParams, flags), ResultType.EXISTING);
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
			GaussianOverlapAnalysis overlap = new GaussianOverlapAnalysis(flags, null, extractSpotParams(params, n), 2);
			overlap.add(extractOtherParams(params, n, npeaks), true);
			double[] overlapData = overlap.getOverlapData();
			return overlapData[1] + params[Gaussian2DFunction.BACKGROUND];
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
				params[j + Gaussian2DFunction.Z_POSITION] = estimatedParams[Gaussian2DFunction.Z_POSITION];
				// Reset the width params if using an astigmatism z-model
				if (fitConfig.getAstigmatismZModel() != null)
				{
					params[j + Gaussian2DFunction.X_SD] = xsd;
					params[j + Gaussian2DFunction.Y_SD] = ysd;
				}
				else
				{
					params[j + Gaussian2DFunction.X_SD] = estimatedParams[Gaussian2DFunction.X_SD];
					params[j + Gaussian2DFunction.Y_SD] = estimatedParams[Gaussian2DFunction.Y_SD];
				}
				params[j + Gaussian2DFunction.ANGLE] = estimatedParams[Gaussian2DFunction.ANGLE];
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

		private void getFittedParams(int i, double[] params, int j)
		{
			float[] fittedParams = fittedNeighbours[i].params;
			params[j + Gaussian2DFunction.SIGNAL] = fittedParams[Gaussian2DFunction.SIGNAL];
			params[j + Gaussian2DFunction.X_POSITION] = fittedParams[Gaussian2DFunction.X_POSITION];
			params[j + Gaussian2DFunction.Y_POSITION] = fittedParams[Gaussian2DFunction.Y_POSITION];
			params[j + Gaussian2DFunction.Z_POSITION] = fittedParams[Gaussian2DFunction.Z_POSITION];
			// Reset the width params if using an astigmatism z-model
			if (fitConfig.getAstigmatismZModel() != null)
			{
				params[j + Gaussian2DFunction.X_SD] = xsd;
				params[j + Gaussian2DFunction.Y_SD] = ysd;
			}
			else
			{
				params[j + Gaussian2DFunction.X_SD] = fittedParams[Gaussian2DFunction.X_SD];
				params[j + Gaussian2DFunction.Y_SD] = fittedParams[Gaussian2DFunction.Y_SD];
			}
			params[j + Gaussian2DFunction.ANGLE] = fittedParams[Gaussian2DFunction.ANGLE];
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
			// The error is now set by the function solver. Not all function solvers can compute 
			// the sum-of-squares so we can no longer update the error to be the independent
			// of the solver.

			//final double r2 = 1 - (gf.getFinalResidualSumOfSquares() / gf.getTotalSumOfSquares());
			//result.setError(r2);
		}

		private double getQAScoreMulti()
		{
			if (qaMulti != null)
				return qaMulti.score;

			// Ensure we have a multi result
			getResultMulti();

			// Note this assumes that this method will be called after a multi fit and that the 
			// residuals were computed.
			final double[] residuals = residualsMulti;
			qaMulti = computeQA((residuals == null) ? null : (FitResult) resultMulti.data, regionBounds, residuals);

			return qaMulti.score;
		}

		public MultiPathFitResult.FitResult getResultDoubletMulti(double residualsThreshold)
		{
			// Fit a multi-fit. If all peaks are valid then precompute them, apart from the central candidate,
			// then fit as a doublet. Validate against an updated neighbourhood using the multi-fit results

			if (computedDoubletMulti)
				return resultDoubletMulti;

			if (residualsThreshold >= 1 || residualsThreshold < 0)
				return null;

			if (getQAScoreMulti() < residualsThreshold)
				return null;

			computedDoubletMulti = true;

			if (resultMulti.status != 0)
				return null;

			// Note this assumes that this method will be called after a multi fit and that the 
			// residuals were computed.
			if (residualsMulti == null)
				return null;

			// Ideally all multi-results must be valid to proceed with doublet fitting. 
			// We do not perform validation of the results. So we assume that the results have 
			// been checked and are valid and continue.

			final PreprocessedPeakResult[] fitResults = resultMulti.results;

			// Get the background for the multi-fit result
			final FitResult multiFitResult = (FitResult) resultMulti.data;
			final double[] fittedParams = multiFitResult.getParameters();
			final float background = (float) fittedParams[Gaussian2DFunction.BACKGROUND];

			// Get the neighbours 
			final CandidateList neighbours = findNeighbours(candidates.get(candidateId)).copy();

			// Exclude the fitted candidate neighbours from the candidate neighbours
			neighbours.removeIf(new CandidateList.Predicate()
			{
				public boolean test(Candidate candidate)
				{
					final int otherId = candidate.index;
					for (int j = 0; j < fitResults.length; j++)
					{
						if (fitResults[j].getCandidateId() == otherId)
						{
							return true;
						}
					}
					return false;
				}
			});

			// Get the fitted neighbours
			final CandidateList peakNeighbours = findPeakNeighbours(candidates.get(candidateId));
			final CandidateList peakNeighbours2 = peakNeighbours.copy();

			// Update with the fitted results from the multi fit
			NEXT_RESULT: for (int j = 0; j < fitResults.length; j++)
			{
				final int otherId = fitResults[j].getCandidateId();
				if (otherId == candidateId)
					// Ignore this as it is the current candidate
					continue;
				// Check if this is already a fitted neighbour
				for (int i = 0; i < peakNeighbours.getSize(); i++)
				{
					if (otherId == peakNeighbours.get(i).index)
					{
						// Options:
						// 1, Change the coordinates to those of the multi-fit
						// 2. Add the new multi-fit result and keep the old result (leaving 2 positions for the neighbour)
						// 3. Leave the position as that originally fitted.

						// Choose option 3 for simplicity. Note that the original fitted coordinates
						// should be more accurate as the fit region was centred around the spot. Also
						// note that due to bounding the fit coordinates (from the multi-fit) of any 
						// existing spot will be limited to a single pixel shift in XY.
						continue NEXT_RESULT;
					}
				}
				// Create a new dummy fitted neighbour
				// (Use similar logic to when we create the actual results in #add(SelectedResult))
				// Note that we do not unfreeze the parameters (i.e. the widths of from astigmatism z-model)
				// since we are only interested in the coordinates.
				final double[] p = fitResults[j].toGaussian2DParameters();
				final float[] params = new float[p.length];
				params[Gaussian2DFunction.BACKGROUND] = background;
				for (int i = 1; i < params.length; i++)
					params[i] = (float) p[i];
				// Store slice results relative to the data frame (not the global bounds)
				// Convert back so that 0,0 is the top left of the data bounds
				params[Gaussian2DFunction.X_POSITION] -= cc.dataBounds.x;
				params[Gaussian2DFunction.Y_POSITION] -= cc.dataBounds.y;
				int x = (int) params[Gaussian2DFunction.X_POSITION];
				int y = (int) params[Gaussian2DFunction.Y_POSITION];
				peakNeighbours2.add(new Candidate(x, y, otherId, params, null, 0, 0, true));
			}

			// Create the precomputed function values. This is the function defined by the 
			// fit result apart from the central candidate.
			int nPeaks = fittedParams.length / Gaussian2DFunction.PARAMETERS_PER_PEAK;
			double[] precomputedFunctionParams = null;
			double[] precomputedFunctionValues = null;
			if (nPeaks > 1)
			{
				// Remove the actual fitted peak. 
				precomputedFunctionParams = extractOtherParams(fittedParams, 0, nPeaks);
				Gaussian2DFunction func = fitConfig.createGaussianFunction(nPeaks - 1, width, height);
				precomputedFunctionValues = new StandardValueProcedure().getValues(func, precomputedFunctionParams);
			}

			// Add any precomputed values already computed. These are from peak outside the fit region.
			if (precomputedFittedNeighboursMulti != null)
			{
				if (precomputedFunctionValues == null)
				{
					precomputedFunctionValues = precomputedFittedNeighboursMulti;
				}
				else
				{
					for (int i = 0; i < precomputedFunctionValues.length; i++)
						precomputedFunctionValues[i] += precomputedFittedNeighboursMulti[i];
				}
			}

			// Build a dummy fit result object for fitting the single candidate to the data.
			// This will be a single evaluation of the function against the data.
			int nFittedParameters = fitConfig.createGaussianFunction(1, width, height).gradientIndices().length;
			int degreesOfFreedom = FastMath.max(region.length - nFittedParameters, 0);
			double[] parameters = Arrays.copyOf(fittedParams, 1 + Gaussian2DFunction.PARAMETERS_PER_PEAK);
			FitResult fitResult = new FitResult(FitStatus.OK, degreesOfFreedom, 0, null, parameters, null, 1,
					nFittedParameters, null, 0, 0);

			// Evaluate the multi fit as if fitted as a single peak. 
			// The resulting function value are used in the doublet fit
			double singleValue = valueMulti;

			//// XXX: These should be the same
			//{
			//	try
			//	{
			//		gf.setComputeResiduals(false);
			//		fitConfig.setPrecomputedFunctionValues(precomputedFunctionValues);
			//		if (!gf.evaluate(region, regionBounds.width, regionBounds.height, 1, parameters))
			//			return null;
			//	}
			//	finally
			//	{
			//		gf.setComputeResiduals(true);
			//		fitConfig.setPrecomputedFunctionValues(null);
			//	}
			//
			//	singleValue = getFitValue();
			//	if (singleValue != valueMulti)
			//		System.err.printf("Not same value: %f != %f\n", singleValue, valueMulti);
			//}

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

			resultDoubletMulti = fitAsDoublet(fitResult, region, precomputedFunctionValues, residualsThreshold,
					neighbours, peakNeighbours2, qaMulti, singleValue, 0);

			//			if (resultMultiDoublet != null && resultMultiDoublet.status == FitStatus.BAD_PARAMETERS.ordinal())
			//			{
			//				System.out.println("Bad params: " + Arrays.toString(parameters));
			//				//gdsc.core.ij.Utils.display("Region", region, width, height);
			//				//gdsc.core.ij.Utils.display("Residuals1", residuals, width, height);
			//			}

			if (resultDoubletMulti != null && resultDoubletMulti.status == 0 && resultDoubletMulti.results != null)
			{
				// The code below builds a combined result for the multi-fit and the primary candidate
				// fitted as a doublet. However this means that validation will be repeated on the spots that have 
				// been validated before. This is a small overhead but allows using the multi-doublet result to 
				// return all the fit results combined.

				// Fitting 2 spots is better than 1 and does not clash with any of the other multi-fit results.
				// Build a combined result with the doublet and the other multi-fit results
				FitResult doubletFitResult = (gdsc.smlm.fitting.FitResult) resultDoubletMulti.getData();

				double[] initialParams1 = doubletFitResult.getInitialParameters();
				double[] params1 = doubletFitResult.getParameters();

				double[] initialParams2 = multiFitResult.getInitialParameters();
				double[] params2 = multiFitResult.getParameters();

				// Create the initial parameters by adding an extra spot
				double[] initialParams = new double[initialParams2.length + Gaussian2DFunction.PARAMETERS_PER_PEAK];
				double[] params = new double[initialParams.length];

				final int srcPos = Gaussian2DFunction.getIndex(1, Gaussian2DFunction.SIGNAL);
				final int destPos = initialParams1.length;
				final int length = initialParams2.length - srcPos;

				System.arraycopy(initialParams1, 0, initialParams, 0, destPos);
				System.arraycopy(initialParams2, srcPos, initialParams, destPos, length);

				System.arraycopy(params1, 0, params, 0, destPos);
				System.arraycopy(params2, srcPos, params, destPos, length);

				final int npeaks = multiFitResult.getNumberOfPeaks() + 1;
				double[] paramDevs = null;
				if (multiFitResult.getParameterDeviations() != null)
				{
					// Recompute the deviations with all the parameters
					paramDevs = new double[params.length];
					// Add the pre-computed function from outside the region. The parameters for this
					// are ignored from the deviations computation.
					fitConfig.setPrecomputedFunctionValues(precomputedFittedNeighboursMulti);
					gf.computeDeviations(region, width, height, npeaks, params, paramDevs);
					fitConfig.setPrecomputedFunctionValues(null);
				}

				// Create all the output results
				PreprocessedPeakResult[] results = new PreprocessedPeakResult[resultDoubletMulti.results.length +
						resultMulti.results.length - 1];

				// We must compute a local background for all the spots
				final double[] frozenParams = params.clone();
				final int flags = GaussianFunctionFactory.freeze(fitConfig.getFunctionFlags(),
						fitConfig.getAstigmatismZModel(), frozenParams);

				int n = 0;
				for (int i = 0; i < resultDoubletMulti.results.length; i++)
				{
					PreprocessedPeakResult r = resultDoubletMulti.results[i];
					results[n] = resultFactory.createPreprocessedPeakResult(r.getCandidateId(), r.getId(),
							initialParams, params, paramDevs, getLocalBackground(n, npeaks, frozenParams, flags),
							(r.isExistingResult()) ? ResultType.EXISTING
									: (r.isNewResult()) ? ResultType.NEW : ResultType.CANDIDATE);
					n++;
				}
				// Ignore the first result (this was removed and fit as the doublet)
				for (int i = 1; i < resultMulti.results.length; i++)
				{
					PreprocessedPeakResult r = resultMulti.results[i];
					// Increment the ID by one since the position in the parameters array is moved to 
					// accommodate 2 preceding peaks and not 1 
					results[n] = resultFactory.createPreprocessedPeakResult(r.getCandidateId(), r.getId() + 1,
							initialParams, params, paramDevs, getLocalBackground(n, npeaks, frozenParams, flags),
							(r.isExistingResult()) ? ResultType.EXISTING
									: (r.isNewResult()) ? ResultType.NEW : ResultType.CANDIDATE);
					n++;
				}

				final int adjust = (fitConfig.isBackgroundFitting()) ? 1 : 0;
				nFittedParameters = npeaks *
						((multiFitResult.getNumberOfFittedParameters() - adjust) / multiFitResult.getNumberOfPeaks()) +
						adjust;
				//@formatter:off
				FitResult mdoubletFitResult =
						doubletFitResult.toBuilder()
						.setDegreesOfFreedom(FastMath.max(region.length - nFittedParameters, 0))
						.setError(doubletFitResult.getError())
						.setInitialParameters(initialParams)
						.setParameters(params)
						.setParameterDeviations(paramDevs)
						.setnPeaks(npeaks)
						.setnFittedParameters(nFittedParameters)
						.setIterations(doubletFitResult.getIterations())
						.setEvaluations(doubletFitResult.getEvaluations())
						.build();
				//@formatter:on

				resultDoubletMulti = createResult(mdoubletFitResult, results);
			}
			else
			{
				// Nothing to return. Do not set to null to allow reporting of the errors
				//resultMultiDoublet = null;
			}

			return resultDoubletMulti;
		}

		public MultiPathFitResult.FitResult getResultSingle()
		{
			if (resultSingle != null)
				return resultSingle;

			// Get each peak that is pre-computed.
			// This is done to compute equivalent deviations to the getResultMulti()
			// by using only those fitted peaks within the region.
			getPrecomputedNeighbourCount();

			// Background of fitted peaks within the region
			double background = 0;
			int backgroundCount = 0;

			// Subtract all fitted neighbours from the region
			if (fittedNeighbourCount != 0)
			{
				// The fitted result will be relative to (0,0) and have a 0.5 pixel offset applied
				final double xOffset = regionBounds.x + 0.5;
				final double yOffset = regionBounds.y + 0.5;

				//Utils.display("Region", region, width, height);

				double[] precomputedFunctionParamsSingle = new double[1 + parametersPerPeak * fittedNeighbourCount];
				for (int i = 0, j = 0; i < fittedNeighbourCount; i++, j += parametersPerPeak)
				{
					// Check if within the region
					if (!precomputed[i])
					{
						background += fittedNeighbours[i].params[Gaussian2DFunction.BACKGROUND];
						backgroundCount++;
					}

					getFittedParams(i, precomputedFunctionParamsSingle, j);

					// Adjust position relative to extracted region
					precomputedFunctionParamsSingle[j + Gaussian2DFunction.X_POSITION] -= xOffset;
					precomputedFunctionParamsSingle[j + Gaussian2DFunction.Y_POSITION] -= yOffset;
				}

				Gaussian2DFunction func = fitConfig.createGaussianFunction(fittedNeighbourCount, width, height);
				precomputedFittedNeighboursSingle = new StandardValueProcedure().getValues(func,
						precomputedFunctionParamsSingle);
			}

			final double[] params = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];

			// Estimate background.
			// Use the multi-fitting background as this uses the background of fitted neighbours if available.
			// Use the fitted peaks within the region or fall-back to the estimate for 
			// a single peak (i.e. low point of the region).
			params[Gaussian2DFunction.BACKGROUND] = (backgroundCount == 0) ? getFittingBackgroundSingle()
					: background / backgroundCount;
			// Store for debugging 
			background = params[Gaussian2DFunction.BACKGROUND];

			final boolean[] amplitudeEstimate = new boolean[1];

			// Re-use an estimate if we have it. Note that this may be quite far from the candidate.
			amplitudeEstimate[0] = getEstimate(candidates.get(candidateId), params, 0, false);
			// If we have no estimate the default will be an amplitude estimate.
			final boolean usingEstimate = !amplitudeEstimate[0];
			if (amplitudeEstimate[0])
			{
				// We can estimate the signal here instead of using the amplitude.
				// Do this when the fitting window covers enough of the Gaussian (e.g. 2.5xSD).
				float signal = 0;
				if (estimateSignal)
				{
					final double oldBackground = params[Gaussian2DFunction.BACKGROUND];
					double sum = 0;
					final int size = width * height;
					for (int i = size; i-- > 0;)
						sum += region[i];
					// Subtract any fitted peaks
					if (precomputedFittedNeighboursSingle != null)
						sum -= Maths.sum(precomputedFittedNeighboursSingle);
					signal = (float) (sum - oldBackground * size);
				}
				if (signal > 0)
				{
					amplitudeEstimate[0] = false;
					params[Gaussian2DFunction.SIGNAL] = signal;
				}
				else
				{
					// Resort to default amplitude estimate. Ensure this is above zero.
					if (params[Gaussian2DFunction.SIGNAL] <= 0)
					{
						// Reset to the single fitting background
						double oldBackground = params[Gaussian2DFunction.BACKGROUND];
						params[Gaussian2DFunction.BACKGROUND] = getFittingBackgroundSingle();
						double backgroundChange = oldBackground - params[Gaussian2DFunction.BACKGROUND];

						params[Gaussian2DFunction.SIGNAL] += backgroundChange;

						if (params[Gaussian2DFunction.SIGNAL] <= 0)
						{
							// Reset to the minimum value in the data.
							oldBackground = params[Gaussian2DFunction.BACKGROUND];
							params[Gaussian2DFunction.BACKGROUND] = getDefaultBackground(region, width, height);
							backgroundChange = oldBackground - params[Gaussian2DFunction.BACKGROUND];

							// Make the amplitude estimate higher by the change in background
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
			}

			// If there were neighbours or we have an estimate 
			// then use off grid pixels to prevent re-estimate of CoM
			// (since the CoM estimate will be skewed by the neighbours, or is not needed)
			if (fittedNeighbourCount != 0 || candidateNeighbourCount != 0 || usingEstimate)
			{
				if ((int) params[Gaussian2DFunction.X_POSITION] == params[Gaussian2DFunction.X_POSITION])
					params[Gaussian2DFunction.X_POSITION] += 0.001;
				if ((int) params[Gaussian2DFunction.Y_POSITION] == params[Gaussian2DFunction.Y_POSITION])
					params[Gaussian2DFunction.Y_POSITION] += 0.001;
			}

			fitConfig.setPrecomputedFunctionValues(precomputedFittedNeighboursSingle);
			final FitResult fitResult = gf.fit(region, width, height, 1, params, amplitudeEstimate,
					params[Gaussian2DFunction.BACKGROUND] == 0);
			fitConfig.setPrecomputedFunctionValues(null);
			valueSingle = getFitValue();
			updateError(fitResult);

			// Ensure the initial parameters are at the candidate position since we may have used an estimate.
			// This will ensure that drift is computed correctly.
			final double[] initialParams = fitResult.getInitialParameters();
			initialParams[Gaussian2DFunction.X_POSITION] = candidates.get(candidateId).x - regionBounds.x;
			initialParams[Gaussian2DFunction.Y_POSITION] = candidates.get(candidateId).y - regionBounds.y;

			// Create the results
			PreprocessedPeakResult[] results = null;
			if (fitResult.getStatus() == FitStatus.OK)
			{
				residualsSingle = gf.getResiduals();

				//				// Debug background estimates
				//				double base = params[0] - fitConfig.getBias();
				//				System.out.printf("[%d] %d %.1f : %.1f  %.2f  %.1f  %.2f  %.1f  %.2f  %.1f  %.2f  %.1f  %.2f\n", slice,
				//						candidateId, params[0], background, (background - params[0]) / base,
				//						getMultiFittingBackground(), (getMultiFittingBackground() - params[0]) / base,
				//						getSingleFittingBackground(), (getSingleFittingBackground() - params[0]) / base,
				//						getDefaultBackground(this.region, width, height),
				//						(getDefaultBackground(this.region, width, height) - params[0]) / base,
				//						getDefaultBackground(region, width, height),
				//						(getDefaultBackground(region, width, height) - params[0]) / base);

				// Debug estimates verses fit.
				// Distinguish those we have estimated using the amplitudeEstimate array.
				// Trivial analysis shows that estimates have a lower relative error than a default initial guess.
				//				int offset = 0;
				//				System.out.printf("[%d] %d C %b %.2f %.2f %.2f %.2f %.2f %.2f\n", slice, candidateId, !amplitudeEstimate[0], 
				//						DoubleEquality.relativeError(initialParams[offset + 1], params[offset + 1]), 
				//						DoubleEquality.relativeError(initialParams[offset + 2], params[offset + 2]), 
				//						DoubleEquality.relativeError(initialParams[offset + 3], params[offset + 3]), 
				//						DoubleEquality.relativeError(initialParams[offset + 4], params[offset + 4]), 
				//						DoubleEquality.relativeError(initialParams[offset + 5], params[offset + 5]), 
				//						DoubleEquality.relativeError(initialParams[offset + 6], params[offset + 6])); 

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

				final CandidateList peakNeighbours = findPeakNeighbours(candidates.get(candidateId));
				if (peakNeighbours.getSize() != 0)
				{
					// Coords for comparison to the real positions
					final float fcx2 = (float) (regionBounds.x + fitParams[Gaussian2DFunction.X_POSITION] + 0.5);
					final float fcy2 = (float) (regionBounds.y + fitParams[Gaussian2DFunction.Y_POSITION] + 0.5);
					float mind2 = (float) distanceToSingleFit2;
					int ii = -1;
					for (int i = 0; i < peakNeighbours.getSize(); i++)
					{
						final float d2 = distance2(fcx2, fcy2,
								peakNeighbours.get(i).params[Gaussian2DFunction.X_POSITION],
								peakNeighbours.get(i).params[Gaussian2DFunction.Y_POSITION]);
						if (mind2 > d2)
						{
							// There is another fitted result that is closer.
							// Note: The fit region is not centred on the other spot so this fit will probably
							// be worse and is discarded (not compared to the existing fit to get the best one). 

							mind2 = d2;
							ii = i;
							otherId = peakNeighbours.get(i).index;
						}
					}
					if (otherId != candidateId)
					{
						if (logger != null)
						{
							logger.info(
									"Bad peak: Fitted coordinates moved closer to another result (%d,%d : x=%.1f,y=%.1f : %.1f,%.1f)",
									candidates.get(candidateId).x, candidates.get(candidateId).y, fcx2, fcy2,
									peakNeighbours.get(ii).params[Gaussian2DFunction.X_POSITION],
									peakNeighbours.get(ii).params[Gaussian2DFunction.Y_POSITION]);
						}
						//System.out.printf("Single drift to another result: [%d,%d] %d\n", slice, candidateId, otherId);						
						resultType = ResultType.EXISTING;

						// Update the initial parameters to the position of the existing result so 
						// that drift is correct for filtering
						initialParams[Gaussian2DFunction.X_POSITION] = peakNeighbours
								.get(ii).params[Gaussian2DFunction.X_POSITION] - cc.fromFitRegionToGlobalX();
						initialParams[Gaussian2DFunction.Y_POSITION] = peakNeighbours
								.get(ii).params[Gaussian2DFunction.Y_POSITION] - cc.fromFitRegionToGlobalY();
					}
				}

				// 4. Check if there is an unfit candidate spot closer than the current candidate.
				//    This represents drift out to fit another unfitted spot.

				if (otherId != candidateId)
				{
					final CandidateList neighbours = findNeighbours(candidates.get(candidateId));
					if (neighbours.getSize() != 0)
					{
						// Position - do not add 0.5 pixel offset to allow distance comparison to integer candidate positions.
						float fcx2 = (float) (regionBounds.x + fitParams[Gaussian2DFunction.X_POSITION]);
						float fcy2 = (float) (regionBounds.y + fitParams[Gaussian2DFunction.Y_POSITION]);
						double mind2 = distanceToSingleFit2;
						for (int j = 0; j < neighbours.getSize(); j++)
						{
							final int id = neighbours.get(j).index;
							if (isFit(id))
								// This will be in the already fitted results instead so ignore...
								continue;
							final double d2 = distance2(fcx2, fcy2, candidates.get(id));
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
										candidates.get(candidateId).x, candidates.get(candidateId).y, fcx2 + 0.5f,
										fcy2 + 0.5f, candidates.get(otherId).x, candidates.get(otherId).y);
							}
							//System.out.printf("Single drift to another candidate: [%d,%d] %d\n", slice, candidateId, otherId);

							// There is another candidate to be fit later that is closer.
							// This may be used as an estimate so we return it as such (i.e we do not ignore it)
							//otherId = candidateId;
							if (otherId > candidateId)
								resultType = ResultType.CANDIDATE;

							// Update the initial parameters to the position of the candidate so 
							// that drift is correct for filtering
							initialParams[Gaussian2DFunction.X_POSITION] = candidates.get(otherId).x - regionBounds.x;
							initialParams[Gaussian2DFunction.Y_POSITION] = candidates.get(otherId).y - regionBounds.y;
						}
					}
				}

				results = new PreprocessedPeakResult[1];

				double[] fitParamDevs = fitResult.getParameterDeviations();
				localBackgroundSingle = 0;

				int npeaks = 1 + fittedNeighbourCount - precomputedFittedNeighbourCount;
				if (npeaks > 1)
				{
					// We must compute a local background using the influence from neighbours.
					// For equivalence with the multi-fit we only include the fits within the region.

					// The fitted result will be relative to (0,0) and have a 0.5 pixel offset applied
					final double xOffset = regionBounds.x + 0.5;
					final double yOffset = regionBounds.y + 0.5;

					functionParamsSingle = new double[1 + parametersPerPeak * npeaks];
					System.arraycopy(fitParams, 0, functionParamsSingle, 0, fitParams.length);
					for (int i = 0, j = parametersPerPeak; i < fittedNeighbourCount; i++)
					{
						// Check if within the region
						if (!precomputed[i])
						{
							getFittedParams(i, functionParamsSingle, j);

							// Adjust position relative to extracted region
							functionParamsSingle[j + Gaussian2DFunction.X_POSITION] -= xOffset;
							functionParamsSingle[j + Gaussian2DFunction.Y_POSITION] -= yOffset;
							j += parametersPerPeak;
						}
					}

					final double[] frozenParams = functionParamsSingle.clone();
					final int flags = GaussianFunctionFactory.freeze(fitConfig.getFunctionFlags(),
							fitConfig.getAstigmatismZModel(), frozenParams);
					localBackgroundSingle = getLocalBackground(0, npeaks, frozenParams, flags);

					if (fitParamDevs != null)
					{
						// Recompute the deviations with all the parameters.
						double[] paramDevs1 = new double[functionParamsSingle.length];
						// These pre-computed values will be those peaks outside the region
						fitConfig.setPrecomputedFunctionValues(getPrecomputedFittedNeighbours());
						if (gf.computeDeviations(region, width, height, npeaks, functionParamsSingle, paramDevs1))
							System.arraycopy(paramDevs1, 0, fitParamDevs, 0, fitParamDevs.length);
						fitConfig.setPrecomputedFunctionValues(null);
					}
				}
				results[0] = resultFactory.createPreprocessedPeakResult(otherId, 0, initialParams, fitParams,
						fitParamDevs, localBackgroundSingle, resultType);
			}

			resultSingle = createResult(fitResult, results);

			return resultSingle;
		}

		private double getFitValue()
		{
			FunctionSolver solver = gf.getFunctionSolver();
			if (solver.getType() == FunctionSolverType.MLE)
				return ((MLEFunctionSolver) solver).getLogLikelihood();
			else if (solver.getType() == FunctionSolverType.WLSE)
				return ((WLSEFunctionSolver) solver).getChiSquared();
			else
				return ((LSEFunctionSolver) solver).getAdjustedCoefficientOfDetermination();
		}

		private String getFitValueName()
		{
			FunctionSolver solver = gf.getFunctionSolver();
			if (solver.getType() == FunctionSolverType.MLE)
				return "Log-likelihood";
			else if (solver.getType() == FunctionSolverType.WLSE)
				return "Chi-Squared";
			else
				return "Adjusted R^2";
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

		private double getQAScoreSingle()
		{
			if (qaSingle != null)
				return qaSingle.score;

			// Ensure we have a single result
			getResultSingle();

			// Note this assumes that this method will be called after a single fit and that the 
			// residuals were computed.
			final double[] residuals = residualsSingle;
			qaSingle = computeQA((residuals == null) ? null : (FitResult) resultSingle.data, regionBounds, residuals);

			return qaSingle.score;
		}

		public MultiPathFitResult.FitResult getResultDoubletSingle(double residualsThreshold)
		{
			if (computedDoubletSingle)
				return resultDoubletSingle;

			if (residualsThreshold >= 1 || residualsThreshold < 0)
				return null;

			if (getQAScoreSingle() < residualsThreshold)
				return null;

			computedDoubletSingle = true;

			final CandidateList neighbours = findNeighbours(candidates.get(candidateId));
			final CandidateList peakNeighbours = findPeakNeighbours(candidates.get(candidateId));

			resultDoubletSingle = fitAsDoublet((FitResult) resultSingle.data, region, precomputedFittedNeighboursSingle,
					residualsThreshold, neighbours, peakNeighbours, qaSingle, valueSingle, localBackgroundSingle);

			// Check if the deviations require updating. 
			// This will be the case if we have stored function parameters for the single fit.
			if (resultDoubletSingle != null && resultDoubletSingle.status == 0 && resultDoubletSingle.results != null &&
					functionParamsSingle != null)
			{
				// Recompute the deviations with all the parameters.
				// For equivalence with the multi-fit we only include the fits within the region.
				FitResult doubletFitResult = (gdsc.smlm.fitting.FitResult) resultDoubletSingle.getData();
				double[] fitParams = doubletFitResult.getParameters();
				double[] fitParamDevs = new double[fitParams.length];

				int npeaks = 2 + fittedNeighbourCount - precomputedFittedNeighbourCount;
				double[] params = new double[1 + parametersPerPeak * npeaks];
				System.arraycopy(fitParams, 0, params, 0, fitParams.length);
				System.arraycopy(functionParamsSingle, 1 + parametersPerPeak, params, fitParams.length,
						params.length - fitParams.length);

				double[] paramDevs = new double[params.length];
				// These pre-computed values will be those peaks outside the region
				fitConfig.setPrecomputedFunctionValues(getPrecomputedFittedNeighbours());
				if (gf.computeDeviations(region, width, height, npeaks, params, paramDevs))
					System.arraycopy(paramDevs, 0, fitParamDevs, 0, fitParamDevs.length);
				fitConfig.setPrecomputedFunctionValues(null);

				// Use the updated deviations
				doubletFitResult = doubletFitResult.toBuilder().setParameterDeviations(fitParamDevs).build();
				resultDoubletSingle = createResult(doubletFitResult, resultDoubletSingle.results);

				double[] initialParams = doubletFitResult.getInitialParameters();
				PreprocessedPeakResult[] results = resultDoubletSingle.results;
				for (int i = 0; i < results.length; i++)
				{
					PreprocessedPeakResult r = results[i];
					results[i] = resultFactory.createPreprocessedPeakResult(r.getCandidateId(), r.getId(),
							initialParams, params, paramDevs, localBackgroundSingle, (r.isExistingResult())
									? ResultType.EXISTING : (r.isNewResult()) ? ResultType.NEW : ResultType.CANDIDATE);
				}
			}

			return resultDoubletSingle;
		}

		/**
		 * Perform quadrant analysis on the residuals.
		 * <p>
		 * Perform quadrant analysis as per rapidSTORM to analyse if the residuals of the the fit are skewed around the
		 * single fit centre. This may indicate the result is actually two spots (a doublet).
		 *
		 * @param fitResult
		 *            the fit result
		 * @param regionBounds
		 *            the region bounds
		 * @param residuals
		 *            the residuals
		 * @return the quadrant analysis
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
		 *
		 * @param fitResult
		 *            the fit result
		 * @param region
		 *            the region
		 * @param precomputedFunctionValues
		 *            the precomputed function values
		 * @param residualsThreshold
		 *            the residuals threshold
		 * @param neighbours
		 *            the neighbours
		 * @param peakNeighbours
		 *            the peak neighbours
		 * @param qa
		 *            the qa object that performed quadrant analysis
		 * @param singleValue
		 *            the objective function value from fitting a single peak
		 * @param localBackground
		 *            the local background
		 * @return the multi path fit result. fit result
		 */
		private MultiPathFitResult.FitResult fitAsDoublet(FitResult fitResult, double[] region,
				double[] precomputedFunctionValues, double residualsThreshold, CandidateList neighbours,
				CandidateList peakNeighbours, QuadrantAnalysis qa, double singleValue, double localBackground)
		{
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
			final double[] doubletParams = new double[1 + 2 * Gaussian2DFunction.PARAMETERS_PER_PEAK];

			// Note: Quadrant analysis sets the positions using 0.5,0.5 as the centre of the pixel.
			// The fitting does not, so subtract 0.5.
			// Note that we set position and signal but leave the z-position and widths to their 
			// standard values since this is 'new' fit.

			doubletParams[Gaussian2DFunction.BACKGROUND] = params[Gaussian2DFunction.BACKGROUND];
			doubletParams[Gaussian2DFunction.SIGNAL] = params[Gaussian2DFunction.SIGNAL] * 0.5;
			doubletParams[Gaussian2DFunction.X_POSITION] = (qa.x1 - 0.5);
			doubletParams[Gaussian2DFunction.Y_POSITION] = (qa.y1 - 0.5);
			doubletParams[Gaussian2DFunction.PARAMETERS_PER_PEAK +
					Gaussian2DFunction.SIGNAL] = params[Gaussian2DFunction.SIGNAL] * 0.5;
			doubletParams[Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_POSITION] = (qa.x2 - 0.5);
			doubletParams[Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.Y_POSITION] = (qa.y2 - 0.5);
			// -+-+-

			// If validation is on in the fit configuration:
			// - Disable checking of position movement since we guessed the 2-peak location.
			// - Disable checking within the fit region as we do that per peak 
			//   (which is better than failing if either peak is outside the region)
			// - Increase the iterations level then reset afterwards.
			// Note: If simple validation is using the precision from the fit deviations
			// then it will be a under-estimate. This is good enough for initial filtering
			// and the deviations for the entire function parameters will be re-computed 
			// for the pre-processed peak results.

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
			fitConfig.setPrecomputedFunctionValues(precomputedFunctionValues);
			final FitResult newFitResult = gf.fit(region, width, height, 2, doubletParams, amplitudeEstimate,
					doubletParams[Gaussian2DFunction.BACKGROUND] == 0);
			fitConfig.setPrecomputedFunctionValues(null);
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

				// TODO - Make the selection criteria for Doublets configurable:
				// MLE - AIC, BIC, LLR
				// WLSE - AIC, BIC, q-values of each chi-square
				// LSE - Adjusted coefficient of determination 

				// Note: Numerical recipes pp 669 uses 0.1 for q-value for weighted least squares fitting
				// This is 0.9 for p!

				double doubleValue = getFitValue();
				final int length = width * height;
				double ic1 = Double.NaN, ic2 = Double.NaN;
				final boolean improvement;

				FunctionSolver solver = gf.getFunctionSolver();
				if (solver.getType() == FunctionSolverType.MLE
				//&& fitConfig.isModelCamera()
				)
				{
					// ------------
					// TODO: Check this is still true as we may need to change the improvement criterion. 
					// Note: the residuals are no longer computed for all solvers.
					// The MLE is only good if we are modelling the camera noise. 
					// The MLE put out by the Poisson model is not better than using the IC from the fit residuals.
					// ------------
					ic1 = Maths.getBayesianInformationCriterion(singleValue, length,
							fitResult.getNumberOfFittedParameters());
					ic2 = Maths.getBayesianInformationCriterion(doubleValue, length,
							newFitResult.getNumberOfFittedParameters());
					// IC should be lower
					improvement = ic2 < ic1;
					if (logger != null)
						logger.info("Model improvement - Log likelihood (IC) : %f (%f) => %f (%f) : %f", singleValue,
								ic1, doubleValue, ic2, ic1 - ic2);
				}
				else if (solver.getType() == FunctionSolverType.WLSE)
				{
					// If using the weighted least squares estimator then we can get the log likelihood from an approximation
					ic1 = Maths.getBayesianInformationCriterionFromResiduals(singleValue, length,
							fitResult.getNumberOfFittedParameters());
					ic2 = Maths.getBayesianInformationCriterionFromResiduals(doubleValue, length,
							newFitResult.getNumberOfFittedParameters());
					// IC should be lower
					improvement = ic2 < ic1;
					if (logger != null)
						logger.info("Model improvement - Chi-squared (IC) : %f (%f) => %f (%f) : %f", singleValue, ic1,
								doubleValue, ic2, ic1 - ic2);
				}
				else if (solver.getType() == FunctionSolverType.LSE)
				{
					// Adjusted r^2 should be higher
					improvement = doubleValue > singleValue;
					if (logger != null)
						logger.info("Model improvement - Adjusted R^2 : %f => %f", singleValue, doubleValue);
				}
				else
				{
					throw new IllegalStateException("Unable to calculate solution improvement");
				}

				if (logger2 != null)
				{
					double[] peakParams = newFitResult.getParameters();
					if (peakParams != null)
					{
						peakParams = Arrays.copyOf(peakParams, peakParams.length);
						int npeaks = peakParams.length / Gaussian2DFunction.PARAMETERS_PER_PEAK;
						for (int i = 0; i < npeaks; i++)
						{
							peakParams[i * Gaussian2DFunction.PARAMETERS_PER_PEAK +
									Gaussian2DFunction.X_POSITION] += 0.5 + regionBounds.x;
							peakParams[i * Gaussian2DFunction.PARAMETERS_PER_PEAK +
									Gaussian2DFunction.Y_POSITION] += 0.5 + regionBounds.y;
						}
					}
					String msg = String.format("Doublet %d [%d,%d] %s (%s) %s [%f -> %f] IC [%f -> %f] = %s\n", slice,
							cc.fromRegionToGlobalX(cx), cc.fromRegionToGlobalY(cy), newFitResult.getStatus(),
							newFitResult.getStatusData(), getFitValueName(), singleValue, doubleValue, ic1, ic2,
							Arrays.toString(peakParams));
					logger2.debug(msg);
				}

				// Check if the predictive power of the model is better with two peaks:
				if (!improvement)
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
					final int offset = n * Gaussian2DFunction.PARAMETERS_PER_PEAK;

					// Ensure the initial parameters are at the candidate position since we may have used an estimate.
					// This will ensure that drift is computed correctly.
					initialParams[Gaussian2DFunction.X_POSITION + offset] = candidates.get(candidateId).x -
							regionBounds.x;
					initialParams[Gaussian2DFunction.Y_POSITION + offset] = candidates.get(candidateId).y -
							regionBounds.y;

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

					double xShift = fitParams[Gaussian2DFunction.X_POSITION + offset] -
							params[Gaussian2DFunction.X_POSITION];
					double yShift = fitParams[Gaussian2DFunction.Y_POSITION + offset] -
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

					// Distance to current candidate 
					xShift = fitParams[Gaussian2DFunction.X_POSITION + offset] -
							initialParams[Gaussian2DFunction.X_POSITION];
					yShift = fitParams[Gaussian2DFunction.Y_POSITION + offset] -
							initialParams[Gaussian2DFunction.Y_POSITION];
					final double distanceToSingleFit2 = xShift * xShift + yShift * yShift;

					// 3. Check we are not closer to a fitted spot. This has already had a chance at 
					//    fitting a doublet so is ignored..

					if (peakNeighbours.getSize() != 0)
					{
						// Coords for comparison to the real positions
						final float fcx2 = (float) (regionBounds.x + fitParams[Gaussian2DFunction.X_POSITION + offset] +
								0.5);
						final float fcy2 = (float) (regionBounds.y + fitParams[Gaussian2DFunction.Y_POSITION + offset] +
								0.5);
						final float d2 = (float) distanceToSingleFit2;
						for (int i = 0; i < peakNeighbours.getSize(); i++)
						{
							if (d2 > distance2(fcx2, fcy2, peakNeighbours.get(i).params[Gaussian2DFunction.X_POSITION],
									peakNeighbours.get(i).params[Gaussian2DFunction.Y_POSITION]))
							{
								if (logger != null)
								{
									logger.info(
											"Bad peak %d: Fitted coordinates moved closer to another result (%d,%d : x=%.1f,y=%.1f : %.1f,%.1f)",
											n, candidates.get(candidateId).x, candidates.get(candidateId).y, fcx2, fcy2,
											peakNeighbours.get(i).params[Gaussian2DFunction.X_POSITION],
											peakNeighbours.get(i).params[Gaussian2DFunction.Y_POSITION]);
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
					if (neighbours.getSize() != 0)
					{
						// Position - do not add 0.5 pixel offset to allow distance comparison to integer candidate positions.
						float fcx2 = (float) (regionBounds.x + fitParams[Gaussian2DFunction.X_POSITION + offset]);
						float fcy2 = (float) (regionBounds.y + fitParams[Gaussian2DFunction.Y_POSITION + offset]);
						double mind2 = distanceToSingleFit2;
						for (int j = 0; j < neighbours.getSize(); j++)
						{
							final int id = neighbours.get(j).index;
							if (isFit(id))
								// This will be in the already fitted results instead so ignore...
								continue;
							final double d2 = distance2(fcx2, fcy2, candidates.get(id));
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
										n, candidates.get(candidateId).x, candidates.get(candidateId).y, fcx2 + 0.5f,
										fcy2 + 0.5f, candidates.get(otherId).x, candidates.get(otherId).y);
							}

							// There is another candidate to be fit later that is closer.
							// This may be used as an estimate so we return it as such (i.e we do not ignore it)

							// Update the initial parameters to the position of the candidate so 
							// that drift is correct for filtering
							initialParams[Gaussian2DFunction.X_POSITION + offset] = candidates.get(otherId).x -
									regionBounds.x;
							initialParams[Gaussian2DFunction.Y_POSITION + offset] = candidates.get(otherId).y -
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
				final PreprocessedPeakResult[] results = new PreprocessedPeakResult[nPeaks];
				for (int i = 0; i < nPeaks; i++)
				{
					// If it is this candidate, or an earlier one that was not fit then this is a new result.
					// Otherwise it is a candidate we will process later
					final ResultType resultType = (candidateIndex[i] <= candidateId) ? ResultType.NEW
							: ResultType.CANDIDATE;
					double[] fitParamDevs = (precomputedFunctionValues == null)
							// For a single fit as a doublet we can use the deviations directly
							? newFitResult.getParameterDeviations()
							// If there was a pre-computed function then the deviations must be recomputed using the neighbours.
							: null;
					results[i] = resultFactory.createPreprocessedPeakResult(candidateIndex[i], position[i],
							initialParams, fitParams, fitParamDevs, localBackground, resultType);
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
						int npeaks = peakParams.length / Gaussian2DFunction.PARAMETERS_PER_PEAK;
						for (int i = 0; i < npeaks; i++)
						{
							peakParams[i * Gaussian2DFunction.PARAMETERS_PER_PEAK +
									Gaussian2DFunction.X_POSITION] += 0.5 + regionBounds.x;
							peakParams[i * Gaussian2DFunction.PARAMETERS_PER_PEAK +
									Gaussian2DFunction.Y_POSITION] += 0.5 + regionBounds.y;
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

		private float distance2(float cx, float cy, float x, float y)
		{
			final float dx = cx - x;
			final float dy = cy - y;
			return dx * dx + dy * dy;
		}
	}

	private double estimateOffsetx, estimateOffsety;

	private void storeEstimate(int i, PreprocessedPeakResult peak, byte filterRank)
	{
		final double[] params = peak.toGaussian2DParameters();
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
		final double dx = candidates.get(i).x - params[Gaussian2DFunction.X_POSITION];
		final double dy = candidates.get(i).y - params[Gaussian2DFunction.Y_POSITION];

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
		final double[] newParams = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
		System.arraycopy(params, n * Gaussian2DFunction.PARAMETERS_PER_PEAK + 1, newParams, 1,
				Gaussian2DFunction.PARAMETERS_PER_PEAK);
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
		final double[] newParams = new double[params.length - Gaussian2DFunction.PARAMETERS_PER_PEAK];
		if (n > 0)
		{
			System.arraycopy(params, 1, newParams, 1, n * Gaussian2DFunction.PARAMETERS_PER_PEAK);
		}
		final int left = nPeaks - (n + 1);
		if (left > 0)
		{
			System.arraycopy(params, (n + 1) * Gaussian2DFunction.PARAMETERS_PER_PEAK + 1, newParams,
					n * Gaussian2DFunction.PARAMETERS_PER_PEAK + 1, left * Gaussian2DFunction.PARAMETERS_PER_PEAK);
		}
		return newParams;
	}

	/**
	 * Join the parameters from the two parameter arrays. The background is taken from the first array. This function
	 * can be used to append the parameters from a pre-computed function onto the parameters from a fit result to build
	 * the entire function parameters.
	 *
	 * @param params1
	 *            the params 1
	 * @param params2
	 *            the params 2
	 * @return the double[]
	 */
	@SuppressWarnings("unused")
	private static double[] joinParams(double[] params1, double[] params2)
	{
		double[] params = new double[params1.length + params2.length - 1];
		System.arraycopy(params1, 0, params, 0, params1.length);
		System.arraycopy(params2, 1, params, params1.length, params2.length - 1);
		return params;
	}

	/**
	 * @return The background estimate when fitting a single peak
	 */
	@SuppressWarnings("unused")
	private float getSingleFittingBackground()
	{
		final float background;
		if (useFittedBackground && fittedBackground.getN() != 0)
		{
			// Use the average background from all results
			background = (float) (fittedBackground.getMean());
		}
		else
		{
			if (this.noise != 0)
			{
				// Initial guess using the noise (assuming all noise is from Poisson background).
				// EMCCD will have increase noise by a factor of sqrt(2)
				CalibrationReader r = new CalibrationReader(fitConfig.getCalibration());
				double gain = (fitConfig.isFitCameraCounts()) ? r.getCountPerPhoton() : 1;
				background = (float) (PeakResultHelper.noiseToLocalBackground(noise, gain, r.isEMCCD()));
			}
			else
			{
				// Initial guess using the data estimator
				background = estimateBackground();
			}
		}
		return background;
	}

	private void resetNeighbours()
	{
		clearGridCache();
		candidateNeighbourCount = 0;
		fittedNeighbourCount = 0;
		allNeighbours = null;
		allFittedNeighbours = null;
	}

	private CandidateList findNeighbours(Candidate candidate)
	{
		if (allNeighbours == null)
		{
			// Using the neighbour grid
			allNeighbours = gridManager.getCandidateNeighbours(candidate);
			allNeighbours.sort();
		}
		return allNeighbours;
	}

	private CandidateList findPeakNeighbours(Candidate candidate)
	{
		if (allFittedNeighbours == null)
		{
			// Using the neighbour grid 
			allFittedNeighbours = gridManager.getFittedNeighbours(candidate.x, candidate.y);
		}
		return allFittedNeighbours;
	}

	/**
	 * Search for any neighbours within a set height of the specified peak that is within the search region bounds.
	 *
	 * @param regionBounds
	 *            the region bounds
	 * @param candidateId
	 *            the candidate index
	 * @param background
	 *            The background in the region
	 * @return The number of neighbours
	 */
	private int findNeighbours(Rectangle regionBounds, int candidateId, float background)
	{
		final int xmin = regionBounds.x;
		final int xmax = xmin + regionBounds.width - 1;
		final int ymin = regionBounds.y;
		final int ymax = ymin + regionBounds.height - 1;

		final Candidate spot = candidates.get(candidateId);

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
		CandidateList neighbours = findNeighbours(spot);
		for (int i = 0; i < neighbours.getSize(); i++)
		{
			final Candidate neighbour = neighbours.get(i);
			if (isFit(neighbour.index))
				continue;
			if (canIgnore(neighbour.x, neighbour.y, xmin, xmax, ymin, ymax, neighbour.intensity, heightThreshold))
				continue;
			candidateNeighbours[candidateNeighbourCount++] = neighbour;
		}

		// XXX Debugging
		//		int c = 0;
		//		// Processing all lower spots.
		//		//for (int i = n + 1; i < candidates.getSize(); i++)
		//		// Processing all spots.
		//		for (int i = 0; i < this.candidates.getSize(); i++)
		//		{
		//			if (i == n || isFit(i))
		//				continue;
		//			if (canIgnore(this.candidates.get(i).x, this.candidates.get(i).y, xmin, xmax, ymin, ymax,
		//					this.candidates.get(i).intensity, heightThreshold))
		//				continue;
		//			//neighbourIndices[c++] = i;
		//			if (neighbourIndices[c++] != i)
		//				throw new RuntimeException("invalid grid neighbours");
		//		}

		// Check all existing maxima. 

		fittedNeighbourCount = 0;
		if (fittedBackground.getN() != 0)
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

			neighbours = findPeakNeighbours(spot);

			for (int i = 0; i < neighbours.getSize(); i++)
			{
				final Candidate neighbour = neighbours.get(i);
				// No height threshold check as this is a validated peak
				final double xw = 2 * neighbour.params[Gaussian2DFunction.X_SD];
				final double yw = 2 * neighbour.params[Gaussian2DFunction.Y_SD];
				final double x = neighbour.params[Gaussian2DFunction.X_POSITION];
				final double y = neighbour.params[Gaussian2DFunction.Y_POSITION];
				if (intersects(x0min, y0min, x0max, y0max, x - xw, y - yw, x + xw, y + yw))
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
	 * Get an estimate of the background level using the median of the image.
	 * 
	 * @return The background level
	 */
	private float estimateBackground()
	{
		createDataEstimator();
		// Use the median
		return dataEstimator.getPercentile(50);
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
				sum += result.getBackground();
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
		createDataEstimator();
		return estimateNoise(dataEstimator, FitProtosHelper.convertNoiseEstimatorMethod(config.getNoiseMethod()));
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
	public static float estimateNoise(float[] data, int width, int height, NoiseEstimatorMethod method)
	{
		// Do the same logic as the non-static method 
		DataEstimator dataEstimator = newDataEstimator(data, width, height);
		return estimateNoise(dataEstimator, FitProtosHelper.convertNoiseEstimatorMethod(method));
	}

	private static float estimateNoise(DataEstimator dataEstimator, NoiseEstimator.Method method)
	{
		// No methods using a background region are good so we just use the global noise estimate
		//if (dataEstimator.isBackgroundRegion())
		//	return dataEstimator.getNoise();		
		return dataEstimator.getNoise(method);
	}

	private void createDataEstimator()
	{
		if (dataEstimator == null)
		{
			int width = job.getBounds().width;
			int height = job.getBounds().height;
			dataEstimator = newDataEstimator(data, width, height);
		}
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
		// This is the total number of results we produce. 
		// Note that although we have may a maximum candidate less than the length
		// of the candidate list, we continue processing candidates if we have an 
		// estimate. This is possibly a candidate that was a good fit but was not labelled
		// as a new result because of drift to another yet-to-be-processed candidate.
		return candidates.getLength();
	}

	/**
	 * Provide the multi-path fit results dynamically
	 */
	private class DynamicMultiPathFitResult extends MultiPathFitResult
	{
		/** The constant for no Quadrant Analysis score */
		private static final double NO_QA_SCORE = -1;

		final ImageExtractor ie, ie2;
		boolean dynamic;
		Rectangle regionBounds;
		double[] region, region2, var_g2;
		CandidateSpotFitter spotFitter;
		FitType fitType;
		boolean isValid;
		@SuppressWarnings("unused")
		int extra = 0;

		public DynamicMultiPathFitResult(ImageExtractor ie, ImageExtractor ie2, boolean dynamic)
		{
			this.frame = FitWorker.this.slice;
			this.width = cc.dataBounds.width;
			this.height = cc.dataBounds.height;
			this.ie = ie;
			this.ie2 = ie2;
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

			// Only provide results if below the max candidate ID or we have a valid estimate
			if (candidateId < candidates.getSize())
				isValid = true;
			else if (isValid(candidateId))
			{
				extra++; // Count these for debugging
				isValid = true;
			}
			else
				isValid = false;

			if (isValid)
			{
				// Set fitting region
				regionBounds = ie.getBoxRegionBounds(candidates.get(candidateId).x, candidates.get(candidateId).y,
						fitting);
				region = ie.crop(regionBounds, region);
				region2 = ie2.crop(regionBounds, region2);

				cc.setRegionBounds(regionBounds);

				// Set up per-pixel noise
				if (cameraModel.isPerPixelModel())
				{
					// Note: The region bounds are relative to the data bounds origin so 
					// convert them to absolute
					Rectangle bounds = new Rectangle(regionBounds);
					bounds.x += cc.dataBounds.x;
					bounds.y += cc.dataBounds.y;
					float[] v = cameraModel.getNormalisedVariance(bounds);
					// Convert to double
					if (var_g2 == null || var_g2.length != v.length)
					{
						var_g2 = SimpleArrayUtils.toDouble(v);
					}
					else
					{
						for (int i = 0; i < v.length; i++)
							var_g2[i] = v[i];
					}
				}

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
				spotFitter = new CandidateSpotFitter(gf, factory, region, region2, var_g2, regionBounds, candidateId);
			}
		}

		@Override
		public FitResult getMultiFitResult()
		{
			FitResult result = super.getMultiFitResult();
			if (result == null && isValid)
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
			if (score == NO_QA_SCORE && isValid)
			{
				score = spotFitter.getQAScoreMulti();
				this.setMultiQAScore(score);
			}
			return score;
		}

		@Override
		public FitResult getMultiDoubletFitResult()
		{
			FitResult result = super.getMultiDoubletFitResult();
			if (result == null && isValid)
			{
				result = spotFitter.getResultDoubletMulti(config.getResidualsThreshold());
				setMultiDoubletFitResult(result);
				fitType.setMultiDoublet(spotFitter.computedDoubletMulti);
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
			if (result == null && isValid)
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
			if (score == NO_QA_SCORE && isValid)
			{
				score = spotFitter.getQAScoreSingle();
				this.setSingleQAScore(score);
			}
			return score;
		}

		@Override
		public FitResult getDoubletFitResult()
		{
			FitResult result = super.getDoubletFitResult();
			if (result == null && isValid)
			{
				result = spotFitter.getResultDoubletSingle(config.getResidualsThreshold());
				setDoubletFitResult(result);
				fitType.setDoublet(spotFitter.computedDoubletSingle);
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
			// This means fitting and then selection of the best filter settings
			// must be iterated until convergence to ensure the fitting+filter is 
			// optimum.

			if (dynamicMultiPathFitResult.isValid)
			{
				// Calling the spot fitter with zero residuals will force the result to be computed if possible
				dynamicMultiPathFitResult.spotFitter.getResultDoubletMulti(0);
				dynamicMultiPathFitResult.spotFitter.getResultDoubletSingle(0);

				// Now update the result if they were previously null
				dynamicMultiPathFitResult.getMultiFitResult();
				dynamicMultiPathFitResult.getMultiQAScore();
				dynamicMultiPathFitResult.getMultiDoubletFitResult();
				dynamicMultiPathFitResult.getSingleFitResult();
				dynamicMultiPathFitResult.getSingleQAScore();
				dynamicMultiPathFitResult.getDoubletFitResult();
			}

			job.setMultiPathFitResult(index, dynamicMultiPathFitResult.copy(false));
		}

		// Send the actual results to the neighbour grid
		if (flushToGrid())
		{
			// Count if there were any new results
			success++;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IMultiPathFitResults#getTotalCandidates()
	 */
	public int getTotalCandidates()
	{
		// This is the total number of candidates Ids we may produce		
		return candidates.getLength();
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

		// Add to the slice results.
		final PreprocessedPeakResult[] results = selectedResult.results;
		if (results == null)
		{
			if (logger != null)
			{
				final int candidateId = dynamicMultiPathFitResult.candidateId;
				//@formatter:off
				logger.info("Not fit %d (%d,%d) %s", candidateId, 
						cc.fromDataToGlobalX(candidates.get(candidateId).x),
						cc.fromDataToGlobalY(candidates.get(candidateId).y),
						(selectedResult.fitResult!=null && selectedResult.fitResult.data!=null) ?
						((FitResult) selectedResult.fitResult.data).getStatus().toString() : "");								
				//@formatter:on
			}
			return;
		}

		final int candidateId = dynamicMultiPathFitResult.candidateId;

		final FitResult fitResult = (FitResult) selectedResult.fitResult.data;

		// The background for each result was the local background. We want the fitted global background
		final float background = (float) fitResult.getParameters()[0];
		final double[] dev = fitResult.getParameterDeviations();

		if (queueSize != 0)
			throw new RuntimeException("There are results queued already!");

		for (int i = 0; i < results.length; i++)
		{
			PreprocessedPeakResult peak = results[i];
			if (peak.isExistingResult())
				continue;
			if (peak.isNewResult())
			{
				final double[] p = peak.toGaussian2DParameters();

				// Store slice results relative to the data frame (not the global bounds)
				// Convert back so that 0,0 is the top left of the data bounds
				p[Gaussian2DFunction.X_POSITION] -= cc.dataBounds.x;
				p[Gaussian2DFunction.Y_POSITION] -= cc.dataBounds.y;

				final float[] params = new float[p.length];
				params[Gaussian2DFunction.BACKGROUND] = background;
				for (int j = 1; j < p.length; j++)
					params[j] = (float) p[j];

				final float[] paramDevs;
				if (dev == null)
				{
					paramDevs = null;
				}
				else
				{
					paramDevs = new float[p.length];
					paramDevs[Gaussian2DFunction.BACKGROUND] = (float) dev[Gaussian2DFunction.BACKGROUND];
					final int offset = peak.getId() * Gaussian2DFunction.PARAMETERS_PER_PEAK;
					for (int j = 1; j < p.length; j++)
						paramDevs[j] = (float) dev[offset + j];
				}

				addSingleResult(peak.getCandidateId(), params, paramDevs, fitResult.getError(), peak.getNoise(),
						peak.getLocationVarianceCRLB());

				if (logger != null)
				{
					// Show the shift, signal and width spread
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
				case FAILED_TO_ESTIMATE_WIDTH:
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
				int npeaks = peakParams.length / Gaussian2DFunction.PARAMETERS_PER_PEAK;
				for (int i = 0; i < npeaks; i++)
				{
					peakParams[i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_POSITION] += cc
							.fromFitRegionToGlobalX();
					peakParams[i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.Y_POSITION] += cc
							.fromFitRegionToGlobalY();
					if (fitConfig.isAngleFitting())
						peakParams[i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.ANGLE] *= 180.0 /
								Math.PI;
				}
			}
			final int x = candidates.get(candidateId).x;
			final int y = candidates.get(candidateId).y;
			logger2.debug("%d:%d [%d,%d] %s (%s) = %s\n", slice, candidateId, cc.fromDataToGlobalX(x),
					cc.fromDataToGlobalY(y), fitResult.getStatus(), fitResult.getStatusData(),
					Arrays.toString(peakParams));
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
		return candidates.get(candidateId).fit;
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
