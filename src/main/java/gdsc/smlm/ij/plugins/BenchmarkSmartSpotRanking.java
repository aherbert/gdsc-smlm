package gdsc.smlm.ij.plugins;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import gdsc.core.ij.Utils;
import gdsc.core.match.ClassificationResult;
import gdsc.core.match.Coordinate;
import gdsc.core.match.FractionClassificationResult;
import gdsc.core.threshold.AutoThreshold;
import gdsc.core.threshold.FloatHistogram;
import gdsc.core.threshold.Histogram;
import gdsc.core.utils.NoiseEstimator.Method;
import gdsc.core.utils.Statistics;
import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.fitting.FitConfiguration;
import gdsc.smlm.ij.plugins.BenchmarkSpotFilter.FilterResult;
import gdsc.smlm.ij.plugins.BenchmarkSpotFilter.ScoredSpot;
import gdsc.smlm.ij.plugins.ResultsMatchCalculator.PeakResultPoint;
import gdsc.smlm.results.MemoryPeakResults;
import ij.IJ;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.text.TextWindow;

/**
 * Attempt to classify the spot candidates into those that do match a result and those that do not. This smart ranking
 * can be used to ensure all the good candidates are processed per frame before the fail limit takes effect.
 */
public class BenchmarkSmartSpotRanking implements PlugIn
{
	private static final String TITLE = "Smart Spot Ranking";

	static FitConfiguration fitConfig;
	private static FitEngineConfiguration config;
	static
	{
		fitConfig = new FitConfiguration();
		config = new FitEngineConfiguration(fitConfig);
		// Set some default fit settings here ...
		// Ensure all candidates are fitted
		config.setFailuresLimit(-1);
		fitConfig.setFitValidation(true);
		fitConfig.setMinPhotons(1); // Do not allow negative photons 
		fitConfig.setCoordinateShiftFactor(0);
		fitConfig.setPrecisionThreshold(0);
		fitConfig.setMinWidthFactor(0);
		fitConfig.setWidthFactor(0);

		fitConfig.setBackgroundFitting(true);
		fitConfig.setMinIterations(0);
		fitConfig.setNoise(0);
		config.setNoiseMethod(Method.QUICK_RESIDUALS_LEAST_MEAN_OF_SQUARES);
	}

	private static AutoThreshold.Method[] thresholdMethods;
	private static boolean[] thresholdMethodOptions;
	private static String[] thresholdMethodNames;
	static
	{
		thresholdMethods = AutoThreshold.Method.values();
		thresholdMethodOptions = new boolean[thresholdMethods.length];
		thresholdMethodNames = new String[thresholdMethods.length];
		for (int i = 0; i < thresholdMethods.length; i++)
		{
			thresholdMethodNames[i] = thresholdMethods[i].name;
			thresholdMethodOptions[i] = true;
		}
		// Turn some off
		// These often fail to converge
		thresholdMethodOptions[AutoThreshold.Method.INTERMODES.ordinal()] = false;
		thresholdMethodOptions[AutoThreshold.Method.MINIMUM.ordinal()] = false;
		// These are slow
		thresholdMethodOptions[AutoThreshold.Method.SHANBHAG.ordinal()] = false;
		thresholdMethodOptions[AutoThreshold.Method.RENYI_ENTROPY.ordinal()] = false;
	}
	private AutoThreshold.Method[] methods = null;

	private static double fractionPositives = 100;
	private static double fractionNegativesAfterAllPositives = 50;
	private static int negativesAfterAllPositives = 10;
	private static boolean selectMethods = true;

	private boolean extraOptions = false;

	private static TextWindow summaryTable = null;

	private MemoryPeakResults results;
	private CreateData.SimulationParameters simulationParameters;

	private static HashMap<Integer, ArrayList<Coordinate>> actualCoordinates = null;
	private static HashMap<Integer, FilterCandidates> filterCandidates;
	private static double fP, fN;
	private static int nP, nN;

	static int lastId = -1, lastFilterId = -1;
	private static double lastFractionPositives = -1;
	private static double lastFractionNegativesAfterAllPositives = -1;
	private static int lastNegativesAfterAllPositives = -1;

	// Allow other plugins to access the results
	static int rankResultsId = 0;
	static HashMap<Integer, RankResults> rankResults;
	static double distanceInPixels;
	static double lowerDistanceInPixels;
	static double candidateTN, candidateFN;

	private class FilterCandidates
	{
		// Integer counts of positives (matches) and negatives
		final int p, n;
		// Double sums of the fractions match score and antiscore 
		final double np, nn;
		final ScoredSpot[] spots;

		public FilterCandidates(int p, int n, double np, double nn, ScoredSpot[] spots)
		{
			this.p = p;
			this.n = n;
			this.np = np;
			this.nn = nn;
			this.spots = spots;
		}
	}

	private class RankResult
	{
		@SuppressWarnings("unused")
		final AutoThreshold.Method m;
		final float t;
		final FractionClassificationResult f;
		final ClassificationResult c;
		/**
		 * Store details about the spots that were accepted
		 */
		@SuppressWarnings("unused")
		final boolean[] good;
		final long time;

		public RankResult(AutoThreshold.Method m, float t, FractionClassificationResult f, ClassificationResult c,
				boolean[] good, long time)
		{
			this.m = m;
			this.t = t;
			this.f = f;
			this.c = c;
			this.good = good;
			this.time = time;
		}
	}

	private class RankResults
	{
		@SuppressWarnings("unused")
		final ScoredSpot[] spots;
		/** Store the z-position of the actual spots for later analysis. Size is the number of actual spots */
		@SuppressWarnings("unused")
		final double[] zPosition;

		ArrayList<RankResult> results = new ArrayList<RankResult>();

		public RankResults(ScoredSpot[] spots, double[] zPosition)
		{
			this.spots = spots;
			this.zPosition = zPosition;
		}

	}

	/**
	 * Used to allow multi-threading of the fitting method
	 */
	private class Worker implements Runnable
	{
		volatile boolean finished = false;
		final BlockingQueue<Integer> jobs;
		final HashMap<Integer, ArrayList<Coordinate>> actualCoordinates;
		final HashMap<Integer, FilterCandidates> filterCandidates;
		final HashMap<Integer, RankResults> results;

		public Worker(BlockingQueue<Integer> jobs, HashMap<Integer, ArrayList<Coordinate>> actualCoordinates,
				HashMap<Integer, FilterCandidates> filterCandidates)
		{
			this.jobs = jobs;
			this.actualCoordinates = actualCoordinates;
			this.filterCandidates = filterCandidates;
			this.results = new HashMap<Integer, RankResults>();
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
					Integer job = jobs.take();
					if (job == null || job.intValue() < 0 || finished)
						break;
					run(job.intValue());
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
			}
		}

		private void run(int frame)
		{
			if (Utils.isInterrupted())
			{
				finished = true;
				return;
			}

			showProgress();

			FilterCandidates candidates = filterCandidates.get(frame);

			// Extract the spot intensities
			final ScoredSpot[] spots = candidates.spots;
			float[] intensity = new float[spots.length];
			for (int i = 0; i < spots.length; i++)
			{
				intensity[i] = spots[i].spot.intensity;
			}

			Coordinate[] actual = ResultsMatchCalculator.getCoordinates(actualCoordinates, frame);
			final double[] zPosition = new double[actual.length];
			for (int i = 0; i < actual.length; i++)
			{
				PeakResultPoint p = (PeakResultPoint) actual[i];
				zPosition[i] = p.peakResult.error;
			}

			RankResults results = new RankResults(spots, zPosition);
			this.results.put(frame, results);

			long t1 = System.nanoTime();
			FloatHistogram histogram = FloatHistogram.buildHistogram(intensity, true);
			// Only compact once
			final int bins = 4096;
			Histogram histogram2 = histogram.compact(bins);
			t1 = System.nanoTime() - t1;

			for (AutoThreshold.Method m : methods)
			{
				long t2 = System.nanoTime();
				final float t = histogram2.getAutoThreshold(m);
				t2 = System.nanoTime() - t2;

				// Score
				final boolean[] good = new boolean[spots.length];
				double tp = 0;
				double fp = 0;
				double tn = 0;
				double fn = 0;
				int itp = 0;
				int ifp = 0;
				int itn = 0;
				int ifn = 0;
				for (int i = 0; i < spots.length; i++)
				{
					if (spots[i].spot.intensity >= t)
					{
						good[i] = true;
						tp += spots[i].getScore();
						fp += spots[i].antiScore();
						if (spots[i].match)
							itp++;
						else
							ifp++;
					}
					else
					{
						fn += spots[i].getScore();
						tn += spots[i].antiScore();
						if (spots[i].match)
							ifn++;
						else
							itn++;
					}
				}

				// Store the results using a copy of the original (to preserve the candidates for repeat analysis)
				results.results.add(new RankResult(m, t, new FractionClassificationResult(tp, fp, tn, fn),
						new ClassificationResult(itp, ifp, itn, ifn), good, t1 + t2));
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		extraOptions = Utils.isExtraOptions();

		simulationParameters = CreateData.simulationParameters;
		if (simulationParameters == null)
		{
			IJ.error(TITLE, "No benchmark spot parameters in memory");
			return;
		}

		results = MemoryPeakResults.getResults(CreateData.CREATE_DATA_IMAGE_TITLE + " (Create Data)");
		if (results == null)
		{
			IJ.error(TITLE, "No benchmark results in memory");
			return;
		}
		if (BenchmarkSpotFilter.filterResults == null)
		{
			IJ.error(TITLE, "No benchmark spot candidates in memory");
			return;
		}
		if (BenchmarkSpotFilter.simulationId != simulationParameters.id)
		{
			IJ.error(TITLE, "Update the benchmark spot candidates for the latest simulation");
			return;
		}

		if (!showDialog())
			return;

		run();
	}

	private boolean showDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage(String.format(
				"Rank candidate spots in the benchmark image created by " + CreateData.TITLE +
						" plugin\nand identified by the " + BenchmarkSpotFilter.TITLE +
						" plugin.\nPSF width = %s nm (Square pixel adjustment = %s nm)\n \nConfigure the fitting:",
				Utils.rounded(simulationParameters.s), Utils.rounded(getSa())));

		gd.addSlider("Fraction_positives", 50, 100, fractionPositives);
		gd.addSlider("Fraction_negatives_after_positives", 0, 100, fractionNegativesAfterAllPositives);
		gd.addSlider("Min_negatives_after_positives", 0, 10, negativesAfterAllPositives);
		gd.addCheckbox("Select_methods", selectMethods);

		// Collect options for fitting that may effect ranking
		final double sa = getSa();
		gd.addNumericField("Initial_StdDev", sa / simulationParameters.a, 3);
		gd.addSlider("Fitting_width", 2, 4.5, config.getFitting());
		gd.addCheckbox("Include_neighbours", config.isIncludeNeighbours());
		gd.addSlider("Neighbour_height", 0.01, 1, config.getNeighbourHeightThreshold());

		// TODO - add options for spot ranking here, e.g. the threshold methods

		// Output options
		//		gd.addCheckbox("Show_score_histograms", showFilterScoreHistograms);
		//		gd.addCheckbox("Show_correlation", showCorrelation);
		//		gd.addCheckbox("Plot_rank_by_intensity", rankByIntensity);
		//		gd.addCheckbox("Save_filter_range", saveFilterRange);

		if (extraOptions)
		{
		}

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		fractionPositives = Math.abs(gd.getNextNumber());
		fractionNegativesAfterAllPositives = Math.abs(gd.getNextNumber());
		negativesAfterAllPositives = (int) Math.abs(gd.getNextNumber());
		selectMethods = gd.getNextBoolean();

		fitConfig.setInitialPeakStdDev(gd.getNextNumber());
		config.setFitting(gd.getNextNumber());
		config.setIncludeNeighbours(gd.getNextBoolean());
		config.setNeighbourHeightThreshold(gd.getNextNumber());

		//		showFilterScoreHistograms = gd.getNextBoolean();
		//		showCorrelation = gd.getNextBoolean();
		//		rankByIntensity = gd.getNextBoolean();
		//		saveFilterRange = gd.getNextBoolean();

		if (extraOptions)
		{
		}

		if (gd.invalidNumber())
			return false;

		int count = 0;
		if (selectMethods)
		{
			methods = new AutoThreshold.Method[thresholdMethods.length];

			gd = new GenericDialog(TITLE);
			gd.addHelp(About.HELP_URL);
			for (int i = 0; i < thresholdMethods.length; i++)
				gd.addCheckbox(thresholdMethodNames[i], thresholdMethodOptions[i]);

			gd.showDialog();

			if (gd.wasCanceled())
				return false;

			for (int i = 0; i < thresholdMethods.length; i++)
			{
				thresholdMethodOptions[i] = gd.getNextBoolean();
				if (thresholdMethodOptions[i])
					methods[count++] = thresholdMethods[i];
			}

			methods = Arrays.copyOf(methods, count);
		}
		else
		{
			// Do them all
			methods = thresholdMethods.clone();
			count = methods.length;
		}

		if (count == 0)
		{
			IJ.error(TITLE, "No threshold methods selected");
			return false;
		}

		return true;
	}

	/** The total progress. */
	int progress, stepProgress, totalProgress;

	/**
	 * Show progress.
	 */
	private synchronized void showProgress()
	{
		if (++progress % stepProgress == 0)
		{
			if (Utils.showStatus("Frame: " + progress + " / " + totalProgress))
				IJ.showProgress(progress, totalProgress);
		}
	}

	private void run()
	{
		// Extract all the results in memory into a list per frame. This can be cached
		boolean refresh = false;
		if (lastId != simulationParameters.id)
		{
			// Do not get integer coordinates
			// The Coordinate objects will be PeakResultPoint objects that store the original PeakResult
			// from the MemoryPeakResults
			actualCoordinates = ResultsMatchCalculator.getCoordinates(results.getResults(), false);
			lastId = simulationParameters.id;
			refresh = true;
		}

		// Extract all the candidates into a list per frame. This can be cached if the settings have not changed
		if (refresh || lastFilterId != BenchmarkSpotFilter.filterResultsId ||
				lastFractionPositives != fractionPositives ||
				lastFractionNegativesAfterAllPositives != fractionNegativesAfterAllPositives ||
				lastNegativesAfterAllPositives != negativesAfterAllPositives)
		{
			filterCandidates = subsetFilterResults(BenchmarkSpotFilter.filterResults);

			lastFilterId = BenchmarkSpotFilter.filterResultsId;
			lastFractionPositives = fractionPositives;
			lastFractionNegativesAfterAllPositives = fractionNegativesAfterAllPositives;
			lastNegativesAfterAllPositives = negativesAfterAllPositives;
		}

		// Create a pool of workers
		final int nThreads = Prefs.getThreads();
		BlockingQueue<Integer> jobs = new ArrayBlockingQueue<Integer>(nThreads * 2);
		List<Worker> workers = new LinkedList<Worker>();
		List<Thread> threads = new LinkedList<Thread>();
		for (int i = 0; i < nThreads; i++)
		{
			Worker worker = new Worker(jobs, actualCoordinates, filterCandidates);
			Thread t = new Thread(worker);
			workers.add(worker);
			threads.add(t);
			t.start();
		}

		// Process the frames
		totalProgress = filterCandidates.size();
		stepProgress = Utils.getProgressInterval(totalProgress);
		progress = 0;
		for (int i : filterCandidates.keySet())
		{
			put(jobs, i);
		}
		// Finish all the worker threads by passing in a null job
		for (int i = 0; i < threads.size(); i++)
		{
			put(jobs, -1);
		}

		// Wait for all to finish
		for (int i = 0; i < threads.size(); i++)
		{
			try
			{
				threads.get(i).join();
			}
			catch (InterruptedException e)
			{
				e.printStackTrace();
			}
		}
		threads.clear();

		IJ.showProgress(1);

		if (Utils.isInterrupted())
		{
			IJ.showStatus("Aborted");
			return;
		}

		IJ.showStatus("Collecting results ...");

		rankResultsId++;
		rankResults = new HashMap<Integer, RankResults>();
		for (Worker w : workers)
		{
			rankResults.putAll(w.results);
		}

		summariseResults(rankResults);

		IJ.showStatus("");
	}

	/**
	 * Extract all the filter candidates in order until the desired number of positives have been reached and the number
	 * of negatives matches the configured parameters.
	 * 
	 * @param filterResults
	 * @return The filter candidates
	 */
	private HashMap<Integer, FilterCandidates> subsetFilterResults(HashMap<Integer, FilterResult> filterResults)
	{
		// Convert fractions from percent 
		final double f1 = Math.min(1, fractionPositives / 100.0);
		final double f2 = fractionNegativesAfterAllPositives / 100.0;

		HashMap<Integer, FilterCandidates> subset = new HashMap<Integer, FilterCandidates>();
		fP = fN = 0;
		nP = nN = 0;
		for (Entry<Integer, FilterResult> result : filterResults.entrySet())
		{
			FilterResult r = result.getValue();

			// Determine the number of positives to find. This score may be fractional.
			fP += r.result.getTP();
			fN += r.result.getFP();

			// Q. Is r.result.getTP() not the same as the total of r.spots[i].match
			// A. Not if we used fractional scoring.

			for (ScoredSpot spot : r.spots)
			{
				if (spot.match)
					nP++;
				else
					nN++;
			}

			// Make the target use the fractional score
			final double targetP = r.result.getTP() * f1;

			// Count the number of positive & negatives
			int p = 0, n = 0;
			double np = 0, nn = 0;

			boolean reachedTarget = false;
			int nAfter = 0;

			int count = 0;
			for (ScoredSpot spot : r.spots)
			{
				count++;
				nn += spot.antiScore();
				if (spot.match)
				{
					np += spot.getScore();
					p++;
					if (!reachedTarget)
					{
						reachedTarget = np >= targetP;
					}
				}
				else
				{
					n++;
					if (reachedTarget)
					{
						nAfter++;
					}
				}

				if (reachedTarget)
				{
					// Check if we have reached both the limits
					if (nAfter >= negativesAfterAllPositives && (double) n / (n + p) >= f2)
						break;
				}
			}

			// Debug
			//System.out.printf("Frame %d : %.1f / (%.1f + %.1f). p=%d, n=%d, after=%d, f=%.1f\n", result.getKey().intValue(),
			//		r.result.getTP(), r.result.getTP(), r.result.getFP(), p, n,
			//		nAfter, (double) n / (n + p));

			subset.put(result.getKey(), new FilterCandidates(p, n, np, nn, Arrays.copyOf(r.spots, count)));
		}
		return subset;
	}

	private void put(BlockingQueue<Integer> jobs, int i)
	{
		try
		{
			jobs.put(i);
		}
		catch (InterruptedException e)
		{
			throw new RuntimeException("Unexpected interruption", e);
		}
	}

	private void summariseResults(HashMap<Integer, RankResults> rankResults)
	{
		createTable();

		// Summarise the ranking results. 
		StringBuilder sb = new StringBuilder(BenchmarkSpotFilter.resultPrefix);

		// nP and nN is the fractional score of the spot candidates 
		addCount(sb, nP + nN);
		addCount(sb, nP);
		addCount(sb, nN);
		addCount(sb, fP);
		addCount(sb, fN);

		double tp = 0;
		double fp = 0;
		int cTP = 0, cFP = 0;
		for (FilterCandidates result : filterCandidates.values())
		{
			tp += result.np;
			fp += result.nn;
			cTP += result.p;
			cFP += result.n;
		}

		//		// This should be the same
		//		double tp2 = 0;
		//		double fp2 = 0;
		//		int cTP2 = 0, cFP2 = 0;
		//		for (RankResults rr : rankResults.values())
		//		{
		//			for (ScoredSpot spot : rr.spots)
		//			{
		//				if (spot.match)
		//					cTP2++;
		//				else
		//					cFP2++;
		//				tp2 += spot.getScore();
		//				fp2 += spot.antiScore();
		//			}
		//		}
		//		if (tp != tp2 || fp != fp2 || cTP != cTP2 || cFP != cFP2)
		//			System.out.println("Error counting");

		// The fraction of positive and negative candidates that were included
		add(sb, (100.0 * cTP) / nP);
		add(sb, (100.0 * cFP) / nN);

		// Add counts of the the candidates
		add(sb, cTP + cFP);
		add(sb, cTP);
		add(sb, cFP);

		// Add fractional counts of the the candidates
		add(sb, tp + fp);
		add(sb, tp);
		add(sb, fp);

		String resultPrefix = sb.toString();

		// TODO

		// Pre-compute the results and have optional sort

		// Add more scoring metrics

		// Add good label to spot candidates and have the benchmark spot filter respect this before applying the fail count limit.
		// This may just involve setting the fail count to zero in all the results that are good, given that spots are in order
		// and we are just thresholding by intensity

		// Allow using the fitted results from benchmark spot fit. Will it make a difference if we fit the candidates (some will fail
		// if weak).
		// Can this be done by allowing the user to select the input (spot candidates or fitted positions)?

		for (int i = 0; i < methods.length; i++)
		{
			tp = 0;
			fp = 0;
			double tn = 0;
			int itp = 0;
			int ifp = 0;
			int itn = 0;
			Statistics s = new Statistics();
			long time = 0;

			for (RankResults rr : rankResults.values())
			{
				RankResult r = rr.results.get(i);
				s.add(r.t);
				time += r.time;
				tp += r.f.getTP();
				fp += r.f.getFP();
				tn += r.f.getTN();
				itp += r.c.getTP();
				ifp += r.c.getFP();
				itn += r.c.getTN();
			}

			sb.setLength(0);
			sb.append(resultPrefix);
			add(sb, methods[i].name);
			sb.append('\t').append(Utils.rounded(s.getMean())).append(" +/- ")
					.append(Utils.rounded(s.getStandardDeviation()));
			add(sb, Utils.timeToString(time / 1e6));

			// TP are all accepted candidates that can be matched to a spot
			// FP are all accepted candidates that cannot be matched to a spot
			// TN are all accepted candidates that cannot be matched to a spot
			// FN = The number of missed spots

			// Raw counts of match or no-match
			addScores(sb, new FractionClassificationResult(itp, ifp, itn, simulationParameters.molecules - itp));

			// Fractional scoring
			addScores(sb, new FractionClassificationResult(tp, fp, tn, simulationParameters.molecules - tp));

			summaryTable.append(sb.toString());
		}
	}

	private static void add(StringBuilder sb, String value)
	{
		sb.append("\t").append(value);
	}

	private static void add(StringBuilder sb, int value)
	{
		sb.append("\t").append(value);
	}

	private static void add(StringBuilder sb, double value)
	{
		add(sb, Utils.rounded(value));
	}

	private static void addCount(StringBuilder sb, double value)
	{
		// Check if the double holds an integer count
		if ((int) value == value)
		{
			sb.append("\t").append((int) value);
		}
		else
		{
			// Otherwise add the counts using at least 2 dp
			if (value > 100)
				sb.append("\t").append(IJ.d2s(value));
			else
				add(sb, Utils.rounded(value));
		}
	}

	private void createTable()
	{
		if (summaryTable == null || !summaryTable.isVisible())
		{
			summaryTable = new TextWindow(TITLE, createHeader(false), "", 1000, 300);
			summaryTable.setVisible(true);
		}
	}

	private String createHeader(boolean extraRecall)
	{
		StringBuilder sb = new StringBuilder(BenchmarkSpotFilter.tablePrefix);
		sb.append("\t");
		sb.append("Spots\t");
		sb.append("nP\t");
		sb.append("nN\t");
		sb.append("fP\t");
		sb.append("fN\t");

		sb.append("% nP\t");
		sb.append("% nN\t");

		sb.append("cTotal\t");
		sb.append("cTP\t");
		sb.append("cFP\t");

		sb.append("cfTotal\t");
		sb.append("cfTP\t");
		sb.append("cfFP\t");

		sb.append("Method\t");
		sb.append("Threshold\t");
		sb.append("Time\t");

		addScoreColumns(sb, null);
		addScoreColumns(sb, "f ");

		return sb.toString();
	}

	private void addScoreColumns(StringBuilder sb, String prefix)
	{
		addScoreColumn(sb, prefix, "tp");
		addScoreColumn(sb, prefix, "fp");
		addScoreColumn(sb, prefix, "tn");
		addScoreColumn(sb, prefix, "fn");
		addScoreColumn(sb, prefix, "Recall");
		addScoreColumn(sb, prefix, "Precision");
		addScoreColumn(sb, prefix, "F1");
		addScoreColumn(sb, prefix, "Jaccard");
		addScoreColumn(sb, prefix, "MCC");

		// TODO - add other scores
	}

	private void addScores(StringBuilder sb, FractionClassificationResult m)
	{
		add(sb, m.getTP());
		add(sb, m.getFP());
		add(sb, m.getTN());
		add(sb, m.getFN());
		add(sb, m.getRecall());
		add(sb, m.getPrecision());
		add(sb, m.getF1Score());
		add(sb, m.getJaccard());
		add(sb, m.getMCC());

		// TODO - add other scores
	}

	private void addScoreColumn(StringBuilder sb, String prefix, String name)
	{
		if (prefix != null)
			sb.append(prefix);
		sb.append(name);
		sb.append("\t");
	}

	private double getSa()
	{
		final double sa = PSFCalculator.squarePixelAdjustment(simulationParameters.s, simulationParameters.a);
		return sa;
	}
}
