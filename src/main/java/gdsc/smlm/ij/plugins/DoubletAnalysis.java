/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
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
package gdsc.smlm.ij.plugins;

import java.awt.Checkbox;
import java.awt.Choice;
import java.awt.Color;
import java.awt.Rectangle;
import java.awt.TextField;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.util.FastMath;

import gdsc.core.ij.IJLogger;
import gdsc.core.ij.Utils;
import gdsc.core.logging.Logger;
import gdsc.core.match.BasePoint;
import gdsc.core.match.Coordinate;
import gdsc.core.match.MatchCalculator;
import gdsc.core.match.PointPair;
import gdsc.core.utils.ImageExtractor;
import gdsc.core.utils.Maths;
import gdsc.core.utils.RampedScore;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.core.utils.StoredDataStatistics;
import gdsc.core.utils.TextUtils;
import gdsc.smlm.data.config.CalibrationReader;
import gdsc.smlm.data.config.CalibrationWriter;
import gdsc.smlm.data.config.FitProtos.NoiseEstimatorMethod;
import gdsc.smlm.data.config.FitProtos.PrecisionMethod;
import gdsc.smlm.data.config.PSFProtosHelper;
import gdsc.smlm.data.config.TemplateProtos.TemplateSettings;
import gdsc.smlm.engine.FitConfiguration;
import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.engine.FitWorker;
import gdsc.smlm.engine.QuadrantAnalysis;
import gdsc.smlm.filters.MaximaSpotFilter;
import gdsc.smlm.filters.Spot;
import gdsc.smlm.fitting.FitResult;
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.FunctionSolver;
// TODO - add support for using the chi-squared distribution to generate a q-value for the fit
//import gdsc.smlm.fitting.WLSEFunctionSolver;
//import gdsc.smlm.fitting.MLEFunctionSolver;
import gdsc.smlm.fitting.FunctionSolverType;
import gdsc.smlm.fitting.Gaussian2DFitter;
import gdsc.smlm.fitting.LSEFunctionSolver;
import gdsc.smlm.fitting.MLEFunctionSolver;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.ij.IJImageSource;
import gdsc.smlm.ij.plugins.ResultsMatchCalculator.PeakResultPoint;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.IJImageConverter;
import gdsc.smlm.results.MemoryPeakResults;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.procedure.TIntProcedure;
import gnu.trove.procedure.TObjectProcedure;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Macro;
import ij.Prefs;
import ij.gui.ExtendedGenericDialog;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.Plot2;
import ij.gui.PlotWindow;
import ij.gui.PointRoi;
import ij.plugin.PlugIn;
import ij.plugin.WindowOrganiser;
import ij.text.TextWindow;

/**
 * Fits spots created by CreateData plugin.
 * <p>
 * Assigns results to filter candidates to determine if spots are either single or doublets or larger clusters. Outputs
 * a table of the single and double fit for each spot with metrics. This can be used to determine the best settings for
 * optimum doublet fitting and filtering.
 */
public class DoubletAnalysis implements PlugIn, ItemListener
{
	/*
	 * Note: 21-Oct-2016
	 *
	 * This plugin may be obsolete now that the BenchmarkSpotFit and BenchmarkFilterAnalysis plugins
	 * can handle singles, multiples and doublets together. This means that the residuals threshold can be
	 * optimised concurrently with the fail count and the filter. It is left within the codebase in case
	 * it is useful in the future.
	 */

	private static final String TITLE = "Doublet Analysis";
	private static FitConfiguration fitConfig, filterFitConfig;
	private static FitEngineConfiguration config;
	private static int lastId = 0;
	static
	{
		config = new FitEngineConfiguration();
		fitConfig = config.getFitConfiguration();

		// Set some default fit settings here ...
		// Ensure all candidates are fitted
		config.setFailuresLimit(-1);
		fitConfig.setSmartFilter(false);
		fitConfig.setDisableSimpleFilter(false);
		fitConfig.setMinPhotons(1); // Do not allow negative photons
		fitConfig.setCoordinateShiftFactor(0); // Disable
		fitConfig.setPrecisionThreshold(0);
		fitConfig.setMinWidthFactor(0);
		fitConfig.setWidthFactor(0);

		fitConfig.setNoise(0);
		config.setNoiseMethod(NoiseEstimatorMethod.QUICK_RESIDUALS_LEAST_MEAN_OF_SQUARES);

		fitConfig.setBackgroundFitting(true);
		fitConfig.setNotSignalFitting(false);
		fitConfig.setComputeDeviations(false);
		fitConfig.setComputeResiduals(true);

		//
		filterFitConfig = new FitConfiguration();
		filterFitConfig.setSmartFilter(false);
		filterFitConfig.setDisableSimpleFilter(false);
		filterFitConfig.setMinPhotons(0);
		filterFitConfig.setCoordinateShiftFactor(0);
		filterFitConfig.setPrecisionThreshold(0);
		filterFitConfig.setMinWidthFactor(0);
		filterFitConfig.setWidthFactor(0);
		filterFitConfig.setPrecisionMethod(PrecisionMethod.MORTENSEN);
	}

	private static boolean useBenchmarkSettings = false;
	private static double iterationIncrease = 1;
	private static boolean ignoreWithNeighbours = false;
	private static boolean showOverlay = false;
	private static boolean showHistograms = false;
	private static boolean showResults = false;
	private static boolean showJaccardPlot = true;
	private static boolean useMaxResiduals = true;
	private static double lowerDistance = 1;
	private static double matchDistance = 1;
	private static double signalFactor = 2;
	private static double lowerSignalFactor = 1;
	private static String[] MATCHING = { "Simple", "By residuals", "By candidate" };
	private static int matching = 0;

	private static boolean analysisUseBenchmarkSettings = false;
	private static double analysisDriftAngle = 45;
	private static double minGap = 0;
	private static boolean analysisShowResults = false;
	private static boolean analysisLogging = false;
	private static String analysisTitle = "";
	private static boolean saveTemplate = false;
	private static String templateFilename = "";

	static
	{
		String currentUsersHomeDir = System.getProperty("user.home");
		templateFilename = currentUsersHomeDir + File.separator + "gdsc.smlm" + File.separator + "template";
	}

	private static String[] SELECTION_CRITERIA = { "R2", "AIC", "BIC", "ML AIC", "ML BIC" };
	private static int selectionCriteria = 4;

	private static TextWindow summaryTable = null, resultsTable = null, analysisTable = null;
	private static ArrayList<DoubletResult> doubletResults;
	private ResidualsScore residualsScore;
	private static ResidualsScore _residualsScoreMax;
	private static ResidualsScore _residualsScoreAv;
	private static int numberOfMolecules;
	private static String analysisPrefix;
	private ImagePlus imp;
	private MemoryPeakResults results;
	private CreateData.SimulationParameters simulationParameters;

	private Choice textPSF;
	private Choice textDataFilterType;
	private Choice textDataFilter;
	private TextField textSmooth;
	private TextField textSearch;
	private TextField textBorder;
	private TextField textFitting;
	private Choice textFitSolver;
	private TextField textMatchDistance;
	private TextField textLowerDistance;
	private TextField textSignalFactor;
	private TextField textLowerFactor;
	private Checkbox cbSmartFilter;
	private TextField textCoordinateShiftFactor;
	private TextField textSignalStrength;
	private TextField textMinPhotons;
	private TextField textPrecisionThreshold;
	private Choice textPrecisionMethod;
	private TextField textMinWidthFactor;
	private TextField textWidthFactor;

	private AtomicInteger ignored = new AtomicInteger(0);

	private static final String[] NAMES = new String[] { "Candidate:N results in candidate",
			"Assigned Result:N results in assigned spot", "Singles:Neighbours", "Doublets:Neighbours",
			"Multiples:Neighbours", "Singles:Almost", "Doublets:Almost", "Multiples:Almost"

	};
	private static final String[] NAMES2 = { "Score n=1", "Score n=2", "Score n=N", "Iter n=1", "Eval n=1", "Iter n>1",
			"Eval n>1" };

	private static boolean[] displayHistograms = new boolean[NAMES.length + NAMES2.length];
	static
	{
		for (int i = 0; i < displayHistograms.length; i++)
			displayHistograms[i] = true;
	}

	private WindowOrganiser windowOrganiser = new WindowOrganiser();

	private class ResidualsScore
	{
		final double[] residuals, jaccard, recall, precision;
		final int maxJaccardIndex;
		double[] bestResiduals = new double[3];

		public ResidualsScore(double[] residuals, double[] jaccard, double[] recall, double[] precision,
				int maxJaccardIndex)
		{
			this.residuals = residuals;
			this.jaccard = jaccard;
			this.recall = recall;
			this.precision = precision;
			this.maxJaccardIndex = maxJaccardIndex;
		}
	}

	/**
	 * Allows plotting the bonus from fitting all spots at a given residuals threshold
	 */
	public class DoubletBonus implements Comparable<DoubletBonus>
	{
		final double rMax, rAv;
		final double tp, fp;
		public double residuals;

		/**
		 * Instantiates a new doublet bonus.
		 *
		 * @param rMax
		 *            the maximum residuals score
		 * @param rAv
		 *            the average residuals score
		 * @param tp
		 *            the additional true positives if this was accepted as a doublet
		 * @param fp
		 *            the additional false positives if this was accepted as a doublet
		 */
		public DoubletBonus(double rMax, double rAv, double tp, double fp)
		{
			this.rMax = rMax;
			this.rAv = rAv;
			this.residuals = rMax;
			this.tp = tp;
			this.fp = fp;
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(DoubletBonus that)
		{
			return Double.compare(this.residuals, that.residuals);
		}

		public void setScore(boolean useMax)
		{
			this.residuals = (useMax) ? rMax : rAv;
		}
	}

	private class ResultCoordinate extends BasePoint implements Comparable<ResultCoordinate>
	{
		final DoubletResult result;
		final int id;

		public ResultCoordinate(DoubletResult result, int id, double x, double y)
		{
			// Add the 0.5 pixel offset
			super((float) (x + 0.5), (float) (y + 0.5));
			this.result = result;
			this.id = id;
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(ResultCoordinate that)
		{
			//if (this.result.spot.intensity > that.result.spot.intensity)
			//	return -1;
			//if (this.result.spot.intensity < that.result.spot.intensity)
			//	return 1;
			//return 0;
			return Integer.compare(this.result.spotIndex, that.result.spotIndex);
		}
	}

	/**
	 * Stores results from single and doublet fitting.
	 */
	public class DoubletResult implements Comparable<DoubletResult>
	{
		final int frame;
		final float noise;
		final Spot spot;
		final int n, c, neighbours, almostNeighbours, spotIndex;
		FitResult fitResult1 = null;
		FitResult fitResult2 = null;
		double sumOfSquares1, sumOfSquares2;
		double ll1, ll2;
		double r1, r2;
		double value1, value2;
		double score1, score2;
		double aic1, aic2, bic1, bic2;
		double maic1, maic2, mbic1, mbic2;
		double[] xshift = new double[2];
		double[] yshift = new double[2];
		double[] a = new double[2];
		double gap;
		int iter1, iter2, eval1, eval2;
		boolean good1, good2, valid, valid2;
		double tp1, fp1, tp2a, fp2a, tp2b, fp2b;

		public DoubletResult(int frame, float noise, Spot spot, int n, int neighbours, int almostNeighbours,
				int spotIndex)
		{
			this.frame = frame;
			this.noise = noise;
			this.spot = spot;
			this.n = n;
			this.c = DoubletAnalysis.getClass(n);
			this.neighbours = neighbours;
			this.almostNeighbours = almostNeighbours;
			this.spotIndex = spotIndex;
			this.tp1 = 0;
			this.fp1 = 1;
			this.tp2a = 0;
			this.fp2a = 1;
			this.tp2b = 0;
			this.fp2b = 1;
		}

		public void addTP1(double score)
		{
			//if (score > 1)
			//	System.out.printf("Bad score %f\n", score);
			if (tp1 != 0)
				System.out.printf("Double counting: %f\n", score);
			tp1 += score;
			fp1 -= score;
		}

		public void addTP2(double score, int id)
		{
			//if (score > 1)
			//	System.out.printf("Bad score %f\n", score);
			if (id == 0)
			{
				tp2a += score;
				fp2a -= score;
			}
			else
			{
				tp2b += score;
				fp2b -= score;
			}
		}

		public double getMaxScore()
		{
			return (score1 > score2) ? score1 : score2;
		}

		public double getAvScore()
		{
			return (score1 + score2) * 0.5;
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(DoubletResult that)
		{
			// Makes the resutls easy to find in the table
			int r = this.frame - that.frame;
			if (r != 0)
				return r;
			r = this.spot.x - that.spot.x;
			if (r != 0)
				return r;
			return this.spot.y - that.spot.y;
			//if (this.spot.intensity > that.spot.intensity)
			//	return -1;
			//if (this.spot.intensity < that.spot.intensity)
			//	return 1;
			//return 0;
		}
	}

	/**
	 * Used to allow multi-threading of the fitting method.
	 */
	private class Worker implements Runnable
	{
		volatile boolean finished = false;
		final BlockingQueue<Integer> jobs;
		final ImageStack stack;
		final TIntObjectHashMap<ArrayList<Coordinate>> actualCoordinates;
		final int fitting;
		final FitConfiguration fitConfig;
		final MaximaSpotFilter spotFilter;
		final Gaussian2DFitter gf;
		final boolean relativeIntensity;
		final double limit;
		final int[] spotHistogram, resultHistogram;
		final int[][] neighbourHistogram;
		final int[][] almostNeighbourHistogram;
		final Overlay o;
		double[] region = null;
		float[] data = null;
		ArrayList<DoubletResult> results = new ArrayList<>();
		int daic = 0, dbic = 0, cic = 0;
		RampedScore rampedScore, signalScore = null;

		/**
		 * Instantiates a new worker.
		 *
		 * @param jobs
		 *            the jobs
		 * @param stack
		 *            the stack
		 * @param actualCoordinates
		 *            the actual coordinates
		 * @param fitConfig
		 *            the fit config
		 * @param o
		 *            the overlay
		 */
		public Worker(BlockingQueue<Integer> jobs, ImageStack stack,
				TIntObjectHashMap<ArrayList<Coordinate>> actualCoordinates, FitConfiguration fitConfig, Overlay o)
		{
			this.jobs = jobs;
			this.stack = stack;
			this.actualCoordinates = actualCoordinates;
			this.fitConfig = fitConfig.clone();
			this.gf = new Gaussian2DFitter(this.fitConfig);
			this.spotFilter = config.createSpotFilter();
			this.relativeIntensity = !spotFilter.isAbsoluteIntensity();

			fitting = config.getFittingWidth();
			// Fit window is 2*fitting+1. The distance limit is thus 0.5 pixel higher than fitting.
			limit = fitting + 0.5;
			spotHistogram = new int[20];
			resultHistogram = new int[spotHistogram.length];
			neighbourHistogram = new int[3][spotHistogram.length];
			almostNeighbourHistogram = new int[3][spotHistogram.length];
			this.o = o;
			rampedScore = new RampedScore(lowerDistance, matchDistance);
			if (signalFactor > 0)
				signalScore = new RampedScore(lowerSignalFactor, signalFactor);
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see java.lang.Runnable#run()
		 */
		@Override
		public void run()
		{
			try
			{
				while (true)
				{
					Integer job = jobs.take();
					if (job == null || job.intValue() < 0)
						break;
					if (!finished)
						// Only run if not finished to allow queue to be emptied
						run(job.intValue());
				}
			}
			catch (InterruptedException e)
			{
				System.out.println(e.toString());
				//throw new RuntimeException(e);
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

			Coordinate[] actual = ResultsMatchCalculator.getCoordinates(actualCoordinates, frame);
			ArrayList<DoubletResult> frameResults = new ArrayList<>(actual.length);

			// Extract the data
			final int maxx = stack.getWidth();
			final int maxy = stack.getHeight();
			data = IJImageConverter.getData(stack.getPixels(frame), maxx, maxy, null, data);

			// Smooth the image and identify spots with a filter
			Spot[] spots = spotFilter.rank(data, maxx, maxy);

			// Match the each actual result to the closest filter candidate.
			// The match must be within the fit window used during fitting, i.e. could the actual
			// result be fit using this candidate.
			int[] matches = new int[actual.length];
			for (int i = 0; i < actual.length; i++)
			{
				double dmin = Double.POSITIVE_INFINITY;
				int match = -1;
				// Get the coordinates, offset to allow for 0.5 to be the centre of the pixel
				double x = actual[i].getX() - 0.5;
				double y = actual[i].getY() - 0.5;
				for (int j = 0; j < spots.length; j++)
				{
					double dx = Math.abs(x - spots[j].x);
					double dy = Math.abs(y - spots[j].y);
					if (dx < limit && dy < limit)
					{
						final double d2 = dx * dx + dy * dy;
						if (dmin > d2)
						{
							dmin = d2;
							match = j;
						}
					}
				}
				matches[i] = match;
			}

			ImageExtractor ie = null;
			float estimatedBackground = 0;
			float noise = 0;

			// Identify single and doublets (and other)
			int singles = 0, doublets = 0, multiples = 0, total = 0, ignored = 0;
			int[] spotMatchCount = new int[spots.length];
			int[] neighbourIndices = new int[spots.length];
			for (int i = 0; i < actual.length; i++)
			{
				if (matches[i] != -1)
				{
					// Count all matches
					int n = 0;
					final int j = matches[i];
					for (int ii = i; ii < matches.length; ii++)
					{
						if (matches[ii] == j)
						{
							n++;
							// Reset to avoid double counting
							matches[ii] = -1;
						}
					}
					switch (n)
					{
						//@formatter:off
						case 1: singles++; break;
						case 2: doublets++; break;
						default: multiples++;
						//@formatter:on
					}

					// Initialise for fitting on first match
					if (ie == null)
					{
						ie = new ImageExtractor(data, maxx, maxy);
						estimatedBackground = estimateBackground(maxx, maxy);
						noise = FitWorker.estimateNoise(data, maxx, maxy, config.getNoiseMethod());
					}

					final Spot spot = spots[j];
					final Rectangle regionBounds = ie.getBoxRegionBounds(spot.x, spot.y, fitting);

					// Count the number of candidates within the fitting window
					// that are potential neighbours,
					// i.e. will fit neighbours be used?
					// It does not matter if the neighbours have a match to a result
					// or not, just that they are present for multiple peak fitting

					final int xmin = regionBounds.x;
					final int xmax = xmin + regionBounds.width - 1;
					final int ymin = regionBounds.y;
					final int ymax = ymin + regionBounds.height - 1;
					final int xmin2 = xmin - fitting;
					final int xmax2 = xmax + fitting;
					final int ymin2 = ymin - fitting;
					final int ymax2 = ymax + fitting;

					final float heightThreshold;
					float background = estimatedBackground;

					if (spot.intensity < background)
						heightThreshold = spot.intensity;
					else
						heightThreshold = (float) ((spot.intensity - background) *
								config.getNeighbourHeightThreshold() + background);

					int neighbourCount = 0;
					int almostNeighbourCount = 0;
					for (int jj = 0; jj < spots.length; jj++)
					{
						if (j == jj)
							continue;
						if (spots[jj].x < xmin2 || spots[jj].x > xmax2 || spots[jj].y < ymin2 || spots[jj].y > ymax2)
							continue;
						if (spots[jj].x < xmin || spots[jj].x > xmax || spots[jj].y < ymin || spots[jj].y > ymax ||
								spots[jj].intensity < heightThreshold)
						{
							almostNeighbourCount++;
						}
						else
						{
							neighbourIndices[neighbourCount++] = jj;
						}
					}

					// Optionally ignore results with neighbours
					if (ignoreWithNeighbours && neighbourCount != 0)
					{
						// TODO - Should we now remove these actual results from the actual[] array.
						// This removes them from any scoring of fit results.
						// For now leave them in as all spot candidates around them are going to be ignored
						// so there should not be any fit results close by.
						ignored += n;
						continue;
					}

					// Store the number of actual results that match to a spot
					addToHistogram(spotHistogram, n, 1);
					// Store the number of results that match to a spot with n results
					addToHistogram(resultHistogram, n, n);
					total += n;
					spotMatchCount[j] = n;

					final int c = DoubletAnalysis.getClass(spotMatchCount[j]);
					addToHistogram(neighbourHistogram[c], neighbourCount, 1);
					addToHistogram(almostNeighbourHistogram[c], almostNeighbourCount, 1);

					// TODO - Fit singles with neighbours?

					// Currently this will only explore how to benchmark the fitting of medium
					// density data with singles and a few doublets with no neighbours. This
					// is fine for low density PALM data but not for high density STORM data.

					// Fit the candidates (as per the FitWorker logic)
					// (Fit even multiple since this is what the FitWorker will do)
					region = ie.crop(regionBounds, region);

					boolean[] amplitudeEstimate = new boolean[1];
					float signal = 0;
					double sum = 0;
					final int width = regionBounds.width;
					final int height = regionBounds.height;
					final int size = width * height;
					for (int k = size; k-- > 0;)
						sum += region[k];
					signal = (float) (sum - background * size);
					if (signal <= 0)
					{
						amplitudeEstimate[0] = true;
						signal = spot.intensity - ((relativeIntensity) ? 0 : background);
						if (signal < 0)
						{
							signal += background;
							background = 0;
						}
					}

					final double[] params = new double[] { background, signal, 0, spot.x - regionBounds.x,
							spot.y - regionBounds.y, 0, 0 };

					final DoubletResult result = new DoubletResult(frame, noise, spot, n, neighbourCount,
							almostNeighbourCount, j);
					result.fitResult1 = gf.fit(region, width, height, 1, params, amplitudeEstimate);
					FunctionSolver f1 = gf.getFunctionSolver();
					result.iter1 = f1.getIterations();
					result.eval1 = f1.getEvaluations();

					// For now only process downstream if the fit was reasonable. This allows a good attempt at doublet fitting.
					result.good1 = goodFit(result.fitResult1, width, height) == 2;

					if (result.good1)
					{
						result.sumOfSquares1 = (f1.getType() == FunctionSolverType.LSE)
								? ((LSEFunctionSolver) f1).getTotalSumOfSquares() : 0;
						result.ll1 = (f1.getType() == FunctionSolverType.MLE)
								? ((MLEFunctionSolver) f1).getLogLikelihood() : 0;
						result.value1 = gf.getValue();

						// Compute residuals and fit as a doublet
						final double[] fitParams = result.fitResult1.getParameters();
						final int cx = (int) Math.round(fitParams[Gaussian2DFunction.X_POSITION]);
						final int cy = (int) Math.round(fitParams[Gaussian2DFunction.Y_POSITION]);
						final double[] residuals = gf.getResiduals();

						QuadrantAnalysis qa = new QuadrantAnalysis();

						// TODO - Also perform quadrant analysis on a new region centred around
						// the fit centre...

						if (qa.quadrantAnalysis(residuals, width, height, cx, cy) && qa.computeDoubletCentres(width,
								height, cx, cy, fitParams[Gaussian2DFunction.X_SD], fitParams[Gaussian2DFunction.Y_SD]))
						{
							result.score1 = qa.score1;
							result.score2 = qa.score2;

							// -+-+-
							// Estimate params using the single fitted peak
							// -+-+-
							final double[] doubletParams = new double[1 + 2 * Gaussian2DFunction.PARAMETERS_PER_PEAK];

							doubletParams[Gaussian2DFunction.BACKGROUND] = fitParams[Gaussian2DFunction.BACKGROUND];
							doubletParams[Gaussian2DFunction.SIGNAL] = fitParams[Gaussian2DFunction.SIGNAL] * 0.5;
							doubletParams[Gaussian2DFunction.X_POSITION] = (float) (qa.x1 - 0.5);
							doubletParams[Gaussian2DFunction.Y_POSITION] = (float) (qa.y1 - 0.5);
							doubletParams[Gaussian2DFunction.PARAMETERS_PER_PEAK +
									Gaussian2DFunction.SIGNAL] = params[Gaussian2DFunction.SIGNAL] * 0.5;
							doubletParams[Gaussian2DFunction.PARAMETERS_PER_PEAK +
									Gaussian2DFunction.X_POSITION] = (float) (qa.x2 - 0.5);
							doubletParams[Gaussian2DFunction.PARAMETERS_PER_PEAK +
									Gaussian2DFunction.Y_POSITION] = (float) (qa.y2 - 0.5);
							// -+-+-

							// Increase the iterations level then reset afterwards.
							final int maxIterations = fitConfig.getMaxIterations();
							final int maxEvaluations = fitConfig.getMaxFunctionEvaluations();
							fitConfig.setMaxIterations((int) (maxIterations *
									FitWorker.ITERATION_INCREASE_FOR_DOUBLETS * iterationIncrease));
							fitConfig.setMaxFunctionEvaluations((int) (maxEvaluations *
									FitWorker.EVALUATION_INCREASE_FOR_DOUBLETS * iterationIncrease));
							gf.setComputeResiduals(false);
							result.fitResult2 = gf.fit(region, width, height, 2, doubletParams, new boolean[2]);
							gf.setComputeResiduals(true);
							fitConfig.setMaxIterations(maxIterations);
							fitConfig.setMaxFunctionEvaluations(maxEvaluations);
							FunctionSolver f2 = gf.getFunctionSolver();
							result.iter2 = f2.getIterations();
							result.eval2 = f2.getEvaluations();
							int r2 = goodFit2(result.fitResult2, width, height);

							// Store all results if we made a fit, even if the fit was not good
							if (r2 != 0)
							{
								result.good2 = r2 == 2;
								result.sumOfSquares2 = (f2.getType() == FunctionSolverType.LSE)
										? ((LSEFunctionSolver) f2).getTotalSumOfSquares() : 0;
								result.ll2 = (f2.getType() == FunctionSolverType.MLE)
										? ((MLEFunctionSolver) f2).getLogLikelihood() : 0;
								result.value2 = gf.getValue();

								final int length = width * height;
								result.aic1 = Maths.getAkaikeInformationCriterionFromResiduals(result.sumOfSquares1,
										length, result.fitResult1.getNumberOfFittedParameters());
								result.aic2 = Maths.getAkaikeInformationCriterionFromResiduals(result.sumOfSquares2,
										length, result.fitResult2.getNumberOfFittedParameters());
								result.bic1 = Maths.getBayesianInformationCriterionFromResiduals(result.sumOfSquares1,
										length, result.fitResult1.getNumberOfFittedParameters());
								result.bic2 = Maths.getBayesianInformationCriterionFromResiduals(result.sumOfSquares2,
										length, result.fitResult2.getNumberOfFittedParameters());
								if (f2.getType() == FunctionSolverType.MLE)
								{
									result.maic1 = Maths.getAkaikeInformationCriterion(result.ll1, length,
											result.fitResult1.getNumberOfFittedParameters());
									result.maic2 = Maths.getAkaikeInformationCriterion(result.ll2, length,
											result.fitResult2.getNumberOfFittedParameters());
									result.mbic1 = Maths.getBayesianInformationCriterion(result.ll1, length,
											result.fitResult1.getNumberOfFittedParameters());
									result.mbic2 = Maths.getBayesianInformationCriterion(result.ll2, length,
											result.fitResult2.getNumberOfFittedParameters());

									// XXX - Debugging: see if the IC computed from the residuals would make a different choice
									// Disable by setting to 1
									if (result.getMaxScore() > 1)
									{
										cic++;
										if (Math.signum(result.aic1 - result.aic2) != Math
												.signum(result.maic1 - result.maic2))
										{
											daic++;
											System.out.printf(
													"AIC difference with residuals [%d] %d,%d : %d  %f vs %f (%.2f)\n",
													frame, spot.x, spot.y, n, Math.signum(result.aic1 - result.aic2),
													Math.signum(result.maic1 - result.maic2), result.getMaxScore());
										}
										if (Math.signum(result.bic1 - result.bic2) != Math
												.signum(result.mbic1 - result.mbic2))
										{
											dbic++;
											System.out.printf(
													"BIC difference with residuals [%d] %d,%d : %d  %f vs %f (%.2f)\n",
													frame, spot.x, spot.y, n, Math.signum(result.bic1 - result.bic2),
													Math.signum(result.mbic1 - result.mbic2), result.getMaxScore());
										}
										if (Double.isInfinite(result.value1) || Double.isInfinite(result.value2))
											System.out.printf("oops\n", result.value1, result.value2);
									}
								}
								else
								{
									result.maic1 = result.aic1;
									result.maic2 = result.aic2;
									result.mbic1 = result.bic1;
									result.mbic2 = result.bic2;
								}
								if (f1.getType() == FunctionSolverType.LSE)
								{
									result.r1 = ((LSEFunctionSolver) f1).getAdjustedCoefficientOfDetermination();
									result.r2 = ((LSEFunctionSolver) f2).getAdjustedCoefficientOfDetermination();
								}

								// Debugging: see if the AIC or BIC ever differ
								//if (Math.signum(result.aic1 - result.aic2) != Math.signum(result.bic1 - result.bic2))
								//	System.out.printf("BIC difference [%d] %d,%d : %d  %f vs %f (%.2f)\n", frame,
								//			spot.x, spot.y, n, Math.signum(result.aic1 - result.aic2),
								//			Math.signum(result.bic1 - result.bic2), result.getMaxScore());

								final double[] newParams = result.fitResult2.getParameters();
								for (int p = 0; p < 2; p++)
								{
									final double xShift = newParams[Gaussian2DFunction.X_POSITION +
											p * Gaussian2DFunction.PARAMETERS_PER_PEAK] -
											params[Gaussian2DFunction.X_POSITION];
									final double yShift = newParams[Gaussian2DFunction.Y_POSITION +
											p * Gaussian2DFunction.PARAMETERS_PER_PEAK] -
											params[Gaussian2DFunction.Y_POSITION];
									result.a[p] = 57.29577951 *
											QuadrantAnalysis.getAngle(qa.vector, new double[] { xShift, yShift });
									result.xshift[p] = xShift / limit;
									result.yshift[p] = yShift / limit;
								}

								// Store the distance between the spots
								final double dx = newParams[Gaussian2DFunction.X_POSITION] -
										newParams[Gaussian2DFunction.PARAMETERS_PER_PEAK +
												Gaussian2DFunction.X_POSITION];
								final double dy = newParams[Gaussian2DFunction.Y_POSITION] -
										newParams[Gaussian2DFunction.PARAMETERS_PER_PEAK +
												Gaussian2DFunction.Y_POSITION];
								result.gap = Math.sqrt(dx * dx + dy * dy);
							}
						}
					}

					// True results, i.e. where there was a choice between selecting fit results of single or doublet
					if (result.good1 && result.good2)
					{
						if (result.neighbours == 0)
						{
							result.valid = true;
							if (result.almostNeighbours == 0)
								result.valid2 = true;
						}
					}

					frameResults.add(result);
				}
			}

			DoubletAnalysis.this.ignored.addAndGet(ignored);

			//System.out.printf("Frame %d, singles=%d, doublets=%d, multi=%d\n", frame, singles, doublets, multiples);
			resultHistogram[0] += actual.length - total - ignored;

			addToOverlay(frame, spots, singles, doublets, multiples, spotMatchCount);

			results.addAll(frameResults);

			// At the end of all the fitting, assign results as true or false positive.

			if (matching == 0)
			{
				// Simple matching based on closest distance.
				// This is valid for comparing the score between residuals=1 (all single) and
				// residuals=0 (all doublets), when all spot candidates are fit.

				ArrayList<ResultCoordinate> f1 = new ArrayList<>();
				ArrayList<ResultCoordinate> f2 = new ArrayList<>();
				for (DoubletResult result : frameResults)
				{
					if (result.good1)
					{
						final Rectangle regionBounds = ie.getBoxRegionBounds(result.spot.x, result.spot.y, fitting);
						double x = result.fitResult1.getParameters()[Gaussian2DFunction.X_POSITION] + regionBounds.x;
						double y = result.fitResult1.getParameters()[Gaussian2DFunction.Y_POSITION] + regionBounds.y;
						f1.add(new ResultCoordinate(result, -1, x, y));
						if (result.good2)
						{
							double x2 = result.fitResult2.getParameters()[Gaussian2DFunction.X_POSITION] +
									regionBounds.x;
							double y2 = result.fitResult2.getParameters()[Gaussian2DFunction.Y_POSITION] +
									regionBounds.y;
							f2.add(new ResultCoordinate(result, 0, x2, y2));
							x2 = result.fitResult2.getParameters()[Gaussian2DFunction.X_POSITION +
									Gaussian2DFunction.PARAMETERS_PER_PEAK] + regionBounds.x;
							y2 = result.fitResult2.getParameters()[Gaussian2DFunction.Y_POSITION +
									Gaussian2DFunction.PARAMETERS_PER_PEAK] + regionBounds.y;
							f2.add(new ResultCoordinate(result, 1, x2, y2));
						}
					}
				}

				if (f1.isEmpty())
					return;

				List<PointPair> pairs = new ArrayList<>();
				MatchCalculator.analyseResults2D(actual, f1.toArray(new ResultCoordinate[f1.size()]), matchDistance,
						null, null, null, pairs);
				for (PointPair pair : pairs)
				{
					ResultCoordinate coord = (ResultCoordinate) pair.getPoint2();
					coord.result.addTP1(getScore(pair.getXYDistance2(), coord, pair.getPoint1()));
				}

				if (f2.isEmpty())
					return;

				// Note: Computing the closest match to all doublets
				// results in a simple analysis since we are comparing all doublets against all singles.
				// There may be a case where a actual coordinate is matched by a bad doublet but also by a
				// good doublet at a higher distance. However when we filter the results to remove bad doublets
				// (low residuals, etc) this result is not scored even though it would also match the better doublet.

				// This may not matter unless the density is high.

				MatchCalculator.analyseResults2D(actual, f2.toArray(new Coordinate[f2.size()]), matchDistance, null,
						null, null, pairs);
				for (PointPair pair : pairs)
				{
					ResultCoordinate coord = (ResultCoordinate) pair.getPoint2();
					coord.result.addTP2(getScore(pair.getXYDistance2(), coord, pair.getPoint1()), coord.id);
				}
			}
			else if (matching == 1)
			{
				// Rank doublets by the residuals score.
				// This is valid for comparing the score between residuals=1 (all singles)
				// and the effect of altering the residuals threshold to allow more doublets.
				// It is not a true effect as doublets with a higher residuals score may not be
				// first spot candidates that are fit.

				// Rank singles by the Candidate spot
				Collections.sort(frameResults, new Comparator<DoubletResult>()
				{
					@Override
					public int compare(DoubletResult o1, DoubletResult o2)
					{
						return o1.spotIndex - o2.spotIndex;
					}
				});

				final double threshold = matchDistance * matchDistance;
				final boolean[] assigned = new boolean[actual.length];

				int count = 0;
				int[] matched = new int[frameResults.size()];
				OUTER: for (int j = 0; j < frameResults.size(); j++)
				{
					DoubletResult result = frameResults.get(j);
					if (result.good1)
					{
						final Rectangle regionBounds = ie.getBoxRegionBounds(result.spot.x, result.spot.y, fitting);
						float x = (float) (result.fitResult1.getParameters()[Gaussian2DFunction.X_POSITION] +
								regionBounds.x);
						float y = (float) (result.fitResult1.getParameters()[Gaussian2DFunction.Y_POSITION] +
								regionBounds.y);

						for (int i = 0; i < actual.length; i++)
						{
							if (assigned[i])
								continue;
							final double d2 = actual[i].distance2(x, y);
							if (d2 <= threshold)
							{
								assigned[i] = true;
								matched[j] = i + 1;
								result.addTP1(getScore(d2, result, -1, actual[i]));
								if (++count == actual.length)
									break OUTER;
								break;
							}
						}
					}
				}

				// Rank the doublets by residuals threshold instead. 1 from the doublet
				// must match the spot that it matched as a single (if still available).
				// The other can match anything else...
				Collections.sort(frameResults, new Comparator<DoubletResult>()
				{
					@Override
					public int compare(DoubletResult o1, DoubletResult o2)
					{
						if (o1.getMaxScore() > o2.getMaxScore())
							return -1;
						if (o1.getMaxScore() < o2.getMaxScore())
							return -1;
						return 0;
					}
				});

				count = 0;
				Arrays.fill(assigned, false);
				OUTER: for (int j = 0; j < frameResults.size(); j++)
				{
					DoubletResult result = frameResults.get(j);
					if (result.good1)
					{
						final Rectangle regionBounds = ie.getBoxRegionBounds(result.spot.x, result.spot.y, fitting);

						if (result.good2)
						{
							float x1 = (float) (result.fitResult2.getParameters()[Gaussian2DFunction.X_POSITION] +
									regionBounds.x);
							float y1 = (float) (result.fitResult2.getParameters()[Gaussian2DFunction.Y_POSITION] +
									regionBounds.y);
							float x2 = (float) (result.fitResult2.getParameters()[Gaussian2DFunction.X_POSITION +
									Gaussian2DFunction.PARAMETERS_PER_PEAK] + regionBounds.x);
							float y2 = (float) (result.fitResult2.getParameters()[Gaussian2DFunction.Y_POSITION +
									Gaussian2DFunction.PARAMETERS_PER_PEAK] + regionBounds.y);

							ResultCoordinate ra = new ResultCoordinate(result, 0, x1, y1);
							ResultCoordinate rb = new ResultCoordinate(result, 1, x2, y2);

							// Q. what did the single match?
							int i = matched[j] - 1;
							if (i != -1 && !assigned[i])
							{
								// One of the doublet pair must match this
								final double d2a = ra.distanceXY2(actual[i]);
								final double d2b = rb.distanceXY2(actual[i]);
								if (d2a < d2b)
								{
									if (d2a <= threshold)
									{
										ra.result.addTP2(getScore(d2a, ra, actual[i]), ra.id);
										if (++count == actual.length)
											break OUTER;
										assigned[i] = true;
										ra = null;
									}
								}
								else
								{
									if (d2b <= threshold)
									{
										rb.result.addTP2(getScore(d2b, rb, actual[i]), rb.id);
										if (++count == actual.length)
											break OUTER;
										assigned[i] = true;
										rb = null;
									}
								}
							}

							// Process the rest of the results
							for (i = 0; i < actual.length; i++)
							{
								if (assigned[i])
									continue;
								if (ra == null)
								{
									final double d2 = rb.distanceXY2(actual[i]);
									if (d2 <= threshold)
									{
										rb.result.addTP2(getScore(d2, rb, actual[i]), rb.id);
										if (++count == actual.length)
											break OUTER;
										assigned[i] = true;
										break;
									}
								}
								else if (rb == null)
								{
									final double d2 = ra.distanceXY2(actual[i]);
									if (d2 <= threshold)
									{
										ra.result.addTP2(getScore(d2, ra, actual[i]), ra.id);
										if (++count == actual.length)
											break OUTER;
										assigned[i] = true;
										break;
									}
								}
								else
								{
									final double d2a = ra.distanceXY2(actual[i]);
									final double d2b = rb.distanceXY2(actual[i]);
									if (d2a < d2b)
									{
										if (d2a <= threshold)
										{
											ra.result.addTP2(getScore(d2a, ra, actual[i]), ra.id);
											if (++count == actual.length)
												break OUTER;
											assigned[i] = true;
											ra = null;
										}
									}
									else
									{
										if (d2b <= threshold)
										{
											rb.result.addTP2(getScore(d2b, rb, actual[i]), rb.id);
											if (++count == actual.length)
												break OUTER;
											assigned[i] = true;
											rb = null;
										}
									}
								}
							}
						}
					}
				}
			}
			else
			{
				// Matching based on spot ranking
				ArrayList<ResultCoordinate> f1 = new ArrayList<>();
				ArrayList<ResultCoordinate> f2 = new ArrayList<>();
				for (DoubletResult result : frameResults)
				{
					if (result.good1)
					{
						final Rectangle regionBounds = ie.getBoxRegionBounds(result.spot.x, result.spot.y, fitting);
						double x = result.fitResult1.getParameters()[Gaussian2DFunction.X_POSITION] + regionBounds.x;
						double y = result.fitResult1.getParameters()[Gaussian2DFunction.Y_POSITION] + regionBounds.y;
						f1.add(new ResultCoordinate(result, -1, x, y));
						if (result.good2)
						{
							double x2 = result.fitResult2.getParameters()[Gaussian2DFunction.X_POSITION] +
									regionBounds.x;
							double y2 = result.fitResult2.getParameters()[Gaussian2DFunction.Y_POSITION] +
									regionBounds.y;
							f2.add(new ResultCoordinate(result, 0, x2, y2));
							x2 = result.fitResult2.getParameters()[Gaussian2DFunction.X_POSITION +
									Gaussian2DFunction.PARAMETERS_PER_PEAK] + regionBounds.x;
							y2 = result.fitResult2.getParameters()[Gaussian2DFunction.Y_POSITION +
									Gaussian2DFunction.PARAMETERS_PER_PEAK] + regionBounds.y;
							f2.add(new ResultCoordinate(result, 1, x2, y2));
						}
					}
				}

				if (f1.isEmpty())
					return;

				Collections.sort(f1);
				Collections.sort(f2);

				// Match the singles to actual coords, rank by the Candidate spot (not distance)
				final double threshold = matchDistance * matchDistance;
				int count = 0;
				boolean[] assigned = new boolean[actual.length];
				OUTER: for (ResultCoordinate r : f1)
				{
					for (int i = 0; i < actual.length; i++)
					{
						if (assigned[i])
							continue;
						final double d2 = r.distanceXY2(actual[i]);
						if (d2 <= threshold)
						{
							r.result.addTP1(getScore(d2, r, actual[i]));
							if (++count == actual.length)
								break OUTER;
							assigned[i] = true;
							break;
						}
					}
				}

				// Match the doublets to actual coords, rank by the Candidate spot (not distance)
				// Process in pairs
				count = 0;
				Arrays.fill(assigned, false);
				OUTER: for (int j = 0; j < f2.size(); j += 2)
				{
					ResultCoordinate ra = f2.get(j);
					ResultCoordinate rb = f2.get(j + 1);

					for (int i = 0; i < actual.length; i++)
					{
						if (assigned[i])
							continue;
						if (ra == null)
						{
							final double d2 = rb.distanceXY2(actual[i]);
							if (d2 <= threshold)
							{
								rb.result.addTP2(getScore(d2, rb, actual[i]), rb.id);
								if (++count == actual.length)
									break OUTER;
								assigned[i] = true;
								break;
							}
						}
						else if (rb == null)
						{
							final double d2 = ra.distanceXY2(actual[i]);
							if (d2 <= threshold)
							{
								ra.result.addTP2(getScore(d2, ra, actual[i]), ra.id);
								if (++count == actual.length)
									break OUTER;
								assigned[i] = true;
								break;
							}
						}
						else
						{
							final double d2a = ra.distanceXY2(actual[i]);
							final double d2b = rb.distanceXY2(actual[i]);
							if (d2a < d2b)
							{
								if (d2a <= threshold)
								{
									ra.result.addTP2(getScore(d2a, ra, actual[i]), ra.id);
									if (++count == actual.length)
										break OUTER;
									assigned[i] = true;
									ra = null;
								}
							}
							else
							{
								if (d2b <= threshold)
								{
									rb.result.addTP2(getScore(d2b, rb, actual[i]), rb.id);
									if (++count == actual.length)
										break OUTER;
									assigned[i] = true;
									rb = null;
								}
							}
						}
					}
				}
			}

		}

		private void addToHistogram(int[] h, int i, int n)
		{
			if (h.length <= i)
			{
				final int newLength = (int) Math.ceil(i * 1.5);
				h = Arrays.copyOf(h, newLength);
			}
			h[i] += n;
		}

		private double getScore(double d2, ResultCoordinate resultCoord, Coordinate coord)
		{
			return getScore(d2, resultCoord.result, resultCoord.id, coord);
		}

		private double getScore(double d2, DoubletResult result, int id, Coordinate coord)
		{
			double matchScore = rampedScore.scoreAndFlatten(d2, 256);
			if (signalScore != null)
			{
				PeakResultPoint p = (PeakResultPoint) coord;
				double s1 = p.peakResult.getIntensity();
				double s2;
				switch (id)
				{
					case -1:
						s2 = result.fitResult1.getParameters()[Gaussian2DFunction.SIGNAL];
						break;

					case 1:
						s2 = result.fitResult2.getParameters()[Gaussian2DFunction.SIGNAL +
								Gaussian2DFunction.PARAMETERS_PER_PEAK];
						break;

					case 0:
					default:
						s2 = result.fitResult2.getParameters()[Gaussian2DFunction.SIGNAL];
				}
				double rsf = s1 / s2;
				double sf = Math.abs((rsf < 1) ? 1 - 1 / rsf : rsf - 1);
				final double fScore = signalScore.scoreAndFlatten(sf, 256);
				matchScore = RampedScore.flatten(matchScore * fScore, 256);
			}
			return matchScore;
		}

		private int goodFit(FitResult fitResult, final int width, final int height)
		{
			if (fitResult == null)
				return 0;
			final double[] params = fitResult.getParameters();
			if (params == null)
				return 0;
			switch (fitResult.getStatus())
			{
				case OK:
					break;

				// The following happen when we are doing validation
				case INSUFFICIENT_SIGNAL:
					return 1; // Failed validation

				case WIDTH_DIVERGED:
				case INSUFFICIENT_PRECISION:
				case COORDINATES_MOVED:
					break; // Ignore these and check again

				default:
					return 0;
			}

			// Do some simple validation

			// Check if centre is within the region
			final double border = FastMath.min(width, height) / 4.0;
			if ((params[Gaussian2DFunction.X_POSITION] < border ||
					params[Gaussian2DFunction.X_POSITION] > width - border) ||
					params[Gaussian2DFunction.Y_POSITION] < border ||
					params[Gaussian2DFunction.Y_POSITION] > height - border)
				return 1;

			// Check the width is reasonable
			final double regionSize = FastMath.max(width, height) * 0.5;
			if (params[Gaussian2DFunction.X_SD] < 0 || params[Gaussian2DFunction.X_SD] > regionSize ||
					params[Gaussian2DFunction.Y_SD] < 0 || params[Gaussian2DFunction.Y_SD] > regionSize)
				return 1;
			return 2;
		}

		private int goodFit2(FitResult fitResult, final int width, final int height)
		{
			if (fitResult == null)
				return 0;
			final double[] params = fitResult.getParameters();
			if (params == null)
				return 0;
			switch (fitResult.getStatus())
			{
				case OK:
					break;

				// The following happen when we are doing validation
				case INSUFFICIENT_SIGNAL:
					return 1; // Failed validation

				case WIDTH_DIVERGED:
				case INSUFFICIENT_PRECISION:
				case COORDINATES_MOVED:
					break; // Ignore these and check again

				default:
					return 0;
			}

			// Do some simple validation

			final double regionSize = FastMath.max(width, height) * 0.5;
			for (int n = 0; n < 2; n++)
			{
				// Check the width is reasonable
				if (params[n * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_SD] < 0 ||
						params[n * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_SD] > regionSize ||
						params[n * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.Y_SD] < 0 ||
						params[n * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.Y_SD] > regionSize)
					return 1;

				// Check if centre is within the region - Border allowing fit slightly outside
				final double borderx = Gaussian2DFunction.SD_TO_HWHM_FACTOR * fitConfig.getInitialXSD();
				final double bordery = Gaussian2DFunction.SD_TO_HWHM_FACTOR * fitConfig.getInitialYSD();
				if ((params[n * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_POSITION] < -borderx ||
						params[n * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_POSITION] > width +
								borderx) ||
						params[n * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.Y_POSITION] < -bordery ||
						params[n * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.Y_POSITION] > height +
								bordery)
				{
					// Perhaps do a check on the quadrant?
					return 1;
				}
			}

			return 2;
		}

		/**
		 * Get an estimate of the background level using the mean of image.
		 *
		 * @param width
		 *            the width
		 * @param height
		 *            the height
		 * @return the float
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
		 * Adds the to overlay.
		 *
		 * @param frame
		 *            the frame
		 * @param spots
		 *            the spots
		 * @param singles
		 *            the singles
		 * @param doublets
		 *            the doublets
		 * @param multiples
		 *            the multiples
		 * @param spotMatchCount
		 *            the spot match count
		 */
		private void addToOverlay(int frame, Spot[] spots, int singles, int doublets, int multiples,
				int[] spotMatchCount)
		{
			if (o != null)
			{
				// Create an output stack with coloured ROI overlay for each n=1, n=2, n=other
				// to check that the doublets are correctly identified.
				final int[] sx = new int[singles];
				final int[] sy = new int[singles];
				final int[] dx = new int[doublets];
				final int[] dy = new int[doublets];
				final int[] mx = new int[multiples];
				final int[] my = new int[multiples];
				final int[] count = new int[3];
				final int[][] coords = new int[][] { sx, dx, mx, sy, dy, my };
				final Color[] color = new Color[] { Color.red, Color.green, Color.blue };
				for (int j = 0; j < spotMatchCount.length; j++)
				{
					final int c = DoubletAnalysis.getClass(spotMatchCount[j]);
					if (c < 0)
						continue;
					coords[c][count[c]] = spots[j].x;
					coords[c + 3][count[c]] = spots[j].y;
					count[c]++;
				}
				for (int c = 0; c < 3; c++)
				{
					final PointRoi roi = new PointRoi(coords[c], coords[c + 3], count[c]);
					roi.setPosition(frame);
					roi.setShowLabels(false);
					roi.setFillColor(color[c]);
					roi.setStrokeColor(color[c]);
					// Overlay uses a vector which is synchronized already
					o.add(roi);
				}
			}
		}

	}

	/**
	 * Gets the class.
	 *
	 * @param n
	 *            the number of results that match to the spot
	 * @return the class (none = -1, single = 0, double = 1, multiple = 2)
	 */
	private static int getClass(int n)
	{
		if (n == 0)
			return -1;
		if (n == 1)
			return 0;
		if (n == 2)
			return 1;
		return 2;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	@Override
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if ("analysis".equals(arg))
		{
			runAnalysis();
		}
		else
		{
			simulationParameters = CreateData.simulationParameters;
			if (simulationParameters == null)
			{
				IJ.error(TITLE, "No simulation parameters in memory");
				return;
			}
			imp = CreateData.getImage();
			if (imp == null)
			{
				IJ.error(TITLE, "No simulation image");
				return;
			}
			results = CreateData.getResults();
			if (results == null)
			{
				IJ.error(TITLE, "No simulation results in memory");
				return;
			}

			if (!showDialog())
				return;

			run();
		}
	}

	/**
	 * Show dialog.
	 *
	 * @return true, if successful
	 */
	private boolean showDialog()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		final double sa = getSa();
		gd.addMessage(
				String.format("Fits the benchmark image created by CreateData plugin.\nPSF width = %s, adjusted = %s",
						Utils.rounded(simulationParameters.s / simulationParameters.a), Utils.rounded(sa)));

		// For each new benchmark width, reset the PSF width to the square pixel adjustment
		if (lastId != simulationParameters.id)
		{
			double w = sa;
			matchDistance = w * Gaussian2DFunction.SD_TO_HWHM_FACTOR;
			lowerDistance = 0.5 * matchDistance;
			fitConfig.setInitialPeakStdDev(w);

			CalibrationWriter cal = new CalibrationWriter(fitConfig.getCalibration());

			cal.setNmPerPixel(simulationParameters.a);
			cal.setCountPerPhoton(simulationParameters.gain);
			cal.setQuantumEfficiency(simulationParameters.qe);
			cal.setExposureTime(100);
			cal.setReadNoise(simulationParameters.readNoise);
			cal.setBias(simulationParameters.bias);
			cal.setCameraType(simulationParameters.cameraType);

			fitConfig.setCalibration(cal.getCalibration());
		}

		// Support for using templates
		String[] templates = ConfigurationTemplate.getTemplateNames(true);
		gd.addChoice("Template", templates, templates[0]);

		// Allow the settings from the benchmark analysis to be used
		gd.addCheckbox("Benchmark_settings", useBenchmarkSettings);

		// Collect options for fitting
		PeakFit.addPSFOptions(gd, fitConfig);
		PeakFit.SimpleFitEngineConfigurationProvider provider = new PeakFit.SimpleFitEngineConfigurationProvider(
				config);
		PeakFit.addDataFilterOptions(gd, provider);
		PeakFit.addSearchOptions(gd, provider);
		PeakFit.addBorderOptions(gd, provider);
		PeakFit.addFittingOptions(gd, provider);

		gd.addChoice("Fit_solver", SettingsManager.getFitSolverNames(), fitConfig.getFitSolver().ordinal());

		gd.addSlider("Iteration_increase", 1, 4.5, iterationIncrease);
		gd.addCheckbox("Ignore_with_neighbours", ignoreWithNeighbours);
		gd.addCheckbox("Show_overlay", showOverlay);
		gd.addCheckbox("Show_histograms", showHistograms);
		gd.addCheckbox("Show_results", showResults);
		gd.addCheckbox("Show_Jaccard_Plot", showJaccardPlot);
		gd.addCheckbox("Use_max_residuals", useMaxResiduals);
		gd.addNumericField("Match_distance", matchDistance, 2);
		gd.addNumericField("Lower_distance", lowerDistance, 2);
		gd.addNumericField("Signal_factor", signalFactor, 2);
		gd.addNumericField("Lower_factor", lowerSignalFactor, 2);
		gd.addChoice("Matching", MATCHING, MATCHING[matching]);

		// Add a mouse listener to the config file field
		if (Utils.isShowGenericDialog())
		{
			Vector<TextField> numerics = gd.getNumericFields();
			Vector<Choice> choices = gd.getChoices();
			int n = 0;
			int ch = 0;

			choices.get(ch++).addItemListener(this);
			Checkbox b = (Checkbox) gd.getCheckboxes().get(0);
			b.addItemListener(this);
			textPSF = choices.get(ch++);
			textDataFilterType = choices.get(ch++);
			textDataFilter = choices.get(ch++);
			textSmooth = numerics.get(n++);
			textSearch = numerics.get(n++);
			textBorder = numerics.get(n++);
			textFitting = numerics.get(n++);
			textFitSolver = choices.get(ch++);
			n++; // Iteration increase
			textMatchDistance = numerics.get(n++);
			textLowerDistance = numerics.get(n++);
			textSignalFactor = numerics.get(n++);
			textLowerFactor = numerics.get(n++);
		}

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		// Ignore the template
		gd.getNextChoice();
		useBenchmarkSettings = gd.getNextBoolean();
		fitConfig.setPSFType(PeakFit.getPSFTypeValues()[gd.getNextChoiceIndex()]);
		config.setDataFilterType(gd.getNextChoiceIndex());
		config.setDataFilter(gd.getNextChoiceIndex(), Math.abs(gd.getNextNumber()), false, 0);
		config.setSearch(gd.getNextNumber());
		config.setBorder(gd.getNextNumber());
		config.setFitting(gd.getNextNumber());
		fitConfig.setFitSolver(gd.getNextChoiceIndex());

		// Avoid stupidness. Note: We are mostly ignoring the validation result and
		// checking the results for the doublets manually.
		fitConfig.setMinPhotons(15); // Realistically we cannot fit lower than this
		// Set the width factors to help establish bounds for bounded fitters
		fitConfig.setMinWidthFactor(1.0 / 10);
		fitConfig.setWidthFactor(10);

		iterationIncrease = gd.getNextNumber();
		ignoreWithNeighbours = gd.getNextBoolean();
		showOverlay = gd.getNextBoolean();
		showHistograms = gd.getNextBoolean();
		showResults = gd.getNextBoolean();
		showJaccardPlot = gd.getNextBoolean();
		useMaxResiduals = gd.getNextBoolean();
		matchDistance = Math.abs(gd.getNextNumber());
		lowerDistance = Math.abs(gd.getNextNumber());
		signalFactor = Math.abs(gd.getNextNumber());
		lowerSignalFactor = Math.abs(gd.getNextNumber());
		matching = gd.getNextChoiceIndex();

		gd.collectOptions();

		if (gd.invalidNumber())
			return false;

		if (lowerDistance > matchDistance)
			lowerDistance = matchDistance;
		if (lowerSignalFactor > signalFactor)
			lowerSignalFactor = signalFactor;

		if (useBenchmarkSettings)
		{
			if (!updateFitConfiguration(config))
				return false;
		}

		boolean configure = true;
		if (useBenchmarkSettings)
		{
			// Only configure the fit solver if not in a macro
			configure = Macro.getOptions() == null;
		}
		if (configure && !PeakFit.configurePSFModel(config))
			return false;
		if (configure && !PeakFit.configureFitSolver(config, IJImageSource.getBounds(imp), null, PeakFit.FLAG_NO_SAVE))
			return false;

		lastId = simulationParameters.id;

		if (showHistograms)
		{
			gd = new ExtendedGenericDialog(TITLE);
			gd.addMessage("Select the histograms to display");

			for (int i = 0; i < NAMES.length; i++)
				gd.addCheckbox(NAMES[i].replace(' ', '_'), displayHistograms[i]);
			for (int i = 0; i < NAMES2.length; i++)
				gd.addCheckbox(NAMES2[i].replace(' ', '_'), displayHistograms[i + NAMES.length]);
			gd.showDialog();
			if (gd.wasCanceled())
				return false;
			for (int i = 0; i < displayHistograms.length; i++)
				displayHistograms[i] = gd.getNextBoolean();
		}

		return true;
	}

	private boolean updateFitConfiguration(FitEngineConfiguration config)
	{
		// Do this first as it sets the initial SD
		if (!BenchmarkSpotFit.updateConfiguration(config))
		{
			IJ.error(TITLE, "Unable to use the benchmark spot fit configuration");
			return false;
		}

		CalibrationWriter cal = new CalibrationWriter(fitConfig.getCalibration());

		cal.setNmPerPixel(simulationParameters.a);
		cal.setCountPerPhoton(simulationParameters.gain);
		cal.setQuantumEfficiency(simulationParameters.qe);
		cal.setExposureTime(100);
		cal.setReadNoise(simulationParameters.readNoise);
		cal.setBias(simulationParameters.bias);
		cal.setCameraType(simulationParameters.cameraType);

		fitConfig.setCalibration(cal.getCalibration());

		if (!BenchmarkSpotFilter.updateConfiguration(config))
		{
			IJ.error(TITLE, "Unable to use the benchmark spot filter configuration");
			return false;
		}
		// Make sure all spots are fit
		config.setFailuresLimit(-1);

		// Get the distance from the filter analysis. This ensures that we compute scores the
		// same as the filter analysis
		if (BenchmarkFilterAnalysis.distanceInPixels > 0)
		{
			matchDistance = BenchmarkFilterAnalysis.distanceInPixels;
			lowerDistance = BenchmarkFilterAnalysis.lowerDistanceInPixels;
			signalFactor = BenchmarkFilterAnalysis.signalFactor;
			lowerSignalFactor = BenchmarkFilterAnalysis.lowerSignalFactor;
		}
		else
		{
			// Use the fit analysis distance if no filter analysis has been run
			matchDistance = BenchmarkSpotFit.distanceInPixels;
			lowerDistance = BenchmarkSpotFit.lowerDistanceInPixels;
			signalFactor = lowerSignalFactor = BenchmarkSpotFit.signalFactor;
		}

		return true;
	}

	/**
	 * Gets the sa.
	 *
	 * @return the sa
	 */
	private double getSa()
	{
		final double sa = PSFCalculator.squarePixelAdjustment(simulationParameters.s, simulationParameters.a) /
				simulationParameters.a;
		return sa;
	}

	/** The total progress. */
	int progress, stepProgress, totalProgress;

	/**
	 * Show progress.
	 */
	private synchronized void showProgress()
	{
		if (progress % stepProgress == 0)
		{
			if (Utils.showStatus("Frame: " + progress + " / " + totalProgress))
				IJ.showProgress(progress, totalProgress);
		}
		progress++;
	}

	/**
	 * Run.
	 */
	private void run()
	{
		doubletResults = null;

		final ImageStack stack = imp.getImageStack();

		// Get the coordinates per frame
		TIntObjectHashMap<ArrayList<Coordinate>> actualCoordinates = ResultsMatchCalculator.getCoordinates(results,
				false);

		final long[] sumCount = new long[1];
		actualCoordinates.forEachValue(new TObjectProcedure<ArrayList<Coordinate>>()
		{
			@Override
			public boolean execute(ArrayList<Coordinate> list)
			{
				sumCount[0] += list.size();
				return true;
			}
		});
		final double density = 1e6 * sumCount[0] / (simulationParameters.a * simulationParameters.a *
				results.getBounds().getWidth() * results.getBounds().getHeight() * actualCoordinates.size());

		// Create a pool of workers
		final int nThreads = Prefs.getThreads();
		final BlockingQueue<Integer> jobs = new ArrayBlockingQueue<>(nThreads * 2);
		List<Worker> workers = new LinkedList<>();
		List<Thread> threads = new LinkedList<>();
		Overlay overlay = (showOverlay) ? new Overlay() : null;
		for (int i = 0; i < nThreads; i++)
		{
			Worker worker = new Worker(jobs, stack, actualCoordinates, fitConfig, overlay);
			Thread t = new Thread(worker);
			workers.add(worker);
			threads.add(t);
			t.start();
		}

		// Fit the frames
		long runTime = System.nanoTime();
		totalProgress = actualCoordinates.size();
		stepProgress = Utils.getProgressInterval(totalProgress);
		progress = 0;
		actualCoordinates.forEachKey(new TIntProcedure()
		{
			@Override
			public boolean execute(int frame)
			{
				put(jobs, frame);
				return true;
			}
		});
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
		threads = null;

		IJ.showProgress(1);
		IJ.showStatus("Collecting results ...");
		runTime = System.nanoTime() - runTime;

		// Collect the results
		int cic = 0, daic = 0, dbic = 0;
		ArrayList<DoubletResult> results = null;
		int maxH = 0, maxH2 = 0, maxH3 = 0;
		for (Worker worker : workers)
		{
			if (results == null)
				results = worker.results;
			else
				results.addAll(worker.results);
			cic += worker.cic;
			daic += worker.daic;
			dbic += worker.dbic;
			maxH = Maths.max(maxH, worker.spotHistogram.length);
			for (int k = 0; k < 3; k++)
			{
				maxH2 = Maths.max(maxH2, worker.neighbourHistogram[k].length);
				maxH3 = Maths.max(maxH3, worker.almostNeighbourHistogram[k].length);
			}
		}
		if (cic > 0)
			System.out.printf("Difference AIC %d, BIC %d, Total %d\n", daic, dbic, cic);
		if (showHistograms)
		{
			double[] spotHistogram = new double[maxH];
			double[] resultHistogram = new double[maxH];
			double[][] neighbourHistogram = new double[3][maxH2];
			double[][] almostNeighbourHistogram = new double[3][maxH3];
			for (Worker worker : workers)
			{
				final int[] h1a = worker.spotHistogram;
				final int[] h1b = worker.resultHistogram;
				for (int j = 0; j < h1a.length; j++)
				{
					spotHistogram[j] += h1a[j];
					resultHistogram[j] += h1b[j];
				}
				final int[][] h2 = worker.neighbourHistogram;
				final int[][] h3 = worker.almostNeighbourHistogram;
				for (int k = 0; k < 3; k++)
				{
					for (int j = 0; j < h2[k].length; j++)
						neighbourHistogram[k][j] += h2[k][j];
					for (int j = 0; j < h3[k].length; j++)
						almostNeighbourHistogram[k][j] += h3[k][j];
				}
			}

			showHistogram(0, spotHistogram);
			showHistogram(1, resultHistogram);
			showHistogram(2, neighbourHistogram[0]);
			showHistogram(3, neighbourHistogram[1]);
			showHistogram(4, neighbourHistogram[2]);
			showHistogram(5, almostNeighbourHistogram[0]);
			showHistogram(6, almostNeighbourHistogram[1]);
			showHistogram(7, almostNeighbourHistogram[2]);
		}
		workers.clear();
		workers = null;

		if (overlay != null)
			imp.setOverlay(overlay);

		MemoryPeakResults.freeMemory();

		Collections.sort(results);
		summariseResults(results, density, runTime);

		windowOrganiser.tile();

		IJ.showStatus("");
	}

	/**
	 * Show histogram.
	 *
	 * @param i
	 *            the i
	 * @param histogram
	 *            the spot histogram
	 */
	private void showHistogram(int i, double[] histogram)
	{
		if (!displayHistograms[i])
			return;
		// Truncate to correct size
		for (int j = histogram.length; j-- > 0;)
			if (histogram[j] != 0)
			{
				histogram = Arrays.copyOf(histogram, j + 1);
				break;
			}
		String[] labels = NAMES[i].split(":");
		Plot2 plot = new Plot2(labels[0], labels[1], "Count");
		double max = Maths.max(histogram);
		plot.setLimits(0, histogram.length, 0, max * 1.05);
		plot.addPoints(SimpleArrayUtils.newArray(histogram.length, 0, 1.0), histogram, Plot2.BAR);
		PlotWindow pw = Utils.display(labels[0], plot);
		if (Utils.isNewWindow())
			windowOrganiser.add(pw.getImagePlus().getID());
	}

	/**
	 * Put.
	 *
	 * @param jobs
	 *            the jobs
	 * @param i
	 *            the i
	 */
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

	/**
	 * Summarise results.
	 *
	 * @param results
	 *            the results
	 * @param density
	 *            the density
	 * @param runTime
	 *            the run time
	 */
	private void summariseResults(ArrayList<DoubletResult> results, double density, long runTime)
	{
		// Store results in memory for later analysis
		doubletResults = results;

		// If we are only assessing results with no neighbour candidates
		// TODO - Count the number of actual results that have no neighbours

		numberOfMolecules = this.results.size() - ignored.get();

		// Store details we want in the analysis table
		StringBuilder sb = new StringBuilder();
		sb.append(Utils.rounded(density)).append('\t');
		sb.append(Utils.rounded(getSa())).append('\t');
		sb.append(config.getFittingWidth()).append('\t');
		sb.append(PSFProtosHelper.getName(fitConfig.getPSFType()));
		sb.append(":").append(PeakFit.getSolverName(fitConfig));
		if (fitConfig.isModelCameraMLE())
		{
			sb.append(":Camera\t");

			// Add details of the noise model for the MLE
			CalibrationReader r = new CalibrationReader(fitConfig.getCalibration());
			sb.append("EM=").append(r.isEMCCD());
			sb.append(":A=").append(Utils.rounded(r.getCountPerElectron()));
			sb.append(":N=").append(Utils.rounded(r.getReadNoise()));
			sb.append('\t');
		}
		else
			sb.append("\t\t");
		analysisPrefix = sb.toString();

		// -=-=-=-=-

		showResults(results, showResults);

		createSummaryTable();

		sb.setLength(0);

		final int n = countN(results);

		// Create the benchmark settings and the fitting settings
		sb.append(numberOfMolecules).append('\t');
		sb.append(n).append('\t');
		sb.append(Utils.rounded(density)).append('\t');
		sb.append(Utils.rounded(simulationParameters.minSignal)).append('\t');
		sb.append(Utils.rounded(simulationParameters.maxSignal)).append('\t');
		sb.append(Utils.rounded(simulationParameters.averageSignal)).append('\t');
		sb.append(Utils.rounded(simulationParameters.s)).append('\t');
		sb.append(Utils.rounded(simulationParameters.a)).append('\t');
		sb.append(Utils.rounded(getSa() * simulationParameters.a)).append('\t');
		sb.append(Utils.rounded(simulationParameters.gain)).append('\t');
		sb.append(Utils.rounded(simulationParameters.readNoise)).append('\t');
		sb.append(Utils.rounded(simulationParameters.b)).append('\t');
		sb.append(Utils.rounded(simulationParameters.noise)).append('\t');
		sb.append(Utils.rounded(simulationParameters.averageSignal / simulationParameters.noise)).append('\t');
		sb.append(config.getFittingWidth()).append('\t');
		sb.append(PSFProtosHelper.getName(fitConfig.getPSFType()));
		sb.append(":").append(PeakFit.getSolverName(fitConfig));
		if (fitConfig.isModelCameraMLE())
		{
			sb.append(":Camera\t");

			// Add details of the noise model for the MLE
			CalibrationReader r = new CalibrationReader(fitConfig.getCalibration());
			sb.append("EM=").append(r.isEMCCD());
			sb.append(":A=").append(Utils.rounded(r.getCountPerElectron()));
			sb.append(":N=").append(Utils.rounded(r.getReadNoise()));
			sb.append('\t');
		}
		else
			sb.append("\t\t");

		// Now output the actual results ...

		// Show histograms as cumulative to avoid problems with bin width
		// Residuals scores
		// Iterations and evaluations where fit was OK

		StoredDataStatistics[] stats = new StoredDataStatistics[NAMES2.length];
		for (int i = 0; i < stats.length; i++)
			stats[i] = new StoredDataStatistics();

		// For Jaccard scoring we need to count the score with no residuals threshold,
		// i.e. Accumulate the score accepting all doublets that were fit
		double tp = 0;
		double fp = 0;

		double bestTp = 0, bestFp = 0;

		ArrayList<DoubletBonus> data = new ArrayList<>(results.size());
		for (DoubletResult result : results)
		{
			final double score = result.getMaxScore();

			// Filter the singles that would be accepted
			if (result.good1)
			{
				// Filter the doublets that would be accepted
				if (result.good2)
				{
					final double tp2 = result.tp2a + result.tp2b;
					final double fp2 = result.fp2a + result.fp2b;
					tp += tp2;
					fp += fp2;

					if (result.tp2a > 0.5)
					{
						bestTp += result.tp2a;
						bestFp += result.fp2a;
					}
					if (result.tp2b > 0.5)
					{
						bestTp += result.tp2b;
						bestFp += result.fp2b;
					}

					// Store this as a doublet bonus
					data.add(new DoubletBonus(score, result.getAvScore(), tp2 - result.tp1, fp2 - result.fp1));
				}
				else
				{
					// No doublet fit so this will always be the single fit result
					tp += result.tp1;
					fp += result.fp1;
					if (result.tp1 > 0.5)
					{
						bestTp += result.tp1;
						bestFp += result.fp1;
					}
				}
			}

			// Build statistics
			final int c = result.c;

			// True results, i.e. where there was a choice between single or doublet
			if (result.valid)
			{
				stats[c].add(score);
			}

			// Of those where the fit was good, summarise the iterations and evaluations
			if (result.good1)
			{
				stats[3].add(result.iter1);
				stats[4].add(result.eval1);
				// Summarise only those which are a valid doublet. We do not really care
				// about the iteration increase for singles that are not doublets.
				if (c != 0 && result.good2)
				{
					stats[5].add(result.iter2);
					stats[6].add(result.eval2);
				}
			}
		}

		// Debug the counts
		//		double tpSingle = 0;
		//		double fpSingle = 0;
		//		double tpDoublet = 0;
		//		double fpDoublet = 0;
		//		int nSingle = 0, nDoublet = 0;
		//		for (DoubletResult result : results)
		//		{
		//			if (result.good1)
		//			{
		//				if (result.good2)
		//				{
		//					tpDoublet += result.tp2;
		//					fpDoublet += result.fp2;
		//					nDoublet++;
		//				}
		//				tpSingle += result.tp1;
		//				fpSingle += result.fp1;
		//				nSingle++;
		//			}
		//		}
		//		System.out.printf("Single %.1f,%.1f (%d) : Doublet %.1f,%.1f (%d)\n", tpSingle, fpSingle, nSingle, tpDoublet, fpDoublet, nDoublet*2);

		// Summarise score for true results
		Percentile p = new Percentile(99);
		for (int c = 0; c < stats.length; c++)
		{
			double[] values = stats[c].getValues();
			// Sorting is need for the percentile and the cumulative histogram so do it once
			Arrays.sort(values);
			sb.append(Utils.rounded(stats[c].getMean())).append("+/-")
					.append(Utils.rounded(stats[c].getStandardDeviation())).append(" (").append(stats[c].getN())
					.append(") ").append(Utils.rounded(p.evaluate(values))).append('\t');

			if (showHistograms && displayHistograms[c + NAMES.length])
				showHistogram(values, NAMES2[c]);
		}

		sb.append(MATCHING[matching]).append('\t');

		// Plot a graph of the additional results we would fit at all score thresholds.
		// This assumes we just pick the the doublet if we fit it (NO FILTERING at all!)

		// Initialise the score for residuals 0
		// Store this as it serves as a baseline for the filtering analysis
		computeScores(data, tp, fp, numberOfMolecules, true);
		_residualsScoreMax = this.residualsScore;
		computeScores(data, tp, fp, numberOfMolecules, false);
		_residualsScoreAv = this.residualsScore;

		residualsScore = (useMaxResiduals) ? _residualsScoreMax : _residualsScoreAv;
		if (showJaccardPlot)
			plotJaccard(residualsScore, null);

		String bestJaccard = Utils.rounded(bestTp / (bestFp + numberOfMolecules)) + '\t';
		analysisPrefix += bestJaccard;

		sb.append(bestJaccard);
		addJaccardScores(sb);

		sb.append('\t').append(Utils.timeToString(runTime / 1000000.0));

		summaryTable.append(sb.toString());
	}

	/**
	 * Compute the total number of true results that all the candidates could have fit, including singles, doublets and
	 * multiples.
	 *
	 * @param results
	 *            the doublet fitting results
	 * @return the total localisation results we could have fit
	 */
	private int countN(ArrayList<DoubletResult> results)
	{
		int n = 0;
		for (DoubletResult r : results)
			n += r.n;
		return n;
	}

	/**
	 * Compute maximum jaccard for all the residuals thresholds.
	 *
	 * @param data
	 *            the data
	 * @param tp
	 *            the true positives at residuals = 0
	 * @param fp
	 *            the false positives at residuals = 0
	 * @param n
	 *            the number of true positives at residuals = 0
	 * @param useMax
	 *            Use the max residuals
	 */
	private void computeScores(ArrayList<DoubletBonus> data, double tp, double fp, int n, boolean useMax)
	{
		// Add data at ends to complete the residuals scale from 0 to 1
		data.add(new DoubletBonus(0, 0, 0, 0));
		data.add(new DoubletBonus(1, 1, 0, 0));

		for (DoubletBonus b : data)
			b.setScore(useMax);
		Collections.sort(data);

		double[] residuals = new double[data.size() + 2];
		double[] jaccard = new double[residuals.length];
		double[] recall = new double[residuals.length];
		double[] precision = new double[residuals.length];
		int maxJaccardIndex = 0;

		int count = 0;
		double last = 0;
		for (DoubletBonus b : data)
		{
			if (last != b.residuals)
			{
				residuals[count] = last;
				jaccard[count] = tp / (n + fp);
				if (tp + fp != 0)
					precision[count] = tp / (tp + fp);
				recall[count] = tp / n;
				if (jaccard[maxJaccardIndex] < jaccard[count])
					maxJaccardIndex = count;
				count++;
			}
			tp -= b.tp;
			fp -= b.fp;
			last = b.residuals;
		}
		residuals[count] = last;
		jaccard[count] = tp / (n + fp);
		if (jaccard[maxJaccardIndex] < jaccard[count])
			maxJaccardIndex = count;
		if (tp + fp != 0)
			precision[count] = tp / (tp + fp);
		recall[count] = tp / n;
		count++;
		residuals = Arrays.copyOf(residuals, count);
		jaccard = Arrays.copyOf(jaccard, count);
		precision = Arrays.copyOf(precision, count);
		recall = Arrays.copyOf(recall, count);
		this.residualsScore = new ResidualsScore(residuals, jaccard, recall, precision, maxJaccardIndex);
	}

	private void plotJaccard(ResidualsScore residualsScore, ResidualsScore residualsScoreReference)
	{
		String title = TITLE + " Score";
		Plot plot = new Plot(title, "Residuals", "Score");
		double max = getMax(0.01, residualsScore);
		if (residualsScoreReference != null)
		{
			max = getMax(max, residualsScore);
		}
		plot.setLimits(0, 1, 0, max * 1.05);
		addLines(plot, residualsScore, 1);
		if (residualsScoreReference != null)
		{
			addLines(plot, residualsScoreReference, 0.5);
		}
		plot.setColor(Color.black);
		plot.addLabel(0, 0,
				String.format("Residuals %s; Jaccard %s (Black); Precision %s (Blue); Recall %s (Red)",
						Utils.rounded(residualsScore.residuals[residualsScore.maxJaccardIndex]),
						Utils.rounded(residualsScore.jaccard[residualsScore.maxJaccardIndex]),
						Utils.rounded(residualsScore.precision[residualsScore.maxJaccardIndex]),
						Utils.rounded(residualsScore.recall[residualsScore.maxJaccardIndex])));
		display(title, plot);
	}

	private double getMax(double max, ResidualsScore residualsScore)
	{
		max = Math.max(max, residualsScore.jaccard[residualsScore.maxJaccardIndex]);
		max = Maths.maxDefault(max, residualsScore.precision);
		max = Maths.maxDefault(max, residualsScore.recall);
		return max;
	}

	private void addLines(Plot plot, ResidualsScore residualsScore, double saturation)
	{
		if (saturation == 1)
		{
			// Colours are the same as the BenchmarkSpotFilter

			plot.setColor(Color.black);
			plot.addPoints(residualsScore.residuals, residualsScore.jaccard, Plot.LINE);
			plot.setColor(Color.blue);
			plot.addPoints(residualsScore.residuals, residualsScore.precision, Plot.LINE);
			plot.setColor(Color.red);
			plot.addPoints(residualsScore.residuals, residualsScore.recall, Plot.LINE);
		}
		else
		{
			plot.setColor(getColor(Color.black, saturation));
			plot.addPoints(residualsScore.residuals, residualsScore.jaccard, Plot.LINE);
			plot.setColor(getColor(Color.blue, saturation));
			plot.addPoints(residualsScore.residuals, residualsScore.precision, Plot.LINE);
			plot.setColor(getColor(Color.red, saturation));
			plot.addPoints(residualsScore.residuals, residualsScore.recall, Plot.LINE);
		}
	}

	private Color getColor(Color color, double saturation)
	{
		int r = color.getRed();
		int g = color.getGreen();
		int b = color.getBlue();
		if (r == 0 && g == 0 && b == 0)
		{
			// Black so make a grey
			return Color.gray;
		}
		float[] hsbvals = Color.RGBtoHSB(r, g, b, null);
		hsbvals[1] *= saturation;
		return Color.getHSBColor(hsbvals[0], hsbvals[1], hsbvals[2]);
	}

	/**
	 * Show a cumulative histogram of the data.
	 *
	 * @param values
	 *            the values
	 * @param xTitle
	 *            The name of plotted statistic
	 */
	public void showHistogram(double[] values, String xTitle)
	{
		double[][] h = Maths.cumulativeHistogram(values, false);

		String title = TITLE + " " + xTitle + " Cumulative";
		Plot2 plot = new Plot2(title, xTitle, "Frequency");
		double xMax = h[0][h[0].length - 1];
		double yMax = h[1][h[1].length - 1];
		plot.setLimits(0, xMax, 0, 1.05 * yMax);
		plot.setColor(Color.blue);
		plot.addPoints(h[0], h[1], Plot2.BAR);
		display(title, plot);
	}

	private void display(String title, Plot plot)
	{
		PlotWindow pw = Utils.display(title, plot);
		if (Utils.isNewWindow())
			windowOrganiser.add(pw.getImagePlus().getID());
	}

	/**
	 * Creates the summary table.
	 */
	private void createSummaryTable()
	{
		if (summaryTable == null || !summaryTable.isVisible())
		{
			summaryTable = new TextWindow(TITLE + " Summary", createSummaryHeader(), "", 1000, 300);
			summaryTable.setVisible(true);
		}
	}

	/**
	 * Creates the summary header.
	 *
	 * @return the string
	 */
	private String createSummaryHeader()
	{
		StringBuilder sb = new StringBuilder();
		sb.append(
				"Molecules\tMatched\tDensity\tminN\tmaxN\tN\ts (nm)\ta (nm)\tsa (nm)\tGain\tReadNoise (ADUs)\tB (photons)\tnoise (ADUs)\tSNR\tWidth\tMethod\tOptions\t");
		for (String name : NAMES2)
			sb.append(name).append('\t');
		sb.append(
				"Matching\tBest J\tJ (r=1)\tMax J\tResiduals\tArea +/-15%\tArea 98%\tMin 98%\tMax 98%\tRange 98%\twMean 98%\tArea >90%\tMin >90%\tMax >90%\tRange >90%\twMean >90%\tRun time");
		return sb.toString();
	}

	/**
	 * Show results.
	 *
	 * @param results
	 *            the results
	 */
	private void showResults(ArrayList<DoubletResult> results, boolean show)
	{
		if (!show)
			return;

		createResultsTable();

		ArrayList<String> list = new ArrayList<>(results.size());
		int flush = 9;
		StringBuilder sb = new StringBuilder();
		for (DoubletResult result : results)
		{
			sb.append(result.frame).append('\t');
			sb.append(result.spot.x).append('\t');
			sb.append(result.spot.y).append('\t');
			sb.append(IJ.d2s(result.spot.intensity, 1)).append('\t');
			sb.append(result.n).append('\t');
			sb.append(result.neighbours).append('\t');
			sb.append(result.almostNeighbours).append('\t');
			sb.append(Utils.rounded(result.score1)).append('\t');
			sb.append(Utils.rounded(result.score2)).append('\t');
			add(sb, result.fitResult1);
			add(sb, result.fitResult2);
			sb.append(IJ.d2s(result.sumOfSquares1, 1)).append('\t');
			sb.append(IJ.d2s(result.sumOfSquares2, 1)).append('\t');
			sb.append(IJ.d2s(result.ll1, 1)).append('\t');
			sb.append(IJ.d2s(result.ll2, 1)).append('\t');
			sb.append(IJ.d2s(result.value1, 1)).append('\t');
			sb.append(IJ.d2s(result.value2, 1)).append('\t');
			sb.append(Utils.rounded(result.r1)).append('\t');
			sb.append(Utils.rounded(result.r2)).append('\t');
			sb.append(Utils.rounded(result.aic1)).append('\t');
			sb.append(Utils.rounded(result.aic2)).append('\t');
			sb.append(Utils.rounded(result.bic1)).append('\t');
			sb.append(Utils.rounded(result.bic2)).append('\t');
			sb.append(Utils.rounded(result.maic1)).append('\t');
			sb.append(Utils.rounded(result.maic2)).append('\t');
			sb.append(Utils.rounded(result.mbic1)).append('\t');
			sb.append(Utils.rounded(result.mbic2)).append('\t');
			sb.append(Utils.rounded(result.a[0])).append('\t');
			sb.append(Utils.rounded(result.a[1])).append('\t');
			sb.append(Utils.rounded(result.gap)).append('\t');
			sb.append(Utils.rounded(result.xshift[0])).append('\t');
			sb.append(Utils.rounded(result.yshift[0])).append('\t');
			sb.append(Utils.rounded(result.xshift[1])).append('\t');
			sb.append(Utils.rounded(result.yshift[1])).append('\t');
			sb.append(result.iter1).append('\t');
			sb.append(result.iter2).append('\t');
			sb.append(result.eval1).append('\t');
			sb.append(result.eval2).append('\t');
			addParams(sb, result.fitResult1);
			addParams(sb, result.fitResult2);

			list.add(sb.toString());
			sb.setLength(0);
			// Flush below 10 lines so ImageJ will layout the columns
			if (--flush == 0)
			{
				resultsTable.getTextPanel().append(list);
				list.clear();
			}
		}

		resultsTable.getTextPanel().append(list);
	}

	/**
	 * Adds the.
	 *
	 * @param sb
	 *            the sb
	 * @param fitResult
	 *            the fit result
	 */
	private void add(StringBuilder sb, FitResult fitResult)
	{
		if (fitResult != null)
		{
			sb.append(fitResult.getStatus()).append('\t');
		}
		else
		{
			sb.append('\t');
		}
	}

	/**
	 * Adds the params.
	 *
	 * @param sb
	 *            the sb
	 * @param fitResult
	 *            the fit result
	 */
	private void addParams(StringBuilder sb, FitResult fitResult)
	{
		if (fitResult != null)
		{
			sb.append(Arrays.toString(fitResult.getParameters())).append('\t');
		}
		else
		{
			sb.append('\t');
		}
	}

	/**
	 * Creates the results table.
	 */
	private void createResultsTable()
	{
		if (resultsTable == null || !resultsTable.isVisible())
		{
			resultsTable = new TextWindow(TITLE + " Results", createResultsHeader(), "", 1000, 300);
			resultsTable.setVisible(true);
		}
	}

	/**
	 * Creates the results header.
	 *
	 * @return the string
	 */
	private String createResultsHeader()
	{
		return "Frame\tx\ty\tI\tn\tneighbours\talmost\tscore1\tscore2\tR1\tR2\tss1\tss2\tll1\tll2\tv1\tv2\tr1\tr2\taic1\taic2\tbic1\tbic2\tmaic1\tmaic2\tmbic1\tmbic2\ta1\ta2\tgap\tx1\ty1\tx2\ty2\ti1\ti2\te1\te2\tparams1\tparams2";
	}

	/**
	 * Run analysis.
	 */
	private void runAnalysis()
	{
		if (doubletResults == null)
		{
			IJ.error(TITLE, "No doublet results in memory");
			return;
		}

		// Ask the user to set filters
		if (!showAnalysisDialog())
			return;

		showResults(doubletResults, analysisShowResults);

		// Store the effect of fitting as a doublet
		ArrayList<DoubletBonus> data = new ArrayList<>(doubletResults.size());
		// True positive and False positives at residuals = 0
		double tp = 0;
		double fp = 0;

		Logger logger = (analysisLogging) ? new IJLogger() : null;

		// Get filters for the single and double fits

		// No coordinate shift for the doublet. We have already done simple checking of the
		// coordinates to get the good=2 flag
		FitConfiguration filterFitConfig2 = filterFitConfig.clone();
		filterFitConfig2.setCoordinateShift(Integer.MAX_VALUE);

		final int size = 2 * config.getFittingWidth() + 1;
		Rectangle regionBounds = new Rectangle(0, 0, size, size);
		final double otherDriftAngle = 180 - analysisDriftAngle;

		// Process all the results
		for (DoubletResult result : doubletResults)
		{
			// Filter the singles that would be accepted
			if (result.good1)
			{
				filterFitConfig.setNoise(result.noise);
				FitStatus fitStatus0 = filterFitConfig.validatePeak(0, result.fitResult1.getInitialParameters(),
						result.fitResult1.getParameters(), result.fitResult1.getParameterDeviations());

				double tp1 = 0, fp1 = 0;
				if (fitStatus0 == FitStatus.OK)
				{
					tp1 = result.tp1;
					fp1 = result.fp1;
				}
				else if (analysisLogging)
					logFailure(logger, 0, result, fitStatus0);

				// Filter the doublets that would be accepted
				// We have already done simple checking to get the good1 = true flag. So accept
				// width diverged spots as OK for a doublet fit
				if ((fitStatus0 == FitStatus.OK || fitStatus0 == FitStatus.WIDTH_DIVERGED) && selectFit(result) &&
						result.good2)
				{
					double tp2 = 0, fp2 = 0;

					// Basic spot criteria (SNR, Photons, width)
					filterFitConfig2.setNoise(result.noise);
					FitStatus fitStatus1 = filterFitConfig2.validatePeak(0, result.fitResult2.getInitialParameters(),
							result.fitResult2.getParameters(), result.fitResult2.getParameterDeviations());
					FitStatus fitStatus2 = filterFitConfig2.validatePeak(1, result.fitResult2.getInitialParameters(),
							result.fitResult2.getParameters(), result.fitResult2.getParameterDeviations());

					// Log basic failures
					boolean[] accept = new boolean[2];
					if (fitStatus1 == FitStatus.OK)
					{
						accept[0] = true;
					}
					else if (analysisLogging)
						logFailure(logger, 1, result, fitStatus1);

					if (fitStatus2 == FitStatus.OK)
					{
						accept[1] = true;
					}
					else if (analysisLogging)
						logFailure(logger, 2, result, fitStatus2);

					// If the basic filters are OK, do some analysis of the doublet.
					// We can filter each spot with criteria such as shift and the angle to the quadrant.
					if (accept[0] || accept[1])
					{
						if (result.gap < minGap)
						{
							accept[0] = accept[1] = false;
							if (analysisLogging)
								logger.info("Reject Doublet (%.2f): Fitted coordinates below min gap (%g<%g)\n",
										result.getMaxScore(), result.gap, minGap);
						}
					}

					if (accept[0] || accept[1])
					{
						// The logic in here will be copied to the FitWorker.quadrantAnalysis routine.
						double[] params = result.fitResult1.getParameters();
						double[] newParams = result.fitResult2.getParameters();

						// Set up for shift filtering
						double shift = filterFitConfig.getCoordinateShift();
						if (shift == 0 || shift == Double.POSITIVE_INFINITY)
						{
							// Allow the shift to span half of the fitted window.
							shift = 0.5 * FastMath.min(regionBounds.width, regionBounds.height);
						}

						// Set an upper limit on the shift that is not too far outside the fit window
						final double maxShiftX, maxShiftY;
						final double factor = Gaussian2DFunction.SD_TO_HWHM_FACTOR;
						if (fitConfig.isXSDFitting())
						{
							// Add the fitted standard deviation to the allowed shift
							maxShiftX = regionBounds.width * 0.5 + factor * params[Gaussian2DFunction.X_SD];
							maxShiftY = regionBounds.height * 0.5 + factor * params[Gaussian2DFunction.Y_SD];
						}
						else
						{
							// Add the configured standard deviation to the allowed shift
							maxShiftX = regionBounds.width * 0.5 + factor * fitConfig.getInitialXSD();
							maxShiftY = regionBounds.height * 0.5 + factor * fitConfig.getInitialYSD();
						}

						for (int n = 0; n < 2; n++)
						{
							if (!accept[n])
								continue;
							accept[n] = false; // Reset

							final double xShift = newParams[Gaussian2DFunction.X_POSITION +
									n * Gaussian2DFunction.PARAMETERS_PER_PEAK] - params[Gaussian2DFunction.X_POSITION];
							final double yShift = newParams[Gaussian2DFunction.Y_POSITION +
									n * Gaussian2DFunction.PARAMETERS_PER_PEAK] - params[Gaussian2DFunction.Y_POSITION];
							if (Math.abs(xShift) > maxShiftX || Math.abs(yShift) > maxShiftY)
							{
								if (analysisLogging)
									logger.info(
											"Reject P%d (%.2f): Fitted coordinates moved outside fit region (x=%g,y=%g)\n",
											n + 1, result.getMaxScore(), xShift, yShift);
								continue;
							}
							if (Math.abs(xShift) > shift || Math.abs(yShift) > shift)
							{
								// Check the domain is OK (the angle is in radians).
								// Allow up to a 45 degree difference to show the shift is along the vector
								if (result.a[n] > analysisDriftAngle && result.a[n] < otherDriftAngle)
								{
									if (analysisLogging)
										logger.info(
												"Reject P%d (%.2f): Fitted coordinates moved into wrong quadrant (x=%g,y=%g,a=%f)",
												n + 1, result.getMaxScore(), xShift, yShift, result.a[n]);
									continue;
								}

								// Note: The FitWorker also checks for drift to another candidate.
							}

							// This is OK
							accept[n] = true;
						}
					}

					if (accept[0])
					{
						tp2 += result.tp2a;
						fp2 += result.fp2a;
					}
					if (accept[1])
					{
						tp2 += result.tp2b;
						fp2 += result.fp2b;
					}

					if (accept[0] || accept[1])
					{
						tp += tp2;
						fp += fp2;

						// Store this as a doublet bonus
						data.add(new DoubletBonus(result.getMaxScore(), result.getAvScore(), tp2 - tp1, fp2 - fp1));
					}
					else
					{
						// No doublet fit so this will always be the single fit result
						tp += tp1;
						fp += fp1;
					}
				}
				else
				{
					// No doublet fit so this will always be the single fit result
					tp += tp1;
					fp += fp1;
				}
			}
		}

		// Compute the max Jaccard
		computeScores(data, tp, fp, numberOfMolecules, useMaxResiduals);

		if (showJaccardPlot)
			plotJaccard(residualsScore, (useMaxResiduals) ? _residualsScoreMax : _residualsScoreAv);

		createAnalysisTable();

		StringBuilder sb = new StringBuilder(analysisPrefix);
		sb.append(analysisTitle).append('\t');
		sb.append((useMaxResiduals) ? "Max" : "Average").append('\t');
		sb.append(SELECTION_CRITERIA[selectionCriteria]).append('\t');
		if (filterFitConfig.isSmartFilter())
		{
			sb.append(filterFitConfig.getSmartFilterName()).append("\t\t\t\t\t\t\t\t");
		}
		else
		{
			sb.append('\t');
			sb.append(filterFitConfig.getCoordinateShiftFactor()).append('\t');
			sb.append(filterFitConfig.getSignalStrength()).append('\t');
			sb.append(filterFitConfig.getMinPhotons()).append('\t');
			sb.append(filterFitConfig.getMinWidthFactor()).append('\t');
			sb.append(filterFitConfig.getMaxWidthFactor()).append('\t');
			sb.append(filterFitConfig.getPrecisionThreshold()).append('\t');
			sb.append(filterFitConfig.getPrecisionMethod()).append('\t');
		}
		sb.append(analysisDriftAngle).append('\t');
		sb.append(minGap).append('\t');

		addJaccardScores(sb);

		analysisTable.append(sb.toString());

		saveTemplate(sb.toString());
	}

	/**
	 * Save PeakFit configuration template using the current benchmark settings.
	 *
	 * @param summary
	 *            the summary
	 */
	private void saveTemplate(String summary)
	{
		if (!saveTemplate)
			return;

		// Start with a clone of the filter settings
		FitEngineConfiguration config = new FitEngineConfiguration();
		FitConfiguration fitConfig = config.getFitConfiguration();
		fitConfig.setFitSettings(filterFitConfig.getFitSettings());

		// Copy settings used during fitting
		updateConfiguration(config);

		// Remove the PSF width to make the template generic
		fitConfig.setInitialPeakStdDev(0);
		fitConfig.setNmPerPixel(0);
		fitConfig.setGain(0);
		fitConfig.setNoise(0);
		// This was done fitting all the results
		config.setFailuresLimit(-1);
		if (useBenchmarkSettings)
		{
			FitEngineConfiguration pConfig = new FitEngineConfiguration();
			// TODO - add option to use latest or the best
			if (BenchmarkFilterAnalysis.updateConfiguration(pConfig, false))
				config.setFailuresLimit(pConfig.getFailuresLimit());
		}

		// Set the residuals
		fitConfig.setComputeResiduals(true);
		// TODO - make the choice of the best residuals configurable
		config.setResidualsThreshold(residualsScore.bestResiduals[2]);

		String filename = BenchmarkFilterAnalysis.getFilename("Template_File", templateFilename);
		if (filename != null)
		{
			templateFilename = filename;
			TemplateSettings.Builder settings = TemplateSettings.newBuilder();
			getNotes(settings, summary);
			settings.setFitEngineSettings(config.getFitEngineSettings());
			if (!SettingsManager.toJSON(settings.build(), filename, SettingsManager.FLAG_SILENT))
				IJ.log("Unable to save the template configuration");
		}
	}

	/**
	 * Updates the given configuration using the latest fitting settings used in benchmarking.
	 *
	 * @param pConfig
	 *            the configuration
	 * @return true, if successful
	 */
	public static boolean updateConfiguration(FitEngineConfiguration pConfig)
	{
		pConfig.mergeFitEngineSettings(config.getFitEngineSettings());
		return true;
	}

	private void getNotes(TemplateSettings.Builder settings, String summary)
	{
		settings.addNotes("Benchmark template");
		if (!TextUtils.isNullOrEmpty(analysisTitle))
			BenchmarkFilterAnalysis.addField(settings, "Doublet Analysis Title", analysisTitle);
		// Add create data settings.
		// Just add the columns and the data from the summary window
		final String header = createAnalysisHeader();
		BenchmarkFilterAnalysis.addField(settings, "Doublet Analysis Summary Fields", header);
		BenchmarkFilterAnalysis.addField(settings, "Doublet Analysis Summary Values", summary);
		// Now pick out key values...
		BenchmarkFilterAnalysis.addKeyFields(settings, header, summary, new String[] { "Density", "s", "Selection",
				"Max J", "Residuals", "Area >90", "Range >90", "wMean >90" });

		// Add any other settings that may be useful in the template
		BenchmarkFilterAnalysis.addField(settings, "Created", BenchmarkFilterAnalysis.getCurrentTimeStamp());
	}

	private void logFailure(Logger logger, int i, DoubletResult result, FitStatus fitStatus)
	{
		logger.info("Reject P%d (%.2f): %s\n", i, result.getMaxScore(), fitStatus.toString());
	}

	private boolean selectFit(DoubletResult result)
	{
		switch (selectionCriteria)
		{
			//@formatter:off
			case 0:	return result.r2    > result.r1;
			case 1:	return result.aic2  < result.aic1;
			case 2:	return result.bic2  < result.bic1;
			case 3:	return result.maic2 < result.maic1;
			case 4:	return result.mbic2 < result.mbic1;
			//@formatter:on
		}
		return false;
	}

	private void addJaccardScores(StringBuilder sb)
	{
		double[] residuals = residualsScore.residuals;
		double[] jaccard = residualsScore.jaccard;
		int maxJaccardIndex = residualsScore.maxJaccardIndex;
		sb.append(Utils.rounded(jaccard[jaccard.length - 1])).append('\t');
		sb.append(Utils.rounded(jaccard[maxJaccardIndex])).append('\t')
				.append(Utils.rounded(residuals[maxJaccardIndex]));

		sb.append('\t').append(Utils.rounded(getArea(residuals, jaccard, maxJaccardIndex, 0.15)));
		//sb.append('\t').append(Utils.rounded(getArea(residuals, jaccard, maxJaccardIndex, 0.3)));
		//sb.append('\t').append(Utils.rounded(getArea(residuals, jaccard, maxJaccardIndex, 1)));

		residualsScore.bestResiduals[0] = residuals[maxJaccardIndex];
		// Find the range that has a Jaccard within a % of the max
		residualsScore.bestResiduals[1] = addRange(sb, residuals, jaccard, maxJaccardIndex,
				jaccard[maxJaccardIndex] * 0.98);
		// Do the same within a fraction of the performance improvement over no residuals
		residualsScore.bestResiduals[2] = addRange(sb, residuals, jaccard, maxJaccardIndex,
				jaccard[jaccard.length - 1] + (jaccard[maxJaccardIndex] - jaccard[jaccard.length - 1]) * 0.9);
	}

	private double addRange(StringBuilder sb, double[] residuals, double[] jaccard, int maxJaccardIndex, double limit)
	{
		double sum = 0;
		double lower = residuals[maxJaccardIndex];
		double upper = residuals[maxJaccardIndex];
		int j = maxJaccardIndex;
		// For weighted mean
		double sumR = 0;
		for (int i = j; i-- > 0;)
		{
			if (jaccard[i] < limit)
			{
				break;
			}
			final double height = (jaccard[j] + jaccard[i]) * 0.5;
			final double width = (residuals[j] - residuals[i]);
			final double area = height * width;
			sum += area;
			j = i;
			lower = residuals[j];
			final double avResiduals = (residuals[j] + residuals[i]) * 0.5;
			sumR += area * avResiduals;
		}
		j = maxJaccardIndex;
		for (int i = j; ++i < jaccard.length;)
		{
			if (jaccard[i] < limit)
			{
				break;
			}
			final double height = (jaccard[j] + jaccard[i]) * 0.5;
			final double width = (residuals[i] - residuals[j]);
			final double area = height * width;
			sum += area;
			j = i;
			upper = residuals[j];
			final double avResiduals = (residuals[j] + residuals[i]) * 0.5;
			sumR += area * avResiduals;
		}
		final double weightedMean;
		if (sum != 0)
		{
			weightedMean = sumR / sum;
			sum /= (upper - lower);
		}
		else
		{
			weightedMean = sum = jaccard[maxJaccardIndex];
		}
		sb.append('\t').append(Utils.rounded(sum)).append('\t').append(Utils.rounded(lower)).append('\t')
				.append(Utils.rounded(upper)).append('\t').append(Utils.rounded(upper - lower)).append('\t')
				.append(Utils.rounded(weightedMean));
		return weightedMean;
	}

	private double getArea(double[] residuals, double[] jaccard, int maxJaccardIndex, double window)
	{
		double sum = 0;
		double lower = Math.max(0, residuals[maxJaccardIndex] - window);
		double upper = Math.min(1, residuals[maxJaccardIndex] + window);
		//		LinearInterpolator li = new LinearInterpolator();
		//		PolynomialSplineFunction fun = li.interpolate(residuals, jaccard);
		//		//TrapezoidIntegrator in = new TrapezoidIntegrator();
		//		SimpsonIntegrator in = new SimpsonIntegrator();
		//
		//		try
		//		{
		//			sum = in.integrate(1000, fun, lower, upper);
		//		}
		//		catch (TooManyEvaluationsException e)
		//		{
		//			System.out.println("Failed to calculate the area: " + e.getMessage());
		int j = maxJaccardIndex;
		for (int i = j; i-- > 0;)
		{
			if (residuals[i] < lower)
			{
				lower = residuals[j];
				break;
			}
			sum += (jaccard[j] + jaccard[i]) * 0.5 * (residuals[j] - residuals[i]);
			j = i;
		}
		j = maxJaccardIndex;
		for (int i = j; ++i < jaccard.length;)
		{
			if (residuals[i] > upper)
			{
				upper = residuals[j];
				break;
			}
			sum += (jaccard[j] + jaccard[i]) * 0.5 * (residuals[i] - residuals[j]);
			j = i;
		}
		//		}
		return sum / (upper - lower);
	}

	/**
	 * Show dialog.
	 *
	 * @return true, if successful
	 */
	private boolean showAnalysisDialog()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		StringBuilder sb = new StringBuilder("Filters the doublet fits and reports the performance increase\n");

		// Show the fitting settings that will effect filters, i.e. fit standard deviation, fit width
		sb.append("SD0 = ").append(Utils.rounded(fitConfig.getInitialXSD())).append("\n");
		//sb.append("SD1 = ").append(Utils.rounded(fitConfig.getInitialPeakStdDev1())).append("\n");
		sb.append("Fit Width = ").append(config.getFittingWidth()).append("\n");

		gd.addMessage(sb.toString());

		// Collect options for filtering
		gd.addChoice("Selection_Criteria", SELECTION_CRITERIA, SELECTION_CRITERIA[selectionCriteria]);

		// Copy the settings used when fitting
		filterFitConfig.setCalibration(fitConfig.getCalibration());
		filterFitConfig.setPSF(fitConfig.getPSF());
		filterFitConfig.setFitSolverSettings(fitConfig.getFitSolverSettings());

		String[] templates = ConfigurationTemplate.getTemplateNames(true);
		gd.addChoice("Template", templates, templates[0]);
		// Allow the settings from the benchmark analysis to be used
		gd.addCheckbox("Benchmark_settings", analysisUseBenchmarkSettings);
		gd.addCheckbox("Smart_filter", fitConfig.isSmartFilter());
		gd.addSlider("Shift_factor", 0.01, 2, filterFitConfig.getCoordinateShiftFactor());
		gd.addNumericField("Signal_strength", filterFitConfig.getSignalStrength(), 2);
		gd.addNumericField("Min_photons", filterFitConfig.getMinPhotons(), 0);
		gd.addSlider("Min_width_factor", 0, 0.99, filterFitConfig.getMinWidthFactor());
		gd.addSlider("Max_width_factor", 1.01, 5, filterFitConfig.getMaxWidthFactor());
		gd.addNumericField("Precision", filterFitConfig.getPrecisionThreshold(), 2);
		gd.addChoice("Precision_method", SettingsManager.getPrecisionMethodNames(),
				filterFitConfig.getPrecisionMethod().ordinal());

		gd.addNumericField("Drift_angle", analysisDriftAngle, 2);
		gd.addNumericField("Min_gap", minGap, 2);

		// Collect display options
		gd.addCheckbox("Show_results", analysisShowResults);
		gd.addCheckbox("Show_Jaccard_Plot", showJaccardPlot);
		gd.addCheckbox("Use_max_residuals", useMaxResiduals);
		gd.addCheckbox("Logging", analysisLogging);
		gd.addStringField("Title", analysisTitle);
		gd.addCheckbox("Save_template", saveTemplate);

		// TODO - Add support for updating a template with a residuals threshold, e.g. from the BenchmarkFilterAnalysis plugin

		// Add a mouse listener to the config file field
		if (Utils.isShowGenericDialog())
		{
			Vector<TextField> numerics = gd.getNumericFields();
			Vector<Checkbox> checkboxes = gd.getCheckboxes();
			Vector<Choice> choices = gd.getChoices();
			int n = 0;
			choices.get(1).addItemListener(this);
			checkboxes.get(0).addItemListener(this);
			cbSmartFilter = checkboxes.get(1);
			textCoordinateShiftFactor = numerics.get(n++);
			textSignalStrength = numerics.get(n++);
			textMinPhotons = numerics.get(n++);
			textMinWidthFactor = numerics.get(n++);
			textWidthFactor = numerics.get(n++);
			textPrecisionThreshold = numerics.get(n++);
			textPrecisionMethod = choices.get(2);
		}

		gd.showDialog();

		if (gd.wasCanceled())
			return false;
		if (gd.invalidNumber())
			return false;

		selectionCriteria = gd.getNextChoiceIndex();
		// Ignore the template
		gd.getNextChoice();
		analysisUseBenchmarkSettings = gd.getNextBoolean();
		fitConfig.setSmartFilter(gd.getNextBoolean());
		filterFitConfig.setCoordinateShiftFactor(gd.getNextNumber());
		filterFitConfig.setSignalStrength(gd.getNextNumber());
		filterFitConfig.setMinPhotons(gd.getNextNumber());
		filterFitConfig.setMinWidthFactor(gd.getNextNumber());
		filterFitConfig.setWidthFactor(gd.getNextNumber());
		filterFitConfig.setPrecisionThreshold(gd.getNextNumber());
		filterFitConfig.setPrecisionMethod(gd.getNextChoiceIndex());
		analysisDriftAngle = gd.getNextNumber();
		minGap = gd.getNextNumber();

		analysisShowResults = gd.getNextBoolean();
		showJaccardPlot = gd.getNextBoolean();
		useMaxResiduals = gd.getNextBoolean();
		analysisLogging = gd.getNextBoolean();
		analysisTitle = gd.getNextString();
		saveTemplate = gd.getNextBoolean();

		if (gd.invalidNumber())
			return false;

		if (analysisUseBenchmarkSettings)
		{
			if (!updateFilterConfiguration(filterFitConfig))
				return false;
		}

		return true;
	}

	private boolean updateFilterConfiguration(FitConfiguration filterFitConfig)
	{
		FitEngineConfiguration c = new FitEngineConfiguration();
		// TODO - add option to use latest or the best
		if (!BenchmarkFilterAnalysis.updateConfiguration(c, false))
		{
			IJ.error(TITLE, "Unable to use the benchmark filter analysis configuration");
			return false;
		}
		filterFitConfig.setFitSettings(c.getFitEngineSettings().getFitSettings());
		return true;
	}

	/**
	 * Creates the analysis table.
	 */
	private void createAnalysisTable()
	{
		if (analysisTable == null || !analysisTable.isVisible())
		{
			analysisTable = new TextWindow(TITLE + " Analysis", createAnalysisHeader(), "", 1200, 300);
			analysisTable.setVisible(true);
		}
	}

	/**
	 * Creates the analysis header.
	 *
	 * @return the string
	 */
	private String createAnalysisHeader()
	{
		return "Density\ts\tWidth\tMethod\tOptions\tBest J\tTitle\tUse residuals\tSelection\tFilter\tShift\tSNR\tPhotons\tMin Width\tWidth\tPrecision\tLocal B\tAngle\tGap\tJ (r=1)\tMax J\tResiduals\tArea +/-15%\tArea 98%\tMin 98%\tMax 98%\tRange 98%\twMean 98%\tArea >90%\tMin >90%\tMax >90%\tRange >90%\twMean >90%";
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see java.awt.event.ItemListener#itemStateChanged(java.awt.event.ItemEvent)
	 */
	@Override
	public void itemStateChanged(ItemEvent e)
	{
		if (e.getSource() instanceof Choice)
		{
			// Update the settings from the template
			Choice choice = (Choice) e.getSource();
			String templateName = choice.getSelectedItem();

			// Get the configuration template
			TemplateSettings template = ConfigurationTemplate.getTemplate(templateName);

			if (textCoordinateShiftFactor != null)
			{
				// Start with a default. This will cause a reset.
				FitConfiguration fitConfig = new FitConfiguration();
				if (template != null)
				{
					fitConfig.setFitSettings(template.getFitEngineSettings().getFitSettings());
				}

				cbSmartFilter.setState(fitConfig.isSmartFilter());
				textCoordinateShiftFactor.setText("" + fitConfig.getCoordinateShiftFactor());
				textSignalStrength.setText("" + fitConfig.getSignalStrength());
				textMinPhotons.setText("" + fitConfig.getMinPhotons());
				textMinWidthFactor.setText("" + fitConfig.getMinWidthFactor());
				textWidthFactor.setText("" + fitConfig.getMaxWidthFactor());
				textPrecisionThreshold.setText("" + fitConfig.getPrecisionThreshold());
				textPrecisionMethod.select(fitConfig.getPrecisionMethod().ordinal());
			}
			else
			{
				if (template != null)
				{
					if (template.hasFitEngineSettings())
					{
						boolean custom = ConfigurationTemplate.isCustomTemplate(templateName);
						FitEngineConfiguration config2 = new FitEngineConfiguration(template.getFitEngineSettings(),
								template.getCalibration(), template.getPsf());
						fitConfig = new FitConfiguration();
						FitConfiguration fitConfig2 = config2.getFitConfiguration();
						if (custom && template.hasPsf())
							textPSF.select(PeakFit.getPSFTypeNames()[fitConfig2.getPSFType().ordinal()]);
						textDataFilterType.select(config2.getDataFilterType().ordinal());
						textDataFilter.select(config2.getDataFilterMethod(0).ordinal());
						textSmooth.setText("" + config2.getDataFilterParameterValue(0));
						textSearch.setText("" + config2.getSearch());
						textBorder.setText("" + config2.getBorder());
						textFitting.setText("" + config2.getFitting());
						textFitSolver.select(SettingsManager.getFitSolverNames()[fitConfig.getFitSolver().ordinal()]);

						// Copy settings not in the dialog for the fit solver
						if (custom)
							fitConfig.setPSF(fitConfig2.getPSF());
						fitConfig.setFitSolverSettings(fitConfig2.getFitSolverSettings());
					}
				}
				else
				{
					// Ignore
				}
			}
		}
		else if (e.getSource() instanceof Checkbox)
		{
			Checkbox checkbox = (Checkbox) e.getSource();
			if (!checkbox.getState())
				return;

			if (textCoordinateShiftFactor != null)
			{
				if (!updateFilterConfiguration(filterFitConfig))
					return;

				cbSmartFilter.setState(filterFitConfig.isSmartFilter());
				textCoordinateShiftFactor.setText("" + filterFitConfig.getCoordinateShiftFactor());
				textSignalStrength.setText("" + filterFitConfig.getSignalStrength());
				textMinPhotons.setText("" + filterFitConfig.getMinPhotons());
				textMinWidthFactor.setText("" + filterFitConfig.getMinWidthFactor());
				textWidthFactor.setText("" + filterFitConfig.getMaxWidthFactor());
				textPrecisionThreshold.setText("" + filterFitConfig.getPrecisionThreshold());
				textPrecisionMethod.select(filterFitConfig.getPrecisionMethod().ordinal());
			}
			else
			{
				if (!updateFitConfiguration(config))
					return;

				textPSF.select(PeakFit.getPSFTypeNames()[fitConfig.getPSFType().ordinal()]);
				textDataFilterType.select(config.getDataFilterType().ordinal());
				textDataFilter.select(config.getDataFilterMethod(0).ordinal());
				textSmooth.setText("" + config.getDataFilterParameterValue(0));
				textSearch.setText("" + config.getSearch());
				textBorder.setText("" + config.getBorder());
				textFitting.setText("" + config.getFitting());
				textFitSolver.select(SettingsManager.getFitSolverNames()[fitConfig.getFitSolver().ordinal()]);
				textMatchDistance.setText("" + matchDistance);
				textLowerDistance.setText("" + lowerDistance);
				textSignalFactor.setText("" + signalFactor);
				textLowerFactor.setText("" + lowerSignalFactor);
			}
		}
	}
}
