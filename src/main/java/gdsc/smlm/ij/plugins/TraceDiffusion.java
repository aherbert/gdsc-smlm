package gdsc.smlm.ij.plugins;

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

import gdsc.smlm.ij.IJTrackProgress;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.settings.ClusteringSettings;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.Trace;
import gdsc.smlm.results.TraceManager;
import gdsc.smlm.utils.Maths;
import gdsc.smlm.utils.Random;
import gdsc.smlm.utils.Sort;
import gdsc.smlm.utils.Statistics;
import gdsc.smlm.utils.StoredDataStatistics;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import ij.plugin.WindowOrganiser;
import ij.text.TextWindow;

import java.awt.Color;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.MultivariateMatrixFunction;
import org.apache.commons.math3.analysis.MultivariateVectorFunction;
import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.exception.TooManyIterationsException;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.OptimizationData;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.PointVectorValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.SimpleValueChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer;
import org.apache.commons.math3.optim.nonlinear.vector.ModelFunction;
import org.apache.commons.math3.optim.nonlinear.vector.ModelFunctionJacobian;
import org.apache.commons.math3.optim.nonlinear.vector.Target;
import org.apache.commons.math3.optim.nonlinear.vector.Weight;
import org.apache.commons.math3.optim.nonlinear.vector.jacobian.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.util.FastMath;

/**
 * Run a tracing algorithm on the peak results to trace molecules across the frames.
 */
public class TraceDiffusion implements PlugIn
{
	private static final String TITLE = "Trace Diffusion";
	private static String inputOption = "";
	private static String header = null;
	private static TextWindow summaryTable = null;

	private static final String[] NAMES = new String[] { "Total Signal", "Signal/Frame", "t-On (s)" };
	private static boolean[] displayHistograms = new boolean[NAMES.length];
	static
	{
		for (int i = 0; i < displayHistograms.length; i++)
			displayHistograms[i] = false;
	}
	private static final int TOTAL_SIGNAL = 0;
	private static final int SIGNAL_PER_FRAME = 1;
	private static final int T_ON = 2;

	private static boolean[] alwaysRemoveOutliers;
	static
	{
		alwaysRemoveOutliers = new boolean[NAMES.length];
		alwaysRemoveOutliers[TOTAL_SIGNAL] = false;
	}

	private static boolean displayMSDHistogram = true;
	private static boolean displayDHistogram = true;

	private static boolean saveTraceDistances = false;
	private static String filename = "";
	private static double minFraction = 0.1;
	private static double minDifference = 2;

	private GlobalSettings globalSettings;
	private ClusteringSettings settings;
	private MemoryPeakResults results;
	// Store exposure time in seconds
	private double exposureTime = 0;
	private double precision;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		if (MemoryPeakResults.countMemorySize() == 0)
		{
			IJ.error(TITLE, "No localisations in memory");
			return;
		}

		if (!showDialog())
			return;

		TraceManager manager = new TraceManager(results);

		// Run the tracing
		manager.setTracker(new IJTrackProgress());
		manager.setDistanceExclusion(settings.distanceExclusion / results.getCalibration().nmPerPixel);
		manager.traceMolecules(settings.distanceThreshold / results.getCalibration().nmPerPixel, 1);
		Trace[] traces = manager.getTraces();

		traces = filterTraces(traces, settings.minimumTraceLength);

		int count = traces.length;
		double D = 0;
		double[][] jdParams = null;
		if (count > 0)
		{
			//--- Save results ---

			// Save the traces to memory
			TraceMolecules.saveResults(results, traces, "Tracks");

			// Sort traces by time to assist the results source in extracting frames sequentially.
			// Do this before saving to assist in debugging using the saved traces file.
			TraceMolecules.sortByTime(traces);

			if (settings.saveTraces)
				TraceMolecules.saveTraces(results, traces, createSettingsComment());

			//--- MSD Analysis ---

			// Conversion constants
			final double px2ToUm2 = results.getCalibration().nmPerPixel * results.getCalibration().nmPerPixel / 1e6;
			final double px2ToUm2PerSecond = px2ToUm2 / exposureTime;

			// Get the maximum trace length
			int length = settings.minimumTraceLength;
			if (!settings.truncate)
			{
				for (Trace trace : traces)
					length = FastMath.max(length, trace.size());
			}

			// Extract the mean-squared distance statistics
			Statistics[] stats = new Statistics[length];
			// Disable sub-sampling
			final boolean subSample = false; // (settings.internalDistances & settings.subSampledDistances);
			for (int i = 0; i < stats.length; i++)
				stats[i] = (subSample) ? new StoredDataStatistics() : new Statistics();

			ArrayList<double[]> distances = (saveTraceDistances) ? new ArrayList<double[]>(traces.length) : null;

			// Store all the jump distances at the specified interval
			StoredDataStatistics jumpDistances = new StoredDataStatistics();
			final int jumpDistanceInterval = settings.jumpDistance;
			final double jdPx2ToUm2PerSecond = px2ToUm2 / (jumpDistanceInterval * exposureTime);

			// Compute squared distances
			StoredDataStatistics msdPerMoleculeAllVsAll = new StoredDataStatistics();
			StoredDataStatistics msdPerMoleculeAdjacent = new StoredDataStatistics();
			for (Trace trace : traces)
			{
				ArrayList<PeakResult> results = trace.getPoints();
				// Sum the MSD and the time
				double sumD = 0, sumD_adjacent = 0;
				int sumT = 0, sumT_adjacent = 0;
				final int traceLength = (settings.truncate) ? settings.minimumTraceLength : trace.size();

				// Get the mean for each time separation
				double[] sum = new double[traceLength + 1];

				// Do the distances to the origin (saving if necessary)
				{
					final float x = results.get(0).getXPosition();
					final float y = results.get(0).getYPosition();
					if (saveTraceDistances)
					{
						double[] msd = new double[traceLength - 1];
						for (int j = 1; j < traceLength; j++)
						{
							final double d = distance2(x, y, results.get(j));
							msd[j - 1] = px2ToUm2 * d;
							final int t = j;
							if (t == 1)
							{
								sumD_adjacent += d;
								sumT_adjacent++;
							}
							else
							{
								sumD += d;
								sumT += t;
							}
							if (t == jumpDistanceInterval)
								jumpDistances.add(jdPx2ToUm2PerSecond * d);
							sum[t] += d;
						}
						distances.add(msd);
					}
					else
					{
						for (int j = 1; j < traceLength; j++)
						{
							final double d = distance2(x, y, results.get(j));
							final int t = j;
							if (t == 1)
							{
								sumD_adjacent += d;
								sumT_adjacent++;
							}
							else
							{
								sumD += d;
								sumT += t;
							}
							if (t == jumpDistanceInterval)
								jumpDistances.add(jdPx2ToUm2PerSecond * d);
							sum[t] += d;
						}
					}
				}

				if (settings.internalDistances)
				{
					// Do the internal distances
					for (int i = 1; i < traceLength; i++)
					{
						final float x = results.get(i).getXPosition();
						final float y = results.get(i).getYPosition();
						for (int j = i + i; j < traceLength; j++)
						{
							final double d = distance2(x, y, results.get(j));
							final int t = j - i;
							if (t == 1)
							{
								sumD_adjacent += d;
								sumT_adjacent++;
							}
							else
							{
								sumD += d;
								sumT += t;
							}
							if (t == jumpDistanceInterval)
								jumpDistances.add(jdPx2ToUm2PerSecond * d);
							sum[t] += d;
						}
					}

					// Add the average distance per time separation to the population
					for (int t = 1; t < traceLength; t++)
					{
						stats[t].add(sum[t] / (traceLength - t));
					}
				}
				else
				{
					// Add the distance per time separation to the population
					for (int t = 1; t < traceLength; t++)
					{
						stats[t].add(sum[t]);
					}
				}

				// Calculate the average displacement for the trace (do not simply use the largest 
				// time separation since this will miss moving molecules that end up at the origin)
				sumD += sumD_adjacent;
				sumT += sumT_adjacent;
				msdPerMoleculeAllVsAll.add(px2ToUm2PerSecond * sumD / sumT);
				msdPerMoleculeAdjacent.add(px2ToUm2PerSecond * sumD_adjacent / sumT_adjacent);
			}

			StoredDataStatistics dStarPerMoleculeAllVsAll = null;
			StoredDataStatistics dStarPerMoleculeAdjacent = null;
			calculatePrecision(traces);
			if (saveTraceDistances || (settings.showHistograms && displayDHistogram))
			{
				dStarPerMoleculeAllVsAll = calculateDiffusionCoefficient(msdPerMoleculeAllVsAll);
				dStarPerMoleculeAdjacent = calculateDiffusionCoefficient(msdPerMoleculeAdjacent);
			}

			if (saveTraceDistances)
			{
				saveTraceDistances(traces.length, distances, msdPerMoleculeAllVsAll, msdPerMoleculeAdjacent,
						dStarPerMoleculeAllVsAll, dStarPerMoleculeAdjacent);
			}

			// Calculate the cumulative jump-distance histogram
			double[][] jdHistogram = Maths.cumulativeHistogram(jumpDistances.getValues(), true);
			// Always show the jump distance histogram
			String jdTitle = TITLE + " Jump Distance";
			Plot jdPlot = new Plot(jdTitle, "Distance (um^2/second)", "Cumulative Probability", jdHistogram[0],
					jdHistogram[1]);
			Utils.display(jdTitle, jdPlot);

			// Plot the per-trace histogram of MSD and D*
			if (settings.showHistograms)
			{
				if (displayMSDHistogram)
				{
					Utils.showHistogram(TITLE, msdPerMoleculeAllVsAll, "MSD/Molecule (all-vs-all)", 0,
							(settings.removeOutliers) ? 1 : 0, settings.histogramBins);
					Utils.showHistogram(TITLE, msdPerMoleculeAdjacent, "MSD/Molecule (adjacent)", 0,
							(settings.removeOutliers) ? 1 : 0, settings.histogramBins);
				}
				if (displayDHistogram)
				{
					Utils.showHistogram(TITLE, dStarPerMoleculeAllVsAll, "D*/Molecule (all-vs-all)", 0,
							(settings.removeOutliers) ? 1 : 0, settings.histogramBins);
					Utils.showHistogram(TITLE, dStarPerMoleculeAdjacent, "D*/Molecule (adjacent)", 0,
							(settings.removeOutliers) ? 1 : 0, settings.histogramBins);
				}
			}

			if (subSample)
			{
				// Extract a representative subset so that the sample size for all lengths is equal
				int size = Integer.MAX_VALUE;
				for (int i = stats.length - 1; i-- > 0;)
					if (stats[i].getN() > 0)
						size = FastMath.min(size, stats[i].getN());
				if (size < Integer.MAX_VALUE)
				{
					Random random = new Random();
					for (int i = stats.length - 1; i-- > 0;)
					{
						double[] data = ((StoredDataStatistics) stats[i]).getValues();
						random.shuffle(data);
						stats[i] = new Statistics(Arrays.copyOf(data, size));
					}
				}
			}

			// Calculate the mean squared distance (MSD)
			double[] x = new double[stats.length];
			double[] y = new double[x.length];
			double[] sd = new double[x.length];
			for (int i = 1; i < stats.length; i++)
			{
				x[i] = i * exposureTime;
				y[i] = stats[i].getMean() * px2ToUm2;
				//sd[i] = stats[i].getStandardDeviation() * px2ToUm2;
				sd[i] = stats[i].getStandardError() * px2ToUm2;
			}

			String title = TITLE + " MSD";
			Plot plot = plotMSD(x, y, sd, title);

			// Fit the MSD using a linear fit that must pass through 0,0
			D = fitMSD(x, y, title, plot);

			// Fit Jump Distance cumulative probability
			jdParams = fitJumpDistance(jumpDistances, jdHistogram, jdTitle, jdPlot);
		}

		summarise(traces, D, jdParams);

		IJ.showStatus(String.format("%d localisations => %d traces", results.size(), count));
	}

	/**
	 * Calculate the average precision of localisation in the traces
	 * 
	 * @param traces
	 */
	private void calculatePrecision(Trace[] traces)
	{
		// Get the average precision of the localisations
		precision = 0;
		final double nmPerPixel = results.getNmPerPixel();
		final double gain = results.getGain();
		final boolean emCCD = results.isEMCCD();
		int n = 0;
		for (Trace trace : traces)
		{
			for (PeakResult r : trace.getPoints())
				precision += r.getPrecision(nmPerPixel, gain, emCCD);
			n += trace.size();
		}
		precision /= n;

		if (precision > 100)
		{
			GenericDialog gd = new GenericDialog(TITLE);
			gd.addMessage("The average precision of the traced results is " + Utils.rounded(precision, 4) +
					" nm.\nPlease verify the precision.");
			gd.addSlider("Precision (nm)", 5, 100, precision);
			gd.showDialog();
			if (!(gd.wasCanceled() || gd.invalidNumber()))
			{
				precision = Math.abs(gd.getNextNumber());
			}
		}
	}

	/**
	 * Calculate the apparent diffusion coefficient of the molecule. This is done by using the mean-squared deviation
	 * between frames divided by the time interval (delta) between frames. This is divided by 4 to produce the diffusion
	 * coefficient from 2D distance analysis. The apparent coefficient is produced by subtracting the squared
	 * localisation error also normalised by the time delta.
	 * <p>
	 * See Uphoff, et al, 2013. Single-molecule DNA repair in live bacteria, PNAS 110, 8063-8068
	 * 
	 * @param msdPerMoleculeAdjacent
	 * @return
	 */
	private StoredDataStatistics calculateDiffusionCoefficient(StoredDataStatistics msdPerMoleculeAdjacent)
	{
		StoredDataStatistics dStarPerMolecule = new StoredDataStatistics();
		final double diffusionCoefficientConversion = 1.0 / 4.0;
		// convert precision to um
		final double error = precision * precision / (1e6 * exposureTime);
		for (double msd : msdPerMoleculeAdjacent.getValues())
		{
			dStarPerMolecule.add(FastMath.max(0, msd * diffusionCoefficientConversion - error));
		}
		return dStarPerMolecule;
	}

	private void saveTraceDistances(int nTraces, ArrayList<double[]> distances, StoredDataStatistics msdPerMolecule,
			StoredDataStatistics msdPerMoleculeAdjacent, StoredDataStatistics dStarPerMolecule,
			StoredDataStatistics dStarPerMoleculeAdjacent)
	{
		String[] path = Utils.decodePath(filename);
		OpenDialog chooser = new OpenDialog("Trace_distances_File", path[0], path[1]);
		if (chooser.getFileName() != null)
		{
			filename = Utils.replaceExtension(chooser.getDirectory() + chooser.getFileName(), "xls");

			OutputStreamWriter out = null;
			try
			{
				FileOutputStream fos = new FileOutputStream(filename);
				out = new OutputStreamWriter(fos, "UTF-8");
				double[] msd = msdPerMolecule.getValues();
				double[] msd2 = msdPerMoleculeAdjacent.getValues();
				double[] dStar = dStarPerMolecule.getValues();
				double[] dStar2 = dStarPerMoleculeAdjacent.getValues();
				out.write(String.format("#%d traces : Precision = %s nm : Exposure time = %s s\n", nTraces,
						Utils.rounded(precision, 4), Utils.rounded(exposureTime, 4)));
				out.write(String
						.format("#TraceId\tMSD all-vs-all (um^2/s)\tMSD adjacent (um^2/s)\tD* all-vs-all(um^2/s)\tD* adjacent(um^2/s)\tDistances (um^2) per %ss ... \n",
								Utils.rounded(exposureTime, 4)));
				for (int i = 0; i < msd.length; i++)
				{
					out.write(Integer.toString(i + 1));
					out.write('\t');
					out.write(Utils.rounded(msd[i], 4));
					out.write('\t');
					out.write(Utils.rounded(msd2[i], 4));
					out.write('\t');
					out.write(Utils.rounded(dStar[i], 4));
					out.write('\t');
					out.write(Utils.rounded(dStar2[i], 4));
					for (double d : distances.get(i))
					{
						out.write('\t');
						out.write(Utils.rounded(d, 4));
					}
					out.write('\n');
				}
			}
			catch (Exception e)
			{
			}
			finally
			{
				if (out != null)
				{
					try
					{
						out.close();
					}
					catch (IOException e)
					{
					}
				}
			}
		}
	}

	private double distance2(final float x, final float y, PeakResult r2)
	{
		final double dx = x - r2.getXPosition();
		final double dy = y - r2.getYPosition();
		return dx * dx + dy * dy;
	}

	/**
	 * Filter traces that are not the minimum length
	 * 
	 * @param traces
	 * @param minimumTraceLength
	 * @return The new traces
	 */
	private Trace[] filterTraces(Trace[] traces, int minimumTraceLength)
	{
		int count = 0;
		for (int i = 0; i < traces.length; i++)
		{
			if (traces[i].size() >= minimumTraceLength)
				traces[count++] = traces[i];
		}

		Utils.log("%d Traces filtered to %d using minimum length %d", traces.length, count, minimumTraceLength);
		return Arrays.copyOf(traces, count);
	}

	private String createSettingsComment()
	{
		return String.format("Molecule tracing : distance-threshold = %f nm", settings.distanceThreshold);
	}

	private void summarise(Trace[] traces, double D, double[][] jdParams)
	{
		IJ.showStatus("Calculating summary ...");

		// Create summary table
		createSummaryTable();

		Statistics[] stats = new Statistics[NAMES.length];
		for (int i = 0; i < stats.length; i++)
		{
			stats[i] = (settings.showHistograms) ? new StoredDataStatistics() : new Statistics();
		}
		for (Trace trace : traces)
		{
			stats[T_ON].add(trace.getOnTime() * exposureTime);
			final double signal = trace.getSignal() / results.getGain();
			stats[TOTAL_SIGNAL].add(signal);
			stats[SIGNAL_PER_FRAME].add(signal / trace.size());
		}

		// Add to the summary table
		StringBuilder sb = new StringBuilder();
		sb.append(results.getName()).append("\t");
		sb.append(Utils.rounded(exposureTime * 1000, 3)).append("\t");
		sb.append(Utils.rounded(settings.distanceThreshold, 3)).append("\t");
		sb.append(Utils.rounded(settings.distanceExclusion, 3)).append("\t");
		sb.append(settings.minimumTraceLength).append("\t");
		sb.append(settings.truncate).append("\t");
		sb.append(settings.internalDistances).append("\t");
		sb.append(settings.fitLength).append("\t");
		sb.append(traces.length).append("\t");
		sb.append(Utils.rounded(D, 4)).append("\t");
		sb.append(Utils.rounded(settings.jumpDistance * exposureTime)).append("\t");
		if (jdParams == null)
		{
			sb.append("\t\t");
		}
		else
		{
			sb.append(format(jdParams[0])).append("\t");
			sb.append(format(jdParams[1])).append("\t");
		}

		for (int i = 0; i < stats.length; i++)
		{
			sb.append(Utils.rounded(stats[i].getMean(), 3)).append("\t");
		}
		if (java.awt.GraphicsEnvironment.isHeadless())
		{
			IJ.log(sb.toString());
			return;
		}
		else
		{
			summaryTable.append(sb.toString());
		}

		if (settings.showHistograms)
		{
			IJ.showStatus("Calculating histograms ...");

			int[] idList = new int[NAMES.length];
			int count = 0;

			boolean requireRetile = false;
			for (int i = 0; i < NAMES.length; i++)
			{
				if (displayHistograms[i])
				{
					idList[count++] = Utils.showHistogram(TITLE, (StoredDataStatistics) stats[i], NAMES[i], 0,
							(settings.removeOutliers || alwaysRemoveOutliers[i]) ? 2 : 0, settings.histogramBins);
					requireRetile = requireRetile || Utils.isNewWindow();
				}
			}

			if (count > 0 && requireRetile)
			{
				idList = Arrays.copyOf(idList, count);
				new WindowOrganiser().tileWindows(idList);
			}
		}

		IJ.showStatus("");
	}

	private String format(double[] jumpD)
	{
		if (jumpD == null || jumpD.length == 0)
			return "";
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < jumpD.length; i++)
		{
			if (i != 0)
				sb.append(", ");
			sb.append(Utils.rounded(jumpD[i], 4));
		}
		return sb.toString();
	}

	private void createSummaryTable()
	{
		if (java.awt.GraphicsEnvironment.isHeadless())
		{
			if (header == null)
			{
				header = createHeader();
				IJ.log(header);
			}
		}
		else
		{
			if (summaryTable == null || !summaryTable.isVisible())
			{
				summaryTable = new TextWindow(TITLE + " Data Summary", createHeader(), "", 800, 300);
				summaryTable.setVisible(true);
			}
		}
	}

	private String createHeader()
	{
		StringBuilder sb = new StringBuilder(
				"Dataset\tExposure time (ms)\tD-threshold (nm)\tEx-threshold (nm)\tMin.Length\tTruncate\tInternal\tFit Length\tTraces\tD (um^2/s)\tJump Distance (s)\tJump D (um^2/s)\tFractions");
		for (int i = 0; i < NAMES.length; i++)
		{
			sb.append("\t").append(NAMES[i]);
		}
		return sb.toString();
	}

	private boolean showDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);

		globalSettings = SettingsManager.loadSettings();
		settings = globalSettings.getClusteringSettings();

		gd.addNumericField("Distance_Threshold (nm)", settings.distanceThreshold, 0);
		gd.addNumericField("Distance_Exclusion (nm)", settings.distanceExclusion, 0);
		gd.addSlider("Min_trace_length", 2, 20, settings.minimumTraceLength);
		gd.addCheckbox("Truncate_traces", settings.truncate);
		gd.addCheckbox("Internal_distances", settings.internalDistances);
		//gd.addCheckbox("Sub-sample_distances", settings.subSampledDistances);
		gd.addSlider("Fit_length", 2, 20, settings.fitLength);
		gd.addSlider("Jump_distance", 1, 20, settings.jumpDistance);
		gd.addSlider("Minimum_difference", 0, 10, minDifference);
		gd.addSlider("Minimum_fraction", 0, 1, minFraction);
		gd.addCheckbox("Save_traces", settings.saveTraces);
		gd.addCheckbox("Save_trace_distances", saveTraceDistances);
		gd.addCheckbox("Show_histograms", settings.showHistograms);

		gd.showDialog();

		if (gd.wasCanceled() || !readDialog(gd))
			return false;

		// Update the settings
		SettingsManager.saveSettings(globalSettings);

		// Load the results
		results = ResultsManager.loadInputResults(inputOption, true);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			IJ.showStatus("");
			return false;
		}

		// Check the results have a calibrated exposure time
		if (results.getCalibration() == null || results.getCalibration().exposureTime <= 0)
		{
			gd = new GenericDialog(TITLE);
			gd.addMessage("Uncalibrated results. Please enter the exposure time for each frame:");
			gd.addNumericField("Exposure_time (ms)", 100, 2);
			gd.showDialog();
			if (gd.wasCanceled() || gd.invalidNumber())
				return false;
			exposureTime = Math.abs(gd.getNextNumber() / 1000);
		}
		else
		{
			exposureTime = results.getCalibration().exposureTime / 1000;
		}

		return true;
	}

	private boolean readDialog(GenericDialog gd)
	{
		inputOption = ResultsManager.getInputSource(gd);
		settings.distanceThreshold = gd.getNextNumber();
		settings.distanceExclusion = Math.abs(gd.getNextNumber());
		settings.minimumTraceLength = (int) Math.abs(gd.getNextNumber());
		settings.truncate = gd.getNextBoolean();
		settings.internalDistances = gd.getNextBoolean();
		//settings.subSampledDistances = gd.getNextBoolean();
		settings.fitLength = (int) Math.abs(gd.getNextNumber());
		settings.jumpDistance = (int) Math.abs(gd.getNextNumber());
		minDifference = Math.abs(gd.getNextNumber());
		minFraction = Math.abs(gd.getNextNumber());
		settings.saveTraces = gd.getNextBoolean();
		saveTraceDistances = gd.getNextBoolean();
		settings.showHistograms = gd.getNextBoolean();

		if (gd.invalidNumber())
			return false;

		if (settings.showHistograms)
		{
			gd = new GenericDialog(TITLE);
			gd.addMessage("Select the histograms to display");
			gd.addCheckbox("Remove_outliers", settings.removeOutliers);
			gd.addNumericField("Histogram_bins", settings.histogramBins, 0);
			for (int i = 0; i < displayHistograms.length; i++)
				gd.addCheckbox(NAMES[i].replace(' ', '_'), displayHistograms[i]);
			gd.addCheckbox("MSD/Molecule", displayMSDHistogram);
			gd.addCheckbox("D*/Molecule", displayDHistogram);
			gd.showDialog();
			if (gd.wasCanceled())
				return false;
			settings.removeOutliers = gd.getNextBoolean();
			settings.histogramBins = (int) Math.abs(gd.getNextNumber());
			for (int i = 0; i < displayHistograms.length; i++)
				displayHistograms[i] = gd.getNextBoolean();
			displayMSDHistogram = gd.getNextBoolean();
			displayDHistogram = gd.getNextBoolean();
		}

		// Check arguments
		try
		{
			Parameters.isAboveZero("Distance threshold", settings.distanceThreshold);
			Parameters.isAbove("Min trace length", settings.minimumTraceLength, 1);
			Parameters.isAboveZero("Histogram bins", settings.histogramBins);
			Parameters.isAbove("Fit length", settings.fitLength, 1);
			Parameters.isAboveZero("Jump distance", settings.jumpDistance);
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		return true;
	}

	private Plot plotMSD(double[] x, double[] y, double[] sd, String title)
	{
		Plot plot = new Plot(title, "Time (s)", "Distance (um^2)", x, y);
		// Set limits before any plotting
		double max = 0;
		for (int i = 1; i < x.length; i++)
		{
			double value = y[i] + sd[i];
			max = FastMath.max(max, value);
		}
		plot.setLimits(0, x[x.length - 1] + exposureTime * 0.5, 0, max);
		plot.setColor(Color.blue);
		for (int i = 1; i < x.length; i++)
		{
			plot.drawLine(x[i], y[i] - sd[i], x[i], y[i] + sd[i]);
		}
		plot.setColor(Color.red);
		Utils.display(title, plot);
		return plot;
	}

	/**
	 * Fit the MSD using a linear fit that must pass through 0,0.
	 * <p>
	 * Update the plot by adding the fit line.
	 * 
	 * @param x
	 * @param y
	 * @param title
	 * @param plot
	 * @return
	 */
	private double fitMSD(double[] x, double[] y, String title, Plot plot)
	{
		double D = 0;
		LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer();
		PointVectorValuePair lvmSolution;
		try
		{
			final LinearFunction function = new LinearFunction(x, y, settings.fitLength);
			double[] parameters = new double[] { function.guess() };
			lvmSolution = optimizer.optimize(new MaxIter(3000), new MaxEval(Integer.MAX_VALUE),
					new ModelFunctionJacobian(new MultivariateMatrixFunction()
					{
						@Override
						public double[][] value(double[] point) throws IllegalArgumentException
						{
							return function.jacobian(point);
						}
					}), new ModelFunction(function), new Target(function.getY()), new Weight(function.getWeights()),
					new InitialGuess(parameters));

			double ss = 0;
			double[] obs = function.getY();
			double[] exp = lvmSolution.getValue();
			for (int i = 0; i < obs.length; i++)
				ss += (obs[i] - exp[i]) * (obs[i] - exp[i]);

			D = lvmSolution.getPoint()[0] / 4;
			Utils.log("Linear fit (%d points) : Gradient = %s, D = %s um^2/s, SS = %f (%d evaluations)", obs.length,
					Utils.rounded(lvmSolution.getPoint()[0], 4), Utils.rounded(D, 4), ss, optimizer.getEvaluations());

			// Add the fit to the plot
			plot.setColor(Color.magenta);
			plot.drawLine(0, 0, x[x.length - 1], x[x.length - 1] * 4 * D);
			Utils.display(title, plot);
		}
		catch (TooManyIterationsException e)
		{
			Utils.log("Failed to fit : Too many iterations (%d)", optimizer.getIterations());
		}
		catch (ConvergenceException e)
		{
			Utils.log("Failed to fit : %s", e.getMessage());
		}
		return D;
	}

	public class LinearFunction implements MultivariateVectorFunction
	{
		double[] x, y;

		public LinearFunction(double[] x, double[] y, int length)
		{
			int to = FastMath.min(x.length, 1 + length);
			this.x = Arrays.copyOfRange(x, 1, to);
			this.y = Arrays.copyOfRange(y, 1, to);
		}

		// Adapted from http://commons.apache.org/proper/commons-math/userguide/optimization.html
		// Use the deprecated API since the new one is not yet documented.

		/**
		 * @return An estimate for the linear gradient
		 */
		public double guess()
		{
			return y[y.length - 1] / x[x.length - 1];
		}

		public double[] getWeights()
		{
			double[] w = new double[x.length];
			Arrays.fill(w, 1);
			return w;
		}

		public double[] getY()
		{
			return y;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see org.apache.commons.math3.analysis.MultivariateVectorFunction#value(double[])
		 */
		public double[] value(double[] variables)
		{
			double[] values = new double[x.length];
			for (int i = 0; i < values.length; i++)
			{
				values[i] = x[i] * variables[0];
			}
			return values;
		}

		double[][] jacobian(double[] variables)
		{
			// Compute the gradients using calculus differentiation:
			// y = ax
			// y' = x
			double[][] jacobian = new double[x.length][variables.length];

			for (int i = 0; i < jacobian.length; ++i)
			{
				jacobian[i][0] = x[i];
			}

			return jacobian;
		}
	}

	/**
	 * Fit the jump distance histogram using a cumulative sum as detailed in
	 * <p>
	 * Update the plot by adding the fit line.
	 * 
	 * @param jumpDistances
	 * @param jdHistogram
	 * @param title
	 * @param plot
	 * @return
	 */
	private double[][] fitJumpDistance(StoredDataStatistics jumpDistances, double[][] jdHistogram, String title,
			Plot plot)
	{
		final double meanDistance = Math.sqrt(jumpDistances.getMean()) * 1e3;
		final double beta = meanDistance / precision;
		Utils.log(
				"Jump Distance analysis : N = %d, Time = %d frames (%s seconds). Mean Distance = %s nm, Precision = %s nm, Beta = %s",
				jumpDistances.getN(), settings.jumpDistance, Utils.rounded(settings.jumpDistance * exposureTime, 4),
				Utils.rounded(meanDistance, 4), Utils.rounded(precision, 4), Utils.rounded(beta, 4));
		int n = 0;
		int N = 10;
		double[] SS = new double[N];
		Arrays.fill(SS, -1);
		double[] ic = new double[N];
		Arrays.fill(ic, Double.POSITIVE_INFINITY);
		double[][] coefficients = new double[N][];
		double[][] fractions = new double[N][];
		double[][] fitParams = new double[N][];
		double bestIC = Double.POSITIVE_INFINITY;
		int best = -1;

		// Guess the D
		final double estimatedD = jumpDistances.getMean() / 4;
		Utils.log("Estimated D = %s um^2/s", Utils.rounded(estimatedD, 4));

		// Fit using a single population model
		LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer();
		try
		{
			final JumpDistanceFunction function = new JumpDistanceFunction(jdHistogram[0], jdHistogram[1], estimatedD);
			PointVectorValuePair lvmSolution = optimizer.optimize(new MaxIter(3000), new MaxEval(Integer.MAX_VALUE),
					new ModelFunctionJacobian(new MultivariateMatrixFunction()
					{
						@Override
						public double[][] value(double[] point) throws IllegalArgumentException
						{
							return function.jacobian(point);
						}
					}), new ModelFunction(function), new Target(function.getY()), new Weight(function.getWeights()),
					new InitialGuess(function.guess()));

			fitParams[n] = lvmSolution.getPointRef();
			SS[n] = calculateSumOfSquares(function.getY(), lvmSolution.getValueRef());
			ic[n] = Maths.getInformationCriterion(SS[n], function.x.length, 1);
			coefficients[n] = fitParams[n];
			fractions[n] = new double[] { 1 };

			Utils.log("Fit Jump distance (N=%d) : D = %s um^2/s, SS = %f, IC = %s (%d evaluations)", n + 1,
					Utils.rounded(fitParams[n][0], 4), SS[n], Utils.rounded(ic[n], 4), optimizer.getEvaluations());

			bestIC = ic[n];
			best = 0;

			addToPlot(function, fitParams[n], jdHistogram, title, plot, Color.magenta);
		}
		catch (TooManyIterationsException e)
		{
			Utils.log("Failed to fit : Too many iterations (%d)", optimizer.getIterations());
		}
		catch (ConvergenceException e)
		{
			Utils.log("Failed to fit : %s", e.getMessage());
		}

		n++;

		// Fit using a mixed population model. 
		// Vary n from 2 to N. Stop when the fit fails or the fit is worse.
		int bestMulti = -1;
		double bestMultiIC = Double.POSITIVE_INFINITY;
		while (n < N)
		{
			// Uses a weighted sum of n exponential functions, each function models a fraction of the particles.
			// An LVM fit cannot restrict the parameters so the fractions do not go below zero.
			// Use the CMEASOptimizer which supports bounded fitting.

			MixedJumpDistanceFunctionMultivariate mixedFunction = new MixedJumpDistanceFunctionMultivariate(
					jdHistogram[0], jdHistogram[1], estimatedD, n + 1);

			double[] lB = mixedFunction.getLowerBounds();
			double[] uB = mixedFunction.getUpperBounds();
			SimpleBounds bounds = new SimpleBounds(lB, uB);

			int maxIterations = 2000;
			double stopFitness = 0; //Double.NEGATIVE_INFINITY;
			boolean isActiveCMA = true;
			int diagonalOnly = 20;
			int checkFeasableCount = 1;
			RandomGenerator random = new Well19937c();
			boolean generateStatistics = false;
			ConvergenceChecker<PointValuePair> checker = new SimpleValueChecker(1e-6, 1e-10);
			// The sigma determines the search range for the variables. It should be 1/3 of the initial search region.
			double[] s = new double[lB.length];
			for (int i = 0; i < s.length; i++)
				s[i] = (uB[i] - lB[i]) / 3;
			OptimizationData sigma = new CMAESOptimizer.Sigma(s);
			OptimizationData popSize = new CMAESOptimizer.PopulationSize((int) (4 + Math.floor(3 * Math
					.log(mixedFunction.x.length))));

			CMAESOptimizer opt = null;
			try
			{
				opt = new CMAESOptimizer(maxIterations, stopFitness, isActiveCMA, diagonalOnly, checkFeasableCount,
						random, generateStatistics, checker);
				PointValuePair constrainedSolution = opt.optimize(new InitialGuess(mixedFunction.guess()),
						new ObjectiveFunction(mixedFunction), GoalType.MINIMIZE, bounds, sigma, popSize, new MaxIter(
								maxIterations), new MaxEval(maxIterations * 2));

				int evaluations = opt.getEvaluations();
				fitParams[n] = constrainedSolution.getPointRef();
				SS[n] = constrainedSolution.getValue();

				// Try and improve using a LVM fit
				final MixedJumpDistanceFunctionGradient mixedFunctionGradient = new MixedJumpDistanceFunctionGradient(
						jdHistogram[0], jdHistogram[1], estimatedD, n + 1);

				PointVectorValuePair lvmSolution;
				try
				{
					lvmSolution = optimizer.optimize(new MaxIter(3000), new MaxEval(Integer.MAX_VALUE),
							new ModelFunctionJacobian(new MultivariateMatrixFunction()
							{
								@Override
								public double[][] value(double[] point) throws IllegalArgumentException
								{
									return mixedFunctionGradient.jacobian(point);
								}
							}), new ModelFunction(mixedFunctionGradient), new Target(mixedFunctionGradient.getY()),
							new Weight(mixedFunctionGradient.getWeights()), new InitialGuess(fitParams[n]));
					double ss = calculateSumOfSquares(mixedFunctionGradient.getY(), lvmSolution.getValue());
					// All fitted parameters must be above zero
					if (ss < SS[n] && Maths.min(lvmSolution.getPoint()) > 0)
					{
						//Utils.log("  Re-fitting improved the SS from %s to %s (-%s%%)", Utils.rounded(SS[n], 4),
						//		Utils.rounded(ss, 4), Utils.rounded(100 * (SS[n] - ss) / SS[n], 4));
						fitParams[n] = lvmSolution.getPoint();
						SS[n] = ss;
						evaluations += optimizer.getEvaluations();
					}
				}
				catch (TooManyIterationsException e)
				{
					//Utils.log("Failed to re-fit : Too many evaluations (%d)", optimizer.getEvaluations());
				}
				catch (ConvergenceException e)
				{
					//Utils.log("Failed to re-fit : %s", e.getMessage());
				}

				// Since the fractions must sum to one we subtract 1 degree of freedom from the number of parameters
				ic[n] = Maths.getInformationCriterion(SS[n], mixedFunction.x.length, fitParams[n].length - 1);

				double[] d = new double[n + 1];
				double[] f = new double[n + 1];
				double sum = 0;
				for (int i = 0; i < d.length; i++)
				{
					f[i] = fitParams[n][i * 2];
					sum += f[i];
					d[i] = fitParams[n][i * 2 + 1];
				}
				for (int i = 0; i < f.length; i++)
					f[i] /= sum;
				// Sort by coefficient size
				sort(d, f);
				coefficients[n] = d;
				fractions[n] = f;

				Utils.log("Fit Jump distance (N=%d) : D = %s um^2/s (%s), SS = %f, IC = %s (%d evaluations)", n + 1,
						format(d), format(f), SS[n], Utils.rounded(ic[n], 4), evaluations);

				boolean valid = true;
				for (int i = 0; i < f.length; i++)
				{
					// Check the fit has fractions above the minimum fraction
					if (f[i] < minFraction)
					{
						Utils.log("Fraction is less than the minimum fraction: %s < %s", Utils.rounded(f[i]),
								Utils.rounded(minFraction));
						valid = false;
						break;
					}
					// Check the coefficients are different
					if (i + 1 < f.length && d[i] / d[i + 1] < minDifference)
					{
						Utils.log("Coefficients are not different: %s / %s = %s < %s", Utils.rounded(d[i]),
								Utils.rounded(d[i + 1]), Utils.rounded(d[i] / d[i + 1]), Utils.rounded(minDifference));
						valid = false;
						break;
					}
				}

				if (!valid)
					break;

				// Store the best model
				if (bestIC > ic[n])
				{
					bestIC = ic[n];
					best = n;
				}

				// Store the best multi model
				if (bestMultiIC < ic[n])
				{
					break;
				}

				bestMultiIC = ic[n];
				bestMulti = n;
			}
			catch (TooManyEvaluationsException e)
			{
				Utils.log("Failed to fit : Too many evaluations (%d)", opt.getEvaluations());
			}

			n++;
		}

		// Add the best fit to the plot and return the parameters.
		if (bestMulti > -1)
		{
			Function function = new MixedJumpDistanceFunctionMultivariate(jdHistogram[0], jdHistogram[1], 0,
					bestMulti + 1);
			addToPlot(function, fitParams[bestMulti], jdHistogram, title, plot, Color.yellow);
		}

		if (best > -1)
		{
			Utils.log("Best fit achieved using %d population%s: D = %s um^2/s, Fractions = %s", best + 1,
					(best == 0) ? "" : "s", format(coefficients[best]), format(fractions[best]));
		}

		return (best > -1) ? new double[][] { coefficients[best], fractions[best] } : null;
	}

	private void sort(double[] d, double[] f)
	{
		// Sort by coefficient size
		int[] indices = new int[f.length];
		for (int i = 0; i < f.length; i++)
		{
			indices[i] = i;
		}
		Sort.sort(indices, d);
		double[] d2 = Arrays.copyOf(d, d.length);
		double[] f2 = Arrays.copyOf(f, f.length);
		for (int i = 0; i < f.length; i++)
		{
			d[i] = d2[indices[i]];
			f[i] = f2[indices[i]];
		}
	}

	private void addToPlot(Function function, double[] params, double[][] jdHistogram, String title, Plot plot,
			Color color)
	{
		final double max = jdHistogram[0][jdHistogram[0].length - 1];
		final int nPoints = 300;
		final double interval = max / nPoints;
		double[] x = new double[nPoints + 1];
		double[] y = new double[nPoints + 1];

		for (int i = 0; i < nPoints; i++)
		{
			x[i] = i * interval;
			y[i] = function.evaluate(x[i], params);
		}
		x[nPoints] = max;
		y[nPoints] = function.evaluate(max, params);

		plot.setColor(color);
		plot.addPoints(x, y, Plot.LINE);
		Utils.display(title, plot);
	}

	private double calculateSumOfSquares(double[] obs, double[] exp)
	{
		double ss = 0;
		for (int i = 0; i < obs.length; i++)
			ss += (obs[i] - exp[i]) * (obs[i] - exp[i]);
		return ss;
	}

	public abstract class Function
	{
		double[] x, y;
		double estimatedD;

		public Function(double[] x, double[] y, double estimatedD)
		{
			this.x = x;
			this.y = y;
			this.estimatedD = estimatedD;
		}

		/**
		 * @return An estimate for the parameters
		 */
		public abstract double[] guess();

		public double[] getWeights()
		{
			double[] w = new double[x.length];
			Arrays.fill(w, 1);
			return w;
		}

		public double[] getX()
		{
			return x;
		}

		public double[] getY()
		{
			return y;
		}

		public abstract double evaluate(double x, double[] parameters);

		public double[][] jacobian(double[] variables)
		{
			double[][] jacobian = new double[x.length][variables.length];

			final double delta = 0.001;
			double[] d = new double[variables.length];
			double[][] variables2 = new double[variables.length][];
			for (int i = 0; i < variables.length; i++)
			{
				d[i] = delta * Math.abs(variables[i]); // Should the delta be changed for each parameter ?
				variables2[i] = Arrays.copyOf(variables, variables.length);
				variables2[i][i] += d[i];
			}
			for (int i = 0; i < jacobian.length; ++i)
			{
				double value = evaluate(x[i], variables);
				for (int j = 0; j < variables.length; j++)
				{
					double value2 = evaluate(x[i], variables2[j]);
					jacobian[i][j] = (value2 - value) / d[j];
				}
			}
			return jacobian;
		}
	}

	public class JumpDistanceFunction extends Function implements MultivariateVectorFunction
	{
		public JumpDistanceFunction(double[] x, double[] y, double estimatedD)
		{
			super(x, y, estimatedD);
		}

		// Adapted from http://commons.apache.org/proper/commons-math/userguide/optimization.html
		// Use the deprecated API since the new one is not yet documented.

		public double[] guess()
		{
			return new double[] { estimatedD };
		}

		public double evaluate(double x, double[] params)
		{
			return 1 - FastMath.exp(-x / (4 * params[0]));
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see org.apache.commons.math3.analysis.MultivariateVectorFunction#value(double[])
		 */
		public double[] value(double[] variables)
		{
			double[] values = new double[x.length];
			final double fourD = 4 * variables[0];
			for (int i = 0; i < values.length; i++)
			{
				values[i] = 1 - FastMath.exp(-x[i] / fourD);
			}
			return values;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see gdsc.smlm.ij.plugins.TraceDiffusion.Function#jacobian(double[])
		 */
		public double[][] jacobian(double[] variables)
		{
			// Compute the gradients using calculus differentiation:
			// y = 1 - a
			// a = exp(b)
			// b = -x / 4D
			//
			// y' = -a'
			// a' = exp(b) * b'
			// b' = -1 * -x / 4D^2 = x / 4D^2
			// y' = -exp(b) * x / 4D^2
			//    = -a * -b / D
			//    = a * b / D
			//    = exp(b) * b / D

			final double d = variables[0];
			final double fourD = 4 * d;
			double[][] jacobian = new double[x.length][variables.length];

			for (int i = 0; i < jacobian.length; ++i)
			{
				final double b = -x[i] / fourD;
				jacobian[i][0] = FastMath.exp(b) * b / d;
			}

			//// Check numerically ...
			//double[][] jacobian2 = super.jacobian(variables);
			//for (int i = 0; i < jacobian.length; i++)
			//{
			//	System.out.printf("dD = %g : %g = %g\n", jacobian[i][0], jacobian2[i][0],
			//			DoubleEquality.relativeError(jacobian[i][0], jacobian2[i][0]));
			//}

			return jacobian;
		}
	}

	public class MixedJumpDistanceFunction extends Function
	{
		int n;

		public MixedJumpDistanceFunction(double[] x, double[] y, double estimatedD, int n)
		{
			super(x, y, estimatedD);
			this.n = n;
		}

		public double[] guess()
		{
			// Store the fraction and then the diffusion coefficient.
			// Q. Should this be modified to set one fraction to always be 1? 
			// Having an actual parameter for fitting will allow the optimisation engine to 
			// adjust the fraction for its diffusion coefficient relative to the others.
			double[] guess = new double[n * 2];
			double d = estimatedD;
			for (int i = 0; i < n; i++)
			{
				// Fraction are all equal
				guess[i * 2] = 1;
				// Diffusion coefficient gets smaller for each fraction
				guess[i * 2 + 1] = d;
				d *= 0.1;
			}
			return guess;
		}

		public double[] getUpperBounds()
		{
			double[] bounds = new double[n * 2];
			for (int i = 0; i < n; i++)
			{
				// Fraction guess is 1 so set the upper limit as 10
				bounds[i * 2] = 10;
				// Diffusion coefficient could be 10x the estimated
				bounds[i * 2 + 1] = estimatedD * 10;
			}
			return bounds;
		}

		public double[] getLowerBounds()
		{
			return new double[n * 2];
		}

		public double evaluate(double x, double[] params)
		{
			double sum = 0;
			double total = 0;
			for (int i = 0; i < n; i++)
			{
				final double f = params[i * 2];
				sum += f * FastMath.exp(-x / (4 * params[i * 2 + 1]));
				total += f;
			}
			return 1 - sum / total;
		}

		public double[] getValue(double[] variables)
		{
			double total = 0;
			for (int i = 0; i < n; i++)
			{
				total += variables[i * 2];
			}

			final double[] fourD = new double[n];
			final double[] f = new double[n];
			for (int i = 0; i < n; i++)
			{
				f[i] = variables[i * 2] / total;
				fourD[i] = 4 * variables[i * 2 + 1];
			}

			double[] values = new double[x.length];
			for (int i = 0; i < values.length; i++)
			{
				double sum = 0;
				for (int j = 0; j < n; j++)
				{
					sum += f[j] * FastMath.exp(-x[i] / (fourD[j]));
				}
				values[i] = 1 - sum;
			}
			return values;
		}
	}

	public class MixedJumpDistanceFunctionGradient extends MixedJumpDistanceFunction implements
			MultivariateVectorFunction
	{
		public MixedJumpDistanceFunctionGradient(double[] x, double[] y, double estimatedD, int n)
		{
			super(x, y, estimatedD, n);
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see org.apache.commons.math3.analysis.MultivariateVectorFunction#value(double[])
		 */
		public double[] value(double[] point) throws IllegalArgumentException
		{
			return getValue(point);
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see gdsc.smlm.ij.plugins.TraceDiffusion.Function#jacobian(double[])
		 */
		public double[][] jacobian(double[] variables)
		{
			// Compute the gradients using calculus differentiation:
			// y = 1 - sum(a)
			// The sum is over n components of the following function
			// a = f * exp(b)
			// b = -x / 4D
			// Each function contributes a fraction f:
			// f = fj / sum_j(f)

			// The gradient is the sum of the individual gradients. The diffusion coefficient is only 
			// used per component. The fraction is used in all, either with the fraction as the 
			// numerator (A) or part of the denominator (B) 
			// E.G. 
			// f(A) = A / (A+B+C)
			// Quotient rule: f = g / h => f' = (g'h - gh') / h^2
			// f'(A) = ((A+B+C) - A) / (A+B+C)^2
			//       = (B+C) / (A+B+C)^2
			//       = (sum(f) - f) / sum(f)^2
			// f'(B) = -A / (A+B+C)^2
			//       = -f / sum(f)^2

			// Differentiate with respect to D is easier:
			// y' = -a'
			// a' = f * exp(b) * b'
			// b' = -1 * -x / 4D^2 = x / 4D^2
			// y' = f * -exp(b) * x / 4D^2
			//    = f * -a * -b / D
			//    = f * a * b / D
			//    = f * exp(b) * b / D

			final double[] fourD = new double[n];
			final double[] f = new double[n];
			double total = 0;
			for (int i = 0; i < n; i++)
			{
				f[i] = variables[i * 2];
				fourD[i] = 4 * variables[i * 2 + 1];
				total += f[i];
			}

			final double[] fraction = new double[n];
			final double[] total_f = new double[n];
			final double[] f_total = new double[n];
			for (int i = 0; i < n; i++)
			{
				fraction[i] = f[i] / total;
				// Because we use y = 1 - sum(a) all coefficients are inverted
				total_f[i] = -1 * (total - f[i]) / (total * total);
				f_total[i] = -1 * -f[i] / (total * total);
			}

			double[][] jacobian = new double[x.length][variables.length];

			double[] b = new double[n];
			for (int i = 0; i < x.length; ++i)
			{
				for (int j = 0; j < n; j++)
					b[j] = -x[i] / fourD[j];

				for (int j = 0; j < n; j++)
				{
					// Gradient for the diffusion coefficient
					jacobian[i][j * 2 + 1] = fraction[j] * FastMath.exp(b[j]) * b[j] / variables[j * 2 + 1];

					// Gradient for the fraction f
					jacobian[i][j * 2] = total_f[j] * FastMath.exp(b[j]);
					for (int k = 0; k < n; k++)
					{
						if (j == k)
							continue;
						jacobian[i][j * 2] += f_total[k] * FastMath.exp(b[k]);
					}
				}
			}

			//// Check numerically ...
			//double[][] jacobian2 = super.jacobian(variables);
			//for (int i = 0; i < jacobian.length; i++)
			//{
			//	StringBuilder sb = new StringBuilder();
			//	for (int j = 0; j < jacobian[i].length; j++)
			//	{
			//		sb.append(" d").append(j).append(" = ").append(jacobian[i][j]).append(" : ")
			//				.append(jacobian2[i][j]).append(" = ")
			//				.append(DoubleEquality.relativeError(jacobian[i][j], jacobian2[i][j]));
			//	}
			//	System.out.println(sb.toString());
			//}

			return jacobian;
		}
	}

	public class MixedJumpDistanceFunctionMultivariate extends MixedJumpDistanceFunction implements
			MultivariateFunction
	{
		public MixedJumpDistanceFunctionMultivariate(double[] x, double[] y, double estimatedD, int n)
		{
			super(x, y, estimatedD, n);
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see org.apache.commons.math3.analysis.MultivariateFunction#value(double[])
		 */
		public double value(double[] parameters)
		{
			double[] obs = getValue(parameters);

			// Optimise the sum of squares
			double ss = 0;
			for (int i = x.length; i-- > 0;)
			{
				double dx = y[i] - obs[i];
				ss += dx * dx;
			}
			return ss;
		}
	}
}
