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

import gdsc.smlm.fitting.JumpDistanceAnalysis;
import gdsc.smlm.fitting.JumpDistanceAnalysis.CurveLogger;
import gdsc.smlm.ij.IJTrackProgress;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.settings.ClusteringSettings;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.IJLogger;
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.Trace;
import gdsc.smlm.results.TraceManager;
import gdsc.smlm.utils.Maths;
import gdsc.smlm.utils.Random;
import gdsc.smlm.utils.Statistics;
import gdsc.smlm.utils.StoredDataStatistics;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.gui.Plot2;
import ij.gui.PlotWindow;
import ij.plugin.PlugIn;
import ij.plugin.WindowOrganiser;
import ij.text.TextWindow;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.analysis.MultivariateMatrixFunction;
import org.apache.commons.math3.analysis.MultivariateVectorFunction;
import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyIterationsException;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.PointVectorValuePair;
import org.apache.commons.math3.optim.nonlinear.vector.ModelFunction;
import org.apache.commons.math3.optim.nonlinear.vector.ModelFunctionJacobian;
import org.apache.commons.math3.optim.nonlinear.vector.Target;
import org.apache.commons.math3.optim.nonlinear.vector.Weight;
import org.apache.commons.math3.optim.nonlinear.vector.jacobian.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.util.FastMath;

/**
 * Run a tracing algorithm on the peak results to trace molecules across the frames.
 */
public class TraceDiffusion implements PlugIn, CurveLogger
{
	private static final String TITLE = "Trace Diffusion";
	private static String inputOption = "";
	private static String header = null;
	private static TextWindow summaryTable = null;

	private static final String[] NAMES = new String[] { "Total Signal", "Signal/Frame", "t-On (s)" };
	private static final boolean[] ROUNDED = new boolean[] { false, false, true };
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
	private static boolean displayTraceLength = false;

	private static boolean saveTraceDistances = false;
	private static boolean saveRawData = false;
	private static String rawDataDirectory = "";
	private boolean directoryChosen = false;
	private static String distancesFilename = "";
	private static double minFraction = 0.1;
	private static double minDifference = 2;
	private static int minN = 1;
	private static int maxN = 5;
	private static boolean debugFitting = false;
	private static boolean multipleInputs = false;
	private static String tracesFilename = "";
	private static String title = "";

	private GlobalSettings globalSettings;
	private ClusteringSettings settings;
	private MemoryPeakResults results;
	private boolean extraOptions;
	private int myMinN = 1;

	// The number of additional datasets
	private int additionalDatasets = 0;

	// Store exposure time in seconds
	private double exposureTime = 0;
	private double precision;

	// Used to tile new plot windows
	private int[] idList = new int[20];
	private int idCount = 0;

	private String jdTitle = TITLE + " Jump Distance";
	private Plot2 jdPlot;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		extraOptions = Utils.isExtraOptions();
		if (MemoryPeakResults.countMemorySize() == 0)
		{
			IJ.error(TITLE, "No localisations in memory");
			return;
		}

		// -=-=-
		// Get the traces:
		// - Trace a single dataset (and store in memory)
		// - Combine trace results held in memory
		// -=-=-

		if (!showTraceDialog())
			return;

		Utils.log(TITLE + "...");
		Trace[] traces = getTraces();

		// -=-=-
		// Analyse the traces
		// -=-=-

		if (!showDialog())
			return;

		int count = traces.length;
		double D = 0;
		int n = 0;
		double[][] jdParams = null;
		if (count > 0)
		{
			//--- MSD Analysis ---

			// Conversion constants
			final double px2ToUm2 = results.getCalibration().nmPerPixel * results.getCalibration().nmPerPixel / 1e6;
			final double px2ToUm2PerSecond = px2ToUm2 / exposureTime;

			// Get the maximum trace length
			int length = settings.minimumTraceLength;
			if (!settings.truncate)
			{
				for (Trace trace : traces)
				{
					if (length < trace.size())
						length = trace.size();
				}
			}

			// Extract the mean-squared distance statistics
			Statistics[] stats = new Statistics[length];
			// Disable sub-sampling
			final boolean subSample = false; // (settings.internalDistances & settings.subSampledDistances);
			for (int i = 0; i < stats.length; i++)
				stats[i] = (subSample) ? new StoredDataStatistics() : new Statistics();

			ArrayList<double[]> distances = (saveTraceDistances || displayTraceLength) ? new ArrayList<double[]>(
					traces.length) : null;

			// Store all the jump distances at the specified interval
			StoredDataStatistics jumpDistances = new StoredDataStatistics();
			final int jumpDistanceInterval = settings.jumpDistance;
			final double jdPx2ToUm2PerSecond = px2ToUm2 / (jumpDistanceInterval * exposureTime);

			// Pre-calculate conversion factors
			double[] convert = null;
			if (settings.msdCorrection)
			{
				convert = new double[length];
				for (int t = 1; t < length; t++)
					convert[t] = JumpDistanceAnalysis.getConversionfactor(t);
			}

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
					if (distances != null)
					{
						double[] msd = new double[traceLength - 1];
						for (int j = 1; j < traceLength; j++)
						{
							final int t = j;
							double d = distance2(x, y, results.get(j));
							if (settings.msdCorrection)
								d *= convert[t];
							msd[j - 1] = px2ToUm2 * d;
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
							final int t = j;
							double d = distance2(x, y, results.get(j));
							if (settings.msdCorrection)
								d *= convert[t];
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
						for (int j = i + 1; j < traceLength; j++)
						{
							final int t = j - i;
							double d = distance2(x, y, results.get(j));
							if (settings.msdCorrection)
								d *= convert[t];
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

			if (displayTraceLength)
			{
				StoredDataStatistics lengths = calculateTraceLengths(distances);
				showHistogram(lengths, "Trace length (um)");
			}

			// Plot the per-trace histogram of MSD and D*
			if (settings.showHistograms)
			{
				if (displayMSDHistogram)
				{
					showHistogram(msdPerMoleculeAllVsAll, "MSD/Molecule (all-vs-all)");
					showHistogram(msdPerMoleculeAdjacent, "MSD/Molecule (adjacent)");
				}
				if (displayDHistogram)
				{
					showHistogram(dStarPerMoleculeAllVsAll, "D*/Molecule (all-vs-all)");
					showHistogram(dStarPerMoleculeAdjacent, "D*/Molecule (adjacent)");
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
			Plot2 plot = plotMSD(x, y, sd, title);

			// Fit the MSD using a linear fit that must pass through 0,0
			D = fitMSD(x, y, title, plot);

			// Jump Distance analysis
			if (saveRawData)
				saveStatistics(jumpDistances, "Jump Distance", "Distance (um^2/second)", false);

			// Calculate the cumulative jump-distance histogram
			double[][] jdHistogram = JumpDistanceAnalysis.cumulativeHistogram(jumpDistances.getValues());

			// Always show the jump distance histogram
			jdTitle = TITLE + " Jump Distance";
			jdPlot = new Plot2(jdTitle, "Distance (um^2/second)", "Cumulative Probability", jdHistogram[0],
					jdHistogram[1]);
			display(jdTitle, jdPlot);

			// Fit Jump Distance cumulative probability
			n = jumpDistances.getN();
			jdParams = fitJumpDistance(jumpDistances, jdHistogram);
		}

		summarise(traces, D, n, jdParams);
	}

	public StoredDataStatistics calculateTraceLengths(ArrayList<double[]> distances)
	{
		StoredDataStatistics lengths = new StoredDataStatistics();
		for (double[] trace : distances)
		{
			double sum = 0;
			for (double d : trace)
			{
				sum += Math.sqrt(d);
			}
			lengths.add(sum);
		}
		return lengths;
	}

	private void display(String title, Plot2 plot)
	{
		PlotWindow pw = Utils.display(title, plot);
		if (Utils.isNewWindow())
			idList[idCount++] = pw.getImagePlus().getID();
	}

	private void showHistogram(StoredDataStatistics stats, String title)
	{
		showHistogram(stats, title, false, false);
	}

	private void showHistogram(StoredDataStatistics stats, String title, boolean alwaysRemoveOutliers, boolean rounded)
	{
		if (saveRawData)
			saveStatistics(stats, title, title, rounded);

		int id = Utils.showHistogram(TITLE, stats, title, 0, (settings.removeOutliers || alwaysRemoveOutliers) ? 1 : 0,
				settings.histogramBins);
		if (Utils.isNewWindow())
			idList[idCount++] = id;
	}

	private void tileNewWindows()
	{
		if (idCount > 0)
		{
			idList = Arrays.copyOf(idList, idCount);
			new WindowOrganiser().tileWindows(idList);
		}
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
		distancesFilename = Utils.getFilename("Trace_Distances_File", distancesFilename);
		if (distancesFilename != null)
		{
			distancesFilename = Utils.replaceExtension(distancesFilename, "xls");

			BufferedWriter out = null;
			try
			{
				FileOutputStream fos = new FileOutputStream(distancesFilename);
				out = new BufferedWriter(new OutputStreamWriter(fos, "UTF-8"));
				double[] msd = msdPerMolecule.getValues();
				double[] msd2 = msdPerMoleculeAdjacent.getValues();
				double[] dStar = dStarPerMolecule.getValues();
				double[] dStar2 = dStarPerMoleculeAdjacent.getValues();
				out.write(String.format("#%d traces : Precision = %s nm : Exposure time = %s s", nTraces,
						Utils.rounded(precision, 4), Utils.rounded(exposureTime, 4)));
				out.newLine();
				out.write(String
						.format("#TraceId\tMSD all-vs-all (um^2/s)\tMSD adjacent (um^2/s)\tD* all-vs-all(um^2/s)\tD* adjacent(um^2/s)\tDistances (um^2) per %ss ... ",
								Utils.rounded(exposureTime, 4)));
				out.newLine();
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
					out.newLine();
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

	private void saveMSD(double[] x, double[] y, double[] se)
	{
		if (!directoryChosen)
			rawDataDirectory = Utils.getDirectory("Data_directory", rawDataDirectory);
		directoryChosen = true;
		if (rawDataDirectory == null)
			return;
		String filename = rawDataDirectory + "MSD.txt";

		BufferedWriter out = null;
		try
		{
			FileOutputStream fos = new FileOutputStream(filename);
			out = new BufferedWriter(new OutputStreamWriter(fos, "UTF-8"));
			out.write("Time (s)\tDistance (um^2)\tS.E.");
			out.newLine();
			for (int i = 0; i < x.length; i++)
			{
				out.write(Utils.rounded(x[i]));
				out.write('\t');
				out.write(Double.toString(y[i]));
				out.write('\t');
				out.write(Double.toString(se[i]));
				out.newLine();
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

	private void saveStatistics(StoredDataStatistics stats, String title, String label, boolean rounded)
	{
		if (!directoryChosen)
			rawDataDirectory = Utils.getDirectory("Data_directory", rawDataDirectory);
		directoryChosen = true;
		if (rawDataDirectory == null)
			return;
		String filename = rawDataDirectory + title.replace("/", " per ").replace("*", "star") + ".txt";

		BufferedWriter out = null;
		try
		{
			FileOutputStream fos = new FileOutputStream(filename);
			out = new BufferedWriter(new OutputStreamWriter(fos, "UTF-8"));
			out.write(label);
			out.newLine();
			double[] data = stats.getValues();
			Arrays.sort(data);
			if (rounded)
			{
				for (double d : data)
				{
					out.write(Utils.rounded(d, 4));
					out.newLine();
				}
			}
			else
			{
				for (double d : data)
				{
					out.write(Double.toString(d));
					out.newLine();
				}
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

	private double distance2(final float x, final float y, PeakResult r2)
	{
		final double dx = x - r2.getXPosition();
		final double dy = y - r2.getYPosition();
		return dx * dx + dy * dy;
	}

	/**
	 * Filter traces that are not the minimum length
	 * 
	 * @param name
	 * @param traces
	 * @param minimumTraceLength
	 * @param ignoreEnds
	 * @return The new traces
	 */
	private Trace[] filterTraces(String name, Trace[] traces, int minimumTraceLength, boolean ignoreEnds)
	{
		final int minLength = (ignoreEnds) ? minimumTraceLength + 2 : minimumTraceLength;
		int count = 0;
		for (int i = 0; i < traces.length; i++)
		{
			if (traces[i].size() >= minLength)
			{
				if (ignoreEnds)
					traces[i].removeEnds();
				traces[count++] = traces[i];
			}
		}

		Utils.log("Filtered results '%s' : %s filtered to %d using minimum length %d (Ignore ends = %b)", name,
				Utils.pleural(traces.length, "trace"), count, minimumTraceLength, ignoreEnds);
		return Arrays.copyOf(traces, count);
	}

	private String createSettingsComment()
	{
		return String.format("Molecule tracing : distance-threshold = %f nm", settings.distanceThreshold);
	}

	private void summarise(Trace[] traces, double D, int n, double[][] jdParams)
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
		StringBuilder sb = new StringBuilder(title);
		sb.append('\t').append(results.getName());
		if (additionalDatasets > 0)
		{
			sb.append(" + ").append(additionalDatasets).append(" others");
		}
		sb.append("\t");
		sb.append(Utils.rounded(exposureTime * 1000, 3)).append("\t");
		sb.append(Utils.rounded(settings.distanceThreshold, 3)).append("\t");
		sb.append(Utils.rounded(settings.distanceExclusion, 3)).append("\t");
		sb.append(settings.minimumTraceLength).append("\t");
		sb.append(settings.ignoreEnds).append("\t");
		sb.append(settings.truncate).append("\t");
		sb.append(settings.internalDistances).append("\t");
		sb.append(settings.fitLength).append("\t");
		sb.append(settings.msdCorrection).append("\t");
		sb.append(settings.mle).append("\t");
		sb.append(traces.length).append("\t");
		sb.append(Utils.rounded(D, 4)).append("\t");
		sb.append(Utils.rounded(settings.jumpDistance * exposureTime)).append("\t");
		sb.append(n).append("\t");
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

			for (int i = 0; i < NAMES.length; i++)
			{
				if (displayHistograms[i])
				{
					showHistogram((StoredDataStatistics) stats[i], NAMES[i], alwaysRemoveOutliers[i], ROUNDED[i]);
				}
			}
		}

		tileNewWindows();

		IJ.showStatus("Finished " + TITLE);
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
				"Title\tDataset\tExposure time (ms)\tD-threshold (nm)\tEx-threshold (nm)\tMin.Length\tIgnoreEnds\tTruncate\tInternal\tFit Length\tCorrection\tMLE\tTraces\tD (um^2/s)\tJump Distance (s)\tN\tJump D (um^2/s)\tFractions");
		for (int i = 0; i < NAMES.length; i++)
		{
			sb.append("\t").append(NAMES[i]);
		}
		return sb.toString();
	}

	private boolean showTraceDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);

		globalSettings = SettingsManager.loadSettings();
		settings = globalSettings.getClusteringSettings();

		gd.addNumericField("Distance_Threshold (nm)", settings.distanceThreshold, 0);
		gd.addNumericField("Distance_Exclusion (nm)", settings.distanceExclusion, 0);
		gd.addSlider("Min_trace_length", 2, 20, settings.minimumTraceLength);
		gd.addCheckbox("Ignore_ends", settings.ignoreEnds);
		gd.addCheckbox("Save_traces", settings.saveTraces);
		gd.addCheckbox("Multiple_inputs", multipleInputs);

		gd.showDialog();

		if (gd.wasCanceled() || !readTraceDialog(gd))
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

	private boolean readTraceDialog(GenericDialog gd)
	{
		inputOption = ResultsManager.getInputSource(gd);
		settings.distanceThreshold = gd.getNextNumber();
		settings.distanceExclusion = Math.abs(gd.getNextNumber());
		settings.minimumTraceLength = (int) Math.abs(gd.getNextNumber());
		settings.ignoreEnds = gd.getNextBoolean();
		settings.saveTraces = gd.getNextBoolean();
		multipleInputs = gd.getNextBoolean();

		if (gd.invalidNumber())
			return false;

		// Check arguments
		try
		{
			Parameters.isAboveZero("Distance threshold", settings.distanceThreshold);
			Parameters.isAbove("Min trace length", settings.minimumTraceLength, 1);
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		return true;
	}

	private Trace[] getTraces()
	{
		ArrayList<MemoryPeakResults> allResults = new ArrayList<MemoryPeakResults>();
		allResults.add(results);

		final double nmPerPixel = results.getCalibration().nmPerPixel;
		if (multipleInputs)
		{
			// Get additional results sets with the same calibration
			MemoryPeakResults r = nextInput(nmPerPixel, allResults);
			while (r != null)
			{
				allResults.add(r);
				r = nextInput(nmPerPixel, allResults);
			}
		}

		ArrayList<Trace> allTraces = new ArrayList<Trace>();
		additionalDatasets = -1;
		for (MemoryPeakResults r : allResults)
		{
			additionalDatasets++;

			TraceManager manager = new TraceManager(r);

			// Run the tracing
			manager.setTracker(new IJTrackProgress());
			manager.setDistanceExclusion(settings.distanceExclusion / nmPerPixel);
			manager.traceMolecules(settings.distanceThreshold / nmPerPixel, 1);
			Trace[] traces = manager.getTraces();

			traces = filterTraces(r.getName(), traces, settings.minimumTraceLength, settings.ignoreEnds);
			allTraces.addAll(Arrays.asList(traces));

			//--- Save results ---
			if (traces.length > 0)
			{
				// Save the traces to memory
				TraceMolecules.saveResults(r, traces, "Tracks");

				if (settings.saveTraces)
				{
					// Sort traces by time to assist the results source in extracting frames sequentially.
					// Do this before saving to assist in debugging using the saved traces file.
					TraceMolecules.sortByTime(traces);
					String newFilename = TraceMolecules.saveTraces(r, traces, createSettingsComment(), tracesFilename,
							additionalDatasets);
					// Only keep the main filename in memory
					if (additionalDatasets == 0)
						tracesFilename = newFilename;
				}
			}
		}

		if (additionalDatasets > 0)
			Utils.log("Multiple inputs provide %d traces", allTraces.size());

		return allTraces.toArray(new Trace[allTraces.size()]);
	}

	private MemoryPeakResults nextInput(double nmPerPixel, ArrayList<MemoryPeakResults> allResults)
	{
		ArrayList<String> source = new ArrayList<String>(3);
		boolean fileInput = false;

		for (MemoryPeakResults r : MemoryPeakResults.getAllResults())
		{
			if (notValid(r, allResults, nmPerPixel))
				continue;
			ResultsManager.addInputSource(source, r, InputSource.MEMORY_SINGLE_FRAME);
		}
		if (source.isEmpty())
			return null;

		String inputOption;

		// if a macro then use the recorder to get the option

		GenericDialog gd = new GenericDialog(TITLE);
		gd.addMessage("Select additional inputs ...\n \nPress cancel to continue with the analysis.");

		// If in macro mode then we must just use the String input field to allow the macro
		// IJ to return the field values from the macro arguments. Using a Choice input
		// will always return a field value.

		String fieldName = "Input" + allResults.size();
		if (IJ.isMacro())
			// Use blank default value so bad macro parameters return nothing
			gd.addStringField(fieldName, "");
		else
			ResultsManager.addInputSourceToDialog(gd, fieldName, "", source, fileInput);

		gd.showDialog();

		if (gd.wasCanceled())
			return null;

		if (IJ.isMacro())
			inputOption = gd.getNextString();
		else
			inputOption = ResultsManager.getInputSource(gd);
		MemoryPeakResults results = ResultsManager.loadInputResults(inputOption, true);
		if (results == null || results.size() == 0)
		{
			return null;
		}

		// Check the results have the same calibrated exposure time and pixel size
		if (results.getCalibration() == null || results.getCalibration().exposureTime / 1000.0 != exposureTime ||
				results.getNmPerPixel() != nmPerPixel)
		{
			return null;
		}

		return results;
	}

	/**
	 * Check if the results are valid for inclusion as additional datasets
	 * 
	 * @param r
	 *            The results
	 * @param allResults
	 *            All the current results
	 * @param nmPerPixel
	 *            The calibrated pixel size of the primary results
	 * @return True if the results are not valid to be included
	 */
	private boolean notValid(MemoryPeakResults r, ArrayList<MemoryPeakResults> allResults, double nmPerPixel)
	{
		// Check the calibration is the same
		if (r.getNmPerPixel() != nmPerPixel)
			return true;
		if (r.getCalibration().exposureTime / 1000.0 != exposureTime)
			return true;
		// Check the results have not already been chosen
		for (MemoryPeakResults r2 : allResults)
		{
			if (r2.getName().equals(r.getName()))
				return true;
		}
		return false;
	}

	private boolean showDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		globalSettings = SettingsManager.loadSettings();
		settings = globalSettings.getClusteringSettings();

		gd.addCheckbox("Truncate_traces", settings.truncate);
		gd.addCheckbox("Internal_distances", settings.internalDistances);
		//gd.addCheckbox("Sub-sample_distances", settings.subSampledDistances);
		gd.addSlider("Fit_length", 2, 20, settings.fitLength);
		gd.addCheckbox("MSD_correction", settings.msdCorrection);
		gd.addCheckbox("Maximum_likelihood", settings.mle);
		gd.addSlider("Fit_restarts", 0, 10, settings.fitRestarts);
		gd.addSlider("Jump_distance", 1, 20, settings.jumpDistance);
		gd.addSlider("Minimum_difference", 0, 10, minDifference);
		gd.addSlider("Minimum_fraction", 0, 1, minFraction);
		if (extraOptions)
			gd.addSlider("Minimum_N", 1, 10, minN);
		gd.addSlider("Maximum_N", 2, 10, maxN);
		gd.addCheckbox("Debug_fitting", debugFitting);
		gd.addCheckbox("Save_trace_distances", saveTraceDistances);
		gd.addCheckbox("Save_raw_data", saveRawData);
		gd.addCheckbox("Show_histograms", settings.showHistograms);
		gd.addStringField("Title", title);

		gd.showDialog();

		if (gd.wasCanceled() || !readDialog(gd))
			return false;

		// Update the settings
		SettingsManager.saveSettings(globalSettings);

		return true;
	}

	private boolean readDialog(GenericDialog gd)
	{
		settings.truncate = gd.getNextBoolean();
		settings.internalDistances = gd.getNextBoolean();
		//settings.subSampledDistances = gd.getNextBoolean();
		settings.fitLength = (int) Math.abs(gd.getNextNumber());
		settings.msdCorrection = gd.getNextBoolean();
		settings.mle = gd.getNextBoolean();
		settings.fitRestarts = (int) Math.abs(gd.getNextNumber());
		settings.jumpDistance = (int) Math.abs(gd.getNextNumber());
		minDifference = Math.abs(gd.getNextNumber());
		minFraction = Math.abs(gd.getNextNumber());
		if (extraOptions)
			myMinN = minN = (int) Math.abs(gd.getNextNumber());
		maxN = (int) Math.abs(gd.getNextNumber());
		debugFitting = gd.getNextBoolean();
		saveTraceDistances = gd.getNextBoolean();
		saveRawData = gd.getNextBoolean();
		settings.showHistograms = gd.getNextBoolean();
		title = gd.getNextString();

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
			gd.addCheckbox("Trace_length", displayTraceLength);
			gd.showDialog();
			if (gd.wasCanceled())
				return false;
			settings.removeOutliers = gd.getNextBoolean();
			settings.histogramBins = (int) Math.abs(gd.getNextNumber());
			for (int i = 0; i < displayHistograms.length; i++)
				displayHistograms[i] = gd.getNextBoolean();
			displayMSDHistogram = gd.getNextBoolean();
			displayDHistogram = gd.getNextBoolean();
			displayTraceLength = gd.getNextBoolean();
		}

		// Check arguments
		try
		{
			Parameters.isAboveZero("Histogram bins", settings.histogramBins);
			Parameters.isAbove("Fit length", settings.fitLength, 1);
			Parameters.isAboveZero("Jump distance", settings.jumpDistance);
			Parameters.isEqualOrAbove("Maximum N", maxN, myMinN);
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		return true;
	}

	private Plot2 plotMSD(double[] x, double[] y, double[] sd, String title)
	{
		if (saveRawData)
			saveMSD(x, y, sd);

		Plot2 plot = new Plot2(title, "Time (s)", "Distance (um^2)", x, y);
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
		display(title, plot);
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
	private double fitMSD(double[] x, double[] y, String title, Plot2 plot)
	{
		// TODO - The Weimann paper (Plos One e64287) fits:
		// MSD(n dt) = 4D n dt + 4s^2
		// n = number of jumps
		// dt = time difference between frames
		// s = localisation precision
		// Thus we should fit an intercept as well.

		double D = 0;
		LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer();
		PointVectorValuePair lvmSolution;
		double ic = 0;
		double gradient = 0, intercept = 0;
		try
		{
			final LinearFunction function = new LinearFunction(x, y, settings.fitLength);
			double[] parameters = new double[] { function.guess() };
			lvmSolution = optimizer.optimize(new MaxIter(3000), new MaxEval(Integer.MAX_VALUE),
					new ModelFunctionJacobian(new MultivariateMatrixFunction()
					{
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

			ic = Maths.getInformationCriterion(ss, obs.length, 1);

			gradient = lvmSolution.getPoint()[0];
			D = gradient / 4;

			Utils.log("Linear fit (%d points) : Gradient = %s, D = %s um^2/s, SS = %f, IC = %f (%d evaluations)",
					obs.length, Utils.rounded(gradient, 4), Utils.rounded(D, 4), ss, ic, optimizer.getEvaluations());
		}
		catch (TooManyIterationsException e)
		{
			Utils.log("Failed to fit : Too many iterations (%d)", optimizer.getIterations());
		}
		catch (ConvergenceException e)
		{
			Utils.log("Failed to fit : %s", e.getMessage());
		}

		try
		{
			final LinearFunctionWithIntercept function = new LinearFunctionWithIntercept(x, y, settings.fitLength);
			lvmSolution = optimizer.optimize(new MaxIter(3000), new MaxEval(Integer.MAX_VALUE),
					new ModelFunctionJacobian(new MultivariateMatrixFunction()
					{
						public double[][] value(double[] point) throws IllegalArgumentException
						{
							return function.jacobian(point);
						}
					}), new ModelFunction(function), new Target(function.getY()), new Weight(function.getWeights()),
					new InitialGuess(function.guess()));

			double ss = 0;
			double[] obs = function.getY();
			double[] exp = lvmSolution.getValue();
			for (int i = 0; i < obs.length; i++)
				ss += (obs[i] - exp[i]) * (obs[i] - exp[i]);

			double ic2 = Maths.getInformationCriterion(ss, obs.length, 1);
			double gradient2 = lvmSolution.getPoint()[0];
			final double s = lvmSolution.getPoint()[1];
			double intercept2 = 4 * s * s;

			if (ic2 < ic || debugFitting)
			{
				// Convert fitted precision in um to nm
				Utils.log(
						"Linear fit with intercept (%d points) : Gradient = %s, Intercept = %s, D = %s um^2/s, precision = %s nm, SS = %f, IC = %f (%d evaluations)",
						obs.length, Utils.rounded(gradient2, 4), Utils.rounded(intercept2, 4),
						Utils.rounded(gradient / 4, 4), Utils.rounded(s * 1000, 4), ss, ic2, optimizer.getEvaluations());
			}

			if (lvmSolution == null || ic2 < ic)
			{
				gradient = gradient2;
				intercept = intercept2;
				D = gradient / 4;
			}

			// Add the fit to the plot
			plot.setColor(Color.magenta);
			plot.drawLine(0, intercept, x[x.length - 1], gradient * x[x.length - 1] + intercept);
			display(title, plot);
		}
		catch (TooManyIterationsException e)
		{
			Utils.log("Failed to fit with intercept : Too many iterations (%d)", optimizer.getIterations());
		}
		catch (ConvergenceException e)
		{
			Utils.log("Failed to fit with intercept : %s", e.getMessage());
		}

		return D;
	}

	public class LinearFunction implements MultivariateVectorFunction
	{
		double[] x, y;
		double[][] jacobian;

		public LinearFunction(double[] x, double[] y, int length)
		{
			int to = FastMath.min(x.length, 1 + length);
			this.x = Arrays.copyOfRange(x, 1, to);
			this.y = Arrays.copyOfRange(y, 1, to);
			jacobian = calculateJacobian();
		}

		// Adapted from http://commons.apache.org/proper/commons-math/userguide/optimization.html

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
			return jacobian;
		}

		double[][] calculateJacobian()
		{
			// Compute the gradients using calculus differentiation:
			// y = ax + c
			// dy_da = x
			double[][] jacobian = new double[x.length][1];

			for (int i = 0; i < jacobian.length; ++i)
			{
				jacobian[i][0] = x[i];
			}

			return jacobian;
		}
	}

	public class LinearFunctionWithIntercept implements MultivariateVectorFunction
	{
		double[] x, y;

		public LinearFunctionWithIntercept(double[] x, double[] y, int length)
		{
			int to = FastMath.min(x.length, 1 + length);
			this.x = Arrays.copyOfRange(x, 1, to);
			this.y = Arrays.copyOfRange(y, 1, to);
		}

		// Adapted from http://commons.apache.org/proper/commons-math/userguide/optimization.html

		/**
		 * @return An estimate for the linear gradient and intercept
		 */
		public double[] guess()
		{
			if (y.length == 1)
				return new double[] { y[0] / x[0], 0 };

			double a = (y[y.length - 1] - y[0]) / (x[x.length - 1] - x[0]);
			// y = ax + 4c^2
			// y = ax + intercept
			// intercept = y - ax
			//           = 4c^2
			double intercept = y[y.length - 1] - a * x[x.length - 1];
			double c = (intercept < 0) ? 0 : Math.sqrt(intercept / 4);

			return new double[] { a, c };
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
			// y = ax + 4c^2
			final double[] values = new double[x.length];
			final double a = variables[0];
			final double intercept = 4 * variables[1] * variables[1];
			for (int i = 0; i < values.length; i++)
			{
				values[i] = a * x[i] + intercept;
			}
			return values;
		}

		double[][] jacobian(double[] variables)
		{
			// Compute the gradients using calculus differentiation:
			// y = ax + 4c^2
			// dy_da = x
			// dy_dc = 8c
			double[][] jacobian = new double[x.length][2];
			final double dy_dc = 8 * variables[1];

			for (int i = 0; i < jacobian.length; ++i)
			{
				jacobian[i][0] = x[i];
				jacobian[i][1] = dy_dc;
			}

			return jacobian;
		}
	}

	/**
	 * Fit the jump distance histogram.
	 * <p>
	 * Update the plot by adding the fit line(s).
	 * 
	 * @param jumpDistances
	 * @param jdHistogram
	 * @return The fitted coefficients and fractions
	 */
	private double[][] fitJumpDistance(StoredDataStatistics jumpDistances, double[][] jdHistogram)
	{
		final double msd = jumpDistances.getMean();
		final double meanDistance = Math.sqrt(msd) * 1e3;
		final double beta = meanDistance / precision;
		Utils.log(
				"Jump Distance analysis : N = %d, Time = %d frames (%s seconds). MSD = %s um^2/second, Mean Distance = %s nm, Precision = %s nm, Beta = %s",
				jumpDistances.getN(), settings.jumpDistance, Utils.rounded(settings.jumpDistance * exposureTime, 4),
				Utils.rounded(msd, 4), Utils.rounded(meanDistance, 4), Utils.rounded(precision, 4),
				Utils.rounded(beta, 4));

		IJLogger logger = new IJLogger(debugFitting, debugFitting);
		JumpDistanceAnalysis jd = new JumpDistanceAnalysis(logger);
		jd.setFitRestarts(settings.fitRestarts);
		jd.setMinFraction(minFraction);
		jd.setMinDifference(minDifference);
		jd.setMinN(myMinN);
		jd.setMaxN(maxN);
		// Update the plot with the fit
		jd.setCurveLogger(this);
		double[][] fit;
		if (settings.mle)
			fit = jd.fitJumpDistancesMLE(jumpDistances.getValues(), jdHistogram);
		else
			fit = jd.fitJumpDistanceHistogram(jumpDistances.getMean(), jdHistogram);
		return fit;
	}

	public int getNumberOfCurvePoints()
	{
		return 300;
	}

	public void saveSinglePopulationCurve(double[][] curve)
	{
		addToJumpDistancePlot(curve[0], curve[1], Color.magenta);
	}

	public void saveMixedPopulationCurve(double[][] curve)
	{
		addToJumpDistancePlot(curve[0], curve[1], Color.yellow);
	}

	private void addToJumpDistancePlot(double[] x, double[] y, Color color)
	{
		jdPlot.setColor(color);
		jdPlot.addPoints(x, y, Plot2.LINE);
		display(jdTitle, jdPlot);
	}
}
