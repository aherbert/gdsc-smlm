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
import gdsc.smlm.results.Calibration;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.Trace;
import gdsc.smlm.results.TraceManager;
import gdsc.smlm.utils.Maths;
import gdsc.smlm.utils.Statistics;
import gdsc.smlm.utils.StoredDataStatistics;
import ij.IJ;
import ij.Macro;
import ij.WindowManager;
import ij.gui.GUI;
import ij.gui.GenericDialog;
import ij.gui.Plot2;
import ij.gui.PlotWindow;
import ij.plugin.PlugIn;
import ij.plugin.WindowOrganiser;
import ij.plugin.frame.Recorder;
import ij.text.TextWindow;

import java.awt.BorderLayout;
import java.awt.Button;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dialog;
import java.awt.FlowLayout;
import java.awt.Frame;
import java.awt.List;
import java.awt.Panel;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Set;

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
	private boolean extraOptions, multiMode;
	private int myMinN = 1;

	// The number of additional datasets
	private int additionalDatasets = 0;

	// Store exposure time in seconds
	private double exposureTime = 0;
	private double precision, beta;
	private double ic = Double.NaN;

	// Used to tile new plot windows
	private int[] idList = new int[20];
	private int idCount = 0;

	private String jdTitle = TITLE + " Jump Distance";
	private Plot2 jdPlot;

	// Used for the macro extensions
	private static double[][] jumpDistanceParameters = null;

	// Used for the multiMode option 
	private static ArrayList<String> selected = new ArrayList<String>();

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		jumpDistanceParameters = null;

		extraOptions = Utils.isExtraOptions();
		if (MemoryPeakResults.countMemorySize() == 0)
		{
			IJ.error(TITLE, "No localisations in memory");
			return;
		}

		ArrayList<MemoryPeakResults> allResults = new ArrayList<MemoryPeakResults>();

		// Option to pick multiple input datasets together using a list box.
		// (Do not support running in macros. The macro options are recorded as if running 
		// using the multiple_inputs option in the trace dialog.)
		if ("multi".equals(arg) && Macro.getOptions() == null)
		{
			if (!showMultiDialog(allResults))
				return;
		}

		// This shows the dialog for selecting trace options
		if (!showTraceDialog(allResults))
			return;

		if (allResults.isEmpty()) // Sense check
			return;

		Utils.log(TITLE + "...");

		// This optionally collects additional datasets then gets the traces:
		// - Trace each single dataset (and store in memory)
		// - Combine trace results held in memory
		Trace[] traces = getTraces(allResults);

		// -=-=-
		// Analyse the traces
		// -=-=-

		// Only show the second dialog if we have traces. 
		// This still allows a zero entry in the results table.
		if (traces.length > 0)
			if (!showDialog())
				return;

		int count = traces.length;
		double[] fitMSDResult = null;
		int n = 0;
		double[][] jdParams = null;
		if (count > 0)
		{
			calculatePrecision(traces, allResults.size() > 1);

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

			// Get the localisation error (4s^2) in um^2
			final double error = (settings.precisionCorrection) ? 4 * precision * precision / 1e6 : 0;
			// Pre-calculate MSD correction factors. This accounts for the fact that the distance moved 
			// in the start/end frames is reduced due to the averaging of the particle location over the 
			// entire frame into a single point. The true MSD may be restored by applying a factor.
			// Note: These are used for the calculation of the diffusion coefficients per molecule and 
			// the MSD passed to the Jump Distance analysis. However the error is not included in the 
			// jump distance analysis so will be subtracted from the fitted D coefficients later.
			final double[] factors;
			if (settings.msdCorrection)
			{
				factors = new double[length];
				for (int t = 1; t < length; t++)
					factors[t] = JumpDistanceAnalysis.getConversionfactor(t);
			}
			else
			{
				factors = Utils.newArray(length, 0.0, 1.0);
			}

			// Extract the mean-squared distance statistics
			Statistics[] stats = new Statistics[length];
			for (int i = 0; i < stats.length; i++)
				stats[i] = new Statistics();

			ArrayList<double[]> distances = (saveTraceDistances || displayTraceLength) ? new ArrayList<double[]>(
					traces.length) : null;

			// Store all the jump distances at the specified interval
			StoredDataStatistics jumpDistances = new StoredDataStatistics();
			final int jumpDistanceInterval = settings.jumpDistance;

			// Compute squared distances
			StoredDataStatistics msdPerMoleculeAllVsAll = new StoredDataStatistics();
			StoredDataStatistics msdPerMoleculeAdjacent = new StoredDataStatistics();
			for (Trace trace : traces)
			{
				ArrayList<PeakResult> results = trace.getPoints();
				// Sum the MSD and the time
				final int traceLength = (settings.truncate) ? settings.minimumTraceLength : trace.size();

				// Get the mean for each time separation
				double[] sumDistance = new double[traceLength + 1];
				double[] sumTime = new double[sumDistance.length];

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
							final double d = distance2(x, y, results.get(j));
							msd[j - 1] = px2ToUm2 * d;
							if (t == jumpDistanceInterval)
								jumpDistances.add(msd[j - 1]);
							sumDistance[t] += d;
							sumTime[t] += t;
						}
						distances.add(msd);
					}
					else
					{
						for (int j = 1; j < traceLength; j++)
						{
							final int t = j;
							final double d = distance2(x, y, results.get(j));
							if (t == jumpDistanceInterval)
								jumpDistances.add(px2ToUm2 * d);
							sumDistance[t] += d;
							sumTime[t] += t;
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
							final double d = distance2(x, y, results.get(j));
							if (t == jumpDistanceInterval)
								jumpDistances.add(px2ToUm2 * d);
							sumDistance[t] += d;
							sumTime[t] += t;
						}
					}

					// Add the average distance per time separation to the population
					for (int t = 1; t < traceLength; t++)
					{
						// Note: (traceLength - t) == count
						stats[t].add(sumDistance[t] / (traceLength - t));
					}
				}
				else
				{
					// Add the distance per time separation to the population
					for (int t = 1; t < traceLength; t++)
					{
						stats[t].add(sumDistance[t]);
					}
				}

				// Fix this for the precision and MSD adjustment.
				// It may be necessary to:
				// - sum the raw distances for each time interval (this is sumDistance[t])
				// - subtract the precision error
				// - apply correction factor for the n-frames to get actual MSD
				// - sum the actual MSD

				double sumD = 0, sumD_adjacent = Math.max(0, sumDistance[1] - error) * factors[1];
				double sumT = 0, sumT_adjacent = sumTime[1];
				for (int t = 1; t < traceLength; t++)
				{
					sumD += Math.max(0, sumDistance[t] - error) * factors[t];
					sumT += sumTime[t];
				}

				// Calculate the average displacement for the trace (do not simply use the largest 
				// time separation since this will miss moving molecules that end up at the origin)

				msdPerMoleculeAllVsAll.add(px2ToUm2PerSecond * sumD / sumT);
				msdPerMoleculeAdjacent.add(px2ToUm2PerSecond * sumD_adjacent / sumT_adjacent);
			}

			StoredDataStatistics dPerMoleculeAllVsAll = null;
			StoredDataStatistics dPerMoleculeAdjacent = null;
			if (saveTraceDistances || (settings.showHistograms && displayDHistogram))
			{
				dPerMoleculeAllVsAll = calculateDiffusionCoefficient(msdPerMoleculeAllVsAll);
				dPerMoleculeAdjacent = calculateDiffusionCoefficient(msdPerMoleculeAdjacent);
			}

			if (saveTraceDistances)
			{
				saveTraceDistances(traces.length, distances, msdPerMoleculeAllVsAll, msdPerMoleculeAdjacent,
						dPerMoleculeAllVsAll, dPerMoleculeAdjacent);
			}

			if (displayTraceLength)
			{
				StoredDataStatistics lengths = calculateTraceLengths(distances);
				showHistogram(lengths, "Trace length (um)");
			}

			// Plot the per-trace histogram of MSD and D
			if (settings.showHistograms)
			{
				if (displayMSDHistogram)
				{
					showHistogram(msdPerMoleculeAllVsAll, "MSD/Molecule (all-vs-all)");
					showHistogram(msdPerMoleculeAdjacent, "MSD/Molecule (adjacent)");
				}
				if (displayDHistogram)
				{
					showHistogram(dPerMoleculeAllVsAll, "D/Molecule (all-vs-all)");
					showHistogram(dPerMoleculeAdjacent, "D/Molecule (adjacent)");
				}
			}

			// Calculate the mean squared distance (MSD)
			double[] x = new double[stats.length];
			double[] y = new double[x.length];
			double[] sd = new double[x.length];
			// Intercept is the 4s^2 (in um^2)
			y[0] = 4 * precision * precision / 1e6;
			for (int i = 1; i < stats.length; i++)
			{
				x[i] = i * exposureTime;
				y[i] = stats[i].getMean() * px2ToUm2;
				//sd[i] = stats[i].getStandardDeviation() * px2ToUm2;
				sd[i] = stats[i].getStandardError() * px2ToUm2;
			}

			String title = TITLE + " MSD";
			Plot2 plot = plotMSD(x, y, sd, title);

			// Fit the MSD using a linear fit
			fitMSDResult = fitMSD(x, y, title, plot);

			// Jump Distance analysis
			if (saveRawData)
				saveStatistics(jumpDistances, "Jump Distance", "Distance (um^2)", false);

			// Calculate the cumulative jump-distance histogram
			double[][] jdHistogram = JumpDistanceAnalysis.cumulativeHistogram(jumpDistances.getValues());

			// Always show the jump distance histogram
			jdTitle = TITLE + " Jump Distance";
			jdPlot = new Plot2(jdTitle, "Distance (um^2)", "Cumulative Probability", jdHistogram[0], jdHistogram[1]);
			display(jdTitle, jdPlot);

			// Fit Jump Distance cumulative probability
			n = jumpDistances.getN();
			jumpDistanceParameters = jdParams = fitJumpDistance(jumpDistances, jdHistogram);
		}

		summarise(traces, fitMSDResult, n, jdParams);
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
	 * @param multi
	 */
	private void calculatePrecision(Trace[] traces, boolean multi)
	{
		// Check the diffusion simulation for a precision
		if (DiffusionRateTest.isSimulated(results.getName()) && !multi)
		{
			precision = DiffusionRateTest.lastSimulatedPrecision;
		}
		else
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
		}

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
	 * Calculate the diffusion coefficient (D) of the molecule. This is done by using the mean-squared deviation
	 * between frames divided by the time interval (delta) between frames. This is divided by 4 to produce the diffusion
	 * coefficient from two-dimensional distance analysis.
	 * <p>
	 * See Uphoff, et al, 2013. Single-molecule DNA repair in live bacteria, PNAS 110, 8063-8068
	 * 
	 * @param msdPerMoleculeAdjacent
	 * @return The D per molecule
	 */
	private StoredDataStatistics calculateDiffusionCoefficient(StoredDataStatistics msdPerMoleculeAdjacent)
	{
		StoredDataStatistics dPerMolecule = new StoredDataStatistics();
		final double diffusionCoefficientConversion = 1.0 / 4.0;
		for (double msd : msdPerMoleculeAdjacent.getValues())
		{
			dPerMolecule.add(msd * diffusionCoefficientConversion);
		}
		return dPerMolecule;
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
						.format("#TraceId\tMSD all-vs-all (um^2/s)\tMSD adjacent (um^2/s)\tD all-vs-all(um^2/s)\tD adjacent(um^2/s)\tDistances (um^2) per %ss ... ",
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

	private void summarise(Trace[] traces, double[] fitMSDResult, int n, double[][] jdParams)
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
		sb.append(settings.precisionCorrection).append("\t");
		sb.append(settings.mle).append("\t");
		sb.append(traces.length).append("\t");
		sb.append(Utils.rounded(precision, 4)).append("\t");
		double D = 0, s = 0;
		if (fitMSDResult != null)
		{
			D = fitMSDResult[0];
			s = fitMSDResult[1];
		}
		sb.append(Utils.rounded(D, 4)).append("\t");
		sb.append(Utils.rounded(s * 1000, 4)).append("\t");
		sb.append(Utils.rounded(settings.jumpDistance * exposureTime)).append("\t");
		sb.append(n).append("\t");
		sb.append(Utils.rounded(beta, 4)).append("\t");
		if (jdParams == null)
		{
			sb.append("\t\t\t");
		}
		else
		{
			sb.append(format(jdParams[0])).append("\t");
			sb.append(format(jdParams[1])).append("\t");
			sb.append(Utils.rounded(ic)).append("\t");
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
		StringBuilder sb = new StringBuilder("Title\tDataset\tExposure time (ms)\tD-threshold (nm)");
		sb.append("\tEx-threshold (nm)\t");
		sb.append("Min.Length\tIgnoreEnds\tTruncate\tInternal\tFit Length");
		sb.append("\tMSD corr.\ts corr.\tMLE\tTraces\ts (nm)\tD (um^2/s)\tfit s (nm)");
		sb.append("\tJump Distance (s)\tN\tBeta\tJump D (um^2/s)\tFractions\tIC");
		for (int i = 0; i < NAMES.length; i++)
		{
			sb.append("\t").append(NAMES[i]);
		}
		return sb.toString();
	}

	private boolean showTraceDialog(ArrayList<MemoryPeakResults> allResults)
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		if (!multiMode)
			ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);

		globalSettings = SettingsManager.loadSettings();
		settings = globalSettings.getClusteringSettings();

		gd.addNumericField("Distance_Threshold (nm)", settings.distanceThreshold, 0);
		gd.addNumericField("Distance_Exclusion (nm)", settings.distanceExclusion, 0);
		gd.addSlider("Min_trace_length", 2, 20, settings.minimumTraceLength);
		gd.addCheckbox("Ignore_ends", settings.ignoreEnds);
		gd.addCheckbox("Save_traces", settings.saveTraces);
		if (!multiMode)
			gd.addCheckbox("Multiple_inputs", multipleInputs);

		gd.showDialog();

		if (gd.wasCanceled() || !readTraceDialog(gd))
			return false;

		// Update the settings
		SettingsManager.saveSettings(globalSettings);

		// Load the results
		if (!multiMode)
		{
			MemoryPeakResults results = ResultsManager.loadInputResults(inputOption, true);
			if (results == null || results.size() == 0)
			{
				IJ.error(TITLE, "No results could be loaded");
				IJ.showStatus("");
				return false;
			}

			if (!checkCalibration(results))
				return false;

			allResults.add(results);
		}

		return true;
	}

	/**
	 * Check the results have a calibrated exposure time and pixel pitch. If not then show a dialog to collect the
	 * calibration.
	 * 
	 * @param results
	 * @return True if calibrated
	 */
	private boolean checkCalibration(MemoryPeakResults results)
	{
		if (results.getCalibration() == null || results.getCalibration().exposureTime <= 0 ||
				results.getCalibration().nmPerPixel <= 0)
		{
			Calibration cal = results.getCalibration();
			if (cal == null)
				cal = new Calibration();

			GenericDialog gd = new GenericDialog(TITLE);
			gd.addMessage("Uncalibrated results! Please enter the calibration:");
			gd.addNumericField("Exposure_time (ms)", cal.exposureTime, 2);
			gd.addNumericField("Pixel_pitch (nm)", cal.nmPerPixel, 2);
			gd.showDialog();
			if (gd.wasCanceled() || gd.invalidNumber())
				return false;
			cal.exposureTime = gd.getNextNumber();
			cal.nmPerPixel = gd.getNextNumber();
			if (cal.exposureTime <= 0 || cal.nmPerPixel <= 0)
				return false;
		}
		return true;
	}

	private boolean readTraceDialog(GenericDialog gd)
	{
		if (!multiMode)
			inputOption = ResultsManager.getInputSource(gd);
		settings.distanceThreshold = gd.getNextNumber();
		settings.distanceExclusion = Math.abs(gd.getNextNumber());
		settings.minimumTraceLength = (int) Math.abs(gd.getNextNumber());
		settings.ignoreEnds = gd.getNextBoolean();
		settings.saveTraces = gd.getNextBoolean();
		if (!multiMode)
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

	private Trace[] getTraces(ArrayList<MemoryPeakResults> allResults)
	{
		this.results = allResults.get(0);

		// Results should be checked for calibration by this point
		exposureTime = results.getCalibration().exposureTime / 1000;

		final double nmPerPixel = results.getCalibration().nmPerPixel;
		if (!multiMode && multipleInputs)
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
		gd.addCheckbox("Precision_correction", settings.precisionCorrection);
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
		settings.precisionCorrection = gd.getNextBoolean();
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
			gd.addCheckbox("D/Molecule", displayDHistogram);
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
	 * @return [D, precision]
	 */
	private double[] fitMSD(double[] x, double[] y, String title, Plot2 plot)
	{
		// The Weimann paper (Plos One e64287) fits:
		// MSD(n dt) = 4D n dt + 4s^2
		// n = number of jumps
		// dt = time difference between frames
		// s = localisation precision
		// Thus we should fit an intercept as well.

		// From the fit D = gradient / (4*exposureTime)

		double D = 0;
		double intercept = 0;
		double precision = 0;

		LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer();
		PointVectorValuePair lvmSolution;
		double ic = 0;

		// Fit with no intercept
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

			double gradient = lvmSolution.getPoint()[0];
			D = gradient / 4;

			Utils.log("Linear fit (%d points) : Gradient = %s, D = %s um^2/s, SS = %s, IC = %s (%d evaluations)",
					obs.length, Utils.rounded(gradient, 4), Utils.rounded(D, 4), Utils.rounded(ss), Utils.rounded(ic),
					optimizer.getEvaluations());
		}
		catch (TooManyIterationsException e)
		{
			Utils.log("Failed to fit : Too many iterations (%d)", optimizer.getIterations());
		}
		catch (ConvergenceException e)
		{
			Utils.log("Failed to fit : %s", e.getMessage());
		}

		// Fit with intercept.
		// Optionally include the intercept (which is the estimated precision).
		boolean fitIntercept = true;
		try
		{
			final LinearFunctionWithIntercept function = new LinearFunctionWithIntercept(x, y, settings.fitLength,
					fitIntercept);
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

			double ic2 = Maths.getInformationCriterion(ss, obs.length, 2);
			double gradient = lvmSolution.getPoint()[0];
			final double s = lvmSolution.getPoint()[1];
			double intercept2 = 4 * s * s;

			if (ic2 < ic || debugFitting)
			{
				// Convert fitted precision in um to nm
				Utils.log(
						"Linear fit with intercept (%d points) : Gradient = %s, Intercept = %s, D = %s um^2/s, precision = %s nm, SS = %s, IC = %s (%d evaluations)",
						obs.length, Utils.rounded(gradient, 4), Utils.rounded(intercept2, 4),
						Utils.rounded(gradient / 4, 4), Utils.rounded(s * 1000, 4), Utils.rounded(ss),
						Utils.rounded(ic2), optimizer.getEvaluations());
			}

			if (lvmSolution == null || ic2 < ic)
			{
				intercept = intercept2;
				D = gradient / 4;
				precision = s;
			}
		}
		catch (TooManyIterationsException e)
		{
			Utils.log("Failed to fit with intercept : Too many iterations (%d)", optimizer.getIterations());
		}
		catch (ConvergenceException e)
		{
			Utils.log("Failed to fit with intercept : %s", e.getMessage());
		}

		if (settings.msdCorrection)
		{
			// Fit with intercept including the MSD correction in the intercept.
			// For the MSD correction we fit including the correction factor (n-1/3)/n:
			// MSD = 4Dt n * (n - 1/3)/n + 4 s^2
			// MSD = 4Dt n - (4Dt) / 3 + 4 s^2
			// i.e. the intercept is allowed to be a small negative.
			try
			{
				// This function fits the jump distance (n) not the time (nt) so update x
				double[] x2 = new double[x.length];
				for (int i = 0; i < x2.length; i++)
					x2[i] = x[i] / exposureTime;

				final LinearFunctionWithMSDCorrectedIntercept function = new LinearFunctionWithMSDCorrectedIntercept(
						x2, y, settings.fitLength, fitIntercept);
				lvmSolution = optimizer.optimize(new MaxIter(3000), new MaxEval(Integer.MAX_VALUE),
						new ModelFunctionJacobian(new MultivariateMatrixFunction()
						{
							public double[][] value(double[] point) throws IllegalArgumentException
							{
								return function.jacobian(point);
							}
						}), new ModelFunction(function), new Target(function.getY()),
						new Weight(function.getWeights()), new InitialGuess(function.guess()));

				double ss = 0;
				double[] obs = function.getY();
				double[] exp = lvmSolution.getValue();
				for (int i = 0; i < obs.length; i++)
					ss += (obs[i] - exp[i]) * (obs[i] - exp[i]);

				double ic2 = Maths.getInformationCriterion(ss, obs.length, 2);
				double gradient = lvmSolution.getPoint()[0];
				final double s = lvmSolution.getPoint()[1];
				double intercept2 = 4 * s * s - gradient / 3;

				// Q. Is this working?
				// Try fixed precision fitting. Is the gradient correct?
				// Revisit all the equations to see if they are wrong.
				// Try adding the x[0] datapoint using the precision.
				// Change the formula to not be linear at x[0] and to just fit the precision, i.e. the intercept2 = 4 * s * s - gradient / 3 is wrong as the 
				// equation is not linear below n=1.

				// Incorporate the exposure time into the gradient to allow comparison to other fits 
				gradient /= exposureTime;

				if (ic2 < ic || debugFitting)
				{
					// Convert fitted precision in um to nm
					Utils.log(
							"Linear fit with MSD corrected intercept (%d points) : Gradient = %s, Intercept = %s, D = %s um^2/s, precision = %s nm, SS = %s, IC = %s (%d evaluations)",
							obs.length, Utils.rounded(gradient, 4), Utils.rounded(intercept2, 4),
							Utils.rounded(gradient / 4, 4), Utils.rounded(s * 1000, 4), Utils.rounded(ss),
							Utils.rounded(ic2), optimizer.getEvaluations());
				}

				if (lvmSolution == null || ic2 < ic)
				{
					intercept = intercept2;
					D = gradient / 4;
					precision = s;
				}
			}
			catch (TooManyIterationsException e)
			{
				Utils.log("Failed to fit with intercept : Too many iterations (%d)", optimizer.getIterations());
			}
			catch (ConvergenceException e)
			{
				Utils.log("Failed to fit with intercept : %s", e.getMessage());
			}
		}

		// Add the fit to the plot
		if (D > 0)
		{
			plot.setColor(Color.magenta);
			plot.drawLine(0, intercept, x[x.length - 1], 4 * D * x[x.length - 1] + intercept);
			display(title, plot);

			checkTraceDistance(D);
		}

		return new double[] { D, precision };
	}

	/**
	 * Check the distance used for tracing covers enough of the cumulative mean-squared distance distribution
	 * 
	 * @param d
	 */
	private void checkTraceDistance(double d)
	{
		double t = exposureTime;
		// Cumul P(r^2) = 1 - exp(-r^2 / 4dt)
		double r = settings.distanceThreshold / 1000;
		double msd = 4 * d * t;
		double p = 1 - FastMath.exp(-r * r / msd);
		Utils.log("Checking trace distance: r = %s nm, D = %s um^2/s, Cumul p(r^2|frame) = %s",
				settings.distanceThreshold, Utils.rounded(d), Utils.rounded(p));
		if (p < 0.95)
		{
			Utils.log("WARNING *** The tracing distance may not be large enough! ***");
		}
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
		final double[] x, y;
		final boolean fitIntercept;

		public LinearFunctionWithIntercept(double[] x, double[] y, int length, boolean fitIntercept)
		{
			this.fitIntercept = fitIntercept;
			int to = FastMath.min(x.length, 1 + length);
			// Optionally include the intercept
			int from = (fitIntercept) ? 0 : 1;
			this.x = Arrays.copyOfRange(x, from, to);
			this.y = Arrays.copyOfRange(y, from, to);
		}

		/**
		 * @return An estimate for the linear gradient and intercept
		 */
		public double[] guess()
		{
			int n1 = (fitIntercept) ? 1 : 0;

			if (y.length == n1 + 1)
				return new double[] { y[n1] / x[n1], 0 };

			double a = (y[y.length - 1] - y[n1]) / (x[x.length - 1] - x[n1]);
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

	public class LinearFunctionWithMSDCorrectedIntercept implements MultivariateVectorFunction
	{
		final double THIRD = 1 / 3.0;
		final double[] x, y;
		final boolean fitIntercept;

		public LinearFunctionWithMSDCorrectedIntercept(double[] x, double[] y, int length, boolean fitIntercept)
		{
			this.fitIntercept = fitIntercept;
			int to = FastMath.min(x.length, 1 + length);
			// Optionally include the intercept
			int from = (fitIntercept) ? 0 : 1;
			this.x = Arrays.copyOfRange(x, from, to);
			this.y = Arrays.copyOfRange(y, from, to);
		}

		/**
		 * @return An estimate for the linear gradient and intercept
		 */
		public double[] guess()
		{
			int n1 = (fitIntercept) ? 1 : 0;

			if (y.length == n1 + 1)
				return new double[] { y[n1] / x[n1], 0 };

			double a = (y[y.length - 1] - y[n1]) / (x[x.length - 1] - x[n1]);
			// y = ax - a/3 + 4c^2
			// y = ax + intercept
			// intercept = y - ax
			//           = 4c^2 - a/3
			// 4c^2 = intercept + a/3
			double intercept = y[y.length - 1] - a * x[x.length - 1];
			intercept += a * THIRD;
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
			// When x>=1:
			// y = ax - a/3 + 4c^2
			// When x==0:
			// y = 4c^2
			final double[] values = new double[x.length];
			final double a = variables[0];
			final double intercept = 4 * variables[1] * variables[1];
			final double error = intercept - a * THIRD;
			int i = 0;
			// Special case for fitting the intercept since the line is not linear below n=1
			if (fitIntercept)
				values[i++] = intercept;
			for (; i < values.length; i++)
			{
				values[i] = a * x[i] + error;
			}
			return values;
		}

		double[][] jacobian(double[] variables)
		{
			// Compute the gradients using calculus differentiation:
			// y = ax - a/3 + 4c^2
			// dy_da = x - 1/3
			// dy_dc = 8c
			double[][] jacobian = new double[x.length][2];
			final double dy_dc = 8 * variables[1];

			for (int i = 0; i < jacobian.length; ++i)
			{
				jacobian[i][0] = x[i] - THIRD;
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
	 *            (in um^2)
	 * @param jdHistogram
	 * @return The fitted coefficients and fractions
	 */
	private double[][] fitJumpDistance(StoredDataStatistics jumpDistances, double[][] jdHistogram)
	{
		final double msd = jumpDistances.getMean();
		final double meanDistance = Math.sqrt(msd) * 1e3;
		// TODO:
		// Q. Should the beta be expressed using the mean-distance or MSD? 
		// Q. Should it be normalised to the frame length. If not then the beta will be invariant on 
		// jump distance length
		beta = meanDistance / precision;
		Utils.log(
				"Jump Distance analysis : N = %d, Time = %d frames (%s seconds). MSD = %s um^2/jump, Mean Distance = %s nm/jump, Precision = %s nm, Beta = %s",
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

		// Set the calibration
		jd.setN(settings.jumpDistance);
		jd.setDeltaT(exposureTime);
		if (settings.precisionCorrection)
			jd.setError(precision, true);
		jd.setMsdCorrection(settings.msdCorrection);

		double[][] fit;
		if (settings.mle)
			fit = jd.fitJumpDistancesMLE(jumpDistances.getValues(), jdHistogram);
		else
			fit = jd.fitJumpDistanceHistogram(jumpDistances.getMean(), jdHistogram);

		// Get the raw fitted D and convert it to a calibrated D*
		if (fit != null)
		{
			fit[0] = jd.calculateApparentDiffusionCoefficient(fit[0]);
			// Check the largest D
			checkTraceDistance(fit[0][0]);
			ic = jd.getInformationCriterion();
		}

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

	/**
	 * Macro extension function.
	 * <p>
	 * Get the number of fitted species from the last call to fit the jump distances.
	 * 
	 * @param args
	 *            0: Double[1] - output the number of species
	 * @return Empty string
	 */
	public static String getNumberOfSpecies(Object[] args)
	{
		int n = 0;
		if (jumpDistanceParameters != null)
		{
			n = jumpDistanceParameters[0].length;
		}
		Double[] array = (Double[]) args[0];
		array[0] = new Double(n);
		return "";
	}

	/**
	 * Macro extension function.
	 * <p>
	 * Get the diffusion coefficient for the requested species from the last call to fit the jump distances.
	 * 
	 * @param args
	 *            0: Double[1] - input the index of the species; 1: Double[1] - output the coefficient
	 * @return Empty string
	 */
	public static String getD(Object[] args)
	{
		double value = 0;
		if (jumpDistanceParameters != null)
		{
			int i = ((Double) args[0]).intValue();
			if (i >= 0 && i < jumpDistanceParameters[0].length)
				value = jumpDistanceParameters[0][i];
		}
		((Double[]) args[1])[0] = new Double(value);
		return "";
	}

	/**
	 * Macro extension function.
	 * <p>
	 * Get the population fraction for the requested species from the last call to fit the jump distances.
	 * 
	 * @param args
	 *            0: Double[1] - input the index of the species; 1: Double[1] - output the population fraction
	 * @return Empty string
	 */
	public static String getF(Object[] args)
	{
		double value = 0;
		if (jumpDistanceParameters != null)
		{
			int i = ((Double) args[0]).intValue();
			if (i >= 0 && i < jumpDistanceParameters[1].length)
				value = jumpDistanceParameters[1][i];
		}
		((Double[]) args[1])[0] = new Double(value);
		return "";
	}

	/**
	 * Macro extension function.
	 * <p>
	 * Get the diffusion coefficient and population fraction for the requested species from the last call to fit the
	 * jump distances.
	 * 
	 * @param args
	 *            0: Double[1] - input the index of the species; 1: Double[1] - output the coefficient; 1: Double[1] -
	 *            output the population fraction
	 * @return Empty string
	 */
	public static String getSpecies(Object[] args)
	{
		double value = 0, value2 = 0;
		;
		if (jumpDistanceParameters != null)
		{
			int i = ((Double) args[0]).intValue();
			if (i >= 0 && i < jumpDistanceParameters[0].length)
			{
				value = jumpDistanceParameters[0][i];
				value2 = jumpDistanceParameters[1][i];
			}
		}
		((Double[]) args[1])[0] = new Double(value);
		((Double[]) args[2])[0] = new Double(value2);
		return "";
	}

	private class MultiDialog extends Dialog implements ActionListener, KeyListener, WindowListener
	{
		private static final long serialVersionUID = -881270633231897572L;

		private Button cancel, okay;
		private boolean wasCanceled;
		private List list;

		public MultiDialog()
		{
			super(WindowManager.getCurrentImage() != null ? (Frame) WindowManager.getCurrentImage().getWindow() : IJ
					.getInstance() != null ? IJ.getInstance() : new Frame(), TITLE, true);
			addKeyListener(this);
			addWindowListener(this);
		}

		public void showDialog()
		{
			add(buildPanels());
			if (IJ.isMacintosh())
				setResizable(false);
			pack();
			GUI.center(this);
			setVisible(true);
			IJ.wait(50); // work around for Sun/WinNT bug
		}

		protected Panel buildPanels()
		{
			Panel p = new Panel();
			BorderLayout layout = new BorderLayout();
			layout.setVgap(3);
			p.setLayout(layout);
			p.add(buildResultsList(), BorderLayout.NORTH, 0);
			p.add(buildButtonPanel(), BorderLayout.CENTER, 1);
			return p;
		}

		protected Component buildResultsList()
		{
			Set<String> names = MemoryPeakResults.getResultNames();
			final int MAX_SIZE = 30;
			int size;
			if (names.size() < MAX_SIZE)
			{
				size = names.size();
			}
			else
			{
				size = MAX_SIZE;
			}
			list = new List(size, true);
			int n = 0;
			for (String name : names)
			{
				list.add(name);
				// Select the same as last time
				if (selected.contains(name))
				{
					list.select(n);
				}
				n++;
			}
			return (Component) list;
		}

		protected Panel buildButtonPanel()
		{
			Panel buttons = new Panel();
			buttons.setLayout(new FlowLayout(FlowLayout.CENTER, 5, 0));
			okay = new Button("OK");
			okay.addActionListener(this);
			okay.addKeyListener(this);
			buttons.add(okay);
			cancel = new Button("Cancel");
			cancel.addActionListener(this);
			cancel.addKeyListener(this);
			buttons.add(cancel);
			return buttons;
		}

		public boolean wasCancelled()
		{
			return wasCanceled;
		}

		@Override
		public void actionPerformed(ActionEvent e)
		{
			Object source = e.getSource();
			if (source == okay || source == cancel)
			{
				wasCanceled = source == cancel;
				dispose();
			}
		}

		@Override
		public void keyTyped(KeyEvent paramKeyEvent)
		{
		}

		@Override
		public void keyPressed(KeyEvent e)
		{
			int keyCode = e.getKeyCode();
			IJ.setKeyDown(keyCode);
			if (keyCode == KeyEvent.VK_ENTER)
			{
				dispose();
			}
			else if (keyCode == KeyEvent.VK_ESCAPE)
			{
				wasCanceled = true;
				dispose();
				IJ.resetEscape();
			}
			else if (keyCode == KeyEvent.VK_W &&
					(e.getModifiers() & Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()) != 0)
			{
				wasCanceled = true;
				dispose();
			}
		}

		@Override
		public void keyReleased(KeyEvent paramKeyEvent)
		{
		}

		public ArrayList<String> getSelectedResults()
		{
			int[] listIndexes = list.getSelectedIndexes();
			selected.clear();

			// Record as if we use the multiple_inputs option in the trace dialog
			if (Recorder.record)
				Recorder.recordOption("Multiple_inputs");

			String name = list.getItem(listIndexes[0]);
			selected.add(name);
			if (Recorder.record)
				Recorder.recordOption("Input", name);

			for (int n = 1; n < listIndexes.length; ++n)
			{
				name = list.getItem(listIndexes[n]);
				selected.add(name);
				if (Recorder.record)
					Recorder.recordOption("Input"+n, name);
			}
			
			return selected;
		}

		@Override
		public void windowClosing(WindowEvent e)
		{
			wasCanceled = true;
			dispose();
		}

		//@formatter:off
	    public void windowActivated(WindowEvent e) {}
	    public void windowOpened(WindowEvent e) {}
	    public void windowClosed(WindowEvent e) {}
	    public void windowIconified(WindowEvent e) {}
	    public void windowDeiconified(WindowEvent e) {}
	    public void windowDeactivated(WindowEvent e) {}
		//@formatter:on
	}

	private boolean showMultiDialog(ArrayList<MemoryPeakResults> allResults)
	{
		multiMode = true;

		// Show a list box containing all the results. This should remember the last set of chosen items.
		MultiDialog md = new MultiDialog();

		md.showDialog();

		if (md.wasCancelled())
			return false;

		for (String name : md.getSelectedResults())
		{
			MemoryPeakResults r = MemoryPeakResults.getResults(name);
			if (r != null)
				allResults.add(r);
		}

		if (allResults.isEmpty())
			return false;

		// Check calibration exists for the first set of results
		if (!checkCalibration(allResults.get(0)))
			return false;

		// Check the calibration is the same for the rest
		Calibration cal = allResults.get(0).getCalibration();
		final double nmPerPixel = cal.nmPerPixel;
		final double exposureTime = cal.exposureTime;
		for (int i = 1; i < allResults.size(); i++)
		{
			MemoryPeakResults results = allResults.get(1);

			if (results.getCalibration() == null || results.getCalibration().exposureTime != exposureTime ||
					results.getNmPerPixel() != nmPerPixel)
			{
				IJ.error(TITLE, "The exposure time and pixel pitch must match across all the results");
				return false;
			}
		}

		return true;
	}
}
