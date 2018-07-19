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
package uk.ac.sussex.gdsc.smlm.ij.plugins;

import java.awt.Rectangle;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.util.FastMath;

import gnu.trove.set.hash.TIntHashSet;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.ExtendedGenericDialog;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.LutLoader;
import ij.plugin.PlugIn;
import ij.plugin.WindowOrganiser;
import ij.process.FloatProcessor;
import ij.text.TextWindow;
import uk.ac.sussex.gdsc.core.clustering.Cluster;
import uk.ac.sussex.gdsc.core.clustering.ClusterPoint;
import uk.ac.sussex.gdsc.core.clustering.ClusteringAlgorithm;
import uk.ac.sussex.gdsc.core.clustering.ClusteringEngine;
import uk.ac.sussex.gdsc.core.data.utils.Converter;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.ij.IJTrackProgress;
import uk.ac.sussex.gdsc.core.ij.Utils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationHelper;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CalibrationOrBuilder;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.ClusteringSettings;
import uk.ac.sussex.gdsc.smlm.data.config.UnitConverterFactory;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.TimeUnit;
import uk.ac.sussex.gdsc.smlm.engine.ParameterisedFitJob;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.Cluster.CentroidMethod;
import uk.ac.sussex.gdsc.smlm.results.ImageSource;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.TextFilePeakResults;
import uk.ac.sussex.gdsc.smlm.results.Trace;
import uk.ac.sussex.gdsc.smlm.results.TraceManager;
import uk.ac.sussex.gdsc.smlm.results.TraceManager.TraceMode;
import uk.ac.sussex.gdsc.smlm.results.count.Counter;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.PrecisionResultProcedure;

/**
 * Run a tracing algorithm on the peak results to trace molecules across the frames.
 */

public class TraceMolecules implements PlugIn
{
	private enum OptimiserPlot
	{
		//@formatter:off
		NONE{ @Override
		public String getName() { return "None"; }},
		NEAREST_NEIGHBOUR{ @Override
		public String getName() { return "Nearest neighbour"; }},
		BILINEAR{ @Override
		public String getName() { return "Bi-linear"; }};
		//@formatter:on

		@Override
		public String toString()
		{
			return getName();
		}

		/**
		 * Gets the name.
		 *
		 * @return the name
		 */
		abstract public String getName();

		public static OptimiserPlot get(int ordinal)
		{
			if (ordinal < 0 || ordinal >= values().length)
				ordinal = 0;
			return values()[ordinal];
		}
	}

	private static TraceMode getTraceMode(int traceMode)
	{
		if (traceMode < 0 || traceMode >= TraceMode.values().length)
			return TraceMode.LATEST_FORERUNNER;
		return TraceMode.values()[traceMode];
	}

	private static ClusteringAlgorithm getClusteringAlgorithm(int clusteringAlgorithm)
	{
		if (clusteringAlgorithm < 0 || clusteringAlgorithm >= ClusteringAlgorithm.values().length)
			return ClusteringAlgorithm.PAIRWISE;
		return ClusteringAlgorithm.values()[clusteringAlgorithm];
	}

	private String TITLE = "Trace or Cluster Molecules";
	private String outputName;
	private static double MIN_BLINKING_RATE = 1; // Should never be <= 0
	private static String inputOption = "";
	private static boolean inputDebugMode = true;
	private static boolean inputOptimiseBlinkingRate = false;

	//private static boolean fitOnlyCentroid = false;
	//private static float distanceThreshold = 1;
	//private static float expansionFactor = 2;
	//private static boolean debugFailures = false;

	private static String header = null;
	private static TextWindow summaryTable = null;

	private static final String[] NAMES = new String[] { "Total Signal", "Signal/Frame", "Blinks", "t-On (s)",
			"t-Off (s)", "Total t-On (s)", "Total t-Off (s)" };
	private static final String[] FILENAMES = new String[] { "total_signal", "signal_per_frame", "blinks", "t_on",
			"t_off", "total_t_on", "total_t_off" };
	private static boolean[] displayHistograms = new boolean[NAMES.length];
	static
	{
		for (int i = 0; i < displayHistograms.length; i++)
			displayHistograms[i] = true;
	}
	private static final int TOTAL_SIGNAL = 0;
	private static final int SIGNAL_PER_FRAME = 1;
	private static final int BLINKS = 2;
	private static final int T_ON = 3;
	private static final int T_OFF = 4;
	private static final int TOTAL_T_ON = 5;
	private static final int TOTAL_T_OFF = 6;

	private static boolean[] integerDisplay;
	static
	{
		integerDisplay = new boolean[NAMES.length];
		integerDisplay[BLINKS] = true;
		// Times are now in fractions of seconds
		//integerDisplay[T_ON] = true;
		//integerDisplay[T_OFF] = true;
		//integerDisplay[TOTAL_T_ON] = true;
		//integerDisplay[TOTAL_T_OFF] = true;
	}
	private static boolean[] alwaysRemoveOutliers;
	static
	{
		alwaysRemoveOutliers = new boolean[NAMES.length];
		alwaysRemoveOutliers[TOTAL_SIGNAL] = false;
	}

	private static String filename = "";

	private ClusteringSettings.Builder settings;
	private MemoryPeakResults results;
	// Store exposure time in seconds
	private double exposureTime = 0;

	// Used for the plotting
	private double[] dThresholds;
	private int[] tThresholds;
	private ArrayList<double[]> zeroCrossingPoints;
	private FloatProcessor fp;
	private Calibration cal;
	// Store the pixel value for the first plotted result
	private int origX, origY;
	private boolean debugMode = false, altKeyDown, optimiseBlinkingRate = false;

	/*
	 * (non-Javadoc)
	 *
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	@Override
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (MemoryPeakResults.isMemoryEmpty())
		{
			IJ.error(TITLE, "No localisations in memory");
			return;
		}
		altKeyDown = Utils.isExtraOptions();

		Trace[] traces = null;
		int totalFiltered = 0;
		if ("cluster".equals(arg))
		{
			// --=-=-=-=-=-
			// Clustering
			// --=-=-=-=-=-
			outputName = "Cluster";

			if (!showClusterDialog())
				return;

			final ClusteringEngine engine = new ClusteringEngine(Prefs.getThreads(),
					getClusteringAlgorithm(settings.getClusteringAlgorithm()), new IJTrackProgress());

			if (settings.getSplitPulses())
			{
				engine.setPulseInterval(settings.getPulseInterval());
				limitTimeThreshold(settings.getPulseInterval());
			}

			final ArrayList<Cluster> clusters = engine.findClusters(convertToClusterPoints(),
					getDistance(settings.getDistanceThreshold(), results.getCalibration()), timeThresholdInFrames());

			if (clusters == null)
			{
				Utils.log("Aborted");
				return;
			}

			traces = convertToTraces(clusters);
		}
		else
		{
			// --=-=-=-=-=-
			// Tracing
			// --=-=-=-=-=-
			outputName = "Trace";

			if (!showDialog())
				return;

			final TraceManager manager = new TraceManager(results);
			manager.setTraceMode(getTraceMode(settings.getTraceMode()));
			manager.setActivationFrameInterval(settings.getPulseInterval());
			manager.setActivationFrameWindow(settings.getPulseWindow());
			manager.setDistanceExclusion(getDistance(settings.getDistanceExclusion(), results.getCalibration()));

			if (settings.getOptimise())
				// Optimise before configuring for a pulse interval
				runOptimiser(manager);

			if (settings.getSplitPulses())
			{
				manager.setPulseInterval(settings.getPulseInterval());
				limitTimeThreshold(settings.getPulseInterval());
			}

			manager.setTracker(new IJTrackProgress());
			manager.traceMolecules(getDistance(settings.getDistanceThreshold(), results.getCalibration()),
					timeThresholdInFrames());
			traces = manager.getTraces();
			totalFiltered = manager.getTotalFiltered();
		}

		// --=-=-=-=-=-
		// Results processing
		// --=-=-=-=-=-

		outputName += (outputName.endsWith("e") ? "" : "e") + "d";
		saveResults(results, traces, outputName);

		// Save singles + single localisations in a trace
		saveCentroidResults(results, getSingles(traces), outputName + " Singles");
		final Trace[] multiTraces = getTraces(traces);
		saveResults(results, multiTraces, outputName + " Multi");

		// Save centroids
		outputName += " Centroids";
		final MemoryPeakResults tracedResults = saveCentroidResults(results, traces, outputName);

		// Save traces separately
		saveCentroidResults(results, multiTraces, outputName + " Multi");

		// Sort traces by time to assist the results source in extracting frames sequentially.
		// Do this before saving to assist in debugging using the saved traces file.
		sortByTime(traces);

		if (settings.getSaveTraces())
			saveTraces(traces);

		summarise(traces, totalFiltered, settings.getDistanceThreshold(), timeThresholdInSeconds());

		IJ.showStatus(String.format("%d localisations => %d traces (%d filtered)", results.size(), tracedResults.size(),
				totalFiltered));

		//// Provide option to refit the traces as single peaks and save to memory
		//if (settings.refitOption)
		//	fitTraces(results, traces);
	}

	private static double getDistance(double distanceThreshold, CalibrationOrBuilder calibration)
	{
		// Convert from NM to native units
		final Converter c = CalibrationHelper.getDistanceConverter(calibration, DistanceUnit.NM);
		return c.convertBack(distanceThreshold);
	}

	/**
	 * Limit the time threshold to the pulse interval duration.
	 *
	 * @param pulseInterval
	 *            the pulse interval
	 */
	private void limitTimeThreshold(int pulseInterval)
	{
		// The pulse interval is in frames.
		// Convert the interval to the correct units.

		double limit;
		if (settings.getTimeUnit() == TimeUnit.FRAME)
			limit = settings.getPulseInterval();
		else
		{
			final TypeConverter<TimeUnit> convert = UnitConverterFactory.createConverter(TimeUnit.FRAME,
					settings.getTimeUnit(), exposureTime);
			limit = convert.convert(pulseInterval);
		}

		if (settings.getTimeThreshold() > limit)
			settings.setTimeThreshold(limit);
	}

	private List<ClusterPoint> convertToClusterPoints()
	{
		return convertToClusterPoints(results);
	}

	/**
	 * Convert a list of peak results into points for the clustering engine.
	 *
	 * @param results
	 *            the results
	 * @return the list of clusters
	 */
	public static List<ClusterPoint> convertToClusterPoints(MemoryPeakResults results)
	{
		final ArrayList<ClusterPoint> points = new ArrayList<>(results.size());
		final Counter counter = new Counter();
		results.forEach(new PeakResultProcedure()
		{
			@Override
			public void execute(PeakResult p)
			{
				points.add(ClusterPoint.newTimeClusterPoint(counter.getAndIncrement(), p.getXPosition(),
						p.getYPosition(), p.getIntensity(), p.getFrame(), p.getEndFrame()));
			}
		});
		return points;
	}

	private Trace[] convertToTraces(ArrayList<Cluster> clusters)
	{
		return convertToTraces(results, clusters);
	}

	/**
	 * Convert the clusters from the clustering engine into traces composed of the original list of peak results.
	 *
	 * @param results
	 *            the results
	 * @param clusters
	 *            the clusters
	 * @return the traces
	 */
	public static Trace[] convertToTraces(MemoryPeakResults results, ArrayList<Cluster> clusters)
	{
		final Trace[] traces = new Trace[clusters.size()];
		int i = 0;
		for (final Cluster cluster : clusters)
		{
			final Trace trace = new Trace();
			trace.setId(i + 1);
			for (ClusterPoint point = cluster.head; point != null; point = point.next)
				// The point Id was the position in the original results array
				trace.add(results.get(point.id));
			traces[i++] = trace;
		}
		return traces;
	}

	/**
	 * Sort traces by time.
	 *
	 * @param traces
	 *            the traces
	 */
	static void sortByTime(Trace[] traces)
	{
		for (final Trace t : traces)
			t.sort();
		Arrays.sort(traces, new Comparator<Trace>()
		{
			@Override
			public int compare(Trace o1, Trace o2)
			{
				return o1.getHead().getFrame() - o2.getHead().getFrame();
			}
		});
	}

	/**
	 * Convert the traces to results.
	 *
	 * @param sourceResults
	 *            the source results
	 * @param traces
	 *            the traces
	 * @param name
	 *            the name
	 * @return the memory peak results
	 */
	static MemoryPeakResults saveResults(MemoryPeakResults sourceResults, Trace[] traces, String name)
	{
		final MemoryPeakResults tracedResults = TraceManager.convertToPeakResults(sourceResults, traces);
		tracedResults.setName(sourceResults.getName() + " " + name);
		MemoryPeakResults.addResults(tracedResults);
		return tracedResults;
	}

	private static MemoryPeakResults saveCentroidResults(MemoryPeakResults sourceResults, Trace[] traces, String name)
	{
		final MemoryPeakResults tracedResults = TraceManager.convertToCentroidPeakResults(sourceResults, traces);
		tracedResults.setName(sourceResults.getName() + " " + name);
		MemoryPeakResults.addResults(tracedResults);
		return tracedResults;
	}

	private static Trace[] getSingles(Trace[] traces)
	{
		final ArrayList<Trace> result = new ArrayList<>();
		for (final Trace t : traces)
			if (t.size() == 1)
				result.add(t);
		return result.toArray(new Trace[result.size()]);
	}

	private static Trace[] getTraces(Trace[] traces)
	{
		final ArrayList<Trace> result = new ArrayList<>();
		for (final Trace t : traces)
			if (t.size() != 1)
				result.add(t);
		return result.toArray(new Trace[result.size()]);
	}

	private void saveTraces(Trace[] traces)
	{
		filename = saveTraces(results, traces, createSettingsComment(), filename, 0);
	}

	/**
	 * Save the traces to the file. A File open dialog is presented and the selected filename returned.
	 * <p>
	 * If the id is above zero then the file open dialog title will have the id appended and the filename is searched
	 * for .[0-9]+. and it is replaced with .id.
	 *
	 * @param sourceResults
	 *            the source results
	 * @param traces
	 *            the traces
	 * @param comment
	 *            the comment
	 * @param filename
	 *            The initial filename
	 * @param id
	 *            The traces id (used if above 0)
	 * @return The select filename (or null)
	 */
	static String saveTraces(MemoryPeakResults sourceResults, Trace[] traces, String comment, String filename, int id)
	{
		IJ.showStatus("Saving traces");
		String title = "Traces_File";
		if (id > 0)
		{
			title += id;
			final String regex = "\\.[0-9]+\\.";
			if (filename.matches(regex))
				filename.replaceAll(regex, "." + (id) + ".");
			else
				Utils.replaceExtension(filename, id + ".xls");
		}
		filename = Utils.getFilename(title, filename);
		if (filename != null)
		{
			filename = Utils.replaceExtension(filename, "xls");

			final boolean showDeviations = sourceResults.hasDeviations();
			// Assume that are results are from a single frame but store the trace ID
			final TextFilePeakResults traceResults = new TextFilePeakResults(filename, showDeviations, false, true);
			traceResults.copySettings(sourceResults);
			traceResults.begin();
			if (!traceResults.isActive())
				IJ.error("Failed to write to file: " + filename);
			else
			{
				traceResults.addComment(comment);
				for (final Trace trace : traces)
					traceResults.addTrace(trace);
				traceResults.end();
			}
		}
		IJ.showStatus("");
		return filename;
	}

	private String createSettingsComment()
	{
		return String.format("Molecule tracing : distance-threshold = %f : time-threshold = %f (%d frames)",
				settings.getDistanceThreshold(), timeThresholdInSeconds(), timeThresholdInFrames());
	}

	private void summarise(Trace[] traces, int filtered, double dThreshold, double tThreshold)
	{
		IJ.showStatus("Calculating summary ...");

		// Create summary table
		createSummaryTable();

		final Statistics[] stats = new Statistics[NAMES.length];
		for (int i = 0; i < stats.length; i++)
			stats[i] = (settings.getShowHistograms() || settings.getSaveTraceData()) ? new StoredDataStatistics()
					: new Statistics();
		int singles = 0;
		for (final Trace trace : traces)
		{
			final int nBlinks = trace.getNBlinks() - 1;
			stats[BLINKS].add(nBlinks);
			final int[] onTimes = trace.getOnTimes();
			final int[] offTimes = trace.getOffTimes();
			double tOn = 0;
			for (final int t : onTimes)
			{
				stats[T_ON].add(t * exposureTime);
				tOn += t * exposureTime;
			}
			stats[TOTAL_T_ON].add(tOn);
			if (offTimes != null)
			{
				double tOff = 0;
				for (final int t : offTimes)
				{
					stats[T_OFF].add(t * exposureTime);
					tOff += t * exposureTime;
				}
				stats[TOTAL_T_OFF].add(tOff);
			}
			final double signal = trace.getSignal() / results.getGain();
			stats[TOTAL_SIGNAL].add(signal);
			stats[SIGNAL_PER_FRAME].add(signal / trace.size());
			if (trace.size() == 1)
				singles++;
		}

		// Add to the summary table
		final StringBuilder sb = new StringBuilder();
		sb.append(results.getName()).append('\t');
		sb.append(outputName.equals("Cluster") ? getClusteringAlgorithm(settings.getClusteringAlgorithm())
				: getTraceMode(settings.getTraceMode())).append('\t');
		sb.append(Utils.rounded(exposureTime * 1000, 3)).append('\t');
		sb.append(Utils.rounded(dThreshold, 3)).append('\t');
		sb.append(Utils.rounded(tThreshold, 3));
		if (settings.getSplitPulses())
			sb.append(" *");
		sb.append('\t');
		sb.append(convertSecondsTotFrames(tThreshold)).append('\t');
		sb.append(traces.length).append('\t');
		sb.append(filtered).append('\t');
		sb.append(singles).append('\t');
		sb.append(traces.length - singles).append('\t');
		for (int i = 0; i < stats.length; i++)
			sb.append(Utils.rounded(stats[i].getMean(), 3)).append('\t');
		if (java.awt.GraphicsEnvironment.isHeadless())
		{
			IJ.log(sb.toString());
			return;
		}
		summaryTable.append(sb.toString());

		if (settings.getShowHistograms())
		{
			IJ.showStatus("Calculating histograms ...");

			int[] idList = new int[NAMES.length];
			int count = 0;

			boolean requireRetile = false;
			for (int i = 0; i < NAMES.length; i++)
				if (displayHistograms[i])
				{
					idList[count++] = Utils.showHistogram(TITLE, (StoredDataStatistics) stats[i], NAMES[i],
							(integerDisplay[i]) ? 1 : 0,
							(settings.getRemoveOutliers() || alwaysRemoveOutliers[i]) ? 2 : 0,
							settings.getHistogramBins());
					requireRetile = requireRetile || Utils.isNewWindow();
				}

			if (count > 0 && requireRetile)
			{
				idList = Arrays.copyOf(idList, count);
				new WindowOrganiser().tileWindows(idList);
			}
		}

		if (settings.getSaveTraceData())
			saveTraceData(stats);

		IJ.showStatus("");
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
		else if (summaryTable == null || !summaryTable.isVisible())
		{
			summaryTable = new TextWindow(TITLE + " Data Summary", createHeader(), "", 800, 300);
			summaryTable.setVisible(true);
		}
	}

	private static String createHeader()
	{
		final StringBuilder sb = new StringBuilder(
				"Dataset\tAlgorithm\tExposure time (ms)\tD-threshold (nm)\tT-threshold (s)\t(Frames)\tMolecules\tFiltered\tSingles\tClusters");
		for (int i = 0; i < NAMES.length; i++)
			sb.append('\t').append(NAMES[i]);
		return sb.toString();
	}

	private void saveTraceData(Statistics[] stats)
	{
		// Get the directory
		IJ.showStatus("Saving trace data");
		final String directory = Utils.getDirectory("Trace_data_directory", settings.getTraceDataDirectory());
		if (directory != null)
		{
			settings.setTraceDataDirectory(directory);
			SettingsManager.writeSettings(settings.build());
			for (int i = 0; i < NAMES.length; i++)
				saveTraceData((StoredDataStatistics) stats[i], NAMES[i], FILENAMES[i]);
		}
		IJ.showStatus("");
	}

	private void saveTraceData(StoredDataStatistics s, String name, String fileSuffix)
	{
		try (BufferedWriter file = new BufferedWriter(
				new FileWriter(settings.getTraceDataDirectory() + TITLE + "." + fileSuffix + ".txt")))
		{
			file.append(name);
			file.newLine();

			for (final double d : s.getValues())
			{
				file.append(Utils.rounded(d, 4));
				file.newLine();
			}
		}
		catch (final Exception e)
		{
			// Q. Add better handling of errors?
			e.printStackTrace();
			IJ.log("Failed to save trace data to results directory: " + settings.getTraceDataDirectory());
		}
	}

	private boolean showDialog()
	{
		TITLE = outputName + " Molecules";
		final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);

		settings = SettingsManager.readClusteringSettings(0).toBuilder();

		gd.addNumericField("Distance_Threshold", settings.getDistanceThreshold(), 2, 6, "nm");
		gd.addNumericField("Distance_Exclusion", settings.getDistanceExclusion(), 2, 6, "nm");
		gd.addNumericField("Time_Threshold", settings.getTimeThreshold(), 2);
		gd.addChoice("Time_unit", SettingsManager.getTimeUnitNames(), settings.getTimeUnit().ordinal());
		final String[] traceModes = SettingsManager.getNames((Object[]) TraceManager.TraceMode.values());
		gd.addChoice("Trace_mode", traceModes, traceModes[getTraceMode(settings.getTraceMode()).ordinal()]);
		gd.addNumericField("Pulse_interval", settings.getPulseInterval(), 0, 6, "Frames");
		gd.addNumericField("Pulse_window", settings.getPulseWindow(), 0, 6, "Frames");
		gd.addCheckbox("Split_pulses", settings.getSplitPulses());
		gd.addCheckbox("Optimise", settings.getOptimise());
		gd.addCheckbox("Save_traces", settings.getSaveTraces());
		gd.addCheckbox("Show_histograms", settings.getShowHistograms());
		gd.addCheckbox("Save_trace_data", settings.getSaveTraceData());
		//gd.addCheckbox("Refit_option", settings.refitOption);
		if (altKeyDown)
			gd.addCheckbox("Debug", inputDebugMode);

		gd.showDialog();

		if (gd.wasCanceled() || !readDialog(gd))
			return false;

		// Update the settings
		SettingsManager.writeSettings(settings.build());

		// Load the results
		results = ResultsManager.loadInputResults(inputOption, true, null, null);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			IJ.showStatus("");
			return false;
		}

		// Store exposure time in seconds
		exposureTime = results.getCalibrationReader().getExposureTime() / 1000;

		return true;
	}

	private boolean readDialog(ExtendedGenericDialog gd)
	{
		inputOption = ResultsManager.getInputSource(gd);
		settings.setDistanceThreshold(gd.getNextNumber());
		settings.setDistanceExclusion(gd.getNextNumber());
		settings.setTimeThreshold(gd.getNextNumber());
		settings.setTimeUnit(SettingsManager.getTimeUnitValues()[gd.getNextChoiceIndex()]);
		settings.setTraceMode(gd.getNextChoiceIndex());
		settings.setPulseInterval((int) gd.getNextNumber());
		settings.setPulseWindow((int) gd.getNextNumber());
		settings.setSplitPulses(gd.getNextBoolean());
		settings.setOptimise(gd.getNextBoolean());
		settings.setSaveTraces(gd.getNextBoolean());
		settings.setShowHistograms(gd.getNextBoolean());
		settings.setSaveTraceData(gd.getNextBoolean());
		//settings.refitOption = gd.getNextBoolean();
		if (altKeyDown)
			debugMode = inputDebugMode = gd.getNextBoolean();

		if (gd.invalidNumber())
			return false;

		if (settings.getShowHistograms())
		{
			gd = new ExtendedGenericDialog(TITLE);
			gd.addMessage("Select the histograms to display");
			gd.addCheckbox("Remove_outliers", settings.getRemoveOutliers());
			gd.addNumericField("Histogram_bins", settings.getHistogramBins(), 0);
			for (int i = 0; i < displayHistograms.length; i++)
				gd.addCheckbox(NAMES[i].replace(' ', '_'), displayHistograms[i]);
			gd.showDialog();
			if (gd.wasCanceled())
				return false;
			settings.setRemoveOutliers(gd.getNextBoolean());
			settings.setHistogramBins((int) Math.abs(gd.getNextNumber()));
			for (int i = 0; i < displayHistograms.length; i++)
				displayHistograms[i] = gd.getNextBoolean();
		}

		// Check arguments
		try
		{
			Parameters.isAboveZero("Distance threshold", settings.getDistanceThreshold());
			Parameters.isAboveZero("Time threshold", settings.getTimeThreshold());
			Parameters.isPositive("Pulse interval", settings.getPulseInterval());
			Parameters.isPositive("Pulse window", settings.getPulseWindow());
			//Parameters.isAboveZero("Histogram bins", settings.getHistogramBins());
		}
		catch (final IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		return true;
	}

	private boolean showClusterDialog()
	{
		TITLE = outputName + " Molecules";
		final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);

		settings = SettingsManager.readClusteringSettings(0).toBuilder();

		gd.addNumericField("Distance_Threshold", settings.getDistanceThreshold(), 2, 6, "nm");
		gd.addNumericField("Time_Threshold", settings.getTimeThreshold(), 2);
		gd.addChoice("Time_unit", SettingsManager.getTimeUnitNames(), settings.getTimeUnit().ordinal());
		final String[] algorithm = SettingsManager.getNames((Object[]) ClusteringAlgorithm.values());
		gd.addChoice("Clustering_algorithm", algorithm,
				algorithm[getClusteringAlgorithm(settings.getClusteringAlgorithm()).ordinal()]);
		gd.addNumericField("Pulse_interval", settings.getPulseInterval(), 0, 6, "Frames");
		gd.addCheckbox("Split_pulses", settings.getSplitPulses());
		gd.addCheckbox("Save_clusters", settings.getSaveTraces());
		gd.addCheckbox("Show_histograms", settings.getShowHistograms());
		gd.addCheckbox("Save_cluster_data", settings.getSaveTraceData());
		gd.addCheckbox("Refit_option", settings.getRefitOption());
		if (altKeyDown)
			gd.addCheckbox("Debug", inputDebugMode);

		gd.showDialog();

		if (gd.wasCanceled() || !readClusterDialog(gd))
			return false;

		// Update the settings
		SettingsManager.writeSettings(settings.build());

		// Load the results
		results = ResultsManager.loadInputResults(inputOption, true, null, null);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			IJ.showStatus("");
			return false;
		}

		// Store exposure time in seconds
		exposureTime = results.getCalibrationReader().getExposureTime() / 1000;

		return true;
	}

	private boolean readClusterDialog(ExtendedGenericDialog gd)
	{
		inputOption = ResultsManager.getInputSource(gd);
		settings.setDistanceThreshold(gd.getNextNumber());
		settings.setTimeThreshold(gd.getNextNumber());
		settings.setTimeUnit(SettingsManager.getTimeUnitValues()[gd.getNextChoiceIndex()]);
		settings.setClusteringAlgorithm(gd.getNextChoiceIndex());
		settings.setPulseInterval((int) gd.getNextNumber());
		settings.setSplitPulses(gd.getNextBoolean());
		settings.setSaveTraces(gd.getNextBoolean());
		settings.setShowHistograms(gd.getNextBoolean());
		settings.setSaveTraceData(gd.getNextBoolean());
		settings.setRefitOption(gd.getNextBoolean());
		if (altKeyDown)
			debugMode = inputDebugMode = gd.getNextBoolean();

		if (gd.invalidNumber())
			return false;

		if (settings.getShowHistograms())
		{
			gd = new ExtendedGenericDialog(TITLE);
			gd.addMessage("Select the histograms to display");
			gd.addCheckbox("Remove_outliers", settings.getRemoveOutliers());
			gd.addNumericField("Histogram_bins", settings.getHistogramBins(), 0);
			for (int i = 0; i < displayHistograms.length; i++)
				gd.addCheckbox(NAMES[i].replace(' ', '_'), displayHistograms[i]);
			gd.showDialog();
			if (gd.wasCanceled())
				return false;
			settings.setRemoveOutliers(gd.getNextBoolean());
			settings.setHistogramBins((int) Math.abs(gd.getNextNumber()));
			for (int i = 0; i < displayHistograms.length; i++)
				displayHistograms[i] = gd.getNextBoolean();
		}

		// Check arguments
		try
		{
			Parameters.isAboveZero("Distance threshold", settings.getDistanceThreshold());
			final ClusteringAlgorithm clusteringAlgorithm = getClusteringAlgorithm(settings.getClusteringAlgorithm());
			if (clusteringAlgorithm == ClusteringAlgorithm.CENTROID_LINKAGE_DISTANCE_PRIORITY ||
					clusteringAlgorithm == ClusteringAlgorithm.CENTROID_LINKAGE_TIME_PRIORITY)
			{
				Parameters.isAboveZero("Time threshold", settings.getTimeThreshold());
				Parameters.isPositive("Pulse interval", settings.getPulseInterval());
			}
			//Parameters.isAboveZero("Histogram bins", settings.getHistogramBins());
		}
		catch (final IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		return true;
	}

	private void runOptimiser(TraceManager manager)
	{
		// Get an estimate of the number of molecules without blinking
		final Statistics stats = new Statistics();
		final double nmPerPixel = this.results.getNmPerPixel();
		final PrecisionResultProcedure pp = new PrecisionResultProcedure(results);
		pp.getPrecision();
		stats.add(pp.precision);
		// Use twice the precision to get the initial distance threshold

		// Use 2.5x sigma as per the PC-PALM protocol in Sengupta, et al (2013) Nature Protocols 8, 345
		final double dEstimate = stats.getMean() * 2.5 / nmPerPixel;
		final int n = manager.traceMolecules(dEstimate, 1);
		//for (double d : new double[] { 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4 })
		//	System.out.printf("d=%.2f, estimate=%d\n", d,
		//			manager.traceMolecules(stats.getMean() * d / this.results.getNmPerPixel(), 1));

		if (!getParameters(n, dEstimate))
			return;

		// TODO - Convert the distance threshold to use nm instead of pixels?
		final List<double[]> results = runTracing(manager, settings.getMinDistanceThreshold(),
				settings.getMaxDistanceThreshold(), settings.getMinTimeThreshold(), settings.getMaxTimeThreshold(),
				settings.getOptimiserSteps());

		// Compute fractional difference from the true value:
		// Use blinking rate directly or the estimated number of molecules
		double nReference;
		int statistic;
		if (optimiseBlinkingRate)
		{
			nReference = settings.getBlinkingRate();
			statistic = 3;
			IJ.log(String.format("Estimating blinking rate: %.2f", nReference));
		}
		else
		{
			nReference = n / settings.getBlinkingRate();
			statistic = 2;
			IJ.log(String.format("Estimating number of molecules: %d / %.2f = %.2f", n, settings.getBlinkingRate(),
					nReference));
		}

		for (final double[] result : results)
			//System.out.printf("%g %g = %g\n", result[0], result[1], result[2]);
			if (optimiseBlinkingRate)
				result[2] = (nReference - result[statistic]) / nReference;
			else
				result[2] = (result[statistic] - nReference) / nReference;

		// Locate the optimal parameters with a fit of the zero contour
		final boolean found = findOptimalParameters(results);

		createPlotResults(results);

		if (!found)
			return;

		// Make fractional difference absolute so that lowest is best
		for (final double[] result : results)
			result[2] = Math.abs(result[2]);

		// Set the optimal thresholds using the lowest value
		double[] best = new double[] { 0, 0, Double.MAX_VALUE };
		for (final double[] result : results)
			if (best[2] > result[2])
				best = result;

		settings.setDistanceThreshold(best[0]);

		// The optimiser works using frames so convert back to the correct units
		final TypeConverter<TimeUnit> convert = UnitConverterFactory.createConverter(TimeUnit.FRAME,
				settings.getTimeUnit(), exposureTime);
		settings.setTimeThreshold(convert.convert(best[1]));

		IJ.log(String.format("Optimal fractional difference @ D-threshold=%g, T-threshold=%f (%d frames)",
				settings.getDistanceThreshold(), timeThresholdInSeconds(), timeThresholdInFrames()));
		SettingsManager.writeSettings(settings.build());
	}

	private boolean getParameters(int n, double d)
	{
		final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE + " Optimiser");

		final String msg = String.format("Estimate %d molecules at d=%f, t=1", n, d);
		IJ.log(msg);
		gd.addMessage(msg);
		gd.addNumericField("Min_Distance_Threshold (px)", settings.getMinDistanceThreshold(), 2);
		gd.addNumericField("Max_Distance_Threshold (px)", settings.getMaxDistanceThreshold(), 2);
		gd.addNumericField("Min_Time_Threshold (frames)", settings.getMinTimeThreshold(), 0);
		gd.addNumericField("Max_Time_Threshold (frames)", settings.getMaxTimeThreshold(), 0);
		gd.addSlider("Steps", 1, 20, settings.getOptimiserSteps());
		gd.addNumericField("Blinking_rate", settings.getBlinkingRate(), 2);
		final String[] plotNames = SettingsManager.getNames((Object[]) OptimiserPlot.values());
		gd.addChoice("Plot", plotNames, plotNames[OptimiserPlot.get(settings.getOptimiserPlot()).ordinal()]);
		if (altKeyDown)
			gd.addCheckbox("Optimise_blinking", inputOptimiseBlinkingRate);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		settings.setMinDistanceThreshold(gd.getNextNumber());
		settings.setMaxDistanceThreshold(gd.getNextNumber());
		settings.setMinTimeThreshold((int) gd.getNextNumber());
		settings.setMaxTimeThreshold((int) gd.getNextNumber());
		settings.setOptimiserSteps((int) gd.getNextNumber());
		settings.setBlinkingRate(gd.getNextNumber());
		settings.setOptimiserPlot(gd.getNextChoiceIndex());
		if (altKeyDown)
			optimiseBlinkingRate = inputOptimiseBlinkingRate = gd.getNextBoolean();

		if (gd.invalidNumber())
			return false;

		if (settings.getMinDistanceThreshold() < 0)
			settings.setMinDistanceThreshold(0);
		if (settings.getMaxDistanceThreshold() < settings.getMinDistanceThreshold())
			settings.setMaxDistanceThreshold(settings.getMinDistanceThreshold());
		if (settings.getMinTimeThreshold() < 0)
			settings.setMinTimeThreshold(0);
		if (settings.getMaxTimeThreshold() < settings.getMinTimeThreshold())
			settings.setMaxTimeThreshold(settings.getMinTimeThreshold());
		if (settings.getOptimiserSteps() < 0)
			settings.setOptimiserSteps(1);
		if (settings.getBlinkingRate() < MIN_BLINKING_RATE)
		{
			IJ.error(gd.getTitle(), "Blinking rate must be above " + MIN_BLINKING_RATE);
			return false;
		}

		if (settings.getMinDistanceThreshold() == settings.getMaxDistanceThreshold() &&
				settings.getMinTimeThreshold() == settings.getMaxTimeThreshold())
		{
			IJ.error(gd.getTitle(), "Nothing to optimise");
			return false;
		}

		SettingsManager.writeSettings(settings.build());

		return true;
	}

	private int convertSecondsTotFrames(double timeInSeconds)
	{
		return (int) Math.round(timeInSeconds / exposureTime);
	}

	private int timeThresholdInFrames()
	{
		return (int) Math.round(timeThresholdIn(TimeUnit.FRAME));
	}

	private double timeThresholdInSeconds()
	{
		return timeThresholdIn(TimeUnit.SECOND);
	}

	private double timeThresholdIn(TimeUnit timeUnit)
	{
		return UnitConverterFactory.createConverter(settings.getTimeUnit(), timeUnit, exposureTime)
				.convert(settings.getTimeThreshold());
	}

	/**
	 * Runs the tracing algorithm using distances and time thresholds between min and max with the configured number
	 * of steps. Steps are spaced using a logarithmic scale.
	 * <p>
	 * Returns a list of [distance,time,N traces]
	 *
	 * @param peakResults
	 *            the peak results
	 * @param minDistanceThreshold
	 *            the min distance threshold
	 * @param maxDistanceThreshold
	 *            the max distance threshold
	 * @param minTimeThreshold
	 *            the min time threshold
	 * @param maxTimeThreshold
	 *            the max time threshold
	 * @param optimiserSteps
	 *            the optimiser steps
	 * @return a list of [distance,time,N traces,blinking rate]
	 */
	public List<double[]> runTracing(MemoryPeakResults peakResults, double minDistanceThreshold,
			double maxDistanceThreshold, int minTimeThreshold, int maxTimeThreshold, int optimiserSteps)
	{
		return runTracing(new TraceManager(peakResults), minDistanceThreshold, maxDistanceThreshold, minTimeThreshold,
				maxTimeThreshold, optimiserSteps);
	}

	/**
	 * Runs the tracing algorithm using distances and time thresholds between min and max with the configured number
	 * of steps. Steps are spaced using a logarithmic scale.
	 * <p>
	 * Returns a list of [distance,time,N traces]
	 *
	 * @param manager
	 *            the manager
	 * @param minDistanceThreshold
	 *            the min distance threshold
	 * @param maxDistanceThreshold
	 *            the max distance threshold
	 * @param minTimeThreshold
	 *            the min time threshold
	 * @param maxTimeThreshold
	 *            the max time threshold
	 * @param optimiserSteps
	 *            the optimiser steps
	 * @return a list of [distance,time,N traces,blinking rate]
	 */
	public List<double[]> runTracing(TraceManager manager, double minDistanceThreshold, double maxDistanceThreshold,
			int minTimeThreshold, int maxTimeThreshold, int optimiserSteps)
	{
		dThresholds = getIntervals(minDistanceThreshold, maxDistanceThreshold, optimiserSteps);
		tThresholds = convert(getIntervals(minTimeThreshold, maxTimeThreshold, optimiserSteps));

		final int total = dThresholds.length * tThresholds.length;
		final ArrayList<double[]> results = new ArrayList<>(total);

		IJ.showStatus("Optimising tracing (" + total + " steps) ...");

		if (debugMode)
			IJ.log("Optimising tracing ...");

		int step = 0;
		for (final double d : dThresholds)
			for (final int t : tThresholds)
			{
				IJ.showProgress(step++, total);
				final int n = manager.traceMolecules(d, t);
				results.add(new double[] { d, t, n, getBlinkingRate(manager.getTraces()) });
				if (debugMode)
					summarise(manager.getTraces(), manager.getTotalFiltered(), d, t);
			}

		if (debugMode)
			IJ.log("-=-=-=-");

		IJ.showStatus("");
		IJ.showProgress(1.0);

		return results;
	}

	private static double getBlinkingRate(Trace[] traces)
	{
		final SummaryStatistics stats = new SummaryStatistics();
		for (final Trace trace : traces)
			stats.addValue(trace.getNBlinks());
		final double blinkingRate = stats.getMean();
		return blinkingRate;
	}

	private static double[] getIntervals(double min, double max, int optimiserSteps)
	{
		if (max < min)
		{
			final double tmp = max;
			max = min;
			min = tmp;
		}
		final double range = max - min;
		if (range == 0)
			return new double[] { min };

		final double[] values = new double[optimiserSteps + ((min != 0) ? 1 : 0)];
		int j = 0;
		if (min != 0)
			values[j++] = min;

		// Build a set of steps from min to max

		// Calculate a factor so that:
		//   f^steps = range + 1
		// => factor^n is in bounds [1:range+1] when n <= steps
		final double f = Math.pow(range + 1, 1.0 / optimiserSteps);

		// Set the first increment, i.e. f^1
		double x = f;

		for (int i = 0; i < optimiserSteps; i++)
		{
			// Set the value starting from min.
			// This is equivalent to: values[i] = min + Math.pow(f, i+1) - 1
			// Note that the bounds is [1:range+1] and so 1 is subtracted
			values[j++] = min + x - 1;
			x *= f;
		}

		return values;
	}

	private static int[] convert(double[] intervals)
	{
		final TIntHashSet set = new TIntHashSet(intervals.length);
		for (final double d : intervals)
			set.add((int) Math.round(d));

		set.remove(0); // Do not allow zero

		final int[] values = set.toArray();
		Arrays.sort(values);
		return values;
	}

	/**
	 * Find the contour that intersects zero on the fractional difference plot.
	 * Find the point on the contour nearest the origin.
	 *
	 * @param results
	 *            the results
	 * @return true, if successful
	 */
	private boolean findOptimalParameters(List<double[]> results)
	{
		// This method only works if there are many results and if the results
		// cover enough of the search space to go from above zero (i.e. not enough traces)
		// to below zero (i.e. too many traces)

		final int maxx = tThresholds.length;
		final int maxy = dThresholds.length;

		// --------
		// Find zero crossings using linear interpolation
		zeroCrossingPoints = new ArrayList<>();
		// --------

		// Pass across all time points
		boolean noZeroCrossingAtT0 = false;
		boolean noZeroCrossingAtTN = false;
		for (int x = 0; x < maxx; x++)
		{
			// Find zero crossings on distance points
			final double[] data = new double[maxy];
			for (int y = 0; y < maxy; y++)
			{
				final int i = y * maxx + x;
				final double[] result = results.get(i);
				data[y] = result[2];
			}
			final double zeroCrossing = findZeroCrossing(data, dThresholds);
			if (zeroCrossing > 0)
				zeroCrossingPoints.add(new double[] { tThresholds[x], zeroCrossing });
			else if (x == 0)
				noZeroCrossingAtT0 = true;
			else if (x == maxx - 1)
				noZeroCrossingAtTN = true;
		}

		// If there were not enough zero crossings then the ranges are wrong
		if (zeroCrossingPoints.size() < 3)
		{
			IJ.log(String.format("Very few zero crossings (%d). Increase the optimisation space",
					zeroCrossingPoints.size()));
			return false;
		}

		// --------
		// Use relative distances to find the zero crossing with the smallest distance from origin
		// and set this as the optimal parameters
		// --------
		double minD = Double.MAX_VALUE;
		final double maxTimeThresholdInFrames = settings.getMaxTimeThreshold();
		// The optimiser works using frames so convert back to the correct units
		final TypeConverter<TimeUnit> convert = UnitConverterFactory.createConverter(TimeUnit.FRAME,
				settings.getTimeUnit(), exposureTime);

		for (final double[] point : zeroCrossingPoints)
		{
			final double dx = point[0] / maxTimeThresholdInFrames;
			final double dy = point[1] / settings.getMaxDistanceThreshold();
			final double d = dx * dx + dy * dy;
			if (d < minD)
			{
				minD = d;
				settings.setDistanceThreshold(point[1]);
				settings.setTimeThreshold(convert.convert(point[0]));
			}
		}

		// --------
		// Add more points to make the plotted line look better when showing the plot.
		// --------

		// Pass across all distance points
		boolean noZeroCrossingAtD0 = false;
		boolean noZeroCrossingAtDN = false;
		final double[] tThresholdsD = toDouble(tThresholds);
		for (int y = 0; y < maxy; y++)
		{
			// Find zero crossings on time points
			final double[] data = new double[maxx];
			for (int x = 0; x < maxx; x++)
			{
				final int i = y * maxx + x;
				final double[] result = results.get(i);
				data[x] = result[2];
			}
			final double zeroCrossing = findZeroCrossing(data, tThresholdsD);
			if (zeroCrossing > 0)
				zeroCrossingPoints.add(new double[] { zeroCrossing, dThresholds[y] });
			else if (y == 0)
				noZeroCrossingAtD0 = true;
			else if (y == maxy - 1)
				noZeroCrossingAtDN = true;
		}

		sortPoints();

		// --------
		// Output a message suggesting if the limits should be updated.
		// --------
		final StringBuilder sb = new StringBuilder();
		boolean reduceTime = false;
		boolean reduceDistance = false;
		if (noZeroCrossingAtDN && settings.getMinTimeThreshold() > 1)
		{
			sb.append(" * No zero crossing at max distance\n");
			reduceTime = true;
		}
		if (noZeroCrossingAtTN && settings.getMinDistanceThreshold() > 0)
		{
			sb.append(" * No zero crossing at max time\n");
			reduceDistance = true;
		}
		if (!noZeroCrossingAtD0 && settings.getMinDistanceThreshold() > 0)
		{
			sb.append(" * Zero crossing at min distance\n");
			reduceDistance = true;
		}
		if (!noZeroCrossingAtT0 && settings.getMinTimeThreshold() > 1)
		{
			sb.append(" * Zero crossing at min time\n");
			reduceTime = true;
		}
		if (reduceTime)
			sb.append(" => Reduce the min time threshold\n");
		if (reduceDistance)
			sb.append(" => Reduce the min distance threshold\n");
		if (sb.length() > 0)
		{
			sb.insert(0, "\nWarning:\n");
			sb.append("\n");
			IJ.log(sb.toString());
		}

		// TODO - Fit a function to the zero crossing points. I am not sure what function
		// is suitable for the asymptotic curve (e.g. 1/x == x^-1), perhaps:
		//   f(x) = a + (bx+c)^n
		// where
		//   n < 0
		//   a = Distance asymptote (equivalent to the distance resolution?)
		//   b = Scaling factor
		//   c = Time asymptote

		//interpolateZeroCrossingPoints();

		return true;
	}

	private static double findZeroCrossing(double[] data, double[] axis)
	{
		if (data[0] < 0)
			return -1;
		for (int i = 1; i < data.length; i++)
			if (data[i] < 0)
			{
				final double fraction = data[i - 1] / (data[i - 1] - data[i]);
				return fraction * axis[i] + (1 - fraction) * axis[i - 1];
			}

		return -1;
	}

	private void sortPoints()
	{
		// Sort by x coord, then y
		Collections.sort(zeroCrossingPoints, new Comparator<double[]>()
		{
			@Override
			public int compare(double[] o1, double[] o2)
			{
				if (o1[0] < o2[0])
					return -1;
				if (o1[0] > o2[0])
					return 1;
				if (o1[1] < o2[1])
					return -1;
				if (o1[1] > o2[1])
					return 1;
				return 0;
			}
		});
	}

	@SuppressWarnings("unused")
	private void interpolateZeroCrossingPoints()
	{
		final double[] x = new double[zeroCrossingPoints.size()];
		final double[] y = new double[zeroCrossingPoints.size()];
		for (int i = 0; i < x.length; i++)
		{
			final double[] point = zeroCrossingPoints.get(i);
			x[i] = point[0];
			y[i] = point[1];
		}
		final PolynomialSplineFunction fx = new SplineInterpolator().interpolate(x, y);
		double minX = x[0];
		final double maxX = x[x.length - 1];
		final double xinc = (maxX - minX) / 50;
		for (minX = minX + xinc; minX < maxX; minX += xinc)
			zeroCrossingPoints.add(new double[] { minX, fx.value(minX) });
		sortPoints();
	}

	/**
	 * Build an image using the values within the results to set X,Y and value.
	 *
	 * @param results
	 *            the results
	 */
	private void createPlotResults(List<double[]> results)
	{
		final int w = 400, h = 400;
		switch (OptimiserPlot.get(settings.getOptimiserPlot()))
		{
			case NONE:
				return;
			case BILINEAR:
				fp = createBilinearPlot(results, w, h);
				break;
			default:
				fp = createNNPlot(results, w, h);
		}

		// Create a calibration to map the pixel position back to distance/time
		cal = new Calibration();
		final double xRange = getRange(settings.getMaxTimeThreshold(), settings.getMinTimeThreshold(), origX, w);
		final double yRange = getRange(settings.getMaxDistanceThreshold(), settings.getMinDistanceThreshold(), origY,
				h);
		cal.pixelWidth = xRange / w;
		cal.pixelHeight = yRange / h;
		cal.xOrigin = origX - settings.getMinTimeThreshold() / cal.pixelWidth;
		cal.yOrigin = origY - settings.getMinDistanceThreshold() / cal.pixelHeight;
		cal.setXUnit("frame");
		cal.setYUnit("pixel");

		showPlot();
	}

	/**
	 * Shows the plot
	 */
	private void showPlot()
	{
		if (OptimiserPlot.get(settings.getOptimiserPlot()) == OptimiserPlot.NONE)
			return;

		// Display the image
		final String title = TITLE + ": | N - N_actual | / N_actual";
		ImagePlus imp = WindowManager.getImage(title);
		if (imp != null)
		{
			fp.setColorModel(imp.getProcessor().getColorModel());
			imp.setProcessor(fp);
		}
		else
		{
			imp = new ImagePlus(title, fp);
			imp.show();
			WindowManager.setTempCurrentImage(imp);
			final LutLoader lut = new LutLoader();
			lut.run("fire");
			WindowManager.setTempCurrentImage(null);
		}

		imp.setCalibration(cal);
		addZeroCrossingPoints(imp);
		imp.updateAndDraw();
	}

	private void addZeroCrossingPoints(ImagePlus imp)
	{
		PolygonRoi roi = null;
		imp.setRoi(roi);
		if (zeroCrossingPoints == null || zeroCrossingPoints.isEmpty())
			return;
		final Calibration cal = imp.getCalibration();
		final int nPoints = zeroCrossingPoints.size();
		final float[] xPoints = new float[nPoints];
		final float[] yPoints = new float[nPoints];
		for (int i = 0; i < nPoints; i++)
		{
			final double[] point = zeroCrossingPoints.get(i);
			// Convert to pixel coordinates.
			xPoints[i] = (float) (cal.xOrigin + (point[0] / cal.pixelWidth));
			yPoints[i] = (float) (cal.yOrigin + (point[1] / cal.pixelHeight));
		}
		roi = new PolygonRoi(xPoints, yPoints, nPoints, Roi.POLYLINE);
		imp.setRoi(roi);
	}

	private FloatProcessor createNNPlot(List<double[]> results, int w, int h)
	{
		final FloatProcessor fp = new FloatProcessor(w, h);

		// Create lookup table that map the tested threshold values to a position in the image
		final int[] xLookup = createLookup(tThresholds, settings.getMinTimeThreshold(), w);
		final int[] yLookup = createLookup(dThresholds, settings.getMinDistanceThreshold(), h);
		origX = (settings.getMinTimeThreshold() != 0) ? xLookup[1] : 0;
		origY = (settings.getMinDistanceThreshold() != 0) ? yLookup[1] : 0;

		final int gridWidth = tThresholds.length;
		final int gridHeight = dThresholds.length;
		for (int y = 0, i = 0; y < gridHeight; y++)
			for (int x = 0; x < gridWidth; x++, i++)
			{
				final int x1 = xLookup[x];
				final int x2 = xLookup[x + 1];
				final int y1 = yLookup[y];
				final int y2 = yLookup[y + 1];
				final double[] result = results.get(i);
				fp.setValue(Math.abs(result[2]));
				fp.setRoi(x1, y1, x2 - x1, y2 - y1);
				fp.fill();
			}
		return fp;
	}

	private static int[] createLookup(int[] values, int min, int scale)
	{
		final double[] newValues = toDouble(values);
		return createLookup(newValues, min, scale);
	}

	private static double[] toDouble(int[] values)
	{
		final double[] newValues = new double[values.length];
		for (int i = 0; i < values.length; i++)
			newValues[i] = values[i];
		return newValues;
	}

	private static int[] createLookup(double[] values, double min, int scale)
	{
		// To allow the lowest result to be plotted, add space at the edge
		// equal to the next interval
		if (min != 0 && values.length > 1)
			min -= values[1] - values[0];

		final int[] lookup = new int[values.length + 1];
		final double range = values[values.length - 1] - min;
		final double scaleFactor = scale / range;
		for (int i = 1; i < values.length; i++)
			lookup[i] = (int) Math.round(scaleFactor * (values[i - 1] - min));
		lookup[values.length] = scale;
		return lookup;
	}

	private FloatProcessor createBilinearPlot(List<double[]> results, int w, int h)
	{
		final FloatProcessor fp = new FloatProcessor(w, h);

		// Create lookup table that map the tested threshold values to a position in the image
		final int[] xLookup = createLookup(tThresholds, settings.getMinTimeThreshold(), w);
		final int[] yLookup = createLookup(dThresholds, settings.getMinDistanceThreshold(), h);
		origX = (settings.getMinTimeThreshold() != 0) ? xLookup[1] : 0;
		origY = (settings.getMinDistanceThreshold() != 0) ? yLookup[1] : 0;

		final int gridWidth = tThresholds.length;
		final int gridHeight = dThresholds.length;
		for (int y = 0, prevY = 0; y < gridHeight; y++)
		{
			for (int x = 0, prevX = 0; x < gridWidth; x++)
			{
				// Get the 4 flanking values
				final double x1y1 = results.get(prevY * gridWidth + prevX)[2];
				final double x1y2 = results.get(y * gridWidth + prevX)[2];
				final double x2y1 = results.get(prevY * gridWidth + x)[2];
				final double x2y2 = results.get(y * gridWidth + x)[2];

				// Pixel range
				final int x1 = xLookup[x];
				final int x2 = xLookup[x + 1];
				final int y1 = yLookup[y];
				final int y2 = yLookup[y + 1];

				final double xRange = x2 - x1;
				final double yRange = y2 - y1;

				for (int yy = y1; yy < y2; yy++)
				{
					final double yFraction = (yy - y1) / yRange;
					for (int xx = x1; xx < x2; xx++)
					{
						// Interpolate
						final double xFraction = (xx - x1) / xRange;
						final double v1 = x1y1 * (1 - xFraction) + x2y1 * xFraction;
						final double v2 = x1y2 * (1 - xFraction) + x2y2 * xFraction;
						final double value = v1 * (1 - yFraction) + v2 * yFraction;
						fp.setf(xx, yy, (float) value);
					}
				}

				prevX = x;
			}
			prevY = y;
		}

		// Convert to absolute for easier visualisation
		final float[] data = (float[]) fp.getPixels();
		for (int i = 0; i < data.length; i++)
			data[i] = Math.abs(data[i]);

		return fp;
	}

	private static double getRange(double max, double min, int orig, int w)
	{
		double r = max - min;
		if (r <= 0)
			r = 1;
		return r * w / (w - orig);
	}

	@SuppressWarnings("unused")
	private void fitTraces(MemoryPeakResults results, Trace[] traces)
	{
		//		// Check if the original image is open and the fit configuration can be extracted
		//		ImageSource source = results.getSource();
		//		if (source == null)
		//			return;
		//		if (!source.open())
		//			return;
		//		FitEngineConfiguration config = (FitEngineConfiguration) XmlUtils.fromXML(results.getConfiguration());
		//		if (config == null)
		//			return;
		//
		//		// Show a dialog asking if the traces should be refit
		//		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		//		gd.addMessage("Do you want to fit the traces as a single peak using a combined image?");
		//
		//		//gd.addCheckbox("Fit_closest_to_centroid", !fitOnlyCentroid);
		//		//gd.addSlider("Distance_threshold", 0.01, 3, distanceThreshold);
		//		gd.addSlider("Expansion_factor", 1, 4.5, expansionFactor);
		//
		//		// Allow fitting settings to be adjusted
		//		FitConfiguration fitConfig = config.getFitConfiguration();
		//		gd.addMessage("--- Gaussian fitting ---");
		//		gd.addChoice("Spot_filter_type", SettingsManager.dataFilterTypeNames,
		//				SettingsManager.dataFilterTypeNames[config.getDataFilterType().ordinal()]);
		//		gd.addChoice("Spot_filter", SettingsManager.dataFilterNames,
		//				SettingsManager.dataFilterNames[config.getDataFilter(0).ordinal()]);
		//		gd.addSlider("Smoothing", 0, 2.5, config.getSmooth(0));
		//		gd.addSlider("Search_width", 0.5, 2.5, config.getSearch());
		//		gd.addSlider("Border", 0.5, 2.5, config.getBorder());
		//		gd.addSlider("Fitting_width", 2, 4.5, config.getFitting());
		//
		//		gd.addChoice("Fit_solver", SettingsManager.fitSolverNames,
		//				SettingsManager.fitSolverNames[fitConfig.getFitSolver().ordinal()]);
		//		gd.addChoice("Fit_function", SettingsManager.fitFunctionNames,
		//				SettingsManager.fitFunctionNames[fitConfig.getFitFunction().ordinal()]);
		//
		//		gd.addChoice("Fit_criteria", SettingsManager.fitCriteriaNames,
		//				SettingsManager.fitCriteriaNames[fitConfig.getFitCriteria().ordinal()]);
		//		gd.addNumericField("Significant_digits", fitConfig.getSignificantDigits(), 0);
		//		gd.addNumericField("Coord_delta", fitConfig.getDelta(), 4);
		//		gd.addNumericField("Lambda", fitConfig.getLambda(), 4);
		//		gd.addNumericField("Max_iterations", fitConfig.getMaxIterations(), 0);
		//		gd.addNumericField("Fail_limit", config.getFailuresLimit(), 0);
		//		gd.addCheckbox("Include_neighbours", config.isIncludeNeighbours());
		//		gd.addSlider("Neighbour_height", 0.01, 1, config.getNeighbourHeightThreshold());
		//		gd.addSlider("Residuals_threshold", 0.01, 1, config.getResidualsThreshold());
		//
		//		//gd.addSlider("Duplicate_distance", 0, 1.5, fitConfig.getDuplicateDistance());
		//
		//		gd.addMessage("--- Peak filtering ---\nDiscard fits that shift; are too low; or expand/contract");
		//
		//		gd.addCheckbox("Smart_filter", fitConfig.isSmartFilter());
		//		gd.addCheckbox("Disable_simple_filter", fitConfig.isDisableSimpleFilter());
		//		gd.addSlider("Shift_factor", 0.01, 2, fitConfig.getCoordinateShiftFactor());
		//		gd.addNumericField("Signal_strength", fitConfig.getSignalStrength(), 2);
		//		gd.addNumericField("Min_photons", fitConfig.getMinPhotons(), 0);
		//		gd.addSlider("Min_width_factor", 0, 0.99, fitConfig.getMinWidthFactor());
		//		gd.addSlider("Width_factor", 1, 4.5, fitConfig.getWidthFactor());
		//		gd.addNumericField("Precision", fitConfig.getPrecisionThreshold(), 2);
		//
		//		gd.addCheckbox("Debug_failures", debugFailures);
		//		gd.showDialog();
		//		if (!gd.wasOKed())
		//		{
		//			source.close();
		//			return;
		//		}
		//
		//		// Get parameters for the fit
		//		//fitOnlyCentroid = !gd.getNextBoolean();
		//		//distanceThreshold = (float) gd.getNextNumber();
		//		expansionFactor = (float) gd.getNextNumber();
		//
		//		config.setDataFilterType(gd.getNextChoiceIndex());
		//		config.setDataFilter(gd.getNextChoiceIndex(), Math.abs(gd.getNextNumber()), 0);
		//		config.setSearch(gd.getNextNumber());
		//		config.setBorder(gd.getNextNumber());
		//		config.setFitting(gd.getNextNumber());
		//		fitConfig.setFitSolver(gd.getNextChoiceIndex());
		//		fitConfig.setFitFunction(gd.getNextChoiceIndex());
		//		fitConfig.setFitCriteria(gd.getNextChoiceIndex());
		//
		//		fitConfig.setSignificantDigits((int) gd.getNextNumber());
		//		fitConfig.setDelta(gd.getNextNumber());
		//		fitConfig.setLambda(gd.getNextNumber());
		//		fitConfig.setMaxIterations((int) gd.getNextNumber());
		//		config.setFailuresLimit((int) gd.getNextNumber());
		//		config.setIncludeNeighbours(gd.getNextBoolean());
		//		config.setNeighbourHeightThreshold(gd.getNextNumber());
		//		config.setResidualsThreshold(gd.getNextNumber());
		//
		//		fitConfig.setSmartFilter(gd.getNextBoolean());
		//		fitConfig.setDisableSimpleFilter(gd.getNextBoolean());
		//		fitConfig.setCoordinateShiftFactor(gd.getNextNumber());
		//		fitConfig.setSignalStrength(gd.getNextNumber());
		//		fitConfig.setMinPhotons(gd.getNextNumber());
		//		fitConfig.setMinWidthFactor(gd.getNextNumber());
		//		fitConfig.setWidthFactor(gd.getNextNumber());
		//		fitConfig.setPrecisionThreshold(gd.getNextNumber());
		//
		//		// Check arguments
		//		try
		//		{
		//			//Parameters.isAboveZero("Distance threshold", distanceThreshold);
		//			Parameters.isAbove("Expansion factor", expansionFactor, 1);
		//			Parameters.isAboveZero("Search_width", config.getSearch());
		//			Parameters.isAboveZero("Fitting_width", config.getFitting());
		//			Parameters.isAboveZero("Significant digits", fitConfig.getSignificantDigits());
		//			Parameters.isAboveZero("Delta", fitConfig.getDelta());
		//			Parameters.isAboveZero("Lambda", fitConfig.getLambda());
		//			Parameters.isAboveZero("Max iterations", fitConfig.getMaxIterations());
		//			Parameters.isPositive("Failures limit", config.getFailuresLimit());
		//			Parameters.isPositive("Neighbour height threshold", config.getNeighbourHeightThreshold());
		//			Parameters.isPositive("Residuals threshold", config.getResidualsThreshold());
		//			Parameters.isPositive("Coordinate Shift factor", fitConfig.getCoordinateShiftFactor());
		//			Parameters.isPositive("Signal strength", fitConfig.getSignalStrength());
		//			Parameters.isPositive("Min photons", fitConfig.getMinPhotons());
		//			Parameters.isPositive("Min width factor", fitConfig.getMinWidthFactor());
		//			Parameters.isPositive("Width factor", fitConfig.getWidthFactor());
		//			Parameters.isPositive("Precision threshold", fitConfig.getPrecisionThreshold());
		//		}
		//		catch (IllegalArgumentException e)
		//		{
		//			IJ.error(TITLE, e.getMessage());
		//			source.close();
		//			return;
		//		}
		//
		//		debugFailures = gd.getNextBoolean();
		//
		//		if (!PeakFit.configureSmartFilter(globalSettings, null))
		//			return;
		//		if (!PeakFit.configureDataFilter(globalSettings, null, false))
		//			return;
		//		if (!PeakFit.configureFitSolver(globalSettings, null, false))
		//			return;
		//
		//		// Adjust settings for a single maxima
		//		config.setIncludeNeighbours(false);
		//		fitConfig.setDuplicateDistance(0);
		//
		//		// Create a fit engine
		//		MemoryPeakResults refitResults = new MemoryPeakResults();
		//		refitResults.copySettings(results);
		//		refitResults.setName(results.getName() + " Trace Fit");
		//		refitResults.setSortAfterEnd(true);
		//		refitResults.begin();
		//		int threadCount = Prefs.getThreads();
		//		PeakResults syncResults = SynchronizedPeakResults.create(refitResults, threadCount);
		//		// No border since we know where the peaks are and we must not miss them due to truncated searching
		//		FitEngine engine = new FitEngine(config, syncResults, threadCount, FitQueue.BLOCKING);
		//
		//		// Either : Only fit the centroid
		//		// or     : Extract a bigger region, allowing all fits to run as normal and then
		//		//          find the correct spot using Euclidian distance.
		//
		//		// Set up the limits
		//		final double stdDev = FastMath.max(fitConfig.getInitialPeakStdDev0(), fitConfig.getInitialPeakStdDev1());
		//		float fitWidth = (float) (stdDev * config.getFitting() * expansionFactor);
		//
		//		IJ.showStatus("Refitting traces ...");
		//
		//		List<JobItem> jobItems = new ArrayList<JobItem>(traces.length);
		//		int singles = 0;
		//		int fitted = 0;
		//		for (int n = 0; n < traces.length; n++)
		//		{
		//			Trace trace = traces[n];
		//
		//			if (n % 32 == 0)
		//				IJ.showProgress(n, traces.length);
		//
		//			// Skip traces with one peak
		//			if (trace.size() == 1)
		//			{
		//				singles++;
		//				// Use the synchronized method to avoid thread clashes with the FitEngine
		//				syncResults.add(trace.getHead());
		//				continue;
		//			}
		//
		//			Rectangle bounds = new Rectangle();
		//			double[] combinedNoise = new double[1];
		//			float[] data = buildCombinedImage(source, trace, fitWidth, bounds, combinedNoise, false);
		//			if (data == null)
		//				continue;
		//
		//			// Fit the combined image
		//			FitParameters params = new FitParameters();
		//			params.noise = (float) combinedNoise[0];
		//			float[] centre = trace.getCentroid();
		//
		//			//if (fitOnlyCentroid)
		//			//{
		//			int newX = (int) Math.round(centre[0]) - bounds.x;
		//			int newY = (int) Math.round(centre[1]) - bounds.y;
		//			params.maxIndices = new int[] { newY * bounds.width + newX };
		//			//}
		//			//else
		//			//{
		//			//	params.filter = new ArrayList<float[]>();
		//			//	params.filter.add(new float[] { centre[0] - bounds.x, centre[1] - bounds.y });
		//			//	params.distanceThreshold = distanceThreshold;
		//			//}
		//
		//			// This is not needed since the bounds are passed using the FitJob
		//			//params.setOffset(new float[] { bounds.x, bounds.y });
		//			int startT = trace.getHead().getFrame();
		//			params.endT = trace.getTail().getFrame();
		//
		//			ParameterisedFitJob job = new ParameterisedFitJob(n, params, startT, data, bounds);
		//			jobItems.add(new JobItem(job, trace, centre));
		//			engine.run(job);
		//			fitted++;
		//		}
		//
		//		engine.end(false);
		//
		//		IJ.showStatus("");
		//		IJ.showProgress(1);
		//
		//		// Check the success ...
		//		FitStatus[] values = FitStatus.values();
		//		int[] statusCount = new int[values.length + 1];
		//		ArrayList<String> names = new ArrayList<String>(Arrays.asList(SettingsManager.getNames((Object[]) values)));
		//		//names.add(String.format("No maxima within %.2f of centroid", distanceThreshold));
		//		int separated = 0;
		//		int success = 0;
		//		final int debugLimit = 3;
		//		for (JobItem jobItem : jobItems)
		//		{
		//			int id = jobItem.getId();
		//			ParameterisedFitJob job = jobItem.job;
		//			Trace trace = jobItem.trace;
		//			int[] indices = job.getIndices();
		//			FitResult fitResult = null;
		//			int status;
		//			if (indices.length < 1)
		//			{
		//				status = values.length;
		//			}
		//			else if (indices.length > 1)
		//			{
		//				//System.out.printf("Multiple fits performed for trace : Job Id = %d\n", id);
		//
		//				// Fits are recorded if (a) they succeeded and were close to the target centroid;
		//				// or (b) if they failed and started close to the target centroid.
		//
		//				// Choose the first OK result. This is all that matters for the success reporting
		//				for (int n = 0; n < indices.length; n++)
		//				{
		//					if (job.getFitResult(n).getStatus() == FitStatus.OK)
		//					{
		//						fitResult = job.getFitResult(n);
		//						break;
		//					}
		//				}
		//				// Otherwise use the closest failure.
		//				if (fitResult == null)
		//				{
		//					final float[] centre = traces[id].getCentroid();
		//					double minD = Double.POSITIVE_INFINITY;
		//					for (int n = 0; n < indices.length; n++)
		//					{
		//						// Since the fit has failed we use the initial parameters.
		//						// Note: This assumes the initial parameters are for a Gaussian 2D function
		//						final double[] params = job.getFitResult(n).getInitialParameters();
		//						final double dx = params[Gaussian2DFunction.X_POSITION] - centre[0];
		//						final double dy = params[Gaussian2DFunction.Y_POSITION] - centre[1];
		//						final double d = dx * dx + dy * dy;
		//						if (minD > d)
		//						{
		//							minD = d;
		//							fitResult = job.getFitResult(n);
		//						}
		//					}
		//				}
		//
		//				status = fitResult.getStatus().ordinal();
		//			}
		//			else
		//			{
		//				fitResult = job.getFitResult(0);
		//				status = fitResult.getStatus().ordinal();
		//			}
		//
		//			// All jobs have only one peak
		//			statusCount[status]++;
		//
		//			// Debug why any fits failed
		//			if (fitResult == null || fitResult.getStatus() != FitStatus.OK)
		//			{
		//				refitResults.addAll(trace.getPoints());
		//				separated += trace.size();
		//
		//				if (debugFailures)
		//				{
		//					FitStatus s = (fitResult == null) ? FitStatus.UNKNOWN : fitResult.getStatus();
		//
		//					// Only display the first n per category to limit the number of images
		//					double[] noise = new double[1];
		//					if (statusCount[status] <= debugLimit)
		//					{
		//						Rectangle bounds = new Rectangle();
		//						buildCombinedImage(source, trace, fitWidth, bounds, noise, true);
		//						float[] centre = trace.getCentroid();
		//						Utils.display(
		//								String.format("Trace %d (n=%d) : x=%f,y=%f", id, trace.size(), centre[0], centre[1]),
		//								slices);
		//
		//						switch (s)
		//						{
		//							case INSUFFICIENT_PRECISION:
		//								float precision = (Float) fitResult.getStatusData();
		//								IJ.log(String.format("Trace %d (n=%d) : %s = %f", id, trace.size(), names.get(status),
		//										precision));
		//								break;
		//							case INSUFFICIENT_SIGNAL:
		//								if (noise[0] == 0)
		//									noise[0] = getCombinedNoise(trace);
		//								float snr = (Float) fitResult.getStatusData();
		//								IJ.log(String.format("Trace %d (n=%d) : %s = %f (noise=%.2f)", id, trace.size(),
		//										names.get(status), snr, noise[0]));
		//								break;
		//							case COORDINATES_MOVED:
		//							case OUTSIDE_FIT_REGION:
		//							case WIDTH_DIVERGED:
		//								float[] shift = (float[]) fitResult.getStatusData();
		//								IJ.log(String.format("Trace %d (n=%d) : %s = %.3f,%.3f", id, trace.size(),
		//										names.get(status), shift[0], shift[1]));
		//								break;
		//							default:
		//								IJ.log(String.format("Trace %d (n=%d) : %s", id, trace.size(), names.get(status)));
		//								break;
		//						}
		//					}
		//				}
		//			}
		//			else
		//			{
		//				success++;
		//
		//				if (debugFailures)
		//				{
		//					// Only display the first n per category to limit the number of images
		//					double[] noise = new double[1];
		//					if (statusCount[status] <= debugLimit)
		//					{
		//						Rectangle bounds = new Rectangle();
		//						buildCombinedImage(source, trace, fitWidth, bounds, noise, true);
		//						float[] centre = trace.getCentroid();
		//						Utils.display(
		//								String.format("Trace %d (n=%d) : x=%f,y=%f", id, trace.size(), centre[0], centre[1]),
		//								slices);
		//					}
		//				}
		//			}
		//		}
		//
		//		IJ.log(String.format("Trace fitting : %d singles : %d / %d fitted : %d separated", singles, success, fitted,
		//				separated));
		//		if (separated > 0)
		//		{
		//			IJ.log("Reasons for fit failure :");
		//			// Start at i=1 to skip FitStatus.OK
		//			for (int i = 1; i < statusCount.length; i++)
		//			{
		//				if (statusCount[i] != 0)
		//					IJ.log("  " + names.get(i) + " = " + statusCount[i]);
		//			}
		//		}
		//
		//		refitResults.end();
		//		MemoryPeakResults.addResults(refitResults);
		//
		//		source.close();
	}

	private ImageStack slices;

	@SuppressWarnings("unused")
	private float[] buildCombinedImage(ImageSource source, Trace trace, float fitWidth, Rectangle bounds,
			double[] combinedNoise, boolean createStack)
	{
		final int w = source.getWidth();
		final int h = source.getHeight();

		// Get the coordinates and the spot bounds
		final float[] centre = trace.getCentroid(CentroidMethod.SIGNAL_WEIGHTED);
		int minX = (int) Math.floor(centre[0] - fitWidth);
		int maxX = (int) Math.ceil(centre[0] + fitWidth);
		int minY = (int) Math.floor(centre[1] - fitWidth);
		int maxY = (int) Math.ceil(centre[1] + fitWidth);

		// Account for crops at the edge of the image
		minX = FastMath.max(0, minX);
		maxX = FastMath.min(w, maxX);
		minY = FastMath.max(0, minY);
		maxY = FastMath.min(h, maxY);

		final int width = maxX - minX;
		final int height = maxY - minY;
		if (width <= 0 || height <= 0)
			// The centre must be outside the image width and height
			return null;
		bounds.x = minX;
		bounds.y = minY;
		bounds.width = width;
		bounds.height = height;

		if (createStack)
			slices = new ImageStack(width, height);

		// Combine the images. Subtract the fitted background to zero the image.
		final float[] data = new float[width * height];
		float sumBackground = 0;
		double noise = 0;
		for (int i = 0; i < trace.size(); i++)
		{
			final PeakResult result = trace.get(i);
			noise += result.getNoise() * result.getNoise();

			final float[] sourceData = source.get(result.getFrame(), bounds);
			final float background = result.getBackground();
			sumBackground += background;
			for (int j = 0; j < data.length; j++)
				data[j] += sourceData[j] - background;
			if (createStack)
				slices.addSlice(new FloatProcessor(width, height, sourceData, null));
		}
		if (createStack)
		{
			// Add a final image that is the average of the individual slices. This allows
			// it to be visualised in the same intensity scale.
			final float[] data2 = Arrays.copyOf(data, data.length);
			final int size = slices.getSize();
			sumBackground /= size;
			for (int i = 0; i < data2.length; i++)
				data2[i] = sumBackground + data2[i] / size;
			slices.addSlice(new FloatProcessor(width, height, data2, null));
		}

		// Combined noise is the sqrt of the sum-of-squares
		combinedNoise[0] = Math.sqrt(noise);

		return data;
	}

	@SuppressWarnings("unused")
	private static double getCombinedNoise(Trace trace)
	{
		double noise = 0;
		for (int i = 0; i < trace.size(); i++)
		{
			final PeakResult result = trace.get(i);
			noise += result.getNoise() * result.getNoise();
		}
		// Combined noise is the sqrt of the sum-of-squares
		return Math.sqrt(noise);
	}

	@SuppressWarnings("unused")
	private class JobItem
	{
		ParameterisedFitJob job;
		Trace trace;

		public JobItem(ParameterisedFitJob job, Trace trace, float[] centre)
		{
			this.job = job;
			this.trace = trace;
		}

		public int getId()
		{
			return job.getId();
		}
	}
}
