package gdsc.smlm.ij.plugins;

import java.awt.Rectangle;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.util.FastMath;

import gdsc.core.clustering.Cluster;
import gdsc.core.clustering.ClusterPoint;
import gdsc.core.clustering.ClusteringAlgorithm;
import gdsc.core.clustering.ClusteringEngine;
import gdsc.core.data.utils.TypeConverter;
import gdsc.core.ij.IJTrackProgress;
import gdsc.core.ij.Utils;
import gdsc.core.utils.Statistics;
import gdsc.core.utils.StoredDataStatistics;
import gdsc.smlm.data.config.UnitConfig.TimeUnit;
import gdsc.smlm.data.config.UnitConverterFactory;
import gdsc.smlm.engine.ParameterisedFitJob;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.settings.ClusteringSettings;
import gdsc.smlm.ij.settings.ClusteringSettingsHelper;
import gdsc.smlm.ij.settings.ClusteringSettingsHelper.OptimiserPlot;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.results.Cluster.CentroidMethod;
import gdsc.smlm.results.Counter;
import gdsc.smlm.results.ImageSource;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.TextFilePeakResults;
import gdsc.smlm.results.Trace;
import gdsc.smlm.results.TraceManager;
import gdsc.smlm.results.procedures.PeakResultProcedure;
import gdsc.smlm.results.procedures.PrecisionResultProcedure;
import gnu.trove.set.hash.TIntHashSet;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.ExtendedGenericDialog;
import ij.gui.PolygonRoi;
import ij.measure.Calibration;
import ij.plugin.LutLoader;
import ij.plugin.PlugIn;
import ij.plugin.WindowOrganiser;
import ij.process.FloatProcessor;
import ij.text.TextWindow;

/**
 * Run a tracing algorithm on the peak results to trace molecules across the frames.
 */
@SuppressWarnings("unused")
public class TraceMolecules implements PlugIn
{
	private String TITLE = "Trace or Cluster Molecules";
	private String outputName;
	private static double MIN_BLINKING_RATE = 1; // Should never be <= 0
	private static String inputOption = "";
	private static boolean inputDebugMode = true;
	private static boolean inputOptimiseBlinkingRate = false;

	//private static boolean fitOnlyCentroid = false;
	//private static float distanceThreshold = 1;
	private static float expansionFactor = 2;
	private static boolean debugFailures = false;

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

	private GlobalSettings globalSettings;
	private ClusteringSettings settings;
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

			ClusteringEngine engine = new ClusteringEngine(Prefs.getThreads(),
					ClusteringSettingsHelper.getClusteringAlgorithm(settings.getClusteringAlgorithm()),
					new IJTrackProgress());

			if (settings.splitPulses)
			{
				engine.setPulseInterval(settings.pulseInterval);
				limitTimeThreshold(settings, settings.pulseInterval);
			}

			ArrayList<Cluster> clusters = engine.findClusters(convertToClusterPoints(),
					settings.distanceThreshold / results.getCalibrationReader().getNmPerPixel(),
					timeInFrames(settings));

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

			TraceManager manager = new TraceManager(results);
			manager.setTraceMode(ClusteringSettingsHelper.getTraceMode(settings.getTraceMode()));
			manager.setActivationFrameInterval(settings.pulseInterval);
			manager.setActivationFrameWindow(settings.pulseWindow);
			manager.setDistanceExclusion(settings.distanceExclusion / results.getCalibrationReader().getNmPerPixel());

			if (settings.optimise)
			{
				// Optimise before configuring for a pulse interval
				runOptimiser(manager);
			}

			if (settings.splitPulses)
			{
				manager.setPulseInterval(settings.pulseInterval);
				limitTimeThreshold(settings, settings.pulseInterval);
			}

			manager.setTracker(new IJTrackProgress());
			manager.traceMolecules(settings.distanceThreshold / results.getCalibrationReader().getNmPerPixel(),
					timeInFrames(settings));
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
		Trace[] multiTraces = getTraces(traces);
		saveResults(results, multiTraces, outputName + " Multi");

		// Save centroids
		outputName += " Centroids";
		MemoryPeakResults tracedResults = saveCentroidResults(results, traces, outputName);

		// Save traces separately
		saveCentroidResults(results, multiTraces, outputName + " Multi");

		// Sort traces by time to assist the results source in extracting frames sequentially.
		// Do this before saving to assist in debugging using the saved traces file.
		sortByTime(traces);

		if (settings.saveTraces)
			saveTraces(traces);

		summarise(traces, totalFiltered, settings.distanceThreshold, timeInSeconds(settings));

		IJ.showStatus(String.format("%d localisations => %d traces (%d filtered)", results.size(), tracedResults.size(),
				totalFiltered));

		//// Provide option to refit the traces as single peaks and save to memory
		//if (settings.refitOption)
		//	fitTraces(results, traces);
	}

	/**
	 * Limit the time threshold to the pulse interval duration.
	 *
	 * @param settings
	 *            the settings
	 * @param pulseInterval
	 *            the pulse interval
	 */
	private void limitTimeThreshold(ClusteringSettings settings, int pulseInterval)
	{
		// The pulse interval is in frames.
		// Convert the interval to the correct units. 

		double limit;
		if (settings.getTimeUnit() == TimeUnit.FRAME)
		{
			limit = settings.pulseInterval;
		}
		else
		{
			TypeConverter<TimeUnit> convert = UnitConverterFactory.createConverter(TimeUnit.FRAME,
					settings.getTimeUnit(), exposureTime);
			limit = convert.convert(pulseInterval);
		}

		if (settings.getTimeThreshold() > limit)
		{
			settings.setTimeThreshold(limit);
		}
	}

	private List<ClusterPoint> convertToClusterPoints()
	{
		return convertToClusterPoints(results);
	}

	/**
	 * Convert a list of peak results into points for the clustering engine
	 * 
	 * @param results
	 * @return
	 */
	public static List<ClusterPoint> convertToClusterPoints(MemoryPeakResults results)
	{
		final ArrayList<ClusterPoint> points = new ArrayList<ClusterPoint>(results.size());
		final Counter counter = new Counter();
		results.forEach(new PeakResultProcedure()
		{
			public void execute(PeakResult p)
			{
				points.add(ClusterPoint.newTimeClusterPoint(counter.getAndIncrement(), p.getXPosition(),
						p.getYPosition(), p.getSignal(), p.getFrame(), p.getEndFrame()));
			}
		});
		return points;
	}

	private Trace[] convertToTraces(ArrayList<Cluster> clusters)
	{
		return convertToTraces(results, clusters);
	}

	/**
	 * Convert the clusters from the clustering engine into traces composed of the original list of peak results
	 * 
	 * @param results
	 * @param clusters
	 * @return
	 */
	public static Trace[] convertToTraces(MemoryPeakResults results, ArrayList<Cluster> clusters)
	{
		Trace[] traces = new Trace[clusters.size()];
		int i = 0;
		for (Cluster cluster : clusters)
		{
			Trace trace = new Trace();
			trace.setId(i + 1);
			for (ClusterPoint point = cluster.head; point != null; point = point.next)
			{
				// The point Id was the position in the original results array
				trace.add(results.get(point.id));
			}
			traces[i++] = trace;
		}
		return traces;
	}

	/**
	 * Sort traces by time
	 * 
	 * @param traces
	 */
	static void sortByTime(Trace[] traces)
	{
		for (Trace t : traces)
			t.sort();
		Arrays.sort(traces, new Comparator<Trace>()
		{
			public int compare(Trace o1, Trace o2)
			{
				return o1.getHead().getFrame() - o2.getHead().getFrame();
			}
		});
	}

	static MemoryPeakResults saveResults(MemoryPeakResults sourceResults, Trace[] traces, String name)
	{
		MemoryPeakResults tracedResults = TraceManager.convertToPeakResults(sourceResults, traces);
		tracedResults.setName(sourceResults.getName() + " " + name);
		MemoryPeakResults.addResults(tracedResults);
		return tracedResults;
	}

	static MemoryPeakResults saveCentroidResults(MemoryPeakResults sourceResults, Trace[] traces, String name)
	{
		MemoryPeakResults tracedResults = TraceManager.convertToCentroidPeakResults(sourceResults, traces);
		tracedResults.setName(sourceResults.getName() + " " + name);
		MemoryPeakResults.addResults(tracedResults);
		return tracedResults;
	}

	private Trace[] getSingles(Trace[] traces)
	{
		ArrayList<Trace> result = new ArrayList<Trace>();
		for (Trace t : traces)
			if (t.size() == 1)
				result.add(t);
		return result.toArray(new Trace[result.size()]);
	}

	private Trace[] getTraces(Trace[] traces)
	{
		ArrayList<Trace> result = new ArrayList<Trace>();
		for (Trace t : traces)
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
	 * @param traces
	 * @param comment
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
			String regex = "\\.[0-9]+\\.";
			if (filename.matches(regex))
				filename.replaceAll(regex, "." + (id) + ".");
			else
				Utils.replaceExtension(filename, id + ".xls");
		}
		filename = Utils.getFilename(title, filename);
		if (filename != null)
		{
			filename = Utils.replaceExtension(filename, "xls");

			boolean showDeviations = sourceResults.hasDeviations();
			// Assume that are results are from a single frame but store the trace ID
			TextFilePeakResults traceResults = new TextFilePeakResults(filename, showDeviations, false, true);
			traceResults.copySettings(sourceResults);
			traceResults.begin();
			if (!traceResults.isActive())
			{
				IJ.error("Failed to write to file: " + filename);
			}
			else
			{
				traceResults.addComment(comment);
				for (Trace trace : traces)
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
				settings.distanceThreshold, timeInSeconds(settings), timeInFrames(settings));
	}

	private void summarise(Trace[] traces, int filtered, double dThreshold, double tThreshold)
	{
		IJ.showStatus("Calculating summary ...");

		// Create summary table
		createSummaryTable();

		Statistics[] stats = new Statistics[NAMES.length];
		for (int i = 0; i < stats.length; i++)
		{
			stats[i] = (settings.showHistograms || settings.saveTraceData) ? new StoredDataStatistics()
					: new Statistics();
		}
		int singles = 0;
		for (Trace trace : traces)
		{
			int nBlinks = trace.getNBlinks() - 1;
			stats[BLINKS].add(nBlinks);
			int[] onTimes = trace.getOnTimes();
			int[] offTimes = trace.getOffTimes();
			double tOn = 0;
			for (int t : onTimes)
			{
				stats[T_ON].add(t * exposureTime);
				tOn += t * exposureTime;
			}
			stats[TOTAL_T_ON].add(tOn);
			if (offTimes != null)
			{
				double tOff = 0;
				for (int t : offTimes)
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
		StringBuilder sb = new StringBuilder();
		sb.append(results.getName()).append('\t');
		sb.append(outputName.equals("Cluster")
				? ClusteringSettingsHelper.getClusteringAlgorithm(settings.getClusteringAlgorithm())
				: ClusteringSettingsHelper.getTraceMode(settings.getTraceMode())).append('\t');
		sb.append(Utils.rounded(exposureTime * 1000, 3)).append('\t');
		sb.append(Utils.rounded(dThreshold, 3)).append('\t');
		sb.append(Utils.rounded(tThreshold, 3));
		if (settings.splitPulses)
			sb.append(" *");
		sb.append('\t');
		sb.append(timeInFrames2(tThreshold)).append('\t');
		sb.append(traces.length).append('\t');
		sb.append(filtered).append('\t');
		sb.append(singles).append('\t');
		sb.append(traces.length - singles).append('\t');
		for (int i = 0; i < stats.length; i++)
		{
			sb.append(Utils.rounded(stats[i].getMean(), 3)).append('\t');
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
					idList[count++] = Utils.showHistogram(TITLE, (StoredDataStatistics) stats[i], NAMES[i],
							(integerDisplay[i]) ? 1 : 0, (settings.removeOutliers || alwaysRemoveOutliers[i]) ? 2 : 0,
							settings.histogramBins);
					requireRetile = requireRetile || Utils.isNewWindow();
				}
			}

			if (count > 0 && requireRetile)
			{
				idList = Arrays.copyOf(idList, count);
				new WindowOrganiser().tileWindows(idList);
			}
		}

		if (settings.saveTraceData)
		{
			saveTraceData(stats);
		}

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
				"Dataset\tAlgorithm\tExposure time (ms)\tD-threshold (nm)\tT-threshold (s)\t(Frames)\tMolecules\tFiltered\tSingles\tClusters");
		for (int i = 0; i < NAMES.length; i++)
		{
			sb.append('\t').append(NAMES[i]);
		}
		return sb.toString();
	}

	private void saveTraceData(Statistics[] stats)
	{
		// Get the directory
		IJ.showStatus("Saving trace data");
		String directory = Utils.getDirectory("Trace_data_directory", settings.traceDataDirectory);
		if (directory != null)
		{
			settings.traceDataDirectory = directory;
			SettingsManager.saveSettings(globalSettings);
			for (int i = 0; i < NAMES.length; i++)
				saveTraceData((StoredDataStatistics) stats[i], NAMES[i], FILENAMES[i]);
		}
		IJ.showStatus("");
	}

	private void saveTraceData(StoredDataStatistics s, String name, String fileSuffix)
	{
		BufferedWriter file = null;
		try
		{
			file = new BufferedWriter(new FileWriter(settings.traceDataDirectory + TITLE + "." + fileSuffix + ".txt"));
			file.append(name);
			file.newLine();

			for (double d : s.getValues())
			{
				file.append(Utils.rounded(d, 4));
				file.newLine();
			}
		}
		catch (Exception e)
		{
			// Q. Add better handling of errors?
			e.printStackTrace();
			IJ.log("Failed to save trace data to results directory: " + settings.traceDataDirectory);
		}
		finally
		{
			if (file != null)
			{
				try
				{
					file.close();
				}
				catch (IOException e)
				{
					e.printStackTrace();
				}
			}
		}
	}

	private boolean showDialog()
	{
		TITLE = outputName + " Molecules";
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);

		globalSettings = SettingsManager.loadSettings();
		settings = globalSettings.getClusteringSettings();

		gd.addNumericField("Distance_Threshold (nm)", settings.distanceThreshold, 2);
		gd.addNumericField("Distance_Exclusion (nm)", settings.distanceExclusion, 2);
		gd.addNumericField("Time_Threshold", settings.getTimeThreshold(), 2);
		gd.addChoice("Time_unit", SettingsManager.getTimeUnitNames(), settings.getTimeUnit().ordinal());
		String[] traceModes = SettingsManager.getNames((Object[]) TraceManager.TraceMode.values());
		gd.addChoice("Trace_mode", traceModes,
				traceModes[ClusteringSettingsHelper.getTraceMode(settings.getTraceMode()).ordinal()]);
		gd.addNumericField("Pulse_interval (frames)", settings.pulseInterval, 0);
		gd.addNumericField("Pulse_window (frames)", settings.pulseWindow, 0);
		gd.addCheckbox("Split_pulses", settings.splitPulses);
		gd.addCheckbox("Optimise", settings.optimise);
		gd.addCheckbox("Save_traces", settings.saveTraces);
		gd.addCheckbox("Show_histograms", settings.showHistograms);
		gd.addCheckbox("Save_trace_data", settings.saveTraceData);
		//gd.addCheckbox("Refit_option", settings.refitOption);
		if (altKeyDown)
		{
			gd.addCheckbox("Debug", inputDebugMode);
		}

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

		// Store exposure time in seconds
		exposureTime = results.getCalibrationReader().getExposureTime() / 1000;

		return true;
	}

	private boolean readDialog(ExtendedGenericDialog gd)
	{
		inputOption = ResultsManager.getInputSource(gd);
		settings.distanceThreshold = gd.getNextNumber();
		settings.distanceExclusion = gd.getNextNumber();
		settings.setTimeThreshold(gd.getNextNumber());
		settings.setTimeUnit(gd.getNextChoiceIndex());
		settings.setTraceMode(gd.getNextChoiceIndex());
		settings.pulseInterval = (int) gd.getNextNumber();
		settings.pulseWindow = (int) gd.getNextNumber();
		settings.splitPulses = gd.getNextBoolean();
		settings.optimise = gd.getNextBoolean();
		settings.saveTraces = gd.getNextBoolean();
		settings.showHistograms = gd.getNextBoolean();
		settings.saveTraceData = gd.getNextBoolean();
		//settings.refitOption = gd.getNextBoolean();
		if (altKeyDown)
		{
			debugMode = inputDebugMode = gd.getNextBoolean();
		}

		if (gd.invalidNumber())
			return false;

		if (settings.showHistograms)
		{
			gd = new ExtendedGenericDialog(TITLE);
			gd.addMessage("Select the histograms to display");
			gd.addCheckbox("Remove_outliers", settings.removeOutliers);
			gd.addNumericField("Histogram_bins", settings.histogramBins, 0);
			for (int i = 0; i < displayHistograms.length; i++)
				gd.addCheckbox(NAMES[i].replace(' ', '_'), displayHistograms[i]);
			gd.showDialog();
			if (gd.wasCanceled())
				return false;
			settings.removeOutliers = gd.getNextBoolean();
			settings.histogramBins = (int) Math.abs(gd.getNextNumber());
			for (int i = 0; i < displayHistograms.length; i++)
				displayHistograms[i] = gd.getNextBoolean();
		}

		// Check arguments
		try
		{
			Parameters.isAboveZero("Distance threshold", settings.distanceThreshold);
			Parameters.isAboveZero("Time threshold", settings.getTimeThreshold());
			Parameters.isPositive("Pulse interval", settings.pulseInterval);
			Parameters.isPositive("Pulse window", settings.pulseWindow);
			Parameters.isAboveZero("Histogram bins", settings.histogramBins);
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		return true;
	}

	private boolean showClusterDialog()
	{
		TITLE = outputName + " Molecules";
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);

		globalSettings = SettingsManager.loadSettings();
		settings = globalSettings.getClusteringSettings();

		gd.addNumericField("Distance_Threshold (nm)", settings.distanceThreshold, 2);
		gd.addNumericField("Time_Threshold", settings.getTimeThreshold(), 2);
		gd.addChoice("Time_unit", SettingsManager.getTimeUnitNames(), settings.getTimeUnit().ordinal());
		String[] algorithm = SettingsManager.getNames((Object[]) ClusteringAlgorithm.values());
		gd.addChoice("Clustering_algorithm", algorithm, algorithm[ClusteringSettingsHelper
				.getClusteringAlgorithm(settings.getClusteringAlgorithm()).ordinal()]);
		gd.addNumericField("Pulse_interval (frames)", settings.pulseInterval, 0);
		gd.addCheckbox("Split_pulses", settings.splitPulses);
		gd.addCheckbox("Save_clusters", settings.saveTraces);
		gd.addCheckbox("Show_histograms", settings.showHistograms);
		gd.addCheckbox("Save_cluster_data", settings.saveTraceData);
		gd.addCheckbox("Refit_option", settings.refitOption);
		if (altKeyDown)
		{
			gd.addCheckbox("Debug", inputDebugMode);
		}

		gd.showDialog();

		if (gd.wasCanceled() || !readClusterDialog(gd))
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

		// Store exposure time in seconds
		exposureTime = results.getCalibrationReader().getExposureTime() / 1000;

		return true;
	}

	private boolean readClusterDialog(ExtendedGenericDialog gd)
	{
		inputOption = ResultsManager.getInputSource(gd);
		settings.distanceThreshold = gd.getNextNumber();
		settings.setTimeThreshold(gd.getNextNumber());
		settings.setTimeUnit(gd.getNextChoiceIndex());
		settings.setClusteringAlgorithm(gd.getNextChoiceIndex());
		settings.pulseInterval = (int) gd.getNextNumber();
		settings.splitPulses = gd.getNextBoolean();
		settings.saveTraces = gd.getNextBoolean();
		settings.showHistograms = gd.getNextBoolean();
		settings.saveTraceData = gd.getNextBoolean();
		settings.refitOption = gd.getNextBoolean();
		if (altKeyDown)
		{
			debugMode = inputDebugMode = gd.getNextBoolean();
		}

		if (gd.invalidNumber())
			return false;

		if (settings.showHistograms)
		{
			gd = new ExtendedGenericDialog(TITLE);
			gd.addMessage("Select the histograms to display");
			gd.addCheckbox("Remove_outliers", settings.removeOutliers);
			gd.addNumericField("Histogram_bins", settings.histogramBins, 0);
			for (int i = 0; i < displayHistograms.length; i++)
				gd.addCheckbox(NAMES[i].replace(' ', '_'), displayHistograms[i]);
			gd.showDialog();
			if (gd.wasCanceled())
				return false;
			settings.removeOutliers = gd.getNextBoolean();
			settings.histogramBins = (int) Math.abs(gd.getNextNumber());
			for (int i = 0; i < displayHistograms.length; i++)
				displayHistograms[i] = gd.getNextBoolean();
		}

		// Check arguments
		try
		{
			Parameters.isAboveZero("Distance threshold", settings.distanceThreshold);
			ClusteringAlgorithm clusteringAlgorithm = ClusteringSettingsHelper
					.getClusteringAlgorithm(settings.getClusteringAlgorithm());
			if (clusteringAlgorithm == ClusteringAlgorithm.CENTROID_LINKAGE_DISTANCE_PRIORITY ||
					clusteringAlgorithm == ClusteringAlgorithm.CENTROID_LINKAGE_TIME_PRIORITY)
			{
				Parameters.isAboveZero("Time threshold", settings.getTimeThreshold());
				Parameters.isPositive("Pulse interval", settings.pulseInterval);
			}
			Parameters.isAboveZero("Histogram bins", settings.histogramBins);
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		return true;
	}

	private void runOptimiser(TraceManager manager)
	{
		// Get an estimate of the number of molecules without blinking
		Statistics stats = new Statistics();
		final double nmPerPixel = this.results.getNmPerPixel();
		PrecisionResultProcedure pp = new PrecisionResultProcedure(results);
		pp.getPrecision();
		stats.add(pp.precision);
		// Use twice the precision to get the initial distance threshold

		// Use 2.5x sigma as per the PC-PALM protocol in Sengupta, et al (2013) Nature Protocols 8, 345
		double dEstimate = stats.getMean() * 2.5 / nmPerPixel;
		int n = manager.traceMolecules(dEstimate, 1);
		//for (double d : new double[] { 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4 })
		//	System.out.printf("d=%.2f, estimate=%d\n", d,
		//			manager.traceMolecules(stats.getMean() * d / this.results.getNmPerPixel(), 1));

		if (!getParameters(n, dEstimate))
			return;

		// TODO - Convert the distance threshold to use nm instead of pixels?		
		List<double[]> results = runTracing(manager, settings.minDistanceThreshold, settings.maxDistanceThreshold,
				settings.minTimeThreshold, settings.maxTimeThreshold, settings.optimiserSteps);

		// Compute fractional difference from the true value:
		// Use blinking rate directly or the estimated number of molecules
		double nReference;
		int statistic;
		if (optimiseBlinkingRate)
		{
			nReference = settings.blinkingRate;
			statistic = 3;
			IJ.log(String.format("Estimating blinking rate: %.2f", nReference));
		}
		else
		{
			nReference = n / settings.blinkingRate;
			statistic = 2;
			IJ.log(String.format("Estimating number of molecules: %d / %.2f = %.2f", n, settings.blinkingRate,
					nReference));
		}

		for (double[] result : results)
		{
			//System.out.printf("%g %g = %g\n", result[0], result[1], result[2]);
			if (optimiseBlinkingRate)
				result[2] = (nReference - result[statistic]) / nReference;
			else
				result[2] = (result[statistic] - nReference) / nReference;
		}

		// Locate the optimal parameters with a fit of the zero contour
		boolean found = findOptimalParameters(results);

		createPlotResults(results);

		if (!found)
			return;

		// Make fractional difference absolute so that lowest is best
		for (double[] result : results)
			result[2] = Math.abs(result[2]);

		// Set the optimal thresholds using the lowest value
		double[] best = new double[] { 0, 0, Double.MAX_VALUE };
		for (double[] result : results)
			if (best[2] > result[2])
				best = result;

		settings.distanceThreshold = best[0];

		// The optimiser works using frames so convert back to the correct units
		TypeConverter<TimeUnit> convert = UnitConverterFactory.createConverter(TimeUnit.FRAME, settings.getTimeUnit(),
				exposureTime);
		settings.setTimeThreshold(convert.convert(best[1]));

		IJ.log(String.format("Optimal fractional difference @ D-threshold=%g, T-threshold=%f (%d frames)",
				settings.distanceThreshold, timeInSeconds(settings), timeInFrames(settings)));
		SettingsManager.saveSettings(globalSettings);
	}

	private boolean getParameters(int n, double d)
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE + " Optimiser");

		String msg = String.format("Estimate %d molecules at d=%f, t=1", n, d);
		IJ.log(msg);
		gd.addMessage(msg);
		gd.addNumericField("Min_Distance_Threshold (px)", settings.minDistanceThreshold, 2);
		gd.addNumericField("Max_Distance_Threshold (px)", settings.maxDistanceThreshold, 2);
		gd.addNumericField("Min_Time_Threshold (frames)", settings.minTimeThreshold, 0);
		gd.addNumericField("Max_Time_Threshold (frames)", settings.maxTimeThreshold, 0);
		gd.addSlider("Steps", 1, 20, settings.optimiserSteps);
		gd.addNumericField("Blinking_rate", settings.blinkingRate, 2);
		String[] plotNames = SettingsManager.getNames((Object[]) ClusteringSettingsHelper.OptimiserPlot.values());
		gd.addChoice("Plot", plotNames,
				plotNames[ClusteringSettingsHelper.getOptimiserPlot(settings.getOptimiserPlot()).ordinal()]);
		if (altKeyDown)
			gd.addCheckbox("Optimise_blinking", inputOptimiseBlinkingRate);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		settings.minDistanceThreshold = gd.getNextNumber();
		settings.maxDistanceThreshold = gd.getNextNumber();
		settings.minTimeThreshold = (int) gd.getNextNumber();
		settings.maxTimeThreshold = (int) gd.getNextNumber();
		settings.optimiserSteps = (int) gd.getNextNumber();
		settings.blinkingRate = gd.getNextNumber();
		settings.setOptimiserPlot(gd.getNextChoiceIndex());
		if (altKeyDown)
		{
			optimiseBlinkingRate = inputOptimiseBlinkingRate = gd.getNextBoolean();
		}

		if (gd.invalidNumber())
			return false;

		if (settings.minDistanceThreshold < 0)
			settings.minDistanceThreshold = 0;
		if (settings.maxDistanceThreshold < settings.minDistanceThreshold)
			settings.maxDistanceThreshold = settings.minDistanceThreshold;
		if (settings.minTimeThreshold < 0)
			settings.minTimeThreshold = 0;
		if (settings.maxTimeThreshold < settings.minTimeThreshold)
			settings.maxTimeThreshold = settings.minTimeThreshold;
		if (settings.optimiserSteps < 0)
			settings.optimiserSteps = 1;
		if (settings.blinkingRate < MIN_BLINKING_RATE)
		{
			IJ.error(gd.getTitle(), "Blinking rate must be above " + MIN_BLINKING_RATE);
			return false;
		}

		if (settings.minDistanceThreshold == settings.maxDistanceThreshold &&
				settings.minTimeThreshold == settings.maxTimeThreshold)
		{
			IJ.error(gd.getTitle(), "Nothing to optimise");
			return false;
		}

		SettingsManager.saveSettings(globalSettings);

		return true;
	}

	private int timeInFrames2(double timeInSeconds)
	{
		return (int) Math.round(timeInSeconds / exposureTime);
	}

	private int timeInFrames(ClusteringSettings settings)
	{
		return (int) Math.round(timeIn(settings, TimeUnit.FRAME));
	}

	private double timeInSeconds(ClusteringSettings settings)
	{
		return timeIn(settings, TimeUnit.SECOND);
	}

	private double timeIn(ClusteringSettings settings, TimeUnit timeUnit)
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
	 * @param minDistanceThreshold
	 * @param maxDistanceThreshold
	 * @param minTimeThreshold
	 * @param maxTimeThreshold
	 * @param optimiserSteps
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
	 * @param minDistanceThreshold
	 * @param maxDistanceThreshold
	 * @param minTimeThreshold
	 * @param maxTimeThreshold
	 * @param optimiserSteps
	 * @return a list of [distance,time,N traces,blinking rate]
	 */
	public List<double[]> runTracing(TraceManager manager, double minDistanceThreshold, double maxDistanceThreshold,
			int minTimeThreshold, int maxTimeThreshold, int optimiserSteps)
	{
		dThresholds = getIntervals(minDistanceThreshold, maxDistanceThreshold, optimiserSteps);
		tThresholds = convert(getIntervals(minTimeThreshold, maxTimeThreshold, optimiserSteps));

		int total = dThresholds.length * tThresholds.length;
		ArrayList<double[]> results = new ArrayList<double[]>(total);

		IJ.showStatus("Optimising tracing (" + total + " steps) ...");

		if (debugMode)
			IJ.log("Optimising tracing ...");

		int step = 0;
		for (double d : dThresholds)
			for (int t : tThresholds)
			{
				IJ.showProgress(step++, total);
				int n = manager.traceMolecules(d, t);
				results.add(new double[] { d, t, n, getBlinkingRate(manager.getTraces()) });
				if (debugMode)
				{
					summarise(manager.getTraces(), manager.getTotalFiltered(), d, t);
				}
			}

		if (debugMode)
			IJ.log("-=-=-=-");

		IJ.showStatus("");
		IJ.showProgress(1.0);

		return results;
	}

	private double getBlinkingRate(Trace[] traces)
	{
		SummaryStatistics stats = new SummaryStatistics();
		for (Trace trace : traces)
			stats.addValue(trace.getNBlinks());
		double blinkingRate = stats.getMean();
		return blinkingRate;
	}

	private double[] getIntervals(double min, double max, int optimiserSteps)
	{
		if (max < min)
		{
			double tmp = max;
			max = min;
			min = tmp;
		}
		double range = max - min;
		if (range == 0)
			return new double[] { min };

		double[] values = new double[optimiserSteps + ((min != 0) ? 1 : 0)];
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

	private int[] convert(double[] intervals)
	{
		TIntHashSet set = new TIntHashSet(intervals.length);
		for (double d : intervals)
			set.add((int) Math.round(d));

		set.remove(0); // Do not allow zero

		int[] values = set.toArray();
		Arrays.sort(values);
		return values;
	}

	/**
	 * Find the contour that intersects zero on the fractional difference plot.
	 * Find the point on the contour nearest the origin.
	 * 
	 * @param results
	 */
	private boolean findOptimalParameters(List<double[]> results)
	{
		// This method only works if there are many results and if the results 
		// cover enough of the search space to go from above zero (i.e. not enough traces)
		// to below zero (i.e. too many traces)

		int maxx = tThresholds.length;
		int maxy = dThresholds.length;

		// --------
		// Find zero crossings using linear interpolation 
		zeroCrossingPoints = new ArrayList<double[]>();
		// --------

		// Pass across all time points
		boolean noZeroCrossingAtT0 = false;
		boolean noZeroCrossingAtTN = false;
		for (int x = 0; x < maxx; x++)
		{
			// Find zero crossings on distance points
			double[] data = new double[maxy];
			for (int y = 0; y < maxy; y++)
			{
				int i = y * maxx + x;
				double[] result = results.get(i);
				data[y] = result[2];
			}
			double zeroCrossing = findZeroCrossing(data, dThresholds);
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
		final double maxTimeThresholdInFrames = settings.maxTimeThreshold;
		// The optimiser works using frames so convert back to the correct units
		TypeConverter<TimeUnit> convert = UnitConverterFactory.createConverter(TimeUnit.FRAME, settings.getTimeUnit(),
				exposureTime);

		for (double[] point : zeroCrossingPoints)
		{
			double dx = point[0] / maxTimeThresholdInFrames;
			double dy = point[1] / settings.maxDistanceThreshold;
			double d = dx * dx + dy * dy;
			if (d < minD)
			{
				minD = d;
				settings.distanceThreshold = point[1];
				settings.setTimeThreshold(convert.convert(point[0]));
			}
		}

		// --------
		// Add more points to make the plotted line look better when showing the plot.
		// --------

		// Pass across all distance points
		boolean noZeroCrossingAtD0 = false;
		boolean noZeroCrossingAtDN = false;
		double[] tThresholdsD = toDouble(tThresholds);
		for (int y = 0; y < maxy; y++)
		{
			// Find zero crossings on time points
			double[] data = new double[maxx];
			for (int x = 0; x < maxx; x++)
			{
				int i = y * maxx + x;
				double[] result = results.get(i);
				data[x] = result[2];
			}
			double zeroCrossing = findZeroCrossing(data, tThresholdsD);
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
		StringBuilder sb = new StringBuilder();
		boolean reduceTime = false;
		boolean reduceDistance = false;
		if (noZeroCrossingAtDN && settings.minTimeThreshold > 1)
		{
			sb.append(" * No zero crossing at max distance\n");
			reduceTime = true;
		}
		if (noZeroCrossingAtTN && settings.minDistanceThreshold > 0)
		{
			sb.append(" * No zero crossing at max time\n");
			reduceDistance = true;
		}
		if (!noZeroCrossingAtD0 && settings.minDistanceThreshold > 0)
		{
			sb.append(" * Zero crossing at min distance\n");
			reduceDistance = true;
		}
		if (!noZeroCrossingAtT0 && settings.minTimeThreshold > 1)
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

	private double findZeroCrossing(double[] data, double[] axis)
	{
		if (data[0] < 0)
			return -1;
		for (int i = 1; i < data.length; i++)
		{
			if (data[i] < 0)
			{
				double fraction = data[i - 1] / (data[i - 1] - data[i]);
				return fraction * axis[i] + (1 - fraction) * axis[i - 1];
			}
		}

		return -1;
	}

	private void sortPoints()
	{
		// Sort by x coord, then y
		Collections.sort(zeroCrossingPoints, new Comparator<double[]>()
		{
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
		double[] x = new double[zeroCrossingPoints.size()];
		double[] y = new double[zeroCrossingPoints.size()];
		for (int i = 0; i < x.length; i++)
		{
			double[] point = zeroCrossingPoints.get(i);
			x[i] = point[0];
			y[i] = point[1];
		}
		PolynomialSplineFunction fx = new SplineInterpolator().interpolate(x, y);
		double minX = x[0];
		double maxX = x[x.length - 1];
		double xinc = (maxX - minX) / 50;
		for (minX = minX + xinc; minX < maxX; minX += xinc)
		{
			zeroCrossingPoints.add(new double[] { minX, fx.value(minX) });
		}
		sortPoints();
	}

	/**
	 * Build an image using the values within the results to set X,Y and value
	 * 
	 * @param results
	 */
	private void createPlotResults(List<double[]> results)
	{
		int w = 400, h = 400;
		switch (ClusteringSettingsHelper.getOptimiserPlot(settings.getOptimiserPlot()))
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
		double xRange = getRange(settings.maxTimeThreshold, settings.minTimeThreshold, origX, w);
		double yRange = getRange(settings.maxDistanceThreshold, settings.minDistanceThreshold, origY, h);
		cal.pixelWidth = xRange / w;
		cal.pixelHeight = yRange / h;
		cal.xOrigin = origX - settings.minTimeThreshold / cal.pixelWidth;
		cal.yOrigin = origY - settings.minDistanceThreshold / cal.pixelHeight;
		cal.setXUnit("frame");
		cal.setYUnit("pixel");

		showPlot();
	}

	/**
	 * Shows the plot
	 */
	private void showPlot()
	{
		if (ClusteringSettingsHelper.getOptimiserPlot(settings.getOptimiserPlot()) == OptimiserPlot.NONE)
			return;

		// Display the image
		String title = TITLE + ": | N - N_actual | / N_actual";
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
			LutLoader lut = new LutLoader();
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
		Calibration cal = imp.getCalibration();
		int nPoints = zeroCrossingPoints.size();
		float[] xPoints = new float[nPoints];
		float[] yPoints = new float[nPoints];
		for (int i = 0; i < nPoints; i++)
		{
			double[] point = zeroCrossingPoints.get(i);
			// Convert to pixel coordinates. 
			xPoints[i] = (float) (cal.xOrigin + (point[0] / cal.pixelWidth));
			yPoints[i] = (float) (cal.yOrigin + (point[1] / cal.pixelHeight));
		}
		roi = new PolygonRoi(xPoints, yPoints, nPoints, PolygonRoi.POLYLINE);
		imp.setRoi(roi);
	}

	private FloatProcessor createNNPlot(List<double[]> results, int w, int h)
	{
		FloatProcessor fp = new FloatProcessor(w, h);

		// Create lookup table that map the tested threshold values to a position in the image
		int[] xLookup = createLookup(tThresholds, settings.minTimeThreshold, w);
		int[] yLookup = createLookup(dThresholds, settings.minDistanceThreshold, h);
		origX = (settings.minTimeThreshold != 0) ? xLookup[1] : 0;
		origY = (settings.minDistanceThreshold != 0) ? yLookup[1] : 0;

		int gridWidth = tThresholds.length;
		int gridHeight = dThresholds.length;
		for (int y = 0, i = 0; y < gridHeight; y++)
		{
			for (int x = 0; x < gridWidth; x++, i++)
			{
				int x1 = xLookup[x];
				int x2 = xLookup[x + 1];
				int y1 = yLookup[y];
				int y2 = yLookup[y + 1];
				double[] result = results.get(i);
				fp.setValue(Math.abs(result[2]));
				fp.setRoi(x1, y1, x2 - x1, y2 - y1);
				fp.fill();
			}
		}
		return fp;
	}

	private int[] createLookup(int[] values, int min, int scale)
	{
		double[] newValues = toDouble(values);
		return createLookup(newValues, min, scale);
	}

	private double[] toDouble(int[] values)
	{
		double[] newValues = new double[values.length];
		for (int i = 0; i < values.length; i++)
		{
			newValues[i] = values[i];
		}
		return newValues;
	}

	private int[] createLookup(double[] values, double min, int scale)
	{
		// To allow the lowest result to be plotted, add space at the edge 
		// equal to the next interval
		if (min != 0 && values.length > 1)
		{
			min -= values[1] - values[0];
		}

		int[] lookup = new int[values.length + 1];
		double range = values[values.length - 1] - min;
		double scaleFactor = scale / range;
		for (int i = 1; i < values.length; i++)
		{
			lookup[i] = (int) Math.round(scaleFactor * (values[i - 1] - min));
		}
		lookup[values.length] = scale;
		return lookup;
	}

	private FloatProcessor createBilinearPlot(List<double[]> results, int w, int h)
	{
		FloatProcessor fp = new FloatProcessor(w, h);

		// Create lookup table that map the tested threshold values to a position in the image
		int[] xLookup = createLookup(tThresholds, settings.minTimeThreshold, w);
		int[] yLookup = createLookup(dThresholds, settings.minDistanceThreshold, h);
		origX = (settings.minTimeThreshold != 0) ? xLookup[1] : 0;
		origY = (settings.minDistanceThreshold != 0) ? yLookup[1] : 0;

		int gridWidth = tThresholds.length;
		int gridHeight = dThresholds.length;
		for (int y = 0, prevY = 0; y < gridHeight; y++)
		{
			for (int x = 0, prevX = 0; x < gridWidth; x++)
			{
				// Get the 4 flanking values
				double x1y1 = results.get(prevY * gridWidth + prevX)[2];
				double x1y2 = results.get(y * gridWidth + prevX)[2];
				double x2y1 = results.get(prevY * gridWidth + x)[2];
				double x2y2 = results.get(y * gridWidth + x)[2];

				// Pixel range
				int x1 = xLookup[x];
				int x2 = xLookup[x + 1];
				int y1 = yLookup[y];
				int y2 = yLookup[y + 1];

				double xRange = x2 - x1;
				double yRange = y2 - y1;

				for (int yy = y1; yy < y2; yy++)
				{
					double yFraction = (yy - y1) / yRange;
					for (int xx = x1; xx < x2; xx++)
					{
						// Interpolate
						double xFraction = (xx - x1) / xRange;
						double v1 = x1y1 * (1 - xFraction) + x2y1 * xFraction;
						double v2 = x1y2 * (1 - xFraction) + x2y2 * xFraction;
						double value = v1 * (1 - yFraction) + v2 * yFraction;
						fp.setf(xx, yy, (float) value);
					}
				}

				prevX = x;
			}
			prevY = y;
		}

		// Convert to absolute for easier visualisation
		float[] data = (float[]) fp.getPixels();
		for (int i = 0; i < data.length; i++)
			data[i] = Math.abs(data[i]);

		return fp;
	}

	private double getRange(double max, double min, int orig, int w)
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
		//		gd.addSlider("Width_factor", 1.01, 5, fitConfig.getWidthFactor());
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
		float[] centre = trace.getCentroid(CentroidMethod.SIGNAL_WEIGHTED);
		int minX = (int) Math.floor(centre[0] - fitWidth);
		int maxX = (int) Math.ceil(centre[0] + fitWidth);
		int minY = (int) Math.floor(centre[1] - fitWidth);
		int maxY = (int) Math.ceil(centre[1] + fitWidth);

		// Account for crops at the edge of the image
		minX = FastMath.max(0, minX);
		maxX = FastMath.min(w, maxX);
		minY = FastMath.max(0, minY);
		maxY = FastMath.min(h, maxY);

		int width = maxX - minX;
		int height = maxY - minY;
		if (width <= 0 || height <= 0)
		{
			// The centre must be outside the image width and height
			return null;
		}
		bounds.x = minX;
		bounds.y = minY;
		bounds.width = width;
		bounds.height = height;

		if (createStack)
			slices = new ImageStack(width, height);

		// Combine the images. Subtract the fitted background to zero the image.
		float[] data = new float[width * height];
		float sumBackground = 0;
		double noise = 0;
		for (PeakResult result : trace.getPoints())
		{
			noise += result.noise * result.noise;

			float[] sourceData = source.get(result.getFrame(), bounds);
			final float background = result.getBackground();
			sumBackground += background;
			for (int i = 0; i < data.length; i++)
			{
				data[i] += sourceData[i] - background;
			}
			if (createStack)
				slices.addSlice(new FloatProcessor(width, height, sourceData, null));
		}
		if (createStack)
		{
			// Add a final image that is the average of the individual slices. This allows
			// it to be visualised in the same intensity scale.
			float[] data2 = Arrays.copyOf(data, data.length);
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
	private double getCombinedNoise(Trace trace)
	{
		double noise = 0;
		for (PeakResult result : trace.getPoints())
		{
			noise += result.noise * result.noise;
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
