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

import java.awt.AWTEvent;
import java.awt.Rectangle;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.regex.Pattern;

import ags.utils.dataStructures.trees.secondGenKD.IntResultHeap;
import gdsc.core.generics.ConcurrentMonoStack;
import gdsc.core.ij.BufferedTextWindow;
import gdsc.core.ij.IJTrackProgress;
import gdsc.core.ij.Utils;
import gdsc.core.logging.Ticker;
import gdsc.core.utils.BooleanArray;
import gdsc.core.utils.BooleanRollingArray;
import gdsc.core.utils.Maths;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.core.utils.Sort;
import gdsc.core.utils.TextUtils;
import gdsc.core.utils.TurboList;

import gdsc.smlm.data.NamedObject;
import gdsc.smlm.data.config.GUIProtos.FailCountManagerSettings;
import gdsc.smlm.data.config.ResultsProtos.ResultsSettings;
import gdsc.smlm.engine.FitEngine;
import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.engine.FitJob.Status;
import gdsc.smlm.engine.FitParameters;
import gdsc.smlm.engine.FitParameters.FitTask;
import gdsc.smlm.engine.ParameterisedFitJob;
import gdsc.smlm.ij.IJImageSource;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.results.ImageSource;
import gdsc.smlm.results.count.ConsecutiveFailCounter;
import gdsc.smlm.results.count.FailCounter;
import gdsc.smlm.results.count.PassRateFailCounter;
import gdsc.smlm.results.count.ResettingFailCounter;
import gdsc.smlm.results.count.RollingWindowFailCounter;
import gdsc.smlm.results.count.WeightedFailCounter;
import gnu.trove.list.array.TByteArrayList;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.ExtendedGenericDialog;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingExtendedGenericDialog;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.plugin.WindowOrganiser;
import ij.text.TextWindow;

/**
 * This plugin handles generation and analysis of fail counts to optimise the stopping criteria for a sequential
 * analysis routine.
 */
public class FailCountManager implements PlugIn
{
	private static final String TITLE = "Fail Count Manager";

	//@formatter:off
	private enum FailCountOption implements NamedObject 
	{
		CREATE_DATA { @Override
		public String getName() { return "Create Data"; } },
		LOAD_DATA { @Override
		public String getName() { return "Load Data"; } },
		SAVE_DATA { @Override
		public String getName() { return "Save Data"; } },
		PLOT_DATA { @Override
		public String getName() { return "Plot Data"; } },
		ANALYSE_DATA { @Override
		public String getName() { return "Analyse Data"; } },
		;

		@Override
		public String getShortName()
		{
			return getName();
		}
		
		public static FailCountOption forOrdinal(int ordinal)
		{
			FailCountOption[] values = FailCountOption.values();
			if (ordinal < 0 || ordinal >= values.length)
					ordinal = 0;
			return values[ordinal]; 
		}
	};
	//@formatter:on
	private static String[] OPTIONS = SettingsManager.getNames((Object[]) FailCountOption.values());

	/**
	 * Hold the fail count data for a single sequential analysis routine
	 */
	private static class FailCountData
	{
		/** The id of the data. */
		public final int id;

		/** The results (pass/fail). */
		private boolean[] results;

		private int maxConsFailCount = -1;

		// These are for plotting so we use float not int
		private float[] candidate = null;
		private float[] consFailCount = null;
		private float[] passCount = null;
		private float[] passRate = null;

		/** The number of results to process before a fail counter is not OK. Used to score a fail counter */
		private int target;
		private float targetPassCount;
		private double maxScore;

		public FailCountData(int id, boolean[] results)
		{
			this.id = id;
			this.results = results;
		}

		public int getMaxConsecutiveFailCount()
		{
			if (maxConsFailCount == -1)
			{
				int consFail = 0;
				int max = 0;
				for (int i = 0; i < results.length; i++)
				{
					if (results[i])
					{
						consFail = 0;
					}
					else
					{
						consFail++;
						if (max < consFail)
							max = consFail;
					}
				}
				maxConsFailCount = max;
			}
			return maxConsFailCount;
		}

		public int getPassCount()
		{
			createData();
			return (int) passCount[passCount.length - 1];
		}

		public int getFailCount()
		{
			return results.length - getPassCount();
		}

		public void createData()
		{
			if (candidate == null)
				initialiseData();
		}

		private synchronized void initialiseData()
		{
			if (candidate != null)
				return;
			int pass = 0;
			int consFail = 0;

			int size = results.length;
			candidate = new float[size];
			passCount = new float[size];
			passRate = new float[size];
			consFailCount = new float[size];
			for (int i = 0; i < size; i++)
			{
				if (results[i])
				{
					pass++;
					consFail = 0;
				}
				else
				{
					consFail++;
					// Only set this when non-zero
					consFailCount[i] = consFail;
				}
				candidate[i] = i + 1;
				passCount[i] = pass;
				passRate[i] = (float) pass / (i + 1);
			}
		}

		public float[] getRollingFailCount(int rollingWindow)
		{
			BooleanRollingArray c = new BooleanRollingArray(rollingWindow);
			int size = results.length;
			float[] failCount = new float[size];
			for (int i = 0; i < size; i++)
			{
				c.add(results[i]);
				failCount[i] = c.getFalseCount();
			}
			return failCount;
		}

		public float[] getWeightedFailCount(int passWeight, int failWeight)
		{
			WeightedFailCounter c = WeightedFailCounter.create(Integer.MAX_VALUE, failWeight, passWeight);
			int size = results.length;
			float[] failCount = new float[size];
			for (int i = 0; i < size; i++)
			{
				c.addResult(results[i]);
				failCount[i] = c.getFailCount();
			}
			return failCount;
		}

		public float[] getResettingFailCount(double resetFraction)
		{
			ResettingFailCounter c = ResettingFailCounter.create(Integer.MAX_VALUE, resetFraction);
			int size = results.length;
			float[] failCount = new float[size];
			for (int i = 0; i < size; i++)
			{
				c.addResult(results[i]);
				failCount[i] = c.getFailCount();
			}
			return failCount;
		}

		public void initialiseAnalysis(double targetPassFraction)
		{
			initialiseData();
			targetPassCount = Math.round(passCount[passCount.length - 1] * targetPassFraction);
			int size = results.length;
			target = 1;
			while (target <= size)
			{
				if (passCount[target - 1] >= targetPassCount)
					break;
				target++;
			}
			maxScore = score(size);
		}

		public double score(FailCounter counter)
		{
			int size = results.length;
			int i = 0;
			while (i < size)
			{
				if (results[i])
					counter.pass();
				else
					counter.fail();
				i++;
				if (!counter.isOK())
					return score(i);
			}
			return maxScore;
		}

		public double score(int n)
		{
			if (n == target)
				return 0; // Perfect
			if (n < target)
			{
				// Penalise stopping too early.
				float remaining = (targetPassCount - passCount[n - 1]) / targetPassCount;
				// This has a score from 0 to 1.
				//return remaining * remaining;
				return remaining;
			}
			else
			{
				// Penalise running too long.
				// Overrun will be above 0.
				float overrun = ((float) (n - target)) / target;

				// This has a score from 0 to Infinity.
				//return overrun * overrun;
				return overrun;

				// This has a score from 0 to Infinity but does not heavily penalise large overrun
				//if (overrun < 1)
				//	 // So the gradient is 1 at x=1, f(x)=0.5*x^2, f'(x)=x
				//	return overrun * overrun / 2.0;
				//else
				//	return overrun;

				// This has a score from 0 to 1. This equally weights under/overrun.
				// However very long overrun is not penalised due to the exponential.
				//  0	0
				//	0.5	0.3934693403
				//	1	0.6321205588
				//	2	0.8646647168
				//	3	0.9502129316
				//	4	0.9816843611
				//	5	0.993262053
				//	10	0.9999546001
				//return 1.0 - FastMath.exp(-overrun);
			}
		}
	}

	private static TurboList<FailCountData> failCountData = new TurboList<FailCountData>(1);
	private static TextWindow resultsWindow = null;

	private FailCountManagerSettings.Builder settings;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	@Override
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		settings = SettingsManager.readFailCountManagerSettings(0).toBuilder();

		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addChoice("Option", OPTIONS, settings.getOption());
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		settings.setOption(gd.getNextChoiceIndex());

		FailCountOption option = FailCountOption.forOrdinal(settings.getOption());
		switch (option)
		{
			case CREATE_DATA:
				createData();
				break;
			case LOAD_DATA:
				loadData();
				break;
			case SAVE_DATA:
				saveData();
				break;
			case PLOT_DATA:
				plotData();
				break;
			case ANALYSE_DATA:
				analyseData();
				break;
			default:
				throw new IllegalStateException("Unknown option: " + option);
		}

		SettingsManager.writeSettings(settings);
	}

	/**
	 * Creates the fail count data by running fitting on the current image.
	 */
	private void createData()
	{
		ImagePlus imp = WindowManager.getCurrentImage();
		if (imp == null)
		{
			IJ.error(TITLE, "No image for fitting");
			return;
		}

		if (!showCreateDataDialog(imp))
			return;

		// Get the current fit configuration
		Configuration c = new Configuration();
		if (!c.showDialog(false))
			return;

		FitEngineConfiguration fitConfig = c.getFitEngineConfiguration();
		// Update stopping criteria.
		fitConfig.resetFailCounter();
		fitConfig.setFailuresLimit(settings.getFailCountLimit());

		ImageSource source = new IJImageSource(imp);
		PeakFit peakFit = new PeakFit(fitConfig, ResultsSettings.getDefaultInstance());
		peakFit.setResultsSuffix("(FailCountAnalysis)");
		if (!peakFit.initialise(source, null, false))
		{
			IJ.error(TITLE, "Failed to initialise the fit engine");
			return;
		}
		FitEngine engine = peakFit.createFitEngine();

		Rectangle bounds = new Rectangle(source.getWidth(), source.getHeight());

		// Run 
		int totalFrames = Math.min(source.getFrames(), settings.getMaxFrames());
		final int step = Utils.getProgressInterval(totalFrames);
		IJ.showProgress(0);
		boolean shutdown = false;
		int slice = 0;
		TurboList<ParameterisedFitJob> jobs = new TurboList<ParameterisedFitJob>(totalFrames);
		while (!shutdown && slice < totalFrames)
		{
			float[] data = source.next();
			if (data == null)
				break;

			if (slice++ % step == 0)
			{
				if (Utils.showStatus("Fitting slice: " + slice + " / " + totalFrames))
					IJ.showProgress(slice, totalFrames);
			}

			ParameterisedFitJob job = createJob(source.getStartFrameNumber(), data, bounds);
			jobs.addf(job);
			engine.run(job);

			shutdown = escapePressed();
		}

		Utils.showStatus("Extracting fail count data");
		engine.end(shutdown);
		IJ.showProgress(1);
		source.close();

		// Extract the fail count data
		TurboList<FailCountData> failCountData = new TurboList<FailCountData>(jobs.size());
		for (int i = 0; i < jobs.size(); i++)
		{
			ParameterisedFitJob job = jobs.getf(i);
			if (job.getStatus() == Status.FINISHED)
			{
				FitParameters fitParams = job.getFitParameters();

				// Find the last success
				boolean[] results = fitParams.pass;
				int end = results.length - 1;
				while (end > 0 && !results[end])
					end--;
				// Add on the configured fail count limit
				end = Math.min(end + 1 + settings.getFailCountLimit(), results.length);
				results = Arrays.copyOf(results, end);
				failCountData.add(new FailCountData(job.getSlice(), results));
			}
		}
		FailCountManager.failCountData = failCountData;
		Utils.showStatus("");

		// Save for the future
		if (settings.getSaveAfterFitting())
			saveData();
	}

	private boolean showCreateDataDialog(ImagePlus imp)
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addMessage(TextUtils.wrap("Run the fit engine on the current image to generate " +
				"pass/fail data for sequential candidates in each frame. A second dialog will " +
				"be shown to check the fit settings.", 80));
		gd.addSlider("Max_frames", 1, imp.getStackSize(), settings.getMaxFrames());
		gd.addNumericField("Fail_count_limit", settings.getFailCountLimit(), 0);
		gd.addCheckbox("Save", settings.getSaveAfterFitting());
		gd.showDialog();
		if (gd.wasCanceled())
			return false;
		settings.setMaxFrames((int) gd.getNextNumber());
		settings.setFailCountLimit((int) gd.getNextNumber());
		settings.setSaveAfterFitting(gd.getNextBoolean());
		try
		{
			Parameters.isAboveZero("Max frames", settings.getMaxFrames());
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}
		return true;
	}

	private ParameterisedFitJob createJob(int startFrame, float[] data, Rectangle bounds)
	{
		FitParameters fitParams = new FitParameters();
		fitParams.fitTask = FitTask.PSF_FITTING;
		// Signal that the fail count should be recorded
		fitParams.pass = new boolean[0];
		return new ParameterisedFitJob(fitParams, startFrame, data, bounds);
	}

	private boolean escapePressed()
	{
		if (IJ.escapePressed())
		{
			IJ.log(TITLE + " stopping ...");
			IJ.beep();
			return true;
		}
		return false;
	}

	/**
	 * Load the data from a file.
	 */
	private void loadData()
	{
		String filename = Utils.getFilename("Fail_count_data_filename", settings.getFilename());
		if (filename == null)
			return;
		settings.setFilename(filename);
		TurboList<FailCountData> failCountData = new TurboList<FailCountData>();

		BufferedReader br = null;
		Pattern pattern = Pattern.compile("[\t, ]+");
		try
		{
			br = new BufferedReader(new FileReader(filename));
			// Ignore the first line
			String line = br.readLine();
			BooleanArray array = new BooleanArray(100);
			int lastId = 0;
			int lastCandidate = 0;
			while ((line = br.readLine()) != null)
			{
				String[] data = pattern.split(line);
				if (data.length != 3)
					throw new IOException("Require 3 fields in the data");

				int id = Integer.parseInt(data[0]);
				if (id < 1)
					throw new IOException("ID must be strictly positive");
				int candidate = Integer.parseInt(data[1]);
				if (candidate < 1)
					throw new IOException("Candidate must be strictly positive");
				boolean ok = guessStatus(data[2]);

				if (lastId != id)
				{
					if (array.size() > 0)
					{
						failCountData.add(new FailCountData(lastId, array.toArray()));
						array.clear();
					}
					if (candidate != 1)
						throw new IOException("Candidate must start at 1");
					lastId = id;
					lastCandidate = candidate - 1; // Ensure continuous
				}
				// Require continuous sequence
				if (candidate - lastCandidate == 1)
				{
					array.add(ok);
					lastCandidate = candidate;
				}
				else
				{
					// Make impossible to add any more for this ID
					lastCandidate = -1;
				}
			}
			// Final ID
			if (array.size() > 0)
				failCountData.add(new FailCountData(lastId, array.toArray()));

			IJ.showMessage(TITLE, "Loaded " + TextUtils.pleural(failCountData.size(), "sequence"));
			FailCountManager.failCountData = failCountData;
		}
		catch (NumberFormatException e)
		{
			IJ.error(TITLE, "Failed to load data:\n" + e.getMessage());
		}
		catch (IOException e)
		{
			IJ.error(TITLE, "Failed to load data:\n" + e.getMessage());
		}
		finally
		{
			if (br != null)
				try
				{
					br.close();
				}
				catch (IOException e)
				{
				}
		}
	}

	/**
	 * Guess the pass/fail status.
	 *
	 * @param string
	 *            the string
	 * @return true, if successful
	 * @throws IOException
	 *             Signals that an I/O exception has occurred.
	 */
	private boolean guessStatus(String string) throws IOException
	{
		int len = string.length();
		if (len < 1)
			return false;

		if (len == 1)
		{
			char c = string.charAt(0);
			if (c == 'y' || c == '1')
				return true;
			if (c == 'n' || c == '0')
				return false;
		}
		else
		{
			string = string.toLowerCase();
			if (string.equals("pass"))
				return true;
			if (string.equals("fail"))
				return false;
			if (string.equals("ok"))
				return true;
		}
		throw new IOException("Unrecognised status: " + string);
	}

	/**
	 * Save the data in memory to file.
	 */
	private void saveData()
	{
		TurboList<FailCountData> failCountData = FailCountManager.failCountData;
		if (failCountData.isEmpty())
		{
			IJ.error(TITLE, "No fail count data in memory");
			return;
		}
		String filename = Utils.getFilename("Fail_count_data_filename", settings.getFilename());
		if (filename == null)
			return;
		settings.setFilename(filename);

		BufferedWriter bw = null;
		try
		{
			bw = new BufferedWriter(new FileWriter(filename));
			bw.write("ID,Candidate,Status");
			bw.newLine();
			for (int i = 0; i < failCountData.size(); i++)
			{
				FailCountData d = failCountData.get(i);
				String prefix = d.id + ",";
				boolean[] pass = d.results;
				for (int j = 0; j < pass.length; j++)
				{
					bw.write(prefix);
					bw.write(Integer.toString(j + 1));
					bw.write(',');
					if (pass[j])
						bw.write('y');
					else
						bw.write('n');
					bw.newLine();
				}
			}
		}
		catch (IOException e)
		{
			IJ.error(TITLE, "Failed to save data:\n" + e.getMessage());
		}
		finally
		{
			if (bw != null)
				try
				{
					bw.close();
				}
				catch (IOException e)
				{
				}
		}
	}

	private static class PlotData
	{
		final int item;
		final int rollingWindow;
		final int passWeight;
		final int failWeight;
		final double resetFraction;
		final boolean fixedXAxis;

		PlotData(int item, boolean fixedXAxis, int rollingWindow, int passWeight, int failWeight, double resetFraction)
		{
			this.item = item;
			this.fixedXAxis = fixedXAxis;
			this.rollingWindow = rollingWindow;
			this.passWeight = passWeight;
			this.failWeight = failWeight;
			this.resetFraction = resetFraction;
		}

		/**
		 * Test if this equals the other object.
		 *
		 * @param that
		 *            the other object
		 * @return true, if successful
		 */
		public boolean equals(PlotData that)
		{
			//@formatter:off
			return that != null && 
					this.item == that.item && 
					this.fixedXAxis == that.fixedXAxis && 
					this.rollingWindow == that.rollingWindow &&
					this.passWeight == that.passWeight &&
					this.failWeight == that.failWeight &&
					this.resetFraction == that.resetFraction
					;
			//@formatter:on
		}

		/**
		 * Checks if is a new item.
		 *
		 * @param that
		 *            the other object
		 * @return true, if is new item
		 */
		public boolean isNewItem(PlotData that)
		{
			if (that == null)
				return true;
			return this.item != that.item;
		}
	}

	private class PlotWorker implements Runnable
	{
		final ConcurrentMonoStack<PlotData> stack;
		final TurboList<FailCountData> failCountData;
		PlotData lastPlotData = null;
		int maxSize = 0;

		PlotWorker(ConcurrentMonoStack<PlotData> stack, TurboList<FailCountData> failCountData)
		{
			this.stack = stack;
			this.failCountData = failCountData;
			for (int i = 0; i < failCountData.size(); i++)
			{
				maxSize = Math.max(maxSize, failCountData.getf(i).results.length);
			}
		}

		@Override
		public void run()
		{
			//while (!Thread.interrupted())
			while (true)
			{
				try
				{
					PlotData plotData = stack.pop();
					if (plotData == null)
						break;
					if (plotData.equals(lastPlotData))
						continue;
					run(plotData);
					lastPlotData = plotData;
				}
				catch (InterruptedException e)
				{
					//Thread.currentThread().interrupt();
					break;
				}
			}
		}

		private void run(PlotData plotData)
		{
			boolean isNew = plotData.isNewItem(lastPlotData) || plotData.fixedXAxis != lastPlotData.fixedXAxis;
			int item = plotData.item - 1; // 0-based index
			if (item < 0 || item >= failCountData.size())
				return;
			FailCountData data = failCountData.get(item);

			data.createData();
			WindowOrganiser wo = new WindowOrganiser();
			if (isNew)
			{
				display(wo, "Pass Count", data.candidate, data.passCount, plotData.fixedXAxis);
				display(wo, "Pass Rate", data.candidate, data.passRate, plotData.fixedXAxis);
				display(wo, "Consecutive Fail Count", data.candidate, data.consFailCount, plotData.fixedXAxis);
			}

			// Only rebuild if changed
			if (isNew || plotData.rollingWindow != lastPlotData.rollingWindow)
			{
				display(wo, "Rolling Fail Count", data.candidate, data.getRollingFailCount(plotData.rollingWindow),
						plotData.fixedXAxis);
			}
			if (isNew || plotData.passWeight != lastPlotData.passWeight ||
					plotData.failWeight != lastPlotData.failWeight)
			{
				display(wo, "Weighted Fail Count", data.candidate,
						data.getWeightedFailCount(plotData.passWeight, plotData.failWeight), plotData.fixedXAxis);
			}
			if (isNew || plotData.resetFraction != lastPlotData.resetFraction)
			{
				display(wo, "Resetting Fail Count", data.candidate, data.getResettingFailCount(plotData.resetFraction),
						plotData.fixedXAxis);
			}

			wo.tile();
		}

		private void display(WindowOrganiser wo, String string, float[] x, float[] y, boolean fixedXAxis)
		{
			String title = TITLE + " " + string;
			Plot plot = new Plot(title, "Candidate", string);
			double maxx = (fixedXAxis) ? maxSize : x[x.length - 1];
			double maxy = Maths.max(y);
			plot.setLimits(x[0], maxx, 0, maxy * 1.05);
			plot.addPoints(x, y, Plot.LINE);
			plot.addLabel(0, 0, "Max = " + maxy);
			Utils.display(title, plot, 0, wo);
		}
	}

	/**
	 * Show an interactive plot of the fail count data.
	 */
	private void plotData()
	{
		TurboList<FailCountData> failCountData = FailCountManager.failCountData;
		if (failCountData.isEmpty())
		{
			IJ.error(TITLE, "No fail count data in memory");
			return;
		}

		// Find max fail count size
		final int max = getMaxConsecutiveFailCount(failCountData);

		final ConcurrentMonoStack<PlotData> stack = new ConcurrentMonoStack<PlotData>();
		new Thread(new PlotWorker(stack, failCountData)).start();

		NonBlockingExtendedGenericDialog gd = new NonBlockingExtendedGenericDialog(TITLE);
		gd.addSlider("Item", 1, failCountData.size(), settings.getPlotItem());
		gd.addCheckbox("Fixed_x_axis", settings.getPlotFixedXAxis());
		gd.addMessage("Rolling Window Fail Count");
		gd.addSlider("Rolling_window", 1, 3 * max, settings.getPlotRollingWindow());
		gd.addMessage("Weighted Fail Count");
		gd.addSlider("Pass_weight", 1, 20, settings.getPlotPassWeight());
		gd.addSlider("Fail_weight", 1, 20, settings.getPlotFailWeight());
		gd.addMessage("Resetting Fail Count");
		gd.addSlider("Reset_fraction", 0.05, 0.95, settings.getPlotResetFraction());
		gd.addDialogListener(new DialogListener()
		{
			@Override
			public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
			{
				int item = (int) gd.getNextNumber();
				boolean fixedXAxis = gd.getNextBoolean();
				int rollingWindow = (int) gd.getNextNumber();
				int passWeight = (int) gd.getNextNumber();
				int failWeight = (int) gd.getNextNumber();
				double resetFraction = gd.getNextNumber();
				settings.setPlotItem(item);
				settings.setPlotRollingWindow(rollingWindow);
				settings.setPlotPassWeight(passWeight);
				settings.setPlotFailWeight(failWeight);
				settings.setPlotResetFraction(resetFraction);
				stack.insert(new PlotData(item, fixedXAxis, rollingWindow, passWeight, failWeight, resetFraction));
				return true;
			}
		});

		gd.hideCancelButton();
		gd.setOKLabel("Close");
		gd.showDialog();

		stack.close(gd.wasCanceled());
	}

	private int getMaxConsecutiveFailCount(TurboList<FailCountData> failCountData)
	{
		int max = 1;
		for (int i = 0; i < failCountData.size(); i++)
		{
			max = Math.max(max, failCountData.getf(i).getMaxConsecutiveFailCount());
		}
		return max;
	}

	private int getMaxFailCount(TurboList<FailCountData> failCountData)
	{
		int max = 1;
		for (int i = 0; i < failCountData.size(); i++)
		{
			max = Math.max(max, failCountData.getf(i).getFailCount());
		}
		return max;
	}

	@SuppressWarnings("unused")
	private int getMaxPassCount(TurboList<FailCountData> failCountData)
	{
		int max = 1;
		for (int i = 0; i < failCountData.size(); i++)
		{
			max = Math.max(max, failCountData.getf(i).getPassCount());
		}
		return max;
	}

	private void analyseData()
	{
		TurboList<FailCountData> failCountData = FailCountManager.failCountData;
		if (failCountData.isEmpty())
		{
			IJ.error(TITLE, "No fail count data in memory");
			return;
		}

		if (!showAnalysisDialog())
			return;

		final int maxCons = getMaxConsecutiveFailCount(failCountData);
		final int maxFail = getMaxFailCount(failCountData);
		//final int maxPass = getMaxPassCount(failCountData);

		// Create a set of fail counters
		final TurboList<FailCounter> counters = new TurboList<FailCounter>();
		TByteArrayList type = new TByteArrayList();
		for (int i = 0; i <= maxCons; i++)
		{
			counters.add(ConsecutiveFailCounter.create(i));
		}
		fill(type, counters, 0);

		// The other counters are user configured.
		// Ideally this would be a search to optimise the best parameters
		// for each counter as any enumeration may be way off the mark.

		// Note that 0 failures in a window can be scored using the consecutive fail counter.
		int max = Math.min(maxFail, settings.getRollingCounterMaxAllowedFailures());
		for (int fail = Maths.min(maxFail,
				Math.max(1, settings.getRollingCounterMinAllowedFailures())); fail <= max; fail++)
		{
			// Note that n-1 failures in window n can be scored using the consecutive fail counter.
			for (int window = Math.max(fail + 2, settings.getRollingCounterMinWindow()); window <= settings
					.getRollingCounterMaxWindow(); window++)
			{
				counters.add(RollingWindowFailCounter.create(fail, window));
			}
			switch (checkCounters(counters))
			{
				case ANALYSE:
					break;
				case CONTINUE:
					break;
				case RETURN:
					return;
				default:
					throw new IllegalStateException();
			}
		}
		fill(type, counters, 1);

		max = Math.min(maxFail, settings.getWeightedCounterMaxAllowedFailures());
		for (int fail = Maths.min(maxFail, settings.getWeightedCounterMinAllowedFailures()); fail <= max; fail++)
		{
			for (int w = settings.getWeightedCounterMinPassDecrement(); w <= settings
					.getWeightedCounterMaxPassDecrement(); w++)
			{
				counters.add(WeightedFailCounter.create(fail, 1, w));
			}
			switch (checkCounters(counters))
			{
				case ANALYSE:
					break;
				case CONTINUE:
					break;
				case RETURN:
					return;
				default:
					throw new IllegalStateException();
			}
		}
		fill(type, counters, 2);

		max = Math.min(maxFail, settings.getResettingCounterMaxAllowedFailures());
		for (int fail = Maths.min(maxFail, settings.getResettingCounterMinAllowedFailures()); fail <= max; fail++)
		{
			for (double f = settings.getResettingCounterMinResetFraction(); f <= settings
					.getResettingCounterMaxResetFraction(); f += settings.getResettingCounterIncResetFraction())
			{
				counters.add(ResettingFailCounter.create(fail, f));
			}
			switch (checkCounters(counters))
			{
				case ANALYSE:
					break;
				case CONTINUE:
					break;
				case RETURN:
					return;
				default:
					throw new IllegalStateException();
			}
		}
		fill(type, counters, 3);

		for (int count = settings.getPassRateCounterMinAllowedCounts(); count <= settings
				.getPassRateCounterMaxAllowedCounts(); count++)
		{
			for (double f = settings.getPassRateCounterMinPassRate(); f <= settings
					.getPassRateCounterMaxPassRate(); f += settings.getPassRateCounterIncPassRate())
			{
				counters.add(PassRateFailCounter.create(count, f));
			}
			switch (checkCounters(counters))
			{
				case ANALYSE:
					break;
				case CONTINUE:
					break;
				case RETURN:
					return;
				default:
					throw new IllegalStateException();
			}
		}
		fill(type, counters, 4);

		counters.trimToSize();

		// Score each of a set of standard fail counters against each frame using how 
		// close they are to the target.
		final double[] score = new double[counters.size()];
		final double targetPassFraction = settings.getTargetPassFraction();

		int nThreads = Prefs.getThreads();
		ExecutorService executor = Executors.newFixedThreadPool(nThreads);
		TurboList<Future<?>> futures = new TurboList<Future<?>>(nThreads);

		final Ticker ticker = Ticker.createStarted(new IJTrackProgress(), failCountData.size(), nThreads > 1);
		IJ.showStatus("Analysing " + TextUtils.pleural(counters.size(), "counter"));
		for (int i = 0; i < failCountData.size(); i++)
		{
			final FailCountData data = failCountData.getf(i);
			futures.add(executor.submit(new Runnable()
			{
				@Override
				public void run()
				{
					if (IJ.escapePressed())
						return;

					// TODO - Ideally this plugin should be run on benchmark data with ground truth.
					// The target could be to ensure all all the correct results are fit 
					// and false positives are excluded from incrementing the pass counter.
					// This could be done by saving the results from a benchmarking scoring
					// plugin to memory as the current dataset.
					data.initialiseAnalysis(targetPassFraction);

					// Score in blocks and then do a synchronized write to the combined score
					Thread t = Thread.currentThread();
					double[] s = new double[8192];
					int i = 0;
					while (i < counters.size())
					{
						if (t.isInterrupted())
							break;
						int block = Math.min(8192, counters.size() - i);
						for (int j = 0; j < block; j++)
						{
							FailCounter counter = counters.getf(i + j).newCounter();
							s[j] = data.score(counter);
						}
						// Write to the combined score
						synchronized (score)
						{
							for (int j = 0; j < block; j++)
							{
								score[i + j] += s[j];
							}
						}
						i += block;
					}
					ticker.tick();
				}
			}));
		}

		Utils.waitForCompletion(futures);
		executor.shutdown();
		IJ.showProgress(1);
		if (IJ.escapePressed())
		{
			IJ.showStatus("");
			IJ.error(TITLE, "Cancelled analysis");
			return;
		}
		IJ.showStatus("Summarising results ...");

		// TODO - check if the top filter is at the bounds of the range
		int minIndex = SimpleArrayUtils.findMinIndex(score);
		Utils.log(TITLE + " Analysis : Best counter = %s (Score = %f)", counters.getf(minIndex).getDescription(),
				score[minIndex]);

		// Show a table of results for the top N for each type
		int topN = Math.min(settings.getTableTopN(), score.length);
		if (topN > 0)
		{
			byte[] types = type.toArray();
			byte maxType = types[types.length - 1];
			createTable();
			for (byte b = 0; b <= maxType; b++)
			{
				int[] indices;
				// Use a heap to avoid a full sort
				IntResultHeap heap = new IntResultHeap(topN);
				for (int i = 0; i < score.length; i++)
					if (types[i] == b)
						heap.addValue(score[i], i);
				if (heap.getSize() == 0)
					continue;
				indices = heap.getData();
				// Ensure sorted
				Sort.sortAscending(indices, score);

				StringBuilder sb = new StringBuilder();
				BufferedTextWindow tw = new BufferedTextWindow(resultsWindow);
				for (int i = 0; i < topN; i++)
				{
					sb.setLength(0);
					int j = indices[i];
					sb.append(i + 1).append('\t');
					sb.append(counters.getf(j).getDescription()).append('\t');
					sb.append(score[j]);
					tw.append(sb.toString());
				}
				tw.flush();
			}
		}

		// TODO - Save the best fail counter to the current fit configuration.

		IJ.showStatus("");
	}

	private void fill(TByteArrayList type, TurboList<FailCounter> counters, int b)
	{
		int n = counters.size() - type.size();
		Utils.log("Type %d = %d", b, n);
		type.fill(type.size(), counters.size(), (byte) b);
	}

	private enum CounterStatus
	{
		CONTINUE, ANALYSE, RETURN;
	}

	private static int maxCounters = 200000;

	private CounterStatus checkCounters(TurboList<FailCounter> counters)
	{
		if (counters.size() > maxCounters)
		{
			GenericDialog gd = new GenericDialog(TITLE);
			gd.addMessage("Too many counters to analyse: " + counters.size());
			gd.addNumericField("Max_counters", maxCounters, 0);
			gd.enableYesNoCancel(" Analyse ", " Continue ");
			gd.showDialog();
			if (gd.wasCanceled())
				return CounterStatus.RETURN;
			if (gd.wasOKed())
				return CounterStatus.ANALYSE;
			int newMaxCounters = (int) gd.getNextNumber();
			if (newMaxCounters <= maxCounters)
			{
				IJ.error(TITLE, "The max counters has not been increased, unable to continue");
				return CounterStatus.RETURN;
			}
			maxCounters = newMaxCounters;
		}
		return CounterStatus.CONTINUE;
	}

	private void createTable()
	{
		if (resultsWindow == null || !resultsWindow.isShowing())
		{
			resultsWindow = new TextWindow(TITLE + " Analysis Results", "Rank\tFail Counter\tScore", "", 600, 400);
		}
	}

	private boolean showAnalysisDialog()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addMessage(TextUtils.wrap("Analysis a set of fail counters on the current pass/fail data.", 80));
		gd.addSlider("Target_pass_fraction", 0.1, 1, settings.getTargetPassFraction());
		gd.addSliderIncludeDefault("Table_top_n", 0, 100, settings.getTableTopN());
		gd.addNumericField("Rolling_counter_min_allowed_failures", settings.getRollingCounterMinAllowedFailures(), 0);
		gd.addNumericField("Rolling_counter_max_allowed_failures", settings.getRollingCounterMaxAllowedFailures(), 0);
		gd.addNumericField("Rolling_counter_min_window", settings.getRollingCounterMinWindow(), 0);
		gd.addNumericField("Rolling_counter_max_window", settings.getRollingCounterMaxWindow(), 0);
		gd.addNumericField("Weighted_counter_min_allowed_failures", settings.getWeightedCounterMinAllowedFailures(), 0);
		gd.addNumericField("Weighted_counter_max_allowed_failures", settings.getWeightedCounterMaxAllowedFailures(), 0);
		gd.addNumericField("Weighted_counter_min_pass_decrement", settings.getWeightedCounterMinPassDecrement(), 0);
		gd.addNumericField("Weighted_counter_max_pass_decrement", settings.getWeightedCounterMaxPassDecrement(), 0);
		gd.addNumericField("Resetting_counter_min_allowed_failures", settings.getResettingCounterMinAllowedFailures(),
				0);
		gd.addNumericField("Resetting_counter_max_allowed_failures", settings.getResettingCounterMaxAllowedFailures(),
				0);
		gd.addNumericField("Resetting_counter_min_pass_decrement", settings.getResettingCounterMinResetFraction(), 2);
		gd.addNumericField("Resetting_counter_max_pass_decrement", settings.getResettingCounterMaxResetFraction(), 2);
		gd.addNumericField("Resetting_counter_inc_pass_decrement", settings.getResettingCounterIncResetFraction(), 2);
		gd.addNumericField("Pass_rate_counter_min_allowed_failures", settings.getPassRateCounterMinAllowedCounts(), 0);
		gd.addNumericField("Pass_rate_counter_max_allowed_failures", settings.getPassRateCounterMaxAllowedCounts(), 0);
		gd.addNumericField("Pass_rate_counter_min_pass_rate", settings.getPassRateCounterMinPassRate(), 3);
		gd.addNumericField("Pass_rate_counter_max_pass_rate", settings.getPassRateCounterMaxPassRate(), 3);
		gd.addNumericField("Pass_rate_counter_inc_pass_rate", settings.getPassRateCounterIncPassRate(), 3);
		gd.showDialog();
		if (gd.wasCanceled())
			return false;
		settings.setTargetPassFraction(gd.getNextNumber());
		settings.setTableTopN((int) gd.getNextNumber());
		settings.setRollingCounterMinAllowedFailures((int) gd.getNextNumber());
		settings.setRollingCounterMaxAllowedFailures((int) gd.getNextNumber());
		settings.setRollingCounterMinWindow((int) gd.getNextNumber());
		settings.setRollingCounterMaxWindow((int) gd.getNextNumber());
		settings.setWeightedCounterMinAllowedFailures((int) gd.getNextNumber());
		settings.setWeightedCounterMaxAllowedFailures((int) gd.getNextNumber());
		settings.setWeightedCounterMinPassDecrement((int) gd.getNextNumber());
		settings.setWeightedCounterMaxPassDecrement((int) gd.getNextNumber());
		settings.setResettingCounterMinAllowedFailures((int) gd.getNextNumber());
		settings.setResettingCounterMaxAllowedFailures((int) gd.getNextNumber());
		settings.setResettingCounterMinResetFraction(gd.getNextNumber());
		settings.setResettingCounterMaxResetFraction(gd.getNextNumber());
		settings.setResettingCounterIncResetFraction(gd.getNextNumber());
		settings.setPassRateCounterMinAllowedCounts((int) gd.getNextNumber());
		settings.setPassRateCounterMaxAllowedCounts((int) gd.getNextNumber());
		settings.setPassRateCounterMinPassRate(gd.getNextNumber());
		settings.setPassRateCounterMaxPassRate(gd.getNextNumber());
		settings.setPassRateCounterIncPassRate(gd.getNextNumber());
		try
		{
			Parameters.isAboveZero("Target pass fraction", settings.getTargetPassFraction());
			Parameters.isAboveZero("Resetting counter inc pass decrement",
					settings.getResettingCounterIncResetFraction());
			Parameters.isAboveZero("Pass rate counter inc pass rate", settings.getPassRateCounterIncPassRate());
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}
		return true;
	}
}
