package gdsc.smlm.ij.plugins;

import java.awt.Rectangle;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.regex.Pattern;

import gdsc.core.ij.Utils;
import gdsc.core.utils.BooleanArray;
import gdsc.core.utils.TextUtils;
import gdsc.core.utils.TurboList;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

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
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.ExtendedGenericDialog;
import ij.plugin.PlugIn;

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
		CREATE_DATA { public String getName() { return "Create Data"; } },
		LOAD_DATA { public String getName() { return "Load Data"; } },
		SAVE_DATA { public String getName() { return "Save Data"; } },
		;

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

		public FailCountData(int id, boolean[] results)
		{
			this.id = id;
			// Find the last success
			int end = results.length - 1;
			while (end > 0 && !results[end])
				end--;
			this.results = Arrays.copyOf(results, end + 1);
		}
	}

	private static TurboList<FailCountData> failCountData = new TurboList<FailCountData>(1);

	private FailCountManagerSettings.Builder settings;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
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
		// XXX - Ensure this disables all advanced stopping criteria.  
		fitConfig.setFailuresLimit(settings.getFailCountLimit());

		ImageSource source = new IJImageSource(imp);
		PeakFit peakFit = new PeakFit(fitConfig, ResultsSettings.getDefaultInstance());
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
				if (Utils.showStatus("Slice: " + slice + " / " + totalFrames))
					IJ.showProgress(slice, totalFrames);
			}

			ParameterisedFitJob job = createJob(source.getStartFrameNumber(), data, bounds);
			jobs.addf(job);
			engine.run(job);

			shutdown = escapePressed();
		}

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
				failCountData.add(new FailCountData(job.getSlice(), fitParams.pass));
			}
		}
		FailCountManager.failCountData = failCountData;

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
}
