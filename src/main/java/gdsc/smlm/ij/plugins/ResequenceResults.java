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

import gdsc.core.ij.IJTrackProgress;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.procedures.PeakResultProcedure;
import gdsc.core.logging.TrackProgress;
import ij.IJ;
import ij.gui.ExtendedGenericDialog;
import ij.plugin.PlugIn;

/**
 * Updates the frame numbers on results that are stored in memory.
 */
public class ResequenceResults implements PlugIn
{
	private static final String TITLE = "Resequence Results";
	private static String inputOption = "";
	private static int start = 1;
	private static int block = 1;
	private static int skip = 0;
	private static boolean logMapping = false;

	/*
	 * (non-)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (MemoryPeakResults.isMemoryEmpty())
		{
			IJ.error(TITLE, "There are no fitting results in memory");
			return;
		}

		if (!showDialog())
			return;

		MemoryPeakResults results = ResultsManager.loadInputResults(inputOption, true);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			return;
		}

		if (resequenceResults(results, start, block, skip, (logMapping) ? new IJTrackProgress() : null))
			IJ.showStatus("Resequenced " + results.getName());
	}

	private boolean showDialog()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Resequence the results in memory (assumed to be continuous from 1).\n" +
				"Describe the regular repeat of the original image:\n" +
				"Start = The first frame that contained the data\n" +
				"Block = The number of continuous frames containing data\n" +
				"Skip = The number of continuous frames to ignore before the next data\n \n" +
				"E.G. 2:9:1 = Data was imaged from frame 2 for 9 frames, 1 frame to ignore, then repeat.");

		ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);
		gd.addNumericField("Start", start, 0);
		gd.addNumericField("Block", block, 0);
		gd.addNumericField("Skip", skip, 0);
		gd.addCheckbox("Log_mapping", logMapping);

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		inputOption = ResultsManager.getInputSource(gd);
		start = (int) gd.getNextNumber();
		block = (int) gd.getNextNumber();
		skip = (int) gd.getNextNumber();
		logMapping = gd.getNextBoolean();

		// Check arguments
		try
		{
			Parameters.isAboveZero("Start", start);
			Parameters.isAboveZero("Block", block);
			Parameters.isPositive("Skip", skip);
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		return true;
	}

	private static class ResequencePeakResultProcedure implements PeakResultProcedure
	{
		int start;
		TrackProgress tracker;

		ResequencePeakResultProcedure(int start, TrackProgress tracker)
		{
			this.start = start;
			this.tracker= tracker;
		}
		
		public void execute(PeakResult r)
		{
			int t = 1; // The current frame in the results
			int mapped = start; // The mapped frame in the results
			int b = 1; // The current block size
			
			boolean print = (tracker != null);
			
			if (t != r.getFrame())
			{
				// Update the mapped position
				while (t < r.getFrame())
				{
					// Move to the next position
					mapped++;

					// Check if this move will make the current block too large 
					if (++b > block)
					{
						// Skip 
						mapped += skip;
						b = 1;
					}

					t++;
				}

				t = r.getFrame();
				print = (tracker != null);
			}

			r.setFrame(mapped);

			if (print)
			{
				print = false;
				tracker.log("Map %d -> %d", t, mapped);
			}
		}
	}
	
	/**
	 * Resequence the results for the original imaging sequence provided. Results are assumed to be continuous from 1.
	 * 
	 * @param results
	 * @param start
	 *            The first frame that contained the data
	 * @param block
	 *            The number of continuous frames containing data
	 * @param skip
	 *            The number of continuous frames to ignore before the next data
	 * @param tracker
	 *            Used to report the mapping
	 * @return
	 */
	private static boolean resequenceResults(MemoryPeakResults results, final int start, final int block,
			final int skip, final TrackProgress tracker)
	{
		if (results == null || results.size() == 0)
			return false;

		results.sort();

		// Assume the results start from frame 1 (or above)
		if (results.getFirstFrame() < 1)
		{
			return false;
		}

		results.forEach(new ResequencePeakResultProcedure(start, tracker));

		return true;
	}
}
