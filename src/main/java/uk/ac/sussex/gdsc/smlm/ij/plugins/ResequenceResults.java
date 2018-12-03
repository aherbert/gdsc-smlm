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

import ij.IJ;
import ij.plugin.PlugIn;
import uk.ac.sussex.gdsc.core.ij.ImageJTrackProgress;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.logging.TrackProgress;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;

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
    @Override
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

        final MemoryPeakResults results = ResultsManager.loadInputResults(inputOption, true, null, null);
        if (results == null || results.size() == 0)
        {
            IJ.error(TITLE, "No results could be loaded");
            return;
        }

        if (resequenceResults(results, start, block, skip, (logMapping) ? new ImageJTrackProgress() : null))
            IJ.showStatus("Resequenced " + results.getName());
    }

    private static boolean showDialog()
    {
        final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
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
        catch (final IllegalArgumentException e)
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
            this.tracker = tracker;
        }

        @Override
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
     *            the results
     * @param start
     *            The first frame that contained the data
     * @param block
     *            The number of continuous frames containing data
     * @param skip
     *            The number of continuous frames to ignore before the next data
     * @param tracker
     *            Used to report the mapping
     * @return true, if successful
     */
    private static boolean resequenceResults(MemoryPeakResults results, final int start, final int block,
            final int skip, final TrackProgress tracker)
    {
        if (results == null || results.size() == 0)
            return false;

        results.sort();

        // Assume the results start from frame 1 (or above)
        if (results.getFirstFrame() < 1)
            return false;

        results.forEach(new ResequencePeakResultProcedure(start, tracker));

        return true;
    }
}
