/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

import uk.ac.sussex.gdsc.core.ij.ImageJTrackProgress;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.TextFilePeakResults;
import uk.ac.sussex.gdsc.smlm.results.Trace;
import uk.ac.sussex.gdsc.smlm.results.TraceManager;

import ij.IJ;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

/**
 * Run a tracing algorithm on the peak results to trace neighbours across the frames.
 */
public class NeighbourAnalysis implements PlugIn {
  private static final String TITLE = "Neighbour Analysis";
  private static String inputOption = "";
  private static double distanceThreshold = 0.6;
  private static int timeThreshold = 1;

  private static String filename = "";

  private MemoryPeakResults results;

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "No localisations in memory");
      return;
    }

    if (!showDialog()) {
      return;
    }

    final TraceManager manager = new TraceManager(results);

    // Run the tracing
    manager.setTracker(new ImageJTrackProgress());
    final Trace[] traces = manager.findNeighbours(distanceThreshold, timeThreshold);

    saveTraces(traces);
  }

  private void saveTraces(Trace[] traces) {
    final String[] path = ImageJUtils.decodePath(filename);
    final OpenDialog chooser = new OpenDialog("Traces_File", path[0], path[1]);
    if (chooser.getFileName() != null) {
      filename = chooser.getDirectory() + chooser.getFileName();

      // Remove extension and replace with .xls
      final int index = filename.lastIndexOf('.');
      if (index > 0) {
        filename = filename.substring(0, index);
      }
      filename += ".xls";

      final boolean showDeviations = results.hasDeviations();
      final TextFilePeakResults traceResults = new TextFilePeakResults(filename, showDeviations);
      traceResults.copySettings(results);
      traceResults.begin();
      if (!traceResults.isActive()) {
        IJ.error(TITLE, "Failed to write to file: " + filename);
        return;
      }
      traceResults.addComment(createSettingsComment());
      for (final Trace trace : traces) {
        traceResults.addCluster(trace); // addTrace(...) does a sort on the results
      }
      traceResults.end();
    }
  }

  private static String createSettingsComment() {
    return String.format("Neighbour tracing : distance-threshold = %f : time-threshold = %d",
        distanceThreshold, timeThreshold);
  }

  private boolean showDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);

    ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);

    gd.addNumericField("Distance_Threshold (px)", distanceThreshold, 4);
    gd.addNumericField("Time_Threshold (frames)", timeThreshold, 0);

    gd.showDialog();

    if (gd.wasCanceled() || !readDialog(gd)) {
      return false;
    }

    // Load the results
    results = ResultsManager.loadInputResults(inputOption, false, DistanceUnit.PIXEL);
    if (results == null || results.size() == 0) {
      IJ.error(TITLE, "No results could be loaded");
      IJ.showStatus("");
      return false;
    }

    return true;
  }

  private static boolean readDialog(ExtendedGenericDialog gd) {
    inputOption = ResultsManager.getInputSource(gd);
    distanceThreshold = gd.getNextNumber();
    timeThreshold = (int) gd.getNextNumber();

    if (distanceThreshold < 0) {
      distanceThreshold = 0;
    }
    if (timeThreshold < 0) {
      timeThreshold = 0;
    }
    if (timeThreshold == 0 && distanceThreshold == 0) {
      IJ.error(TITLE, "No thresholds specified");
      return false;
    }
    return !gd.invalidNumber();
  }
}
