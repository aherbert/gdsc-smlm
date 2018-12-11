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

import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;

import ij.IJ;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;

/**
 * Allows results held in memory to be renamed.
 */
public class RenameResults implements PlugIn {
  private static final String TITLE = "Rename Results";

  private String renameText = "";

  /*
   * (non-)
   *
   * @see ij.plugin.PlugIn#run(java.lang.String)
   */
  @Override
  public void run(String arg) {
    SMLMUsageTracker.recordPlugin(this.getClass(), arg);

    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "There are no fitting results in memory");
      return;
    }

    if (!showDialog()) {
      return;
    }

    IJ.showStatus("Renamed " + renameResults() + " result sets");
  }

  private boolean showDialog() {
    final GenericDialog gd = new GenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);

    gd.addMessage("To rename the results in memory update the second name field as desired.\n"
        + "(Note the semi-colon at the end of the line is needed for macro recording.)");

    final StringBuilder sb = new StringBuilder();
    for (final String name : MemoryPeakResults.getResultNames()) {
      sb.append(name).append(" = ").append(name).append(";\n");
    }

    gd.addTextAreas(sb.toString(), null, 20, 80);

    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }

    renameText = gd.getNextText();

    return true;
  }

  private int renameResults() {
    // Start with the original names for the mapping from old to new
    final HashMap<String, String> mappedNames = new HashMap<>();
    for (final String name : MemoryPeakResults.getResultNames()) {
      mappedNames.put(name, name);
    }

    // Get the new names
    final String[] lines = renameText.split("[;\n]");
    for (final String line : lines) {
      final String[] fields = line.split("[=]");
      if (fields.length == 2) {
        final String oldName = fields[0].trim();
        final String newName = fields[1].trim();

        if (!mappedNames.containsKey(oldName)) {
          IJ.error(TITLE, "An unknown original name has been specified: " + oldName);
          return 0;
        }

        if (oldName.equals(newName)) {
          // No update required
          continue;
        }

        mappedNames.put(oldName, newName);
      }
    }

    // Check the new names are unique
    final Set<String> newNames = new HashSet<>();
    for (final String newName : mappedNames.values()) {
      if (newNames.contains(newName)) {
        IJ.error(TITLE, "A duplicate new name has been specified: " + newName);
        return 0;
      }
      newNames.add(newName);
    }

    // Rename
    final List<MemoryPeakResults> renamedResults = new LinkedList<>();
    for (final Entry<String, String> entry : mappedNames.entrySet()) {
      if (entry.getKey().equals(entry.getValue())) {
        continue;
      }
      // Remove from memory and store in a list
      final MemoryPeakResults results = MemoryPeakResults.removeResults(entry.getKey());
      if (results != null) {
        results.setName(entry.getValue());
        renamedResults.add(results);
      }
    }

    // Add back to memory
    for (final MemoryPeakResults results : renamedResults) {
      MemoryPeakResults.addResults(results);
    }

    return renamedResults.size();
  }
}
