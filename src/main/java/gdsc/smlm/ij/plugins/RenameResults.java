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

import gdsc.smlm.results.MemoryPeakResults;
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
 * Filters PeakFit results that are stored in memory using the configured filters.
 */
public class RenameResults implements PlugIn
{
	private static final String TITLE = "Rename Results";

	private String renameText = "";

	/*
	 * (non-)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		if (MemoryPeakResults.countMemorySize() == 0)
		{
			IJ.error(TITLE, "There are no fitting results in memory");
			return;
		}

		if (!showDialog())
			return;

		IJ.showStatus("Renamed " + renameResults() + " result sets");
	}

	private boolean showDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("To rename the results in memory update the second name field as desired.\n"
				+ "(Note the semi-colon at the end of the line is needed for macro recording.)");

		StringBuilder sb = new StringBuilder();
		for (String name : MemoryPeakResults.getResultNames())
		{
			sb.append(name).append(" = ").append(name).append(";\n");
		}

		gd.addTextAreas(sb.toString(), null, 20, 80);

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		renameText = gd.getNextText();

		return true;
	}

	private int renameResults()
	{
		// Start with the original names for the mapping from old to new
		HashMap<String, String> mappedNames = new HashMap<String, String>();
		for (String name : MemoryPeakResults.getResultNames())
			mappedNames.put(name, name);
		
		// Get the new names
		String[] lines = renameText.split("\n");
		for (String line : lines)
		{
			String[] fields = line.split("[=;]");
			if (fields.length == 2)
			{
				String oldName = fields[0].trim();
				String newName = fields[1].trim();

				if (!mappedNames.containsKey(oldName))
				{
					IJ.error(TITLE, "An unknown original name has been specified: " + oldName);
					return 0;
				}

				if (oldName.equals(newName))
					// No update required
					continue;
				
				mappedNames.put(oldName, newName);
			}
		}

		// Check the new names are unique
		Set<String> newNames = new HashSet<String>();
		for (String newName : mappedNames.values())
		{
			if (newNames.contains(newName))
			{
				IJ.error(TITLE, "A duplicate new name has been specified: " + newName);
				return 0;
			}
			newNames.add(newName);
		}
		
		// Rename
		List<MemoryPeakResults> renamedResults = new LinkedList<MemoryPeakResults>();
		for (Entry<String, String> entry : mappedNames.entrySet())
		{
			if (entry.getKey().equals(entry.getValue()))
				continue;
			// Remove from memory and store in a list
			MemoryPeakResults results = MemoryPeakResults.removeResults(entry.getKey());
			if (results != null)
			{
				results.setName(entry.getValue());
				renamedResults.add(results);
			}
		}
		
		// Add back to memory
		for (MemoryPeakResults results : renamedResults)
			MemoryPeakResults.addResults(results);
		
		return renamedResults.size();
	}
}
