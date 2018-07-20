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

import com.google.protobuf.Message;

import ij.IJ;
import ij.Prefs;
import ij.gui.ExtendedGenericDialog;
import ij.plugin.PlugIn;
import uk.ac.sussex.gdsc.core.ij.Utils;
import uk.ac.sussex.gdsc.smlm.ij.settings.Constants;
import uk.ac.sussex.gdsc.smlm.results.PeakResultsReader;
import uk.ac.sussex.gdsc.smlm.utils.XmlUtils;

/**
 * This plugin allows the header to be displayed from a PeakFit results file.
 */
public class ShowResultsHeader implements PlugIn
{
	private static String TITLE = "Show Results Header";

	private static String inputFilename = Prefs.get(Constants.inputFilename, "");
	private static boolean raw = false;

	/*
	 * (non-Javadoc)
	 *
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	@Override
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addMessage("Show the results header in the ImageJ log");
		gd.addFilenameField("Filename", inputFilename, 30);
		gd.addCheckbox("Raw", raw);

		gd.showDialog();
		if (gd.wasCanceled())
			return;

		inputFilename = gd.getNextString();
		raw = gd.getNextBoolean();

		Prefs.set(Constants.inputFilename, inputFilename);

		final PeakResultsReader reader = new PeakResultsReader(inputFilename);
		final String header = reader.getHeader();
		if (header == null)
		{
			IJ.error(TITLE, "No header found in file: " + inputFilename);
			return;
		}
		if (raw)
		{
			// The ImageJ TextPanel class correctly stores lines with tab characters.
			// However when it is drawn in ij.text.TextCanvas using java.awt.Graphics.drawChars(...)
			// the instance of this class is sun.java2d.SunGraphics2D which omits '\t' chars.
			// This may be a problem specific to the Linux JRE.
			// TODO - Find out if this is a Linux specific bug.

			// Output the raw text. This preserves the tabs in the Cut/Copy commands.
			IJ.log(header);
			// Replace tabs by 4 spaces:
			//IJ.log(header.replace("\t", "    "));
			return;
		}
		// Output what information we can extract
		boolean found = false;
		found |= show("Format", reader.getFormat().toString());
		found |= show("Name", reader.getName());
		found |= show("Bounds", reader.getBounds());
		found |= show("Source", reader.getSource());
		found |= show("Calibration", reader.getCalibration());
		found |= show("PSF", reader.getPSF());
		found |= show("Configuration", reader.getConfiguration());
		if (!found)
			IJ.error(TITLE, "No header information found in file: " + inputFilename);
	}

	private static boolean show(String title, Object data)
	{
		if (data == null)
			return false;
		String text = (data instanceof String) ? (String) data : XmlUtils.toXML(data);
		if (text.startsWith("{"))
			text = uk.ac.sussex.gdsc.smlm.utils.JSONUtils.simplify(text);
		else if (text.startsWith("<"))
			text = uk.ac.sussex.gdsc.core.utils.XmlUtils.prettyPrintXml(text);
		Utils.log("%s: %s", title, text);
		return true;
	}

	private static boolean show(String title, Message data)
	{
		if (data == null)
			return false;
		Utils.log("%s:\n%s", title, data.toString());
		return true;
	}
}
