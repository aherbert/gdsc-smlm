package gdsc.smlm.ij.plugins;

import gdsc.core.ij.Utils;
import gdsc.smlm.ij.settings.Constants;
import gdsc.smlm.results.PeakResultsReader;
import gdsc.smlm.utils.XmlUtils;
import ij.IJ;
import ij.Prefs;
import ij.gui.ExtendedGenericDialog;
import ij.plugin.PlugIn;

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
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addMessage("Show the results header in the ImageJ log");
		gd.addFilenameField("Filename", inputFilename, 30);
		gd.addCheckbox("Raw", raw);

		gd.showDialog();
		if (gd.wasCanceled())
			return;

		inputFilename = gd.getNextString();
		raw = gd.getNextBoolean();

		Prefs.set(Constants.inputFilename, inputFilename);

		PeakResultsReader reader = new PeakResultsReader(inputFilename);
		String header = reader.getHeader();
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
		found |= show("Calibration", reader.getCalibration());
		found |= show("Configuration", reader.getConfiguration());
		if (!found)
			IJ.error(TITLE, "No header information found in file: " + inputFilename);
	}

	private boolean show(String title, Object data)
	{
		if (data == null)
			return false;
		String text = (data instanceof String) ? (String) data : XmlUtils.toXML(data);
		if (text.startsWith("<"))
			text = XmlUtils.prettyPrintXml(text);
		Utils.log("%s: %s", title, text);
		return true;
	}
}
