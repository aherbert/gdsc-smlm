package gdsc.smlm.ij.plugins;

import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;

import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

/**
 * Filters PeakFit results that are stored in memory using various fit criteria.
 */
public class CropResults implements PlugIn
{
	private static final String TITLE = "Crop Results";
	private static String inputOption = "";
	private static double border = 0;
	private static double x = 0, y = 0, width = 0, height = 0;
	private static boolean selectRegion;

	private MemoryPeakResults results;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (MemoryPeakResults.countMemorySize() == 0)
		{
			IJ.error(TITLE, "There are no fitting results in memory");
			return;
		}

		// Show a dialog allowing the results set to be filtered
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addMessage("Select a dataset to crop");
		ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		inputOption = ResultsManager.getInputSource(gd);
		results = ResultsManager.loadInputResults(inputOption, false);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			IJ.showStatus("");
			return;
		}

		if (!showDialog())
			return;

		cropResults();
	}

	private boolean showDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		Rectangle bounds = results.getBounds(true);

		gd.addMessage(String.format("x=%d,y=%d,w=%d,h=%d", bounds.x, bounds.y, bounds.width, bounds.height));
		gd.addNumericField("Border", border, 2);
		gd.addCheckbox("Select_region", selectRegion);
		gd.addNumericField("X", x, 2);
		gd.addNumericField("Y", y, 2);
		gd.addNumericField("Width", width, 2);
		gd.addNumericField("Height", height, 2);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		border = Math.max(0, gd.getNextNumber());
		selectRegion = gd.getNextBoolean();
		x = gd.getNextNumber();
		y = gd.getNextNumber();
		width = Math.max(0, gd.getNextNumber());
		height = Math.max(0, gd.getNextNumber());

		return true;
	}

	/**
	 * Apply the filters to the data
	 */
	private void cropResults()
	{
		MemoryPeakResults newResults = new MemoryPeakResults();
		newResults.copySettings(results);
		newResults.setName(results.getName() + " Cropped");
		MemoryPeakResults.addResults(newResults);

		// These bounds are integer. But this is because the results are meant to come from an image.
		Rectangle bounds = results.getBounds(true);

		// The crop bounds can be floating point...
		
		// Border
		double xx = bounds.x + border;
		double yy = bounds.y + border;
		double w = Math.max(0, bounds.width - 2 * border);
		double h = Math.max(0, bounds.height - 2 * border);
		Rectangle2D borderBounds = new Rectangle2D.Double(xx, yy, w, h);

		// Bounding box
		if (selectRegion)
		{
			Rectangle2D boxBounds = new Rectangle2D.Double(x, y, width, height);
			borderBounds = borderBounds.createIntersection(boxBounds);
		}

		if (borderBounds.getWidth() > 0 && borderBounds.getHeight() > 0)
		{
			for (PeakResult result : results.getResults())
			{
				if (borderBounds.contains(result.getXPosition(), result.getYPosition()))
					newResults.add(result);
			}
		}

		IJ.showStatus(newResults.size() + " Cropped localisations");
	}
}
