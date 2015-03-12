package ij.gui;

import ij.measure.Measurements;
import ij.plugin.filter.Analyzer;

/**
 * Extension of the ij.gui.PlotWindow class to add functionality
 */
public class PlotWindow2 extends PlotWindow
{
	private static final long serialVersionUID = 5935603633626914545L;

	private static int precision = Analyzer.getPrecision();
	private static boolean scientific;
	static
	{
		int measurements = Analyzer.getMeasurements();
		scientific = (measurements & Measurements.SCIENTIFIC_NOTATION) != 0;
	}
	private static boolean update = false;
	private static long time = 0;
	private static boolean lock = false;

	/**
	 * Construct a plot window.
	 * This method throws an exception on some platforms since the super constructor is package private.
	 * @param plot
	 */
	PlotWindow2(Plot plot)
	{
		super(plot);
	}

	@Deprecated
	public PlotWindow2(String title, String xLabel, String yLabel, double[] xValues, double[] yValues)
	{
		super(title, xLabel, yLabel, xValues, yValues);
	}

	@Deprecated
	public PlotWindow2(String title, String xLabel, String yLabel, float[] xValues, float[] yValues)
	{
		super(title, xLabel, yLabel, xValues, yValues);
	}

	private boolean askForPrecision;

	@Override
	public void saveAsText()
	{
		requireNewPrecision();
		super.saveAsText();
	}

	@Override
	public void showList()
	{
		requireNewPrecision();
		super.showList();
	}

	@Override
	public void copyToClipboard()
	{
		requireNewPrecision();
		super.copyToClipboard();
	}

	/**
	 * Set a flag to show a dialog to ask for the precision if one second has elapsed since the last time
	 */
	private void requireNewPrecision()
	{
		askForPrecision = System.currentTimeMillis() > time + 1000;
	}

	/**
	 * Override the getPrecision function to allow the user to select the precision for the numbers
	 * 
	 * @param values
	 * @return
	 */
	@Override
	public int getPrecision(float[] values)
	{
		// Use a simple lock to ensure no two threads ask at the same time
		if (askForPrecision && lock())
		{
			GenericDialog gd = new GenericDialog("Plot precision");
			gd.addSlider("Plot_precision", 0, 9, precision);
			gd.addCheckbox("Scientific_notation", scientific);
			gd.addCheckbox("Update_preferences", update);
			gd.showDialog();
			if (!gd.wasCanceled())
			{
				int p = (int) gd.getNextNumber();
				scientific = gd.getNextBoolean();
				update = gd.getNextBoolean();
				if (!gd.invalidNumber())
				{
					precision = Math.max(0, Math.min(p, 9));

					if (update)
					{
						Analyzer.setPrecision(p);
						Analyzer.setMeasurement(Analyzer.SCIENTIFIC_NOTATION, scientific);
					}
				}
			}
			time = System.currentTimeMillis();
			askForPrecision = false;
			lock = false;
		}

		int setDigits = precision;
		boolean scientificNotation = scientific;
		int minDecimalPlaces = 4;
		if (scientificNotation)
		{
			if (setDigits < minDecimalPlaces)
				setDigits = minDecimalPlaces;
			return -setDigits;
		}
		int digits = minDecimalPlaces;
		if (setDigits > digits)
			digits = setDigits;
		boolean realValues = false;
		for (int i = 0; i < values.length; i++)
		{
			if ((int) values[i] != values[i])
			{
				realValues = true;
				break;
			}
		}
		if (!realValues)
			digits = 0;
		return digits;
	}

	private synchronized boolean lock()
	{
		if (!lock)
		{
			return lock = true;
		}
		return false;
	}
}
