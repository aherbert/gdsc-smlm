package ij.gui;

import ij.measure.Measurements;
import ij.plugin.filter.Analyzer;

/**
 * Extension of the ij.gui.PlotWindow class to add functionality
 */
public class SuperPlotWindow extends PlotWindow
{
	private static final long serialVersionUID = 5935603633626914545L;
	
	private static int precision = Analyzer.getPrecision();
	private static boolean scientific;
	static {
		int measurements = Analyzer.getMeasurements();
		scientific = (measurements & Measurements.SCIENTIFIC_NOTATION) != 0;
	}
	private static boolean update = false;
	private static long time = 0;
	private static boolean lock = false;

	public SuperPlotWindow(Plot plot)
	{
		super(plot);
	}

	@Deprecated
	public SuperPlotWindow(String title, String xLabel, String yLabel, double[] xValues, double[] yValues)
	{
		super(title, xLabel, yLabel, xValues, yValues);
	}

	@Deprecated
	public SuperPlotWindow(String title, String xLabel, String yLabel, float[] xValues, float[] yValues)
	{
		super(title, xLabel, yLabel, xValues, yValues);
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
		// Show a dialog to ask for the precision if one second has elapsed
		if (System.currentTimeMillis() > time + 1000 && lock())
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
