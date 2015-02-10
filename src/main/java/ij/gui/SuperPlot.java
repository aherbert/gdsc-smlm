package ij.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.macro.Interpreter;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import ij.util.Tools;

import java.lang.reflect.Field;

/**
 * Extension of the ij.gui.Plot class to add functionality
 */
public class SuperPlot extends Plot
{
	/** Draw a bar plot */
	public static final int BAR = 999;

	public SuperPlot(String title, String xLabel, String yLabel, float[] xValues, float[] yValues)
	{
		super(title, xLabel, yLabel, xValues, yValues);
	}

	public SuperPlot(String title, String xLabel, String yLabel, double[] xValues, double[] yValues)
	{
		super(title, xLabel, yLabel, xValues, yValues);
	}

	public SuperPlot(String dummy, String title, String xLabel, String yLabel, float[] xValues, float[] yValues)
	{
		super(title, xLabel, yLabel, xValues, yValues);
	}

	public SuperPlot(String title, String xLabel, String yLabel)
	{
		super(title, xLabel, yLabel, (float[]) null, (float[]) null);
	}

	public SuperPlot(String title, String xLabel, String yLabel, int flags)
	{
		super(title, xLabel, yLabel, (float[]) null, (float[]) null, flags);
	}

	public SuperPlot(String title, String xLabel, String yLabel, float[] xValues, float[] yValues, int flags)
	{
		super(title, xLabel, yLabel, xValues, yValues, flags);
	}

	public SuperPlot(String title, String xLabel, String yLabel, double[] xValues, double[] yValues, int flags)
	{
		super(title, xLabel, yLabel, xValues, yValues, flags);
	}

	/**
	 * Adds a set of points to the plot, adds a curve if shape is set to LINE or bars if the shape is BAR
	 * 
	 * @param x
	 *            the x-coodinates
	 * @param y
	 *            the y-coodinates
	 * @param shape
	 *            CIRCLE, X, BOX, TRIANGLE, CROSS, DOT, LINE or BAR
	 */
	@Override
	public void addPoints(float[] x, float[] y, int shape)
	{
		if (xValues == null)
		{
			// Set the limits if this is the first set of data. The limits are usually set in the constructor
			// but we may want to not pass in the values to the constructor and then immediately call 
			// addPoints(x, y, SuperPlot.BAR)
			double[] a = Tools.getMinMax(x);
			xMin = a[0];
			xMax = a[1];
			a = Tools.getMinMax(y);
			yMin = a[0];
			yMax = a[1];
		}
		
		// This only works if the addPoints super method ignores the BAR option but still store the values
		setup();
		switch (shape)
		{
			case BAR:
				ImageProcessor ip = getProcessor();
				int flags = 0;
				try
				{
					Field f = Plot.class.getDeclaredField("flags");
					f.setAccessible(true);
					flags = f.getInt(this);
				}
				catch (Exception e)
				{
					// Ignore
				}
				float[] xValues = createHistogramAxis(x);
				float[] yValues = createHistogramValues(y);
				drawFloatPolyline(ip, ((flags & X_LOG_NUMBERS) != 0) ? arrayToLog(xValues) : xValues,
						((flags & Y_LOG_NUMBERS) != 0) ? arrayToLog(yValues) : yValues, xValues.length);
				break;
		}

		super.addPoints(x, y, shape);
	}

	private float[] arrayToLog(float[] val)
	{
		float[] newVal = new float[val.length];
		for (int i = 0; i < newVal.length; i++)
			newVal[i] = (float) (Math.log10(val[i]));
		return newVal;
	}

	/**
	 * For the provided histogram x-axis bins, produce an x-axis for plotting. This functions doubles up the histogram
	 * x-positions to allow plotting a square line profile using the ImageJ plot command.
	 * 
	 * @param histogramX
	 * @return
	 */
	public static float[] createHistogramAxis(float[] histogramX)
	{
		float[] axis = new float[histogramX.length * 2 + 2];
		int index = 0;
		for (int i = 0; i < histogramX.length; ++i)
		{
			axis[index++] = histogramX[i];
			axis[index++] = histogramX[i];
		}
		if (histogramX.length > 0)
		{
			float dx = (histogramX.length == 1) ? 1 : (histogramX[1] - histogramX[0]);
			axis[index++] = histogramX[histogramX.length - 1] + dx;
			axis[index++] = histogramX[histogramX.length - 1] + dx;
		}
		return axis;
	}

	/**
	 * For the provided histogram y-axis values, produce a y-axis for plotting. This functions doubles up the histogram
	 * values to allow plotting a square line profile using the ImageJ plot command.
	 * 
	 * @param histogramY
	 * @return
	 */
	public static float[] createHistogramValues(float[] histogramY)
	{
		float[] axis = new float[histogramY.length * 2 + 2];

		int index = 0;
		axis[index++] = 0;
		for (int i = 0; i < histogramY.length; ++i)
		{
			axis[index++] = histogramY[i];
			axis[index++] = histogramY[i];
		}
		return axis;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.gui.Plot#show()
	 */
	@Override
	public PlotWindow show()
	{
		draw();
		ImageProcessor ip = getProcessor();
		if (Prefs.useInvertingLut && (ip instanceof ByteProcessor) && !Interpreter.isBatchMode() &&
				IJ.getInstance() != null)
		{
			ip.invertLut();
			ip.invert();
		}
		if ((IJ.macroRunning() && IJ.getInstance() == null) || Interpreter.isBatchMode())
		{
			String title = "";
			try
			{
				Field f = Plot.class.getDeclaredField("title");
				f.setAccessible(true);
				title = f.get(this).toString();
			}
			catch (Exception e)
			{
				// Ignore
			}
			ImagePlus imp = new ImagePlus(title, ip);
			WindowManager.setTempCurrentImage(imp);
			imp.setProperty("XValues", xValues); //Allows values to be retrieved by 
			imp.setProperty("YValues", yValues); // by Plot.getValues() macro function
			Interpreter.addBatchModeImage(imp);
			return null;
		}
		ImageWindow.centerNextImage();
		SuperPlotWindow pw = new SuperPlotWindow(this);
		ImagePlus imp = pw.getImagePlus();
		if (IJ.isMacro() && imp != null) // wait for plot to be displayed
			IJ.selectWindow(imp.getID());
		return pw;
	}
}
