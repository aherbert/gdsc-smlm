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
import java.lang.reflect.Method;

/**
 * Extension of the ij.gui.Plot class to add functionality
 */
public class Plot2 extends Plot
{
	/** Draw a bar plot */
	public static final int BAR = 999;
	
	private static boolean failedOverride = false;

	public Plot2(String title, String xLabel, String yLabel, float[] xValues, float[] yValues)
	{
		super(title, xLabel, yLabel, xValues, yValues);
	}

	public Plot2(String title, String xLabel, String yLabel, double[] xValues, double[] yValues)
	{
		super(title, xLabel, yLabel, xValues, yValues);
	}

	public Plot2(String dummy, String title, String xLabel, String yLabel, float[] xValues, float[] yValues)
	{
		super(title, xLabel, yLabel, xValues, yValues);
	}

	public Plot2(String title, String xLabel, String yLabel)
	{
		super(title, xLabel, yLabel, (float[]) null, (float[]) null);
	}

	public Plot2(String title, String xLabel, String yLabel, int flags)
	{
		super(title, xLabel, yLabel, (float[]) null, (float[]) null, flags);
	}

	public Plot2(String title, String xLabel, String yLabel, float[] xValues, float[] yValues, int flags)
	{
		super(title, xLabel, yLabel, xValues, yValues, flags);
	}

	public Plot2(String title, String xLabel, String yLabel, double[] xValues, double[] yValues, int flags)
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
		// Override to allow a Bar plot. If this fails due to an exception then a line plot will be used.
		if (failedOverride)
		{
			if (shape == Plot2.BAR)
				shape = Plot.LINE;
			super.addPoints(x, y, shape);
			return;
		}

		// Set the limits if this is the first set of data. The limits are usually set in the constructor
		// but we may want to not pass in the values to the constructor and then immediately call 
		// addPoints(x, y, Plot2.BAR)
		boolean setLimits = false;
		try
		{
			// Check for any values
			Field f = Plot.class.getDeclaredField("xValues");
			f.setAccessible(true);
			Object value = f.get(this);
			setLimits = (value == null);
		}
		catch (Throwable e)
		{
			// Ignore
		}
		if (setLimits)
		{
			try
			{
				double[] a = Tools.getMinMax(x);
				setFieldValue("xMin", a[0]);
				setFieldValue("xMax", a[1]);
				a = Tools.getMinMax(y);
				setFieldValue("yMin", a[0]);
				setFieldValue("yMax", a[1]);
			}
			catch (Exception e)
			{
				// Ignore
			}
		}

		// This only works if the addPoints super method ignores the BAR option but still store the values
		try
		{
			Method m = Plot.class.getDeclaredMethod("setup");
			m.setAccessible(true);
			m.invoke(this);

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
					catch (Throwable e)
					{
						// Ignore
					}
					float[] xValues = createHistogramAxis(x);
					float[] yValues = createHistogramValues(y);

					m = Plot.class.getDeclaredMethod("drawFloatPolyline", ImageProcessor.class, float[].class,
							float[].class, int.class);
					m.setAccessible(true);

					if ((flags & X_LOG_NUMBERS) != 0)
						xValues = arrayToLog(xValues);
					if ((flags & Y_LOG_NUMBERS) != 0)
						yValues = arrayToLog(yValues);

					m.invoke(this, ip, xValues, yValues, xValues.length);
					break;
			}
		}
		catch (Throwable e)
		{
			// Revert to drawing a line
			failedOverride = true;
			shape = Plot.LINE;
		}

		super.addPoints(x, y, shape);
	}

	private void setFieldValue(String name, Object value) throws Exception
	{
		Field f = Plot.class.getDeclaredField(name);
		f.setAccessible(true);
		f.set(this, value);
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
		// Override to show a PlotWindow2 object
		if (failedOverride)
			return super.show();
			
		try
		{
			Method m = Plot.class.getDeclaredMethod("draw");
			m.setAccessible(true);
			m.invoke(this);

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
				catch (Throwable e)
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
			// This may throw an IllegalAccessError on some platforms
			PlotWindow2 pw = new PlotWindow2(this);
			ImagePlus imp = pw.getImagePlus();
			if (IJ.isMacro() && imp != null) // wait for plot to be displayed
				IJ.selectWindow(imp.getID());
			return pw;
		}
		catch (Throwable e)
		{
			// Ignore
			failedOverride = true;
		}

		return super.show();
	}
}
