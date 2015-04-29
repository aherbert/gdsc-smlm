package gdsc.smlm.ij.utils;

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

import gdsc.smlm.utils.Maths;
import gdsc.smlm.utils.StoredDataStatistics;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Macro;
import ij.WindowManager;
import ij.gui.ImageWindow;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.gui.ProgressBar;
import ij.gui.Plot2;
import ij.io.DirectoryChooser;
import ij.io.OpenDialog;
import ij.plugin.frame.Recorder;
import ij.process.ImageProcessor;
import ij.text.TextWindow;

import java.awt.Frame;
import java.awt.Window;
import java.lang.reflect.Field;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import org.apache.commons.math3.util.FastMath;

/**
 * Contains helper functions
 */
public class Utils
{
	private static boolean newWindow = false;

	/**
	 * Splits a full path into the directory and filename
	 * 
	 * @param path
	 * @return directory and filename
	 */
	public static String[] decodePath(String path)
	{
		String[] result = new String[2];
		if (path == null)
			path = "";
		int i = path.lastIndexOf('/');
		if (i == -1)
			i = path.lastIndexOf('\\');
		if (i > 0)
		{
			result[0] = path.substring(0, i + 1);
			result[1] = path.substring(i + 1);
		}
		else
		{
			result[0] = OpenDialog.getDefaultDirectory();
			result[1] = path;
		}
		return result;
	}

	/**
	 * Round the double to the specified significant digits
	 * 
	 * @param d
	 * @param significantDigits
	 * @return
	 */
	public static String rounded(double d, int significantDigits)
	{
		if (Double.isInfinite(d) || Double.isNaN(d))
			return "" + d;
		BigDecimal bd = new BigDecimal(d);
		bd = bd.round(new MathContext(significantDigits));
		return "" + bd.doubleValue();
	}

	/**
	 * Round the double to 4 significant digits
	 * 
	 * @param d
	 * @return
	 */
	public static String rounded(double d)
	{
		return rounded(d, 4);
	}

	/**
	 * Show the image. Replace a currently open image with the specified title or else create a new image.
	 * 
	 * @param title
	 * @param ip
	 * @return the
	 */
	public static ImagePlus display(String title, ImageProcessor ip)
	{
		newWindow = false;
		ImagePlus imp = WindowManager.getImage(title);
		if (imp == null)
		{
			imp = new ImagePlus(title, ip);
			imp.show();
			newWindow = true;
		}
		else
		{
			imp.setProcessor(ip);
		}
		return imp;
	}

	/**
	 * Show the image. Replace a currently open image with the specified title or else create a new image.
	 * 
	 * @param title
	 * @param slices
	 * @return the image
	 */
	public static ImagePlus display(String title, ImageStack slices)
	{
		newWindow = false;
		ImagePlus imp = WindowManager.getImage(title);
		if (imp == null)
		{
			imp = new ImagePlus(title, slices);
			imp.show();
			newWindow = true;
		}
		else
		{
			imp.setStack(slices);
		}
		return imp;
	}

	/**
	 * Show the plot. Replace a currently open plot with the specified title or else create a new plot window.
	 * 
	 * @param title
	 * @param plot
	 * @return the plot window
	 */
	public static PlotWindow display(String title, Plot plot)
	{
		newWindow = false;
		Frame plotWindow = null;
		int[] wList = WindowManager.getIDList();
		int len = wList != null ? wList.length : 0;
		for (int i = 0; i < len; i++)
		{
			ImagePlus imp = WindowManager.getImage(wList[i]);
			if (imp != null && imp.getWindow() instanceof PlotWindow)
			{
				if (imp.getTitle().equals(title))
				{
					plotWindow = imp.getWindow();
					break;
				}
			}
		}
		PlotWindow p;
		if (plotWindow == null)
		{
			p = plot.show();
			newWindow = true;
		}
		else
		{
			p = (PlotWindow) plotWindow;
			p.drawPlot(plot);
		}
		return p;
	}

	/**
	 * Close the named window
	 * 
	 * @param name
	 */
	public static void close(String name)
	{
		Window w = WindowManager.getWindow(name);
		if (w != null)
		{
			if (w instanceof ImageWindow)
				((ImageWindow) w).close();
			else if (w instanceof Frame)
			{
				WindowManager.removeWindow(w);
				w.dispose();
			}
		}
	}

	/**
	 * Calculate a histogram given the provided data
	 * 
	 * @param data
	 * @param numBins
	 *            The number of histogram bins between min and max
	 * @return The histogram as a pair of arrays: { value[], frequency[] }
	 */
	public static float[][] calcHistogram(float[] data, int numBins)
	{
		float min = Float.POSITIVE_INFINITY;
		float max = Float.NEGATIVE_INFINITY;
		for (float f : data)
		{
			if (min > f)
				min = f;
			if (max < f)
				max = f;
		}
		return calcHistogram(data, min, max, numBins);
	}

	/**
	 * Calculate a histogram given the provided data.
	 * <p>
	 * The histogram will create the specified number of bins to accommodate all data between the minimum and maximum
	 * inclusive. The number of bins must be above one so that min and max are in different bins. If min and max are the
	 * same then the number of bins is set to 1.
	 * 
	 * @param data
	 * @param min
	 *            The minimum value to include (inclusive)
	 * @param max
	 *            The maximum value to include (inclusive)
	 * @param numBins
	 *            The number of histogram bins between min and max (must be above one)
	 * @return The histogram as a pair of arrays: { value[], frequency[] }
	 */
	public static float[][] calcHistogram(float[] data, double min, double max, int numBins)
	{
		// Parameter check
		if (numBins < 2)
			numBins = 2;
		if (max < min)
		{
			double tmp = max;
			max = min;
			min = tmp;
		}
		final double binSize;
		if (max == min)
		{
			numBins = 1;
			binSize = 1;
		}
		else
		{
			binSize = (max - min) / (numBins - 1);
		}

		final float[] value = new float[numBins];
		final float[] frequency = new float[numBins];

		for (int i = 0; i < numBins; i++)
		{
			value[i] = (float) (min + i * binSize);
		}

		for (double d : data)
		{
			int bin = (int) ((d - min) / binSize);
			if (bin < 0)
			{ /* this data is smaller than min */
			}
			else if (bin >= numBins)
			{ /* this data point is bigger than max */
			}
			else
			{
				frequency[bin]++;
			}
		}

		return new float[][] { value, frequency };
	}

	/**
	 * Calculate a histogram given the provided data
	 * 
	 * @param data
	 * @param numBins
	 *            The number of histogram bins between min and max
	 * @return The histogram as a pair of arrays: { value[], frequency[] }
	 */
	public static double[][] calcHistogram(double[] data, int numBins)
	{
		double min = Double.POSITIVE_INFINITY;
		double max = Double.NEGATIVE_INFINITY;
		for (double f : data)
		{
			if (min > f)
				min = f;
			if (max < f)
				max = f;
		}
		return calcHistogram(data, min, max, numBins);
	}

	/**
	 * Calculate a histogram given the provided data.
	 * <p>
	 * The histogram will create the specified number of bins to accommodate all data between the minimum and maximum
	 * inclusive. The number of bins must be above one so that min and max are in different bins. If min and max are the
	 * same then the number of bins is set to 1.
	 * 
	 * @param data
	 * @param min
	 *            The minimum value to include (inclusive)
	 * @param max
	 *            The maximum value to include (inclusive)
	 * @param numBins
	 *            The number of histogram bins between min and max (must be above one)
	 * @return The histogram as a pair of arrays: { value[], frequency[] }
	 */
	public static double[][] calcHistogram(double[] data, double min, double max, int numBins)
	{
		// Parameter check
		if (numBins < 2)
			numBins = 2;
		if (max < min)
		{
			double tmp = max;
			max = min;
			min = tmp;
		}
		final double binSize;
		if (max == min)
		{
			numBins = 1;
			binSize = 1;
		}
		else
		{
			binSize = (max - min) / (numBins - 1);
		}

		final double[] value = new double[numBins];
		final double[] frequency = new double[numBins];

		for (int i = 0; i < numBins; i++)
		{
			value[i] = (min + i * binSize);
		}

		for (double d : data)
		{
			int bin = (int) ((d - min) / binSize);
			if (bin < 0)
			{ /* this data is smaller than min */
			}
			else if (bin >= numBins)
			{ /* this data point is bigger than max */
			}
			else
			{
				frequency[bin]++;
			}
		}

		return new double[][] { value, frequency };
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

	/**
	 * For the provided histogram x-axis bins, produce an x-axis for plotting. This functions doubles up the histogram
	 * x-positions to allow plotting a square line profile using the ImageJ plot command.
	 * 
	 * @param histogramX
	 * @return
	 */
	public static double[] createHistogramAxis(double[] histogramX)
	{
		double[] axis = new double[histogramX.length * 2 + 2];
		int index = 0;
		for (int i = 0; i < histogramX.length; ++i)
		{
			axis[index++] = histogramX[i];
			axis[index++] = histogramX[i];
		}
		if (histogramX.length > 0)
		{
			double dx = (histogramX.length == 1) ? 1 : (histogramX[1] - histogramX[0]);
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
	public static double[] createHistogramValues(double[] histogramY)
	{
		double[] axis = new double[histogramY.length * 2 + 2];

		int index = 0;
		axis[index++] = 0;
		for (int i = 0; i < histogramY.length; ++i)
		{
			axis[index++] = histogramY[i];
			axis[index++] = histogramY[i];
		}
		return axis;
	}

	/**
	 * Return the histogram statistics
	 * 
	 * @param x
	 *            Histogram values
	 * @param y
	 *            Histogram counts
	 * @return Array containing: { mean, standard deviation }
	 */
	public static double[] getHistogramStatistics(float[] x, float[] y)
	{
		// Get the average
		double n = 0;
		double sum = 0.0;
		double sum2 = 0.0;
		for (int i = 0; i < x.length; i++)
		{
			if (y[i] > 0)
			{
				float count = y[i];
				float value = x[i];
				n += count;
				sum += value * count;
				sum2 += (value * value) * count;
			}
		}
		double av = sum / n;

		// Get the Std.Dev
		double stdDev;
		if (n > 0)
		{
			double d = n;
			stdDev = (d * sum2 - sum * sum) / d;
			if (stdDev > 0.0)
				stdDev = Math.sqrt(stdDev / (d - 1.0));
			else
				stdDev = 0.0;
		}
		else
			stdDev = 0.0;

		return new double[] { av, stdDev };
	}

	/**
	 * Logs a message to the ImageJ log
	 * 
	 * @param format
	 * @param args
	 */
	public static void log(String format, Object... args)
	{
		IJ.log(String.format(format, args));
	}

	/**
	 * Check if the escape key has been pressed. Show a status aborted message if true.
	 * 
	 * @return True if aborted
	 */
	public static boolean isInterrupted()
	{
		if (IJ.escapePressed())
		{
			IJ.beep();
			IJ.showStatus("Aborted");
			return true;
		}
		return false;
	}

	/**
	 * Show a histogram of the data
	 * 
	 * @param title
	 *            The title to prepend to the plot name
	 * @param stats
	 * @param name
	 *            The name of plotted statistic
	 * @param minWidth
	 *            The minimum bin width to use (e.g. set to 1 for integer values)
	 * @param removeOutliers
	 *            Remove outliers. 1 - 1.5x IQR. 2 - remove top 2%.
	 * @param bins
	 *            The number of bins to use
	 * @return The histogram window ID
	 */
	public static int showHistogram(String title, StoredDataStatistics stats, String name, double minWidth,
			int removeOutliers, int bins)
	{
		return showHistogram(title, stats, name, minWidth, removeOutliers, bins, true, null);
	}

	/**
	 * Show a histogram of the data
	 * 
	 * @param title
	 *            The title to prepend to the plot name
	 * @param stats
	 * @param name
	 *            The name of plotted statistic
	 * @param minWidth
	 *            The minimum bin width to use (e.g. set to 1 for integer values)
	 * @param removeOutliers
	 *            Remove outliers. 1 - 1.5x IQR. 2 - remove top 2%.
	 * @param bins
	 *            The number of bins to use
	 * @param label
	 *            The label to add
	 * @return The histogram window ID
	 */
	public static int showHistogram(String title, StoredDataStatistics stats, String name, double minWidth,
			int removeOutliers, int bins, String label)
	{
		return showHistogram(title, stats, name, minWidth, removeOutliers, bins, true, label);
	}

	/**
	 * Show a histogram of the data
	 * 
	 * @param title
	 *            The title to prepend to the plot name
	 * @param stats
	 * @param name
	 *            The name of plotted statistic
	 * @param minWidth
	 *            The minimum bin width to use (e.g. set to 1 for integer values)
	 * @param removeOutliers
	 *            Remove outliers. 1 - 1.5x IQR. 2 - remove top 2%.
	 * @param bins
	 *            The number of bins to use
	 * @param barChart
	 *            Use a bar chart, else plot non-zero bin counts as a line plot
	 * @param label
	 *            The label to add
	 * @return The histogram window ID
	 */
	public static int showHistogram(String title, StoredDataStatistics stats, String name, double minWidth,
			int removeOutliers, int bins, boolean barChart, String label)
	{
		double[] values = stats.getValues();
		if (values == null || values.length < 2)
			return 0;
		double[] limits = Maths.limits(values);
		double yMin = limits[0];
		double yMax = limits[1];

		switch (removeOutliers)
		{
			case 1:
				// Get the inter quartile range
				double lower = stats.getStatistics().getPercentile(25);
				double upper = stats.getStatistics().getPercentile(75);
				double iqr = upper - lower;
				yMin = FastMath.max(lower - iqr, yMin);
				yMax = FastMath.min(upper + iqr, yMax);
				break;

			case 2:
				// Remove top 2%
				yMax = stats.getStatistics().getPercentile(98);
				break;

		}

		if (minWidth > 0)
		{
			double binSize = (yMax - yMin) / ((bins < 2) ? 1 : bins - 1);
			if (binSize < minWidth)
			{
				bins = (int) ((yMax - yMin) / minWidth) + 1;
				//yMax = bins * minWidth + yMin;
			}
		}
		//		else
		//		{
		//			// Calculate the resolution, i.e. the smallest gap between data points
		//			double resolution = Double.POSITIVE_INFINITY;
		//			for (int i=1; i<values.length; i++)
		//			{
		//				if (values[i-1] != values[i])
		//				{
		//					if (resolution > values[i] - values[i-1])
		//						resolution = values[i] - values[i-1];
		//				}
		//			}
		//			
		//			// Set the number of bins as the most needed to separate the data points. 
		//			// This prevents gaps xin the histogram
		//			if (resolution != Double.POSITIVE_INFINITY)
		//			{
		//				int numBins = 1 + (int)((yMax - yMin) / resolution);
		//				if (bins > numBins)
		//					bins = numBins;
		//			}
		//		}

		title += " " + name;

		double[][] hist = Utils.calcHistogram(values, yMin, yMax, bins);

		double[] xValues, yValues;

		if (barChart)
		{
			// Standard histogram
			xValues = hist[0]; //Utils.createHistogramAxis(hist[0]);
			yValues = hist[1]; //Utils.createHistogramValues(hist[1]);
		}
		else
		{
			// Line plot of non-zero values
			int c = 0;
			xValues = new double[hist[0].length];
			yValues = new double[xValues.length];
			for (int i = 0; i < xValues.length; i++)
			{
				if (hist[1][i] != 0)
				{
					xValues[c] = hist[0][i];
					yValues[c] = hist[1][i];
					c++;
				}
			}
			xValues = Arrays.copyOf(xValues, c);
			yValues = Arrays.copyOf(yValues, c);
		}

		Plot2 plot = new Plot2(title, name, "Frequency");
		if (xValues.length > 0)
		{
			double dx = 0;
			if (barChart)
				dx = (xValues.length == 1) ? 1 : (xValues[1] - xValues[0]);
			double xMax = xValues[xValues.length - 1] + dx;
			double xPadding = 0.05 * (xMax - xValues[0]);
			plot.setLimits(xValues[0] - xPadding, xMax + xPadding, 0,
					Maths.max(yValues) * 1.05);
		}
		plot.addPoints(xValues, yValues, (barChart) ? Plot2.BAR : Plot.LINE);
		if (label != null)
			plot.addLabel(0, 0, label);
		PlotWindow window = Utils.display(title, plot);
		return window.getImagePlus().getID();
	}

	/**
	 * @return True is the last call to display created a new window
	 */
	public static boolean isNewWindow()
	{
		return newWindow;
	}

	private static ProgressBar progressBar = null;

	/**
	 * Use reflection to replace the progress bar with null
	 * 
	 * @param showProgress
	 *            Set to true to disable the progress bar
	 */
	public static void setShowProgress(boolean showProgress)
	{
		if (progressBar == null)
		{
			progressBar = IJ.getInstance().getProgressBar();
		}

		ProgressBar newProgressBar;
		if (showProgress)
		{
			newProgressBar = progressBar;
		}
		else
		{
			newProgressBar = null;
		}

		try
		{
			Field f = IJ.class.getDeclaredField("progressBar");
			f.setAccessible(true);
			f.set(IJ.class, newProgressBar);
		}
		catch (Exception e)
		{
			// Ignore
		}
	}

	/**
	 * Convert time in milliseconds into a nice string
	 * 
	 * @param time
	 * @return
	 */
	public static String timeToString(double time)
	{
		String units = " ms";
		if (time > 1000) // 1 second
		{
			time /= 1000;
			units = " s";

			if (time > 180) // 3 minutes
			{
				time /= 60;
				units = " min";
			}
		}
		return Utils.rounded(time, 4) + units;
	}

	/**
	 * Replace the filename extension with the specified extension
	 * 
	 * @param filename
	 * @param extension
	 * @return the new filename
	 */
	public static String replaceExtension(String filename, String extension)
	{
		if (filename != null)
		{
			int index = filename.lastIndexOf('.');
			if (index > 0)
			{
				filename = filename.substring(0, index);
			}
			filename += (extension.startsWith(".")) ? extension : "." + extension;
		}
		return filename;
	}

	/**
	 * Check if the current window has the given headings, refreshing the headings if necessary.
	 * Only works if the window is showing.
	 * 
	 * @param textWindow
	 * @param headings
	 * @param preserve
	 *            Preserve the current data (note that is may not match the new headings)
	 * @return True if the window headings were changed
	 */
	public static boolean refreshHeadings(TextWindow textWindow, String headings, boolean preserve)
	{
		if (textWindow != null && textWindow.isShowing())
		{
			if (!textWindow.getTextPanel().getColumnHeadings().equals(headings))
			{
				StringBuffer sb = new StringBuffer();
				if (preserve)
					for (int i = 0; i < textWindow.getTextPanel().getLineCount(); i++)
						sb.append(textWindow.getTextPanel().getLine(i)).append("\n");

				textWindow.getTextPanel().setColumnHeadings(headings);

				if (preserve)
					textWindow.append(sb.toString());

				return true;
			}
		}
		return false;
	}

	/**
	 * Create and fill an array
	 * 
	 * @param length
	 *            The length of the array
	 * @param start
	 *            The start
	 * @param increment
	 *            The increment
	 * @return The new array
	 */
	public static double[] newArray(int length, double start, double increment)
	{
		double[] data = new double[length];
		for (int i = 0; i < length; i++, start += increment)
			data[i] = start;
		return data;
	}

	/**
	 * Create and fill an array
	 * 
	 * @param length
	 *            The length of the array
	 * @param start
	 *            The start
	 * @param increment
	 *            The increment
	 * @return The new array
	 */
	public static int[] newArray(int length, int start, int increment)
	{
		int[] data = new int[length];
		for (int i = 0; i < length; i++, start += increment)
			data[i] = start;
		return data;
	}

	/**
	 * Waits for all threads to complete computation.
	 * 
	 * @param futures
	 */
	public static void waitForCompletion(List<Future<?>> futures)
	{
		try
		{
			for (Future<?> f : futures)
			{
				f.get();
			}
		}
		catch (ExecutionException ex)
		{
			ex.printStackTrace();
		}
		catch (InterruptedException e)
		{
			e.printStackTrace();
		}
	}

	/**
	 * Open a directory selection dialog using the given title (and optionally the default directory)
	 * 
	 * @param title
	 *            The dialog title
	 * @param directory
	 *            The default directory to start in
	 * @return The directory (or null if the dialog is cancelled)
	 */
	public static String getDirectory(String title, String directory)
	{
		String defaultDir = OpenDialog.getDefaultDirectory();
		if (directory != null && directory.length() > 0)
			OpenDialog.setDefaultDirectory(directory);
		DirectoryChooser chooser = new DirectoryChooser(title);
		directory = chooser.getDirectory();
		OpenDialog.setDefaultDirectory(defaultDir);
		return directory;
	}

	/**
	 * Open a file selection dialog using the given title (and optionally the default path)
	 * 
	 * @param title
	 *            The dialog title
	 * @param filename
	 *            The default path to start with
	 * @return The path (or null if the dialog is cancelled)
	 */
	public static String getFilename(String title, String filename)
	{
		String[] path = Utils.decodePath(filename);
		OpenDialog chooser = new OpenDialog(title, path[0], path[1]);
		if (chooser.getFileName() != null)
		{
			return chooser.getDirectory() + chooser.getFileName();
		}
		return null;
	}

	/**
	 * Determine if the plugin is running with extra options. Checks for the ImageJ shift or alt key down properties. If
	 * running in a macro then searches the options string for the 'extraoptions' flag.
	 * <p>
	 * If the extra options are required then adds the 'extraoptions' flag to the macro recorder options.
	 * 
	 * @return True if extra options are required
	 */
	public static boolean isExtraOptions()
	{
		final String EXTRA = "extraoptions";
		boolean extraOptions = IJ.altKeyDown() || IJ.shiftKeyDown();
		if (!extraOptions && IJ.isMacro())
		{
			extraOptions = (Macro.getOptions() != null && Macro.getOptions().contains(EXTRA));
		}
		if (extraOptions)
			Recorder.recordOption(EXTRA);
		return extraOptions;
	}

	/**
	 * Convert the input array to a double
	 * 
	 * @param a
	 * @return The new array
	 */
	public static double[] toDouble(float[] a)
	{
		if (a == null)
			return null;
		double[] b = new double[a.length];
		for (int i = 0; i < a.length; i++)
			b[i] = a[i];
		return b;
	}

	/**
	 * Convert the input array to a float
	 * 
	 * @param a
	 * @return The new array
	 */
	public static float[] toFloat(double[] a)
	{
		if (a == null)
			return null;
		float[] b = new float[a.length];
		for (int i = 0; i < a.length; i++)
			b[i] = (float) a[i];
		return b;
	}
}
