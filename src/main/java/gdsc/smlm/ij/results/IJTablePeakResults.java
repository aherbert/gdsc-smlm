package gdsc.smlm.ij.results;

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

import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.ij.utils.CoordinateProvider;
import gdsc.smlm.ij.utils.ImageROIPainter;
import gdsc.smlm.results.Calibration;
import gdsc.smlm.results.PeakResult;
import ij.WindowManager;
import ij.text.TextPanel;
import ij.text.TextWindow;

import java.awt.Frame;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;

/**
 * Saves the fit results to an ImageJ results table.
 * <p>
 * The table supports mouse click events to draw the selected coordinates on the original source image using the
 * ImageROIPainter.
 */
public class IJTablePeakResults extends IJAbstractPeakResults implements CoordinateProvider
{
	// Store the ROI painters that have been attached to TextPanels so they can be updated
	// with a new image source
	private static HashMap<TextPanel, ImageROIPainter> map = new HashMap<TextPanel, ImageROIPainter>();

	private boolean showDeviations = true;
	private boolean showEndFrame = false;
	private boolean clearAtStart = false;
	private boolean showCalibratedValues = false;
	private String peakIdColumnName = "Peak";
	private String source = null;
	private String sourceText = null;
	private String tableTitle = "Fit Results";
	private TextWindow resultsWindow;
	private double gain = 1;
	private double nmPerPixel = 1;
	private boolean addCounter = false;

	private int indexT = -1, indexX = -1, indexY = -1;

	private int size = 0;

	public IJTablePeakResults(boolean showDeviations)
	{
		this.showDeviations = showDeviations;
	}

	public IJTablePeakResults(boolean showDeviations, String source)
	{
		this.showDeviations = showDeviations;
		this.source = source;
	}

	public IJTablePeakResults(boolean showDeviations, String source, boolean clearAtStart)
	{
		this.showDeviations = showDeviations;
		this.source = source;
		this.clearAtStart = clearAtStart;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.PeakResults#begin()
	 */
	public void begin()
	{
		createSourceText();
		createResultsWindow();
		if (clearAtStart)
		{
			resultsWindow.getTextPanel().clear();
		}
		if (showCalibratedValues)
		{
			Calibration cal = getCalibration();
			if (cal != null)
			{
				gain = cal.gain;
				nmPerPixel = cal.nmPerPixel;
			}
			else
			{
				gain = 1;
				nmPerPixel = 1;
			}
		}
	}

	/**
	 * Create the result window (if it is not available)
	 */
	private void createResultsWindow()
	{
		String header = createResultsHeader();

		ImageROIPainter roiPainter = null;
		for (Frame f : WindowManager.getNonImageWindows())
		{
			if (f != null && tableTitle.equals(f.getTitle()) && f instanceof TextWindow)
			{
				resultsWindow = (TextWindow) f;

				// Check if the existing table matches the desired header
				String currentHeader = resultsWindow.getTextPanel().getColumnHeadings();
				if (!currentHeader.startsWith(header))
				{
					resultsWindow = null;
					continue;
				}

				roiPainter = map.get(resultsWindow.getTextPanel());
			}
		}

		if (resultsWindow == null || !resultsWindow.isShowing())
		{
			resultsWindow = new TextWindow(tableTitle, header, "", 800, 300);
			roiPainter = new ImageROIPainter(resultsWindow.getTextPanel(), "", this);

			// The ROI painter adds itself to the TextPanel as a mouse listener. However
			// the TextPanel addMouseListener() adds to the private TextCanvas object so it 
			// cannot be retrieved. Store the painter in a global lookup table.
			map.put(resultsWindow.getTextPanel(), roiPainter);
		}

		if (roiPainter != null && getSource() != null)
		{
			roiPainter.setTitle(getSource().getOriginal().getName());

			// Update the coordinate provider (avoids memory leaks with old objects lying around)
			roiPainter.setCoordProvider(this);

			// Get the headings for extracting the coordinates 
			String[] headings = resultsWindow.getTextPanel().getColumnHeadings().split("\t");
			for (int i = 0; i < headings.length; i++)
			{
				if (headings[i].equals(peakIdColumnName))
				{
					indexT = i;
					continue;
				}
				if (headings[i].equals("X"))
				{
					indexX = i;
					continue;
				}
				if (headings[i].equals("Y"))
				{
					indexY = i;
					continue;
				}
			}
		}
	}

	private String createResultsHeader()
	{
		StringBuilder sb = new StringBuilder();
		if (addCounter)
			sb.append("#\t");
		if (sourceText != null)
			sb.append("Source\t");
		sb.append(peakIdColumnName);
		if (showEndFrame)
			sb.append("\tEnd ").append(peakIdColumnName);
		sb.append("\torigX");
		sb.append("\torigY");
		sb.append("\torigValue");
		sb.append("\tError");
		sb.append("\tNoise");
		sb.append("\tSNR");
		sb.append("\tBackground");
		addDeviation(sb);
		sb.append("\tSignal");
		addDeviation(sb);
		sb.append("\tAngle");
		addDeviation(sb);
		sb.append("\tX");
		addDeviation(sb);
		sb.append("\tY");
		addDeviation(sb);
		sb.append("\tX SD");
		addDeviation(sb);
		sb.append("\tY SD");
		addDeviation(sb);
		if (this.calibration != null)
		{
			sb.append("\tPrecision (nm)");
		}
		return sb.toString();
	}

	private void createSourceText()
	{
		StringBuilder sb = new StringBuilder();
		if (source != null)
			sb.append(source);
		else if (getSource() != null)
			sb.append(getSource().getName());
		if (bounds != null)
		{
			if (sb.length() > 0)
				sb.append(": ");
			sb.append(getBoundsString());
		}
		if (sb.length() > 0)
		{
			sb.append("\t");
			sourceText = sb.toString();
		}
	}

	private void addDeviation(StringBuilder sb)
	{
		if (showDeviations)
			sb.append("\t+/-");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.AbstractPeakResults#add(int, int, int, float, double, float, float[], float[])
	 */
	public void add(int peak, int origX, int origY, float origValue, double error, float noise, float[] params,
			float[] paramsDev)
	{
		addPeak(peak, peak, origX, origY, origValue, error, noise, params, paramsDev);
	}

	private void addPeak(int peak, int endFrame, int origX, int origY, float origValue, double error, float noise,
			float[] params, float[] paramsDev)
	{
		float precision = 0;
		if (this.calibration != null)
		{
			final double s = (params[Gaussian2DFunction.X_SD] + params[Gaussian2DFunction.Y_SD]) * 0.5 *
					calibration.nmPerPixel;
			precision = (float) PeakResult.getPrecision(calibration.nmPerPixel, s, params[Gaussian2DFunction.SIGNAL] /
					calibration.gain, noise / calibration.gain, calibration.emCCD);
		}
		final float snr = (noise > 0) ? params[Gaussian2DFunction.SIGNAL] / noise : 0;
		if (showCalibratedValues)
		{
			// Do not calibrate the original values
			//origX *= nmPerPixel;
			//origY *= nmPerPixel;
			//origValue /= gain;
			noise /= gain;
			params = Arrays.copyOf(params, params.length);
			params[Gaussian2DFunction.SIGNAL] /= gain;
			params[Gaussian2DFunction.BACKGROUND] = (float) ((params[Gaussian2DFunction.BACKGROUND] - calibration.bias) / gain);
			params[Gaussian2DFunction.X_POSITION] *= nmPerPixel;
			params[Gaussian2DFunction.X_SD] *= nmPerPixel;
			params[Gaussian2DFunction.Y_POSITION] *= nmPerPixel;
			params[Gaussian2DFunction.Y_SD] *= nmPerPixel;
			if (paramsDev != null)
			{
				paramsDev = Arrays.copyOf(paramsDev, paramsDev.length);
				paramsDev[Gaussian2DFunction.SIGNAL] /= gain;
				paramsDev[Gaussian2DFunction.BACKGROUND] /= gain;
				paramsDev[Gaussian2DFunction.X_POSITION] *= nmPerPixel;
				paramsDev[Gaussian2DFunction.X_SD] *= nmPerPixel;
				paramsDev[Gaussian2DFunction.Y_POSITION] *= nmPerPixel;
				paramsDev[Gaussian2DFunction.Y_SD] *= nmPerPixel;
			}
		}
		if (showDeviations)
		{
			if (paramsDev == null)
				paramsDev = new float[7];
			addResult(peak, endFrame, origX, origY, origValue, error, noise, snr, params[0], paramsDev[0], params[1],
					paramsDev[1], params[2], paramsDev[2], params[3], paramsDev[3], params[4], paramsDev[4], params[5],
					paramsDev[5], params[6], paramsDev[6], precision);
		}
		else
		{
			addResult(peak, endFrame, origX, origY, origValue, error, noise, snr, params[0], params[1], params[2],
					params[3], params[4], params[5], params[6], precision);
		}
	}

	private void addResult(int peak, int endPeak, float origX, float origY, float origValue, double error,
			float... args)
	{
		StringBuilder sb = new StringBuilder();
		if (addCounter)
			sb.append(size + 1).append("\t");
		if (sourceText != null)
			sb.append(sourceText);
		// Do not calibrate the original values		
		//if (showCalibratedValues)
		//	sb.append(peak).append(String.format("\t%g", origX)).append(String.format("\t%g", origY));
		//else
		sb.append(peak);
		if (showEndFrame)
			sb.append("\t").append(endPeak);
		sb.append(String.format("\t%.0f", origX)).append(String.format("\t%.0f", origY));
		sb.append(String.format("\t%g", origValue));
		sb.append(String.format("\t%g", error));
		for (float f : args)
			sb.append(String.format("\t%g", f));

		addResult(sb.toString());
	}

	private synchronized void addResult(String result)
	{
		size++;
		resultsWindow.append(result);
	}

	public void addAll(Collection<PeakResult> results)
	{
		for (PeakResult result : results)
			addPeak(result.peak, result.getEndFrame(), result.origX, result.origY, result.origValue, result.error,
					result.noise, result.params, result.paramsStdDev);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.PeakResults#size()
	 */
	public int size()
	{
		return size;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.PeakResults#end()
	 */
	public void end()
	{
		// Nothing to do
	}

	/**
	 * @return the name of the peak column
	 */
	public String getPeakIdColumnName()
	{
		return peakIdColumnName;
	}

	/**
	 * @param peakIdColumnName
	 *            the name of the peak column
	 */
	public void setPeakIdColumnName(String peakIdColumnName)
	{
		this.peakIdColumnName = peakIdColumnName;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#isActive()
	 */
	public boolean isActive()
	{
		return resultsWindow != null && resultsWindow.isShowing();
	}

	/**
	 * @return the table title
	 */
	public String getTableTitle()
	{
		return tableTitle;
	}

	/**
	 * Use to set the title of the table. If an existing table exists with the same title then it will be appended,
	 * otherwise a new table is created.
	 * 
	 * @param tableTitle
	 *            the table title
	 */
	public void setTableTitle(String tableTitle)
	{
		if (tableTitle != null && tableTitle.length() > 0)
			this.tableTitle = tableTitle;
	}

	/**
	 * @return True if the deviations of the parameters should be shown
	 */
	public boolean isShowDeviations()
	{
		return showDeviations;
	}

	/**
	 * @param showDeviations
	 *            True if the deviations of the parameters should be shown
	 */
	public void setShowDeviations(boolean showDeviations)
	{
		this.showDeviations = showDeviations;
	}

	/**
	 * @return True if the table should be cleared in {@link #begin()}
	 */
	public boolean isClearAtStart()
	{
		return clearAtStart;
	}

	/**
	 * @param clearAtStart
	 *            True if the table should be cleared in {@link #begin()}
	 */
	public void setClearAtStart(boolean clearAtStart)
	{
		this.clearAtStart = clearAtStart;
	}

	/**
	 * @return True if the calibration shoudl be used to adjust the table values
	 */
	public boolean isShowCalibratedValues()
	{
		return showCalibratedValues;
	}

	/**
	 * @param showCalibratedValues
	 *            True if the calibration shoudl be used to adjust the table values
	 */
	public void setShowCalibratedValues(boolean showCalibratedValues)
	{
		this.showCalibratedValues = showCalibratedValues;
	}

	/**
	 * @return the addCounter
	 */
	public boolean isAddCounter()
	{
		return addCounter;
	}

	/**
	 * @param addCounter
	 *            the addCounter to set
	 */
	public void setAddCounter(boolean addCounter)
	{
		this.addCounter = addCounter;
	}

	/**
	 * @return the resultsWindow
	 */
	public TextWindow getResultsWindow()
	{
		return resultsWindow;
	}

	/**
	 * @param line
	 * @return
	 */
	public double[] getCoordinates(String line)
	{
		// Extract the startT and x,y coordinates from the PeakResult line
		String[] fields = line.split("\t");
		try
		{
			int startT = Integer.valueOf(fields[indexT]);
			double x = Double.valueOf(fields[indexX]);
			double y = Double.valueOf(fields[indexY]);
			if (showCalibratedValues)
			{
				x /= nmPerPixel;
				y /= nmPerPixel;
			}
			return new double[] { startT, x, y };
		}
		catch (ArrayIndexOutOfBoundsException e)
		{
			// Will happen if any index is still at the default of -1 or if there are not enough fields
		}
		catch (NumberFormatException e)
		{
			// In case any field is not a number
		}
		return null;
	}

	/**
	 * @return If true show the results end frame in the table
	 */
	public boolean isShowEndFrame()
	{
		return showEndFrame;
	}

	/**
	 * @param showEndFrame
	 *            If true show the results end frame in the table
	 */
	public void setShowEndFrame(boolean showEndFrame)
	{
		this.showEndFrame = showEndFrame;
	}
}
