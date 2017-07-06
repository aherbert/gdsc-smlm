package gdsc.smlm.ij.plugins;

import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.FastMath;

import gdsc.core.ij.Utils;
import gdsc.core.utils.StoredDataStatistics;
import gdsc.smlm.data.config.UnitConfig.DistanceUnit;
import gdsc.smlm.data.config.UnitConfig.IntensityUnit;
import gdsc.smlm.engine.FitConfiguration;

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

import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.results.IJTablePeakResults;
import gdsc.smlm.results.Counter;
import gdsc.smlm.results.ImageSource;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.procedures.HeightResultProcedure;
import gdsc.smlm.results.procedures.PeakResultProcedure;
import gdsc.smlm.results.procedures.PrecisionResultProcedure;
import gdsc.smlm.results.procedures.StandardResultProcedure;
import gdsc.smlm.results.procedures.WidthResultProcedure;
import gdsc.smlm.results.procedures.XYResultProcedure;
import gdsc.smlm.utils.XmlUtils;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.ExtendedGenericDialog;
import ij.gui.Plot2;
import ij.gui.PointRoi;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.text.TextPanel;

/**
 * Extract the spots from the original image into a stack, ordering the spots by various rankings.
 */
public class SpotInspector implements PlugIn, MouseListener
{
	private static final String TITLE = "Spot Inspector";

	private static String inputOption = "";
	private static String[] SORT_ORDER = new String[] { "SNR", "Precision", "Amplitude", "Signal", "Error",
			"Original Value", "X SD", "Y SD", "Width factor", "Shift" };
	private static int sortOrderIndex = 1;
	private static int radius = 5;
	private static boolean showCalibratedValues = true;
	private static boolean plotScore = true;
	private static boolean plotHistogram = true;
	private static int histogramBins = 100;
	private static boolean removeOutliers = true;

	private MemoryPeakResults results;
	private TextPanel textPanel;
	private List<PeakResultRank> rankedResults;

	private static int currentId = 0;
	private int id;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (MemoryPeakResults.isMemoryEmpty())
		{
			IJ.error(TITLE, "No localisations in memory");
			return;
		}

		if (!showDialog())
			return;

		// Load the results
		results = ResultsManager.loadInputResults(inputOption, false, DistanceUnit.PIXEL, null);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			IJ.showStatus("");
			return;
		}

		// Check if the original image is open
		ImageSource source = results.getSource();
		if (source == null)
		{
			IJ.error(TITLE, "Unknown original source image");
			return;
		}
		source = source.getOriginal();
		if (!source.open())
		{
			IJ.error(TITLE, "Cannot open original source image: " + source.toString());
			return;
		}
		final float stdDevMax = getStandardDeviation(results);
		if (stdDevMax < 0)
		{
			// TODO - Add dialog to get the initial peak width
			IJ.error(TITLE, "Fitting configuration (for initial peak width) is not available");
			return;
		}

		// Rank spots
		rankedResults = new ArrayList<PeakResultRank>(results.size());

		// Data for the sorting
		final PrecisionResultProcedure pp;
		if (sortOrderIndex == 1)
		{
			pp = new PrecisionResultProcedure(results);
			pp.getPrecision();
		}
		else
			pp = null;

		// Build procedures to get:
		// Shift = position in pixels - originXY
		final StandardResultProcedure sp;
		if (sortOrderIndex == 9)
		{
			sp = new StandardResultProcedure(results, DistanceUnit.PIXEL);
			sp.getXYR();
		}
		else
			sp = null;

		// SD = gaussian widths only for Gaussian PSFs
		final WidthResultProcedure wp;
		if (sortOrderIndex >= 6 && sortOrderIndex <= 8)
		{
			wp = new WidthResultProcedure(results, DistanceUnit.PIXEL);
			wp.getWxWy();
		}
		else
			wp = null;

		// Amplitude for Gaussian PSFs
		final HeightResultProcedure hp;
		if (sortOrderIndex == 2)
		{
			hp = new HeightResultProcedure(results, IntensityUnit.PHOTON);
			hp.getH();
		}
		else
			hp = null;

		final Counter c = new Counter();
		results.forEach(new PeakResultProcedure()
		{
			public void execute(PeakResult r)
			{
				float[] score = getScore(r, c.getAndIncrement(), pp, sp, wp, hp, stdDevMax);
				rankedResults.add(new PeakResultRank(r, score[0], score[1]));
			}
		});
		Collections.sort(rankedResults);

		// Prepare results table
		IJTablePeakResults table = new IJTablePeakResults(false, results.getName(), true);
		table.copySettings(results);
		table.setTableTitle(TITLE);
		table.setAddCounter(true);
		table.setShowZ(results.is3D());
		// TODO - Add to settings
		table.setShowFittingData(true);
		table.setShowNoiseData(true);
		
		if (showCalibratedValues)
		{
			table.setDistanceUnit(DistanceUnit.NM);
			table.setIntensityUnit(IntensityUnit.PHOTON);
		}
		else
		{
			table.setDistanceUnit(DistanceUnit.PIXEL);
			table.setIntensityUnit(IntensityUnit.COUNT);
		}
		table.begin();

		// Add a mouse listener to jump to the frame for the clicked line
		textPanel = table.getResultsWindow().getTextPanel();

		// We must ignore old instances of this class from the mouse listeners
		id = ++currentId;
		textPanel.addMouseListener(this);

		// Add results to the table
		int n = 0;
		for (PeakResultRank rank : rankedResults)
		{
			rank.rank = n++;
			table.add(rank.peakResult);
		}
		table.end();

		if (plotScore || plotHistogram)
		{
			// Get values for the plots
			float[] xValues = null, yValues = null;
			double yMin, yMax;

			int spotNumber = 0;
			xValues = new float[rankedResults.size()];
			yValues = new float[xValues.length];
			for (PeakResultRank rank : rankedResults)
			{
				xValues[spotNumber] = spotNumber + 1;
				yValues[spotNumber++] = recoverScore(rank.score);
			}

			// Set the min and max y-values using 1.5 x IQR 
			DescriptiveStatistics stats = new DescriptiveStatistics();
			for (float v : yValues)
				stats.addValue(v);
			if (removeOutliers)
			{
				double lower = stats.getPercentile(25);
				double upper = stats.getPercentile(75);
				double iqr = upper - lower;

				yMin = FastMath.max(lower - iqr, stats.getMin());
				yMax = FastMath.min(upper + iqr, stats.getMax());

				IJ.log(String.format("Data range: %f - %f. Plotting 1.5x IQR: %f - %f", stats.getMin(), stats.getMax(),
						yMin, yMax));
			}
			else
			{
				yMin = stats.getMin();
				yMax = stats.getMax();

				IJ.log(String.format("Data range: %f - %f", yMin, yMax));
			}

			plotScore(xValues, yValues, yMin, yMax);
			plotHistogram(yValues, yMin, yMax);
		}

		// Extract spots into a stack
		final int w = source.getWidth();
		final int h = source.getHeight();
		final int size = 2 * radius + 1;
		ImageStack spots = new ImageStack(size, size, rankedResults.size());

		// To assist the extraction of data from the image source, process them in time order to allow 
		// frame caching. Then set the appropriate slice in the result stack
		Collections.sort(rankedResults, new Comparator<PeakResultRank>()
		{
			public int compare(PeakResultRank o1, PeakResultRank o2)
			{
				if (o1.peakResult.getFrame() < o2.peakResult.getFrame())
					return -1;
				if (o1.peakResult.getFrame() > o2.peakResult.getFrame())
					return 1;
				return 0;
			}
		});

		for (PeakResultRank rank : rankedResults)
		{
			PeakResult r = rank.peakResult;

			// Extract image
			// Note that the coordinates are relative to the middle of the pixel (0.5 offset)
			// so do not round but simply convert to int
			final int x = (int) (r.getXPosition());
			final int y = (int) (r.getYPosition());

			// Extract a region but crop to the image bounds
			int minX = x - radius;
			int minY = y - radius;
			int maxX = FastMath.min(x + radius + 1, w);
			int maxY = FastMath.min(y + radius + 1, h);

			int padX = 0, padY = 0;
			if (minX < 0)
			{
				padX = -minX;
				minX = 0;
			}
			if (minY < 0)
			{
				padY = -minY;
				minY = 0;
			}
			int sizeX = maxX - minX;
			int sizeY = maxY - minY;

			float[] data = source.get(r.getFrame(), new Rectangle(minX, minY, sizeX, sizeY));
			// Prevent errors with missing data
			if (data == null)
				data = new float[sizeX * sizeY];
			ImageProcessor spotIp = new FloatProcessor(sizeX, sizeY, data, null);

			// Pad if necessary, i.e. the crop is too small for the stack
			if (padX > 0 || padY > 0 || sizeX < size || sizeY < size)
			{
				ImageProcessor spotIp2 = spotIp.createProcessor(size, size);
				spotIp2.insert(spotIp, padX, padY);
				spotIp = spotIp2;
			}
			int slice = rank.rank + 1;
			spots.setPixels(spotIp.getPixels(), slice);
			spots.setSliceLabel(Utils.rounded(rank.originalScore), slice);
		}

		source.close();

		ImagePlus imp = Utils.display(TITLE, spots);
		imp.setRoi((PointRoi) null);

		// Make bigger		
		for (int i = 10; i-- > 0;)
			imp.getWindow().getCanvas().zoomIn(imp.getWidth() / 2, imp.getHeight() / 2);
	}

	private float getStandardDeviation(MemoryPeakResults results2)
	{
		// Standard deviation is only needed for the width filtering
		if (sortOrderIndex != 8)
			return 0;
		// XXX - This should be refactored out to the results configuration (i.e. store the PSF type)
		FitEngineConfiguration config = (FitEngineConfiguration) XmlUtils.fromXML(results.getConfiguration());
		if (config == null || config.getFitConfiguration() == null)
		{
			return -1;
		}
		FitConfiguration fitConfig = config.getFitConfiguration();
		float stdDevMax = (float) ((fitConfig.getInitialXSD() > 0) ? fitConfig.getInitialXSD() : 1);
		if (fitConfig.getInitialYSD() > 0)
			stdDevMax = (float) FastMath.max(fitConfig.getInitialYSD(), stdDevMax);
		return stdDevMax;
	}

	private void plotScore(float[] xValues, float[] yValues, double yMin, double yMax)
	{
		if (plotScore)
		{
			String title = TITLE + " Score";
			Plot2 plot = new Plot2(title, "Rank", SORT_ORDER[sortOrderIndex], xValues, yValues);
			plot.setLimits(1, xValues.length, yMin, yMax);

			Utils.display(title, plot);
		}
	}

	private void plotHistogram(float[] data, double yMin, double yMax)
	{
		if (plotHistogram)
		{
			String title = TITLE + " Histogram";
			Utils.showHistogram(title, new StoredDataStatistics(data), SORT_ORDER[sortOrderIndex], 0,
					(removeOutliers) ? 1 : 0, histogramBins);
		}
	}

	private float[] getScore(PeakResult r, int i, PrecisionResultProcedure pp, StandardResultProcedure sp,
			WidthResultProcedure wp, HeightResultProcedure hp, float stdDevMax)
	{
		// Return score so high is better
		float score;
		boolean negative = false;
		switch (sortOrderIndex)
		{
			case 9: // Shift
				// We do not have the original centroid so use the original X/Y
				score = FastMath.max(sp.x[i] - r.origX + 0.5f, sp.y[i] - r.origY + 0.5f);
				negative = true;
				break;
			case 8: // Width factor
				score = getFactor(FastMath.max(wp.wx[i], wp.wy[i]), stdDevMax);
				negative = true;
				break;
			case 7:
				score = wp.wy[i];
				negative = true;
				break;
			case 6:
				score = wp.wx[i];
				negative = true;
				break;
			case 5: // Original value
				score = (float) r.origValue;
				break;
			case 4: // Error
				score = (float) r.error;
				negative = true;
				break;
			case 3: // Signal
				score = (float) (r.getSignal());
				break;
			case 2: // Amplitude
				score = hp.h[i];
				break;
			case 1: // Precision
				score = (float) pp.precision[i];
				negative = true;
				break;
			default: // SNR
				score = r.getSignal() / r.noise;
		}
		return new float[] { (negative) ? -score : score, score };
	}

	private float recoverScore(float score)
	{
		// Reset the sign of the score
		switch (sortOrderIndex)
		{
			case 9: // Shift
				return -score;
			case 8: // Width factor
				return -score;
			case 7:
				return -score;
			case 6:
				return -score;
			case 5: // Original value
				return score;
			case 4: // Error
				return -score;
			case 3: // Signal
				return score;
			case 2: // Amplitude
				return score;
			case 1: // Precision
				return -score;
			default: // SNR
				return score;
		}
	}

	/**
	 * Get the relative change factor between f and g
	 * 
	 * @param f
	 * @param g
	 * @return
	 */
	private static float getFactor(float f, float g)
	{
		if (f > g)
			return f / g;
		return g / f;
	}

	private class PeakResultRank implements Comparable<PeakResultRank>
	{
		int rank;
		PeakResult peakResult;
		float score, originalScore;

		public PeakResultRank(PeakResult r, float s, float original)
		{
			peakResult = r;
			score = s;
			originalScore = original;
		}

		public int compareTo(PeakResultRank o)
		{
			// High is better
			if (score > o.score)
				return -1;
			if (score < o.score)
				return 1;
			return 0;
		}
	}

	private boolean showDialog()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);

		gd.addChoice("Ranking", SORT_ORDER, SORT_ORDER[sortOrderIndex]);
		gd.addSlider("Radius", 1, 15, radius);
		gd.addCheckbox("Calibrated_table", showCalibratedValues);
		gd.addCheckbox("Plot_score", plotScore);
		gd.addCheckbox("Plot_histogram", plotHistogram);
		gd.addNumericField("Histogram_bins", histogramBins, 0);
		gd.addCheckbox("Remove_outliers", removeOutliers);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		inputOption = ResultsManager.getInputSource(gd);
		sortOrderIndex = gd.getNextChoiceIndex();
		radius = (int) gd.getNextNumber();
		showCalibratedValues = gd.getNextBoolean();
		plotScore = gd.getNextBoolean();
		plotHistogram = gd.getNextBoolean();
		histogramBins = (int) gd.getNextNumber();
		removeOutliers = gd.getNextBoolean();

		// Check arguments
		try
		{
			Parameters.isAboveZero("Radius", radius);
			Parameters.isAbove("Histogram bins", histogramBins, 1);
		}
		catch (IllegalArgumentException ex)
		{
			IJ.error(TITLE, ex.getMessage());
			return false;
		}

		return true;
	}

	public void mouseClicked(MouseEvent e)
	{
		if (id != currentId)
			return;
		// Show the result that was double clicked in the result table
		if (e.getClickCount() > 1)
		{
			int rank = textPanel.getSelectionStart() + 1;

			// Show the spot that was double clicked
			ImagePlus imp = WindowManager.getImage(TITLE);
			if (imp != null && rank > 0 && rank <= imp.getStackSize())
			{
				imp.setSlice(rank);
				if (imp.getWindow() != null)
					imp.getWindow().toFront();

				// Add an ROI to to the image containing all the spots in that part of the frame
				PeakResult r = rankedResults.get(rank - 1).peakResult;

				final int x = (int) (r.getXPosition());
				final int y = (int) (r.getYPosition());

				// Find bounds
				final int minX = x - radius;
				final int minY = y - radius;
				final int maxX = x + radius + 1;
				final int maxY = y + radius + 1;

				// Create ROIs
				final ArrayList<float[]> spots = new ArrayList<float[]>();
				results.forEach(DistanceUnit.PIXEL, new XYResultProcedure()
				{
					public void executeXY(float x, float y)
					{
						if (x > minX && x < maxX && y > minY && y < maxY)
						{
							// Use only unique points
							final float xPosition = x - minX;
							final float yPosition = y - minY;
							if (!contains(spots, xPosition, yPosition))
								spots.add(new float[] { xPosition, yPosition });
						}
					}
				});

				int points = spots.size();
				float[] ox = new float[points];
				float[] oy = new float[points];
				for (int i = 0; i < points; i++)
				{
					ox[i] = spots.get(i)[0];
					oy[i] = spots.get(i)[1];
				}
				imp.setRoi(new PointRoi(ox, oy, points));
			}
		}
	}

	private boolean contains(ArrayList<float[]> spots, float xPosition, float yPosition)
	{
		for (float[] data : spots)
			if (data[0] == xPosition && data[1] == yPosition)
				return true;
		return false;
	}

	public void mousePressed(MouseEvent e)
	{
		// TODO Auto-generated method stub

	}

	public void mouseReleased(MouseEvent e)
	{
		// TODO Auto-generated method stub

	}

	public void mouseEntered(MouseEvent e)
	{
		// TODO Auto-generated method stub

	}

	public void mouseExited(MouseEvent e)
	{
		// TODO Auto-generated method stub

	}

}
