package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.Color;
import java.util.Arrays;

import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import gdsc.core.data.utils.TypeConverter;
import gdsc.core.ij.Utils;
import gdsc.core.utils.Maths;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.core.utils.SimpleLock;
import gdsc.core.utils.Sort;
import gdsc.core.utils.Statistics;
import gdsc.core.utils.StoredData;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.data.config.UnitProtos.TimeUnit;
import gdsc.smlm.fitting.JumpDistanceAnalysis;

/*----------------------------------------------------------------------------- 
 * GDSC Plugins for ImageJ
 * 
 * Copyright (C) 2011 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.results.IdPeakResultComparator;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.procedures.PeakResultProcedure;
import gdsc.smlm.results.procedures.PrecisionResultProcedure;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import ij.IJ;
import ij.gui.DialogListener;
import ij.gui.ExtendedGenericDialog;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingExtendedGenericDialog;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.plugin.PlugIn;
import ij.plugin.WindowOrganiser;

/**
 * Analyses the track lengths of traced data
 */
public class TraceLengthAnalysis implements PlugIn, DialogListener, PeakResultProcedure
{
	private static final String TITLE = "Trace Length Analysis";
	private static final float WIDTH = 0.3f;

	private static String inputOption = "";
	private static double dThreshold = 0;
	private static boolean normalise = false;

	private TypeConverter<DistanceUnit> distanceConverter;
	private TypeConverter<TimeUnit> timeConverter;
	private SimpleLock lock = new SimpleLock();
	private double _msdThreshold = -1;
	private boolean _normalise = false;
	private int _index;
	private double error = 0;
	private double[] d; // MSD of trace
	private int[] length; // Length of trace
	private int[] id; // trace id
	private double minX, maxX;
	private int[] h1, h2;
	private float[] x1, x2, y1, y2;

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
		MemoryPeakResults results = ResultsManager.loadInputResults(inputOption, false, null, null);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			return;
		}

		try
		{
			distanceConverter = results.getDistanceConverter(DistanceUnit.UM);
			timeConverter = results.getTimeConverter(TimeUnit.SECOND);
		}
		catch (Exception e)
		{
			IJ.error(TITLE, "Cannot convert units to um or seconds: " + e.getMessage());
			return;
		}

		// Get the localisation error (4s^2) in raw units^2
		double precision = 0;
		try
		{
			PrecisionResultProcedure p = new PrecisionResultProcedure(results);
			p.getLSEPrecision();

			// Precision in nm using the median
			precision = new Percentile().evaluate(p.precision, 50);
			// Maths.sum(p.precision) / p.precision.length;
			double rawPrecision = distanceConverter.convertBack(precision / 1e3); // Convert from nm to um to raw units
			// Get the localisation error (4s^2) in units^2
			error = 4 * rawPrecision * rawPrecision;
		}
		catch (Exception e)
		{
			Utils.log(TITLE + " - Unable to compute precision: " + e.getMessage());
		}

		// Analyse the track lengths
		results = results.copy();
		results.sort(IdPeakResultComparator.INSTANCE);
		// Ensure the first result triggers an id change
		lastid = results.getFirst().getId() - 1;
		results.forEach(this);
		store(); // For the final track
		d = dList.toArray();
		length = lengthList.toArray();
		id = idList.toArray();
		int[] limits = Maths.limits(length);
		minX = limits[0] - 1;
		maxX = limits[1] + 1;
		h1 = new int[limits[1] + 1];
		h2 = new int[h1.length];
		x1 = createHistogramAxis(h1.length, -WIDTH);
		x2 = createHistogramAxis(h1.length, 0f);
		y1 = new float[x1.length];
		y2 = new float[x1.length];

		// Sort by MSD
		int[] indices = SimpleArrayUtils.newArray(d.length, 0, 1);
		Sort.sortAscending(indices, d, false);
		double[] d2 = d.clone();
		int[] length2 = length.clone();
		int[] id2 = id.clone();
		for (int i = 0; i < indices.length; i++)
		{
			d[i] = d2[indices[i]];
			length[i] = length2[indices[i]];
			id[i] = id2[indices[i]];
		}

		// Interactive analysis
		final NonBlockingExtendedGenericDialog gd = new NonBlockingExtendedGenericDialog(TITLE);
		gd.addMessage(String.format("Split traces into fixed or moving using the track diffusion coefficient (D).\n" +
				"Localistion error has been subtracted from jumps (%s nm).", Utils.rounded(precision)));
		Statistics s = new Statistics(d);
		double av = s.getMean();
		String msg = String.format("Average D per track = %s um^2/s", Utils.rounded(av));
		gd.addMessage(msg);
		// Histogram the diffusion coefficients
		WindowOrganiser wo = new WindowOrganiser();
		int id = Utils.showHistogram("Trace diffusion coefficient", new StoredData(d), "D (um^2/s)", 0, 1, 0, msg);
		if (Utils.isNewWindow())
			wo.add(id);
		double min = Utils.xValues[0];
		double max = Utils.xValues[Utils.xValues.length - 1];
		// see if we can build a nice slider range from the histogram limits
		if (max - min < 5)
			// Because sliders are used when the range is <5 and floating point
			gd.addSlider("D_threshold", min, max, dThreshold);
		else
			gd.addNumericField("D_threshold", dThreshold, 2, 6, "um^2/s");
		gd.addCheckbox("Normalise", normalise);
		gd.addDialogListener(this);
		if (Utils.isShowGenericDialog())
		{
			draw(wo);
			wo.tile();
		}
		gd.setOKLabel("Save datasets");
		gd.setCancelLabel("Close");
		gd.showDialog();

		if (gd.wasCanceled())
			return;

		// Sort by ID
		PeakResult[] list = results.toArray();
		Arrays.sort(list, new IdPeakResultComparator());

		createResults(results, "Fixed", 0, _index, list);
		createResults(results, "Moving", _index, d.length, list);
	}

	private void createResults(MemoryPeakResults results, String suffix, int from, int to, PeakResult[] list)
	{
		MemoryPeakResults out = new MemoryPeakResults();
		out.copySettings(results);
		out.setName(results.getName() + " " + suffix);

		// Sort target ids
		int[] target = Arrays.copyOfRange(id, from, to);
		Arrays.sort(target);

		for (int i = 0, j = 0; i < list.length && j < target.length;)
		{
			int nextId = target[j++];
			// Move forward
			while (i < list.length && list[i].getId() < nextId)
				i++;
			// Write out
			while (i < list.length && list[i].getId() == nextId)
				out.add(list[i++]);
		}

		MemoryPeakResults.addResults(out);
	}

	private boolean showDialog()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addMessage("Analyse the track length of traced data");
		ResultsManager.addInput(gd, "Input", inputOption, InputSource.MEMORY_CLUSTERED);
		gd.showDialog();
		if (gd.wasCanceled())
			return false;
		inputOption = ResultsManager.getInputSource(gd);
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.gui.DialogListener#dialogItemChanged(ij.gui.GenericDialog, java.awt.AWTEvent)
	 */
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
	{
		dThreshold = gd.getNextNumber();
		normalise = gd.getNextBoolean();
		update();
		return true;
	}

	private void draw(WindowOrganiser wo)
	{
		_msdThreshold = dThreshold;
		_normalise = normalise;

		// Find the index in the MSD array
		int index = Arrays.binarySearch(d, _msdThreshold);
		if (index < 0)
			index = -index - 1;
		_index = index;

		// Histogram the distributions
		computeHistogram(0, index, length, h1);
		computeHistogram(index, length.length, length, h2);
		int sum1 = (int) Maths.sum(h1);
		int sum2 = (int) Maths.sum(h2);

		float max1 = createHistogramValues(h1, (_normalise) ? sum1 : 1, y1);
		float max2 = createHistogramValues(h2, (_normalise) ? sum2 : 2, y2);

		String title = "Trace length distribution";
		Plot plot = new Plot(title, "Length", "Frequency");
		plot.setLimits(minX, maxX, 0, Math.max(max1, max2) * 1.05);
		plot.setColor(Color.red);
		plot.addPoints(x1, y1, Plot.LINE);
		plot.setColor(Color.blue);
		plot.addPoints(x2, y2, Plot.LINE);
		plot.setColor(Color.black);
		double p = 100.0 * sum1 / (sum1 + sum2);
		plot.addLabel(0, 0, String.format("Fixed (red) = %d (%s%%), Moving (blue) = %d (%s%%)", sum1, Utils.rounded(p),
				sum2, Utils.rounded(100 - p)));
		PlotWindow pw = Utils.display(title, plot, Utils.NO_TO_FRONT);
		if (wo != null)
		{
			// First call with the window organiser put at the front
			pw.toFront();
			if (Utils.isNewWindow())
				wo.add(pw);
		}
	}

	private static void computeHistogram(int i, int end, int[] length, int[] h)
	{
		Arrays.fill(h, 0);
		while (i < end)
		{
			h[length[i++]]++;
		}
	}

	/**
	 * For the provided histogram x-axis bins, produce an x-axis for plotting. This functions doubles up the histogram
	 * x-positions to allow plotting a square line profile using the ImageJ plot command.
	 * 
	 * @param length
	 * @return
	 */
	public static float[] createHistogramAxis(int length, float offset)
	{
		float[] axis = new float[length * 4];
		int index = 0;
		float offset2 = offset + WIDTH;
		for (int i = 0; i < length; ++i)
		{
			axis[index++] = i + offset;
			axis[index++] = i + offset;
			axis[index++] = i + offset2;
			axis[index++] = i + offset2;
		}
		return axis;
	}

	/**
	 * For the provided histogram y-axis values, produce a y-axis for plotting. This functions doubles up the histogram
	 * values to allow plotting a square line profile using the ImageJ plot command.
	 *
	 * @param histogramY
	 *            the histogram Y
	 * @param norm
	 *            the normalisation
	 * @param axis
	 *            the axis
	 * @return the float
	 */
	public static float createHistogramValues(int[] histogramY, double norm, float[] axis)
	{
		// Assume axis[0] and axis[axis.length-1] == 0
		int index = 1;
		float max = 0;
		for (int i = 0; i < histogramY.length; ++i)
		{
			float v = (float) (histogramY[i] / norm);
			axis[index++] = v;
			axis[index++] = v;
			if (max < v)
				max = v;
			index += 2;
		}
		return max;
	}

	private void update()
	{
		if (lock.acquire())
		{
			// Run in a new thread to allow the GUI to continue updating
			new Thread(new Runnable()
			{
				public void run()
				{
					try
					{
						// Continue while the parameter is changing
						//@formatter:off
							while (
									_msdThreshold != dThreshold ||
									_normalise != normalise
									)
							//@formatter:on
						{
							draw(null);
						}
					}
					finally
					{
						// Ensure the running flag is reset
						lock.release();
					}
				}
			}).start();
		}
	}

	private int lastid = -1;
	private float lastx, lasty;
	private int startFrame, lastFrame, n;
	private double sumSquared;
	private TDoubleArrayList dList = new TDoubleArrayList();
	private TIntArrayList lengthList = new TIntArrayList();
	private TIntArrayList idList = new TIntArrayList();

	public void execute(PeakResult peakResult)
	{
		int id = peakResult.getId();
		int frame = peakResult.getFrame();
		float x = peakResult.getXPosition();
		float y = peakResult.getYPosition();
		if (lastid != id)
		{
			store();
			lastid = id;
			startFrame = frame;
			sumSquared = 0;
			n = 0;
		}
		else
		{
			// Compute the jump
			int jump = frame - lastFrame;
			// Get the raw distance but subtract the expected localisation error
			double d2 = Math.max(0, Maths.distance2(lastx, lasty, x, y) - error);
			// We expect the Mean Squared Distance (MSD) to scale linearly 
			// with time so just weight each jump by the time gap.
			// However we apply a correction factor for diffusion with frames.
			sumSquared += JumpDistanceAnalysis.convertObservedToActual(d2, jump);
			n += jump;
		}
		lastFrame = frame;
		lastx = x;
		lasty = y;
	}

	private void store()
	{
		if (lastid == 0 || n == 0)
			return;
		int length = lastFrame - startFrame;
		double msd = sumSquared / n;
		// 4D = MSD => D = MSD / 4
		// Mean squared distance is in raw units squared.
		// Convert twice since this is a squared distance then from frames to seconds.
		// Use the convertBack() since this is a divide in the final units: um^2/s
		dList.add(timeConverter.convertBack(distanceConverter.convert(distanceConverter.convert(msd))) / 4.0);
		lengthList.add(length);
		idList.add(lastid);
	}
}
