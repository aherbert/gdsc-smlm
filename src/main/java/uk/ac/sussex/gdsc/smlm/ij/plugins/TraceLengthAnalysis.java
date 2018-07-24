/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package uk.ac.sussex.gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.Color;
import java.util.Arrays;

import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import ij.IJ;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.plugin.PlugIn;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.ij.Utils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.NonBlockingExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.utils.Maths;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleLock;
import uk.ac.sussex.gdsc.core.utils.Sort;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.StoredData;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.TimeUnit;
import uk.ac.sussex.gdsc.smlm.fitting.JumpDistanceAnalysis;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.PrecisionResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.sort.IdFramePeakResultComparator;

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
	private final SimpleLock lock = new SimpleLock();
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
	@Override
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
		catch (final Exception e)
		{
			IJ.error(TITLE, "Cannot convert units to um or seconds: " + e.getMessage());
			return;
		}

		// Get the localisation error (4s^2) in raw units^2
		double precision = 0;
		try
		{
			final PrecisionResultProcedure p = new PrecisionResultProcedure(results);
			p.getPrecision();

			// Precision in nm using the median
			precision = new Percentile().evaluate(p.precision, 50);
			// Maths.sum(p.precision) / p.precision.length;
			final double rawPrecision = distanceConverter.convertBack(precision / 1e3); // Convert from nm to um to raw units
			// Get the localisation error (4s^2) in units^2
			error = 4 * rawPrecision * rawPrecision;
		}
		catch (final Exception e)
		{
			Utils.log(TITLE + " - Unable to compute precision: " + e.getMessage());
		}

		// Analyse the track lengths
		results = results.copy();
		results.sort(IdFramePeakResultComparator.INSTANCE);
		// Ensure the first result triggers an id change
		lastid = results.getFirst().getId() - 1;
		results.forEach(this);
		store(); // For the final track
		d = dList.toArray();
		length = lengthList.toArray();
		id = idList.toArray();
		final int[] limits = Maths.limits(length);
		minX = limits[0] - 1;
		maxX = limits[1] + 1;
		h1 = new int[limits[1] + 1];
		h2 = new int[h1.length];
		x1 = createHistogramAxis(h1.length, -WIDTH);
		x2 = createHistogramAxis(h1.length, 0f);
		y1 = new float[x1.length];
		y2 = new float[x1.length];

		// Sort by MSD
		final int[] indices = SimpleArrayUtils.newArray(d.length, 0, 1);
		Sort.sortAscending(indices, d, false);
		final double[] d2 = d.clone();
		final int[] length2 = length.clone();
		final int[] id2 = id.clone();
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
		final Statistics s = new Statistics(d);
		final double av = s.getMean();
		final String msg = String.format("Average D per track = %s um^2/s", Utils.rounded(av));
		gd.addMessage(msg);
		// Histogram the diffusion coefficients
		final WindowOrganiser wo = new WindowOrganiser();
		final int id = Utils.showHistogram("Trace diffusion coefficient", new StoredData(d), "D (um^2/s)", 0, 1, 0, msg);
		if (Utils.isNewWindow())
			wo.add(id);
		final double min = Utils.xValues[0];
		final double max = Utils.xValues[Utils.xValues.length - 1];
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
		final PeakResult[] list = results.toArray();
		Arrays.sort(list, new IdFramePeakResultComparator());

		createResults(results, "Fixed", 0, _index, list);
		createResults(results, "Moving", _index, d.length, list);
	}

	private void createResults(MemoryPeakResults results, String suffix, int from, int to, PeakResult[] list)
	{
		final MemoryPeakResults out = new MemoryPeakResults();
		out.copySettings(results);
		out.setName(results.getName() + " " + suffix);

		// Sort target ids
		Arrays.sort(id, from, to);

		for (int i = 0; i < list.length && from < to;)
		{
			final int nextId = id[from++];
			// Move forward
			while (i < list.length && list[i].getId() < nextId)
				i++;
			// Write out
			while (i < list.length && list[i].getId() == nextId)
				out.add(list[i++]);
		}

		MemoryPeakResults.addResults(out);
	}

	private static boolean showDialog()
	{
		final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
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
	@Override
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
		final int sum1 = (int) Maths.sum(h1);
		final int sum2 = (int) Maths.sum(h2);

		final float max1 = createHistogramValues(h1, (_normalise) ? sum1 : 1, y1);
		final float max2 = createHistogramValues(h2, (_normalise) ? sum2 : 2, y2);

		final String title = "Trace length distribution";
		final Plot plot = new Plot(title, "Length", "Frequency");
		plot.setLimits(minX, maxX, 0, Math.max(max1, max2) * 1.05);
		plot.setColor(Color.red);
		plot.addPoints(x1, y1, Plot.LINE);
		plot.setColor(Color.blue);
		plot.addPoints(x2, y2, Plot.LINE);
		plot.setColor(Color.black);
		final double p = 100.0 * sum1 / (sum1 + sum2);
		plot.addLabel(0, 0, String.format("Fixed (red) = %d (%s%%), Moving (blue) = %d (%s%%)", sum1, Utils.rounded(p),
				sum2, Utils.rounded(100 - p)));
		final PlotWindow pw = Utils.display(title, plot, Utils.NO_TO_FRONT);
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
			h[length[i++]]++;
	}

	/**
	 * For the provided histogram x-axis bins, produce an x-axis for plotting. This functions doubles up the histogram
	 * x-positions to allow plotting a square line profile using the ImageJ plot command.
	 *
	 * @param length
	 *            the length
	 * @param offset
	 *            the offset
	 * @return the x-axis
	 */
	public static float[] createHistogramAxis(int length, float offset)
	{
		final float[] axis = new float[length * 4];
		int index = 0;
		final float offset2 = offset + WIDTH;
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
			final float v = (float) (histogramY[i] / norm);
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
			// Run in a new thread to allow the GUI to continue updating
			new Thread(new Runnable()
			{
				@Override
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
								draw(null);
					}
					finally
					{
						// Ensure the running flag is reset
						lock.release();
					}
				}
			}).start();
	}

	private int lastid = -1;
	private float lastx, lasty;
	private int startFrame, lastFrame, n;
	private double sumSquared;
	private final TDoubleArrayList dList = new TDoubleArrayList();
	private final TIntArrayList lengthList = new TIntArrayList();
	private final TIntArrayList idList = new TIntArrayList();

	@Override
	public void execute(PeakResult peakResult)
	{
		final int id = peakResult.getId();
		final int frame = peakResult.getFrame();
		final float x = peakResult.getXPosition();
		final float y = peakResult.getYPosition();
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
			final int jump = frame - lastFrame;
			// Get the raw distance but subtract the expected localisation error
			final double d2 = Math.max(0, Maths.distance2(lastx, lasty, x, y) - error);
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
		final int length = lastFrame - startFrame;
		final double msd = sumSquared / n;
		// 4D = MSD => D = MSD / 4
		// Mean squared distance is in raw units squared.
		// Convert twice since this is a squared distance then from frames to seconds.
		// Use the convertBack() since this is a divide in the final units: um^2/s
		dList.add(timeConverter.convertBack(distanceConverter.convert(distanceConverter.convert(msd))) / 4.0);
		lengthList.add(length);
		idList.add(lastid);
	}
}
