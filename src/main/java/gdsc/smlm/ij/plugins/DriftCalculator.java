package gdsc.smlm.ij.plugins;

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

import java.awt.Color;
import java.awt.Point;
import java.awt.Rectangle;
import java.util.ArrayList;

import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import java.util.Arrays;
import java.util.Collections;

import gdsc.smlm.fitting.function.Gaussian2DFunction;
import gdsc.smlm.fitting.utils.DoubleEquality;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import ij.util.Tools;

/**
 * Calculates drift using the feducial markers within ROI added to the ROI manager.
 * <p>
 * Adapted from the drift calculation method in QuickPALM.
 */
public class DriftCalculator implements PlugIn
{
	private static String TITLE = "Drift Calculator";

	private static String inputOption = "";
	private static double relativeError = 0.1;
	private static double smoothing = 0.25;
	private static int iterations = 1;
	private static boolean plotDrift = true;
	private static boolean updateResults = false;

	private static PlotWindow plotx = null;
	private static PlotWindow ploty = null;

	private double[] lastdx;
	private double[] lastdy;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		// Require some fit results and selected regions
		if (MemoryPeakResults.countMemorySize() == 0)
		{
			IJ.error(TITLE, "There are no fitting results in memory");
			return;
		}
		Roi[] rois = getRois();
		if (rois == null)
			return;

		if (!showDialog())
			return;

		MemoryPeakResults results = ResultsManager.loadInputResults(inputOption, false);
		exec(results, rois, relativeError, smoothing, iterations);
	}

	private Roi[] getRois()
	{
		RoiManager rmanager = RoiManager.getInstance();
		if (rmanager == null || rmanager.getCount() == 0)
		{
			IJ.error("Add ROIs to use in drift correction to the RoiManager first (select region then press [t]).");
			return null;
		}
		return rmanager.getRoisAsArray();
	}

	private boolean showDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Compute the drift using marked ROIs.");
		ResultsManager.addInput(gd, inputOption, InputSource.Memory);
		gd.addMessage("Stopping criteria");
		gd.addNumericField("Relative_error", relativeError, 3);
		gd.addMessage("loess smoothing parameters");
		gd.addSlider("Smoothing", 0.001, 1, smoothing);
		gd.addSlider("Iterations", 1, 10, iterations);
		gd.addCheckbox("Plot_drift", plotDrift);
		gd.addMessage("Apply drift correction to in-memory results");
		gd.addCheckbox("Update_results", updateResults);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		inputOption = ResultsManager.getInputSource(gd);
		relativeError = gd.getNextNumber();
		smoothing = gd.getNextNumber();
		iterations = (int) gd.getNextNumber();
		plotDrift = gd.getNextBoolean();
		updateResults = gd.getNextBoolean();

		// Check arguments
		try
		{
			Parameters.isAboveZero("Relative error", relativeError);
			Parameters.isPositive("Smoothing", smoothing);
			Parameters.isPositive("Iterations", iterations);
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}
		
		return true;
	}

	private void exec(MemoryPeakResults results, Roi[] rois, double relativeError, double smoothWindow, int iterations)
	{
		int[] limits = findTimeLimits(results);
		Spot[][] roiSpots = findSpots(results, rois);

		// Check we have enough data
		if (roiSpots.length == 0)
		{
			IJ.error("No peak fit results in the selected ROIs");
			return;
		}

		// Convert smoothing window to a fraction of the image size
		double smoothing = 0;
		if (smoothWindow > 0 && smoothWindow <= 1)
		{
			smoothing = smoothWindow;
		}

		double[] dx = new double[limits[1] + 1];
		double[] dy = new double[dx.length];

		double[] sum = new double[dx.length];
		double[] weights = calculateWeights(roiSpots, dx.length, sum);
		lastdx = null;
		double change = calculateDrift(roiSpots, weights, sum, dx, dy, smoothing, iterations);
		if (change == Double.NaN)
			return;
		double error = relativeError + 1;
		IJ.log(String.format("Drift Calculator : Initial drift %g", change));

		for (int i = 1; i <= 100 && error > relativeError && change > 1e-16; i++)
		{
			double lastChange = change;
			change = calculateDrift(roiSpots, weights, sum, dx, dy, smoothing, iterations);
			if (change == Double.NaN)
				return;
			if (i > 1)
			{
				error = DoubleEquality.relativeError(change, lastChange);
				IJ.log(String.format("Iteration %d : Total change %g : Relative change %g", i, change, error));
			}
			else
				IJ.log(String.format("Iteration %d : Total change %g", i, change));
		}

		// Interpolate values for all time limits
		double[][] values = extractValues(weights, limits[0], limits[1], dx, dy);
		PolynomialSplineFunction fx = new SplineInterpolator().interpolate(values[0], values[1]);
		PolynomialSplineFunction fy = new SplineInterpolator().interpolate(values[0], values[2]);

		// Interpolator can only create missing values within the range provided by the input values.
		// The two ends have to be extrapolated.
		// TODO: Perform extrapolation. Currently the end values are used.

		// Find first drift point
		int startT = limits[0], endT = limits[1];
		for (int t = startT; t <= endT; t++)
		{
			if (weights[t] != 0)
			{
				// Use same value to start
				startT = t;
				final double d = dx[t];
				while (t-- > 0)
					dx[t] = d;
				break;
			}
		}
		// Find last drift point
		for (int t = endT; t > startT; t--)
		{
			if (weights[t] != 0)
			{
				// Use same value to end
				endT = t;
				final double d = dx[t];
				while (++t < dx.length)
					dx[t] = d;
				break;
			}
		}

		for (int t = startT; t <= endT; t++)
		{
			if (weights[t] == 0)
			{
				dx[t] = fx.value(t);
				dy[t] = fy.value(t);
			}
		}

		if (plotDrift)
			plotDrift(limits, weights, dx, dy);

		if (updateResults)
			for (PeakResult r : results)
			{
				r.params[Gaussian2DFunction.X_POSITION] += dx[r.peak];
				r.params[Gaussian2DFunction.Y_POSITION] += dy[r.peak];
			}
	}

	private int[] findTimeLimits(MemoryPeakResults results)
	{
		int min = Integer.MAX_VALUE;
		int max = 0;
		for (PeakResult r : results)
		{
			if (min > r.peak)
				min = r.peak;
			if (max < r.peak)
				max = r.peak;
		}
		return new int[] { min, max };
	}

	/**
	 * Build a list of the points that are within each roi
	 * 
	 * @param results
	 * @param rois
	 * @return
	 */
	private Spot[][] findSpots(MemoryPeakResults results, Roi[] rois)
	{
		ArrayList<Spot[]> roiSpots = new ArrayList<Spot[]>(rois.length);
		for (int i = 0; i < rois.length; i++)
		{
			Spot[] spots = findSpots(results, rois[i].getBounds());
			if (spots.length > 0)
				roiSpots.add(spots);
		}
		return roiSpots.toArray(new Spot[roiSpots.size()][]);
	}

	private Spot[] findSpots(MemoryPeakResults results, Rectangle bounds)
	{
		ArrayList<Spot> list = new ArrayList<Spot>(100);
		int maxx = bounds.x + bounds.width;
		int maxy = bounds.y + bounds.height;

		// Find spots within the ROI
		for (PeakResult r : results)
		{
			float x = r.params[Gaussian2DFunction.X_POSITION];
			float y = r.params[Gaussian2DFunction.Y_POSITION];
			if (x > bounds.x && x < maxx && y > bounds.y && y < maxy)
			{
				list.add(new Spot(r.peak, x, y, r.getSignal()));
			}
		}

		// For each frame pick the strongest spot
		Collections.sort(list);

		ArrayList<Spot> newList = new ArrayList<Spot>(list.size());

		int currentT = -1;
		for (Spot spot : list)
		{
			if (currentT != spot.t)
			{
				newList.add(spot);
			}
			currentT = spot.t;
		}

		return newList.toArray(new Spot[newList.size()]);
	}

	/**
	 * For each ROI calculate the sum of the spot intensity. Also compute the sum of the intensity for each time point.
	 * 
	 * @param roiSpots
	 * @param dx
	 *            The total number of timepoints
	 * @param sum
	 *            The sum of the intensity for each ROI
	 * @return The sum of the intensity for each time point.
	 */
	private double[] calculateWeights(Spot[][] roiSpots, int timepoints, double[] sum)
	{
		double[] weights = new double[timepoints];
		for (int i = 0; i < roiSpots.length; i++)
		{
			for (Spot s : roiSpots[i])
			{
				weights[s.t] += s.s;
				sum[i] += s.s;
			}
		}
		return weights;
	}

	/**
	 * Calculate the drift as displacement of each spot from the centre-of-mass. Update the current drift parameters.
	 * 
	 * @param roiSpots
	 * @param dx
	 * @param dy
	 * @param iterations
	 *            Iterations for loess smoothing
	 * @param smoothing
	 *            loess smoothing fraction
	 * @return The total update to the drift parameters (Euclidian distance)
	 */
	private double calculateDrift(Spot[][] roiSpots, double[] weights, double[] sum, double[] dx, double[] dy,
			double smoothing, int iterations)
	{
		double[] newDx = new double[dx.length];
		double[] newDy = new double[dy.length];

		// For each ROI
		for (int i = 0; i < roiSpots.length; i++)
		{
			// Calculate centre-of-mass using the current position (coord + drift)
			double cx = 0, cy = 0;
			for (Spot s : roiSpots[i])
			{
				cx += s.s * (s.x + dx[s.t]);
				cy += s.s * (s.y + dy[s.t]);
			}
			cx /= sum[i];
			cy /= sum[i];

			// Calculate update to the drift as centre-of-mass minus the current position (coord + drift)
			for (Spot s : roiSpots[i])
			{
				newDx[s.t] += s.s * (cx - (s.x + dx[s.t]));
				newDy[s.t] += s.s * (cy - (s.y + dy[s.t]));
			}
		}

		// Normalise
		for (int t = 0; t < dx.length; t++)
		{
			if (weights[t] != 0)
			{
				newDx[t] /= weights[t];
				newDy[t] /= weights[t];
			}

			// New drift = previous drift + update
			newDx[t] += dx[t];
			newDy[t] += dy[t];
		}

		// Store the pure drift values for plotting
		//if (lastdx == null)
		//{
		lastdx = Arrays.copyOf(newDx, newDx.length);
		lastdy = Arrays.copyOf(newDy, newDy.length);
		//}

		// Perform smoothing
		if (smoothing > 0)
		{
			double[][] values = extractValues(weights, 0, dx.length - 1, newDx, newDy);

			// Smooth
			LoessInterpolator loess = new LoessInterpolator(smoothing, iterations);
			values[1] = loess.smooth(values[0], values[1]);
			values[2] = loess.smooth(values[0], values[2]);

			// Add back
			int n = 0;
			for (int t = 0; t < dx.length; t++)
			{
				if (weights[t] != 0)
				{
					newDx[t] = values[1][n];
					newDy[t] = values[2][n];
					n++;
				}
			}

			for (int i = 0; i < newDx.length; i++)
			{
				if (Double.isNaN(newDx[i]))
				{
					IJ.log(String.format("Loess smoothing created bad X-estimate at point %d/%d", i, newDx.length));
					return Double.NaN;
				}
				if (Double.isNaN(newDy[i]))
				{
					IJ.log(String.format("Loess smoothing created bad Y-estimate at point %d/%d", i, newDx.length));
					return Double.NaN;
				}
			}
		}

		// Average drift correction should be zero
		normalise(newDx);
		normalise(newDy);

		// Calculate change and update the input drift parameters
		double change = 0;
		for (int t = 0; t < dx.length; t++)
		{
			if (weights[t] != 0)
			{
				double d1 = dx[t] - newDx[t];
				double d2 = dy[t] - newDy[t];
				change += Math.sqrt(d1 * d1 + d2 * d2);
				dx[t] = newDx[t];
				dy[t] = newDy[t];
			}
		}
		return change;
	}

	private void normalise(double[] data)
	{
		double av1 = 0;
		for (int i = 0; i < data.length; i++)
			av1 += data[i];
		av1 /= data.length;

		for (int i = 0; i < data.length; i++)
			data[i] -= av1;
	}

	/**
	 * For all indices between min and max, if the data array is not zero then add the index and the values from array 1
	 * and 2 to the output.
	 * 
	 * @param data
	 * @param minT
	 * @param maxT
	 * @param array1
	 * @param array2
	 * @return Array of [index][array1][array2]
	 */
	private double[][] extractValues(double[] data, int minT, int maxT, double[] array1, double[] array2)
	{
		// Extract data points for smoothing
		int timepoints = maxT - minT + 1;
		double[][] values = new double[3][timepoints];
		int n = 0;
		for (int t = minT; t <= maxT; t++)
		{
			if (data[t] != 0)
			{
				values[0][n] = t;
				values[1][n] = array1[t];
				values[2][n] = array2[t];
				n++;
			}
		}
		values[0] = Arrays.copyOf(values[0], n);
		values[1] = Arrays.copyOf(values[1], n);
		values[2] = Arrays.copyOf(values[2], n);
		return values;
	}

	private void plotDrift(int[] limits, double[] weights, double[] dx, double[] dy)
	{
		// Build an array of timepoints from the min to the max
		double[] completeT = new double[limits[1] + 1];
		for (int i = limits[0]; i < completeT.length; i++)
			completeT[i] = i;

		// Drift should be centred around zero
		normalise(lastdx);
		normalise(lastdy);

		// Extract the interpolated points and the original drift
		double[][] interpolated = extractValues(completeT, limits[0], limits[1], dx, dy);
		double[][] original = extractValues(weights, limits[0], limits[1], lastdx, lastdy);

		plotx = plotDrift(plotx, null, interpolated, original, "Drift X", 1);
		ploty = plotDrift(ploty, plotx, interpolated, original, "Drift Y", 2);
	}

	private PlotWindow plotDrift(PlotWindow src, PlotWindow parent, double[][] interpolated, double[][] original,
			String name, int index)
	{
		// Store previous location
		Point location = null;
		if (src != null && src.isShowing())
		{
			location = src.getLocation();
			src.close();
		}

		// Create plot
		double[] a = Tools.getMinMax(original[0]);
		double[] b = Tools.getMinMax(original[index]);

		Plot plot = new Plot(name, "Frame", "Drift (px)", (float[]) null, (float[]) null);
		plot.setLimits(a[0], a[1], b[0], b[1]);
		plot.setColor(new Color(0, 0, 155)); // De-saturated blue
		plot.addPoints(original[0], original[index], Plot.CROSS);
		plot.setColor(java.awt.Color.RED);
		plot.addPoints(interpolated[0], interpolated[index], Plot.LINE);
		src = plot.show();

		// Set the location
		if (location != null)
		{
			src.setLocation(location);
		}
		else if (parent != null)
		{
			location = parent.getLocation();
			location.y += parent.getHeight();
			src.setLocation(location);
		}

		return src;
	}

	private class Spot implements Comparable<Spot>
	{
		int t;
		double x;
		double y;
		double s; // signal

		public Spot(int t, double x, double y, double s)
		{
			this.t = t;
			this.x = x;
			this.y = y;
			this.s = s;
		}

		public int compareTo(Spot that)
		{
			// Sort in time order
			if (this.t == that.t)
			{
				// ... then signal
				if (this.s > that.s)
					return -1;
				if (this.s < that.s)
					return 1;
				return 0;
			}
			return this.t - that.t;
		}
	}
}
