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

import gdsc.smlm.fitting.function.Gaussian2DFunction;
import gdsc.smlm.fitting.utils.DoubleEquality;
import gdsc.smlm.ij.IJTrackProgress;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.results.IJImagePeakResults;
import gdsc.smlm.ij.results.ImagePeakResultsFactory;
import gdsc.smlm.ij.results.ResultsImage;
import gdsc.smlm.ij.results.ResultsMode;
import gdsc.smlm.ij.utils.AlignImagesFFT;
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.TrackProgress;
import gdsc.smlm.utils.Maths;
import ij.IJ;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import ij.process.Blitter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.Color;
import java.awt.Point;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

/**
 * Calculates drift in localisation results. Can use the feducial markers within ROI added to the ROI manager or by
 * aligning N consecutive frames with the overall image.
 */
public class DriftCalculator implements PlugIn
{
	private static String TITLE = "Drift Calculator";
	private static String[] METHODS = new String[] { "Image Alignment", "Marked ROIs" };
	private static int method = 0;

	private static String inputOption = "";
	private static int maxIterations = 10;
	private static double relativeError = 0.1;
	private static double smoothing = 0.25;
	private static int iterations = 1;
	private static boolean plotDrift = true;
	private static boolean updateResults = false;

	// Parameters to control the image alignment algorithm
	private static int frames = 500;
	private static int minimimLocalisations = 50;
	private static String[] SIZES = new String[] { "128", "256", "512", "1024", "2048" };
	private static String reconstructionSize = SIZES[1];

	private static PlotWindow plotx = null;
	private static PlotWindow ploty = null;

	private int interpolationStart, interpolationEnd;
	private double[] lastdx;
	private double[] lastdy;

	private TrackProgress tracker = new IJTrackProgress();

	// Used to multi-thread the image alignment
	private ExecutorService threadPool = null;
	private int progressCounter = 0;
	private int totalCounter = 0;

	private synchronized void incrementProgress()
	{
		tracker.progress(++progressCounter, totalCounter);
	}

	/**
	 * Use a runnable to allow multi-threaded operation. Input parameters
	 * that are manipulated should have synchronized methods.
	 */
	private class ImageWorker implements Runnable
	{
		AlignImagesFFT aligner;
		ImageProcessor ip;
		int t;
		Rectangle alignBounds;
		List<double[]> alignments;

		public ImageWorker(AlignImagesFFT aligner, ImageProcessor ip, int t, Rectangle alignBounds,
				List<double[]> alignments)
		{
			this.aligner = aligner;
			this.ip = ip;
			this.t = t;
			this.alignBounds = alignBounds;
			this.alignments = alignments;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Runnable#run()
		 */
		public void run()
		{
			incrementProgress();
			double[] result = aligner.align(ip, AlignImagesFFT.WindowMethod.Tukey, alignBounds,
					AlignImagesFFT.SubPixelMethod.Cubic);
			// Create a result for failures
			if (result == null)
				result = new double[] { Double.NaN, Double.NaN, t };
			// Store the time point with the result
			result[2] = t;
			alignments.add(result);
		}
	}

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

		if (!showDialog(rois))
			return;

		MemoryPeakResults results = ResultsManager.loadInputResults(inputOption, false);
		if (results.size() < 2)
		{
			IJ.error(TITLE, "There are not enough fitting results for drift correction");
			return;
		}
		double[][] drift = null;
		int[] limits = findTimeLimits(results);
		switch (method)
		{
			case 1:
				drift = calculateUsingMarkers(results, limits, rois, relativeError, smoothing, iterations);
				break;
			default:
				if (!showImageDialog())
					return;
				drift = calculateUsingFrames(results, limits, Integer.parseInt(reconstructionSize));
		}

		if (drift == null)
			return;

		Utils.log("Drift correction interpolated for frames [%d - %d] of [%d - %d] (%s%%)", interpolationStart,
				interpolationEnd, limits[0], limits[1],
				Utils.rounded((100.0 * (interpolationEnd - interpolationStart + 1)) / (limits[1] - limits[0] + 1)));

		GenericDialog gd = new GenericDialog(TITLE);
		gd.addMessage("Apply drift correction to in-memory results?");
		gd.addCheckbox("Update_results", updateResults);
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		if (updateResults = gd.getNextBoolean())
		{
			Utils.log("Applying drift correction");
			final double[] dx = drift[0];
			final double[] dy = drift[1];
			for (PeakResult r : results)
			{
				r.params[Gaussian2DFunction.X_POSITION] += dx[r.peak];
				r.params[Gaussian2DFunction.Y_POSITION] += dy[r.peak];
			}
		}
	}

	private Roi[] getRois()
	{
		RoiManager rmanager = RoiManager.getInstance();
		if (rmanager == null || rmanager.getCount() == 0)
		{
			IJ.log("To use feducial markers for drift correction, add ROIs to the RoiManager (select a region then press [t]).");
			return null;
		}
		return rmanager.getRoisAsArray();
	}

	private boolean showDialog(Roi[] rois)
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Compute the drift in localisation results");
		ResultsManager.addInput(gd, inputOption, InputSource.Memory);
		String[] items = METHODS;
		if (rois == null)
		{
			items = Arrays.copyOf(METHODS, 1);
			method = 0;
		}
		gd.addChoice("Method", items, items[method]);
		gd.addSlider("Max_iterations", 0, 100, maxIterations);
		gd.addMessage("Stopping criteria");
		gd.addNumericField("Relative_error", relativeError, 3);
		gd.addMessage("loess smoothing parameters");
		gd.addSlider("Smoothing", 0.001, 1, smoothing);
		gd.addSlider("Iterations", 1, 10, iterations);
		gd.addCheckbox("Plot_drift", plotDrift);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		inputOption = ResultsManager.getInputSource(gd);
		method = gd.getNextChoiceIndex();
		maxIterations = (int) gd.getNextNumber();
		relativeError = gd.getNextNumber();
		smoothing = gd.getNextNumber();
		iterations = (int) gd.getNextNumber();
		plotDrift = gd.getNextBoolean();

		// Check arguments
		try
		{
			Parameters.isAboveZero("Relative error", relativeError);
			Parameters.isPositive("Smoothing", smoothing);
			Parameters.isEqualOrBelow("Smoothing", smoothing, 1);
			Parameters.isPositive("Iterations", iterations);
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		return true;
	}

	private boolean showImageDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Compute the drift in localisation results using sub-image alignment");
		gd.addNumericField("Frames", frames, 0);
		gd.addSlider("Minimum_localisations", 10, 50, minimimLocalisations);
		gd.addChoice("FFT size", SIZES, reconstructionSize);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		frames = (int) gd.getNextNumber();
		minimimLocalisations = (int) gd.getNextNumber();
		reconstructionSize = gd.getNextChoice();

		// Check arguments
		try
		{
			Parameters.isAboveZero("Frames", frames);
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		return true;
	}

	/**
	 * Calculates drift using the feducial markers within ROI.
	 * <p>
	 * Adapted from the drift calculation method in QuickPALM.
	 * 
	 * @param results
	 * @param limits
	 * @param rois
	 * @param relativeError
	 * @param smoothWindow
	 * @param iterations
	 * @return the drift { dx[], dy[] }
	 */
	private double[][] calculateUsingMarkers(MemoryPeakResults results, int[] limits, Roi[] rois, double relativeError,
			double smoothWindow, int iterations)
	{
		Spot[][] roiSpots = findSpots(results, rois);

		// Check we have enough data
		if (roiSpots.length == 0)
		{
			IJ.error("No peak fit results in the selected ROIs");
			return null;
		}

		double smoothing = 0.2;
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
		if (Double.isNaN(change))
			return null;
		double error = relativeError + 1;
		IJ.log(String.format("Drift Calculator : Initial drift %g", change));

		for (int i = 1; i <= maxIterations && error > relativeError && change > 1e-16; i++)
		{
			double lastChange = change;
			change = calculateDrift(roiSpots, weights, sum, dx, dy, smoothing, iterations);
			if (Double.isNaN(change))
				return null;
			if (i > 1)
			{
				error = DoubleEquality.relativeError(change, lastChange);
				IJ.log(String.format("Iteration %d : Total change %g : Relative change %g", i, change, error));
			}
			else
				IJ.log(String.format("Iteration %d : Total change %g", i, change));
		}

		interpolate(dx, dy, weights);

		if (plotDrift)
			plotDrift(limits, weights, dx, dy);

		return new double[][] { dx, dy };
	}

	private static boolean smooth(double[] newDx, double[] newDy, double[] originalDriftTimePoints, double smoothing,
			int iterations)
	{
		double[][] values = extractValues(originalDriftTimePoints, 0, newDx.length - 1, newDx, newDy);

		// Smooth
		LoessInterpolator loess = new LoessInterpolator(smoothing, iterations);
		values[1] = loess.smooth(values[0], values[1]);
		values[2] = loess.smooth(values[0], values[2]);

		// Add back
		int n = 0;
		for (int t = 0; t < newDx.length; t++)
		{
			if (originalDriftTimePoints[t] != 0)
			{
				newDx[t] = values[1][n];
				newDy[t] = values[2][n];
				n++;

				if (Double.isNaN(newDx[t]))
				{
					IJ.log(String.format("Loess smoothing created bad X-estimate at point %d/%d", t, newDx.length));
					return false;
				}
				if (Double.isNaN(newDy[t]))
				{
					IJ.log(String.format("Loess smoothing created bad Y-estimate at point %d/%d", t, newDx.length));
					return false;
				}
			}
		}

		return true;
	}

	private void interpolate(double[] dx, double[] dy, double[] originalDriftTimePoints)
	{
		// Interpolator can only create missing values within the range provided by the input values.
		// The two ends have to be extrapolated.
		// TODO: Perform extrapolation. Currently the end values are used.

		// Find end points
		int startT = 0;
		while (originalDriftTimePoints[startT] == 0)
			startT++;
		int endT = originalDriftTimePoints.length - 1;
		while (originalDriftTimePoints[endT] == 0)
			endT--;

		// Extrapolate using a constant value
		for (int t = startT; t-- > 0;)
		{
			dx[t] = dx[startT];
			dy[t] = dy[startT];
		}
		for (int t = endT; ++t < dx.length;)
		{
			dx[t] = dx[endT];
			dy[t] = dy[endT];
		}

		double[][] values = extractValues(originalDriftTimePoints, startT, endT, dx, dy);
		PolynomialSplineFunction fx = new SplineInterpolator().interpolate(values[0], values[1]);
		PolynomialSplineFunction fy = new SplineInterpolator().interpolate(values[0], values[2]);

		for (int t = startT; t <= endT; t++)
		{
			if (originalDriftTimePoints[t] == 0)
			{
				dx[t] = fx.value(t);
				dy[t] = fy.value(t);
			}
		}

		this.interpolationStart = startT;
		this.interpolationEnd = endT;
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
		lastdx = Arrays.copyOf(newDx, newDx.length);
		lastdy = Arrays.copyOf(newDy, newDy.length);

		// Perform smoothing
		if (smoothing > 0)
		{
			if (!smooth(newDx, newDy, weights, smoothing, iterations))
				return Double.NaN;
		}

		// Average drift correction for the calculated points should be zero to allow change comparison
		normalise(newDx, weights);
		normalise(newDy, weights);

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
			// When weights == 0 any shift we calculate is later ignored since the points are interpolated
			//dx[t] = newDx[t];
			//dy[t] = newDy[t];
		}
		return change;
	}

	/**
	 * Normalise the data so that the points identified by non-zeros in the toProcess array have a centre of mass of
	 * zero. The shift is calculated on a subset of the points but applied to all points.
	 * 
	 * @param data
	 * @param toProcess
	 */
	private void normalise(double[] data, double[] toProcess)
	{
		double av1 = 0;
		int count = 0;
		for (int i = 0; i < data.length; i++)
		{
			if (toProcess[i] != 0)
			{
				av1 += data[i];
				count++;
			}

		}
		av1 /= count;

		for (int i = 0; i < data.length; i++)
			//if (toProcess[i] != 0)
			data[i] -= av1;
	}

	@SuppressWarnings("unused")
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
	private static double[][] extractValues(double[] data, int minT, int maxT, double[] array1, double[] array2)
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

	private void plotDrift(int[] limits, double[] originalDriftTimePoints, double[] dx, double[] dy)
	{
		// Build an array of timepoints from the min to the max
		double[] completeT = new double[limits[1] + 1];
		for (int i = limits[0]; i < completeT.length; i++)
			completeT[i] = i;

		// Drift should be centred around zero for the calculated points to produce a fair plot
		normalise(lastdx, originalDriftTimePoints);
		normalise(lastdy, originalDriftTimePoints);

		// Extract the interpolated points and the original drift
		double[][] interpolated = extractValues(dx, limits[0], limits[1], dx, dy);
		double[][] original = extractValues(originalDriftTimePoints, limits[0], limits[1], lastdx, lastdy);

		plotx = plotDrift(plotx, null, interpolated, original, "Drift X", 1);
		ploty = plotDrift(ploty, plotx, interpolated, original, "Drift Y", 2);
	}

	private PlotWindow plotDrift(PlotWindow src, PlotWindow parent, double[][] interpolated, double[][] original,
			String name, int index)
	{
		// Create plot
		double[] a = Maths.limits(interpolated[0]);
		double[] b = Maths.limits(original[index]);
		b = Maths.limits(b, interpolated[index]);

		Plot plot = new Plot(name, "Frame", "Drift (px)", (float[]) null, (float[]) null);
		plot.setLimits(a[0], a[1], b[0], b[1]);
		plot.setColor(new Color(0, 0, 155)); // De-saturated blue
		plot.addPoints(original[0], original[index], Plot.CROSS);
		plot.setColor(java.awt.Color.RED);
		plot.addPoints(interpolated[0], interpolated[index], Plot.LINE);
		src = Utils.display(name, plot);

		if (Utils.isNewWindow() && parent != null)
		{
			Point location = parent.getLocation();
			location.y += parent.getHeight();
			src.setLocation(location);
		}

		return src;
	}

	/**
	 * Calculates drift using images from N consecutive frames aligned to the overall image.
	 * 
	 * @param results
	 * @params limits
	 * @param reconstructionSize
	 * @return the drift { dx[], dy[] }
	 */
	private double[][] calculateUsingFrames(MemoryPeakResults results, int[] limits, int reconstructionSize)
	{
		double[] dx = new double[limits[1] + 1];
		double[] dy = new double[dx.length];

		// Extract the localisations into blocks of N consecutive frames
		ArrayList<ArrayList<PeakResult>> blocks = new ArrayList<ArrayList<PeakResult>>();
		results.sort();
		List<PeakResult> peakResults = results.getResults();
		int t = 0;
		ArrayList<PeakResult> nextBlock = null;
		for (PeakResult r : peakResults)
		{
			if (r.peak > t)
			{
				while (r.peak > t)
					t += frames;
				// To avoid blocks without many results only create a new block if the min size has been met
				if (nextBlock == null || nextBlock.size() >= minimimLocalisations)
					nextBlock = new ArrayList<PeakResult>();
				blocks.add(nextBlock);
			}
			nextBlock.add(r);
		}

		if (blocks.size() < 2)
		{
			tracker.log("ERROR : Require at least 2 images for drift calculation");
			return null;
		}

		// Check the final block has enough localisations
		if (nextBlock.size() < minimimLocalisations)
		{
			blocks.remove(blocks.size() - 1);
			ArrayList<PeakResult> combinedBlock = blocks.get(blocks.size() - 1);
			combinedBlock.addAll(nextBlock);

			if (blocks.size() < 2)
			{
				tracker.log("ERROR : Require at least 2 images for drift calculation");
				return null;
			}
		}

		// Find the average time point for each block
		int[] blockT = new int[blocks.size()];
		t = 0;
		for (ArrayList<PeakResult> block : blocks)
		{
			long sum = 0;
			for (PeakResult r : block)
			{
				sum += r.peak;
			}
			blockT[t++] = (int) (sum / block.size());
		}

		// Calculate a scale to use when constructing the images for alignment
		Rectangle bounds = results.getBounds(true);
		float scale = (reconstructionSize - 1f) / Math.max(bounds.width, bounds.height);

		threadPool = Executors.newFixedThreadPool(Prefs.getThreads());

		double[] originalDriftTimePoints = getOriginalDriftTimePoints(dx, blockT);
		lastdx = null;
		double change = calculateDrift(blocks, blockT, bounds, scale, dx, dy, originalDriftTimePoints, smoothing,
				iterations);
		if (Double.isNaN(change))
			return null;
		if (plotDrift)
			plotDrift(limits, originalDriftTimePoints, dx, dy);
		double error = relativeError + 1;
		IJ.log(String.format("Drift Calculator : Initial drift %g", change));

		for (int i = 1; i <= maxIterations && error > relativeError && change > 1e-16; i++)
		{
			double lastChange = change;
			change = calculateDrift(blocks, blockT, bounds, scale, dx, dy, originalDriftTimePoints, smoothing,
					iterations);
			if (Double.isNaN(change))
				return null;
			if (plotDrift)
				plotDrift(limits, originalDriftTimePoints, dx, dy);
			if (i > 1)
			{
				error = DoubleEquality.relativeError(change, lastChange);
				IJ.log(String.format("Iteration %d : Total change %g : Relative change %g", i, change, error));
			}
			else
				IJ.log(String.format("Iteration %d : Total change %g", i, change));
		}

		if (plotDrift)
			plotDrift(limits, originalDriftTimePoints, dx, dy);

		return new double[][] { dx, dy };
	}

	/**
	 * Create an array to show the time-point of the original calculated drift alignment
	 * 
	 * @param dx
	 *            The drift array
	 * @param timepoints
	 *            Array of timepoints for which there is a drift calculation
	 * @return array matching dx length with non-zero values for each identified timepoint
	 */
	private double[] getOriginalDriftTimePoints(double[] dx, int[] timepoints)
	{
		double[] originalDriftTimePoints = new double[dx.length];
		for (int i = 0; i < timepoints.length; i++)
			originalDriftTimePoints[timepoints[i]] = 1;
		return originalDriftTimePoints;
	}

	/**
	 * Calculate the drift by aligning N consecutive frames with the overall image. Update the current drift parameters.
	 * 
	 * @param blocks
	 * @param blockT
	 * @param bounds
	 * @param scale
	 * @param dx
	 * @param dy
	 * @param smoothing2
	 * @param iterations2
	 * @return
	 */
	private double calculateDrift(ArrayList<ArrayList<PeakResult>> blocks, int[] blockT, Rectangle bounds, float scale,
			double[] dx, double[] dy, double[] originalDriftTimePoints, double smoothing, int iterations)
	{
		double[] newDx = new double[dx.length];
		double[] newDy = new double[dy.length];

		// Construct images using the current drift
		tracker.status("Constructing images");

		// Built an image for each block of results.
		final ImageProcessor[] blockIp = new ImageProcessor[blocks.size()];
		for (int i = 0; i < blocks.size(); i++)
		{
			tracker.progress(i, blocks.size());
			IJImagePeakResults blockImage = newImage(bounds, scale);
			for (PeakResult r : blocks.get(i))
			{
				float[] params = Arrays.copyOf(r.params, 7);
				params[Gaussian2DFunction.X_POSITION] += dx[r.peak];
				params[Gaussian2DFunction.Y_POSITION] += dy[r.peak];
				blockImage.add(r.peak, 0, 0, 0f, 0d, 0f, params, null);
			}
			blockIp[i] = getImage(blockImage);
		}
		tracker.progress(1);

		// Build an image with all results.
		FloatProcessor allIp = new FloatProcessor(blockIp[0].getWidth(), blockIp[0].getHeight());
		for (ImageProcessor ip : blockIp)
			allIp.copyBits(ip, 0, 0, Blitter.ADD);

		// Align
		tracker.status("Aligning images");
		final AlignImagesFFT aligner = new AlignImagesFFT();
		aligner.init(allIp, AlignImagesFFT.WindowMethod.None, false);
		final Rectangle alignBounds = AlignImagesFFT.createHalfMaxBounds(allIp.getWidth(), allIp.getHeight(),
				allIp.getWidth(), allIp.getHeight());

		List<double[]> alignments = Collections.synchronizedList(new ArrayList<double[]>(blocks.size()));
		List<Future<?>> futures = new LinkedList<Future<?>>();
		progressCounter = 0;
		totalCounter = blocks.size();

		for (int i = 0; i < blocks.size(); i++)
		{
			futures.add(threadPool.submit(new ImageWorker(aligner, blockIp[i], blockT[i], alignBounds, alignments)));
		}
		waitForCompletion(futures);
		tracker.progress(1);

		// Used to flag when an alignment has failed
		originalDriftTimePoints = Arrays.copyOf(originalDriftTimePoints, originalDriftTimePoints.length);
		int ok = 0;
		for (double[] result : alignments)
		{
			int t = (int) result[2];
			if (Double.isNaN(result[0]))
			{
				// TODO: How to ignore bad alignments? 
				// Only do smoothing where there was an alignment?
				originalDriftTimePoints[t] = 0;
				tracker.log("WARNING : Unable to align image for time %d to the overall projection", t);
			}
			else
			{
				ok++;
				newDx[t] = result[0] / scale;
				newDy[t] = result[1] / scale;
				// New drift = update + previous drift
				newDx[t] += dx[t];
				newDy[t] += dy[t];
			}
		}

		if (ok < 2)
		{
			tracker.log("ERROR : Unable to align more than 1 image to the overall projection");
			return Double.NaN;
		}

		// Store the pure drift values for plotting
		lastdx = Arrays.copyOf(newDx, newDx.length);
		lastdy = Arrays.copyOf(newDy, newDy.length);

		// Perform smoothing
		if (smoothing > 0)
		{
			tracker.status("Smoothing drift");
			if (!smooth(newDx, newDy, originalDriftTimePoints, smoothing, iterations))
				return Double.NaN;
		}

		// Interpolate values for all time limits
		tracker.status("Interpolating drift");
		interpolate(newDx, newDy, originalDriftTimePoints);

		// Average drift correction for the calculated points should be zero to allow change comparison
		normalise(newDx, originalDriftTimePoints);
		normalise(newDy, originalDriftTimePoints);

		// Calculate change and update the input drift parameters
		double change = 0;
		for (int t = 0; t < dx.length; t++)
		{
			if (originalDriftTimePoints[t] != 0)
			{
				double d1 = dx[t] - newDx[t];
				double d2 = dy[t] - newDy[t];
				change += Math.sqrt(d1 * d1 + d2 * d2);
			}
			// Update all points since interpolation has already been done
			dx[t] = newDx[t];
			dy[t] = newDy[t];
		}

		tracker.status("");
		return change;
	}

	/**
	 * Waits for all threads to complete computation.
	 * 
	 * @param futures
	 */
	private void waitForCompletion(List<Future<?>> futures)
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

	private IJImagePeakResults newImage(Rectangle bounds, float imageScale)
	{
		IJImagePeakResults image = ImagePeakResultsFactory.createPeakResultsImage(ResultsImage.SIGNAL_INTENSITY, true,
				false, "", bounds, 100, 1, imageScale, 0, ResultsMode.ADD);
		image.setDisplayImage(false);
		image.begin();
		return image;
	}

	private ImageProcessor getImage(IJImagePeakResults imageResults)
	{
		imageResults.end();
		return imageResults.getImagePlus().getProcessor();
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
