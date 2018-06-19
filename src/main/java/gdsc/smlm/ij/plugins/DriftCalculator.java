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
package gdsc.smlm.ij.plugins;

import java.awt.Color;
import java.awt.Point;
import java.awt.Rectangle;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.InputMismatchException;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.NoSuchElementException;
import java.util.Scanner;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.regex.Pattern;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.util.FastMath;

import gdsc.core.ij.AlignImagesFFT;
import gdsc.core.ij.AlignImagesFFT.SubPixelMethod;
import gdsc.core.ij.AlignImagesFFT.WindowMethod;
import gdsc.core.ij.IJTrackProgress;
import gdsc.core.ij.Utils;
import gdsc.core.logging.TrackProgress;
import gdsc.core.utils.Maths;
import gdsc.core.utils.TurboList;
import gdsc.core.utils.UnicodeReader;
import gdsc.smlm.data.config.ResultsProtos.ResultsImageMode;
import gdsc.smlm.data.config.ResultsProtos.ResultsImageType;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.results.IJImagePeakResults;
import gdsc.smlm.ij.results.ImagePeakResultsFactory;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.count.Counter;
import gdsc.smlm.results.procedures.PeakResultProcedure;
import gdsc.smlm.results.procedures.StandardResultProcedure;
import gdsc.smlm.results.procedures.XYRResultProcedure;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.ExtendedGenericDialog;
import ij.gui.Plot;
import ij.gui.Plot2;
import ij.gui.PlotWindow;
import ij.gui.Roi;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import ij.process.Blitter;
import ij.process.FHT;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/**
 * Calculates drift in localisation results. Can use the feducial markers within ROI added to the ROI manager or by
 * aligning N consecutive frames with the overall image.
 */
public class DriftCalculator implements PlugIn
{
	private static String TITLE = "Drift Calculator";

	private static String driftFilename = "";
	private static final String SUB_IMAGE_ALIGNMENT = "Localisation Sub-Images";
	private static final String DRIFT_FILE = "Drift File";
	private static final String STACK_ALIGNMENT = "Reference Stack Alignment";
	private static final String MARKED_ROIS = "Marked ROIs";
	private static String method = "";
	private static String[] UPDATE_METHODS = new String[] { "None", "Update", "New dataset", "New truncated dataset" };
	private static int updateMethod = 0;

	private static String inputOption = "";
	private static int maxIterations = 50;
	private static double relativeError = 0.01;
	private static double smoothing = 0.25;
	private static boolean limitSmoothing = true;
	private static int minSmoothingPoints = 10;
	private static int maxSmoothingPoints = 50;
	private static int iterations = 1;
	private static boolean plotDrift = true;
	private static boolean saveDrift = false;

	// Parameters to control the image alignment algorithm
	private static int frames = 2000;
	private static int minimimLocalisations = 50;
	private static String[] SIZES = new String[] { "128", "256", "512", "1024", "2048" };
	private static String reconstructionSize = SIZES[1];
	private static String stackTitle = "";
	private static int startFrame = 1;
	private static int frameSpacing = 1;
	private static int interpolationMethod = ImageProcessor.BILINEAR;
	private static SubPixelMethod subPixelMethod = AlignImagesFFT.SubPixelMethod.CUBIC;

	private static PlotWindow plotx = null;
	private static PlotWindow ploty = null;

	private int interpolationStart, interpolationEnd;
	private double[] calculatedTimepoints;
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
	 * Align images to the reference initialised in the given aligner
	 */
	private class ImageAligner implements Runnable
	{
		AlignImagesFFT aligner;
		ImageProcessor[] ip;
		int[] t;
		Rectangle alignBounds;
		List<double[]> alignments;
		int from, to;

		public ImageAligner(AlignImagesFFT aligner, ImageProcessor[] ip, int[] t, Rectangle alignBounds,
				List<double[]> alignments, int from, int to)
		{
			this.aligner = aligner;
			this.ip = ip;
			this.t = t;
			this.alignBounds = alignBounds;
			this.alignments = alignments;
			this.from = from;
			this.to = to;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Runnable#run()
		 */
		@Override
		public void run()
		{
			for (int i = from; i < to && i < ip.length; i++)
			{
				incrementProgress();
				// Window method is ignored since the image processor is already an FHT image
				double[] result = aligner.align(ip[i], WindowMethod.TUKEY, alignBounds, subPixelMethod);
				// Create a result for failures
				if (result == null)
					result = new double[] { Double.NaN, Double.NaN, t[i] };
				// Store the time point with the result
				result[2] = t[i];
				alignments.add(result);
			}
		}
	}

	/**
	 * Duplicate and translate images
	 */
	private class ImageTranslator implements Runnable
	{
		ImageProcessor[] images, ip;
		double[] dx, dy;
		int from, to;

		public ImageTranslator(ImageProcessor[] images, ImageProcessor[] ip, double dx[], double dy[], int from, int to)
		{
			this.images = images;
			this.ip = ip;
			this.dx = dx;
			this.dy = dy;
			this.from = from;
			this.to = to;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Runnable#run()
		 */
		@Override
		public void run()
		{
			for (int i = from; i < to && i < ip.length; i++)
			{
				incrementProgress();
				ip[i] = images[i].duplicate();
				if (dx[i] != 0 || dy[i] != 0)
				{
					ip[i].setInterpolationMethod(interpolationMethod);
					ip[i].translate(dx[i], dy[i]);
				}
			}
		}
	}

	/**
	 * Creates an image reconstruction from the provided localisations
	 */
	private class ImageBuilder implements Runnable
	{
		ArrayList<Localisation> localisations;
		ImageProcessor[] images;
		int i;
		Rectangle bounds;
		float scale;
		double[] dx, dy;

		public ImageBuilder(ArrayList<Localisation> localisations, ImageProcessor[] images, int i, Rectangle bounds,
				float scale, double[] dx, double[] dy)
		{
			this.localisations = localisations;
			this.images = images;
			this.i = i;
			this.bounds = bounds;
			this.scale = scale;
			this.dx = dx;
			this.dy = dy;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Runnable#run()
		 */
		@Override
		public void run()
		{
			incrementProgress();
			IJImagePeakResults blockImage = newImage(bounds, scale);
			for (Localisation r : localisations)
			{
				blockImage.add(r.t, (float) (r.x + dx[r.t]), (float) (r.y + dy[r.t]), r.s);
			}
			images[i] = getImage(blockImage);
		}
	}

	/**
	 * Prepare the slices in a stack for image correlation.
	 */
	private class ImageFHTInitialiser implements Runnable
	{
		ImageStack stack;
		ImageProcessor[] images;
		AlignImagesFFT aligner;
		FHT[] fhtImages;
		int from, to;

		public ImageFHTInitialiser(ImageStack stack, ImageProcessor[] images, AlignImagesFFT aligner, FHT[] fhtImages,
				int from, int to)
		{
			this.stack = stack;
			this.images = images;
			this.aligner = aligner;
			this.fhtImages = fhtImages;
			this.from = from;
			this.to = to;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Runnable#run()
		 */
		@Override
		public void run()
		{
			for (int i = from; i < to && i < images.length; i++)
			{
				incrementProgress();
				images[i] = stack.getProcessor(i + 1);
				AlignImagesFFT.applyWindowSeparable(images[i], WindowMethod.TUKEY);
				fhtImages[i] = aligner.transformTarget(images[i], WindowMethod.NONE);
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	@Override
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		// Require some fit results and selected regions
		if (MemoryPeakResults.isMemoryEmpty())
		{
			IJ.error(TITLE, "There are no fitting results in memory");
			return;
		}
		Roi[] rois = getRois();
		String[] stackTitles = createStackImageList();

		if (!showDialog(rois, stackTitles))
			return;

		MemoryPeakResults results = ResultsManager.loadInputResults(inputOption, false, DistanceUnit.PIXEL, null);
		if (results == null || results.size() < 2)
		{
			IJ.error(TITLE, "There are not enough fitting results for drift correction");
			return;
		}
		double[][] drift = null;
		int[] limits = findTimeLimits(results);
		if (method.equals(MARKED_ROIS))
		{
			drift = calculateUsingMarkers(results, limits, rois);
		}
		else if (method.equals(STACK_ALIGNMENT))
		{
			ImageStack stack = showStackDialog(stackTitles);
			if (stack == null)
				return;
			drift = calculateUsingImageStack(stack, limits);
		}
		else if (method.equals(DRIFT_FILE))
		{
			drift = calculateUsingDriftFile(limits);
		}
		else
		{
			if (!showSubImageDialog())
				return;
			drift = calculateUsingFrames(results, limits, Integer.parseInt(reconstructionSize));
		}

		if (drift == null)
			return;

		Utils.log("Drift correction interpolated for frames [%d - %d] of [%d - %d] (%s%%)", interpolationStart,
				interpolationEnd, limits[0], limits[1],
				Utils.rounded((100.0 * (interpolationEnd - interpolationStart + 1)) / (limits[1] - limits[0] + 1)));

		applyDriftCorrection(results, drift);
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

	private boolean showDialog(Roi[] rois, String[] stackTitles)
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Correct the drift in localisation results");
		ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);
		ArrayList<String> methods = new ArrayList<String>(4);
		methods.add(SUB_IMAGE_ALIGNMENT);
		methods.add(DRIFT_FILE);
		if (rois != null)
			methods.add(MARKED_ROIS);
		if (stackTitles != null)
			methods.add(STACK_ALIGNMENT);
		String[] items = methods.toArray(new String[methods.size()]);
		gd.addChoice("Method", items, method);
		gd.addMessage("Stopping criteria");
		gd.addSlider("Max_iterations", 0, 100, maxIterations);
		gd.addNumericField("Relative_error", relativeError, 3);
		gd.addMessage("LOESS smoothing parameters");
		gd.addSlider("Smoothing", 0.001, 1, smoothing);
		gd.addCheckbox("Limit_smoothing", limitSmoothing);
		gd.addSlider("Min_smoothing_points", 5, 50, minSmoothingPoints);
		gd.addSlider("Max_smoothing_points", 5, 50, maxSmoothingPoints);
		gd.addSlider("Smoothing_iterations", 1, 10, iterations);
		gd.addCheckbox("Plot_drift", plotDrift);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		inputOption = ResultsManager.getInputSource(gd);
		method = gd.getNextChoice();
		maxIterations = (int) gd.getNextNumber();
		relativeError = gd.getNextNumber();
		smoothing = gd.getNextNumber();
		limitSmoothing = gd.getNextBoolean();
		minSmoothingPoints = (int) gd.getNextNumber();
		maxSmoothingPoints = (int) gd.getNextNumber();
		iterations = (int) gd.getNextNumber();
		plotDrift = gd.getNextBoolean();

		// Check arguments
		try
		{
			Parameters.isPositive("Max iterations", maxIterations);
			Parameters.isAboveZero("Relative error", relativeError);
			Parameters.isPositive("Smoothing", smoothing);
			if (limitSmoothing)
			{
				Parameters.isEqualOrAbove("Min smoothing points", minSmoothingPoints, 3);
				Parameters.isEqualOrAbove("Max smoothing points", maxSmoothingPoints, 3);
				Parameters.isEqualOrAbove("Max smoothing points", maxSmoothingPoints, minSmoothingPoints);
			}
			Parameters.isEqualOrBelow("Smoothing", smoothing, 1);
			Parameters.isPositive("Smoothing iterations", iterations);
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		return true;
	}

	private boolean showSubImageDialog()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Compute the drift using localisation sub-image alignment");
		gd.addNumericField("Frames", frames, 0);
		gd.addSlider("Minimum_localisations", 10, 50, minimimLocalisations);
		gd.addChoice("FFT size", SIZES, reconstructionSize);
		String[] methods = new String[] { AlignImagesFFT.SubPixelMethod.CUBIC.toString(),
				AlignImagesFFT.SubPixelMethod.GAUSSIAN.toString() };
		gd.addChoice("Sub-pixel_method", methods, subPixelMethod.toString());

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		frames = (int) gd.getNextNumber();
		minimimLocalisations = (int) gd.getNextNumber();
		reconstructionSize = gd.getNextChoice();
		subPixelMethod = (gd.getNextChoiceIndex() == 0) ? AlignImagesFFT.SubPixelMethod.CUBIC
				: AlignImagesFFT.SubPixelMethod.GAUSSIAN;

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

	private ImageStack showStackDialog(String[] stackTitles)
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Compute the drift using a reference stack alignment");

		gd.addChoice("Stack_image", stackTitles, stackTitle);
		gd.addMessage("Frame = previous + spacing");
		gd.addNumericField("Start_frame", startFrame, 0);
		gd.addSlider("Frame_spacing", 1, 20, frameSpacing);
		String[] methods = ImageProcessor.getInterpolationMethods();
		gd.addChoice("Interpolation_method", methods, methods[interpolationMethod]);
		methods = new String[] { AlignImagesFFT.SubPixelMethod.CUBIC.toString(),
				AlignImagesFFT.SubPixelMethod.GAUSSIAN.toString() };
		gd.addChoice("Sub-pixel_method", methods, subPixelMethod.toString());
		gd.showDialog();

		if (gd.wasCanceled())
			return null;

		stackTitle = gd.getNextChoice();
		startFrame = (int) gd.getNextNumber();
		frameSpacing = (int) gd.getNextNumber();
		interpolationMethod = gd.getNextChoiceIndex();
		subPixelMethod = (gd.getNextChoiceIndex() == 0) ? AlignImagesFFT.SubPixelMethod.CUBIC
				: AlignImagesFFT.SubPixelMethod.GAUSSIAN;

		try
		{
			Parameters.isAboveZero("Start frame", startFrame);
			Parameters.isAboveZero("Frame spacing", frameSpacing);
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return null;
		}

		ImagePlus imp = WindowManager.getImage(stackTitle);
		if (imp != null && imp.getStackSize() > 1)
			return imp.getImageStack();
		return null;
	}

	/**
	 * Build a list of suitable stack images
	 * 
	 * @return
	 */
	private String[] createStackImageList()
	{
		int[] idList = WindowManager.getIDList();
		if (idList != null)
		{
			String[] list = new String[idList.length];
			int count = 0;
			for (int id : idList)
			{
				ImagePlus imp = WindowManager.getImage(id);
				if (imp != null && imp.getStackSize() > 1)
				{
					list[count++] = imp.getTitle();
				}
			}
			return Arrays.copyOf(list, count);
		}
		return null;
	}

	private void applyDriftCorrection(MemoryPeakResults results, double[][] drift)
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addMessage("Apply drift correction to in-memory results?");
		gd.addChoice("Update_method", UPDATE_METHODS, UPDATE_METHODS[updateMethod]);
		// Option to save the drift unless it was loaded from file
		if (method != DRIFT_FILE)
			gd.addCheckbox("Save_drift", saveDrift);
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		updateMethod = gd.getNextChoiceIndex();
		if (method != DRIFT_FILE)
		{
			saveDrift = gd.getNextBoolean();
			saveDrift(calculatedTimepoints, lastdx, lastdy);
		}
		if (updateMethod == 0)
			return;

		final double[] dx = drift[0];
		final double[] dy = drift[1];

		// Note: We can use the raw procedure on the results because we requested 
		// the results were in pixels

		if (updateMethod == 1)
		{
			// Update the results in memory
			Utils.log("Applying drift correction to the results set: " + results.getName());
			results.forEach(new PeakResultProcedure()
			{
				@Override
				public void execute(PeakResult r)
				{
					r.setXPosition((float) (r.getXPosition() + dx[r.getFrame()]));
					r.setYPosition((float) (r.getYPosition() + dy[r.getFrame()]));
				}
			});
		}
		else
		{
			// Create a new set of results
			final MemoryPeakResults newResults = new MemoryPeakResults(results.size());
			newResults.copySettings(results);
			newResults.setName(results.getName() + " (Corrected)");
			MemoryPeakResults.addResults(newResults);
			final boolean truncate = updateMethod == 3;
			Utils.log("Creating %sdrift corrected results set: " + newResults.getName(),
					(truncate) ? "truncated " : "");
			results.forEach(new PeakResultProcedure()
			{
				@Override
				public void execute(PeakResult r)
				{
					if (truncate)
					{
						if (r.getFrame() < interpolationStart || r.getFrame() > interpolationEnd)
							return;
					}
					r.setXPosition((float) (r.getXPosition() + dx[r.getFrame()]));
					r.setYPosition((float) (r.getYPosition() + dy[r.getFrame()]));
					newResults.add(r);
				}
			});
		}
	}

	/**
	 * Calculates drift using the feducial markers within ROI.
	 * <p>
	 * Adapted from the drift calculation method in QuickPALM.
	 * 
	 * @param results
	 * @param limits
	 * @param rois
	 * @return the drift { dx[], dy[] }
	 */
	private double[][] calculateUsingMarkers(MemoryPeakResults results, int[] limits, Roi[] rois)
	{
		Spot[][] roiSpots = findSpots(results, rois, limits);

		// Check we have enough data
		if (roiSpots.length == 0)
		{
			IJ.error("No peak fit results in the selected ROIs");
			return null;
		}

		double[] dx = new double[limits[1] + 1];
		double[] dy = new double[dx.length];

		double[] sum = new double[roiSpots.length];
		double[] weights = calculateWeights(roiSpots, dx.length, sum);

		double smoothing = updateSmoothingParameter(weights);

		lastdx = null;
		double change = calculateDriftUsingMarkers(roiSpots, weights, sum, dx, dy, smoothing, iterations);
		if (Double.isNaN(change) || tracker.isEnded())
			return null;
		Utils.log("Drift Calculator : Initial drift " + Utils.rounded(change));

		for (int i = 1; i <= maxIterations; i++)
		{
			change = calculateDriftUsingMarkers(roiSpots, weights, sum, dx, dy, smoothing, iterations);
			if (Double.isNaN(change))
				return null;

			if (converged(i, change, getTotalDrift(dx, dy, weights)))
				break;
		}

		if (tracker.isEnded())
			return null;

		interpolate(dx, dy, weights);

		plotDrift(limits, dx, dy);
		saveDrift(weights, dx, dy);

		return new double[][] { dx, dy };
	}

	/**
	 * Update the smoothing parameter using the upper and lower limits for the number of points to use for smoothing
	 * 
	 * @param data
	 *            The data to be smoothed
	 * @return The updated smoothing parameter
	 */
	private double updateSmoothingParameter(double[] data)
	{
		if (!limitSmoothing)
			return smoothing;

		int n = countNonZeroValues(data);

		int bandwidthInPoints = (int) (smoothing * n);

		// Check the bounds for the smoothing
		int original = bandwidthInPoints;
		if (minSmoothingPoints > 0)
		{
			bandwidthInPoints = FastMath.max(bandwidthInPoints, minSmoothingPoints);
		}
		if (maxSmoothingPoints > 0)
		{
			bandwidthInPoints = FastMath.min(bandwidthInPoints, maxSmoothingPoints);
		}

		double newSmoothing = (double) bandwidthInPoints / n;
		if (original != bandwidthInPoints)
			Utils.log("Updated smoothing parameter for %d data points to %s (%d smoothing points)", n,
					Utils.rounded(newSmoothing), bandwidthInPoints);

		return newSmoothing;
	}

	/**
	 * Count the number of points where the data array is not zero
	 * 
	 * @param data
	 * @return the number of points
	 */
	private int countNonZeroValues(double[] data)
	{
		int n = 0;
		for (double d : data)
		{
			if (d != 0)
				n++;
		}
		return n;
	}

	private double getTotalDrift(double[] dx, double[] dy, double[] originalDriftTimePoints)
	{
		double totalDrift = 0;
		for (int t = 0; t < dx.length; t++)
		{
			if (originalDriftTimePoints[t] != 0)
			{
				totalDrift += Math.sqrt(dx[t] * dx[t] + dy[t] * dy[t]);
			}
		}
		return totalDrift;
	}

	private boolean converged(int iteration, double change, double totalDrift)
	{
		double error = change / totalDrift;
		Utils.log("Iteration %d : Drift %s : Total change %s : Relative change %s", iteration,
				Utils.rounded(totalDrift), Utils.rounded(change), Utils.rounded(error));
		if (error < relativeError || change < 1e-16)
			return true;
		if (tracker.isEnded())
		{
			Utils.log("WARNING : Drift calculation was interrupted");
			return true;
		}
		return false;
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
					Utils.log("ERROR : Loess smoothing created bad X-estimate at point %d/%d", t, newDx.length);
					return false;
				}
				if (Double.isNaN(newDy[t]))
				{
					Utils.log("ERROR : Loess smoothing created bad Y-estimate at point %d/%d", t, newDx.length);
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

		PolynomialSplineFunction fx, fy;
		if (values[0].length < 3)
		{
			fx = new LinearInterpolator().interpolate(values[0], values[1]);
			fy = new LinearInterpolator().interpolate(values[0], values[2]);
		}
		else
		{
			fx = new SplineInterpolator().interpolate(values[0], values[1]);
			fy = new SplineInterpolator().interpolate(values[0], values[2]);
		}

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
		StandardResultProcedure sp = new StandardResultProcedure(results);
		sp.getT();
		return Maths.limits(sp.frame);
	}

	/**
	 * Build a list of the points that are within each roi
	 * 
	 * @param results
	 * @param rois
	 * @param limits
	 * @return
	 */
	private Spot[][] findSpots(MemoryPeakResults results, Roi[] rois, int[] limits)
	{
		ArrayList<Spot[]> roiSpots = new ArrayList<Spot[]>(rois.length);
		for (int i = 0; i < rois.length; i++)
		{
			Spot[] spots = findSpots(results, rois[i].getBounds(), limits);
			if (spots.length > 0)
				roiSpots.add(spots);
		}
		return roiSpots.toArray(new Spot[roiSpots.size()][]);
	}

	private Spot[] findSpots(MemoryPeakResults results, Rectangle bounds, int[] limits)
	{
		final TurboList<Spot> list = new TurboList<Spot>(limits[1] - limits[0] + 1);
		final float minx = bounds.x;
		final float miny = bounds.y;
		final float maxx = bounds.x + bounds.width;
		final float maxy = bounds.y + bounds.height;

		// Find spots within the ROI
		results.forEach(DistanceUnit.PIXEL, new XYRResultProcedure()
		{
			@Override
			public void executeXYR(float x, float y, PeakResult r)
			{
				if (x > minx && x < maxx)
				{
					if (y > miny && y < maxy)
						list.add(new Spot(r.getFrame(), x, y, r.getIntensity()));
				}
			}
		});

		// For each frame pick the strongest spot
		Collections.sort(list);

		TurboList<Spot> newList = new TurboList<Spot>(list.size());

		int currentT = -1;
		for (Spot spot : list)
		{
			if (currentT != spot.t)
			{
				newList.add(spot);
				currentT = spot.t;
			}
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
	private double calculateDriftUsingMarkers(Spot[][] roiSpots, double[] weights, double[] sum, double[] dx,
			double[] dy, double smoothing, int iterations)
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
		calculatedTimepoints = Arrays.copyOf(weights, weights.length);
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

	private void plotDrift(int[] limits, double[] dx, double[] dy)
	{
		if (!plotDrift)
			return;

		// Build an array of timepoints from the min to the max
		double[] completeT = new double[limits[1] + 1];
		for (int i = limits[0]; i < completeT.length; i++)
			completeT[i] = i;

		// Drift should be centred around zero for the calculated points to produce a fair plot
		normalise(lastdx, calculatedTimepoints);
		normalise(lastdy, calculatedTimepoints);

		// Extract the interpolated points and the original drift
		double[][] interpolated = extractValues(dx, limits[0], limits[1], dx, dy);
		double[][] original = extractValues(calculatedTimepoints, limits[0], limits[1], lastdx, lastdy);

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

		Plot2 plot = new Plot2(name, "Frame", "Drift (px)", (float[]) null, (float[]) null);
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
	 * Saves the T,X,Y values to file for all t in the originalDriftTimePoints array which are not zero.
	 * 
	 * @param originalDriftTimePoints
	 * @param dx
	 * @param dy
	 */
	private void saveDrift(double[] originalDriftTimePoints, double[] dx, double[] dy)
	{
		if (!saveDrift)
			return;
		if (!getDriftFilename())
			return;
		BufferedWriter out = null;
		try
		{
			out = new BufferedWriter(new FileWriter(driftFilename));
			out.write("Time\tX\tY\n");
			for (int t = 0; t < dx.length; t++)
			{
				if (originalDriftTimePoints[t] != 0)
				{
					out.write(String.format("%d\t%f\t%f\n", t, dx[t], dy[t]));
				}
			}
			Utils.log("Saved calculated drift to file: " + driftFilename);
		}
		catch (IOException e)
		{
		}
		finally
		{
			try
			{
				out.close();
			}
			catch (IOException e)
			{
			}
		}
	}

	/**
	 * Calculates drift using T,X,Y records read from a file.
	 * 
	 * @param limits
	 * @return the drift { dx[], dy[] }
	 */
	private double[][] calculateUsingDriftFile(int[] limits)
	{
		// Read drift TXY from file
		calculatedTimepoints = new double[limits[1] + 1];
		lastdx = new double[calculatedTimepoints.length];
		lastdy = new double[calculatedTimepoints.length];

		if (!getDriftFilename())
			return null;

		if (readDriftFile(limits) < 2)
		{
			Utils.log("ERROR : Not enough drift points within the time limits %d - %d", limits[0], limits[1]);
			return null;
		}

		double[] dx = Arrays.copyOf(lastdx, lastdx.length);
		double[] dy = Arrays.copyOf(lastdy, lastdy.length);

		double smoothing = updateSmoothingParameter(calculatedTimepoints);

		// Perform smoothing
		if (smoothing > 0)
		{
			if (!smooth(dx, dy, calculatedTimepoints, smoothing, iterations))
				return null;
		}

		// Average drift correction for the calculated points should be zero
		normalise(dx, calculatedTimepoints);
		normalise(dy, calculatedTimepoints);

		interpolate(dx, dy, calculatedTimepoints);

		plotDrift(limits, dx, dy);

		return new double[][] { dx, dy };
	}

	private boolean getDriftFilename()
	{
		String[] path = Utils.decodePath(driftFilename);
		OpenDialog chooser = new OpenDialog("Drift_file", path[0], path[1]);
		if (chooser.getFileName() == null)
			return false;
		driftFilename = chooser.getDirectory() + chooser.getFileName();
		Utils.replaceExtension(driftFilename, "tsv");
		return true;
	}

	/**
	 * Read the drift file storing the T,X,Y into the class level calculatedTimepoints, lastdx and lastdy
	 * arrays. Ignore any records where T is outside the limits.
	 * 
	 * @param limits
	 * @return The number of records read
	 */
	private int readDriftFile(int[] limits)
	{
		int ok = 0;
		BufferedReader input = null;
		try
		{
			FileInputStream fis = new FileInputStream(driftFilename);
			input = new BufferedReader(new UnicodeReader(fis, null));

			String line;
			Pattern pattern = Pattern.compile("[\t, ]+");
			while ((line = input.readLine()) != null)
			{
				if (line.length() == 0)
					continue;
				if (Character.isDigit(line.charAt(0)))
				{
					try
					{
						Scanner scanner = new Scanner(line);
						scanner.useDelimiter(pattern);
						scanner.useLocale(Locale.US);
						final int t = scanner.nextInt();
						if (t < limits[0] || t > limits[1])
							continue;
						final double x = scanner.nextDouble();
						final double y = scanner.nextDouble();
						calculatedTimepoints[t] = ++ok;
						lastdx[t] = x;
						lastdy[t] = y;
						scanner.close();
					}
					catch (InputMismatchException e)
					{
					}
					catch (NoSuchElementException e)
					{
					}
				}
			}
		}
		catch (IOException e)
		{
			// ignore
		}
		finally
		{
			try
			{
				if (input != null)
					input.close();
			}
			catch (IOException e)
			{
				// Ignore
			}
		}
		return ok;
	}

	private class BlockPeakResultProcedure implements PeakResultProcedure
	{
		final ArrayList<ArrayList<Localisation>> blocks = new ArrayList<ArrayList<Localisation>>();
		ArrayList<Localisation> nextBlock = null;
		final Counter counter = new Counter();

		@Override
		public void execute(PeakResult r)
		{
			if (r.getFrame() > counter.getCount())
			{
				while (r.getFrame() > counter.getCount())
					counter.increment(frames);
				// To avoid blocks without many results only create a new block if the min size has been met
				if (nextBlock == null || nextBlock.size() >= minimimLocalisations)
					nextBlock = new ArrayList<Localisation>();
				blocks.add(nextBlock);
			}
			nextBlock.add(new Localisation(r.getFrame(), r.getXPosition(), r.getYPosition(), r.getIntensity()));
		}
	}

	/**
	 * Calculates drift using images from N consecutive frames aligned to the overall image.
	 * 
	 * @param results
	 * @param limits
	 * @param reconstructionSize
	 * @return the drift { dx[], dy[] }
	 */
	private double[][] calculateUsingFrames(MemoryPeakResults results, int[] limits, int reconstructionSize)
	{
		// Extract the localisations into blocks of N consecutive frames
		BlockPeakResultProcedure p = new BlockPeakResultProcedure();
		results.sort();
		results.forEach(new BlockPeakResultProcedure());

		final ArrayList<ArrayList<Localisation>> blocks = p.blocks;
		ArrayList<Localisation> nextBlock = p.nextBlock;

		if (blocks.size() < 2)
		{
			tracker.log("ERROR : Require at least 2 images for drift calculation");
			return null;
		}

		// Check the final block has enough localisations
		if (nextBlock.size() < minimimLocalisations)
		{
			blocks.remove(blocks.size() - 1);
			if (blocks.size() < 2)
			{
				tracker.log("ERROR : Require at least 2 images for drift calculation");
				return null;
			}

			ArrayList<Localisation> combinedBlock = blocks.get(blocks.size() - 1);
			combinedBlock.addAll(nextBlock);
		}

		// Find the average time point for each block
		int[] blockT = new int[blocks.size()];
		int t = 0;
		for (ArrayList<Localisation> block : blocks)
		{
			long sum = 0;
			for (Localisation r : block)
			{
				sum += r.t;
			}
			blockT[t++] = (int) (sum / block.size());
		}

		// Calculate a scale to use when constructing the images for alignment
		Rectangle bounds = results.getBounds(true);
		float scale = (reconstructionSize - 1f) / FastMath.max(bounds.width, bounds.height);

		threadPool = Executors.newFixedThreadPool(Prefs.getThreads());

		double[] dx = new double[limits[1] + 1];
		double[] dy = new double[dx.length];

		double[] originalDriftTimePoints = getOriginalDriftTimePoints(dx, blockT);
		lastdx = null;

		double smoothing = updateSmoothingParameter(originalDriftTimePoints);

		double change = calculateDriftUsingFrames(blocks, blockT, bounds, scale, dx, dy, originalDriftTimePoints,
				smoothing, iterations);
		if (Double.isNaN(change) || tracker.isEnded())
			return null;

		plotDrift(limits, dx, dy);
		Utils.log("Drift Calculator : Initial drift " + Utils.rounded(change));

		for (int i = 1; i <= maxIterations; i++)
		{
			change = calculateDriftUsingFrames(blocks, blockT, bounds, scale, dx, dy, originalDriftTimePoints,
					smoothing, iterations);
			if (Double.isNaN(change))
				return null;

			plotDrift(limits, dx, dy);

			if (converged(i, change, getTotalDrift(dx, dy, originalDriftTimePoints)))
				break;
		}

		if (tracker.isEnded())
			return null;

		plotDrift(limits, dx, dy);

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
	 * @param smoothing
	 * @param iterations
	 * @return
	 */
	private double calculateDriftUsingFrames(ArrayList<ArrayList<Localisation>> blocks, int[] blockT, Rectangle bounds,
			float scale, double[] dx, double[] dy, double[] originalDriftTimePoints, double smoothing, int iterations)
	{
		// Construct images using the current drift
		tracker.status("Constructing images");

		// Built an image for each block of results.
		final ImageProcessor[] images = new ImageProcessor[blocks.size()];

		List<Future<?>> futures = new LinkedList<Future<?>>();
		progressCounter = 0;
		totalCounter = images.length * 2;

		for (int i = 0; i < images.length; i++)
		{
			futures.add(threadPool.submit(new ImageBuilder(blocks.get(i), images, i, bounds, scale, dx, dy)));
		}
		Utils.waitForCompletion(futures);

		for (int i = 0; i < blocks.size(); i++)
		{
			tracker.progress(i, blocks.size());
			IJImagePeakResults blockImage = newImage(bounds, scale);
			for (Localisation r : blocks.get(i))
			{
				blockImage.add(r.t, (float) (r.x + dx[r.t]), (float) (r.y + dy[r.t]), r.s);
			}
			images[i] = getImage(blockImage);
		}

		// Build an image with all results.
		FloatProcessor allIp = new FloatProcessor(images[0].getWidth(), images[0].getHeight());
		for (ImageProcessor ip : images)
			allIp.copyBits(ip, 0, 0, Blitter.ADD);

		return calculateDrift(blockT, scale, dx, dy, originalDriftTimePoints, smoothing, iterations, images, allIp,
				true);
	}

	/**
	 * Calculate the drift of images to the reference image. Update the current drift parameters.
	 * 
	 * @param imageT
	 *            The frame number for each image
	 * @param scale
	 *            The image scale (used to adjust the drift to the correct size)
	 * @param dx
	 *            The X drift
	 * @param dy
	 *            The Y drift
	 * @param originalDriftTimePoints
	 *            Non-zero when the frame number refers to an aligned image frame
	 * @param smoothing
	 *            LOESS smoothing parameter
	 * @param iterations
	 *            LOESS iterations parameter
	 * @param images
	 *            The images to align
	 * @param reference
	 *            The reference image
	 * @param includeCurrentDrift
	 *            Set to true if the input images already have the current drift applied. The new drift will be added to
	 *            the current drift.
	 * @return
	 */
	private double calculateDrift(int[] imageT, float scale, double[] dx, double[] dy, double[] originalDriftTimePoints,
			double smoothing, int iterations, final ImageProcessor[] images, FloatProcessor reference,
			boolean includeCurrentDrift)
	{
		// Align
		tracker.status("Aligning images");
		final AlignImagesFFT aligner = new AlignImagesFFT();
		aligner.init(reference, WindowMethod.NONE, false);
		final Rectangle alignBounds = AlignImagesFFT.createHalfMaxBounds(reference.getWidth(), reference.getHeight(),
				reference.getWidth(), reference.getHeight());

		List<double[]> alignments = Collections.synchronizedList(new ArrayList<double[]>(images.length));
		List<Future<?>> futures = new LinkedList<Future<?>>();

		int imagesPerThread = getImagesPerThread(images);
		for (int i = 0; i < images.length; i += imagesPerThread)
		{
			futures.add(threadPool.submit(
					new ImageAligner(aligner, images, imageT, alignBounds, alignments, i, i + imagesPerThread)));
		}
		Utils.waitForCompletion(futures);
		tracker.progress(1);

		// Used to flag when an alignment has failed
		originalDriftTimePoints = Arrays.copyOf(originalDriftTimePoints, originalDriftTimePoints.length);
		double[] newDx = new double[dx.length];
		double[] newDy = new double[dy.length];
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
				if (includeCurrentDrift)
				{
					// New drift = update + previous drift
					newDx[t] += dx[t];
					newDy[t] += dy[t];
				}
			}
		}

		if (ok < 2)
		{
			tracker.log("ERROR : Unable to align more than 1 image to the overall projection");
			return Double.NaN;
		}

		// Store the pure drift values for plotting
		calculatedTimepoints = Arrays.copyOf(originalDriftTimePoints, originalDriftTimePoints.length);
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
	 * Get the number of images that should be processed on each thread
	 * 
	 * @param images
	 *            The list of images
	 * @return The images per thread
	 */
	private int getImagesPerThread(final ImageProcessor[] images)
	{
		return FastMath.max(1, (int) Math.round((double) images.length / Prefs.getThreads()));
	}

	private IJImagePeakResults newImage(Rectangle bounds, float imageScale)
	{
		IJImagePeakResults image = ImagePeakResultsFactory.createPeakResultsImage(ResultsImageType.DRAW_INTENSITY, true,
				false, "", bounds, 100, 1, imageScale, 0, ResultsImageMode.IMAGE_ADD);
		image.setDisplayImage(false);
		image.begin();
		return image;
	}

	private ImageProcessor getImage(IJImagePeakResults imageResults)
	{
		imageResults.end();
		return imageResults.getImagePlus().getProcessor();
	}

	/**
	 * Calculates drift using images from a reference stack aligned to the overall z-projection.
	 * 
	 * @param stack
	 * 
	 * @param limits
	 * @return the drift { dx[], dy[] }
	 */
	private double[][] calculateUsingImageStack(ImageStack stack, int[] limits)
	{
		// Update the limits using the stack size
		int upperT = startFrame + frameSpacing * (stack.getSize() - 1);
		limits[1] = FastMath.max(limits[1], upperT);

		// TODO - Truncate the stack if there are far too many frames for the localisation limits

		tracker.status("Constructing images");

		threadPool = Executors.newFixedThreadPool(Prefs.getThreads());

		// Built an image and FHT image for each slice
		final ImageProcessor[] images = new ImageProcessor[stack.getSize()];
		final FHT[] fhtImages = new FHT[stack.getSize()];

		List<Future<?>> futures = new LinkedList<Future<?>>();
		progressCounter = 0;
		totalCounter = images.length;

		int imagesPerThread = getImagesPerThread(images);
		final AlignImagesFFT aligner = new AlignImagesFFT();
		FloatProcessor referenceIp = stack.getProcessor(1).toFloat(0, null);
		// We do not care about the window method because this processor will not 
		// actually be used for alignment, it is a reference for the FHT size		
		aligner.init(referenceIp, WindowMethod.NONE, false);
		for (int i = 0; i < images.length; i += imagesPerThread)
		{
			futures.add(threadPool
					.submit(new ImageFHTInitialiser(stack, images, aligner, fhtImages, i, i + imagesPerThread)));
		}
		Utils.waitForCompletion(futures);
		tracker.progress(1);

		if (tracker.isEnded())
			return null;

		double[] dx = new double[limits[1] + 1];
		double[] dy = new double[dx.length];

		double[] originalDriftTimePoints = new double[dx.length];
		int[] blockT = new int[stack.getSize()];
		for (int i = 0, t = startFrame; i < stack.getSize(); i++, t += frameSpacing)
		{
			originalDriftTimePoints[t] = 1;
			blockT[i] = t;
		}

		double smoothing = updateSmoothingParameter(originalDriftTimePoints);

		lastdx = null;
		// For the first iteration calculate drift to the first image in the stack
		// (since the average projection may have a large drift blurring the image)
		double change = calculateDriftUsingImageStack(referenceIp, images, fhtImages, blockT, dx, dy,
				originalDriftTimePoints, smoothing, iterations);
		if (Double.isNaN(change) || tracker.isEnded())
			return null;

		plotDrift(limits, dx, dy);
		Utils.log("Drift Calculator : Initial drift " + Utils.rounded(change));

		for (int i = 1; i <= maxIterations; i++)
		{
			change = calculateDriftUsingImageStack(null, images, fhtImages, blockT, dx, dy, originalDriftTimePoints,
					smoothing, iterations);
			if (Double.isNaN(change))
				return null;

			plotDrift(limits, dx, dy);

			if (converged(i, change, getTotalDrift(dx, dy, originalDriftTimePoints)))
				break;
		}

		if (tracker.isEnded())
			return null;

		plotDrift(limits, dx, dy);

		return new double[][] { dx, dy };
	}

	/**
	 * Calculate the drift of images to the reference image. If no reference is provided then produce a combined
	 * z-projection. Update the current drift parameters.
	 * 
	 * @param reference
	 * @param images
	 *            The images to align
	 * @param fhtImages
	 *            The images to align (pre-transformed to a FHT)
	 * @param blockT
	 *            The frame number for each image
	 * @param dx
	 *            The X drift
	 * @param dy
	 *            The Y drift
	 * @param originalDriftTimePoints
	 *            Non-zero when the frame number refers to an aligned image frame
	 * @param smoothing
	 * @param iterations
	 * @return The change in the drift (NaN is an error occurred)
	 */
	private double calculateDriftUsingImageStack(FloatProcessor reference, ImageProcessor[] images, FHT[] fhtImages,
			int[] blockT, double[] dx, double[] dy, double[] originalDriftTimePoints, double smoothing, int iterations)
	{
		progressCounter = 0;
		totalCounter = images.length;

		if (reference == null)
		{
			// Construct images using the current drift
			tracker.status("Constructing reference image");

			// Built an image using the current drift
			List<Future<?>> futures = new LinkedList<Future<?>>();
			totalCounter = images.length * 2;

			final ImageProcessor[] blockIp = new ImageProcessor[images.length];
			double[] threadDx = new double[images.length];
			double[] threadDy = new double[images.length];
			for (int i = 0; i < images.length; i++)
			{
				threadDx[i] = dx[blockT[i]];
				threadDy[i] = dy[blockT[i]];
			}

			int imagesPerThread = getImagesPerThread(images);
			for (int i = 0; i < images.length; i += imagesPerThread)
			{
				futures.add(threadPool
						.submit(new ImageTranslator(images, blockIp, threadDx, threadDy, i, i + imagesPerThread)));
			}
			Utils.waitForCompletion(futures);

			// Build an image with all results.
			reference = new FloatProcessor(blockIp[0].getWidth(), blockIp[0].getHeight());
			for (ImageProcessor ip : blockIp)
			{
				reference.copyBits(ip, 0, 0, Blitter.ADD);
			}
		}

		// Ensure the reference is windowed
		AlignImagesFFT.applyWindowSeparable(reference, WindowMethod.TUKEY);

		return calculateDrift(blockT, 1f, dx, dy, originalDriftTimePoints, smoothing, iterations, fhtImages, reference,
				false);
	}

	/**
	 * Used to precalculate the localisation signal and store it with T,X,Y values with double precision
	 */
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

		@Override
		public int compareTo(Spot that)
		{
			// Sort in time order
			if (this.t == that.t)
			{
				// ... then signal
				return Double.compare(that.s, this.s);
			}
			return (this.t < that.t) ? -1 : 1;
		}
	}

	/**
	 * Used to precalculate the localisation signal and store it with T,X,Y values
	 */
	private class Localisation
	{
		int t;
		float x;
		float y;
		float s; // signal

		public Localisation(int t, float x, float y, float s)
		{
			this.t = t;
			this.x = x;
			this.y = y;
			this.s = s;
		}
	}
}
