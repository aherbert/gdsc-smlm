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
import java.awt.Rectangle;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.math3.random.HaltonSequenceGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import gdsc.core.clustering.DensityManager;
import gdsc.core.ij.Utils;
import gdsc.smlm.data.config.ResultsProtos.ResultsImageMode;
import gdsc.smlm.data.config.ResultsProtos.ResultsImageType;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.results.IJImagePeakResults;
import gdsc.smlm.ij.results.ImagePeakResultsFactory;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.procedures.StandardResultProcedure;
import gdsc.smlm.results.procedures.XYRResultProcedure;
import gnu.trove.list.array.TDoubleArrayList;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.ExtendedGenericDialog;
import ij.gui.Plot2;
import ij.plugin.PlugIn;
import ij.plugin.frame.Recorder;

/**
 * Produces an image on localisation using their density
 */
public class DensityImage implements PlugIn
{
	private static String TITLE = "Density Image";
	private static String inputOption = "";
	private static float radius = 1.5f;
	private static boolean chooseRoi = false;
	private static String roiImage = "";
	private static boolean adjustForBorder = true;
	private static int imageScale = 2;
	private static boolean cumulativeImage = false;
	private static boolean useSquareApproximation = false;
	private static int resolution = 10;

	private static String[] ScoreMethods = new String[] { "Density", "Ripley's K", "Ripley's K / Area", "Ripley's L",
			"Ripley's L - r", "Ripley's L / r", "Ripley's (L - r) / r" };
	private static int scoreMethodIndex = 0;

	private static boolean filterLocalisations = true;
	private static double filterThreshold = 0;
	private static boolean computeRipleysPlot = false;

	private static double minR = 0.2;
	private static double maxR = 3;
	private static double incrementR = 0.2;
	private static boolean confidenceIntervals = false;

	private Rectangle roiBounds;
	private double scaledRoiMinX, scaledRoiMaxX, scaledRoiMinY, scaledRoiMaxY;
	private int roiImageWidth, roiImageHeight;

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
		int size = MemoryPeakResults.countMemorySize();
		if (size == 0)
		{
			IJ.error(TITLE, "There are no fitting results in memory");
			return;
		}

		if (!showDialog())
			return;

		MemoryPeakResults results = ResultsManager.loadInputResults(inputOption, false, DistanceUnit.PIXEL, null);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			IJ.showStatus("");
			return;
		}

		boolean[] isWithin = new boolean[1];
		results = cropWithBorder(results, isWithin);
		if (results.size() == 0)
		{
			IJ.error(TITLE, "No results within the crop region");
			IJ.showStatus("");
			return;
		}

		long start = System.currentTimeMillis();
		IJ.showStatus("Calculating density ...");

		boolean useAdjustment = adjustForBorder && !isWithin[0];

		DensityManager dm = createDensityManager(results);
		int[] density = null;
		if (useSquareApproximation)
			density = dm.calculateSquareDensity(radius, resolution, useAdjustment);
		else
			density = dm.calculateDensity(radius, useAdjustment);

		density = cropBorder(results, density);

		// Convert to float
		ScoreCalculator calc = createCalculator(results);
		float[] densityScore = calc.calculate(density);

		int filtered = plotResults(results, densityScore, calc);

		logDensityResults(results, density, radius, filtered);

		if (computeRipleysPlot)
			computeRipleysPlot(results);

		double seconds = (System.currentTimeMillis() - start) / 1000.0;
		IJ.showStatus(TITLE + " complete : " + seconds + "s");
	}

	private static DensityManager createDensityManager(MemoryPeakResults results)
	{
		if (results == null || results.size() == 0)
			throw new IllegalArgumentException("Results are null or empty");

		StandardResultProcedure sp = new StandardResultProcedure(results, DistanceUnit.PIXEL);
		sp.getXY();
		return new DensityManager(sp.x, sp.y, results.getBounds());
	}

	private ScoreCalculator createCalculator(MemoryPeakResults results)
	{
		switch (scoreMethodIndex)
		{
			case 1:
				// Ripley's K (Density / av. density)
				return new KScoreCalculator(results, 0);

			case 2:
				// Ripley's K / area
				return new KScoreCalculator(results, 1);

			case 3:
				// Ripley's L
				return new LScoreCalculator(results, 0);

			case 4:
				// Ripley's L - r
				return new LScoreCalculator(results, 1);

			case 5:
				// Ripley's L / r
				return new LScoreCalculator(results, 2);

			case 6:
				// Ripley's (L - r) / r
				return new LScoreCalculator(results, 3);

			case 0:
			default:
				return new DensityScoreCalculator(results);
		}
	}

	private interface ScoreCalculator
	{
		/**
		 * Get the density score for the input density counts.
		 *
		 * @param density
		 *            the density
		 * @return the float[]
		 */
		float[] calculate(int[] density);

		/**
		 * Get the score threshold for filtering results using the configured filter threshold.
		 *
		 * @return the threshold
		 */
		float getThreshold();
	}

	private class DensityScoreCalculator implements ScoreCalculator
	{
		MemoryPeakResults results;

		public DensityScoreCalculator(MemoryPeakResults results)
		{
			this.results = results;
		}

		@Override
		public float[] calculate(int[] density)
		{
			float[] score = new float[density.length];
			for (int i = 0; i < score.length; i++)
				score[i] = density[i];
			return score;
		}

		protected float getAverageDensity()
		{
			Rectangle bounds = results.getBounds();
			float area = bounds.width * bounds.height;
			return results.size() / area;
		}

		protected float getRegionArea()
		{
			return radius * radius * ((useSquareApproximation) ? 4 : (float) Math.PI);
		}

		@Override
		public float getThreshold()
		{
			float expected = getAverageDensity() * getRegionArea();
			return (float) (expected * filterThreshold);
		}
	}

	private class KScoreCalculator extends DensityScoreCalculator
	{
		int mode;

		public KScoreCalculator(MemoryPeakResults results, int mode)
		{
			super(results);
			this.mode = mode;
		}

		@Override
		public float[] calculate(int[] density)
		{
			float[] score = new float[density.length];
			// K(r)
			float regionDivisor = getAverageDensity();
			if (mode == 1)
				// K(r) / area
				regionDivisor *= getRegionArea();
			for (int i = 0; i < score.length; i++)
			{
				score[i] = density[i] / regionDivisor;
			}
			return score;
		}

		@Override
		public float getThreshold()
		{
			// Note: K(r) ~ Area
			// Since K(r) should be equal to the area to make the filter threshold scale appropriately
			// we adjust the threshold by the area
			if (mode == 0)
				return (float) filterThreshold * getRegionArea();
			// K(r) / area == 1
			// => no adjustment as this is radius scale independent
			return (float) filterThreshold;
		}
	}

	private class LScoreCalculator extends KScoreCalculator
	{
		public LScoreCalculator(MemoryPeakResults results, int mode)
		{
			super(results, mode);
		}

		@Override
		public float[] calculate(int[] density)
		{
			// Compute a normalised variance stabilised per particle L-score.
			// As in: Scarselli, et al.
			// Cell type-specific β2-adrenergic receptor clusters identified using PALM
			// microscopy are not lipid raft related, but depend on actin cytoskeleton integrity.
			// J Biol Chem. 2012 May 11;287(20):16768-80
			// Note:
			// I have re-arranged the score to be:
			//   Li(r) = Math.sqrt((Sample density / Average density) / pi) - r
			// This should be above zero if the density around the spot is higher than the average sample density.

			float[] score = new float[density.length];
			float regionDivisor = getAverageDensity() * ((useSquareApproximation) ? 4 : (float) Math.PI);
			for (int i = 0; i < score.length; i++)
			{
				// L(r)
				score[i] = (float) Math.sqrt(density[i] / regionDivisor);
			}
			if (mode == 1 || mode == 3)
			{
				// L(r) - r
				// (L(r) - r) / r
				for (int i = 0; i < score.length; i++)
					score[i] -= radius;
			}
			if (mode == 2 || mode == 3)
			{
				// L(r) / r
				// (L(r) - r) / r
				for (int i = 0; i < score.length; i++)
					score[i] /= radius;
			}
			return score;
		}

		@Override
		public float getThreshold()
		{
			// Note:
			// L(r) is proportional to radius
			// K(r) is proportional to area
			// To make the filtered results the same to the K(r) function we could use the
			// sqrt of the filterThreshold

			double threshold = filterThreshold;
			//double threshold = Math.sqrt(filterThreshold);

			// Note: L(r) ~ r
			// Since L(r) should be equal to the radius to make the filter threshold scale appropriately
			// we adjust the threshold by the radius
			if (mode == 0)
				return (float) threshold * radius;
			// L(r) - r == 0
			if (mode == 1)
				return (float) threshold * radius - radius;
			// L(r) - r  / r == 0
			if (mode == 3)
				return (float) (threshold * radius - radius) / radius;
			// L(r) / r == 1
			// => no adjustment as this is radius scale independent
			return (float) threshold;
		}
	}

	/**
	 * Crop the results to the ROI. Add a border using the sampling radius so that counts do not have to be approximated
	 * (i.e. all spots around the edge of the ROI will have a complete image to sample from). The results are modified
	 * in place.
	 *
	 * @param results
	 * @param isWithin
	 *            Set to true if the added border is within the original bounds (i.e. no adjustment for missing counts
	 *            is required)
	 * @return
	 */
	private MemoryPeakResults cropWithBorder(MemoryPeakResults results, boolean[] isWithin)
	{
		isWithin[0] = false;
		if (roiBounds == null)
			return results;

		// Adjust bounds relative to input results image:
		// Use the ROI relative to the frame the ROI is drawn on.
		// Map those fractional coordinates back to the original data bounds.
		Rectangle bounds = results.getBounds();
		double xscale = (double) roiImageWidth / bounds.width;
		double yscale = (double) roiImageHeight / bounds.height;

		// Compute relative to the results bounds (if present)
		scaledRoiMinX = bounds.x + roiBounds.x / xscale;
		scaledRoiMaxX = scaledRoiMinX + roiBounds.width / xscale;
		scaledRoiMinY = bounds.y + roiBounds.y / yscale;
		scaledRoiMaxY = scaledRoiMinY + roiBounds.height / yscale;

		// Allow for the border
		final float minX = (int) (scaledRoiMinX - radius);
		final float maxX = (int) Math.ceil(scaledRoiMaxX + radius);
		final float minY = (int) (scaledRoiMinY - radius);
		final float maxY = (int) Math.ceil(scaledRoiMaxY + radius);

		// Create a new set of results within the bounds
		final MemoryPeakResults newResults = new MemoryPeakResults();
		newResults.begin();
		results.forEach(DistanceUnit.PIXEL, new XYRResultProcedure()
		{
			@Override
			public void executeXYR(float x, float y, PeakResult result)
			{
				if (x >= minX && x <= maxX && y >= minY && y <= maxY)
					newResults.add(result);
			}
		});
		newResults.end();
		newResults.copySettings(results);
		newResults.setBounds(new Rectangle((int) minX, (int) minY, (int) (maxX - minX), (int) (maxY - minY)));
		isWithin[0] = minX >= bounds.x && minY >= bounds.y && maxX <= (bounds.x + bounds.width) &&
				maxY <= (bounds.y + bounds.height);
		return newResults;
	}

	/**
	 * Remove any results which fall in the radius border added around the ROI. Results are modified in place and a new
	 * density array is returned.
	 *
	 * @param results
	 * @param density
	 * @return
	 */
	private int[] cropBorder(MemoryPeakResults results, int[] density)
	{
		if (roiBounds == null)
			return density;
		final float minX = (int) (scaledRoiMinX);
		final float maxX = (int) Math.ceil(scaledRoiMaxX);
		final float minY = (int) (scaledRoiMinY);
		final float maxY = (int) Math.ceil(scaledRoiMaxY);
		// Clone the results then add back those that are within the bounds
		PeakResult[] peakResults = results.toArray();
		results.begin();
		int count = 0;
		for (int i = 0; i < peakResults.length; i++)
		{
			PeakResult peakResult = peakResults[i];
			float x = peakResult.getXPosition();
			float y = peakResult.getYPosition();
			if (x < minX || x > maxX || y < minY || y > maxY)
				continue;
			results.add(peakResult);
			density[count++] = density[i];
		}
		results.end();
		results.setBounds(new Rectangle((int) minX, (int) minY, (int) (maxX - minX), (int) (maxY - minY)));
		return Arrays.copyOf(density, count);
	}

	/**
	 * Output a log message of the results including the average density for localisations and the expected average.
	 *
	 * @param results
	 *            the results
	 * @param density
	 *            the density
	 * @param radius
	 *            the radius
	 * @param filtered
	 *            the filtered
	 * @return the summary statistics
	 */
	private static SummaryStatistics logDensityResults(MemoryPeakResults results, int[] density, float radius,
			int filtered)
	{
		float region = (float) (radius * radius * ((useSquareApproximation) ? 4 : Math.PI));

		Rectangle bounds = results.getBounds();
		float area = bounds.width * bounds.height;
		float expected = results.size() * region / area;
		SummaryStatistics summary = new SummaryStatistics();

		for (int i = 0; i < results.size(); i++)
		{
			summary.addValue(density[i]);
		}

		DensityManager dm = createDensityManager(results);

		// Compute this using the input density scores since the radius is the same.
		final double l = (useSquareApproximation) ? dm.ripleysLFunction(radius) : dm.ripleysLFunction(density, radius);

		String msg = String.format("Density %s : N=%d, %.0fpx : Radius=%s : L(r) - r = %s : E = %s, Obs = %s (%sx)",
				results.getName(), summary.getN(), area, rounded(radius), rounded(l - radius), rounded(expected),
				rounded(summary.getMean()), rounded(summary.getMean() / expected));
		if (filterLocalisations)
			msg += String.format(" : Filtered=%d (%s%%)", filtered, rounded(filtered * 100.0 / density.length));
		IJ.log(msg);

		return summary;
	}

	private static String rounded(double d)
	{
		return Utils.rounded(d, 3);
	}

	/**
	 * Draw an image of the density for each localisation. Optionally filter results below a threshold.
	 *
	 * @param results
	 *            the results
	 * @param density
	 *            the density
	 * @param scoreCalculator
	 *            the score calculator
	 * @return the number of localisations drawn
	 */
	private static int plotResults(MemoryPeakResults results, float[] density, ScoreCalculator scoreCalculator)
	{
		// Filter results using 5x higher than average density of the sample in a 150nm radius:
		// Annibale, et al (2011). Identification of clustering artifacts in photoactivated localization microscopy.
		// Nature Methods, 8, pp527-528
		MemoryPeakResults newResults = null;
		float densityThreshold = Float.NEGATIVE_INFINITY; // No filtering
		if (filterLocalisations)
		{
			densityThreshold = scoreCalculator.getThreshold();
			newResults = new MemoryPeakResults();
			newResults.copySettings(results);
			newResults.setName(results.getName() + " Density Filter");
		}

		// Draw an image - Use error so that a floating point value can be used on a single pixel
		IJImagePeakResults image = ImagePeakResultsFactory.createPeakResultsImage(ResultsImageType.DRAW_INTENSITY,
				false, false, results.getName() + " Density", results.getBounds(), results.getNmPerPixel(),
				results.getGain(), imageScale, 0,
				(cumulativeImage) ? ResultsImageMode.IMAGE_ADD : ResultsImageMode.IMAGE_MAX);
		image.setDisplayFlags(image.getDisplayFlags() | IJImagePeakResults.DISPLAY_NEGATIVES);
		image.setLutName("grays");
		image.begin();
		StandardResultProcedure sp = new StandardResultProcedure(newResults, DistanceUnit.PIXEL);
		sp.getXYR();
		for (int i = 0; i < density.length; i++)
		{
			if (density[i] < densityThreshold)
				continue;
			image.add(sp.x[i], sp.y[i], density[i]);
			if (newResults != null)
				newResults.add(sp.peakResults[i]);
		}
		image.end();

		// Add to memory
		if (newResults != null && newResults.size() > 0)
			MemoryPeakResults.addResults(newResults);

		return image.size();
	}

	private boolean showDialog()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		// Build a list of all images with a region ROI
		List<String> titles = new LinkedList<>();
		if (WindowManager.getWindowCount() > 0)
		{
			for (int imageID : WindowManager.getIDList())
			{
				ImagePlus imp = WindowManager.getImage(imageID);
				if (imp != null && imp.getRoi() != null && imp.getRoi().isArea())
					titles.add(imp.getTitle());
			}
		}

		gd.addMessage("Show an image using the localisation density");

		ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);

		gd.addNumericField("Radius", radius, 3);
		if (!titles.isEmpty())
			gd.addCheckbox((titles.size() == 1) ? "Use_ROI" : "Choose_ROI", chooseRoi);
		gd.addCheckbox("Adjust_for_border", adjustForBorder);
		gd.addSlider("Image_Scale", 1, 15, imageScale);
		gd.addCheckbox("Cumulative_image", cumulativeImage);

		gd.addCheckbox("Use_square_approx", useSquareApproximation);
		gd.addNumericField("Square_resolution", resolution, 0);
		gd.addChoice("Score", ScoreMethods, ScoreMethods[scoreMethodIndex]);

		gd.addMessage(
				"Filter localisations using the L-score / Relative density.\nFiltered results will be added to memory:");
		gd.addCheckbox("Filter_localisations", filterLocalisations);
		gd.addNumericField("Filter_threshold", filterThreshold, 2);

		gd.addCheckbox("Compute_Ripleys_L_plot", computeRipleysPlot);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		inputOption = ResultsManager.getInputSource(gd);

		radius = (float) gd.getNextNumber();
		if (!titles.isEmpty())
			chooseRoi = gd.getNextBoolean();
		adjustForBorder = gd.getNextBoolean();
		imageScale = (int) gd.getNextNumber();
		cumulativeImage = gd.getNextBoolean();

		useSquareApproximation = gd.getNextBoolean();
		resolution = (int) gd.getNextNumber();
		scoreMethodIndex = gd.getNextChoiceIndex();

		filterLocalisations = gd.getNextBoolean();
		filterThreshold = gd.getNextNumber();

		computeRipleysPlot = gd.getNextBoolean();

		// Check arguments
		try
		{
			Parameters.isAboveZero("Radius", radius);
			Parameters.isAboveZero("Image scale", imageScale);
			Parameters.isAboveZero("Resolution", resolution);
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		if (!titles.isEmpty() && chooseRoi)
		{
			if (titles.size() == 1)
			{
				roiImage = titles.get(0);
				Recorder.recordOption("Image", roiImage);
			}
			else
			{
				String[] items = titles.toArray(new String[titles.size()]);
				gd = new ExtendedGenericDialog(TITLE);
				gd.addMessage("Select the source image for the ROI");
				gd.addChoice("Image", items, roiImage);
				gd.showDialog();
				if (gd.wasCanceled())
					return false;
				roiImage = gd.getNextChoice();
			}
			ImagePlus imp = WindowManager.getImage(roiImage);

			roiBounds = imp.getRoi().getBounds();
			roiImageWidth = imp.getWidth();
			roiImageHeight = imp.getHeight();
		}
		else
		{
			roiBounds = null;
		}

		return true;
	}

	/**
	 * Compute the Ripley's L-function for user selected radii and show it on a plot.
	 *
	 * @param results
	 */
	private void computeRipleysPlot(MemoryPeakResults results)
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addMessage("Compute Ripley's L(r) - r plot");
		gd.addNumericField("Min_radius", minR, 2);
		gd.addNumericField("Max_radius", maxR, 2);
		gd.addNumericField("Increment", incrementR, 2);
		gd.addCheckbox("Confidence_intervals", confidenceIntervals);

		gd.showDialog();
		if (gd.wasCanceled())
			return;
		minR = gd.getNextNumber();
		maxR = gd.getNextNumber();
		incrementR = gd.getNextNumber();
		confidenceIntervals = gd.getNextBoolean();

		if (minR > maxR || incrementR < 0 || gd.invalidNumber())
		{
			IJ.error(TITLE, "Invalid radius parameters");
			return;
		}

		DensityManager dm = createDensityManager(results);
		double[][] values = calculateLScores(dm);

		// 99% confidence intervals
		final int iterations = (confidenceIntervals) ? 99 : 0;
		double[] upper = null;
		double[] lower = null;
		Rectangle bounds = results.getBounds();
		// Use a uniform distribution for the coordinates
		HaltonSequenceGenerator dist = new HaltonSequenceGenerator(2);
		dist.skipTo(new Well19937c(System.currentTimeMillis() + System.identityHashCode(this)).nextInt());
		for (int i = 0; i < iterations; i++)
		{
			IJ.showProgress(i, iterations);
			IJ.showStatus(String.format("L-score confidence interval %d / %d", i + 1, iterations));

			// Randomise coordinates
			float[] x = new float[results.size()];
			float[] y = new float[x.length];
			for (int j = x.length; j-- > 0;)
			{
				final double[] d = dist.nextVector();
				x[j] = (float) (d[0] * bounds.width);
				y[j] = (float) (d[1] * bounds.height);
			}
			double[][] values2 = calculateLScores(new DensityManager(x, y, bounds));
			if (upper == null || lower == null)
			{
				upper = values2[1];
				lower = upper.clone();
			}
			else
			{
				for (int m = upper.length; m-- > 0;)
				{
					if (upper[m] < values2[1][m])
						upper[m] = values2[1][m];
					if (lower[m] > values2[1][m])
						lower[m] = values2[1][m];
				}
			}
		}

		String title = results.getName() + " Ripley's (L(r) - r) / r";
		Plot2 plot = new Plot2(title, "Radius", "(L(r) - r) / r", values[0], values[1]);
		// Get the limits
		double yMin = min(0, values[1]);
		double yMax = max(0, values[1]);
		if (iterations > 0)
		{
			yMin = min(yMin, lower);
			yMax = max(yMax, upper);
		}
		plot.setLimits(0, values[0][values[0].length - 1], yMin, yMax);
		if (iterations > 0)
		{
			plot.setColor(Color.BLUE);
			plot.addPoints(values[0], upper, 1);
			plot.setColor(Color.RED);
			plot.addPoints(values[0], lower, 1);
			plot.setColor(Color.BLACK);
		}
		Utils.display(title, plot);
	}

	private static double min(double min, double[] data)
	{
		for (double d : data)
			if (min > d)
				min = d;
		return min;
	}

	private static double max(double max, double[] data)
	{
		for (double d : data)
			if (max < d)
				max = d;
		return max;
	}

	private static double[][] calculateLScores(DensityManager dm)
	{
		TDoubleArrayList x = new TDoubleArrayList();
		TDoubleArrayList y = new TDoubleArrayList();
		x.add(0.0);
		y.add(0.0);

		for (double r = minR; r < maxR; r += incrementR)
		{
			double l = dm.ripleysLFunction(r);
			x.add(r);
			double score = (r > 0) ? (l - r) / r : 0;
			y.add(score);
		}

		double[][] values = new double[2][];
		values[0] = x.toArray();
		values[1] = y.toArray();
		return values;
	}
}
