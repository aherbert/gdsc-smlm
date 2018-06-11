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
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import gdsc.core.ij.Utils;
import gdsc.core.utils.DoubleData;
import gdsc.core.utils.Maths;
import gdsc.core.utils.RollingArray;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.core.utils.Statistics;
import gdsc.core.utils.StoredData;
import gdsc.core.utils.StoredDataStatistics;
import gdsc.smlm.data.config.CalibrationHelper;
import gdsc.smlm.data.config.CreateDataSettingsHelper;
import gdsc.smlm.data.config.GUIProtos.CreateDataSettings;
import gdsc.smlm.data.config.PSFProtos.PSFType;
import gdsc.smlm.data.config.PSFHelper;

import gdsc.smlm.fitting.JumpDistanceAnalysis;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.model.DiffusionType;
import gdsc.smlm.model.ImageModel;
import gdsc.smlm.model.MoleculeModel;
import gdsc.smlm.model.SphericalDistribution;
import gdsc.smlm.results.ExtendedPeakResult;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.Plot2;
import ij.gui.PlotWindow;
import ij.plugin.LutLoader;
import ij.plugin.PlugIn;
import ij.plugin.WindowOrganiser;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import ij.text.TextWindow;

/**
 * Move a set of molecules and calculates the diffusion rate. Uses settings from the CreateData
 * plugin so that the diffusion should be equivalent.
 * 
 */
public class DiffusionRateTest implements PlugIn
{
	private static final String TITLE = "Diffusion Rate Test";
	private static TextWindow msdTable = null;

	// Used to allow other plugins to detect if a dataset is simulated
	static double lastSimulatedPrecision = 0;
	static String[] lastSimulatedDataset = new String[2];

	static boolean isSimulated(String name)
	{
		for (String name2 : lastSimulatedDataset)
			if (name.equals(name2))
				return true;
		return false;
	}

	/**
	 * Used to aggregate points into results
	 */
	public class Point
	{
		public int id;
		public double x, y;

		/**
		 * Create a cluster point
		 * 
		 * @param id
		 * @param x
		 * @param y
		 */
		public Point(int id, double x, double y)
		{
			this.id = id;
			this.x = x;
			this.y = y;
		}

		public Point(int id, double[] xyz)
		{
			this(id, xyz[0], xyz[1]);
		}

		public double distance2(Point p)
		{
			final double dx = x - p.x;
			final double dy = y - p.y;
			return dx * dx + dy * dy;
		}

		public double distance2(Point p, double error, RandomGenerator rand)
		{
			final double dx = (x + rand.nextGaussian() * error) - (p.x + rand.nextGaussian() * error);
			final double dy = (y + rand.nextGaussian() * error) - (p.y + rand.nextGaussian() * error);
			return dx * dx + dy * dy;
		}
	}

	private CreateDataSettings.Builder settings;
	private static boolean useConfinement = false;
	private static int confinementAttempts = 5;
	private static int fitN = 20;
	private static boolean showDiffusionExample = false;
	private static double magnification = 5;
	private static int aggregateSteps = 10;
	private static int msdAnalysisSteps = 0;
	private static double precision = 0;
	private int myAggregateSteps = 0;
	private int myMsdAnalysisSteps = 0;
	private boolean extraOptions = false;
	private double conversionFactor;
	private double myPrecision = 0;

	private int[] idList = new int[12];
	private int idCount = 0;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	@Override
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (IJ.controlKeyDown())
		{
			simpleTest();
			return;
		}

		extraOptions = Utils.isExtraOptions();

		if (!showDialog())
			return;

		lastSimulatedDataset[0] = lastSimulatedDataset[1] = "";
		lastSimulatedPrecision = 0;

		final int totalSteps = (int) Math.ceil(settings.getSeconds() * settings.getStepsPerSecond());

		conversionFactor = 1000000.0 / (settings.getPixelPitch() * settings.getPixelPitch());

		// Diffusion rate is um^2/sec. Convert to pixels per simulation frame.
		final double diffusionRateInPixelsPerSecond = settings.getDiffusionRate() * conversionFactor;
		final double diffusionRateInPixelsPerStep = diffusionRateInPixelsPerSecond / settings.getStepsPerSecond();
		final double precisionInPixels = myPrecision / settings.getPixelPitch();
		final boolean addError = myPrecision != 0;

		Utils.log(TITLE + " : D = %s um^2/sec, Precision = %s nm", Utils.rounded(settings.getDiffusionRate(), 4),
				Utils.rounded(myPrecision, 4));
		Utils.log("Mean-displacement per dimension = %s nm/sec",
				Utils.rounded(1e3 * ImageModel.getRandomMoveDistance(settings.getDiffusionRate()), 4));
		if (extraOptions)
			Utils.log("Step size = %s, precision = %s",
					Utils.rounded(ImageModel.getRandomMoveDistance(diffusionRateInPixelsPerStep)),
					Utils.rounded(precisionInPixels));

		// Convert diffusion co-efficient into the standard deviation for the random walk
		DiffusionType diffusionType = CreateDataSettingsHelper.getDiffusionType(settings.getDiffusionType());
		final double diffusionSigma = (diffusionType == DiffusionType.LINEAR_WALK)
				// Q. What should this be? At the moment just do 1D diffusion on a random vector
				? ImageModel.getRandomMoveDistance(diffusionRateInPixelsPerStep)
				: ImageModel.getRandomMoveDistance(diffusionRateInPixelsPerStep);
		Utils.log("Simulation step-size = %s nm", Utils.rounded(settings.getPixelPitch() * diffusionSigma, 4));

		// Move the molecules and get the diffusion rate
		IJ.showStatus("Simulating ...");
		final long start = System.nanoTime();
		final long seed = System.currentTimeMillis() + System.identityHashCode(this);
		RandomGenerator[] random = new RandomGenerator[3];
		RandomGenerator[] random2 = new RandomGenerator[3];
		for (int i = 0; i < 3; i++)
		{
			random[i] = new Well19937c(seed + i * 12436);
			random2[i] = new Well19937c(seed + i * 678678 + 3);
		}
		Statistics[] stats2D = new Statistics[totalSteps];
		Statistics[] stats3D = new Statistics[totalSteps];
		StoredDataStatistics jumpDistances2D = new StoredDataStatistics(totalSteps);
		StoredDataStatistics jumpDistances3D = new StoredDataStatistics(totalSteps);
		for (int j = 0; j < totalSteps; j++)
		{
			stats2D[j] = new Statistics();
			stats3D[j] = new Statistics();
		}
		SphericalDistribution dist = new SphericalDistribution(
				settings.getConfinementRadius() / settings.getPixelPitch());
		Statistics asymptote = new Statistics();

		// Save results to memory
		MemoryPeakResults results = new MemoryPeakResults(totalSteps);
		results.setCalibration(
				CalibrationHelper.create(settings.getPixelPitch(), 1, 1000.0 / settings.getStepsPerSecond()));
		results.setName(TITLE);
		results.setPSF(PSFHelper.create(PSFType.CUSTOM));
		int peak = 0;
		// Store raw coordinates
		ArrayList<Point> points = new ArrayList<Point>(totalSteps);

		StoredData totalJumpDistances1D = new StoredData(settings.getParticles());
		StoredData totalJumpDistances2D = new StoredData(settings.getParticles());
		StoredData totalJumpDistances3D = new StoredData(settings.getParticles());

		for (int i = 0; i < settings.getParticles(); i++)
		{
			if (i % 16 == 0)
			{
				IJ.showProgress(i, settings.getParticles());
				if (Utils.isInterrupted())
					return;
			}

			peak++; // Increment the frame so that tracing analysis can distinguish traces
			double[] origin = new double[3];
			final int id = i + 1;
			MoleculeModel m = new MoleculeModel(id, origin.clone());
			if (addError)
				origin = addError(origin, precisionInPixels, random);
			if (useConfinement)
			{
				// Note: When using confinement the average displacement should asymptote
				// at the average distance of a point from the centre of a ball. This is 3r/4.
				// See: http://answers.yahoo.com/question/index?qid=20090131162630AAMTUfM
				// The equivalent in 2D is 2r/3. However although we are plotting 2D distance
				// this is a projection of the 3D position onto the plane and so the particles
				// will not be evenly spread (there will be clustering at centre caused by the
				// poles)
				final double[] axis = (diffusionType == DiffusionType.LINEAR_WALK) ? nextVector() : null;
				for (int j = 0; j < totalSteps; j++)
				{
					double[] xyz = m.getCoordinates();
					double[] originalXyz = xyz.clone();
					for (int n = confinementAttempts; n-- > 0;)
					{
						if (diffusionType == DiffusionType.GRID_WALK)
							m.walk(diffusionSigma, random);
						else if (diffusionType == DiffusionType.LINEAR_WALK)
							m.slide(diffusionSigma, axis, random[0]);
						else
							m.move(diffusionSigma, random);

						if (!dist.isWithin(m.getCoordinates()))
						{
							// Reset position
							for (int k = 0; k < 3; k++)
								xyz[k] = originalXyz[k];
						}
						else
						{
							// The move was allowed
							break;
						}
					}

					points.add(new Point(id, xyz));

					if (addError)
						xyz = addError(xyz, precisionInPixels, random2);

					peak = record(xyz, id, peak, stats2D[j], stats3D[j], jumpDistances2D, jumpDistances3D, origin,
							results);
				}
				asymptote.add(distance(m.getCoordinates()));
			}
			else
			{
				if (diffusionType == DiffusionType.GRID_WALK)
				{
					for (int j = 0; j < totalSteps; j++)
					{
						m.walk(diffusionSigma, random);
						double[] xyz = m.getCoordinates();
						points.add(new Point(id, xyz));
						if (addError)
							xyz = addError(xyz, precisionInPixels, random2);
						peak = record(xyz, id, peak, stats2D[j], stats3D[j], jumpDistances2D, jumpDistances3D, origin,
								results);
					}
				}
				else if (diffusionType == DiffusionType.LINEAR_WALK)
				{
					final double[] axis = nextVector();
					for (int j = 0; j < totalSteps; j++)
					{
						m.slide(diffusionSigma, axis, random[0]);
						double[] xyz = m.getCoordinates();
						points.add(new Point(id, xyz));
						if (addError)
							xyz = addError(xyz, precisionInPixels, random2);
						peak = record(xyz, id, peak, stats2D[j], stats3D[j], jumpDistances2D, jumpDistances3D, origin,
								results);
					}
				}
				else
				{
					for (int j = 0; j < totalSteps; j++)
					{
						m.move(diffusionSigma, random);
						double[] xyz = m.getCoordinates();
						points.add(new Point(id, xyz));
						if (addError)
							xyz = addError(xyz, precisionInPixels, random2);
						peak = record(xyz, id, peak, stats2D[j], stats3D[j], jumpDistances2D, jumpDistances3D, origin,
								results);
					}
				}
			}

			// Debug: record all the particles so they can be analysed
			// System.out.printf("%f %f %f\n", m.getX(), m.getY(), m.getZ());
			final double[] xyz = m.getCoordinates();
			double d2 = 0;
			totalJumpDistances1D.add(d2 = xyz[0] * xyz[0]);
			totalJumpDistances2D.add(d2 += xyz[1] * xyz[1]);
			totalJumpDistances3D.add(d2 += xyz[2] * xyz[2]);
		}
		final double time = (System.nanoTime() - start) / 1000000.0;
		IJ.showProgress(1);

		MemoryPeakResults.addResults(results);
		lastSimulatedDataset[0] = results.getName();
		lastSimulatedPrecision = myPrecision;

		// Convert pixels^2/step to um^2/sec
		final double msd2D = (jumpDistances2D.getMean() / conversionFactor) /
				(results.getCalibrationReader().getExposureTime() / 1000);
		final double msd3D = (jumpDistances3D.getMean() / conversionFactor) /
				(results.getCalibrationReader().getExposureTime() / 1000);
		Utils.log(
				"Raw data D=%s um^2/s, Precision = %s nm, N=%d, step=%s s, mean2D=%s um^2, MSD 2D = %s um^2/s, mean3D=%s um^2, MSD 3D = %s um^2/s",
				Utils.rounded(settings.getDiffusionRate()), Utils.rounded(myPrecision), jumpDistances2D.getN(),
				Utils.rounded(results.getCalibrationReader().getExposureTime() / 1000),
				Utils.rounded(jumpDistances2D.getMean() / conversionFactor), Utils.rounded(msd2D),
				Utils.rounded(jumpDistances3D.getMean() / conversionFactor), Utils.rounded(msd3D));

		aggregateIntoFrames(points, addError, precisionInPixels, random2);

		IJ.showStatus("Analysing results ...");

		if (showDiffusionExample)
		{
			showExample(totalSteps, diffusionSigma, random);
		}

		// Plot a graph of mean squared distance
		double[] xValues = new double[stats2D.length];
		double[] yValues2D = new double[stats2D.length];
		double[] yValues3D = new double[stats3D.length];
		double[] upper2D = new double[stats2D.length];
		double[] lower2D = new double[stats2D.length];
		double[] upper3D = new double[stats3D.length];
		double[] lower3D = new double[stats3D.length];

		SimpleRegression r2D = new SimpleRegression(false);
		SimpleRegression r3D = new SimpleRegression(false);

		final int firstN = (useConfinement) ? fitN : totalSteps;
		for (int j = 0; j < totalSteps; j++)
		{
			// Convert steps to seconds
			xValues[j] = (j + 1) / settings.getStepsPerSecond();

			// Convert values in pixels^2 to um^2
			final double mean2D = stats2D[j].getMean() / conversionFactor;
			final double mean3D = stats3D[j].getMean() / conversionFactor;
			final double sd2D = stats2D[j].getStandardDeviation() / conversionFactor;
			final double sd3D = stats3D[j].getStandardDeviation() / conversionFactor;
			yValues2D[j] = mean2D;
			yValues3D[j] = mean3D;
			upper2D[j] = mean2D + sd2D;
			lower2D[j] = mean2D - sd2D;
			upper3D[j] = mean3D + sd3D;
			lower3D[j] = mean3D - sd3D;

			if (j < firstN)
			{
				r2D.addData(xValues[j], yValues2D[j]);
				r3D.addData(xValues[j], yValues3D[j]);
			}
		}

		// TODO - Fit using the equation for 2D confined diffusion:
		// MSD = 4s^2 + R^2 (1 - 0.99e^(-1.84^2 Dt / R^2)
		// s = localisation precision
		// R = confinement radius
		// D = 2D diffusion coefficient
		// t = time

		final PolynomialFunction fitted2D, fitted3D;
		if (r2D.getN() > 0)
		{
			// Do linear regression to get diffusion rate

			final double[] best2D = new double[] { r2D.getIntercept(), r2D.getSlope() };
			fitted2D = new PolynomialFunction(best2D);

			final double[] best3D = new double[] { r3D.getIntercept(), r3D.getSlope() };
			fitted3D = new PolynomialFunction(best3D);

			// For 2D diffusion: d^2 = 4D
			// where: d^2 = mean-square displacement

			double D = best2D[1] / 4.0;
			String msg = "2D Diffusion rate = " + Utils.rounded(D, 4) + " um^2 / sec (" + Utils.timeToString(time) +
					")";
			IJ.showStatus(msg);
			Utils.log(msg);

			D = best3D[1] / 6.0;
			Utils.log("3D Diffusion rate = " + Utils.rounded(D, 4) + " um^2 / sec (" + Utils.timeToString(time) + ")");
		}
		else
		{
			fitted2D = fitted3D = null;
		}

		// Create plots
		plotMSD(totalSteps, xValues, yValues2D, lower2D, upper2D, fitted2D, 2);
		plotMSD(totalSteps, xValues, yValues3D, lower3D, upper3D, fitted3D, 3);

		plotJumpDistances(TITLE, jumpDistances2D, 2, 1);
		plotJumpDistances(TITLE, jumpDistances3D, 3, 1);

		// Show the total jump length for debugging
		//plotJumpDistances(TITLE + " total", totalJumpDistances1D, 1, totalSteps);
		//plotJumpDistances(TITLE + " total", totalJumpDistances2D, 2, totalSteps);
		//plotJumpDistances(TITLE + " total", totalJumpDistances3D, 3, totalSteps);

		if (idCount > 0)
			new WindowOrganiser().tileWindows(idList);

		if (useConfinement)
			Utils.log("3D asymptote distance = %s nm (expected %.2f)",
					Utils.rounded(asymptote.getMean() * settings.getPixelPitch(), 4),
					3 * settings.getConfinementRadius() / 4);
	}

	/**
	 * Plot the MSD.
	 *
	 * @param totalSteps
	 *            the total steps
	 * @param xValues
	 *            the x values (the time)
	 * @param yValues
	 *            the y values (MSD)
	 * @param lower
	 *            the lower bounds (mean-SD)
	 * @param upper
	 *            the upper bounds (mean+SD)
	 * @param fitted
	 *            the fitted line
	 * @param dimensions
	 *            the number of dimensions for the jumps
	 */
	private void plotMSD(int totalSteps, double[] xValues, double[] yValues, double[] lower, double[] upper,
			PolynomialFunction fitted, int dimensions)
	{
		// TODO Auto-generated method stub
		String title = TITLE + " " + dimensions + "D";
		Plot2 plot = new Plot2(title, "Time (seconds)", "Mean-squared Distance (um^2)", xValues, yValues);
		double[] limits = Maths.limits(upper);
		limits = Maths.limits(limits, lower);
		plot.setLimits(0, totalSteps / settings.getStepsPerSecond(), limits[0], limits[1]);
		plot.setColor(Color.blue);
		plot.addPoints(xValues, lower, Plot.LINE);
		plot.addPoints(xValues, upper, Plot.LINE);
		if (fitted != null)
		{
			plot.setColor(Color.red);
			plot.addPoints(new double[] { xValues[0], xValues[xValues.length - 1] },
					new double[] { fitted.value(xValues[0]), fitted.value(xValues[xValues.length - 1]) }, Plot.LINE);
		}
		plot.setColor(Color.black);

		PlotWindow pw1 = Utils.display(title, plot);
		if (Utils.isNewWindow())
			idList[idCount++] = pw1.getImagePlus().getID();
	}

	/**
	 * Plot a cumulative histogram and standard histogram of the jump distances.
	 *
	 * @param title
	 *            the title
	 * @param jumpDistances
	 *            the jump distances
	 * @param dimensions
	 *            the number of dimensions for the jumps
	 * @param steps
	 *            the steps
	 */
	private void plotJumpDistances(String title, DoubleData jumpDistances, int dimensions, int steps)
	{
		// Cumulative histogram
		// --------------------
		final double factor = conversionFactor;
		double[] values = jumpDistances.values();
		for (int i = 0; i < values.length; i++)
			values[i] /= factor;
		String title2 = title + " Cumulative Jump Distance " + dimensions + "D";
		double[][] jdHistogram = JumpDistanceAnalysis.cumulativeHistogram(values);
		DiffusionType diffusionType = CreateDataSettingsHelper.getDiffusionType(settings.getDiffusionType());
		if (diffusionType == DiffusionType.GRID_WALK)
		{
			// In this case with a large simulation size the jumps are all 
			// the same distance so the histogram is a single step. Check the plot
			// range will be handled by ImageJ otherwise pad it out a bit.
			double[] x = jdHistogram[0];
			double[] y = jdHistogram[1];
			if (x[x.length - 1] - x[0] < 0.01)
			{
				double[] x2 = new double[x.length + 3];
				double[] y2 = new double[y.length + 3];
				System.arraycopy(x, 0, x2, 2, x.length);
				System.arraycopy(y, 0, y2, 2, y.length);
				x2[0] = x[0] - 0.1;
				x2[1] = x[0];
				x2[x2.length - 1] = x[x.length - 1] + 0.1;
				y2[0] = 0;
				y2[1] = 0;
				y2[y2.length - 1] = 1;
				jdHistogram[0] = x2;
				jdHistogram[1] = y2;

				// Add some artificial points to allow the plot to be drawn
				values = Arrays.copyOf(values, values.length + 2);
				values[values.length - 2] = x[0] - 0.1;
				values[values.length - 1] = x[x.length - 1] + 0.1;
			}
		}
		Plot2 jdPlot = new Plot2(title2, "Distance (um^2)", "Cumulative Probability", jdHistogram[0], jdHistogram[1]);
		PlotWindow pw2 = Utils.display(title2, jdPlot);
		if (Utils.isNewWindow())
			idList[idCount++] = pw2.getImagePlus().getID();

		// This is the Chi-squared distribution: The sum of the squares of k independent
		// standard normal random variables with k = dimensions. It is a special case of
		// the gamma distribution. If the normals have non-unit variance the distribution 
		// is scaled.
		// Chi       ~ Gamma(k/2, 2)      // using the scale parameterisation of the gamma
		// s^2 * Chi ~ Gamma(k/2, 2*s^2)
		// So if s^2 = 2D:
		// 2D * Chi  ~ Gamma(k/2, 4D)
		double estimatedD = steps * settings.getDiffusionRate() / settings.getStepsPerSecond();
		if (myPrecision > 0)
			estimatedD += myPrecision * myPrecision / 1e6;
		double max = Maths.max(values);
		double[] x = SimpleArrayUtils.newArray(1000, 0, max / 1000);
		double k = dimensions / 2.0;
		double mean = 4 * estimatedD;

		GammaDistribution dist = new GammaDistribution(k, mean);

		double[] y = new double[x.length];
		for (int i = 0; i < x.length; i++)
			y[i] = dist.cumulativeProbability(x[i]);

		jdPlot.setColor(Color.red);
		jdPlot.addPoints(x, y, Plot.LINE);
		Utils.display(title2, jdPlot);

		// Histogram
		// ---------
		title2 = title + " Jump " + dimensions + "D";
		StoredDataStatistics jumpDistances2 = new StoredDataStatistics(values);
		int plotId = Utils.showHistogram(title2, jumpDistances2, "Distance (um^2)", 0, 0,
				Math.max(20, values.length / 1000));
		if (Utils.isNewWindow())
			idList[idCount++] = plotId;

		// Recompute the expected function
		for (int i = 0; i < x.length; i++)
			y[i] = dist.density(x[i]);

		// Scale to have the same area
		if (Utils.xValues.length > 1)
		{
			final double area1 = jumpDistances2.getN() * (Utils.xValues[1] - Utils.xValues[0]);
			final double area2 = dist.cumulativeProbability(x[x.length - 1]);
			final double scale = area1 / area2;
			for (int i = 0; i < y.length; i++)
				y[i] *= scale;
		}
		jdPlot = Utils.plot;
		jdPlot.setColor(Color.red);
		jdPlot.addPoints(x, y, Plot.LINE);
		Utils.display(WindowManager.getImage(plotId).getTitle(), jdPlot);
	}

	/**
	 * Add a random Gaussian XY shift using the specified precision
	 * 
	 * @param xyz
	 * @param precision
	 * @param random
	 * @return The new xyz
	 */
	private double[] addError(double[] xyz, double precision, RandomGenerator[] random)
	{
		final double[] xy = xyz.clone();
		for (int i = 0; i < 2; i++)
		{
			final double shift = random[i].nextGaussian() * precision;
			xy[i] += shift;
		}
		return xy;
	}

	private int record(double[] xyz, int id, int peak, Statistics stats2D, Statistics stats3D,
			StoredDataStatistics jumpDistances2D, StoredDataStatistics jumpDistances3D, double[] origin,
			MemoryPeakResults results)
	{
		final double dx = xyz[0] - origin[0];
		final double dy = xyz[1] - origin[1];
		final double dz = xyz[2] - origin[2];
		final double jump2D = dx * dx + dy * dy;
		jumpDistances2D.add(jump2D);
		jumpDistances3D.add(jump2D + dz * dz);
		for (int i = 0; i < 3; i++)
			origin[i] = xyz[i];

		final double d2 = squared2D(xyz);
		stats2D.add(d2);
		stats3D.add(d2 + xyz[2] * xyz[2]);

		final float[] params = PeakResult.createParams(0, 10f, (float) xyz[0], (float) xyz[1], (float) xyz[2]);
		final float noise = 0.1f;
		results.add(new ExtendedPeakResult(peak, (int) params[PeakResult.X], (int) params[PeakResult.Y], 10, 0, noise,
				0, params, null, peak, id));
		return ++peak;
	}

	/**
	 * Get the squared distance from the origin in 2D (using XY coordinates)
	 * 
	 * @param coordinates
	 * @return
	 */
	private double squared2D(double[] coordinates)
	{
		return coordinates[0] * coordinates[0] + coordinates[1] * coordinates[1];
	}

	/**
	 * Get the distance from the origin in 3D
	 * 
	 * @param coordinates
	 * @return
	 */
	private double distance(double[] coordinates)
	{
		return Math.sqrt(
				coordinates[0] * coordinates[0] + coordinates[1] * coordinates[1] + coordinates[2] * coordinates[2]);
	}

	private boolean showDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);

		settings = SettingsManager.readCreateDataSettings(0).toBuilder();

		if (settings.getStepsPerSecond() < 1)
			settings.setStepsPerSecond(1);

		gd.addNumericField("Pixel_pitch (nm)", settings.getPixelPitch(), 2);
		gd.addNumericField("Seconds", settings.getSeconds(), 1);
		gd.addSlider("Steps_per_second", 1, 15, settings.getStepsPerSecond());
		if (extraOptions)
		{
			gd.addSlider("Aggregate_steps", 2, 20, aggregateSteps);
			gd.addNumericField("MSD_analysis_steps", msdAnalysisSteps, 0);
		}
		gd.addNumericField("Particles", settings.getParticles(), 0);
		gd.addNumericField("Diffusion_rate (um^2/sec)", settings.getDiffusionRate(), 2);
		if (extraOptions)
			gd.addNumericField("Precision (nm)", precision, 2);
		String[] diffusionTypes = SettingsManager.getNames((Object[]) DiffusionType.values());
		gd.addChoice("Diffusion_type", diffusionTypes, diffusionTypes[settings.getDiffusionType()]);
		gd.addCheckbox("Use_confinement", useConfinement);
		gd.addSlider("Confinement_attempts", 1, 20, confinementAttempts);
		gd.addSlider("Confinement_radius (nm)", 0, 3000, settings.getConfinementRadius());
		gd.addSlider("Fit_N", 5, 20, fitN);
		gd.addCheckbox("Show_example", showDiffusionExample);
		gd.addSlider("Magnification", 1, 10, magnification);

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		settings.setPixelPitch(Math.abs(gd.getNextNumber()));
		settings.setSeconds(Math.abs(gd.getNextNumber()));
		settings.setStepsPerSecond(Math.abs(gd.getNextNumber()));
		if (extraOptions)
		{
			myAggregateSteps = aggregateSteps = Math.abs((int) gd.getNextNumber());
			myMsdAnalysisSteps = msdAnalysisSteps = Math.abs((int) gd.getNextNumber());
		}
		settings.setParticles(Math.abs((int) gd.getNextNumber()));
		settings.setDiffusionRate(Math.abs(gd.getNextNumber()));
		if (extraOptions)
			myPrecision = precision = Math.abs(gd.getNextNumber());
		settings.setDiffusionType(gd.getNextChoiceIndex());
		useConfinement = gd.getNextBoolean();
		confinementAttempts = Math.abs((int) gd.getNextNumber());
		settings.setConfinementRadius(Math.abs(gd.getNextNumber()));
		fitN = Math.abs((int) gd.getNextNumber());
		showDiffusionExample = gd.getNextBoolean();
		magnification = gd.getNextNumber();

		// Save before validation so that the current values are preserved.
		SettingsManager.writeSettings(settings.build());

		// Check arguments
		try
		{
			Parameters.isAboveZero("Pixel Pitch", settings.getPixelPitch());
			Parameters.isAboveZero("Seconds", settings.getSeconds());
			Parameters.isAboveZero("Steps per second", settings.getStepsPerSecond());
			Parameters.isAboveZero("Particles", settings.getParticles());
			Parameters.isPositive("Diffusion rate", settings.getDiffusionRate());
			Parameters.isAboveZero("Magnification", magnification);
			Parameters.isAboveZero("Confinement attempts", confinementAttempts);
			Parameters.isAboveZero("Fit N", fitN);
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		if (settings.getDiffusionRate() == 0)
			IJ.error(TITLE, "Warning : Diffusion rate is zero");

		if (gd.invalidNumber())
			return false;

		SettingsManager.writeSettings(settings.build());

		return true;
	}

	private RandomGenerator rg = null;

	private double[] nextVector()
	{
		if (rg == null)
			rg = new Well19937c();

		// Allow dimensions to be changed for testing
		double l = 0;
		final double[] v = new double[3];
		final int size = 3;
		final int dim = 3; // Normalise over a different size 
		while (l == 0)
		{
			for (int i = 0; i < size; i++)
			{
				v[i] = rg.nextGaussian();
			}
			for (int i = 0; i < dim; i++)
			{
				l += v[i] * v[i];
			}
		}
		l = Math.sqrt(l);
		for (int i = 0; i < size; i++)
			v[i] /= l;
		return v;
	}

	private void showExample(int totalSteps, double diffusionSigma, RandomGenerator[] random)
	{
		MoleculeModel m = new MoleculeModel(0, new double[3]);
		float[] xValues = new float[totalSteps];
		float[] x = new float[totalSteps];
		float[] y = new float[totalSteps];
		DiffusionType diffusionType = CreateDataSettingsHelper.getDiffusionType(settings.getDiffusionType());
		final double[] axis = (diffusionType == DiffusionType.LINEAR_WALK) ? nextVector() : null;
		for (int j = 0; j < totalSteps; j++)
		{
			if (diffusionType == DiffusionType.GRID_WALK)
				m.walk(diffusionSigma, random);
			else if (diffusionType == DiffusionType.LINEAR_WALK)
				m.slide(diffusionSigma, axis, random[0]);
			else
				m.move(diffusionSigma, random);
			x[j] = (float) (m.getX());
			y[j] = (float) (m.getY());
			xValues[j] = (float) ((j + 1) / settings.getStepsPerSecond());
		}

		// Plot x and y coords on a timeline
		String title = TITLE + " example coordinates";
		Plot2 plot = new Plot2(title, "Time (seconds)", "Distance (um)");
		float[] xUm = convertToUm(x);
		float[] yUm = convertToUm(y);
		float[] limits = Maths.limits(xUm);
		limits = Maths.limits(limits, yUm);
		plot.setLimits(0, totalSteps / settings.getStepsPerSecond(), limits[0], limits[1]);
		plot.setColor(Color.red);
		plot.addPoints(xValues, xUm, Plot.LINE);
		plot.setColor(Color.blue);
		plot.addPoints(xValues, yUm, Plot.LINE);

		Utils.display(title, plot);

		// Scale up and draw 2D position
		for (int j = 0; j < totalSteps; j++)
		{
			x[j] *= magnification;
			y[j] *= magnification;
		}
		float[] limitsx = getLimits(x);
		float[] limitsy = getLimits(y);

		int width = (int) (limitsx[1] - limitsx[0]);
		int height = (int) (limitsy[1] - limitsy[0]);

		// Ensure we draw something, even it is a simple dot at the centre for no diffusion
		if (width == 0)
		{
			width = (int) (32 * magnification);
			limitsx[0] = -width / 2;
		}
		if (height == 0)
		{
			height = (int) (32 * magnification);
			limitsy[0] = -height / 2;
		}

		ImageProcessor ip = new ByteProcessor(width, height);

		// Adjust x and y using the minimum to centre
		x[0] -= limitsx[0];
		y[0] -= limitsy[0];

		for (int j = 1; j < totalSteps; j++)
		{
			// Adjust x and y using the minimum to centre
			x[j] -= limitsx[0];
			y[j] -= limitsy[0];

			// Draw a line
			ip.setColor(32 + (223 * j) / (totalSteps - 1));
			ip.drawLine(round(x[j - 1]), round(y[j - 1]), round(x[j]), round(y[j]));
		}
		// Draw the final position
		ip.putPixel(round(x[totalSteps - 1]), round(y[totalSteps - 1]), 255);

		ImagePlus imp = Utils.display(TITLE + " example", ip);

		// Apply the fire lookup table
		WindowManager.setTempCurrentImage(imp);
		LutLoader lut = new LutLoader();
		lut.run("fire");
		WindowManager.setTempCurrentImage(null);

	}

	private float[] convertToUm(float[] x)
	{
		final float factor = (float) (settings.getPixelPitch() / 1000.0);
		float[] newX = new float[x.length];
		for (int i = 0; i < x.length; i++)
			newX[i] = x[i] * factor;
		return newX;
	}

	private int round(float f)
	{
		return Math.round(f);
	}

	private float[] getLimits(float[] x)
	{
		float[] limits = Maths.limits(x);
		limits[0] = (float) Math.floor(limits[0]);
		limits[1] = (float) Math.ceil(limits[1]);
		return limits;
	}

	private void aggregateIntoFrames(ArrayList<Point> points, boolean addError, double precisionInPixels,
			RandomGenerator[] random)
	{
		if (myAggregateSteps < 1)
			return;

		MemoryPeakResults results = new MemoryPeakResults(points.size() / myAggregateSteps);
		results.setCalibration(CalibrationHelper.create(settings.getPixelPitch(), 1,
				myAggregateSteps * 1000.0 / settings.getStepsPerSecond()));
		results.setName(TITLE + " Aggregated");
		results.setPSF(PSFHelper.create(PSFType.CUSTOM));
		MemoryPeakResults.addResults(results);
		lastSimulatedDataset[1] = results.getName();
		int id = 0;
		int peak = 1;
		int n = 0;
		double cx = 0, cy = 0;
		// Get the mean square distance
		double sum = 0;
		int count = 0;
		PeakResult last = null;
		for (Point result : points)
		{
			final boolean newId = result.id != id;
			if (n >= myAggregateSteps || newId)
			{
				if (n != 0)
				{
					double[] xyz = new double[] { cx / n, cy / n };
					if (addError)
						xyz = addError(xyz, precisionInPixels, random);
					final float[] params = PeakResult.createParams(0, n, (float) xyz[0], (float) xyz[1], 0);
					final float noise = 0.1f;
					PeakResult r = new ExtendedPeakResult(peak, (int) params[PeakResult.X], (int) params[PeakResult.Y],
							n, 0, noise, 0, params, null, peak, id);
					results.add(r);
					if (last != null)
					{
						sum += last.distance2(r);
						count++;
					}
					last = r;
					n = 0;
					cx = cy = 0;
					peak++;
				}
				if (newId)
				{
					peak++; // Increment the frame so that tracing analysis can distinguish traces
					last = null;
					id = result.id;
				}
			}
			n++;
			cx += result.x;
			cy += result.y;
		}

		// Final peak
		if (n != 0)
		{
			double[] xyz = new double[] { cx / n, cy / n };
			if (addError)
				xyz = addError(xyz, precisionInPixels, random);
			final float[] params = PeakResult.createParams(0, n, (float) xyz[0], (float) xyz[1], 0);
			final float noise = 0.1f;
			PeakResult r = new ExtendedPeakResult(peak, (int) params[PeakResult.X], (int) params[PeakResult.Y], n, 0,
					noise, 0, params, null, peak, id);
			results.add(r);
			if (last != null)
			{
				sum += last.distance2(r);
				count++;
			}
		}

		// MSD in pixels^2 / frame
		double msd = sum / count;
		// Convert to um^2/second
		Utils.log("Aggregated data D=%s um^2/s, Precision=%s nm, N=%d, step=%s s, mean=%s um^2, MSD = %s um^2/s",
				Utils.rounded(settings.getDiffusionRate()), Utils.rounded(myPrecision), count,
				Utils.rounded(results.getCalibrationReader().getExposureTime() / 1000),
				Utils.rounded(msd / conversionFactor),
				Utils.rounded((msd / conversionFactor) / (results.getCalibrationReader().getExposureTime() / 1000)));

		msdAnalysis(points);
	}

	/**
	 * Tabulate the observed MSD for different jump distances
	 * 
	 * @param points
	 */
	private void msdAnalysis(ArrayList<Point> points)
	{
		if (myMsdAnalysisSteps == 0)
			return;

		IJ.showStatus("MSD analysis ...");
		IJ.showProgress(1, myMsdAnalysisSteps);

		// This will only be fast if the list is an array
		Point[] list = points.toArray(new Point[points.size()]);

		// Compute the base MSD
		Point origin = new Point(0, 0, 0);
		double sum = origin.distance2(list[0]);
		int count = 1;
		for (int i = 1; i < list.length; i++)
		{
			Point last = list[i - 1];
			Point current = list[i];

			if (last.id == current.id)
			{
				sum += last.distance2(current);
			}
			else
			{
				sum += origin.distance2(current);
			}
			count++;
		}
		createMsdTable((sum / count) * settings.getStepsPerSecond() / conversionFactor);

		// Create a new set of points that have coordinates that 
		// are the rolling average over the number of aggregate steps
		RollingArray x = new RollingArray(aggregateSteps);
		RollingArray y = new RollingArray(aggregateSteps);

		int id = 0;
		int length = 0;
		for (Point p : points)
		{
			if (p.id != id)
			{
				x.clear();
				y.clear();
			}
			id = p.id;
			x.add(p.x);
			y.add(p.y);
			// Only create a point if the full aggregation size is reached
			if (x.isFull())
			{
				list[length++] = new Point(id, x.getAverage(), y.getAverage());
			}
		}

		// Q - is this useful?
		final double p = myPrecision / settings.getPixelPitch();
		final long seed = System.currentTimeMillis() + System.identityHashCode(this);
		RandomGenerator rand = new Well19937c(seed);

		final int totalSteps = (int) Math.ceil(settings.getSeconds() * settings.getStepsPerSecond() - aggregateSteps);
		final int limit = Math.min(totalSteps, myMsdAnalysisSteps);
		final int interval = Utils.getProgressInterval(limit);
		final ArrayList<String> results = new ArrayList<String>(totalSteps);
		for (int step = 1; step <= myMsdAnalysisSteps; step++)
		{
			if (step % interval == 0)
				IJ.showProgress(step, limit);

			sum = 0;
			count = 0;
			for (int i = step; i < length; i++)
			{
				Point last = list[i - step];
				Point current = list[i];

				if (last.id == current.id)
				{
					if (p == 0)
					{
						sum += last.distance2(current);
						count++;
					}
					else
					{
						// This can be varied but the effect on the output with only 1 loop 
						// is the same if enough samples are present
						for (int ii = 1; ii-- > 0;)
						{
							sum += last.distance2(current, p, rand);
							count++;
						}
					}
				}
			}
			if (count == 0)
				break;
			results.add(addResult(step, sum, count));

			// Flush to auto-space the columns
			if (step == 9)
			{
				msdTable.getTextPanel().append(results);
				results.clear();
			}
		}
		msdTable.getTextPanel().append(results);

		IJ.showProgress(1);
	}

	private void createMsdTable(double baseMsd)
	{
		String header = createHeader(baseMsd);
		if (msdTable == null || !msdTable.isVisible())
		{
			msdTable = new TextWindow("MSD Analysis", header, "", 800, 300);
			msdTable.setVisible(true);
		}
		else
		{
			//msdTable.getTextPanel().clear();
		}
	}

	private String prefix;
	private double exposureTime;

	private String createHeader(double baseMsd)
	{
		final double apparentD = baseMsd / 4;
		StringBuilder sb = new StringBuilder();
		sb.append(settings.getDiffusionRate()).append('\t');
		sb.append(myPrecision).append('\t');
		sb.append(Utils.rounded(apparentD)).append('\t');
		sb.append(Utils.rounded(1.0 / settings.getStepsPerSecond())).append('\t');
		sb.append(myAggregateSteps).append('\t');
		// Exposure time is the aggregated frame time 
		exposureTime = myAggregateSteps / settings.getStepsPerSecond();
		sb.append(Utils.rounded(exposureTime)).append('\t');
		prefix = sb.toString();
		return "D (um^2/s)\tPrecision (nm)\tDsim (um^2/s)\tStep (s)\tResolution\tFrame (s)\tt (s)\tn\tN\tMSD (um^2)\tD (um^2/s)";
	}

	private String addResult(int step, double sum, int count)
	{
		StringBuilder sb = new StringBuilder();
		// Exposure time is the aggregated frame time 
		final double msd = (sum / count) / conversionFactor;
		// Jump distance separation is the number of steps
		final double t = step / settings.getStepsPerSecond();
		sb.append(prefix);
		sb.append(Utils.rounded(t)).append('\t');
		sb.append(Utils.rounded(t / exposureTime)).append('\t');
		sb.append(count).append('\t');
		// Not rounded to preserve precision 
		sb.append(msd).append('\t');
		sb.append(msd / (4 * t));
		return sb.toString();
	}

	private static double simpleD = 0.5;
	private static int simpleSteps = 1;
	private static int simpleParticles = 10000;
	private static boolean linearDiffusion = false;
	private static String simpleDir = null;

	/**
	 * Perform a simple diffusion test. This can be used to understand the distributions that are generated during
	 * 3D diffusion.
	 */
	private void simpleTest()
	{
		if (!showSimpleDialog())
			return;

		StoredDataStatistics[] stats2 = new StoredDataStatistics[3];
		StoredDataStatistics[] stats = new StoredDataStatistics[3];
		RandomGenerator[] random = new RandomGenerator[3];
		final long seed = System.currentTimeMillis() + System.identityHashCode(this);
		for (int i = 0; i < 3; i++)
		{
			stats2[i] = new StoredDataStatistics(simpleParticles);
			stats[i] = new StoredDataStatistics(simpleParticles);
			random[i] = new Well19937c(seed + i);
		}

		final double scale = Math.sqrt(2 * simpleD);
		final int report = Math.max(1, simpleParticles / 200);
		for (int particle = 0; particle < simpleParticles; particle++)
		{
			if (particle % report == 0)
				IJ.showProgress(particle, simpleParticles);
			double[] xyz = new double[3];
			if (linearDiffusion)
			{
				double[] dir = nextVector();
				for (int step = 0; step < simpleSteps; step++)
				{
					final double d = ((random[1].nextDouble() > 0.5) ? -1 : 1) * random[0].nextGaussian();
					for (int i = 0; i < 3; i++)
					{
						xyz[i] += dir[i] * d;
					}
				}
			}
			else
			{
				for (int step = 0; step < simpleSteps; step++)
				{
					for (int i = 0; i < 3; i++)
					{
						xyz[i] += random[i].nextGaussian();
					}
				}
			}
			for (int i = 0; i < 3; i++)
				xyz[i] *= scale;
			double msd = 0;
			for (int i = 0; i < 3; i++)
			{
				msd += xyz[i] * xyz[i];
				stats2[i].add(msd);
				// Store the actual distances
				stats[i].add(xyz[i]);
			}
		}
		IJ.showProgress(1);

		for (int i = 0; i < 3; i++)
		{
			plotJumpDistances(TITLE, stats2[i], i + 1);
			// Save stats to file for fitting
			save(stats2[i], i + 1, "msd");
			save(stats[i], i + 1, "d");
		}

		if (idCount > 0)
			new WindowOrganiser().tileWindows(idList);
	}

	private boolean showSimpleDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);

		gd.addNumericField("D", simpleD, 2);
		gd.addNumericField("Steps", simpleSteps, 0);
		gd.addNumericField("Particles", simpleParticles, 0);
		gd.addCheckbox("Linear_diffusion", linearDiffusion);
		gd.addStringField("Directory", simpleDir, 30);

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		simpleD = gd.getNextNumber();
		simpleSteps = (int) gd.getNextNumber();
		simpleParticles = (int) gd.getNextNumber();
		linearDiffusion = gd.getNextBoolean();
		simpleDir = gd.getNextString();
		if (!new File(simpleDir).exists())
			simpleDir = null;

		return true;
	}

	/**
	 * Plot a cumulative histogram and standard histogram of the jump distances.
	 *
	 * @param title
	 *            the title
	 * @param jumpDistances
	 *            the jump distances
	 * @param dimensions
	 *            the number of dimensions for the jumps
	 * @param steps
	 *            the steps
	 */
	private void plotJumpDistances(String title, DoubleData jumpDistances, int dimensions)
	{
		// Cumulative histogram
		// --------------------
		double[] values = jumpDistances.values();
		String title2 = title + " Cumulative Jump Distance " + dimensions + "D";
		double[][] jdHistogram = JumpDistanceAnalysis.cumulativeHistogram(values);
		Plot2 jdPlot = new Plot2(title2, "Distance (um^2)", "Cumulative Probability", jdHistogram[0], jdHistogram[1]);
		PlotWindow pw2 = Utils.display(title2, jdPlot);
		if (Utils.isNewWindow())
			idList[idCount++] = pw2.getImagePlus().getID();

		// Plot the expected function
		// This is the Chi-squared distribution: The sum of the squares of k independent
		// standard normal random variables with k = dimensions. It is a special case of
		// the gamma distribution. If the normals have non-unit variance the distribution 
		// is scaled.
		// Chi       ~ Gamma(k/2, 2)      // using the scale parameterisation of the gamma
		// s^2 * Chi ~ Gamma(k/2, 2*s^2)
		// So if s^2 = 2D:
		// 2D * Chi  ~ Gamma(k/2, 4D)
		double estimatedD = simpleD * simpleSteps;
		double max = Maths.max(values);
		double[] x = SimpleArrayUtils.newArray(1000, 0, max / 1000);
		double k = dimensions / 2.0;
		double mean = 4 * estimatedD;

		GammaDistribution dist = new GammaDistribution(k, mean);

		double[] y = new double[x.length];
		for (int i = 0; i < x.length; i++)
			y[i] = dist.cumulativeProbability(x[i]);

		jdPlot.setColor(Color.red);
		jdPlot.addPoints(x, y, Plot.LINE);
		Utils.display(title2, jdPlot);

		// Histogram
		// ---------
		title2 = title + " Jump " + dimensions + "D";
		int plotId = Utils.showHistogram(title2, jumpDistances, "Distance (um^2)", 0, 0,
				Math.max(20, values.length / 1000));
		if (Utils.isNewWindow())
			idList[idCount++] = plotId;

		// Recompute the expected function
		for (int i = 0; i < x.length; i++)
			y[i] = dist.density(x[i]);

		// Scale to have the same area
		if (Utils.xValues.length > 1)
		{
			final double area1 = jumpDistances.size() * (Utils.xValues[1] - Utils.xValues[0]);
			final double area2 = dist.cumulativeProbability(x[x.length - 1]);
			final double scaleFactor = area1 / area2;
			for (int i = 0; i < y.length; i++)
				y[i] *= scaleFactor;
		}
		jdPlot = Utils.plot;
		jdPlot.setColor(Color.red);
		jdPlot.addPoints(x, y, Plot.LINE);
		Utils.display(WindowManager.getImage(plotId).getTitle(), jdPlot);
	}

	private void save(StoredDataStatistics storedDataStatistics, int dimensions, String prefix)
	{
		if (simpleDir == null)
			return;
		File file = new File(simpleDir, prefix + dimensions + "d.txt");
		OutputStreamWriter out = null;
		final String newLine = System.getProperty("line.separator");
		try
		{
			FileOutputStream fos = new FileOutputStream(file);
			out = new OutputStreamWriter(fos, "UTF-8");
			for (double d : storedDataStatistics.getValues())
			{
				out.write(Double.toString(d));
				out.write(newLine);
			}
		}
		catch (Exception e)
		{
			e.printStackTrace(); // Show the error
		}
		finally
		{
			if (out != null)
			{
				try
				{
					out.close();
				}
				catch (IOException e)
				{
					e.printStackTrace();
				}
			}
		}
	}
}
