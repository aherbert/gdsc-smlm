package gdsc.smlm.ij.plugins.pcpalm;

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

import gdsc.smlm.ij.IJTrackProgress;
import gdsc.smlm.ij.plugins.About;
import gdsc.smlm.ij.plugins.Parameters;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.results.Calibration;
import gdsc.smlm.results.DensityManager;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.clustering.ClusterPoint;
import gdsc.smlm.results.clustering.Cluster;
import gdsc.smlm.results.clustering.ClusteringAlgorithm;
import gdsc.smlm.results.clustering.ClusteringEngine;
import gdsc.smlm.utils.Maths;
import gdsc.smlm.utils.StoredDataStatistics;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.analysis.DifferentiableMultivariateVectorFunction;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.MultivariateMatrixFunction;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.OptimizationData;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.SimpleValueChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer;
import org.apache.commons.math3.optimization.PointVectorValuePair;
import org.apache.commons.math3.optimization.general.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;

/**
 * Find clusters of molecules using a partial centroid-linkage hierarchical clustering algorithm.
 * <p>
 * Points are added to the nearest cluster if they are below the distance threshold to the cluster centroid. The cluster
 * centroid is updated. All points above the cluster distance threshold remain as single molecules.
 * <p>
 * The purpose is to join colocalising molecules into clusters.
 * <p>
 * See Puchnar, et al (2013). Counting molecules in single organelles with superresolution microscopy allows tracking of
 * the endosome maturation trajectory. PNAS. doi:10.1073/pnas.1309676110
 */
public class PCPALMClusters implements PlugInFilter
{
	static String TITLE = "PC-PALM Clusters";

	private static double distance = 50;
	private static ClusteringAlgorithm sClusteringAlgorithm = ClusteringAlgorithm.Pairwise;
	private static int maxN = 0;
	private static boolean showCumulativeHistogram = false;

	private ClusteringAlgorithm clusteringAlgorithm = ClusteringAlgorithm.Pairwise;
	private boolean logging = false;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	public int setup(String arg, ImagePlus imp)
	{
		if (PCPALMMolecules.molecules == null || PCPALMMolecules.molecules.size() < 2)
		{
			error("Require a set of molecules for analysis.\n" + "Please create a binary molecule image using " +
					PCPALMMolecules.TITLE);
			return DONE;
		}
		if (!showDialog())
			return DONE;

		PCPALMMolecules.logSpacer();
		Utils.log(TITLE);
		PCPALMMolecules.logSpacer();

		PCPALMAnalysis analysis = new PCPALMAnalysis();
		ArrayList<Molecule> molecules = analysis.cropToRoi(imp);

		if (molecules.size() < 2)
		{
			error("No results within the crop region");
			return DONE;
		}

		Utils.log("Using %d molecules (Density = %s um^-2)", molecules.size(),
				Utils.rounded(molecules.size() / analysis.croppedArea));
		long start = System.currentTimeMillis();

		long s1 = System.nanoTime();
		this.logging = true;

		ClusteringEngine engine = new ClusteringEngine();
		engine.setClusteringAlgorithm(clusteringAlgorithm);
		engine.setTracker(new IJTrackProgress());
		ArrayList<Cluster> clusters = engine.findClusters(convertToPoint(molecules), distance);

		if (clusters == null)
		{
			Utils.log("Aborted");
			return DONE;
		}
		Utils.log("Finished : %d total clusters (%s ms)", clusters.size(),
				Utils.rounded((System.nanoTime() - s1) / 1e6));

		// Display a histogram of the cluster sizes
		StoredDataStatistics data = new StoredDataStatistics();
		for (Cluster c : clusters)
			data.add(c.n);
		String title = TITLE + " Molecules/cluster";

		// Ensure the bin width is never less than 1
		float yMax = (int) Math.ceil(Maths.max(data.getFloatValues()));
		int nBins = (int) (yMax + 1);
		float[][] hist = Utils.calcHistogram(data.getFloatValues(), 0, yMax, nBins);

		// Create the axes
		float[] xValues = Utils.createHistogramAxis(hist[0]);
		float[] yValues = Utils.createHistogramValues(hist[1]);

		// Plot
		yMax = Maths.max(yValues);
		Plot plot = new Plot(title, "Molecules/cluster", "Frequency", xValues, yValues);
		if (xValues.length > 0)
		{
			double xPadding = 0.05 * (xValues[xValues.length - 1] - xValues[0]);
			plot.setLimits(xValues[0] - xPadding, xValues[xValues.length - 1] + xPadding, 0, yMax * 1.05);
		}
		Utils.display(title, plot);

		double[] fitParameters = fitBinomial(data);
		if (fitParameters != null)
		{
			// Add the binomial to the histogram
			int n = (int) fitParameters[0];
			double p = fitParameters[1];

			Utils.log("Optimal fit : N=%d, p=%s", n, Utils.rounded(p));

			// Calculate the estimate number of clusters from the observed molecules:
			// Actual = (Observed / p-value) / N
			final double actual = (molecules.size() / p) / n;
			Utils.log("Estimated number of clusters : (%d / %s) / %d = %s", molecules.size(), Utils.rounded(p), n,
					Utils.rounded(actual));

			double[] x = new double[n + 2];
			double[] y = new double[n + 2];

			BinomialDistribution dist = new BinomialDistribution(n, p);
			for (int i = 0; i <= n; i++)
			{
				x[i] = i + 0.5;
				y[i] = dist.probability(i) * data.getN();
			}
			x[n + 1] = n + 1.5;
			y[n + 1] = 0;

			// Redraw the plot since the limits may have changed
			plot = new Plot(title, "Molecules/cluster", "Frequency", xValues, yValues);
			double xPadding = 0.05 * (xValues[xValues.length - 1] - xValues[0]);
			plot.setLimits(xValues[0] - xPadding, xValues[xValues.length - 1] + xPadding, 0,
					Maths.maxDefault(yMax, y) * 1.05);
			plot.setColor(Color.magenta);
			plot.addPoints(x, y, Plot.LINE);
			plot.addPoints(x, y, Plot.CIRCLE);
			plot.setColor(Color.black);
			Utils.display(title, plot);
		}

		// Save cluster centroids to a results set in memory. Then they can be plotted.
		MemoryPeakResults results = new MemoryPeakResults(clusters.size());
		results.setName(TITLE);
		// Set an arbitrary calibration so that the lifetime of the results is stored in the exposure time
		// The results will be handled as a single mega-frame containing all localisation. 
		results.setCalibration(new Calibration(100, 1, PCPALMMolecules.seconds * 1000));
		// Make the standard deviation such that the Gaussian volume will be 95% at the distance threshold
		final float sd = (float) (distance / 1.959964);
		for (Cluster c : clusters)
		{
			results.add(new PeakResult((float) c.x, (float) c.y, sd, c.n));
		}
		MemoryPeakResults.addResults(results);

		// TODO - What saving options should we have?
		// Save the clusters to file?

		double seconds = (System.currentTimeMillis() - start) / 1000.0;
		String msg = TITLE + " complete : " + seconds + "s";
		IJ.showStatus(msg);
		Utils.log(msg);
		return DONE;
	}

	private List<ClusterPoint> convertToPoint(ArrayList<Molecule> molecules)
	{
		ArrayList<ClusterPoint> points = new ArrayList<ClusterPoint>(molecules.size());
		int id = 0;
		for (Molecule m : molecules)
		{
			points.add(new ClusterPoint(id++, m.x, m.y));
		}
		return points;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
	 */
	public void run(ImageProcessor ip)
	{
		// Do nothing		
	}

	private boolean showDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Find clusters using centroid-linkage clustering.");

		gd.addNumericField("Distance (nm)", distance, 0);
		String[] names = SettingsManager.getNames((Object[]) ClusteringAlgorithm.values());
		gd.addChoice("Algorithm", names, names[sClusteringAlgorithm.ordinal()]);
		gd.addSlider("Max_N", 0, 10, maxN);
		gd.addCheckbox("Show_cumulative_histogram", showCumulativeHistogram);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		distance = gd.getNextNumber();
		clusteringAlgorithm = sClusteringAlgorithm = ClusteringAlgorithm.values()[gd.getNextChoiceIndex()];
		maxN = (int) Math.abs(gd.getNextNumber());
		showCumulativeHistogram = gd.getNextBoolean();

		// Check arguments
		try
		{
			Parameters.isAboveZero("Distance", distance);
		}
		catch (IllegalArgumentException ex)
		{
			error(ex.getMessage());
			return false;
		}

		return true;
	}

	private void error(String message)
	{
		Utils.log("ERROR : " + message);
		IJ.error(TITLE, message);
	}

	private double[] fitBinomial(StoredDataStatistics data)
	{
		double[][] histogram = Maths.cumulativeHistogram(data.getValues(), true);

		int n = 0;
		double bestSS = Double.POSITIVE_INFINITY;
		double[] parameters = null;
		int worse = 0;
		int N = (int) histogram[0][histogram[0].length - 1];
		if (maxN > 0 && N > maxN)
			N = maxN;
		final double mean = data.getMean();

		Utils.log("Mean cluster size = %s", Utils.rounded(mean));
		Utils.log("Fitting cumulative Binomial distribution excluding X=0 (zero size clusters cannot be observed)");

		// Plot the cumulative histogram
		String title = TITLE + " Cumulative Distribution";
		Plot plot = null;
		if (showCumulativeHistogram)
		{
			plot = new Plot(title, "N", "Cumulative Probability (ex. X=0)", histogram[0], histogram[1]);
			plot.setLimits(0, N, 0, 1.05);
			plot.addPoints(histogram[0], histogram[1], Plot.CIRCLE);
			Utils.display(title, plot);
		}

		// Since varying the N should be done in integer steps do this
		// for n=1,2,3,... until the SS peaks then falls off (is worse then the best 
		// score several times in succession)
		while (n++ < N)
		{
			PointValuePair solution = fitBinomial(histogram, n, mean);
			if (solution == null)
				break;

			double p = solution.getPointRef()[0];

			Utils.log("Fitted Binomial distribution : N=%d, p=%s. SS=%g", n, Utils.rounded(p), solution.getValue());

			if (bestSS > solution.getValue())
			{
				bestSS = solution.getValue();
				parameters = new double[] { n, p };
				worse = 0;
			}
			else if (bestSS < Double.POSITIVE_INFINITY)
			{
				if (++worse >= 3)
					break;
			}

			if (showCumulativeHistogram)
				addToPlot(n, p, title, plot, new Color((float) n / N, 0, 1f - (float) n / N));
		}

		// Add best it in magenta
		if (showCumulativeHistogram && parameters != null)
			addToPlot((int) parameters[0], parameters[1], title, plot, Color.magenta);

		return parameters;
	}

	private void addToPlot(int n, double p, String title, Plot plot, Color color)
	{
		double[] x = new double[n + 1];
		double[] y = new double[n + 1];

		BinomialDistribution dist = new BinomialDistribution(n, p);

		// Normalise to 1 excluding the x=0 point
		double total = 0;
		for (int i = 1; i <= n; i++)
		{
			total += dist.probability(i);
		}

		double cumul = 0;
		for (int i = 1; i <= n; i++)
		{
			cumul += dist.probability(i) / total;
			x[i] = i;
			y[i] = cumul;
		}

		plot.setColor(color);
		plot.addPoints(x, y, Plot.LINE);
		//plot.addPoints(x, y, Plot.CIRCLE);
		Utils.display(title, plot);
	}

	private PointValuePair fitBinomial(double[][] histogram, int n, double mean)
	{
		BinomialModelFunction function = new BinomialModelFunction(histogram, n);

		// The model is only fitting the probability p
		// For a binomial n*p = mean => p = mean/n
		double[] initialSolution = new double[] { Math.min(mean / n, 1) };

		double[] lB = new double[1];
		double[] uB = new double[] { 1 };
		SimpleBounds bounds = new SimpleBounds(lB, uB);

		// Fit
		// CMAESOptimizer or BOBYQAOptimizer support bounds

		// CMAESOptimiser based on Matlab code:
		// https://www.lri.fr/~hansen/cmaes.m
		// Take the defaults from the Matlab documentation
		int maxIterations = 2000;
		double stopFitness = 0; //Double.NEGATIVE_INFINITY;
		boolean isActiveCMA = true;
		int diagonalOnly = 0;
		int checkFeasableCount = 1;
		RandomGenerator random = new Well19937c();
		boolean generateStatistics = false;
		ConvergenceChecker<PointValuePair> checker = new SimpleValueChecker(1e-6, 1e-10);
		// The sigma determines the search range for the variables. It should be 1/3 of the initial search region.
		OptimizationData sigma = new CMAESOptimizer.Sigma(new double[] { (uB[0] - lB[0]) / 3 });
		OptimizationData popSize = new CMAESOptimizer.PopulationSize((int) (4 + Math.floor(3 * Math
				.log(function.p.length))));

		try
		{
			CMAESOptimizer opt = new CMAESOptimizer(maxIterations, stopFitness, isActiveCMA, diagonalOnly,
					checkFeasableCount, random, generateStatistics, checker);
			PointValuePair solution = opt.optimize(new InitialGuess(initialSolution), new ObjectiveFunction(function),
					GoalType.MINIMIZE, bounds, sigma, popSize, new MaxIter(maxIterations), new MaxEval(
							maxIterations * 2));
			if (solution == null)
				return null;

			// Improve fit with a gradient based LVM optimizer
			LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer();
			try
			{
				BinomialModelFunctionGradient gradientFunction = new BinomialModelFunctionGradient(histogram, n);
				PointVectorValuePair lvmSolution = optimizer.optimize(3000, gradientFunction, gradientFunction.p,
						gradientFunction.getWeights(), solution.getPoint());

				double ss = 0;
				double[] obs = gradientFunction.p;
				double[] exp = lvmSolution.getValue();
				for (int i = 0; i < obs.length; i++)
					ss += (obs[i] - exp[i]) * (obs[i] - exp[i]);
				// Check the pValue is valid since the LVM is not bounded.
				double pValue = lvmSolution.getPointRef()[0];
				if (ss < solution.getValue() && pValue <= 1 && pValue >= 0)
				{
					//Utils.log("Re-fitting improved the SS from %s to %s (-%s%%)",
					//		Utils.rounded(solution.getValue(), 4), Utils.rounded(ss, 4),
					//		Utils.rounded(100 * (solution.getValue() - ss) / solution.getValue(), 4));
					return new PointValuePair(lvmSolution.getPoint(), ss);
				}
			}
			catch (TooManyEvaluationsException e)
			{
				Utils.log("Failed to re-fit: Too many evaluations (%d)", optimizer.getEvaluations());
			}
			catch (ConvergenceException e)
			{
				Utils.log("Failed to re-fit: %s", e.getMessage());
			}
			catch (Exception e)
			{
				// Ignore this ...
			}

			return solution;
		}
		catch (Exception e)
		{
			Utils.log("Failed to fit Binomial distribution with N=%d : %s", n, e.getMessage());
		}
		return null;
	}

	/**
	 * Evaluates the cumulative binomial probability distribution. Assumes the
	 * input data is a cumulative histogram from 0 to N in integer increments.
	 */
	public class BinomialModel
	{
		int trials;
		double[] p;

		public BinomialModel(double[][] histogram, int trials)
		{
			this.trials = trials;

			double[] nValues = histogram[0];
			double[] pValues = histogram[1];
			int N = (int) nValues[nValues.length - 1];
			p = new double[N + 1];

			// Pad the histogram out for any missing values between 0 and N
			for (int i = 1; i < nValues.length; i++)
			{
				int j = (int) nValues[i - 1];
				int k = (int) nValues[i];
				for (int ii = j; ii < k; ii++)
					p[ii] = pValues[i - 1];
			}
			p[N] = pValues[pValues.length - 1];
		}

		/**
		 * Get the cumulative probability function for the input pValue ignoring the X=0 data point
		 * 
		 * @param pValue
		 * @return
		 */
		public double[] getP(double pValue)
		{
			BinomialDistribution dist = new BinomialDistribution(trials, pValue);

			double cumul = 0;

			// Ignore x=0 since we cannot see a zero size cluster.
			// This is done by re-normalising the cumulative probability excluding x=0 
			// to match the input curve.
			double[] p2 = new double[p.length];
			for (int i = 1; i < p.length; i++)
			{
				cumul += dist.probability(i);
				p2[i] = cumul;
			}

			for (int i = 1; i < p.length; i++)
			{
				p2[i] /= cumul;
			}

			return p2;
		}
	}

	/**
	 * Allow optimisation using Apache Commons Math 3 MultivariateFunction optimisers
	 */
	public class BinomialModelFunction extends BinomialModel implements MultivariateFunction
	{
		public BinomialModelFunction(double[][] histogram, int trials)
		{
			super(histogram, trials);
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see org.apache.commons.math3.analysis.MultivariateFunction#value(double[])
		 */
		public double value(double[] parameters)
		{
			double[] p2 = getP(parameters[0]);
			double ss = 0;
			for (int i = 1; i < p.length; i++)
			{
				final double dx = p[i] - p2[i];
				ss += dx * dx;
			}
			return ss;
		}
	}

	/**
	 * Allow optimisation using Apache Commons Math 3 MultivariateFunction optimisers
	 */
	public class BinomialModelFunctionGradient extends BinomialModel implements
			DifferentiableMultivariateVectorFunction
	{
		public BinomialModelFunctionGradient(double[][] histogram, int trials)
		{
			super(histogram, trials);

			// We could ignore the first p value as it is always zero:
			//p = Arrays.copyOfRange(p, 1, p.length);
			// BUT then we would have to override the getP() method since this has 
			// an offset of 1 and assumes the index of p is X.
		}

		public double[] getWeights()
		{
			double[] w = new double[p.length];
			Arrays.fill(w, 1);
			return w;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see org.apache.commons.math3.analysis.MultivariateFunction#value(double[])
		 */
		public double[] value(double[] point) throws IllegalArgumentException
		{
			return getP(point[0]);
		}

		public MultivariateMatrixFunction jacobian()
		{
			return new MultivariateMatrixFunction()
			{
				public double[][] value(double[] variables)
				{
					return jacobian(variables);
				}
			};
		}

		double[][] jacobian(double[] variables)
		{
			// Compute the gradients using numerical differentiation
			final double pValue = variables[0];
			double[][] jacobian = new double[p.length][1];

			final double delta = 0.001 * pValue;
			double[] p2 = getP(pValue);
			double[] p3 = getP(pValue + delta);

			for (int i = 0; i < jacobian.length; ++i)
			{
				jacobian[i][0] = (p3[i] - p2[i]) / delta;
			}
			return jacobian;
		}
	}
}
