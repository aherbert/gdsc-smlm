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

import gdsc.smlm.fitting.BinomialFitter;
import gdsc.smlm.ij.IJTrackProgress;
import gdsc.smlm.ij.plugins.About;
import gdsc.smlm.ij.plugins.Parameters;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.IJLogger;
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.results.Calibration;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.clustering.Cluster;
import gdsc.smlm.results.clustering.ClusterPoint;
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
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.optim.PointValuePair;

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

		int mode = 3;

		// TODO - Depending on the mode we need to add the fictional n=0 data point (in blue)

		double[] fitParameters = fitBinomial(data, mode);
		if (fitParameters != null)
		{
			// Add the binomial to the histogram
			int n = (int) fitParameters[0];
			double p = fitParameters[1];

			Utils.log("Optimal fit : N=%d, p=%s", n, Utils.rounded(p));

			BinomialDistribution dist = new BinomialDistribution(n, p);

			// Depending on the mode either a standard or a zero-truncated binomial was fitted.
			// pi is the adjustment factor for the probability density.
			double pi = (mode == 3) ? 1 / (1 - dist.probability(0)) : 1;

			// Calculate the estimate number of clusters from the observed molecules:
			// Actual = (Observed / p-value) / N
			final double actual = (molecules.size() / p) / n;
			Utils.log("Estimated number of clusters : (%d / %s) / %d = %s", molecules.size(), Utils.rounded(p), n,
					Utils.rounded(actual));
			// This must be adjusted if fitting a zero truncated distribution.
			Utils.log("Estimated total number of clusters : %s * (%d / %s) / %d = %s", Utils.rounded(pi),
					molecules.size(), Utils.rounded(p), n, Utils.rounded(pi * actual));

			double[] x = new double[n + 2];
			double[] y = new double[n + 2];

			for (int i = 0; i <= n; i++)
			{
				x[i] = i + 0.5;
				y[i] = dist.probability(i) * data.getN() * pi;
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

	private double[] fitBinomial(StoredDataStatistics data, int mode)
	{
		double[] cumulativeHistogram = BinomialFitter.getHistogram(data.getValues(), true);

		double bestSS = Double.POSITIVE_INFINITY;
		double[] parameters = null;
		int worse = 0;
		int N = (int) cumulativeHistogram.length - 1;
		double[] values = Utils.newArray(cumulativeHistogram.length, 0.0, 1.0);
		int minN = 1;
		if (maxN > 0 && N > maxN)
			N = maxN;
		final double mean = data.getMean();

		String name = (mode == 3) ? "Zero-truncated Binomial distribution" : "Binomial distribution";

		Utils.log("Mean cluster size = %s", Utils.rounded(mean));
		Utils.log("Fitting cumulative " + name);

		// Plot the cumulative histogram
		String title = TITLE + " Cumulative Distribution";
		Plot plot = null;
		if (showCumulativeHistogram)
		{
			plot = new Plot(title, "N", "Cumulative Probability", values, cumulativeHistogram);
			plot.setLimits(0, N, 0, 1.05);
			plot.addPoints(values, cumulativeHistogram, Plot.CIRCLE);
			Utils.display(title, plot);
		}

		// We need the original histogram, not the cumulative histogram
		double[] histogram = new double[cumulativeHistogram.length];
		for (int i = histogram.length; i-- > 1;)
		{
			histogram[i] = cumulativeHistogram[i] - cumulativeHistogram[i - 1];
		}

		// Since varying the N should be done in integer steps do this
		// for n=1,2,3,... until the SS peaks then falls off (is worse then the best 
		// score several times in succession)
		BinomialFitter bf = new BinomialFitter(new IJLogger());
		bf.setMaximumLikelihood(false);
		for (int n = minN; n <= N; n++)
		{
			// Optionally iterate this if we are estimating the n=0 data using an initial p value.
			// TODO - Add options for the mode input and the initial p-value
			PointValuePair solution = bf.fitBinomial(histogram, n, mean, mode, 0);
			if (solution == null)
				break;

			double p = solution.getPointRef()[0];

			Utils.log("Fitted %s : N=%d, p=%s. SS=%g", name, n, Utils.rounded(p), solution.getValue());

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
				addToPlot(n, p, mode, title, plot, new Color((float) n / N, 0, 1f - (float) n / N));
		}

		// Add best it in magenta
		if (showCumulativeHistogram && parameters != null)
			addToPlot((int) parameters[0], parameters[1], mode, title, plot, Color.magenta);

		return parameters;
	}

	private void addToPlot(int n, double p, int mode, String title, Plot plot, Color color)
	{
		double[] x = new double[n + 1];
		double[] y = new double[n + 1];

		BinomialDistribution dist = new BinomialDistribution(n, p);

		int startIndex = (mode == 3) ? 1 : 0;

		// Normalise optionally excluding the x=0 point
		double total = 1;
		if (startIndex > 0)
			total -= dist.probability(0);

		double cumul = 0;
		for (int i = startIndex; i <= n; i++)
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
}
