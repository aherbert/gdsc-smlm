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
import gdsc.smlm.utils.UnicodeReader;
import ij.IJ;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.InputMismatchException;
import java.util.List;
import java.util.Locale;
import java.util.NoSuchElementException;
import java.util.Scanner;
import java.util.regex.Pattern;

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
public class PCPALMClusters implements PlugIn
{
	private class HistogramData
	{
		double[][] histogram;
		int frames;
		double area;
		String units;

		public HistogramData(double[][] h, int f, double a, String u)
		{
			histogram = h;
			frames = f;
			area = a;
			units = u;
		}
		
		public boolean isCalibrated()
		{
			return frames > 0 && area > 0;
		}
	}

	static String TITLE = "PC-PALM Clusters";

	private static int runMode = 0;
	private static double distance = 50;
	private static ClusteringAlgorithm sClusteringAlgorithm = ClusteringAlgorithm.ClosestParticle;
	private static int minN = 1;
	private static int maxN = 0;
	private static boolean maximumLikelihood = false;
	private static boolean showCumulativeHistogram = false;
	private static boolean multiThread = true;
	private static boolean sWeightedClustering = false;
	private static boolean saveHistogram = false;
	private static String histogramFile = "";
	private static boolean calibrateHistogram = false;
	private static int frames = 1;
	private static double area = 1;
	private static String[] UNITS = { "pixels^2", "um^2" };
	private static String units = UNITS[0];

	private ClusteringAlgorithm clusteringAlgorithm = ClusteringAlgorithm.ClosestParticle;
	private boolean weightedClustering = false;

	private boolean fileInput = false;

	// Set by the loadHistogram method
	String frequencyTitle = "";

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	@Override
	public void run(String arg)
	{
		if (!showDialog())
			return;

		PCPALMMolecules.logSpacer();
		Utils.log(TITLE);
		PCPALMMolecules.logSpacer();
		long start = System.currentTimeMillis();

		// Create a histogram of the cluster sizes
		int observed = 1;
		double count;
		double[] cumulativeHistogram;
		double mean;
		String title = TITLE + " Molecules/cluster";
		String xTitle = "Molecules/cluster";
		String yTitle = "Frequency";
		float[] xValues, yValues;

		if (fileInput)
		{
			// Load the histogram from file			
			HistogramData histogramData = loadHistogram(histogramFile);
			if (histogramData == null)
				return;
			double[][] hist = histogramData.histogram;

			// Create the data required for fitting and plotting
			xValues = convert(Utils.createHistogramAxis(hist[0]));
			yValues = convert(Utils.createHistogramValues(hist[1]));
			double sum = 0;
			count = 0;
			for (int i = 0; i < hist[0].length; i++)
			{
				count += hist[1][i];
				sum += hist[1][i] * i;
			}
			mean = sum / count;
			cumulativeHistogram = new double[hist[0].length];
			final double norm = 1.0 / count;
			count = 0;
			for (int i = 0; i < hist[0].length; i++)
			{
				count += hist[1][i];
				cumulativeHistogram[i] = count * norm;
			}
		}
		else
		{
			// Perform clustering analysis to generate the histogram of cluster sizes
			PCPALMAnalysis analysis = new PCPALMAnalysis();
			ArrayList<Molecule> molecules = analysis.cropToRoi(WindowManager.getCurrentImage());

			if (molecules.size() < 2)
			{
				error("No results within the crop region");
				return;
			}

			Utils.log("Using %d molecules (Density = %s um^-2)", molecules.size(),
					Utils.rounded(molecules.size() / analysis.croppedArea));

			long s1 = System.nanoTime();
			ClusteringEngine engine = new ClusteringEngine();
			engine.setClusteringAlgorithm(clusteringAlgorithm);
			engine.setTracker(new IJTrackProgress());
			if (multiThread)
				engine.setThreadCount(Prefs.getThreads());
			IJ.showStatus("Clustering ...");
			ArrayList<Cluster> clusters = engine.findClusters(convertToPoint(molecules), distance);
			IJ.showStatus("");

			if (clusters == null)
			{
				Utils.log("Aborted");
				return;
			}
			observed = molecules.size();
			Utils.log("Finished : %d total clusters (%s ms)", clusters.size(),
					Utils.rounded((System.nanoTime() - s1) / 1e6));

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

			// Get the data for fitting
			StoredDataStatistics data = new StoredDataStatistics();
			for (Cluster c : clusters)
				data.add(c.n);
			cumulativeHistogram = BinomialFitter.getHistogram(data.getValues(), true);
			mean = data.getMean();

			float[] values = data.getFloatValues();
			float yMax = (int) Math.ceil(Maths.max(values));
			int nBins = (int) (yMax + 1);
			float[][] hist = Utils.calcHistogram(values, 0, yMax, nBins);

			saveHistogram(hist, xTitle, yTitle);

			// Create the axes
			xValues = Utils.createHistogramAxis(hist[0]);
			yValues = Utils.createHistogramValues(hist[1]);
			count = clusters.size();
		}

		// TODO - option to load a noise histogram and subtract it before fitting.
		// Either the clustering algorithm was run and the calibration was provided
		// Or the histogram was loaded and we can detect if it is calibrated
		// If the histogram is calibrated then ask the user if they wish to load a calibrated noise histogram
		
		// Updating the histogram will require:
		// - replotting the noise adjusted histogram (leave the raw one)
		// - the count to be updated
		
		// TODO make both routine above return a HistogramData object. This allows common processing from here.
		
		
		// Plot the (noise adjusted) histogram
		float yMax = Maths.max(yValues);
		Plot plot = new Plot(title, xTitle, yTitle, xValues, yValues);
		if (xValues.length > 0)
		{
			double xPadding = 0.05 * (xValues[xValues.length - 1] - xValues[0]);
			plot.setLimits(xValues[0] - xPadding, xValues[xValues.length - 1] + xPadding, 0, yMax * 1.05);
		}
		Utils.display(title, plot);
		
		// Fit the histogram
		double[] fitParameters = fitBinomial(cumulativeHistogram, mean);
		if (fitParameters != null)
		{
			// Add the binomial to the histogram
			int n = (int) fitParameters[0];
			double p = fitParameters[1];

			Utils.log("Optimal fit : N=%d, p=%s", n, Utils.rounded(p));

			BinomialDistribution dist = new BinomialDistribution(n, p);

			// A zero-truncated binomial was fitted.
			// pi is the adjustment factor for the probability density.
			double pi = 1 / (1 - dist.probability(0));

			if (!fileInput)
			{
				// Calculate the estimated number of clusters from the observed molecules:
				// Actual = (Observed / p-value) / N
				final double actual = (observed / p) / n;
				Utils.log("Estimated number of clusters : (%d / %s) / %d = %s", (int) observed, Utils.rounded(p), n,
						Utils.rounded(actual));
			}

			double[] x = new double[n + 2];
			double[] y = new double[n + 2];

			// Scale the values to match those on the histogram
			final double normalisingFactor = count * pi;
			for (int i = 0; i <= n; i++)
			{
				x[i] = i + 0.5;
				y[i] = dist.probability(i) * normalisingFactor;
			}
			x[n + 1] = n + 1.5;
			y[n + 1] = 0;

			// Redraw the plot since the limits may have changed
			plot = new Plot(title, xTitle, yTitle, xValues, yValues);
			double xPadding = 0.05 * (xValues[xValues.length - 1] - xValues[0]);
			plot.setLimits(xValues[0] - xPadding, xValues[xValues.length - 1] + xPadding, 0,
					Maths.maxDefault(yMax, y) * 1.05);
			plot.setColor(Color.magenta);
			plot.addPoints(x, y, Plot.LINE);
			plot.addPoints(x, y, Plot.CIRCLE);
			plot.setColor(Color.black);
			Utils.display(title, plot);
		}

		double seconds = (System.currentTimeMillis() - start) / 1000.0;
		String msg = TITLE + " complete : " + seconds + "s";
		IJ.showStatus(msg);
		Utils.log(msg);
		return;
	}

	private float[] convert(double[] values)
	{
		float[] data = new float[values.length];
		for (int i = 0; i < data.length; i++)
			data[i] = (float) values[i];
		return data;
	}

	private void saveHistogram(float[][] hist, String xTitle, String yTitle)
	{
		if (!saveHistogram)
			return;
		histogramFile = Utils.getFilename("Histogram_file", histogramFile);
		if (histogramFile != null)
		{
			float[] values = hist[1];
			//if (normaliseHistogram)
			//{
			//	values = Arrays.copyOf(values, values.length);
			//	yTitle += "/(frame*" + units + ")";
			//	final double normalisingFactor = 1.0 / (frames * area);
			//	for (int i = 0; i < hist[1].length; i++)
			//		values[i] *= normalisingFactor;
			//}

			histogramFile = Utils.replaceExtension(histogramFile, "tsv");
			BufferedWriter output = null;
			try
			{
				output = new BufferedWriter(new FileWriter(histogramFile));
				if (calibrateHistogram)
				{
					output.write(String.format("Frames  %d", frames));
					output.newLine();
					output.write(String.format("Area    %f", area));
					output.newLine();
					output.write(String.format("Units   %s", units));
					output.newLine();
				}
				output.write(xTitle);
				output.write("\t");
				output.write(yTitle);
				output.newLine();
				for (int i = 0; i < hist[0].length; i++)
				{
					output.write(String.format("%d\t%s", (int) hist[0][i], Utils.rounded(values[i])));
					output.newLine();
				}
			}
			catch (Exception e)
			{
				e.printStackTrace();
				IJ.log("Failed to save histogram to file: " + histogramFile);
			}
			finally
			{
				if (output != null)
				{
					try
					{
						output.close();
					}
					catch (IOException e)
					{
						e.printStackTrace();
					}
				}
			}
		}
	}

	private HistogramData loadHistogram(String filename)
	{
		frequencyTitle = "";
		BufferedReader input = null;
		try
		{
			int f = 0;
			double a = 0;
			String u = "";

			FileInputStream fis = new FileInputStream(filename);
			input = new BufferedReader(new UnicodeReader(fis, null));

			String line;
			int count = 0;

			ArrayList<double[]> data = new ArrayList<double[]>();

			// Read the header and store the calibration if present
			while ((line = input.readLine()) != null)
			{
				count++;
				if (line.length() == 0)
					continue;
				if (Character.isDigit(line.charAt(0)))
					// This is the first record
					break;
				String[] fields = line.split("[\t, ]+");
				if (fields[0].equalsIgnoreCase("frames"))
					f = Integer.parseInt(fields[1]);
				if (fields[0].equalsIgnoreCase("area"))
					a = Double.parseDouble(fields[1]);
				if (fields[0].equalsIgnoreCase("units"))
					u = fields[1];
			}

			final Pattern pattern = Pattern.compile("[\t, ]+");
			while (line != null)
			{
				if (line.length() == 0)
					continue;
				if (!Character.isDigit(line.charAt(0)))
					continue;

				// Extract the first 2 fields
				Scanner scanner = new Scanner(line);
				scanner.useLocale(Locale.US);
				scanner.useDelimiter(pattern);

				try
				{
					int molecules = scanner.nextInt();
					double frequency = scanner.nextDouble();

					// Check for duplicates
					for (double[] d : data)
					{
						if (d[0] == molecules)
						{
							error("Duplicate molecules field on line " + count);
							return null;
						}
					}

					data.add(new double[] { molecules, frequency });
				}
				catch (InputMismatchException e)
				{
					error("Incorrect fields on line " + count);
					return null;
				}
				catch (NoSuchElementException e)
				{
					error("Incorrect fields on line " + count);
					return null;
				}
				finally
				{
					scanner.close();
				}

				// Get the next line
				line = input.readLine();
				count++;
			}

			if (data.isEmpty())
			{
				error("No data in file " + filename);
				return null;
			}

			// Create a contiguous histogram from zero
			int maxN = 0;
			for (double[] d : data)
			{
				if (maxN < d[0])
					maxN = (int) d[0];
			}

			double[][] hist = new double[2][maxN + 1];
			for (int n = 0; n <= maxN; n++)
			{
				hist[0][n] = n;
				for (double[] d : data)
				{
					if (n == d[0])
						hist[1][n] = d[1];
				}
			}
			return new HistogramData(hist, f, a, u);
		}
		catch (IOException e)
		{
			IJ.error(TITLE, "Unable to read from file " + filename);
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
		return null;
	}

	private List<ClusterPoint> convertToPoint(ArrayList<Molecule> molecules)
	{
		ArrayList<ClusterPoint> points = new ArrayList<ClusterPoint>(molecules.size());
		int id = 0;
		for (Molecule m : molecules)
		{
			points.add(ClusterPoint.newClusterPoint(id++, m.x, m.y, (weightedClustering) ? m.photons : 1));
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
		if (PCPALMMolecules.molecules == null || PCPALMMolecules.molecules.size() < 2)
		{
			Utils.log(TITLE + " defaulting to File mode");
			fileInput = true;
		}
		else
		{
			GenericDialog gd = new GenericDialog(TITLE);
			String[] items = { "Clustering", "File" };

			gd.addMessage("Fit a Binomial distribution to a histogram of cluster sizes.\n \nSelect the method to generate the histogram:");
			gd.addChoice("Method", items, items[runMode]);
			gd.showDialog();
			if (gd.wasCanceled())
				return false;
			runMode = gd.getNextChoiceIndex();
			fileInput = (runMode == 1);
		}

		if (fileInput)
		{
			if ((histogramFile = Utils.getFilename("Histogram_file", histogramFile)) == null)
				return false;
		}

		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		// Check if the molecules have weights
		boolean haveWeights = false;
		if (!fileInput)
		{
			haveWeights = checkForWeights();

			gd.addMessage("Find clusters using centroid-linkage clustering.");

			gd.addNumericField("Distance (nm)", distance, 0);
			String[] names = SettingsManager.getNames((Object[]) ClusteringAlgorithm.values());
			gd.addChoice("Algorithm", names, names[sClusteringAlgorithm.ordinal()]);
			gd.addCheckbox("Multi_thread", multiThread);
			if (haveWeights)
				gd.addCheckbox("Weighted_clustering", sWeightedClustering);
		}

		gd.addSlider("Min_N", 1, 10, minN);
		gd.addSlider("Max_N", 0, 10, maxN);
		gd.addCheckbox("Show_cumulative_histogram", showCumulativeHistogram);
		gd.addCheckbox("Maximum_likelihood", maximumLikelihood);

		if (!fileInput)
		{
			gd.addCheckbox("Save_histogram", saveHistogram);
			gd.addMessage("Histogram calibration (optional)");
			gd.addCheckbox("Calibrate_histogram", calibrateHistogram);
			gd.addNumericField("Frames", frames, 0);
			gd.addNumericField("Area", area, 2);
			gd.addChoice("Units", UNITS, units);
		}

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		if (!fileInput)
		{
			distance = gd.getNextNumber();
			clusteringAlgorithm = sClusteringAlgorithm = ClusteringAlgorithm.values()[gd.getNextChoiceIndex()];
			multiThread = gd.getNextBoolean();
			if (haveWeights)
				weightedClustering = sWeightedClustering = gd.getNextBoolean();
		}
		minN = (int) Math.abs(gd.getNextNumber());
		maxN = (int) Math.abs(gd.getNextNumber());
		showCumulativeHistogram = gd.getNextBoolean();
		maximumLikelihood = gd.getNextBoolean();
		if (!fileInput)
		{
			saveHistogram = gd.getNextBoolean();
			calibrateHistogram = gd.getNextBoolean();
			frames = (int) Math.abs(gd.getNextNumber());
			area = Math.abs(gd.getNextNumber());
			units = gd.getNextChoice();
		}

		// Check arguments
		try
		{
			Parameters.isAboveZero("Min N", minN);
			if (!fileInput)
			{
				Parameters.isAboveZero("Distance", distance);
				Parameters.isAboveZero("Frames", frames);
				Parameters.isAboveZero("Area", area);
			}
		}
		catch (IllegalArgumentException ex)
		{
			error(ex.getMessage());
			return false;
		}

		return true;
	}

	private boolean checkForWeights()
	{
		for (Molecule m : PCPALMMolecules.molecules)
			if (m.photons <= 0)
				return false;
		return true;
	}

	private void error(String message)
	{
		Utils.log("ERROR : " + message);
		IJ.error(TITLE, message);
	}

	private double[] fitBinomial(double[] cumulativeHistogram, double mean)
	{
		double bestSS = Double.POSITIVE_INFINITY;
		double[] parameters = null;
		int worse = 0;
		int N = (int) cumulativeHistogram.length - 1;
		double[] values = Utils.newArray(cumulativeHistogram.length, 0.0, 1.0);
		int min = minN;
		if (maxN > 0 && N > maxN)
			N = maxN;
		if (min > N)
			min = N;

		String name = "Zero-truncated Binomial distribution";

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
		bf.setMaximumLikelihood(maximumLikelihood);
		for (int n = min; n <= N; n++)
		{
			PointValuePair solution = bf.fitBinomial(histogram, mean, n, true);
			if (solution == null)
				continue;

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

		int startIndex = 1;

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
