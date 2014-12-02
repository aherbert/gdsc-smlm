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

import gdsc.smlm.ij.plugins.About;
import gdsc.smlm.ij.plugins.Parameters;
import gdsc.smlm.ij.utils.LoggingOptimiserFunction;
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.utils.Maths;
import gdsc.smlm.utils.Statistics;
import gdsc.smlm.utils.UnicodeReader;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import ij.text.TextWindow;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.InputMismatchException;
import java.util.NoSuchElementException;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.MultivariateMatrixFunction;
import org.apache.commons.math3.analysis.MultivariateVectorFunction;
import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.exception.TooManyIterationsException;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.optim.GradientChecker;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.OptimizationData;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.PointVectorValuePair;
import org.apache.commons.math3.optim.PositionChecker;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.SimpleValueChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunctionGradient;
import org.apache.commons.math3.optim.nonlinear.scalar.gradient.BFGSOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer;
import org.apache.commons.math3.optim.nonlinear.vector.ModelFunction;
import org.apache.commons.math3.optim.nonlinear.vector.ModelFunctionJacobian;
import org.apache.commons.math3.optim.nonlinear.vector.Target;
import org.apache.commons.math3.optim.nonlinear.vector.Weight;
import org.apache.commons.math3.optim.nonlinear.vector.jacobian.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well44497b;
import org.apache.commons.math3.util.FastMath;

/**
 * Use the PC-PALM protocol to fit correlation curve(s) using the random or clustered model.
 * <p>
 * See Sengupta, et al (2013). Quantifying spatial resolution in point-localisation superresolution images using pair
 * correlation analysis. Nature Protocols 8, pp345-354.
 */
public class PCPALMFitting implements PlugIn
{
	static String TITLE = "PC-PALM Fitting";

	private static String INPUT_FROM_FILE = "Load from file";
	private static String INPUT_PREVIOUS = "Re-use previous curve";
	private static String INPUT_ANALYSIS = "Select PC-PALM Analysis results";
	private static String HEADER_PEAK_DENSITY = "Peak density (um^-2)";
	private static String HEADER_SPATIAL_DOMAIN = "Spatial domain";

	private static String inputOption = "";
	private static double correlationDistance = 800; // nm
	private static double estimatedPrecision = -1;
	private static double copiedEstimatedPrecision = -1;
	private static double blinkingRate = -1;
	private static double copiedBlinkingRate = -1;
	private static boolean showErrorBars = false;
	private static boolean fitClusteredModels = false;
	private static boolean saveCorrelationCurve = false;
	private static String inputFilename = "";
	private static String outputFilename = "";
	private static int fitRestarts = 3;
	private static boolean useLSE = false;
	private static double fitAboveEstimatedPrecision = 0;
	private static double fittingTolerance = 0; // Zero to ignore
	private static double gr_protein_threshold = 1.5;

	private RandomModelFunction randomModel;
	private ClusteredModelFunctionGradient clusteredModel;
	private EmulsionModelFunctionGradient emulsionModel;

	private int boundedEvaluations;

	// Information criterion of models
	private double ic1, ic2, ic3;
	private boolean valid1, valid2;

	// Used for the results table
	private static TextWindow resultsTable = null;

	private boolean doneHeader = false;

	private int offset = 0;
	private double[][] gr;
	private double peakDensity;
	private boolean spatialDomain;

	// Save the input for the analysis
	static double[][] previous_gr;
	private static double previous_peakDensity;
	private static boolean previous_spatialDomain;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		//		if (PCPALMAnalysis.results.isEmpty())
		//		{
		//			IJ.error(TITLE, "Require a set of correlation curves for analysis.\n" +
		//					"Please create a g(r) curve using " + PCPALMAnalysis.TITLE);
		//			return;
		//		}
		if (!showDialog())
			return;

		long start = System.currentTimeMillis();
		header();

		analyse();

		double seconds = (System.currentTimeMillis() - start) / 1000.0;
		String msg = TITLE + " complete : " + seconds + "s";
		IJ.showStatus(msg);
		log(msg);
	}

	private void header()
	{
		if (!doneHeader)
		{
			doneHeader = true;
			PCPALMMolecules.logSpacer();
			log(TITLE);
			PCPALMMolecules.logSpacer();
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
	 */
	public void run(ImageProcessor ip)
	{
		// Nothing to do
	}

	private boolean showDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		// Build a list of results to use in the analysis
		if (!getCorrelationResults())
			return false;

		if (spatialDomain)
		{
			// Spatial domain results are just combined to a curve
			// Add option to save the results curve
			gd.addMessage("Options:");
			gd.addCheckbox("Save_correlation_curve", saveCorrelationCurve);
			gd.showDialog();
			if (gd.wasCanceled())
				return false;
			saveCorrelationCurve = gd.getNextBoolean();
			return true;
		}

		if (estimatedPrecision < 0 || copiedEstimatedPrecision != PCPALMMolecules.sigmaS)
		{
			copiedEstimatedPrecision = estimatedPrecision = PCPALMMolecules.sigmaS;
		}
		if (blinkingRate < 0 || copiedBlinkingRate != PCPALMAnalysis.blinkingRate)
		{
			copiedBlinkingRate = blinkingRate = PCPALMAnalysis.blinkingRate;
		}

		gd.addMessage("Analyse clusters using Pair Correlation.");

		gd.addNumericField("Estimated_precision", estimatedPrecision, 2);
		gd.addNumericField("Blinking_rate", blinkingRate, 2);
		gd.addCheckbox("Show_error_bars", showErrorBars);
		gd.addSlider("Fit_restarts", 0, 5, fitRestarts);
		gd.addCheckbox("Refit_using_LSE", useLSE);
		gd.addSlider("Fit_above_estimated_precision", 0, 2.5, fitAboveEstimatedPrecision);
		gd.addSlider("Fitting_tolerance", 0, 200, fittingTolerance);
		gd.addSlider("gr_random_threshold", 1, 2.5, gr_protein_threshold);
		gd.addCheckbox("Fit_clustered_models", fitClusteredModels);
		gd.addCheckbox("Save_correlation_curve", saveCorrelationCurve);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		estimatedPrecision = gd.getNextNumber();
		blinkingRate = gd.getNextNumber();
		showErrorBars = gd.getNextBoolean();
		fitRestarts = (int) Math.abs(gd.getNextNumber());
		useLSE = gd.getNextBoolean();
		fitAboveEstimatedPrecision = Math.abs(gd.getNextNumber());
		fittingTolerance = Math.abs(gd.getNextNumber());
		gr_protein_threshold = gd.getNextNumber();
		fitClusteredModels = gd.getNextBoolean();
		saveCorrelationCurve = gd.getNextBoolean();

		// Check arguments
		try
		{
			Parameters.isAbove("Correlation distance", correlationDistance, 1);
			Parameters.isAbove("Estimated precision", estimatedPrecision, 0);
			Parameters.isAbove("Blinking_rate", blinkingRate, 0);
			Parameters.isAbove("gr random threshold", gr_protein_threshold, 1);
		}
		catch (IllegalArgumentException ex)
		{
			IJ.error(TITLE, ex.getMessage());
			return false;
		}

		return true;
	}

	private boolean getCorrelationResults()
	{
		// Option to:
		// - load a correlation curve
		// - use previous results (if available)
		// - select a set of analysis results (if available)
		String[] options = new String[] { INPUT_FROM_FILE, "", "" };
		int count = 1;
		if (previous_gr != null)
			options[count++] = INPUT_PREVIOUS;
		if (!PCPALMAnalysis.results.isEmpty())
			options[count++] = INPUT_ANALYSIS;

		options = Arrays.copyOf(options, count);
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addMessage("Select the source for the correlation curve");
		gd.addChoice("Input", options, inputOption);
		gd.showDialog();
		if (gd.wasCanceled())
			return false;
		inputOption = gd.getNextChoice();

		if (inputOption.equals(INPUT_PREVIOUS))
		{
			// In the case of a macro the previous results may be null
			if (previous_gr == null)
				return false;

			gr = previous_gr;
			peakDensity = previous_peakDensity;
			spatialDomain = previous_spatialDomain;
			return true;
		}
		else if (inputOption.equals(INPUT_FROM_FILE))
		{
			return loadCorrelationCurve();
		}

		// Fill the results list with analysis results from PCPALM Analysis
		ArrayList<CorrelationResult> results = new ArrayList<CorrelationResult>();
		if (!selectAnalysisResults(results))
			return false;

		// We have some results. Convert them to the format used for fitting.

		header();
		log("Computing combined pair correlation curve (%d datasets)", results.size());

		spatialDomain = results.get(0).spatialDomain;

		// Get average peak density
		peakDensity = 0;

		int size = 0;
		for (CorrelationResult r : results)
		{
			peakDensity += r.peakDensity;
			size = FastMath.max(size, r.gr[0].length);
		}
		peakDensity /= results.size();

		// Combine all selected g(r) curves
		gr = combineCurves(results, size);

		return true;
	}

	private boolean selectAnalysisResults(ArrayList<CorrelationResult> results)
	{
		// If no results then fail
		if (PCPALMAnalysis.results.isEmpty())
			return false;

		// If only one result then use that
		if (PCPALMAnalysis.results.size() == 1)
		{
			results.add(PCPALMAnalysis.results.get(0));
			return true;
		}

		// Otherwise build a set of matched analysis results
		try
		{
			while (selectNextCorrelation(results))
				;
		}
		catch (Exception e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		// Remove bad results from the dataset.
		ArrayList<CorrelationResult> newResults = new ArrayList<CorrelationResult>(results.size());
		for (int i = 0; i < results.size(); i++)
		{
			CorrelationResult r = results.get(i);
			// If the PC-PALM Analysis has been done on too few molecules then the g(r) curve will be bad 
			if (r.n < 10)
			{
				header();
				log("Excluding dataset ID %d - Too few unique points (%f)", r.id, r.n);
				continue;
			}
			// If the PC-PALM Analysis has a g(r) curve all below 1 then it is not valid
			int offset = r.spatialDomain ? 0 : 1;
			double max = 0;
			for (int j = offset; j < r.gr[1].length; j++)
			{
				if (max < r.gr[1][j])
					max = r.gr[1][j];
			}
			if (max < 1)
			{
				header();
				log("Excluding dataset ID %d - g(r) curve is always below 1 (max = %f)", r.id, max);
				continue;
			}
			newResults.add(r);
		}

		results = newResults;
		return !results.isEmpty();
	}

	private boolean selectNextCorrelation(ArrayList<CorrelationResult> results)
	{
		ArrayList<String> titles = buildTitlesList(results);

		// Show a dialog allowing the user to select an input image
		if (titles.isEmpty())
			return false;

		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);
		gd.addMessage("Select the next correlation curve\nFrequency domain curves are identified with *");
		int n = (results.size() + 1);

		// If in macro mode then we must just use the String input field to allow the macro
		// IJ to return the field values from the macro arguments. Using a Choice input
		// will always return a field value.

		if (IJ.isMacro())
			gd.addStringField("R_" + n, titles.get(0));
		else
			gd.addChoice("R_" + n, titles.toArray(new String[titles.size()]), "");

		gd.addMessage("Cancel to finish");
		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		String title;
		if (IJ.isMacro())
			title = gd.getNextString();
		else
			title = gd.getNextChoice();

		// Check the correlation exists. If not then exit. This is mainly relevant for Macro mode since
		// the loop will continue otherwise since the titles list is not empty.
		String[] fields = title.split("\\*?:");
		try
		{
			int id = Integer.parseInt(fields[0]);
			for (CorrelationResult r : PCPALMAnalysis.results)
			{
				if (r.id == id)
				{
					results.add(r);
					return true;
				}
			}
		}
		catch (NumberFormatException e)
		{
		}
		return false;
	}

	private ArrayList<String> buildTitlesList(ArrayList<CorrelationResult> results)
	{
		// Make all subsequent results match the same nmPerPixel limit
		double nmPerPixel = 0;
		boolean spatialDomain = false;
		boolean filter = false;
		if (!results.isEmpty())
		{
			filter = true;
			nmPerPixel = results.get(0).nmPerPixel;
			spatialDomain = results.get(0).spatialDomain;
		}

		ArrayList<String> titles = new ArrayList<String>();
		for (CorrelationResult r : PCPALMAnalysis.results)
		{
			if (alreadySelected(results, r) ||
					(filter && (r.nmPerPixel != nmPerPixel || r.spatialDomain != spatialDomain)))
				continue;
			titles.add(String.format("%d%s: %s (%s nm/px)", r.id, (r.spatialDomain) ? "" : "*", r.source.getName(),
					Utils.rounded(r.nmPerPixel, 3)));
		}
		return titles;
	}

	private boolean alreadySelected(ArrayList<CorrelationResult> results, CorrelationResult r)
	{
		for (CorrelationResult r2 : results)
			if (r.id == r2.id)
				return true;
		return false;
	}

	/**
	 * Log a message to the IJ log window
	 * 
	 * @param format
	 * @param args
	 */
	private static void log(String format, Object... args)
	{
		IJ.log(String.format(format, args));
	}

	/**
	 * Perform the PC Analysis
	 * <p>
	 * Spatial domain results can just be combined to an average curve.
	 * <p>
	 * Frequency domain results can be fit using the g(r) model.
	 */
	private void analyse()
	{
		previous_gr = gr;
		previous_peakDensity = peakDensity;
		previous_spatialDomain = spatialDomain;

		String axisTitle;
		if (spatialDomain)
		{
			offset = 0;
			axisTitle = "molecules/um^2";
		}
		else
		{
			// Ignore the r=0 value by starting with an offset if necessary
			offset = (gr[0][0] == 0) ? 1 : 0;
			axisTitle = "g(r)";
		}
		String title = TITLE + " " + axisTitle;
		Plot plot = PCPALMAnalysis.plotCorrelation(gr, offset, title, axisTitle, spatialDomain, showErrorBars);

		if (spatialDomain)
		{
			saveCorrelationCurve(gr);
			log("Created correlation curve from the spatial domain (Plot title = " + title + ")");
			return;
		}

		// -------------
		// Model fitting for g(r) correlation curves
		// -------------
		log("Fitting g(r) correlation curve from the frequency domain");
		log("Average peak density = %s um^-2. Blinking estimate = %s", Utils.rounded(peakDensity, 4),
				Utils.rounded(blinkingRate, 4));

		createResultsTable();

		// Get the protein density in nm^2. 
		peakDensity /= 1e6;

		// Use the blinking rate estimate to estimate the density
		// (factors in the over-counting of the same molecules)
		double proteinDensity = peakDensity / blinkingRate;

		ArrayList<double[]> curves = new ArrayList<double[]>();

		// Fit the g(r) curve for r>0 to equation 2
		Color color = Color.red;
		String resultColour = "Red";
		double[] parameters = fitRandomModel(gr, estimatedPrecision, proteinDensity, resultColour);
		if (parameters != null)
		{
			log("  Plot %s: Over-counting estimate = %s", randomModel.getName(),
					Utils.rounded(peakDensity / parameters[1], 4));
			log("  Plot %s == %s", randomModel.getName(), resultColour.toString());
			plot.setColor(color);
			plot.addPoints(randomModel.getX(), randomModel.value(parameters), Plot.LINE);
			addNonFittedPoints(plot, gr, randomModel, parameters);
			Utils.display(title, plot);
			if (saveCorrelationCurve)
				curves.add(extractCurve(gr, randomModel, parameters));
		}

		// Fit the clustered models if the random model fails or if chosen as an option
		if (!valid1 || fitClusteredModels)
		{
			// Fit the g(r) curve for r>0 to equation 3
			color = Color.blue;
			resultColour = "Blue";
			parameters = fitClusteredModel(gr, estimatedPrecision, proteinDensity, resultColour);

			if (parameters != null)
			{
				log("  Plot %s: Over-counting estimate = %s", clusteredModel.getName(),
						Utils.rounded(peakDensity / parameters[1], 4));
				log("  Plot %s == %s, ", clusteredModel.getName(), resultColour.toString());
				plot.setColor(color);
				plot.addPoints(clusteredModel.getX(), clusteredModel.value(parameters), Plot.LINE);
				addNonFittedPoints(plot, gr, clusteredModel, parameters);
				Utils.display(title, plot);
				if (saveCorrelationCurve)
					curves.add(extractCurve(gr, clusteredModel, parameters));
			}

			// Fit to an emulsion model for a distribution confined to circles
			color = Color.magenta;
			resultColour = "Magenta";
			parameters = fitEmulsionModel(gr, estimatedPrecision, proteinDensity, resultColour);

			if (parameters != null)
			{
				log("  Plot %s: Over-counting estimate = %s", emulsionModel.getName(),
						Utils.rounded(peakDensity / parameters[1], 4));
				log("  Plot %s == %s", emulsionModel.getName(), resultColour.toString());
				plot.setColor(color);
				plot.addPoints(emulsionModel.getX(), emulsionModel.value(parameters), Plot.LINE);
				addNonFittedPoints(plot, gr, emulsionModel, parameters);
				Utils.display(title, plot);
				if (saveCorrelationCurve)
					curves.add(extractCurve(gr, emulsionModel, parameters));
			}
		}

		saveCorrelationCurve(gr, curves.toArray(new double[0][0]));
	}

	private void addNonFittedPoints(Plot plot, double[][] gr, BaseModelFunction model, double[] parameters)
	{
		double[] x = new double[gr[0].length];
		double[] y = new double[x.length];
		int j = 0;
		for (int i = offset; i < gr[0].length; i++)
		{
			// Output points that were not fitted
			if (gr[0][i] < randomModel.getX()[0])
			{
				x[j] = gr[0][i];
				y[j] = model.evaluate(gr[0][i], parameters);
				j++;
			}
		}
		x = Arrays.copyOf(x, j);
		y = Arrays.copyOf(y, j);
		plot.addPoints(x, y, Plot.CIRCLE);
	}

	private double[] extractCurve(double[][] gr, BaseModelFunction model, double[] parameters)
	{
		double[] y = new double[gr[0].length - offset];
		for (int i = offset; i < gr[0].length; i++)
		{
			y[i - offset] = model.evaluate(gr[0][i], parameters);
		}
		return y;
	}

	private double[][] combineCurves(ArrayList<CorrelationResult> results, int maxSize)
	{
		double[][] gr = new double[3][maxSize];
		Statistics[] gr_ = new Statistics[maxSize];
		for (int i = 0; i < maxSize; i++)
			gr_[i] = new Statistics();

		for (CorrelationResult r : results)
		{
			for (int i = 0; i < r.gr[0].length; i++)
			{
				gr[0][i] = r.gr[0][i]; // All scales should be the same so over-write is OK

				// Note that sometimes the analysis generates values that are very bad (e.g. if too 
				// few points were analysed). Perhaps we should exclude outliers for each distance interval.

				// NaN values can be generated so ignore them
				if (!Double.isNaN(r.gr[1][i]))
					gr_[i].add(r.gr[1][i]);
			}
		}
		for (int i = 0; i < maxSize; i++)
		{
			gr[1][i] = gr_[i].getMean();
			gr[2][i] = gr_[i].getStandardError();
		}
		return gr;
	}

	private void saveCorrelationCurve(double[][] gr, double[]... curves)
	{
		if (!saveCorrelationCurve)
			return;
		outputFilename = Utils.getFilename("Output_Correlation_File", outputFilename);
		if (outputFilename != null)
		{
			outputFilename = Utils.replaceExtension(outputFilename, "xls");

			BufferedWriter output = null;
			try
			{
				output = new BufferedWriter(new FileWriter(outputFilename));
				writeHeader(output, HEADER_PEAK_DENSITY, Double.toString(previous_peakDensity));
				writeHeader(output, HEADER_SPATIAL_DOMAIN, Boolean.toString(previous_spatialDomain));
				output.write("#r\tg(r)\tS.E.");
				for (int j = 0; j < curves.length; j++)
					output.write(String.format("\tModel %d", j + 1));
				output.newLine();
				// Ignore the r=0 value by starting with an offset if necessary
				for (int i = offset; i < gr[0].length; i++)
				{
					output.write(String.format("%f\t%f\t%f", gr[0][i], gr[1][i], gr[2][i]));
					for (int j = 0; j < curves.length; j++)
						output.write(String.format("\t%f", curves[j][i - offset]));
					output.newLine();
				}
			}
			catch (Exception e)
			{
				// Q. Add better handling of errors?
				e.printStackTrace();
				IJ.log("Failed to save correlation curve to file: " + outputFilename);
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

	private void writeHeader(BufferedWriter output, String header, String value) throws IOException
	{
		output.write("#");
		output.write(header);
		output.write(" = ");
		output.write(value);
		output.newLine();
	}

	/**
	 * Load a correlation curve from file. Will set the global gr, peakDensity and spatialDomain variables. If the data
	 * fails to be loaded then the method will return false.
	 * 
	 * @return True if loaded
	 */
	private boolean loadCorrelationCurve()
	{
		inputFilename = Utils.getFilename("Input_Correlation_File", inputFilename);
		if (inputFilename == null)
			return false;

		// Set the analysis variables
		boolean spatialDomainSet = false;
		boolean peakDensitySet = false;

		BufferedReader input = null;
		try
		{
			FileInputStream fis = new FileInputStream(inputFilename);
			input = new BufferedReader(new UnicodeReader(fis, null));

			String line;
			int count = 0;

			Pattern pattern = Pattern.compile("#([^=]+) = ([^ ]+)");

			// Read the header
			while ((line = input.readLine()) != null)
			{
				count++;

				if (line.length() == 0)
					continue;
				if (line.charAt(0) != '#')
				{
					// This is the first record
					break;
				}

				// This is a header line. Extract the key-value pair
				Matcher match = pattern.matcher(line);
				if (match.find())
				{
					if (match.group(1).equals(HEADER_SPATIAL_DOMAIN))
					{
						// Do not use Boolean.parseBoolean because this will not fail if the field is 
						// neither true/false - it only return true for a match to true
						spatialDomainSet = true;
						if (match.group(2).equalsIgnoreCase("true"))
							spatialDomain = true;
						else if (match.group(2).equalsIgnoreCase("false"))
							spatialDomain = false;
						else
							// We want to know if the field is not true/false
							spatialDomainSet = false;
					}
					else if (match.group(1).equals(HEADER_PEAK_DENSITY))
					{
						try
						{
							peakDensity = Double.parseDouble(match.group(2));
							peakDensitySet = true;
						}
						catch (NumberFormatException e)
						{
							// Ignore this.
						}
					}
				}
			}

			if (!peakDensitySet)
			{
				IJ.error(TITLE, "No valid " + HEADER_PEAK_DENSITY + " record in file " + inputFilename);
				return false;
			}
			if (!spatialDomainSet)
			{
				IJ.error(TITLE, "No valid " + HEADER_SPATIAL_DOMAIN + " record in file " + inputFilename);
				return false;
			}

			// Read the data: gr[0][i], gr[1][i], gr[2][i]
			ArrayList<double[]> data = new ArrayList<double[]>();
			while (line != null)
			{
				if (line.length() == 0)
					continue;
				if (line.charAt(0) == '#')
					continue;

				// Extract the first 3 fields
				Scanner scanner = new Scanner(line);
				scanner.useDelimiter("[\t ,]+");

				double r, g;
				try
				{
					r = scanner.nextDouble();
					g = scanner.nextDouble();
				}
				catch (InputMismatchException e)
				{
					IJ.error(TITLE, "Incorrect fields on line " + count);
					scanner.close();
					return false;
				}
				catch (NoSuchElementException e)
				{
					IJ.error(TITLE, "Incorrect fields on line " + count);
					scanner.close();
					return false;
				}
				// Allow the file to be missing the curve error. This is only used for plotting anyway.
				double error = 0;
				try
				{
					error = scanner.nextDouble();
				}
				catch (InputMismatchException e)
				{
				}
				catch (NoSuchElementException e)
				{
				}
				scanner.close();
				data.add(new double[] { r, g, error });

				// Read the next line
				line = input.readLine();
				count++;
			}

			if (data.isEmpty())
			{
				IJ.error(TITLE, "No data in file " + inputFilename);
				return false;
			}

			gr = new double[3][data.size()];
			for (int i = 0; i < data.size(); i++)
			{
				final double[] d = data.get(i);
				gr[0][i] = d[0];
				gr[1][i] = d[1];
				gr[2][i] = d[2];
			}

		}
		catch (IOException e)
		{
			IJ.error(TITLE, "Unable to read from file " + inputFilename);
			return false;
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

		return true;
	}

	/**
	 * Fits the correlation curve with r>0 to the random model using the estimated density and precision. Parameters
	 * must be fit within a tolerance of the starting values.
	 * 
	 * @param gr
	 * @param sigmaS
	 *            The estimated precision
	 * @param proteinDensity
	 *            The estimate protein density
	 * @return The fitted parameters [precision, density]
	 */
	private double[] fitRandomModel(double[][] gr, double sigmaS, double proteinDensity, String resultColour)
	{
		final RandomModelFunction myFunction = new RandomModelFunction();
		randomModel = myFunction;

		log("Fitting %s: Estimated precision = %f nm, estimated protein density = %g um^-2", randomModel.getName(),
				sigmaS, proteinDensity * 1e6);

		randomModel.setLogging(true);

		for (int i = offset; i < gr[0].length; i++)
		{
			// Only fit the curve above the estimated resolution (points below it will be subject to error)
			if (gr[0][i] > sigmaS * fitAboveEstimatedPrecision)
				randomModel.addPoint(gr[0][i], gr[1][i]);
		}
		LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer();

		PointVectorValuePair optimum;
		try
		{
			optimum = optimizer.optimize(new MaxIter(3000), new MaxEval(Integer.MAX_VALUE), new ModelFunctionJacobian(
					new MultivariateMatrixFunction()
					{
						public double[][] value(double[] point) throws IllegalArgumentException
						{
							return myFunction.jacobian(point);
						}
					}), new ModelFunction(myFunction), new Target(myFunction.getY()),
					new Weight(myFunction.getWeights()), new InitialGuess(new double[] { sigmaS, proteinDensity }));
		}
		catch (TooManyIterationsException e)
		{
			log("Failed to fit %s: Too many iterations (%d)", randomModel.getName(), optimizer.getIterations());
			return null;
		}
		catch (ConvergenceException e)
		{
			log("Failed to fit %s: %s", randomModel.getName(), e.getMessage());
			return null;
		}

		randomModel.setLogging(false);

		double[] parameters = optimum.getPoint();
		// Ensure the width is positive
		parameters[0] = Math.abs(parameters[0]);

		double ss = 0;
		double[] obs = randomModel.getY();
		double[] exp = optimum.getValue();
		for (int i = 0; i < obs.length; i++)
			ss += (obs[i] - exp[i]) * (obs[i] - exp[i]);
		ic1 = Maths.getInformationCriterion(ss, randomModel.size(), parameters.length);

		final double fitSigmaS = parameters[0];
		final double fitProteinDensity = parameters[1];

		// Check the fitted parameters are within tolerance of the initial estimates
		double e1 = parameterDrift(sigmaS, fitSigmaS);
		double e2 = parameterDrift(proteinDensity, fitProteinDensity);

		log("  %s fit: SS = %f. cAIC = %f. %d evaluations", randomModel.getName(), ss, ic1, optimizer.getEvaluations());
		log("  %s parameters:", randomModel.getName());
		log("    Average precision = %s nm (%s%%)", Utils.rounded(fitSigmaS, 4), Utils.rounded(e1, 4));
		log("    Average protein density = %s um^-2 (%s%%)", Utils.rounded(fitProteinDensity * 1e6, 4),
				Utils.rounded(e2, 4));

		valid1 = true;
		if (fittingTolerance > 0 && (Math.abs(e1) > fittingTolerance || Math.abs(e2) > fittingTolerance))
		{
			log("  Failed to fit %s within tolerance (%s%%): Average precision = %f nm (%s%%), average protein density = %g um^-2 (%s%%)",
					randomModel.getName(), Utils.rounded(fittingTolerance, 4), fitSigmaS, Utils.rounded(e1, 4),
					fitProteinDensity * 1e6, Utils.rounded(e2, 4));
			valid1 = false;
		}

		if (valid1)
		{
			// ---------
			// TODO - My data does not comply with this criteria. 
			// This could be due to the PC-PALM Molecule code limiting the nmPerPixel to fit the images in memory 
			// thus removing correlations at small r.
			// It could also be due to the nature of the random simulations being 3D not 2D membranes 
			// as per the PC-PALM paper. 
			// ---------
			// Evaluate g(r)protein where:
			// g(r)peaks = g(r)protein + g(r)stoch
			// g(r)peaks ~ 1           + g(r)stoch
			// Verify g(r)protein should be <1.5 for all r>0
			double[] gr_stoch = randomModel.value(parameters);
			double[] gr_peaks = randomModel.getY();
			double[] gr_ = randomModel.getX();

			//SummaryStatistics stats = new SummaryStatistics();
			for (int i = 0; i < gr_peaks.length; i++)
			{
				// Only evaluate above the fitted average precision 
				if (gr_[i] < fitSigmaS)
					continue;

				// Note the RandomModelFunction evaluates g(r)stoch + 1;
				double gr_protein_i = gr_peaks[i] - (gr_stoch[i] - 1);

				if (gr_protein_i > gr_protein_threshold)
				{
					// Failed fit
					log("  Failed to fit %s: g(r)protein %s > %s @ r=%s", randomModel.getName(),
							Utils.rounded(gr_protein_i, 4), Utils.rounded(gr_protein_threshold, 4),
							Utils.rounded(gr_[i], 4));
					valid1 = false;
				}
				//stats.addValue(gr_i);
				//System.out.printf("g(r)protein @ %f = %f\n", gr[0][i], gr_protein_i);
			}
		}

		addResult(randomModel.getName(), resultColour, valid1, fitSigmaS, fitProteinDensity, 0, 0, 0, 0, ic1);

		return parameters;
	}

	private double parameterDrift(double start, double end)
	{
		if (end < start)
		{
			return -100 * (start - end) / end;
		}
		return 100 * (end - start) / start;
	}

	/**
	 * Fits the correlation curve with r>0 to the clustered model using the estimated density and precision. Parameters
	 * must be fit within a tolerance of the starting values.
	 * 
	 * @param gr
	 * @param sigmaS
	 *            The estimated precision
	 * @param proteinDensity
	 *            The estimated protein density
	 * @return The fitted parameters [precision, density, clusterRadius, clusterDensity]
	 */
	private double[] fitClusteredModel(double[][] gr, double sigmaS, double proteinDensity, String resultColour)
	{
		final ClusteredModelFunctionGradient myFunction = new ClusteredModelFunctionGradient();
		clusteredModel = myFunction;
		log("Fitting %s: Estimated precision = %f nm, estimated protein density = %g um^-2", clusteredModel.getName(),
				sigmaS, proteinDensity * 1e6);

		clusteredModel.setLogging(true);
		for (int i = offset; i < gr[0].length; i++)
		{
			// Only fit the curve above the estimated resolution (points below it will be subject to error)
			if (gr[0][i] > sigmaS * fitAboveEstimatedPrecision)
				clusteredModel.addPoint(gr[0][i], gr[1][i]);
		}

		double[] parameters;
		// The model is: sigma, density, range, amplitude
		double[] initialSolution = new double[] { sigmaS, proteinDensity, sigmaS * 5, 1 };
		int evaluations = 0;

		// Constrain the fitting to be close to the estimated precision (sigmaS) and protein density.
		// LVM fitting does not support constrained fitting so use a bounded optimiser.
		SumOfSquaresModelFunction clusteredModelMulti = new SumOfSquaresModelFunction(clusteredModel);
		double[] x = clusteredModelMulti.x;

		// Put some bounds around the initial guess. Use the fitting tolerance (in %) if provided.
		double limit = (fittingTolerance > 0) ? 1 + fittingTolerance / 100 : 2;
		double[] lB = new double[] { initialSolution[0] / limit, initialSolution[1] / limit, 0, 0 };
		// The amplitude and range should not extend beyond the limits of the g(r) curve.
		double[] uB = new double[] { initialSolution[0] * limit, initialSolution[1] * limit, Maths.max(x),
				Maths.max(gr[1]) };
		log("Fitting %s using a bounded search: %s < precision < %s & %s < density < %s", clusteredModel.getName(),
				Utils.rounded(lB[0], 4), Utils.rounded(uB[0], 4), Utils.rounded(lB[1] * 1e6, 4),
				Utils.rounded(uB[1] * 1e6, 4));

		PointValuePair constrainedSolution = runBoundedOptimiser(gr, initialSolution, lB, uB, clusteredModelMulti);

		if (constrainedSolution == null)
			return null;

		parameters = constrainedSolution.getPointRef();
		evaluations = boundedEvaluations;

		// Refit using a LVM  
		if (useLSE)
		{
			log("Re-fitting %s using a gradient optimisation", clusteredModel.getName());
			LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer();
			PointVectorValuePair lvmSolution;
			try
			{
				lvmSolution = optimizer.optimize(new MaxIter(3000), new MaxEval(Integer.MAX_VALUE),
						new ModelFunctionJacobian(new MultivariateMatrixFunction()
						{
							public double[][] value(double[] point) throws IllegalArgumentException
							{
								return myFunction.jacobian(point);
							}
						}), new ModelFunction(myFunction), new Target(myFunction.getY()),
						new Weight(myFunction.getWeights()), new InitialGuess(parameters));
				evaluations += optimizer.getEvaluations();

				double ss = 0;
				double[] obs = clusteredModel.getY();
				double[] exp = lvmSolution.getValue();
				for (int i = 0; i < obs.length; i++)
					ss += (obs[i] - exp[i]) * (obs[i] - exp[i]);
				if (ss < constrainedSolution.getValue())
				{
					log("Re-fitting %s improved the SS from %s to %s (-%s%%)", clusteredModel.getName(), Utils.rounded(
							constrainedSolution.getValue(), 4), Utils.rounded(ss, 4), Utils.rounded(100 *
							(constrainedSolution.getValue() - ss) / constrainedSolution.getValue(), 4));
					parameters = lvmSolution.getPoint();
				}
			}
			catch (TooManyIterationsException e)
			{
				log("Failed to re-fit %s: Too many iterations (%d)", clusteredModel.getName(),
						optimizer.getIterations());
			}
			catch (ConvergenceException e)
			{
				log("Failed to re-fit %s: %s", clusteredModel.getName(), e.getMessage());
			}
		}

		clusteredModel.setLogging(false);

		// Ensure the width is positive
		parameters[0] = Math.abs(parameters[0]);
		//parameters[2] = Math.abs(parameters[2]);

		double ss = 0;
		double[] obs = clusteredModel.getY();
		double[] exp = clusteredModel.value(parameters);
		for (int i = 0; i < obs.length; i++)
			ss += (obs[i] - exp[i]) * (obs[i] - exp[i]);
		ic2 = Maths.getInformationCriterion(ss, clusteredModel.size(), parameters.length);

		final double fitSigmaS = parameters[0];
		final double fitProteinDensity = parameters[1];
		final double domainRadius = parameters[2]; //The radius of the cluster domain
		final double domainDensity = parameters[3]; //The density of the cluster domain

		// This is from the PC-PALM paper. However that paper fits the g(r)protein exponential convolved in 2D
		// with the g(r)PSF. In this method we have just fit the exponential
		final double nCluster = 2 * domainDensity * Math.PI * domainRadius * domainRadius * fitProteinDensity;

		double e1 = parameterDrift(sigmaS, fitSigmaS);
		double e2 = parameterDrift(proteinDensity, fitProteinDensity);

		log("  %s fit: SS = %f. cAIC = %f. %d evaluations", clusteredModel.getName(), ss, ic2, evaluations);
		log("  %s parameters:", clusteredModel.getName());
		log("    Average precision = %s nm (%s%%)", Utils.rounded(fitSigmaS, 4), Utils.rounded(e1, 4));
		log("    Average protein density = %s um^-2 (%s%%)", Utils.rounded(fitProteinDensity * 1e6, 4),
				Utils.rounded(e2, 4));
		log("    Domain radius = %s nm", Utils.rounded(domainRadius, 4));
		log("    Domain density = %s", Utils.rounded(domainDensity, 4));
		log("    nCluster = %s", Utils.rounded(nCluster, 4));

		// Check the fitted parameters are within tolerance of the initial estimates
		valid2 = true;
		if (fittingTolerance > 0 && (Math.abs(e1) > fittingTolerance || Math.abs(e2) > fittingTolerance))
		{
			log("  Failed to fit %s within tolerance (%s%%): Average precision = %f nm (%s%%), average protein density = %g um^-2 (%s%%)",
					clusteredModel.getName(), Utils.rounded(fittingTolerance, 4), fitSigmaS, Utils.rounded(e1, 4),
					fitProteinDensity * 1e6, Utils.rounded(e2, 4));
			valid2 = false;
		}

		// Check extra parameters. Domain radius should be higher than the precision. Density should be positive
		if (domainRadius < fitSigmaS)
		{
			log("  Failed to fit %s: Domain radius is smaller than the average precision (%s < %s)",
					clusteredModel.getName(), Utils.rounded(domainRadius, 4), Utils.rounded(fitSigmaS, 4));
			valid2 = false;
		}
		if (domainDensity < 0)
		{
			log("  Failed to fit %s: Domain density is negative (%s)", clusteredModel.getName(),
					Utils.rounded(domainDensity, 4));
			valid2 = false;
		}

		if (ic2 > ic1)
		{
			log("  Failed to fit %s - Information Criterion has increased %s%%", clusteredModel.getName(),
					Utils.rounded((100 * (ic2 - ic1) / ic1), 4));
			valid2 = false;
		}

		addResult(clusteredModel.getName(), resultColour, valid2, fitSigmaS, fitProteinDensity, domainRadius,
				domainDensity, nCluster, 0, ic2);

		return parameters;
	}

	private PointValuePair runBoundedOptimiser(double[][] gr, double[] initialSolution, double[] lB, double[] uB,
			SumOfSquaresModelFunction function)
	{
		// Create the functions to optimise
		ObjectiveFunction objective = new ObjectiveFunction(new SumOfSquaresMultivariateFunction(function));
		ObjectiveFunctionGradient gradient = new ObjectiveFunctionGradient(new SumOfSquaresMultivariateVectorFunction(
				function));

		final boolean debug = false;

		// Try a BFGS optimiser since this will produce a deterministic solution and can respect bounds.
		PointValuePair optimum = null;
		boundedEvaluations = 0;
		final MaxEval maxEvaluations = new MaxEval(2000);
		MultivariateOptimizer opt = null;
		for (int iteration = 0; iteration <= fitRestarts; iteration++)
		{
			try
			{
				opt = new BFGSOptimizer();
				final double relativeThreshold = 1e-6;
				final double absoluteThreshold = 1e-10;

				// Configure maximum step length for each dimension using the bounds
				double[] stepLength = new double[lB.length];
				for (int i = 0; i < stepLength.length; i++)
					stepLength[i] = (uB[i] - lB[i]) * 0.3333333;

				// The GoalType is always minimise so no need to pass this in
				optimum = opt.optimize(maxEvaluations, gradient, objective, 
						new InitialGuess((optimum==null) ? initialSolution : optimum.getPointRef()),
						new SimpleBounds(lB, uB), new GradientChecker(relativeThreshold,
								absoluteThreshold), new PositionChecker(relativeThreshold,
								absoluteThreshold), new BFGSOptimizer.StepLength(stepLength));
				if (debug)
					System.out.printf("BFGS Iter %d = %g (%d)\n", iteration, optimum.getValue(), opt.getEvaluations());
			}
			catch (TooManyEvaluationsException e)
			{
				break; // No need to restart
			}
			catch (RuntimeException e)
			{
				break; // No need to restart
			}
			finally
			{
				boundedEvaluations += opt.getEvaluations();
			}
		}

		// Try a CMAES optimiser which is non-deterministic. To overcome this we perform restarts.

		// CMAESOptimiser based on Matlab code:
		// https://www.lri.fr/~hansen/cmaes.m
		// Take the defaults from the Matlab documentation
		double stopFitness = 0; //Double.NEGATIVE_INFINITY;
		boolean isActiveCMA = true;
		int diagonalOnly = 0;
		int checkFeasableCount = 1;
		RandomGenerator random = new Well44497b(); //Well19937c();
		boolean generateStatistics = false;
		ConvergenceChecker<PointValuePair> checker = new SimpleValueChecker(1e-6, 1e-10);
		// The sigma determines the search range for the variables. It should be 1/3 of the initial search region.
		double[] range = new double[lB.length];
		for (int i = 0; i < lB.length; i++)
			range[i] = (uB[i] - lB[i]) / 3;
		OptimizationData sigma = new CMAESOptimizer.Sigma(range);
		OptimizationData popSize = new CMAESOptimizer.PopulationSize((int) (4 + Math.floor(3 * Math
				.log(initialSolution.length))));
		SimpleBounds bounds = new SimpleBounds(lB, uB);

		// Restart the optimiser several times and take the best answer.
		for (int iteration = 0; iteration <= fitRestarts; iteration++)
		{
			try
			{
				// Start from the initial solution
				opt = new CMAESOptimizer(maxEvaluations.getMaxEval(), stopFitness, isActiveCMA, diagonalOnly,
						checkFeasableCount, random, generateStatistics, checker);
				PointValuePair constrainedSolution = opt.optimize(new InitialGuess(initialSolution), objective,
						GoalType.MINIMIZE, bounds, sigma, popSize, maxEvaluations);
				if (debug)
					System.out.printf("CMAES Iter %d initial = %g (%d)\n", iteration, constrainedSolution.getValue(),
							opt.getEvaluations());
				boundedEvaluations += opt.getEvaluations();
				if (optimum == null || constrainedSolution.getValue() < optimum.getValue())
				{
					optimum = constrainedSolution;
				}
			}
			catch (TooManyEvaluationsException e)
			{
			}
			catch (TooManyIterationsException e)
			{
			}
			finally
			{
				boundedEvaluations += maxEvaluations.getMaxEval();
			}
			if (optimum == null)
				continue;
			try
			{
				// Also restart from the current optimum
				opt = new CMAESOptimizer(maxEvaluations.getMaxEval(), stopFitness, isActiveCMA, diagonalOnly,
						checkFeasableCount, random, generateStatistics, checker);
				PointValuePair constrainedSolution = opt.optimize(new InitialGuess(optimum.getPointRef()), objective,
						GoalType.MINIMIZE, bounds, sigma, popSize, maxEvaluations);
				if (debug)
					System.out.printf("CMAES Iter %d restart = %g (%d)\n", iteration, constrainedSolution.getValue(),
							opt.getEvaluations());
				if (constrainedSolution.getValue() < optimum.getValue())
				{
					optimum = constrainedSolution;
				}
			}
			catch (TooManyEvaluationsException e)
			{
			}
			catch (TooManyIterationsException e)
			{
			}
			finally
			{
				boundedEvaluations += maxEvaluations.getMaxEval();
			}
		}
		return optimum;
	}

	/**
	 * Fits the correlation curve with r>0 to the clustered model using the estimated density and precision. Parameters
	 * must be fit within a tolerance of the starting values.
	 * 
	 * @param gr
	 * @param sigmaS
	 *            The estimated precision
	 * @param proteinDensity
	 *            The estimated protein density
	 * @return The fitted parameters [precision, density, clusterRadius, clusterDensity]
	 */
	private double[] fitEmulsionModel(double[][] gr, double sigmaS, double proteinDensity, String resultColour)
	{
		final EmulsionModelFunctionGradient myFunction = new EmulsionModelFunctionGradient();
		emulsionModel = myFunction;
		log("Fitting %s: Estimated precision = %f nm, estimated protein density = %g um^-2", emulsionModel.getName(),
				sigmaS, proteinDensity * 1e6);

		emulsionModel.setLogging(true);
		for (int i = offset; i < gr[0].length; i++)
		{
			// Only fit the curve above the estimated resolution (points below it will be subject to error)
			if (gr[0][i] > sigmaS * fitAboveEstimatedPrecision)
				emulsionModel.addPoint(gr[0][i], gr[1][i]);
		}

		double[] parameters;
		// The model is: sigma, density, range, amplitude, alpha
		double[] initialSolution = new double[] { sigmaS, proteinDensity, sigmaS * 5, 1, sigmaS * 5 };
		int evaluations = 0;

		// Constrain the fitting to be close to the estimated precision (sigmaS) and protein density.
		// LVM fitting does not support constrained fitting so use a bounded optimiser.
		SumOfSquaresModelFunction emulsionModelMulti = new SumOfSquaresModelFunction(emulsionModel);
		double[] x = emulsionModelMulti.x;
		double[] y = emulsionModelMulti.y;

		// Range should be equal to the first time the g(r) curve crosses 1
		for (int i = 0; i < x.length; i++)
			if (y[i] < 1)
			{
				initialSolution[4] = initialSolution[2] = (i > 0) ? (x[i - 1] + x[i]) * 0.5 : x[i];
				break;
			}

		// Put some bounds around the initial guess. Use the fitting tolerance (in %) if provided.
		double limit = (fittingTolerance > 0) ? 1 + fittingTolerance / 100 : 2;
		double[] lB = new double[] { initialSolution[0] / limit, initialSolution[1] / limit, 0, 0, 0 };
		// The amplitude and range should not extend beyond the limits of the g(r) curve.
		// TODO - Find out the expected range for the alpha parameter.  
		double[] uB = new double[] { initialSolution[0] * limit, initialSolution[1] * limit, Maths.max(x),
				Maths.max(gr[1]), Maths.max(x) * 2 };
		log("Fitting %s using a bounded search: %s < precision < %s & %s < density < %s", emulsionModel.getName(),
				Utils.rounded(lB[0], 4), Utils.rounded(uB[0], 4), Utils.rounded(lB[1] * 1e6, 4),
				Utils.rounded(uB[1] * 1e6, 4));

		PointValuePair constrainedSolution = runBoundedOptimiser(gr, initialSolution, lB, uB, emulsionModelMulti);

		if (constrainedSolution == null)
			return null;

		parameters = constrainedSolution.getPointRef();
		evaluations = boundedEvaluations;

		// Refit using a LVM  
		if (useLSE)
		{
			log("Re-fitting %s using a gradient optimisation", emulsionModel.getName());
			LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer();
			PointVectorValuePair lvmSolution;
			try
			{
				lvmSolution = optimizer.optimize(new MaxIter(3000), new MaxEval(Integer.MAX_VALUE),
						new ModelFunctionJacobian(new MultivariateMatrixFunction()
						{
							public double[][] value(double[] point) throws IllegalArgumentException
							{
								return myFunction.jacobian(point);
							}
						}), new ModelFunction(myFunction), new Target(myFunction.getY()),
						new Weight(myFunction.getWeights()), new InitialGuess(parameters));
				evaluations += optimizer.getEvaluations();

				double ss = 0;
				double[] obs = emulsionModel.getY();
				double[] exp = lvmSolution.getValue();
				for (int i = 0; i < obs.length; i++)
					ss += (obs[i] - exp[i]) * (obs[i] - exp[i]);
				if (ss < constrainedSolution.getValue())
				{
					log("Re-fitting %s improved the SS from %s to %s (-%s%%)", emulsionModel.getName(), Utils.rounded(
							constrainedSolution.getValue(), 4), Utils.rounded(ss, 4), Utils.rounded(100 *
							(constrainedSolution.getValue() - ss) / constrainedSolution.getValue(), 4));
					parameters = lvmSolution.getPoint();
				}
			}
			catch (TooManyIterationsException e)
			{
				log("Failed to re-fit %s: Too many iterations (%d)", emulsionModel.getName(), optimizer.getIterations());
			}
			catch (ConvergenceException e)
			{
				log("Failed to re-fit %s: %s", emulsionModel.getName(), e.getMessage());
			}
		}

		emulsionModel.setLogging(false);

		// Ensure the width is positive
		parameters[0] = Math.abs(parameters[0]);
		//parameters[2] = Math.abs(parameters[2]);

		double ss = 0;
		double[] obs = emulsionModel.getY();
		double[] exp = emulsionModel.value(parameters);
		for (int i = 0; i < obs.length; i++)
			ss += (obs[i] - exp[i]) * (obs[i] - exp[i]);
		ic3 = Maths.getInformationCriterion(ss, emulsionModel.size(), parameters.length);

		final double fitSigmaS = parameters[0];
		final double fitProteinDensity = parameters[1];
		final double domainRadius = parameters[2]; //The radius of the cluster domain
		final double domainDensity = parameters[3]; //The density of the cluster domain
		final double coherence = parameters[4]; //The coherence length between circles

		// This is from the PC-PALM paper. It may not be correct for the emulsion model.
		final double nCluster = 2 * domainDensity * Math.PI * domainRadius * domainRadius * fitProteinDensity;

		double e1 = parameterDrift(sigmaS, fitSigmaS);
		double e2 = parameterDrift(proteinDensity, fitProteinDensity);

		log("  %s fit: SS = %f. cAIC = %f. %d evaluations", emulsionModel.getName(), ss, ic3, evaluations);
		log("  %s parameters:", emulsionModel.getName());
		log("    Average precision = %s nm (%s%%)", Utils.rounded(fitSigmaS, 4), Utils.rounded(e1, 4));
		log("    Average protein density = %s um^-2 (%s%%)", Utils.rounded(fitProteinDensity * 1e6, 4),
				Utils.rounded(e2, 4));
		log("    Domain radius = %s nm", Utils.rounded(domainRadius, 4));
		log("    Domain density = %s", Utils.rounded(domainDensity, 4));
		log("    Domain coherence = %s", Utils.rounded(coherence, 4));
		log("    nCluster = %s", Utils.rounded(nCluster, 4));

		// Check the fitted parameters are within tolerance of the initial estimates
		valid2 = true;
		if (fittingTolerance > 0 && (Math.abs(e1) > fittingTolerance || Math.abs(e2) > fittingTolerance))
		{
			log("  Failed to fit %s within tolerance (%s%%): Average precision = %f nm (%s%%), average protein density = %g um^-2 (%s%%)",
					emulsionModel.getName(), Utils.rounded(fittingTolerance, 4), fitSigmaS, Utils.rounded(e1, 4),
					fitProteinDensity * 1e6, Utils.rounded(e2, 4));
			valid2 = false;
		}

		// Check extra parameters. Domain radius should be higher than the precision. Density should be positive
		if (domainRadius < fitSigmaS)
		{
			log("  Failed to fit %s: Domain radius is smaller than the average precision (%s < %s)",
					emulsionModel.getName(), Utils.rounded(domainRadius, 4), Utils.rounded(fitSigmaS, 4));
			valid2 = false;
		}
		if (domainDensity < 0)
		{
			log("  Failed to fit %s: Domain density is negative (%s)", emulsionModel.getName(),
					Utils.rounded(domainDensity, 4));
			valid2 = false;
		}

		if (ic3 > ic1)
		{
			log("  Failed to fit %s - Information Criterion has increased %s%%", emulsionModel.getName(),
					Utils.rounded((100 * (ic3 - ic1) / ic1), 4));
			valid2 = false;
		}

		addResult(emulsionModel.getName(), resultColour, valid2, fitSigmaS, fitProteinDensity, domainRadius,
				domainDensity, nCluster, coherence, ic3);

		return parameters;
	}

	/**
	 * Abstract base model function class for common functionality.
	 */
	public abstract class BaseModelFunction extends LoggingOptimiserFunction
	{
		public BaseModelFunction(String name)
		{
			super(name);
		}

		/**
		 * Evaluate the correlation function
		 * 
		 * @param r
		 *            The correlation radius
		 * @param parameters
		 *            The parameters
		 * @return
		 */
		public abstract double evaluate(double r, final double[] parameters);

		/**
		 * Evaluate the jacobian of the correlation function for all data points (see
		 * {@link #addData(double[], double[])})
		 * 
		 * @param parameters
		 * @return The jacobian
		 */
		public abstract double[][] jacobian(double[] parameters);

		/**
		 * Get the value of the function for all data points corresponding to the last call to
		 * {@link #jacobian(double[])}
		 * 
		 * @return The corresponding value
		 */
		public abstract double[] getValue();
	}

	/**
	 * Allow optimisation using Apache Commons Math 3 Gradient Optimiser
	 * <p>
	 * g(r)peaks = g(r)stoch + 1
	 * <p>
	 * where
	 * <p>
	 * g(r)stoch = (1/4*pi*s^2*p) * exp(-r^2/4s^2)
	 * <p>
	 * s = average single molecule positional uncertainty (precision)
	 * <p>
	 * p = average protein density
	 */
	public class RandomModelFunction extends BaseModelFunction implements MultivariateVectorFunction
	{
		double[] lastValue = null;

		public RandomModelFunction()
		{
			super("Random Model");
		}

		@Override
		public double[] getValue()
		{
			return lastValue;
		}

		// Adapted from http://commons.apache.org/proper/commons-math/userguide/optimization.html
		// Use the deprecated API since the new one is not yet documented.

		public double[][] jacobian(double[] variables)
		{
			// Compute the gradients using calculus differentiation
			final double sigma = variables[0];
			final double density = variables[1];
			final double[][] jacobian = new double[x.size()][2];
			lastValue = new double[x.size()];

			for (int i = 0; i < jacobian.length; ++i)
			{
				final double r = this.x.get(i);

				final double a = 1.0 / (4 * Math.PI * density * sigma * sigma);
				final double b = -r * r / (4 * sigma * sigma);
				final double c = FastMath.exp(b);

				// value  = a * c
				lastValue[i] = a * c;

				// Differentiate with respect to sigma:
				// value' = a' * c + a * c'  [ Product rule ]
				// c = FastMath.exp(b)
				// c' = b' * FastMath.exp(b)     [ Chain rule ]
				// value' = a' * c + a * b' * c
				jacobian[i][0] = (-2 * a / sigma) * c + a * (-2 * b / sigma) * c;

				// Differentiate with respect to density:
				// c' = 0 since density does not feature in c
				// => value' = a' * c
				jacobian[i][1] = (-a / density) * c;
			}

			//// Check numerically ...
			//double[][] jacobian2 = jacobian2(variables);
			//for (int i = 0; i < jacobian.length; i++)
			//{
			//	System.out.printf("dSigma = %g : %g = %g. dDensity = %g : %g = %g\n", jacobian[i][0], jacobian2[i][0],
			//			DoubleEquality.relativeError(jacobian[i][0], jacobian2[i][0]), jacobian[i][1], jacobian2[i][1],
			//			DoubleEquality.relativeError(jacobian[i][1], jacobian2[i][1]));
			//}

			return jacobian;
		}

		@SuppressWarnings("unused")
		private double[][] jacobian2(double[] variables)
		{
			// Compute the gradients using numerical differentiation
			final double sigma = variables[0];
			final double density = variables[1];
			final double[][] jacobian = new double[x.size()][2];
			lastValue = new double[x.size()];

			final double delta = 0.001;
			final double[][] d = new double[variables.length][variables.length];
			for (int i = 0; i < variables.length; i++)
				d[i][i] = delta * Math.abs(variables[i]); // Should the delta be changed for each parameter ?
			for (int i = 0; i < jacobian.length; ++i)
			{
				final double r = this.x.get(i);
				final double value = lastValue[i] = evaluate(r, sigma, density);
				for (int j = 0; j < variables.length; j++)
				{
					final double value2 = evaluate(r, sigma + d[0][j], density + d[1][j]);
					jacobian[i][j] = (value2 - value) / d[j][j];
				}
			}
			return jacobian;
		}

		/**
		 * Evaluate the correlation function
		 * 
		 * @param r
		 *            The correlation radius
		 * @param sigma
		 *            Average precision
		 * @param density
		 *            Average protein density
		 * @return
		 */
		public double evaluate(double r, final double sigma, final double density)
		{
			return (1.0 / (4 * Math.PI * density * sigma * sigma)) * FastMath.exp(-r * r / (4 * sigma * sigma)) + 1;
		}

		/**
		 * Evaluate the correlation function
		 * 
		 * @param r
		 *            The correlation radius
		 * @param parameters
		 *            The parameters
		 * @return
		 */
		public double evaluate(double r, final double[] parameters)
		{
			return evaluate(r, parameters[0], parameters[1]);
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see org.apache.commons.math3.analysis.MultivariateVectorFunction#value(double[])
		 */
		public double[] value(double[] variables)
		{
			increment();
			double[] values = new double[x.size()];
			for (int i = 0; i < values.length; i++)
			{
				values[i] = evaluate(x.get(i), variables[0], variables[1]);
			}
			return values;
		}
	}

	/**
	 * Base implementation of the PC-PALM clustered model. This is used to fit g(r) curves of membrane proteins which
	 * appear to be distributed as per a fluctuations model.
	 * <p>
	 * g(r)peaks = g(r)stoch + g(r)protein
	 * <p>
	 * where
	 * <p>
	 * g(r)stoch = (1/4*pi*s^2*p) * exp(-r^2/4s^2)
	 * <p>
	 * s = average single molecule positional uncertainty (precision)<br>
	 * p = average protein density
	 * <p>
	 * g(r)protein = (A*exp(-r/l)+1) conv g(r)PSF
	 * <p>
	 * A = proportional to density of proteins in the cluster<br>
	 * l = proportional to length of the cluster<br>
	 * conv = a convolution operation
	 * <p>
	 * g(r)PSF = (1/4*pi*s^2) * exp(-r^2/4s^2)
	 * <p>
	 * Note: The clustered model described in the PLoS One paper models g(r)protein using the exponential directly, i.e.
	 * there is no convolution !!!
	 */
	public abstract class ClusteredModelFunction extends BaseModelFunction
	{
		double[] lastValue = null;

		public ClusteredModelFunction()
		{
			super("Clustered Model");
		}

		@Override
		public double[] getValue()
		{
			return lastValue;
		}

		// Adapted from http://commons.apache.org/proper/commons-math/userguide/optimization.html
		// Use the deprecated API since the new one is not yet documented.

		public double[][] jacobian(double[] variables)
		{
			// Compute the gradients using calculus differentiation
			final double sigma = variables[0];
			final double density = variables[1];
			final double range = variables[2];
			final double amplitude = variables[3];
			final double[][] jacobian = new double[x.size()][variables.length];
			lastValue = new double[x.size()];

			for (int i = 0; i < jacobian.length; ++i)
			{
				final double r = this.x.get(i);

				final double a = 1.0 / (4 * Math.PI * density * sigma * sigma);
				final double b = -r * r / (4 * sigma * sigma);
				final double c = FastMath.exp(b);

				final double d = -r / range;
				final double e = FastMath.exp(d);

				// value  = a * c + 
				//          amplitude * e + 1
				lastValue[i] = a * c + amplitude * e + 1;

				// Differentiate with respect to sigma:
				// value' = a' * c + a * c'  [ Product rule ]
				// c = FastMath.exp(b)
				// c' = b' * FastMath.exp(b)     [ Chain rule ]
				// value' = a' * c + a * b' * c
				jacobian[i][0] = (-2 * a / sigma) * c + a * (-2 * b / sigma) * c;

				// Differentiate with respect to density:
				// c' = 0 since density does not feature in c
				// => value' = a' * c
				jacobian[i][1] = (-a / density) * c;

				// Differentiate with respect to range:
				// value' = amplitude * e'
				// e = FastMath.exp(d)
				// e' = d' * FastMath.exp(d)     [ Chain rule ]
				jacobian[i][2] = amplitude * (-1 * d / range) * e;

				// Differentiate with respect to amplitude:
				jacobian[i][3] = e;
			}

			//// Check numerically ...
			//double[][] jacobian2 = jacobian2(variables);
			//for (int i = 0; i < jacobian.length; i++)
			//{
			//	System.out.printf("dSigma = %g : %g = %g. dDensity = %g : %g = %g. dRange = %g : %g = %g. dAmplitude = %g : %g = %g\n", 
			//			jacobian[i][0], jacobian2[i][0], DoubleEquality.relativeError(jacobian[i][0], jacobian2[i][0]), 
			//			jacobian[i][1], jacobian2[i][1], DoubleEquality.relativeError(jacobian[i][1], jacobian2[i][1]), 
			//			jacobian[i][2], jacobian2[i][2], DoubleEquality.relativeError(jacobian[i][2], jacobian2[i][2]), 
			//			jacobian[i][3], jacobian2[i][3], DoubleEquality.relativeError(jacobian[i][3], jacobian2[i][3])
			//			);
			//}

			return jacobian;
		}

		double[][] jacobian2(double[] variables)
		{
			// Compute the gradients using numerical differentiation
			final double sigma = variables[0];
			final double density = variables[1];
			final double range = variables[2];
			final double amplitude = variables[3];
			final double[][] jacobian = new double[x.size()][variables.length];
			lastValue = new double[x.size()];

			final double delta = 0.001;
			double[][] d = new double[variables.length][variables.length];
			for (int i = 0; i < variables.length; i++)
				d[i][i] = delta * Math.abs(variables[i]); // Should the delta be changed for each parameter ?
			for (int i = 0; i < jacobian.length; ++i)
			{
				final double r = this.x.get(i);
				final double value = lastValue[i] = evaluate(r, sigma, density, range, amplitude);
				for (int j = 0; j < variables.length; j++)
				{
					final double value2 = evaluate(r, sigma + d[0][j], density + d[1][j], range + d[2][j], amplitude +
							d[3][j]);
					jacobian[i][j] = (value2 - value) / d[j][j];
				}
			}
			return jacobian;
		}

		/**
		 * Evaluate the correlation function
		 * 
		 * @param r
		 *            The correlation radius
		 * @param sigma
		 *            Average precision
		 * @param density
		 *            Average protein density
		 * @param range
		 *            Range of the cluster
		 * @param amplitude
		 *            Amplitude of the cluster
		 * @return
		 */
		public double evaluate(double r, final double sigma, final double density, final double range,
				final double amplitude)
		{
			final double gr_stoch = (1.0 / (4 * Math.PI * density * sigma * sigma)) *
					FastMath.exp(-r * r / (4 * sigma * sigma));
			final double gr_protein = amplitude * FastMath.exp(-r / range) + 1;
			return gr_stoch + gr_protein;
		}

		/**
		 * Evaluate the correlation function
		 * 
		 * @param r
		 *            The correlation radius
		 * @param parameters
		 *            The parameters
		 * @return
		 */
		public double evaluate(double r, final double[] parameters)
		{
			return evaluate(r, parameters[0], parameters[1], parameters[2], parameters[3]);
		}
	}

	/**
	 * Allow optimisation using Apache Commons Math 3 Gradient Optimiser
	 */
	public class ClusteredModelFunctionGradient extends ClusteredModelFunction implements MultivariateVectorFunction
	{
		// Adapted from http://commons.apache.org/proper/commons-math/userguide/optimization.html
		// Use the deprecated API since the new one is not yet documented.

		/*
		 * (non-Javadoc)
		 * 
		 * @see org.apache.commons.math3.analysis.MultivariateVectorFunction#value(double[])
		 */
		public double[] value(double[] variables)
		{
			increment();
			double[] values = new double[x.size()];
			for (int i = 0; i < values.length; i++)
			{
				values[i] = evaluate(x.get(i), variables[0], variables[1], variables[2], variables[3]);
			}
			return values;
		}
	}

	public class SumOfSquaresModelFunction
	{
		BaseModelFunction f;
		double[] x, y;

		// Cache the value
		double[] lastParameters;
		double lastSS;

		public SumOfSquaresModelFunction(BaseModelFunction f)
		{
			this.f = f;
			x = f.getX();
			y = f.getY();
		}

		public double evaluate(double[] parameters)
		{
			if (sameVariables(parameters))
				return lastSS;

			lastParameters = null;

			double ss = 0;
			for (int i = x.length; i-- > 0;)
			{
				final double dx = f.evaluate(x[i], parameters) - y[i];
				ss += dx * dx;
			}
			return ss;
		}

		/**
		 * Check if the variable match those last used for computation of the value
		 * 
		 * @param parameters
		 * @return True if the variables are the same
		 */
		private boolean sameVariables(double[] parameters)
		{
			if (lastParameters != null)
			{
				for (int i = 0; i < parameters.length; i++)
					if (parameters[i] != lastParameters[i])
						return false;
				return true;
			}
			return false;
		}

		public double[] gradient(double[] parameters)
		{
			// We can compute the jacobian for all the functions.
			// To get the gradient for the SS we need:
			// f(x) = (g(x) - y)^2
			// f'(x) = 2 * (g(x) - y) * g'(x)

			final double[][] jacobian = f.jacobian(parameters);
			final double[] gx = f.getValue();
			lastSS = 0;
			lastParameters = parameters.clone();

			final double[] gradient = new double[parameters.length];
			for (int i = 0; i < x.length; i++)
			{
				final double dx = gx[i] - y[i];
				lastSS += dx * dx;
				final double twodx = 2 * dx;
				for (int j = 0; j < gradient.length; j++)
				{
					final double g1 = twodx * jacobian[i][j];
					gradient[j] += g1;

					//// Check this is correct
					//final double[] p = parameters.clone();
					//final double h = 0.01 * p[j]; // p[j] is unlikely to be zero
					//p[j] += h;
					//final double dx1 = f.evaluate(x[i], p) - y[i];
					//final double ss1 = dx1 * dx1;
					//double delta = p[j];
					//p[j] -= h;
					//delta -= p[j];
					//final double dx2 = f.evaluate(x[i], p) - y[i];
					//final double ss2 = dx2 * dx2;
					//final double g2 = (ss1 - ss2) / delta;
					//if (!gdsc.smlm.fitting.utils.DoubleEquality.almostEqualRelativeOrAbsolute(g1, g2, 1e-2, 1e-5))
					//	System.out.printf("[%d][%d] %f == %f\n", i, j, g1, g2);
				}

			}
			return gradient;
		}
	}

	public class SumOfSquaresMultivariateFunction implements MultivariateFunction
	{
		SumOfSquaresModelFunction f;

		public SumOfSquaresMultivariateFunction(SumOfSquaresModelFunction f)
		{
			this.f = f;
		}

		public double value(double[] point)
		{
			return f.evaluate(point);
		}
	}

	public class SumOfSquaresMultivariateVectorFunction implements MultivariateVectorFunction
	{
		SumOfSquaresModelFunction f;

		public SumOfSquaresMultivariateVectorFunction(SumOfSquaresModelFunction f)
		{
			this.f = f;
		}

		public double[] value(double[] point) throws IllegalArgumentException
		{
			return f.gradient(point);
		}
	}

	/**
	 * Base implementation of the emulsion clustered model. This model assumes a random distribution of non-overlapping
	 * circles in 2D. The molecules can be located at any position within the circles.
	 * <p>
	 * g(r)peaks = g(r)stoch + g(r)protein
	 * <p>
	 * where
	 * <p>
	 * g(r)stoch = (1/4*pi*s^2*p) * exp(-r^2/4s^2)
	 * <p>
	 * s = average single molecule positional uncertainty (precision)<br>
	 * p = average protein density
	 * <p>
	 * g(r)protein = (A*exp(-r/alpha)*cos(pi*r/(2*r0))+1)
	 * <p>
	 * A = proportional to density of proteins in the cluster<br>
	 * alpha = measure of the coherence length between circles<br>
	 * r0 = Average circle radius
	 * <p>
	 * Note: Described in figure 3 of Veatch, et al (2012) Plos One, e31457
	 */
	public abstract class EmulsionModelFunction extends BaseModelFunction
	{
		double[] lastValue = null;

		public EmulsionModelFunction()
		{
			super("Emulsion Clustered Model");
		}

		@Override
		public double[] getValue()
		{
			return lastValue;
		}

		// Adapted from http://commons.apache.org/proper/commons-math/userguide/optimization.html
		// Use the deprecated API since the new one is not yet documented.

		public double[][] jacobian(double[] variables)
		{
			// Compute the gradients using calculus differentiation
			final double sigma = variables[0];
			final double density = variables[1];
			final double range = variables[2];
			final double amplitude = variables[3];
			final double alpha = variables[4];
			final double[][] jacobian = new double[x.size()][variables.length];
			lastValue = new double[x.size()];

			for (int i = 0; i < jacobian.length; ++i)
			{
				final double r = this.x.get(i);

				final double a = 1.0 / (4 * Math.PI * density * sigma * sigma);
				final double b = -r * r / (4 * sigma * sigma);
				final double c = FastMath.exp(b);

				final double d = -r / alpha;
				final double e = FastMath.exp(d);
				final double f = 0.5 * Math.PI * r / range;
				final double g = Math.cos(f);

				// value  = a * c + 
				//          amplitude * e * g + 1
				lastValue[i] = a * c + amplitude * e * g + 1;

				// Differentiate with respect to sigma:
				// value' = a' * c + a * c'  [ Product rule ]
				// c = FastMath.exp(b)
				// c' = b' * FastMath.exp(b)     [ Chain rule ]
				// value' = a' * c + a * b' * c
				jacobian[i][0] = (-2 * a / sigma) * c + a * (-2 * b / sigma) * c;

				// Differentiate with respect to density:
				// c' = 0 since density does not feature in c
				// => value' = a' * c
				jacobian[i][1] = (-a / density) * c;

				// Differentiate with respect to range:
				// value' = amplitude * e * g'
				// g = Math.cos(f)
				// g' = f' * -Math.sin(f)    [ Chain rule ]
				//jacobian[i][2] = amplitude * e * (-f / range) * -Math.sin(f);
				jacobian[i][2] = amplitude * e * (f / range) * Math.sin(f);

				// Differentiate with respect to amplitude:
				jacobian[i][3] = e * g;

				// Differentiate with respect to alpha:
				// value' = amplitude * e' * g
				// e = FastMath.exp(d)
				// e' = d' * FastMath.exp(d)     [ Chain rule ]
				// e' = d' * e               
				jacobian[i][4] = amplitude * (-1 * d / alpha) * e * g;
			}

			//// Check numerically ...
			//double[][] jacobian2 = jacobian2(variables);
			//for (int i = 0; i < jacobian.length; i++)
			//{
			//	System.out.printf("dSigma = %g : %g = %g. dDensity = %g : %g = %g. dRange = %g : %g = %g. dAmplitude = %g : %g = %g. dAlpha = %g : %g = %g\n",
			//					jacobian[i][0], jacobian2[i][0], DoubleEquality.relativeError(jacobian[i][0], jacobian2[i][0]), 
			//					jacobian[i][1], jacobian2[i][1], DoubleEquality.relativeError(jacobian[i][1], jacobian2[i][1]),								
			//					jacobian[i][2], jacobian2[i][2], DoubleEquality.relativeError(jacobian[i][2], jacobian2[i][2]), 
			//					jacobian[i][3],	jacobian2[i][3], DoubleEquality.relativeError(jacobian[i][3], jacobian2[i][3]),
			//					jacobian[i][4], jacobian2[i][4], DoubleEquality.relativeError(jacobian[i][4], jacobian2[i][4]));
			//}

			return jacobian;
		}

		double[][] jacobian2(double[] variables)
		{
			// Compute the gradients using numerical differentiation
			final double sigma = variables[0];
			final double density = variables[1];
			final double range = variables[2];
			final double amplitude = variables[3];
			final double alpha = variables[4];
			final double[][] jacobian = new double[x.size()][variables.length];
			lastValue = new double[x.size()];

			final double delta = 0.001;
			final double[][] d = new double[variables.length][variables.length];
			for (int i = 0; i < variables.length; i++)
				d[i][i] = delta * Math.abs(variables[i]); // Should the delta be changed for each parameter ?
			for (int i = 0; i < jacobian.length; ++i)
			{
				final double r = this.x.get(i);
				final double value = lastValue[i] = evaluate(r, sigma, density, range, amplitude, alpha);
				for (int j = 0; j < variables.length; j++)
				{
					final double value2 = evaluate(r, sigma + d[0][j], density + d[1][j], range + d[2][j], amplitude +
							d[3][j], alpha + d[4][j]);
					jacobian[i][j] = (value2 - value) / d[j][j];
				}
			}
			return jacobian;
		}

		/**
		 * Evaluate the correlation function
		 * 
		 * @param r
		 *            The correlation radius
		 * @param sigma
		 *            Average precision
		 * @param density
		 *            Average protein density
		 * @param range
		 *            Average circle radius
		 * @param amplitude
		 *            Amplitude of the cluster
		 * @param alpha
		 *            Measure of the coherence length between circles
		 * @return
		 */
		public double evaluate(double r, final double sigma, final double density, final double range,
				final double amplitude, final double alpha)
		{
			final double gr_stoch = (1.0 / (4 * Math.PI * density * sigma * sigma)) *
					FastMath.exp(-r * r / (4 * sigma * sigma));
			final double gr_protein = amplitude * FastMath.exp(-r / alpha) * Math.cos(0.5 * Math.PI * r / range) + 1;
			return gr_stoch + gr_protein;
		}

		/**
		 * Evaluate the correlation function
		 * 
		 * @param r
		 *            The correlation radius
		 * @param parameters
		 *            The parameters
		 * @return
		 */
		public double evaluate(double r, final double[] parameters)
		{
			return evaluate(r, parameters[0], parameters[1], parameters[2], parameters[3], parameters[4]);
		}
	}

	/**
	 * Allow optimisation using Apache Commons Math 3 Gradient Optimiser
	 */
	public class EmulsionModelFunctionGradient extends EmulsionModelFunction implements MultivariateVectorFunction
	{
		/*
		 * (non-Javadoc)
		 * 
		 * @see org.apache.commons.math3.analysis.MultivariateVectorFunction#value(double[])
		 */
		public double[] value(double[] variables)
		{
			increment();
			double[] values = new double[x.size()];
			for (int i = 0; i < values.length; i++)
			{
				values[i] = evaluate(x.get(i), variables[0], variables[1], variables[2], variables[3], variables[4]);
			}
			return values;
		}
	}

	private void createResultsTable()
	{
		if (resultsTable == null || !resultsTable.isVisible())
		{
			StringBuilder sb = new StringBuilder();
			sb.append("Model\t");
			sb.append("Colour\t");
			sb.append("Valid\t");
			sb.append("Precision (nm)\t");
			sb.append("Density (um^-2)\t");
			sb.append("Domain Radius (nm)\t");
			// TODO - Find out the units of the domain density
			sb.append("Domain Density\t");
			sb.append("N-cluster\t");
			sb.append("Coherence\t");
			sb.append("cAIC\t");
			resultsTable = new TextWindow(TITLE, sb.toString(), (String) null, 800, 300);
		}
	}

	private void addResult(String model, String resultColour, boolean valid, double precision, double density,
			double domainRadius, double domainDensity, double nCluster, double coherence, double ic)
	{
		StringBuilder sb = new StringBuilder();
		sb.append(model).append("\t");
		sb.append(resultColour.toString()).append("\t");
		sb.append(valid).append("\t");
		sb.append(Utils.rounded(precision, 4)).append("\t");
		sb.append(Utils.rounded(density * 1e6, 4)).append("\t");
		sb.append(getString(domainRadius)).append("\t");
		sb.append(getString(domainDensity)).append("\t");
		sb.append(getString(nCluster)).append("\t");
		sb.append(getString(coherence)).append("\t");
		sb.append(Utils.rounded(ic, 4)).append("\t");
		resultsTable.append(sb.toString());
	}

	private String getString(double value)
	{
		return (value == 0) ? "-" : Utils.rounded(value, 4);
	}
}
