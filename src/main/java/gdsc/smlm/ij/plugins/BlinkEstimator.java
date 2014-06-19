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

import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.plugins.pcpalm.Molecule;
import gdsc.smlm.ij.plugins.pcpalm.PCPALMMolecules;
import gdsc.smlm.ij.utils.LoggingOptimiserFunction;
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.TraceManager;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.PlugIn;

import java.awt.Color;
import java.util.ArrayList;

import org.apache.commons.math3.analysis.DifferentiableMultivariateVectorFunction;
import org.apache.commons.math3.analysis.MultivariateMatrixFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optimization.PointVectorValuePair;
import org.apache.commons.math3.optimization.general.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.util.Precision;

/**
 * Estimates the flourophore blinking rate from a set of localisations.
 * <p>
 * Uses the method of Annibale, et al (2011). Quantitative Photo Activated Localization Microscopy: Unravelling the
 * Effect of Photoblinking. PLoS ONE 6, e22678.
 */
public class BlinkEstimator implements PlugIn
{
	private static String TITLE = "Blink Estimator";

	private static String inputOption = "";
	private static int s_maxDarkTime = 80;
	private static boolean s_relativeDistance = true;
	private static int histogramBins = 50;
	private static boolean showHistogram = true;
	private static double s_searchDistance = 2.5;
	private static int s_nFittedPoints = 5;
	private static int rangeFittedPoints = 0;
	private static boolean fitIntercept = true;
	private static boolean s_timeAtLowerBound = false;

	private BlinkingFunction blinkingModel;
	private double r2;
	private double adjustedR2;

	public double msPerFrame;
	public int maxDarkTime = s_maxDarkTime;
	public boolean relativeDistance = s_relativeDistance;
	public double searchDistance = s_searchDistance;
	public int nFittedPoints = s_nFittedPoints;
	public boolean timeAtLowerBound = s_timeAtLowerBound;
	public boolean showPlots = false;

	private double[] parameters = null;
	private boolean increaseNFittedPoints = false;

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

		if (!showDialog())
			return;

		MemoryPeakResults results = ResultsManager.loadInputResults(inputOption, true);
		msPerFrame = results.getCalibration().exposureTime;
		Utils.log("%s: %d localisations", TITLE, results.size());

		showPlots = true;
		if (rangeFittedPoints > 0)
		{
			computeFitCurves(results, true);
		}
		else
		{
			computeBlinkingRate(results, true);
		}
	}

	private boolean showDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Compute the blinking rate by fitting counts to dark-time.\nSee Annibale et al (2011) PLos ONE 6, e22678.");
		ResultsManager.addInput(gd, inputOption, InputSource.Memory);

		gd.addNumericField("Max_dark_time (frames)", s_maxDarkTime, 0);
		gd.addCheckbox("Relative_distance", s_relativeDistance);
		gd.addNumericField("Histogram_bins", histogramBins, 0);
		gd.addCheckbox("Show_histogram", showHistogram);
		gd.addSlider("Search_distance", 0.5, 5, s_searchDistance);
		gd.addSlider("Fitted_points", 4, 15, s_nFittedPoints);
		gd.addSlider("Range_of_fitted_points", 0, 15, rangeFittedPoints);
		gd.addCheckbox("Time_at_lower_bound", s_timeAtLowerBound);
		//gd.addCheckbox("Fit_intercept", fitIntercept);
		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		inputOption = gd.getNextChoice();
		maxDarkTime = s_maxDarkTime = (int) gd.getNextNumber();
		relativeDistance = s_relativeDistance = gd.getNextBoolean();
		histogramBins = (int) gd.getNextNumber();
		showHistogram = gd.getNextBoolean();
		searchDistance = s_searchDistance = gd.getNextNumber();
		nFittedPoints = s_nFittedPoints = (int) gd.getNextNumber();
		rangeFittedPoints = (int) gd.getNextNumber();
		timeAtLowerBound = s_timeAtLowerBound = gd.getNextBoolean();
		//fitIntercept = gd.getNextBoolean();

		// Check arguments
		try
		{
			Parameters.isAbove("Max dark time", maxDarkTime, 3);
			Parameters.isAbove("Histogram bins", histogramBins, 1);
			Parameters.isAboveZero("Search distance", searchDistance);
			Parameters.isAbove("n-Fitted points", nFittedPoints, 3);
			Parameters.isPositive("Range of fitted points", rangeFittedPoints);
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		return true;
	}

	private void computeFitCurves(MemoryPeakResults results, boolean verbose)
	{
		// Calculate the counts verses dark time curve
		double[] Ntd = calculateCounts(results, maxDarkTime, searchDistance, relativeDistance, verbose);
		double[] td = calculateTd(Ntd);

		Ntd = shift(Ntd);
		td = shift(td);

		// Fit curve 
		double[] nPoints = new double[rangeFittedPoints + 1];
		double[][] parameters = new double[3][nPoints.length];
		double[] r2 = new double[rangeFittedPoints + 1];
		double[] adjustedR2 = new double[rangeFittedPoints + 1];
		for (int n = 0; n <= rangeFittedPoints; n++)
		{
			nPoints[n] = n + nFittedPoints;
			double[] p = fit(td, Ntd, (int) nPoints[n], false);
			if (p == null)
			{
				// Leave as empty in the output plots 
				continue;
			}
			for (int i = 0; i < p.length; i++)
				parameters[i][n] = p[i];
			r2[n] = this.r2;
			adjustedR2[n] = this.adjustedR2;
		}

		// Plot
		plot("Fitted points", "N", nPoints, parameters[0]);
		plot("Fitted points", "nBlinks", nPoints, parameters[1]);
		plot("Fitted points", "tOff", nPoints, parameters[2]);
		plot("Fitted points", "Adjusted R^2", nPoints, adjustedR2);
	}

	/**
	 * Remove the first element of the array. Return the rest of the array
	 * 
	 * @param d
	 * @return
	 */
	private double[] shift(double[] d)
	{
		if (fitIntercept)
			return d;
		double[] d2 = new double[d.length - 1];
		System.arraycopy(d, 1, d2, 0, d2.length);
		return d2;
	}

	private void plot(String xAxisTitle, String yAxisTitle, double[] x, double[] y)
	{
		String title = TITLE + " " + yAxisTitle;
		Plot plot = new Plot(title, xAxisTitle, yAxisTitle, x, y);
		Utils.display(title, plot);
	}

	double computeBlinkingRate(MemoryPeakResults results, boolean verbose)
	{
		parameters = null;
		increaseNFittedPoints = false;

		// Calculate the counts verses dark time curve
		double[] Ntd = calculateCounts(results, maxDarkTime, searchDistance, relativeDistance, verbose);
		double[] td = calculateTd(Ntd);

		if (verbose)
			Utils.log("  Estimate %.0f molecules at td = %.0f ms", Ntd[0], td[0]);

		Ntd = shift(Ntd);
		td = shift(td);

		// Fit curve
		parameters = fit(td, Ntd, nFittedPoints, verbose);
		if (parameters == null)
			return 0;

		// Display
		if (showPlots)
		{
			String title = TITLE + " Molecule Counts";
			Plot plot = new Plot(title, "td (ms)", "Count", td, Ntd);
			Utils.display(title, plot);

			plot.setColor(Color.red);
			plot.addPoints(blinkingModel.getX(), blinkingModel.value(parameters), Plot.CIRCLE);

			// Add the rest that is not fitted
			double[] xOther = new double[td.length - blinkingModel.size()];
			double[] yOther = new double[xOther.length];
			for (int i = 0, t = blinkingModel.size(); i < xOther.length; i++, t++)
			{
				xOther[i] = td[t];
				yOther[i] = blinkingModel.evaluate(td[t], parameters);
			}

			plot.setColor(Color.blue);
			plot.addPoints(xOther, yOther, Plot.CROSS);
			Utils.display(title, plot);
		}

		// Check if the fitted curve asymptotes above the real curve
		if (blinkingModel.evaluate(td[Ntd.length - 1], parameters) < Ntd[Ntd.length - 1])
		{
			if (verbose)
			{
				Utils.log("  *** Warning ***");
				Utils.log("  Fitted curve does not asymptote above real curve. Increase the number of fitted points to sample more of the overcounting regime");
				Utils.log("  ***************");
			}
			increaseNFittedPoints = true;
		}

		// Blinking rate is 1 + nBlinks
		double blinkingRate = 1 + parameters[1];
		if (verbose)
			Utils.log("  Blinking rate = %s", Utils.rounded(blinkingRate, 4));
		return blinkingRate;
	}

	public double computeBlinkingRate(MemoryPeakResults results)
	{
		return computeBlinkingRate(results, false);
	}

	/**
	 * Calculate the counts of molecules using different dark times. The distance threshold for molecule tracing will be
	 * absolute or relative. If relative it is set using the average precision multiplied by the search distance.
	 * <p>
	 * Note that index 0 corresponds to a t-threshold of 1 in the tracing algorithm, i.e. adjacent frames in the
	 * sequence. This is equivalent to a dark time of (up to) the frame acquisition rate, i.e. the molecule is not
	 * allowed to blink.
	 * 
	 * @param results
	 * @param maxDarkTime
	 * @param searchDistance
	 * @param relativeDistance
	 * @param verbose
	 *            Output log messages
	 * @return the counts of molecules
	 */
	public double[] calculateCounts(MemoryPeakResults results, int maxDarkTime, double searchDistance,
			boolean relativeDistance, boolean verbose)
	{
		double distanceThreshold;
		if (relativeDistance)
		{
			double averagePrecision = calculateAveragePrecision(results, verbose);
			distanceThreshold = averagePrecision * searchDistance / results.getNmPerPixel();
			if (verbose)
				Utils.log("Average precision = %f, Distance threshold = %f px", averagePrecision, distanceThreshold);
		}
		else
		{
			distanceThreshold = searchDistance;
			Utils.log("Distance threshold = %f px", distanceThreshold);
		}

		double[] Ntd = new double[maxDarkTime + 1];

		TraceManager tm = new TraceManager(results);
		IJ.showStatus("Computing counts ...");
		for (int td = 0; td <= maxDarkTime; td++)
		{
			IJ.showProgress(td, maxDarkTime);
			Ntd[td] = tm.traceMolecules(distanceThreshold, td + 1);
		}
		IJ.showProgress(1);
		IJ.showStatus("");

		return Ntd;
	}

	/**
	 * Calculate the dark time corresponding to the molecule counts.
	 * <p>
	 * Note that index 0 corresponds to a t-threshold of 1 in the tracing algorithm, i.e. adjacent frames in the
	 * sequence. This is equivalent to a dark time of (up to) the frame acquisition rate, i.e. the molecule is not
	 * allowed to blink.
	 * <p>
	 * The returned Td values are the lower bounds of the dark time, i.e. t-threshold 1 equals 0 dark frames (0ms),
	 * t-threshold 2 equals 1 dark frame (n ms per frame), etc. This behaviour can be changed by setting the
	 * {@link #timeAtLowerBound} flag to false. Then the time will reflect the upper bounds of the dark time, i.e.
	 * t-threshold 1 equals 0 dark frames (n ms per frame), t-threshold 2 equals 1 dark frame (2n ms per frame), etc.
	 * 
	 * @param Ntd
	 * @return
	 */
	public double[] calculateTd(double[] Ntd)
	{
		double[] td = new double[Ntd.length];
		for (int t = 0; t < td.length; t++)
		{
			if (timeAtLowerBound)
			{
				// Using the lower bounds of the dark time allows the blink estimator to predict the sampled blinks
				// statistic produced by the Create Data plugin.
				td[t] = t * msPerFrame;
			}
			else
			{
				// Adjust for the number of frames that the molecule is allowed to be in a dark state
				td[t] = (t + 1) * msPerFrame;
			}
		}
		return td;
	}

	private double calculateAveragePrecision(MemoryPeakResults results, boolean verbose)
	{
		PCPALMMolecules fitter = new PCPALMMolecules();
		ArrayList<Molecule> molecules = fitter.extractLocalisations(results);
		String title = (verbose) ? TITLE + " Localisation Precision" : null;

		double fittedAverage = fitter.calculateAveragePrecision(molecules, title, histogramBins, true, true);

		// Sense check the precision
		if (fittedAverage < 5 || fittedAverage > 60)
		{
			GenericDialog gd = new GenericDialog(TITLE);
			gd.addMessage("Estimated precision is not within expected bounds.\nPlease enter an estimate:");
			gd.addSlider("Precision", 5, 60, fittedAverage);
			gd.showDialog();
			if (!gd.wasCanceled())
				fittedAverage = gd.getNextNumber();
		}

		// The fitter does checks for a good fit to the histogram so just return the value
		return fittedAverage;
	}

	/**
	 * Fit the dark time to counts of molecules curve. Only use the first n fitted points.
	 * <p>
	 * Calculates:<br/>
	 * N = The number of photoblinking molecules in the sample<br/>
	 * nBlink = The average number of blinks per flourophore<br/>
	 * tOff = The off-time
	 * 
	 * @param td
	 *            The dark time
	 * @param ntd
	 *            The counts of molecules
	 * @param nFittedPoints
	 * @param log
	 *            Write the fitting results to the ImageJ log window
	 * @return The fitted parameters [N, nBlink, tOff], or null if no fit was possible
	 */
	@SuppressWarnings("deprecation")
	public double[] fit(double[] td, double[] ntd, int nFittedPoints, boolean log)
	{
		blinkingModel = new BlinkingFunction();
		blinkingModel.setLogging(true);
		for (int i = 0; i < nFittedPoints; i++)
			blinkingModel.addPoint(td[i], ntd[i]);

		// Different convergence thresholds seem to have no effect on the resulting fit, only the number of
		// iterations for convergence
		double initialStepBoundFactor = 100;
		double costRelativeTolerance = 1e-6;
		double parRelativeTolerance = 1e-6;
		double orthoTolerance = 1e-6;
		double threshold = Precision.SAFE_MIN;
		LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer(initialStepBoundFactor,
				costRelativeTolerance, parRelativeTolerance, orthoTolerance, threshold);
		try
		{
			double[] obs = blinkingModel.getY();

			PointVectorValuePair optimum = optimizer.optimize(1000, blinkingModel, obs, blinkingModel.getWeights(),
					new double[] { ntd[0], 0.1, td[1] });

			blinkingModel.setLogging(false);

			double[] parameters = optimum.getPoint();

			double[] exp = optimum.getValue();
			double mean = 0;
			for (double d : obs)
				mean += d;
			mean /= obs.length;
			double ssResiduals = 0, ssTotal = 0;
			for (int i = 0; i < obs.length; i++)
			{
				ssResiduals += (obs[i] - exp[i]) * (obs[i] - exp[i]);
				ssTotal += (obs[i] - mean) * (obs[i] - mean);
			}

			r2 = 1 - ssResiduals / ssTotal;
			adjustedR2 = getAdjustedCoefficientOfDetermination(ssResiduals, ssTotal, obs.length, parameters.length);

			if (log)
			{
				Utils.log("  Fit %d points. R^2 = %s. Adjusted R^2 = %s", obs.length, Utils.rounded(r2, 4),
						Utils.rounded(adjustedR2, 4));
				Utils.log("  N=%s, nBlink=%s, tOff=%s (%s frames)", Utils.rounded(parameters[0], 4),
						Utils.rounded(parameters[1], 4), Utils.rounded(parameters[2], 4),
						Utils.rounded(parameters[2] / msPerFrame, 4));
			}

			return parameters;
		}
		catch (TooManyEvaluationsException e)
		{
			if (log)
				Utils.log("  Failed to fit %d points", blinkingModel.size());
			return null;
		}
	}

	public double getNMolecules()
	{
		if (parameters != null)
			return parameters[0];
		return 0;
	}

	public double getNBlinks()
	{
		if (parameters != null)
			return parameters[1];
		return 0;
	}

	public double getTOff()
	{
		if (parameters != null)
			return parameters[2];
		return 0;
	}

	/**
	 * @param ssResiduals
	 *            Sum of squared residuals from the model
	 * @param ssTotal
	 *            SStotal is the sum of the squared differences from the mean of the dependent variable (total sum of
	 *            squares)
	 * @param n
	 *            Number of observations
	 * @param d
	 *            Number of parameters in the model
	 * @return
	 */
	private double getAdjustedCoefficientOfDetermination(double ssResiduals, double ssTotal, int n, int d)
	{
		if (n - d - 1 <= 0)
			return 1 - (ssResiduals / ssTotal);
		return 1 - (ssResiduals / ssTotal) * ((n - 1) / (n - d - 1));
	}

	/**
	 * @return the coefficient of determination of the previous fit
	 */
	public double getR2()
	{
		return r2;
	}

	/**
	 * @return the adjusted coefficient of determination of the previous fit
	 */
	public double getAdjustedR2()
	{
		return adjustedR2;
	}

	/**
	 * @return the increaseNFittedPoints
	 */
	public boolean isIncreaseNFittedPoints()
	{
		return increaseNFittedPoints;
	}

	/**
	 * Allow optimisation using Apache Commons Math 3 Gradient Optimiser
	 * <p>
	 * N(td) = N . (1 + nBlink . exp((1-td)/tOff)
	 * <p>
	 * where
	 * <p>
	 * N(td) = The number of calculated molecules at different dark times (td)<br/>
	 * N = The number of photoblinking molecules in the sample<br/>
	 * nBlink = The average number of blinks per flourophore<br/>
	 * td = The dark time<br/>
	 * tOff = The off-time<br/>
	 */
	@SuppressWarnings("deprecation")
	public class BlinkingFunction extends LoggingOptimiserFunction implements DifferentiableMultivariateVectorFunction
	{
		public BlinkingFunction()
		{
			super("Blinking Model");
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see gdsc.smlm.ij.plugins.OptimiserFunction#getWeights()
		 */
		public double[] getWeights()
		{
			// Bias the early values
			//double[] w = new double[y.size()];
			//for (int i = 0; i < w.length; i++)
			//	w[i] = w.length - i;
			//return w;
			return super.getWeights();
		}

		// Adapted from http://commons.apache.org/proper/commons-math/userguide/optimization.html
		// Use the deprecated API since the new one is not yet documented.

		private double[][] jacobian(double[] variables)
		{
			// Compute the gradients using calculus differentiation
			final double N = variables[0];
			final double nBlink = variables[1];
			final double tOff = variables[2];
			double[][] jacobian = new double[x.size()][variables.length];

			for (int i = 0; i < jacobian.length; ++i)
			{
				double td = this.x.get(i);

				final double a = (1 - td) / tOff;
				final double b = Math.exp(a);

				// value  = N * (1 + nBlink * b)
				//        = N + N * nBlink * exp(a)

				// Differentiate with respect to N:
				jacobian[i][0] = 1 + nBlink * b;

				// Differentiate with respect to nBlink:
				jacobian[i][1] = N * b;

				// Differentiate with respect to tOff:
				jacobian[i][2] = N * nBlink * b * -a / tOff;
			}

			//// Check numerically ...
			//double[][] jacobian2 = jacobian2(variables);
			//for (int i = 0; i < jacobian.length; i++)
			//{
			//	System.out.printf("N = %g : %g = %g. nBlink = %g : %g = %g. tOff = %g : %g = %g\n", jacobian[i][0],
			//			jacobian2[i][0], DoubleEquality.relativeError(jacobian[i][0], jacobian2[i][0]), jacobian[i][1],
			//			jacobian2[i][1], DoubleEquality.relativeError(jacobian[i][1], jacobian2[i][1]), jacobian[i][2],
			//			jacobian2[i][2], DoubleEquality.relativeError(jacobian[i][2], jacobian2[i][2]));
			//}

			return jacobian;
		}

		@SuppressWarnings("unused")
		private double[][] jacobian2(double[] variables)
		{
			// Compute the gradients using numerical differentiation
			final double N = variables[0];
			final double nBlink = variables[1];
			final double tOff = variables[2];
			double[][] jacobian = new double[x.size()][variables.length];

			final double delta = 0.001;
			double[][] d = new double[variables.length][variables.length];
			for (int i = 0; i < variables.length; i++)
				d[i][i] = delta * Math.abs(variables[i]); // Should the delta be changed for each parameter ?
			for (int i = 0; i < jacobian.length; ++i)
			{
				double r = this.x.get(i);
				double value = evaluate(r, N, nBlink, tOff);
				for (int j = 0; j < variables.length; j++)
				{
					double value2 = evaluate(r, N + d[0][j], nBlink + d[1][j], tOff + d[2][j]);
					jacobian[i][j] = (value2 - value) / d[j][j];
				}
			}
			return jacobian;
		}

		/**
		 * Evaluate the function
		 * 
		 * @param td
		 *            The dark time
		 * @param N
		 *            The number of molecules in the sample
		 * @param nBlink
		 *            The blinking rate
		 * @param tOff
		 *            The off-time
		 * @return
		 */
		public double evaluate(double td, double N, double nBlink, double tOff)
		{
			return N * (1.0 + nBlink * Math.exp((1 - td) / tOff));
		}

		public double evaluate(double td, double[] parameters)
		{
			return evaluate(td, parameters[0], parameters[1], parameters[2]);
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
				values[i] = evaluate(x.get(i), variables[0], variables[1], variables[2]);
			}
			return values;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see org.apache.commons.math3.analysis.DifferentiableMultivariateVectorFunction#jacobian()
		 */
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
	}
}
