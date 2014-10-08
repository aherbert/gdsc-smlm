package gdsc.smlm.fitting.nonlinear;

import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.FunctionSolver;
import gdsc.smlm.fitting.function.NonLinearFunction;
import gdsc.smlm.fitting.linear.EJMLLinearSolver;
import gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator;
import gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator3;
import gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator4;
import gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator5;
import gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator6;
import gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator7;
import gdsc.smlm.fitting.nonlinear.gradient.GradientCalculatorFactory;
import gdsc.smlm.fitting.utils.DoubleEquality;

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
 * 
 * This is an adaption of the C-code contained in the CcpNmr Analysis Program:
 *   CCPN website (http://www.ccpn.ac.uk/). 
 * The CCPN code was based on Numerical Recipes. 
 *---------------------------------------------------------------------------*/

/**
 * Uses Levenberg-Marquardt method to fit a nonlinear model with coefficients (a) for a
 * set of data points (x, y).
 */
public class NonLinearFit implements FunctionSolver
{
	private static final int SUM_OF_SQUARES_BEST = 0;
	private static final int SUM_OF_SQUARES_OLD = 1;
	private static final int SUM_OF_SQUARES_NEW = 2;

	private EJMLLinearSolver solver = new EJMLLinearSolver();
	private GradientCalculator calculator;

	private double[] beta = new double[0];
	private double[] da;
	private float[] ap = new float[0];

	private double[][] covar;
	private double[][] alpha;
	private double initialLambda = 0.01;
	private double lambda;
	private double[] sumOfSquaresWorking;

	private double totalSumOfSquares;
	private double initialResidualSumOfSquares;
	private double finalResidualSumOfSquares;
	private int numberOfFittedParameters;
	private int numberOfFittedPoints;

	/**
	 * Default constructor
	 */
	public NonLinearFit()
	{
		init(3, 1e-10f);
	}

	/**
	 * Default constructor
	 * 
	 * @param significantDigits
	 *            Validate the Levenberg-Marquardt fit solution to the specified number of significant digits
	 * @param maxAbsoluteError
	 *            Validate the Levenberg-Marquardt fit solution using the specified maximum absolute error
	 */
	public NonLinearFit(int significantDigits, float maxAbsoluteError)
	{
		init(significantDigits, maxAbsoluteError);
	}

	private void init(int significantDigits, float maxAbsoluteError)
	{
		solver.setEqual(new DoubleEquality(significantDigits, maxAbsoluteError));
	}

	private boolean nonLinearModel(int n, float[] y, float[] a, NonLinearFunction func, boolean initialStage)
	{
		// The NonLinearFunction evaluates a function with parameters a but only computes the gradient
		// for m <= a.length parameters. The parameters can be accessed using the gradientIndices() method.  

		int[] gradientIndices = func.gradientIndices();
		int m = gradientIndices.length;

		if (initialStage)
		{
			numberOfFittedParameters = m;
			numberOfFittedPoints = n;
			finalResidualSumOfSquares = 0;

			lambda = initialLambda;
			for (int j = a.length; j-- > 0;)
				ap[j] = a[j];
			sumOfSquaresWorking[SUM_OF_SQUARES_BEST] = calculator.findLinearised(n, y, a, alpha, beta, func);
			initialResidualSumOfSquares = sumOfSquaresWorking[SUM_OF_SQUARES_BEST];
			totalSumOfSquares = getSumOfSquares(n, y);
			if (calculator.isNaNGradients())
			{
				//System.out.println("Bad initial gradients");
				return false;
			}
		}

		// Set previous using the current best fit result we have
		sumOfSquaresWorking[SUM_OF_SQUARES_OLD] = sumOfSquaresWorking[SUM_OF_SQUARES_BEST];

		for (int i = m; i-- > 0;)
		{
			// If the gradient vector is very small set to zero so that this is ignored
			// TODO - At what level should gradients be ignored (i.e. the parameter has no effect?)
			// Note that analysis on a test dataset showed no difference in results. Those that are caught 
			// for bad gradients must therefore go on to fail on peak filtering criteria. At least this
			// gives the option of not filtering.
			da[i] = (Math.abs(beta[i]) < 1e-16) ? 0 : beta[i];
			for (int j = m; j-- > 0;)
				covar[i][j] = alpha[i][j];
			covar[i][i] *= 1 + lambda;
		}

		// Solve the gradient equation A x = b:
		// A = Hessian matrix (covar)
		// x = Parameter shift (output da) 
		// b = Gradient vector (beta : input as da)
		if (!solver.solveWithZeros(covar, da))
			return false;

		// Update the parameters. Ensure to use the gradient indices to update the correct parameters
		for (int j = m; j-- > 0;)
			ap[gradientIndices[j]] = (float) (a[gradientIndices[j]] + da[j]);

		sumOfSquaresWorking[SUM_OF_SQUARES_NEW] = calculator.findLinearised(n, y, ap, covar, da, func);

		if (calculator.isNaNGradients())
		{
			//System.out.println("Bad working gradients");
			return false; // Stop now
			//lambda *= 10.0; // Allow to continue
		}
		else if (sumOfSquaresWorking[SUM_OF_SQUARES_NEW] < sumOfSquaresWorking[SUM_OF_SQUARES_OLD])
		{
			lambda *= 0.1;

			for (int i = 0; i < m; i++)
				for (int j = m; j-- > 0;)
					alpha[i][j] = covar[i][j];

			for (int j = m; j-- > 0;)
			{
				beta[j] = da[j];
			}
			for (int j = a.length; j-- > 0;)
			{
				a[j] = ap[j];
			}

			sumOfSquaresWorking[SUM_OF_SQUARES_BEST] = sumOfSquaresWorking[SUM_OF_SQUARES_NEW];
		}
		else
		{
			lambda *= 10.0;
		}

		return true;
	}

	private double getSumOfSquares(final int n, float[] y)
	{
		double sx = 0, ssx = 0;
		for (int i = n; i-- > 0;)
		{
			sx += y[i];
			ssx += y[i] * y[i];
		}
		final double sumOfSquares = ssx - (sx * sx) / (n);
		return sumOfSquares;
	}

	private FitStatus doFit(int n, float[] y, float[] y_fit, float[] a, float[] a_dev, double[] error,
			NonLinearFunction func, StoppingCriteria sc, double noise)
	{
		int[] gradientIndices = func.gradientIndices();
		int nparams = gradientIndices.length;

		sc.initialise(a);
		if (!nonLinearModel(n, y, a, func, true))
			return (calculator.isNaNGradients()) ? FitStatus.INVALID_GRADIENTS_IN_NON_LINEAR_MODEL
					: FitStatus.SINGULAR_NON_LINEAR_MODEL;
		sc.evaluate(sumOfSquaresWorking[SUM_OF_SQUARES_OLD], sumOfSquaresWorking[SUM_OF_SQUARES_NEW], a);

		while (sc.areNotSatisfied())
		{
			if (!nonLinearModel(n, y, a, func, false))
				return (calculator.isNaNGradients()) ? FitStatus.INVALID_GRADIENTS_IN_NON_LINEAR_MODEL
						: FitStatus.SINGULAR_NON_LINEAR_MODEL;

			sc.evaluate(sumOfSquaresWorking[SUM_OF_SQUARES_OLD], sumOfSquaresWorking[SUM_OF_SQUARES_NEW], a);
		}

		if (!sc.areAchieved())
			return FitStatus.FAILED_TO_CONVERGE;

		if (a_dev != null)
		{
			// This is used to calculate the parameter covariance matrix.
			// Solve the gradient matrix corresponding to the best Chi-squared 
			// stored in alpha and beta. 
			if (!solver.solveWithZeros(alpha, beta))
				return FitStatus.SINGULAR_NON_LINEAR_SOLUTION;

			if (!solver.invert(covar))
				return FitStatus.SINGULAR_NON_LINEAR_SOLUTION;

			for (int i = 0; i < nparams; i++)
				a_dev[gradientIndices[i]] = (float) Math.sqrt(Math.max(covar[i][i], 0));
		}

		if (y_fit != null)
		{
			for (int i = 0; i < n; i++)
				y_fit[i] = func.eval(i);
		}

		finalResidualSumOfSquares = sumOfSquaresWorking[SUM_OF_SQUARES_BEST];
		error[0] = getError(finalResidualSumOfSquares, noise, n, numberOfFittedParameters);

		return FitStatus.OK;
	}

	/**
	 * Compute the error
	 * @param residualSumOfSquares
	 * @param noise
	 * @param numberOfFittedPoints
	 * @param numberOfFittedParameters
	 * @return
	 */
	public static double getError(double residualSumOfSquares, double noise, int numberOfFittedPoints,
			int numberOfFittedParameters)
	{
		double error = residualSumOfSquares;

		// Divide by the uncertainty in the individual measurements yi to get the chi-squared
		if (noise > 0)
		{
			error /= numberOfFittedPoints * noise * noise;
		}

		// This updates the chi-squared value to the average error for a single fitted
		// point using the degrees of freedom (N-M)?
		// Note: This matches the mean squared error output from the MatLab fitting code.
		// If a noise estimate was provided for individual measurements then this will be the
		// reduced chi-square (see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2892436/)
		if (numberOfFittedPoints > numberOfFittedParameters)
			error /= (numberOfFittedPoints - numberOfFittedParameters);
		else
			error = 0;

		return error;
	}

	/**
	 * Uses Levenberg-Marquardt method to fit a nonlinear model with coefficients (a) for a
	 * set of data points (x, y).
	 * <p>
	 * It is assumed that the data points x[i] corresponding to y[i] are consecutive integers from zero.
	 * 
	 * @param n
	 *            The number of points to fit, n <= y.length (allows input data y to be used as a buffer)
	 * @param y
	 *            Set of n data points to fit (input)
	 * @param y_fit
	 *            Fitted data points (output)
	 * @param a
	 *            Set of m coefficients (input/output)
	 * @param a_dev
	 *            Standard deviation of the set of m coefficients (output)
	 * @param error
	 *            Output parameter. The Mean-Squared Error (MSE) for the fit if noise is 0. If noise is provided then
	 *            this will be applied to create a reduced chi-square measure.
	 * @param func
	 *            Non-linear fitting function
	 * @param sc
	 *            The stopping criteria for the function
	 * @param noise
	 *            Estimate of the noise in the individual measurements
	 * @return The fit status
	 */
	public FitStatus fit(final int n, final float[] y, final float[] y_fit, final float[] a, final float[] a_dev,
			final double[] error, final NonLinearFunction func, final StoppingCriteria sc, final double noise)
	{
		final int nparams = func.gradientIndices().length;

		// Create dynamically for the parameter sizes
		calculator = GradientCalculatorFactory.newCalculator(nparams);

		// Initialise storage. 
		// Note that covar and da are passed to EJMLLinerSolver and so must be the correct size. 
		beta = new double[nparams];
		da = new double[nparams];
		covar = new double[nparams][nparams];
		alpha = new double[nparams][nparams];
		ap = new float[a.length];

		// Store the { best, previous, new } sum-of-squares values 
		sumOfSquaresWorking = new double[3];

		final FitStatus result = doFit(n, y, y_fit, a, a_dev, error, func, sc, noise);

		return result;
	}

	/**
	 * Used for debugging
	 * 
	 * @param format
	 * @param o
	 */
	void printf(String format, Object... o)
	{
		System.out.printf(format, o);
	}

	/**
	 * @param initialLambda
	 *            the initial lambda for the Levenberg-Marquardt fitting routine
	 */
	public void setInitialLambda(double initialLambda)
	{
		this.initialLambda = initialLambda;
	}

	/**
	 * @return the initialLambda
	 */
	public double getInitialLambda()
	{
		return initialLambda;
	}

	/**
	 * @return the totalSumOfSquares
	 */
	public double getTotalSumOfSquares()
	{
		return totalSumOfSquares;
	}

	/**
	 * @return the initialResidualSumOfSquares
	 */
	public double getInitialResidualSumOfSquares()
	{
		return initialResidualSumOfSquares;
	}

	/**
	 * @return the numberOfFittedParameters
	 */
	public int getNumberOfFittedParameters()
	{
		return numberOfFittedParameters;
	}

	/**
	 * @return the numberOfFittedPoints
	 */
	public int getNumberOfFittedPoints()
	{
		return numberOfFittedPoints;
	}

	private NonLinearFunction func;
	private StoppingCriteria sc;

	/**
	 * Set the non-linear function for the {@link #fit(int, float[], float[], float[], float[], double[], double)}
	 * method
	 * 
	 * @param sc
	 */
	public void setNonLinearFunction(NonLinearFunction func)
	{
		this.func = func;
	}

	/**
	 * Set the stopping criteria for the {@link #fit(int, float[], float[], float[], float[], double[], double)} method
	 * 
	 * @param sc
	 */
	public void setStoppingCriteria(StoppingCriteria sc)
	{
		this.sc = sc;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#fit(int, float[], float[], float[], float[], double[], double)
	 */
	public FitStatus fit(int n, float[] y, float[] y_fit, float[] a, float[] a_dev, double[] error, double noise)
	{
		return fit(n, y, y_fit, a, a_dev, error, func, sc, noise);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#getResidualSumOfSquares()
	 */
	public double getResidualSumOfSquares()
	{
		return finalResidualSumOfSquares;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#getIterations()
	 */
	public int getIterations()
	{
		return sc.getIteration();
	}
}
