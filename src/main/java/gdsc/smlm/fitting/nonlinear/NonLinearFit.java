package gdsc.smlm.fitting.nonlinear;

import java.util.Arrays;

import gdsc.core.utils.DoubleEquality;
import gdsc.smlm.fitting.FisherInformationMatrix;
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.FunctionSolverType;
import gdsc.smlm.fitting.MLEFunctionSolver;
import gdsc.smlm.fitting.WLSEFunctionSolver;
import gdsc.smlm.fitting.linear.EJMLLinearSolver;
import gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator;
import gdsc.smlm.fitting.nonlinear.gradient.GradientCalculatorFactory;
import gdsc.smlm.fitting.nonlinear.gradient.MLEGradientCalculator;
import gdsc.smlm.fitting.nonlinear.stop.ErrorStoppingCriteria;
import gdsc.smlm.function.ChiSquaredDistributionTable;
import gdsc.smlm.function.GradientFunction;
import gdsc.smlm.function.NonLinearFunction;
import gdsc.smlm.function.PoissonCalculator;

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
 * <p>
 * The MLEFunctionSolver is supported if the flag to use a Poisson MLE model is set. If the function supports weights
 * then the WLSEFunctionSolver is supported. The default implementation supports the LSEFunctionSolver.
 */
public class NonLinearFit extends LSEBaseFunctionSolver implements MLEFunctionSolver, WLSEFunctionSolver
{
	protected static final int SUM_OF_SQUARES_BEST = 0;
	protected static final int SUM_OF_SQUARES_OLD = 1;
	protected static final int SUM_OF_SQUARES_NEW = 2;

	protected EJMLLinearSolver solver = new EJMLLinearSolver();
	protected GradientCalculator calculator;
	protected StoppingCriteria sc;

	protected double[] beta = new double[0];
	protected double[] da;
	protected double[] ap = new double[0];

	protected double[][] covar;
	protected double[][] alpha;
	protected double initialLambda = 0.01;
	protected double lambda;
	protected double[] sumOfSquaresWorking;

	protected double initialResidualSumOfSquares;

	protected NonLinearFunction func;
	protected double[] lastY_fit;
	protected double ll = Double.NaN;

	/**
	 * Default constructor
	 * 
	 * @param func
	 *            The function to fit
	 */
	public NonLinearFit(NonLinearFunction func)
	{
		this(func, null);
	}

	/**
	 * Default constructor
	 * 
	 * @param func
	 *            The function to fit
	 * @param sc
	 *            The stopping criteria
	 */
	public NonLinearFit(NonLinearFunction func, StoppingCriteria sc)
	{
		this(func, sc, 1e-3, 1e-10);
	}

	/**
	 * Default constructor
	 * 
	 * @param func
	 *            The function to fit
	 * @param sc
	 *            The stopping criteria
	 * @param maxRelativeError
	 *            Validate the Levenberg-Marquardt fit solution using the specified maximum relative error
	 * @param maxAbsoluteError
	 *            Validate the Levenberg-Marquardt fit solution using the specified maximum absolute error
	 */
	public NonLinearFit(NonLinearFunction func, StoppingCriteria sc, double maxRelativeError, double maxAbsoluteError)
	{
		super(func);
		this.func = func;
		init(sc, maxRelativeError, maxAbsoluteError);
	}

	@Override
	protected void preProcess()
	{
		super.preProcess();
		ll = Double.NaN;
	}

	private void init(StoppingCriteria sc, double maxRelativeError, double maxAbsoluteError)
	{
		setStoppingCriteria(sc);
		solver.setEqual(new DoubleEquality(maxRelativeError, maxAbsoluteError));
	}

	protected boolean nonLinearModel(int n, double[] y, double[] a, boolean initialStage)
	{
		// The NonLinearFunction evaluates a function with parameters a but only computes the gradient
		// for m <= a.length parameters. The parameters can be accessed using the gradientIndices() method.  

		final int[] gradientIndices = f.gradientIndices();
		final int m = gradientIndices.length;

		if (initialStage)
		{
			lambda = initialLambda;
			for (int j = a.length; j-- > 0;)
				ap[j] = a[j];
			sumOfSquaresWorking[SUM_OF_SQUARES_BEST] = calculator.findLinearised(n, y, a, alpha, beta, func);
			initialResidualSumOfSquares = sumOfSquaresWorking[SUM_OF_SQUARES_BEST];
			if (calculator.isNaNGradients())
			{
				//System.out.println("Bad initial gradients");
				return false;
			}
		}

		// Set previous using the current best fit result we have
		sumOfSquaresWorking[SUM_OF_SQUARES_OLD] = sumOfSquaresWorking[SUM_OF_SQUARES_BEST];

		// Solve the gradient equation A x = b:
		// A = Hessian matrix (alpha)
		// x = Parameter shift (output da) 
		// b = Gradient vector (beta)
		if (!solve(a, m))
			return false;

		// Update the parameters. Ensure to use the gradient indices to update the correct parameters
		updateFitParameters(a, gradientIndices, m, da, ap);

		sumOfSquaresWorking[SUM_OF_SQUARES_NEW] = calculator.findLinearised(n, y, ap, covar, da, func);

		if (calculator.isNaNGradients())
		{
			//System.out.println("Bad working gradients");
			return false; // Stop now
			//lambda *= 10.0; // Allow to continue
		}
		else if (sumOfSquaresWorking[SUM_OF_SQUARES_NEW] < sumOfSquaresWorking[SUM_OF_SQUARES_OLD])
		{
			accepted(a, ap, m);
		}
		else
		{
			increaseLambda();
		}

		return true;
	}

	/**
	 * Called when there was a successful improvement of the fit. The lambda parameter should be reduced.
	 *
	 * @param a
	 *            The old fit parameters
	 * @param ap
	 *            The new fit parameters
	 * @param m
	 *            the number of fitted parameters (matches gradient indicies length)
	 */
	protected void accepted(double[] a, double[] ap, int m)
	{
		decreaseLambda();

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

	protected void decreaseLambda()
	{
		lambda *= 0.1;
	}

	protected void increaseLambda()
	{
		lambda *= 10.0;
	}

	/**
	 * Solve the gradient equation A x = b: *
	 * 
	 * <pre>
	 * A = Hessian matrix (alpha)
	 * x = Parameter shift (output da)
	 * b = Gradient vector (beta)
	 * </pre>
	 * 
	 * The Hessian and gradient parameter from the current best scoring parameter set are assumed to be in alpha and
	 * beta. The lambda parameter is used to weight the diagonal of the Hessian.
	 *
	 * @param a
	 *            the current fit parameters
	 * @param m
	 *            the number of fit parameters
	 * @return true, if successful
	 */
	protected boolean solve(double[] a, final int m)
	{
		createLinearProblem(m);
		return solve(covar, da);
	}

	/**
	 * Creates the linear problem.
	 * <p>
	 * The Hessian and gradient parameter from the current best scoring parameter set are assumed to be in alpha and
	 * beta. These are copied into the working variables covar and da. The lambda parameter is used to weight the
	 * diagonal of the Hessian.
	 *
	 * @param m
	 *            the number of fit parameters
	 */
	protected void createLinearProblem(final int m)
	{
		for (int i = m; i-- > 0;)
		{
			da[i] = beta[i];
			for (int j = m; j-- > 0;)
				covar[i][j] = alpha[i][j];
			covar[i][i] *= (1 + lambda);
		}
	}

	/**
	 * Solves (one) linear equation, a x = b
	 * <p>
	 * On input have a[n][n], b[n]. On output b replaced by x[n].
	 * <p>
	 * Note: Any zero elements in b are not solved.
	 * 
	 * @return False if the equation is singular (no solution)
	 */
	protected boolean solve(double[][] a, double[] b)
	{
		// 04-May-2017
		// EJMLLinearSolver was updated to use the pseudoInverse to create a solution
		// when the matrix is singular. Thus we no longer check for zeros.		

		// If the gradient vector is very small set to zero so that this is ignored.

		// TODO - At what level should gradients be ignored (i.e. the parameter has no effect?).
		// Note that analysis on a test dataset showed no difference in results. Those that are caught 
		// for bad gradients must therefore go on to fail on peak filtering criteria. At least this
		// gives the option of not filtering.
		//for (int i = b.length; i-- > 0;)
		//	if (Math.abs(b[i]) < 1e-16)
		//		b[i] = 0;

		// TODO
		// Q. Do we need a better qr decomposition that uses the largest Eigen column first. 
		// There is a version from Apache commons math.
		// We could assess the magnitude of each value in the gradient vector and rearrange.

		return solver.solve(a, b);
	}

	/**
	 * Update the fit parameters. Note that not all parameters are fit and therefore the gradients indices are used to
	 * map the fit parameters to the parameters array.
	 * <p>
	 * This method can be overridden to provide bounded update to the parameters.
	 *
	 * @param a
	 *            the current fit parameters
	 * @param gradientIndices
	 *            the gradient indices (maps the fit parameter index to the parameter array)
	 * @param m
	 *            the number of fit parameters
	 * @param da
	 *            the parameter shift
	 * @param ap
	 *            the new fit parameters
	 */
	protected void updateFitParameters(double[] a, int[] gradientIndices, int m, double[] da, double[] ap)
	{
		for (int j = m; j-- > 0;)
			ap[gradientIndices[j]] = a[gradientIndices[j]] + da[j];
	}

	private FitStatus doFit(int n, double[] y, double[] y_fit, double[] a, double[] a_dev, StoppingCriteria sc)
	{
		final int[] gradientIndices = f.gradientIndices();

		sc.initialise(a);
		if (!nonLinearModel(n, y, a, true))
			return (calculator.isNaNGradients()) ? FitStatus.INVALID_GRADIENTS : FitStatus.SINGULAR_NON_LINEAR_MODEL;
		sc.evaluate(sumOfSquaresWorking[SUM_OF_SQUARES_OLD], sumOfSquaresWorking[SUM_OF_SQUARES_NEW], a);

		while (sc.areNotSatisfied())
		{
			if (!nonLinearModel(n, y, a, false))
				return (calculator.isNaNGradients()) ? FitStatus.INVALID_GRADIENTS
						: FitStatus.SINGULAR_NON_LINEAR_MODEL;

			sc.evaluate(sumOfSquaresWorking[SUM_OF_SQUARES_OLD], sumOfSquaresWorking[SUM_OF_SQUARES_NEW], a);
		}

		if (!sc.areAchieved())
		{
			if (sc.getIteration() >= sc.getMaximumIterations())
				return FitStatus.TOO_MANY_ITERATIONS;
			return FitStatus.FAILED_TO_CONVERGE;
		}

		if (a_dev != null)
		{
			// Copy the diagonal from the current solution alpha (so we do not 
			// recompute it if inversion fails) 
			final double[] I = new double[beta.length];
			for (int i = I.length; i-- > 0;)
				I[i] = alpha[i][i];
			if (!computeDeviations(a_dev))
			{
				// Matrix inversion failed. In order to return a solution 
				// return the reciprocal of the diagonal of the Fisher information 
				// for a loose bound on the limit 
				//final double[] I = calculator.fisherInformationDiagonal(n, a, f);
				Arrays.fill(a_dev, 0);
				for (int i = gradientIndices.length; i-- > 0;)
					a_dev[gradientIndices[i]] = FisherInformationMatrix.reciprocal(I[i]);
			}
		}

		value = sumOfSquaresWorking[SUM_OF_SQUARES_BEST];

		// Compute fitted data points
		if (y_fit != null)
		{
			for (int i = 0; i < n; i++)
				y_fit[i] = func.eval(i);
		}

		return FitStatus.OK;
	}

	/**
	 * Compute the parameter deviations using the covariance matrix of the solution
	 *
	 * @param a_dev
	 *            the a dev
	 * @return true, if successful
	 */
	private boolean computeDeviations(double[] a_dev)
	{
		if (isMLE())
		{
			// The Hessian matrix refers to the log-likelihood ratio.
			// Compute and invert a matrix related to the Poisson log-likelihood.
			// This assumes this does achieve the maximum likelihood estimate for a 
			// Poisson process.
			MLEGradientCalculator c = (MLEGradientCalculator) calculator;
			double[][] I = c.fisherInformationMatrix(lastY.length, null, func);
			
			// Use a dedicated solver optimised for inverting the matrix diagonal 
			FisherInformationMatrix m = new FisherInformationMatrix(I);
			m.setEqual(solver.getEqual());
			setDeviations(a_dev, m.crlb(true));
		}
		else
		{
			// Call the invert method directly on alpha 
			if (!solver.invertSymmPosDef(alpha))
				return false;

			setDeviationsFromMatrix(a_dev, alpha);
		}
		return true;
	}

	/**
	 * Uses Levenberg-Marquardt method to fit a nonlinear model with coefficients (a) for a
	 * set of data points (x, y).
	 * <p>
	 * It is assumed that the data points x[i] corresponding to y[i] are consecutive integers from zero.
	 *
	 * @param y
	 *            Set of n data points to fit (input)
	 * @param y_fit
	 *            Fitted data points (output)
	 * @param a
	 *            Set of m coefficients (input/output)
	 * @param a_dev
	 *            Standard deviation of the set of m coefficients (output)
	 * @return The fit status
	 */
	public FitStatus computeFit(double[] y, double[] y_fit, final double[] a, final double[] a_dev)
	{
		int n = y.length;
		final int nparams = f.gradientIndices().length;

		// Create dynamically for the parameter sizes
		calculator = GradientCalculatorFactory.newCalculator(nparams, isMLE());

		// Initialise storage. 
		// Note that covar and da are passed to EJMLLinerSolver and so must be the correct size. 
		beta = new double[nparams];
		da = new double[nparams];
		covar = new double[nparams][nparams];
		alpha = new double[nparams][nparams];
		ap = new double[a.length];

		// Store the { best, previous, new } sum-of-squares values 
		sumOfSquaresWorking = new double[3];

		boolean copyYfit = false;
		if (isMLE())
		{
			// We must have positive data
			y = ensurePositive(n, y);

			// Store the function values for use in computing the log likelihood
			lastY = y;
			if (y_fit == null)
			{
				// Re-use space
				if (lastY_fit == null || lastY_fit.length < y.length)
					lastY_fit = new double[y.length];
				y_fit = lastY_fit;
				// We will not need to copy y_fit later since lastY_fit is used direct
				copyYfit = false;
			}
		}

		final FitStatus result = doFit(n, y, y_fit, a, a_dev, sc);
		this.evaluations = this.iterations = sc.getIteration();

		if (isMLE())
		{
			// Ensure we have a private copy of the the y_fit since the any calling
			// code may modify it
			if (copyYfit)
			{
				if (lastY_fit == null || lastY_fit.length < y.length)
					lastY_fit = new double[y.length];
				System.arraycopy(y_fit, 0, lastY_fit, 0, y.length);
			}
		}

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
	 * @return the initialResidualSumOfSquares
	 */
	public double getInitialResidualSumOfSquares()
	{
		return initialResidualSumOfSquares;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#setGradientFunction(gdsc.smlm.function.GradientFunction)
	 */
	public void setGradientFunction(GradientFunction f)
	{
		super.setGradientFunction(f);
		if (!(f instanceof NonLinearFunction))
			throw new IllegalArgumentException("Function must be a NonLinearFunction");
	}

	/**
	 * Set the stopping criteria for the {@link #fit(int, double[], double[], double[], double[], double[], double)}
	 * method
	 * 
	 * @param sc
	 */
	public void setStoppingCriteria(StoppingCriteria sc)
	{
		if (sc == null)
			sc = new ErrorStoppingCriteria();
		this.sc = sc;
	}

	/**
	 * Checks if set to perform Maximum Likelihood Estimation assuming Poisson model.
	 *
	 * @return true if is set to perform MLE
	 */
	public boolean isMLE()
	{
		return type == FunctionSolverType.MLE;
	}

	/**
	 * Sets to true to perform Maximum Likelihood Estimation assuming Poisson model.
	 * <p>
	 * This modifies the standard LVM as described in Laurence & Chromy (2010) Efficient maximum likelihood estimator.
	 * Nature Methods 7, 338-339. The input data must be Poisson distributed for this to be relevant.
	 *
	 * @param mle
	 *            true to perform Maximum Likelihood Estimation
	 */
	public void setMLE(boolean mle)
	{
		if (mle)
			type = FunctionSolverType.MLE;
		else
		{
			type = (func.canComputeWeights()) ? FunctionSolverType.WLSE : FunctionSolverType.LSE;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#computeValue(double[], double[], double[])
	 */
	@Override
	public boolean computeValue(double[] y, double[] y_fit, double[] a)
	{
		final int n = y.length;
		final int nparams = f.gradientIndices().length;

		// Create dynamically for the parameter sizes
		calculator = GradientCalculatorFactory.newCalculator(nparams, isMLE());

		if (isMLE())
		{
			// We must have positive data
			y = ensurePositive(n, y);

			// Store the function values for use in computing the log likelihood
			lastY = y;
			if (y_fit == null)
			{
				// Re-use space
				if (lastY_fit == null || lastY_fit.length < y.length)
					lastY_fit = new double[y.length];
				y_fit = lastY_fit;
			}
			else
			{
				lastY_fit = y_fit;
			}
		}

		value = calculator.findLinearised(n, y, y_fit, a, func);

		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.LSEBaseFunctionSolver#getTotalSumOfSquares()
	 */
	@Override
	public double getTotalSumOfSquares()
	{
		if (type == FunctionSolverType.LSE)
			return super.getTotalSumOfSquares();
		throw new IllegalStateException();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.WLSEFunctionSolver#getChiSquared()
	 */
	public double getChiSquared()
	{
		if (type == FunctionSolverType.WLSE)
			// The weighted MLE will produce the chi-squared
			return value;
		throw new IllegalStateException();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.MLEFunctionSolver#getLogLikelihood()
	 */
	public double getLogLikelihood()
	{
		if (type == FunctionSolverType.MLE && lastY != null)
		{
			// The MLE version directly computes the log-likelihood ratio.
			// We must compute the log likelihood for a Poisson MLE.
			if (Double.isNaN(ll))
				ll = PoissonCalculator.logLikelihood(lastY_fit, lastY);
			return ll;
		}
		throw new IllegalStateException();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.MLEFunctionSolver#getLogLikelihoodRatio()
	 */
	public double getLogLikelihoodRatio()
	{
		if (type == FunctionSolverType.MLE)
			// The MLE version directly computes the log-likelihood ratio
			return value;
		throw new IllegalStateException();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.MLEFunctionSolver#getQ()
	 */
	public double getQ()
	{
		if (type == FunctionSolverType.MLE)
			// Value will be the log-likelihood ratio for the MLE.
			// Wilks theorum states the LLR approaches the chi-squared distribution for large n.
			return ChiSquaredDistributionTable.computeQValue(value,
					getNumberOfFittedPoints() - getNumberOfFittedParameters());
		if (type == FunctionSolverType.WLSE)
			// Value will be the Chi-squared
			return ChiSquaredDistributionTable.computeQValue(value,
					getNumberOfFittedPoints() - getNumberOfFittedParameters());
		throw new IllegalStateException();
	}
}
