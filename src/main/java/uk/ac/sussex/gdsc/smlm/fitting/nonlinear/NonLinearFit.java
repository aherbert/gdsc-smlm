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
package uk.ac.sussex.gdsc.smlm.fitting.nonlinear;

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.smlm.fitting.FisherInformationMatrix;
import uk.ac.sussex.gdsc.smlm.fitting.FitStatus;
import uk.ac.sussex.gdsc.smlm.fitting.FunctionSolverType;
import uk.ac.sussex.gdsc.smlm.fitting.MLEFunctionSolver;
import uk.ac.sussex.gdsc.smlm.fitting.WLSEFunctionSolver;
import uk.ac.sussex.gdsc.smlm.fitting.linear.EJMLLinearSolver;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.GradientCalculatorFactory;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.stop.ErrorStoppingCriteria;
import uk.ac.sussex.gdsc.smlm.function.ChiSquaredDistributionTable;
import uk.ac.sussex.gdsc.smlm.function.GradientFunction;
import uk.ac.sussex.gdsc.smlm.function.NonLinearFunction;
import uk.ac.sussex.gdsc.smlm.function.PoissonCalculator;

/**
 * Uses Levenberg-Marquardt method to fit a nonlinear model with coefficients (a) for a
 * set of data points (x, y).
 * <p>
 * The MLEFunctionSolver is supported if the flag to use a Poisson MLE model is set. If the function supports weights
 * then the WLSEFunctionSolver is supported. The default implementation supports the LSEFunctionSolver.
 */
public class NonLinearFit extends LSEBaseFunctionSolver implements MLEFunctionSolver, WLSEFunctionSolver
{
	/** Index for the best sum-of-squares in {@link #sumOfSquaresWorking} */
	protected static final int SUM_OF_SQUARES_BEST = 0;
	/** Index for the new sum-of-squares in {@link #sumOfSquaresWorking} */
	protected static final int SUM_OF_SQUARES_OLD = 1;
	/** Index for the previous sum-of-squares in {@link #sumOfSquaresWorking} */
	protected static final int SUM_OF_SQUARES_NEW = 2;

	/** The solver. */
	protected EJMLLinearSolver solver = new EJMLLinearSolver();
	/** The calculator. */
	protected GradientCalculator calculator;
	/** The stopping criteria. */
	protected StoppingCriteria sc;

	/** The beta. */
	protected double[] beta = new double[0];

	/** Working space for beta. */
	protected double[] da;

	/** The updated parameters a.
	 * This is equal to the current parameters a plus the solution x to A x = b
	 * with A = {@link #covar} and b = gradient vector (beta). */
	protected double[] ap = new double[0];

	/**
	 * Working space for the alpha matrix.
	 * This is equal to {@link #alpha} with the diagonal scaled by 1 + {@link #lambda}.
	 */
	protected double[][] covar;

	/** The alpha matrix. */
	protected double[][] alpha;

	/** The initial lambda. */
	protected double initialLambda = 0.01;

	/** The working lambda. */
	protected double lambda;

	/** The sum of squares results. */
	protected double[] sumOfSquaresWorking;

	/** The initial residual sum of squares. */
	protected double initialResidualSumOfSquares;

	/** The function. */
	protected NonLinearFunction func;

	/** The y data from the last successful fit. */
	protected double[] lastyFit;

	/** The log-likelihood. Used for the MLE LVM algorithm. */
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

	/**
	 * Compute a step for the LVM non linear model.
	 *
	 * @param n
	 *            the number of data points
	 * @param y
	 *            the data
	 * @param a
	 *            the parameters
	 * @param initialStage
	 *            the initial stage flag
	 * @return true, if successful
	 */
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
				//System.out.println("Bad initial gradients");
				return false;
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
			//System.out.println("Bad working gradients");
			return false; // Stop now
			//lambda *= 10.0; // Allow to continue
		else if (sumOfSquaresWorking[SUM_OF_SQUARES_NEW] < sumOfSquaresWorking[SUM_OF_SQUARES_OLD])
			accepted(a, ap, m);
		else
			increaseLambda();

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
			beta[j] = da[j];
		for (int j = a.length; j-- > 0;)
			a[j] = ap[j];
		sumOfSquaresWorking[SUM_OF_SQUARES_BEST] = sumOfSquaresWorking[SUM_OF_SQUARES_NEW];
	}

	/**
	 * Decrease lambda. Call this when the solution to the matrix improved the score.
	 */
	protected void decreaseLambda()
	{
		lambda *= 0.1;
	}

	/**
	 * Increase lambda. Call this when the solution to the matrix do not improve the score, or the matrix had no
	 * solution.
	 */
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
	 * @param a
	 *            the a
	 * @param b
	 *            the b
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

	private FitStatus doFit(int n, double[] y, double[] yFit, double[] a, double[] aDev, StoppingCriteria sc)
	{
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

		if (aDev != null)
			if (!computeDeviations(n, y, aDev))
				return FitStatus.SINGULAR_NON_LINEAR_SOLUTION;

		value = sumOfSquaresWorking[SUM_OF_SQUARES_BEST];

		// Compute fitted data points
		if (yFit != null)
			for (int i = 0; i < n; i++)
				yFit[i] = func.eval(i);

		return FitStatus.OK;
	}

	/**
	 * Compute the parameter deviations using the covariance matrix of the solution.
	 *
	 * @param n
	 *            the n
	 * @param y
	 *            the y
	 * @param aDev
	 *            the parameter deviations
	 * @return true, if successful
	 */
	private boolean computeDeviations(int n, double[] y, double[] aDev)
	{
		if (isMLE())
		{
			// The Hessian matrix refers to the log-likelihood ratio.
			// Compute and invert a matrix related to the Poisson log-likelihood.
			// This assumes this does achieve the maximum likelihood estimate for a
			// Poisson process.
			final double[][] I = calculator.fisherInformationMatrix(n, null, func);
			if (calculator.isNaNGradients())
				throw new FunctionSolverException(FitStatus.INVALID_GRADIENTS);

			// Use a dedicated solver optimised for inverting the matrix diagonal.
			final FisherInformationMatrix m = new FisherInformationMatrix(I);
			setDeviations(aDev, m);
			return true;
		}
		final double[] covar = calculator.variance(n, null, func);
		if (covar != null)
		{
			setDeviations(aDev, covar);
			return true;
		}
		return false;
	}

	/**
	 * Uses Levenberg-Marquardt method to fit a nonlinear model with coefficients (a) for a
	 * set of data points (x, y).
	 * <p>
	 * It is assumed that the data points x[i] corresponding to y[i] are consecutive integers from zero.
	 *
	 * @param y
	 *            Set of n data points to fit (input)
	 * @param yFit
	 *            Fitted data points (output)
	 * @param a
	 *            Set of m coefficients (input/output)
	 * @param aDev
	 *            Standard deviation of the set of m coefficients (output)
	 * @return The fit status
	 */
	@Override
	public FitStatus computeFit(double[] y, double[] yFit, final double[] a, final double[] aDev)
	{
		final int n = y.length;
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
			if (yFit == null)
			{
				// Re-use space
				if (lastyFit == null || lastyFit.length < y.length)
					lastyFit = new double[y.length];
				yFit = lastyFit;
				// We will not need to copy yFit later since lastyFit is used direct
				copyYfit = false;
			}
		}

		final FitStatus result = doFit(n, y, yFit, a, aDev, sc);
		this.evaluations = this.iterations = sc.getIteration();

		if (isMLE())
			// Ensure we have a private copy of the the yFit since the any calling
			// code may modify it
			if (copyYfit)
			{
				if (lastyFit == null || lastyFit.length < y.length)
					lastyFit = new double[y.length];
				System.arraycopy(yFit, 0, lastyFit, 0, y.length);
			}

		return result;
	}

	/**
	 * Used for debugging.
	 *
	 * @param format
	 *            the format
	 * @param o
	 *            the o
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
	 * @see uk.ac.sussex.gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#setGradientFunction(uk.ac.sussex.gdsc.smlm.function.GradientFunction)
	 */
	@Override
	public void setGradientFunction(GradientFunction f)
	{
		super.setGradientFunction(f);
		if (!(f instanceof NonLinearFunction))
			throw new IllegalArgumentException("Function must be a NonLinearFunction");
		func = (NonLinearFunction) f;
	}

	/**
	 * Set the stopping criteria for the {@link #fit(double[], double[], double[], double[])}
	 * method.
	 *
	 * @param sc
	 *            the new stopping criteria
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
		return getType() == FunctionSolverType.MLE;
	}

	/**
	 * Sets to true to perform Maximum Likelihood Estimation assuming Poisson model.
	 * <p>
	 * This modifies the standard LVM as described in Laurence &amp; Chromy (2010) Efficient maximum likelihood estimator.
	 * Nature Methods 7, 338-339. The input data must be Poisson distributed for this to be relevant.
	 *
	 * @param mle
	 *            true to perform Maximum Likelihood Estimation
	 */
	public void setMLE(boolean mle)
	{
		if (mle)
			setType(FunctionSolverType.MLE);
		else
			setType((func.canComputeWeights()) ? FunctionSolverType.WLSE : FunctionSolverType.LSE);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#computeValue(double[], double[], double[])
	 */
	@Override
	public boolean computeValue(double[] y, double[] yFit, double[] a)
	{
		final int n = y.length;

		// Create dynamically for the parameter sizes
		calculator = GradientCalculatorFactory.newCalculator(f.getNumberOfGradients(), isMLE());

		if (isMLE())
		{
			// We must have positive data
			y = ensurePositive(n, y);

			// Store the function values for use in computing the log likelihood
			lastY = y;
			if (yFit == null)
			{
				// Re-use space
				if (lastyFit == null || lastyFit.length < y.length)
					lastyFit = new double[y.length];
				yFit = lastyFit;
			}
			else
				lastyFit = yFit;
		}

		value = calculator.findLinearised(n, y, yFit, a, func);

		return true;
	}

	@Override
	public boolean computeDeviations(double[] y, double[] a, double[] aDev)
	{
		calculator = GradientCalculatorFactory.newCalculator(f.getNumberOfGradients(), isMLE());

		if (isMLE())
			return super.computeDeviations(y, a, aDev);

		// LSE computation
		final double[] covar = calculator.variance(y.length, a, func);
		if (covar != null)
		{
			setDeviations(aDev, covar);
			return true;
		}
		return false;
	}

	@Override
	protected FisherInformationMatrix computeFisherInformationMatrix(double[] y, double[] a)
	{
		// Compute and invert a matrix related to the Poisson log-likelihood.
		// This assumes this does achieve the maximum likelihood estimate for a
		// Poisson process.
		final double[][] I = calculator.fisherInformationMatrix(y.length, a, func);
		if (calculator.isNaNGradients())
			throw new FunctionSolverException(FitStatus.INVALID_GRADIENTS);
		return new FisherInformationMatrix(I);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.fitting.nonlinear.LSEBaseFunctionSolver#getTotalSumOfSquares()
	 */
	@Override
	public double getTotalSumOfSquares()
	{
		if (getType() == FunctionSolverType.LSE)
			return super.getTotalSumOfSquares();
		throw new IllegalStateException();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.fitting.WLSEFunctionSolver#getChiSquared()
	 */
	@Override
	public double getChiSquared()
	{
		if (getType() == FunctionSolverType.WLSE)
			// The weighted LSE will produce the chi-squared
			return value;
		throw new IllegalStateException();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.fitting.MLEFunctionSolver#getLogLikelihood()
	 */
	@Override
	public double getLogLikelihood()
	{
		if (getType() == FunctionSolverType.MLE && lastY != null)
		{
			// The MLE version directly computes the log-likelihood ratio.
			// We must compute the log likelihood for a Poisson MLE.
			if (Double.isNaN(ll))
				ll = PoissonCalculator.fastLogLikelihood(lastyFit, lastY);
			return ll;
		}
		throw new IllegalStateException();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.fitting.MLEFunctionSolver#getLogLikelihoodRatio()
	 */
	@Override
	public double getLogLikelihoodRatio()
	{
		if (getType() == FunctionSolverType.MLE)
			// The MLE version directly computes the log-likelihood ratio
			return value;
		throw new IllegalStateException();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.fitting.MLEFunctionSolver#getQ()
	 */
	@Override
	public double getQ()
	{
		if (getType() == FunctionSolverType.MLE)
			// Value will be the log-likelihood ratio for the MLE.
			// Wilks theorum states the LLR approaches the chi-squared distribution for large n.
			return ChiSquaredDistributionTable.computeQValue(value,
					getNumberOfFittedPoints() - getNumberOfFittedParameters());
		if (getType() == FunctionSolverType.WLSE)
			// Value will be the Chi-squared
			return ChiSquaredDistributionTable.computeQValue(value,
					getNumberOfFittedPoints() - getNumberOfFittedParameters());
		throw new IllegalStateException();
	}
}
