package gdsc.smlm.fitting.nonlinear;

import gdsc.smlm.fitting.FisherInformationMatrix;
import gdsc.smlm.fitting.FunctionSolverType;
import gdsc.smlm.fitting.MLEFunctionSolver;
import gdsc.smlm.fitting.nonlinear.gradient.LVMGradientProcedure;
import gdsc.smlm.fitting.nonlinear.gradient.MLELVMGradientProcedureFactory;
import gdsc.smlm.fitting.nonlinear.gradient.PoissonGradientProcedure;
import gdsc.smlm.fitting.nonlinear.gradient.PoissonGradientProcedureFactory;
import gdsc.smlm.function.ChiSquaredDistributionTable;
import gdsc.smlm.function.Gradient1Function;
import gdsc.smlm.function.PoissonCalculator;
import gdsc.smlm.function.PrecomputedGradient1Function;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Uses the Levenberg-Marquardt method to fit a gradient function with coefficients (a) using maximum likelihood
 * estimation.
 * <p>
 * This solver computes a modified Chi-squared expression to perform Maximum Likelihood Estimation assuming Poisson
 * model. See Laurence & Chromy (2010) Efficient maximum likelihood estimator. Nature Methods 7, 338-339. The input data
 * must be Poisson distributed for this to be relevant.
 * <p>
 * Per observation variances can be added to both the target x data and the function value to optimise the LLR sCMOS
 * function for Poisson data. See Huang et al, (2015). Video-rate nanoscopy using sCMOS camera–specific single-molecule
 * localization algorithms.
 * Nature Methods 10, 653–658.
 */
public class MLELVMSteppingFunctionSolver extends LVMSteppingFunctionSolver implements MLEFunctionSolver
{
	/** The last Y fit. */
	protected double[] lastY_fit;

	/** The ll. */
	protected double ll = Double.NaN;

	/** The per observation variances. This is not null if fitting using the method of Huang, et al (2015). */
	protected double[] w;
	/**
	 * The gradient function used by the procedure. This may be wrapped to add the per observation variances if fitting
	 * using the LLR sCMOS method of Huang, et al (2015).
	 */
	protected Gradient1Function f1;

	/**
	 * Create a new stepping function solver.
	 *
	 * @param f
	 *            the function
	 * @throws NullPointerException
	 *             if the function is null
	 */
	public MLELVMSteppingFunctionSolver(Gradient1Function f)
	{
		super(FunctionSolverType.LSE, f);
	}

	/**
	 * Create a new stepping function solver.
	 *
	 * @param f
	 *            the function
	 * @param maxRelativeError
	 *            Validate the Levenberg-Marquardt fit solution using the specified maximum relative error
	 * @param maxAbsoluteError
	 *            Validate the Levenberg-Marquardt fit solution using the specified maximum absolute error
	 * @throws NullPointerException
	 *             if the function is null
	 */
	public MLELVMSteppingFunctionSolver(Gradient1Function f, double maxRelativeError, double maxAbsoluteError)
	{
		super(FunctionSolverType.LSE, f, maxRelativeError, maxAbsoluteError);
	}

	/**
	 * Create a new stepping function solver.
	 *
	 * @param f
	 *            the function
	 * @param tc
	 *            the tolerance checker
	 * @param bounds
	 *            the bounds
	 * @throws NullPointerException
	 *             if the function or tolerance checker is null
	 */
	public MLELVMSteppingFunctionSolver(Gradient1Function f, ToleranceChecker tc, ParameterBounds bounds)
	{
		super(FunctionSolverType.LSE, f, tc, bounds);
	}

	/**
	 * Create a new stepping function solver.
	 *
	 * @param f
	 *            the function
	 * @param tc
	 *            the tolerance checker
	 * @param bounds
	 *            the bounds
	 * @param maxRelativeError
	 *            Validate the Levenberg-Marquardt fit solution using the specified maximum relative error
	 * @param maxAbsoluteError
	 *            Validate the Levenberg-Marquardt fit solution using the specified maximum absolute error
	 * @throws NullPointerException
	 *             if the function or tolerance checker is null
	 */
	public MLELVMSteppingFunctionSolver(Gradient1Function f, ToleranceChecker tc, ParameterBounds bounds,
			double maxRelativeError, double maxAbsoluteError)
	{
		super(FunctionSolverType.LSE, f, tc, bounds, maxRelativeError, maxAbsoluteError);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#preProcess()
	 */
	@Override
	protected void preProcess()
	{
		ll = Double.NaN;
		lastY_fit = null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.LVMSteppingFunctionSolver#prepareY(double[])
	 */
	@Override
	protected double[] prepareY(double[] y)
	{
		// We can handle per-observation variances as detailed in
		// Huang, et al. (2015) by simply adding the variances to the target data.

		final int n = y.length;
		w = getWeights(n);
		if (w != null)
		{
			final double[] x = new double[n];
			for (int i = 0; i < n; i++)
			{
				// Also ensure the input y is positive
				x[i] = (y[i] > 0) ? y[i] + w[i] : w[i];
			}
			return x;
		}
		else
			return ensurePositive(y);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.LVMSteppingFunctionSolver#createGradientProcedure(double[])
	 */
	@Override
	protected LVMGradientProcedure createGradientProcedure(double[] y)
	{
		// We can handle per-observation variances as detailed in
		// Huang, et al. (2015) by simply adding the variances to the computed value.
		f1 = (Gradient1Function) f;
		if (w != null)
		{
			f1 = PrecomputedGradient1Function.wrapGradient1Function(f1, w);
		}
		return MLELVMGradientProcedureFactory.create(y, f1);
	}

	@Override
	protected double computeFitValue(double[] a)
	{
		// Cache this
		lastA = a;
		return super.computeFitValue(a);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.LVMSteppingFunctionSolver#computeValues(double[])
	 */
	@Override
	protected void computeValues(double[] y_fit)
	{
		super.computeValues(y_fit);

		// Cache the values to compute the log-likelihood
		int size = f1.size();
		if (lastY_fit == null)
			lastY_fit = new double[size];
		System.arraycopy(y_fit, 0, lastY_fit, 0, size);
		
		if (w != null)
		{
			// The function was wrapped to add the per-observation variances
			// to the computed value, these must be subtracted to get the actual value
			for (int i = 0, n = w.length; i < n; i++)
			{
				y_fit[i] -= w[i];
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.SteppingFunctionSolver#computeFisherInformationMatrix()
	 */
	@Override
	protected FisherInformationMatrix computeFisherInformationMatrix()
	{
		// The Hessian matrix refers to the log-likelihood ratio.
		// Compute and invert a matrix related to the Poisson log-likelihood.
		// This assumes this does achieve the maximum likelihood estimate for a 
		// Poisson process.
		PoissonGradientProcedure p = PoissonGradientProcedureFactory.create(f1);
		p.computeFisherInformation(lastA);
		p.getLinear(walpha);
		return new FisherInformationMatrix(walpha, beta.length);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.MLEFunctionSolver#getLogLikelihood()
	 */
	public double getLogLikelihood()
	{
		if (Double.isNaN(ll))
		{
			if (lastY_fit == null)
			{
				// Evaluate the function values if necessary
				int size = f1.size();
				lastY_fit = new double[size];
				super.computeValues(lastY_fit);
			}

			ll = PoissonCalculator.fastLogLikelihood(lastY_fit, lastY);
		}
		return ll;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.MLEFunctionSolver#getLogLikelihoodRatio()
	 */
	public double getLogLikelihoodRatio()
	{
		// This method computes the log-likelihood ratio directly
		return value;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.MLEFunctionSolver#getQ()
	 */
	public double getQ()
	{
		// Wilks theorum states the LLR approaches the chi-squared distribution for large n.
		return ChiSquaredDistributionTable.computeQValue(getLogLikelihoodRatio(),
				getNumberOfFittedPoints() - getNumberOfFittedParameters());
	}
}
