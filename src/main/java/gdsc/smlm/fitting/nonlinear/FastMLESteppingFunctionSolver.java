package gdsc.smlm.fitting.nonlinear;

import gdsc.smlm.fitting.FisherInformationMatrix;
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.FunctionSolverType;
import gdsc.smlm.fitting.MLEFunctionSolver;
import gdsc.smlm.fitting.nonlinear.gradient.FastMLEGradient2Procedure;
import gdsc.smlm.fitting.nonlinear.gradient.FastMLEGradient2ProcedureFactory;
import gdsc.smlm.fitting.nonlinear.gradient.PoissonGradientProcedure;
import gdsc.smlm.fitting.nonlinear.gradient.PoissonGradientProcedureFactory;
import gdsc.smlm.function.ChiSquaredDistributionTable;
import gdsc.smlm.function.Gradient1Function;
import gdsc.smlm.function.Gradient2Function;

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
 * Uses the Fast MLE method to fit a gradient function with coefficients (a).
 * <p>
 * Calculates the Newton-Raphson update vector for a Poisson process using the first and second partial derivatives.
 * <p>
 * Ref: Smith et al, (2010). Fast, single-molecule localisation that achieves theoretically minimum uncertainty.
 * Nature Methods 7, 373-375 (supplementary note), Eq. 12.
 */
public class FastMLESteppingFunctionSolver extends SteppingFunctionSolver implements MLEFunctionSolver
{
	/** The log-likelihood. */
	protected double ll = Double.NaN;
	/** The log-likelihood ratio. */
	protected double llr = Double.NaN;

	/** The gradient procedure. */
	protected FastMLEGradient2Procedure gradientProcedure;

	/**
	 * Create a new stepping function solver.
	 *
	 * @param f
	 *            the function
	 * @param maxRelativeError
	 *            the max relative error
	 * @param maxAbsoluteError
	 *            the max absolute error
	 * @throws NullPointerException
	 *             if the function is null
	 */
	public FastMLESteppingFunctionSolver(Gradient2Function f, double maxRelativeError, double maxAbsoluteError)
	{
		this(f, new ToleranceChecker(maxRelativeError, maxAbsoluteError), null);
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
	public FastMLESteppingFunctionSolver(Gradient2Function f, ToleranceChecker tc, ParameterBounds bounds)
	{
		super(FunctionSolverType.MLE, f, tc, bounds);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#preProcess()
	 */
	@Override
	protected void preProcess()
	{
		ll = llr = Double.NaN;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.SteppingFunctionSolver#prepareFitValue(double[], double[])
	 */
	@Override
	protected double[] prepareFitValue(double[] y, double[] a)
	{
		// Ensure the gradient procedure is created
		y = prepareY(y);
		gradientProcedure = createGradientProcedure(y);
		// Ensure maximisation
		tc.setMinimiseValue(false);
		return y;
	}

	/**
	 * Prepare Y for the gradient procedure by ensuring positive values.
	 *
	 * @param y
	 *            the y
	 * @return the new y
	 */
	protected double[] prepareY(double[] y)
	{
		return ensurePositive(y);
	}

	/**
	 * Creates the gradient procedure.
	 *
	 * @param y
	 *            the y
	 * @return the newton raphson gradient 2 procedure
	 */
	protected FastMLEGradient2Procedure createGradientProcedure(double[] y)
	{
		return FastMLEGradient2ProcedureFactory.create(y, (Gradient2Function) f);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.SteppingFunctionSolver#computeFitValue(double[])
	 */
	@Override
	protected double computeFitValue(double[] a)
	{
		gradientProcedure.computeUpdate(a);

		if (gradientProcedure.isNaNGradients())
			throw new FunctionSolverException(FitStatus.INVALID_GRADIENTS);

		// Log-likelihood only needs to be computed if the tolerance checker 
		// is testing the value
		if (tc.checkValue)
			ll = gradientProcedure.computeLogLikelihood();

		return ll;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.SteppingFunctionSolver#computeStep(double[])
	 */
	@Override
	protected void computeStep(double[] step)
	{
		gradientProcedure.getUpdate(step);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.SteppingFunctionSolver#accept(double, double[], double, double[])
	 */
	@Override
	protected boolean accept(double currentValue, double[] a, double newValue, double[] newA)
	{
		// TODO - Extend the method to implement a combination of Newton-Raphson and Bisection 
		// (see Numerical Recipes in C++, 2nd Ed, page 370, function rtsafe)

		// Always accept the step. The Smith, et al (2010) paper used 10 steps until
		// convergence, with no apparent checking of the log-likelihood value or parameters.
		// The Newton-Raphson method converges fast but does require a good initial
		// estimate for the parameters.
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.SteppingFunctionSolver#prepareFunctionValue(double[], double[])
	 */
	@Override
	protected double[] prepareFunctionValue(double[] y, double[] a)
	{
		y = prepareY(y);
		gradientProcedure = createGradientProcedure(y);
		return y;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.SteppingFunctionSolver#computeFunctionValue(double[], double[])
	 */
	@Override
	protected double computeFunctionValue(double[] y_fit, double[] a)
	{
		ll = gradientProcedure.computeLogLikelihood(a);
		if (y_fit != null)
			System.arraycopy(gradientProcedure.u, 0, y_fit, 0, gradientProcedure.u.length);
		return ll;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.SteppingFunctionSolver#computeValues(double[])
	 */
	@Override
	protected void computeValues(double[] y_fit)
	{
		// Used the cached values
		System.arraycopy(gradientProcedure.u, 0, y_fit, 0, gradientProcedure.u.length);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.SteppingFunctionSolver#computeFisherInformationMatrix()
	 */
	@Override
	protected FisherInformationMatrix computeFisherInformationMatrix()
	{
		// The fisher information is that for a Poisson process
		PoissonGradientProcedure p = PoissonGradientProcedureFactory.create((Gradient1Function) f);
		p.computeFisherInformation(null);
		return new FisherInformationMatrix(p.getLinear(), gradientProcedure.n);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#getValue()
	 */
	@Override
	public double getValue()
	{
		// Override this to return the log likelihood since the value may not 
		// actually be computed during computeFitValue(double[])
		return getLogLikelihood();
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
			ll = gradientProcedure.computeLogLikelihood();
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
		if (Double.isNaN(llr))
		{
			llr = gradientProcedure.computeLogLikelihoodRatio();
		}
		return llr;
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
