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
import gdsc.smlm.function.Gradient2Function;
import gdsc.smlm.function.PrecomputedGradient2Function;

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
 * <p>
 * Ref: Huang et al, (2015). Video-rate nanoscopy using sCMOS camera–specific single-molecule localization algorithms.
 * Nature Methods 10, 653–658.
 */
public class FastMLESteppingFunctionSolver extends SteppingFunctionSolver implements MLEFunctionSolver
{
	/** The log-likelihood. */
	protected double ll = Double.NaN;
	/** Flag if the log-likelihood is the pseudo log-likelihood. */
	protected boolean isPseudoLogLikelihood = false;
	/** The log-likelihood ratio. */
	protected double llr = Double.NaN;

	/** The per observation variances. This is not null if fitting using the method of Huang, et al (2015). */
	protected double[] w;
	/**
	 * The gradient function used by the procedure. This may be wrapped to add the per observation variances if fitting
	 * using the method of Huang, et al (2015).
	 */
	protected Gradient2Function f2;
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
		isPseudoLogLikelihood = false;
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

	/**
	 * Creates the gradient procedure.
	 *
	 * @param y
	 *            the y
	 * @return the newton raphson gradient 2 procedure
	 */
	protected FastMLEGradient2Procedure createGradientProcedure(double[] y)
	{
		// We can handle per-observation variances as detailed in
		// Huang, et al. (2015) by simply adding the variances to the computed value.
		f2 = (Gradient2Function) f;
		if (w != null)
		{
			f2 = PrecomputedGradient2Function.wrapGradient2Function(f2, w);
		}
		return FastMLEGradient2ProcedureFactory.create(y, f2);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.SteppingFunctionSolver#computeFitValue(double[])
	 */
	@Override
	protected double computeFitValue(double[] a)
	{
		gradientProcedure.computeSecondDerivative(a);

		if (gradientProcedure.isNaNGradients())
			throw new FunctionSolverException(FitStatus.INVALID_GRADIENTS);

		// Log-likelihood only needs to be computed if the tolerance checker 
		// is testing the value. Use the Pseudo log-likelihood for speed.
		if (tc.checkValue)
		{
			ll = gradientProcedure.computePseudoLogLikelihood();
			isPseudoLogLikelihood = true;
		}

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
		final double[] d1 = gradientProcedure.d1;
		final double[] d2 = gradientProcedure.d2;

		// Simple Newton-Raphson update step as per Smith et al, (2010), SI Eq. 13:
		// parameter -> new parameter + delta
		// => new parameter = parameter - delta  
		for (int i = 0; i < step.length; i++)
			step[i] = -d1[i] / d2[i];
		
		// TODO - Extend this method by implementing Line Search and Backtracking
		// (see Numerical Recipes in C++, 2nd Ed, page 388-389, function lnsrch).
		// This can still be done in the context of the stepping function solver.
		// Adjustments:
		// - Always ensure we compute the function value. This is needed to determine
		// if the step computed a better value.
		// - The first call to computeFitValue() computes derivatives. The value is obtained from 
		// computePseudoLogLikelihood() 
		// - All subsequent calls to computeFitValue() only call computeValue() and implement
		// backtracking if the value does not improve with the full Newton step.
		// - When a suitable value is achieved then this should be returned from computeFitValue()
		// - The next call to compute step must evaluate the derivatives.
		// - The function evaluations counter should be appropriately incremented
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.SteppingFunctionSolver#accept(double, double[], double, double[])
	 */
	@Override
	protected boolean accept(double currentValue, double[] a, double newValue, double[] newA)
	{
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
		isPseudoLogLikelihood = false;
		if (y_fit != null)
			copyFunctionValue(y_fit);
		return ll;
	}

	/**
	 * Copy the function value into the y_fit array.
	 *
	 * @param y_fit
	 *            the function values
	 */
	private void copyFunctionValue(double[] y_fit)
	{
		final double[] u = gradientProcedure.u;
		if (w != null)
		{
			// The function was wrapped to add the per-observation variances
			// to the computed value, these must be subtracted to get the actual value
			for (int i = 0, n = u.length; i < n; i++)
			{
				y_fit[i] = u[i] - w[i];
			}
		}
		else
			System.arraycopy(u, 0, y_fit, 0, u.length);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.SteppingFunctionSolver#computeValues(double[])
	 */
	@Override
	protected void computeValues(double[] y_fit)
	{
		copyFunctionValue(y_fit);
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
		PoissonGradientProcedure p = PoissonGradientProcedureFactory.create(f2);
		p.computeFisherInformation(null); // Assume preinitialised function
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
			isPseudoLogLikelihood = false;
		}
		else if (isPseudoLogLikelihood)
		{
			isPseudoLogLikelihood = false;
			ll -= gradientProcedure.computeLogXFactorialTerm();
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
			llr = gradientProcedure.computeLogLikelihoodRatio(getLogLikelihood());
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
