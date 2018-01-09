package gdsc.smlm.fitting.nonlinear;

import java.util.Arrays;

import org.ejml.data.DenseMatrix64F;

import gdsc.core.utils.SimpleArrayUtils;
import gdsc.core.utils.Sort;
import gdsc.smlm.data.NamedObject;
import gdsc.smlm.fitting.FisherInformationMatrix;
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.FunctionSolverType;
import gdsc.smlm.fitting.MLEFunctionSolver;
import gdsc.smlm.fitting.linear.EJMLLinearSolver;
import gdsc.smlm.fitting.nonlinear.gradient.FastMLEGradient2Procedure;
import gdsc.smlm.fitting.nonlinear.gradient.FastMLEGradient2ProcedureFactory;
import gdsc.smlm.fitting.nonlinear.gradient.FastMLEJacobianGradient2Procedure;
import gdsc.smlm.fitting.nonlinear.gradient.PoissonGradientProcedure;
import gdsc.smlm.fitting.nonlinear.gradient.PoissonGradientProcedureFactory;
import gdsc.smlm.function.ChiSquaredDistributionTable;
import gdsc.smlm.function.ExtendedGradient2Function;
import gdsc.smlm.function.Gradient2Function;
import gdsc.smlm.function.PrecomputedExtendedGradient2Function;
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
	/**
	 * Define the method to use when the line search direction is not in the same direction as
	 * that defined by the first derivative gradient.
	 */
	public enum LineSearchMethod implements NamedObject
	{
		/**
		 * Do nothing to handle the incorrect orientation. The default solver action is taken. This may cause the search
		 * to take an invalid move or it may error.
		 */
		NONE("None"),

		/**
		 * Ignore any search direction that is in the opposite direction to the first derivative gradient.
		 */
		IGNORE("Ignore"),

		/**
		 * Progressively ignore any search direction that is in the opposite direction to the first derivative gradient.
		 * Do this in order of the magnitude of the error
		 */
		PARTIAL_IGNORE("Partial ignore");
		
		private final String name;

		private LineSearchMethod(String name)
		{
			this.name = name;
		}

		@Override
		public String toString()
		{
			return name;
		}

		public String getName()
		{
			return name;
		}

		public String getShortName()
		{
			return name;
		}
	}

	protected LineSearchMethod lineSearchMethod = LineSearchMethod.NONE;

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
	/** The Jacobian gradient procedure. */
	protected FastMLEJacobianGradient2Procedure jacobianGradientProcedure;
	/** The jacobian. */
	protected double[] jacobian = null;

	protected EJMLLinearSolver solver = null;

	public static final double DEFAULT_MAX_RELATIVE_ERROR = LVMSteppingFunctionSolver.DEFAULT_MAX_RELATIVE_ERROR;
	public static final double DEFAULT_MAX_ABSOLUTE_ERROR = LVMSteppingFunctionSolver.DEFAULT_MAX_ABSOLUTE_ERROR;

	protected double[] aOld, searchDirection;
	protected boolean firstEvaluation;

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

	/**
	 * Enable computing the Newton-Raphson step using a full Jacobian solution. This requires computation of the
	 * Jacobian of second order partial derivatives with respect to parameters [i,j] and inversion using matrix
	 * decomposition.
	 *
	 * @param enable
	 *            Set to true to enable
	 * @deprecated The computation of the step using the full Jacobian is invalid
	 */
	@Deprecated
	void enableJacobianSolution(boolean enable)
	{
		// This method is defined at the package level for JUnit testing. It is not public
		// as the method does not work.		
		enableJacobianSolution(enable, DEFAULT_MAX_RELATIVE_ERROR, DEFAULT_MAX_ABSOLUTE_ERROR);
	}

	/**
	 * Enable computing the Newton-Raphson step using a full Jacobian solution. This requires computation of the
	 * Jacobian of second order partial derivatives with respect to parameters [i,j] and inversion using matrix
	 * decomposition.
	 *
	 * @param enable
	 *            Set to true to enable
	 * @param maxRelativeError
	 *            Validate the Jacobian solution using the specified maximum relative error
	 * @param maxAbsoluteError
	 *            Validate the Jacobian solution using the specified maximum absolute error
	 * @deprecated The computation of the step using the full Jacobian is invalid
	 */
	@Deprecated
	void enableJacobianSolution(boolean enable, double maxRelativeError, double maxAbsoluteError)
	{
		// This method is defined at the package level for JUnit testing. It is not public
		// as the method does not work.		
		if (enable)
		{
			if (!(f instanceof ExtendedGradient2Function))
				throw new IllegalStateException("Jacobian requires an " + ExtendedGradient2Function.class.getName());
			solver = new EJMLLinearSolver(maxRelativeError, maxAbsoluteError);
		}
		else
			solver = null;
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
		firstEvaluation = true;
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
		if (solver != null)
		{
			if (w != null)
			{
				f2 = PrecomputedExtendedGradient2Function.wrapExtendedGradient2Function((ExtendedGradient2Function) f,
						w);
			}
			jacobian = new double[f2.size()];
			return jacobianGradientProcedure = new FastMLEJacobianGradient2Procedure(y, (ExtendedGradient2Function) f2);
		}
		else
		{
			if (w != null)
			{
				f2 = PrecomputedGradient2Function.wrapGradient2Function(f2, w);
			}
			return FastMLEGradient2ProcedureFactory.create(y, f2);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.SteppingFunctionSolver#computeFitValue(double[])
	 */
	@Override
	protected double computeFitValue(double[] a)
	{
		if (lineSearchMethod != LineSearchMethod.NONE)
		{
			// The code below will adjust the search direction
			if (firstEvaluation)
			{
				// For the first evaluation we store the old value and initialise
				firstEvaluation = false;
				aOld = a.clone();
				searchDirection = new double[a.length];
			}
			else
			{
				// All subsequent calls to computeFitValue() must check the search direction
				for (int i = 0; i < searchDirection.length; i++)
					// Configure the search direction with the full Newton step
					searchDirection[i] = a[i] - aOld[i];
				
				double[] gradient = gradientProcedure.d1;
				final int[] gradientIndices = f.gradientIndices();
				
				double slope = 0.0;
				for (int i = 0; i < gradient.length; i++)
					slope += gradient[i] * searchDirection[gradientIndices[i]];

				if (slope <= 0.0)
				{
					// The slope is invalid so update the position by removing bad 
					// search direction components
					
					switch (lineSearchMethod)
					{
						case IGNORE:
							// Ignore any search direction that is in the opposite direction to the
							// first derivative gradient.
							slope = 0.0;
							for (int i = 0; i < gradient.length; i++)
							{
								double slopeComponent = gradient[i] * searchDirection[gradientIndices[i]];
								if (slopeComponent < 0)
								{
									// Ignore this component
									a[gradientIndices[i]] = aOld[gradientIndices[i]];
								}
								else
								{
									slope += slopeComponent;
								}
							}
							if (slope == 0)
							{
								// No move so just set converged
								tc.setConverged();
								return ll;
								//throw new FunctionSolverException(FitStatus.LINE_SEARCH_ERROR, "No slope");
							}
							break;
							
						case PARTIAL_IGNORE:
							// Progressively ignore any search direction that is in the opposite direction to 
							// the first derivative gradient. Do this in order of the magnitude of the error
							double[] slopeComponents = new double[gradient.length];
							for (int i = 0; i < slopeComponents.length; i++)
								slopeComponents[i] = gradient[i] * searchDirection[gradientIndices[i]];
							int[] indices = SimpleArrayUtils.newArray(slopeComponents.length, 0, 1);
							Sort.sortAscending(indices, slopeComponents);
							int j = 0;
							while (slope <= 0 && j < slopeComponents.length && slopeComponents[indices[j]] <= 0)
							{
								int i = indices[j];
								// Ignore this component
								slope -= slopeComponents[i];
								a[gradientIndices[i]] = aOld[gradientIndices[i]];
								j++;
							}
							if (j == slopeComponents.length)
							{
								// No move so just set converged
								tc.setConverged();
								return ll;
								//throw new FunctionSolverException(FitStatus.LINE_SEARCH_ERROR, "No slope");
							}
							break;
							
						default:
							throw new IllegalStateException("Unknown line search method: " + lineSearchMethod);
					}
				}
			}
		}

		computeGradients(a);

		// Log-likelihood only needs to be computed if the tolerance checker 
		// is testing the value. Use the Pseudo log-likelihood for speed.
		if (tc.checkValue)
		{
			ll = gradientProcedure.computePseudoLogLikelihood();
			isPseudoLogLikelihood = true;
		}

		return ll;
	}

	/**
	 * Compute the gradients for the Newton step using the gradient procedure.
	 *
	 * @param a
	 *            the funtion parameters
	 */
	protected void computeGradients(double[] a)
	{
		if (solver != null)
			jacobianGradientProcedure.computeJacobian(a);
		else
			gradientProcedure.computeSecondDerivative(a);

		if (gradientProcedure.isNaNGradients())
			throw new FunctionSolverException(FitStatus.INVALID_GRADIENTS);
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

		if (solver != null)
		{
			// Solve the Jacobian. This is an implementation of the Newton-Raphson method
			// for systems of non-linear equations (see Numerical Recipes in C++, 2nd Ed, section 9.6)
			// XXX This does not work.
			// The the first order derivatives "are not n independent, arbitrary functions,
			// rather they obey so-called integrability conditions that are highly restrictive".
			// This code is deprecated and may be removed.
			for (int i = 0; i < step.length; i++)
				step[i] = -d1[i];
			jacobianGradientProcedure.getJacobianLinear(jacobian);
			DenseMatrix64F m = DenseMatrix64F.wrap(d1.length, d1.length, jacobian);
			System.out.println(m.toString());
			System.out.println(Arrays.toString(d2));
			if (solver.solve(jacobian, step))
			{
				// XXX - debug the difference
				double[] step2 = new double[d1.length];
				for (int i = 0; i < step.length; i++)
					step2[i] = -d1[i] / d2[i];
				System.out.printf("[%d] Jacobian Step %s vs %s\n", tc.getIterations(), Arrays.toString(step),
						Arrays.toString(step2));
				return;
			}
		}

		// Simple Newton-Raphson update step as per Smith et al, (2010), SI Eq. 13:
		// parameter -> new parameter + delta
		// => new parameter = parameter - delta  
		for (int i = 0; i < step.length; i++)
			step[i] = -d1[i] / d2[i];
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
	protected double computeFunctionValue(double[] yFit, double[] a)
	{
		ll = gradientProcedure.computeLogLikelihood(a);
		isPseudoLogLikelihood = false;
		if (yFit != null)
			copyFunctionValue(yFit);
		return ll;
	}

	/**
	 * Copy the function value into the yFit array.
	 *
	 * @param yFit
	 *            the function values
	 */
	private void copyFunctionValue(double[] yFit)
	{
		final double[] u = gradientProcedure.u;
		if (w != null)
		{
			// The function was wrapped to add the per-observation variances
			// to the computed value, these must be subtracted to get the actual value
			for (int i = 0, n = u.length; i < n; i++)
			{
				yFit[i] = u[i] - w[i];
			}
		}
		else
			System.arraycopy(u, 0, yFit, 0, u.length);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.SteppingFunctionSolver#computeValues(double[])
	 */
	@Override
	protected void computeValues(double[] yFit)
	{
		copyFunctionValue(yFit);
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

	public LineSearchMethod getLineSearchMethod()
	{
		return lineSearchMethod;
	}

	public void setLineSearchMethod(LineSearchMethod lineSearchMethod)
	{
		this.lineSearchMethod = lineSearchMethod;
	}
}
