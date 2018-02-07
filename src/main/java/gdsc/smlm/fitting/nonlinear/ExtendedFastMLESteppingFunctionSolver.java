package gdsc.smlm.fitting.nonlinear;

import java.util.Arrays;

import org.ejml.data.DenseMatrix64F;

import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.linear.EJMLLinearSolver;
import gdsc.smlm.fitting.nonlinear.gradient.FastMLEGradient2Procedure;
import gdsc.smlm.fitting.nonlinear.gradient.FastMLEGradient2ProcedureFactory;
import gdsc.smlm.fitting.nonlinear.gradient.FastMLEJacobianGradient2Procedure;
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
public class ExtendedFastMLESteppingFunctionSolver extends FastMLESteppingFunctionSolver
{
	/** The Jacobian gradient procedure. */
	protected FastMLEJacobianGradient2Procedure jacobianGradientProcedure;
	/** The jacobian. */
	protected double[] jacobian = null;

	protected EJMLLinearSolver solver = null;

	public static final double DEFAULT_MAX_RELATIVE_ERROR = LVMSteppingFunctionSolver.DEFAULT_MAX_RELATIVE_ERROR;
	public static final double DEFAULT_MAX_ABSOLUTE_ERROR = LVMSteppingFunctionSolver.DEFAULT_MAX_ABSOLUTE_ERROR;

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
	public ExtendedFastMLESteppingFunctionSolver(Gradient2Function f, double maxRelativeError, double maxAbsoluteError)
	{
		super(f, maxRelativeError, maxAbsoluteError);
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
	public ExtendedFastMLESteppingFunctionSolver(Gradient2Function f, ToleranceChecker tc, ParameterBounds bounds)
	{
		super(f, tc, bounds);
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
}
