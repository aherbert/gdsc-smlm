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
package gdsc.smlm.fitting.nonlinear;

import gdsc.smlm.fitting.FisherInformationMatrix;
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.FunctionSolverType;
import gdsc.smlm.fitting.WLSEFunctionSolver;
import gdsc.smlm.fitting.nonlinear.gradient.LVMGradientProcedure;
import gdsc.smlm.fitting.nonlinear.gradient.WLSQLVMGradientProcedureFactory;
import gdsc.smlm.fitting.nonlinear.gradient.WPoissonGradientProcedure;
import gdsc.smlm.fitting.nonlinear.gradient.WPoissonGradientProcedureFactory;
import gdsc.smlm.function.ChiSquaredDistributionTable;
import gdsc.smlm.function.Gradient1Function;

/**
 * Uses the Levenberg-Marquardt method to fit a gradient function with coefficients (a) using weighted least squares
 * estimation.
 * <p>
 * This solver computes a modified Chi-squared expression to perform Weighted Least Squares Estimation assuming a
 * Poisson model with a Gaussian noise component. The weight per observation is equal to 1/[variance + max(y, 0) + 1].
 * <p>
 * See Ruisheng, et al (2017) Algorithmic corrections for localization microscopy with sCMOS cameras - characterisation
 * of a computationally efficient localization approach. Optical Express 25, Issue 10, pp 11701-11716.
 */
public class WLSELVMSteppingFunctionSolver extends LVMSteppingFunctionSolver implements WLSEFunctionSolver
{
	/**
	 * Create a new stepping function solver.
	 *
	 * @param f
	 *            the function
	 * @throws NullPointerException
	 *             if the function is null
	 */
	public WLSELVMSteppingFunctionSolver(Gradient1Function f)
	{
		super(FunctionSolverType.WLSE, f);
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
	public WLSELVMSteppingFunctionSolver(Gradient1Function f, double maxRelativeError, double maxAbsoluteError)
	{
		super(FunctionSolverType.WLSE, f, maxRelativeError, maxAbsoluteError);
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
	public WLSELVMSteppingFunctionSolver(Gradient1Function f, ToleranceChecker tc, ParameterBounds bounds)
	{
		super(FunctionSolverType.WLSE, f, tc, bounds);
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
	public WLSELVMSteppingFunctionSolver(Gradient1Function f, ToleranceChecker tc, ParameterBounds bounds,
			double maxRelativeError, double maxAbsoluteError)
	{
		super(FunctionSolverType.WLSE, f, tc, bounds, maxRelativeError, maxAbsoluteError);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.LVMSteppingFunctionSolver#createGradientProcedure(double[])
	 */
	@Override
	protected LVMGradientProcedure createGradientProcedure(double[] y)
	{
		return WLSQLVMGradientProcedureFactory.create(y, getWeights(y.length), (Gradient1Function) f);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.SteppingFunctionSolver#computeFisherInformationMatrix()
	 */
	@Override
	protected FisherInformationMatrix computeFisherInformationMatrix(double[] yFit)
	{
		// Get the values if necessary
		if (yFit != null && yFit.length == ((Gradient1Function) f).size())
			computeValues(yFit);

		// TODO. Check if these deviations are correct.
		// Note: Huang et al (2015) compute:
		// Iab = 1 / (uk + var/g^2) * duda * dudb
		// with uk the expected photon count.
		// This will compute: 
		// Iab = 1 / (xk + 1.0 + var/g^2) * duda * dudb
		// with xk the observed photon count.

		// The last Hessian matrix should be stored in the working alpha.
		return new FisherInformationMatrix(walpha, beta.length);
	}

	@Override
	protected FisherInformationMatrix computeFunctionFisherInformationMatrix(double[] y, double[] a)
	{
		// Compute using the scaled Hessian as per the above method.
		// Use a dedicated procedure that omits computing beta.
		WPoissonGradientProcedure p = WPoissonGradientProcedureFactory.create(y, getWeights(y.length),
				(Gradient1Function) f);
		p.computeFisherInformation(a);
		if (p.isNaNGradients())
			throw new FunctionSolverException(FitStatus.INVALID_GRADIENTS);
		return new FisherInformationMatrix(p.getLinear(), p.n);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.WLSEFunctionSolver#getChiSquared()
	 */
	@Override
	public double getChiSquared()
	{
		return value;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.WLSEFunctionSolver#getQ()
	 */
	@Override
	public double getQ()
	{
		// Value will be the Chi-squared
		return ChiSquaredDistributionTable.computeQValue(value,
				getNumberOfFittedPoints() - getNumberOfFittedParameters());
	}
}
