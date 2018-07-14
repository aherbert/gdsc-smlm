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

import gdsc.core.utils.BitFlags;
import gdsc.smlm.fitting.FisherInformationMatrix;
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.FunctionSolverType;
import gdsc.smlm.function.Gradient1Function;
import gdsc.smlm.function.Gradient1FunctionStore;
import gdsc.smlm.function.GradientFunction;
import gdsc.smlm.function.ValueFunction;
import gdsc.smlm.function.ValueProcedure;

/**
 * Abstract class for FunctionSolvers that use update steps to the current parameters.
 */
public abstract class SteppingFunctionSolver extends BaseFunctionSolver
{
	/**
	 * Simple class to allow the values to be computed.
	 */
	private static class SimpleValueProcedure implements ValueProcedure
	{
		int i = 0;
		double[] yFit;

		SimpleValueProcedure(double[] yFit)
		{
			this.yFit = yFit;
		}

		@Override
		public void execute(double value)
		{
			yFit[i++] = value;
		}
	}

	protected int[] gradientIndices;
	protected final ToleranceChecker tc;
	protected final ParameterBounds bounds;
	private double[] weights = null;

	/**
	 * Create a new stepping function solver
	 *
	 * @param type
	 *            the type
	 * @param f
	 *            the function
	 * @throws NullPointerException
	 *             if the function is null
	 */
	public SteppingFunctionSolver(FunctionSolverType type, Gradient1Function f)
	{
		this(type, f, new ToleranceChecker(1e-3, 1e-6), null);
	}

	/**
	 * Create a new stepping function solver
	 *
	 * @param type
	 *            the type
	 * @param f
	 *            the function
	 * @param tc
	 *            the tolerance checker
	 * @param bounds
	 *            the bounds
	 * @throws NullPointerException
	 *             if the function or tolerance checker is null
	 * @throws IllegalArgumentException
	 *             if the bounds are not constructed with the same gradient function
	 */
	public SteppingFunctionSolver(FunctionSolverType type, Gradient1Function f, ToleranceChecker tc,
			ParameterBounds bounds)
	{
		super(type, f);
		if (tc == null)
			throw new NullPointerException("Null tolerance checker");
		this.tc = tc;
		if (bounds == null)
			bounds = new ParameterBounds(f);
		else if (bounds.getGradientFunction() != f)
			throw new IllegalArgumentException("Bounds must be constructed with the same gradient function");
		this.bounds = bounds;
	}

	/**
	 * Compute fit.
	 *
	 * @param y
	 *            the y
	 * @param yFit
	 *            the yFit
	 * @param a
	 *            the a
	 * @param aDev
	 *            the aDev
	 * @return the fit status
	 */
	@Override
	protected FitStatus computeFit(double[] y, double[] yFit, double[] a, double[] aDev)
	{
		// Lay out a simple iteration loop for a stepping solver.
		// The sub-class must compute the next step.
		// This class handles attenuation of the step.
		// The sub-class determines if the step is accepted or rejected.

		gradientIndices = f.gradientIndices();
		final double[] step = new double[gradientIndices.length];
		final double[] newA = a.clone();

		// Initialise for fitting
		bounds.initialise();
		tc.reset();
		String name = this.getClass().getSimpleName();

		try
		{
			lastY = prepareFitValue(y, a);

			// First evaluation
			double currentValue = computeFitValue(a);
			log("%s Value [%s] = %s : %s\n", name, tc.getIterations(), currentValue, a);

			int status = 0;
			while (true)
			{
				// Compute next step
				computeStep(step);
				log("%s Step [%s] = %s\n", name, tc.getIterations(), step);

				// Apply bounds to the step
				bounds.applyBounds(a, step, newA);

				// Evaluate
				double newValue = computeFitValue(newA);
				log("%s Value [%s] = %s : %s\n", name, tc.getIterations(), newValue, newA);

				// Check stopping criteria
				status = tc.converged(currentValue, a, newValue, newA);
				log("%s Status [%s] = %s\n", name, tc.getIterations(), status);
				if (status != 0)
				{
					value = newValue;
					System.arraycopy(newA, 0, a, 0, a.length);
					break;
				}

				// Check if the step was an improvement
				if (accept(currentValue, a, newValue, newA))
				{
					log("%s Accepted [%s]\n", name, tc.getIterations());
					currentValue = newValue;
					System.arraycopy(newA, 0, a, 0, a.length);
					bounds.accepted(a, newA);
				}
			}

			log("%s End [%s] = %s\n", name, tc.getIterations(), status);

			if (BitFlags.anySet(status, ToleranceChecker.STATUS_CONVERGED))
			{
				log("%s Converged [%s]\n", name, tc.getIterations());
				// A solver may compute both at the same time...
				if (aDev != null)
					computeDeviationsAndValues(aDev, yFit);
				else if (yFit != null)
					computeValues(yFit);
				return FitStatus.OK;
			}

			// Check the iterations
			if (BitFlags.areSet(status, ToleranceChecker.STATUS_MAX_ITERATIONS))
				return FitStatus.TOO_MANY_ITERATIONS;

			// We should not reach here unless we missed something
			return FitStatus.FAILED_TO_CONVERGE;
		}
		catch (FunctionSolverException e)
		{
			// XXX - debugging
			String msg = e.getMessage();
			if (msg != null)
				System.out.printf("%s failed: %s - %s\n", getClass().getSimpleName(), e.fitStatus.getName(), msg);
			else
				System.out.printf("%s failed: %s\n", getClass().getSimpleName(), e.fitStatus.getName());
			return e.fitStatus;
		}
		finally
		{
			iterations = tc.getIterations();
			// Allow subclasses to increment this
			if (evaluations == 0)
				evaluations = iterations;
		}
	}

	/**
	 * Log progress from the solver.
	 *
	 * @param format
	 *            the format
	 * @param args
	 *            the args
	 */
	private void log(String format, Object... args)
	{
		//		// Convert arrays to a single string
		//		for (int i=0; i<args.length; i++)
		//			if (args[i] instanceof double[])
		//				args[i] = java.util.Arrays.toString((double[])args[i]);
		//		System.out.printf(format, args);
	}

	/**
	 * Prepare y for fitting, e.g. ensure strictly positive values.
	 *
	 * @param y
	 *            the y
	 * @param a
	 *            the parameters
	 * @return the new y
	 */
	protected abstract double[] prepareFitValue(double[] y, double[] a);

	/**
	 * Compute the fit value using the parameters. The y data is the same as that passed to
	 * {@link #prepareFitValue(double[], double[])}.
	 * <p>
	 * This method is followed by a call to {@link #computeStep(double[])} so the step could be pre-computed here.
	 *
	 * @param a
	 *            the parameters
	 * @return the fit value
	 */
	protected abstract double computeFitValue(double[] a);

	/**
	 * Compute the update step for the current parameters.
	 *
	 * @param step
	 *            the step
	 */
	protected abstract void computeStep(double[] step);

	/**
	 * Determine if the step should be accepted. If accepted then the current parameters and function value are updated
	 * and any bounds on the step size may be updated.
	 * <p>
	 * Note that although this class handles convergence on the value/parameters it is left to the sub-class to
	 * determine if each step should be accepted.
	 *
	 * @param currentValue
	 *            the current value
	 * @param a
	 *            the current parameters
	 * @param newValue
	 *            the new value
	 * @param newA
	 *            the new parameters
	 * @return true, if successful
	 */
	protected abstract boolean accept(double currentValue, double[] a, double newValue, double[] newA);

	/**
	 * Compute the deviations for the parameters a from the last call to
	 * {@link #computeFitValue(double[])}. Optionally store the function values.
	 *
	 * @param aDev
	 *            the parameter deviations
	 * @param yFit
	 *            the y fit (may be null)
	 */
	protected void computeDeviationsAndValues(double[] aDev, double[] yFit)
	{
		// Use a dedicated solver optimised for inverting the matrix diagonal.
		// The last Hessian matrix should be stored in the working alpha.
		final FisherInformationMatrix m = computeFisherInformationMatrix(yFit);

		setDeviations(aDev, m);
	}

	/**
	 * Compute the Fisher Information matrix for the parameters a from the last call to
	 * {@link #computeFitValue(double[])}. This can be used to set the covariances for each of the fitted
	 * parameters.
	 * <p>
	 * Alternatively a sub-class can override {@link #computeDeviationsAndValues(double[], double[])} directly and
	 * provide a dummy implementation of this function as it will not be used, e.g. throw an exception.
	 *
	 * @param yFit
	 *            the y fit (may be null)
	 * @return the Fisher Information matrix
	 */
	protected abstract FisherInformationMatrix computeFisherInformationMatrix(double[] yFit);

	/**
	 * Compute the function y-values using the y and parameters a from the last call to
	 * {@link #computeFitValue(double[])}.
	 * <p>
	 * Utility method to compute the function values using the preinitialised function.
	 * Sub-classes may override this if they have cached the function values from the
	 * last execution of a forEach procedure.
	 * <p>
	 * The base gradient function is used. If sub-classes wrap the function (e.g. with per-observation weights) then
	 * these will be omitted.
	 *
	 * @param yFit
	 *            the y fit values
	 */
	protected void computeValues(double[] yFit)
	{
		ValueFunction f = (ValueFunction) this.f;
		f.forEach(new SimpleValueProcedure(yFit));
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#computeValue(double[], double[], double[])
	 */
	@Override
	protected boolean computeValue(double[] y, double[] yFit, double[] a)
	{
		// If the yFit array is not null then wrap the gradient function.
		// Compute the value and the wrapper will store the values appropriately.
		// Then reset the gradient function.

		// Note: If a sub class wraps the function with weights
		// then the weights will not be stored in the function value.
		// Only the value produced by the original function is stored:
		// Wrapped (+weights) < FunctionStore < Function

		// However if the base function is already wrapped then this will occur:
		// Wrapped (+weights) < FunctionStore < Wrapped (+precomputed) < Function

		gradientIndices = f.gradientIndices();
		if (yFit != null && yFit.length == ((Gradient1Function) f).size())
		{
			GradientFunction tmp = f;
			f = new Gradient1FunctionStore((Gradient1Function) f, yFit, null);
			lastY = prepareFunctionValue(y, a);
			value = computeFunctionValue(a);
			f = tmp;
		}
		else
		{
			lastY = prepareFunctionValue(y, a);
			value = computeFunctionValue(a);
		}
		return true;
	}

	/**
	 * Prepare y for computing the function value, e.g. ensure strictly positive values.
	 *
	 * @param y
	 *            the y
	 * @param a
	 *            the parameters
	 * @return the new y
	 */
	protected abstract double[] prepareFunctionValue(double[] y, double[] a);

	/**
	 * Compute the function value. The y data is the same as that passed to
	 * {@link #prepareFunctionValue(double[], double[])}
	 *
	 * @param a
	 *            the parameters
	 * @return the function value
	 */
	protected abstract double computeFunctionValue(double[] a);

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#computeFisherInformationMatrix(double[], double[])
	 */
	@Override
	protected FisherInformationMatrix computeFisherInformationMatrix(double[] y, double[] a)
	{
		gradientIndices = f.gradientIndices();
		y = prepareFunctionFisherInformationMatrix(y, a);
		return computeFunctionFisherInformationMatrix(y, a);
	}

	/**
	 * Prepare y for computing the Fisher information matrix, e.g. ensure strictly positive values.
	 *
	 * @param y
	 *            the y
	 * @param a
	 *            the parameters
	 * @return the new y
	 */
	protected abstract double[] prepareFunctionFisherInformationMatrix(double[] y, double[] a);

	/**
	 * Compute the Fisher information matrix.
	 *
	 * @param y
	 *            the y
	 * @param a
	 *            the parameters
	 * @return the Fisher Information matrix
	 */
	protected abstract FisherInformationMatrix computeFunctionFisherInformationMatrix(double[] y, double[] a);

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.fitting.FunctionSolver#isBounded()
	 */
	@Override
	public boolean isBounded()
	{
		// Bounds are tighter than constraints and we support those
		return true;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.fitting.FunctionSolver#setBounds(double[], double[])
	 */
	@Override
	public void setBounds(double[] lower, double[] upper)
	{
		bounds.setBounds(lower, upper);
	}

	/**
	 * Warning: If the function is changed then the clamp values may require updating. However setting a new function
	 * does not set the clamp values to null to allow caching when the clamp values are unchanged, e.g. evaluation of a
	 * different function in the same parameter space.
	 * <p>
	 * Setting a new function removes the current bounds.
	 *
	 * @param f
	 *            the new gradient function
	 * @see gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#setGradientFunction(gdsc.smlm.function.GradientFunction)
	 */
	@Override
	public void setGradientFunction(GradientFunction f)
	{
		super.setGradientFunction(f);
		bounds.setGradientFunction(f);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#isWeighted()
	 */
	@Override
	public boolean isWeighted()
	{
		return true;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#setWeights(double[])
	 */
	@Override
	public void setWeights(double[] weights)
	{
		this.weights = weights;
	}

	/**
	 * Gets the weights for observations of size n, e.g. the per observation variance term.
	 *
	 * @param n
	 *            the size
	 * @return the weights
	 */
	public double[] getWeights(int n)
	{
		return (weights == null || weights.length != n) ? null : weights;
	}
}
