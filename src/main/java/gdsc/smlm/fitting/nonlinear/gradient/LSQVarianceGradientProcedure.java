package gdsc.smlm.fitting.nonlinear.gradient;

import java.util.Arrays;

import gdsc.smlm.fitting.linear.EJMLLinearSolver;
import gdsc.smlm.function.Gradient1Function;
import gdsc.smlm.function.Gradient1Procedure;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2018 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Compute the variance of the parameters of the function assuming a least squares fit of a Poisson process.
 * <p>
 * Uses the Mortensen formula (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 25.
 */
public class LSQVarianceGradientProcedure implements Gradient1Procedure
{
	public static final int STATUS_OK = 0;
	public static final int STATUS_BAD_GRADIENTS = 1;
	public static final int STATUS_FAILED_INVERSION = 2;

	protected final Gradient1Function func;

	/**
	 * The number of gradients
	 */
	public final int n;
	/** The working size. */
	private final int size;
	/**
	 * Working space for I = sum_i { Ei,a * Ei,b }
	 */
	protected final double[] I;
	/**
	 * Working space for E = sum_i { Ei * Ei,a * Ei,b }
	 */
	protected final double[] E;

	/** The solver. */
	protected final EJMLLinearSolver solver;

	/** The variance. */
	public final double[] variance;

	/**
	 * Instantiates a new LSQ variance gradient procedure.
	 *
	 * @param func
	 *            Gradient function
	 */
	public LSQVarianceGradientProcedure(final Gradient1Function func)
	{
		this(func, EJMLLinearSolver.createForInversion(1e-2));
	}

	/**
	 * Instantiates a new LSQ variance gradient procedure.
	 *
	 * @param func
	 *            Gradient function
	 * @param solver
	 *            the solver
	 * @throws IllegalArgumentException
	 *             if the solver is null
	 */
	public LSQVarianceGradientProcedure(final Gradient1Function func, EJMLLinearSolver solver)
			throws IllegalArgumentException
	{
		if (solver == null)
			throw new IllegalArgumentException("The solver cannot be null");

		this.func = func;
		this.n = func.getNumberOfGradients();
		size = n * (n + 1) / 2;

		I = new double[func.size()];
		E = new double[I.length];
		variance = new double[n];

		this.solver = solver;
	}

	/**
	 * Evaluate the function and compute the variance of the parameters.
	 *
	 * @param a
	 *            Set of coefficients for the function (if null the function must be pre-initialised)
	 * @return the result
	 */
	public int variance(final double[] a)
	{
		initialise();
		if (a != null)
			func.initialise1(a);
		func.forEach((Gradient1Procedure) this);
		if (finish())
			return STATUS_BAD_GRADIENTS;
		if (!solver.invert(I, n))
			return STATUS_FAILED_INVERSION;
		computeVariance();
		return STATUS_OK;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient1Procedure#execute(double, double[])
	 */
	public void execute(final double Ei, double[] Eix)
	{
		for (int a = 0; a < n; a++)
		{
			for (int b = 0, j = a * n; b <= a; b++, j++)
			{
				double v = Eix[a] * Eix[b];
				I[j] += v;
				E[j] += Ei * v;
			}
		}
	}

	/**
	 * Initialise for the computation using the first order gradients.
	 */
	protected void initialise()
	{
		Arrays.fill(I, 0);
		Arrays.fill(E, 0);
		Arrays.fill(variance, 0);
	}

	/**
	 * Finish the computation using the first order gradients. Check the gradients are OK then generate symmetric matrix
	 * for I and E.
	 *
	 * @return true, if the gradient computation failed (e.g. NaN gradients
	 */
	protected boolean finish()
	{
		for (int i = 0; i < size; i++)
			if (Double.isNaN(I[i]))
				return true;
		// Generate symmetric matrix
		for (int a = 0; a < n; a++)
			for (int b = 0; b < a; b++)
			{
				int j = a * n + b;
				int k = b * n + a;
				//	System.out.printf("I[%d] = I[%d];\n", k, j);
				//	System.out.printf("E[%d] = E[%d];\n", k, j);
				I[k] = I[j];
				E[k] = E[j];
			}
		return false;
	}

	/**
	 * Compute the variance using the inverted I and E matrices.
	 */
	protected void computeVariance()
	{
		for (int a = 0; a < n; a++)
		{
			// Note: b==a as we only do the diagonal
			double v = 0;
			for (int ap = 0; ap < n; ap++)
			{
				for (int bp = 0; bp < n; bp++)
				{
					v += I[a * n + ap] * E[ap * n + bp] * I[bp * n + a];
				}
			}
			variance[a] = v;
		}
	}
}
