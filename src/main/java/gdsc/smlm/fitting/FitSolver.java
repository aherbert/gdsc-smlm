package gdsc.smlm.fitting;

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
 *---------------------------------------------------------------------------*/

/**
 * Define the fitting solver
 */
public enum FitSolver
{
	/**
	 * Custom Levenberg-Marquardt least-squares fitting.
	 * <p>
	 * Uses the Hessian matrix with a Newton optimisation method that requires inversion of the Hessian.
	 */
	LVM("LVM"),
	/**
	 * Custom Levenberg-Marquardt least-squares fitting with weights. The weights require a function that provides the
	 * expected variance for each data point. Without weights the results match the LVM method. 
	 * <p>
	 * Uses the Hessian matrix with a Newton optimisation method that requires inversion of the Hessian.
	 */
	LVM_WEIGHTED("LVM (Weighted)"),
	/**
	 * Apache Commons Math LVM least-squares fitting
	 * <p>
	 * Uses the Jacobian matrix with a quasi-Newton optimisation (that approximates the inverted Hessian).
	 */
	LVM_QUASI_NEWTON("LVM (Quasi-Newton)"),
	/**
	 * Maximum Likelihood Estimator
	 * <p>
	 * Uses a Poisson noise model for the probability density function of the data.
	 */
	MLE("Maximum Likelihood");

	private String name;

	private FitSolver(String name)
	{
		this.name = name;
	}

	@Override
	public String toString()
	{
		return name;
	}

}