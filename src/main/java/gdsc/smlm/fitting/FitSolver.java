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
	LVM("Least Squares Estimator", "LSE"),
	/**
	 * Custom Levenberg-Marquardt least-squares fitting with weights. The weights require a function that provides the
	 * expected variance for each data point. Without weights the results match the LVM method.
	 * <p>
	 * Uses the Hessian matrix with a Newton optimisation method that requires inversion of the Hessian.
	 */
	LVM_WEIGHTED("Weighted Least Squares Estimator", "WLSE"),
	/**
	 * Apache Commons Math Levenberg-Marquardt least-squares fitting.
	 * <p>
	 * Uses the Jacobian matrix with a quasi-Newton optimisation (that approximates the inverted Hessian). Quasi-Newton
	 * methods should avoid problems with a Hessian that cannot be inverted, e.g. in the case of round-off error
	 * introduced by vastly different magnitudes in the gradients.
	 */
	LVM_QUASI_NEWTON("Least Squares Estimator (Quasi-Newton)", "LSEqn"),
	/**
	 * Maximum Likelihood Estimator
	 * <p>
	 * Uses a Poisson noise model for the probability density function of the data.
	 */
	MLE("Maximum Likelihood Estimator", "MLE");

	private String name;
	private String shortName;

	private FitSolver(String name, String shortName)
	{
		this.name = name;
		this.shortName = shortName;
	}

	@Override
	public String toString()
	{
		return name + " (" + shortName + ")";
	}
	
	public String getName()
	{
		return name;
	}
	
	public String getShortName()
	{
		return shortName;
	}
}