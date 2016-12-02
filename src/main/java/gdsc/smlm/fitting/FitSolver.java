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
	//@formatter:off
	/**
	 * Custom Levenberg-Marquardt least-squares fitting.
	 * <p>
	 * Uses the Hessian matrix with a Newton optimisation method that requires inversion of the Hessian.
	 */
	LVM{ 
		public String getName() { return "Least Squares Estimator"; } 
		public String getShortName() { return "LSE"; }},
	/**
	 * Custom Levenberg-Marquardt least-squares fitting.
	 * <p>
	 * Uses the Hessian matrix with a Newton optimisation method that requires inversion of the Hessian.
	 * <p>
	 * Parameters can be bounded using a hard-stop limit.
	 */
	BOUNDED_LVM{ 
		public String getName() { return "Bounded Least Squares Estimator"; } 
		public String getShortName() { return "BLSE"; }},
	/**
	 * Custom Levenberg-Marquardt least-squares fitting for Poisson data using the method of Laurence & Chromy (2010) Nature Methods 7, 338-339.
	 * <p>
	 * Uses the Hessian matrix with a Newton optimisation method that requires inversion of the Hessian.
	 * <p>
	 * Parameters are bounded using a hard-stop limit to prevent negative function values being produced.
	 */
	LVM_MLE{ 
		public String getName() { return "Bounded LVM Maximum Likelihood Estimator"; } 
		public String getShortName() { return "LVM MLE"; }},
	/**
	 * Custom Levenberg-Marquardt least-squares fitting with weights. The weights require a function that provides the
	 * expected variance for each data point. Without weights the results match the LVM method.
	 * <p>
	 * Uses the Hessian matrix with a Newton optimisation method that requires inversion of the Hessian.
	 */
	LVM_WEIGHTED{ 
		public String getName() { return "Weighted Least Squares Estimator"; } 
		public String getShortName() { return "WLSE"; }},
	/**
	 * Custom Levenberg-Marquardt least-squares fitting with weights. The weights require a function that provides the
	 * expected variance for each data point. Without weights the results match the LVM method.
	 * <p>
	 * Uses the Hessian matrix with a Newton optimisation method that requires inversion of the Hessian.
	 * <p>
	 * Parameters can be bounded using a hard-stop limit.
	 */
	BOUNDED_LVM_WEIGHTED{ 
		public String getName() { return "Bounded Weighted Least Squares Estimator"; } 
		public String getShortName() { return "BWLSE"; }},
	/**
	 * Apache Commons Math Levenberg-Marquardt least-squares fitting.
	 * <p>
	 * Uses the Jacobian matrix with a quasi-Newton optimisation (that approximates the inverted Hessian). Quasi-Newton
	 * methods should avoid problems with a Hessian that cannot be inverted, e.g. in the case of round-off error
	 * introduced by vastly different magnitudes in the gradients.
	 */
	LVM_QUASI_NEWTON{ 
		public String getName() { return "Least Squares Estimator (Quasi-Newton)"; } 
		public String getShortName() { return "LSEqn"; }},
	/**
	 * Maximum Likelihood Estimator
	 * <p>
	 * Uses a Poisson noise model for the probability density function of the data.
	 */
	MLE{ 
		public String getName() { return "Maximum Likelihood Estimator"; } 
		public String getShortName() { return "MLE"; }};
	//@formatter:on

	@Override
	public String toString()
	{
		return getName() + " (" + getShortName() + ")";
	}

	/**
	 * Gets the name.
	 *
	 * @return the name
	 */
	abstract public String getName();

	/**
	 * Gets the short name.
	 *
	 * @return the short name
	 */
	abstract public String getShortName();

}