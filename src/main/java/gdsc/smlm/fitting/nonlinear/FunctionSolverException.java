package gdsc.smlm.fitting.nonlinear;

import gdsc.smlm.fitting.FitStatus;

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
 * Exception to throw if a function solver failed.
 */
public class FunctionSolverException extends RuntimeException
{
	/** The Constant serialVersionUID. */
	private static final long serialVersionUID = -5131234527135746186L;

	/** The fit status. This indicates why the solver failed. */
	public final FitStatus fitStatus;

	/**
	 * Instantiates a new function solver exception.
	 *
	 * @param fitStatus the fit status
	 */
	public FunctionSolverException(FitStatus fitStatus)
	{
		super();
		this.fitStatus = fitStatus;
	}

	/**
	 * Instantiates a new function solver exception.
	 *
	 * @param fitStatus the fit status
	 * @param message the message
	 */
	public FunctionSolverException(FitStatus fitStatus, String message)
	{
		super(message);
		this.fitStatus = fitStatus;
	}

	/**
	 * Instantiates a new function solver exception.
	 *
	 * @param fitStatus the fit status
	 * @param message the message
	 * @param cause the cause
	 */
	public FunctionSolverException(FitStatus fitStatus, String message, Throwable cause)
	{
		super(message, cause);
		this.fitStatus = fitStatus;
	}

	/**
	 * Instantiates a new function solver exception.
	 *
	 * @param fitStatus the fit status
	 * @param cause the cause
	 */
	public FunctionSolverException(FitStatus fitStatus, Throwable cause)
	{
		super(cause);
		this.fitStatus = fitStatus;
	}
}
