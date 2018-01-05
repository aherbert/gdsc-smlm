package gdsc.smlm.results.count;

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
 * Interface for the evaluation of whether to stop a sequential analysis.
 */
public interface FailCounter
{
	/**
	 * Gets the description of the fail counter.
	 *
	 * @return the description (including any parameter values)
	 */
	public String getDescription();
	
	/**
	 * Called when the most recent event passed.
	 */
	public void pass();

	/**
	 * Called when the n most recent events passed.
	 * <p>
	 * This method can be used when a series of events are known to pass.
	 *
	 * @param n the n
	 */
	public void pass(int n);

	/**
	 * Called when the most recent event failed. It is expected that the result of {@link #isOK()} may change after
	 * calling this method.
	 */
	public void fail();

	/**
	 * Called when the n most recent event failed. It is expected that the result of {@link #isOK()} may change after
	 * calling this method.
	 * <p>
	 * This method can be used when a series of events are known to fail.
	 *
	 * @param n the n
	 */
	public void fail(int n);

	/**
	 * Checks if it is ok to continue the analysis. This is set to false when the analysis should stop.
	 *
	 * @return true, if is ok to continue
	 */
	public boolean isOK();

	/**
	 * Create a duplicate fail counter reset to the initialised state.
	 *
	 * @return the fail counter
	 */
	public FailCounter newCounter();

	/**
	 * Reset the counter.
	 */
	public void reset();
}