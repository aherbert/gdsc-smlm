package gdsc.smlm.ij.plugins;

/*----------------------------------------------------------------------------- 
 * GDSC Plugins for ImageJ
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
 * A worker to update results based on new settings.
 * 
 * @author Alex Herbert
 */
public abstract class WorkflowWorker<S, R>
{
	/**
	 * Compare the settings and return false if any settings that the work depends on have changed
	 * <p>
	 * Both objects will not be null.
	 *
	 * @param current
	 *            the current
	 * @param previous
	 *            the previous
	 * @return true if settings have changed
	 */
	public abstract boolean equalSettings(S current, S previous);

	/**
	 * Compare the results and return false if any results that the work depends on have changed
	 * <p>
	 * Either object could be null (if no results have yet been generated for this work).
	 *
	 * @param current
	 *            the current
	 * @param previous
	 *            the previous
	 * @return true if results have changed
	 */
	public abstract boolean equalResults(R current, R previous);

	/**
	 * Creates the results.
	 *
	 * @param settings
	 *            the settings
	 * @param results
	 *            the results
	 * @return the results
	 */
	public abstract R createResults(S settings, R results);

	/**
	 * Called when there are new results in the current work. This can be used to reset the worker before
	 * {@link #createResult(Work)} is called.
	 */
	protected void newResults()
	{
	}
}
