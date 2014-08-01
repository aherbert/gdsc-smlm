package gdsc.smlm.results;

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
 * Track the progress of processing results
 */
public interface TrackProgress
{
	/**
	 * Specify progress as a fraction
	 * 
	 * @param fraction
	 */
	public void progress(double fraction);

	/**
	 * Specify progress as the position relative to the total
	 * 
	 * @param position
	 * @param total
	 */
	public void progress(long position, long total);

	/**
	 * Logs a message on the progress
	 * 
	 * @param format
	 * @param args
	 */
	public void log(String format, Object... args);

	/**
	 * Sets the status on the progress
	 * 
	 * @param format
	 * @param args
	 */
	public void status(String format, Object... args);

	/**
	 * Return true if the processing of results is active
	 * 
	 * @return
	 */
	public boolean stop();
}
