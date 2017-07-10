package gdsc.smlm.results;

import gdsc.smlm.data.config.PSFHelper;

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
 * Saves the fit results to file
 */
public abstract class SMLMFilePeakResults extends FilePeakResults
{
	private boolean showDeviations = true;
	private boolean showEndFrame = false;
	private boolean showId = false;
	protected String peakIdColumnName = "Frame";

	public SMLMFilePeakResults(String filename)
	{
		super(filename);
	}

	public SMLMFilePeakResults(String filename, boolean showDeviations)
	{
		super(filename);
		this.showDeviations = showDeviations;
	}

	public SMLMFilePeakResults(String filename, boolean showDeviations, boolean showEndFrame)
	{
		super(filename);
		this.showDeviations = showDeviations;
		this.showEndFrame = showEndFrame;
	}

	public SMLMFilePeakResults(String filename, boolean showDeviations, boolean showEndFrame, boolean showId)
	{
		super(filename);
		this.showDeviations = showDeviations;
		this.showEndFrame = showEndFrame;
		this.showId = showId;
	}

	/**
	 * Check size of the parameter arrays. The array size must match the number of parameters
	 *
	 * @param numberOfParams
	 *            the number of params
	 * @param params
	 *            the params
	 * @param paramsStdDevs
	 *            the params std devs
	 * @throws IllegalArgumentException
	 *             if the parameter arrays are the wrong size
	 */
	protected static void checkSize(int numberOfParams, float[] params, float[] paramsStdDevs)
			throws IllegalArgumentException
	{
		checkSize(numberOfParams, params);
		checkSize(numberOfParams, paramsStdDevs);
	}

	/**
	 * Check size of the parameter arrays. The array size must match the number of parameters
	 *
	 * @param numberOfParams
	 *            the number of params
	 * @param a
	 *            the a
	 * @throws IllegalArgumentException
	 *             if the parameter array is the wrong size
	 */
	protected static void checkSize(int numberOfParams, float[] a) throws IllegalArgumentException
	{
		if (a.length < numberOfParams)
			throw new IllegalArgumentException("Incorrect number of parameters " + a.length);
	}

	/**
	 * @return A line containing the file format version
	 */
	protected String getVersion()
	{
		StringBuilder sb = new StringBuilder();
		sb.append(isBinary() ? "Binary" : "Text");
		sb.append(".");
		sb.append(isShowDeviations() ? "D1" : "D0");
		sb.append(".E");
		int extended = isShowEndFrame() ? 1 : 0;
		extended += isShowId() ? 2 : 0;
		sb.append(extended);
		// Version 1 had signal and amplitude in the results. It did not have the version string.
		// .V2 = Version 2 has only signal in the results
		// .V3 = Version 3 has an improved calibration header
		sb.append(".V3");
		return sb.toString();
	}

	/**
	 * @return the name of the peak column
	 */
	public String getPeakIdColumnName()
	{
		return peakIdColumnName;
	}

	/**
	 * @param peakIdColumnName
	 *            the name of the peak column
	 */
	public void setPeakIdColumnName(String peakIdColumnName)
	{
		this.peakIdColumnName = peakIdColumnName;
	}

	/**
	 * @return True if the records contain the parameter deviations
	 */
	public boolean isShowDeviations()
	{
		return showDeviations;
	}

	/**
	 * @return True if the records contain the result end frame
	 */
	public boolean isShowEndFrame()
	{
		return showEndFrame;
	}

	/**
	 * @return True if the records contain a result Id
	 */
	public boolean isShowId()
	{
		return showId;
	}
}
