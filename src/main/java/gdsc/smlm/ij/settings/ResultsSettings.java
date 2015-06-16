package gdsc.smlm.ij.settings;

import gdsc.smlm.ij.results.ResultsImage;
import gdsc.smlm.ij.results.ResultsTable;

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
 * Contain the settings for the Results
 */
public class ResultsSettings implements Cloneable
{
	public boolean logProgress = false;
	public boolean showDeviations = false;
	private ResultsImage resultsImage = ResultsImage.LOCALISATIONS;
	public boolean weightedImage = true;
	public boolean equalisedImage = true;
	public double precision = 0.3;
	public double imageScale = 8;
	public int imageRollingWindow = 0;
	public String resultsDirectory = null;
	public String resultsFilename = null;
	public boolean binaryResults = false;
	public boolean resultsInMemory = false;
	private ResultsTable resultsTable = ResultsTable.NONE;

	/**
	 * @return the resultsImage
	 */
	public ResultsImage getResultsImage()
	{
		return resultsImage;
	}

	/**
	 * @param resultsImage
	 *            the resultsImage to set
	 */
	public void setResultsImage(ResultsImage resultsImage)
	{
		this.resultsImage = resultsImage;
	}

	/**
	 * @param resultsImage
	 *            the resultsImage to set
	 */
	public void setResultsImage(int resultsImage)
	{
		if (resultsImage >= 0 && resultsImage < ResultsImage.values().length)
		{
			setResultsImage(ResultsImage.values()[resultsImage]);
		}
	}

	/**
	 * @return the resultsTable
	 */
	public ResultsTable getResultsTable()
	{
		return resultsTable;
	}

	/**
	 * @param resultsTable
	 *            the resultsTable to set
	 */
	public void setResultsTable(ResultsTable resultsTable)
	{
		this.resultsTable = resultsTable;
	}

	/**
	 * @param resultsTable
	 *            the resultsImage to set
	 */
	public void setResultsTable(int resultsTable)
	{
		if (resultsTable >= 0 && resultsTable < ResultsTable.values().length)
		{
			setResultsTable(ResultsTable.values()[resultsTable]);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public ResultsSettings clone()
	{
		ResultsSettings c;
		try
		{
			c = (ResultsSettings) super.clone();
		}
		catch (CloneNotSupportedException e)
		{
			return null;
		}
		return c;
	}
}
