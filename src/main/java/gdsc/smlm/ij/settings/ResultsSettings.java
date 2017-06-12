package gdsc.smlm.ij.settings;

import gdsc.smlm.data.config.SMLMSettings.AngleUnit;
import gdsc.smlm.data.config.SMLMSettings.DistanceUnit;
import gdsc.smlm.data.config.SMLMSettings.IntensityUnit;
import gdsc.smlm.ij.results.ResultsFileFormat;
import gdsc.smlm.ij.results.ResultsImage;

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
	private ResultsImage resultsImage;
	public boolean weightedImage = true;
	public boolean equalisedImage = true;
	public double precision = 0.3;
	public double imageScale = 8;
	public int imageRollingWindow = 0;
	public String resultsDirectory = null;

	public String resultsFilename = null;
	private ResultsFileFormat resultsFileFormat;
	private DistanceUnit fileDistanceUnit;
	private IntensityUnit fileIntensityUnit;
	private AngleUnit fileAngleUnit;
	public boolean fileComputePrecision = false;

	public boolean resultsInMemory = true;

	public boolean showResultsTable = false;
	private DistanceUnit tableDistanceUnit;
	private IntensityUnit tableIntensityUnit;
	private AngleUnit tableAngleUnit;
	public boolean tableComputePrecision = false;
	public int tableRoundingPrecision = 4;

	public ResultsSettings()
	{
		initialiseState();
	}

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
	 * @return the resultsFileFormat
	 */
	public ResultsFileFormat getResultsFileFormat()
	{
		return resultsFileFormat;
	}

	/**
	 * @param resultsFileFormat
	 *            the resultsFileFormat to set
	 */
	public void setResultsFileFormat(ResultsFileFormat resultsFileFormat)
	{
		this.resultsFileFormat = resultsFileFormat;
	}

	/**
	 * @param resultsFileFormat
	 *            the resultsFileFormat to set
	 */
	public void setResultsFileFormat(int resultsFileFormat)
	{
		if (resultsFileFormat >= 0 && resultsFileFormat < ResultsFileFormat.values().length)
		{
			setResultsFileFormat(ResultsFileFormat.values()[resultsFileFormat]);
		}
	}

	/**
	 * @return the fileDistanceUnit
	 */
	public DistanceUnit getFileDistanceUnit()
	{
		return fileDistanceUnit;
	}

	/**
	 * @param fileDistanceUnit
	 *            the fileDistanceUnit to set
	 */
	public void setFileDistanceUnit(DistanceUnit fileDistanceUnit)
	{
		this.fileDistanceUnit = fileDistanceUnit;
	}

	/**
	 * @param fileDistanceUnit
	 *            the fileDistanceUnit to set
	 */
	public void setFileDistanceUnit(int fileDistanceUnit)
	{
		if (fileDistanceUnit >= 0 && fileDistanceUnit < DistanceUnit.values().length)
		{
			setFileDistanceUnit(DistanceUnit.values()[fileDistanceUnit]);
		}
	}

	/**
	 * @return the fileIntensityUnit
	 */
	public IntensityUnit getFileIntensityUnit()
	{
		return fileIntensityUnit;
	}

	/**
	 * @param fileIntensityUnit
	 *            the fileIntensityUnit to set
	 */
	public void setFileIntensityUnit(IntensityUnit fileIntensityUnit)
	{
		this.fileIntensityUnit = fileIntensityUnit;
	}

	/**
	 * @param fileIntensityUnit
	 *            the fileIntensityUnit to set
	 */
	public void setFileIntensityUnit(int fileIntensityUnit)
	{
		if (fileIntensityUnit >= 0 && fileIntensityUnit < IntensityUnit.values().length)
		{
			setFileIntensityUnit(IntensityUnit.values()[fileIntensityUnit]);
		}
	}

	/**
	 * @return the fileAngleUnit
	 */
	public AngleUnit getFileAngleUnit()
	{
		return fileAngleUnit;
	}

	/**
	 * @param fileAngleUnit
	 *            the fileAngleUnit to set
	 */
	public void setFileAngleUnit(AngleUnit fileAngleUnit)
	{
		this.fileAngleUnit = fileAngleUnit;
	}

	/**
	 * @param fileAngleUnit
	 *            the fileAngleUnit to set
	 */
	public void setFileAngleUnit(int fileAngleUnit)
	{
		if (fileAngleUnit >= 0 && fileAngleUnit < AngleUnit.values().length)
		{
			setFileAngleUnit(AngleUnit.values()[fileAngleUnit]);
		}
	}

	/**
	 * @return the tableDistanceUnit
	 */
	public DistanceUnit getTableDistanceUnit()
	{
		return tableDistanceUnit;
	}

	/**
	 * @param tableDistanceUnit
	 *            the tableDistanceUnit to set
	 */
	public void setTableDistanceUnit(DistanceUnit tableDistanceUnit)
	{
		this.tableDistanceUnit = tableDistanceUnit;
	}

	/**
	 * @param tableDistanceUnit
	 *            the tableDistanceUnit to set
	 */
	public void setTableDistanceUnit(int tableDistanceUnit)
	{
		if (tableDistanceUnit >= 0 && tableDistanceUnit < DistanceUnit.values().length)
		{
			setTableDistanceUnit(DistanceUnit.values()[tableDistanceUnit]);
		}
	}

	/**
	 * @return the tableIntensityUnit
	 */
	public IntensityUnit getTableIntensityUnit()
	{
		return tableIntensityUnit;
	}

	/**
	 * @param tableIntensityUnit
	 *            the tableIntensityUnit to set
	 */
	public void setTableIntensityUnit(IntensityUnit tableIntensityUnit)
	{
		this.tableIntensityUnit = tableIntensityUnit;
	}

	/**
	 * @param tableIntensityUnit
	 *            the tableIntensityUnit to set
	 */
	public void setTableIntensityUnit(int tableIntensityUnit)
	{
		if (tableIntensityUnit >= 0 && tableIntensityUnit < IntensityUnit.values().length)
		{
			setTableIntensityUnit(IntensityUnit.values()[tableIntensityUnit]);
		}
	}

	/**
	 * @return the tableAngleUnit
	 */
	public AngleUnit getTableAngleUnit()
	{
		return tableAngleUnit;
	}

	/**
	 * @param tableAngleUnit
	 *            the tableAngleUnit to set
	 */
	public void setTableAngleUnit(AngleUnit tableAngleUnit)
	{
		this.tableAngleUnit = tableAngleUnit;
	}

	/**
	 * @param tableAngleUnit
	 *            the tableAngleUnit to set
	 */
	public void setTableAngleUnit(int tableAngleUnit)
	{
		if (tableAngleUnit >= 0 && tableAngleUnit < AngleUnit.values().length)
		{
			setTableAngleUnit(AngleUnit.values()[tableAngleUnit]);
		}
	}

	/**
	 * Ensure that the internal state of the object is initialised. This is used after deserialisation since some state
	 * is not saved but restored from other property values.
	 */
	public void initialiseState()
	{
		if (resultsImage == null)
			resultsImage = ResultsImage.LOCALISATIONS;
		if (resultsFileFormat == null)
			resultsFileFormat = ResultsFileFormat.GDSC_TEXT;
		if (fileDistanceUnit == null)
			fileDistanceUnit = DistanceUnit.PIXEL; // Legacy default
		if (fileIntensityUnit == null)
			fileIntensityUnit = IntensityUnit.PHOTON;
		if (fileAngleUnit == null)
			fileAngleUnit = AngleUnit.RADIAN;
		if (tableDistanceUnit == null)
			tableDistanceUnit = DistanceUnit.PIXEL; // Legacy default
		if (tableIntensityUnit == null)
			tableIntensityUnit = IntensityUnit.PHOTON;
		if (tableAngleUnit == null)
			tableAngleUnit = AngleUnit.DEGREE;
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
