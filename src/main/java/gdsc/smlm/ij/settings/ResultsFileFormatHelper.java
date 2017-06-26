package gdsc.smlm.ij.settings;

import gdsc.smlm.data.config.ResultsConfig.ResultsFileFormat;

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

public class ResultsFileFormatHelper
{
	/**
	 * Gets the name.
	 *
	 * @param resultsFileFormat
	 *            the results file format
	 * @return the name
	 */
	public static String getName(ResultsFileFormat resultsFileFormat)
	{
		switch (resultsFileFormat)
		{
			case BINARY:
				return "Binary";
			case FILE_NONE:
				return "None";
			case MALK:
				return "MALK (Molecular Accuracy Localisation Keep)";
			case TEXT:
				return "Text";
			case TSF:
				return "TSF (Tagged Spot File)";
			case UNRECOGNIZED:
			default:
				return "Unknown";
		}
	}

	/**
	 * Gets the extension.
	 *
	 * @param resultsFileFormat
	 *            the results file format
	 * @return the extension
	 */
	public static String getExtension(ResultsFileFormat resultsFileFormat)
	{
		switch (resultsFileFormat)
		{
			case BINARY:
				return "bin";
			case FILE_NONE:
				return "";
			case MALK:
				return "txt";
			case TEXT:
				return "xls";
			case TSF:
				return "tsf";
			case UNRECOGNIZED:
			default:
				return "";
		}
	}

	/**
	 * Checks if is gdsc.
	 *
	 * @param resultsFileFormat
	 *            the results file format
	 * @return true, if is gdsc
	 */
	public static boolean isGDSC(ResultsFileFormat resultsFileFormat)
	{
		switch (resultsFileFormat)
		{
			case BINARY:
			case TEXT:
				return true;
			default:
				return false;
		}
	}
}
