package gdsc.smlm.data.config;

import gdsc.smlm.data.config.ResultsProtos.ResultsFileFormat;
import gdsc.smlm.data.config.ResultsProtos.ResultsImageType;
import gdsc.smlm.data.config.ResultsProtos.ResultsSettings;

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
 * Contains helper functions for the ResultsProtos class.
 */
public class ResultsProtosHelper
{
	/** The default ResultsSettings */
	public static final ResultsSettings defaultResultsSettings;
	static
	{
		ResultsSettings.Builder builder = ResultsSettings.newBuilder();
		builder.getResultsImageSettingsBuilder().setWeighted(true);
		builder.getResultsImageSettingsBuilder().setEqualised(true);
		builder.getResultsImageSettingsBuilder().setAveragePrecision(30);
		builder.getResultsImageSettingsBuilder().setScale(1);
		builder.getResultsTableSettingsBuilder().setRoundingPrecision(4);
		builder.getResultsInMemorySettingsBuilder().setInMemory(true);
		defaultResultsSettings = builder.build();
	}	
	
	/**
	 * Gets the name.
	 *
	 * @param value
	 *            the results file format
	 * @return the name
	 */
	public static String getName(ResultsFileFormat value)
	{
		switch (value)
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
				return "Unknown";
			default:
				throw new IllegalStateException("Unknown name: " + value);
		}
	}

	/**
	 * Gets the extension.
	 *
	 * @param value
	 *            the results file format
	 * @return the extension
	 */
	public static String getExtension(ResultsFileFormat value)
	{
		switch (value)
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
	 * @param value
	 *            the results file format
	 * @return true, if is gdsc
	 */
	public static boolean isGDSC(ResultsFileFormat value)
	{
		switch (value)
		{
			case BINARY:
			case TEXT:
				return true;
			default:
				return false;
		}
	}

	/**
	 * Gets the name.
	 *
	 * @param value
	 *            the results image
	 * @return the name
	 */
	public static String getName(ResultsImageType value)
	{
		switch (value)
		{
			case DRAW_FITTED_PSF:
				return "Fitted PSF";
			case DRAW_FIT_ERROR:
				return "Fit Error";
			case DRAW_FRAME_NUMBER:
				return "Frame number";
			case DRAW_INTENSITY:
				return "Intensity";
			case DRAW_INTENSITY_AVERAGE_PRECISION:
				return "Intensity (width=av.precision)";
			case DRAW_INTENSITY_PRECISION:
				return "Intensity (width=precision)";
			case DRAW_LOCALISATIONS:
				return "Localisations";
			case DRAW_LOCALISATIONS_AVERAGE_PRECISION:
				return "Localisations (width=av.precision)";
			case DRAW_LOCALISATIONS_PRECISION:
				return "Localisations (width=precision)";
			case DRAW_NONE:
				return "None";
			case DRAW_Z_POSITION:
				return "Z position";
			case UNRECOGNIZED:
				return "Unknown";
			default:
				throw new IllegalStateException("Unknown name: " + value);
		}
	}
}
