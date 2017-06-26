package gdsc.smlm.ij.settings;

import gdsc.smlm.data.config.ResultsConfig.ResultsImageType;

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

public class ResultsImageTypeHelper
{

	/**
	 * Gets the name.
	 *
	 * @param resultsImage
	 *            the results image
	 * @return the name
	 */
	public static String getName(ResultsImageType resultsImage)
	{
		switch (resultsImage)
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
			case UNRECOGNIZED:
			default:
				return "Unknown";
		}
	}
}
