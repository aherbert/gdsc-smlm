package gdsc.smlm.ij.results;

import java.awt.Rectangle;

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

public class ImagePeakResultsFactory
{
	/**
	 * Create a PeakResults image using the specified parameters
	 * 
	 * @param resultsImage
	 *            The type of image
	 * @param weighted
	 *            Flag to indicate is the values should be bilinearly weighted on surrounding 4 pixels (applied when
	 *            plotting localisations at a single point)
	 * @param equalised
	 *            Flag to indicate if the image should have histogram equalisation applied
	 * @param title
	 *            The title of the image
	 * @param bounds
	 *            Define the bounding rectangle of the result coordinates
	 * @param nmPerPixel
	 *            The results scale in nanometers per pixel
	 * @param gain
	 *            The results gain
	 * @param imageScale
	 *            Define the scale of the image relative to the bounding rectangle
	 * @param precision
	 *            For average precision plots this parameter specifies the fixed width of the PSF (in nm).
	 *            If less than zero then defaults to the nmPerPixel value.
	 * @param mode
	 *            The mode for showing consecutive results in the same pixel location
	 * @return The PeakResults image
	 */
	public static IJImagePeakResults createPeakResultsImage(ResultsImage resultsImage, boolean weighted,
			boolean equalised, String title, Rectangle bounds, double nmPerPixel, double gain, double imageScale,
			double precision, ResultsMode mode)
	{
		IJImagePeakResults image;
		switch (resultsImage)
		{
			case PSF:
			case SIGNAL_PRECISION:
			case LOCALISATIONS_PRECISION:
			case SIGNAL_AV_PRECISION:
			case LOCALISATIONS_AV_PRECISION:
				// Special case for full PSF image
				PSFImagePeakResults image2 = new PSFImagePeakResults(title, bounds, (float) imageScale);
				if (resultsImage == ResultsImage.SIGNAL_AV_PRECISION ||
						resultsImage == ResultsImage.LOCALISATIONS_AV_PRECISION)
				{
					// Fixed width display (in pixels)
					image2.setWidth((float) (((precision <= 0) ? nmPerPixel : precision) / nmPerPixel));
				}
				else if (resultsImage == ResultsImage.SIGNAL_PRECISION ||
						resultsImage == ResultsImage.LOCALISATIONS_PRECISION)
				{
					image2.setCalculatedPrecision(true);
				}
				image = image2;
				break;

			default:
				image = new IJImagePeakResults(title, bounds, (float) imageScale);
		}
		int flags = 0;

		switch (resultsImage)
		{
			case SIGNAL_INTENSITY:
			case SIGNAL_PRECISION:
			case PSF:
				flags |= IJImagePeakResults.DISPLAY_SIGNAL;
				break;
			case FRAME_NUMBER:
				flags |= IJImagePeakResults.DISPLAY_PEAK;
				break;
			case ERROR:
				flags |= IJImagePeakResults.DISPLAY_ERROR;
				break;
			default:
				// Nothing to do for the other cases
				break;
		}

		switch (mode)
		{
			case MAX:
				flags |= IJImagePeakResults.DISPLAY_MAX;
				break;
			case REPLACE:
				flags |= IJImagePeakResults.DISPLAY_REPLACE;
				break;
			default:
				// Nothing to do for the other cases
				break;
		}

		if (weighted)
		{
			flags |= IJImagePeakResults.DISPLAY_WEIGHTED;
		}
		if (equalised)
		{
			flags |= IJImagePeakResults.DISPLAY_EQUALIZED;
		}
		image.setDisplayFlags(flags);
		image.setCalibration(nmPerPixel, gain);
		return image;
	}
}
