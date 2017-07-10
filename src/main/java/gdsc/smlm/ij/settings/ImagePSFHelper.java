package gdsc.smlm.ij.settings;

import java.util.HashMap;

import gdsc.smlm.data.config.PSFProtos.ImagePSF;
import gdsc.smlm.data.config.PSFProtos.ImagePSFOrBuilder;

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
 * Contain helper functions for the ImagePOSF class
 */
public class ImagePSFHelper
{
	/**
	 * Convert the ImagePSF to a string.
	 *
	 * @param imagePsf
	 *            the image psf
	 * @return the string
	 */
	public static String toString(ImagePSFOrBuilder imagePsf)
	{
		return SettingsManager.toJSON(imagePsf, SettingsManager.FLAG_JSON_WHITESPACE);
	}

	/**
	 * Get the ImagePSF from a string.
	 *
	 * @param string
	 *            the string
	 * @return the image PSF
	 */
	public static ImagePSF fromString(String string)
	{
		ImagePSF.Builder builder = ImagePSF.newBuilder();
		if (SettingsManager.fromJSON(string, builder))
			return builder.build();
		return null;
	}

	/**
	 * Creates the ImagePSF.
	 *
	 * @param centreImage
	 *            the centre image
	 * @param pixelSize
	 *            the pixel size
	 * @param pixelDepth
	 *            the pixel depth
	 * @param imageCount
	 *            the image count
	 * @param fwhm
	 *            the fwhm
	 * @return the image PSF
	 */
	public static ImagePSF create(int centreImage, double pixelSize, double pixelDepth, int imageCount, double fwhm)
	{
		return create(centreImage, pixelSize, pixelDepth, imageCount, fwhm, null);
	}

	/**
	 * Creates the ImagePSF.
	 *
	 * @param centreImage
	 *            the centre image
	 * @param pixelSize
	 *            the pixel size
	 * @param pixelDepth
	 *            the pixel depth
	 * @param imageCount
	 *            the image count
	 * @param fwhm
	 *            the fwhm
	 * @param notes
	 *            the notes
	 * @return the image PSF
	 */
	public static ImagePSF create(int centreImage, double pixelSize, double pixelDepth, int imageCount, double fwhm,
			HashMap<String, String> notes)
	{
		ImagePSF.Builder builder = ImagePSF.newBuilder();
		builder.setCentreImage(centreImage);
		builder.setPixelSize(pixelSize);
		builder.setPixelDepth(pixelDepth);
		builder.setImageCount(imageCount);
		builder.setFwhm(fwhm);
		if (notes != null)
			builder.putAllNotes(notes);
		return builder.build();
	}
}
