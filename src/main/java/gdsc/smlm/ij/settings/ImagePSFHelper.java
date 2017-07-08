package gdsc.smlm.ij.settings;

import gdsc.smlm.utils.XmlUtils;

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
	public static String toString(ImagePSF imagePsf)
	{
		return XmlUtils.toXML(imagePsf);
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
		Object o = XmlUtils.fromXML(string);
		if (o != null && o instanceof ImagePSF)
			return (ImagePSF) o;
		return null;
	}
}
