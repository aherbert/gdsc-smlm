package gdsc.smlm.ij.settings;

import java.util.HashMap;

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
 * Contain the settings for the PSF Creator plugin
 * @deprecated This will be replaced with a proto object
 */
@Deprecated
public class ImagePSF
{
	private int zCentre = 1;
	private double nmPerPixel = 100;
	private double nmPerSlice = 20;
	private int nImages = 1;
	private double fwhm = 1;
	private HashMap<String, String> note;
	private HashMap<Integer,Offset> offset;

	public ImagePSF(int zCentre, double nmPerPixel, double nmPerSlice, int nImages, double fwhm)
	{
		this(zCentre, nmPerPixel, nmPerSlice, nImages, fwhm, null);
	}

	/**
	 * Notes should be added as "key: value\n"
	 * 
	 * @param zCentre
	 * @param nmPerPixel
	 * @param nmPerSlice
	 * @param nImages
	 * @param fwhm
	 * @param note
	 *            String of notes using the "key: value\n" format
	 */
	public ImagePSF(int zCentre, double nmPerPixel, double nmPerSlice, int nImages, double fwhm,
			HashMap<String, String> note)
	{
		this.setCentreImage(zCentre);
		this.setPixelSize(nmPerPixel);
		this.setPixelDepth(nmPerSlice);
		this.setImageCount(nImages);
		this.setFwhm(fwhm);
		this.note = (note == null) ? new HashMap<String, String>() : note;
	}

	/**
	 * Add a note to the note string. Notes should be added as "key: value\n". An existing note with the same key will
	 * be replaced.
	 * 
	 * @param key
	 * @param value
	 */
	public void addNote(String key, String value)
	{
		note.put(key, value);
	}

	public int getCentreImage()
	{
		return zCentre;
	}

	public void setCentreImage(int zCentre)
	{
		this.zCentre = zCentre;
	}

	public double getPixelSize()
	{
		return nmPerPixel;
	}

	public void setPixelSize(double nmPerPixel)
	{
		this.nmPerPixel = nmPerPixel;
	}

	public double getPixelDepth()
	{
		return nmPerSlice;
	}

	public void setPixelDepth(double nmPerSlice)
	{
		this.nmPerSlice = nmPerSlice;
	}

	public int getImageCount()
	{
		return nImages;
	}

	public void setImageCount(int nImages)
	{
		this.nImages = nImages;
	}

	public double getFwhm()
	{
		return fwhm;
	}

	public void setFwhm(double fwhm)
	{
		this.fwhm = fwhm;
	}

	public HashMap<Integer,Offset> getOffset()
	{
		return offset;
	}

	public void setOffset(HashMap<Integer,Offset> offset)
	{
		this.offset = offset;
	}
}
