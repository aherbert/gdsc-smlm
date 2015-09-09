package gdsc.smlm.ij.settings;

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
 */
public class PSFSettings
{
	public int zCentre = 1;
	public double nmPerPixel = 100;
	public double nmPerSlice = 20;
	public int nImages = 1;
	public double fwhm = 1;
	/**
	 * Notes should be added as "key: value\n"
	 */
	public String note = null;
	public PSFOffset[] offset;

	public PSFSettings(int zCentre, double nmPerPixel, double nmPerSlice, int nImages, double fwhm)
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
	public PSFSettings(int zCentre, double nmPerPixel, double nmPerSlice, int nImages, double fwhm, String note)
	{
		this.zCentre = zCentre;
		this.nmPerPixel = nmPerPixel;
		this.nmPerSlice = nmPerSlice;
		this.nImages = nImages;
		this.fwhm = fwhm;
		this.note = note;
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
		StringBuilder sb = new StringBuilder();
		boolean found = false;
		if (note != null)
		{
			String[] lines = note.split("\n");
			for (int i = 0; i < lines.length; i++)
			{
				if (lines[i].startsWith(key))
				{
					// Only update the first matching key, others are removed
					if (!found)
						sb.append(key).append(": ").append(value).append("\n");
					found = true;
				}
				else
					sb.append(lines[i]).append("\n");
			}
		}
		if (!found)
			sb.append(key).append(": ").append(value).append("\n");
		note = sb.toString();
	}
}
