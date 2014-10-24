package gdsc.smlm.engine;

import java.util.List;

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
 * Specifies additional parameters for the job.
 */
public class FitParameters
{
	public enum FitTask
	{
		PSF_FITTING("PSF Fitting"),
		MAXIMA_IDENITIFICATION("Maxima Identification");

		private String name;

		private FitTask(String name)
		{
			this.name = name;
		}

		@Override
		public String toString()
		{
			return name;
		}
	}	
	
	private float[] offset = null;

	/**
	 * The noise for the image data
	 */
	public float noise = Float.NaN;
	/**
	 * The maxima to fit within the data
	 */
	public int[] maxIndices = null;
	/**
	 * The background for the image data
	 */
	public float background = Float.NaN;
	/**
	 * Only maxima within the distance threshold to these coordinates will be included in the results.
	 */
	public List<float[]> filter = null;
	/**
	 * The distance threshold to use when checking if fitted peaks match the desired results. 
	 */
	public float distanceThreshold = 1;
	/**
	 * The task to perform
	 */
	public FitTask fitTask = FitTask.PSF_FITTING;
	/**
	 * The frame acquisition end time for the input data. Used when data represents multiple frames.  
	 */
	public int endT = -1;

	/**
	 * The offset to apply to all fitted results
	 */
	public float[] getOffset()
	{
		return offset;
	}

	/**
	 * Set the X,Y offset. Must be an array of length 2 otherwise the offset is reset to null
	 * 
	 * @param offset
	 */
	public void setOffset(float[] offset)
	{
		if (offset != null)
		{
			if (offset.length != 2)
				offset = null;
		}
		this.offset = offset;
	}
}
