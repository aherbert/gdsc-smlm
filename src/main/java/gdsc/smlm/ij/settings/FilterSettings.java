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
 * Contain the settings for the result filtering
 */
public class FilterSettings
{
	public float maxDrift = 0;
	public float minSignal = 0;
	public float minSNR = 0;
	public double maxPrecision = 0;
	public float maxWidth = 0;
	public float minWidth = 0;
	public String maskTitle = "";
	public String freeFilter = "";
	public String filterTemplate = "";
	public String filterAnalysisDirectory = "";
	public String filterSetFilename = "";
}
