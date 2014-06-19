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

import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.fitting.FitConfiguration;

/**
 * Contain the configuration for a single run of the batch fitting plugin
 */
public class BatchRun
{
	public String image;
	public FitEngineConfiguration fitEngineConfiguration = null;
	
	public BatchRun()
	{
		fitEngineConfiguration = new FitEngineConfiguration(new FitConfiguration());
	}
	public BatchRun(String image, FitEngineConfiguration fitEngineConfiguration)
	{
		this.image = image;
		this.fitEngineConfiguration = fitEngineConfiguration;
	}
}
