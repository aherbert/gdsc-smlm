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

/**
 * Contain the settings for the gdsc.fitting package
 * 
 * @deprecated The settings have been moved to the gdsc.smlm.data.config package.
 */
@Deprecated
public class GlobalSettings
{
	private String notes = null;
	private FitEngineConfiguration fitEngineConfiguration = null;
	private PSFEstimatorSettings psfEstimatorSettings = null;
	private PSFCalculatorSettings psfCalculatorSettings = null;
	private FilterSettings filterSettings = null;
	private ClusteringSettings clusteringSettings = null;
	private CreateDataSettings createDataSettings = null;
	private OPTICSSettings opticsSettings = null;

	/**
	 * @return the notes
	 */
	public String getNotes()
	{
		return notes;
	}

	/**
	 * @param notes
	 *            the notes to set
	 */
	public void setNotes(String notes)
	{
		this.notes = notes;
	}

	/**
	 * @return the fitEngineConfiguration
	 */
	public FitEngineConfiguration getFitEngineConfiguration()
	{
		if (fitEngineConfiguration == null)
			fitEngineConfiguration = new FitEngineConfiguration();
		return fitEngineConfiguration;
	}

	/**
	 * @return the psfEstimatorSettings
	 */
	public PSFEstimatorSettings getPsfEstimatorSettings()
	{
		if (psfEstimatorSettings == null)
			psfEstimatorSettings = new PSFEstimatorSettings();
		return psfEstimatorSettings;
	}

	/**
	 * @return the psfCalculatorSettings
	 */
	public PSFCalculatorSettings getPsfCalculatorSettings()
	{
		if (psfCalculatorSettings == null)
			psfCalculatorSettings = new PSFCalculatorSettings();
		return psfCalculatorSettings;
	}

	/**
	 * @return the filterSettings
	 */
	public FilterSettings getFilterSettings()
	{
		if (filterSettings == null)
			filterSettings = new FilterSettings();
		return filterSettings;
	}

	/**
	 * @return the clusteringSettings
	 */
	public ClusteringSettings getClusteringSettings()
	{
		if (clusteringSettings == null)
			clusteringSettings = new ClusteringSettings();
		return clusteringSettings;
	}

	/**
	 * @return the createDataSettings
	 */
	public CreateDataSettings getCreateDataSettings()
	{
		if (createDataSettings == null)
			createDataSettings = new CreateDataSettings();
		return createDataSettings;
	}

	/**
	 * @return the opticsSettings
	 */
	public OPTICSSettings getOPTICSSettings()
	{
		if (opticsSettings == null)
			opticsSettings = new OPTICSSettings();
		return opticsSettings;
	}

	/**
	 * @param config
	 */
	public void setFitEngineConfiguration(FitEngineConfiguration config)
	{
		this.fitEngineConfiguration = config;
	}

	/**
	 * Check the setting is not currently null. If the setting is null then a call to the get() method will initialise a
	 * default object.
	 * 
	 * @return true if the setting is not null
	 */
	public boolean isFitEngineConfiguration()
	{
		return (fitEngineConfiguration != null);
	}

	/**
	 * Check the setting is not currently null. If the setting is null then a call to the get() method will initialise a
	 * default object.
	 * 
	 * @return true if the setting is not null
	 */
	public boolean isPsfEstimatorSettings()
	{
		return (psfEstimatorSettings != null);
	}

	/**
	 * Check the setting is not currently null. If the setting is null then a call to the get() method will initialise a
	 * default object.
	 * 
	 * @return true if the setting is not null
	 */
	public boolean isPsfCalculatorSettings()
	{
		return (psfCalculatorSettings != null);
	}

	/**
	 * Check the setting is not currently null. If the setting is null then a call to the get() method will initialise a
	 * default object.
	 * 
	 * @return true if the setting is not null
	 */
	public boolean isFilterSettings()
	{
		return (filterSettings != null);
	}

	/**
	 * Check the setting is not currently null. If the setting is null then a call to the get() method will initialise a
	 * default object.
	 * 
	 * @return true if the setting is not null
	 */
	public boolean isClusteringSettings()
	{
		return (clusteringSettings != null);
	}

	/**
	 * Check the setting is not currently null. If the setting is null then a call to the get() method will initialise a
	 * default object.
	 * 
	 * @return true if the setting is not null
	 */
	public boolean isCreateDataSettings()
	{
		return (createDataSettings != null);
	}

	/**
	 * Check the setting is not currently null. If the setting is null then a call to the get() method will initialise a
	 * default object.
	 * 
	 * @return true if the setting is not null
	 */
	public boolean isOPTICSSettings()
	{
		return (opticsSettings != null);
	}
}
