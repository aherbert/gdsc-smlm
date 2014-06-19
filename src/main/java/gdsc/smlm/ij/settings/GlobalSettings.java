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
import gdsc.smlm.results.Calibration;

/**
 * Contain the settings for the gdsc.fitting package
 */
public class GlobalSettings
{
	private Calibration calibration = null;
	private FitEngineConfiguration fitEngineConfiguration = null;
	private PSFEstimatorSettings psfEstimatorSettings = null;
	private PSFCalculatorSettings psfCalculatorSettings = null;
	private ResultsSettings resultsSettings = null;
	private FilterSettings filterSettings = null;
	private ClusteringSettings clusteringSettings = null; 
	private CreateDataSettings createDataSettings = null;

	/**
	 * @return the calibration
	 */
	public Calibration getCalibration()
	{
		if (calibration == null)
			calibration = new Calibration();
		return calibration;
	}

	/**
	 * @return the fitEngineConfiguration
	 */
	public FitEngineConfiguration getFitEngineConfiguration()
	{
		if (fitEngineConfiguration == null)
			fitEngineConfiguration = new FitEngineConfiguration(new FitConfiguration());
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
	 * @return the resultsSettings
	 */
	public ResultsSettings getResultsSettings()
	{
		if (resultsSettings == null)
			resultsSettings = new ResultsSettings();
		return resultsSettings;
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
	 * @param config
	 */
	public void setFitEngineConfiguration(FitEngineConfiguration config)
	{
		this.fitEngineConfiguration = config;
	}

	/**
	 * @param resultsSettings
	 */
	public void setResultsSettings(ResultsSettings resultsSettings)
	{
		this.resultsSettings = resultsSettings;
	}
}
