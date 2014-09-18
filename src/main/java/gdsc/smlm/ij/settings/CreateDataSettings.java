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
 * Contain the settings for the Create Data plugin
 */
public class CreateDataSettings
{
	public int size = 512;
	public double depth = 3000;
	public int seconds = 100;
	public int exposureTime = 100;
	public int stepsPerSecond = 10;
	public String illumination = "";
	public String backgroundImage = "";
	public String distribution = "";
	public String distributionMask = "";
	public double distributionMaskSliceDepth = 25;
	public double background = 1;
	private double emGain = 255;
	private double cameraGain = 0.1557;
	private double quantumEfficiency = 0.95;
	private double totalGain = 0;
	public double readNoise = 46;
	public int bias = 500;
	public int particles = 300;
	/**
	 * Photons are modelled using an emission rate (photon emission is constant when the fluorophore is on) sampled from
	 * a distribution. The rate (average of the distribution) is correlated to the total on time. Observations on real
	 * data indicate the correlation is around -0.35, i.e. shorter bursts are brighter.
	 */
	public int photonsPerSecond = 1000;

	/**
	 * For a simple localisation model the photons are randomly selected between photonsPerSecond and
	 * photonsPerSecondMaximum
	 */
	public int photonsPerSecondMaximum = 2000;
	/**
	 * Set to true to use a custom distribution for the random photon emission. The default is to use a gamma
	 * distribution with the defined shape parameter.
	 */
	public boolean customPhotonDistribution = false;
	public String photonDistribution = "";
	public double photonShape = 2.5;
	public double correlation = -0.35;
	public String psfModel = "";
	public String psfImageName = "";
	public double wavelength = 561;
	public double numericalAperture = 1.4;
	public double pixelPitch = 107;
	public double density = 1;
	public double diffusionRate = 0;
	public boolean compoundMolecules = false;
	public String compoundText = "";
	public boolean rotateInitialOrientation = false;
	public boolean rotateDuringSimulation = false;
	public boolean rotate2D = false;
	public boolean useGridWalk = true;
	public double fixedFraction = 0;
	public String confinement = "";
	public String confinementMask = "";
	public double confinementMaskSliceDepth = 25;
	public double confinementRadius = 10;
	public int pulseInterval = 0;
	public double pulseRatio = 100;
	/**
	 * Average t-On in milliseconds
	 */
	public double tOn = 40;
	/**
	 * Average t-Off for the short dark state in milliseconds
	 */
	public double tOffShort = 25;
	/**
	 * Average t-Off for the long dark state in milliseconds
	 */
	public double tOffLong = 631;
	/**
	 * Average number of short blinks, i.e. when transitioning from the on state to the short dark state
	 */
	public double nBlinksShort = 6.3;
	/**
	 * Average number of long blinks, i.e. transitions to the long dark state
	 */
	public double nBlinksLong = 1.8;
	/**
	 * Set to true to use a geometric distribution for the nBlinks. Default is to use a Poisson.
	 */
	public boolean nBlinksGeometricDistribution = false;

	public double minPhotons = 30;
	public double minSNRt1 = 0;
	public double minSNRtN = 0;
	public boolean saveImage = false;
	public boolean saveImageResults = false;
	public boolean saveLocalisations = false;
	public boolean saveFluorophores = false;
	public String imageFilename = "";
	public String imageResultsFilename = "";
	public String localisationsFilename = "";
	public String fluorophoresFilename = "";

	public int cellSize = 32;
	public double probabilityBinary = 0.1;
	public double minBinaryDistance = 0;
	public double maxBinaryDistance = 30;

	public boolean showHistograms = false;
	public boolean chooseHistograms = false;
	public int histogramBins = 100;
	public boolean removeOutliers = false;
	public float densityRadius = 1000;

	/**
	 * Get the total gain (the EM-gain multiplied by the camera gain). If either the gain is disabled then it will be
	 * ignored. A total gain of zero means no gain is applied.
	 * 
	 * @return
	 */
	public double getTotalGain()
	{
		return totalGain;
	}

	/**
	 * @return the emGain
	 */
	public double getEmGain()
	{
		return emGain;
	}

	/**
	 * EM-gain cannot be below 1. If so it is set to zero and disabled.
	 * 
	 * @param emGain
	 *            the emGain to set
	 */
	public void setEmGain(double emGain)
	{
		if (emGain < 1)
			emGain = 0;
		this.emGain = emGain;
		setTotalGain();
	}

	/**
	 * @return the cameraGain
	 */
	public double getCameraGain()
	{
		return cameraGain;
	}

	/**
	 * Camera gain cannot be below 0. If so it is set to zero and disabled.
	 * 
	 * @param cameraGain
	 *            the cameraGain to set
	 */
	public void setCameraGain(double cameraGain)
	{
		if (cameraGain < 0)
			cameraGain = 0;
		this.cameraGain = cameraGain;
		setTotalGain();
	}

	/**
	 * @return the quantumEfficiency
	 */
	public double getQuantumEfficiency()
	{
		return quantumEfficiency;
	}

	/**
	 * QE cannot be below 0 or above 1. If so it is set to one and disabled.
	 * 
	 * @param quantumEfficiency
	 *            the quantumEfficiency to set
	 */
	public void setQuantumEfficiency(double quantumEfficiency)
	{
		if (quantumEfficiency < 0 || quantumEfficiency > 1)
			quantumEfficiency = 1;
		this.quantumEfficiency = quantumEfficiency;
		setTotalGain();
	}

	/**
	 * Ensure that the internal state of the object is initialised. This is used after deserialisation since some state
	 * is not saved but restored from other property values.
	 */
	public void initialiseState()
	{
		setTotalGain();
	}

	/**
	 * Set the total gain (the EM-gain multiplied by the camera gain multiplied by the quantum efficiency). If either
	 * the gain is disabled then it will be ignored. A total gain of zero means no gain is applied.
	 */
	private void setTotalGain()
	{
		totalGain = 0;
		if (cameraGain > 0)
		{
			totalGain = (emGain > 0) ? emGain * cameraGain : cameraGain;
		}
		else if (emGain > 0)
			totalGain = emGain;
		totalGain *= quantumEfficiency;
	}
}
