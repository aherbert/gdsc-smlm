package gdsc.smlm.data.config;

import java.util.ArrayList;
import java.util.List;

import gdsc.smlm.data.config.SMLMSettings.PSF;
import gdsc.smlm.data.config.SMLMSettings.PSFParameter;
import gdsc.smlm.data.config.SMLMSettings.PSFParameterUnit;
import gdsc.smlm.data.config.SMLMSettings.PSFType;

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
 * Contains helper functions for the PSF class.
 */
public class PSFHelper
{
	private final PSF.Builder psfBuilder;

	/**
	 * Instantiates a new psf helper.
	 *
	 * @param psf
	 *            the psf
	 * @throws IllegalArgumentException
	 *             if the psf is null
	 */
	public PSFHelper(PSF psf) throws IllegalArgumentException
	{
		if (psf == null)
			throw new IllegalArgumentException("PSF is null");
		psfBuilder = psf.toBuilder();
	}

	/**
	 * Gets the psf.
	 *
	 * @return the psf
	 */
	public PSF getPSF()
	{
		return psfBuilder.build();
	}

	/**
	 * Checks if is a Gaussian 2D PSF.
	 *
	 * @param psf
	 *            the psf
	 * @return true, if is a Gaussian 2D PSF
	 * @throws ConfigurationException
	 *             if the psf is null
	 */
	public static boolean isGaussian2D(PSF psf) throws ConfigurationException
	{
		if (psf == null)
			throw new ConfigurationException("psf is null");
		switch (psf.getPsfType())
		{
			case OneAxisGaussian2D:
			case AstigmaticGaussian2D:
			case TwoAxisAndThetaGaussian2D:
			case TwoAxisGaussian2D:
				return true;
			case Custom:
			case UNRECOGNIZED:
			default:
				break;
		}
		return false;
	}

	/**
	 * Gets the Gaussian 2D x-width and y-width indices for the additional parameters.
	 *
	 * @param psf
	 *            the psf
	 * @return the Gaussian 2D x-width and y-width indices for the additional parameters.
	 * @throws ConfigurationException
	 *             if the psf is null, or not a Gaussian 2D function
	 */
	public static int[] getGaussian2DWxWyIndices(PSF psf) throws ConfigurationException
	{
		if (psf == null)
			throw new ConfigurationException("psf is null");
		switch (psf.getPsfType())
		{
			case OneAxisGaussian2D:
				return new int[] { 0, 0 };
			case AstigmaticGaussian2D:
			case TwoAxisAndThetaGaussian2D:
			case TwoAxisGaussian2D:
				return new int[] { 0, 1 };
			case Custom:
			case UNRECOGNIZED:
			default:
				break;
		}
		throw new ConfigurationException("psf is not Gaussian2D");
	}

	private static final ArrayList<PSFParameter> sxParameters, sxsyParameters, sxsyaParameters;
	static
	{
		PSFParameter.Builder builder = PSFParameter.newBuilder();

		builder.setName("S");
		builder.setUnit(PSFParameterUnit.DISTANCE);
		sxParameters = new ArrayList<PSFParameter>(1);
		sxParameters.add(builder.build());

		sxsyParameters = new ArrayList<PSFParameter>(2);
		builder.setName("Sx");
		sxsyParameters.add(builder.build());
		builder.setName("Sy");
		sxsyParameters.add(builder.build());

		sxsyaParameters = new ArrayList<PSFParameter>(sxsyParameters);
		builder.setName("Angle");
		builder.setUnit(PSFParameterUnit.ANGLE);
		sxsyaParameters.add(builder.build());
		sxsyaParameters.trimToSize();
	}

	/**
	 * Gets the parameters for the PSF. If the PSF is a standard Gaussian2D function then the configured parameter list
	 * is only used if it has the correct size and parameter types, otherwise a default list is returned.
	 *
	 * @param psf
	 *            the psf
	 * @return the parameters
	 * @throws ConfigurationException
	 *             if the psf is null, or not recognised
	 */
	public static List<PSFParameter> getParameters(PSF psf) throws ConfigurationException
	{
		if (psf == null)
			throw new ConfigurationException("psf is null");

		List<PSFParameter> list = psf.getParameterList();
		switch (psf.getPsfType())
		{
			case OneAxisGaussian2D:
				return checkParameters(sxParameters, list);

			case AstigmaticGaussian2D:
			case TwoAxisGaussian2D:
				return checkParameters(sxsyParameters, list);

			case TwoAxisAndThetaGaussian2D:
				return checkParameters(sxsyaParameters, list);

			case UNRECOGNIZED:
				throw new ConfigurationException("psf is not recognised");

			case Custom:
			default:
				break;
		}

		return list;
	}

	/**
	 * Check the list has the same size and parameter type as the default and return it, otherwise return the default.
	 * This allows the names to be customised but not the types.
	 *
	 * @param defaultList
	 *            the default list
	 * @param list
	 *            the list
	 * @return the list (or the default)
	 */
	private static List<PSFParameter> checkParameters(ArrayList<PSFParameter> defaultList, List<PSFParameter> list)
	{
		if (list != null && list.size() == defaultList.size())
		{
			// Check each
			int i = 0;
			for (PSFParameter p : list)
			{
				if (p.getUnit() != defaultList.get(i++).getUnit())
					return defaultList;
			}
			return list;
		}
		return defaultList;
	}
	
	/**
	 * Creates the PSF builder.
	 *
	 * @param psfType the PSF type
	 * @return the PSF builder
	 */
	public static PSF.Builder createBuilder(PSFType psfType)
	{
		PSF.Builder builder = PSF.newBuilder();
		builder.setPsfType(psfType);
		return builder;
	}
	
	/**
	 * Creates the PSF.
	 *
	 * @param psfType the PSF type
	 * @return the PSF
	 */
	public static PSF create(PSFType psfType)
	{
		return createBuilder(psfType).build();
	}
}
