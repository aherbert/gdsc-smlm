package gdsc.smlm.data.config;

import java.util.List;

import gdsc.smlm.data.config.PSFProtos.PSF;
import gdsc.smlm.data.config.PSFProtos.PSFOrBuilder;
import gdsc.smlm.data.config.PSFProtos.PSFParameter;
import gdsc.smlm.data.config.PSFProtos.PSFParameterUnit;
import gdsc.smlm.data.config.PSFProtos.PSFType;
import gdsc.smlm.results.PeakResult;

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
	/** The index in the Gaussian 2D PSF parameters for the first axis standard deviation (Sx). */
	public static final int INDEX_SX = 0;
	/** The index in the Gaussian 2D PSF parameters for the second axis standard deviation (Sy). */
	public static final int INDEX_SY = 1;
	/** The index in the Gaussian 2D PSF parameters for the theta rotation angle. */
	public static final int INDEX_THETA = 2;

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
	public static boolean isGaussian2D(PSFOrBuilder psf) throws ConfigurationException
	{
		if (psf == null)
			throw new ConfigurationException("psf is null");
		switch (psf.getPsfType())
		{
			case ONE_AXIS_GAUSSIAN_2D:
			case ASTIGMATIC_GAUSSIAN_2D:
			case TWO_AXIS_AND_THETA_GAUSSIAN_2D:
			case TWO_AXIS_GAUSSIAN_2D:
				return true;
			case CUSTOM:
			case UNRECOGNIZED:
			default:
				break;
		}
		return false;
	}

	/**
	 * Gets the Gaussian 2D x-width and y-width indices for the PeakResult parameters.
	 * <p>
	 * Note that the indices can be used directly with the PeakResult parameters array as they have been adjusted using
	 * an offset of PeakResult.STANDARD_PARAMETERS.
	 *
	 * @param psf
	 *            the psf
	 * @return the Gaussian 2D x-width and y-width indices for the PeakResult parameters.
	 * @throws ConfigurationException
	 *             if the psf is null, or not a Gaussian 2D function
	 */
	public static int[] getGaussian2DWxWyIndices(PSFOrBuilder psf) throws ConfigurationException
	{
		if (psf == null)
			throw new ConfigurationException("psf is null");
		switch (psf.getPsfType())
		{
			case ONE_AXIS_GAUSSIAN_2D:
				return new int[] { PeakResult.STANDARD_PARAMETERS + INDEX_SX,
						PeakResult.STANDARD_PARAMETERS + INDEX_SX };
			case ASTIGMATIC_GAUSSIAN_2D:
			case TWO_AXIS_AND_THETA_GAUSSIAN_2D:
			case TWO_AXIS_GAUSSIAN_2D:
				return new int[] { PeakResult.STANDARD_PARAMETERS + INDEX_SX,
						PeakResult.STANDARD_PARAMETERS + INDEX_SY };
			case CUSTOM:
			case UNRECOGNIZED:
			default:
				break;
		}
		throw new ConfigurationException("psf is not Gaussian2D");
	}

	/**
	 * Gets the Gaussian 2D x-width and y-width for the PSF parameters.
	 *
	 * @param psf
	 *            the psf
	 * @return the Gaussian 2D x-width and y-width for the PSF parameters.
	 * @throws ConfigurationException
	 *             if the psf is null, or not a Gaussian 2D function
	 */
	public static double[] getGaussian2DWxWy(PSFOrBuilder psf) throws ConfigurationException
	{
		if (psf == null)
			throw new ConfigurationException("psf is null");
		double sx, sy;
		switch (psf.getPsfType())
		{
			case ONE_AXIS_GAUSSIAN_2D:
				sx = getParameterValue(psf, INDEX_SX, 1);
				return new double[] { sx, sx };
			case ASTIGMATIC_GAUSSIAN_2D:
			case TWO_AXIS_AND_THETA_GAUSSIAN_2D:
			case TWO_AXIS_GAUSSIAN_2D:
				sx = getParameterValue(psf, INDEX_SX, 1);
				sy = getParameterValue(psf, INDEX_SY, 1);
				return new double[] { sx, sy };
			case CUSTOM:
			case UNRECOGNIZED:
			default:
				break;
		}
		throw new ConfigurationException("psf is not Gaussian2D");
	}

	private static double getParameterValue(PSFOrBuilder psf, int i, double defaultValue)
	{
		if (psf.getParametersCount() > i)
		{
			double v = psf.getParameters(i).getValue();
			if (v > 0)
				return v;
		}
		return defaultValue;
	}


	/**
	 * Gets the Gaussian 2D x-width for the PSF parameters.
	 *
	 * @param psf
	 *            the psf
	 * @return the Gaussian 2D x-width for the PSF parameters.
	 * @throws ConfigurationException
	 *             if the psf is null, or not a one-axis Gaussian 2D function
	 */
	public static double getGaussian2DWx(PSFOrBuilder psf) throws ConfigurationException
	{
		if (psf == null)
			throw new ConfigurationException("psf is null");
		switch (psf.getPsfType())
		{
			case ONE_AXIS_GAUSSIAN_2D:
				return getParameterValue(psf, INDEX_SX, 1);
			case ASTIGMATIC_GAUSSIAN_2D:
			case TWO_AXIS_AND_THETA_GAUSSIAN_2D:
			case TWO_AXIS_GAUSSIAN_2D:
				// These are not supported 
				throw new ConfigurationException("Gaussian2D psf is not a one-axis Gaussian");
			case CUSTOM:
			case UNRECOGNIZED:
			default:
				break;
		}
		throw new ConfigurationException("psf is not Gaussian2D");
	}
	
	/**
	 * Gets the Gaussian 2D angle index for the PeakResult parameters.
	 * <p>
	 * Note that the index can be used directly with the PeakResult parameters array as they have been adjusted using
	 * an offset of PeakResult.STANDARD_PARAMETERS.
	 *
	 * @param psf
	 *            the psf
	 * @return the Gaussian 2D x-width and y-width indices for the PeakResult parameters.
	 * @throws ConfigurationException
	 *             if the psf is null, or not a rotated two axis Gaussian 2D function
	 */
	public static int getGaussian2DAngleIndex(PSFOrBuilder psf) throws ConfigurationException
	{
		if (psf == null)
			throw new ConfigurationException("psf is null");
		switch (psf.getPsfType())
		{
			case TWO_AXIS_AND_THETA_GAUSSIAN_2D:
				return PeakResult.STANDARD_PARAMETERS + INDEX_THETA;
			default:
				break;
		}
		throw new ConfigurationException("psf is not a rotated two axis Gaussian2D");
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
	public static List<PSFParameter> getParameters(PSFOrBuilder psf) throws ConfigurationException
	{
		if (psf == null)
			throw new ConfigurationException("psf is null");

		List<PSFParameter> list = psf.getParametersList();
		switch (psf.getPsfType())
		{
			case ONE_AXIS_GAUSSIAN_2D:
				return checkParameters(PSFProtosHelper.defaultOneAxisGaussian2DPSF.getParametersList(), list);

			case ASTIGMATIC_GAUSSIAN_2D:
			case TWO_AXIS_GAUSSIAN_2D:
				return checkParameters(PSFProtosHelper.defaultTwoAxisGaussian2DPSF.getParametersList(), list);

			case TWO_AXIS_AND_THETA_GAUSSIAN_2D:
				return checkParameters(PSFProtosHelper.defaultTwoAxisAndThetaGaussian2DPSF.getParametersList(), list);

			case UNRECOGNIZED:
				throw new ConfigurationException("psf is not recognised");

			case CUSTOM:
			default:
				break;
		}

		return list;
	}

	/**
	 * Gets the count of parameters for the PSF.
	 *
	 * @param psf
	 *            the psf
	 * @return the parameter count
	 * @throws ConfigurationException
	 *             if the psf is null, or not recognised
	 */
	public static int getParameterCount(PSFOrBuilder psf) throws ConfigurationException
	{
		if (psf == null)
			throw new ConfigurationException("psf is null");

		switch (psf.getPsfType())
		{
			case ONE_AXIS_GAUSSIAN_2D:
				return 1;

			case ASTIGMATIC_GAUSSIAN_2D:
			case TWO_AXIS_GAUSSIAN_2D:
				return 2;

			case TWO_AXIS_AND_THETA_GAUSSIAN_2D:
				return 3;

			case UNRECOGNIZED:
				throw new ConfigurationException("psf is not recognised");

			case CUSTOM:
			default:
				return psf.getParametersCount();
		}
	}

	/**
	 * Gets the count of parameters for the PSF.
	 *
	 * @param psf
	 *            the psf
	 * @return the parameter count (may be zero if the PSF is not configured)
	 */
	public static int getParameterCountSafe(PSFOrBuilder psf)
	{
		try
		{
			return getParameterCount(psf);
		}
		catch (ConfigurationException e)
		{
			return 0;
		}
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
	private static List<PSFParameter> checkParameters(List<PSFParameter> defaultList, List<PSFParameter> list)
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
	 * @param psfType
	 *            the PSF type
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
	 * @param psfType
	 *            the PSF type
	 * @return the PSF
	 */
	public static PSF create(PSFType psfType)
	{
		return createBuilder(psfType).build();
	}

	/**
	 * Checks if is a 3D.
	 *
	 * @param psf
	 *            the psf
	 * @return true, if is 3D
	 */
	public static boolean is3D(PSFOrBuilder psf)
	{
		if (psf != null)
		{
			return psf.getPsfType() == PSFType.ASTIGMATIC_GAUSSIAN_2D;
		}
		return false;
	}

	/**
	 * Checks if is two axis gaussian 2 D.
	 *
	 * @param psfType
	 *            the psf type
	 * @return true, if is two axis gaussian 2 D
	 */
	public static boolean isTwoAxisGaussian2D(PSFType psfType)
	{
		switch (psfType)
		{
			case ASTIGMATIC_GAUSSIAN_2D:
			case TWO_AXIS_GAUSSIAN_2D:
			case TWO_AXIS_AND_THETA_GAUSSIAN_2D:
				return true;
			default:
				return false;
		}
	}

	/**
	 * Checks for angle parameters.
	 *
	 * @param psf
	 *            the psf
	 * @return true, if successful
	 */
	public static boolean hasAngleParameters(PSFOrBuilder psf)
	{
		if (psf != null)
		{
			for (PSFParameter p : getParameters(psf))
				if (p.getUnit() == PSFParameterUnit.ANGLE)
					return true;
		}
		return false;
	}
}
