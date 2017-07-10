package gdsc.smlm.data.config;

import gdsc.smlm.data.config.PSFProtos.PSF;
import gdsc.smlm.data.config.PSFProtos.PSFParameter;
import gdsc.smlm.data.config.PSFProtos.PSFParameterUnit;
import gdsc.smlm.data.config.PSFProtos.PSFType;

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
 * Contains helper functions for the PSFProtos class.
 */
public class PSFProtosHelper
{
	/** The default one-axis Gaussian 2D PSF */
	public static final PSF defaultOneAxisGaussian2DPSF;
	/** The default two-axis Gaussian 2D PSF */
	public static final PSF defaultTwoAxisGaussian2DPSF;
	/** The default two-axis and theta Gaussian 2D PSF */
	public static final PSF defaultTwoAxisAndThetaGaussian2DPSF;

	static
	{
		PSFParameter.Builder paramBuilder = PSFParameter.newBuilder();
		PSF.Builder builder = PSF.newBuilder();

		builder.setPsfType(PSFType.ONE_AXIS_GAUSSIAN_2D);
		paramBuilder.setName("S");
		paramBuilder.setValue(1);
		paramBuilder.setUnit(PSFParameterUnit.DISTANCE);
		builder.addParameters(paramBuilder.build());
		defaultOneAxisGaussian2DPSF = builder.build();

		builder.clear();
		builder.setPsfType(PSFType.TWO_AXIS_GAUSSIAN_2D);
		paramBuilder.setName("Sx");
		builder.addParameters(paramBuilder.build());
		paramBuilder.setName("Sy");
		builder.addParameters(paramBuilder.build());
		defaultTwoAxisGaussian2DPSF = builder.build();

		builder.setPsfType(PSFType.TWO_AXIS_AND_THETA_GAUSSIAN_2D);
		paramBuilder.setName("Angle");
		paramBuilder.setUnit(PSFParameterUnit.ANGLE);
		paramBuilder.setValue(0);
		builder.addParameters(paramBuilder.build());
		defaultTwoAxisAndThetaGaussian2DPSF = builder.build();
	}

	/**
	 * Gets the name.
	 *
	 * @param value
	 *            the value
	 * @return the name
	 */
	public static String getName(PSFType value)
	{
		switch (value)
		{
			case ASTIGMATIC_GAUSSIAN_2D:
				return "Astigmatic Gaussian 2D";
			case CUSTOM:
				return "Custom";
			case ONE_AXIS_GAUSSIAN_2D:
				return "Circular Gaussian 2D";
			case PSF_TYPE_NA:
				return "NA";
			case TWO_AXIS_AND_THETA_GAUSSIAN_2D:
				return "Rotating Elliptical Gaussian 2D";
			case TWO_AXIS_GAUSSIAN_2D:
				return "Elliptical Gaussian 2D";
			case UNRECOGNIZED:
				return "Unknown";
			default:
				throw new IllegalStateException("Unknown name: " + value);
		}
	}

	/**
	 * Gets the default PSF.
	 *
	 * @param value
	 *            the value
	 * @return the default PSF
	 */
	public static PSF getDefaultPSF(PSFType value)
	{
		switch (value)
		{
			case ONE_AXIS_GAUSSIAN_2D:
				return defaultOneAxisGaussian2DPSF;
			case TWO_AXIS_GAUSSIAN_2D:
			case ASTIGMATIC_GAUSSIAN_2D:
				return defaultTwoAxisGaussian2DPSF;
			case TWO_AXIS_AND_THETA_GAUSSIAN_2D:
				return defaultTwoAxisAndThetaGaussian2DPSF;

			case CUSTOM:
				return PSF.getDefaultInstance();
			case PSF_TYPE_NA:
				return PSF.getDefaultInstance();
			case UNRECOGNIZED:
			default:
				throw new IllegalStateException("No default PSF for type: " + value);
		}
	}
}
