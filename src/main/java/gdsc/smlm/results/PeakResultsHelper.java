package gdsc.smlm.results;

import java.util.ArrayList;

import gdsc.core.data.utils.Converter;
import gdsc.core.data.utils.IdentityTypeConverter;
import gdsc.core.data.utils.TypeConverter;
import gdsc.core.ij.Utils;
import gdsc.core.utils.TurboList;
import gdsc.smlm.data.config.ConfigurationException;
import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.data.config.UnitHelper;
import gdsc.smlm.data.config.SMLMSettings.AngleUnit;
import gdsc.smlm.data.config.SMLMSettings.DistanceUnit;
import gdsc.smlm.data.config.SMLMSettings.IntensityUnit;
import gdsc.smlm.data.config.SMLMSettings.PSF;
import gdsc.smlm.data.config.SMLMSettings.PSFParameter;
import gdsc.smlm.data.config.SMLMSettings.PSFParameterUnit;

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
 * Helper class for converting peak results
 */
public class PeakResultsHelper
{
	private Calibration calibration;
	private PSF psf;

	public PeakResultsHelper(Calibration calibration, PSF psf)
	{
		this.calibration = calibration;
		this.psf = psf;
	}

	/** The intensity unit. */
	private IntensityUnit intensityUnit;

	/** The background converter. */
	private TypeConverter<IntensityUnit> backgroundConverter = null;
	/** The intensity converter. */
	private TypeConverter<IntensityUnit> intensityConverter = null;

	/**
	 * Gets the intensity unit.
	 *
	 * @return the intensity unit
	 */
	public IntensityUnit getIntensityUnit()
	{
		return intensityUnit;
	}

	/**
	 * Sets the intensity unit.
	 *
	 * @param intensityUnit
	 *            the new intensity unit
	 */
	public void setIntensityUnit(IntensityUnit intensityUnit)
	{
		backgroundConverter = intensityConverter = null;
		this.intensityUnit = intensityUnit;
	}

	/**
	 * Gets the background converter for the configured units. If the calibration is null then an identity converter is
	 * returned.
	 *
	 * @return the background converter
	 */
	public TypeConverter<IntensityUnit> getBackgroundConverter()
	{
		if (backgroundConverter == null)
		{
			if (calibration == null)
			{
				backgroundConverter = intensityConverter = new IdentityTypeConverter<IntensityUnit>(null);
			}
			else
			{
				ArrayList<TypeConverter<IntensityUnit>> list = calibration.getDualIntensityConverterSafe(intensityUnit);
				intensityConverter = list.get(0);
				backgroundConverter = list.get(1);
			}
		}
		return backgroundConverter;
	}

	/**
	 * Checks for intensity converter.
	 *
	 * @return true, if successful
	 */
	public boolean hasIntensityConverter()
	{
		return intensityConverter != null;
	}

	/**
	 * Gets the intensity converter for the configured units. If the calibration is null then an identity converter is
	 * returned.
	 *
	 * @return the intensity converter
	 */
	public TypeConverter<IntensityUnit> getIntensityConverter()
	{
		if (intensityConverter == null)
		{
			getBackgroundConverter();
		}
		return intensityConverter;
	}

	/** The distance unit. */
	private DistanceUnit distanceUnit;

	/** The distance converter. */
	private TypeConverter<DistanceUnit> distanceConverter = null;

	/**
	 * Gets the distance unit.
	 *
	 * @return the distance unit
	 */
	public DistanceUnit getDistanceUnit()
	{
		return distanceUnit;
	}

	/**
	 * Sets the distance unit.
	 *
	 * @param distanceUnit
	 *            the new distance unit
	 */
	public void setDistanceUnit(DistanceUnit distanceUnit)
	{
		distanceConverter = null;
		this.distanceUnit = distanceUnit;
	}

	/**
	 * Checks for distance converter.
	 *
	 * @return true, if successful
	 */
	public boolean hasDistanceConverter()
	{
		return distanceConverter != null;
	}

	/**
	 * Gets the distance converter for the configured units. If the calibration is null then an identity converter is
	 * returned.
	 *
	 * @return the distance converter
	 */
	public TypeConverter<DistanceUnit> getDistanceConverter()
	{
		if (distanceConverter == null)
		{
			distanceConverter = (calibration == null) ? new IdentityTypeConverter<DistanceUnit>(null)
					: calibration.getDistanceConverterSafe(distanceUnit);
		}
		return distanceConverter;
	}

	/** The angle unit. */
	private AngleUnit angleUnit;

	/** The angle converter. */
	private TypeConverter<AngleUnit> angleConverter = null;

	/**
	 * Gets the angle unit.
	 *
	 * @return the angle unit
	 */
	public AngleUnit getAngleUnit()
	{
		return angleUnit;
	}

	/**
	 * Sets the angle unit.
	 *
	 * @param angleUnit
	 *            the new angle unit
	 */
	public void setAngleUnit(AngleUnit angleUnit)
	{
		angleConverter = null;
		this.angleUnit = angleUnit;
	}

	/**
	 * Checks for angle converter.
	 *
	 * @return true, if successful
	 */
	public boolean hasAngleConverter()
	{
		return angleConverter != null;
	}

	/**
	 * Gets the angle converter for the configured units. If the calibration is null then an identity converter is
	 * returned.
	 *
	 * @return the angle converter
	 */
	public TypeConverter<AngleUnit> getAngleConverter()
	{
		if (angleConverter == null)
		{
			angleConverter = (calibration == null) ? new IdentityTypeConverter<AngleUnit>(null)
					: calibration.getAngleConverterSafe(angleUnit);
		}
		return angleConverter;
	}

	/**
	 * Gets the converters for the peak results parameters. This includes the standard parameters and any additional
	 * parameters defined in the PSF. If a parameter unit type is undefined then an identity converter is created.
	 *
	 * @return the converters
	 */
	public Converter[] getConverters()
	{
		TurboList<Converter> list = new TurboList<Converter>(5);

		getBackgroundConverter();
		getDistanceConverter();

		list.add(backgroundConverter);
		list.add(intensityConverter);
		list.add(distanceConverter);
		list.add(distanceConverter);
		list.add(distanceConverter);
		if (psf != null)
		{
			try
			{
				for (PSFParameter p : PSFHelper.getParameters(psf))
				{
					switch (p.getUnit())
					{
						case DISTANCE:
							list.add(distanceConverter);
							break;
						case INTENSITY:
							list.add(intensityConverter);
							break;
						case ANGLE:
							list.add(getAngleConverter());
							break;
						default:
							list.add(new IdentityTypeConverter<PSFParameterUnit>(p.getUnit()));
					}
				}
			}
			catch (ConfigurationException e)
			{
			}
		}
		return list.toArray(new Converter[list.size()]);
	}

	/**
	 * Gets the names for the peak results parameters. This includes the standard parameters and any additional
	 * parameters defined in the PSF. If a parameter name is undefined then unknown is returned.
	 *
	 * @return the converters
	 */
	public String[] getNames()
	{
		TurboList<String> list = new TurboList<String>(5);

		list.add("Background");
		list.add("Intensity");
		list.add("X");
		list.add("Y");
		list.add("Z");
		if (psf != null)
		{
			try
			{
				for (PSFParameter p : PSFHelper.getParameters(psf))
				{
					String name = p.getName();
					list.add(Utils.isNullOrEmpty(name) ? "unknown" : name);
				}
			}
			catch (ConfigurationException e)
			{
			}
		}
		return list.toArray(new String[list.size()]);
	}

	/**
	 * Gets the unit names for the peak results parameters. This includes the standard parameters and any additional
	 * parameters defined in the PSF. If a parameter unit is undefined then an empty string is returned.
	 *
	 * @return the converters
	 */
	public String[] getUnitNames()
	{
		TurboList<String> list = new TurboList<String>(5);

		getBackgroundConverter();
		getDistanceConverter();

		String intensityUnit = (intensityConverter.to() != null) ? UnitHelper.getShortName(intensityConverter.to())
				: "";
		String distanceUnit = (distanceConverter.to() != null) ? UnitHelper.getShortName(distanceConverter.to()) : "";
		String angleUnit = null;

		list.add(intensityUnit);
		list.add(intensityUnit);
		list.add(distanceUnit);
		list.add(distanceUnit);
		list.add(distanceUnit);
		if (psf != null)
		{
			try
			{
				for (PSFParameter p : PSFHelper.getParameters(psf))
				{
					switch (p.getUnit())
					{
						case DISTANCE:
							list.add(distanceUnit);
							break;
						case INTENSITY:
							list.add(intensityUnit);
							break;
						case ANGLE:
							list.add(angleUnit = getAngleUnit(angleUnit));
							break;
						default:
							list.add("");
					}
				}
			}
			catch (ConfigurationException e)
			{
			}
		}
		return list.toArray(new String[list.size()]);
	}

	private String getAngleUnit(String angleUnit)
	{
		if (angleUnit == null)
		{
			getAngleConverter();
			angleUnit = (angleConverter.to() != null) ? angleUnit = UnitHelper.getShortName(angleConverter.to()) : "";

		}
		return angleUnit;
	}

	/**
	 * Gets the calibration, updated with the current output units of the converters.
	 *
	 * @return the calibration
	 */
	public Calibration getCalibration()
	{
		if (calibration == null)
			return null;

		// Clone the calibration as it may change
		Calibration calibration = this.calibration.clone();

		if (hasIntensityConverter())
			calibration.setIntensityUnit(getIntensityConverter().to());
		if (hasDistanceConverter())
			calibration.setDistanceUnit(getDistanceConverter().to());
		if (hasAngleConverter())
			calibration.setAngleUnit(getAngleConverter().to());

		return calibration;
	}
}
