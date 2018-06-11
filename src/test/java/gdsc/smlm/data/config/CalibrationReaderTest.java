/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 * 
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package gdsc.smlm.data.config;

import org.junit.Assert;
import org.junit.Test;

import gdsc.core.data.utils.TypeConverter;
import gdsc.smlm.data.config.CalibrationProtos.AngleCalibration;
import gdsc.smlm.data.config.CalibrationProtos.Calibration;
import gdsc.smlm.data.config.CalibrationProtos.DistanceCalibration;
import gdsc.smlm.data.config.CalibrationProtos.IntensityCalibration;
import gdsc.smlm.data.config.UnitProtos.AngleUnit;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.data.config.UnitProtos.IntensityUnit;

public class CalibrationReaderTest
{
	double nmPerPixel = 104.5;
	double gain = 45;
	double bias = 100;

	@Test
	public void canGetDistanceConverter()
	{
		Calibration.Builder builder = Calibration.newBuilder();
		DistanceCalibration.Builder distanceBuilder = builder.getDistanceCalibrationBuilder();

		distanceBuilder.setNmPerPixel(nmPerPixel);
		distanceBuilder.setDistanceUnit(DistanceUnit.PIXEL);

		Calibration c = builder.build();

		CalibrationReader reader = new CalibrationReader(c);

		TypeConverter<DistanceUnit> distanceConverter = reader.getDistanceConverter(DistanceUnit.NM);
		Assert.assertEquals(distanceConverter.from(), DistanceUnit.PIXEL);
		Assert.assertEquals(distanceConverter.to(), DistanceUnit.NM);

		TypeConverter<DistanceUnit> distanceConverter2 = CalibrationHelper.getDistanceConverter(c, DistanceUnit.NM);
		Assert.assertEquals(distanceConverter2.from(), DistanceUnit.PIXEL);
		Assert.assertEquals(distanceConverter2.to(), DistanceUnit.NM);

		Assert.assertEquals(distanceConverter.getFunction(), distanceConverter2.getFunction());
	}

	@Test
	public void canGetIntensityConverter()
	{
		Calibration.Builder builder = Calibration.newBuilder();
		IntensityCalibration.Builder intensityBuilder = builder.getIntensityCalibrationBuilder();

		intensityBuilder.setIntensityUnit(IntensityUnit.PHOTON);
		intensityBuilder.setCountPerPhoton(gain);

		Calibration c = builder.build();

		CalibrationReader reader = new CalibrationReader(c);

		TypeConverter<IntensityUnit> intensityConverter = reader.getIntensityConverter(IntensityUnit.COUNT);
		Assert.assertEquals(intensityConverter.from(), IntensityUnit.PHOTON);
		Assert.assertEquals(intensityConverter.to(), IntensityUnit.COUNT);

		TypeConverter<IntensityUnit> intensityConverter2 = CalibrationHelper.getIntensityConverter(c,
				IntensityUnit.COUNT);
		Assert.assertEquals(intensityConverter2.from(), IntensityUnit.PHOTON);
		Assert.assertEquals(intensityConverter2.to(), IntensityUnit.COUNT);

		Assert.assertEquals(intensityConverter.getFunction(), intensityConverter2.getFunction());
	}

	@Test
	public void canGetAngleConverter()
	{
		Calibration.Builder builder = Calibration.newBuilder();
		AngleCalibration.Builder psfBuilder = builder.getAngleCalibrationBuilder();

		psfBuilder.setAngleUnit(AngleUnit.RADIAN);

		Calibration c = builder.build();

		CalibrationReader reader = new CalibrationReader(c);

		TypeConverter<AngleUnit> angleConverter = reader.getAngleConverter(AngleUnit.DEGREE);
		Assert.assertEquals(angleConverter.from(), AngleUnit.RADIAN);
		Assert.assertEquals(angleConverter.to(), AngleUnit.DEGREE);

		TypeConverter<AngleUnit> angleConverter2 = CalibrationHelper.getAngleConverter(c, AngleUnit.DEGREE);
		Assert.assertEquals(angleConverter2.from(), AngleUnit.RADIAN);
		Assert.assertEquals(angleConverter2.to(), AngleUnit.DEGREE);

		Assert.assertEquals(angleConverter.getFunction(), angleConverter2.getFunction());
	}
}
