package gdsc.smlm.data.config;

import org.junit.Assert;
import org.junit.Test;

import gdsc.core.data.utils.TypeConverter;
import gdsc.smlm.data.config.SMLMSettings.AngleUnit;
import gdsc.smlm.data.config.SMLMSettings.Calibration;
import gdsc.smlm.data.config.SMLMSettings.DistanceCalibration;
import gdsc.smlm.data.config.SMLMSettings.DistanceUnit;
import gdsc.smlm.data.config.SMLMSettings.IntensityCalibration;
import gdsc.smlm.data.config.SMLMSettings.IntensityUnit;
import gdsc.smlm.data.config.SMLMSettings.PSFCalibration;

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
		intensityBuilder.setGain(gain);

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
		PSFCalibration.Builder psfBuilder = builder.getPsfCalibrationBuilder();

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
