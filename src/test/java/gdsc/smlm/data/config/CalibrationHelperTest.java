package gdsc.smlm.data.config;

import java.util.ArrayList;

import org.junit.Assert;
import org.junit.Test;

import gdsc.core.data.utils.TypeConverter;
import gdsc.smlm.data.config.SMLMSettings.AngleUnit;
import gdsc.smlm.data.config.SMLMSettings.Calibration;
import gdsc.smlm.data.config.SMLMSettings.CameraCalibration;
import gdsc.smlm.data.config.SMLMSettings.DistanceCalibration;
import gdsc.smlm.data.config.SMLMSettings.DistanceUnit;
import gdsc.smlm.data.config.SMLMSettings.IntensityCalibration;
import gdsc.smlm.data.config.SMLMSettings.IntensityUnit;
import gdsc.smlm.data.config.SMLMSettings.PSFCalibration;

public class CalibrationHelperTest
{
	@Test
	public void canUpdateCalibration()
	{
		Calibration.Builder builder = Calibration.newBuilder();
		DistanceCalibration.Builder distanceBuilder = builder.getDistanceCalibrationBuilder();
		IntensityCalibration.Builder intensityBuilder = builder.getIntensityCalibrationBuilder();
		PSFCalibration.Builder psfBuilder = builder.getPsfCalibrationBuilder();
		CameraCalibration.Builder cameraBuilder = builder.getCameraCalibrationBuilder();

		double nmPerPixel = 104.5;
		double gain = 45;
		double bias = 100;

		distanceBuilder.setNmPerPixel(nmPerPixel);
		distanceBuilder.setUnit(DistanceUnit.PIXEL);

		intensityBuilder.setUnit(IntensityUnit.PHOTON);
		intensityBuilder.setGain(gain);

		psfBuilder.setAngleUnit(AngleUnit.RADIAN);
		cameraBuilder.setBias(bias);

		Calibration c = builder.build();

		CalibrationHelper helper = new CalibrationHelper(c);

		TypeConverter<DistanceUnit> distanceConverter = helper.getDistanceConverter(DistanceUnit.NM);
		Assert.assertEquals(distanceConverter.from(), DistanceUnit.PIXEL);
		Assert.assertEquals(distanceConverter.to(), DistanceUnit.NM);
		Assert.assertEquals(helper.getCalibration().getDistanceCalibration().getUnit(), DistanceUnit.NM);

		TypeConverter<DistanceUnit> distanceConverter2 = CalibrationHelper.getDistanceConverter(c, DistanceUnit.NM);
		Assert.assertEquals(distanceConverter2.from(), DistanceUnit.PIXEL);
		Assert.assertEquals(distanceConverter2.to(), DistanceUnit.NM);
		
		Assert.assertEquals(distanceConverter.getFunction(), distanceConverter2.getFunction());

		ArrayList<TypeConverter<IntensityUnit>> list = helper.getIntensityConverter(IntensityUnit.COUNT);
		Assert.assertEquals(list.get(0).from(), IntensityUnit.PHOTON);
		Assert.assertEquals(list.get(0).to(), IntensityUnit.COUNT);
		Assert.assertEquals(list.get(1).from(), IntensityUnit.PHOTON);
		Assert.assertEquals(list.get(1).to(), IntensityUnit.COUNT);
		Assert.assertEquals(helper.getCalibration().getIntensityCalibration().getUnit(), IntensityUnit.COUNT);

		ArrayList<TypeConverter<IntensityUnit>> list2 = CalibrationHelper.getIntensityConverter(c, IntensityUnit.COUNT);
		Assert.assertEquals(list2.get(0).from(), IntensityUnit.PHOTON);
		Assert.assertEquals(list2.get(0).to(), IntensityUnit.COUNT);
		Assert.assertEquals(list2.get(1).from(), IntensityUnit.PHOTON);
		Assert.assertEquals(list2.get(1).to(), IntensityUnit.COUNT);

		Assert.assertEquals(list.get(0).getFunction(), list2.get(0).getFunction());
		Assert.assertEquals(list.get(1).getFunction(), list2.get(1).getFunction());

		TypeConverter<AngleUnit> angleConverter = helper.getAngleConverter(AngleUnit.DEGREE);
		Assert.assertEquals(angleConverter.from(), AngleUnit.RADIAN);
		Assert.assertEquals(angleConverter.to(), AngleUnit.DEGREE);
		Assert.assertEquals(helper.getCalibration().getPsfCalibration().getAngleUnit(), AngleUnit.DEGREE);

		TypeConverter<AngleUnit> angleConverter2 = CalibrationHelper.getAngleConverter(c, AngleUnit.DEGREE);
		Assert.assertEquals(angleConverter2.from(), AngleUnit.RADIAN);
		Assert.assertEquals(angleConverter2.to(), AngleUnit.DEGREE);
		
		Assert.assertEquals(angleConverter.getFunction(), angleConverter2.getFunction());
	}
}
