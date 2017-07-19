package gdsc.smlm.data.config;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

import gdsc.smlm.data.config.UnitProtos.AngleUnit;
import gdsc.smlm.data.config.CalibrationProtos.CameraType;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.data.config.UnitProtos.IntensityUnit;

public class CalibrationWriterTest
{
	@Test
	public void canWrite()
	{
		RandomGenerator r = new Well19937c();
		for (int i = 0; i < 100; i++)
			canWrite(r);
	}

	private void canWrite(RandomGenerator r)
	{
		double amplification = 1 + r.nextDouble();
		double bias = 1 + r.nextDouble();
		double exposureTime = 1 + r.nextDouble();
		double gain = 1 + r.nextDouble();
		double nmPerPixel = 1 + r.nextDouble();
		double readNoise = 1 + r.nextDouble();
		AngleUnit angleUnit = AngleUnit.values()[r.nextInt(AngleUnit.values().length - 1)];
		CameraType cameraType = CameraType.values()[r.nextInt(CameraType.values().length - 1)];
		DistanceUnit distanceUnit = DistanceUnit.values()[r.nextInt(DistanceUnit.values().length - 1)];
		IntensityUnit intensityUnit = IntensityUnit.values()[r.nextInt(IntensityUnit.values().length - 1)];

		CalibrationWriter writer = new CalibrationWriter();

		Assert.assertEquals(writer.getCountPerElectron(), 0, 0);
		Assert.assertEquals(writer.getBias(), 0, 0);
		Assert.assertEquals(writer.getExposureTime(), 0, 0);
		Assert.assertEquals(writer.getCountPerPhoton(), 0, 0);
		Assert.assertEquals(writer.getNmPerPixel(), 0, 0);
		Assert.assertEquals(writer.getReadNoise(), 0, 0);
		Assert.assertFalse(writer.hasCountPerElectron());
		Assert.assertFalse(writer.hasBias());
		Assert.assertFalse(writer.hasExposureTime());
		Assert.assertFalse(writer.hasCountPerPhoton());
		Assert.assertFalse(writer.hasNmPerPixel());
		Assert.assertFalse(writer.hasReadNoise());
		Assert.assertEquals(writer.getAngleUnit(), AngleUnit.ANGLE_UNIT_NA);
		Assert.assertEquals(writer.getCameraType(), CameraType.CAMERA_TYPE_NA);
		Assert.assertEquals(writer.getDistanceUnit(), DistanceUnit.DISTANCE_UNIT_NA);
		Assert.assertEquals(writer.getIntensityUnit(), IntensityUnit.INTENSITY_UNIT_NA);

		writer.setCountPerElectron(amplification);
		writer.setBias(bias);
		writer.setExposureTime(exposureTime);
		writer.setCountPerPhoton(gain);
		writer.setNmPerPixel(nmPerPixel);
		writer.setReadNoise(readNoise);
		writer.setAngleUnit(angleUnit);
		writer.setCameraType(cameraType);
		writer.setDistanceUnit(distanceUnit);
		writer.setIntensityUnit(intensityUnit);

		Assert.assertEquals(writer.getCountPerElectron(), amplification, 0);
		Assert.assertEquals(writer.getBias(), bias, 0);
		Assert.assertEquals(writer.getExposureTime(), exposureTime, 0);
		Assert.assertEquals(writer.getCountPerPhoton(), gain, 0);
		Assert.assertEquals(writer.getNmPerPixel(), nmPerPixel, 0);
		Assert.assertEquals(writer.getReadNoise(), readNoise, 0);
		Assert.assertTrue(writer.hasCountPerElectron());
		Assert.assertTrue(writer.hasBias());
		Assert.assertTrue(writer.hasExposureTime());
		Assert.assertTrue(writer.hasCountPerPhoton());
		Assert.assertTrue(writer.hasNmPerPixel());
		Assert.assertTrue(writer.hasReadNoise());
		Assert.assertEquals(writer.getAngleUnit(), angleUnit);
		Assert.assertEquals(writer.getCameraType(), cameraType);
		Assert.assertEquals(writer.getDistanceUnit(), distanceUnit);
		Assert.assertEquals(writer.getIntensityUnit(), intensityUnit);

		CalibrationReader reader = new CalibrationReader(writer.getCalibration());

		Assert.assertEquals(reader.getCountPerElectron(), amplification, 0);
		Assert.assertEquals(reader.getBias(), bias, 0);
		Assert.assertEquals(reader.getExposureTime(), exposureTime, 0);
		Assert.assertEquals(reader.getCountPerPhoton(), gain, 0);
		Assert.assertEquals(reader.getNmPerPixel(), nmPerPixel, 0);
		Assert.assertEquals(reader.getReadNoise(), readNoise, 0);
		Assert.assertTrue(reader.hasCountPerElectron());
		Assert.assertTrue(reader.hasBias());
		Assert.assertTrue(reader.hasExposureTime());
		Assert.assertTrue(reader.hasCountPerPhoton());
		Assert.assertTrue(reader.hasNmPerPixel());
		Assert.assertTrue(reader.hasReadNoise());
		Assert.assertEquals(reader.getAngleUnit(), angleUnit);
		Assert.assertEquals(reader.getCameraType(), cameraType);
		Assert.assertEquals(reader.getDistanceUnit(), distanceUnit);
		Assert.assertEquals(reader.getIntensityUnit(), intensityUnit);
	}
}
