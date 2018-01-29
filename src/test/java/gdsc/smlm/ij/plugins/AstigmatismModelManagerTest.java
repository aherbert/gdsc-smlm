package gdsc.smlm.ij.plugins;

import org.junit.Assert;
import org.junit.Test;

import gdsc.core.data.utils.TypeConverter;
import gdsc.smlm.data.config.PSFProtos.AstigmatismModel;
import gdsc.smlm.data.config.UnitConverterFactory;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.function.gaussian.AstigmatismZModel;

public class AstigmatismModelManagerTest
{
	@Test
	public void canConvertModel()
	{
		DistanceUnit[] unit = new DistanceUnit[] { DistanceUnit.PIXEL, DistanceUnit.NM, DistanceUnit.UM };
		for (int i = 0; i < unit.length; i++)
			for (int j = 0; j < unit.length; j++)
				canConvertModel(unit[i], unit[j]);
	}

	private void canConvertModel(DistanceUnit zDistanceUnit, DistanceUnit sDistanceUnit)
	{
		// Use a reasonable z-depth function from the Smith, et al (2010) paper (page 377)
		double sx = 1.08;
		double sy = 1.01;
		double gamma = 0.389;
		double d = 0.531;
		double Ax = -0.0708;
		double Bx = -0.073;
		double Ay = 0.164;
		double By = 0.0417;
		double nmPerPixel = 100;

		//Ax = Ay = 0;
		//Bx = By = 0;

		AstigmatismModel.Builder builder = AstigmatismModel.newBuilder();
		builder.setGamma(gamma);
		builder.setD(d);
		builder.setS0X(sx);
		builder.setAx(Ax);
		builder.setBx(Bx);
		builder.setS0Y(sy);
		builder.setAy(Ay);
		builder.setBy(By);
		builder.setZDistanceUnit(DistanceUnit.UM);
		builder.setSDistanceUnit(DistanceUnit.PIXEL);
		builder.setNmPerPixel(nmPerPixel);

		AstigmatismModel model1 = builder.build();
		AstigmatismModel model2 = AstigmatismModelManager.convert(model1, zDistanceUnit, sDistanceUnit);

		AstigmatismZModel m1 = AstigmatismModelManager.create(model1);
		AstigmatismZModel m2 = AstigmatismModelManager.create(model2);

		TypeConverter<DistanceUnit> zc = UnitConverterFactory.createConverter(DistanceUnit.UM, zDistanceUnit,
				nmPerPixel);
		TypeConverter<DistanceUnit> sc = UnitConverterFactory.createConverter(DistanceUnit.PIXEL, sDistanceUnit,
				nmPerPixel);

		for (double z = -0.5; z <= 0.5; z += 0.1)
		{
			double e = sc.convert(m1.getSx(z));
			double o = m2.getSx(zc.convert(z));
			//System.out.printf("%f vs %f\n", e, o);
			Assert.assertEquals(e, o, e * 1e-8);
		}
	}
}
