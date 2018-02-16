package gdsc.smlm.data.config;

import org.junit.Assert;
import org.junit.Test;

import gdsc.smlm.data.config.PSFProtos.AstigmatismModel;
import gdsc.smlm.data.config.PSFProtos.PSF;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;

public class PSFProtosHelperTest
{
	@Test
	public void canConvertAstigmatismModel()
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
		
		DistanceUnit zDistanceUnit = DistanceUnit.UM;
		DistanceUnit sDistanceUnit = DistanceUnit.PIXEL;

		AstigmatismModel.Builder builder = AstigmatismModel.newBuilder();
		builder.setGamma(gamma);
		builder.setD(d);
		builder.setS0X(sx);
		builder.setAx(Ax);
		builder.setBx(Bx);
		builder.setS0Y(sy);
		builder.setAy(Ay);
		builder.setBy(By);
		builder.setZDistanceUnit(zDistanceUnit);
		builder.setSDistanceUnit(sDistanceUnit);
		builder.setNmPerPixel(nmPerPixel);
		
		AstigmatismModel model1 = builder.build();
		PSF psf = PSFProtosHelper.createPSF(model1, zDistanceUnit, sDistanceUnit);
		//System.out.println(psf);

		AstigmatismModel model2 = PSFProtosHelper.createModel(psf, zDistanceUnit, sDistanceUnit, nmPerPixel);
		Assert.assertEquals(model1, model2);
	}
}
