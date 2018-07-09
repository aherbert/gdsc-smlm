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

import gdsc.smlm.data.config.PSFProtos.AstigmatismModel;
import gdsc.smlm.data.config.PSFProtos.PSF;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;

@SuppressWarnings({ "javadoc" })
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
