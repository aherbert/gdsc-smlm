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

import com.google.protobuf.InvalidProtocolBufferException;
import com.google.protobuf.TextFormat;
import com.google.protobuf.TextFormat.ParseException;
import com.google.protobuf.util.JsonFormat;
import com.google.protobuf.util.JsonFormat.Printer;

import gdsc.smlm.data.config.PSFProtos.PSFParameter;
import gdsc.smlm.data.config.PSFProtos.PSFParameterUnit;
import gdsc.smlm.data.config.PSFProtos.PSFType;
import gdsc.smlm.utils.JSONUtils;
import gdsc.test.TestSettings;

@SuppressWarnings({ "javadoc" })
public class PSFProtosTest
{
	@Test
	public void canWriteAndReadString() throws ParseException, InvalidProtocolBufferException
	{
		PSFProtos.PSF.Builder psfBuilder = PSFProtos.PSF.newBuilder();
		PSFParameter.Builder psfParamBuilder = PSFProtos.PSFParameter.newBuilder();
		psfBuilder.setPsfType(PSFType.TWO_AXIS_AND_THETA_GAUSSIAN_2D);
		psfParamBuilder.setName("X\"SD");
		psfParamBuilder.setUnit(PSFParameterUnit.DISTANCE);
		psfParamBuilder.setValue(1.1);
		psfBuilder.addParameters(psfParamBuilder);
		psfParamBuilder.setName("Y SD");
		psfParamBuilder.setUnit(PSFParameterUnit.DISTANCE);
		psfParamBuilder.setValue(1.2);
		psfBuilder.addParameters(psfParamBuilder);
		psfParamBuilder.setName("'Angle");
		psfParamBuilder.setUnit(PSFParameterUnit.ANGLE);
		psfParamBuilder.clearValue();
		psfBuilder.addParameters(psfParamBuilder);

		// Standard string
		String e = psfBuilder.toString();
		PSFProtos.PSF psf = psfBuilder.build();
		String o = psf.toString();
		TestSettings.debugln(o);
		Assert.assertEquals(e, o);

		psfBuilder.clear();
		TextFormat.merge(o, psfBuilder);
		Assert.assertTrue("Merge string", psf.equals(psfBuilder.build()));

		// Short string
		String o2 = TextFormat.shortDebugString(psf);
		TestSettings.debugln(o2);

		psfBuilder.clear();
		TextFormat.merge(o2, psfBuilder);
		Assert.assertTrue("Merge short string", psf.equals(psfBuilder.build()));

		// JSON
		Printer printer = JsonFormat.printer().omittingInsignificantWhitespace();
		String json = printer.print(psf);
		TestSettings.debugln(json);
		json = JSONUtils.simplify(json);
		TestSettings.debugln(json);

		psfBuilder.clear();
		JsonFormat.parser().merge(json, psfBuilder);
		Assert.assertTrue("Merge JSON", psf.equals(psfBuilder.build()));
	}
}
