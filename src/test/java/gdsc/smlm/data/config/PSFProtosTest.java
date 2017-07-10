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
		System.out.printf(o);
		Assert.assertEquals(e, o);

		psfBuilder.clear();
		TextFormat.merge(o, psfBuilder);
		Assert.assertTrue("Merge string", psf.equals(psfBuilder.build()));

		// Short string
		String o2 = TextFormat.shortDebugString(psf);
		System.out.println(o2);

		psfBuilder.clear();
		TextFormat.merge(o2, psfBuilder);
		Assert.assertTrue("Merge short string", psf.equals(psfBuilder.build()));

		// JSON
		Printer printer = JsonFormat.printer().omittingInsignificantWhitespace();
		String json = printer.print(psf);
		System.out.println(json);
		json = JSONUtils.simplify(json);
		System.out.println(json);

		psfBuilder.clear();
		JsonFormat.parser().merge(json, psfBuilder);
		Assert.assertTrue("Merge JSON", psf.equals(psfBuilder.build()));
	}
}
