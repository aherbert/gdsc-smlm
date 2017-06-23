package gdsc.smlm.data.config;

import org.junit.Assert;
import org.junit.Test;

import com.google.protobuf.InvalidProtocolBufferException;
import com.google.protobuf.util.JsonFormat;
import com.google.protobuf.util.JsonFormat.Printer;

import gdsc.smlm.data.config.PSFConfig.PSFParameter;
import gdsc.smlm.data.config.PSFConfig.PSFParameterUnit;
import gdsc.smlm.data.config.PSFConfig.PSFType;

public class PSFConfigTest
{
	@Test
	public void canWriteAndReadString()
	{
		PSFConfig.PSF.Builder psfBuilder = PSFConfig.PSF.newBuilder();
		PSFParameter.Builder psfParamBuilder = PSFConfig.PSFParameter.newBuilder();
		psfBuilder.setPsfType(PSFType.TWO_AXIS_AND_THETA_GAUSSIAN_2D);
		psfParamBuilder.setName("X SD");
		psfParamBuilder.setUnit(PSFParameterUnit.DISTANCE);
		psfBuilder.addParameter(psfParamBuilder);
		psfParamBuilder.setName("Y SD");
		psfParamBuilder.setUnit(PSFParameterUnit.DISTANCE);
		psfBuilder.addParameter(psfParamBuilder);
		psfParamBuilder.setName("Angle");
		psfParamBuilder.setUnit(PSFParameterUnit.ANGLE);
		psfBuilder.addParameter(psfParamBuilder);
		//psfBuilder.addParameterName("Y SD");
		String e = psfBuilder.toString();
		PSFConfig.PSF psf = psfBuilder.build();
		String o = psf.toString();
		Assert.assertEquals(e, o);
		//psf.getParameterName(0);
		
		// Standard string
		//System.out.printf(o);
		
		try
		{
			// JSON
			Printer printer = JsonFormat.printer().omittingInsignificantWhitespace();
			String json = printer.print(psf);
			//System.out.println(json);
			
			psfBuilder.clear();
			JsonFormat.parser().merge(json, psfBuilder);
			Assert.assertEquals(e, psfBuilder.toString());
		}
		catch (InvalidProtocolBufferException e1)
		{
			// This should be OK
			e1.printStackTrace();
		}
	}
}
