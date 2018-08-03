package uk.ac.sussex.gdsc.smlm.data.config;

import java.util.logging.Level;
import java.util.logging.Logger;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import com.google.protobuf.InvalidProtocolBufferException;
import com.google.protobuf.TextFormat;
import com.google.protobuf.TextFormat.ParseException;
import com.google.protobuf.util.JsonFormat;
import com.google.protobuf.util.JsonFormat.Printer;

import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFParameter;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFParameterUnit;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.utils.JSONUtils;
import uk.ac.sussex.gdsc.test.TestLog;

@SuppressWarnings({ "javadoc" })
public class PSFProtosTest
{
    @Test
    public void canWriteAndReadString() throws ParseException, InvalidProtocolBufferException
    {
        final Logger logger = Logger.getLogger(PSFProtosTest.class.getName());
        final Level logLevel = Level.FINE;
        final PSFProtos.PSF.Builder psfBuilder = PSFProtos.PSF.newBuilder();
        final PSFParameter.Builder psfParamBuilder = PSFProtos.PSFParameter.newBuilder();
        psfBuilder.setPsfType(PSFType.CUSTOM);
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

        // Test weird characters
        psfParamBuilder.setName("other-char");
        psfBuilder.addParameters(psfParamBuilder);
        psfParamBuilder.setName("other{char");
        psfBuilder.addParameters(psfParamBuilder);
        psfParamBuilder.setName("other_char");
        psfBuilder.addParameters(psfParamBuilder);

        // Standard string
        final String e = psfBuilder.toString();
        final PSFProtos.PSF psf = psfBuilder.build();
        final String o = psf.toString();
        logger.log(TestLog.getRecord(logLevel, o));
        Assertions.assertEquals(e, o);

        psfBuilder.clear();
        TextFormat.merge(o, psfBuilder);
        Assertions.assertTrue(psf.equals(psfBuilder.build()), "Merge string");

        // Short string
        final String o2 = TextFormat.shortDebugString(psf);
        logger.log(TestLog.getRecord(logLevel, o2));

        psfBuilder.clear();
        TextFormat.merge(o2, psfBuilder);
        Assertions.assertTrue(psf.equals(psfBuilder.build()), "Merge short string");

        // JSON
        final Printer printer = JsonFormat.printer().omittingInsignificantWhitespace();
        String json = printer.print(psf);
        logger.log(TestLog.getRecord(logLevel, json));
        json = JSONUtils.simplify(json);
        logger.log(TestLog.getRecord(logLevel, json));

        psfBuilder.clear();
        JsonFormat.parser().merge(json, psfBuilder);
        Assertions.assertTrue(psf.equals(psfBuilder.build()), "Merge JSON");
    }
}
