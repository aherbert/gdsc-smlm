/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.data.config;

import com.google.protobuf.InvalidProtocolBufferException;
import com.google.protobuf.TextFormat;
import com.google.protobuf.TextFormat.ParseException;
import com.google.protobuf.util.JsonFormat;
import com.google.protobuf.util.JsonFormat.Printer;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFParameter;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFParameterUnit;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.utils.JsonUtils;
import uk.ac.sussex.gdsc.test.utils.TestLogging;
import uk.ac.sussex.gdsc.test.utils.TestLogging.TestLevel;

@SuppressWarnings({"javadoc"})
class PsfProtosTest {
  @Test
  void canWriteAndReadString() throws ParseException, InvalidProtocolBufferException {
    final Logger logger = Logger.getLogger(PsfProtosTest.class.getName());
    final Level logLevel = TestLevel.TEST_INFO;
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
    logger.log(TestLogging.getRecord(logLevel, o));
    Assertions.assertEquals(e, o);

    psfBuilder.clear();
    TextFormat.merge(o, psfBuilder);
    Assertions.assertTrue(psf.equals(psfBuilder.build()), "Merge string");

    // Short string
    final String o2 = TextFormat.shortDebugString(psf);
    logger.log(TestLogging.getRecord(logLevel, o2));

    psfBuilder.clear();
    TextFormat.merge(o2, psfBuilder);
    Assertions.assertTrue(psf.equals(psfBuilder.build()), "Merge short string");

    // JSON
    final Printer printer = JsonFormat.printer().omittingInsignificantWhitespace();
    String json = printer.print(psf);
    logger.log(TestLogging.getRecord(logLevel, json));
    json = JsonUtils.simplify(json);
    logger.log(TestLogging.getRecord(logLevel, json));

    psfBuilder.clear();
    JsonFormat.parser().merge(json, psfBuilder);
    Assertions.assertTrue(psf.equals(psfBuilder.build()), "Merge JSON");
  }
}
