/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.ij.plugins;

import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.AstigmatismModel;
import uk.ac.sussex.gdsc.smlm.data.config.UnitConverterUtils;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.function.gaussian.AstigmatismZModel;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;

import org.junit.jupiter.api.Test;

@SuppressWarnings({"javadoc"})
public class AstigmatismModelManagerTest {
  @Test
  public void canConvertModel() {
    final DistanceUnit[] unit =
        new DistanceUnit[] {DistanceUnit.PIXEL, DistanceUnit.NM, DistanceUnit.UM};
    for (int i = 0; i < unit.length; i++) {
      for (int j = 0; j < unit.length; j++) {
        canConvertModel(unit[i], unit[j]);
      }
    }
  }

  private static void canConvertModel(DistanceUnit zDistanceUnit, DistanceUnit sdDistanceUnit) {
    // Use a reasonable z-depth function from the Smith, et al (2010) paper (page 377)
    final double sx = 1.08;
    final double sy = 1.01;
    final double gamma = 0.389;
    final double d = 0.531;
    final double Ax = -0.0708;
    final double Bx = -0.073;
    final double Ay = 0.164;
    final double By = 0.0417;
    final double nmPerPixel = 100;

    // Ax = Ay = 0;
    // Bx = By = 0;

    final AstigmatismModel.Builder builder = AstigmatismModel.newBuilder();
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

    final AstigmatismModel model1 = builder.build();
    final AstigmatismModel model2 =
        AstigmatismModelManager.convert(model1, zDistanceUnit, sdDistanceUnit);

    final AstigmatismZModel m1 = AstigmatismModelManager.create(model1);
    final AstigmatismZModel m2 = AstigmatismModelManager.create(model2);

    final TypeConverter<DistanceUnit> zc =
        UnitConverterUtils.createConverter(DistanceUnit.UM, zDistanceUnit, nmPerPixel);
    final TypeConverter<DistanceUnit> sc =
        UnitConverterUtils.createConverter(DistanceUnit.PIXEL, sdDistanceUnit, nmPerPixel);

    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-8, 0);
    for (double z = -0.5; z <= 0.5; z += 0.1) {
      final double e = sc.convert(m1.getSx(z));
      final double o = m2.getSx(zc.convert(z));
      // logger.fine(FunctionUtils.getSupplier("%f vs %f", e, o);
      TestAssertions.assertTest(e, o, predicate);
    }
  }
}
