/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;

@SuppressWarnings({"unchecked", "javadoc"})
class IntensityUnitTest {
  @Test
  void canConvert() {
    final double offset = 120;
    final double countPerPhoton = 45.5;
    for (int photon = 1; photon < 100; photon++) {
      //@formatter:off
      check(offset, countPerPhoton,
        new ExpectedUnit<>(IntensityUnit.COUNT, offset + photon * countPerPhoton),
        new ExpectedUnit<>(IntensityUnit.PHOTON, photon)
      );
      check(0, countPerPhoton,
            new ExpectedUnit<>(IntensityUnit.COUNT, photon * countPerPhoton),
            new ExpectedUnit<>(IntensityUnit.PHOTON, photon)
      );
      check(countPerPhoton,
            new ExpectedUnit<>(IntensityUnit.COUNT, photon * countPerPhoton),
            new ExpectedUnit<>(IntensityUnit.PHOTON, photon)
      );
      //@formatter:on
    }
  }

  private static void check(double offset, double countPerPhoton,
      ExpectedUnit<IntensityUnit>... expectedUnits) {
    final int n = expectedUnits.length;
    TypeConverter<IntensityUnit> conv;
    for (int i = 0; i < n; i++) {
      final IntensityUnit u1 = expectedUnits[i].unit;
      final double v1 = expectedUnits[i].value;
      for (int j = 0; j < n; j++) {
        final IntensityUnit u2 = expectedUnits[j].unit;
        conv = UnitConverterUtils.createConverter(u1, u2, offset, countPerPhoton);
        final double o = conv.convert(v1);
        Assertions.assertEquals(expectedUnits[j].value, o, 1e-5, () -> u1 + " to " + u2);
      }
    }
  }

  private static void check(double countPerPhoton, ExpectedUnit<IntensityUnit>... expectedUnits) {
    final int n = expectedUnits.length;
    TypeConverter<IntensityUnit> conv;
    for (int i = 0; i < n; i++) {
      final IntensityUnit u1 = expectedUnits[i].unit;
      final double v1 = expectedUnits[i].value;
      for (int j = 0; j < n; j++) {
        final IntensityUnit u2 = expectedUnits[j].unit;
        conv = UnitConverterUtils.createConverter(u1, u2, countPerPhoton);
        final double o = conv.convert(v1);
        Assertions.assertEquals(expectedUnits[j].value, o, 1e-5, () -> u1 + " to " + u2);
      }
    }
  }
}
