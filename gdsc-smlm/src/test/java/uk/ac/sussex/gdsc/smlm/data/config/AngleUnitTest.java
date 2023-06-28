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

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit;

@SuppressWarnings({"unchecked", "javadoc"})
class AngleUnitTest {
  @Test
  void canConvert() {
    for (int a = -360; a <= 360; a++) {
      //@formatter:off
      check(
        new ExpectedUnit<>(AngleUnit.DEGREE, a),
        new ExpectedUnit<>(AngleUnit.RADIAN, Math.toRadians(a))
      );
      //@formatter:on
    }
  }

  private static void check(ExpectedUnit<AngleUnit>... expectedUnits) {
    final int n = expectedUnits.length;
    TypeConverter<AngleUnit> conv;
    for (int i = 0; i < n; i++) {
      final AngleUnit u1 = expectedUnits[i].unit;
      final double v1 = expectedUnits[i].value;
      for (int j = 0; j < n; j++) {
        final AngleUnit u2 = expectedUnits[j].unit;
        conv = UnitConverterUtils.createConverter(u1, u2);
        final double o = conv.convert(v1);
        Assertions.assertEquals(expectedUnits[j].value, o, 1e-5, () -> u1 + " to " + u2);
      }
    }
  }
}
