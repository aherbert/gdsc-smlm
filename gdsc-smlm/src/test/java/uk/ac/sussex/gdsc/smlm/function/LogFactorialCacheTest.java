/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.function;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

@SuppressWarnings({"javadoc"})
class LogFactorialCacheTest {
  @Test
  void testInvalidRangeThrows() {
    final LogFactorialCache c1 = new LogFactorialCache(20);
    Assertions.assertThrows(IllegalArgumentException.class, () -> c1.ensureRange(10, 9));
  }

  @Test
  void testLogFactorial() {
    final int size = 200;
    final LogFactorialCache c1 = new LogFactorialCache();
    final LogFactorialCache c2 = new LogFactorialCache(size);
    final LogFactorialCache c3 = c1.withN(size);
    final LogFactorialCache c4 = c2.withN(0);

    final LogFactorialCache c5 = new LogFactorialCache();
    c5.increaseMaxN(size / 2);

    final LogFactorialCache c6 = new LogFactorialCache(size / 2);
    c6.ensureRange(0, size / 2);
    c6.increaseMaxN(size);
    c6.increaseMaxN(size);

    final LogFactorialCache c7 = c6.withN(size / 2);
    c7.ensureRange(0, size);
    c7.ensureRange(0, size);

    Assertions.assertEquals(size, c2.getMaxN());
    Assertions.assertEquals(size, c3.getMaxN());
    Assertions.assertEquals(size / 2, c5.getMaxN());
    Assertions.assertEquals(size, c6.getMaxN());
    Assertions.assertEquals(size, c7.getMaxN());

    for (int i = 0; i <= size; i++) {
      final int n = i;
      final double expected = LogFactorial.value(n);
      if (n <= c1.getMaxN()) {
        Assertions.assertEquals(expected, c1.getLogFactorial(n));
      } else {
        Assertions.assertThrows(ArrayIndexOutOfBoundsException.class, () -> c1.getLogFactorial(n));
      }
      Assertions.assertEquals(expected, c2.getLogFactorial(n));
      Assertions.assertEquals(expected, c2.getLogFactorial(n));
      Assertions.assertEquals(expected, c3.getLogFactorial(n));
      if (n <= c4.getMaxN()) {
        Assertions.assertEquals(expected, c4.getLogFactorial(n));
      } else {
        Assertions.assertThrows(ArrayIndexOutOfBoundsException.class, () -> c4.getLogFactorial(n));
      }
      if (n <= c5.getMaxN()) {
        Assertions.assertEquals(expected, c5.getLogFactorial(n));
      } else {
        Assertions.assertThrows(ArrayIndexOutOfBoundsException.class, () -> c5.getLogFactorial(n));
      }
      Assertions.assertEquals(n <= size / 2 ? expected : 0, c6.getLogFactorialUnsafe(n));
      Assertions.assertEquals(expected, c7.getLogFactorialUnsafe(n));
    }
  }
}
