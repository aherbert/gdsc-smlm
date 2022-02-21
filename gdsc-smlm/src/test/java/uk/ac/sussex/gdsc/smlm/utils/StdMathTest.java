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

package uk.ac.sussex.gdsc.smlm.utils;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

@SuppressWarnings({"javadoc"})
class StdMathTest {
  @ParameterizedTest
  @ValueSource(doubles = {0, Double.MIN_VALUE, Double.MIN_NORMAL, 1e-6, 1e-3, 1, 1.23, Math.PI,
      4.53e2, 789.57e6, Double.MAX_VALUE, Double.POSITIVE_INFINITY, Double.NaN})
  void testExp(double x) {
    Assertions.assertEquals(Math.exp(x), StdMath.exp(x));
    Assertions.assertEquals(Math.exp(-x), StdMath.exp(-x));
  }
}
