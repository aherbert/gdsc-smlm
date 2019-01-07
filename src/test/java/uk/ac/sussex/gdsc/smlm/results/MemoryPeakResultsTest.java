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

package uk.ac.sussex.gdsc.smlm.results;

import uk.ac.sussex.gdsc.core.utils.MemoryUtils;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

@SuppressWarnings({"javadoc"})
public class MemoryPeakResultsTest {
  // Note: This test may not be very robust. It may fail depending on JVM platform.

  @Test
  public void canEstimatePeakResultMemorySize() {
    final int parameters = PeakResult.STANDARD_PARAMETERS;
    final long size = MemoryUtils.measureSize(10000,
        () -> new PeakResult(0, 0, 0, 0, 0, 0, 0, new float[parameters], null));
    Assertions.assertEquals(size, MemoryPeakResults.PEAK_RESULT_SIZE, size * 0.1);
  }

  @Test
  public void canEstimatePeakResultMemorySizeWithDeviations() {
    final int parameters = PeakResult.STANDARD_PARAMETERS;
    final long size = MemoryUtils.measureSize(10000,
        () -> new PeakResult(0, 0, 0, 0, 0, 0, 0, new float[parameters], new float[parameters]));
    Assertions.assertEquals(size, MemoryPeakResults.PEAK_RESULT_SIZE_WITH_DEVIATIONS, size * 0.1);
  }
}
