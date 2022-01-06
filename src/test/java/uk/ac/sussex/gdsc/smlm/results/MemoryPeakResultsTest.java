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

package uk.ac.sussex.gdsc.smlm.results;

import java.util.function.Supplier;
import java.util.logging.Logger;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.MemoryUtils;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;

@SuppressWarnings({"javadoc"})
class MemoryPeakResultsTest {
  private static Logger logger = Logger.getLogger(MemoryPeakResultsTest.class.getName());

  /** The margin of error when checking size against the fixed values. */
  private static final double SIZE_TOLERANCE = 0.2;

  @Test
  void canEstimatePeakResultMemorySize() {
    checkSize(MemoryPeakResults.PEAK_RESULT_SIZE,
        () -> new PeakResult(0, 0, 0, 0, 0, 0, 0, new float[PeakResult.STANDARD_PARAMETERS], null));
  }

  @Test
  void canEstimatePeakResultMemorySizeWithDeviations() {
    checkSize(MemoryPeakResults.PEAK_RESULT_SIZE_WITH_DEVIATIONS,
        () -> new PeakResult(0, 0, 0, 0, 0, 0, 0, new float[PeakResult.STANDARD_PARAMETERS],
            new float[PeakResult.STANDARD_PARAMETERS]));
  }

  /**
   * Check the size of the object.
   *
   * <p>Note: This test may not be very robust. It may fail depending on JVM platform so it just
   * logs an error but does not assert.
   *
   * @param expected the expected size
   * @param supplier the supplier of objects
   */
  private static void checkSize(int expected, Supplier<Object> supplier) {
    final long actual = MemoryUtils.measureSize(10000, supplier);
    final double error = DoubleEquality.relativeError(actual, expected);
    // This is flaky so do not assert the test
    // Assertions.assertEquals(size, expected, Math.abs(size) * SIZE_TOLERANCE);
    logger.log(TestLogUtils.getResultRecord(error < SIZE_TOLERANCE,
        "Memory expected=%d : measured=%d : error=%f", expected, actual, error));
  }

  @Test
  void canCopyAndAssignZeroIds() {
    final MemoryPeakResults results = new MemoryPeakResults(5);
    results.add(new PeakResult(0, 1, 20));
    results.add(new IdPeakResult(2, 3, 21, 99));
    results.add(new IdPeakResult(4, 5, 22, 1));
    results.add(new PeakResult(6, 7, 23));

    final MemoryPeakResults results2 = results.copyAndAssignZeroIds();
    assertResult(results2.get(0), 0, 1, 20, -1);
    assertResult(results2.get(1), 2, 3, 21, 99);
    assertResult(results2.get(2), 4, 5, 22, 1);
    assertResult(results2.get(3), 6, 7, 23, -2);
  }

  private static void assertResult(PeakResult result, float x, float y, float intensity, int id) {
    Assertions.assertEquals(x, result.getXPosition());
    Assertions.assertEquals(y, result.getYPosition());
    Assertions.assertEquals(intensity, result.getIntensity());
    Assertions.assertEquals(id, result.getId());
  }
}
