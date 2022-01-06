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

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

@SuppressWarnings({"javadoc"})
class PeakResultHelperTest {
  @Test
  void canConvertLocalBackgroundToNoise() {
    final double gain = 6;

    final double[] photons = {0, 1, 2, 4, 10, 50, 100};

    // CCD
    for (final double p : photons) {
      // Assuming a Poisson distribution N photons should have a noise of sqrt(N).
      // However the input and output are in ADU counts so we apply the gain.
      final double n = PeakResultHelper.localBackgroundToNoise(p * gain, gain, false);
      Assertions.assertEquals(Math.sqrt(p) * gain, n, () -> "CCD " + p);
    }

    // EM-CCD
    for (final double p : photons) {
      // Assuming a Poisson distribution N photons should have a noise of sqrt(N * 2)
      // (due to the EM-CCD noise factor of 2).
      // However the input and output are in ADU counts so we apply the gain.
      final double n = PeakResultHelper.localBackgroundToNoise(p * gain, gain, true);
      Assertions.assertEquals(Math.sqrt(2 * p) * gain, n, () -> "EM-CCD " + p);
    }
  }

  @Test
  void canConvertLocalBackgroundToNoiseAndBack() {
    final double gain = 6;

    final double[] photons = {0, 1, 2, 4, 10, 50, 100};

    for (final boolean emCcd : new boolean[] {false, true}) {
      for (final double p : photons) {
        final double b = p * gain;
        final double n = PeakResultHelper.localBackgroundToNoise(b, gain, emCcd);
        final double b2 = PeakResultHelper.noiseToLocalBackground(n, gain, emCcd);
        Assertions.assertEquals(b, b2, 1e-6, () -> emCcd + " " + p);
      }
    }
  }
}
