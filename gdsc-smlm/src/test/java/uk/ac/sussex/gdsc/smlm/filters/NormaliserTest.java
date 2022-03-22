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

package uk.ac.sussex.gdsc.smlm.filters;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngFactory;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;

@SuppressWarnings({"javadoc"})
class NormaliserTest extends AbstractFilterTest {
  @SeededTest
  void nonNormaliserCanCopyToOutDataWithBorder(RandomSeed seed) {
    final UniformRandomProvider rg = RngFactory.create(seed.get());

    for (final int width : primes) {
      for (final int height : primes) {
        final float[] data = createData(rg, width, height);
        for (final int boxSize : boxSizes) {
          // logger.fine(() -> String.format("%dx%d : border=%d", width, height, boxSize);

          // Assume fixed normaliser works
          final FixedNormaliser n = new FixedNormaliser(1);
          final NonNormaliser nn = new NonNormaliser();
          final float[] e = new float[data.length];
          final float[] o = new float[data.length];
          n.normalise(data, e, width, height, boxSize);
          nn.normalise(data, o, width, height, boxSize);
          Assertions.assertArrayEquals(o, e,
              () -> String.format("%dx%d : border=%d", width, height, boxSize));
        }
      }
    }
  }
}
