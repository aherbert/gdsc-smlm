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

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.smlm.results.count.Counter;

@SuppressWarnings({"javadoc"})
public class CounterTest {
  @Test
  public void canIncrementAndGet() {
    Counter counter = new Counter();
    Assertions.assertEquals(1, counter.incrementAndGet());
    Assertions.assertEquals(2, counter.incrementAndGet());
    Assertions.assertEquals(3, counter.incrementAndGet());

    counter = new Counter();
    Assertions.assertEquals(10, counter.incrementAndGet(10));
    Assertions.assertEquals(20, counter.incrementAndGet(10));
    Assertions.assertEquals(30, counter.incrementAndGet(10));
  }

  @Test
  public void canGetAndIncrement() {
    Counter counter = new Counter();
    Assertions.assertEquals(0, counter.getAndIncrement());
    Assertions.assertEquals(1, counter.getAndIncrement());
    Assertions.assertEquals(2, counter.getAndIncrement());

    counter = new Counter();
    Assertions.assertEquals(0, counter.getAndIncrement(10));
    Assertions.assertEquals(10, counter.getAndIncrement(10));
    Assertions.assertEquals(20, counter.getAndIncrement(10));
  }

  @Test
  public void canDecrementAndGet() {
    Counter counter = new Counter();
    Assertions.assertEquals(-1, counter.decrementAndGet());
    Assertions.assertEquals(-2, counter.decrementAndGet());
    Assertions.assertEquals(-3, counter.decrementAndGet());

    counter = new Counter();
    Assertions.assertEquals(-10, counter.decrementAndGet(10));
    Assertions.assertEquals(-20, counter.decrementAndGet(10));
    Assertions.assertEquals(-30, counter.decrementAndGet(10));
  }

  @Test
  public void canGetAndDecrement() {
    Counter counter = new Counter();
    Assertions.assertEquals(0, counter.getAndDecrement());
    Assertions.assertEquals(-1, counter.getAndDecrement());
    Assertions.assertEquals(-2, counter.getAndDecrement());

    counter = new Counter();
    Assertions.assertEquals(0, counter.getAndDecrement(10));
    Assertions.assertEquals(-10, counter.getAndDecrement(10));
    Assertions.assertEquals(-20, counter.getAndDecrement(10));
  }
}
