/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
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

@SuppressWarnings({ "javadoc" })
public class CounterTest
{
    @Test
    public void canIncrementAndGet()
    {
        Counter c = new Counter();
        Assertions.assertEquals(1, c.incrementAndGet());
        Assertions.assertEquals(2, c.incrementAndGet());
        Assertions.assertEquals(3, c.incrementAndGet());

        c = new Counter();
        Assertions.assertEquals(10, c.incrementAndGet(10));
        Assertions.assertEquals(20, c.incrementAndGet(10));
        Assertions.assertEquals(30, c.incrementAndGet(10));
    }

    @Test
    public void canGetAndIncrement()
    {
        Counter c = new Counter();
        Assertions.assertEquals(0, c.getAndIncrement());
        Assertions.assertEquals(1, c.getAndIncrement());
        Assertions.assertEquals(2, c.getAndIncrement());

        c = new Counter();
        Assertions.assertEquals(0, c.getAndIncrement(10));
        Assertions.assertEquals(10, c.getAndIncrement(10));
        Assertions.assertEquals(20, c.getAndIncrement(10));
    }

    @Test
    public void canDecrementAndGet()
    {
        Counter c = new Counter();
        Assertions.assertEquals(-1, c.decrementAndGet());
        Assertions.assertEquals(-2, c.decrementAndGet());
        Assertions.assertEquals(-3, c.decrementAndGet());

        c = new Counter();
        Assertions.assertEquals(-10, c.decrementAndGet(10));
        Assertions.assertEquals(-20, c.decrementAndGet(10));
        Assertions.assertEquals(-30, c.decrementAndGet(10));
    }

    @Test
    public void canGetAndDecrement()
    {
        Counter c = new Counter();
        Assertions.assertEquals(0, c.getAndDecrement());
        Assertions.assertEquals(-1, c.getAndDecrement());
        Assertions.assertEquals(-2, c.getAndDecrement());

        c = new Counter();
        Assertions.assertEquals(0, c.getAndDecrement(10));
        Assertions.assertEquals(-10, c.getAndDecrement(10));
        Assertions.assertEquals(-20, c.getAndDecrement(10));
    }
}
