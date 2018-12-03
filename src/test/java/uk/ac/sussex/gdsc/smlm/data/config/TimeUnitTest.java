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
package uk.ac.sussex.gdsc.smlm.data.config;

import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.TimeUnit;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

@SuppressWarnings({ "unchecked", "javadoc" })
public class TimeUnitTest
{
    @Test
    public void canConvert()
    {
        final double msPerFrame = 35;
        for (int frame = 1; frame < 10; frame++)
            //@formatter:off
    		check(msPerFrame,
    			new ExpectedUnit<>(TimeUnit.FRAME, frame),
    			new ExpectedUnit<>(TimeUnit.SECOND, frame * msPerFrame / 1e3),
    			new ExpectedUnit<>(TimeUnit.MILLISECOND, frame * msPerFrame)
    			);
    		//@formatter:on
    }

    private static void check(double msPerFrame, ExpectedUnit<TimeUnit>... expectedUnits)
    {
        final int n = expectedUnits.length;
        TypeConverter<TimeUnit> c;
        for (int i = 0; i < n; i++)
        {
            final TimeUnit u1 = expectedUnits[i].u;
            final double v1 = expectedUnits[i].value;
            for (int j = 0; j < n; j++)
            {
                final TimeUnit u2 = expectedUnits[j].u;
                c = UnitConverterFactory.createConverter(u1, u2, msPerFrame);
                final double o = c.convert(v1);
                Assertions.assertEquals(expectedUnits[j].value, o, 1e-5, () -> u1 + " to " + u2);
            }
        }
    }
}
