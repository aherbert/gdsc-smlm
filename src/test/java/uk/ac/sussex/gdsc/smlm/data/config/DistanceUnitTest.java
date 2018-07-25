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

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;

@SuppressWarnings({ "unchecked", "javadoc" })
public class DistanceUnitTest
{
	@Test
	public void canConvert()
	{
		final double nmPerPixel = 104.5;
		for (int pixel = 1; pixel < 10; pixel++)
			//@formatter:off
    		check(nmPerPixel,
    			new ExpectedUnit<>(DistanceUnit.PIXEL, pixel),
    			new ExpectedUnit<>(DistanceUnit.UM, pixel * nmPerPixel / 1e3),
    			new ExpectedUnit<>(DistanceUnit.NM, pixel * nmPerPixel)
    			);
    		//@formatter:on
	}

	private static void check(double nmPerPixel, ExpectedUnit<DistanceUnit>... expectedUnits)
	{
		final int n = expectedUnits.length;
		TypeConverter<DistanceUnit> c;
		for (int i = 0; i < n; i++)
		{
			final DistanceUnit u1 = expectedUnits[i].u;
			final double v1 = expectedUnits[i].value;
			for (int j = 0; j < n; j++)
			{
				final DistanceUnit u2 = expectedUnits[j].u;
				c = UnitConverterFactory.createConverter(u1, u2, nmPerPixel);
				final double o = c.convert(v1);
				Assertions.assertEquals(expectedUnits[j].value, o, 1e-5, () -> u1 + " to " + u2);
			}
		}
	}
}
