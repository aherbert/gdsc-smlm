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
package gdsc.smlm.data.config;

import org.junit.Assert;
import org.junit.Test;

import gdsc.core.data.utils.TypeConverter;
import gdsc.smlm.data.config.UnitProtos.IntensityUnit;

@SuppressWarnings({ "unchecked", "javadoc" })
public class IntensityUnitTest
{
	@Test
	public void canConvert()
	{
		double offset = 120;
		double countPerPhoton = 45.5;
		for (int photon = 1; photon < 100; photon++)
		{
			//@formatter:off
    		check(offset, countPerPhoton,
    			new ExpectedUnit<>(IntensityUnit.COUNT, offset + photon * countPerPhoton),
    			new ExpectedUnit<>(IntensityUnit.PHOTON, photon)
    			);
    		check(0, countPerPhoton,
        			new ExpectedUnit<>(IntensityUnit.COUNT, photon * countPerPhoton),
        			new ExpectedUnit<>(IntensityUnit.PHOTON, photon)
        			);
    		check(countPerPhoton,
        			new ExpectedUnit<>(IntensityUnit.COUNT, photon * countPerPhoton),
        			new ExpectedUnit<>(IntensityUnit.PHOTON, photon)
        			);
    		//@formatter:on
		}
	}

	private void check(double offset, double countPerPhoton, ExpectedUnit<IntensityUnit>... expectedUnits)
	{
		int n = expectedUnits.length;
		TypeConverter<IntensityUnit> c;
		for (int i = 0; i < n; i++)
		{
			IntensityUnit u1 = expectedUnits[i].u;
			double v1 = expectedUnits[i].value;
			for (int j = 0; j < n; j++)
			{
				IntensityUnit u2 = expectedUnits[j].u;
				c = UnitConverterFactory.createConverter(u1, u2, offset, countPerPhoton);
				double o = c.convert(v1);
				Assert.assertEquals(u1 + " to " + u2, expectedUnits[j].value, o, 1e-5);
			}
		}
	}

	private void check(double countPerPhoton, ExpectedUnit<IntensityUnit>... expectedUnits)
	{
		int n = expectedUnits.length;
		TypeConverter<IntensityUnit> c;
		for (int i = 0; i < n; i++)
		{
			IntensityUnit u1 = expectedUnits[i].u;
			double v1 = expectedUnits[i].value;
			for (int j = 0; j < n; j++)
			{
				IntensityUnit u2 = expectedUnits[j].u;
				c = UnitConverterFactory.createConverter(u1, u2, countPerPhoton);
				double o = c.convert(v1);
				Assert.assertEquals(u1 + " to " + u2, expectedUnits[j].value, o, 1e-5);
			}
		}
	}
}
