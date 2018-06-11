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
import gdsc.smlm.data.config.UnitProtos.AngleUnit;

@SuppressWarnings("unchecked")
public class AngleUnitTest
{
	@Test
	public void canConvert()
	{
		double degToRad = Math.PI / 180.0;
		for (int a = -360; a <= 360; a++)
		{
			//@formatter:off
    		check(degToRad, 
    			new ExpectedUnit<AngleUnit>(AngleUnit.DEGREE, a), 
    			new ExpectedUnit<AngleUnit>(AngleUnit.RADIAN, a * degToRad)
    			);
    		//@formatter:on
		}
	}

	private void check(double degToRad, ExpectedUnit<AngleUnit>... expectedUnits)
	{
		int n = expectedUnits.length;
		TypeConverter<AngleUnit> c;
		for (int i = 0; i < n; i++)
		{
			AngleUnit u1 = expectedUnits[i].u;
			double v1 = expectedUnits[i].value;
			for (int j = 0; j < n; j++)
			{
				AngleUnit u2 = expectedUnits[j].u;
				c = UnitConverterFactory.createConverter(u1, u2);
				double o = c.convert(v1);
				Assert.assertEquals(u1 + " to " + u2, expectedUnits[j].value, o, 1e-5);
			}
		}
	}
}
